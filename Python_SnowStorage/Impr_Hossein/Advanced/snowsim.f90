! ============================================================
!  gfortran -g -fcheck=all -fbacktrace -Wall snowsim.f90 dopri5.f90 YBER_ODEPACK.f90 ODEPACK_MODULES.f90 -o snowsim -lopenblas 
!  gfortran snowsim.f90 dopri5.f90 YBER_ODEPACK.f90 ODEPACK_MODULES.f90 -o snowsim -lopenblas
! ============================================================
module enhanced_features_module
    implicit none

    type InsulationParameters
        double precision :: Hi, k_dry, k_sat, n_k
        double precision :: W_sat, W_field
        double precision :: alpha_dry, alpha_wet, n_alpha
        double precision :: delta_k_age, tau_k_years
        double precision :: delta_alpha_age, tau_alpha_years
        double precision :: zeta0, gamma_H, gamma_W
        double precision :: beta_w, K_E, K_D
        double precision :: Lv, rho_w, c_w, Tfreeze
        double precision :: rho_air, C_E, U10, P0
    end type InsulationParameters
    
    type InsulationState
        double precision :: W, age_days             ! State variables
        double precision :: k_eff, alpha_eff, f_sat ! Diagnostic outputs
    end type InsulationState
    
    type Forcing
        double precision :: Isolar, Prain, T_rain, RH, Ta
    end type Forcing
    
    type LayerProperties
        double precision :: LWC              ! Liquid water content [-]
        double precision :: ice_fraction     ! Ice volume fraction [-]
        double precision :: height           ! Layer height [m]
    end type LayerProperties
    
    type CSVData
        integer :: n_points
        double precision, allocatable :: time(:)      ! Hours from start
        double precision, allocatable :: temp(:)      ! Air temperature [°C]
        double precision, allocatable :: wind(:)      ! Wind speed [m/s]
        double precision, allocatable :: precip(:)    ! Precipitation [m/h]
        double precision, allocatable :: solar(:)     ! Solar irradiance [W/m²]
        double precision, allocatable :: rh(:)        ! Relative humidity [%]
    end type CSVData
    
    ! Module variables for insulation state persistence
    double precision, allocatable :: T_ins_stored(:)
    logical :: T_ins_initialized = .false.

contains

    pure function interpolate_data(data_vec, t_query, dt_data_val) result(interp_val)
        implicit none
        double precision, intent(in) :: data_vec(:)
        double precision, intent(in) :: t_query
        double precision, intent(in) :: dt_data_val
        double precision             :: interp_val
        
        double precision :: idx_float, frac
        integer :: idx_low, idx_high, n
        
        n = size(data_vec)
        
        ! Calculate float index (0-based like Python)
        idx_float = t_query / dt_data_val
        
        ! Calculate integer bounds (convert to 1-based Fortran indexing)
        idx_low = int(floor(idx_float)) + 1
        idx_high = int(ceiling(idx_float)) + 1
        
        ! CLAMPING: ensure indices are between 1 and n
        if (idx_low < 1) idx_low = 1
        if (idx_high < 1) idx_high = 1
        if (idx_low > n) idx_low = n
        if (idx_high > n) idx_high = n
        
        ! Handle exact matches
        if (idx_low == idx_high) then
            interp_val = data_vec(idx_low)
            return
        end if
        
        ! Linear interpolation
        frac = idx_float - floor(idx_float)
        interp_val = data_vec(idx_low) * (1.0_8 - frac) + data_vec(idx_high) * frac
                        
    end function interpolate_data
    
end module enhanced_features_module

! ============================================================
!  CSV DATA READING SUBROUTINE
! ============================================================
subroutine read_csv_data(filename, csv_data, status)
    use enhanced_features_module
    implicit none
    
    character(len=*), intent(in) :: filename
    type(CSVData), intent(out) :: csv_data
    integer, intent(out) :: status
    
    integer :: unit_num, ios, n_lines, i
    character(len=500) :: line_buffer
    character(len=30)  :: dummy_timestamp ! To hold the "2023-04-01T00:00" string
    
    status = 0
    unit_num = 15 ! Changed to a safer unit number
    
    ! 1. First pass: count valid data lines
    open(unit=unit_num, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
        status = -1
        return
    end if
    
    read(unit_num, '(A)', iostat=ios) line_buffer  ! Skip header line
    n_lines = 0
    do
        ! Read the whole line as a string to count it safely
        read(unit_num, '(A)', iostat=ios) line_buffer
        if (ios /= 0) exit
        if (len_trim(line_buffer) > 0) n_lines = n_lines + 1
    end do
    rewind(unit_num)
    read(unit_num, '(A)') line_buffer ! Skip header again
    
    ! 2. Allocate arrays
    csv_data%n_points = n_lines
    if (n_lines == 0) then
        status = -3
        close(unit_num)
        return
    end if
    
    allocate(csv_data%time(n_lines), csv_data%temp(n_lines), &
             csv_data%wind(n_lines), csv_data%precip(n_lines), &
             csv_data%solar(n_lines), csv_data%rh(n_lines))
    
    ! 3. Second pass: read numeric data
    do i = 1, n_lines
        ! Use a comma-separated read: first field is string, others are numbers
        read(unit_num, *, iostat=ios) dummy_timestamp, &
             csv_data%temp(i), csv_data%wind(i), csv_data%precip(i), &
             csv_data%solar(i), csv_data%rh(i)
             
        if (ios /= 0) then
            print *, "Error reading line ", i+1, " content: ", dummy_timestamp
            status = -2
            close(unit_num)
            return
        end if
        
        ! Map time as a simple hour counter (1, 2, 3...)
        csv_data%time(i) = dble(i-1) 
    end do
    
    close(unit_num)
end subroutine read_csv_data


! ============================================================
!  REFREEZING SUBROUTINE (no dt parameter - removed)
! ============================================================
subroutine refreezing_layer(T_layer, LWC_layer, ice_frac, dz_layer, Tfreeze, rho_i, rho_w, &
                            c_s, c_w, Lf, new_T, new_LWC, new_ice_frac, refrozen_mass)
    implicit none
    
    double precision, intent(in) :: T_layer, LWC_layer, ice_frac, dz_layer
    double precision, intent(in) :: Tfreeze, rho_i, rho_w, c_s, c_w, Lf
    double precision, intent(out) :: new_T, new_LWC, new_ice_frac, refrozen_mass
    
    double precision :: dT_max, dtheta_w_max, dtheta_w, dtheta_i, dT
    
    ! No refreezing if at or above freezing, or no liquid water
    if ((T_layer >= Tfreeze) .or. (LWC_layer <= 0.0D0)) then
        new_T = T_layer
        new_LWC = LWC_layer
        new_ice_frac = ice_frac
        refrozen_mass = 0.0D0
        return
    end if
    
    ! Calculate maximum refreezing based on cold content
    dT_max = T_layer - Tfreeze
    dtheta_w_max = -(dT_max * (ice_frac * rho_i * c_s + LWC_layer * rho_w * c_w)) / (rho_w * Lf)
    dtheta_w = min(LWC_layer, dtheta_w_max)
    dtheta_i = (rho_w / rho_i) * dtheta_w
    dT = (dtheta_w * rho_w * Lf) / (ice_frac * rho_i * c_s + LWC_layer * rho_w * c_w)
    
    new_T = T_layer + dT
    new_LWC = LWC_layer - dtheta_w
    new_ice_frac = ice_frac + dtheta_i
    refrozen_mass = dtheta_w * dz_layer * rho_w  ! Include layer thickness
    
end subroutine refreezing_layer

! ============================================================
!  PERCOLATION SUBROUTINE (Bucket method) - returns mass flux
! ============================================================
subroutine percolate_water(LWC_array, heights, theta_e, n_layers, rho_w, runoff)
    implicit none
    
    integer, intent(in) :: n_layers
    double precision, dimension(n_layers), intent(inout) :: LWC_array
    double precision, dimension(n_layers), intent(in) :: heights
    double precision, intent(in) :: theta_e, rho_w
    double precision, intent(out) :: runoff
    
    integer :: i
    double precision :: excess, excess_mass
    
    runoff = 0.0D0
    
    ! Percolate through layers
    do i = 1, n_layers - 1
        if (LWC_array(i) > theta_e) then
            excess = LWC_array(i) - theta_e
            LWC_array(i) = theta_e
            excess_mass = excess * heights(i)
            LWC_array(i+1) = LWC_array(i+1) + excess_mass / heights(i+1)
        end if
    end do
    
    ! Bottom layer runoff - return mass flux [kg/m²]
    if (LWC_array(n_layers) > theta_e) then
        excess = LWC_array(n_layers) - theta_e
        LWC_array(n_layers) = theta_e
        runoff = excess * heights(n_layers) * rho_w  ! [kg/m²]
    end if
    
end subroutine percolate_water

! ============================================================
!  TRIDIAGONAL MATRIX ALGORITHM (TDMA / Thomas Algorithm)
! ============================================================
subroutine solve_tdma(a, b, c, d, n, x)
    implicit none
    
    integer, intent(in) :: n
    double precision, dimension(n), intent(in) :: a, b, c, d
    double precision, dimension(n), intent(out) :: x
    
    double precision, dimension(n) :: c_prime, d_prime
    double precision :: denom
    integer :: i
    
    ! Forward sweep
    c_prime(1) = c(1) / b(1)
    d_prime(1) = d(1) / b(1)
    
    do i = 2, n
        denom = b(i) - a(i) * c_prime(i-1)
        c_prime(i) = c(i) / denom
        d_prime(i) = (d(i) - a(i) * d_prime(i-1)) / denom
    end do
    
    ! Backward substitution
    x(n) = d_prime(n)
    do i = n-1, 1, -1
        x(i) = d_prime(i) - c_prime(i) * x(i+1)
    end do
    
end subroutine solve_tdma

! ============================================================
!  MULTI-LAYER INSULATION RESISTANCE CALCULATOR
! ============================================================
subroutine compute_insulation_resistance_multilayer(T_surf, T_snow, k_eff, Hi, &
                                                    N_ins, D_ins, dt_substep, &
                                                    R_ins, T_profile)
    use enhanced_features_module
    implicit none
    
    double precision, intent(in) :: T_surf, T_snow, k_eff, Hi, D_ins, dt_substep
    integer, intent(in) :: N_ins
    double precision, intent(out) :: R_ins
    double precision, dimension(N_ins+1), intent(out) :: T_profile
    
    ! Local variables
    double precision :: dx, dFo, q_flux
    integer :: n_substeps, i, j
    double precision, dimension(N_ins+1) :: a_vec, b_vec, c_vec, d_vec, T_n
    
    dx = Hi / dble(N_ins)
    dFo = D_ins * dt_substep / (dx**2)
    n_substeps = int(3600.0D0 / dt_substep)
    
    ! Initialize on first call
    if (.not. T_ins_initialized) then
        if (allocated(T_ins_stored)) deallocate(T_ins_stored)
        allocate(T_ins_stored(N_ins+1))
        do i = 1, N_ins+1
            T_ins_stored(i) = T_surf + (T_snow - T_surf) * dble(i-1) / dble(N_ins)
        end do
        T_ins_initialized = .true.
    end if
    
    T_n = T_ins_stored
    
    ! Time integration loop
    do i = 1, n_substeps
        ! Build tridiagonal system
        a_vec = 0.0D0
        b_vec = 0.0D0
        c_vec = 0.0D0
        d_vec = 0.0D0
        
        ! Boundary condition at outer surface (node 1)
        b_vec(1) = 1.0D0
        c_vec(1) = 0.0D0
        d_vec(1) = T_surf
        
        ! Interior nodes
        do j = 2, N_ins
            a_vec(j) = -dFo
            b_vec(j) = 1.0D0 + 2.0D0 * dFo
            c_vec(j) = -dFo
            d_vec(j) = T_n(j)
        end do
        
        ! Boundary condition at inner surface (node N_ins+1)
        a_vec(N_ins+1) = 0.0D0
        b_vec(N_ins+1) = 1.0D0
        c_vec(N_ins+1) = 0.0D0
        d_vec(N_ins+1) = T_snow
        
        ! Solve tridiagonal system
        call solve_tdma(a_vec, b_vec, c_vec, d_vec, N_ins+1, T_n)
    end do
    
    ! Store for next call
    T_ins_stored = T_n
    T_profile = T_n
    
    ! Compute effective resistance from heat flux
    q_flux = k_eff * (T_n(1) - T_n(2)) / dx
    
    if (abs(T_surf - T_snow) > 0.01D0) then
        if (abs(q_flux) > 1.0D-6) then
            R_ins = (T_surf - T_snow) / q_flux
        else
            R_ins = Hi / k_eff
        end if
    else
        R_ins = Hi / k_eff
    end if
    
end subroutine compute_insulation_resistance_multilayer

! ============================================================
!  MODIFIED INSULATION STEP WITH MULTI-LAYER OPTION
! ============================================================
subroutine insulation_step_enhanced(state_in, forc, p, delta_t, USE_MULTILAYER, &
                                   N_ins, D_ins, T_snow_surf, &
                                   R_ins, q_solar, q_rain_snow, q_evap, state_out)
    use enhanced_features_module
    implicit none
    
    ! Inputs
    type(InsulationState), intent(inout) :: state_in
    type(Forcing), intent(in) :: forc
    type(InsulationParameters), intent(in) :: p
    double precision, intent(in) :: delta_t, T_snow_surf, D_ins
    logical, intent(in) :: USE_MULTILAYER
    integer, intent(in) :: N_ins
    
    ! Outputs
    double precision, intent(out) :: R_ins, q_solar, q_rain_snow, q_evap
    type(InsulationState), intent(out) :: state_out
    
    ! Local variables
    double precision :: W, age_days, f, age_yr
    double precision :: k_moist, k_age_factor, k_eff
    double precision :: alpha_moist, alpha_age, alpha_eff
    double precision :: eta_rain, P_in_mass, zeta_rain
    double precision :: Tc_s, Tc_a, e_sat_s, e_sat_a, e_surf, e_air, VPD
    double precision :: E0, f_breath, E, D, dWdt, W_new, age_days_new
    double precision :: T_surf_approx
    double precision, dimension(N_ins+1) :: T_profile
    
    ! Extract state variables
    W = state_in%W
    age_days = state_in%age_days
    
    ! 1) Moisture & age factors
    f = max(0.0D0, min(1.0D0, W / p%W_sat))
    age_yr = age_days / 365.0D0
    
    ! Conductivity
    k_moist = p%k_dry + (p%k_sat - p%k_dry) * (f**p%n_k)
    k_age_factor = 1.0D0 + p%delta_k_age * (1.0D0 - exp(-age_yr / p%tau_k_years))
    k_eff = k_moist * k_age_factor
    
    ! Compute resistance (multi-layer or simple)
    if (USE_MULTILAYER) then
        T_surf_approx = forc%Ta  ! Approximation
        call compute_insulation_resistance_multilayer(T_surf_approx, T_snow_surf, &
                                                     k_eff, p%Hi, N_ins, D_ins, 10.0D0, &
                                                     R_ins, T_profile)
    else
        R_ins = p%Hi / k_eff
    end if
    
    ! Absorptivity
    alpha_moist = p%alpha_dry + (p%alpha_wet - p%alpha_dry) * (f**p%n_alpha)
    alpha_age = alpha_moist + p%delta_alpha_age * (1.0D0 - exp(-age_yr / p%tau_alpha_years))
    alpha_eff = max(0.0D0, min(1.0D0, alpha_age))
    
    ! 2) Solar absorption
    q_solar = alpha_eff * forc%Isolar
    
    ! 3) Rain infiltration & heat
    eta_rain = max(0.0D0, 1.0D0 - f)
    P_in_mass = eta_rain * p%rho_w * forc%Prain
    
    zeta_rain = p%zeta0 * exp(-p%gamma_H * p%Hi) * exp(-p%gamma_W * f)
    q_rain_snow = zeta_rain * p%rho_w * p%c_w * forc%Prain * (forc%T_rain - p%Tfreeze)
    
    ! 4) Evaporation
    Tc_s = p%Tfreeze - 273.15D0
    Tc_a = forc%Ta - 273.15D0
    
    e_sat_s = 611.0D0 * exp(17.27D0 * Tc_s / (Tc_s + 237.3D0))
    e_sat_a = 611.0D0 * exp(17.27D0 * Tc_a / (Tc_a + 237.3D0))
    
    e_surf = e_sat_s
    e_air = forc%RH * e_sat_a
    
    VPD = max(0.0D0, e_surf - e_air)
    
    E0 = p%rho_air * p%C_E * p%U10 * VPD / p%P0
    f_breath = exp(-p%beta_w * f)
    E = E0 * f_breath
    
    q_evap = -p%Lv * E
    
    ! 5) Drainage & moisture update
    D = p%K_D * max(0.0D0, W - p%W_field)
    dWdt = P_in_mass - E - D
    W_new = W + delta_t * dWdt
    W_new = max(0.0D0, min(p%W_sat, W_new))
    
    age_days_new = age_days + delta_t / 86400.0D0
    
    ! Update output state
    state_out = state_in
    state_out%W = W_new
    state_out%age_days = age_days_new
    state_out%k_eff = k_eff
    state_out%alpha_eff = alpha_eff
    state_out%f_sat = f
    
end subroutine insulation_step_enhanced

program main_enhanced
    use enhanced_features_module
    implicit none
    
    type(InsulationParameters) :: InsPar
    type(InsulationState)      :: InsState
    type(CSVData)              :: met_data
    type(Forcing)              :: forc
    
    ! EXTERNAL function declarations
    double precision, external :: ground_flux
    
    ! Control flags
    logical :: USE_ADVANCED_INSULATION, USE_REFREEZING, USE_PERCOLATION
    logical :: USE_MULTILAYER_INSULATION, USE_CSV_DATA
    
    ! Variable declarations
    integer :: Nt, k, status, i, n_hours
    double precision :: sigma, Lf, rho_i, rho_s, c_s, rho_w, c_w, Tfreeze
    double precision :: Hs, Ns, dz_s, k_snow, Hi, k_i_base, alpha_const
    double precision :: eta_rain_const, Hg_ins, kg_ins, h_conv, epsilon
    double precision :: T_mean, h_rad, h_eff, R_eff, R_layer, R_12, R_23
    double precision :: R_g_ins, R_3g, Tg, Cs_layer, T1_init, T2_init, T3_init
    double precision, dimension(3) :: Tempe, R_nm, T_new
    double precision, dimension(3) :: LWC, ice_fractions, layer_heights
    double precision :: theta_e
    
    ! Multi-layer insulation parameters
    integer :: N_ins
    double precision :: dz_ins, rho_dry, moist_cont, rho_wet, c_dry, c_wet, D_ins
    
    ! Time integration
    double precision, allocatable, dimension(:) :: t_vec
    double precision :: t0, tf, dt, t, t_mid, dt_data
    
    ! History arrays
    double precision, dimension(:,:), allocatable :: T_hist, LWC_hist
    double precision, dimension(:), allocatable :: qnet_surf_hist, qa_hist, qsolar_hist
    double precision, dimension(:), allocatable :: qrain_hist, qevap_hist, qground_hist
    double precision, dimension(:), allocatable :: Ta_hist, Isolar_hist, Prain_hist
    double precision, dimension(:), allocatable :: melt_rate_hist, refrozen_hist, runoff_hist
    double precision, dimension(:), allocatable :: W_hist, k_eff_hist, alpha_hist
    double precision, dimension(:), allocatable :: Rins_hist, fsat_hist
    
    ! Flux and energy variables
    double precision :: R_ins, q_solar_ins, q_rain_snow, q_evap_val
    double precision :: R_a2s, q_a_mid, q_surf_mid, q_ground_mid
    double precision :: E_melt, dE_melt, dM_melt, surface_melt_water
    double precision :: total_refrozen, refrozen_mass, runoff
    double precision :: Ta_C, Ta_K, Isolar, Prain_mh, Prain, wind_speed, RH_pct, RH_frac
    
    ! Progress tracking
    integer :: pct
    integer, dimension(4) :: progress_points
    
    ! Solver choice
    INTEGER :: SOLVER_CHOICE
    
    ! --- NEW: Weather Data for Common Block ---
    INTEGER, PARAMETER :: MAX_WEATHER = 100000 
    DOUBLE PRECISION :: temp_vec_com(MAX_WEATHER)
    DOUBLE PRECISION :: dt_data_com
    INTEGER          :: itemp_size_com

    ! --- UPDATED: Common block for LSODA ---
    DOUBLE PRECISION :: R_eff_com, R_ins_com, q_solar_com, q_rain_com, Tg_com
    DOUBLE PRECISION :: q_evap_com, R_12_com, R_23_com, R_3g_com, Cs_layer_com
    
    ! The order in this list must match the RHS and the Integrator exactly!
    COMMON /SNOW_PARAMS/ R_eff_com, R_ins_com, q_solar_com, q_rain_com, q_evap_com, &
                        R_12_com, R_23_com, R_3g_com, Cs_layer_com, Tg_com, &
                        temp_vec_com, dt_data_com, itemp_size_com
    
    ! ============================================================
    !  Enhanced Multi-layer Snow Storage RC Model with Real Data
    ! ============================================================
    
    print *, "============================================================"
    print *, "Enhanced Snow Storage RC Model with Real Data"
    print *, "============================================================"
    
    ! Control flags
    USE_ADVANCED_INSULATION = .TRUE.
    USE_REFREEZING = .TRUE.
    USE_PERCOLATION = .TRUE.
    USE_MULTILAYER_INSULATION = .TRUE.
    USE_CSV_DATA = .TRUE.
    
    ! Solver choice
    SOLVER_CHOICE = 2  ! 1=RK4, 2=DOPRI5, 3=LSODA
    
    ! ============================================================
    ! Physical constants
    ! ============================================================
    sigma   = 5.670374419D-8
    Lf      = 3.34D5
    rho_i   = 917.0D0
    rho_s   = 400.0D0
    c_s     = 2100.0D0
    rho_w   = 1000.0D0
    c_w     = 4180.0D0
    Tfreeze = 273.15D0
    
    ! ============================================================
    ! Snow & insulation geometry
    ! ============================================================
    Hs   = 2.0D0
    Ns   = 3.0D0
    dz_s = Hs / Ns
    
    k_snow = 0.25D0
    Hi     = 0.10D0  ! Matching Python: 0.1m insulation
    
    ! Multi-layer insulation setup
    if (USE_MULTILAYER_INSULATION) then
        N_ins = 20
        dz_ins = Hi / dble(N_ins)
    else
        N_ins = 1
        dz_ins = Hi
    end if
    
    ! Insulation material properties
    k_i_base   = 0.32D0
    rho_dry    = 100.0D0
    moist_cont = 50.0D0
    rho_wet    = rho_dry + moist_cont/100.0D0*1000.0D0
    c_dry      = 0.99D3
    c_wet      = (1.0D0 - moist_cont/100.0D0)*c_dry + moist_cont/100.0D0*c_w
    D_ins      = k_i_base / (c_wet * rho_wet)
    
    ! Simple insulation parameters
    alpha_const    = 0.80D0
    eta_rain_const = 1.0D0
    
    ! Ground insulation
    Hg_ins = 0.3D0
    kg_ins = 0.04D0
    
    ! ============================================================
    ! Surface heat transfer
    ! ============================================================
    h_conv  = 8.0D0
    epsilon = 0.95D0
    T_mean  = 273.15D0 + 3.0D0
    h_rad   = 4.0D0 * epsilon * sigma * T_mean**3
    h_eff   = h_conv + h_rad
    R_eff   = 1.0D0 / h_eff
    
    ! ============================================================
    ! Thermal resistances
    ! ============================================================
    R_layer = dz_s / k_snow
    R_12    = R_layer
    R_23    = R_layer
    R_g_ins = Hg_ins / kg_ins
    R_3g    = R_layer + R_g_ins
    R_nm    = [R_12, R_23, R_3g]
    
    ! ============================================================
    ! Ground & snow capacity
    ! ============================================================
    Tg       = 273.15D0 + 2.0D0
    Cs_layer = rho_s * c_s * dz_s
    
    ! ============================================================
    ! Initial conditions
    ! ============================================================
    T1_init = 273.15D0 - 2.0D0
    T2_init = 273.15D0 - 4.0D0
    T3_init = 273.15D0 - 6.0D0
    Tempe = [T1_init, T2_init, T3_init]
    
    ! Layer properties for refreezing/percolation
    LWC = [0.0D0, 0.0D0, 0.0D0]
    theta_e = 0.03D0
    ice_fractions = [0.4D0, 0.4D0, 0.4D0]
    layer_heights = [dz_s, dz_s, dz_s]
    
    ! ============================================================
    ! Advanced insulation parameters
    ! ============================================================
    if (USE_ADVANCED_INSULATION) then
        InsPar = InsulationParameters( &
            Hi = Hi, &
            k_dry = 0.06D0, &
            k_sat = 0.30D0, &
            n_k = 1.5D0, &
            W_sat = 30.0D0, &
            W_field = 10.0D0, &
            alpha_dry = 0.10D0, &
            alpha_wet = 0.25D0, &
            n_alpha = 1.0D0, &
            delta_k_age = 0.5D0, &
            tau_k_years = 2.0D0, &
            delta_alpha_age = 0.05D0, &
            tau_alpha_years = 2.0D0, &
            zeta0 = 0.3D0, &
            gamma_H = 0.5D0, &
            gamma_W = 2.0D0, &
            beta_w = 3.0D0, &
            K_E = 1.0D-5, &
            K_D = 5.0D-6, &
            Lv = 2.5D6, &
            rho_w = rho_w, &
            c_w = c_w, &
            Tfreeze = Tfreeze, &
            rho_air = 1.2D0, &
            C_E = 1.3D-3, &
            U10 = 2.0D0, &
            P0 = 101325.0D0 &
        )
        
        InsState = InsulationState(&
            W = 5.0D0, &
            age_days = 0.0D0, &
            k_eff = 0.D0, &
            alpha_eff = 0.D0, &
            f_sat = 0.D0 &
        )
    else
        InsState = InsulationState( &
            W = 0.0D0, &
            age_days = 0.0D0, &
            k_eff = 0.D0, &
            alpha_eff = 0.D0, &
            f_sat = 0.D0 &
        )
    end if
    
    ! ============================================================
    ! Load meteorological data
    ! ============================================================
    print *, ""
    !if (USE_CSV_DATA) then
        print *, "Loading meteorological data from DATA.csv..."
        call read_csv_data('DATA.csv', met_data, status)
        if (status /= 0) then
            print *, "FATAL ERROR: Could not read DATA.csv. File may be missing or formatted incorrectly."
            stop
        end if
        
    !    if (status /= 0) then
    !        print *, "ERROR: Could not read DATA.csv (status = ", status, ")"
    !        print *, "Falling back to synthetic forcing functions."
    !        USE_CSV_DATA = .FALSE.
    !    else
    !        print '(A, I0, A)', " Loaded ", met_data%n_points, " hourly data points"
    !        if (allocated(met_data%time)) then
    !            print '(A, F0.1, A, F0.1)', " Period: ", met_data%time(1), &
    !                                       " to ", met_data%time(met_data%n_points)
    !        end if
    !    end if
    !end if
    
    ! ============================================================
    ! Time integration settings
    ! ============================================================
    t0 = 0.0D0
    dt = 600.0D0      ! 10 minutes
    dt_data = 3600.0D0  ! CSV data timestep (1 hour)
    
    !if (USE_CSV_DATA) then
        ! Simulate for length of available data
    n_hours = met_data%n_points
    tf = dble(n_hours) * dt_data
    !else
    !    ! Default simulation length
    !    tf = 120.0D0 * 24.0D0 * 3600.0D0  ! 120 days
    !    n_hours = int(tf / dt_data)
    !end if
    
    ! Use arange function directly (no pre-allocation needed)
    ! Remove the manual 'allocate(t_vec(Nt))' line entirely
    t_vec = arange(t0, tf, dt)
    Nt = size(t_vec)
    
    print *, ""
    print *, "Simulation settings:"
    print '(A, I0, A, F6.1, A)', "  Duration: ", n_hours, " hours (", &
                                  dble(n_hours)/24.0D0, " days)"
    print '(A, F8.1, A, F5.1, A)', "  Time step: ", dt, " s (", dt/60.0D0, " min)"
    print '(A, I0)', "  Total steps: ", Nt
    print '(A, L, A, I0, A)', "  Multi-layer insulation: ", USE_MULTILAYER_INSULATION, &
                              " (", N_ins, " layers)"
    print '(A, L)', "  Refreezing: ", USE_REFREEZING
    print '(A, L)', "  Percolation: ", USE_PERCOLATION
    
    SELECT CASE (SOLVER_CHOICE)
    CASE (1)
        print *, " Solver: RK4"
    CASE (2)
        print *, " Solver: DOPRI5"
    CASE (3)
        print *, " Solver: LSODA"
    END SELECT
    
    ! ============================================================
    ! Allocate history arrays
    ! ============================================================
    allocate(T_hist(Nt, 3), LWC_hist(Nt, 3))
    T_hist = 0.0D0
    LWC_hist = 0.0D0
    T_hist(1, :) = Tempe
    LWC_hist(1, :) = LWC
    
    allocate(qnet_surf_hist(Nt), qa_hist(Nt), qsolar_hist(Nt))
    allocate(qrain_hist(Nt), qevap_hist(Nt), qground_hist(Nt))
    qnet_surf_hist = 0.0D0
    qa_hist = 0.0D0
    qsolar_hist = 0.0D0
    qrain_hist = 0.0D0
    qevap_hist = 0.0D0
    qground_hist = 0.0D0
    
    allocate(Ta_hist(Nt), Isolar_hist(Nt), Prain_hist(Nt))
    Ta_hist = 0.0D0
    Isolar_hist = 0.0D0
    Prain_hist = 0.0D0
    
    allocate(melt_rate_hist(Nt), refrozen_hist(Nt), runoff_hist(Nt))
    melt_rate_hist = 0.0D0
    refrozen_hist = 0.0D0
    runoff_hist = 0.0D0
    E_melt = 0.0D0
    
    allocate(W_hist(Nt), k_eff_hist(Nt), alpha_hist(Nt))
    allocate(Rins_hist(Nt), fsat_hist(Nt))
    W_hist = 0.0D0
    k_eff_hist = 0.0D0
    alpha_hist = 0.0D0
    Rins_hist = 0.0D0
    fsat_hist = 0.0D0
    
    ! ============================================================
    ! Main time loop
    ! ============================================================
    print *, ""
    print *, "Running simulation..."
    
    ! Progress tracking
    progress_points = [int(Nt * 0.25), int(Nt * 0.5), int(Nt * 0.75), Nt]
    
    do k = 1, Nt-1
        t = t_vec(k)
        t_mid = t + dt / 2.0D0
        
        ! Show progress at specific intervals
        ! We use max(1, ...) to prevent division by zero for very short runs
        if (mod(k, max(1, Nt/4)) == 0 .or. k == Nt-1) then
            ! dble(k)/dble(Nt) ensures floating point division
            ! nint() rounds 24.999 to 25
            pct = nint(100.0D0 * dble(k) / dble(Nt))
            
            ! Force the final step to show 100%
            if (k == Nt-1) pct = 100
            
            print '(A, I3, A)', "  Progress: ", pct, "%"
        end if
        
        ! ========================================================
        ! Get forcing data (CSV or synthetic)
        ! ========================================================
        !if (USE_CSV_DATA) then
            ! Interpolate from CSV data
            Ta_C = interpolate_data(met_data%temp, t_mid, dt_data)
            Isolar = interpolate_data(met_data%solar, t_mid, dt_data)
            Prain_mh = interpolate_data(met_data%precip, t_mid, dt_data)
            wind_speed = interpolate_data(met_data%wind, t_mid, dt_data)
            RH_pct = interpolate_data(met_data%rh, t_mid, dt_data)
            
            ! Convert units
            Ta_K = Ta_C + 273.15D0
            Prain = Prain_mh / 3600.0D0  ! [m/h] -> [m/s]
            RH_frac = RH_pct / 100.0D0
            
            ! Build forcing structure
            forc%Isolar = Isolar
            forc%Prain = Prain
            forc%T_rain = Ta_K  ! Assume rain at air temperature
            forc%RH = RH_frac
            forc%Ta = Ta_K
        !else
        !    ! Use synthetic forcing functions
        !    call set_forcing(t_mid, forc)
        !end if
        
        ! ========================================================
        ! Update insulation properties
        ! ========================================================
        if (USE_ADVANCED_INSULATION) then
            call insulation_step_enhanced(InsState, forc, InsPar, dt, &
                                         USE_MULTILAYER_INSULATION, &
                                         N_ins, D_ins, Tempe(1), &
                                         R_ins, q_solar_ins, q_rain_snow, q_evap_val, &
                                         InsState)
            
            W_hist(k) = InsState%W
            k_eff_hist(k) = InsState%k_eff
            alpha_hist(k) = InsState%alpha_eff
            Rins_hist(k) = R_ins
            fsat_hist(k) = InsState%f_sat
        else
            R_ins = Hi / k_i_base
            q_solar_ins = alpha_const * forc%Isolar
            q_rain_snow = eta_rain_const * rho_w * c_w * &
                         forc%Prain * (forc%T_rain - Tfreeze)
            q_evap_val = 0.0D0
        end if
        
        R_a2s = R_eff + R_ins
        
        ! ========================================================
        ! Integrate temperature (choose solver)
        ! ========================================================
        SELECT CASE (SOLVER_CHOICE)
        CASE (1)
            ! RK4
            call step_RK4(t, dt, Tempe, R_eff, R_ins, q_solar_ins, q_rain_snow, &
                         q_evap_val, met_data%temp, n_hours, dt_data, & 
                         R_nm, Cs_layer, Tg, T_new)
        CASE (2)
            ! DOPRI5
            call integrate_DOPRI5(t, dt, Tempe, R_ins, R_eff, q_solar_ins, &
                                q_rain_snow, q_evap_val, R_nm, Cs_layer, &
                                Tg, met_data%temp, n_hours, dt_data, T_new)
        CASE (3)
            ! LSODA - set common block
            R_eff_com = R_eff
            R_ins_com = R_ins
            q_solar_com = q_solar_ins
            q_rain_com = q_rain_snow
            q_evap_com = q_evap_val
            R_12_com = R_nm(1)
            R_23_com = R_nm(2)
            R_3g_com = R_nm(3)
            Cs_layer_com = Cs_layer
            Tg_com = Tg
            
            call integrate_LSODA(t, dt, Tempe, R_ins, R_eff, q_solar_ins, &
                                q_rain_snow, q_evap_val, R_nm, Cs_layer, &
                                Tg, met_data%temp, n_hours, dt_data, T_new)
        END SELECT
        
        ! ========================================================
        ! Calculate fluxes
        ! ========================================================
        q_a_mid = (forc%Ta - Tempe(1)) / R_a2s
        q_surf_mid = q_a_mid + q_solar_ins + q_rain_snow + q_evap_val
        q_ground_mid = ground_flux(Tempe(3), Tg, R_3g)
        
        ! ========================================================
        ! Refreezing in each layer
        ! ========================================================
        total_refrozen = 0.0D0
        if (USE_REFREEZING) then
            do i = 1, 3
                call refreezing_layer(T_new(i), LWC(i), ice_fractions(i), layer_heights(i), &
                                     Tfreeze, rho_i, rho_w, c_s, c_w, Lf, &
                                     T_new(i), LWC(i), ice_fractions(i), refrozen_mass)
                total_refrozen = total_refrozen + refrozen_mass
            end do
        end if
        refrozen_hist(k) = total_refrozen
        
        ! ========================================================
        ! Surface melting
        ! ========================================================
        surface_melt_water = 0.0D0
        if (T_new(1) > Tfreeze) then
            T_new(1) = Tfreeze
            if (q_surf_mid > 0.0D0) then
                dE_melt = q_surf_mid * dt
                dM_melt = dE_melt / (rho_i * Lf)
                surface_melt_water = dM_melt  ! [m w.e.]
                LWC(1) = LWC(1) + surface_melt_water / layer_heights(1)
            else
                dE_melt = 0.0D0
                dM_melt = 0.0D0
            end if
            E_melt = E_melt + dE_melt
            melt_rate_hist(k) = dM_melt / dt
        else
            melt_rate_hist(k) = 0.0D0
        end if
        
        ! ========================================================
        ! Percolation
        ! ========================================================
        runoff = 0.0D0
        if (USE_PERCOLATION) then
            call percolate_water(LWC, layer_heights, theta_e, 3, rho_w, runoff)
        end if
        runoff_hist(k) = runoff
        
        ! ========================================================
        ! Update and store
        ! ========================================================
        Tempe = T_new
        T_hist(k+1, :) = Tempe
        LWC_hist(k+1, :) = LWC
        
        qnet_surf_hist(k) = q_surf_mid
        qa_hist(k) = q_a_mid
        qsolar_hist(k) = q_solar_ins
        qrain_hist(k) = q_rain_snow
        qevap_hist(k) = q_evap_val
        qground_hist(k) = q_ground_mid
        
        Ta_hist(k) = forc%Ta
        Isolar_hist(k) = forc%Isolar
        Prain_hist(k) = forc%Prain
    end do
    
    ! Final values
    !if (USE_CSV_DATA) then
        Ta_hist(Nt) = interpolate_data(met_data%temp, t_vec(Nt), dt_data) + 273.15D0
        Isolar_hist(Nt) = interpolate_data(met_data%solar, t_vec(Nt), dt_data)
        Prain_hist(Nt) = interpolate_data(met_data%precip, t_vec(Nt), dt_data) / 3600.0D0
    !else
    !    call set_forcing(t_vec(Nt), forc)
    !    Ta_hist(Nt) = forc%Ta
    !    Isolar_hist(Nt) = forc%Isolar
    !    Prain_hist(Nt) = forc%Prain
    !end if
    
    if (USE_ADVANCED_INSULATION) then
        W_hist(Nt) = InsState%W
        k_eff_hist(Nt) = InsState%k_eff
        alpha_hist(Nt) = InsState%alpha_eff
        Rins_hist(Nt) = R_ins
        fsat_hist(Nt) = InsState%f_sat
    end if
    
    !print *, "  Progress: 100%"
    print *, ""
    print *, "Simulation complete!"
    
    ! ============================================================
    ! Enhanced energy diagnostics
    ! ============================================================
    call compute_energy_diagnostics_enhanced(t_vec, Nt, qa_hist, qsolar_hist, &
                                            qrain_hist, qevap_hist, qground_hist, &
                                            T_hist, Cs_layer, E_melt, rho_i, Lf, &
                                            refrozen_hist, runoff_hist, n_hours, &
                                            USE_ADVANCED_INSULATION, USE_MULTILAYER_INSULATION, &
                                            USE_REFREEZING, USE_PERCOLATION, &
                                            N_ins, SOLVER_CHOICE)
    
    ! ============================================================
    ! Cleanup
    ! ============================================================
    !if (USE_CSV_DATA) then
    !    if (allocated(met_data%time))   deallocate(met_data%time)
    !    if (allocated(met_data%temp))   deallocate(met_data%temp)
    !    if (allocated(met_data%wind))   deallocate(met_data%wind)
    !    if (allocated(met_data%precip)) deallocate(met_data%precip)
    !    if (allocated(met_data%solar))  deallocate(met_data%solar)
    !    if (allocated(met_data%rh))     deallocate(met_data%rh)
    !end if
    
    deallocate(t_vec)
    deallocate(T_hist, LWC_hist)
    deallocate(qnet_surf_hist, qa_hist, qsolar_hist)
    deallocate(qrain_hist, qevap_hist, qground_hist)
    deallocate(Ta_hist, Isolar_hist, Prain_hist)
    deallocate(melt_rate_hist, refrozen_hist, runoff_hist)
    deallocate(W_hist, k_eff_hist, alpha_hist)
    deallocate(Rins_hist, fsat_hist)
    deallocate(met_data%time)
    deallocate(met_data%temp)
    deallocate(met_data%wind)
    deallocate(met_data%precip)
    deallocate(met_data%solar)
    deallocate(met_data%rh)

contains

    ! ============================================================
    !  NUMPY STYLE ARANGE FUNCTION
    ! ============================================================
    function arange(start, stop, step) result(arr)
        implicit none
        double precision, intent(in) :: start, stop, step
        double precision, allocatable :: arr(:)
        integer :: n, m
        
        ! Match Python's np.arange behavior: [start, stop)
        ! Does NOT include stop
        if (abs(step) < tiny(1.0D0)) then
            n = 1
        else
            ! Python: n = ceil((stop - start) / step) but excludes endpoint
            n = floor((stop - start) / step)
            if (n < 0) n = 0
        end if
        
        n = max(0, n)
        if (n == 0) then
            allocate(arr(0))
            return
        end if
        
        allocate(arr(n))
        
        do m = 1, n
            arr(m) = start + (m-1)*step
        end do
    end function arange

end program main_enhanced

! ============================================================
!  Integration routines for discrete data
! ============================================================

FUNCTION trapz_discrete(y, x, n) RESULT(integral)
    !======================================================================
    ! Trapezoidal rule for discrete data arrays
    ! y: array of function values
    ! x: array of corresponding x values
    ! n: size of arrays
    !======================================================================
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    DOUBLE PRECISION, DIMENSION(n), INTENT(IN) :: y, x
    DOUBLE PRECISION :: integral
    INTEGER :: i
    
    integral = 0.0D0
    do i = 1, n-1
        integral = integral + 0.5D0 * (y(i+1) + y(i)) * (x(i+1) - x(i))
    end do
END FUNCTION trapz_discrete

FUNCTION simpson_discrete(y, x, n) RESULT(integral)
    !======================================================================
    ! Simpson's 1/3 rule for discrete data arrays (composite)
    ! Handles both even and odd number of points
    ! y: array of function values
    ! x: array of corresponding x values  
    ! n: size of arrays
    !======================================================================
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    DOUBLE PRECISION, DIMENSION(n), INTENT(IN) :: y, x
    DOUBLE PRECISION :: integral
    INTEGER :: i
    DOUBLE PRECISION :: h, sum_odd, sum_even
    
    ! External functions
    DOUBLE PRECISION, EXTERNAL :: trapz_discrete
    
    if (n < 3) then
        ! Fall back to trapezoidal for too few points
        integral = trapz_discrete(y, x, n)
        return
    end if
    
    ! Check if spacing is uniform (required for standard Simpson's rule)
    h = x(2) - x(1)
    
    ! For uniform spacing and odd number of points (even number of intervals)
    if (mod(n, 2) == 1) then
        sum_odd = 0.0D0
        sum_even = 0.0D0
        
        do i = 2, n-1, 2
            sum_odd = sum_odd + y(i)
        end do
        
        do i = 3, n-2, 2
            sum_even = sum_even + y(i)
        end do
        
        integral = (h / 3.0D0) * (y(1) + 4.0D0*sum_odd + 2.0D0*sum_even + y(n))
    else
        ! For even number of points: use Simpson's on n-1 points + trapezoidal on last interval
        sum_odd = 0.0D0
        sum_even = 0.0D0
        
        do i = 2, n-2, 2
            sum_odd = sum_odd + y(i)
        end do
        
        do i = 3, n-3, 2
            sum_even = sum_even + y(i)
        end do
        
        integral = (h / 3.0D0) * (y(1) + 4.0D0*sum_odd + 2.0D0*sum_even + y(n-1))
        integral = integral + 0.5D0 * h * (y(n-1) + y(n))
    end if
END FUNCTION simpson_discrete

! ============================================================
!  Energy diagnostics subroutine
! ============================================================
subroutine compute_energy_diagnostics_enhanced(t_vec, Nt, qa_hist, qsolar_hist, &
                                              qrain_hist, qevap_hist, qground_hist, &
                                              T_hist, Cs_layer, E_melt, rho_i, Lf, &
                                              refrozen_hist, runoff_hist, n_hours, &
                                              USE_ADVANCED, USE_MULTILAYER, &
                                              USE_REFREEZE, USE_PERCOL, &
                                              N_ins, SOLVER_CHOICE)
    IMPLICIT NONE
    ! Arguments
    INTEGER, INTENT(IN) :: Nt, SOLVER_CHOICE, n_hours, N_ins
    DOUBLE PRECISION, DIMENSION(Nt), INTENT(IN) :: t_vec, qa_hist, qsolar_hist
    DOUBLE PRECISION, DIMENSION(Nt), INTENT(IN) :: refrozen_hist, runoff_hist
    DOUBLE PRECISION, DIMENSION(Nt), INTENT(IN) :: qrain_hist, qevap_hist, qground_hist
    DOUBLE PRECISION, DIMENSION(Nt, 3), INTENT(IN) :: T_hist
    DOUBLE PRECISION, INTENT(IN) :: Cs_layer, E_melt, rho_i, Lf
    LOGICAL, INTENT(IN) :: USE_ADVANCED, USE_MULTILAYER, USE_REFREEZE, USE_PERCOL
    
    ! Local Variables
    DOUBLE PRECISION :: E_a, E_solar, E_rain, E_evap, E_g, E_refrozen
    DOUBLE PRECISION :: E_total_in, E_snow_change, E_balance, M_melt, M_runoff
    
    ! External function for integration
    DOUBLE PRECISION, EXTERNAL :: simpson_discrete
    
    ! 1. Integrate energy fluxes using Simpson's rule [J/m2]
    E_a     = simpson_discrete(qa_hist, t_vec, Nt)
    E_solar = simpson_discrete(qsolar_hist, t_vec, Nt)
    E_rain  = simpson_discrete(qrain_hist, t_vec, Nt)
    E_evap  = simpson_discrete(qevap_hist, t_vec, Nt)
    E_g     = simpson_discrete(qground_hist, t_vec, Nt)
    
    E_total_in = E_a + E_solar + E_rain + E_evap + E_g
    
    ! 2. Energy storage and Phase changes
    ! Change in internal energy: Cs * (T_final_sum - T_initial_sum)
    E_snow_change = Cs_layer * (sum(T_hist(Nt, :)) - sum(T_hist(1, :)))
    
    ! Refreezing energy: sum of refrozen mass [kg/m2] * Latent heat [J/kg]
    E_refrozen = sum(refrozen_hist) * Lf
    
    ! Energy Balance: In - (Storage + Melt + Refreeze)
    E_balance = E_total_in - (E_snow_change + E_melt + E_refrozen)
    
    ! 3. Mass Balance
    M_melt = E_melt / (rho_i * Lf)
    M_runoff = sum(runoff_hist)
    
    ! --- Print Results (Formatted to match Python output) ---
    print *, ""
    print *, "============================================================"
    print *, "Energy Balance Diagnostics"
    print *, "============================================================"

    SELECT CASE (SOLVER_CHOICE)
    CASE (1)
        print *, "  Solver: RK4"
    CASE (2)
        print *, "  Solver: DOPRI5"
    CASE (3)
        print *, "  Solver: LSODA"
    END SELECT

    print '(A, L1)', "  Advanced insulation:    ", USE_ADVANCED
    print '(A, L1, A, I0, A)', "  Multi-layer insulation: ", USE_MULTILAYER, " (", N_ins, " layers)"
    print '(A, L1)', "  Refreezing:             ", USE_REFREEZE
    print '(A, L1)', "  Percolation:            ", USE_PERCOL

    print *, ""
    print *, "Energy fluxes [MJ/m²]:"
    print '(A, F10.3)', "  Air convection:     ", E_a / 1.0D6
    print '(A, F10.3)', "  Solar radiation:    ", E_solar / 1.0D6
    print '(A, F10.3)', "  Rain heat:          ", E_rain / 1.0D6
    print '(A, F10.3)', "  Evaporation:        ", E_evap / 1.0D6
    print '(A, F10.3)', "  Ground heat:        ", E_g / 1.0D6
    print *, "  -----------------------------"
    print '(A, F10.3)', "  Total input:        ", E_total_in / 1.0D6

    print *, ""
    print *, "Energy storage/losses [MJ/m²]:"
    print '(A, F10.3)', "  Snow temperature:   ", E_snow_change / 1.0D6
    print '(A, F10.3)', "  Melting:            ", E_melt / 1.0D6
    print '(A, F10.3)', "  Refreezing:         ", E_refrozen / 1.0D6
    print *, "  -----------------------------"
    print '(A, F10.3)', "  Energy residual:    ", E_balance / 1.0D6

    print *, ""
    print *, "Mass balance:"
    print '(A, F10.3, A)', "  Total melt:         ", M_melt, " m w.e."
    print '(A, F10.3, A)', "  Total runoff:       ", M_runoff, " kg/m²"
    ! Calculate melt rate in mm/day (M_melt * 1000 to get mm)
    print '(A, F10.3, A)', "  Melt rate (avg):    ", (M_melt / (real(n_hours)/24.0D0)) * 1000.0D0, " mm/day"
    print *, "============================================================"

END SUBROUTINE compute_energy_diagnostics_enhanced

! ============================================================
!  Ground flux
! ============================================================
double precision function ground_flux(T3, Tg, R_3g) result(flux)
    ! Ground -> bottom snow layer flux [W/m^2].
    implicit none
    double precision, intent(in) :: T3, Tg, R_3g

    flux = (Tg - T3) / R_3g

end function ground_flux

!subroutine set_forcing(t_mid, forc)
!    use enhanced_features_module
!    implicit none
!    
!    double precision :: t_mid
!    type(Forcing) :: forc
!    double precision, external :: Ta_fun, Isolar_fun, Prain_fun
!    
!    forc%Isolar = Isolar_fun(t_mid)
!    forc%Prain  = Prain_fun(t_mid)
!    forc%T_rain = 273.15D0 + 3.0D0
!    forc%RH     = 0.80D0
!    forc%Ta     = Ta_fun(t_mid)
!end subroutine set_forcing

! ============================================================
!  dT/dt for snow layers
! ============================================================
subroutine dTdt(t, Tv, R_a2s, q_solar, q_rain, q_evap, &
            temp_vec, itemp_size, dt_data, R_nm, Cs_layer, Tg, dT)
    ! dT/dt for snow layers
    ! Tv = [T1, T2, T3]
    ! R_nm = [R_12, R_23, R_3g]
    ! Returns: dT = [dT1, dT2, dT3]
    use enhanced_features_module
    implicit none
    double precision, intent(in) :: t, dt_data
    double precision, dimension(3), intent(in) :: Tv
    double precision, dimension(3), intent(in) :: R_nm
    double precision, intent(in) :: R_a2s, q_solar, q_rain, q_evap
    integer, intent(in) :: itemp_size
    double precision, dimension(itemp_size), intent(in) :: temp_vec
    double precision, intent(in) :: Cs_layer, Tg
    double precision, dimension(3), intent(out) :: dT

    double precision :: T1, T2, T3, Ta, q_a, q_surf, q_12, dT1, R_3g
    double precision :: q_21, q_23, dT2, q_32, q_3g, dT3, R_12, R_23
    
    T1 = Tv(1)
    T2 = Tv(2)
    T3 = Tv(3)
    Ta = interpolate_data(temp_vec, t, dt_data) + 273.15

    R_12 = R_nm(1)
    R_23 = R_nm(2)
    R_3g = R_nm(3)

    q_a    = (Ta - T1) / R_a2s
    q_surf = q_a + q_solar + q_rain + q_evap

    ! layer 1
    q_12 = (T2 - T1) / R_12
    dT1  = (q_surf + q_12) / Cs_layer

    ! layer 2
    q_21 = (T1 - T2) / R_12
    q_23 = (T3 - T2) / R_23
    dT2  = (q_21 + q_23) / Cs_layer

    ! layer 3
    q_32 = (T2 - T3) / R_23
    q_3g = (Tg - T3) / R_3g
    dT3  = (q_32 + q_3g) / Cs_layer

    dT = [dT1, dT2, dT3]

end subroutine dTdt

SUBROUTINE step_RK4(t, dt, T_current, R_eff, R_ins, q_solar, q_rain, &
                   q_evap, temp_vec, itemp_size, dt_data, &
                   R_nm, Cs_layer, Tg, T_new)
    IMPLICIT NONE
    ! --- Input Scalars (Read-only) ---
    DOUBLE PRECISION, INTENT(IN) :: t, dt, R_eff, R_ins, q_solar, q_rain, q_evap
    DOUBLE PRECISION, INTENT(IN) :: Cs_layer, Tg, dt_data
    INTEGER, INTENT(IN)          :: itemp_size 
    
    ! --- Input Arrays (Read-only) ---
    DOUBLE PRECISION, DIMENSION(itemp_size), INTENT(IN) :: temp_vec
    DOUBLE PRECISION, DIMENSION(3), INTENT(IN)          :: T_current, R_nm
    
    ! --- Output Array (Written by this routine) ---
    DOUBLE PRECISION, DIMENSION(3), INTENT(OUT)         :: T_new
    
    ! --- Local Variables (Internal use only) ---
    DOUBLE PRECISION :: R_a2s
    DOUBLE PRECISION, DIMENSION(3) :: k1_vec, k2_vec, k3_vec, k4_vec
    
    R_a2s = R_eff + R_ins
    
    call dTdt(t, T_current, R_a2s, q_solar, q_rain, q_evap, &
             temp_vec, itemp_size, dt_data, R_nm, Cs_layer, Tg, k1_vec)
    call dTdt(t + dt/2.0D0, T_current + dt*k1_vec/2.0D0, &
             R_a2s, q_solar, q_rain, q_evap, &
             temp_vec, itemp_size, dt_data, R_nm, Cs_layer, Tg, k2_vec)
    call dTdt(t + dt/2.0D0, T_current + dt*k2_vec/2.0D0, &
             R_a2s, q_solar, q_rain, q_evap, &
             temp_vec, itemp_size, dt_data, R_nm, Cs_layer, Tg, k3_vec)
    call dTdt(t + dt, T_current + dt*k3_vec, &
             R_a2s, q_solar, q_rain, q_evap, &
             temp_vec, itemp_size, dt_data, R_nm, Cs_layer, Tg, k4_vec)
    
    T_new = T_current + (dt/6.0D0) * (k1_vec + 2.0D0*k2_vec + 2.0D0*k3_vec + k4_vec)
END SUBROUTINE step_RK4

! ------------------------------------------------------------
!  OPTION 2: DOPRI5 (Adaptive RK, explicit, no Jacobian)
! ------------------------------------------------------------
SUBROUTINE integrate_DOPRI5(t, dt, T_current, R_ins, R_eff, q_solar, q_rain, &
                            q_evap, R_nm, Cs_layer, Tg, temp_vec, itemp_size, &
                            dt_data, T_new)
    IMPLICIT NONE
    DOUBLE PRECISION :: t, dt, R_ins, R_eff, q_solar, q_rain, q_evap
    DOUBLE PRECISION :: Cs_layer, Tg, dt_data
    INTEGER :: itemp_size
    DOUBLE PRECISION, DIMENSION(itemp_size) :: temp_vec
    DOUBLE PRECISION, DIMENSION(3) :: T_current, R_nm, T_new
    
    ! DOPRI5 parameters
    INTEGER, PARAMETER :: N = 3
    ! INCREASE THESE: Give the solver plenty of breathing room
    INTEGER, PARAMETER :: LWORK = 200 
    INTEGER, PARAMETER :: LIWORK = 100
    
    DOUBLE PRECISION :: WORK(LWORK)
    INTEGER :: IWORK(LIWORK)
    
    ! RPAR needs to hold the 11 scalars + the meteorological vector
    DOUBLE PRECISION :: RPAR(11 + itemp_size) 
    INTEGER :: IPAR(10) ! Declare as an array, not a single integer
    
    INTEGER :: ITOL, IOUT, IDID
    DOUBLE PRECISION :: RTOL, ATOL, XEND
    DOUBLE PRECISION :: Y(N)
    
    EXTERNAL :: snow_rhs_dopri5, solout_dummy
    
    ! Initialize work arrays to zero to prevent garbage checks
    WORK = 0.0D0
    IWORK = 0
    
    ! 1. Pack Scalars into RPAR
    RPAR(1)  = R_eff
    RPAR(2)  = R_ins
    RPAR(3)  = q_solar
    RPAR(4)  = q_rain
    RPAR(5)  = q_evap
    RPAR(6)  = R_nm(1)
    RPAR(7)  = R_nm(2)
    RPAR(8)  = R_nm(3)
    RPAR(9)  = Cs_layer
    RPAR(10) = Tg
    RPAR(11) = dt_data
    
    ! 2. Pack the Vector into RPAR starting at index 12
    RPAR(12 : 11 + itemp_size) = temp_vec
    
    ! 3. Set IPAR array
    IPAR(1) = itemp_size ! Store the size in the first slot
    
    Y = T_current
    ITOL = 0      ! Scalar tolerances
    RTOL = 1.0D-8 ! Tighten tolerances slightly for stability
    ATOL = 1.0D-8
    XEND = t + dt
    IOUT = 0      ! No dense output
    
    ! IMPORTANT: Check if solver reached the end
    CALL DOPRI5(N, snow_rhs_dopri5, t, Y, XEND, &
                RTOL, ATOL, ITOL, &
                solout_dummy, IOUT, &
                WORK, LWORK, IWORK, LIWORK, RPAR, IPAR, IDID)
    
    if (IDID < 0) then
        print *, "DOPRI5 Warning: Solver returned IDID = ", IDID
    end if
    
    T_new = Y
END SUBROUTINE integrate_DOPRI5

SUBROUTINE snow_rhs_dopri5(N, T_TIME, Y, F, RPAR, IPAR)
    use enhanced_features_module
    IMPLICIT NONE
    INTEGER :: N, IPAR(*)            ! Explicitly an array
    DOUBLE PRECISION :: T_TIME, Y(N), F(N), RPAR(*)
    
    ! Internal variables
    DOUBLE PRECISION :: R_eff, R_ins, R_a2s, q_solar, q_rain, q_evap
    DOUBLE PRECISION :: R_12, R_23, R_3g, Cs_layer, Tg, dt_data, Ta
    INTEGER :: v_size
    
    ! --- UNPACKING ---
    ! Get the integer size from the Integer bucket
    v_size   = IPAR(1) 
    
    ! Get the scalars from the Real bucket
    R_eff    = RPAR(1)
    R_ins    = RPAR(2)
    q_solar  = RPAR(3)
    q_rain   = RPAR(4)
    q_evap   = RPAR(5)
    R_12     = RPAR(6)
    R_23     = RPAR(7)
    R_3g     = RPAR(8)
    Cs_layer = RPAR(9)
    Tg       = RPAR(10)
    dt_data  = RPAR(11)
    
    ! Ta interpolation:
    ! Use v_size to correctly slice the array starting at index 12
    Ta = interpolate_data(RPAR(12 : 11 + v_size), T_TIME, dt_data) + 273.15
    
    R_a2s = R_eff + R_ins
    
    ! Physics (same as before)
    F(1) = ((Ta - Y(1))/R_a2s + q_solar + q_rain + q_evap + (Y(2)-Y(1))/R_12) / Cs_layer
    F(2) = ((Y(1)-Y(2))/R_12 + (Y(3)-Y(2))/R_23) / Cs_layer
    F(3) = ((Y(2)-Y(3))/R_23 + (Tg-Y(3))/R_3g) / Cs_layer
END SUBROUTINE snow_rhs_dopri5

! Dummy output routine (required by DOPRI5)
SUBROUTINE solout_dummy(NR, XOLD, X, Y, N, CON, ICOMP, ND, RPAR, IPAR, IRTRN)
    IMPLICIT NONE
    INTEGER :: NR, N, ND, ICOMP(ND), IPAR, IRTRN
    DOUBLE PRECISION :: XOLD, X, Y(N), CON(*), RPAR(*)
    IRTRN = 0
    ! Avoid warnings
    IF (NR < 0 .OR. XOLD < 0.0D0 .OR. X < 0.0D0) CONTINUE
    IF (N < 0 .OR. ND < 0 .OR. IPAR < 0) CONTINUE
    IF (Y(1) < -999.0D0 .OR. CON(1) < -999.0D0) CONTINUE
    IF (ICOMP(1) < 0 .OR. RPAR(1) < -999.0D0) CONTINUE
END SUBROUTINE solout_dummy

! ------------------------------------------------------------
!  OPTION 2: LSODA (Auto stiff/non-stiff switching)
! ------------------------------------------------------------
SUBROUTINE integrate_LSODA(t, dt, T_current, R_ins, R_eff, q_solar, q_rain, &
                           q_evap, R_nm, Cs_layer, Tg, temp_vec, itemp_size, &
                           dt_data, T_new)
    USE odepack_interface
    USE odepack_common
    IMPLICIT NONE
    
    ! Input arguments
    DOUBLE PRECISION :: t, dt, R_ins, R_eff, q_solar, q_rain, q_evap
    DOUBLE PRECISION :: Cs_layer, Tg, dt_data
    INTEGER :: itemp_size
    DOUBLE PRECISION, DIMENSION(itemp_size) :: temp_vec
    DOUBLE PRECISION, DIMENSION(3) :: T_current, R_nm, T_new
    
    ! --- COMMON BLOCK FOR PARAMETERS AND WEATHER ---
    ! We must use a fixed size for the weather array here
    INTEGER, PARAMETER :: MAX_WEATHER = 100000 
    DOUBLE PRECISION :: C_R_eff, C_R_ins, C_q_solar, C_q_rain, C_q_evap
    DOUBLE PRECISION :: C_R_12, C_R_23, C_R_3g, C_Cs_layer, C_Tg
    DOUBLE PRECISION :: C_temp_vec(MAX_WEATHER), C_dt_data
    INTEGER :: C_itemp_size
    COMMON /SNOW_PARAMS/ C_R_eff, C_R_ins, C_q_solar, C_q_rain, C_q_evap, &
                        C_R_12, C_R_23, C_R_3g, C_Cs_layer, C_Tg, &
                        C_temp_vec, C_dt_data, C_itemp_size
    
    ! LSODA parameters
    INTEGER, PARAMETER :: NEQ = 3
    INTEGER :: ITOL, ITASK, ISTATE, IOPT, JT, LRW, LIW
    DOUBLE PRECISION :: RTOL, ATOL(NEQ), TOUT
    DOUBLE PRECISION :: Y(NEQ)
    DOUBLE PRECISION, ALLOCATABLE :: RWORK(:)
    INTEGER, ALLOCATABLE :: IWORK(:)
    TYPE(odepack_common_data), TARGET :: common_data
    
    EXTERNAL :: snow_rhs_lsoda, jac_dummy
    
    ! Pack data into the COMMON block
    C_R_eff    = R_eff
    C_R_ins    = R_ins
    C_q_solar  = q_solar
    C_q_rain   = q_rain
    C_q_evap   = q_evap
    C_R_12     = R_nm(1)
    C_R_23     = R_nm(2)
    C_R_3g     = R_nm(3)
    C_Cs_layer = Cs_layer
    C_Tg       = Tg
    C_dt_data  = dt_data
    C_itemp_size = itemp_size
    C_temp_vec(1:itemp_size) = temp_vec(1:itemp_size)
    
    ! LSODA setup
    LRW = MAX(20 + 16*NEQ, 22 + 9*NEQ + NEQ*NEQ)
    LIW = 20 + NEQ
    ALLOCATE(RWORK(LRW), IWORK(LIW))
    
    Y = T_current
    IWORK = 0
    RWORK = 0.0D0
    common_data%ierr = 0
    ITOL = 2
    RTOL = 1.0D-8
    ATOL = 1.0D-8
    ITASK = 1
    ISTATE = 1
    IOPT = 0
    JT = 2
    TOUT = t + dt
    
    CALL DLSODA(snow_rhs_lsoda, NEQ, Y, t, TOUT, ITOL, RTOL, ATOL, &
                ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, &
                jac_dummy, JT, common_data)
    
    IF (ISTATE < 0) THEN
        WRITE(*,'(A,I3)') 'DLSODA error: ISTATE = ', ISTATE
        STOP
    END IF
    
    T_new = Y
    DEALLOCATE(RWORK, IWORK)
END SUBROUTINE integrate_LSODA

! Right-hand side for LSODA
SUBROUTINE snow_rhs_lsoda(NEQ, T_TIME, Y, YDOT, common_data)
    USE odepack_common
    use enhanced_features_module
    IMPLICIT NONE
    INTEGER :: NEQ
    DOUBLE PRECISION :: T_TIME, Y(NEQ), YDOT(NEQ)
    TYPE(odepack_common_data) :: common_data
    
    ! --- MATCHING COMMON BLOCK ---
    INTEGER, PARAMETER :: MAX_WEATHER = 100000
    DOUBLE PRECISION :: R_eff, R_ins, q_solar, q_rain, q_evap
    DOUBLE PRECISION :: R_12, R_23, R_3g, Cs_layer, Tg
    DOUBLE PRECISION :: temp_vec(MAX_WEATHER), dt_data
    INTEGER :: itemp_size
    COMMON /SNOW_PARAMS/ R_eff, R_ins, q_solar, q_rain, q_evap, &
                        R_12, R_23, R_3g, Cs_layer, Tg, &
                        temp_vec, dt_data, itemp_size
    
    DOUBLE PRECISION :: R_a2s, Ta, q_a, q_surf, q_12, q_21, q_23, q_32, q_3g
    !DOUBLE PRECISION, EXTERNAL :: interpolate_data
    
    R_a2s = R_eff + R_ins
    
    ! Ta is now interpolated directly from the weather vector in COMMON
    Ta = interpolate_data(temp_vec(1:itemp_size), T_TIME, dt_data) + 273.15
    
    ! Physics calculation
    q_a = (Ta - Y(1)) / R_a2s
    q_surf = q_a + q_solar + q_rain + q_evap
    
    q_12 = (Y(2) - Y(1)) / R_12
    YDOT(1) = (q_surf + q_12) / Cs_layer
    
    q_21 = (Y(1) - Y(2)) / R_12
    q_23 = (Y(3) - Y(2)) / R_23
    YDOT(2) = (q_21 + q_23) / Cs_layer
    
    q_32 = (Y(2) - Y(3)) / R_23
    q_3g = (Tg - Y(3)) / R_3g
    YDOT(3) = (q_32 + q_3g) / Cs_layer
    
    common_data%ierr = 0
END SUBROUTINE snow_rhs_lsoda

! Dummy Jacobian for LSODA (JT=2 means internally generated)
SUBROUTINE jac_dummy(NEQ, T, Y, ML, MU, PD, NROWPD, common_data)
    USE odepack_common
    IMPLICIT NONE
    INTEGER :: NEQ, ML, MU, NROWPD
    DOUBLE PRECISION :: T, Y(NEQ), PD(NROWPD, NEQ)
    TYPE(odepack_common_data) :: common_data
    
    PD = 0.0D0
    common_data%ierr = 0
    
    ! Avoid warnings
    IF (NEQ < 0 .OR. ML < 0 .OR. MU < 0 .OR. NROWPD < 0) CONTINUE
    IF (T < -999.0D0 .OR. Y(1) < -999.0D0) CONTINUE
END SUBROUTINE jac_dummy

! ============================================================
! Enhanced energy diagnostics with refreezing and runoff
! ============================================================
!SUBROUTINE compute_energy_diagnostics_enhanced(t_vec, Nt, qa_hist, qsolar_hist, &
!                                              qrain_hist, qevap_hist, qground_hist, &
!                                              T_hist, Cs_layer, E_melt, rho_i, Lf, &
!                                              refrozen_hist, runoff_hist, n_hours, &
!                                              USE_ADVANCED, USE_MULTILAYER, &
!                                              USE_REFREEZE, USE_PERCOL, &
!                                              N_ins, SOLVER_CHOICE)
!    IMPLICIT NONE
!    INTEGER, INTENT(IN) :: Nt, SOLVER_CHOICE, n_hours, N_ins
!    DOUBLE PRECISION, DIMENSION(Nt), INTENT(IN) :: t_vec, qa_hist, qsolar_hist
!    DOUBLE PRECISION, DIMENSION(Nt), INTENT(IN) :: qrain_hist, qevap_hist, qground_hist
!    DOUBLE PRECISION, DIMENSION(Nt), INTENT(IN) :: refrozen_hist, runoff_hist
!    DOUBLE PRECISION, DIMENSION(Nt, 3), INTENT(IN) :: T_hist
!    DOUBLE PRECISION, INTENT(IN) :: Cs_layer, E_melt, rho_i, Lf
!    LOGICAL, INTENT(IN) :: USE_ADVANCED, USE_MULTILAYER, USE_REFREEZE, USE_PERCOL
!    
!    DOUBLE PRECISION :: E_a, E_solar, E_rain, E_evap, E_g
!    DOUBLE PRECISION :: E_total_in, E_snow_change, E_refrozen, E_balance
!    DOUBLE PRECISION :: M_melt, M_runoff, melt_rate_avg
!    
!    DOUBLE PRECISION, EXTERNAL ::

! ============================================================
!  EXAMPLE: Integration into main time loop
! ============================================================
! This shows how to use the enhanced features in your main program
!
! Add to variable declarations in main program:
!   use enhanced_features_module
!   type(CSVData) :: met_data
!   type(LayerProperties), dimension(3) :: snow_layers
!   logical :: USE_CSV_DATA, USE_REFREEZING, USE_PERCOLATION, USE_MULTILAYER
!   double precision :: theta_e, total_refrozen, runoff
!   double precision, dimension(3) :: LWC, ice_fractions, layer_heights
!   double precision, allocatable :: refrozen_hist(:), runoff_hist(:), LWC_hist(:,:)
!   integer :: N_ins_layers
!   double precision :: D_insulation, k_i_base
!
! Initialization (before time loop):
!   USE_CSV_DATA = .TRUE.
!   USE_REFREEZING = .TRUE.
!   USE_PERCOLATION = .TRUE.
!   USE_MULTILAYER = .TRUE.
!   
!   ! Initialize layer properties
!   theta_e = 0.03D0  ! Field capacity
!   LWC = [0.0D0, 0.0D0, 0.0D0]
!   ice_fractions = [0.4D0, 0.4D0, 0.4D0]
!   layer_heights = [dz_s, dz_s, dz_s]
!   
!   ! Multi-layer insulation parameters
!   N_ins_layers = 20
!   k_i_base = 0.32D0
!   rho_dry = 100.0D0
!   moist_cont = 50.0D0
!   rho_wet = rho_dry + moist_cont/100.0D0*1000.0D0
!   c_dry = 0.99D3
!   c_wet = (1.0D0 - moist_cont/100.0D0)*c_dry + moist_cont/100.0D0*c_w
!   D_insulation = k_i_base / (c_wet * rho_wet)
!   
!   ! Allocate history arrays
!   allocate(refrozen_hist(Nt), runoff_hist(Nt), LWC_hist(Nt, 3))
!   refrozen_hist = 0.0D0
!   runoff_hist = 0.0D0
!   LWC_hist = 0.0D0
!   LWC_hist(1, :) = LWC
!   
!   ! Read CSV data if enabled
!   if (USE_CSV_DATA) then
!       call read_csv_data('DATA.csv', met_data, status)
!       if (status /= 0) then
!           print *, "Error reading CSV file. Using synthetic forcing."
!           USE_CSV_DATA = .FALSE.
!       else
!           print *, "Loaded ", met_data%n_points, " hourly data points"
!       end if
!   end if
!
! Inside time loop (after computing T_new, before updating Tempe):
!   
!   ! Get forcing data
!   if (USE_CSV_DATA) then
!       Ta_C = interpolate_data(met_data%temp, met_data%n_points, t_mid, 3600.0D0)
!       forc%Ta = Ta_C + 273.15D0
!       forc%Isolar = interpolate_data(met_data%solar, met_data%n_points, t_mid, 3600.0D0)
!       forc%Prain = interpolate_data(met_data%precip, met_data%n_points, t_mid, 3600.0D0) / 3600.0D0
!       wind_speed = interpolate_data(met_data%wind, met_data%n_points, t_mid, 3600.0D0)
!       RH_pct = interpolate_data(met_data%rh, met_data%n_points, t_mid, 3600.0D0)
!       forc%RH = RH_pct / 100.0D0
!   else
!       ! Use synthetic forcing functions
!       forc%Ta = Ta_fun(t_mid)
!       forc%Isolar = Isolar_fun(t_mid)
!       forc%Prain = Prain_fun(t_mid)
!       forc%RH = 0.80D0
!   end if
!   
!   ! Call enhanced insulation step
!   if (USE_ADVANCED_INSULATION) then
!       call insulation_step_enhanced(InsState, forc, InsPar, dt, USE_MULTILAYER, &
!                                    N_ins_layers, D_insulation, T_new(1), &
!                                    R_ins, q_solar_ins, q_rain_snow, q_evap, InsState)
!   end if
!   
!   ! Refreezing in each layer
!   total_refrozen = 0.0D0
!   if (USE_REFREEZING) then
!       do i = 1, 3
!           call refreezing_layer(T_new(i), LWC(i), ice_fractions(i), layer_heights(i), &
!                                Tfreeze, rho_i, rho_w, c_s, c_w, Lf, &
!                                T_new(i), LWC(i), ice_fractions(i), refrozen_mass)
!           total_refrozen = total_refrozen + refrozen_mass
!       end do
!   end if
!   refrozen_hist(k) = total_refrozen
!   
!   ! Surface melting (add melt water to LWC)
!   if (T_new(1) > Tfreeze) then
!       T_new(1) = Tfreeze
!       if (q_surf_mid > 0.0D0) then
!           dE_melt = q_surf_mid * dt
!           dM_melt = dE_melt / (rho_i * Lf)
!           surface_melt_water = dM_melt  ! [m w.e.]
!           LWC(1) = LWC(1) + surface_melt_water / layer_heights(1)
!       else
!           dE_melt = 0.0D0
!           dM_melt = 0.0D0
!       end if
!       E_melt = E_melt + dE_melt
!       melt_rate_hist(k) = dM_melt / dt
!   else
!       melt_rate_hist(k) = 0.0D0
!   end if
!   
!   ! Percolation
!   if (USE_PERCOLATION) then
!       call percolate_water(LWC, layer_heights, theta_e, 3, rho_w, runoff)
!       runoff_hist(k) = runoff
!   end if
!   
!   ! Store LWC history
!   LWC_hist(k+1, :) = LWC
!
! After time loop - Enhanced energy diagnostics:
!   E_refrozen = sum(refrozen_hist) * Lf
!   E_balance = E_total_in - (E_snow_change + E_melt + E_refrozen)
!   M_runoff = sum(runoff_hist)
!   
!   print *, ""
!   print '(A, ES14.6, A)', " E_refreezing    = ", E_refrozen, " J/m²"
!   print '(A, F10.3, A)', " Total runoff    = ", M_runoff, " kg/m²"
!   
!   ! Cleanup
!   if (USE_CSV_DATA) then
!       deallocate(met_data%time, met_data%temp, met_data%wind)
!       deallocate(met_data%precip, met_data%solar, met_data%rh)
!   end if
!   deallocate(refrozen_hist, runoff_hist, LWC_hist)