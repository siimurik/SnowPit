program main
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
        double precision :: W, age_days            ! State variables
        double precision :: k_eff, alpha_eff, f_sat  ! Diagnostic outputs
    end type InsulationState
    
    type(InsulationParameters) :: InsPar
    type(InsulationState)      :: InsState

    type Forcing
        double precision :: Isolar, Prain, T_rain, RH, Ta
    end type Forcing
    
    ! EXTERNAL function declarations
    !double precision, external :: arange

    ! Variable declarations
    logical :: USE_ADVANCED_INSULATION
    integer :: Nt
    double precision :: sigma, Lf, rho_i, rho_s, c_s, rho_w, c_w, Tfreeze
    double precision :: Hs, Ns, dz_s, k_snow, Hi, k_i_const, alpha_const
    double precision :: eta_rain_const, Hg_ins, kg_ins, h_conv, epsilon
    double precision :: T_mean, h_rad, h_eff, R_eff, R_layer, R_12, R_23
    double precision :: R_g_ins, R_3g, Tg, Cs_layer, T1_init, T2_init, T3_init
    double precision, dimension(3) :: T, R_nm
    double precision :: t0, tf, delta_t
    double precision, allocatable, dimension(:) :: t_vec

    ! Allocate 2D array for temperature history: Nt x 3
    double precision, dimension(:,:), allocatable :: T_hist

    ! Allocate 1D arrays for various histories
    double precision, dimension(:), allocatable :: qnet_surf_hist
    double precision, dimension(:), allocatable :: qa_hist
    double precision, dimension(:), allocatable :: qsolar_hist
    double precision, dimension(:), allocatable :: qrain_hist
    double precision, dimension(:), allocatable :: qevap_hist
    double precision, dimension(:), allocatable :: qground_hist

    double precision, dimension(:), allocatable :: Ta_hist
    double precision, dimension(:), allocatable :: Isolar_hist
    double precision, dimension(:), allocatable :: Prain_hist

    double precision, dimension(:), allocatable :: melt_rate_hist
    double precision :: E_melt

    ! Insulation diagnostics
    double precision, dimension(:), allocatable :: W_hist
    double precision, dimension(:), allocatable :: k_eff_hist
    double precision, dimension(:), allocatable :: alpha_hist
    double precision, dimension(:), allocatable :: Rins_hist
    double precision, dimension(:), allocatable :: fsat_hist
    

    ! ============================================================
    !  3-layer Snow Storage RC Model with switchable insulation
    ! ============================================================

    USE_ADVANCED_INSULATION = .TRUE.  ! True: moisture model, False: constant R

    ! ---------- Physical constants ----------
    sigma   = 5.670374419D-8      ! Stefan-Boltzmann [W/m^2 K^4]
    Lf      = 3.34d5              ! Latent heat of fusion [J/kg]
    rho_i   = 917.0D0             ! Ice density [kg/m^3]
    rho_s   = 400.0D0             ! Snow bulk density [kg/m^3]
    c_s     = 2100.0D0            ! Snow specific heat [J/kg K]
    rho_w   = 1000.0D0            ! Water density [kg/m^3]
    c_w     = 4180.0D0            ! Water specific heat [J/kg K]
    Tfreeze = 273.15D0            ! 0°C [K]

    ! ---------- Snow & insulation ----------
    Hs   = 2.0D0                  ! total snow thickness [m]
    Ns   = 3.0D0                  ! number of snow layers
    dz_s = Hs / Ns                ! thickness per snow layer [m]

    k_snow = 0.25D0               ! snow conductivity [W/mK]
    Hi     = 0.6D0                ! insulation thickness [m]

    ! Simple (constant) insulation parameters
    k_i_const     = 0.06D0        ! [W/mK]
    alpha_const   = 0.10D0        ! [-]
    eta_rain_const = 1.0D0        ! fraction of rain heat reaching snow

    ! Ground insulation
    Hg_ins = 0.3D0                ! [m]
    kg_ins = 0.04D0               ! [W/mK]

    ! ---------- Surface HTC (air-side, conv + LW) ----------
    h_conv  = 8.0D0               ! convective coefficient [W/m^2K]
    epsilon = 0.95D0
    T_mean  = 273.15D0 + 3.0D0    ! nominal mean temperature [K]
    h_rad   = 4.0D0 * epsilon * sigma * T_mean**3
    h_eff   = h_conv + h_rad

    R_eff = 1.0D0 / h_eff         ! air-side resistance [m^2K/W]

    ! ---------- Thermal resistances in snow & ground ----------
    R_layer = dz_s / k_snow
    R_12    = R_layer
    R_23    = R_layer

    R_g_ins = Hg_ins / kg_ins
    R_3g    = R_layer + R_g_ins
    R_nm = [R_12, R_23, R_3g]

    ! ---------- Ground & snow capacity ----------
    Tg        = 273.15D0 + 2.0D0
    Cs_layer  = rho_s * c_s * dz_s   ! [J/(m^2 K)] per snow layer

    ! ---------- Initial snow temperatures ----------
    T1_init = 273.15D0 - 2.0D0
    T2_init = 273.15D0 - 4.0D0
    T3_init = 273.15D0 - 6.0D0
    T = [T1_init, T2_init, T3_init]

    ! ============================================================
    !  Advanced insulation parameters
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

            Lv = 2.5d6, &        
            rho_w = rho_w, &
            c_w = c_w, &
            Tfreeze = Tfreeze, &

            ! better evaporation model
            rho_air = 1.2D0, &      ! [kg/m^3]
            C_E = 1.3D-3, &         ! bulk transfer coeff
            U10 = 2.0D0, &          ! [m/s]
            P0 = 101325.0D0 &       ! [Pa]
        )
        
        InsState = InsulationState(&
            W = 5.0D0, &            ! [kg/m^2]
            age_days = 0.0D0, &     ! days
            k_eff = 0.D0, &
            alpha_eff = 0.D0, &
            f_sat = 0.D0 &
            )
    else
        ! Set InsPar as not used
        InsState = InsulationState( &
            W = 0.0D0, &            ! [kg/m^2]
            age_days = 0.0D0, &     ! days
            k_eff = 0.D0, &
            alpha_eff = 0.D0, &
            f_sat = 0.D0 &
            )
    end if

    ! Optional: Print some values to verify
    print *, "Physical constants:"
    print '(A, F12.8)', " sigma = ", sigma
    print '(A, F12.4)', " rho_w = ", rho_w
    print '(A, F12.2)', " c_w = ", c_w
    print '(A, F12.6)', " Tfreeze = ", Tfreeze
        
    if (USE_ADVANCED_INSULATION) then
        print *, "Advanced insulation:"
        print '(A, F8.4)', " InsPar%k_dry = ", InsPar%k_dry
        print '(A, F8.4)', " InsPar%k_sat = ", InsPar%k_sat
        print '(A, F8.2)', " InsState%W = ", InsState%W
    end if

    ! ============================================================
    !  Time integration settings
    ! ============================================================
    t0 = 0.0D0
    tf = 120.0D0 * 24.0D0 * 3600.0D0   ! 120 days
    delta_t = 600.0D0                   ! 10 minutes
    
    ! Use arange function directly (no pre-allocation needed)
    t_vec = arange(t0, tf + delta_t, delta_t)
    Nt = size(t_vec)  ! Get actual size from function
    
    allocate(T_hist(Nt, 3))
    T_hist = 0.0D0
    T_hist(1, :) = T  ! Fortran is 1-indexed, not 0-indexed like Python

    ! Histories
    allocate(qnet_surf_hist(Nt), qa_hist(Nt), qsolar_hist(Nt), &
            qrain_hist(Nt), qevap_hist(Nt), qground_hist(Nt))

    ! Initialize all arrays to zero
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

    allocate(melt_rate_hist(Nt))
    melt_rate_hist = 0.0D0
    E_melt = 0.0D0

    allocate(W_hist(Nt), k_eff_hist(Nt), alpha_hist(Nt), &
            Rins_hist(Nt), fsat_hist(Nt))
    W_hist = 0.0D0
    k_eff_hist = 0.0D0
    alpha_hist = 0.0D0
    Rins_hist = 0.0D0
    fsat_hist = 0.0D0

    ! Free memory
    deallocate(t_vec)
    deallocate(T_hist)
    deallocate(qnet_surf_hist, qa_hist, qsolar_hist)
    deallocate(qrain_hist, qevap_hist, qground_hist)
    deallocate(Ta_hist, Isolar_hist, Prain_hist)
    deallocate(melt_rate_hist)
    deallocate(W_hist, k_eff_hist, alpha_hist)
    deallocate(Rins_hist, fsat_hist)

contains

    ! ============================================================
    !  Insulation step (advanced model)
    ! ============================================================
    subroutine insulation_step(state_in, forc, p, dt, R_ins, q_solar, &
                            q_rain_snow, q_evap, state_out)
        ! Advanced insulation model with moisture, age, rain, and evaporation
        implicit none
        
        ! Inputs
        type(InsulationState), intent(in) :: state_in
        type(Forcing), intent(in) :: forc
        type(InsulationParameters), intent(in) :: p
        double precision, intent(in) :: dt
        
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
        
        ! Extract state variables
        W = state_in%W
        age_days = state_in%age_days
        
        ! 1) Moisture & age factors
        f = max(0.0D0, min(1.0D0, W / p%W_sat))  ! saturation fraction [0,1]
        age_yr = age_days / 365.0D0
        
        ! conductivity: k(W, age)
        k_moist = p%k_dry + (p%k_sat - p%k_dry) * (f**p%n_k)
        k_age_factor = 1.0D0 + p%delta_k_age * (1.0D0 - exp(-age_yr / p%tau_k_years))
        k_eff = k_moist * k_age_factor
        R_ins = p%Hi / k_eff
        
        ! absorptivity: alpha(W, age)
        alpha_moist = p%alpha_dry + (p%alpha_wet - p%alpha_dry) * (f**p%n_alpha)
        alpha_age = alpha_moist + p%delta_alpha_age * (1.0D0 - exp(-age_yr / p%tau_alpha_years))
        alpha_eff = max(0.0D0, min(1.0D0, alpha_age))  ! clip to [0,1]
        
        ! 2) Solar absorption
        q_solar = alpha_eff * forc%Isolar
        
        ! 3) Rain infiltration & heat to snow
        eta_rain = max(0.0D0, 1.0D0 - f)
        P_in_mass = eta_rain * p%rho_w * forc%Prain   ! [kg/m^2/s]
        
        zeta_rain = p%zeta0 * exp(-p%gamma_H * p%Hi) * exp(-p%gamma_W * f)
        q_rain_snow = zeta_rain * p%rho_w * p%c_w * forc%Prain * (forc%T_rain - p%Tfreeze)
        
        ! 4) Evaporation (bulk aerodynamic model)
        Tc_s = p%Tfreeze - 273.15D0      ! ~0°C
        Tc_a = forc%Ta - 273.15D0
        
        ! saturation vapour pressures [Pa] - Tetens formula
        e_sat_s = 611.0D0 * exp(17.27D0 * Tc_s / (Tc_s + 237.3D0))
        e_sat_a = 611.0D0 * exp(17.27D0 * Tc_a / (Tc_a + 237.3D0))
        
        e_surf = e_sat_s           ! saturated surface
        e_air = forc%RH * e_sat_a
        
        VPD = max(0.0D0, e_surf - e_air)    ! Vapour Pressure Deficit [Pa]
        
        ! bulk aerodynamic evaporation [kg/m^2/s]
        E0 = p%rho_air * p%C_E * p%U10 * VPD / p%P0
        
        f_breath = exp(-p%beta_w * f)
        E = E0 * f_breath
        
        q_evap = -p%Lv * E      ! [W/m^2] (negative = cooling)
        
        ! 5) Drainage & moisture update
        D = p%K_D * max(0.0D0, W - p%W_field)
        dWdt = P_in_mass - E - D
        W_new = W + dt * dWdt
        W_new = max(0.0D0, min(p%W_sat, W_new))  ! clip to [0, W_sat]
        
        age_days_new = age_days + dt / 86400.0D0
        
        ! Update output state
        state_out = state_in  ! copy all fields
        state_out%W = W_new
        state_out%age_days = age_days_new
        
        ! Store diagnostics
        state_out%k_eff = k_eff
        state_out%alpha_eff = alpha_eff
        state_out%f_sat = f
        
    end subroutine insulation_step

    function arange(start, stop, step) result(arr)
        ! Fortran version of numpy.arange
        implicit none
        double precision, intent(in) :: start, stop, step
        double precision, allocatable :: arr(:)
        
        integer :: n, i
        
        ! Calculate number of elements
        if (abs(step) < tiny(1.0D0)) then
            n = 1
        else
            n = floor((stop - start) / step) + 1
        end if
        
        ! Python excludes stop if reached exactly
        if (step > 0.0D0 .and. start + (n-1)*step >= stop) then
            n = n - 1
        else if (step < 0.0D0 .and. start + (n-1)*step <= stop) then
            n = n - 1
        end if
        
        ! Ensure at least 1 element
        n = max(1, n)
        
        ! Now allocate the result array
        allocate(arr(n))
        
        ! Fill array
        do i = 1, n
            arr(i) = start + (i-1)*step
        end do
        
    end function arange

end program main

! ============================================================
!  Forcing functions
! ============================================================
double precision function Ta_fun(t) result(T_air)
    ! Air temperature [k] with daily cycle
    ! elemental: can be called with scalar or array arguments
    implicit none
    double precision, intent(in) :: t
    double precision, parameter :: day = 24.D0 * 3600.D0
    double precision, parameter :: Tmean = 273.15D0 + 3.D0
    double precision, parameter :: Tamp  = 6.D0
    double precision, parameter :: pi = 4.D0*atan(1.D0)

    T_air = Tmean + Tamp * sin(2.D0 * pi * (t / day))
end function Ta_fun

elemental double precision function Isolar_fun(t) result(I_solar)
    ! Diurnal solar [W/m^2], half-sine from 6h to 18h.
    implicit none
    double precision, intent(in) :: t
    double precision, parameter :: DAY_SECONDS = 86400.0D0      ! 24*3600
    double precision, parameter :: DAWN = 21600.0D0             ! 6*3600
    double precision, parameter :: DUSK = 64800.0D0             ! 18*3600  
    double precision, parameter :: PI = acos(-1.0D0)
    double precision :: tau
    
    tau = mod(t, DAY_SECONDS)
    
    if (tau < DAWN .or. tau > DUSK) then
        I_solar = 0.0D0
    else
        I_solar = 600.0D0 * sin(PI * (tau - DAWN) / (DUSK - DAWN))
    end if
end function Isolar_fun

elemental double precision function Prain_fun(t) result(Prain)
    ! Rain rate [m/s]: 2h event every 5 days.
    implicit none
    double precision, intent(in) :: t
    double precision, parameter :: day = 24.D0 * 3600.D0
    double precision, parameter :: period = 5.D0 * day
    double precision :: tau
    double precision, parameter :: EVENT_DURATION = 2.0 * 3600.0  ! 2 hours in seconds

    tau = mod(t, period)

    if (tau < EVENT_DURATION) then
        Prain = 5.D-7
    else 
        Prain = 0.D0
    end if
end function Prain_fun

! ============================================================
!  Ground flux
! ============================================================
elemental double precision function ground_flux(T3, Tg, R_3g) result(flux)
    ! Ground -> bottom snow layer flux [W/m^2].
    implicit none
    double precision, intent(in) :: T3, Tg, R_3g

    flux = (Tg - T3) / R_3g

end function ground_flux

! ============================================================
!  dT/dt for snow layers
! ============================================================
function dTdt(t, Tv, R_a2s, q_solar, q_rain, q_evap, &
                        R_nm, Cs_layer, Tg) result(dT)
    ! dT/dt for snow layers
    ! Tv = [T1, T2, T3]
    ! R_nm = [R_12, R_23, R_3g]
    ! Returns: dT = [dT1, dT2, dT3]
    implicit none
    double precision, intent(in) :: t
    double precision, dimension(3), intent(in) :: Tv
    double precision, dimension(3), intent(in) :: R_nm
    double precision, intent(in) :: R_a2s, q_solar, q_rain, q_evap
    double precision, intent(in) :: Cs_layer, Tg
    double precision, dimension(3) :: dT

    double precision, external :: Ta_fun
    double precision :: T1, T2, T3, Ta, q_a, q_surf, q_12, dT1, R_3g
    double precision :: q_21, q_23, dT2, q_32, q_3g, dT3, R_12, R_23
    
    T1 = Tv(1)
    T2 = Tv(2)
    T3 = Tv(3)
    Ta = Ta_fun(t)

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

end function dTdt

