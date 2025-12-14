! Compile with DOPRI5 subroutines:
! gfortran -O2 snowRC_adv.f90 dopri5.f90 YBER_ODEPACK.f90 ODEPACK_MODULES.f90 -o snow -lopenblas
! gfortran snowRC_adv.f90 dopri5.f90 YBER_ODEPACK.f90 ODEPACK_MODULES.f90 lapack.f90 lapackc.f90 dc_lapack.f90 -o snow 

module types_module
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
    
end module types_module

program main
    use types_module
    implicit none
    
    type(InsulationParameters) :: InsPar
    type(InsulationState)      :: InsState
    
    ! EXTERNAL function declarations
    !double precision, external :: arange
    double precision, external :: Isolar_fun, Prain_fun, Ta_fun, ground_flux

    ! Variable declarations
    logical :: USE_ADVANCED_INSULATION
    integer :: Nt, k
    double precision :: sigma, Lf, rho_i, rho_s, c_s, rho_w, c_w, Tfreeze
    double precision :: Hs, Ns, dz_s, k_snow, Hi, k_i_const, alpha_const
    double precision :: eta_rain_const, Hg_ins, kg_ins, h_conv, epsilon
    double precision :: T_mean, h_rad, h_eff, R_eff, R_layer, R_12, R_23
    double precision :: R_g_ins, R_3g, Tg, Cs_layer, T1_init, T2_init, T3_init
    double precision, dimension(3) :: Temp, R_nm
    double precision, allocatable, dimension(:) :: t_vec
    double precision :: t0, tf, dt, t, t_mid
    type(Forcing) :: forc
    double precision :: R_ins, q_solar_ins, q_rain_snow, q_evap_val
    !double precision :: R_a2s
    !double precision, dimension(3) :: k1_vec, k2_vec, k3_vec, k4_vec, 
    double precision, dimension(3) :: T_new
    !double precision :: Ta_mid, T1_mid, 
    double precision :: q_a_mid, q_surf_mid, q_ground_mid
    !double precision :: dE_melt, dM_melt

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

    INTEGER :: SOLVER_CHOICE 

    double precision :: R_eff_com, R_ins_com, q_solar_com, q_rain_com, Tg_com
    double precision :: q_evap_com, R_12_com, R_23_com, R_3g_com, Cs_layer_com

    COMMON /SNOW_PARAMS/ R_eff_com, R_ins_com, q_solar_com, q_rain_com, q_evap_com, &
                        R_12_com, R_23_com, R_3g_com, Cs_layer_com, Tg_com
    

    ! ============================================================
    !  3-layer Snow Storage RC Model with switchable insulation
    ! ============================================================

    USE_ADVANCED_INSULATION = .TRUE.  ! True: moisture model, False: constant R

    ! ============================================================
    
    print *, "Chosen integration method:"
    SOLVER_CHOICE = 3  ! 1=RK4, 2=DOPRI5, 3=LSODA
    SELECT CASE (SOLVER_CHOICE)
    CASE (1)
        print *, "RK4 - FOURTH ORDER RUNGE-KUTTA"
    CASE (2)
        print *, "DOPRI5 - DORMAND & PRINCE EXPLICIT RUNGE-KUTTA"
        print *, "         METHOD OF ORDER (4)5"
    CASE (3)
        print *, "LSODA -"
        print *, "LIVERMORE SOLVER FOR ORDINARY DIFFERENTIAL EQUATIONS, WITH"
        print *, "AUTOMATIC METHOD SWITCHING FOR STIFF AND NONSTIFF PROBLEMS."
    END SELECT

    ! ============================================================

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
    Temp = [T1_init, T2_init, T3_init]

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
    !print *, "Physical constants:"
    !print '(A, F12.8)', " sigma = ", sigma
    !print '(A, F12.4)', " rho_w = ", rho_w
    !print '(A, F12.2)', " c_w = ", c_w
    !print '(A, F12.6)', " Tfreeze = ", Tfreeze
    !    
    !if (USE_ADVANCED_INSULATION) then
    !    print *, "Advanced insulation:"
    !    print '(A, F8.4)', " InsPar%k_dry = ", InsPar%k_dry
    !    print '(A, F8.4)', " InsPar%k_sat = ", InsPar%k_sat
    !    print '(A, F8.2)', " InsState%W = ", InsState%W
    !end if

    ! ============================================================
    !  Time integration settings
    ! ============================================================
    t0 = 0.0D0
    tf = 120.0D0 * 24.0D0 * 3600.0D0    ! 120 days
    dt = 600.0D0                        ! 10 minutes
    
    ! Use arange function directly (no pre-allocation needed)
    Nt = floor((tf - t0) / dt) + 1
    allocate(t_vec(Nt))
    t_vec = arange(t0, tf + dt, dt)
    Nt = size(t_vec)  ! Get actual size from function
    
    allocate(T_hist(Nt, 3))
    T_hist       = 0.0D0
    T_hist(1, :) = Temp  ! Fortran is 1-indexed, not 0-indexed like Python

    ! Histories
    allocate(qnet_surf_hist(Nt), qa_hist(Nt), qsolar_hist(Nt), &
            qrain_hist(Nt), qevap_hist(Nt), qground_hist(Nt))

    ! Initialize all arrays to zero
    qnet_surf_hist = 0.0D0
    qa_hist        = 0.0D0
    qsolar_hist    = 0.0D0
    qrain_hist     = 0.0D0
    qevap_hist     = 0.0D0
    qground_hist   = 0.0D0

    allocate(Ta_hist(Nt), Isolar_hist(Nt), Prain_hist(Nt))
    Ta_hist     = 0.0D0
    Isolar_hist = 0.0D0
    Prain_hist  = 0.0D0

    allocate(melt_rate_hist(Nt))
    melt_rate_hist = 0.0D0
    E_melt         = 0.0D0

    allocate(W_hist(Nt), k_eff_hist(Nt), alpha_hist(Nt), &
            Rins_hist(Nt), fsat_hist(Nt))
    W_hist     = 0.0D0
    k_eff_hist = 0.0D0
    alpha_hist = 0.0D0
    Rins_hist  = 0.0D0
    fsat_hist  = 0.0D0

    ! ============================================================
    !  Time loop (RK4)
    ! ============================================================
    do k = 1, Nt-1
        t = t_vec(k)
        t_mid = t + dt / 2.0D0
        
        ! --- forcing at mid-step ---
        call set_forcing(t_mid, forc)

        ! --- insulation: choose simple or advanced model ---
        ! Compute insulation (advanced or simple)
        call get_insulation_properties(USE_ADVANCED_INSULATION, InsState, forc, InsPar, &
                                    dt, Hi, k_i_const, alpha_const, eta_rain_const, &
                                    rho_w, c_w, Tfreeze, &
                                    R_ins, q_solar_ins, q_rain_snow, q_evap_val)
        
        !! total air → surface snow resistance
        ! Choose integration method
        SELECT CASE (SOLVER_CHOICE)
        CASE (1)
            ! Original RK4
            call step_RK4(t, dt, Temp, R_eff, R_ins, q_solar_ins, q_rain_snow, &
                        q_evap_val, R_nm, Cs_layer, Tg, T_new)
        CASE (2)
            ! Adaptive DOPRI5
            call integrate_DOPRI5(t, dt, Temp, R_ins, R_eff, q_solar_ins, &
                                q_rain_snow, q_evap_val, R_nm, Cs_layer, Tg, T_new)
        CASE (3)
            ! Adaptive DOPRI5
            ! LSODA with auto stiff/non-stiff switching
            ! First set common block for LSODA

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
            
            call integrate_LSODA(t, dt, Temp, R_ins, R_eff, q_solar_ins, &
                                q_rain_snow, q_evap_val, R_nm, Cs_layer, Tg, T_new)
        
        END SELECT

        ! Compute fluxes and apply melting
        call process_surface_energy(t_mid, Temp, T_new, R_eff, R_ins, &
                                q_solar_ins, q_rain_snow, q_evap_val, &
                                Tg, R_3g, Tfreeze, dt, rho_i, Lf, &
                                q_a_mid, q_surf_mid, q_ground_mid, &
                                E_melt, melt_rate_hist(k))
        
        ! Update and store
        Temp = T_new

        ! --- store histories ---
        !Temp = T_new
        T_hist(k+1, :) = Temp
        
        qnet_surf_hist(k) = q_surf_mid
        qa_hist(k)        = q_a_mid
        qsolar_hist(k)    = q_solar_ins
        qrain_hist(k)     = q_rain_snow
        qevap_hist(k)     = q_evap_val
        qground_hist(k)   = q_ground_mid
        
        Ta_hist(k)     = Ta_fun(t)
        Isolar_hist(k) = Isolar_fun(t)
        Prain_hist(k)  = Prain_fun(t)

    end do
    
    ! last forcing values
    Ta_hist(Nt)     = Ta_fun(t_vec(Nt))
    Isolar_hist(Nt) = Isolar_fun(t_vec(Nt))
    Prain_hist(Nt)  = Prain_fun(t_vec(Nt))
    
    if (USE_ADVANCED_INSULATION) then
        W_hist(Nt)     = InsState%W
        k_eff_hist(Nt) = InsState%k_eff
        alpha_hist(Nt) = InsState%alpha_eff
        Rins_hist(Nt)  = R_ins
        fsat_hist(Nt)  = InsState%f_sat
    end if

    ! BROKEN!!!
    !call save_final_timestep(Nt, t_vec, USE_ADVANCED_INSULATION, InsState, R_ins, &
    !                    Ta_hist, Isolar_hist, Prain_hist, &
    !                    W_hist, k_eff_hist, alpha_hist, Rins_hist, fsat_hist)
    
    !print *, " W_hist     = ", W_hist(Nt)
    !print *, " k_eff_hist = ", k_eff_hist(Nt)
    !print *, " alpha_hist = ", alpha_hist(Nt)
    !print *, " Rins_hist  = ", Rins_hist(Nt)
    !print *, " fsat_hist  = ", fsat_hist(Nt)
    
    ! ============================================================
    !  Energy diagnostics
    ! ============================================================
    call compute_energy_diagnostics(t_vec, Nt, qa_hist, qsolar_hist, &
                                qrain_hist, qevap_hist, qground_hist, &
                                T_hist, Cs_layer, E_melt, rho_i, Lf, &
                                USE_ADVANCED_INSULATION, SOLVER_CHOICE)


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

double precision function Isolar_fun(t) result(I_solar)
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

double precision function Prain_fun(t) result(Prain)
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
double precision function ground_flux(T3, Tg, R_3g) result(flux)
    ! Ground -> bottom snow layer flux [W/m^2].
    implicit none
    double precision, intent(in) :: T3, Tg, R_3g

    flux = (Tg - T3) / R_3g

end function ground_flux

! ============================================================
!  dT/dt for snow layers
! ============================================================
subroutine dTdt(t, Tv, R_a2s, q_solar, q_rain, q_evap, &
                        R_nm, Cs_layer, Tg, dT)
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
    double precision, dimension(3), intent(out) :: dT

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

end subroutine dTdt

! ============================================================
!  Insulation step (advanced model)
! ============================================================
subroutine insulation_step(state_in, forc, p, delta_t, R_ins, q_solar, &
                        q_rain_snow, q_evap, state_out)
    ! Advanced insulation model with moisture, age, rain, and evaporation
    use types_module
    implicit none
    
    ! Inputs
    type(InsulationState), intent(inout) :: state_in
    type(Forcing), intent(in) :: forc
    type(InsulationParameters), intent(in) :: p
    double precision, intent(in) :: delta_t
    
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
    W_new = W + delta_t * dWdt
    W_new = max(0.0D0, min(p%W_sat, W_new))  ! clip to [0, W_sat]
    
    age_days_new = age_days + delta_t / 86400.0D0
    
    ! Update output state
    state_out = state_in  ! copy all fields
    state_out%W = W_new
    state_out%age_days = age_days_new
    
    ! Store diagnostics
    state_out%k_eff = k_eff
    state_out%alpha_eff = alpha_eff
    state_out%f_sat = f
    
end subroutine insulation_step

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

FUNCTION romberg_discrete(y, x, n, max_order) RESULT(integral)
    !======================================================================
    ! Romberg integration for discrete uniformly-spaced data
    ! Uses Richardson extrapolation on trapezoidal rule
    ! 
    ! y: array of function values
    ! x: array of corresponding x values (must be uniformly spaced)
    ! n: size of arrays
    ! max_order: maximum order of extrapolation (typically 4-6)
    !======================================================================
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, max_order
    DOUBLE PRECISION, DIMENSION(n), INTENT(IN) :: y, x
    DOUBLE PRECISION :: integral
    
    INTEGER, PARAMETER :: MAXORDER = 10
    DOUBLE PRECISION, DIMENSION(MAXORDER, MAXORDER) :: R
    INTEGER :: i, j, k, m, step_size
    DOUBLE PRECISION :: h, sum_val
    INTEGER :: actual_order
    
    ! Limit order based on available points
    actual_order = min(max_order, MAXORDER)
    actual_order = min(actual_order, int(log(real(n-1))/log(2.0)) + 1)
    
    if (n < 2) then
        integral = 0.0D0
        return
    end if
    
    h = x(2) - x(1)  ! Assume uniform spacing
    
    ! R(1,1): Trapezoidal rule with all points
    R(1, 1) = 0.5D0 * h * (y(1) + y(n))
    do i = 2, n-1
        R(1, 1) = R(1, 1) + h * y(i)
    end do
    
    ! Build Romberg table using subset of points
    do i = 2, actual_order
        ! Use every 2^(i-1) points
        step_size = 2**(i-1)
        
        if (1 + (n-1)/step_size < 2) exit  ! Not enough points
        
        sum_val = 0.0D0
        m = 0
        do k = 1, n, step_size
            if (k == 1 .or. k == n) then
                sum_val = sum_val + 0.5D0 * y(k)
            else
                sum_val = sum_val + y(k)
            end if
            m = m + 1
        end do
        
        R(i, 1) = sum_val * h * step_size
        
        ! Richardson extrapolation
        do j = 2, i
            R(i, j) = R(i, j-1) + (R(i, j-1) - R(i-1, j-1)) / (4.0D0**(j-1) - 1.0D0)
        end do
    end do
    
    ! Return the highest order estimate
    integral = R(actual_order, actual_order)
END FUNCTION romberg_discrete

! ============================================================
!  Energy diagnostics section (replaces Python trapz calls)
! ============================================================
SUBROUTINE compute_energy_diagnostics(t_vec, Nt, qa_hist, qsolar_hist, &
                                     qrain_hist, qevap_hist, qground_hist, &
                                     T_hist, Cs_layer, E_melt, rho_i, Lf, &
                                     USE_ADVANCED, SOLVER_CHOICE)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: Nt, SOLVER_CHOICE
    DOUBLE PRECISION, DIMENSION(Nt), INTENT(IN) :: t_vec, qa_hist, qsolar_hist
    DOUBLE PRECISION, DIMENSION(Nt), INTENT(IN) :: qrain_hist, qevap_hist, qground_hist
    DOUBLE PRECISION, DIMENSION(Nt, 3), INTENT(IN) :: T_hist
    DOUBLE PRECISION, INTENT(IN) :: Cs_layer, E_melt, rho_i, Lf
    LOGICAL, INTENT(IN) :: USE_ADVANCED
    
    DOUBLE PRECISION :: E_a, E_solar, E_rain, E_evap, E_g
    DOUBLE PRECISION :: E_total_in, E_snow_change, E_balance, M_melt
    
    ! External functions
    DOUBLE PRECISION, EXTERNAL :: simpson_discrete
    
    ! Integrate energy fluxes using Simpson's rule
    E_a     = simpson_discrete(qa_hist, t_vec, Nt)
    E_solar = simpson_discrete(qsolar_hist, t_vec, Nt)
    E_rain  = simpson_discrete(qrain_hist, t_vec, Nt)
    E_evap  = simpson_discrete(qevap_hist, t_vec, Nt)
    E_g     = simpson_discrete(qground_hist, t_vec, Nt)
    
    E_total_in = E_a + E_solar + E_rain + E_evap + E_g
    
    ! Change in snow internal energy (all 3 layers)
    E_snow_change = Cs_layer * (sum(T_hist(Nt, :)) - sum(T_hist(1, :)))
    
    E_balance = E_total_in - (E_snow_change + E_melt)
    M_melt = E_melt / (rho_i * Lf)
    
    ! Print results
    print *, ""
    print *, "============================================================"
    SELECT CASE (SOLVER_CHOICE)
    CASE (1)
        print *, "Energy Diagnostics (RK4)"
    CASE (2)
        print *, "Energy Diagnostics (DOPRI5)"
    CASE (3)
        print *, "Energy Diagnostics (LSODA)"
    END SELECT

    print *, "============================================================"
    print '(A, L)', " USE_ADVANCED_INSULATION = ", USE_ADVANCED
    print *, ""
    print '(A, ES14.6, A)', " E_air           = ", E_a, " J/m²"
    print '(A, ES14.6, A)', " E_solar         = ", E_solar, " J/m²"
    print '(A, ES14.6, A)', " E_rain          = ", E_rain, " J/m²"
    print '(A, ES14.6, A)', " E_evap          = ", E_evap, " J/m²"
    print '(A, ES14.6, A)', " E_ground        = ", E_g, " J/m²"
    print '(A, ES14.6, A)', " E_total_in      = ", E_total_in, " J/m²"
    print *, ""
    print '(A, ES14.6, A)', " E_snow_change   = ", E_snow_change, " J/m²"
    print '(A, ES14.6, A)', " E_melt          = ", E_melt, " J/m²"
    print '(A, ES14.6, A)', " Energy residual = ", E_balance, " J/m²"
    print *, ""
    print '(A, F10.6, A)', " Total melted snow thickness = ", M_melt, " m"
    print *, "============================================================"
END SUBROUTINE compute_energy_diagnostics

subroutine set_forcing(t_mid, forc)
    use types_module
    implicit none
    
    double precision :: t_mid
    type(Forcing) :: forc
    double precision, external :: Ta_fun, Isolar_fun, Prain_fun
    
    forc%Isolar = Isolar_fun(t_mid)
    forc%Prain  = Prain_fun(t_mid)
    forc%T_rain = 273.15D0 + 3.0D0
    forc%RH     = 0.80D0
    forc%Ta     = Ta_fun(t_mid)
end subroutine set_forcing

SUBROUTINE get_insulation_properties(USE_ADVANCED, InsState, forc, InsPar, &
                                     dt, Hi, k_i_const, alpha_const, eta_rain_const, &
                                     rho_w, c_w, Tfreeze, &
                                     R_ins, q_solar, q_rain_snow, q_evap)
    use types_module
    IMPLICIT NONE
    LOGICAL :: USE_ADVANCED
    TYPE(InsulationState) :: InsState
    TYPE(Forcing) :: forc
    TYPE(InsulationParameters) :: InsPar
    DOUBLE PRECISION :: dt, Hi, k_i_const, alpha_const, eta_rain_const
    DOUBLE PRECISION :: rho_w, c_w, Tfreeze
    DOUBLE PRECISION :: R_ins, q_solar, q_rain_snow, q_evap
    
    if (USE_ADVANCED) then
        call insulation_step(InsState, forc, InsPar, dt, &
                           R_ins, q_solar, q_rain_snow, q_evap, InsState)
    else
        R_ins = Hi / k_i_const
        q_solar = alpha_const * forc%Isolar
        q_rain_snow = eta_rain_const * rho_w * c_w * &
                     forc%Prain * (forc%T_rain - Tfreeze)
        q_evap = 0.0D0
    end if
END SUBROUTINE get_insulation_properties

SUBROUTINE step_RK4(t, dt, T_current, R_eff, R_ins, q_solar, q_rain, &
                   q_evap, R_nm, Cs_layer, Tg, T_new)
    IMPLICIT NONE
    DOUBLE PRECISION :: t, dt, R_eff, R_ins, q_solar, q_rain, q_evap
    DOUBLE PRECISION :: Cs_layer, Tg, R_a2s
    DOUBLE PRECISION, DIMENSION(3) :: T_current, R_nm, T_new
    DOUBLE PRECISION, DIMENSION(3) :: k1_vec, k2_vec, k3_vec, k4_vec
    
    R_a2s = R_eff + R_ins
    
    call dTdt(t, T_current, R_a2s, q_solar, q_rain, q_evap, &
             R_nm, Cs_layer, Tg, k1_vec)
    call dTdt(t + dt/2.0D0, T_current + dt*k1_vec/2.0D0, &
             R_a2s, q_solar, q_rain, q_evap, R_nm, Cs_layer, Tg, k2_vec)
    call dTdt(t + dt/2.0D0, T_current + dt*k2_vec/2.0D0, &
             R_a2s, q_solar, q_rain, q_evap, R_nm, Cs_layer, Tg, k3_vec)
    call dTdt(t + dt, T_current + dt*k3_vec, &
             R_a2s, q_solar, q_rain, q_evap, R_nm, Cs_layer, Tg, k4_vec)
    
    T_new = T_current + (dt/6.0D0) * (k1_vec + 2.0D0*k2_vec + 2.0D0*k3_vec + k4_vec)
END SUBROUTINE step_RK4

SUBROUTINE process_surface_energy(t_mid, T_old, T_new, R_eff, R_ins, &
                                  q_solar, q_rain, q_evap, Tg, R_3g, &
                                  Tfreeze, dt, rho_i, Lf, &
                                  q_a, q_surf, q_ground, E_melt, melt_rate)
    IMPLICIT NONE
    DOUBLE PRECISION :: t_mid, R_eff, R_ins, q_solar, q_rain, q_evap
    DOUBLE PRECISION :: Tg, R_3g, Tfreeze, dt, rho_i, Lf
    DOUBLE PRECISION :: q_a, q_surf, q_ground, E_melt, melt_rate
    DOUBLE PRECISION, DIMENSION(3) :: T_old, T_new
    DOUBLE PRECISION :: Ta_mid, T1_mid, R_a2s, dE_melt, dM_melt
    DOUBLE PRECISION, EXTERNAL :: Ta_fun, ground_flux
    
    R_a2s = R_eff + R_ins
    Ta_mid = Ta_fun(t_mid)
    T1_mid = T_old(1)
    
    q_a = (Ta_mid - T1_mid) / R_a2s
    q_surf = q_a + q_solar + q_rain + q_evap
    q_ground = ground_flux(T_old(3), Tg, R_3g)
    
    ! Apply melting clamp
    if (T_new(1) > Tfreeze) then
        T_new(1) = Tfreeze
        if (q_surf > 0.0D0) then
            dE_melt = q_surf * dt
            dM_melt = dE_melt / (rho_i * Lf)
        else
            dE_melt = 0.0D0
            dM_melt = 0.0D0
        end if
        E_melt = E_melt + dE_melt
        melt_rate = dM_melt / dt
    else
        melt_rate = 0.0D0
    end if
END SUBROUTINE process_surface_energy

! ============================================================
!  INTEGRATION OPTIONS: DOPRI5 and LSODA
!  Place these subroutines in the CONTAINS section
! ============================================================

! ------------------------------------------------------------
!  OPTION 1: DOPRI5 (Adaptive RK, explicit, no Jacobian)
! ------------------------------------------------------------
SUBROUTINE integrate_DOPRI5(t, dt, T_current, R_ins, R_eff, q_solar, q_rain, &
                            q_evap, R_nm, Cs_layer, Tg, T_new)
    IMPLICIT NONE
    DOUBLE PRECISION :: t, dt, R_ins, R_eff, q_solar, q_rain, q_evap
    DOUBLE PRECISION :: Cs_layer, Tg
    DOUBLE PRECISION, DIMENSION(3) :: T_current, R_nm, T_new
    
    ! DOPRI5 parameters
    INTEGER, PARAMETER :: N = 3              ! System dimension
    INTEGER, PARAMETER :: LWORK = 8*N + 21
    INTEGER, PARAMETER :: LIWORK = 21
    DOUBLE PRECISION :: WORK(LWORK), RPAR(20)
    INTEGER :: IWORK(LIWORK), ITOL, IOUT, IDID, IPAR
    DOUBLE PRECISION :: RTOL, ATOL, XEND
    DOUBLE PRECISION :: Y(N)
    
    EXTERNAL :: snow_rhs_dopri5, solout_dummy
    
    ! Pack parameters into RPAR
    RPAR(1)  = R_eff
    RPAR(2)  = R_ins
    RPAR(3)  = q_solar
    RPAR(4)  = q_rain
    RPAR(5)  = q_evap
    RPAR(6)  = R_nm(1)  ! R_12
    RPAR(7)  = R_nm(2)  ! R_23
    RPAR(8)  = R_nm(3)  ! R_3g
    RPAR(9)  = Cs_layer
    RPAR(10) = Tg
    
    ! Initial conditions
    Y = T_current
    
    ! Solver settings
    ITOL = 0                    ! Scalar tolerances
    RTOL = 1.0D-4              ! Relative tolerance
    ATOL = 1.0D-3              ! Absolute tolerance (0.001 K)
    IOUT = 0                   ! No dense output needed
    XEND = t + dt              ! End time
    
    ! Clear work arrays
    WORK = 0.0D0
    IWORK = 0
    
    ! Call DOPRI5
    CALL DOPRI5(N, snow_rhs_dopri5, t, Y, XEND, &
                RTOL, ATOL, ITOL, &
                solout_dummy, IOUT, &
                WORK, LWORK, IWORK, LIWORK, RPAR, IPAR, IDID)
    
    ! Check for errors
    IF (IDID < 0) THEN
        WRITE(*,'(A,I3)') 'DOPRI5 error: IDID = ', IDID
        STOP
    END IF
    
    ! Return solution
    T_new = Y
END SUBROUTINE integrate_DOPRI5

! Right-hand side for DOPRI5
SUBROUTINE snow_rhs_dopri5(N, T_TIME, Y, F, RPAR, IPAR)
    IMPLICIT NONE
    INTEGER :: N, IPAR
    DOUBLE PRECISION :: T_TIME, Y(N), F(N), RPAR(*)
    DOUBLE PRECISION :: R_eff, R_ins, R_a2s, q_solar, q_rain, q_evap
    DOUBLE PRECISION :: R_12, R_23, R_3g, Cs_layer, Tg
    DOUBLE PRECISION :: T1, T2, T3, Ta, q_a, q_surf
    DOUBLE PRECISION :: q_12, q_21, q_23, q_32, q_3g
    DOUBLE PRECISION, EXTERNAL :: Ta_fun
    
    ! Unpack parameters
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
    
    R_a2s = R_eff + R_ins
    
    T1 = Y(1)
    T2 = Y(2)
    T3 = Y(3)
    Ta = Ta_fun(T_TIME)
    
    ! Surface energy balance
    q_a = (Ta - T1) / R_a2s
    q_surf = q_a + q_solar + q_rain + q_evap
    
    ! Layer 1 (surface)
    q_12 = (T2 - T1) / R_12
    F(1) = (q_surf + q_12) / Cs_layer
    
    ! Layer 2 (middle)
    q_21 = (T1 - T2) / R_12
    q_23 = (T3 - T2) / R_23
    F(2) = (q_21 + q_23) / Cs_layer
    
    ! Layer 3 (bottom)
    q_32 = (T2 - T3) / R_23
    q_3g = (Tg - T3) / R_3g
    F(3) = (q_32 + q_3g) / Cs_layer
    
    ! Avoid unused variable warning
    IF (IPAR < -999) WRITE(*,*) IPAR
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
                           q_evap, R_nm, Cs_layer, Tg, T_new)
    USE odepack_interface
    USE odepack_common
    IMPLICIT NONE
    
    DOUBLE PRECISION :: t, dt, R_ins, R_eff, q_solar, q_rain, q_evap
    DOUBLE PRECISION :: Cs_layer, Tg
    DOUBLE PRECISION, DIMENSION(3) :: T_current, R_nm, T_new
    
    ! LSODA parameters
    INTEGER, PARAMETER :: NEQ = 3
    INTEGER :: ITOL, ITASK, ISTATE, IOPT, JT, LRW, LIW, IPAR
    DOUBLE PRECISION :: RTOL, ATOL(NEQ), TOUT
    DOUBLE PRECISION :: Y(NEQ), RPAR(20)
    INTEGER, ALLOCATABLE :: IWORK(:)
    DOUBLE PRECISION, ALLOCATABLE :: RWORK(:)
    TYPE(odepack_common_data), TARGET :: common_data
    
    EXTERNAL :: snow_rhs_lsoda, jac_dummy
    
    ! Calculate work array sizes for LSODA
    LRW = MAX(20 + 16*NEQ, 22 + 9*NEQ + NEQ*NEQ)
    LIW = 20 + NEQ
    
    ALLOCATE(RWORK(LRW), IWORK(LIW))
    
    ! Pack parameters into RPAR
    RPAR(1)  = R_eff
    RPAR(2)  = R_ins
    RPAR(3)  = q_solar
    RPAR(4)  = q_rain
    RPAR(5)  = q_evap
    RPAR(6)  = R_nm(1)  ! R_12
    RPAR(7)  = R_nm(2)  ! R_23
    RPAR(8)  = R_nm(3)  ! R_3g
    RPAR(9)  = Cs_layer
    RPAR(10) = Tg
    
    ! Initial conditions
    Y = T_current
    
    ! Initialize arrays
    IWORK = 0
    RWORK = 0.0D0
    common_data%ierr = 0
    
    ! Solver settings
    ITOL = 2                    ! Array of absolute tolerances
    RTOL = 1.0D-5              ! Relative tolerance
    ATOL = 1.0D-4              ! Absolute tolerance (0.0001 K)
    ITASK = 1                  ! Normal computation to TOUT
    ISTATE = 1                 ! First call
    IOPT = 0                   ! No optional inputs
    JT = 2                     ! Internally generated Jacobian
    IPAR = 0
    
    TOUT = t + dt
    
    ! Call DLSODA
    CALL DLSODA(snow_rhs_lsoda, NEQ, Y, t, TOUT, ITOL, RTOL, ATOL, &
                ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, &
                jac_dummy, JT, common_data)
    
    ! Check for errors
    IF (ISTATE < 0) THEN
        WRITE(*,'(A,I3)') 'DLSODA error: ISTATE = ', ISTATE
        SELECT CASE (ISTATE)
        CASE (-1)
            WRITE(*,'(A)') 'Excess work done. Try increasing MXSTEP.'
        CASE (-2)
            WRITE(*,'(A)') 'Excess accuracy requested.'
        CASE (-3)
            WRITE(*,'(A)') 'Illegal input detected.'
        CASE (-4)
            WRITE(*,'(A)') 'Repeated error test failures.'
        CASE (-5)
            WRITE(*,'(A)') 'Repeated convergence failures.'
        CASE (-6)
            WRITE(*,'(A)') 'Error weight became zero.'
        END SELECT
        STOP
    END IF
    
    ! Return solution
    T_new = Y
    
    DEALLOCATE(RWORK, IWORK)
END SUBROUTINE integrate_LSODA

! Right-hand side for LSODA
SUBROUTINE snow_rhs_lsoda(NEQ, T_TIME, Y, YDOT, common_data)
    USE odepack_common
    IMPLICIT NONE
    INTEGER :: NEQ
    DOUBLE PRECISION :: T_TIME, Y(NEQ), YDOT(NEQ)
    TYPE(odepack_common_data) :: common_data
    
    ! Note: RPAR needs to be passed via common block or module for LSODA
    ! For simplicity, we'll access global variables from the main program
    ! In production code, use a module to store these parameters
    
    DOUBLE PRECISION :: R_eff, R_ins, R_a2s, q_solar, q_rain, q_evap
    DOUBLE PRECISION :: R_12, R_23, R_3g, Cs_layer, Tg
    DOUBLE PRECISION :: T1, T2, T3, Ta, q_a, q_surf
    DOUBLE PRECISION :: q_12, q_21, q_23, q_32, q_3g
    DOUBLE PRECISION, EXTERNAL :: Ta_fun
    COMMON /SNOW_PARAMS/ R_eff, R_ins, q_solar, q_rain, q_evap, &
                        R_12, R_23, R_3g, Cs_layer, Tg
    
    R_a2s = R_eff + R_ins
    
    T1 = Y(1)
    T2 = Y(2)
    T3 = Y(3)
    Ta = Ta_fun(T_TIME)
    
    ! Surface energy balance
    q_a = (Ta - T1) / R_a2s
    q_surf = q_a + q_solar + q_rain + q_evap
    
    ! Layer 1 (surface)
    q_12 = (T2 - T1) / R_12
    YDOT(1) = (q_surf + q_12) / Cs_layer
    
    ! Layer 2 (middle)
    q_21 = (T1 - T2) / R_12
    q_23 = (T3 - T2) / R_23
    YDOT(2) = (q_21 + q_23) / Cs_layer
    
    ! Layer 3 (bottom)
    q_32 = (T2 - T3) / R_23
    q_3g = (Tg - T3) / R_3g
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
