import numpy as np
import matplotlib.pyplot as plt
from numba import njit
import csv
from datetime import datetime

# ============================================================
#  Enhanced Multi-layer Snow Storage RC Model with Real Data
#  Improvements:
#  1. Real meteorological data from CSV
#  2. Multi-layer insulation option
#  3. Refreezing with cold content
#  4. Percolation with bucket method
#  5. Advanced insulation model with moisture and age
# ------------------------------------------------------------
# New features:
#   * Fixed huge bottleneck of multilayer insulation solver 
#     using Numba JIT compilation.
#   * Function for multilayer insulation now comes with transient 
#     Robin BCs at both boundaries. Outer boundary is dynamic 
#     based on wind and environment. Inner boundary is static 
#     based on snow contact.
#   * Bottom ground flux now uses Robin BC with ground HTC.
# ============================================================

USE_ADVANCED_INSULATION = True
USE_REFREEZING = True
USE_PERCOLATION = True
USE_MULTILAYER_INSULATION = True

# ---------- Physical constants ----------
sigma   = 5.670374419e-8      # Stefan-Boltzmann [W/m^2 K^4]
Lf      = 3.34e5              # Latent heat of fusion [J/kg]
rho_i   = 917.0               # Ice density [kg/m^3]
rho_s   = 400.0               # Snow bulk density [kg/m^3]
c_s     = 2100.0              # Snow specific heat [J/kg K]
rho_w   = 1000.0              # Water density [kg/m^3]
c_w     = 4180.0              # Water specific heat [J/kg K]
Tfreeze = 273.15              # 0°C [K]

# ---------- Snow & insulation ----------
Hs   = 2.0                    # total snow thickness [m]
Ns   = 3                      # number of snow layers
dz_s = Hs / Ns                # thickness per snow layer [m]

k_snow = 0.25                 # snow conductivity [W/mK]
Hi     = 0.20                 # insulation thickness [m]

# Multi-layer insulation setup
if USE_MULTILAYER_INSULATION:
    N_ins = 20                # number of insulation layers
    dz_ins = Hi / N_ins * 0.1 # thickness per insulation layer [m]
else:
    N_ins = 1
    dz_ins = Hi

# Insulation material properties
k_i_base = 0.32               # base thermal conductivity [W/mK]
rho_dry = 100.0               # dry density [kg/m^3]
moist_cont = 50.0             # moisture content [%]
rho_wet = rho_dry + moist_cont/100.0*1000  # wet density [kg/m^3]
c_dry = 0.99e3                # dry specific heat [J/(kg*K)]
c_wet = (1.0 - moist_cont/100.0)*c_dry + moist_cont/100.0*c_w  # wet specific heat
D_ins = k_i_base / (c_wet * rho_wet)  # thermal diffusivity [m^2/s]

# Ground boundary condition - Robin BC
h_ground = 3.0  # Ground heat transfer coefficient [W/m²K]
                # Typical range: 2-5 W/(m²K) for soil interface
                # Reference: NREL/TP-550-33954 (Deru, 2003)
                # Lower values = better insulated ground
                # Higher values = more conductive ground/higher water table

Tg_deep = 273.15 + 2.0  # Deep ground temperature [K]
                         # This represents the undisturbed ground temperature
                         # at sufficient depth (typically 2-3m)

# Simple (constant) insulation parameters
alpha_const   = 0.80          # solar absorptivity
eta_rain_const = 1.0          # fraction of rain heat reaching snow

# Ground insulation
Hg_ins = 0.3                  # [m]
kg_ins = 0.04                 # [W/mK]

# ---------- Surface HTC (air-side, conv + LW) ----------
h_conv  = 8.0                 # convective coefficient [W/m^2K]
epsilon = 0.95
T_mean  = 273.15 + 3.0        # nominal mean temperature [K]
h_rad   = 4.0 * epsilon * sigma * T_mean**3
h_eff   = h_conv + h_rad

R_eff = 1.0 / h_eff           # air-side resistance [m^2K/W]

# ---------- Thermal resistances in snow & ground ----------
R_layer = dz_s / k_snow
R_12    = R_layer
R_23    = R_layer

R_g_ins = Hg_ins / kg_ins
R_3g    = R_layer + R_g_ins

# ---------- Ground & snow capacity ----------
Tg        = 273.15 + 2.0
Cs_layer  = rho_s * c_s * dz_s   # [J/(m^2 K)] per snow layer

# ---------- Initial snow temperatures ----------
T1_init = 273.15 - 2.0
T2_init = 273.15 - 4.0
T3_init = 273.15 - 6.0
T = np.array([T1_init, T2_init, T3_init], dtype=float)

# ---------- Layer properties for refreezing/percolation ----------
LWC = np.array([0.0, 0.0, 0.0], dtype=float)
theta_e = 0.04  # Irreducible water content

# Ice fractions (assuming constant for now)
ice_fractions = np.array([0.4, 0.4, 0.4])
heights = np.array([dz_s, dz_s, dz_s])

# ============================================================
#  CSV Data Loading Functions
# ============================================================
def read_csv_data(filename):
    """Read meteorological data from CSV file."""
    data = {
        'time': [],
        'temp': [],      # Air temperature [°C]
        'wind': [],      # Wind speed [m/s]
        'precip': [],    # Precipitation [m/h]
        'solar': [],     # Global solar irradiance [W/m²]
        'rh': []         # Relative humidity [%]
    }
    
    with open(filename, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                data['time'].append(row['Time'])
                data['temp'].append(float(row['Temp_C']))
                data['wind'].append(float(row['Air_Vel_m/s']))
                data['precip'].append(float(row['Prec_m/h']))
                data['solar'].append(float(row['Glo_Sol_Ir_W/m2']))
                data['rh'].append(float(row['RH_%']))
            except (ValueError, KeyError) as e:
                print(f"Warning: Skipping row due to error: {e}")
                continue
    
    return data

def interpolate_data(data_vec, t_query, dt_data=3600.0):
    """
    Interpolate hourly data to arbitrary time points.
    
    Args:
        data_vec: List of data values (hourly)
        t_query: Query time [s]
        dt_data: Time step of input data [s] (default 3600s = 1 hour)
    
    Returns:
        Interpolated value
    """
    idx = t_query / dt_data
    idx_low = int(np.floor(idx))
    idx_high = int(np.ceil(idx))
    
    # Boundary handling
    if idx_low < 0:
        return data_vec[0]
    if idx_high >= len(data_vec):
        return data_vec[-1]
    if idx_low == idx_high:
        return data_vec[idx_low]
    
    # Linear interpolation
    frac = idx - idx_low
    return data_vec[idx_low] * (1 - frac) + data_vec[idx_high] * frac

# ============================================================
#  Advanced insulation parameters
# ============================================================
if USE_ADVANCED_INSULATION:
    InsPar = {
        "Hi":     Hi,
        "k_dry":  0.06,
        "k_sat":  0.30,
        "n_k":    1.5,

        "W_sat":   30.0,
        "W_field": 10.0,

        "alpha_dry": 0.10,
        "alpha_wet": 0.25,
        "n_alpha":   1.0,

        "delta_k_age":     0.5,
        "tau_k_years":     2.0,
        "delta_alpha_age": 0.05,
        "tau_alpha_years": 2.0,

        "zeta0":   0.3,
        "gamma_H": 0.5,
        "gamma_W": 2.0,

        "beta_w": 3.0,
        "K_E":    1e-5,
        "K_D":    5e-6,

        "Lv":      2.5e6,
        "rho_w":   rho_w,
        "c_w":     c_w,
        "Tfreeze": Tfreeze,

        "rho_air": 1.2,
        "C_E":     1.3e-3,
        "U10":     2.0,
        "P0":      101325.0
    }

    InsState = {
        "W":        5.0,
        "age_days": 0.0
    }
else:
    InsPar = None
    InsState = {
        "W":        0.0,
        "age_days": 0.0
    }

# ============================================================
#  Multi-layer insulation thermal solver (TDMA)
# ============================================================
@njit
def solve_tdma(a, b, c, d, n):
    """Tridiagonal matrix algorithm (Thomas algorithm)."""
    c_prime = np.zeros(n)
    d_prime = np.zeros(n)
    x = np.zeros(n)

    c_prime[0] = c[0] / b[0]
    d_prime[0] = d[0] / b[0]

    for i in range(1, n):
        denom = b[i] - a[i] * c_prime[i - 1]
        c_prime[i] = c[i] / denom
        d_prime[i] = (d[i] - a[i] * d_prime[i - 1]) / denom

    x[-1] = d_prime[-1]
    for i in range(n - 2, -1, -1):
        x[i] = d_prime[i] - c_prime[i] * x[i + 1]

    return x

@njit
def compute_insulation_resistance_multilayer_core(T_env, h_out, T_snow, h_in, k_eff, dt_main, T_n_init, dt_substep=10.0):
    """
    Core computation for insulation temperature profile using Robin Boundary Conditions.
    This version takes the initial temperature profile as input and returns the updated one.
    
    Args:
        T_env: Equivalent environmental temperature (Sol-Air temp) [K]
        h_out: External Heat Transfer Coefficient [W/m²K] (from wind speed)
        T_snow: Temperature of the snow directly below insulation [K]
        h_in: Contact conductance/HTC between insulation and snow [W/m²K]
        k_eff: Effective thermal conductivity of insulation [W/mK]
        dt_main: Main simulation time step [s] (required, typically 600s)
        T_n_init: Initial temperature profile [K] (N_ins + 1 nodes)
        dt_substep: Sub-time step for the thermal solver [s] (default 10s)
        
    Returns: 
        R_ins: Effective resistance (calculated from resulting flux)
        T_n: New temperature profile
    """
    # Calculate number of substeps to cover one main time step
    n_substeps = int(dt_main / dt_substep)
    
    dx = dz_ins  # Cell size [m]
    
    T_n = T_n_init.copy()
    
    # Arrays for TDMA
    N_nodes = N_ins + 1
    a = np.zeros(N_nodes)
    b = np.zeros(N_nodes)
    c = np.zeros(N_nodes)
    d = np.zeros(N_nodes)
    
    # Time-stepping loop
    for _ in range(n_substeps):
        # Fourier number for this substep
        dFo = D_ins * dt_substep / dx**2
        
        # Biot numbers (dimensionless)
        dBio_o = h_out * dx / k_eff
        dBio_i = h_in * dx / k_eff
        
        # ============================================================
        # Outer boundary (node 0): Robin BC with external environment
        # ============================================================
        b[0] = (1.0 + 2.0 * dFo + 2.0 * dFo * dBio_o)
        c[0] = -2.0 * dFo
        d[0] = T_n[0] + 2.0 * dFo * dBio_o * T_env
        
        # ============================================================
        # Interior nodes: Standard implicit scheme
        # ============================================================
        for j in range(1, N_nodes - 1):
            a[j] = -dFo
            b[j] = (1.0 + 2.0 * dFo)
            c[j] = -dFo
            d[j] = T_n[j]
        
        # ============================================================
        # Inner boundary (node -1): Robin BC with snow contact
        # ============================================================
        a[-1] = -2.0 * dFo
        b[-1] = (1.0 + 2.0 * dFo + 2.0 * dFo * dBio_i)
        c[-1] = 0.0  # No node beyond this
        d[-1] = T_n[-1] + 2.0 * dFo * dBio_i * T_snow
        
        # Solve the tridiagonal system
        T_n = solve_tdma(a, b, c, d, N_nodes)
    
    # ============================================================
    # Calculate effective resistance from heat flux
    # ============================================================
    # Heat flux at the snow interface (inner boundary)
    # q = h_in * (T_ins_bottom - T_snow)
    q_into_snow = h_in * (T_n[-1] - T_snow)
    
    # Effective resistance: R = ΔT / q
    delta_T_total = T_env - T_snow
    
    if abs(q_into_snow) > 1e-9 and abs(delta_T_total) > 0.01:
        R_ins = delta_T_total / q_into_snow
    else:
        # Fallback to simple conduction if flux is very small
        R_ins = Hi / k_eff
    
    return R_ins, T_n


def compute_insulation_resistance_multilayer(T_env, h_out, T_snow, h_in, k_eff, dt_main, dt_substep=10.0):
    """
    Wrapper function that manages state for the insulation temperature profile.
    """
    global T_ins_profile
    
    if not USE_MULTILAYER_INSULATION:
        return Hi / k_eff, None
    
    # Initialize profile if first call
    if T_ins_profile is None:
        T_ins_profile = np.linspace(T_env, T_snow, N_ins + 1)
    
    # Call the JIT-compiled core function
    R_ins, T_ins_profile = compute_insulation_resistance_multilayer_core(
        T_env, h_out, T_snow, h_in, k_eff, dt_main, T_ins_profile, dt_substep
    )
    
    return R_ins, T_ins_profile

# ============================================================
#  Refreezing function
# ============================================================
@njit
def refreezing_layer(T_layer, LWC_layer, ice_frac):
    """Refreeze water in a single layer based on cold content."""
    if (T_layer >= Tfreeze) or (LWC_layer <= 0.0):
        return T_layer, LWC_layer, ice_frac, 0.0
    
    dT_max = T_layer - Tfreeze
    dtheta_w_max = -(dT_max * (ice_frac * rho_i * c_s + LWC_layer * rho_w * c_w)) / (rho_w * Lf)
    dtheta_w = min(LWC_layer, dtheta_w_max)
    dtheta_i = (rho_w / rho_i) * dtheta_w
    dT = (dtheta_w * rho_w * Lf) / (ice_frac * rho_i * c_s + LWC_layer * rho_w * c_w)
    
    new_T = T_layer + dT
    new_LWC = LWC_layer - dtheta_w
    new_ice_frac = ice_frac + dtheta_i
    refrozen_mass = dtheta_w * dz_s * rho_w
    
    return new_T, new_LWC, new_ice_frac, refrozen_mass

# ============================================================
#  Percolation function
# ============================================================
@njit
def percolate_water(LWC_array, heights, theta_e):
    """Percolate water through layers using bucket method."""
    n_layers = len(LWC_array)
    new_LWC = LWC_array.copy()
    
    for i in range(n_layers - 1):
        if new_LWC[i] > theta_e:
            excess = new_LWC[i] - theta_e
            new_LWC[i] = theta_e
            excess_mass = excess * heights[i]
            new_LWC[i+1] += excess_mass / heights[i+1]
    
    runoff = 0.0
    if new_LWC[n_layers-1] > theta_e:
        excess = new_LWC[n_layers-1] - theta_e
        new_LWC[n_layers-1] = theta_e
        runoff = excess * heights[n_layers-1] * rho_w
    
    return new_LWC, runoff

# ============================================================
#  Ground flux
# ============================================================
def ground_flux(T3):
    """Ground -> bottom snow layer flux [W/m^2]."""
    return (Tg - T3) / R_3g

def ground_flux_robin_bc(T3, Tg_deep, h_ground, k_soil=1.5, L_soil=1.0):
    """
    Ground flux with combined conduction and interface resistance.
    
    Args:
        T3: Temperature of bottom snow layer [K]
        Tg_deep: Deep ground temperature [K]
        h_ground: Ground interface heat transfer coefficient [W/m²K]
        k_soil: Soil thermal conductivity [W/mK] (typical: 0.5-2.5)
        L_soil: Effective soil layer thickness [m] (typical: 0.5-2.0)
    
    Returns:
        q_ground: Heat flux from ground to snow [W/m²]

    Links:
        - https://www.cableizer.com/documentation/h_tr/
        - https://docs.nrel.gov/docs/fy03osti/33954.pdf
    """
    # Combined resistance: R_total = R_conduction + R_interface
    R_cond = L_soil / k_soil  # Conduction resistance through soil
    R_interface = 1.0 / h_ground  # Interface resistance
    R_total = R_cond + R_interface
    
    q_ground = (Tg_deep - T3) / R_total
    return q_ground

# ============================================================
#  dT/dt for snow layers
# ============================================================
def dTdt(t, Tv, stepPar, met_data, dt_data):
    """
    Tv = [T1, T2, T3]
    stepPar: dict with R_a2s, q_solar, q_rain, q_evap
    """
    T1, T2, T3 = Tv
    
    # Get Ta at the specific time point
    Ta = interpolate_data(met_data['temp'], t, dt_data) + 273.15
    
    q_a = (Ta - T1) / stepPar["R_a2s"]
    q_surf = q_a + stepPar["q_solar"] + stepPar["q_rain"] + stepPar["q_evap"]

    q_12 = (T2 - T1) / R_12
    dT1  = (q_surf + q_12) / Cs_layer

    q_21 = (T1 - T2) / R_12
    q_23 = (T3 - T2) / R_23
    dT2  = (q_21 + q_23) / Cs_layer

    q_32 = (T2 - T3) / R_23
    #q_3g = (Tg - T3) / R_3g # Fixed ground flux
    q_3g = ground_flux_robin_bc(T3, Tg_deep, h_ground) # Robin BC ground flux
    dT3  = (q_32 + q_3g) / Cs_layer

    return np.array([dT1, dT2, dT3])

@njit
def compute_h_out(wind_speed):
    """
    Calculate external heat transfer coefficient based on wind speed.
    
    Args:
        wind_speed: Wind speed [m/s]
    
    Returns:
        h_out: External heat transfer coefficient [W/m²K]
    """
    if wind_speed <= 5.0:
        return 6.0 + 4.0 * wind_speed
    else:
        return 7.41 * (wind_speed**0.78)

# ============================================================
#  Insulation step (advanced model) - UPDATED
# ============================================================
def insulation_step(state_in, forc, p, dt):
    """Advanced insulation model with moisture, age, rain, and evaporation."""
    W        = state_in["W"]
    age_days = state_in["age_days"]

    f      = np.clip(W / p["W_sat"], 0.0, 1.0)
    age_yr = age_days / 365.0

    k_moist     = p["k_dry"] + (p["k_sat"] - p["k_dry"]) * (f**p["n_k"])
    k_age_factor = 1.0 + p["delta_k_age"] * (1.0 - np.exp(-age_yr / p["tau_k_years"]))
    k_eff       = k_moist * k_age_factor

    alpha_moist = p["alpha_dry"] + (p["alpha_wet"] - p["alpha_dry"]) * (f**p["n_alpha"])
    alpha_age   = alpha_moist + p["delta_alpha_age"] * (1.0 - np.exp(-age_yr / p["tau_alpha_years"]))
    alpha_eff   = np.clip(alpha_age, 0.0, 1.0)

    q_solar = alpha_eff * forc["Isolar"]

    eta_rain   = max(0.0, 1.0 - f)
    P_in_mass  = eta_rain * p["rho_w"] * forc["Prain"]

    zeta_rain   = p["zeta0"] * np.exp(-p["gamma_H"] * p["Hi"]) * np.exp(-p["gamma_W"] * f)
    q_rain_snow = zeta_rain * p["rho_w"] * p["c_w"] * forc["Prain"] * (forc["T_rain"] - p["Tfreeze"])

    Tc_s = p["Tfreeze"] - 273.15
    Tc_a = forc["Ta"] - 273.15

    e_sat_s = 611.0 * np.exp(17.27 * Tc_s / (Tc_s + 237.3))
    e_sat_a = 611.0 * np.exp(17.27 * Tc_a / (Tc_a + 237.3))

    e_surf = e_sat_s
    e_air  = forc["RH"] * e_sat_a

    VPD = max(0.0, e_surf - e_air)

    E0 = p["rho_air"] * p["C_E"] * p["U10"] * VPD / p["P0"]
    f_breath = np.exp(-p["beta_w"] * f)
    E = E0 * f_breath

    q_evap = -p["Lv"] * E

    D     = p["K_D"] * max(0.0, W - p["W_field"])
    dWdt  = P_in_mass - E - D
    W_new = W + dt * dWdt
    W_new = np.clip(W_new, 0.0, p["W_sat"])

    age_days_new = age_days + dt / 86400.0

    state_out = dict(state_in)
    state_out["W"]        = W_new
    state_out["age_days"] = age_days_new
    state_out["k_eff"]    = k_eff
    state_out["alpha_eff"] = alpha_eff
    state_out["f_sat"]    = f
    
    # Compute multi-layer resistance if enabled
    T_surf_approx = forc["Ta"]  # Approximation for surface temperature
    T_snow_approx = T[0]  # Top snow layer temperature
    h_out = forc["h_out"]  # Get h_out from forcing
    h_in = 99.75  # Internal heat transfer coefficient [W/m²K]
    
    R_ins, T_profile = compute_insulation_resistance_multilayer(
        T_surf_approx, 
        h_out, 
        T_snow_approx, 
        h_in, 
        k_eff,
        dt_main=dt
    )
    
    return R_ins, q_solar, q_rain_snow, q_evap, state_out

# ============================================================
#  Main simulation
# ============================================================
def main():
    global T, LWC, ice_fractions, InsState, T_ins_profile

    # Reset insulation profile for clean start
    T_ins_profile = None
    
    print("="*60)
    print("Enhanced Snow Storage RC Model with Real Data")
    print("="*60)
    
    # Load real meteorological data
    print("\nLoading meteorological data from DATA.csv...")
    try:
        met_data = read_csv_data('DATA.csv')
        print(f"Loaded {len(met_data['temp'])} hourly data points")
        print(f"Period: {met_data['time'][0]} to {met_data['time'][-1]}")
    except FileNotFoundError:
        print("ERROR: DATA.csv not found!")
        print("Please ensure DATA.csv is in the same directory as this script.")
        return
    
    # Time integration settings
    t0 = 0.0
    dt = 600.0  # 10 min (much finer than 1-hour data)
    dt_data = 3600.0  # Input data time step [s]
    
    # Simulate for the length of available data
    n_hours = len(met_data['temp'])
    tf = n_hours * dt_data
    t_vec = np.arange(t0, tf, dt)
    Nt = len(t_vec)
    
    print(f"\nSimulation settings:")
    print(f"  Duration: {n_hours} hours ({n_hours/24:.1f} days)")
    print(f"  Time step: {dt} s ({dt/60:.1f} min)")
    print(f"  Total steps: {Nt}")
    print(f"  Multi-layer insulation: {USE_MULTILAYER_INSULATION} ({N_ins} layers)")
    print(f"  Refreezing: {USE_REFREEZING}")
    print(f"  Percolation: {USE_PERCOLATION}")
    
    # Initialize history arrays
    T_hist = np.zeros((Nt, 3))
    T_hist[0, :] = T
    LWC_hist = np.zeros((Nt, 3))
    LWC_hist[0, :] = LWC
    
    qnet_surf_hist = np.zeros(Nt)
    qa_hist        = np.zeros(Nt)
    qsolar_hist    = np.zeros(Nt)
    qrain_hist     = np.zeros(Nt)
    qevap_hist     = np.zeros(Nt)
    qground_hist   = np.zeros(Nt)
    
    Ta_hist     = np.zeros(Nt)
    Isolar_hist = np.zeros(Nt)
    Prain_hist  = np.zeros(Nt)
    
    melt_rate_hist = np.zeros(Nt)
    refrozen_hist  = np.zeros(Nt)
    runoff_hist    = np.zeros(Nt)
    E_melt = 0.0
    
    W_hist     = np.zeros(Nt)
    k_eff_hist = np.zeros(Nt)
    alpha_hist = np.zeros(Nt)
    Rins_hist  = np.zeros(Nt)
    fsat_hist  = np.zeros(Nt)
    
    # Progress tracking
    progress_points = [int(Nt * p) for p in [0.25, 0.5, 0.75, 1.0]]
    
    # Main time loop
    print("\nRunning simulation...")
    for k in range(Nt-1):
        t = t_vec[k]
        t_mid = t + dt/2.0
        
        # Show progress
        if k in progress_points:
            pct = int(100 * k / Nt)
            print(f"  Progress: {pct}%")
        
        # Interpolate forcing data from hourly measurements
        Ta_C = interpolate_data(met_data['temp'], t_mid, dt_data)
        Isolar = interpolate_data(met_data['solar'], t_mid, dt_data)
        Prain_mh = interpolate_data(met_data['precip'], t_mid, dt_data)
        wind = interpolate_data(met_data['wind'], t_mid, dt_data)
        RH_pct = interpolate_data(met_data['rh'], t_mid, dt_data)
        
        # Convert units
        Ta_K = Ta_C + 273.15
        Prain = Prain_mh / 3600.0  # [m/h] -> [m/s]
        RH_frac = RH_pct / 100.0

        h_out = compute_h_out(wind)

        forc = {
            "Isolar": Isolar,
            "Prain":  Prain,
            "T_rain": Ta_K,
            "RH":     RH_frac,
            "Ta":     Ta_K,
            "U10":    wind,
            "h_out":  h_out 
        }
        
        # Update insulation
        if USE_ADVANCED_INSULATION:
            R_ins, q_solar_ins, q_rain_snow, q_evap, InsState = \
                insulation_step(InsState, forc, InsPar, dt)
            
            W_hist[k]      = InsState["W"]
            k_eff_hist[k]  = InsState["k_eff"]
            alpha_hist[k]  = InsState["alpha_eff"]
            Rins_hist[k]   = R_ins
            fsat_hist[k]   = InsState["f_sat"]
        else:
            R_ins       = Hi / k_i_base
            q_solar_ins = alpha_const * forc["Isolar"]
            q_rain_snow = eta_rain_const * rho_w * c_w * \
                          forc["Prain"] * (forc["T_rain"] - Tfreeze)
            q_evap      = 0.0
        
        R_a2s = R_eff + R_ins
        
        stepPar = {
            "R_a2s":   R_a2s,
            "q_solar": q_solar_ins,
            "q_rain":  q_rain_snow,
            "q_evap":  q_evap,
            "Ta":      Ta_K
        }
        
        # RK4 integration
        k1_vec = dTdt(t, T, stepPar, met_data, dt_data)
        k2_vec = dTdt(t + dt/2.0, T + dt*k1_vec/2.0, stepPar, met_data, dt_data)
        k3_vec = dTdt(t + dt/2.0, T + dt*k2_vec/2.0, stepPar, met_data, dt_data)
        k4_vec = dTdt(t + dt, T + dt*k3_vec, stepPar, met_data, dt_data)
        T_new = T + (dt/6.0) * (k1_vec + 2*k2_vec + 2*k3_vec + k4_vec)
        
        # Calculate fluxes
        T1_mid   = T[0]
        q_a_mid  = (Ta_K - T1_mid) / R_a2s
        q_surf_mid = q_a_mid + q_solar_ins + q_rain_snow + q_evap
        #q_ground_mid = ground_flux(T[2])    # Dirichlet BC ground flux
        q_ground_mid = ground_flux_robin_bc(T[2], Tg_deep, h_ground) # Robin BC ground flux

        
        # Refreezing
        total_refrozen = 0.0
        if USE_REFREEZING:
            for i in range(3):
                T_new[i], LWC[i], ice_fractions[i], refrozen = \
                    refreezing_layer(T_new[i], LWC[i], ice_fractions[i])
                total_refrozen += refrozen
        
        refrozen_hist[k] = total_refrozen
        
        # Surface melting
        surface_melt_water = 0.0
        if T_new[0] > Tfreeze:
            T_new[0] = Tfreeze
            if q_surf_mid > 0.0:
                dE_melt = q_surf_mid * dt
                dM_melt = dE_melt / (rho_i * Lf)
                surface_melt_water = dM_melt
            else:
                dE_melt = 0.0
                dM_melt = 0.0
            E_melt += dE_melt
            melt_rate_hist[k] = dM_melt / dt
            LWC[0] += surface_melt_water / dz_s
        else:
            melt_rate_hist[k] = 0.0
        
        # Percolation
        runoff = 0.0
        if USE_PERCOLATION:
            LWC, runoff = percolate_water(LWC, heights, theta_e)
        
        runoff_hist[k] = runoff
        
        # Store results
        T = T_new
        T_hist[k+1, :] = T
        LWC_hist[k+1, :] = LWC
        
        qnet_surf_hist[k] = q_surf_mid
        qa_hist[k]        = q_a_mid
        qsolar_hist[k]    = q_solar_ins
        qrain_hist[k]     = q_rain_snow
        qevap_hist[k]     = q_evap
        qground_hist[k]   = q_ground_mid
        
        Ta_hist[k]     = Ta_K
        Isolar_hist[k] = Isolar
        Prain_hist[k]  = Prain
    
    # Final values
    Ta_hist[-1] = interpolate_data(met_data['temp'], t_vec[-1], dt_data) + 273.15
    Isolar_hist[-1] = interpolate_data(met_data['solar'], t_vec[-1], dt_data)
    Prain_hist[-1] = interpolate_data(met_data['precip'], t_vec[-1], dt_data) / 3600.0
    
    if USE_ADVANCED_INSULATION:
        W_hist[-1]      = InsState["W"]
        k_eff_hist[-1]  = InsState["k_eff"]
        alpha_hist[-1]  = InsState["alpha_eff"]
        Rins_hist[-1]   = R_ins
        fsat_hist[-1]   = InsState["f_sat"]
    
    print("  Progress: 100%")
    print("\nSimulation complete!")
    
    # ============================================================
    #  Energy diagnostics
    # ============================================================
    print("\n" + "="*60)
    print("Energy Balance Diagnostics")
    print("="*60)
    
    E_a     = np.trapezoid(qa_hist,      t_vec)
    E_solar = np.trapezoid(qsolar_hist,  t_vec)
    E_rain  = np.trapezoid(qrain_hist,   t_vec)
    E_evap  = np.trapezoid(qevap_hist,   t_vec)
    E_g     = np.trapezoid(qground_hist, t_vec)
    E_total_in = E_a + E_solar + E_rain + E_evap + E_g
    
    E_snow_change = Cs_layer * np.sum(T_hist[-1,:] - T_hist[0,:])
    E_refrozen = np.sum(refrozen_hist) * Lf
    E_balance = E_total_in - (E_snow_change + E_melt + E_refrozen)
    M_melt = E_melt / (rho_i * Lf)
    M_runoff = np.sum(runoff_hist)
    
    print(f"\nConfiguration:")
    print(f"  Advanced insulation:    {USE_ADVANCED_INSULATION}")
    print(f"  Multi-layer insulation: {USE_MULTILAYER_INSULATION} ({N_ins} layers)")
    print(f"  Refreezing:             {USE_REFREEZING}")
    print(f"  Percolation:            {USE_PERCOLATION}")
    
    print(f"\nEnergy fluxes [MJ/m²]:")
    print(f"  Air convection:     {E_a/1e6:>10.3f}")
    print(f"  Solar radiation:    {E_solar/1e6:>10.3f}")
    print(f"  Rain heat:          {E_rain/1e6:>10.3f}")
    print(f"  Evaporation:        {E_evap/1e6:>10.3f}")
    print(f"  Ground heat:        {E_g/1e6:>10.3f}")
    print(f"  ─────────────────────────────")
    print(f"  Total input:        {E_total_in/1e6:>10.3f}")
    
    print(f"\nEnergy storage/losses [MJ/m²]:")
    print(f"  Snow temperature:   {E_snow_change/1e6:>10.3f}")
    print(f"  Melting:            {E_melt/1e6:>10.3f}")
    print(f"  Refreezing:         {E_refrozen/1e6:>10.3f}")
    print(f"  ─────────────────────────────")
    print(f"  Energy residual:    {E_balance/1e6:>10.3f}")
    
    print(f"\nMass balance:")
    print(f"  Total melt:         {M_melt:>10.3f} m w.e.")
    print(f"  Total runoff:       {M_runoff:>10.3f} kg/m²")
    print(f"  Melt rate (avg):    {M_melt/(n_hours/24)*1000:>10.3f} mm/day")
    
    # ============================================================
    #  Plots
    # ============================================================
    print("\nGenerating plots...")
    days = t_vec / (24.0 * 3600.0)
    
    fig = plt.figure(figsize=(14, 10))
    
    # Temperature evolution
    ax1 = plt.subplot(4, 2, 1)
    ax1.plot(days, Ta_hist - 273.15, 'k--', linewidth=1, label='Air temp', alpha=0.7)
    ax1.plot(days, T_hist[:,0] - 273.15, label='T1 (surface)', linewidth=1.5)
    ax1.plot(days, T_hist[:,1] - 273.15, label='T2 (middle)', linewidth=1.5)
    ax1.plot(days, T_hist[:,2] - 273.15, label='T3 (bottom)', linewidth=1.5)
    ax1.axhline(y=0, color='gray', linestyle=':', alpha=0.5)
    ax1.set_ylabel('Temperature [°C]')
    ax1.legend(loc='best', fontsize=8)
    ax1.grid(True, alpha=0.3)
    ax1.set_title('Snow Layer Temperatures')
    
    # Solar radiation
    ax2 = plt.subplot(4, 2, 2)
    ax2.fill_between(days, 0, Isolar_hist, alpha=0.5)
    ax2.set_ylabel('Solar [W/m²]')
    ax2.grid(True, alpha=0.3)
    ax2.set_title('Solar Radiation Input')
    
    # Liquid water content
    ax3 = plt.subplot(4, 2, 3)
    ax3.plot(days, LWC_hist[:,0], label='Layer 1 (surface)', linewidth=1.5)
    ax3.plot(days, LWC_hist[:,1], label='Layer 2 (middle)', linewidth=1.5)
    ax3.plot(days, LWC_hist[:,2], label='Layer 3 (bottom)', linewidth=1.5)
    ax3.axhline(y=theta_e, color='r', linestyle='--', label='Field capacity', linewidth=1)
    ax3.set_ylabel('LWC [-]')
    ax3.legend(loc='best', fontsize=8)
    ax3.grid(True, alpha=0.3)
    ax3.set_title('Liquid Water Content')
    
    # Precipitation
    ax4 = plt.subplot(4, 2, 4)
    ax4.fill_between(days, 0, Prain_hist * 1000.0 * 3600.0, alpha=0.5, color='blue')
    ax4.set_ylabel('Precip [mm/h]')
    ax4.grid(True, alpha=0.3)
    ax4.set_title('Precipitation')
    
    # Melt and runoff
    ax5 = plt.subplot(4, 2, 5)
    cumulative_melt = np.cumsum(melt_rate_hist * dt)
    cumulative_runoff = np.cumsum(runoff_hist)
    ax5.plot(days, cumulative_melt, label='Cumulative melt', linewidth=2)
    ax5.plot(days, cumulative_runoff / rho_w, label='Cumulative runoff', linewidth=2)
    ax5.set_ylabel('Water [m w.e.]')
    ax5.set_xlabel('Time [days]')
    ax5.legend(loc='best', fontsize=8)
    ax5.grid(True, alpha=0.3)
    ax5.set_title('Melt and Runoff')
    
    # Heat fluxes
    ax6 = plt.subplot(4, 2, 6)
    ax6.plot(days, qsolar_hist, label='Solar', linewidth=1, alpha=0.7)
    ax6.plot(days, qrain_hist, label='Rain', linewidth=1, alpha=0.7)
    ax6.plot(days, qevap_hist, label='Evaporation', linewidth=1, alpha=0.7)
    ax6.plot(days, qa_hist, label='Air convection', linewidth=1, alpha=0.7)
    ax6.set_ylabel('Heat flux [W/m²]')
    ax6.set_xlabel('Time [days]')
    ax6.legend(loc='best', fontsize=8)
    ax6.grid(True, alpha=0.3)
    ax6.set_title('Surface Heat Fluxes')
    
    # Insulation properties (if advanced model)
    if USE_ADVANCED_INSULATION:
        ax7 = plt.subplot(4, 2, 7)
        ax7_twin = ax7.twinx()
        ln1 = ax7.plot(days, k_eff_hist, 'b-', label='k_eff', linewidth=1.5)
        ln2 = ax7_twin.plot(days, alpha_hist, 'r-', label='α_eff', linewidth=1.5)
        ax7.set_ylabel('k_eff [W/mK]', color='b')
        ax7_twin.set_ylabel('α_eff [-]', color='r')
        ax7.tick_params(axis='y', labelcolor='b')
        ax7_twin.tick_params(axis='y', labelcolor='r')
        lns = ln1 + ln2
        labs = [l.get_label() for l in lns]
        ax7.legend(lns, labs, loc='best', fontsize=8)
        ax7.grid(True, alpha=0.3)
        ax7.set_xlabel('Time [days]')
        ax7.set_title('Insulation Properties')
        
        ax8 = plt.subplot(4, 2, 8)
        ax8.plot(days, W_hist, linewidth=1.5, color='steelblue')
        ax8.set_ylabel('Moisture [kg/m²]')
        ax8.set_xlabel('Time [days]')
        ax8.grid(True, alpha=0.3)
        ax8.set_title('Insulation Moisture Content')
    
    plt.tight_layout()
    plt.savefig('snow_storage_simulation_v2.png', dpi=150, bbox_inches='tight')
    print("  Saved plot: snow_storage_simulation_v2.png")
    plt.show()
    
    print("\n" + "="*60)
    print("Analysis complete!")
    print("="*60)

if __name__ == "__main__":
    main()
