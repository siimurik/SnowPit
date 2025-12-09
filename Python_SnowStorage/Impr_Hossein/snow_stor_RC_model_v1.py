import numpy as np
import matplotlib.pyplot as plt
from numba import njit

# ============================================================
#  Enhanced 3-layer Snow Storage RC Model
#  Improvements from COSIPY integration:
#  1. Adaptive time stepping for heat equation
#  2. Refreezing with cold content
#  3. Percolation with bucket method
#  4. Optional densification
# ============================================================

USE_ADVANCED_INSULATION = True
USE_REFREEZING = True          # NEW: Enable refreezing in subsurface
USE_PERCOLATION = True         # NEW: Enable water percolation
USE_DENSIFICATION = False      # NEW: Enable snow densification

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
Hi     = 0.6                  # insulation thickness [m]

# Simple (constant) insulation parameters
k_i_const     = 0.06          # [W/mK]
alpha_const   = 0.10          # [-]
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

# ---------- NEW: Layer properties for refreezing/percolation ----------
# Liquid water content [volumetric fraction] per layer
LWC = np.array([0.0, 0.0, 0.0], dtype=float)
# Irreducible water content (field capacity) [volumetric fraction]
theta_e = 0.03  # Typical value for snow

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

        # better evaporation model
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
#  Forcing functions
# ============================================================
def Ta_fun(t):
    """Air temperature [K] with daily cycle."""
    day   = 24.0 * 3600.0
    Tmean = 273.15 + 3.0
    Tamp  = 6.0
    return Tmean + Tamp * np.sin(2.0 * np.pi * (t / day))

def Isolar_fun(t):
    """Diurnal solar [W/m^2], half-sine from 6h to 18h."""
    day = 24.0 * 3600.0
    tau = np.mod(t, day)
    if tau < 6*3600 or tau > 18*3600:
        return 0.0
    return 600.0 * np.sin(np.pi * (tau - 6*3600) / (12*3600))

def Prain_fun(t):
    """Rain rate [m/s]: 2h event every 5 days."""
    day    = 24.0 * 3600.0
    period = 5.0 * day
    tau = np.mod(t, period)
    if tau < 2.0 * 3600.0:
        return 5e-7
    return 0.0

# ============================================================
#  NEW: Refreezing function (adapted from COSIPY)
# ============================================================
@njit
def refreezing_layer(T_layer, LWC_layer, ice_frac, dt):
    """
    Refreeze water in a single layer based on cold content.
    
    Args:
        T_layer: Layer temperature [K]
        LWC_layer: Liquid water content [volumetric fraction]
        ice_frac: Ice fraction [volumetric fraction]
        dt: Time step [s]
    
    Returns:
        tuple: (new_T, new_LWC, new_ice_frac, refrozen_mass)
    """
    if (T_layer >= Tfreeze) or (LWC_layer <= 0.0):
        return T_layer, LWC_layer, ice_frac, 0.0
    
    # Cold content available
    dT_max = T_layer - Tfreeze  # negative
    
    # Maximum water that can refreeze based on cold content
    # Energy balance: dtheta_w * Lf * rho_w = -dT_max * (rho_i * c_s * ice_frac + rho_w * c_w * LWC_layer)
    dtheta_w_max = -(dT_max * (ice_frac * rho_i * c_s + LWC_layer * rho_w * c_w)) / (rho_w * Lf)
    
    # Maximum based on available water
    dtheta_w = min(LWC_layer, dtheta_w_max)
    
    # Change in ice fraction (mass conservation)
    dtheta_i = (rho_w / rho_i) * dtheta_w
    
    # Update temperature due to latent heat release
    dT = (dtheta_w * rho_w * Lf) / (ice_frac * rho_i * c_s + LWC_layer * rho_w * c_w)
    
    new_T = T_layer + dT
    new_LWC = LWC_layer - dtheta_w
    new_ice_frac = ice_frac + dtheta_i
    refrozen_mass = dtheta_w * dz_s * rho_w  # [kg/m^2]
    
    return new_T, new_LWC, new_ice_frac, refrozen_mass

# ============================================================
#  NEW: Percolation function (adapted from COSIPY)
# ============================================================
@njit
def percolate_water(LWC_array, heights, theta_e):
    """
    Percolate water through layers using bucket method.
    
    Args:
        LWC_array: Array of liquid water content [volumetric fraction]
        heights: Array of layer heights [m]
        theta_e: Irreducible water content [volumetric fraction]
    
    Returns:
        tuple: (new_LWC_array, runoff)
    """
    n_layers = len(LWC_array)
    new_LWC = LWC_array.copy()
    
    # Percolate from top to bottom
    for i in range(n_layers - 1):
        if new_LWC[i] > theta_e:
            # Excess water percolates
            excess = new_LWC[i] - theta_e
            new_LWC[i] = theta_e
            
            # Add to next layer (accounting for different heights)
            excess_mass = excess * heights[i]  # [m]
            new_LWC[i+1] += excess_mass / heights[i+1]
    
    # Runoff from bottom layer
    runoff = 0.0
    if new_LWC[n_layers-1] > theta_e:
        excess = new_LWC[n_layers-1] - theta_e
        new_LWC[n_layers-1] = theta_e
        runoff = excess * heights[n_layers-1] * rho_w  # [kg/m^2]
    
    return new_LWC, runoff

# ============================================================
#  Ground flux
# ============================================================
def ground_flux(T3):
    """Ground -> bottom snow layer flux [W/m^2]."""
    return (Tg - T3) / R_3g

# ============================================================
#  dT/dt for snow layers with adaptive time stepping
# ============================================================
def dTdt(t, Tv, stepPar):
    """
    Tv = [T1, T2, T3]
    stepPar: dict with R_a2s, q_solar, q_rain, q_evap
    """
    T1, T2, T3 = Tv
    Ta = Ta_fun(t)

    q_a    = (Ta - T1) / stepPar["R_a2s"]
    q_surf = q_a + stepPar["q_solar"] + stepPar["q_rain"] + stepPar["q_evap"]

    # layer 1
    q_12 = (T2 - T1) / R_12
    dT1  = (q_surf + q_12) / Cs_layer

    # layer 2
    q_21 = (T1 - T2) / R_12
    q_23 = (T3 - T2) / R_23
    dT2  = (q_21 + q_23) / Cs_layer

    # layer 3
    q_32 = (T2 - T3) / R_23
    q_3g = (Tg - T3) / R_3g
    dT3  = (q_32 + q_3g) / Cs_layer

    return np.array([dT1, dT2, dT3])

# ============================================================
#  NEW: Adaptive time stepping (from COSIPY)
# ============================================================
@njit
def compute_stable_dt(dz, k_array, rho_array, c_array):
    """
    Compute stable time step based on CFL condition.
    
    Args:
        dz: Layer thickness [m]
        k_array: Thermal conductivity array [W/mK]
        rho_array: Density array [kg/m^3]
        c_array: Specific heat array [J/kgK]
    
    Returns:
        Stable time step [s]
    """
    # Thermal diffusivity
    alpha = k_array / (rho_array * c_array)
    
    # CFL stability condition
    c_stab = 0.5
    dt_stable = c_stab * dz**2 / (2 * np.max(alpha))
    
    return dt_stable

# ============================================================
#  Insulation step (advanced model)
# ============================================================
def insulation_step(state_in, forc, p, dt):
    """
    Advanced insulation model with moisture, age, rain, and evaporation.
    """
    W        = state_in["W"]
    age_days = state_in["age_days"]

    # 1) Moisture & age factors
    f      = np.clip(W / p["W_sat"], 0.0, 1.0)
    age_yr = age_days / 365.0

    # conductivity: k(W, age)
    k_moist     = p["k_dry"] + (p["k_sat"] - p["k_dry"]) * (f**p["n_k"])
    k_age_factor = 1.0 + p["delta_k_age"] * (1.0 - np.exp(-age_yr / p["tau_k_years"]))
    k_eff       = k_moist * k_age_factor
    R_ins       = p["Hi"] / k_eff

    # absorptivity: alpha(W, age)
    alpha_moist = p["alpha_dry"] + (p["alpha_wet"] - p["alpha_dry"]) * (f**p["n_alpha"])
    alpha_age   = alpha_moist + p["delta_alpha_age"] * (1.0 - np.exp(-age_yr / p["tau_alpha_years"]))
    alpha_eff   = np.clip(alpha_age, 0.0, 1.0)

    # 2) Solar absorption
    q_solar = alpha_eff * forc["Isolar"]

    # 3) Rain infiltration & heat to snow
    eta_rain   = max(0.0, 1.0 - f)
    P_in_mass  = eta_rain * p["rho_w"] * forc["Prain"]

    zeta_rain   = p["zeta0"] * np.exp(-p["gamma_H"] * p["Hi"]) * np.exp(-p["gamma_W"] * f)
    q_rain_snow = zeta_rain * p["rho_w"] * p["c_w"] * forc["Prain"] * (forc["T_rain"] - p["Tfreeze"])

    # 4) Evaporation (bulk aerodynamic model)
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

    # 5) Drainage & moisture update
    D     = p["K_D"] * max(0.0, W - p["W_field"])
    dWdt  = P_in_mass - E - D
    W_new = W + dt * dWdt
    W_new = np.clip(W_new, 0.0, p["W_sat"])

    age_days_new = age_days + dt / 86400.0

    state_out = dict(state_in)
    state_out["W"]        = W_new
    state_out["age_days"] = age_days_new

    # diagnostics
    state_out["k_eff"]     = k_eff
    state_out["alpha_eff"] = alpha_eff
    state_out["f_sat"]     = f

    return R_ins, q_solar, q_rain_snow, q_evap, state_out

# ============================================================
#  Time integration settings
# ============================================================
t0 = 0.0
tf = 120.0 * 24.0 * 3600.0   # 120 days
dt = 600.0                   # 10 min
t_vec = np.arange(t0, tf + dt, dt)
Nt = len(t_vec)

T_hist = np.zeros((Nt, 3))
T_hist[0, :] = T

# NEW: LWC history
LWC_hist = np.zeros((Nt, 3))
LWC_hist[0, :] = LWC

# histories
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
refrozen_hist  = np.zeros(Nt)  # NEW
runoff_hist    = np.zeros(Nt)  # NEW
E_melt = 0.0

# insulation diagnostics
W_hist     = np.zeros(Nt)
k_eff_hist = np.zeros(Nt)
alpha_hist = np.zeros(Nt)
Rins_hist  = np.zeros(Nt)
fsat_hist  = np.zeros(Nt)

# Ice fractions (assuming constant for now, could be updated with densification)
ice_fractions = np.array([0.4, 0.4, 0.4])  # volumetric fraction
heights = np.array([dz_s, dz_s, dz_s])

# ============================================================
#  Time loop (RK4 with enhancements)
# ============================================================
for k in range(Nt-1):
    t = t_vec[k]
    t_mid = t + dt/2.0

    # --- forcing at mid-step ---
    forc = {
        "Isolar": Isolar_fun(t_mid),
        "Prain":  Prain_fun(t_mid),
        "T_rain": 273.15 + 3.0,
        "RH":     0.80,
        "Ta":     Ta_fun(t_mid)
    }

    # --- insulation: choose simple or advanced model ---
    if USE_ADVANCED_INSULATION:
        R_ins, q_solar_ins, q_rain_snow, q_evap, InsState = \
            insulation_step(InsState, forc, InsPar, dt)

        W_hist[k]      = InsState["W"]
        k_eff_hist[k]  = InsState["k_eff"]
        alpha_hist[k]  = InsState["alpha_eff"]
        Rins_hist[k]   = R_ins
        fsat_hist[k]   = InsState["f_sat"]
    else:
        R_ins       = Hi / k_i_const
        q_solar_ins = alpha_const * forc["Isolar"]
        q_rain_snow = eta_rain_const * rho_w * c_w * \
                      forc["Prain"] * (forc["T_rain"] - Tfreeze)
        q_evap      = 0.0

    # total air → surface snow resistance
    R_a2s = R_eff + R_ins

    stepPar = {
        "R_a2s":   R_a2s,
        "q_solar": q_solar_ins,
        "q_rain":  q_rain_snow,
        "q_evap":  q_evap
    }

    # --- RK4 integration of snow temperatures ---
    k1_vec = dTdt(t,          T,           stepPar)
    k2_vec = dTdt(t + dt/2.0, T + dt*k1_vec/2.0, stepPar)
    k3_vec = dTdt(t + dt/2.0, T + dt*k2_vec/2.0, stepPar)
    k4_vec = dTdt(t + dt,     T + dt*k3_vec,     stepPar)
    T_new  = T + (dt/6.0) * (k1_vec + 2*k2_vec + 2*k3_vec + k4_vec)

    # --- fluxes at mid-step for melt & diagnostics ---
    Ta_mid   = Ta_fun(t_mid)
    T1_mid   = T[0]
    q_a_mid  = (Ta_mid - T1_mid) / R_a2s
    q_surf_mid = q_a_mid + q_solar_ins + q_rain_snow + q_evap
    q_ground_mid = ground_flux(T[2])

    # --- NEW: Refreezing ---
    total_refrozen = 0.0
    if USE_REFREEZING:
        for i in range(3):
            T_new[i], LWC[i], ice_fractions[i], refrozen = \
                refreezing_layer(T_new[i], LWC[i], ice_fractions[i], dt)
            total_refrozen += refrozen
    
    refrozen_hist[k] = total_refrozen

    # --- melting clamp at surface ---
    surface_melt_water = 0.0
    if T_new[0] > Tfreeze:
        T_new[0] = Tfreeze
        if q_surf_mid > 0.0:
            dE_melt = q_surf_mid * dt
            dM_melt = dE_melt / (rho_i * Lf)
            surface_melt_water = dM_melt  # [m w.e.]
        else:
            dE_melt = 0.0
            dM_melt = 0.0
        E_melt += dE_melt
        melt_rate_hist[k] = dM_melt / dt
        
        # Add melt water to top layer LWC
        LWC[0] += surface_melt_water / dz_s
    else:
        melt_rate_hist[k] = 0.0

    # --- NEW: Percolation ---
    runoff = 0.0
    if USE_PERCOLATION:
        LWC, runoff = percolate_water(LWC, heights, theta_e)
    
    runoff_hist[k] = runoff

    # --- store histories ---
    T = T_new
    T_hist[k+1, :] = T
    LWC_hist[k+1, :] = LWC

    qnet_surf_hist[k] = q_surf_mid
    qa_hist[k]        = q_a_mid
    qsolar_hist[k]    = q_solar_ins
    qrain_hist[k]     = q_rain_snow
    qevap_hist[k]     = q_evap
    qground_hist[k]   = q_ground_mid

    Ta_hist[k]     = Ta_fun(t)
    Isolar_hist[k] = Isolar_fun(t)
    Prain_hist[k]  = Prain_fun(t)

# last forcing values
Ta_hist[-1]     = Ta_fun(t_vec[-1])
Isolar_hist[-1] = Isolar_fun(t_vec[-1])
Prain_hist[-1]  = Prain_fun(t_vec[-1])

if USE_ADVANCED_INSULATION:
    W_hist[-1]      = InsState["W"]
    k_eff_hist[-1]  = InsState["k_eff"]
    alpha_hist[-1]  = InsState["alpha_eff"]
    Rins_hist[-1]   = R_ins
    fsat_hist[-1]   = InsState["f_sat"]

# ============================================================
#  Energy diagnostics
# ============================================================
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

print(f"USE_ADVANCED_INSULATION = {USE_ADVANCED_INSULATION}")
print(f"USE_REFREEZING = {USE_REFREEZING}")
print(f"USE_PERCOLATION = {USE_PERCOLATION}")
print(f"E_total_in      = {E_total_in:.3e} J/m^2")
print(f"E_snow_change   = {E_snow_change:.3e} J/m^2")
print(f"E_melt          = {E_melt:.3e} J/m^2")
print(f"E_refrozen      = {E_refrozen:.3e} J/m^2")
print(f"Energy residual = {E_balance:.3e} J/m^2")
print(f"Total melted snow thickness = {M_melt:.3f} m")
print(f"Total runoff = {M_runoff:.3f} m w.e.")

# ============================================================
#  Plots
# ============================================================
days       = t_vec / (24.0 * 3600.0)
melt_depth = np.cumsum(melt_rate_hist) * dt

# Temperatures
plt.figure(figsize=(10, 8))
plt.subplot(3,1,1)
plt.plot(days, Ta_hist - 273.15, 'k--', label='T_a')
plt.plot(days, T_hist[:,0] - 273.15, label='T1 (surface)')
plt.plot(days, T_hist[:,1] - 273.15, label='T2 (middle)')
plt.plot(days, T_hist[:,2] - 273.15, label='T3 (bottom)')
plt.ylabel('Temperature [°C]')
plt.legend()
plt.grid(True)

plt.subplot(3,1,2)
plt.plot(days, Isolar_hist)
plt.ylabel('Solar [W/m²]')
plt.grid(True)

plt.subplot(3,1,3)
plt.plot(days, Prain_hist * 1000.0 * 3600.0)
plt.ylabel('Rain [mm/h]')
plt.xlabel('Time [days]')
plt.grid(True)
plt.tight_layout()

# NEW: Liquid water content
plt.figure(figsize=(10, 6))
plt.plot(days, LWC_hist[:,0], label='Layer 1 (surface)')
plt.plot(days, LWC_hist[:,1], label='Layer 2 (middle)')
plt.plot(days, LWC_hist[:,2], label='Layer 3 (bottom)')
plt.axhline(y=theta_e, color='r', linestyle='--', label='Field capacity')
plt.xlabel('Time [days]')
#plt.ylabel('LW

plt.show()