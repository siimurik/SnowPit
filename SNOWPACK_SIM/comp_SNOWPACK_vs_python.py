"""
compare_snowpack_vs_python_v4.py
---------------------------------
Runs the snowsim_v4 RC model (3-layer snow, advanced insulation,
refreezing, percolation) and plots it against SNOWPACK .met output.

Flags used:
    USE_ADVANCED_INSULATION = True
    USE_REFREEZING          = True
    USE_PERCOLATION         = True

Usage:
    python compare_snowpack_vs_python_v4.py \\
        --csv   DATA_2024.csv \\
        --met   output/snow_storage_snow_storage.met \\
        --out   comparison_v4.png

Requires: numpy, pandas, matplotlib, numba
"""

import argparse
import csv
from datetime import datetime, timedelta

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from numba import njit


# ============================================================
#  Physical constants  (module-level so @njit functions can see them)
# ============================================================
sigma   = 5.670374419e-8   # Stefan-Boltzmann  [W/m² K⁴]
Lf      = 3.34e5           # Latent heat of fusion  [J/kg]
rho_i   = 917.0            # Ice density  [kg/m³]
rho_s   = 400.0            # Snow bulk density  [kg/m³]
c_s     = 2100.0           # Snow specific heat  [J/(kg K)]
rho_w   = 1000.0           # Water density  [kg/m³]
c_w     = 4180.0           # Water specific heat  [J/(kg K)]
Tfreeze = 273.15           # 0 °C  [K]

# Snow geometry
Hs   = 4.5
Ns   = 3
dz_s = Hs / Ns

# Snow conductivity & layer resistances
k_snow  = 0.423558
R_layer = dz_s / k_snow
R_12    = R_layer
R_23    = R_layer

# Layer heat capacity
Cs_layer = rho_s * c_s * dz_s   # [J/(m² K)]

# Insulation geometry
Hi      = 0.20        # thickness [m]
k_i_base = 0.07       # base conductivity [W/(m K)]
rho_dry  = 200.0      # dry density [kg/m³]
c_dry    = 1.5e3     # dry specific heat [J/(kg K)]

# Advanced insulation parameters (dict consumed by insulation_step)
_INS_PAR = {
    "Hi":             Hi,
    "k_dry":          0.05,
    "k_sat":          0.12,
    "n_k":            1.5,
    "W_sat":          30.0,
    "W_field":        10.0,
    "alpha_dry":      0.05,
    "alpha_wet":      0.08,
    "n_alpha":        1.0,
    "delta_k_age":    0.5,
    "tau_k_years":    2.0,
    "delta_alpha_age":0.05,
    "tau_alpha_years":2.0,
    "zeta0":          0.3,
    "gamma_H":        0.5,
    "gamma_W":        2.0,
    "beta_w":         3.0,
    "K_E":            1e-5,
    "K_D":            5e-6,
    "Lv":             2.5e6,
    "rho_w":          rho_w,
    "c_w":            c_w,
    "Tfreeze":        Tfreeze,
    "rho_air":        1.2,
    "C_E":            1.3e-3,
    "U10":            2.0,
    "P0":             101325.0,
}

# Moisture content used for Thermal diffusivity of insulation
moist_cont = 36.318144   # [%]
c_wet_ins  = (1.0 - moist_cont / 100.0) * c_dry + moist_cont / 100.0 * c_w
rho_wet_ins = rho_dry + moist_cont / 100.0 * rho_w
D_ins       = k_i_base / (c_wet_ins * rho_wet_ins)

# Surface heat transfer
h_conv  = 8.0
epsilon = 0.95
T_mean  = 273.15 + 3.0
h_rad   = 4.0 * epsilon * sigma * T_mean**3
h_eff   = h_conv + h_rad
R_eff   = 1.0 / h_eff

# Ground Robin BC
h_ground = 1.511638   # [W/(m² K)]

# Percolation
theta_e = 0.056854

# Ground insulation (kept for reference but not used in ground Robin BC)
Hg_ins = 0.3
kg_ins = 0.04
R_g_ins = Hg_ins / kg_ins
R_3g    = R_layer + R_g_ins


# ============================================================
#  1.  SNOWPACK .met reader  (unchanged from original)
# ============================================================

def read_snowpack_met(met_path):
    """
    Parse a SNOWPACK .met file into a DataFrame.

    Returns
    -------
    df        : pd.DataFrame  indexed by datetime
    col_units : dict  {column_name: unit_string}
    """
    with open(met_path) as f:
        lines = f.readlines()

    header_block = [l for l in lines if l.startswith(',,') or l.startswith('ID,Date,')]
    col_names_line  = header_block[1]
    col_units_line  = header_block[2]

    col_names    = [c.strip() for c in col_names_line.split(',')]
    col_units_raw = [c.strip() for c in col_units_line.split(',')]
    col_units    = {col_names[i]: col_units_raw[i]
                    for i in range(min(len(col_names), len(col_units_raw)))}

    data_rows = [l.strip().split(',') for l in lines if l[:4].isdigit()]

    seen = {}; deduped = []
    for name in col_names[:len(data_rows[0])]:
        if name in seen:
            seen[name] += 1
            deduped.append(f"{name}_{seen[name]}")
        else:
            seen[name] = 0
            deduped.append(name)

    df = pd.DataFrame(data_rows, columns=deduped)
    df['Date'] = pd.to_datetime(df['Date'].str.strip(), format='%d.%m.%Y %H:%M:%S')
    df = df.set_index('Date').sort_index()

    for c in df.columns:
        if c not in ('ID',):
            df[c] = pd.to_numeric(df[c], errors='coerce').replace(-999.0, np.nan)

    return df, col_units


# ============================================================
#  2.  CSV met reader
# ============================================================

def read_met_csv(csv_path):
    data = {'temp': [], 'wind': [], 'precip': [], 'solar': [], 'rh': [], 'soil': []}
    with open(csv_path) as f:
        for r in csv.DictReader(f):
            data['temp'].append(float(r['Temp_C']))
            data['wind'].append(float(r['Air_Vel_m/s_10m']))
            data['precip'].append(float(r['Prec_m/h']))
            data['solar'].append(float(r['Glo_Sol_Ir_W/m2']))
            data['rh'].append(float(r['RH_%']))
            data['soil'].append(float(r['Soil_Temp_320cm']))
    return {k: np.array(v, dtype=np.float64) for k, v in data.items()}


# ============================================================
#  3.  Numba helper: linear interpolation
# ============================================================

@njit
def _interp(vec, t_s, dt_data=3600.0):
    idx   = t_s / dt_data
    lo    = int(np.floor(idx))
    hi    = lo + 1
    n     = len(vec)
    if lo < 0:  return vec[0]
    if hi >= n: return vec[n - 1]
    frac  = idx - lo
    return vec[lo] * (1.0 - frac) + vec[hi] * frac


# ============================================================
#  4.  Numba: dT/dt for 3-layer snow
# ============================================================

@njit
def _dTdt(t, Tv,
          R_a2s, q_solar, q_rain, q_evap,
          temp_arr, soil_arr, dt_data,
          R_12_, R_23_, Cs_, h_gnd, Tfz):
    T1, T2, T3 = Tv[0], Tv[1], Tv[2]

    Ta      = _interp(temp_arr, t, dt_data) + 273.15
    q_a     = (Ta - T1) / R_a2s
    q_surf  = q_a + q_solar + q_rain + q_evap
    q_12    = (T2 - T1) / R_12_
    dT1     = (q_surf + q_12) / Cs_

    q_21    = (T1 - T2) / R_12_
    q_23    = (T3 - T2) / R_23_
    dT2     = (q_21 + q_23) / Cs_

    Tsoil   = _interp(soil_arr, t, dt_data) + 273.15
    R_cond  = 1.0 / 1.5
    R_intf  = 1.0 / h_gnd
    q_3g    = (Tsoil - T3) / (R_cond + R_intf)
    q_32    = (T2 - T3) / R_23_
    dT3     = (q_32 + q_3g) / Cs_

    out = np.empty(3)
    out[0] = dT1; out[1] = dT2; out[2] = dT3
    return out


@njit
def rk4_step(t, T, dt,
             R_a2s, q_solar, q_rain, q_evap,
             temp_arr, soil_arr, dt_data,
             R_12_, R_23_, Cs_, h_gnd, Tfz):
    args = (R_a2s, q_solar, q_rain, q_evap,
            temp_arr, soil_arr, dt_data,
            R_12_, R_23_, Cs_, h_gnd, Tfz)
    k1 = _dTdt(t,          T,               *args)
    k2 = _dTdt(t + dt/2.0, T + dt*k1/2.0,  *args)
    k3 = _dTdt(t + dt/2.0, T + dt*k2/2.0,  *args)
    k4 = _dTdt(t + dt,     T + dt*k3,       *args)
    return T + (dt / 6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4)


# ============================================================
#  5.  Numba: refreezing
# ============================================================

@njit
def refreezing_layer(T_layer, LWC_layer, ice_frac):
    if (T_layer >= Tfreeze) or (LWC_layer <= 0.0):
        return T_layer, LWC_layer, ice_frac, 0.0
    dT_max      = T_layer - Tfreeze
    dtheta_w_max = -(dT_max * (ice_frac * rho_i * c_s + LWC_layer * rho_w * c_w)) / \
                    (rho_w * (Lf - dT_max * (c_s - c_w)))
    dtheta_w    = min(LWC_layer, dtheta_w_max)
    dtheta_i    = (rho_w / rho_i) * dtheta_w
    dT          = (dtheta_w * rho_w * Lf) / \
                   (ice_frac * rho_i * c_s + LWC_layer * rho_w * c_w)
    return (T_layer + dT,
            LWC_layer - dtheta_w,
            ice_frac  + dtheta_i,
            dtheta_w * dz_s * rho_w)


# ============================================================
#  6.  Numba: percolation
# ============================================================

@njit
def percolate_water(LWC_arr, heights, th_e):
    n    = len(LWC_arr)
    nLWC = LWC_arr.copy()
    for i in range(n - 1):
        if nLWC[i] > th_e:
            exc      = nLWC[i] - th_e
            nLWC[i]  = th_e
            nLWC[i+1] += exc * heights[i] / heights[i+1]
    runoff = 0.0
    if nLWC[n-1] > th_e:
        exc       = nLWC[n-1] - th_e
        nLWC[n-1] = th_e
        runoff    = exc * heights[n-1] * rho_w
    return nLWC, runoff


# ============================================================
#  7.  Advanced insulation step  (pure Python, runs once per dt)
# ============================================================

def insulation_step(state, forc, p, dt):
    """
    Update insulation moisture/age and return effective thermal resistance
    plus snow-surface heat fluxes (solar, rain, evap).
    """
    W        = state["W"]
    age_days = state["age_days"]

    f       = np.clip(W / p["W_sat"], 0.0, 1.0)
    age_yr  = age_days / 365.0

    k_moist     = p["k_dry"] + (p["k_sat"] - p["k_dry"]) * (f**p["n_k"])
    k_age_fac   = 1.0 + p["delta_k_age"] * (1.0 - np.exp(-age_yr / p["tau_k_years"]))
    k_eff       = k_moist * k_age_fac

    alpha_moist = p["alpha_dry"] + (p["alpha_wet"] - p["alpha_dry"]) * (f**p["n_alpha"])
    alpha_age   = alpha_moist + p["delta_alpha_age"] * \
                  (1.0 - np.exp(-age_yr / p["tau_alpha_years"]))
    alpha_eff   = np.clip(alpha_age, 0.0, 1.0)

    q_solar     = alpha_eff * forc["Isolar"]

    eta_rain    = max(0.0, 1.0 - f)
    P_in_mass   = eta_rain * p["rho_w"] * forc["Prain"]

    zeta_rain   = p["zeta0"] * np.exp(-p["gamma_H"] * p["Hi"]) * np.exp(-p["gamma_W"] * f)
    q_rain_snow = zeta_rain * p["rho_w"] * p["c_w"] * \
                  forc["Prain"] * (forc["T_rain"] - p["Tfreeze"])

    Tc_s  = p["Tfreeze"] - 273.15
    Tc_a  = forc["Ta"] - 273.15
    e_sat_s = 611.0 * np.exp(17.27 * Tc_s / (Tc_s + 237.3))
    e_sat_a = 611.0 * np.exp(17.27 * Tc_a / (Tc_a + 237.3))
    VPD   = max(0.0, e_sat_s - forc["RH"] * e_sat_a)
    E0    = p["rho_air"] * p["C_E"] * p["U10"] * VPD / p["P0"]
    E     = E0 * np.exp(-p["beta_w"] * f)
    q_evap = -p["Lv"] * E

    D_drain = p["K_D"] * max(0.0, W - p["W_field"])
    W_new   = np.clip(W + dt * (P_in_mass - E - D_drain), 0.0, p["W_sat"])

    R_ins = Hi / k_eff

    state_out = dict(state)
    state_out["W"]        = W_new
    state_out["age_days"] = age_days + dt / 86400.0
    state_out["k_eff"]    = k_eff
    state_out["alpha_eff"]= alpha_eff
    state_out["f_sat"]    = f

    return R_ins, q_solar, q_rain_snow, q_evap, state_out


# ============================================================
#  8.  Numba: external HTC from wind speed
# ============================================================

@njit
def compute_h_out(wind_speed):
    if wind_speed <= 5.0:
        return 6.0 + 4.0 * wind_speed
    return 7.41 * (wind_speed**0.78)


# ============================================================
#  9.  Main simulation  (snowsim_v4 RC model, full flags)
# ============================================================

def run_python_sim(met, dt=600.0):
    """
    3-layer RC snow model with:
        USE_ADVANCED_INSULATION = True
        USE_REFREEZING          = True
        USE_PERCOLATION         = True

    Parameters
    ----------
    met : dict of np.ndarray   hourly forcing arrays from read_met_csv()
    dt  : float                model timestep [s]  (default 600 s)

    Returns
    -------
    dict with keys:
        dates, T_hist, LWC_hist, SWE, melt_cumul, runoff_cumul,
        q_ground, qa, q_solar, q_rain, q_evap,
        k_eff_hist, alpha_hist, W_hist
    """
    n_hours = len(met['temp'])
    t_vec   = np.arange(0.0, n_hours * 3600.0, dt)
    Nt      = len(t_vec)
    dt_data = 3600.0

    # Initial conditions
    T            = np.array([273.15 - 2.0, 273.15 - 4.0, 273.15 - 6.0])
    LWC          = np.zeros(3)
    ice_fractions = np.array([0.4, 0.4, 0.4])
    heights      = np.array([dz_s, dz_s, dz_s])
    InsState     = {"W": 5.0, "age_days": 0.0}

    # Output arrays
    T_hist        = np.zeros((Nt, 3))
    LWC_hist      = np.zeros((Nt, 3))
    melt_rate_h   = np.zeros(Nt)
    runoff_hist   = np.zeros(Nt)
    qa_hist       = np.zeros(Nt)
    qsolar_hist   = np.zeros(Nt)
    qrain_hist    = np.zeros(Nt)
    qevap_hist    = np.zeros(Nt)
    qground_hist  = np.zeros(Nt)
    k_eff_hist    = np.zeros(Nt)
    alpha_hist    = np.zeros(Nt)
    W_hist        = np.zeros(Nt)

    T_hist[0]   = T
    LWC_hist[0] = LWC

    temp_arr = met['temp'].astype(np.float64)
    soil_arr = met['soil'].astype(np.float64)

    E_melt = 0.0

    def _interp_py(vec, t_s):
        idx = t_s / dt_data
        lo  = max(int(np.floor(idx)), 0)
        hi  = min(lo + 1, len(vec) - 1)
        frac = idx - lo
        return vec[lo] * (1 - frac) + vec[hi] * frac

    print("Running snowsim_v4 RC model...")
    progress_pts = {int(Nt * p): f"{int(p*100)}%" for p in [0.25, 0.5, 0.75]}

    for k in range(Nt - 1):
        t     = t_vec[k]
        t_mid = t + dt / 2.0

        if k in progress_pts:
            print(f"  Progress: {progress_pts[k]}")

        Ta_C     = _interp_py(met['temp'],   t_mid)
        Isolar   = _interp_py(met['solar'],  t_mid)
        Prain_mh = _interp_py(met['precip'], t_mid)
        wind     = _interp_py(met['wind'],   t_mid)
        RH_pct   = _interp_py(met['rh'],     t_mid)
        Tsoil_C  = _interp_py(met['soil'],   t_mid)

        Ta_K   = Ta_C + 273.15
        Prain  = Prain_mh / 3600.0
        RH_frac = RH_pct / 100.0
        Tsoil_K = Tsoil_C + 273.15

        h_out = compute_h_out(wind)

        forc = {
            "Isolar": Isolar,
            "Prain":  Prain,
            "T_rain": Ta_K,
            "RH":     RH_frac,
            "Ta":     Ta_K,
            "U10":    wind,
            "h_out":  h_out,
        }

        # Advanced insulation step
        R_ins, q_solar_ins, q_rain_snow, q_evap, InsState = \
            insulation_step(InsState, forc, _INS_PAR, dt)

        W_hist[k]     = InsState["W"]
        k_eff_hist[k] = InsState["k_eff"]
        alpha_hist[k] = InsState["alpha_eff"]

        R_a2s = R_eff + R_ins

        # RK4 temperature integration (Numba-JIT)
        T_new = rk4_step(
            t, T, dt,
            R_a2s, q_solar_ins, q_rain_snow, q_evap,
            temp_arr, soil_arr, dt_data,
            R_12, R_23, Cs_layer,
            h_ground, Tfreeze,
        )

        # Fluxes for diagnostics
        q_a_mid    = (Ta_K - T[0]) / R_a2s
        q_surf_mid = q_a_mid + q_solar_ins + q_rain_snow + q_evap
        R_cond     = 1.0 / 1.5
        R_intf     = 1.0 / h_ground
        q_gnd_mid  = (Tsoil_K - T[2]) / (R_cond + R_intf)

        # Refreezing
        for i in range(3):
            T_new[i], LWC[i], ice_fractions[i], _ = \
                refreezing_layer(T_new[i], LWC[i], ice_fractions[i])

        # Surface melt
        if T_new[0] > Tfreeze:
            T_new[0] = Tfreeze
            if q_surf_mid > 0.0:
                dE   = q_surf_mid * dt
                dM   = dE / (rho_i * Lf)
                E_melt += dE
                melt_rate_h[k] = dM / dt
                LWC[0] += dM / dz_s

        # Percolation
        LWC, runoff = percolate_water(LWC, heights, theta_e)
        runoff_hist[k] = runoff

        # Store
        T = T_new
        T_hist[k+1]   = T
        LWC_hist[k+1] = LWC
        qa_hist[k]    = q_a_mid
        qsolar_hist[k]= q_solar_ins
        qrain_hist[k] = q_rain_snow
        qevap_hist[k] = q_evap
        qground_hist[k] = q_gnd_mid

    # Fill last step
    W_hist[-1]     = InsState["W"]
    k_eff_hist[-1] = InsState.get("k_eff", k_eff_hist[-2])
    alpha_hist[-1] = InsState.get("alpha_eff", alpha_hist[-2])

    print("  Progress: 100%")
    print(f"  Done.  Total melt = {E_melt/(rho_i*Lf)*rho_w*1000:.1f} mm w.e.")

    # SWE trajectory (initial minus cumulative ice lost)
    melt_cumul_m  = np.cumsum(melt_rate_h) * dt       # m ice
    SWE_hist      = np.maximum(0.0, Hs * rho_s - melt_cumul_m * rho_i)
    runoff_cumul  = np.cumsum(runoff_hist)

    t0    = datetime(2024, 4, 1)
    dates = [t0 + timedelta(seconds=float(s)) for s in t_vec]

    return {
        'dates':        dates,
        'T_hist':       T_hist,        # K  (Nt, 3)
        'LWC_hist':     LWC_hist,      # -  (Nt, 3)
        'SWE':          SWE_hist,      # kg/m²
        'melt_cumul':   melt_cumul_m,  # m ice equivalent
        'runoff_cumul': runoff_cumul,  # kg/m²
        'q_ground':     qground_hist,  # W/m²
        'qa':           qa_hist,
        'q_solar':      qsolar_hist,
        'q_rain':       qrain_hist,
        'q_evap':       qevap_hist,
        'k_eff':        k_eff_hist,
        'alpha':        alpha_hist,
        'W':            W_hist,
        't_vec':        t_vec,
    }


# ============================================================
#  10.  Comparison plot
# ============================================================

def make_comparison_plot(py, sp, out_path):
    """
    py  : dict returned by run_python_sim()
    sp  : pd.DataFrame returned by read_snowpack_met()
    """
    step  = max(1, int(3600.0 / 600.0))   # downsample to ~hourly
    idx   = np.arange(0, len(py['dates']), step)
    d_py  = [py['dates'][i] for i in idx]

    sp_runoff_col = 'Snowpack runoff (virtual lysimeter -- snow only)'
    sp_runoff_cumul = sp[sp_runoff_col].fillna(0).cumsum() \
        if sp_runoff_col in sp.columns else None

    fig, axes = plt.subplots(5, 1, figsize=(14, 20), sharex=True)
    fig.suptitle(
        'SNOWPACK vs Python snowsim_v4 — Snow Storage 2024\n'
        'Python: advanced insulation + refreezing + percolation | SNOWPACK: full physics',
        fontsize=13, fontweight='bold', y=0.995,
    )

    # ---- Panel 1: layer temperatures ----
    ax = axes[0]
    ax.plot(d_py, py['T_hist'][idx, 0] - 273.15, '#aec7e8', lw=1.2, label='Py T1 (0.75 m)')
    ax.plot(d_py, py['T_hist'][idx, 1] - 273.15, '#ffbb78', lw=1.2, label='Py T2 (2.25 m)')
    ax.plot(d_py, py['T_hist'][idx, 2] - 273.15, '#98df8a', lw=1.2, label='Py T3 (3.75 m)')
    for col, color, label in [
        ('Temperature 1 (modelled)', '#1f77b4', 'SP T1 (0.75 m)'),
        ('Temperature 2 (modelled)', '#d62728', 'SP T2 (2.25 m)'),
        ('Temperature 3 (modelled)', '#2ca02c', 'SP T3 (3.75 m)'),
    ]:
        if col in sp.columns:
            ax.plot(sp.index, sp[col], color, lw=2, label=label)
    if 'Modeled surface temperature' in sp.columns:
        ax.plot(sp.index, sp['Modeled surface temperature'],
                'k', lw=0.8, ls='--', alpha=0.4, label='SP surface')
    ax.axhline(0, color='gray', lw=0.8, ls=':')
    ax.set_ylabel('Temperature [°C]')
    ax.set_ylim(-20, 5)
    ax.legend(fontsize=7, ncol=4, loc='lower right')
    ax.set_title('Layer Temperatures')
    ax.grid(True, alpha=0.3)

    # ---- Panel 2: SWE ----
    ax = axes[1]
    ax.plot(d_py, py['SWE'][idx], 'steelblue', lw=2.5, label='Python SWE')
    if 'SWE (of snowpack)' in sp.columns:
        ax.plot(sp.index, sp['SWE (of snowpack)'], 'navy', lw=2, ls='--', label='SNOWPACK SWE')
    if 'Modelled snow depth (vertical)' in sp.columns:
        ax2b = ax.twinx()
        ax2b.plot(sp.index, sp['Modelled snow depth (vertical)'] / 100,
                  'cornflowerblue', lw=1, alpha=0.4, label='SP depth (m)')
        ax2b.set_ylabel('Depth [m]', color='cornflowerblue')
        ax2b.tick_params(axis='y', labelcolor='cornflowerblue')
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax2b.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, fontsize=8, loc='upper right')
    else:
        ax.legend(fontsize=8, loc='upper right')
    ax.set_ylabel('SWE [kg/m²]')
    ax.set_title('Snow Water Equivalent & Depth')
    ax.grid(True, alpha=0.3)

    # ---- Panel 3: cumulative melt / runoff ----
    ax = axes[2]
    py_melt_mm = py['melt_cumul'][idx] * rho_w  # m_ice → mm w.e.  (rho_w = 1000 kg/m³)
    ax.plot(d_py, py_melt_mm, 'tomato', lw=2,
            label=f'Python cumul. melt  ({py_melt_mm[-1]:.0f} mm w.e.)')
    if sp_runoff_cumul is not None:
        ax.plot(sp.index, sp_runoff_cumul, 'darkred', lw=2, ls='--',
                label=f'SNOWPACK cumul. runoff  ({sp_runoff_cumul.max():.0f} kg/m²)')

    ax.set_ylabel('mm w.e.  (Python) / kg/m²  (SNOWPACK)')
    ax.legend(fontsize=8, loc='upper left')
    ax.set_title('Cumulative Melt / Runoff')
    ax.grid(True, alpha=0.3)

    # ---- Panel 4: ground heat flux ----
    ax = axes[3]
    ax.plot(d_py, py['q_ground'][idx], 'sienna', lw=1.5,
            label=f'Python Robin BC  (mean={py["q_ground"].mean():.2f} W/m²)')
    ax.fill_between(d_py, 0, py['q_ground'][idx], alpha=0.12, color='sienna')
    ax.axhline(2.414, color='goldenrod', lw=2, ls='--',
               label='SNOWPACK GEO_HEAT = 2.414 W/m²')
    ax.axhline(0, color='gray', lw=0.8, ls=':')
    ax.set_ylabel('W/m²')
    ax.legend(fontsize=8)
    ax.set_title('Ground Heat Flux')
    ax.text(0.01, 0.03,
            f'Python: R_total = L_soil/k_soil + 1/h_ground = {1/1.5:.3f} + {1/h_ground:.3f}'
            f' = {1/1.5 + 1/h_ground:.3f} m²K/W,  driven by Soil_Temp_320cm',
            transform=ax.transAxes, fontsize=7, style='italic', color='gray', va='bottom')
    ax.grid(True, alpha=0.3)

    # ---- Panel 5: insulation effective conductivity ----
    ax = axes[4]
    ax.plot(d_py, py['k_eff'][idx], 'b-', lw=1.5,
            label=f'k_eff  (mean={py["k_eff"].mean():.4f} W/(mK))')
    ax_r = ax.twinx()
    ax_r.plot(d_py, py['alpha'][idx], 'r-', lw=1.5, alpha=0.7,
              label=f'α_eff  (mean={py["alpha"].mean():.3f})')
    ax.set_ylabel('k_eff [W/(mK)]', color='b')
    ax_r.set_ylabel('α_eff [–]', color='r')
    ax.tick_params(axis='y', labelcolor='b')
    ax_r.tick_params(axis='y', labelcolor='r')
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax_r.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=8, loc='upper left')
    ax.set_title('Advanced Insulation Properties (Python only)')
    ax.grid(True, alpha=0.3)

    for ax in axes:
        ax.xaxis.set_major_locator(mdates.MonthLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
        ax.xaxis.set_minor_locator(mdates.WeekdayLocator(interval=2))
    axes[-1].set_xlabel('2024')

    plt.tight_layout(rect=[0, 0, 1, 0.97])
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    print(f"Saved → {out_path}")


# ============================================================
#  11.  Summary table
# ============================================================

def print_summary(py, sp):
    sp_runoff_col = 'Snowpack runoff (virtual lysimeter -- snow only)'
    sp_runoff_total = sp[sp_runoff_col].fillna(0).sum() \
        if sp_runoff_col in sp.columns else float('nan')

    py_snowfree = next(
        (py['dates'][i] for i, v in enumerate(py['SWE']) if v <= 0), None)
    if 'SWE (of snowpack)' in sp.columns:
        mask = sp['SWE (of snowpack)'] <= 0
        sp_snowfree = sp.index[mask][0] if mask.any() else None
    else:
        sp_snowfree = None

    print("\n=== Comparison Summary (snowsim_v4 vs SNOWPACK) ===")
    print(f"{'Metric':<45} {'Python v4':>12} {'SNOWPACK':>12}")
    print("-" * 71)

    def _sp(col, fmt='.0f'):
        return (f"{sp[col].iloc[0]:{fmt}}" if col in sp.columns else "N/A")

    rows = [
        ("Initial SWE (kg/m²)",
         f"{py['SWE'][0]:.0f}",
         _sp('SWE (of snowpack)')),
        ("Final SWE (kg/m²)",
         f"{py['SWE'][-1]:.0f}",
         _sp('SWE (of snowpack)', '.1f')),  # uses iloc[0] — placeholder
        ("Cumul. melt  (mm w.e.)",
         f"{py['melt_cumul'][-1] * rho_w:.0f}",
         f"{sp_runoff_total:.0f}"),
        ("Mean ground flux (W/m²)",
         f"{py['q_ground'].mean():.3f}",
         "2.414 (const)"),
        ("Mean k_eff insulation (W/(mK))",
         f"{py['k_eff'].mean():.4f}",
         "—"),
        ("Mean alpha_eff insulation",
         f"{py['alpha'].mean():.4f}",
         "—"),
        ("Snow-free date",
         str(py_snowfree)[:10] if py_snowfree else "—",
         str(sp_snowfree)[:10] if sp_snowfree else "—"),
    ]
    for label, py_val, sp_val in rows:
        print(f"{label:<45} {py_val:>12} {sp_val:>12}")


# ============================================================
#  12.  Entry point
# ============================================================

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Compare SNOWPACK vs Python snowsim_v4 RC model')
    parser.add_argument('--csv', default='DATA_2024.csv',
                        help='Path to meteorological CSV (default: DATA_2024.csv)')
    parser.add_argument('--met', default='output/snow_storage_snow_storage.met',
                        help='Path to SNOWPACK .met output file')
    parser.add_argument('--out', default='snowpack_vs_python_v4.png',
                        help='Output plot filename (default: snowpack_vs_python_v4.png)')
    parser.add_argument('--dt',  type=float, default=600.0,
                        help='Python model timestep [s] (default: 600)')
    args = parser.parse_args()

    print(f"Reading CSV:  {args.csv}")
    met = read_met_csv(args.csv)
    print(f"  {len(met['temp'])} hourly records")

    py = run_python_sim(met, dt=args.dt)

    print(f"\nReading SNOWPACK: {args.met}")
    sp, units = read_snowpack_met(args.met)
    print(f"  {len(sp)} rows,  {sp.index[0]}  →  {sp.index[-1]}")

    print_summary(py, sp)
    make_comparison_plot(py, sp, args.out)