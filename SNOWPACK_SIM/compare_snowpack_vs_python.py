"""
compare_snowpack_vs_python.py
------------------------------
Runs a simplified version of snowsim_v3.py and plots it against
SNOWPACK .met output side by side.

Usage:
    python compare_snowpack_vs_python.py \
        --csv   DATA_2024.csv \
        --met   output/snow_storage_snow_storage.met \
        --out   comparison.png

Requires: numpy, pandas, matplotlib
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


# ============================================================
#  1.  SNOWPACK .met reader
# ============================================================

def read_snowpack_met(met_path):
    """
    Parse a SNOWPACK .met file into a DataFrame.

    The file format is:
        [STATION_PARAMETERS]  ... metadata ...
        [HEADER]
        ,,1,2,3,...            <- column numbers (ignored)
        ID,Date,Name1,...      <- column names
        ,,unit1,...            <- units (ignored)
        [DATA]
        0203,dd.mm.yyyy HH:MM:SS, val1, val2, ...

    Returns
    -------
    df : pd.DataFrame  indexed by datetime, columns named from the header row.
    col_units : dict   {column_name: unit_string}
    """
    with open(met_path) as f:
        lines = f.readlines()

    # Find the three header lines (start with ',,', 'ID,', ',,')
    header_block = [l for l in lines if l.startswith(',,') or l.startswith('ID,Date,')]
    col_numbers_line = header_block[0]   # ,,1,2,3,...
    col_names_line   = header_block[1]   # ID,Date,Name,...
    col_units_line   = header_block[2]   # ,,W m-2,...

    col_names = [c.strip() for c in col_names_line.split(',')]
    col_units_raw = [c.strip() for c in col_units_line.split(',')]
    col_units = {col_names[i]: col_units_raw[i]
                 for i in range(min(len(col_names), len(col_units_raw)))}

    # Data rows start with the station code (e.g. '0203,')
    data_rows = [l.strip().split(',') for l in lines if l[:4].isdigit()]

    # Deduplicate column names (SNOWPACK uses '-' for unnamed columns)
    seen = {}
    deduped = []
    for name in col_names[:len(data_rows[0])]:
        if name in seen:
            seen[name] += 1
            deduped.append(f"{name}_{seen[name]}")
        else:
            seen[name] = 0
            deduped.append(name)

    df = pd.DataFrame(data_rows, columns=deduped)

    # Parse datetime (format: "01.04.2024 00:00:00")
    df['Date'] = pd.to_datetime(df['Date'].str.strip(), format='%d.%m.%Y %H:%M:%S')
    df = df.set_index('Date').sort_index()

    # Convert all value columns to numeric, replacing -999 nodata with NaN
    for c in df.columns:
        if c not in ('ID',):
            df[c] = pd.to_numeric(df[c], errors='coerce').replace(-999.0, np.nan)

    return df, col_units


# ============================================================
#  2.  Simplified Python snow model  (snowsim_v3 core, no Numba)
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
    return {k: np.array(v) for k, v in data.items()}


def interp(vec, t_s, dt_data=3600.0):
    """Linear interpolation of hourly data vector at time t_s [seconds]."""
    idx = t_s / dt_data
    lo = int(np.floor(idx))
    hi = min(lo + 1, len(vec) - 1)
    lo = max(lo, 0)
    frac = idx - lo
    return vec[lo] * (1 - frac) + vec[hi] * frac


def run_python_sim(met, dt=600.0):
    """
    Simplified 3-layer RC snow model from snowsim_v3.py.
    Advanced insulation and Numba removed; constant insulation resistance used.

    USE_ADVANCED_INSULATION = False
    USE_REFREEZING = False
    USE_PERCOLATION = True
    USE_MULTILAYER_INSULATION = False

    Parameters
    ----------
    met : dict of np.ndarray   hourly forcing arrays
    dt  : float                model timestep [s]

    Returns
    -------
    dict with keys: dates, T_hist, SWE, melt_cumul, runoff_cumul, q_ground
    """
    # Physical constants
    Lf = 3.34e5; rho_i = 917.0; rho_s = 400.0
    rho_w = 1000.0; c_w = 4180.0; c_s = 2100.0; Tfreeze = 273.15
    sigma = 5.670374419e-8

    # Snow geometry
    Hs = 4.5; Ns = 3; dz_s = Hs / Ns
    k_snow = 0.731665
    R_layer = dz_s / k_snow

    # Insulation (constant, simplified)
    Hi = 0.20; k_i_base = 0.32
    R_ins = Hi / k_i_base

    # Surface heat transfer
    h_conv = 8.0; epsilon = 0.95; T_mean = 276.15
    h_rad = 4 * epsilon * sigma * T_mean**3
    R_eff = 1.0 / (h_conv + h_rad)
    R_a2s = R_eff + R_ins

    # Ground Robin BC  (R = L/k_soil + 1/h_ground)
    h_ground = 3.675021
    R_ground = 1.0 / 1.5 + 1.0 / h_ground

    # Layer heat capacity
    Cs = rho_s * c_s * dz_s

    # Liquid water
    theta_e = 0.055119
    alpha = 0.80         # constant solar absorptivity (no advanced insulation)

    # Time vector
    n_hours = len(met['temp'])
    t_vec = np.arange(0.0, n_hours * 3600.0, dt)
    Nt = len(t_vec)

    # Initial conditions  (matches .sno file)
    T = np.array([271.15, 269.15, 267.15])
    LWC = np.zeros(3)

    # Output arrays
    T_hist       = np.zeros((Nt, 3))
    SWE_hist     = np.zeros(Nt)
    melt_cumul   = np.zeros(Nt)
    runoff_cumul = np.zeros(Nt)
    q_ground_hist = np.zeros(Nt)

    T_hist[0] = T
    SWE_hist[0] = Hs * rho_s
    cum_melt = 0.0
    cum_runoff = 0.0

    def dTdt(t, Tv):
        T1, T2, T3 = Tv
        Ta = interp(met['temp'], t) + 273.15
        Is = interp(met['solar'], t)
        P  = interp(met['precip'], t) / 3600.0   # m/s
        q_a    = (Ta - T1) / R_a2s
        q_sol  = alpha * Is
        q_rain = rho_w * c_w * P * (Ta - Tfreeze)
        q_surf = q_a + q_sol + q_rain
        Tsoil  = interp(met['soil'], t) + 273.15
        q_3g   = (Tsoil - T3) / R_ground
        dT1 = (q_surf + (T2 - T1) / R_layer) / Cs
        dT2 = ((T1 - T2) / R_layer + (T3 - T2) / R_layer) / Cs
        dT3 = ((T2 - T3) / R_layer + q_3g) / Cs
        return np.array([dT1, dT2, dT3])

    for k in range(Nt - 1):
        t = t_vec[k]; tm = t + dt / 2.0

        # RK4 temperature update
        k1 = dTdt(t, T)
        k2 = dTdt(tm, T + dt * k1 / 2)
        k3 = dTdt(tm, T + dt * k2 / 2)
        k4 = dTdt(t + dt, T + dt * k3)
        Tn = T + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)

        # Surface melt
        Ta_K = interp(met['temp'], tm) + 273.15
        Is   = interp(met['solar'], tm)
        P    = interp(met['precip'], tm) / 3600.0
        q_a  = (Ta_K - T[0]) / R_a2s
        q_surf = q_a + alpha * Is + rho_w * c_w * P * (Ta_K - Tfreeze)
        if Tn[0] >= Tfreeze and q_surf > 0:
            Tn[0] = Tfreeze
            dm = q_surf * dt / (rho_i * Lf)
            cum_melt += dm
            LWC[0] += dm / dz_s

        SWE = max(0.0, Hs * rho_s - cum_melt * rho_i)

        # Percolation (bucket)
        for i in range(2):
            if LWC[i] > theta_e:
                exc = LWC[i] - theta_e; LWC[i] = theta_e; LWC[i+1] += exc
        runoff = 0.0
        if LWC[2] > theta_e:
            exc = LWC[2] - theta_e; LWC[2] = theta_e
            runoff = exc * dz_s * rho_w
            cum_runoff += runoff

        # Ground flux for diagnostics
        Tsoil_K = interp(met['soil'], tm) + 273.15
        q_g = (Tsoil_K - Tn[2]) / R_ground

        T = Tn
        T_hist[k+1]        = T
        SWE_hist[k+1]      = SWE
        melt_cumul[k+1]    = cum_melt
        runoff_cumul[k+1]  = cum_runoff
        q_ground_hist[k+1] = q_g

    t0 = datetime(2024, 4, 1)
    dates = [t0 + timedelta(seconds=float(s)) for s in t_vec]

    return {
        'dates':        dates,
        'T_hist':       T_hist,         # K, shape (Nt, 3)
        'SWE':          SWE_hist,       # kg/m²
        'melt_cumul':   melt_cumul,     # m ice equivalent
        'runoff_cumul': runoff_cumul,   # kg/m²
        'q_ground':     q_ground_hist,  # W/m²
    }


# ============================================================
#  3.  Plot
# ============================================================

def make_comparison_plot(py, sp, out_path):
    """
    py  : dict returned by run_python_sim()
    sp  : pd.DataFrame returned by read_snowpack_met()
    """
    # Downsample Python output to hourly for cleaner lines
    step = 6   # 6 × 10 min = 1 hour
    idx  = np.arange(0, len(py['dates']), step)
    d_py = [py['dates'][i] for i in idx]

    # SNOWPACK cumulative runoff (column 38 is a rate per output interval → cumsum)
    sp_runoff_cumul = sp['Snowpack runoff (virtual lysimeter -- snow only)'].fillna(0).cumsum()

    fig, axes = plt.subplots(4, 1, figsize=(14, 16), sharex=True)
    fig.suptitle(
        'SNOWPACK vs Python snowsim_v3 — Snow Storage Apr–Aug 2024\n'
        'Python: simplified constant insulation | SNOWPACK: full physics',
        fontsize=13, fontweight='bold', y=0.99
    )

    # ---- Panel 1: layer temperatures ----
    ax = axes[0]
    # Python (light shades)
    ax.plot(d_py, py['T_hist'][idx, 0] - 273.15, '#aec7e8', lw=1.2, label='Py T1 (0.75 m)')
    ax.plot(d_py, py['T_hist'][idx, 1] - 273.15, '#ffbb78', lw=1.2, label='Py T2 (2.25 m)')
    ax.plot(d_py, py['T_hist'][idx, 2] - 273.15, '#98df8a', lw=1.2, label='Py T3 (3.75 m)')
    # SNOWPACK (bold)
    ax.plot(sp.index, sp['Temperature 1 (modelled)'], '#1f77b4', lw=2, label='SP T1 (0.75 m)')
    ax.plot(sp.index, sp['Temperature 2 (modelled)'], '#d62728', lw=2, label='SP T2 (2.25 m)')
    ax.plot(sp.index, sp['Temperature 3 (modelled)'], '#2ca02c', lw=2, label='SP T3 (3.75 m)')
    ax.plot(sp.index, sp['Modeled surface temperature'], 'k', lw=0.8, ls='--', alpha=0.4,
            label='SP surface')
    ax.axhline(0, color='gray', lw=0.8, ls=':')
    ax.set_ylabel('Temperature [°C]'); ax.set_ylim(-20, 5)
    ax.legend(fontsize=7, ncol=4, loc='lower right')
    ax.set_title('Layer Temperatures (SP values go NaN once depth sensor exits melting snowpack)')
    ax.grid(True, alpha=0.3)

    # ---- Panel 2: SWE + snow depth ----
    ax = axes[1]
    ax.plot(d_py, py['SWE'][idx], 'steelblue', lw=2.5, label='Python SWE')
    ax.plot(sp.index, sp['SWE (of snowpack)'], 'navy', lw=2, ls='--', label='SNOWPACK SWE')
    ax2b = ax.twinx()
    ax2b.plot(sp.index, sp['Modelled snow depth (vertical)'] / 100,
              'cornflowerblue', lw=1, alpha=0.4, label='SP depth (m)')
    ax.set_ylabel('SWE [kg/m²]')
    ax2b.set_ylabel('Depth [m]', color='cornflowerblue')
    ax2b.tick_params(axis='y', labelcolor='cornflowerblue')
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2b.get_legend_handles_labels()
    ax.legend(h1+h2, l1+l2, fontsize=8, loc='upper right')
    ax.set_title('Snow Water Equivalent & Depth')
    ax.grid(True, alpha=0.3)

    # ---- Panel 3: cumulative melt / runoff ----
    ax = axes[2]
    ax.plot(d_py, py['melt_cumul'][idx] * 917,
            'tomato', lw=2, label='Python cumul. melt [kg/m²]')
    ax.plot(sp.index, sp_runoff_cumul,
            'darkred', lw=2, ls='--',
            label=f'SNOWPACK cumul. runoff [kg/m²]  (total={sp_runoff_cumul.max():.0f})')
    ax.set_ylabel('kg/m²')
    ax.legend(fontsize=8, loc='upper left')
    ax.set_title('Cumulative Melt / Runoff')
    ax.grid(True, alpha=0.3)

    # ---- Panel 4: ground heat flux ----
    ax = axes[3]
    ax.plot(d_py, py['q_ground'][idx], 'sienna', lw=1.5,
            label=f'Python Robin BC  (mean={py["q_ground"].mean():.1f} W/m²)')
    ax.fill_between(d_py, 0, py['q_ground'][idx], alpha=0.12, color='sienna')
    ax.axhline(0.06, color='goldenrod', lw=2, ls='--',
               label='SNOWPACK GEO_HEAT=0.06 W/m²  (SNP_SOIL=false)')
    ax.axhline(0, color='gray', lw=0.8, ls=':')
    ax.set_ylabel('W/m²'); ax.set_ylim(-1, 12)
    ax.legend(fontsize=8)
    ax.set_title('Ground Heat Flux  →  set SNP_SOIL=true to activate TSG-driven BC in SNOWPACK')
    ax.text(0.01, 0.05,
            'Python: R_total = L_soil/k_soil + 1/h_ground = 0.667+0.272 = 0.939 m²K/W, '
            'driven by Soil_Temp_320cm\n'
            'SNOWPACK: constant GEO_HEAT because SNP_SOIL=false bypasses TSG column',
            transform=ax.transAxes, fontsize=7, style='italic', color='gray', va='bottom')
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
#  4.  Summary table
# ============================================================

def print_summary(py, sp):
    sp_runoff_total = sp['Snowpack runoff (virtual lysimeter -- snow only)'].fillna(0).sum()
    py_snowfree = next(
        (py['dates'][i] for i, v in enumerate(py['SWE']) if v <= 0), None)
    sp_snowfree = sp[sp['SWE (of snowpack)'] <= 0].index[0] \
        if (sp['SWE (of snowpack)'] <= 0).any() else None

    print("\n=== Comparison Summary ===")
    print(f"{'Metric':<40} {'Python':>12} {'SNOWPACK':>12}")
    print("-" * 66)
    rows = [
        ("Initial SWE (kg/m²)",      f"{py['SWE'][0]:.0f}",           f"{sp['SWE (of snowpack)'].iloc[0]:.0f}"),
        ("Final SWE (kg/m²)",        f"{py['SWE'][-1]:.0f}",          f"{sp['SWE (of snowpack)'].iloc[-1]:.1f}"),
        ("Cumul. melt/runoff (kg/m²)",f"{py['melt_cumul'][-1]*917:.0f}", f"{sp_runoff_total:.0f}"),
        ("Min T at 0.75 m (°C)",     f"{(py['T_hist'][:,0]-273.15).min():.2f}",
                                      f"{sp['Temperature 1 (modelled)'].min():.2f}"),
        ("Mean ground flux (W/m²)",  f"{py['q_ground'].mean():.2f}",  "0.06 (const)"),
        ("Snow-free date",           str(py_snowfree)[:10] if py_snowfree else "—",
                                      str(sp_snowfree)[:10] if sp_snowfree else "—"),
    ]
    for label, py_val, sp_val in rows:
        print(f"{label:<40} {py_val:>12} {sp_val:>12}")


# ============================================================
#  Main
# ============================================================

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compare SNOWPACK vs Python snow model')
    parser.add_argument('--csv', default='DATA_2024.csv',
                        help='Path to meteorological CSV (default: DATA_2024.csv)')
    parser.add_argument('--met', default='output/snow_storage_snow_storage.met',
                        help='Path to SNOWPACK .met output file')
    parser.add_argument('--out', default='snowpack_vs_python.png',
                        help='Output plot filename (default: snowpack_vs_python.png)')
    parser.add_argument('--dt',  type=float, default=600.0,
                        help='Python model timestep in seconds (default: 600)')
    args = parser.parse_args()

    print(f"Reading CSV:  {args.csv}")
    met = read_met_csv(args.csv)
    print(f"  {len(met['temp'])} hourly records")

    print(f"Running Python sim (dt={args.dt}s)...")
    py = run_python_sim(met, dt=args.dt)
    print(f"  Done. SWE {py['SWE'][0]:.0f} → {py['SWE'][-1]:.0f} kg/m²")

    print(f"Reading SNOWPACK: {args.met}")
    sp, units = read_snowpack_met(args.met)
    print(f"  {len(sp)} rows, {sp.index[0]} to {sp.index[-1]}")

    print_summary(py, sp)
    make_comparison_plot(py, sp, args.out)
