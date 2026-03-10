"""
csv_to_smet.py
--------------
Convert DATA_2024.csv to a SNOWPACK-compatible .smet file.

Usage:
    python csv_to_smet.py                                      # defaults
  > python csv_to_smet.py --csv DATA_2024.csv --out input/snow_storage.smet
    python csv_to_smet.py --csv DATA_2024.csv --out input/snow_storage.smet \\
        --lat 59.39805 --lon 24.60277 --alt 33.16 \\
        --station-id snow_storage --station-name "Snow Storage Site" \\
        --tz 0 --step 10

What this script does
---------------------
1. Reads the CSV (hourly, m/h precipitation).
2. Converts units to SNOWPACK conventions:
       TA   [K]         = Temp_C + 273.15
       RH   [0–1]       = RH_% / 100
       VW   [m/s]       = Air_Vel_m/s_10m   (no change)
       ISWR [W/m²]      = Glo_Sol_Ir_W/m2   (no change)
       PSUM [kg/m²]     = Prec_m/h * 1000   (m/h → kg/m²/h = mm/h)
                          NOTE: kept as hourly accumulation so the
                          SNOWPACK accumulate resampler with period=<step>
                          can distribute it correctly.
       TSG  [K]         = Soil_Temp_320cm + 273.15
3. Prepends one anchor row (first CSV row, timestamp shifted back by
   one hour, PSUM forced to 0) so that SNOWPACK's accumulate resampler
   has a prior data point at the very first model step.
4. Writes the SMET file.

The anchor row is the reason the current .smet starts at 2024-03-31T23:00
even though the CSV starts at 2024-04-01T00:00.  It is intentional and
required — without it SNOWPACK reports "missing precipitation" on the
first timestep because the accumulate resampler needs two consecutive
timestamps to compute an increment.
"""

import argparse
import csv
from datetime import datetime, timedelta


# ── unit conversion ───────────────────────────────────────────────────────────

def to_smet_row(row: dict) -> dict:
    """Convert one CSV row to SMET field values (all units, no formatting yet)."""
    return {
        'TA':   float(row['Temp_C']) + 273.15,
        'RH':   float(row['RH_%']) / 100.0,
        'VW':   float(row['Air_Vel_m/s_10m']),
        'ISWR': float(row['Glo_Sol_Ir_W/m2']),
        # m/h × 1000 → kg/m²/h  (= mm/h; hourly accumulation for accumulate resampler)
        'PSUM': float(row['Prec_m/h']) * 1000.0,
        'TSG':  float(row['Soil_Temp_320cm']) + 273.15,
    }


def format_row(timestamp: str, vals: dict) -> str:
    return (f"{timestamp} "
            f"{vals['TA']:.2f} "
            f"{vals['RH']:.3f} "
            f"{vals['VW']:.2f} "
            f"{vals['ISWR']:.1f} "
            f"{vals['PSUM']:.5f} "
            f"{vals['TSG']:.2f}")


# ── main ──────────────────────────────────────────────────────────────────────

def convert(csv_path, out_path, lat, lon, alt, station_id, station_name, tz):

    # Read CSV
    rows = []
    with open(csv_path, newline='') as f:
        for r in csv.DictReader(f):
            rows.append(r)

    if not rows:
        raise ValueError(f"No data rows found in {csv_path}")

    print(f"Read {len(rows)} rows from {csv_path}")
    print(f"  Period: {rows[0]['Time']}  →  {rows[-1]['Time']}")

    # Build SMET data lines
    data_lines = []

    # ── anchor row ────────────────────────────────────────────────────────────
    # Duplicate of the first real row but timestamped one hour earlier and with
    # PSUM=0.  This gives the accumulate resampler a "previous" value so it can
    # compute an increment at the very first model timestep.
    first_ts = datetime.strptime(rows[0]['Time'], '%Y-%m-%dT%H:%M')
    anchor_ts = (first_ts - timedelta(hours=1)).strftime('%Y-%m-%dT%H:%M')
    anchor_vals = to_smet_row(rows[0])
    anchor_vals['PSUM'] = 0.0   # no precip before data starts
    data_lines.append(format_row(anchor_ts, anchor_vals))

    # ── real data ─────────────────────────────────────────────────────────────
    for r in rows:
        data_lines.append(format_row(r['Time'], to_smet_row(r)))

    # ── write SMET ────────────────────────────────────────────────────────────
    header = f"""\
SMET 1.1 ASCII
[HEADER]
station_id       = {station_id}
station_name     = {station_name}
latitude         = {lat}
longitude        = {lon}
altitude         = {alt}
nodata           = -999
tz               = {tz}
source           = {csv_path}
fields           = timestamp TA RH VW ISWR PSUM TSG
[DATA]"""

    with open(out_path, 'w') as f:
        f.write(header + '\n')
        f.write('\n'.join(data_lines) + '\n')

    print(f"Written {len(data_lines)} rows (1 anchor + {len(rows)} data) → {out_path}")
    print(f"  First row: {data_lines[0][:60]}")
    print(f"  Second row: {data_lines[1][:60]}")
    print(f"  Last row:  {data_lines[-1][:60]}")

    # Quick sanity check
    non_zero_psum = sum(1 for l in data_lines[1:] if float(l.split()[5]) > 0)
    print(f"  Non-zero PSUM rows: {non_zero_psum}")


# ── CLI ───────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Convert DATA_2024.csv to a SNOWPACK .smet file')

    parser.add_argument('--csv',          default='DATA_2024.csv',
                        help='Input CSV file (default: DATA_2024.csv)')
    parser.add_argument('--out',          default='input/snow_storage.smet',
                        help='Output .smet file (default: input/snow_storage.smet)')

    # Site metadata
    parser.add_argument('--lat',          type=float, default=59.39805555555556,
                        help='Latitude [decimal degrees]')
    parser.add_argument('--lon',          type=float, default=24.602777777777778,
                        help='Longitude [decimal degrees]')
    parser.add_argument('--alt',          type=float, default=33.16,
                        help='Altitude [m a.s.l.]')
    parser.add_argument('--station-id',   default='snow_storage',
                        help='station_id header field')
    parser.add_argument('--station-name', default='Snow Storage Site',
                        help='station_name header field')
    parser.add_argument('--tz',           type=int,   default=0,
                        help='Time zone offset from UTC (default: 0)')

    args = parser.parse_args()

    convert(
        csv_path     = args.csv,
        out_path     = args.out,
        lat          = args.lat,
        lon          = args.lon,
        alt          = args.alt,
        station_id   = args.station_id,
        station_name = args.station_name,
        tz           = args.tz,
    )