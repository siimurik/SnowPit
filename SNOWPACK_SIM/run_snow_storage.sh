#!/usr/bin/env bash
# run_snow_storage.sh
# Run SNOWPACK simulation for snow_storage configuration
# Replicates snowsim_v3.py melting conditions using DATA_2024.csv

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# ── Output directory ──────────────────────────────────────────
mkdir -p output

# ── Locate SNOWPACK binary ────────────────────────────────────
SNOWPACK_BIN=""
for candidate in \
    snowpack \
    /usr/local/bin/snowpack \
    /opt/snowpack/bin/snowpack \
    "$HOME/.local/bin/snowpack"; do
    if command -v "$candidate" &>/dev/null 2>&1; then
        SNOWPACK_BIN="$candidate"
        break
    fi
done

if [ -z "$SNOWPACK_BIN" ]; then
    echo "ERROR: 'snowpack' binary not found."
    echo "  Install from: https://gitlabext.wsl.ch/snow-models/snowpack"
    echo "  Or via conda:  conda install -c conda-forge snowpack"
    exit 1
fi

echo "=============================================="
echo " SNOWPACK snow_storage simulation"
echo " Replicating snowsim_v3.py melting conditions"
echo " Period: 2024-04-01 to 2024-08-31"
echo "=============================================="
echo "Using binary: $SNOWPACK_BIN"
echo ""

# ── Run ───────────────────────────────────────────────────────
"$SNOWPACK_BIN" -c snow_storage.ini -e 2024-08-31T23:00

echo ""
echo "Simulation complete. Results written to: output/"
echo ""
echo "Output files:"
ls -lh output/ 2>/dev/null || echo "  (none yet)"
