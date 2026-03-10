```
snowpack -c snow_storage.ini -e 2024-08-31T23:00 2>&1 | head -40
```

**`snow_storage.smet`** — generated from all 3,672 rows of DATA_2024.csv (Apr–Aug 2024), with fields:
- `TA` (K), `RH` (0–1), `VW` (m/s), `ISWR` (W/m²), `PSUM` (kg/m²/h from m/h×1000), `TSG` (soil temp in K from `Soil_Temp_320cm`)

**`snow_storage.sno`** — 3-layer initial profile matching your Python setup:
- Hs=4.5m (3×1.5m), ρ=400 kg/m³, ice fraction≈0.436, LWC=0
- Layer temperatures: T1=−2°C, T2=−4°C, T3=−6°C

**`snow_storage.ini`** — key parameter mappings:
| Python | SNOWPACK |
|---|---|
| `dt=600s` | `CALCULATION_STEP_LENGTH=10.0` |
| `h_conv=8.0, ε=0.95` | `ROUGHNESS_LENGTH=0.0025`, `EMISSIVITY_SNOW=0.95` |
| `Hi=0.20m, k_dry=0.06` | `CRUST_THICKNESS/CONDUCTIVITY` |
| `theta_e=0.055119`, bucket percolation | `RESIDUAL_WATER_CONTENT`, `LB_COND_WATERFLUX=BUCKET` |
| Robin BC ground, `h_ground=3.675`, TSG | `BOTTOM_BOUNDARY_CONDITION=SOIL_TEMP`, `k_soil=1.5` |

**One thing to verify:** the `latitude`/`longitude`/`altitude` in both `.smet` and `.sno` are placeholders (60.0°N, 25.0°E, 100m) — replace with your actual site coordinates so SNOWPACK's radiation geometry is correct.
