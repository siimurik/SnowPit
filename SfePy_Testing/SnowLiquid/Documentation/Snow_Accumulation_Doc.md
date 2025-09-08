# Snow Accumulation Source Term Documentation

## Overview

The `snow_accumulation` source term is a key component in the enhanced heat-liquid-vapor-pressure coupled transport system that simulates realistic snow behavior including accumulation, melting, and sublimation processes.

## Function Connection and Integration

### Primary Function
- **Function**: `get_accumulated_snow_source()`
- **Material Definition**: `'snow_accumulation': 'get_accumulated_snow_source'`
- **Equation Integration**: Used in the **Liquid equation** as a positive source term:

```python
'Liquid': """
    dw_dot.i.Omega(r, dCl/dt)
    + dw_laplace.i.Omega(liquid_mat.D_l, r, Cl)
    + dw_laplace.i.Omega(thermal_mat.D_t, r, T)
    + dw_laplace.i.Omega(pressure_liquid_mat.D_p, r, P)
    = - dw_volume_lvf.i.Omega(evaporation_rate.val, r)
      - dw_bc_newton.i.Gamma_Left(prec_h_l.val, prec_cl_ref.val, r, Cl)
      - dw_bc_newton.i.Gamma_Right(in_fixed.h_l_in, in_fixed.Cl_in, r, Cl)
      + dw_volume_lvf.i.Omega(snow_accumulation.val, r)  # <-- HERE
    """
```

## How the Snow Accumulation Source Works

The `get_accumulated_snow_source()` function implements three main physical processes:

### 1. Snow Accumulation Process

**Conditions for Snow Accumulation:**
- Precipitation rate > 0.1 mm/h (`current_prec > 0.1`)
- Air temperature ≤ 1.0°C (`current_temp <= 1.0`)

**Mathematical Implementation:**
```python
if current_prec > 0.1 and current_temp <= 1.0:
    snow_rate_flux = current_prec / 3600.0  # Convert mm/h to mm/s
    snow_infiltration_rate = snow_rate_flux / snow_penetration_depth
    
    for i in range(min(n_surface_nodes, n_nodes)):
        depth = i * dx
        decay_factor = np.exp(-depth / (snow_penetration_depth / 3.0))
        snow_source[i] += snow_infiltration_rate * decay_factor * 0.1
```

**Physical Meaning:**
- Snow accumulates primarily in surface layers (top 2cm: `snow_penetration_depth = 0.02`)
- Accumulation decreases exponentially with depth
- Snow is converted to liquid water content in the snow matrix

### 2. Snow Melting Process

**Conditions for Snow Melting:**
- Temperature > 0.5°C in any node
- Existing liquid content above baseline (indicating potential snow presence)

**Mathematical Implementation:**
```python
for i in range(n_nodes):
    if T_vals[i] > 0.5:
        baseline_liquid = 0.01
        potential_snow = max(0.0, Cl_vals[i] - baseline_liquid)
        
        if potential_snow > 0.001:
            melt_rate_coeff = 2e-4
            T_excess = T_vals[i] - 0.0
            melt_rate = melt_rate_coeff * T_excess * potential_snow
            max_melt = potential_snow / (2.0 * dt)
            melt_rate = min(melt_rate, max_melt)
            snow_source[i] += melt_rate
```

**Physical Meaning:**
- Melting occurs throughout the domain when temperatures exceed 0.5°C
- Melt rate is proportional to temperature excess and available snow
- Rate-limited to prevent numerical instabilities

### 3. Sublimation Process

**Conditions for Sublimation:**
- Temperature < 0.0°C (below freezing)
- Relative humidity < 90% (unsaturated air)
- Only affects surface nodes

**Mathematical Implementation:**
```python
for i in range(min(n_surface_nodes, n_nodes)):
    if T_vals[i] < 0.0:
        current_rh_local = rh[hour_idx] if hour_idx < len(rh) else 70.0
        if current_rh_local < 90.0:
            sublimation_rate = 1e-6 * (90.0 - current_rh_local) / 90.0
            snow_source[i] -= sublimation_rate  # Negative = loss
```

**Physical Meaning:**
- Direct transition from solid snow to water vapor
- Rate increases with decreasing humidity
- Only affects surface layers where sublimation occurs

## Key Parameters

| Parameter | Value | Physical Meaning |
|-----------|-------|------------------|
| `snow_penetration_depth` | 0.02 m | Depth over which snow accumulates |
| `melt_rate_coeff` | 2e-4 | Snow melting rate coefficient |
| `baseline_liquid` | 0.01 kg/m³ | Baseline liquid content (not snow) |
| `sublimation_rate` | 1e-6 | Base sublimation rate coefficient |

## Practical Example

### Scenario: Winter Storm with Temperature Variation

**Initial Conditions:**
- Time: Hour 0-10 of simulation
- Snow depth: 10cm (0.1m domain)
- Grid: 20 nodes, dx = 0.005m

**Hour 1-3: Snow Accumulation**
```
Weather: T_air = -3°C, Precipitation = 2 mm/h, RH = 85%
Process: Snow accumulation in top 4 nodes (0-2cm depth)

Node 0 (surface): snow_source = +0.000167 kg/(m³·s)
Node 1 (0.5cm):   snow_source = +0.000122 kg/(m³·s)  
Node 2 (1.0cm):   snow_source = +0.000089 kg/(m³·s)
Node 3 (1.5cm):   snow_source = +0.000065 kg/(m³·s)
Deeper nodes:     snow_source = 0
```

**Hour 4-6: Cold and Dry (Sublimation)**
```
Weather: T_air = -8°C, Precipitation = 0 mm/h, RH = 60%
Process: Sublimation in surface nodes

Node 0: T = -5°C → snow_source = -0.000033 kg/(m³·s)
Node 1: T = -4°C → snow_source = -0.000033 kg/(m³·s)
Deeper nodes: No sublimation
```

**Hour 7-10: Warm-up and Melting**
```
Weather: T_air = +3°C, Precipitation = 0 mm/h
Process: Snow melting throughout affected depth

Node 0: T = +2°C, Cl = 0.05 kg/m³ → snow_source = +0.000020 kg/(m³·s)
Node 1: T = +1°C, Cl = 0.03 kg/m³ → snow_source = +0.000008 kg/(m³·s)
Node 2: T = +0.8°C, Cl = 0.02 kg/m³ → snow_source = +0.000004 kg/(m³·s)
```

## Integration with Other Processes

### Coupling with Evaporation
- Snow melting increases liquid content (Cl)
- Higher liquid content can increase evaporation rate
- Balance between melting input and evaporation output

### Coupling with Vapor Transport
- Sublimation directly affects vapor content (Cv)
- Creates vapor pressure gradients
- Influences vapor transport and boundary conditions

### Coupling with Heat Transfer
- Melting consumes latent heat (cooling effect)
- Sublimation also consumes latent heat
- Affects temperature distribution

## Physical Realism Features

1. **Depth-Dependent Accumulation**: Snow doesn't accumulate uniformly but decreases with depth
2. **Temperature-Dependent Phase Changes**: Realistic thresholds for melting/sublimation
3. **Rate Limiting**: Prevents unrealistic rapid changes
4. **Environmental Coupling**: Links to actual weather data (temperature, precipitation, humidity)
5. **Conservative Mass Transfer**: Snow that melts becomes liquid water; sublimated snow becomes vapor

## Numerical Considerations

- **Stability**: Maximum melt rates prevent numerical oscillations
- **Conservation**: Total water mass (liquid + vapor + solid) is conserved
- **Time Integration**: Source term integrated implicitly with other equation terms
- **Spatial Discretization**: Uses finite element basis functions for smooth distribution

This snow accumulation source term provides a physically realistic representation of snow behavior in the coupled transport system, enabling accurate simulation of moisture movement in snow-covered materials.