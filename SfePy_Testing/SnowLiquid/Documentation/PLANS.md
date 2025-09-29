### 31.4 MOISTURE MOVEMENT AND DRYING BEHAVIOR

Several researches have been carried out to understand the mechanism of moisture movement in clay during drying. Newitt et al. [12] and Wakabayashi [13] investigated the moisture movement in clay by liquid diffusion and vapor diffusion, which affect the drying characteristics particularly the falling rate. They concluded that the liquid diffusion dominates the movement until about $20\%$-dry basis in moisture content for stoneware clay and $30\%$ for the mixture of $80\%$ Kibushi clay and $20\%$ feldspar. Wakabayashi [14] also evaluated the effective moisture diffusion coefficient of some sorts of clay such as Kibushi, Gairome, stoneware, feldspar, and their mixtures. The effective diffusion coefficient is available for the brief description of the moisture movement behavior. The effective diffusion coefficient D can be defined by

$$
N = -\frac{D}{v_0} \frac{\mathrm{d} w}{\mathrm{d}x}
$$

They correlated $D$ as

$$
D = \frac{\alpha K \gamma^2 \rho_{\mathrm{p}} w^3 e^{-Kw}}{k \mu S^2 (1+\gamma w)}
$$

where $k$ is Carman’s constant ( $= 5$), $N$ is the mass flux of moisture, $S$ is the specific surface area of clay, $v_0$ is the specific volume of dry clay, $w$ is the moisture content with dry base, $g$ is the specific gravity of a clay particle, $\rho_{\mathrm{p}}$ is the density of a clay particle, and $m$ is the viscosity of liquid water. Constants $\alpha$ and $K$ are given empirically by $2.02 \cdot 10^{11} \; \mathrm{Pa}$ and $24.9$ kg-solid/kg-water, respectively. 

To understand the mechanism further, a microscopic investigation on the moisture content in clay is required. The osmotic suction potential was introduced as the driving force of moisture movement, described in Section 31.3, and was successively applied to the prediction of moisture movement in wet clay [15,16]. The theoretical analysis on the twodimensional moisture transfer of cylindrical clay was performed taking into account the effects of both osmotic suction and strain–stress caused by the shrinkage [17,18]. However, only the transient mass-transfer equation was analyzed, assuming a constant drying rate on the external surface of the clay. Comini and Lewis [19] developed the finite element method to solve simultaneous heat and moisture transfer for noncompressible porous media for a complicated geometry of the axial symmetry such as electric insulators. The three-dimensional problem on heat and moisture transfer, involving the drying shrinkage, was analyzed by Hasatani et al. [20]. Their model was not only limited to surface evaporation period but also applied to geometries more complicated than a simple slab shape. In the future, it is expected that drying kinetics would be experimentally studied and the simple analysis be made available for the entire drying periods. 

Drying is a macroscopical phenomenon involving simultaneous heat and mass transfer. Suppose that heat conduction and moisture diffusion are dominant during the overall transfer process in a homogeneous medium. The multidimensional conservation equations for an anisotropic medium can then be expressed as

$$
\begin{gathered}
c_{\mathrm{p}} \rho_{\mathrm{m}} \frac{\partial T}{\partial t}=-\nabla \mathbf{J}_{\mathbf{h}}+\dot{q}+\Delta H_{\mathrm{V}} \varepsilon_{\mathrm{L}} \frac{\partial C}{\partial t} \\
\frac{\partial C}{\partial t}=-\nabla \mathbf{J}_{\mathbf{m}}
\end{gathered}
$$

The heat flux $\mathbf{J_h}$ and mass flux $\mathbf{J_m}$ in Equation 31.3 and Equation 31.4 are given by

$$
\begin{array}{r}
\mathbf{J}_{\mathbf{h}}=-k_{\mathrm{t}} \nabla T-k_{\mathrm{c}} \nabla C-k_{\mathrm{p}} \nabla p \\
\mathbf{J}_{\mathbf{m}}=-D_{\mathrm{W}} \nabla C-D_{\mathrm{t}} \nabla T-D_{\mathrm{p}} \nabla p
\end{array}
$$

The motion of water is influenced by the water pressure, which is induced by the capillary or osmotic suction resulting from evaporation in porous media. The water flows in the relatively larger pores saturated with water rather than by diffusion process. In such a case, the pressure diffusion term dominates the mass-transfer rate in Equation 31.6. The Soret effect described by the temperature diffusion term is generally considered small compared with other terms.
Then Equation 31.6 can be simplified as

$$
\mathbf{J_m} = -D_{\mathbf{p}} \nabla p
$$

The liquid flow rate in the porous media is given by Darcy’s law:

$$
\mathbf{J}_w = - \frac{k_s}{\mu} \nabla p
$$

where $k_s$ denotes the permeability. As the mass-transfer rate should equal the water-flow rate if the vapor flow is ignored, the pressure diffusivity $D_{\mathbf{p}}$ is derived from both Equation 31.8 and Equation 31.7:

$$
D_{\mathrm{p}} = \frac{k_s}{\mu}
$$

Hence Equation 31.4 can be rewritten from Equation 31.7 and Equation 31.9, assuming uniform properties in the body,

$$
\frac{\partial C}{\partial t} = \frac{k_s}{\mu} \nabla^2 p
$$

The capillary pressure of water is determined by a balance of the interfacial energies among the three phases: solid, liquid, and vapor, the wetting angle of liquid–solid, and the radii of the pores, which are dependent on the pore structure in the medium and the amount of water existing in the pores. The osmotic pressure is dependent on the pore structure and the liquid–solid interface. Therefore, one must predict the pressure $p$ statistically from the pore structure distribution model. In order to simplify the model, alternatively, the mass-transfer equation is often expressed by introducing the moisture content w as the driving force for moisture transfer,

$$
\frac{\partial w}{\partial t} = \nabla \cdot \left(  \frac{D}{v_0} \nabla w \right)
$$

The initial and boundary conditions must be specified
depending on the drying system and the surrounding
atmosphere to which a medium is exposed.


--- 

Based on the reference material from the *Handbook of Industrial Drying*, we can construct a **system of coupled differential equations** that account for **heat transfer, liquid moisture transport, and vapor diffusion** in porous media (e.g., clay, snow, or insulation materials). Below is the proposed system, incorporating key mechanisms from the chapter:

---

### **Governing Equations for Coupled Heat and Moisture Transport**

#### **1. Heat Transfer Equation (Energy Conservation)**  
Incorporates conduction, latent heat effects, and convective coupling with moisture:  
$$
\rho_m c_p \frac{\partial T}{\partial t} = \nabla \cdot \left( k_t \nabla T \right) + \Delta H_v \varepsilon_L \frac{\partial C}{\partial t} + \dot{q}
$$  
**Terms**:  
- $ \rho_m $: Bulk density of wet material (solid + moisture)  
- $ c_p $: Specific heat capacity  
- $ k_t $: Thermal conductivity  
- $ \Delta H_v $: Latent heat of vaporization  
- $ \varepsilon_L $: Liquid-phase volume fraction  
- $ \dot{q} $: External heat source (e.g., solar radiation)  

---

#### **2. Mass Transfer Equation (Moisture Conservation)**  
Combines **liquid diffusion**, **vapor diffusion**, and **capillary-driven flow**:  
$$
\frac{\partial C}{\partial t} = \nabla \cdot \left( D_w \nabla C + D_t \nabla T + \frac{k_s}{\mu} \nabla p \right)
$$  
**Simplified Form (Dominant Capillary Flow)**:  
$$
\frac{\partial C}{\partial t} = \frac{k_s}{\mu} \nabla^2 p
$$  
**Terms**:  
- $ C $: Moisture concentration (kg/m³)  
- $ D_w $: Moisture diffusivity (liquid phase)  
- $ D_t $: Thermal moisture diffusivity (Soret effect)  
- $ k_s $: Permeability of porous medium  
- $ \mu $: Dynamic viscosity of water  
- $ p $: Capillary/osmotic pressure (driving force)  

---

I want to solve this system of partial differential equations:

$$
\begin{aligned}
&\text{(Heat)} \quad \rho_m c_p \frac{\partial T}{\partial t} = \nabla \cdot (k_t \nabla T) + L_v \left( \frac{\partial C_l}{\partial t} + \nabla \cdot \mathbf{J_v} \right) + \dot{q} \\
&\text{(Liquid)} \quad \frac{\partial C_l}{\partial t} = \nabla \cdot \left( D_l \nabla C_l + \frac{k_s}{\mu} (\nabla p_c + \rho_l g \nabla z) \right) - \dot{m}_{evap} \\
&\text{(Vapor)} \quad \frac{\partial C_v}{\partial t} = \nabla \cdot (D_v \nabla \rho_v + \rho_v \mathbf{vg}) + \dot{m}_{evap} \\
&\text{(Air)} \quad \phi \frac{\partial \rho_a}{\partial t} = \nabla \cdot \left( \rho_a \frac{k_g}{\mu_g} \nabla p_a \right)
\end{aligned}
$$

where:
- $ \dot{m}_{evap} $ = evaporation/condensation rate,
- $ \mathbf{v_g} = -\frac{k_g}{\mu_g} \nabla p_a $ = gas velocity,
- $ \phi $ = porosity.
The vapor flux is given by:

$$
\mathbf{J_v} = \underbrace{-D_v \nabla \rho_v}_{\text{Vapor diffusion}} + \underbrace{\rho_v \mathbf{vg}}_{\text{Vapor advection}}
$$

currently I have solved this part SfePy
$$
\left\lbrace
\begin{aligned}
&\frac{\partial T}{\partial t}(x, t) = \alpha \frac{\partial^2 T}{\partial x^2}(x, t) && \text{(Heat equation)} \\
&\text{Initial condition:} \\
&T(x, 0) = T_0(x) && \forall x \in [0, L] \\
&\text{Boundary conditions:} \\
&-k \frac{\partial T}{\partial x}(0, t) = h_o \big(T(0, t) - T_o(t)\big) && \text{(Robin BC at $x=0$)} \\
&-k \frac{\partial T}{\partial x}(L, t) = h_i \big(T(L, t) - T_i(t)\big) && \text{(Robin BC at $x=L$)} \\
&\text{where:} \\
&\alpha = \frac{k}{\rho c_p} \quad \text{(Thermal diffusivity)}
\end{aligned}
\right.
$$

It is named test_fem_snow.py

I need to start adding more terms. what can I start adding step by step to make the project manageable. I have some examples which I plan to use but do not know how to implement correctly. Here is one example where the use coupled equations named laplace_coupling_lcbcs.py which could be used as a reference to expand the equation.

---

Looking at your current heat equation implementation and the target system, I can suggest a step-by-step approach to gradually build up to the full coupled system. Let me analyze what you have and propose a manageable progression.

## Current State Analysis

Your `test_fem_snow.py` successfully implements a 1D transient heat equation with Robin boundary conditions. This is a solid foundation.

## Suggested Step-by-Step Implementation Plan

### Step 1: Add Source Terms to Heat Equation
Start by adding the latent heat term and heat source to your existing heat equation:

```python
# Modified heat equation with source terms
equations = {
    'Heat': """dw_dot.i.Omega(mat.rho_cp, s, dT/dt)
             + dw_laplace.i.Omega(mat.lam, s, T)
             = - dw_bc_newton.i.Gamma_Left(h_out_dyn.val, T_out_dyn.val, s, T)
               - dw_bc_newton.i.Gamma_Right(in_fixed.h_in, in_fixed.T_in, s, T)
               + dw_volume_lvf.i.Omega(mat.L_v, s, dCl_dt)
               + dw_volume_lvf.i.Omega(mat.q_dot, s, 1)
             """
}
```

This requires adding:
- `L_v` (latent heat of vaporization) to materials
- `dCl_dt` (time derivative of liquid content) - start with a simple prescribed function
- `q_dot` (heat source term) - can be zero initially

### Step 2: Add Liquid Water Transport (Simplified)
Add the liquid water equation as a second field, but start with diffusion only:

```python
fields = {
    'temperature': ('real', 1, 'Omega', 1),
    'liquid': ('real', 1, 'Omega', 1),  # New field for liquid content
}

variables = {
    'T': ('unknown field', 'temperature', 0, 1),
    's': ('test field', 'temperature', 'T'),
    'Cl': ('unknown field', 'liquid', 1, 1),    # New liquid variable
    'r': ('test field', 'liquid', 'Cl'),        # New test function
}

equations = {
    'Heat': """dw_dot.i.Omega(mat.rho_cp, s, dT/dt)
             + dw_laplace.i.Omega(mat.lam, s, T)
             = - dw_bc_newton.i.Gamma_Left(h_out_dyn.val, T_out_dyn.val, s, T)
               - dw_bc_newton.i.Gamma_Right(in_fixed.h_in, in_fixed.T_in, s, T)
               + dw_volume_lvf.i.Omega(mat.L_v, s, dCl/dt)""",
               
    'Liquid': """dw_dot.i.Omega(mat.one, r, dCl/dt)
               + dw_laplace.i.Omega(mat.D_l, r, Cl)
               = dw_volume_lvf.i.Omega(mat.m_evap, r, 1)"""
}
```

### Step 3: Add Coupling Between Heat and Liquid
Introduce temperature-dependent evaporation:

```python
def get_evap_rate(ts, coors, mode=None, **kwargs):
    """Temperature-dependent evaporation rate"""
    if mode != 'qp':
        return {}
    
    # Simple temperature-dependent model
    # In reality, this would depend on T, Cl, and vapor pressure
    val = nm.zeros((coors.shape[0], 1, 1), dtype=nm.float64)
    # Add your evaporation model here
    
    return {'val': val}
```

### Step 4: Add Vapor Transport (Simplified)
Add vapor as a third field with diffusion only initially:

```python
fields = {
    'temperature': ('real', 1, 'Omega', 1),
    'liquid': ('real', 1, 'Omega', 1),
    'vapor': ('real', 1, 'Omega', 1),  # New vapor field
}

equations = {
    # ... existing equations ...
    'Vapor': """dw_dot.i.Omega(mat.one, w, dCv/dt)
              + dw_laplace.i.Omega(mat.D_v, w, Cv)
              = - dw_volume_lvf.i.Omega(mat.m_evap, w, 1)"""
}
```

### Step 5: Add Air Transport
Finally, add the air equation (this can be simpler in 1D):

### Step 6: Add Advection Terms
Once you have the basic diffusion working, add the advection terms using `dw_convect` terms.

## Immediate Next Step Implementation

Here's a concrete first step you can implement right now:## Key Changes in Step 1:

1. **Added latent heat term**: `dw_volume_lvf.i.Omega(mat.L_v, liquid_rate.val, s)` represents the latent heat release/absorption due to phase change
2. **Added heat source term**: `dw_volume_lvf.i.Omega(heat_source.val, s)` for any additional heat sources
3. **Prescribed liquid change rate**: `get_liquid_change_rate()` provides a simple test function for dCl/dt

## Testing This Step:

1. Create a simple CSV file `t_o_and_h_o.csv` with temperature and heat transfer coefficient data
2. Run this version and observe how the latent heat term affects the temperature evolution
3. Compare results with your original version to see the impact

## Next Steps After This Works:

Once you verify this step works correctly, you can move to Step 2 by:

1. Adding the liquid field as a second unknown
2. Implementing the liquid transport equation 
3. Coupling the evaporation rate between heat and liquid equations

This approach lets you:
- **Test incrementally** - each step builds on the previous working version
- **Debug easily** - isolate problems to the newly added components
- **Understand physics** - see how each term affects the solution
- **Maintain working code** - always have a functioning baseline

Would you like me to help you implement Step 2 (adding the liquid field) once you've tested this first enhancement?

### 1. Enhance Evaporation/Condensation Modeling (\dot{m}_{evap})
The current implementation uses a simplistic empirical rate: $\dot{m}_{evap} = k_{evap} \cdot C_l \cdot \max(0, T - T_{ref})$. This assumes **evaporation only occurs above a reference temperature and scales linearly with liquid content, without considering vapor dynamics or thermodynamic equilibrium**.

- **Improvements for Realism**:
  - Adopt a physics-based model for $\dot{m}_{evap}$, such as one based on the Hertz-Knudsen equation or vapor pressure deficit: $\dot{m}_{evap} = \beta (p_{v,sat}(T) - p_v)$, where $\beta$ is a mass transfer coefficient, $p_{v,sat}$ is saturation vapor pressure (e.g., using Clausius-Clapeyron relation), and $p_v$ is local vapor pressure derived from $\rho_v$ and ideal gas law.
  - Incorporate porosity ($\phi$) and specific surface area (e.g., for snow grains, ~210 m²/m³ as noted in comments) to scale evaporation with available interface area: $\dot{m}_{evap} \propto a_{int} (p_{v,sat} - p_v)$, where $a_{int}$ is the ice-liquid-vapor interface area.
  - Add condensation when p_v > p_{v,sat} (negative \dot{m}_{evap}), enabling refreezing or frost formation in cold/dry conditions.
  - Temperature dependence: Make k_evap (or \beta) vary with T, as diffusion and viscosity change (e.g., Arrhenius form: k_evap \propto e^{-E_a / (R T)}).

- **Implementation Tips**:
  - In `get_evaporation_rate`, compute p_{v,sat} using an empirical formula like p_{v,sat} = 611.2 \exp(17.67 T / (T + 243.5)) (Tetens formula, T in °C).
  - Require solving for \rho_v from the vapor equation (see below), coupling it to evaporation.
  - Validate against snowpack data (e.g., from SNOWPACK model), where evaporation rates peak during melt periods.

### 2. Incorporate Vapor Transport and Advection (\mathbf{J_v} and Vapor Equation)
The code omits the vapor equation and flux \mathbf{J_v} = -D_v \nabla \rho_v + \rho_v \mathbf{v_g}, simplifying latent heat to only -L_v \cdot \dot{m}_{evap} instead of L_v (\partial C_l / \partial t + \nabla \cdot \mathbf{J_v}).

- **Improvements for Realism**:
  - Add the vapor conservation equation: \partial C_v / \partial t = \nabla \cdot (D_v \nabla \rho_v + \rho_v \mathbf{v_g}) + \dot{m}_{evap}, where C_v = \phi (1 - S_l) \rho_v (S_l = liquid saturation = C_l / (\phi \rho_l)).
  - Explicitly compute and diverge \mathbf{J_v} in the heat equation for vapor advection effects, which are crucial in porous media like snow where vapor diffusion drives sublimation/deposition.
  - Include gas velocity \mathbf{v_g} = - (k_g / \mu_g) \nabla p_a from Darcy's law, coupling to the air equation.
  - Account for porosity (\phi): Scale vapor diffusion D_v = D_{v0} \phi^{4/3} (1 - S_l)^{10/3} (tortuosity correction, e.g., Millington-Quirk model).

- **Implementation Tips**:
  - Introduce a new field/variable for \rho_v in SfePy (similar to 'liquid' for C_l).
  - Add terms like `dw_div_grad.i.Omega(D_v, s, rho_v)` for diffusion and an advection term using `dw_convect.i.Omega(v_g, s, rho_v)`.
  - Boundary conditions: Set vapor BC based on ambient humidity, e.g., \rho_{v,bc} = RH \cdot \rho_{v,sat}(T_o) (RH from data).
  - This will capture realistic vapor-driven heat pipes in snow, where sublimation at warm bases and deposition at cold surfaces occur.

### 3. Add Capillary and Gravity Effects in Liquid Transport
Liquid transport is limited to diffusion (\nabla \cdot (D_l \nabla C_l)), ignoring capillary pressure and gravity: (k_s / \mu) (\nabla p_c + \rho_l g \nabla z).

- **Improvements for Realism**:
  - Implement Richards equation form: \partial C_l / \partial t = \nabla \cdot [ (k_s / \mu) (\nabla p_c + \rho_l g \nabla z) ] - \dot{m}_{evap}, where p_c = p_c(S_l) (e.g., van Genuchten model: p_c = -1/\alpha [S_l^{-1/m} - 1]^{1/n}).
  - Parameterize for snow: Low permeability k_s (~10^{-10} m²), gravity dominant in wet snow (downward drainage), capillary for retention in dry snow.
  - Make D_l moisture-dependent: D_l(S_l) = D_{l0} S_l^\gamma (e.g., \gamma=3 for unsaturated flow).

- **Implementation Tips**:
  - Add gravity as a source term: `dw_volume_lvf.i.Omega( (k_s / \mu) \rho_l g , r)` (for vertical 1D, z increasing upward).
  - For capillary: Introduce nonlinear diffusion or use a pressure-based formulation (solve for p_c instead of C_l, with C_l = C_l(p_c)).
  - Use snow-specific retention curves (e.g., from Yamaguchi et al., 2010: residual saturation ~0.07, van Genuchten parameters \alpha=3.5 m^{-1}, n=4.3).

### 4. Include Air Flow and Pressure Dynamics (Air Equation)
The air equation is absent, so gas velocity \mathbf{v_g} isn't computed, neglecting advection in vapor and potential pressure-driven flows.

- **Improvements for Realism**:
  - Add \phi \partial \rho_a / \partial t = \nabla \cdot (\rho_a (k_g / \mu_g) \nabla p_a), with \mathbf{v_g} feeding into vapor advection.
  - Couple to evaporation: \dot{m}_{evap} produces vapor, affecting p_a if assuming ideal gas mixture (p_a = p_g - p_v).
  - For snow, include wind-pumping effects at boundaries: Time-varying p_a gradients from ambient wind.

- **Implementation Tips**:
  - New field for p_a or \rho_a.
  - Compute \mathbf{v_g} in a material function and use in advection terms.
  - Boundary: Set p_a = atmospheric pressure + wind-induced fluctuations (e.g., sinusoidal based on data).

### 5. Refine Boundary Conditions and Environmental Coupling
Outer BCs use Robin type with humidity-to-liquid conversion, but the conversion is linear and ad-hoc, ignoring saturation vapor pressure.

- **Improvements for Realism**:
  - Enhance `humidity_to_liquid_content`: Use sorption isotherms for snow (e.g., hysteresis in wetting/drying), incorporating p_{v,sat}(T): C_l,eq = C_{max} [ RH \exp( -p_c / (\rho_l R_v T) ) ]^{1/\lambda}.
  - Add precipitation: In `get_cl_o`, include rain/snowmelt influx as a flux BC: J_l,in = v_rain \rho_l + q_rain (sensible heat from rain).
  - Inner BC: If modeling building insulation, couple to indoor vapor diffusion.

- **Implementation Tips**:
  - Integrate Psat_WV as noted in comments: p_{sat} = 610.78 \exp(17.27 T / (T + 237.3)).
  - Add solar radiation: q_dot = (1 - albedo) SWR + LWR - emitted (Stefan-Boltzmann), using data for t_o.

### 6. Account for Material Property Variations and Phase Changes
Properties like k_t, D_l, c_p are constant; no explicit ice phase.

- **Improvements for Realism**:
  - Temperature/moisture dependence: k_t = k_dry + b C_l (b~2.5 for snow), \rho c_p varying with ice/liquid fractions.
  - Explicit freezing: If T < 0, convert liquid to ice with latent heat release, reducing permeability.
  - Porosity evolution: \phi decreases with densification from melt/refreeze.

- **Implementation Tips**:
  - Make materials state-dependent in SfePy (e.g., via functions).
  - Track ice fraction as another variable.

### 7. Numerical and Validation Enhancements
- **Nonlinear Coupling**: Current Newton solver may need tuning (e.g., tighter eps, line search) for strong evap coupling.
- **Dimensionality**: Extend to 2D/3D for lateral flows in snowpacks.
- **Validation**: Compare to field data (e.g., SNOTEL sites) or models like CROCUS/SNOWPACK.
- **Efficiency**: Use adaptive time-stepping for fast events (e.g., rain onset).

```python
Liquid --[evaporation]--> Vapor --[diffusion]--> Boundaries
   ↓                        ↑
[Heat removal]         [Heat added during condensation]
```

To convert the given strong form equations to weak integral form and express them in SfePy syntax, we'll break down each term. The system consists of two equations: one for heat transfer and one for mass transfer.

### Given Strong Form Equations:
1. Heat equation:
   $$
   c_{\mathrm{p}} \rho_{\mathrm{m}} \frac{\partial T}{\partial t} = -\nabla \cdot \mathbf{J}_{\mathbf{h}} + \dot{q} + \Delta H_{\mathrm{V}} \varepsilon_{\mathrm{L}} \frac{\partial C}{\partial t}
   $$
2. Mass equation:
   $$
   \frac{\partial C}{\partial t} = -\nabla \cdot \mathbf{J}_{\mathbf{m}}
   $$

With fluxes:
   $$
   \mathbf{J}_{\mathbf{h}} = -k_{\mathrm{t}} \nabla T - k_{\mathrm{c}} \nabla C - k_{\mathrm{p}} \nabla p
   $$
   $$
   \mathbf{J}_{\mathbf{m}} = -D_{\mathrm{W}} \nabla C - D_{\mathrm{t}} \nabla T - D_{\mathrm{p}} \nabla p
   $$

### Step 1: Write the weak form for each equation
We multiply each equation by a test function and integrate over the domain $\Omega$. For the heat equation, use test function $s_T$; for the mass equation, use test function $s_C$.

#### Heat equation weak form:
Multiply by $s_T$ and integrate:
$$
\int_{\Omega} s_T \left( c_p \rho_m \frac{\partial T}{\partial t} \right) = \int_{\Omega} s_T \left( -\nabla \cdot \mathbf{J_h} + \dot{q} + \Delta H_V \varepsilon_L \frac{\partial C}{\partial t} \right)
$$

Apply integration by parts to the flux divergence term:
$$
\int_{\Omega} s_T (-\nabla \cdot \mathbf{J_h}) = \int_{\Omega} \nabla s_T \cdot \mathbf{J_h} - \int_{\partial \Omega} s_T (\mathbf{J_h} \cdot \mathbf{n})
$$
Typically, the boundary term is handled by Neumann conditions. So:
$$
\int_{\Omega} s_T c_p \rho_m \frac{\partial T}{\partial t} = \int_{\Omega} \nabla s_T \cdot \mathbf{J_h} - \int_{\partial \Omega} s_T (\mathbf{J_h} \cdot \mathbf{n}) + \int_{\Omega} s_T \dot{q} + \int_{\Omega} s_T \Delta H_V \varepsilon_L \frac{\partial C}{\partial t}
$$

Rearrange:
$$
\int_{\Omega} s_T c_p \rho_m \frac{\partial T}{\partial t} - \int_{\Omega} \nabla s_T \cdot \mathbf{J_h} - \int_{\Omega} s_T \Delta H_V \varepsilon_L \frac{\partial C}{\partial t} = \int_{\Omega} s_T \dot{q} - \int_{\partial \Omega} s_T (\mathbf{J_h} \cdot \mathbf{n})
$$

#### Mass equation weak form:
Multiply by $s_C$ and integrate:
$$
\int_{\Omega} s_C \frac{\partial C}{\partial t} = \int_{\Omega} s_C (-\nabla \cdot \mathbf{J_m})
$$

Apply integration by parts:
$$
\int_{\Omega} s_C \frac{\partial C}{\partial t} = \int_{\Omega} \nabla s_C \cdot \mathbf{J_m} - \int_{\partial \Omega} s_C (\mathbf{J_m} \cdot \mathbf{n})
$$

So:
$$
\int_{\Omega} s_C \frac{\partial C}{\partial t} - \int_{\Omega} \nabla s_C \cdot \mathbf{J_m} = - \int_{\partial \Omega} s_C (\mathbf{J_m} \cdot \mathbf{n})
$$

### Step 2: Substitute the flux expressions
Substitute $\mathbf{J_h}$ and $\mathbf{J_m}$:

$$
\mathbf{J_h} = -k_t \nabla T - k_c \nabla C - k_p \nabla p
$$
$$
\mathbf{J_m} = -D_W \nabla C - D_t \nabla T - D_p \nabla p
$$

So the weak forms become:

#### Heat equation:
$$
\int_{\Omega} s_T c_p \rho_m \frac{\partial T}{\partial t} 
- \int_{\Omega} \nabla s_T \cdot \left( -k_t \nabla T - k_c \nabla C - k_p \nabla p \right)
- \int_{\Omega} s_T \Delta H_V \varepsilon_L \frac{\partial C}{\partial t}
= \int_{\Omega} s_T \dot{q} - \int_{\partial \Omega} s_T (\mathbf{J_h} \cdot \mathbf{n})
$$

Which expands to:
$$
\int_{\Omega} s_T c_p \rho_m \frac{\partial T}{\partial t}
+ \int_{\Omega} k_t \nabla s_T \cdot \nabla T
+ \int_{\Omega} k_c \nabla s_T \cdot \nabla C
+ \int_{\Omega} k_p \nabla s_T \cdot \nabla p
- \int_{\Omega} s_T \Delta H_V \varepsilon_L \frac{\partial C}{\partial t}
= \int_{\Omega} s_T \dot{q} - \int_{\partial \Omega} s_T (\mathbf{J_h} \cdot \mathbf{n})
$$

#### Mass equation:
$$
\int_{\Omega} s_C \frac{\partial C}{\partial t}
- \int_{\Omega} \nabla s_C \cdot \left( -D_W \nabla C - D_t \nabla T - D_p \nabla p \right)
= - \int_{\partial \Omega} s_C (\mathbf{J_m} \cdot \mathbf{n})
$$

Which expands to:
$$
\int_{\Omega} s_C \frac{\partial C}{\partial t}
+ \int_{\Omega} D_W \nabla s_C \cdot \nabla C
+ \int_{\Omega} D_t \nabla s_C \cdot \nabla T
+ \int_{\Omega} D_p \nabla s_C \cdot \nabla p
= - \int_{\partial \Omega} s_C (\mathbf{J_m} \cdot \mathbf{n})
$$

### Step 3: Write in SfePy syntax
In SfePy, weak forms are expressed using terms like `dw_dot` (for time derivatives), `dw_laplace` (for diffusion), `dw_diffusion` (for general diffusion with coefficient matrix), and `dw_volume_integrate` (for source terms). Boundary fluxes are handled with `dw_surface_integrate`.

Assume:
- $T$ and $C$ are unknown fields.
- $p$ might be given or another unknown. For simplicity, we assume $p$ is given (or solved elsewhere).

#### Heat equation weak form in SfePy syntax:
Terms:
1. $\int_{\Omega} s_T c_p \rho_m \frac{\partial T}{\partial t}$ → `dw_dot.i.Omega(heat_capacity, s_T, dT/dt)`
   - `heat_capacity = c_p * rho_m`
2. $\int_{\Omega} k_t \nabla s_T \cdot \nabla T$ → `dw_laplace.i.Omega(thermal_conductivity, s_T, T)`
   - `thermal_conductivity = k_t`
3. $\int_{\Omega} k_c \nabla s_T \cdot \nabla C$ → `dw_diffusion.i.Omega(coupling_kc, s_T, C)`
   - `coupling_kc = k_c` (scalar) or matrix if anisotropic.
4. $\int_{\Omega} k_p \nabla s_T \cdot \nabla p$ → `dw_diffusion.i.Omega(kp_coeff, s_T, p)`
   - This term requires $p$ as a given function.
5. $- \int_{\Omega} s_T \Delta H_V \varepsilon_L \frac{\partial C}{\partial t}$ → `- dw_dot.i.Omega(enthalpy_coeff, s_T, dC/dt)`
   - `enthalpy_coeff = Delta_HV * epsilon_L`
6. $\int_{\Omega} s_T \dot{q}$ → `dw_volume_integrate.i.Omega(heat_source, s_T)`
   - `heat_source = dot_q`
7. $- \int_{\partial \Omega} s_T (\mathbf{J_h} \cdot \mathbf{n})$ → `- dw_surface_integrate.i.Gamma(heat_flux, s_T)`
   - This is typically used for Neumann BCs. If the flux is given, use this.

#### Mass equation weak form in SfePy syntax:
1. $\int_{\Omega} s_C \frac{\partial C}{\partial t}$ → `dw_dot.i.Omega(one, s_C, dC/dt)`
   - Here `one` is a coefficient equal to 1.
2. $\int_{\Omega} D_W \nabla s_C \cdot \nabla C$ → `dw_laplace.i.Omega(DW, s_C, C)`
   - `DW = D_W`
3. $\int_{\Omega} D_t \nabla s_C \cdot \nabla T$ → `dw_diffusion.i.Omega(Dt, s_C, T)`
   - `Dt = D_t`
4. $\int_{\Omega} D_p \nabla s_C \cdot \nabla p$ → `dw_diffusion.i.Omega(Dp, s_C, p)`
5. $- \int_{\partial \Omega} s_C (\mathbf{J_m} \cdot \mathbf{n})$ → `- dw_surface_integrate.i.Gamma(mass_flux, s_C)`

### Final Weak Form in SfePy Syntax:
For the heat equation:
```
dw_dot.i.Omega(heat_capacity, s_T, dT/dt)
+ dw_laplace.i.Omega(thermal_conductivity, s_T, T)
+ dw_diffusion.i.Omega(kc, s_T, C)
+ dw_diffusion.i.Omega(kp, s_T, p)
- dw_dot.i.Omega(enthalpy_coeff, s_T, dC/dt)
= dw_volume_integrate.i.Omega(heat_source, s_T)
- dw_surface_integrate.i.Gamma(heat_flux, s_T)
```

For the mass equation:
```
dw_dot.i.Omega(one, s_C, dC/dt)
+ dw_laplace.i.Omega(DW, s_C, C)
+ dw_diffusion.i.Omega(Dt, s_C, T)
+ dw_diffusion.i.Omega(Dp, s_C, p)
= - dw_surface_integrate.i.Gamma(mass_flux, s_C)
```

### Notes:
- The terms involving $\nabla p$ assume $p$ is known. If $p$ is also unknown, you need an additional equation.
- The boundary terms (`dw_surface_integrate`) are for Neumann conditions. If there are Dirichlet conditions, they are enforced directly on the field.
- The coefficients (like `heat_capacity`, `thermal_conductivity`, etc.) must be defined in the material.
- The `i.Omega` and `i.Gamma` refer to integration over the domain and boundary, respectively. The exact region name (e.g., `Omega`, `Gamma_N`) depends on your mesh regions.

This gives you the weak form and the corresponding SfePy syntax for your coupled heat and mass transfer system.

# Enhanced Heat-Liquid-Vapor Coupled Transport System
## Corrected Governing Equations

### Physical Model
Three-field coupled system based on clay drying research with enhanced physics for snow/insulation applications.

### Mathematical Formulation

**Heat Equation (Temperature T):**
```
ρcp ∂T/∂t = λ ∇²T + k_c ∇²Cl
```

**Liquid Equation (Liquid Content Cl):**
```
∂Cl/∂t = D_eff(Cl) ∇²Cl + D_t ∇²T + ∇·(D_p ∇p) - m_evap
```

**Vapor Equation (Vapor Content Cv):**
```
∂Cv/∂t = D_v ∇²Cv + m_evap
```

### Enhanced Coupling Mechanisms

1. **Heat-Moisture Coupling**: `k_c ∇²Cl` in heat equation
   - Physical basis: Moisture gradients create internal heat sources/sinks
   - Implementation: `dw_laplace.i.Omega(heat_coupling_mat.k_c, s, Cl)`

2. **Thermal Diffusion (Soret Effect)**: `D_t ∇²T` in liquid equation
   - Physical basis: Temperature gradients drive moisture movement
   - Implementation: `dw_laplace.i.Omega(thermal_mat.D_t, r, T)`

3. **Pressure-Driven Flow**: `∇·(D_p ∇p)` in liquid equation
   - Physical basis: Capillary and osmotic pressures drive liquid transport
   - Pressure: `p = p_capillary(Cl, T) + p_osmotic(Cl, T) - p_atm`
   - Implementation: Pre-calculated divergence via `dw_volume_lvf.i.Omega(pressure_gradient.val, r)`

4. **Shrinkage-Modified Diffusivity**: `D_eff(Cl) = D_l × (1 - α_shrink × Cl)`
   - Physical basis: Higher moisture content reduces effective diffusivity
   - Implementation: Variable diffusivity in `dw_laplace.i.Omega(shrinkage_diffusivity.val, r, Cl)`

5. **Enhanced Evaporation**: 
   ```
   m_evap = k_evap × f_humidity × f_temperature × f_liquid × f_vapor
   ```
   - Physical basis: Multi-factor evaporation model including ambient humidity
   - Implementation: `dw_volume_lvf.i.Omega(evaporation_rate.val, r/w)`

### SfePy Implementation

```python
equations = {
    'Heat': """dw_dot.i.Omega(mat.rho_cp, s, dT/dt)
             + dw_laplace.i.Omega(mat.lam, s, T)                    # λ ∇²T
             + dw_laplace.i.Omega(heat_coupling_mat.k_c, s, Cl)     # k_c ∇²Cl
             = - dw_bc_newton.i.Gamma_Left(h_out_dyn.val, T_out_dyn.val, s, T)
               - dw_bc_newton.i.Gamma_Right(in_fixed.h_in, in_fixed.T_in, s, T)
             """,
             
    'Liquid': """dw_dot.i.Omega(r, dCl/dt)
               + dw_laplace.i.Omega(shrinkage_diffusivity.val, r, Cl) # D_eff ∇²Cl
               + dw_laplace.i.Omega(thermal_mat.D_t, r, T)            # D_t ∇²T
               = - dw_volume_lvf.i.Omega(evaporation_rate.val, r)     # -m_evap
                 - dw_bc_newton.i.Gamma_Left(out_fixed.h_l_out, cl_out_dyn.val, r, Cl)
                 - dw_bc_newton.i.Gamma_Right(in_fixed.h_l_in, in_fixed.Cl_in, r, Cl)
                 + dw_volume_lvf.i.Omega(pressure_gradient.val, r)    # ∇·(D_p ∇p)
               """,
               
    'Vapor': """dw_dot.i.Omega(w, dCv/dt)
              + dw_laplace.i.Omega(vapor_mat.D_v, w, Cv)             # D_v ∇²Cv
              = + dw_volume_lvf.i.Omega(evaporation_rate.val, w)     # +m_evap
                - dw_bc_newton.i.Gamma_Left(h_v_out_dyn.val, cv_out_dyn.val, w, Cv)
                - dw_bc_newton.i.Gamma_Right(in_fixed.h_v_in, in_fixed.Cv_in, w, Cv)
              """
}
```

### Material Properties

```python
materials = {
    'mat': ({'lam': λ, 'rho_cp': ρcp},),                    # Heat capacity and conductivity
    'heat_coupling_mat': ({'k_c': k_c},),                  # Heat-moisture coupling
    'thermal_mat': ({'D_t': D_t},),                        # Thermal diffusion (Soret)
    'vapor_mat': ({'D_v': D_v},),                          # Vapor diffusion
    'shrinkage_diffusivity': 'get_shrinkage_modified_diffusivity', # Variable liquid diffusion
    'pressure_gradient': 'get_pressure_gradient_source',    # Pressure-driven transport
    'evaporation_rate': 'get_evaporation_rate',            # Phase change
}
```

### Boundary Conditions

**Left Boundary (x = 0, Outer Surface):**
- Heat: `-λ ∂T/∂x = h_o(t) [T - T_o(t)]` (Robin BC with time-dependent convection)
- Liquid: `-D_eff ∂Cl/∂x = h_l_o [Cl - Cl_o(t)]` (Moisture exchange with environment)
- Vapor: `-D_v ∂Cv/∂x = h_v_o(t) [Cv - Cv_o(t)]` (Vapor exchange based on humidity)

**Right Boundary (x = d_ins, Inner Surface):**
- Heat: `-λ ∂T/∂x = h_i [T - T_i]` (Fixed interior conditions)
- Liquid: `-D_eff ∂Cl/∂x = h_l_i [Cl - Cl_i]` (Interior moisture conditions)
- Vapor: `-D_v ∂Cv/∂x = h_v_i [Cv - Cv_i]` (Interior vapor conditions)

### Key Corrections Made

1. **Eliminated duplicate terms**: Removed conflicting `dw_volume_lvf` terms that duplicated `dw_laplace` physics
2. **Proper coupling implementation**: Used `dw_laplace` for all diffusion-type coupling terms
3. **Consistent flux formulation**: Based on clay research flux equations
4. **Maintained pressure handling**: Kept pre-calculated pressure gradient as volume source since pressure field not explicitly solved

This formulation properly represents the coupled transport physics while maintaining numerical tractability within the three-field framework.

Looking at this discontinuous Galerkin example, **your understanding and implementation are correct**. Let me explain why:

## How `dw_volume_lvf` Works as a Divergence Operator

In the DG example, you can see:
```python
"- dw_volume_lvf.i.Omega(g.val, v)"
```

Where `g.val` comes from `source_fun` that returns a pre-calculated source term. This is exactly what you're doing with your pressure gradient.

## The Key Insight

`dw_volume_lvf.i.Omega(f, v)` represents the weak form integral:
```
∫_Ω f(x) * v(x) dΩ
```

When `f(x)` is a **pre-calculated divergence** (like `∇·(D_p ∇p)`), then this term effectively adds the divergence to your equation.

## Your Implementation is Correct

Your approach:
```python
+ dw_volume_lvf.i.Omega(pressure_gradient.val, r)  # ∇·(D_p ∇p)
```

Where `pressure_gradient.val` returns the pre-computed `∇·(D_p ∇p)`, is mathematically sound because:

1. **You calculate the full divergence** in `get_pressure_gradient_source()`
2. **You apply it as a volume source** using `dw_volume_lvf`
3. **This gives the same weak form** as if you had solved for pressure explicitly

## Why This Works

The mathematical equivalence is:
```
Method 1 (explicit pressure field):
∫ D_p ∇p · ∇v dΩ  (using dw_laplace)

Method 2 (your approach):
∫ [∇·(D_p ∇p)] · v dΩ  (using dw_volume_lvf with pre-calculated divergence)
```

Both represent the same physics, just with different numerical implementation strategies.

## Validation from the Example

The DG example confirms this pattern:
- Complex differential operators can be pre-calculated
- Volume source terms (`dw_volume_lvf`) can represent these operators
- This is a legitimate numerical technique

## Your Implementation Status

Your current approach is **mathematically correct** and **numerically valid**. The only improvement would be ensuring your finite difference calculation in `get_pressure_gradient_source()` is accurate, but the overall framework is sound.

The DG example validates that `dw_volume_lvf` can indeed act as a divergence operator when you pre-calculate the divergence properly.

Of course. These are all **Robin boundary conditions** (also called mixed boundary conditions of the third kind). They describe a flux that is proportional to the difference between the internal field variable and an external value.

Here they are rewritten in a more formal LaTeX notation.

---

### Left Boundary ($x = 0$, Outer Surface)

At this boundary, the flux of each quantity is proportional to the difference between its value at the surface and an external value.

**Heat Transfer:**
$$
-k_t \frac{\partial T}{\partial x} \bigg|_{x=0} = h_{t,o} \left[ T(0, t) - T_o(t) \right]
$$

**Liquid Transport:**
$$
-D_l \frac{\partial C_l}{\partial x} \bigg|_{x=0} = h_{l,o} \left[ C_l(0, t) - C_{l,o}(t) \right]
$$

**Vapor Transport:**
$$
-D_v \frac{\partial C_v}{\partial x} \bigg|_{x=0} = h_{v,o} \left[ C_v(0, t) - C_{v,o}(t) \right]
$$

**Pressure:**
$$
-k_p \frac{\partial P}{\partial x} \bigg|_{x=0} = h_{p,o} \left[ P(0, t) - P_{\text{atm}} \right]
$$

---

### Right Boundary ($x = d_{\text{ins}}$, Inner Surface)

The same form applies at the right boundary, but with different proportionality constants and external values.

**Heat Transfer:**
$$
-k_t \frac{\partial T}{\partial x} \bigg|_{x=d_{\text{ins}}} = h_{t,i} \left[ T(d_{\text{ins}}, t) - T_i \right]
$$

**Liquid Transport:**
$$
-D_l \frac{\partial C_l}{\partial x} \bigg|_{x=d_{\text{ins}}} = h_{l,i} \left[ C_l(d_{\text{ins}}, t) - C_{l,i} \right]
$$

**Vapor Transport:**
$$
-D_v \frac{\partial C_v}{\partial x} \bigg|_{x=d_{\text{ins}}} = h_{v,i} \left[ C_v(d_{\text{ins}}, t) - C_{v,i} \right]
$$

**Pressure:**
$$
-D_p \frac{\partial P}{\partial x} \bigg|_{x=d_{\text{ins}}} = h_{p,i} \left[ P(d_{\text{ins}}, t) - P_i \right]
$$

---

### Summary

All these equations are of the canonical Robin/mixed BC form:
$$
-\alpha \frac{\partial \psi}{\partial n} = \beta (\psi - \psi_{\text{external}})
$$
where $\alpha$ is a transport coefficient (e.g., $k_t, D_l$), $\frac{\partial \psi}{\partial n}$ is the derivative in the direction of the outward-facing normal, $\beta$ is a transfer coefficient (e.g., $h$), and $\psi_{\text{external}}$ is the external driving value.

Of course. Here are the boundary conditions in LaTeX format.

### Left Boundary ($x = 0$, Outer Surface)

**Heat:**
$$
-\lambda \frac{\partial T}{\partial x} \bigg|_{x=0} = h_o \left( T - T_o \right) + \alpha \cdot I_{\text{solar}}
$$

**Liquid:**
$$
-D_l \frac{\partial C_l}{\partial x} \bigg|_{x=0} = m_{\text{rain}}(t)
$$

**Vapor:**
$$
-D_v \frac{\partial C_v}{\partial x} \bigg|_{x=0} = k_{\text{mass}} \left( \rho_{v,\text{sat}} - \rho_v \right)
$$

**Pressure:**
$$
-k_p \frac{\partial P}{\partial x} \bigg|_{x=0} = h_p \left( P - P_{\text{atm}} \right)
$$

---

### Right Boundary ($x = d_{\text{ins}}$, Inner Surface)

**Heat:**
$$
-\lambda \frac{\partial T}{\partial x} \bigg|_{x=d_{\text{ins}}} = h_i \left( T - 0 \right)
$$

**Liquid:**
$$
-D_l \frac{\partial C_l}{\partial x} \bigg|_{x=d_{\text{ins}}} = k_{\text{drain}} \cdot C_l
$$

**Vapor:**
$$
\frac{\partial C_v}{\partial x} \bigg|_{x=d_{\text{ins}}} = 0
$$

**Pressure:**
$$
P(d_{\text{ins}}, t) = P_{\text{atm}} + \rho_{\text{snow}} \cdot g \cdot d_{\text{ins}}
$$

## Boundary Conditions
For the external drying surfaces of the sample, the boundary conditions are assumed to be of the following form

$$
\begin{aligned}
& \left.\mathbf{J}_{\mathbf{w}}\right|_{x=0^{+}} \cdot \hat{\mathbf{n}}=h_{\mathrm{m}} c M_{\mathrm{v}} \ln \left(\frac{1-x_{\infty}}{1-\left.x_{\mathrm{v}}\right|_{x=0}}\right) \\
& \left.\mathbf{J}_{\mathbf{e}}\right|_{x=0^{+}} \cdot \hat{\mathbf{n}}=h\left(\left.T\right|_{x=0}-T_{\infty}\right) \\
& \left.P_{\mathrm{g}}\right|_{x=0^{+}}=P_{\mathrm{atm}}
\end{aligned}
$$

where $\mathbf{J}_{\mathbf{w}}$ and $\mathbf{J}_{\mathbf{e}}$ represent the fluxes of total moisture and total enthalpy at the boundary, respectively, and x
denotes the normal position from the boundary in the external medium.