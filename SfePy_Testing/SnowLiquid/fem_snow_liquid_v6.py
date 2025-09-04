r"""
==================================================================
Enhanced Heat-Liquid-Vapor Coupled Transport System with 
Additional Moisture Movement Terms Based on Clay Research
------------------------------------------------------------------
Three-field coupled system with enhanced physics
------------------------------------------------------------------
GOVERNING EQUATIONS:
Heat Equation (Temperature T) ,Liquid Equation (Liquid Content Cl),
Vapor Equation (Vapor Content Cv).

    Heat:   ρcp ∂T/∂t = λ ∂²T/∂x² + k_c ∂Cl/∂x
    Liquid: ∂Cl/∂t = D_eff(Cl) ∂²Cl/∂x² + D_p ∂²p/∂x² + D_t ∂²T/∂x² - m_evap
    Vapor:  ∂Cv/∂t = D_v ∂²Cv/∂x² + m_evap
------------------------------------------------------------------
ENHANCED COUPLING MECHANISMS:
1. Pressure-driven liquid flow: 
    p = p_capillary(Cl, T) + p_osmotic(Cl, T) - p_atm
  where:
    - p_capillary: 2σcos(θ)/r_eff(Cl) with temperature-dependent 
      surface tension
    - p_osmotic: Van't Hoff approximation π = cRT
2. Soret effect (thermal diffusion): D_t ∂T/∂x drives moisture 
   movement
3. Heat-moisture coupling: k_c ∂Cl/∂x contributes to heat flux
4. Shrinkage-modified diffusivity: 
    D_eff(Cl) = D_l * (1 - α_shrink * Cl)
5. Enhanced evaporation: m_evap = k_evap * (1 - RH) with 
   humidity dependence
------------------------------------------------------------------
BOUNDARY CONDITIONS:

Left Boundary (x = 0, Outer Surface):
    Heat:        -λ ∂T/∂x = h_o(t) [T - T_o(t)]
    Liquid:      -D_eff ∂Cl/∂x = h_l_o [Cl - Cl_o(t)]  
    Vapor:       -D_v ∂Cv/∂x = h_v_o [Cv - Cv_o(t)]

Right Boundary (x = d_ins, Inner Surface):
    Heat:        -λ ∂T/∂x = h_i [T - T_i]
    Liquid:      -D_eff ∂Cl/∂x = h_l_i [Cl - Cl_i]
    Vapor:       -D_v ∂Cv/∂x = h_v_i [Cv - Cv_i]

Time-dependent outer boundary conditions:
    h_o(t), T_o(t), Cl_o(t), Cv_o(t) from environmental data

Fixed inner boundary conditions:
    h_i = 99.75, T_i = 0°C, Cl_i = 0.01, Cv_i = 0.005
------------------------------------------------------------------
MATERIAL PROPERTIES:
    λ = 0.32 W/(mK), ρcp = ρ_wet * c_wet, D_l = 1e-8 m²/s, D_v = 2e-5 m²/s
    D_p = 1e-9 m²/(Pa·s), D_t = 5e-10 m²/(K·s), k_c = 0.1 W·s/(m·K·kg)
------------------------------------------------------------------
CODE EXECUTION INSTRUCTIONS:

Run this script with:
    sfepy-run fem_snow_liquid_v5.py

For running in parallel, install:
    pip install mpi4py

Run on mulptiple cores with:
    mpirun -n 4 sfepy-run --app=bvp-mM --debug-mpi fem_snow_liquid_v5.py
==================================================================
"""

###############################################################################
###############################################################################
###############################################################################

# TODO! Fix coors.shape problem!!!!

def interpolate_nodes_to_qp_1d(node_values, n_qp):
    """
    Interpolate nodal values to quadrature points for 1D linear elements.
    
    For 1D linear elements, quadrature points are typically at element centers.
    
    Parameters:
    -----------
    node_values : array, shape (n_nodes,)
        Values at mesh nodes
    n_qp : int
        Number of quadrature points (should equal number of elements)
        
    Returns:
    --------
    qp_values : array, shape (n_qp,)
        Values interpolated to quadrature points
    """
    n_nodes = len(node_values)
    n_elements = n_nodes - 1
    
    if n_qp != n_elements:
        raise ValueError(f"Expected {n_elements} quadrature points, got {n_qp}")
    
    qp_values = nm.zeros(n_qp)
    
    # For 1D linear elements with quadrature point at element center
    # QP value = 0.5 * (left_node + right_node)
    for i in range(n_elements):
        qp_values[i] = 0.5 * (node_values[i] + node_values[i + 1])
    
    return qp_values

def get_evaporation_rate_corrected(ts, coors, mode=None, equations=None, term=None, 
                                 problem=None, **kwargs):
    """
    Corrected evaporation rate function with proper node-to-QP interpolation.
    """
    if mode != 'qp' or coors is None:
        return {}
    
    # Get current relative humidity from environmental data
    hour_idx = min(int(ts.time / 3600), len(rh) - 1)
    current_rh = rh[hour_idx]
    humidity_factor = max(0, 1.0 - current_rh/100.0)
    
    # Get current state variables (these are at NODES)
    variables = problem.get_variables()
    T_var = variables['T']
    Cl_var = variables['Cl']
    Cv_var = variables['Cv']
    
    # Get state vectors at nodes
    T_vals_nodes = T_var.get_state_in_region(problem.domain.regions['Omega'])
    Cl_vals_nodes = Cl_var.get_state_in_region(problem.domain.regions['Omega'])
    Cv_vals_nodes = Cv_var.get_state_in_region(problem.domain.regions['Omega'])
    
    n_qp = coors.shape[0]  # Number of quadrature points
    
    print(f"Debug: Nodes: {len(T_vals_nodes)}, QP: {n_qp}")
    
    try:
        # Interpolate from nodes to quadrature points
        T_vals = interpolate_nodes_to_qp_1d(T_vals_nodes, n_qp)
        Cl_vals = interpolate_nodes_to_qp_1d(Cl_vals_nodes, n_qp)  
        Cv_vals = interpolate_nodes_to_qp_1d(Cv_vals_nodes, n_qp)
        
        print(f"Interpolation successful: T range = [{T_vals.min():.3f}, {T_vals.max():.3f}]")
        
    except Exception as e:
        print(f"Interpolation failed: {e}")
        print("Using element-wise averages as fallback")
        
        # Fallback: use averages (your current approach)
        T_vals = nm.full(n_qp, T_vals_nodes.mean())
        Cl_vals = nm.full(n_qp, Cl_vals_nodes.mean()) 
        Cv_vals = nm.full(n_qp, Cv_vals_nodes.mean())
    
    # Calculate evaporation factors (now at quadrature points)
    T_factor = nm.maximum(0, nm.minimum(1, (T_vals - T_ref)/(T_evap_max - T_ref)))
    Cl_factor = nm.maximum(0, (Cl_vals - Cl_min)/Cl_min)
    
    T_K = T_vals + 273.15
    Psat = Psat_WV(T_K) * 100.0
    M_wv = 0.018
    R_gas = 8.314
    rho_v_sat = Psat * M_wv / (R_gas * T_K)
    vapor_factor = nm.maximum(0, 1.0 - Cv_vals / rho_v_sat)
    
    # Combined evaporation rate
    evap_rate = k_evap * humidity_factor * T_factor * Cl_factor * vapor_factor
    
    print(f"Evaporation rate range: [{evap_rate.min():.2e}, {evap_rate.max():.2e}]")
    
    # Should now have correct shape
    assert len(evap_rate) == coors.shape[0], f"Shape mismatch: {len(evap_rate)} != {coors.shape[0]}"
    
    val = evap_rate.reshape((coors.shape[0], 1, 1))
    return {'val': val}

def get_pressure_gradient_source_corrected(ts, coors, mode=None, equations=None, term=None,
                                         problem=None, **kwargs):
    """
    Corrected pressure gradient function with proper interpolation.
    """
    if mode != 'qp' or coors is None:
        return {}
    
    try:
        if problem is not None and hasattr(problem, 'get_variables'):
            variables = problem.get_variables()
            
            if 'T' in variables and 'Cl' in variables:
                T_var = variables['T']
                Cl_var = variables['Cl']
                
                # Get nodal values
                T_vals_nodes = T_var.get_state_in_region(problem.domain.regions['Omega'])
                Cl_vals_nodes = Cl_var.get_state_in_region(problem.domain.regions['Omega'])
                
                n_qp = coors.shape[0]
                
                # Interpolate to quadrature points
                T_vals = interpolate_nodes_to_qp_1d(T_vals_nodes, n_qp)
                Cl_vals = interpolate_nodes_to_qp_1d(Cl_vals_nodes, n_qp)
                
                # Calculate pressures at QP
                p_cap = calculate_capillary_pressure(Cl_vals, T_vals)
                p_osm = calculate_osmotic_pressure(Cl_vals, T_vals)
                total_pressure = p_cap + p_osm - p_atm
                
                # Pressure gradient effect
                pressure_effect = D_p * total_pressure / (mu_water * dx**2)
                
            else:
                pressure_effect = nm.zeros(coors.shape[0])
        else:
            pressure_effect = nm.zeros(coors.shape[0])
            
    except Exception as e:
        print(f"Pressure calculation failed: {e}")
        pressure_effect = nm.zeros(coors.shape[0])
    
    val = pressure_effect.reshape((coors.shape[0], 1, 1))
    return {'val': val}

###############################################################################
###############################################################################
###############################################################################

from __future__ import absolute_import
import numpy as nm
from sfepy.discrete.fem import Mesh
from sfepy.discrete.fem.meshio import UserMeshIO
import csv
import os



# Physical parameters
d_ins = 0.1       # Insulation thickness [m]
dx = 0.005        # Cell size [m]
n_el = int(d_ins / dx)  # Number of elements
nodes = n_el + 1        # Number of nodes

# Enhanced material properties
lam_i = 0.32      # Thermal conductivity [W/(mK)]
c_wet = 2.59E03   # J/(kg*K)
rho_dry = 100.0   # kg/m^3
moist_cont = 50.0 # %
rho_wet = rho_dry + moist_cont*10.0 # kg/m^3
D = lam_i/(c_wet * rho_wet) # m^2/s

# Enhanced transport coefficients
D_l = 1e-8        # Liquid diffusion coefficient [m^2/s]
D_v = 2e-5        # Vapor diffusion coefficient [m^2/s]
D_p = 1e-9        # Pressure diffusion coefficient [m^2/(Pa*s)]
D_t = 5e-10       # Thermal diffusion coefficient [m^2/(K*s)] - Soret effect
k_c = 0.1         # Heat-moisture coupling coefficient [W*s/(m*K*kg)]

# Pressure and suction parameters
k_s = 1e-15       # Permeability [m^2]
mu_water = 1e-3   # Dynamic viscosity of water [Pa*s]
p_atm = 101325    # Atmospheric pressure [Pa]
gamma_surface = 0.072  # Surface tension [N/m]
contact_angle = 60     # Contact angle [degrees]

# Additional physical constants
L_v = 2.45e6      # Latent heat of vaporization [J/kg]
R_v = 461.5       # Specific gas constant for water vapor [J/(kg*K)]
alpha_shrink = 0.01   # Shrinkage coefficient [1/kg_water]

# Enhanced evaporation model parameters
k_evap = 1e-7     # Evaporation rate constant [1/s]
T_ref = 0.0       # Reference temperature [°C]
Cl_min = 0.001    # Minimum liquid content for evaporation [kg/m³]
T_evap_max = 10.0 # Maximum temperature for linear evaporation model [°C]

def Psat_WV(T_K):
    """Water vapour saturation pressure (same as before)"""
    Tc = 647.096  # Critical temperature, K
    Pc = 220640   # Critical pressure, hPa
    C1 = -7.85951783
    C2 = 1.84408259
    C3 = -11.7866497
    C4 = 22.6807411
    C5 = -15.9618719
    C6 = 1.80122502
    teta = 1 - T_K / Tc
    x = Tc / T_K * (C1*teta + C2*teta**1.5 + C3*teta**3 + C4*teta**3.5 + \
                    C5*teta**4 + C6*teta**7.5)
    x = nm.exp(x) * Pc
    return x

def calculate_capillary_pressure(Cl, T):
    """
    Calculate capillary pressure based on moisture content and temperature
    Based on Young-Laplace equation and pore structure
    """
    # Simplified model: p_cap = f(moisture_content, temperature)
    # Higher moisture -> lower capillary pressure (larger effective pore radius)
    Cl_norm = nm.maximum(0.001, Cl / 0.1)  # Normalize to prevent division by zero
    
    # Temperature effect on surface tension
    T_celsius = T
    sigma_temp = gamma_surface * (1 - 0.0015 * T_celsius)  # Approximate temperature dependence
    
    # Effective capillary radius (increases with moisture content)
    r_eff = 1e-6 * (Cl_norm ** 0.3)  # Effective radius [m]
    
    # Young-Laplace equation: ΔP = 2σcos(θ)/r
    cos_theta = nm.cos(nm.radians(contact_angle))
    p_cap = 2 * sigma_temp * cos_theta / r_eff
    
    return p_cap

def calculate_osmotic_pressure(Cl, T):
    """
    Calculate osmotic pressure based on moisture content
    Simplified model for clay-water interaction
    """
    # Van't Hoff equation approximation for dilute solutions
    # π = cRT where c is concentration
    T_K = T + 273.15
    concentration = Cl / 1000.0  # Approximate concentration [mol/L]
    pi_osmotic = concentration * R_v * T_K * 0.001  # Convert to Pa
    
    return pi_osmotic

# NOTE! UNUSED!!!!
def get_pressure_gradient_source(ts, coors, mode=None, equations=None, term=None,
                                problem=None, **kwargs):
    """
    Calculate pressure-driven moisture flux term: ∇·(D_p ∇p)
    where p = p_capillary + p_osmotic
    """
    if mode != 'qp' or coors is None:
        return {}
    
    try:
        if problem is not None and hasattr(problem, 'get_variables'):
            variables = problem.get_variables()
            
            if 'T' in variables and 'Cl' in variables:
                T_var = variables['T']
                Cl_var = variables['Cl']
                
                T_vals = T_var.get_state_in_region(problem.domain.regions['Omega'])
                Cl_vals = Cl_var.get_state_in_region(problem.domain.regions['Omega'])
                
                # Calculate pressures
                p_cap = calculate_capillary_pressure(Cl_vals, T_vals)
                p_osm = calculate_osmotic_pressure(Cl_vals, T_vals)
                total_pressure = p_cap + p_osm - p_atm
                
                # Pressure gradient effect (simplified as pressure magnitude for volume term)
                pressure_effect = D_p * total_pressure / (mu_water * dx**2)
                
                print("Pressure effect", pressure_effect, pressure_effect.shape)

                if len(pressure_effect) != coors.shape[0]:
                    pressure_effect = nm.full(coors.shape[0], pressure_effect.mean())
            else:
                pressure_effect = nm.zeros(coors.shape[0])
        else:
            pressure_effect = nm.zeros(coors.shape[0])
            
    except Exception as e:
        pressure_effect = nm.zeros(coors.shape[0])
    
    val = pressure_effect.reshape((coors.shape[0], 1, 1))
    return {'val': val}

def get_soret_effect(ts, coors, mode=None, equations=None, term=None,
                    problem=None, **kwargs):
    """
    Calculate Soret effect (thermal diffusion): D_t * ∇T
    This represents moisture movement due to temperature gradients
    """
    if mode != 'qp' or coors is None:
        return {}
    
    try:
        if problem is not None and hasattr(problem, 'get_variables'):
            variables = problem.get_variables()
            
            if 'T' in variables:
                T_var = variables['T']
                T_vals = T_var.get_state_in_region(problem.domain.regions['Omega'])
                
                # Approximate temperature gradient (simplified for 1D)
                if len(T_vals) > 1:
                    T_grad = nm.gradient(T_vals) / dx
                    # Thermal diffusion effect
                    soret_effect = D_t * T_grad
                    
                    if len(soret_effect) != coors.shape[0]:
                        soret_effect = nm.full(coors.shape[0], soret_effect.mean())
                else:
                    soret_effect = nm.zeros(coors.shape[0])
            else:
                soret_effect = nm.zeros(coors.shape[0])
        else:
            soret_effect = nm.zeros(coors.shape[0])
            
    except Exception as e:
        soret_effect = nm.zeros(coors.shape[0])
    
    val = soret_effect.reshape((coors.shape[0], 1, 1))
    return {'val': val}

def get_heat_moisture_coupling(ts, coors, mode=None, equations=None, term=None,
                             problem=None, **kwargs):
    """
    Heat-moisture coupling term: k_c * ∇C
    This represents heat flux due to moisture gradients
    """
    if mode != 'qp' or coors is None:
        return {}
    
    try:
        if problem is not None and hasattr(problem, 'get_variables'):
            variables = problem.get_variables()
            
            if 'Cl' in variables:
                Cl_var = variables['Cl']
                Cl_vals = Cl_var.get_state_in_region(problem.domain.regions['Omega'])
                
                # Approximate moisture gradient (simplified for 1D)
                if len(Cl_vals) > 1:
                    Cl_grad = nm.gradient(Cl_vals) / dx
                    # Heat-moisture coupling effect
                    coupling_effect = k_c * Cl_grad / (rho_wet * c_wet)  # Normalize
                    print("coupling effect", coupling_effect, coupling_effect.shape)
                    if len(coupling_effect) != coors.shape[0]:
                        coupling_effect = nm.full(coors.shape[0], coupling_effect.mean())
                else:
                    coupling_effect = nm.zeros(coors.shape[0])
            else:
                coupling_effect = nm.zeros(coors.shape[0])
        else:
            coupling_effect = nm.zeros(coors.shape[0])
            
    except Exception as e:
        coupling_effect = nm.zeros(coors.shape[0])
    
    val = coupling_effect.reshape((coors.shape[0], 1, 1))
    return {'val': val}

def get_heat_moist_coup(ts, coors, mode=None, equations=None, term=None,
                             problem=None, **kwargs):
    """
    Heat-moisture coupling term: k_c * ∇C
    This represents heat flux due to moisture gradients
    """
    if mode != 'qp' or coors is None:
        return {}
    
    # Get variables
    variables = problem.get_variables()
    Cl_var = variables['Cl']
    Cl_nodal = Cl_var()  # Get all 21 nodal values
    
    # Calculate gradient at NODES (where we have enough points)
    nodal_coors = nm.linspace(0.0, d_ins, 21)
    dx_nodal = nodal_coors[1] - nodal_coors[0]  # Actual nodal spacing
    
    Cl_grad_nodal = nm.gradient(Cl_nodal) / dx_nodal  # Gradient at 21 nodes
    
    # Now interpolate gradient to quadrature points
    Cl_grad_qp = 0.5 * (Cl_grad_nodal[:-1] + Cl_grad_nodal[1:])  # 20 values
    
    # Heat-moisture coupling effect
    coupling_effect = k_c * Cl_grad_qp / (rho_wet * c_wet)
    
    val = coupling_effect.reshape((coors.shape[0], 1, 1))
    return {'val': val}

def get_shrinkage_modified_diffusivity(ts, coors, mode=None, equations=None, 
                                       term=None, problem=None, **kwargs):
                                    
    """
    Shrinkage-modified diffusion coefficient
    D_eff = D_l * (1 - α_shrink * w)
    """
    if mode != 'qp' or coors is None:
        return {}
    
    try:
        if problem is not None and hasattr(problem, 'get_variables'):
            variables = problem.get_variables()
            
            if 'Cl' in variables:
                Cl_var = variables['Cl']
                Cl_vals = Cl_var.get_state_in_region(problem.domain.regions['Omega'])
                
                # Calculate effective diffusivity considering shrinkage
                # Higher moisture -> more shrinkage -> reduced diffusivity
                shrinkage_factor = nm.maximum(0.1, 1.0 - alpha_shrink * Cl_vals)
                D_eff = D_l * shrinkage_factor
                
                if len(D_eff) != coors.shape[0]:
                    D_eff = nm.full(coors.shape[0], D_eff.mean())
            else:
                D_eff = nm.full(coors.shape[0], D_l)
        else:
            D_eff = nm.full(coors.shape[0], D_l)
            
    except Exception as e:
        D_eff = nm.full(coors.shape[0], D_l)
    
    val = D_eff.reshape((coors.shape[0], 1, 1))
    return {'val': val}

def get_shrink_mod_diff(ts, coors, mode=None, equations=None,
                       term=None, problem=None, **kwargs):
    if mode != 'qp' or coors is None:
        return {}
    
    # Get nodal values
    variables = problem.get_variables()
    Cl_var = variables['Cl']
    Cl_nodal = Cl_var()  # 21 values at nodes
    
    # INTERPOLATE to quadrature points
    # For 1D linear elements with 1 quadrature point per element:
    # The quadrature point is at the center, so average the two nodes
    Cl_qp = 0.5 * (Cl_nodal[:-1] + Cl_nodal[1:])  # Now 20 values!
    print("Cl_qp shape", Cl_qp.shape)
    
    shrinkage_factor = nm.maximum(0.1, 1.0 - alpha_shrink * Cl_qp)
    D_eff = D_l * shrinkage_factor  # Now 20 values - CORRECT!
    print("D_eff shape", D_eff.shape)

    # NOW the shapes match: 20 quadrature points = 20 material values
    val = D_eff.reshape((coors.shape[0], 1, 1))  # (20, 1, 1)
    return {'val': val}

# Keep all existing functions (read_input_data, get_h_o, etc.) from original code
def read_input_data(filename="DATA.csv"):
    """Read air temperature, air velocity, precipitation, global solar irradation 
    and relative humidity data from CSV file."""
    airTemp, airVel, prec, gloSolIr, relHum = [], [], [], [], []
    with open('DATA.csv', 'r') as csvfile:
        csvreader = csv.reader(csvfile)
        next(csvreader, None)
        for row in csvreader:
            airTemp.append(float(row[0].strip()))
            airVel.append(float(row[1].strip()))
            prec.append(float(row[2].strip()))
            gloSolIr.append(float(row[3].strip()))
            relHum.append(float(row[4].strip())) 
    return airTemp, airVel, prec, gloSolIr, relHum

# Keep existing boundary condition functions (get_h_o, get_t_o, etc.)
def get_h_o(ts, coors, mode=None, **kwargs):
    """Time-dependent heat transfer coefficient."""
    if mode != 'qp' or coors is None: return {}
    hour_idx = min(int(ts.time / 3600), len(h_o) - 1)
    val = nm.full((coors.shape[0], 1, 1), h_o[hour_idx], dtype=nm.float64)
    return {'val': val}

def get_h_v_o(ts, coors, mode=None, **kwargs):
    """Time-dependent vapor heat transfer coefficient."""
    if mode != 'qp' or coors is None: return {}
    hour_idx = min(int(ts.time / 3600), len(h_v_o) - 1)
    val = nm.full((coors.shape[0], 1, 1), h_v_o[hour_idx], dtype=nm.float64)
    return {'val': val}

def get_t_o(ts, coors, mode=None, **kwargs):
    """Time-dependent sol-air temperature."""
    if mode != 'qp' or coors is None: return {}
    hour_idx = min(int(ts.time / 3600), len(t_o) - 1)
    val = nm.full((coors.shape[0], 1, 1), t_o[hour_idx], dtype=nm.float64)
    return {'val': val}

# Keep other existing functions...
def get_cl_o(ts, coors, mode=None, **kwargs):
    """Time-dependent outer liquid reference content."""
    if mode != 'qp' or coors is None: return {}
    hour_idx = min(int(ts.time / 3600), len(rh) - 1)
    current_rh = rh[hour_idx] if hour_idx < len(rh) else 70.0
    humidity_factor = current_rh / 100.0
    cl_adjusted = cl_o * (0.5 + 0.5 * humidity_factor)
    val = nm.full((coors.shape[0], 1, 1), cl_adjusted, dtype=nm.float64)
    return {'val': val}

def get_cv_o_old(ts, coors, mode=None, **kwargs):
    """Time-dependent outer vapor reference content."""
    if mode != 'qp' or coors is None: return {}
    hour_idx = min(int(ts.time / 3600), len(rh) - 1)
    current_rh = rh[hour_idx] if hour_idx < len(rh) else 70.0
    current_t = t_o[hour_idx] if hour_idx < len(t_o) else -5.0
    T_K = current_t + 273.15
    try:
        Psat = Psat_WV(T_K)
        rho_v_sat = Psat * 100 * 0.018 / (8.314 * T_K)
        rho_v_actual = rho_v_sat * current_rh / 100.0
        cv_adjusted = max(cv_o * 0.1, min(cv_o * 2.0, rho_v_actual))
    except:
        cv_adjusted = cv_o * (current_rh / 100.0)
    val = nm.full((coors.shape[0], 1, 1), cv_adjusted, dtype=nm.float64)
    return {'val': val}

def get_cv_i(ts, coors, mode=None, **kwargs):
    """Time-dependent inner vapor reference content."""
    if mode != 'qp' or coors is None: return {}
    val = nm.full((coors.shape[0], 1, 1), cv_i, dtype=nm.float64)
    return {'val': val}

def get_vapor_content_from_RH(T_celsius, RH_percent):
    """
    Calculate vapor content from temperature and relative humidity.
    
    Uses ideal gas law: ρ_v = (P_v * M_v) / (R * T)
    where P_v = RH * P_sat
    
    Parameters:
    -----------
    T_celsius : float or array
        Temperature in Celsius
    RH_percent : float or array  
        Relative humidity in percent
        
    Returns:
    --------
    vapor_content : float or array
        Vapor content in kg/m³
    """
    T_K = T_celsius + 273.15
    
    # Saturation pressure
    P_sat = Psat_WV(T_K)  # Pa
    
    # Actual vapor pressure
    P_vapor = P_sat * RH_percent / 100.0  # Pa
    
    # Ideal gas law: ρ = P*M/(R*T)
    M_water = 0.018015  # kg/mol, molar mass of water
    R_gas = 8.314       # J/(mol·K), universal gas constant
    
    vapor_content = (P_vapor * M_water) / (R_gas * T_K)  # kg/m³
    
    return vapor_content

def get_cv_o(ts, coors, mode=None, **kwargs):
    """
    Time-dependent outer vapor reference content based on environmental 
    RH and temperature.
    
    This replaces the existing get_cv_o() function with proper physics.
    """
    if mode != 'qp' or coors is None:
        return {}
    
    hour_idx = min(int(ts.time / 3600), len(rh) - 1)
    
    # Get environmental conditions
    current_rh = rh[hour_idx] if hour_idx < len(rh) else 70.0
    current_t = t_o[hour_idx] if hour_idx < len(t_o) else -5.0
    
    try:
        # Calculate vapor content from thermodynamics
        cv_physical = get_vapor_content_from_RH(current_t, current_rh)
        
        # Apply reasonable bounds to prevent numerical issues
        # But keep them physically meaningful
        cv_min = 1e-6   # Very dry conditions (practically zero)
        cv_max = 0.05   # Very humid conditions (near saturation at 30°C)
        
        cv_adjusted = nm.clip(cv_physical, cv_min, cv_max)
        
    except Exception as e:
        # Fallback to simple approximation if calculation fails
        print(f"Warning: Vapor content calculation failed, using approximation: {e}")
        cv_adjusted = 0.005 * (current_rh / 100.0)
    
    val = nm.full((coors.shape[0], 1, 1), cv_adjusted, dtype=nm.float64)
    return {'val': val}

def get_rh_for_vapor_bc(ts, coors, mode=None, **kwargs):
    """
    Time-dependent relative humidity for direct use in boundary conditions.
    
    This can be used if you want to formulate the BC differently.
    """
    if mode != 'qp' or coors is None:
        return {}
    
    hour_idx = min(int(ts.time / 3600), len(rh) - 1)
    current_rh = rh[hour_idx] if hour_idx < len(rh) else 70.0
    
    val = nm.full((coors.shape[0], 1, 1), current_rh, dtype=nm.float64)
    return {'val': val}

def calculate_vapor_bc_coefficient(T_celsius, RH_percent):
    """
    Calculate the vapor boundary condition coefficient for Robin BC.
    
    For Robin BC: -D_v ∂Cv/∂x = h_v (Cv - Cv_ref)
    
    The coefficient h_v should be related to mass transfer coefficient
    from convective mass transfer theory.
    
    Parameters:
    -----------
    T_celsius : float
        Temperature in Celsius
    RH_percent : float
        Relative humidity in percent
        
    Returns:
    --------
    h_v_eff : float
        Effective vapor transfer coefficient [m/s]
    """
    # Base mass transfer coefficient (can be related to heat transfer coefficient)
    # Analogy: h_m = h_h / (ρ * cp * Le^(2/3)) where Le is Lewis number
    
    h_base = 2e-5  # Base mass transfer coefficient [m/s]
    
    # Temperature effect (higher T → higher molecular motion → higher h_v)
    #T_factor = 1 + 0.01 * T_celsius
    T_factor = [1 + 0.01 * T for T in T_celsius]
    
    # Humidity effect (higher RH → lower driving force → need higher h_v for same flux)
    # This accounts for the nonlinear nature of vapor transport
    #RH_factor = 1 + 0.5 * (RH_percent / 100.0)
    RH_factor = [1 + 0.5 * (RH / 100.0) for RH in RH_percent]

    
    h_v_eff = [h_base * T_f * RH_f for T_f, RH_f in 
               zip(T_factor, RH_factor)]
    
    return h_v_eff

def get_evaporation_rate(ts, coors, mode=None, equations=None, term=None, 
                        problem=None, **kwargs):
    """
    Enhanced temperature and liquid-content dependent evaporation rate.
    Includes humidity effects using SfePy API.
    
    Evaporation model:
    m_evap = k_evap * f_humidity * f_temperature * f_liquid * f_vapor
    
    where:
    - f_humidity = (1 - RH/100) accounts for ambient humidity
    - f_temperature = max(0, (T - T_ref)/(T_evap_max - T_ref)) 
    - f_liquid = max(0, (Cl - Cl_min)/Cl_min) accounts for available liquid
    - f_vapor = 1 - Cv/rho_v_sat limits evaporation near saturation
    """
    if mode != 'qp' or coors is None:
        return {}
    
    # Get current relative humidity from environmental data
    hour_idx = min(int(ts.time / 3600), len(rh) - 1)
    current_rh = rh[hour_idx]
    humidity_factor = max(0, 1.0 - current_rh/100.0)
    print("humidity_factor", humidity_factor)
    
    # Get current state variables
    variables = problem.get_variables()
    T_var = variables['T']
    Cl_var = variables['Cl']
    Cv_var = variables['Cv']
    
    # Get state vectors
    T_vals = T_var.get_state_in_region(problem.domain.regions['Omega'])
    Cl_vals = Cl_var.get_state_in_region(problem.domain.regions['Omega'])
    Cv_vals = Cv_var.get_state_in_region(problem.domain.regions['Omega'])
    
    # Temperature factor: linear increase from T_ref to T_evap_max
    # Essentially, when T < 0, no evaporation can occur.
    T_factor = nm.maximum(0, nm.minimum(1, (T_vals - T_ref)/(T_evap_max - T_ref)))
    print("T_factor", T_factor) # zeros

    # Liquid availability factor
    Cl_factor = nm.maximum(0, (Cl_vals - Cl_min)/Cl_min)
    print("Cl factor", Cl_factor)
    
    # Vapor saturation factor (reduces evaporation when vapor content is high)
    T_K = T_vals + 273.15
    Psat = Psat_WV(T_K) * 100.0  # hPa to Pa
    M_wv = 0.018     # kg/mol, molar mass of water vapor
    R_gas = 8.314    # J/(mol*K)
    rho_v_sat = Psat * M_wv / (R_gas * T_K)  # kg/m³
    print("rho_v_sat", rho_v_sat)

    vapor_factor = nm.maximum(0, 1.0 - Cv_vals / rho_v_sat)
    print("vapor factor", vapor_factor)
    
    # Combined evaporation rate
    evap_rate = k_evap * humidity_factor * T_factor * Cl_factor * vapor_factor
    
    print("\nEvaporation rate", evap_rate, evap_rate.shape)

    # Ensure we have the right shape for quadrature points
    if len(evap_rate) != coors.shape[0]:
        print(f"WARNING: Dimension mismatch! Nodes: {len(evap_rate)}, QP: {coors.shape[0]}")
        print("Using average evaporation rate as fallback")
        evap_rate = nm.full(coors.shape[0], evap_rate.mean())
    
    val = evap_rate.reshape((coors.shape[0], 1, 1))
    return {'val': val}

def get_evap_rate(ts, coors, mode=None, equations=None, term=None, 
                        problem=None, **kwargs):
    """
    Enhanced temperature and liquid-content dependent evaporation rate.
    Includes humidity effects using SfePy API.
    
    Evaporation model:
    m_evap = k_evap * f_humidity * f_temperature * f_liquid * f_vapor
    
    where:
    - f_humidity = (1 - RH/100) accounts for ambient humidity
    - f_temperature = max(0, (T - T_ref)/(T_evap_max - T_ref)) 
    - f_liquid = max(0, (Cl - Cl_min)/Cl_min) accounts for available liquid
    - f_vapor = 1 - Cv/rho_v_sat limits evaporation near saturation
    """
    if mode != 'qp' or coors is None:
        return {}
    
    # Get current relative humidity from environmental data
    hour_idx = min(int(ts.time / 3600), len(rh) - 1)
    current_rh = rh[hour_idx]
    humidity_factor = max(0, 1.0 - current_rh/100.0)
    
    # Get current state variables
    variables = problem.get_variables()
    T_var = variables['T']
    Cl_var = variables['Cl']
    Cv_var = variables['Cv']
    
    # Get state vectors
    T_vals = T_var.get_state_in_region(problem.domain.regions['Omega'])   # 21 points
    Cl_vals = Cl_var.get_state_in_region(problem.domain.regions['Omega']) # 21 points
    Cv_vals = Cv_var.get_state_in_region(problem.domain.regions['Omega']) # 21 points
    
    # INTERPOLATE to quadrature points (average of adjacent nodes)
    T_qp = 0.5 * (T_vals[:-1] + T_vals[1:])  # 20 values
    Cl_qp = 0.5 * (Cl_vals[:-1] + Cl_vals[1:])  # 20 values  
    Cv_qp = 0.5 * (Cv_vals[:-1] + Cv_vals[1:])  # 20 values
    
    # Now do calculations with 20 QUADRATURE point values
    T_factor = nm.maximum(0, nm.minimum(1, (T_qp - T_ref)/(T_evap_max - T_ref)))
    Cl_factor = nm.maximum(0, (Cl_qp - Cl_min)/Cl_min)
    
    T_K = T_qp + 273.15
    Psat = Psat_WV(T_K) * 100.0
    rho_v_sat = Psat * 0.018 / (8.314 * T_K)
    vapor_factor = nm.maximum(0, 1.0 - Cv_qp / rho_v_sat)
    
    evap_rate = k_evap * humidity_factor * T_factor * Cl_factor * vapor_factor  # 20 values
    
    print("\nEvaporation rate", evap_rate, evap_rate.shape)

    # Ensure we have the right shape for quadrature points
    if len(evap_rate) != coors.shape[0]:
        print(f"WARNING: Dimension mismatch! Nodes: {len(evap_rate)}, QP: {coors.shape[0]}")
        print("Using average evaporation rate as fallback")
        evap_rate = nm.full(coors.shape[0], evap_rate.mean())
    
    val = evap_rate.reshape((coors.shape[0], 1, 1))
    return {'val': val}

def calculate_enhanced_evaporation(Cl_vals, Cv_vals, T_vals, rh_current):
    """
    Enhanced evaporation model considering all physical effects.
    
    Physical basis:
    m_evap = k_evap * f_kinetic * f_equilibrium * f_transport
    
    where:
    - f_kinetic: Kinetic theory (temperature dependence)  
    - f_equilibrium: Phase equilibrium (humidity dependence)
    - f_transport: Transport limitations (pore structure, diffusion)
    
    Parameters:
    -----------
    Cl_vals : array
        Liquid moisture content [kg/m³]
    Cv_vals : array
        Vapor moisture content [kg/m³]
    T_vals : array
        Temperature [°C]
    rh_current : float
        Current relative humidity [%]
        
    Returns:
    --------
    evap_rate : array
        Evaporation rate [kg/(m³·s)]
    """
    T_K = T_vals + 273.15
    
    # Kinetic factor: Hertz-Knudsen equation influence
    f_kinetic = nm.exp((T_vals - 0.0) / 20.0)  # Exponential temperature dependence
    
    # Equilibrium factor: Driving force for phase change
    Psat = Psat_WV(T_K) * 100.0  # Pa
    M_wv = 0.018  # kg/mol
    R = 8.314     # J/(mol·K)
    rho_v_sat = Psat * M_wv / (R * T_K)  # Saturation vapor density
    
    rho_v_ambient = rho_v_sat * rh_current / 100.0
    f_equilibrium = nm.maximum(0, (rho_v_sat - Cv_vals) / rho_v_sat)
    
    # Transport factor: Pore structure and diffusion limitations
    # Available liquid water
    f_liquid = nm.maximum(0, (Cl_vals - 0.001) / nm.maximum(0.001, Cl_vals))
    
    # Pore connectivity (decreases as moisture decreases)
    connectivity = nm.minimum(1.0, Cl_vals / 0.05)  # Normalized to typical max moisture
    
    # Diffusion limitation in partially saturated media
    tortuosity = 1.0 + 2.0 * (1.0 - connectivity)
    f_transport = connectivity / tortuosity
    
    # Combined evaporation rate
    evap_rate = k_evap * f_kinetic * f_equilibrium * f_transport * f_liquid
    
    return evap_rate

def mesh_hook(mesh, mode):
    """Generate the 1D mesh."""
    if mode == 'read':
        coors = nm.linspace(0.0, d_ins, nodes).reshape((nodes, 1))
        print("OG coors shape", coors.shape)
        conn = nm.arange(nodes, dtype=nm.int32).repeat(2)[1:-1].reshape((-1, 2))
        mat_ids = nm.zeros(nodes - 1, dtype=nm.int32)
        descs = ['1_2']
        mesh = Mesh.from_data('heat_liquid_1d', coors, None, [conn], [mat_ids], descs)
        return mesh
    elif mode == 'write':
        pass

def calculate_pressure_driven_flux(Cl_vals, T_vals, dx):
    """
    Calculate pressure-driven moisture flux using Darcy's law.
    
    Physical basis:
    q = -(k/μ) * ∇p where k is permeability, μ is viscosity
    
    For moisture transport:
    ∂Cl/∂t = ∇·(D_p ∇p) where D_p = k/(μ φ) and φ is porosity
    
    Parameters:
    -----------
    Cl_vals : array
        Liquid moisture content [kg/m³]
    T_vals : array
        Temperature [°C]
    dx : float
        Grid spacing [m]
    
    Returns:
    --------
    flux_divergence : array
        Divergence of pressure-driven flux [kg/(m³·s)]
    """
    n_nodes = len(Cl_vals)
    flux_divergence = nm.zeros_like(Cl_vals)
    
    # Calculate total pressure at each node
    p_cap = calculate_capillary_pressure(Cl_vals, T_vals)
    p_osm = calculate_osmotic_pressure(Cl_vals, T_vals)
    p_total = p_cap + p_osm
    
    # Calculate moisture-dependent effective diffusivity
    # Higher moisture content → higher effective permeability
    phi_eff = 0.3 + 0.2 * (Cl_vals / nm.max(Cl_vals))  # Effective porosity
    mu_eff = 1e-3 * (1 + 0.1 * Cl_vals)  # Viscosity increases with moisture
    D_p_eff = k_s / (mu_eff * phi_eff)  # Effective pressure diffusivity
    
    # Calculate flux at cell faces using upwind scheme
    for i in range(n_nodes):
        if i == 0:
            # Forward difference at left boundary
            if n_nodes > 1:
                dp_dx = (p_total[1] - p_total[0]) / dx
                D_avg = 0.5 * (D_p_eff[0] + D_p_eff[1])
                flux_right = -D_avg * dp_dx
                flux_divergence[i] = -flux_right / dx
        elif i == n_nodes - 1:
            # Backward difference at right boundary
            dp_dx = (p_total[i] - p_total[i-1]) / dx
            D_avg = 0.5 * (D_p_eff[i-1] + D_p_eff[i])
            flux_left = -D_avg * dp_dx
            flux_divergence[i] = flux_left / dx
        else:
            # Central differences for interior nodes
            # Left face flux
            dp_dx_left = (p_total[i] - p_total[i-1]) / dx
            D_avg_left = 0.5 * (D_p_eff[i-1] + D_p_eff[i])
            flux_left = -D_avg_left * dp_dx_left
            
            # Right face flux
            dp_dx_right = (p_total[i+1] - p_total[i]) / dx
            D_avg_right = 0.5 * (D_p_eff[i] + D_p_eff[i+1])
            flux_right = -D_avg_right * dp_dx_right
            
            # Flux divergence
            flux_divergence[i] = (flux_right - flux_left) / dx
    
    return flux_divergence

def calculate_soret_effect_proper(Cl_vals, T_vals, dx):
    """
    Calculate Soret effect (thermal diffusion) with proper physics.
    
    Physical basis:
    The Soret effect causes mass flux due to temperature gradients:
    J_soret = -ρ * D * S_T * ∇T
    
    where S_T is the Soret coefficient [1/K]
    
    For moisture transport:
    ∂Cl/∂t = ∇·(D_T * Cl * ∇T / T) 
    
    where D_T is thermal diffusion coefficient and the Cl/T term
    accounts for concentration dependence.
    
    Parameters:
    -----------
    Cl_vals : array
        Liquid moisture content [kg/m³]
    T_vals : array
        Temperature [°C]
    dx : float
        Grid spacing [m]
    
    Returns:
    --------
    soret_source : array
        Soret effect contribution to moisture equation [kg/(m³·s)]
    """
    n_nodes = len(Cl_vals)
    soret_source = nm.zeros_like(Cl_vals)
    
    # Convert to Kelvin for proper thermodynamic calculations
    T_K = T_vals + 273.15
    
    # Temperature-dependent Soret coefficient
    # Typically negative for water (moves from hot to cold)
    S_T = -1e-3 * (1 + 0.001 * T_vals)  # [1/K]
    
    # Moisture-dependent diffusivity enhancement
    D_T_eff = D_t * (1 + 0.5 * Cl_vals / nm.max(Cl_vals))
    
    for i in range(n_nodes):
        if i == 0 and n_nodes > 1:
            # Forward difference
            dT_dx = (T_K[1] - T_K[0]) / dx
            Cl_avg = 0.5 * (Cl_vals[0] + Cl_vals[1])
            T_avg = 0.5 * (T_K[0] + T_K[1])
            D_avg = 0.5 * (D_T_eff[0] + D_T_eff[1])
            S_avg = 0.5 * (S_T[0] + S_T[1])
            
            # Soret flux = D_T * S_T * Cl * ∇T / T
            soret_flux = D_avg * S_avg * Cl_avg * dT_dx / T_avg
            soret_source[i] = -soret_flux / dx
            
        elif i == n_nodes - 1 and n_nodes > 1:
            # Backward difference
            dT_dx = (T_K[i] - T_K[i-1]) / dx
            Cl_avg = 0.5 * (Cl_vals[i-1] + Cl_vals[i])
            T_avg = 0.5 * (T_K[i-1] + T_K[i])
            D_avg = 0.5 * (D_T_eff[i-1] + D_T_eff[i])
            S_avg = 0.5 * (S_T[i-1] + S_T[i])
            
            soret_flux = D_avg * S_avg * Cl_avg * dT_dx / T_avg
            soret_source[i] = soret_flux / dx
            
        elif n_nodes > 2:
            # Central differences with flux form
            # Left face
            dT_dx_left = (T_K[i] - T_K[i-1]) / dx
            Cl_left = 0.5 * (Cl_vals[i-1] + Cl_vals[i])
            T_left = 0.5 * (T_K[i-1] + T_K[i])
            D_left = 0.5 * (D_T_eff[i-1] + D_T_eff[i])
            S_left = 0.5 * (S_T[i-1] + S_T[i])
            flux_left = D_left * S_left * Cl_left * dT_dx_left / T_left
            
            # Right face  
            dT_dx_right = (T_K[i+1] - T_K[i]) / dx
            Cl_right = 0.5 * (Cl_vals[i] + Cl_vals[i+1])
            T_right = 0.5 * (T_K[i] + T_K[i+1])
            D_right = 0.5 * (D_T_eff[i] + D_T_eff[i+1])
            S_right = 0.5 * (S_T[i] + S_T[i+1])
            flux_right = D_right * S_right * Cl_right * dT_dx_right / T_right
            
            # Divergence of flux
            soret_source[i] = (flux_right - flux_left) / dx
    
    return soret_source

def calculate_dufour_effect(Cl_vals, T_vals, dx):
    """
    Calculate Dufour effect - heat flux due to mass gradients.
    
    This is the reciprocal of Soret effect (Onsager reciprocal relations).
    
    Physical basis:
    J_q_dufour = -L_TD * ∇(μ_m/T) ≈ -D_TD(T) * ∇Cl / T
    
    where:
    - L_TD is the cross-coupling coefficient (Onsager)
    - μ_m is the chemical potential of moisture
    - D_TD is the temperature-dependent Dufour coefficient
    
    Parameters:
    -----------
    Cl_vals : array
        Liquid moisture content [kg/m³]
    T_vals : array  
        Temperature [°C]
    dx : float
        Grid spacing [m]
        
    Returns:
    --------
    dufour_source : array
        Dufour effect contribution to heat equation [W/m³]
    """
    n_nodes = len(Cl_vals)
    dufour_source = nm.zeros_like(Cl_vals)
    
    # Convert to Kelvin for proper thermodynamic treatment
    T_K = T_vals + 273.15
    
    # Temperature-dependent Dufour coefficient
    # Related to Soret coefficient by Onsager reciprocal relations: L_TD = L_DT
    # D_TD(T) = k_c * T * f(T) where f(T) accounts for temperature dependence
    T_ref = 273.15  # Reference temperature [K]
    D_TD_base = k_c  # Base heat-mass coupling coefficient [W·s/(m·K·kg)]
    
    # Temperature-dependent enhancement factor
    # Physical basis: stronger coupling at higher temperatures due to increased molecular motion
    temp_factor = (T_K / T_ref) * nm.exp(-500.0 / T_K)  # Arrhenius-type dependence
    D_TD_eff = D_TD_base * temp_factor
    
    for i in range(n_nodes):
        if i == 0 and n_nodes > 1:
            # Forward difference
            dCl_dx = (Cl_vals[1] - Cl_vals[0]) / dx
            # Average properties at the interface
            T_avg = 0.5 * (T_K[0] + T_K[1])
            D_avg = 0.5 * (D_TD_eff[0] + D_TD_eff[1])
            # Dufour heat flux = D_TD * ∇Cl / T (thermodynamically consistent)
            dufour_source[i] = D_avg * dCl_dx / T_avg
            
        elif i == n_nodes - 1 and n_nodes > 1:
            # Backward difference
            dCl_dx = (Cl_vals[i] - Cl_vals[i-1]) / dx  
            T_avg = 0.5 * (T_K[i-1] + T_K[i])
            D_avg = 0.5 * (D_TD_eff[i-1] + D_TD_eff[i])
            dufour_source[i] = D_avg * dCl_dx / T_avg
            
        elif n_nodes > 2:
            # Central differences with proper averaging
            dCl_dx = (Cl_vals[i+1] - Cl_vals[i-1]) / (2 * dx)
            # Use local temperature and coefficient
            dufour_source[i] = D_TD_eff[i] * dCl_dx / T_K[i]
    
    return dufour_source

###############################################################################
#                                MAIN SECTION                                 #
###############################################################################

# Read data and calculate boundary conditions (same as before)
airTemp, airVel, prec, gloSolIr, rh = read_input_data()

# Calculation of equivalent sol-air temperature (in °C)
h = 22.7; alpha = 0.8; T_cor_fact = 4.0
t_o = [alpha * glob_solir / h + air_temp - T_cor_fact for 
       glob_solir, air_temp in zip(gloSolIr, airTemp)]

# Calculation of convective heat transfer coefficient 
h_o = [6.0 + 4.0*vel if vel <= 5.0 else 7.41*(vel**0.78) for vel in airVel]

# Boundary conditions (same as before)
t_i = 0.0; h_i = 99.75
cl_i = 0.01; cl_o = 0.02; h_l_i = 1e-6; h_l_o = 2e-6
cv_i = 0.005; cv_o = 0.008; h_v_i = 1e-5; #h_v_o = 2e-5
h_v_o = calculate_vapor_bc_coefficient(airTemp, rh)

# Mesh is generated at runtime by mesh_hook function - allows parameterized 
# 1D domain
# (length = d_ins, resolution = dx) without needing external mesh files
filename_mesh = UserMeshIO(mesh_hook)

# Enhanced materials with new transport coefficients
materials = {
    'mat': ({'lam': lam_i, 'rho_cp': rho_wet * c_wet},),
    'liquid_mat': ({'D_l': D_l},),    # Liquid diffusion
    'vapor_mat': ({'D_v': D_v},),     # Vapor diffusion
    'pressure_mat': ({'D_p': D_p},),  # Pressure diffusion
    'thermal_mat': ({'D_t': D_t},),   # Thermal diffusion (Soret)
    'coupling_mat': ({'k_c': k_c},),  # Heat-moisture coupling
    
    # Dynamic boundary conditions
    'h_out_dyn': 'get_h_o',
    'h_v_out_dyn': 'get_h_v_o',
    'T_out_dyn': 'get_t_o', 
    'cl_out_dyn': 'get_cl_o',  
    'cv_out_dyn': 'get_cv_o',  
    'cv_in_dyn': 'get_cv_i',   
    
    # Fixed boundary conditions
    'in_fixed': ({'h_in': h_i, 'T_in': t_i, 'h_l_in': h_l_i, 
                  'Cl_in': cl_i, 'h_v_in': h_v_i, 'Cv_in': cv_i},),
    'out_fixed': ({'h_l_out': h_l_o
                   #, 'h_v_out': h_v_o
                   },),
    
    # Enhanced transport terms
    'pressure_gradient': 'get_pressure_gradient_source',          
    'soret_effect': 'get_soret_effect',                           
    'heat_moisture_coupling': 'get_heat_moisture_coupling',       
    'shrinkage_diffusivity': 'get_shrinkage_modified_diffusivity',
    
    # Existing terms
    'evaporation_rate': 'get_evaporation_rate',
}

functions = {
    'get_h_o': (get_h_o,),
    'get_h_v_o': (get_h_v_o,),  # Dynamic vapor convection coeff based on RH
    'get_t_o': (get_t_o,),
    'get_cl_o': (get_cl_o,),      # Added missing liquid boundary function
    'get_cv_o': (get_cv_o,),      # Added vapor boundary function  
    'get_cv_i': (get_cv_i,),      # Added inner vapor boundary function
    'get_evaporation_rate': (get_evaporation_rate,),
    
    # NEW FUNCTIONS
    'get_pressure_gradient_source': (get_pressure_gradient_source,),
    'get_soret_effect': (get_soret_effect,),
    'get_heat_moisture_coupling': (get_heat_moisture_coupling,),
    'get_shrinkage_modified_diffusivity': (get_shrink_mod_diff,),
}

regions = {
    'Omega': 'all',
    'Gamma_Left': ('vertices in (x < 0.00001)', 'facet'),   # Outer boundary
    'Gamma_Right': ('vertices in (x > %f)' % (d_ins - 0.00001), 'facet'), # Inner
}

fields = {
    'temperature': ('real', 1, 'Omega', 1),
    'liquid': ('real', 1, 'Omega', 1),
    'vapor': ('real', 1, 'Omega', 1),
}

variables = {
    'T': ('unknown field', 'temperature', 0, 1),
    's': ('test field', 'temperature', 'T'),
    'Cl': ('unknown field', 'liquid', 1, 1),
    'r': ('test field', 'liquid', 'Cl'),
    'Cv': ('unknown field', 'vapor', 2, 1),
    'w': ('test field', 'vapor', 'Cv'),
}

integrals = {'i': 1,}

# ENHANCED EQUATIONS with additional physics terms
equations = {
    'Heat': """dw_dot.i.Omega(mat.rho_cp, s, dT/dt)
             + dw_laplace.i.Omega(mat.lam, s, T)
             = - dw_bc_newton.i.Gamma_Left(h_out_dyn.val, T_out_dyn.val, s, T)
               - dw_bc_newton.i.Gamma_Right(in_fixed.h_in, in_fixed.T_in, s, T)
               + dw_volume_lvf.i.Omega(heat_moisture_coupling.val, s)
             """,
             
    'Liquid': """dw_dot.i.Omega(r, dCl/dt)
               + dw_laplace.i.Omega(shrinkage_diffusivity.val, r, Cl)
               = - dw_volume_lvf.i.Omega(evaporation_rate.val, r)
                 - dw_bc_newton.i.Gamma_Left(out_fixed.h_l_out, cl_out_dyn.val, r, Cl)
                 - dw_bc_newton.i.Gamma_Right(in_fixed.h_l_in, in_fixed.Cl_in, r, Cl)
                 + dw_volume_lvf.i.Omega(pressure_gradient.val, r)
                 + dw_volume_lvf.i.Omega(soret_effect.val, r)
               """,
               
    'Vapor': """dw_dot.i.Omega(w, dCv/dt)
              + dw_laplace.i.Omega(vapor_mat.D_v, w, Cv)
              = + dw_volume_lvf.i.Omega(evaporation_rate.val, w)
                - dw_bc_newton.i.Gamma_Left(h_v_out_dyn.val, cv_out_dyn.val, w, Cv)
                - dw_bc_newton.i.Gamma_Right(in_fixed.h_v_in, in_fixed.Cv_in, w, Cv)
              """
}

# Boundary conditions (empty but required)
ebcs = {}

# Initial conditions
ics = {
    'ic_T': ('Omega', {'T.0': 0.0}),
    'ic_Cl': ('Omega', {'Cl.0': 0.015}),
    'ic_Cv': ('Omega', {'Cv.0': 0.006}),
}

# Time and solver parameters
nr_of_hours = 2
stop = nr_of_hours * 3600.0
dt = 10.0
nr_of_steps = int(stop/dt)


solvers = {
    'ls': ('ls.scipy_direct', {}),
    'newton': ('nls.newton', {
        'i_max': 25,      # Increased for enhanced nonlinear system
        'eps_a': 1e-8,
        'eps_r': 1e-6,
        'is_linear': False,
    }),
    'ts': ('ts.simple', {
        't0': 0.0,
        't1': stop,
        'dt': dt,
        'n_step': nr_of_steps,
        'verbose': 1,
    }),
}

def save_enhanced_results(out, problem, state, extend=False):
    """Enhanced results saving with additional physics terms."""
    import os
    
    filename = os.path.join(problem.conf.options['output_dir'], "enhanced_moisture_results_v6.csv")
    header = [
        "Time (s)", "Node Index", "Position (m)", 
        "Temperature (°C)", "Liquid Content (kg/m³)", "Vapor Content (kg/m³)",
        "Capillary Pressure (Pa)", "Osmotic Pressure (Pa)", 
        "Pressure Gradient Effect", "Soret Effect", "Dufour Effect", "Shrinkage Factor"
    ]
    
    # Get state
    state_vec = state.get_state()
    coors = problem.fields['temperature'].get_coor()
    n_dof = problem.fields['temperature'].n_nod
    
    if len(state_vec) >= 3 * n_dof:
        T_vals = state_vec[:n_dof]
        Cl_vals = state_vec[n_dof:2*n_dof]
        Cv_vals = state_vec[2*n_dof:3*n_dof]
    else:
        T_vals = state_vec[:n_dof]
        Cl_vals = state_vec[n_dof:] if len(state_vec) > n_dof else nm.zeros(n_dof)
        Cv_vals = nm.zeros(n_dof)
    
    # Calculate enhanced physics terms
    p_cap = calculate_capillary_pressure(Cl_vals, T_vals)
    p_osm = calculate_osmotic_pressure(Cl_vals, T_vals)
    shrinkage_factors = nm.maximum(0.1, 1.0 - alpha_shrink * Cl_vals)
    pressure_flux = calculate_pressure_driven_flux(Cl_vals, T_vals, dx)
    soret_source = calculate_soret_effect_proper(Cl_vals, T_vals, dx)
    dufour_source = calculate_dufour_effect(Cl_vals, T_vals, dx)  
    
    # Write results every hour
    if problem.ts.time % 3600 < problem.ts.dt:
        file_exists = os.path.exists(filename) and os.path.getsize(filename) > 0
        if not file_exists:
            with open(filename, 'w', newline='') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow(header)
        
        with open(filename, 'a', newline='') as csvfile:
            writer = csv.writer(csvfile)
            for i, coord in enumerate(coors):
                writer.writerow([
                    problem.ts.time, i, coord[0], T_vals[i], Cl_vals[i], Cv_vals[i],
                    p_cap[i], p_osm[i], 
                    pressure_flux[i],
                    soret_source[i],
                    dufour_source[i],
                    shrinkage_factors[i]
                ])
        
        print(f"Hour {int(problem.ts.time/3600)}: Enhanced physics simulation")
    
    return out

options = {
    'nls': 'newton',
    'ls': 'ls', 
    'ts': 'ts',
    'save_times': [3600*i for i in range(1, nr_of_hours + 1)],
    'post_process_hook': save_enhanced_results,
    'output_dir': './output_snow_step6',
}