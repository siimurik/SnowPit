r"""
==================================================================
Enhanced Heat-Liquid-Vapor Coupled Transport System with 
Additional Moisture Movement Terms Based on Clay Research
------------------------------------------------------------------
Enhanced system includes:
1. Pressure-driven liquid flow (Darcy's law)
2. Osmotic/capillary pressure effects
3. Soret effect (thermal diffusion)
4. Cross-coupling between temperature and moisture
5. Shrinkage effects on transport properties

Based on equations from moisture movement research:
- Pressure diffusion: J_m = -D_p * ∇p
- Thermal diffusion: J_m = -D_t * ∇T  (Soret effect)
- Heat-moisture coupling: J_h = -k_c * ∇C

==================================================================
"""
from __future__ import absolute_import
import numpy as nm
import math
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
    x = Tc / T_K * (C1*teta + C2*teta**1.5 + C3*teta**3 + C4*teta**3.5 + C5*teta**4 + C6*teta**7.5)
    x = math.exp(x) * Pc
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

def get_shrinkage_modified_diffusivity(ts, coors, mode=None, equations=None, term=None,
                                     problem=None, **kwargs):
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

# Keep all existing functions (read_input_data, get_h_o, etc.) from original code
def read_input_data(filename="DATA.csv"):
    """Read air temperature, air velocity, precipitation, global solar irradation and relative humidity data from CSV file."""
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

# Read data and calculate boundary conditions (same as before)
airTemp, airVel, prec, gloSolIr, rh = read_input_data()

h = 22.7; alpha = 0.8; T_cor_fact = 4.0
t_o = [alpha * glob_solir / h + air_temp - T_cor_fact for glob_solir, air_temp in zip(gloSolIr, airTemp)]
h_o = [6.0 + 4.0*vel if vel <= 5.0 else 7.41*(vel**0.78) for vel in airVel]

# Boundary conditions (same as before)
t_i = 0.0; h_i = 99.75
cl_i = 0.01; cl_o = 0.02; h_l_i = 1e-6; h_l_o = 2e-6
cv_i = 0.005; cv_o = 0.008; h_v_i = 1e-5; h_v_o = 2e-5

# Keep existing boundary condition functions (get_h_o, get_t_o, etc.)
def get_h_o(ts, coors, mode=None, **kwargs):
    if mode != 'qp' or coors is None: return {}
    hour_idx = min(int(ts.time / 3600), len(h_o) - 1)
    val = nm.full((coors.shape[0], 1, 1), h_o[hour_idx], dtype=nm.float64)
    return {'val': val}

def get_t_o(ts, coors, mode=None, **kwargs):
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

def get_cv_o(ts, coors, mode=None, **kwargs):
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

def get_evaporation_rate(ts, coors, mode=None, equations=None, term=None, problem=None, **kwargs):
    """Enhanced evaporation rate (same as before)"""
    if mode != 'qp' or coors is None: return {}
    # ... (same implementation as before)
    hour_idx = min(int(ts.time / 3600), len(rh) - 1)
    current_rh = rh[hour_idx] if hour_idx < len(rh) else 70.0
    humidity_factor = max(0, 1.0 - current_rh/100.0)
    evap_rate = nm.full(coors.shape[0], k_evap * humidity_factor * 0.5)
    val = evap_rate.reshape((coors.shape[0], 1, 1))
    return {'val': val}

def mesh_hook(mesh, mode):
    """Generate the 1D mesh."""
    if mode == 'read':
        coors = nm.linspace(0.0, d_ins, nodes).reshape((nodes, 1))
        conn = nm.arange(nodes, dtype=nm.int32).repeat(2)[1:-1].reshape((-1, 2))
        mat_ids = nm.zeros(nodes - 1, dtype=nm.int32)
        descs = ['1_2']
        mesh = Mesh.from_data('heat_liquid_1d', coors, None, [conn], [mat_ids], descs)
        return mesh
    elif mode == 'write':
        pass

filename_mesh = UserMeshIO(mesh_hook)

# Enhanced materials with new transport coefficients
materials = {
    'mat': ({'lam': lam_i, 'rho_cp': rho_wet * c_wet},),
    'liquid_mat': ({'D_l': D_l},),
    'vapor_mat': ({'D_v': D_v},),
    'pressure_mat': ({'D_p': D_p},),  # NEW: Pressure diffusion
    'thermal_mat': ({'D_t': D_t},),   # NEW: Thermal diffusion (Soret)
    'coupling_mat': ({'k_c': k_c},),  # NEW: Heat-moisture coupling
    
    # Dynamic boundary conditions
    'h_out_dyn': 'get_h_o',
    'T_out_dyn': 'get_t_o', 
    'cl_out_dyn': 'get_cl_o',  # Added missing liquid boundary function
    'cv_out_dyn': 'get_cv_o',  # Added vapor boundary function
    'cv_in_dyn': 'get_cv_i',   # Added inner vapor boundary function
    
    # Fixed boundary conditions
    'in_fixed': ({'h_in': h_i, 'T_in': t_i, 'h_l_in': h_l_i, 'Cl_in': cl_i, 'h_v_in': h_v_i, 'Cv_in': cv_i},),
    'out_fixed': ({'h_l_out': h_l_o, 'h_v_out': h_v_o},),
    
    # Enhanced transport terms
    'pressure_gradient': 'get_pressure_gradient_source',      # NEW
    'soret_effect': 'get_soret_effect',                      # NEW
    'heat_moisture_coupling': 'get_heat_moisture_coupling',   # NEW
    'shrinkage_diffusivity': 'get_shrinkage_modified_diffusivity',  # NEW
    
    # Existing terms
    'evaporation_rate': 'get_evaporation_rate',
}

functions = {
    'get_h_o': (get_h_o,),
    'get_t_o': (get_t_o,),
    'get_cl_o': (get_cl_o,),      # Added missing liquid boundary function
    'get_cv_o': (get_cv_o,),      # Added vapor boundary function  
    'get_cv_i': (get_cv_i,),      # Added inner vapor boundary function
    'get_evaporation_rate': (get_evaporation_rate,),
    
    # NEW FUNCTIONS
    'get_pressure_gradient_source': (get_pressure_gradient_source,),
    'get_soret_effect': (get_soret_effect,),
    'get_heat_moisture_coupling': (get_heat_moisture_coupling,),
    'get_shrinkage_modified_diffusivity': (get_shrinkage_modified_diffusivity,),
}

regions = {
    'Omega': 'all',
    'Gamma_Left': ('vertices in (x < 0.00001)', 'facet'),
    'Gamma_Right': ('vertices in (x > %f)' % (d_ins - 0.00001), 'facet'), 
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
                - dw_bc_newton.i.Gamma_Left(out_fixed.h_v_out, cv_out_dyn.val, w, Cv)
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
        'n_step': int(stop/dt),
        'verbose': 1,
    }),
}

def save_enhanced_results(out, problem, state, extend=False):
    """Enhanced results saving with additional physics terms."""
    import os
    
    filename = os.path.join(problem.conf.options['output_dir'], "enhanced_moisture_results.csv")
    header = [
        "Time (s)", "Node Index", "Position (m)", 
        "Temperature (°C)", "Liquid Content (kg/m³)", "Vapor Content (kg/m³)",
        "Capillary Pressure (Pa)", "Osmotic Pressure (Pa)", 
        "Pressure Gradient Effect", "Soret Effect", "Shrinkage Factor"
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
                    p_cap[i], p_osm[i], 0.0, 0.0, shrinkage_factors[i]  # Simplified output
                ])
        
        print(f"Hour {int(problem.ts.time/3600)}: Enhanced physics simulation")
    
    return out

options = {
    'nls': 'newton',
    'ls': 'ls', 
    'ts': 'ts',
    'save_times': [3600*i for i in range(1, nr_of_hours + 1)],
    'post_process_hook': save_enhanced_results,
    'output_dir': './output_enhanced_moisture',
}