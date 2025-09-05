r"""
==================================================================
Enhanced Heat-Liquid-Vapor-Pressure Coupled Transport System with 
Additional Moisture Movement Terms for Snow
------------------------------------------------------------------
Four-field coupled system with enhanced physics including 
realistic pressure
------------------------------------------------------------------
GOVERNING EQUATIONS:
Heat Equation (Temperature T), Liquid Equation (Liquid Content Cl),
Vapor Equation (Vapor Content Cv), Pressure Equation (Pore Pressure P).

    Heat:     ρcp ∂T/∂t = λ ∇²T + k_c ∇²Cl + ΔH εᴸ ∂Cl/∂t + k_p ∇²P
    Liquid:   ∂Cl/∂t = D_eff(Cl) ∇²Cl + D_t ∇²T + D_p ∇²P - m_evap
    Vapor:    ∂Cv/∂t = D_v ∇²Cv + D_t_v ∇²T + D_p_v ∇²P + m_evap
    Pressure: ∂P/∂t = D_p_eff ∇²P + α_T ∇²T + α_L ∇²Cl + α_V ∇²Cv

Physical meaning of the last term: 
  Pressure changes due to thermal expansion, liquid filling pores 
  and vapor pressure effects.
------------------------------------------------------------------
BOUNDARY CONDITIONS:

Left Boundary (x = 0, Outer Surface):
    Heat:        -λ ∂T/∂x = h_o(t) [T - T_o(t)]
    Liquid:      -D_eff ∂Cl/∂x = h_l_o [Cl - Cl_o(t)]  
    Vapor:       -D_v ∂Cv/∂x = h_v_o [Cv - Cv_o(t)]
    Pressure:    -k_p ∂P/∂x = h_p_o [P - P_atm] (atmospheric pressure BC)

Right Boundary (x = d_ins, Inner Surface):
    Heat:        -λ ∂T/∂x = h_i [T - T_i]
    Liquid:      -D_eff ∂Cl/∂x = h_l_i [Cl - Cl_i]
    Vapor:       -D_v ∂Cv/∂x = h_v_i [Cv - Cv_i]
    Pressure:    -D_p ∂P/∂x = h_p_i [P - P_i]
==================================================================
"""

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

# Pressure and suction parameters
k_s = 1e-10       # Permeability [m^2]
phi = 0.4         # Snow porosity
mu_water = 1e-3   # Dynamic viscosity of water [Pa*s]
p_atm = 101325    # Atmospheric pressure [Pa]
gamma_surface = 0.075  # Surface tension [N/m]
contact_angle = 5      # Contact angle [degrees]

# Enhanced transport coefficients
D_l = 1e-6        # Liquid diffusion coefficient [m^2/s]
D_v = 2e-5        # Vapor diffusion coefficient [m^2/s]
D_p = k_s / (mu_water * phi) # Pressure diffusion coefficient [m^2/(Pa*s)]
D_t = 1e-10       # Thermal diffusion coefficient [m^2/(K*s)] - Soret effect
k_c = 0.1         # Heat-moisture coupling coefficient [W*s/(m*K*kg)]
k_p = 1e-9        # Pressure-heat-moisture coefficient

# NEW: Pressure field coefficients
D_p_eff = 1e-8    # Effective pressure diffusion [m^2/s]
alpha_T = 1e3     # Thermal expansion coefficient for pressure [Pa/K]
alpha_L = 5e4     # Liquid-pressure coupling [Pa*m^3/kg]
alpha_V = 2e4     # Vapor-pressure coupling [Pa*m^3/kg]
h_p_o = 1e-6      # Pressure transfer coefficient at outer boundary [m/s]
h_p_i = 5e-7      # Base pressure transfer coefficient [m/s]
P_base = p_atm + 1000  # Base pressure (slightly above atmospheric) [Pa]

# Pressure-moisture coupling coefficients
D_t_v = 5e-11     # Thermal diffusion for vapor [m^2/(K*s)]
D_p_v = 1e-11     # Pressure diffusion for vapor [m^2/(Pa*s)]

# Additional physical constants
L_v = 2.45e6      # Latent heat of vaporization [J/kg]
L_s = 2.83e6      # Latent heat of sublimation [J/kg]
R_v = 461.5       # Specific gas constant for water vapor [J/(kg*K)]

# Enhanced evaporation model parameters
k_evap = 1e-7     # Evaporation rate constant [1/s]
T_ref = 0.0       # Reference temperature [°C]
Cl_min = 0.001    # Minimum liquid content for evaporation [kg/m³]
T_evap_max = 10.0 # Maximum temperature for linear evaporation model [°C]

def Psat_WV(T_K):
    """Water vapour saturation pressure"""
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

def get_p_atm(ts, coors, mode=None, **kwargs):
    """Atmospheric pressure boundary condition"""
    if mode != 'qp' or coors is None:
        return {}
    val = nm.full((coors.shape[0], 1, 1), p_atm, dtype=nm.float64)
    return {'val': val}

def get_h_p_o(ts, coors, mode=None, **kwargs):
    """Pressure transfer coefficient at outer boundary"""
    if mode != 'qp' or coors is None:
        return {}
    val = nm.full((coors.shape[0], 1, 1), h_p_o, dtype=nm.float64)
    return {'val': val}

def get_heat_moist_coup(ts, coors, mode=None, equations=None, term=None,
                             problem=None, **kwargs):
    """Heat-moisture coupling term: k_c * ∇C"""
    if mode != 'qp' or coors is None:
        return {}
    
    variables = problem.get_variables()
    Cl_var = variables['Cl']
    Cl_vals = Cl_var.get_state_in_region(problem.domain.regions['Omega'])
    Cl_vals_flat = Cl_vals.flatten()
    
    Cl_grad = nm.gradient(Cl_vals_flat) / dx
    coupling_effect = k_c * Cl_grad / (rho_wet * c_wet)
    
    val = coupling_effect.reshape((coupling_effect.shape[0], 1, 1))
    return {'val': val}

def read_input_data(filename="DATA.csv"):
    """Read environmental data from CSV file"""
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

def get_h_o(ts, coors, mode=None, **kwargs):
    """Time-dependent heat transfer coefficient"""
    if mode != 'qp' or coors is None: return {}
    hour_idx = min(int(ts.time / 3600), len(h_o) - 1)
    val = nm.full((coors.shape[0], 1, 1), h_o[hour_idx], dtype=nm.float64)
    return {'val': val}

def get_h_v_o(ts, coors, mode=None, **kwargs):
    """Time-dependent vapor heat transfer coefficient"""
    if mode != 'qp' or coors is None: return {}
    hour_idx = min(int(ts.time / 3600), len(h_v_o) - 1)
    val = nm.full((coors.shape[0], 1, 1), h_v_o[hour_idx], dtype=nm.float64)
    return {'val': val}

def get_t_o(ts, coors, mode=None, **kwargs):
    """Time-dependent sol-air temperature"""
    if mode != 'qp' or coors is None: return {}
    hour_idx = min(int(ts.time / 3600), len(t_o) - 1)
    val = nm.full((coors.shape[0], 1, 1), t_o[hour_idx], dtype=nm.float64)
    return {'val': val}

def get_cl_o(ts, coors, mode=None, **kwargs):
    """Time-dependent outer liquid reference content"""
    if mode != 'qp' or coors is None: return {}
    hour_idx = min(int(ts.time / 3600), len(rh) - 1)
    current_rh = rh[hour_idx] if hour_idx < len(rh) else 70.0
    humidity_factor = current_rh / 100.0
    cl_adjusted = cl_o * (0.5 + 0.5 * humidity_factor)
    val = nm.full((coors.shape[0], 1, 1), cl_adjusted, dtype=nm.float64)
    return {'val': val}

def get_vapor_content_from_RH(T_celsius, RH_percent):
    """Calculate vapor content from temperature and relative humidity"""
    T_K = T_celsius + 273.15
    P_sat = Psat_WV(T_K)  # Pa
    P_vapor = P_sat * RH_percent / 100.0  # Pa
    M_water = 0.018015  # kg/mol
    R_gas = 8.314       # J/(mol·K)
    vapor_content = (P_vapor * M_water) / (R_gas * T_K)  # kg/m³
    return vapor_content

def get_cv_o(ts, coors, mode=None, **kwargs):
    """Time-dependent outer vapor reference content"""
    if mode != 'qp' or coors is None:
        return {}
    
    hour_idx = min(int(ts.time / 3600), len(rh) - 1)
    current_rh = rh[hour_idx] if hour_idx < len(rh) else 70.0
    current_t = t_o[hour_idx] if hour_idx < len(t_o) else -5.0
    
    try:
        cv_physical = get_vapor_content_from_RH(current_t, current_rh)
        cv_min = 1e-6
        cv_max = 0.05
        cv_adjusted = nm.clip(cv_physical, cv_min, cv_max)
    except Exception as e:
        print(f"Warning: Vapor content calculation failed: {e}")
        cv_adjusted = 0.005 * (current_rh / 100.0)
    
    val = nm.full((coors.shape[0], 1, 1), cv_adjusted, dtype=nm.float64)
    return {'val': val}

def get_cv_i(ts, coors, mode=None, **kwargs):
    """Time-dependent inner vapor reference content"""
    if mode != 'qp' or coors is None: return {}
    val = nm.full((coors.shape[0], 1, 1), cv_i, dtype=nm.float64)
    return {'val': val}

def calculate_vapor_bc_coefficient(T_celsius, RH_percent):
    """Calculate vapor boundary condition coefficient"""
    h_base = 2e-5
    T_factor = [1 + 0.01 * T for T in T_celsius]
    RH_factor = [1 + 0.5 * (RH / 100.0) for RH in RH_percent]
    h_v_eff = [h_base * T_f * RH_f for T_f, RH_f in zip(T_factor, RH_factor)]
    return h_v_eff

def get_evaporation_rate(ts, coors, mode=None, equations=None, term=None, 
                         problem=None, **kwargs):
    """Stabilized evaporation rate"""
    if mode != 'qp' or coors is None:
        return {}

    hour_idx = min(int(ts.time / 3600), len(rh) - 1)
    current_rh = rh[hour_idx]
    
    variables = problem.get_variables()
    T_vals = variables['T'].get_state_in_region(problem.domain.regions['Omega']).flatten()
    Cl_vals = variables['Cl'].get_state_in_region(problem.domain.regions['Omega']).flatten()
    Cv_vals = variables['Cv'].get_state_in_region(problem.domain.regions['Omega']).flatten()

    humidity_factor = nm.tanh(2 * (1.0 - current_rh / 100.0))
    T_factor = 0.5 * (1 + nm.tanh((T_vals - T_ref) / 5.0))
    Cl_factor = nm.tanh(5 * nm.maximum(0.0, Cl_vals - Cl_min))
    
    T_K = T_vals + 273.15
    Psat = nm.maximum(100.0, Psat_WV(T_K) * 100.0)
    M_wv = 0.018
    R_gas = 8.314
    rho_v_sat = Psat * M_wv / (R_gas * T_K)
    vapor_factor = nm.tanh(2 * nm.maximum(0.0, 1.0 - Cv_vals / rho_v_sat))

    evap_rate = k_evap * 0.1 * humidity_factor * T_factor * Cl_factor * vapor_factor
    
    val = evap_rate.reshape((evap_rate.shape[0], 1, 1))
    return {'val': val}

def get_prec_enhanced_cl_ref(ts, coors, mode=None, **kwargs):
    """Precipitation-enhanced reference liquid concentration"""
    if mode != 'qp' or coors is None:
        return {}
    
    hour_idx = min(int(ts.time / 3600), len(prec) - 1)
    current_prec = prec[hour_idx] if hour_idx < len(prec) else 0.0
    current_temp = t_o[hour_idx] if hour_idx < len(t_o) else -5.0
    current_rh = rh[hour_idx] if hour_idx < len(rh) else 70.0
    
    base_cl_ref = 0.02 * (current_rh / 100.0)
    
    if current_prec > 0.1:
        if current_temp >= 2.0:
            rain_intensity_factor = min(2.0, current_prec / 5.0)
            prec_cl_ref = 0.15 * rain_intensity_factor
        elif current_temp <= -2.0:
            snow_intensity_factor = min(1.5, current_prec / 3.0)
            prec_cl_ref = 0.05 * snow_intensity_factor
        else:
            rain_fraction = (current_temp + 2.0) / 4.0
            rain_contribution = 0.15 * min(2.0, current_prec / 5.0) * rain_fraction
            snow_contribution = 0.05 * min(1.5, current_prec / 3.0) * (1.0 - rain_fraction)
            prec_cl_ref = rain_contribution + snow_contribution
        
        enhanced_cl_ref = max(base_cl_ref, prec_cl_ref)
        enhanced_cl_ref = min(enhanced_cl_ref, 0.3)
    else:
        enhanced_cl_ref = base_cl_ref
    
    val = nm.full((coors.shape[0], 1, 1), enhanced_cl_ref, dtype=nm.float64)
    return {'val': val}

def get_prec_enhanced_h_l(ts, coors, mode=None, **kwargs):
    """Precipitation-enhanced mass transfer coefficient"""
    if mode != 'qp' or coors is None:
        return {}
    
    hour_idx = min(int(ts.time / 3600), len(prec) - 1)
    current_prec = prec[hour_idx] if hour_idx < len(prec) else 0.0
    current_temp = t_o[hour_idx] if hour_idx < len(t_o) else -5.0
    current_vel = airVel[hour_idx] if hour_idx < len(airVel) else 2.0
    current_rh = rh[hour_idx] if hour_idx < len(rh) else 70.0
    
    base_h_l = 2e-6
    temp_factor = nm.exp(0.03 * current_temp)
    wind_factor = 1.0 + 0.5 * current_vel
    humidity_factor = 1.0 + 0.3 * (current_rh / 100.0)
    
    if current_prec > 0.1:
        if current_temp >= 2.0:
            rain_enhancement = 2.0 + 3.0 * min(1.0, current_prec / 10.0)
        elif current_temp <= -2.0:
            snow_enhancement = 0.7 + 1.5 * min(1.0, current_prec / 5.0)
        else:
            rain_fraction = (current_temp + 2.0) / 4.0
            rain_enh = 2.0 + 3.0 * min(1.0, current_prec / 10.0)
            snow_enh = 0.7 + 1.5 * min(1.0, current_prec / 5.0)
            rain_enhancement = rain_fraction * rain_enh + (1.0 - rain_fraction) * snow_enh
        prec_factor = rain_enhancement
    else:
        prec_factor = 1.0
    
    enhanced_h_l = base_h_l * temp_factor * wind_factor * humidity_factor * prec_factor
    enhanced_h_l = max(enhanced_h_l, 1e-8)
    enhanced_h_l = min(enhanced_h_l, 1e-3)
    
    val = nm.full((coors.shape[0], 1, 1), enhanced_h_l, dtype=nm.float64)
    return {'val': val}

def get_accumulated_snow_source(ts, coors, mode=None, equations=None, 
                               term=None, problem=None, **kwargs):
    """Snow accumulation and melting source term"""
    if mode != 'qp' or coors is None:
        return {}
    
    hour_idx = min(int(ts.time / 3600), len(prec) - 1)
    current_prec = prec[hour_idx] if hour_idx < len(prec) else 0.0
    current_temp = t_o[hour_idx] if hour_idx < len(t_o) else -5.0
    
    variables = problem.get_variables()
    T_vals = variables['T'].get_state_in_region(problem.domain.regions['Omega']).flatten()
    Cl_vals = variables['Cl'].get_state_in_region(problem.domain.regions['Omega']).flatten()
    
    n_nodes = len(T_vals)
    snow_source = nm.zeros(n_nodes)
    
    snow_penetration_depth = 0.02
    n_surface_nodes = max(1, int(snow_penetration_depth / dx))
    
    # Snow accumulation
    if current_prec > 0.1 and current_temp <= 1.0:
        snow_rate_flux = current_prec / 3600.0
        snow_infiltration_rate = snow_rate_flux / snow_penetration_depth
        
        for i in range(min(n_surface_nodes, n_nodes)):
            depth = i * dx
            decay_factor = nm.exp(-depth / (snow_penetration_depth / 3.0))
            snow_source[i] += snow_infiltration_rate * decay_factor * 0.1
    
    # Snow melting
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
    
    # Sublimation
    for i in range(min(n_surface_nodes, n_nodes)):
        if T_vals[i] < 0.0:
            current_rh_local = rh[hour_idx] if hour_idx < len(rh) else 70.0
            if current_rh_local < 90.0:
                sublimation_rate = 1e-6 * (90.0 - current_rh_local) / 90.0
                snow_source[i] -= sublimation_rate
    
    val = snow_source.reshape((n_nodes, 1, 1))
    return {'val': val}

def get_L_v_eps_L(ts, coors, mode=None, equations=None, term=None,
                 problem=None, **kwargs):
    """Calculate L_v * ε_L = L_v * (Cl / ρ_water)"""
    if mode != 'qp' or coors is None:
        return {}
    
    variables = problem.get_variables()
    Cl_var = variables['Cl']
    Cl_vals = Cl_var.get_state_in_region(problem.domain.regions['Omega'])
    Cl_vals_flat = Cl_vals.flatten()
    
    # Temperature-dependent latent heat
    L = nm.where((variables['T'].get_state_in_region(problem.domain.regions['Omega']).flatten()) < 0.0,
                    L_s, L_v)
    eps_L = Cl_vals_flat / 1000.0
    val = (L * eps_L).reshape((-1,1,1))
    return {'val': val}

def mesh_hook(mesh, mode):
    """Generate the 1D mesh"""
    if mode == 'read':
        coors = nm.linspace(0.0, d_ins, nodes).reshape((nodes, 1))
        conn = nm.arange(nodes, dtype=nm.int32).repeat(2)[1:-1].reshape((-1, 2))
        mat_ids = nm.zeros(nodes - 1, dtype=nm.int32)
        descs = ['1_2']
        mesh = Mesh.from_data('heat_liquid_1d', coors, None, [conn], [mat_ids], descs)
        return mesh
    elif mode == 'write':
        pass

###############################################################################
#                                MAIN SECTION                                 #
###############################################################################

# Read data and calculate boundary conditions
airTemp, airVel, prec, gloSolIr, rh = read_input_data()

# Calculation of equivalent sol-air temperature (in °C)
h = 22.7; alpha = 0.8; T_cor_fact = 4.0
t_o = [alpha * glob_solir / h + air_temp - T_cor_fact for 
       glob_solir, air_temp in zip(gloSolIr, airTemp)]

# Calculation of convective heat transfer coefficient 
h_o = [6.0 + 4.0*vel if vel <= 5.0 else 7.41*(vel**0.78) for vel in airVel]

# Boundary conditions
t_i = 0.0; h_i = 99.75
cl_i = 0.01; cl_o = 0.02; h_l_i = 1e-6; h_l_o = 2e-6
cv_i = 0.005; cv_o = 0.008; h_v_i = 1e-5
h_v_o = calculate_vapor_bc_coefficient(airTemp, rh)

# Mesh generation
filename_mesh = UserMeshIO(mesh_hook)

# Enhanced materials with pressure field components
materials = {
    'mat': ({'lam': lam_i, 'rho_cp': rho_wet * c_wet},),
    'liquid_mat': ({'D_l': D_l},),
    'vapor_mat': ({'D_v': D_v},),
    'pressure_mat': ({'D_p_eff': D_p_eff},),  # Effective pressure diffusion
    'thermal_mat': ({'D_t': D_t},),
    'coupling_mat': ({'k_c': k_c},),
    'pressure_heat_mat': ({'k_p': k_p},),  # Pressure-heat coupling
    
    # Pressure coupling coefficients
    'pressure_thermal': ({'alpha_T': alpha_T},),
    'pressure_liquid': ({'alpha_L': alpha_L},),
    'pressure_vapor': ({'alpha_V': alpha_V},),
    
    # Enhanced transport coefficients
    'thermal_vapor_mat': ({'D_t_v': D_t_v},),
    'pressure_liquid_mat': ({'D_p': D_p},),
    'pressure_vapor_mat': ({'D_p_v': D_p_v},),
    
    # Dynamic boundary conditions
    'h_out_dyn': 'get_h_o',
    'h_v_out_dyn': 'get_h_v_o',
    'T_out_dyn': 'get_t_o', 
    'cl_out_dyn': 'get_cl_o',  
    'cv_out_dyn': 'get_cv_o',  
    'cv_in_dyn': 'get_cv_i',
    'p_atm_dyn': 'get_p_atm',
    'h_p_out_dyn': 'get_h_p_o',
    'prec_cl_ref': 'get_prec_enhanced_cl_ref',
    'prec_h_l': 'get_prec_enhanced_h_l',
    
    # Fixed boundary conditions
    'in_fixed': ({'h_in': h_i, 'T_in': t_i, 'h_l_in': h_l_i, 
                  'Cl_in': cl_i, 'h_v_in': h_v_i, 'Cv_in': cv_i, 
                  'P_in': P_base, 'h_p_in': h_p_i},),
    'out_fixed': ({'h_l_out': h_l_o},),
    
    # Enhanced transport terms
    'evaporation_rate': 'get_evaporation_rate',
    'latent_heat_coeff': 'get_L_v_eps_L',
    'snow_accumulation': 'get_accumulated_snow_source',
    'heat_moisture_coupling': 'get_heat_moist_coup',
}

functions = {
    'get_h_o': (get_h_o,),
    'get_h_v_o': (get_h_v_o,),
    'get_t_o': (get_t_o,),
    'get_cl_o': (get_cl_o,),
    'get_cv_o': (get_cv_o,),
    'get_cv_i': (get_cv_i,),
    'get_p_atm': (get_p_atm,),
    'get_h_p_o': (get_h_p_o,),
    'get_evaporation_rate': (get_evaporation_rate,),
    'get_heat_moist_coup': (get_heat_moist_coup,),
    'get_L_v_eps_L': (get_L_v_eps_L,),
    'get_prec_enhanced_cl_ref': (get_prec_enhanced_cl_ref,),
    'get_prec_enhanced_h_l': (get_prec_enhanced_h_l,),
    'get_accumulated_snow_source': (get_accumulated_snow_source,),
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
    'pressure': ('real', 1, 'Omega', 1),
}

variables = {
    'T': ('unknown field', 'temperature', 0, 1),
    's': ('test field', 'temperature', 'T'),
    'Cl': ('unknown field', 'liquid', 1, 1),
    'r': ('test field', 'liquid', 'Cl'),
    'Cv': ('unknown field', 'vapor', 2, 1),
    'w': ('test field', 'vapor', 'Cv'),
    'P': ('unknown field', 'pressure', 3, 1),
    'q': ('test field', 'pressure', 'P'),
}

integrals = {'i': 1,}

# ENHANCED EQUATIONS with realistic pressure field
equations = {
    'Heat': """
            dw_dot.i.Omega(mat.rho_cp, s, dT/dt)
            + dw_laplace.i.Omega(mat.lam, s, T)
            + dw_laplace.i.Omega(coupling_mat.k_c, s, Cl)
            + dw_laplace.i.Omega(pressure_heat_mat.k_p, s, P)
            + dw_dot.i.Omega(latent_heat_coeff.val, s, dCl/dt)
            = - dw_bc_newton.i.Gamma_Left(h_out_dyn.val, T_out_dyn.val, s, T)
              - dw_bc_newton.i.Gamma_Right(in_fixed.h_in, in_fixed.T_in, s, T)
            """,
             
    'Liquid': """
            dw_dot.i.Omega(r, dCl/dt)
            + dw_laplace.i.Omega(liquid_mat.D_l, r, Cl)
            + dw_laplace.i.Omega(thermal_mat.D_t, r, T)
            + dw_laplace.i.Omega(pressure_liquid_mat.D_p, r, P)
            = - dw_volume_lvf.i.Omega(evaporation_rate.val, r)
              - dw_bc_newton.i.Gamma_Left(prec_h_l.val, prec_cl_ref.val, r, Cl)
              - dw_bc_newton.i.Gamma_Right(in_fixed.h_l_in, in_fixed.Cl_in, r, Cl)
              + dw_volume_lvf.i.Omega(snow_accumulation.val, r)
            """,
               
    'Vapor': """
            dw_dot.i.Omega(w, dCv/dt)
            + dw_laplace.i.Omega(vapor_mat.D_v, w, Cv)
            + dw_laplace.i.Omega(thermal_vapor_mat.D_t_v, r, T)
            + dw_laplace.i.Omega(pressure_vapor_mat.D_p_v, r, P)
            = + dw_volume_lvf.i.Omega(evaporation_rate.val, w)
              - dw_bc_newton.i.Gamma_Left(h_v_out_dyn.val, cv_out_dyn.val, w, Cv)
              - dw_bc_newton.i.Gamma_Right(in_fixed.h_v_in, in_fixed.Cv_in, w, Cv)
            """,
    
    'Pressure': """
            dw_dot.i.Omega(q, dP/dt)
            + dw_laplace.i.Omega(pressure_mat.D_p_eff, q, P)
            + dw_laplace.i.Omega(pressure_thermal.alpha_T, q, T)
            + dw_laplace.i.Omega(pressure_liquid.alpha_L, q, Cl)
            + dw_laplace.i.Omega(pressure_vapor.alpha_V, q, Cv)
            = - dw_bc_newton.i.Gamma_Left(h_p_out_dyn.val, p_atm_dyn.val, q, P)
              - dw_bc_newton.i.Gamma_Right(in_fixed.h_p_in, in_fixed.P_in, q, P)
            """
}

# Boundary conditions (empty but required)
ebcs = {}

# Initial conditions
ics = {
    'ic_T': ('Omega', {'T.0': 0.0}),
    'ic_Cl': ('Omega', {'Cl.0': 0.015}),
    'ic_Cv': ('Omega', {'Cv.0': 0.006}),
    'ic_P': ('Omega', {'P.0': P_base}),  # Initialize pressure to base value
}

# Time and solver parameters
nr_of_hours = 2
stop = nr_of_hours * 3600.0
dt = 10.0
nr_of_steps = int(stop/dt)

solvers = {
    'ls': ('ls.scipy_direct', {}),
    'newton': ('nls.newton', {
        'i_max': 30,      # Increased for enhanced nonlinear system
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
    """Enhanced results saving with pressure field"""
    
    filename = os.path.join(problem.conf.options['output_dir'], "enhanced_moisture_pressure_results_v8.csv")
    header = [
        "Time (s)", "Node Index", "Position (m)", 
        "Temperature (°C)", "Liquid Content (kg/m³)", "Vapor Content (kg/m³)",
        "Pressure (Pa)", "Pressure_rel (Pa)"
    ]
    
    # Get state
    state_vec = state.get_state()
    coors = problem.fields['temperature'].get_coor()
    n_dof = problem.fields['temperature'].n_nod
    
    T_vals = state_vec[:n_dof]
    Cl_vals = state_vec[n_dof:2*n_dof]
    Cv_vals = state_vec[2*n_dof:3*n_dof]
    P_vals = state_vec[3*n_dof:4*n_dof]
    P_rel_vals = P_vals - p_atm  # Relative pressure (gauge pressure)
    
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
                    problem.ts.time, i, coord[0], T_vals[i], Cl_vals[i], 
                    Cv_vals[i], P_vals[i], P_rel_vals[i]
                ])
        
        print(f"Hour {int(problem.ts.time/3600)}: Enhanced physics with pressure simulation")
        print(f"  Pressure range: {P_vals.min():.1f} - {P_vals.max():.1f} Pa")
        print(f"  Relative pressure range: {P_rel_vals.min():.1f} - {P_rel_vals.max():.1f} Pa")
    
    return out

options = {
    'nls': 'newton',
    'ls': 'ls', 
    'ts': 'ts',
    'save_times': [3600*i for i in range(1, nr_of_hours + 1)],
    'post_process_hook': save_enhanced_results,
    'output_dir': './output_snow_pressure_v8',
}