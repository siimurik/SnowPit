r"""
==================================================================
Step 4: Heat-Liquid-Vapor Coupled Transport System
------------------------------------------------------------------
Three-field coupled system:
    Heat:   ∂T/∂t = α ∂²T/∂x² + L_v * m_evap/(ρ*cp) + q_dot/(ρ*cp)
    Liquid: ∂Cl/∂t = D_l * ∂²Cl/∂x² - m_evap
    Vapor:  ∂Cv/∂t = D_v * ∂²Cv/∂x² + m_evap

where m_evap = f(T, Cl, Cv) couples all three equations
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

# Material properties
lam_i = 0.32      # Thermal conductivity [W/(mK)]
c_wet = 2.59E03   # J/(kg*K)
rho_dry = 100.0   # kg/m^3
moist_cont = 50.0 # %
rho_wet = rho_dry + moist_cont*10.0 # kg/m^3
D = lam_i/(c_wet * rho_wet) # m^2/s

# Enhanced parameters for liquid and vapor transport
D_l = 1e-8        # Liquid diffusion coefficient [m^2/s]
D_v = 2e-5        # Vapor diffusion coefficient [m^2/s] (much faster than liquid)
L_v = 2.45e6      # Latent heat of vaporization [J/kg]
q_dot = 0.0       # Heat source [W/m^3]

# Enhanced evaporation model parameters
k_evap = 1e-7     # Increased evaporation rate constant [1/s]
T_ref = 0.0       # Reference temperature [°C]
Cl_min = 0.001    # Minimum liquid content for evaporation [kg/m³]
T_evap_max = 10.0 # Maximum temperature for linear evaporation model [°C]

def Psat_WV(T_K):
    """
    Water vapour saturation pressure
    W. Wagner and A. Pruß:" The IAPWS Formulation 1995 for the
    Thermodynamic Properties of Ordinary Water Substance for General and Scientific Use ",
    Journal of Physical and Chemical Reference Data, June 2002 ,Volume 31, Issue 2, pp.
    387535)
    Returns Saturation vapor pressure (hPa)
    """
    Tc = 647.096  # Critical temperature, K
    Pc = 220640   # Critical pressure, hPa
    C1 = -7.85951783
    C2 = 1.84408259
    C3 = -11.7866497
    C4 = 22.6807411
    C5 = -15.9618719
    C6 = 1.80122502
    teta = 1 - T_K / Tc
    x = Tc / T_K * (C1 * teta + C2 * teta ** 1.5 + C3 * teta ** 3 + C4 * teta ** 3.5 + C5 * teta ** 4 + C6 * teta ** 7.5)
    x = math.exp(x) * Pc
    return x

def read_temp_and_hcoeff_from_csv(filename="t_rh_to_ho.csv"):
    """Read temperature and heat transfer coefficient data from CSV file."""
    time, rh, t_o, h_o = [], [], [], []
    with open(filename, 'r') as csvfile:
        csvreader = csv.reader(csvfile)
        for row in csvreader:
            time.append(float(row[0].strip()))  # time in seconds
            rh.append(float(row[1].strip()))    # relative humidity %
            t_o.append(float(row[2].strip()))   # outer layer temperature
            h_o.append(float(row[3].strip()))   # outer layer convection 

    return time, rh, t_o, h_o

# Boundary conditions
t_i = 0.0         # Inner temperature [°C]
h_i = 99.75       # Inner heat transfer coefficient [W/m2K]

# Liquid boundary conditions (Robin type)
cl_i = 0.01       # Inner reference liquid content [kg/m^3]
cl_o = 0.02       # Outer reference liquid content [kg/m^3] 
h_l_i = 1e-6      # Inner liquid mass transfer coefficient [m/s]
h_l_o = 2e-6      # Outer liquid mass transfer coefficient [m/s]

# Vapor boundary conditions (Robin type)
cv_i = 0.005      # Inner reference vapor content [kg/m^3]
cv_o = 0.008      # Outer reference vapor content [kg/m^3]
h_v_i = 1e-5      # Inner vapor mass transfer coefficient [m/s]  
h_v_o = 2e-5      # Outer vapor mass transfer coefficient [m/s]

time, rh, t_o, h_o = read_temp_and_hcoeff_from_csv()

def mesh_hook(mesh, mode):
    """Generate the 1D mesh."""
    if mode == 'read':
        coors = nm.linspace(0.0, d_ins, nodes).reshape((nodes, 1))
        conn = nm.arange(nodes, dtype=nm.int32).repeat(2)[1:-1].reshape((-1, 2))
        mat_ids = nm.zeros(nodes - 1, dtype=nm.int32)
        descs = ['1_2']

        mesh = Mesh.from_data('heat_liquid_1d', coors, None,
                             [conn], [mat_ids], descs)
        return mesh
    elif mode == 'write':
        pass

def get_h_o(ts, coors, mode=None, **kwargs):
    """Time-dependent heat transfer coefficient."""
    if mode != 'qp' or coors is None:
        return {}

    hour_idx = min(int(ts.time / 3600), len(h_o) - 1)
    val = nm.full((coors.shape[0], 1, 1), h_o[hour_idx], dtype=nm.float64)

    return {'val': val}

def get_t_o(ts, coors, mode=None, **kwargs):
    """Time-dependent ambient temperature."""
    if mode != 'qp' or coors is None:
        return {}

    hour_idx = min(int(ts.time / 3600), len(t_o) - 1)
    val = nm.full((coors.shape[0], 1, 1), t_o[hour_idx], dtype=nm.float64)

    return {'val': val}

def get_rh(ts, coors, mode=None, **kwargs):
    """Time-dependent relative humidity."""
    if mode != 'qp' or coors is None:
        return {}

    hour_idx = min(int(ts.time / 3600), len(rh) - 1)
    val = nm.full((coors.shape[0], 1, 1), rh[hour_idx], dtype=nm.float64)

    return {'val': val}

def get_cv_o(ts, coors, mode=None, **kwargs):
    """
    Time-dependent outer vapor reference content.
    Based on relative humidity and saturation vapor pressure.
    """
    if mode != 'qp' or coors is None:
        return {}

    # Get current environmental conditions
    hour_idx = min(int(ts.time / 3600), len(rh) - 1)
    current_rh = rh[hour_idx] if hour_idx < len(rh) else 70.0
    current_t = t_o[hour_idx] if hour_idx < len(t_o) else -5.0
    
    # Calculate saturation vapor pressure and actual vapor content
    T_K = current_t + 273.15  # Convert to Kelvin
    try:
        Psat = Psat_WV(T_K)  # Saturation pressure in hPa
        # Convert to vapor density (approximate ideal gas law)
        # ρ_v = P * M / (R * T) where M ≈ 0.018 kg/mol, R = 8.314 J/(mol*K)
        rho_v_sat = Psat * 100 * 0.018 / (8.314 * T_K)  # kg/m³
        rho_v_actual = rho_v_sat * current_rh / 100.0
        cv_adjusted = max(cv_o * 0.1, min(cv_o * 2.0, rho_v_actual))
    except:
        # Fallback if saturation pressure calculation fails
        cv_adjusted = cv_o * (current_rh / 100.0)
    
    val = nm.full((coors.shape[0], 1, 1), cv_adjusted, dtype=nm.float64)

    return {'val': val}
def get_cl_o(ts, coors, mode=None, **kwargs):
    """
    Time-dependent outer liquid reference content.
    Now includes humidity effects for better physical realism.
    """
    if mode != 'qp' or coors is None:
        return {}

    # Get current relative humidity and temperature
    hour_idx = min(int(ts.time / 3600), len(rh) - 1)
    current_rh = rh[hour_idx]
    current_t = t_o[hour_idx]
    
    # Adjust reference liquid content based on humidity
    # Higher humidity -> higher reference liquid content
    cl_base = cl_o
    humidity_factor = current_rh / 100.0  # Convert % to fraction
    cl_adjusted = cl_base * (0.5 + 0.5 * humidity_factor)  # Range: 0.5*cl_o to cl_o
    
    val = nm.full((coors.shape[0], 1, 1), cl_adjusted, dtype=nm.float64)

    return {'val': val}

def get_cv_i(ts, coors, mode=None, **kwargs):
    """Time-dependent inner vapor reference content."""
    if mode != 'qp' or coors is None:
        return {}

    val = nm.full((coors.shape[0], 1, 1), cv_i, dtype=nm.float64)
    return {'val': val}

def get_evaporation_rate(ts, coors, mode=None, equations=None, term=None, 
                        problem=None, **kwargs):
    """
    Enhanced temperature and liquid-content dependent evaporation rate.
    Now includes humidity effects and more robust state access using SfePy API.
    
    Evaporation model:
    m_evap = k_evap * f_humidity * f_temperature * f_liquid
    
    where:
    - f_humidity = (1 - RH/100) accounts for ambient humidity
    - f_temperature = max(0, (T - T_ref)/(T_evap_max - T_ref)) 
    - f_liquid = max(0, (Cl - Cl_min)/Cl_min) accounts for available liquid
    """
    if mode != 'qp' or coors is None:
        return {}
    
    try:
        # Get current relative humidity
        hour_idx = min(int(ts.time / 3600), len(rh) - 1)
        current_rh = rh[hour_idx]
        humidity_factor = max(0, 1.0 - current_rh/100.0)  # 0 at 100% RH, 1 at 0% RH
        
        # Get current state from problem using correct SfePy API
        if problem is not None and hasattr(problem, 'get_variables'):
            variables = problem.get_variables()
            
            # Get current variable values
            if 'T' in variables and 'Cl' in variables:
                T_var = variables['T']
                Cl_var = variables['Cl']
                
                # Get the state vectors
                T_vals = T_var.get_state_in_region(problem.domain.regions['Omega'])
                Cl_vals = Cl_var.get_state_in_region(problem.domain.regions['Omega'])
                
                # Also get vapor content if available
                Cv_vals = None
                if 'Cv' in variables:
                    Cv_var = variables['Cv']
                    Cv_vals = Cv_var.get_state_in_region(problem.domain.regions['Omega'])
                
                # Temperature factor: linear increase from T_ref to T_evap_max
                T_factor = nm.maximum(0, nm.minimum(1, (T_vals - T_ref)/(T_evap_max - T_ref)))
                
                # Liquid availability factor
                Cl_factor = nm.maximum(0, (Cl_vals - Cl_min)/Cl_min)
                
                # Vapor saturation factor (reduces evaporation when vapor content is high)
                if Cv_vals is not None:
                    # Improved: Use local saturation vapor density via Psat_WV and ideal gas law
                    T_K = T_vals + 273.15
                    try:
                        Psat = Psat_WV(T_K) * 100.0  # hPa to Pa
                        # Ideal gas law: rho_v_sat = Psat * M / (R * T)
                        M_wv = 0.018     # kg/mol, molar mass of water vapor
                        R_gas = 8.314    # J/(mol*K)
                        rho_v_sat = Psat * M_wv / (R_gas * T_K)  # kg/m³
                        # Vapor factor: how close is Cv to saturation (limits evaporation)
                        vapor_factor = nm.maximum(0, 1.0 - Cv_vals / rho_v_sat)
                    except Exception as e:
                        vapor_factor = nm.ones_like(T_vals)
                else:
                    vapor_factor = nm.ones_like(T_vals)
                
                # Combined evaporation rate
                evap_rate = k_evap * humidity_factor * T_factor * Cl_factor * vapor_factor
                
                # Ensure we have the right shape for quadrature points
                if len(evap_rate) != coors.shape[0]:
                    # Simple interpolation/averaging if dimensions don't match
                    evap_rate = nm.full(coors.shape[0], evap_rate.mean())
            else:
                # Variables not available, use fallback
                evap_rate = nm.full(coors.shape[0], k_evap * humidity_factor * 0.5 * 1.0)
            
        else:
            # Problem not available or no variables method, use fallback
            evap_rate = nm.full(coors.shape[0], k_evap * humidity_factor * 0.5 * 1.0)
            
    except Exception as e:
        # Robust fallback - don't print warnings every iteration to avoid spam
        if not hasattr(get_evaporation_rate, 'warning_printed'):
            print(f"Note: Using simplified evaporation model - {e}")
            get_evaporation_rate.warning_printed = True
            
        hour_idx = min(int(ts.time / 3600), len(rh) - 1)
        current_rh = rh[hour_idx] if hour_idx < len(rh) else 70.0
        humidity_factor = max(0, 1.0 - current_rh/100.0)
        evap_rate = nm.full(coors.shape[0], k_evap * humidity_factor * 0.5)
    
    val = evap_rate.reshape((coors.shape[0], 1, 1))
    return {'val': val}

def get_latent_heat_source(ts, coors, mode=None, equations=None, term=None,
                          problem=None, **kwargs):
    """
    Latent heat source: L_v * ∂Cl/∂t ≈ -L_v * m_evap
    This represents the cooling effect due to evaporation.
    """
    if mode != 'qp' or coors is None:
        return {}
    
    # Get evaporation rate using the same robust approach
    evap_data = get_evaporation_rate(ts, coors, mode, equations, term, problem, **kwargs)
    evap_rate = evap_data['val'].flatten()
    
    # Latent heat source = -L_v * evaporation_rate
    # Negative because evaporation removes liquid (cooling effect)
    latent_source = -L_v * evap_rate / (rho_wet * c_wet)  # Normalized by heat capacity
    
    val = latent_source.reshape((coors.shape[0], 1, 1))
    return {'val': val}

def get_heat_source(ts, coors, mode=None, **kwargs):
    """Heat source term."""
    if mode != 'qp' or coors is None:
        return {}
    
    val = nm.full((coors.shape[0], 1, 1), q_dot, dtype=nm.float64)
    return {'val': val}

filename_mesh = UserMeshIO(mesh_hook)

materials = {
    'mat': ({'lam': lam_i, 'rho_cp': rho_wet * c_wet},),
    'liquid_mat': ({'D_l': D_l},),
    'vapor_mat': ({'D_v': D_v},),  # Added vapor diffusion material
    'h_out_dyn': 'get_h_o',
    'T_out_dyn': 'get_t_o', 
    'Cl_out_dyn': 'get_cl_o',
    'Cv_out_dyn': 'get_cv_o',  # Added vapor boundary condition
    'Cv_in_dyn': 'get_cv_i',   # Added inner vapor boundary condition
    'rh_dyn': 'get_rh',
    'in_fixed': ({'h_in': h_i, 'T_in': t_i, 'h_l_in': h_l_i, 'Cl_in': cl_i, 
                  'h_v_in': h_v_i},),  # Added vapor transfer coeff
    'out_fixed': ({'h_l_out': h_l_o, 'h_v_out': h_v_o},),  # Added vapor transfer coeff
    'evaporation_rate': 'get_evaporation_rate',
    'latent_heat_source': 'get_latent_heat_source',
    'heat_source': 'get_heat_source',
}

functions = {
    'get_h_o': (get_h_o,),
    'get_t_o': (get_t_o,),
    'get_cl_o': (get_cl_o,),
    'get_cv_o': (get_cv_o,),  # Added vapor boundary function
    'get_cv_i': (get_cv_i,),  # Added inner vapor boundary function
    'get_rh': (get_rh,),
    'get_evaporation_rate': (get_evaporation_rate,),
    'get_latent_heat_source': (get_latent_heat_source,),
    'get_heat_source': (get_heat_source,),
}

regions = {
    'Omega': 'all',
    'Gamma_Left': ('vertices in (x < 0.00001)', 'facet'),
    'Gamma_Right': ('vertices in (x > %f)' % (d_ins - 0.00001), 'facet'), 
}

# Three fields: temperature, liquid content, and vapor content
fields = {
    'temperature': ('real', 1, 'Omega', 1),
    'liquid': ('real', 1, 'Omega', 1),
    'vapor': ('real', 1, 'Omega', 1), # Added vapor field
}
# Four variables with history for time derivatives
variables = {
    'T': ('unknown field', 'temperature', 0, 1),
    's': ('test field', 'temperature', 'T'),
    'Cl': ('unknown field', 'liquid', 1, 1),
    'r': ('test field', 'liquid', 'Cl'),
    'Cv': ('unknown field', 'vapor', 2, 1),  # Fixed: history size = 2, order = 1
    'w': ('test field', 'vapor', 'Cv') # Added vapor test field
}

# Boundary conditions
ebcs = {
    # Robin BCs are handled in the equations
}

integrals = {
    'i': 1,
}

# Three-field coupled equations: Heat-Liquid-Vapor
equations = {
    'Heat': """dw_dot.i.Omega(mat.rho_cp, s, dT/dt)
             + dw_laplace.i.Omega(mat.lam, s, T)
             = - dw_bc_newton.i.Gamma_Left(h_out_dyn.val, T_out_dyn.val, s, T)
               - dw_bc_newton.i.Gamma_Right(in_fixed.h_in, in_fixed.T_in, s, T)
               + dw_volume_lvf.i.Omega(latent_heat_source.val, s)
               + dw_volume_lvf.i.Omega(heat_source.val, s)
             """,
             
    'Liquid': """dw_dot.i.Omega(r, dCl/dt)
               + dw_laplace.i.Omega(liquid_mat.D_l, r, Cl)
               = - dw_volume_lvf.i.Omega(evaporation_rate.val, r)
                 - dw_bc_newton.i.Gamma_Left(out_fixed.h_l_out, Cl_out_dyn.val, r, Cl)
                 - dw_bc_newton.i.Gamma_Right(in_fixed.h_l_in, in_fixed.Cl_in, r, Cl)
               """,
               
    'Vapor': """dw_dot.i.Omega(w, dCv/dt)
              + dw_laplace.i.Omega(vapor_mat.D_v, w, Cv)
              = + dw_volume_lvf.i.Omega(evaporation_rate.val, w)
                - dw_bc_newton.i.Gamma_Left(out_fixed.h_v_out, Cv_out_dyn.val, w, Cv)
                - dw_bc_newton.i.Gamma_Right(in_fixed.h_v_in, Cv_in_dyn.val, w, Cv)
              """
}

# Initial conditions for all three fields
ics = {
    'ic_T': ('Omega', {'T.0': 0.0}),
    'ic_Cl': ('Omega', {'Cl.0': 0.015}),   # Initial liquid content
    'ic_Cv': ('Omega', {'Cv.0': 0.006}),   # Initial vapor content
}


# Time parameters
nr_hour = len(t_o)
start = 0.0
nr_of_hours = 2
stop = nr_of_hours * 3600.0
dt = 10.0
nr_of_steps = int(stop/dt)

solvers = {
    'ls': ('ls.scipy_direct', {}),
    'newton': ('nls.newton', {
        'i_max': 20,      # Increased iterations for nonlinear system
        'eps_a': 1e-8,
        'eps_r': 1e-6,    # Added relative tolerance
        'is_linear': False,
    }),
    'ts': ('ts.simple', {
        't0': start,
        't1': stop,
        'dt': dt,
        'n_step': nr_of_steps,
        'verbose': 1,
    }),
}

    
def save_coupled_results(out, problem, state, extend=False):
    """Save temperature, liquid content, and vapor content with diagnostics."""
    import os

    filename = os.path.join(problem.conf.options['output_dir'], "coupled_results_step4.csv")
    header = [
        "Time (s)", "Node Index", "Position (m)", 
        "Temperature (°C)", "Liquid Content (kg/m³)", "Vapor Content (kg/m³)",
        "Evaporation Rate (kg/m³s)", "Latent Heat Effect (W/m³)",
        "RH (%)", "T_outdoor (°C)"
    ]

    # Get the full state vector
    state_vec = state.get_state()
    coors = problem.fields['temperature'].get_coor()
    n_dof = problem.fields['temperature'].n_nod

    # Handle different state vector configurations
    if len(state_vec) == 3 * n_dof:
        T_vals = state_vec[:n_dof]
        Cl_vals = state_vec[n_dof:2*n_dof]
        Cv_vals = state_vec[2*n_dof:3*n_dof]
    elif len(state_vec) == 2 * n_dof:
        T_vals = state_vec[:n_dof]
        Cl_vals = state_vec[n_dof:2*n_dof]
        Cv_vals = nm.zeros(n_dof)
    else:
        T_vals = state_vec[:n_dof]
        Cl_vals = nm.zeros(n_dof)
        Cv_vals = nm.zeros(n_dof)

    # Calculate evaporation rate for diagnostics
    try:
        evap_data = get_evaporation_rate(problem.ts, coors[:, 0], 'qp', problem=problem)
        evap_rates = evap_data['val'].flatten()
    except:
        evap_rates = nm.zeros(len(T_vals))

    # Calculate latent heat effect
    try:
        latent_data = get_latent_heat_source(problem.ts, coors[:, 0], 'qp', problem=problem)
        latent_effects = latent_data['val'].flatten()
    except:
        latent_effects = nm.zeros(len(T_vals))

    # Get current environmental conditions
    hour_idx = min(int(problem.ts.time / 3600), len(rh) - 1)
    current_rh = rh[hour_idx] if hour_idx < len(rh) else 70.0
    current_t_o = t_o[hour_idx] if hour_idx < len(t_o) else -5.0

    # Write header if file does not exist or is empty
    file_needs_header = not os.path.exists(filename) or os.path.getsize(filename) == 0
    if file_needs_header:
        with open(filename, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(header)

    # Append data for each node at full-hour marks
    if problem.ts.time % 3600 < problem.ts.dt:  
        with open(filename, 'a', newline='') as csvfile:
            writer = csv.writer(csvfile)
            for i, coord in enumerate(coors):
                writer.writerow([
                    problem.ts.time, i, coord[0], T_vals[i], Cl_vals[i], Cv_vals[i],
                    evap_rates[i], latent_effects[i], current_rh, current_t_o
                ])

    # Print summary every hour for monitoring
    if problem.ts.time % 3600 < problem.ts.dt:
        print(f"\nHour {int(problem.ts.time/3600)}: T_avg={T_vals.mean():.2f}°C, "
              f"Cl_avg={Cl_vals.mean():.4f} kg/m³, "
              f"Cv_avg={Cv_vals.mean():.4f} kg/m³, "
              f"Total_evap={evap_rates.sum():.2e} kg/m³s, "
              f"RH={current_rh:.1f}%")

    return out

options = {
    'nls': 'newton',
    'ls': 'ls',
    'ts': 'ts',
    'save_times': [3600*i for i in range(1, nr_of_hours + 1)],
    'post_process_hook': save_coupled_results,
    'output_dir': './output_snow_step4',
}