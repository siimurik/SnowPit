r"""
==================================================================
Step 2: Coupled Heat-Liquid Transport System
------------------------------------------------------------------
Solves the coupled system:
    Heat:   ∂T/∂t = α ∂²T/∂x² + L_v * ∂Cl/∂t / (ρ*cp) + q_dot / (ρ*cp)
    Liquid: ∂Cl/∂t = D_l ∂²Cl/∂x² - m_evap

where m_evap = f(T, Cl) couples the equations
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

# Material properties
lam_i = 0.32      # Thermal conductivity [W/(mK)]
c_wet = 2.59E03   # J/(kg*K)
rho_dry = 100.0   # kg/m^3
moist_cont = 50.0 # %
rho_wet = rho_dry + moist_cont*10.0 # kg/m^3
D = lam_i/(c_wet * rho_wet) # m^2/s

# New parameters for liquid transport
D_l = 1e-8        # Liquid diffusion coefficient [m^2/s]
L_v = 2.45e6      # Latent heat of vaporization [J/kg]
q_dot = 0.0       # Heat source [W/m^3]

# Evaporation model parameters
k_evap = 1e-8     # Evaporation rate constant [1/s]
T_ref = 0.0       # Reference temperature [°C]

def read_temp_and_hcoeff_from_csv(filename="t_o_and_h_o.csv"):
    """Read temperature and heat transfer coefficient data from CSV file."""
    t_o, h_o = [], []
    with open(filename, 'r') as csvfile:
        csvreader = csv.reader(csvfile)
        for row in csvreader:
            if len(row) >= 2:
                t_o.append(float(row[0].strip()))
                h_o.append(float(row[1].strip()))
    return t_o, h_o

# Boundary conditions
t_i = 0.0         # Inner temperature [°C]
h_i = 99.75       # Inner heat transfer coefficient [W/m2K]

# Liquid boundary conditions (Robin type)
cl_i = 0.01       # Inner reference liquid content [kg/m^3]
cl_o = 0.02       # Outer reference liquid content [kg/m^3] # TODO! After one hourly data use the moisture from the previous hour
h_l_i = 1e-6      # Inner liquid mass transfer coefficient [m/s]
h_l_o = 2e-6      # Outer liquid mass transfer coefficient [m/s]

t_o, h_o = read_temp_and_hcoeff_from_csv()

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

def get_cl_o(ts, coors, mode=None, **kwargs):
    """
    Time-dependent outer liquid reference content.
    Could vary with humidity, weather conditions, etc.
    For now, keep constant but structure allows for variation.
    """
    if mode != 'qp' or coors is None:
        return {}

    # Could make this time-dependent based on weather data
    # For now, use constant reference value
    val = nm.full((coors.shape[0], 1, 1), cl_o, dtype=nm.float64)

    return {'val': val}

def get_evaporation_rate(ts, coors, mode=None, equations=None, term=None, 
                        problem=None, **kwargs):
    """
    Temperature and liquid-content dependent evaporation rate.
    m_evap = k_evap * Cl * max(0, T - T_ref)
    
    This couples the heat and liquid equations.
    """
    if mode != 'qp' or coors is None:
        return {}
    
    # Get current state - we need both T and Cl fields
    if equations is not None and hasattr(equations, 'variables'):
        try:
            # Get current temperature and liquid content
            state = problem.get_current_state()
            
            # Extract temperature values (first field)
            n_dof_T = problem.fields['temperature'].n_nod
            T_vals = state[:n_dof_T]
            
            # Extract liquid content values (second field)  
            Cl_vals = state[n_dof_T:n_dof_T*2]
            
            # Interpolate to quadrature points if needed
            # For simplicity, assume we're at nodes
            if len(T_vals) == coors.shape[0]:
                T_qp = T_vals
                Cl_qp = Cl_vals
            else:
                # More complex interpolation would be needed for real quadrature points
                T_qp = nm.full(coors.shape[0], T_vals.mean())
                Cl_qp = nm.full(coors.shape[0], Cl_vals.mean())
            
            # Evaporation model: only evaporate if T > T_ref and Cl > 0
            evap_rate = k_evap * nm.maximum(0, Cl_qp) * nm.maximum(0, T_qp - T_ref)
            
        except:
            # Fallback: use simple constant rate
            evap_rate = nm.full(coors.shape[0], k_evap * 0.01 * 1.0)
    else:
        # Fallback: use simple constant rate  
        evap_rate = nm.full(coors.shape[0], k_evap * 0.01 * 1.0)
    
    val = evap_rate.reshape((coors.shape[0], 1, 1))
    return {'val': val}

def get_latent_heat_source(ts, coors, mode=None, equations=None, term=None,
                          problem=None, **kwargs):
    """
    Latent heat source: L_v * ∂Cl/∂t
    This will be computed from the evaporation rate.
    """
    if mode != 'qp' or coors is None:
        return {}
    
    # Get evaporation rate
    evap_data = get_evaporation_rate(ts, coors, mode, equations, term, problem, **kwargs)
    evap_rate = evap_data['val'].flatten()
    
    # Latent heat source = -L_v * evaporation_rate
    # Negative because evaporation removes liquid (cooling effect)
    latent_source = -L_v * evap_rate
    
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
    'h_out_dyn': 'get_h_o',
    'T_out_dyn': 'get_t_o', 
    'Cl_out_dyn': 'get_cl_o',
    'in_fixed': ({'h_in': h_i, 'T_in': t_i, 'h_l_in': h_l_i, 'Cl_in': cl_i},),
    'out_fixed': ({'h_l_out': h_l_o},),
    'evaporation_rate': 'get_evaporation_rate',
    'latent_heat_source': 'get_latent_heat_source',
    'heat_source': 'get_heat_source',
}

functions = {
    'get_h_o': (get_h_o,),
    'get_t_o': (get_t_o,),
    'get_cl_o': (get_cl_o,),
    'get_evaporation_rate': (get_evaporation_rate,),
    'get_latent_heat_source': (get_latent_heat_source,),
    'get_heat_source': (get_heat_source,),
}

regions = {
    'Omega': 'all',
    'Gamma_Left': ('vertices in (x < 0.00001)', 'facet'),
    'Gamma_Right': ('vertices in (x > %f)' % (d_ins - 0.00001), 'facet'), 
}

# Two fields: temperature and liquid content
fields = {
    'temperature': ('real', 1, 'Omega', 1),
    'liquid': ('real', 1, 'Omega', 1),
}

# Four variables: T, test_T, Cl, test_Cl
# Note: Added history=1 to enable time derivatives
variables = {
    'T': ('unknown field', 'temperature', 0, 1),
    's': ('test field', 'temperature', 'T'),
    'Cl': ('unknown field', 'liquid', 1, 1), 
    'r': ('test field', 'liquid', 'Cl'),
}

# Boundary conditions - Now using Robin BCs for liquid transport
ebcs = {
    # No more fixed liquid content - using Robin BCs instead
}

integrals = {
    'i': 1,
}

# Coupled equations with Robin BCs for both temperature and liquid
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
               """
}

# Initial conditions for both fields
ics = {
    'ic_T': ('Omega', {'T.0': 0.0}),
    'ic_Cl': ('Omega', {'Cl.0': 0.015}),  # Initial liquid content
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
        'i_max': 15,
        'eps_a': 1e-8,
        'is_linear': False,  # Now nonlinear due to coupling!
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
    """Save both temperature and liquid content at full-hour intervals."""
    filename = os.path.join(problem.conf.options['output_dir'], "coupled_results_step2.csv")
    
    # Get the full state vector
    state_vec = state.get_state()
    coors = problem.fields['temperature'].get_coor()
    
    # Split state vector into temperature and liquid parts
    n_dof = problem.fields['temperature'].n_nod
    T_vals = state_vec[:n_dof]
    Cl_vals = state_vec[n_dof:2*n_dof]
    
    # Write headers only for the first timestep
    if problem.ts.step == 0:
        with open(filename, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["Time (s)", "Node Index", "Position (m)", 
                           "Temperature (°C)", "Liquid Content (kg/m³)"])

    # Save data at full-hour marks
    if problem.ts.time % 3600 < problem.ts.dt:  
        with open(filename, 'a', newline='') as csvfile:
            writer = csv.writer(csvfile)
            for i, coord in enumerate(coors):
                writer.writerow([problem.ts.time, i, coord[0], T_vals[i], Cl_vals[i]])

    return out

options = {
    'nls': 'newton',
    'ls': 'ls',
    'ts': 'ts',
    'save_times': [3600*i for i in range(1, nr_of_hours + 1)],
    'post_process_hook': save_coupled_results,
    'output_dir': './output_snow_step2',
}