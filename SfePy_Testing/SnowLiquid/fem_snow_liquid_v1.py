r"""
==================================================================
Step 1: Heat equation with latent heat source term
------------------------------------------------------------------
Extends the basic heat equation to include:
    ∂T/∂t = α ∂²T/∂x² + L_v * ∂Cl/∂t / (ρ*cp) + q_dot / (ρ*cp)

This is the first step toward the full coupled system.
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

# New parameters for latent heat
L_v = 2.45e6      # Latent heat of vaporization [J/kg]
q_dot = 0.0       # Heat source [W/m^3] - start with zero

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
t_o, h_o = read_temp_and_hcoeff_from_csv()

def mesh_hook(mesh, mode):
    """Generate the 1D mesh."""
    if mode == 'read':
        coors = nm.linspace(0.0, d_ins, nodes).reshape((nodes, 1))
        conn = nm.arange(nodes, dtype=nm.int32).repeat(2)[1:-1].reshape((-1, 2))
        mat_ids = nm.zeros(nodes - 1, dtype=nm.int32)
        descs = ['1_2']

        mesh = Mesh.from_data('heat_1d', coors, None,
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

def get_latent_heat_source(ts, coors, mode=None, **kwargs):
    """
    Latent heat source term: L_v * dCl/dt
    Combines the latent heat and liquid change rate into a single source term.
    """
    if mode != 'qp' or coors is None:
        return {}
    
    # Simple sinusoidal variation for testing
    # Positive = condensation (releases heat)
    # Negative = evaporation (absorbs heat)
    dCl_dt = 1e-6 * nm.sin(2 * nm.pi * ts.time / 3600.0)  # kg/(m^3*s)
    latent_source = L_v * dCl_dt  # W/m^3
    val = nm.full((coors.shape[0], 1, 1), latent_source, dtype=nm.float64)
    
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
    'h_out_dyn': 'get_h_o',
    'T_out_dyn': 'get_t_o', 
    'in_fixed': ({'h_in': h_i, 'T_in': t_i},),
    'latent_heat_source': 'get_latent_heat_source',  # Combined L_v * dCl/dt
    'heat_source': 'get_heat_source',                # q_dot source
}

functions = {
    'get_h_o': (get_h_o,),
    'get_t_o': (get_t_o,),
    'get_latent_heat_source': (get_latent_heat_source,),
    'get_heat_source': (get_heat_source,),
}

regions = {
    'Omega': 'all',
    'Gamma_Left': ('vertices in (x < 0.00001)', 'facet'),
    'Gamma_Right': ('vertices in (x > %f)' % (d_ins - 0.00001), 'facet'), 
}

fields = {
    'temperature': ('real', 1, 'Omega', 1),
}

variables = {
    'T': ('unknown field', 'temperature', 0, 1),
    's': ('test field', 'temperature', 'T'),
}

ebcs = {}

integrals = {
    'i': 1,
}

# Enhanced heat equation with source terms
equations = {
    'Heat': """dw_dot.i.Omega(mat.rho_cp, s, dT/dt)
             + dw_laplace.i.Omega(mat.lam, s, T)
             = - dw_bc_newton.i.Gamma_Left(h_out_dyn.val, T_out_dyn.val, s, T)
               - dw_bc_newton.i.Gamma_Right(in_fixed.h_in, in_fixed.T_in, s, T)
               + dw_volume_lvf.i.Omega(latent_heat_source.val, s)
               + dw_volume_lvf.i.Omega(heat_source.val, s)
             """
}

ics = {
    'ic': ('Omega', {'T.0': 0.0}),
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
        'i_max': 10,
        'eps_a': 1e-8,
        'is_linear': True,
    }),
    'ts': ('ts.simple', {
        't0': start,
        't1': stop,
        'dt': dt,
        'n_step': nr_of_steps,
        'verbose': 1,
    }),
}

def save_temperature_results(out, problem, state, extend=False):
    """Save temperature values at full-hour intervals."""
    filename = os.path.join(problem.conf.options['output_dir'], "temperature_results_step1.csv")
    
    T_var = state.get_state()
    coors = problem.fields['temperature'].get_coor()

    if T_var.shape[0] != len(coors):  
        raise ValueError(f"Mismatch between state vector size ({T_var.shape[0]}) and number of nodes ({len(coors)})")

    if problem.ts.step == 0:
        with open(filename, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["Time (s)", "Node Index", "Position (m)", "Temperature (°C)"])

    if problem.ts.time % 3600 < problem.ts.dt:  
        with open(filename, 'a', newline='') as csvfile:
            writer = csv.writer(csvfile)
            for i, coord in enumerate(coors):
                writer.writerow([problem.ts.time, i, coord[0], T_var[i]])

    return out

options = {
    'nls': 'newton',
    'ls': 'ls',
    'ts': 'ts',
    'save_times': [3600*i for i in range(1, nr_of_hours + 1)],
    'post_process_hook': save_temperature_results,
    'output_dir': './output_snow_step1',
}