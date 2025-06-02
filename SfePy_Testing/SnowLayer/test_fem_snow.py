r"""
Transient 1D heat equation for snow layer with time-dependent BCs.
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
rho_wet = rho_dry + moist_cont/100.0*1000 # kg/m^3
D = lam_i/(c_wet * rho_wet) # m^2/s

# Boundary conditions
t_i = 0.0         # Inner temperature [°C]
h_i = 99.75       # Inner heat transfer coefficient [W/m2K]

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

t_o, h_o = read_temp_and_hcoeff_from_csv()
nr_hour = len(t_o)  # Number of hours

# Time parameters
hours_to_simulate = 10 # Simulating for 10 hours
total_time = hours_to_simulate * 3600  # 10 hours in seconds

# Ensure enough timesteps while keeping control over total steps
max_timesteps = 500  # Set a reasonable cap (higher than 50)

dt = 3600.0  # Set dt to 1 hour (3600 seconds)
n_step = hours_to_simulate  # Set number of steps equal to number of hours

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
    """
    Time-dependent heat transfer coefficient values (h_o).
    """
    if mode != 'qp' or coors is None:
        return {}

    hour_idx = min(int(ts.time / 3600), len(h_o) - 1)  # Ensure valid index
    val = nm.full((coors.shape[0], 1, 1), h_o[hour_idx], dtype=nm.float64)  # Correct shape

    return {'val': val}

def get_t_o(ts, coors, mode=None, **kwargs):
    """
    Time-dependent ambient temperature values (t_o).
    """
    if mode != 'qp' or coors is None:
        return {}

    hour_idx = min(int(ts.time / 3600), len(t_o) - 1)  # Ensure valid index
    val = nm.full((coors.shape[0], 1, 1), t_o[hour_idx], dtype=nm.float64)  # Correct shape

    return {'val': val}


filename_mesh = UserMeshIO(mesh_hook)

materials = {
    'mat': ({'lam': lam_i, 'rho_cp': rho_wet * c_wet},),
    'heat_loss_h': 'get_h_o',
    'heat_loss_T': 'get_t_o',
}


functions = {
    'get_h_o': (get_h_o,),
    'get_t_o': (get_t_o,),
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

# Remove EBCs if using Newton BCs
ebcs = {
    't2': ('Gamma_Right', {'T.0': 0.0}),  # Constant temperature on right
}

integrals = {
    'i': 2,
}

equations = {
    'Heat': """dw_dot.i.Omega(mat.rho_cp, s, dT/dt)
             + dw_laplace.i.Omega(mat.lam, s, T)
             = dw_bc_newton.i.Gamma_Left(heat_loss_h.val, heat_loss_T.val, s, T)"""
}



ics = {
    'ic': ('Omega', {'T.0': 0.0}),
}

solvers = {
    'ls': ('ls.scipy_direct', {}),
    'newton': ('nls.newton', {
        'i_max': 1,
        'eps_a': 1e-10,
        'is_linear': True,
    }),
    'ts': ('ts.simple', {
        't0': 0.0,
        't1': total_time,  # 10 hours total
        'dt': dt,  # 1 hour step size
        'n_step': n_step,  # 10 total steps
        'verbose': True,
    }),
}


def save_temperature_results(out, problem, state, extend=False):
    """ Save temperature values at each timestep into a CSV file. """
    filename = os.path.join(problem.conf.options['output_dir'], "temperature_results.csv")
    
    # Get the full state vector instead of trying to index it as a dictionary
    T_var = state.get_state()  # Returns the full solution vector
    coors = problem.fields['temperature'].get_coor()  # Get mesh node coordinates

    # Ensure correct indexing based on the shape of T_var
    if T_var.shape[0] != len(coors):  # Check if dimensions match
        raise ValueError(f"Mismatch between state vector size ({T_var.shape[0]}) and number of nodes ({len(coors)})")

    # Write headers only for the first timestep
    if problem.ts.step == 0:
        with open(filename, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["Time (s)", "Node Index", "Position (m)", "Temperature (°C)"])

    # Append results for each timestep
    with open(filename, 'a', newline='') as csvfile:
        writer = csv.writer(csvfile)
        for i, coord in enumerate(coors):
            writer.writerow([problem.ts.time, i, coord[0], T_var[i]])

    return out  # Ensure compatibility with Sfepy execution


options = {
    'nls': 'newton',
    'ls': 'ls',
    'ts': 'ts',
    'save_times': 'all',
    'post_process_hook': save_temperature_results,  # Runs after each timestep
    'output_dir': './output_snow',
}