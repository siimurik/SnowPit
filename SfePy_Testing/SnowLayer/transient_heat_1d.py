r"""
Transient 1D heat equation with time-dependent boundary conditions.

This replicates the functionality of the original transient1D() function.
"""
import numpy as nm
from sfepy.discrete.fem import Mesh
from sfepy.discrete.fem.meshio import UserMeshIO
import csv

# Parameters from your original implementation
d_ins = 0.1  # m
lam_i = 0.32  # W/(mK)
D = 8.33e-8  # m²/s (thermal diffusivity)
h_i = 99.75  # W/m2K (inner heat transfer coefficient)
t_i = 0.0  # °C (inner temperature)

# Read time-dependent boundary conditions from CSV
# Assuming CSV has columns: hour, t_o, h_o
# Read boundary condition data from CSV
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
dt = 10.0  # Time step in seconds
nh = int(3600 / dt)  # Number of time steps per hour
total_time = nr_hour * 3600  # Total simulation time in seconds
n_step = nr_hour * nh  # Total number of time steps

# Mesh parameters
dx = 0.005  # Cell size, m
n_el = int(d_ins / dx)  # Number of elements
nodes = n_el + 1  # Number of nodes

def mesh_hook(mesh, mode):
    """
    Generate the 1D mesh.
    """
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

filename_mesh = UserMeshIO(mesh_hook)

def get_bc_values(ts, coor, mode=None, region_name=None, **kwargs):
    """
    Time-dependent boundary condition values.
    """
    if mode != 'qp': return {}
    
    hour = ts.time / 3600  # Current hour
    hour_idx = min(int(hour), nr_hour - 1)  # Index for current hour
    
    nqp = coor.shape[0]
    val = nm.zeros((nqp, 1, 2), dtype=nm.float64)  # h and T_inf
    
    if region_name == 'Gamma_Left':
        # Left boundary (outer surface)
        val[:, 0, 0] = h_o[hour_idx]  # h value
        val[:, 0, 1] = t_o[hour_idx]   # T_inf value
    elif region_name == 'Gamma_Right':
        # Right boundary (inner surface)
        val[:, 0, 0] = h_i  # h value
        val[:, 0, 1] = t_i   # T_inf value
    
    return {'val': val}

materials = {
    'mat': ({'lam': lam_i, 'rho_cp': lam_i/D},),
    'bc': 'get_bc_values',
}

functions = {
    'get_bc_values': (get_bc_values,),
}

regions = {
    'Omega': 'all',
    'Gamma_Left': ('vertices in (x < 0.00001)', 'facet'),
    'Gamma_Right': ('vertices in (x > %f)' % (d_ins - 1e-5), 'facet'),
}

fields = {
    'temperature': ('real', 1, 'Omega', 1),
}

variables = {
    'T': ('unknown field', 'temperature', 0, 1),  # history=1 for time-dependent
    's': ('test field', 'temperature', 'T'),
}

ebcs = {
}

integrals = {
    'i': 2,
}

equations = {
    'Heat': """dw_dot.i.Omega(mat.rho_cp, s, dT/dt)
             + dw_laplace.i.Omega(mat.lam, s, T)
             = dw_bc_newton.i.Gamma_Left(bc.val, s, T)
             + dw_bc_newton.i.Gamma_Right(bc.val, s, T)"""
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
        't1': total_time,
        'dt': dt,
        'n_step': n_step,
        'verbose': True,
    }),
}

options = {
    'nls': 'newton',
    'ls': 'ls',
    'ts': 'ts',
    'save_times': 'all',
    'output_dir': './output_heat_1d',
}