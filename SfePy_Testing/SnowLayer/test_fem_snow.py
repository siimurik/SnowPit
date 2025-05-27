r"""
Transient 1D heat equation for snow layer with time-dependent BCs.
"""
from __future__ import absolute_import
import numpy as nm
from sfepy.discrete.fem import Mesh
from sfepy.discrete.fem.meshio import UserMeshIO
import csv

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
dt = 10.0  # Time step in seconds
#nh = int(3600 / dt)  # Number of time steps per hour
nh = 1
# Adjust time parameters to focus on the first 10 hours
hours_to_simulate = 10
total_time = hours_to_simulate * 3600  # 10 hours in seconds
n_step = hours_to_simulate * nh  # Only timesteps for 10 hours


def step_hook(pb, ts, variables):
    """
    Extract temperature at the leftmost and rightmost nodes and store it.
    """
    T_field = pb.get_variables()['T']
    coors = pb.fields['temperature'].get_coor()

    num_nodes = coors.shape[0]  # Get number of nodes in mesh
    print(f"Mesh contains {num_nodes} nodes.")  # Debugging

    left_index = 0  # Always first node
    right_index = num_nodes - 1  # Last valid node

    T_left = T_field.data[left_index]  # Bottom layer temperature
    T_right = T_field.data[right_index]  # Top layer temperature

    filename = "temperature_probe.csv"
    with open(filename, mode='a', newline='') as file:
        writer = csv.writer(file)
        writer.writerow([ts.time, T_left, T_right])

    print(f"Saved timestep {ts.time:.2f}s: Bottom Temp = {T_left}, Top Temp = {T_right}")


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
        't1': total_time,
        'dt': dt,
        'n_step': n_step,
        'verbose': True,
    }),
}

options = {
    'step_hook': 'step_hook',  # Calls the function at every timestep
    'nls': 'newton',
    'ls': 'ls',
    'ts': 'ts',
    'save_times': 'all',
    'output_dir': './output_snow',
}