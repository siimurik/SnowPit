r"""
==================================================================
Transient 1D heat equation for snow layer with time-dependent BCs.
------------------------------------------------------------------
Solves the 1D transient heat equation:
    ∂T/∂t = α ∂²T/∂x²
with boundary conditions:
    At x=0: -k∂T/∂x = h_o(T - T_o)
    At x=L: -k∂T/∂x = h_i(T - T_i)
where:
    α = k/(ρc_p) is thermal diffusivity
------------------------------------------------------------------
For uploading to Github:
    git add test_fem_snow.py 
    git commit -m "text here"
    git push origin main
------------------------------------------------------------------
Run this script with:
    sfepy-run transient_heat_1d.py
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
rho_wet = rho_dry + moist_cont*10.0 # /100.0*1000 # kg/m^3
D = lam_i/(c_wet * rho_wet) # m^2/s

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
    Retrieves the time-dependent heat transfer coefficient (`h_o`) for a given time step.

    This function returns a uniform heat transfer coefficient at full-hour intervals.
    The same value remains applied throughout the hour until the next update.

    Parameters:
        ts (TimeStepper): Time step object containing the current simulation time in seconds.
        coors (array): Coordinates of mesh points (required for function compatibility).
        mode (str, optional): Specifies the evaluation mode (must be 'qp' for query points).
        **kwargs: Additional arguments (not used here).

    Returns:
        dict: A dictionary with the key `'val'` containing an array of heat transfer coefficient values
              that remain constant for the duration of each hour.

    Example:
        - At `ts.time = 0s` → Uses `h_o[0]`
        - At `ts.time = 3600s` → Uses `h_o[1]`
        - At `ts.time = 5400s` → Still uses `h_o[1]` (until next full hour)
    """
    if mode != 'qp' or coors is None:
        return {}

    hour_idx = min(int(ts.time / 3600), len(h_o) - 1)  # Select correct hourly index
    val = nm.full((coors.shape[0], 1, 1), h_o[hour_idx], dtype=nm.float64)  # Apply constant value

    return {'val': val}

def get_t_o(ts, coors, mode=None, **kwargs):
    """
    Retrieves the time-dependent ambient temperature (`t_o`) for a given time step.

    Similar to `get_h_o()`, this function applies ambient temperature values at 
    full-hour marks. The temperature remains unchanged until the next hour.

    Parameters:
        ts (TimeStepper): Time step object containing the current simulation time in seconds.
        coors (array): Coordinates of mesh points (required for function compatibility).
        mode (str, optional): Specifies the evaluation mode (must be 'qp' for query points).
        **kwargs: Additional arguments (not used here).

    Returns:
        dict: A dictionary with the key `'val'` containing an array of ambient temperature values
              that remain constant for the duration of each hour.

    Example:
        - At `ts.time = 0s` → Uses `t_o[0]`
        - At `ts.time = 3600s` → Uses `t_o[1]`
        - At `ts.time = 5400s` → Still uses `t_o[1]` (until next full hour)
    """
    if mode != 'qp' or coors is None:
        return {}

    hour_idx = min(int(ts.time / 3600), len(t_o) - 1)  # Select correct hourly index
    val = nm.full((coors.shape[0], 1, 1), t_o[hour_idx], dtype=nm.float64)  # Apply constant value

    return {'val': val}


###############################################################################
#                                MAIN SECTION                                 #
###############################################################################

# Read data and calculate boundary conditions
airTemp, airVel, prec, gloSolIr, rh = read_input_data()

# Boundary conditions
t_i = 0.0 # Inner layer temp
h_i = 99.75 # Inner layer heat transfer coefficient 

# Calculation of equivalent sol-air temperature (in °C)
h = 22.7; alpha = 0.8; T_cor_fact = 4.0
t_o = [alpha * glob_solir / h + air_temp - T_cor_fact for 
       glob_solir, air_temp in zip(gloSolIr, airTemp)]

# Calculation of convective heat transfer coefficient 
h_o = [6.0 + 4.0*vel if vel <= 5.0 else 7.41*(vel**0.78) for vel in airVel]

# Mesh generation
filename_mesh = UserMeshIO(mesh_hook)

materials = {
    'mat': ({'lam': lam_i, 'rho_cp': rho_wet * c_wet},),
    'h_out_dyn': 'get_h_o',  # Left boundary (variable)
    'T_out_dyn': 'get_t_o',  # Left boundary (variable)
    'in_fixed': ({'h_in': h_i, 'T_in': t_i},),  # Right boundary (fixed values)
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
ebcs = {}
#ebcs = {
#    't2': ('Gamma_Right', {'T.0': 0.0}),  # Constant temperature on right
#}

integrals = {
    'i': 1,
}

equations = {
    'Heat': """dw_dot.i.Omega(mat.rho_cp, s, dT/dt)
             + dw_laplace.i.Omega(mat.lam, s, T)
             = - dw_bc_newton.i.Gamma_Left(h_out_dyn.val, T_out_dyn.val, s, T)
               - dw_bc_newton.i.Gamma_Right(in_fixed.h_in, in_fixed.T_in, s, T)
             """
}


ics = {
    'ic': ('Omega', {'T.0': 0.0}),
}

# Time parameters - match original function's resolution
nr_hour = len(t_o)  # Number of total hours
start = 0.0
nr_of_hours = 10
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
        't1': stop,             # Hours in seconds
        'dt': dt,               # Step size of 10 seconds
        'n_step': nr_of_steps,  # Ensure 360 steps are taken for every hour # has precedence over dt!
        'verbose': 1,
    }),
}


def save_temperature_results(out, problem, state, extend=False):
    """ Save temperature values only at full-hour intervals (every 3600s) into a CSV file. """
    filename = os.path.join(problem.conf.options['output_dir'], "temperature_results.csv")
    
    # Get the full state vector
    T_var = state.get_state()  # Returns the full solution vector
    coors = problem.fields['temperature'].get_coor()  # Get mesh node coordinates

    # Ensure correct indexing based on the shape of T_var
    if T_var.shape[0] != len(coors):  
        raise ValueError(f"Mismatch between state vector size ({T_var.shape[0]}) and number of nodes ({len(coors)})")

    # Write headers only for the first timestep
    if problem.ts.step == 0:
        with open(filename, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["Time (s)", "Node Index", "Position (m)", "Temperature (°C)"])

    # **Only save data at full-hour marks (3600s, 7200s, ...)**
    if problem.ts.time % 3600 < problem.ts.dt:  
        with open(filename, 'a', newline='') as csvfile:
            writer = csv.writer(csvfile)
            for i, coord in enumerate(coors):
                writer.writerow([problem.ts.time, i, coord[0], T_var[i]])

    return out  # Ensure compatibility with Sfepy execution


options = {
    'nls': 'newton',
    'ls': 'ls',
    'ts': 'ts',
    'save_times': [3600*i for i in range(1, nr_of_hours + 1)],
    'post_process_hook': save_temperature_results,  # Runs after each timestep
    'output_dir': './top_snow_layer_1d',
}