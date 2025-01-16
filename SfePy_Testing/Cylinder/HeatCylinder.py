from __future__ import absolute_import
from sfepy import data_dir
import numpy as nm
import os

filename_mesh = data_dir + '/meshes/3d/cylinder.mesh'

# Time settings
t0 = 0.0        # Start time (s)
t1 = 10.0       # End time (s)
n_step = 100    # Number of time steps

# Define materials
materials = {'m': ({'D': 0.0001},)}

# Define regions for the cylinder
regions = {
    'Omega': 'all',  # Entire domain
    'Bottom': ('vertices in (x < 0.00001)', 'facet'),  # Bottom
    'Top': ('vertices in (x > 0.099999)', 'facet')  # Top
}

# Define fields
fields = {'temperature': ('real', 1, # 1: The number of components of the field. Here, it is a scalar field (1 component).
                          'Omega', 1)} # The approximation order of the field. This specifies the polynomial order used in the finite element method.

# Define variables
variables = {
    'T': ('unknown field', 'temperature', 
            1,  # The initial value or guess for the unknown field. Here, it is set to 1. 
            1), # The history type, which is used for time-dependent problems. This indicates that the variable has a history of one previous time step.
    's': ('test field', 'temperature', 'T')
}

# Define initial conditions
def get_ic(coor, ic):
    return nm.zeros_like(coor[:, 0])

# Define initial conditions
#def get_ic(coor, ic):
#    top_temp    = -2.0
#    bottom_temp =  2.0
#    x = -coor[:, 0] # neg bc x-axis points up; "0" - x ("1" - y; "2" - z)
#    x_min, x_max = x.min(), x.max()
#    gradient = (top_temp - bottom_temp) / (x_max - x_min)
#    return top_temp + gradient * (x - x_min)

# External temperature function
def external_temperature(t):
    # Simulate daily temperature cycle
    T_sun = 20.0 + 10.0 * nm.sin(2 * nm.pi * t / 24.0)  # Radiation from the sun
    T_sensible = 5.0 * nm.sin(2 * nm.pi * t / 24.0)  # Sensible heat
    T_latent = 2.0 * nm.sin(2 * nm.pi * t / 24.0)  # Latent heat
    T_rain = 1.0 * nm.sin(2 * nm.pi * t / 24.0)  # Heat from the rain
    return T_sun + T_sensible + T_latent + T_rain

# Robin boundary condition function
def robin_bc(ts, coor, bc, problem, equations=None, term=None, mode=None, region=None, ig=None):
    T_ext = external_temperature(ts.time)
    h = 0.1  # Heat transfer coefficient
    k = 0.0001  # Thermal conductivity of the material
    T = problem.get_variables()['T'].get_state_in_region(region)
    return h * (T - T_ext) / k


# Robin boundary conditions
bcs = {
    'robin_top': ('Top', {'T.0': robin_bc}),
}

# Essential boundary conditions (Dirichlet)
ebcs = {
    'fixed_bottom': ('Bottom', {'T.0': 2.0}),
    'fixed_top': ('Top', {'T.0': -2.0}),
}

# Define integrals
integrals = {'i': 2}

# Define equations
equations = {
    'transient-diffusion': """
        dw_dot.i.Omega(s, dT/dt) +
        dw_laplace.i.Omega(m.D, s, T) = 0
    """
}

# Solvers
solvers = {
    'ls': ('ls.scipy_direct', {}),
    'newton': ('nls.newton', {
        'i_max': 10,
        'eps_a': 1e-10,
    }),
    'ts': ('ts.simple', {
        't0': t0,
        't1': t1,
        'n_step': n_step,
    }),
}

# Options
options = {
    'nls': 'newton',
    'ls': 'ls',
    'ts': 'ts',
    'save_times': 'all',
    'output_dir': './output',
}

# Ensure output directory exists
os.makedirs(options['output_dir'], exist_ok=True)

# Initial conditions
functions = {'get_ic': (get_ic,), 'robin': (robin_bc,)}
ics = {'ic': ('Omega', {'T.0': 'get_ic'})}