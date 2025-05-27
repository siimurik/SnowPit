"""
FEM is very sensitive to boundary conditions. If bottom is surrounded by cold 
temperature, heat won't progress at all. The same goes for the top which is 
surrounded by cold temperatures.
"""
from __future__ import absolute_import
from sfepy import data_dir
import numpy as nm

filename_mesh = data_dir + '/meshes/3d/cube_medium_hexa.mesh'

# Time settings.
t0 = 0.0
t1 = 10.0  # Final time.
n_step = 100  # Number of time steps.

# Define materials
materials = {
    'm': ({'D': 0.01, 'v': [[0.1], [0.0], [0.0]]},),
}

# Define regions
regions = {
    'Omega': 'all',
    'Left': ('vertices in (x < -0.499)', 'facet'),
    'Right': ('vertices in (x > 0.499)', 'facet'),
    'Front': ('vertices in (y > 0.499)', 'facet'),
    'Back': ('vertices in (y < -0.499)', 'facet'),
    'Bottom': ('vertices in (z < -0.499)', 'facet'),
    'Top': ('vertices in (z > 0.499)', 'facet'),
}

# Define fields
fields = {
    'temperature': ('real', 1, 'Omega', 1),
}

# Define variables
variables = {
    'T': ('unknown field', 'temperature', 0, 1),
    's': ('test field', 'temperature', 'T'),
}


# Define temperatures based on depth.
def soil_temp(depth):
    surface_temp = 0.0  # Surface temperature in Celsius
    temp_gradient = 4.0  # Temperature difference in Celsius per unit depth
    return surface_temp + temp_gradient * depth

# Apply soil temperature to initial condition.
def get_ic(coor, ic):
    depth = -coor[:, 2]  # Assuming z-axis points upwards
    return soil_temp(depth)

# Robin boundary condition for the top face.
ambient_temp = 4.0
heat_transfer_coef = 10.0
heat_flux = 300.0/273.15

def robin(ts, coor, bc=None, problem=None, **kwargs):
    T_inf = kwargs.get('T_inf', ambient_temp)  # Ambient temperature
    h = kwargs.get('h', heat_transfer_coef)  # Heat transfer coefficient
    q = kwargs.get('q', heat_flux)  # Solar radiation heat flux
    return h * (coor[:, 2] - T_inf) + q  # Returns the values required for the LCBC

# Boundary conditions.
ground_temp = 4.0
ebcs = {
    #'fixed_left': ('Left', {'T.0': 5.0}),
    #'fixed_right': ('Right', {'T.0': 5.0}),
    'fixed_bottom': ('Bottom', {'T.0': ground_temp}),
}

bcs = {
    'top': ('Top', {'grad_t': (robin, {'T_inf': ambient_temp, 
                                       'h': heat_transfer_coef, 
                                       'q': heat_flux
    })}),
}


# Initial conditions.
functions = {
    'get_ic': (get_ic,),
    'robin': (robin,),
}

ics = {
    'ic': ('Omega', {'T.0': 'get_ic'}),
}

# Define integrals
integrals = {
    'i': 3,
}

# Define equations
equations = {
    'advection-diffusion': """
        dw_dot.i.Omega(s, dT/dt)
      + dw_advect_div_free.i.Omega(m.v, s, T)
      + dw_laplace.i.Omega(m.D, s, T)
      = 0
    """
}


# Solvers.
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


# Options.
options = {
    'nls': 'newton',
    'ls': 'ls',
    'ts': 'ts',
    'save_times': 'all',
    'output_dir': './output',
}