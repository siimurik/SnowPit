r"""
Transient advection-diffusion equation with divergence-free advection velocity.

Find :math:`T` such that:

.. math::
    \int_{\Omega} s \pdiff{T}{t}
    + \int_{\Omega} s \nabla \cdot \left(\ul{v} T \right)
    + \int_{\Omega} D \nabla s \cdot \nabla T
    = 0
    \;, \quad \forall s \;.

View the results using::

  sfepy-view cube_medium_hexa.*.vtk -f T:wu 1:vw
"""

from __future__ import absolute_import
from sfepy import data_dir
import numpy as nm

#filename_mesh = data_dir + '/meshes/3d/cube_medium_hexa.mesh'
filename_mesh = "C:\\Users\\sipuga\\Documents\\SnowStorageSolvers\\SfePy_Testing\\SnowPit\\cubeTetra.mesh"

# Time settings
t0 = 0.0
t1 = 10.0
n_step = 100

# Define materials
materials = {
    'm': ({'D': 1.0, 'v': [[0.0], [0.0], [0.5]]},),
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

# Define essential boundary conditions (Dirichlet)
soil_temp = 0.0  # °C
bottom_temp = 4.0  # °C

ebcs = {
    'T_left': ('Left', {'T.0': soil_temp}),
    'T_right': ('Right', {'T.0': soil_temp}),
    'T_front': ('Front', {'T.0': soil_temp}),
    'T_back': ('Back', {'T.0': soil_temp}),
    'T_bottom': ('Bottom', {'T.0': bottom_temp}),
}

# Define Robin boundary condition on the top face
robin_bc = {
    'T_top': ('Top', {'T.0': 'get_robin_value'}),
}

# Define functions for Robin boundary condition
# heat_source_value = 10.0  # Heat source value for Robin BC
# def get_robin_value(ts, coors, **kwargs):
#     return heat_source_value + 10.0 * ts.time  # Example: linear increase with time
# 
# # Define Robin boundary condition on the top face
# robin_bc = {
#     'T_top': ('Top', {'T.0': 'get_robin_value'}),
# }
# 
# functions = {
#     'get_robin_value': (get_robin_value,),
# }


# Robin boundary condition for the top face.
ambient_temp = 4.0
heat_transfer_coef = 10.0
heat_flux = 300.0/273.15

def robin(ts, coor, bc=None, problem=None, **kwargs):
    T_inf = kwargs.get('T_inf', ambient_temp)  # Ambient temperature
    h = kwargs.get('h', heat_transfer_coef)  # Heat transfer coefficient
    q = kwargs.get('q', heat_flux)  # Solar radiation heat flux
    return h * (coor[:, 2] - T_inf) + q  # Returns the values required for

bcs = {
    'top': ('Top', {'grad_t': (robin, {'T_inf': ambient_temp, 
                                       'h': heat_transfer_coef, 
                                       'q': heat_flux
    })}),
}

# Initial conditions.
functions = {
    'robin': (robin,),
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