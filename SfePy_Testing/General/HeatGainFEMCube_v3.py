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

# Accompanying mesh from SfePy
filename_mesh = data_dir + '/meshes/3d/cube_medium_hexa.mesh'

# Homemade mesh 
#filename_mesh = "C:\\Users\\sipuga\\Documents\\SnowStorageSolvers\\SfePy_Testing\\SnowPit\\cubeTetra.mesh"

# Time settings
t0 = 0.0        # Start time (s)
t1 = 10.0       # End time (s)
n_step = 100    # Number of time steps

# Define materials
materials = {
    'm': ({'D': 0.01, # Diffusion coefficient (m²/s)
           'v': [[0.0], [0.0], [0.0035]]},),    # convection speed (m/s)
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
    'temperature': ('real', 1,   # The number of components (1 means scalar field, e.g., temperature).
                    'Omega', 1), # The approximation order (linear approximation).
}

# Define variables
variables = {
    'T': ('unknown field', 'temperature', 
          0,    # The time derivative order (0 means a steady-state problem; 1 would mean a transient problem).
          1),   # The approximation order (consistent with the field’s order).
    's': # The test field (used in the weak formulation of the problem). 
         # It represents the same temperature field and is paired with 'T'.
         ('test field', 'temperature', 'T'),
}

# Define temperatures based on depth.
def soil_temp(depth):
    surface_temp = 0.0  # Surface temperature (°C)
    temp_gradient = 4.0  # Temperature difference (°C/m)
    return surface_temp + temp_gradient * depth

# Apply soil temperature to initial condition.
def get_ic(coor, ic):
    depth = -coor[:, 2]  # Assuming z-axis points upwards
    return soil_temp(depth)

# Robin boundary condition for the top face.
ambient_temp = 4000.0  # Ambient temperature (°C)
heat_transfer_coef = 10000.0  # Base heat transfer coefficient (W/m²·°C)
heat_flux = 30000.0  # Solar radiation heat flux (W/m²)
layer_thickness = 0.000005  # Thickness of the insulating layer (m)
layer_conductivity = 11111.1  # Thermal conductivity of the insulating layer (W/m·°C)

def robin(ts, coor, bc=None, problem=None, **kwargs):
    T_inf = kwargs.get('T_inf', ambient_temp)  # Ambient temperature
    h_base = kwargs.get('h', heat_transfer_coef)  # Base heat transfer coefficient
    q = kwargs.get('q', heat_flux)  # Solar radiation heat flux
    
    # Effective heat transfer coefficient accounting for insulation
    d = kwargs.get('thickness', layer_thickness)  # Thickness of the layer
    k = kwargs.get('conductivity', layer_conductivity)  # Thermal conductivity
    h_eff = k / (d + k / h_base)
    
    # Robin condition: q + h_eff * (T - T_inf)
    return h_eff * (coor[:, 2] - T_inf) + q


bcs = {
    'top': ('Top', {'grad_t': (robin, {'T_inf': ambient_temp, 
                                       'h': heat_transfer_coef, 
                                       'q': heat_flux
    })}),
}

#bcs = {
#    'robin': ('Top', {'T': 'robin'}),
#}


### Define essential boundary conditions (Dirichlet)
# Set walls at a fixed value (False)
add_gradient = True    # <--- add temp gradient to sides if set to 'True'
if add_gradient == False:
    soil_temp = 0.0     # °C
    bottom_temp = 4.0   # °C
    ebcs = {    # essential boundary conditions
    'T_left': ('Left', {'T.0': soil_temp}),
    'T_right': ('Right', {'T.0': soil_temp}),
    'T_front': ('Front', {'T.0': soil_temp}),
    'T_back': ('Back', {'T.0': soil_temp}),
    #'T_top': ('Top', {'T.0': bottom_temp}),
    'T_bottom': ('Bottom', {'T.0': bottom_temp}),
    }
    # Initial conditions.
    functions = {
    'robin': (robin,),
    }
else:   # add gradient if True
    # Example values for materials: 'm': ({'D': 0.01, 'v': [[0.0], [0.0], [0.035]]},)
    # Boundary conditions.
    ground_temp = 4.0   # °C
    ebcs = {
        #'fixed_left': ('Left', {'T.0': 5.0}),
        #'fixed_right': ('Right', {'T.0': 5.0}),
        'fixed_bottom': ('Bottom', {'T.0': ground_temp}),
    }
    # Initial conditions.
    functions = {
        'get_ic': (get_ic,),
        'robin': (robin,),
    }
    ics = {
    'ic': ('Omega', {'T.0': 'get_ic'}),
    }
# End of if block

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

