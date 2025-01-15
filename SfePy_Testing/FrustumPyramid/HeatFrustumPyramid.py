r"""
Transient diffusion equation:

    \int_{\Omega} s \frac{\partial T}{\partial t} + 
    \int_{\Omega} D \nabla s \cdot \nabla T = 0 \;, 
    \quad \forall s \in H^1(\Omega).

"""

from __future__ import absolute_import
from sfepy import data_dir
import numpy as nm
import os

# Homemade mesh 
filename_mesh = "C:\\Users\\sipuga\\Documents\\SnowStorageSolvers\\SfePy_Testing\\SnowPit\\FrustumPyramid\\frustum_pyramid.mesh"

# Time settings
t0 = 0.0        # Start time (s)
t1 = 10.0       # End time (s)
n_step = 100    # Number of time steps

# Define materials
materials = {
    'm': ({'D': 0.1,  # Diffusion coefficient (m²/s)
    },),    
}

# Define regions for the frustum pyramid
regions = {
    'Omega': 'all',  # Entire domain
    'Left': ('vertices in (x < -2.0)', 'facet'),  # Left boundary
    'Right': ('vertices in (x > 2.0)', 'facet'),  # Right boundary
    'Front': ('vertices in (y > 4.0)', 'facet'),  # Front boundary
    'Back': ('vertices in (y < -4.0)', 'facet'),  # Back boundary
    'Bottom': ('vertices in (z < 0.001)', 'facet'),  # Bottom boundary
    'Top': ('vertices in (z > 3.999)', 'facet'),  # Top boundary
}

# Define fields
fields = {
    'temperature': ('real', 1,   # The number of components (1 means scalar field, e.g., temperature).
                    'Omega', 1), # The approximation order (linear approximation).
}

# Define variables
variables = {
    'T': (  'unknown field', 'temperature', 
            1,    # The time derivative order (0 means a steady-state problem; 1 would mean a transient problem).
            1),   # The approximation order (consistent with the field’s order).
            's': ('test field', 'temperature', 'T'),  # The test field
}

# Define temperatures based on depth.
def soil_temp(depth):
    surface_temp  = 0.0  # Surface temperature (°C)
    temp_gradient = 4.0  # Temperature difference (°C/m)
    return surface_temp + temp_gradient * depth

# Apply soil temperature to initial condition.
def get_ic(coor, ic):
    depth = -coor[:, 2]  # Assuming z-axis points upwards
    return soil_temp(depth)

# Robin boundary condition for the top face.
ambient_temp = 4.0  # Ambient temperature (°C)
heat_transfer_coef = 10.0  # Base heat transfer coefficient (W/m²·°C)
heat_flux = 300.0  # Solar radiation heat flux (W/m²)
layer_thickness = 0.1  # Thickness of the insulating layer (m)
layer_conductivity = 0.1  # Thermal conductivity of the insulating layer (W/m·°C)

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

# Boundary conditions for the top face using Robin condition.
bcs = {
    'top': ('Top', {'grad_t': (robin, {'T_inf': ambient_temp, 
                                       'h': heat_transfer_coef, 
                                       'q': heat_flux,
                                       'thickness': layer_thickness,
                                       'conductivity': layer_conductivity
    })}),
}

# Essential boundary conditions (Dirichlet)
add_gradient = True  # Flag to add temperature gradient to sides if set to 'True'

if add_gradient == False:
    soil_temp = 0.0     # °C
    bottom_temp = 4.0   # °C
    ebcs = {    # essential boundary conditions
        'T_left': ('Left', {'T.0': soil_temp}),
        'T_right': ('Right', {'T.0': soil_temp}),
        'T_front': ('Front', {'T.0': soil_temp}),
        'T_back': ('Back', {'T.0': soil_temp}),
        'T_bottom': ('Bottom', {'T.0': bottom_temp}),
    }
    # Initial conditions.
    functions = {
        'robin': (robin,),
    }
else:   # add gradient if True
    ground_temp = 4.0   # °C
    ebcs = {
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

# Define integrals
integrals = {
    'i': 3,  # Third-order quadrature
}

# Define equations
equations = {
    'transient-diffusion': """
        dw_dot.i.Omega(s, dT/dt)
      + dw_laplace.i.Omega(m.D, s, T)
      = 0
    """
}

# Solvers.
solvers = {
    'ls': ('ls.scipy_direct', {}),  # Linear solver
    'newton': ('nls.newton', {  # Non-linear solver
        'i_max': 10,  # Maximum number of iterations
        'eps_a': 1e-10,  # Absolute tolerance
    }),
    'ts': ('ts.simple', {  # Time-stepping solver
        't0': t0,  # Start time
        't1': t1,  # End time
        'n_step': n_step,  # Number of time steps
    }),
}

# Options.
options = {
    'nls': 'newton',  # Non-linear solver
    'ls': 'ls',  # Linear solver
    'ts': 'ts',  # Time-stepping solver
    'save_times': 'all',  # Save results at all time steps
    'output_dir': './output',  # Directory to save the output
}

# Ensure output directory exists
os.makedirs(options['output_dir'], exist_ok=True)
