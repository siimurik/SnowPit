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

t0 = 0.0
t1 = 0.1
n_step = 101

h = 10.0  # W/m2/K
T0 = -2.0 # °C
materials = {
    'flux' : ({'val' : 25.0},),
    'm': ({'D': 0.01,},),  
    'heat_loss': ({ 'h_bot': -h, 'T_bot_inf': T0,
                    'h_top': -h, 'T_top_inf': T0},)
}

field_1 = {
    'name': 'temperature',
    'dtype': 'real',
    'shape': (1,),
    'region': 'Omega',
    'approx_order': 1,
}

variable_1 = {
    'name': 'T',
    'kind': 'unknown field',
    'field': 'temperature',
    'order': 0,
    'history': 1,
}

variable_2 = {
    'name': 's',
    'kind': 'test field',
    'field': 'temperature',
    'dual': 'T',  
}

# Define regions for the frustum pyramid
regions = {
    'Omega': 'all',  # Entire domain
    'Left': ('vertices in (x < -1.999)', 'facet'),  # Left boundary
    'Right': ('vertices in (x > 1.999)', 'facet'),  # Right boundary
    'Front': ('vertices in (y > 3.999)', 'facet'),  # Front boundary
    'Back': ('vertices in (y < -3.999)', 'facet'),  # Back boundary
    'Bottom': ('vertices in (z < 0.001)', 'facet'),  # Bottom boundary
    'Top': ('vertices in (z > 3.999)', 'facet'),  # Top boundary
}

ebcs = {
    'T1': ('Bottom', {'T.0': 2.0}),
}

ics = {
    'ic': ('Omega', {'T.0': T0}),
}

integrals = {
    'i' : 2,
}

equations = {
    'Temperature': """
    dw_dot.i.Omega( s, dT/dt ) + dw_laplace.i.Omega( m.D, s, T ) =
    + dw_integrate.i.Top(flux.val, s)
    + dw_bc_newton.i.Top(heat_loss.h_top, heat_loss.T_top_inf, s, T)
    """
}

solver_0 = {
    'name': 'ls',
    'kind': 'ls.scipy_direct',
    'use_presolve': True,
}

solver_1 = {
    'name': 'newton',
    'kind': 'nls.newton',

    'i_max': 1,
    'eps_a': 1e-10,
    'eps_r': 1.0,
    'macheps': 1e-16,
    'lin_red': 1e-2,  # Linear system error < (eps_a * lin_red).
    'ls_red': 0.1,
    'ls_red_warp': 0.001,
    'ls_on': 1.1,
    'ls_min': 1e-5,
    'check': 0,
    'delta': 1e-6,
    'is_linear': True,
}

solver_2 = {
    'name': 'ts',
    'kind': 'ts.simple',
    't0': t0,
    't1': t1,
    'dt': None,
    'n_step': n_step,  # has precedence over dt!
    'verbose': 1,
}

options = {
    'nls': 'newton',
    'ls': 'ls',
    'ts': 'ts',
    'save_times': 'all',
    'output_dir': './output_pyramid',  # Directory to save the output
}
