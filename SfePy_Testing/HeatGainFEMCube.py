from sfepy import data_dir
import numpy as nm

filename_mesh = data_dir + '/meshes/3d/cube_medium_tetra.mesh'

# Time settings.
t0 = 0.0
t1 = 10.0  # Final time.
n_step = 100  # Number of time steps.

# Material properties.
material_1 = {
    'name': 'coef',
    'values': {'val': 0.01},  # Thermal conductivity.
    'kind': 'stationary',  # Thermal properties are time-independent.
}

# Field definition for temperature.
field_1 = {
    'name': 'temperature',
    'dtype': 'real',
    'shape': (1,),
    'region': 'Omega',
    'approx_order': 1,
}

# Unknown and test variables.
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

# Regions in the mesh.
regions = {
    'Omega': 'all',
    'Left': ('vertices in (x < -0.499)', 'facet'),
    'Right': ('vertices in (x > 0.499)', 'facet'),
    'Bottom': ('vertices in (z < -0.499)', 'facet'),
    'Top': ('vertices in (z > 0.499)', 'facet'),
}

# Updated Boundary Conditions reflecting soil temperatures.
soil_temp = 5.0  # Mid-range soil temperature between 0 and 8 degrees Celsius
ebcs = {
    'fixed_left': ('Left', {'T.all': soil_temp}),
    'fixed_right': ('Right', {'T.all': soil_temp}),
    'fixed_bottom': ('Bottom', {'T.all': soil_temp}),
}

# Robin boundary condition for the top face.
ambientTemp = 4.0
heatTransferCoef = 10.0
heatFlux = 300.0
def robin(ts, coor, bc=None, problem=None, **kwargs):
    T_inf = kwargs.get('T_inf', ambientTemp)  # Ambient temperature in Celsius
    h = kwargs.get('h', heatTransferCoef)  # Heat transfer coefficient
    q = kwargs.get('q', heatFlux)  # Heat flux due to solar radiation
    return h * (coor[:, 2] - T_inf) + q

bcs = {
    'top': ('Top', {'grad_t': (robin, {'T_inf': ambientTemp, 'h': heatTransferCoef, 'q': heatFlux})}),
}

# Initial condition in Celsius.
def get_ic(coor, ic):
    """Set initial temperature to -5°C everywhere."""
    return nm.full(coor.shape[0], -5.0)

functions = {
    'get_ic' : (get_ic,),
    'robin': (robin,),
}

ics = {
    'ic': ('Omega', {'T.0': 'get_ic'}),
}

# Integral definition for weak formulation.
integral_1 = {
    'name': 'i',
    'order': 1,
}

# PDE definition.
equations = {
    'HeatTransfer':
    """dw_dot.i.Omega(s, dT/dt)
     + dw_laplace.i.Omega(coef.val, s, T) = 0"""
}

# Linear and nonlinear solvers.
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
    'lin_red': 1e-2,
    'is_linear': True,
}

# Time-stepping solver.
solver_2 = {
    'name': 'ts',
    'kind': 'ts.simple',
    't0': t0,
    't1': t1,
    'n_step': n_step,
    'verbose': 1,
}

# Solver options.
options = {
    'nls': 'newton',
    'ls': 'ls',
    'ts': 'ts',
    'save_times': 'all',
    'output_dir': './output',  # Adjust as needed.
}
