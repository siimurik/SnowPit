from __future__ import absolute_import
from sfepy import data_dir

# Mesh filename
filename_mesh = data_dir + '/meshes/3d/cube_medium_hexa.mesh'

############# Heat Equation (Conduction and Convection).

# Material properties (e.g., thermal conductivity for water or ice)
material_1 = {
    'name' : 'thermal_conductivity',
    'values' : {'val' : 0.6},  # Example for ice, adjust as needed for your material
}

# Define the temperature field
field_1 = {
    'name' : 'temperature',
    'dtype' : 'real',
    'shape' : (1,),
    'region' : 'Omega',
    'approx_order' : 1,
}

# Region for the cube
region_1000 = {
    'name' : 'Omega',
    'select' : 'cells of group 0',
}

# Integral order (first-order quadrature)
integral_1 = {
    'name' : 'i',
    'order' : 1,
}

# Solver options
solver_0 = {
    'name' : 'ls',
    'kind' : 'ls.scipy_direct',
}

# Variables (temperature and test functions)
variable_1 = {
    'name' : 'T',
    'kind' : 'unknown field',
    'field' : 'temperature',
    'order' : 0,
}

variable_2 = {
    'name' : 's',
    'kind' : 'test field',
    'field' : 'temperature',
    'dual' : 'T',
}

# Define the regions for boundary conditions
region_0 = {
    'name' : 'Surface',
    'select' : 'vertices of surface',
    'kind' : 'facet',
}

region_1 = {
    'name' : 'Bottom',
    'select' : 'vertices in (z < -0.4999999)',
    'kind' : 'facet',
}

region_2 = {
    'name' : 'Top',
    'select' : 'vertices in (z > 0.4999999)',
    'kind' : 'facet',
}

region_03 = {
    'name' : 'Left',
    'select' : 'vertices in (x < -0.4999999)',
    'kind' : 'facet',
}

# Dirichlet Boundary Conditions for 5 faces (constant temperature)
ebc_1 = {
    'name' : 'T0',
    'region' : 'Surface',
    'dofs' : {'T.0' : -3.0},  # Example temperature
}

ebc_4 = {
    'name' : 'T1',
    'region' : 'Top',
    'dofs' : {'T.0' : 1.0},  # Example temperature
}

ebc_3 = {
    'name' : 'T2',
    'region' : 'Bottom',
    'dofs' : {'T.0' : -1.0},  # Example temperature
}

ebc_2 = {
    'name' : 'T3',
    'region' : 'Left',
    'dofs' : {'T.0' : 2.0},  # Example temperature
}

# Robin Boundary Condition for the top face
robin_bc = {
    'name' : 'robin_bc',
    'region' : 'Top',
    'dofs' : {'T.0' : (0.0, 5.0)},  # Example for Robin condition (thermal resistance)
}

# Define the equation for the heat equation (convection and conduction)
equations = {
    'heat_equation' : """dw_laplace.i.Omega( thermal_conductivity.val, s, T ) = 0""",
}

# Solver for the system
solver_1 = {
    'name'      : 'newton',
    'kind'      : 'nls.newton',
    'i_max'     : 1,
    'eps_a'     : 1e-10,
    'eps_r'     : 1.0,
    'macheps'   : 1e-16,
    'lin_red'   : 1e-2,
    'ls_red'    : 0.1,
    'check'     : 0,
}

# To use the material properties, thermal boundary conditions, and solve
