from __future__ import absolute_import
from sfepy.base.base import Struct
import numpy as nm
from sfepy import data_dir

filename_mesh = data_dir + '/meshes/3d/cube_medium_hexa.mesh'

# Define material properties for cold material (water or ice) and the top boundary with temperature fluctuations
material_1 = {
    'name': 'coef',
    'values': {'val': 1.0},
}

# Define the top boundary where the Robin condition will be applied, with temperature fluctuation source term
def temperature_fluctuation_source_term(ts, coors, mode=None, **kwargs):
    if mode == 'qp':  # apply to quadrature points
        nqp, dim = coors.shape
        amplitude = 10.0  # example amplitude of temperature fluctuation in degrees Celsius
        period = 24 * 60 * 60  # period of fluctuation (24 hours in seconds)
        time = ts.time
        fluctuation = amplitude * nm.sin(2 * nm.pi * time / period)  # simple sinusoidal fluctuation
        source = fluctuation * nm.ones((nqp, 1), dtype=nm.float64)
        return {'source_term': source}
    return {}

# Define the regions for the cube (for boundaries and volume)
regions = {
    'Omega': 'cells of group 0',  # interior of the cube
    'Top': 'cells in (z > 0.4999999)',  # the top surface of the cube
    'Bottom': 'cells in (z < -0.4999999)',  # bottom surface
    'Surface': 'vertices of surface',  # all surfaces
}

# Define the field for temperature (scalar field)
fields = {
    'temperature': ('real', 1, 'Omega', 1),  # temperature field
}

# Define the variables for temperature
variables = {
    'T': ('unknown field', 'temperature'),  # unknown temperature field
    's': ('test field', 'temperature', 'T'),  # test field for Galerkin method
}

# Boundary conditions (Dirichlet for most sides, Robin for the top)
ebcs = {
    'T0': ('Surface', {'T.0': -3.0}),  # bottom surface fixed temperature
    'T1': ('Top', {'T.0': 1.0}),  # Dirichlet for the top
    'T2': ('Bottom', {'T.0': -1.0}),  # Dirichlet for the bottom
    'T3': ('Left', {'T.0': 2.0}),  # Dirichlet for the left face
}

# Equation for Laplace problem with the source term for the top boundary
equations = {
    'temperature_equation': """dw_laplace.i.Omega( coef.val, s, T )
                               + dw_integrate.5.Top(source_term, s) = 0""",
}

# Set up the solvers and solver options
solvers = {
    'ls': ('ls.scipy_direct', {}),
    'newton': ('nls.newton',
               {'i_max': 1,
                'eps_a': 1e-6,
                'eps_r': 1.0}),
}

# Define the material function to include source terms
functions = {
    'temperature_fluctuation_source_term': (temperature_fluctuation_source_term,),
}

# Solver setup
options = {
    'post_process_hook': 'postproc',
}

def postproc(out, pb, state, extend=False):
    # Extract the temperature fluctuation results for visualization
    fluctuation = pb.evaluate('ev_integrate_mat.5.Top(source_term, T)', mode='el_avg')
    out['fluctuation'] = Struct(name='output_data',
                                mode='cell',
                                data=fluctuation.reshape(fluctuation.shape[0], 1, 1, 1),
                                dofs=None)
    return out
