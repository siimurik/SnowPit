from sfepy import data_dir
import numpy as nm

#filename_mesh = data_dir + '/meshes/3d/cube_big_tetra.mesh'
filename_mesh = data_dir + '/meshes/3d/cube_medium_tetra.mesh'


# Time settings.
t0 = 0.0
t1 = 10.0
n_step = 200

# Material properties. 
material_1 = {
    'name': 'coef',
    'values': {'val': 0.01}, # Thermal conductivity
    'kind': 'stationary',    # Thermal properties are time-independent.
}

material_2 = {
    'name': 'm',
    'values': {'D': 0.1, 
               'v': [[1.0], [0.0], [0.0]]}, # Diffusion coefficient
    'kind': 'stationary', 
}

# List of materials.
materials = [
    material_1,  # Thermal conductivity material.
    material_2   # Advection-diffusion properties.
]

# Fields.
fields = {
    'temperature': ('real', 1, 'Omega', 1),
}

# Variables.
variables = {
    'T': ('unknown field', 'temperature', 0, 1),
    's': ('test field', 'temperature', 'T'),
}

# Regions.
regions = {
    'Omega': 'all',
    'Left': ('vertices in (x < -0.499)', 'facet'),
    'Right': ('vertices in (x > 0.499)', 'facet'),
    'Front': ('vertices in (y > 0.499)', 'facet'),
    'Back': ('vertices in (y < -0.499)', 'facet'),
    'Bottom': ('vertices in (z < -0.499)', 'facet'),
    'Top': ('vertices in (z > 0.499)', 'facet'),
}

# Boundary Conditions.
soil_temp   = -5.0 # °C
bottom_temp = 4.0 # °C
ebcs = {
    'fixed_left': ('Left', {'T.all':   soil_temp}),
    'fixed_right': ('Right', {'T.all': soil_temp}),
    'fixed_front': ('Front', {'T.all': soil_temp}),
    'fixed_back': ('Back', {'T.all':   soil_temp}),
    'fixed_bottom': ('Bottom', {'T.all': bottom_temp}),
}

ambient_temp = 12.0 # °C
heat_transf_coef = 300.0/273.15 # Heat transfer coefficient W/(m^2*°C).
def robin(ts, coor, bc=None, problem=None, **kwargs): 
    T_inf = kwargs.get('T_inf', ambient_temp)  # Ambient temperature °C
    h = kwargs.get('h',heat_transf_coef) 
    return -h * (coor[:, 2] - T_inf)

bcs = {
    'robin_top': ('Top', 
        {'T.robin': (robin, 
            {
                'T_inf': ambient_temp, 
                'h': heat_transf_coef
            })
        }),
}

# Initial Conditions.
def get_ic(coor, ic=None):
    return nm.full(coor.shape[0], soil_temp) # °C

functions = {
    #'get_ic': (get_ic,),
    'robin': (robin,),
}

#ics = {
#    'ic': ('Omega', {'T.0': 'get_ic'}),
#}

# Integral.
integrals = {
    'i': 2,
}

# TODO: Try convection instead of advection
# The positive sign indicates that the convection process is contributing 
# to the transport of heat. It represents the added effect of the bulk fluid 
# movement, which can carry heat along with it, potentially leading to an 
# increase in temperature in some regions.

# Equations.
equations = {
    'advection_diffusion':
    """
       dw_dot.i.Omega(s, dT/dt)
     + dw_advect_div_free.i.Omega(m.v, s, T)
     + dw_laplace.i.Omega(m.D, s, T)    
     = 0
    """
}
# "- dw_laplace.i.Omega(m.D, s, T)": The negative sign is essential 
# because heat flows from higher to lower temperature regions, 
# which aligns with the concept of reducing temperature gradients

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
