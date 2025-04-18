import numpy as nm
from sfepy import data_dir
from sfepy.discrete.probes import LineProbe
import matplotlib.pyplot as plt
import os

filename_mesh = data_dir + '/meshes/3d/cylinder.mesh'

t0 = 0.0
t1 = 0.1
n_step = 101

nominal_heat_flux = 6.36e5
alpha = 0.25
t_start = nm.array([0., 20., 40.])  # times when heating starts (seconds)
t_stop = nm.array([10., 30., 50.])  # times when heating stops (seconds)

h = 10.0  # W/m2/K
T0 = -2.0  # °C
materials = {
    'flux': ({'val': 25.0},),
    'm': ({'D': 0.01,},),
    'heat_loss': ({'h_bot': -h, 'T_bot_inf': T0,
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

regions = {
    'Omega': 'all',
    'Gamma_Left': ('vertices in (x < 0.00001)', 'facet'),
    'Gamma_Right': ('vertices in (x > 0.099999)', 'facet'),
}

ebcs = {
    'T1': ('Gamma_Left', {'T.0': 2.0}),
}

def get_ic(coor, ic):
    return -2.0

def get_flux_value(ts, coors, mode=None, **kwargs):
    """Defines heat flux as a function of time."""
    if mode == 'qp':
        shape = (coors.shape[0], 1, 1)
        if nm.any((ts.time >= t_start) & (ts.time <= t_stop)):
            flux = alpha * nominal_heat_flux
        else:
            flux = 0.
        val = flux * nm.ones(shape, dtype=nm.float64)
        return {'val': val}

functions = {
    'get_flux_value': (get_flux_value,),
}

ics = {
    'ic': ('Omega', {'T.0': T0}),
}

integral_1 = {
    'name': 'i',
    'order': 2,
}

equations = {
    'Temperature': """
    dw_dot.i.Omega( s, dT/dt ) + dw_laplace.i.Omega( m.D, s, T ) =
    + dw_integrate.i.Gamma_Right(flux.val, s)
    + dw_bc_newton.i.Gamma_Right(heat_loss.h_top, heat_loss.T_top_inf, s, T)
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

# Function to generate probe for temperature along the cylinder axis
def gen_probe():
    p0, p1 = nm.array([0., 0., 0.]), nm.array([0.1, 0., 0.])
    line_probe = LineProbe(p0, p1, n_point=100, share_geometry=True)
    return line_probe

line_probe = gen_probe()
probe_results = []
heat_flux_data = []

# Function to save heat flux at Gamma_Right boundary
def save_heat_flux(ts, problem, state):
    flux = problem.evaluate('ev_integrate.i.Gamma_Right(flux.val)')
    heat_flux_data.append((ts.time, flux))

def step_hook(pb, ts, variables):
    T_field = pb.get_variables()['T']
    pars, vals = line_probe(T_field)
    probe_results.append(vals)
    save_heat_flux(ts, pb, variables)

def post_process_hook(out, pb, state, extend=False):
    ts = pb.ts
    if ts.step == ts.n_step - 1:
        fig, (ax1, ax2) = plt.subplots(nrows=2)
        temperature_image = nm.array(probe_results).squeeze()
        m = ax1.imshow(temperature_image.T, origin='lower', aspect='auto')
        ax1.set_xlabel("time step")
        ax1.set_ylabel("distance along cylinder height (x-axis)")
        fig.colorbar(m, ax=ax1, label="temperature")
        ax2.plot(temperature_image.T[0], label="Gamma_Left")
        ax2.plot(temperature_image.T[-1], label="Gamma_Right")
        ax2.set_xlabel("time step")
        ax2.set_ylabel("temperature (°C)")
        ax2.legend()
        
        # Delta_T plot
        delta_T = temperature_image.T[-1] - temperature_image.T[0]
        fig, ax3 = plt.subplots()
        ax3.plot(delta_T, label="Delta_T")
        ax3.set_xlabel("distance along cylinder height (x-axis)")
        ax3.set_ylabel("Delta_T (°C)")
        ax3.legend()
        fig.savefig(os.path.join(pb.output_dir, 'delta_T_plot.png'), bbox_inches='tight')
        
        # Heat flux plot
        heat_flux_data_np = nm.array(heat_flux_data)
        fig, ax4 = plt.subplots()
        ax4.plot(heat_flux_data_np[:, 0], heat_flux_data_np[:, 1], label="Heat Flux")
        ax4.set_xlabel("time (s)")
        ax4.set_ylabel("Heat Flux (W/m^2)")
        ax4.legend()
        fig.savefig(os.path.join(pb.output_dir, 'heat_flux_plot.png'), bbox_inches='tight')
        
        fig.tight_layout()
        fig.savefig(os.path.join(pb.output_dir, 'heat_probe_time_evolution.png'),
                    bbox_inches='tight')
    return out

options = {
    'step_hook': 'step_hook',
    'post_process_hook': 'post_process_hook',
    'nls': 'newton',
    'ls': 'ls',
    'ts': 'ts',
    'save_times': 'all',
    'output_dir': './output_robin',  # Directory to save the output
}