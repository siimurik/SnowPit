import csv
import numpy as np
from sfepy.discrete.fem import Mesh, FEDomain, Field
from sfepy.discrete import Problem, Variables, Integral, Equation, Equations
from sfepy.terms import Term
from sfepy.solvers.ts import TimeStepper
from sfepy.solvers.nls import Newton
from sfepy.solvers.ls import ScipyDirect


def read_temp_and_hcoeff_from_csv(filename="t_o_and_h_o.csv"):
    """
    Read temperature and heat transfer coefficient data from CSV file.
    
    Parameters:
        filename (str): Path to CSV file (default: "t_o_and_h_o.csv")
        
    Returns:
        tuple: (t_o, h_o) where both are lists of floats
    """
    t_o = []
    h_o = []
    
    with open(filename, 'r') as csvfile:
        csvreader = csv.reader(csvfile)
        for row in csvreader:
            if len(row) >= 2:  # Ensure there are at least 2 columns
                try:
                    t_o.append(float(row[0].strip()))  # First column is temperature
                    h_o.append(float(row[1].strip()))  # Second column is h coefficient
                except ValueError:
                    print(f"Warning: Could not convert row {row} to numbers. Skipping.")
    
    return t_o, h_o

def transient1D_sfepy(t_o, h_o, d_ins, lam_i, D, dx=0.005, dt=10.0, h_i=99.75):
    """
    Temperature distribution, transient 1D with Robin BCs using SfePy's TimeStepper.
    
    Parameters:
        t_o (list): Outer temperature for each hour [°C]
        h_o (list): Outer heat transfer coefficient for each hour [W/m²K]
        d_ins (float): Insulation thickness [m]
        lam_i (float): Thermal conductivity [W/mK]
        D (float): Thermal diffusivity [m²/s]
        dx (float): Approximate element size [m]
        dt (float): Time step [s]
        h_i (float): Inner heat transfer coefficient [W/m²K]
        
    Returns:
        ndarray: Temperature distribution with shape (n_nodes, n_hours)
    """
    t_i = 0.0  # Inner temperature [°C]
    nr_hour = len(t_o)  # Number of hours
    total_time = 3600 * nr_hour  # Total simulation time [s]
    nh = int(3600 / dt)  # Number of time steps per hour
    
    # Create 1D mesh
    n_el = max(2, int(round(d_ins / dx)))  # Ensure at least 2 elements
    vertices = np.linspace(0, d_ins, n_el + 1)[:, None]  # Nodes as 2D array
    cells = np.array([[i, i+1] for i in range(n_el)])  # Line elements
    mesh = Mesh.from_data('1d_mesh', vertices, None, [cells], None, '1_2')
    
    # Create domain and regions
    domain = FEDomain('domain', mesh)
    omega = domain.create_region('Omega', 'all')
    gamma1 = domain.create_region('Gamma1', 'vertices in (x < 1e-10)', 'facet')  # Outer surface
    gamma2 = domain.create_region('Gamma2', 'vertices in (x > %f - 1e-10)' % d_ins, 'facet')  # Inner surface
    
    # Create field for temperature
    field = Field.from_args('temperature', np.float64, 'scalar', omega, approx_order=1)
    
    # Create variables
    t_var = Variables.from_conf([{
        'name': 't',
        'kind': 'unknown field',
        'field': 'temperature',
        'order': 0,
    }, {
        'name': 'v',
        'kind': 'test field',
        'field': 'temperature',
        'dual': 't',
    }], fields={'temperature': field})
    
    # Material properties (thermal conductivity and volumetric heat capacity)
    rho_c = lam_i / D  # ρc = λ/D [J/m³K]
    materials = {
        'mat': ({
            'lam': lam_i, 
            'rho_c': rho_c
        },),
        'bc': ({
            'h_o': h_o[0], 
            't_o': t_o[0], 
            'h_i': h_i, 
            't_i': t_i
        },),  # Initial BC values
    }
    
    # Create integral
    integral = Integral('i', order=2)
    
    # Define terms
    t1 = Term.new('dw_laplace(mat.lam, v, t)', integral, omega, mat=materials['mat'])
    t2 = Term.new('dw_surface_dot(bc.h_o, v, t, Gamma1)', integral, gamma1, bc=materials['bc'])
    t3 = Term.new('dw_surface_dot(bc.h_i, v, t, Gamma2)', integral, gamma2, bc=materials['bc'])
    t4 = Term.new('dw_volume_dot(mat.rho_c, v, t)', integral, omega, mat=materials['mat'])
    t5 = Term.new('dw_surface_integrate(bc.h_o, bc.t_o, v, Gamma1)', integral, gamma1, bc=materials['bc'])
    t6 = Term.new('dw_surface_integrate(bc.h_i, bc.t_i, v, Gamma2)', integral, gamma2, bc=materials['bc'])
    t7 = Term.new('dw_volume_dot(mat.rho_c, v, t[-1])', integral, omega, mat=materials['mat'])
    
    # Create equation
    eq = Equation('temperature', 
                  t1 + (dt**-1) * t4 + t2 + t3 - t5 - t6 - (dt**-1) * t7)
    eqs = Equations([eq])
    
    # Setup problem
    pb = Problem('heat_conduction', equations=eqs)
    pb.set_variables(t_var)
    pb.time_update(materials)
    
    # Set initial condition
    pb.set_ics({'t': {'Omega': 0.0}})
    
    # Configure solvers
    ls = ScipyDirect({})
    nls = Newton({}, lin_solver=ls)
    pb.set_solver(nls)
    
    # Initialize results storage
    n_nodes = vertices.shape[0]
    T_nh = np.zeros((n_nodes, nr_hour))
    
    # Time stepping
    tss = TimeStepper(dt=dt, t0=0.0, t1=total_time)
    current_hour = -1
    
    # Store initial state
    state_old = pb.create_state()
    state_old.set_default(0.0)
    
    for ts in tss:
        # Check if we moved to the next hour
        new_hour = int(ts.time // 3600)
        if new_hour != current_hour and new_hour < nr_hour:
            current_hour = new_hour
            # Update boundary condition parameters for current hour
            materials['bc'] = ({
                'h_o': h_o[current_hour], 
                't_o': t_o[current_hour], 
                'h_i': h_i, 
                't_i': t_i
            },)
            pb.time_update(materials)
        
        # Update time stepping data
        pb.update_time_stepper(ts)
        
        # Solve the time step
        state = pb.solve(state0=state_old)
        
        # Store results at each full hour
        if ts.step % nh == 0 and current_hour < nr_hour and current_hour >= 0:
            T_nh[:, current_hour] = state()['t']
        
        # Update old state for next time step
        state_old = state.copy()
    
    return T_nh

def main():
    # Read the data from CSV
    try:
        t_o, h_o = read_temp_and_hcoeff_from_csv()
        print(f"Read {len(t_o)} data points")
        print(f"First 5 temperatures: {t_o[:5]}")
        print(f"First 5 h coefficients: {h_o[:5]}")

        # Use with your SfePy function
        results = transient1D_sfepy(t_o, h_o, d_ins=0.1, lam_i=0.32, D=1.5e-7)
        print(f"Results shape: {results.shape}")
        
    except FileNotFoundError:
        print("CSV file not found. Please ensure 't_o_and_h_o.csv' exists.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    main()