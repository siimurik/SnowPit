import numpy as np
import csv
from sfepy.discrete.fem import Mesh, FEDomain, Field
from sfepy.discrete import Problem
from sfepy.solvers.ts import TimeStepper

def read_temp_and_hcoeff_from_csv(filename="t_o_and_h_o.csv"):
    """Read temperature and heat transfer coefficient data from CSV file."""
    t_o, h_o = [], []
    with open(filename, 'r') as csvfile:
        csvreader = csv.reader(csvfile)
        for row in csvreader:
            if len(row) >= 2:
                t_o.append(float(row[0].strip()))
                h_o.append(float(row[1].strip()))
    return t_o, h_o

def transient1D_sfepy(t_o, h_o, d_ins, lam_i, D, dx=0.005, dt=10.0, h_i=99.75):
    """
    Working version for SfePy 2025.1
    Returns temperature distribution with shape (n_nodes, n_hours).
    """
    # Parameters and initialization
    t_i = 0.0  # Inner temperature [°C]
    nr_hour = len(t_o)
    total_time = 3600 * nr_hour
    nh = int(3600 / dt)
    
    # Create proper 1D mesh - updated for SfePy 2025.1
    n_el = max(2, int(round(d_ins / dx)))
    vertices = np.linspace(0, d_ins, n_el + 1)[:, None]  # Nodes as 2D array
    cells = np.column_stack((np.arange(n_el), np.arange(1, n_el+1)))
    
    # Correct mesh creation syntax for current SfePy version
    mesh = Mesh.from_data('1d_mesh', 
                         vertices, 
                         None, 
                         [cells], 
                         [('1_2',)],  # Note the tuple format
                         '1_1')       # Geometric dimension
    
    # Create domain and regions
    domain = FEDomain('domain', mesh)
    omega = domain.create_region('Omega', 'all')
    gamma1 = domain.create_region('Gamma1', 'vertices in (x < 1e-10)', 'facet')
    gamma2 = domain.create_region('Gamma2', f'vertices in (x > {d_ins - 1e-10})', 'facet')
    
    # Create field and problem
    field = Field.from_args('temperature', np.float64, 'scalar', omega, approx_order=1)
    pb = Problem('heat_conduction')
    pb.domain = domain
    
    # Variables and initial condition
    t_var = pb.create_variables(['t'])
    pb.set_ics(ics={'t': {'Omega': 0.0}})
    
    # Material properties
    rho_c = lam_i / D  # ρc = λ/D [J/m³K]
    pb.materials = {
        'mat': ({'lam': lam_i, 'rho_c': rho_c},),
        'bc': ({'h_o': h_o[0], 't_o': t_o[0], 'h_i': h_i, 't_i': t_i},),
    }
    
    # Time-dependent BC update
    def update_bc(pb, hour):
        if hour < nr_hour:
            pb.materials['bc'].set_data({'h_o': h_o[hour], 't_o': t_o[hour]})
    
    # Define equation
    equation = """
    dw_laplace(mat.lam, t, v) + 
    dw_surface_robin(bc.h_o, t, v, Gamma1) + 
    dw_surface_robin(bc.h_i, t, v, Gamma2) + 
    dw_volume_dot(mat.rho_c/ts.dt, t, v) = 
    dw_surface_robin(bc.h_o, bc.t_o, v, Gamma1) + 
    dw_surface_robin(bc.h_i, bc.t_i, v, Gamma2) + 
    dw_volume_dot(mat.rho_c/ts.dt, t[-1], v)
    """
    pb.set_equations({'temperature': equation})
    
    # Results storage
    n_nodes = vertices.shape[0]
    T_nh = np.zeros((n_nodes, nr_hour))
    
    # Time stepping
    tss = TimeStepper(dt=dt, t0=0.0, t1=total_time)
    current_hour = -1
    
    for ts in tss:
        new_hour = int(ts.time // 3600)
        if new_hour != current_hour:
            current_hour = new_hour
            update_bc(pb, current_hour)
        
        state = pb.solve()
        
        if ts.step % nh == 0 and current_hour < nr_hour:
            T_nh[:, current_hour] = state['t']
    
    return T_nh

def main():
    try:
        # Read input data
        t_o, h_o = read_temp_and_hcoeff_from_csv()
        print(f"Read {len(t_o)} data points")
        print(f"First 5 temperatures: {t_o[:5]}")
        print(f"First 5 h coefficients: {h_o[:5]}")
        
        # Run simulation with adjusted parameters
        results = transient1D_sfepy(t_o, h_o, 
                                   d_ins=0.1, 
                                   lam_i=0.32, 
                                   D=1.5e-7,
                                   dx=0.01,  # Slightly larger dx for stability
                                   dt=60.0)  # Larger time step
        
        print("Simulation completed successfully!")
        print(f"Results shape: {results.shape}")
        
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        raise  # Re-raise to see full traceback

if __name__ == "__main__":
    main()