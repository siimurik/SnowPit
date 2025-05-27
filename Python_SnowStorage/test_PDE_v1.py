import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.integrate import solve_ivp
from scipy import sparse
from scipy.sparse.linalg import spsolve

# ----- SIMULATION PARAMETERS -----
# Spatial parameters
L = 0.5  # Snow layer thickness [m]
nx = 50  # Number of grid points
dx = L / (nx - 1)  # Grid spacing [m]
x = np.linspace(0, L, nx)  # Spatial grid [m]

# Temporal parameters
hours = 24  # Simulation duration [h]
t_end = hours * 3600  # Simulation end time [s]
dt_output = 1800  # Output time step [s]
t_eval = np.arange(0, t_end + dt_output, dt_output)  # Times to evaluate solution

# ----- MATERIAL PROPERTIES -----
# General properties
rho_i = 917    # Density of ice [kg/m³]
rho_l = 1000   # Density of water [kg/m³]
rho_v = 0.7    # Reference vapor density [kg/m³]
g = 9.81       # Gravitational acceleration [m/s²]

# Initial snow properties (vary with depth)
def get_initial_snow_properties():
    # Snow density varies with depth (increases with depth)
    porosity = 0.7 - 0.3 * (x/L)  # Porosity decreases with depth
    rho_s = (1 - porosity) * rho_i  # Snow density [kg/m³]
    
    # Derived properties based on density
    # Thermal conductivity [W/(m·K)] - Using Calonne et al. (2011) empirical relation
    k_t = 0.024 - 1.23e-4 * rho_s + 2.5e-6 * rho_s**2
    
    # Specific heat capacity [J/(kg·K)]
    c_p = 2100  # Approximate value for ice/snow
    
    # Liquid water parameters
    k_s = 1e-10 * (1 + 5 * (x/L))  # Hydraulic conductivity [m²], increases with depth
    D_l = 2e-9 * np.ones_like(x)   # Liquid diffusivity [m²/s]
    mu = 1.8e-3  # Water viscosity [Pa·s]
    
    # Vapor transport parameters
    D_v = 2.2e-5 * porosity  # Vapor diffusivity [m²/s], varies with porosity
    
    # Latent heat
    L_v = 2.5e6  # Latent heat of vaporization [J/kg]
    
    return {
        'rho_s': rho_s,
        'porosity': porosity,
        'k_t': k_t,
        'c_p': c_p,
        'k_s': k_s,
        'D_l': D_l,
        'mu': mu,
        'D_v': D_v,
        'L_v': L_v,
        'rho_m': rho_s  # Effective moisture density (approx. as snow density)
    }

# ----- EXTERNAL CONDITIONS (HOURLY) -----
def get_external_conditions(t):
    """Return external conditions based on time in seconds"""
    hour = (t / 3600) % 24  # Convert to hour of day
    
    # Temperature with diurnal cycle [°C]
    T_ext = -5 + 10 * np.sin(2 * np.pi * (hour - 10) / 24)
    
    # Solar radiation [W/m²]
    if 6 <= hour <= 18:  # Daytime
        G_solar = 500 * np.sin(np.pi * (hour - 6) / 12)
    else:  # Nighttime
        G_solar = 0
    
    # Atmospheric humidity [kg/m³]
    rho_v_ext = 0.003 + 0.001 * np.sin(2 * np.pi * (hour - 12) / 24)
    
    # Rain flux [kg/(m²·s)]
    rain_flux = 0
    if 3 <= hour <= 5:  # Light snow/rain in early morning
        rain_flux = 1e-5
    
    # Atmospheric pressure [Pa]
    p_atm = 101325
    
    return {
        'T_ext': T_ext + 273.15,  # Convert to Kelvin
        'G_solar': G_solar,
        'rho_v_ext': rho_v_ext,
        'rain_flux': rain_flux,
        'p_atm': p_atm
    }

# ----- BOUNDARY CONDITIONS -----
# Heat transfer coefficients
h_ext = 15    # External convective heat transfer coefficient [W/(m²·K)]
h_int = 5     # Internal heat transfer coefficient [W/(m²·K)]
alpha = 0.7   # Solar absorption coefficient

# Moisture transfer coefficient
h_m = 0.01    # Mass transfer coefficient [m/s]

# Internal temperature [K]
T_int = 273.15 - 1  # Slightly below freezing

# ----- PHASE CHANGE MODEL -----
def calculate_evaporation_rate(T, C_l, C_v, props):
    """Calculate evaporation rate based on current state"""
    # Saturation vapor density at current temperature [kg/m³]
    T_C = T - 273.15  # Convert to Celsius
    rho_v_sat = 0.611 * np.exp(17.27 * T_C / (T_C + 237.3)) / (461.5 * T) * 1000
    
    # Evaporation rate proportional to departure from equilibrium
    k_evap = 1e-4  # Evaporation rate coefficient [1/s]
    
    # Evaporation is limited by available liquid water
    m_evap = k_evap * (rho_v_sat - C_v) * (C_l > 0)
    
    # Ensure we don't evaporate more water than available
    m_evap = np.minimum(m_evap, C_l / dt_output)
    return np.maximum(m_evap, 0)  # No condensation in this simple model

# ----- INITIAL CONDITIONS -----
def initial_conditions():
    # Get material properties
    props = get_initial_snow_properties()
    
    # Initial temperature profile [K]
    T_0 = 273.15 - 10 + 8 * (1 - x/L)  # Colder at surface, warmer at bottom
    
    # Initial liquid water content [kg/m³]
    C_l_0 = 5 * np.exp(-10 * x/L)  # More moisture near surface
    
    # Initial vapor content [kg/m³]
    # Based on temperature-dependent saturation vapor density
    T_C = T_0 - 273.15
    C_v_0 = 0.8 * (0.611 * np.exp(17.27 * T_C / (T_C + 237.3)) / (461.5 * T_0) * 1000)
    
    # Combine into single state vector
    y0 = np.concatenate([T_0, C_l_0, C_v_0])
    
    return y0, props

# ----- FINITE DIFFERENCE MATRICES -----
def create_fd_matrices():
    # Create standard second derivative operator (center points only)
    main_diag = -2 * np.ones(nx)
    off_diag = np.ones(nx-1)
    diagonals = [off_diag, main_diag, off_diag]
    offsets = [-1, 0, 1]
    D2 = sparse.diags(diagonals, offsets, shape=(nx, nx), format='csr')
    D2 = D2 / dx**2
    
    # Create standard first derivative operator (center points only)
    main_diag_D1 = np.zeros(nx)
    off_diag_upper = 0.5 * np.ones(nx-1)
    off_diag_lower = -0.5 * np.ones(nx-1)
    diagonals_D1 = [off_diag_lower, main_diag_D1, off_diag_upper]
    D1 = sparse.diags(diagonals_D1, offsets, shape=(nx, nx), format='csr')
    D1 = D1 / dx
    
    # Convert to CSR format for easier manipulation
    D2_array = D2.toarray()
    D1_array = D1.toarray()
    
    # Manually modify boundary rows for boundary conditions
    # Second derivative operator
    D2_array[0, 0:3] = np.array([2, -2, 0]) / dx**2  # Forward difference at boundary
    D2_array[-1, -3:] = np.array([0, -2, 2]) / dx**2  # Backward difference at boundary
    
    # First derivative operator
    D1_array[0, 0:3] = np.array([-3, 4, -1]) / (2*dx)  # One-sided 2nd order at boundary
    D1_array[-1, -3:] = np.array([1, -4, 3]) / (2*dx)  # One-sided 2nd order at boundary
    
    # Convert back to sparse format
    D2 = sparse.csr_matrix(D2_array)
    D1 = sparse.csr_matrix(D1_array)
    
    return D1, D2

# ----- ODE SYSTEM -----
def snow_pde_system(t, y, D1_array, D2_array):
    """
    Convert the PDE system to ODEs for scipy.integrate.solve_ivp
    t: time
    y: state vector [T, C_l, C_v]
    D1_array, D2_array: finite difference matrices (as dense arrays)
    """
    # Get current material properties
    props = get_initial_snow_properties()  # In a real model, these might be state-dependent
    
    # Get external conditions
    ext = get_external_conditions(t)
    
    # Extract state variables
    T = y[:nx]
    C_l = y[nx:2*nx]
    C_v = y[2*nx:3*nx]
    
    # Ensure non-negative values
    C_l = np.maximum(C_l, 0)
    C_v = np.maximum(C_v, 0)
    
    # Calculate evaporation rate
    m_evap = calculate_evaporation_rate(T, C_l, C_v, props)
    
    # ---- Heat Transport Equation ----
    # Interior points
    k_t_avg = 0.5 * (props['k_t'][:-1] + props['k_t'][1:])  # Average thermal conductivity
    k_t_avg = np.append(k_t_avg, props['k_t'][-1])
    k_t_avg = np.insert(k_t_avg, 0, props['k_t'][0])
    
    # Apply spatial operators
    dT_dt = np.zeros_like(T)
    
    # Internal heat diffusion
    for i in range(1, nx-1):
        dT_dt[i] = (1 / (props['rho_m'][i] * props['c_p'])) * (
            (k_t_avg[i+1] * (T[i+1] - T[i]) - k_t_avg[i] * (T[i] - T[i-1])) / dx**2
        )
        
        # Add latent heat term
        dT_dt[i] += props['L_v'] * m_evap[i] / (props['rho_m'][i] * props['c_p'])
    
    # Boundary conditions for heat
    # Surface (x=0)
    dT_dt[0] = (1 / (props['rho_m'][0] * props['c_p'])) * (
        h_ext * (ext['T_ext'] - T[0]) + alpha * ext['G_solar']
    ) / dx
    
    # Interior boundary (x=L)
    dT_dt[-1] = (1 / (props['rho_m'][-1] * props['c_p'])) * (
        h_int * (T_int - T[-1])
    ) / dx
    
    # ---- Liquid Moisture Transport Equation ----
    dC_l_dt = np.zeros_like(C_l)
    
    # Interior liquid transport
    for i in range(1, nx-1):
        # Diffusion term
        dC_l_dt[i] = (
            (props['D_l'][i+1] * (C_l[i+1] - C_l[i]) - props['D_l'][i-1] * (C_l[i] - C_l[i-1])) / dx**2
        )
        
        # Gravity term - simplified as dz/dx = 1 for vertical orientation
        dC_l_dt[i] += (
            (props['k_s'][i] / props['mu']) * rho_l * g
        ) / dx
        
        # Phase change term
        dC_l_dt[i] -= m_evap[i]
    
    # Boundary conditions for liquid moisture
    # Surface (x=0): rain input
    dC_l_dt[0] = ext['rain_flux'] / dx
    
    # Interior boundary (x=L): no flux
    dC_l_dt[-1] = 0
    
    # ---- Vapor Transport Equation ----
    dC_v_dt = np.zeros_like(C_v)
    
    # Interior vapor diffusion
    for i in range(1, nx-1):
        dC_v_dt[i] = (
            (props['D_v'][i+1] * (C_v[i+1] - C_v[i]) - props['D_v'][i-1] * (C_v[i] - C_v[i-1])) / dx**2
        )
        
        # Add evaporation source term
        dC_v_dt[i] += m_evap[i]
    
    # Boundary conditions for vapor
    # Surface (x=0): vapor exchange with atmosphere
    dC_v_dt[0] = h_m * (ext['rho_v_ext'] - C_v[0]) / dx
    
    # Interior boundary (x=L): no flux
    dC_v_dt[-1] = 0
    
    # Combine into single derivative vector
    dy_dt = np.concatenate([dT_dt, dC_l_dt, dC_v_dt])
    
    return dy_dt

# ----- SOLVE THE SYSTEM -----
def solve_snow_system():
    print("Setting up initial conditions...")
    # Get initial conditions
    y0, props = initial_conditions()
    
    # Create finite difference matrices
    D1, D2 = create_fd_matrices()
    
    # Convert to arrays for the solver
    D1_array = D1.toarray()
    D2_array = D2.toarray()
    
    print(f"Solving system for {hours} hours...")
    # Solve the system
    solution = solve_ivp(
        lambda t, y: snow_pde_system(t, y, D1_array, D2_array),
        [0, t_end],
        y0,
        method='BDF',  # Backward Differentiation Formula for stiff systems
        t_eval=t_eval,
        rtol=1e-4,
        atol=1e-6
    )
    
    print(f"Solution complete! Shape: {solution.y.shape}")
    return solution, props

# ----- VISUALIZE RESULTS -----
def visualize_results(solution, props):
    # Extract solution components
    T = solution.y[:nx, :]
    C_l = solution.y[nx:2*nx, :]
    C_v = solution.y[2*nx:3*nx, :]
    
    # Convert time to hours for plotting
    t_hours = solution.t / 3600
    
    # Plot temperature profiles at different times
    plt.figure(figsize=(12, 8))
    
    # Temperature subplot
    plt.subplot(3, 1, 1)
    times_to_plot = [0, 6, 12, 18, 23]  # Hours to plot
    for i, hour in enumerate(times_to_plot):
        idx = np.argmin(np.abs(t_hours - hour))
        plt.plot(x, T[:, idx] - 273.15, label=f'{hour}h')
    
    plt.xlabel('Depth (m)')
    plt.ylabel('Temperature (°C)')
    plt.title('Temperature Profile in Snow')
    plt.legend()
    plt.grid(True)
    
    # Liquid moisture subplot
    plt.subplot(3, 1, 2)
    for i, hour in enumerate(times_to_plot):
        idx = np.argmin(np.abs(t_hours - hour))
        plt.plot(x, C_l[:, idx], label=f'{hour}h')
    
    plt.xlabel('Depth (m)')
    plt.ylabel('Liquid Water Content (kg/m³)')
    plt.title('Liquid Water Distribution')
    plt.legend()
    plt.grid(True)
    
    # Vapor subplot
    plt.subplot(3, 1, 3)
    for i, hour in enumerate(times_to_plot):
        idx = np.argmin(np.abs(t_hours - hour))
        plt.plot(x, C_v[:, idx] * 1000, label=f'{hour}h')  # Convert to g/m³ for better visibility
    
    plt.xlabel('Depth (m)')
    plt.ylabel('Vapor Content (g/m³)')
    plt.title('Water Vapor Distribution')
    plt.legend()
    plt.grid(True)
    
    plt.tight_layout()
    plt.savefig('snow_profiles.png', dpi=300)
    plt.show()
    
    # Plot time evolution at different depths
    plt.figure(figsize=(12, 8))
    
    # Depths to plot
    depths = [0, 0.1, 0.25, 0.45]  # in meters
    depth_indices = [np.argmin(np.abs(x - d)) for d in depths]
    
    # Temperature subplot
    plt.subplot(3, 1, 1)
    for i, idx in enumerate(depth_indices):
        plt.plot(t_hours, T[idx, :] - 273.15, label=f'{x[idx]:.2f}m')
    
    # Also plot external temperature
    ext_temp = np.array([get_external_conditions(t)['T_ext'] for t in solution.t]) - 273.15
    plt.plot(t_hours, ext_temp, 'k--', label='External')
    
    plt.xlabel('Time (hours)')
    plt.ylabel('Temperature (°C)')
    plt.title('Temperature Evolution')
    plt.legend()
    plt.grid(True)
    
    # Liquid moisture subplot
    plt.subplot(3, 1, 2)
    for i, idx in enumerate(depth_indices):
        plt.plot(t_hours, C_l[idx, :], label=f'{x[idx]:.2f}m')
    
    plt.xlabel('Time (hours)')
    plt.ylabel('Liquid Water Content (kg/m³)')
    plt.title('Liquid Water Evolution')
    plt.legend()
    plt.grid(True)
    
    # Vapor subplot
    plt.subplot(3, 1, 3)
    for i, idx in enumerate(depth_indices):
        plt.plot(t_hours, C_v[idx, :] * 1000, label=f'{x[idx]:.2f}m')
    
    # Also plot external vapor density
    ext_vapor = np.array([get_external_conditions(t)['rho_v_ext'] for t in solution.t]) * 1000
    plt.plot(t_hours, ext_vapor, 'k--', label='External')
    
    plt.xlabel('Time (hours)')
    plt.ylabel('Vapor Content (g/m³)')
    plt.title('Water Vapor Evolution')
    plt.legend()
    plt.grid(True)
    
    plt.tight_layout()
    plt.savefig('snow_evolution.png', dpi=300)
    plt.show()

# ----- MAIN EXECUTION -----
if __name__ == "__main__":
    solution, props = solve_snow_system()
    visualize_results(solution, props)