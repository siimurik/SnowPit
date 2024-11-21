import numpy as np
import numba as nb
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Parameters 
L, W, D = 1, 1, 1                # Pit dimensions in meters
dx, dy, dz = 0.1, 0.1, 0.1      # Spatial step
dt = 0.1                         # Reduced time step (change this as needed)
k, rho, c_p = 0.3, 400.0, 2100.0  # Thermal conductivity, density, specific heat
alpha = k / (rho * c_p)         # Thermal diffusivity
nx, ny, nz = int(L/dx), int(W/dy), int(D/dz)  # Number of grid points

# CFL Condition Check
CFL_x = alpha * dt / dx**2
CFL_y = alpha * dt / dy**2
CFL_z = alpha * dt / dz**2

if CFL_x <= 0.5 and CFL_y <= 0.5 and CFL_z <= 0.5:
    print("CFL condition is satisfied.")
else:
    print("CFL condition is NOT satisfied.")

# Initialize temperature field
T_initial = np.full((nx, ny, nz), -10.0)  # Initial snow temperature (e.g., -10°C)
T_ground = 5.0  # Ground temperature (°C)

# Number of timesteps
n_steps = 1000  

@nb.njit
def time_step_with_bc(T, alpha, dx, dy, dz, dt, n_steps, T_ground, T_ext_values, k_wood, h_wood, d_wood):
    T_new = np.empty_like(T)
    nx, ny, nz = T.shape

    for step in range(n_steps):
        T_ext = T_ext_values[step]

        for i in range(1, nx-1):
            for j in range(1, ny-1):
                for k in range(1, nz-1):
                    T_new[i, j, k] = T[i, j, k] + alpha * dt * (
                        (T[i+1, j, k] - 2*T[i, j, k] + T[i-1, j, k]) / dx**2 +
                        (T[i, j+1, k] - 2*T[i, j, k] + T[i, j-1, k]) / dy**2 +
                        (T[i, j, k+1] - 2*T[i, j, k] + T[i, j, k-1]) / dz**2
                    )

        # Apply boundary conditions
        T_new[0, :, :] = T_ground  # Left
        T_new[-1, :, :] = T_ground  # Right
        T_new[:, 0, :] = T_ground  # Front
        T_new[:, -1, :] = T_ground  # Back
        T_new[:, :, 0] = T_ground  # Bottom

        # 2. Robin BC for the top
        for i in range(nx):
            for j in range(ny):
                dT_dz_top = (T[i, j, -1] - T[i, j, -2]) / dz
                T_new[i, j, -1] = T[i, j, -1] + dt * (
                    -k_wood * dT_dz_top / d_wood + h_wood * (T_ext - T[i, j, -1])
                )

        # Smoothing step (can be adjusted or removed)
        #for i in range(nx):
        #    for j in range(ny):
        #        # Smooth bottom
        #        T_new[i, j, 0] = 0.5 * (T_new[i, j, 0] + T[i, j, 1])  
        #        # Smooth top
        #        T_new[i, j, -1] = 0.5 * (T_new[i, j, -1] + T[i, j, -2])  
        #        
        #        # Smooth left edge
        #        if i == 0:  # Left edge
        #            T_new[i, j, :] = 0.5 * (T_new[i, j, :] + T_new[i + 1, j, :])
        #        
        #        # Smooth right edge
        #        if i == nx - 1:  # Right edge
        #            T_new[i, j, :] = 0.5 * (T_new[i, j, :] + T_new[i - 1, j, :])
        #        
        #        # Smooth front edge
        #        if j == 0:  # Front edge
        #            T_new[i, j, :] = 0.5 * (T_new[i, j, :] + T_new[i, j + 1, :])
        #        
        #        # Smooth back edge
        #        if j == ny - 1:  # Back edge
        #            T_new[i, j, :] = 0.5 * (T_new[i, j, :] + T_new[i, j - 1, :])

        # Swap arrays for the next iteration
        T, T_new = T_new, T

        # Debugging output: print max and min temperature
        if step % 100 == 0:  # Print every 100 steps
            print("Step: ", step, "Max temp", np.round(np.max(T), 4), "Min temp", np.round(np.min(T), 4))

    return T

# Simulation parameters
k_wood = 0.15   # Thermal conductivity of woodchips (W/m·K)
h_wood = 10.0   # Heat transfer coefficient (W/m²·K)
d_wood = 0.2    # Thickness of woodchip layer (m)

# External temperature as a sinusoidal function (e.g., diurnal variation)
time = np.arange(0, n_steps * dt, dt)  # Time array
T_ext_values = 15.0 + 10.0 * np.sin(2 * np.pi * time / (24 * 3600))  # Increased amplitude

# Initialize temperature field (snow at -10°C)
T_initial = np.full((nx, ny, nz), -10.0)

# Run the simulation
T_final = time_step_with_bc(
    T_initial, alpha, dx, dy, dz, dt, n_steps,
    T_ground, T_ext_values, k_wood, h_wood, d_wood
)

# Debugging output: Final maximum and minimum temperatures
print(f"Final Max Temp: {np.max(T_final)}, Final Min Temp: {np.min(T_final)}")

# Create a meshgrid for the 3D plot
X, Y = np.meshgrid(np.linspace(0, L, nx), np.linspace(0, W, ny))
Z_mid = T_final[:, :, nz // 2]  # Take a slice at the midpoint of the depth

# Create a 3D surface plot
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(X, Y, Z_mid, cmap='viridis', edgecolor='none')

# Customize the plot
ax.set_title('Temperature Distribution at Mid Depth')
ax.set_xlabel('Width (m)')
ax.set_ylabel('Length (m)')
ax.set_zlabel('Temperature (°C)')
fig.colorbar(surf, ax=ax, label='Temperature (°C)')

plt.show()
