import os
import h5py
from fipy import CellVariable, Grid3D, Viewer, TransientTerm, DiffusionTerm
from fipy.tools import numerix

# Define the mesh
nx, ny, nz = 30, 30, 30
dx, dy, dz = 0.1, 0.1, 0.1
mesh = Grid3D(nx=nx, ny=ny, nz=nz, dx=dx, dy=dy, dz=dz)

# Define the variable
temperature = CellVariable(name="temperature", mesh=mesh, value=273.15)  # starting at 0°C

# Define the coefficients
k_soil = 1.0  # Thermal conductivity of soil (W/m·K)
k_insulation = 0.1  # Thermal conductivity of insulation (W/m·K)
rho = 1000  # Density (kg/m³)
c = 1000  # Specific heat capacity (J/kg·K)

alpha_soil = k_soil / (rho * c)

# Define the PDE
eq = TransientTerm() == DiffusionTerm(coeff=alpha_soil)

# Set initial conditions
temperature.setValue(273.15)  # All at 0°C initially
temperature.setValue(253.15, where=(mesh.x < dx * 10) & (mesh.y < dy * 10) & (mesh.z < dz * 10))  # Cube interior at -20°C

# Set boundary conditions (example)
temperature.constrain(293.15, mesh.facesTop)  # Top face at 20°C due to sun warming

# Define time-stepping parameters
timeStepDuration = 0.1
steps = 100  # Number of timesteps

# Create output folder
output_folder = "output"
os.makedirs(output_folder, exist_ok=True)

# HDF5 file setup
hdf5_path = os.path.join(output_folder, "temperature_data.h5")
with h5py.File(hdf5_path, "w") as h5file:
    # Create a dataset to store all timesteps
    dset = h5file.create_dataset(
        "temperature", 
        shape=(steps, nx, ny, nz), 
        dtype='float64', 
        compression="gzip"  # Optional: compress the data
    )

    # Solve the equation over time
    for step in range(steps):
        eq.solve(var=temperature, dt=timeStepDuration)
        
        # Save the current timestep's data to the HDF5 dataset
        dset[step, :, :, :] = temperature.value.reshape((nx, ny, nz))
        print(f"Saved timestep {step} to HDF5 file")

print(f"All timesteps saved to {hdf5_path}")

# Display results (this would typically require a 3D visualization tool)
Viewer(vars=temperature).plot()
