from fipy import CellVariable, Grid3D, Viewer, TransientTerm, DiffusionTerm
from fipy.tools import numerix
import numpy as np

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

# Solve the equation over time
timeStepDuration = 0.1
steps = 100
for step in range(steps):
    eq.solve(var=temperature, dt=timeStepDuration)

# Display results (this would typically require a 3D visualization tool)
Viewer(vars=temperature).plot()
