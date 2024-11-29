from fipy import CellVariable, Grid1D, DiffusionTerm

# 1D grid
mesh = Grid1D(dx=1.0, nx=10)

# Variable
temperature = CellVariable(name="Temperature", mesh=mesh, value=0.0)

# Heat equation
eq = DiffusionTerm(coeff=1.0) == 0

# Solve
temperature.constrain(100, where=mesh.facesLeft)  # Left boundary condition
temperature.constrain(0, where=mesh.facesRight)   # Right boundary condition
eq.solve(var=temperature)

print(temperature)
