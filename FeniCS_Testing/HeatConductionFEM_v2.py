from fenics import *
import numpy as np

# Parameters
L = 1.0  # Length of the cube (1m x 1m x 1m)
T = 300  # Initial temperature of the cube (Kelvin)
heat_source_strength = 1.0  # Strength of the heat source (arbitrary units)

# Create mesh and define function space
mesh = BoxMesh(Point(0, 0, 0), Point(L, L, L), 10, 10, 10)  # 10 divisions along each axis
V = FunctionSpace(mesh, 'P', 1)  # Linear elements

# Dirichlet boundary condition on 5 faces
class BottomAndSidesBoundary(SubDomain):
    def inside(self, x, on_boundary):
        # Apply Dirichlet condition on x[2] == 0 (bottom), x[0] == 0 or x[0] == L (sides)
        return on_boundary and (x[2] == 0 or x[0] == 0 or x[1] == 0 or x[0] == L or x[1] == L)

# Robin boundary condition on the top face
#class TopBoundary(SubDomain):
#    def inside(self, x, on_boundary):
#        return near(x[2], L)  # Top face
#
#top_bc = RobinBC(V, k=1.0, T_inf=300, boundary=TopBoundary())


# Create the boundary condition (for 5 faces)
bc1 = DirichletBC(V, T, BottomAndSidesBoundary())

# Define function and variational problem
u = TrialFunction(V)
v = TestFunction(V)

# Initial temperature field (use a spatially-varying function as initial condition)
u_0 = Expression('300 + 100*exp(-10*(pow(x[0]-0.5, 2) + pow(x[1]-0.5, 2) + pow(x[2]-0.5, 2)))', degree=2)
u_n = interpolate(u_0, V)

# Time-stepping parameters
dt = 0.1
t_end = 10.0
t = 0.0

# Define the weak form of the heat equation (with a heat source term)
a = u*v*dx + dt*dot(grad(u), grad(v))*dx
L = u_n*v*dx + dt*Constant(heat_source_strength)*v*dx  # Add heat source term

# Output file for visualization
file = File('heat_cubefile.pvd')  # Output file for visualization
file << u_n  # Output initial condition

# Time-stepping loop
while t < t_end:
    solve(a == L, u_n, [bc1])  # Solve for the next time step

    # Update the solution and output it
    file << u_n  # Output current solution

    u_n.assign(u_n)  # Update the solution for the next time step
    t += dt

