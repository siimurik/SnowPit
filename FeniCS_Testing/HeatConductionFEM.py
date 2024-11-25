from fenics import *
import numpy as np
import matplotlib.pyplot as plt

# Import mesh from XML file
mesh = Mesh('outputTetMesh.xml')

# Define function space
V = FunctionSpace(mesh, 'P', 1)

# Define Dirichlet boundary condition for 5 faces
class DirichletBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and (
            near(x[0], -1, tol=1e-6) or near(x[1], -1, tol=1e-6) or near(x[2], -1, tol=1e-6) or
            near(x[0], 1, tol=1e-6) or near(x[1], 1, tol=1e-6)
        )

# Create boundary condition object
dirichlet_boundary = DirichletBoundary()
boundaries = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
boundaries.set_all(0)
dirichlet_boundary.mark(boundaries, 1)

# Define Dirichlet BC value and apply
u_D = Constant(0)
bc_dirichlet = DirichletBC(V, u_D, boundaries, 1)

# Define Robin boundary condition for the top insulated with woodchips
class RobinBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[2], 1, tol=1e-6) and on_boundary

# Mark Robin boundary
robin_boundary = RobinBoundary()
robin_boundary.mark(boundaries, 2)

# Print boundary markings for debugging
print("Boundary Markings:")
for i in range(boundaries.size()):
    print(f"Facet {i}: {boundaries[i]}")

alpha = Constant(1.0)  # Coefficient for Robin boundary condition
beta = Constant(0.0)  # Coefficient for Robin boundary condition
ds = Measure('ds', domain=mesh, subdomain_data=boundaries)

# Combine boundary conditions
bc_combined = [bc_dirichlet]

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(0)  # Source term
a = dot(grad(u), grad(v)) * dx + alpha * u * v * ds(2)  # Robin BC term
L = f * v * dx + beta * v * ds(2)  # Source term

# Compute solution
u = Function(V)
solve(a == L, u, bc_combined)

# Plot solution
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Extract mesh points for plotting
x = mesh.coordinates()[:, 0]
y = mesh.coordinates()[:, 1]
z = u.compute_vertex_values(mesh)  # Get the solution values at mesh points

# Create a scatter plot for the solution
sc = ax.scatter(x, y, z, c=z, cmap='viridis')  # Use a colormap for better visualization

# Set labels
ax.set_xlabel('X-axis')
ax.set_ylabel('Y-axis')
ax.set_zlabel('Solution u')

# Show color bar
plt.colorbar(sc, ax=ax, label='Solution value')

# Show plot
plt.show()

# Plot boundaries for debugging
plot(boundaries, title='Boundary Markers')
plt.show()
