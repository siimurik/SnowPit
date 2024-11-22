"""
Solves the heat equation PDE:

                2         
∂              ∂          
──(u(t, x)) = ───(u(t, x))  +  f
∂t              2         
              ∂x          

in a rod with x ∈ [0.0, 1.0]

with homogeneous Dirichlet Boundary Conditions (u(0.0) = u(1.0) = 0.0)

and a heat source f

using a sine initial condition:

      1 |                      ...........                      
        |                 .....           .....                 
        |              ...                     ...              
        |            ..                           ..            
0.55555 |----------..-------------------------------..----------
        |       ...                                   ...       
        |     ..                                         ..     
        |   ..                                             ..   
        | ..                                                 .. 
      0 |_______________________________________________________
         0                          0.5                        1
"""

import numpy as np
import matplotlib.pyplot as plt
from mpi4py import MPI
from dolfinx import mesh, fem
from ufl import TrialFunction, TestFunction, lhs, rhs, grad, dot, dx
from ufl.finiteelement import FiniteElement

# Create a mesh
n_elements = 32
domain = mesh.create_interval(MPI.COMM_WORLD, n_elements, [0.0, 1.0])

# Define a Function Space
# Create a UFL finite element and FunctionSpace
element = FiniteElement("CG", domain.ufl_cell(), 1)
lagrange_polynomial_space_first_order = fem.FunctionSpace(domain, element)

# The value of the solution on the boundary
u_on_boundary = fem.Constant(domain, 0.0)

# Define boundary conditions
boundary_dofs = fem.locate_dofs_geometrical(
    lagrange_polynomial_space_first_order,
    lambda x: np.isclose(x[0], 0.0) | np.isclose(x[0], 1.0),
)
boundary_condition = fem.dirichletbc(u_on_boundary, boundary_dofs)

# The initial condition, u(t=0, x) = sin(pi * x)
u_old = fem.Function(lagrange_polynomial_space_first_order)
u_old.interpolate(lambda x: np.sin(np.pi * x[0]))

# Plot the initial condition
plt.figure()
x = lagrange_polynomial_space_first_order.tabulate_dof_coordinates()
plt.plot(x[:, 0], u_old.x.array, label="t=0.0")

# Time-stepping parameters
time_step_length = 0.1
n_time_steps = 5

# The forcing term (heat source)
heat_source = fem.Constant(domain, 0.0)

# Define the variational problem
u_trial = TrialFunction(lagrange_polynomial_space_first_order)
v_test = TestFunction(lagrange_polynomial_space_first_order)

a = (
    u_trial * v_test * dx
    + time_step_length * dot(grad(u_trial), grad(v_test)) * dx
)
L = (
    u_old * v_test * dx
    + time_step_length * heat_source * v_test * dx
)

# Prepare for solving
u_solution = fem.Function(lagrange_polynomial_space_first_order)

# Time-stepping
for i in range(n_time_steps):
    # Assemble the linear system
    problem = fem.petsc.LinearProblem(a, L, bcs=[boundary_condition])
    u_solution = problem.solve()

    # Update the solution for the next step
    u_old.x.array[:] = u_solution.x.array

    # Plot the current solution
    plt.plot(x[:, 0], u_solution.x.array, label=f"t={(i+1)*time_step_length:1.1f}")

# Final plot
plt.legend()
plt.title("Heat Conduction in a Rod with Homogeneous Dirichlet BC")
plt.xlabel("x position")
plt.ylabel("Temperature")
plt.grid()
plt.show()
