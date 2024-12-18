# Transient Advection-Diffusion Equation Solver

This project solves the transient advection-diffusion equation with a divergence-free advection velocity using the finite element method (FEM) with SfePy.

## Equation

The equation being solved is:

$$
\int_{\Omega} s \frac{\partial T}{\partial t} + \int_{\Omega} D \nabla s \cdot \nabla T = 0 \;, \quad \forall s \in H^1(\Omega).
$$

where:
- $T $ is the temperature.
- $D $ is the diffusion coefficient.
- $s $ is the test function.
- $\Omega$ is the domain over which the problem is defined (aka a volume in space, specifically a cube in our case).

## Boundary Conditions

- **Robin Boundary Condition**: Applied on the top face with ambient temperature, heat transfer coefficient, and solar radiation heat flux.
- **Dirichlet Boundary Conditions**: Applied on the sides and bottom face with specified temperatures.

## Specifications

- **Mesh**: Uses a 3D cube mesh (`cube_medium_hexa.mesh`).
- **Time Settings**: 
  - Start time: 0.0 s
  - End time: 10.0 s
  - Number of time steps: 100
- **Materials**: 
  - Diffusion coefficient $D $: 0.01 m²/s
- **Regions**: 
  - `Omega`: Entire domain
  - `Left`, `Right`, `Front`, `Back`, `Bottom`, `Top`: Facets of the cube
- **Fields**: 
  - `temperature`: Scalar field with linear approximation
- **Variables**: 
  - `T`: Unknown temperature field
  - `s`: Test field
- **Initial Conditions**: 
  - Temperature based on depth
- **Boundary Conditions**: 
  - Robin condition on the top face
  - Dirichlet conditions on the sides and bottom face
- **Solvers**: 
  - Linear solver: `ls.scipy_direct`
  - Nonlinear solver: `nls.newton`
  - Time-stepping solver: `ts.simple`

## Running the Code

To run the simulation, execute the script with SfePy. View the results using:

```sh
sfepy-view cube_medium_hexa.*.vtk -f T:wu 1:vw
```

## Dependencies

- SfePy
- NumPy
- JAX

I understand it can seem quite complex at first! Let's break it down into simpler terms.

### What is the Weak Form?

The **weak form** of an equation is a way of rewriting a differential equation so that it can be solved more easily using numerical methods like the finite element method (FEM). Here's a step-by-step explanation:

1. **Differential Equation**: This is the original form of the equation, which involves derivatives. For example, the heat equation involves the second derivative of temperature with respect to space.

2. **Integral Form**: To convert the differential equation into an integral form, we multiply it by a test function $s$and integrate over the domain $\Omega$. This process is called "testing" the equation.

3. **Integration by Parts**: This step reduces the order of the derivatives, making the equation easier to handle numerically. It also naturally incorporates boundary conditions.

### Why is it Necessary?

- **Numerical Stability**: The weak form often provides better numerical stability and convergence properties.
- **Flexibility**: It can handle complex geometries and boundary conditions more effectively.
- **Implementation**: It simplifies the implementation of boundary conditions and the assembly of the system of equations.

### What Makes it "Weak"?

The term "weak" doesn't imply that the formulation is lazy or inferior. Instead, it means that the solution is not required to satisfy the differential equation at every point (which is the "strong" form). Instead, it only needs to satisfy the equation in an "averaged" sense over the domain. This relaxation makes it easier to find approximate solutions.

### Example

Consider the heat equation in one dimension:

$$
-\frac{d^2 T}{dx^2} = f
$$

In the weak form, we multiply by a test function $s $and integrate:

$$
\int_{\Omega} s \left( -\frac{d^2 T}{dx^2} \right) \, dx = \int_{\Omega} s f \, dx
$$

Using integration by parts, we get:

$$
\int_{\Omega} \frac{d s}{dx} \frac{d T}{dx} \, dx = \int_{\Omega} s f \, dx
$$

This is the weak form of the heat equation. It involves lower-order derivatives, making it more suitable for numerical methods.

### Summary

- **Weak Form**: A reformulation of a differential equation that is easier to solve numerically.
- **Necessary**: Provides numerical stability, flexibility, and easier implementation.
- **"Weak"**: Indicates that the solution satisfies the equation in an averaged sense, not pointwise.

I hope this makes it clearer! Feel free to ask if you have more questions.

Great observation! The minus sign disappears due to the process of integration by parts. Let's go through it step-by-step:

### Original Equation

We start with the differential equation:

$$
-\frac{d^2 T}{dx^2} = f
$$

### Multiplying by a Test Function

We multiply both sides by a test function $s $and integrate over the domain $\Omega $:

$$
\int_{\Omega} s \left( -\frac{d^2 T}{dx^2} \right) \, dx = \int_{\Omega} s f \, dx
$$

### Integration by Parts

Integration by parts is used to move the derivative from $T $to $s $. The formula for integration by parts is:

$$
\int u \frac{dv}{dx} \, dx = uv \bigg|_{\partial \Omega} - \int v \frac{du}{dx} \, dx
$$

In our case, let $u = s $and $dv = \frac{d^2 T}{dx^2} dx $. Applying integration by parts:

$$
\int_{\Omega} s \left( -\frac{d^2 T}{dx^2} \right) \, dx = - \left[ s \frac{dT}{dx} \right]_{\partial \Omega} + \int_{\Omega} \frac{ds}{dx} \frac{dT}{dx} \, dx
$$

### Boundary Terms

The term $\left[ s \frac{dT}{dx} \right]_{\partial \Omega} $represents the boundary terms. If $s $is chosen such that it vanishes on the boundary (common in FEM), this term disappears.

### Resulting Weak Form

We are left with:

$$
\int_{\Omega} \frac{ds}{dx} \frac{dT}{dx} \, dx = \int_{\Omega} s f \, dx
$$

The minus sign disappears because it is absorbed during the integration by parts process. The derivative of $T $is transferred to $s $, and the sign change is a result of this transfer.

### Summary

- The minus sign disappears due to integration by parts.
- Integration by parts transfers the derivative from $T $to $s $, changing the sign in the process.