# Documentation

Download FVTool for MATLAB from this [link](https://se.mathworks.com/matlabcentral/fileexchange/46637-a-simple-finite-volume-solver-for-matlab?s_tid=FX_rc3_behav). 

For the code to work, first run `FVToolStartUp.m` from the FVTool library. 

$$
\begin{cases}
\alpha \frac{\partial T}{\partial t} + \nabla \cdot \left(-D \nabla T\right) + \nabla \cdot (uT) = 0 & \text{in the domain}, \\
T = T_{\text{boundary}} & \text{on Dirichlet boundaries}, \\
-D \nabla T \cdot \mathbf{n} + h (T - T_{\infty}) = q & \text{on Robin boundaries}.
\end{cases}
$$

In this system:
- $\alpha \frac{\partial T}{\partial t}$ is the transient term, representing the change in temperature over time.
- $\nabla \cdot \left(-D \nabla T\right)$ is the diffusion term, representing the spread of temperature due to diffusion.
- $\nabla \cdot (uT)$ is the convection term, representing the transport of temperature due to fluid flow.
- $T = T_{\text{boundary}}$ represents the Dirichlet boundary conditions, where $T_{\text{boundary}}$ is the prescribed temperature on these boundaries.
- $-D \nabla T \cdot \mathbf{n} + h (T - T_{\infty}) = q$ represents the Robin boundary conditions, where $\mathbf{n}$ is the normal vector to the boundary, $h$ is the heat transfer coefficient, $T_{\infty}$ is the external temperature (e.g., heat source), and $q$ is the heat flux.

## Further explanation

Sure! Let's break down the equation:

$$
\alpha \frac{\partial T}{\partial t} + \nabla \cdot \left(-D \nabla T\right) = 0
$$

This is a form of the **transient heat conduction equation** (also known as the **heat equation**) with diffusion. Here's what each term represents:

### Terms in the Equation

1. **$\alpha \frac{\partial T}{\partial t}$**:
   - $\alpha$: This is a constant that represents the thermal diffusivity of the material. It combines the effects of thermal conductivity, density, and specific heat capacity.
   - $\frac{\partial T}{\partial t}$: This is the partial derivative of temperature $T$ with respect to time $t$. It represents the rate of change of temperature over time.

2. **$\nabla \cdot \left(-D \nabla T\right)$**:
   - $\nabla T$: This is the gradient of the temperature field $T$. It represents the rate and direction of change of temperature in space.
   - $D$: This is the diffusion coefficient, which measures how quickly heat diffuses through the material.
   - $-D \nabla T$: This term represents the heat flux due to diffusion. The negative sign indicates that heat flows from regions of higher temperature to regions of lower temperature.
   - $\nabla \cdot \left(-D \nabla T\right)$: This is the divergence of the heat flux. It represents the net rate of heat flow out of a point in the material.

### Physical Interpretation

The equation describes how temperature $T$ changes over time and space within a material due to heat conduction. 

- The term $\alpha \frac{\partial T}{\partial t}$ represents the change in temperature over time.
- The term $\nabla \cdot \left(-D \nabla T\right)$ represents the spatial distribution of heat flow due to diffusion.

### Why is it Important?

This equation is fundamental in heat transfer analysis. It helps predict how temperature evolves in a material over time, which is crucial for applications in engineering, physics, and environmental science.

### Summary

- **$\alpha \frac{\partial T}{\partial t}$**: Temporal change in temperature.
- **$\nabla \cdot \left(-D \nabla T\right)$**: Spatial distribution of heat flow.
- **Equation**: Describes transient heat conduction in a material.

I hope this explanation helps! Let me know if you have any more questions.


# Divegence theorem

The general divergence theorem, also known as Gauss's theorem, relates the flux of a vector field through a closed surface to the divergence of the field inside the surface.

$$
\int_{\partial \Omega} \mathbf{F} \cdot \mathbf{n} \, dS = \int_{\Omega} \nabla \cdot \mathbf{F} \, dV
$$

where:
- $ \Omega $ is the volume (domain) enclosed by the surface $ \partial \Omega $.
- $ \mathbf{F} $ is the vector field.
- $ \mathbf{n} $ is the outward-pointing unit normal vector on the surface $ \partial \Omega $.
- $ dS $ is the surface element.
- $ dV $ is the volume element.
- $ \nabla \cdot \mathbf{F} $ is the divergence of the vector field $ \mathbf{F} $.

This theorem essentially states that the total flux of $ \mathbf{F} $ through the boundary $ \partial \Omega $ is equal to the integral of the divergence of $ \mathbf{F} $ over the volume $ \Omega $.
