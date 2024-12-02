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

This adjusted system combines all terms to represent the temperature distribution with the specified boundary conditions. If you need further modifications or have any other questions, feel free to ask!