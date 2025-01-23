# Heat Transfer Simulation in a Frustum Cone Using SfePy

This repository contains an example of solving a heat transfer problem in a frustum cone geometry using [SfePy](https://sfepy.org), a finite element analysis library. The simulation models transient heat conduction with boundary conditions representing convective heat loss and a fixed temperature at the bottom.

---

## Problem Definition

The heat transfer problem is governed by the following equations:

### Governing Equation
$$
\begin{aligned}
\frac{\partial T}{\partial t} &= \alpha \nabla^2 T  && \forall x \in \Omega, \forall t \in [t_0, t_1],
\end{aligned}
$$
where:
- $T$ is the temperature in the domain ( $^\circ\text{C}$ ),
- $\alpha$ is the thermal diffusivity ( $\text{m}^2/\text{s}$ ),
- $\Omega$ is the domain of the frustum cone.

### Boundary Conditions
1. **Top Surface ($\Gamma_{\text{top}}$)**  
   Convective heat transfer and heat flux at the top boundary:
   $$
   \begin{aligned}
   -\alpha \nabla T \cdot n &= q(t) + h_{\text{top}}(T - T_{\text{top\_inf}}) && \forall x \in \Gamma_{\text{top}},
   \end{aligned}
   $$
   where:
   - $q(t)$ is the heat flux ($\text{W/m}^2$),
   - $h_{\text{top}}$ is the heat transfer coefficient ($\text{W/m}^2/\text{K}$),
   - $T_{\text{top\_inf}}$ is the ambient temperature at the top ( $^\circ\text{C}$ ).

2. **Bottom Surface ($\Gamma_{\text{bottom}}$)**  
   Fixed temperature boundary condition:
   $$
   \begin{aligned}
   T &= 2.0 && \forall x \in \Gamma_{\text{bottom}}.
   \end{aligned}
   $$

3. **Side Surfaces ($\Gamma_{\text{side}}$)**  
   Insulated boundary condition (no heat flux):
   $$
   \begin{aligned}
   \nabla T \cdot n &= 0 && \forall x \in \Gamma_{\text{side}}.
   \end{aligned}
   $$

### Initial Condition
The initial temperature distribution is defined as:
$$
\begin{aligned}
T(x, 0) &= T_0 && \forall x \in \Omega,
\end{aligned}
$$
where $T_0 = -2.0  ^\circ\text{C}$.

# Thermal diffusion coefficient

The diffusion coefficient in the heat equation, often denoted as $ D $ or $ \alpha $, is a measure of how quickly heat spreads through a material. It is determined by the material's properties and can be calculated using the formula:

$$
\alpha = \frac{k}{\rho c_p}
$$

where:
- $ k $ is the thermal conductivity of the material,
- $ \rho $ is the density of the material,
- $ c_p $ is the specific heat capacity at constant pressure[1](https://www.tec-science.com/thermodynamics/heat/heat-equation-diffusion-equation/)[2](https://www.uni-muenster.de/imperia/md/content/physik_tp/lectures/ws2016-2017/num_methods_i/heat.pdf).

This coefficient appears in the heat equation, which in one dimension is written as:

$$
\frac{\partial T}{\partial t} = \alpha \frac{\partial^2 T}{\partial x^2}
$$

Here, $ T $ is the temperature, $ t $ is time, and $ x $ is the spatial coordinate[2](https://www.uni-muenster.de/imperia/md/content/physik_tp/lectures/ws2016-2017/num_methods_i/heat.pdf).

# Mass Diffusion vs Thermal Diffusion

The diffusion coefficient is sometimes denoted as $ D $ and sometimes as $ \alpha $ depending on the context and the specific type of diffusion being described.

- **$ D $**: This symbol is commonly used in the context of **mass diffusion** (e.g., the diffusion of particles, molecules, or gases). It represents the diffusivity or diffusion coefficient in Fick's laws of diffusion.

- **$ \alpha $**: This symbol is typically used in the context of **thermal diffusion** (heat conduction). It represents the thermal diffusivity in the heat equation.

Both symbols essentially describe how quickly something spreads out, but they are used in different fields to avoid confusion. $ D $ is used for mass diffusion, while $ \alpha $ is used for thermal diffusion.

Does that help clarify things? Feel free to ask if you have more questions!
---

## Simulation Setup

### Geometry
The computational domain is a frustum cone, described by the mesh file located at:
```
C:\Users\sipuga\Documents\SnowStorageSolvers\SfePy_Testing\SnowPit\FrustumCone\frustum_cone.mesh
```

### Material Properties
- Thermal diffusivity: $\alpha = 0.01 \, \text{m}^2/\text{s}$.

### Boundary Heat Loss
- Heat transfer coefficient: $h = 10.0 \, \text{W/m}^2/\text{K}$,
- Ambient temperature: $T_{\text{inf}} = -2.0  ^\circ\text{C}$.

### Numerical Parameters
- Time interval: $t \in [0.0, 0.1] \, \text{s}$,
- Number of time steps: $n = 101$.

### Regions
- **Omega**: Entire domain,
- **Gamma\_Top**, **Gamma\_Bottom**, **Gamma\_Side**: Named boundaries for applying specific conditions.

---

## How to Run

1. Install [SfePy](https://sfepy.org/doc-devel/installation.html) and its dependencies.
2. Save the script and the mesh file in your working directory.
3. Execute the script using:
   ```bash
   sfepy-run heat_frustum_cone.py
   ```

# Documentation

Here’s a refined and more correct version of the **first formulation**, ensuring consistency between the governing equations, boundary conditions, and notation:

$$
\begin{aligned}
\frac{\partial T}{\partial t} &= \alpha \nabla^2 T  && \forall (x,y,z) \in \Omega, \forall t \in [t_0, t_1], \\
-\alpha \nabla T \cdot n &= q(t) + h_{top}(T - T_{\text{top\_inf}}) && \forall (x,y,z) \in \Gamma_{\text{top}}, \\
T &= 2.0 && \forall (x,y,z) \in \Gamma_{\text{bottom}}, \\
\nabla T \cdot n &= 0 && \forall (x,y,z) \in \Gamma_{\text{side}}.
\end{aligned}
$$

---

### Changes Made:
1. **Sign Consistency**:
   - The heat diffusion equation $ \frac{\partial T}{\partial t} = \alpha \nabla^2 T $ is written in its standard form (positive $\alpha \nabla^2 T$).
   - The boundary flux condition on $\Gamma_{\text{top}}$ is expressed as $-\alpha \nabla T \cdot n$, which is the outward diffusive heat flux.

2. **Boundary Flux on $\Gamma_{\text{top}}$**:
   - The term $-\alpha \nabla T \cdot n = q(t) + h_{top}(T - T_{\text{top\_inf}})$ is now explicitly consistent with the heat flux condition. If (\alpha \nabla T \cdot n) is negative, heat is entering the object (inward flux). Here:
     - $q(t)$ represents any additional flux (e.g., external source/sink).
     - $h_{top}(T - T_{\text{top\_inf}})$ models convective heat transfer.

3. **General Clarity**:
   - The domain ($\Omega$) and boundary ($\Gamma$) are clearly distinguished, with $n$ as the outward normal vector for boundary conditions.
   - Notation for $T_{\text{top\_inf}}$ reflects a consistent subscript style for clarity.

---

### Interpretation:
- The governing equation describes heat diffusion within the domain ($\Omega$).
- At the top boundary ($\Gamma_{\text{top}}$), the heat flux includes contributions from an external source $q(t)$ and convective heat transfer.
- Fixed temperature and insulated boundaries are specified at $\Gamma_{\text{bottom}}$ and $\Gamma_{\text{side}}$, respectively.

---

# Derivation of the weak form

To derive the weak integral form of the given problem step by step, we follow these steps:

---

### **1. Multiply the PDE by a test function and integrate over the domain**
Introduce a test function $s(x)$ from a suitable function space, and multiply both sides of the governing equation by $s$:
$$
s \frac{\partial T}{\partial t} = s \alpha \nabla^2 T \quad \text{in } \Omega.
$$
Integrate over the domain $\Omega$:
$$
\int_\Omega s \frac{\partial T}{\partial t} \, d\Omega = \int_\Omega s \alpha \nabla^2 T \, d\Omega.
$$

---

### **2. Handle the time derivative term**
The term involving the time derivative can be rewritten as:
$$
\int_\Omega s \frac{\partial T}{\partial t} \, d\Omega.
$$
This is the time-dependent contribution and remains as-is in the weak form.

---

### **3. Apply integration by parts to the second-order spatial term**
For the Laplacian term $\nabla^2 T$, use the vector calculus identity and integration by parts:
$$
\int_\Omega s \alpha \nabla^2 T \, d\Omega = -\int_\Omega \alpha \nabla s \cdot \nabla T \, d\Omega + \int_{\partial \Omega} s \alpha \nabla T \cdot n \, d\Gamma,
$$
where $n$ is the outward unit normal on the boundary $\partial \Omega$.

- The first term $-\int_\Omega \alpha \nabla s \cdot \nabla T \, d\Omega$ becomes part of the weak formulation.
- The second term $\int_{\partial \Omega} s \alpha \nabla T \cdot n \, d\Gamma$ involves boundary contributions.

---

### **4. Apply boundary conditions**
Substitute the boundary conditions to simplify the boundary integral term:
1. On $\Gamma_{\text{top}}$, $-\alpha \nabla T \cdot n = q(t) + h_{\text{top}}(T - T_{\text{top\_inf}})$:
   $$
   \int_{\Gamma_{\text{top}}} s (-\alpha \nabla T \cdot n) \, d\Gamma = \int_{\Gamma_{\text{top}}} s \big(q(t) + h_{\text{top}}(T - T_{\text{top\_inf}})\big) \, d\Gamma.
   $$

2. On $\Gamma_{\text{bottom}}$, $T = 2.0$:
   This is a Dirichlet condition. In the weak form, the test function $s$ must satisfy $s = 0$ on $\Gamma_{\text{bottom}}$, so this term does not contribute directly to the weak form.

3. On $\Gamma_{\text{side}}$, $\nabla T \cdot n = 0$:
   This implies no flux, so the corresponding integral over $\Gamma_{\text{side}}$ vanishes.

---

### **5. Combine terms into the weak integral form**
The weak form is obtained by combining all terms:
$$
\int_\Omega s \frac{\partial T}{\partial t} \, d\Omega 
- \int_\Omega \alpha \nabla s \cdot \nabla T \, d\Omega
+ \int_{\Gamma_{\text{top}}} s \big(q(t) + h_{\text{top}}(T - T_{\text{top\_inf}})\big) \, d\Gamma = 0.
$$

---

### **6. Incorporate initial and boundary conditions**
- Initial condition: $T(x, 0) = T_0$ in $\Omega$.
- Dirichlet boundary condition ($\Gamma_{\text{bottom}}$): Enforced by modifying the function space for $T$, where $T$ satisfies the fixed value $T = 2.0$ on $\Gamma_{\text{bottom}}$.

---

### **Final Weak Form**
The weak form of the problem is:
$$
\begin{aligned}
\int_\Omega s \frac{\partial T}{\partial t} \, d\Omega 
- \int_\Omega \alpha \nabla s \cdot \nabla T \, d\Omega 
+ \int_{\Gamma_{\text{top}}} s \big(q(t) + h_{\text{top}}(T - T_{\text{top\_inf}})\big) \, d\Gamma &= 0, \\
T &= 2.0 && \text{on } \Gamma_{\text{bottom}}, \\
T(x, 0) &= T_0 && \text{in } \Omega.
\end{aligned}
$$
This form is suitable for numerical approximation using finite elements or other variational methods.

You're on the right track conceptually, but there are some clarifications to consider when implementing Neumann boundary conditions for all the regions. Here’s how to reason about it and why some terms may not be strictly necessary.

---

# **Neumann Boundary Conditions**

1. **Adding Neumann Terms**:
   Neumann boundary conditions specify the flux $q_n = -\alpha \nabla T \cdot n$ on a given boundary. In the weak form, these terms appear as integrals over the respective boundary:
   $$
   \int_{\Gamma} q_n \, s \, d\Gamma.
   $$
   If the flux is zero (homogeneous Neumann), this integral contributes nothing to the equation, so you **do not need to explicitly add terms like `dw_integrate.i.Gamma_X(0.0, s)`**. These terms are mathematically redundant because the integral evaluates to zero.

2. **Including Non-Zero Neumann Flux**:
   For boundaries where the flux is non-zero, like the top boundary ($q_{\text{top}} = \text{flux.val}$), you should include terms like:
   ```python
   + dw_integrate.i.Gamma_Top(flux.val, s)
   ```
   This is already correct in your current implementation.

3. **Robin/Convective Terms**:
   For boundaries with heat loss (Robin condition), such as $\Gamma_{\text{top}}$, you correctly include the `dw_bc_newton` term:
   ```python
   + dw_bc_newton.i.Gamma_Top(heat_loss.h_top, heat_loss.T_top_inf, s, T)
   ```

4. **Homogeneous Neumann Boundaries**:
   For the left, right, and bottom boundaries, if the flux is zero (homogeneous Neumann), **you do not need to explicitly include `dw_integrate` terms**. These boundaries are naturally enforced by the weak form unless explicitly overridden.

---

### **Proposed Equation**
If your goal is to specify:
- Non-zero flux on `Gamma_Top` (Neumann condition and/or Robin condition),
- Zero flux on `Gamma_Left`, `Gamma_Right`, and `Gamma_Bottom` (homogeneous Neumann),

The correct equation should be:
```python
equations = {
    'Temperature': """
    dw_dot.i.Omega( s, dT/dt ) + dw_laplace.i.Omega( m.\alpha, s, T ) =
    + dw_integrate.i.Gamma_Top(flux.val, s)
    + dw_bc_newton.i.Gamma_Top(heat_loss.h_top, heat_loss.T_top_inf, s, T)
    """
}
```

### **Why Simplified Terms Aren't Needed**
Including terms like:
```python
+ dw_integrate.i.Gamma_Left(0.0, s)
+ dw_integrate.i.Gamma_Right(0.0, s)
+ dw_integrate.i.Gamma_Bottom(0.0, s)
```
doesn’t add anything to the equation because:
1. The flux is zero, so the integral evaluates to zero.
2. Homogeneous Neumann conditions are naturally enforced in the weak form by the absence of boundary terms.

---

### **Final Note**
To summarize:
- Explicitly add terms for non-zero flux or Robin conditions.
- Skip adding terms for homogeneous Neumann boundaries unless you explicitly want to document their presence (though it is unnecessary for numerical accuracy).