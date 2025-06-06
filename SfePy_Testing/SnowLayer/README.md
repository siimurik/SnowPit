
# Weak Formulation of the 1D Transient Heat Equation

We start with the **strong form** of the PDE and boundary conditions:

$$
\left\lbrace
\begin{aligned}
&\frac{\partial T}{\partial t}(x, t) = \alpha \frac{\partial^2 T}{\partial x^2}(x, t) && \text{(Heat equation)} \\
&\text{Initial condition:} \\
&T(x, 0) = T_0(x) && \forall x \in [0, L] \\
&\text{Boundary conditions:} \\
&-k \frac{\partial T}{\partial x}(0, t) = h_o \big(T(0, t) - T_o(t)\big) && \text{(Robin BC at $x=0$)} \\
&-k \frac{\partial T}{\partial x}(L, t) = h_i \big(T(L, t) - T_i(t)\big) && \text{(Robin BC at $x=L$)} \\
&\text{where:} \\
&\alpha = \frac{k}{\rho c_p} \quad \text{(Thermal diffusivity)}
\end{aligned}
\right.
$$

---

### Step 1: Rewrite the PDE in divergence form
First, express the heat equation in terms of heat flux $ q = -k \frac{\partial T}{\partial x} $:

$$
\rho c_p \frac{\partial T}{\partial t} = \frac{\partial}{\partial x} \left( k \frac{\partial T}{\partial x} \right)
$$

---

### Step 2: Multiply by a test function $ s(x) $ and integrate
We introduce a test function $ s(x) \in H^1([0, L]) $ (i.e., $ s $ is smooth and vanishes where Dirichlet conditions are imposed, but here we have Robin BCs, so no strict vanishing is required). Multiply the PDE by $ s $ and integrate over the domain $ \Omega = [0, L] $:

$$
\int_0^L \rho c_p \frac{\partial T}{\partial t} s \, dx = \int_0^L \frac{\partial}{\partial x} \left( k \frac{\partial T}{\partial x} \right) s \, dx
$$

---

### Step 3: Apply integration by parts to the right-hand side
Using integration by parts on the right-hand side:

$$
\int_0^L \frac{\partial}{\partial x} \left( k \frac{\partial T}{\partial x} \right) s \, dx = \left[ k \frac{\partial T}{\partial x} s \right]_0^L - \int_0^L k \frac{\partial T}{\partial x} \frac{\partial s}{\partial x} \, dx
$$

Substitute the boundary conditions at $ x = 0 $ and $ x = L $:

$$
\left[ k \frac{\partial T}{\partial x} s \right]_0^L = -h_o (T(0, t) - T_o(t)) s(0) + h_i (T(L, t) - T_i(t)) s(L)
$$

---

### Step 4: Combine terms into the weak form
Now, substitute back into the equation:

$$
\int_0^L \rho c_p \frac{\partial T}{\partial t} s \, dx + \int_0^L k \frac{\partial T}{\partial x} \frac{\partial s}{\partial x} \, dx = -h_o (T(0, t) - T_o(t)) s(0) + h_i (T(L, t) - T_i(t)) s(L)
$$

---

### Final Weak Formulation
The **weak form** is: Find $ T(x, t) $ such that for all test functions $ s(x) $,

$$
\boxed{
\begin{aligned}
\int_0^L \rho c_p \frac{\partial T}{\partial t} s \, dx + \int_0^L k \frac{\partial T}{\partial x} \frac{\partial s}{\partial x} \, dx + h_o T(0, t) s(0) - h_i T(L, t) s(L) \\
= h_o T_o(t) s(0) - h_i T_i(t) s(L)
\end{aligned}
}
$$

### Explanation of Terms:
1. **Time derivative term**:  
   $$
   \int_0^L \rho c_p \frac{\partial T}{\partial t} s \, dx
   $$  
   Represents the transient (time-dependent) heat storage.

2. **Diffusion term**:  
   $$
   \int_0^L k \frac{\partial T}{\partial x} \frac{\partial s}{\partial x} \, dx
   $$  
   Represents heat conduction (spatial derivatives).

3. **Boundary terms (Robin conditions)**:  
   $$
   h_o T(0, t) s(0) - h_i T(L, t) s(L)
   $$  
   Account for convective heat transfer at the boundaries.

4. **Source terms from BCs**:  
   $$
   h_o T_o(t) s(0) - h_i T_i(t) s(L)
   $$  
   Represent the external thermal forcing (ambient temperatures $ T_o $ and $ T_i $).

---

### Comparison with Your Reference Example
Your reference weak form (from `sfepy`) is more general (for 2D/3D and multiple boundary types). The above is its **1D specialization** for:
- Robin BCs at both ends (no Neumann/Dirichlet).
- No internal heat sources.
- Constant material properties (but easily extended to $ k(x) $, $ \rho(x) $, etc.).

This matches the structure of the original PDE system but in integral form for finite element analysis.

---

# Weak Formulation in `sfepy`-like Syntax
```python
'Heat': r"""
dw_dot.i.Omega(mat.rho_cp, s, dT/dt)
+ dw_laplace.i.Omega(mat.k, s, T)
= - dw_bc_newton.i.Gamma_Left(h_o, T_o, s, T)
  - dw_bc_newton.i.Gamma_Right(h_i, T_i, s, T)
"""
```

### Explanation of Terms:
1. **`dw_dot.i.Omega(mat.rho_cp, s, dT/dt)`**  
   - Represents the time-derivative term:  
     $$
     \int_0^L \rho c_p \frac{\partial T}{\partial t} s \, dx
     $$  
   - `mat.rho_cp` = $\rho c_p$ (density × specific heat).  

2. **`dw_laplace.i.Omega(mat.k, s, T)`**  
   - Represents the diffusion term:  
     $$
     \int_0^L k \frac{\partial T}{\partial x} \frac{\partial s}{\partial x} \, dx
     $$  
   - `mat.k` = $k$ (thermal conductivity).  

3. **`dw_bc_newton.i.Gamma_Left(h_o, T_o, s, T)`**  
   - Implements the Robin BC at $x = 0$:  
     $$
     -h_o (T(0, t) - T_o) s(0)
     $$  
   - `h_o` = outer heat transfer coefficient, `T_o` = ambient temperature.  

4. **`dw_bc_newton.i.Gamma_Right(h_i, T_i, s, T)`**  
   - Implements the Robin BC at $x = L$:  
     $$
     -h_i (T(L, t) - T_i) s(L)
     $$  
   - `h_i` = inner heat transfer coefficient, `T_i` = fixed inner temperature.  

---

### Key Adjustments to Match `sfepy`:
- **Naming Conventions**:  
  - `Gamma_Left` → Boundary at $x = 0$.  
  - `Gamma_Right` → Boundary at $x = L$.  
- **Sign Convention**:  
  - `sfepy` typically moves boundary terms to the right-hand side (hence the `-` signs).  
- **Material Parameters**:  
  - `mat.rho_cp` and `mat.k` are material properties (could also be fields for spatially varying cases).  
- **Dynamic vs. Fixed Inputs**:  
  - `h_out_dyn.val` and `T_out_dyn.val` in your example suggest time-dependent inputs, replaced here with generic `h_o` and `T_o` for clarity.  

---

### Full Correspondence with Mathematical Weak Form:
The `sfepy` syntax directly mirrors the weak form derived earlier:
$$
\int_0^L \rho c_p \frac{\partial T}{\partial t} s \, dx + \int_0^L k \frac{\partial T}{\partial x} \frac{\partial s}{\partial x} \, dx = -h_o (T(0, t) - T_o) s(0) - h_i (T(L, t) - T_i) s(L)
$$

### Notes:
1. **Time Discretization**:  
   In `sfepy`, `dT/dt` would be handled by the time-stepping solver (e.g., implicit Euler).  
2. **Boundary Terms**:  
   The `dw_bc_newton` terms assume `sfepy`’s built-in handling of Robin BCs. For custom BCs, you might use `dw_surface_integrate.i.Gamma_Left(h_o * (T - T_o) * s)`.  
3. **Extensions**:  
   - For **spatially varying** $k(x)$ or $h(x)$, replace `mat.k` with a field.  
   - For **Neumann/Dirichlet BCs**, use `dw_bc_dirichlet` or `dw_bc_neumann`.  

---

# **Mathematical equivalency**

### **Form 1 (Compact Form)**
$$
\boxed{
\int_0^L \rho c_p \frac{\partial T}{\partial t} s \, dx + \int_0^L k \frac{\partial T}{\partial x} \frac{\partial s}{\partial x} \, dx = -h_o (T(0, t) - T_o) s(0) - h_i (T(L, t) - T_i) s(L)
}
$$

#### Key Features:
1. **Right-hand side (RHS)** groups all boundary terms.  
2. **Sign convention**:  
   - The negative sign on the RHS comes from moving the boundary terms from the left-hand side (LHS) during integration by parts.  
   - Physically, this represents heat flux **into** the domain (consistent with `sfepy`'s syntax).  

---

### **Form 2 (Expanded Form)**
$$
\boxed{
\begin{aligned}
\int_0^L \rho c_p \frac{\partial T}{\partial t} s \, dx + \int_0^L k \frac{\partial T}{\partial x} \frac{\partial s}{\partial x} \, dx &+ h_o T(0, t) s(0) - h_i T(L, t) s(L) \\
&= h_o T_o(t) s(0) - h_i T_i(t) s(L)
\end{aligned}
}
$$

#### Key Features:
1. **LHS** explicitly separates:  
   - The diffusion term ($\int k \nabla T \cdot \nabla s$).  
   - The boundary terms ($h_o T(0,t)s(0)$, $-h_i T(L,t)s(L)$).  
2. **RHS** contains only the known external forcing terms ($h_o T_o s(0)$, $-h_i T_i s(L)$).  

---

# **Why Are They Equivalent?**
1. **Expand Form 1's RHS**:  
   $$
   -h_o (T(0,t) - T_o)s(0) - h_i (T(L,t) - T_i)s(L) = -h_o T(0,t)s(0) + h_o T_o s(0) - h_i T(L,t)s(L) + h_i T_i s(L)
   $$  
   Rearrange:  
   $$
   -h_o T(0,t)s(0) - h_i T(L,t)s(L) + h_o T_o s(0) + h_i T_i s(L)
   $$  

2. **Move the boundary terms to the LHS**:  
   $$
   \text{LHS} + h_o T(0,t)s(0) - h_i T(L,t)s(L) = h_o T_o s(0) - h_i T_i s(L)
   $$  
   This matches **Form 2 exactly**.

---

### **Key Differences:**
| Feature               | Form 1 (Compact)                          | Form 2 (Expanded)                          |
|-----------------------|------------------------------------------|--------------------------------------------|
| **Boundary Terms**    | On RHS (grouped with forcing)            | Split between LHS (unknown $T$) and RHS (known $T_o, T_i$) |
| **Physical Meaning**  | Net heat flux into domain                | Explicit separation of reactive ($h T$) and forcing ($h T_o$) terms |
| **Implementation**    | Used in `sfepy` (e.g., `dw_bc_newton`)  | Common in theoretical derivations          |

---

### **Which Form to Use?**
1. **For `sfepy`/FEM codes**:  
   Use **Form 1** (matches `dw_bc_newton` syntax):  
   ```python
   'Heat': """
   dw_dot.i.Omega(mat.rho_cp, s, dT/dt) 
   + dw_laplace.i.Omega(mat.k, s, T) 
   = - dw_bc_newton.i.Gamma_Left(h_o, T_o, s, T) 
     - dw_bc_newton.i.Gamma_Right(h_i, T_i, s, T)
   """
   ```

2. **For theoretical clarity**:  
   Use **Form 2** (clearly separates unknown $T$ and known $T_o, T_i$).

---

### **Summary**
Both forms are correct, but:  
- **Form 1** is better for computational implementation (directly maps to FEM weak forms).  
- **Form 2** is better for analytical derivations (explicitly splits terms).  

The choice depends on whether you’re writing code or doing math by hand!

The apparent discrepancy arises from **dimensionality** and **how boundary integrals are represented** in 1D vs. higher dimensions. Let me clarify why the Robin BC terms lose their integral signs in 1D and how this connects to the general weak form.

---

### **1. General Weak Form (2D/3D)**
In higher dimensions, boundary terms are **surface integrals** over $\Gamma$:
$$
\int_\Gamma \underbrace{(-k \nabla T \cdot n)}_{\text{Boundary flux}} s \, dS
$$
For Robin BCs ($-k \nabla T \cdot n = h(T - T_\infty)$), this becomes:
$$
\int_\Gamma h(T - T_\infty) s \, dS
$$
- The integral is necessary because $\Gamma$ is a **region** (surface/curve).  
- Example from your `sfepy` reference:
  $$
  \int_{\Gamma_{source}} h(T - T_\infty) s \, dS + \int_{\Gamma_{plate}} h(T - T_\infty) s \, dS
  $$

---

### **2. 1D Weak Form (No Boundary Integrals)**
In 1D, the "boundary" reduces to **two points** ($x = 0$ and $x = L$). The boundary "integral" simplifies to **point evaluation**:
$$
\int_{\Gamma} (-k \nabla T \cdot n) s \, dS \rightarrow \left[ (-k \partial_x T) s \right]_0^L
$$
- For Robin BCs at $x = 0$ and $x = L$:
  $$
  \left[ (-k \partial_x T) s \right]_0^L = -h_o(T(0,t) - T_o)s(0) - h_i(T(L,t) - T_i)s(L)
  $$
- **No integrals appear** because:
  - In 1D, $\Gamma = \{0, L\}$ is a zero-dimensional set (points).  
  - The "integral" over $\Gamma$ is just a **sum of evaluations** at these points.

---

### **Why the Difference?**
| Dimension | Boundary Integral Representation          | Example                                      |
|-----------|------------------------------------------|---------------------------------------------|
| **2D/3D** | $\int_\Gamma h(T - T_\infty) s \, dS$  | $\int_{\Gamma_{source}} h(T - T_\infty) s$ |
| **1D**    | $h(T - T_\infty) s \big|_{x=0}^L$      | $-h_o(T(0,t) - T_o)s(0) - h_i(T(L,t) - T_i)s(L)$ |

---

### **Key Implications**
1. **Mathematical Consistency**:  
   - Both forms are correct. The 1D case is a specialization where boundary integrals collapse to point evaluations.  
   - In FEM codes like `sfepy`, the 1D version would still use `dw_bc_robin`, but internally, it handles the point evaluation.

2. **Your Numerical Instability**:  
   - If you mistakenly applied **2D-style integrals** to 1D (e.g., multiplying by a "boundary length"), the scaling would be wrong, causing energy imbalance.  
   - The **signs** must still be negative for Robin BCs to ensure heat flows **out** when $T > T_\infty$.

3. **Spectral/Finite Difference Methods**:  
   - In methods that don’t use weak forms, Robin BCs are enforced directly as equations at boundary points (no integrals).  

---

### **How to Reconcile the Two Forms**
#### **General Weak Form (N-D)**
$$
\int_\Omega \rho c_p \frac{\partial T}{\partial t} s \, dV + \int_\Omega k \nabla T \cdot \nabla s \, dV = - \int_\Gamma h(T - T_\infty) s \, dS
$$

#### **1D Weak Form**
$$
\int_0^L \rho c_p \frac{\partial T}{\partial t} s \, dx + \int_0^L k \partial_x T \partial_x s \, dx = -h_o(T(0,t) - T_o)s(0) - h_i(T(L,t) - T_i)s(L)
$$

#### **Connection**
- In 1D, $\int_\Gamma \rightarrow \sum_{x \in \{0, L\}}$.  
- The **integral signs vanish** because the boundary is a discrete set.

---

### **Practical Advice**
1. **For 1D Codes**:  
   - Implement Robin BCs as **point conditions** (no integrals).  
   - Ensure the negative sign is preserved:
     ```python
     rhs = -h_left * (T[0] - T_o) * s[0] - h_right * (T[-1] - T_i) * s[-1]
     ```

2. **For `sfepy`/FEM**:  
   - Use `dw_bc_robin` (it automatically handles dimensionality).  
   - The integral is hidden in the assembly process.

3. **Debugging Tip**:  
   - For 1D, compare with the analytical solution for a rod with Robin BCs:
     $$
     T(x) = T_\infty + (T_0 - T_\infty) e^{-h x / k}
     $$
     If your numerical solution diverges, check:
     - **Signs** on Robin terms.  
     - **Scaling** (no integrals in 1D!).  

---

### **Summary**
- **In 2D/3D**: Robin BCs are surface integrals ($\int_\Gamma$).  
- **In 1D**: They reduce to point evaluations (**no integrals**).  
- **Signs are critical**: Negative for correct energy outflow.  
- Your instability likely came from either:
  1. Missing negative signs, or  
  2. Incorrectly treating 1D BCs as integrals.  

The two forms are consistent—the 1D version is just the degenerate case of the general N-D weak form.