%{
============================================================================================================
Transient diffusion equation
PDE and boundary conditions
    $$
    \begin{cases}
    \alpha \frac{\partial T}{\partial t} + \nabla \cdot \left(-D \nabla T\right) = 0 & \text{in the domain}, \\
    T = T_{\text{boundary}} & \text{on Dirichlet boundaries}, \\
    -D \nabla T \cdot \mathbf{n} + h (T - T_{\infty}) = q & \text{on Robin boundaries}.
    \end{cases}
    $$
============================================================================================================
%}

clc
clear

%% Define the domain and create a mesh structure
L = 50;  % domain length
Nx = 11; % number of cells
m = createMesh3D(Nx, Nx, Nx, L, L, L);

%% Create the boundary condition structure
BC = createBC(m); % all Neumann boundary condition structure
BC.left.a(:) = 0; BC.left.b(:) = 1; BC.left.c(:) = 0; % Dirichlet on the left
BC.right.a(:) = 0; BC.right.b(:) = 1; BC.right.c(:) = 0; % Dirichlet on the right
BC.front.a(:) = 0; BC.front.b(:) = 1; BC.front.c(:) = 0; % Dirichlet on the front
BC.back.a(:) = 0; BC.back.b(:) = 1; BC.back.c(:) = 0; % Dirichlet on the back
BC.bottom.a(:) = 0; BC.bottom.b(:) = 1; BC.bottom.c(:) = 4; % Dirichlet on the bottom (constant temperature)

%% Heated isolated boundary on top
% Define the thickness of the insulating layer (in meters)
insulation_thickness = 0.1; % example thickness

% Define the thermal conductivity of the insulating material (W/m·K)
k_insulation = 0.04; % example value for typical insulation material

% Calculate the heat transfer coefficient (h = k / thickness)
h_insulation = k_insulation / insulation_thickness;

% Robin boundary condition on the top (insulated with heat source)
heat_source_value = 10; % example heat source value representing solar heat
BC.top.a(:) = 1; % Neumann condition: true 
BC.top.b(:) = h_insulation; 
BC.top.c(:) = heat_source_value; % Robin boundary condition on top

%% Define the transfer coefficients
D_val = 0.1;
D = createCellVariable(m, D_val);
alfa = createCellVariable(m, 1);
%u = createFaceVariable(m, [0,0,0.5]);

%% Define initial values with optional temperature gradient
add_grad = true;   % <--- add temp gradient to sides if set to 'true'
if add_grad == false
    c_init = 0; % starting temperature of 0 on all faces
else
    c_init = get_ic(m); % starting temperature with gradient
end

c_old = createCellVariable(m, c_init, BC); % initial values
c = c_old; % assign the old value of the cells to the current values

%% Loop
Dave = harmonicMean(D);
Mdiff = diffusionTerm(Dave);
[Mbc, RHSbc] = boundaryCondition(BC);
FL = fluxLimiter('Superbee');
%Mconv = convectionTvdTerm(u, c, FL);
dt = 1; % time step
final_t = 100;

for t = dt:dt:final_t
    [M_trans, RHS_trans] = transientTerm(c_old, dt, alfa);
    M = M_trans - Mdiff + Mbc;% + Mconv;
    RHS = RHS_trans + RHSbc;
    c = solvePDE(m, M, RHS);
    c_old = c;
    figure(1); 
    visualizeCells(c); 
    colorbar; % add colorbar for better visualization
    clim([0 4]); % set color axis limits
    drawnow;
end

%% Visualization
figure(1); 
visualizeCells(c); 
colorbar; % add colorbar for better visualization
clim([0 4]); % set color axis limits

%% Temperature Gradient Function
function temp = side_temp(y)
    surface_temp = 0.0;  % Temperature at the top (°C)
    bottom_temp = 4.0;   % Temperature at the bottom (°C)
    height = 50;        % Height of the domain
    temp = surface_temp + (bottom_temp - surface_temp) * (y / height);
end

%% Initial Condition Function
function c_init = get_ic(m)
    % Get the coordinates of the cell centers
    [~, Y, ~] = ndgrid(m.cellcenters.x, m.cellcenters.y, m.cellcenters.z);
    c_init = side_temp(Y);
end