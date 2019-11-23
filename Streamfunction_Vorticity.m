%% SIMPLE

% Code written by Jay Mehta (July 2019).

%% Clear Everything
clear all
clc
rmpath('Source')
rmpath('Source/dst_idst')

%% Add Paths
addpath('Source')
addpath('Source/dst_idst')

%% Set up the problem domain and the problem object

% Set up the parameters that have to be passed

params = flow_parameters_init;
domain = domain_parameters_init;

% Domain
Nx = 64;
Ny = 64;
domain.Nx = Nx;
domain.Ny = Ny;

x_range = [0 1];
y_range = [0 1];
domain.x_range = x_range;
domain.y_range = y_range;

dx = (x_range(2)-x_range(1))/Nx;
params.dx = dx;
dy = (y_range(2)-y_range(1))/Ny;

%% Setting up the Domain of the Problem

[X_n, Y_n] = DomainSetup(params,domain,"node");
[X_e_x, Y_e_x] = DomainSetup(params,domain,"xe");
[X_e_y, Y_e_y] = DomainSetup(params,domain,"ye");
[X_c,Y_c] = DomainSetup(params,domain,"cell");

domain.X_n = X_n;
domain.Y_n = Y_n;
domain.X_e_x = X_e_x;
domain.Y_e_x = Y_e_x;
domain.X_e_y = X_e_y;
domain.Y_e_y = Y_e_y;
domain.X_c = X_c;
domain.Y_c = Y_c;

%% Initializing the flow parameters

U = 1;
V = 0;
params.U = U;

Re = 400;
nu = U / Re;
params.nu = nu;

%% Pre-Allocate Memory

gamma = NodeData(Nx,Ny);
sf = NodeData(Nx,Ny);
velocity = EdgeData(Nx,Ny);

a_p = EdgeData(Nx,Ny);
a_e = EdgeData(Nx,Ny);
a_w = EdgeData(Nx,Ny);
a_n = EdgeData(Nx,Ny);
a_s = EdgeData(Nx,Ny);
d = EdgeData(Nx,Ny);

% A_p = CellData(Nx,Ny);
% A_e = CellData(Nx,Ny);
% A_w = CellData(Nx,Ny);
% A_n = CellData(Nx,Ny);
% A_s = CellData(Nx,Ny);

% p = CellData(Nx,Ny);
% p_prime = CellData(Nx,Ny);
% 
% rhs = CellData(Nx,Ny);

% rf_p = 5e-1;
rf_g = 5e-1;

max_gs_iter = 5;

conv_x = 1;
conv_y = 1;
tol = 1e-6;

%% Apply Boundary Conditions

uB = zeros(1,Nx+1);
uT = ones(1,Nx+1);
uL = zeros(1,Ny+2);
uR = zeros(1,Ny+2);

vB = zeros(1,Nx+2);
vT = zeros(1,Nx+1);
vL = zeros(1,Ny+1);
vR = zeros(1,Ny+1);

bc = bc_params(uL,uR,uB,uT,vL,vR,vB,vT);

velocity = apply_bc(bc,velocity);

%% Solving

for iter = 1:100
    
    % Boundary conditions for vorticity are nuanced.
    for i = 2:Nx
        gamma.x(i,1) = -2/dx*velocity.x(i,2);
        gamma.x(i,Ny+1) = -2/dx*(uT(i) - velocity.x(i,Ny+1));
    end
    for j = 2:Ny
        gamma.x(1,j) = 2/dx*velocity.y(2,j);
        gamma.x(Nx+1,j) = -2/dx*velocity.y(Nx+1,j);
    end
    
    gamma_old = gamma;
    
    velocity_x_c = interpol(velocity,CellData(Nx,Ny),1);
    velocity_u_v = interpol(velocity_x_c,EdgeData(Nx,Ny),2); % U in V's location
    velocity_y_c = interpol(velocity,NodeData(Nx,Ny),2);
    velocity_v_u = interpol(velocity_y_c,EdgeData(Nx,Ny),1); % V in U's location
    
    for i = 2:Nx
        for j = 2:Ny
            a_p.x(i,j) = (0.5 * dx * (velocity_u_v.y(i+1,j) - velocity_u_v.y(i,j))...
                        + 0.5 * dx * (velocity_v_u.x(i,j) - velocity_v_u.x(i,j-1))...
                        + 4*nu)/rf_g;
            a_e.x(i,j) = -0.5 * dx * velocity_u_v.y(i+1,j) + nu;
            a_w.x(i,j) = 0.5 * dx * velocity_u_v.y(i,j) + nu;
            a_n.x(i,j) = -0.5 * dx * velocity_v_u.x(i,j) + nu;
            a_s.x(i,j) = 0.5 * dx * velocity_v_u.x(i,j-1) + nu;
            
        end
    end
    
    gamma_star = gamma_old; % Why not use a zero-prediction field?
    
    for gs_iter = 1:max_gs_iter
        for i = 2:Nx
            for j = 2:Ny
                gamma_star.x(i,j) = (a_e.x(i,j) * gamma_star.x(i+1,j)...
                    + a_w.x(i,j) * gamma_star.x(i-1,j)...
                    + a_n.x(i,j) * gamma_star.x(i,j+1)...
                    + a_s.x(i,j) * gamma_star.x(i,j-1)...
                    + (1-rf_g) * a_p.x(i,j) * gamma_old.x(i,j)) / a_p.x(i,j);
            end
        end
    end

    gs_conv(1) = max(max(abs(gamma_star.x-gamma_old.x)));        
    gs_conv = max(gs_conv);
    gamma = gamma_star;
    
    for i = 2:Nx
        gamma.x(i,1) = -2/dx*velocity.x(i,2);
        gamma.x(i,Ny+1) = -2/dx*(uT(i) - velocity.x(i,Ny+1));
    end
    for j = 2:Ny
        gamma.x(1,j) = 2/dx*velocity.y(2,j);
        gamma.x(Nx+1,j) = -2/dx*velocity.y(Nx+1,j);
    end
    gamma.x = -gamma.x*dx^2;
    
    sf = smoothing(sf,gamma,"sor","tol",1e-3);
    sf.x(1,:) = 0;
    sf.x(end,:) = 0;
    sf.x(:,1) = 0;
    sf.x(:,end) = 0;
    
    gamma.x = -gamma.x/dx^2;
    
    velocity = curl_2(sf);
    velocity.x = velocity.x/dx;
    velocity.y = velocity.y/dy;
    velocity = apply_bc(bc,velocity);
end