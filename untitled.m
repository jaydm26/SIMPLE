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

%% Create the L^-1 operator using Lattice Green's function (Liska, Colonius, 2016)

% g_hat = L_inv(domain,"node");

%% Setting up the Object

R = 0.5;
params.char_L = R;

xc = 0;
yc = 0;
body_function = @(x,y) (x-xc)^2 + (y-yc)^2 - (R)^2;

N_theta = floor(2*pi*R/dx);
d_theta = 2*pi/N_theta;

theta_range = 0:d_theta:2*pi-d_theta;

xi = zeros(length(theta_range),1);
eta = zeros(length(theta_range),1);

for i = 1:length(theta_range)
    xi(i) = xc + R*cos(theta_range(i));
    eta(i) = yc + R*sin(theta_range(i));
end

%% Initializing the flow parameters

U = 1;
V = 0;
params.U = U;

Re = 100;
% nu = U * R / Re;
nu = 0.01;
params.nu = nu;

Co = 2.5e-1;
Fo = 5;
dt = min([Fo * dx^2/nu,Co*dx]);
params.dt = dt;

Fo = nu * dt/dx^2;
Co = dt/dx;

params.Fo = Fo;
params.Co = Co;

t_steady = (2*R)^2/nu;
tf = t_steady;
time_range = 0:dt:tf;

%% Pre-Allocate Memory

velocity = EdgeData(Nx,Ny); % Velocity Field
pressure = CellData(Nx,Ny);
Fx = zeros(length(xi),1);
Fy = zeros(length(eta),1);
ub = zeros(length(xi),1); % X-component of Velocity on the body
vb = zeros(length(eta),1); % Y-component of Velocity on the body
Drag = 0;
Lift = 0;

a_p = EdgeData(Nx,Ny);
a_e = EdgeData(Nx,Ny);
a_w = EdgeData(Nx,Ny);
a_n = EdgeData(Nx,Ny);
a_s = EdgeData(Nx,Ny);
d = EdgeData(Nx,Ny);

A_p = CellData(Nx,Ny);
A_e = CellData(Nx,Ny);
A_w = CellData(Nx,Ny);
A_n = CellData(Nx,Ny);
A_s = CellData(Nx,Ny);

p = CellData(Nx,Ny);

rf_p = 7.5e-1;
rf_v = 7.5e-1;

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

while conv_x > tol || conv_y > tol
    
    velocity_old = velocity;
    velocity_x_c = interpol(velocity_old,CellData(Nx,Ny),1);
    velocity_x_n = interpol(velocity_old,CellData(Nx,Ny),1);
    velocity_y_c = interpol(velocity_old,CellData(Nx,Ny),2);
    velocity_y_n = interpol(velocity_old,CellData(Nx,Ny),2);
    
    for i = 2:Nx
        for j = 2:Ny+1
            a_p.x(i,j) = (0.5 * dx * (velocity_x_c.x(i+1,j) - velocity_x_c.x(i,j))...
                        + 0.5 * dx * (velocity_y_n.x(i,j) - velocity_y_n.x(i,j-1))...
                        + 4*nu)/rf_v;
            a_e.x(i,j) = -0.5 * dx * velocity_x_c.x(i+1,j) + nu;
            a_w.x(i,j) = 0.5 * dx* velocity_x_c.x(i,j) + nu;
            a_n.x(i,j) = -0.5 * dx * velocity_y_n.x(i,j) + nu;
            a_s.x(i,j) = 0.5 * dx * velocity_y_n.x(i,j-1) + nu;
            
            d.x(i,j) = dx * 1/a_p.x(i,j);
        end
    end
    
    for i = 2:Nx+1
        for j = 2:Ny
            a_p.y(i,j) = (0.5 * dx * (velocity_y_c.x(i,j+1) - velocity_y_c.x(i,j))...
                        + 0.5 * dx * (velocity_x_n.x(i,j) - velocity_x_n.x(i-1,j))...
                        + 4*nu)/rf_v;
            a_e.y(i,j) = -0.5 * dx * velocity_x_n.x(i,j) + nu;
            a_w.y(i,j) = 0.5 * dx * velocity_x_n.x(i-1,j) + nu;
            a_n.y(i,j) = -0.5 * dx * velocity_y_c.x(i,j+1) + nu;
            a_s.y(i,j) = 0.5 * dx * velocity_y_c.x(i,j) + nu;
            
            d.y(i,j) = dx * 1/a_p.y(i,j);
        end
    end
    
    velocity_star = velocity_old;
    
    for gs_iter = 1:max_gs_iter
        
        for i = 2:Nx
            for j = 2:Ny+1
                velocity_star.x(i,j) = (a_e.x(i,j) * velocity_star.x(i+1,j)...
                    + a_w.x(i,j) * velocity_star.x(i-1,j)...
                    + a_n.x(i,j) * velocity_star.x(i,j+1)...
                    + a_s.x(i,j) * velocity_star.x(i,j-1)...
                    - dx * (p.x(i+1,j) - p.x(i,j))...
                    + (1-rf_v) * a_p.x(i,j) * velocity_old.x(i,j)) / a_p.x(i,j);
            end
        end
        
        
        
        for i = 2:Nx+1
            for j = 2:Ny
                velocity_star.y(i,j) = (a_e.y(i,j) * velocity_star.y(i+1,j)...
                    + a_w.y(i,j) * velocity_star.y(i-1,j)...
                    + a_n.y(i,j) * velocity_star.y(i,j+1)...
                    + a_s.y(i,j) * velocity_star.y(i,j-1)...
                    - dx * (p.x(i,j+1) - p.x(i,j))...
                    + (1-rf_v) * a_p.y(i,j) * velocity_old.y(i,j)) / a_p.y(i,j);
            end
        end
        
        velocity_star = apply_bc(bc,velocity_star);
        
        gs_conv(1) = max(max(abs(velocity_star.x-velocity_old.x)));
        gs_conv(2) = max(max(abs(velocity_star.y-velocity_old.y)));
        
        gs_conv = max(gs_conv);
        
    end
    
    for i = 2:Nx+1
        for j = 2:Ny+1
            A_p.x(i,j) = d.x(i,j) + d.x(i-1,j) + d.y(i,j) + d.y(i,j-1);
            A_e.x(i,j) = d.x(i,j);
            A_w.x(i,j) = d.x(i-1,j);
            A_n.x(i,j) = d.y(i,j);
            A_s.x(i,j) = d.y(i,j-1);
            
            rhs = CellData(Nx,Ny);
            rhs.x(i,j) = -velocity_star.x(i,j) + velocity_star.x(i-1,j)...
                - velocity_star.y(i,j) + velocity_star.y(i,j-1);
        end
    end
    
    p_prime = CellData(Nx,Ny);
    
    for i = 2:Nx+1
        for j = 2:Ny+1
            p_prime.x(i,j) = (A_e.x(i,j) * p_prime.x(i+1,j)...
                + A_w.x(i,j) * p_prime.x(i-1,j)...
                + A_n.x(i,j) * p_prime.x(i,j+1)...
                + A_s.x(i,j) * p_prime.x(i,j-1)...
                + rhs.x(i,j)) / A_p.x(i,j);
        end
    end
    
    p.x = p.x + rf_p * p_prime.x;
    
    velocity_prime = grad(p_prime);
    velocity_prime.x = -velocity_prime.x .* d.x;
    velocity_prime.y = -velocity_prime.y .* d.y;
    
    velocity.x = velocity_star.x + velocity_prime.x;
    velocity.y = velocity_star.y + velocity_prime.y;
    
    conv_x = norm(velocity_prime.x(2:Nx,2:Ny+1))/norm(velocity.x(2:Nx,2:Ny+1));
    conv_y = norm(velocity_prime.y(2:Nx+1,2:Ny))/norm(velocity.y(2:Nx+1,2:Ny));
end

%% Post Processing

gamma = curl_2(velocity);
gamma.x = -gamma.x;

sf = smoothing(NodeData(Nx,Ny),gamma,'sor','tol',1e-3);
