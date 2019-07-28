clear all
clc

Nx= 3;
rho= 1.0;
visc= 0.01;
alphaP= 0.5;
alphaV= 0.5;
maxIter= 10000;  % Maximum number of SIMPLE iterations
maxGSiter= 100;
u_TOP= 1.0;
u_BOTTOM= 0.0;
v_LEFT= 0.0;
v_RIGHT= 0.0;
Ny = Nx;
h = 1/Nx;

u      = zeros(Nx+1,Ny+2);   % x velocity values
uStar  = zeros(Nx+1,Ny+2);   % Temporary x velocity values
uPrime = zeros(Nx+1,Ny+2);   % x velocity correction
v      = zeros(Nx+2,Ny+1);   % y velocity values
vStar  = zeros(Nx+2,Ny+1);   % Temporary y velocity values
vPrime = zeros(Nx+2,Ny+1);
p      = zeros(Nx+2,Ny+2);
pOld   = zeros(Nx+2,Ny+2);
pPrime = zeros(Nx+2,Ny+2);
dU = zeros(Nx+1,Ny+2);
dV = zeros(Nx+2,Ny+1);

u(:,:) = 0.0; v(:,:) = 0.0; p(:,:) = 0.0;
uOld = u;
vOld = v;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start the SIMPLE iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = 1
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Update ghost cell velocities using linear interpolation
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:Nx+1
  uOld(i,1)    = 2*u_BOTTOM - uOld(i,2);
  uOld(i,Ny+2) = 2*u_TOP    - uOld(i,Ny+1);
end
for j = 1:Ny+1
  vOld(1,j)    = 2*v_LEFT  - vOld(2,j);
  vOld(Nx+2,j) = 2*v_RIGHT - vOld(Nx+1,j);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% STEP 1a. Solve the x momentum equation (Eqn. 1a of Handout 11) using 
% Gauss-Seidel. The result is stored in uStar. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First setup the coefficients
for i = 2:Nx   % These i and j are the indices of the u-cell
  for j = 2:Ny+1
    mdot_e =  rho * h * 0.5 * (uOld(i,j) + uOld(i+1,j));
    mdot_w = -rho * h * 0.5 * (uOld(i,j) + uOld(i-1,j));
    mdot_n =  rho * h * 0.5 * (vOld(i,j) + vOld(i+1,j));
    mdot_s = -rho * h * 0.5 * (vOld(i,j-1) + vOld(i+1,j-1));
% Center cell coefficient with applied under-relaxation
    au_C(i,j) = (0.5 * (mdot_e + mdot_w + mdot_n + mdot_s) + 4*visc) / alphaV;
    au_E(i,j) = - 0.5 * mdot_e + visc;
    au_W(i,j) = - 0.5 * mdot_w + visc;
    au_N(i,j) = - 0.5 * mdot_n + visc;
    au_S(i,j) = - 0.5 * mdot_s + visc;
    % East neighbor coefficient
    % West neighbor coefficient
    % North neighbor coefficient
    % South neighbor coefficient
    dU(i,j) = h / au_C(i,j);
  end
end

uStar = uOld;

for iter = 1:maxGSiter
    for i = 2:Nx   % These i and j are the indices of the u-cell
        for j = 2:Ny+1
            uStar(i,j) = (au_E(i,j) * uStar(i+1,j) + au_W(i,j) * uStar(i-1,j) ...
                + au_N(i,j) * uStar(i,j+1) + au_S(i,j) * uStar(i,j-1) ...
                - h * (pOld(i+1,j) - pOld(i,j)) ...
                + (1-alphaV) * au_C(i,j) * uOld(i,j) ) / au_C(i,j);
        end
    end
end

for i = 2:Nx+1   % These i and j are the indices of the v-cell
  for j = 2:Ny
    mdot_e =  rho * h * 0.5 * (uOld(i,j) + uOld(i,j+1));
    mdot_w = -rho * h * 0.5 * (uOld(i-1,j) + uOld(i-1,j+1));
    mdot_n =  rho * h * 0.5 * (vOld(i,j) + vOld(i,j+1));
    mdot_s = -rho * h * 0.5 * (vOld(i,j) + vOld(i,j-1));
% Center cell coefficient with applied under-relaxation
    av_C(i,j) = (0.5 * (mdot_e + mdot_w + mdot_n + mdot_s) + 4*visc) / alphaV;
    
    av_E(i,j) = - 0.5 * mdot_e + visc;
    av_W(i,j) = - 0.5 * mdot_w + visc;
    av_N(i,j) = - 0.5 * mdot_n + visc;
    av_S(i,j) = - 0.5 * mdot_s + visc;
    
    dV(i,j) = h / av_C(i,j);
  end
end

vStar = vOld;

for iter = 1:maxGSiter
    for i = 2:Nx+1 % These i and j are the indices of the v-cell
        for j = 2:Ny
            vStar(i,j) = (av_E(i,j) * vStar(i+1,j) + av_W(i,j) * vStar(i-1,j) ...
                + av_N(i,j) * vStar(i,j+1) + av_S(i,j) * vStar(i,j-1) ...
                - h * (pOld(i,j+1) - pOld(i,j)) ...
                + (1-alphaV) * av_C(i,j) * vOld(i,j) ) / av_C(i,j);
        end
    end
end

% First setup the coefficients
for i = 2:Nx+1   % These i and j are the indices of the p-cell
  for j = 2:Ny+1
    aPC_C(i,j) = dU(i,j) + dU(i-1,j) + dV(i,j) + dV(i,j-1);  % Center cell coefficient
    aPC_E(i,j) = dU(i,j);
    aPC_W(i,j) = dU(i-1,j);
    aPC_N(i,j) = dV(i,j);
    aPC_S(i,j) = dV(i,j-1);
    bPC(i,j) = -uStar(i,j) + uStar(i-1,j) - vStar(i,j) + vStar(i,j-1);  % RHS value
  end
end

i = 2;
j = 2;
aPC_C(i,j) = 1.0;  % Center cell coefficient
aPC_E(i,j) = 0.0;  % East neighbor coefficient
aPC_W(i,j) = 0.0;  % West neighbor coefficient
aPC_N(i,j) = 0.0;  % North neighbor coefficient
aPC_S(i,j) = 0.0;  % South neighbor coefficient
bPC(i,j)   = 0.0;  % RHS value

% Gauss-Seidel loop
pPrime(:,:) = 0;   % Initialize pPrime's to zero.
for iter = 1:maxGSiter
  for i = 2:Nx+1
    for j = 2:Ny+1
      pPrime(i,j) = (aPC_E(i,j) * pPrime(i+1,j) + aPC_W(i,j) * pPrime(i-1,j) ...
                   + aPC_N(i,j) * pPrime(i,j+1) + aPC_S(i,j) * pPrime(i,j-1) ...
                   + bPC(i,j) ) / aPC_C(i,j);
    end
  end
end

for i = 2:Nx
  for j = 2:Ny+1
    uPrime(i,j) = -(pPrime(i+1,j) - pPrime(i,j)) * dU(i,j);
  end
end

for i = 2:Nx+1
  for j = 2:Ny
     vPrime(i,j) = -(pPrime(i,j+1) - pPrime(i,j)) * dV(i,j);
  end
end

for i = 2:Nx+1
  for j = 2:Ny+1
    p(i,j) = pOld(i,j) + alphaP * pPrime(i,j);
  end
end

for i = 2:Nx
  for j = 2:Ny+1
     u(i,j) = uStar(i,j) + uPrime(i,j);
  end
end

for i = 2:Nx+1
  for j = 2:Ny
     v(i,j) = vStar(i,j) + vPrime(i,j);
  end
end

% uCenter = u(Nx/2+1, Ny/2+1);
% vCenter = v(Nx/2+1, Ny/2+1);

uOld = u;
vOld = v;
pOld = p;

end