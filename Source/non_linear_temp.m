function nltt = non_linear_temp(params,velocity,T)
    %NON_LINEAR_TEMP Computes the non-linear term duT/dx + dvT/dy present in the 
    % Energy equation.
    %
    % nltt = non_linear_temp(params,velocity,T)
    %
    % Variable lookup:
    %
    % velocity: Velocity field (EdgeData).
    %
    % T: Temperature field (CellData).
    %
    % domain: domain parameters.
    %
    % nltt: computed non-linear terms.
    %
    % Created by Jay Mehta (18 July 2019)
    
    %% Four terms to evaluate for d(uu) and d(vu) for X-direction
    
    Nx = T.size(1);
    Ny = T.size(2);
    dx = params.dx;
    dy = params.dx;
    
    nltt = CellData(Nx,Ny);
    
%     w = min(1.2 * params.dt * max(max(max(abs(velocity.x))),...
%         max(max(abs(velocity.y)))),1);
    w = 1;

    T1 = interpol(T,EdgeData(Nx,Ny),1);
    T1_upwind = upwinding(T,EdgeData(Nx,Ny),1);
    uT = EdgeData(Nx,Ny);
    uT.x = velocity.x .* (T1.x - w * T1_upwind.x);
    duTdx = div(uT,CellData(Nx,Ny),1);
    duTdx.x = duTdx.x/dx;
    
    T2 = interpol(T,EdgeData(Nx,Ny),2);
    T2_upwind = upwinding(T,EdgeData(Nx,Ny),2);
    vT = EdgeData(Nx,Ny);
    vT.y = velocity.y .* (T2.y - w * T2_upwind.y);
    dvTdy = div(vT,CellData(Nx,Ny),2);
    dvTdy.x = dvTdy.x/dy;
    
    nltt.x = duTdx.x + dvTdy.x;
    nltt.x(:,1) = 0;
    nltt.x(:,end) = 0;
    nltt.x(1,:) = 0;
    nltt.x(end,:) = 0;
end