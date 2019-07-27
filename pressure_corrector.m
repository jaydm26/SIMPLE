function pressure = pressure_corrector(params,velocity,velocity_star)
    
    Nx = velocity.size(1);
    Ny = velocity.size(2);
    
    dx = params.dx;
    nu = params.nu;
    
    rhs = CellData(Nx,Ny);
    pressure = CellData(Nx,Ny);
    
    velocity_x_n = interpol(velocity,NodeData(Nx,Ny),1);
    velocity_y_n = interpol(velocity,NodeData(Nx,Ny),2);
    velocity_x_c = interpol(velocity,CellData(Nx,Ny),1);
    velocity_y_c = interpol(velocity,CellData(Nx,Ny),2);
    
    div_velocity_star_x = div(velocity_star,CellData(Nx,Ny),1);
    div_velocity_star_y = div(velocity_star,CellData(Nx,Ny),2);
    
    for i = 1:Nx+1
        for j = 1:Ny+2
            rhs.x(i,j) = rhs.x(i,j) - div_velocity_star_x.x(i,j); % U^*(i,j) - U^*(i+1,j)
        end
    end
        
    for i = 1:Nx+2
        for j = 1:Ny+2
            rhs.x(i,j) = rhs.x(i,j) - div_velocity_star_y.x(i,j); % V^*(i,j) - V^*(i,j+1)
        end
    end
    
     %% Round 1
    
    div_velocity_x_c = div(velocity_x_c,EdgeData(Nx,Ny),1);
    div_velocity_y_n = div(velocity_y_n,EdgeData(Nx,Ny),2);
    A = zeros(Nx,Nx);
    for j = 2:Ny+1
        for i = 2:Nx+1
            % Need to calculate aN,aE,aS,aW for each (i,j).
            % In the first round, since we are solving in the X-direction,
            % we need to obttain aE and aW.
            aE = 0.5/dx * div_velocity_x_c.x(i,j) + 0.5/dx * div_velocity_y_n.x(i,j) + 4*nu/dx^2;
            aW = 0.5/dx * div_velocity_x_c.x(i-1,j) + 0.5/dx * div_velocity_y_n.x(i-1,j) + 4*nu/dx^2;
            A(i-1,i-1) = 1/dx * (1/aE + 1/aW);
            if i == 2
                A(i,i-1) = -1/dx * 1/aE; % East
            elseif i == Nx+1
                A(i-2,i-1) = -1/dx * 1/aW; % West
            else
                A(i,i-1) = -1/dx * 1/aE; % East
                A(i-2,i-1) = -1/dx * 1/aW; % West
            end
        end
    
        % Now construct the AX = B problem
        
        rhs.x(2,j) = rhs.x(2,j) + 1/dx * 1/(0.5/dx * div_velocity_x_c.x(1,j) + 0.5/dx * div_velocity_y_n.x(1,j) + 4*nu/dx^2) * pressure.x(1,j);
        rhs.x(Nx+1,j) = rhs.x(Nx+1,j) + 1/dx * 1/(0.5/dx * div_velocity_x_c.x(Nx+1,j) + 0.5/dx * div_velocity_y_n.x(Nx+1,j) + 4*nu/dx^2) * pressure.x(Nx+2,j);
        B = rhs.x(2:Nx+1,j);
        a = zeros(length(A)-1,1);
        b = zeros(length(A),1);
        c = zeros(length(A)-1,1); 
        for i = 1:length(A)
            if i == 1
                b(i,1) = A(i,i);
                c(i,1) = A(i,i+1);
            elseif i == length(A)
                a(i-1,1) = A(i,i-1);
                b(i,1) = A(i,i);
            else
                a(i-1,1) = A(i,i-1);
                b(i,1) = A(i,i);
                c(i,1) = A(i,i+1);
            end
        end
        pressure.x(2:Nx+1,j) = trisolve(a,b,c,B,'reg');
    end
    
    
    %% Round 2
    
    rhs.x = pressure.x';
    pressure.x = pressure.x';
    velocity_y_c.x = velocity_y_c.x';
    velocity_x_n.x = velocity_x_n.x';
    div_velocity_y_c = div(velocity_y_c,EdgeData(Nx,Ny),1);
    div_velocity_x_n = div(velocity_x_n,EdgeData(Nx,Ny),2);
    
    A = zeros(Ny,Ny);
    for j = 2:Nx+1
        for i = 2:Ny+1
            % Need to calculate aN,aE,aS,aW for each (i,j).
            % In the second round, since we are solving in the Y-direction,
            % we need to obttain aN and aS.
            aN = 0.5/dx * div_velocity_y_c.x(i,j) + 0.5/dx * div_velocity_x_n.x(i,j) + 4*nu/dx^2;
            aS = 0.5/dx * div_velocity_y_c.x(i-1,j) + 0.5/dx * div_velocity_x_n.x(i-1,j) + 4*nu/dx^2;
            A(i-1,i-1) = 1/dx * (1/aN + 1/aS);
            if i == 2
                A(i,i-1) = -1/dx * 1/aN; % North
            elseif i == Ny+1
                A(i-2,i-1) = -1/dx * 1/aS; % South
            else
                A(i,i-1) = -1/dx * 1/aN; % North
                A(i-2,i-1) = -1/dx * 1/aS; % South
            end
        end
        
        rhs.x(2,j) = rhs.x(2,j) + 1/dx * 1/(0.5/dx * div_velocity_y_c.x(1,j) + 0.5/dx * div_velocity_x_n.x(1,j) + 4*nu/dx^2) * pressure.x(1,j);
        rhs.x(Ny+1,j) = rhs.x(Ny+1,j) + 1/dx * 1/(0.5/dx * div_velocity_y_c.x(Ny+1,j) + 0.5/dx * div_velocity_x_n.x(Ny+1,j) + 4*nu/dx^2) * pressure.x(Ny+2,j);
        B = rhs.x(2:Ny+1,j);
        a = zeros(length(A)-1,1);
        b = zeros(length(A),1);
        c = zeros(length(A)-1,1);
        for i = 1:length(A)
            if i == 1
                b(i,1) = A(i,i);
                c(i,1) = A(i,i+1);
            elseif i == length(A)
                a(i-1,1) = A(i,i-1);
                b(i,1) = A(i,i);
            else
                a(i-1,1) = A(i,i-1);
                b(i,1) = A(i,i);
                c(i,1) = A(i,i+1);
            end
        end
        pressure.x(2:Ny+1,j) = trisolve(a,b,c,B,'reg');
    end

    %% Undoing the interchange and transpose
    
    pressure.x = pressure.x';
    
    pressure.x(1,2:Ny+1) = pressure.x(2,2:Ny+1);
    pressure.x(end,2:Ny+1) = pressure.x(end-1,2:Ny+1);
    pressure.x(2:Nx+1,1) = pressure.x(2:Nx+1,2);
    pressure.x(2:Nx+1,end) = pressure.x(2:Nx+1,end-1);
    
    pressure.x(1,1) = 0;
    pressure.x(1,end) = 0;
    pressure.x(end,1) = 0;
    pressure.x(end,end) = 0;
    
end
