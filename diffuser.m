function velocity = diffuser(params,bc,rhs,velocity)
    
    Nx = velocity.size(1);
    Ny = velocity.size(2);
    
    dx = params.dx;
    nu = params.nu;
    
    velocity_bc = EdgeData(Nx,Ny);
    velocity_bc = apply_bc(bc,velocity_bc);
    
    %% For X-direction (U-velocity)
    
    A = zeros(Nx-1,Nx-1);
    for j = 2:Ny+1
        for i = 2:Nx
            A(i-1,i-1) = 2*nu/dx^2; % Center
            if i == 2
                A(i,i-1) = -nu/dx^2; % East
            elseif i == Nx
                A(i-2,i-1) = -nu/dx^2; % West
            else
                A(i,i-1) = -nu/dx^2 ; % East
                A(i-2,i-1) = -nu/dx^2; % West
            end
        end
    
    % Now constructing the AX=B problem.

        rhs.x(2,j) = rhs.x(2,j) + (nu/dx^2) * velocity_bc.x(1,j);
        rhs.x(Nx,j) = rhs.x(Nx,j) + (nu/dx^2) * velocity_bc.x(Nx+1,j);
        B = rhs.x(2:Nx,j);
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
        velocity.x(2:Nx,j) = trisolve(a,b,c,B,'reg');
    end

%% For X-direction (V-velocity)

% This is the off-direction one. So the diagonal terms will have to be 
% modified.
    
    A = zeros(Nx,Nx);
    for j = 2:Ny
        for i = 2:Nx+1
            A(i-1,i-1) = 2*nu/dx^2; % F(East - West) + D(East-West)
            if i == 2
                A(i,i-1) = -nu/dx^2; % east face
                A(i-1,i-1) = A(i-1,i-1) + nu/dx^2; % Modified term
            elseif i == Nx+1
                A(i-2,i-1) = -nu/dx^2; % West
                A(i-1,i-1) = A(i-1,i-1) + nu/dx^2; % Modified term
            else
                A(i,i-1) = -nu/dx^2; % East
                A(i-2,i-1) = -nu/dx^2; % West
            end
        end
    
    % Now constructing the AX=B problem.

        rhs.y(2,j) = rhs.y(2,j) + (nu/dx^2) * velocity_bc.y(1,j);
        rhs.y(Nx+1,j) = rhs.y(Nx+1,j) + (nu/dx^2) * velocity_bc.y(Nx+2,j);
        B = rhs.y(2:Nx+1,j);
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
            velocity.y(2:Nx+1,j) = trisolve(a,b,c,B,'reg');
    end

    %% Round 2 (X & Y are interchanged and transposed)
    
    rhs.x = velocity.y';
    rhs.y = velocity.x';
    temp = velocity;
    velocity.x = temp.y';
    velocity.y = temp.x';
    temp = velocity_bc;
    velocity_bc.x = temp.y';
    velocity_bc.y = temp.x';
    
    %% For Y-direction (V-velocity)
    A = zeros(Ny-1,Ny-1);
    for j = 2:Nx+1
        for i = 2:Ny
            A(i-1,i-1) = 2*nu/dx^2; % Middle
            if i == 2
                A(i,i-1) = -nu/dx^2; % East
            elseif i == Ny
                A(i-2,i-1) = -nu/dx^2; % West
            else
                A(i,i-1) = -nu/dx^2; % East
                A(i-2,i-1) = -nu/dx^2; % West
            end
        end
        
        rhs.x(2,j) = rhs.x(2,j) + (nu/dx^2) * velocity_bc.x(1,j);
        rhs.x(Ny,j) = rhs.x(Ny,j) + (nu/dx^2) * velocity_bc.x(Ny+1,j);
        B = rhs.x(2:Ny,j);
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
        velocity.x(2:Ny,j) = trisolve(a,b,c,B,'reg');
    end

%% For Y-direction (U-velocity)

% This is the off-direction one. So the diagonal terms will have to be 
% modified
    
    A = zeros(Ny,Ny);
    for j = 2:Nx
        for i = 2:Ny+1
            A(i-1,i-1) = 2*nu/dx^2; % Middle
            if i == 2
                A(i,i-1) = -nu/dx^2; % East
                A(i-1,i-1) = A(i-1,i-1) + nu/dx^2; % Modified Term
            elseif i == Nx+1
                A(i-2,i-1) = -nu/dx^2;  % West
                A(i-1,i-1) = A(i-1,i-1) + nu/dx^2; % Modified Term
            else
                A(i,i-1) = -nu/dx^2; % East
                A(i-2,i-1) = -nu/dx^2; % West
            end
        end
        
        rhs.y(2,j) = rhs.y(2,j) + (nu/dx^2) * velocity_bc.y(1,j);
        rhs.y(Ny+1,j) = rhs.y(Ny+1) + (nu/dx^2) * velocity_bc.y(Ny+2,j);
        B = rhs.y(2:Ny+1,j);
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
            velocity.y(2:Ny+1,j) = trisolve(a,b,c,B,'reg');
    end
    
%% Undoing the interchange and transpose
    
    velocity_temp = velocity;
    velocity.x = velocity_temp.y';
    velocity.y = velocity_temp.x';
    
    velocity = apply_bc(bc,velocity);
    
end