function velocity = velocity_corrector(params,bc,velocity,velocity_star,pressure_prime)
    
    Nx = velocity.size(1);
    Ny = velocity.size(2);
    
    dx = params.dx;
    nu = params.nu;
    
    velocity_x_n = interpol(velocity,NodeData(Nx,Ny),1);
    velocity_y_n = interpol(velocity,NodeData(Nx,Ny),2);
    velocity_x_c = interpol(velocity,CellData(Nx,Ny),1);
    velocity_y_c = interpol(velocity,CellData(Nx,Ny),2);
    
    div_velocity_x_c = div(velocity_x_c,EdgeData(Nx,Ny),1);
    div_velocity_y_n = div(velocity_y_n,EdgeData(Nx,Ny),2);
    div_velocity_y_c = div(velocity_y_c,EdgeData(Nx,Ny),2);
    div_velocity_x_n = div(velocity_x_n,EdgeData(Nx,Ny),1);
    grad_pressure_prime = grad(pressure_prime);
    
    for i = 2:Nx
        for j = 2:Ny+1
            a = 0.5/dx * div_velocity_x_c.x(i,j) + 0.5/dx * div_velocity_y_n.x(i,j) + 4 * nu/dx^2;
            velocity.x(i,j) = velocity_star.x(i,j) - 1/dx * 1/a * grad_pressure_prime.x(i,j);
        end
    end
    
    for i = 2:Nx+1
        for j = 2:Ny
            a = 0.5/dx * div_velocity_y_c.y(i,j) + 0.5/dx * div_velocity_x_n.y(i,j) + 4 * nu/dx^2;
            velocity.y(i,j) = velocity_star.y(i,j) - 1/dx * 1/a * grad_pressure_prime.y(i,j);
        end
    end
    
    velocity = apply_bc(bc,velocity);
end
