function Tnxt = ExpEulTemporalSolver(delta, ts, te, T, Lx, Ly)
    
    Nx = (size(T,2)-2); % Total pts in x: Nx + 2
    Ny = (size(T,1)-2); % Total pts in y: Ny + 2
    dx = Lx/(Nx+1); % x-space step
    dy = Ly/(Ny+1); % y-space step
    
    c = 1 - 2*delta*((dx)^-2 + (dy)^-2);    % Coefficient for center term
    h = delta/(dx)^2;   % Coefficient for left, right terms ('horizontal')
    v = delta/(dy)^2;   % Coefficient for top, bottom terms ('vertical')
                
    time = ts;
    Tnxt = T;
    
    while time < te
        T = Tnxt;   % Create matrix for the next timestep
        for j = 2 : Ny+1     
            for i = 2 : Nx+1        
                Tnxt(j,i) = c*T(j,i) + h*(T(j-1,i) + T(j+1,i)) + v*(T(j,i-1) + T(j,i+1));
            end
        end       
        time = time + delta;
    end
end
