function Tnxt = ImpEulTemporalSolver(delta, ts, te, T, Lx, Ly)
    
    Nx = (size(T,2)-2); % Total pts in x: Nx + 2
    Ny = (size(T,1)-2); % Total pts in y: Ny + 2
    dx = Lx/(size(T,2)-1); % x-space step
    dy = Ly/(size(T,1)-1); % y-space step
    
    c = 1 + 2*delta*((dx)^-2 + (dy)^-2);    % Coefficient for center term
    h = delta/(dx)^2;   % Coefficient for left, right terms ('horizontal')
    v = delta/(dy)^2;   % Coefficient for top, bottom terms ('vertical')
    
    time = ts;
    Tnxt = T;
    while time < te
        T = Tnxt;   % Create matrix for the next timestep
        res = 1;
        while res > 1e-6
            for i = 2: Nx+1
                for j = 2: Ny+1
                    Tnxt(j,i) = (1/c)*(T(j,i) + h*(Tnxt(j-1,i) + Tnxt(j+1,i)) + v*(Tnxt(j,i-1) + Tnxt(j,i+1)));
                end
            end
            res = 0;
            for i = 2: Nx+1
                for j = 2: Ny+1;        
                    res = res + (c*Tnxt(j,i) - h*(Tnxt(j-1,i) + Tnxt(j+1,i)) - v*(Tnxt(j,i-1) + Tnxt(j,i+1)) - T(j,i))^2; % Calculate Residual
                end
            end
            res = sqrt(res/(Nx*Ny))  % Root-mean
        end
        time = time + delta;
    end
end
