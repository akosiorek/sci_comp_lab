function [T, time, error] = GaussSeidelSolver(b, Nx, Ny)
% b has length (Nx+2 * Ny+2)
% Solves matrix T (Nx+2 x Ny+2) using Gauss-Seidel method

    tic
    T = zeros(Ny+2, Nx+2);  % Create Solution Matrix, automatically set boundaries to zero, initial guess to zero
    hx = 1/(1+Nx); hy = 1/(1+Ny);   % Discretization levels
    res = 1;  % Define residual parameter
    while res > 1e-4    % Check for convergence
        res = 0;    % Set residual at beginning of iteration to zero
        for j = 2: Ny+1
            for i = 2: Nx+1;
                T(j, i) = 0.5*(1/(1/hx^2 + 1/hy^2)) * ((T(j,i-1) + T(j,i+1))/(hx^2)...
                    + (T(j-1,i) + T(j+1,i))/(hy^2) - b((Nx+2)*(j-1)+i));    % solve for each point                
            end
        end
        for j = 2: Ny+1
            for i = 2: Nx+1;        
                res = res + (b((Nx+2)*(j-1)+i) - (T(j,i-1) + T(j,i+1))/(hx^2)...
                    - (T(j-1,i) + T(j+1,i))/(hy^2) + 2*(1/hx^2 + 1/hy^2)*T(j,i))^2; % Calculate residual
            end
        end
        res = sqrt(res/(Nx*Ny));  % Root-mean
    end
    error = 0;    % res is now used to store error
    for j = 2: Ny+1
        for i = 2: Nx+1
            error = error + (T(j,i) - sin(pi*(i-1)*hx)*sin(pi*(j-1)*hy))^2;
        end
    end
    error = sqrt(error/(Nx*Ny));
    time = toc;
end
