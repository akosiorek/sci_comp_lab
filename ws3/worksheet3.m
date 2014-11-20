function worksheet3()
   close, clear, clc
   
    Nx = 27;
    Ny = 27;
    
    hx = 1 / (Nx + 1);
    hy = 1 / (Ny + 1);
    
    gaussSeidlTolerance = 10^-4;
    
    
    exact_solution = sin(pi * (0:hy:1))' * sin(pi * (0:hx:1));
    
    
    d2T = Rhs(hy:hy:1-hy, hx:hx:1-hx);
    direct_solution = DirectSolver(Ny, Nx, hy, hx, d2T, exact_solution);
    iterative_solution = GaussSeidl(Ny, Nx, hy, hx, d2T, gaussSeidlTolerance, exact_solution);
    
    subplot(1, 3, 1)
    contourf(exact_solution);    
    title('exact');
    
    subplot(1, 3, 2)
    contourf(direct_solution);    
    title('direct');
    
    subplot(1, 3, 3)
    contourf(iterative_solution);
    title('iterative');
    
end



function coeffs = CoefficientMatrix(Nx, Ny, hx, hy)
   % CoefficientMatrix for the heat equation linear system
   
    [a, b, c] = Coeffs(hx, hy);
    template = [a b zeros(1, Nx - 2) c];
    coeffs = toeplitz([template zeros(1, Nx * Ny - length(template))]);
end

function [a, b, c] = Coeffs(hx, hy) 

    c = hx ^ -2;
    b = hy ^ -2;
    a = -2 * (b + c);    
end

function rhs = Rhs(y, x) 
    rhs = - 2 * pi * pi * sin(pi * y)' * sin(pi * x);
end 

function solution = DirectSolver(Nx, Ny, hx, hy, rhs, exact_solution)
    
    dTCoeffs = CoefficientMatrix(Nx, Ny, hx, hy);   
    d2T = reshape(rhs, [numel(rhs), 1]);
    T = dTCoeffs \ d2T;     
    solution = padarray(reshape(T, Ny, Nx), [1 1]);
    disp({'Direct norm = ',sqrt(mean(mean((solution - exact_solution).^2)))});
end

function solution = GaussSeidl(Ny, Nx, hy, hx, rhs, tol, exact_solution) 
    
    [a, b, c] = Coeffs(hx, hy);
    
    b = b / -a;
    c = c / -a; 
    rhs = rhs / a;
    x = zeros(Ny + 2, Nx + 2);
    
    residual = realmax;
    iter = 0;
    N = (Ny + 2) * (Nx + 2);
    while residual > tol
        for i = 2 : Ny + 1
            for j = 2 : Nx + 1
                x(i, j) = b * (x(i, j-1) + x(i, j+1)) + c * (x(i-1, j) + x(i+1, j)) + rhs(i-1, j-1);
            end
        end        
            
        residual = sqrt(mean(mean((x - exact_solution).^2)));
        iter = iter + 1;
    end
    
    disp({'GaussSeidl norm = ', residual});
    solution = x;
end















