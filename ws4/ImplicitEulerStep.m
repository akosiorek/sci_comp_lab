function [T, runTime] = ImplicitEulerStep(oldT, dt)

    
    s = size(oldT) - 2;    %   Assume size(oldT) = (Nx + 2, Ny + 2)
    Nx = s(1);
    Ny = s(2);
    T = double(oldT);  % create a solution matrix by copying the old solution and taking it as a starting point
    hx = 1/(1+Nx);  % step sizes
    hy = 1/(1+Ny);
    
    % define weights so as not to compute them on each iteration of the triple loop below
    hx2 = hx^-2;
    hy2 = hy^-2;
    invGlobalWeight = (1 + 2 * dt * (hx2 + hy2));
    hx2dt = hx2 * dt;
    hy2dt = hy2 * dt;
    globalWeight = 1 / invGlobalWeight;
    xWeight = globalWeight * hx2dt;
    yWeight = globalWeight * hy2dt;
    
    res = 1;  % Define residual parameter
    tic
    %iter = 0;
    while res > 10^-6    % Check for convergence
        res = 0;    % Set residual at beginning of iteration to zero
        for j = 2 : Ny + 1
            for i = 2 : Nx + 1;
                T(j, i) = globalWeight * oldT(j, i) + xWeight * (T(j, i - 1) + T(j, i + 1))...
                    + yWeight * (T(j - 1, i) + T(j + 1, i));
            end
        end
        for j = 2 : Ny + 1
            for i = 2 : Nx + 1;                
                res = res + (oldT(j, i) + hx2dt * (T(j, i - 1) + T(j, i + 1))...
                    + hy2dt * (T(j - 1, i) + T(j + 1, i)) - invGlobalWeight * T(j, i)) ^ 2;
            end
        end
        res = sqrt(res/(Nx*Ny));  % Root-mean
       %iter = iter + 1;
    end   
    runTime = toc;  
end