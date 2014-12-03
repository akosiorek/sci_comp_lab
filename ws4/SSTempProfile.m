function Tinf = SSTempProfile(~, ~, ~, T, Lx, Ly)

    Nx = (size(T,2)-2); % Total pts in x: Nx + 2
    Ny = (size(T,1)-2); % Total pts in y: Ny + 2
    dx = Lx/(size(T,2)-1);  % x-space step
    dy = Ly/(size(T,1)-1);  % y-space step
    b = zeros(Nx*Ny, 1);
    A = CreateMatrixA(Nx, Ny, dx, dy);
    Temp = A\b;
    Tinf = zeros(Nx+2, Ny+2);
    for i = 1: Nx
        for j = 1: Ny
            Tinf(j+1, i+1) = Temp((i-1)*Nx+j);
        end
    end
end
