function worksheet3()
   
    Nx = 10;
    Ny = 10;
    
    coeffs  = CoefficientMatrix(Nx, Ny)
    
    
end



function coeffs = CoefficientMatrix(Nx, Ny)
   % CoefficientMatrix for the heat equation linear system
   
    template = [-4 1 0 1];
    coeffs = toeplitz([template zeros(1, Nx * Ny - 4)]);
end