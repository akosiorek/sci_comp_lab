%Scientific Computing Lab - Computational Science and Engineering
%Worksheet 2 Submission (11.11.2014)

%Group 8:
%Joshi, Saumitra Vinay
%Kosiorek, Adam
%Perez de Alba Ortiz, Alberto

function [] = worksheet3()
     close all
     clear
     clc

     Nx = [7, 15, 31, 63, 127]; Ny=Nx; %Arrays of Nx,Ny values and names
     N= {'7' , '15' , '31', '63', '127'};
     
     solvers = {@FullMatrixSolver, @SparseMatrixSolver,@GaussSeidelSolverS}; %Arrays of solver functions and names
     solversNames = {'direct solution with full matrix', 'direct solution with sparse matrix', 'iterative solution with Gauss-Seidel'};
     
     solutions = cell(length(solvers),length(Nx)); %Arrays for storing solutions, runtimes, storage requirements and errors
     runTimes = zeros(length(solvers),length(Nx));
     storageReqs = zeros(length(solvers),length(Nx));
     errors = zeros(length(solvers),length(Nx));
    
    
    for i=1:length(Nx); %Loop for calculating solutions
            
        A=CreateMatrixA(Nx(i),Ny(i)); %Create matrix A and vector b corresponding to PDE
        b=CreateVectorb(Nx(i),Ny(i));
        GSS=0; %Indicate if Gauss-Seidel Solver is going to be used (it demands changes on b: including boundary values)
        
        for j=1:length(solvers)
            
            if (i==length(Nx) && j~=length(solvers) ) 
                %Do not calculate last Nx,Ny solution for direct solvers
            elseif (j==length(solvers))
                GSS=1; %Using Gauss-Seidel Solver
                [solutions{j,i}, runTimes(j,i), storageReqs(j,i)] = solvers{j}(A,b,Nx(i),Ny(i)); %Calling corresponding solver function
                errors(j,i)=CalculateError(Nx(i),Ny(i),solutions{j,i},GSS); %Calculating corresponding error
            elseif (j~=length(solvers))
                [solutions{j,i}, runTimes(j,i), storageReqs(j,i)] = solvers{j}(A,b,Nx(i),Ny(i)); %Calling corresponding solver function
                errors(j,i)=CalculateError(Nx(i),Ny(i),solutions{j,i},GSS); %Calculating corresponding error
            end
                
        end
  
    end

        
    GSS=0; %Indicate if Gauss-Seidel Solver is going to be used (it demands changes on A and b)
    for i=1:length(solvers) %Loop for plots and result tables.
        index=1; %Index for subplots
        
        for j=1:length(Nx)-1
        
            if (i==length(solvers))
                GSS=1; %Using Gauss-Seidel Solver
            end
            
            figure (i)
           
            PlotResults(solutions{i,j}, Nx(j), Ny(j),index,N{j},GSS) %Plot results
            suptitle(solversNames{i}); %Subplots title
            
            index=index+2; %Update subplot index
            
            
        end
        
            %Display result table
            Nt= {'     7\t' , '      15\t' , '31\t', '  63\t', '   127\t'};
            fprintf(strcat('\t\t', solversNames{i}, '\n\n'));
            fprintf('Nx, Ny\t'); fprintf(strcat('\t',Nt{1:end-1},'\n\n'));
            fprintf('runtime \t'); format short; disp(runTimes(i,1:end-1));
            fprintf('storage\t '); format short; disp(storageReqs(i,1:end-1));
    end
    
    %Calculate error reduction
    errorRed=zeros(length(Nx));
    for j = 2: length(errors)    
        errorRed(end,j) = (errors(end, j - 1) / errors(end, j));
    end
    
    %Display error reduction table for GSS
    fprintf(strcat('\t\t', solversNames{end}, '\n\n'))
    fprintf('Nx, Ny\t'); fprintf(strcat('\t',Nt{:},'\n\n'));
    fprintf('error\t\t'); format short; disp(errors(end,:));
    fprintf('error red.\t'); format short; disp(errorRed(end,:)); 
    
end

function [A] = CreateMatrixA(Nx, Ny)

    A = eye(Nx * Ny) * - 2 * ( (Nx + 1)^2 + (Ny + 1)^2); %Set A size 

    for iy=0:Ny-1 %Loops for calculating each value 
        for ix=1:Nx  clc

            if (ix~=1) %Conditions based of PDE form
                A( (Nx * iy) + ix , (Nx * iy) + ix - 1) = (Nx + 1)^2;
            end
            if (ix~=Nx )
                A( (Nx * iy) + ix , (Nx * iy) + ix + 1) = (Nx + 1)^2;
            end

            if ( iy>=1 )
                A( (Nx * iy) + ix , (Nx * iy) + ix - Nx)=(Ny + 1)^2;
            end
            if (iy~=Nx-1)
                A( (Nx * iy) + ix , (Nx * iy) + ix + Nx) = (Ny + 1)^2;
            end

        end
    end
end

function [b] = CreateVectorb(Nx,Ny)

    b = zeros(Nx * Ny,1); %Set b size
    index=1; %set increasing index value
    for iy=1:Ny %Loop for filling values of each element
        for ix=1:Nx
            b(index)= -2 * pi^2 * sin(pi *  ix/(Nx + 1) ) * sin (pi *  iy/(Ny + 1) ); %Exact solution
            index=index+1;
        end
    end

end

function [T, runTime, storageReq] = FullMatrixSolver(A,b,Nx,Ny)

    tic 

    T = A\b; %Solve

    runTime = toc; %Calculate runtime

    storageReq = numel(A) + numel(b) + numel(T); %Calculate storage requirements

end

function [T, runTime, storageReq] = SparseMatrixSolver(A,b,Nx,Ny)

    S = sparse(A); %transform A to sparse matri

    tic

    T = S\b; %Solve

    runTime = toc; %Calculate runtime

    storageReq = nnz(S) + numel(b) + numel(T); %Calculate storage requirements

end

function [] = PlotResults(t, Nx, Ny, index,N, GSS)

  if (GSS~=1) %if plot is for GSS, then Vector2Matrix is not needed
      T=Vector2Matrix(t, Nx,Nx);
  else
      T=t;
  end
  x = 0:1/(Nx+1):1; %Set x and y
  y = 0:1/(Ny+1):1;
  subplot(4, 2, index); %Assign subplot position
  index=index+1;
  mesh(x, y, T); %Surface plot
  title(strcat('Surface plot for Nx=Ny= ', N)); %Add title, labels...
  xlabel('x');
  ylabel('y');
  zlabel('T');
  subplot(4, 2, index);
  contour(x, y, T); %Contour plot
  title(strcat('Countour plot for Nx=Ny= ', N));
  xlabel('x');
  ylabel('y');
  zlabel('T');
  
end

function [T] = Vector2Matrix(t, Nx, Ny)

    T=zeros(Nx + 2, Ny + 2); %Set matrix T size
    index=1;
    for iy=2:Ny+1 %Loops for assigning vector elements to matrix elements.
        for ix=2:Nx+1
            T(iy, ix)= t(index);
            index=index+1;
        end
    end
end

function [e] = CalculateError(Nx,Ny,T,GSS)

if (GSS~=1)
    T=Vector2Matrix(T,Nx,Ny); %if plot is for GSS, then Vector2Matrix is not needed
end
hx = 1/(1+Nx); hy = 1/(1+Ny); %calculate hx,y
e=0;

    for j = 2: Ny+1 %Loops for calculating each elements difference
        for i = 2: Nx+1
            e = e + (T(j,i) - sin(pi*(i-1)*hx)*sin(pi*(j-1)*hy))^2;
             
        end
    end
    e = sqrt(e/(Nx*Ny)); %Calculate error
end

function [T, runTime, storageReq] = GaussSeidelSolverS(A, b, Nx, Ny)

bBorders=zeros((Nx+2)*(Ny+2),1);

for i=1:Nx
    bBorders((Nx+2)*i+2:(Nx+2)*i+2+Nx-1)=b(1+Nx*(i-1):i*Nx);  
end

% b has length (Nx+2 * Ny+2)
% Solves matrix T (Nx+2 x Ny+2) using Gauss-Seidel method
    
    tic
    T = zeros(Ny+2, Nx+2);  % Create Solution Matrix, automatically set boundaries to zero, initial guess to zero
    hx = 1/(1+Nx); hy = 1/(1+Ny);   % Discretization levels
    res = 1;  % Define residual parameter
    while res > 10^-4    % Check for convergence
        res = 0;    % Set residual at beginning of iteration to zero
        for j = 2: Ny+1
            for i = 2: Nx+1;
                T(j, i) = 0.5*(1/(1/hx^2 + 1/hy^2)) * ((T(j,i-1) + T(j,i+1))/(hx^2)...
                    + (T(j-1,i) + T(j+1,i))/(hy^2) - bBorders((Nx+2)*(j-1)+i));    % solve for each point                
            end
        end
        for j = 2: Ny+1
            for i = 2: Nx+1;        
                res = res + (bBorders((Nx+2)*(j-1)+i) - (T(j,i-1) + T(j,i+1))/(hx^2)...
                    - (T(j-1,i) + T(j+1,i))/(hy^2) + 2*(1/hx^2 + 1/hy^2)*T(j,i))^2; % Calculate residual
            end
        end
        res = sqrt(res/(Nx*Ny));  % Root-mean
    end
   
    runTime = toc;
    storageReq= numel(T) + numel(b);
    
end
