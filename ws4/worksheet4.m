function [] = worksheet4()

    close all
    clear
    clc
    Lx=1; Ly=Lx;
    Nx = [3, 7, 15, 31]; Ny = Nx;   % Arrays of Nx, Ny values
    N = {'3', '7', '15', '31'}; % N: string version
    deltas = 2.^(-4: -1: -10);  % Array of timesteps
    dels = {'1/64', '1/128', '1/256', '1/512', '1/1024', '1/2048', '1/4096'};   % deltas: string version
    
    dispTimes = [1/8, 2/8, 3/8, 4/8];   % Create array of display times
    dT = {'1/8', '2/8', '3/8', '4/8'};   % Display times: string version
    
    figH = 1: length(dispTimes);    % Create figure handles
    for i = 1: length(figH)
        figure(figH(i));    % Create figures
        set(figH(i), 'Visible', 'off'); % Hide them temporarily
        set(figH(i), 'Position', [25, 25, 1300, 800]);  % Set size
        set(gcf, 'color', 'white');
    end
    
    %Create matrix to mark stable cases
    stabilityExp=ones(length(Nx),length(deltas));
    stabilityImp=ones(length(Nx),1);
    
    fprintf('Solving Explicit Euler scheme ... \n\n');
    for i = 1: length(Nx)
        fprintf(strcat('Solving for Nx=Ny=', N{i}, ' ... \n\n' ));
        for j = 1: length(deltas)
            T = [zeros(1,Nx(i)+2); zeros(Nx(i),1),ones(Nx(i),Nx(i)), zeros(Nx(i),1); zeros(1,Nx(i)+2)]; %Initialize T matrix 
            for k=1:length(dispTimes)
                T = ExpEulTemporalSolver(deltas(j), dispTimes(k)-1/8, dispTimes(k), T, Lx, Ly);  % Call Explicit Euler Temporal Solver and update T
                x = 0: 1/(Nx(i)+1): 1;   % Set x and y
                y = 0: 1/(Ny(i)+1): 1;
                
                figure(figH(k)); %Set new figure
                rows=length(Nx); %Set subplot rows and columns
                columns=length(deltas);
                handles=zeros(rows,columns);
                handles(i,j) = subplot(rows,columns,j+(i-1)*columns);
                mesh(x, y, T); %Surface plot
                pbaspect([1 1 1]);
                xlabel('x');
                ylabel('y');
                zlabel('T');           
                
                %Subplot column labels
                if i == 1
                    text(0.25, 1.25, strcat('dt=', dels{j}), 'Units', 'normalized', 'FontSize', 12);
                end
                %Subplot row labels
                if j == 1
                    text(-1.1, 0.5, strcat('Nx=Ny=', N{i}), 'Units', 'normalized', 'FontSize', 12);
                end
                
                set(figH(k), 'Visible', 'off'); % Hide figure temporarily
                
                %Indicate if case is stable 
                %*considering that unstability causes oscillations that make T < 0
                if (any(T(:) < 0))
                stabilityExp(i,j)=0;
                end
              
            end
        end
    end
    
    fprintf('Solving Implicit Euler scheme ... \n\n');
    for i = 1: length(Nx)
        fprintf(strcat('Solving for Nx=Ny=', N{i}, ' ... \n\n' ));
        j = 1;
            T = [zeros(1,Nx(i)+2); zeros(Nx(i),1),ones(Nx(i),Nx(i)), zeros(Nx(i),1); zeros(1,Nx(i)+2)]; %Initialize T matrix 
            for k=1:length(dispTimes)
                T = ImpEulTemporalSolver(deltas(j), dispTimes(k)-1/8, dispTimes(k), T, Lx, Ly);  % Call Explicit Euler Temporal Solver and update T
                x = 0: 1/(Nx(i)+1): 1;   % Set x and y
                y = 0: 1/(Ny(i)+1): 1;
                
                figure(4+k); %Set new figure
                rows=length(Nx); %Set subplot rows and columns
                columns=1;
                handles=zeros(rows,columns);
                handles(i,j) = subplot(rows,columns,j+(i-1)*columns);
                mesh(x, y, T); %Surface plot
                pbaspect([1 1 1]);
                xlabel('x');
                ylabel('y');
                zlabel('T');           
                
                %Subplot column labels
                if i == 1
                    text(0.25, 1.25, strcat('dt=', dels{j}), 'Units', 'normalized', 'FontSize', 12);
                end
                %Subplot row labels
                if j == 1
                    text(-1.75, 0.5, strcat('Nx=Ny=', N{i}), 'Units', 'normalized', 'FontSize', 12);
                end
                
                set(figH(k), 'Visible', 'off'); % Hide figure temporarily
                
                %Indicate if case is stable 
                %*considering that unstability causes oscillations that make T < 0
                if (any(T(:) < 0))
                stabilityImp(i,j)=0;
                end
        end
    end
    
    for i = 1: length(figH)
        figure(figH(i));    % Activate figure i
        suptitle(strcat('State of Solutions of Explicit Euler scheme at t = ', dT(i)));  % Set title
    end
    
    for i=1:4
        figure(4+i);
        suptitle(strcat('State of Solutions of Implicit Euler scheme at t = ', dT(i)));
    end
    
    %Display stability matrix
    dtlabels='dt=1/64 dt=1/128 dt=1/256 dt=1/512 dt=1/1024 dt=1/2048 dt=1/4096';
    Nlabels='Nx=Ny=3 Nx=Ny=7 Nx=Ny=15 Nx=Ny=31';
    
    printmat(stabilityExp,'Stable Cases for Explicit Euler Scheme', Nlabels, dtlabels)
    printmat(stabilityImp,'Stable Cases for Implicit Euler Scheme', Nlabels, 'dt=1/64')
end
