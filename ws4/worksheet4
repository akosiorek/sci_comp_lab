function [] = worksheet4()

    close all
    clear
    clc

    Lx = 1; Ly = 1; % Set dimensions of domain
    Nx = [3, 7, 15, 31]; Ny = Nx;   % Arrays of Nx, Ny values
    %N = {'3', '7', '15', '31'}; % N: string version
    deltas = 2.^(-4: -1: -10);  % Array of timesteps
    %dels = {'1/64', '1/128', '1/256', '1/512', '1/1024', '1/2048', '1/4096'};   % deltas: string version
    
    dispTimes = [1/8, 2/8, 3/8, 4/8];   % Create array of display times
    dT = {'1/8', '2/8', '3/8', '4/8'};   % Display times: string version
    
    figH = 1: length(dispTimes);    % Create figure handles
    for i = 1: length(figH)
        figure(figH(i));    % Create figures
        set(figH(i), 'Visible', 'off'); % Hide them temporarily
        set(figH(i), 'Position', [25, 25, 1300, 800]);  % Set size
        title(strcat('State of Solutions at t = ', dT(i)));  % Set title
        set(gcf, 'color', 'white');
    end
    
    for i = 1: length(Nx)
        for j = 1: length(deltas)
            T = ones(Ny(i)+2, Nx(i)+2); % Create T-matrix
            T(1:Ny(i)+2, 1) = 0; T(1:Ny(i)+2, Nx(i)+2) = 0; T(1, 1:Nx(i)+2) = 0; T(Ny(i)+2, 1:Nx(i)+2) = 0; % Initializa T-matrix
            time = 0; k = 1;  % Initialize time to zero, set counter for display control
            while time < dispTimes(length(dispTimes))
                T = ExpEulTemporalSolver(deltas(j), time, dispTimes(k), T, Lx, Ly);  % Call Explicit Euler Temporal Solver and update T
                figure(figH(k));    % Activate figure for appropriate dispTime
                index = j + (i-1)*length(deltas); % Set appropriate index in subplot
                PlotResults(T, Nx(i), Ny(i), index, length(Nx), length(deltas));
                time = dispTimes(k);
                k = k + 1;
            end
        end
    end
end
