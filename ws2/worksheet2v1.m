function worksheet2v1()
    
    tic, clear, clc
    funs = {@ExpEuler, @Heun, @ImpEuler, @AdamMoulton, @AdamMoultonL1, @AdamMoultonL2};
    names = {'Explicit Euler method (q = 1)', 'Method of Heun (q = 2)', 'Implicit Euler (q = 1)', 'Adam-Moulton Method (q = 2)', 'Adam-Moulton Linerization 1 (q = 2)', 'Adam-Moulton Linerization 2 (q = 2)'};  % Method names
    deltas = [1, 0.5, 0.25, 0.125, 0.0625, 0.03125]; % Array of timesteps
    colours = {'-r', '-b', '-g', '-m', '-k', '-c', '--k'}; % Colours for graphs
    
    p = @(t) (200 ./ ( 20 - 10 * exp(-7*t)));  % Exact solution of p(t)
    p_dot = @(p) (p .* ( 7 - p .* 0.7));   % Gradient of p(t) at point p
    g = '(0.7*delta)*p^2 + (1 - 7*delta)*p - x0';
    g_dot = '1.4*delta*p + 1 - 7*delta';
    e_NR = 1e-4; % Newton Rhapson Convergence Limit
    
    p0 = 20; min = 0; max = 5;   % Initial value, Start time, End time
    
    errors = zeros(length(funs), length(deltas));   % Array of errors (wrt exact solution)
    error_approx = zeros(size(deltas));             % Array of approx. errors (wrt solution at lowest timestep)
    approxes = cell(length(funs), length(deltas));  % Array of approximations (solutions)
    
    for i = 1: length(funs) % Method Loop
        fprintf(strcat('\t\t', names{i}, '\n\n'));   % Display active method name
        figure(i);   % Assign method handle 'i' to created figure
        hold on; grid; % Allows all plots per method on the same graph, Adds grid
    
        best_error = realmax; best_approx = 0; % Slots for lowest error and corresponding approximate result
        order = zeros(1, length(deltas)); % Array for storing order
   
        for j = 1: length(deltas)    % Timestep loop
            delta = deltas(j);   % bin for current timestep (apart from for legibility purposes, is this needed?)
            precise = p(min: delta: max);    % Caluclate precise solution values at points falling under current timestep
            approxes{i, j} = funs{i}(p_dot, p0, min, max, delta, g, g_dot, e_NR);    % Store resultant approximation array in the cell for method i and stepsize j by calling function for method i
            errors(i, j) = Error(precise, approxes{i, j}, delta);    % Store error in similar slot of 'errors'
            if errors(i, j) < best_error
                best_error = errors(i, j);   % Update minimum error slot-value for current method if new error is lower than 'best_error' so far
                best_approx = approxes{i, j};    % Similarly update best approximation for current method
            end 
            x_vals = min: delta: max;  % Create array for x-axis
            plot(x_vals, approxes{i, j}, colours{j});    % Add approximation for step-size j on held plot in jth colour
        end  % End of timestep loop
   
        plot(x_vals, precise, colours{j+1}); % Add the precise solution on the held plot in 'j+1'th colour
        legend('dt=1', 'dt=0.5', 'dt=0.25', 'dt=0.125', 'dt=0.0625', 'dt=0.03125', 'dt=precise');  % Add legend
        xlabel('time'); ylabel('value'); title(names(i)); % Add x-label, y-label, title (title from array of names)
        hold off     % Free figure for next iteration
        saveas(i, names{i}); % Save figure using array of names
    
        for j = 2: length(deltas)    % Loop for order. Note: Order calculation not possible for t = 0.03125
            order(j) = (errors(i, j - 1) / errors(i, j));
        end
    
        for j = 1: length(deltas)    % Loop for approximate errors
            error_approx(j) = Error(best_approx, approxes{i, j}, deltas(j));
        end
   
        % Presentation of Tabular Results
        format rat, fprintf('delta\t   '); disp(deltas)    % Display timestep sizes
        format shortEng, fprintf('error\t    '); disp(errors(i, :)); % Display errors
        fprintf('error red.'); disp(order);  % Display order
        fprintf('error app.  '); disp(error_approx); % Display approximate errors
        fprintf('\n\n'); % Spacing to visually differentiate b/w outputs of different methods
    end     % End of method loop

    fprintf('Runtime (seconds): '), disp(toc)   % Display runtime
end

function output = ExpEuler(fun_dot, f0, min, max, delta, ~, ~, ~)   % Explicit Euler Solver
    
    N = Steps(max, min, delta); % Call steps for determining number of data points
    output = [f0, zeros(1, N-1)];
    for i = 2: N
        f = output(i - 1);
        output(i) = f + delta * fun_dot(f); % Call fun_dot for derivative
    end
end

function output = Heun(fun_dot, f0, min, max, delta, ~, ~, ~)    % Method of Heun Solver
    
    N = Steps(max, min, delta); % Call steps for determining number of data points
    output = [f0, zeros(1, N-1)];
    for i = 2: N
        f = output(i - 1);
        f_1dot = fun_dot(f);
        f_2dot = fun_dot(f + delta * f_1dot);
        output(i) = f + 0.5 * delta * (f_1dot + f_2dot);
    end
end

function output = ImpEuler(~, f0, min, max, delta, g, gd, e)   % Implicit Euler Solver
    
    N = Steps(max, min, delta); % Call steps for determining number of data points
    output = [f0, zeros(1, N-1)];
    for i = 2: N
        f = output(i - 1);
        output(i) = NR_solver(g, gd, f, e, delta); % NR-Solver for the next point
    end
end

function output = AdamMoulton(fun_dot, f0, min, max, delta, g, gd, e)
    
    N = Steps(max, min, delta); % Call steps for determining number of data points
    output = [f0, zeros(1, N-1)];
    for i = 2: N
        f = output(i - 1);
        fnxt = NR_solver(g, gd, f, e, delta); % NR-Solver for the next point
        output(i) = f + delta * 0.5 * (fun_dot(f) + fun_dot(fnxt)); % Call fun_dot for derivative
    end
end

function output = AdamMoultonL1(fun_dot, f0, min, max, delta, ~, ~, ~)

    N = Steps(max, min, delta); % Call steps for determining number of data points
    output = [f0, zeros(1, N-1)];
    for i = 2: N
        f = output(i - 1);
        fnxt = ((7*delta/20)*f^2 - (7*delta+1)*f)/(7*delta*f/20 + 1);
        output(i) = f + delta * 0.5 * (fun_dot(f) + fun_dot(fnxt));% + (7*(1 - 0.1*fnxt)*f)); % Call fun_dot for derivative
    end
    disp(output)
end

function output = AdamMoultonL2(fun_dot, f0, min, max, delta, ~, ~, ~)
    
    N = Steps(max, min, delta); % Call steps for determining number of data points
    output = [f0, zeros(1, N-1)];
    for i = 2: N
        f = output(i - 1);
        fnxt = ((7*delta/20)*f^2 - (7*0.5*delta+1)*f)/(7*delta*f/20 + 1 - 7*delta/2);
        output(i) = f + delta * 0.5 * (fun_dot(f) + fun_dot(fnxt));% + (7*(1 - 0.1*f)*fnxt)); % Call fun_dot for derivative
    end
end

function x = NR_solver(g, gd, x0, e, delta) % Newton-Rhapson Solver: Determines the next point

    if ((7*delta - 1)^2 + 4*(0.7*delta*x0)) < 0 % Check if solution exists
        x = NaN;
    else
        y = @(p) (eval(g));
        yd = @(p) (eval(gd));
        % dt, although highlighted by MatLab, is used in g, gd.
        x = x0;
        while (abs(y(x)) > e)
            t = x;
            x = t - (y(t)/yd(t));
        end
    end
end

function e = Error(ref, approx, delta)  % Error calculator: Reference, Approxmation, Timestep

    if (length(ref) == length(approx)) %for same timestep size
        e = sqrt((delta / 5) .* sum((approx - ref).^2)); %error calculation
    else    % for different timestep sizes
        longer = ref;
        shorter = approx;
        if length(ref) < length(approx) %Create arrays for adjustment
            longer = approx;
            shorter = ref;
        end
        factor = (length(longer) - 1) / (length(shorter) - 1); % set adjustment factor
        adjusted_longer = zeros(size(shorter));
        for i = 1: length(shorter) %loop for inserting values into adjusted array
            adjusted_longer(i) = longer((i - 1) * factor + 1);
        end
        e = sqrt((delta / 5) .* sum((adjusted_longer - shorter).^2)); %error calculation
    end
end

function stable = Stability(exactSolution, approximation)
%   Stability: Checks if a given vector is a stable approximation of an
%   exact solution

    errors = abs(approximation - exactSolution);
    trendLine = polyfit(1:length(approximation), errors, 1);
    stable = trendLine(1) > 0;
end

function n = Steps(max, min, delta) % Function to calculate number of points for given delta
    n = (max - min) / delta + 1; 
end
