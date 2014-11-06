function worksheet1()
clc


p = @(t) 10 ./ ( 1 + 9 * exp(-t));
p_dot = @(p) p .* ( 1 - p ./ 10);

p0 = 1;
min = 0;
max = 5;
deltas = [1, 0.5, 0.25, 0.125];



funs = {@euler, @heun, @runge_kutta};
names = {'explicit Euler method (q = 1)', 'method of Heun (q = 2)', 'Runge-Kutta method (q = 4)'};
colours = {'-r', '-b', '-g', '-m'};

errors = zeros(length(funs), length(deltas));
approxes = cell(length(funs), length(deltas));


for i = 1:length(funs)
    fprintf(strcat('\t\t', names{i}, '\n\n'));
    
    
   figure_label = figure(i);    
   hold on
   grid
   
   best_error = realmax;
   best_approx = 0;
   
   for j = 1:length(deltas)
       delta = deltas(j);
       precise = p(min:delta:max);
       approxes{i, j} = funs{i}(p_dot, p0, min, max, delta);
       errors(i, j) = error(precise, approxes{i, j}, delta);
       if errors(i, j) < best_error
           best_error = errors(i, j);
           best_approx = approxes{i, j};
       end
       
       x_vals = min:delta:max;
       plot(x_vals, approxes{i, j}, colours{j});
   end
   
   order = zeros(1, length(deltas));
   order(1) = 0;
   for j = 2:length(deltas)
       order(j) = log2(errors(i, j - 1) / errors(i, j));
   end
   
   error_approx = zeros(size(deltas));
   for j = 1:length(deltas)
       error_approx(j) = error(best_approx, approxes{i, j}, deltas(j));
   end
   
   
   
   plot(x_vals, precise, '-k');
   legend('1', '0.5', '0.25', '0.125', 'precise');
   
    xlabel('time');
    ylabel('value');
    title(names(i));
    hold off
    
    saveas(figure_label, names{i});    
    

    format short
    fprintf('delta\t   '); disp(deltas)
    format shortEng
    fprintf('error\t    '); disp(errors(i, :));
    fprintf('error red.'); disp(order);
    fprintf('error app.  '); disp(error_approx);
    fprintf('\n\n');
end


end

function output = euler(fun_dot, f0, min, max, delta)
    N = steps(max, min, delta);
    output = zeros(1, N);
    output(1) = f0;
    for i = 2:N
        f = output(i - 1);
        output(i) = f + delta * fun_dot(f);
    end
end

function output = heun(fun_dot, f0, min, max, delta)
    N = steps(max, min, delta);
    output = zeros(1, N);
    output(1) = f0;
    for i = 2:N
        f = output(i - 1);
        f_1dot = fun_dot(f);
        f_2dot = fun_dot(f + delta * f_1dot);
        output(i) = f + 0.5 * delta * (f_1dot + f_2dot);
    end
end

function output = runge_kutta(fun_dot, f0, min, max, delta)
    N = steps(max, min, delta);
    output = zeros(1, N);
    output(1) = f0;
    for i = 2:N
        f = output(i - 1);
        f_1dot = fun_dot(f);
        f_2dot = fun_dot(f + 0.5 * delta * f_1dot);
        f_3dot = fun_dot(f + 0.5 * delta * f_2dot);
        f_4dot = fun_dot(f + delta * f_3dot);
        output(i) = f + delta * (f_1dot + 2 * (f_2dot + f_3dot) + f_4dot) / 6;
    end
end

function e = error(precise, approx, delta)

    if (length(precise) == length(approx))
        e = sqrt((delta / 5) .* sum((approx - precise).^2));
    
    else % in case we compare inputs calculated with different deltas (they have different sizes)
        longer = precise;
        shorter = approx;
        
        if length(precise) < length(approx)
            longer = approx;
            shorter = precise;
        end
        
        factor = (length(longer) - 1) / (length(shorter) - 1);
        adjusted_longer = zeros(size(shorter));
        
        for i = 1:length(shorter)
            adjusted_longer(i) = longer((i - 1) * factor + 1);
        end
        
        e = sqrt((delta / 5) .* sum((adjusted_longer - shorter).^2));   
        end
end

function n = steps(max, min, delta) 
    n = (max - min) / delta + 1; 
end