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


for i = 1:length(funs)
    fprintf(strcat('\t\t', names{i}, '\n\n'));
    
    
   figure_label = figure(i);    
   hold on
   grid
   for j = 1:length(deltas)
       delta = deltas(j);
       precise = p(min:delta:max);
       approx = funs{i}(p_dot, p0, min, max, delta);
       errors(i, j) = error(precise, approx, delta);
       
       x_vals = min:delta:max;
       plot(x_vals, approx, colours{j});
   end
   plot(x_vals, precise, '-k');
   legend('1', '0.5', '0.25', '0.125', 'precise');
   
    xlabel('time');
    ylabel('value');
    title(names(i));
    hold off
    
    saveas(figure_label, names{i});

    format short
    fprintf('delta'); disp(deltas)
    format shortEng
    fprintf('error'); disp(errors(i, :));
    fprintf('error red.\n');
    fprintf('error app.\n');
    fprintf('\n\n');
end

order = zeros(3, 3);
for i = 2:length(errors)
    order(:, i - 1) = log2(errors(:, i - 1) ./ errors(:, i));
end
mean(order, 2)


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
    e = sqrt((delta / 5) .* sum((approx - precise).^2));
end

function n = steps(max, min, delta) 
    n = (max - min) / delta + 1; 
end