function result = sci_comp_ex1()

p = @(t) 10 ./ ( 1 + 9 * exp(-t));
p_dot = @(p) p .* ( 1 - p ./ 10);


min = 0;
max = 1;
delta = 0.01;
p0 = p(min);

p_precise = p(min:delta:max);
p_euler = euler(p_dot, p0, min, max, delta)


end

function output = euler(fun_dot, f0, min, max, delta)

    steps = (max - min) / delta;

    output = zeros(1, steps);
    output(1) = f0;
    for i = 2:steps

        output(i) = output(i - 1) + delta * fun_dot(min + i * delta);  

    end

end