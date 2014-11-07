function stable = stability(exactSolution, approximation)
%   STABILITY Checks if a given vector is a stable approximation of an
%   exact solution
%
%   stable = stability(exactSolution, approxmiation)
%

    errors = abs(approximation - exactSolution);
    trendLine = polyfit(1:length(approximation), errors, 1)
    stable = trendLine(1) > 0;
end

function output = sign_changes(x, n)
    signs = sign(x - means(x, n));
    
    for i=2:length(signs)
        if(abs(signs(i - 1) - signs(i)) < eps)
            signs(i - 1) = 0;
        else
            signs(i - 1) = 1;
        end
    end
    changes = cumsum(signs);
    output = [changes(n+1:length(changes))] - changes(1:length(changes) - n); 
    output = [zeros(1, n) output];
    
end

function output = means(x, n)

    output = zeros(size(x));
    output(n) = sum(output(1:n));
    for i = n+1:length(x)
        output(i) = output(i - 1) + x(i) - x(i - n);
    end
    output = output / n;    
end