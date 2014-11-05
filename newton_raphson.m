function root = newton_raphson(x0, max_iter, eps, fun, dfun)

    if nargin < 5
        syms x;
        dfun_symbolic(x) = diff(fun(x));
        dfun = @(y) double(subs(dfun_symbolic, x, y));
    end
    
    x1 = x0;
    error = realmax;
    
    for i = 1:max_iter
        
       x1 = x0;
       x0 = x0 - fun(x0) / dfun(x0);
       %error = fun(x0);
       error = abs(x1 - x0);
       if(error < eps)
           break;
       end
    end
    
    root = x0;
end