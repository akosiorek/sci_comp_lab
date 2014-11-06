function root = bisection(fun, a, b, min_interval, eps) 

    while(abs(a - b) > min_interval && abs(fun(a)) > eps && abs(fun(b)) > eps)
        
        c = (b + a) / 2;
        fun(c);
        disp([a, b, c]);
        if(fun(a) * fun(c) < 0)
            b = c;
        else
            a = c;
        end
        
    end
    
    if(abs(fun(a)) > abs(fun(b)))
        root = b;
    else
        root = a;
    end
    
    if(root > eps)
        root = NaN;
    end
end