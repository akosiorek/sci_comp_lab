%------NewtonMethod--------
% Inputs are timeStepSize, NewtonStopCriteria, p_i and methodName
% p_i and p_j stand respectively for solutions n-th and n+1-th
% iteration and maxIteration stand for iteration # and maximum iterations allowed
% methodName defines which functions G(p) and G'(p) to use.

function [p_j] = NewtonMethod(timeStepSize, NewtonStopCriteria,p_i,methodName)
    difference_ij = realmax; %difference subjected to NewtonStopcriteria
    maxIteration = 1 / NewtonStopCriteria; %Maximum number of iterations
    p_j = p_i; %first root guess
    iteration = 0;
    while ((difference_ij > NewtonStopCriteria) && (iteration<maxIteration)) %evaluate restrictions
        p_j = p_j - G(timeStepSize,p_j,p_i, methodName) / Gp(timeStepSize,p_j,p_i,methodName); %call functions to evaluate next root value.
        difference_ij = abs(p_j - p_i); %evaluate difference for NewtonStopcriteria
        iteration = iteration + 1; % count iterations
    end
end

%G(p) for NewtonMethod
function [G] = G(timeStepSize,p_j,p_i, methodName)
    if (strcmp(methodName,'ImplicitEuler') || strcmp(methodName,'AdamsMoulton'))
        G=p_j-timeStepSize*p_dot(p_j)-p_i; %evaluates depending of current method
    elseif strcmp(methodName,'AdamsMoultonLinearization1')
        G=p_i + (timeStepSize / 2)*(7 * (1 - p_i/10) * p_i + 7 * (1 - p_j/10) * p_i) ;
    elseif strcmp(methodName,'AdamsMoultonLinearization2')
        G=p_i + (timeStepSize / 2)*(7 * (1 - p_i/10) * p_i + 7 * (1 - p_i/10) * p_j) ;
    end
end

%G'(p) for NewtonMethod
function [Gp] = Gp(timeStepSize,p_j,p_i, methodName)
    if (strcmp(methodName,'ImplicitEuler') || strcmp(methodName,'AdamsMoulton'))
        Gp=1-timeStepSize*p_ddot(p_j);  %evaluates depending of current method
    elseif strcmp(methodName,'AdamsMoultonLinearization1')
        Gp=(timeStepSize / 2) * (7 * (p_i / 10));
    elseif strcmp(methodName,'AdamsMoultonLinearization2')
        Gp=(timeStepSize / 2) * (7 * (1 - p_i / 10));
    end
end

%For evaluations (these two should be defined in the main code):
function p_dot = p_dot(p)
    p_dot = 7.*(1-p./10).*p; % dp/dt at point p
end

function p_ddot = p_ddot(p)
    p_ddot = (7-p.*14./10); % d^2p/dt^2 at point p
end

