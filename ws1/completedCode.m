function tutorial1()
% Solves the ODE for given initial conditon at different timesteps. Prints
% comparative images for all three methods, displays E, E-tilda, Order.
% Accounts for non-equal divisions by dt.
clear, clc

f = '(1 - (y/10))*y'; y0 = 1; tf = 5;
Rs(4, 3) = struct('Y', 0.0, 'E', 0.0, 'Et', 0.0, 'P', 0.0, 'xaxis', 0.0);
colors = ['-r', '-g', '-m', '-k'];

for mt = 1: 3
    switch mt, case 1, disp('Explicit Euler'), case 2, disp('Method of Heun'), case 3, disp('Fourth-order RK'), end
    disp('E, E-tilda');
    ord = 0;
    for j = 4: -1: 1
        dt = 1/2^(j-1);
        [Rs(j, mt).Y, Rs(j, mt).E, Rs(j, mt).Et, Rs(j, mt).P] = solver(f, y0, dt, tf, mt, Rs);
        if j < 4, ord = ord + (Rs(j, mt).E*(Rs(j+1, mt).E^-1)); end    % Calculate Order
        switch j, case 4, disp('dt = 0.125'), case 3, disp('dt = 0.250'), case 2, disp('dt = 0.500'), case 1, disp('dt = 1.000'), end
        disp(Rs(j, mt).E); disp(Rs(j, mt).Et);
        Rs(j, mt).xaxis = 0: dt: tf;
        plot(Rs(j, mt).xaxis, Rs(j, mt).Y, colors(j)); hold on;
    end
    disp('Order of Growth: ');  disp(round(log2(round(0.3333*ord))));
    plot(Rs(4, mt).xaxis, Rs(4, mt).P); %Add the accurate curve to the plot
    xlabel('t'); ylabel('y'); legend('dt=0.125', 'dt=0.250', 'dt=0.500', 'dt=1.000', 'Precise');
    h = figure; saveas(h, num2str(mt));
    hold off;
end
    

function [y, E, Et, yex] = solver(f, y0, dt, tf, ch, Rs)
% This function offers a choice of three mts to solve the equation:
%                            dy/dt = f(y)
% with inputs f (function, y0 (initial value), dt (timestep), tf (final
% time) and ch (choice of mt). Enter ch = 1 for explicit Euler, ch = 2
% for Heun's mt and ch = 3 for Fourth-order Runge Kutta.

if rem(tf, dt) == 0, n = tf/dt;
else n = floor(tf/dt)+1; end

gr = @(y)(eval(f)); % Read text-fromat function as a math expression
y = [y0, zeros(1, n)]; % Define output scoop
yex = [y0, zeros(1, n)];
d1 = 0; d2 = 0; k = dt/0.125;

switch ch
    case 1
        % Explicit Euler
        for i = 2: n
            y(i) = y(i-1) + dt*(gr(y(i-1)));
            yex(i) = 10/(1+9*exp(-(i-1)*dt));
            d1 = d1 + (yex(i) - y(i))^2;
            if k ~= 1
                d2 = d2 + (Rs(4, ch).Y(k*(i-1)+1) - y(i))^2;
            end
        end
        y(i+1) = y(i) + (tf - dt*(n-1))*gr(y(i)); % Calculate y after last timestep (which is ~= dt)
        yex(i+1) = 10/(1+9*exp(-tf));
        d1 = d1 + (yex(i+1) - y(i+1))^2;
        if k ~= 1
            d2 = d2 + (Rs(4, ch).Y(k*(i)+1) - y(i+1))^2;
        end
        
    case 2
        % Heun's mt
        for i = 2: n
            y(i) = y(i-1) + dt*0.5*(gr(y(i-1)) + gr(y(i-1)+dt*gr(y(i-1))));
            yex(i) = 10/(1+9*exp(-(i-1)*dt));
            d1 = d1 + (yex(i) - y(i))^2;
            if k ~= 1
                d2 = d2 + (Rs(4, ch).Y(k*(i-1)+1) - y(i))^2;
            end
        end
        y(i+1) = y(i) + (tf - dt*(n-1))*0.5*(gr(y(i)) + gr(y(i)+(tf - dt*(n-1))*gr(y(i)))); % Calculate y after last timestep (which is ~= dt)
        yex(i+1) = 10/(1+9*exp(-tf));
        d1 = d1 + (yex(i) - y(i))^2;
        if k ~= 1
            d2 = d2 + (Rs(4, ch).Y(k*(i)+1) - y(i+1))^2;
        end
        
    case 3
        % Fourth-order Runge Kutta
        for i = 2: n
            y1 = gr(y(i-1)); y2 = gr(y(i-1) + (dt/2)*y1); y3 = gr(y(i-1) + (dt/2)*y2); y4 = gr(y(i-1) + dt*y3);
            y(i) = y(i-1) + dt*(1/6)*(y1 + 2*(y2+y3) + y4);
            yex(i) = 10/(1+9*exp(-(i-1)*dt));
            d1 = d1 + (yex(i) - y(i))^2;
            if k ~= 1
                d2 = d2 + (Rs(4, ch).Y(k*(i-1)+1) - y(i))^2;
            end
        end
        dtn = (tf - dt*(n-1)); % Modified timestep for last iteration
        y1 = gr(y(i)); y2 = gr(y(i) + (dtn/2)*y1); y3 = gr(y(i) + (dtn/2)*y2); y4 = gr(y(i) + dtn*y3);
        y(i+1) = y(i) + dt*(1/6)*(y1 + 2*(y2+y3) + y4); % Calculate y after last timestep (which is ~= dt)
        yex(i+1) = 10/(1+9*exp(-tf));
        d1 = d1 + (yex(i) - y(i))^2;
        if k ~= 1
            d2 = d2 + (Rs(4, ch).Y(k*(i)+1) - y(i+1))^2;
        end
end

E = sqrt((dt/5)*d1);
Et = sqrt((dt/5)*d2);
