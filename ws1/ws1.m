function [] = ws1()

t_e=5;
Dt=[1,1/2,1/4,1/8];

P_1=[0:Dt(1):t_e; pt(0:Dt(1):t_e);expliciteuler(Dt(1),t_e);heun(Dt(1),t_e);rungekutta(Dt(1),t_e)];
P_2=[0:Dt(2):t_e; pt(0:Dt(2):t_e);expliciteuler(Dt(2),t_e);heun(Dt(2),t_e);rungekutta(Dt(2),t_e)];
P_3=[0:Dt(3):t_e; pt(0:Dt(3):t_e);expliciteuler(Dt(3),t_e);heun(Dt(3),t_e);rungekutta(Dt(3),t_e)];
P_4=[0:Dt(4):t_e; pt(0:Dt(4):t_e);expliciteuler(Dt(4),t_e);heun(Dt(4),t_e);rungekutta(Dt(4),t_e)];

E_exact_1=(Dt(1)/5)*sum((P_1(3:5,:)-[P_1(2,:);P_1(2,:);P_1(2,:)]).^2,2);
E_exact_2=(Dt(2)/5)*sum((P_2(3:5,:)-[P_2(2,:);P_2(2,:);P_2(2,:)]).^2,2);
E_exact_3=(Dt(3)/5)*sum((P_3(3:5,:)-[P_3(2,:);P_3(2,:);P_3(2,:)]).^2,2);
E_exact_4=(Dt(4)/5)*sum((P_4(3:5,:)-[P_4(2,:);P_4(2,:);P_4(2,:)]).^2,2);

E_exact= [E_exact_1, E_exact_2, E_exact_3, E_exact_4];

E_red= [zeros(3,1),E_exact_2./E_exact_1];

fprintf('Explicit Euler method (q=1)\n');
format rat
fprintf('Dt'); disp(Dt)
format short
fprintf('error'); disp(E_exact(1,:))
fprintf('error red.');
fprintf('error app.'); 

fprintf('method of Heun (q=2)\n');
format rat
fprintf('Dt'); disp(Dt)
format short
fprintf('error'); disp(E_exact(2,:))
fprintf('error red.'); 
fprintf('error app.'); 

fprintf('Runge-Kutta method (q=4)\n');
format rat
fprintf('Dt'); disp(Dt)
format short
fprintf('error'); disp(E_exact(2,:))
fprintf('error red.'); 
fprintf('error app.');


figure(1);
hold on
grid on
title('Explicit Euler method (q=1)');
xlabel('t');
ylabel('p(t)');
plot(P_1(1,:),P_1(3,:)); 
plot(P_2(1,:),P_2(3,:)); 
plot(P_3(1,:),P_3(3,:)); 
plot(P_4(1,:),P_4(3,:)); 
plot_p_t(Dt(4),t_e)
legend('Dt=1', 'Dt=1/2', 'Dt=1/4', 'Dt=1/8', 'Exact')

figure(2);
hold on
grid on
title('method of Heun (q=2)');
xlabel('t');
ylabel('p(t)');
plot(P_1(1,:),P_1(4,:)); 
plot(P_2(1,:),P_2(4,:)); 
plot(P_3(1,:),P_3(4,:)); 
plot(P_4(1,:),P_4(4,:)); 
plot_p_t(Dt(4),t_e)
legend('Dt=1', 'Dt=1/2', 'Dt=1/4', 'Dt=1/8', 'Exact')

figure(3);
hold on
grid on
title('Runge-Kutta method (q=4)');
xlabel('t');
ylabel('p(t)');
plot(P_1(1,:),P_1(5,:)); 
plot(P_2(1,:),P_2(5,:)); 
plot(P_3(1,:),P_3(5,:)); 
plot(P_4(1,:),P_4(5,:)); 
plot_p_t(Dt(4),t_e)
legend('Dt=1', 'Dt=1/2', 'Dt=1/4', 'Dt=1/8', 'Exact')

end

function [] = plot_p_t(Dt,t_e)

t=0:Dt:t_e;
p=10./(1+9.*exp(-t));
plot(t,p)

end

function [pt]=pt(t)

pt=10./(1+9*exp(-t));

end
function [dpdt]=dpdt(p)

dpdt=(1-p/10)*p;

end

function [P] = expliciteuler(Dt,t_e)
p=1;
P=0:Dt:t_e;
n=1;
for t=0:Dt:t_e
P(n)=p;
p=p+Dt*dpdt(p);
n=n+1;
end
end

function [P] = heun(Dt,t_e)       
p=1;
P=0:Dt:t_e;
n=1;
for t=0:Dt:t_e
P(n)=p;
dpdt_1=dpdt(p);
dpdt_2=dpdt(p+Dt*dpdt(p));
p=p+Dt*(1/2)*(dpdt_1+dpdt_2);
n=n+1;
end
end

function [P] = rungekutta(Dt,t_e)

p=1;
P=0:Dt:t_e;
n=1;
for t=0:Dt:t_e
P(n)=p;
dpdt_1=dpdt(p);
dpdt_2=dpdt(p+Dt*(1/2)*dpdt_1);
dpdt_3=dpdt(p+Dt*(1/2)*dpdt_2);
dpdt_4=dpdt(p+Dt*dpdt_3);
p=p+Dt*(1/6)*(dpdt_1+2*dpdt_2+2*dpdt_3+dpdt_4);

n=n+1;
end
end
    
