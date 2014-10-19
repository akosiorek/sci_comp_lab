function [] = ws1()

t_e=5;
Dt=[1,1/2,1/4,1/8];

plot_p_t(Dt(end),t_e)

R=[Dt; ones(3,4)];
for k=1:4
R(2,k)=expliciteuler(Dt(k),t_e);
R(3,k)=heun(Dt(k),t_e);
R(4,k)=rungekutta(Dt(k),t_e);
end

printmat(R,'Worksheet 1: Results','Dt Euler Heun RK', '1 1/2 1/4 1/8');

V=10/(1+9*exp(-5));

E=[Dt; ones(3,4)];

for j=1:4
for k=2:4
    E(k,j)=sqrt((E(1,j)/5)*(R(k,j)-V)^2);
end
end

printmat(E,'Worksheet 1: Errors','Dt Exact Euler Heun RK', '1 1/2 1/4 1/8');

Euler=[Dt;E(2,:);zeros(2,4)];
Heun=[Dt;E(3,:);zeros(2,4)];
RK=[Dt;E(4,:);zeros(2,4)];

printmat(Euler,'Explicit Euler method (q=1)','Dt error error_red error_app','1 1/2 1/4 1/8') 
printmat(Heun,'method of Heun (q=2)','Dt error error_red error_app','1 1/2 1/4 1/8') 
printmat(Euler,'Runge-Kutta method (q=4)','Dt error error_red error_app','1 1/2 1/4 1/8') 

end

function [] = plot_p_t(Dt,t_e)

t=0:Dt:t_e;
p=10./(1+9.*exp(-t));
plot(p)
title('p(t)=10/[1+9exp(-t)]');
xlabel('t');
ylabel('p(t)');
grid on

end

function [dpdt]=dpdt(p)

dpdt=(1-p/10)*p;

end

function [p] = expliciteuler(Dt,t_e)

p=1;
for t=0:Dt:t_e
p=p+Dt*dpdt(p);
end

end

function [p] = heun(Dt,t_e)    
    
p=1;
for t=0:Dt:t_e
dpdt_1=dpdt(p);
dpdt_2=dpdt(p+Dt*dpdt(p));
p=p+Dt*(1/2)*(dpdt_1+dpdt_2);
end

end

function [p] = rungekutta(Dt,t_e)

p=1;
for t=0:Dt:t_e
dpdt_1=dpdt(p);
dpdt_2=dpdt(p+Dt*(1/2)*dpdt_1);
dpdt_3=dpdt(p+Dt*(1/2)*dpdt_2);
dpdt_4=dpdt(p+Dt*dpdt_3);
p=p+Dt*(1/6)*(dpdt_1+2*dpdt_2+2*dpdt_3+dpdt_4);
end

end
    
