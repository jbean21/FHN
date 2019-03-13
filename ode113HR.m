% Solves HR using ode113

global I; 
I = 0;

tspan = [0 50];
x0 = [0 0 0];

[t,x] = ode113(@(t,x) odeHR(t,x), tspan, x0);
plot(t,x(:,1),'-o',t,x(:,2),'-.',t,x(:,3));
str = sprintf('ode113, HindmarshRose (I=%d)', I);
title(str);
xlabel('Time, $t$ /s', 'Interpreter', 'latex');
ylabel('Value (dimensionless)', 'Interpreter', 'latex');
thelegend27 = legend('x','y', 'z','Interpreter','latex');
thelegend27.FontSize = 14;



function dxdt = odeHR(t,x)
dxdt = zeros(3,1);
global I;
a = 1;
b = 3;
c = 1;
d = 5;
r = 0.001;
s = 4;
x_r = -1.6;
dxdt(1) = x(2) + b*(x(1)).^2 - a*(x(1)).^3 - x(3) +I;
dxdt(2) = c-d*(x(1)).^2 - x(2);
dxdt(3) = r*(s*(x(1)-x_r)-x(3));
end

