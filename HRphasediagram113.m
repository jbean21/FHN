% Solves HR using ode113

global I a b c d r s x_r; 
I = 0;
a = 1;
b = 3;
c = 1;
d = 5;
r = 0.001;
s = 4;
x_r = -1.6;

tspan = [0 50];
x0 = [0 0 0];

[t,values] = ode113(@(t,x) odeHR(t,x), tspan, x0);
x = values(:,1);
y = values(:,2);
z = values(:,3);

dxarray = dx(x,y,z);
dyarray = dy(x,y,z);
dzarray = dz(x,y,z);

%plot(t,x(:,1),'-o',t,x(:,2),'-.',t,x(:,3));
p = plot(values(:,1),dxarray,'-',values(:,2),dyarray,':',values(:,3),dzarray, '-.');
set(p, 'linewidth', 3);
str = sprintf('ode113, HindmarshRose (I=%d), Phase Diagram', I);
title(str);
xlabel('Variable', 'Interpreter', 'latex');
ylabel('Derivative', 'Interpreter', 'latex');
thelegend27 = legend('x','y', 'z','Interpreter','latex');
thelegend27.FontSize = 14;

function HR = odeHR(t,x)
global I a b c d r s x_r;
HR = zeros(3,1);
HR(1) = x(2) + b*(x(1)).^2 - a*(x(1)).^3 - x(3) +I;
HR(2) = c-d*(x(1)).^2 - x(2);
HR(3) = r*(s*(x(1)-x_r)-x(3));
end

function dxdt = dx(x,y,z)
global I a b;
dxdt = y + b*(x).^2 - a*(x).^3 - z +I;
end

function dydt = dy(x,y,z)
global c d;
dydt = c-d*(x).^2 - y;
end

function dzdt = dz(x,y,z)
global r s x_r;
dzdt = r*(s*(x-x_r)-z);
end