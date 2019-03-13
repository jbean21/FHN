% *** Program that solves the HR equations using the fourth-order Runge Kutta method ***
% Adapted from https://math.okstate.edu/people/yqwang/teaching/math4513_fall11/Notes/rungekutta.pdf

function rungekutta
h = 0.001;
t = 0;
x0 = 0.0;
y0 = 0.0;
z0 = 0.0;

global I;
I = 0;

xarray = [];
yarray = [];
zarray = [];
tarray = [];
fprintf("Step 0: t = %12.8f, w = %12.8f\n", t, y0);
for i=1:50000
    k1x = h*dxdt(t,x0,y0,z0);
    k1y = h*dydt(t,x0,y0,z0);
    k1z = h*dzdt(t,x0,y0,z0);
    
    k2x = h*dxdt(t+h/2, x0+k1x/2, y0, z0);
    k2y = h*dydt(t+h/2, x0, y0+k1y/2, z0);
    k2z = h*dzdt(t+h/2, x0, y0, z0+k1z/2);
    
    k3x = h*dxdt(t+h/2, x0+k2x/2, y0, z0);
    k3y = h*dydt(t+h/2, x0, y0+k2y/2, z0);
    k3z = h*dzdt(t+h/2, x0, y0, z0+k2z/2);
    
    k4x = h*dxdt(t+h, x0+k3x, y0, z0);
    k4y = h*dydt(t+h, x0, y0+k3y, z0);
    k4z = h*dzdt(t+h, x0, y0, z0+k3z);
    
    x0 = x0 + (k1x+2*k2x+2*k3x+k4x)/6;
    y0 = y0 + (k1y+2*k2y+2*k3y+k4y)/6;
    z0 = z0 + (k1z+2*k2z+2*k3z+k4z)/6;
    
    xarray = [xarray, x0];
    yarray = [yarray, y0];
    zarray = [zarray, z0];
    tarray = [tarray, t];
    
    t = t + h;
    
    %fprintf("Step %d: t = %6.4f, x = %18.15f\n, y = %18.15f\n, z = %18.15f\n", i, t, x0, y0, z0);
end
    
p = plot(tarray, xarray, '-', tarray, yarray, ':', tarray, zarray, '-.','MarkerSize',0.5);
%p = plot(tarray, xarray, '-', tarray, yarray, ':','MarkerSize',0.5);
set(p(1), 'linewidth', 3);
set(p(2), 'linewidth', 3);
set(p(3), 'linewidth', 3);
str = sprintf('Runge Kutta, HindmarshRose (I=%d)', I);
title(str);
xlabel('Time $t$ /s', 'Interpreter','latex');
ylabel('Function Value (dimensionless)', 'Interpreter','latex');
thelegend27 = legend('$x(t)$','$y(t)$', '$z(t)$', 'Interpreter','latex', 'location', 'southwest');
thelegend27.FontSize = 14;

function eqn1 = dxdt(t,x,y,z)
global I;
% Change I at the top of this file
a=1;
b=3;
%eqn1 = y-ax^3+bx^2+I-z;
eqn1 = y-a*x.^3+b*x.^2+I-z;

function eqn2 = dydt(t,x,y,z)
c=1;
d=5;
eqn2 = c-d*x.^2-y;

function eqn3 = dzdt(t,x,y,z)
r=0.001;
s=4;
x1=-1.6;
eqn3 = r*(s*(x-x1)-z);