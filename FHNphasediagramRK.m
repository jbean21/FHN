
function rungekutta
h = 0.001;
t = 0;
v0 = 0;
w0 = 0.0;


global I;
I = 0.5;

varray = [];
warray = [];
tarray = [];

dvarray = [];
dwarray = [];


fprintf("Step 0: t = %12.8f, w = %12.8f\n", t, w0);
for i=1:50000
    
    dvarray = [dvarray, dvdt(t,v0,w0)];
    dwarray = [dwarray, dwdt(t,v0,w0)];
   
    
    k1x = h*dvdt(t,v0,w0);
    k1y = h*dwdt(t,v0,w0);
    
    
    k2x = h*dvdt(t+h/2, v0+k1x/2, w0);
    k2y = h*dwdt(t+h/2, v0, w0+k1y/2);
    
    
    k3x = h*dvdt(t+h/2, v0+k2x/2, w0);
    k3y = h*dwdt(t+h/2, v0, w0+k2y/2);
   
    
    k4x = h*dvdt(t+h, v0+k3x, w0);
    k4y = h*dwdt(t+h, v0, w0+k3y);
   
    
    v0 = v0 + (k1x+2*k2x+2*k3x+k4x)/6;
    w0 = w0 + (k1y+2*k2y+2*k3y+k4y)/6;
    
    
    varray = [varray, v0];
    warray = [warray, w0];
    tarray = [tarray, t];
    
    t = t + h;
    
    %fprintf("Step %d: t = %6.4f, x = %18.15f\n, y = %18.15f\n, z = %18.15f\n", i, t, x0, y0, z0);
end
    
%p = plot(tarray, xarray, '-', tarray, yarray, ':', tarray, zarray, '-.','MarkerSize',0.5);
%p = plot(tarray, xarray, '-', tarray, yarray, ':','MarkerSize',0.5);
p = plot(varray,dvarray,'-',warray,dwarray,':');
set(p(1), 'linewidth', 3);
set(p(2), 'linewidth', 3);

str = sprintf('Runge Kutta, FitzHugh-Nagumo (I=%d)', I);
title(str);
xlabel('Function value', 'Interpreter','latex');
ylabel('Derivative value)', 'Interpreter','latex');
%thelegend27 = legend('$v_{dot}(t)$','$w_{dot}(t)$', 'Interpreter','latex', 'location', 'southwest');
%thelegend27.FontSize = 14;

function eqn1 = dvdt(t,v,w)
global I;
% Change I at the top of this file

eqn1 = v - (1/3)*(v.^3) - w +I;

function eqn2 = dwdt(t,v,w)
a = 0.7;
b = 0.8;
T = 12.5;
eqn2 = (1/T)*(v+a-(b*w));

