
clear all;

tspan = [0 200];
y0 = [0 0];
y1 = [-1 -1];

global I;
I = 0;

N = 1000;

result = zeros(2, 1);

%plotter(tspan, y0, N, 0.01);
%y0 = [0 0];
%plotter(tspan, y0, N, -0.01);
%y0 = [10 -2];
%plotter(tspan, y0, N, 0.1);
%y0 = [-10 2];
%plotter(tspan, y0, N, -0.1);
%y0 = [10 -4];
%plotter(tspan, y0, N, 0.01);
%y0 = [-10 4];
%plotter(tspan, y0, N, -0.01);
[t,y] = ode113(@(t,y) odefcn(t,y), tspan, y0);
[t,x] = ode113(@(t,y) odefcn(t,y), tspan, y1);

plot(y(:,1), y(:,2));
hold on
plot(x(:,1), x(:,2));
hold off

str = sprintf('Phase Diagram, FitzHugh-Nagumo (I=%d), ode113', I);
title(str);
xlabel('$v$', 'Interpreter','latex');
ylabel('$w$', 'Interpreter','latex');
xlim([-3 3]);
ylim([-2 2]);
%xlim([-3 3]);
%ylim([-2 2]);


function dxdt = odefcn(t,x)

dxdt = zeros(2, 1);

a = 0.7;
b = 0.8;
T = 3.0;
I_ext = 0.5;
dxdt(1) = T*(x(1) - (1/3)*(x(1).^3) - x(2) + I_ext);
dxdt(2) = (1/T)*(x(1)+a-b*x(2));
end

%function p = plotter(tspan, y0, N, incr)
%for j = 1:1:N
%    [t,y] = ode113(@(t,y) odefcn(t,y), tspan, y0);
%    plot(y(:,1), y(:,2));
%    hold on
%    y0(1) = y0(1) + incr;
%end
%end

