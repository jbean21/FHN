
tspan = [0 10];
y0 = [0 0];

global I;
I = 0.5;

N = 100;

result = zeros(2, 1);

plotter(tspan, y0, N, 0.01);
y0 = [0 0];
plotter(tspan, y0, N, -0.01);
y0 = [10 -2];
plotter(tspan, y0, N, 0.1);
y0 = [-10 2];
plotter(tspan, y0, N, -0.1);
y0 = [10 -4];
plotter(tspan, y0, N, 0.01);
y0 = [-10 4];
plotter(tspan, y0, N, -0.01);

str = sprintf('Runge Kutta, FitzHugh-Nagumo (I=%d)', I);
title(str);
xlabel('v', 'Interpreter','latex');
ylabel('$\dot{v}$', 'Interpreter','latex');
xlim([-3 3]);
ylim([-2 2]);


function dxdt = v(t,x)

dxdt = zeros(2, 1);

a = 0.7;
b = 0.8;
T = 12.5;
I_ext = 0.5;
dxdt(1) = x(1) - (1/3)*(x(1).^3) - x(2) + I_ext;
dxdt(2) = (1/T)*(x(1)+a-b*x(2));
end


function dvdt = vdot(t, v, w)
I_ext = 0.5;
dvdt = v - (1/3)*(v.^3) - w + I_ext;
end

function dwdt = wdot(t, x)
a = 0.7;
b = 0.8;
T = 12.5;
I_ext = 0.5;
dwdt = (1/T)*(x(1)+a-b*x(2));
end

function p = plotter(tspan, y0, N, incr)
vl = [];
wl = [];
for j = 1:1:N
    [t,y] = ode113(@(t,y) v(t,y), tspan, y0);
    y0(1) = y0(1) + incr;
    plot(y(:,1), vdot(t, y(:,1), y(:,2)));
    hold on
end
hold off
    
end


