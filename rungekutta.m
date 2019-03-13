% *** Program that solves an ODE using the fourth-order Runge Kutta method ***
% Adapted from https://math.okstate.edu/people/yqwang/teaching/math4513_fall11/Notes/rungekutta.pdf

function rungekutta
h = 0.0001;
t = 0;
y0 = 0.01;
yarray = [];
tarray = [];
fprintf("Step 0: t = %12.8f, w = %12.8f\n", t, y0);
for i=1:50000
    k1 = h*f(t,y0);
    k2 = h*f(t+h/2, y0+k1/2);
    k3 = h*f(t+h/2, y0+k2/2);
    k4 = h*f(t+h, y0+k3);
    y0 = y0 + (k1+2*k2+2*k3+k4)/6;
    yarray = [yarray, y0];
    tarray = [tarray, t];
    t = t + h;
    fprintf("Step %d: t = %6.4f, w = %18.15f\n", i, t, y0);
    end
plot(tarray, yarray);
title("Runge Kutta, equation dydt = y-t^2+1");
xlabel('Time $t$ /s', 'Interpreter','latex');
ylabel('Function $y(t)$ (dimensionless)', 'Interpreter','latex');

function eqn1 = f(t,y)
eqn1 = y-t^2+1;