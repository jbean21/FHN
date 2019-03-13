
tspan = [0 100];
x0 = [0 0];
opts = odeset('Reltol',1e-13,'AbsTol',1e-14,'Stats','on');
[t,x] = ode113(@(t,x) odeFN(t,x), tspan, x0, opts);
figure(1)
plot(t,x(:,1),'-o');



function dxdt = odeFN(t,x)
dxdt = zeros(2,1);

a = 0.75;%.7;
b = .8;%.8;
T = 3;

% Change this value to test for hopf bifurcation
% It will not be exact due to error from the ODE solver
I = 1.475;

dxdt(1) = T*(x(1) - (1/3)*(x(1).^3) - x(2) + I);
dxdt(2) = (1/T)*(x(1)+a-b*x(2));

end