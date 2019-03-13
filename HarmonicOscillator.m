% *** FORWORD ***
% The reason for plotting the harmonic oscillator with a specific ang.
% frequency is so that we can determine the order of the time that will be
% required when performing Shampine Gordon on the Hindmarsh Rose equations.
% This is the reason Sandra was talking about units and labelling the axes
% of our graphs.

%Question: What about the constant of integration?

% Initialise variables.
% Note: arrays in MATLAB do not use commas.
tspan = [0 10];
y0=[0 0.01];
omega = 2;

% Solves the ODE and plots the result.
[t,y] = ode45(@(t,y) odefcn(t,y,omega), tspan, y0);
plot(t,y(:,1),'-o',t,y(:,2),'-.')
title('Simple Pendulum');
xlabel('Time $t$ /s', 'Interpreter','latex');
ylabel('Angular Displacement $\theta$ /rad', 'Interpreter','latex');
thelegend27 = legend('$\theta(t)$','$\dot{\theta}(t)$', 'Interpreter','latex');
thelegend27.FontSize = 14;

% The following function (called odefcn()) defines a second order ODE in
% two constituent parts, dydt(1) and dydt(2). This linearisation is needed
% for ode45 to work. The parts are presented in a column vector (of size
% 2x1) which is initialised to (0) by the zeros() function.
%                              (0)

function dydt= odefcn(t,y,omega)
dydt = zeros(2,1);
dydt(1) = y(2);
dydt(2) = -(omega^2)*sin(y(1));
end
