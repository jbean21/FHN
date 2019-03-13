% *** FORWORD ***
% The reason for plotting the harmonic oscillator with a specific ang.
% frequency is so that we can determine the order of the time that will be
% required when performing Shampine Gordon on the Hindmarsh Rose equations.
% This is the reason Sandra was talking about units and labelling the axes
% of our graphs.

%Question: What about the constant of integration?

% Initialise variables.
% Note: arrays in MATLAB do not use commas.
tspan = [0 20];
y0=[0 0];
omega = 2;
N = 300;
for j = 1:1:N
    [t,y] = ode45(@(t,y) odefcn(t,y,omega), tspan, y0);
    plot(y(:,1),y(:,2));
    hold on
    y0(1) = y0(1) + 0.04;
end
y0(1) = 0;
for j = 1:1:N
    [t,y] = ode45(@(t,y) odefcn(t,y,omega), tspan, y0);
    plot(y(:,1),y(:,2));
    hold on
    y0(1) = y0(1) - 0.04;
end
y0(1) = -20;
y0(2) = 2;
for j = 1:1:N
    [t,y] = ode45(@(t,y) odefcn(t,y,omega), tspan, y0);
    plot(y(:,1),y(:,2));
    hold on
    y0(1) = y0(1) - 0.1;
end
y0(1) = 20;
y0(2) = -2;
for j = 1:1:N
    [t,y] = ode45(@(t,y) odefcn(t,y,omega), tspan, y0);
    plot(y(:,1),y(:,2));
    hold on
    y0(1) = y0(1) + 0.1;
end
y0(1) = 20;
y0(2) = -4;
for j = 1:1:N
    [t,y] = ode45(@(t,y) odefcn(t,y,omega), tspan, y0);
    plot(y(:,1),y(:,2));
    hold on
    y0(1) = y0(1) + 0.1;
end
y0(1) = -20;
y0(2) = 4;
for j = 1:1:N
    [t,y] = ode45(@(t,y) odefcn(t,y,omega), tspan, y0);
    plot(y(:,1),y(:,2));
    hold on
    y0(1) = y0(1) - 0.1;
end
hold off

title('Harmonic Oscillator Phase Diagram');
xlabel('Angular Displacement $\theta$ /rad', 'Interpreter','latex');
ylabel('Angular Frequency $\dot{\theta}$ /rad s$^{-1}$', 'Interpreter','latex');
xlim([-10 10]);
ylim([-5 5]);
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