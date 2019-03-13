% *** FORWORD ***
% The reason for plotting the harmonic oscillator with a specific ang.
% frequency is so that we can determine the order of the time that will be
% required when performing Shampine Gordon on the Hindmarsh Rose equations.
% This is the reason Sandra was talking about units and labelling the axes
% of our graphs.

%Question: What about the constant of integration?

% Initialise variables.
% Note: arrays in MATLAB do not use commas.
tspan = [0 50];
y0=[0 0.001];
omega = 2;
t_0 = 0;
t_end = 10;

% Solves the ODE and plots the result.
[t,y] = ode45(@(t,y) odefcn(t,y,omega), tspan, y0);

Fs = 500;            % Sampling frequency
T = 1/Fs;             % Sampling period
L = 500;             % Signal Duration (ms)
t = (0:L-1)*T;        % Time vector

% Take FFT
Y = fft(y);

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
plot(f,P1)
%ax = gca;
%ax.YAxis.Exponent = 0;
title({'Single-Sided Amplitude Spectrum of a Simple Pendulum' ''})
xlabel('Frequency /Hz')
%ylabel('Power')

%plot(t,y(:,1),'-o',t,y(:,2),'-.')
%title('Simple Pendulum');
%xlabel('Time $t$ /s', 'Interpreter','latex');
%ylabel('Angular Displacement $\theta$ /rad', 'Interpreter','latex');
%thelegend27 = legend('$\theta(t)$','$\dot{\theta}(t)$', 'Interpreter','latex');
%thelegend27.FontSize = 14;

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
