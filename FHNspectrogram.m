% https://community.plm.automation.siemens.com/t5/Testing-Knowledge-Base/Window-Types-Hanning-Flattop-Uniform-Tukey-and-Exponential/ta-p/445063
% https://community.plm.automation.siemens.com/t5/Testing-Knowledge-Base/Windows-and-Spectral-Leakage/ta-p/432760
tspan = 0:1e-3:200;
x0 = [0 0 0 0];

[t,x] = ode113(@(t,x) odeFN(t,x), tspan, x0);
%figure(1);
%plot(fft(x(:,2)));

str = sprintf('ode113, Fitzhugh-Nagumo (I=0.5)');
title(str);
xlabel('Time, $t$ /s', 'Interpreter', 'latex');
ylabel('Value (dimensionless)', 'Interpreter', 'latex');
%thelegend27 = legend('x','y','Interpreter','latex');
%thelegend27.FontSize = 14;

%figure(2);
%hold on

%plot(fft(x(:,1)));
%length(x(:,1))
spectrogram(x(:,1), 30000, 29000, 'yaxis');%, hamming(1000, 'periodic'));
%hold off

function dxdt = odeFN(t,x)
dxdt = zeros(4,1);

a = 0;%.75;
b = 0;%.8;
c = 3;
I_ext = 0.1;
k = 0.133;

dxdt(1) = c*(x(1) - (1/3)*x(1)^3 + x(2) + I_ext + k*(x(3) - x(1))); %x(1) - (1/3)*x(1)^3 - x(2) + I_ext + alpha*(x(3) - x(1));
dxdt(2) = (-1/c)*(x(1) - a + b*x(2));%(1/T)*(x(1) + a - b*x(2));
dxdt(3) = c*(x(3) - (1/3)*x(3)^3 + x(4) + k*(x(1) - x(3)));
dxdt(4) = (-1/c)*(x(3) - a + b*x(4));%(1/T)*(x(3) + a - b*x(4));



%dxdt(1) = x(1) - (1/3)*(x(1).^3) - x(2) + I_ext;
%dxdt(2) = (1/T)*(x(1)+a-b*x(2));
%dxdt(3) = 0;
%dxdt(4) = 0;
end