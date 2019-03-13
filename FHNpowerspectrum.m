% https://community.plm.automation.siemens.com/t5/Testing-Knowledge-Base/Digital-Signal-Processing-Sampling-Rates-Bandwidth-Spectral/ta-p/402991
ti = 0;
tf = 200;
tspan = ti:1e-4:tf;
y0=[0 0 0 0];

[t,y] = ode113(@(t,y) coupled(t,y), tspan, y0);

freq1 = 0:(1/tf):(1/tf)*(length(y(:,1))-1);
ps = fft(y(:,1));
ps1 = ps.*conj(ps);
ps1 = ps1/max(ps1);

% Check is ps1 is even in length. If not, make it even so future division
% by 2 returns and integer value.
if mod(length(ps1), 2) ~= 0
    ps1 = ps1(1:length(ps1)-1);
end

figure(1);
plot(freq1(1:length(ps1)/2),ps1(1:length(ps1)/2));

title({'Single-Sided Power Spectrum of a FHN neruron' ''})
xlabel('Frequency $f$ /Hz', 'interpreter', 'latex')
ylabel('Normalised Power', 'interpreter', 'latex')
xlim([-0.01 0.5]);
ylim([-0.1 1.2]);

function dxdt = odefcn(t,x)

dxdt = zeros(2, 1);

a = 0;%.4;
b = 0;%.8;
T = 10;
I_ext = 0.1;
dxdt(1) = x(1) - (1/3)*(x(1).^3) - x(2) + I_ext;
dxdt(2) = (1/T)*(x(1)+a-b*x(2));
end

function [dxdt] = coupled(t,x)
dxdt = zeros(4,1);

a = 0.75;
b = 0.8;
c = 3.0;
I_ext = -0.58;
k=0.113;

% Bi-directional Coupling
dxdt(1) = c*(x(1) - (1/3)*x(1)^3 + x(2) + I_ext + k*(x(3)-x(1))); %x(1) - (1/3)*x(1)^3 - x(2) + I_ext + alpha*(x(3) - x(1));
dxdt(2) = (-1/c)*(x(1) - a + b*x(2));%(1/T)*(x(1) + a - b*x(2));
dxdt(3) = c*(x(3) - (1/3)*x(3)^3 + x(4) + k*(x(1) - x(3)));
dxdt(4) = (-1/c)*(x(3) - a + b*x(4));%(1/T)*(x(3) + a - b*x(4));

end
