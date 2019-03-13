
tspan = [0 50];
x0 = [0.5 0];

I=-1

[t,x] = ode113(@(t,x) odeFN(t,x), tspan, x0);
plot(t,x(:,1),'-o',t,x(:,2),'-.');
str = sprintf('ode113, Fitzhugh-Nagumo (I=%d)', I);
title(str);
xlabel('Time, $t$ /s', 'Interpreter', 'latex');
ylabel('Value (dimensionless)', 'Interpreter', 'latex');
thelegend27 = legend('x','y','Interpreter','latex');
thelegend27.FontSize = 14;
%xlim([0 30]);


function dxdt = odeFN(t,x)
dxdt = zeros(2,1);

a = 0.7;
b = 0.8;
T = 3;
I_ext = 0.1;
dxdt(1) = T*(x(1) - (1/3)*(x(1).^3) - x(2) - 0);
dxdt(2) = (1/T)*(x(1));

end