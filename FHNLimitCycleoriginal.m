tspan = 0:1e-4:15
y0 = [0.02653 -0.08815];
y1 = [0 1];
a = 0.75;
b = 0.8;
I_ext = 0.5;
 
[t,y] = ode113(@(t,y) odefcn(t,y), tspan, y0);
%[t,x] = ode113(@(t,x) odefcn(t,x), tspan, y1);
 
v = linspace(-10,10);
w = linspace(-10,10);
 
vnull = -(1/3)*(v.^3) + v + I_ext;
wnull = v+I_ext;
 
 
plot(y(:,1), y(:,2),'r');
hold on
%plot(x(:,1), x(:,2),'g');
plot(v,vnull,'k');
plot(v,wnull,'b');
hold off
 
str = sprintf('Phase Diagram, FitzHugh-Nagumo (I=%d), ode113', I_ext);
title(str);
xlabel('$v$', 'Interpreter','latex');
ylabel('$w$', 'Interpreter','latex');
xlim([-2.5 2.5]);
ylim([-0.6 1.8]);
 
 
function dxdt = odefcn(t,x)
 
dxdt = zeros(2, 1);
 
a = 0.75;
b = 0.8;
T = 3.0;
I_ext = 0.5;
dxdt(1) = T*(x(1) - (1/3)*(x(1).^3) - x(2) + I_ext);
dxdt(2) = (1/T)*(x(1)+a-b*x(2));
end