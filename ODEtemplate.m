tspan = [0 5];
y0=[0 0.01];
A=1;
B=2;

[t,y] = ode45(@(t,y) odefcn(t,y,A,B), tspan, y0);

plot(t,y(:,1),'-o',t,y(:,2),'-.')
%plot(t,y,'-o')
%title('Solution of van der Pol Equation (\mu = 1) with ODE45');
%xlabel('Time t');
%ylabel('Solution y');
%legend('y_1','y_2')

function dydt = odefcn(t,y,A,B)
dydt = zeros(2,1);
dydt(1) = y(2);
dydt(2) = (A/B)*t.*y(1);
end