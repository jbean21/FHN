
tspan = 0:1e-3:1000;
x0 = [0 0 0 0];

a = 0.75;
b = 0.8;
c = 3.0;
I_ext = -0.58;
k = 0.113;

cross = 0;

[t,x] = ode113(@(t,x) odeFN(t,x), tspan, x0);
figure(1);
plot(t,x(:,1), 'b');
hold on;
plot(t,x(:,3), 'r');
%xlim([0 50]);

%figure(2);
%plot(x(:,1),x(:,2),'g');
%hold on
%plot(x(:,3),x(:,4),'r');
%xlim([-3 3]);
%ylim([-1 2]);
%str = sprintf('ode113, Two Coupled Neurones, Fitzhugh-Nagumo (I=)');
%title(str);
%xlabel('$v$', 'Interpreter', 'latex');
%ylabel('$w$', 'Interpreter', 'latex');
%thelegend27 = legend('x','y','Interpreter','latex');
%thelegend27.FontSize = 14;
%hold off

taus = [];
timepoints = [];
timepoints = [timepoints 0];
threshold = 1e-2;

for i=2:length(x(:,1))-1
    
    %dxdt = (x(i+1,1) - x(i-1,1)) / (2*(t(i+1) - t(i)));
    %abs(x(i,1)-cross)
    if x(i,1) < cross + threshold && cross - threshold < x(i,1)
        timepoints = [timepoints t(i)];
    end
end
timepoints(1) = [];
newtime = [];
i = 1;

while i < length(timepoints)
    counter = 0;
    value = timepoints(i);
    if timepoints(i+1) - value < 1
        for j = i:length(timepoints)-1
            if timepoints(j) - value < 1
                counter = counter + 1;
            end
        end
    
        c = i + round(counter/2);
        if c < length(timepoints) && counter ~= 0
            newtime = [newtime timepoints(c)];
            i = i + counter;
        else
            newtime = [newtime timepoints(i)];
            i = i + 1;
        end
    end
end
        
newtime

for i=1:length(newtime)-2
    if mod(i, 2) == 0
        taus = [taus newtime(i+2)-newtime(i)];
    end
end

%for i=1:length(timepoints)-1
%    taus = [taus timepoints(i+1) - timepoints(i)];
%end
taus
figure(2);
for i = 1:length(taus)-1
        scatter(taus(i), taus(i+1), 20, [0.5 0.5 0.5], 'filled');
        hold on
end
hold off
box on
xlabel('Time /s')
ylabel('Time /s')



function [dxdt] = odeFN(t,x)
dxdt = zeros(4,1);

a = 0.75;
b = 0.8;
c = 3.0;
I_ext = -0.58;
k = 0.113;

% Bi-directional Coupling
dxdt(1) = c*(x(1) - (1/3)*x(1)^3 + x(2) + I_ext + k*(x(3) - x(1))); %x(1) - (1/3)*x(1)^3 - x(2) + I_ext + alpha*(x(3) - x(1));
dxdt(2) = (-1/c)*(x(1) - a + b*x(2));%(1/T)*(x(1) + a - b*x(2));
dxdt(3) = c*(x(3) - (1/3)*x(3)^3 + x(4) + k*(x(1) - x(3)));
dxdt(4) = (-1/c)*(x(3) - a + b*x(4));%(1/T)*(x(3) + a - b*x(4));

end
