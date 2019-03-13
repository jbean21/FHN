
tspan = 0:1e-2:2000;
global I_ext
global a
global b
global gamma
a = -0.01;
b = 0.8;
I_ext = 0.000001;
gamma = 4

%For 3 diagrams, set a=+/- and change y0 to position 1 or 3 on line 27

fixedpoint_v = roots([(-1) (a+1) (-a-b/gamma) (I_ext)]) %Correct
fixedpoint_w = []

fixedpoint_v(imag(fixedpoint_v)~=0) = []

for i=1:length(fixedpoint_v)
    z = fixedpoint_v(i)
    if isreal(z)
        fixedpoint_w = [fixedpoint_w b*z/gamma]
    end
end


y0 = [fixedpoint_v(1), fixedpoint_w(1)]
%y0 = [0 0]
[t,y] = ode113(@(t,y) odefcn(t,y), tspan, y0);
%[t,x] = ode113(@(t,x) odefcn(t,x), tspan, y1);

v = linspace(-10,10,1000);
%w = linspace(-10,10,1000);

stationarypoints_v = roots([-3 2*(a+1) -a])
stationarypoints_w = stationarypoints_v.*(a-stationarypoints_v).*(stationarypoints_v-1)+I_ext;

vnull = v-1/3*v.^3+I_ext;
wnull1 = b*v/gamma;
wnull2 = v/44;
wnull3 = v/100;
%vnull = v.*(a-v).*(v-1);
%wnull = grad*v;


plot(y(:,1), y(:,2),'-.k','linewidth', 3);
hold on
%plot(x(:,1), x(:,2),'g');
plot(v,vnull,'k','linewidth', 5);
%hold on
plot(v,wnull1,'k','linewidth', 5);
%plot(v,wnull2,'-.k', 'Linewidth', 1.2);
%plot(v,wnull3,'k', 'Linewidth', 1.2);
scatter(fixedpoint_v(1),fixedpoint_w(1), 500, 'k', 'filled')
if length(fixedpoint_v)==2
    scatter(fixedpoint_v(2),fixedpoint_w(2), 500, 'k', 'filled')
end
if length(fixedpoint_v)==3
    scatter(fixedpoint_v(3),fixedpoint_w(3), 500, 'k', 'filled')
    scatter(fixedpoint_v(2),fixedpoint_w(2), 500, 'k', 'filled', 'Markerfacecolor', 'w', 'markeredgecolor', 'k')
end
%if length(stationarypoints_v) == 1
%    scatter(stationarypoints_v(1),stationarypoints_w(1), 80, 'b', 'filled')
%end
%if length(stationarypoints_v) == 2
%    scatter(stationarypoints_v(1),stationarypoints_w(1), 80, 'b', 'filled')
%    scatter(stationarypoints_v(2),stationarypoints_w(2), 80, 'b', 'filled')
%end
%plot(0.8509,0.01934, '.k', 'markersize', 18)
hold off

%str = sprintf('Phase Diagram, FitzHugh-Nagumo (I=%d), ode113', I_ext);
%title(str);
xlabel('$v(t)$', 'Interpreter','latex');
ylabel('$w(t)$', 'Interpreter','latex');
set(gca, 'FontSize', 50, 'FontName', 'Times New Roman','linewidth', 5);
xlim([-5 5]);
ylim([-2 2]);
%xlim([-2.5 2.25]);
%ylim([-1.25 1.25]);


function dxdt = odefcn2(t,x)

dxdt = zeros(2, 1);

a = 0.7;
b = 0.8;
T = 3.0;
I_ext = 0;
dxdt(1) = T*(x(1)-1/3*x(1)^3 - x(2) + I_ext);
dxdt(2) = (1/T)*(x(1)+a-b*x(2));
end

function dxdt = odefcn(t,x)
 
dxdt = zeros(2, 1);
global b
global gamma
T = 2;
global I_ext
global a
dxdt(1) = T*(x(1)-1/3*x(1)^3 - x(2) + I_ext);
dxdt(2) = (1/T)*(x(1)+a-b*x(2));
end
