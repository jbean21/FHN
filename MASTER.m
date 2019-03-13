ti = 0;
tf = 90;
c = 3.0; % frequency
transient_time = 1/c;
tspan = ti:1e-4:tf;
x0 = [0 0 0 0];

path = "";

[t,x] = ode113(@(t,x) odeFN(t,x,c), tspan, x0);

% Remove transient behaviour
[array, index] = min(abs(t-transient_time));
t = t(index:end);
x = [x(index:end,1) x(index:end,2) x(index:end,3) x(index:end,4)];

% Normalise the time array
t = t/c;

%% ====== Power Spectrum ====== %%
% https://community.plm.automation.siemens.com/t5/Testing-Knowledge-Base/Digital-Signal-Processing-Sampling-Rates-Bandwidth-Spectral/ta-p/402991

freq1 = 0:(1/tf):(1/tf)*(length(x(:,1))-1);
ps = fft(x(:,1));
ps1 = ps.*conj(ps);
ps1 = ps1/max(ps1);

% Check is ps1 is even in length. If not, make it even so future division
% by 2 returns and integer value.
if mod(length(ps1), 2) ~= 0
    ps1 = ps1(1:length(ps1)-1);
end

fig = figure(1);
plot(freq1(1:length(ps1)/2),ps1(1:length(ps1)/2));

%title({'Single-Sided Power Spectrum of a FHN neruron' ''})
xlabel('Normalised Frequency', 'interpreter', 'latex')
ylabel('Normalised Power', 'interpreter', 'latex')
xlim([-0.01 0.5]);
ylim([-0.1 1.2]);

box on
set(gca, 'FontSize',[14],'FontName','Times New Roman');%, 'FontWeight', 'bold');
saveas(fig,path+'PowerSpectrum.fig');
saveas(fig,path+'PowerSpectrum.jpg');

%% ====== Amplitude Plot ====== %%
fig = figure(2);
plot(t,x(:,1),t,x(:,3));
xlabel('Normalised Time', 'Interpreter', 'latex');
ylabel('$v_1, v_2$', 'Interpreter', 'latex');

box on
set(gca, 'FontSize',[14],'FontName','Times New Roman');%, 'FontWeight', 'bold');
saveas(fig,path+'CoupledFHN.fig');
saveas(fig,path+'CoupledFHN.jpg');

%% ====== Phase Diagram ====== %%
fig = figure(3);
plot(x(:,1),x(:,2),'k');
hold on
plot(x(:,3),x(:,4),'r');
hold off

xlim([-3 3]);
ylim([-1 2]);
xlabel('$v$', 'Interpreter', 'latex');
ylabel('$w$', 'Interpreter', 'latex');

box on
set(gca, 'FontSize',[14],'FontName','Times New Roman');%, 'FontWeight', 'bold');
saveas(fig,path+'PhaseDiagram.fig');
saveas(fig,path+'PhaseDiagram.jpg');



%% ====== Delay Plot ====== %%

cross = 0;

taus = [];
timepoints = [];
timepoints = [timepoints 0];
threshold = 1e-2;

for i=2:length(x(:,1))-1
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
        
%newtime

for i=1:length(newtime)-2
    if mod(i, 2) == 0
        taus = [taus newtime(i+2)-newtime(i)];
    end
end

%taus

N = 1; % N takes values of 1 to length(taus) - you choose

fig = figure(4);
for i = 1:length(taus)-1
    if i+N < length(taus)
        scatter(taus(i), taus(i+N), 20, [0.5 0.5 0.5], 'filled');
        hold on
    end
end
hold off

xlabel('$\tau_1$', 'Interpreter', 'latex')
ylabel('$\tau_2$', 'Interpreter', 'latex')

box on
set(gca,'FontSize',[14],'FontName','Times New Roman');%, 'FontWeight', 'bold');
saveas(fig,path+'DelayPlot.fig');
saveas(fig,path+'DelayPlot.jpg');


%% ====== Phi1-Phi2 ====== %%

ft1 = fft(x(:,1));
ft2 = fft(x(:,3));

phi1 = wrapTo2Pi(angle(ft1));
phi2 = wrapTo2Pi(angle(ft2));

dphi = phi2-phi1;

fig = figure(5);
plot(freq1(1:length(dphi)),dphi);
xlabel('Normalised Frequency')
ylabel('Phase Difference, $\phi_2-\phi1$', 'Interpreter', 'latex')

box on
set(gca, 'FontSize',14,'FontName','Times New Roman');%, 'FontWeight', 'bold');
saveas(fig,path+'PhaseDifference.fig');
saveas(fig,path+'PhaseDifference.jpg');


%% ====== Power Spectrum varying coupling strength ====== %%
k=0;

c = {};
counter = 1;

while k<1
[t,x] = ode113(@(t,x) coupled(t,x,k), tspan, x0);
freq1 = 0:(1/tf):(1/tf)*(length(x(:,1))-1);

ps1 = fft(x(:,1));
ps1v1 = ps1.*conj(ps1);
ps1v1 = ps1v1/max(ps1v1);
phi1 = angle(ps1);

ps2 = fft(x(:,3));
ps1v2 = ps2.*conj(ps2);
ps1v2 = ps1v2/max(ps1v2);
phi2 = angle(ps2);

dphi = phi2-phi1;

% Check is ps1 is even in length. If not, make it even so future division
% by 2 returns and integer value.
if mod(length(ps1v1), 2) ~= 0
    ps1v1 = ps1v1(1:length(ps1v1)-1);
end

if mod(length(ps1v2), 2) ~= 0
    ps1v2 = ps1v2(1:length(ps1v2)-1);
end

c{counter} = ps1v1;
counter = counter + 1;
k = k + 0.2;
end

fig = figure(6);
plot(freq1(1:length(ps1v1)/2),c{1}(1:length(ps1v1)/2),freq1(1:length(ps1v1)/2),c{2}(1:length(ps1v1)/2),freq1(1:length(ps1v1)/2),c{3}(1:length(ps1v1)/2),freq1(1:length(ps1v1)/2),c{4}(1:length(ps1v1)/2),freq1(1:length(ps1v1)/2),c{5}(1:length(ps1v1)/2));
ylim([-0.1  1.2]);
xlim([-0.01 1]);
xlabel('Normalised Frequency')
ylabel('Normalized Power')

box on
set(gca, 'FontSize',[14],'FontName','Times New Roman');%, 'FontWeight', 'bold');
saveas(fig,path+'PowerSpectrumVaryK.fig');
saveas(fig,path+'PowerSpectrumVaryK.jpg');



%% Coupled Equations
function [dxdt] = odeFN(t,x,c)
dxdt = zeros(4,1);

a = 0.75;
b = 0.8;
I_ext = -0.58;
k = 0.113;

% Reverse unidirectional Coupling
dxdt(1) = c*(x(1) - (1/3)*x(1)^3 + x(2) + I_ext); %x(1) - (1/3)*x(1)^3 - x(2) + I_ext + alpha*(x(3) - x(1));
dxdt(2) = (-1/c)*(x(1) - a + b*x(2));%(1/T)*(x(1) + a - b*x(2));
dxdt(3) = c*(x(3) - (1/3)*x(3)^3 + x(4) + k*(x(1) - x(3)));
dxdt(4) = (-1/c)*(x(3) - a + b*x(4));%(1/T)*(x(3) + a - b*x(4));

end

function [dxdt] = coupled(t,x,k)
dxdt = zeros(4,1);

a = 0.75;
b = 0.8;
c = 3.0;
I_ext = -0.58;

% Reverse unidirectional Coupling
dxdt(1) = c*(x(1) - (1/3)*x(1)^3 + x(2) + I_ext); %x(1) - (1/3)*x(1)^3 - x(2) + I_ext + alpha*(x(3) - x(1));
dxdt(2) = (-1/c)*(x(1) - a + b*x(2));%(1/T)*(x(1) + a - b*x(2));
dxdt(3) = c*(x(3) - (1/3)*x(3)^3 + x(4) + k*(x(1) - x(3)));
dxdt(4) = (-1/c)*(x(3) - a + b*x(4));%(1/T)*(x(3) + a - b*x(4));

end
