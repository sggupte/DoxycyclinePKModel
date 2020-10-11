%% Workspace initiation
clear, format short e, figure(1), clf

%% Establishing constants
numDays = 86;
Const = [.11,.593,75,18,0.0055,0.29,0.5445,1.73,11.3,38.8];

tspan = 0:1:numDays*24*2;   %set time from 0 to number of Days       

%% Define i(t)
% Dosing Strategy 1
i_t = zeros(1, numDays*24*2+1);      % 100 mg every 24 hours
for i = 0:(numDays)
    i_t(i*24+1) = 100;
end

% Dosing Strategy 2
i_t2 = zeros(1, numDays*24*2+1);      % 50 mg twice a day
for i = 0:numDays
    i_t2(i*24+1) = 50;
    i_t2(i*24+12+1) = 50;  
end

% Dosing Strategy 3
i_t3 = zeros(1, numDays*24*2+1);      % 200 mg once every other day 
for i = 0:numDays/2
    i_t3(i*48+1) = 200;
end

%% Solving ODE system (before the pathogen)
yinit = [0,0,0];            % Set initial values equal to 0
dayIntroduced = 2;          % Set time we introduce the parasite
i_t = i_t2;                  % Choose your dosing strategy

tspan1 = tspan(1:24*dayIntroduced);
t_i = tspan1;
i_ta = i_t(1:24*dayIntroduced);

DiffFileName = 'myode1';
DE = eval(sprintf('@(t, y, C, t_i, i_ta) %s(t,y,C, t_i, i_ta)', DiffFileName));
[tout1, yout1] = ode45(@(t,y) DE(t,y,Const,t_i,i_ta), tspan1, yinit);

%% Solving Second ODE system (introducing the pathogen)
last = length(yout1(:,1));
initParaNum = 1000;         % Set initial number of parasites

%Sets initial conditions of second set of equations as the last value of 
%the first set of equations 
yinit = [yout1(last,1),yout1(last,2),yout1(last,3),initParaNum]; 
tspan2 = tspan(24*dayIntroduced+1:end);
t_i = tspan2;
i_tb = i_t(24*dayIntroduced+1:end);

DiffFileName = 'myode2';
DE = eval(sprintf('@(t, y, C, t_i, i_tb) %s(t,y,C, t_i, i_tb)', DiffFileName));
[tout2, yout2] = ode45(@(t,y) DE(t,y,Const,t_i,i_tb), tspan2, yinit);

% Add 0 padding to yout1
size_yout = size(yout1);
yout1 = [yout1,zeros(size_yout(1),1)];

tout = [tout1',tout2']';
yout = [yout1',yout2']';
%% Plots
yyaxis left; %% Creates axis measured in mg
plot(tout,yout(:,1),'k-',tout,yout(:,2),'b-',tout,yout(:,3),'r-'); hold on
xlabel('Time (hours)')
ylabel('Amount of Doxycycline in Compartment (mg)')
xlim([1,4129])
ylim([0,140])

% Calculate Values for Table 2
max_liver = max(yout(:,3))*10^(-3)/1071;
max_blood = max(yout(:,2))*10^(-3)/5000;
mean_liver = mean(yout(:,3))*10^(-3)/1071;
mean_blood = mean(yout(:,2))*10^(-3)/5000;

yyaxis right %% Creates axis measured in cells
plot(tout,yout(:,4),'Color',[0,0.5,0],'LineWidth',3);
ylabel('Number of Parasite Cells');
legend('GI','Blood','Liver','Parasite');
ylim([0,5000]);
xlim([0,3900]);

ax = gca;
set(gca,'YColor',[0,0.5,0]);
