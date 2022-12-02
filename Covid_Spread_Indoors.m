clc;
close all;
clear;
%% CIVE 60005 Coursework
% Controlling outdoor and indoor spread of COVID-19
% 01704926
% Ching Wai Wong

% Constants
G=1/(60*60); %Quanta emmision rate of breathing
T_indoor=17+273; %indoor temperature
T_outdoor=7+273; %outdoor temperature
g=9.81;

%% Q2a maximum number of people without temperature dropping
No=linspace(1,30,30); % Intial guess of 30 people
d=linspace(0.1,0.5,10); % intial guess of d up to 0.5m
t=[0 3600];% Period of 1 hr in seconds
bi=(T_indoor-T_outdoor)*g/T_outdoor;
X=[bi 0 0]; %initial conditions buoyancy concentration probability

% Optimize No and d
R=zeros(length(No),length(d));
for i=1:length(No)
    for j=1:length(d)
        [t,x]=ode45(@(t,x) model(t,x,No(i),d(j),G),t,X);
        
        %Constraints T>17 and p<0.01
        if  min(x(:,1))>=bi &&  max(x(:,3))<0.01
            R(i,j)=d(j);
        else
            R(i,j)=NaN;
        end
    end
end

%Find maximum value of N
dmax=max(R,[],'all');
[row,column]=find(R==dmax);
Nomax=max(row);
fprintf('The max number of people is %g and the depth of window is %g m \n',Nomax,dmax);


% Find the temperature
[t,x]=ode45(@(t,x) model(t,x,Nomax,dmax,G),t,X);
Tmax=((x(:,1)*T_outdoor/g)+T_outdoor)-273;


% plot for concentration probability and temperature against time
figure;
subplot(3,1,1);
plot(t/3600,x(:,2));
xlabel('Time [hr]');
ylabel('Concentration');
ylim([0 0.0015])
title('Concentration against time');

hold on;
subplot(3,1,2);
plot(t/3600,x(:,3));
xlabel('Time [hr]');
ylabel('Probability');
yline(0.01,'--r','LineWidth',1.5);
ylim([0 0.03]);
title('Probability against time');

subplot(3,1,3);
plot(t/3600,Tmax);
xlabel('Time [hr]');
ylabel('Temperature [C]');
ylim([16.5 17.5]);
yline(17,'--r','LineWidth',1.5);
title('Temperature against time');
%% Q2b maximum number of people when time of 0.5 hr is left
No_b=0; % All occupants leave the space
G_b=0; % No occupants therefore no source
t_b=[0 1800]; % 0.5 hr
b_b=((15+273)-T_outdoor)*g/T_outdoor; % Initial condition for buoyancy (guess)
c_b=0.00015; % Initial condition for concentration (guess)
p_b=0; % Initial condition for probabilty (guess)

% Optimize No_b and d
R_b=zeros(length(No),length(d));
for i=1:length(No)
    for j=1:length(d)
        % From t=0 to t=1 hrs (full classroom)
        X=[b_b c_b p_b];
        [t,x]=ode45(@(t,x) model(t,x,No(i),d(j),G),t,X);
        
        % From t=1 to t=1.5 hrs (empty classroom)
        X_b=[x(end,1); x(end,2); x(end,3)];% Boundary conditions from t[0 1]
        [t_b,x_b]=ode45(@(t_b,x_b) model(t_b,x_b,No_b,d(j),G_b),t_b,X_b);
        
        % Repeat until b(0)=b(1.5) and c(0)=c(1.5)
        while abs(x_b(end,1)-x(1,1))>0 && abs(x_b(end,2)-x(1,2))>0
            b_b=x_b(end,1);
            c_b=x_b(end,2);
            X=[b_b c_b p_b];
            [t,x]=ode45(@(t,x) model(t,x,No(i),d(j),G),t,X);
            X_b=[x(end,1); x(end,2); x(end,3)];
            [t_b,x_b]=ode45(@(t_b,x_b) model(t_b,x_b,No_b,d(j),G_b),t_b,X_b);
            %Constraints T>17 and p<0.01
            if min(x(:,1))>=bi  &&  max(x(:,3))<0.01
                R_b(i,j)=d(j);
            else
                R_b(i,j)=NaN;            
            end
        end

    end
end

%Find maximum value of N
dmax_b=max(R_b,[],'all');
[row,column]=find(R_b==dmax_b);
Nomax_b=max(row);
%fprintf('The max number of people is %g and the depth of window is %g m \n',Nomax_b,dmax_b);

%Find the combined result for t[0 1.5]
X=[b_b c_b p_b];
[t,k]=ode45(@(t_b,x_b) model(t_b,x_b,Nomax_b,dmax_b,G),t,X);
X_b=[k(end,1); k(end,2); k(end,3)];
[t_b,k_b]=ode45(@(t_b,k_b) model(t_b,k_b,No_b,dmax_b,G_b),t_b,X_b);
% Repeat until b(0)=b(1.5) and c(0)=c(1.5)
while abs(k_b(end,1)-k(1,1))>0 && abs(k_b(end,2)-k(1,2))>0
    b_b=k_b(end,1);
    c_b=k_b(end,2);
    X=[b_b c_b p_b];
    [t,k]=ode45(@(t_b,x_b) model(t_b,x_b,Nomax_b,dmax_b,G),t,X);
    X_b=[k(end,1); k(end,2); k(end,3)];
    [t_b,k_b]=ode45(@(t_b,k_b) model(t_b,k_b,No_b,dmax_b,G_b),t_b,X_b);
end

T1_b=-273+((k(:,1)*T_outdoor/g)+T_outdoor);
T5_b=-273+((k_b(:,1)*T_outdoor/g)+T_outdoor);
Temp_ans=[T1_b;T5_b(2:end)]; % Compile temperature
t_ans=[t;t_b(2:end)+3600]/3600; % Compile time
b_ans=[k(:,1);k_b((2:end),1)]; % Compile buoyancy
c_ans=[k(:,2);k_b((2:end),2)]; % Compile concentration
p_ans=[k(:,3);k_b((2:end),3)]; % Compile probability

% Plot concentration probability and temperature against time
figure;
subplot(3,1,1);
plot(t_ans,c_ans);
xlabel('Time [hr]');
ylabel('Concentration');
ylim([0 0.002])
title('Concentration against time');

hold on;
subplot(3,1,2);
plot(t_ans,p_ans);
xlabel('Time [hr]');
ylabel('Probability');
yline(0.01,'--r','LineWidth',1.5);
ylim([0 0.03]);
title('Probability against time');

subplot(3,1,3);
plot(t_ans,Temp_ans);
xlabel('Time [hr]');
ylabel('Temperature [C]');
ylim([12 30]);
yline(17,'--r','LineWidth',1.5);
title('Temperature against time');
%% Question 2c
% See report 
%% function
%% Q2a
function fun=model(t,y,No,d,G)
%fun is the buoyancy concentration probability
%t is the time duration
%y is a vector buoyancy concentration and probability
%No is the number of people in the room
%d is the depth of the window

% Constants
V=200; %Volume
w=6; %Width
r=0.0002; %Volume flux of repiration per person
Fp=0.0028; %Equivalent of 100W heat load perperson
Fc=0.028; %Equivalent of 1000W heat load from central heating
F_total=Fp*No+Fc;%Total load
k=0.5;

Q=(k*w*sqrt(y(1)*d^3))/3;
fun(1,1)=(F_total-Q*y(1))/V; %db/dt
fun(2,1)=(G/V)-(Q*y(2)/V); %dc/dt
fun(3,1)=r*(No-1)*(1-y(3))*y(2); %dp/dt

end
