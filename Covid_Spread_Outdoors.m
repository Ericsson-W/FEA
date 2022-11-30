clc;
close all;
clear;
%% CIVE 60005 Coursework
% Controlling outdoor and indoor spread of COVID-19
% 01704926
% Ching Wai Wong

% Constants
DT=1; %Diffusion constant
z=0; %Ground Level
z_h=1.75; %Average head height
M_h=1.75; %Source height of virus
t_d=(0:3*60*60); %Concert duration is 3 hours in seconds
G=98*4/3600; %Quanta emmision rate ,people likely to be singing loudly in a concert
U=3; %Wind speed
% Note it is multiplied by 4 people and converted to per second

%%Q1a calculate the temporal evoution of virus concentration at head height
%at platforms 1,3,7,13 for conditions of no wind
% Note x and y coordinates are measured from platform 13 (the source)
% Platform 1
CP1=Conc_no_wind(-8,8,z_h,M_h,DT,t_d,G);
% Platform 3
CP3=Conc_no_wind(0,8,z_h,M_h,DT,t_d,G);
% Platform 7
CP7=Conc_no_wind(-4,4,z_h,M_h,DT,t_d,G);
% Platform 13
CP13=Conc_no_wind(1,0,z_h,M_h,DT,t_d,G);

% Plot concentrations a function of time
figure;
plot(t_d,CP1,'g',t_d,CP3,'r',t_d,CP7,'b',t_d,CP13,'k','LineWidth',2);
hold on;
legend('Platform 1','Platform 3','Platform 7','Platform 13');
xlabel('Time[s]');
ylabel('Concentration [m^-3]');
title('Concentration Profile at Various Platforms');

%% Q1b Spatial distributions of virusconcentration at no wind and wind with
% speed of 3ms-1 blowing from south-west
% Measuring spread of up to 10m from platform 13
x=linspace(-10,10,100);
y=linspace(-10,10,100);
[X,Y]=meshgrid(x,y);

C_nowind=Conc_wind(x,y,z,M_h,DT,t_d(end),G,0);
C_wind=Conc_wind(x,y,z,M_h,DT,t_d(end),G,U);

% Correct the solution so it is not infinte concentration in the centre
C_max_nowind=Conc_wind(1,0,z,M_h,DT,t_d(end),G,0);
C_max_wind=Conc_wind(1,0,z,M_h,DT,t_d(end),G,U);
for i=1:size(C_nowind,1)
    for j=1:size(C_nowind,1)
        if C_nowind(i,j)>C_max_nowind
           C_nowind(i,j)=C_max_nowind;
        elseif C_wind(i,j)>C_max_wind
           C_wind(i,j)=C_max_wind;
        end
    end
end


% Plot the concentration over area
figure;
contourf(X,Y,C_nowind);
title('Spatial Distribution with no wind (Diffusion Only)');
xlabel('x [m]');
ylabel('y [m]');
colourbar=colorbar;
colourbar.Label.String='concentration [m^-3]';
grid on;
axis equal tight;


figure;
contourf(X,Y,C_wind);
title('Spatial Distribution with wind (Diffusion and Advection)');
xlabel('x [m]');
ylabel('y [m]');
colourbar=colorbar;
colourbar.Label.String='concentration [m^-3]';
grid on;
axis equal tight

%% Q1c Spatial distributions of virus witt disinfectant applied
C_nowind_disinfectant=Conc_disinfectant(x,y,z_h,M_h,DT,t_d(end),G,0);
C_wind_disinfectant=Conc_disinfectant(x,y,z_h,M_h,DT,t_d(end),G,U);

% Correct the solution to avoid singularity
Cmax1=((G/(4*pi*DT*sqrt(1^2+(z_h-M_h)^2))*erfc(sqrt(1^2+(z_h-M_h)^2)/sqrt(4*DT*t_d(end)))-(G/(4*pi*DT*sqrt(1^2 + (z_h+M_h)^2)))*erfc(sqrt(1^2 + (z_h+M_h)^2)/sqrt(4*DT*t_d(end)))));
D=C_nowind_disinfectant>Cmax1;
C_nowind_disinfectant(D)=Cmax1;

Cmax2=Conc_disinfectant(1,0,z,M_h,DT,t_d(end),G,3);
Cmax2 = ((G/(4*pi*DT*sqrt(1^2+(z_h-M_h)^2))*exp(-(U/2*DT)*(sqrt(1^2+(z_h-M_h)^2-1))-(G/(4*pi*DT*sqrt(1^2+(z_h+M_h)^2))*exp(-(U/2*DT)*(sqrt(1^2+(z_h+M_h)^2-1)))))));
D = C_wind_disinfectant>Cmax2;
C_wind_disinfectant(D)=Cmax2;

% Plot concentration over area
figure;
contourf(X,Y,C_nowind_disinfectant);
title('Spatial Distribution with disinfectant and no wind (Diffusion)');
xlabel('x [m]');
ylabel('y [m]');
colourbar=colorbar;colourbar.Label.String='concentration [m^-3]';
grid on;
axis equal tight;

figure;
contourf(X,Y,C_wind_disinfectant);
title('Spatial Distribution with disinfectant and wind (Diffusion and Advection)');
xlabel('x [m]');
ylabel('y [m]');
colourbar=colorbar;colourbar.Label.String='concentration [m^-3]';
grid on;
axis equal tight
%% Q1d Vertical distribution of concetration in front of centre of stage
xs=0; 
ys=15; 
zs=(0:0.1:5.0); 

% Preset sizes for each case
Cs_nowind=zeros(1,length(zs)); %no wind and no disinfectant
Cs_wind=zeros(1,length(zs)); %wind with no disinfectant
Cs_nowind_disinfectant=zeros(1,length(zs));%no wind with disinfectant
Cs_wind_disinfectant=zeros(1,length(zs)); %wind with disinfectant

zs_real=zeros(1,length(zs));%real source height 

zs_image=zeros(1,length(zs));%image source height

% iterate for values of z
for i=1:length(zs)
     zs_real(i)=zs(i)-M_h;
     zs_image(i)=zs(i)+M_h;
     rs_real=sqrt(xs^2 + ys^2 + zs_real.^2);
     rs_image=sqrt(xs^2 + ys^2 + zs_image.^2);

     %no wind and no disinfectant
     Cs_nowind(i)=((G/(4*pi*DT*rs_real(i)))*erfc(rs_real(i)/sqrt(4 *DT*t_d(end)))+(G/(4*pi*DT*rs_image(i)))*erfc(rs_image(i)/sqrt(4*DT*t_d(end))));

     %wind with no disinfectant
     Cs_wind(i)=((G/(4*pi*DT*rs_real(i)))*exp(-(U/2*DT)*(rs_real(i)-ys))+(G/(4*pi*DT*rs_image(i)))*exp(-(U/2*DT)*(rs_image(i)-ys)));

     %no wind with disinfectant
     Cs_nowind_disinfectant(i)=((G/(4*pi*DT*rs_real(i)))*erfc(rs_real(i)/sqrt(4*DT*t_d(end)))-(G/(4*pi*DT*rs_image(i)))*erfc(rs_image(i)/sqrt(4*DT*t_d(end))));

     %wind with disinfectant
     Cs_wind_disinfectant(i)=((G/(4*pi*DT*rs_real(i)))*exp(-(U/2*DT)*(rs_real(i)-ys))-(G/(4*pi*DT*rs_image(i)))*exp(-(U/2*DT)*(rs_image(i)-ys)));
end

% Combined plots of concentration along the height for each case
figure;
plot(Cs_nowind,zs,'k',Cs_wind,zs,'b',Cs_nowind_disinfectant,zs,'r',Cs_wind_disinfectant,zs,'g','LineWiDTh',1)
legend('No wind and no disinfectant','Wind with no disinfectant','No wind with disinfectant','Wind with disinfectant')
title('Vertical distribution in centre of stage');
ylabel('Height from ground level [m]');
xlabel('Concentration [m^-3]');

%% Question 1e Raising each platform to 2m 
% Now z is increashed by 2
ze=z+2;

C_nowind2m=zeros(length(y),length(x)); %no wind no disinfectnat
C_wind2m=zeros(length(y),length(x)); %wind no disinfectnat
C_nowind_disinfectant2m=zeros(length(y),length(x));% no wind with disinfectnat
C_wind_disinfectant2m=zeros(length(y),length(x)); % wind and disinfectant
for i=1:length(x)
     for j=1:length(y)
         z_real=ze-M_h-2;
         z_image=ze+M_h+2; 
         
         r_real=sqrt(x(i)^2 + y(j)^2 + z_real^2);
         r_image=sqrt(x(i)^2 + y(j)^2 + z_image^2);
         
         % no wind
         C_nowind2m(i,j)=((G/(4*pi*DT*r_real))*erfc(r_real/sqrt(4*DT*t_d(end)))+(G/(4*pi*DT*r_image))*erfc(r_image/sqrt(4*DT*t_d(end))));
         % wind
         C_wind2m(i,j)=((G/(4*pi*DT*r_real))*exp(-(U/2*DT)*(r_real-x(i))))+(G/(4*pi*DT*r_image))*exp(-(U/2*DT)*(r_image-x(i)));
         % no wind with disinfectant
         C_nowind_disinfectant2m(i,j)=((G/(4*pi*DT*r_real))*erfc(r_real/sqrt(4*DT*t_d(end)))-(G/(4*pi*DT*r_image))*erfc(r_image/sqrt(4*DT*t_d(end))));
         % wind with disinfectant
         C_wind_disinfectant2m(i,j)=((G/(4*pi*DT*r_real))*exp(-(U/2*DT)*(r_real-x(i)))-(G/(4*pi*DT*r_image))*exp(-(U/2*DT)*(r_image-x(i))));
     
     end         
 end


% Plot concentration over area
figure;
contourf(X,Y,C_nowind2m);
title('Spatial Distribution of raised platform with no wind (Diffusion Only)');
xlabel('x [m]');
ylabel('y [m]');
colourbar=colorbar;
colourbar.Label.String='concentration [m^-3]';
grid on;
axis equal tight;


figure;
contourf(X,Y,C_wind2m);
title('Spatial Distribution of raised platform with wind (Diffusion and Advection)');
xlabel('x [m]');
ylabel('y [m]');
colourbar=colorbar;
colourbar.Label.String='concentration [m^-3]';
grid on;
axis equal tight

figure;
contourf(X,Y,C_nowind_disinfectant2m);
title('Spatial Distribution of raised platform disinfectant and no wind');
xlabel('x [m]');
ylabel('y [m]');
colourbar=colorbar;
colourbar.Label.String='concentration [m^-3]';
grid on;
axis equal tight

figure;
contourf(X,Y,C_wind_disinfectant2m);
title('Spatial Distribution of raised platform with wind and disinfectant');
xlabel('x [m]');
ylabel('y [m]');
colourbar=colorbar;
colourbar.Label.String='concentration [m^-3]';
grid on;
axis equal tight
%% Functions
%% Q1a
function[C]=Conc_no_wind(x,y,z,h,DT,t,G)
% x y z coordinates from source
% h head height
% DT diffusion constant
%t duration considered
%G quanta emmision rate

C=zeros(1,size(t,2));
for i=1:size(t,2)
    z_real=z+h;
    z_image=z-h;
    r_real=(x^2+y^2+z_real^2)^0.5;
    r_image=(x^2+y^2+z_image^2)^0.5;
    % Combining Real and Image source for concentration profile
    C(i)=(G/(4*pi*DT*r_real))*erfc(r_real/(4*DT*t(i)).^0.5)+((G/(4*pi*DT*r_image))*erfc(r_image/(4*DT*t(i)).^0.5));
end
end

%% Q1b
function[C]=Conc_wind(x,y,z,h,DT,t,G,U)
% x y z coordinates from source
% h head height
% DT diffusion constant
%t duration considered
%G quanta emmision rate
%U wind speed
C=zeros(length(y),length(x));

for i=1:length(x)
    for j=1:length(y)
        z_real=z+h;
        z_image=z-h;
        r_real=sqrt(x(i)^2+y(j)^2+z_real^2);
        r_image=sqrt(x(i)^2+y(j)^2+z_image^2);
        if U==0
            C(i,j)=((G/(4*pi*DT*r_real))*erfc(r_real/(4*DT*t).^0.5))+((G/(4*pi*DT*r_image))*erfc(r_image/(4*DT*t).^0.5));
        else
            C(i,j)=(G/(4*pi*DT*r_real))*exp(-(U/(2*DT))*(r_real-x(i)))+(G/(4*pi*DT*r_image))*exp(-(U/(2*DT))*(r_image-x(i)));
        end
    end
end
end

%% Q1c
function[C]=Conc_disinfectant(x,y,z,h,DT,t,G,U)
% x y z coordinates from source
% h head height
% DT diffusion constant
%t duration considered
%G quanta emmision rate
%U wind speed
C=zeros(length(y),length(x));
for i = 1:length(x)
     for j = 1:length(y)
         z_real = z-h;
         z_image = z+h;
         r_real = sqrt(x(i)^2 + y(j)^2 + z_real^2);
         r_image = sqrt(x(i)^2 + y(j)^2 + z_image^2);
         if U==0
            C(i,j)=((G/(4*pi*DT*r_real))*erfc(r_real/sqrt(4*DT*t))-(G/(4*pi*DT*r_image))*erfc(r_image/sqrt(4*DT*t)));
         else
            C(i,j)=((G/(4*pi*DT*r_real))*exp(-(U/2*DT)*(r_real-x(i)))-(G/(4*pi*DT*r_image))*exp(-(U/2*DT)*(r_image-x(i))));
         end
     end
end
end
