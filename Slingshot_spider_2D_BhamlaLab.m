clear all
close all
clc

%% Parameters
%Spider
pars.ms = 1.6e-6;%OK-Excel Sheet and Order of magnitude Estimation                                                       %Spider Mass (kg)
ds = 1.3e-3;% From Beast - TRC7 ~1.3mm                                                   %Spider Diameter (m)
pars.f=1;
pars.Et = 0.35e10;  
pars.Er = 0.45e10;
% Tension line
dt0 = 1e-6;% Tension Line Diameter                                                             
pars.At = (pi/4)*dt0^2; %Cross-sectional Area (m^2)                                                                                              

% Radial line
dr = 1e-6;                                                                %Diameter (m)
pars.Ar = (pi/4)*dr^2;                                                      %Cross-sectional Area (m^2)                                                             %Elastic Modulus (Pa)
pars.NumRadLine =8; % There is a ratio between Rad lines to Capture lines 

%Initial Conditions
lt0=4e-2;
pars.lte=1.5*lt0;%Assume lt0=abs(x0) 

lr0 = 3.5e-2; %(Estimated from video)-Initial stretched radial length  %Stretched Length (m)
% lr0 = 2.5e-2;
phir0=30;% From TRC7
% y0=-lr0*sind(phir0);
% x0=1e-3;
% pars.lre=sqrt(-y0^2+(lr0)^2);
pars.lre=0.97*lr0;
% y0=-sqrt(lr0^2-pars.lre^2);
y0=-sqrt(lr0^2-pars.lre^2)*0.1;

y0=-0.01236;
y0=-0.012;

x0=-0.0355e-2;
x0=-0.045e-2;
% x0=4*x0;

                                                         
%Damping force terms
SingleRadial=4.5e-2;
SingleCapLen=6.50e-3;% Top Capture: 10 mm but close to camera (~8mm), Lower capture is 2.66 mm, Average is ~5mm

% SingleRadial=3.5e-2;
% SingleCapLen=7e-3;% Top Capture: 10 mm but close to camera (~8mm), Lower capture is 2.66 mm, Average is ~

NumCap=13;
Lenr = pars.NumRadLine*SingleRadial;% calculate from Images - Function of N - https://www.dropbox.com/work/CHBE-Bhamla-Research/Projects/Slingshot%20Spider/Slingshot%20Spider%20Manuscripts/Invited%20Modeling%20Manuscript/Experimental/Image%20Analysis/Web?preview=EpeirotypusWeb-Fig44-45.psd                                                                 %Combined Length of All Radial Silk (m)
Lenc = pars.NumRadLine*SingleCapLen*NumCap;% From Images                                                                 %Combined Length of All Capture Silk (m)
Len = Lenr + Lenc;

muair=1.85e-5;                                                               %viscosity of air (Ns/m^2)
Dwr = 4*pi*muair*Lenr; 
Dwc = 4*pi*muair*Lenc; 

% Spider drag
rhoair=1;
vmaxspider=4;%Oder of magnitude 1m/s
Re=ds*vmaxspider*rhoair/muair;
Cd=1.25;    
Ds=Cd*(pi*(0.5*ds)^2)/2;

% Damping coefficient
pars.cwr = Dwr; 
pars.cwc = Dwc;
pars.cs = Ds;

%% Load Field Data Displacement and time
%%{


load('TCR71A_Kinematics.mat')
xa=X_spider;
ya=Y_spider+abs(Y_spider(1));
ta=time;

maxya_loc=find(ya==max(ya));

%}
%% Simulation
time = linspace(0,0.095,10000);
BC=[x0 0 y0 0];                                                                  %input displacement or initial displacement of the web
[t,dynamics] = ode45(@slingshotspider2D,time,BC,[],pars);

% Kinematics from Simulation
xs = dynamics(:,1);%Spider x- displacement
vxs = dynamics(:,2);%Spider x- velocity
ys = dynamics(:,3);%Spider y- displacement
vys = dynamics(:,4);%Spider y- velocity

%% Forces

lr=sqrt(ys.^2+(pars.lre-xs).^2);
sinphir= ys./(lr);
cosphir=(pars.lre-xs)./lr;

lt=sqrt(xs.^2+(pars.lte+ys).^2);
sinphit=xs./(lt);
cosphit=(pars.lte-ys)/lt;

delta_lr=lr-pars.lre;
delta_lt=lt-pars.lte;

%% Spring Constant and Damping
Kt=(pars.Et*pars.At)/pars.lte;
Kr=(pars.NumRadLine*pars.Er*pars.Ar)/((pars.lre));
Cw=pars.cwr+pars.cwc;
Cs=pars.cs;

%Inertia
ms=pars.ms;

%% Forces
lr1=sqrt(ys.^2+(pars.lre+xs).^2);
lr2=sqrt(ys.^2+(pars.lre-xs).^2);

sinphir1= ys./(lr1);
cosphir1=(pars.lre+xs)./lr1;

sinphir2= ys./(lr2);
cosphir2=(pars.lre-xs)./lr2;

lt=sqrt(xs.^2+(pars.lte+ys).^2);
sinphit=xs./(lt);
cosphit=(pars.lte-ys)./lt;

% Delta l
delta_lr1=lr1-pars.lre;
delta_lr2=lr2-pars.lre;
delta_lt=lt-pars.lte;

%% Spring Constant and Damping
Kt=(pars.Et*pars.At)/pars.lte;
Kr=(pars.NumRadLine*pars.Er*pars.Ar)/((pars.lre));
Cw=pars.cwr+pars.cwc;
Cs=pars.cs;

%Inertia
ms=pars.ms;

% Forces
Fr1x=-0.08*0.5.*Kr.*delta_lr1.*cosphir1.*heaviside(lr1-pars.lre);
Fr1y=-0.5.*Kr.*delta_lr1.*sinphir1.*heaviside(lr1-pars.lre);

Fr2x=0.08*0.5.*Kr.*delta_lr2.*cosphir2.*heaviside(lr2-pars.lre);
Fr2y=-0.5.*Kr.*delta_lr2.*sinphir2.*heaviside(lr2-pars.lre);

%tension
Ftx=Kt.*delta_lt.*sinphit.*heaviside(lt-pars.lte);
Fty=-Kt.*delta_lt.*cosphit.*heaviside(lt-pars.lte);
%drag
Cwx=-sinphir1*Cw.*vxs+Cs.*vxs.^2;
Cwy=-cosphir1*Cw.*vys+Cs.*vys.^2;
Cwx=0.5*Cwx;
Cwy=0.5*Cwy;

blue=[0.3 0.5 1];

%% Plotting
i=0;
i=i+1;
figure(i)
subplot(2,2,1)
plot(t*1e3, (xs)*1e3,'color',[0.3 0.3 0.3],'linewidth',3,'displayname','Model')
hold on
plot(ta*1e3, (xa+abs(xa(1)))*1e3,'o','markerfacecolor',[0.2549 0.4118 0.8824],'markeredgecolor',[0.2549 0.4118 0.8824],'markersize',2,'linewidth',2,'displayname','Field Data')
xlabel('Time, t (ms)')
ylabel('x (mm)')
set(gca,'linewidth',1,'fontsize',18)
xlim([0 0.08]*1e3)

subplot(2,2,2)
plot(t*1e3, smooth(ys+abs(ys(1)))*1e3,'color','k','linewidth',1.5,'displayname','Model')
hold on
plot(ta*1e3, ya*1e3,'color','b','linewidth',1.5,'displayname','Model')
xlabel('Time, t (ms)')
ylabel('Displacement, y (mm)')
set(gca,'linewidth',1,'fontsize',12)
xlim([0 0.08]*1e3)

i=i+1;
figure(i)

plot(xs*1e3, ys*1e3,'-','color',[0.3 0.3 0.3],'linewidth',2,'displayname','Model')
hold on
plot(0, 0,'o','markerfacecolor','k','markeredgecolor','k','markersize',8)
plot(x0*1e3, y0*1e3,'o','markerfacecolor',[1 1 1],'markeredgecolor','r','markersize',8,'linewidth',3)
xlabel('X position (mm)')
ylabel('Y position (mm)')
set(gca,'linewidth',1,'fontsize',14)
xlim([-1 1])
ylim([-0.02 0.005]*1e3)
width=350;
height=350;
set(gcf,'position',[10,10,width,height]) 
i=i+1
figure(i)
plot((xa+abs(xa(end)))*1e3, (ya-abs(ya(end)))*1e3,':','color',[0.2549 0.4118 0.8824],'linewidth',3)
hold on
plot((xa(1)+abs(xa(end)))*1e3, (ya(1)-abs(ya(end)))*1e3,'o','markerfacecolor',[1 1 1],'markeredgecolor','r','markersize',8,'linewidth',3)
plot((xa(end)+abs(xa(end)))*1e3, (ya(end)-abs(ya(end)))*1e3,'o','markerfacecolor','k','markeredgecolor','k','markersize',8)
xlabel('X position (mm)')
ylabel('Y position (mm)')
set(gca,'linewidth',1,'fontsize',14)
xlim([-1e-3 1e-3]*1e3)
ylim([-0.02 0.005]*1e3)
width=350;
height=350;
set(gcf,'position',[10,10,width,height]) 
% figure(2)
% plot(t, Fr1x,'--','color','r','linewidth',1.5,'displayname','Fr1x')
% hold on
% plot(t, Fr2x,'-','color','r','linewidth',1.5,'displayname','Fr2x')
% plot(t, Fr1y,'--','color','b','linewidth',1.5,'displayname','Fr1y')
% plot(t, Fr2y,'-','color','b','linewidth',1.5,'displayname','Fr2y')
% plot(t, Ftx,'--','color','k','linewidth',1.5,'displayname','Ftx')
% plot(t, Fty,'-','color','k','linewidth',1.5,'displayname','Fty')
% plot(t, Cwx,'--','color',[0.5 0.5 0.5],'linewidth',1.5,'displayname','Fdragx')
% plot(t, Cwy,'-','color',[0.5 0.5 0.5],'linewidth',1.5,'displayname','Fdragy')
% legend
% xlabel('Time (s)')
% ylabel('Forces (N)')
% set(gca,'linewidth',1,'fontsize',12)

i=i+1
figure(i)


blue=[0.3 0.5 1];

hold on
plot(ta*1e3, ya*1e3,'o','markerfacecolor',[0.2549 0.4118 0.8824],'markeredgecolor',[0.2549 0.4118 0.8824],'markersize',3,'linewidth',3,'displayname','Field Data')
plot(t*1e3, (ys+abs(ys(1)))*1e3,'color',[0.3 0.3 0.3],'linewidth',3,'displayname','Model')
xlabel('Time, t (ms)')
ylabel('Displacement, y (mm)')
legend('Location','NorthEast')
set(gca,'linewidth',1,'fontsize',16)
xlim([0 0.08]*1e3)
purple=[0.55 0.17 0.9];
ylim([0 20])
box on
width=1000;
height=450;
set(gcf,'position',[10,10,width,height])


i=i+1
figure(i)
green=[0 0.5 0];
orange=[1 0.5 0];

plot(t*1e3, Fr1x,':','color',orange,'linewidth',3,'displayname','Fr1x')
hold on
plot(t*1e3, Fr2x,':','color',green,'linewidth',3,'displayname','Fr2x')
plot(t*1e3, Fr1y,'--','color',orange,'linewidth',3,'displayname','Fr1y')
plot(t*1e3, Fr2y,'--','color',green,'linewidth',3,'displayname','Fr2y')
plot(t*1e3, Ftx,':','color',purple,'linewidth',3,'displayname','Ftx')
plot(t*1e3, Fty,'-','color',purple,'linewidth',3,'displayname','Fty')
plot(t*1e3, Cwx,':','color',[0.3 0.3 0.3],'linewidth',3,'displayname','Fdragx')
plot(t*1e3, Cwy,'-','color',[0.3 0.3 0.3],'linewidth',3,'displayname','Fdragy')

legend('Location','NorthEast')
xlabel('Time, t (ms)')
ylabel('Forces (N)')
set(gca,'linewidth',1,'fontsize',16)
xlim([0 0.08]*1e3)

width=1000;
height=450;
set(gcf,'position',[10,10,width,height])

% filename=sprintf('C_%.02f_x05mm.mat',100);
% save(filename,'t','xs','ys')

function dynamics=slingshotspider2D(t,bc,pars)

%Initial condition
x=bc(1);
xdot=bc(2);
y=bc(3);
ydot=bc(4);

%% Geometry
% Internal parameters
lr1=sqrt(y^2+(pars.lre+x)^2);
lr2=sqrt(y^2+(pars.lre-x)^2);

sinphir1= y/(lr1);
cosphir1=(pars.lre+x)/lr1;

sinphir2= y/(lr2);
cosphir2=(pars.lre-x)/lr2;

lt=sqrt(x^2+(pars.lte+y)^2);
sinphit=x/(lt);
cosphit=(pars.lte-y)/lt;

% Delta l
delta_lr1=lr1-pars.lre;
delta_lr2=lr2-pars.lre;
delta_lt=lt-pars.lte;

%% Spring Constant and Damping
Kt=(pars.Et*pars.At)/pars.lte;
Kr=(pars.NumRadLine*pars.Er*pars.Ar)/((pars.lre));

Cw=pars.cwr+pars.cwc;
Cs=pars.cs;

%Inertia
ms=pars.ms;

% Forces
Fr1x=0.5*Kr*delta_lr1*cosphir1*heaviside(lr1-pars.lre);
Fr1y=0.5*Kr*delta_lr1*sinphir1*heaviside(lr1-pars.lre);

Fr2x=0.5*Kr*delta_lr2*cosphir2*heaviside(lr2-pars.lre);
Fr2y=0.5*Kr*delta_lr2*sinphir2*heaviside(lr2-pars.lre);

Ftx=Kt*delta_lt*sinphit*heaviside(lt-pars.lte);
Fty=Kt*delta_lt*cosphit*heaviside(lt-pars.lte);

Cwx=(Cw*xdot+Cs*xdot^2)*(sinphir1+sinphir2)*0.5;
Cwy=(Cw*ydot+Cs*ydot^2)*(cosphir1+cosphir2)*0.5;

Cwx=0.5*Cwx;% Damping in Half
Cwy=0.5*Cwy;% Damping in Half



%% Equation
dxdt=xdot;
dxdotdt=(-Ftx-0.08*(Fr1x-Fr2x)-Cwx)/ms;
dydt=ydot;
dydotdt=(-Fty-(Fr1y+Fr2y)-Cwy)/ms;

%% Output
dynamics=[dxdt;dxdotdt;dydt;dydotdt];
end



