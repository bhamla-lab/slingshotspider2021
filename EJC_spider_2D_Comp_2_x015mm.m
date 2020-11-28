clear all
close all
clc

%% Parameters
%Spider
pars.ms = 1.6e-6;%OK-Excel Sheet and Order of magnitude Estimation                                                       %Spider Mass (kg)
ds = 1.3e-3;% From Beast - TRC7 ~1.3mm                                                   %Spider Diameter (m)

pars.Et = 0.4e10;  

% Tension line
dt0 = 1e-6;% Tension Line Diameter                                                             
pars.At = (pi/4)*dt0^2; %Cross-sectional Area (m^2)                                                                                              

% Radial line
dr = 1e-6;                                                                %Diameter (m)
pars.Ar = (pi/4)*dr^2;                                                      %Cross-sectional Area (m^2)
pars.Er = pars.Et;                                                             %Elastic Modulus (Pa)
pars.NumRadLine =8; % There is a ratio between Rad lines to Capture lines 

%Initial Conditions
lt0=4e-2;
pars.lte=1.5*lt0;%Assume lt0=abs(x0) 

lr0 = 2.5e-2; %(Estimated from video)-Initial stretched radial length  %Stretched Length (m)
% lr0 = 2.5e-2;
phir0=30;% From TRC7
% y0=-lr0*sind(phir0);
% x0=1e-3;
% pars.lre=sqrt(-y0^2+(lr0)^2);
pars.lre=0.94*lr0;
% y0=-sqrt(lr0^2-pars.lre^2);
y0=-sqrt(lr0^2-pars.lre^2)*0.9;

y0=-0.01236;
x0=-0.0355e-2;
x0=-0.045e-2;
x0=15*x0;

                                                         
%Damping force terms
SingleRadial=5e-2;
SingleCapLen=6e-3;% Top Capture: 10 mm but close to camera (~8mm), Lower capture is 2.66 mm, Average is ~5mm
NumCap=10;
Lenr = pars.NumRadLine*SingleRadial;% calculate from Images - Function of N - https://www.dropbox.com/work/CHBE-Bhamla-Research/Projects/Slingshot%20Spider/Slingshot%20Spider%20Manuscripts/Invited%20Modeling%20Manuscript/Experimental/Image%20Analysis/Web?preview=EpeirotypusWeb-Fig44-45.psd                                                                 %Combined Length of All Radial Silk (m)
Lenc = (pars.NumRadLine)*SingleCapLen*NumCap;% From Images                                                                 %Combined Length of All Capture Silk (m)
Len = Lenr + Lenc;

muair=1.85e-5;                                                               %viscosity of air (Ns/m^2)
Dwr = 4*pi*muair*Lenr; 
Dwc = 4*pi*muair*Lenc; 

%Spider drag
rhoair=1;
vmaxspider=4;%Oder of magnitude 1m/s
Re=ds*vmaxspider*rhoair/muair;
Cd=24/Re;
Ds=Cd*(pi*(0.5*ds)^2)/2;

%Damping coefficient
pars.cwr = Dwr; 
pars.cwc = Dwc;
pars.cs = Ds;

%% Load Field Data Displacement and time
%%{


load('TCR71A_Kinematics.mat')
xa=X_spider;
ya=Y_spider+abs(Y_spider(1));
ta=time;

% load('TCR71B_Kinematics.mat')
% xb=X_spider;
% yb=Y_spider+abs(Y_spider(1));
% tb=time;
% 
% load('TCR71C_Kinematics.mat')
% xc=X_spider;
% yc=Y_spider+abs(Y_spider(1));
% tc=time;
% 
% load('TCR71E_Kinematics.mat')
% xe=X_spider;
% ye=Y_spider+abs(Y_spider(1));
% te=time;
% 
% load('TCR71F_Kinematics.mat')
% xf=X_spider;
% yf=Y_spider+abs(Y_spider(1));
% tf=time;


maxya_loc=find(ya==max(ya));
% maxyb_loc=find(yb==max(yb));
% maxyc_loc=find(yc==max(yc));
% maxye_loc=find(ye==max(ye));
% maxyf_loc=find(yf==max(yf));


% tmax=max([ta(maxya_loc) tb(maxyb_loc) tc(maxyc_loc) te(maxye_loc) tf(maxyf_loc)]);
% ta=ta+(tmax-ta(maxya_loc));
% tb=tb+(tmax-tb(maxyb_loc));
% tc=tc+(tmax-tc(maxyc_loc));
% te=te+(tmax-te(maxye_loc));
% tf=tf+(tmax-tf(maxyf_loc));

%}
%% Simulation
time = linspace(0,0.1,10000);
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
Kr=Kr*0.42;
Kt=1.25*Kt;
% Forces
Fr1x=-0.125*0.5.*Kr.*delta_lr1.*cosphir1.*heaviside(lr1-pars.lre);
Fr1y=-0.5.*Kr.*delta_lr1.*sinphir1.*heaviside(lr1-pars.lre);

Fr2x=0.125*0.5.*Kr.*delta_lr2.*cosphir2.*heaviside(lr2-pars.lre);
Fr2y=-0.5.*Kr.*delta_lr2.*sinphir2.*heaviside(lr2-pars.lre);

%tension
Ftx=Kt.*delta_lt.*sinphit.*heaviside(lt-pars.lte);
Fty=-Kt.*delta_lt.*cosphit.*heaviside(lt-pars.lte);
%drag
Cwx=-0.05*Cw.*vxs+Cs.*vxs.^2;
Cwy=-0.7*Cw.*vys+Cs.*vys.^2;

blue=[0.3 0.5 1];

%% Plotting
i=0;
i=i+1;
figure(i)
subplot(2,2,1)
plot(t*1e3, (xs)*1e3,'color',[0.3 0.3 0.3],'linewidth',1.5,'displayname','Model')
hold on
plot(ta*1e3, (xa+abs(xa(1)))*1e3,'o','markerfacecolor',[0.2549 0.4118 0.8824],'markeredgecolor',[0.2549 0.4118 0.8824],'markersize',1.75,'linewidth',1.75,'displayname','Field Data')
xlabel('Time, t (ms)')
ylabel('Displacement, x (mm)')
set(gca,'linewidth',1,'fontsize',12)
xlim([0 0.08]*1e3)

subplot(2,2,2)
plot(t*1e3, smooth(ys+abs(ys(1)))*1e3,'color','k','linewidth',1.5,'displayname','Model')
hold on
plot(ta*1e3, ya*1e3,'color','b','linewidth',1.5,'displayname','Model')
xlabel('Time, t (ms)')
ylabel('Displacement, y (mm)')
set(gca,'linewidth',1,'fontsize',12)
xlim([0 0.08]*1e3)

subplot(2,2,3)
title('Model')
plot(xs*1e3, ys*1e3,'-','color',[0.3 0.3 0.3],'linewidth',2,'displayname','Model')
hold on
plot(0, 0,'o','markerfacecolor','k','markeredgecolor','k','markersize',8)
plot(x0*1e3, y0*1e3,'o','markerfacecolor',[1 1 1],'markeredgecolor','r','markersize',8,'linewidth',2)
xlabel('X position (mm)')
ylabel('Y position (mm)')
set(gca,'linewidth',1,'fontsize',12)
xlim([-1 1])
ylim([-0.02 0.005]*1e3)
grid 

subplot(2,2,4)
title('Field')
plot((xa+abs(xa(end)))*1e3, (ya-abs(ya(end)))*1e3,':','color',[0.2549 0.4118 0.8824],'linewidth',2)
hold on
plot((xa(1)+abs(xa(end)))*1e3, (ya(1)-abs(ya(end)))*1e3,'o','markerfacecolor',[1 1 1],'markeredgecolor','r','markersize',8,'linewidth',2)
plot((xa(end)+abs(xa(end)))*1e3, (ya(end)-abs(ya(end)))*1e3,'o','markerfacecolor','k','markeredgecolor','k','markersize',8)
xlabel('X position (mm)')
ylabel('Y position (mm)')
set(gca,'linewidth',1,'fontsize',12)
xlim([-1e-3 1e-3]*1e3)
ylim([-0.02 0.005]*1e3)
grid 


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


figure(3)

subplot(2,1,1)
blue=[0.3 0.5 1];

hold on
plot(ta*1e3, ya*1e3,'o','markerfacecolor',[0.2549 0.4118 0.8824],'markeredgecolor',[0.2549 0.4118 0.8824],'markersize',1.75,'linewidth',1.75,'displayname','Field Data')
plot(t*1e3, (ys+abs(ys(1)))*1e3,'color',[0.3 0.3 0.3],'linewidth',2,'displayname','Model')
xlabel('Time, t (ms)')
ylabel('Displacement, y (mm)')
legend('Location','NorthEastOutside')
set(gca,'linewidth',1,'fontsize',12)
xlim([0 0.08]*1e3)
purple=[0.55 0.17 0.9];
box on

green=[0 0.5 0];
orange=[1 0.5 0];
subplot(2,1,2)
plot(t*1e3, Fr1x,':','color',orange,'linewidth',2,'displayname','Fr1x')
hold on
plot(t*1e3, Fr2x,':','color',green,'linewidth',2,'displayname','Fr2x')
plot(t*1e3, Fr1y,'--','color',orange,'linewidth',2,'displayname','Fr1y')
plot(t*1e3, Fr2y,'--','color',green,'linewidth',2,'displayname','Fr2y')
plot(t*1e3, Ftx,':','color',purple,'linewidth',2,'displayname','Ftx')
plot(t*1e3, Fty,'-','color',purple,'linewidth',2,'displayname','Fty')
plot(t*1e3, Cwx,':','color',[0.3 0.3 0.3],'linewidth',2,'displayname','Fdragx')
plot(t*1e3, Cwy,'-','color',[0.3 0.3 0.3],'linewidth',2,'displayname','Fdragy')

legend('Location','NorthEastOutside')
xlabel('Time, t (ms)')
ylabel('Forces (N)')
set(gca,'linewidth',1,'fontsize',12)
xlim([0 0.08]*1e3)

save('Kt_1.25_x015mm.mat','t','xs','ys')

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
Kr=0.25*Kr;
Kt=1.6*Kt;
% Forces
Fr1x=0.5*Kr*delta_lr1*cosphir1*heaviside(lr1-pars.lre);
Fr1y=0.5*Kr*delta_lr1*sinphir1*heaviside(lr1-pars.lre);

Fr2x=0.5*Kr*delta_lr2*cosphir2*heaviside(lr2-pars.lre);
Fr2y=0.5*Kr*delta_lr2*sinphir2*heaviside(lr2-pars.lre);

Ftx=Kt*delta_lt*sinphit*heaviside(lt-pars.lte);
Fty=Kt*delta_lt*cosphit*heaviside(lt-pars.lte);

% Cwx=0.8*Cw*xdot+0*Cs*xdot^2;
Cwx=0.01*Cw*xdot+Cs*xdot^2;
Cwy=0.7*Cw*ydot+Cs*ydot^2;


%% Equation
dxdt=xdot;
dxdotdt=(-Ftx-0.075*(Fr1x-Fr2x)-Cwx)/ms;
dydt=ydot;
dydotdt=(-Fty-(Fr1y+Fr2y)-Cwy)/ms;

% Fr1x=0.5*Kr*delta_lr*cosphir;
% Fr2x=0.5*Kr*delta_lr*cosphir;
% 
% Fr1y=0.5*Kr*delta_lr*sinphir;
% Fr2y=0.5*Kr*delta_lr*sinphir;
% 
% dxdt=xdot;
% dxdotdt=(-Fr1x*heaviside(x)+Fr2x*heaviside(-x))/ms;
% dydt=ydot;
% dydotdt=(-Kt*delta_lt*heaviside(y)*cosphit-Fr1y-Fr2y-Cw*ydot-Cs*ydot^2)/ms;

%% Output
dynamics=[dxdt;dxdotdt;dydt;dydotdt];
end



