% INITIALIZE ===========================================================

close all
clear

% INPUT PARAMETERS ==================================================

fMHz=2000;  % frequency in MHz


% SECONDARY PARAMETERS ==============================================

f=fMHz*1e6;
lambdac=300/fMHz;
kc=2*pi/lambdac;


% GEOMETRIC INPUTS =================================================

elev=30;     % elevation (deg)
azim=20;     % azimuth  (deg)
D=10;        % distance from MS to building
hRx=1.5;      % MS antenna hieght

% bulding

hB=13;       % bulding height (m)
x1=0;        % beguinning og building
x2=20;       % end of buidling

% mobile route ========================================================

stepxRx=lambdac;
xRx=[-50:stepxRx:50];

% Intersection point of ray with screen =============================

x0=D*tand(azim)+xRx;
y0=D*tand(elev)/cosd(azim)+hRx;
d2=(y0-hRx)/sind(elev);


% Fresnel radius

R1=sqrt(lambdac*abs(d2));

% Integration limits ==================================================

% vertical 

v11=sqrt(2)*(0-y0)/R1;
v21=inf;

v12=sqrt(2)*(hB-y0)/R1;
v22=inf;

v13=sqrt(2)*(0-y0)/R1;
v23=inf;


% horizontal ======================================

u11=-inf;
u21=sqrt(2)*(x1-x0)/R1;

figure,plot(xRx,u21,'k')
xlabel('Traveled distance (m)')
ylabel('u21 parameter')

u12=sqrt(2)*(x1-x0)/R1;
u22=sqrt(2)*(x2-x0)/R1;

figure,plot(xRx,u12,'k',xRx,u22,'k:')
xlabel('Traveled distance (m)')
ylabel('u12 & u22 parameters')
legend('u12','u22')

u13=sqrt(2)*(x2-x0)/R1;
u23=inf;

figure,plot(xRx,u13,'k')
xlabel('Traveled distance (m)')
ylabel('u13 parameter')



% Normalized field strength ============================================

Enormalized=j/2*(...
    ((mfun('FresnelC',u21)-mfun('FresnelC',u11))-j*(mfun('FresnelS',u21)-mfun('FresnelS',u11))).*...
     ((mfun('FresnelC',v21)-mfun('FresnelC',v11))-j*((mfun('FresnelS',v21)-mfun('FresnelS',v11))))...
     +...
    ((mfun('FresnelC',u22)-mfun('FresnelC',u12))-j*(mfun('FresnelS',u22)-mfun('FresnelS',u12))).*...
     ((mfun('FresnelC',v22)-mfun('FresnelC',v12))-j*((mfun('FresnelS',v22)-mfun('FresnelS',v12))))...
     +...
     ((mfun('FresnelC',u23)-mfun('FresnelC',u13))-j*(mfun('FresnelS',u23)-mfun('FresnelS',u13))).*...
     ((mfun('FresnelC',v23)-mfun('FresnelC',v13))-j*((mfun('FresnelS',v23)-mfun('FresnelS',v13))))...
     );


BuidlingOutline=[min(xRx),0; x1,0; x1,hB; x2,hB; x2,0; max(xRx),0];
 
figure,plot(xRx,20*log10(abs(Enormalized)),'k',BuidlingOutline(:,1),BuidlingOutline(:,2),'k:')
xlabel('Traveled distance (m)')
ylabel('Received field strength relative to direct ray (dB)')

figure,plot(xRx,abs(Enormalized),'k',BuidlingOutline(:,1),BuidlingOutline(:,2),'k:')
xlabel('Traveled distance (m)')
ylabel('Received field strength relative to direct ray (lienar units)')

figure,plot(xRx,abs(Enormalized),'k')
xlabel('Traveled distance (m)')
ylabel('Received field strength relative to direct ray (linear units)')
