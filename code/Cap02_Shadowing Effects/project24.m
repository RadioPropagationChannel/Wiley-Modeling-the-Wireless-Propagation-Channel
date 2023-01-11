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

elev_a=30;     % elevation (deg)
azim_a=10;     % azimuth  (deg)

elev_b=30;
azim_b=0;

elev_c=30;
azim_c=-30;

D=10;        % distance from MS to building
hRx=1.5;      % MS antenna hieght

% buldings

hB1=13;       % bulding height (m)
x1=0;        % beginning of building
x2=20;       % end of building

hB2=10;       % bulding 2 height (m)
x3=30;        % beginning of building 2
x4=45;        % end of building 2

hB3=15;       % bulding 3 height (m)
x5=50;        % beginning of building 3
x6=70;        % end of building 3
 
hB4=12;       % bulding 4 height (m)
x7=80;        % beginning of building 4
x8=90;        % end of building 4

% mobile route ========================================================

stepxRx=lambdac;
xRx=[-20:stepxRx:110];

% Building outline  

BuidlingOutline=[min(xRx),0; x1,0; x1,hB1; x2,hB1; x2,0; ...
    x3,0; x3,hB2; x4,hB2; x4,0; x5,0; x5,hB3; x6,hB3; x6,0; ...
    x7,0; x7,hB4; x8,hB4; x8,0; max(xRx),0];

% Intersection point of ray with screen =============================

x0a=D*tand(azim_a)+xRx;
y0a=D*tand(elev_a)/cosd(azim_a)+hRx;
d2a=(y0a-hRx)/sind(elev_a);

x0b=D*tand(azim_b)+xRx;
y0b=D*tand(elev_b)/cosd(azim_b)+hRx;
d2b=(y0b-hRx)/sind(elev_b);

x0c=D*tand(azim_c)+xRx;
y0c=D*tand(elev_c)/cosd(azim_c)+hRx;
d2c=(y0c-hRx)/sind(elev_c);

% Fresnel radius

R1a=sqrt(lambdac*abs(d2a));
R1b=sqrt(lambdac*abs(d2b)); 
R1c=sqrt(lambdac*abs(d2c));

% Integration limits ==================================================

% vertical 
% link a
v11a=sqrt(2)*(0-y0a)/R1a;
v21a=inf;

v12a=sqrt(2)*(hB1-y0a)/R1a;
v22a=inf;

v13a=sqrt(2)*(0-y0a)/R1a;
v23a=inf;

v14a=sqrt(2)*(hB2-y0a)/R1a;
v24a=inf;

v15a=sqrt(2)*(0-y0a)/R1a;
v25a=inf;

v16a=sqrt(2)*(hB3-y0a)/R1a;
v26a=inf;

v17a=sqrt(2)*(0-y0a)/R1a;
v27a=inf;

v18a=sqrt(2)*(hB4-y0a)/R1a;
v28a=inf;

v19a=sqrt(2)*(0-y0a)/R1a;
v29a=inf;


% link b
v11b=sqrt(2)*(0-y0b)/R1b;
v21b=inf;

v12b=sqrt(2)*(hB1-y0b)/R1b;
v22b=inf;

v13b=sqrt(2)*(0-y0b)/R1b;
v23b=inf;

v14b=sqrt(2)*(hB2-y0b)/R1b;
v24b=inf;

v15b=sqrt(2)*(0-y0b)/R1b;
v25b=inf;

v16b=sqrt(2)*(hB3-y0b)/R1b;
v26b=inf;

v17b=sqrt(2)*(0-y0b)/R1b;
v27b=inf;

v18b=sqrt(2)*(hB4-y0b)/R1b;
v28b=inf;

v19b=sqrt(2)*(0-y0b)/R1b;
v29b=inf;

% link c
v11c=sqrt(2)*(0-y0c)/R1c;
v21c=inf;

v12c=sqrt(2)*(hB1-y0c)/R1c;
v22c=inf;

v13c=sqrt(2)*(0-y0c)/R1c;
v23c=inf;

v14c=sqrt(2)*(hB2-y0c)/R1c;
v24c=inf;

v15c=sqrt(2)*(0-y0c)/R1c;
v25c=inf;

v16c=sqrt(2)*(hB3-y0c)/R1c;
v26c=inf;

v17c=sqrt(2)*(0-y0c)/R1c;
v27c=inf;

v18c=sqrt(2)*(hB4-y0c)/R1c;
v28c=inf;

v19c=sqrt(2)*(0-y0c)/R1c;
v29c=inf;


% horizontal ======================================

% link a
u11a=-inf;
u21a=sqrt(2)*(x1-x0a)/R1a;
% figure,plot(xRx,u21a)

u12a=sqrt(2)*(x1-x0a)/R1a;
u22a=sqrt(2)*(x2-x0a)/R1a;
% figure,plot(xRx,u12a,'r',xRx,u22a,'b')

u13a=sqrt(2)*(x2-x0a)/R1a;
u23a=sqrt(2)*(x3-x0a)/R1a;
% figure,plot(xRx,u13a,'r',xRx,u23a,'b')

u14a=sqrt(2)*(x3-x0a)/R1a;
u24a=sqrt(2)*(x4-x0a)/R1a;
% figure,plot(xRx,u14a,'r',xRx,u24a,'b')

u15a=sqrt(2)*(x4-x0a)/R1a;
u25a=sqrt(2)*(x5-x0a)/R1a;
% figure,plot(xRx,u15a,'r',xRx,u25a,'b')

u16a=sqrt(2)*(x5-x0a)/R1a;
u26a=sqrt(2)*(x6-x0a)/R1a;
% figure,plot(xRx,u16a,'r',xRx,u26a,'b')

u17a=sqrt(2)*(x6-x0a)/R1a;
u27a=sqrt(2)*(x7-x0a)/R1a;
% figure,plot(xRx,u17a,'r',xRx,u27a,'b')

u18a=sqrt(2)*(x7-x0a)/R1a;
u28a=sqrt(2)*(x8-x0a)/R1a;
% figure,plot(xRx,u18a,'r',xRx,u28a,'b')

u19a=sqrt(2)*(x8-x0a)/R1a;
u29a=inf;
% figure,plot(xRx,u19a)

% link b
u11b=-inf;
u21b=sqrt(2)*(x1-x0b)/R1b;
% figure,plot(xRx,u21b)

u12b=sqrt(2)*(x1-x0b)/R1b;
u22b=sqrt(2)*(x2-x0b)/R1b;
% figure,plot(xRx,u12b,'r',xRx,u22b,'b')

u13b=sqrt(2)*(x2-x0b)/R1b;
u23b=sqrt(2)*(x3-x0b)/R1b;
% figure,plot(xRx,u13b,'r',xRx,u23b,'b')

u14b=sqrt(2)*(x3-x0b)/R1b;
u24b=sqrt(2)*(x4-x0b)/R1b;
% figure,plot(xRx,u14b,'r',xRx,u24b,'b')

u15b=sqrt(2)*(x4-x0b)/R1b;
u25b=sqrt(2)*(x5-x0b)/R1b;
% figure,plot(xRx,u15b,'r',xRx,u25b,'b')

u16b=sqrt(2)*(x5-x0b)/R1b;
u26b=sqrt(2)*(x6-x0b)/R1b;
% figure,plot(xRx,u16b,'r',xRx,u26b,'b')

u17b=sqrt(2)*(x6-x0b)/R1b;
u27b=sqrt(2)*(x7-x0b)/R1b;
% figure,plot(xRx,u17b,'r',xRx,u27b,'b')

u18b=sqrt(2)*(x7-x0b)/R1b;
u28b=sqrt(2)*(x8-x0b)/R1b;
% figure,plot(xRx,u18b,'r',xRx,u28b,'b')

u19b=sqrt(2)*(x8-x0b)/R1b;
u29b=inf;
% figure,plot(xRx,u19b)


% link c
u11c=-inf;
u21c=sqrt(2)*(x1-x0c)/R1c;
% figure,plot(xRx,u21c)

u12c=sqrt(2)*(x1-x0c)/R1c;
u22c=sqrt(2)*(x2-x0c)/R1c;
% figure,plot(xRx,u12c,'r',xRx,u22c,'b')

u13c=sqrt(2)*(x2-x0c)/R1c;
u23c=sqrt(2)*(x3-x0c)/R1c;
% figure,plot(xRx,u13c,'r',xRx,u23c,'b')

u14c=sqrt(2)*(x3-x0c)/R1c;
u24c=sqrt(2)*(x4-x0c)/R1c;
% figure,plot(xRx,u14c,'r',xRx,u24c,'b')

u15c=sqrt(2)*(x4-x0c)/R1c;
u25c=sqrt(2)*(x5-x0c)/R1c;
% figure,plot(xRx,u15c,'r',xRx,u25c,'b')

u16c=sqrt(2)*(x5-x0c)/R1c;
u26c=sqrt(2)*(x6-x0c)/R1c;
% figure,plot(xRx,u16c,'r',xRx,u26c,'b')

u17c=sqrt(2)*(x6-x0c)/R1c;
u27c=sqrt(2)*(x7-x0c)/R1c;
% figure,plot(xRx,u17c,'r',xRx,u27c,'b')

u18c=sqrt(2)*(x7-x0c)/R1c;
u28c=sqrt(2)*(x8-x0c)/R1c;
% figure,plot(xRx,u18c,'r',xRx,u28c,'b')

u19c=sqrt(2)*(x8-x0c)/R1c;
u29c=inf;
% figure,plot(xRx,u19c)


% Normalized field strength ==========================================

% link a

Enormalized_a=j/2*(...
    ((mfun('FresnelC',u21a)-mfun('FresnelC',u11a))-j*(mfun('FresnelS',u21a)-mfun('FresnelS',u11a))).*...
     ((mfun('FresnelC',v21a)-mfun('FresnelC',v11a))-j*((mfun('FresnelS',v21a)-mfun('FresnelS',v11a))))...
     +...
    ((mfun('FresnelC',u22a)-mfun('FresnelC',u12a))-j*(mfun('FresnelS',u22a)-mfun('FresnelS',u12a))).*...
     ((mfun('FresnelC',v22a)-mfun('FresnelC',v12a))-j*((mfun('FresnelS',v22a)-mfun('FresnelS',v12a))))...
     +...
     ((mfun('FresnelC',u23a)-mfun('FresnelC',u13a))-j*(mfun('FresnelS',u23a)-mfun('FresnelS',u13a))).*...
     ((mfun('FresnelC',v23a)-mfun('FresnelC',v13a))-j*((mfun('FresnelS',v23a)-mfun('FresnelS',v13a))))...
     +...
     ((mfun('FresnelC',u24a)-mfun('FresnelC',u14a))-j*(mfun('FresnelS',u24a)-mfun('FresnelS',u14a))).*...
     ((mfun('FresnelC',v24a)-mfun('FresnelC',v14a))-j*((mfun('FresnelS',v24a)-mfun('FresnelS',v14a))))...
     +...
     ((mfun('FresnelC',u25a)-mfun('FresnelC',u15a))-j*(mfun('FresnelS',u25a)-mfun('FresnelS',u15a))).*...
     ((mfun('FresnelC',v25a)-mfun('FresnelC',v15a))-j*((mfun('FresnelS',v25a)-mfun('FresnelS',v15a))))...
     +...
     ((mfun('FresnelC',u26a)-mfun('FresnelC',u16a))-j*(mfun('FresnelS',u26a)-mfun('FresnelS',u16a))).*...
     ((mfun('FresnelC',v26a)-mfun('FresnelC',v16a))-j*((mfun('FresnelS',v26a)-mfun('FresnelS',v16a))))...
     +...
     ((mfun('FresnelC',u27a)-mfun('FresnelC',u17a))-j*(mfun('FresnelS',u27a)-mfun('FresnelS',u17a))).*...
     ((mfun('FresnelC',v27a)-mfun('FresnelC',v17a))-j*((mfun('FresnelS',v27a)-mfun('FresnelS',v17a))))...
     +...
     ((mfun('FresnelC',u28a)-mfun('FresnelC',u18a))-j*(mfun('FresnelS',u28a)-mfun('FresnelS',u18a))).*...
     ((mfun('FresnelC',v28a)-mfun('FresnelC',v18a))-j*((mfun('FresnelS',v28a)-mfun('FresnelS',v18a))))...
     +...
     ((mfun('FresnelC',u29a)-mfun('FresnelC',u19a))-j*(mfun('FresnelS',u29a)-mfun('FresnelS',u19a))).*...
     ((mfun('FresnelC',v29a)-mfun('FresnelC',v19a))-j*((mfun('FresnelS',v29a)-mfun('FresnelS',v19a))))...
     );


figure,plot(xRx,20*log10(abs(Enormalized_a)),'k',BuidlingOutline(:,1),BuidlingOutline(:,2),'k:')
xlabel('Traveled distance (m)')
ylabel('Received field strength relative to direct ray (dB)')
title(['Elevation ' num2str(elev_a) ', Orientation ' num2str(azim_a)])

figure,plot(xRx,abs(Enormalized_a),'k',BuidlingOutline(:,1),BuidlingOutline(:,2),'k:')
xlabel('Traveled distance (m)')
ylabel('Received field strength relative to direct ray (linear units)')
title(['Elevation ' num2str(elev_a) ', Orientation ' num2str(azim_a)])

figure,plot(xRx,abs(Enormalized_a),'k')
xlabel('Traveled distance (m)')
ylabel('Received field strength relative to direct ray (linear units)')
title(['Elevation ' num2str(elev_a) ', Orientation ' num2str(azim_a)])

% link b
Enormalized_b=j/2*(...
    ((mfun('FresnelC',u21b)-mfun('FresnelC',u11b))-j*(mfun('FresnelS',u21b)-mfun('FresnelS',u11b))).*...
     ((mfun('FresnelC',v21b)-mfun('FresnelC',v11b))-j*((mfun('FresnelS',v21b)-mfun('FresnelS',v11b))))...
     +...
    ((mfun('FresnelC',u22b)-mfun('FresnelC',u12b))-j*(mfun('FresnelS',u22b)-mfun('FresnelS',u12b))).*...
     ((mfun('FresnelC',v22b)-mfun('FresnelC',v12b))-j*((mfun('FresnelS',v22b)-mfun('FresnelS',v12b))))...
     +...
     ((mfun('FresnelC',u23b)-mfun('FresnelC',u13b))-j*(mfun('FresnelS',u23b)-mfun('FresnelS',u13b))).*...
     ((mfun('FresnelC',v23b)-mfun('FresnelC',v13b))-j*((mfun('FresnelS',v23b)-mfun('FresnelS',v13b))))...
     +...
     ((mfun('FresnelC',u24b)-mfun('FresnelC',u14b))-j*(mfun('FresnelS',u24b)-mfun('FresnelS',u14b))).*...
     ((mfun('FresnelC',v24b)-mfun('FresnelC',v14b))-j*((mfun('FresnelS',v24b)-mfun('FresnelS',v14b))))...
     +...
     ((mfun('FresnelC',u25b)-mfun('FresnelC',u15b))-j*(mfun('FresnelS',u25b)-mfun('FresnelS',u15b))).*...
     ((mfun('FresnelC',v25b)-mfun('FresnelC',v15b))-j*((mfun('FresnelS',v25b)-mfun('FresnelS',v15b))))...
     +...
     ((mfun('FresnelC',u26b)-mfun('FresnelC',u16b))-j*(mfun('FresnelS',u26b)-mfun('FresnelS',u16b))).*...
     ((mfun('FresnelC',v26b)-mfun('FresnelC',v16b))-j*((mfun('FresnelS',v26b)-mfun('FresnelS',v16b))))...
     +...
     ((mfun('FresnelC',u27b)-mfun('FresnelC',u17b))-j*(mfun('FresnelS',u27b)-mfun('FresnelS',u17b))).*...
     ((mfun('FresnelC',v27b)-mfun('FresnelC',v17b))-j*((mfun('FresnelS',v27b)-mfun('FresnelS',v17b))))...
     +...
     ((mfun('FresnelC',u28b)-mfun('FresnelC',u18b))-j*(mfun('FresnelS',u28b)-mfun('FresnelS',u18b))).*...
     ((mfun('FresnelC',v28b)-mfun('FresnelC',v18b))-j*((mfun('FresnelS',v28b)-mfun('FresnelS',v18b))))...
     +...
     ((mfun('FresnelC',u29b)-mfun('FresnelC',u19b))-j*(mfun('FresnelS',u29b)-mfun('FresnelS',u19b))).*...
     ((mfun('FresnelC',v29b)-mfun('FresnelC',v19b))-j*((mfun('FresnelS',v29b)-mfun('FresnelS',v19b))))...
     );


figure,plot(xRx,20*log10(abs(Enormalized_b)),'k',BuidlingOutline(:,1),BuidlingOutline(:,2),'k:')
xlabel('Traveled distance (m)')
ylabel('Received field strength relative to direct ray (dB)')
title(['Elevation ' num2str(elev_b) ', Orientation ' num2str(azim_b)])

figure,plot(xRx,abs(Enormalized_b),'k',BuidlingOutline(:,1),BuidlingOutline(:,2),'k:')
xlabel('Traveled distance (m)')
ylabel('Received field strength relative to direct ray (linear units)')
title(['Elevation ' num2str(elev_b) ', Orientation ' num2str(azim_b)])

figure,plot(xRx,abs(Enormalized_b),'k')
xlabel('Traveled distance (m)')
ylabel('Received field strength relative to direct ray (linear units)')
title(['Elevation ' num2str(elev_b) ', Orientation ' num2str(azim_b)])

% link c
Enormalized_c=j/2*(...
    ((mfun('FresnelC',u21c)-mfun('FresnelC',u11c))-j*(mfun('FresnelS',u21c)-mfun('FresnelS',u11c))).*...
     ((mfun('FresnelC',v21c)-mfun('FresnelC',v11c))-j*((mfun('FresnelS',v21c)-mfun('FresnelS',v11c))))...
     +...
    ((mfun('FresnelC',u22c)-mfun('FresnelC',u12c))-j*(mfun('FresnelS',u22c)-mfun('FresnelS',u12c))).*...
     ((mfun('FresnelC',v22c)-mfun('FresnelC',v12c))-j*((mfun('FresnelS',v22c)-mfun('FresnelS',v12c))))...
     +...
     ((mfun('FresnelC',u23c)-mfun('FresnelC',u13c))-j*(mfun('FresnelS',u23c)-mfun('FresnelS',u13c))).*...
     ((mfun('FresnelC',v23c)-mfun('FresnelC',v13c))-j*((mfun('FresnelS',v23c)-mfun('FresnelS',v13c))))...
     +...
     ((mfun('FresnelC',u24c)-mfun('FresnelC',u14c))-j*(mfun('FresnelS',u24c)-mfun('FresnelS',u14c))).*...
     ((mfun('FresnelC',v24c)-mfun('FresnelC',v14c))-j*((mfun('FresnelS',v24c)-mfun('FresnelS',v14c))))...
     +...
     ((mfun('FresnelC',u25c)-mfun('FresnelC',u15c))-j*(mfun('FresnelS',u25c)-mfun('FresnelS',u15c))).*...
     ((mfun('FresnelC',v25c)-mfun('FresnelC',v15c))-j*((mfun('FresnelS',v25c)-mfun('FresnelS',v15c))))...
     +...
     ((mfun('FresnelC',u26c)-mfun('FresnelC',u16c))-j*(mfun('FresnelS',u26c)-mfun('FresnelS',u16c))).*...
     ((mfun('FresnelC',v26c)-mfun('FresnelC',v16c))-j*((mfun('FresnelS',v26c)-mfun('FresnelS',v16c))))...
     +...
     ((mfun('FresnelC',u27c)-mfun('FresnelC',u17c))-j*(mfun('FresnelS',u27c)-mfun('FresnelS',u17c))).*...
     ((mfun('FresnelC',v27c)-mfun('FresnelC',v17c))-j*((mfun('FresnelS',v27c)-mfun('FresnelS',v17c))))...
     +...
     ((mfun('FresnelC',u28c)-mfun('FresnelC',u18c))-j*(mfun('FresnelS',u28c)-mfun('FresnelS',u18c))).*...
     ((mfun('FresnelC',v28c)-mfun('FresnelC',v18c))-j*((mfun('FresnelS',v28c)-mfun('FresnelS',v18c))))...
     +...
     ((mfun('FresnelC',u29c)-mfun('FresnelC',u19c))-j*(mfun('FresnelS',u29c)-mfun('FresnelS',u19c))).*...
     ((mfun('FresnelC',v29c)-mfun('FresnelC',v19c))-j*((mfun('FresnelS',v29c)-mfun('FresnelS',v19c))))...
     );


figure,plot(xRx,20*log10(abs(Enormalized_c)),'k',BuidlingOutline(:,1),BuidlingOutline(:,2),'k:')
xlabel('Traveled distance (m)')
ylabel('Received field strength relative to direct ray (dB)')
title(['Elevation ' num2str(elev_c) ', Orientation ' num2str(azim_c)])

figure,plot(xRx,abs(Enormalized_c),'k',BuidlingOutline(:,1),BuidlingOutline(:,2),'k:')
xlabel('Traveled distance (m)')
ylabel('Received field strength relative to direct ray (linear units)')
title(['Elevation ' num2str(elev_c) ', Orientation ' num2str(azim_c)])

figure,plot(xRx,abs(Enormalized_c),'k')
xlabel('Traveled distance (m)')
ylabel('Received field strength relative to direct ray (linear units)')
title(['Elevation ' num2str(elev_c) ', Orientation ' num2str(azim_c)])


% Correlation ==============================================================

LagAxis=([1:2*length(Enormalized_a)-1]-length(Enormalized_a)).*stepxRx;

cor_ab=xcorr(abs(Enormalized_a)-mean(abs(Enormalized_a)),abs(Enormalized_b)-mean(abs(Enormalized_b)),'coeff');

figure,plot(LagAxis,cor_ab,'k')
xlabel('Distance lag (m)')
ylabel('Cross-correlation between links a and b')
grid

corrcoef((abs(Enormalized_a)-mean(abs(Enormalized_a))),(abs(Enormalized_b)-mean(abs(Enormalized_b))))

cor_ac=xcorr(abs(Enormalized_a)-mean(abs(Enormalized_a)),abs(Enormalized_c)-mean(abs(Enormalized_c)),'coeff');

figure,plot(LagAxis,cor_ac,'k')
xlabel('Distance lag (m)')
ylabel('Cross-correlation between links a and c')
grid

corrcoef((abs(Enormalized_a)-mean(abs(Enormalized_a))),(abs(Enormalized_c)-mean(abs(Enormalized_c))))

cor_bc=xcorr(abs(Enormalized_b)-mean(abs(Enormalized_b)),abs(Enormalized_c)-mean(abs(Enormalized_c)),'coeff');
figure,plot(LagAxis,cor_bc,'k')
xlabel('Distance lag (m)')
ylabel('Cross-correlation between links b and c')
grid

corrcoef((abs(Enormalized_b)-mean(abs(Enormalized_b))),(abs(Enormalized_c)-mean(abs(Enormalized_c))))


