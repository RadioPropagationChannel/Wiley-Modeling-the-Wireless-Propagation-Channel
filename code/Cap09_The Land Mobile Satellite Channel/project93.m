%
% project93
%
% INITIALIZE ===========================================================

close all
clear

% INPUT PARAMETERS ==================================================

fMHz=2000;       % frequency in MHz
avPower=-15;     % sigma^2  diffuse multipath power rel to LOS in dB
V= 10;           % MS speed in m/s
F= 4;            % sampling fraction of wavelength

% =============== Butterworth filter parameters =========================
Wp=0.1;        
Ws=0.3;
Rp=3; %dB
Rs=40; %dB

% SECONDARY PARAMETERS ==============================================

f=fMHz*1e6;
lambdac=300/fMHz;
kc=2*pi/lambdac;
sigma=sqrt(0.5*10.^(avPower/10));     % Rayleigh parameter


stepxRx=lambdac/F;                % sampling spacing in mobile route 
ts=stepxRx/V;                     % sample time spacing 
fs=1/ts;                          % sampling frequency

% GEOMETRIC INPUTS =================================================

elev_a=30;     % elevation (deg)
azim_a=0;     % azimuth  (deg)

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

xRx=[-20:stepxRx:110];

% Building outline  

BuidlingOutline=[min(xRx),0; x1,0; x1,hB1; x2,hB1; x2,0; ...
    x3,0; x3,hB2; x4,hB2; x4,0; x5,0; x5,hB3; x6,hB3; x6,0; ...
    x7,0; x7,hB4; x8,hB4; x8,0; max(xRx),0];

% Intersection point of ray with screen =============================

x0a=D*tand(azim_a)+xRx;
y0a=D*tand(elev_a)/cosd(azim_a)+hRx;
d2a=(y0a-hRx)/sind(elev_a);

% Fresnel radius

R1a=sqrt(lambdac*abs(d2a));

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



% =========================================================================
% Introducing diffuse multipath

leng=length(xRx);           % length of series of diffuse multipath to be simulated
t_axis=[0:leng-1]*ts;       % create time axis 

ry=rayleigh(sigma,leng);    % FUNCION for generating a Rayleigh series
ry_mod=abs(ry);
t_axis=[0:leng-1]*ts;
figure,plot(t_axis,20*log10(ry_mod),'k')
axis([0 max(t_axis) -50 10])
xlabel('Time (s)'), ylabel('Relative signal level (dB)')
title('Rayleigh fading'),grid

% ===============  Filtered Rayleigh series ========================

[ryfaux, B, A]=filtersignal(ry,Wp,Ws,Rp,Rs);
[H,fre]=freqz(B,A,512,fs);      % For computinf the filter's gain 
% Calculation of filter gain
     [h,T]=impz(B,A);
     h2=h.^2;
     gainFaux=sqrt(cumsum(h2));
     gainF=gainFaux(length(h))
     ryf=ryfaux/gainF;
     ryf=ryf';                   % traspose to adapt to Enormalized_a arrangement
% ........go on .........
ryf_mod=abs(ryf);                        
figure,plot(t_axis*V+xRx(1),20*log10(ryf_mod),'k')
% axis([0 max(t_axis*V) -50 10])
xlabel('Traveled (m)')
ylabel('Reletive signal level (dB)')
title('Filtered Rayleigh fading'), grid

% =======================================================================

% add diffraction modeled series and diffuse scattering

Enormalized_aWithDiffuse=Enormalized_a+ryf;

figure,plot(xRx,20*log10(abs(Enormalized_aWithDiffuse)),'k',BuidlingOutline(:,1),BuidlingOutline(:,2),'k:')
xlabel('Traveled distance (m)')
ylabel('Received field strength relative to direct ray (dB)')
title(['Elevation ' num2str(elev_a) ', Orientation ' num2str(azim_a)])

figure,plot(xRx,abs(Enormalized_aWithDiffuse),'k',BuidlingOutline(:,1),BuidlingOutline(:,2),'k:')
xlabel('Traveled distance (m)')
ylabel('Received field strength relative to direct ray (linear units)')
title(['Elevation ' num2str(elev_a) ', Orientation ' num2str(azim_a)])

figure,plot(xRx,abs(Enormalized_aWithDiffuse),'k')
xlabel('Traveled distance (m)')
ylabel('Received field strength relative to direct ray (linear units)')
title(['Elevation ' num2str(elev_a) ', Orientation ' num2str(azim_a)])

