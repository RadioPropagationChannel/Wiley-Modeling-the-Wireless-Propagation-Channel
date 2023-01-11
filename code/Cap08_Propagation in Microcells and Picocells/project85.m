% project85

% INITIALIZE ===========================================================
clc
close all
clear

% INPUT PARAMETERS ==================================================

fMHz=2000;  % frequency in MHz


% SECONDARY PARAMETERS ==============================================

f=fMHz*1e6;
lambdac=300/fMHz;
kc=2*pi/lambdac;


% GEOMETRIC INPUTS =================================================

xt=100;   % m
yt=50;      % m   (0 m indicates same height as screen)
zt=50

stepRx = lambdac/8;
xr=-10;                % m
yr=[-10:stepRx:10]';    % m   SAMPLING AT THE RECEIVER SIDE
zr=0;

% WINDOW =================================

ya1=-0.5;      % beginning of aperture
ya2=0.5;
za1=-0.4
za2=0.4;


% Point O

y0=-xt*(yr-yt)./(xr-xt)+yt;
z0=-xt*(zr-zt)./(xr-xt)+zt

% distances

d1=sqrt(xt^2+(yt-y0).^2+(zt-z0).^2);
d2=sqrt(xr^2+(yr-y0).^2+(zr-z0).^2);

% Radius of 1st Fresnel zone at window

R1=sqrt(lambdac.*d1.*d2./(d1+d2)); % Radious of 1st Fresnel zone

% diffraction parameters

u1=sqrt(2)*(ya1-y0)./R1;
u2=sqrt(2)*(ya2-y0)./R1;

v1=sqrt(2)*(za1-z0)./R1;
v2=sqrt(2)*(za2-z0)./R1;

Enormalized=j/2.*((mfun('FresnelC',u2)-mfun('FresnelC',u1))-j*(mfun('FresnelS',u2)-mfun('FresnelS',u1))).*...
    ((mfun('FresnelC',v2)-mfun('FresnelC',v1))-j*(mfun('FresnelS',v2)-mfun('FresnelS',v1)));

figure,plot(yr,20*log10(abs(Enormalized)),'k','LineWidth',2)
xlabel('Receiver position (m)')
ylabel('Relative field strength (dB)')
title('Building entry through window')
grid


