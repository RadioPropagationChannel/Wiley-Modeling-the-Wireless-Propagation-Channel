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

xt=-1000;   % m
yt=0;      % m   (0 m indicates same height as screen)

stepRx = 1;
xr=1000;                % m
yr=[-50:stepRx:50]';    % m   SAMPLING AT THE RECEIVER SIDE


% SAMPLING POINTS ALONG THE APERTURE =================================

xa=0;      % beginning of aperture
ya=0;      % beginning of aperture

R1=sqrt(lambdac*abs(xr)*abs(xt)/(abs(xr)+abs(xt))); % Fresnel radius


yy=((yr-yt)./(xr-xt)).*(xa-xt)+yt;   % intersection with screen
hh=ya-yy;                            % obstruction
v=sqrt(2)*hh/R1;                     % normalized obstruction

Enormalized=(1-j)*j/2.*((0.5-mfun('FresnelC',v))-j*(0.5-mfun('FresnelS',v)));

figure,plot(hh,20*log10(abs(Enormalized)),'k','LineWidth',2)
xlabel('Obstruction, h (m)')
ylabel('Field strength relative to the direct signal (dB)')
title('Knife-edge model')
grid

figure,plot(v,20*log10(abs(Enormalized)),'k','LineWidth',2)
xlabel('Normalized obstruction, v (m)')
ylabel('Field strength relative to the direct signal (dB)')
title('Knife-edge model')
grid

figure,plot(yr,20*log10(abs(Enormalized)),'k','LineWidth',2)
xlabel('Receive antenna height (m)')
ylabel('Field strength relative to the direct signal (dB)')
title('Knife-edge model')
grid


