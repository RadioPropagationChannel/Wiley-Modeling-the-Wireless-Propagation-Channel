%
% project25
%
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

xt=-5000;   % m
yt=100;      % m   

stepRx = 0.1;
xr=5000;                % m
yr=[50:stepRx:500]';    % m   SAMPLING AT THE RECEIVER SIDE


% SAMPLING POINTS ALONG THE APERTURE =================================

xa=0;       % beginning of aperture
ya=100;      % beginning of aperture

R1=sqrt(lambdac*abs(xr)*abs(xt)/(abs(xr)+abs(xt))); % Fresnel radius

% 1st contribution Tx-Rx

yy=((yr-yt)./(xr-xt)).*(xa-xt)+yt;   % intersection with screen
hh=ya-yy;                            % obstruction
v=sqrt(2)*hh/R1;                     % normalized obstruction
figure,plot(yr,v,'k', 'LineWidth',2)
title('Normalized obstruction, v, for link Tx-Rx')
xlabel('Rx antenna height (m)')
ylabel('Normalized obstruction, v')

Enormalized1=(1-j)*j/2.*((0.5-mfun('FresnelC',v))-j*(0.5-mfun('FresnelS',v)));

figure,plot(yr,20*log10(abs(Enormalized1)),'k', 'LineWidth',2)
xlabel('Rx antenna height (m)')
title('Tx-Rx contribution. Relative field strength')
ylabel('Relative field strength (dB)')

figure,plot(yr,abs(Enormalized1),'k', 'LineWidth',2)
xlabel('Rx antenna height (m)')
title('Tx-Rx contribution. Relative field strength')
ylabel('Relative field strength (lin. units)')

% 2nd contribution Tx'-Rx

yy=((yr+yt)./(xr-xt)).*(xa-xt)-yt;   % intersection with screen
hh=ya-yy;                            % obstruction
v=sqrt(2)*hh/R1;                     % normalized obstruction
figure,plot(yr,v,'k', 'LineWidth',2)
title('Normalized obstruction, v, for link Image Tx-Rx')
xlabel('Rx antenna height (m)')
ylabel('Normalized obstruction, v')


Enormalized2=(1-j)*j/2.*((0.5-mfun('FresnelC',v))-j*(0.5-mfun('FresnelS',v)));

figure,plot(yr,20*log10(abs(Enormalized2)),'k', 'LineWidth',2)
xlabel('Rx antenna height (m)')
title('ImageTx-Rx contribution. Relative field strength')
ylabel('Relative field strength (dB)')

figure,plot(yr,abs(Enormalized2),'k', 'LineWidth',2)
xlabel('Rx antenna height (m)')
title('ImageTx-Rx contribution. Relative field strength')
ylabel('Relative field strength (lin. units)')

% 3rd contribution Tx-Rx'

yy=((-yr-yt)./(xr-xt)).*(xa-xt)+yt;   % intersection with screen
hh=ya-yy;                             % obstruction

v=sqrt(2)*hh/R1;                      % normalized obstruction
figure,plot(yr,v,'k', 'LineWidth',2)
title('Normalized obstruction, v, for link Tx-Image Rx')
xlabel('Rx antenna height (m)')
ylabel('Normalized obstruction, v')

Enormalized3=(1-j)*j/2.*((0.5-mfun('FresnelC',v))-j*(0.5-mfun('FresnelS',v)));

figure,plot(yr,20*log10(abs(Enormalized3)),'k', 'LineWidth',2)
xlabel('Rx antenna height (m)')
title('Tx-Image Rx contribution. Relative field strength')
ylabel('Relative field strength (dB)')

figure,plot(yr,abs(Enormalized3),'k', 'LineWidth',2)
xlabel('Rx antenna height (m)')
title('Tx-Image Rx contribution. Relative field strength')
ylabel('Relative field strength (lin units)')

% 4th contribution Tx'-Rx'================================================

yy=((-yr+yt)./(xr-xt)).*(xa-xt)-yt;   % intersection with screen
hh=ya-yy;                            % obstruction
v=sqrt(2)*hh/R1;                     % normalized obstruction
figure,plot(yr,v,'k', 'LineWidth',2)
title('Normalized obstruction, v, for link Image Tx-Image Rx')
xlabel('Rx antenna height (m)')
ylabel('Normalized obstruction, v')

Enormalized4=(1-j)*j/2.*((0.5-mfun('FresnelC',v))-j*(0.5-mfun('FresnelS',v)));

figure,plot(yr,20*log10(abs(Enormalized4)),'k', 'LineWidth',2)
xlabel('Rx antenna height (m)')
title('Image Tx-Image Rx contribution. Relative field strength')
ylabel('Relative field strength (dB)')

figure,plot(yr,abs(Enormalized4),'k', 'LineWidth',2)
xlabel('Rx antenna height (m)')
title('Image Tx-Image Rx contribution. Relative field strength')
ylabel('Relative field strength (lin. units)')

% Combine ================================================================

% Phase correction
d_DD=sqrt(((xr-xt)^2)+((yt-yr).^2));
d_RD=sqrt(((xr-xt)^2)+((-yt-yr).^2));
d_DR=sqrt(((xr-xt)^2)+((yt+yr).^2));
d_RR=sqrt(((xr-xt)^2)+((yr-yt).^2));

Phase_Cor2=exp(-j*kc*(d_RD-d_DD));
Phase_Cor3=exp(-j*kc*(d_DR-d_DD));
Phase_Cor4=exp(-j*kc*(d_RR-d_DD));

Enormalized2_p=Enormalized2.*Phase_Cor2;
Enormalized3_p=Enormalized3.*Phase_Cor3;
Enormalized4_p=Enormalized4.*Phase_Cor4;

% ===============================================================

Enormalized=Enormalized1-Enormalized2_p-Enormalized3_p+Enormalized4_p;


figure,plot(yr,20*log10(abs(Enormalized)),'k', 'LineWidth',2)
xlabel('Rx antenna height (m)')
title('Overall relative field strength')
ylabel('Relative field strength (dB)')

figure,plot(yr,abs(Enormalized),'k', 'LineWidth',2)
xlabel('Rx antenna height (m)')
title('Overall relative field strength')
ylabel('Relative field strength (lin. units)')

figure,plot(yr,20*log10(abs(Enormalized)),'k',yr,20*log10(abs(Enormalized1)),'k:',...
    yr,20*log10(abs(Enormalized2)),'k-.',yr,20*log10(abs(Enormalized3)),'k--',...
    yr,20*log10(abs(Enormalized4)),'k:')
xlabel('Rx antenna height (m)')
title('Relative field strength')
ylabel('Relative field strength (dB)')
legend('Overall','Tx-Rx path','ImageTx-Rx path','Tx-ImageRx path','ImageTx-ImageRx path') 
