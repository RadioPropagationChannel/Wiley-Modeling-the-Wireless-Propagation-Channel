%
% project21
%
% INITIALIZE ===========================================================
close all
clear
clc
% INPUT PARAMETERS ==================================================

fMHz=100;  % frequency in MHz


% SECONDARY PARAMETERS ==============================================

f=fMHz*1e6;
lambdac=300/fMHz;
kc=2*pi/lambdac;

stepAperture=lambdac/8;
maxAperture=300;      % end of window (m)
Na=200;               % Beguinning of linear decaying part of window (m)
psi=100;              % edge heigth (m)


% GEOMETRIC INPUTS =================================================

xt=-1000;   % m
yt=100;     % m   (0 m indicates same height as screen)

stepRx = 1;
xr=1000;                % m
yr=[0:stepRx:250]';    % m   SAMPLING AT THE RECEIVER SIDE

DistTxRx=(xr-xt)+((yr-yt).^2)/(2*(xr-xt));

% SAMPLING POINTS ALONG THE APERTURE =================================
xa=0;
ya=[psi:stepAperture:maxAperture];

% Window along the aperture.
w=triang_win(2*length(find(ya>Na)));
wa=[ones(1,length((find(ya<=Na)))) w(floor(length(w)/2)+1:end)'];

% Tx-Aperture side calculations ==========================================
DistTxApertureX=xa-xt;
DistTxApertureY=abs(ya-yt);
RTxAperture=DistTxApertureX+((DistTxApertureY.^2)/(2*DistTxApertureX));

Efs1=exp(-j*kc*RTxAperture)/DistTxApertureX.*wa;

% Aperture-Rx side calculations ==========================================
DistApertureRxX=xr-xa;
DistApertureRxY=abs(repmat(yr,1,length(ya))-repmat(ya,length(yr),1));
RApertureRx=DistApertureRxX+((DistApertureRxY.^2)/(2*DistApertureRxX));

Fd=sqrt(kc*(xa-xt)/(2*pi*j*(xr-xt)*xr));
Ed=Fd*exp(-j*kc*RApertureRx).*repmat(Efs1,length(yr),1)*stepAperture;
Ed=sum(Ed,2);

Efs2=exp(-j*kc*DistTxRx)/(xr-xt);

figure,plot(yr,abs(Ed./Efs2),'k')
xlabel('z(m)');
ylabel('Normalized field strength')
grid

% R1=sqrt(lambdac*abs(xr)*abs(xt)/(abs(xr)+abs(xt)));
% xx=0;
% yy=((yr-yt)./(xr-xt)).*(xx-xt)+yt;
% hh=-yy;
% v=sqrt(2)*hh/R1;
% figure,plot(v,abs(Ed./Efs2))
% ylabel('dB');
% xlabel('v');