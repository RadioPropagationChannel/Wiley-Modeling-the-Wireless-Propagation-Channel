%
% project262
%
% INITIALIZE ===========================================================
close all
clear
clc
% INPUT PARAMETERS ==================================================

fMHz=2e3;  % frequency in MHz


% SECONDARY PARAMETERS ==============================================
f=fMHz*1e6;
lambdac=300/fMHz;
kc=2*pi/lambdac;

stepAperture=lambdac/2;
maxAperture=300;
Na=200;
psi1=100;   % Flat-topped obstacle height in meters.


% GEOMETRIC INPUTS =================================================

xt=0;   % m
yt=100;      % m   (0 m indicates same height as screen)

stepRx = 1;
xr=3000;    % m
yr=[0:stepRx:150]';    % m   SAMPLING AT THE RECEIVER SIDE


% Fresnel reflection coefficient.
ro=-1;

% SAMPLING POINTS ALONG THE APERTURE =================================
xa1=1000; % There are 2 apertures.
xa2=2000;
ya=[psi1:stepAperture:maxAperture]; 

% Plot of the geometry.
figure,plot(repmat(xt,1,length(0:yt)),0:yt,'g',xt,yt,'.g'),hold on
plot(repmat(xa1,1,length(0:psi1)),0:psi1,'LineWidth',2);
plot(xa1:xa2,repmat(psi1,1,length(xa1:xa2)),'b','LineWidth',2);
plot(repmat(xa2,1,length(0:psi1)),0:psi1,'b','LineWidth',2);
plot(repmat(xr,1,length(yr)),yr,'.g')
hold off
ylabel('Height (m)')
xlabel('Distance (m)');
axis([xt-100 xr+100 0 yr(end)+10])

% Triangular window along the aperture.
w=triang_win(2*length(find(ya>Na)));
wa=[ones(1,length((find(ya<=Na)))) w(floor(length(w)/2)+1:end)'];

% Tx-Aperture1 side calculations ==========================================
DistTxApertureX=xa1-xt;
DistTxApertureY=ya-yt;
RTxAperture=DistTxApertureX+((DistTxApertureY.^2)/(2*DistTxApertureX));

Efs1=exp(-j*kc*RTxAperture)/DistTxApertureX.*wa;

% Aperture1-Aperture2 side calculations ===================================
% Direct component.
DistAperture12X=xa2-xa1;
DistAperture12Y=repmat(ya',1,length(ya))-repmat(ya,length(ya),1);
RAperture12=DistAperture12X+((DistAperture12Y.^2)/(2*DistAperture12X));

Fd12=sqrt(kc*(xa1-xt)/(2*pi*j*(xa2-xt)*(xa2-xa1)));
E12=Fd12*exp(-j*kc*RAperture12).*repmat(Efs1,length(ya),1)*stepAperture;
E12=sum(E12,2).*wa';

% Reflected component.
hAperture12=ya-psi1;
RAperture12R=RAperture12+(2*repmat(hAperture12,length(ya),1).*repmat(hAperture12',1,length(ya))/DistAperture12X);

E12R=Fd12*exp(-j*kc*RAperture12R).*repmat(Efs1,length(ya),1)*stepAperture*ro;
E12R=sum(E12R,2).*wa';

% The field in the second aperture is given by the sum of the direct and
% the reflected components.
Ed12=E12+E12R;

% Free space field calculation.
ind_rx=find(ya<=yr(end));
DistTx2=(xa2-xt)+(((ya(ind_rx)'-yt).^2)/(2*(xa2-xt)));
Efs2=exp(-j*kc*DistTx2)/(xa2-xt);

% Field relative to free space plot in the obstacle 2 due to the obstacle 1
% and the reflection.
figure,plot(ya(ind_rx),abs(Ed12(ind_rx)./Efs2),'k'),grid
xlabel('z(m)');
ylabel('Normalized field strength')
title('Top of screen # 2')

% Aperture2-Rx side calculations ==========================================
DistAperture2RxX=xr-xa2;
DistAperture2RxY=repmat(yr,1,length(ya))-repmat(ya,length(yr),1);
RAperture2Rx=DistAperture2RxX+((DistAperture2RxY.^2)/(2*DistAperture2RxX));

Fd2rx=sqrt(kc*(xa2-xt)/(2*pi*j*(xr-xt)*(xr-xa2)));
EdRx=Fd2rx*exp(-j*kc*RAperture2Rx).*repmat(Ed12.',length(yr),1)*stepAperture;
EdRx=sum(EdRx,2);

% Free space field calculation.
DistTxRx=(xr-xt)+(((yr-yt).^2)/(2*(xr-xt)));
Efs=exp(-j*kc*DistTxRx)/(xr-xt);

% Plots
figure,plot(yr,abs(EdRx./Efs),'k'),grid
xlabel('z(m)');
ylabel('Normalized field strength')
title('After flat-topped obstacle');