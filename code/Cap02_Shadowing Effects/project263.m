% INITIALIZE ===========================================================

close all
clear

% Double wedge case.
% INPUT PARAMETERS ==================================================

fMHz=2e3;  % frequency in MHz

% SECONDARY PARAMETERS ==============================================
f=fMHz*1e6;
lambdac=300/fMHz;
kc=2*pi/lambdac;

stepAperture=lambdac;
maxAperture=320;
Na=220;
psi1=100;
psi2=120;


% GEOMETRIC INPUTS =================================================

xt=0;   % m
yt=100; % m   

stepRx = 1;
xr=4000;    % m
yr=[0:stepRx:150]';    % m   SAMPLING AT THE RECEIVER SIDE


% Fresnel reflection coefficient.
ro=-1;

% SAMPLING POINTS ALONG THE APERTURE =================================
xa1=1000; % There are 2 apertures.
xv=2000;
xa2=3000;
ya1=[psi1:stepAperture:maxAperture]; 
ya2=[0:stepAperture:maxAperture]; 
ya3=[psi2:stepAperture:maxAperture];

% Profile.
figure,plot(repmat(xt,1,length(0:yt)),0:yt,'g',xt,yt,'.g'),hold on
yprof=[(xt:xa1-1)*(psi1/(xa1-xt)) ((xa1:xv-1)-xa1)*(-psi1/(xv-xa1))+psi1 ((xv:xa2-1)-xv)*(psi2/(xa2-xv)) ((xa2:xr)-xa2)*(-psi2/(xr-xa2))+psi2];
plot(xt:xr,yprof,'LineWidth',2);
plot(repmat(xr,1,length(yr)),yr,'.g'),hold off
title('Geometry');
xlabel('Distance (m)');
ylabel('Height (m)');
axis([xt-100 xr+100 0 yr(end)+10])


% Window along the aperture.
w1=triang_win(2*length(find(ya1>Na)));
wa1=[ones(1,length((find(ya1<=Na)))) w1(floor(length(w1)/2)+1:end)'];
w2=triang_win(2*length(find(ya2>Na)));
wa2=[ones(1,length((find(ya2<=Na)))) w2(floor(length(w2)/2)+1:end)'];
w3=triang_win(2*length(find(ya3>Na)));
wa3=[ones(1,length((find(ya3<=Na)))) w3(floor(length(w3)/2)+1:end)'];



% Tx-Aperture1 side calculations ==========================================
DistTxApertureX=xa1-xt;

% Direct component.
DistTxApertureY=ya1-yt;
RTxAperture=DistTxApertureX+((DistTxApertureY.^2)/(2*DistTxApertureX));

Efs1_D=exp(-j*kc*RTxAperture)/DistTxApertureX.*wa1;

% Reflected component.
h1TxAperture=yt;
h2TxAperture=ya1-psi1;
RTxApertureR=RTxAperture+((2*h1TxAperture.*h2TxAperture)/DistTxApertureX);

Efs1_R=ro*exp(-j*kc*RTxApertureR)/DistTxApertureX.*wa1;

% Total field in section 1.
Efs1=Efs1_D+Efs1_R;


% Aperture1-valley side calculations ======================================
DistAperture1vX=xv-xa1;

% Direct component.
DistAperture1vY=repmat(ya2',1,length(ya1))-repmat(ya1,length(ya2),1);
RAperture1v=DistAperture1vX+((DistAperture1vY.^2)/(2*DistAperture1vX));

Fd1v=sqrt(kc*(xa1-xt)/(2*pi*j*(xv-xt)*(xv-xa1)));
E1vD=Fd1v*exp(-j*kc*RAperture1v).*repmat(Efs1,length(ya2),1)*stepAperture;
E1vD=sum(E1vD,2).*wa2';

% Reflected component.
h1Aperture1v=ya1-psi1;
h2Aperture1v=ya2;
RAperture1vR=RAperture1v+((2*repmat(h1Aperture1v,length(ya2),1).*repmat(h2Aperture1v',1,length(ya1)))/DistAperture1vX);

E1vR=Fd1v*exp(-j*kc*RAperture1vR).*repmat(Efs1,length(ya2),1)*stepAperture*ro;
E1vR=sum(E1vR,2).*wa2';

% Total field in section 2.
Ed1v=E1vD+E1vR;

% Field plot in the valley due to the obstacle 1 and the reflection.
ind_rx=find(ya2<=yr(end));

% Free space field calculation.
DistTxv=(xv-xt)+(((ya2(ind_rx)'-yt).^2)/(2*(xv-xt)));
Efsv=exp(-j*kc*DistTxv)/(xv-xt);

% Plots
figure,plot(ya2(ind_rx),abs(Ed1v(ind_rx)./Efsv),'k'),grid
xlabel('z(m)');
ylabel('Normalized field strength')
title('After first obstacle')

% figure,plot(ya2,20*log10(abs(Ed1v))),grid

% Valley-Aperture2 side calculations ======================================
DistAperturev2X=xa2-xv;

% Direct component.
DistAperturev2Y=repmat(ya3',1,length(ya2))-repmat(ya2,length(ya3),1);
RAperturev2=DistAperturev2X+((DistAperturev2Y.^2)/(2*DistAperturev2X));

Fdv2=sqrt(kc*(xv-xt)/(2*pi*j*(xa2-xt)*(xa2-xv)));
Ev2=Fdv2*exp(-j*kc*RAperturev2).*repmat(Ed1v.',length(ya3),1)*stepAperture;
Ev2=sum(Ev2,2).*wa3';

% Reflected component.
h1Aperturev2=ya2;
h2Aperturev2=ya3-psi2;
RAperturev2R=RAperturev2+((2*repmat(h1Aperturev2,length(ya3),1).*repmat(h2Aperturev2',1,length(ya2)))/DistAperturev2X);

Ev2R=Fdv2*exp(-j*kc*RAperturev2R).*repmat(Ed1v.',length(ya3),1)*stepAperture*ro;
Ev2R=sum(Ev2R,2).*wa3';

% The field in the second aperture is given by the sum of the direct and
% the reflected components.
Edv2=Ev2+Ev2R;

% Field plot in the valley due to the obstacle 1 and the reflection.
ind_rx=find(ya3<=yr(end));

% Free space field calculation.
DistTx2=(xa2-xt)+(((ya3(ind_rx)-yt).^2)/(2*(xa2-xt)));
Efs2=exp(-j*kc*DistTx2)/(xa2-xt);

% Plots
figure,plot(ya3(ind_rx),abs(Edv2(ind_rx)./Efs2.'),'k'),grid
xlabel('z (m)');
ylabel('Normalized field strength')
title('On top of second obstacle')


% Aperture2-Rx side calculations ==========================================
DistAperture2RxX=xr-xa2;

% Direct component
DistAperture2RxY=repmat(yr,1,length(ya3))-repmat(ya3,length(yr),1);
RAperture2Rx=DistAperture2RxX+((DistAperture2RxY.^2)/(2*DistAperture2RxX));

Fd2rx=sqrt(kc*(xa2-xt)/(2*pi*j*(xr-xt)*(xr-xa2)));
ERx=Fd2rx*exp(-j*kc*RAperture2Rx).*repmat(Edv2.',length(yr),1)*stepAperture;
ERx=sum(ERx,2);

% Reflected component.
h1Aperture2Rx=ya3-psi2;
h2Aperture2Rx=yr;
RAperture2RxR=RAperture2Rx+((2*repmat(h1Aperture2Rx,length(yr),1).*repmat(h2Aperture2Rx,1,length(ya3)))/DistAperture2RxX);

ERxR=Fd2rx*exp(-j*kc*RAperture2RxR).*repmat(Edv2.',length(yr),1)*stepAperture*ro;
ERxR=sum(ERxR,2);

% Total field.
EdRx=ERx+ERxR;

% Free space field calculation.
DistTxRx=(xr-xt)+(((yr-yt).^2)/(2*(xr-xt)));
Efs=exp(-j*kc*DistTxRx)/(xr-xt);

% Plots
figure,plot(yr,abs(EdRx./Efs),'k'),grid
xlabel('z(m)');
ylabel('Normalized field strength')
title('After second Obstacle');