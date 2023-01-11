%
% project83
%
%=======================================================================
clear
close all
clc
%=======================================================================
f=2000           % Frequency in MHz
ptx=1;           % Transmitter power in Watts
gtx=1;           % Transmitter gain in linear units (1 is isotropic) 

%=======================================================================
lambdac=300/f;     % wavelength (m)
kc=2*pi/lambdac;   % propagation constant


%=======================================================================
zt=10;                  % transmitter height
xt=0;

x1=10;                  % one side of street building 1 position
x2=-10;        % other side of steet building 2 position 

zr=1.5;                 % receiver height
xrStep=lambdac;         % MS route along x-axis (m).
xrIni=3;
xr=[x1+xrIni:xrStep:1e4]; 
xr1=-6;

yt=0;         
yr=45;
yStep=lambdac;             
yw=[0:yStep:1e4];       % Wall coordinates
yr1=[1:yStep:1e4];      % Route along the y-axis (m).

y3=30;                  % Aperture y-axis coordinates.
y4=50;      
ya=[y3;y4];

z3=0;                   % Aperture z-axis coordinates.
z4=Inf;
za=[z3;z4];


% Plot of the geometry 
xend=100;
yend=100;
figure,plot(xt,yt,'^k','MarkerFaceColor','k'),hold on
plot(xr(find(xr<=xend)),yr*ones(length(find(xr<=xend)),1),'.k')
plot(xr1*ones(length(find(yr1<=yend)),1),yr1(find(yr1<=yend)),'.k')
plot(x1*ones(length(0:yStep:y3),1),0:yStep:y3,'k','LineWidth',2)
plot(x1*ones(length(y4:yStep:yend),1),y4:yStep:yend,'k','LineWidth',2);
plot(x2*ones(length(0:yStep:yend),1),0:yStep:yend,'k','LineWidth',2);
plot(x1:xrStep:xend,y3*ones(length(x1:xrStep:xend),1),'k','LineWidth',2);
plot(x1:xrStep:xend,y4*ones(length(x1:xrStep:xend),1),'k','LineWidth',2),hold off
title('Geometry');
xlabel('x(m)');
ylabel('y(m)');
axis([2*x2 xend yt-(yend/10) yend])


% Reflection coefficients ===============================================
RG=-1;
RW=-0.5; 

% Route along the y-axis. ================================================
imagezt=-zt;          % for Ground reflected wave
imagext1=2*x1;        % for wall 1 reflected wave (right wall)
imagext2=2*x2;        % for Wall 2 reflected wave (left wall)

dTxRx=sqrt((xr1-xt).^2+(yr1-yt).^2+(zr-zt).^2);   % Direct ray (1)
dimageTxRx=sqrt((xr1-xt).^2+(yr1-yt).^2+(zr-imagezt).^2);  % Ground reflected ray (2)

dimageTx1Rx=sqrt((xr1-imagext1).^2+(yr1-yt).^2+(zr-zt).^2); % building 1 reflected ray (3)
dimageTx2Rx=sqrt((xr1-imagext2).^2+(yr1-yt).^2+(zr-zt).^2); % building 2 reflected ray (4)

E0=exp(-j*kc*dTxRx)./dTxRx;  % direct ray, free space field strength in dB

EG=RG*exp(-j*kc*dimageTxRx)./dimageTxRx;
EW1=RW*exp(-j*kc*dimageTx1Rx)./dimageTx1Rx;
EW2=RW*exp(-j*kc*dimageTx2Rx)./dimageTx2Rx;

Et=(E0+EG+EW1+EW1);  % E field 4 ray model V/m.

FdB=20*log10(abs(Et./E0));
figure,plot(yr1,FdB,'k','LineWidth',2);
 xlabel('Traveled distance(m)');
title('Rx signal(dB). Route along the y-axis');


% Free-space losses.
Ly=32.5+20*log10(f)+20*log10(yr1*1e-3);

figure,plot(yr1,-Ly+FdB,'k','LineWidth',2);
 xlabel('Traveled distance(m)');
title('Rx signal with free space losses (dB). Route along the y-axis')

% Route along the aperture (x-axis).
% Fresnel-Kirchhoff diffraction parameter.
Fe=x1*(xr-x1)./xr;              % Focal distance (x-axis).
alpha=sqrt(2./(lambdac*Fe));


%=====================================================================
% Imaginary sources coordinates.
imagezr=-zr;               % for Ground reflected wave
imageyr3=2*y3-yr;           % for wall 3.
imageyr4=2*y4-yr;           % for wall 4.

% Path distances for each case.
DistTxRx=sqrt((xr-xt).^2+(yr-yt).^2+(zr-zt).^2);                % Direct ray distance.
DistTxRxG=sqrt((xr-xt).^2+(yr-yt).^2+(imagezr-zt).^2);          % Reflected ray in the ground of the alley.
DistTxRx3=sqrt((xr-xt).^2+(imageyr3-yt).^2+(zr-zt).^2);         % Reflected ray in the wall 3.
DistTxRx4=sqrt((xr-xt).^2+(imageyr4-yt).^2+(zr-zt).^2);         % Reflected ray in the wall 4.

DistTxGRx=sqrt((xr-xt).^2+(yr-yt).^2+(zr-imagezt).^2);          % Reflected ray in the ground of the main street.
DistTxGRxG=sqrt((xr-xt).^2+(yr-yt).^2+(imagezr-imagezt).^2);    % Reflected ray in the ground of the main street and alley.
DistTxGRx3=sqrt((xr-xt).^2+(imageyr3-yt).^2+(zr-imagezt).^2);   % Reflected ray in the ground of the main street and in the wall 3.
DistTxGRx4=sqrt((xr-xt).^2+(imageyr4-yt).^2+(zr-imagezt).^2);   % Reflected ray in the ground of the main street and in the wall 4.

DistTx2Rx=sqrt((xr-imagext2).^2+(yr-yt).^2+(zr-zt).^2);          % Reflected ray in the wall 2 .
DistTx2RxG=sqrt((xr-imagext2).^2+(yr-yt).^2+(imagezr-zt).^2);    % Reflected ray in the wall 2 and in the ground of the alley.
DistTx2Rx3=sqrt((xr-imagext2).^2+(imageyr3-yt).^2+(zr-zt).^2);   % Reflected ray in the wall 2 and in the wall 3.
DistTx2Rx4=sqrt((xr-imagext2).^2+(imageyr4-yt).^2+(zr-zt).^2);   % Reflected ray in the wall 2 and in the wall 4.


% Intersection point of the tx-rx path into the screen.
OTxRx=repmat([xt;yt;zt],1,length(xr))+x1.*[ones(1,length(xr));(yt-yr)./(xt-xr);(zt-zr)./(xt-xr)];                   % Direct ray.
OTxRxG=repmat([xt;yt;zt],1,length(xr))+x1.*[ones(1,length(xr));(yt-yr)./(xt-xr);(zt-imagezr)./(xt-xr)];             % Alley ground reflected ray.
OTxRx3=repmat([xt;yt;zt],1,length(xr))+x1.*[ones(1,length(xr));(yt-imageyr3)./(xt-xr);(zt-zr)./(xt-xr)];            % Alley wall 3 reflected ray.
OTxRx4=repmat([xt;yt;zt],1,length(xr))+x1.*[ones(1,length(xr));(yt-imageyr4)./(xt-xr);(zt-zr)./(xt-xr)];            % Alley wall 4 reflected ray.

OTxGRx=repmat([xt;yt;imagezt],1,length(xr))+x1.*[ones(1,length(xr));(yt-yr)./(xt-xr);(imagezt-zr)./(xt-xr)];        % Main street ground reflected ray.
OTxGRxG=repmat([xt;yt;imagezt],1,length(xr))+x1.*[ones(1,length(xr));(yt-yr)./(xt-xr);(imagezt-imagezr)./(xt-xr)];  % Main street and alley ground reflected ray.
OTxGRx3=repmat([xt;yt;imagezt],1,length(xr))+x1.*[ones(1,length(xr));(yt-imageyr3)./(xt-xr);(imagezt-zr)./(xt-xr)]; % Main street ground and alley wall 3 reflected ray.
OTxGRx4=repmat([xt;yt;imagezt],1,length(xr))+x1.*[ones(1,length(xr));(yt-imageyr4)./(xt-xr);(imagezt-zr)./(xt-xr)]; % Main street ground and alley wall 4 reflected ray.

OTx2Rx=repmat([imagext2;yt;zt],1,length(xr))-(imagext2-x1).*[ones(1,length(xr));(yt-yr)./(imagext2-xr);(zt-zr)./(imagext2-xr)];         % Main street wall 2 reflected ray.
OTx2RxG=repmat([imagext2;yt;zt],1,length(xr))-(imagext2-x1).*[ones(1,length(xr));(yt-yr)./(imagext2-xr);(zt-imagezr)./(imagext2-xr)];   % Main street wall 2 and alley ground reflected ray.
OTx2Rx3=repmat([imagext2;yt;zt],1,length(xr))-(imagext2-x1).*[ones(1,length(xr));(yt-imageyr3)./(imagext2-xr);(zt-zr)./(imagext2-xr)];  % Main street wall 2 and wall 3 reflected ray.
OTx2Rx4=repmat([imagext2;yt;zt],1,length(xr))-(imagext2-x1).*[ones(1,length(xr));(yt-imageyr4)./(imagext2-xr);(zt-zr)./(imagext2-xr)];  % Main street wall 2 and wall 4 reflected ray.


% Aperture coordinates referred to the intersection point.
yaOTxRx=repmat(ya,1,length(xr))-repmat(OTxRx(2,:),length(ya),1);
zaOTxRx=z3-OTxRx(3,:);

yaOTxRxG=repmat(ya,1,length(xr))-repmat(OTxRxG(2,:),length(ya),1);
zaOTxRxG=z3-OTxRxG(3,:);

yaOTxRx3=repmat(ya,1,length(xr))-repmat(OTxRx3(2,:),length(ya),1);
zaOTxRx3=z3-OTxRx3(3,:);

yaOTxRx4=repmat(ya,1,length(xr))-repmat(OTxRx4(2,:),length(ya),1);
zaOTxRx4=z3-OTxRx4(3,:);

yaOTxGRx=repmat(ya,1,length(xr))-repmat(OTxGRx(2,:),length(ya),1);
zaOTxGRx=z3-OTxGRx(3,:);

yaOTxGRxG=repmat(ya,1,length(xr))-repmat(OTxGRxG(2,:),length(ya),1);
zaOTxGRxG=z3-OTxGRxG(3,:);

yaOTxGRx3=repmat(ya,1,length(xr))-repmat(OTxGRx3(2,:),length(ya),1);
zaOTxGRx3=z3-OTxGRx3(3,:);

yaOTxGRx4=repmat(ya,1,length(xr))-repmat(OTxGRx4(2,:),length(ya),1);
zaOTxGRx4=z3-OTxGRx4(3,:);

yaOTx2Rx=repmat(ya,1,length(xr))-repmat(OTx2Rx(2,:),length(ya),1);
zaOTx2Rx=z3-OTx2Rx(3,:);

yaOTx2RxG=repmat(ya,1,length(xr))-repmat(OTx2RxG(2,:),length(ya),1);
zaOTx2RxG=z3-OTx2RxG(3,:);

yaOTx2Rx3=repmat(ya,1,length(xr))-repmat(OTx2Rx3(2,:),length(ya),1);
zaOTx2Rx3=z3-OTx2Rx3(3,:);

yaOTx2Rx4=repmat(ya,1,length(xr))-repmat(OTx2Rx4(2,:),length(ya),1);
zaOTx2Rx4=z3-OTx2Rx4(3,:);


% Integration variables u (y-axis) and v (z-axis).
uTxRx=repmat(alpha,length(ya),1).*yaOTxRx;
vTxRx=alpha.*zaOTxRx;

uTxRxG=repmat(alpha,length(ya),1).*yaOTxRxG;
vTxRxG=alpha.*zaOTxRxG;

uTxRx3=repmat(alpha,length(ya),1).*yaOTxRx3;
vTxRx3=alpha.*zaOTxRx3;

uTxRx4=repmat(alpha,length(ya),1).*yaOTxRx4;
vTxRx4=alpha.*zaOTxRx4;

uTxGRx=repmat(alpha,length(ya),1).*yaOTxGRx;
vTxGRx=alpha.*zaOTxGRx;

uTxGRxG=repmat(alpha,length(ya),1).*yaOTxGRxG;
vTxGRxG=alpha.*zaOTxGRxG;

uTxGRx3=repmat(alpha,length(ya),1).*yaOTxGRx3;
vTxGRx3=alpha.*zaOTxGRx3;

uTxGRx4=repmat(alpha,length(ya),1).*yaOTxGRx4;
vTxGRx4=alpha.*zaOTxGRx4;

uTx2Rx=repmat(alpha,length(ya),1).*yaOTx2Rx;
vTx2Rx=alpha.*zaOTx2Rx;

uTx2RxG=repmat(alpha,length(ya),1).*yaOTx2RxG;
vTx2RxG=alpha.*zaOTx2RxG;

uTx2Rx3=repmat(alpha,length(ya),1).*yaOTx2Rx3;
vTx2Rx3=alpha.*zaOTx2Rx3;

uTx2Rx4=repmat(alpha,length(ya),1).*yaOTx2Rx4;
vTx2Rx4=alpha.*zaOTx2Rx4;

% Fresnel integrals for each axis.
[CuTxRx,SuTxRx]=Fresnel_integrals(uTxRx);
[CvTxRx,SvTxRx]=Fresnel_integrals(vTxRx);

[CuTxRxG,SuTxRxG]=Fresnel_integrals(uTxRxG);
[CvTxRxG,SvTxRxG]=Fresnel_integrals(vTxRxG);

[CuTxRx3,SuTxRx3]=Fresnel_integrals(uTxRx3);
[CvTxRx3,SvTxRx3]=Fresnel_integrals(vTxRx3);

[CuTxRx4,SuTxRx4]=Fresnel_integrals(uTxRx4);
[CvTxRx4,SvTxRx4]=Fresnel_integrals(vTxRx4);

[CuTxGRx,SuTxGRx]=Fresnel_integrals(uTxGRx);
[CvTxGRx,SvTxGRx]=Fresnel_integrals(vTxGRx);

[CuTxGRxG,SuTxGRxG]=Fresnel_integrals(uTxGRxG);
[CvTxGRxG,SvTxGRxG]=Fresnel_integrals(vTxGRxG);

[CuTxGRx3,SuTxGRx3]=Fresnel_integrals(uTxGRx3);
[CvTxGRx3,SvTxGRx3]=Fresnel_integrals(vTxGRx3);

[CuTxGRx4,SuTxGRx4]=Fresnel_integrals(uTxGRx4);
[CvTxGRx4,SvTxGRx4]=Fresnel_integrals(vTxGRx4);

[CuTx2Rx,SuTx2Rx]=Fresnel_integrals(uTx2Rx);
[CvTx2Rx,SvTx2Rx]=Fresnel_integrals(vTx2Rx);

[CuTx2RxG,SuTx2RxG]=Fresnel_integrals(uTx2RxG);
[CvTx2RxG,SvTx2RxG]=Fresnel_integrals(vTx2RxG);

[CuTx2Rx3,SuTx2Rx3]=Fresnel_integrals(uTx2Rx3);
[CvTx2Rx3,SvTx2Rx3]=Fresnel_integrals(vTx2Rx3);

[CuTx2Rx4,SuTx2Rx4]=Fresnel_integrals(uTx2Rx4);
[CvTx2Rx4,SvTx2Rx4]=Fresnel_integrals(vTx2Rx4);

% Aperture integral.
BTxRx=([CuTxRx(2,:)-CuTxRx(1,:)]-j*[SuTxRx(2,:)-SuTxRx(1,:)]).*([0.5-CvTxRx]-j*[0.5-SvTxRx]);
BTxRxG=([CuTxRxG(2,:)-CuTxRxG(1,:)]-j*[SuTxRxG(2,:)-SuTxRxG(1,:)]).*([0.5-CvTxRxG]-j*[0.5-SvTxRxG]);
BTxRx3=([CuTxRx3(2,:)-CuTxRx3(1,:)]-j*[SuTxRx3(2,:)-SuTxRx3(1,:)]).*([0.5-CvTxRx3]-j*[0.5-SvTxRx3]);
BTxRx4=([CuTxRx4(2,:)-CuTxRx4(1,:)]-j*[SuTxRx4(2,:)-SuTxRx4(1,:)]).*([0.5-CvTxRx4]-j*[0.5-SvTxRx4]);

BTxGRx=([CuTxGRx(2,:)-CuTxGRx(1,:)]-j*[SuTxGRx(2,:)-SuTxGRx(1,:)]).*([0.5-CvTxGRx]-j*[0.5-SvTxGRx]);
BTxGRxG=([CuTxGRxG(2,:)-CuTxGRxG(1,:)]-j*[SuTxGRxG(2,:)-SuTxGRxG(1,:)]).*([0.5-CvTxGRxG]-j*[0.5-SvTxGRxG]);
BTxGRx3=([CuTxGRx3(2,:)-CuTxGRx3(1,:)]-j*[SuTxGRx3(2,:)-SuTxGRx3(1,:)]).*([0.5-CvTxGRx3]-j*[0.5-SvTxGRx3]);
BTxGRx4=([CuTxGRx4(2,:)-CuTxGRx4(1,:)]-j*[SuTxGRx4(2,:)-SuTxGRx4(1,:)]).*([0.5-CvTxGRx4]-j*[0.5-SvTxGRx4]);

BTx2Rx=([CuTx2Rx(2,:)-CuTx2Rx(1,:)]-j*[SuTx2Rx(2,:)-SuTx2Rx(1,:)]).*([0.5-CvTx2Rx]-j*[0.5-SvTx2Rx]);
BTx2RxG=([CuTx2RxG(2,:)-CuTx2RxG(1,:)]-j*[SuTx2RxG(2,:)-SuTx2RxG(1,:)]).*([0.5-CvTx2RxG]-j*[0.5-SvTx2RxG]);
BTx2Rx3=([CuTx2Rx3(2,:)-CuTx2Rx3(1,:)]-j*[SuTx2Rx3(2,:)-SuTx2Rx3(1,:)]).*([0.5-CvTx2Rx3]-j*[0.5-SvTx2Rx3]);
BTx2Rx4=([CuTx2Rx4(2,:)-CuTx2Rx4(1,:)]-j*[SuTx2Rx4(2,:)-SuTx2Rx4(1,:)]).*([0.5-CvTx2Rx4]-j*[0.5-SvTx2Rx4]);


% Phase correction (The Fresnel integrals don't consider the phase
% differences between the direct signal ray and the reflected one).
pcTxRxG=exp(-j*kc*(DistTxRxG-DistTxRx));
pcTxRx3=exp(-j*kc*(DistTxRx3-DistTxRx));
pcTxRx4=exp(-j*kc*(DistTxRx4-DistTxRx));

pcTxGRx=exp(-j*kc*(DistTxGRx-DistTxRx));
pcTxGRxG=exp(-j*kc*(DistTxGRxG-DistTxRx));
pcTxGRx3=exp(-j*kc*(DistTxGRx3-DistTxRx));
pcTxGRx4=exp(-j*kc*(DistTxGRx4-DistTxRx));

pcTx2Rx=exp(-j*kc*(DistTx2Rx-DistTxRx));
pcTx2RxG=exp(-j*kc*(DistTx2RxG-DistTxRx));
pcTx2Rx3=exp(-j*kc*(DistTx2Rx3-DistTxRx));
pcTx2Rx4=exp(-j*kc*(DistTx2Rx4-DistTxRx));


 % Diffraction factor (in complex numbers).
 FTxRx=(j/2)*BTxRx;
 FTxRxG=(j/2)*BTxRxG*RG.*pcTxRxG;
 FTxRx3=(j/2)*BTxRx3*RW.*pcTxRx3;
 FTxRx4=(j/2)*BTxRx4*RW.*pcTxRx4;
 
 FTxGRx=(j/2)*BTxGRx*RG.*pcTxGRx;
 FTxGRxG=(j/2)*BTxGRxG*(RG^2).*pcTxGRxG;
 FTxGRx3=(j/2)*BTxGRx3*RG*RW.*pcTxGRx3;
 FTxGRx4=(j/2)*BTxGRx4*RG*RW.*pcTxGRx4;

 FTx2Rx=(j/2)*BTx2Rx*RW.*pcTx2Rx;
 FTx2RxG=(j/2)*BTx2RxG*RW*RG.*pcTx2RxG;
 FTx2Rx3=(j/2)*BTx2Rx3*(RW^2).*pcTx2Rx3;
 FTx2Rx4=(j/2)*BTx2Rx4*(RW^2).*pcTx2Rx4;
 
 % Total Diffraction factor.
 FT=FTxRx+FTxRxG+FTxRx3+FTxRx4+FTxGRx+FTxGRxG+FTxGRx3+FTxGRx4+...
     FTx2Rx+FTx2RxG+FTx2Rx3+FTx2Rx4;
 
 
 % Total Diffraction factor plot.
 figure,plot(xr,20*log10(abs(FT)),'k','LineWidth',2);
  xlabel('Traveled distance(m)');
 title('Rx signal (dB). Route along x-axis');
 ylabel('Inverse of excess loss, -Lexcess (dB)')

 
 % Free-space losses.
 Lx=32.5+20*log10(f)+20*log10(xr*1e-3);

 figure,plot(xr,-Lx+20*log10(abs(FT)),'k','LineWidth',2);
 xlabel('Traveled distance(m)');
 ylabel('Inverse of path loss, -L (dB)')
 title('Rx signal with free space losses (dB). Route along x-axis');
 
 
 % Rx signal (the 2 routes).
 figure,plot(yr1,FdB,':k',[xr1:xr(1)-1 xr]-xr1,[NaN(1,length(xr1:xr(1)-1)) 20*log10(abs(FT))],'k','LineWidth',2);
 xlabel('Traveled distance(m)');
 title('Rx signal (dB)');
 legend('y-axis route','x-axis route');
  ylabel('Relative signal level (dB)');
 
 figure,plot(yr1,-Ly+FdB,':k',[xr1:xr(1)-1 xr]-xr1,[NaN(1,length(xr1:xr(1)-1)) -Lx+20*log10(abs(FT))],'k','LineWidth',2);
 xlabel('Traveled distance(m)');
 title('Rx signal (dB)');
 legend('y-axis route','x-axis route');
 ylabel('Inverse of path loss, -L (dB)')