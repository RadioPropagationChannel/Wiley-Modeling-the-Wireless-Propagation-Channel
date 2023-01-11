%
% project82
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
zt=10;        % transmitter height
xt=0;
zr=1.5        % receiver height
xr=-6; 

yt=0;         % MS route along y-axis
yrStep=0.5;
yr=[1:yrStep:100000];

x1=10;        % one side of street buidling 1 position
x2=-10;        % other side of steet building 2 position 

% Reflection coefficients ===============================================
RG=-1;
RW=-0.5;

%=====================================================================
imagezt=-zt;          % for Ground reflected wave
imagext1=2*x1;        % for Wall 1 reflected wave
imagext2=2*x2;        % for Wall 2 reflected wave
imagext12=-(2*abs(x1)+2*abs(x2));  % for Wall 1 and Wall 2 wave
imagext21=2*abs(x1)+2*abs(x2);     % for Wall 2 and Wall 1 wave

dTxRx=sqrt((xr-xt).^2+(yr-yt).^2+(zr-zt).^2);   % Direct ray (1)
dimageTxRx=sqrt((xr-xt).^2+(yr-yt).^2+(zr-imagezt).^2);  % Ground reflected ray (2)

dimageTx1Rx=sqrt((xr-imagext1).^2+(yr-yt).^2+(zr-zt).^2); % building 1 reflected ray (3)
dimageTx2Rx=sqrt((xr-imagext2).^2+(yr-yt).^2+(zr-zt).^2); % buildign2 reflected ray (4)

dimageTx12Rx=sqrt((xr-imagext12).^2+(yr-yt).^2+(zr-zt).^2); % building 1 and 2 reflected ray (5)
dimageTx21Rx=sqrt((xr-imagext21).^2+(yr-yt).^2+(zr-zt).^2); % building 2 and 1 reflected ray (6)

%====================================================================
% figure,plot(yr,dTxRx,'r',yr,dimageTxRx,'g',yr,dimageTx1Rx,'k',yr,dimageTx2Rx,'y')
% 
% figure,plot(yr,dTxRx-dimageTxRx,'g',yr,dTxRx-dimageTx1Rx,'r',yr,dTxRx-dimageTx2Rx,'b')


%========================================================================

E0=sqrt(60*ptx*gtx)*exp(-j*kc*dTxRx)./dTxRx;  % direct ray, free space field strength in V/m

% E field 2 ray model (direct ray and ground reflection) V/m
E2=sqrt(60*ptx*gtx)*(exp(-j*kc*dTxRx)./dTxRx+...
    RG*exp(-j*kc*dimageTxRx)./dimageTxRx);     

% E field 4 ray model (direct ray, ground, wall 1 and wall 2 reflections) V/m
E4=sqrt(60*ptx*gtx)*(exp(-j*kc*dTxRx)./dTxRx+...
    RG*exp(-j*kc*dimageTxRx)./dimageTxRx...
    +RW*exp(-j*kc*dimageTx1Rx)./dimageTx1Rx...
    +RW*exp(-j*kc*dimageTx2Rx)./dimageTx2Rx);

% E field 6 ray model (direct ray,wall 1 and 2, wall 2 and 1, ground, wall 1 
% and wall 2 reflections) V/m
E6=sqrt(60*ptx*gtx)*(exp(-j*kc*dTxRx)./dTxRx+...
    RG*exp(-j*kc*dimageTxRx)./dimageTxRx...
    +RW*exp(-j*kc*dimageTx1Rx)./dimageTx1Rx...
    +RW*exp(-j*kc*dimageTx2Rx)./dimageTx2Rx...
    +(RW.^2)*exp(-j*kc*dimageTx12Rx)./dimageTx12Rx...
    +(RW.^2)*exp(-j*kc*dimageTx21Rx)./dimageTx21Rx);

% field strenghts are in V/m, now they are converted to dBuV/m
figure,semilogx(yr,20*log10(abs(E0)*1e6),'k:',yr,20*log10(abs(E2)*1e6),'k--',yr,20*log10(abs(E4)*1e6),'k','LineWidth',2)
xlabel('Distance from Tx (m)')
ylabel('Received field strength (dB \muV/m)')
title('Recived field: free space, 2-ray model and 4-ray model')
legend('free space','2-ray model','4-ray model')


figure,semilogx(yr,20*log10(abs(E0)*1e6),'k:',yr,20*log10(abs(E4)*1e6),'k--',yr,20*log10(abs(E6)*1e6),'k','LineWidth',2)
xlabel('Distance from Tx (m)')
ylabel('Received field strength (dB \muV/m)')
title('Recived field: free space, 4-ray model and 6-ray model')
legend('free space','4-ray model','6-ray model')


figure,semilogx(yr,20*log10(abs(E2./E0)),'k--',yr,20*log10(abs(E4./E0)),'k','LineWidth',2);
xlabel('Distance from Tx (m)')
ylabel('Relative received field strength (dB)')
title('Relative recived field for 2-ray model and 4-ray model')
legend('2-ray model','4-ray model')


figure,semilogx(yr,20*log10(abs(E4./E0)),'k--',yr,20*log10(abs(E6./E0)),'k','LineWidth',2);
xlabel('Distance from Tx (m)')
ylabel('Relative received field strength (dB)')
title('Relative recived field for 4-ray model and 6-ray model')
legend('4-ray model','6-ray model')

% Power Delay Profile (PDP).
c=3e8;
yr=100; % (m) Static receiver

dTxRx=sqrt((xr-xt).^2+(yr-yt).^2+(zr-zt).^2);   % Direct ray (1)
dimageTxRx=sqrt((xr-xt).^2+(yr-yt).^2+(zr-imagezt).^2);  % Ground reflected ray (2)
dimageTx1Rx=sqrt((xr-imagext1).^2+(yr-yt).^2+(zr-zt).^2); % building 1 reflected ray (3)
dimageTx2Rx=sqrt((xr-imagext2).^2+(yr-yt).^2+(zr-zt).^2); % buildign2 reflected ray (4)
dimageTx12Rx=sqrt((xr-imagext12).^2+(yr-yt).^2+(zr-zt).^2); % building 1 and 2 reflected ray (5)
dimageTx21Rx=sqrt((xr-imagext21).^2+(yr-yt).^2+(zr-zt).^2); % building 2 and 1 reflected ray (6)

delayImageTxRx=(dimageTxRx-dTxRx)/c;
delayImageTx1Rx=(dimageTx1Rx-dTxRx)/c;
delayImageTx2Rx=(dimageTx2Rx-dTxRx)/c;
delayImageTx12Rx=(dimageTx12Rx-dTxRx)/c;
delayImageTx21Rx=(dimageTx21Rx-dTxRx)/c;

E0=sqrt(60*ptx*gtx)*exp(-j*kc*dTxRx)./dTxRx;  % direct ray, free space field strength in V/m.
EimageTxRx=sqrt(60*ptx*gtx)*(RG*exp(-j*kc*dimageTxRx)./dimageTxRx);  % Ground reflected ray field strength in V/m.
EimageTx1Rx=sqrt(60*ptx*gtx)*(RW*exp(-j*kc*dimageTx1Rx)./dimageTx1Rx); % building 1 reflected ray field strength in V/m.
EimageTx2Rx=sqrt(60*ptx*gtx)*(RW*exp(-j*kc*dimageTx2Rx)./dimageTx2Rx); % building 2 reflected ray field strength in V/m.
EimageTx12Rx=sqrt(60*ptx*gtx)*((RW.^2)*exp(-j*kc*dimageTx12Rx)./dimageTx12Rx); % building 1 and 2 reflected ray field strength in V/m.
EimageTx21Rx=sqrt(60*ptx*gtx)*((RW.^2)*exp(-j*kc*dimageTx21Rx)./dimageTx21Rx); % building 2 and 1 reflected ray field strength in V/m.

PimageTxRx=20*log10(abs(EimageTxRx./E0));        % Ground reflected ray power.
PimageTx1Rx=20*log10(abs(EimageTx1Rx./E0));      % Building 1 reflected ray power.
PimageTx2Rx=20*log10(abs(EimageTx2Rx./E0));      % Building 2 reflected ray power.
PimageTx12Rx=20*log10(abs(EimageTx12Rx./E0));    % Building 1 and 2 reflected ray power.
PimageTx21Rx=20*log10(abs(EimageTx21Rx./E0));    % Building 2 and 1 reflected ray power.

% 2-rays model PDP.
delay2=[0 delayImageTxRx];
power2=[0 PimageTxRx];

figure,stem2D(delay2,power2,floor(power2(end)))
xlabel('Excess delay (s)')
ylabel('Relative level (dB) 2-ray model')

% 4-rays model PDP.
delay4=[0 delayImageTxRx delayImageTx1Rx delayImageTx2Rx];
power4=[0 PimageTxRx PimageTx1Rx PimageTx2Rx];

figure,stem2D(delay4,power4,floor(min(power4))-5)
xlabel('Excess delay (s)')
ylabel('Relative level (dB) 4-ray model')

% 6-rays model PDP.
delay6=[0 delayImageTxRx delayImageTx1Rx delayImageTx2Rx delayImageTx12Rx delayImageTx21Rx];
power6=[0 PimageTxRx PimageTx1Rx PimageTx2Rx PimageTx12Rx PimageTx21Rx];

figure,stem2D(delay6,power6,floor(min(power6))-5)
xlabel('Excess delay (s)')
ylabel('Relative level (dB) 6-ray model')

% PDP parameters;
[D2,S2]=PDPparameters(delay2,power2);
[D4,S4]=PDPparameters(delay4,power4);
[D6,S6]=PDPparameters(delay6,power6);

fprintf('\n2-rays model\n Excess delay \t %g\n Delay spread \t %g\n\n',D2,S2);
fprintf('4-rays model\n Excess delay \t %g\n Delay spread \t %g\n\n',D4,S4);
fprintf('6-rays model\n Excess delay \t %g\n Delay spread \t %g\n\n',D6,S6);