%
% project811
%
%=======================================================================
clear
close all
clc
%=======================================================================
f=2000;          % Frequency in MHz
ptx=1;           % Transmitter power in Watts
gtx=1;           % Transmitter gain in linear units (1 is isotropic) 

%======================================================================
lambdac=300/f;     % wavelength (m)
kc=2*pi/lambdac;   % propagation constant

%=======================================================================
ht=10;        % transmitter height
hr=1.5        % receiver height

xt=0;

xrStep=1;
xr=[1:xrStep:100000];

% Reflection coefficients ===============================================
RG=-1;

% ray path lengths =======================================================
imageht=-ht;

dTxRx=sqrt((xt-xr).^2+(ht-hr).^2);
dimageTxRx=sqrt((xt-xr).^2+(imageht-hr).^2);

% Field strengths =======================================================
E0=sqrt(60*ptx*gtx)*exp(-j*kc*dTxRx)./dTxRx;      % in V/m

E=sqrt(60*ptx*gtx)*(exp(-j*kc*dTxRx)./dTxRx...
    +RG*exp(-j*kc*dimageTxRx)./dimageTxRx);

% field strenghts are in V/m, now they are converted to dBuV/m
figure,semilogx(xr,20*log10(abs(E0)*1e6),'k:',xr,20*log10(abs(E)*1e6),'k', 'LineWidth',2)
xlabel('Distance from Tx (m)')
ylabel('Received field strength (dB \muV/m)')
title('Recived field: free space, 2-ray model')
legend('free space','2-ray model')

figure,semilogx(xr,20*log10(abs(E./E0)),'k', 'LineWidth',2)
xlabel('Distance from Tx (m)')
ylabel('Relative received field strength (dB \muV/m)')
title('Relative recived field strength. 2-ray model')




