%
% project812
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


xt=0;
xr=200;       % MS at a fixed distance 

hrStep=0.1;            
hr=[0.1:hrStep:10];        % receiver height

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
figure,semilogy(20*log10(abs(E0)*1e6),hr,'k:',20*log10(abs(E)*1e6),hr,'k','LineWidth',2)
ylabel('Rx antenna height (m)')
xlabel('Received field strength (dB \muV/m)')
title('Recived field: free space, 2-ray model')
legend('free space','2-ray model')

