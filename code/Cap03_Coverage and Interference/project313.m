%
% project313
%
%========================================================================
clear
close all
clc
%========================================================================
fMHz=2000;     % Frequency in MHz
Lcorr=30;       % Correlation distance in meters 

EIRP=30;        % EIRP dBm
AA=100;        % Loss at 1 km 
n=3.6;          % Propagation exponent
SS=7;           % Location variability in dB

%=======================================================================
lambdac=300/fMHz;   % wavelength in m
ds=lambdac;
Distance= 30;       % simulated distance (km)
InterpRate=round(Lcorr/ds)
Lcorr=InterpRate*ds       % slightly correct Lcorr to make it a multiple of ds

Nsamples=round((Distance-1)*1000/Lcorr);

% Create Gaussian series ================================================

d_axis1=1+[0:Nsamples-1]'*Lcorr/1000;

veryslowVars=EIRP-AA-10*n*log10(d_axis1);
slowVars=randn(Nsamples,1)*SS;
R=veryslowVars+slowVars;

figure,plot(d_axis1,veryslowVars,'k',d_axis1,R,'k:')
xlabel('Traversed distance (km)')
ylabel('Slow signal variations (dBm)')
ylabel('Slow variations autocorrelation')


d_axis2=(1000+([1:1/InterpRate:Nsamples]-1)*Lcorr)/1000;

Rinterpolated=interp1(d_axis1,R,d_axis2,'spline');
 
figure,plot(d_axis1,veryslowVars,'k',d_axis2,Rinterpolated,'k:')
xlabel('Traversed distance (km)')
ylabel('Slow signal variations (dBm)')
ylabel('Slow variations autocorrelation')
