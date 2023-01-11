%
% project312    (spline interpolation)
%
%========================================================================
clear
close all
clc
%========================================================================
fMHz=2000;     % Frequency in MHz
Lcorr=30;       % Correlation distance in meters 
MM=-70;          % Larger area mean in dBm
SS=7;          % Larger area std or location variability in dB

ThresHold=-80;

%=======================================================================
lambdac=300/fMHz;   % wavelength in m
ds=lambdac;
Nsamples=200;
InterpRate=round(Lcorr/ds)
Lcorr=InterpRate*ds       % slightly correct Lcorr to make it a multiple of ds

% =====================================================

% Create Gaussian series ================================================

R=randn(Nsamples,1);
d_axis1=[1:Nsamples]*Lcorr-Lcorr;

d_axis2=([1:1/InterpRate:Nsamples]-1)*Lcorr;

Rinterpolated=interp1(d_axis1,R,d_axis2,'spline');
Rinterpolated=Rinterpolated*SS+MM;

mean(Rinterpolated)
std(Rinterpolated)

plot(d_axis1,R*SS+MM,'k^',d_axis2,Rinterpolated,'k','LineWidth',2)
xlabel('Traversed distance (m)')
ylabel('Slow signal variations (dBm)')
ylabel('Slow variations autocorrelation')


figure,plot([d_axis2(1); d_axis2(end)],[ThresHold; ThresHold],'k:', d_axis2,Rinterpolated,'k','LineWidth',2)
xlabel('Traversed distance (m)')
ylabel('Slow signal variations (dBm)')
ylabel('Slow variations autocorrelation')

% evaluate coverage probability
Nbelow=find(Rinterpolated<ThresHold);
disp('Coverage probability')
CovProb=1-(length(Nbelow)/length(Rinterpolated))

% normalized Gaussian parameters for calculating theoretical area coverage
k=(ThresHold-MM)/SS;
ThCovProb=0.5*erfc(k/sqrt(2))

    