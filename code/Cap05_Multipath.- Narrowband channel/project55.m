%
% p5oject55
%
%=============================== RESET ==================================
clear
close all
clc
% basic inputs ==========================================================

fc=2000;     % MHz  Carrier frequency
F=100;        % sampling rate: fraction of wave length
V=10;        %  m/s MS1 speed 
NFFT=128;    % Number of points in FFT
Nsamples=1000 % Number of samples 
avPower=-20;   % 0.5*sigma^2  RMSsquared value

% geometry inputs ========================================================

dBS=1000;     
angleBS=135;
BSx=dBS*cosd(angleBS) % location of transmitter (BS) x-coordinate
BSy=dBS*sind(angleBS)  % location of transmitter (BS) y-coordinate


% locations of point scatterers =========================================

SC=[100 100
    -100 50
    -40 30
    100 70
    -70 -80
    -30 -60
    5 120
    -40 110
    0 -110
    -60 30
    50 -60
    -80 45
    -45 -80];

SCx=SC(:,1);
SCy=SC(:,2);

NSC=length(SCx); % Number of scatterers; 

figure,plot(SCx,SCy,'k*', BSx,BSy,'k^'), hold on

% indirect parameters ===================================================
lambdac=300/fc;    % m wavelength
Dx=lambdac/F;      % m sampling spacing 
ts=Dx/V;          % s time sampling interval
fs=1/ts;           % Hz sampling frequency
kc=2*pi/lambdac;   % propagation constant

a=sqrt(10.^(avPower/10)/NSC)  % magnitude of echoes
sigma=sqrt(0.5*10.^(avPower/10))     % Rayleigh parameter

fm=V/lambdac;       % max Doppler shift

timeaxis=ts.*[0:Nsamples-1];

MS0=-V*timeaxis(end)/2;        % initial location of receiver (MS) x-coordinate

MSx=MS0+V.*timeaxis;  % MS route along x-axis
MSy=zeros(Nsamples);  % MS route along x-axis (y=0)
plot(MSx,MSy,'k')
xlabel('Distance (m)');
ylabel('Distance (m)');

MINx=min(min(BSx, SCx))-100;
MAXx=max(max(BSx, SCx))+100;
MINy=min(min(min(BSy, SCy)))-100;
MAXy=max(max(max(BSy, SCy)))+100;
axis([MINx MAXx MINy MAXy])


% calculate distance matrix =============================================
distBSSC=sqrt((BSx-SCx).^2+(BSy-SCy).^2);

distBSSCext=repmat(distBSSC,1,Nsamples);

distSCMS=zeros(NSC,Nsamples);
for ii=1:Nsamples
    distSCMS(:,ii)=sqrt((SCx-MSx(ii)).^2+SCy.^2);
end

distBSSCMS=distBSSCext+distSCMS;

% calculate complex envelope ===========================================
ray=a*exp(-j*kc*distBSSCMS);
r=sum(ray);

% plot amplitude and phase =============================================
figure,plot(timeaxis,abs(r),'k')
xlabel('Time (s)')
ylabel('Magnitude of complex envelope')

% Autocorrelation of magnitude |r| with time/distance =====================

rho=xcorr(abs(r)-mean(abs(r)),'coeff');
timelags=([0:length(rho)-1]-length(rho)/2+1)*ts
distancelags=([0:length(rho)-1]-length(rho)/2+1)*Dx
rhotheoretical=besselj(0,timelags*2*pi*fm).^2;
figure,plot(timelags,rho,'k:',timelags,rhotheoretical,'k')
xlabel('Antenna spacing (s)')
ylabel('Autocorrelation coefficient')
title('Autocorrelation. Time domain') 
legend('Simulated','Theoretical')

figure,plot(distancelags,rho,'k:',distancelags,rhotheoretical,'k')
xlabel('Antenna spacing (m)')
ylabel('Autocorrelation coefficient')
title('Autocorrelation. Space domain') 
legend('Simulated','Theoretical')

% =================== Switch diversity ==================================
% ==== at 1/8 lambda

offset0125=round(F/8);
series1=abs(r(1:length(r)-offset0125));
series2=abs(r(offset0125+1:length(r)));
divseries0125=max(series1,series2);
figure,plot(timeaxis(1:length(series1)),20*log10(series1),'k:', ...
    timeaxis(1:length(series1)),20*log10(series2),'k-.',...
    timeaxis(1:length(series1)), 20*log10(divseries0125),'k')
xlabel('Time (s)')
ylabel('Relative received signal level (dB)')
title('Selection diversity. Spacing 1/8 \lambda. S1, S2 and Combination of S1 and S2')

% ==== at 1/4 lambda

offset025=round(F/4);
series1=abs(r(1:length(r)-offset025));
series2=abs(r(offset025+1:length(r)));
divseries025=max(series1,series2);
figure,plot(timeaxis(1:length(series1)),20*log10(series1),'k:', ...
    timeaxis(1:length(series1)),20*log10(series2),'k-.',...
    timeaxis(1:length(series1)), 20*log10(divseries025),'k')
xlabel('Time (s)')
ylabel('Relative received signal level (dB)')
title('Selection diversity. Spacing 1/4 \lambda. S1, S2 and Combination of S1 and S2')

% ==== at 1/2 lambda

offset05=round(F/2);
series1=abs(r(1:length(r)-offset05));
series2=abs(r(offset05+1:length(r)));
divseries05=max(series1,series2);
figure,plot(timeaxis(1:length(series1)),20*log10(series1),'k:', ...
    timeaxis(1:length(series1)),20*log10(series2),'k-.',...
    timeaxis(1:length(series1)), 20*log10(divseries05),'k')
xlabel('Time (s)')
ylabel('Relative received signal level (dB)')
title('Selection diversity. Spacing 1/2 \lambda. S1, S2 and Combination of S1 and S2')

% ==== at 1 lambda

offset1=F;
series1=abs(r(1:length(r)-offset1));
series2=abs(r(offset1+1:length(r)));
divseries1=max(series1,series2);
figure,plot(timeaxis(1:length(series1)),20*log10(series1),'k:', ...
    timeaxis(1:length(series1)),20*log10(series2),'k-.',...
    timeaxis(1:length(series1)), 20*log10(divseries1),'k')
xlabel('Time (s)')
ylabel('Relative received signal level (dB)')
title('Selection diversity. Spacing \lambda. S1, S2 and Combination of S1 and S2')

% ======= calculating diversity gain using CDF ====================

[CDFx,CDFy]=fCDF(20*log10(abs(r)));
[CDFx0125,CDFy0125]=fCDF(20*log10(divseries0125));
[CDFx025,CDFy025]=fCDF(20*log10(divseries025));
[CDFx05,CDFy05]=fCDF(20*log10(divseries05));
[CDFx1,CDFy1]=fCDF(20*log10(divseries1));
figure, semilogy(CDFx,CDFy,'k', CDFx0125,CDFy0125,'k:', CDFx025,CDFy025,'k-.',CDFx05,CDFy05,'k--',CDFx1,CDFy1,'k.')
xlabel('Relative received signal level (dB)')
ylabel('Probability the abscissa is not exceeded')
legend('w/o diversity','\lambda/8 separation','\lambda/4 separation','\lambda/2 separation','\lambda separation', 'Location', 'SouthEast')

