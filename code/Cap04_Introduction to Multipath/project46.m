%
% project46
%
% =======================================================================
clear
close all
clc
% basic inputs ==========================================================

fc=2000;     % MHz  Carrier frequency
F=16;        % sampling rate: fraction of wave length
V=10;        %  m/s MS1 speed 
NFFT=128;    % Number of points in FFT
Nsamples=100000 % Number of samples 

% geometry inputs ========================================================

BSx=-1000; % location of transmitter (BS) x-coordinate
BSy=1000;  % location of transmitter (BS) y-coordinate
MS0=0;        % initial location of receiver (MS) x-coordinate

% locations of point scatterers =========================================

SC=[100 100
    -100 50
    -40 30
    100 70
    -70 -80
    -30 -60
    0 120
    -40 110
    0 -110
    -60 30
    50 -60
    -80 45
    -45 -80];

SCx=SC(:,1);
SCy=SC(:,2);

NSC=length(SCx); % Number of scatterers; 

figure,plot(SCx,SCy,'*', BSx,BSy,'^'), hold on

% indirect parameters ===================================================
lambdac=300/fc;    % m wavelength
Dx=lambdac/F;      % m sampling spacing 
ts=Dx/V;          % s time sampling interval
fs=1/ts;           % Hz sampling frequency
kc=2*pi/lambdac;   % propagation constant

timeaxis=ts.*[0:Nsamples-1];

MSx=MS0;  % MS is stationary at orign 
MSy=0; ;  % MS is stationary at orign 
plot(MSx,MSy,'+r')
xlabel('Distance (m)');
ylabel('Distance (m)');

MINx=min(min(BSx, SCx))-100;
MAXx=max(max(BSx, SCx))+100;
MINy=min(min(min(BSy, SCy)))-100;
MAXy=max(max(max(BSy, SCy)))+100;
axis([MINx MAXx MINy MAXy])


% calculate distances ==================================================
distBSSC=sqrt((BSx-SCx).^2+(BSy-SCy).^2);

distSCMS=sqrt((SCx-MSx).^2+SCy.^2);

distBSSCMS=distBSSC+distSCMS;

% introduce slowly varuing phases ==========================================

tSlopes=[0.01 
    -0.2 
    -0.03 
    0.1 
    0.6 
    -0.25 
    0.1 
    0.2 
    -0.13 
    0.9 
    -0.2 
    -0.4 
    0.3]; 

SCphases=tSlopes*timeaxis;
figure,plot(timeaxis,SCphases)
xlabel('Time (s)')
ylabel('Scatterer phases (rad)')

% calculate complex envelope ===========================================

distBSSCMSext=repmat(distBSSCMS,1,Nsamples); 

ray=exp(-j*kc*distBSSCMSext+j*SCphases);
r=sum(ray);

% plot amplitude and phase =============================================
figure,plot(timeaxis,abs(r))
xlabel('Time (s)')
ylabel('Magnitude of complex envelope')

figure,plot(timeaxis,unwrap(angle(r)))
xlabel('Time (s)')
ylabel('Phase of complex envelope (rad)')

% plot normalized spectrum in dB =======================================
spectrumr=fftshift((abs(fft(r,NFFT))).^2);
freqaxis=[0:NFFT-1]*fs/NFFT-fs/2;
figure,plot(freqaxis,10*log10(spectrumr)-max(10*log10(spectrumr)))
xlabel('Doppler shift (Hz)')
ylabel('Normalized frquency response (dB)')

% Calculate CDF ========================================================
[CDFx,CDFy]=fCDF(abs(r));

sigma=sqrt(NSC/2);
[CDFyTH]=RayleighCDF(sigma,CDFx)
figure,plot(CDFx,CDFy,CDFx,CDFyTH)
legend('Simulated CDF', 'Theoretical CDF', 'Location', 'SouthEast');
xlabel('Magnitude complex envelope');
ylabel('Probability the abscissa is not exceeded');






