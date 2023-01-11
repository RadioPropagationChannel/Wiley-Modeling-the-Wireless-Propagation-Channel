%
% project413
%
% ======================================================================
clear
close all
clc
% basic inputs ==========================================================

fc=2000;     % MHz  Carrier frequency
F=16;        % sampling rate: fraction of wave length
V=10;        %  m/s MS1 speed 
NFFT=128;    % Number of points in FFT
Nsamples=100 % Number of samples 

% geometry inputs ========================================================

dBS=1000;     % distance of BS to origin 
alpha = 120;   % degree. Angle of BS-MS with MS route 

% inidirect gemeotric parameters ========================================

BSx=dBS*cosd(alpha);  % loc of BS x-coord
BSy=dBS*sind(alpha);  % loc of BS y-coord

% indirect parameters ===================================================
lambdac=300/fc;    % m wavelength
Dx=lambdac/F;      % m sampling spacing 
ts=Dx/V;          % s time sampling interval
fs=1/ts;           % Hz sampling frequency
kc=2*pi/lambdac;   % propagation constant

timeaxis=ts.*[0:Nsamples];     % s  elapsed time axis
disaxis=Dx.*[0:Nsamples];      % n  traveled distance axis 

MSx=V.*timeaxis;    % MS route sampling points

% radio path length ==================================================== 

distBSMS=sqrt((BSx-MSx).^2+(BSy).^2);

figure,plot(timeaxis,distBSMS)
xlabel('Time (s)');
ylabel('Distance between transmitter (BS) and mobile receiver (MS) in meters');

% complex envelope: amplitude and phase ================================
r=exp(-j*kc.*distBSMS);

figure,plot(disaxis,abs(r))
xlabel('Traveled distance (m)')
ylabel('Magnitude of complex envelope')

figure,plot(timeaxis,abs(r))
xlabel('Time (s)')
ylabel('Magnitude of complex envelope')

figure,plot(disaxis,unwrap(angle(r)))
xlabel('Traveled distance (m)')
ylabel('Absolute phase of complex envelope (rad)')

figure,plot(timeaxis,unwrap(angle(r)))
xlabel('Time (s)')
ylabel('Absolute phase of complex envelope (rad)')

figure,plot(disaxis,angle(r))
xlabel('Traveled distance (m)')
ylabel('Modulo-\pi phase of complex envelope (rad)')
axis([0 1 -3.5 3.5])

% complex envelope spectrum ============================================
spectrumr=fftshift((abs(fft(r,NFFT))).^2);
freqaxis=[0:NFFT-1]*fs/NFFT-fs/2;
figure,plot(freqaxis,10*log10(spectrumr)-max(10*log10(spectrumr)))
xlabel('Doppler shift (Hz)')
ylabel('Normalized frquency response (dB)')


