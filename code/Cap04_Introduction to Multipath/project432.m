%
% project432     (alpha 90)
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
Nsamples=100 % Number of samples 

% geometry inputs ========================================================

BSx=-1000; % location of transmitter (BS) x-coordinate
BSy=1000;  % location of transmitter (BS) y-coordinate
MS0=0;        % initial location of receiver (MS) x-coordinate
SC1x=100;    % location of scatterer-1 x-coordinate
SC1y=0;      % location of scatterer-1 y-coordinate

alpha2 = 90;       % degree
distSC2=100; % m distance from origin (first route point) to sc2

% inidirect gemeotric parameters ========================================

SC2x=distSC2*cosd(alpha2);  % loc of scatterer-1 x-coord
SC2y=distSC2*sind(alpha2);  % loc of scatterer-1 y-coord


% indirect parameters ===================================================
lambdac=300/fc;    % m wavelength
Dx=lambdac/F;      % m sampling spacing 
ts=Dx/V;          % s time sampling interval
fs=1/ts;           % Hz sampling frequency
kc=2*pi/lambdac;   % propagation constant

timeaxis=ts.*[0:Nsamples];

MSx=MS0+V.*timeaxis;

distBSSC1MS=sqrt((BSx-SC1x).^2+(BSy-SC1y).^2)+sqrt((SC1x-MSx).^2+SC1y.^2)
distBSSC2MS=sqrt((BSx-SC2x).^2+(BSy-SC2y).^2)+sqrt((SC2x-MSx).^2+SC2y.^2)
figure,plot(timeaxis,distBSSC1MS,timeaxis,distBSSC2MS)
xlabel('Time (s)');
ylabel('Distance (m)');
legend('Distance between BS, scatterer 1 and MS', 'Distance between BS, scatterer 2 and MS', 'Location', 'Best');


ray1=exp(-j*kc.*distBSSC1MS);
ray2=exp(-j*kc.*distBSSC2MS);
figure,plot(timeaxis,unwrap(angle(ray1)),timeaxis,unwrap(angle(ray2)))
xlabel('Time (s)');
ylabel('Phase (rad)');
legend('Ray of scatterer 1', 'Ray of scatterer 2', 'Location', 'Best');

r=ray1+ray2;
figure,plot(timeaxis,abs(r))
xlabel('Time (s)')
ylabel('Magnitude of complex envelope')

figure,plot(timeaxis,unwrap(angle(r)))
xlabel('Time (s)')
ylabel('Phase of complex envelope (rad)')

spectrumr=fftshift((abs(fft(r,NFFT))).^2);
freqaxis=[0:NFFT-1]*fs/NFFT-fs/2;
figure,plot(freqaxis,10*log10(spectrumr)-max(10*log10(spectrumr)))
xlabel('Doppler shift (Hz)')
ylabel('Normalized frquency response (dB)')


