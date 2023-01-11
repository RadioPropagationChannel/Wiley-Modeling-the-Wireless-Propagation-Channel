%
% project421
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
BSy=0;  % location of transmitter (BS) y-coordinate
MS0=0;        % initial location of receiver (MS) x-coordinate
MSy=0;

alpha1 = 0;       % degree
distSC1=1000; % m distance from origin (first route point) to sc2

% inidirect gemeotric parameters ========================================

SC1x=distSC1*cosd(alpha1);  % loc of scatterer-1 x-coord
SC1y=distSC1*sind(alpha1);  % loc of scatterer-1 y-coord


% indirect parameters ===================================================
lambdac=300/fc;    % m wavelength
Dx=lambdac/F;      % m sampling spacing 
ts=Dx/V;          % s time sampling interval
fs=1/ts;           % Hz sampling frequency
kc=2*pi/lambdac;   % propagation constant

timeaxis=ts.*[0:Nsamples];

MSx=MS0+V.*timeaxis;

distBSMS=sqrt((BSx-MSx).^2+(BSy-MSy).^2);
distBSSC1MS=sqrt((BSx-SC1x).^2+(BSy-SC1y).^2)+sqrt((SC1x-MSx).^2+SC1y.^2);
figure,plot(timeaxis,distBSMS,timeaxis,distBSSC1MS)
axisFig = axis;
axis([axisFig(1) axisFig(2) axisFig(3)-100 axisFig(4)+100]);
xlabel('Time (s)');
ylabel('Distance (m)');
legend('Distance between BS and MS', 'Distance between BS, scatterer and MS', 'Location', 'Best');

ray1=exp(-j*kc.*distBSMS);
ray2=-exp(-j*kc.*distBSSC1MS);   % including ref coefficient -1
figure,plot(timeaxis,unwrap(angle(ray1)),timeaxis,unwrap(angle(ray2)))
xlabel('Time (s)');
ylabel('Phase (rad)');
legend('Direct ray', 'Reflected ray', 'Location', 'Best');

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


