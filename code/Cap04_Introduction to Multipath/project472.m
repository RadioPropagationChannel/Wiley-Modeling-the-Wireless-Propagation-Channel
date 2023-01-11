%
% project472
%
% =======================================================================
clear
close all
clc
% basic inputs ==========================================================

fc=2000;   % MHz  Carrier frequency
F=16;       % sampling rate: fraction of wave length
V1=5;      %  m/s MS1 speed 
V2=10;      % m/s  MS2 speed
V3=0;       % m/s scatterer speed
NFFT=128;   % Number of points in FFT
Nsamples=100 % Number of samples 

% geometry inputs ========================================================

trasmitter0=-500;  % initial location of transmitter (MS1) x-coordinate
receiver0=0;        % initial location of receiver (MS2) x-coordinate
scatterer0=500;     % initial location of scatterer x-coordinate

% indirect parameters ===================================================
lambdac=300/fc;   % m wavelength
Dx=lambdac/F;     % m sampling spacing 
ts=Dx/V2;         % s time sampling interval
fs=1/ts;           % Hz sampling frequency
kc=2*pi/lambdac;  % propagation constant

timeaxis=ts.*[0:Nsamples];

receiver=receiver0+V2.*timeaxis;
transmitter=trasmitter0+V1.*timeaxis;
scatterer=scatterer0+V3.*timeaxis;
figure,plot(timeaxis,receiver,timeaxis,transmitter,timeaxis,scatterer)
axisFig=axis;
axis([axisFig(1) axisFig(2) axisFig(3)-100 axisFig(4)+100]);
xlabel('Time (s)');
ylabel('Traveled distance (m)');
legend('Receiver', 'Transmitter', 'Scatterer', 'Location', 'Best');


distanceMS1MS2=abs(transmitter-receiver);
distanceMS1ScattererMS2=abs(scatterer-transmitter)+abs(receiver-scatterer);
figure,plot(timeaxis,distanceMS1MS2,timeaxis,distanceMS1ScattererMS2)
axisFig=axis;
axis([axisFig(1) axisFig(2) axisFig(3)-100 axisFig(4)+100]);
xlabel('Time (s)');
ylabel('Distance (m)');
legend('Distance between MS1 and MS2', 'Distance between MS1, scatterer and MS2', 'Location', 'Best');


ray0=exp(-j*kc.*distanceMS1MS2);
ray1=exp(-j*kc.*distanceMS1ScattererMS2);
figure,plot(timeaxis,unwrap(angle(ray0)),timeaxis,unwrap(angle(ray1)))
xlabel('Time (s)');
ylabel('Phase (rad)');
legend('Direct ray', 'Reflected ray', 'Location', 'Best');


r=ray0+ray1;
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


