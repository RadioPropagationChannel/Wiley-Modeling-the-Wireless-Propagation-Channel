%
% project513 width +/-60 degree, direction 90 degree
%
% Initialize ============================================================
clear
close all
clc
% basic inputs ==========================================================

fc=2000;         % MHz  Carrier frequency
F=16;             % sampling rate: fraction of wave length
V=10;            %  m/s MS1 speed 
NFFT=1024;       % Number of points in FFT
Nsamples=1000;   % Number of samples 
NSC=1000;        % Number of scatterers
avPower=-20;     % sigma^2  Raverage power

% geometry inputs ========================================================

dBS=1000;     
angleBS=135;
BSx=dBS*cosd(angleBS) % location of transmitter (BS) x-coordinate
BSy=dBS*sind(angleBS)  % location of transmitter (BS) y-coordinate

% locations of point scatterers =========================================

minalpha=30;     % -60
maxalpha=150;    % +60


D=200;                        % radius from origin
% alpha=rand(NSC,1)*360;        % random draw of angles of arrival
alpha=rand(NSC,1)*(maxalpha-minalpha)+minalpha;        % random draw of angles of arrival

SCx=D.*cosd(alpha);
SCy=D.*sind(alpha);

figure,plot(SCx,SCy,'*', BSx,BSy,'^'), hold on

% indirect parameters ===================================================
lambdac=300/fc;    % m wavelength
Dx=lambdac/F;      % m sampling spacing 
ts=Dx/V;          % s time sampling interval
fs=1/ts;           % Hz sampling frequency
kc=2*pi/lambdac;   % propagation constant

a=sqrt(10.^(avPower/10)/NSC)  % magnitude of echoes
sigma=sqrt(0.5*10.^(avPower/10))     % Rayleigh parameter

fm=V/lambdac                % max Doppler shift


timeaxis=ts.*[0:Nsamples-1];

MS0=-V*timeaxis(end)/2;        % initial location of receiver (MS) x-coordinate

MSx=MS0+V.*timeaxis;  % MS route along x-axis
MSy=zeros(Nsamples,1);  % MS route along x-axis (y=0)
plot(MSx,MSy,'r')
xlabel('Distance (m)');
ylabel('Distance (m)');

MINx=min(min(BSx, SCx))-100;
MAXx=max(max(BSx, SCx))+100;
MINy=min(min(min(BSy, SCy)))-100;
MAXy=max(max(max(BSy, SCy)))+100;
axis([MINx MAXx MINy MAXy])
axis equal

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
figure,plot(timeaxis,abs(r))
xlabel('Time (s)')
ylabel('Magnitude of complex envelope')

figure,plot(timeaxis,20*log10(abs(r)))
xlabel('Time (s)')
ylabel('Magnitude of complex envelope (dB)')

figure,plot(timeaxis,angle(r))
xlabel('Time (s)')
ylabel('Modulo-\pi phase of complex envelope (rad)')

figure,plot(timeaxis,unwrap(angle(r)))
xlabel('Time (s)')
ylabel('Phase of complex envelope (rad)')

% plot normalized RF spectrum ===============================
spectrumr=fftshift((abs(fft(r,NFFT))).^2);
freqaxis=[0:NFFT-1]*fs/NFFT-fs/2;
%figure,plot(freqaxis,10*log10(spectrumr)-max(10*log10(spectrumr)))
figure,plot(freqaxis,spectrumr)
xlabel('Doppler shift (Hz)')
ylabel('Complex envelope spectrum');


