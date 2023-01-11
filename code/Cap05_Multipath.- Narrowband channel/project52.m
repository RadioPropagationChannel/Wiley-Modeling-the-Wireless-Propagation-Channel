% RICE CASE

clear
close all

% basic inputs ==========================================================

fc=2000;     % MHz  Carrier frequency
F=16;        % sampling rate: fraction of wave length
V=10;        %  m/s MS1 speed 
NFFT=1024;    % Number of points in FFT
Nsamples=10000; % Number of samples 
NSC=100;       % Number of scatterers
avPower=-20;   % Average local power
A0=0;           % Level of direct signal 

% geometry inputs ========================================================

dBS=1000;     
angleBS=15;
BSx=dBS*cosd(angleBS) % location of transmitter (BS) x-coordinate
BSy=dBS*sind(angleBS)  % location of transmitter (BS) y-coordinate

% locations of point scatterers =========================================

D=200;                        % radius from origin
alpha=rand(NSC,1)*360;        % random draw of angles of arrival

SCx=D.*cosd(alpha);
SCy=D.*sind(alpha);

figure,plot(SCx,SCy,'k*', BSx,BSy,'k^'), hold on

% indirect parameters ===================================================
lambdac=300/fc;    % m wavelength
Dx=lambdac/F;      % m sampling spacing 
ts=Dx/V;          % s time sampling interval
fs=1/ts;           % Hz sampling frequency
kc=2*pi/lambdac;   % propagation constant

a=sqrt(10.^(avPower/10)/NSC)  % magnitude of echoes
sigma=sqrt(0.5*10.^(avPower/10))     % Rayleigh parameter

a0=10^(A0/20);                 % Magnitude of direct signal
K=10*log10(a0^2/2/sigma^2)     % Rice K-Factor (dB)
fm=V/lambdac                % max Doppler shift (Hz)


timeaxis=ts.*[0:Nsamples-1];

MS0=-V*timeaxis(end)/2;        % initial location of receiver (MS) x-coordinate

MSx=MS0+V.*timeaxis;  % MS route along x-axis
MSy=zeros(Nsamples,1);  % MS route along x-axis (y=0)
plot(MSx,MSy,'r')
xlabel('Propagation scenario. Distance (m)')
ylabel('Propagation scenario. Distance (m)')
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

% Direct signal distance vector ======================================== 
distBSMS=sqrt((BSx-MSx').^2+(BSy-MSy).^2);

% calculate complex envelope ===========================================
ray=a*exp(-j*kc*distBSSCMS);
r=sum(ray);                     % without direct signal

r=r+a0*exp(-j*kc*distBSMS');

% plot amplitude and phase =============================================
figure,plot(timeaxis,abs(r),'k')
xlabel('Time (s)')
ylabel('Magnitude of complex envelope')

figure,plot(timeaxis,20*log10(abs(r)),'k')
xlabel('Time (s)')
ylabel('Magnitude of complex envelope (dB)')
axisRicedB=axis;
axis([axisRicedB(1) axisRicedB(2) -25 5 ])

figure,plot(timeaxis,angle(r),'k')
xlabel('Time (s)')
ylabel('Modulo-\pi phase of complex envelope (rad)')

figure,plot(timeaxis,unwrap(angle(r)),'k')
xlabel('Time (s)')
ylabel('Phase of complex envelope (rad)')

% plot normalized RF spectrum =======================================
spectrumr=fftshift((abs(fft(r,NFFT))).^2);
freqaxis=[0:NFFT-1]*fs/NFFT-fs/2;
figure,plot(freqaxis,10*log10(spectrumr)-max(10*log10(spectrumr)),'k')
xlabel('Doppler shift (Hz)')
ylabel('Spectrum of complex envelope (dB/max)')


% plot normalized BB spectrum (magnitude of envelope) ===================
spectrumrBF=fftshift((fft(abs(r)-mean(abs(r)),NFFT)).^2);
freqaxis=[0:NFFT-1]*fs/NFFT-fs/2;
figure,plot(freqaxis,10*log10(spectrumrBF)-max(10*log10(spectrumrBF)),'k')
xlabel('Frequency (Hz)')
ylabel('Frequency response (dB/max)')



% Calculate CDF of simulated magnitude and theoretical ==============
[CDFx,CDFy]=fCDF(abs(r));

yrct=[];
for ii=0:0.01:ceil(max(abs(r))),
    yrct=[yrct riceCDF(a0,sigma,ii)];
end
xrct=[0:0.01:ceil(max(abs(r)))];
figure, plot(CDFx,CDFy,'r*',xrct,yrct,'k')
axis([0 2 0 1])
grid
ylabel('Probability the abscissa is not exceeded')
xlabel('Normalized signal magnitude')
legend('Simulated','Theoretical')
