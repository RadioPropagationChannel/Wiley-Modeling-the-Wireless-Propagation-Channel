%
% project45 
%
% =======================================================================
clear
close all
clc
% basic inputs ==========================================================

fc=200;     % MHz  Carrier frequency
F=50;        % sampling rate: fraction of wave length
V=10;        %  m/s MS1 speed 
NFFT=128;    % Number of points in FFT
Nmeters=3;   % Number of meters in the distance axis.

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

figure,plot(SCx,SCy,'*', BSx,BSy,'^'), hold on

% indirect parameters ===================================================
lambdac=300/fc;    % m wavelength
Dx=lambdac/F;      % m sampling spacing 
ts=Dx/V;           % s time sampling interval
fs=1/ts;           % Hz sampling frequency
kc=2*pi/lambdac;   % propagation constant
Nsamples=Nmeters/Dx; % Number of samples

% timeaxis=ts.*[0:Nsamples-1];
distanceaxis=Dx.*[0:Nsamples-1];

% MSx=MS0+V.*timeaxis;  % MS route along x-axis
MSx=MS0+distanceaxis;  % MS route along x-axis
MSy=repmat(distanceaxis',1,length(distanceaxis));     % MS routes along x-axis for different y values

for m=1:length(MSy)
    plot(MSx,MSy(m,:),'r')
end
xlabel('Distance (m)');
ylabel('Distance (m)');

MINx=min(min(BSx, SCx))-100;
MAXx=max(max(BSx, SCx))+100;
MINy=min(min(min(BSy, SCy)))-100;
MAXy=max(max(max(BSy, SCy)))+100;
axis([MINx MAXx MINy MAXy])


% calculate distance matrix =============================================
distBSSC=sqrt((BSx-SCx).^2+(BSy-SCy).^2);

distBSSCext=repmat(distBSSC,[1 Nsamples Nsamples]);

distSCMS=zeros(NSC,Nsamples,Nsamples);

for jj=1:Nsamples
    for ii=1:Nsamples
        distSCMS(:,ii,jj)=sqrt((SCx-MSx(ii)).^2+(SCy-MSy(jj,1)).^2);
    end
end

distBSSCMS=distBSSCext+distSCMS;

% calculate complex envelope ===========================================
ray=exp(-j*kc*distBSSCMS);
ra=sum(ray);
r(:,:)=ra(1,:,:);

figure,surf(distanceaxis,distanceaxis,abs(r));
xlabel('Traveled distance (m)');
ylabel('Traveled distance (m)');
zlabel('Magnitude of complex envelope');

% % plot amplitude and phase =============================================
% figure,plot(timeaxis,abs(r(:,1)))
% xlabel('Time (s)')
% ylabel('Magnitude of complex envelope')
% 
% figure,plot(distanceaxis,abs(r(:,1)))
% xlabel('Distance (m)')
% ylabel('Magnitude of complex envelope')
% 
% figure,plot(timeaxis,unwrap(angle(r(:,1))))
% xlabel('Time (s)')
% ylabel('Phase of complex envelope (rad)')
% 
% % plot normalized spectrum in dB =======================================
% spectrumr=fftshift((abs(fft(r(:,1),NFFT))).^2);
% freqaxis=[0:NFFT-1]*fs/NFFT-fs/2;
% figure,plot(freqaxis,10*log10(spectrumr)-max(10*log10(spectrumr)))
% xlabel('Doppler shift (Hz)')
% ylabel('Normalized frquency response (dB)')
% 
% % Calculate CDF ========================================================
% [CDFx,CDFy]=fCDF(abs(r(:,1)));
% 
% sigma=sqrt(NSC/2);
% [CDFyTH]=RayleighCDF(sigma,CDFx)
% figure,plot(CDFx,CDFy,CDFx,CDFyTH)