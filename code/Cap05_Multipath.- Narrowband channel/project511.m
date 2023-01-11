%
% project511 Clarke's model
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

minalpha=0;
maxalpha=360;


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
ylabel('Absolute phase of complex envelope (rad)')

% plot normalized RF spectrum ===============================
spectrumr=fftshift((abs(fft(r,NFFT))).^2);
freqaxis=[0:NFFT-1]*fs/NFFT-fs/2;
%figure,plot(freqaxis,10*log10(spectrumr)-max(10*log10(spectrumr)))
figure,plot(freqaxis,spectrumr)
xlabel('Doppler shift (Hz)')
ylabel('Complex envelope Spectrum');


% plot normalized BB spectrum (magnitude of envelope) ===================
spectrumrBF=fftshift(abs((fft(abs(r)-mean(abs(r)),NFFT))).^2);
freqaxis=[0:NFFT-1]*fs/NFFT-fs/2;

figure,plot(freqaxis,10*log10(spectrumrBF)-max(10*log10(spectrumrBF)))
% figure,plot(freqaxis,spectrumrBF)
xlabel('Frequency (Hz)')
ylabel('Spectrum of magnitude of complex envelope');

% calculate autocorrelation function ===================================
correl=xcorr(abs(r)-mean(abs(r)),'coeff');
correlaxis=([0:length(correl)-1]-length(correl)/2+1)*Dx;

correlTheoretical=besselj(0,correlaxis*2*pi/lambdac);

figure,plot(correlaxis,correl,'k',correlaxis,correlTheoretical,'k:')
axis([0 1 -0.35 1])
xlabel('Spacing (m)')
ylabel('Complex envelope autocorrelation coefficient')
legend('Simulation','Theoretical')


% Calculate CDF of magnitude =======================================
[CDFx,CDFy]=fCDF(abs(r));

[CDFyTH]=RayleighCDF(sigma,CDFx);
figure,plot(CDFx,CDFy,CDFx,CDFyTH)
xlabel('Signal magnitude')
ylabel('Probability the abscissa is not exceeded')
legend('Simulation','Theoretical')
aaux=axis;
axis([aaux(1) aaux(2)  0 1])

% Calculate histogram of phase =======================================

[histy,histx]=hist(angle(r),20);
histy=histy/length(angle(r));     % normalize histigram for total probability equal one
figure,bar(histx,histy)
xlabel('Phase (degree)');
ylabel('Probability');


% Normalize signal amplitude wrt to rms value

rmsr=sqrt(std(abs(r))^2+mean(abs(r))^2);      % calculate rms value of signal magnitide
rrho=abs(r)/rmsr;                             % normalize magnitide wrt rms value
figure,plot(timeaxis,20*log10(rrho))   % plot normalized signal envelope wrt rms value (dB)
xlabel('Time (s)')
ylabel('Magnitude of complex envelope (dB/RMS)')

% average fade curation, afd ====================================
[axisafd,afd]=afduration(rrho,ts);            % returns afd with abscissas in dB/RMS and ordintes in s
afdinwavelengths = afd*V/lambdac;             % for converting afd in s to wavelengths

afdtheoretical= (exp(10.^(axisafd./10))-1)./(sqrt(2*pi)*10.^(axisafd./20));

figure, semilogy(axisafd,afdinwavelengths,axisafd,afdtheoretical); %inwavelengths)
xlabel('Signal level (dB/RMS)')
ylabel('Average fade durations (wavelengths)')
title('Average fade durations , afd ')
axis([min(axisafd) max(axisafd)  0.01 10])


% Distribution of fade durations wrt lev dB ====================

lev=[-25 -20 -15 -10 -5 0 5];


% figure
% hold on
% for ii=1:length(lev)
%     [durationsatlevel]=durationoffades(rrho,ts,lev(ii));
%     durationsatlevel=durationsatlevel*V/lambdac;         % convert to wavelenghts
%     [x,y]=fCDF(durationsatlevel);
%     plot(x,y)
% end
% 
% ylabel('Probability the abscissa is not exceeded')
% w=['Duration (s) of fades below ' num2str(lev) ' dB/RMS in wavelengths'];
% xlabel(w)
% z=['Distribution of durationes at  ' num2str(lev) ' dB/RMS'];
% title(z) 
% axisdurs=axis; axis([axisdurs(1) axisdurs(2) 0 1])
% hold off

% level crossing rate calculations =====================================

[axislcr,lcr]=lcrate(rrho,ts);
lcr=lcr*lambdac/V;                % convert to crossings per wavelength

lcrtheoretical=sqrt(2*pi)*10.^(axislcr./20).*exp(-10.^(axislcr./10));
figure,semilogy(axislcr,lcr,axislcr,lcrtheoretical)
xlabel('Signal level (dB/RMS)')
ylabel('Level crossing rate, crossings/wavelenght')
title('Level crossing rate')
axis([min(axislcr) max(axislcr) 0.01 2])

% random FM ===========================================================

randomFM=diff(angle(r))/ts;
randomFMnorm=randomFM/(2*pi*fm);
axisrandomFM=[0:length(randomFM)-1]*ts;

figure
subplot(3,1,1)
plot(timeaxis(1:1000),20*log10(abs(r(1:1000))),'k')
ylabel('Level (dB/LOS)')
subplot(3,1,2)
plot(timeaxis(1:1000),angle(r(1:1000)),'k')
ylabel('Phase (Rad.)')
subplot(3,1,3)
plot(axisrandomFM(1:999), randomFM(1:999),'k')
ylabel('Random FM')
xlabel('Time (s)')
% histogram of random FM
[y,x]=hist(randomFMnorm,300);
y=y/length(randomFMnorm);
figure,bar(x,y)
axishistrandomFM=axis;
% axis([-10000 10000 axishistrandomFM(3) axishistrandomFM(4)])
xlabel('Random FM (Hz)');
ylabel('Probability');

% Cumulative distribution of random FM
[CDFx,CDFy]=fCDF(randomFMnorm);
figure,plot(CDFx,CDFy)
% axis([-10000 10000 0 1])
xlabel('Random FM (Hz)');
ylabel('Probability the abscissa is not exceeded');
