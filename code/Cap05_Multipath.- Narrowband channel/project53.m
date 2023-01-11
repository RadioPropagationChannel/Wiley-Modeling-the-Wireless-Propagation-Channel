%
% project53
%
%========================================================================
clear
close all
clc
% basic inputs ==========================================================

fc=2000;       % MHz  Carrier frequency
F=16;          % sampling rate: fraction of wave length
V=10;          %  m/s MS1 speed 
NFFT=1024;     % Number of points in FFT
Nsamples=1000; % Number of samples 
avPower=-20;   % 0.5*sigma^2  RMSsquared value
N0= 8;        % number of sinusoidal generators -1
N=2*(2*N0+1);  % equivalent number of rays

% indirect parameters ===================================================
lambdac=300/fc;                 % m wavelength
Dx=lambdac/F;                   % m sampling spacing 
ts=Dx/V;                        % s time sampling interval
fs=1/ts;                        % Hz sampling frequency
kc=2*pi/lambdac;                % propagation constant

a=sqrt(2*10.^(avPower/10)/N)    % magnitude of echoes
sigma=sqrt(0.5*10.^(avPower/10))

fm=V/lambdac;                    % max Doppler shift

timeaxis=ts.*[0:Nsamples-1];

% generator ======================================================= 

nn=[1:N0]';
dopplerfreqs=fm*cos(2*pi*nn/N);

% magnitudedopplerfreqs=ones(length(dopplerfreqs),1)*2;
% magnitudedopplerfreqs=[magnitudedopplerfreqs; sqrt(2)];
% figure,stem([dopplerfreqs; fm ],magnitudedopplerfreqs)

beta=pi*nn/(N0+1);   % phases and amplifier gains
alpha=pi/4;             % phases and amplifier gains

betaaux=repmat(beta,1,length(timeaxis));
Iaux=2*cos(betaaux).*cos(2*pi*dopplerfreqs*timeaxis);
Qaux=2*sin(betaaux).*cos(2*pi*dopplerfreqs*timeaxis);

Iaux2=(1/sqrt(2))*2*cos(alpha).*cos(2*pi*fm*timeaxis);
Qaux2=(1/sqrt(2))*2*sin(alpha).*cos(2*pi*fm*timeaxis);

Iaux=[Iaux; Iaux2]; 
Qaux=[Qaux; Qaux2]; 

% plotting the amplified amplitudes of sines and cosines

aI=2*cos(beta);
aI=[aI; sqrt(2)*cos(alpha)];

aQ=2*sin(beta);
aQ=[aQ; sqrt(2)*sin(alpha)];

figure,stem([dopplerfreqs; fm ],aI,'k')
xlabel('Generator frequencies, I-rail (Hz)') 
ylabel('Magnitude')
figure,stem([dopplerfreqs; fm ],aQ,'k')
xlabel('Generator frequencies, Q-rail (Hz)') 
ylabel('Magnitude')
figure,stem([dopplerfreqs; fm ],sqrt(aI.^2+aQ.^2),'k')
xlabel('Generator frequencies I&Q rails combined (Hz)') 
ylabel('Magnitude')

% complex envelope =============================================

I=sum(Iaux);      % sum all rays on in-phase rail
Q=sum(Qaux);      % sum all rays on quadrature rail
r=(I+j.*Q);       % complex envelope 
r=r*a;            % and normalize complex envelope 

% plot amplitude and phase =============================================
figure,plot(timeaxis,abs(r),'k')
xlabel('Time (s)')
ylabel('Magnitude of complex envelope')

figure,plot(timeaxis,20*log10(abs(r)),'k')
xlabel('Time (s)')
ylabel('Magnitude of complex envelope (dB)')

figure,plot(timeaxis,angle(r),'k')
xlabel('Time (s)')
ylabel('Modulo-\pi phase of complex envelope (rad)')

figure,plot(timeaxis,unwrap(angle(r)),'k')
xlabel('Time (s)')
ylabel('Absolute phase of complex envelope (rad)')

% plot normalized RF spectrum ===============================
spectrumr=fftshift((abs(fft(r,NFFT))).^2);
freqaxis=[0:NFFT-1]*fs/NFFT-fs/2;
figure,plot(freqaxis,spectrumr,'k')
xlabel('Doppler shift (Hz)')
ylabel('Frequency response (dB/max)')
title('Doppler spectrum')

% calculate autocorrelation function ===================================
correl=xcorr(r,'coeff');
correl=real(correl);
correlaxis=([0:length(correl)-1]-length(correl)/2+1)*ts;
correltheoretical=(besselj(0,2*pi*fm*correlaxis));

figure,plot(correlaxis,correl,'k:',correlaxis,real(correltheoretical),'k')
legend('Simulation','Theoretical')
xlabel('Time lag (s)')
ylabel('Correlation coefficient')

% calculae correlation between in-phase and quadrature sereis
xcorrel=xcorr(real(r),imag(r),'coeff');
correlaxis=([0:length(correl)-1]-length(correl)/2+1)*ts;
figure,plot(correlaxis,xcorrel,'k')
xlabel('Time lag (s)')
ylabel('Cross-correlation coefficient')
title('I v. Q cross-correlation coefficient') 

% Calculate distributions of in-phase and quadrature rails
[CDFxI,CDFyI]=fCDF(real(r));
[CDFxQ,CDFyQ]=fCDF(imag(r));
[CDFyTheoretical]=GaussianCDF(0,sigma,CDFxI);
figure,plot(CDFxI,CDFyTheoretical,'k',CDFxI,CDFyI,'+',CDFxQ,CDFyQ,'^')
xlabel('I/Q rail value in linear units')
ylabel('Probability the abscissa is not exceeded')
title('CDFs of in-phase and Quadrature components') 
legend('Gaussian','I-Rail','Q-Rail')


% Calculate CDF of magnitude =======================================
[CDFx,CDFy]=fCDF(abs(r));

[CDFyTH]=RayleighCDF(sigma,CDFx);
figure,plot(CDFx,CDFy,'k:',CDFx,CDFyTH,'k')
xlabel('Magnitude of complex signal')
ylabel('Probability the abscissa is not exceeded')
title('Simulated and Rayleigh theoretical CDFs') 
legend('Simulated','Rayleigh')

% Calculate histogram of phase =======================================

[histy,histx]=hist(angle(r),20);
histy=histy/length(angle(r));     % normalize histigram for total probability equal one
figure,bar(histx,histy,'k')
xlabel('Phase (rad)')
ylabel('Probability')
title('Phase histogram')

