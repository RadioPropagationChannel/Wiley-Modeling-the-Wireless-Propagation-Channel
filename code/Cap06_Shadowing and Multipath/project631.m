%
% project631
%
% Initialization ========================================================

clear all, close all

% Common parameters ======================================================

Nsamples=100000;         % approx.  number of samples to be generated 
f= 2000e6;             % carrier frequency (Hz)
V=10;                  % MS speed (m/s)
Lcorr=9.0;             % Slow variations. Correlation distance (m) 
F=4;                    % Sampling: Fraction of wavelength 
NFFT= 256;              % No. of FFT points

% Model parameters MP,S,M  (LOO) =======================================

MP=-10;                         % rms squared value of diffuse mulitpath dB
M=-8.0;                        % mean in dB
S= 3.0;                        % std in dB

sigma=sqrt(0.5*10^(MP/10));    % convert to linear units 

AoA=180;                       % angle of arrival of direct signal

% Secondary parameters ===================================================

lambdac=3e8/f;                  % wavelength (m)
kc=2*pi/lambdac;                % wave number  
fm=V/lambdac;                   % Max Doppler
ts=(lambdac/F)/V;               % sampling spacing (s)
fs=1/ts;                       % sampling freq. (Hz)

% Randomize seed using seconds field in computer clock ==================

clk = clock;
randn('state',sum(100*clk(6)));

% Gen.  slow variations. Lower rail in diagram 
% Uncorrelated samples are assumed to be spaced Lcorr (m) ===============

samplesLcorr=Lcorr/(lambdac/F);      % No. of samples within Lcorr
samplesLcorr=round(samplesLcorr);   % Integer number: modifies slightly Lcorr

Nslowsamples=Nsamples/samplesLcorr;     % No of slow var. samples 
Nslowsamples=round(Nslowsamples);      % Integer number 

Nsamples=Nslowsamples*samplesLcorr;    % New total number of samples


% Direct signal amplitude ==============================================

slow=randn(Nslowsamples,1);      % uncorrelated slow variations
A=(slow*S)+M;         % Normal distr.: Mean M and std S (ampl. dir. signal) 
a=10.^(A/20);         %convenrt to linear units (log-normal) 

% Interpolate slow direct signal variations =============================

x=[0:Nslowsamples-1]*Lcorr;     % axis in m (samples spaced Lcorr m)
x2=[0:Nslowsamples*samplesLcorr-1]*Lcorr/samplesLcorr;


G2=interp1(x,A,x2,'spline');         % Interpolated amplitude for A 
figure,plot(x2,G2,'r',x,A,'b')
xlabel('Traveled distance (m)')
ylabel('Signal level (dB/LOS)')
title('Slow signal variartions: before and after interpolation') 
legend('Slow variations after interpolation', 'Slow variations before interpolation', 'Location', 'Best');


G2ampl=10.^(G2/20);

% Direct signal's phase ===========================================

phi=kc.*[0:Nsamples-1]*lambdac./F*cosd(AoA);      % phase inccreases linearly with taveled distance 
phase_ini=2*pi*rand;            % initial phase (drawn randmoly)
phi=phi+phase_ini;               % overall phase of dir signal

%figure,plot([0:Nsamples-1]*lambdac./F,phi)

% overall complex direct signal ========================================

G2filt=G2ampl.*exp(j*phi);

% calculate Direct signal spectrum =======================================

SG2filt=fftshift(abs(fft(G2filt,NFFT)).^2);
freqaxis=[0:NFFT-1]*(fs/NFFT)-fs/2;
figure,plot(freqaxis,10*log10(SG2filt)-max(10*log10(SG2filt)))
title('Doppler spectrum, Direct signal')
xlabel('Frequency (Hz)')
ylabel('Normalized Doppler spectrum (dB)')
        
% Generate fast variations % upper rail =================================

I=randn(1,Nsamples);        % In-phase component  
Q=randn(1,Nsamples);        % Quadrature component 

% % JAKES filter ===========================================================

freqstep=5;
[filtertimedomain]=jakes(freqstep,fm,fs);
close(4:6)

% LOO time series generator  ============================================

% RAYLEIGH part. Upper rail in schematic diagram ========================

% Rayleigh before filtering (Jakes)
ray=(I+i*Q);

% Fiter and multiply by sigma. Upper rail ================================

ray_filt=sigma*conv(ray,filtertimedomain);

% Remove end samples after convolution ===================================

len=length(filtertimedomain);

ray_filt=ray_filt(len:end-len/2);
x2=x2(1:end-len/2);
G2filt=G2filt(1:end-len/2);
Nsamples=Nsamples-len/2;

figure,plot(x2,abs(ray_filt))
xlabel('Traveled distance (m)')
ylabel('Signal level (Linear units)')
title('Fast signal variartions after filtering') 

% Doppler spectrum of Rayleigh series ====================================

Sray_filt=fftshift(abs(fft(ray_filt,NFFT)).^2);
freqaxis=[0:NFFT-1]*(fs/NFFT)-fs/2;
figure,plot(freqaxis,10*log10(Sray_filt)-max(10*log10(Sray_filt)))
title('Doppler spectrum, Rayleigh series')
xlabel('Frequency (Hz)')
ylabel('Normalized Doppler spectrum (dB)')


% Building LOO series ==================================================

rLOO=ray_filt+G2filt;

figure,plot([0:Nsamples-1]*lambdac/F,20*log10(abs(rLOO)),'g',...
    [0:Nsamples-1]*lambdac/F,20*log10(abs(G2filt)),'r')
xlabel('Traveled distance (m)')
ylabel('Signal level (dB/LOS)')
title('Overall signal variations and slow signal variartions') 

% spectrum

SrLOO=fftshift(abs(fft(rLOO,NFFT)).^2);
freqaxis=[0:NFFT-1]*(fs/NFFT)-fs/2;
figure,plot(freqaxis,10*log10(SrLOO)-max(10*log10(SrLOO)))
title('Doppler spectrum, Loo series')
xlabel('Frequency (Hz)')
ylabel('Normalized Doppler spectrum (dB)')

% Calculate CDF ======================================================

% we remove samples at the end because interpolation process produces 
% unwanted peaks at the end of series

removesamples=round(Lcorr/(lambdac/F))

[xLOO,yLOO]=fCDF(abs(rLOO(1:end-removesamples)));
[xLN,yLN]=fCDF(abs(G2filt(1:end-removesamples)));
[xray,yray]=fCDF(abs(ray_filt(1:end-removesamples)));
figure,plot(20*log10(xLOO),yLOO,'g',20*log10(xLN),yLN,'r-.',20*log10(xray),yray,'k:')
title('CDF: Overall Loo, slow-lognormal and fast series') 
xlabel('Signal level (dB/LOS)')
ylabel('Porbability the abscissa is not exceeded')
legend('Loo','Lognormal','Rayleigh')

% Plotting the phase ================================================

figure, plot([0:Nsamples-1]*lambdac/F,unwrap(angle(rLOO)),'g',...
    [0:Nsamples-1]*lambdac/F,unwrap(angle(G2filt)),'r-.',...
    [0:Nsamples-1]*lambdac/F,unwrap(angle(ray_filt)),'k:');
title('Absolute phase (Rad.)');
xlabel('Traveled distance (m)');
ylabel('Phase (Rad.)')
legend('Loo','Direct ray','Diffuse multipath')


