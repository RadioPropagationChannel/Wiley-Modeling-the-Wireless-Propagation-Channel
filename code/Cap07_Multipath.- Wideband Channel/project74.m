%
% project 74     % cost 205 model for BAD URBAN case
%
% =======================================================================
clear
close all
clc
% =======================================================================
Nsamples=1000;         % number of samples to be generated 
fMHz=900;               % Freq. in MHz
V=10;                   % mobile speed in m/s
F=8;                   % Sampling: Fraction of wavelength 
NFFT= 128;              % No. of FFT points
freqstep=1;           % freq step in defining Doppler filters

% ========================================================================
lambdac=300/fMHz;              % wavelength
kc=2*pi/lambdac;               % wave number  
fm=V/lambdac;                  % Max Doppler
ts=(lambdac/F)/V;              % sampling spacing (s)
fs=1/ts;                       % sampling freq. (Hz)
timeaxis=[0:Nsamples-1]*ts;    % time axis
% =======================================================================
% BU Bad Urban
% Tap #  Del(us)  Pwr(dB)  Doppler
%  1     0          -7      CLASS
%  2     0.2        -3      CLASS
%  3     0.4        -1      CLASS 
%  4     0.8         0      GAUS1
%  5     1.6        -2      GAUS1
%  6     2.2        -6      GAUS2
%  7     3.2        -7      GAUS2
%  8     5.0        -1      GAUS2
%  9     6.0        -2      GAUS2
% 10     7.2        -7      GAUS2
% 11     8.2        -10     GAUS2
% 12    10.0        -15     GAUS2

TapPWR=[-7; -3; -1; 0; -2; -6; -7; -1; -2; -7; -10; -15]; 
TapNo=[1; 2; 3; 4; 5; 6; 7; 8; 9; 10; 11; 12];
TAU=[0; 0.2; 0.4; 0.8; 1.6; 2.2; 3.2; 5.0; 6.0; 7.2; 8.2; 10.0];    
TDL=[];
DOP=[];

figure,stem2D(TAU,TapPWR,-40)
auxy=axis;
axis([auxy(1)-1 auxy(2)+1 auxy(3) auxy(4)+10])
xlabel('Excess delay (us)')
ylabel('Relative level (dB)')
title('COST 207 Bad Urban TDL model')

% =======================================================================
% Tap #  Del(us)  Pwr(dB)  Doppler
%  1     0          -7      CLASS
% =======================================================================

I=randn(1,Nsamples);        % In-phase component  
Q=randn(1,Nsamples);        % Quadrature component 

[filtertimedomain]=jakes(freqstep,fm,fs);

ray=sqrt(0.5)*(I+j*Q);
ray=ray*10^(0.1*TapPWR(1));
ray_filt=conv(ray,filtertimedomain);

len=length(filtertimedomain);  % remove extra samples after convolution
ray_filt=ray_filt(len:length(ray_filt)-len-1);
timeaxis=[0:Nsamples-len-2]*ts;    % time axis

TDL=[TDL ray_filt'];
DOP=[DOP fftshift(abs(fft(ray_filt,NFFT)').^2)];

% =======================================================================
% Tap #  Del(us)  Pwr(dB)  Doppler
%  2     0.2        -3      CLASS
% =======================================================================

I=randn(1,Nsamples);        % In-phase component  
Q=randn(1,Nsamples);        % Quadrature component 

[filtertimedomain]=jakes(freqstep,fm,fs);

ray=sqrt(0.5)*(I+j*Q);
ray=ray*10^(0.1*TapPWR(2));
ray_filt=conv(ray,filtertimedomain);

len=length(filtertimedomain);  % remove extra samples after convolution
ray_filt=ray_filt(len:length(ray_filt)-len-1);

TDL=[TDL ray_filt'];
DOP=[DOP fftshift(abs(fft(ray_filt,NFFT)').^2)];

% =======================================================================
% Tap #  Del(us)  Pwr(dB)  Doppler
%  3     0.4        -1      CLASS 
% =======================================================================

I=randn(1,Nsamples);        % In-phase component  
Q=randn(1,Nsamples);        % Quadrature component 

[filtertimedomain]=jakes(freqstep,fm,fs);

ray=sqrt(0.5)*(I+j*Q);
ray=ray*10^(0.1*TapPWR(3));
ray_filt=conv(ray,filtertimedomain);

len=length(filtertimedomain);  % remove extra samples after convolution
ray_filt=ray_filt(len:length(ray_filt)-len-1);

TDL=[TDL ray_filt'];
DOP=[DOP fftshift(abs(fft(ray_filt,NFFT)').^2)];

% ========================================================================
% Tap #  Del(us)  Pwr(dB)  Doppler
%  4     0.8         0      GAUS1
% ========================================================================

I=randn(1,Nsamples);        % In-phase component  
Q=randn(1,Nsamples);        % Quadrature component 

fr11=-0.8*fm;
fr12=0.05*fm;
fr21=0.4*fm;
fr22=0.1*fm;
A=-10;             % attenuation of second component wrt first 

[filtertimedomain]=gaussCOST207(freqstep,fm,fs,fr11,fr12,fr21,fr22,A);

ray=sqrt(0.5)*(I+j*Q);
ray=ray*10^(0.1*TapPWR(4));
ray_filt=conv(ray,filtertimedomain);

len=length(filtertimedomain);  % remove extra samples after convolution
ray_filt=ray_filt(len:length(ray_filt)-len-1);

TDL=[TDL ray_filt'];
DOP=[DOP fftshift(abs(fft(ray_filt,NFFT)').^2)];

% ========================================================================
% Tap #  Del(us)  Pwr(dB)  Doppler
%  5     1.6        -2      GAUS1
% ========================================================================

I=randn(1,Nsamples);        % In-phase component  
Q=randn(1,Nsamples);        % Quadrature component 

fr11=-0.8*fm;
fr12=0.05*fm;
fr21=0.4*fm;
fr22=0.1*fm;
A=-10;             % attenuation of second component wrt first 

[filtertimedomain]=gaussCOST207(freqstep,fm,fs,fr11,fr12,fr21,fr22,A);

close, close, close 

ray=sqrt(0.5)*(I+j*Q);
ray=ray*10^(0.1*TapPWR(5));
ray_filt=conv(ray,filtertimedomain);

len=length(filtertimedomain);  % remove extra samples after convolution
ray_filt=ray_filt(len:length(ray_filt)-len-1);

TDL=[TDL ray_filt'];
DOP=[DOP fftshift(abs(fft(ray_filt,NFFT)').^2)];

% ========================================================================
% Tap #  Del(us)  Pwr(dB)  Doppler
%  6     2.2        -6      GAUS2
% ========================================================================

I=randn(1,Nsamples);        % In-phase component  
Q=randn(1,Nsamples);        % Quadrature component 

fr1=0.7*fm;
fr2=0.1*fm;
fr1=-0.4*fm;
fr2=0.15*fm;
A=-15;             % attenuation of second component wrt first 

[filtertimedomain]=gaussCOST207(freqstep,fm,fs,fr11,fr12,fr21,fr22,A);

ray=sqrt(0.5)*(I+j*Q);
ray=ray*10^(0.1*TapPWR(6));
ray_filt=conv(ray,filtertimedomain);

len=length(filtertimedomain);  % remove extra samples after convolution
ray_filt=ray_filt(len:length(ray_filt)-len-1);

TDL=[TDL ray_filt'];
DOP=[DOP fftshift(abs(fft(ray_filt,NFFT)').^2)];

% ========================================================================
% Tap #  Del(us)  Pwr(dB)  Doppler
%  7     3.2        -7      GAUS2
% ========================================================================

I=randn(1,Nsamples);        % In-phase component  
Q=randn(1,Nsamples);        % Quadrature component 

fr1=0.7*fm;
fr2=0.1*fm;
fr1=-0.4*fm;
fr2=0.15*fm;
A=-15;             % attenuation of second component wrt first 

[filtertimedomain]=gaussCOST207(freqstep,fm,fs,fr11,fr12,fr21,fr22,A);

close, close, close 

ray=sqrt(0.5)*(I+j*Q);
ray=ray*10^(0.1*TapPWR(7));
ray_filt=conv(ray,filtertimedomain);

len=length(filtertimedomain);  % remove extra samples after convolution
ray_filt=ray_filt(len:length(ray_filt)-len-1);

TDL=[TDL ray_filt'];
DOP=[DOP fftshift(abs(fft(ray_filt,NFFT)').^2)];

% ========================================================================
% Tap #  Del(us)  Pwr(dB)  Doppler
%  8     5.0        -1      GAUS2
% ========================================================================

I=randn(1,Nsamples);        % In-phase component  
Q=randn(1,Nsamples);        % Quadrature component 

fr1=0.7*fm;
fr2=0.1*fm;
fr1=-0.4*fm;
fr2=0.15*fm;
A=-15;             % attenuation of second component wrt first 

[filtertimedomain]=gaussCOST207(freqstep,fm,fs,fr11,fr12,fr21,fr22,A);

close, close, close 

ray=sqrt(0.5)*(I+j*Q);
ray=ray*10^(0.1*TapPWR(8));
ray_filt=conv(ray,filtertimedomain);

len=length(filtertimedomain);  % remove extra samples after convolution
ray_filt=ray_filt(len:length(ray_filt)-len-1); 

TDL=[TDL ray_filt'];
DOP=[DOP fftshift(abs(fft(ray_filt,NFFT)').^2)];

% ========================================================================
% Tap #  Del(us)  Pwr(dB)  Doppler
%  9     6.0        -2      GAUS2
% ========================================================================

I=randn(1,Nsamples);        % In-phase component  
Q=randn(1,Nsamples);        % Quadrature component 

fr1=0.7*fm;
fr2=0.1*fm;
fr1=-0.4*fm;
fr2=0.15*fm;
A=-15;             % attenuation of second component wrt first 

[filtertimedomain]=gaussCOST207(freqstep,fm,fs,fr11,fr12,fr21,fr22,A);

close, close, close 

ray=sqrt(0.5)*(I+j*Q);
ray=ray*10^(0.1*TapPWR(9));
ray_filt=conv(ray,filtertimedomain);

len=length(filtertimedomain);  % remove extra samples after convolution
ray_filt=ray_filt(len:length(ray_filt)-len-1);

TDL=[TDL ray_filt'];
DOP=[DOP fftshift(abs(fft(ray_filt,NFFT)').^2)];

% ========================================================================
% Tap #  Del(us)  Pwr(dB)  Doppler
% 10     7.2        -7      GAUS2
% ========================================================================

I=randn(1,Nsamples);        % In-phase component  
Q=randn(1,Nsamples);        % Quadrature component 

fr1=0.7*fm;
fr2=0.1*fm;
fr1=-0.4*fm;
fr2=0.15*fm;
A=-15;             % attenuation of second component wrt first 

[filtertimedomain]=gaussCOST207(freqstep,fm,fs,fr11,fr12,fr21,fr22,A);

close, close, close 

ray=sqrt(0.5)*(I+j*Q);
ray=ray*10^(0.1*TapPWR(10));
ray_filt=conv(ray,filtertimedomain);

len=length(filtertimedomain);  % remove extra samples after convolution
ray_filt=ray_filt(len:length(ray_filt)-len-1);

TDL=[TDL ray_filt'];
DOP=[DOP fftshift(abs(fft(ray_filt,NFFT)').^2)];

% ========================================================================
% Tap #  Del(us)  Pwr(dB)  Doppler
% 11     8.2        -10     GAUS2
% ========================================================================

I=randn(1,Nsamples);        % In-phase component  
Q=randn(1,Nsamples);        % Quadrature component 

fr1=0.7*fm;
fr2=0.1*fm;
fr1=-0.4*fm;
fr2=0.15*fm;
A=-15;             % attenuation of second component wrt first 

[filtertimedomain]=gaussCOST207(freqstep,fm,fs,fr11,fr12,fr21,fr22,A);

close, close, close 

ray=sqrt(0.5)*(I+j*Q);
ray=ray*10^(0.1*TapPWR(11));
ray_filt=conv(ray,filtertimedomain);

len=length(filtertimedomain);  % remove extra samples after convolution
ray_filt=ray_filt(len:length(ray_filt)-len-1);

TDL=[TDL ray_filt'];
DOP=[DOP fftshift(abs(fft(ray_filt,NFFT)').^2)];

% ========================================================================
% Tap #  Del(us)  Pwr(dB)  Doppler
% 12    10.0        -15     GAUS2
% ========================================================================

I=randn(1,Nsamples);        % In-phase component  
Q=randn(1,Nsamples);        % Quadrature component 

fr1=0.7*fm;
fr2=0.1*fm;
fr1=-0.4*fm;
fr2=0.15*fm;
A=-15;             % attenuation of second component wrt first 

[filtertimedomain]=gaussCOST207(freqstep,fm,fs,fr11,fr12,fr21,fr22,A);

close, close, close 

ray=sqrt(0.5)*(I+j*Q);
ray=ray*10^(0.1*TapPWR(12));
ray_filt=conv(ray,filtertimedomain);

len=length(filtertimedomain);  % remove extra samples after convolution
ray_filt=ray_filt(len:length(ray_filt)-len-1);

TDL=[TDL ray_filt'];
DOP=[DOP fftshift(abs(fft(ray_filt,NFFT)').^2)];

% =======================================================================
% Plotting TDL time-series
% =======================================================================
otheraxis=ones(length(timeaxis),1);

figure, hold
for ii=1:12
    plot3(otheraxis.*TapNo(ii),timeaxis, abs(TDL(:,ii)),'k')
end
grid
view(3)
title('TDL time-series')
xlabel('Tap No.')
zlabel('Relative signal level (lin.units)')
ylabel('Time (s)')

figure, hold
for ii=1:12
    plot3(otheraxis.*TAU(ii),timeaxis, 20*log10(abs(TDL(:,ii))),'k')
end
grid
view(3)
title('TDL time-series')
xlabel('Tap No.')
zlabel('Relative signal level (dB)')
ylabel('Time (s)')

% =======================================================================
% Plotting TDL Doppler spectra 
% =======================================================================
freqaxis=[0:NFFT-1]*(fs/NFFT)-fs/2;
otheraxis2=ones(length(freqaxis),1);

figure, hold
for ii=1:12
    plot3(otheraxis2.*TapNo(ii),freqaxis, DOP(:,ii),'k')
end
grid
view(3)

title('TDL Doppler spectra')
xlabel('Tap No.')
zlabel('Doppler spectrum')
ylabel('Doppler frequency (Hz)')

