%
% project103
%
% =======================================================================
% Initialize 
% =======================================================================
clear
close all
clc
% =======================================================================
% basic inputs 
% =======================================================================
fc=2000;         % MHz  Carrier frequency
F=8;             % sampling rate: fraction of wave length
V=10;            %  m/s MS1 speed
Nsamples=500;   % Number of samples

% =======================================================================
% indirect parameters 
% =======================================================================
lambdac=300/fc;    % m wavelength
Dx=lambdac/F;      % m sampling spacing
ts=Dx/V;          % s time sampling interval
fs=1/ts;           % Hz sampling frequency
kc=2*pi/lambdac;   % propagation constant
fm=V/lambdac                % max Doppler shift
timeaxis=ts.*[0:Nsamples-1];
N_tx=3
N_rx=3

% =====================================================================

load RMIMO

% ======================================================================
% Cholesky decomposition
% ======================================================================
C=chol(RMIMO);
C=C';                 % Pick lower-triangular matrix

% ======================================================================
% Butterworth Doppler filter parameters 
% ======================================================================
Wp=0.5*2*fm/fs;        
Ws=1*2*fm/fs;
Rp=3; %dB
Rs=40; %dB

% =======================================================================  
% Filter Rayleigh series 
% =======================================================================

rsim=cell(N_tx,N_rx);
rsimfilt=cell(N_tx,N_rx);
for ii=1:N_tx
    for jj=1:N_rx
        rsim{ii,jj}=sqrt(0.5)*(randn(Nsamples,1)+j.*randn(Nsamples,1));
        % Multiply by sqrt(0.5) to have 
        % unit variance complex Gaussian process
        [auxx, B, A]=filtersignal(rsim{ii,jj},Wp,Ws,Rp,Rs);        
        [h,T]=impz(B,A);
        gainF=sqrt(sum(h.^2));
        rsimfilt{ii,jj}=auxx/gainF; 
    end
end

% ======================================================================= 
% Plot Doppler filter resposne 
% =======================================================================
% Calculation of filter gain

[H,fre]=freqz(B,A,512,fs);

figure, plot(fre,20*log10(abs(H)))
xlabel('Frequency (Hz)')
ylabel('Magnitude of filter response (dB)')
title('Doppler filter')
auxx=axis;
axis([0 fs/2 auxx(3) 0])
grid



% convert cell in vector

Hsimfilt=zeros(Nsamples,N_tx*N_rx);

col=0;
figure,hold
for ii=1:N_tx
    for jj=1:N_rx
        col=col+1;
        Hsimfilt(:,col)=rsimfilt{ii,jj}; 
        plot(timeaxis,20*log10(abs(Hsimfilt(:,col))))
    end
end
xlabel('Time (s)')
ylabel('Normalized, Doppler shaped signals (dB)')


Hsimfilt=Hsimfilt';        % need to traspose to multiply by C
Hsimfiltcor=C*Hsimfilt;

figure,plot(timeaxis,20*log10(abs(Hsimfiltcor)))
xlabel('Time (s)')
ylabel('Normalized, Doppler shaped MIMO signals (dB)')


