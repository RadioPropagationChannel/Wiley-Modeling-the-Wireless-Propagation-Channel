%
% project322  (Power control speed 56 m/s)
%
% Initialization ========================================================

clear all, close all

% Common parameters ======================================================

Nsamples=100050;         % approx.  number of samples to be generated 
f= 900e6;             % carrier frequency (Hz)
V=56;                  % MS speed (m/s) 
Lcorr=9.0;             % Slow variations. Correlation distance (m) 
F=500;                    % Sampling: Fraction of wavelength 
NFFT= 256;              % No. of FFT points

% Model parameters S,M  (SUZUKI) =======================================

	M=-13;               % Mean value of sigma in dB
	S= 6.0;              % std of sigma in dB
    
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


% Variations of sigma ==============================================

slow=randn(Nslowsamples+1,1);      % uncorrelated slow variations (one more sample to avoid problems in interoplation
A=(slow*S)+M;         % Normal distr.: Mean M and std S  
a=10.^(A/20);         %convenrt to linear units (log-normal) 

% Interpolate slow variations of sigma =============================

% x=[0:Nslowsamples-1]*Lcorr;     % axis in m (samples spaced Lcorr m)
% x2=[0:Nslowsamples*samplesLcorr-1]*Lcorr/samplesLcorr;

x=[0:Nslowsamples]*Lcorr;     % axis in m (samples spaced Lcorr m)
x2=[0:Nslowsamples*samplesLcorr]*Lcorr/samplesLcorr;

G2=interp1(x,A,x2,'spline');         % Interpolated amplitude for A 

% we use one extra sample in x to avoid interpolation problems at end
% points. Now we discard it and the corresponding interpolated samples 
x=x(1:Nslowsamples);
A=A(1:Nslowsamples);
x2=x2(1:Nslowsamples*samplesLcorr+1);
G2=G2(1:Nslowsamples*samplesLcorr+1);

% new corrected number of samples
Nsamples=Nslowsamples*samplesLcorr+1

figure,plot(x2,G2,'r',x,A,'b')
xlabel('Traveled distance (m)')
ylabel('Signal level (dB/LOS)')
title('Slow signal variartions: before and after interpolation') 
legend('Slow variations after interpolation', 'Slow variations before interpolation', 'Location', 'Best');

G2filt=10.^(G2/20);


% Generate fast variations % upper rail =================================

I=randn(1,Nsamples);        % In-phase component  
Q=randn(1,Nsamples);        % Quadrature component 

% JAKES filter ===========================================================

freqstep=5;
[filtertimedomain]=jakeswofigs(freqstep,fm,fs);


% Suzuki time series generator  ============================================

% RAYLEIGH part. Upper rail in schematic diagram ========================

% Rayleigh before filtering (Jakes)

ray=(I+j*Q);


% Fiter. Upper rail ================================

ray_filt=conv(ray,filtertimedomain);


% Remove end samples after convolution ===================================
len=length(filtertimedomain);
ray_filt=ray_filt(len:end-len/2);
x2=x2(1:end-len/2);
G2filt=G2filt(1:end-len/2);
Nsamples=Nsamples-len/2;

% Building Suzuki series ==================================================

rSUZ=ray_filt.*G2filt;

figure,plot([0:Nsamples-1]*lambdac/F,20*log10(abs(rSUZ)),'g',...
    [0:Nsamples-1]*lambdac/F,20*log10(abs(G2filt)),'r')
xlabel('Traveled distance (m)')
ylabel('Signal level (dB/LOS)')
title('Overall signal variations and slow signal variartions') 

taxisSUZ=[0:Nsamples-1]*ts;
figure,plot(taxisSUZ,20*log10(abs(rSUZ)),'g',...
    taxisSUZ,20*log10(abs(G2filt)),'r')
xlabel('time (s)');
ylabel('Signal level (dB/LOS)')
title('Overall signal variations and slow signal variartions') 

% Calculate CDF ======================================================

[xSUZ,ySUZ]=fCDF(abs(rSUZ));
[xLN,yLN]=fCDF(abs(G2filt));
[xray,yray]=fCDF(abs(ray_filt));

%==========================================================================
% Mobile side.
% Sliding window.
wint=1/1500;               % Window size in s.
winLen=round(wint/ts);     % Window size in samples.

threshold=median(20*log10(abs(rSUZ)));  % power threshold in dB.

numWin=ceil(Nsamples/winLen);
control=0;
rSUZpc=[];
sigControl=[];
for win=1:numWin
    lowWin=(win-1)*winLen+1;
    upWin=win*winLen;
    if upWin>Nsamples
        upWin=Nsamples;
    end
    
    signal=20*log10(abs(rSUZ(lowWin:upWin)))+control;
    rSUZpc=[rSUZpc signal];
    sigControl=[sigControl control*ones(1,length(lowWin:upWin))];
    
    winMean=mean(signal);
    if winMean>threshold
        control=control-1;
    elseif winMean<threshold
        control=control+1;
    end
end

figure,plot(taxisSUZ,20*log10(abs(rSUZ)),'g',taxisSUZ,sigControl,'r',taxisSUZ,rSUZpc)
xlabel('time (s)');
ylabel('Normalized received power (dB)')
legend('Received signal without PC','Control signal','Power Controlled signal','Location','SouthWest');

% Distribution of the power controlled signal.
[xSUZpc,ySUZpc]=fCDF(rSUZpc);
[xcontrol,ycontrol]=fCDF(sigControl);

figure,plot(20*log10(xSUZ),ySUZ,'g',xSUZpc,ySUZpc),grid
title('CDFs of received signal with and without Power Control') 
xlabel('Signal level (dB/LOS)')
ylabel('Probability the abscissa is not exceeded')
legend('Received signal','Power controlled signal','Location','NorthWest')
axis_fig=axis;
axis([axis_fig(1) axis_fig(2) 0 1])

stdSUZ=std(20*log10(rSUZ)) 
stdSUZpc=std(rSUZpc)