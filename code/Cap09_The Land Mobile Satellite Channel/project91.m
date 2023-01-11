%
% project91
%
% =======================================================================
clear
close all
clc
% =======================================================================
fMHz=1540;         % frequency in MHz
LFrame=1;          % min duration of one state in m 

RouteLength=100;   % Simulated route length in m
lambdac=300/fMHz;  % wavelength in m
F=6;               % Sampling fraction of wavelength
% =======================================================================
sigmaRayl=0.2;     % Sigma Rayleigh 
sigmaRice=0.15     % Sigma Rice
aR=1;            % a Rice             
ds=lambdac/F;    % sampling spacing

InterpRate=round(LFrame/ds);
ds=LFrame/InterpRate  % change slightly ds so that LFrame is a multiple of the new ds

% ======================================================================
% 
 fs=1/ds;       % sampling frequency in cycles/m 
 
% =============== Butterworth filter parameters =========================
Wp=0.09;        
Ws=0.16;
Rp=3; %dB
Rs=50; %dB

% =======================================================================
% Markov channel parameters
 P=[0.95 0.05
    0.1  0.9];

% ========================================================================

CurrentState=1;
StateSeries=[CurrentState];
d_axis=[0];
SigmaSeries=[sigmaRice];        % because CurrentState = 1

NoDraws=round(RouteLength/LFrame);

for ii=1:NoDraws
    drawState=rand(1,1);
    if CurrentState ==1,
        if  drawState<=P(1,1),
            StateSeries=[StateSeries; 1];
            SigmaSeries=[SigmaSeries; sigmaRice];
            d_axis=[d_axis; ii];
            CurrentState=1;
        else
            StateSeries=[StateSeries; 0];
            d_axis=[d_axis; ii];
            SigmaSeries=[SigmaSeries; sigmaRayl];
            CurrentState=0;
        end
    else
        if  drawState<=P(2,1),
            StateSeries=[StateSeries; 1];
            SigmaSeries=[SigmaSeries; sigmaRice];
            d_axis=[d_axis; ii];
            CurrentState=1;
        else
            StateSeries=[StateSeries; 0];
            SigmaSeries=[SigmaSeries; sigmaRayl];
            d_axis=[d_axis; ii];
            CurrentState=0;
        end
    end
end
d_axis=d_axis*LFrame;            % convert axis from samples to meters

% figure,plot(d_axis,StateSeries)
% aa=axis;
% axis([aa(1) aa(2) -0.5 1.5])

% figure,plot(d_axis,SigmaSeries)

% =======================================================================

% Change sampling rate to fraction F of the wavelength

InterpStateSeries=[];
InterpSigmaSeries=[];
for ii=1:length(StateSeries)
    if StateSeries(ii)==0, 
        InterpStateSeries=[InterpStateSeries; StateSeries(ii); zeros(InterpRate-1,1)];
        InterpSigmaSeries=[InterpSigmaSeries; ones(InterpRate,1)*sigmaRayl];
    else
        InterpStateSeries=[InterpStateSeries; StateSeries(ii); ones(InterpRate-1,1)];
        InterpSigmaSeries=[InterpSigmaSeries; ones(InterpRate,1)*sigmaRice];
    end
end

d_axisInterp=[0:length(InterpStateSeries)-1]*ds;
figure,plot(d_axisInterp,InterpStateSeries,'k') 
aa=axis;
axis([aa(1) aa(2) -0.5 1.5])
xlabel('Traveled distance (m)')
ylabel('State series')

figure,plot(d_axisInterp,InterpSigmaSeries,'k') 
xlabel('Traveled distance (m)')
ylabel('Sigma. Multipath parameter')
aa=axis;
axis([aa(1) aa(2) 0.0 0.5])


% Running mean low-pass filtering ===================================================

lengthwindow=InterpRate;       % Samples per state frame
averagingwindow=ones(lengthwindow,1)/lengthwindow;
FiltInterpStateSeries=conv(averagingwindow,InterpStateSeries);
FiltInterpSigmaSeries=conv(averagingwindow,InterpSigmaSeries);

FiltInterpStateSeries=FiltInterpStateSeries(lengthwindow:end);   % discard samples after convolution 
FiltInterpSigmaSeries=FiltInterpSigmaSeries(lengthwindow:end);   % discard samples after convolution 

figure, plot(d_axisInterp,FiltInterpStateSeries,'k')
aa=axis;
axis([aa(1) aa(2) -0.5 1.5])
xlabel('Traveled distance (m)')
ylabel('State series')

figure, plot(d_axisInterp,FiltInterpSigmaSeries,'k')
xlabel('Traveled distance (m)')
ylabel('Sigma. Multipath parameter')

% =======================================================================
% GENERATE FAST VARIATIONS

r=rayleigh(FiltInterpSigmaSeries,length(FiltInterpStateSeries));

% ===============  Filtered Rayliegh series ========================

[rFilt, B, A]=filtersignal(r,Wp,Ws,Rp,Rs);
[H,fre]=freqz(B,A,512,fs);      % For computing the filter's gain 
% Calculation of filter gain
     [h,T]=impz(B,A);
     gainF=sqrt(sum(h.^2));
     rFilt=rFilt/gainF;
% ........go on .........

figure,plot(d_axisInterp,20*log10(abs(rFilt)),'k')
aa=axis;
axis([aa(1) aa(2) -50 10])
xlabel('Traveled distance (m)')
ylabel('Multipath signal level (dB)')

% =============== Filter resposne ========================================

[H,fre]=freqz(B,A,512,fs);

figure, plot(fre,20*log10(abs(H)),'k')
xlabel('Frequency (Cycles/m)')
ylabel('Magnitude of filter response (dB)')
auxx=axis;
axis([0 fs/2 auxx(3) 0])
grid

% ========================================================================
TotalSeries=abs(FiltInterpStateSeries+r);
figure,plot(d_axisInterp,TotalSeries,'k')
aa=axis;
axis([aa(1) aa(2) -0.5 1.5])
xlabel('Traveled distance (m)')
ylabel('Overall signal (lin. units)')

figure,plot(d_axisInterp,20*log10(TotalSeries),'k')
xlabel('Traveled distance (m)')
ylabel('Overall signal level (dB)')

 
