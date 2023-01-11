%
% project714    (includes scenario editor)
%
% Initialize ============================================================
clear
close all
clc
warning off MATLAB:divideByZero
% basic inputs ==========================================================

fc=2000;         % MHz  Carrier frequency
F=8;            % sampling rate: fraction of wave length
V=10;            %  m/s MS1 speed 
NFFT=1024;       % Number of points in FFT
Nsamples=200;    % Number of samples 
avPower=-20;     % sigma^2  Raverage power

% =======================================================================
% Frequency inputs

step_f=0.01;                % Freq axis step MHz
faxis=[1997.5:step_f:2002.5];    % Freq axis in  MHz

% geometry inputs ========================================================

dBS=1000;     
angleBS=130;
BSx=dBS*cosd(angleBS) % location of transmitter (BS) x-coordinate
BSy=dBS*sind(angleBS)  % location of transmitter (BS) y-coordinate

% locations of point scatterers =========================================

fig=figure;
plot(BSx,BSy,'k^'), hold on

% indirect parameters ===================================================

lambdac=300/fc;    % m wavelength
Dx=lambdac/F;      % m sampling spacing 
ts=Dx/V;           % s time sampling interval
fs=1/ts;           % Hz sampling frequency
kc=2*pi/lambdac;   % propagation constant
fm=V/lambdac       % max Doppler shift
cc=3e8;            % speed of light
timeaxis=ts.*[0:Nsamples-1];

% ========================================================================
MS0=-V*timeaxis(end)/2;   % initial location of receiver (MS) x-coordinate

MSx=MS0+V.*timeaxis;    % MS route along x-axis
MSy=zeros(Nsamples,1);  % MS route along x-axis (y=0)
plot(MSx,MSy,'k','LineWidth',5)

MINx=min(min(BSx,MSx))-1500;
MAXx=max(max(BSx,MSx))+1500;
MINy=min(min(BSy,MSy))-1500;
MAXy=max(max(BSy,MSy))+1500;
axis([MINx MAXx MINy MAXy])
plot([0 0],[MINy MAXy], 'k:')
plot([MINx MAXx],[0 0], 'k:')

% SCENARIO EDITOR =======================================================
% placing point-scatterers in propagation scenario.

[SCx,SCy] = getpts(fig);
NSC=length(SCx);
plot(SCx,SCy,'k+');

% =======================================================================
% calculate distance matrix 
% =======================================================================

distBSSC=sqrt((BSx-SCx).^2+(BSy-SCy).^2);

distBSSCext=repmat(distBSSC,1,Nsamples);

distSCMS=zeros(NSC,Nsamples);
for ii=1:Nsamples
    distSCMS(:,ii)=sqrt((SCx-MSx(ii)).^2+SCy.^2);
end

distBSSCMS=distBSSCext+distSCMS;

% ======================================================================

distBSMS1=sqrt((BSx-MSx(1)).^2+(BSy-MSy(1)).^2);

a=(distBSMS1./distBSSC(:)).*(1./distSCMS(:,1));
DeltaPower=avPower-10*log10(sum(a.^2));
deltaa=10.^(DeltaPower/20)
a=deltaa*a



% ========================================================================
% frequency response for faxis 

FreqResp=zeros(Nsamples,length(faxis));

for k1=1:Nsamples               % scan route points
    for k2=1:length(faxis)      % scan frequencies
        for k3=1:NSC
            wl2=300/faxis(k2);
            FreqResp(k1,k2)=[FreqResp(k1,k2) + a(k3)*exp(-j*(2*pi/wl2)*distBSSCMS(k3,k1))];
        end
    end
end


figure;mesh(faxis,timeaxis,abs(FreqResp))
ylabel('Time (s)')
xlabel('Frequency (MHz)')
zlabel('Level (l.u.)')
title('Time.varying frequency response')

figure;mesh(faxis,timeaxis,20*log10(abs(FreqResp)))
ylabel('Time (s)')
xlabel('Frequency (MHz)')
zlabel('Level (dB)')
title('Time-varying frequency response')

figure;plot(timeaxis,20*log10(abs(FreqResp(:,1))),'k')
xlabel('time (s)')
ylabel('level (dB)')
title('Complex envelope magnitude central frequency')

figure;plot(faxis,20*log10(abs(FreqResp(1,:))),'k')
xlabel('frequency (MHz)')
ylabel('level (dB)')
title('Frequency response for first route point')

figure;plot(faxis,angle(FreqResp(1,:)),'k')
xlabel('frequency (MHz)')
ylabel('phase (Rad)')
title('Frequency response for first route point')

%=======================================================================
% Impulse response through IFFT

ImpResp=zeros(Nsamples,length(faxis));
for k4=1:Nsamples
    ImpResp(k4,:)=ifft(FreqResp(k4,:));
end
taumax=1/(step_f.*1e6);
step_tau=taumax/(length(faxis)-1);
step_tau=taumax/(length(faxis));               %<--------------?????
tauaxis=[0:length(faxis)-1].*step_tau;

figure;mesh(tauaxis,MSx,abs(ImpResp))
xlabel('delay (s)')
ylabel('route point (m)')
zlabel('level (l.u.)')
title('Impulse response variable with time')

figure;plot(tauaxis,abs(ImpResp(1,:)),'k')
xlabel('delay (s)')
ylabel('level (l.u.)')
title('Impulse response for the first route point')

figure;plot(tauaxis,abs(ImpResp(1,:)),'k')
xlabel('delay (s)')
ylabel('level (l.u.)')
title('Impulse response for the first route point')
aaa=axis;
hold on
for ii=1:length(distBSSCMS(:,1)),
    plot([distBSSCMS(ii,1)/cc; distBSSCMS(ii,1)/cc],[0; aaa(4)],'k:')
end



