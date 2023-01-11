% =======================================================================
% project732    (includes scenario editor)
% =======================================================================
% Initialize ============================================================
clear
close all
clc
warning off 
% basic inputs ==========================================================

fc=2000;         % MHz  Carrier frequency
F=4;             % sampling rate: fraction of wave length
V=10;            % m/s MS1 speed 
NFFT=64;         % Number of points in FFT
Nsamples=400;    % Number of route samples 
avPower=-20;     % sigma^2  Raverage power
delaystep=1e-7   % delay discretization setep in s
step_f=0.01;     % Freq axis step MHz

% geometry inputs ========================================================

dBS=5000;     
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

%========================================================================
% axes
% =======================================================================
timeaxis=ts.*[0:Nsamples-1];
Doppleraxis=([0:NFFT-1]-NFFT/2)*(fs/(NFFT-1));
faxis=[1999:step_f:2001];    % Freq axis in  MHz
% DELAY AXIS DEPENDS ON MAX DELAY, SET LATER 

% ========================================================================
MS0=-V*timeaxis(end)/2;   % initial location of receiver (MS) x-coordinate

MSx=MS0+V.*timeaxis;    % MS route along x-axis
MSy=zeros(Nsamples,1)';  % MS route along x-axis (y=0)
plot(MSx,MSy,'k','LineWidth',5)

MINx=min(min(BSx,MSx))-1000;
MAXx=max(max(BSx,MSx))+1000;
MINy=min(min(BSy,MSy))-1000;
MAXy=max(max(BSy,MSy))+1000;
axis([MINx MAXx MINy MAXy])
plot([0 0],[MINy MAXy], 'k:')
plot([MINx MAXx],[0 0], 'k:')

%=========================================================================
% SCENARIO EDITOR 
% ========================================================================
% placing point-scatterers in propagation scenario.

[SCx,SCy] = getpts(fig);
NSC=length(SCx);
plot(SCx,SCy,'k+');
xlabel('Distance (m)')
ylabel('Distance (m)')

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

distBSMS1aux=sqrt((BSx-MSx).^2+(BSy-MSy).^2);   
distBSMS1=min(min(distBSMS1aux));               % Ref distance is min BSMS dist 

% a=(distBSMS1./distBSSC(:)).*(distBSMS1./distSCMS(:,1));
a=(distBSMS1./sqrt(distBSSC(:))).*(distBSMS1./sqrt(distSCMS(:,1)));  % <-----

DeltaPower=avPower-10*log10(sum(a.^2));
deltaa=10.^(DeltaPower/20);             % to achieve reference power
a=deltaa*a;

% =====================================================================
% Define time-varying complex magnitudes of point scatterer contributions 
% amplitudes remain constant while phases change

aa=zeros(NSC,Nsamples);     % create variable 

for k1=1:Nsamples           % scan route points
    for k2=1:NSC            % scan scatterers
        aa(k2,k1)=a(k2)*exp(-j*kc*distBSSCMS(k2,k1));  % time-varying phase
    end
end

% ======================================================================

distBSSCMS1=distBSSCMS-distBSMS1;     % set a new refernece for delays wrt to 
DelaysNormalized=distBSSCMS1/cc;      % arrival of direct ray, here assumed 
                                      % to be totally blocked 
% DelaysNormalized=distBSSCMS/cc; 
                                      
DelaysNormalized=round(DelaysNormalized/delaystep);  % quantify delays delaystep (s)

auxx=size(DelaysNormalized);

auxx2=max(max(DelaysNormalized))+1;    % to include 0 ns delay
ImpulseResponse=zeros(auxx(2),auxx2);     % Create delay profile with step delaystep (s)

for jj=1:auxx(2)           % scan route locations
    for ii=1:auxx(1)       % scan scatterers
        indexx=DelaysNormalized(ii,jj)+1;               %<--------- ????
        ImpulseResponse(jj,indexx)=ImpulseResponse(jj,indexx)+aa(ii,jj);     
                                    % put in corresponding delay bin 
                                    % complex amplitude of delta
    end
end

axisdelayprofile=[0:auxx2-1];   % axis in delaystep units
figure, hold
for ii=1:Nsamples
 stem(axisdelayprofile*delaystep*1e6,abs(ImpulseResponse(ii,:)))  
 % delays in us
 % accumulate deltas with time and delay on same plot
end
xlabel('Delay (\mus)')
ylabel('Relative signal level (lin.units)')
title('Absoultue value of time varying impulse response, h(\tau;t)')

figure,stem3(axisdelayprofile*delaystep*1e6,timeaxis,abs(ImpulseResponse))
% delays in us
view(3)
xlabel('Delay (\mus)')
ylabel('Time (s)')
zlabel('Relative signal level (lin.units)')
title('Absoultue value of time varying impulse response, h(\tau;t)')
% ========================================================================
% frequency response for faxis (Freq axis)
% ========================================================================

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
title('Time-varying frequency response')

figure;mesh(faxis,timeaxis,20*log10(abs(FreqResp)))
ylabel('Time (s)')
xlabel('Frequency (MHz)')
zlabel('Level (dB)')
title('Time-varying frequency response')

figure;plot(faxis,20*log10(abs(FreqResp(1,:))),'k')
xlabel('frequency (MHz)')
ylabel('level (dB)')
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
title('Time-varying impulse response. Magnitude')

figure;plot(tauaxis,abs(ImpResp(1,:)),'k')
xlabel('delay (s)')
zlabel('level (l.u.)')
title('Impulse response for the first route point. Magnitude')

% =======================================================================
% Compute scattering matrix
% =======================================================================

% ImpulseResponse(jj,kk)   jj is route ppoint and kk is delay bin

auxz=size(ImpulseResponse);
ScatMat=zeros(NFFT,auxz(2));

for kk=1:auxz(2)           % scan delay bins 
    ScatMat(:,kk)=fftshift(abs(fft(ImpulseResponse(:,kk),NFFT)).^2);
end
figure,surf(axisdelayprofile*delaystep*1e6,Doppleraxis,ScatMat)
 shading interp
 colormap(jet) 
xlabel('Delay (\mus)')
ylabel('Doppler (Hz)')
zlabel('Level (lin. units)')
title('Channel scattering matrix')


ScatMatdB=10*log10(ScatMat);
Floor=-40; 
[ww,zz]=find(ScatMatdB==-inf);
ScatMatdB(ww,zz)=Floor;

figure,surf(axisdelayprofile*delaystep*1e6,Doppleraxis,ScatMatdB)
 shading interp
 colormap(jet) 
xlabel('Delay (\mus)')
ylabel('Doppler (Hz)')
zlabel('Level (dB)')
title('Channel scattering matrix')

% =======================================================================
% Tapped delay line
% =======================================================================

apdp=abs(ImpulseResponse(1,:)).^2;
APDP=10*log10(apdp)
minAPDP=-60;
maxAPDP=max(APDP);
maxDel=max(axisdelayprofile*delaystep*1e6);
figure, stem2D(axisdelayprofile*delaystep*1e6,APDP,minAPDP) 
axis([-0.5 maxDel+0.5 minAPDP maxAPDP+10])

xlabel('Delay (\mus)')
ylabel('Relative signal level (dB)')
title('Averaged power delay profile')
 
