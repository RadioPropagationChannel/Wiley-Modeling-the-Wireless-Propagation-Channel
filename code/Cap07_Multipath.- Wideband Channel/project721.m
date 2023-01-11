% =======================================================================
% project721    
% =======================================================================
% Initialize ============================================================
clear
close all
clc
warning off 
% basic inputs ==========================================================

fc=2000;         % MHz  Carrier frequency
F=8;            % sampling rate: fraction of wave length
V=10;            %  m/s MS1 speed 
NFFT=1024;       % Number of points in FFT
Nsamples=200;    % Number of samples 
avPower=-20;     % sigma^2  Raverage power

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
MSy=zeros(Nsamples,1)';  % MS route along x-axis (y=0)
plot(MSx,MSy,'k','LineWidth',5)

MINx=min(min(BSx,MSx))-500;
MAXx=max(max(BSx,MSx))+500;
MINy=min(min(BSy,MSy))-500;
MAXy=max(max(BSy,MSy))+500;
axis([MINx MAXx MINy MAXy])
plot([0 0],[MINy MAXy], 'k:')
plot([MINx MAXx],[0 0], 'k:')

%=========================================================================
% SCENARIO EDITOR 
% ========================================================================
% placing point-scatterers in propagation scenario.

SCx =[ 74.85; 44.55; -12.25; -42.55; -91.79; -42.55;... 
   -103.15; -57.70; -31.19; 10.46; 40.76; 139.23;...
   63.49; -16.04; -76.64; -106.94; -65.27; 6.68;  ...
   33.19; -27.40; 93.79;  78.64; 116.51; 82.42; ...
   55.91; -12.25; -122.09; -148.60; -141.02; -122.09; ...
   -4.68; 71.06; 139.23; 135.45; 108.93; 25.61;...
  -61.49; -137.23];

SCy =[ 60.28; 44.78; 44.78; 29.29; 80.93; 86.09; ...
   -22.34; -27.50; -58.48; -63.65; -42.99; -58.48; ...
  -94.63; -130.78; -94.63; 44.78; 122.24; 173.88; ...
  101.59; 153.22; -22.34; -110.12; 3.47; 122.24; ...
  163.55; 189.37; 111.91; 18.96; -58.48; -120.45; ...
  -151.43; -141.11; 34.46; 117.08; 189.37; 225.52; ...
  225.52; 168.72];

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

% set scatterer levels

a=(distBSMS1./distBSSC(:)).*(1./distSCMS(:,1));
DeltaPower=avPower-10*log10(sum(a.^2));
deltaa=10.^(DeltaPower/20);             % to achieve reference power
a=deltaa*a;

% =====================================================================
% Define time-varying complex magnitudes of point scatterer contributions 
% amplitudes remain constant while phases change

aa=zeros(Nsamples,NSC);     % create variable 

for k1=1:Nsamples           % scan route points
    for k2=1:NSC            % scan scatterers
        aa(k1,k2)=a(k2)*exp(-j*kc*distBSSCMS(k2,k1));  % time-varying phase
    end
end
 
% ======================================================================

distBSSCMS1=distBSSCMS-distBSMS1;     % set a new refernece for delays wrt to 
DelaysNormalized=distBSSCMS1/cc;      % arrival of direct ray, here assumed 
                                      % to be totally blocked 

DelaysNormalized=round(DelaysNormalized*1e9);  % convert quantify delays to 1 ns
%figure,mesh(DelaysNormalized)
auxx=size(DelaysNormalized);

auxx2=max(max(DelaysNormalized))+1;    % to include 0 ns delay
DelayProfile=zeros(auxx(2),auxx2);     % Create delay profile with step of 1 ns

for jj=1:auxx(2)           % scan route locations
    for ii=1:auxx(1)       % scan scatterers
        indexx=DelaysNormalized(ii,jj)+1;               %<--------- ????
        DelayProfile(jj,indexx)=DelayProfile(jj,indexx)+aa(jj,ii);     
                                    % put in corresponding delay bin 
                                    % complex amplitude of delta
    end
end

axisdelayprofile=[0:auxx2-1];   % axis in ns
figure, hold
for ii=1:Nsamples
 stem(axisdelayprofile,abs(DelayProfile(ii,:)))  % accumulate deltas with time and delay on same plot
end
xlabel('Excess delays (ns)')
ylabel('Several ideal channel responses')

% ======================================================================= 
% Parameters of ideal PDP

% PDP parameters: D, S 
%========================================================================
D=sum(abs(DelayProfile(1,:)).*axisdelayprofile)/sum(abs(DelayProfile(1,:)))  
% 1st route point
S=sqrt(sum(abs(DelayProfile(1,:)).*(axisdelayprofile-D).^2)/sum(abs(DelayProfile(1,:))))
%=======================================================================
% Simulate channel sounding
%=======================================================================
% build sounding pulse

SequenceLength=511;      % length of sounding sequence
tau_sampl_spacing=1;     % in ns
tau_chip=100;            % in ns   1/Chip rate of sounding sequence
SamplesInImpulse=round(tau_chip/tau_sampl_spacing);
tau_chip_mod=SamplesInImpulse*tau_sampl_spacing;

SlopePulse=SequenceLength/(SamplesInImpulse+1);

SoundImpulse=[0:SamplesInImpulse]*SlopePulse;

SoundImpulseAux=[SamplesInImpulse-1:-1:0]*SlopePulse;
SoundImpulse=[SoundImpulse SoundImpulseAux];

% figure, plot([0:length(SoundImpulse)-1],SoundImpulse)

% Normilize sound impulse to unit energy
ESoundImpulse=sum(SoundImpulse.^2);
SoundImpulse=SoundImpulse/sqrt(ESoundImpulse);
figure, plot([0:length(SoundImpulse)-1],SoundImpulse)
xlabel('Samples')
ylabel('Normalized sounding pulse')

figure,hold
pdp=[];
for ii=1:Nsamples,
    pdpaux=conv(SoundImpulse,DelayProfile(ii,:));
    auxx3=length(pdpaux);
    plot([0:auxx3-1]-tau_chip,abs(pdpaux)) 
    pdp=[pdp, pdpaux'];
    %pause
end
xlabel('Delay (ns)')
ylabel('Convolutions of sounding pulse and series of imp. responses')

% figure, plot([0:auxx3-1]-tau_chip,10*log10(abs(pdpaux)))


figure,surf(timeaxis, [0:auxx3-1]-tau_chip,abs(pdp))    
shading interp, colormap(jet) 
zlabel('Instantaneous PDP series (lin. units)')
ylabel('Excess delay (ns)')
xlabel('Time (s)')

figure,surf(timeaxis, [0:auxx3-1]-tau_chip,10*log10(abs(pdp)))
shading interp, colormap(jet) 
zlabel('Instantaneous PDP series (dB)')
ylabel('Excess delay (ns)')
xlabel('Time (s)')

% =======================================================================
% Compute averaged power delay profile APDP from instantaneous PDP
% =======================================================================
apdp=sum(abs(pdp'));    % trapose pdp since sum operates row-wise
apdp=apdp/Nsamples;      % 

figure, plot([0:auxx3-1]-tau_chip,apdp)
ylabel('Average PDP (lin.units)')
xlabel('Excess delay (ns)')

figure, plot([0:auxx3-1]-tau_chip,10*log10(apdp))
ylabel('Average PDP (dB)')
xlabel('Excess delay (ns)')
% =======================================================================
% Compute parameters of APDP
% =======================================================================
D=sum(apdp.*[0:auxx3-1])/sum(apdp)
S=sqrt(sum(apdp.*([0:auxx3-1]-D).^2)/sum(apdp))
