%
% project1022   % MIMO scenario. MS moving. Small antenna separation
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
Nsamples=100;   % Number of samples
NSC=100;         % Number of scatterers
avPower=0;     % sigma^2  Raverage power

% =======================================================================
% indirect parameters 
% =======================================================================
lambdac=300/fc;    % m wavelength
Dx=lambdac/F;      % m sampling spacing
ts=Dx/V;          % s time sampling interval
fs=1/ts;           % Hz sampling frequency
kc=2*pi/lambdac;   % propagation constant
a=sqrt(10.^(avPower/10)/NSC)  % magnitude of echoes
fm=V/lambdac                % max Doppler shift
timeaxis=ts.*[0:Nsamples-1];
% =======================================================================
% geometric inputs 
% =======================================================================

dBS=500;
angleBS=180;
BSx=dBS*cosd(angleBS) % location of transmitter (BS) x-coordinate
BSy=dBS*sind(angleBS)  % location of transmitter (BS) y-coordinate

fig=figure;    % Open scenario plot
hold on

% =======================================================================
% Transmitter geometry
% =======================================================================
N_tx=3;             % Transmitter antennas number.
delta_tx=lambdac/16;   % Distance between adjacent transmitter antennas (m).
epsilon=0;         % Angle between Y-axis and the transmitter antennas axis (deg).

BS_epaxis=(-(N_tx-1)/2:(N_tx-1)/2)*delta_tx;
BSxi=(BS_epaxis*sind(180-epsilon))+BSx;
BSyi=(BS_epaxis*cosd(180-epsilon))+BSy;
plot(BSxi,BSyi,'k^')

% Mobile geometry
N_rx=3;             % No. of MS antennas 
delta_rx=lambdac/16;   % Distance between adjacent mobile antennas (m).

MSyi=(-(N_rx-1)/2:(N_rx-1)/2)*delta_rx;

%===================================================

MS0=-V*timeaxis(end)/2;        % initial location of receiver (MS) x-coordinate

MSx=MS0+V.*timeaxis;  % MS route along x-axis
MSy=zeros(Nsamples,1);  % MS route along x-axis (y=0)
plot(MSx,MSy,'r')
%===================================================

plot(repmat(MS0,1,length(MSyi)),MSyi,'r.')
%===================================================
MINx=min(min([BSxi MSx]))-200;
MAXx=max(max([BSxi MSx]))+200;
MINy=min(min(min([BSyi MSy'])))-200;
MAXy=max(max(max([BSyi MSy'])))+200;
axis([MINx MAXx MINy MAXy])
axis equal

% locations of point scatterers =========================================

minalpha=0;
maxalpha=360;

D=200;                        % radius from origin
alpha=rand(NSC,1)*(maxalpha-minalpha)+minalpha;        % random draw of angles of arrival

SCx=D.*cosd(alpha);
SCy=D.*sind(alpha);

plot(SCx,SCy,'*')
xlabel('Distance (m)');
ylabel('Distance (m)');

% =======================================================================
% calculate distance matrix 

distBSSC=cell(N_tx,1);
distBSSCext=cell(N_tx,1);
for ii=1:N_tx
    distBSSC{ii}=sqrt((BSxi(ii)-SCx).^2+(BSyi(ii)-SCy).^2);;
    distBSSCext{ii}=repmat(distBSSC{ii},[1 Nsamples]);
end

dist_BSCMS=cell(1,N_rx);
for ii=1:N_rx
    distSCMS{ii}=sqrt((repmat(SCx,1,Nsamples)-repmat(MSx,NSC,1)).^2+(repmat(SCy,1,Nsamples)-MSyi(ii)).^2);
end

distBSSCMS=cell(N_tx,N_rx);
for ii=1:N_tx
    for jj=1:N_rx
        distBSSCMS{ii,jj}=distBSSCext{ii}+distSCMS{jj};
    end
end

% =======================================================================
% calculate complex envelope 
% =======================================================================
ray=cell(N_tx,N_rx);
r=cell(N_tx,N_rx);
figure,hold
for ii=1:N_tx
    for jj=1:N_rx
        ray{ii,jj}=a*exp(-j*kc*distBSSCMS{ii,jj});
        r{ii,jj}=sum(ray{ii,jj});
        plot(timeaxis,20*log10(abs(r{ii,jj})),'k')        
    end
end
xlabel('Time (s)')
ylabel('Magnitude of complex envelope (dB)')
title('All transmitters and receivers')

% =======================================================================
% convert cell in matrix
% =======================================================================
H=zeros(N_tx,N_rx,Nsamples);

for ii=1:N_tx
    for jj=1:N_rx
        H(ii,jj,:)=r{ii,jj}; 
    end
end

% calculate eigenvalues

Neigens=min(N_tx,N_rx);
eigens=zeros(Neigens,Nsamples);
for ii=1:Nsamples
    eigens(:,ii)=svd(H(:,:,ii));
end

eigens=eigens.^2;  % before they were singular values, now eigenvalues

figure,plot(timeaxis,10*log10(eigens),'k')
xlabel('Time (s)')
ylabel('Eigenvalues (dB)')


CDFx=[];
CDFy=[];
for ii=1:min(N_tx,N_rx)
    [x,y]=fCDF(eigens(ii,:));
    CDFx=[CDFx, x'];
    CDFy=[CDFy, y'];
end
figure,semilogy(10*log10(CDFx),CDFy)
xlabel('Eigenvalues (dB)')
ylabel('Probability the abscissa is not exceeded')

% ======================================================================
% claculate capacity time-series with equal power assingment to all models 

SNR=20;   % Signal to noise ratio in dB
snr=10^(0.1*SNR);

CSISO=log2(1+snr.*abs(r{1,1}).^2);
CMIMO=log2(1+snr.*eigens./Neigens);

figure,plot(timeaxis,CMIMO,'k:',timeaxis,sum(CMIMO),'k',timeaxis,CSISO,'k.-')
xlabel('Time (s)')
ylabel('Capacity (b/s/Hz)')
legend('MIMO channels ','Overall MIMO','SISO', 'Location', 'Best')



[xMIMO,yMIMO]=fCDF(sum(CMIMO))
[xSISO,ySISO]=fCDF(CSISO)

figure,semilogy(xMIMO,yMIMO,'k',xSISO,ySISO,'k.-')
xlabel('Capacity (b/s/Hz)')
ylabel('Probability the abscissa is not exceeded')
legend('MIMO','SISO', 'Location', 'Best')

% ======================================================================
% calculate RMIMO

RMIMO=zeros(N_tx*N_rx,N_tx*N_rx);
row=1;
col=1;
for ii=1:N_tx
    for jj=1:N_rx
        for kk=1:N_tx
            for mm=1:N_rx
                auxx=corrcoef(r{ii,jj},r{kk,mm});
                RMIMO(row,col)=auxx(1,2);
                col=col+1;
            end
        end
        row=row+1;
        col=1;
    end
end

% save RMIMO RMIMO

% ========================================================================
% BS side correlations
% ========================================================================

RBS=zeros(N_tx,N_tx);
for ii=1:N_tx
    for kk=1:N_tx
        auxx=corrcoef(r{ii,1},r{kk,1});
        RBS(ii,kk)=auxx(1,2);
    end
end
RBS
% ========================================================================
% MS side correlations
% ========================================================================

RMS=zeros(N_rx,N_rx);
for ii=1:N_rx
    for kk=1:N_rx
        auxx=corrcoef(r{1,ii},r{1,kk});
        RMS(ii,kk)=auxx(1,2);
    end
end
RMS

% =======================================================================
RMIMOkron=kron(RBS,RMS)

abs(RMIMO)-abs(RMIMOkron)

abs(RMIMO-RMIMOkron)

