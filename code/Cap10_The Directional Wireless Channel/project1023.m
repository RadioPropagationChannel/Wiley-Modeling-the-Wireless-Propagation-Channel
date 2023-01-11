%
% project1023   % MIMO scenario. MS static. Large antenna separation
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
 
NSC=100;         % Number of scatterers
avPower=0;     % sigma^2  Raverage power

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nreal=100;      % Number of realizations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% =======================================================================
% indirect parameters 
% =======================================================================
lambdac=300/fc;    % m wavelength
kc=2*pi/lambdac;   % propagation constant
a=sqrt(10.^(avPower/10)/NSC)  % magnitude of echoes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
realaxis=1:Nreal;    % realization axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
delta_tx=lambdac;   % Distance between adjacent transmitter antennas (m).
epsilon=0;         % Angle between Y-axis and the transmitter antennas axis (deg).

BS_epaxis=(-(N_tx-1)/2:(N_tx-1)/2)*delta_tx;
BSxi=(BS_epaxis*sind(180-epsilon))+BSx;
BSyi=(BS_epaxis*cosd(180-epsilon))+BSy;
plot(BSxi,BSyi,'k^')

% Mobile geometry
N_rx=3;             % No. of MS antennas 
delta_rx=lambdac;   % Distance between adjacent mobile antennas (m).

MSyi=(-(N_rx-1)/2:(N_rx-1)/2)*delta_rx;

%===================================================
MSx=0;
MSy=0;
plot(MSx,MSy,'r')
%===================================================

plot(repmat(MSx,1,length(MSyi)),MSyi,'r.')
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
alpha=rand(NSC,1)*(maxalpha-minalpha)+minalpha;  % random draw of angles of arrival

SCx=D.*cosd(alpha);
SCy=D.*sind(alpha);

plot(SCx,SCy,'*')
xlabel('Distance (m)');
ylabel('Distance (m)');

% =======================================================================
% calculate distance matrix 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
distBSSC=zeros(N_tx,1,NSC);
for ii=1:N_tx
    distBSSC(ii,1,:)=sqrt((BSxi(ii)-SCx).^2+(BSyi(ii)-SCy).^2);
end
distBSSC=repmat(distBSSC,[1 N_rx 1]);

distSCMS=zeros(1,N_rx,NSC);
for ii=1:N_rx
    distSCMS(1,ii,:)=sqrt((SCx-MSx).^2+(SCy-MSyi(ii)).^2);
end
distSCMS=repmat(distSCMS,[N_tx 1 1]);

distBSSCMS=distBSSC+distSCMS;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

% =======================================================================
% calculate complex envelope 
% =======================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r=zeros(N_tx,N_rx,Nreal);
for ii=1:Nreal
    % Random phases for each transmitter, receiver and scatterer.
    phi=[];
    for ll=1:N_tx
        for jj=1:N_rx
            for kk=1:NSC
                phi(ll,jj,kk)=rand(1,1)*2*pi;
            end
        end
    end
  
    ray=a*exp(-j*(kc*distBSSCMS-phi));
    r(:,:,ii)=sum(ray,3);
end


cadColors=['b  ';'g  ';'r  ';'c  ';'m  ';'y  ';'k  ';'b--';'g--'];
figure,hold on
for ii=1:N_tx
    for jj=1:N_rx
        dataplot=zeros(1,Nreal);
        dataplot(1,:)=abs(r(ii,jj,:));
        plot(realaxis,dataplot,cadColors(N_tx*(ii-1)+jj,:));
    end
end
hold off
xlabel('Realization')
ylabel('Magnitude of complex envelope (dB)')
title('All transmitters and receivers')
legend('Tx1-Rx1','Tx1-Rx2','Tx1-Rx3','Tx2-Rx1','Tx2-Rx2','Tx2-Rx3',...
    'Tx3-Rx1','Tx3-Rx2','Tx3-Rx3','Location','Best');

H=r;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% calculate singular values

Neigens=min(N_tx,N_rx);
eigens=zeros(Neigens,Nreal);
for ii=1:Nreal
    eigens(:,ii)=svd(H(:,:,ii));
end

eigens=eigens.^2;  % before they were singular values, now eigenvalues

figure,plot(realaxis,10*log10(eigens),'k')
xlabel('Realization')
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

CSISO(1,:)=log2(1+snr.*abs(r(1,1,:)).^2);
CMIMO=log2(1+snr.*eigens./Neigens);

% figure,plot(realaxis,CMIMO,realaxis,sum(CMIMO),realaxis,CSISO,'.-')
figure,plot(realaxis,CMIMO,'k:',realaxis,sum(CMIMO),'k',realaxis,CSISO,'k.-')
xlabel('Realization')
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
                auxx=corrcoef(r(ii,jj,:),r(kk,mm,:));
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
        auxx=corrcoef(r(ii,1,:),r(kk,1,:));
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
        auxx=corrcoef(r(1,ii,:),r(1,kk,:));
        RMS(ii,kk)=auxx(1,2);
    end
end
RMS

% =======================================================================
RMIMOkron=kron(RBS,RMS)

abs(RMIMO)-abs(RMIMOkron)

abs(RMIMO-RMIMOkron)

