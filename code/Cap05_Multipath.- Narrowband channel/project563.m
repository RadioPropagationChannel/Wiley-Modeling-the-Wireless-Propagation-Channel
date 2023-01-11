%
% p5oject562     % Epsiolon 15 degree
%
%=============================== RESET ==================================
clear
close all
clc
% basic inputs ==========================================================

fc=2000;                  % MHz  Carrier frequency
F=16;                     % sampling rate: fraction of wave length
V=10;                     %  m/s MS1 speed 
NFFT=128;                 % Number of points in FFT
Nsamples=1000             % Number of samples 
NSC = 100                 % number of point scatterers
avPower=-20;              % sigma^2  Average power 

lambdac=300/fc;    % m wavelength

% ================== BS diversity parameters ===========================

epsilon =15;        % (degree) angle between baseline between antennas and link direction
P = 1;               % separation between BS antennas in number of wavelengths 
Delta = P*lambdac;   % m actual separation between antennas
Nsims=50;            % number of simulations

% ======================== geometry inputs ===========================

dBS=20000;              % distance BS origin (m) beginning mobile route
BSx=dBS*cosd(90)  % location of transmitter (BS) x-coordinate
BSy=dBS*sind(90)  % location of transmitter (BS) y-coordinate

% Create a ring of point scatterers =========================================

D=100;                        % radius from origin
alpha=rand(NSC,1)*360;        % random draw of angles of arrival

SCx=D.*cosd(alpha);
SCy=D.*sind(alpha);

figure,plot(SCx,SCy,'k*', BSx,BSy,'k^'), hold on

% indirect parameters ===================================================

Dx=lambdac/F;      % m sampling spacing 
ts=Dx/V;          % s time sampling interval
fs=1/ts;           % Hz sampling frequency
kc=2*pi/lambdac;   % propagation constant

a=sqrt(10.^(avPower/10)/NSC)  % magnitude of echoes
sigma=sqrt(0.5*10.^(avPower/10))     % Rayleigh parameter

fm=V/lambdac;                % max Doppler shift

timeaxis=ts.*[0:Nsamples-1];

MS0=-V*timeaxis(end)/2;        % initial location of receiver (MS) x-coordinate

MSx=MS0+V.*timeaxis;  % MS route along x-axis
MSy=zeros(Nsamples);  % MS route along x-axis (y=0)
plot(MSx,MSy,'k')

MINx=min(min(BSx, SCx))-100;
MAXx=max(max(BSx, SCx))+100;
MINy=min(min(min(BSy, SCy)))-100;
MAXy=max(max(max(BSy, SCy)))+100;
axis([MINx MAXx MINy MAXy])
axis equal
xlabel('Propagation scenario. Dimensions in meters')
ylabel('Propagation scenario. Dimensions in meters')

m1=[];
m2=[];
divr=[];

for ii=0:Nsims

    % ============ create antenna separations =================
    
    BSx1=BSx+(ii*Delta/2)*sind(epsilon);
    BSx2=BSx-(ii*Delta/2)*sind(epsilon);
    
    BSy1=BSy+(ii*Delta/2)*cosd(epsilon);
    BSy2=BSy-(ii*Delta/2)*cosd(epsilon);
    
    % calculate distance matrix =============================================
    distBSSC1=sqrt((BSx1-SCx).^2+(BSy1-SCy).^2);
    distBSSC2=sqrt((BSx2-SCx).^2+(BSy2-SCy).^2);

    distBSSCext1=repmat(distBSSC1,1,Nsamples);
    distBSSCext2=repmat(distBSSC2,1,Nsamples);

    distSCMS=zeros(NSC,Nsamples);
    for ii=1:Nsamples
        distSCMS(:,ii)=sqrt((SCx-MSx(ii)).^2+SCy.^2);
    end

    distBSSCMS1=distBSSCext1+distSCMS;
    distBSSCMS2=distBSSCext2+distSCMS;

    % calculate complex envelope ===========================================
    ray1=a*exp(-j*kc*distBSSCMS1);
    ray2=a*exp(-j*kc*distBSSCMS2);

    r1=sum(ray1); m1aux=abs(r1);
    r2=sum(ray2); m2aux=abs(r2);
    
    divraux=max(m1aux,m2aux);

    m1=[m1 m1aux'];
    m2=[m2 m2aux'];
    divr=[divr divraux'];
    
    % =====================================================================

end  % of for loop


% ======= calculating diversity gain using CDF ====================


CDFxx1=[];CDFxx2=[];CDFxx3=[];
CDFyy1=[];CDFyy2=[];CDFyy3=[];

for ii=0:Nsims
    [CDFx1,CDFy1]=fCDF(20*log10(m1(:,ii+1)));
    [CDFx2,CDFy2]=fCDF(20*log10(m2(:,ii+1)));
    [CDFx3,CDFy3]=fCDF(20*log10(divr(:,ii+1)));

    CDFxx1=[CDFxx1 CDFx1']; CDFxx2=[CDFxx2 CDFx2']; CDFxx3=[CDFxx3 CDFx3'];
    CDFyy1=[CDFyy1 CDFy1']; CDFyy2=[CDFyy2 CDFy2']; CDFyy3=[CDFyy3 CDFy3'];

end



figure,semilogy(CDFxx1(:,1),CDFyy1(:,1), CDFxx3,CDFyy3)
xlabel('Relative received signal level (dB)')
ylabel('Probability the abscissa is not exceeded')

% ========== evolution of cross-correlation coeff ==============

rho=[];
for ii=1:Nsims
    rhoaux=xcorr(m1(:,ii)-mean(m1(:,ii)),m2(:,ii)-mean(m2(:,ii)),'coeff');
    rho=[rho rhoaux(length(m1))];
end


% Theoretical values of cross correlation

kk=D/dBS;
z1=kc*kk*sind(epsilon).*Delta*[0:Nsims-1];
z2=0.5*kk^2*kc*sqrt(1-(3*cosd(epsilon)^2/4)).*Delta*[0:Nsims-1];
rhotheoretical=besselj(0,z1).^2.*besselj(0,z2).^2;

figure,plot(Delta*[0:Nsims-1]/lambdac,rhotheoretical,'k',Delta*[0:Nsims-1]/lambdac,rho,'k.')
xlabel('Antenna separation (wavelengths)')
ylabel('Cross-correlation coefficient at BS')
legend('Theoretical', 'Simulated', 'Location', 'SouthWest');
