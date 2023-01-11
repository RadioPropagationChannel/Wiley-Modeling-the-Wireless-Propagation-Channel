%
% project311     (autoregresive filter)
%
%========================================================================
clear
close all
clc
%========================================================================

Lcorr=10;       % Correlation distance in meters 
MM=-80;          % Larger area mean in dB
SS=5;          % Larger area std or location variability in dB

%=======================================================================

ds=1;               % m   
Nsamples=200;       % Number of samples to be generated  

% Fiter parameters =====================================================
cc=exp(-ds/Lcorr);
bb=SS*sqrt(1-cc^2);

A=[1 -cc];
B=1;

% Create Gaussian series ================================================
R=randn(Nsamples,1);
Rfiltered=filter(B,A,R)*bb;
Rfiltered=Rfiltered+MM;


d_axis=[0:Nsamples-1]*ds;
plot(d_axis,R+MM,'k:',d_axis,Rfiltered,'k')
xlabel('Traversed distance (m)')
ylabel('Slow signal variations (dBm)')
legend('Uncorrelated Gaussian','Correlated Gaussian ')


mean_Rfiltered=mean(Rfiltered)

std_Rfiltered=std(Rfiltered)


Rfwithout=(Rfiltered-mean_Rfiltered)/std(Rfiltered);


Rfcorr=xcorr(Rfwithout,'coeff');
Rcorr=xcorr(R,'coeff');
figure,plot([-Nsamples+1:Nsamples-1]*ds,Rcorr,'k:',[-Nsamples+1:Nsamples-1]*ds,Rfcorr,'k'),grid
ylabel('Autocorrelation coefficient')
xlabel('Sample spacing (m)')
legend('Uncorrelated Gaussian','Correlated Gaussian')