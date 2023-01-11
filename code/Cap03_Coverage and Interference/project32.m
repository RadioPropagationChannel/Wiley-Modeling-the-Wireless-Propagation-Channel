%
% project32    (two cross correlated series)
%
%========================================================================
clear
close all
clc
%========================================================================
fMHz=2000;     % Frequency in MHz
Lcorr=30;       % Correlation distance in meters 
MM1=-70;          % Larger area mean in dBm
SS1=9;          % Larger area std or location variability in dB
MM2=-70;          % Larger area mean in dBm
SS2=9;          % Larger area std or location variability in dB

ThresHold=-80;

corr=[1 0.9
      0.9  1];
%=======================================================================
lambdac=300/fMHz;   % wavelength in m
ds=lambdac;
Nsamples=4000;
InterpRate=round(Lcorr/ds)
Lcorr=InterpRate*ds       % slightly correct Lcorr to make it a multiple of ds

% Create Gaussian series ================================================
R1=randn(1,Nsamples);
R2=randn(1,Nsamples);

disp('Cross-correlation coefficient')
corrcoef(R1,R2)

d_axis1=[1:Nsamples]*Lcorr-Lcorr;

d_axis2=([1:1/InterpRate:Nsamples]-1)*Lcorr;

Rinterpolated1=interp1(d_axis1,R1,d_axis2,'spline');
Rinterpolated2=interp1(d_axis1,R2,d_axis2,'spline');

Rinterpolated1=Rinterpolated1*SS1+MM1;
Rinterpolated2=Rinterpolated2*SS2+MM2;

plot(d_axis2,Rinterpolated1,'k:',d_axis2,Rinterpolated2,'k','LineWidth',2)
xlabel('Traversed distance (m)')
ylabel('Slow signal variations (dBm)')
ylabel('Slow variations autocorrelation')

disp('Cross-correlation coefficient')
corrcoef((Rinterpolated1-mean(Rinterpolated1))/std(Rinterpolated1),...
    (Rinterpolated2-mean(Rinterpolated2))/std(Rinterpolated2))

plot([d_axis2(1); d_axis2(end)],[ThresHold; ThresHold],'k:',d_axis2,Rinterpolated1,'k:',...
    d_axis2,Rinterpolated2,'k--',d_axis2,max(Rinterpolated1,Rinterpolated2),'k','LineWidth',2)
xlabel('Traversed distance (m)')
ylabel('Slow signal variations (dBm)')
ylabel('Slow variations autocorrelation')


% evaluate coverage probability =========================================

Nbelow1=find(Rinterpolated1<ThresHold);
disp('Coverage probability Signal 1')
Nabove1=find(Rinterpolated1>=ThresHold);
CovProb1=length(Nabove1)/length(Rinterpolated1)

Nbelow2=find(Rinterpolated2<ThresHold);
disp('Coverage probability Signal 2')
Nabove2=find(Rinterpolated2>=ThresHold);
CovProb2=length(Nabove2)/length(Rinterpolated2)

Nbelow12=find(Rinterpolated1<ThresHold | Rinterpolated2<ThresHold );
disp('Coverage probability with diversity Signals 1&2')
Nabove12=find(Rinterpolated1>=ThresHold | Rinterpolated2>=ThresHold);
CovProb12=length(Nabove12)/length(Rinterpolated1)

% Theoretical coverage probability 

kk1=(ThresHold-MM1)/SS1;
a1=0.5*erfc(kk1/sqrt(2))
u1=1-a1;

kk2=(ThresHold-MM2)/SS2;
a2=0.5*erfc(kk2/sqrt(2))
u2=1-a2;

aTuncor=1-u1*u2

% ======================================================================
% we reuse R1 and R2 from above

%corrcoef(R1,R2)

cholcorr=chol(corr)';
Rcorr=cholcorr*[R1;R2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%corrcoef(Rcorr(1,:),Rcorr(2,:))

Rinterpolated1=interp1(d_axis1,Rcorr(1,:),d_axis2,'spline');
Rinterpolated2=interp1(d_axis1,Rcorr(2,:),d_axis2,'spline');


Rinterpolated1=Rinterpolated1*SS1+MM1;
Rinterpolated2=Rinterpolated2*SS2+MM2;

figure,plot(d_axis2,Rinterpolated1,'k:',d_axis2,Rinterpolated2,'k','LineWidth',2)
xlabel('Traversed distance (m)')
ylabel('Slow signal variations (dBm)')
ylabel('Slow variations autocorrelation')

disp('Cross-correlation coefficient')
corrcoef((Rinterpolated1-mean(Rinterpolated1))/std(Rinterpolated1),...
    (Rinterpolated2-mean(Rinterpolated2))/std(Rinterpolated2))

plot([d_axis2(1); d_axis2(end)],[ThresHold; ThresHold],'k:',d_axis2,Rinterpolated1,'k:',...
    d_axis2,Rinterpolated2,'k--',d_axis2,max(Rinterpolated1,Rinterpolated2),'k','LineWidth',2)
xlabel('Traversed distance (m)')
ylabel('Slow signal variations (dBm)')
ylabel('Slow variations autocorrelation')


% evaluate coverage probability ========================================

Nbelow1=find(Rinterpolated1<ThresHold);
disp('Coverage probability Signal 1')
Nabove1=find(Rinterpolated1>=ThresHold);
CovProb1=length(Nabove1)/length(Rinterpolated1)

Nbelow2=find(Rinterpolated2<ThresHold);
disp('Coverage probability Signal 2')
Nabove2=find(Rinterpolated2>=ThresHold);
CovProb2=length(Nabove2)/length(Rinterpolated2)

Nbelow12=find(Rinterpolated1<ThresHold | Rinterpolated2<ThresHold );
disp('Coverage probability with diversity Signals 1&2')
Nabove12=find(Rinterpolated1>=ThresHold | Rinterpolated2>=ThresHold);
CovProb12=length(Nabove12)/length(Rinterpolated1)

%========================================================================
% evaluate joint coverage probability theoretically
% in case of partial correlation

rho=corr(1,2);    % rho12 correlation coefficient
aTcorr=1-(rho*sqrt(u1*(1-u1))*sqrt(u2*(1-u2))+u1*u2)


