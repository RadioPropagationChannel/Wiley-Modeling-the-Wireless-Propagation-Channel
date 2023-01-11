%
% project11
%
% INITIALIZE =========================================================
clear 
close all
clc
% ====================================================================

Nbins=10;

% ==================================================================== 
load series11      

timeaxis=series11(:,1);       % timeaxis in s
P=series11(:,2);              % P in dBm

figure,plot(timeaxis,P,'k')
title('Received power')
ylabel('Received power (dBm)')
xlabel('Elapsed time (s)')

p=10.^(P/10);   % now p is in mW
p=p/1000;       % now p is in W

figure,plot(timeaxis,p,'k')
title('Received power')
ylabel('Received power (W)')
xlabel('Elapsed time (s)')

v=sqrt(2*50*p);

figure,plot(timeaxis,v,'k')
title('Received voltage')
ylabel('Received voltage (v)')
xlabel('Elapsed time (s)')

sigma=mean(v)/1.25

vnorm=v/sigma;    % normalize vortage wrt modal value

% Again compute sigma of normalized v (vnorm)

sigmavnorm=mean(vnorm)/1.25

figure,plot(timeaxis,vnorm,'k')
title('Normalized received voltage')
ylabel('Normalized received voltage (v/mode(v))')
xlabel('Elapsed time (s)')

[CDFx,CDFy]=fCDF(vnorm);

[CDFyTheoretical]=RayleighCDF(sigmavnorm,CDFx);
figure, plot(CDFx,CDFyTheoretical,'k',CDFx,CDFy,'k:') 
title('Sample CDF')
ylabel('Probability the abscissa is not exceeded')
xlabel('Normalized received voltage (v/mode(v))')
legend('Theoretical','Sample')

% Chi-Square test 1st try ===============================================

[HvnormY,HvnormX]=hist(vnorm,Nbins);

HvnormYtheoretical=RayleighHIST(HvnormX,sigmavnorm);
HvnormYtheoretical=HvnormYtheoretical*length(vnorm);
figure,bar(HvnormX,HvnormYtheoretical,0.8,'y')
hold on
bar(HvnormX,HvnormY,0.2,'r')
hold off
title('Histogam of theoretical Rayleigh and of series')
legend('Theoretical','Sample')
xlabel('Normalized signal level')
ylabel('Bin frecuencies')

% Chi-Square parameter ================================================

D2=sum((HvnormYtheoretical-HvnormY).^2./HvnormYtheoretical)


% Chi-Square test 2nd try ===============================================

RR=BINSequalprobRayleigh(Nbins,1)

HvnormY=histc(vnorm,RR)
HvnormY=HvnormY(1:Nbins)


HvnormYtheoretical=1/Nbins;
HvnormYtheoretical=HvnormYtheoretical*length(vnorm);

% showing results =====================================================

theoreticalfreq=ones(Nbins,1).*HvnormYtheoretical;
D2partial=(HvnormYtheoretical-HvnormY).^2./HvnormYtheoretical;

format bank
tableresults=[RR(1:end-1), RR(2:end), HvnormY,theoreticalfreq,D2partial]
format short

% figure,bar(HvnormX,HvnormYtheoretical,0.8,'y')
% hold on
% bar(HvnormX,HvnormY,0.2,'r')
% hold off
% title('Histogam of theoretical Rayleigh and of series')

% Chi-Square parameter ================================================

D2=sum((HvnormYtheoretical-HvnormY).^2./HvnormYtheoretical)


%alpha=1-gammainc(0.5*chi2,0.5*df); significance level


