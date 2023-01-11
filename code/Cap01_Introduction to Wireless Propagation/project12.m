%
% project12
%
% INITIALIZE =========================================================
clear 
close all
clc
% ====================================================================

lambdac=300/2000;             % wavelength for a 2 GHz carrier

NofWavelengths=10;   
SamplesperWavelength=4
lengthwindow=NofWavelengths*SamplesperWavelength
windowinmeters=NofWavelengths*lambdac

load series12                 % Suzuki series

% distanceaxis=series12(1:10000,1);   % distance axis in m (sampled wavelength/4)
% P=series12(1:10000,2);              % P in dBm

distanceaxis=series12(:,1);   % distance axis in m (sampled wavelength/4)
P=series12(:,2);              % P in dBm


figure,plot(distanceaxis,P,'k')
title('Received power')
ylabel('Received power (dBm)')
xlabel('Traveled distance (m)')

p=10.^(P/10);   % now p is in mW
p=p/1000;       % now p is in W

figure,plot(distanceaxis,p,'k')
title('Received power')
ylabel('Received power (W)')
xlabel('Traveled distance (m)')

v=sqrt(2*50*p);                % compute the voltage

figure,plot(distanceaxis,v,'k')
title('Received voltage')
ylabel('Received voltage (v)')
xlabel('Traveled distance (m)')

% Low-pass filtering ===================================================
disp('Filtering overall signal')
averagingwindow=ones(1,lengthwindow)/lengthwindow;
vfiltered=conv(averagingwindow,v);
vfiltered=vfiltered(1:length(v));   % discard last samples after convolution 
disp('Endo of filtering')

% ======================================================================
vfast=v./vfiltered;    % remove slow variations and nomalize

% PLOTTING ============================================================
figure,plot(distanceaxis,vfiltered,'k')
title('Received voltage. Slow variations')
ylabel('Received voltage. Slow variations (v)')
xlabel('Traveled distance (m)')

figure,plot(distanceaxis,20*log10(vfiltered*1e6),'k')
title('Received voltage. Slow variations')
ylabel('Received voltage. Slow variations (dB  \muV)')
xlabel('Traveled distance (m)')

figure,plot(distanceaxis,vfast,'k')
axis([0 distanceaxis(end) 0 5])
title('Received normalized voltage. Fast variations')
ylabel('Received voltage. Fast variations (lin. units)')
xlabel('Traveled distance (m)')

figure,plot(distanceaxis,20*log10(vfast),'k')
title('Received normalized voltage. Fast variations')
ylabel('Received voltage. Fast variations (dB)')
xlabel('Traveled distance (m)')

figure,plot(distanceaxis,20*log10(v*1e6),'k:')
hold on
plot(distanceaxis,20*log10(vfiltered*1e6),'k','LineWidth',2)
hold off
title('Received voltage. Overall and Slow variations')
ylabel('Received voltage. Overall and Slow variations (dB  \muV)')
xlabel('Traveled distance (m)')
title('ZOOM IN TO SEE OVERALL AND SLOW VARIATIONS')

% STATISTICAL ANALYSIS ==================================================
disp('Mean of fast variations'), mean(vfast(length(averagingwindow):end))
disp('Standard deviation of slow variations (dB)'),SS=std(20*log10(vfiltered(length(averagingwindow):end)*1e6))
disp('Mean of slow variations (dB)'),MM=mean(20*log10(vfiltered(length(averagingwindow):end)*1e6))

sigma=mean(vfast(length(averagingwindow):end))/1.25;
[CDFxfast,CDFyfast]=fCDF(vfast(length(averagingwindow):end))
[CDFyTheoreticalfast]=RayleighCDF(sigma,CDFxfast)
figure,plot(CDFxfast,CDFyfast,'k:',CDFxfast,CDFyTheoreticalfast,'k')
title('CDFs of measured fast variations and corresponding theoretical distribution')
xlabel('Normalized voltage')
ylabel('Probability the abscissa is not exceeded')
legend('Measured fast variations','Theoretical fast variations')

[CDFxslow,CDFyslow]=fCDF(20*log10(vfiltered(length(averagingwindow):end)*1e6));
[CDFxoverall,CDFyoverall]=fCDF(20*log10(v(length(averagingwindow):end)*1e6));

figure,plot(CDFxslow,CDFyslow,'k',CDFxoverall,CDFyoverall,'k:')
title('Measured CDFs of Overall and Slow variations')
xlabel('Voltage (dB  \muV)')
ylabel('Probability the abscissa is not exceeded')
legend('Slow variations','Overall variations')

[CDFyslowtheoretical]=GaussianCDF(MM,SS,CDFxslow);
figure,plot(CDFxslow,CDFyslow,'k.-',CDFxslow,CDFyslowtheoretical,'k')
title('CDFs of measured slow variations and theoretical distribution')
xlabel('Voltage (dB  \muV)')
ylabel('Probability the abscissa is not exceeded')
legend('Measured slow variations','Theoretical slow variations')
