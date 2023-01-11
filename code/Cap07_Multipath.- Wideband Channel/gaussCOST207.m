% ========================================================================
% Gaussian filter 
% ========================================================================

function [filtertimedomain]=gaussCOST207(freqstep,fm,fs,fr11,fr12,fr21,fr22,A)

ts=1/fs;

freq1=[freqstep/2:freqstep:round(fs/2)]';  % create frequency axis
freq2=-freq1;
freq2=sort(freq2);
freq=[freq2; freq1];


% Frequency domain
H=exp(-(freq-fr11).^2./(2*fr12.^2))+(10^(0.2*A))*exp(-(freq-fr21).^2./(2*fr22.^2));

% figure,stem(freq,H)
figure,plot(freq,H,'k')
title('Gaussian filter')
xlabel('Frequency (Hz)')
ylabel('Frequency response, H(f) (linear units)')

    
% Convert filter to time domain ========================================
 
timeaxisfilter=([0:length(freq)-1]-round(length(freq)/2))*ts;
filtertimedomain=sqrt(2*pi)*fr12*exp(-2*(pi*fr12*timeaxisfilter).^2).*exp(j*2*pi*fr11.*timeaxisfilter)+...
    (10^(0.2*A))*sqrt(2*pi)*fr22*exp(-2*(pi*fr22*timeaxisfilter).^2).*exp(j*2*pi*fr21.*timeaxisfilter);


filtEnergy=sum(abs(filtertimedomain).^2);
filtertimedomain=filtertimedomain/sqrt(filtEnergy);    % the energy of the filter is one

figure,stem(timeaxisfilter, abs(filtertimedomain),'k')
title('Gaussian filter')
xlabel('Time (s)')
ylabel('Impulse response, h(t), magnitude (linear units)')

figure,plot(timeaxisfilter, unwrap(angle(filtertimedomain)),'k')
title('Gaussian filter')
xlabel('Time (s)')
ylabel('Impulse response, h(t), phase (linear units)')

