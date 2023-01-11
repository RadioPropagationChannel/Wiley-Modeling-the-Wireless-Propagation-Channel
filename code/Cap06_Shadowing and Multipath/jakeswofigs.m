function [filtertimedomain]=jakeswofigs(freqstep,fm,fs)

% JAKES filter without figures

ts=1/fs;
A=[];
A2=[];

%freq=[0:freqstep:round(fs)]-fs/2;

freq=[freqstep/2:freqstep:round(fs/2)]';
freq2=-freq;
freq2=sort(freq2);

for ii=1:length(freq)
    if abs(freq(ii))< fm-2,
        auxA=1/sqrt(1-(freq(ii)/fm).^2);
        A=[A; auxA];
        A2=[auxA; A2];
    else
        A=[A; 0]; 
        A2=[0; A2];
    end
end

freq=[freq2; freq];
A = [A2; A];

% figure,stem(freq,A)
% title('Jackes filter')
% xlabel('Frequency (Hz)')
% ylabel('Frequency response, H(f) (linear units)')
%     
% figure,plot(freq,A)
% title('Jakes filter. Shifted');

a=ifft(A);
filter_wrap=fftshift(fft(a));        % JACKES Filter
filter_wrap=real(filter_wrap);       % Imaginary part is extremely small

shiftedfreq=freq+fs/2;
% figure, plot(shiftedfreq,filter_wrap)
% title('Jakes filter');

% Convert filter to time domain ========================================

filtertimedomain=fftshift(ifft(filter_wrap));        
filtertimedomain=real(filtertimedomain);
filtEnergy=sum(filtertimedomain.^2);
filtertimedomain=filtertimedomain/sqrt(filtEnergy);

timeaxisfilter=[0:length(filtertimedomain)-1]*ts;

% figure,stem(timeaxisfilter, filtertimedomain)
% title('Jackes filter')
% xlabel('Time (s)')
% ylabel('Impulse response, h(t) (linear units)')
