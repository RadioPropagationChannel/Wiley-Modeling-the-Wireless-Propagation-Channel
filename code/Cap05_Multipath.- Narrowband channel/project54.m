%
%project54
%
% =============== RESET ==========================================
clear 
close all
clc
% =============== INPUTS ===================================
a=1;
sigma=0.12;    %0.2236   %0.3;
fs=120;         % Hz sampling frequency 
ts=1/fs;       % s sample spacing 
% =============== Butterworth filter parameters =========================
Wp=0.1;             % Actual freq is fs*Wp/2        
Ws=0.3;
Rp=3; %dB
Rs=40; %dB
% =============== RESULTS ================================
leng=10000; % Simulated series length 
lev=-15;   % RMS level
nfftpts=1024; %255;  % No. of FFT points

% =============== RAYLEIGH SERIES ===============================

ry=rayleigh(sigma,leng);    % FUNCION
ry_mod=abs(ry);
t_axis=[0:leng-1]*ts;
figure,plot(t_axis,20*log10(ry_mod),'k')
axis([0 max(t_axis) -50 10])
xlabel('Time (s)')
ylabel('Relative signal level (dB)')
title('Rayleigh fading')
grid
% =============== RICE SERIES ===============

rc=rice(a,sigma,leng);
rc_mod=abs(rc);
figure, plot(t_axis,20*log10(rc_mod),'k')
axis([0 max(t_axis) -50 10])
xlabel('Time (s)')
ylabel('Relative signal level (dB)')
title('Rice fading')
grid

% =============== CDF Rayleigh ===============

[xry,yry]=fCDF(ry_mod);            
[yryt]=RayleighCDF(sigma,xry);
figure, plot(xry,yry,'k:',xry,yryt,'k')
axis([0 2 0 1])
grid
xlabel('Magnitude of complex envelope')
ylabel('Probability the abscissa is not exceeded')
legend('Simulated','Rayleigh theoretical')
title('Simulated and Rayleigh theoretical CDFs')

% =============== CDF Rice ===============

[xrc,yrc]=fCDF(rc_mod);
yrct=[];
for i=0:0.01:ceil(max(rc_mod)),
    yrct=[yrct riceCDF(a,sigma,i)];
end
xrct=[0:0.01:ceil(max(rc_mod))];
figure, plot(xrc,yrc,'k:',xrct,yrct,'k')
axis([0 2 0 1])
grid
xlabel('Magnitude of complex envelope')
ylabel('Probability the abscissa is not exceeded')
legend('Simulated','Rice theoretical')
title('Simulated and Rice theoretical CDFs')
% ===============  Filtered Rayleigh series ========================

[ryfaux, B, A]=filtersignal(ry,Wp,Ws,Rp,Rs);
[H,fre]=freqz(B,A,512,fs);      % For computinf the filter's gain 
% Calculation of filter gain
     [h,T]=impz(B,A);
     h2=h.^2;
     gainFaux=sqrt(cumsum(h2));
     gainF=gainFaux(length(h))
     ryf=ryfaux/gainF;
% ........go on .........
ryf_mod=abs(ryf);                        
figure,plot(t_axis,20*log10(ryf_mod),'k')
axis([10 max(t_axis) -50 10])
xlabel('Time (s)')
ylabel('Reletive signal level (dB)')
title('Filtered Rayleigh fading')
grid
% =============== Filter resposne ========================================

[H,fre]=freqz(B,A,512,fs);
figure, plot(fre,abs(H),'k')
axis([0 fs/2 0 1])
xlabel('Frequence (Hz)')
ylabel('Magnitude of frequency response (lin. units)')
title('Doppler filter')
grid

figure, plot(fre,20*log10(abs(H)),'k')
xlabel('Frequence (Hz)')
ylabel('Magnitude of frequency response (dB)')
title('Doppler filter')
auxx=axis;
axis([0 fs/2 auxx(3) 0])
grid
% =============== Filtered Rice series ========================

rcf=(real(ryf)+a)+j*imag(ryf);
rcf_mod=abs(rcf);
figure,plot(t_axis,20*log10(rcf_mod),'k')
axis([10 max(t_axis) -50 10])
xlabel('Time (s)')
ylabel('Reletive signal level (dB)')
title('Filtered Rice fading')
grid
% =============== Filtered and unfiltered CDFs. Rayleigh case ===========

[xry,yry]=fCDF(ry_mod);            
[xryf,yryf]=fCDF(abs(ryf)); 
figure,plot(xry,yry,'k:',xryf,yryf,'k')
axis([0 2 0 1])
grid
xlabel('Magnitude of complex envelope')
ylabel('Probability the abscissa is not exceeded')
legend('Before filtering','After filtering')
title('Rayleigh series before and after filtering')

% ================ Filtered and unfiltered CDFs. Rice case ==============

[xrc,yrc]=fCDF(rc_mod);
[xrcf,yrcf]=fCDF(abs(rcf));
figure, plot(xrc,yrc,'k:',xrcf,yrcf,'k')
axis([0 2 0 1])
grid
xlabel('Magnitude of complex envelope')
ylabel('Probability the abscissa is not exceeded')
legend('Before filtering','After filtering')
title('Rice series before and after filtering')

% =================== Unfiltererd Rayleigh series spectrum===================

[sry,Fry]=spectrumsignal(ry,nfftpts,fs);
figure,plot(Fry,10*log10((abs(sry)).^2)-max(10*log10((abs(sry)).^2)),'k')                              
xlabel('Frequency (Hz)')
ylabel('Frequency response (dB/Max)')
title('Rayleigh Doppler spectrum before filtering')
grid
% =================== Unfiltererd Rice series spectrum ===================

[src,Frc]=spectrumsignal(rc,nfftpts,fs);
figure,plot(Frc,10*log10((abs(src)).^2)-max(10*log10((abs(src)).^2)),'k')                             
xlabel('Frequency (Hz)')
ylabel('Frequency response (dB/Max)')
title('Rice Doppler spectrum before filtering')
grid
%=============== Spectrum of filtered Rayleigh series ===========

[sryf,Fryf]=spectrumsignal(ryf,nfftpts,fs);
figure, plot(Fryf,10*log10((abs(sryf)).^2)-max(10*log10((abs(sryf)).^2)),'k')                              
xlabel('Frequency (Hz)')
ylabel('Frequency response (dB/Max)')
title('Filtered Rayleigh Doppler spectrum')
grid
%================== Spectrum of filtered Rice series ===================

% rcf=rcf-mean(rcf);                %<-------------- dudas
[srcf,Frcf]=spectrumsignal(rcf,nfftpts,fs);
figure, plot(Frcf,10*log10((abs(srcf)).^2)-max(10*log10((abs(srcf)).^2)),'k')                             
xlabel('Frequency (Hz)')
ylabel('Frequency response (dB/Max)')
title('Filtered Rice Doppler spectrum')
grid
% ===================== Rayleigh afd ====================================

[ejexafd,afd]=afduration(ryf_mod,ts);
figure, plot(ejexafd,afd,'k')
xlabel('Level (dB/LOS)')
ylabel('afd (s)')
title('Average fade durations, afd ')
grid

% % ====================== Filter spectrum =================================
% 
% [sryf,Fryf]=spectrumsignal(ryf,nfftpts,fs);
% figure,plot(Fryf,10*log10((abs(sryf)).^2)-max(10*log10((abs(sryf)).^2)),'k')
% axis([0 max(Fryf) -200 0])
% xlabel('Frequency (Hz)')
% ylabel('Frequency response (dB/Max)')
% title('Filtered Rayleigh Doppler spectrum')
% grid
% ======================= Filter impulse resposne =======================

figure,plot(T*ts,h,'k')
xlabel('Time (s)')
title('Doppler filter. Impulse response')
grid