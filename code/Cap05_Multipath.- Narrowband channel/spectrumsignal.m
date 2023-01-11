function [s,F]=spectrumsignal(r,nfftpts,fs)
% S=20*log10(abs(fftshift(fft(r,nfftpts))));
% S=S-max(S);
s=fftshift(fft(r,nfftpts));
pasofrec=fs/nfftpts;
F=([0:nfftpts-1]-(nfftpts/2))*pasofrec;
