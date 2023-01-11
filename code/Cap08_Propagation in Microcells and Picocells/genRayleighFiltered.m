function [t_axis,ryf]=genRayleighFiltered(leng, fs,sigma,Wp,Ws,Rp,Rs);
%
%
% ======================================================================
 
ts=1/fs;       % sample spacing (s)

% =============== RAYLEIGH SERIES ===============================

ry=rayleigh(sigma,leng);    % FUNCION
ry_mod=abs(ry);
t_axis=[0:leng-1]*ts;

% ===============  Filtered Rayleigh series ========================

[ryfaux, B, A]=filtersignal(ry,Wp,Ws,Rp,Rs);
[H,fre]=freqz(B,A,512,fs);      % For computing the filter's gain 
% Calculation of filter gain
     [h,T]=impz(B,A);
     h2=h.^2;
     gainFaux=sqrt(cumsum(h2));
     gainF=gainFaux(length(h));
%--------------------------------------------------------- 
ryf=ryfaux/gainF;



