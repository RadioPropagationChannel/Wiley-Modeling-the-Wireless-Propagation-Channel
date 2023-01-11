%
% Intro 83      % slab
%
%=======================================================================
clear
close all
clc
%=======================================================================
ang_inc=10;       % Incident angle (degrees)
frec=2.0e9;       % Frequency (Hz)
%Medium 1===============================================================
er0=1;            % Relative permitivity medium 1.
sigma0=0;         % Conductivity medium 1(S/m).  
%Medium 2===============================================================
er1=2.5;          % Relative permitivity medium 2.
sigma1=0.1;       % Conductivity medium 2(S/m).   
width=1e-1;       % Thickness of slab (m).
%======================================================================
%Generalised reflection and penetration coefficient as a function of
%frequency
mRpar_g=[];
mRper_g=[];
mTpar_g=[];
mTper_g=[];
mfrec=[1:0.01:9]*1e9;

for indice=1:length(mfrec)
    
    k0=cprop(mfrec(indice),er0,sigma0);
    k1=cprop(mfrec(indice),er1,sigma1);
    epp=ep(mfrec(indice),er1,sigma1);
    [rpar,rper]=rfresnel(ang_inc*pi/180,epp);
    [s,d]=s_d(ang_inc,er1,width);
    [Rpar_g,Rper_g]=rfresnel_g(rpar,rper,k0,k1,d,s,ang_inc*pi/180);
    [tpar_g,tper_g]=tfresnel_g(rpar,rper,k0,k1,d,s,ang_inc*pi/180);
    mRpar_g=[mRpar_g;Rpar_g];
    mRper_g=[mRper_g;Rper_g];
    mTpar_g=[mTpar_g;tpar_g];
    mTper_g=[mTper_g;tper_g];
  
end
%Plot Generalised reflection and penetration coefficient 
figure
plot(mfrec,abs(mRpar_g),'-k')
hold on
plot(mfrec,abs(mRper_g),'--k')
hold on
plot(mfrec,abs(mTpar_g),':k')
hold on
plot(mfrec,abs(mTper_g),'-.k')
xlabel('Frequency (Hz)')
title('Generalized reflection coefficient as a function of frequency')
ylabel('\midCoefficient\mid')
legend('Reflection C. Par','Reflection C. Per','Transmission C. Par'...
    ,'Transmission C. Per')
%Generalised reflection and penetration coefficient as a function of
%thickness
mwidth=0:0.01:1;
mRpar_g=[];
mRper_g=[];
mTpar_g=[];
mTper_g=[];
for indice=1:length(mwidth)
    k0=cprop(frec,er0,sigma0);
    k1=cprop(frec,er1,sigma1);
    epp=ep(frec,er1,sigma1);
    [rpar,rper]=rfresnel(ang_inc*pi/180,epp);
    [s,d]=s_d(ang_inc,er1,mwidth(indice));
    [Rpar_g,Rper_g]=rfresnel_g(rpar,rper,k0,k1,d,s,ang_inc*pi/180);
    [tpar_g,tper_g]=tfresnel_g(rpar,rper,k0,k1,d,s,ang_inc*pi/180);
    mRpar_g=[mRpar_g;Rpar_g];
    mRper_g=[mRper_g;Rper_g];
    mTpar_g=[mTpar_g;tpar_g];
    mTper_g=[mTper_g;tper_g];
  
end
%Plot Generalised reflection and penetration coefficient 
figure
plot(mwidth,abs(mRpar_g),'-k')
hold on
plot(mwidth,abs(mRper_g),'--k')
hold on
plot(mwidth,abs(mTpar_g),':k')
hold on
plot(mwidth,abs(mTper_g),'-.k')
xlabel('Thickness (m)')
title('Generalized reflection coefficient as a function of thickness')
ylabel('\midCoefficient\mid')
legend('Reflection C. Par','Reflection C. Per','Transmission C. Par'...
    ,'Transmission C. Per')


%Plot Generalised reflection and penetration coefficient in dB
figure
plot(mwidth,20*log10(abs(mRpar_g)),'-k')
hold on
plot(mwidth,20*log10(abs(mRper_g)),'--k')
hold on
plot(mwidth,20*log10(abs(mTpar_g)),':k')
hold on
plot(mwidth,20*log10(abs(mTper_g)),'-.k')
xlabel('Thickness (m)')
title('Generalized reflection coefficient as a function of thickness')
ylabel('\midCoefficient\mid')
legend('Reflection C. Par','Reflection C. Per','Transmission C. Par'...
    ,'Transmission C. Per')