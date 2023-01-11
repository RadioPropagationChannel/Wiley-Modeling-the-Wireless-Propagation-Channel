%
% Intro 81
%
%=======================================================================
clear
close all
clc
%=======================================================================
f=2000             % Frequency in MHz

% Very dry Ground ======================================================
er1=3;             % Relative permitivity.
sigma1=7e-4;       % Conductivity (S/m).

% Wet Ground ===========================================================
er2=15;             % Relative permitivity.
sigma2=4e-1;        % Conductivity (S/m).
%=======================================================================
% Reflection and Transmission coefficients Very dry Ground==============
Mrpar=[];
Mrper=[];
Mtpar=[];
Mtper=[];
mang_inc=0:0.01:90;
epp=ep(f*1e6,er1,sigma1);
for index=1:length(mang_inc)
    
    [rpar,rper]=rfresnel(mang_inc(index)*pi/180,epp);
    [tpar,tper]=tfresnel(mang_inc(index)*pi/180,epp);
    Mrpar=[Mrpar;rpar];
    Mrper=[Mrper;rper];
    Mtpar=[Mtpar;tpar];
    Mtper=[Mtper;tper];
    
end
figure
plot(mang_inc,abs(Mrpar),'-k','LineWidth',2);
hold on
plot(mang_inc,abs(Mrper),'--k','LineWidth',2);
hold on
plot(mang_inc,abs(Mtpar),':k','LineWidth',2);
hold on
plot(mang_inc,abs(Mtper),'-.k','LineWidth',2);
xlabel('Incident angle \theta (degrees) ')
ylabel('\midCoefficient\mid')
grid
title('Reflection and Transmission coefficients-Very dry ground')
legend('Reflection C. Par','Reflection C. Per','Transmission C. Par'...
    ,'Transmission C. Per')

% Reflection and Transmission coefficients Wet Ground==================
Mrpar=[];
Mrper=[];
Mtpar=[];
Mtper=[];
epp=ep(f*1e6,er2,sigma2);
for index=1:length(mang_inc)
    
    [rpar,rper]=rfresnel(mang_inc(index)*pi/180,epp);
    [tpar,tper]=tfresnel(mang_inc(index)*pi/180,epp);
    Mrpar=[Mrpar;rpar];
    Mrper=[Mrper;rper];
    Mtpar=[Mtpar;tpar];
    Mtper=[Mtper;tper];
    
end

figure
plot(mang_inc,abs(Mrpar),'-k','LineWidth',2);
hold on
plot(mang_inc,abs(Mrper),'--k','LineWidth',2);
hold on
plot(mang_inc,abs(Mtpar),':k','LineWidth',2);
hold on
plot(mang_inc,abs(Mtper),'-.k','LineWidth',2);
xlabel('Incident angle \theta (degrees) ')
ylabel('\midCoefficient\mid')
grid
title('Reflection and Transmission coefficients-Wet ground')
legend('Reflection C. Par','Reflection C. Per','Transmission C. Par',...
    'Transmission C. Per')

clear Mrpar Mrper Mtpar Mtper
