%
% intro32   (Hata model)
%
%========================================================================
close all
clear
clc
%========================================================================
fMHz=2000;        % freq in MHz
d=[0.1:0.1:10];   % distance in km
hm=1.5;           % mobile antenna height in m
ht=30;            % base station antenna height in m
Cm_p=0;   
Cm_g=3;   

am_p=(1.1*log10(fMHz)-0.7)*hm-(1.56*log10(fMHz)-0.8);

if (fMHz<=200)
    am_g=8.29*((log10(1.54*hm))^2)-1.1;    % am=8.29*((log10(1.54)*hm)^2)-1.1 ??????
elseif (fMHz>400)
    am_g=3.2*((log10(11.75*hm))^2)-4.97;   % am=3.2*((log10(11.75)*hm)^2)-4.97 ?????
end


if (fMHz<1500 && fMHz>150)
    % Hata basic loss (Urban environment).
    AA=69.55+26.16*log10(fMHz)-13.82*log10(ht);
    AA_g=AA-am_g;
    AA_p=AA-am_p;
elseif (fMHz>=1500 && fMHz<=2000)
    % Hata-COST231 basic loss.
    AA=46.3+33.9*log10(fMHz)-13.82*log10(ht);
    AA_g=AA-am_g+Cm_g;
    AA_p=AA-am_p+Cm_p;
else
    fprintf('Wrong frequency. Correct frequency range (150,2000] MHz');
    return
end

% Losses due to the distance.
BB=44.9-6.55*log10(ht);
nn=BB/10;
L_g=AA_g+BB*log10(d);
L_p=AA_p+BB*log10(d);

% Suburban zone loss.
Ls=L_p-2*((log10(fMHz/28))^2)-5.4;

% Rural zone loss.
Lr=L_p-4.78*((log10(fMHz))^2)+18.33*log10(fMHz)-40.94;

figure,plot(d,-L_g,d,-L_p,'r',d,-Ls,'g',d,-Lr,'c');
legend('Urban. Big city','Urban. Small-medium city','Suburban','Rural');
title('Hata model');
xlabel('d (km)');
ylabel('Attenuation (dB)');