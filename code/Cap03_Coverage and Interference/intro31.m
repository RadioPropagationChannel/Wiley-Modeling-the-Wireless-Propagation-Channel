%
% intro31
%
%=======================================================================
clear
close all
clc
%======================================================================


sigmaL_over_n=[0.1:0.1:8];

b=10*log10(exp(1))./(sigmaL_over_n*sqrt(2));

% Pobability range for the plots
Px0=[.5 .55 .6 .65 .7 .75 .8 .85 .9 .95];
a=erfinv(1-2*Px0);

figure,hold on
for m=1:length(a)
    Fu=0.5*(1-erf(a(m))+exp((1-(2*a(m)*b))./(b.^2)).*(1-erf((1-(a(m)*b))./b)));
    plot(sigmaL_over_n,Fu)
end
hold off
axis([sigmaL_over_n(1) sigmaL_over_n(end) 0.5 1])
xlabel('Location variability-propagation exponent ratio, \sigma_L/n');
ylabel('Area coverage probability, F_u ');