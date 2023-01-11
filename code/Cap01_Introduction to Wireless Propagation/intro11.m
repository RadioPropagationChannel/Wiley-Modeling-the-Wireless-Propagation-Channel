%
% intro11
%
%=======================================================================
clear
close all
clc
%=======================================================================
xaxis=[0:0.1:4];
sigma=1;
%=======================================================================
pdf=Rayleighpdf(sigma,xaxis);
CDF=RayleighCDF(sigma,xaxis);

modalvalue=1;
standarddeviation=0.665;
meanvalue=1.25
RMSvalue=1.41
medianvalue=1.18;

figure,plot(xaxis,pdf,'k',xaxis,CDF,'k','LineWidth',2)
hold on
plot([modalvalue modalvalue],[0 1],'k:')
plot([standarddeviation standarddeviation],[0 1],'k:')
plot([RMSvalue RMSvalue],[0 1],'k:')
plot([medianvalue medianvalue],[0 1],'k:')
plot([meanvalue meanvalue],[0 1],'k:')

hold off
xlabel('Random variable, r')
ylabel('pdf and CDF')
title('Rayleigh distribution')



