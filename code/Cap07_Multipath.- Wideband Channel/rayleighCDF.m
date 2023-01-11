function [x,y]=rayleighCDF(sigma);
x=[0:0.01:10];
y=1-exp(-(x.^2)./(2*sigma^2));