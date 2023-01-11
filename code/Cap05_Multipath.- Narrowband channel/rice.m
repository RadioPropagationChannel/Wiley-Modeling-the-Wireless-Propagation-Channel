function [r]=rice(a,sigma,NSamples);
% Generates a Ricean time series of length NSamples 
% paramters a and sigma

ii=randn(NSamples,1)*sigma+a;
qq=randn(NSamples,1)*sigma;
r=ii+j*qq;