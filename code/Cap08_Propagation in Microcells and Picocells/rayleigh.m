function [r]=rayleigh(sigma,NSamples)
% generates a Rayliegh time series of lenght NSamples and modal value sigma 

ii=randn(NSamples,1).*sigma;
qq=randn(NSamples,1).*sigma;
r=ii+j*qq;
