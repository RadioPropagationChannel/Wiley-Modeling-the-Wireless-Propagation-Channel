function k=cprop(frec,er,sigma)
% DESCRIPTION:
%  The function calculates the propagation constant
% INPUTS
%  frec: frequency
%  er: relative permeability 
%  sigma: conductivity
% OUTPUT
%  k: propagation constant

mu0=4*pi*1e-7;
e0=1e-9/(36*pi);

w=2*pi*frec;
e=er*e0;

alfa=w*sqrt(mu0*e/2)*sqrt(sqrt(1+(sigma/(w*e))^2)-1);
beta=w*sqrt(mu0*e/2)*sqrt(sqrt(1+(sigma/(w*e))^2)+1);

k=alfa+j*beta;