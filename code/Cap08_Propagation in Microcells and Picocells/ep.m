function ep=ep(f,er,sigma)
% DESCRIPTION:
% This function calculates one parameter for the Fresnel function
% INPUTS
%  f: frequency
%  er: relative permeability 
%  sigma: conductivity 
% OUTPUT
%  ep: parameter for the Fresnel function
 
 

e0=1e-9/(36*pi);

w=2*pi*f;
e=er*e0;

ep=e/e0-j*sigma/(w*e0);