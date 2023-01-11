function [s,d]=s_d(tita_i,er,len)
% DESCRIPTION
%  This function calculates parameters for the transmission coefficients
% INPUTS
%  tita_i: incident angle
%  er: relative permeability
%  len: thickness
% OUTPUTS
%  d: difference of paths between reflected or transmitted rays 
%  s: distance of the ray inside the wall


% la er de la formula de 'd' tiene que estar elevada al cuadrado.
% la variable len es el grosor.

tita_i=tita_i*pi/180;
s=len/sqrt(1-(sin(tita_i)^2)/er);
if ((tita_i < 1e-3)&(tita_i > -1e-3))
    d=0;
else
    d=2*len/sqrt((er/sin(tita_i))^2-1);
end