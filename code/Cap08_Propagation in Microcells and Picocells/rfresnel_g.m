function [rpar_g,rper_g]=rfresnel_g(rpar,rper,k0,k1,d,s,tita_i)
% DESCRIPTION
%  This function calculates the generalized reflection coefficients (parallel and perpendicular)
% INPUTS
%  rpar: parallel reflection coefficient
%  rper: perpendicular reflection coefficient
%  k0: first medium propagation constant
%  k1: second medium propagation constant
%  d: difference of paths between reflected or transmitted rays 
%  s: distance of the ray inside the wall
%  tita_i: incident angle
% OUTPUTS
%  rpar_g: general parallel reflection coefficient
%  rper_g: general perpendicular reflection coefficient
j=sqrt(-1);
alfa=real(k1);
beta=imag(k1);


c=exp(-2*alfa*s+j*(k0*d*sin(tita_i)-2*beta*s));
num=(1-rpar^2)*c;
den=1-rpar^2*c;
rpar_g=rpar*(1-num/den);

num2=(1-rper^2)*c;
den2=1-rper^2*c;
rper_g=rper*(1-num2/den2);