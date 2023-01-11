function [tpar_g,tper_g]=tfresnel_g(rpar,rper,k0,k1,d,s,tita_i)
% DESCRIPTION:
%  This function calculates the generalized transmitted coefficients (parallel and perpendicular)
% INPUTS
%  rpar: parallel reflection coefficient
%  rper: perpendicular reflection coefficient
%  k0: first medium propagation constant
%  k1: second medium propagation constant
%  d: difference of paths between reflected or transmitted rays 
%  s: distance of the ray inside the wall
%  tita_i: incident angle
% OUTPUT
%  tpar_g: general parallel transmission coefficient
%  tper_g: general perpendicular transmission coefficient
j=sqrt(-1);
alfa=real(k1);
beta=imag(k1);

c=exp(-1*(alfa*s+j*beta*s));
num=(1-rpar^2)*c;
den=1-rpar^2*(c^2*exp(j*k0*d*sin(tita_i)));
tpar_g=(num/den);


num2=(1-rper^2)*c;
den2=1-rper^2*(c^2*exp(j*k0*d*sin(tita_i)));
tper_g=(num2/den2);

