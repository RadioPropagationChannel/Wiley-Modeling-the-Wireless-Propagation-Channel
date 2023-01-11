function [tpar,tper]=tfresnel(ang_inc,ep)
% DESCRIPTION
%  This function calculates the reflection coefficients (parallel and perpendicular)
% INPUTS
%  ang_inc: incident angle
%  ep: parameter for the Fresnel function
% OUTPUT
%  tpar: parallel penetration coefficient
%  tper: perpendicular penetration coefficient


tper=2.*cos(ang_inc)./(cos(ang_inc)+sqrt(ep-sin(ang_inc).^2));
tpar=2.*sqrt(ep).*cos(ang_inc)./(ep.*cos(ang_inc)+sqrt(ep-sin(ang_inc).^2));
