function [rpar,rper]=rfresnel(ang_inc,ep)
% DESCRIPTION
%  This function calculates the reflection coefficients (parallel and perpendicular)
% INPUTS
%  ang_inc: incident angle
%  ep: parameter for the Fresnel function
% OUTPUT
%  rpar: parallel reflection coefficient
%  rper: perpendicular reflection coefficient
 
rper=(cos(ang_inc)-sqrt(ep-(sin(ang_inc))^2))/(cos(ang_inc)+sqrt(ep-(sin(ang_inc))^2));
rpar=((ep*cos(ang_inc)-sqrt(ep-(sin(ang_inc))^2))/(ep*cos(ang_inc)+sqrt(ep-(sin(ang_inc))^2)));
