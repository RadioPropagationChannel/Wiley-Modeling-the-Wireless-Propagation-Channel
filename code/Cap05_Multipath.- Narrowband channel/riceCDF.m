function [y]=riceCDF(a,sigma,x);
y=quad('ricetheoretical',0,x,1e-6,0,a,sigma);