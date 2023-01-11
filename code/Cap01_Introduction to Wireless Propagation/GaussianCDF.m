function [CDFy]=GaussianCDF(M,S,CDFx)

k=(CDFx-M)/S;

CDFy=1-0.5*erfc(k/sqrt(2));