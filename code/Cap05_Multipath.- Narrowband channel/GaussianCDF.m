function [CDFy]=gaussianCDF(mean,sigma,axisCDFx)
CDFy=1-0.5*erfc((axisCDFx-mean)/sigma/sqrt(2));
