function [CDFy]=RayleighCDF(sigma,axisCDFx)
CDFy=1-exp(-(axisCDFx.^2)./(2*sigma.^2));
