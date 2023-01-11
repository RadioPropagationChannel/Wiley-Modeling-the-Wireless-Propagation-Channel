function [pdfy]=Rayleighpdf(sigma,axispdfx)
pdfy=(axispdfx./sigma).*exp(-(axispdfx.^2)./(2*sigma.^2));
