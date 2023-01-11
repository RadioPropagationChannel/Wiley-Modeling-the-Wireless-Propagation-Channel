function [D,S]=PDPparameters(delays,powers)

D=sum(delays.*powers)./sum(powers);
RMS2=sum(delays.^2.*powers)./sum(powers);
S=sqrt(RMS2-D.^2);