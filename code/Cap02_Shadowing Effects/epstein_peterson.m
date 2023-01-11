function [Ld,nu] = epstein_peterson(pos,h,n,lambda)
% 
% function [Ld,nu] = epstein_peterson(pos,h,n,lambda)
% It calculates the field for each path of the Epstein & Peterson model.
%
% pos   :   transmitter, obstacles and receiver x-axis positions vector (m).
% h     :   transmitter, obstacles and receiver heights vector (m).
% n     :   obstacles number (from 1 to 3).
% lambda:   wavelength (m).
%
% Ld    :   total losses given by Epstein & Peterson model (dB). 
% nu    :   difraction parameter.

if n<1
    Ld=0;
    nu=NaN;
    return;
end

Ldi=[];
for m=1:n
    d1=pos(m+1)-pos(m);
    d2=pos(m+2)-pos(m+1);
    hobs=h(m+1)-((d1*(h(m+2)-h(m))/(pos(m+2)-pos(m)))+h(m));
    nu=hobs*sqrt((2/lambda)*(d1+d2)/(d1*d2));
   
    if (nu<-0.7)
        Ldi(m)=0;
    else
        Ldi(m)=6.9+20*log10(sqrt(((nu-0.1)^2)+1)+nu-0.1);
    end
end
Ld=sum(Ldi);