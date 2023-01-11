function [C,S]= Fresnel_integrals(x)
% 
% function [C,S]= Fresnel_integrals(x)
% It calculates the Fresnel integrals using the error function.
% 
% x :  abscissa.
% C :  Fresnel cosine for x.
% S :  Fresnel sine for x.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z=sqrt(pi)*(1-j)*x/2;
complex=((1+j)/2)*erfz(z);
C=real(complex);
S=imag(complex);