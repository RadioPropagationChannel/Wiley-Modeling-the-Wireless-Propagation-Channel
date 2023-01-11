function [hry1t]=rayHIST(hry2,sigma);

salto=hry2(2)-hry2(1);
extremo1=hry2-salto/2;
extremo2=hry2+salto/2;

hry1t=exp(-(extremo1.^2)./(2*sigma^2))-exp(-(extremo2.^2)./(2*sigma^2));
