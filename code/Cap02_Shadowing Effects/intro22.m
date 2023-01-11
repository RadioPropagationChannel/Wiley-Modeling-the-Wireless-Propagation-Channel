% Comparation between knife edge difraction losses and its approximation
close all;
clear;

v=-3:.1:3;

ind=find(v>=-0.7);

Ld=zeros(1,length(v));
Ld(ind)=6.9+20*log10(sqrt(((v(ind)-0.1).^2)+1)+v(ind)-0.1);

L_fresnel=20*log10(abs((1+j)/2.*((0.5-mfun('FresnelC',v))-j*(0.5-mfun('FresnelS',v)))));

figure,plot(v,L_fresnel,'k',v,-Ld,'k:','lineWidth',2),legend('Fresnel','Approximation')
xlabel('\nu parameter')
ylabel('Relative field strength level (dB)')