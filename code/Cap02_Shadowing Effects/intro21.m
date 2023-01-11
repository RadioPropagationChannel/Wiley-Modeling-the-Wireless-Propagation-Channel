% INITIALIZE ===========================================================

close all
clear

%  Fresnel integrals ==================================================

v=[-5:0.01:5];
C=mfun('FresnelC',v);
S=mfun('FresnelS',v);

figure, plot(v, C, 'k',v,S,'k--','LineWidth',2)
grid
title('Fresnel Cosine and Sine integrals')
xlabel('Input variable, v')
ylabel('Cosine and Sine integrals')
legend('C(v)','S(v)')

figure,plot(C,S,'k','LineWidth',2)
axis([-1 1 -1 1])
axis square
grid
xlabel('C(v)')
ylabel('S(v)')
title('Cornu spiral')

% Knife-edge diffraction

KE=0.5*(1+j)*((0.5-C)-j*(0.5-S));

figure,plot(v,abs(KE),'k','LineWidth',2)
grid
title('Knife-edge diffraction')
ylabel('Received signal rel. to direct signal (linear units)')
xlabel('Normalized obstruction, v')


figure,plot(v,20*log10(abs(KE)),'k','LineWidth',2)
grid
title('Knife-edge diffraction')
ylabel('Received signal rel. to direct signal (dB)')
xlabel('Normalized obstruction, v')


figure,plot(v,angle(KE),'k','LineWidth',2)
title('Knife-edge diffraction')
ylabel('Phase (Rad)')
xlabel('Normalized obstruction, v')
