

% Comparison with theoretical CDF ==================================00

[F_r_axis,F_r]=suzukitheoretical(M,S);

figure,plot(20*log10(xSUZ),ySUZ,'g.-',20*log10(F_r_axis),F_r,'r')
title('CDF: Suzuki series and theoretical') 
xlabel('Signal level (dB/LOS)')
ylabel('Porbability the abscissa is not exceeded')
legend('Time-series','Theoretical')

figure,plot((xSUZ),ySUZ,'g.-',(F_r_axis),F_r,'r')
title('CDF: Suzuki series and theoretical') 
xlabel('Signal level (linear units)')
ylabel('Porbability the abscissa is not exceeded')
legend('Time-series','Theoretical')
