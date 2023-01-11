

% Comparison with theoretical CDF ==================================00

[F_r_axis,F_r]=lootheoretical(M,S,MP)

figure,plot(20*log10(xLOO),yLOO,'g.-',20*log10(F_r_axis),F_r,'r')
title('CDF: Loo series and theoretical') 
xlabel('Signal level (dB/LOS)')
ylabel('Porbability the abscissa is not exceeded')
legend('Time-series','Theoretical')

figure,plot((xLOO),yLOO,'g.-',(F_r_axis),F_r,'r')
title('CDF: Loo series and theoretical') 
xlabel('Signal level (linear units)')
ylabel('Porbability the abscissa is not exceeded')
legend('Time-series','Theoretical')
