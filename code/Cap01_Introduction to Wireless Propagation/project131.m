%
% project131
%
% INITIALIZE =========================================================
clear 
close all
clc
% ====================================================================

load series131                 % Rayleigh case 

timeaxis=series131(:,1);       % timeaxis in s
I=real(series131(:,2));        % signal levels are referred to the
Q=imag(series131(:,2));        % direct's signal level

% ======================================================================
figure,plot(timeaxis,I,'k',timeaxis,Q,'k.-')
title('Normalized received signal, real part and imaginary parts')
ylabel('Normalized received signal, real and imaginary parts (lin. units)')
xlabel('Elapsed time (s)')
legend('Real part','Imginary part')

figure,plot(I,Q,'k')
title('In-phase v. quadrature plot. Rayleigh case')
ylabel('Quadrature component')
xlabel('In-phase component')
axis([-1 1 -1 1])

figure,plot(timeaxis,abs(I+j.*Q),'k')
title('Normalized received signal, magnitude')
ylabel('Normalized received signal, magnitude (lin. units)')
xlabel('Elapsed time (s)')

figure,plot(timeaxis,20*log10(abs(I+j.*Q)),'k')
title('Normalized received signal, magnitude')
ylabel('Normalized received signal, magnitude (dB)')
xlabel('Elapsed time (s)')

figure,plot(timeaxis,angle(I+j.*Q),'k')
title('Received signal, modulo-\pi phase')
ylabel('Received signal, modulo-\pi phase (radians)')
xlabel('Elapsed time (s)')

figure,plot(timeaxis,unwrap(angle(I+j.*Q)),'k')
title('Received signal, absolute phase')
ylabel('Received signal, phase (radians)')
xlabel('Elapsed time (s)')


