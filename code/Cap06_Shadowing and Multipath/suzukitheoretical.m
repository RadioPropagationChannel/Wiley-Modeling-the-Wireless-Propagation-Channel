function      [F_r_axis,F_r]=suzukitheoretical(M,S)

warning off

% calculate pdf

f_r=[];
f_r_axis=[];
for r=0:.01:6
    p = @(sigma) 8.686*r/(S*sqrt(2*pi)).*(1./sigma.^3).*exp(-(20.*log10(sigma)-M).^2./(2*S^2)).*exp(-r^2./(2.*sigma.^2)); 
    aux=quad(p,0,3);
    if r>0.5 & aux<1e-4,break,end
    f_r=[f_r aux];
    f_r_axis=[f_r_axis r];
end

figure, plot(f_r_axis,f_r)
xlabel('signal amplitude (lin. units)')
title('Suzuki pdf');

% Calculate CDF 

F_r=[];
F_r_axis=[];

for kk=0:0.01:6
%     kk
    p = @(sigma,r) 8.686*r/(S*sqrt(2*pi)).*(1./sigma.^3).*exp(-(20.*log10(sigma)-M).^2./(2*S^2)).*exp(-r^2./(2.*sigma.^2)); 
    aux=dblquad(p,0,3,0,kk);
    if 1-aux < 1e-4 & kk > 0.5, break,end 
    F_r=[F_r aux];
    F_r_axis=[F_r_axis kk];
end

figure, plot(F_r_axis,F_r)
title('Suzuki CDF');
ylabel('Probability the abscissa is not exceeded');
xlabel('Amplitude (linear units)');

figure, plot(20*log10(F_r_axis),F_r)
title('Suzuki CDF');
ylabel('Probability the abscissa is not exceeded');
xlabel('Amplitude (dB/LOS)');

warning on