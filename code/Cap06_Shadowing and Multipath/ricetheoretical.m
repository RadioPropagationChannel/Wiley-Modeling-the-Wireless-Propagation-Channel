function [F_r_axis,F_r]=ricetheoretical(A,MP)

warning off

sigma=sqrt(10^(MP/10)/2);
a=10^(A/20);

f_r=[];
f_r_axis=[];


for r=0:0.001:5
    aux = r/(sigma^2)*exp(-(r.^2+a.^2)./(2*sigma^2)).*besseli(0,r.*a/sigma^2);
    if r>0.5 & aux<1e-4, break, end
    f_r=[f_r aux];
    f_r_axis=[f_r_axis r];
end

figure,plot(f_r_axis,f_r)
title('Rice pdf');

% Rice CDF 

F_r=[]; 
F_r_axis=[];

for kk=0:0.001:3
    p= @(r) r./(sigma^2).*exp(-(r.^2+a.^2)./(2.*sigma^2)).*besseli(0,r.*a./sigma^2);
    aux=quad(p,0,kk);
    if 1-aux < 1e-4 & kk > 0.5, break,end 
    F_r=[F_r aux];
    F_r_axis=[F_r_axis kk];
end

figure,plot(F_r_axis,F_r)
title('Rice CDF');
ylabel('probability the abscissa is not exceeded');
xlabel('Amplitude (linear units)');

warning on