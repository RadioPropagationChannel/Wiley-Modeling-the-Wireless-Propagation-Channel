function [F_r_axis,F_r]=lootheoretical(M,S,MP)
format long
warning off

sigma=sqrt(10^(MP/10)/2);


% 3D function.
X=[];
tam_a=10;
tam_r=10;


ii=1;
jj=1;
for r=0:.1:tam_r
    for a=0:.1:tam_a
        X(jj,ii)=1/a*exp(-(20*log10(a)-M)^2/(2*S^2))*exp(-(r^2+a^2)/(2*sigma^2))*besseli(0,r*a/sigma^2);
        ii=ii+1;
    end
    ii=1;
    jj=jj+1;
end
r=0:.1:tam_a;
a=0:.1:tam_r;

figure, surf(a,r,X);
[ro,co]=find(X>=10^-3);
title('Probability density function of Loo distribution')
xlabel('a');
ylabel('r');
zlabel('Probability')

% Funcion Distribucion Loo
paso=(r(max(ro))-r(min(ro)))/50;
F_r_axis=[r(min(ro)):paso:r(max(ro))];
lon=length(F_r_axis);
F_r=ones(1,lon);
j=1;
for k=r(min(ro)):paso:r(max(ro))
   p = @(a,r) 8.686*r/(sigma^2*S*sqrt(2*pi)).*1./a.*exp(-(20.*log10(a)-M).^2./(2*S^2)).*exp(-(r.^2+a.^2)./(2*sigma^2)).*besseli(0,r.*a/sigma^2);
   F_r(j)=dblquad(p,a(min(co)),a(max(co)),r(min(ro)),k);
   j=j+1;
end

figure,hold on
plot(F_r_axis,F_r,'r')
xlabel('r');
ylabel('Probability the abscissa is not exceeded');
title('CDF of Loo distribution');

warning on
format short