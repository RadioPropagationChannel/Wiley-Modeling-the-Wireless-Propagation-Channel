%
% project34    (multiple interfererers)
%
%========================================================================
close all
clear
clc
%========================================================================
Ncases=5000;    % Montercarlo. Number of draws 
%========================================================================

MC=-80;
MI1=-110;
MI2=-120;
MI3=-125;
MI4=-108;
MI5=-90;

SC=8;
SI1=7;
SI2=6;
SI3=8;
SI4=10;
SI5=7;
%======================================================================

rho=[1.0 0.5 0.7 0.0 0.0 0.4
     0.5 1.0 0.0 0.0 0.0 0.3
     0.7 0.0 1.0 0.0 0.0 0.0
     0.0 0.0 0.0 1.0 0.5 0.0
     0.0 0.0 0.0 0.5 1.0 0.3
     0.4 0.3 0.0 0.0 0.3 1.0];

%=======================================================================
C=randn(1,Ncases)*SC+MC;
I1=randn(1,Ncases)*SI1+MI1;
I2=randn(1,Ncases)*SI2+MI2;
I3=randn(1,Ncases)*SI3+MI3;
I4=randn(1,Ncases)*SI4+MI4;
I5=randn(1,Ncases)*SI5+MI5;
%=======================================================================
i1=10.^(I1/10);
i2=10.^(I2/10);
i3=10.^(I3/10);
i4=10.^(I4/10);
i5=10.^(I5/10);
iTot=i1+i2+i3+i4+i5;
ITot=10*log10(iTot);
%========================================================================
CIR=C-ITot;
figure,plot([1:Ncases],CIR,'k',[1:Ncases], C, 'k--',[1:Ncases], ITot,'k:')
xlabel('Draw number')
ylabel('Level (dBm) uncorrelated case')
legend('CIR','C','ITot');

figure,plot([1:Ncases],ITot,'k',[1:Ncases], I1, 'k:',[1:Ncases], I2,'k:',...
    [1:Ncases], I3, 'k:',[1:Ncases], I4,'k:',[1:Ncases], I5, 'k:')
xlabel('Draw number')
ylabel('Level (dBm) uncorrelated case')
legend('ITot', 'I1','I2','I3','I4','I5');

message='CIR mean: ';
disp([message num2str(mean(CIR)) ' dB. Uncorrelated case'])
message='CIR std: ';
disp([message num2str(std(CIR)) ' dB. Uncorrelated case'])

[CDFxu,CDFyu]=fCDF(CIR);
[CDFyuTH]=GaussianCDF(mean(CIR),std(CIR),CDFxu);

figure,plot(CDFxu,CDFyu,'k:',CDFxu,CDFyuTH,'k')
xlabel('CIR (dB)')
ylabel('Probabilty of not exceeding the abscissa')
title('Simulated and theoretical CDF of CIR. Uncorrelated case')
legend('Simulated CIR','Theoretical CIR')
%=========================================================================

C=(C-MC)/SC;
I1=(I1-MI1)/SI1;
I2=(I2-MI2)/SI2;
I3=(I3-MI3)/SI3;
I4=(I4-MI4)/SI4;
I5=(I5-MI5)/SI5;


cholcorr=chol(rho)';
CIcorr=cholcorr*[C;I1;I2;I3;I4;I5];

C=CIcorr(1,:)*SC+MC;
I1=CIcorr(2,:)*SI1+MI1;
I2=CIcorr(3,:)*SI2+MI2;
I3=CIcorr(4,:)*SI2+MI3;
I4=CIcorr(5,:)*SI4+MI4;
I5=CIcorr(6,:)*SI5+MI5;

%=========================================================================
i1=10.^(I1/10);
i2=10.^(I2/10);
i3=10.^(I3/10);
i4=10.^(I4/10);
i5=10.^(I5/10);
iTot=i1+i2+i3+i4+i5;
ITot=10*log10(iTot);

mean(ITot)
std(ITot)

%========================================================================
CIR=C-ITot;
figure,plot([1:Ncases],CIR,'k',[1:Ncases], C, 'k--',[1:Ncases], ITot,'k:')
xlabel('Draw number')
ylabel('Level (dBm) correlated case')
legend('CIR','C','ITot');

figure,plot([1:Ncases],ITot,'k',[1:Ncases], I1, 'k:',[1:Ncases], I2,'k:',...
    [1:Ncases], I3, 'k:',[1:Ncases], I4,'k:',[1:Ncases], I5, 'k:')
xlabel('Draw number')
ylabel('Level (dBm) correlated case')
legend('ITot', 'I1','I2','I3','I4','I5');

message='CIR mean: ';
disp([message num2str(mean(CIR)) ' dB. Correlated case'])
message='CIR std: ';
disp([message num2str(std(CIR)) ' dB. Correlated case'])

[CDFxc,CDFyc]=fCDF(CIR);
[CDFycTH]=GaussianCDF(mean(CIR),std(CIR),CDFxc);

figure,plot(CDFxc,CDFyc,'k',CDFxc,CDFycTH,'k:')
xlabel('CIR (dB)')
ylabel('Probabilty of not exceeding the abscissa')
title('Simulated and theoretical CDF of CIR. Correlated case')
legend('Simulated CIR','Theoretical CIR')

figure,plot(CDFxc,CDFyc,'k',CDFxu,CDFyu,'k:')
xlabel('CIR (dB)')
ylabel('Probabilty of not exceeding the abscissa')

legend('CIR correlated case','CIR uncorrelated case')