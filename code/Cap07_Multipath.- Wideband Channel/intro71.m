%
% intro71
%
%========================================================================
clear
close all
clc
%========================================================================
% PDP definition 
%========================================================================
tau=[0 1e3 2e3 5e3]; % 
tau=[0  0.1e3  0.20e3  0.35e3  0.45e3]   % Delays in ns
PDP=[-20 -10 -10 0]; %
PDP=[0 -3 -30 -20 -10]   % Powers in dB

%=========================================================================
% convert powers in dB to linear units
%=========================================================================
pdp=10.^(PDP/10)

figure, stem2D(tau,PDP,-100)
auxx=axis;
axis([auxx(1)-10 auxx(2)+10 auxx(3)-10 auxx(4)+10 ])
xlabel('Excess delay (ns)')
ylabel('Relative level (dB)')

figure, stem(tau,pdp)
auxx=axis;
axis([auxx(1)-0.1 auxx(2)+0.1 auxx(3)-0.1 auxx(4)+0.1 ])
xlabel('Excess delay (ns)')
ylabel('Relative level (linear units)')

%=========================================================================
% PDP parameters: D, S 
%========================================================================
D=sum(pdp.*tau)/sum(pdp)
S=sqrt(sum(pdp.*(tau-D).^2)/sum(pdp))

%=======================================================================
% Now the frequency spaced correlation function, R(Df)
%========================================================================
stepDf=10000;
MaxDf=2e7;
Df=[0:stepDf:MaxDf]';


R=[];
for ii=0:stepDf:MaxDf
    auxR=sum(pdp.*exp(-j*2*pi*ii*tau*1e-9));
    R=[R; auxR];
end
figure,plot(Df,abs(R)/max(abs(R)),'k')
xlabel('Frequency spacing (Hz)')
ylabel('Frequency-spaced correlation')
