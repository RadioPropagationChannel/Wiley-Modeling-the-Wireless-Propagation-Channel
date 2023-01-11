%
% project314     (two BS and handover)
%
%========================================================================
clear
close all
clc
%========================================================================
fMHz=2000;     % Frequency in MHz
Lcorr=30;       % Correlation distance in meters 

EIRP=30;        % EIRP dBm
AA=100;        % Loss at 1 km 
n=3.6;          % Propagation exponent
SS=7;           % Location variability in dB

%=======================================================================
%lambdac=300/fMHz;   % wavelength in m
ds=1;               % sample spacing in m  
Distance= 6;       % simulated distance (km)
InterpRate=round(Lcorr/ds);   % Interpolation rate 
Lcorr=InterpRate*ds;       % slightly correct Lcorr to make it a multiple of ds

Nsamples=round(Distance*1000/Lcorr);

BSseparation=Distance;   % separation between BSs in km
HO_ThresHold=-85;  % Handover threshild dBm
HO_Margin=15;      % Handover margin dB

% Create Gaussian series ================================================

d_axis1=[0:Nsamples-1]'*Lcorr/1000;

warning off
veryslowVars=EIRP-AA-10*n*log10(d_axis1);
veryslowVars(1)=veryslowVars(2);        % replace veryslowVars(1)=inf
warning on

slowVars=randn(Nsamples,1)*SS;             % 1st draw for BS1
R1=veryslowVars+slowVars;                  % Received power from BS1 

slowVars=randn(Nsamples,1)*SS;             % 1st draw for BS2
R2=veryslowVars+slowVars;                  % Received power from BS2 


%=======================================================================

d_axis2=([0:1/InterpRate:Nsamples])*Lcorr/1000;

Rinterpolated1=interp1(d_axis1,R1,d_axis2,'spline');
Rinterpolated2=interp1(d_axis1,R2,d_axis2,'spline');

d_axis2BS1=d_axis2;                      % Distance axis for BS1 signal (same as before)
d_axis2BS2=max(d_axis2BS1)-d_axis2BS1;   % Distance axis for BS2 signal 
d_axis2BS2=BSseparation-d_axis2BS1;   % Distance axis for BS2 signal 

figure,plot(d_axis2BS1,Rinterpolated1,'k',d_axis2BS2,Rinterpolated2,'k:', ...
    [d_axis2BS2(end) d_axis2BS1(end)],[HO_ThresHold HO_ThresHold],'k--', ...
    [d_axis2BS2(end) d_axis2BS1(end)],[HO_ThresHold+HO_Margin HO_ThresHold+HO_Margin],'k:')
xlabel('Traversed distance (km)')
ylabel('Received signal (dBm)')


% ========================================================================
% HO algorithm 1
% ========================================================================
BSseries=[1];           % Serving Base Station series
currentBS=1;            % starting with MS hooked up to BS1
RHO=Rinterpolated1(1);  % MS recieving BS1 signal
switchBS=0;             % Swith 

indexBS1=length(Rinterpolated1);
indexBS2=length(Rinterpolated2); 

% for ii=2:length(Rinterpolated1)-1, 
for ii=2:length(Rinterpolated1),     
    if currentBS==1 & switchBS==0, 
        if Rinterpolated1(ii)> HO_ThresHold, 
             RHO=[RHO; Rinterpolated1(ii)];
             currentBS=1;
             BSseries=[BSseries; currentBS];
             switchBS=1; 
        else if Rinterpolated1(ii) < Rinterpolated2(indexBS2-ii+1),
             RHO=[RHO; Rinterpolated2(indexBS2-ii+1)];
             currentBS=2;
             BSseries=[BSseries; currentBS];
             switchBS=1; 
            else
             RHO=[RHO; Rinterpolated1(ii)];
             currentBS=1;
             BSseries=[BSseries; currentBS];
             switchBS=1; 
            end
        end
    end
    if currentBS==2 & switchBS==0, 
        if Rinterpolated2(indexBS2-ii+1)> HO_ThresHold, 
             RHO=[RHO; Rinterpolated2(indexBS2-ii+1)];
             currentBS=2;
             BSseries=[BSseries; currentBS];
             switchBS=1; 
        else if Rinterpolated2(indexBS2-ii+1) < Rinterpolated1(ii),
             RHO=[RHO; Rinterpolated1(ii)];
             currentBS=1;
             BSseries=[BSseries; currentBS];
             switchBS=1; 
            else
             RHO=[RHO; Rinterpolated2(indexBS2-ii+1)];
             currentBS=2;
             BSseries=[BSseries; currentBS];
             switchBS=1; 
            end
        end
    end
    switchBS=0;       % reset swicth BS flag
end

figure,plot(d_axis2BS1(1:indexBS1),BSseries(1:indexBS1),'k')
a=axis;axis([a(1) a(2) 0 3])
xlabel('Traversed distance (m)')
ylabel('No.of BS handling the call')

figure,plot(d_axis2BS1(1:indexBS1),RHO(1:indexBS1),'k',...
    d_axis2BS1(1:indexBS1),Rinterpolated1(1:indexBS1),'k:',...
    d_axis2BS2(1:indexBS2),Rinterpolated2(1:indexBS2),'k--')
xlabel('Traversed distance (km)')
ylabel('Received signal (dBm)')

% =======================================================================
% HO algorithm 2   Including HO Margin
%========================================================================
BSseries=[1];           % Serving Base Station series
currentBS=1;            % starting with MS hooked up to BS1
RHO=Rinterpolated1(1);  % MS recieving BS1 signal
switchBS=0;             % Swith 

indexBS1=length(Rinterpolated1);
indexBS2=length(Rinterpolated2); 

% for ii=2:length(Rinterpolated1)-1, 
for ii=2:length(Rinterpolated1),     
    if currentBS==1 & switchBS==0, 
        if Rinterpolated1(ii)> HO_ThresHold, 
             RHO=[RHO; Rinterpolated1(ii)];
             currentBS=1;
             BSseries=[BSseries; currentBS];
             switchBS=1; 
        else if Rinterpolated1(ii)+HO_Margin < Rinterpolated2(indexBS2-ii+1),
             RHO=[RHO; Rinterpolated2(indexBS2-ii+1)];
             currentBS=2;
             BSseries=[BSseries; currentBS];
             switchBS=1; 
            else
             RHO=[RHO; Rinterpolated1(ii)];
             currentBS=1;
             BSseries=[BSseries; currentBS];
             switchBS=1; 
            end
        end
    end
    if currentBS==2 & switchBS==0, 
        if Rinterpolated2(indexBS2-ii+1)> HO_ThresHold, 
             RHO=[RHO; Rinterpolated2(indexBS2-ii+1)];
             currentBS=2;
             BSseries=[BSseries; currentBS];
             switchBS=1; 
        else if Rinterpolated2(indexBS2-ii+1)+HO_Margin < Rinterpolated1(ii),
             RHO=[RHO; Rinterpolated1(ii)];
             currentBS=1;
             BSseries=[BSseries; currentBS];
             switchBS=1; 
            else
             RHO=[RHO; Rinterpolated2(indexBS2-ii+1)];
             currentBS=2;
             BSseries=[BSseries; currentBS];
             switchBS=1; 
            end
        end
    end
    switchBS=0;       % reset swicth BS flag
end

figure,plot(d_axis2BS1(1:indexBS1),BSseries(1:indexBS1),'k')
a=axis;axis([a(1) a(2) 0 3])
xlabel('Traversed distance (m)')
ylabel('No.of BS handling the call')

figure,plot(d_axis2BS1(1:indexBS1),RHO(1:indexBS1),'k',...
    d_axis2BS1(1:indexBS1),Rinterpolated1(1:indexBS1),'k:',...
    d_axis2BS2(1:indexBS2),Rinterpolated2(1:indexBS2),'k--')
xlabel('Traversed distance (km)')
ylabel('Received signal (dBm)')