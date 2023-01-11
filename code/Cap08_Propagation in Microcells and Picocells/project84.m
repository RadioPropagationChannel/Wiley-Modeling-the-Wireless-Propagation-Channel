%
% project84            Saleh and Valenzuela model
%
%=======================================================================
clear
close all
clc
%=======================================================================
fMHz=1500;  % frequency MHz
dd=30;      % distance in m
nn=4;       % propagagation exponent 
EIRP=1;      % EIRP in dBW

% S&V model parameters

lambda=5;     % ns^-1
LAMBDA=25;   % ns^-1
gamma=29;       % ns
GAMMA=60;       % ns

% time-series generator and Butterworth filter parameters 
fs=20;      % Sampling frequency in Hz 
Wp=0.01;     % normalized freqquency wrt fs/2. Passband
Ws=0.2;     % normalized freqquency wrt fs/2. Stopband
Rp=3;       % Passband spec in dB
Rs=40;      % Stopband spec in dB
leng=500;  % Simulated series length 

%=======================================================================
lambdac=300/fMHz;


% Free spece loss and received power at 1 m from Tx

Lfs1m=32.4+20*log10(0.001)+20*log10(fMHz);  % free space loss at 1 m distance
Pr1m=EIRP-Lfs1m;           % received power at 1 m assuming Gr=0 dBi;
Pr=Pr1m-10*nn*log10(dd)    % received power under an nn propoagation law
Lfsddm=32.4+20*log10(dd/1000)+20*log10(fMHz); % free space loss at dd m
Prfs=EIRP-Lfsddm           % received power under fs conditions (dBW)
prRelativedB=Pr-Prfs       % received power relative to free space (dB)


% =======================================================================
%GENERATE CLUSTER ARRIVALS
% =======================================================================
NClusterArrivals=100;   % this is a possible maximum, not used in the end

TT=genExponential(LAMBDA,NClusterArrivals);    % generate exponential inter-arrivals
TT=[0;TT];          % first arrival is at 0 ns
TClusters=cumsum(TT);      
prClusters=ones(length(TT),1); % we start off assuming all cluster powers are equal with value one 

figure,stem(TClusters,prClusters)
xlabel('Times of arrival of clusters (ns)')
ylabel('Uniform cluster powers')

% Cluster power decay according to GAMMA

prClusters=prClusters.*exp(-TClusters/GAMMA);   % apply power decay rate

DeltaPower=prRelativedB-10*log10(sum(prClusters))  % compute difference with objective 
prClusters=prClusters*10.^(DeltaPower/10);         % correct powers to fit objective

figure,stem2D(TClusters,10*log10(prClusters),-80)
xlabel('Times of arrival of clusters (ns)')
ylabel('Cluster powers (dB)')
grid

figure,bar([1:length(prClusters)],10*log10(cumsum(prClusters)))
xlabel('Number of clusters accumulated');
ylabel('Accumulated relative power (dB)');

% select number of significant clusters

stopT=find(10*log10(cumsum(prClusters))>prRelativedB-0.5);  % within 0.5 dB
stopT=min(stopT)          % select first


TClusters=TClusters(1:stopT);
prClusters=prClusters(1:stopT);

figure,stem2D(TClusters,10*log10(prClusters),-80)
xlabel('Times of arrival of clusters (ns)')
ylabel('Cluster powers (dB)')
grid

figure,stem(TClusters,prClusters)
xlabel('Times of arrival of clusters (ns)')
ylabel('Cluster powers')
grid

%=======================================================================
% now the same process has to be done again for each cluster
% within each cluster we need to draw the ray delays and average powers

NClusters=length(TClusters)

prClustersdB=10*log10(prClusters);
NRayArrivals=100;

TimeTotalRays=[];
PowerTotalRays=[];
 
for ii=1: NClusters
   
    %GENERATE RAY ARRIVALS FOR CLUSTER ii

    TT=genExponential(lambda,NRayArrivals);
    TT=[0;TT];          % first arrival is at 0 ns
    TRays=cumsum(TT);      
    prRays=ones(length(TT),1);   

    % Cluster power decay according to gamma
    prRays=prRays.*exp(-TRays/gamma);

    DeltaPower=prClustersdB(ii)-10*log10(sum(prRays));
    prRays=prRays*10.^(DeltaPower/10);

    % select number of significant rays

    stopT=find(10*log10(cumsum(prRays))>prClustersdB(ii)-0.5);  % within 0.5 dB
    stopT=min(stopT);          % select first

    TRays=TRays(1:stopT);
    prRays=prRays(1:stopT);
    
%     TRays=TRays+TClusters(ii);
    TimeTotalRays=[TimeTotalRays; TRays+TClusters(ii)];
    PowerTotalRays=[PowerTotalRays; prRays];
      
end   % finish all clusters

figure,stem2D(TimeTotalRays,10*log10(PowerTotalRays),-80)
xlabel('Times of arrival of rays (ns)')
ylabel('Ray powers (dB)')
grid

figure,stem(TimeTotalRays,PowerTotalRays)
xlabel('Times of arrival of rays (ns)')
ylabel('Ray powers')
grid

% calculate average delay DD and delay spread SS

disp('average delay DD and delay spread SS')
[DD,SS]=PDPparameters(TimeTotalRays,PowerTotalRays)

% =====================================================================
sigmaRay=0.5*sqrt(PowerTotalRays);
Ntaps=length(PowerTotalRays)
TDLseries=[];
for ii=1:Ntaps
    [t_axis,ryf]=genRayleighFiltered(leng,fs,sigmaRay(ii),Wp,Ws,Rp,Rs);
    TDLseries=[TDLseries ryf];
end
otheraxis=ones(length(t_axis),1);

figure, hold
for ii=1:Ntaps
     plot3(otheraxis.*TimeTotalRays(ii),t_axis, 20*log10(abs(TDLseries(:,ii))),'k')
end
grid
view(3)
ylabel('Time (s)')
zlabel('Reletive signal level (dB)')
xlabel('Excess delay (ns)')


