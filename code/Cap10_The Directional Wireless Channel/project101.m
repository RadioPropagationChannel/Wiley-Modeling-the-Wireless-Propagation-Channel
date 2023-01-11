%
% project101     Generate synthetic point scatterer scenario   
%                based on extended Saleh and Valenzuela model
%
%=======================================================================
clear
close all
clc
%=======================================================================
fMHz=1500;  % frequency MHz
dd=50;      % distance in m
nn=4;       % propagagation exponent    
EIRP=1      % EIRP in dBW                  
cc=3e8;     % speed of light
% NRaysCluster=10;    % No of rays per cluster
sigmaLaplace=5;   % Laplace standard deviation in degree

% S & V model parameters

lambda=5;     % ns^-1
LAMBDA=25;   % ns^-1
gamma=29;       % ns
GAMMA=60       % ns


%=======================================================================
% secondary parameters
lambdac=300/fMHz;
dirRayDelay=1e9*dd/cc;   % in ns


% Free spece loss and received power at 1 m from Tx

Lfs1m=32.4+20*log10(0.001)+20*log10(fMHz);  % free space loss at 1 m distance
Pr1m=EIRP-Lfs1m;           % received power at 1 m assuming Gr=0 dBi;
Pr=Pr1m-10*nn*log10(dd)    % received power under an nn propoagation law
Lfsddm=32.4+20*log10(dd/1000)+20*log10(fMHz); % free space loss at dd m
Prfs=EIRP-Lfsddm           % received power under fs conditions (dBW)
prRelativedB=Pr-Prfs       % received power relative to free space (dB)

%========================================================================

NClusterArrivals=100;

%========================================================================
%GENERATE CLUSTER ARRIVALS
% =======================================================================

TT=genExponential(LAMBDA,NClusterArrivals);
TT=[0;TT];          % first arrival is at 0 ns
TClusters=cumsum(TT);      
prClusters=ones(length(TT),1);   

figure,stem(TClusters,prClusters,'k')
xlabel('Times of arrival of clusters (ns)')
ylabel('Uniform cluster powers')

% Cluster power decay according to GAMMA
prClusters=prClusters.*exp(-TClusters/GAMMA);

DeltaPower=prRelativedB-10*log10(sum(prClusters))
prClusters=prClusters*10.^(DeltaPower/10);

figure,stem2D(TClusters,10*log10(prClusters),-80)
xlabel('Times of arrival of clusters (ns)')
ylabel('Cluster powers')
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

% =======================================================================
% up to now we have our intermediate resuls: 
% TClusters contains delays refferred to first arrival (not direct ray)
% prClusters
% the distance between BS and MS is dd
% we have to calculate for each delay the corresponding ellipse parameters
% aa and bb

TClustersAbsolute=TClusters+dirRayDelay;

NClusters=length(TClusters)
xCluster=[];
yCluster=[];
aaCluster=[];
bbCluster=[];
eccentricityCluster=[];
thetaCluster=[];

figscenario=figure, hold 
plot(-dd,0,'k^',dd,0,'k^')
for ii=1:NClusters                % generate cluster nominal locations
    aa=cc*TClustersAbsolute(ii)*1e-9/2+dd;
    bb=sqrt(aa^2-dd^2);
    eccentricity=dd/aa;
    theta=rand(1,1)*360;
    xx=(dd+aa*cosd(theta))/(1+eccentricity*cosd(theta));
    yy=(aa-eccentricity*xx)*sind(theta);
    
    aaCluster=[aaCluster; aa];     % storing obtained values for later
    bbCluster=[bbCluster; bb];
    eccentricityCluster=[eccentricityCluster; eccentricity];
    thetaCluster=[thetaCluster; theta];
    xCluster=[xCluster; xx];
    yCluster=[yCluster; yy];
    
    plot(xx,yy,'k*')               % Plot point-scatterer
    
    thetaall=[0:360];
    xxall=(dd+aa.*cosd(thetaall))./(1+eccentricity.*cosd(thetaall))
    yyall=(aa-eccentricity.*xxall).*sind(thetaall)
    plot(xxall,yyall,'k:')         % Plot ellipse
end

% ======================================================================

%=======================================================================
% now this has to be done again for each cluster
% within each cluster we need to draw the ray delays and average powers
 
prClustersdB=10*log10(prClusters);
NRayArrivals=100;

TimeTotalRays=[];
PowerTotalRays=[];
RaysInCluster=[];     % to keep track of how many rays there are in each cluster

for ii=1: NClusters    %stopT      % ii points consecutively at all clusters
   
    %GENERATE RAY ARRIVALS WITH CLUSTER ii

    TT=genExponential(lambda,NRayArrivals);
    TT=[0;TT];          % first arrival is at 0 ns
    TRays=cumsum(TT);      
    prRays=ones(length(TT),1);   

    % Cluster power decay according to gamma
    prRays=prRays.*exp(-TRays/gamma);          % < .........?????

    DeltaPower=prClustersdB(ii)-10*log10(sum(prRays));
    prRays=prRays*10.^(DeltaPower/10);

    % select number of significant rays

    stopT=find(10*log10(cumsum(prRays))>prClustersdB(ii)-0.5);  % within 0.5 dB
    stopT=min(stopT)          % select first
    RaysInCluster=[RaysInCluster; stopT];

    TRays=TRays(1:stopT);
    prRays=prRays(1:stopT);
    
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

% =======================================================================
% up to now we have our intermediate resuls in 
% TimeTotalRays contains delays refferred to first arrival (not direct ray)
% PowerTotalRays
% the distance between BS and MS is dd
% we have to calculate for each delay the corresponding ellipse parameters
% aa and bb

TimeTotalRaysAbsolute=TimeTotalRays+dirRayDelay;   % absolute delays in ns

% =======================================================================
% generation of actual point scatterer locations  
% =======================================================================
 
xPointScats=[];
yPointScats=[];
AoA=[];

for ii=1:NClusters
    for jj=1:RaysInCluster(ii)
        DeltaTheta=genLaplacian(sigmaLaplace);   % draw Laplace deviation from nominal cluster angle
        newTheta=thetaCluster(ii)+DeltaTheta;
          
        aa=cc*TimeTotalRaysAbsolute(sum(RaysInCluster(1:ii-1))+jj)*1e-9/2+dd;     % delays are in ns are converted to s
        bb=sqrt(aa^2-dd^2);
        eccentricity=dd/aa;

        xx=(dd+aa*cosd(newTheta))/(1+eccentricity*cosd(newTheta));
        yy=(aa-eccentricity*xx)*sind(newTheta); 
                
        AoA=[AoA; newTheta];
        xPointScats=[xPointScats; xx];
        yPointScats=[yPointScats; yy];
    end
end

figure(figscenario)
plot(xPointScats,yPointScats, 'ko')
xlabel('Distance (m)');
ylabel('Distance (m)');
axis equal


figure,stem3(AoA,TimeTotalRays,PowerTotalRays)
xlabel('AoA (deg)')
ylabel('time (ns)')
zlabel('power');