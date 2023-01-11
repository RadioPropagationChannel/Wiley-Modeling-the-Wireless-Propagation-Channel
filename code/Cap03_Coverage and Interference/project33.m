%
% project33    
%
%=========================================================================
close all
clear
clc
%=========================================================================
Lcorr=30;   % correlation distance in meters
ls=5;       % sampling spacing in meters 
SL=10;       % location variability in dB

maxradius=500;  % maximum radious of simulated field
AA=137;      % Loss at 1 km (propagation model) in dB
nn=3.5;     % Propagation exponent
EIRP=0;     % EIRI in dBm
Gr=0;       % Receiver gain in dBi

ThresHold=-110;   % Threshold in dBm
fringe=100; % converage contour
%=========================================================================
interprate=round(Lcorr/ls);
Lcorr=ls*interprate;
Nsamples=round(maxradius/Lcorr);

%=========================================================================
map1=randn(2*Nsamples,2*Nsamples)*SL;
xaxis=[-Nsamples:Nsamples-1]*Lcorr;
yaxis=[-Nsamples:Nsamples-1]*Lcorr;
figure,mesh(xaxis,yaxis,map1)
xlabel('distance (m)')
ylabel('distance (m)')

% Interpolated distances
xaxisinterp=[-Nsamples:1/interprate:Nsamples-1]*Lcorr;
yaxisinterp=[-Nsamples:1/interprate:Nsamples-1]*Lcorr;

map1interp=interp2(xaxis,yaxis,map1,xaxisinterp,yaxisinterp','spline');
figure,mesh(xaxisinterp,yaxisinterp,map1interp)
xlabel('distance (m)')
ylabel('distance (m)')

% calculate distances 
map2interp=zeros(length(xaxisinterp),length(yaxisinterp));
for ii=1:length(xaxisinterp),
    for jj=1:length(yaxisinterp),
        map2interp(jj,ii)=sqrt(xaxisinterp(ii).^2+yaxisinterp(jj).^2);
    end
end
[y0,x0]=find(map2interp==0);    % locate if distance is zero
map2interp(y0,x0)=1;            % change distance equal zero to avoid log(0)
figure,mesh(xaxisinterp, yaxisinterp,map2interp)
xlabel('distance (m)')
ylabel('distance (m)')

map3interp=EIRP-AA-10*nn*log10(map2interp/1000)+Gr;   % distance in m need km in formula
figure,mesh(xaxisinterp,yaxisinterp,map3interp)
xlabel('distance (m)')
ylabel('distance (m)')

map4interp=map1interp+map3interp;
figure,mesh(xaxisinterp,yaxisinterp,map4interp)
xlabel('distance (m)')
ylabel('distance (m)')

map5interp=map4interp-ThresHold;
figure,mesh(map5interp)


[xx,yy]=find(map5interp>0);
map6interp=zeros(length(xaxisinterp),length(yaxisinterp));
for ii=1:length(xx)
    indexx=xx(ii);
    indexy=yy(ii);
    map6interp(indexx,indexy)=1;
end
figure,contour(xaxisinterp, yaxisinterp, map6interp)
xlabel('distance (m)')
ylabel('distance (m)')
title('Coverage cotours')

% area coverage calculation
coverageCount=0;
[xx,yy]=find(map2interp<fringe);
for ii=1:length(xx)
    indexx=xx(ii);
    indexy=yy(ii);
    if map6interp(indexx,indexy)==1,
        coverageCount=coverageCount+1;
    end
end
AreaCovProb=coverageCount/length(xx);

message='Area coverage probability: ';
disp([message num2str(AreaCovProb)])

% fringe coverage calculation
coverageCount=0;
[xx,yy]=find(map2interp<fringe+5 & map2interp>fringe-5);
length(xx)
for ii=1:length(xx)
    indexx=xx(ii);
    indexy=yy(ii);
    if map6interp(indexx,indexy)==1,
        coverageCount=coverageCount+1;
    end
end
FringeCovProb=coverageCount/length(xx);

message='Fringe coverage probability: ';
disp([message num2str(FringeCovProb)])

