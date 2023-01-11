% INITIALIZE ===========================================================

close all
clear

% INPUT PARAMETERS ==================================================

fMHz=2e3;  % frequency in MHz


% SECONDARY PARAMETERS ==============================================
f=fMHz*1e6;
lambdac=300/fMHz;
kc=2*pi/lambdac;

stepAperture=lambdac/2;
maxAperture=300; 
Na=200;
psi=100; % Obstacles height (All have the same height).

% GEOMETRIC INPUTS =================================================

xt=0;   % m
yt=100;      % m   (0 m indicates same height as screen)

stepRx = 1;
yr=[0:stepRx:150]';    % m   SAMPLING AT THE RECEIVER SIDE



% SAMPLING POINTS ALONG THE APERTURE =================================
xa=[1000 2000 3000 4000 5000]; % There are 4 obstacles, the last one point is the reception point for the case with 4 obtacles.
ya=[psi:stepAperture:maxAperture]; % All obstacles have the same height.

% Plot of the geometry.
figure,plot(repmat(xt,1,length(0:yt)),0:yt,'g',xt,yt,'.g'),hold on
for m=1:length(xa)-1
    plot(repmat(xa(m),1,length(0:psi)),0:psi,'LineWidth',2);
end
plot(repmat(xa(end),1,length(yr)),yr,'.g')
hold off
ylabel('Height (m)')
xlabel('Distance (m)');
axis([xt-100 xa(end)+100 0 yr(end)+10])



% Window along the aperture.
w=triang_win(2*length(find(ya>Na)));
wa=[ones(1,length((find(ya<=Na)))) w(floor(length(w)/2)+1:end)'];

% Tx-Aperture1 side calculations ==========================================
DistTxApertureX=xa(1)-xt;
DistTxApertureY=ya-yt;
% DistTxAperture=sqrt(DistTxApertureX.^2+DistTxApertureY.^2);

RTxAperture=DistTxApertureX+((DistTxApertureY.^2)/(2*DistTxApertureX));

Efs1=exp(-j*kc*RTxAperture)/DistTxApertureX.*wa;
Edant=Efs1.';

for m=1:length(xa)-1
    % Rx
    xr=xa(m+1);                % m
    
    % Aperture-Rx side calculations ==========================================
    DistApertureRxX=xr-xa(m);
    
    DistApertureijY=repmat(ya',1,length(ya))-repmat(ya,length(ya),1);
    DistApertureRxY=repmat(yr,1,length(ya))-repmat(ya,length(yr),1);

    DistTxRx=(xr-xt)+((yr-yt).^2)/(2*(xr-xt));
    
    RApertureij=DistApertureRxX+((DistApertureijY.^2)/(2*DistApertureRxX));
    RApertureRx=DistApertureRxX+((DistApertureRxY.^2)/(2*DistApertureRxX));
    
    Efs=exp(-j*kc*DistTxRx)/(xr-xt);
    
    Fd=sqrt(kc*(xa(m)-xt)/(2*pi*j*(xr-xt)*(xr-xa(m))));
    Ed=Fd*exp(-j*kc*RApertureRx).*repmat(Edant.',length(yr),1)*stepAperture;
    Ed=sum(Ed,2);
  
    Edant=Fd*exp(-j*kc*RApertureij).*repmat(Edant.',length(ya),1)*stepAperture;
    Edant=sum(Edant,2).*wa';
    
    figure,plot(yr,abs(Ed./Efs),'k')
    title(['Obstacle # ' num2str(m)]);
    ylabel('Normalized field strength')
    xlabel('z(m)');
end
