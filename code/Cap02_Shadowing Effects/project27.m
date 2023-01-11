close all;
clear;

% Data
fc=2e8;         % (Hz)
c=3e8;          % (m/s)
lambdac=c/fc;    % (m)
kc=2*pi/lambdac;

% Profile geometry.
slope=[0.01 0.01 0.01 0.01];
height=[10 15 20 15];

% Transmitter geometry.
ht=2;           % Transmit antenna height (m)

xt=10;
yt=ht+height(1);
profy=-fliplr(cumsum(ones(1,xt-1)))*slope(1)+height(1);
profy=[profy height(1):-slope(1):0];
xp=xt;

figure,plot(repmat(xt,1,2),[height(1) yt],'g',xt,yt,'g.')
title('Profile')

% Profile generation.

for m=2:length(slope)
    profy=[profy 0:slope(m):height(m)];
    xp=[xp length(profy)];
    profy=[profy height(m)-slope(m):-slope(m):0];
end
profx=0:length(profy)-1;

hold on,plot(profx,profy)


% Receiver geometry.
hr=1.5;      % Receive antenna heigth (m)

% Shadow areas.
txp_slope=(height(2:end)-yt)./(xp(2:end)-xt);

xsh=[];
txp_total=zeros(1,length(profx));
for m=2:length(height)
    x=xp(m):profx(end);
    
    % LOS line between the present peak and the next one.
    txp_line=height(m)+((x-xp(m))*txp_slope(m-1));
    
    height_dr=txp_line-profy(x);
    
    ind_shadow=find(height_dr>=hr);
 
    xsh=[xsh setdiff(x(ind_shadow),xsh)];
    txp_aux=zeros(1,length(profx));
    txp_aux(x(ind_shadow))=txp_line(ind_shadow);
    
    txp_total=max(txp_total,txp_aux);
end

xsh=sort(xsh);
ind_zone=find(diff(txp_total(xsh))>1);
ind_zone=[0 ind_zone];

for m=1:length(ind_zone)
    if m==length(ind_zone)
        plot(repmat(xsh(ind_zone(m)+1),1,2),[profy(xsh(ind_zone(m)+1)) txp_total(xsh(ind_zone(m)+1))],'r',xsh(ind_zone(m)+1:end),...
            txp_total(xsh(ind_zone(m)+1:end)),'r');
    else
        plot(repmat(xsh(ind_zone(m)+1),1,2),[profy(xsh(ind_zone(m)+1)) txp_total(xsh(ind_zone(m)+1))],'r',xsh(ind_zone(m)+1:ind_zone(m+1)),...
            txp_total(xsh(ind_zone(m)+1:ind_zone(m+1))),'r',repmat(xsh(ind_zone(m+1)),1,2),[profy(xsh(ind_zone(m+1))) txp_total(xsh(ind_zone(m+1)))],'r');
    end
end
xlabel('Distance(m)')
ylabel('Height (m)')

% Receiver locations.
xr=xt+1:profx(end); 
yr=profy(xr)+hr;

% Receiver locations that are and aren't in shadow areas.
% Direct visibility zones.
dr=sqrt(((xr-xt).^2)+((yr(xr-xt)-yt).^2));


% Free space loss
Lfs=32.5+20*log10(fc*1e-6)+20*log10(dr*1e-3);

hold off
figure,plot(xr,-Lfs,'k:')

% Plane earth loss
Lpe=120-20*log10(ht*hr)+40*log10(dr*1e-3);
hold on,plot(xr,-Lpe,'k-.')

% Diffraction loss due to the profile.
for m=1:length(xr)
    pos=[xp(find(xp<xr(m))) xr(m)];
    h=[height(find(xp<xr(m))) yr(m)];
    
    Ld=epstein_peterson(pos,h,length(pos)-2,lambdac);
    
    L(m)=max(Lfs(m),Lpe(m))+Ld;
end

plot(xr,-L,'k')
xlabel('Distance (m)')
ylabel('Inverse of path loss (dB)')
title('Total loss (dB)')
legend('Free space','Plane earth','Total');