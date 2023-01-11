%
% project92
%
% =======================================================================
clear
close all
clc
% ========================================================================
fMHz=2000;    % frequency in MHz

BldHeight=10;   % building heigths in m 
StWidth=15;   % Street Width in m
Dy=3;         % Separation from center of street in m
hR=1.5;       % MS antenna height in m
% ======================================================================

MKA=atand(BldHeight/(StWidth-Dy));
message='Masking angle : '
disp([message num2str(MKA) ' degree'])

% =======================================================================
load newSatMatrix

% =======================================================================

[a,b]=size(newSatMatrix);
NoSats=b/3;


SatConnections=zeros(a,NoSats);

for ii=1:NoSats     % go through all satellites
    for jj=1:a                    % go through all time instants
        if newSatMatrix(jj,1+3*(ii-1))> 0,     % coverage is possible only for elevations > 0 deg

            if newSatMatrix(jj,2+3*(ii-1))>= 0 & newSatMatrix(jj,2+3*(ii-1))<= 90
                % 1st quadrant
                elev=newSatMatrix(jj,1+3*(ii-1));
                azim=newSatMatrix(jj,2+3*(ii-1));
                hSat=(StWidth-Dy)/cosd(azim)*tand(elev)+hR;
                if hSat > BldHeight,
                    SatConnections(jj,ii)=1;
                else
                    SatConnections(jj,ii)=0;
                end
            end                        % end of 1st quadrant

            if newSatMatrix(jj,2+3*(ii-1))> 90 & newSatMatrix(jj,2+3*(ii-1))<= 180
                % 2nd quadrant
                elev=newSatMatrix(jj,1+3*(ii-1));
                azim=180-newSatMatrix(jj,2+3*(ii-1));
                hSat=(StWidth+Dy)/cosd(azim)*tand(elev)+hR;
                if hSat > BldHeight,
                    SatConnections(jj,ii)=1;
                else
                    SatConnections(jj,ii)=0;
                end
            end                        % end of 2nd quadrant

            if newSatMatrix(jj,2+3*(ii-1))> 180 & newSatMatrix(jj,2+3*(ii-1))< 270
                % 3rd quadrant
                elev=newSatMatrix(jj,1+3*(ii-1));
                azim=newSatMatrix(jj,2+3*(ii-1))-180;
                hSat=(StWidth+Dy)/cosd(azim)*tand(elev)+hR;
                if hSat > BldHeight,
                    SatConnections(jj,ii)=1;
                else
                    SatConnections(jj,ii)=0;
                end
            end                       % end of 3rd quadrant

            if newSatMatrix(jj,2+3*(ii-1))>= 270 & newSatMatrix(jj,2+3*(ii-1))< 360
                % 4th quadrant
                elev=newSatMatrix(jj,1+3*(ii-1));
                azim=360-newSatMatrix(jj,2+3*(ii-1));
                hSat=(StWidth-Dy)/cosd(azim)*tand(elev)+hR;
                if hSat > BldHeight,
                    SatConnections(jj,ii)=1;
                else
                    SatConnections(jj,ii)=0;
                end
            end                       % end of 4th quadrant

        else
            SatConnections(jj,ii)=0;          % no coverage elev < 0 degreee
        end
        %         pause
    end     % end of jj for (time)
end         % end of ii for (sats)

for ii=1:NoSats
    figure,plot(SatConnections(:,ii)*90,'k')
    hold
    plot(newSatMatrix(:,1+3*(ii-1)),'k:')
    title(['Satellite # ' num2str(ii)])
    xlabel('time (s)')
    ylabel('Sat. elevation (degree) / ON-OFF coverage')
    axis([1 a -0.5 100])
    figure,plot(newSatMatrix(:,2+3*(ii-1)),'k:')
    hold
    plot(SatConnections(:,ii)*90,'k')
    title(['Satellite # ' num2str(ii)])
    xlabel('time (s)')
    ylabel('Sat. azimuth (degree)')
end

figure, hold on
for ii=1:NoSats
    plot(SatConnections(:,ii))
    axis([1 a -0.5 1.5])
end
xlabel('time (s)')
ylabel('ON-OFF coverage')

figure,pcolor(SatConnections')
ylabel('satellite #')
xlabel('time (s)')
title('Satellite availabilities v. time')
shading flat

figure
for ii=1:NoSats
    out=mmpolar(newSatMatrix(:,2+3*(ii-1))*pi/180,-newSatMatrix(:,1+3*(ii-1)),'k',...
        'Style','compass','tticksign','+','ttickdelta',10,'RTickAngle',5, 'Rlimit',[-90 0] );
    P.RtickLabel={'80','60','40','20','0'};
    mmpolar(P)
    hold on
end
text(0.9,-0.05,'Elevation','Rotation',270)
text(0.15,1.2,'Azimuth','Rotation',270)

