%
% project94
%
% =======================================================================
clear
close all
clc
% ========================================================================
fMHz=2000;    % frequency in MHz
lambdac=300/fMHz;   % wavelength in m
ts=1;    % dampling interval in s
% =======================================================================
load newSatMatrix

% =======================================================================

[a,b]=size(newSatMatrix);
NoSats=b/3;

SatPhases=zeros(a,NoSats);
for ii=1:NoSats                   % go through all satellites
    for jj=1:a                    % go through all time instants
        if newSatMatrix(jj,1+3*(ii-1))> 0,     % coverage is possible only for elevations > 0 deg
            SatPhases(jj,ii)=-2*pi*newSatMatrix(jj,3+3*(ii-1))*1000/lambdac;
        else
           SatPhases(jj,ii)=nan;          % no coverage elev < 0 degreee 
        end
    end     % end of jj for (time)
end         % end of ii for (sats)

for ii=1:NoSats
    figure,plot(SatPhases(:,ii),'k')
    xlabel('Time (s)')
    ylabel('Phase (Rad)')
    title(['Satellite # ' num2str(ii)])
end

SatDopplers=1/(2*pi)*diff(SatPhases);

for ii=1:NoSats
    figure,plot(SatDopplers(:,ii),'k')
    xlabel('Time (s)')
    ylabel('Doppler shift (Hz)')
    title(['Satellite # ' num2str(ii)])
end