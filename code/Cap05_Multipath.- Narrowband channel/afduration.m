function [axisafd,afd]=afduration(series,ts)
% returns afd with abscissas in dB/RMS and ordintes in s

SERIES=20*log10(series);  % SERIES is in dB and series is in linear units

afd=[];
durations=[];
currentduration=0;
axisafd=[];
for lev=ceil(min(SERIES))+1:floor(max(SERIES))-1
    afdaux=find(SERIES<=lev);
    countafd=0;
    for ii=2:length(afdaux),
        if afdaux(ii-1)+1 ~= afdaux(ii), 
            countafd=countafd+1; 
        end 
    end 
    if countafd ~= 0,
        afdaux2=length(afdaux)/countafd;
        afd=[afd afdaux2];
        axisafd=[axisafd lev];
    end  
end

afd=afd*ts;

    
 