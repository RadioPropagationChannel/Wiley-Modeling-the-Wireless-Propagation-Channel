function [axislcr,lcr]=lcrate(series,ts)
SERIES=20*log10(series);  % SERIES is in dB and series is in linear units

crossings=[];
axislcr=[ceil(min(SERIES))+1:floor(max(SERIES))-1];

for lev=ceil(min(SERIES))+1:floor(max(SERIES))-1
    SERIESaux=SERIES-lev; 
    countlcr=0;
    for ii=2:length(SERIES)
        if SERIESaux(ii)>0 & SERIESaux(ii-1)<=0,
            countlcr=countlcr+1;
        end
    end
    crossings=[crossings countlcr];
end
durationobservation=length(SERIES)*ts;    
lcr=crossings/durationobservation;