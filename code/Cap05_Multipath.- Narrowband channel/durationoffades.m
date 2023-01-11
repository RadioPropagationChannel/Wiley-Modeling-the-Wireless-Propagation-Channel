function [durationsatlevel]=durationoffades(series,ts,lev1);

% returns vector of all durations below specified level

SERIES=20*log10(series);

durationsatlevel=[];
duraux=find(SERIES<=lev1);
currentdur=1;
    for ii=2:length(duraux),
        if duraux(ii-1)+1 == duraux(ii),
            currentdur=currentdur+1; 
        else
            durationsatlevel=[durationsatlevel; currentdur];
            currentdur=1;
        end     
    end 

durationsatlevel=durationsatlevel*ts;



