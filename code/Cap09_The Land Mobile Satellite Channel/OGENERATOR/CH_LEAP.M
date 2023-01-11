
function [leap] = ch_leap (year)
%
% Checks if the year is a leap-year
%
% Inputs:
%	year	: year 
%
% Outputs:
% 	leap	: 1 if leap year, 0 otherwise
%
%
% By: Manuel Bentosinos-Rico - University of Vigo, Spain
% Release: 6/98



leap = 0;
if (round (year/100) == (year/100))
    if (round (year/400) == (year/400))
        leap = 1;
    end
elseif (round (year/4) == (year/4))
	leap = 1;
end
