function [jd, day] = jul_date (year,month, day, hours, mins, seconds)
%
% Calculation of Julian date corresponding to (h,m,s) UT 
% on month day, year
% The format of the year is four digits
%
% Inputs:
% 	year, month, day, hour, mins, seconds: Start time (format: 1997, 30, 1, 14, 26, 34)
% Outputd:
%	jd	: julian date corresponding to the given day
%	day	: universal time day
%
%
% By: Manuel Bentosinos-Rico - University of Vigo, Spain
% Release: 6/98




% The reference year is 1986
x = year - 1986;       

jd = x * 365;
[sal]=ch_leap(year-1);
if (sal == 1)
    jd = jd + floor(x/4) + 1;   % A day per leap-year is added
else
    jd = jd + floor(x/4);
end
ndays_feb = 28;   % Number of days of February (non-leap year)
[sal]=ch_leap(year-1);
if (sal == 1)
    ndays_feb = 29;
end

if (month == 2)
    day = day + 31;
elseif (month == 3)
    day = day + ndays_feb + 31;
elseif (month == 4)
    day = day + ndays_feb + 62;
elseif (month == 5)
    day = day + ndays_feb + 92;
elseif (month == 6)
    day = day + ndays_feb + 123;
elseif (month == 7)
    day = day + ndays_feb + 153;
elseif (month == 8)
    day = day + ndays_feb + 184;
elseif (month == 9)
    day = day + ndays_feb + 215;
elseif (month == 10)
    day = day + ndays_feb + 245;
elseif (month == 11)
    day = day + ndays_feb + 276; 
elseif (month == 12)
    day = day + ndays_feb + 306; 
end

% Calculation of UT
univ_time = (hours + mins/60 + seconds/3600) / 24;
day = day + univ_time;
jd = jd + day;

% Adding to the Julian date for Jan 0.0, 1986
jd = jd + 2446430.5;