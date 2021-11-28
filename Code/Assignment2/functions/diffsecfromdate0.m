function [t] = diffsecfromdate0(date,date0)


% diffsecfromdate0.m - computes the difference in seconds between two dates.
%
% PROTOTYPE:
%   [t] = diffsecfromdate0(date,date0) 
%
% INPUT:
%   date  [6]     Second date as 6-element vector of numbers
%                 [year, month, day, hour, minute, second].
%
%   date0 [6]     First date as 6-element vector of numbers
%                 [year, month, day, hour, minute, second].
%
% OUTPUT:
%   t     [1x1]   Time value of dates' gap                          [s]
%    
% CONTRIBUTORS:
%   Andrea Bersani
%   Giovanni Chiarolla
%   Jacopo Fabbri
%   Matteo Manicaglia
%
% VERSIONS:
%   2021-1: Last version

year0 = date0(1);
month0 = date0(2);
day0 = date0(3);
hour0 = date0(4);
minute0 = date0(5);
second0 = date0(6);

year = date(1);
month = date(2);
day = date(3);
hour = date(4);
minute = date(5);
second = date(6);

if year ~= year0
    ydiff = year-year0;
    daysYvect = zeros(1,ydiff);
    for kY = 1:ydiff
        yearC = year0+kY-1;
        if mod(yearC,4) ~= 0 || (mod(yearC,100) == 0 && mod(yearC,400) ~= 0)
            daysY = 365;
        else
            daysY = 366;
        end
        daysYvect(kY) = daysY;
    end
else
    daysYvect = 0;
end

if mod(year,4) ~= 0 || (mod(year,100) == 0 && mod(year,400) ~= 0)
    nfeb = 28;
else
    nfeb = 29;
end

daysMvect = [31,nfeb,31,30,31,30,31,31,30,31,30,31];
for kM = 1:12
    if month == kM
        daysM = sum(daysMvect(1:kM-1));
    end
end
for kM0 = 1:12
    if month0 == kM0
        daysM0 = sum(daysMvect(1:kM0-1));
    end
end
daysMdiff = daysM-daysM0;

ddiff = day - day0;
hdiff = hour - hour0;
mdiff = minute - minute0;
sdiff = second - second0;

t = (((sum(daysYvect(:))+daysMdiff+ddiff)*24+hdiff)*60+mdiff)*60+sdiff;
end