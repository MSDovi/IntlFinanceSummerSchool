function OUT = Date2Num(dates,frequency)
% =======================================================================
% Convert a cell array of dates into a double array. The first period of 
% the year is coded as a round number. For quarterly data the step change 
% is 0.25, so that 1984Q2 => 1984.25; For  monthly data the step change is 
% 0.833, so that 1984M2 => 1984.0833.
% =======================================================================
% OUT = Date2Num(DATA)
% -----------------------------------------------------------------------
% INPUTS 
%	- DATA: a (T x 1) vector with dates in cell format 
% -----------------------------------------------------------------------
% OPTIONAL INPUT
%   - frequency : monthly ('m'), quarterly ('q') [default] or annual ('y') 
% -----------------------------------------------------------------------
% OUTPUT
%	- OUT:  a (T x 1) vector with dates in numeric format
% =======================================================================
% EXAMPLE 
%   - Convert to double a monthly cell array:
%       dates = {'1992M1';'1992M2';'1992M3';'1992M4';'1992M5';};
%       OUT = Date2Num(dates,'m')
% RELATED
%   - DatesCreate, DatesCount
% =======================================================================
% Ambrogio Cesa Bianchi, December 2016
% ambrogiocesabianchi@gmail.com
% -----------------------------------------------------------------------

% Check inputs
if ~exist('frequency','var')
    frequency = 'q';
end
[n,c] = size(dates);

% Convert cell array into double array
if strcmp(frequency,'q')
    OUT = nan(n,c);
    for i=1:n
        for j=1:c
            aux = dates{i,j};
            year = str2double(aux(1:4));
            period = str2double(aux(6));
            OUT(i,j) = year + period/4 -0.25;
        end
    end
elseif strcmp(frequency,'m')
    OUT = nan(n,c);
    for i=1:n
        for j=1:c
            aux = dates{i,j};
            year = str2double(aux(1:4));
            period = str2double(aux(6));
            OUT(i,j) = year + period/12 -0.0833;
        end
    end
end