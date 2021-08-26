function opt = FigFontOption(fsize)
% =======================================================================
% Inputs for FigFont
% =======================================================================
% OPTIONAL INPUT
%   - fsize: size of the font [default 12]
% =========================================================================
% Ambrogio Cesa Bianchi, March 2015
% ambrogio.cesabianchi@gmail.com


%% CHECK INPUT
% =======================================================================
if ~exist('fsize','var')
    fsize = 12;
end

%% SET OPTIONS
% =======================================================================
opt.fsize   = fsize;
opt.fname   = 'Helvetica';
opt.fweight = 'light';
