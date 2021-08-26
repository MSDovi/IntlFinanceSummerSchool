function SRout = SR(VAR,R,VARopt)
% =======================================================================
% Compute IRs, VDs, and HDs for a VAR model estimated with VARmodel and 
% identified with sign restrictions
% =======================================================================
% SRout = SR(VAR,R,VARopt)
% -----------------------------------------------------------------------
% INPUT
%   - VAR: structure, result of VARmodel function
%   - R: 3-D matrix containing the sign restrictions (nvar,3,nshocks)
%       described below
%   - VARopt: options of the VAR (see VARopt from VARmodel)
% -----------------------------------------------------------------------
% OUTPUT
%   - SRout
%       * IRall : 4-D matrix of IRs  (nsteps,nvar,nshocks,ndraws)
%       * VDall : 4-D matrix of VDs (nsteps,nvar,nshocks,ndraws)
%       * HDall : 4-D matrix of HDs (nsteps,nvar,nshocks,ndraws)
%       * Ball  : 4-D matrix of Bs (nvar,nvar,nshocks,ndraws)
%       * IRmed : median of IRall
%       * VDmed : median of VDall
%       * Bmed  : median of Ball
%       * HDmed : median of HDall
%       * IRinf : (100-pctg)/2 percentile of IRall
%       * VDinf : (100-pctg)/2 percentile of VDall
%       * IRinf : (100-pctg)/2 percentile of HDall
%       * IRsup : 100 - (100-pctg)/2 percentile of IRall
%       * VDsup : 100 - (100-pctg)/2 percentile of VDall
%       * HDsup : 100 - (100-pctg)/2 percentile of HDall
% =======================================================================
% Ambrogio Cesa Bianchi, March 2020
% ambrogio.cesabianchi@gmail.com

% Note. This code follows the notation as in the lecture notes available at
% https://sites.google.com/site/ambropo/MatlabCodes

% The R matrix is a 3-D matrix with dimension (nvar,3,nshocks). Assume you
% have 3 variables and you want to identify just one shock. Then:
% 
%              from        to          sign
%  R(:,:,1) = [ 1           4           1          % VAR1
%               1           4          -1          % VAR2
%               0           0           0];        % VAR3
% 
%   - The first column defines the first period from which the restriction is imposed 
%   - The second column defines the last period to which the restriction is imposed 
%   - The last column defines the sign of the restriction: positive (1) or negative (-1)
%   - To leave unrestricted set all three columns to zero
% 
% In the above example we set the following restrictions. VAR1 must respond
% with positive sign from period 1 to period 4; VAR2 must respond with 
% negative sign from period 1 to period 4; VAR3 is left unrestricted
%
% An additional shock could be defined as follows:
%              from        to          sign
%  R(:,:,2) = [ 1           4           1          % VAR1
%               1           4           1          % VAR2
%               1           4          -1];        % VAR3



%% Check inputs
%==========================================================================
if ~exist('R','var')
    error('You have not provided sign restrictions (R)')
end

if ~exist('VARopt','var')
    error('You need to provide VAR options (VARopt from VARmodel)');
end


%% Retrieve parameters and preallocate variables
%==========================================================================
nvar    = VAR.nvar;
nvar_ex = VAR.nvar_ex;
nsteps  = VARopt.nsteps;
ndraws  = VARopt.ndraws;
nobs    = VAR.nobs;
nlag    = VAR.nlag;
pctg    = VARopt.pctg;

% Determines the number of shocks, of sign restrictions (nsr) and zero 
% restrictions (nzr)
nshocks = size(R,3);
nzr = 0; 
for ss=1:nshocks
    if sum(sum(R(:,:,ss)))==0
        nzr = nzr +1;
    end    
end
nsr = nshocks - nzr;

% Initialize empty matrix for the IR draws
IRstore  = nan(nsteps,nvar,nvar,ndraws); 
VDstore = nan(nsteps,nvar,nvar,ndraws); 
HDstore_shock  = nan(nobs+nlag,nvar,nvar,ndraws); 
HDstore_init   = nan(nobs+nlag,nvar,ndraws); 
HDstore_const  = nan(nobs+nlag,nvar,ndraws); 
HDstore_trend  = nan(nobs+nlag,nvar,ndraws); 
HDstore_trend2 = nan(nobs+nlag,nvar,ndraws); 
HDstore_endo   = nan(nobs+nlag,nvar,ndraws); 
if nvar_ex>0; HDstore_exo = nan(nobs+nlag,nvar,nvar_ex,ndraws); end
Bstore = nan(nvar,nvar,ndraws); 


%% Sign restriction routine
%==========================================================================
jj = 0; % accepted draws
kk = 0; % total draws
ww = 1; % index for printing on screen
while jj < ndraws
    kk = kk+1;

    % Draw F and sigma from the posterior and set up VAR_draw
    [sigma_draw, Ft_draw] = VARdrawpost(VAR);
    VAR_draw = VAR;
    VAR_draw.Ft = Ft_draw;
    VAR_draw.sigma = sigma_draw;
    VARopt.ident = 'sr';
    
    % Reshape the R matrix to get the "SIGN" matrix (variables x shocks)
    % NOTE: does not support yet restrictions over multiple horizons
	SIGN = zeros(size(R,1)); % square matrix of size equal to num of variables
    for ss=1:size(R,3)  % for each shock, assing the signs to "SIGN"
        SIGN(:,ss) = R(:,3,ss);
    end

    % Compute rotated B matrix
    B = SignRestrictions(sigma_draw,[],SIGN); %B = B';
    VAR_draw.BfromSR = B; % Update VAR_draw with the rotated B matrix 
    jj = jj+1 ;

    % Compute and store IR, VD, HD
    [aux_irf, VAR_draw] = VARir(VAR_draw,VARopt); 
    IRstore(:,:,:,jj)  = aux_irf;
    aux_fevd = VARvd(VAR_draw,VARopt);
    VDstore(:,:,:,jj)  = aux_fevd;
    aux_hd = VARhd(VAR_draw);
    HDstore_shock(:,:,:,jj) = aux_hd.shock;
    HDstore_init(:,:,jj)    = aux_hd.init;
    HDstore_const(:,:,jj)   = aux_hd.const;
    HDstore_trend(:,:,jj)   = aux_hd.trend;
    HDstore_trend2(:,:,jj)  = aux_hd.trend2;
    HDstore_endo(:,:,jj)    = aux_hd.endo;
    if nvar_ex>0; HDstore_exo(:,:,jj) = aux_hd.exo; end

    % Store B
    Bstore(:,:,jj) = VAR_draw.B;
    
    % Display number of loops
    if jj==10*ww
        disp(['Loop: ' num2str(jj) ' / ' num2str(kk) ' draws']);
        ww=ww+1;
    end

end
disp('-- Done!');
disp(' ');


%% Store results
%==========================================================================
% Store all accepted IRs and VDs
SRout.IRall  = IRstore;
SRout.VDall = VDstore;
SRout.Ball = Bstore;

% Compute and save median impulse response
SRout.IRmed(:,:,:)  = median(IRstore,4);
SRout.VDmed(:,:,:) = median(VDstore,4);
SRout.Bmed(:,:,:) = median(Bstore,3);
SRout.HDmed.shock(:,:,:)  = median(HDstore_shock,4);
SRout.HDmed.init(:,:,:)   = median(HDstore_init,3);
SRout.HDmed.const(:,:,:)  = median(HDstore_const,3);
SRout.HDmed.trend(:,:,:)  = median(HDstore_trend,3);
SRout.HDmed.trend2(:,:,:) = median(HDstore_trend2,3);
SRout.HDmed.endo(:,:,:)   = median(HDstore_endo,3);
if nvar_ex>0; SRout.HDmed.exo(:,:,:) = median(HDstore_exo,3); end

% Compute lower and upper bounds
pctg_inf = (100-pctg)/2; 
pctg_sup = 100 - (100-pctg)/2;

aux = prctile(IRstore,[pctg_inf pctg_sup],4);
SRout.IRinf = aux(:,:,:,1);
SRout.IRsup = aux(:,:,:,2);

aux = prctile(VDstore,[pctg_inf pctg_sup],4);
SRout.VDinf = aux(:,:,:,1);
SRout.VDsup = aux(:,:,:,2);

aux = prctile(HDstore_shock,[pctg_inf pctg_sup],4);
SRout.HDinf.shock = aux(:,:,:,1);
SRout.HDsup.shock = aux(:,:,:,2);

aux = prctile(HDstore_init,[pctg_inf pctg_sup],3);
SRout.HDinf.init = aux(:,:,1);
SRout.HDsup.init = aux(:,:,2);

aux = prctile(HDstore_const,[pctg_inf pctg_sup],3);
SRout.HDinf.const = aux(:,:,1);
SRout.HDsup.const = aux(:,:,2);

aux = prctile(HDstore_trend,[pctg_inf pctg_sup],3);
SRout.HDinf.trend = aux(:,:,1);
SRout.HDsup.trend = aux(:,:,2);

aux = prctile(HDstore_trend2,[pctg_inf pctg_sup],3);
SRout.HDinf.trend2 = aux(:,:,1);
SRout.HDsup.trend2 = aux(:,:,2);

aux = prctile(HDstore_endo,[pctg_inf pctg_sup],3);
SRout.HDinf.endo = aux(:,:,1);
SRout.HDsup.endo = aux(:,:,2);

if nvar_ex>0
    aux = prctile(HDstore_exo,[pctg_inf pctg_sup],3);
    SRout.HDinf.exo = aux(:,:,1);
    SRout.HDsup.exo = aux(:,:,2);
end

clear aux
