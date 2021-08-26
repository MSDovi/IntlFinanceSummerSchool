warning off
close all; 
clear all; clc; 
format long 

%% USER DEFINED SETTINGS
%--------------------------------------------------------------------------
% Select whether to estimate the model
f_estim = 0;

% Select parameters to loop over ***will only work if f_estim=0***
models = {
          'mod1', 'eta_f', '0';...
          'mod2', 'eta_f', '0.5';...
          'mod3', 'eta_f', '1';...
          };

% Select shock, available shocks are listed below in the switch block
shock = {'eR_star'};
switch [shock{1}]
    case 'eR_star'
        f_shock = 1;
    case 'etheta_star'
        f_shock = 2;
    case 'eA_star'
        f_shock = 3;
    case 'etau_k'
        f_shock = 4;
    case 'etau_fk'
        f_shock = 5;
    case 'etau_bg'
        f_shock = 6;
end

% Select options for estimation 
if f_estim == 1
    addpath('IRFmatch Codes');
    % horizons to match IRF
    horizon_est  = 16; 
    % variables to match IRF
    var_est(1,1:2) = {'exp_spr_star_est', 'US EBP (bp)'};
    var_est(2,1:2) = {'pi_star_ann_est', 'US INFL (pp)'};
    var_est(3,1:2) = {'realGDP_star_est', 'US RGDP (%)'};
    var_est(4,1:2) = {'realGDP_est', 'RGDP (%)'};
    var_est(5,1:2) = {'pi_ann_est', 'INFL (pp)'};
    var_est(6,1:2) = {'realexp_est', 'EXP (%)'};
    var_est(7,1:2) = {'r_ann_est', 'IS (pp)'};
    var_est(8,1:2) = {'realexr_est', 'RXUSD (%)'};
    var_est(9,1:2) = {'eqp_est', 'EQUITY (%)'}; 
    % set priors
    vprior(1,1:4) = {'gamma_star', 'beta_pdf', 'gamma_star', '0.15'};
    vprior(2,1:4) = {'h', 'normal_pdf', 'h', '0.3'};
    vprior(3,1:4) = {'sigma', 'normal_pdf', 'sigma', '1'};
    vprior(4,1:4) = {'lev_star_bar', 'normal_pdf', 'lev_star_bar', '1'};
    vprior(5,1:4) = {'lev_bar', 'normal_pdf', 'lev_bar', '1'};
    vprior(6,1:4) = {'gamma_b', 'normal_pdf', 'gamma_b', '0.1'};
    vprior(7,1:4) = {'phik', 'normal_pdf', 'phik', '0.5'};
    vprior(8,1:4) = {'eta_f', 'beta_pdf', 'eta_f', '0.2'};
    vprior(9,1:4) = {'xi', 'beta_pdf', 'xi', '0.1'};
    vprior(10,1:4) = {'xi_star', 'beta_pdf', 'xi_star', '0.1'};
    vprior(11,1:4) = {'thetapi', 'normal_pdf', 'thetapi', '10'};
    vprior(12,1:4) = {'thetapi_star', 'normal_pdf', 'thetapi_star', '10'};

    % process the above settings
    nplots = size(var_est,1);
    fid_est_stoch = fopen('estim_stoch.mod','w');
    var_est_str = [' '];
    for i = 1 : nplots
        var_est_str=strcat(var_est_str, strcat(var_est{i,1}, {' '}));
    end
    fprintf(fid_est_stoch,'%s\n',['stoch_simul(order=1,irf=', num2str(horizon_est), ') ', var_est_str{1}, ';']);
    fclose(fid_est_stoch);
    fid_est_prior = fopen('estim_priors.mod','w');
    for i = 1 : size(vprior,1)
        fprintf(fid_est_prior,'%s\n',[vprior{i,1} ', ' vprior{i,2} ', ' vprior{i,3} ', ' vprior{i,4} ';']);
    end
    fclose(fid_est_prior);
end

% If estimation is selected over-write calibrated model comparison
if f_estim==1; models={}; end
nmodels = size(models,1);

% Generate macro-processor command for dynare file
fid_op = fopen('do_options.mod','w');
fprintf(fid_op,'%s\n',['@#define f_shock = ', num2str(f_shock)]);
fprintf(fid_op,'%s\n',['@#define f_estim=', num2str(f_estim)]);
fclose(fid_op);


%% SOLVE THE MODEL IN DYNARE
%--------------------------------------------------------------------------
if nmodels == 0
    fid = fopen('RECALIB_AQ.mod','w');
    fclose(fid);
    dynare CFR_Model_AQ noclearall
    close all
    % save IRF
    IR = oo_.irfs;
else
    for ii=1:nmodels
        % set parameters
        baseline_param_AQ; % give baseline parameterization...
        eval([models{ii,2} '=' models{ii,3} ';']) % ... but over-write the parameter you're looping over
        % solve for targets
        calibrate_model_AQ
        % save in CALIB
        fid = fopen('RECALIB_AQ.mod','w');
        % Parameter you're looping over:
        fprintf(fid,'%s\n',[models{ii,2} '=' models{ii,3} ';']);
        % Parameters used to match calibration targets:
        fprintf(fid,'%s\n',['gamma_b=', num2str(gamma_b,16), ';']);
        fprintf(fid,'%s\n',['theta=', num2str(theta,16), ';']);
        fprintf(fid,'%s\n',['theta_star_bar=', num2str(theta_star,16), ';']);
        fprintf(fid,'%s\n',['xib=', num2str(xib,16), ';']);
        fprintf(fid,'%s\n',['xib_star=', num2str(xib_star,16), ';']);
        fprintf(fid,'%s\n',['gdp_bar=', num2str(gdp_bar,16), ';']);
        fprintf(fid,'%s\n',['g_bar=', num2str(g_bar,16), ';']);
        fprintf(fid,'%s\n',['g_star_bar=', num2str(g_star_bar,16), ';']);
        fclose(fid);
        % run dynare
        dynare CFR_Model_AQ noclearall
        close all
        %save IRF
        IR.(models{ii,1}) = oo_.irfs;
    end
end


%% CHARTS AND RESULTS
%--------------------------------------------------------------------------
if f_estim == 0
    %% Select variables to table and plot (dynare name, plot name, multiplying factor)
    vselect(1,1:3)   = {'y','Output',100};
    vselect(2,1:3)   = {'c','Consumption',100};
    vselect(3,1:3)   = {'inv','Investment',100};
    vselect(4,1:3)   = {'l','Hours',100};
    vselect(5,1:3)   = {'pi_ann','CPI inflation',10000};
    vselect(6,1:3)   = {'r_ann','NIR',10000};
    vselect(7,1:3)   = {'ner','NER',100};
    vselect(8,1:3)   = {'rer','RER',100};
    vselect(9,1:3)   = {'q','Q',100};
    vselect(10,1:3)  = {'exp_spr_star','US Credit Spread',10000};
    vselect(11,1:3)  = {'k','Capital',100};
    vselect(12,1:3)  = {'pih_ann','PPI inflation',100};
    vselect(13,1:3)  = {'b_star_h','Interbank lending',100};
    vselect(14,1:3)  = {'yh_star','Export',100};
    vselect(15,1:3)  = {'yf','Import',100};

    %% Create Table of Moments
    varnames = char(vselect{:,2});
    for i = 1:size(varnames,1) 
        varidx(i)=strmatch(vselect{i,1},var_list_,'exact'); % get the position index of variables in the output file oo_
        std(i,:) = sqrt(oo_.var(varidx(i),varidx(i)));
    end
    aux = [oo_.mean(varidx) std];
    in.cnames = char('Mean','Std. Dev.');
    in.rnames = char('Variable', varnames);
    disp('-----------')
    mprint(aux,in)

    %% IRFs
    addpath(genpath('.../VARToolbox'))
    nsteps = 40;
    % Initialize plots
    nplots = size(vselect,1);
    FigSize(30,18); 
    row = 3; col = 6;
    IRF = nan(nsteps,length(vselect),nmodels);
    style = {'-',':','--','-.'};
    % Loop over variables
    for ii=1:nplots
        nameirf = {[vselect{ii,1} '_' shock{1}]};
        subplot(row,col,ii)
        % Loop over models
        if nmodels == 0
            h = plot(1:nsteps,vselect{ii,3}*IR.(nameirf{1})(1:nsteps),'LineWidth',2.5,'Color',cmap(1),'LineStyle',style{1}); hold on
        else
            for jj=1:nmodels %nmodels
               h(jj) = plot(1:nsteps,vselect{ii,3}*IR.(models{jj,1}).(nameirf{1})(1:nsteps),'LineWidth',2.5,'Color',cmap(jj),'LineStyle',style{jj}); hold on
            end
        end
        plot(1:nsteps,zeros(nsteps),'LineWidth',.5,'Color',rgb('black')); hold on
        set(gca,'xLim',[1 40],'xTick', [10 20 30 40]);
        if ii>col*(row-1) ; xlabel('Quarters'); end
        title(vselect{ii,2})
        if vselect{ii,3}==100; ylabel('Percent'); else; ylabel('Basis points'); end
        axis tight; box on
        set(gca,'Layer','top')
    end
    % Create legends
    if nmodels
        for kk = 1 : nmodels
          legendtag(kk) = {['\' models{kk,2} '=' models{kk,3}]};
        end
        legend(legendtag,'Location','SouthEast')
    end
%     saveFigure(['IRF_' shock{1}],1)
    clear vselect
    rmpath(genpath('.../VARToolbox'))

else % if f_estim == 1, Plot VAR IRFs and compare them with calibrated DSGE 
    
    FigSize   
    for jj = 1:nplots
        nameirf = {[var_est{jj,1} '_' shock{1}]};
        subplot(2,5,jj)
        H = PlotSwathe(IRFmatch_data.irf(:,jj),squeeze([IRFmatch_data.low(:,jj),IRFmatch_data.upp(:,jj)])); hold on; hg = plot([CALIB.oo_.irfs.(nameirf{1})(1:horizon_est)'],'--','LineWidth',1.5,'Color',cmap(2)); 
        title(var_est{jj,2}); legend([H.patch hg],{'VAR';'Calibrated DSGE'})
    end

    % Plot VAR IRFs and compare them with estimated DSGE
    FigSize
    for qq = 1:nplots
        nameirf = {[var_est{qq,1} '_' shock{1}]};
        subplot(2,5,qq)
        H = PlotSwathe(IRFmatch_data.irf(:,qq),squeeze([IRFmatch_data.low(:,qq),IRFmatch_data.upp(:,qq)])); hold on; hg = plot([ESTIM.oo_.irfs.(nameirf{1})(1:horizon_est)'],'--','LineWidth',1.5,'Color',cmap(2)); 
        title(var_est{qq,2}); legend([H.patch hg],{'VAR';'Estimated DSGE'})
    end

    % Print estimates and compare then with calibrated/priors
    calib_names = char(vprior{:,1});
    for i = 1:size(vprior,1) 
        calib_names_index(i)=strmatch(vprior{i,1},CALIB.M_.param_names,'exact'); % get the position index of variables in the output file oo_
    end
    calib_values = CALIB.M_.params(calib_names_index);
    aux = [calib_values xparam oo_.prior.mean];
    in.cnames = char('calibrated','posterior','prior');
    in.rnames = char(' ', calib_names);
    disp('-----------')
    disp('Estimation results:')
    mprint(aux,in)

end


