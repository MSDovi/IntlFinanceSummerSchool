
    @#include "estim_stoch.mod"    

    // Save stuff for later
    CALIB.M_ = M_;
    CALIB.oo_ = oo_;
    save('CALIB','CALIB')

    // IRF matching settings
    %--------------------------------------------------------------------------
    // Priors (center priors on true VAR parameters) ***
    estimated_params;
    @#include "estim_priors.mod"
    end;
             
    // Declare observed variable
    varobs r_star_ann;

    // Settings
    options_.irfs_match_estimation = 1;                                        % IRF Matching

    // Settings
    load IRFmatch_data;
     % Action: Back up the standard deviation from IRFmatch_data here

    for i = 1:size(StdNorm,1)                                                  % Replace any 0 standard deviations in StdNorm with an arbitrary value
        for j = 1:size(StdNorm,2)
            if StdNorm(i,j) == 0
             StdNorm(i,j) = 0.001;
            end
        end
    end

    % Define Var_IRF, invVarNorm, and logdetVarNorm here
    
    M_.Var_IRF         = Var_IRF;                                              % set model settings to align with data
    M_.invVarNorm      = invVarNorm;
    M_.logdetVarNorm   = logdetVarNorm;
    M_.horizon_est     = horizon_est;
    M_.invVarNorm      = invVarNorm;
    M_.mod_var_list    = {var_est{:,1}};
    M_.mod_shock_list  = shock;
    M_.SR              = [];       


    // Estimation
    %--------------------------------------------------------------------------
    estimation(datafile='r_star_ann.mat',diffuse_filter,mh_replic=0,mh_nblocks=2,mh_jscale=0.8,mode_compute = 7,plot_priors = 1);
    xparam = get_posterior_parameters('mode',M_,estim_params_,oo_,options_);
    M_ = set_all_parameters(xparam,estim_params_,M_);
    @#include "estim_stoch.mod"


    ESTIM.M_ = M_;
    ESTIM.oo_ = oo_;
    save('ESTIM','ESTIM');