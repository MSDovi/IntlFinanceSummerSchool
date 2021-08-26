clear;
%% Setup
addpath('ProblemSet') %change this to 'Solutions' to compare your results to the solutions
country = 'US'; %US, Canada, Sweden, UK
iv_type = 'FF4'; %FF4 (three-month futures), FF1 (current month futures)
p = 12; %number of lags
horizon = 25; %number of horizons for IRF
scaling = 0.1924518; %scaling in GK--can be set to anything
boot_reps = 2500;
sig_level = 0.1;
 
%% Load data and construct matrices
data_raw = readtable(char(append('Data/', country, '_', iv_type, '.xlsx')), 'VariableNamingRule', 'preserve');
Y = table2array(data_raw(:, 1:end-1));
n = size(Y, 2); 
Z = table2array(data_raw(:, end));
 
%% Calculate reduced-form VAR
rs = reduced_form(Y, p);
 
%% Identify column of structural covariance matrix
rs = ivsvar(rs, Z, horizon, scaling);
Z = rs.Z; %shortened IV

%% Calculate confidence bands based on Wild bootstrap
T = rs.T;
T_Z = rs.T_Z;

irf_boot = zeros(horizon, n, boot_reps);
for rep = 1:boot_reps
    rand_index = ((rand(1, T)<.5)*2 - 1)';
    res_boot = rs.res .* (rand_index * ones(1, n));
    Y_star = zeros(p + T, n);
    Y_star(1:p, :) = rs.Y(1:p, :);
    for i = p + 1 : p + T
        Y_lagged = Y_star(i-p:i-1, :);
        Y_lagged = flip(Y_lagged)';
        Y_lagged = Y_lagged(:);
        Y_star(i, :) = Y_lagged'*rs.rf_coeff(2:(p*n)+1, :) + rs.rf_coeff(1, :) +  res_boot(i - p, :);
    end
    Z_star = Z.*rand_index(T - T_Z + 1:end);
    boot_rs = reduced_form(Y_star, p);
    boot_rs = ivsvar(boot_rs, Z_star, horizon, scaling);
    irf_boot(:, :, rep) = boot_rs.irf;
end
        
%% Create plots
irf = rs.irf;
irf_lower = quantile(irf_boot, sig_level, 3);
irf_upper = quantile(irf_boot, 1 - sig_level, 3);

for i=1:n
    subplot(3, 2,  i);
    hold on
    plot(1:1:horizon, irf(:, i), '-k');
    plot(1:1:horizon, irf_lower(:, i), '--k');
    plot(1:1:horizon, irf_upper(:, i), '--k');
    title(data_raw.Properties.VariableNames(i));
    yline(0);
end
