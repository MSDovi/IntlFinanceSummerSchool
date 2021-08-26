function rs = ivsvar(rs, Z, horizon, scaling)

%% Unpack struct and save other useful quantities
rf_res = rs.res;
T = rs.T;
T_Z_start = sum(isnan(Z)) + 1;
Z = Z(T_Z_start: end);
T_Z = size(Z, 1);
Z = [ones(T_Z, 1), Z]; 
p = rs.p;
n = rs.n;
rf_coeff = rs.rf_coeff;


%% Identification of structural covariance matrix (up to scale)
% To Do:
% 1. Extract the residuals corresponding to the period for which IV is
%    available (use code in section above)
% 2. Compute TSLS coefficient for n = 2, ..., 6 ==> save under tsls_coeff


irf_coeffs = [scaling tsls_coeff];

%% Compute IRFs
irf = zeros(p + horizon, n);
irf(p + 1, :) = irf_coeffs;
for i = 2:horizon
    lv = irf(i:p + i - 1, :);
    lv = flip(lv)';
    lv = lv(:);
    irf(p + i, :) = lv'*rf_coeff(2:n*p+1, :);
end
irf = irf(p + 1: p + horizon, :);

rs.Z = Z(:, 2); %remove the constant
rs.irf = irf;
rs.tsls_coeff = tsls_coeff;
rs.T_Z = T_Z;


end

