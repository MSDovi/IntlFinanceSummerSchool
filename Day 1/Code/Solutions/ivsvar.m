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
res_reduced = rf_res(T - T_Z + 1: T, :);
pi_hat = (Z'*Z)\(Z'*res_reduced(:, 1));
res_fit = [ones(T_Z, 1) Z*pi_hat];
tsls_coeff = (res_fit'*res_fit)\res_fit'*res_reduced(:, 2:end);
tsls_coeff = scaling * tsls_coeff(2:end, :); %remove constant coefficient
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

