function rs = reduced_form(Y, p)

%% Compute reduced-form VAR
% To Do:
% 1. Compute RHS lagged VAR matrix                ==> save under X
% 2. Remove NAs from X and Y                      ==> overwrite X and Y
% 2. Compute reduced-form coefficients of the VAR ==> save under rf_coeff
% 3. Compute number of effective observations     ==> save under T
% 4. Compute reduced-form residuals               ==> save under rf_res






rs.Y = Y;
rs.X = X;
rs.rf_coeff = rf_coeff;
rs.res = rf_res;
[rs.T rs.n] = size(Y);
rs.p = p;
rs.T = T;


end

