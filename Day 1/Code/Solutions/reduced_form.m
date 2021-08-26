function rs = reduced_form(Y, p)

%% Compute reduced-form VAR
X = [lagmatrix(Y, [1:1:p])];
X = X(p + 1: end, :);
Y = Y(p + 1: end, :);
T = size(Y, 1);
X = [ones(T, 1), X];
rf_coeff = (X'*X)\(X'*Y);
rf_res = Y - X*rf_coeff;

rs.Y = Y;
rs.X = X;
rs.rf_coeff = rf_coeff;
rs.res = rf_res;
[rs.T rs.n] = size(Y);
rs.p = p;
rs.T = T;



end

