%% Calculate the Steady State

l=log(l_bar);
rk=log(rk_bar);
q=log(1);
r=log(r_bar);
pi=0;
pm=log((rho-1)/rho);
ph=log(1);
pf=log(1);
oph_pcp=log(1);
oph_pcp_star=log(1);
oph_lcp=log(1);
oph_lcp_star=log(1);
dh=log(1);
dh_star=log(1);
A=log(1);
pih=0;
pif=0;
lev=log(lev_bar);
x=log(x_bar);

%l_star=log(l_star_bar);
rk_star=log(rk_star_bar);
q_star=log(1);
r_star=log(r_star_bar);
rb_star=r_star;
pi_star=0;
pm_star=log((rho_star-1)/rho_star);
ph_star=log(1);
pf_star=log(1);
opf_star=log(1);
df_star=log(1);
A_star=log(1);
pih_star=0;
pif_star=0;
lev_star=log(lev_star_bar);
theta_star=theta_star_bar;

tau_k=0;
tau_bg=0;

rer=log(1);
ner=log(1);

rr=r;
w=log((exp(A)*exp(pm)*(alfa^alfa)*((1-alfa)^(1-alfa))/((exp(rk)-1+delta)^alfa))^(1/(1-alfa)));
ow=w;
k=log((alfa/(1-alfa))*exp(w)*exp(l)/(exp(rk)-1+delta));
inv=log(delta*exp(k));
y=log(exp(w)*exp(l)/(1-alfa)/exp(pm));
%l=log((1-alfa)*exp(pm)*exp(y)/exp(w));
z=k;
b_star=log(exp(x)*exp(z));

rb_star=rk_star;
rr_star=r_star; 
w_star=log((exp(A_star)*exp(pm_star)*(alfa_star^alfa_star)*((1-alfa_star)^(1-alfa_star))/((exp(rk_star)-1+delta_star)^alfa_star))^(1/(1-alfa_star)));
y_star=log(((1-a)/(a-a_star)*exp(y)-(1-exp(rb_star))*exp(b_star))/(a_star/(a-a_star)*(1-n)/n));
l_star=log(exp(y_star)*(1-alfa_star)*exp(pm_star)/exp(w_star));
ow_star=w_star;
k_star=log((alfa_star/(1-alfa_star))*exp(w_star)*exp(l_star)/(exp(rk_star)-1+delta_star));
inv_star=log(delta_star*exp(k_star));
z_star=k_star;

psi_star=log(exp(lev_star)*theta_star);
mu0_star=log(beta_star*(1-omega_star+omega_star*exp(psi_star))*(exp(rk_star)-exp(r_star)));
mu1_star=log(beta_star*(1-omega_star+omega_star*exp(psi_star))*exp(r_star));
kappa=log(theta*(1+gamma_b/2*(exp(x)^2))*exp(lev));
muk=log(beta*(1-omega+omega*exp(kappa))*(exp(rk)-exp(r)*(1+tau_k)));
mud=log(beta*(1-omega+omega*exp(kappa))*exp(r));
mub=log(beta*(1-omega+omega*exp(kappa))*(exp(r)*(1+tau_k)-exp(rb_star)*(1+tau_k)/(1-tau_bg)));
mu=log(exp(mub)/exp(muk));
nw=log(1/exp(lev)*exp(z));
d=log((1+tau_k)*exp(z)-exp(nw)-exp(b_star)*(1-tau_bg));
d_star=log((1-1/exp(lev_star))*(exp(z_star)+exp(b_star)));
nw_star=log(1/exp(lev_star)*(exp(z_star)+exp(b_star)));

yd=log(((1-a_star)*exp(y)-a_star*(1-n)/n*exp(y_star))/(a-a_star));
yd_star=log((a*exp(y_star)-(1-a)*n/(1-n)*exp(y))/(a-a_star));
yh=log(a*exp(yd));
yf=log((1-a)*(n/(1-n))*exp(yd));
yh_star=log(a_star*((1-n)/n)*exp(yd_star));
yf_star=log((1-a_star)*exp(yd_star));
lambda=log(chi*(exp(l)^zeta)*nul/(exp(w)*(nul-1)));
lambda_star=log(chi_star*(exp(l_star)^zeta_star)*nul_star/(exp(w_star)*(nul_star-1)));
c=log((exp(lambda)^(-1/sigma))/(1-h));
c_star=log((exp(lambda_star)^(-1/sigma_star))/(1-h_star));
xw1=log((exp(l)^(1+zeta))*(exp(w)^(nul*zeta))/(1-beta*xiw));
xw2=log(exp(lambda)*exp(l)/(1-beta*xiw));
x1_pcp=log(exp(pm)*exp(y)*exp(lambda)/(1-beta*xi));
x2_pcp=log(exp(y)*exp(lambda)/(1-beta*xi));
x1_lcp=log(exp(pm)*exp(yh)*exp(lambda)/(1-beta*xi));
x2_lcp=log(exp(yh)*exp(lambda)/(1-beta*xi));
x3_lcp=log(exp(pm)*exp(yh_star)*exp(lambda)/(1-beta*xi));
x4_lcp=log(exp(yh_star)*exp(lambda)/(1-beta*xi));
xw1_star=log((exp(l_star)^(1+zeta_star))*(exp(w_star)^(nul_star*zeta_star))/(1-beta_star*xiw_star));
xw2_star=log(exp(lambda_star)*exp(l_star)/(1-beta_star*xiw_star));
x3_star=log(exp(pm_star)*exp(y_star)*exp(lambda_star)/(1-beta_star*xi_star));
x4_star=log(exp(y_star)*exp(lambda_star)/(1-beta_star*xi_star));
g=log(exp(yd)-exp(c)-exp(inv));
g_star=log(exp(yd_star)-exp(c_star)-exp(inv_star));
Deltah=log(1);
delta1=log(1);
delta2=log(1);
tot=log(1);

%% Calculate Parameters from Targets

gamma_b=2*exp(mu)/((exp(x)^2)*exp(mu)+2*exp(x))

theta1=(beta*(1-omega)*exp(r)+exp(lev)*beta*(1-omega)*(exp(rk)-exp(r))+exp(lev)*beta*exp(x)*(1-omega)*(exp(r)-exp(rb_star)));
theta2=exp(lev)*((gamma_b/2)*(exp(x)^2)+1)*(1-beta*omega*exp(r)-exp(lev)*beta*omega*(exp(rk)-exp(r))-exp(x)*exp(lev)*beta*omega*(exp(r)-exp(rb_star)));
theta=theta1/theta2

% thetah=beta*(1-omega)*((exp(rkh)-exp(r))*exp(levh)+exp(r))/(exp(levh)-beta*omega*(exp(levh)^2)*(exp(rkh)-exp(r))-exp(levh)*beta*omega*exp(r))
% thetaf=beta*(1-omega)*((exp(rkf)-exp(rb_star))*exp(levf)+exp(rb_star))/(exp(levf)-beta*omega*(exp(levf)^2)*(exp(rkf)-exp(rb_star))-exp(levf)*beta*omega*exp(rb_star))
theta_star=beta_star*(1-omega_star)*((exp(rk_star)-exp(r_star))*exp(lev_star)+exp(r_star))/(exp(lev_star)-beta_star*omega_star*(exp(lev_star)^2)*(exp(rk_star)-exp(r_star))-exp(lev_star)*beta_star*omega_star*exp(r_star))

xib=(exp(nw)+omega*exp(r)*exp(d)+omega*exp(rb_star)*exp(b_star))/(exp(rk)*exp(k))-omega
xib_star=(exp(nw_star)+omega_star*exp(r_star)*exp(d_star))/(exp(rk_star)*exp(k_star)+exp(rb_star)*exp(b_star))-omega_star

% chi=exp(w)*(nul-1)*exp(lambda)/((exp(l)^zeta)*nul)
% chi_star=exp(w_star)*(nul_star-1)*exp(lambda_star)/((exp(l_star)^zeta)*nul_star)

gdp_bar=exp(y)
g_bar=exp(g)
g_star_bar=exp(g_star)

% clearvars -except thetah thetaf theta_star xibh xibf xib_star g_bar g_star_bar gdp_bar