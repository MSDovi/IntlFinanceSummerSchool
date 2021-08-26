% Calibration

% home parameters
h=0.65;
sigma=1;
chi=24;
zeta=1;
beta=1.03^(-0.25);
openness=0.4;
n=0.1;
a=1-(1-n)*openness;
epsilon=1.5;
rho=6;
nul=21;
alfa=0.36;
phik=1.5;
delta=0.02;
xi=0.6;
xiw=0.65;
rhoA=0.946;
rhott=0.75;
eta=3/4;
eta_f=0.5;
xib=6.588920742448545e-04;
theta=0.547347639854640;
gamma_b=1.065535593466435;
omega=0.97;
gamma=0.75;
thetapi=1.5;
thetaxr=0;
rhog=0;
varphi_k=0;
varphi_bg=0;

% foreign parameters
h_star=0.65;
sigma_star=1;
chi_star=24;
zeta_star=1;
beta_star=1.01^(-0.25);
a_star=n*openness;
epsilon_star=1.5;
rho_star=6;
nul_star=21;
alfa_star=0.36;
phik_star=1.5;
delta_star=0.02;
xi_star=0.6;
xiw_star=0.65;
rhoA_star=0.946;
rhott_star=0.75;
xib_star=0.001914433567820;
omega_star=0.97;
gamma_star=0.75;
thetapi_star=1.5;
rhog_star=0;

% calibration targets
l_bar=0.3;
rk_bar=1.045^0.25;
rk_star_bar=1.025^0.25;
lev_bar=5;
x_bar=0.3;
lev_star_bar=5;

% steady state values
r_bar=1/beta;
r_star_bar=1/beta_star;
theta_star_bar=0.505356397803559;  % theta_star becomes a variable subject to exogenous process
g_bar=0.125731678202725;
g_star_bar=0.005969010980998;
gdp_bar=1.074235242915810;