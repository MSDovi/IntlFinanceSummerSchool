var

// variables of Home country

c // consumption
l // hours worked
w // real wage
ow // optimal wage setting
xw1 // optimal wage setting factor
xw2 // optimal wage setting factor
lambda // marginal utility of consumption
k // capital
rk // rental rate on capital
q // price for capital
inv // investment
r // gross nominal deposit rate
pi // CPI inflation rate
yd // total demand
yh // home goods demand
yf // foreign goods demand
y // output
pm // price of the intermediate goods
ph // relative price of the home good bundle over the CPI index
pf // relative price of the foreign good bundle over the CPI index
oph_pcp // PCP optimal pricing of home goods in home over the home PPI
oph_pcp_star // PCP optimal pricing of home goods in foreign over the foreigh PPI
oph_lcp // LCP optimal pricing of home goods in home over the home PPI
oph_lcp_star // LCP optimal pricing of home goods in foreign over the foreigh PPI
x1_pcp // pricing factor
x2_pcp // pricing factor
x1_lcp // pricing factor
x2_lcp // pricing factor
x3_lcp // pricing factor
x4_lcp // pricing factor
dh // price dispersion for home products
dh_star // price dispersion for home products
A // H-type productivity
pih // PPI inflation of home products
//pif // PPI inflation of foreign products
g // government purchases
nw // net worth of the home currency bank
lev // leverage ratio of the home currency bank
mub // relative cost advantage of borrowing abroad 
muk // expected credit spread
mud // return factor
mu //ratio of mub over muk
kappa // expected total net worth over current net worth, H-type bank
x   //ratio of foreign currency debt over total assets 
d // deposit in the H-type bank
z // claims issued by the intermediate goods producer
rr // real interest rate
Omega // bankers' discounting factor
Thetax // proportion of divertible funds in home banks

// variables of Foreign country

c_star // consumption
l_star // hours worked
w_star // real wage
ow_star // optimal wage setting
xw1_star // optimal wage setting factor
xw2_star // optimal wage setting factor
lambda_star // marginal utility of consumption
k_star // capital
rk_star // rental rate on capital
q_star // discounting factor for investment
inv_star // investment
r_star // gross nominal interest rate
pi_star // CPI inflation rate
yd_star // total demand
yh_star // home goods demand
yf_star // foreign goods demand
y_star // output
pm_star // price of the intermediate goods
ph_star // relative price of the home good bundle over the CPI index
pf_star // relative price of the foreign good bundle over the CPI index
opf_star // optimal pricing of foreign goods in foreign over the foreign PPI
x3_star // pricing factor
x4_star // pricing factor
df_star // price dispersion for foreign products
A_star // productivity
pih_star // PPI inflation of home products
pif_star // PPI inflation of foreign products
g_star // government purchases
d_star // deposit in the bank
nw_star // net worth of the bank
lev_star // leverage ratio of the bank
theta_star // strength of the IC for banks
mu0_star // expected credit spread
mu1_star // return factor
b_star // interbank lending to home F-type banks
psi_star // expected net worth over current net worth
rb_star // interbank lending rate
z_star // claims issued by the intermediate good producer
rr_star // real interest rate

// exchange rates

ner // nominal exchange rate
rer // real exchange rate
tot // P_F/P_H
Deltah // price passthrough ner*P_H^*/P_H
delta1 // factor for price passthrough
delta2 // factor for price passthrough

// report variables

exp_spr_star // expected bank premium spread
spr_star // bank premium spread
r_ann // annualised nominal interest rate
r_star_ann // annualised nominal interest rate
pi_ann // annualised CPI inflation
pih_ann // annualised PPI inflation
lev_mkt_star // market value of leverage
b_star_h // interbank lending in Home currency
ca // current account
exp_spr_star_est // expected bank premium spread scaled for estimation
pi_star_ann_est // annualised PPI inflation scaled for estimation
pi_ann_est  // annualised CPI inflation scaled for estimation
r_ann_est // annualised nominal interest rate scaled for estimation
eqp_est // equity price scaled for estimation
realGDP_est // real Home GDP scaled for estimation
realGDP_star_est // real Foreign GDP scaled for estimation
realexr_est // real exchange rate scaled for estimation
realexp_est; // real export scaled for estimation

varexo

eA // productivity shock
eA_star // productivity shock
eR // monetary policy shock
eR_star //monetary policy shock
etheta_star // theta shock
eG // government spending shock
eG_star; // government spending shock

parameters

// home parameters
h // habit formation
sigma // risk aversion
chi // relative weight on the disutility of labour
zeta // inverse of frisch elasticity of labour supply
beta // discounting coefficient
openness // degree of openness
a // home bias
n // population ratio
epsilon // elasticity of substitution between home and foreign goods
rho // elasticity of substitution across intermediate goods
nul // differentiated wage markup
alfa // share of capital in the production of intermediate goods
phik // cost of investment adjustment
delta // depreciation rate
xi // price rigidity
xiw // wage rigidity
rhoA // persistence of productivity shocks
rhott // persistence of the shock on theta
eta // proportion of home currency firms
eta_f // proportion of PCP firms
xib // proportion of assets brought by new bankers
theta // proportion of divertible funds
gamma_b // increase in the proportion of divertible funds with the increase in proportion of foreign currency denominated debt
omega // survival rate of bankers
gamma // smoothing coeffient in Taylor rule
thetapi // response parameter to inflation in Taylor rule
thetaxr // response parameter to XR in Taylor rule
rhog // persistence of government spending shock
varphi_k // persistence of macroprudential tax shock
varphi_bg // persistence of macroprudential tax shock

// foreign parameters
h_star // habit formation
sigma_star // risk aversion
chi_star // relative weight on the disutility of labour
zeta_star // inverse of frisch elasticity of labour supply
beta_star // discounting coefficient
a_star // home bias
epsilon_star // elasticity of substitution between home and foreign goods
rho_star // elasticity of substitution across intermediate goods
nul_star // differentiated wage markup
alfa_star // share of capital in the production of intermediate goods
phik_star // cost of investment adjustment
delta_star // depreciation rate
xi_star // price rigidity
xiw_star // wage rigidity
rhoA_star // persistence of productivity shocks
rhott_star // persistence of the shock on theta
xib_star // proportion of assets brought by new bankers
omega_star  // steady state proportion of divertible funds
gamma_star // smoothing coeffient in Taylor rule
thetapi_star // response parameter to inflation in Taylor rule
rhog_star // persistence of government spending shock

// calibration targets
l_bar // steady state labour hours
rk_bar // steady state return for capital
rk_star_bar // steady state return for capital
lev_bar // steady state leverage ratio
x_bar // steady state proportion of foreign currency denominated debt
lev_star_bar // steady state leverage ratio

// steady state values
r_bar // steady state home real interest rate
r_star_bar //steady state foreign interest rate
theta_star_bar // steady state theta
g_bar // steady state government purchase
g_star_bar // steady state government purchase
gdp_bar; // steady state Foreign GDP

@#include "baseline_param_AQ.m"
@#include "RECALIB_AQ.mod"

model;

// Households, Home

(exp(c)-h*exp(c(-1)))^(-sigma)=exp(lambda); // FOC for consumption
exp(lambda)=beta*exp(r)*exp(lambda(+1))/(1+pi(+1)); // FOC for deposit holdings

// Households, Foreign

(exp(c_star)-h_star*exp(c_star(-1)))^(-sigma_star)=exp(lambda_star); // FOC for consumption
exp(lambda_star)=beta_star*exp(r_star)*exp(lambda_star(+1))/(1+pi_star(+1)); // FOC for deposit holdings

// International Risk Sharing

exp(yd)+(exp(rb_star(-1))/(1+pi_star))*exp(b_star(-1))*exp(rer)=exp(b_star)*exp(rer)+exp(ph)*exp(yh)+exp(ph_star)*exp(rer)*exp(yh_star); // budget constraint for Home households
exp(rer)/exp(rer(-1))=exp(ner)/exp(ner(-1))*(1+pi_star)/(1+pi); //definition of ner

// The Labour Market, Home

exp(ow)^(1+nul*zeta)=chi*nul/(nul-1)*exp(xw1)/exp(xw2); // optimal wage setting
exp(xw1)=(exp(l)^(1+zeta))*(exp(w)^(zeta*nul))+beta*xiw*((exp(w(+1))/exp(w))^nul)*((1+pi(+1))^(nul*(1+zeta)))*exp(xw1(+1)); // optimal wage setting factor
exp(xw2)=exp(l)*exp(lambda)+beta*xiw*((exp(w(+1))/exp(w)*(1+pi(+1)))^nul)*exp(xw2(+1))/(1+pi(+1)); // optimal wage setting factor
exp(w)=((1-xiw)*(exp(ow)^(1-nul))+xiw*((exp(w(-1))/(1+pi))^(1-nul)))^(1/(1-nul)); // law of motion for real wage

// The Labour Market, Foreign

exp(ow_star)^(1+nul_star*zeta_star)=chi_star*nul_star/(nul_star-1)*exp(xw1_star)/exp(xw2_star); // optimal wage setting
exp(xw1_star)=(exp(l_star)^(1+zeta_star))*(exp(w_star)^(zeta_star*nul_star))+beta_star*xiw_star*((exp(w_star(+1))/exp(w_star))^nul_star)*((1+pi_star(+1))^(nul_star*(1+zeta_star)))*exp(xw1_star(+1)); // optimal wage setting factor
exp(xw2_star)=exp(l_star)*exp(lambda_star)+beta_star*xiw_star*((exp(w_star(+1))/exp(w_star)*(1+pi_star(+1)))^nul_star)*exp(xw2_star(+1))/(1+pi_star(+1)); // optimal wage setting factor
exp(w_star)=((1-xiw_star)*(exp(ow_star)^(1-nul_star))+xiw_star*((exp(w_star(-1))/(1+pi_star))^(1-nul_star)))^(1/(1-nul_star)); // law of motion for real wage

// Local Banks

exp(lev)=exp(mud)/(exp(Thetax)-exp(muk)-exp(mub)*exp(x)); // optimal leverage of the bank
exp(lev)=exp(q)*exp(z)/exp(nw); // definition of the leverage
exp(mub)=exp(Omega)*(exp(r)/(1+pi(+1))-exp(rb_star)/(1+pi_star(+1))*exp(rer(+1))/exp(rer)); // definition of mub
exp(muk)=exp(Omega)*(exp(rk(+1))-exp(r)/(1+pi(+1))); // definition of muk
exp(mud)=exp(Omega)*exp(r)/(1+pi(+1)); // definition of mud
exp(q)*exp(z)=exp(nw)+exp(d)+exp(rer)*exp(b_star); // balance constraint
exp(nw)=(omega+xib)*exp(rk)*exp(q(-1))*exp(z(-1))-omega*(exp(r(-1))/(1+pi))*exp(d(-1))-omega*(exp(rb_star(-1))/(1+pi_star))*exp(b_star(-1))*exp(rer); // evolution of net worth
exp(kappa)=exp(Thetax)*exp(lev); // definition of kappa
exp(x)=exp(b_star)*exp(rer)/(exp(q)*exp(z)); //definition of x
exp(x)=(sqrt(1+(2/gamma_b)*(exp(mu)^2))-1)/exp(mu); //optimal proportion of interbank borrowing
exp(mu)=exp(mub)/exp(muk); //definition of mu
exp(Omega)=beta*(exp(lambda(+1))/exp(lambda))*(1-omega+omega*exp(kappa(+1))); //definition of Omega
exp(Thetax)=theta*(1+gamma_b/2*(exp(x)^2)); //definition of Thetax

// Global Banks

exp(lev_star)=exp(mu1_star)/(theta_star-exp(mu0_star)); // optimal leverage of the bank
exp(lev_star)=(exp(q_star)*exp(z_star)+exp(b_star))/exp(nw_star); // definition of the leverage
exp(mu0_star)=beta_star*(exp(lambda_star(+1))/exp(lambda_star))*(1-omega_star+omega_star*exp(psi_star(+1)))*(exp(rk_star(+1))-exp(r_star)/(1+pi_star(+1))); // definition of mu0
exp(mu1_star)=beta_star*(exp(lambda_star(+1))/exp(lambda_star))*(1-omega_star+omega_star*exp(psi_star(+1)))*(exp(r_star)/(1+pi_star(+1))); // definition of mu1
beta_star*(exp(lambda_star(+1))/exp(lambda_star))*(1-omega_star+omega_star*exp(psi_star(+1)))*exp(rk_star(+1))=beta_star*(exp(lambda_star(+1))/exp(lambda_star))*(1-omega_star+omega_star*exp(psi_star(+1)))*exp(rb_star)/(1+pi_star(+1)); // optimal rate for interbank lending
exp(q_star)*exp(z_star)+exp(b_star)=exp(nw_star)+exp(d_star); // balance constraint
exp(nw_star)=(omega_star+xib_star)*(exp(rk_star)*exp(q_star(-1))*exp(z_star(-1))+exp(rb_star(-1))/(1+pi_star)*exp(b_star(-1)))-omega_star*(exp(r_star(-1))/(1+pi_star))*exp(d_star(-1)); // evolution of net worth
exp(psi_star)=theta_star*exp(lev_star); // definition of psi
theta_star=(theta_star_bar^(1-rhott))*(theta_star(-1)^rhott)/exp(etheta_star); // exogenous evolution process of theta

// Intermediate Goods Packagers, Home

exp(pm)=((exp(q(-1))*exp(rk)-(1-delta)*exp(q))^alfa)*(exp(w)^(1-alfa))/(exp(A)*(alfa^alfa)*((1-alfa)^(1-alfa))); // determine pm
exp(A)=(exp(A(-1))^rhoA)/exp(eA); // law of motion for technology
exp(w)=(1-alfa)*exp(pm)*(exp(dh)*exp(yh)+exp(dh_star)*exp(yh_star))/exp(l); // FOC for labour
exp(rk)=(alfa*exp(pm)*(exp(dh)*exp(yh)+exp(dh_star)*exp(yh_star))/exp(k(-1))+(1-delta)*exp(q))/exp(q(-1)); // FOC for capital
exp(z)=exp(k); // exchange of claims and capital

// Intermediate Goods Producers, Foreign

exp(pm_star)=((exp(q_star(-1))*exp(rk_star)-(1-delta_star)*exp(q_star))^alfa_star)*(exp(w_star)^(1-alfa_star))/(exp(A_star)*(alfa_star^alfa_star)*((1-alfa_star)^(1-alfa_star))); // determine pm
exp(A_star)=(exp(A_star(-1))^rhoA_star)/exp(eA_star); // law of motion for technology
exp(w_star)=(1-alfa_star)*exp(pm_star)*exp(df_star)*exp(y_star)/exp(l_star); // FOC for labour
exp(rk_star)=(alfa_star*exp(pm_star)*exp(df_star)*exp(y_star)/exp(k_star(-1))+(1-delta_star)*exp(q_star))/exp(q_star(-1)); // FOC for capital
exp(z_star)=exp(k_star); // exchange of claims and capital

// Capital Producers, Home

exp(q)=1+(phik/2)*((exp(inv)/exp(inv(-1))-1)^2)+phik*(exp(inv)/exp(inv(-1))-1)*(exp(inv)/exp(inv(-1)))-phik*beta*(exp(lambda(+1))/exp(lambda))*(exp(inv(+1))/exp(inv)-1)*((exp(inv(+1))^2)/(exp(inv)^2)); // FOC for investment
exp(k)=exp(inv)+(1-delta)*exp(k(-1)); // capital accumulation

// Capital Producers, Foreign

exp(q_star)=1+(phik_star/2)*((exp(inv_star)/exp(inv_star(-1))-1)^2)+phik_star*(exp(inv_star)/exp(inv_star(-1))-1)*(exp(inv_star)/exp(inv_star(-1)))-phik_star*beta_star*(exp(lambda_star(+1))/exp(lambda_star))*(exp(inv_star(+1))/exp(inv_star)-1)*((exp(inv_star(+1))^2)/(exp(inv_star)^2)); // FOC for investment
exp(k_star)=exp(inv_star)+(1-delta_star)*exp(k_star(-1)); // capital accumulation

// demand aggregation

exp(yh)=a*(exp(ph)^(-epsilon))*exp(yd); // demand for home products
exp(yf)=(1-a)*(n/(1-n))*(exp(pf)^(-epsilon))*exp(yd); // demand for foreign products
exp(y)=exp(yh)+exp(yh_star); // total output, home
exp(yf_star)=(1-a_star)*(exp(pf_star)^(-epsilon_star))*exp(yd_star); // demand for foreign products
exp(yh_star)=a_star*((1-n)/n)*(exp(ph_star)^(-epsilon_star))*exp(yd_star); // demand for home products
exp(y_star)=exp(yf)+exp(yf_star); // total demand

// PCP Final Goods Producers, Home

exp(oph_pcp)=(rho/(rho-1))*(exp(x1_pcp)/exp(x2_pcp)); // optimal pricing for home goods
exp(x1_pcp)=exp(y)*exp(pm)*exp(lambda)+xi*beta*((1+pih(+1))^rho)*exp(x1_pcp(+1)); // definition of x1
exp(x2_pcp)=exp(y)*exp(lambda)*exp(ph)+xi*beta*((1+pih(+1))^(rho-1))*exp(x2_pcp(+1)); // definition of x2
exp(oph_pcp_star)=exp(ph)/(exp(rer)*exp(ph_star))*exp(oph_pcp); // law of one price from PCP firms

// LCP Final Goods Producers, Home

exp(oph_lcp)=(rho/(rho-1))*(exp(x1_lcp)/exp(x2_lcp)); // optimal pricing for home goods
exp(x1_lcp)=exp(yh)*exp(pm)*exp(lambda)+xi*beta*((1+pih(+1))^rho)*exp(x1_lcp(+1)); // definition of x1
exp(x2_lcp)=exp(yh)*exp(lambda)*exp(ph)+xi*beta*((1+pih(+1))^(rho-1))*exp(x2_lcp(+1)); // definition of x2
exp(oph_lcp_star)*exp(rer)=(rho/(rho-1))*(exp(x3_lcp)/exp(x4_lcp)); // optimal pricing for home goods
exp(x3_lcp)=exp(yh_star)*exp(pm)*exp(lambda)+xi*beta*((1+pih_star(+1))^rho)*exp(x3_lcp(+1)); // definition of x3
exp(x4_lcp)=exp(yh_star)*exp(lambda)*exp(ph_star)+xi*beta*(exp(ner(+1))/exp(ner))*((1+pih_star(+1))^(rho-1))*(1+pi_star(+1))*exp(x4_lcp(+1))/(1+pi(+1)); // definition of x4

// demand, price evolution and dispersions for home producers

xi*(1+pih)^(rho-1)+(1-xi)*(eta_f*(exp(oph_pcp)^(1-rho))+(1-eta_f)*(exp(oph_lcp)^(1-rho)))=1; // law of motion for prices
xi*((1+pih_star)^(rho-1))+(1-xi)*(eta_f*(exp(oph_pcp_star)^(1-rho))+(1-eta_f)*(exp(oph_lcp_star)^(1-rho)))=1; // law of motion for prices
exp(dh)=(1-xi)*(eta_f*(exp(oph_pcp)^(-rho))+(1-eta_f)*(exp(oph_lcp)^(-rho)))+xi*((1+pih)^rho)*exp(dh(-1)); // law of motion for price dispersion
exp(dh_star)=(1-xi)*(eta_f*(exp(oph_pcp_star)^(-rho))+(1-eta_f)*(exp(oph_lcp_star)^(-rho)))+xi*((1+pih_star)^rho)*exp(dh_star(-1)); // law of motion for price dispersion

// Final Goods Producers, Foreign

exp(opf_star)=(rho_star/(rho_star-1))*(exp(x3_star)/exp(x4_star)); // optimal pricing for foreign goods
exp(x3_star)=exp(y_star)*exp(pm_star)*exp(lambda_star)+xi_star*beta_star*((1+pif_star(+1))^rho_star)*exp(x3_star(+1)); // definition of x3_star
exp(x4_star)=exp(y_star)*exp(lambda_star)*exp(pf_star)+xi_star*beta_star*((1+pif_star(+1))^(rho_star-1))*exp(x4_star(+1)); // definition of x4_star
xi_star*(1+pif_star)^(rho_star-1)+(1-xi_star)*(exp(opf_star)^(1-rho_star))=1; // law of motion for prices
exp(df_star)=(1-xi_star)*(exp(opf_star)^(-rho_star))+xi_star*((1+pif_star)^rho_star)*exp(df_star(-1)); // law of motion for price dispersion

// relative prices and inflation

(1+pih)/(1+pi)=exp(ph)/exp(ph(-1)); //definition of inflation
(1+pif_star)/(1+pi_star)=exp(pf_star)/exp(pf_star(-1)); //definition of inflation
a*(exp(ph)^(1-epsilon))+(1-a)*(exp(pf)^(1-epsilon))=1; //links of relative prices
a_star*(exp(ph_star)^(1-epsilon_star))+(1-a_star)*(exp(pf_star)^(1-epsilon_star))=1; //links of relative prices
exp(ph)^(epsilon-1)=a+(1-a)*exp(tot)^(1-epsilon); //definition of ph
exp(pf_star)^(epsilon_star-1)=a_star*(exp(Deltah)/exp(tot))^(1-epsilon_star)+(1-a_star); //definition of pf_star
exp(rer)^(1-epsilon)=(a_star*(exp(Deltah))^(1-epsilon)+(1-a_star)*(exp(tot))^(1-epsilon))/(a+(1-a)*(exp(tot)^(1-epsilon))); // definition of tot
(exp(Deltah))=eta_f+(1-eta_f)*exp(rer)*exp(delta1)/exp(delta2); // price passthrough
exp(delta1)^(1-rho_star)=xi*(exp(delta1(-1))/(1+pi_star))^(1-rho_star)+(1-xi)*(exp(oph_lcp_star)*exp(ph_star))^(1-rho_star); //definition of delta1
exp(delta2)^(1-rho)=xi*(exp(delta2(-1))/(1+pi))^(1-rho)+(1-xi)*(exp(oph_lcp)*exp(ph))^(1-rho); // definition of delta2

// Monetary Policy and Government Purchases

exp(g)/g_bar=((exp(g(-1))/g_bar)^rhog)*exp(eG); // exogenous government purchase
exp(g_star)/g_star_bar=((exp(g_star(-1))/g_star_bar)^rhog_star)*exp(eG_star); // exogenous government purchase
exp(r)=((exp(r(-1)))^gamma)*((r_bar*((1+pi)^thetapi*(exp(ner)/exp(ner(-1)))^thetaxr))^(1-gamma))*exp(eR); // Taylor Rule Home
exp(r_star)=((exp(r_star(-1)))^gamma_star)*((r_star_bar*((1+pi_star)^thetapi_star))^(1-gamma_star))/exp(eR_star); // Taylor Rule Foreign
exp(rr)=exp(r(-1))/(1+pi); // definition of real interest rate
exp(rr_star)=exp(r_star(-1))/(1+pi_star); // definition of real interest rate

// Market Clearing

exp(yd)=exp(c)+exp(inv)+exp(g)+(phik/2)*((exp(inv)/exp(inv(-1))-1)^2)*exp(inv); // resource constraint for Home
exp(yd_star)=exp(c_star)+exp(inv_star)+exp(g_star)+(phik_star/2)*((exp(inv_star)/exp(inv_star(-1))-1)^2)*exp(inv_star); // resource constraint for Foreign

// variables for presentation
r_ann=(exp(r))^4;  // definition of annualised interest rate
r_star_ann=(exp(r_star))^4;  // definition of annualised interest rate
pi_ann=(1+pi)^4-1;  // definition of annualised CPI inflation rate
pih_ann=(1+pih)^4-1;  // definition of annualised PPI inflation rate
ca=(-exp(b_star)+exp(b_star(-1)))*exp(rer)/gdp_bar; // definition of current account over steady state GDP
spr_star=(exp(rk_star)^4-exp(r_star(-1))^4/(1+pi_star)^4); // definition of bank premium of Foreign banks
exp_spr_star=(exp(rk_star(+1))^4-exp(r_star)^4/(1+pi_star(+1))^4); // definition of bank premium of Foreign banks
exp(lev_mkt_star)=(exp(q_star)*exp(z_star(-1))+exp(b_star(-1)))/exp(nw_star(-1)); // market value of leverage
exp(b_star_h)=exp(b_star)*exp(rer); // interbank lending in Home currency
eqp_est=q*100;
realGDP_est=y*100;
realexr_est=rer*100;
realGDP_star_est=y_star*100;
realexp_est=yh_star*100;
exp_spr_star_est=(exp(rk_star(+1))^4-exp(r_star)^4/(1+pi_star(+1))^4)*10000;
pi_star_ann_est=((1+pi_star)^4-1)*100; 
pi_ann_est=((1+pi)^4-1)*100;  
r_ann_est=((exp(r))^4)*100;
end;

@#include "do_options.mod"

shocks;
@#if f_shock == 1 && f_estim == 1
    var eR_star; stderr 0.0015;
@#endif
@#if f_shock == 1 && f_estim == 0
    var eR_star; stderr 0.0025;
@#endif
@#if f_shock == 2
    var etheta_star; stderr 0.01;
@#endif
@#if f_shock == 3
    var eA_star; stderr 0.01;
@#endif
@#if f_shock == 4
    var etau_k; stderr 0.0001;
@#endif
@#if f_shock == 5
    var etau_fk; stderr 0.0001;
@#endif
@#if f_shock == 6
    var etau_bg; stderr 0.0001;
@#endif
end;

steady_state_model;  // Calculate the steady state
@#include "steadystate_AQ.mod"
end;

resid;
steady;
check;

@#if f_estim == 0
    stoch_simul(order=1, nograph, irf=40) y c inv l r_ann pi_ann pih_ann rer ca q k b_star_h ner r_star_ann lev_star exp_spr_star q_star k_star nw_star yh_star yf lev_mkt_star;
@#else
    @#include "IRFmatch.mod"
@#endif