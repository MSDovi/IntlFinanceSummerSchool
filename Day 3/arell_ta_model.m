%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DESCRIPTION:
% This code solves the Arellano (2008) soveregin default model on a fine grid for output.
% Code by Sergio de Ferra and Enrico Mallucci.


% Variable spread1993 is the simulated time-series for spreads from 1993Q1 
% to 2001Q3 according to the model. Variable spread1993 is the input 
% variable that is used to compute the time-varying volatility of interest
% rate spreads.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;

% Parameters set independently of targets
gamma_c = 2;            % RRA in consumer preferences

% Intertemporal Parameters 
beta     =  0.95282;          % Subjective discount factor (Cuadra, Sanchez, Sapriza = .97) Arellano: .953
lambda   = .282;          % Probability of re-entering sovereign debt contract after exclusion, Arellano: .282

def_cost_asymm = 1; % = 2 if default cost is quadratic; 1 if default cost a la arellano; 0 if aguiar gopinath

wc_par_symm    = .969;          % Welfare cost of def't: Fraction of output available in default state. Not piecewise.
wc_par_asymm   = .969;          % Welfare cost of def't: Fraction of output available in default state. Piecewise.
cost1 = -0.195;
cost2 = 0.145;
    
% r-process parameters
mu_r    = (1.017)-1;

% y-process parameters
sigma_ey = .025;%.006;         % Standard deviation of innovation of y. Arellano = 0.025
rho_y    = .945;%.98;          % AR coefficient for y. Arellano = 0.945
    
% Number of y-gridpoints
N_y     = 556; 

% Number of s.d. spanned by processes. If N_variable = 1, degenerate distribution at mean
int_y       = 2*(N_y~=1) + eps*(N_y==1);    



%% Discretize Process for y and generate transition matrix
% Tauchen (1986) method, as implemented in code by Martin Floden

mu= -.5*sigma_ey^2/(1-rho_y^2);
rho = rho_y;
sigma = sigma_ey;

m=int_y; 
Z     = zeros(N_y,1);
Zprob = zeros(N_y,N_y);
a     = (1-rho)*mu;

Z(N_y)  = m * sqrt(sigma^2 / (1 - rho^2));
Z(1)  = -Z(N_y);
zstep = (Z(N_y) - Z(1)) / (N_y - 1);

for i=2:(N_y-1)
    Z(i) = Z(1) + zstep * (i - 1);
end 

Z = Z + a / (1-rho);

for j = 1:N_y
    for k = 1:N_y
        if k == 1
            Zprob(j,k) = 0.5 * erfc(-((Z(1) - a - rho * Z(j) + zstep / 2) / sigma)/sqrt(2));
        elseif k == N_y
            Zprob(j,k) = 1 - 0.5 * erfc(-((Z(N_y) - a - rho * Z(j) - zstep / 2) / sigma)/sqrt(2));
        else
            Zprob(j,k) = 0.5 * erfc(-((Z(k) - a - rho * Z(j) + zstep / 2) / sigma)/sqrt(2)) - ...
                         0.5 * erfc(-((Z(k) - a - rho * Z(j) - zstep / 2) / sigma)/sqrt(2));
        end
    end
end

% Vector of log(y)-values
ly_vec =Z;

% Transition matrix
P_y = Zprob;

% Vector of (y)-values
y_vec = exp(ly_vec);

% Mean of output process
y_mean = P_y^100000*y_vec;
y_mean = y_mean(floor(N_y/2)+1);

% Make yearly interest rates into quarterly debt prices, vector sized (N_y,1)
qrf_vec = (1 + repmat(mu_r,[N_y,1])).^-(1/4); 

%% Define utility under autarky

if def_cost_asymm==1
c_aut = y_vec;    
c_aut(y_vec>wc_par_asymm*y_mean) = wc_par_asymm*y_mean;
elseif def_cost_asymm==0
c_aut = wc_par_symm*y_vec;
elseif def_cost_asymm==2
c_aut = y_vec - max(0,cost1.*y_vec+cost2.*y_vec.^2);
end

util_aut = 1/(1-gamma_c)*c_aut.^(1-gamma_c);

%% Setup dynamic sovereign default problem

% Parameters for iteration
damp_v        = .8;          % dampening in vfi
damp_q        = .8;          % dampening in q-update
maxiter_v   = 500;
maxiter_q   = 600;
tol_v       = 1e-4;
tol_q       = 1e-8;
log_diff_q      = [];

iter_v  = 1;
iter_q  = 1;
diff_v  = 7;
diff_q  = 7;

loadguesses     = 0;    % 1 if loading guesses from external file

% Define grid for government debt ------

single_space    = 0;    % 1 if grid is evenly spaced, 0 if 2-steps grid

b_g_min   = -0.05;
b_g_max   = 1.18;
N_b_g       = round(456);

if single_space
    
    b_g_vec = linspace(b_g_min,b_g_max,N_b_g);
    
else
    
    b_g_mid = .3;
    N_b_g_mid = floor(.85*N_b_g);
    
    
low_bound =b_g_min;
mid_point =b_g_mid;
upp_bound =b_g_max;
N_1 =N_b_g_mid;
N_tot =N_b_g;

if N_1 > 0 && N_1 < N_tot && ( (low_bound < mid_point && mid_point < upp_bound) || (low_bound > mid_point && mid_point > upp_bound) )

y_1 = linspace(low_bound,mid_point,N_1);

y_2 = linspace(mid_point+(upp_bound-mid_point)/(N_tot-N_1),upp_bound,N_tot-N_1);

b_g_vec = [y_1 y_2];

else
    
b_g_vec = linspace(low_bound,upp_bound,N_tot);

end

end

% Ensure presence of 0 on b-vector
[~,i_b_g_zero] = min(abs(b_g_vec));
b_g_vec(i_b_g_zero) = 0;

% 3-dimensional grid for state-debt. Dim: x,b,b'
b_state_grid = repmat(b_g_vec,[N_y,1,N_b_g]);

% 2-dimensional grid for choice-debt. Dim: x,b'
b_choice_2grid = repmat(b_g_vec,[N_y,1]);

% 3-dimensional grid for y-state. Dim: x,b,b'
y_state_grid = repmat(y_vec,[1,N_b_g,N_b_g]);

% -------------------------------------------------
% Start from scratch or load guesses
loadguess = true;

% Pre-allocate grids
if loadguess == true
    
    guesses = load('arell_fine_replication');
    
    q_g = guesses.q_g;
    v_guess = guesses.v_guess;
    v_bad_guess = guesses.v_bad_guess ;
    i_b_max_new = guesses.i_b_g_pf ;
    def_new= guesses.def_pf;
    
else

    
    % Form initial guesses
    q_g = qrf_vec*ones(1,N_b_g);
    v_guess = zeros(N_y,N_b_g);
    v_bad_guess = zeros(N_y,1);
    i_b_max_new = NaN(N_y,N_b_g);
    def_new= NaN(N_y,N_b_g);
end


%% Sovereign default model solution by iteration

while diff_q > tol_q && iter_q < maxiter_q


while diff_v > tol_v && iter_v<maxiter_v
    
    % Expected continuation value
    e_v_guess = P_y*v_guess;
    e_v_guess_3grid = permute(repmat(e_v_guess,[1 1 N_b_g]),[1 3 2]);
    
    % value if re-entering contract with no debt
    v_badgood_guess = v_guess(:,i_b_g_zero);
    e_v_badgood_guess = P_y*v_badgood_guess;
    
    % Expected value of remaining in exclusion state
    e_v_bad_guess = P_y*v_bad_guess;
    
    % Resources from borrowing
    borr_rev_choice = q_g.*b_choice_2grid;
    borr_rev_choice_3grid = permute(repmat(borr_rev_choice,[1,1,N_b_g]),[1 3 2]); % Make grid with dimension z,b,b'
    
    % Consumption implied by state and choice for b'. Une a guess for price index. Dim: z,b,b'
    cons_choice = y_state_grid - b_state_grid + borr_rev_choice_3grid;
    cons_choice(cons_choice<eps) = eps;
    
    % Period-utility implied by state and choice for b'. Dim: z,b,b'
    util_choice = 1/(1-gamma_c)*cons_choice.^(1-gamma_c);
        
    % Form maximand in borrowing choice
    borrower_maximand = util_choice + beta*e_v_guess_3grid;
    
    % Store old choice for borrowing
    i_b_max_old = i_b_max_new;
    
    % Choice for borrowing and update value function
    [v_good_guess_new, i_b_max_new] = max(borrower_maximand,[],3);
        
    
    % Evaluate value function in default state
    v_bad_guess_new = util_aut + beta*(lambda*e_v_badgood_guess + (1-lambda)*e_v_bad_guess);
    
    %Store old default policy function
    def_old = def_new;
    
    % Evaluate policy function for default
    def_new = repmat(v_bad_guess_new,[1,N_b_g])>v_good_guess_new;
    
    
    % Evaluate value function including discrete choice for default
    v_guess_new = max(repmat(v_bad_guess_new,[1,N_b_g]),v_good_guess_new);
    
    % Store old value-function guesses
    v_old = v_guess;
    v_bad_old = v_bad_guess;
    
    % Update value-function guesses
    v_guess = damp_v*v_guess_new + (1-damp_v)*v_old;
    v_bad_guess = damp_v*v_bad_guess_new + (1-damp_v)*v_bad_old;
    
    % Evaluate change in v
    diff_v = max(abs(v_guess_new(:)-v_old(:)))
    
    iter_v = iter_v+1;
    
end

% Default probability conditional on debt choice
e_def = P_y*def_new;

% Compute new sovereign debt price
q_g_new = repmat(qrf_vec,[1,N_b_g]).*(1-e_def);

% Store old sovereign debt price
q_g_old = q_g;

% Evaluate change in q_g
diff_q = max(abs(q_g_new(:)-q_g_old(:)))

% Update sovereign debt price
q_g = damp_q*q_g_new + (1-damp_q)*q_g_old;



iter_q

iter_q = iter_q + 1;

iter_v

% Reset counter and diff for inner-v loop
diff_v = 7;
iter_v = 1;


end


%% Define PF's for interesting variables, and reshape in a reasonable way

% Value functions
v       = v_guess;
v_good  = v_good_guess_new;
v_bad   = v_bad_guess;

%Default
def_pf = def_new;

% Debt - index
i_b_g_pf = i_b_max_new;
i_b_g_pf(def_pf) = i_b_g_zero;

% Policy function for q_g
q_g_pf = NaN(N_y,N_b_g);
for i_x = 1:N_y
q_g_pf(i_x,:) = (q_g(i_x,i_b_g_pf(i_x,:)));
end

q_g_pf(def_pf) = NaN;

% PF for q_g, with zeros where there would be default
q_g_pf_zeros = q_g_pf;
q_g_pf_zeros(isnan(q_g_pf)) = 0;

% Policy function for r_g
r_g_pf = q_g_pf.^-1-1;

% PF for r_g, with zeros where there would be default
r_g_pf_zeros = r_g_pf;
r_g_pf_zeros(isnan(q_g_pf)) = 0;

% Policy function for b_g
b_g_pf = b_g_vec(i_b_g_pf);
b_g_pf(def_pf)  = NaN;

% Policy function for CA
ca_pf = -(b_g_vec(i_b_g_pf) - repmat(b_g_vec,[N_y,1]));
ca_pf(def_pf)  = NaN;


% Policy function for TB
tb_pf = -(q_g_pf_zeros.*b_g_vec(i_b_g_pf) - repmat(b_g_vec,[N_y,1]));
tb_pf(def_pf)  = 0; % Impose balanced trade when defaulting


% Reshape policy functions into grid
% q_g (function)
q_g_grid =reshape(q_g,[N_y,N_b_g]);

% q_g (policy function)
q_g_pf_grid     =reshape(q_g_pf,[N_y,N_b_g]);

% Default
def_pf_grid     =reshape(def_pf,[N_y,N_b_g]);

% CA
ca_pf_grid      =reshape(ca_pf,[N_y,N_b_g]);

% Debt
b_g_pf_grid =reshape(b_g_pf,[N_y,N_b_g]);



%% Load empirical data for Argentina


[data,~,~]=xlsread('data_arell_ta.xlsx', 'data','B56:G91');
spreads_data = data(:,6);
nx_y_data =data(:,5);

[y_data,~,~]=xlsread('data_arell_ta.xlsx', 'data','E17:E96');

[~,year_labels,~]=xlsread('data_arell_ta.xlsx', 'data','A56:A91');


gdp=detrend(log(y_data));
gdp = gdp(40:end-5,:); %Only consider shocks between 1993Q1 and 2001Q3
i_y_path=zeros(size(gdp,1),1);

% Find path of i_y that closely replicates Argentine data
for i=1:size(gdp,1)
    prova=gdp(i,1);
    error =(prova-ly_vec).^2;
    [~, i_y_close] = min(error);
    i_y_path(i,1)=i_y_close;

end




%% Simulate model

% Number of simulated periods
T_sim = 5000;

% Initial value of exogenous variable
init_x_sim = sub2ind([1 N_y],floor(N_y/2)+1);

% Number of simulations
J=10;

% Pre-allocate for speed
i_b_g_sim       = NaN(T_sim,J); % Sovereign debt series (index)
def_sim         = NaN(T_sim,J); % Default series
q_g_sim         = NaN(T_sim,J); % Sovereign debt price series
r_g_sim         = NaN(T_sim,J); % Sovereign debt int.rate series
y_sim           = NaN(T_sim,J); % Output series
tb_sim          = NaN(T_sim,J); % Trade balance series
ca_sim          = NaN(T_sim,J); % Current account series


i_x_sim = NaN(T_sim,J);
access_sim      = NaN(T_sim,J); % Access to int'l financial markets

i_b_g_sim(1,:)    = i_b_g_zero.*ones(1,J);
access_sim(1,:)   = ones(1,J);

for j=1:J

    
% Simulated shock process - - - - - -- - - - -- - - - - - - -- - - - - - 
% Uses code from markov.m file, from Ljungqvist-Sargent textbook
    
T=P_y;
n=T_sim-size(gdp,1);
s0=init_x_sim;
V=1:N_y;


[r, c]=size(T);

if r ~= c
  disp('error using markov function');
  disp('transition matrix must be square');
  return;
end

for k=1:r
  if sum(T(k,:)) ~= 1
    disp('error using markov function')
    disp(['row ',num2str(k),' does not sum to one']);
    disp(['normalizing row ',num2str(k),'']);
    T(k,:)=T(k,:)/sum(T(k,:));
  end
end
[v1, v2]=size(V);
if v1 ~= 1 ||v2 ~=r
  disp('error using markov function');
  disp(['state value vector V must be 1 x ',num2str(r),''])
  if v2 == 1 &&v2 == r
    disp('transposing state valuation vector');
    V=V';
  else
    return;
  end  
end
if s0 < 1 ||s0 > r
  disp(['initial state ',num2str(s0),' is out of range']);
  disp(['initial state defaulting to 1']);
  s0=1;
end

X=rand(n-1,1);
s=zeros(r,1);
s(s0)=1;
cum=T*triu(ones(size(T)));

for k=1:length(X)
  state(:,k)=s;
  ppi=[0 s'*cum];
  s=((X(k)<=ppi(2:r+1)).*(X(k)>ppi(1:r)))';
end
i_x_sim1=(V*state)';


i_x_sim(:,j) = [init_x_sim; i_x_sim1; i_y_path];


% Simulated redemption process
redem_sim = (rand(T_sim+1,1)<lambda);
    
y_sim(:,j) = y_vec(i_x_sim(:,j));

% End of simulated shock process - - - - - -- - - - -- - - - - - - -- - - -


% Simulated process for endogenous variables
for t = 1:T_sim
    
    if access_sim(t,j) % Access to international financial markets
        
        % Choose default/repayment
        def_sim(t,j) = def_pf(i_x_sim(t,j),i_b_g_sim(t,j));
        
        % Lose next-period market access if defaulting and not redempted tomorrow
        access_sim(t+1,j) = (1-def_sim(t,j)*(1-redem_sim(t+1)));
        
        % Choose next-period debt (=0 if defaulting)
        i_b_g_sim(t+1,j) = i_b_g_pf(i_x_sim(t,j),i_b_g_sim(t,j));
        
        % Additional fiscal variables
        q_g_sim(t,j) = q_g_pf(i_x_sim(t,j),i_b_g_sim(t,j));
        
        tb_sim(t,j) = tb_pf(i_x_sim(t,j),i_b_g_sim(t,j));
        
        ca_sim(t,j) = ca_pf(i_x_sim(t,j),i_b_g_sim(t,j));
        
        r_g_sim(t,j) = r_g_pf(i_x_sim(t,j),i_b_g_sim(t,j));

        
    else    % Financial autarky
        
        def_sim(t,j) = 0;
        access_sim(t+1,j) = redem_sim(t+1);
        
        i_b_g_sim(t+1,j) = i_b_g_zero;
        
        tb_sim(t,j) = 0;
%         
    end
    
end

% Sovereign debt
b_g_sim(:,j) = b_g_vec(i_b_g_sim(:,j));

end

r_g_sim = q_g_sim.^-4-1;

spread_sim = 10000*(r_g_sim-mu_r);

%% Compare simulated with empirical data

set(groot, 'defaultTextInterpreter','latex'); 
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(0, 'DefaultAxesFontSize', 18);




f = figure;
f.Position = [1,1,800,450];
subplot(3,1,1)
hold on;
plot(mean(y_sim(end-size(gdp,1)+1:end,:)'));
plot(exp(gdp),'--r')
xticks(1:4:36)
xticklabels(year_labels(1:4:36))
legend('Model','Data')
title('Detrended Log GDP')


subplot(3,1,2)
hold on;
plot(nanmean(spread_sim(end-size(gdp,1)+1:end,:)'));
plot(spreads_data*100,'--r')
xticks(1:4:36)
xticklabels(year_labels(1:4:36))
title('Sovereign Spreads')

subplot(3,1,3)
hold on;
plot(100*mean(tb_sim(end-size(gdp,1)+1:end-1,:)./y_sim(end-size(gdp,1)+1:end-1,:),2))
plot(nx_y_data*100,'--r')
xticks(1:4:36)
xticklabels(year_labels(1:4:36))
title('Trade Balance to Output Ratio')

%% Figures for lecture

% Default / Repyament Policy Function

f = figure;
% f.Position = [1,1,1000,450];
mesh(b_g_vec/4,y_vec,def_pf)
xlim([b_g_vec(1)/4-.0001,.2])
ylim([y_vec(1)-.001,y_vec(end)])
xlabel('Debt')
ylabel('Output')
view( 0 , 90)
map = [0, 0, 0.3;
  255/256, 51/256, 0];
colormap(map)
cbar = colorbar;
cbar.Ticks = [.25 .75] ; %Create 8 ticks from zero to 1
cbar.TickLabels = {'Repay','Default'}
cbar.TickLabelInterpreter  = 'latex';
cbar.TickDirection ='out';
title('Choice of Default/Repayment')
%% Cost of borrowing

f = figure;
f.Position = [1,1,800,450];
hold on;
plot(b_g_vec/4,q_g(floor(N_y/2)+1,:),'Linewidth',2)
plot(b_g_vec/4,q_g(95,:),'Linewidth',2)
xlim([b_g_vec(1)/4-.001,.1])
title('Price of Debt, $q(y,b^\prime)$')
xlabel('Debt Issued, $b^\prime$')
ylabel('$q$')
legend('Mean Output','Low Output')