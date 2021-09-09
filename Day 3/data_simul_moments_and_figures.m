%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Description:
%This file contains the codes that produce graph and moments reported in the paper

%Instructions:
%Only run this code if DataModel.xlsx,longjie_replication.mat, shortjie_replication.mat, 
%    and long_quadjie_replication.mat  are  up to date. 


%If databases are not up-to-date you need to perform the following three
%steps before running this code:
% 1: Generate databases Databases/longjie_replication.mat, Databases/shortjie_replication.mat, 
%    and Databases/long_quadjie_replication.mat using the codes esd_tvv_longjie.m, 
%    esd_tvv_shortjie.m
%2: Generate parameter estimates and error estimates for the volatility process
%   of interes rate spreads, using the fortran code
%3: Store parameter and error estimated from fortran in the excel file
%   Databases/DataModel.xlsx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Moments for the volatility of interest rates from the DATA
clear all
%Load errors for the process of the volatility of interest rate spreads obtained from the econometric model
[errors,~,~]=xlsread('P:\Research\20160628_default_volatility\Codes\Replication Codes\Databases\DataModel.xlsx', 'ARGspreadERRORS');
%Load Parameter estimates for process of the volatility of interest rate spreads, obtained from the econometric model
[parameters,~,~]=xlsread('P:\Research\20160628_default_volatility\Codes\Replication Codes\Databases\DataModel.xlsx', 'ARGspreadPARAM','A1:D2001');


%Load empirical data
[data,~,~]=xlsread('P:\Research\20160628_default_volatility\Codes\Replication Codes\Databases\DataModel.xlsx', 'data','B56:G91');
spreads = data(:,6);
spread = spreads*100;
r=data(:,6);
c=data(:,1);
y=data(:,4);
nx=data(:,5);

clear gdp data

%Prepare data
r = r*100;
gdp = detrend(log(y));
c = detrend(log(c));
% nx =[nan(39,1);nx];

%-------------------------------------------
% Cyclical Properties in the data
%-------------------------------------------
%These are the moments reported in Table 9
disp('sigma_c/sigma_y')
mean(std(c)./std(gdp))

disp('sigma_ny/sigma_y')
nxy = nx./(y/100000);
mean(std(nxy)./std(gdp))

disp('rho(c,y)')
corrcoef(c,gdp)

disp('rho(nxy,y)')
corrcoef(nxy,gdp)

disp('rho(spread,y)')
corrcoef(spread(isnan(spread)==0),gdp(isnan(spread)==0))

%-------------------------------------------
% Compute the Stochastic Volatility
%-------------------------------------------

%Prepare variables
par = mean(parameters);

% this is the length of the process
err=errors';
err=err(:);
err(isnan(err)==1) = [];
size_data = size(err,1)/2000;
err=reshape(err,2000,size_data)';

begin =1;
sigma_t=zeros(size_data-1,1);
for j= 1:2000
%Compute time varying volatility
sigma_t(1,j) = par(1,4)+par(1,3)*err(1,j);
for t=2:size_data
sigma_t(t,j) = (1-par(1,1))*par(1,4) +par(1,1).*sigma_t(t-1,j)+par(1,3)*err(t,j);
end
sigma_t = (sigma_t);
sigmatj=sigma_t(:,j);

%Compute contemporaneous correlations
corr_spreads1=corrcoef(sigma_t(begin:end,j),spread(begin:end,1));
corr_c1=corrcoef(sigma_t(begin:end,j),c(begin:end,1));
corr_y1=corrcoef(sigma_t(begin:end,j),gdp(begin:end,1));
corr_nx1=corrcoef(sigma_t(begin:end,j),nx(begin:end,1));
corr_spreads(j,1)=corr_spreads1(1,2);
corr_c(j,1)=corr_c1(1,2);
corr_y(j,1)=corr_y1(1,2);
corr_nx(j,1)=corr_nx1(1,2);

%Compute correlation between volatility and income at tifferent time
%horizons
corr_y11=corrcoef(sigma_t(begin:end-1,j),gdp(begin+1:end,1));
corr_y12=corrcoef(sigma_t(begin:end-2,j),gdp(begin+2:end,1));
corr_y13=corrcoef(sigma_t(begin:end-3,j),gdp(begin+3:end,1));
corr_y14=corrcoef(sigma_t(begin:end-4,j),gdp(begin+4:end,1));
corr_y15=corrcoef(sigma_t(begin:end-5,j),gdp(begin+5:end,1));
corr_y1_1=corrcoef(sigma_t(begin+1:end,j),gdp(begin:end-1,1));
corr_y21=corrcoef(sigma_t(begin+2:end,j),gdp(begin:end-2,1));
corr_y31=corrcoef(sigma_t(begin+3:end,j),gdp(begin:end-3,1));
corr_y41=corrcoef(sigma_t(begin+4:end,j),gdp(begin:end-4,1));
corr_y51=corrcoef(sigma_t(begin+5:end,j),gdp(begin:end-5,1));


corr_y_t(j,6)=corr_y1(1,2);
corr_y_t(j,7)=corr_y11(1,2);
corr_y_t(j,8)=corr_y12(1,2);
corr_y_t(j,9)=corr_y13(1,2);
corr_y_t(j,10)=corr_y14(1,2);
corr_y_t(j,11)=corr_y15(1,2);
corr_y_t(j,5)=corr_y1_1(1,2);
corr_y_t(j,4)=corr_y21(1,2);
corr_y_t(j,3)=corr_y31(1,2);
corr_y_t(j,2)=corr_y41(1,2);
corr_y_t(j,1)=corr_y51(1,2);
end

%---------------------------------------------------------------------
%These are the moments reported in the first column of Table 2 and 5
%---------------------------------------------------------------------
disp('corr_spreads')
corr_spreads(corr_spreads>mean(corr_spreads)+2.5*std(corr_spreads)|corr_spreads<mean(corr_spreads)-2.5*std(corr_spreads))=[];
disp(mean(corr_spreads))
disp(std(corr_spreads))

disp('corr_y')
corr_y(corr_y>mean(corr_y)+2.5*std(corr_y)|corr_y<mean(corr_y)-2.5*std(corr_y))=[];
disp(mean(corr_y))
disp(std(corr_y))

disp('corr_c')
corr_c(corr_c>mean(corr_c)+2.5*std(corr_c)|corr_c<mean(corr_c)-2.5*std(corr_c))=[];
disp(mean(corr_c))
disp(std(corr_c))

disp('corr_nx')
corr_nx(corr_nx>mean(corr_nx)+2.5*std(corr_nx)|corr_nx<mean(corr_nx)-2.5*std(corr_nx))=[];
disp(mean(corr_nx))
disp(std(corr_nx))

%---------------------------------------------------------------------
%These are the correlations plotted in Figure 1 and 8 for Argentina 
%--------------------------------------------------------------------------

disp('corr_y_t')
disp(mean(corr_y_t))
disp(std(corr_y_t))

%Compute the average path of the volatility of interest rates from the
%different extractions
sigma_t_data=mean(sigma_t,2);

%Compute the average correlations
corr_spreads_fin=mean(corr_spreads);
corr_y_fin=mean(corr_y);
corr_c_fin=mean(corr_c);
corr_nx_fin=nanmean(corr_nx);

spreads1 = spread*ones(1,80);
spreads_data = spreads1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Moments for the volatility of interest rate spreads from the BASELINE MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Keep some of the data variables
clearvars -except sigma_t_data error_data spreads_data error_gdp_data sigma_t y_data corr_y_t

%Load errors for the process of the volatility of interest rate spreads
%obtained from the econometric model using simulated series from the
%baseline model
[errors,~,~]=xlsread('P:\Research\20160628_default_volatility\Codes\Replication Codes\Databases\DataModel.xlsx', 'SimulationsModelERROR');
%Load parameters for the process of the volatility of interest rate spreads
%obtained from the econometric model using simulated series from the
%baseline model
[parameters,~,~]=xlsread('P:\Research\20160628_default_volatility\Codes\Replication Codes\Databases\DataModel.xlsx', 'SimulationsModelPARAM','A1:D2001');


%Load simulated time series obtained from the model
[data,~,~]=xlsread('P:\Research\20160628_default_volatility\Codes\Replication Codes\Databases\DataModel.xlsx', 'BaselineModelSim','A2:D37');
spread = data(:,1);
y = data(:,2);
c = data(:,3);
nx = data(:,4);

%Prepare variables
gdp = log(y); %GDP process from simulations is already detrended and logged too
c = log(c);
par = mean(parameters);
clear data

%---------------------------------------------------------------------
% Compute the Stochastic Volatility
%---------------------------------------------------------------------

err=errors';
err=err(:);
err(isnan(err)==1) = [];
size_data = size(err,1)/2000;
err=reshape(err,2000,size_data)';

clear corr_spreads corr_c corr_y corr_nx
begin =1;
sigma_t=zeros(size_data,1);
for j= 1:2000
%Compute time varying volatility
sigma_t(1,j) = par(1,4)+par(1,3)*err(1,j);
for t=2:size_data
sigma_t(t,j) = (1-par(1,1))*par(1,4) +par(1,1).*sigma_t(t-1,j)+par(1,3)*err(t,j);
end
sigma_t = (sigma_t);
sigmatj=sigma_t(:,j);

%Compute contemporaneous correlations
corr_spreads1=corrcoef(sigma_t(begin:end,j),spread(begin:end-1,1));
corr_c1=corrcoef(sigma_t(begin:end,j),c(begin:end-1,1));
corr_y1=corrcoef(sigma_t(begin:end,j),gdp(begin:end-1,1));
corr_nx1=corrcoef(sigma_t(begin:end,j),nx(begin:end-1,1));
corr_spreads(j,1)=corr_spreads1(1,2);
corr_c(j,1)=corr_c1(1,2);
corr_y(j,1)=corr_y1(1,2);
corr_nx(j,1)=corr_nx1(1,2);


%Compute correlation between the volatility of interest rate spreads and
%income at different time horizons
corr_y11=corrcoef(sigma_t(begin:end-1,j),gdp(begin+1:end-1,1));
corr_y12=corrcoef(sigma_t(begin:end-2,j),gdp(begin+2:end-1,1));
corr_y13=corrcoef(sigma_t(begin:end-3,j),gdp(begin+3:end-1,1));
corr_y14=corrcoef(sigma_t(begin:end-4,j),gdp(begin+4:end-1,1));
corr_y15=corrcoef(sigma_t(begin:end-5,j),gdp(begin+5:end-1,1));
corr_y1_1=corrcoef(sigma_t(begin+1:end,j),gdp(begin:end-2,1));
corr_y21=corrcoef(sigma_t(begin+2:end,j),gdp(begin:end-3,1));
corr_y31=corrcoef(sigma_t(begin+3:end,j),gdp(begin:end-4,1));
corr_y41=corrcoef(sigma_t(begin+4:end,j),gdp(begin:end-5,1));
corr_y51=corrcoef(sigma_t(begin+5:end,j),gdp(begin:end-6,1));
corr_y_t_model(j,6)=corr_y1(1,2);
corr_y_t_model(j,7)=corr_y11(1,2);
corr_y_t_model(j,8)=corr_y12(1,2);
corr_y_t_model(j,9)=corr_y13(1,2);
corr_y_t_model(j,10)=corr_y14(1,2);
corr_y_t_model(j,11)=corr_y15(1,2);
corr_y_t_model(j,5)=corr_y1_1(1,2);
corr_y_t_model(j,4)=corr_y21(1,2);
corr_y_t_model(j,3)=corr_y31(1,2);
corr_y_t_model(j,2)=corr_y41(1,2);
corr_y_t_model(j,1)=corr_y51(1,2);
end

%--------------------------------------------------------------
%These are the moments reported in the second column of Table 5
%--------------------------------------------------------------
disp('corr_spreads')
corr_spreads(corr_spreads>mean(corr_spreads)+2.5*std(corr_spreads)|corr_spreads<mean(corr_spreads)-2.5*std(corr_spreads))=[];
disp(mean(corr_spreads))
disp(std(corr_spreads))

disp('corr_y')
corr_y(corr_y>mean(corr_y)+2.5*std(corr_y)|corr_y<mean(corr_y)-2.5*std(corr_y))=[];
disp(mean(corr_y))
disp(std(corr_y))

disp('corr_c')
corr_c(corr_c>mean(corr_c)+2.5*std(corr_c)|corr_c<mean(corr_c)-2.5*std(corr_c))=[];
disp(mean(corr_c))
disp(std(corr_c))

disp('corr_nx')
corr_nx(corr_nx>mean(corr_nx)+2.5*std(corr_nx)|corr_nx<mean(corr_nx)-2.5*std(corr_nx))=[];
disp(mean(corr_nx))
disp(std(corr_nx))

%--------------------------------------------------------------
%These are the correlations plotted in Figure  8 for Argentina 
%--------------------------------------------------------------------------
disp('corr_y_t_model')
disp(mean(corr_y_t_model))
disp(std(corr_y_t_model))

corr_spreads_fin=mean(corr_spreads);
corr_c_fin=mean(corr_c);
corr_y_fin=mean(corr_y);
corr_nx_fin=mean(corr_nx);

%Compute the average path of interest rate volatility from simulations
sigma_t_moments=mean(sigma_t,2);






%% Figure 8 in the paper
% Plot Figure for the correlations at different leads and lag in Argentina
%----------------------------------------
close all
figure;
hFig = figure;
date = [-5:1:5]';
plot(date, mean(corr_y_t,1),'b','LineWidth',1.5)
hold on
% plot(date, corr_tfp_sigma_data_CIP,'--b','LineWidth',1)
% plot(date, corr_tfp_sigma_data_CIM,'--b','LineWidth',1)
plot(date, mean(corr_y_t_model,1),'-.r','LineWidth',1.5)
xlim([-5.25 5.25])
xlabel('Leads and Lags','FontSize',14);
ylabel('\rho (vol_{t},y_{t+lead})','FontSize',14);
legend  ({'Data', 'Model'},'Location','Best','FontSize',14)
legend BOXOFF
set(gcf, 'PaperPosition', [0 0 8 6.5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [8 6.5]); %Set the paper to have width 5 and height 5.
set(gca,'XTick',-5:1:5);
print(hFig,'-dpdf','P:\Research\20160628_default_volatility\Codes\Replication Codes\Output\tcorrelationJIE.pdf','-opengl')




%% Figure 5 in the paper 
% Scatter plot for the correlation between the level of spreads and their volatility
%--------------------------------------------------------------------------
close all
figure;
hFig = figure;
box on;
spreads1 = spread*ones(1,80);
scatter(spreads_data(:,10),sigma_t_data,'+','r')
% lsline
% k = lsline;
% set(k(1),'color','b')
xlabel('Spreads','FontSize',14);
ylabel('Spreads Volatility','FontSize',14);
% title('Volatility and Spreads` Level')
hold on
scatter(spreads1(1:end-1,10),mean(sigma_t(1:end,:),2)-0.1,'b')
h = lsline;
set(h(1),'color','b','linewidth',2)
set(h(2),'color','r','linewidth',2,'Linestyle','--')
% xlim([200 1800])
% ylim([-6.05 -4])
legend  ({'Data','Model'},'Location','Best','FontSize',14)
legend BOXOFF
set(gcf, 'PaperPosition', [0 0 8 4.5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [8 4.5]); %Set the paper to have width 5 and height 5.
print(hFig,'-dpdf','P:\Research\20160628_default_volatility\Codes\Replication Codes\Output\scatter1993JIE.pdf','-opengl')


%% %Figure 4 in the paper 
%Plots the time-series of the spreads 1993Q1-2001Q3
%---------------------------------------------------


clearvars -except sigma_t_data sigma_t_moments  spreads_data spread

[spreadsimul,~,~]=xlsread('P:\Research\20160628_default_volatility\Codes\Replication Codes\Databases\DataModel.xlsx', 'AlternativeModelSim','B9:C43');


spreddata=spreads_data(:,1);
spreddata = [spreddata; 2860 ];
spredarelllt=spread;
spredarellfine=spreadsimul(:,1);
spredquadlt=spreadsimul(:,2);


close all
figure(1);
spredarellfine(isnan(spredarellfine)==1)=inf;
hFig = figure(1);
date = [1993:0.25:2002]';
plot(date, spreddata(1:end,1),'k','LineWidth',2.5)
hold on
% plot(date, spredarell,'--b','LineWidth',1.5)
plot(date(1:end-1,1), spredarelllt(1:end,1),'--bd','LineWidth',1.5)
% plot(date, spredarelllt,'-d','LineWidth',1.5)
% plot(date, spredquadlt,'--*','LineWidth',1.5)
xlabel('Years','FontSize',14);
ylabel('Spreads','FontSize',14);
legend  ({'Data','Model'},'Location','Best','FontSize',14)
legend BOXOFF
set(gcf, 'PaperPosition', [0 0 8 6.5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [8 6.5]); %Set the paper to have width 5 and height 5.
print(hFig,'-dpdf','P:\Research\20160628_default_volatility\Codes\Replication Codes\Output\spreadsseriesJIE.pdf','-opengl')

%% Figure 6 in the paper 
%Plots the time-series of spreads in different versions of the model
%--------------------------------------------------------------------------
close all
figure(2);
hFig = figure(2);
date = [1993:0.25:2001.75]';
plot(date, spreddata(1:end-1,1),'k','LineWidth',2.5)
hold on
% plot(date, spredarell(1:end-1,1),'--b','LineWidth',1.5)
plot(date(1:end-1), spredarellfine(1:end,1),':r*','LineWidth',1)
plot(date(1:end-1), spredquadlt(1:end,1),'-.gs','LineWidth',1)
xlabel('Years','FontSize',14);
ylabel('Spreads','FontSize',14);
legend  ({'Data', 'ST-Bonds', 'LT-Bonds & Quad. Costs'},'Location','Best','FontSize',14)
legend BOXOFF
set(gcf, 'PaperPosition', [0 0 8 6.5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [8 6.5]); %Set the paper to have width 5 and height 5.
print(hFig,'-dpdf','P:\Research\20160628_default_volatility\Codes\Replication Codes\Output\spreadsseries1JIE.pdf','-opengl')


%% Figure 7 in the paper
%Plots the time-series of the volatility *of interest rate spreads in different versions of the model
%----------------------------------------

clear all
[parameters,~,~]=xlsread('P:\Research\20160628_default_volatility\Codes\Replication Codes\Databases\DataModel.xlsx', 'paramSummary','B2:E5');
[errors,~,~]=xlsread('P:\Research\20160628_default_volatility\Codes\Replication Codes\Databases\DataModel.xlsx', 'errorsummary','A1:L24000');

    errdata=errors(:,1:3)';
    errdata=errdata(:);
    errdata(isnan(errdata)==1) = [];
    size_data = size(errdata,1)/2000
    errdata=reshape(errdata,2000,size_data)';

    errarellfine=errors(:,4:6)';
    errarellfine=errarellfine(:);
    errarellfine(isnan(errarellfine)==1) = [];
    size_arellfine = size(errarellfine,1)/2000
    errarellfine=reshape(errarellfine,2000,size_arellfine)';

    errarelllt=errors(:,7:9)';
    errarelllt=errarelllt(:);
    errarelllt(isnan(errarelllt)==1) = [];
    size_arelllt = size(errarelllt,1)/2000
    errarelllt=reshape(errarelllt,2000,size_arelllt)';

    errquadlt=errors(:,7:9)';
    errquadlt=errquadlt(:);
    errquadlt(isnan(errquadlt)==1) = [];
    size_quadlt = size(errquadlt,1)/2000
    errquadlt=reshape(errquadlt,2000,size_quadlt)';

    
begin =1;
sigma_t_data=zeros(size_data,1);
%Compute time varying spreads
for j= 1:2000
sigma_t_data(1,j) = parameters(1,4)+parameters(1,3)*errdata(1,j);
for t=2:size_data
sigma_t_data(t,j) = (1-parameters(1,1))*parameters(1,4) +parameters(1,1).*sigma_t_data(t-1,j)+parameters(1,3)*errdata(t,j);
end
end
sigma_t_data=mean(sigma_t_data,2);
sigma_t_data(end,:) =[];

begin =1;
sigma_t_arellfine=zeros(size_arellfine,1);
%Compute time varying spreads
for j= 1:2000
sigma_t_arellfine(1,j) = parameters(2,4)+parameters(2,3)*errarellfine(1,j);
for t=2:size_arellfine
sigma_t_arellfine(t,j) = (1-parameters(2,1))*parameters(2,4) +parameters(2,1).*sigma_t_arellfine(t-1,j)+parameters(2,3)*errarellfine(t,j);
end
end
sigma_t_arellfine=mean(sigma_t_arellfine,2);
% sigma_t_arellfine = sigma_t_arellfine(2:end,1);

begin =1;
sigma_t_arelllt=zeros(size_arelllt,1);
%Compute time varying spreads
for j= 1:2000
sigma_t_arelllt(1,j) = parameters(3,4)+parameters(3,3)*errarelllt(1,j);
for t=2:size_arelllt
sigma_t_arelllt(t,j) = (1-parameters(3,1))*parameters(3,4) +parameters(3,1).*sigma_t_arelllt(t-1,j)+parameters(3,3)*errarelllt(t,j);
end
end
sigma_t_arelllt=mean(sigma_t_arelllt,2);
% sigma_t_arelllt = sigma_t_arelllt(2:end,1);

   
begin =1;
sigma_t_quadlt=zeros(size_quadlt,1);
%Compute time varying spreads
for j= 1:2000
sigma_t_quadlt(1,j) = parameters(4,4)+parameters(4,3)*errquadlt(1,j);
for t=2:size_quadlt
sigma_t_quadlt(t,j) = (1-parameters(4,1))*parameters(4,4) +parameters(4,1).*sigma_t_quadlt(t-1,j)+parameters(4,3)*errquadlt(t,j);
end
end
sigma_t_quadlt=mean(sigma_t_quadlt,2);
% sigma_t_quadlt = sigma_t_quadlt(2:end,1);

% sigma_t_arellfine(end+1:end+2,:)=NaN;
sigma_t_arelllt(end+1,:)=NaN;
sigma_t_quadlt(end+1,:)=NaN;

close all
figure(1);
hFig = figure(1);
date = [1993:0.25:2001.5]';
plot(date, sigma_t_data,'k','LineWidth',2.5)
hold on
% plot(date, spredarell,'--b','LineWidth',1.5)
plot(date, sigma_t_arelllt(1:end-1,1),'--bd','LineWidth',1.5)
plot(date, sigma_t_arellfine,':r*','LineWidth',1)
plot(date, sigma_t_quadlt(1:end-1,1),'-.gs','LineWidth',1)
xlabel('Years','FontSize',14);
ylabel('Time Varying Volatility of Spreads','FontSize',14);
legend  ({'Data', 'LT-Bonds', 'ST-Bonds', 'LT-Bonds & Quad. Costs'},'Location','North','FontSize',14)
legend BOXOFF
set(gcf, 'PaperPosition', [0 0 8 6.5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [8 6.5]); %Set the paper to have width 5 and height 5.
print(hFig,'-dpdf','P:\Research\20160628_default_volatility\Codes\Replication Codes\Output\spreadsvol1JIE.pdf','-opengl')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figures and moments for Latin American countries


%----------------------------------------
%Brazil
%----------------------------------------

clear all
%Load Data obtained 
[errors,~,~]=xlsread('P:\Research\20160628_default_volatility\Codes\Replication Codes\Databases\DataModel.xlsx', 'BRAspreadERRORS');
[parameters,~,~]=xlsread('P:\Research\20160628_default_volatility\Codes\Replication Codes\Databases\DataModel.xlsx', 'BRAspreadPARAM','A1:D2001');
[data,~,~]=xlsread('P:\Research\20160628_default_volatility\Codes\Replication Codes\Databases\DataModel.xlsx', 'Brazil','B2:E56');
spreads = data(:,4);
c = data(:,2);
y = data(:,1);
nx = data(:,3);
clear gdp
y_data=y;
gdp = detrend(log(y));
c = detrend(log(c));

%Prepare variables
par = mean(parameters);

 % this is the length of the process
err=errors';
err=err(:);
err(isnan(err)==1) = [];
size_data = size(err,1)/2000;
err=reshape(err,2000,size_data)';
% err= median(err,2);

begin =1;
sigma_t=zeros(size_data-1,1);
%Compute time varying spreads
for j= 1:2000
%Compute volatility
sigma_t(1,j) = par(1,4)+par(1,3)*err(1,j);
for t=2:size_data
sigma_t(t,j) = (1-par(1,1))*par(1,4) +par(1,1).*sigma_t(t-1,j)+par(1,3)*err(t,j);
end

%Compute correlations 
sigma_t = (sigma_t);
sigmatj=sigma_t(:,j);
corr_spreads1=corrcoef(sigma_t(begin:end,j),spreads(begin:end,1));
corr_c1=corrcoef(sigma_t(begin:end,j),c(begin:end,1));
corr_y1=corrcoef(sigma_t(begin:end,j),gdp(begin:end,1));
corr_y11=corrcoef(sigma_t(begin:end-1,j),gdp(begin+1:end,1));
corr_y12=corrcoef(sigma_t(begin:end-2,j),gdp(begin+2:end,1));
corr_y13=corrcoef(sigma_t(begin:end-3,j),gdp(begin+3:end,1));
corr_y14=corrcoef(sigma_t(begin:end-4,j),gdp(begin+4:end,1));
corr_y15=corrcoef(sigma_t(begin:end-5,j),gdp(begin+5:end,1));
corr_y1_1=corrcoef(sigma_t(begin+1:end,j),gdp(begin:end-1,1));
corr_y21=corrcoef(sigma_t(begin+2:end,j),gdp(begin:end-2,1));
corr_y31=corrcoef(sigma_t(begin+3:end,j),gdp(begin:end-3,1));
corr_y41=corrcoef(sigma_t(begin+4:end,j),gdp(begin:end-4,1));
corr_y51=corrcoef(sigma_t(begin+5:end,j),gdp(begin:end-5,1));

corr_nx1=corrcoef(sigma_t(begin:end,j),nx(begin:end,1));

corr_spreads(j,1)=corr_spreads1(1,2);
corr_c(j,1)=corr_c1(1,2);
corr_y(j,1)=corr_y1(1,2);
corr_y_t(j,6)=corr_y1(1,2);
corr_y_t(j,7)=corr_y11(1,2);
corr_y_t(j,8)=corr_y12(1,2);
corr_y_t(j,9)=corr_y13(1,2);
corr_y_t(j,10)=corr_y14(1,2);
corr_y_t(j,11)=corr_y15(1,2);
corr_y_t(j,5)=corr_y1_1(1,2);
corr_y_t(j,4)=corr_y21(1,2);
corr_y_t(j,3)=corr_y31(1,2);
corr_y_t(j,2)=corr_y41(1,2);
corr_y_t(j,1)=corr_y51(1,2);
corr_nx(j,1)=corr_nx1(1,2);
end
%contemporaneous correlations reported in Table 2
disp('corr_spreads')
corr_spreads(corr_spreads>mean(corr_spreads)+2.5*std(corr_spreads)|corr_spreads<mean(corr_spreads)-2.5*std(corr_spreads))=[];
disp(mean(corr_spreads))
disp(std(corr_spreads))

disp('corr_y')
corr_y(corr_y>mean(corr_y)+2.5*std(corr_y)|corr_y<mean(corr_y)-2.5*std(corr_y))=[];
disp(mean(corr_y))
disp(std(corr_y))

disp('corr_c')
corr_c(corr_c>mean(corr_c)+2.5*std(corr_c)|corr_c<mean(corr_c)-2.5*std(corr_c))=[];
disp(mean(corr_c))
disp(std(corr_c))

disp('corr_nx')
corr_nx(corr_nx>mean(corr_nx)+2.5*std(corr_nx)|corr_nx<mean(corr_nx)-2.5*std(corr_nx))=[];
disp(mean(corr_nx))
disp(std(corr_nx))


gdp_data_bz = detrend(log(y_data));
sigma_t_data_bz=mean(sigma_t,2);

%Correaltion between volatility and income at different time horizons
corr_tfp_sigma_data_bz =mean(corr_y_t);


%----------------------------------------
%Ecuador
%----------------------------------------
clearvars -except  corr_tfp_sigma_data_bz

%Load Data obtained 
[errors,~,~]=xlsread('P:\Research\20160628_default_volatility\Codes\Replication Codes\Databases\DataModel.xlsx', 'ECUspreadERRORS');
[parameters,~,~]=xlsread('P:\Research\20160628_default_volatility\Codes\Replication Codes\Databases\DataModel.xlsx', 'ECUspreadPARAM','A1:D2001');
[data,~,~]=xlsread('P:\Research\20160628_default_volatility\Codes\Replication Codes\Databases\DataModel.xlsx', 'Ecuador','B2:E56');
spreads = data(:,4);
c = data(:,2);
y = data(:,1);
nx = data(:,3);
clear gdp
y_data=y;
gdp = detrend(log(y));
c = detrend(log(c));

%Prepare variables
par = mean(parameters);

 % this is the length of the process

err=errors';
err=err(:);
err(isnan(err)==1) = [];
size_data = size(err,1)/2000;
err=reshape(err,2000,size_data)';
% err= median(err,2);

clear corr_spreads corr_c corr_y corr_nx
begin =1;
sigma_t=zeros(size_data-1,1);
%Compute time varying spreads
for j= 1:2000
sigma_t(1,j) = par(1,4)+par(1,3)*err(1,j);
for t=2:size_data
sigma_t(t,j) = (1-par(1,1))*par(1,4) +par(1,1).*sigma_t(t-1,j)+par(1,3)*err(t,j);
end
sigma_t = (sigma_t);
sigmatj=sigma_t(:,j);
corr_spreads1=corrcoef(sigma_t(begin:end,j),spreads(begin:end,1));
corr_c1=corrcoef(sigma_t(begin:end,j),c(begin:end,1));
corr_y1=corrcoef(sigma_t(begin:end,j),gdp(begin:end,1));
corr_y11=corrcoef(sigma_t(begin:end-1,j),gdp(begin+1:end,1));
corr_y12=corrcoef(sigma_t(begin:end-2,j),gdp(begin+2:end,1));
corr_y13=corrcoef(sigma_t(begin:end-3,j),gdp(begin+3:end,1));
corr_y14=corrcoef(sigma_t(begin:end-4,j),gdp(begin+4:end,1));
corr_y15=corrcoef(sigma_t(begin:end-5,j),gdp(begin+5:end,1));
corr_y1_1=corrcoef(sigma_t(begin+1:end,j),gdp(begin:end-1,1));
corr_y21=corrcoef(sigma_t(begin+2:end,j),gdp(begin:end-2,1));
corr_y31=corrcoef(sigma_t(begin+3:end,j),gdp(begin:end-3,1));
corr_y41=corrcoef(sigma_t(begin+4:end,j),gdp(begin:end-4,1));
corr_y51=corrcoef(sigma_t(begin+5:end,j),gdp(begin:end-5,1));

corr_nx1=corrcoef(sigma_t(begin:end,j),nx(begin:end,1));

corr_spreads(j,1)=corr_spreads1(1,2);
corr_c(j,1)=corr_c1(1,2);
corr_y(j,1)=corr_y1(1,2);
corr_y_t(j,6)=corr_y1(1,2);
corr_y_t(j,7)=corr_y11(1,2);
corr_y_t(j,8)=corr_y12(1,2);
corr_y_t(j,9)=corr_y13(1,2);
corr_y_t(j,10)=corr_y14(1,2);
corr_y_t(j,11)=corr_y15(1,2);
corr_y_t(j,5)=corr_y1_1(1,2);
corr_y_t(j,4)=corr_y21(1,2);
corr_y_t(j,3)=corr_y31(1,2);
corr_y_t(j,2)=corr_y41(1,2);
corr_y_t(j,1)=corr_y51(1,2);
corr_nx(j,1)=corr_nx1(1,2);
end

disp('corr_spreads')
corr_spreads(corr_spreads>mean(corr_spreads)+2.5*std(corr_spreads)|corr_spreads<mean(corr_spreads)-2.5*std(corr_spreads))=[];
disp(mean(corr_spreads))
disp(std(corr_spreads))

disp('corr_y')
corr_y(corr_y>mean(corr_y)+2.5*std(corr_y)|corr_y<mean(corr_y)-2.5*std(corr_y))=[];
disp(mean(corr_y))
disp(std(corr_y))

disp('corr_c')
corr_c(corr_c>mean(corr_c)+2.5*std(corr_c)|corr_c<mean(corr_c)-2.5*std(corr_c))=[];
disp(mean(corr_c))
disp(std(corr_c))

disp('corr_nx')
corr_nx(corr_nx>mean(corr_nx)+2.5*std(corr_nx)|corr_nx<mean(corr_nx)-2.5*std(corr_nx))=[];
disp(mean(corr_nx))
disp(std(corr_nx))

corr_tfp_sigma_data_ec =mean(corr_y_t);



%----------------------------------------
%Venezuela
%----------------------------------------

clearvars -except  corr_tfp_sigma_data_bz corr_tfp_sigma_data_ec
%Load Data Obtained from fortran
[errors,~,~]=xlsread('P:\Research\20160628_default_volatility\Codes\Replication Codes\Databases\DataModel.xlsx', 'VNZspreadERRORS');
[parameters,~,~]=xlsread('P:\Research\20160628_default_volatility\Codes\Replication Codes\Databases\DataModel.xlsx', 'VNZspreadPARAM','A1:D2001');
[data,~,~]=xlsread('P:\Research\20160628_default_volatility\Codes\Replication Codes\Databases\DataModel.xlsx', 'Venezuela','B2:E58');
spreads = data(:,4);
c = data(:,2);
y = data(:,1);
nx = data(:,3);
clear gdp
y_data=y;
gdp = detrend(log(y(14:end,1)));
c = detrend(log(c(18:end,1)));

%Prepare variables
par = mean(parameters);

err=errors';
err=err(:);
err(isnan(err)==1) = [];
size_data = size(err,1)/2000;
err=reshape(err,2000,size_data)';
% err= median(err,2);

clear corr_spreads corr_c corr_y corr_nx
begin =1;
sigma_t=zeros(size_data-1,1);
%Compute time varying spreads
for j= 1:2000
sigma_t(1,j) = par(1,4)+par(1,3)*err(1,j);
for t=2:size_data
sigma_t(t,j) = (1-par(1,1))*par(1,4) +par(1,1).*sigma_t(t-1,j)+par(1,3)*err(t,j);
end
sigma_t = (sigma_t);
sigmatj=sigma_t(:,j);
corr_spreads1=corrcoef(sigma_t(begin:end,j),spreads(begin:end,1));
corr_c1=corrcoef(sigma_t(18:end,j),c(begin:end,1));
corr_y1=corrcoef(sigma_t(14:end,j),gdp(begin:end,1));
corr_y11=corrcoef(sigma_t(14:end-1,j),gdp(begin+1:end,1));
corr_y12=corrcoef(sigma_t(14:end-2,j),gdp(begin+2:end,1));
corr_y13=corrcoef(sigma_t(14:end-3,j),gdp(begin+3:end,1));
corr_y14=corrcoef(sigma_t(14:end-4,j),gdp(begin+4:end,1));
corr_y15=corrcoef(sigma_t(14:end-5,j),gdp(begin+5:end,1));
corr_y1_1=corrcoef(sigma_t(14+1:end,j),gdp(begin:end-1,1));
corr_y21=corrcoef(sigma_t(14+2:end,j),gdp(begin:end-2,1));
corr_y31=corrcoef(sigma_t(14+3:end,j),gdp(begin:end-3,1));
corr_y41=corrcoef(sigma_t(14+4:end,j),gdp(begin:end-4,1));
corr_y51=corrcoef(sigma_t(14+5:end,j),gdp(begin:end-5,1));

corr_nx1=corrcoef(sigma_t(2:end,j),nx(begin+1:end,1));

corr_spreads(j,1)=corr_spreads1(1,2);
corr_c(j,1)=corr_c1(1,2);
corr_y(j,1)=corr_y1(1,2);
corr_y_t(j,6)=corr_y1(1,2);
corr_y_t(j,7)=corr_y11(1,2);
corr_y_t(j,8)=corr_y12(1,2);
corr_y_t(j,9)=corr_y13(1,2);
corr_y_t(j,10)=corr_y14(1,2);
corr_y_t(j,11)=corr_y15(1,2);
corr_y_t(j,5)=corr_y1_1(1,2);
corr_y_t(j,4)=corr_y21(1,2);
corr_y_t(j,3)=corr_y31(1,2);
corr_y_t(j,2)=corr_y41(1,2);
corr_y_t(j,1)=corr_y51(1,2);
corr_nx(j,1)=corr_nx1(1,2);
end

disp('corr_spreads')
corr_spreads(corr_spreads>mean(corr_spreads)+2.5*std(corr_spreads)|corr_spreads<mean(corr_spreads)-2.5*std(corr_spreads))=[];
disp(mean(corr_spreads))
disp(std(corr_spreads))

disp('corr_y')
corr_y(corr_y>mean(corr_y)+2.5*std(corr_y)|corr_y<mean(corr_y)-2.5*std(corr_y))=[];
disp(mean(corr_y))
disp(std(corr_y))

disp('corr_c')
corr_c(corr_c>mean(corr_c)+2.5*std(corr_c)|corr_c<mean(corr_c)-2.5*std(corr_c))=[];
disp(mean(corr_c))
disp(std(corr_c))

disp('corr_nx')
corr_nx(corr_nx>mean(corr_nx)+2.5*std(corr_nx)|corr_nx<mean(corr_nx)-2.5*std(corr_nx))=[];
disp(nanmean(corr_nx))
disp(std(corr_nx))

corr_tfp_sigma_data_ve =mean(corr_y_t);


%----------------------------------------
%Argentina
%----------------------------------------

clearvars -except  corr_tfp_sigma_data_bz corr_tfp_sigma_data_ec  corr_tfp_sigma_data_ve
%Load Data Obtained from fortran
[errors,~,~]=xlsread('P:\Research\20160628_default_volatility\Codes\Replication Codes\Databases\DataModel.xlsx', 'ARGspreadERRORS');
[parameters,~,~]=xlsread('P:\Research\20160628_default_volatility\Codes\Replication Codes\Databases\DataModel.xlsx', 'ARGspreadPARAM','A1:D2001');
[spreads,~,~]=xlsread('P:\Research\20160628_default_volatility\Codes\Replication Codes\Databases\DataModel.xlsx', 'data','B56:G91');
%Load empirical data
[data,~,~]=xlsread('P:\Research\20160628_default_volatility\Codes\Replication Codes\Databases\DataModel.xlsx', 'data','B56:G91');
spreads = data(:,6);
spread = spreads*100;
r=data(:,6);
c=data(:,1);
y=data(:,4);
nx=data(:,5);
gdp = detrend(log(y));
c = detrend(log(c));
%Prepare variables
par = mean(parameters);

 % this is the length of the process
err=errors';
err=err(:);
err(isnan(err)==1) = [];
size_data = size(err,1)/2000;
err=reshape(err,2000,size_data)';
% err= median(err,2);

clear corr_spreads corr_c corr_y corr_nx corr_y_t
begin =1;
sigma_t=zeros(size_data-1,1);
%Compute time varying spreads
for j= 1:2000
sigma_t(1,j) = par(1,4)+par(1,3)*err(1,j);
for t=2:size_data
sigma_t(t,j) = (1-par(1,1))*par(1,4) +par(1,1).*sigma_t(t-1,j)+par(1,3)*err(t,j);
end
sigma_t = (sigma_t);
sigmatj=sigma_t(:,j);
corr_spreads1=corrcoef(sigma_t(begin:end,j),spread(begin:end,1));
corr_c1=corrcoef(sigma_t(begin:end,j),c(begin:end,1));
corr_y1=corrcoef(sigma_t(begin:end,j),gdp(begin:end,1));
corr_y11=corrcoef(sigma_t(begin:end-1,j),gdp(begin+1:end,1));
corr_y12=corrcoef(sigma_t(begin:end-2,j),gdp(begin+2:end,1));
corr_y13=corrcoef(sigma_t(begin:end-3,j),gdp(begin+3:end,1));
corr_y14=corrcoef(sigma_t(begin:end-4,j),gdp(begin+4:end,1));
corr_y15=corrcoef(sigma_t(begin:end-5,j),gdp(begin+5:end,1));
corr_y1_1=corrcoef(sigma_t(begin+1:end,j),gdp(begin:end-1,1));
corr_y21=corrcoef(sigma_t(begin+2:end,j),gdp(begin:end-2,1));
corr_y31=corrcoef(sigma_t(begin+3:end,j),gdp(begin:end-3,1));
corr_y41=corrcoef(sigma_t(begin+4:end,j),gdp(begin:end-4,1));
corr_y51=corrcoef(sigma_t(begin+5:end,j),gdp(begin:end-5,1));

corr_nx1=corrcoef(sigma_t(begin:end,j),nx(begin:end,1));

corr_spreads(j,1)=corr_spreads1(1,2);
corr_c(j,1)=corr_c1(1,2);
corr_y(j,1)=corr_y1(1,2);
corr_y_t(j,6)=corr_y1(1,2);
corr_y_t(j,7)=corr_y11(1,2);
corr_y_t(j,8)=corr_y12(1,2);
corr_y_t(j,9)=corr_y13(1,2);
corr_y_t(j,10)=corr_y14(1,2);
corr_y_t(j,11)=corr_y15(1,2);
corr_y_t(j,5)=corr_y1_1(1,2);
corr_y_t(j,4)=corr_y21(1,2);
corr_y_t(j,3)=corr_y31(1,2);
corr_y_t(j,2)=corr_y41(1,2);
corr_y_t(j,1)=corr_y51(1,2);
corr_nx(j,1)=corr_nx1(1,2);
end

disp('corr_spreads')
corr_spreads(corr_spreads>mean(corr_spreads)+2.5*std(corr_spreads)|corr_spreads<mean(corr_spreads)-2.5*std(corr_spreads))=[];
disp(mean(corr_spreads))
disp(std(corr_spreads))

disp('corr_y')
corr_y(corr_y>mean(corr_y)+2.5*std(corr_y)|corr_y<mean(corr_y)-2.5*std(corr_y))=[];
disp(mean(corr_y))
disp(std(corr_y))

disp('corr_c')
corr_c(corr_c>mean(corr_c)+2.5*std(corr_c)|corr_c<mean(corr_c)-2.5*std(corr_c))=[];
disp(mean(corr_c))
disp(std(corr_c))

disp('corr_nx')
corr_nx(corr_nx>mean(corr_nx)+2.5*std(corr_nx)|corr_nx<mean(corr_nx)-2.5*std(corr_nx))=[];
disp(mean(corr_nx))
disp(std(corr_nx))

disp('corr_y_t')
% corr_nx(corr_y_t>mean(corr_y_t)+2.5*std(corr_y_t)|corr_y_t<mean(corr_y_t)-2.5*std(corr_y_t))=[];
disp(mean(corr_y_t))
disp(std(corr_y_t))

corr_tfp_sigma_data_ar =mean(corr_y_t);


%-------------------------------------------------------------------
%Figure 1
%Correlation between volatility and income at different time horizons
%-------------------------------------------------------------------

close all
figure;
hFig = figure;
date = [-5:1:5]';
plot(date, corr_tfp_sigma_data_ar,':bs','LineWidth',1.5)
hold on
plot(date, corr_tfp_sigma_data_bz,'-ro','LineWidth',1.5)
hold on
plot(date, corr_tfp_sigma_data_ec,'--k^','LineWidth',1.5)
hold on
plot(date, corr_tfp_sigma_data_ve,'-g*','LineWidth',1.5)
xlim([-5.25 5.25])
ylim([-0.2 0.1])
xlabel('Leads and Lags','FontSize',14);
ylabel('\rho (vol_{t},y_{t+lead})','FontSize',14);
legend  ({'Argentina','Brazil', 'Ecuador', 'Venezuela'},'Location','Best','FontSize',14)
legend BOXOFF
set(gcf, 'PaperPosition', [0 0 8 6.5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [8 6.5]); %Set the paper to have width 5 and height 5.
set(gca,'XTick',-5:1:5);
print(hFig,'-dpdf','P:\Research\20160628_default_volatility\Codes\Replication Codes\Output\tcorrelationJIEothercounties.pdf','-opengl')



%% Conditional interest rate spreads


%-------------------------------------------------------------------
%Figure 2
%Conditional interest rate spreads-Baseline Model
%-------------------------------------------------------------------
clear all
load('P:\Research\20160628_default_volatility\Codes\Replication Codes\Databases\longjie_replication.mat')
i_b_g_mean =84;
i_b_g_high = 89;

close all
figure;
hFig=figure;
hold on;
box on;
plot(y_vec_3sh,10000*(1*(((q_g_pf(:,i_b_g_zero-30)).^-1-delta)-(mu_r))),'linewidth',2,'Linestyle',':')
plot(y_vec_3sh,10000*(1*(((q_g_pf(:,i_b_g_mean)).^-1-delta)-(mu_r))),'k','linewidth',2,'Linestyle','-.')
plot(y_vec_3sh,10000*(1*(((q_g_pf(:,i_b_g_high)).^-1-delta)-(mu_r))),'r','linewidth',2)

% axis([.85 1.17 -100 1500])
xlab = xlabel('Endowment','FontSize',15);
ylab = ylabel('Spread (bp)','FontSize',15);
% tit = title(Spreads policy function');

legend({'No Debt', 'Mean debt','High debt'},'FontSize',15);
legend BOXOFF
set(gcf, 'PaperPosition', [0 0 8 6.5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [8 6.5]); %Set the paper to have width 5 and height 5.
% set(findall(gcf,'type','text'),'FontSize',14,'fontWeight','bold') 
print(hFig,'-dpdf','P:\Research\20160628_default_volatility\Codes\Replication Codes\Output\policy1993jie.pdf','-opengl')

%-------------------------------------------------------------------
%Figure 9
%Conditional interest rate spreads-ST debt Model 
%-------------------------------------------------------------------
clear all
load('P:\Research\20160628_default_volatility\Codes\Replication Codes\Databases\shortjie_replication.mat')
i_b_g_mean =141;
i_b_g_high = 269;

close all
figure;
hFig=figure;
hold on;
box on;
plot(y_vec_3sh,10000*(1*(((q_g_pf(:,i_b_g_zero-30)).^-4-1)-(mu_r))),'Color',[1 0 0],'linewidth',2,'Linestyle',':')
plot(y_vec_3sh,10000*(1*(((q_g_pf(:,i_b_g_mean)).^-4-1)-(mu_r))),'Color',[0.8 0.3 0.0],'linewidth',2,'Linestyle','-.')
plot(y_vec_3sh,10000*(1*(((q_g_pf(:,i_b_g_high)).^-4-1)-(mu_r))),'Color',[0.6 0.0 0.3],'linewidth',2)

% axis([.85 1.17 -100 1500])
xlab = xlabel('Endowment','FontSize',15);
ylab = ylabel('Spread (bp)','FontSize',15);
% tit = title(Spreads policy function');

legend({'No Debt', 'Mean debt','High debt'},'FontSize',15);
legend BOXOFF
set(gcf, 'PaperPosition', [0 0 8 6.5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [8 6.5]); %Set the paper to have width 5 and height 5.
print(hFig,'-dpdf','P:\Research\20160628_default_volatility\Codes\Replication Codes\Output\policy1993arelljie.pdf','-opengl')



%-------------------------------------------------------------------
%Figure 10
%Conditional interest rate spreads-LT debt Model with quadratic costs
%-------------------------------------------------------------------
clear all
load('P:\Research\20160628_default_volatility\Codes\Replication Codes\Databases\long_quadjie_replication.mat')
i_b_g_mean =64;
i_b_g_high = 72;

close all
figure;
hFig=figure;
hold on;
box on;
plot(y_vec_3sh,10000*(1*(((q_g_pf(:,i_b_g_zero-30)).^-1-delta)-(mu_r))),'Color',[0 1 0],'linewidth',2,'Linestyle',':')
plot(y_vec_3sh,10000*(1*(((q_g_pf(:,i_b_g_mean)).^-1-delta)-(mu_r))),'Color',[0.3 0.8 0.0],'linewidth',2,'Linestyle','-.')
plot(y_vec_3sh,10000*(1*(((q_g_pf(:,i_b_g_high)).^-1-delta)-(mu_r))),'Color',[0 0.6 0.3],'linewidth',2)

% axis([.85 1.17 -100 1500])
xlab = xlabel('Endowment','FontSize',15);
ylab = ylabel('Spread (bp)','FontSize',15);
% tit = title(Spreads policy function');

legend({'No Debt', 'Mean debt','High debt'},'FontSize',15);
legend BOXOFF
set(gcf, 'PaperPosition', [0 0 8 6.5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [8 6.5]); %Set the paper to have width 5 and height 5.
print(hFig,'-dpdf','P:\Research\20160628_default_volatility\Codes\Replication Codes\Output\policy1993quadjie.pdf','-opengl')


%% Non-Linearity: Debt Issuance Policy Function, Unconditional and Conditional Debt Pricing Function

%-------------------------------------------------------------------
%Figure 3
%-------------------------------------------------------------------
% Non-Linear policy function for interest rate

load('P:\Research\20160628_default_volatility\Codes\Replication Codes\Databases\longjie_replication.mat')



mu_r    = .017;%.001;

set(0,'defaultTextInterpreter','latex');
set(0,'defaultAxesTickLabelInterpreter','latex','defaultAxesFontSize',20)
set(0,'DefaultLegendInterpreter','latex')



% Specify debt-state level
% i_b_g_mean = i_b_g_zero+50;
b_g_mean = 0.025;0.0053;
[~, i_b_g_mean] = min(abs(b_g_vec-b_g_mean));

% Arbitrary output levels
y_high = 1.07;
y_low = .93;%1.03;
[~, i_y_high] = min(abs(y_vec_3sh-y_high));
[~, i_y_low] = min(abs(y_vec_3sh-y_low));

% Axis margins
y_axlim = [.97*y_low,1.03*y_high];[.98,1.1];
b_axlim = [0.0235,.0265];
q_axlim = [min(q_g_pf(:,i_b_g_mean)),max(q_g_pf(:,i_b_g_mean))];


figg = figure;

figg.Position = [10, 10, 900,900];
% Colors
color_order =  get(gca,'colororder');
color_1 = color_order(1,:);
color_2 = color_order(2,:);
color_3 = color_order(3,:);



subplot(2,1,1)


hold on;
box on;
% Borrowing P.F.
plot(y_vec_3sh,b_g_pf(:,i_b_g_mean),'linewidth',3)

% Vertical/Horiz.line line y_high
plot([y_vec_3sh(i_y_high)-eps y_vec_3sh(i_y_high)+eps],...
    [0,b_g_pf(i_y_high,i_b_g_mean)],':k')
plot([min(y_vec_3sh) max(y_vec_3sh)],...
    [b_g_pf(i_y_high,i_b_g_mean),b_g_pf(i_y_high,i_b_g_mean)],':k')

% Vertical/Horiz. line y_low
plot([y_vec_3sh(i_y_low)-eps y_vec_3sh(i_y_low)+eps],...
    [0,b_g_pf(i_y_low,i_b_g_mean)],':k')
plot([min(y_vec_3sh) max(y_vec_3sh)],...
    [b_g_pf(i_y_low,i_b_g_mean),b_g_pf(i_y_low,i_b_g_mean)],':k')



xlim(y_axlim)
ylim(b_axlim)

xlabel('$y$')
ylabel('$\left|b^\prime\right|$')


xticks([y_vec_3sh(i_y_low) y_vec_3sh(i_y_high)])
xticklabels({'$y_{L}$','$y_{H}$'})
yticks([b_g_pf(i_y_low,i_b_g_mean) b_g_pf(i_y_high,i_b_g_mean)]) 
yticklabels({'$b^\prime_{L}$','$b^\prime_{H}$'})

legend('$b^\prime(y,b)$','location','SouthEast')
legend boxoff;

title({'Panel A: Debt Issuance Policy Function'})

subplot(2,1,2)
hold on;
box on;


% Q- Policy function
plot(y_vec_3sh,q_g_pf(:,i_b_g_mean),'linewidth',3,'color',color_1)

% Price functions varying y, keeping b_prime constant
plot(y_vec_3sh,q_g(:,i_b_g_pf(i_y_high,i_b_g_mean)),'linewidth',2,'color',color_2,'linestyle','--')
plot(y_vec_3sh,q_g(:,i_b_g_pf(i_y_low,i_b_g_mean)),'linewidth',2,'color',color_3,'linestyle','-.')



% Vertical/Horiz.line line y_high
plot([y_vec_3sh(i_y_high)-eps y_vec_3sh(i_y_high)+eps],...
    [0,q_g_pf(i_y_high,i_b_g_mean)],':k')
plot([min(y_vec_3sh) max(y_vec_3sh)],...
    [q_g_pf(i_y_high,i_b_g_mean),q_g_pf(i_y_high,i_b_g_mean)],':k')

% Vertical/Horiz. line y_low
plot([y_vec_3sh(i_y_low)-eps y_vec_3sh(i_y_low)+eps],...
    [0,q_g_pf(i_y_low,i_b_g_mean)],':k')
plot([min(y_vec_3sh) max(y_vec_3sh)],...
    [q_g_pf(i_y_low,i_b_g_mean),q_g_pf(i_y_low,i_b_g_mean)],':k')

% Horiz. Line q_HL
plot([min(y_vec_3sh) max(y_vec_3sh)],...
    [q_g(i_y_low,i_b_g_pf(i_y_high,i_b_g_mean)),q_g(i_y_low,i_b_g_pf(i_y_high,i_b_g_mean))],':k')

xlim(y_axlim)
ylim(q_axlim)   

xlabel('$y$')
ylabel('$q$')

xticks([y_vec_3sh(i_y_low) y_vec_3sh(i_y_high)])
xticklabels({'$y_{L}$','$y_{H}$'})
yticks([q_g(i_y_low,i_b_g_pf(i_y_high,i_b_g_mean))...
q_g_pf(i_y_low,i_b_g_mean) q_g_pf(i_y_high,i_b_g_mean)]) 
yticklabels({'$q_{HL}$','$q_{L}$','$q_{H}$'})

legend('$q^c(y,b)$','$q(y,b^\prime_H)$','$q(y,b^\prime_L)$','location','East') %,'Color','w','EdgeColor','w'
legend boxoff;


title({'Panel B: Conditional and Unconditional Debt Pricing Functions'})
print(hFig,'-dpdf','P:\Research\20160628_default_volatility\Codes\Replication Codes\Output\nonlinearity_2p_coarse-eps-converted-to.pdf.pdf','-opengl')
