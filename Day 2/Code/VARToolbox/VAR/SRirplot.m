function SRirplot(IR,VARopt,INF,SUP)
% =======================================================================
% Plot the IRs computed with SR (sign restriction procedure)
% =======================================================================
% SRirplot(IR,VARopt,INF,SUP)
% -----------------------------------------------------------------------
% INPUT
%   - IR(:,:,:): matrix with periods, variable, shock
%   - VARopt: options of the VAR (from VARmodel and SR)
% -----------------------------------------------------------------------
% OPTIONAL INPUT
%   - INF: lower error band
%   - SUP: upper error band
% =======================================================================
% Ambrogio Cesa Bianchi, March 2015
% ambrogio.cesabianchi@gmail.com


%% Check inputs
%===============================================
if ~exist('VARopt','var')
    error('You need to provide VAR options (VARopt from VARmodel)');
end
% If there is VARopt check that vnames and snames are not empty
vnames = VARopt.vnames;
snames = VARopt.snames;
if isempty(vnames)
    error('You need to add label for endogenous variables in VARopt');
end
if isempty(snames)
    error('You need to add label for shocks in VARopt');
end


%% Define some parameters
%===============================================
filename = [VARopt.figname 'IR_SR_'];
quality = VARopt.quality;
suptitle = VARopt.suptitle;
pick = VARopt.pick;

% Initialize IR matrix
nshocks = length(snames); [nsteps, nvars, ~] = size(IR);

% If one shock is chosen, set the right value for nshocks
if pick<0 || pick>nvars
    error('The selected shock is non valid')
else
    if pick==0
        pick=1;
    else
        nshocks = pick;
    end
end

% Define the rows and columns for the subplots
row = round(sqrt(nvars));
col = ceil(sqrt(nvars));

% Define a timeline
steps = 1:1:nsteps;
x_axis = zeros(1,nsteps);


%% Plot
%================================================
FigSize(20,14)
for jj=pick:nshocks                
    for ii=1:nvars
        subplot(row,col,ii);
        plot(steps,IR(:,ii,jj),'LineStyle','-','Color',[0.01 0.09 0.44],'LineWidth',2);
        hold on
        if exist('INF','var') && exist('SUP','var')
            PlotSwathe(IR(:,ii,jj),[INF(:,ii,jj) SUP(:,ii,jj)],cmap(1));
        end
        xlim([1 nsteps]);
        plot(x_axis,'k','LineWidth',0.5); set(gca,'Layer','top');
        title([vnames{ii} ' to ' snames{jj}], 'FontWeight','bold','FontSize',10); 
    end
    % Save
    FigName = [filename num2str(jj)];
    if quality 
        if suptitle==1
            Alphabet = char('a'+(1:nshocks)-1);
            SupTitle([Alphabet(jj) ') IR to a shock to '  snames{jj}])
        end
        set(gcf, 'Color', 'w');
        export_fig(FigName,'-pdf','-png','-painters')
    else
        print('-dpng','-r100',FigName);
        print('-dpdf','-r100',FigName);
    end
    clf('reset');
end

close all
