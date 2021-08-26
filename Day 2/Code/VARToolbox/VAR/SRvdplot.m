function SRvdplot(VD,VARopt,INF,SUP)
% =======================================================================
% Plot the VDs computed with SR (sign restriction procedure)
% =======================================================================
% SRfevdplot(VD,VARopt,INF,SUP)
% -----------------------------------------------------------------------
% INPUT
%   - VD(:,:,:): matrix with 't' steps, the VD due to 'j' shock for 
%       'k' variable
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
filename = [VARopt.figname 'VD_SR'];
quality = VARopt.quality;
suptitle = VARopt.suptitle;
pick = VARopt.pick;

% Initialize VD matrix
nshocks = length(snames); [nsteps, nvars, ~] = size(VD);

% If one variable is chosen, set the right value for nvars
if pick<0 || pick>nvars
    error('The selected variable is non valid')
else
    if pick==0
        pick=1;
    else
        nvars = pick;
    end
end


% Define the rows and columns for the subplots
row = round(sqrt(nshocks));
col = ceil(sqrt(nshocks));

% Define a timeline
steps = 1:1:nsteps;
x_axis = zeros(1,nsteps);



%% Plot
%=========
FigSize
for ii=pick:nvars
    for jj=1:nshocks
        subplot(row,col,jj);
        plot(steps,VD(:,jj,ii),'LineStyle','-','Color',[0.01 0.09 0.44],'LineWidth',2);
        hold on
        if exist('INF','var') && exist('SUP','var')
            PlotSwathe(VD(:,jj,ii),[INF(:,jj,ii) SUP(:,jj,ii)],cmap(1));
        end
        xlim([1 nsteps]); ylim([0 100]);
        plot(x_axis,'k','LineWidth',0.5); set(gca,'Layer','top');
        title([vnames{ii} ' to ' snames{jj}], 'FontWeight','bold','FontSize',10); 
    end
    % Save
    FigName = [filename num2str(ii)];
    if quality 
        if suptitle==1
            Alphabet = char('a'+(1:nvars)-1);
            SupTitle([Alphabet(ii) ') VD of '  vnames{ii}])
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
