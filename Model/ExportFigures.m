%% Resimulate figures
clear;
saveFigures = true;
%% Figure 1
% only the bottom two panels
f = figure(101);clf;
pCa = 11;
drawFig1 = true;
RunCombinedModel;
f = gcf;

if saveFigures
    saveas(f, '../Figures/Figure1_BC', 'png');
end
%% Figure 2
% exportgraphics(f,'../Figures/RepreRamps.png','Resolution',150)
f = figure(2);clf;
PlotFig2AndDie = true;
addpath ../DataProcessing/
try
    AverageRamps

    if saveFigures
        saveas(f, '../Figures/Figure2', 'png');
    end

catch e
    disp('Skipping figure 2, whichi requires postprocessed data. Start running LoadBakersExpPassiveCa.m on raw data.')
end


%% Figure 3
% FigFitDecayOverlay
f = figure(3);clf;
aspect = 1.5;
f.Position = [300 200 7.2*96 7.2*96/aspect];
x = [3.7242    0.2039    4.8357];
% loads data cleaned and averaged in AverageRamps
load("../data/pca11data.mat");
[c rampShift] = evalPowerFit(x, Farr, Tarr, true, [], false);
leg = gcf().Children(1);
leg.Position = [0.4519    0.8571    0.5096    0.1359];

if saveFigures
    saveas(f, '../Figures/Figure3', 'png');
end

%% Figure 4
clearvars -except saveFigures
f = figure(4); clf;
pCa = 11;

drawPlots = true;
plotDetailedPlots = true;
plotInSeparateFigure = true;
RunCombinedModel;

if saveFigures
    saveas(f, '../Figures/Figure4', 'png');
end

%% Figure 5
% states of relaxed
clearvars -except saveFigures
statesFig = figure(5); clf;
placeHolderFig = figure(105);clf;
pCa = 11;

drawPlots = true;
plotDetailedPlots = false;
plotInSeparateFigure = true;
drawAllStates = 3;
rampSet = [3];

RunCombinedModel;

if saveFigures
    saveas(statesFig, '../Figures/Figure5', 'png');
end

%% Figure 6
% Model, just with high pca
clearvars -except saveFigures
f = figure(6); clf;
pCa = 4.51;
% rampSet = [4];
drawPlots = true;
plotDetailedPlots = true;
plotInSeparateFigure = true;
RunCombinedModel;

if saveFigures
    saveas(f, '../Figures/Figure6', 'png');
end

%% Figure 7

% states of high ca
clearvars -except saveFigures
statesFig = figure(7); clf;
placeHolderFig = figure(105);clf;
pCa = 4.51;

drawPlots = true;
plotDetailedPlots = false;
plotInSeparateFigure = true;
drawAllStates = 3;
rampSet = [3];

RunCombinedModel;

if saveFigures
    saveas(statesFig, '../Figures/Figure7', 'png');
end

%% Figure 8
clearvars -except saveFigures

pcax = [4.51, 5.5, 5.75, 6, 6.2, 11];   
ylims = [60, 54, 44, 32, 22, 16];

cf = figure(80085); clf;

aspect = 1.5;
set(cf, 'Position', [500  300  7.2*96 7.2*96/aspect]);
tiledlayout('flow', 'TileSpacing','Compact')
for i_pca = 1:size(pcax, 2)
    pCa = pcax(i_pca);
    
    nexttile;hold on;
    plotInSeparateFigure = false;
    % rampSet = [4];
    clear params;
    RunCombinedModel;
 
    % leg = gca().Legend;
    % leg.NumColumns = 1;
    title(sprintf('pCa = %g', pcax(i_pca)), Interpreter="latex");
    xlim([1e-2 2e2]);     ylim([0, ylims(i_pca)]);
    xticks([1e-2 1e0 1e2]);
    legend(gca(), 'off');
    switch (i_pca)
        case {1, 2, 3}
        % top three
        % set(gca,'Xticklabel',[]);
        otherwise
            xlabel('t (s)', Interpreter="latex");
    end
    switch (i_pca)
        case {1, 4}
        % left two keep ylabel
        otherwise
        ylabel("")
    end

    set(gca, 'FontSize', 12);
end

if saveFigures
    saveas(cf, '../Figures/Figure8', 'png');
end

%% Figure 9
clearvars -except saveFigures

% identified param vals on a Hill curve
paramNames = {'\F_{ss}', 'n_ss' , 'k_p'    , 'n_p'    , 'k_d'   , 'n_d'    , '\alpha_U'  , 'n_U'    , '\mu'    , '\delta_U'    , 'k_{A}'   , 'k_{D}'   };
paramSel = [7 8 3 11];

% SAme same as in RunCombinedModel
paramSet = [...
      5.19       12.8       4345       2.37      4e+04       2.74  8.658e+05      5.797      0.678      0.165   0.005381      0.383 
      5.19       12.8       3998       2.37      4e+04       2.74  2.837e+06      6.211      0.678      0.165   0.005213      0.383 
      5.19       12.8       3109       2.37      4e+04       2.74  3.184e+06      6.526      0.678      0.165   0.002699      0.383 
      5.19       12.8       1623       2.37      4e+04       2.74  1.807e+07       7.96      0.678      0.165  3.261e-05      0.383 
      5.19       12.8      887.7       2.37      4e+04       2.74  1.876e+07      8.505      0.678      0.165  3.478e-06      0.383 
      5.19       12.8      512.3       2.37      4e+04       2.74  2.668e+07      9.035      0.678      0.165        NaN        NaN 
    ];
pcax = [4.51, 5.5, 5.75, 6, 6.2, 11];

cf = figure(9);clf;
tiledlayout('flow', TileSpacing='compact');

modNames_modSel = paramNames(paramSel);
modNames_units = {"-", "s^{-1}", "s^{-1}", "s^{-1}"}
for i_ms = 1:length(paramSel)
    nexttile();hold on;
    fitHill(paramSet, paramSel(i_ms), pcax);
    ylabel(sprintf('$%s (%s)$', modNames_modSel{i_ms}, modNames_units{i_ms}), Interpreter="latex");
    % title(sprintf('Fitting param %s', modNames{modSel(i_ms)}));
    set(gca,'TickLabelInterpreter','latex')
end
fontsize(12, 'points');
aspect = 1.5;

set(cf, 'Position', [500  300  7.2*96 7.2*96/aspect]);
% set(cf, 'Position', [500  300  3.5*96 3.5*96/aspect]);
K = [5.9076, 5.9544, 5.7493, 5.8618];
K_mean = round(mean(K), 2)
K_std = round(std(K), 2)

if saveFigures
    saveas(cf, '../Figures/Figure9', 'png');
end
%% Figure 10 - hysteresis - choose relaxed or high calcium? Can add refolding too.
clearvars -except saveFigures
cf = figure(10);clf;
drawPlots = true;

% one second period
rampSet = [3];
    
plotInSeparateFigure = true;
% pCa 4.5
params = [5.19       12.8       4345       2.37      4e+04       2.74  8.658e+05      5.797      0.678      0.165   0.005381      0.383];
pCa = 4.51;
% pCa 11
% params = [5.19       12.8      512.3       2.37      4e+04       2.74  2.668e+07      9.035      0.678      0.165        NaN        NaN];
% pCa = 11;

% params(9) = 0.07;
simtype = 'sin';
params(9) = 5e-1;
RunCombinedModel;


if saveFigures
    saveas(cf, '../Figures/Figure10', 'png');
end

%% Figure 11 - PEVK knockout
clearvars -except saveFigures
cf = figure(11);clf;
fontsize(12, 'points');

pCa = 4.51;
% standard best fit
params = [ 5.19       12.8       4345       2.37      4e+04       2.74  8.658e+05      5.797      0.678      0.165   0.005381      0.383 ];
% knock out the PEVK attachment
params(11) = 1e-9;

% default plot
drawPlots = true;
plotDetailedPlots = false;
plotInSeparateFigure = false;
RunCombinedModel;
fontsize(12, 'points');

% Loop through each line and change it to dashed
lines = findobj(gca, 'Type', 'Line');
for i = 1:length(lines)
    set(lines(i), 'LineStyle', '-.');
end

% NOW - optimized for PEVK knockout
params = [5.19, 12.8, 7.56e+03, 2.75, 5.4e+04, 2.63, 1.87e+05, 5.85, 1.31, 0.165, 1e-09, 0.383, ];
drawPlots = false;
RunCombinedModel;

% Add to plot
for j = max(rampSet):-1:1
    plot(Time{j},Force{j},':', 'linewidth',2, 'Color', 'k'); 
end

% get rid of extra legends
leg = legend; 
leg.Interpreter = "tex"
leg.String{2} = "PEVK knockout"; leg.String{3} = "PEVK knockout - reoptimized";
leg.String(4:end) = []; leg.NumColumns = 1;

% general
fontsize(12, 'points');
aspect = 1.5;set(cf, 'Position', [500  300  7.2*96 7.2*96/aspect]);

if saveFigures
    saveas(cf, '../Figures/Figure11', 'png');
end
