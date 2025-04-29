%% Resimulate figures
clear;
saveFigures = true;
%% Figure 1
% only the bottom two panels
f = figure(101);clf;
clearvars -except saveFigures
pCa = 11;
drawFig1 = true;
RunCombinedModel;
f = gcf;

add_panel_labels(["", "(d)", "(c)"]);
 
if saveFigures
    saveas(gcf, '../Figures/Figure1_BC', 'png');
    saveas(gcf, '../Figures/Figure1_BC', 'fig');
end

%% Figure 2
% exportgraphics(f,'../Figures/RepreRamps.png','Resolution',150)
f = figure(2);clf;
PlotFig2AndDie = true;
addpath ../DataProcessing/
try
    AverageRamps
    %%
    add_panel_labels(["(b)", "(a)"], 15);

    if saveFigures
        saveas(f, '../Figures/Figure2', 'png');
        saveas(f, '../Figures/Figure2', 'fig');
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

add_panel_labels(["(b)", "", "(a)", "", "(c)"], 35);

if saveFigures
    f = gcf;
    saveas(f, '../Figures/Figure3', 'png');
    saveas(f, '../Figures/Figure3', 'fig');
end

%% Figure 4
clearvars -except saveFigures
f = figure(4); clf;
pCa = 11;
alphaF_0 = 0.0;
drawPlots = true;
plotDetailedPlots = true;
plotInSeparateFigure = true;
RunCombinedModel;

add_panel_labels(["(f)", "(e)", "(d)", "(c)", "(b)", "(a)"], 15);

if saveFigures
    saveas(f, '../Figures/Figure4', 'png');
    saveas(f, '../Figures/Figure4', 'fig');    
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

% Move the labels manually
add_panel_labels(["(b)", "", "", "", "(a)"], 15);

if saveFigures
    statesFig = gcf;
    saveas(statesFig, '../Figures/Figure5', 'png');
    saveas(statesFig, '../Figures/Figure5', 'fig');
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

add_panel_labels(["(f)", "(e)", "(d)", "(c)", "(b)", "(a)"], 15);

if saveFigures
    saveas(f, '../Figures/Figure6', 'png');
    saveas(f, '../Figures/Figure6', 'fig');
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

%% Move the labels manually
add_panel_labels(["(c)", "", "", "","","","","(b)","(a)",], 15);
%%
if saveFigures
    statesFig = gcf;
    saveas(statesFig, '../Figures/Figure7', 'png');
    saveas(statesFig, '../Figures/Figure7', 'fig');
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

add_panel_labels(["(f)", "(e)", "(d)", "(c)", "(b)", "(a)"], 15);


if saveFigures
    cf = gcf;
    saveas(cf, '../Figures/Figure8', 'png');
    saveas(cf, '../Figures/Figure8', 'fig');
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
   % {'\alpha_U'}    {'n_U'}    {'k_p'}    {'k_{A}'}
modNames_units = {"s$^{-1}$", "-", "kPa", "s$^{-1}$"}
for i_ms = 1:length(paramSel)
    nexttile();hold on;
    fitHill(paramSet, paramSel(i_ms), pcax);
    ylabel(sprintf('$%s$ (%s)', modNames_modSel{i_ms}, modNames_units{i_ms}), Interpreter="latex");
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

add_panel_labels(["(d)", "(c)", "(b)", "(a)"], 15, [], [-25, 20]);

if saveFigures
    saveas(cf, '../Figures/Figure9', 'png');
    saveas(cf, '../Figures/Figure9', 'fig');
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

add_panel_labels(["(b)", "(c)", "(a)"], 15);


if saveFigures
    cf = gcf;
    saveas(cf, '../Figures/Figure10', 'png');
    saveas(cf, '../Figures/Figure10', 'fig');
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
xlabel('$t$ (s)', Interpreter='latex')
if saveFigures
    saveas(cf, '../Figures/Figure11', 'png');
    saveas(cf, '../Figures/Figure11', 'fig');
end

return

%% Graphical abstract figure - comparison relaxed / high ca in a single plot

clearvars -except saveFigures
cf = figure(101); clf;hold on;
pCa = 11;
alphaF_0 = 0.0;
drawPlots = true;
plotDetailedPlots = false;
plotInSeparateFigure = false;
RunCombinedModel;

Time11 = Time; Force11 = Force;

pCa = 4.51;
clear params;
RunCombinedModel;
Time4 = Time;
Force4 = Force;

%% get any other line fort he legend - have to manually select that first though
hNew = gco; 

lgd = legend;
existingHandles = lgd.PlotChildren;
existingLabels = lgd.String;

newHandles = [(existingHandles); hNew]; newLabels = [existingLabels(:); {'Relaxed'}];
legend(newHandles, newLabels); xlabel('t (s)')

if saveFigures
    saveas(cf, '../Figures/Figure_GraphicalAbstract', 'png');
    saveas(cf, '../Figures/Figure_GraphicalAbstract', 'fig');
end

%% Figure supporting figure 01 - Proof the rates are stabilized at higher Ca
% we need more force to have the same transition rate

% Run Relaxed
clearvars -except saveFigures
pCa = 11;
RunCombinedModel;
Fp_11 = Fp; RU_11 = RU;
clearvars -except Fp_11 RU_11 saveFigures;

% Run high-Ca
pCa = 4.51;
RunCombinedModel;
Fp_4 = Fp; RU_4 = RU;

% plot it
cf = figure(102);clf; nexttile; hold on;
for i_ng = 1:Ng
    plot(Fp_11(:, i_ng), RU_11(:, i_ng), '-', 'DisplayName', sprintf('No Ca: State %d → %d', i_ng, i_ng+1))
end

for i_ng = 1:Ng
    plot(Fp_4(:, i_ng), RU_4(:, i_ng), '--', 'DisplayName', sprintf('High Ca: State %d → %d', i_ng, i_ng+1))
end

xlabel('$\Theta$ (kPa)', Interpreter='latex'); ylabel(['$U_{n\rightarrow n+1}(n, s)$  (s$^{-1}$)'], Interpreter='latex');
legend show; leg = legend;leg.Location = "eastoutside";

if saveFigures
    saveas(cf, '../Figures/Figure_S01_rates_vs_force', 'png');
    saveas(cf, '../Figures/Figure_S01_rates_vs_force', 'fig');
end


%% functions

function labelSubplots(figHandle, labels)
% Adds labels (A, B, C, ...) to subplots or UI panels in a figure.
% figHandle: handle to the figure
% labels: optional cell array of custom labels. If omitted, uses 'A', 'B', ...

if nargin < 2
    labels = arrayfun(@(k) char('A'+(k-1)), 1:numel(findall(figHandle, 'type', 'axes')), 'UniformOutput', false);
end

axList = findall(figHandle, 'type', 'axes');
axList = flipud(axList);  % axes are stored in reverse creation order

for k = 1:min(numel(axList), numel(labels))
    ax = axList(k);
    label = labels{k};

    % Slightly shrink axes to make room
    pos = get(ax, 'Position');
    newPos = pos;
    newPos(2) = pos(2) + 0.03;   % shift up a bit
    newPos(4) = pos(4) - 0.03;   % shrink height
    set(ax, 'Position', newPos);

    % Add annotation text label in figure space
    axUnits = get(ax, 'Units');
    set(ax, 'Units', 'normalized');
    pos = get(ax, 'Position');
    set(ax, 'Units', axUnits);

    x = pos(1) - 0.02;
    y = pos(2) + pos(4) - 0.02;
    annotation(figHandle, 'textbox', [x y 0.05 0.05], ...
        'String', label, 'FontWeight', 'bold', ...
        'EdgeColor', 'none', 'FontSize', 14, ...
        'VerticalAlignment', 'top');
end
end

function add_panel_labels(labelList, shrinkPixels, skipLabel, shift)
    % Adds panel labels (A, B, C, ...) to current figure's axes or tiles
    % labelList: optional cell array of custom labels {'A','B','C',...}
    % shrinkPixels: vertical space (in pixels) to subtract from each axis for label
    % skipLabel: logical array indicating which axes to skip labeling (but still shrink)

    if nargin < 1 || isempty(labelList)
        labelList = arrayfun(@(x) char(64 + x), 1:26, 'UniformOutput', false); % Default: 'A', 'B', ...
    end
    if nargin < 2
        shrinkPixels = 10; % Default shrink amount in pixels
    end

    fig = gcf;
    figPos = getpixelposition(fig); % [left bottom width height] in pixels
    figHeight = figPos(4);
    figWidth = figPos(3);

    allAxes = findall(fig, 'Type', 'axes');
    % Ignore legends, colorbars etc.
    allAxes = allAxes(~arrayfun(@(a) isa(a, 'matlab.graphics.illustration.Legend') || ...
                                  strcmp(get(a, 'Tag'), 'Colorbar'), allAxes));

    % Sort axes by position (top to bottom, left to right)
    axPos = arrayfun(@(a) get(a, 'Position'), allAxes, 'UniformOutput', false);
    axPos = cat(1, axPos{:});
    [A, idx] = sortrows(axPos, [2,1]);  % Top to bottom, then left to right
    allAxes = allAxes(idx);
    labelList = labelList(idx);    
    
    if nargin < 3 || isempty(skipLabel)
         skipLabel = cellfun(@isempty, labelList);
    end

    if nargin < 4
         shift = [0, 0];
    end

    % Convert pixel height to normalized units
    shrinkNorm = shrinkPixels / figHeight;

    for i = 1:length(allAxes)
        ax = allAxes(i);
        pos = get(ax, 'Position');

        % Shrink height by fixed normalized amount
        newPos = pos;
        newPos(4) = max(pos(4) - shrinkNorm, 0); % Don't let it go negative

        set(ax, 'Position', newPos);

        if ~skipLabel(i)
            % Add label in the freed vertical space
            labelHeightNorm = shrinkPixels / figHeight;
            labelBox = [newPos(1) - 10 / figWidth + shift(1)/figWidth, newPos(2) + newPos(4) + shift(2)/figHeight, 0.03, labelHeightNorm];

            annotation('textbox', labelBox, ...
                'String', labelList{i}, ...
                'FontWeight', 'bold', ...
                'FontSize', 14, ...
                'EdgeColor', 'none', ...
                'VerticalAlignment', 'bottom', ...
                'HorizontalAlignment', 'left');
        end
    end
end
