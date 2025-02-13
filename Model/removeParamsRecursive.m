%Parameters to be potentially excluded one-by-one initially
modNames = {'k_p(NoCa)', 'k_d', 'n_p', 'n_U', 'n_d', '\alpha_U', 'k_{A}', 'k_{D}', 'k_p', '\Theta_{ss}', 'b', 'c', 'd', '\mu', 'alphaF_0','k_{PEVK,A} (low Ca)', 'k_{PEVK,D} (low Ca)', 'Lref', 'delU', ...
        '\alpha_U', 'Fss_pCa', 'kd_pCa', 'kDf'};
mod4 = [433, 4e+04, 2.37, 6.05, 2.74, 1.82e+06, 0.00601, 0.383, 4.37e+03, 5.19, 0, 12.8, 0.0039, 0.678, 0, NaN, NaN, 1, 0.165, NaN, NaN, 4e+04, 0, ];
mod6_2 = [433, 4e+04, 2.02, 7.26, 2.62, 2.37e+06, 4.47e-05, 4.32, 476, 5.19, 0, 12.8, 0.0039, 0.678, 0, NaN, NaN, 1, 0.165, NaN, NaN, 2.89e+04, 0, ];
mod11 = [432.5020, 40000, 2.3741    8.8100    2.5458, 8360000, 0.0171    0.6954  838.5174    5.1935    0.0000,  12.8000    0.0039    0.6780         0  NaN       NaN    1.0000    0.1646       NaN      NaN  40000           0];
% adjusted to pCa format
mod11_ = mod11;mod11_(9) = mod11(1);
resultOptim = load('pCasResultOptim').resultOptim;

%% optim the fit
% isolateRunCombinedModel([], [], mod4, 4.4, true)
% isolateRunCombinedModel([], [], mod11_, 11, true)

% for i = 1:length(resultOptim)
%     resultOptim{i}.modSel = validCombos{i}.modSel;
% end
% save resultOptimValidCombos resultOptim
numOpt = 3;
options = optimset('Display','none', 'TolFun', 1e-2, 'Algorithm','sqp', 'TolX', 0.01, 'PlotFcns', @optimplotfval, 'MaxIter', 100);

for i = [2 3 5]
    modSel = resultOptim{i}.modSel;
    resultOptim{i}.mod11 = mod6_2;
    evalLin = @(curmod) isolateRunCombinedModel(curmod, modSel, resultOptim{i}.mod11,  [11]);

    for ijk = 1:numOpt
        [x c] = fminsearch(evalLin, resultOptim{i}.mod11(modSel), options);
        resultOptim{i}.mod11(modSel) = x;
        resultOptim{i}.mod11_c = c;
    end
end
%%
for i = 1:length(resultOptim)
    modSel = resultOptim{i}.modSel;
    resultOptim{i}.mod6_2 = mod6_2;
    evalLin = @(curmod) isolateRunCombinedModel(curmod, modSel, resultOptim{i}.mod6_2,  [6.2]);

    for ijk = 1:numOpt
        [x c] = fminsearch(evalLin, resultOptim{i}.mod6_2(modSel), options);
        resultOptim{i}.mod6_2(modSel) = x;
        resultOptim{i}.mod6_2_c = c;
    end
end

%%
% Specify your initial set of selected parameters
modSel = [3 4 5 6 7 8 9];
% mod4 = rand(1, 23);
% modSel = [3 4 5 6];
if ~exist('combinations', 'var')
    combinations = {};
end

tic
[bestModSel, bestOptParams, bestCost, combinations] = findOptimalParameters(modSel, mod4, mod6_2, .4, combinations);
toc
%%
combinations = load('ohwow').combinations;
% Create arrays to hold the lengths of modSel and the costs
lengths = cellfun(@(c) numel(c.modSel), combinations);
costs = cellfun(@(c) c.cost, combinations);

% Combine lengths and costs into a matrix for sorting
sortMatrix = [lengths(:), costs(:)];

% Sort first by the length (ascending) and then by cost (ascending)
[~, sortedIndices] = sortrows(sortMatrix, [1, 2]);

% Re-arrange the combinations array based on sorted indices
combinations = combinations(sortedIndices);

% Output sorted combinations (optional)
for i = 1:length(combinations)
    fprintf('Combination %d: ModSel = [%s], NumElements = %d, Cost = %.3f\n', ...
        i, sprintf('%d ', combinations{i}.modSel), lengths(sortedIndices(i)), costs(sortedIndices(i)));
end

selectedCombo = combinations{16};
figure(2);clf;
isolateRunCombinedModel(selectedCombo.params, selectedCombo.modSel, mod4, 11, true)

%% Test a combination of parameters optim from pCa 4.4 to 11
selectedCombo = combinations{16};

% get rid of the attached ones
bool_exclModSel = ~ismember(selectedCombo.modSel, [7 8], 'legacy');
params = selectedCombo.params(bool_exclModSel);
exclModSel = selectedCombo.modSel(bool_exclModSel);

% test
mod11_ = mod4;
isolateRunCombinedModel(params, exclModSel, mod4, 11, true)

% optim
evalLin = @(curmod) isolateRunCombinedModel(curmod, exclModSel, mod4, [11]);
[x c] = fminsearch(evalLin, params, options);
mod11_(exclModSel) = x;
% run with all ramps
isolateRunCombinedModel([], [], mod11_, 11, true, [2 3 4])

% Output results
% fprintf('Optimal parameter set: %s\n', mat2str(bestModSel));
% fprintf('Optimized parameters: %s\n', mat2str(bestOptParams));
% fprintf('Achieved cost: %.3f\n', bestCost);
%% Test and list valid combos
maxValidCost = 0.4;
params = mod4;
% Create arrays to hold the lengths of modSel and the costs
lengths = cellfun(@(c) numel(c.modSel), combinations);
costs = cellfun(@(c) c.cost, combinations);

validCombos = combinations(costs < maxValidCost & lengths <= 4)

% Output sorted combinations (optional)
for i = 1:length(validCombos)
    fprintf('Combination %d: ModSel = [%s], NumElements = %d, Cost = %.3f\n', ...
        i, sprintf('%d ', validCombos{i}.modSel), length(validCombos{i}.modSel), validCombos{i}.cost);
end

params(validCombos{end}.modSel) = validCombos{end}.params;
isolateRunCombinedModel([], [], params, 6.2, true)
%% Optim all params for the valid combos
if ~exist('resultOptim', 'var')
    resultOptim = {};
end
numOpt = 4;
options = optimset('Display','none', 'TolFun', 1e-2, 'Algorithm','sqp', 'TolX', 0.01, 'PlotFcns', @optimplotfval, 'MaxIter', 50);


for i = 1:length(validCombos)

    modSel = validCombos{i}.modSel;

    if length(resultOptim) < i || ~isfield(resultOptim{i}, 'mod6')
        resultOptim{i}.mod6 = mod6_2;
    end
    evalLin = @(curmod) isolateRunCombinedModel(curmod, modSel, mod4, [6]);

    for ijk = 1:numOpt
        [x c] = fminsearch(evalLin, resultOptim{i}.mod6(modSel), options);
        resultOptim{i}.mod6(modSel) = x;
        resultOptim{i}.mod6_c = c;
    end

    if ~isfield(resultOptim{i}, 'mod5_8')
        resultOptim{i}.mod5_8 = resultOptim{i}.mod6;
    end
    evalLin = @(curmod) isolateRunCombinedModel(curmod, modSel, mod4, [5.8]);
    for ijk = 1:numOpt
        [x c] = fminsearch(evalLin, resultOptim{i}.mod5_8(modSel), options);
        resultOptim{i}.mod5_8(modSel) = x;
        resultOptim{i}.mod5_8_c = c;
    end    
    
    if ~isfield(resultOptim{i}, 'mod5_5')
        resultOptim{i}.mod5_5 = resultOptim{i}.mod5_8;
    end
    evalLin = @(curmod) isolateRunCombinedModel(curmod, modSel, mod4,  [5.5]);
    for ijk = 1:numOpt
        [x c] = fminsearch(evalLin, resultOptim{i}.mod5_5(modSel), options);
        resultOptim{i}.mod5_5(modSel) = x;
        resultOptim{i}.mod5_5_c = c;
    end    
end

% save pCasResultOptim resultOptim

%% PLOT valid combos and parameter values

for ires = 1:length(resultOptim)
sqsum = sum([resultOptim{ires}.mod5_5_c, resultOptim{ires}.mod5_8_c, resultOptim{ires}.mod6_c].^2);
    fprintf('%g: sqsum %g. Costs 5.5: %g, Costs 5.8: %g, Costs 6: %g \n', ires, sqsum, ...
        resultOptim{ires}.mod5_5_c, resultOptim{ires}.mod5_8_c, resultOptim{ires}.mod6_c);
        

    if sqsum > 1
        continue;
    end
% figure(ires + 1000);clf;
modSet = [mod4;resultOptim{ires}.mod5_5;resultOptim{ires}.mod5_8;resultOptim{ires}.mod6;mod6_2;resultOptim{ires}.mod11];
pcax = [4.4, 5.5, 5.75, 6, 6.2, 11];
% tiledlayout(1,size(modSel, 2));
modSel = resultOptim{ires}.modSel;
for i = 1:size(modSel, 2)
    nexttile();
    plot(-pcax, modSet(:, modSel(i)), 'x-', LineWidth=3);
    title(modNames{modSel(i)});
end



end

%% Try to optimize param per param

modSet = [mod4;mod4;mod4;mod4;mod4;mod11];
pcax = [4.4, 5.5, 5.75, 6, 6.2, 11];

% first round of optim
% modSet = [mod4;mod4;mod11];
% pcax = [4.4, 5.5, 11];
pcaSel = [1 2]

% Second round of optim - we already have those
disp(['params = [' sprintf('%1.3g, ', modSet(1, :)) '];'])
disp(['params = [' sprintf('%1.3g, ', modSet(2, :)) '];'])
disp(['params = [' sprintf('%1.3g, ', modSet(3, :)) '];'])
mod4 = [433, 4e+04, 2.37, 5.78, 2.74, 8.71e+05, 0.00564, 0.383, 4.31e+03, 5.19, 0, 12.8, 0.0039, 0.678, 0, NaN, NaN, 1, 0.165, NaN, NaN, 4e+04, 0, ];
mod5_5 = [433, 4e+04, 2.37, 6.27, 2.74, 3e+06, 0.00567, 0.383, 3.78e+03, 5.19, 0, 12.8, 0.0039, 0.678, 0, NaN, NaN, 1, 0.165, NaN, NaN, 4e+04, 0, ];
mod11 = [433, 4e+04, 2.37, 8.38, 2.55, 4.3e+06, 0.0171, 0.695, 442, 5.19, 0, 12.8, 0.0039, 0.678, 0, NaN, NaN, 1, 0.165, NaN, NaN, 4e+04, 0, ];
% init for second round - now we use matrix of all to prevent
% non-monotonic, but subselct only pca to optim
modSet = [mod4;mod5_5;mod5_5;mod5_5;mod5_5;mod11];
pcax = [4.4, 5.5, 5.75, 6, 6.2, 11];
%% get the other bound right
pcaSel = [5 6];
% tweak pca11
modSet(6, [7 8]) = [NaN, NaN];
% fix pca4 attachment
modSet(1, 7) = NaN;
modSet(5, 7) = 5e-3;
%% round 3 - get the middle right
pcaSel = [2 3 4];

%%
pcaSel = [1];
selectedCombo = combinations{16};
modSel = selectedCombo.modSel;
% modSel = [1]
optMods = modSet(pcaSel, modSel);
optModsInLine = reshape(optMods, [], 1);
% isolateRunCombinedModelAllpCas(optMods, modSel, modSet, pcax, pcaSel)
%
opts = optimset('Display','none', 'TolFun', 1e-2, 'Algorithm','sqp', 'TolX', 0.01, 'PlotFcns', @optimplotfval, 'MaxIter', 500);

for ijk = 1:4
    evalLin = @(curmod) isolateRunCombinedModelAllpCas(curmod, modSel, modSet, pcax, pcaSel)
    optModsInLine = fminsearch(evalLin, optModsInLine, opts);
end

optMods = reshape(optModsInLine, length(pcaSel), length(modSel));
% modSet = [mod4;mod4;mod4;mod4;mod4;mod4];
% pcax = [4.4, 5.5, 5.75, 6, 6.2, 11];
modSet(pcaSel, modSel) = optMods;
% modSet(end, :) = mod11_

% save optimModSetAll modSet 
%% test fit
% pcaSel = [5 6];
isolateRunCombinedModelAllpCas(optModsInLine, modSel, modSet, pcax, pcaSel, true)
%% test just monotonicity

isolateRunCombinedModelAllpCas([], modSel, modSet, pcax, [])
%% Run all ramps and combinations - Test modSet
pcax = [4.4, 5.5, 5.75, 6, 6.2, 11];   
cf = figure(80085); clf;
aspect = 1.5;
set(cf, 'Position', [500  300  7.2*96 7.2*96/aspect]);

modSel = [4 6 7 9];
ylims = [60, 54, 44, 32, 22, 16];
tiledlayout('flow', 'TileSpacing','Compact')
for i_pca = 1:size(modSet, 1)
    
    % i_pca = 5;
    % modSet(5, 7) = 1e-3;
    % modSet(5, 4) = modSet(6, 4);
    % modSet(1, 7) = modSet(2, 7) ;
    % modSet(6, 5) = modSet(5, 5) ;
    params = modSet(i_pca, :);
    nexttile;hold on;
    isolateRunCombinedModel([], [], params, pcax(i_pca), true, [1 2 3 4])
    leg = gca().Legend;
    leg.NumColumns = 1;
    title(sprintf('pCa %g', pcax(i_pca)), Interpreter="latex");
    xlim([1e-2 2e2]);
    ylim([0, ylims(i_pca)]);
    xticks([1e-2 1e0 1e2]);
    switch (i_pca)
        case {1, 2, 3}
        % top three
        % set(gca,'Xticklabel',[]);
        otherwise
            xlabel('Time (s)', Interpreter="latex");
    end
    switch (i_pca)
        case {1, 4}
        % left two keep ylabel
        otherwise
        ylabel("")
    end

    set(gca, 'FontSize', 12);
end

%% plot interesting parameter values into one plot
figure(1000);clf;
% selectedParams = [3 4 5 6 7 8 9];
selectedParams = modSel;
% selectedParams = 1:23;
% tiledlayout(4,1);
tiledlayout('flow');
for i_selpar = 1:length(selectedParams)
    nexttile(i_selpar);hold on;
    % for ires = [2 3 5]
        % modSet = [mod4;resultOptim{ires}.mod5_5;resultOptim{ires}.mod5_8;resultOptim{ires}.mod6;mod6_2;resultOptim{ires}.mod11];
        % pcax = [4.4, 5.5, 5.75, 6, 6.2, 11];        
        plot(-pcax, modSet(:, selectedParams(i_selpar)), 'x:', LineWidth=3);
    % end
end
    % legend(modNames())
%% Show the resulting matrix
for i = 1:size(modSet,1)
    fprintf('%10.4g ', modSet(i,:));  % Prints each row with 4 significant figures
    fprintf('\n');  % New line after each row
end

%%
ms = [...
       433      4e+04       2.37      5.797       2.74  8.658e+05   0.005381      0.383       4345       5.19          0       12.8          0      0.678          0        NaN        NaN          1      0.165        NaN        NaN      4e+04          0 
       433      4e+04       2.37      6.211       2.74  2.837e+06   0.005213      0.383       3998       5.19          0       12.8          0      0.678          0        NaN        NaN          1      0.165        NaN        NaN      4e+04          0 
       433      4e+04       2.37      6.526       2.74  3.184e+06   0.002699      0.383       3109       5.19          0       12.8          0      0.678          0        NaN        NaN          1      0.165        NaN        NaN      4e+04          0 
       433      4e+04       2.37       7.96       2.74  1.807e+07  3.261e-05      0.383       1623       5.19          0       12.8          0      0.678          0        NaN        NaN          1      0.165        NaN        NaN      4e+04          0 
       433      4e+04       2.37      8.505       2.74  1.876e+07  3.478e-06      0.383      887.7       5.19          0       12.8          0      0.678          0        NaN        NaN          1      0.165        NaN        NaN      4e+04          0 
       433      4e+04       2.37      9.035       2.74  2.668e+07        NaN        NaN      512.3       5.19          0       12.8          0      0.678          0        NaN        NaN          1      0.165        NaN        NaN      4e+04          0 
    ];
%% PEVK KO

ms_PevkKO = ms(1, :);
ms_PevkKO(7) = 1e-12;
isolateRunCombinedModel([], [], ms_PevkKO, [4.4], true, [1 2 3 4])

%% experiment
ms_mu0 = ms(end, :);
% ms_mu0(1, 14) = 5e-1;
isolateRunCombinedModel([], [], ms_mu0, [11], true, [4])
%%
ms_mu0 = ms(1, :);
% ms_mu0(1, 14) = 5e-1;
isolateRunCombinedModel([], [], ms_mu0, [4.4], true, [3])

%% Fig 10 - sin refolding experiment
drawPlots = true;
rampSet = [3];

mod = ms(1, :);    
    
plotInSeparateFigure = true;

simtype = 'sin';
pCa = 11;
mod(15) = 2;
RunCombinedModel;


%% fit params a func
cf = figure(8008135);clf;
fitHill(modSet, modSel(3), pcax)
tiledlayout('flow', TileSpacing='compact');
% tiledlayout(2, 2)
% tiledlayout(1, 4)
modNames_modSel = modNames(modSel);
modNames_units = {"-", "s^{-1}", "s^{-1}", "s^{-1}"}
for i_ms = 1:length(modSel)
    nexttile();hold on;
    fitHill(modSet, modSel(i_ms), pcax);
    ylabel(sprintf('$%s (%s)$', modNames_modSel{i_ms}, modNames_units{i_ms}), Interpreter="latex");
    % title(sprintf('Fitting param %s', modNames{modSel(i_ms)}));
    set(gca,'TickLabelInterpreter','latex')
end
fontsize(12, 'points');
aspect = 1.5;

set(cf, 'Position', [500  300  7.2*96 7.2*96/aspect]);
% set(cf, 'Position', [500  300  3.5*96 3.5*96/aspect]);
K = [5.9076, 5.9544, 5.7493, 5.8618];
K_mean = mean(K)
K_std = std(K)
%% Functions below
function [bestModSel, bestParams, bestCost, combinations] = findOptimalParameters(modSel, params, optimizedParams, costThreshold, combinations)
    modNames = {'k_p(NoCa)', 'k_d', 'n_p', 'n_U', 'n_d', 'alphaU', 'k_{PEVK,A}', 'k_{PEVK,D}', 'k_p(highCa)', 'Fss', 'b', 'c', 'd', 'mu', 'alphaF_0','k_{PEVK,A} (low Ca)', 'k_{PEVK,D} (low Ca)', 'Lref', 'delU', ...
        'AlphaU_pCa', 'Fss_pCa', 'kd_pCa', 'kDf'};

    % Internal recursive function
    function [localBestModSel, localBestParams, localBestCost, level] = exploreSubset(currentModSel, currentParams, level)
        
        % Create a logical array to check if each combination's modSel matches the target modSel
        isFound = any(cellfun(@(c) isfield(c, 'modSel') && isequal(sort(c.modSel), sort(currentModSel)), combinations));

        
        % If no parameters left, return high cost        
        if isempty(currentModSel) || isFound
            localBestModSel = [];
            localBestParams = [];
            localBestCost = inf;
            return;
        end

        % Evaluate cost for current subset
        % currentParams = params(currentModSel);

        % Optimization options
        opts = optimset('Display','none', 'TolFun', 1e-2, 'Algorithm','sqp', 'TolX', 0.01, 'PlotFcns', @optimplotfval, 'MaxIter', 150);

        currentCost = inf;
        lastTestCost = inf;
        for ijk = 1:4
        fprintf('Level %g, Optimizing for %s, search %g.. ', level, sprintf('%d ', currentModSel), ijk)
            if currentCost < costThreshold
                fprintf('... skipping, already good enough. \n');
                break;
            elseif lastTestCost < inf && lastTestCost < currentCost*1.05
                fprintf('... skipping 3rd optim, not enough change.\n');
                break;
            end
            % save one last
           lastTestCost = currentCost;
           [currentParams, currentCost] = fminsearch(@(p) isolateRunCombinedModel(p, currentModSel, params, 6.2), currentParams, opts);
        end 
        level = level + 1;
        
        
        combinations{end+1} = struct('modSel', currentModSel, 'params', currentParams, 'cost', currentCost);
        if currentCost < costThreshold
            % Store the valid combination and its cost
            
            fprintf('ok, cost = %.3f\n', currentCost);
        else
            localBestModSel = [];
            localBestParams = [];
            localBestCost = inf;
            fprintf('fail.\n');
            return;
        end

        % Initialize the best found in this scope
        localBestModSel = currentModSel;
        localBestParams = currentParams;
        localBestCost = currentCost;

        % Explore further by removing each parameter recursively
        for i = 1:length(currentModSel)
            
            reducedModSel = currentModSel;
            fprintf('Reducing by %d\n', reducedModSel(i));
            reducedModSel(i) = [];

            initParams = currentParams;
            initParams(i) = [];

            [candidateModSel, candidateParams, candidateCost, level] = exploreSubset(reducedModSel, initParams, level);
            
            % Update best set if found a suitable candidate
            if candidateCost < costThreshold && length(candidateModSel) < length(localBestModSel)
                localBestModSel = candidateModSel;
                localBestParams = candidateParams;
                localBestCost = candidateCost;
            end
        end
    end

    % Start recursion
    [bestModSel, bestParams, bestCost, level] = exploreSubset(modSel, optimizedParams(modSel), 0);
    fprintf('Level %g\n', level);

    % Output all valid combinations at the end of the function
    fprintf('Valid parameter combinations under the threshold:\n');
    for k = 1:length(combinations)
        fprintf('Combination %d: ModSel = [%s], Cost = %.3f\n', ...
            k, sprintf('%d ', combinations{k}.modSel), combinations{k}.cost);
    end
    % save ohwow;
end

function cost = isolateRunCombinedModelAllpCas(optModsInLine, modSel, modSet, pCas, pcaSel, drawPlots)

pcax = pCas(pcaSel);
optMods = reshape(optModsInLine, length(pcax), length(modSel));
% modSet = [mod4;mod4;mod4;mod4;mod4;mod4];
% pcax = [4.4, 5.5, 5.75, 6, 6.2, 11];

modSet(pcaSel, modSel) = optMods;

if nargin < 6
    drawPlots = false;
end
    rampSet = [4];
total_cost = 0;
    for i_pca = 1:length(pcax)
        mod = modSet(i_pca, :);
        mod(modSel) = optMods(i_pca, :);
        pCa = pcax(i_pca);
        RunCombinedModel;
        total_cost(i_pca) = cost;
    end
% enforce monotonicity over whole modSet
lambda = 1000;

% c_mono = lambda * min(...
%     sum(max(0, optMods(1:end-1) - optMods(2:end))), ... % Non-decreasing penalty
%     sum(max(0, optMods(2:end) - optMods(1:end-1)))); % non increasing penalty
% 
% Compute differences between consecutive elements for each column
differences = diff(modSet./modSet(1, :));
differences = differences(:, modSel);

% Compute non-decreasing penalties: positive values imply non-monotonicity
non_decreasing_penalties = sum(max(0, -differences));

% Compute non-increasing penalties: positive values imply non-monotonicity
non_increasing_penalties = sum(max(0, differences));

% Compute the minimum penalties for each column
column_penalties = min(non_decreasing_penalties, non_increasing_penalties);

% Total penalty (sum of penalties for all columns, scaled by lambda)
c_mono = lambda * sum(column_penalties);

    cost = sum(total_cost) + c_mono; 
end


function cost = isolateRunCombinedModel(optMods, modSel, params, pCa, drawPlots, rampSet)
% just to isolate the script, so the variables can't intervene
% figure(pCa*10);    

if nargin < 5
    drawPlots = false;
    rampSet = [4];
elseif nargin < 6 
    rampSet = [4];
end

    mod = params;
    mod(modSel) = optMods;
    % cost = sum((optMods/10 - 2).^2);
    plotInSeparateFigure = false;
    RunCombinedModel;
end


function fitHill(modSet, modSel, pcax)
cutXaxis = true;
% Example data points (replace with actual values)
x_data = pcax;
y_data = modSet(:, modSel)';  % Dependent variable

filter = ~isnan(y_data);
% y_data = y_data(filter);
% x_data = x_data(filter);
y_data(~filter) = 0;

%% Choose initial parameter guesses [A, K, n]
A0 = max(y_data);     % Approximate max value
A_0 = min(y_data);
K0 = 6; %median(x_data);  % Midpoint estimate
n0 = 100;               % Initial Hill coefficient

params0 = [A0, K0, n0, A_0];

% Define Hill function models
hill_rising = @(p, x) (p(1) * x.^p(3)) ./ (p(2).^p(3) + x.^p(3)) + A_0;  % Rising Hill function
hill_decreasing = @(p, x) p(1) ./ (1 + (x./p(2)).^p(3) + A_0);           % Decreasing Hill function

% Perform curve fitting using nonlinear least squares
options = optimset('Display', 'none', 'TolFun', 1e-6);
params_rising = lsqcurvefit(hill_rising, params0, x_data, y_data, [], [], options);
params_decreasing = lsqcurvefit(hill_decreasing, params0, x_data, y_data, [], [], options);

% Plot results
if cutXaxis
    x_data(end) = 8;
end

% Generate fine x-values for plotting
x_fine = linspace(min(x_data)*0.95, max(x_data)*1.05, 100);

% Evaluate fitted models
y_fit_rising = hill_rising(params_rising, x_fine);
y_fit_decreasing = hill_decreasing(params_decreasing, x_fine);

% Compute residuals
residual_rising = norm(y_data - hill_rising(params_rising, x_data));
residual_decreasing = norm(y_data - hill_decreasing(params_decreasing, x_data));

% Determine best fit model
if residual_rising < residual_decreasing
    best_fit = 'Rising Hill Function';
    best_params = params_rising;
    y_best_fit = y_fit_rising;
else
    best_fit = 'Decreasing Hill Function';
    best_params = params_decreasing;
    y_best_fit = y_fit_decreasing;
end

scatter(x_data, y_data, 80, 'kx','linewidth', 3, 'DisplayName', 'Data Points');
% plot(x_fine, y_fit_rising, 'b-', 'LineWidth', 2, 'DisplayName', 'Hill (Rising)');
% plot(x_fine, y_fit_decreasing, 'g--', 'LineWidth', 2, 'DisplayName', 'Hill (Decreasing)');
plot(x_fine, y_best_fit, 'k-', 'LineWidth', 1.5, 'DisplayName', ['Best Fit: ', best_fit]);

if cutXaxis
    % Set custom x-ticks with a break
    xticks([4 5 6 7 8]);  % Include 7 and 11 for a break
    xticklabels({'-4', '-5', '-6', '...', '11'}); % Add '...' as break
end
%%
legend(gca(), 'off');
set(gca,'TickLabelInterpreter','latex')
xlabel('pCa', Interpreter='latex');
% ylabel('Param value');

% title(['Hill Function Fit - Best Model: ', best_fit]);

% Display best-fit parameters
disp(['Best Fit Model: ', best_fit]);
disp(['A = ', num2str(best_params(1))]);
disp(['K = ', num2str(best_params(2))]);
disp(['n = ', num2str(best_params(3))]);

end