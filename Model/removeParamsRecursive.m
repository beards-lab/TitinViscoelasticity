%Parameters to be potentially excluded one-by-one initially
modNames = {'k_p(NoCa)', 'k_d', 'n_p', 'n_U', 'n_d', 'alphaU', 'k_{PEVK,A}', 'k_{PEVK,D}', 'k_p(highCa)', 'Fss', 'b', 'c', 'd', 'mu', 'alphaF_0','k_{PEVK,A} (low Ca)', 'k_{PEVK,D} (low Ca)', 'Lref', 'delU', ...
        'AlphaU_pCa', 'Fss_pCa', 'kd_pCa', 'kDf'};
mod4 = [433, 4e+04, 2.37, 6.05, 2.74, 1.82e+06, 0.00601, 0.383, 4.37e+03, 5.19, 0, 12.8, 0.0039, 0.678, 0, NaN, NaN, 1, 0.165, NaN, NaN, 4e+04, 0, ];
mod6_2 = [433, 4e+04, 2.02, 7.26, 2.62, 2.37e+06, 4.47e-05, 4.32, 476, 5.19, 0, 12.8, 0.0039, 0.678, 0, NaN, NaN, 1, 0.165, NaN, NaN, 2.89e+04, 0, ];

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

% isolateRunCombinedModel(bestOptParams, bestModSel, params, 6.2)

% Output results
fprintf('Optimal parameter set: %s\n', mat2str(bestModSel));
fprintf('Optimized parameters: %s\n', mat2str(bestOptParams));
fprintf('Achieved cost: %.3f\n', bestCost);
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

%%
for ires = 1:length(resultOptim)
sqsum = sum([resultOptim{ires}.mod5_5_c, resultOptim{ires}.mod5_8_c, resultOptim{ires}.mod6_c].^2);
    fprintf('Costs 5.5: %g, Costs 5.8: %g, Costs 6: %g, sum %g, sqsum %g \n',...
        resultOptim{ires}.mod5_5_c, resultOptim{ires}.mod5_8_c, resultOptim{ires}.mod6_c, ...
        sum([resultOptim{ires}.mod5_5_c, resultOptim{ires}.mod5_8_c, resultOptim{ires}.mod6_c]), ...
        sqsum);

    if sqsum > 1
        continue;
    end
figure(ires + 1000);clf;
modSet = [mod4;resultOptim{ires}.mod5_5;resultOptim{ires}.mod5_8;resultOptim{ires}.mod6;mod6_2, mod11];
pcax = [4.4, 5.5, 5.75, 6, 6.2, 11];
% tiledlayout(1,size(modSel, 2));
modSel = validCombos{ires}.modSel;
for i = 1:size(modSel, 2)
    nexttile();
    plot(-pcax, modSet(:, modSel(i)), 'x-', LineWidth=3);
    title(modNames{modSel(i)});
end



end

%% function [bestModSel, bestParams, bestCost] = findOptimalParameters(modSel, params, costThreshold)
% 
%     % Internal recursive function
%     function [localBestModSel, localBestParams, localBestCost, level] = exploreSubset(currentModSel, level)
%         % If no parameters left, return high cost
%         if isempty(currentModSel)
%             localBestModSel = [];
%             localBestParams = [];
%             localBestCost = inf;
%             return;
%         end
% 
%         % Evaluate cost for current subset
%         currentParams = params(currentModSel);
%         % options = optimset('Display', 'off');
%         opts = optimset('Display','none', 'TolFun', 5e-3, 'Algorithm','sqp', 'TolX', 0.01, 'PlotFcns', @optimplotfval, 'MaxIter', 500);
% 
% 
%         fprintf('Level %g, Optimizng for %s.. ', level, sprintf('%d ', currentModSel))
%         [currentParams, currentCost] = fminsearch(@(p) isolateRunCombinedModel(p, currentModSel, params, 6.2), currentParams, opts);
%         level = level + 1;
% 
%         if currentCost < costThreshold
%             localBestModSel = currentModSel;
%             localBestParams = currentParams;
%             localBestCost = currentCost;
%             fprintf('ok, go on.\n ');
%         else
%             localBestModSel = [];
%             localBestParams = [];
%             localBestCost = inf;
%             % doesn't work here anymore
%             fprintf('fail. \n ');
%             return;
%         end
% 
%         % Explore further by removing each parameter recursively
%         for i = 1:length(currentModSel)            
%             reducedModSel = currentModSel;
%             fprintf('Reducing by %d \n', reducedModSel(i));
%             reducedModSel(i) = [];            
%             [candidateModSel, candidateParams, candidateCost, level] = exploreSubset(reducedModSel, level);
% 
%             % Update best set if found better candidate
%             if candidateCost < costThreshold && length(candidateModSel) < length(localBestModSel)
%                 localBestModSel = candidateModSel;
%                 localBestParams = candidateParams;
%                 localBestCost = candidateCost;
%             end
%         end
%     end
% 
%     % Start recursion
%     [bestModSel, bestParams, bestCost, level] = exploreSubset(modSel, 0);
%     fprintf('Level %g', level)
% 
% end
%%
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

%%

function cost = isolateRunCombinedModel(optMods, modSel, params, pCa, drawPlots)
% just to isolate the script, so the variables can't intervene
% figure(pCa*10);    
if nargin < 5
    drawPlots = false;
end
    rampSet = [4];
    mod = params;
    mod(modSel) = optMods;
    % cost = sum((optMods/10 - 2).^2);
    RunCombinedModel
end
