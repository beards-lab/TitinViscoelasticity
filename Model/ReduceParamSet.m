% We try to fit pCa 6.2 with reduced set of parameters, keeping the pCa 4
% instead
modNames = {'k_p(NoCa)', 'k_d', 'n_p', 'n_U', 'n_d', 'alphaU', 'k_{PEVK,A}', 'k_{PEVK,D}', 'k_p(highCa)', 'Fss', 'b', 'c', 'd', 'mu', 'alphaF_0','k_{PEVK,A} (low Ca)', 'k_{PEVK,D} (low Ca)', 'Lref', 'delU', ...
        'AlphaU_pCa', 'Fss_pCa', 'kd_pCa', 'kDf'};
mod4 = [433, 4e+04, 2.37, 6.05, 2.74, 1.82e+06, 0.00601, 0.383, 4.37e+03, 5.19, 0, 12.8, 0.0039, 0.678, 0, NaN, NaN, 1, 0.165, NaN, NaN, 4e+04, 0, ];
mod6_2 = [433, 4e+04, 2.02, 7.26, 2.62, 2.37e+06, 4.47e-05, 4.32, 476, 5.19, 0, 12.8, 0.0039, 0.678, 0, NaN, NaN, 1, 0.165, NaN, NaN, 2.89e+04, 0, ];

% Specify your initial set of selected parameters
modSel = [3 4 5 6 7 8 9];

% Define a threshold for the cost function
costThreshold = 0.4;

% Parameters to be potentially excluded one-by-one initially
removalCandidates = modSel;

%init
bestModSel = modSel;
bestCost = inf;
bestParams = mod6_2;
foundImprovement = true;
level = 0;
while foundImprovement
    foundImprovement = false;
    for i = 1:length(removalCandidates)
        % Try removing the current candidate
        testModSel = removalCandidates;
        testModSel(i) = [];
        % testOptParams = bestOptParams(testModSel);
        
        % set the other to different set value
        testParams = bestParams;
        testParams(removalCandidates(i)) = mod4(removalCandidates(i));

        % init the search set
        testOptParams = testParams(testModSel);
                
        % Function handle for the cost function
        costFunctionHandle = @(optMods) isolateRunCombinedModel(optMods, testModSel, testParams, 6.2);

        % costFunctionHandle(testOptParams)

        % Optimize for the reduced parameters
        opts = optimset('Display','none', 'TolFun', 5e-3, 'Algorithm','sqp', 'TolX', 0.01, 'PlotFcns', @optimplotfval, 'MaxIter', 100);

        testCost = inf;
        lastTestCost = inf;
        for ijk = 1:3
            fprintf('Iterating mod %s, level %g, search %g \n', modNames{modSel(i)}, level, ijk);
            if testCost < costThreshold
                fprintf('... skipping, already good enough.');
                continue
            elseif lastTestCost < inf && lastTestCost < testCost*1.05
                fprintf('... skipping 3rd optim, not enough change.');
                continue;
            end
            % save one last
           lastTestCost = testCost;

            [testOptParams, testCost] = fminsearch(costFunctionHandle, testOptParams, opts);
        end
        
        % Evaluate cost with reduced parameter set
        if testCost < costThreshold
            bestModSel = testModSel;
            bestParams = testParams;
            bestParams(testModSel) = testOptParams;
            bestCost = testCost;
            foundImprovement = true;
            removalCandidates(i) = []; % Update candidate list
            level = level + 1;
            break; % Restart the loop with the updated candidate list
        end
    end
end

% Output the simplified parameter set
fprintf('Final reduced parameter set: %s\n', mat2str(bestModSel));
fprintf('Final optimized parameters: %s\n', mat2str(bestParams));
fprintf('Final achieved cost: %.3f\n', bestCost);
save ReduceParams2
%%
figure(80085)
tsts = mod6_2;
% tsts(finalModSel) = bestOptParams(finalModSel);
isolateRunCombinedModel([], [], tsts, 6.2)

%%


function cost = isolateRunCombinedModel(optMods, modSel, params, pCa)
% just to isolate the script, so the variables can't intervene
% figure(pCa*10);    
    drawPlots = false;
    rampSet = [4];
    mod = params;
    mod(modSel) = optMods;
    RunCombinedModel;
end
