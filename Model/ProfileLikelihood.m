% Initial guess for parameters
mod4 = [433, 4e+04, 2.37, 6.05, 2.74, 1.82e+06, 0.00601, 0.383, 4.37e+03, 5.19, 0, 12.8, 0.0039, 0.678, 0, NaN, NaN, 1, 0.165, NaN, NaN, 4e+04, 0, ];

mod6_2 = [433, 4e+04, 2.02, 7.26, 2.62, 2.37e+06, 4.47e-05, 4.32, 476, 5.19, 0, 12.8, 0.0039, 0.678, 0, NaN, NaN, 1, 0.165, NaN, NaN, 2.89e+04, 0, ];

% params = mod4;
params = mod6_2;
% pCa = 4.4;
pCa = 6.2;

modSel = [3 4 5 6 7 8 9];
% modSel = [3 4 5];

selMods = params(modSel);

% Options for optimization
opts = optimset('Display','iter', 'TolFun', 5e-3, 'Algorithm','sqp', 'TolX', 0.01, 'PlotFcns', @optimplotfval, 'MaxIter', 100);
% opts = optimset('Display', 'off');

% Pre-define storage for profile likelihoods
numSteps = 8; % Number of steps for each parameter profile

if ~exist('paramProfiles', 'var')
    paramProfiles = cell(length(selMods), 1);
end

% Loop over each parameter to create profiles
for i = 1:length(selMods)
    % Store the fixed parameter profile
    profileValues = linspace(selMods(i) * 0.5, selMods(i) * 1.5, numSteps);

    if ~isfield(paramProfiles{i}, 'likelihood')
        profileLikelihood = inf(size(profileValues));
        lastProfileLikelihood = inf(size(profileValues));
    else
        profileLikelihood = paramProfiles{i}.likelihood;
        lastProfileLikelihood = paramProfiles{i}.lastProfileLikelihood;
    end
    

    for j = 1:numSteps
        if profileLikelihood(j) < inf
            fprintf('Skipping mod %s, step %g, already done...\n', modNames{modSel(i)}, j);
            continue;
        end
        % Fix the i-th parameter
        paramsFixed = params;
        paramsFixed(modSel(i)) = profileValues(j);
        
        % Create a vector of parameters to optimize        
        modSelDif = setdiff(modSel, modSel(i));

        reducedObjective = @(p) isolateRunCombinedModel(p, modSelDif, paramsFixed, pCa); 
        
        % Initial guess for reduced dimensions
        if j == 1 ||  ~isfield(paramProfiles{i}, 'params') || isempty(paramProfiles{i}.params) ...
            || all(size(paramProfiles{i}.params) ~= [numSteps, length(modSelDif)])
            reducedParams0 = params(modSelDif);
        else
            % start with last iter to speed up
            reducedParams0 = paramProfiles{i}.params(j-1, :);
        end
        

        % Optimization
        for ijk = 1:3
            fprintf('Iterating mod %s, step %g, search %g \n', modNames{modSel(i)}, j, ijk);
            if lastProfileLikelihood(j) < 0.3
                fprintf('... skipping, already good enough.');
                continue
            elseif lastProfileLikelihood(j) < inf && lastProfileLikelihood(j) < profileLikelihood(j)*1.05;
                fprintf('... skipping 3rd optim, not enough change.');
                continue;
            end
            % save one last
           lastProfileLikelihood(j) = profileLikelihood(j);
           [reducedParams0, profileLikelihood(j)] = fminsearch(reducedObjective, reducedParams0, opts);
        end
        % Store the profiles
        reducedParams(j, :) = reducedParams0;
        paramProfiles{i} = struct('values', profileValues, 'likelihood', profileLikelihood, 'lastProfileLikelihood', lastProfileLikelihood, 'params', reducedParams);        
    end
    
    % Store the profiles
    % paramProfiles{i} = struct('values', profileValues, 'likelihood', profileLikelihood, 'likelihood', lastProfileLikelihood, 'params', reducedParams0);
end

% Plot profile likelihoods
figure;
for i = 1:length(selMods)
    nexttile;
    plot(paramProfiles{i}.values, paramProfiles{i}.likelihood, '-o');
    xlabel(['Parameter ', num2str(i)]);
    ylabel('Cost');
    title(['Profile Likelihood for Parameter ', modNames(modSel(i))]);
end
save(sprintf('profilelikelihood_%g.mat', pCa), "paramProfiles");
%% load profiles

profile4 = load(sprintf('profilelikelihood_%g.mat', 4.4)).paramProfiles;
profile6_2 = load(sprintf('profilelikelihood_%g.mat', 6.2)).paramProfiles;

% Plot profile likelihoods
figure;
for i = 1:length(selMods)
    nexttile;
    plot(profile4{i}.values, profile4{i}.likelihood, '-x', profile6_2{i}.values, profile6_2{i}.likelihood,'-o');
    xlabel(['Parameter ', num2str(i)]);
    ylabel('Cost');
    title(['Profile Likelihood for Parameter ', modNames(modSel(i))]);
end


%% Plot all pCas
pCas = [4.4, 5.5, 5.8, 6, 6.2, 11];
figure(80085);clf;
colors = lines(6);
tile_semilogx = nexttile();
hold on;
for ipCa = 1:length(pCas)
    if pCas(ipCa) < 11
        filename = sprintf('..\\Data\\AvgMava_pCa%0.1f_%gs.csv', pCas(ipCa), 0.1);     
    else
        filename = sprintf('..\\Data\\AvgRelaxedMavaSet_%gs.csv', 0.1);     
    end
    datatables{ipCa} = readtable(filename);

    hd{ipCa} = errorbar(datatables{ipCa}.Time-2,datatables{ipCa}.F,datatables{ipCa}.SD, '-', LineWidth=2, Color=colors(ipCa, :), CapSize=0);
    set([hd{ipCa}.Bar, hd{ipCa}.Line], 'ColorType', 'truecoloralpha', 'ColorData', [hd{ipCa}.Line.ColorData(1:3); 255*0.4])
    % set(h.Cap, 'EdgeColorType', 'truecoloralpha', 'EdgeColorData', [h.Cap.EdgeColorData(1:3); 255*alpha])
    ylabel('$\Theta$ (kPa)', Interpreter='latex')
end
tile_semilogx.XScale='log';




%% FUNCTIONS

function cost = isolateRunCombinedModel(optMods, modSel, params, pCa)
% just to isolate the script, so the variables can't intervene
% figure(pCa*10);    
drawPlots = false;
    rampSet = [4];
    mod = params;
    mod(modSel) = optMods;
    RunCombinedModel;
end
