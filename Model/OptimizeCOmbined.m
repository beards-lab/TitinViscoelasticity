close all
clear;
modNames = {'k_p(NoCa)', 'k_d', 'n_p', 'n_U', 'n_d', 'alphaU', 'k_{PEVK,A}', 'k_{PEVK,D}', 'k_p(highCa)', 'Fss', 'b', 'c', 'd', 'mu', 'alphaF_0','k_{PEVK,A} (low Ca)', 'k_{PEVK,D} (low Ca)', 'Lref', 'delU', ...
        'AlphaU_pCa', 'Fss_pCa', 'kd_pCa', 'kDf'};
runOptim = false;

%% Baseline
% Nice fit including pCa4.4 0.1s, predicting rest
params = [468, 3.83e+04, 2.3, 9, 2.33, 8.36e+06, 4.98, 84.9, 1.73e+03, 4.89, 1.01e-08, 12.8, 0.00389, 0.678, 0, NaN, NaN, 1, 0.175, NaN, NaN, 5.04e+04, 0, ];
modSel = 1:length(params);

tic
evalCombined(params(modSel), params, modSel, [11 4.4])
toc
% modSel = [1:6 7 8 9 14 20 22]

%% save baseline
m = params;
evalCombined(m, params, 1:length(mod), [11 4.4])
f = figure(144)
exportgraphics(f,sprintf('../Figures/Fig_pCa%g.png', 4.4),'Resolution',150)
f = figure(210)
exportgraphics(f,sprintf('../Figures/Fig_pCa%g.png', 11),'Resolution',150)

%% Knockouts - no PEVK attachment
m = params;
m(7) = 0;% PEVK attachment
m(8) = 1e9; % PEVK detachment
evalCombined(m, params, 1:length(mod), [4.4])
f = figure(144)
exportgraphics(f,sprintf('../Figures/FigKnockoutPEVK_pCa%g.png', 4.4),'Resolution',150)

%% Knockouts - less viscosity
m = params;
m(14) = 1e-1;% 
evalCombined(m, params, 1:length(mod), [11])
f = figure(144)
exportgraphics(f,sprintf('../Figures/FigKnockoutVisc_pCa%g.png', 4.4),'Resolution',150)
f = figure(210)
exportgraphics(f,sprintf('../Figures/FigKnockoutVisc_pCa%g.png', 11),'Resolution',150)

%% Choose one optimalization or run them one after another
%% test in GA
if runOptim
    
    % parpool
    ga_Opts = optimoptions('ga', ...
        'PopulationSize',64, ...            % 250
        'Display','iter', ...
        'MaxStallGenerations',8, ...  % 10
        'UseParallel',true);
    modSel = 1:14;
    % Ng = length(modSel);
    ub = 100*ones(1, length(mod));
    ub([3, 4, 5, 10, 14]) = [10 10 10, 10, 10];
    lb = 0.001*ones(1, length(mod));
    lb([3, 4, 5, 10, 14]) = [.1 .1 .1 .1 .1];
    evalLogCombined = @(logMod) evalCombined(10.^logMod, mod, modSel);
    init = log10(mod(modSel));
    
    [p_OptimGA,Res_OptimGA,~,~,FinPopGA,FinScoreGA] = ...
        ga(evalLogCombined,length(modSel), ...
        [],[],[],[],...
        log10(lb(modSel)),log10(ub(modSel)),[],ga_Opts);
    
    mod(modSel) = 10.^p_OptimGA;
    % use fminserach afterwards
end
%% Optimization - linear parameter space
if runOptim
    options = optimset('Display','iter', 'TolFun', 1e-3, 'Algorithm','sqp', 'TolX', 0.01, 'PlotFcns', @optimplotfval, 'MaxIter', 500);
    modSel = [2 3 4 5 6 7 8 9 10 11 12 13]; % selects modifiers to optimize for
    
    init = params(modSel);
    evalLin = @(optMods) evalCombined(optMods, params, modSel, [11])
    x = fminsearch(evalLin, init, options);
    params(modSel) = x;
end
%% optim in log param space
if runOptim

    modSel = [2 3 4 5 6 7 8 9 10 11 12 13]; % selects modifiers to optimize for
    
    init = max(-10, log10(params(modSel)));
    evalLogCombined = @(logMod) evalCombined(10.^logMod, params, modSel, [4.4]);
    x = fminsearch(evalLogCombined, init, options);
    params(modSel) = 10.^x; 
    % list params
    a = [modSel; params(modSel)]; sprintf('%d: %1.3g\n', a(:))
    % list params for save
    disp(['params = [' sprintf('%1.3g, ', params(:)) '];'])
end
%% Optimize for middle pCa's
if runOptim

    params(22) = NaN;
    mod5_5 = params;mod5_8 = params;mod6=params;
    modSel = [7 9];
    
    evalLin = @(curmod) evalCombined([curmod params(2)*params(9)/params(1)], mod5_5, [modSel 22], [5.5]);
    x = fminsearch(evalLin, mod5_5(modSel), options);
    mod5_5(modSel) = x;
    
    evalLin = @(curmod) evalCombined([curmod params(2)*params(9)/params(1)], mod5_8, [modSel 22], [5.75]);
    x = fminsearch(evalLin, mod5_5(modSel), options);
    mod5_8(modSel) = x;
    
    evalLin = @(curmod) evalCombined([curmod params(2)*params(9)/params(1)], mod6, [modSel 22], [6]);
    x = fminsearch(evalLin, mod5_8(modSel), options);
    mod6(modSel) = x;
end
%% Optimized for middle pCa's - saved last iteration
params(22) = NaN;
modSel = [7 9];
mod5_5 = params;mod5_8 = params;mod6 = params;
mod5_5(modSel) = [0.000165, 1.2e+03];
mod5_8(modSel) = [0.000188, 1.04e+03];
mod6(modSel) = [0.000238, 396];

%% Running all
modSel = 1:22;
evalCombined(mod5_5, params, modSel, [5.5]);
evalCombined(mod5_8, params, modSel, [5.75]);
evalCombined(mod6, params, modSel, [6]);

%% Compare params for different pCa's
mod = params;
mod11 = params;
mod11(modSel) = [1e-6 mod(1)]
modSet = [mod;mod5_5;mod5_8;mod6;mod11];
pcax = [4.4, 5.5, 5.75, 6, 11];
figure(40);clf;
tiledlayout(1,4);
for i = 1:size(modSel, 2)
    nexttile();
    semilogy(-pcax, modSet(:, modSel(i)), 'x-');
    legend(modNames{modSel(i)});
end

disp(['mod5_5 = [' sprintf('%1.3g, ', mod5_5(modSel)) '];'])
disp(['mod5_8 = [' sprintf('%1.3g, ', mod5_8(modSel)) '];'])
disp(['mod6 = [' sprintf('%1.3g, ', mod6(modSel)) '];'])

%%
function totalCost = evalCombined(optMods, mod, modSel, pCas)

    if nargin < 4
        % evaluate pCa 4.4 and 10 by default
        pCas = [4.4 10];
    end
    %normal - optimizing for all
    % modSel = 1:15;

    % modSel = [1 2 3 5 6 10];
    % optimizing only subset of mods
    % mod = [1.1697    1.0418    0.9774    0.9737    0.9858    1.0265    0.9403    1.0837    0.9889    0.8988 1 1 1];
    % mod = [1.16970000000000	0.928400000000000	0.977400000000000	1.02340000000000	1.01370000000000	1.10320000000000	0.937900000000000	1.19500000000000	0.909900000000000	0.898800000000000	1	1	1 1 1];
    % mod = [0.0165    0.7035    0.4403    1.0234    1.0077    0.5754    0.9379    1.1950    0.9099    0.8988    1.0000    1.1421    1.4792    1.1156    2.9834];
    % mod = [0.0185    0.8479    0.4307    1.0234    1.0326 0.5971    0.9379    1.1950    0.9099    0.8988 1.0000    1.4450    0.7510    1.2811    2.7365];
    % for absolute average
    % mod = [0.0343    0.7101    0.4731    1.0234    1.0916    1.9353    0.9379 1.1950    0.9099    0.8988    0.5952    2.0416    0.7510    1.2811 4.1891];
    % optim for -log10 weighting
    % mod = [0.928351480405723	0.928351480405723	1.01367539052550	1.02336567490158	1.01367539052550	1.10319213342611	0.937882365838957	1.19500150970587	0.909890571615859 1 1 1 1];
    % reset the search
    % mod = ones(1, 15);
    % pCa 11 decay with reduced mu
    % mod =  [0.0299    0.8084    0.4798    1.0000    1.0862    1.6973    1.0000    1.0000    1.0000    1.2983    1.0000    1.0000    1.0000    0.1    0.50000];
    % modSel = [1 2 3 5 6 10];
    

    % % mod([1:4 6:10 13]) = optMods;
    % modSel = [11, 12, 13];
    % modSel = [1 2 3 5 6 11 12 15];
    % modSel = [7, 8, 9];
    % modSel = [1 2 3 5 6 10];
    
    % store the init
    baseMods = mod;
    mod(modSel) = optMods;

    drawPlots = true;
    % drawPlots = false;
    totalCost = 0;

    if drawPlots && 1
        figInd = 100;
        try 
            set(groot,'CurrentFigure',figInd); % replace figure(indFig) without stealing the focus
        catch 
            f = figure(figInd); % if indFig is not an exsting figure it creates it (and steal the focus)
            f.Name = 'Parameter modifiers';
        end      
        clf;
        % plotParams(baseMods, [-1, 4]);hold on;
        plotParams(mod./baseMods, [-1, 1]);hold off;
    end


    %% pCa 4
    if ismember(4.4, pCas)
        pCa = 4.4;
        plotOnBackground(drawPlots, pCa);
        % RunCombinedModel;
        cost = isolateRunCombinedModel(mod, pCa, drawPlots);
        totalCost = totalCost + cost;
    end
    %% pCa 5.5
    if ismember(5.5, pCas)
        pCa = 5.5;
        plotOnBackground(drawPlots, pCa);
        % RunCombinedModel;
        cost = isolateRunCombinedModel(mod, pCa, drawPlots);
        totalCost = totalCost + cost;
    end    
    %% pCa 5.75
    if ismember(5.75, pCas)
        pCa = 5.75;
        plotOnBackground(drawPlots, pCa);
        % RunCombinedModel;
        cost = isolateRunCombinedModel(mod, pCa, drawPlots);
        totalCost = totalCost + cost;
    end

    %% pCa 6
    if ismember(6, pCas)
        pCa = 6;
        plotOnBackground(drawPlots, pCa);
        % RunCombinedModel;
        cost = isolateRunCombinedModel(mod, pCa, drawPlots);
        totalCost = totalCost + cost;
    end

    %% no Ca, but evaluate PEVK binding
    if ismember(10, pCas)    
        pCa = 10;
        plotOnBackground(drawPlots, pCa);
        % RunCombinedModel;
        cost = isolateRunCombinedModel(mod, pCa, drawPlots);
        totalCost = totalCost + cost*10;
    end
    %% no Ca, PEVK binding not evaluated
    if ismember(11, pCas)    
        pCa = 11;
        plotOnBackground(drawPlots, pCa);
        % RunCombinedModel;
        cost = isolateRunCombinedModel(mod, pCa, drawPlots);
        totalCost = totalCost + cost*10;
    end    
% return
end 

function plotOnBackground(drawPlots, pCa, cols)
    if nargin < 3
        cols = 2;
    end
    figInd = 100 + round(pCa*10);
    if drawPlots
        try 
            set(groot,'CurrentFigure',figInd); % replace figure(indFig) without stealing the focus
        catch 
            f = figure(figInd); % if indFig is not an exsting figure it creates it (and steal the focus)
            f.Name = ['pCa ' num2str(pCa)];
            if cols == 1
                aspect = 2;
                % normal size of 2-col figure on page is 7.2 inches
                % matlab's pixel is 1/96 of an inch
                f.Position = [300 200 3.5*96 3.5*96/aspect];
            else
                aspect = 1.5;
                f.Position = [300 200 7.2*96 7.2*96/aspect];
            end
        end
    end
end


function cost = isolateRunCombinedModel(mod, pCa, drawPlots)
% just to isolate the script, so the variables can't intervene
    % drawPlots = true;
    RunCombinedModel;
end

function plotParams(mod, mm, resetGca)

    % modNames = {'k_p(NoCa)', 'k_d', 'n_p', 'n_U', 'n_d', 'alphaU', 'k_{PEVK,A}', 'k_{PEVK,D}', 'k_p(highCa)', 'Fss', 'b', 'c', 'd', 'mu', 'alphaF_0','k_{PEVK,A} (low Ca)', 'k_{PEVK,D} (low Ca)', 'Lref', 'delU'};
    modNames = {'k_p(NoCa)', 'k_d', 'n_p', 'n_U', 'n_d', 'alphaU', 'k_{PEVK,A}', 'k_{PEVK,D}', 'k_p(highCa)', 'Fss', 'b', 'c', 'd', 'mu', 'alphaF_0','k_{PEVK,A} (low Ca)', 'k_{PEVK,D} (low Ca)', 'Lref', 'delU', ...
        'AlphaU_pCa', 'Fss_pCa', 'kd_pCa', 'kDf'};

    indxs = strcat(string(1:length(modNames)), ':');
    modNames = strcat(indxs, modNames);

    if nargin < 2 || isempty(mm)
        % try
        %     mima = get(gca, 'RLim');
        %     mi = mima(1);ma=mima(2);
        % catch
            mi = min(floor(log10(mod)));
            ma = mi + 3;
        % end
    else
        mi = mm(1);
        ma = mm(2);
    end
    n = length(mod);
    
    polarplot(linspace(2*pi/n, 2*pi, n), max(1e-3, log10(mod) - mi), 'x-', LineWidth=2);
    
    if ~strcmp('degrees', get(gca, 'ThetaAxisUnits')) || (nargin >= 3 && resetGca)
        % degrees means it is not adjusted yet 
        % or the param has not been provided 
        % or we want to reset
        return
    end
    rng = mi:ma;
    set(gca, 'ThetaAxisUnits', 'radians', ...
        'thetatick', linspace(2*pi/n, 2*pi, n), 'thetaticklabel', modNames, ...
        'Rlim', [min(rng) max(rng)]-mi, 'RTick', rng - mi, 'RTickLabel', 10.^rng);
end

function [c0, cost_sap, cost_sam] = runSa(pCa, mod, saSet)
    drawPlots = false;
    cost_sap = []; % SA plus
    cost_sam = []; % SA minus

    % pCa = 4;
    if drawPlots
        figure(100);
    end
    cost = isolateRunCombinedModel(mod, pCa, drawPlots);
    c0 = cost;
    %
    % saSet = [1:9 12:15];
    % saSet = [1:8 13];
    % saSet = [14 15];
    
    SAFact = 1.05;
    for i_m = saSet
        mod(i_m) = mod(i_m)*SAFact;
        fprintf('Mod %g is up to %g..', i_m, mod(i_m));
        % figure(i_m)
        cost = isolateRunCombinedModel(mod, pCa, drawPlots);
        cost_sap(i_m) = cost;
        
        mod(i_m) = mod(i_m)/SAFact/SAFact;
        fprintf('costing %1.4e€ and down to %g...', cost, mod(i_m));
        cost = isolateRunCombinedModel(mod, pCa, drawPlots);
        cost_sam(i_m) = cost;
        mod(i_m) = mod(i_m)*SAFact;
        fprintf('costing %1.4e€. \n', cost);
    end
end