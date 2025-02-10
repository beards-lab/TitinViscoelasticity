rerunFitting = false;

%% Low calcium data fit
f = figure(3);
aspect = 1.5;
f.Position = [300 200 7.2*96 7.2*96/aspect];

% Averaged and 
% for iarr = 1:length(arr)
%     Farr{iarr} = dsc{ids, arr(iarr)}.datatableZDCorr.F;
%     Tarr{iarr} = dsc{ids, arr(iarr)}.datatableZDCorr.t - 10;
% 
%     i_cutoff = find(Farr{iarr} > 4 & Tarr{iarr} > 50, 1, 'last');
%     Farr{iarr} = Farr{iarr}(1:i_cutoff);
%     Tarr{iarr} = Tarr{iarr}(1:i_cutoff);
% end
% Farr{1} = [];
% Tarr{1} = [];

% Averaged data
load ../data/pca11data.mat
% Load from averaged simulation input data
rds = [100, 10, 1, 0.1];

% linear residuals, offset, c = 125, logC = 4.9
x = [3.7242    0.2039    4.8357];

% linear residuals, no offset, c = 264, logC - 0.89
% x = [8.7614    0.0790 0];

% log residuals, optimizes to no offset, logC = 0.83, cLin - 289
% x = [8.7132    0.0755    0]

if rerunFitting
    init = x;
    fitfunOpt = @(x) evalPowerFit(x, Farr, Tarr, false);
    x = fminsearch(fitfunOpt, init, options)
end

[c rampShift] = evalPowerFit(x, Farr, Tarr, true, [], false);

%% model data - relaxed
f = figure(5);
aspect = 1.5;
f.Position = [300 200 7.2*96 7.2*96/aspect];
load ..\data\pca11modeldata.mat
% load pca11modeldataDoubleStates.mat

% model fit
x = [4.7976    0.2392    4.8212];

rws_s = 4;rws_l = rws_s*4;rws_loglog = 6;
tiledlayout(3, rws_loglog+rws_l, "Padding","compact", "TileSpacing","compact");
tile_semilogx = nexttile(1, [2, rws_l]);
tile_loglog = nexttile(rws_l + 1, [3, rws_loglog]);
nexttile((rws_loglog+rws_l)*2 + 1, [1 rws_s]);
nexttile((rws_loglog+rws_l)*2 + 1 + rws_s, [1 rws_s]);
nexttile((rws_loglog+rws_l)*2 + 1 + 2*rws_s, [1 rws_s]);
nexttile((rws_loglog+rws_l)*2 + 1 + 3*rws_s, [1 rws_s]);
axes(tile_loglog);
[c rspca] = evalPowerFit(x, Farr, Tarr, 'loglogOnly', [], false)
% exportgraphics(f,'../Figures/FigDecayOverlayModelRelaxed.png','Resolution',150)
% Farr = Force;Tarr = Time;
%% pCa model
Farr = {};
Tarr = {};
for i_rd = 1:length(rds)
    % tb = readtable(['..\Data\AvgRelaxedMavaSet_' num2str(rds(i_rd)) 's.csv']);
    tb = readtable(['..\Data\AvgMavaSetpCa4.4_' num2str(rds(i_rd)) 's.csv']);
    Farr{i_rd} = tb.F;
    Tarr{i_rd} = tb.Time - tb.Time(1); % start at zero
end

% pcadata = load('..\pca4.4modeldata.mat');
% Farr = pcadata.Farr;Tarr = pcadata.Tarr;

%% pCa model
f = figure(4);clf;
f.Position = [300 200 7.2*96 7.2*96/aspect];
x = [17.9381    0.2339    4.3359];
[c rspca] = evalPowerFit(x, Farr, Tarr, 'loglogOnly', [], true)

%%

% Ca does not affect parallel stiffness - we use that from teh relaxed
x0 = 3.7242;

% fit the whole Ca
init = [x0    0.2768    5.5530];

% best all tail: c = 28.1
% init = [6.3546    0.2157    x0  -12.2928];

% best tail when relaxing the x0: c = 23. Optically much worse fit though
% init = [35.2299    1.2706    7.3737  -12.2928]

% best tail for fast(0.1 and 1)
% Tarr{1} = {};Farr{1} = {};
% Tarr{2} = {};Farr{2} = {};
% init = [10.2511    0.3704    4.8400   -9.7544]

% best tail for slow (10 and 100)
% Tarr{3} = {};Farr{3} = {};
% Tarr{4} = {};Farr{4} = {};
% init = [4.6950    0.1268    4.8400  -15.7143];

% best head does not work. Pretty bad
% init = [6.4546    0.357    x0  8.2928];


rerunFitting = true;

if rerunFitting
    options = optimset('Display','iter', 'TolFun', 1e-4, 'Algorithm','sqp', 'UseParallel', true, ...
        'TolX', 0.0001, 'PlotFcns', @optimplotfval, 'MaxIter', 150);

    pcaFitFunFixB = @(x)evalPowerFit([x(1), x(2), x0, x(3)], Farr, Tarr, false, [], true);
    x = fminsearch(pcaFitFunFixB, [init(1) init(2) init(4)], options);
    init([1, 2, 4]) = x;
end

[c rampShift] = evalPowerFit(init, Farr, Tarr, true, [], true);
c
%% init([1 3]) = x;
figure(101); 
% init = [x(1), init(2), x(2)];
x = [5.2540    0.1328    3.1201];
x = [5.502540    0.18328    3.1201]
[c rs] = evalPowerFit(x, Farr, Tarr, true, [], false)

%%
% free power exponent
% init = x;
pcaFitFun = @(x)evalPowerFit(x, Farr, Tarr, true, [], true);
x = fminsearch(pcaFitFun, init, options);
init = x;

x = init
c = evalPowerFit(x, Farr, Tarr, true, [], true)


% pcaFitFun = @(x)evalPowerFit(x, Farr, Tarr, false, rampShift, true);
% x = fminsearch(pcaFitFun, init, options)
% x = init;

%% For all iters
figure(1);clf;
% Define the parameters for each setting, including ids and arr
settings = {
    1, [2 3 4 5];  % relaxed
    1, [9 8 7 6];  % relaxed reversed
    3, [2 3 4 5];  % PNB and MAVA, no Ca    
    3, [7 8 9 11]; % max Ca, PNB and MAVA
    4, [2 3 4 5];  % 10C - no Ca, just PNB Mava
    4, [7 8 9 11]; % 10C - max Ca, PNB and MAVA
    5, [2 3 4 5];  % extracted relaxed
    6, [2 3 4 5]   % extracted activated
};

% Generate noise
noise_percentage = 0.05;

% Parameters
siginfs = 0:.25:6;  % A vector with 7 elements
num_iterations = size(settings, 1);  % Number of different sets of ids and arr

% Initialize cell arrays to store results for each iteration
all_c = cell(num_iterations, 1);
all_res = cell(num_iterations, 1);

% Loop over the number of iterations
for iter = 1:num_iterations
    ids = settings{iter, 1};
    arr = settings{iter, 2};
    
    % Assume Farr and Tarr are generated based on ids and arr
    % Sample data for demonstration; replace with actual logic for generating Farr and Tarr
clear Farr Tarr;
    for iarr = 1:length(arr)
        Farr{iarr} = dsc{ids, arr(iarr)}.datatableZDCorr.F;
        Tarr{iarr} = dsc{ids, arr(iarr)}.datatableZDCorr.t - 10;
    
        i_cutoff = find(Farr{iarr} > 4 & Tarr{iarr} > 50, 1, 'last');
        Farr{iarr} = Farr{iarr}(1:i_cutoff);
        Tarr{iarr} = Tarr{iarr}(1:i_cutoff);
    end

    % Initialize variables for storing results
    c = zeros(1, length(siginfs));
    res = zeros(length(siginfs), 2);
    init = [5, 0.1];  % Initial guess, adjust as needed

    % Optimization loop
    for isig = 1:length(siginfs)
        options = optimset('Display', 'iter', 'TolFun', 1e-4, 'TolX', 0.0001, 'MaxIter', 100);
        
        % add noise - about 10%
        for iarr = 1:length(arr)
            noise = noise_percentage * randn(size(Farr{iarr})) .* Farr{iarr};
            % Add the noise to the signal
            Farr_noi{iarr} = Farr{iarr} + noise;
        end
        
        fitfunOpt = @(sigInit) evalPowerFit([sigInit, siginfs(isig)], Farr_noi, Tarr, false);
        x = fminsearch(fitfunOpt, init, options);
        % c(isig) = fitfunOpt(x);  % This stores the objective value
        figure(2)
        c(isig) = evalPowerFit([x, siginfs(isig)], Farr_noi, Tarr, true);
        figure(1)
        res(isig, :) = x;  % This stores the optimized parameters
    end
    
    % Store the results for the current iteration
    all_c{iter} = c;
    all_res{iter} = res;
end

% save('sigmaLandscapeLog', 'all_c', 'all_res', 'siginfs')
save('sigmaLandscapeNoiseLog', 'all_c', 'all_res', 'siginfs')
%% Plot results for each iteration
    % load sigmaLandscapeLog.mat;
    % load sigmaLandscapeNoise.mat
    load sigmaLandscapeNoiseLog.mat
    
    figure(2);clf;hold on;
    num_iterations = [1 2 3 5]; % NO Ca
    % title(['No Ca']);

    % figure(2);clf;hold on;
    % num_iterations = [4 6]; % Ca
    % title(['Max Ca']);

    % figure(3);clf;hold on;
    % num_iterations = [7 8]; % Extracted
    % title(['Extracted']);

    % num_iterations = [1 2]; % NO Ca
    
    sig_cut_off = 20;

clear leg
for i = 1:length(num_iterations)
    iter = num_iterations(i);
    leg(i) = plot(siginfs, all_res{iter}(:, 2), LineWidth=2);
    
    % relative
    % scatter(siginfs, all_res{iter}(:, 2), 40, all_c{iter}/min(all_c{iter}), 'filled');
    
    % absolute
    scatter(siginfs(1:end), all_res{iter}(1:end, 2), 40, min(sig_cut_off, all_c{iter}(1:end)), 'filled');

    xlabel('siginfs');
    ylabel('res(:, 2)');
    
    colorbar; % Add colorbar to indicate the value of c
    grid on;
    % settings(iter)
end

% Convert numeric array to cell array of strings
legend_labels = arrayfun(@(x) num2str(x), num_iterations, 'UniformOutput', false);

legend(leg, legend_labels, Location="northwest")