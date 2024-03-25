
% Data daved by save('pca4dataAdj60s.mat', 'Tarr', 'Farr') from
% AverageRamps[Ca].m
% load pca4data.mat
% load pca4dataAdj.mat
% load pca4dataAdj60s.mat
load pca11data.mat


fitfun = @(init) evalPowerFit(init, Farr, Tarr, true);
figure(2);clf;
% best fit for pca 11
x = [4.4271    0.2121    4.8964]

% best fit for pca 4.4 - 12.7 t^-1.16, using only first 10s
% x = [12.7405    1.1640   15.9898];

% best fit for pca4.4: limit the Fss to that of pca 11, using only first 10s.
% Surprisingly, the time constant is nearly same!
% note: run the AverageRampsCa to get the Ca Farr
% x = [18.8677    0.2162   4.8964];

% experiment: free Fss, using only 0.1 and 1s ramps for 10s to avoid
% messing with the remaining force buildup. The optimizer pushes the Fss
% high
% x = [16.9548    0.2781    7.0702];

% experiment: fixed Fss, using only 0.1 and 1s ramps for 10s. No big
% difference from the all ramps fit.
% x = [19.0119    0.2407    4.8964];

% Best fit for pca4.4 with Frem correction, no assumptions
% x = [10.4953    0.5869   9];
% x = [10.4953    0.21   11];

% pCa 4.4 with Frem correction, for only 0.1 and 1s ramps for 30s
% x = [11.2095    0.6120   12.2428]

% pCa 4.4 with Frem correction, assuming Fss
% x= [17.0511    0.2176    4.8964];

% % pCa 4.4, optimizing just the tail
% x = [4.3209    0.2100   13.6246];

% limit pCa time constant to the relaxed time constant
% x = [8.6226    0.2100   11.3698];
% no Ca
aspect = 1.5;
f = figure(2);
% normal size of 2-col figure on page is 7.2 inches
% matlab's pixel is 1/96 of an inch
f.Position = [300 200 7.2*96 7.2*96/aspect];
x = [4.4271    0.2121    4.8964];
[c rampShift] = fitfun(x)
f = gcf();
exportgraphics(f,'../Figures/FigDecayOverlay.png','Resolution',150)
% saveas(f, 'Figures/FigDecayOverlaypCa4.png')
% exportgraphics(f,'Figures/FigDecayOverlaypCa4.4_7.2.png','Resolution',300)
% saveas(f, 'Figures/FigDecayOverlaypCa4Corr2.png')
%% pCa data
aspect = 1.5;
% rampShift = [5.3980    0.8234    0.2223   0.0100];
pcadata = load('../pCa4dataNoAdj60sFremCorr.mat');
% pcadata = load('../pca4data60sNoFremCorr.mat');
% pcadata = load('..\pca4.4modeldata.mat');
% pcadata.FarrCorr = pcadata.Farr;pcadata.TarrCorr = pcadata.Tarr;
% Farr = pcadata.FarrCorr;Tarr = pcadata.Tarr;
Farr = pcadata.Farr;Tarr = pcadata.Tarr;
% x = [4.1240    0.2121   12.0286];
% load('pCa4dataNoAdj60sFremCorrShifted.mat')
% x = [4.0648    0.2121   12.0505];
f = figure(4);
f.Position = [300 200 7.2*96 7.2*96/aspect];
% keep the b, fit a and Tss
% x = [4.1237    0.2121   12.0289];
% x = [4.1233    0.2121   12.0290];
% x = [16.2753    0.2990    6.5963];
x = [18.5295    0.2121    3.4365];
% x = [20.8050    0.1941    4.5377]; % incomplete fit somehow
% % x = [11.3305    0.3950   13.4757]; % alternative
% alternative: fit all, incl. b: uncorrected data
% x = [11.4470    0.3895   13.3954]; % best fit in 0-10 zone
x = [11.3325    0.3949   13.4736]; % even better
% alternative: fit all, incl. b: corrected data
% x = [13.2146    0.3448    8.8122]
[c rspca] = evalPowerFit(x, Farr, Tarr, true, [], true)
f = gcf();
% exportgraphics(f,'../Figures/FigDecayOverlaypCa4.png','Resolution',150)
% exportgraphics(f,'../Figures/FigDecayOverlayModelpCa4.png','Resolution',150)
%% model - relaxed
f = figure(5);
aspect = 1.5;
f.Position = [300 200 7.2*96 7.2*96/aspect];
load ..\pca11modeldata.mat
% load pca11modeldataDoubleStates.mat
% data fit
x = [4.4271    0.2121    4.8964];
% model fit
% x = [2.2393    0.5408    6.8351];
% retuned model fit with 1e3 alphaU mod
% x = [2.0713    0.5512    6.8491];
% tail-optimized fit using Dan's parameter space
% x = [5.2749    0.3628    5.4011]
x = [3.3117    0.2878    5.8467]
x = [3.3117    0.2878 4.8964]
x = [4.7976    0.2392    4.8212];
%
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
f = figure(4);clf;
f.Position = [300 200 7.2*96 7.2*96/aspect];
pcadata = load('..\pca4.4modeldata.mat');
Farr = pcadata.Farr;Tarr = pcadata.Tarr;
x = [17.9381    0.2339    4.3359];
[c rspca] = evalPowerFit(x, Farr, Tarr, 'loglogOnly', [], true)

%%
options = optimset('Display','iter', 'TolFun', 1e-4, 'Algorithm','sqp', 'UseParallel', true, ...
    'TolX', 0.0001, 'PlotFcns', @optimplotfval, 'MaxIter', 150);
% init(1:3) = [0 0 0]
%% fit the no Ca
init = x;
fitfunOpt = @(init) evalPowerFit(init, Farr, Tarr, false);
x = fminsearch(fitfunOpt, init, options)

[c rampShift] = evalPowerFit(x, Farr, Tarr, true, [], false);
%% fit pCa data

% fixed power exponent b
% init = [x(1) x(3)]
% pcaFitFunFixB = @(x)evalPowerFit([x(1), init(2), x(2)], Farr, Tarr, false, [], true);
% x = fminsearch(pcaFitFunFixB, [init(1) init(3)], options);
% init([1 3]) = x;

% free power exponent
init = x;
pcaFitFun = @(x)evalPowerFit(x, Farr, Tarr, false, [], true);
x = fminsearch(pcaFitFun, init, options);
init = x;

x = init
c = evalPowerFit(x, Farr, Tarr, true, [], true)


% pcaFitFun = @(x)evalPowerFit(x, Farr, Tarr, false, rampShift, true);
% x = fminsearch(pcaFitFun, init, options)
% x = init;


