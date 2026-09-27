%% FitFirstStretch.m
% ENTRY POINT. Stretch-hold (+ 0 ms restretch) fitting of the 2025-11-21
% refolding protocol, before any refolding parameters are fitted.
%
% Data (10 kHz per-protocol exports, NOT the 1 kHz 03_/05_Log_* files, which
% have 2-3 samples on the ~2.8 ms ramp and under-read the peak by ~25%):
%   SS{1}  first stretch, ensemble mean of the 8 *_refolding_<cond>.txt repeats
%   SS{2}  5ms_Ramp_<cond>.txt (5.7 ms ramp)     SS{3}  0.1s_Ramp_<cond>.txt
%   (active also SS{4} 1s_Ramp_Active.txt)
%   SR     release + restretch at 40 s of 0ms_refolding_<cond>.txt (0.2 ms gap)
% Model: simStretchHold.m (same dXdT, grid and states as RunCombinedModel,
%   driven by the measured length trace), optional structural variants
%   (strain-dependent drag mu(s), parallel viscosity eta, Maxwell element,
%   slip-bond detachment, smooth refolding - see its header), force seen
%   through a 2nd-order transducer (sensorFilter.m, f0 ~1.4-1.9 kHz), compared
%   on log-spaced bins (loadStretchHold.m / loadRestretch.m, modelBinned.m).
% Cost: sum over datasets of mean(((model - data)/max(data)).^2) per bin
%   (costStretchHold.m; resVariant.m for lsqnonlin; runFitJob.m for batch).
%
% Workspace config: cond ('Relax' | 'Active'), resultFile (.mat with p, sensor).

if ~exist('cond', 'var'), cond = 'Relax'; end
dd = '..\Data\2025 11 21 Export/';
rf = {'0ms', '5ms', '10ms', '50ms', '100ms', '1s', '10s', '30s'};

SS = {};
SS{1} = loadStretchHold(strcat(dd, rf, ['_refolding_' cond '.txt']), [], 200);
SS{2} = loadStretchHold([dd '5ms_Ramp_' cond '.txt'], [], 200);
SS{3} = loadStretchHold([dd '0.1s_Ramp_' cond '.txt'], [], 200);
if strcmp(cond, 'Active')
    SS{4} = loadStretchHold([dd '1s_Ramp_Active.txt'], [], 200);
end
SR = loadRestretch([dd '0ms_refolding_' cond '.txt']);
SJ = [SS, {SR}];

%% Parameter set to show
% result files hold p (params), sensor ([f0 zeta]) and note
if ~exist('resultFile', 'var')
    resultFile = sprintf('fitStretchHold_%s_best.mat', cond);
end
R = load(resultFile);
pEval = R.p; senEval = R.sensor;
fprintf('%s: %s\n', resultFile, R.note);
[costAll, partsAll] = costStretchHold(pEval, SJ, [], senEval);
fprintf('%s: total cost %.5f, per dataset %s\n', resultFile, costAll, mat2str(partsAll(:, 1)', 3));

%% Figures: first stretch evaluation and restretch
SSe = SS; %#ok<NASGU> (used by ShowFirstStretchEval)
cd ../Figures; ShowFirstStretchEval; cd ../Model;
[FbR, tfR, FfR] = modelBinned(pEval, SR, senEval);
f = figure(773); clf; f.Position = [100 100 1100 420]; tiledlayout(1, 2, 'TileSpacing', 'compact');
nexttile; hold on; box on; w = SR.t < 8e-3;
plot(1e3*SR.t(w), SR.F(w), '.', 'Color', [.6 .6 .6]);
plot(1e3*SR.tb(SR.tb < 8e-3), SR.Fb(SR.tb < 8e-3), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 3);
plot(1e3*(tfR - SR.tRs), FfR, 'b-', 'LineWidth', 1.5);
xlim([-2 8]); xlabel('t from restretch onset (ms)'); ylabel('\Theta (kPa)');
legend('data 10 kHz', 'data bins', 'model + sensor'); title(sprintf('%s restretch (0 ms gap)', cond));
nexttile; hold on; box on; pos = SR.tb > 0;
plot(SR.tb(pos), SR.Fb(pos), 'ko', 'MarkerSize', 3); plot(SR.tb(pos), FbR(pos), 'b-', 'LineWidth', 1.5);
set(gca, 'XScale', 'log'); xlabel('t from restretch onset (s)'); title(sprintf('restretch cost %.5f', partsAll(end, 1)));
exportgraphics(f, sprintf('../Figures/Restretch_%s.png', cond), 'Resolution', 110);
