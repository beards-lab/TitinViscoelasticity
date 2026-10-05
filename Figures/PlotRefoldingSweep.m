%% PlotRefoldingSweep.m
% ENTRY POINT. Double-ramp (stretch - release - slack gap - restretch)
% refolding protocol of 2025-11-21, simulated as ONE continuous sequence with
% the current best stretch-hold fit (fitStretchHold_<cond>_best.mat, seen
% through its own fitted transducer, sensorFilter.m); only the refolding
% rate alphaF_0 is swept (Model/SweepRefolding.m -> Model/refoldSweep.mat).
%  Fig 1  whole sequence per slack duration: full time course (linear time)
%         + zoom on the first stretch and on the restretch (same axes).
%  Fig 2  peaks vs slack duration: restretch peak (kPa), restretch / first
%         peak, and force 1 s into the restretch hold / same for the first.
% Peaks: data and model-through-sensor both smoothed over 0.7 ms (one
% transducer period) so ringing does not set the peak.
% Workspace config (optional): cond ('Relax' or 'Active').

addpath('../Model');
if ~exist('cond', 'var'), cond = 'Relax'; end
if ~exist('R', 'var') || ~isfield(R, 'Fr'), R = load('../Model/refoldSweep.mat'); end
c = find(strcmp(R.conds, cond));
for g = 1:numel(R.rf)                % raw traces are not stored in the sweep
    if ~isfield(R.SS{c, g}, 'F')
        Sr = loadProtocol(['../Data/2025 11 21 Export/' R.rf{g} '_refolding_' cond '.txt']);
        R.SS{c, g}.t = Sr.t; R.SS{c, g}.F = Sr.F; R.SS{c, g}.Lraw = Sr.Lraw;
    end
end
nG = numel(R.rf); nA = numel(R.aGrid);
cData = [.55 .55 .55]; cRef = [.80 .80 .80]; cA0 = [42 120 214]/255;
cBest = [235 104 52]/255; cFib = [.15 .15 .15];
seqBlue = [205 226 251; 158 197 244; 109 167 236; 57 135 229; 37 106 191; 24 79 149; 13 54 107]/255;
aCol = interp1(linspace(0, 1, size(seqBlue, 1)), seqBlue, linspace(0, 1, nA));
aLab = @(a) sprintf('%.3g', R.aGrid(a));
gapLab = {'0 ms', '5 ms', '10 ms', '50 ms', '100 ms', '1 s', '10 s', '30 s'};
sm = @(y, dt) movmean(y, max(1, round(7e-4/dt)));      % 0.7 ms smoothing

%% Per gap x alpha: restretch cost, peaks, 1 s force
cRs = nan(nG, nA); pk1M = cRs; pkRM = cRs; f1M = cRs; fRM = cRs;
pk1D = nan(nG, 1); pkRD = pk1D; f1D = pk1D; fRD = pk1D; ringD = pk1D; ringM = pk1D;
[bb, ab] = butter(3, [1100 2300]/(0.5*1e4));
for g = 1:nG
    S = R.SS{c, g}; kR = numel(S.ev) - 1; dtd = median(diff(S.t));
    rs = S.evBin == kR; sc = max(S.Fb(S.evBin == 1));
    wR = find(cellfun(@(t) t(1) <= S.ev(kR) + 1e-4 && t(end) >= S.ev(kR), R.tf{c, g, 1}), 1, 'last');
    pkD = @(t0) max(sm(S.F(S.t >= t0 & S.t <= t0 + 0.02), dtd));
    pk1D(g) = pkD(S.ev(1)); pkRD(g) = pkD(S.ev(kR));
    at1 = @(k) find(S.evBin == k & S.tb - S.ev(k) >= 1, 1);
    f1D(g) = S.Fb(at1(1)); fRD(g) = S.Fb(at1(kR));
    for a = 1:nA
        M = R.Fb{c, g, a}; tf = R.tf{c, g, a}; Ff = R.Ff{c, g, a};
        cRs(g, a) = mean(((M(rs) - S.Fb(rs))/sc).^2);
        pkM = @(w, t0) max(sm(Ff{w}(tf{w} >= t0 & tf{w} <= t0 + 0.02), 1e-5));
        pk1M(g, a) = pkM(1, S.ev(1)); pkRM(g, a) = pkM(wR, S.ev(kR));
        f1M(g, a) = M(at1(1)); fRM(g, a) = M(at1(kR));
    end
end
[~, aBest] = min(sum(cRs, 1));
fprintf('%s: best common alphaF_0 = %s /s (restretch cost %.4g -> %.4g)\n', cond, aLab(aBest), ...
    sum(cRs(:, 1)), sum(cRs(:, aBest)));
fprintf('  gap     peak1 d/m   restretch peak d | a=0 | best    ratio d | a=0 | best\n');
for g = 1:nG
    fprintf('  %6s  %5.1f/%5.1f   %5.1f | %5.1f | %5.1f    %.2f | %.2f | %.2f\n', gapLab{g}, pk1D(g), ...
        pk1M(g, aBest), pkRD(g), pkRM(g, 1), pkRM(g, aBest), pkRD(g)/pk1D(g), ...
        pkRM(g, 1)/pk1M(g, 1), pkRM(g, aBest)/pk1M(g, aBest));
end

% ringing (1.1-2.3 kHz band) from 5 ms before to 20 ms after restretch onset:
% data vs best model seen through the sensor
for g = 1:nG
    S = R.SS{c, g}; kR = numel(S.ev) - 1;
    tf = R.tf{c, g, aBest}; Ff = R.Ff{c, g, aBest};
    w = find(cellfun(@(t) t(1) <= S.ev(kR) + 1e-4 && t(end) >= S.ev(kR), tf), 1, 'last');
    in = S.t >= S.ev(kR) - 5e-3 & S.t <= S.ev(kR) + 20e-3;
    ringD(g) = rms(filtfilt(bb, ab, S.F(in)));
    ringM(g) = rms(filtfilt(bb, ab, interp1(tf{w}, Ff{w}, S.t(in), 'linear', 'extrap')));
end
fprintf('  ringing rms 1.1-2.3 kHz around the restretch, data / model+sensor (kPa):\n');
for g = 1:nG, fprintf('    %6s  %.2f / %.2f\n', gapLab{g}, ringD(g), ringM(g)); end

%% Fig 1: the whole sequence for one slack duration (gap 1 s)
g = 6; S = R.SS{c, g}; kR = numel(S.ev) - 1;
f = figure(901); clf; f.Color = 'w'; f.Position = [40 40 1300 820];
T = tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(T, 1, [1 2]); plot(ax, S.t, S.Lraw, 'k', 'LineWidth', 1.2);
set(ax, 'XLim', [-1 72], 'YLim', [0.75 1.22], 'TickDir', 'out', 'FontSize', 9, 'Box', 'off');
ylabel(ax, 'L / L_0'); title(ax, sprintf('imposed length, slack gap %s', gapLab{g}), 'FontWeight', 'normal');
text(ax, [0.5 30.5 31.5 70.5], [1.2 0.9 1.2 0.84], {'stretch 0.95\rightarrow1.175', 'release to 0.95', ...
    'restretch', 'final release'}, 'FontSize', 8, 'Color', cData, 'VerticalAlignment', 'top');
ax = nexttile(T, 3, [1 2]); hold(ax, 'on');
plot(ax, S.tb, S.Fb, '-', 'Color', cData, 'LineWidth', 1);
plot(ax, S.tb, R.Fb{c, g, 1}, '-', 'Color', cA0, 'LineWidth', 1.4);
plot(ax, S.tb, R.Fb{c, g, aBest}, '-', 'Color', cBest, 'LineWidth', 1.4);
set(ax, 'XLim', [-1 72], 'YLim', [-2 1.1*max(S.Fb)], 'TickDir', 'out', 'FontSize', 9, 'Box', 'off');
xlabel(ax, 'time from first stretch onset (s)'); ylabel(ax, '\Theta (kPa)');
title(ax, 'force, whole sequence (data: bin means)', 'FontWeight', 'normal');
legend(ax, {'data', 'model, no refolding (\alpha_F^0 = 0)', sprintf('model, \\alpha_F^0 = %s s^{-1}', aLab(aBest))}, ...
    'Location', 'northeast', 'Box', 'off');
zk = [1 kR]; zn = {'first stretch', 'restretch'};
for q = 1:2
    ax = nexttile(T, 4 + q); hold(ax, 'on');
    t0 = S.ev(zk(q)); in = S.t >= t0 - 2e-3 & S.t <= t0 + 0.03;
    plot(ax, 1e3*(S.t(in) - t0), S.F(in), '-', 'Color', cData, 'LineWidth', 0.8);
    [tm, y0, yb, yr] = modelAround(R, c, aBest, g, zk(q), 0.03);
    plot(ax, 1e3*tm, y0, '-', 'Color', cA0, 'LineWidth', 1.4);
    plot(ax, 1e3*tm, yb, '-', 'Color', cBest, 'LineWidth', 1.4);
    plot(ax, 1e3*tm, yr, '--', 'Color', cFib, 'LineWidth', 1);
    set(ax, 'XLim', [-2 30], 'YLim', [-2 1.1*max(S.Fb)], 'TickDir', 'out', 'FontSize', 9, 'Box', 'off');
    xlabel(ax, sprintf('time from %s onset (ms)', zn{q})); ylabel(ax, '\Theta (kPa)');
    title(ax, sprintf('zoom: %s', zn{q}), 'FontWeight', 'normal');
    if q == 2
        legend(ax, {'data (10 kHz)', '\alpha_F^0 = 0, through sensor', ...
            sprintf('\\alpha_F^0 = %s, through sensor', aLab(aBest)), ...
            sprintf('\\alpha_F^0 = %s, fibre force (no sensor)', aLab(aBest))}, 'Location', 'northeast', 'Box', 'off');
    end
end
title(T, sprintf('%s: the double-ramp protocol simulated as one sequence (best stretch-hold fit, only \\alpha_F^0 varied)', cond), 'FontSize', 12);
exportgraphics(f, sprintf('RefoldSweep_sequence_%s.png', cond), 'Resolution', 150);

%% Fig 2: restretch for every slack duration
f = figure(902); clf; f.Color = 'w'; f.Position = [60 40 1300 640];
T = tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
yMax = 1.1*max(pk1D);
for g = 1:nG
    S = R.SS{c, g}; kR = numel(S.ev) - 1;
    ax = nexttile(T); hold(ax, 'on');
    i1 = S.t >= S.ev(1) - 2e-3 & S.t <= S.ev(1) + 0.03;
    iR = S.t >= S.ev(kR) - 2e-3 & S.t <= S.ev(kR) + 0.03;
    plot(ax, 1e3*(S.t(i1) - S.ev(1)), sm(S.F(i1), 1e-4), '-', 'Color', cRef, 'LineWidth', 1);
    plot(ax, 1e3*(S.t(iR) - S.ev(kR)), S.F(iR), '-', 'Color', cData, 'LineWidth', 0.8);
    [tm, y0, yb] = modelAround(R, c, aBest, g, kR, 0.03);
    plot(ax, 1e3*tm, y0, '-', 'Color', cA0, 'LineWidth', 1.4);
    plot(ax, 1e3*tm, yb, '-', 'Color', cBest, 'LineWidth', 1.4);
    set(ax, 'XLim', [-2 30], 'YLim', [-2 yMax], 'TickDir', 'out', 'FontSize', 9, 'Box', 'off');
    title(ax, sprintf('slack gap %s', gapLab{g}), 'FontWeight', 'normal');
    if mod(g, 4) == 1, ylabel(ax, '\Theta (kPa)'); else, ax.YTickLabel = []; end
    if g > 4, xlabel(ax, 'time from restretch onset (ms)'); end
end
lg = legend(nexttile(T, 1), {'first stretch, same file (0.7 ms mean)', 'restretch data', ...
    'model, \alpha_F^0 = 0', sprintf('model, \\alpha_F^0 = %s s^{-1}', aLab(aBest))}, ...
    'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 9); lg.Layout.Tile = 'north';
title(T, sprintf('%s: restretch after each slack gap (model seen through the sensor)', cond), 'FontSize', 12);
exportgraphics(f, sprintf('RefoldSweep_restretch_%s.png', cond), 'Resolution', 150);

%% Fig 3: peaks vs slack duration
f = figure(903); clf; f.Color = 'w'; f.Position = [60 60 1350 430];
T = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
xg = R.gap; xg(1) = 1e-3;                   % 0 ms plotted at 1 ms
Q = {pkRD, pkRM, 'restretch peak (kPa)'; pkRD./pk1D, pkRM./pk1M, 'restretch peak / first peak'; ...
     fRD./f1D, fRM./f1M, 'force 1 s into hold: restretch / first'};
for m = 1:3
    ax = nexttile(T); hold(ax, 'on');
    D = Q{m, 1}; M = Q{m, 2};
    for a = 2:nA, plot(ax, xg, M(:, a), '-', 'Color', aCol(a, :), 'LineWidth', 0.8, 'HandleVisibility', 'off'); end
    if m == 1, yline(ax, mean(pk1D), ':', 'first-stretch peak (data)', 'Color', cData, 'FontSize', 7, 'HandleVisibility', 'off'); end
    if m > 1, yline(ax, 1, ':', 'full recovery', 'Color', cData, 'FontSize', 7, 'HandleVisibility', 'off'); end
    plot(ax, xg, M(:, 1), '-', 'Color', cA0, 'LineWidth', 2);
    plot(ax, xg, M(:, aBest), '-', 'Color', cBest, 'LineWidth', 2);
    plot(ax, xg, D, 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'w', 'MarkerSize', 8);
    set(ax, 'XScale', 'log', 'XLim', [7e-4 50], 'XTick', [1e-3 1e-2 1e-1 1 10], ...
        'XTickLabel', {'0', '10 ms', '100 ms', '1 s', '10 s'}, 'TickDir', 'out', 'FontSize', 9, ...
        'YGrid', 'on', 'GridAlpha', 0.12);
    xlabel(ax, 'slack gap at 0.95 L_0'); title(ax, Q{m, 3}, 'FontWeight', 'normal');
end
lg = legend(ax, {'model, \alpha_F^0 = 0', sprintf('model, best \\alpha_F^0 = %s s^{-1}', aLab(aBest)), 'data'}, ...
    'Location', 'southeast', 'FontSize', 8); lg.Box = 'off';
title(T, sprintf(['%s: recovery vs slack duration; thin blue = other \\alpha_F^0 (%s ... %s s^{-1}, ' ...
    'darker = faster refolding)'], cond, aLab(2), aLab(nA)), 'FontSize', 11);
exportgraphics(f, sprintf('RefoldSweep_peaks_%s.png', cond), 'Resolution', 150);

function [t, y0, yb, yr] = modelAround(R, c, aBest, g, k, t1)
% model force around movement k of gap g: fine sensor window, then log bins;
% y0 = alphaF_0 = 0 and yb = aBest through the sensor, yr = aBest fibre force
S = R.SS{c, g}; t0 = S.ev(k);
tf = R.tf{c, g, aBest};
w = find(cellfun(@(t) t(1) <= t0 + 1e-4 && t(end) >= t0, tf), 1, 'last');
iw = tf{w} >= t0 - 2e-3 & tf{w} <= t0 + t1;
bi = S.tb > tf{w}(end) & S.tb <= t0 + t1;
t = [tf{w}(iw); S.tb(bi)] - t0;
y0 = [R.Ff{c, g, 1}{w}(iw); R.Fb{c, g, 1}(bi)];
yb = [R.Ff{c, g, aBest}{w}(iw); R.Fb{c, g, aBest}(bi)];
yr = [R.Fr{c, g, aBest}{w}(iw); nan(nnz(bi), 1)];
end
