%% PlotRefoldingVariants.m
% ENTRY POINT. Compares the refolding functions on the whole relaxed
% double-ramp protocol:
%   single rate  alphaF_0, on/off mask       (Model/refoldSweep.mat)
%   #1 nDep      alphaF_0*((n+1)/Ng)^gammaF  (Model/refoldVariants.mat)
%   #3 fGate     alphaF_0*exp(-Fp/F_R)       (Model/refoldVariants.mat)
%  Fig 1  restretch peak / first peak vs slack gap: data, and per value of
%         the second parameter the model at its best alphaF_0.
%  Fig 2  restretch time courses, best of each refolding function.
%  Fig 3  restretch cost maps (alphaF_0 x second parameter).
% Cost = restretch bins only, ((model - data)/first-stretch peak)^2 averaged
% per gap and summed over the 8 gaps (as in PlotRefoldingSweep.m).

addpath('../Model');
dd = '../Data/2025 11 21 Export/';
if ~exist('RV', 'var'), RV = load('../Model/refoldVariants.mat'); end
if ~exist('R1', 'var'), R1 = load('../Model/refoldSweep.mat'); end
rf = RV.rf; nG = numel(rf);
gapLab = {'0 ms', '5 ms', '10 ms', '50 ms', '100 ms', '1 s', '10 s', '30 s'};
xg = [0 5e-3 10e-3 50e-3 0.1 1 10 30]; xg(1) = 1e-3;
S = RV.SS;
for g = 1:nG
    Sr = loadProtocol([dd rf{g} '_refolding_Relax.txt']);
    S{g}.t = Sr.t; S{g}.F = Sr.F;
end
cData = [.55 .55 .55]; cRef = [.80 .80 .80];
cV = [42 120 214; 235 104 52; 27 175 122]/255;   % single, #1, #3
ramp = @(c, n) interp1([0 1], [0.82 + 0.18*c; 0.55*c], linspace(0, 1, n));

%% Data metrics
pk1D = nan(nG, 1); pkRD = pk1D; sm = @(y, dt) movmean(y, max(1, round(7e-4/dt)));
for g = 1:nG
    kR = numel(S{g}.ev) - 1; dtd = median(diff(S{g}.t));
    pkD = @(t0) max(sm(S{g}.F(S{g}.t >= t0 & S{g}.t <= t0 + 0.02), dtd));
    pk1D(g) = pkD(S{g}.ev(1)); pkRD(g) = pkD(S{g}.ev(kR));
end
rD = pkRD./pk1D;

%% Model metrics: M.cost (alpha x val), M.ratio (gap x alpha x val)
c1 = find(strcmp(R1.conds, 'Relax'));
M = struct();
M.single = metrics(squeeze(R1.Fb(c1, :, :)), squeeze(R1.tf(c1, :, :)), squeeze(R1.Ff(c1, :, :)), S);
M.single.def = struct('name', '', 'vals', NaN, 'aGrid', R1.aGrid);
for v = RV.variants
    X = RV.res.(v{1});
    M.(v{1}) = metrics(X.Fb, X.tf, X.Ff, S);
    M.(v{1}).def = X.def;
end
names = {'single', 'nDep', 'fGate'};
titles = {'single rate, on/off mask', '#1 n-dependent  \alpha_F^0((n+1)/N_g)^{\gamma_F}', ...
    '#3 force gate  \alpha_F^0 exp(-F_p/F_R)'};
best = struct();
for m = 1:3
    Q = M.(names{m}); [cmin, k] = min(Q.cost(:)); [ia, iv] = ind2sub(size(Q.cost), k);
    best.(names{m}) = [ia iv];
    fprintf('%-6s best cost %.4f at alphaF_0 = %.3g', names{m}, cmin, Q.def.aGrid(ia));
    if ~isnan(Q.def.vals(1)), fprintf(', %s = %g', Q.def.name, Q.def.vals(iv)); end
    fprintf(' | ratio %s\n', sprintf('%.2f ', Q.ratio(:, ia, iv)));
end
fprintf('data                                         | ratio %s\n', sprintf('%.2f ', rD));

%% Fig 1: recovery vs slack gap
f = figure(911); clf; f.Color = 'w'; f.Position = [40 60 1400 440];
T = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for m = 1:3
    Q = M.(names{m}); nV = numel(Q.def.vals); col = ramp(cV(m, :), max(nV, 2));
    ax = nexttile(T); hold(ax, 'on'); lab = {};
    for iv = 1:nV
        [cmin, ia] = min(Q.cost(:, iv));
        isB = iv == best.(names{m})(2);
        plot(ax, xg, Q.ratio(:, ia, iv), '-', 'Color', col(min(iv, end), :), 'LineWidth', 1 + 1.2*isB);
        if isnan(Q.def.vals(1))
            lab{end+1} = sprintf('\\alpha_F^0 = %.3g, cost %.3f', Q.def.aGrid(ia), cmin); %#ok<SAGROW>
        else
            lab{end+1} = sprintf('%s = %g (\\alpha_F^0 %.3g), cost %.3f', Q.def.name, ...
                Q.def.vals(iv), Q.def.aGrid(ia), cmin); %#ok<SAGROW>
        end
    end
    plot(ax, xg, rD, 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'w', 'MarkerSize', 8);
    yline(ax, 1, ':', 'Color', cData, 'HandleVisibility', 'off');
    set(ax, 'XScale', 'log', 'XLim', [7e-4 50], 'YLim', [0.4 1.05], 'XTick', [1e-3 1e-2 1e-1 1 10], ...
        'XTickLabel', {'0', '10 ms', '100 ms', '1 s', '10 s'}, 'TickDir', 'out', 'FontSize', 9, ...
        'YGrid', 'on', 'GridAlpha', 0.12);
    xlabel(ax, 'slack gap at 0.95 L_0'); if m == 1, ylabel(ax, 'restretch peak / first peak'); end
    title(ax, titles{m}, 'FontWeight', 'normal');
    lg = legend(ax, [lab, {'data'}], 'Location', 'southeast', 'FontSize', 7); lg.Box = 'off';
end
title(T, 'Relaxed: recovery vs slack gap; each line = best \alpha_F^0 for that second-parameter value (thick = overall best)', 'FontSize', 11);
exportgraphics(f, 'RefoldVariants_peaks.png', 'Resolution', 150);

%% Fig 2: restretch time courses, best of each function
f = figure(912); clf; f.Color = 'w'; f.Position = [60 40 1400 640];
T = tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
src = {R1, RV.res.nDep, RV.res.fGate};
for g = 1:nG
    kR = numel(S{g}.ev) - 1; t0 = S{g}.ev(kR);
    ax = nexttile(T); hold(ax, 'on');
    i1 = S{g}.t >= S{g}.ev(1) - 2e-3 & S{g}.t <= S{g}.ev(1) + 0.03;
    iR = S{g}.t >= t0 - 2e-3 & S{g}.t <= t0 + 0.03;
    plot(ax, 1e3*(S{g}.t(i1) - S{g}.ev(1)), sm(S{g}.F(i1), 1e-4), '-', 'Color', cRef, 'LineWidth', 1);
    plot(ax, 1e3*(S{g}.t(iR) - t0), S{g}.F(iR), '-', 'Color', cData, 'LineWidth', 0.8);
    for m = 1:3
        b = best.(names{m});
        if m == 1
            Fb = R1.Fb{c1, g, b(1)}; tf = R1.tf{c1, g, b(1)}; Ff = R1.Ff{c1, g, b(1)};
        else
            Fb = src{m}.Fb{g, b(1), b(2)}; tf = src{m}.tf{g, b(1), b(2)}; Ff = src{m}.Ff{g, b(1), b(2)};
        end
        w = find(cellfun(@(t) t(1) <= t0 + 1e-4 && t(end) >= t0, tf), 1, 'last');
        iw = tf{w} >= t0 - 2e-3 & tf{w} <= t0 + 0.03;
        bi = S{g}.tb > tf{w}(end) & S{g}.tb <= t0 + 0.03;
        plot(ax, 1e3*([tf{w}(iw); S{g}.tb(bi)] - t0), [Ff{w}(iw); Fb(bi)], '-', 'Color', cV(m, :), 'LineWidth', 1.3);
    end
    set(ax, 'XLim', [-2 30], 'YLim', [-2 1.1*max(pk1D)], 'TickDir', 'out', 'FontSize', 9, 'Box', 'off');
    title(ax, sprintf('slack gap %s', gapLab{g}), 'FontWeight', 'normal');
    if mod(g, 4) == 1, ylabel(ax, '\Theta (kPa)'); else, ax.YTickLabel = []; end
    if g > 4, xlabel(ax, 'time from restretch onset (ms)'); end
end
b1 = best.single; b2 = best.nDep; b3 = best.fGate;
lg = legend(nexttile(T, 1), {'first stretch, same file', 'restretch data', ...
    sprintf('single rate (\\alpha_F^0 %.3g)', M.single.def.aGrid(b1(1))), ...
    sprintf('#1 n-dependent (\\alpha_F^0 %.3g, \\gamma_F %g)', M.nDep.def.aGrid(b2(1)), M.nDep.def.vals(b2(2))), ...
    sprintf('#3 force gate (\\alpha_F^0 %.3g, F_R %g)', M.fGate.def.aGrid(b3(1)), M.fGate.def.vals(b3(2)))}, ...
    'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 9); lg.Layout.Tile = 'north';
title(T, 'Relaxed: restretch after each slack gap, best of each refolding function (through the sensor)', 'FontSize', 12);
exportgraphics(f, 'RefoldVariants_restretch.png', 'Resolution', 150);

%% Fig 3: cost maps
f = figure(913); clf; f.Color = 'w'; f.Position = [80 80 1100 420];
T = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
seqBlue = [205 226 251; 158 197 244; 109 167 236; 57 135 229; 37 106 191; 24 79 149; 13 54 107]/255;
cmap = flipud(interp1(linspace(0, 1, 7), seqBlue, linspace(0, 1, 64)));
lim = log10([min([M.nDep.cost(:); M.fGate.cost(:)]), min(M.single.cost(:))*3]);
for m = 2:3
    Q = M.(names{m}); ax = nexttile(T);
    imagesc(ax, log10(Q.cost')); hold(ax, 'on');
    [~, ia] = min(Q.cost, [], 1);
    plot(ax, ia, 1:numel(Q.def.vals), 'w.-', 'MarkerSize', 12);
    plot(ax, best.(names{m})(1), best.(names{m})(2), 'o', 'MarkerFaceColor', cV(m, :), 'MarkerEdgeColor', 'w', 'MarkerSize', 10);
    colormap(ax, cmap); caxis(ax, lim);
    set(ax, 'YDir', 'normal', 'XTick', 1:3:numel(Q.def.aGrid), ...
        'XTickLabel', arrayfun(@(a) sprintf('%.3g', a), Q.def.aGrid(1:3:end), 'UniformOutput', false), ...
        'YTick', 1:numel(Q.def.vals), 'YTickLabel', string(Q.def.vals), 'FontSize', 9);
    xlabel(ax, '\alpha_F^0 (s^{-1})'); ylabel(ax, Q.def.name); title(ax, titles{m}, 'FontWeight', 'normal');
end
cb = colorbar(ax); cb.Label.String = 'log_{10} restretch cost';
title(T, sprintf('Relaxed: restretch cost (single-rate best = %.3f; white = best \\alpha_F^0 per row)', min(M.single.cost(:))), 'FontSize', 11);
exportgraphics(f, 'RefoldVariants_cost.png', 'Resolution', 150);

function Q = metrics(Fb, tf, Ff, S)
% Fb, tf, Ff: gap x alpha (x val) cells; cost summed over gaps, peak ratios
sz = size(Fb); if numel(sz) == 2, sz(3) = 1; end
Q.cost = zeros(sz(2), sz(3)); Q.ratio = nan(sz);
sm = @(y) movmean(y, 70);                           % 0.7 ms on the 1e-5 grid
for g = 1:sz(1)
    kR = numel(S{g}.ev) - 1; rs = S{g}.evBin == kR;
    sc = max(S{g}.Fb(S{g}.evBin == 1));
    for a = 1:sz(2)
        for v = 1:sz(3)
            M = Fb{g, a, v}; t = tf{g, a, v}; F = Ff{g, a, v};
            Q.cost(a, v) = Q.cost(a, v) + mean(((M(rs) - S{g}.Fb(rs))/sc).^2);
            wR = find(cellfun(@(x) x(1) <= S{g}.ev(kR) + 1e-4 && x(end) >= S{g}.ev(kR), t), 1, 'last');
            pk = @(w, t0) max(sm(F{w}(t{w} >= t0 & t{w} <= t0 + 0.02)));
            Q.ratio(g, a, v) = pk(wR, S{g}.ev(kR))/pk(1, S{g}.ev(1));
        end
    end
end
end
