%% PlotRefit.m
% ENTRY POINT. Before vs after the joint refit (Model/RefitRefolding.m) on the
% whole double-ramp protocol and the stretch-hold ramps.
%   before   current best stretch-hold fit, no refolding (Model/refoldSweep.mat)
%   sweep    best 2-D sweep point, mechanics fixed (Relax: #1 n-dependent,
%            refoldVariants.mat; Active: attachment cycling / slack release,
%            activeAttach.mat) - both with their fit's own sensor
%   refit    joint refit, mechanics + refolding (+ attachment), shared sensor
%  Fig 1  recovery vs slack gap + the ramps used in the fit
%  Fig 2  restretch after every slack gap
% Workspace config (optional): cond ('Relax' | 'Active').

addpath('../Model');
if ~exist('cond', 'var'), cond = 'Relax'; end
dd = '../Data/2025 11 21 Export/';
F = load(sprintf('../Model/refit_%s_nDep.mat', cond));
R1 = load('../Model/refoldSweep.mat'); c1 = find(strcmp(R1.conds, cond));
rf = F.rf; nG = numel(rf);
gapLab = {'0 ms', '5 ms', '10 ms', '50 ms', '100 ms', '1 s', '10 s', '30 s'};
xg = [1e-3 5e-3 10e-3 50e-3 0.1 1 10 30];
S = F.SP;
for g = 1:nG
    Sr = loadProtocol([dd rf{g} '_refolding_' cond '.txt']); S{g}.t = Sr.t; S{g}.F = Sr.F;
end
cData = [.55 .55 .55]; cRef = [.80 .80 .80];
cB = [42 120 214]/255; cS = [27 175 122]/255; cR = [235 104 52]/255;

% model sets: {Fb, tf, Ff} per gap
mdl = cell(3, 1);
mdl{1} = cellfun(@(k) {R1.Fb{c1, k, 1}, R1.tf{c1, k, 1}, R1.Ff{c1, k, 1}}, num2cell(1:nG), 'UniformOutput', false);
if strcmp(cond, 'Relax')
    V = load('../Model/refoldVariants.mat'); X = V.res.nDep;
    cst = zeros(size(X.Fb, 2), size(X.Fb, 3));
    for g = 1:nG, for a = 1:size(cst, 1), for v = 1:size(cst, 2)
        cst(a, v) = cst(a, v) + rcost(X.Fb{g, a, v}, S{g}); end, end, end
    [~, k] = min(cst(:)); [ia, iv] = ind2sub(size(cst), k);
    mdl{2} = cellfun(@(g) {X.Fb{g, ia, iv}, X.tf{g, ia, iv}, X.Ff{g, ia, iv}}, num2cell(1:nG), 'UniformOutput', false);
    sLab = sprintf('sweep: #1 refolding (\\alpha_F^0 %.3g, \\gamma_F %g), mechanics fixed', X.def.aGrid(ia), X.def.vals(iv));
else
    A = load('../Model/activeAttach.mat'); sz = size(A.Fb); cst = zeros(sz(2:end));
    for j = 1:numel(A.Fb)
        [g, a, d, r] = ind2sub(sz, j); cst(a, d, r) = cst(a, d, r) + rcost(A.Fb{j}, S{g});
    end
    [~, k] = min(cst(:)); [a, d, r] = ind2sub(size(cst), k);
    mdl{2} = cellfun(@(g) {A.Fb{g, a, d, r}, A.tf{g, a, d, r}, A.Ff{g, a, d, r}}, num2cell(1:nG), 'UniformOutput', false);
    sLab = sprintf('sweep: high attachment, k_A %.3g, k_{Dslack} %g, relaxed #1 refolding', A.kAg(a), A.kDsg(d));
end
mdl{3} = cellfun(@(g) {F.FbP{2, g}, F.tfP{2, g}, F.FfP{2, g}}, num2cell(1:nG), 'UniformOutput', false);
labs = {'before: best stretch-hold fit, no refolding', sLab, 'refit: mechanics + refolding jointly, shared sensor'};
cols = {cB, cS, cR};

%% Metrics
sm = @(y, dt) movmean(y, max(1, round(7e-4/dt)));
rD = nan(nG, 1); rM = nan(nG, 3); cM = zeros(1, 3);
for g = 1:nG
    kR = numel(S{g}.ev) - 1; dtd = median(diff(S{g}.t));
    pk = @(t0) max(sm(S{g}.F(S{g}.t >= t0 & S{g}.t <= t0 + 0.02), dtd));
    rD(g) = pk(S{g}.ev(kR))/pk(S{g}.ev(1));
    for m = 1:3
        [Fb, tf, Ff] = mdl{m}{g}{:};
        cM(m) = cM(m) + rcost(Fb, S{g});
        rM(g, m) = peakF(tf, Ff, S{g}.ev(kR))/peakF(tf, Ff, S{g}.ev(1));
    end
end
rampC = @(m) sum(cellfun(@(Fb, Sq) mean(((Fb - Sq.Fb)/max(Sq.Fb)).^2), F.FbR(m, :), F.SR));
fprintf('%s restretch cost: before %.4f | sweep %.4f | refit %.4f;  ramps cost: refit start %.4f -> refit %.4f\n', ...
    cond, cM, rampC(1), rampC(2));
fprintf('  %s\n', F.note);
fprintf('  ratio data  %s\n', sprintf('%.2f ', rD));
for m = 1:3, fprintf('  ratio m%d    %s\n', m, sprintf('%.2f ', rM(:, m))); end
pn = {'Fss', 'n_ss', 'kp', 'np', 'kd', 'nd', 'alphaU', 'nU', 'mu', 'delU', 'kA', 'kD', 'alphaF_0'};
pn(14:28) = {'mu1', 'm_mu', 'Fbeta', 'signedFd', 'clipFd', 'eta', 'wSlack', 'cComp', 'muRec', 'kM', 'tauM', 'kDf', 'F_R', 'kDslack', 'gammaF'};
fprintf('  free parameters, start -> refit:\n');
for i = F.spec.idx, fprintf('    %-9s %10.4g -> %10.4g\n', pn{i}, F.p0(i), F.p(i)); end

%% Fig 1: recovery vs gap + ramps
nR = numel(F.SR);
f = figure(921); clf; f.Color = 'w'; f.Position = [40 60 1450 520];
T = tiledlayout(2, 2 + nR, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(T, 1, [2 2]); hold(ax, 'on');
for m = 1:3, plot(ax, xg, rM(:, m), '-', 'Color', cols{m}, 'LineWidth', 2); end
plot(ax, xg, rD, 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'w', 'MarkerSize', 8);
yline(ax, 1, ':', 'Color', cData, 'HandleVisibility', 'off');
set(ax, 'XScale', 'log', 'XLim', [7e-4 50], 'YLim', [0 1.05], 'XTick', [1e-3 1e-2 1e-1 1 10], ...
    'XTickLabel', {'0', '10 ms', '100 ms', '1 s', '10 s'}, 'TickDir', 'out', 'FontSize', 9, 'YGrid', 'on', 'GridAlpha', 0.12);
xlabel(ax, 'slack gap at 0.95 L_0'); ylabel(ax, 'restretch peak / first peak');
title(ax, sprintf('recovery; restretch cost %.3f | %.3f | %.3f', cM), 'FontWeight', 'normal');
lg = legend(ax, [labs, {'data'}], 'Location', 'southoutside', 'FontSize', 8); lg.Box = 'off';
rampLab = {'first stretch (2.8 ms, 8-repeat pool)', '5.7 ms ramp', '100 ms ramp', '1 s ramp'};
for q = 1:nR
    Sq = F.SR{q};
    for row = 1:2      % top: 0-60 ms linear, bottom: whole hold, log time
        ax = nexttile(T, (row - 1)*(2 + nR) + 2 + q); hold(ax, 'on');
        plot(ax, Sq.tb, Sq.Fb, '.', 'Color', cData, 'MarkerSize', 5);
        plot(ax, Sq.tb, F.FbR{1, q}, '-', 'Color', cB, 'LineWidth', 1.2);
        plot(ax, Sq.tb, F.FbR{2, q}, '-', 'Color', cR, 'LineWidth', 1.2);
        if row == 1
            tEnd = [0.03 0.04 0.2 1.5]; set(ax, 'XLim', [0 tEnd(q)]);
            title(ax, rampLab{q}, 'FontWeight', 'normal');
        else
            set(ax, 'XScale', 'log', 'XLim', [1e-3 30]); xlabel(ax, 'time from onset (s)');
        end
        set(ax, 'TickDir', 'out', 'FontSize', 8, 'Box', 'off');
        if q == 1, ylabel(ax, '\Theta (kPa)'); end
    end
end
title(T, sprintf(['%s: refit on ramps + whole double-ramp; ramps: grey data, blue refit start, ' ...
    'orange refit (ramps cost %.4f \\rightarrow %.4f)'], cond, rampC(1), rampC(2)), 'FontSize', 11);
exportgraphics(f, sprintf('Refit_%s_recovery.png', cond), 'Resolution', 150);

%% Fig 2: restretch after every gap
f = figure(922); clf; f.Color = 'w'; f.Position = [60 40 1400 640];
T = tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
yMax = 1.15*max(cellfun(@(Sg) max(Sg.Fb(Sg.evBin == 1)), S));
for g = 1:nG
    kR = numel(S{g}.ev) - 1; t0 = S{g}.ev(kR);
    ax = nexttile(T); hold(ax, 'on');
    i1 = S{g}.t >= S{g}.ev(1) - 2e-3 & S{g}.t <= S{g}.ev(1) + 0.03;
    iR = S{g}.t >= t0 - 2e-3 & S{g}.t <= t0 + 0.03;
    plot(ax, 1e3*(S{g}.t(i1) - S{g}.ev(1)), sm(S{g}.F(i1), 1e-4), '-', 'Color', cRef, 'LineWidth', 1);
    plot(ax, 1e3*(S{g}.t(iR) - t0), S{g}.F(iR), '-', 'Color', cData, 'LineWidth', 0.8);
    for m = 1:3
        [Fb, tf, Ff] = mdl{m}{g}{:};
        w = find(cellfun(@(t) t(1) <= t0 + 1e-4 && t(end) >= t0, tf), 1, 'last');
        iw = tf{w} >= t0 - 2e-3 & tf{w} <= t0 + 0.03; bi = S{g}.tb > tf{w}(end) & S{g}.tb <= t0 + 0.03;
        plot(ax, 1e3*([tf{w}(iw); S{g}.tb(bi)] - t0), [Ff{w}(iw); Fb(bi)], '-', 'Color', cols{m}, 'LineWidth', 1.3);
    end
    set(ax, 'XLim', [-2 30], 'YLim', [-2 yMax], 'TickDir', 'out', 'FontSize', 9, 'Box', 'off');
    title(ax, sprintf('slack gap %s', gapLab{g}), 'FontWeight', 'normal');
    if mod(g, 4) == 1, ylabel(ax, '\Theta (kPa)'); else, ax.YTickLabel = []; end
    if g > 4, xlabel(ax, 'time from restretch onset (ms)'); end
end
lg = legend(nexttile(T, 1), [{'first stretch, same file', 'restretch data'}, labs], ...
    'NumColumns', 3, 'Box', 'off', 'FontSize', 8); lg.Layout.Tile = 'north';
title(T, sprintf('%s: restretch after each slack gap, before / sweep / refit', cond), 'FontSize', 12);
exportgraphics(f, sprintf('Refit_%s_restretch.png', cond), 'Resolution', 150);

function c = rcost(Fb, S)
kR = numel(S.ev) - 1; rs = S.evBin == kR; sc = max(S.Fb(S.evBin == 1));
c = mean(((Fb(rs) - S.Fb(rs))/sc).^2);
end
function pk = peakF(tf, Ff, t0)
w = find(cellfun(@(t) t(1) <= t0 + 1e-4 && t(end) >= t0, tf), 1, 'last');
pk = max(movmean(Ff{w}(tf{w} >= t0 & tf{w} <= t0 + 0.02), 70));
end
