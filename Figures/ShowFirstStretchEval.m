%% ShowFirstStretchEval.m
% ENTRY POINT (plotting). Shows how a parameter set is evaluated against the
% relaxed first stretch-hold (see Model/FitFirstStretch.m):
%   A  raw 10 kHz samples (8-trial ensemble), log-bin means, model force
%      (true) and model seen through the 2nd-order force transducer
%   B  the same on log time over the whole 30 s hold
%   C  normalized bin residuals - each bin is one term of the cost
%   D  the 5.7 ms and 100 ms ramps of the same preparation (also in the cost)
% Workspace: pEval (params), senEval ([f0 zeta]), SSe (loadStretchHold sets),
% cond ('Relax' | 'Active', optional, for titles/file names).

addpath('../Model');
if ~exist('cond', 'var'), cond = 'Relax'; end
q = 1; S = SSe{q};
[Fb, tf, Ff] = modelBinned(pEval, S, senEval);
Ftrue = simStretchHold(pEval, tf, S.ramp);
[cost, parts, Fm] = costStretchHold(pEval, SSe, [], senEval);

f = figure(772); clf; f.Position = [80 80 1350 720];
tl = tiledlayout(2, 3, 'TileSpacing', 'compact');

% A: fast zoom
nexttile(1, [1 2]); hold on; box on;
w = S.t <= 8e-3;
plot(1e3*S.t(w), S.F(w), '.', 'Color', [0.65 0.65 0.65], 'MarkerSize', 8);
wb = S.tb <= 8e-3;
plot(1e3*S.tb(wb), S.Fb(wb), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
plot(1e3*tf, Ftrue, '--', 'Color', [0 0.45 0.74], 'LineWidth', 1.2);
plot(1e3*tf, Ff, '-', 'Color', [0 0.45 0.74], 'LineWidth', 1.8);
plot(1e3*S.tb(wb), Fb(wb), 's', 'Color', [0.85 0.33 0.1], 'MarkerSize', 6, 'LineWidth', 1.2);
ylabel('\Theta - baseline (kPa)');
yyaxis right; plot(1e3*S.t(w), S.Lraw(w), 'k:'); ylabel('L (L_0)'); yyaxis left;
xlim([0 8]); xlabel('t from ramp onset (ms)');
legend('data: 8-trial mean, 10 kHz samples', 'data: log-bin mean', 'model force', ...
    sprintf('model through sensor (f_0 %.0f Hz, \\zeta %.2f)', senEval), 'model: log-bin mean', 'length', ...
    'Location', 'northeast');
title('A  2.8 ms stretch, ramp and ringing');

% B: log time, whole hold
nexttile(3); hold on; box on;
plot(S.tb, S.Fb, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 3);
plot(S.tb, Fb, '-', 'Color', [0.85 0.33 0.1], 'LineWidth', 1.5);
set(gca, 'XScale', 'log'); xlim([1e-4 30]); xlabel('t (s)'); ylabel('\Theta (kPa)');
legend('data bins', 'model bins'); title('B  log bins over the 30 s hold');

% C: residuals per bin
nexttile(4, [1 2]); hold on; box on;
res = (Fb - S.Fb)/max(S.Fb);
stem(S.tb, res, 'filled', 'MarkerSize', 3, 'Color', [0.85 0.33 0.1]);
yline(0, 'k-');
yline([-1 1]*S.noise/sqrt(max(1, median(S.nb)))/max(S.Fb), 'k:');
set(gca, 'XScale', 'log'); xlim([1e-4 30]);
xlabel('t (s)'); ylabel('(model - data)/max(data)');
title(sprintf('C  bin residuals: cost_{2.8ms} = mean(res^2) = %.5f  (%d bins, ~%d bins/decade)', ...
    parts(1, 1), numel(S.Fb), round(numel(S.Fb)/log10(30/S.tb(1)))));

% D: other ramps
nexttile(6); hold on; box on;
cl = [0.47 0.67 0.19; 0.49 0.18 0.56];
for k = 2:3
    plot(SSe{k}.tb, SSe{k}.Fb, 'o', 'Color', cl(k-1, :), 'MarkerSize', 3);
    plot(SSe{k}.tb, Fm{k}, '-', 'Color', cl(k-1, :), 'LineWidth', 1.5);
end
set(gca, 'XScale', 'log'); xlim([1e-4 30]); xlabel('t (s)'); ylabel('\Theta (kPa)');
legend('5.7 ms data', 'model', '100 ms data', 'model', 'Location', 'northeast');
title(sprintf('D  cost %.5f / %.5f', parts(2, 1), parts(3, 1)));

title(tl, sprintf('%s first stretch-hold, cost (3 sets) %.5f   params %s', cond, cost, mat2str(pEval(1:10), 3)), 'FontSize', 10);
exportgraphics(f, sprintf('FirstStretchEval_%s.png', cond), 'Resolution', 110);
