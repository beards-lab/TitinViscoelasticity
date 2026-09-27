%% CompareAttachmentRamps.m
% ENTRY POINT (plotting). High Ca (pCa 4.51), ramps only (2.8 ms pooled first
% stretch, 5.7 ms, 100 ms, 1 s; restretch excluded - it needs refolding).
% Compares the low-attachment fit (joint stage 1, seed 1) with the
% high-attachment fit refined by SweepHighAttach4.m (4 parameters), and shows
% the attached fraction during the 2.8 ms stretch-hold.
% Output: CompareAttachment_highCa_ramps.png

addpath('../Model'); cd ../Model
dd = '..\Data\2025 11 21 Export\'; rf = {'0ms','5ms','10ms','50ms','100ms','1s','10s','30s'};
pad = @(p) [p, NaN(1, 27 - numel(p))]; sRef = 0.2;
Rb = load('fitStretchHold_Relax_best.mat'); Ab = load('fitStretchHold_Active_best.mat');
W = load('jointFit_localtest\sweep4_highAttach.mat');            % refined high attachment
[~, pH, ~, sH] = jointParams(W.th4, W.J);
pRelStiff = pad(Rb.p); pRelStiff([11 12 25 27]) = NaN;          % low attachment (seed 1)
pA1 = pRelStiff; pA1(3) = Ab.p(3)*sRef^Ab.p(4)/sRef^pRelStiff(4); pA1([7 8 14]) = Ab.p([7 8 14]);
pA1([11 12 25 27]) = [0.005*4.1, 0.995*4.1, 0.05, 1];
Jl = W.J; Jl.pRel0 = pRelStiff; Jl.pAct0 = pA1;
ckL = load('jointFit_localtest\ckpt_seed1_stage1.mat'); [~, pL, ~, sL] = jointParams(ckL.th, Jl);
SS = {loadStretchHold(strcat(dd, rf, '_refolding_Active.txt'), [], 200), loadStretchHold([dd '5ms_Ramp_Active.txt'], [], 200), ...
      loadStretchHold([dd '0.1s_Ramp_Active.txt'], [], 200), loadStretchHold([dd '1s_Ramp_Active.txt'], [], 200)};
tit = {'2.8 ms first stretch (8-trial pool)', '5.7 ms ramp', '100 ms ramp', '1 s ramp'};
P = {pH, pL}; Sn = {sH, sL};
Fb = cell(4, 2); Tf = Fb; Ff = Fb;
if isempty(gcp('nocreate')), parpool('Processes'); end
parfor k = 1:8
    [q, m] = ind2sub([4 2], k);
    [fb, tf, ff] = modelBinned(P{m}, SS{q}, Sn{m}); Fb{k} = fb; Tf{k} = tf; Ff{k} = ff;
end
tl_ = [0; logspace(-4, log10(29), 140)'];
[~, ~, faH] = simStretchHold(pH, tl_, SS{1}.ramp); [~, ~, faL] = simStretchHold(pL, tl_, SS{1}.ramp);
cd ../Figures

colH = [0 0.45 0.74]; colL = [0.85 0.33 0.1];
f = figure; f.Position = [60 60 1500 560];
tl = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile; hold on; box on; S = SS{1}; w = S.t <= 8e-3;
plot(1e3*S.t(w), S.F(w), '.', 'Color', [0.7 0.7 0.7]);
plot(1e3*S.tb(S.tb <= 8e-3), S.Fb(S.tb <= 8e-3), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 3);
plot(1e3*Tf{1,1}, Ff{1,1}, '-', 'Color', colH, 'LineWidth', 1.5); plot(1e3*Tf{1,2}, Ff{1,2}, '-', 'Color', colL, 'LineWidth', 1.5);
xlim([0 8]); xlabel('t from ramp onset (ms)'); ylabel('\Theta (kPa)'); title('2.8 ms stretch, first 8 ms');
cH = 0; cL = 0;
for q = 1:4
    nexttile; hold on; box on; S = SS{q};
    plot(S.tb, S.Fb, 'ko', 'MarkerSize', 3);
    plot(S.tb, Fb{q,1}, '-', 'Color', colH, 'LineWidth', 1.5); plot(S.tb, Fb{q,2}, '-', 'Color', colL, 'LineWidth', 1.5);
    set(gca, 'XScale', 'log'); xlim([1e-4 30]); xlabel('t (s)'); ylabel('\Theta (kPa)');
    c1 = mean(((Fb{q,1} - S.Fb)/max(S.Fb)).^2); c2 = mean(((Fb{q,2} - S.Fb)/max(S.Fb)).^2); cH = cH + c1; cL = cL + c2;
    title(sprintf('%s   high %.4f | low %.4f', tit{q}, c1, c2));
    if q == 1, legend('data (log bins)', 'high attachment (refined)', 'low attachment', 'Location', 'northeast'); end
end
nexttile; hold on; box on;
semilogx(tl_(2:end), 100*faH(2:end), '-', 'Color', colH, 'LineWidth', 1.5); semilogx(tl_(2:end), 100*faL(2:end), '-', 'Color', colL, 'LineWidth', 1.5);
set(gca, 'XScale', 'log'); xlim([1e-4 30]); ylim([0 100]); xlabel('t (s)'); ylabel('attached (%)');
title(sprintf('attached fraction, 2.8 ms stretch (high: min %.1f%%)', 100*min(faH)));
title(tl, sprintf('High Ca, ramps only: total cost high attachment %.4f, low attachment %.4f', cH, cL));
exportgraphics(f, 'CompareAttachment_highCa_ramps.png', 'Resolution', 110);
fprintf('high: cost %.5f | fA rest %.5f, min %.4f at %.2f ms, at 29 s %.4f | r %.3g kDf %.3g mu1 %.3g (relaxed %.3g)\n', cH, faH(1), min(faH), 1e3*tl_(find(faH == min(faH), 1)), faH(end), pH(11)+pH(12), pH(25), pH(14), W.J.pRel0(14));
fprintf('low : cost %.5f | fA rest %.4f\n', cL, faL(1));
fprintf('peaks data %.2f | high %.2f | low %.2f ; F(25 s) data %.2f high %.2f low %.2f\n', max(SS{1}.Fb), max(Fb{1,1}), max(Fb{1,2}), interp1(SS{1}.tb, SS{1}.Fb, 25), interp1(SS{1}.tb, Fb{1,1}, 25), interp1(SS{1}.tb, Fb{1,2}, 25));
