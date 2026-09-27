%% CompareAttachment.m
% ENTRY POINT (plotting). High- vs low-attachment hypotheses against the
% 2025-11-21 stretch-hold data, one figure per Ca level:
%   low Ca  (pCa 11):  relaxed mechanics each hypothesis needs
%                      (high attachment: soft distal element; low: stiff)
%   high Ca (pCa 4.51): stage-1 joint fits from OptimizeStretchHoldJoint.m
%                      (Model/jointFit_localtest/ckpt_seed{3,1}_stage1.mat)
% Panels: first 8 ms of the 2.8 ms stretch, then every dataset on log time.
% Output: CompareAttachment_lowCa.png, CompareAttachment_highCa.png

addpath('../Model');
cd ../Model
dd = '..\Data\2025 11 21 Export\';
rf = {'0ms', '5ms', '10ms', '50ms', '100ms', '1s', '10s', '30s'};
pad = @(p) [p, NaN(1, 27 - numel(p))];
sRef = 0.2;

%% parameter sets
Rb = load('fitStretchHold_Relax_best.mat');       % stiff distal (low attachment)
Rs = load('fitStretchHold_Relax_softDistal.mat');  % soft distal (high attachment)
Ab = load('fitStretchHold_Active_best.mat');
pRelStiff = pad(Rb.p); pRelStiff([11 12 25 27]) = NaN;
pRelSoft  = pad(Rs.p); pRelSoft([11 12 25 27]) = NaN;
J = struct('sharedIdx', [9 3 4 5 6 7 8 1 2 10 14 15 19 23 24], ...
    'sharedLog', logical([1 1 0 1 0 1 0 0 0 0 1 0 1 1 1]), 'refKp', true, 'refAU', true, ...
    'caNames', {{'kpX', 'nU', 'aUX', 'mu1X', 'fA', 'r', 'kDf', 'kDslack'}}, ...
    'senRel', Rb.sensor, 'senAct', Ab.sensor);
% high attachment (seed 3)
Jh = J; Jh.pRel0 = pRelSoft; Jh.pAct0 = setAct(pRelSoft, 0.95, 1000, 0.3, 3000);
ck = load('jointFit_localtest\ckpt_seed3_stage1.mat');
[pRelH, pActH, senRelH, senActH] = jointParams(ck.th, Jh);
% low attachment (seed 1)
pA1 = pRelStiff; pA1(3) = Ab.p(3)*sRef^Ab.p(4)/sRef^pRelStiff(4); pA1([7 8 14]) = Ab.p([7 8 14]);
Jl = J; Jl.pRel0 = pRelStiff; Jl.pAct0 = setAct(pA1, 0.005, 4.1, 0.05, 1);
ck = load('jointFit_localtest\ckpt_seed1_stage1.mat');
[pRelL, pActL, senRelL, senActL] = jointParams(ck.th, Jl);

%% data
SSrel = {loadStretchHold(strcat(dd, rf, '_refolding_Relax.txt'), [], 200), ...
         loadStretchHold([dd '5ms_Ramp_Relax.txt'], [], 200), ...
         loadStretchHold([dd '0.1s_Ramp_Relax.txt'], [], 200), ...
         loadRestretch([dd '0ms_refolding_Relax.txt'])};
SSact = {loadStretchHold(strcat(dd, rf, '_refolding_Active.txt'), [], 200), ...
         loadStretchHold([dd '5ms_Ramp_Active.txt'], [], 200), ...
         loadStretchHold([dd '0.1s_Ramp_Active.txt'], [], 200), ...
         loadStretchHold([dd '1s_Ramp_Active.txt'], [], 200), ...
         loadRestretch([dd '0ms_refolding_Active.txt'])};
titRel = {'2.8 ms first stretch', '5.7 ms ramp', '100 ms ramp', 'restretch (0 ms gap)'};
titAct = {'2.8 ms first stretch', '5.7 ms ramp', '100 ms ramp', '1 s ramp', 'restretch (0 ms gap)'};

%% model evaluations (in parallel)
jobs = {};
for q = 1:4, jobs(end+1, :) = {SSrel{q}, pRelH, senRelH}; jobs(end+1, :) = {SSrel{q}, pRelL, senRelL}; end %#ok<SAGROW>
for q = 1:5, jobs(end+1, :) = {SSact{q}, pActH, senActH}; jobs(end+1, :) = {SSact{q}, pActL, senActL}; end %#ok<SAGROW>
Fb = cell(size(jobs, 1), 1); Tf = Fb; Ff = Fb;
if isempty(gcp('nocreate')), parpool('Processes'); end
parfor k = 1:size(jobs, 1)
    [Fb{k}, tf, ff] = modelBinned(jobs{k, 2}, jobs{k, 1}, jobs{k, 3});
    Tf{k} = tf; Ff{k} = ff;
end
cd ../Figures

%% figures
nr = 4; na = 5;
drawCond(SSrel, titRel, Fb(1:2*nr), Tf(1:2*nr), Ff(1:2*nr), 'Low Ca (pCa 11)', 'CompareAttachment_lowCa.png');
drawCond(SSact, titAct, Fb(2*nr+1:end), Tf(2*nr+1:end), Ff(2*nr+1:end), 'High Ca (pCa 4.51)', 'CompareAttachment_highCa.png');
fprintf('attached at rest: high %.3f, low %.4f\n', pActH(11)/(pActH(11)+pActH(12)), pActL(11)/(pActL(11)+pActL(12)));

function drawCond(SS, tit, Fb, Tf, Ff, figTitle, fname)
colH = [0 0.45 0.74]; colL = [0.85 0.33 0.1];
n = numel(SS);
cost = @(k, q) mean(((Fb{k} - SS{q}.Fb)/max(SS{q}.Fb)).^2);
f = figure; f.Position = [60 60 1500 560];
tl = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
% fast window of the first stretch
nexttile; hold on; box on; S = SS{1}; w = S.t <= 8e-3;
plot(1e3*S.t(w), S.F(w), '.', 'Color', [0.7 0.7 0.7]);
plot(1e3*S.tb(S.tb <= 8e-3), S.Fb(S.tb <= 8e-3), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 3);
plot(1e3*Tf{1}, Ff{1}, '-', 'Color', colH, 'LineWidth', 1.5);
plot(1e3*Tf{2}, Ff{2}, '-', 'Color', colL, 'LineWidth', 1.5);
xlim([0 8]); xlabel('t from ramp onset (ms)'); ylabel('\Theta (kPa)'); title('2.8 ms stretch, first 8 ms');
ch = 0; cl = 0;
for q = 1:n
    nexttile; hold on; box on; S = SS{q}; pos = S.tb > 0;
    plot(S.tb(pos), S.Fb(pos), 'ko', 'MarkerSize', 3);
    plot(S.tb(pos), Fb{2*q-1}(pos), '-', 'Color', colH, 'LineWidth', 1.5);
    plot(S.tb(pos), Fb{2*q}(pos), '-', 'Color', colL, 'LineWidth', 1.5);
    set(gca, 'XScale', 'log'); xlim([1e-4 30]); xlabel('t (s)'); ylabel('\Theta (kPa)');
    c1 = cost(2*q-1, q); c2 = cost(2*q, q); ch = ch + c1; cl = cl + c2;
    title(sprintf('%s   high %.4f | low %.4f', tit{q}, c1, c2));
    if q == 1, legend('data (log bins)', 'high attachment', 'low attachment', 'Location', 'northeast'); end
end
title(tl, sprintf('%s: total cost high attachment %.4f, low attachment %.4f', figTitle, ch, cl));
exportgraphics(f, fname, 'Resolution', 110);
end

function p = setAct(p, fA, r, kDf, kDs)
p(11) = fA*r; p(12) = (1 - fA)*r; p(25) = kDf; p(27) = kDs;
end
