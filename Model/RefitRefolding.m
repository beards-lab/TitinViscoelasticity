%% RefitRefolding.m
% ENTRY POINT. Joint refit of the stretch-hold mechanics AND the n-dependent
% refolding (#1: n+1 -> n at alphaF_0*((n+1)/Ng)^gammaF with the on/off mask,
% simStretchHold index 13 and 28) on
%   ramps      first stretch (8-repeat pool), 5.7 ms, 100 ms (active: + 1 s)
%   protocols  whole double-ramp sequence of the 8 slack gaps (loadProtocol)
% Start, relaxed: fitStretchHold_Relax_best.mat + the best (alphaF_0, gammaF)
% of SweepRefoldingVariants.m; free: mechanics set of relaxedMaxwell_setup.
% Start, active: best point of SweepActiveAttachment.m (high attachment with
% slack release / reattachment); free: the same mechanics + kA kD kDf kDslack.
% Both: + log alphaF_0 + gammaF.
% Sensor FIXED at the shared transducer [1560 Hz, zeta 0.09]
% (DataProcessing/AnalyzeRinging.m), so it cannot trade off with mu/eta.
% Workspace config: cond ('Relax' | 'Active'), maxIter (default 15), dryRun
% (true: print the start cost and stop).
% Run as: batch('RefitRefolding', 'Pool', n, 'Workspace', struct('cond', ...));
% checkpoint refit_<cond>_ckpt.mat (a rerun resumes from it), result
% refit_<cond>_nDep.mat (p, sensor, note + model traces for the figures).

if ~exist('cond', 'var'), cond = 'Relax'; end
if ~exist('maxIter', 'var'), maxIter = 15; end
dd = '..\Data\2025 11 21 Export\';
rf = {'0ms', '5ms', '10ms', '50ms', '100ms', '1s', '10s', '30s'};
senFix = [1560 0.09];
wProt = 0.5;                                  % weight of each protocol dataset

% datasets
SR = {loadStretchHold(strcat(dd, rf, ['_refolding_' cond '.txt']), [], 200), ...
      loadStretchHold([dd '5ms_Ramp_' cond '.txt'], [], 200), ...
      loadStretchHold([dd '0.1s_Ramp_' cond '.txt'], [], 200)};
if strcmp(cond, 'Active'), SR{end+1} = loadStretchHold([dd '1s_Ramp_Active.txt'], [], 200); end
SP = cellfun(@(r) loadProtocol([dd r '_refolding_' cond '.txt']), rf, 'UniformOutput', false);
SS = [SR, SP];

% start point and free set
if strcmp(cond, 'Relax')
    B = load('fitStretchHold_Relax_best.mat');
    U = load('relaxedMaxwell_setup.mat'); spec = U.specM; lb = U.lbM(1:end-2); ub = U.ubM(1:end-2);
    [aS, gS] = bestNDep(load('refoldVariants.mat'));
    p0 = [B.p, NaN(1, 28 - numel(B.p))]; p0(13) = aS; p0(28) = gS;
    odeOpts = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
else
    % high-attachment model with slack release / reattachment: best point of
    % SweepActiveAttachment.m; mechanics set as relaxed + kA kD kDf kDslack
    U = load('relaxedMaxwell_setup.mat'); spec = U.specM; lb = U.lbM(1:end-2); ub = U.ubM(1:end-2);
    spec.idx = [spec.idx, 11 12 25 27]; spec.islog = [spec.islog, true(1, 4)];
    lb = [lb, log([1e-5 1e-5 1e-4 1e-2])]; ub = [ub, log([1e7 1e7 1e3 1e6])];
    p0 = bestAttach(load('activeAttach.mat'));
    if ~(p0(13) > 0)                          % refolding off at the sweep best:
        [p0(13), p0(28)] = bestNDep(load('refoldVariants.mat'));  % start at relaxed #1
    end
    p0(27) = max(p0(27), 1e-2);
    aS = p0(13); gS = p0(28);
    odeOpts = odeset('RelTol', 1e-5, 'AbsTol', 1e-5);
end
spec.pBase = p0; spec.sensor = senFix; spec.odeOpts = odeOpts;
spec.idx = [spec.idx, 13, 28]; spec.islog = [spec.islog, true, false];
spec.wts = [ones(1, numel(SR)), wProt*ones(1, numel(SP))];
lb = [lb, log(1e-3), 0]; ub = [ub, log(1e5), 12];
th = resVariant([], spec, SS); senTh = th(end-1:end); th0 = th(1:end-2);
th0 = min(max(th0, lb + 1e-6), ub - 1e-6);
res = @(x) resVariant([x, senTh], spec, SS);

ckf = sprintf('refit_%s_ckpt.mat', cond);
if exist(ckf, 'file')
    ck = load(ckf); x0 = ck.x;
    fprintf('resuming from checkpoint (iteration %d, resnorm %.5f)\n', ck.it, ck.res);
else
    x0 = th0;
end
r0 = res(th0); c0 = sum(r0.^2);
fprintf('%s start: alphaF_0 %.3g, gammaF %g, cost %.5f (sensor fixed [%g %g])\n', cond, aS, gS, c0, senFix);
if exist('dryRun', 'var') && dryRun, return; end   % start cost only
opt = optimoptions('lsqnonlin', 'Display', 'iter', 'UseParallel', true, ...
    'FiniteDifferenceStepSize', 1e-2, 'MaxIterations', maxIter, ...
    'FunctionTolerance', 1e-6, 'StepTolerance', 1e-5, ...
    'OutputFcn', @(x, ov, st) ckptSave(x, ov, st, ckf));
[x, rn, ~, flag, out] = lsqnonlin(res, x0, lb, ub, opt);
[~, p, sensor] = resVariant([x, senTh], spec, SS);
fprintf('%s refit: cost %.5f -> %.5f, flag %d, %d iterations; alphaF_0 %.3g, gammaF %.3g\n', ...
    cond, c0, rn, flag, out.iterations, p(13), p(28));

% model traces for the figures (start and refit, same fixed sensor)
P = {p0, p};
FbR = cell(2, numel(SR)); FbP = cell(2, numel(SP)); tfP = FbP; FfP = FbP; FrP = FbP;
for m = 1:2
    for q = 1:numel(SR), FbR{m, q} = modelBinned(P{m}, SR{q}, sensor, odeOpts); end
    parfor q = 1:numel(SP)
        [FbP{m, q}, tfP{m, q}, FfP{m, q}, FrP{m, q}] = modelBinned(P{m}, SP{q}, sensor, odeOpts);
    end
end
for q = 1:numel(SR), SR{q} = rmfield(SR{q}, {'A', 'tSim', 'ramp', 't', 'F', 'Lraw'}); end
for q = 1:numel(SP), SP{q} = rmfield(SP{q}, {'A', 'tSim', 'ramp', 't', 'F', 'Lraw'}); end
note = sprintf(['%s, joint refit mechanics + n-dependent refolding (#1) on ramps + 8 double-ramp ' ...
    'protocols (weight %.2g each), sensor fixed [%g %g], cost %.5f -> %.5f'], cond, wProt, senFix, c0, rn);
save(sprintf('refit_%s_nDep.mat', cond), 'p', 'p0', 'sensor', 'note', 'x', 'th0', 'rn', 'c0', 'flag', ...
    'out', 'spec', 'FbR', 'FbP', 'tfP', 'FfP', 'FrP', 'SR', 'SP', 'rf');

function [a, g] = bestNDep(V)
% best (alphaF_0, gammaF) of the n-dependent sweep: restretch cost over gaps
X = V.res.nDep; S = V.SS;
cost = zeros(numel(X.def.aGrid), numel(X.def.vals));
for gi = 1:numel(S)
    kR = numel(S{gi}.ev) - 1; rs = S{gi}.evBin == kR; sc = max(S{gi}.Fb(S{gi}.evBin == 1));
    for ai = 1:size(cost, 1)
        for vi = 1:size(cost, 2)
            cost(ai, vi) = cost(ai, vi) + mean(((X.Fb{gi, ai, vi}(rs) - S{gi}.Fb(rs))/sc).^2);
        end
    end
end
[~, k] = min(cost(:)); [ia, iv] = ind2sub(size(cost), k);
a = X.def.aGrid(ia); g = X.def.vals(iv);
end

function p = bestAttach(A)
% best point of SweepActiveAttachment.m (restretch cost over the gaps)
S = A.SS; sz = size(A.Fb); cost = zeros(sz(2:end));
for j = 1:numel(A.Fb)
    [gi, ai, di, ri] = ind2sub(sz, j);
    kR = numel(S{gi}.ev) - 1; rs = S{gi}.evBin == kR; sc = max(S{gi}.Fb(S{gi}.evBin == 1));
    cost(ai, di, ri) = cost(ai, di, ri) + mean(((A.Fb{j}(rs) - S{gi}.Fb(rs))/sc).^2);
end
[~, k] = min(cost(:)); [ai, di, ri] = ind2sub(size(cost), k);
p = A.p0; p(11) = A.kAg(ai); p(12) = A.kAg(ai)*(1 - A.fA)/A.fA;
p(27) = A.kDsg(di); p(13) = A.refold(ri, 1); p(28) = A.refold(ri, 2);
end

function stop = ckptSave(x, ov, st, f)
stop = false;
if strcmp(st, 'iter'), it = ov.iteration; res = ov.resnorm; save(f, 'x', 'it', 'res'); end
end
