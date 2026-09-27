%% OptimizeStretchHoldJoint.m
% ENTRY POINT. Joint low/high-Ca fit of the 2025-11-21 stretch-hold data, for
% a many-core desktop. Relaxed (pCa 11) and active (pCa 4.51) share all
% mechanical parameters; only the Ca-specific set differs, as in the
% published fit (kp, alphaU, nU, kA) plus the attachment kinetics needed by
% the new data (attached fraction, cycling rate, slip-bond and slack
% detachment). Relaxed data alone cannot separate the proximal and distal
% chain stiffness; the attached (active) chains load the distal element
% directly, which pins it down - hence one joint fit rather than two.
%
% Datasets (10 kHz exports in ..\Data\2025 11 21 Export\):
%   relaxed: first stretch (8-trial ensemble of *_refolding_Relax.txt),
%            5ms_Ramp_Relax, 0.1s_Ramp_Relax, restretch of 0ms_refolding_Relax
%   active:  first stretch (8-trial ensemble), 5ms_Ramp_Active,
%            0.1s_Ramp_Active, 1s_Ramp_Active, restretch of 0ms_refolding_Active
% Model and cost: see FitFirstStretch.m (simStretchHold.m with dXdTvar.m,
%   sensorFilter.m, log bins, cost = sum_q w_q*mean(((model-data)/max(data))^2)).
%
% Seeds (multi-start over the attachment hypothesis):
%   1 lowAttach   ~0.5 % attached at rest (as the published fit, ~1.4 %),
%                 Ca-specific kp, alphaU, nU, mu1 from the independent active fit
%   2 midAttach    50 % attached, cycling 100/s
%   3 highAttach   95 % attached, cycling 1000/s, soft distal element
% Each seed: stage 1 fits the Ca-specific set + active sensor with the shared
% mechanics fixed, stage 2 fits everything jointly.
%
% HOW TO RUN (MATLAB R2023a+, Optimization + Parallel Computing Toolboxes):
%   cd Model; nWorkers = 30; OptimizeStretchHoldJoint
%   - one residual evaluation (9 datasets, ODE tol 1e-6) takes ~3-5 min on
%     one core; lsqnonlin evaluates nFree+1 of them per iteration in
%     parallel, so ~30 workers give ~1 evaluation time per iteration
%     (~2-3 h per seed for 40 iterations). The script prints its own estimate.
%   - dryRun = true only evaluates and prints the seed costs (and timing).
%   - it checkpoints every iteration to outDir; re-running resumes from the
%     latest checkpoint of each stage (resume = true).
% Workspace config (defaults below): nWorkers, seedsToRun, maxIterWarm,
%   maxIterJoint, odeTol, wRel, wAct, caNames, outDir, resume, dryRun, writeBest.

if ~exist('nWorkers', 'var'),     nWorkers = max(1, feature('numcores') - 1); end
if ~exist('seedsToRun', 'var'),   seedsToRun = 1:3; end
if ~exist('maxIterWarm', 'var'),  maxIterWarm = 10; end
if ~exist('maxIterJoint', 'var'), maxIterJoint = 40; end
if ~exist('odeTol', 'var'),       odeTol = 1e-6; end   % 1e-4 is too noisy for FD gradients
if ~exist('wRel', 'var'),         wRel = [1 1 1 1]; end
if ~exist('wAct', 'var'),         wAct = [1 1 1 1 1]; end
if ~exist('caNames', 'var'),      caNames = {'kpX', 'nU', 'aUX', 'mu1X', 'fA', 'r', 'kDf', 'kDslack'}; end
if ~exist('outDir', 'var'),       outDir = 'jointFit_results'; end
if ~exist('resume', 'var'),       resume = true; end
if ~exist('dryRun', 'var'),       dryRun = false; end
if ~exist('writeBest', 'var'),    writeBest = true; end  % write fitStretchHold_*_joint.mat
if ~exist(outDir, 'dir'), mkdir(outDir); end

%% Data
dd = '..\Data\2025 11 21 Export\';
rf = {'0ms', '5ms', '10ms', '50ms', '100ms', '1s', '10s', '30s'};
SSrel = {loadStretchHold(strcat(dd, rf, '_refolding_Relax.txt'), [], 200), ...
         loadStretchHold([dd '5ms_Ramp_Relax.txt'], [], 200), ...
         loadStretchHold([dd '0.1s_Ramp_Relax.txt'], [], 200), ...
         loadRestretch([dd '0ms_refolding_Relax.txt'])};
SSact = {loadStretchHold(strcat(dd, rf, '_refolding_Active.txt'), [], 200), ...
         loadStretchHold([dd '5ms_Ramp_Active.txt'], [], 200), ...
         loadStretchHold([dd '0.1s_Ramp_Active.txt'], [], 200), ...
         loadStretchHold([dd '1s_Ramp_Active.txt'], [], 200), ...
         loadRestretch([dd '0ms_refolding_Active.txt'])};
namesRel = {'rel 2.8ms', 'rel 5.7ms', 'rel 100ms', 'rel restretch'};
namesAct = {'act 2.8ms', 'act 5.7ms', 'act 100ms', 'act 1s', 'act restretch'};

%% Parameter layout (params indices as in simStretchHold.m, 27 entries)
%            mu0 Fp0  np  kd  nd  U0  nU  Fss n_ss delU mu1  m  eta kM tauM
sharedIdx = [  9   3   4   5   6   7   8   1   2   10  14  15  19  23  24];
sharedLog = logical([1 1 0 1 0 1 0 0 0 0 1 0 1 1 1]);
lbSh = [log(1e-4) log(0.1) 1 log(10)  1 log(1e-2)  2 0.5  2 0.05 log(1e-3)  1 log(1e-5) log(0.01) log(1e-4)];
ubSh = [log(1)    log(100) 9 log(1e9) 9 log(1e6)  25 6   30 0.5  log(30)   30 log(1)    log(100)  log(1)];
caLb = struct('kpX', log(0.01), 'nU', 2,  'aUX', log(1e-4), 'mu1X', log(1e-2), 'fA', -12, 'r', log(1e-3), 'kDf', log(1e-4), 'kDslack', log(1e-2));
caUb = struct('kpX', log(1e3),  'nU', 25, 'aUX', log(1e4),  'mu1X', log(1e3),  'fA',  12, 'r', log(1e6),  'kDf', log(100),  'kDslack', log(1e6));
senLb = [log(800) 0.05]; senUb = [log(4000) 1];
lb = [lbSh, cellfun(@(n) caLb.(n), caNames), senLb, senLb];
ub = [ubSh, cellfun(@(n) caUb.(n), caNames), senUb, senUb];
nSh = numel(sharedIdx); nCa = numel(caNames);

%% Seeds
Rb = load('fitStretchHold_Relax_best.mat');            % stiff distal element
Ab = load('fitStretchHold_Active_best.mat');           % independent active fit
pad = @(p) [p, NaN(1, 27 - numel(p))];
pRelStiff = pad(Rb.p); pRelStiff([11 12 25 27]) = NaN;
if exist('fitStretchHold_Relax_softDistal.mat', 'file')
    Rs = load('fitStretchHold_Relax_softDistal.mat'); pRelSoft = pad(Rs.p);  % relaxed refit, soft distal
else
    pRelSoft = pRelStiff; pRelSoft(5) = 40/0.225^pRelStiff(6);
end
pRelSoft([11 12 25 27]) = NaN;
% low attachment: Ca-specific set from the independent active fit, transferred
% at the reference strain (that fit has its own np), so the folded-chain
% force and unfolding rate at s = 0.2 um match it
sRef = 0.2; pA1 = pRelStiff;
pA1(3) = Ab.p(3)*sRef^Ab.p(4)/sRef^pRelStiff(4);         % same Fp0(sRef), shared np
pA1([7 8 14]) = Ab.p([7 8 14]);
pA1 = mkAct(pA1, 0.005, 4.1, 0.05, 1);
seeds = struct('name', {'lowAttach', 'midAttach', 'highAttach'}, ...
    'pRel', {pRelStiff, pRelSoft, pRelSoft}, ...
    'pAct', {pA1, mkAct(pRelSoft, 0.5, 100, 0.1, 1000), mkAct(pRelSoft, 0.95, 1000, 0.3, 3000)}, ...
    'senRel', {Rb.sensor, Rb.sensor, Rb.sensor}, 'senAct', {Ab.sensor, Ab.sensor, Ab.sensor});

J0 = struct('sharedIdx', sharedIdx, 'sharedLog', sharedLog, 'refKp', true, 'refAU', true, ...
    'caNames', {caNames}, 'wRel', wRel, 'wAct', wAct, 'odeOpts', odeset('RelTol', odeTol, 'AbsTol', odeTol));

%% Pool and timing
if isempty(gcp('nocreate')), parpool('Processes', nWorkers); end
pool = gcp;
fprintf('%s: %d workers, seeds %s\n', datetime('now'), pool.NumWorkers, mat2str(seedsToRun));

results = struct([]);
for si = seedsToRun
    sd = seeds(si);
    J = J0; J.pRel0 = sd.pRel; J.pAct0 = sd.pAct; J.senRel = sd.senRel; J.senAct = sd.senAct;
    th0 = jointParams([], J);
    th0 = min(max(th0, lb + 1e-6), ub - 1e-6);
    J.thetaFix = th0;
    J.free = true(size(th0));
    t1 = tic; r0 = resJoint(th0, J, SSrel, SSact); tEval = toc(t1);
    fprintf('\n== seed %d %s: initial cost %.5f (one evaluation %.0f s; ~%.0f min per joint iteration)\n', ...
        si, sd.name, sum(r0.^2), tEval, ceil((numel(th0) + 1)/pool.NumWorkers)*tEval/60);
    if dryRun, continue; end

    % stage 1: Ca-specific + active sensor, shared mechanics fixed
    free1 = false(size(th0)); free1(nSh+1:nSh+nCa) = true; free1(end-1:end) = true;
    J1 = J; J1.wRel(:) = 0;              % relaxed predictions do not change in stage 1
    th1 = runStage(th0, free1, J1, SSrel, SSact, lb, ub, maxIterWarm, ...
        fullfile(outDir, sprintf('ckpt_seed%d_stage1.mat', si)), resume);
    % stage 2: everything
    th2 = runStage(th1, true(size(th0)), J, SSrel, SSact, lb, ub, maxIterJoint, ...
        fullfile(outDir, sprintf('ckpt_seed%d_stage2.mat', si)), resume);

    [pRel, pAct, senRel, senAct] = jointParams(th2, J);
    J.thetaFix = th2; J.free = true(size(th2));
    r = resJoint(th2, J, SSrel, SSact);
    nbs = cellfun(@(S) numel(S.Fb), [SSrel, SSact]); ce = cumsum([0 nbs]);
    parts = arrayfun(@(q) sum(r(ce(q)+1:ce(q+1)).^2), 1:numel(nbs));
    res = struct('seed', sd.name, 'theta', th2, 'pRel', pRel, 'pAct', pAct, 'senRel', senRel, ...
        'senAct', senAct, 'cost', sum(r.^2), 'parts', parts, 'fA', pAct(11)/(pAct(11) + pAct(12)), 'J', J);
    save(fullfile(outDir, sprintf('result_seed%d_%s.mat', si, sd.name)), '-struct', 'res');
    results = [results, res]; %#ok<AGROW>
    fprintf('seed %d %s: cost %.5f, resting attached fraction %.3f\n', si, sd.name, res.cost, res.fA);
    fprintf('   %s\n', strjoin(compose('%s %.4f', string([namesRel, namesAct])', parts'), ', '));
end

%% Best seed -> result files readable by FitFirstStretch.m
if ~isempty(results) && writeBest
    [~, ib] = min([results.cost]); b = results(ib);
    p = b.pRel; sensor = b.senRel; note = sprintf('joint fit, seed %s, total cost %.5f (relaxed part)', b.seed, b.cost);
    save('fitStretchHold_Relax_joint.mat', 'p', 'sensor', 'note');
    p = b.pAct; sensor = b.senAct; note = sprintf('joint fit, seed %s, total cost %.5f, attached at rest %.3f', b.seed, b.cost, b.fA);
    save('fitStretchHold_Active_joint.mat', 'p', 'sensor', 'note');
    fprintf('\nBest: seed %s, cost %.5f. Plot with: cond = ''Active''; resultFile = ''fitStretchHold_Active_joint.mat''; FitFirstStretch\n', b.seed, b.cost);
    for k = 1:numel(results)
        fprintf('  %-10s cost %.5f  attached at rest %.3f\n', results(k).seed, results(k).cost, results(k).fA);
    end
end

%% ---------------------------------------------------------------------
function p = mkAct(p, fA, r, kDf, kDs)
% Active seed from a relaxed set: resting attached fraction fA, cycling rate
% r = kA + kD, slip-bond sensitivity kDf, slack detachment rate kDs.
p(11) = fA*r; p(12) = (1 - fA)*r; p(25) = kDf; p(27) = kDs;
end

function th = runStage(th, free, J, SSrel, SSact, lb, ub, maxIter, ckpt, resume)
% One lsqnonlin stage on the free subset of theta, resuming from ckpt.
if resume && exist(ckpt, 'file')
    ck = load(ckpt);
    if isfield(ck, 'done') && ck.done
        fprintf('   %s already finished (cost %.5f)\n', ckpt, ck.resnorm); th = ck.th; return;
    end
    th = ck.th; fprintf('   resuming %s at iteration %d (cost %.5f)\n', ckpt, ck.iteration, ck.resnorm);
end
J.thetaFix = th; J.free = free;
opt = optimoptions('lsqnonlin', 'Display', 'iter', 'UseParallel', true, ...
    'FiniteDifferenceStepSize', 1e-2, 'MaxIterations', maxIter, ...
    'FunctionTolerance', 1e-6, 'StepTolerance', 1e-5, ...
    'OutputFcn', @(x, ov, st) saveCkpt(x, ov, st, th, free, ckpt));
x = lsqnonlin(@(x) resJoint(x, J, SSrel, SSact), th(free), lb(free), ub(free), opt);
th(free) = x;
ck.th = th; ck.resnorm = sum(resJoint(x, J, SSrel, SSact).^2); ck.iteration = maxIter; ck.done = true; %#ok<STRNU>
save(ckpt, '-struct', 'ck');
end

function stop = saveCkpt(x, ov, state, th, free, f)
stop = false;
if strcmp(state, 'iter')
    th(free) = x;
    ck.th = th; ck.resnorm = ov.resnorm; ck.iteration = ov.iteration; ck.done = false; ck.time = datetime('now'); %#ok<STRNU>
    save(f, '-struct', 'ck');
end
end
