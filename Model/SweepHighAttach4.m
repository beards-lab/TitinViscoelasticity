%% SweepHighAttach4.m
% ENTRY POINT. Four-parameter refinement of the high-attachment active model
% on the ramps only (2.8 ms pooled first stretch, 5.7 ms, 100 ms, 1 s; the
% restretch is left out because it needs refolding). Starts from stage 1 of
% OptimizeStretchHoldJoint.m seed 3 (jointFit_localtest/ckpt_seed3_stage1.mat);
% shared mechanics (soft distal relaxed set) and the sensor stay fixed.
% Free: mu1X, fA, r, kDf (chosen by a one-step sensitivity pass).
% Run as: batch('SweepHighAttach4', 'Pool', 8); checkpoints to
% jointFit_localtest/ckpt_sweep4.mat, result in jointFit_localtest/sweep4_highAttach.mat.

dd = '..\Data\2025 11 21 Export\';
rf = {'0ms', '5ms', '10ms', '50ms', '100ms', '1s', '10s', '30s'};
pad = @(p) [p, NaN(1, 27 - numel(p))];
Rs = load('fitStretchHold_Relax_softDistal.mat'); Rb = load('fitStretchHold_Relax_best.mat'); Ab = load('fitStretchHold_Active_best.mat');
pRelSoft = pad(Rs.p); pRelSoft([11 12 25 27]) = NaN;
pAct0 = pRelSoft; fA = 0.95; r = 1000; pAct0([11 12 25 27]) = [fA*r, (1-fA)*r, 0.3, 3000];
J = struct('sharedIdx', [9 3 4 5 6 7 8 1 2 10 14 15 19 23 24], ...
    'sharedLog', logical([1 1 0 1 0 1 0 0 0 0 1 0 1 1 1]), 'refKp', true, 'refAU', true, ...
    'caNames', {{'kpX', 'nU', 'aUX', 'mu1X', 'fA', 'r', 'kDf', 'kDslack'}}, ...
    'senRel', Rb.sensor, 'senAct', Ab.sensor, 'pRel0', pRelSoft, 'pAct0', pAct0, ...
    'wRel', [], 'wAct', [1 1 1 1], 'odeOpts', odeset('RelTol', 1e-6, 'AbsTol', 1e-6));
SSact = {loadStretchHold(strcat(dd, rf, '_refolding_Active.txt'), [], 200), ...
         loadStretchHold([dd '5ms_Ramp_Active.txt'], [], 200), ...
         loadStretchHold([dd '0.1s_Ramp_Active.txt'], [], 200), ...
         loadStretchHold([dd '1s_Ramp_Active.txt'], [], 200)};
SSrel = {};
ck = load('jointFit_localtest\ckpt_seed3_stage1.mat');
th0 = ck.th; caIdx = 15 + (1:8);
free = false(size(th0)); free(caIdx([4 5 6 7])) = true;         % mu1X fA r kDf
lb = [log(1e-2) -12 log(1e-3) log(1e-4)]; ub = [log(1e3) 12 log(1e6) log(100)];
J.thetaFix = th0; J.free = free;
x0 = min(max(th0(free), lb + 1e-6), ub - 1e-6);
c0 = sum(resJoint(x0, J, SSrel, SSact).^2);
fprintf('start cost (4 ramps) %.5f\n', c0);
ckf = 'jointFit_localtest\ckpt_sweep4.mat';
opt = optimoptions('lsqnonlin', 'Display', 'iter', 'UseParallel', true, 'FiniteDifferenceStepSize', 1e-2, ...
    'MaxIterations', 20, 'FunctionTolerance', 1e-6, 'StepTolerance', 1e-5, ...
    'OutputFcn', @(x, ov, st) ckptSave(x, ov, st, ckf));
[x, rn, ~, flag, out] = lsqnonlin(@(x) resJoint(x, J, SSrel, SSact), x0, lb, ub, opt);
th4 = th0; th4(free) = x;
[~, pA4, ~, sA4] = jointParams(th4, J);
save('jointFit_localtest\sweep4_highAttach.mat', 'th4', 'pA4', 'sA4', 'rn', 'c0', 'flag', 'out', 'J');
fprintf('cost %.5f -> %.5f, flag %d, %d iterations\n', c0, rn, flag, out.iterations);
fprintf('fA %.5f  r %.4g  kDf %.4g  mu1 %.4g\n', pA4(11)/(pA4(11)+pA4(12)), pA4(11)+pA4(12), pA4(25), pA4(14));

function stop = ckptSave(x, ov, st, f)
stop = false;
if strcmp(st, 'iter'), it = ov.iteration; res = ov.resnorm; save(f, 'x', 'it', 'res'); end
end
