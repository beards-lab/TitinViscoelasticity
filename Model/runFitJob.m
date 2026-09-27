function out = runFitJob(spec, SS, th0, lb, ub, maxIter, ckptFile)
% runFitJob  lsqnonlin on resVariant, for use with batch(... 'Pool', n) so a
% long fit runs asynchronously. Utility (not an entry point) for FitFirstStretch.m.
% ckptFile (optional): .mat file updated every iteration with the current
% point (th, resnorm, iteration), so a running job can be inspected.
if nargin < 7, ckptFile = ''; end
opt = optimoptions('lsqnonlin', 'Display', 'iter', 'UseParallel', true, ...
    'FiniteDifferenceStepSize', 1e-2, 'MaxIterations', maxIter, ...
    'FunctionTolerance', 1e-6, 'StepTolerance', 1e-5);
if ~isempty(ckptFile)
    opt = optimoptions(opt, 'OutputFcn', @(x, ov, st) saveCkpt(x, ov, st, ckptFile));
end
t = tic;
[out.th, out.res, out.r, out.flag, out.output] = lsqnonlin(@(th) resVariant(th, spec, SS), th0, lb, ub, opt);
[~, out.p, out.sensor] = resVariant(out.th, spec, SS);
out.time = toc(t);
end

function stop = saveCkpt(x, ov, state, f)
stop = false;
if strcmp(state, 'iter')
    ck.th = x; ck.resnorm = ov.resnorm; ck.iteration = ov.iteration; ck.time = datetime('now'); %#ok<STRNU>
    save(f, '-struct', 'ck');
end
end
