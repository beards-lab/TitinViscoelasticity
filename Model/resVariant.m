function [r, params, sensor] = resVariant(theta, spec, SS)
% resVariant  Generic residual for lsqnonlin over any subset of params.
% Utility (not an entry point) for FitFirstStretch.m.
%   spec.pBase   full parameter vector (see simStretchHold.m, up to 16 entries)
%   spec.idx     indices of params that theta(1:numel(idx)) sets
%   spec.islog   logical, same size as idx: theta is log(param)
%   spec.refKp   if true, the value at index 3 is Fp0(sRef) = kp*sRef^np
%   spec.refAU   if true, the value at index 7 is U0(sRef) = alphaU*sRef^nU
%   spec.wts     dataset weights (default ones)
%   spec.odeOpts optional odeset for the fit (default: simStretchHold's 1e-4)
%   theta(end-1:end) = [log f0, zeta] of the force transducer (sensorFilter.m)
% Called with theta = [] it returns the theta for spec.pBase and spec.sensor.

sRef = 0.2;
p = spec.pBase;
if isempty(theta)
    v = p(spec.idx);
    if spec.refKp, v(spec.idx == 3) = p(3)*sRef^p(4); end
    if spec.refAU, v(spec.idx == 7) = p(7)*sRef^p(8); end
    v(spec.islog) = log(v(spec.islog));
    r = [v, log(spec.sensor(1)), spec.sensor(2)];
    return;
end
n = numel(spec.idx);
v = theta(1:n);
v(spec.islog) = exp(v(spec.islog));
p(spec.idx) = v;
if spec.refKp, p(3) = p(3)/sRef^p(4); end
if spec.refAU, p(7) = p(7)/sRef^p(8); end
params = p;
sensor = [exp(theta(n+1)) theta(n+2)];

wts = ones(1, numel(SS));
if isfield(spec, 'wts') && ~isempty(spec.wts), wts = spec.wts; end
odeOpts = [];
if isfield(spec, 'odeOpts'), odeOpts = spec.odeOpts; end
r = cell(numel(SS), 1);
for q = 1:numel(SS)
    nb = numel(SS{q}.Fb);
    Fb = modelBinned(params, SS{q}, sensor, odeOpts);
    if any(~isfinite(Fb))
        Fb = 10*max(SS{q}.Fb)*ones(nb, 1);
    end
    r{q} = sqrt(wts(q)/nb)*(Fb - SS{q}.Fb)/max(SS{q}.Fb);
end
r = vertcat(r{:});
end
