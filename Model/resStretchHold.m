function [r, params] = resStretchHold(theta, pBase, SS, wts, peakW)
% resStretchHold  Residual vector for lsqnonlin on the first stretch-hold.
% Utility (not an entry point) for FitFirstStretch.m.
%
% theta is a decorrelated, log-scaled version of the pCa 11 parameters:
%   theta = [log mu, log Fp0(sRef), np, log kd, nd, log U0(sRef), nU, Fss, delU]
% optionally followed by [log f0, zeta] of the force transducer (sensorFilter.m)
% and n_ss (parallel element exponent)
% with Fp0(sRef) = kp*sRef^np (folded proximal force at sRef) and
% U0(sRef) = alphaU*sRef^nU (first unfolding rate at sRef). In the raw
% params (kp,np) and (alphaU,nU) are nearly collinear, which makes the
% simplex crawl along a narrow valley.
% Called with theta = [] and pBase = params it returns theta instead.

sRef = 0.2;
if isempty(theta)
    p = pBase;
    r = [log(p(9)), log(p(3)*sRef^p(4)), p(4), log(p(5)), p(6), ...
         log(p(7)*sRef^p(8)), p(8), p(1), p(10)];
    return;
end
if nargin < 4 || isempty(wts), wts = ones(1, numel(SS)); end
if nargin < 5 || isempty(peakW), peakW = 0; end
params = pBase;
params(9)  = exp(theta(1));
params(4)  = theta(3);
params(3)  = exp(theta(2))/sRef^theta(3);
params(5)  = exp(theta(4));
params(6)  = theta(5);
params(8)  = theta(7);
params(7)  = exp(theta(6))/sRef^theta(7);
params(1)  = theta(8);
params(10) = theta(9);
sensor = [];
if numel(theta) >= 11
    sensor = [exp(theta(10)) theta(11)];
end
if numel(theta) >= 12
    params(2) = theta(12);
end

r = cell(numel(SS), 1);
for q = 1:numel(SS)
    nb = numel(SS{q}.Fb);
    Fb = modelBinned(params, SS{q}, sensor);
    if any(~isfinite(Fb))
        Fb = 10*max(SS{q}.Fb)*ones(nb, 1); % smooth-ish penalty, keeps r finite
    end
    r{q} = sqrt(wts(q)/nb)*(Fb - SS{q}.Fb)/max(SS{q}.Fb);
    if peakW > 0 % optional extra weight on the (binned) peak value
        r{q}(end+1) = sqrt(peakW*wts(q))*(max(Fb) - max(SS{q}.Fb))/max(SS{q}.Fb);
    end
end
r = vertcat(r{:});
end
