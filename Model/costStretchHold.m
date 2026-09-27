function [cost, parts, Fm] = costStretchHold(params, SS, wts, sensor)
% costStretchHold  Log-binned cost of simStretchHold against loadStretchHold
% datasets. Utility (not an entry point) for FitFirstStretch.m.
%   cost  = sum_q wts(q) * mean(((A*Fmodel - Fb)/max(Fb)).^2)
%   parts = per-dataset [cost, model bin-peak, data bin-peak]

if nargin < 3 || isempty(wts), wts = ones(1, numel(SS)); end
if nargin < 4, sensor = []; end % [f0 zeta] transducer, see sensorFilter.m
parts = nan(numel(SS), 3);
Fm = cell(1, numel(SS));
cost = 0;
for q = 1:numel(SS)
    if wts(q) == 0, continue; end
    Fb = modelBinned(params, SS{q}, sensor);
    if any(~isfinite(Fb))
        cost = 1e3; return;
    end
    Fm{q} = Fb;
    c = mean(((Fb - SS{q}.Fb)/max(SS{q}.Fb)).^2);
    parts(q, :) = [c, max(Fb), max(SS{q}.Fb)];
    cost = cost + wts(q)*c;
end
end
