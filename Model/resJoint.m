function r = resJoint(x, J, SSrel, SSact)
% resJoint  Residual vector of a joint low/high-Ca stretch-hold fit, for
% lsqnonlin. Utility (not an entry point) for OptimizeStretchHoldJoint.m.
%   x      free entries of theta (J.free, logical mask over theta); the
%          others are taken from J.thetaFix
%   SSrel  relaxed datasets, SSact active datasets (loadStretchHold /
%          loadRestretch); J.wRel, J.wAct their weights; J.odeOpts solver options
% Same per-dataset cost as costStretchHold.m: sum(r.^2) = sum_q w_q*mean(res_q.^2).

theta = J.thetaFix;
theta(J.free) = x;
[pRel, pAct, senRel, senAct] = jointParams(theta, J);
SS = [SSrel, SSact];
P = [repmat({pRel}, 1, numel(SSrel)), repmat({pAct}, 1, numel(SSact))];
Sen = [repmat({senRel}, 1, numel(SSrel)), repmat({senAct}, 1, numel(SSact))];
w = [J.wRel, J.wAct];
r = cell(numel(SS), 1);
for q = 1:numel(SS)
    nb = numel(SS{q}.Fb);
    if w(q) == 0
        r{q} = zeros(0, 1); continue;
    end
    Fb = modelBinned(P{q}, SS{q}, Sen{q}, J.odeOpts);
    if any(~isfinite(Fb))
        Fb = 10*max(SS{q}.Fb)*ones(nb, 1);  % keeps r finite for the solver
    end
    r{q} = sqrt(w(q)/nb)*(Fb - SS{q}.Fb)/max(SS{q}.Fb);
end
r = vertcat(r{:});
end
