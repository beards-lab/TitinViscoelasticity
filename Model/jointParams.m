function [pRel, pAct, senRel, senAct] = jointParams(theta, J)
% jointParams  Map one theta vector to the relaxed and active parameter sets
% of a joint low/high-Ca fit. Utility (not an entry point) for
% OptimizeStretchHoldJoint.m; params layout as in simStretchHold.m.
%
% theta = [shared | Ca-specific | log f0_rel, zeta_rel, log f0_act, zeta_act]
%   shared      J.sharedIdx (J.sharedLog: log-scaled). With J.refKp / J.refAU
%               the values at 3 / 7 are Fp0(sRef) = kp*sRef^np and
%               U0(sRef) = alphaU*sRef^nU (decorrelated from np / nU).
%   Ca-specific J.caNames, applied to a copy of the relaxed set:
%     'kpX'     log multiplier of kp          (published: kp is Ca-dependent)
%     'nU'      nU at high Ca                 (published: nU is Ca-dependent)
%     'aUX'     log multiplier of U0(sRef)    (published: alphaU is Ca-dependent)
%     'mu1X'    log multiplier of the strain-dependent drag mu1
%     'fA'      logit of the resting attached fraction kA/(kA+kD)
%     'r'       log of the attach/detach cycling rate kA + kD (1/s)
%     'kDf'     log slip-bond force sensitivity (1/kPa), detach kD*(1 + kDf*Fd)
%     'kDslack' log extra detachment rate of attached chains with a slack
%               distal segment (1/s)
% Called with theta = [] it returns theta for J.pRel0 / J.pAct0 / sensors.

sRef = 0.2;
if isempty(theta)
    pR = J.pRel0; pA = J.pAct0;
    v = pR(J.sharedIdx);
    if J.refKp, v(J.sharedIdx == 3) = pR(3)*sRef^pR(4); end
    if J.refAU, v(J.sharedIdx == 7) = pR(7)*sRef^pR(8); end
    v(J.sharedLog) = log(v(J.sharedLog));
    c = zeros(1, numel(J.caNames));
    for i = 1:numel(J.caNames)
        switch J.caNames{i}
            case 'kpX',     c(i) = log(pA(3)/pR(3));
            case 'nU',      c(i) = pA(8);
            case 'aUX',     c(i) = log((pA(7)*sRef^pA(8))/(pR(7)*sRef^pR(8)));
            case 'mu1X',    c(i) = log(pA(14)/pR(14));
            case 'fA',      f = pA(11)/(pA(11) + pA(12)); c(i) = log(f/(1 - f));
            case 'r',       c(i) = log(pA(11) + pA(12));
            case 'kDf',     c(i) = log(pA(25));
            case 'kDslack', c(i) = log(pA(27));
        end
    end
    pRel = [v, c, log(J.senRel(1)), J.senRel(2), log(J.senAct(1)), J.senAct(2)];
    return;
end

n = numel(J.sharedIdx);
pR = J.pRel0;
v = theta(1:n);
v(J.sharedLog) = exp(v(J.sharedLog));
pR(J.sharedIdx) = v;
if J.refKp, pR(3) = pR(3)/sRef^pR(4); end
if J.refAU, pR(7) = pR(7)/sRef^pR(8); end
pR([11 12 25 27]) = NaN;                 % no PEVK attachment at pCa 11

pA = pR;
pA([11 12 25 27]) = J.pAct0([11 12 25 27]);
c = theta(n+1:n+numel(J.caNames));
U0 = pR(7)*sRef^pR(8);                   % relaxed unfolding rate at sRef
fA = NaN; r = NaN;
for i = 1:numel(J.caNames)               % nU before aUX so U0(sRef) is kept
    switch J.caNames{i}
        case 'kpX',     pA(3) = pR(3)*exp(c(i));
        case 'nU',      pA(8) = c(i);
        case 'mu1X',    pA(14) = pR(14)*exp(c(i));
        case 'fA',      fA = 1/(1 + exp(-c(i)));
        case 'r',       r = exp(c(i));
        case 'kDf',     pA(25) = exp(c(i));
        case 'kDslack', pA(27) = exp(c(i));
    end
end
iaU = find(strcmp(J.caNames, 'aUX'), 1);
if isempty(iaU)
    pA(7) = U0/sRef^pA(8);
else
    pA(7) = U0*exp(c(iaU))/sRef^pA(8);
end
if isnan(fA), fA = J.pAct0(11)/(J.pAct0(11) + J.pAct0(12)); end
if isnan(r),  r = J.pAct0(11) + J.pAct0(12); end
pA(11) = fA*r; pA(12) = (1 - fA)*r;

k = n + numel(J.caNames);
senRel = [exp(theta(k+1)) theta(k+2)];
senAct = [exp(theta(k+3)) theta(k+4)];
pRel = pR; pAct = pA;
end
