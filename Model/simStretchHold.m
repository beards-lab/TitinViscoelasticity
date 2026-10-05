function [F, L, fAtt] = simStretchHold(params, tOut, ramp, opts, Nx)
% simStretchHold  Lean single stretch-hold; relaxed (kA, kD NaN) or active
% (PEVK attachment states, kA/kD pre-equilibrated at rest).
% Utility (not an entry point) used by FitFirstStretch.m. Same model as
% RunCombinedModel.m / dXdT.m (Nx = 15, Ng = 10, ode15s), but
%   - driven by the measured length trace (ramp.Vfun, ramp.tEnd) instead of
%     an ideal constant-velocity ramp - the peak of a ~2.5 ms stretch is very
%     sensitive to the true motor trajectory,
%   - returns force on the fixed output grid tOut (s, 0 = ramp onset), so the
%     cost is a smooth function of params (no dependence on ode15s steps).
%
% params: [Fss n_ss kp np kd nd alphaU nU mu delU kA kD alphaF_0 ...
%          mu1 m_mu Fbeta]
%         (kA, kD NaN = relaxed; alphaF_0 optional, default 0)
% Optional structural variants (NaN or absent = original model):
%   mu1, m_mu  strain-dependent viscosity mu(s) = mu + mu1*(s/Lmax)^m_mu
%   signedFd   (index 17) 1 = distal force output signed like in dXdT
%   clipFd     (index 18) 1 = distal element tension-only in the dynamics too
%   eta        (index 19) parallel fibre-level viscosity, F += eta*dL/dt
%   wSlack     (index 20) mu1 term uses (s - w*slack(n)) instead of s, w in [0,1]
%   cComp      (index 21) distal compressive stiffness factor in [0,1], used in
%              dynamics and force output (1 = dXdT dynamics)
%   muRec      (index 22) drag of recoiling chains (Vp < 0); NaN = same mu
%   kM, tauM   (index 23:24) parallel linear Maxwell element (spring kM, time tauM)
%   kDf        (index 25) slip-bond detachment kD*(1 + kDf*Fd) (active only)
%   F_R        (index 26) smooth refolding alphaF_0*exp(-Fp(s,n+1)/F_R) instead
%              of the on/off mask of dXdT (continuous in all parameters)
%   kDslack    (index 27) extra detachment rate of attached chains with a slack
%              distal segment (L < s), active only
%   gammaF     (index 28) n-dependent refolding alphaF_0*((n+1)/Ng)^gammaF for
%              n+1 -> n (combines with the on/off mask, or with F_R)
%   Fbeta      Bell-type force-dependent unfolding
%              U(n->n+1) = alphaU*exp(Fp(s,n)/Fbeta)*(Ng-n), nU unused

if nargin < 4 || isempty(opts)
    opts = odeset('RelTol', 1e-4, 'AbsTol', 1e-4); % as in RunCombinedModel
end
if nargin < 5 || isempty(Nx)
    Nx = 15; % as in RunCombinedModel; larger Nx = finer strain grid
end
F = nan(size(tOut)); L = nan(size(tOut)); fAtt = nan(size(tOut));
if any(params(~isnan(params)) < 0)
    return;
end

Lmax = 0.225; Ng = 10; L_0 = 1;
ds = Lmax/(Nx-1); s = (0:Nx-1)'.*ds;

Fss = params(1); n_ss = params(2); kp = params(3); np = params(4);
kd = params(5); nd = params(6); alphaU = params(7); nU = params(8);
mu = params(9); delU = params(10)/Ng;
alphaF_0 = 0;
if numel(params) >= 13 && ~isnan(params(13))
    alphaF_0 = params(13);
end

slack = (0:Ng).*delU;
Fp = kp*(max(0, s-slack)/L_0).^np;
RU = alphaU*((max(0, s-slack(1:Ng))/L_0).^nU).*(ones(Nx,1).*(Ng - (0:Ng-1)));
if numel(params) >= 16 && ~isnan(params(16))
    RU = alphaU*exp(Fp(:, 1:Ng)/params(16)).*(ones(Nx,1).*(Ng - (0:Ng-1)));
end
if numel(params) >= 15 && all(~isnan(params(14:15)))
    if numel(params) >= 20 && ~isnan(params(20)) && params(20) > 0
        % drag set by the proximal extension beyond a fraction w of the
        % unfolded slack (w = 1: s - slack(n); w = 0: s)
        mu = mu + params(14)*(max(0, s - params(20)*slack)/Lmax).^params(15); % Nx x (Ng+1)
    else
        mu = mu + params(14)*(s/Lmax).^params(15); % Nx x 1 column
    end
end
cComp = 1;                                       % distal compression factor
if numel(params) >= 18 && params(18) == 1, cComp = 0; end % tension-only
if numel(params) >= 21 && ~isnan(params(21)), cComp = params(21); end
if numel(params) >= 26 && ~isnan(params(26)) && alphaF_0 > 0
    % smooth refolding n+1 -> n suppressed by the proximal force of state n+1
    alphaF_0 = alphaF_0*exp(-Fp(:, 2:Ng+1)/params(26));   % Nx x Ng
end
if numel(params) >= 28 && ~isnan(params(28)) && any(alphaF_0(:) > 0)
    % n-dependent refolding: n+1 -> n at alphaF_0*((n+1)/Ng)^gammaF (fastest
    % for chains with many unfolded domains); on/off mask kept unless F_R set
    fac = ((1:Ng)/Ng).^params(28);
    if isscalar(alphaF_0)
        alphaF_0 = alphaF_0*fac.*ones(Nx, 1);
        alphaF_0(Fp(:, 2:Ng+1) > 0) = 0;          % same mask as dXdT
    else
        alphaF_0 = alphaF_0.*fac;
    end
end
kDf = 0;                                         % force-dependent detachment
if numel(params) >= 25 && ~isnan(params(25)), kDf = params(25); end
muRec = NaN;                                     % recoil drag (Vp < 0)
if numel(params) >= 22, muRec = params(22); end
kDslack = 0;                                     % detachment of slack attached chains
if numel(params) >= 27 && ~isnan(params(27)), kDslack = params(27); end
isAct = numel(params) >= 12 && all(~isnan(params(11:12)));
if ~isscalar(mu) || cComp ~= 1 || ~isnan(muRec) || kDf > 0 || ~isscalar(alphaF_0) || kDslack > 0 ...
        || (isAct && any(alphaF_0(:) > 0))   % dXdT refolds attached chains with the pu flux
    odefun = @(t, x, varargin) dXdTvar(t, x, varargin{:}, cComp, muRec, kDslack);
else
    odefun = @dXdT;
end

% active (pCa < 11): PEVK attachment states pa, pre-equilibrated at rest
active = numel(params) >= 12 && all(~isnan(params(11:12)));
kA = NaN; kD = NaN;
pu = zeros(Nx, Ng+1); pu(1,1) = 1/ds;
if active
    kA = params(11); kD = params(12);
    fA = kA/(kA + kD);
    pa = pu*fA; pu = pu*(1 - fA);
    x0 = [pu(:); pa(:); 0];
else
    x0 = [pu(:); 0];
end

% sparsity of the ODE Jacobian (speed only; same solution and tolerances):
% p(s,n) couples to s+-1 (sliding), n+-1 (unfolding/refolding), its
% attached/unattached twin, and the length L (last state)
if isempty(odeget(opts, 'JPattern'))
    nh = Nx*(Ng+1); nb = numel(x0) - 1;
    [ii, jj] = ndgrid(1:Nx, 1:Ng+1); k = ii + Nx*(jj - 1);
    nbr = {k, k(max(ii-1, 1) + Nx*(jj-1)), k(min(ii+1, Nx) + Nx*(jj-1)), ...
           k(ii + Nx*(max(jj-1, 1) - 1)), k(ii + Nx*(min(jj+1, Ng+1) - 1))};
    rI = []; cI = [];
    for b = 0:(nb/nh - 1)
        for q = 1:numel(nbr), rI = [rI; k(:) + b*nh]; cI = [cI; nbr{q}(:) + b*nh]; end %#ok<AGROW>
    end
    if nb > nh                                   % attach/detach twins
        rI = [rI; k(:); k(:) + nh]; cI = [cI; k(:) + nh; k(:)];
    end
    rI = [rI; (1:nb)']; cI = [cI; (nb+1)*ones(nb, 1)];
    opts = odeset(opts, 'JPattern', sparse(rI, cI, 1, nb+1, nb+1) ~= 0);
end
tOut = tOut(:);
% sections: {[t1 t2], velocity, pinned length at t2 (NaN = none)}
if isfield(ramp, 'seg')      % general protocol (e.g. stretch-hold-release-restretch)
    sections = ramp.seg;
else                         % single stretch-hold
    sections = {[0 ramp.tEnd], ramp.Vfun, Lmax; [ramp.tEnd tOut(end)], 0, NaN};
end
X = nan(numel(tOut), numel(x0));
Vout = zeros(numel(tOut), 1);            % imposed velocity, for the eta term
lastwarn('');
for k = 1:size(sections, 1)
    tr = sections{k, 1};
    in = tOut > tr(1) & tOut <= tr(2);
    if k == 1
        in = in | tOut == 0;
    end
    tspan = unique([tr(1); tOut(in); tr(2)]);
    [tt, xx] = ode15s(odefun, tspan, x0, opts, Nx, Ng, ds, kA, kD, kd, Fp, RU, ...
        alphaF_0, mu, L_0, nd, kDf, sections{k, 2});
    if tt(end) < tr(2)
        return; % integration failure
    end
    X(in, :) = interp1(tt, xx, tOut(in));
    if isa(sections{k, 2}, 'function_handle')
        Vout(in) = sections{k, 2}(tOut(in));
    end
    x0 = xx(end, :)';
    if ~isnan(sections{k, 3})
        x0(end) = sections{k, 3}; % pin the hold length exactly
    end
end

L = X(:, end);
P = X(:, 1:end-1);                           % rows: time, cols: (s, n[, pa])
if active                                    % attached fraction of all chains
    nh = Nx*(Ng+1);
    fAtt = sum(P(:, nh+1:end), 2)./sum(P, 2);
end
if numel(params) >= 21 && ~isnan(params(21))
    % distal force output consistent with the dynamics' compressive branch
    Fd = kd*(max(0, (L - s')/L_0).^nd - params(21)*max(0, (s' - L)/L_0).^nd);
elseif numel(params) >= 17 && params(17) == 1
    % signed distal force, consistent with dXdT (compression allowed)
    Fd = kd*sign(L - s').*abs((L - s')/L_0).^nd;
else
    Fd = kd*max(0, (L - s')/L_0).^nd;         % time x Nx (as RunCombinedModel)
end
F = ds*sum(Fd.*reshape(sum(reshape(P, size(P, 1), Nx, []), 3), [], Nx), 2) ...
    + Fss/Lmax^n_ss*max(L, 0).^n_ss;
if numel(params) >= 19 && ~isnan(params(19))
    F = F + params(19)*Vout;              % parallel Newtonian viscosity
end
if numel(params) >= 24 && all(~isnan(params(23:24)))
    % parallel linear Maxwell element: dFM/dt = kM*V - FM/tauM (V imposed)
    kM = params(23); tauM = params(24); FM = 0; FMout = zeros(numel(tOut), 1);
    for k = 1:size(sections, 1)
        tr = sections{k, 1}; in = tOut >= tr(1) & tOut <= tr(2);
        Vk = sections{k, 2};
        if isa(Vk, 'function_handle')
            ts = unique([tr(1); tOut(in); tr(2)]);
            [tm, fm] = ode45(@(t, y) kM*Vk(t) - y/tauM, ts, FM, ...
                odeset('RelTol', 1e-6, 'AbsTol', 1e-8, 'MaxStep', 1e-4));
            if numel(ts) == 2, fm = interp1(tm, fm, ts); tm = ts; end
            FMout(in) = interp1(tm, fm, tOut(in));
            FM = fm(end);
        else                                % hold: exponential decay
            FMout(in) = FM*exp(-(tOut(in) - tr(1))/tauM);
            FM = FM*exp(-(tr(2) - tr(1))/tauM);
        end
    end
    F = F + FMout;
end
end
