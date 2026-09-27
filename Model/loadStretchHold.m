function S = loadStretchHold(file, tWin, nBins)
% loadStretchHold  Load a 10 kHz ASI export and prepare the first stretch-hold
% (0.95 -> 1.175) for fitting. Utility (not an entry point) for FitFirstStretch.m.
%
%   S.t, S.F   raw time (s, 0 = ramp onset) and baseline-subtracted force (kPa)
%   S.ramp     .Vfun (model velocity, um/s, from the measured length trace),
%              .tEnd (end of the ramp section)
%   S.tSim     time grid the model must be evaluated on
%   S.A        sparse bin-averaging matrix: S.A*F(S.tSim) -> log-binned force
%   S.tb, S.Fb log-binned (geometric-mean) times and data force
%
% Log binning replaces the uniform sampling: equal weight per time decade,
% so the ~ms peak is not swamped by ~1e5 hold samples, and bin averaging of
% both data and model damps the ~1.5 kHz transducer ringing after the peak.

if nargin < 2 || isempty(tWin), tWin = 29.9; end
if nargin < 3 || isempty(nBins), nBins = 100; end
Lmax = 0.225;

% file may be a cell array of repeats of the same protocol: they are sampled
% phase-locked to the stimulus, so they are averaged sample by sample
% (baselines removed per trial below via the averaged pre-ramp window)
if ischar(file), file = {file}; end
for k = 1:numel(file)
    ek = readtable(file{k}, 'VariableNamingRule', 'preserve');
    ek.Properties.VariableNames = {'Time', 'L', 'F', 'SL'};
    if k == 1
        e = ek; Fsum = ek.F; Lsum = ek.L; n = height(ek);
    else
        n = min(n, height(ek));
        Fsum = Fsum(1:n) + ek.F(1:n); Lsum = Lsum(1:n) + ek.L(1:n);
    end
end
e = e(1:n, :); e.F = Fsum/numel(file); e.L = Lsum/numel(file);
t = e.Time/1000;

% ramp onset: first sample moving faster than 1 L0/s near t = 10 s
v = gradient(movmean(e.L, 5), t);
i0 = find(abs(v) > 1 & t > 9 & t < 11, 1);
if isempty(i0) % slow ramp (< 1 L0/s): first 0.5 % level crossing
    L0pre = mean(e.L(t > 1 & t < 9.9));
    Lplat = max(movmean(e.L(t > 9.9 & t < 25), 101));
    i0 = find(t > 9.9 & e.L > L0pre + 0.005*(Lplat - L0pre), 1);
end
t0 = t(i0);
pre  = t > 1 & t < t0 - 0.05;
Lpre = mean(e.L(pre));
Lplat = max(movmean(e.L(t > t0 & t < t0 + 15), 101));
t99 = t(find(t > t0 & e.L >= Lpre + 0.99*(Lplat - Lpre), 1));
Lpost = mean(e.L(t > t99 + 0.2 & t < t99 + 1));
S.F0 = mean(e.F(pre));
S.noise = std(e.F(pre));

keep = t >= t0 & t <= t0 + tWin;
S.t = t(keep) - t0;
S.F = e.F(keep) - S.F0;
S.Lraw = e.L(keep);

% model length = measured length mapped onto 0..Lmax
Lm = min(max((movmean(S.Lraw, 3) - Lpre)/(Lpost - Lpre)*Lmax, 0), Lmax);
Lm(1) = 0;
iEnd = find(Lm >= 0.998*Lmax, 1);
S.ramp.tEnd = S.t(iEnd);
pp = pchip(S.t(1:iEnd), Lm(1:iEnd));
[br, cf] = unmkpp(pp);
ppd = mkpp(br, cf(:, 1:3).*[3 2 1]);
S.ramp.Vfun = @(tq) ppval(ppd, min(max(tq, br(1)), br(end)));
S.ramp.Lfun = @(tq) ppval(pp, min(max(tq, br(1)), br(end)));

% log bins over (dt, tWin]; drop empty bins
dt = median(diff(S.t));
edges = logspace(log10(dt/2), log10(tWin), nBins + 1);
[~, ~, bin] = histcounts(S.t, edges);
ok = bin > 0;
[ub, ~, g] = unique(bin(ok));
nb = numel(ub);
% fine region: evaluate the model at every sample; beyond it at the bin mean
tFine = S.ramp.tEnd + 0.02;
tOk = S.t(ok);
fineBins = unique(g(tOk <= tFine));
fineMask = ismember(g, fineBins);
coarseBins = setdiff((1:nb)', fineBins);
tMean = accumarray(g, tOk, [nb 1], @mean);
tS = [tOk(fineMask); tMean(coarseBins)];
rowBin = [g(fineMask); coarseBins];
[S.tSim, ~, ic] = unique(tS);
cnt = accumarray(rowBin, 1, [nb 1]);
S.A = sparse(rowBin, ic, 1./cnt(rowBin), nb, numel(S.tSim));
S.Fb = accumarray(g, S.F(ok), [nb 1], @mean);
S.tb = accumarray(g, S.t(ok), [nb 1], @(x) exp(mean(log(x))));
S.nb = accumarray(g, 1, [nb 1]);
S.file = file;
end
