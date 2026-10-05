function S = loadProtocol(file, tStop, nBins)
% loadProtocol  Prepare the WHOLE 2025-11-21 refolding protocol of one 10 kHz
% *_refolding_*.txt export for simStretchHold/modelBinned: first stretch
% (10 s), 30 s hold, release to 0.95 L0 (40 s), gap, restretch, hold, final
% release to 0.80 L0 (80 s, clipped to model L = 0). Utility (not an entry
% point) for SweepRefolding.m. Same conventions as loadRestretch.m:
%   sim time 0 = onset of the first stretch; the model is driven by the
%   measured length during every movement (S.ramp.seg), holds are pinned.
%   S.tb, S.Fb, S.A, S.tSim: every sample during a movement (+20 ms), then
%   nBins log bins in the time since that movement's onset, up to the next.
%   S.ev: movement onsets (sim time), S.evEnd: their ends,
%   S.evBin: per bin, index of the movement it follows.
if nargin < 2 || isempty(tStop), tStop = 72; end   % sim time (80 s = 70)
if nargin < 3 || isempty(nBins), nBins = 60; end
Lmax = 0.225;

e = readtable(file, 'VariableNamingRule', 'preserve');
e.Properties.VariableNames = {'Time', 'L', 'F', 'SL'};
t = e.Time/1000;
v = gradient(movmean(e.L, 5), t);
i0 = find(abs(v) > 1 & t > 9 & t < 11, 1);
t0 = t(i0);
pre = t > 1 & t < t0 - 0.05;
Lpre = mean(e.L(pre));
Lpost = mean(e.L(t > t0 + 0.2 & t < t0 + 1));
S.F0 = mean(e.F(pre)); S.noise = std(e.F(pre));
tt = t - t0;
dt = median(diff(tt));
Lm = min(max((movmean(e.L, 3) - Lpre)/(Lpost - Lpre)*Lmax, 0), Lmax);

% movement episodes: |v| > 0.5 L0/s, pieces < 2 ms apart merged
mv = abs(v) > 0.5 & tt >= -dt & tt <= tStop;
d = diff([0; mv; 0]); on = find(d == 1); off = find(d == -1) - 1;
k = 1;
while k < numel(on)
    if tt(on(k+1)) - tt(off(k)) < 2e-3
        off(k) = off(k+1); on(k+1) = []; off(k+1) = [];
    else
        k = k + 1;
    end
end
seg = cell(0, 3); tPrev = 0; ev = zeros(numel(on), 1); evEnd = ev;
for k = 1:numel(on)
    a = on(k) - 2; b = off(k) + 5;
    if k == 1, a = find(tt >= 0, 1); end
    bEnd = b + round(0.05/dt);                 % level after the movement,
    if k < numel(on), bEnd = min(bEnd, on(k+1) - 3); end  % before the next one
    plat = median(Lm(b:max(b, bEnd)));
    if plat > 0.98*Lmax, plat = Lmax; elseif plat < 0.02*Lmax, plat = 0; end
    Lk = Lm(a:b); Lk(end) = plat;
    if k == 1, Lk(1) = 0; end
    Vk = pchipVel(tt(a:b), Lk);
    if tt(a) > tPrev
        seg(end+1, :) = {[tPrev tt(a)], 0, NaN}; %#ok<AGROW>
    end
    seg(end+1, :) = {[tt(a) tt(b)], Vk, plat}; %#ok<AGROW>
    tPrev = tt(b); ev(k) = tt(a); evEnd(k) = tt(b);
end
seg(end+1, :) = {[tPrev tStop], 0, NaN};
S.ramp.seg = seg;
S.ramp.tEnd = evEnd(1);
S.ev = ev; S.evEnd = evEnd;
win = [max(ev - 1e-3, 0), evEnd + 15e-3];   % sensor windows, overlaps merged
k = 1;
while k < size(win, 1)
    if win(k+1, 1) <= win(k, 2)
        win(k, 2) = max(win(k, 2), win(k+1, 2)); win(k+1, :) = [];
    else
        k = k + 1;
    end
end
S.filtWin = win;

% bins: all samples during each movement + 20 ms, log bins after
keep = tt >= 0 & tt <= tStop;
S.t = tt(keep); S.F = e.F(keep) - S.F0; S.Lraw = e.L(keep);
bin = zeros(size(S.t)); fine = false(size(S.t)); nb = 0;
for k = 1:numel(ev)
    tNext = tStop; if k < numel(ev), tNext = ev(k+1); end
    in = S.t >= ev(k) & S.t < tNext;
    if k == numel(ev), in = in | S.t == tStop; end
    tF = min(evEnd(k) + 0.02, tNext);          % short gaps: all fine
    f = in & S.t <= tF;                        % one bin per sample
    idx = find(f);
    bin(idx) = nb + (1:numel(idx))'; fine(idx) = true;
    nb = nb + numel(idx);
    r = in & ~f;
    if ~any(r), continue; end
    tau = S.t(r) - ev(k);
    edges = logspace(log10(tF - ev(k)), log10(tNext - ev(k) + dt), nBins + 1);
    [~, ~, bk] = histcounts(tau, edges);
    [ub, ~, g] = unique(bk(bk > 0));
    ri = find(r); ri = ri(bk > 0);
    bin(ri) = nb + g;
    nb = nb + numel(ub);
end
use = bin > 0; gU = bin(use); tU = S.t(use); fU = fine(use);
coarse = setdiff((1:nb)', unique(gU(fU)));
tMean = accumarray(gU, tU, [nb 1], @mean);
tS = [tU(fU); tMean(coarse)];
rowBin = [gU(fU); coarse];
[S.tSim, ~, ic] = unique(tS);
cnt = accumarray(rowBin, 1, [nb 1]);
S.A = sparse(rowBin, ic, 1./cnt(rowBin), nb, numel(S.tSim));
S.Fb = accumarray(gU, S.F(use), [nb 1], @mean);
S.tb = tMean;
S.nb = accumarray(gU, 1, [nb 1]);
% movement each bin belongs to (bins are numbered in time order)
evBin = zeros(nb, 1);
for k = 1:numel(ev), evBin(S.tb >= ev(k)) = k; end
S.evBin = evBin;
S.file = file;
end

function Vfun = pchipVel(tk, Lk)
pp = pchip(tk, Lk);
[br, cf] = unmkpp(pp);
ppd = mkpp(br, cf(:, 1:3).*[3 2 1]);
Vfun = @(tq) ppval(ppd, min(max(tq, br(1)), br(end)));
end
