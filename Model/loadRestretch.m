function S = loadRestretch(file, tWin, nBins)
% loadRestretch  Prepare the release + restretch at ~40 s of a 10 kHz
% *_refolding_*.txt export (stretch at 10 s, 30 s hold, release, gap,
% restretch, hold) for simStretchHold/modelBinned. Utility (not an entry
% point) for FitFirstStretch.m. Same conventions as loadStretchHold.m:
%   sim time 0 = onset of the FIRST stretch; the model is driven by the
%   measured length through the whole protocol (S.ramp.seg);
%   S.tb, S.Fb, S.A, S.tSim are log bins in time since RESTRETCH onset
%   (S.tRs, sim time), plus one bin per sample during the release.
if nargin < 2 || isempty(tWin), tWin = 29.9; end
if nargin < 3 || isempty(nBins), nBins = 200; end
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
Lm = min(max((movmean(e.L, 3) - Lpre)/(Lpost - Lpre)*Lmax, 0), Lmax);

% first stretch
w1 = tt >= 0 & tt < 0.1;
i1 = find(w1 & Lm >= 0.998*Lmax, 1);
tE1 = tt(i1);
k1 = find(tt >= 0, 1):i1;
Lk = Lm(k1); Lk(1) = 0;
[V1, ~] = pchipVel(tt(k1), Lk);

% release + gap + restretch (release ~30 s after the first stretch)
mv = find(abs(v) > 1 & tt > 29 & tt < 31 + 40);
iR0 = mv(1) - 2;                                  % release start
wr = iR0:(iR0 + round(2/median(diff(tt))));        % search the next 2 s
[~, iMin] = min(Lm(wr)); iMin = wr(1) + iMin - 1;  % bottom of the release
iEnd = iMin - 1 + find(Lm(iMin:end) >= 0.998*Lmax, 1);
kr = iR0:iEnd;
Lk = Lm(kr); Lk(end) = Lmax;
[V2, L2] = pchipVel(tt(kr), Lk);
% restretch onset: last sample before the rise out of the release minimum
iRs = iMin - 1 + find(Lm(iMin:iEnd) > Lm(iMin) + 0.002*Lmax, 1) - 1;
S.tRs = tt(iRs); S.tRel = tt(iR0); S.tEndRs = tt(iEnd);
S.gap = tt(iRs) - tt(find(Lm(iR0:iMin) <= Lm(iMin) + 0.002*Lmax, 1) + iR0 - 1);
S.ramp.seg = {[0 tE1], V1, Lmax; [tE1 tt(iR0)], 0, NaN; ...
              [tt(iR0) tt(iEnd)], V2, Lmax; [tt(iEnd) S.tRs + tWin], 0, NaN};
S.ramp.tEnd = tE1;
S.ramp.Lfun = L2;
S.filtWin = [tt(iR0) - 1e-3, tt(iEnd) + 15e-3];

% data window around the restretch, relative time tau
keep = tt >= S.tRel - 1e-3 & tt <= S.tRs + tWin;
S.t = tt(keep) - S.tRs;               % relative to restretch onset
S.F = e.F(keep) - S.F0;
S.Lraw = e.L(keep);
dt = median(diff(S.t));
edges = logspace(log10(dt/2), log10(tWin), nBins + 1);
[~, ~, bin] = histcounts(S.t, edges);
nPre = nnz(S.t < dt/2);                % release samples: one bin each
bin(S.t < dt/2) = 0;
ok = bin > 0;
[~, ~, g] = unique(bin(ok));
g = g + nPre;
gAll = zeros(size(S.t)); gAll(S.t < dt/2) = 1:nPre; gAll(ok) = g;
use = gAll > 0; gU = gAll(use); tU = S.t(use);
nb = max(gU);
fineBins = unique(gU(tU <= S.tEndRs - S.tRs + 0.02));
fineMask = ismember(gU, fineBins);
coarse = setdiff((1:nb)', fineBins);
tMean = accumarray(gU, tU, [nb 1], @mean);
tS = [tU(fineMask); tMean(coarse)];
rowBin = [gU(fineMask); coarse];
[tSimRel, ~, ic] = unique(tS);
cnt = accumarray(rowBin, 1, [nb 1]);
S.A = sparse(rowBin, ic, 1./cnt(rowBin), nb, numel(tSimRel));
S.tSim = tSimRel + S.tRs;             % sim time
S.Fb = accumarray(gU, S.F(use), [nb 1], @mean);
S.tb = accumarray(gU, tU, [nb 1], @mean); % relative, can be negative
S.nb = accumarray(gU, 1, [nb 1]);
S.file = file;
end

function [Vfun, Lfun] = pchipVel(tk, Lk)
pp = pchip(tk, Lk);
[br, cf] = unmkpp(pp);
ppd = mkpp(br, cf(:, 1:3).*[3 2 1]);
Vfun = @(tq) ppval(ppd, min(max(tq, br(1)), br(end)));
Lfun = @(tq) ppval(pp, min(max(tq, br(1)), br(end)));
end
