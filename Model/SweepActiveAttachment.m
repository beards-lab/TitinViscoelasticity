%% SweepActiveAttachment.m
% ENTRY POINT. Titin-actin attachment cycling / release at slack on the whole
% ACTIVE double-ramp protocol. Start: high-attachment active model
% (jointFit_localtest/sweep4_highAttach.mat: ~100 % attached, slip-bond kDf,
% slack release kDslack) with its own sensor. Swept:
%   kA        reattachment rate (kD scaled with it, attached fraction fixed)
%   kDslack   release rate of attached chains whose distal segment is slack
%   refold    off, or the relaxed #1 best (alphaF_0*((n+1)/Ng)^gammaF): titin
%             domain refolding taken as shared between relaxed and active
% Run as: batch('SweepActiveAttachment', 'Pool', 8) from Model/; parts in
% activeAttach_parts/ (restartable), result activeAttach.mat.

dd = '..\Data\2025 11 21 Export\';
rf = {'0ms', '5ms', '10ms', '50ms', '100ms', '1s', '10s', '30s'};
kAg = [1e1 3e1 1e2 3e2 1e3 1e4 8.53e5];
kDsg = [0 1e2 1e5];
H = load('jointFit_localtest\sweep4_highAttach.mat');
p0 = [H.pA4, NaN(1, 28 - numel(H.pA4))]; sen = H.sA4;
fA = p0(11)/(p0(11) + p0(12));
Rv = load('refoldVariants.mat'); [aR, gR] = bestNDepRelax(Rv);
refold = [0 NaN; aR gR];                       % [alphaF_0 gammaF] per level
o = odeset('RelTol', 1e-5, 'AbsTol', 1e-5);    % as in the active fits
SS = cellfun(@(r) loadProtocol([dd r '_refolding_Active.txt']), rf, 'UniformOutput', false);

[gi, ai, di, ri] = ndgrid(1:numel(rf), 1:numel(kAg), 1:numel(kDsg), 1:2);
pd = 'activeAttach_parts'; if ~exist(pd, 'dir'), mkdir(pd); end
pf = @(j) fullfile(pd, sprintf('a%d_d%d_r%d_%s.mat', ai(j), di(j), ri(j), rf{gi(j)}));
todo = find(~arrayfun(@(j) exist(pf(j), 'file') == 2, 1:numel(gi)));
fprintf('%d of %d simulations to run (relaxed refolding %.3g /s, gamma %g)\n', numel(todo), numel(gi), aR, gR);
parfor q = 1:numel(todo)
    j = todo(q); p = p0;
    p(11) = kAg(ai(j)); p(12) = kAg(ai(j))*(1 - fA)/fA;
    p(27) = kDsg(di(j)); p(13) = refold(ri(j), 1); p(28) = refold(ri(j), 2);
    [Fb, tf, Ff, Fr] = modelBinned(p, SS{gi(j)}, sen, o);
    parsave(pf(j), Fb, tf, Ff, Fr);
end
Fb = cell(size(gi)); tf = Fb; Ff = Fb; Fr = Fb;
for j = 1:numel(gi)
    Y = load(pf(j)); Fb{j} = Y.Fb; tf{j} = Y.tf; Ff{j} = Y.Ff; Fr{j} = Y.Fr;
end
for k = 1:numel(SS), SS{k} = rmfield(SS{k}, {'A', 'tSim', 'ramp', 't', 'F', 'Lraw'}); end
save('activeAttach.mat', 'Fb', 'tf', 'Ff', 'Fr', 'SS', 'p0', 'sen', 'fA', 'kAg', 'kDsg', 'refold', 'rf');

function [a, g] = bestNDepRelax(V)
X = V.res.nDep; S = V.SS;
cost = zeros(numel(X.def.aGrid), numel(X.def.vals));
for gi = 1:numel(S)
    kR = numel(S{gi}.ev) - 1; rs = S{gi}.evBin == kR; sc = max(S{gi}.Fb(S{gi}.evBin == 1));
    for ai = 1:size(cost, 1)
        for vi = 1:size(cost, 2)
            cost(ai, vi) = cost(ai, vi) + mean(((X.Fb{gi, ai, vi}(rs) - S{gi}.Fb(rs))/sc).^2);
        end
    end
end
[~, k] = min(cost(:)); [ia, iv] = ind2sub(size(cost), k);
a = X.def.aGrid(ia); g = X.def.vals(iv);
end

function parsave(f, Fb, tf, Ff, Fr)
save(f, 'Fb', 'tf', 'Ff', 'Fr');
end
