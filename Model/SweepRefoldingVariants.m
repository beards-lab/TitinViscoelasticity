%% SweepRefoldingVariants.m
% ENTRY POINT. Two modified refolding functions on the whole double-ramp
% protocol (as SweepRefolding.m: best stretch-hold fit + own
% sensor, all other parameters fixed), each a 2-D sweep:
%   'nDep'  n-dependent rate, on/off mask kept (simStretchHold index 28):
%           n+1 -> n at alphaF_0*((n+1)/Ng)^gammaF      grid alphaF_0 x gammaF
%   'fGate' smooth force gate instead of the mask (index 26):
%           n+1 -> n at alphaF_0*exp(-Fp(s,n+1)/F_R)    grid alphaF_0 x F_R
% Workspace config (optional): cond ('Relax' default | 'Active'), variants
% (cellstr, default both).
% Run as: batch('SweepRefoldingVariants', 'Pool', 8) from Model/;
% (active: batch(..., 'Workspace', struct('cond', 'Active'))); per-simulation
% results in refoldVariants[_Active]_parts/ (restartable), collected in
% refoldVariants[_Active].mat, plotted by Figures/PlotRefoldingVariants.m.

if ~exist('cond', 'var'), cond = 'Relax'; end
if ~exist('variants', 'var'), variants = {'nDep', 'fGate'}; end
tag = ''; if ~strcmp(cond, 'Relax'), tag = ['_' cond]; end
dd = '..\Data\2025 11 21 Export\';
rf = {'0ms', '5ms', '10ms', '50ms', '100ms', '1s', '10s', '30s'};
V.nDep  = struct('idx', 28, 'name', '\gamma_F', 'vals', [0.5 1 2 3 4 6], ...
                 'aGrid', logspace(-1, 4, 16));
V.fGate = struct('idx', 26, 'name', 'F_R (kPa)', 'vals', [0.01 0.03 0.1 0.3 1 3 10], ...
                 'aGrid', logspace(-2, 3, 16));

B = load(['fitStretchHold_' cond '_best.mat']);
p0 = [B.p, NaN(1, 28 - numel(B.p))]; sen = B.sensor;
tol = struct('Relax', 1e-6, 'Active', 1e-5);        % as in the fits' notes
o = odeset('RelTol', tol.(cond), 'AbsTol', tol.(cond));
SS = cellfun(@(r) loadProtocol([dd r '_refolding_' cond '.txt']), rf, 'UniformOutput', false);

pd = ['refoldVariants' tag '_parts'];
if ~exist(pd, 'dir'), mkdir(pd); end
jobs = zeros(0, 4);                                  % [variant, value, alpha, gap]
for v = 1:numel(variants)
    X = V.(variants{v});
    [gi, ai, vi] = ndgrid(1:numel(rf), 1:numel(X.aGrid), 1:numel(X.vals));
    jobs = [jobs; repmat(v, numel(gi), 1), vi(:), ai(:), gi(:)]; %#ok<AGROW>
end
pf = @(j) fullfile(pd, sprintf('%s_v%02d_a%02d_%s.mat', variants{jobs(j, 1)}, ...
    jobs(j, 2), jobs(j, 3), rf{jobs(j, 4)}));
todo = find(~arrayfun(@(j) exist(pf(j), 'file') == 2, 1:size(jobs, 1)));
fprintf('%d of %d simulations to run\n', numel(todo), size(jobs, 1));
Vs = cellfun(@(n) V.(n), variants);
parfor q = 1:numel(todo)
    j = todo(q); X = Vs(jobs(j, 1));
    p = p0; p(13) = X.aGrid(jobs(j, 3)); p(X.idx) = X.vals(jobs(j, 2));
    [Fb, tf, Ff, Fr] = modelBinned(p, SS{jobs(j, 4)}, sen, o);
    parsave(pf(j), Fb, tf, Ff, Fr);
end

res = struct();
for v = 1:numel(variants)
    X = V.(variants{v});
    sz = [numel(rf), numel(X.aGrid), numel(X.vals)];
    Fb = cell(sz); tf = Fb; Ff = Fb; Fr = Fb;
    for j = find(jobs(:, 1) == v)'
        Y = load(pf(j)); k = sub2ind(sz, jobs(j, 4), jobs(j, 3), jobs(j, 2));
        Fb{k} = Y.Fb; tf{k} = Y.tf; Ff{k} = Y.Ff; Fr{k} = Y.Fr;
    end
    res.(variants{v}) = struct('def', X, 'Fb', {Fb}, 'tf', {tf}, 'Ff', {Ff}, 'Fr', {Fr});
end
for k = 1:numel(SS)
    SS{k} = rmfield(SS{k}, {'A', 'tSim', 'ramp', 't', 'F', 'Lraw'});
end
save(['refoldVariants' tag '.mat'], 'res', 'SS', 'p0', 'sen', 'rf', 'variants', 'cond');

function parsave(f, Fb, tf, Ff, Fr)
save(f, 'Fb', 'tf', 'Ff', 'Fr');
end
