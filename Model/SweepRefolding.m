%% SweepRefolding.m
% ENTRY POINT. Runs the whole 2025-11-21 refolding protocol (first stretch,
% 30 s hold, release, gap 0 ms..30 s, restretch, hold, final release;
% Model/loadProtocol.m) with the current best stretch-hold fits
% (fitStretchHold_<cond>_best.mat, each seen through its own fitted sensor,
% sensorFilter.m) and sweeps ONLY the refolding rate alphaF_0 (params(13);
% dXdT on/off refolding: n+1 -> n only while the proximal segment of state
% n+1 is slack). Everything else, incl. ODE tolerances of the fits, is fixed.
% Workspace config (optional): conds, aGrid.
% Run as: batch('SweepRefolding', 'Pool', 8) from Model/;
% per-simulation results in refoldSweep_parts/ (a rerun skips finished ones),
% collected in refoldSweep.mat, plotted by Figures/PlotRefoldingSweep.m.

if ~exist('conds', 'var'), conds = {'Relax', 'Active'}; end
if ~exist('aGrid', 'var'), aGrid = [0, logspace(-3, 3, 13)]; end
dd = '..\Data\2025 11 21 Export\';
rf = {'0ms', '5ms', '10ms', '50ms', '100ms', '1s', '10s', '30s'};
gap = [0 5e-3 10e-3 50e-3 0.1 1 10 30];
tol = struct('Relax', 1e-6, 'Active', 1e-5);    % as in the fits' notes

SS = cell(numel(conds), numel(rf)); P = cell(1, numel(conds)); sen = P;
for c = 1:numel(conds)
    B = load(['fitStretchHold_' conds{c} '_best.mat']);
    P{c} = B.p; sen{c} = B.sensor;
    for g = 1:numel(rf)
        SS{c, g} = loadProtocol([dd rf{g} '_refolding_' conds{c} '.txt']);
    end
end

[ci, gi, ai] = ndgrid(1:numel(conds), 1:numel(rf), 1:numel(aGrid));
nJ = numel(ci);
pd = 'refoldSweep_parts';     % one file per simulation: restartable
if ~exist(pd, 'dir'), mkdir(pd); end
pf = @(j) fullfile(pd, sprintf('%s_%s_a%02d.mat', conds{ci(j)}, rf{gi(j)}, ai(j)));
todo = find(~arrayfun(@(j) exist(pf(j), 'file') == 2, 1:nJ));
fprintf('%d of %d simulations to run\n', numel(todo), nJ);
parfor q = 1:numel(todo)
    j = todo(q);
    p = P{ci(j)}; p(13) = aGrid(ai(j));
    o = odeset('RelTol', tol.(conds{ci(j)}), 'AbsTol', tol.(conds{ci(j)}));
    [Fb, tf, Ff, Fr] = modelBinned(p, SS{ci(j), gi(j)}, sen{ci(j)}, o);
    parsave(pf(j), Fb, tf, Ff, Fr);
end
Fb = cell(size(ci)); tf = Fb; Ff = Fb; Fr = Fb;
for j = 1:nJ
    X = load(pf(j)); Fb{j} = X.Fb; tf{j} = X.tf; Ff{j} = X.Ff; Fr{j} = X.Fr;
end
for k = 1:numel(SS)   % keep bins and events; raw traces are reloaded for plots
    SS{k} = rmfield(SS{k}, {'A', 'tSim', 'ramp', 't', 'F', 'Lraw'});
end
save('refoldSweep.mat', 'Fb', 'tf', 'Ff', 'Fr', 'SS', 'P', 'sen', 'aGrid', 'conds', 'rf', 'gap');

function parsave(f, Fb, tf, Ff, Fr)
save(f, 'Fb', 'tf', 'Ff', 'Fr');
end
