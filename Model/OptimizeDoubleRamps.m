%% OptimizeDoubleRamps.m
% fminsearch optimization for double-ramp refolding data.
% Follows the pattern of OptimizeCOmbined.m but targets the doubleramp2
% protocol (ramp->slack->ramp) from Data/2025 11 21 Export/.
%
% Sequential strategy:
%   Phase 1 - pCa=11 (low-Ca/relaxed): fit all mechanical params + alphaF_0
%   Phase 2 - add pCa=4.51 (high-Ca/active): fit Ca-specific modifiers
%   Phase 3 - fine-tune refolding: mu, kA, kD, alphaF_0 for both conditions

close all
clear;

paramNames = {'\F_{ss}', 'n_{ss}', 'k_p', 'n_p', 'k_d', 'n_d', ...
              '\alpha_U', 'n_U', '\mu', '\delta_U', 'k_A', 'k_D', '\alpha_{F,0}'};

% pCa=11 baseline (row 6 of RunCombinedModel paramSet) + alphaF_0 initial guess
% alphaF_0 MUST be > 0 so cycles 2-8 can refold and produce peaks
params_low  = [3.884	15.19	199.752	5.43	31856	2.798	8.151E+06	13.168	0.0685	0.197	NaN	NaN	1733.86];

% pCa=4.51 baseline (row 1 of RunCombinedModel paramSet) + alphaF_0 initial guess
params_high = [3.884	15.19	224.948	5.43	31856	2.798	8.359E+06	13.4455	0.063	0.197	0.007139	0.201	1866.12];

options = optimset('Display', 'iter', 'TolFun', 1e-3, 'TolX', 0.01, ...
                   'PlotFcns', @optimplotfval, 'MaxIter', 100);

%% Quick eval to verify finite cost before optimizing
disp('--- Baseline cost (pCa=11) ---')
evalDoubleRamps([], params_low, params_high, [], [11])
disp('--- Baseline cost (pCa=4.51) ---')
evalDoubleRamps([], params_low, params_high, [], [4.51])

%% Phase 1: pCa=11, fit mechanical params + alphaF_0
% kA(11) and kD(12) are NaN for pCa=11 so they are excluded
runOptim1 = true;
if runOptim1
    modSel = [1 2 3 4 5 6 7 8 9 10 13];   % all non-NaN params + alphaF_0
    pCas   = [11];
    init   = params_low(modSel);
    evalLin = @(x) evalDoubleRamps(x, params_low, params_high, modSel, pCas);
    x = fminsearch(evalLin, init, options);
    params_low(modSel) = x;
    comment = 'Phase 1: pCa=11 low-Ca double-ramp, mechanical params + alphaF_0';
    save optresDoubleRamps_ph1 params_low params_high modSel comment
end

%% Phase 2: pCa=4.51 + pCa=11, Ca-specific modifiers on params_high
% params_low is kept fixed from Phase 1
runOptim2 = false;
if runOptim2
    modSel = [3 7 8 11 12];   % kp, alphaU, nU, kA, kD  (Ca-dependent)
    %pCas   = [4.51 11];
    pCas   = [4.51];
    init   = params_high(modSel);
    evalLin = @(x) evalDoubleRamps(x, params_low, params_high, modSel, pCas);
    x = fminsearch(evalLin, init, options);
    params_high(modSel) = x;
    comment = 'Phase 2: pCa=4.51+11, Ca-specific kp/alphaU/nU/kA/kD';
    save optresDoubleRamps_ph2 params_low params_high modSel comment
end

%% Phase 3: refolding params shared between Ca conditions
runOptim3 = false;
if runOptim3
    modSel = [9 11 12 13];   % mu, kA, kD, alphaF_0
    pCas   = [4.51 11];
    init   = params_high(modSel);
    evalLin = @(x) evalDoubleRamps(x, params_low, params_high, modSel, pCas);
    x = fminsearch(evalLin, init, options);
    params_high(modSel) = x;
    params_low(modSel)  = x;   % mu and alphaF_0 are Ca-independent
    comment = 'Phase 3: refolding params mu/kA/kD/alphaF_0 both Ca conditions';
    save optresDoubleRamps_ph3 params_low params_high modSel comment
end

%% Print current params
fprintf('\n--- params_low (pCa=11) ---\n')
for i = 1:length(paramNames)
    fprintf('%2d) %-16s = %g\n', i, paramNames{i}, params_low(i));
end
fprintf('\n--- params_high (pCa=4.51) ---\n')
for i = 1:length(paramNames)
    fprintf('%2d) %-16s = %g\n', i, paramNames{i}, params_high(i));
end

%% =========================================================================
%  Functions
%% =========================================================================

function totalCost = evalDoubleRamps(optMods, params_low, params_high, modSel, pCas)
% Top-level objective. Applies optMods to both param sets at modSel indices,
% then sums cost over requested Ca conditions.

    if ~isempty(modSel) && ~isempty(optMods)
        params_low(modSel)  = optMods;
        params_high(modSel) = optMods;
    end

    totalCost = 0;

    if ismember(11, pCas)
        cost = isolateRunDoubleRamps(params_low, 11);
        fprintf('  pCa=11   cost = %g\n', cost);
        totalCost = totalCost + cost;
    end

    if ismember(4.51, pCas) || ismember(4.4, pCas)
        cost = isolateRunDoubleRamps(params_high, 4.51);
        fprintf('  pCa=4.51 cost = %g\n', cost);
        totalCost = totalCost + cost;
    end
end

function cost = isolateRunDoubleRamps(params, pCa)
% Isolated wrapper: calls RunCombinedModel with doubleramp2 simtype,
% then computes a peak-based cost (first and second peaks per slack duration).
%
% Key design choices:
%   - Force all vectors to column (:) before arithmetic to avoid
%     outer-product when row x column subtraction creates a matrix.
%   - Use thresholds relative to data amplitude so the cost stays
%     well-behaved while fminsearch explores parameter space.
%   - Filter data peaks to positive values (relaxed data has a negative
%     baseline from drift/offset that creates spurious negative peaks).

    drawPlots = false;
    rampSet   = [1];
    alphaF_0  = params(13);

    if pCa >= 10
        simtype = 'velocitytable_doubleramp2_relaxed';
    else
        simtype = 'velocitytable_doubleramp2_active';
    end

    RunCombinedModel;   % → Force{1}, Time{1}, datatables{1}

    % Guard: simulation failed to produce output
    if isempty(Force) || isempty(Force{1})
        cost = 1e6;
        return;
    end

    % Force column vectors — Force{1} is a row vector (built by indexing in a
    % loop), while datatables{1}.F is a column vector from readtable.
    % Mixing orientations causes row-column outer-product instead of
    % element-wise subtraction, making cost a matrix instead of scalar.
    sim_force = Force{1}(:);
    sim_time  = Time{1}(:);

    % Guard: negative params or diverged solution
    if any(isnan(sim_force)) || max(sim_force) <= 0
        cost = 1e6;
        return;
    end

    % Thresholds relative to sim amplitude for robust peak detection
    sim_prom = 0.05 * max(sim_force);

    sim_peaks = findpeaks(sim_force, sim_time, ...
        'MinPeakWidth',       0.5e-3, ...
        'MaxPeakWidth',       200,    ...
        'MinPeakProminence',  sim_prom, ...
        'MinPeakDistance',    20);
    sim_peaks = sim_peaks(:);   % ensure column

    if length(sim_peaks) < 2
        cost = 1e6;
        return;
    end
    sim_first  = sim_peaks(1:2:end);
    sim_second = sim_peaks(2:2:end);

    % Data peaks (datatables{1} loaded by RunCombinedModel)
    if isempty(datatables) || isempty(datatables{1})
        cost = 1e6;
        return;
    end
    data_force = datatables{1}.F(:);
    data_time  = datatables{1}.Time(:);

    % Threshold relative to data amplitude
    data_prom = 0.2 * max(data_force);

    data_peaks = findpeaks(data_force, data_time, ...
        'MinPeakWidth',       1.5e-3, ...
        'MaxPeakWidth',       5,      ...
        'MinPeakProminence',  data_prom, ...
        'MinPeakDistance',    20);
    data_peaks = data_peaks(:);

    % Relaxed data has a negative DC offset — keep only positive peaks
    data_peaks = data_peaks(data_peaks > 0);

    if length(data_peaks) < 2
        cost = 1e6;
        return;
    end
    data_first  = data_peaks(1:2:end);
    data_second = data_peaks(2:2:end);

    % Normalized peak cost — scale matched to En{j} in RunCombinedModel
    n   = min([length(sim_first), length(data_first), ...
               length(sim_second), length(data_second)]);
    if n < 1
        cost = 1e6;
        return;
    end
    ref  = max(data_first);
    err1 = (sim_first(1:n)  - data_first(1:n) ) / ref;
    err2 = (sim_second(1:n) - data_second(1:n)) / ref;
    cost = 1e3 * mean(err1.^2 + err2.^2);

    if ~isfinite(cost)
        cost = 1e6;
    end
end
