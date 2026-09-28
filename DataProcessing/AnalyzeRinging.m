%% AnalyzeRinging.m
% ENTRY POINT. Identifies the force-transducer resonance that rings after the
% ~2.8 ms first stretch of the 2025-11-21 refolding protocol, and shows the
% pooled force with the transducer response removed (Model/deconvSensor.m).
%
% Pooling: the 8 *_refolding_<cond>.txt repeats share the first stretch
% (0.95 -> 1.175 L0), are phase-locked, and are averaged sample by sample.
% Sensor model = mass-spring-damper transducer (Model/sensorFilter.m), one
% [f0 zeta] shared by all conditions, chosen so the deconvolved force has
% minimal 1.1-2.3 kHz energy 0.2-9 ms after the peak.
% Result (2026-09-27): f0 ~1570 Hz, zeta ~0.10 in relaxed, active and at
% 3/10/30 C (ring energy -94..-99 %); the frequency does not follow fibre
% stiffness, so the fibre is not part of the resonator.
%
% Workspace config: conds (cell, default {'Relax', 'Active'}).

addpath('../Model');
if ~exist('conds', 'var'), conds = {'Relax', 'Active'}; end
dd = '../Data/2025 11 21 Export/';
rf = {'0ms', '5ms', '10ms', '50ms', '100ms', '1s', '10s', '30s'};
dt = 1e-4; tg = (-0.02:dt:0.1)';

%% Pool the first stretch of the 8 repeats
Fp = zeros(numel(tg), numel(conds)); Lp = Fp;
for c = 1:numel(conds)
    for r = 1:numel(rf)
        e = readtable([dd rf{r} '_refolding_' conds{c} '.txt'], 'VariableNamingRule', 'preserve');
        e.Properties.VariableNames = {'Time', 'L', 'F', 'SL'};
        t = e.Time/1000;
        v = gradient(movmean(e.L, 5), t);
        t0 = t(find(abs(v) > 1 & t > 9 & t < 11, 1));
        F0 = mean(e.F(t > t0 - 0.3 & t < t0 - 0.05));
        Fp(:, c) = Fp(:, c) + interp1(t - t0, e.F - F0, tg)/numel(rf);
        Lp(:, c) = Lp(:, c) + interp1(t - t0, e.L, tg)/numel(rf);
    end
end

%% Identify the shared transducer [f0 zeta]
[bb, ab] = butter(3, [1100 2300]/(0.5/dt));
[~, iPk] = max(Fp(tg < 0.01, :));
ringE = @(x) sum(arrayfun(@(c) sum(filtfilt(bb, ab, ...
    deconvSensor(Fp(:, c), dt, x(1), x(2))).^2 .* ...
    (tg > tg(iPk(c)) + 2e-4 & tg < tg(iPk(c)) + 9e-3)), 1:numel(conds)));
sen = fminsearch(ringE, [1600 0.1]);
fprintf('shared transducer: f0 = %.0f Hz, zeta = %.3f\n', sen);
raw = sum(arrayfun(@(c) sum(filtfilt(bb, ab, Fp(:, c)).^2 .* ...
    (tg > tg(iPk(c)) + 2e-4 & tg < tg(iPk(c)) + 9e-3)), 1:numel(conds)));
fprintf('ring energy removed: %.0f %%\n', 100*(1 - ringE(sen)/raw));

%% Measured vs deconvolved
f = figure(801); clf; f.Position = [100 100 1200 420];
tiledlayout(1, numel(conds), 'TileSpacing', 'compact');
for c = 1:numel(conds)
    Fd = deconvSensor(Fp(:, c), dt, sen(1), sen(2));
    v = gradient(Lp(:, c), dt);
    [pk, i] = max(Fd(tg < 0.01));
    fprintf('%s: peak measured %.1f, deconvolved %.1f kPa; drop in 0.5 ms after peak %.1f kPa\n', ...
        conds{c}, max(Fp(tg < 0.01, c)), pk, pk - Fd(i + round(5e-4/dt)));
    nexttile; hold on; box on; grid on;
    plot(1e3*tg, Fp(:, c), 'Color', [.6 .6 .6]);
    plot(1e3*tg, Fd, 'k', 'LineWidth', 1.5);
    plot(1e3*tg, v/max(v)*0.3*pk, 'g');
    xlim([-1 12]); xlabel('t from stretch onset (ms)'); ylabel('\Theta (kPa)');
    title(sprintf('%s, 8-repeat mean', conds{c}));
    % round trip: the transducer applied to the deconvolved force gives the
    % measured force back (minus > 3 kHz content cut by deconvSensor)
    Frt = sensorFilter(Fd, dt, sen(1), sen(2));
    plot(1e3*tg, Frt, 'r--');
    fprintf('%s: round trip rms(measured - sensorFilter(deconvolved)) %.2f kPa (noise %.2f)\n', ...
        conds{c}, rms(Fp(tg > -1e-3 & tg < 0.012, c) - Frt(tg > -1e-3 & tg < 0.012)), std(Fp(tg < -1e-3, c)));
    legend('measured', sprintf('deconvolved (f_0 %.0f Hz, \\zeta %.2f)', sen), ...
        'velocity (scaled)', 'deconvolved \rightarrow sensorFilter');
end

%% Model vs data: conv(model) vs raw data  <=>  model vs deconv(data)
% Best parameter sets (Model/fitStretchHold_<cond>_best.mat) were fitted with
% their own free sensor; here the same model force is seen through (a) that
% sensor, (b) the shared sensor identified above, (c) no sensor, compared to
% the deconvolved data. Fitting (b) to raw data and (c) to deconvolved data
% have the same residual up to the filter H^-1, which amplifies noise near
% and above f0 - so fits keep using (b) (modelBinned.m).
w = tg > -1e-3 & tg < 0.012;          % ramp + ringing window
f = figure(802); clf; f.Position = [100 100 1200 750];
tiledlayout(2, numel(conds), 'TileSpacing', 'compact');
for c = 1:numel(conds)
    R = load(['../Model/fitStretchHold_' conds{c} '_best.mat']);
    S = loadStretchHold(strcat(dd, rf, ['_refolding_' conds{c} '.txt']), [], 200);
    dtf = 1e-5; tf = (0:dtf:0.02)';
    Fm = simStretchHold(R.p, tf, S.ramp, odeset('RelTol', 1e-6, 'AbsTol', 1e-6));
    % model is at rest before onset (0 kPa); resample on the data grid
    toGrid = @(y) interp1([-0.03; tf], [0; y], tg, 'linear', y(end));
    Mraw = toGrid(Fm);
    Mfit = toGrid(sensorFilter(Fm, dtf, R.sensor(1), R.sensor(2)));
    Msen = toGrid(sensorFilter(Fm, dtf, sen(1), sen(2)));
    Fd = deconvSensor(Fp(:, c), dt, sen(1), sen(2));
    e = @(a, b) rms(a(w) - b(w));
    fprintf(['%s (0-12 ms rms, kPa): raw vs model+fit sensor [%.0f %.3f] %.2f | ' ...
        'raw vs model+shared sensor %.2f | deconv vs model %.2f | raw vs model %.2f\n'], ...
        conds{c}, R.sensor, e(Fp(:, c), Mfit), e(Fp(:, c), Msen), e(Fd, Mraw), e(Fp(:, c), Mraw));
    nexttile(c); hold on; box on; grid on;
    plot(1e3*tg, Fp(:, c), 'Color', [.6 .6 .6]);
    plot(1e3*tg, Mfit, 'b', 'LineWidth', 1.2); plot(1e3*tg, Msen, 'r', 'LineWidth', 1.2);
    xlim([-1 12]); ylabel('\Theta (kPa)'); title(sprintf('%s: raw data vs model through sensor', conds{c}));
    legend('measured', sprintf('model + fit sensor [%.0f %.2f]', R.sensor), ...
        sprintf('model + shared sensor [%.0f %.2f]', sen));
    nexttile(numel(conds) + c); hold on; box on; grid on;
    plot(1e3*tg, Fd, 'k', 'LineWidth', 1.2); plot(1e3*tg, Mraw, 'm', 'LineWidth', 1.2);
    xlim([-1 12]); xlabel('t from stretch onset (ms)'); ylabel('\Theta (kPa)');
    title(sprintf('%s: deconvolved data vs model (no sensor)', conds{c}));
    legend('deconvolved data', 'model force');
end
