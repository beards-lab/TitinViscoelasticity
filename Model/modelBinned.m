function [Fb, tf, Ff] = modelBinned(params, S, sensor, odeOpts)
% modelBinned  Model force on the log bins of a loadStretchHold dataset S,
% optionally seen through the force transducer (sensor = [f0 zeta], see
% sensorFilter.m). Utility for FitFirstStretch.m.
%   S.filtWin  optional rows [t1 t2] of fast windows where the sensor
%              dynamics matter (default: [0, end of ramp + 15 ms])
%   tf, Ff     uniform fine-grid model force (filtered) in those windows
if nargin < 3, sensor = []; end
if nargin < 4, odeOpts = []; end % [] = simStretchHold default (1e-4)
tf = []; Ff = [];
if isempty(sensor)
    F = simStretchHold(params, S.tSim, S.ramp, odeOpts);
else
    dtf = 1e-5;
    if isfield(S, 'filtWin')
        win = S.filtWin;
    else
        win = [0, S.ramp.tEnd + 15e-3];  % ringing has died out well before
    end
    tw = cell(size(win, 1), 1);
    for w = 1:size(win, 1)
        tw{w} = (win(w, 1):dtf:win(w, 2))';
    end
    tAll = unique([vertcat(tw{:}); S.tSim]);
    Fall = simStretchHold(params, tAll, S.ramp, odeOpts);
    if any(~isfinite(Fall))
        Fb = nan(size(S.Fb)); return;
    end
    F = interp1(tAll, Fall, S.tSim);
    tf = cell(size(tw)); Ff = tf;
    for w = 1:numel(tw)
        tf{w} = tw{w};
        Ff{w} = sensorFilter(interp1(tAll, Fall, tw{w}), dtf, sensor(1), sensor(2));
        in = S.tSim >= tw{w}(1) & S.tSim <= tw{w}(end);
        F(in) = interp1(tw{w}, Ff{w}, S.tSim(in));
    end
    if numel(tf) == 1
        tf = tf{1}; Ff = Ff{1};
    end
end
Fb = S.A*F;
end
