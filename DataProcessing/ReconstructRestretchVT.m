i = 5;

% plot(AllFiles{i}.t, 2*AllFiles{i}.L);
clf;
hold on;
plot(AllFiles{i}.t, AllFiles{i}.L);
plot(AllFiles{i}.t(1:end -1), diff(AllFiles{i}.L));
plot(AllFiles{i}.t(1:end -2), diff(diff(AllFiles{i}.L)));

%%
% Assumptions: 
% L is your actuator signal (displacement/length)

% Active
vtmode = 'active';i = 4; x_act = 7.85 + [-5 75 85 65 -14 35 65 50]*1e-3;

% Relaxed
vtmode = 'relaxed';i = 5; x_act = 0;

t = AllFiles{i}.t/1000; % s
L = AllFiles{i}.L(t<2000);
t = t(t<2000);

dt = 1e-3;
% dt is your sampling time

% 1. Calculate the derivative (velocity)
v = diff(L) ./ dt; 
v = [v(1); v]; % Pad to keep length N
clf;hold on;
plot(t, L-1)
plot(t, [diff(L);0]/dt);
% plot(t, [diff(diff(L));0;0]/dt);
v1 = 15;
v2 = 45;
v3 = 45 + 0*37.55;
v4 = 0.015;



slack_durs = [0 5 10 20 50 100 200 500]*1e-3;
% slack_durs = [0 5 10 20]*1e-3;

offs = x_act + [0 0 49 173 223 293 319 407]*1e-3;

% slack_durs = [];
velocitytable = [0 ,0;
                14.150, -v1;
                14.150 + 0.15/v1, 0;
                24.120 + 0,     v4;
                24.120 + 0.15/v4, 0;
];
d = 0;
for i_slack_dur = 1:length(slack_durs)
    % cs = sum(slack_durs(1:max(1, i_slack_dur-4)));
    offset = (171.944)*(i_slack_dur-1) + offs(i_slack_dur);
    % offs(i_slack_dur) = offset;
    d = slack_durs(i_slack_dur);
md = 0;
vt = [
    108.416,        v2;
    108.416 + 0.2255/v2, 0;
    138.416, -v3;
    138.416 + 0.2255/v3, 0;
    d + 138.416 + 0.2255/v3, v3;
    d + 138.416 + 2*0.2255/v3, 0;
    md + 168.416, -v3;
    md + 168.416 + 0.3755/v3, 0;
    md + 188.416, v4;
    md + 188.416 + 0.15/v4, 0;
    md + 200, 0;    ];
    vt(:, 1) = vt(:, 1) + offset;
    velocitytable = [velocitytable;vt];
end

% --- Input Parameters ---
L0 = 0.95; 
time = velocitytable(:, 1);
velocity = velocitytable(:, 2);
% 2. Calculate displacement (dL) for each segment
% We use the velocity at the START of the interval
dT = diff(time);
dL = velocity(1:end-1) .* dT;

% 3. Integrate to find Position (L)
% We pad with a 0 at the start to align with L0
lengthtable = L0 + [0; cumsum(dL)];
velocitytable = [velocitytable,lengthtable];

% --- Plotting ---
figure(2);clf;
% subplot(2,1,1);

plot(time, lengthtable, 'r-o', 'MarkerFaceColor', 'r', 'MarkerSize', 4); hold on;
L_interp = interp1(t, L, time);
plot(time, L_interp, 'b-o', 'MarkerFaceColor', 'r', 'MarkerSize', 4); hold on;
plot(t, L, 'k', 'LineWidth', 1);
ylabel('Position L(t)');
title('Reconstructed Position (Actuator Signal)');
grid on;
% figure(3); clf;
% plot(time, [0 diff(Length)])

% save the vt
if strcmp(vtmode,'active')
    fn = '..\data\velocitytable_doubleramp_active.csv';
    
elseif strcmp(vtmode,'relaxed')
    fn = '..\data\velocitytable_doubleramp_relaxed.csv';
    writematrix(velocitytable, '..\data\velocitytable_doubleramp_relaxed.csv');
else
    disp("I dunno, exiting")
    return;
end
writematrix(["Time" "Velocity" "ML"], fn)
writematrix(velocitytable, fn, WriteMode="append");

return;
%%
% 2. Define a threshold to distinguish noise from actual movement
% 10% of the expected constant velocity is usually a safe bet
threshold = 0.5 * max(abs(v));

% 3. Identify movement states
%  1 = Increasing Ramp
% -1 = Decreasing Ramp
%  0 = Steady / Plateau
states = zeros(size(v));
states(v > threshold) = 1;
states(v < -threshold) = -1;

% 4. Find the timing of state changes
changes = diff(states);

% Find indices where movement starts/ends
% Start of Up-Ramp: state goes from 0 to 1
up_start_idx = find(changes == 1);
% End of Up-Ramp: state goes from 1 to 0
up_end_idx   = find(changes == -1 & states(1:end-1) == 1);

% Start of Down-Ramp: state goes from 0 to -1
down_start_idx = find(changes == -1 & states(1:end-1) == 0);
% End of Down-Ramp: state goes from -1 to 0
down_end_idx   = find(changes == 1 & states(1:end-1) == -1);

% Convert indices to time
t_up_start = t(up_start_idx);
t_down_start = t(down_start_idx);

%%

% --- Setup and Noise Handling ---
L_smooth = movmean(L, 1); % Light smoothing
v_raw = [0; diff(L_smooth) ./ diff(t)]; % Calculate velocity

% Define threshold (adjust based on noise floor)
v_threshold = 0.5 * median(abs(v_raw(abs(v_raw) > 0))); 

% --- Logic to Detect Ramp Regions ---
is_up = v_raw > v_threshold;
is_down = v_raw < -v_threshold;

% Find transitions (1 = start/end)
up_starts = find(diff([0; is_up]) == 1);
up_ends   = find(diff([is_up; 0]) == -1);
down_starts = find(diff([0; is_down]) == 1);
down_ends   = find(diff([is_down; 0]) == -1);

% --- Visualization Check ---
figure;
subplot(2,1,1);
plot(t, L, 'k', 'LineWidth', 1.5); hold on;
plot(t(up_starts), L(up_starts), 'go', 'DisplayName', 'Start Up');
plot(t(down_starts), L(down_starts), 'ro', 'DisplayName', 'Start Down');
title('Actuator Signal (L) with Detected Ramp Timings');
legend; grid on;

subplot(2,1,2);
plot(t, v_raw, 'b');
title('Calculated Velocity (v)');
ylabel('v = dL/dt'); grid on;

%% Prefilter the data first
datatable = readtable('..\Data\2025 09 19 Export\06 Log Double Ramps Relax PNB Mava.txt');    
datatable.Properties.VariableNames = {'Time', 'L','F'};
datatable.Time = datatables{1}.Time /1000;% convert to ms

% get the slacks
i_slackZones = find(velocitytable(:, 2) == 0 & velocitytable(:, 3) < 0.9);
slackZones = [i_slackZones, i_slackZones+1];
