vtmode = 'active';
% vtmode = 'relaxed';
%% signal load

if strcmp(vtmode,'active')
    datatable = readtable('../data/2025 11 21 Export/05_Log_Active_Refolding.txt', 'filetype', 'text', 'NumHeaderLines',4);
elseif strcmp(vtmode,'relaxed')
    datatable = readtable('../data/2025 11 21 Export/03_Log_Relax_Refolding.txt', 'filetype', 'text', 'NumHeaderLines',4);
end

if length(datatable.Properties.VariableNames) == 3
    datatable.Properties.VariableNames = {'t', 'L','F'};
elseif length(datatable.Properties.VariableNames) == 4
    datatable.Properties.VariableNames = {'t', 'L','F', 'SL'};
else
    disp('Wat?')
end

t = datatable.t/1000;
L = datatable.L;
F = datatable.F;
    
%% --- 1. Signal Preparation ---

V = diff(L) ./ diff(t);
t_v = t(1:end-1);

% Smooth velocity slightly to handle noise (window of 5-10 samples)
V_smooth = movmean(V, 5); 

% --- 2. Separate Peak Finding ---
% Adjust 'MinPeakHeight' and 'MinPeakDistance' based on your noise levels
% Positive Jumps (v1 and v3)
[pks_pos, locs_pos] = findpeaks(V_smooth, ...
    'MinPeakHeight', 20, ... 
    'MinPeakDistance', 50); 

% Negative Jumps (v2 and v4) - we use -V_smooth to find dips as peaks
[pks_neg, locs_neg] = findpeaks(-V_smooth, ...
    'MinPeakHeight', 30, ...
    'MinPeakDistance', 50);

% --- 3. Map to specific positions ---
% Assuming 5 repeats, we expect 10 pos peaks and 10 neg peaks total
i_p1 = locs_pos(1:2:end); 
i_p3 = locs_pos(2:2:end);
i_p2 = locs_neg(1:2:end);
i_p4 = locs_neg(2:2:end);

% --- 4. Visualization (Verification Step) ---
figure('Color', 'w', 'Name', 'Peak Detection Verification');
ax(1) = subplot(2,1,1);
plot(t, L, 'k'); hold on;
plot(t_v(i_p1), L(i_p1), 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Start v1');
plot(t_v(i_p3), L(i_p3), 'mo', 'MarkerFaceColor', 'm', 'DisplayName', 'Start v3');
plot(t_v(i_p2), L(i_p2), 'bo', 'MarkerFaceColor', 'b', 'DisplayName', 'Start v2');
plot(t_v(i_p4), L(i_p4), 'co', 'MarkerFaceColor', 'c', 'DisplayName', 'Start v4');
ylabel('Length (L)'); legend('show'); grid on;
title('Position Data with Detected Jump Starts');

ax(2) = subplot(2,1,2);
plot(t_v, V, 'Color', [0.7 0.7 0.7]); hold on;
plot(t_v, V_smooth, 'r', 'LineWidth', 1);
ylabel('Velocity (V)'); grid on;
title('Smoothed Velocity for Peak Detection');

linkaxes(ax, 'x');

%% --- 2. Fill Velocity Table with Loops ---

num_repeats = length(i_p1);
% baseline to pos1
v1 = 100;
% pos1 to pos2
v2 = -150;
% pos2 to pos3
v3 = 100;
% pos3 to pos4
v4 = -40;
% pos4 to pos5
v5 = 0.02;

baseline = 0.95;
pos1 = 1.175;
pos2 = baseline;
pos3 = pos1;
pos4 = 0.8;
pos5 = baseline;

vt = [t(i_p1(1)) - 10, 0];

for i = 1:num_repeats
    % Movement 1: Baseline to Pos1
    t1_start = t(i_p1(i));
    dt1 = abs(baseline-pos1) / abs(v1);
    vt = [vt; t1_start, v1; t1_start + dt1, 0];
    
    % Movement 2: Pos1 to Pos2
    t2_start = t(i_p2(i));
    dt2 = abs(pos2 - pos1) / abs(v2);
    vt = [vt; t2_start, v2; t2_start + dt2, 0];
    
    % Movement 3: Pos2 to Pos3
    t3_start = t(i_p3(i));
    dt3 = abs(pos3 - pos2) / abs(v3);
    vt = [vt; t3_start, v3; t3_start + dt3, 0];
    
    % Movement 4: Pos3 to Pos4
    t4_start = t(i_p4(i));
    dt4 = abs(pos4 - pos3) / abs(v4);
    vt = [vt; t4_start, v4; t4_start + dt4, 0];
    
    % Movement 5: Pos4 to Pos5 (Slow drift)
    % Since this is v5=0.05, we trigger it shortly after p4 ends or at a fixed offset
    t5_start = t4_start + dt4 + 10; 
    dt5 = abs(pos5 - pos4) / abs(v5);
    vt = [vt; t5_start, v5; t5_start + dt5, 0];
end

% Sort table by time (just in case)
% vt = sortrows(vt, 1);

% --- 2. Corrected Reconstruction ---
% Use 'previous' to treat the table as a command sequence (step function)
velocitytable = vt;
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

% Plotting
figure;
subplot(2,1,1);
plot(t, L, 'k', time, lengthtable, 'r--'); title('Position Reconstruction');
legend('Original', 'From Table');
subplot(2,1,2);
plot(t(1:end-1), V); hold on;
plot(vt(:,1), vt(:,2), 'ro'); title('Velocity Table Check');

%% Save

% save the vt
if strcmp(vtmode,'active')
    fn = '..\data\velocitytable_doubleramp2_active.csv';
    
elseif strcmp(vtmode,'relaxed')
    fn = '..\data\velocitytable_doubleramp2_relaxed.csv';
else
    disp("I dunno, exiting")
    return;
end
writematrix(["Time" "Velocity" "ML"], fn)
writematrix(velocitytable, fn, WriteMode="append");
