%% Get fmaxes - for each dataset there is a Fmax
% expects _dataset_, leg_names

% index: dataset, 
% measurement, sequence;
figure(1);clf;hold on;
fmax_pointers = [2, 1;2 1;2, 1;2, 1;2, 1;2,1];
for i = 1:size(dataset, 2)
    if fmax_pointers(i, 1) == 0
        % no fmax measured
        dataset_maxF(i) = NaN;
        continue;
    end
    dtst = dataset{i}.dsc{fmax_pointers(i, 1), fmax_pointers(i, 2)};
    rmp = dtst.datatable;

    smoothrmp = useMovAvg(rmp);



    % subplot(122);hold on;
    % plot(rmp.t, rmp.L);
    % subplot(121);hold on;
    pfm(i) = plot(rmp.t, rmp.F);

    %  the pCa starts about halfway the dataset
    i_start = find(rmp.t > 100, 1, 'first');
    i_firstSteps = find(rmp.L(i_start:end) > 1.001, 1, 'first') - 5 + i_start;
    i_step_dwn = min(find(rmp.L(i_start:end) < 0.95, 1, 'last') +5 + i_start, length(rmp.L));
    smoothrmp.F(i_firstSteps:i_step_dwn) = 0;
    % % we use avg now, no need for this
    % i_firstSteps = length(rmp.t)
    % rmp.t([i_start, i_firstSteps])

    % plot(rmp.t(i_start:i_firstSteps), rmp.F(i_start:i_firstSteps), 'x')
    [maxF i_max] = max(smoothrmp.F(i_start:end));
    dataset_maxF(i) = maxF;
    plot(smoothrmp.t, smoothrmp.F, '--');
    plot(rmp.t(i_max+i_start), maxF, 'rx', LineWidth=4, MarkerSize=12)
end
legend(pfm, leg_names);

%%

function rmpOut = useMovAvg(rmp)

% MATLAB Script: Moving Window Average for Signal in Table
% Assumes rmp is a table with time values in rmp.t and signal values in rmp.F.

% Parameters
window_size = .5; % Define the window size in units of rmp.t (e.g., seconds)

% Check if the table contains required columns
if ~all(ismember({'t', 'F'}, rmp.Properties.VariableNames))
    error('The table rmp must contain columns named ''t'' and ''F''.');
end

% Extract time and signal from the table
time = rmp.t;
signal = rmp.F;

% Initialize the averaged signal
smoothed_signal = zeros(size(signal));

% Apply the moving window averaging
for i = 1:length(signal)
    % Define the window range based on time
    start_time = time(i) - window_size / 2;
    end_time = time(i) + window_size / 2;
    
    % Find indices within the time window
    indices_in_window = find(time >= start_time & time <= end_time);
    
    % Compute the average in the window
    smoothed_signal(i) = mean(signal(indices_in_window));
end

rmpOut.t = time;rmpOut.F = smoothed_signal;

return
%% Plot the original and smoothed signal
figure;
plot(time, signal, 'b-', 'DisplayName', 'Original Signal'); hold on;
plot(time, smoothed_signal, 'r-', 'DisplayName', 'Smoothed Signal');
legend;
xlabel('Time');
ylabel('Signal');
title('Moving Window Averaging');
grid on;
end