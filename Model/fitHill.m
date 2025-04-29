function fitHill(modSet, modSel, pcax)
cutXaxis = true;
% Example data points (replace with actual values)
x_data = pcax;
y_data = modSet(:, modSel)';  % Dependent variable

filter = ~isnan(y_data);
% y_data = y_data(filter);
% x_data = x_data(filter);
y_data(~filter) = 0;

%% Choose initial parameter guesses [A, K, n]
A0 = max(y_data);     % Approximate max value
A_0 = min(y_data);
K0 = 6; %median(x_data);  % Midpoint estimate
n0 = 100;               % Initial Hill coefficient

params0 = [A0, K0, n0, A_0];

% Define Hill function models
hill_rising = @(p, x) (p(1) * x.^p(3)) ./ (p(2).^p(3) + x.^p(3)) + A_0;  % Rising Hill function
hill_decreasing = @(p, x) p(1) ./ (1 + (x./p(2)).^p(3) + A_0);           % Decreasing Hill function

% Perform curve fitting using nonlinear least squares
options = optimset('Display', 'none', 'TolFun', 1e-6);
params_rising = lsqcurvefit(hill_rising, params0, x_data, y_data, [], [], options);
params_decreasing = lsqcurvefit(hill_decreasing, params0, x_data, y_data, [], [], options);

% Plot results
if cutXaxis
    x_data(end) = 8;
end

% Generate fine x-values for plotting
x_fine = linspace(min(x_data)*0.95, max(x_data)*1.05, 100);

% Evaluate fitted models
y_fit_rising = hill_rising(params_rising, x_fine);
y_fit_decreasing = hill_decreasing(params_decreasing, x_fine);

% Compute residuals
residual_rising = norm(y_data - hill_rising(params_rising, x_data));
residual_decreasing = norm(y_data - hill_decreasing(params_decreasing, x_data));

% Determine best fit model
if residual_rising < residual_decreasing
    best_fit = 'Rising Hill Function';
    best_params = params_rising;
    y_best_fit = y_fit_rising;
else
    best_fit = 'Decreasing Hill Function';
    best_params = params_decreasing;
    y_best_fit = y_fit_decreasing;
end

scatter(x_data, y_data, 80, 'kx','linewidth', 3, 'DisplayName', 'Data Points');
% plot(x_fine, y_fit_rising, 'b-', 'LineWidth', 2, 'DisplayName', 'Hill (Rising)');
% plot(x_fine, y_fit_decreasing, 'g--', 'LineWidth', 2, 'DisplayName', 'Hill (Decreasing)');
plot(x_fine, y_best_fit, 'k-', 'LineWidth', 0.5, 'DisplayName', ['Best Fit: ', best_fit]);

if cutXaxis
    % Set custom x-ticks with a break
    xticks([4 5 6 7 8]);  % Include 7 and 11 for a break
    xticklabels({'-4', '-5', '-6', '...', '11'}); % Add '...' as break
end
%%
legend(gca(), 'off');
set(gca,'TickLabelInterpreter','latex')
xlabel('pCa', Interpreter='latex');
% ylabel('Param value');

% title(['Hill Function Fit - Best Model: ', best_fit]);

% Display best-fit parameters
disp(['Best Fit Model: ', best_fit]);
disp(['A = ', num2str(best_params(1))]);
disp(['K = ', num2str(best_params(2))]);
disp(['n = ', num2str(best_params(3))]);

end