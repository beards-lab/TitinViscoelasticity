% processes the relaxed data
% clear -except dataset;

% % using only the final dataset
% dataset{1} = load('DataStruct20230919.mat');
% dataset{2} = load('DataStruct20230927.mat');
% dataset{3} = load('DataStruct20230928.mat');
% dataset{4} = load('DataStruct20231027.mat');
% dataset{5} = load('DataStruct20231102.mat');
% dataset{6} = load('DataStruct20231107.mat');
% isMale = [1, 1, 0, 0, 1, 0];

% using only the MAVA dataset
addpath ../data
dataset{1} = load('DataStruct20241010.mat');
dataset{2} = load('DataStruct20241121.mat');
dataset{3} = load('DataStruct20241212.mat');
dataset{4} = load('DataStruct20241217.mat');
dataset{5} = load('DataStruct20241219.mat');
dataset{6} = load('DataStruct20241220.mat');
isMale = [1, 0, 0, 0, 1, 1];
leg_names = ["20241010 M","20241121 F","20241212 F","20241217 F","20241219 M","20241220 M"];
%% Get Fmax
AverageRamps_findFmax;

%% Figure representative ramps
if exist('PlotFig2AndDie', 'var') && PlotFig2AndDie
i_dtst = 1;
ts = [dataset{i_dtst}.dsc{1, 1}.datatable.t; ...
    dataset{i_dtst}.dsc{1, 1}.datatable.t(end) + dataset{i_dtst}.dsc{2, 1}.datatable.t ;...
    dataset{i_dtst}.dsc{1, 1}.datatable.t(end) + dataset{i_dtst}.dsc{2, 1}.datatable.t(end)  + dataset{i_dtst}.dsc{3, 1}.datatable.t];
Ls = [dataset{i_dtst}.dsc{1, 1}.datatable.L; ...
    dataset{i_dtst}.dsc{2, 1}.datatable.L ;...
    dataset{i_dtst}.dsc{3, 1}.datatable.L];
Fs = [dataset{i_dtst}.dsc{1, 1}.datatable.F; ...
    dataset{i_dtst}.dsc{2, 1}.datatable.F ;...
    dataset{i_dtst}.dsc{3, 1}.datatable.F];


ts = ts(1:189864); Ls = Ls(1:189864); Fs = Fs(1:189864);
f = gcf();clf;
% tiledlayout('flow', 'TileSpacing','none')
% nexttile;
aspect = 3;
% normal size of 2-col figure on page is 7.2 inches, half is 3.5 inch
% matlab's pixel is 1/96 of an inch
f.Position = [300 200 7.2*96 7.2*96/aspect];
yl1 = [0.75 1.2];
yl2 = [-2.5 15];
gc = axes('Position', [0.1 0.2 0.39 0.75], 'Box','on', 'BoxStyle','full');

% PNB area
to = 2905; % time offset panel A
% time offset panel A - relaxed fast to slow area
to = 906.85;
plot(ts - to, Ls, 'k--'); ylim(yl1);
ylabel('Muscle length ($L/L_0$)', 'Interpreter','latex')
yyaxis right; plot(ts-to, Fs, 'k-', LineWidth=2);
% PNB area
% xlim([0, 795]);ylim([-5 15]);
% fast to slow area
xlim([-50, 1000]);ylim(yl2);


% leg = legend('ML', 'T', 'Interpreter','latex', 'FontSize',12, NumColumns=2, Location='northwest');
% leg.Position = leg.Position + [0 0.07 0 0];
% ylabel('Tension (kPa)', 'Interpreter','latex')
xlabel('$t$ (s)', 'Interpreter','latex');
fontsize(12, 'points');
gc = gca; 
gc.YAxis(2).Color = [0 0 0];
% gc.YAxis(2).Visible = "off";
gc.YAxis(2).TickLabel = [];
gc.XTick = [0 200 400 600 800];


gc = axes('Position', [0.51 0.2 0.39 0.75], 'Box','on', 'BoxStyle','full');

% co = -2.75;
co = 0
plot(ts-to +co, Ls, 'k--');ylim(yl1);
% ylabel('Muscle length ($L/L_0$)', 'Interpreter','latex')
yyaxis right; plot(ts-to + co, Fs, 'k-', LineWidth=2);
% PNB area
% xlim([0, 795]);ylim([-5 15]);
% fast to slow area
xlim([-0.05, 0.5]);ylim(yl2);

leg = legend('ML ($L/L_0$)', '$\Theta$ (kPa)', 'Interpreter','latex', 'FontSize',12, NumColumns=1, Location='southeast');
leg.Position = leg.Position + [0 0.07 0 0];
xlabel('$t$ (s)', 'Interpreter','latex');
ylabel('\textbf{Stress (kPa)}', 'Interpreter','latex')
fontsize(12, 'points');
gc = gca; 
gc.YAxis(2).Color = [0 0 0];
% gc.YAxis(1).Visible = "off";
gc.YAxis(1).TickLabel = [];
% gc.Position = gc.Position + [0.05 0 0 0];
gc.XTick = [0 0.25 0.5];

% exportgraphics(f,'../Figures/RepreRamps.png','Resolution',150)
    return;
end

%% cell for each ramp
normalize = false;
% ramp durations - indexes the _relaxed_
rds = [100, 10, 1, 0.1];

% Only final dataset with the new protocol
relaxed{1} = [1, 1, 9;2, 1, 6;3, 1, 6;4,1,9;5, 1, 9;6, 1, 9];
relaxed{2} = [1, 1, 8;2, 1, 7;3, 1, 7;4,1,8;5, 1, 8;6, 1, 8];
relaxed{3} = [1, 1, 7;2, 1, 8;3, 1, 8;4,1,7;5, 1, 7;6, 1, 7];
relaxed{4} = [1, 1, 6;2, 1, 9;3, 1, 9;4,1,6;5, 1, 6;6, 1, 6];

% peaks = relaxed;

% Only with decay 60s+ - use for identifying the decay
relaxed{1} = [2, 1, 6;3, 1, 6;4, 1, 9;5,1,9;6, 1, 9];
relaxed{2} = [2, 1, 7;3, 1, 7;4, 1, 8;5,1,8;6, 1, 8];
relaxed{3} = [2, 1, 8;3, 1, 8;4, 1, 7;5,1,7;6, 1, 7];
relaxed{4} = [2, 1, 9;3, 1, 9;4, 1, 6;5,1,6;6, 1, 6];

% % with PNB
% relaxed{1} = [1, 1, 2;2, 1, 2;2, 1, 9;3, 1, 6;4, 1, 6;5,1,9];
% relaxed{2} = [1, 1, 3;2, 1, 3;2, 1, 8;3, 1, 7;4, 1, 7;5,1,8];
% relaxed{3} = [1, 1, 4;2, 1, 4;2, 1, 7;3, 1, 8;4, 1, 8;5,1,7];
% relaxed{4} = [1, 1, 5;2, 1, 5;2, 1, 6;3, 1, 9;4, 1, 9;5,1,6];

% PNB and MAVA datasets
% relaxed{i_rds} = [i_dataset, i_condition, i_ramp]
dsName = 'AvgRelaxedMAVASet';
relaxed{4} = [1, 1, 6;2, 1, 6;3, 1, 6;4, 1, 6;5, 1, 6;6, 1, 6;];
relaxed{3} = [1, 1, 7;2, 1, 7;3, 1, 7;4, 1, 7;5, 1, 7;6, 1, 7;];
relaxed{2} = [1, 1, 8;2, 1, 8;3, 1, 8;4, 1, 8;5, 1, 8;6, 1, 8;];
relaxed{1} = [1, 1, 9;2, 1, 9;3, 1, 9;4, 1, 9;5, 1, 9;6, 1, 9;];

dsName = 'AvgRelaxedMAVASetBckwd';
relaxed{1} = [1, 1, 2;2, 1, 2;3, 1, 2;4, 1, 2;5, 1, 2;6, 1, 2;];
relaxed{2} = [1, 1, 3;2, 1, 3;3, 1, 3;4, 1, 3;5, 1, 3;6, 1, 3;];
relaxed{3} = [1, 1, 4;2, 1, 4;3, 1, 4;4, 1, 4;5, 1, 4;6, 1, 4;];
relaxed{4} = [1, 1, 5;2, 1, 5;3, 1, 5;4, 1, 5;5, 1, 5;6, 1, 5;];
 
% dsName = 'AvgRelaxedMAVASet3'; % from the PNB+Mava treated series
% relaxed{1} = [1, 3, 2;2, 3, 2;3, 3, 2;4, 3, 2;5, 3, 2;6, 3, 2;];
% relaxed{2} = [1, 3, 3;2, 3, 3;3, 3, 3;4, 3, 3;5, 3, 3;6, 3, 3;];
% relaxed{3} = [1, 3, 4;2, 3, 4;3, 3, 4;4, 3, 4;5, 3, 4;6, 3, 4;];
% relaxed{4} = [1, 3, 5;2, 3, 5;3, 3, 5;4, 3, 5;5, 3, 5;6, 3, 5;];

% %%
% dtst = dataset{2}.dsc;
% dtst{1, 1}.datasetTitle
% rmp = dtst{1, 1}.datatable;

%
f = figure(2);clf;hold on;
aspect = 1.5;
% normal size of 2-col figure on page is 7.2 inches, half is 3.5 inch
% matlab's pixel is 1/96 of an inch
f.Position = [300 200 3.5*96 3.5*96/aspect];

clear peaks peaks_norm;
if ~normalize
    dsName = [dsName '_NotNorm'];
end
for i_rds = 1:length(rds)
    % sp =subplot(4, 4, (i_rds-1)*4 +  (1:2));cla;
    sp =subplot(4, 1, i_rds);cla;
    % excluding 100s
    % if i_rds ~= 1 % skip plotting of 100s
    %     sp =subplot(3, 4, (i_rds-2)*4 +  (1:2));cla;
    % end
    outF = [];sum_squared_diff = []; n = 1;clear leg;
    rampSet = relaxed{i_rds};
    clin = lines(size(rampSet, 1)+1);

    for i_logtrace = 1:size(rampSet, 1)
        % experiment dataset - whole measurement session
        eds = dataset{rampSet(i_logtrace, 1)}.dsc;
        % dataset - particular conditions
        dtst = eds{rampSet(i_logtrace, 2), rampSet(i_logtrace, 3)};
        leg{i_logtrace} = [dtst.folder ':' dtst.datasetTitle];
        rmp = dtst.datatableZDCorr;
        i_0 = find(rmp.t >= 10, 1);
        % i_end = find(rmp.t >= 10+28+rds(i_rds), 1);
        dt = (rmp.t(end) - rmp.t(1))/(length(rmp.t)-1);
        i_end = find(rmp.L > 1.15, 1, 'last') - 0.5/dt;
        % cut out
        rmp = rmp(i_0:i_end, :);
        rmp.t = rmp.t - 10;
        
        % base rebase not base
        if ~normalize
        base_rel = 1;
        else
        % base on absolute peak
        % base_rel = max(rmp.F);
        % base on steady state
        % base_rel = rmp.F(end);
        % base on max relaxed peak
        % base_rel = peaks_relaxed(4, i_logtrace);
        % base on Fmax, if available (in AverageRampsCa)
        base_rel = dataset_maxF(rampSet(i_logtrace, 1));
        % base on staeady state of the slowest ramp (in AverageRampsCa)
        % base_rel = dataset_ssF(i_logtrace);
        end

        if isnan(base_rel) || base_rel == 0
            continue;
        end

        F = rmp.F/base_rel;
        % F = rmp.F;
        
        % if isMale(i_logtrace)
        if true || i_rds ~= 1 % skip 100s
            semilogx(rmp.t, F, '-', LineWidth=0.5, Color=[clin(i_logtrace, :), 0.1]);hold on;
        end
        % else
        %     semilogx(rmp.t, F, ':', LineWidth=1.5);hold on;
        % end
        
        % base on obsolute peak
        peaks(i_rds, i_logtrace) = max(rmp.F);
        % normalization        
        peaks_norm(i_rds, i_logtrace) = max(F);

        % base on peak over steady state
        % peaks(i_rds, i_logtrace) = max(rmp.F) - base_rel;
        if isempty(outF)
            outF = F;
            Fmax = base_rel;
            n = 1;
            sum_squared_diff = outF*0;
        else
            L = min(length(outF), length(rmp.F));
            n = n + 1;
            % shrink to shortest one
            outF = outF(1:L) + (F(1:L) - outF(1:L))/n;
            outT = rmp.t(1:L);
            % avg the peaks to rescale back
            Fmax = Fmax + (base_rel - Fmax)/n;
             
			% estimate the standard error out of sum squared
            % tested it works the same on-line as using final averaged outF
            sum_squared_diff = sum_squared_diff(1:L) + (F(1:L) - outF(1:L)).^2;
        end    

        %%
        % fitrg = rmp.L > 1.1 & rmp.t < dtst.rd;
        % rmpFitrg = rmp(fitrg, :);
        % % fitfun = @(a, b, c, d, x) min(a*max(x+d, 0).^(b) + c +0*d, 1e2);
        % fitfun = @(a, b, x) a.*(x + b);
        % 
        % [ae goodness] = fit(rmpFitrg.L, rmpFitrg.F,fitfun, 'StartPoint',[1, 1]);
        % as(i_rds, i_logtrace) = ae.a;

        % figure(7);
        % semilogx(rmp.t + rampShift(i_rds), rmp.F, '-', Color=[clin(i_logtrace, :), 0.1]);hold on;
        % figure(3);        
    end

%% resample and save
    
    if outT(end) < 200
        % resample up the peak and for the tail separately, up to 60s
        % t_s = [linspace(0, 1, 20)*rds(i_rds), ...                
        %        rds(i_rds) + linspace(0, min(60, outT(end) - rds(i_rds)), 40*4)];
        % t_s = [linspace(0, 1, 20)*rds(i_rds), ...                
        %        rds(i_rds) + logspace(log10(1e-3), log10(min(60, outT(end) - rds(i_rds))), 40)];
        t_s = [linspace(0, 1, 20)*rds(i_rds), ...                
               logspace(log10(rds(i_rds)), log10(rds(i_rds) + min(60, outT(end) - rds(i_rds))), 80)];
        

    else
        % resample up the peak and for the tail separately, up to 300s
        % t_s = [linspace(0, 1, 20)*rds(i_rds), ...                
        %        rds(i_rds) + linspace(0, min(300, outT(end) - rds(i_rds)), 40*4*5)];
        % t_s = [linspace(0, 1, 20)*rds(i_rds), ...                
        %        rds(i_rds) + logspace(log10(1e-3), log10(min(300, outT(end) - rds(i_rds))), 50)];
        t_s = [linspace(0, 1, 20)*rds(i_rds), ...                
               logspace(log10(rds(i_rds)), log10(rds(i_rds) + min(300, outT(end) - rds(i_rds))), 100)];

        % I do not know how to get 50 samples with the extension, it just
        % overlaps in the plot somehow
        % clf;
        % plot(1:length(t_s), t_s, 1:length(t_s2), t_s2)
    end
    %% resample log equally
    % t_s = [logspace(log10(1e-3), log10(outT(end)), 40)];
    % remove consequent duplicates at joints
    t_s = t_s(~[false t_s(2:end) == t_s(1:end-1)]);
    SD = sqrt(sum_squared_diff / (n * (n - 1)));
    % Force, Length and SD interpolation (relative to Fmax)
    FLSDint = interp1(outT, [outF, rmp.L(1:L), SD], t_s, "pchip", 'extrap');
    
    tab_rmpAvg = table(t_s' + 2, FLSDint(:, 2), FLSDint(:, 1)*Fmax, FLSDint(:, 3)*Fmax);
    tab_rmpAvg.Properties.VariableNames = {'Time', 'L', 'F', 'SD'};
    writetable(tab_rmpAvg, ['../data/' dsName '_' num2str(rds(i_rds)) 's.csv']);
%% plot the AVG
    % semilogx(outT, outF, 'k-');hold on;
    % semilogx(t_s, FLSDint(:, 1) + FLSDint(:, 3), '--', Color=[clin(end, :), 1], LineWidth=2);
    % semilogx(t_s, FLSDint(:, 1) - FLSDint(:, 3), '--', Color=[clin(end, :), 1], LineWidth=2);
    % plot(t_s, FLSDint(:, 1), '-|', Color=[clin(end, :), 0.5], LineWidth=2)
    
    % Fill the area between upper and lower bounds to show the confidence interval
    % fill([outT', fliplr(outT')], [(outF - SE*1.96)', fliplr((outF + SE*1.96)')], 'b');

    errorbar(t_s,FLSDint(:, 1),FLSDint(:, 3), '-', LineWidth=2, Color=clin(end, :))
    fontsize(14, 'points');
    % pos = sp.Position;
    % sp.Position = [0.1 pos(2),0.35, pos(4)];
    xlim([1e-2 160])    
    yl = ylim;
    ylabel('$T_{rel}$', 'Interpreter','latex');
    yyaxis right; ylim(yl*Fmax);
    ylabel('T (kPa)', 'Interpreter','latex');
    % plot(rmp.t(1:L), outF(1:L) + (rmp.F(1:L) - outF(1:L))/n)
    leg{length(leg) +1} = 'Averaged';
    title(['Ramp ' num2str(rds(i_rds)) 's'])

    % if i_rds == 1
        % x = 0.1;y = 0.03;
        % legend(leg, 'Interpreter','none', 'Position', ...
        %     [sp.Position(1) + sp.Position(3) + x, sp.Position(2) - y, 1 - sp.Position(1) - sp.Position(3) - x, sp.Position(4)])
    xlim([1e-2 300]);
    if true || i_rds == 4
        xlabel('t (s)')
        xticks([0.1 10])
    else
        xticks([]);
    end   
    
    Tarr{i_rds} = outT;
    Farr{i_rds} = outF*Fmax;

end
save(['pca11data' dsName '.mat'], "Farr", "Tarr", "peaks");
peaks11 = peaks;
%%
aspect = 2;
% f = figure(4);f.Position = [300 200 7.2*96 7.2*96/aspect];clf;

% gc = axes('Position', [0.1 0.1 0.39 0.8], 'Box','on', 'BoxStyle','full');
% gc = axes('Position', [0.6 0.1 0.39 0.8], 'Box','on', 'BoxStyle','full', 'Color','r');

% subplot(4, 4, [3 16]);cla;
% remove zero peaks as NaNs to fix the average - only when something was missing
peaks(peaks == 0) = NaN;


% for x axis we go from fastest to slowest
x_ax = length(rds):-1:1;


% colors
clin = lines(size(peaks, 2)+1);
% BW only
clin = repmat([0 0 0], [size(peaks, 2)+1, 1]);
% blue only
clin = repmat([0    0.4470    0.7410], [size(peaks, 2)+1, 1]);
% Red only
% clin = repmat([0.8500    0.3250    0.0980], [size(peaks, 2)+1, 1]);

boxplot(fliplr(peaks'), PlotStyle="traditional", Notch="off", Colors=clin(end, :));hold on;

clear lp;
for i_pk = 1:size(peaks, 2)
    if all(isnan(peaks(:, i_pk))) || all (peaks(:, i_pk) == 0)
        % it was just a placeholder
        leg{i_pk} = '';
        continue;
    end    
    % set(gca, 'colororderindex', 1);
    lp = plot(x_ax', peaks(:, i_pk), '.--', 'MarkerSize',12, LineWidth=1.5, Color=clin(i_pk, :));
    hold on;    
    % semilogx(rds, as(:, i_pk), 'x:', 'MarkerSize',12, LineWidth=0.5, Color=clin(i_pk, :));hold on;
end
% set(gca, 'YAxisLocation', 'right');
% just for the legend
lp(length(lp)+1) = plot(NaN, NaN, 's--',LineWidth=2, Color=clin(end, :), MarkerSize=5);
% validLeg = ~cellfun(@isempty,leg);
% %# remove empty cells
% leg_cleared = leg(validLeg);

% leg_cleared = {'20230919 M','20230927 M','20230928 F','20231027 F','20231102 M','20231107 F', 'Averaged'};
% 5 longer datasets
% leg_cleared = {'20230927 M','20230928 F','20231027 F','20231102 M','20231107 F', 'Averaged'};
% summed up
leg_cleared = {'Individual data', 'Averaged data'};
legend(lp(end-1:end), leg_cleared, 'Interpreter','none', 'AutoUpdate','off', 'Location','best')

plot(x_ax, nanmean(peaks, 2)', '-', LineWidth=2, Color=clin(end, :))
plot(x_ax, nanmean(peaks, 2)', 's', LineWidth=3, Color=clin(end, :), MarkerSize=5)


xlabel('$t_r$ (s)', Interpreter='latex');
% ylabel('$T_{rel}$ (kPa)', Interpreter='latex')
% yl = ylim();yyaxis right;ylim(yl*Fmax);
ylabel('$\Theta$ (kPa)', Interpreter='latex');
set(gca, 'XTickLabel', {'0.1', '1', '10', '100'});
% g = gca();g.YAxis(2).Color = [0 0 0];
ylim([0 inf])
set(gca, 'FontSize', 12);
title('Peak stress - relaxed')
% exportgraphics(f,'Figures/AvgpPeaksRelaxed.png','Resolution',150)
%%

% Data (example values - replace with your own)
%% peak sum up - run BEFORE AverageRampsCa
d = load("pca11dataAvgRelaxedMAVASet.mat");
% d = load("pca11dataAvgRelaxedMAVASet_NotNorm.mat");
% d = load("pca11dataAvgRelaxedMAVASet3.mat");
% d = load("pca11dataAvgRelaxedMAVASet3_NotNorm.mat");
% d = load("pca11dataAvgRelaxedMAVASetBckwd.mat" );

peaks_1  = flipud(d.peaks);

d = load("pca11dataAvgRelaxedMAVASet3.mat");
% d = load("pca11dataAvgRelaxedMAVASet3_NotNorm.mat");
% d = load("pca11dataAvgRelaxedMAVASetBckwd.mat" );
% d = load("pca11dataAvgRelaxedMAVASetBckwd_NotNorm.mat" );

legtext = {'Relaxed, preconditioned', 'Relaxed, PNB+Mava'};
% legtext = {'Relaxed, not preconditioned', 'Relaxed, PNB+Mava', };

peaks_2 = flipud(d.peaks);
categories = {'0.1', '1', '10', '100'};



% Initialize figure
figure(333);clf;
hold on;

% Colors
colors = lines(2);  % 2 distinct colors


% Prepare for individual line plots across categories
x_vals = 1:length(categories);
markers = ["o", "s", "d", "^", "v", "h"];

% Plot individual lines for Data1 and Data2 across categories
for subj = 1:nSubjects
    y_data1 = peaks_1(:, subj);  % 1 value per category
    y_data2 = peaks_2(:, subj);
    
    plot(x_vals - 0.15, y_data1, markers(subj), 'Color', [0.5 0.5 1], 'MarkerFaceColor', [0.5 0.5 1]);
    plot(x_vals + 0.15, y_data2, markers(subj), 'Color', [1 0.5 0.5], 'MarkerFaceColor', [1 0.5 0.5]);
end


% Plot boxcharts
for i = 1:numel(categories)
    % Get data
    y1 = peaks_1(i, :);
    y2 = peaks_2(i, :);
    
    % Box positions (group i with offsets for D1 and D2)
    x1 = repmat(i - 0.15, 1, length(y1));  % shift left
    x2 = repmat(i + 0.15, 1, length(y2));  % shift right
    
    % Plot boxcharts
    b1 = boxchart(x1, y1, 'BoxFaceColor', colors(1,:), 'BoxWidth', 0.25);
    b2 = boxchart(x2, y2, 'BoxFaceColor', colors(2,:), 'BoxWidth', 0.25);
    
    % Add only first to legend
    if i == 1
        set(get(get(b1,'Annotation'),'LegendInformation'), 'IconDisplayStyle', 'on');
        set(get(get(b2,'Annotation'),'LegendInformation'), 'IconDisplayStyle', 'on');
    else
        set(get(get(b1,'Annotation'),'LegendInformation'), 'IconDisplayStyle', 'off');
        set(get(get(b2,'Annotation'),'LegendInformation'), 'IconDisplayStyle', 'off');
    end
        % Statistical test (paired)
    [~, p] = ttest(y1, y2);
    
    % Determine annotation
    if p < 0.001
        stars = '***';
    elseif p < 0.01
        stars = '**';
    elseif p < 0.05
        stars = '*';
    else
        stars = '';
    end
    % Plot asterisk annotation above the boxes
    maxY = max([y1, y2]);
    
    % text(i, maxY + 0.8, sprintf('%.4f %s', p, stars), 'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
end

% Axes settings
xlim([0.5, numel(categories) + 0.5]);
xticks(1:numel(categories));
xticklabels(categories);
xlabel('Ramp duration (s)');
ylabel('Peak Value (kPa)');
title('Grouped Boxplots comparison');
legend(legtext, 'Location', 'Best');
grid on;
