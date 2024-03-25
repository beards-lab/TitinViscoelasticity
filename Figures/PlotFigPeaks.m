% Plots comparison of peak heights for all ramps and for intermediate
% calcium levels

pCas = [4.4, 5.5, 5.75, 6, 11];

pCaPeaksModel = nan(length(pCas), 4);
pCaPeaksData = nan(length(pCas), 4);

for i_pCa = 1:length(pCas)
    pCa = pCas(i_pCa);
    model_data = load(sprintf('..\\pca%gmodeldata.mat', pCa));
    
    for j_rds = 1:length(model_data.Farr)
        if isempty(model_data.Farr{j_rds})
            continue;
        end
        
        pCaPeaksModel(i_pCa, j_rds) = max(model_data.Farr{j_rds});
        if isinf(pCa) || pCa >= 10
          %   % hack - the no-Ca noPNB experiments had higher ramps
          %   datatable = readtable(['..\Data\bakers_passiveStretch_' num2str(rds(i_rd)*1000) 'ms.csv']);
          %   datatable.Properties.VariableNames = {'Time'  'ML'  'F'  'SL'};
          %   datatables{i_rd} = datatable;
          % elseif isnan(pCa)
              % newest format of experiments    
            datatables{j_rds} = readtable(['..\Data\AvgRelaxed_' num2str(rds(j_rds)) 's.csv']);
        else
            datatables{j_rds} = readtable(['..\Data\AvgpCa' num2str(pCa) '_' num2str(rds(j_rds)) 's.csv']);
        end
        [m, i_m] = max(datatables{j_rds}.F);
        pCaPeaksData(i_pCa, j_rds) = m;
        pCaPeaksData_SD(i_pCa, j_rds) = datatables{j_rds}.SD(i_m);
    end
end
%%
f = figure(232);clf;
aspect = 2;
f.Position = [300 200 7.2*96 7.2*96/aspect];
% tiledlayout('flow', TileSpacing='compact');
lw = 1;
ms = 6;
% tl1 = nexttile;
tl1 = axes('Position', [0.1300    0.1790    0.3400    0.7460]);
fontsize(12, 'points');hold on;

% for i = [1 5]
i = 1;
pld = errorbar(rds, pCaPeaksData(i, :), pCaPeaksData_SD(i, :), 'ok', LineWidth=lw);
pl1 = semilogx(rds, pCaPeaksModel(i, :), 'ks--', 'MarkerFaceColor','k', LineWidth=lw, MarkerSize = ms)
i = 5;
errorbar(rds, pCaPeaksData(i, :), pCaPeaksData_SD(i, :), 'ok', LineWidth=lw);
pl2 = semilogx(rds, pCaPeaksModel(i, :), 'ks--', 'MarkerFaceColor','k', LineWidth=lw, MarkerSize = ms)
% only for the elgend
plm = plot(NaN, NaN, 'ks', 'MarkerFaceColor','k', LineWidth=lw, MarkerSize = ms);

% end
tl1.XScale = 'log';
set(gca, 'Xtick',fliplr(rds));
ylabel('$\Theta$ (kPa)', Interpreter='latex');
ylim([0, 50]);
xlabel('$t_r$ (s)', Interpreter='latex');
% legend([pld, pl1, pl2], 'Data', 'pCa 4.4 (model)', 'Relaxed (model)')
% legend([pl1, pl2], 'pCa 4.4 (model)', 'Relaxed (model)')
legend([pld, plm], 'data', 'model')
text(3, 32, 'pCa 4.51', Interpreter='latex');
text(0.2, 18, 'pCa 11', Interpreter='latex');


% tl2 = nexttile;fontsize(12, 'points');
w = 0.08;s = 0.015;
tl2 = axes('Position', [0.5650    0.1790    w    0.7460]);

fontsize(12, 'points');hold on;
pld = errorbar(-pCasP, pCaPeaksData(:, i_rds),pCaPeaksData_SD(:, i_rds), 'ok', LineWidth=lw);
plot(-pCasP, pCaPeaksModel(:, i_rds), 'ks--', 'MarkerFaceColor','k',LineWidth=2, MarkerSize = ms)
ylim([0, 50]);
d = 1;
xlim([-11 - d -11 + d])
set(gca, 'XTick',-[11]);
set(gca, 'Xticklabels', ["11"]);
box off;
ylabel('$\Theta$ (kPa)', Interpreter='latex');

tl3 = axes('Position', [0.5650 + w + s    0.1790    0.3400 - w - s    0.7460]);
hold on;
ylim([0, 50]);
xlabel('pCa', Interpreter='latex')
fontsize(12, 'points')
xlim([-6.5 -4])
pCasP = pCas;
pCasP(1) = 4.51;
% pCasP(end) = 8;

% for i_rds = [4]
    pld = errorbar(-pCasP, pCaPeaksData(:, i_rds),pCaPeaksData_SD(:, i_rds), 'ok', LineWidth=lw);
    plot(-pCasP, pCaPeaksModel(:, i_rds), 'ks--', 'MarkerFaceColor','k',LineWidth=2, MarkerSize = ms)
    % pl1 = plot(-pCasP(end), pCaPeaksModel(end), 'ks:', 'MarkerFaceColor','k',LineWidth=2, MarkerSize = ms);
% end
% legend([pld pl1], 'Data', 'Ramp 0.1s (model)', 'Location', 'southwest')
text(-5.7, 22, '$t_r$ = 0.1 s', Interpreter='latex');
set(gca, 'XTick',-[7:-1:4]);
set(gca, 'Xticklabels', ["" "6" "5" "4"]);
set(gca, 'YTick', []);
% xlim([-8.5 -4])
fontsize(12, 'points')
tl1.Position
tl2.Position
% exportgraphics(f,'../Figures/FigModelPeaks.png','Resolution',150)
% exportgraphics(f,'../Figures/FigModelPeaks.eps')

