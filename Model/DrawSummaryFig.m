clearvars -except saveFigures
f = figure(3);clf;

f.Position = [300.0000  200.0000  691.2000  279.4000];


%
tile1 = subplot(221);
clear params
pCa = 11;
rampSet = [4];
drawPlots = true;
plotDetailedPlots = false;
plotInSeparateFigure = false;
RunCombinedModel;
legend('Location','southeast');
title('pCa = 11', Interpreter='latex');
xlabel('t (s)', 'Interpreter','latex');
ylabel('Stress (kPa)', Interpreter='latex')
%
tile2 = subplot(223);
RunCombinedModel;
legend('off')
ylabel('')
%
tile3 = subplot(222);
clear params
pCa = 4.51;
rampSet = [4];
drawPlots = true;
plotDetailedPlots = false;
plotInSeparateFigure = false;
RunCombinedModel;
xlabel('t (s)', 'Interpreter','latex');
title('pCa = 4.51', Interpreter='latex');
legend('off')
ylabel('Stress (kPa)', Interpreter='latex')
%
tile4 = subplot(224);
RunCombinedModel;
legend('off');
ylabel('')

fontsize(12, 'points')
%%
set(tile1, 'XScale', 'log', 'XLim', [1e-2 160], 'Ylim', [0 20],  'Position', [0.1169    0.150    0.35   0.77]);
set(tile2, 'XScale', 'linear', 'XLim', [0 0.25], 'Ylim', [0 20], 'Position', [0.295   0.55    0.1613    0.35]);

set(tile3, 'XScale', 'log', 'XLim', [1e-2 160], 'Ylim', [0 60], 'Position', [0.5880    0.15    0.35    0.77]);
set(tile4, 'XScale', 'linear', 'XLim', [0 0.25], 'Ylim', [0 60], 'Position', [0.77    0.55    0.1613    0.35]);

