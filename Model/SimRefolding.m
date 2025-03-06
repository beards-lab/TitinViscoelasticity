% Test reforlding

clearvars -except saveFigures
cf = figure(1);clf;
drawPlots = true;

% one second period
rampSet = [3];
    
plotInSeparateFigure = false;

% pCa 4.5
params = [5.19       12.8       4345       2.37      4e+04       2.74  8.658e+05      5.797      0.678      0.165   0.005381      0.383];
pCa = 4.51;
% pCa 11
params = [5.19       12.8      512.3       2.37      4e+04       2.74  2.668e+07      9.035      0.678      0.165        NaN        NaN];
pCa = 11;

% params(9) = 0.07;
simtype = 'ramp';

alphaF_0 = .1;
% params(1) = 2.5;
RunCombinedModel;

figure(2);clf;
% params(9) = 0.07;
simtype = 'sin';
RunCombinedModel;