% Run slack experiment


simtype = 'velocitytable_slack';

loaddata = load('../data/bakers_slack8mM_all.mat');

velocitytable = loaddata.velocitytable;
velocitytable(:, 3) = velocitytable(:, 4)/2;
velocitytable = velocitytable(:, 1:3);

vtb = array2table(velocitytable );    
vtb.Properties.VariableNames = {'Time', 'Velocity','ML'};
writetable(vtb, ['../data/' simtype '.csv']);
% vtb

drawPlots = false;
rampSet = [1];
alphaF_0 = 100;

%% first round
pCa = 11;
clear params;
RunCombinedModel;
pca12Result.t = Time{1};
pca12Result.F = Force{1};

%% second round
pCa = 4.51;
clear params;
RunCombinedModel;
pca4Result.t = Time{1};
pca4Result.F = Force{1};    


%% postprocess and visualize
tx = 1000;
t0 = datatables{1}.Time(1)*tx;


F_data{1} = interp1(datatables{1}.Time, datatables{1}.F, Time{1}); % total force interpolated

figure(1);clf;
L_data{1} = interp1(datatables{1}.Time, datatables{1}.L, Time{1}); % total force interpolated
plot(Time{j}*tx -t0, Length{j} + 0.95, Time{j}*tx -t0, L_data{1}/2); 

figure(2);clf;
plot(Time{1}*tx -t0, F_data{1}, '-',pca12Result.t*tx -t0, pca12Result.F, '-', pca4Result.t*tx -t0, pca4Result.F, '-',lineWidth = 2);hold on;

legend('Slack protocol Force', 'Titin VE - pCa 11', 'Titin VE - pCa 4.5')
xlabel('time (ms)')
ylabel('Force (kPa)');