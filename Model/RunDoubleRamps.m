% RunDoubleRamps
% simtype = 'velocitytable_doubleramp_relaxed.csv';
simtype = 'velocitytable_doubleramp2_relaxed.csv';
pCa = 11;

% simtype = 'velocitytable_doubleramp_active.csv';
simtype = 'velocitytable_doubleramp2_active.csv';
pCa = 4.51;

% simtype = 'ramp';
clear params;
alphaF_0 = 0;

rampSet = [1 2 3 4];
rampSet = [1];
drawAllStates = 1;
% params(9 ) =  0.6780*0.1;
compareFig = false;
RunCombinedModel;
%%
figure(161); nexttile(5);set(gca, 'XScale','linear', 'XLim', [0 inf])
%%
params(9 ) =  0.6780*1;
% compareFig = false;
alphaF_0 = 10;
RunCombinedModel;
figure(161); nexttile(5);set(gca, 'XScale','linear', 'XLim', [0 inf]);
figure(161); nexttile(9);set(gca, 'XScale','linear', 'XLim', [0 inf]);
yyaxis right;
plot(Time{1}, (Length{1}+0.95));
% %%
% figure(202);clf;
% nexttile;
% plot(datatables{1}.Time, datatables{1}.L-0.95,Time{1}, Length{1});
% nexttile;
% plot(Time{1}, Force{1}, datatables{1}.Time, datatables{1}.F, LineWidth=2);

%% remove the trendline

data_force = datatables{1}.F;
data_length = datatables{1}.L-0.95;
data_time = datatables{1}.Time;
i_zeropoints = data_length < -0.1;
data_force0 = data_force(i_zeropoints);
data_force0Time = datatables{1}.Time(i_zeropoints);
data_force_corr = data_force - sf(data_time);
sim_time = Time{1};
sim_force = Force{1};

figure(296);clf;

% show the cleared trend in the data
% nexttile;hold on;
% plot(data_time, data_force);
% plot(data_force0Time, data_force0);
% sf = fit(data_force0Time, data_force0, 'poly1');
% plot(data_time, sf(data_time), LineWidth=2);


nexttile;hold on;

plot(data_time, data_force_corr);
plot(sim_time, sim_force, LineWidth=2);
[data_forcepeaks, i_data_forcepeaks] = findpeaks(data_force_corr,data_time,'MinPeakWidth',1.5e-3,'MaxPeakWidth',1, ...
              'MinPeakProminence',5, 'Annotate','extents','MinPeakDistance',25);
sim_forcepeaks = findpeaks(sim_force,sim_time,'MinPeakWidth',.5e-3,'MaxPeakWidth',100, ...
              'MinPeakProminence',2,'MinPeakDistance',25, 'Annotate','extents');

% just for the annotation
findpeaks(sim_force,sim_time,'MinPeakWidth',.5e-3,'MaxPeakWidth',100, ...
              'MinPeakProminence',2,'MinPeakDistance',25, 'Annotate','extents');

% data_forcepeaks = data_forcepeaks/max(data_forcepeaks);
% sim_forcepeaks = sim_forcepeaks/max(sim_forcepeaks);

first_peaks = data_forcepeaks(1:2:end);
second_peaks = data_forcepeaks(2:2:end);
sim_first_peaks = sim_forcepeaks(1:2:end);
sim_second_peaks = sim_forcepeaks(2:2:end);

slack_durs = [0 5 10 20 50 100 200 500]*1e-3 + 10e-3;
% sim_slack_durs = [0 5 10 20]*1e-3 + 10e-3;
sim_slack_durs = slack_durs;

nexttile();
semilogx(slack_durs, first_peaks, 'ks-',slack_durs,second_peaks, 'kv--', LineWidth=2); 
hold on;
semilogx(sim_slack_durs, sim_first_peaks, 'rx-',sim_slack_durs,sim_second_peaks, 'r+--', LineWidth=2); 
legend('First peak', 'Second peak','SIM: First peak', 'SIM: Second peak');

xlabel('Refoldind duration (s)')



%%

    datatables{1} = readtable('..\Data\2025 11 21 Export/05_Log_Active_Refolding.txt');    
    datatables{1}.Properties.VariableNames = {'Time', 'L','F', 'SL'};
    datatables{1}.Time = datatables{1}.Time /1000;% convert to ms
