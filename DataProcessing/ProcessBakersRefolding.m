%% Titin refolding loading and processing

S1 = dir('../data/2025 09 19 Export');
S1 = dir('../data/2025 11 21 Export');

S1 = S1(~[S1.isdir]);
[~,idx] = sort({S1.name});
S1 = S1(idx);

% mergedTables = [struct2table(S1);struct2table(S2)];
% S = table2struct(mergedTables);
S = S1;
%%
% figure(101);clf;hold on;
skipPlots = false;

dsc = cell(0); % dataset structure cell array
clear AllFiles;
processList = 1:length(S);
processList = [7, 8, 9, 10, 11, 12, 13, 14];
for i = processList
    fprintf('Processing %d:%s..\n', i, S(i).name);
    ds = struct();
    ds.filename = S(i).name;
    folders = split(S(i).folder, '\');
    ds.folder = folders{end};

% continue;
    datatable = readtable([S(i).folder '/' S(i).name], 'filetype', 'text', 'NumHeaderLines',4);
    if length(datatable.Properties.VariableNames) == 3
        datatable.Properties.VariableNames = {'t', 'L','F'};
    elseif length(datatable.Properties.VariableNames) == 4
        datatable.Properties.VariableNames = {'t', 'L','F', 'SL'};
    else
        disp('Wat?')
    end

    ds.t = datatable.t;
    ds.L = datatable.L;
    ds.F = datatable.F;

    if length(datatable.Properties.VariableNames) == 4
        ds.SL = datatable.SL;
    else
        ds.SL = [];
    end
AllFiles{i} = ds;
end

%%
for i = 1:length(S)
    fprintf('%d: %s\n', i, S(i).name);
end
%%
close all;
isel = [1 2 3 4 5];
isel = [7 8 9 10 11]
% isel = [32, 34, 31, 33];

% for i = 1:length(S)
for i = isel

% i = 9;
  
    figure(100+i);
    clf;
    
    nexttile;plot(AllFiles{i}.t, AllFiles{i}.F);
    title(sprintf('%s', AllFiles{i}.filename));
    xlabel('$t$ (ms)', Interpreter='latex');    ylabel('F (kPa)', Interpreter='latex');
    
    if ~isempty(AllFiles{i}.SL)
        nexttile;plot(AllFiles{i}.t, 2*AllFiles{i}.L, AllFiles{i}.t, AllFiles{i}.SL);
        xlabel('$t$ (ms)', Interpreter='latex');    ylabel('2*ML (-), SL (um)', Interpreter='latex');
        % nexttile;plot(AllFiles{i}.L, AllFiles{i}.SL);
    else
        nexttile;plot(AllFiles{i}.t, 2*AllFiles{i}.L);
        xlabel('$t$ (ms)', Interpreter='latex');    ylabel('2*ML (-)', Interpreter='latex');
        
    end

 
    % rng = Time{j} < 20;
    t = AllFiles{i}.t';
    x = AllFiles{i}.L';
    y = AllFiles{i}.F'; 
    % y2 = Force{j} - Force_par{j};
    z = zeros(size(t));
    col = linspace(0, t(end), length(t)); 
    lw = 1;

    nexttile;
    surface([x;x],[y;y],[z;z],[col;col],...
            'facecol','no',...
            'edgecol','interp',...
            'linew',lw);
    % clim([0 50 100]) 
    colormap(turbo);    
    ylabel('$F$ (kPa)', Interpreter='latex');     xlabel('$L$ ($\mu$m)', Interpreter='latex');
    cb = colorbar;  title(cb, 't (s)');
    disp('')

    nexttile;plot(t, y, t, [0 diff(x)./diff(t)*1000]);
    % nexttile;plot(y, [0 diff(x)./diff(t)*1000]);
    % nexttile;plot(t, y, t, x*10);
end

%%
clf;
drawPlots = true;

% one second period 
% rampSet = [3]; 
    
plotInSeparateFigure = true;
% pCa 4.5
params = [5.19       12.8       4345       2.37      4e+04       2.74  8.658e+05      5.797      0.678      0.165   0.005381      0.383];
pCa = 4.51;
% pCa 11
% params = [5.19       12.8      512.3       2.37      4e+04       2.74  2.668e+07      9.035      0.678      0.165        NaN        NaN, 1];
params(9) = 0.1;
% alphaF_0 = 1;
% pCa = 11;
rampSet = 1;

% params(9) = 0.07;
simtype = 'sin0_25';
simtype = 'sin2_5';
% params(9) = 5e-1;

RunCombinedModel;
%% Are the protocols the same?
figure(202);clf;
% plot(AllFiles{9}.t, AllFiles{9}.L,AllFiles{11}.t, AllFiles{11}.L);
n9 = length(AllFiles{9}.L);
n11 = length(AllFiles{11}.L);
plot(1:n9, AllFiles{9}.L,1:n11, AllFiles{11}.L);
plot(AllFiles{11}.t(1:end-1)/1000, diff(AllFiles{11}.L)./diff(AllFiles{11}.t/1000))
L = AllFiles{11}.L;
t = AllFiles{11}.t/1000;


% baseline to pos1
v1 = 100;
% pos1 to pos2
v2 = -150;
% pos2 to pos3
v3 = 100;
% pos3 to pos4
v4 = -40;
% pos4 to pos5
v5 = 0.025;

baseline = 0.95;
pos1 = 1.175;
pos2 = baseline;
pos3 = pos1;
pos4 = 0.8;
pos5 = baseline;
% get indexes of baseline to pos1 and pos3 based on the velocity dX/dt +
% conversion from ms
i_p1p3 = find(diff(diff(L)./diff(t) > 30));
i_p1 = i_p1p3(1:2:end);
i_p3 = i_p1p3(2:2:end);
i_p2p4 = find(diff(diff(L)./diff(t) < -30));
i_p2 = i_p2p4(1:2:end);
i_p4 = i_p2p4(2:2:end);
i_p5 = i_p2p4 + 1e4;

hold on; plot(i_p1p3, L(i_p1p3), '*')

% create velocitytable - time, velocity
i = 1;
vt = [t(i_p1(i)), v1
    t(i_p1(i)) + (pos1-pos2)/v1, 0
    t(i_p2(i)), v2
    t(i_p2(i)) + (pos2-pos3)/v2, 0
    ]
%%
i = 9;
figure(i)
nexttile(1);hold on;plot(AllFiles{i}.t/1000, AllFiles{i}.L);
nexttile(2);hold on;plot(AllFiles{i}.L, AllFiles{i}.F);
nexttile(3);hold on;plot(AllFiles{i}.t/1000, AllFiles{i}.F);