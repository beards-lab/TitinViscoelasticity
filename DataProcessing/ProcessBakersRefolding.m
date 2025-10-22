%% Titin refolding loading and processing

S1 = dir('../data/2025 09 19 Export');

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
for i = 1:length(S)
    fprintf('Processing %d..\n', i)
    ds = struct();
    ds.filename = S(i).name;
    folders = split(S(i).folder, '\');
    ds.folder = folders{end};


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
 

for i = 6:length(S)
% i = 9;
 
    clf;
    
    nexttile;plot(AllFiles{i}.t, AllFiles{i}.F);
    title(sprintf('%s', AllFiles{i}.filename));
    xlabel('$t$ (s)', Interpreter='latex');    ylabel('F (kPa)', Interpreter='latex');
    
    if ~isempty(AllFiles{i}.SL)
        nexttile;plot(AllFiles{i}.t, 2*AllFiles{i}.L, AllFiles{i}.t, AllFiles{i}.SL);
        nexttile;plot(AllFiles{i}.L, AllFiles{i}.SL);
    else
        nexttile;plot(AllFiles{i}.t/1000, AllFiles{i}.L);
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

    nexttile;
    plot(t, y, t, [0 diff(x)./diff(t)*1000]);
    nexttile;plot(y, [0 diff(x)./diff(t)*1000]);
    % nexttile;plot(t, y, t, x*10);
end

%%
drawPlots = true;

% one second period 
% rampSet = [3]; 
    
plotInSeparateFigure = true;
% pCa 4.5
% params = [5.19       12.8       4345       2.37      4e+04       2.74  8.658e+05      5.797      0.678      0.165   0.005381      0.383];
% pCa = 4.51;
% pCa 11
params = [5.19       12.8      512.3       2.37      4e+04       2.74  2.668e+07      9.035      0.678      0.165        NaN        NaN, 1];
params(9) = 0.1;
alphaF_0 = 1;
pCa = 11;
rampSet = 1;

% params(9) = 0.07;
simtype = 'sin0_25';
% params(9) = 5e-1;

RunCombinedModel;
i = 23;
nexttile(1);hold on;plot(AllFiles{i}.t/1000, AllFiles{i}.L);
nexttile(2);hold on;plot(AllFiles{i}.L, AllFiles{i}.F);
nexttile(3);hold on;plot(AllFiles{i}.t/1000, AllFiles{i}.F);