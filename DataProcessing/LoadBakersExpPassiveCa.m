%% read new Anthony's dataset
% expected format:
% SS_NN[_RD][_TYPE][_holdXXXs][_CONDS].txt
% SS - zero padded dataset order
% NN - zero padded order, e.g. 01. Log has 01, if present
% RD - ramp duration, ending with 's', e.g. 100s or 0.1s or 'Log' for log
% TYPE: 
%  'Relax', 
%  'RelaxS2F' - relax slow to fast, 
%  'RelaxF2S' - fast to slow, 
%  'RelaxHold' - ramp hold 300s
%  'pCa4.4', 
%  'pCa5.5' etc.
% CONDS - optionally any other conditions or control sequence. 
%  'PNB' - it was in a PNB solution
%  'NoZDR' - not suitable for zero-drift removal, i.e. interpolation
%  between two slack regions. Zero-value removal is preferred

% holdXXXs : optional, specifying different ramp hold time, e.g. hold300s
close all;
clear;clc;

% rds = [0.1 1 10 100];
% rds = fliplr([0.1 1 10 100]);

% ramp height 0.95 - 1.175, i.e. 1.9 - 2.35um
% S1 = dir('../data/PassiveCaSrc2/20230518_renamed');
% S1 = dir('../data/PassiveCaSrc2/20230919');
% S1 = dir('../data/PassiveCaSrc2/20230927');
% S1 = dir('../data/PassiveCaSrc2/20230928');
% S1 = dir('../data/PassiveCaSrc2/20231027');
% S1 = dir('../data/PassiveCaSrc2/20231102');
% S1 = dir('../data/PassiveCaSrc2/20231107');
% S1 = dir('../data/PassiveCaSrc2/20240705');
% S1 = dir('../data/PassiveCaSrc2/20241010');
% S1 = dir('../data/PassiveCaSrc2/20241121');
% S1 = dir('../data/PassiveCaSrc2/20241212');
S1 = dir('../data/PassiveCaSrc2/20241217');
S1 = dir('../data/PassiveCaSrc2/20241219');
S1 = dir('../data/PassiveCaSrc2/20241220');

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
for i = 1:length(S)
    ds = struct();
    ds.filename = S(i).name;
    folders = split(S(i).folder, '\');
    ds.folder = folders{end};


    % categorized based on name
    np = split(S(i).name, {'_', '.txt'}); % name part
    if length(np) < 3
        % unsupported format, we want at least three elements
        continue;
    end
    ds.set = str2double(np{1}); % dataset set. A set is one series of ramps
    ds.order = str2double(np{2}); % dataset order. One piece of ramp.
    if isnan(ds.set) || isnan(ds.order)
        % unsupported format
        continue;
    end
    
    ds.rd = str2double(np{3}(1:end-1)); % ramp duration
    ds.isLog = strcmp(np{3}, 'Log');
    if ds.isLog
        ds.order = 1;
    end
    ds.type = np{4};
    ds.pCa = str2double(np{4}(4:end));
    % optional format: hold100s, otherwise normal hold time
    % holdTime = str2double(np{5}(5:end-1));
    if length(np) > 5 && length(np{5}) > 5 && ...
        contains(np{5}, 'hold') && ~isnan(str2double(np{5}(5:end-1)))
        ds.holdTime = str2double(np{5}(5:end-1));
    elseif ds.isLog
        ds.holdTime = NaN;
    else
        % revert to normal
        ds.holdTime = 60;
    end

    % zero drift removal
    if any(contains(np, 'NoZDR'))
        ds.ZDR = false;
    else
        ds.ZDR = true;
    end

    % cut out the title    
    ds.datasetTitle = sprintf('%d:%s', ds.set, ds.type);
    
    if ds.isLog
        ds.datasetLegend = sprintf('%s log', ds.type);
    else
        ds.datasetLegend = sprintf('%s:%0.1fs', ds.type, ds.rd);
    end

    % log
    fprintf('Loading %s (%s)...', ds.datasetLegend, ds.filename);

    % read the data
    datatable = readtable([S(i).folder '/' S(i).name], 'filetype', 'text', 'NumHeaderLines',4);
    if length(datatable.Properties.VariableNames) == 3
        datatable.Properties.VariableNames = {'t', 'L','F'};
    elseif length(datatable.Properties.VariableNames) == 4
        datatable.Properties.VariableNames = {'t', 'L','F', 'SL'};
    else
        disp('Wat?')
    end
    % datatable.F = movmean(datatable.F, [8 8]);
    ds.datatable = datatable;
    % sampling frequency. The set starts with 0, thus subtracking one extra element
    ds.fs = (length(datatable.t) - 1)/datatable.t(end); 
    fprintf('Done. Loaded %1.2fs at %0.1f Hz.\n', datatable.t(end), ds.fs);

    % the zero drift is not yet applied!
    % i_ss = datatable.t > datatable.t(end) - 46.0 & datatable.t < datatable.t(end) - 45.5;% index of steady states
    % ss_cur = mean(datatable.F(i_ss));    
    % ds.peaks = max(datatable.F);        
    % ds.ss = ss_cur;

    dsc{ds.set, ds.order} = ds;

    if skipPlots
        continue
    end
    figure(ds.set);
    set(gcf, 'Position',  [769.8000   41.8000  766.4000 740.8000]);

    subplot(211);hold on;plot(datatable.t, datatable.L);
    xlabel('Time (s)');ylabel('Muscle length (L/L0)')
    subplot(212);hold on;
    plot(datatable.t, datatable.F);
    xlabel('Time (s)');ylabel('Tension (kPa)')
    % plot(datatable.t, F);
    % plot(datatable.t(i_ss), repmat(ss_cur, [sum(i_ss), 1])+2, 'LineWidth',2);
    title(ds.datasetLegend, 'Interpreter','None');
    % ylim([0 10])
    % pause
end

%% convert to struct
% fieldnames(dsc{1})
% dsa = cell2mat(dsc)


%% Filter out zero drift in a separate pass
for i_logtrace = 1:size(dsc, 1)
    %% take the log. Log is number 1
    % i_logtrace = 3
    ds = dsc{i_logtrace, 1};

    % is slack? ML below 0.85 definitely is
    is_slack = ds.datatable.L < 0.80; 
    % allow first 2s to be considered as slack too. Usually over zero though
    % is_slack = is_slack | ds.datatable.t < 2; % is slack?
    % is_slack(1) = 0; % we need the diff to be 0 at the beggining

    zdt = ds.datatable.t(is_slack); % zero drift time
    zdF = ds.datatable.F(is_slack); % zero drift Force

    % Force drift as a function of time
    
    % calculate averages of zones for fitting first
    zeroStartZones = find(diff(is_slack) > 0);
    zeroEndZones = find(diff(is_slack) < 0);
    zeroZone_Favg = zeros(1, length(zeroEndZones));
    zeroZone_tavg = zeros(1, length(zeroEndZones));
    for i_zones = 1:length(zeroStartZones)
        zeroZone_Favg(i_zones) = mean(ds.datatable.F(zeroStartZones(i_zones):zeroEndZones(i_zones)));
        zeroZone_tavg(i_zones) = mean(ds.datatable.t(zeroStartZones(i_zones):zeroEndZones(i_zones)));
    end            
    
    % if no slack found
    if length(zeroZone_tavg) == 0
        zeroZone_tavg = [ds.datatable.t(1) ds.datatable.t(end)];
        zeroZone_Favg = [0, 0];
    elseif length(zeroZone_tavg) == 1
        zeroZone_tavg = [ds.datatable.t(1) zeroZone_tavg ds.datatable.t(end)];
        zeroZone_Favg = [zeroZone_Favg, zeroZone_Favg, zeroZone_Favg];
    end

    zeroDriftType = 'piecewiseSpline';
    switch zeroDriftType
        case 'linear'
            f_Fdt = @(a, b, c, x)0*a.*x.^2 +b.*x + c;        
            [ae be] = fit(zdt, zdF, f_Fdt, 'StartPoint', [1e-6, 1e-4, 1]);        
            % plot(timebase, y_exp(ae.a, ae.b, ae.c, timebase), '--', 'Linewidth', 2);
            f_zd = @(t)f_Fdt(ae.a, ae.b, ae.c, t); % zero drift
        case 'quadratic'
            f_Fdt = @(a, b, c, x)a.*x.^2 +b.*x + c;        
            [ae be] = fit(zdt, zdF, f_Fdt, 'StartPoint', [1e-6, 1e-4, 1]);        
            % plot(timebase, y_exp(ae.a, ae.b, ae.c, timebase), '--', 'Linewidth', 2);
            f_zd = @(t)f_Fdt(ae.a, ae.b, ae.c, t); % zero drift
        case 'piecewiseLinear'
            % new drift analysis - spline between locations
            [FO goodness] = fit(zeroZone_tavg', zeroZone_Favg', 'linearinterp');
            % [FO goodness] = fit(zeroZone_tavg', zeroZone_Favg', 'smoothingspline');
            f_zd = @(t)FO(t);
        case 'piecewiseSpline'
            p = pchip(zeroZone_tavg', zeroZone_Favg');
            f_zd = @(t) ppval(p, t);
    end
    % zdS = FO(ds.datatable.t); % zerodrift spline

    figure(i_logtrace); clf; 
    % plot(curve, ds.datatable.t, ds.datatable.F);hold on;
    %     plot(ds.datatable.t(is_slack), ds.datatable.F(is_slack));
    %     ylim([min(ds.datatable.F), max(ds.datatable.F)])
    subplot(222)
    plot(ds.datatable.t, ds.datatable.F, ...
        zdt, zdF, ...
        ds.datatable.t, f_zd(ds.datatable.t));
    yyaxis right;    
    plot(ds.datatable.t, ds.datatable.L, ...
        ds.datatable.t(is_slack), ds.datatable.L(is_slack));
    legend('Raw Force reading', 'Zero regions', ['Zero drift (' zeroDriftType ')']);
    %% xcorr: semimanual data gathering
    %{
    % ignored, results saved
    i_logtrace = 3;
    ful = dsc{i_logtrace, 1}.datatable;
    i = 15; % ramp order in dataset
    srchrng = (i-1)*0.2e4 - 0.1e4 + (1:0.2e4); % correlation range - xcorr does good peaks, but can't really find the one
    % srchrng = (3e0:1.6e3); % correlation range - xcorr does good peaks, but can't really find the one
    seg = dsc{i_logtrace, i}.datatable; 

    % take common sampling rate at 10/s
    Fs = 10;
    timebase = 0:1/Fs:ful.t(end);
    F_fulI = interp1(ful.t, ful.F, 0:1/Fs:ful.t(end));
    F_segI = interp1(seg.t, seg.F, 0:1/Fs:seg.t(end));

    [c, d] = xcorr(F_fulI, F_segI);
    c = c(d>0); d = d(d>0);    % limit only positive shift
    [~, i_xcm] = max(c(srchrng));
    i_xcm = i_xcm + srchrng(1);
    clf;hold on;
    plot(d, c/max(c), 1:length(F_fulI), F_fulI, i_xcm + (1:length(F_segI)), F_segI, 'Linewidth', 2);
    fprintf('Shift for i:%d, t:%0.3f (%s)\n', dsc{i_logtrace, i}.order, timebase(i_xcm), dsc{i_logtrace, i}.filename)
    %}
%% manually assigned shift from 2 to end of all ramps
    ramp_shift_array{i_logtrace} = nan(size(dsc, 2), 1);
    switch dsc{i_logtrace, 1}.folder
        case '20230919'
            ramp_shift_array{1} = [186.8,426.9,576.9,716.9,952.9,1092.9,1232.9,1383,1721.9];
        case '20230927'
            ramp_shift_array{1} = [106.600,276.500,446.600,626.600,896.500,1166.600,1346.500,1516.600];
            ramp_shift_array{2} = [80, 140];
            ramp_shift_array{3} = [134.600, 404.550, 584.500, 754.500, 939.760,1030.500,1300.500,1480.500,1650.500,2066.600,2241.500,2416.500,2591.500];
        case '20230928'
            ramp_shift_array{1} = [99.600, 269.600, 439.600, 619.700, 889.700,1159.700,1339.600,1509.700];
            ramp_shift_array{3} = [108.700, 378.500, 558.600, 728.500, 909.800,1004.500,1274.500,1454.500,1624.500,2040.600,2215.600,2390.600,2565.600];
        case {'20231027', '20231102', '20231107'} 
            ramp_shift_array{1} = [107.500,377.500,557.500,727.500,900.100,1308.500,1478.500,1658.600];
            ramp_shift_array{2} = [68, 139.100];
            ramp_shift_array{3} = [121.300,391.200,571.200,741.200,929.600,1021.200,1291.200,1471.200,1641.200,1811.200,2227.200,2402.200,2577.100,2752.200];
        case {'20240705', '20241010', '20241121', '20241212', '20241217', '20241219', '20241220'}
            ramp_shift_array{1} = [107.500,377.500,557.500,727.500,900.100,1308.500,1478.500,1658.600];
            ramp_shift_array{2} = [68-11.9, 139.100 - 8.9];
            ramp_shift_array{3} = [121.300,391.200,571.200,741.200,929.600-6.38,1021.200,1291.200,1471.200,1641.200,1811.200,2227.200,2402.200,2577.100,2752.200] - 8;
        otherwise
            error('Unknown dataset, check this at line ~251');
    end

    % Identify position of individual ramps
    if ~isnan(ramp_shift_array{i_logtrace})
        % this is weird and has been identified semi-manually
        % ramp_order = [2,3,4,5,6,7,8,9];
        ramp_shift = ramp_shift_array{i_logtrace};
    elseif length(find(diff(is_slack) > 0)) >= 4
        % We use start of the slack, which beggins 4s from the end of indi ramp
        % of the individual cut-out
        i_slackStart = find(diff(is_slack) > 0);
        t_slackStart = (ds.datatable.t(i_slackStart));    
        % take last 4, some traces have some weird beggining
        t_slackStart = t_slackStart([end-3:end]);
        for i = 1:size(t_slackStart, 1)
            ramp_shift(i) = t_slackStart(i) - (dsc{i_logtrace, i + 1}.datatable.t(end) - 36.4);
        end
    else
        % no correction, perhaps no need for a correction
        for i_ramp = 1:size(dsc, 2)-1
            if isempty(dsc{i_logtrace, i_ramp +1})
                % no more ramps here
                break;
            end
            dsc{i_logtrace, i_ramp +1}.datatableZDCorr = dsc{i_logtrace, i_ramp +1}.datatable;
            dsc{i_logtrace, i_ramp +1}.ramp_shift = 0;
        end
        continue;
    end
    % compare to the 100s, 10s, 1s and 0.1s ramp-up cutouts    
    figure(i_logtrace);subplot(221);cla;hold on;
    zeroDrift = cell(1, 5);
    zdr = [];nozdr = [];

    for i_ramp = 1:size(dsc, 2)-1
        if isempty(dsc{i_logtrace, i_ramp +1})
            % no more ramps here
            break;
        end
        rmpdt = dsc{i_logtrace, i_ramp +1}.datatable;
        % TODO pair to global slack zones of the log instead!
        t_slack = rmpdt.t(end) - 31.4 + [0, 8];
        i_slack = rmpdt.t > t_slack(1) & rmpdt.t < t_slack(2);
        avg_Fslack(i_ramp) = mean(rmpdt.F(i_slack));
        plot(rmpdt.t + ramp_shift(i_ramp), rmpdt.L)
        plot(rmpdt.t(i_slack) + ramp_shift(i_ramp), rmpdt.L(i_slack), '*--')
    end
    plot(ds.datatable.t, ds.datatable.L, ':')
    legend('Muscle length', 'Zero regions', 'Whole trace')
    
    figure(i_logtrace);subplot(212);cla;hold on;
    for i_ramp = 1:size(dsc, 2) - 1
        if isempty(dsc{i_logtrace, i_ramp +1})
            % no more ramps here
            break;
        end
        rmp = dsc{i_logtrace, i_ramp + 1};
        rmpdt = rmp.datatable;

        if rmp.ZDR
            % zero drift fitted
            zdr = plot(rmpdt.t + ramp_shift(i_ramp), rmpdt.F - f_zd(rmpdt.t + ramp_shift(i_ramp)), '-', LineWidth=0.5);
            % save the force corrected for the zero drift
            rmpdt.F = rmpdt.F - f_zd(rmpdt.t + ramp_shift(i_ramp));        
        else
            % compare with simple zero shift - e.g. we get oinly single true zer inbetween the baths
            nozdr = plot(rmpdt.t + ramp_shift(i_ramp), rmpdt.F - avg_Fslack(i_ramp), ':', LineWidth=1.5);
            % save the force corrected for the zero shift
            rmpdt.F = rmpdt.F - avg_Fslack(i_ramp);        
        end
        
        dsc{i_logtrace, i_ramp +1}.datatableZDCorr = rmpdt;
        dsc{i_logtrace, i_ramp +1}.ramp_shift = ramp_shift(i_ramp);
    end
    plot(dsc{i_logtrace, 1}.datatable.t, dsc{i_logtrace, 1}.datatable.F - f_zd(dsc{i_logtrace, 1}.datatable.t), 'k:')
    % plot(ds.datatable.t(is_slack), ds.datatable.F(is_slack), '.', 'LineWidth',2)
    % plot(ds.datatable.t, f_Fdt(ae.a, ae.b, ae.c, ds.datatable.t))
    % plot(ds.datatable.t, ds.datatable.F, ':')
    legend([zdr, nozdr], [zeroDriftType ' removal'], 'Zero-value removal')
    axis tight;
end

%% Find peaks and SS
for i_logtrace = 1:size(dsc, 1)
    figure(100 + i_logtrace);clf;hold on;
    for i_ramp = 1:size(dsc, 2) -1
        if isempty(dsc{i_logtrace, i_ramp+1})
            break;
        end
        % use zeroDrift corrected!
        rmpdt = dsc{i_logtrace, i_ramp + 1}.datatableZDCorr;
        rs = dsc{i_logtrace, i_ramp + 1}.ramp_shift;

        
        % i_ss = (rmpdt.t > rmpdt.t(end) - 41.5 - 5 & rmpdt.t < rmpdt.t(end) - 40.5 -5);% indexes of steady states for 20230518
        i_ss = (rmpdt.t > rmpdt.t(end) - 41.5 & rmpdt.t < rmpdt.t(end) - 40.5);% indexes of steady states
        ss_cur = mean(rmpdt.F(i_ss));    
        [peak peakPos] = max(rmpdt.F);
        dsc{i_logtrace, i_ramp + 1}.peak = peak;
        dsc{i_logtrace, i_ramp + 1}.ss = ss_cur;

        
        plot(rmpdt.t + rs, rmpdt.F, ...
            rmpdt.t(i_ss) + rs, rmpdt.F(i_ss), 'x', ...
            rmpdt.t(peakPos) + rs, peak, 'x');

        title(['Peaks and steady state:' dsc{i_logtrace, i_ramp+1}.datasetTitle])

    end
end
%% Save the data file

save(['DataStruct' dsc{1, 1}.folder], 'dsc')

%% plot All the peaks and ss

% rds = [0.1 1 10 100];
figure(21);clf; 
markers = 'sd<^>vox.ph*';
colors = lines(size(dsc, 1));

peaks = nan(size(dsc));
rds = nan(size(dsc));
pCas = nan(size(dsc));
for i_logtrace = 1:size(dsc, 1)
    for i_ramp = 1:size(dsc, 2) -1
        if isempty(dsc{i_logtrace, i_ramp+1})
            break;
        end
        peaks(i_logtrace,i_ramp) = dsc{i_logtrace, i_ramp+1}.peak;
        rds(i_logtrace,i_ramp) = dsc{i_logtrace, i_ramp+1}.rd;
        pCas(i_logtrace,i_ramp) = dsc{i_logtrace, i_ramp+1}.pCa;
    end
    % subplot(211);

    semilogx(rds(i_logtrace, :), peaks(i_logtrace, :), [markers(i_logtrace) '-'], 'LineWidth',2, 'Color', colors(i_logtrace, :));
    legnames{i_logtrace} = dsc{i_logtrace, i_ramp}.datasetTitle;
    hold on;    
end
title('All peaks and SS')
legend(legnames, 'AutoUpdate',false)

% separate loop to separate the legends
ss = zeros(size(rds));
for i_logtrace = 1:size(dsc, 1)
    for i_ramp = 1:size(dsc, 2) -1
        if isempty(dsc{i_logtrace, i_ramp+1})
            break;
        end
        ss(i_logtrace, i_ramp) = dsc{i_logtrace, i_ramp+1}.ss;
        % rds(i_ramp) = dsc{i_logtrace, i_ramp+1}.rd;
    end
    % subplot(211);
    semilogx(rds(i_logtrace, :), ss(i_logtrace, :), [markers(i_logtrace) '--'], 'LineWidth',1, 'Color', colors(i_logtrace, :));
    hold on;    
end

%% compare fast to slow with slow to fast
figure(23);clf;hold on;
colors = lines(3);
rmpsF2S = zeros(0, 2);rmpsS2F = zeros(0, 2);rmpsPNB = zeros(0, 2);
pf2s=[]; ps2f=[]; ppnb = [];pf = [];
% ramp shift by dictionary
rs = dictionary([100, 10, 1, 0.1], [0, 200, 300, 400]);
% order of the ramp in measurements - first come first noted
order = '';

% relax only now
for i_logtrace = 1:size(dsc, 1)
    for i_ramp = 1:size(dsc, 2) -1
        rmp = dsc{i_logtrace, i_ramp+1};
        if isempty(rmp) 
            continue;
        elseif strcmp(rmp.type, 'RelaxF2S')
            if isempty(order)
                order = 'F2S';
            end
            rmpsF2S = [rmpsF2S; rmp.rd, rmp.peak, rmp.ss];
            pf2s = plot(rmp.datatableZDCorr.t + rs(rmp.rd), movmean(rmp.datatableZDCorr.F, [32 32]), '--', 'Color',colors(1, :));
        elseif strcmp(rmp.type, 'RelaxS2F')
            if isempty(order)
                order = 'S2F';
            end
            rmpsS2F = [rmpsS2F; rmp.rd, rmp.peak, rmp.ss];
            ps2f = plot(rmp.datatableZDCorr.t + rs(rmp.rd), movmean(rmp.datatableZDCorr.F, [32 32]), '--', 'Color',colors(2, :));
        elseif strcmp(rmp.type, 'PNBRelax') || strcmp(rmp.type, 'pCa11')
            rmpsPNB = [rmpsPNB; rmp.rd, rmp.peak, rmp.ss];
            ppnb = plot(rmp.datatableZDCorr.t + rs(rmp.rd), movmean(rmp.datatableZDCorr.F, [32 32]), '--', 'Color',colors(3, :));        
        % elseif strcmp(rmp.type, 'pCa11')
        %     rmpsPNB = [rmpsPNB; rmp.rd, rmp.peak, rmp.ss];
        %     ppnb = plot(rmp.datatableZDCorr.t + rs(rmp.rd), movmean(rmp.datatableZDCorr.F, [32 32]), '--', 'Color',colors(3, :));                    
        end
    end
end
title('Ramp order differences - relaxed')    
if strcmp(order, 'F2S')
    legend([pf2s, ps2f, ppnb], 'Fast to slow (first)', 'Slow to fast ', 'PNB slow to fast', 'Location','northeast')
    axes('position', [0.15 0.62, 0.28, 0.3]);cla;
    semilogx(rmpsF2S(:, 1), rmpsF2S(:, 2),  's-', 'Color', colors(1, :), 'LineWidth',2);
    hold on;
    semilogx(rmpsS2F(:, 1), rmpsS2F(:, 2), 'd-', 'Color', colors(2, :), 'LineWidth',2);
    semilogx(rmpsPNB(:, 1), rmpsPNB(:, 2),  's-', 'Color', colors(3, :), 'LineWidth',2);
    legend('fast to slow (first)', 'Slow to fast', 'PNB slow to fast')

elseif strcmp(order, 'S2F') && ~isempty(rmpsF2S)
    legend([pf2s, ps2f, ppnb], 'Fast to slow', 'Slow to fast (first)', 'PNB slow to fast', 'Location','northeast')
    axes('position', [0.15 0.62, 0.28, 0.3]);cla;
    semilogx(rmpsF2S(:, 1), rmpsF2S(:, 2),  's-', 'Color', colors(1, :), 'LineWidth',2);
    hold on;
    semilogx(rmpsS2F(:, 1), rmpsS2F(:, 2), 'd-', 'Color', colors(2, :), 'LineWidth',2);
    semilogx(rmpsPNB(:, 1), rmpsPNB(:, 2),  's-', 'Color', colors(3, :), 'LineWidth',2);
    legend('fast to slow', 'Slow to fast (first)', 'PNB slow to fast')

else
    % F2S is probably empty
    legend([ps2f, ppnb], 'Slow to fast (first)', 'PNB slow to fast', 'Location','northeast')
    axes('position', [0.15 0.62, 0.28, 0.3]);cla;    
    semilogx(rmpsS2F(:, 1), rmpsS2F(:, 2), 'd-', 'Color', colors(2, :), 'LineWidth',2);
    hold on;
    semilogx(rmpsPNB(:, 1), rmpsPNB(:, 2),  's-', 'Color', colors(3, :), 'LineWidth',2);
    legend('Slow to fast', 'PNB slow to fast')
end 

% semilogx(rmpsS2F(:, 1), rmpsS2F(:, 2), 'd-', rmpsF2S(:, 1), rmpsF2S(:, 2),  's-', 'LineWidth',2)
title('Peaks', 'Position',[1, max(rmpsS2F(:, 2)), 0])

%% End of data processing

return;

%% serial stiffness to sarcomere?

colors = lines(size(dsc, 2));
styles = {'-', '--', ':', 's-', 'd--', 'o:', '<-', '^--', '>:'};
figure(41);clf;
leg = cell(0);
for i_logtrace = 1:size(dsc,1)
% i_logtrace = 1
    for i_ramp = 1:size(dsc, 2)-1
        if isempty(dsc{i_logtrace, i_ramp+1})
            break;
        end
    dtst = dsc{i_logtrace, i_ramp+1}; % dataset
    rmp = dtst.datatableZDCorr;
    % limit to ramp plateau
    % wind = rmp.L > 1;
    wind = find(rmp.L > 1.17 & rmp.SL > 2); % window of interest
    wind = wind(1:100:end);
    if ~contains(rmp.Properties.VariableNames, 'SL')
        continue;
    end
    subplot(221);
    semilogx(rmp.t(wind), rmp.SL(wind), styles{i_logtrace}, 'Color', colors(i_ramp, :), LineWidth=2);
    xlabel('t');ylabel('SL');title('SL in time');
    hold on;
    subplot(223);
    semilogx(rmp.t(wind), rmp.F(wind), styles{i_logtrace}, 'Color', colors(i_ramp, :), LineWidth=2);
    xlabel('t');ylabel('F');title('F in time');
    hold on;
    subplot(122);
    plot(rmp.SL(wind), rmp.F(wind), styles{i_logtrace}, 'Color', colors(i_ramp, :), LineWidth=2);
    xlabel('SL');ylabel('F');title('F to SL relation');
    hold on;
    leg{length(leg)+1} =dtst.datasetLegend;
    end
    % pause;
end
legend(leg)


%% Fit the onset
for i_logtrace = [1, 2, 3, 4, 5]
%%
    figure(10 + i_logtrace);clf;hold on;
    colors = lines(size(dsc, 2));
    as = []; bs = []; cs = []; ds = []; rd = []; p = []; pCa = [];

    for i_ramp = 1:size(dsc, 2)-1
        if isempty(dsc{i_logtrace, i_ramp+1})
            break;
        elseif isnan(dsc{i_logtrace, i_ramp+1}.rd)
            continue;
        end
        dtst = dsc{i_logtrace, i_ramp+1}; % dataset
        rmp = dtst.datatableZDCorr;
        offset = 10 + dtst.rd/4;

        fitrg = rmp.L > 1.1 & rmp.t < 10 + dtst.rd;
        rmpFitrg = rmp(fitrg, :);
        % fitfun = @(a, b, c, d, x) min(a*max(x+d, 0).^(b) + c +0*d, 1e2);
        fitfun = @(a, b, c, d, x) a.*(x + b) + a*b*c*d*0;

        [ae goodness] = fit(rmpFitrg.L, rmpFitrg.F,fitfun, 'StartPoint',[1, 1, 1, 0.01]);
        as(i_ramp) = ae.a;
        bs(i_ramp) = ae.b;
        cs(i_ramp) = ae.c;
        ds(i_ramp) = ae.d;
        rd(i_ramp) = dtst.rd;
        pCa(i_ramp) = dtst.pCa;
        % clf;hold on;
        p(i_ramp) = plot(rmpFitrg.L, rmpFitrg.F, '-','Linewidth', 0.5, 'Color', colors(i_ramp, :));
        plot(rmpFitrg.L, fitfun(ae.a, ae.b,ae.c, ae.d, rmpFitrg.L), '-','Linewidth', 4, 'Color', colors(i_ramp, :)*0.8);        
        text(rmpFitrg.L(end), rmpFitrg.F(end), sprintf('a = %1.2e\nb = %1.2e', ae.a, ae.b))
        legends{i_ramp} = sprintf('%s: a = %2.1d with RMSE = %3.1e', dtst.datasetLegend, ae.a, goodness.rmse);
    end
    title(sprintf('Fit for a.*exp(x*b) + c'), 'Interpreter', 'none')
    legend(p(p~=0), legends(p~=0))
    axes('position', [0.15, 0.6, 0.25, 0.3])
    % semilogx(rd, as, 'x-', rd, bs, 'v-', rd, cs, '^-');hold on;
    semilogx(rd, as/max(as), 's:', rds(i_logtrace, :), peaks(i_logtrace, :)/max(peaks(i_logtrace, :)), '^:', 'Linewidth', 2);
    legend('Stiffness Slope (Norm.)', 'Peaks (normalized)');

    axes('position', [0.45, 0.6, 0.25, 0.3])
    plot(-pCa, as/max(as), 's:', -pCas(i_logtrace, :), peaks(i_logtrace, :)/max(peaks(i_logtrace, :)), '^:', 'Linewidth', 2);
    legend('Stiffness Slope (Norm.)', 'Peaks (normalized)');
    
end

%% Compare extracted time-courses
figure(24);clf; hold on;
mult = 1.5;
for i = [1 3 6 7]
    for j = 1:4
        subplot(2,2,j);hold on;
        % Scale equally
        % plot(timecourses{i, j}(:, 1), movmean(timecourses{i, j}(:, 2), [16*j^1.5*2]))
        
        % compare scaling the other experiment set
        if i < 6 % first experiment series
            plot(timecourses{i, j}(:, 1), movmean(timecourses{i, j}(:, 2), [8]))
        else
            plot(timecourses{i, j}(:, 1), mult*movmean(timecourses{i, j}(:, 2), [8]))
        end
        title(sprintf('Ramp %0.1fs', rds(5-j)));
        xlim([10 40 + rds(5-j)])
        ylim([0 40])
    end        
end
subplot(2,2,1);legend(legnames{[1 3 6 7]}, 'Interpreter', 'None');



%% Resample and save

intrs = [2 3 4]
pCas = [11, 6, 4];
% rds = [100, 10, 1, 0.1];
% ts_d = 
clf; hold on;
for i = 1:length(intrs)
    for i_rd = 1:length(rds)
        %%
        clf; hold on;
        rd = rds(i_rd);
        datatable_cur = timecourses{intrs(i), i_rd};

        % correct for zero drift
        datatable_cur.F = datatable_cur.F - f_zd{i_rd};
        
        % cut out the ramp and decay
        % validIds = datatable_cur(:, 1) >= 8 & datatable_cur(:, 1) < rd + 10 + 30;
        validIds = datatable_cur.t >= 8 & datatable_cur.t < rd + 10 + 30;
        datatable_cur = datatable_cur(validIds, :);
        % shift the time base
        datatable_cur(:, 1) =  datatable_cur(:, 1) - 8;


        plot(datatable_cur.t,datatable_cur.F, '--');
        % sample times for datas; baseline, ramp-up, peak decay, long tail
        dwnsmpl = [0:0.5:2, 2 + (0:max(2,rd/30):rd), 2 + rd + (0:max(2,rd/30):min(30,rd)), (2 + 2*rd):0.5:(30 + rd + 2 - 1.0)];

        datatable_cur.Properties.VariableNames = {'Time', 'L', 'F', 'SL'};
        i_data = interp1(datatable_cur.Time, 1:length(datatable_cur.Time), dwnsmpl, 'nearest', 'extrap');
        
        dsf = 8;scaleF = 1;
        datatable_cur.F = movmean(datatable_cur.F*scaleF,[dsf/2 dsf/2]); % force filtered
        % datatable_cur = 
        
        plot(dwnsmpl, datatable_cur.F(i_data), 'x-', 'Linewidth', 2);
        % 
        % plot(dwnsmpl, zeros(1, length(dwnsmpl)), '|-', 'Linewidth', 2);
        % DownSampleAndSplit(datatable, dwnsmpl, [], 2.0, 1, 1.5, filename, 0, 1);

        saveAs = ['bakers_passiveStretch_pCa' num2str(pCas(i)) '_' num2str(rds(i_rd)*1000) 'ms'];
        writetable(datatable_cur, ['data/PassiveCa_3/' saveAs '.csv']);
        disp(['Saved as ' saveAs ' and csv'])
    end
end


%% fft?
