% processes the Ca data, expects peaks from normal
% clearvars -except dataset;

% 
% dataset{1} = load('DataStruct20230518_renamed.mat');
dataset{1} = load('DataStruct20230919.mat');
dataset{2} = load('DataStruct20230927.mat');
dataset{3} = load('DataStruct20230928.mat');
dataset{4} = load('DataStruct20231027.mat');
dataset{5} = load('DataStruct20231102.mat');
dataset{6} = load('DataStruct20231107.mat');
%% Get fmaxes - for each dataset there is a Fmax
% index: dataset, 
% measurement, sequence;
figure(1);clf;hold on;
fmax_pointers = [8, 1;2 1;2, 1;2, 1;2, 1;2,1];
for i = 1:size(dataset, 2)
    if fmax_pointers(i, 1) == 0
        % no fmax measured
        dataset_maxF(i) = NaN;
        continue;
    end
    dtst = dataset{i}.dsc{fmax_pointers(i, 1), fmax_pointers(i, 2)};
    rmp = dtst.datatable;
    % subplot(122);hold on;
    % plot(rmp.t, rmp.L);
    % subplot(121);hold on;
    pfm(i) = plot(rmp.t, rmp.F);

    %  the pCa starts about halfway the dataset
    i_start = find(rmp.t > 100, 1, 'first');
    i_firstSteps = find(rmp.L(i_start:end) > 1.001, 1, 'first') - 5 + i_start;
    % rmp.t([i_start, i_firstSteps])

    % plot(rmp.t(i_start:i_firstSteps), rmp.F(i_start:i_firstSteps), 'x')
    [maxF i_max] = max(rmp.F(i_start:i_firstSteps));
    dataset_maxF(i) = maxF;
    plot(rmp.t(i_max+i_start), maxF, 'rx', LineWidth=4, MarkerSize=12)
end
legend(pfm, '20230919 M','20230927 M','20230928 F','20231027 F','20231102 M','20231107 F');
%% Get the ss F of the slowest active ramp
ss_pointers = [2, 2;3,7;3,7;3,7;3,7;3,7];
dataset_ssF = [];
clf;hold on;
for i = 1:size(dataset, 2)
    if fmax_pointers(i, 1) == 0
        % no fmax measured
        dataset_ssF(i) = NaN;
        continue;
    end
    dttbl = dataset{i}.dsc{ss_pointers(i, 1), ss_pointers(i, 2)}.datatable;
    i_ss = find(dttbl.L > max(dttbl.L)*0.999, 40, 'last') - 10;
    dataset_ssF(i) = mean(dttbl.F(i_ss));
    pl(i) = plot(dttbl.t, dttbl.F); plot(dttbl.t(i_ss), dttbl.F(i_ss), '|');
end
legend(pl, '20230919 M','20230927 M','20230928 F','20231027 F','20231102 M','20231107 F');
%% cell for each ramp
% ramp durations
rds = [100, 10, 1, 0.1];
% rds = [1, 0.1];
% dataset, logtrace, ramp for each ramp duration
% only final dataset
pCa = 4.4;
% pCa = 5.5;
% pCa = 5.75;
% pCa = 6;
% pCa = 6.25;
limit60sPlus = true;

switch pCa
    case 4.4
        pCaData{1} = [1, 2, 2;2, 3, 7;3, 3, 7;4,3, 7;5, 3, 7;6, 3, 7];
        pCaData{2} = [1, 2, 3;2, 3, 8;3, 3, 8;4,3, 8;5, 3, 8;6, 3, 8];
        pCaData{3} = [1, 2, 4;2, 3, 9;3, 3, 9;4,3, 9;5, 3, 9;6, 3, 9];
        pCaData{4} = [1, 2, 5;2, 3,10;3, 3,10;4,3,10;5, 3,10;6, 3,10];
        % indexes of active PNB 300s hold datasets
        i_pca100ms_hold = [2 5; 3 10; 3 10;3 11; 3 11; 3 11];    
        dsName = 'AvgpCa4.4';
    case 5.5
        %% 5.5 - only 100ms in latter experiments
        pCaData{1} = [];pCaData{2} = [];pCaData{3} = [];
        pCaData{4} = [1,3,5;2,3,11;3,3,11;4,3,12;5,3,12;6,3,12];
        i_pca100ms_hold = [3 5; 3 11; 3 11;3 12; 3 12; 3 12];
        dsName = 'AvgpCa5.5';
    case 5.75
        %% 5.75 - only 100ms in latter experiments
        pCaData{1} = [];pCaData{2} = [];pCaData{3} = [];
        pCaData{4} = [1,4,5;2,3,12;3,3,12;4,3,13;5,3,13;6,3,13];
        i_pca100ms_hold = [4 5; 3 12; 3 12;3 13; 3 13; 3 13];
        dsName = 'AvgpCa5.75';
    case 6.0
        %% 6.0 - only 100ms in latter experiments
        pCaData{1} = [];pCaData{2} = [];pCaData{3} = [];
        pCaData{4} = [1,5,5;2,3,13;3,3,13;4,3,14;5,3,14;6,3,14];
        i_pca100ms_hold = [5 5; 3 13; 3 13;3 14; 3 14; 3 14];
        dsName = 'AvgpCa6';
    case 6.25
        %% 6.25 - only 100ms in latter experiments
        pCaData{1} = [];pCaData{2} = [];pCaData{3} = [];
        pCaData{4} = [1,6,5;2,3,14;3,3,14;4,3,15;5,3,15;6,3,15];
        i_pca100ms_hold = [5 5; 3 14; 3 14;3 15; 3 15; 3 15];
        dsName = 'AvgpCa6.25';
end

%% Get corresponding relaxed for scaling options
% Only final dataset with the new protocol
relaxed{1} = [1, 1, 9;2, 1, 6;3, 1, 6;4,1,9;5, 1, 9;6, 1, 9];
relaxed{2} = [1, 1, 8;2, 1, 7;3, 1, 7;4,1,8;5, 1, 8;6, 1, 8];
relaxed{3} = [1, 1, 7;2, 1, 8;3, 1, 8;4,1,7;5, 1, 7;6, 1, 7];
relaxed{4} = [1, 1, 6;2, 1, 9;3, 1, 9;4,1,6;5, 1, 6;6, 1, 6];

%% Adjustment of remaining active force
% figure(6); clf; 
% clf;
% 
% % linear
% % FremFun = @(t, rd, FremMax) min(FremMax(1), (t < rd).*FremMax(1).*(rd-t)./rd) + ... % 
% %             (t >= rd).*FremMax(2).*(t-rd)/(t(end)*0 + 300 -rd);
% % first order
% % k = .015;
% % s = 1;
% 
% % exponential
% % FremFun = @(t, rd, FremMax) min(FremMax(1), (t < rd).*((1-s).*FremMax(1) + s.*FremMax(1).*(rd-t)./rd)) +... % 
% %             (t >= rd).*(...
% %             s*FremMax(2).*...
% %             (1-exp(-k*(t-rd))) +...
% %             (1-s).*FremMax(2)); % exponential approximation
% % s.*FremMax(2)/(1-exp(-k*(t(end) - rd))).*... % scaling to hit the FremMax at the end
% 
% % get the pCa 11 long hold for decay fitting 
% % row = dataset#, col 1 measurement set, col 2 ramp index
% i_Relax_100ms_hold300 = [1 10;1 9;1 9;1 6;1 6;1 6];
% % indexes of validation datasets
% % same or one step slower
% i_validation = [2 4;3 9;3 9;3 10;3 10;3 10];
% % 10s ramp up
% i_validation = [2 3;3 8;3 8;3 8;3 8;3 8];
% i_validation = [];
% tiledlayout(2,5);
% cl = lines(7);
% 
% for i_dtst = 0:length(i_pca100ms_hold)
%     if i_dtst == 0
%         % representative fit
%         nexttile([2 2]);
%         i_dtst = 5;
%         repre = true;
%     else
%         nexttile();
%         repre = false;
%     end
%     dtst = dataset{i_dtst}.dsc{i_pca100ms_hold(i_dtst, 1), i_pca100ms_hold(i_dtst, 2)};
%     dtst_relaxed = dataset{i_dtst}.dsc{i_Relax_100ms_hold300(i_dtst, 1), i_Relax_100ms_hold300(i_dtst, 2)};
%     rmp = dtst.datatableZDCorr;
%     rmpRel = dtst_relaxed.datatableZDCorr;
% 
% 
%     i_0 = find(rmp.t >= 10, 1);
%     % i_end = find(rmp.t >= 10+28+rds(i_rds), 1);
%     dt = (rmp.t(end) - rmp.t(1))/(length(rmp.t)-1);
%     % i_end = find(rmp.L > 1.15, 1, 'last') - 0.5/dt;
%     % int_avg
%     rmp.F = movmean(rmp.F, [8 8]);
%     rmp.t = rmp.t - 10;
%     rmpRes = rmp(1:100:end, :);
%     rmpRelRes = rmpRel(1:100:end, :);
%     rmpRelRes.t = rmpRelRes.t - 10;
%     fitRelRng = rmpRelRes.t > dtst_relaxed.rd & rmpRelRes.t < 300;
% 
%     % fit the relaxed decay first
%     f_fitRelax = @(a, b, c, x) a*(x-dtst_relaxed.rd).^(-b) + c;
%     [rae rgood] = fit(rmpRelRes.t(fitRelRng), rmpRelRes.F(fitRelRng), f_fitRelax, 'StartPoint', [1, 0.01, 5], 'Lower', [0, 0, 0]);    
% 
%     % plot(rmpRel.t-10, rmpRel.F, rmpRelRes.t(fitRelRng), rmpRelRes.F(fitRelRng),...
%     %     rmpRelRes.t(fitRelRng), f_fitRelax(rae.a, rae.b, rae.c, rmpRelRes.t(fitRelRng)))
% 
%     fprintf('Dataset %d: %0.1f(t)^(-%0.2f) + %0.2f\n', i_dtst, rae.a, rae.b, rae.c);    
%     %% remaining force at the beginning and at the end
%     FremMax = [mean(rmp.F(1:i_0-1/dt)), mean(rmp.F(length(rmp.F)-1/dt:end))];
% 
%     % identify the exponential coefficient by fitting
%     t_fit_start = 30;    
%     % we have hold, we can identify the remaining force rise
%     % hold time segmentation not fully reliable
%     % fitrng = rmpRes.t > t_fit_start & rmpRes.t < dtst.holdTime;
%     if t_fit_start > rmpRes.t(end) - 40 - 5
%         % 30s hold only - we basically ignore the rampup and fit a straight line
%         fitRng = rmpRes.t > t_fit_start;
%     else
%         fitRng = rmpRes.t > t_fit_start & rmpRes.t < rmpRes.t(end) - 40;
%     end
% 
%     f_k = @(a, b, c, x) a*(1-exp(-b*x)) + c;
%     f_fit = @(a, b, c, x) f_k(a, b, c, x) + f_fitRelax(rae.a, rae.b, rae.c, x);
%     [ae good] = fit(rmpRes.t(fitRng), rmpRes.F(fitRng), f_fit, 'StartPoint', [1, 0.01, 5], 'Lower', [0, 0, 0], 'Upper', [100, 0.1, 100]);    
% 
%     % Match that by the fit function
%     FremFitShift = FremMax(2) - f_k(ae.a, ae.b, 0, rmpRes.t(end));
%     % minRemForce = f_fit(ae.a, ae.b, FremFitShift, dtst.rd);
%     FremFun = @(t, rd) (t <= 0)*FremMax(1) + ...
%             (t > 0 & t < rd).*max(-10, min(FremMax(1), (f_k(ae.a, ae.b, FremFitShift, rd) - FremMax(1))/rd*t + FremMax(1))) +... % 
%             (t >= rd).*(...
%             f_k(ae.a, ae.b, FremFitShift, t)); % exponential approximation
% 
%     % plot the raw force
%     pltF = plot(rmpRes.t, rmpRes.F, '-', Color=cl(1, :));
%     hold on;
%     % show the fits
%     x_ax = 0.1:0.1:350;
%     if repre
%         pltFitRelax = plot(x_ax, f_fitRelax(rae.a, rae.b, rae.c, x_ax), ':', LineWidth=2, Color=cl(2, :));
%         pltFitRemF = plot(x_ax, f_k(ae.a, ae.b, ae.c, x_ax), ':', LineWidth=2, Color=cl(3, :));
%     end
%     pltFit = plot(x_ax, f_fit(ae.a, ae.b, ae.c, x_ax), ':', LineWidth=4, Color=cl(4, :));
%     xp = rmpRes.t(fitRng);
%     pltFitZone = fill([xp(1) xp(end) xp(end) xp(1)], [0 0 50 50], [1 0.8 0.8], 'FaceAlpha',0.4, EdgeColor='none');
%     % show the remaining force approximation
%     pltApprx = plot(rmp.t, FremFun(rmp.t, dtst.rd), ':', LineWidth=3, Color=cl(5, :));
% 
%     % show the corrected dataset
%     pltCorr = plot(rmpRes.t, rmpRes.F - FremFun(rmpRes.t, dtst.rd), '-', Color=cl(6, :));    
%     % approximation zones
%     % pltB = plot(rmpRes.t(fitrng), rmpRes.F(fitrng), 'x')
%     % pltB = plot(rmp.t(1:i_0-1/dt), rmp.F(1:i_0-1/dt), rmp.t(length(rmp.F)-1/dt:end), rmp.F(length(rmp.F)-1/dt:end), LineWidth=2, Color=[0.8500    0.3250    0.0980]);
%     xlim([rmp.t(1), rmp.t(end)]);
%     if repre
%         legend([pltF, pltFitZone, pltFit, pltFitRelax, pltFitRemF, pltApprx, pltCorr], ...
%             {'Measured data', 'Fitting zone', 'Data fit combined', 'Fit of relaxed decay', 'Fit of F_{rem}', ...
%             'Approximation of F_{rem}', 'Corrected data'}, 'AutoUpdate','off')
%         title({sprintf('Representative remaining force correction for pCa %0.2f', pCa), sprintf('(dataset %d)', i_dtst)});
%         xlabel('Time (s)');ylabel('Tension (kPa)')
%     else
%         title(sprintf('Dataset %d correction', i_dtst));
%     end
%     FremFunArr{i_dtst} = FremFun;
% 
%     %% repeat for the validation ramp, using the same parametrization of the FremFun
%     if length(i_validation) >= i_dtst
%         dtst = dataset{i_dtst}.dsc{i_validation(i_dtst, 1), i_validation(i_dtst, 2)};
%         rmp = dtst.datatableZDCorr;
%         set(gca,'ColorOrderIndex',1)
%         F = movmean(rmp.F, [16 16]);
%         rmp.t = rmp.t - 10;
%         % plot the raw force
%         pltF = plot(rmp.t, F, ':');
%         % show the remaining force approximation
%         pltApprx = plot(rmp.t, FremFun(rmp.t, dtst.rd), '--');
%         % show the corrected dataset
%         pltCorr = plot(rmp.t, F - FremFun(rmp.t, dtst.rd), '-');    
%     end
%     set(gca, 'FontSize', 14);
%     ylim([-5, 50])
%     % %% Test fit the corrected one
%     % FcBarca = rmpRes.F - FremFun(rmpRes.t, dtst.rd);
%     % f_fitRelaxt = @(a, b, c, x) rae.a*(x-dtst_relaxed.rd).^(-rae.b) + c + a*b*c*0;
%     % [raet rgoodt] = fit(rmpRes.t(fitRng)-dtst.rd, FcBarca(fitRng), f_fitRelaxt, 'StartPoint', [1, 0.01, 5], 'Lower', [0, 0, 0]);    
%     % % figure(21);clf;
%     % plot(rmpRes.t, FcBarca);
%     % hold on;
%     % % plot(rmpRes.t(fitRng), FcBarca(fitRng));
%     % x = 0.1:0.1:300;
%     % plot(x+dtst.rd, f_fitRelaxt(rae.a, rae.b, raet.c, x), '-', 'Linewidth', 4);
%     % 
%     % plot(x+dtst.rd, f_fitRelaxt(raet.a, raet.b, raet.c, x), '-', 'Linewidth', 4);
% end

%% Average all data
% figure(5);clf;

figure(4);clf;hold on;

clear sp peaks peaks_norm leg;

for i_rds = 1:length(rds)
    sp{i_rds} =subplot(4, 4, (i_rds-1)*4 +  (1:2));cla;
    % % excluding 100s
    % if i_rds ~= 1 % skip plotting of 100s
    %     sp = subplot(3, 4, (i_rds-2)*4 +  (1:2));cla;
    % end
    outF = [];sum_squared_diff = []; n = 1;clear leg;
    % clear output arr
    Farr{i_rds} = outF;
    rampSet = pCaData{i_rds};
    if isempty(rampSet)
        continue;
    end
    clin = lines(size(rampSet, 1)+1);

    for i_logtrace = 1:size(rampSet, 1)
        if limit60sPlus && rampSet(i_logtrace, 1) == 1                 
            % we want only decays 60s+            
            continue;
        end

        % if rampSet(i_logtrace, 1) == 0
        %     % placeholder to match the relaxed peaks, skipping            
        %     continue;
        % end
        % experiment dataset - whole measurement session
        eds = dataset{rampSet(i_logtrace, 1)}.dsc;
        % dataset - particular conditions
        dtst = eds{rampSet(i_logtrace, 2), rampSet(i_logtrace, 3)};
        rmp = dtst.datatableZDCorr;
        i_0 = find(rmp.t >= 10, 1);
        % i_end = find(rmp.t >= 10+28+rds(i_rds), 1);
        dt = (rmp.t(end) - rmp.t(1))/(length(rmp.t)-1);
        i_end = find(rmp.L > 1.15, 1, 'last') - 0.5/dt;
        
        % trim the ramp - start to end
        rmp = rmp(1:i_end, :);
        % rmp = rmp(i_0:i_end, :);
        rmp.t = rmp.t - 10;
        % i_0 = 1;

        
        
        %% filter out the remaining force
        % saving the original 
        rmp.Fraw = rmp.F;
        % rmp.F = rmp.F - FremFunArr{i_logtrace}(rmp.t, dtst.rd);            
        
        % base on nothing - just absolute peak values
        base_rel = 1;
        % base on absolute peak of each ramp
        % base_rel = max(rmp.F);
        % base on absolute peak of fastest ramp
        % base_rel = ...
        %     dataset{pCa4_4Data{4}(i_logtrace, 1)}.dsc{pCa4_4Data{4}(i_logtrace, 2), pCa4_4Data{4}(i_logtrace, 3)}.peak;
        % base on steady state
        % base_rel = rmp.F(end);
        % base on relaxed peaks for each ramp
        % base_rel = eds{relaxed{i_rds}(i_logtrace, 2),...
                        % relaxed{i_rds}(i_logtrace, 3)}.peak;
        % base on relaxed peaks for fastest ramp only
        % base_rel = eds{relaxed{4}(i_logtrace, 2),...
                        % relaxed{4}(i_logtrace, 3)}.peak;
        % base on Fmax
        base_rel = dataset_maxF(rampSet(i_logtrace, 1));
        % base on staeady state of the slowest ramp (in AverageRampsCa)
        % base_rel = dataset_ssF(i_logtrace);

        if isnan(base_rel) || base_rel == 0
            continue;
        end

        % adjusted for Fremaining, normalized
        F = rmp.F/base_rel;
        % just raw force, normalized
        % F = rmp.Fraw/base_rel;
        % just raw force
        % F = rmp.F;

        % prepare legend
        leg{i_logtrace} = [dtst.folder ':' dtst.datasetTitle];        
        % semilogx(rmp.t, movmean(F, [16 16]), ':', Color=clin(i_logtrace, :));hold on;
        semilogx(rmp.t, F, '-', Color=[clin(i_logtrace, :), 0.05]);hold on;
        % semilogx(rmp.t, FremFunArr{i_logtrace}(rmp.t, dtst.rd)/base_rel, '--', Color=clin(i_logtrace, :));
        % semilogx(dtst.datatableZDCorr.t - 10, dtst.datatableZDCorr.F, '-', Color=clin(i_logtrace, :));hold on;
        % semilogx(rmp.t, FremFun(rmp.t, dtst.rd), '--', Color=clin(i_logtrace, :));
        set(gca, 'FontSize', 14);
        % base on absolute peak
        peaks(i_rds, i_logtrace) = max(F*base_rel);
        % base on peak over steady state
        % peaks(i_rds, i_logtrace) = max(rmp.F) - base_rel;
        
        % normalization        
        peaks_norm(i_rds, i_logtrace) = max(F);

        if isempty(outF)
            outF = F;
            outFr = rmp.Fraw;
            outL = rmp.L;
            Fmax = base_rel;
            n = 1;
            sum_squared_diff = outF*0;
        else
            L = min(length(outF), length(rmp.F));
            n = n + 1;
            % shrink to shortest one
            outF  = outF(1:L)  + (F(1:L)     - outF(1:L))/n;
            outFr = outFr(1:L) + (rmp.Fraw(1:L) - outFr(1:L))/n;
            outT = rmp.t(1:L);
            outL  = outL(1:L)  + (rmp.L(1:L)     - outL(1:L))/n;
            % avg the peaks to rescale back
            Fmax = Fmax + (base_rel - Fmax)/n;

            % estimate the standard error out of sum squared
            sum_squared_diff = sum_squared_diff(1:L) + (F(1:L) - outF(1:L)).^2;            
        end    

        %% getting stiffnes based on linear fit
        % fitrg = rmp.L > 1.1 & rmp.t < dtst.rd;
        % rmpFitrg = rmp(fitrg, :);
        % % fitfun = @(a, b, c, d, x) min(a*max(x+d, 0).^(b) + c +0*d, 1e2);
        % fitfun = @(a, b, x) a.*(x + b);
        % 
        % [ae goodness] = fit(rmpFitrg.L, rmpFitrg.F,fitfun, 'StartPoint',[1, 1]);
        % as(i_rds, i_logtrace) = ae.a;
        % figure(5);
        % semilogx(rmp.t + rampShift(i_rds), rmp.Fraw, '-', Color=[clin(i_logtrace, :), 0.1]);hold on;
        % figure(4);                
    end
    
    % figure(5);
    % rampShift = [-92, -8.6, -0.58, 0];
    % semilogx(outT + rampShift(i_rds), outFr, 'k-', LineWidth=1);hold on;
    % xlim([1e-2 1e2])
    % figure(4);
    % Tarr{i_rds} = t_s;
    % Farr{i_rds} = interp1(outT, outFr, t_s, "pchip", 'extrap');
    Tarr{i_rds} = outT;
    Larr{i_rds} = outL;
    % Farr{i_rds} = outFr;
    Farr{i_rds} = outF*Fmax;
    FmaxArr{i_rds} = Fmax;
    SDarr{i_rds} = sqrt(sum_squared_diff / (n * (n - 1)));
end

%% Frem correction of the averaged
cl = lines(6);
% using the AverageRamps.m - just for the figures here
pca11data = load('pca11data.mat');
FarrRel = pca11data.Farr; TarrRel = pca11data.Tarr;
% load pca4data60sNoFremCorr.mat
% identified from fitting the no Ca ramps FitDecayOverlay.m
rampShift = [5.3980    0.8234    0.2223   0.0100];
xRel = [4.4271    0.2121    4.8964];
ref60s = NaN;
f = figure(9); clf;
aspect = 1.5;
f.Position = [300 200 7.2*96 7.2*96/aspect];
tiledlayout(1, 4, "Padding","compact", "TileSpacing","tight")
% drawPlots = [4 3 1];
for i_rds = [4 3 2 1]
    if i_rds == 0
        i_rds = 4;
        nexttile(1, [1 1])
        title('Ramp: 0.1s')
    else
        nexttile();
    end
    clear rmp rmpRel rmpRelRes;
    if isempty(Farr{i_rds})
        FarrCorr{i_rds} = [];
        continue;
    end
    rmp = table(); rmpRel = table();

    rmp.F = movmean(Farr{i_rds}, [0 0]);
    rmp.t = Tarr{i_rds} - rds(i_rds);
    rmpRes = rmp(1:50:end, :);
    rmpRel.F = movmean(FarrRel{i_rds}, [8 8]);
    rmpRel.t = TarrRel{i_rds} - rds(i_rds);
    rmpRelRes = rmpRel(1:50:end, :);
    fitRelRng = rmpRelRes.t > 0;

    %% fit the relaxed decay first
    f_fitRelax = @(x) xRel(1)*max(1e-6, (x + rampShift(i_rds))).^(-xRel(2)) + xRel(3);
    % [rae rgood] = fit(rmpRelRes.t(fitRelRng), rmpRelRes.F(fitRelRng), f_fitRelax, 'StartPoint', [10, 0.1, 5, 0.001], 'Lower', [0, 0, 0, 0]);    

    % fprintf('Dataset %d: %0.1f(t)^(-%0.2f) + %0.2f\n', i_rds, rae.a, rae.b, rae.c);    
    %% remaining force at the beginning, assuming it is the same at the end
    % sel = rmp.t < - rds(i_rds) & rmp.t > -1 - rds(i_rds);
    % plot(rmp.t(sel), rmp.F(sel))
    FremMax = mean(rmp.F(rmp.t < - rds(i_rds) & rmp.t > -1 - rds(i_rds)));

    % align the Frem0?
    % if isnan(ref60s)
    %     ref60s = FremMax;
    % end
    % shiftRel = ref60s - FremMax;
    % rmpRes.F = rmpRes.F + shiftRel;
    % FremMax = FremMax + shiftRel;

    %% identify the exponential coefficient by fitting
    t_fit_start = 15;    
    fitRng = rmpRes.t > t_fit_start & rmpRes.t <= rmpRes.t(end);
    
    f_k = @(a, b, c, x) a*(1-exp(-b*x)) + c;
    f_fit = @(a, b, c, x) f_k(a, b, c, x) + f_fitRelax(x);
    [ae good] = fit(rmpRes.t(fitRng), rmpRes.F(fitRng), f_fit, 'StartPoint', [1, 0.01, 5], 'Lower', [0, 0, 0], 'Upper', [1000, 0.1, 100]);    
    fprintf('> fit goodness: %0.3f, Frem time constant: %0.4f \n', good.rmse, ae.b);
    
    % what is the shift, if we fix the end - the force rise must coincide
    % with the FremMax
    FremFitShift = FremMax - f_k(ae.a, ae.b, 0, rmpRes.t(end));
    FremFun = @(t, rd) (t <= -rd)*FremMax + ...
            (t > -rd & t < 0).*((f_k(ae.a, ae.b, FremFitShift, 0) - FremMax)/rd*t + FremFitShift) +... % 
            (t >= 0).*(...
            f_k(ae.a, ae.b, FremFitShift, t)); % exponential approximation

    
    % show the fits
    cla;
    x_ax = 0.1:0.1:rmpRes.t(end);

    pltF = plot(rmpRes.t, rmpRes.F, '-', Color=cl(1, :),LineWidth=1 ); hold on;        
    % pltFit = plot(x_ax, f_fit(ae.a, ae.b, ae.c, x_ax), ':', LineWidth=4, Color=cl(1, :)*0.8);
    % pltRelax = plot(rmpRelRes.t, rmpRelRes.F, Color=cl(2, :));
    % pltFitRelax = plot(x_ax, f_fitRelax(x_ax), ':', LineWidth=3, Color=cl(2, :));
    pltFitRemF = plot(x_ax, f_k(ae.a, ae.b, ae.c, x_ax), ':', LineWidth=3, Color=cl(3, :));
    
    xp = rmpRes.t(fitRng);
    % pltFitZone = fill([xp(1) xp(end) xp(end) xp(1)], [-10 -10 60 60], [1 0.8 0.8], 'FaceAlpha',0.24, EdgeColor='none');
    
    % show the remaining force approximation
    pltApprx = plot(rmpRes.t, FremFun(rmpRes.t, rds(i_rds)), '-.', LineWidth=3, Color=cl(4, :));
    
    % show the corrected dataset
    % clf;
    pltCorr = loglog(rmpRes.t, rmpRes.F - FremFun(rmpRes.t, rds(i_rds)), '-', Color=cl(5, :), LineWidth=3);    
    % pltFrem0 = plot(-1 - rds(i_rds), FremMax, 'x', MarkerSize=10, LineWidth=3, Color=cl(4, :));
    
    %
    % approximation zones
    xlim([rmpRes.t(1), rmpRes.t(end)]);
    %
    FarrCorr{i_rds} = rmp.F - FremFun(rmp.t, rds(i_rds));
    TarrCorr{i_rds} = rmp.t + rds(i_rds);
    title(sprintf('Ramp %g s',rds(i_rds)), 'interpreter', 'LaTex');


    ylim([-2 48]);
    if i_rds < 4
        yticks([]);
    else
        yticks([0 10 20 30 40]);
        ylabel('$\Theta$ (kPa)', Interpreter='latex')
    end
    xlabel("$t - t_r$ (s)", Interpreter="latex")
    if i_rds == 1
        % leg = legend([pltF, pltFit, pltApprx,...
        %     pltRelax, pltFitRelax, pltFrem0,  ...
        %     pltCorr, pltFitRemF, pltFitZone],...
        %     "$\Theta_A$", "$\Theta_F$", '$\Theta^*$', ...
        %     "Relaxed stress", "$\sim C_1 (t-t_r)^{-\alpha}$","$\Theta_\infty$",...
        %     "$\Theta_{corr} = \Theta_A - \Theta^*$", "$\sim C_2(1-e^{-b(t-t_r)})$", "$t - t_r > 15$",...             
        %     Location="northeast", NumColumns=3, Interpreter="latex");
        leg = legend([pltF, pltApprx,...            
            pltCorr, pltFitRemF ],...
            "$\Theta_A$", '$\Theta^*$', ...
            "$\Theta_{corr} = \Theta_A - \Theta^*$", "$\sim C_2(1-e^{-b(t-t_r)})$", ...             
            Location="northeast", NumColumns=2, Interpreter="latex");
        
        leg.ItemTokenSize = [12 12]
    end
end
fontsize(12, "points")
f = gcf();
exportgraphics(f,'Figures/Frem_Correction.png','Resolution',150)


%%
% Farr = FarrCorr;
% Tarr = TarrCorr;
% save('pCa4dataNoAdj60sFremCorrShifted.mat', 'Tarr','Farr')
if pCa == 4.4
    save('pCa4dataNoAdj60sFremCorr.mat', 'Tarr','FarrCorr', 'Farr')
end

%% resample and save and plot the AVG
for i_rds = [4 3 2 1]
    if isempty(FarrCorr{i_rds})
        continue;
    end
    t = TarrCorr{i_rds};
    t_ignore = 0.001; 
    % resample up the peak and for the tail separately
    % t_s = [linspace(0, 1, 20)*(rmp.t(i_0) + rds(i_rds) - t_ignore) ...
    %     logspace(log10(rmp.t(i_0) + rds(i_rds) + t_ignore), log10(rmp.t(i_0) + rds(i_rds) + 10), 20) ...
    %     logspace(log10(rmp.t(i_0) + rds(i_rds) + 10), log10(outT(end)), 5)];

    % resample up the peak and for the tail separately
    t_s = [linspace(0, 1, 20)*(rds(i_rds) - t_ignore) ...
        logspace(log10(rds(i_rds)), log10(rds(i_rds) + min(60, t(end) - rds(i_rds))), 80)];
            % rds(i_rds) + linspace(t_ignore, min(60, t(end) - rds(i_rds)), 40)];...        
           % rds(i_rds) + logspace(log10(t_ignore), log10(min(60, t(end) - rds(i_rds))), 40)];

    % resample log equally
    % t_s = [logspace(log10(1e-3), log10(outT(end)), 40)];

    % remove consequent duplicates at joints
    t_s = t_s(~[false t_s(2:end) == t_s(1:end-1)]);
    
    % force and length interpolation
    % save denormalized
    FLSDint = interp1(t, [FarrCorr{i_rds}, Larr{i_rds}, SDarr{i_rds}*FmaxArr{i_rds}], t_s, "pchip", 'extrap');
    figure(7);hold on;plot(t_s, FLSDint(:, 1), 's-');
    tab_rmpAvg = table(t_s' + 2, FLSDint(:, 2), FLSDint(:, 1), FLSDint(:, 3));
    tab_rmpAvg.Properties.VariableNames = {'Time', 'L', 'F', 'SD'};
    writetable(tab_rmpAvg, ['data/' dsName '_' num2str(rds(i_rds)) 's.csv']);
%% plot the AVG

    % Fill the area between upper and lower bounds to show the confidence interval
    % fill([outT', fliplr(outT')], [(outF*Fmax - SD*1.96)', fliplr((outF*Fmax + SD*1.96)')], 'b');

    % semilogx(t_s, FLSDint(:, 1) + FLSDint(:, 3), '-', Color=[clin(end, :), 1], LineWidth=2);
    % semilogx(t_s, FLSDint(:, 1) - FLSDint(:, 3), '-', Color=[clin(end, :), 1], LineWidth=2);
    % semilogx(t_s, FLSDint(:, 1), '-|', LineWidth=2, Color=clin(end, :))
    
    % plot normalized
    axes(sp{i_rds})
    errorbar(t_s,FLSDint(:, 1)/FmaxArr{i_rds},FLSDint(:, 3)/FmaxArr{i_rds}, '-', LineWidth=2, Color=clin(end, :))
    xlim([1e-2 2e2])
    yl = ylim;
    ylabel('\itT_{m,rel}');
    yyaxis right; ylim(yl*Fmax);
    ylabel('\itT_m (kPa)');
    yyaxis left;
    % plot(rmp.t(1:L), outF(1:L) + (rmp.F(1:L) - outF(1:L))/n)
    leg{length(leg) +1} = 'Averaged';
    title(['Ramp ' num2str(rds(i_rds)) 's'])
    if i_rds == 1
        % x = 0.1;y = 0.03;
        % legend(leg, 'Interpreter','none', 'Position', ...
        %     [sp.Position(1) + sp.Position(3) + x, sp.Position(2) - y, 1 - sp.Position(1) - sp.Position(3) - x, sp.Position(4)])
    elseif i_rds == 4
        xlabel('Time (s)')
    end    
end
peaks44 = peaks;
%% peak sum up figure
peaks = peaks44;
% normtype = 'Norm to highest peak'
% normtype = 'Norm to steady state'
% normtype = 'Norm to relaxed peaks for each ramp'
% normtype = 'Norm to highest relaxed peak'
% normtype = 'Norm to Fmax'
% subplot(4, 4, [3 16]);cla;
aspect = 2;
f = figure(3); f.Position = [800 200 7.2*96 7.2*96/aspect];

gc = axes('Position', [0.6 0.1 0.39 0.8], 'Box','on', 'BoxStyle','full', 'Color','r');

% figure(11);
% subplot(1, 5, 5);cla reset;
% remove zero peaks as NaNs to fix the average - only when something was missing
peaks(peaks == 0) = NaN;

% for x axis we go from fastest to slowest
x_ax = length(rds):-1:1;


boxplot(fliplr(peaks'), PlotStyle="traditional", Notch="off", Color='k');hold on;

% coefficient of variance
cv = nanstd(peaks_norm')./nanmean(peaks_norm');
fprintf('Coefficient of variance: %1.3f \n', sum(cv));
% colors
clin = lines(size(peaks, 2)+1);
% BW only
clin = repmat([0 0 0], [size(peaks, 2)+1, 1]);
if limit60sPlus
    % shift by one bc we skipped
    clin(2:end, :) = clin(1:end-1, :);
end
for i_pk = 1:size(peaks_norm, 2)
    if all(isnan(peaks_norm(:, i_pk)))
        % it was just a placeholder
        leg{i_pk} = '';
        continue;
    end
    % set(gca, 'colororderindex', 1);
   
    
    plot(x_ax, peaks(:, i_pk), '.:', 'MarkerSize',12, LineWidth=1.5, Color=clin(i_pk, :));hold on;    
    % semilogx(rds, peaks_relaxed(:, i_pk), 'x-', 'MarkerSize',12, LineWidth=0.5, Color=clin(i_pk, :));hold on;
    % semilogx(rds, as(:, i_pk), 'x:', 'MarkerSize',12, LineWidth=0.5, Color=clin(i_pk, :));hold on;
end

% just for the legend
plot(NaN, NaN, 's-',LineWidth=2, Color=clin(end, :), MarkerSize=5)
%# find non-empty legend
% validLeg = ~cellfun(@isempty,leg);
% %# remove empty cells
% leg_cleared = leg(validLeg);
leg_cleared = {'20230919 M','20230927 M','20230928 F','20231027 F','20231102 M','20231107 F', 'Averaged'};
leg_cleared = {'20230927 M','20230928 F','20231027 F','20231102 M','20231107 F', 'Averaged'};

% legend(leg_cleared, 'Interpreter','none', 'AutoUpdate','off', 'Location','northeast')

plot(x_ax, nanmean(peaks, 2)', '-', LineWidth=2, Color=clin(end, :))
plot(x_ax, nanmean(peaks, 2)', 's', LineWidth=3, Color=clin(end, :), MarkerSize=5)

% plot(x_ax, nanmedian(peaks_norm, 2)', '--', LineWidth=3, Color=clin(end, :))
% plot(x_ax, nanmedian(peaks_norm, 2)', '_', LineWidth=5, Color=clin(end, :), MarkerSize=12)


xlabel('$t_r$ (s)', Interpreter='latex');
% xlim([0.08, 150])
% ylabel('$T_{rel}$ (kPa)', Interpreter='latex')
% ylim([0 1]);
% yl = ylim();
% yyaxis right;ylim(yl*Fmax);
ylabel('$\Theta$ (kPa)', Interpreter='latex');
% semilogx(rds, mean(peaks, 2)', '_', LineWidth=3, MarkerSize=12)
% set(gca, 'XTick', fliplr(rds));
set(gca, 'XTickLabel', {'0.1', '1', '10', '100'});
g = gca();
% g.YAxis(2).Color = [0 0 0];
ylim([0 inf])
set(gca, 'FontSize', 12);
% title(sprintf('Peak stress - pCa %0.2g', pCa));
title(sprintf('Peak stress - pCa %g', 4.51));
% exportgraphics(f,sprintf('Figures/AvgpPeakspCa%0.2g.png', pCa),'Resolution',150)
exportgraphics(f,sprintf('Figures/AvgpPeaks.png'),'Resolution',150)