% Runs the passive model - could be run separately in a clear environment
% or by OptimizeCOmbined.m

%% PLoting and visualization options
if ~exist('drawPlots', 'var')
    drawPlots = true;
end
if ~exist('drawAllStates', 'var')
    % indicates the ramp for which we plot the states
    drawAllStates = 0;
end
if ~exist('plotDetailedPlots', 'var')
    plotDetailedPlots = false;
end
if ~exist('drawFig1', 'var')
    drawFig1 = false;
end
if ~exist('exportRun', 'var')
    exportRun = false;
end
if ~exist('rerunFitting', 'var')
    rerunFitting = false;
end
if ~exist('plotInSeparateFigure', 'var')
    plotInSeparateFigure = true;
end

%% model inputs settings
if ~exist('pCa', 'var')
    % pCa is one of these: [4.51, 5.5, 5.75, 6, 6.2, 11]
    pCa = 11;    
end
if ~exist('params', 'var')    
    paramNames = {'\F_{ss}', 'n_ss' , 'k_p'    , 'n_p'    , 'k_d'   , 'n_d'    , '\alpha_U'  , 'n_U'    , '\mu'    , '\delta_U'    , 'k_{A}'   , 'k_{D}'   };

    % Nice fit of everything    
    paramSet = [...
      5.19       12.8       4345       2.37      4e+04       2.74  8.658e+05      5.797      0.678      0.165   0.005381      0.383 
      5.19       12.8       3998       2.37      4e+04       2.74  2.837e+06      6.211      0.678      0.165   0.005213      0.383 
      5.19       12.8       3109       2.37      4e+04       2.74  3.184e+06      6.526      0.678      0.165   0.002699      0.383 
      5.19       12.8       1623       2.37      4e+04       2.74  1.807e+07       7.96      0.678      0.165  3.261e-05      0.383 
      5.19       12.8      887.7       2.37      4e+04       2.74  1.876e+07      8.505      0.678      0.165  3.478e-06      0.383 
      5.19       12.8      512.3       2.37      4e+04       2.74  2.668e+07      9.035      0.678      0.165        NaN        NaN 
    ];
   pcax = [4.51, 5.5, 5.75, 6, 6.2, 11];
   i_pcax = find(pCa == pcax);
   params = paramSet(i_pcax, :);
end
if ~exist('rampSet', 'var')
    rampSet = [1 2 3 4];
    % rampSet = [4]; % fastest only
end
if ~exist('simtype', 'var')
    % simtype = 'sin';
    % simtype = 'rampbeat';
    simtype = 'ramp';
end

if ~exist('datasetname', 'var')
    datasetname = 'pnbmava';
    % datasetname = 'StitchingRelaxed_1_1.csv';
end
if drawAllStates && ~exist('time_snaps', 'var')
    % decide for timepoints
    % fixed time or fraction of ramp durations?
    % time_snaps = [0, 0.1, 1, 10, 30, 40, 100]
    % time_snaps = [0, rds(j), rds(j) + 30, 60, 120, 160];
    if strcmp(simtype, 'ramp')
        time_snaps = [1e-3, 1, 2, 10];
    elseif strcmp(simtype, 'sin')
        time_snaps = [Tc/2, Tc, 3*Tc/2, 2*Tc];
    end
end

if ~exist('alphaF_0', 'var')
    if length(params) >= 13
        alphaF_0 = params(13); % refolding constant - refolding enabled
    else
        alphaF_0 = 0; % refolding constant - refolding disabled
    end
end

%% preapre environment
load SoHot.mat;
% Cool colorbar creation, stored in mat file though
%     colors = [
%         1, 1, 1; % White for zero
%         0, 0.5, 1; % Blue (complementary to orange)
%         0, 0, 1; % Dark blue
%         0, 0, 0; % Black for max
%     ];
%     % Interpolate the colors
%     SoCool = interp1(linspace(0, 1, size(colors, 1)), colors, linspace(0, 1, 256));
% save SoCool SoCool;
load SoCool.mat;
datasetname = 'pnbmava';
colors = gray(5);

pCa_noEffect = 11; % at 11 we do not consider PEVK attachment



clear Force
clear Time
clear Length
clear outStruct;
if any(params < 0) 
    cost = inf;
    return;
end

figInd = get(groot,'CurrentFigure'); % replace figure(indFig) later without stealing the focus

%% load data files
% pCa 11 - load Relaxed, do not run the extended, PEVK attachment model
% pCa 10 - load Relaxed, run the PEVK attachment model
% pCa < 10 - load AvgpCa dataset, Ca effect in place

% rds = fliplr([0.02 0.1, 1, 10 100]);
rds = fliplr([0.1, 1, 10, 100]);
% rds = fliplr([0.1, 10]);
for i_rd = 1:length(rds)
  if isinf(pCa) || pCa >= 10
      if datasetname == "pnbonly"
          filename = ['..\Data\AvgRelaxed_' num2str(rds(i_rd)) 's.csv'];
      elseif datasetname == "pnbmava"
          filename = ['..\Data\AvgRelaxedMavaSet_' num2str(rds(i_rd)) 's.csv'];
      end

  else
      % file names are somewhat different
      switch (pCa)
          case 4.51
            file_pCa = 4.4;
          case 5.75
            file_pCa = 5.8;
          otherwise
            file_pCa = pCa;
      end
      if datasetname == "pnbonly"
        filename = ['..\Data\AvgpCa' num2str(file_pCa) '_' num2str(rds(i_rd)) 's.csv'];
      elseif datasetname == "pnbmava"
        filename = sprintf('..\\Data\\AvgMava_pCa%0.1f_%gs.csv', file_pCa, rds(i_rd));     
      end
  end

  if exist([filename], "file")
      datatables{i_rd} = readtable(filename);
  else
      datatables{i_rd} = [];
  end
  
end

%% Model numerical parameters

% half-sarcomere ramp height
Lmax = 1.175 - 0.95; % Lmax = 0.225;
Nx   = 15;          % number of space steps
ds   = 1*(Lmax)/(Nx-1);      % space step size
% use this to make sure the ds is big wnough
% ds   = 1*(Lmax + 0.015)/(Nx-1);
s  = (0:1:Nx-1)'.*ds; % strain vector
Ng = 10; 

%% model input adjustable parameters
% new params assignments
Fss  = params(1 );		% Fss   Steady state level
n_ss = params(2 );		% n_ss  Steady state level exponent
kp   = params(3 );		% kp    proximal chain force constantkS   =
np   = params(4 );		% np    proximal chain force exponent
kd   = params(5 );		% kd    proximal chain force constant high Cabist
nd   = params(6 );		% nd    distal chain force exponent
alphaU  = params(7 );		% a_U   chain unfolding rate constant
nU   = params(8 );		% nU    unfolding rate exponent
mu   = params(9 );		% mu    small enough not to affect the result
delU = params(10)/Ng;  % delU 
kA   = params(11);     % kA   
kD   = params(12); 	% kD    PEVK detachment rate

kDf = 0; % universal modifier
L_0  = 1; % reference sarcomere length (um)

% Calculate proximal globular chain force Fp(s,n) for every strain and
% value. 
slack = (0:Ng).*delU;
% Fp = kp*(max(0,s-slack)/Lref).^(np); 
Fp = kp*(max(0,s-slack)/L_0).^(np); 
% Calculate the globular chain folding/unfolding probability transition
% rates
% RU = alphaU*(max(0,s-slack(1:Ng))).^nU; % unfolding rates from state n to (n+1)
RU = alphaU*((max(0,s-slack(1:Ng))/L_0).^nU).*(ones(Nx,1).*(Ng - (0:Ng-1))); % unfolding rates from state n to (n+1)
% clf;mesh(RU)
% drawFig1 = false;
%% visualizing the Force plot - fig 1A&B
if drawFig1
    % Fp = kp*(max(0,s-slack)/Lref).^(np); 

    % do not plot zeros
    % Fp(Fp < 1e-1) = NaN;
    % sx = s(1):0.1:s(end);
    sx = s;
    % Fp = interp1(s, Fp(:, :), sx);
    % f = figure(3); 
    f = gcf();
    set(gcf, 'Position', [500  300  7.2*96 3.5*96])
    gc = axes('Position', [0.1, 0.2, 0.35, 0.7]);
    pl1 = plot(sx, Fp(:, 1), 'k-', 'linewidth', 2); hold on;
    pl2 = plot(sx, Fp(:, 2), 'k-.', 'linewidth', 2); hold on;
    pl3 = plot(sx, Fp(:, 3), 'k--', 'linewidth', 2); hold on;
    plN = plot(repmat(sx, [1, Ng-3]), Fp(:, 4:Ng), ':', 'linewidth', 2, 'Color', 'k'); 
    
    % arrow
    % t0 = 5;
    % plot(s([14, 24]), [t0 t0], 'k');scatter(s(14), [t0], 'k>', 'filled');
    % sx_20 = zeros(1, Ng);
    % for n = 1:Ng
    %     i_s = find(Fp(:, n) > 0, 1);
    %     sx_20(n) = interp1(Fp(i_s:end, n), s(i_s:end), t0);
    % end
    % scatter(sx_20, repmat(t0, [1 Ng]), 'kx', linewidth=1.5)   
    legend([pl1, pl2, pl3, plN(1)], '$\sigma_p(0,s)$', '$\sigma_p(1,s)$', '$\sigma_p(2,s)$', '$\sigma_p(3..10,s)$', 'Location', 'Northwest', Interpreter='latex');
    
    set(gca, 'FontSize', 12)
    aspect = 1.5;
    
    ylim([0 inf]);xlim([0 inf])
    xlabel('$s$ ($\mu$m)', Interpreter='latex');ylabel('$\sigma_p(n,s)$ (kPa)', Interpreter='latex');
    % visualizing the unfolding rate = fig 1C
    % clf;hold on;
    axes('position', [0.55 0.2 0.35, 0.7]);hold on;
    
    plot(s, RU(:, 1), 'k-','linewidth', 2);
    plot(s, RU(:, 2), 'k-.','linewidth', 2);
    plot(s, RU(:, 3), 'k--','linewidth', 2);
    plot(repmat(s, [1 Ng-3]), RU(:, 4:end), 'k:','linewidth', 2);
    
    plot(repmat(s(end), [Ng, 1])', RU(end, :)', 'ks', 'linewidth', 2);
    for n = 1:3
        text(s(end) + 0.01, RU(end, n), ['$U_{' num2str(n)  '\rightarrow ' num2str(n+1) '}$'], 'Fontsize', 12, Interpreter='latex')
    end
    text(s(end) + 0.01, RU(end, 4), '...', 'FontSize',14, 'FontWeight','bold')
    xlabel('$s$ ($\mu$m)', Interpreter='latex');ylabel(['$U_{n\rightarrow n+1}(n, s)$  (s$^{-1}$)'], Interpreter='latex');
    gc = gca;
    set(gca, 'FontSize', 12);
    aspect = 1.5;
    % set(gcf, 'Position', [500  300  3.5*96 3.5*96/aspect])
    % gc.Position = gc.Position + [0 0 -0.1 0];
    set(gca, 'TickLength',[0.025 0.025]);
    
    % inset 
    axes('position', [0.62 0.55 0.21 0.35], 'YAxisLocation','right')
    plot(1:Ng, RU(end, :), 'ks-', 'linewidth', 2);
    xlabel('$n$', Interpreter='latex')
    text(4, max(RU(end, 1))*0.7, ["$U_{n\rightarrow n+1}$ at" sprintf("%0.3f $\\mu$m", max(s))], 'FontSize',14, 'FontWeight','normal', Interpreter='latex')
    set(gca, 'FontSize', 12);
    % set(gcf, 'Position', [500  240  400  300])
    % set(gca, 'YAxisLocation', 'right');
    set(gca, 'YTick',[0 200 400], 'XTick', [0 5 10], 'TickLength',[0.05 0.025]);
    exportgraphics(f,'../Figures/ModelRates2.png','Resolution',150)
    exportgraphics(f,'../Figures/ModelRates.eps','Resolution',150)

    % reconstruct Fp again without the NaN's
    % Fp = kp*(max(0,s-slack)/Lref).^(np); 
    % Fp(isnan(Fp() == NaN) = 0;
    return;
end
%% Folding rate design and visualization
% RF = alphaF*(max(0,delU-s))  % folding rates from state n+1 to n            
% RF = alphaF_0 + alphaF*(slack(1:Ng) - s).*(slack(1:Ng) > s);  % folding rates from state n+1 to n            
RF = alphaF_0;
% RF = 0.1 + alphaF*(slack(1:Ng) > s);  % folding rates from state n+1 to n            
% mesh(RF);view(3)

% clf;hold on;
% plot(repmat(s, [1 11]), RF(:, 1:11), 's-', 'linewidth', 2);
% xlabel('Strain (um)');ylabel('Transition rate');
% legend('U_0', 'U_{2\rightarrow1}', 'U_{3\rightarrow2}', '...')

% RF = 0;
%
% Initial state
PU = zeros(1,Ng+1); % initial unfolded probabilities for un-attached rectifier state
PA = zeros(1,Ng+1); % initial unfolded probabilities for attached rectifier state

pu = zeros(Nx,1)*PU;
pa = zeros(Nx,1)*PA;
pu(1,1) = 1/ds; 

if pCa >= pCa_noEffect
    % no Ca effect assumed
    x0 = reshape(pu,[(Ng+1)*Nx,1]);
else
    % might have some Ca effect
    x0 = reshape([pu, pa],[2*(Ng+1)*Nx,1]);
end
x0 = [x0; 0]; 

% Vlist = [1 10 100 1000 5000]*Lmax/100; %  half-sarcomere velocity
% reducing number of ranges
% Vlist = [10 100 1000]*Lmax/100; %  half-sarcomere velocity (um/s)
Vlist = Lmax./rds;
Force = cell(1, 5); 
Time = cell(1, 5); 
Length = cell(1, 5); 
for j = rampSet
  if isempty(datatables{j})
      % fprintf('Skipping pCa %0.2f %0.0fs dataset\n', pCa, rds(j))
      continue;
  end
  % tic
  V = Vlist(j); % ramp velocity

  Tend_ramp = Lmax/V; % length of ramp

  % testing tolerances, all with +/- same total cost
  % abstol 1e-6 7.7s,  1e-3 3.6s, 1e-1 2.7s and 1e1 3.1s
  % reltol 1e-3 (normal) 7.7s, 1e-1 6.8s and 8s for 1e-5 
    
  % test solvers - do not touch ode15s
  % ode15s 7.5s, ode45 53s, ode23s 173s, ode23t 9.7s

  % test sparse matrix - takes forever
  % S = [reshape(triu(ones(size(pu))),[(Ng+1)*Nx,1]); 1];
  % Sp = S*S';  % opts = odeset('JPattern',Sp');

  % [t0,x0] = ode15s(@dXdT,[-100:1:0],   x0,opts,Nx,Ng,ds,kA,kD,kd,Fp,RU,RF,mu,Lref,nd,0);
  % [t1,x1] = ode15s(@dXdT,[0 Tend_ramp],x0(end,:),opts,Nx,Ng,ds,kA,kD,kd,Fp,RU,RF,mu,Lref,nd,V);
  % 
  
  % whole decay till the bitter end
  % [t2,x2] = ode15s(@dXdT,[Tend_ramp 200],x1(end,:),opts,Nx,Ng,ds,kA,kD,kS,Fc,RU,RF,mu,Ls0,nS,0);
  % limited decay
  % [t2,x2] = ode15s(@dXdT,[Tend_ramp Tend_ramp + 40],x1(end,:),[],Nx,Ng,ds,kA,kD,kd,Fp,RU,RF,mu,Lref,nd,0);
  % only ramp up, no decay  
  % x2 = [];t2 = []; 

  % t = [t1(1:end); t2(2:end)];% prevent overlap at tend_ramp
  % x = [x1; x2(2:end, :)];

  
  if strcmp(simtype, 'ramp')
      %% normal
      times = [-100, 0;0 Tend_ramp;Tend_ramp Tend_ramp + 200];
      % times = [-100, 0;0 Tend_ramp];
      velocities = {0 V 0};
      L0 = 0;
  elseif strcmp(simtype, 'sin')
    %% sinusoidal driving
      % cycle time
      Tc = rds(j);
      times = [-100, 0;0 10*Tc];  
      positions = @(t)-Lmax/2*cos(2*pi*t/Tc) + Lmax/2;
      % differentiating positions
      V = @(t)2/2*Lmax*pi/Tc *sin(2*pi*t/Tc);
    
      % syms fpos(t);   % fpos(t) = @(t)Lmax*sin(2*pi*t/Tc);
    
      velocities = {0, V};
      L0 = 0;%Lmax/2;
  elseif strcmp(simtype, 'refolding')
    % provide times and velocities
  elseif strncmp(simtype, 'velocitytable_', 14)
    % simtype = 'velocitytable_relaxed.csv'
    vtb = readtable(['../data/' simtype ]);
    times = [vtb.Time(1:end-1),vtb.Time(2:end)];
    velocities = num2cell(vtb.Velocity);
    L0 = 0.05; % that is 0.95 + 0.05 = 1
  end
  %% repeated - refolding
  % times = [-100, 0;0 Tend_ramp;Tend_ramp Tend_ramp + 40;... % normal ramp-up
  %     Tend_ramp + 40 Tend_ramp + 41;... % rampdown in 1s
  %     Tend_ramp + 41 Tend_ramp + 41 + 60;... % hold down - wait for refold
  %     Tend_ramp + 41 + 60 Tend_ramp + 41 + 60 + Tend_ramp;... ramp-up
  %     Tend_ramp + 41 + 60 + Tend_ramp Tend_ramp + 41 + 60 + Tend_ramp + 40]; % final hold
  % velocities = [0 V 0 ...
  %     -Lmax/1 0 V 0];
  %% 10 heartbeats
  % times = [-100, 0];velocities = [0];
  % HR = 600;Tc = 60/HR;
  % sf = 0.2; % systole fraction
  % df = 0.4; % diastole fraction
  % Lmax = (2.2 - 1.9)/2;
  % Vc = Lmax/(Tc*sf);Vr = Lmax/(Tc*df);
  % for b = 1:10
  %   % relaxation - blowing up
  %   times = [times; times(end) times(end)+Lmax/Vr];
  %   velocities = [velocities Vr];
  %   % hold
  %   times = [times; times(end) times(end)+Tc*(1-sf - df)];
  %   velocities = [velocities 0];
  %   % contraction
  %   times = [times; times(end) times(end) + Lmax/Vc];
  %   velocities = [velocities -Vc];
  %   drawAllStates = true;
  % end
%%
  pu = zeros(Nx,1)*PU;
  pa = zeros(Nx,1)*PA;
  pu(1,1) = 1/ds; 
  if pCa >= pCa_noEffect
    x0 = reshape(pu,[(Ng+1)*Nx,1]);
  else
    x0 = reshape([pu, pa],[2*(Ng+1)*Nx,1]);
  end
  x0 = [x0; L0(1)]; 
  % opts = odeset('RelTol',1e-1, 'AbsTol',1e-1);          
  opts = odeset('RelTol',1e-2, 'AbsTol',1e-2);
  % assert(length(times) == length(velocities), 'Must be same length')
  t = []; x = [];
  for i_section = 1:size(times, 1)
      [t1,x1] = ode15s(@dXdT,times(i_section, :),x0,opts,Nx,Ng,ds,kA,kD,kd,Fp,RU,RF,mu,L_0,nd,kDf,velocities{i_section});
      t = [t; t1(2:end)];
      x = [x; x1(2:end, :)];
      % prep the init vector again
      x0 = x1(end,:);
      if length(L0) >=i_section
          x0(end) = L0(i_section);
      end
  end
  % validRng = t >= 0;
  % t = t(validRng);x = x(validRng, :);
  %% ramp downm, wait and up again
  % Vdown = Lmax./0.1;
  % [t3,x3] = ode15s(@dXdT, t2(end)+[0 0.1],x2(end,:),[],Nx,Ng,ds,kA,kD,kd,Fp,RU,RF,mu,Lref,nd,-Vdown);
  % [t4,x4] = ode15s(@dXdT, t3(end)+[0 100],x3(end,:),[],Nx,Ng,ds,kA,kD,kd,Fp,RU,RF,mu,Lref,nd,0);
  % [t5,x5] = ode15s(@dXdT, t4(end)+[0 rds(j)],x4(end,:),[],Nx,Ng,ds,kA,kD,kd,Fp,RU,RF,mu,Lref,nd,V);
  % [t6,x6] = ode15s(@dXdT, t5(end)+[0 10],x5(end,:),[],Nx,Ng,ds,kA,kD,kd,Fp,RU,RF,mu,Lref,nd,0);
  % 
  % 
  % t = [t1(1:end); t2(2:end); t3(2:end); t4(2:end); t5(2:end); t6(2:end)];% prevent overlap at tend_ramp
  % x = [x1; x2(2:end, :); x3(2:end, :); x4(2:end, :); x5(2:end, :); x6(2:end, :)];
%%
Time{j} = t;
states{j} = [];states_a{j} = [];    strains{j} = []; i_time_snaps = [];
maxPu = 0; maxPa = 0;
% save pca4statesenv
% load pca4statesenv
% drawAllStates = 1;
% save relaxstatesenv
% load relaxstatesenv
% drawAllStates = true;
    if drawAllStates == j
        % save current figure
        g = gcf;
        % open up a new one
        layout_x = 2 + (pCa < 10);
        if ~exist('statesFig')
            statesFig = figure(50+pCa*10+j); clf;
        else
            figure(statesFig);
        end
        
        tiledlayout(layout_x, 4, 'TileSpacing','compact', Padding='loose');
        fontsize(12, 'points')
        aspect = 2;
        % normal size of 2-col figure on page is 7.2 inches
        % matlab's pixel is 1/96 of an inch
        statesFig.Position = [300 200 7.2*96 7.2*96/aspect];

        % decide for timepoints
        % fixed time or fraction of ramp durations?
        % time_snaps = [0, 0.1, 1, 10, 30, 40, 100]
        % time_snaps = [0, rds(j), rds(j) + 30, 60, 120, 160];
        if strcmp(simtype, 'ramp')
            time_snaps = [1e-3, rds(j), 0*1 + 2*rds(j), max(100, rds(j)*10)];
        elseif strcmp(simtype, 'sin')
            time_snaps = [Tc/2, Tc, 3*Tc/2, 2*Tc];
        end
        % i_time_snaps = find(t > time_snaps)
    
        % disable
        % time_snaps = [];
        i_time_snaps = [];
        for i = 1:length(time_snaps)
            if time_snaps(i) > t
                break;
            end
            i_time_snaps(i) = find(t>=time_snaps(i), 1, 'first');    
        end
    end

  for i = 1:length(t)
    xi = x(i,:);
    Length{j}(i) = xi(end);
    pu = reshape( xi(1:(Ng+1)*Nx), [Nx,Ng+1]);
    if pCa >= pCa_noEffect
        pa = 0;
    else
        pa = reshape( xi((Ng+1)*Nx+1:2*(Ng+1)*Nx), [Nx,Ng+1]);
    end
    Fd = kd*max(0,(Length{j}(i) - s)/L_0).^nd; 
    Force_pa{j} = ds*sum(sum(Fd.*pa ));
    Force{j}(i) =  ds*sum(sum(Fd.*pu )) + Force_pa{j};
    states{j}(i, 1:Ng+1) = sum(pu);
    states_a{j}(i, 1:Ng+1) = sum(pa);
    strains{j}(i, 1:Nx) = sum(pu, 2);

    if drawAllStates == j
        if any(ismember(i_time_snaps, i))
            % sum of all states, converting to %
            ss = (sum(pu(:)) + sum(pa(:)))*ds;
            maxPu = max(maxPu, sum(pu(:))*ds*100);
            maxPa = max(maxPa, sum(pa(:))*ds*100);
            i_snap = find(i_time_snaps == i, 1, 'first');


            % next plot
            tl = nexttile(i_snap);            
            %% 3D plot
            % surf(0:Ng, s, pa*ds*100, 'EdgeColor', 'none');hold on;
            % % contour(0:Ng, s, pu*ds*100, 'EdgeColor','red');hold on;
            % 
            % % surf(s, 0:Ng, Fp')
            % colormap(1-gray);
            % % xlim([0, s(end)]);
            % % shading(gca, 'interp')
            % % view(90, -90); 
            % clim([0 maxPu]);
            % % xlabel('s'); ylabel('State');
            % % title(sprintf('U (%fs), S= %0.1f', t(i), sum(pu(:))));
            % if i_snap == 1
            %     xlabel('# State');
            % elseif i_snap == length(i_time_snaps)
            %     cb = colorbar;
            %     title(cb, '%')
            % end
            % ylabel('$s$ ($\mu$m)', Interpreter='latex');
            %% 2D shaded plot - unattached
            imagesc(s, 0:Ng, pu'*ds*100);hold on; box on;
            set(gca,'YDir','normal');
            mesh([s  - ds/2; s(end) + ds/2], (0:Ng+1) - 0.5, zeros(size(pu') + [1 1]), LineStyle='-', EdgeColor=[1 1 1]*0.9, FaceColor='none')
            contour(s, 0:Ng, pu'*ds*100, [0.01 0.01], 'k-', LineWidth=2)

            % surface(s, 0:Ng, pu'*ds*100, 'EdgeColor',[1 1 1]*0.9);hold on;box on;
            % % [x,y,z]=meshgrid(s,0:Ng,-1);
            % % surface(x, y, z)
            % h = 10;lw = 0.5;
            % line([s(1) s(1)], [0, Ng], [h h], Color='black', linewidth=lw)
            % line([s(end) s(end)], [0, Ng], [h h], Color='black', linewidth=lw)
            % line([s(1) s(end)], [0, 0], [h h], Color='black', linewidth=lw)
            % line([s(1) s(end)], [Ng, Ng], [h h], Color='black', linewidth=lw)
            % imagesc(s, 0:Ng, pu'*ds*100);
            % surf(s, 0:Ng, Fp')
            colormap(tl, SoHot);
            
            % xlim([0, s(end)]);
            % shading(gca, 'interp')
            % view(90, -90); 
            clim([0 round(maxPu)]);
            % xlabel('s'); ylabel('State');
            % title(sprintf('U (%fs), S= %0.1f', t(i), sum(pu(:))));
            if i_snap == 1
                ylabel('$n$', Interpreter='latex');
            elseif i_snap == length(i_time_snaps)
                cb = colorbar;
                cb.Ticks = unique([0 50 round(maxPu)]);
                title(cb, 'Probability of unattached $p_u$ (\%)',Interpreter='latex');
                cb.Title.HorizontalAlignment = 'right';
                set(tl, "YTick", [])
            else
                set(tl, "YTick", [])
            end
            if pCa >= 10
                % xlabel is done by attached panels
                xlabel('$s$ ($\mu$m)', Interpreter='latex');
            else
                ga = gca;
                ga.XTick = [];
            end 
            fontsize(12, 'points');
            set(tl, 'TickLabelInterpreter', 'latex');
            %% 2D shaded plot - attached
            if pCa < 10
                aspect = 1.4; % make the plots rectangular for pca as well
                statesFig.Position = [300 200 7.2*96 7.2*96/aspect];

                tl = nexttile(i_snap + 4);
                % surface(s, 0:Ng, pa'*ds*100, 'EdgeColor',[1 1 1]*0.9);hold on;box on;                
                
                imagesc(s, 0:Ng, pa'*ds*100);hold on; box on;
                set(gca,'YDir','normal');
                mesh([s  - ds/2; s(end) + ds/2], (0:Ng+1) - 0.5, zeros(size(pa') + [1 1]), LineStyle='-', EdgeColor=[1 1 1]*0.9, FaceColor='none')
                contour(s, 0:Ng, pa'*ds*100, [0.01 0.01], 'k-', LineWidth=2)

                % line([s(1) s(1)], [0, Ng], [h h], Color='black', linewidth=lw)
                % line([s(end) s(end)], [0, Ng], [h h], Color='black', linewidth=lw)
                % line([s(1) s(end)], [0, 0], [h h], Color='black', linewidth=lw)
                % line([s(1) s(end)], [Ng, Ng], [h h], Color='black', linewidth=lw)

                % surf(s, 0:Ng, Fp')
                % colormap(1-gray);
                % colormap(SoHot);
                colormap(tl, SoCool);
                % xlim([0, s(end)]);
                % shading(gca, 'interp')
                % view(90, -90); 
                clim([0 round(maxPa, 1)]);
                % xlabel('s'); ylabel('State');
                % title(sprintf('U (%fs), S= %0.1f', t(i), sum(pu(:))));
                if i_snap == 1
                    ylabel('$n$', Interpreter='latex')
                elseif i_snap == length(i_time_snaps)
                    cb = colorbar;
                    cb.Ticks = [0 round(maxPa, 1)];
                    title(cb, 'Probability of attached $p_a$ (\%)',Interpreter='latex')
                    cb.Title.HorizontalAlignment = 'right';
                    set(tl, "YTick", [])
                else
                    set(tl, "YTick", [])
                end
                xlabel('$s$ ($\mu$m)', Interpreter='latex');
                box on;
                fontsize(12, 'points');
                set(tl, 'TickLabelInterpreter', 'latex');
            end
            %% lollipop plot
            % % identify leading states
            % lead_i = sum(pu(:, 1:Ng+1)) > 10;
            % % Disabling unimportant states
            % pup = pu;pup(pup < 1e-1) = NaN;
            % pap = pa;pap(pap < 1e-1) = NaN;
            % 
            % phpa = stem(repmat(s, [1 sum(lead_i)]), pup(:, lead_i)/ss, ':', 'LineWidth', 1.5, MarkerSize=7);
            % if length(pa) > 1
            %     hold on;
            %     ca = gca;
            %     ca.ColorOrderIndex = 1;
            %     phpu = stem(repmat(s, [1 sum(lead_i)]), pap(:, lead_i)*10/ss, ':', 'filled', LineWidth=1.5, MarkerSize=5);
            %     if true || i_snap == 1
            %     leg = legend([phpa,phpu], [strcat("$S_U$ ", string(find(lead_i == 1))), ...
            %         strcat("$S_A$ ", string(find(lead_i == 1)), " (x10)")],'Location', 'best', NumColumns=2, Interpreter='latex');                
            %     else
            %         leg = legend(phpa, strcat("State ", string(find(lead_i == 1))), 'Location', 'best');
            %     end
            % else
            %     leg = legend(phpa, strcat("State ", string(find(lead_i == 1))), 'Location', 'best');
            % end
            % xlabel('$s$ ($\mu$m)', Interpreter='latex');          
            % if nt.Position(1) < 0.4
            %     % left aligned
            %     ylabel('$P_S$ (\%)', Interpreter='latex')
            % end
            % title(leg, sprintf('t = %g (%g)', round(t(i), 3), ss));
            % leg.ItemTokenSize = [10 10];
            % axis([0, s(end), 0, ceil(max(pu(:))/10)*10]);
            %%
            
            
            

            
            % ignore pa for a moment
            % if ~(pa == 0)
            %     subplot(2, length(i_time_snaps), i_snap + length(i_time_snaps));
            %     surf(s, 0:Ng -1, pa, 'EdgeColor','none');
            %     view(90, -90); xlabel('s'); ylabel('State');
            %     title(sprintf('A (%fs)', t(i)))
            % end
        end
    end    

    % if t(i) >= rds(j) && t(max(1, i-1)) < rds(j)
    %     % at the peak
    %     outStruct{j, 1}.pa = pa;
    %     outStruct{j, 1}.pu = pu;
    % elseif t(i) >= rds(j)*2 && t(max(1, i-1)) < rds(j)*2
    %     % after t_rd after the peak
    %     outStruct{j, 2}.pa = pa;
    %     outStruct{j, 2}.pu = pu;
    % elseif t(i) >= rds(j) + 20 && t(max(1, i-1)) < rds(j) + 20
    %     % after 20s after the peak
    %     outStruct{j, 3}.pa = pa;
    %     outStruct{j, 3}.pu = pu;
    % end
    end
% Continuous states in time: time vs. state occupancy plot
% if drawAllStates
%     figure(60+j); clf;
%     set(gcf, 'Name', sprintf('States for pCa %d at %0.1f ramp', pCa, rds(j)) );
%     subplot(131);
%     h = surf(0:Ng, Time{j}, states{j}, 'EdgeColor','none');
%     shading(gca, 'interp')
%     view(90, -90)
%     ylabel('Time (s)'); xlabel('State occupancy (#)')
%     title('States pu');
%     colorbar;
% 
%     subplot(132);
%     plot(t, Force{j}, LineWidth=2);
%     subplot(132);
%     surf(0:Ng, Time{j}, states_a{j});
%     shading(gca, 'interp')
%     view(90, -90)
%     ylabel('Time (s)'); xlabel('State occupancy (#)')
%     title('States pa');
%     colorbar;
% 
%     subplot(133);
%     surf((1:Nx)*ds, Time{j}, strains{j});
%     shading(gca, 'interp')
%     view(90, -90)
%     ylabel('Time (s)'); xlabel('Strains (um)')
%     title('Strains');
% 
%     colorbar;
% 
%     % set to preset figure
%     figure(g);
% end

%

    
    % show state occupation at the end of the ramp
    % clf;
    % i = round(size(x2, 1)*3/4);
    % xi = x2(i,:);    
    % pu = reshape( xi(1:(Ng+1)*Nx), [Nx,Ng+1]);
    % surf(pu)
    % xlabel('Ng');ylabel('Nx');zlabel('State probability');    
    % hold on;
    % title(['Ramp ' num2str(Tend_ramp) 's, V = ' num2str(V) ' ML/s at t ' num2str(t2(i))]);


  % add parallel static force
  % Force{j} = Force{j} + mod(10)*3;
  % Force{j} = Force{j} + 130*(Length{j}).^4;
  % Force{j} = Force{j} + 0.0045*(1.6 + 2*Length{j}).^7;

  % decay offset is set, now use a nonlinear func to fit the ramp onset
  % a*(-b + Lmax).^c + d = 1.2716*3;
  % a*(-b + Lmax).^c = 1.2716*3 - d;
  
  
  % apply constraints
  if n_ss < 0 
      cost = inf;
      return;
  end
  
  % calculate a, so that the max value is the same    
  k_ss = (Fss)/((Lmax)^n_ss);  
  Force_par{j} = k_ss*max(Length{j}, 0).^n_ss;
  % calc force
  Force{j} = Force{j} + Force_par{j}; 
    
  if drawAllStates == j
        nexttile([1 4]);
        semilogx(Time{j}, Force{j}, 'k-');hold on;
        scatter(Time{j}(i_time_snaps), Force{j}(i_time_snaps), 'ko', 'filled');
        xlim([1e-3 max(1300, rds(j)*10) + 30]);
        xlim([1e-3 200]);
        % nexttile([1 1]);
        % loglog(Time{j}, Force{j});hold on;
        % scatter(Time{j}(i_time_snaps), Force{j}(i_time_snaps), 'o', 'filled');
        % xlim([1e-3 inf])
        xlabel('$t$ (s)', Interpreter='latex');
        legend('Tension', 'Insets', 'Location','northwest');
        ylabel('$\Theta$ (kPa)', Interpreter='latex')
        fontsize(12, 'points');
        set(tl, 'TickLabelInterpreter', 'latex');
        % exportgraphics(f,sprintf('../Figures/States%g_%gs.png', pCa, rds(j)),'Resolution',150)
    end

  if any(isnan(Force{j}))
      disp('error');
      cost = inf;
      return;
  end

  % Force{j} = Force{j} + (0.55e6)*0.225^8*mod(10); 


  % Fss = mod(10)*3; % optimized previously
  % Force{j} = Force{j} + Fss;

    % ttoc = toc;
  % fprintf("Ramp %0.1fs takes %1.1fs \n", Tend_ramp, ttoc)
  %% plot the beating
  % g = gcf;
  % figure(2000+pCa*10);clf;
  % plot(Time{j}, Force{j}, Time{j}, Force{j} - Force_par{j});
  % legend('Tension total', 'Tension viscous (kPa)');
  % title(sprintf('pCa %0.1f, HR %0.0f bpm', pCa, HR));
  % figure(g);

  if strcmp(simtype, 'sin')  
  %% Plot sinusoidal outcome
    aspect = 2;
    % figure(900 + j*10 + round(pCa));clf;    
    clf;tiledlayout(2,2, TileSpacing="compact");
    % figure(900935);
    % clf; tiledlayout('flow');
    set(gcf, 'Position', [500  300  7.2*96 7.2*96/aspect])
    rng = Time{j} < 20;
    t = Time{j}(rng)';
    x = Length{j}(rng) + 0.95;
    y = Force{j}(rng);
    y2 = Force{j} - Force_par{j};
    z = zeros(size(t));
    col = linspace(0, t(end), length(t)); 
    lw = 1;

    % if exist('tl1', 'var')
        nexttile(1);hold on;
    % else
        % tl1 = nexttile;
    % end
    plot(t, x, 'k', linewidth = 1);
    % plot(t, x, t, y);legend('Length', 'Force');
    xlabel('$t$ (s)', Interpreter='latex');    ylabel('$L$ ($\mu$m)', Interpreter='latex');
    ylim([0.95, 1.2])

    % if exist('tl2', 'var')
        % nexttile(2);hold on;
    % else
        tl2 = nexttile(2, [2 1]);
    % end

    surface([x;x],[y;y],[z;z],[col;col],...
            'facecol','no',...
            'edgecol','interp',...
            'linew',lw);
    % clim([0 50 100])
    colormap(turbo);    
    ylabel('$\Theta$ (kPa)', Interpreter='latex');     xlabel('$L$ ($\mu$m)', Interpreter='latex');
    cb = colorbar;  title(cb, 't (s)');

    % if exist('tl3', 'var')
        nexttile(3);
    % else    
        % tl3 = nexttile;
    % end
    % plot(x_(t), x, x_(t), y);legend('Length', 'Force');
    plot(t, y, 'k', linewidth = lw);
    xlabel('$t$ (s)', Interpreter='latex');
    ylabel('$\Theta$ (kPa)', Interpreter='latex');
    fontsize(12, 'points');
    % exportgraphics(gcf,sprintf('../Figures/FigHysteresis%g_%gs.png', pCa, rds(j)),'Resolution',150)
  elseif strncmp(simtype, 'velocitytable_', 14)
      
      %%
      set(groot,'CurrentFigure',figInd); % replace figure(indFig) without stealing the focus
      cf = clf;
      F_data{1} = interp1(datatables{1}.Time, datatables{1}.F, Time{1}); % total force interpolated
      nexttile;
      L_data{1} = interp1(datatables{1}.Time, datatables{1}.L, Time{1}); % total force interpolated
      plot(Time{j}, Length{j} + 0.95, Time{j}, L_data{1});
      
      nexttile;
      plot(Time{1}, F_data{1}, '-',Time{1}, Force{1}, '--',lineWidth = 2 );
      % plot(Time{j}, Force{j}, datatables{1}.Time, datatables{1}.F, Time{1}, F_data{1});

      err =  nansum((F_data{1} - Force{1}').^2);

      % findpeaks(datatables{1}.F, datatables{1}.Time, MinPeakDistance=120);
        [py, px] = findpeaks(datatables{1}.F, datatables{1}.Time, MinPeakDistance=120);
        % hold on;
        sel = [2 3 4 5 6 8 9 11];
        
        % plot(px(sel), py(sel), 's', LineWidth=4);
        
        [pm] = findpeaks(Force{1}, Time{1}, MinPeakDistance=140);
        
        if length(pm) == length(sel)
            peak_err = sum(((py(sel) - pm').^2))*3e3;
        else
            peak_err = 1e5;
        end

      cost = err + peak_err;
  end
end


if strcmp(simtype, 'sin')
    % for sin type its the end
    return
end
% Get error for the whole ramp-up and decay
t_endFreeware = zeros(1, 5); % time when we start counting the costs
% alternatively, get the error from decay only
% t_endFreeware  = Lmax./Vlist + 2;

%% Evaluating all ramps at once
En = cell(1, length(rampSet));
Es = cell(0);
for j = rampSet
    datatable_cur = datatables{j};
    if isempty(datatable_cur)
        continue;
    end
    inds = find(datatable_cur.Time >= t_endFreeware(j));
    datatable_cur = datatable_cur(inds, :);
    t_int{j} = datatable_cur.Time - 2;
    % t_int{j} = datatable_cur.Time;
    Ftot_int{j} = interp1(Time{j}, Force{j}, t_int{j}); % total force interpolated

    % weighting to fit semilogx
    x = -log10(max(rds(j), t_int{j}));
    % center and scale
    w = (x+abs(min(x)))/mean(x+abs(min(x)));
    
    % semilogx(t_int{j}, x/sum(x));
    % plot(t_int{j}, w);hold on;

    % no weighing, already in the data
    w = 1;
%%
    Es{j} = w.*((Ftot_int{j} - datatable_cur.F)/max(Ftot_int{j})).^2; % error set
    Es{j} = Es{j}(~isnan(Es{j})); % zero outside bounds
    En{j} = 1e3*sum(Es{j})/length(Es{j}); % normalized error
    
    PeakData(j, 1) = rds(j);
    PeakData(j, 2) = max(datatable_cur.F);
    PeakModel(j) = max(Ftot_int{j});

end

PeakModel = nan(1, length(PeakData));
for j = 1:length(PeakModel)
    m = max(Force{j});
    if ~isempty(m)
        PeakModel(j) = m;
    end
end

Ep = nansum(((PeakData(:, 2) - PeakModel(:))/max(Ftot_int{j})).^2);
% discarding peak fit for high Ca's
if pCa < 10
    % Ep = 0;
end

cost = Ep*100 + sum([En{1:end}], 'all');

if exist('drawPlots', 'var') && ~drawPlots
    return;
end

Tarr = t_int;Farr = Ftot_int;
if exportRun
    % save data for decay overlay loglog plot
    % Tarr = t_int;Farr = Ftot_int(1:4);
    % if pCa < 10
    %     save pcagraphingenv.mat
    % else
    %     save relaxgraphingenv.mat
    % end
    save(sprintf('..\\pca%gmodeldata.mat', pCa), 'Tarr', 'Farr')
    % save('..\pca11modeldataDoubleStates.mat', 'Tarr', 'Farr')
end

% tic
aspect = 2;
if plotInSeparateFigure
    set(groot,'CurrentFigure',figInd); % replace figure(indFig) without stealing the focus
    cf = clf;
    % cf = gcf;
     
    % normal size of 2-col figure on page is 7.2 inches
    % matlab's pixel is 1/96 of an inch
    f.Position = [300 200 7.2*96 7.2*96/aspect];
    
    % f.Position = [300 200 7.2*96 7.2*96/aspect];
    % colors = lines(max(rampSet)+1); colors(1:end-1, :) = colors(2:end, :);
    
    % colors(1:end-1, :) = colors(2:end, :);
    fs = 12;
    % tl = tiledlayout(3, 4, 'TileSpacing', 'compact');
    rws_s = 4;rws_l = rws_s*4;rws_loglog = 6;
    tl = tiledlayout(3, rws_loglog+rws_l, "Padding","compact", "TileSpacing","compact");
    tile_semilogx = nexttile(1, [2, rws_l]);
    tile_loglog = nexttile(rws_l + 1, [3, rws_loglog]);
    tile_r(1) = nexttile((rws_loglog+rws_l)*2 + 1, [1 rws_s]);
    tile_r(2) = nexttile((rws_loglog+rws_l)*2 + 1 + rws_s, [1 rws_s]);
    tile_r(3) = nexttile((rws_loglog+rws_l)*2 + 1 + 2*rws_s, [1 rws_s]);
    tile_r(4) = nexttile((rws_loglog+rws_l)*2 + 1 + 3*rws_s, [1 rws_s]);
    tile_positions = [tile_semilogx.Position;tile_loglog.Position;...
        tile_r(1).Position;tile_r(2).Position;tile_r(3).Position;tile_r(4).Position];
    clf;
    
    reportCosts = false;
    if reportCosts 
        title(tl, sprintf('pCa %g costs %g', pCa, round(cost, 3)));    
    end
    % tile_semilogx = nexttile(1, [2, rws_l]);hold on;
    % tile_semilogx = subplot(3, rws_l + rws_loglog, );hold on;
    tile_semilogx = axes('Position', tile_positions(1, :));hold on;
else
    f.Position = [300 200 7.2*96 7.2*96/aspect];
    tile_semilogx = gca();
    hold on;
end
ym = 0;
for j = max(rampSet):-1:1
    if isempty(Force{j})
        continue;
    end
    if ~exist('compareFig', 'var') || ~compareFig
        hd{j} = errorbar(datatables{j}.Time-2,datatables{j}.F,datatables{j}.SD, '-', LineWidth=2, Color=colors(3, :), CapSize=0);
        set([hd{j}.Bar, hd{j}.Line], 'ColorType', 'truecoloralpha', 'ColorData', [hd{j}.Line.ColorData(1:3); 255*0.4])
        hm{j} = semilogx(Time{j},Force{j},'-', 'linewidth',2.5, 'Color', colors(1, :)); 
        
        % if ~exist('ignoreLegends', 'var') || ~ignoreLegends
        hl = legend([hd{max(rampSet)} hm{max(rampSet)}], ...
        'Data',...    
        'Model',...
        NumColumns=3);
        % end
    
        hl.ItemTokenSize = [30, 20];
        hl.Box = 'off';

    else
        hm{j} = semilogx(Time{j},Force{j},'r-', 'linewidth',1.5); 
    end
    
    % hph{j} = plot(nan, nan, 'x',Color=[1 1 1]); % just a placeholder
    % set(h.Cap, 'EdgeColorType', 'truecoloralpha', 'EdgeColorData', [h.Cap.EdgeColorData(1:3); 255*alpha])
    ym = max([ym Force{j}, datatables{j}.F']);
    set(tile_semilogx, 'TickLength', [0.0125 0.05]);
    set(tile_semilogx, 'TickLabelInterpreter', 'latex');
    ylabel('$\Theta$ (kPa)', Interpreter='latex')
    
    switch (rds(j))
        case 0.1
            % fastest
            tp = [rds(j)*0.8, max(Force{j})*0.92];
        case {1 10 100}
            tp = [rds(j)*1.3, max(Force{j})*1.1 + 0.5];
    end
        
    if plotInSeparateFigure
        % too small otherwise
        text(tp(1), tp(2), ...
        sprintf('$t_r$ = %g s', rds(j)),...
        'horizontalAlignment', 'right', VerticalAlignment='bottom', Interpreter='latex');    
    end
end
xlim([1e-2, 160])
ylim([0 ceil((ym)/5)*5])
yl = ylim();
tile_semilogx.XScale='log';
box on;
% hl = legend([hph{1} hph{2} hph{3} hph{4} hd{4} hd{3} hd{2} hd{1} hm{4} hm{3} hm{2} hm{1}], ...
%     't_r 0.1 s:', '  t_r 1 s:',' t_r 10 s:','t_r 100 s:',...
%     'Data', 'Data', 'Data', 'Data', ...    
%     'Model', 'Model', 'Model', 'Model',...
%     NumColumns=3);

if ~plotDetailedPlots || ~plotInSeparateFigure
    return;
end

% title('Model fit', Interpreter='latex')

    % 'MData Ramp-up 0.1 s', 'MData Ramp-up 0.1 s', 'MData Ramp-up 0.1 s', 'MData Ramp-up 0.1 s'...
    % )
    %%

% tile_semilogx.Position = tile_positions(1, :).*[1 1 1 1]
%% nexttile;
% semilogx(PeakData(:, 1), PeakData(:, 2), 'ko', LineWidth=2);hold on;
% semilogx(PeakData(:, 1), PeakModel, 'x', 'MarkerEdgeColor', [1 1 1]*0.5, LineWidth=2, MarkerSize=8);
% axis([1e-1 1e2 0 ym])
% semilogx(PeakData(:, 1), PeakModel, '--', Color=[1 1 1]*0.5, LineWidth=1);
% hl = legend('Data', 'Model', 'Location', 'best')
% title(hl, 'Peaks')
% ylabel('Tension (kPa)', Interpreter='latex')
%% loglog plot
% tile_loglog = nexttile(rws_l + 1, [3, rws_loglog]);
tile_loglog = axes('Position', tile_positions(2, :).*[1.05 1.8 1 1] + [0 0 0 -0.0440] );box on;
% tile_loglog.Position = tile_positions(2, :).*[1.05 1.8 1 1] + [0 0 0 -0.0440] 
% toc
% return;
%% best fit from FigFitDecayOverlay
 options = optimset('Display','iter', 'TolFun', 1e-4, 'Algorithm','sqp', 'UseParallel', true, ...
        'TolX', 0.0001, 'PlotFcns', @optimplotfval, 'MaxIter', 150);
 addpath ../DataProcessing

 if pCa > 10
    % x = [4.7976    0.2392    4.8212];
    x = [3.7242    0.2039    4.8357]; % data
    x = [3.4642    0.2413    5.0916]; % model refit
    x = [3.7364    0.2187    4.8357]; % model refit using the same Theta_inf
    [c rspca] = evalPowerFit(x, Force, Time, 'loglogOnly', [], false);
    
    if rerunFitting && c > 3
        %%
        disp('Rerunning power law fit...');
        % init = x;
        % fitfunOpt = @(x) evalPowerFit(x, Farr, Tarr, false);
        % x = fminsearch(fitfunOpt, init, options)
        % [c rspca] = evalPowerFit(x, Force, Time, 'loglogOnly', [], false);
        %% keep the Theta_inf
        init = x([1, 2]);
        fitfunOpt = @(p) evalPowerFit([p x(3)], Farr, Tarr, false);
        p = fminsearch(fitfunOpt, init, options)
        [c rspca] = evalPowerFit([p x(3)], Force, Time, 'loglogOnly', [], false);
    end

else
    % best overall fit
    % x = [15.1141    0.4346    4.5143];
    Phi_inf = 5.0506; % assumed from relaxed
    x = [6.4651    0.2383  Phi_inf -30];
    Phi_inf = 4.8357; % assumed from data
    x = [6.3528    0.2136  Phi_inf -30];
    % PEVK KO fit close enough
    % x = [6.1750    0.2272    5.0506  -39.1079];
    [c rspca] = evalPowerFit(x, Farr, Tarr, 'loglogOnly', [], true);

    if rerunFitting && c > 0.04
        %%
        init = x([1 2]);
        
        fitfunOpt = @(opt) evalPowerFit([opt(1) opt(2) x(3) x(4)], Farr, Tarr, false);
        opt = fminsearch(fitfunOpt, init, options)
        x([1 2]) = opt;
        [c rspca] = evalPowerFit(x, Farr, Tarr, 'loglogOnly', [], true);
        %% keep the theta_inf

    end
end


%%
for j = max(rampSet):-1:1
    if isempty(Force{j})
        continue;
    end
    % sp = nexttile((rws_loglog+rws_l)*2 + 1 + (4-j)*rws_s, [1 rws_s]);
    sp = axes('Position', tile_positions(7-j, :).*[1 1.8 1.2 0.8] + (4-j).*[-0.008 0 0 0]);box on;
    hold on;
    if reportCosts 
        title(sp, sprintf('Ramp %g costs %g', rds(j), round(En{j}, 3)));    
    end
    h = errorbar(datatables{j}.Time-2,datatables{j}.F,datatables{j}.SD, '-', LineWidth=2, Color=colors(3, :), CapSize=0);
    set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h.Line.ColorData(1:3); 255])
    plot(Time{j},Force{j},'-', 'linewidth',2, 'Color', 'k'); 
    xlim([0 min(rds(j)*3, 160)]);
    xl = xlim();
    ylim(yl);
    % ylim([0, inf]);
    % yl = ylim();
    % yl(2) = ceil(yl(2)/5)*5;
    % yticks([0 20 40]);
    xticks([0 rds(j), rds(j)*2])
    if j == 4
        ylabel('$\Theta$ (kPa)', Interpreter='latex')
    else
        yticklabels([]);
    end
    xlabel('$t$ (s)', Interpreter='latex', HorizontalAlignment='left');
    set(sp, 'TickLength', [0.05 0.05]);
    set(sp, 'TickLabelInterpreter', 'latex')
    switch (rds(j))
        case 0.1
            % fastest at the bottom
    tit = text(xl(2), 3, sprintf('$t_r$ = %g ', rds(j)), 'Interpreter', 'latex', ...
        HorizontalAlignment='right', VerticalAlignment='bottom');
        otherwise
    tit = text(xl(2), yl(2), sprintf('$t_r$ = %g ', rds(j)), 'Interpreter', 'latex', ...
        HorizontalAlignment='right', VerticalAlignment='top');
    end
    % tit.Position = tit.Position + [0 max(Force{j}) 0];
    % pos = sp.Position;
    % sp.Position = pos + [0 0 0.1 0];
    fontsize(12, 'points');
end    


aspect = 1.5;
set(cf, 'Position', [500  300  7.2*96 7.2*96/aspect]);
if exportRun
    exportgraphics(cf,sprintf('../Figures/ModelFitpCa%g.png', pCa),'Resolution',150)
    exportgraphics(cf,sprintf('../Figures/ModelFitpCa%g.eps', pCa))
end
% toc
