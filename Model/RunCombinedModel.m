% Assuming mod = ones(10,1)
% assuming pCa = Inf
load SoHot.mat;
% simtype = 'sin';
simtype = 'ramp';
% simtype = 'rampbeat';

% rampSet = [1];
% rampSet = [2 4];
% rampSet = [3]; % nly 100ms
rampSet = [1 2 3 4];
% rampSet = [4];

% pCa 11 - load Relaxed, do not run the extended, PEVK attachment model
% pCa 10 - load Relaxed, run the PEVK attachment model
% pCa < 10 - load AvgpCa dataset, Ca effect in place
clear Force
clear Time
clear Length
clear outStruct;
if any(mod < 0) 
    cost = inf;
    return;
end
drawAllStates = false;
drawFig1 = false;
exportRun = false;

figInd = get(groot,'CurrentFigure'); % replace figure(indFig) later without stealing the focus

% rds = fliplr([0.02 0.1, 1, 10 100]);
rds = fliplr([0.1, 1, 10, 100]);
% rds = fliplr([0.1, 10]);
for i_rd = 1:length(rds)
  if isinf(pCa) || pCa >= 10
  %   % hack - the no-Ca noPNB experiments had higher ramps
  %   datatable = readtable(['..\Data\bakers_passiveStretch_' num2str(rds(i_rd)*1000) 'ms.csv']);
  %   datatable.Properties.VariableNames = {'Time'  'ML'  'F'  'SL'};
  %   datatables{i_rd} = datatable;
  % elseif isnan(pCa)
      % newest format of experiments    
    datatables{i_rd} = readtable(['..\Data\AvgRelaxed_' num2str(rds(i_rd)) 's.csv']);
  else
    if exist(['..\Data\AvgpCa' num2str(pCa) '_' num2str(rds(i_rd)) 's.csv'], "file")
        datatables{i_rd} = readtable(['..\Data\AvgpCa' num2str(pCa) '_' num2str(rds(i_rd)) 's.csv']);
    else
        datatables{i_rd} = [];
    end
  end
  % else
  %   % new format for pCa experiments
  %   datatables{i_rd} = readtable(['..\Data\PassiveCa_2\bakers_passiveStretch_pCa' num2str(pCa) '_' num2str(1000*rds(i_rd)) 'ms.csv']);
  % end
end
% if isinf(pCa)
%     % hack - the no-Ca noPNB experiments had higher ramps
%   Lmax = 0.4;
% else    
    % follow-up experiments had lower ramps
    Lmax = 1.175 - 0.95;
% end

% Ls0  = 0.10*mod(13);
% Nx   = 25;          % number of space steps
% ds   = (0.36-Ls0)/(Nx-1);      % space step size
% s  = (0:1:Nx-1)'.*ds; % strain vector
% Ng  = 20;            % number of glubules on globular chain
% delU = 0.0125*mod(1);

% half-sarcomere ramp height
% Lmax = 0.225;
% Nx   = 25;          % number of space steps
Nx   = 25;          % number of space steps
ds   = 1*(Lmax)/(Nx-1);      % space step size
% use this to make sure the ds is big wnough
% ds   = 1*(Lmax + 0.015)/(Nx-1);
s  = (0:1:Nx-1)'.*ds; % strain vector
Ng = 15; 
% so all unfolded make mod(19) = Ng*delU slack, i.e. ~0.15 um. individual delU is about 0.01
delU = mod(19)/Ng;



% propose a function to kA = f(pCa)
% kA   = 1*mod(14);
% kD   = 1*mod(15);
% kC   = 103.33*mod(2);
% kS   = 300*mod(3);         % series element spring constant
% alphaU = 2000*mod(4);       % chain unfolding rate constant
% alphaF = 1*mod(5);
% nC = 1.77*mod(6);
% nS = 2.56*mod(7);
% nU = 4*mod(8);
% mu = 2.44*mod(9); 

% g0 = ones(1,11);
% g0 = mod;
% mod = 0;

kp   = mod(1); % ok      % proximal chain force constant
kd   = mod(2); % ok       % distal chain force constant

if ~isnan(mod(16))
    kA   = mod(16); % PEVK attachment rate
else
    kA   = mod(7);
end

if ~isnan(mod(17))
    kD   = mod(17); % PEVK detachment rate
else
    kD   = mod(8); % PEVK detachment rate
end

alphaU = mod(6);         % chain unfolding rate constant
% Fss = 3.2470*mod(10); % Parallel steady state force
% Fss = 4.89;
Fss = mod(10);
if pCa < 10
    kp   = mod(9);      % proximal chain force constantkS   = g0(2)*14122;        % distal chain force constant
    kA   = mod(7);
    kD   = mod(8); % PEVK detachment rate
    if ~isnan(mod(20))
        % cant exceed the no Ca unfolding rate!
        alphaU = min(alphaU, mod(20));         % chain unfolding rate constant
    end
    if ~isnan(mod(21))
        Fss = mod(21);
    end
    if ~isnan(mod(22))
        kd   = mod(22);      % proximal chain force constant high Cabist
    else
        kd = mod(2)*mod(9)/mod(1); % scaled to kp_ca in the same ratio as kd_relaxed to kp_relaxed
    end
        
end
kDf = mod(23);
alphaF = 0; % chain folding rate constant - not implemented yet
np = mod(3); % proximal chain force exponent
nd = mod(5); % distal chain force exponent
nU = mod(4); % unfolding rate exponent
nF = 1; % folding rate exponent (not implemented yet)
mu = mod(14); % small enough not to affect the result
L_0  = mod(18); % reference sarcomere length (um)
alphaF = 0;
alphaF_0 = mod(15);

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
    f = figure(3); clf;
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
    legend([pl1, pl2, pl3, plN(1)], '$\sigma_p(0,s)$', '$\sigma_p(1,s)$', '$\sigma_p(2,s)$', '$\sigma_p(3..14,s)$', 'Location', 'Northwest', Interpreter='latex');
    
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

if pCa >= 11 
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
      fprintf('Skipping pCa %0.2f %0.0fs dataset\n', pCa, rds(j))
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
      times = [-100, 0;0 Tend_ramp;Tend_ramp Tend_ramp + 10000];
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
      x_ = 0:Tc/100:Tc;
      plot(x_, positions(x_), '-', x_, V(x_), x_, cumsum(V(x_)).*[1 diff(x_)], '--');
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
  if pCa >= 11
    x0 = reshape(pu,[(Ng+1)*Nx,1]);
  else
    x0 = reshape([pu, pa],[2*(Ng+1)*Nx,1]);
  end
  x0 = [x0; L0]; 
  opts = odeset('RelTol',1e-3, 'AbsTol',1e-2);          
  % assert(length(times) == length(velocities), 'Must be same length')
  t = []; x = [];
  for i_section = 1:size(times, 1)
      [t1,x1] = ode15s(@dXdT,times(i_section, :),x0,opts,Nx,Ng,ds,kA,kD,kd,Fp,RU,RF,mu,L_0,nd,kDf,velocities{i_section});
      t = [t; t1(2:end)];
      x = [x; x1(2:end, :)];
      % prep the init vector again
      x0 = x1(end,:);
  end
  validRng = t >= 0;
  t = t(validRng);x = x(validRng, :);

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
    if drawAllStates
        % save current figure
        g = gcf;
        % open up a new one
        layout_x = 2 + (pCa < 10);
        f = figure(50+pCa*10+j); clf; tiledlayout(layout_x, 4, 'TileSpacing','compact', Padding='loose');
        fontsize(12, 'points')
        aspect = 2;
        % normal size of 2-col figure on page is 7.2 inches
        % matlab's pixel is 1/96 of an inch
        f.Position = [300 200 7.2*96 7.2*96/aspect];

        % decide for timepoints
        % fixed time or fraction of ramp durations?
        % time_snaps = [0, 0.1, 1, 10, 30, 40, 100]
        % time_snaps = [0, rds(j), rds(j) + 30, 60, 120, 160];
        if strcmp(simtype, 'ramp')
            time_snaps = [1e-3, rds(j), 0*1 + 2*rds(j), max(1000, rds(j)*10)];
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
    if pCa >= 11
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

    if drawAllStates
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
            colormap(SoHot);
            
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
                cb.Ticks = [0 50 round(maxPu)];
                title(cb, 'Unatt %')
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
            fontsize(12, 'points')
            %% 2D shaded plot - attached
            if pCa < 10
                aspect = 1.4; % make the plots rectangular for pca as well
                f.Position = [300 200 7.2*96 7.2*96/aspect];

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
                colormap(SoHot);
                % xlim([0, s(end)]);
                % shading(gca, 'interp')
                % view(90, -90); 
                clim([0 round(maxPa)]);
                % xlabel('s'); ylabel('State');
                % title(sprintf('U (%fs), S= %0.1f', t(i), sum(pu(:))));
                if i_snap == 1
                    ylabel('$n$', Interpreter='latex')
                elseif i_snap == length(i_time_snaps)
                    cb = colorbar;
                    cb.Ticks = [0 round(maxPa)];
                    title(cb, 'Att %')
                    set(tl, "YTick", [])
                else
                    set(tl, "YTick", [])
                end
                xlabel('$s$ ($\mu$m)', Interpreter='latex');
                box on;
                fontsize(12, 'points')
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
  
  b = mod(11);
  c = mod(12);
  d = mod(13);
  % apply constraints
  if b < 0 || c <= 0 || d < 0 
      cost = inf;
      return;
  end
  
  % calculate a, so that the max value is the same    
  a = (Fss - d)/((Lmax -b)^c);
  Force_par{j} = a*max(Length{j} - b, 0).^c + d;
  % calc force
  Force{j} = Force{j} + Force_par{j}; 
    
  if drawAllStates
        nexttile([1 4]);
        semilogx(Time{j}, Force{j}, 'k-');hold on;
        scatter(Time{j}(i_time_snaps), Force{j}(i_time_snaps), 'ko', 'filled');
        xlim([1e-3 max(1300, rds(j)*10) + 30]);
        % nexttile([1 1]);
        % loglog(Time{j}, Force{j});hold on;
        % scatter(Time{j}(i_time_snaps), Force{j}(i_time_snaps), 'o', 'filled');
        % xlim([1e-3 inf])
        xlabel('$t$ (s)', Interpreter='latex');
        legend('Tension', 'Insets', 'Location','northwest');
        ylabel('$\Theta$ (kPa)', Interpreter='latex')
        fontsize(12, 'points')
        exportgraphics(f,sprintf('../Figures/States%g_%gs.png', pCa, rds(j)),'Resolution',150)
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
    figure(900 + j*10 + round(pCa));clf;    tiledlayout(2,2, TileSpacing="compact");
    set(gcf, 'Position', [500  300  7.2*96 7.2*96/aspect])
    rng = Time{j} < 20;
    t = Time{j}(rng)';
    x = Length{j}(rng) + 0.95;
    y = Force{j}(rng);
    y2 = Force{j} - Force_par{j};
    z = zeros(size(t));
    col = linspace(0, t(end), length(t));   
    x_ = @(t) t - floor(t/Tc)*Tc;
    lw = 1

    nexttile;
    plot(t, x, 'k', linewidth = lw);
    % plot(t, x, t, y);legend('Length', 'Force');
    xlabel('$t$ (s)', Interpreter='latex');    ylabel('$L$ ($\mu$m)', Interpreter='latex');
    ylim([0.95, 1.2])

    nexttile(2, [2 1]);
    surface([x;x],[y;y],[z;z],[col;col],...
            'facecol','no',...
            'edgecol','interp',...
            'linew',lw);
    clim([0 50 100])
    colormap(turbo);    
    ylabel('$\Theta$ (kPa)', Interpreter='latex');     xlabel('$L$ ($\mu$m)', Interpreter='latex');
    cb = colorbar;  cb. title(cb, 't (s)');

    nexttile;
    % plot(x_(t), x, x_(t), y);legend('Length', 'Force');
    plot(t, y, 'k', linewidth = lw);
    xlabel('$t$ (s)', Interpreter='latex');
    ylabel('$\Theta$ (kPa)', Interpreter='latex');
    fontsize(12, 'points');
    exportgraphics(gcf,sprintf('../Figures/FigHysteresis%g_%gs.png', pCa, rds(j)),'Resolution',150)
  end
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
    Ep = 0;
end

cost = Ep*100 + sum([En{1:end}], 'all');

if exist('drawPlots', 'var') && ~drawPlots
    return;
end

Tarr = t_int;Farr = Ftot_int(1:4);
if exportRun
    % save data for decay overlay loglog plot
    Tarr = t_int;Farr = Ftot_int(1:4);
    % if pCa < 10
    %     save pcagraphingenv.mat
    % else
    %     save relaxgraphingenv.mat
    % end
    % save(sprintf('..\\pca%gmodeldata.mat', pCa), 'Tarr', 'Farr')
    % save('..\pca11modeldataDoubleStates.mat', 'Tarr', 'Farr')
end

try
% tic
set(groot,'CurrentFigure',figInd); % replace figure(indFig) without stealing the focus
cf = clf;

aspect = 2;
% normal size of 2-col figure on page is 7.2 inches
% matlab's pixel is 1/96 of an inch
% f.Position = [300 200 7.2*96 7.2*96/aspect];

% f.Position = [300 200 7.2*96 7.2*96/aspect];
% colors = lines(max(rampSet)+1); colors(1:end-1, :) = colors(2:end, :);
colors = gray(5);
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

ym = 0;
for j = max(rampSet):-1:1
    if isempty(Force{j})
        continue;
    end
    hd{j} = errorbar(datatables{j}.Time-2,datatables{j}.F,datatables{j}.SD, '-', LineWidth=2, Color=colors(3, :), CapSize=0);
    set([hd{j}.Bar, hd{j}.Line], 'ColorType', 'truecoloralpha', 'ColorData', [hd{j}.Line.ColorData(1:3); 255*0.4])
    hm{j} = semilogx(Time{j},Force{j},'-', 'linewidth',2, 'Color', 'k'); 
    hph{j} = plot(nan, nan, 'x',Color=[1 1 1]); % just a placeholder
    % set(h.Cap, 'EdgeColorType', 'truecoloralpha', 'EdgeColorData', [h.Cap.EdgeColorData(1:3); 255*alpha])
    ym = max([ym Force{j}, datatables{j}.F']);
    set(tile_semilogx, 'TickLength', [0.0125 0.05]);
    set(tile_semilogx, 'TickLabelInterpreter', 'latex');
    ylabel('$\Theta$ (kPa)', Interpreter='latex')
    
    if rds(j) > 1
        tp = [rds(j)*1.4, max(Force{j})*1.2];
    else
        tp = [rds(j)*0.95, max(Force{j})*0.95];
    end
    text(tp(1), tp(2), ...
        sprintf('$t_r$ = %g s', rds(j)),...
        'horizontalAlignment', 'right', VerticalAlignment='bottom', Interpreter='latex');    
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
hl = legend([hd{1} hm{1}], ...
    'Data',...    
    'Model',...
    NumColumns=3);

hl.ItemTokenSize = [30, 20];
hl.Box = 'off';

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
return;
%% best fit from FigFitDecayOverlay
if pCa > 10
    x = [4.7976    0.2392    4.8212];
    [c rspca] = evalPowerFit(x, Force, Time, 'loglogOnly', [], false);
else
    x = [17.9381    0.2339    4.3359];
    [c rspca] = evalPowerFit(x, Force, Time, 'loglogOnly', [], true);
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
    xl = xlim;
    ylim(yl)
    yticks([0 20 40]);
    xticks([0 rds(j), rds(j)*2])
    if j == 4
        ylabel('$\Theta$ (kPa)', Interpreter='latex')
    else
        yticklabels([]);
    end
    xlabel('$t$ (s)', Interpreter='latex', HorizontalAlignment='left');
    set(sp, 'TickLength', [0.05 0.05]);
    set(sp, 'TickLabelInterpreter', 'latex')
    tit = text(xl(2), min(yl(2)-5, max(Force{j})), sprintf('$t_r$ = %g', rds(j)), 'Interpreter', 'latex', ...
        HorizontalAlignment='right', VerticalAlignment='bottom');
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
return
%%
% zoomIns = [0 200 0 10;...
%            0 20 0 15;...
%            0 2 0 15;...
%            0 .2 0 15;...
%            0 .04 0 15;...
%            ];

clf;
colors = lines(max(rampSet)+1);

% max out of all
ym = ceil( max(cell2mat(Force)) / 5 ) * 5;

% prepare in advance so that it wont draw over my inset
tile_semilogx = subplot(212);hold on;
pos = get(tile_semilogx, 'Position');
ylabel('Tension (kPa)')
xlabel('Time (s)')
set(gca,'Fontsize',14)
title(sprintf('Force response to %.2g ML ramp-up at pCa=%g, costing %1.4e€', Lmax, pCa, cost), 'Parent',tile_semilogx);
set(tile_semilogx, 'XLim', [-1 300]);
set(tile_semilogx, 'YLim', [0 ym*1])

% shift of peaks to have the same tail - just guessed
% shift = [-94, -7.2, -0.35, 0];
% based on pCa 11 shift in data
shift = [5.4, 0.82, 0.22, 0.01] - [100 10 1 0.1];
ymaxScale = 0;
for j = max(rampSet):-1:1
    if isempty(Force{j})
        continue;
    end
% figure(j); clf; axes('position',[0.15 0.15 0.8 0.80]); hold on; box on;
    % subplot(1, 3, j);hold on;
%% primary plot - semilog
    subplot(221)
    semilogx(datatables{j}.Time-2,datatables{j}.F,'-','linewidth',2, 'Color', [colors(j+1, :), 0.3]);
    hold on;
    
    yyaxis right;
    semilogx(t_int{j},cumsum(Es{j}),'--','linewidth',2, 'Color', [colors(j+1, :), 0.3]);
    yyaxis left;
    
    semilogx(Time{j},Force{j},'-', 'linewidth',1, 'Color', colors(j+1, :)*0.8); 
    % semilogx(t_int{j},Es{j},':', 'linewidth',1, 'Color', colors(j+1, :)*0.9); 
    axis([1e-2, 3e2, 0, ym]);  
    set(gca,'Fontsize',14)
    title('Tension response to muscle length ramp-up')
    xlabel('Time (s)')
    ylabel('Tension (kPa)')

%% other view - shifted to see the tail overlap

    subplot(222)
    % Estimating the true offset: Fss = C*(Tss)^-alpha + Fss_true;
    % Fss = Force_par{j}(end); % "steady state" at the end 
    % tss = Time{j}(end) - rds(j);
    % Fss_true = Fss - (4.22*tss^-0.21);
    Fss_true = Force_par{j}(end);
    Fss_true = Force{j}(end) - Force_par{j}(end);
% shift(j) = 0;
    loglog(datatables{j}.Time-2 + shift(j),datatables{j}.F - Fss_true,'-','linewidth',2, 'Color', [colors(j+1, :), 0.3]);
    hold on;
    loglog(Time{j} + shift(j),Force{j} - Fss_true,'-', 'linewidth',1, 'Color', colors(j+1, :)*0.8); 
    % plot(Time{j} + shift(j),Force{j},styles{j}, 'linewidth',1, 'Color', colors(j+1, :)*0.8); 
    
    % semilogx(t_int{j},Es{j},':', 'linewidth',1, 'Color', colors(j+1, :)*0.9); 
    % axis([1e-2, 1e2, 0, ym]);  
    xlim([1e-2, 3e2]);
    ymaxScale = max(ymaxScale, max(Force{j} - Fss_true));
    yminScale = (Force{j}(end) - Fss_true); % this should be around the same
    ylim([max(1e-2, 0.8*yminScale), 1.2*ymaxScale])
    if j == 1
        % only after the last one
        legend('Ramp 10s (Data)', 'Ramp 10s (Model)', 'Ramp 1s (Data)', 'Ramp 1s (Model)', 'Ramp 0.1s (Data)', 'Ramp 0.1s (Model)');
    end
    % legend('Ramp 10s, shifted by -8.6s', 'Ramp 1s, shifted by -0.78', 'Ramp 0.1s' );

    set(gca,'Fontsize',14)
    title('Tension response to muscle length ramp-up: shifted peaks')   

%% secondary plots - timebase. Need to cut out

    plot(datatables{j}.Time-2,datatables{j}.F,'-','linewidth',2, 'Color', [colors(j+1, :), 0.15], Parent=tile_semilogx);
    plot(Time{j},Force{j},'-', 'linewidth',1, 'Color', colors(j+1, :)*0.8, Parent=tile_semilogx);
    plot(t_int{j},Es{j},'--|','linewidth',2, 'Color', [colors(j+1, :), 0.3], Parent=tile_semilogx);

    % plot(t_int{j},Ftot_int{j},'r','linewidth',2.5);
    % max out of all

    
    ylabel('Tension (kPa)')
    xlabel('Time (s)')
    % set(gca,'Xtick',0:50:200)
    
    % zoom-in inset
    
    w = pos(3)*0.18; h = pos(4)*0.5;
    x = pos(1) + pos(3) - (j*1.3 - 0.3)*w; y = pos(2) + pos(4) - h;
    axes('Position',[x, y, w, h]);hold on;
    % axes('position',[0.5 0.5 0.4 0.4]); hold on; box on;
    plot(datatables{j}.Time-2,datatables{j}.F,'-','linewidth',3, 'Color', [colors(j+1, :), 0.15]);
    plot(t_int{j},Es{j},'--','linewidth',2, 'Color', [colors(j+1, :), 0.3]);
    % plot(Time{j},Force{j}, 'r:', 'linewidth',2); 
    plot(t_int{j},Ftot_int{j},'-','linewidth',2, 'Color', colors(j+1, :)*0.9);
    ym_inset = ceil( max(Force{j}) / 10 ) * 10;
    axis([0, rds(j)*2, 0, ym_inset]);  
    set(gca,'Fontsize',14)    
    
end
% cla;
pos1 = get(subplot(221), 'Position');
w = pos1(3)*0.45; h = pos1(4)*0.4;
x = pos1(1) + pos1(3) - w; y = pos1(2) + pos1(4) - h;
axes('Position',[x, y, w, h]);
semilogx(PeakData(:, 1), PeakData(:, 2), 'ko', LineWidth=2);hold on;
semilogx(PeakData(:, 1), PeakModel, 'x', 'MarkerEdgeColor', [1 1 1]*0.5, LineWidth=2, MarkerSize=8);
axis([1e-1 1e2 0 max([ym PeakData(:, 2)'])])
semilogx(PeakData(:, 1), PeakModel, '--', Color=[1 1 1]*0.5, LineWidth=1);
legend('Peaks (Data)', 'Peaks (Model)', 'Location', 'southeast')
%%
h = annotation('textbox', [0.07 0.95 0 0], 'String', 'A)', 'FitBoxToText', false, 'FontSize', 32, 'FontWeight','bold');
h = annotation('textbox', [0.5 0.95 0 0], 'String', 'B)', 'FitBoxToText', false, 'FontSize', 32, 'FontWeight','bold');
h = annotation('textbox', [0.07 0.5 0 0], 'String', 'C)', 'FitBoxToText', false, 'FontSize', 32, 'FontWeight','bold');
%%
catch e
    disp(e.message)
end
%%
return;
Tss =  (0.50e6)*0.225^8;
figure(11); 
loglog(Time{4}-0.1,Force{4}-Tss, ...
       Time{3}-1.0,Force{3}-Tss, ...
       Time{2}-10,Force{2}-Tss, Time{1}-10,Force{1}-Tss); 
grid
save('../modeltesting2.mat', 'Time', 'Force', 'Tss')
% fig = gcf;
% set(gcf, 'Position', [50 50 1200 700])
% saveas(fig, ['..\Figures\Fig_' fig.Name], 'png')
return;

%% Overlap plots

%% Draw plots
figure(1001);clf;hold on;legend()
% plot n.1: 
colororder(jet(Ng));
for n = 1:1:Ng
    
    plot(outStruct{1, 1}.pu(:, n), 'x-', LineWidth=2)
    plot(outStruct{1, 2}.pu(:, n), 'x--', LineWidth=2)
    plot(outStruct{1, 3}.pu(:, n), 'x:', LineWidth=2)
end

% plot u to s for different N

% Visualize states in time?
% / Passive resting fit - both semilog and linear?
% / Shift the peaks so the tails overlap?
% / Fit for maximal CA
% Values of the Ca sensitive params - kA, KD?, kd
