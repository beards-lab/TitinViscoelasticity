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

%% Test while recreating the relaxed rptocol
clear;

figure(2); clf;
datasetname = 'StitchingRelaxed_1_1_red.csv';
pCa = 11;
simtype = 'velocitytable_relaxed.csv';
% simtype = 'velocitytable_relaxed_limhold.csv';
alphaF_0 = 1;
% drawFig1 = true;
drawAllStates = 1;
time_snaps = [270, 360, 898, 1307];
RunCombinedModel;
cost
% hold on;
% plot(Time{1}, Force_par{1}, '--', LineWidth=1)
% set(gca(), "XScale", 'linear');xlim([-inf inf]);ylim([-inf inf]);

% figure(12); plot(states{1})
%% nexttile;plot(Time{1}, states{1})
figure(10);clf;hold on;
num_states = size(states{1}, 2);
cmap = lines(num_states); % Generate 15 distinct colors
for i = 1:num_states
    plot(Time{1}, states{1}(:, i), 'Color', cmap(i,:), 'LineWidth', 1.5);
end
legend

%% resim sin
alphaF_0 = 1;
simtype = 'sin';
pCa = 11;
rampSet = [4 3];
RunCombinedModel;
%% Construct the whole protocol
% clf;
Lmax = 0.225;
% mouse HR = 300
% f_hr = 300/60;
% T = 0.2;
Rates = [0.1, 0.2 0.4];
% pos_slacks = @(t) (t < 0.1).*t*(-(Lmax+0.15)./0.1) + (t >= 0.1 & t < 20)*(-0.15) + 0.15*10*(t>=20 & t < 20 +10).*(t-20) + (t >= 30)*0;
% Define the piecewise function
% syms t Tc;
% Tc = 0.1;
clear t Tc
numbeats = 30;

% t_off = numbeats*Tc;

% sin
% segment = piecewise( ...
%     t < 120 - numbeats*Tc, 0, ...
%     t < 120, (-Lmax/2 * sin(2*pi*(t-120- numbeats*Tc)/Tc + pi/2) + Lmax/2), ......
%     t < 120 + 0.1, - 0.15/0.1 * (t-120), ... % Linear drop to -0.15 in 0.1s
%     t < 120 + 0.1 + 19.9, -0.15, ...               % Stay at -0.15 for 20s
%     t < 120 + 0.1 + 19.9 + 10, -0.15 + (0.15 / 10) * (t - (120 + 0.1 + 19.9)), ... % Rise to 0 in 10s   
%     t >= 120 + 0.1 + 19.9 + 10, 0 ...                   % Stay at 0 after 30.1s
% );

segment = @(t, Tc) ...  
    (t < 120) .* 0 + ...  
    ((t >= 120) & (t < 120 + numbeats*Tc)) .* (-Lmax/2 * sin(2*pi*(t-120)/Tc + pi/2) + Lmax/2) + ...  
    ((t >= 120 + numbeats*Tc) & (t < 120 + numbeats*Tc + 0.1)) .* (- 0.15/0.1 * (t-(120 + numbeats*Tc))) + ...  
    ((t >= 120 + numbeats*Tc + 0.1) & (t < 120 + numbeats*Tc + 0.1 + 19.9)) .* (-0.15) + ...  
    ((t >= 120 + numbeats*Tc + 0.1 + 19.9) & (t < 120 + numbeats*Tc + 0.1 + 19.9 + 10)) .* (-0.15 + (0.15 / 10) * (t - (120 + numbeats*Tc + 0.1 + 19.9))) + ...  
    (t >= 120 + numbeats*Tc + 0.1 + 19.9 + 10) .* 0;  

segment_dt = @(t, Tc) ...
    (t < 120) .* 0 + ...  % t < 120
    ((t >= 120) & (t < 120 + numbeats*Tc)) .* (-Lmax * pi / Tc * cos(2*pi*(t-120)/Tc + pi/2)) + ...  % 120 <= t < 120 + numbeats*Tc
    ((t >= 120 + numbeats*Tc) & (t < 120 + numbeats*Tc + 0.1)) .* (-0.15 / 0.1) + ...  % 120 + numbeats*Tc <= t < 120 + numbeats*Tc + 0.1
    ((t >= 120 + numbeats*Tc + 0.1) & (t < 120 + numbeats*Tc + 0.1 + 19.9)) .* 0 + ...  % 120 + numbeats*Tc + 0.1 <= t < 120 + numbeats*Tc + 0.1 + 19.9
    ((t >= 120 + numbeats*Tc + 0.1 + 19.9) & (t < 120 + numbeats*Tc + 0.1 + 19.9 + 10)) .* (0.15 / 10) + ...  % 120 + numbeats*Tc + 0.1 + 19.9 <= t < 120 + numbeats*Tc + 0.1 + 19.9 + 10
    (t >= 120 + numbeats*Tc + 0.1 + 19.9 + 10) .* 0;  % t >= 120 + numbeats*Tc + 0.1 + 19.9 + 10



% % sawtooth
% Define the function  
% segment = @(t, Tc) ...  
%     (t < 120) .* 0 + ...  
%     ((t >= 120) & (t < 120 + numbeats*Tc)) .* (Lmax/2 * sawtooth(2*pi*(t-120)/Tc, 0.5) + Lmax/2) + ...  
%     ((t >= 120 + numbeats*Tc) & (t < 120 + numbeats*Tc + 0.1)) .* (- 0.15/0.1 * (t-120 - numbeats*Tc)) + ...  
%     ((t >= 120+ numbeats*Tc + 0.1) & (t < 120+ numbeats*Tc + 0.1 + 19.9)) .* (-0.15) + ...  
%     ((t >= 120+ numbeats*Tc + 0.1 + 19.9) & (t < 120 + numbeats*Tc + 0.1 + 19.9 + 10)) .* (-0.15 + (0.15 / 10) * (t - (120 + 0.1 + 19.9 + numbeats*Tc))) + ...  
%     (t >= 120+ numbeats*Tc + 0.1 + 19.9 + 10) .* 0;  
% 
% segment_dt = @(t, Tc) ...
%     (t < 120) .* 0 + ...  % Derivative is 0 for t < 120
%     ((t >= 120) & (t < 120 + numbeats*Tc)) .* (Lmax/2 * (2*pi/Tc) * (sawtooth(2*pi*(t-120)/Tc, 0.5) - circshift(sawtooth(2*pi*(t-120)/Tc, 0.5), 1))) + ...  % Derivative of the sawtooth
%     ((t >= 120 + numbeats*Tc) & (t < 120 + numbeats*Tc + 0.1)) .* (-0.15/0.1) + ...  % Derivative for the linear piece
%     ((t >= 120+ numbeats*Tc + 0.1) & (t < 120+ numbeats*Tc + 0.1 + 19.9)) .* 0 + ...  % Constant segment
%     ((t >= 120+ numbeats*Tc + 0.1 + 19.9) & (t < 120 + numbeats*Tc + 0.1 + 19.9 + 10)) .* (0.15 / 10) + ...  % Linear segment
%     (t >= 120+ numbeats*Tc + 0.1 + 19.9 + 10) .* 0;  % After this point


% f = matlabFunction(pos_slacks, 'File', 'pieceofshit.m')
t_vals = linspace(-100, 170, 153000); % Generate time points
% y_vals = double(subs(segment, {t, Tc}, {t_vals, 0.1})); % Evaluate function at each t
% clf;hold on;
% plot(t_vals, segment(t_vals, 0.1), '-|', t_vals, segment_dt(t_vals, 0.1));
hold on; plot(t_vals, segment(t_vals, 0.2), '-|', t_vals, segment_dt(t_vals, 0.2));
% plot(t_vals, segment(t_vals, 0.4), '-|', t_vals, segment_dt(t_vals, 1));
%%
% clf;
% f_d = @(Tc_) 150;% double(30.1 + numbeats*Tc_ + 120);
% duratios = [-f_d(1) + numbeats*1,f_d(1)- numbeats*1, f_d(0.1), f_d(0.2), f_d(0.4)];
duratios = [-100, 100, 30 + 0.01, 120, numbeats*0.1 +   30, 120, numbeats*0.2 +   30, 120, numbeats*1 +   30];
cumtimes = cumsum(duratios);
times = [cumtimes(1:end-1);cumtimes(2:end)]';
velocities = { 0,...
                @(t_)segment_dt(t_ + 120+numbeats, 1),... % slack only
                @(t_)segment_dt(t_- cumtimes(3), 0.1),... % cut in two so the integrator does not skip at max step
                @(t_)segment_dt(t_- cumtimes(3), 0.1),...
                @(t_)segment_dt(t_- cumtimes(5), 0.2),...
                @(t_)segment_dt(t_- cumtimes(5), 0.2),...
                @(t_)segment_dt(t_- cumtimes(7), 1),...                
                @(t_)segment_dt(t_- cumtimes(7), 1),...                
                };
L0 = [0 0 0 0 0 0 0 0 0];

simtype = 'refolding';
pCa = 4.51;
% pCa = 11;
rampSet = [4];
% set new params based on pCa
clear params;
tic
RunCombinedModel;
toc

figure(161);clf;
nexttile(1); plot(Time{4}, Length{4}, '-|');
nexttile(2); plot(Time{4}, Force{4});
%%
figure(3);clf;
plot(Length{4})
%%

plot(Length{4}(3000:13000), Force{4}(3000:13000));
%%
f1 = @(t_) (t_ < 0 & t_ > -f_d(1) + numbeats*1).*pieceofhit(1, t_ + f_d(1)) + ... % use the slack tail only 
          (t_ < f_d(0.1) & t_ > 0) .* pieceofhit(0.1, t_) + ... 
          (t_ < f_d(0.2)*2 & t_ > f_d(0.1)) .* pieceofhit(0.2, t_ -f_d(0.1)) + ...
          (t_ < f_d(0.4)*3 & t_ > f_d(0.2)*2) .*pieceofhit(0.4, t_ -f_d(0.2)*2);


% test
% tic
v_vals = f1(t_vals);
% toc
% v_vals = [];
% tic 
% for i = 1:length(t_vals)
%     v_vals(i) = f1(t_vals(i));
% end
% toc
% 
% % Numerical integration to compute position (using cumulative trapezoidal rule)
p_vals = cumtrapz(t_vals, v_vals);  % Integrating velocity to get position
hold on; plot(t_vals, p_vals + double(subs(segment, {'t', 'Tc'}, [0 0.1])), '-');
% plot(t_vals, v_vals)
%%
simtype = 'refolding';
times = [-200 350];
velocities = {f1};
pCa = 11;
RunCombinedModel;

figure(161);clf;
nexttile;plot(Time{4}, Length{4});
nexttile; plot(Time{4}, Force{4});

%%
% positions = @(t) (t < 60) .* (-Lmax/2 * cos(2*pi*0.1*t/Tc) + Lmax/2) + ...
%                  ((t >= 60) & (t < 120)) .* 0 + ...
%                  ((t >= 120) & (t < 180)) .* (-Lmax/2 * cos(2*pi*0.2*(t-120)/Tc) + Lmax/2) + ...
%                  ((t >= 180) & (t < 240)) .* 0 + ...
%                  ((t >= 240) & (t < 300)) .* (-Lmax/2 * cos(2*pi*0.4*(t-240)/Tc) + Lmax/2);
% 
% 
%       times = [-100, 0;0 10*Tc];
% 
%       positions = @(t)-Lmax/2*cos(2*pi*t/Tc) + Lmax/2;
%       % differentiating positions
%       V = @(t)2/2*Lmax*pi/Tc *sin(2*pi*t/Tc);
% 
%       % syms fpos(t);   % fpos(t) = @(t)Lmax*sin(2*pi*t/Tc);
% 
%       velocities = {0, V};
%       L0 = 0;%Lmax/2;  

      t = 0:0.01:100;
      plot(t, pos_slacks(t))
