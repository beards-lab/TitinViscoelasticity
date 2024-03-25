x0 = ones(1, 10);

% % param modifiers
% 1  delU
% 2  kC
% 3  kS
% 4  alphaU
% 5  alphaF
% 6  nC
% 7  nS
% 8  nU
% 9 mu
% 10 Fss
% // 10 11 12 - params for ramp-up
% 13 Ls0
% 14 kA
% 15 kD
modNames = {'delU', 'kC', 'kS', 'alphaU', 'alphaF', 'nC', 'nS', 'nU', 'mu', 'b', 'c', 'd/Fps', 'Ls0', 'kA', 'kD', 'N/A'};
x0 = [ 1.0057    1.5981    1.3605    0.0917    1.4696    1.1553    1.0838    0.6758    1.1417]; % first shot

x0 = [ 0.9973    2.9945    2.8851    0.1308    -0.8659    1.3337    1.2446    0.7553     0.9829]; % negative param

% now, fixing the offset to 2.0 at 0.4 offset
x0 = [    1.0394    1.1020    0.8118    0.9860   0.7941    0.9136    0.9657    1.0193    1.4571]; % not really wokring

% optimizing with variable offset, not evaluating ramp-up, cost = 300.7
x0 = [0.9522    1.1349    1.0334 0.9123    0.4409    1.0839    0.9946    0.8851    1.1345    1.2716];

% optimized for all params, excl. Fss
x0 = [0.9522    1.1349    1.0334 0.9123    0.4409    1.0839    0.9946    0.8851    1.1345   0.6833    1.4778    0.8111 1 1];

% optimizing everything at once as a last step
x0 = [0.9736, 1.1519, 1.0295, 0.9306, 0.4439, 1.0758, 1.0011, 0.8704, 1.0829, 0.6774, 1.4252, 0.8157, 1.0302, 1.0118];

x0 = [1.0504    1.2978    0.9544    1.0226    0.4784    1.0429    0.9584    0.8259    0.8629    0.7200    1.3634    0.8411    0.9660    1.0150 1 1];

% optim for -log10 weighting
mod = [1.16970000000000	0.928400000000000	0.977400000000000	1.02340000000000	1.01370000000000	1.10320000000000	0.937900000000000	1.19500000000000	0.909900000000000	0.898800000000000	1	1	1];
%% set data set
% pCa = Inf; % First round of passive ramp-ups data, 200s long, dML 0.4
% pCa = 11; % Ca&PNB experiments, resting ramp-ups, 30s decay, dML = 0.225

mod = [0.9736, 1.1519, 1.0295, 0.9306, 0.4439, 1.0758, 1.0011, 0.8704, 1.0829, 0.6774, 1.4252, 0.8157, 1.0302, 1.0118 1 1 1 1];
mod = [1.0504    1.2978    0.9544    1.0226    0.4784    1.0429    0.9584    0.8259    0.8629    0.7200    1.3634    0.8411    0.9660    1.0150 1 1];


% test for pCa 11
mod(1) = 2.2;
mod(2) = 4.4;
mod(3) = 2.2;
% mod(5) = 2;
mod(6) = 1.03;
% mod(8) = 1.2;
mod(13) = -2.8;
%%
% mod = [2.0398    1.3113    3.8942    1.3500    0.4784 0.7398    0.8176    0.7869    0.8629    0.7200 1.3634    1   -2.8000    1.0150    0.6382 -0.5199];
% mod = [1.1592    1.0379    0.9763    0.9779    1.1237    1.0935    0.9365    1.0882    0.9846    0.8931];
% mod = [1.1697    1.0418    0.9774    0.9737    0.9858    1.0265    0.9403    1.0837    0.9889    0.8988, 1, 1, 1];
% mod = [1.16970000000000	0.928400000000000	0.977400000000000	1.02340000000000	1.01370000000000	1.10320000000000	0.937900000000000	1.19500000000000	0.909900000000000	0.898800000000000	1	1	1 1 1];
% mod = [0.0185    0.8479    0.4307    1.0234    1.0326 0.5971    0.9379    1.1950    0.9099    0.8988 1.0000    1.4450    0.7510    1.2811    2.7365];
tic
saSet = 1:22;
[c0_11, cost11_sap, cost11_sam] = runSa(11, mod, saSet);
[c0_4, cost4_sap, cost4_sam] = runSa(4.4, mod, saSet);
% normalize
cost11_sam = cost11_sam/c0_11;
cost11_sap = cost11_sap/c0_11;
cost4_sam = cost4_sam/c0_4;
cost4_sap = cost4_sap/c0_4;
toc
%% plot the result
% modNames = {'k_p(NoCa)', 'k_d', 'n_p', 'n_U', 'n_d', 'alphaU', 'k_{PEVK,A}', 'k_{PEVK,D}', 'k_p(highCa)', 'Fss', 'b', 'c', 'd', 'mu', 'alphaF_0'};

c0_11 = 1;
c0_4 = 1;

% identify tradeoffs
eps = (c0_4 + c0_11)*1e-3;
% better result
bett = (cost11_sap + cost4_sap + eps < c0_11 + c0_4) ...
    | (cost11_sam + cost4_sam + eps < c0_11 + c0_4);
betts = strings(1, length(cost4_sap));betts(bett) = "$$$";
% tradeoffs between low and high Ca - suggesting some dependency?
tdoCa = (cost11_sap + eps < c0_11 & cost4_sap > c0_4 + eps) | (cost11_sam + eps < c0_11 & cost4_sam > c0_4 + eps) ...
    | (cost11_sap > c0_11 + eps & cost4_sap + eps < c0_4) | (cost11_sam > c0_11 + eps & cost4_sam + eps < c0_4);
tdoCas = strings(1, length(cost4_sap));tdoCas(tdoCa) = "*";

figure(101);clf; 
b = bar([cost11_sap; cost4_sap; zeros(size(cost4_sap)); cost11_sam;cost4_sam]'); hold on;
title('Sensitivity Anal')
plot([0.5, length(cost4_sap)], [c0_11 c0_11], 'b--')
plot([0.5, length(cost4_sap)], [c0_4 c0_4], 'm--')
title('Grouped sensitivity analysis (* indicates Ca trade-off, $$$ potential of improvement)')
legend('10xpCa11+','pCa4+','', '10xpCa11-','pCa4-','pCa11_0','pCa4_0')
set(gca,'XTick',saSet);
set(gca,'XTickLabels',strcat(string(1:length(cost4_sap)), ':', modNames(1:length(cost4_sap)),tdoCas, betts));
set(gca, "FontSize", 14)
%%
% close all
tic
modNames = {'k_p(NoCa)', 'k_d', 'n_p', 'n_U', 'n_d', 'alphaU', 'k_{PEVK,A}', 'k_{PEVK,D}', 'k_p(highCa)', 'Fss', 'b', 'c', 'd', 'mu', 'alphaF_0'};

% optim for -log10 weighting
% mod = [1.16970000000000	0.928400000000000	0.977400000000000	1.02340000000000	1.01370000000000	1.10320000000000	0.937900000000000	1.19500000000000	0.909900000000000	0.898800000000000	1	1	1 1 1];
% reoptim on avg data
% mod = [0.0165    0.7035    0.4403    1.0234    1.0077    0.5754    0.9379    1.1950    0.9099    0.8988    1.0000    1.1421    1.4792    1.1156    2.9834];
% reoptim on avg data incl. 20231102
% mod = [0.0185    0.8479    0.4307    1.0234    1.0326 0.5971    0.9379    1.1950    0.9099    0.8988 1.0000    1.4450    0.7510    1.2811    2.7365];
% modSel = [2 3 4 6 7 8 9 12 14 15];
% mod(16) = 0;
% evalCombined(mod(modSel))
% mod = [1.1592    1.0379    0.9763    0.9779    1.1237    1.0935    0.9365    1.0882    0.9846    0.8931];
% modSel = [7, 9];
% mod_pca6 = [0.1657, 0.3895];
% mod = [1.1697    1.0418    0.9774    0.9737    0.9858    1.0265    0.9403    1.0837    0.9889    0.8988 1 1 1];
% for absolute average - relaxed
% mod = [0.0343    0.7101    0.4731    1.0234    1.0916    1.9353    0.9379 1.1950    0.9099    0.8988    0.5952    2.0416    0.7510    1.2811 4.1891];
% for absolute average - pca 4.4
% mod = [0.0343    0.7101    0.4731    1.0234    1.0916    1.9353    1.5266 0.8291    0.0356    0.8988    0.5952    2.0416    0.7510    1.2811 4.1891];
% for absolute average excluding the ramp-up
% mod = [0.0142    0.4204    0.4267    1.0234    0.7567 0.3025    0.1535    1.0318    0.0133    0.8988 0.8679    4.5429    0.7510    1.2811    1.6208];
%% optim for filtering out remaining force (attempt 01)
mod = [0.0231    0.2275    0.4870    1.0234    0.7909   0.2929    0.1113    0.2652    0.0218    0.8988    0.7426    1.8401    0.7510    1.2811    1.6663];
mod = [0.0231    0.2275    0.4870    1.0234    0.7909   0.2929    0.1113    0.2652    0.0218    0.8988    0.7426    1.8401    0.7510    5.2811    1.6663];
figure(99);polarplot(linspace(2*pi/15, 2*pi, 15), mod, 'x-');evalCombined(mod, mod, 1:15)
% optim for pCa 11 only, from original estimates, without ramp up
% modSel = [1 2 3 5 6 10 15];
% mod = 0.0278    0.7953    0.4768    1.0000    1.0770    1.7512    1.0000    1.0000    1.0000    1.3607    1.0000    1.0000    1.0000    1.0000   -0.1149
% fixing the refolding to non-negative
% mod =[0.0299    0.8084    0.4798    1.0000    1.0862    1.6973    1.0000    1.0000    1.0000    1.2983    1.0000    1.0000    1.0000    1.0000    0.50000];
% figure(99);polarplot(linspace(2*pi/15, 2*pi, 15), mod, 'x-');evalCombined(mod, mod, 1:15)
% retuned for reduced mu
% mod =  [0.0287    0.7984    0.4766    1.0000    1.0715    1.8054    1.0000    1.0000    1.0000    1.3254    1.0000    1.0000    1.0000    1.0000    0.5000];
% figure(99);polarplot(linspace(2*pi/15, 2*pi, 15), mod, 'x-');evalCombined(mod, mod, 1:15)
% Retuned for ramp-ups, strange bump in 10s
% mod =  [0.0293    0.8574    0.4871    1.0000    0.9528    1.7817    1.0000    1.0000    1.0000    1.3369    3.6372    0.2425    0.0030    0.1000    0.5000];
% figure(99);polarplot(linspace(2*pi/15, 2*pi, 15), mod, 'x-');evalCombined(mod, mod, 1:15)
%% retuned for pCa 11 and 4.4, from scratch except for ramp-up and fixed mu and alphaR0 
% Candidate 1
tic
modSel = [1 2 3 4 5 6 7 8 9 10]; mod = [0.4852    0.2070    1.0403    1.1617    0.7393    1.2668    1.3151    1.5592    1.1949    1.6778    3.6372    0.2425    0.0030    0.1000    0.5000];
figure(99);polarplot(linspace(2*pi/15, 2*pi, 15), mod, 'x-');evalCombined(mod, [mod 1 1], 1:15)
toc
%% yet another set
mod =  [0.0299    0.8084    0.4798    1.0000    1.0862    1.6973    1.0000    1.0000    1.0000    1.2983    1.0000    1.0000    1.0000    1.0000    0.50000];
figure(99);polarplot(linspace(2*pi/15, 2*pi, 15), mod, 'x-');evalCombined(mod, mod, 1:15)
%% retuned for reduced mu
modSel = [1 2 3 4 5 6 7 8 9 10]; mod =  [0.0287    0.7984    0.4766    1.0000    1.0715    1.8054    1.0000    1.0000    1.0000    1.3254    1.0000    1.0000    1.0000    1.0000    0.5000];
figure(99);polarplot(linspace(2*pi/15, 2*pi, 15), mod, 'x-');evalCombined(mod, mod, 1:15)
%% retuned for both pCas in log space
% candidate 2
mod = [0.4804    0.1794    0.9727    1.6729 0.7547   1.1479e+03 0.0245    0.1301    0.6788    1.4584 3.6372    0.2425    0.0030    0.1000 0.5000];
figure(99);polarplot(linspace(2*pi/15, 2*pi, 15), mod, 'x-');evalCombined(mod, mod, 1:15)
%% retuned for pCa 4 in log space only, incl. ramp-up
mod = [0.4804    0.0661    0.9735    2.0109    0.5917 1.6109e+05 0.1694    0.2900    0.7105    1.4348    0.3534    0.1366    0.0002 0.1000    0.5000];
figure(99);polarplot(linspace(2*pi/15, 2*pi, 15), mod, 'x-');evalCombined(mod, mod, 1:15)
%% retuned for pCa 4 in log space only, w/o. the ramp-up
mod = [0.4804    0.0667    0.9711    1.9971    0.5672 2.0654e+05 0.0919    0.3252    0.6417    2.0434    0.3652    0.0866    0.0002 0.1000    0.5000];
figure(99);polarplot(linspace(2*pi/15, 2*pi, 15), mod, 'x-');evalCombined(mod, [mod 1 1], 1:15);
%% Tuned for data readjusted for remaining force and scaled by Fmax - with 1.4e4
mod = [0.6108    0.0578    1.0426    1.1110    0.5432    1.3586    0.9776    1.3503 1.2387    1.8919    3.6372    0.2425    0.0030    0.7421    0.4327    0.9964 1.0091];
figure(99);plotParams(mod);evalCombined(mod, mod, 1:17)
%% FineTuning for b, c and d
mod = [0.6108    0.0578    1.0426    1.1110    0.5432    1.3586    0.9776    1.3503    1.2387    1.8919    0.001    1.0877    0.001    0.7421    0.4327    0.9964    1.0091];
figure(99);plotParams(mod);evalCombined(mod, mod, 1:17)
%% Retuned for Fmax normalized data without Frem adjustment, incl. mu
mod = [0.4847    0.0514    1.0389    1.1202    0.5064    1.3293    0.9757    1.3062    1.3521    1.9492    0.0003    1.0841    0.0007    0.7551    0.6101    0.9964    1.0090];
% mod = [0.4777    0.0503    1.0404    1.1199    0.5045    1.3213     0.9756    1.2949    1.3656    1.9521    0.0003    1.0843    0.0008    0.7595    0.6170    0.9964    1.0089];
figure(99);plotParams(mod);evalCombined(mod, mod, 1:17)
%% unfinished finetuning for Fmax normalized data without Frem with free PEVK_A_Relaxed, costs 3e4
mod = [0.2651    0.0171    1.0557    1.1643    0.3515    1.4502    0.9884    1.1832    1.0980    2.0724    0.0000    0.8893    0.0027    0.6528    0.6272    0.9964    1.0068];
figure(99);plotParams(mod);evalCombined(mod, mod, 1:17)
%% Fmax normalized data with power law decay correction, started with manual ramp onset with 1.6e4
mod = [0.6798    0.0944    1.0921    1.2789    0.7231   1.4513    0.1460    0.1229    1.2924    1.5949    1.0000    0.5000         0    0.8211    0.9438    0.6262    0.6327];
figure(99);plotParams(mod);evalCombined(mod, mod, 1:17)
%% additionally, allowing optim the delU and Lref
mod = [0.2723    0.0474    1.0668    1.3219    0.4870    1.4287    0.1309    0.3002    0.5014    1.8149    1.0000    0.5000     0    0.5101    0.9050    0.8296    1.2636    1.0114    0.9628];
figure(99);plotParams(mod);evalCombined(mod, mod, 1:17)
%% running a little longer
mod = [0.2211    0.0356    0.9894    1.2774    0.4615    1.7193    0.1453    0.2687    0.3892    1.8224    1.0000    0.5000         0    0.4032    0.9260    0.8455    1.3289    1.0145    0.9811];
figure(99);plotParams(mod);evalCombined(mod, mod, 1:17)
%% and longer
mod = [0.1675    0.0321    0.9585    1.2787     0.4273    1.4061    0.0986    0.3406    0.2290    1.5851    0.4737    1.1618    0.0031    0.3933    1.4191    0.4482    1.6593    0.9094    0.6945];
figure(99);plotParams(mod);evalCombined(mod, mod, 1:17)

%% and longer until 6.3e3
clf
mod = [0.1403    0.0493    0.9278    1.2795    0.4052    1.1240    0.0933    0.6228    0.1716    1.5952    0.0008    1.4658    0.0023    0.4388    2.7788    0.6790    2.9644    0.9875    0.8030 1.1240];
plotParams(mod, [-2, 1]);hold on;
%% add alphaU_LowCa, optim at 0.1 and 10s only with 4.3e3, but 7.25e3 for all three and 9.6593e+03 for all four
clf;
mod = [0.1382    0.0875    0.9117    1.3178    0.4802    1.1331    0.1095    0.7093    0.1545    1.5725    0.0000    1.4225    0.0027    0.5243    2.5843    0.5718    2.7233    1.0000    1.0000    1.1331];
plotParams(mod, [-2, 3]);hold on;
evalCombined(mod, mod, 1:length(mod), [4.4 10])
%% reduced set - no low Ca PEVK attachment
% optim at 0.1 and 10s only with 3.8e3, but 7.6e3 for all three and 9.46e3 for all 4
mod = [0.2880    0.0763    0.9277    1.7014    0.5508  632.7487    0.2108    0.5573    0.2802    1.5806    1.0005    1.3387    0.1604    0.5548    1.1785    NaN    NaN    1.0000    1.0000  415.4241];
evalCombined(mod, mod, 1:length(mod), [4.4 11])
%% starting at reduced, but extending with non Ca PEVK
mod = [0.2520    0.0793    0.9654    1.7583    0.5718  612.9789    0.2057    0.4547    0.2984    1.6962    0.9569    1.3798    0.1658    0.5642    1.1109   0.1956    0.6015 1.0000    1.0000  454.1430];
% mod(setdiff(1:20, modSel)) = NaN;
plotParams(mod, [-2, 3]);hold on;
legend('Extended set', 'Reduced set')
evalCombined(mod, mod, 1:length(mod), [4.4 10])

%% test evaluate Ca sensitivity importance - REDUCED
% reduced candidate
% mod = [0.2880    0.0763    0.9277    1.7014    0.5508  632.7487    0.2108    0.5573    0.2802    1.5806    1.0005    1.3387    0.1604    0.5548    1.1785    NaN    NaN    1.0000    1.0000  415.4241 1];
evalCombined(mod, mod, 1:length(mod), [4.4])
f = figure(144)
saveas(f, 'rc_base.png')
%% no Ca dep on kp
m = mod;
m(9) = mod(1); % the multiplier in RunCombinedModel
evalCombined(m, mod, 1:length(mod), [4.4])
f = figure(144)
saveas(f, 'rc_no_kp.png')
%% no Ca dep on kd
m = mod;
m(22) = mod(2); % the multiplier in RunCombinedModel
evalCombined(m, mod, 1:length(mod), [4.4])
f = figure(144)
saveas(f, 'rc_no_kd.png')
%% no PEVK attachment
m=mod;
m(7) = 1e-6;% PEVK attachment
m(8) = 1e6; % PEVK detachment
evalCombined(m, mod, 1:length(mod), [4.4])
f = figure(144)
saveas(f, 'rc_no_PEVK.png')
%% no alphaU dep
m = mod;
m(20) = mod(6);% PEVK attachment
evalCombined(m, mod, 1:length(mod), [4.4])
f = figure(144)
title('rc_no_alphaU', 'Interpreter','none')
saveas(f, 'rc_no_alphaU.png')
%% test evaluate Ca sensitivity importance - EXTENDED
% extended candidate
mod = [0.1382    0.0875    0.9117    1.3178    0.4802    1.1331    0.1095    0.7093    0.1545    1.5725    0.0000    1.4225    0.0027    0.5243    2.5843    0.5718    2.7233    1.0000    1.0000    1.1331];
evalCombined(mod, mod, 1:length(mod), [4.4])
f = figure(144)
saveas(f, 'ec_base.png')
%% no Ca dep on kp
m = mod;
m(9) = mod(1)/4.78; % the multiplier in RunCombinedModel
evalCombined(m, mod, 1:length(mod), [4.4])
f = figure(144)
saveas(f, 'ec_no_kp.png')
%% no PEVK attachment
m = params;
m(7) = 1e-9;% PEVK attachment
m(8) = 1e9; % PEVK detachment
evalCombined(m, params, 1:length(mod), [4.4])
f = figure(144)
% saveas(f, 'ec_no_PEVK.png')


%% no alphaU dep
m = mod;
m(20) = mod(6);% PEVK attachment
evalCombined(m, mod, 1:length(mod), [4.4])
f = figure(144)
title('ec_no_alphaU', 'Interpreter','none')
saveas(f, 'ec_no_alphaU.png')
%% no PEVK at 0 Ca
m = mod;
m(7) = 1e-6;% PEVK attachment
m(8) = 1e6; % PEVK detachment
evalCombined(m, mod, 1:length(mod), [10])
f = figure(200)
title('ec_no_PEVK_at_pCa10', 'Interpreter','none')
saveas(f, 'ec_no_PEVK_at_pCa10.png')
%% Can I reproduce the fit with NO PEVK?
m = mod;
m = [m ones(1, 22-length(m))];
% evalCombined(m, [mod 1 1], 1:22, [11])

% disable PEVK
% m(7) = 1e-6; m(8) = 1;
m(16) = 1e-6;m(17) = 1e6;
m(16) = NaN; m(17) = NaN;
% increase the unfolding rate
% m(6) = 1.5;
m(7) = mod(7)*0.9;
m(8) = mod(8)*1.2;
% kp
% m(1) = mod(1)*0.5;
% kd
% m(2) = mod(2)*0.5;
m(20) = NaN;
m(21) = NaN;
% kd stiffening
m(22) = mod(2)*5*2.2;
% kp stiffening
m(9) = mod(1)*2.2;
evalCombined(m, [mod 1 1], 1:22, [11])
modSel = [1:10 14 22];
%%
return
tic
% modSel = [1 2 3 5 6 10];
% default is 403
% modSel = [7, 9];
mod = [mod ones(1, 20 - length(mod))]
modSel = 1:length(mod);
%%
tic
% mod(10) = 1; mod(11) = 1;mod(12) = 0.5;mod(13) = 0;mod(1) = 0.8;mod(2) = 0.12;
evalCombined(mod(modSel), mod, modSel, [11])
% evalCombined(mod_pca6)
% evalCombined(mod)
% evalCombined([1 1 1 1 1 1])
% rampSet = [2 3 4] % 18s
% rampSet = [3] % 6s
% rampSet = [4] % 3s
toc
%% Knockouts
m = params;
evalCombined(m, params, 1:length(mod), [11 4.4])
f = figure(144)
exportgraphics(f,sprintf('../Figures/FigKnockoutBase_pCa%g.png', 4.4),'Resolution',150)
f = figure(210)
exportgraphics(f,sprintf('../Figures/FigKnockoutBase_pCa%g.png', 11),'Resolution',150)


%% no PEVK attachment
m = params;
m(7) = 0;% PEVK attachment
m(8) = 1e9; % PEVK detachment
evalCombined(m, params, 1:length(mod), [4.4])
f = figure(144)
exportgraphics(f,sprintf('../Figures/FigKnockoutPEVK_pCa%g.png', 4.4),'Resolution',150)

%% little viscosity attachment
m = params;
m(14) = 1e-1;% 
evalCombined(m, params, 1:length(mod), [11])
%%
f = figure(144)
exportgraphics(f,sprintf('../Figures/FigKnockoutVisc_pCa%g.png', 4.4),'Resolution',150)
f = figure(210)
exportgraphics(f,sprintf('../Figures/FigKnockoutVisc_pCa%g.png', 11),'Resolution',150)


%% All four ramps costing 1e4
mod2params = [10293*0.7, 14122, 3.27, 6.0, 3.25, (8.4137e5)/0.7, 16.44, 14.977, 10203*4.78*0.7, 3.2470, 0.05, 7, 1, 1, 0.05, 0.1*16.44, 14.977, 0.9, 0.1*15 , (8.4137e5)*0.7, 3.2470, 14122, 0.1];
mod2params = [10203*0.7, 14122, 3.27, 6.0, 3.25, (8.4137e5)*0.7, 16.44, 14.977, 10203*4.78*0.7, 3.2470, 0.05, 7, 1, 1, 0.05, 0.1*16.44, 14.977, 0.9, 0.0125*14, 0.7*8.4137e5, 3.2470, 14122, 0.1];
% optimized Fss
mod = [0.1360    0.4811    0.8499    1.3379    0.6881    2.9849    0.0849    1.1904    0.1302    1.6995         0    1.4243    0.0027    0.7668    2.4840       NaN       NaN    1.0000    1.0000       NaN       NaN    0.4930];


% Fss fixed at data fitted 4.9
% mod = [0.1455    0.4841    0.8456    1.3729    0.6897    2.8866    0.0849    1.1904    0.1302    1.5060         0    1.4243    0.0027    0.7660    2.4840       NaN       NaN    1.0000    1.0000       NaN       NaN   0.4930];

% Fss fixed at 4.89 and baked into the params, 
% reoptimized at modSel = [1    2     3     4     5     6    14] for pCa 11
% mod = [0.2806    9.3110    0.9803    1.9501    0.9566 1.0098e+03 0.0849    1.1904    0.1302    1.5060         0    1.4243    0.0027    0.7076    2.4840       NaN       NaN    1.0000    1.0000       NaN       NaN    0.4930];

% reoptim for pca 11 AND 4.4 with 2.3e4 for all ramps
% mod = [0.2729    6.9247    0.9816    1.9550    0.9273  908.8536    0.0297    1.1195    0.1996    1.5060         0    1.4243    0.0027    0.8796    2.5534       NaN       NaN    1.0000    1.0000       NaN       NaN    4.9791];

% optim all ramps, pCa 11 AND 4.4 with 1.9e4 in 58s with modSel = [1     2     3     4     5     6     7     9    14    22];
% mod = [0.266    6.682    0.971    1.936    0.924    1033.800    0.021    1.119    0.212    1.506    0.000    1.424    0.003    0.946    2.553    NaN    NaN    1.000    1.000    NaN    NaN    5.310]

% optim for better tail only
% mod = [1.0777   16.5334    0.5891    0.6831    2.2466    1.4398    1.0000    1.0000    1.0000    1.0000    1.0000    1.0102    0.9400    1.0181    1.0000    1.0000    1.0000    1.0000    0.9985    1.0000    1.0000    1.0000];
mod(23) = 0;
mod(15) = 0;
params = mod.*mod2params;
params(10) = 4.89;
tic
evalCombined(params(modSel), params, modSel, [11])
toc
%% reinint to match Dan's model
% mod = ones(22, 1); modSel = [1:6 13 14];

% best fit to dan's init params
modSel = [1:6 12 14]; mod = [1.0906   16.8733    0.5964    0.6912    2.1513    1.4117    1.0000    1.0000    1.0000    1.0000    1.0000    1.0004    0.9400    1.0222    1.0000    1.0000    1.0000    1.0000    0.9985    1.0000    1.0000    1.0000];
mod = [0.1152   11.2004    0.4241    1.0295    1.8167    3.3966    1.0000 1.0000    1.0000    1.0000    1.0000   15.7907    1.0000   41.1979 1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000];
% resampling data in logspace to fit the peak better
mod = [0.1227    5.7413    0.3624    0.8583    1.6670    4.4089    1.0000   1.0000    1.0000    1.0000    1.0000   23.7214    1.0000   47.2949 1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000];
mod = [0.1035    6.5217    0.3626    0.8185    1.7903    4.7326    1.0000    1.0000    1.0000    1.3824    1.0000    1.0000    1.0000   84.2683    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000];
% mod = [0.0973    2.9994    0.3890    1.0250    1.5260    6.4763    1.0000    1.0000    1.0000    1.0056    0.9924    3.6707    0.9998  146.8509    1.0046    1.0075    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000];

% modSel = [1     2     3     4     5     6    14]; mod = [0.1499    6.4878    0.3322    0.7361    1.6326    4.6280    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    0.9201    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000];
% optim relaxed WITH PEVK, resampling the data log since START of the ramp - does not sample the
% peak as dense for longer ramps
% modSel = [1     2     3     4     5     6    12    14    15    16]; 
% mod = [0.1360    4.9380    0.4368    1.0218    1.6463    6.4469    1.0000    1.0000    1.0000    1.0000    1.0000    2.8594    1.0000   13.8852   1.0412    1.0079    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000];
% throwing in couple more params, still pCa 10 in log space from ramp begin
% modSel = [1     2     3     4     5     6    10    11    12    13    14   15    16];
% testing pCa
% mod([20 21]) = NaN; modSel = [7 8 9 22];
% mod([16, 17]) = NaN;
% solo pCa 
modSel = [7     8     9    22]; mod = [0.1035    6.5217    0.3626    0.8185    1.7903    4.7326    2.1097    1.7771    0.3800    1.3824    1.0000    1.0000    1.0000   84.2683    1.0000    1.0000    1.0000    1.0000    1.0000       NaN       NaN   19.3337];
% combined
modSel = [1     2     3     4     5     6     7     8     9    22]; mod = [0.0780    2.7953    0.4019    0.9791    1.5029    4.4448    2.6626    2.1066    0.2870    1.3824    1.0000    1.0000    1.0000   84.2683   1.0000    1.0000    1.0000    1.0000    1.0000       NaN       NaN    8.3950];
% actually pretty good combined. Final candidate?
mod = [0.0921    2.2736    0.3899    0.9107    1.4463    6.2490    2.6767    2.3335    0.3324    1.3824    1.0000    1.0000    1.0000   84.2683    1.0000    1.0000    1.0000    1.0000    1.0000       NaN       NaN    6.4155 0];
% reopt for Ca only, incl. the Ca dependency
% modSel = [7     8     9    22 23]; mod = [0.0921    2.2736    0.3899    0.9107    1.4463    6.2490    3.5190    1.7304    0.3318    1.3824    1.0000    1.0000    1.0000   84.2683    1.0000    1.0000    1.0000    1.0000    1.0000       NaN       NaN    7.6587    0.8273];
tic
% knock off of PEVK, reoptim. alphaU did not help here
% mod([7 8 9 22]) = [1e-6 1e6 0.3282   11.2732];
% add alphaU, reoptim it all
modSel = [1     2     3     4     5     6     9    20 22];
% % mod([7 8 9 20 22]) = [1e-6 1e6 0.3282   mod(20)/0.7 11.2732];
% mod(modSel) = [0.1106    2.0271    0.4126    0.9117    1.4129    6.7698    0.3801   14.6443    9.5714];

% no unfolding of bound states costs 39.85
% mod = [0.0944    2.5913    0.4081     0.9095    1.4934    2.6794   90.0356  373.9408    0.2881    1.3824    1.0000    1.0000    1.0000   89.8580    1.0000    1.0000    1.0000    1.0000    1.0000       NaN       NaN   11.5157         0]
% mod(12) = 1;mod(14) = 1;
% mod(10) = 0.1;mod(15) = 10;
% get init from Dan's codes
% modSel = [7 8 9 22]; mod(modSel) = [5 1 4 4].*mod([16 17 1 2]);
% mod([20 21]) = NaN;
% figure(144)
% mod(23 - length(mod)) = NaN;
mod(15) = 0;
% mod(3) = 0.5;
% tic
% mod([7 8]) = [1e-6 1e6];
% mod([9 22]) =  [ 0.3282      11.2732];
% mod([7 8]) = [100/200 500/50];
% mod(7) = 40;
% mod(8) = 1000;
% mod(20) = NaN;
mod2param = [600, 500, 2, 6, 2, 1e6, 200, 50, 600, 3.2842, 0.01, 8, 0.01, 0.2, 1e-3, 200, 50, 0.9, 0.01*15, (8.4137e5)*0.7, 3.2470, 500, 0.1];
% create abs params and fix Lref in Dan's codes
params = mod.*mod2param; params([1 2 6]) = params([1 2 6]).*[(Lref^params(3)), (Lref^params(5)), (Lref^params(4))];
% Nice fit including pCa4.4 0.1s, predicting rest
params = [367, 3e+04, 2.3, 9, 2.33, 3.24e+06, 4.98, 84.9, 1.36e+03, 4.89, 1.01e-08, 12.8, 0.00389, 0.678, 0, NaN, NaN, 0.9, 0.175, NaN, NaN, 3.94e+04, 0, ];

convertLref = ones(size(params));
% convert alphaU and kd and kp, to Lref = L_0 = 1 mum
% kp = 1 9 np = 3, kd = 2 22 nd = 5, alphaU = 6 20 nu = 4
kpcorr = 1/(params(18)^params(3)); kdcorr = 1/(params(18)^params(5)); alphaUCorr = 1/(params(18)^params(4));
convertLref([1 2 6 9 18 20 22]) = [kpcorr kdcorr alphaUCorr kpcorr 1/params(18) alphaUCorr kdcorr];
params = params.*convertLref;

tic
mod(14) = 0.1;
evalCombined(params(modSel), params, modSel, [11])
toc
% modSel = [1:6 7 8 9 14 20 22]

%% 
tic
% modSel = [1:6    11    12    14]; params = [367, 2.39e+04, 2.3, 9, 2.33, 3.24e+06, 1.4, 17.8, 4.44e+03, 4.89, 1.2e-10, 17.6, 0.0027, 0.792, 0, NaN, NaN, 0.9, 0.175, NaN, NaN, 6.96e+03, 0, ];
% modSel = [11 12 13 14];
% params([7 8 9 22 23]) = [5 60 1380 3e4 0];
% modSel = [7 8 22];
m = params;
m(14) = 1;
evalCombined(m(modSel), m, modSel, [11])
% evalCombined(params(modSel), params, modSel, [11])
toc
%%
for i = 1:length(params)
    fprintf('%d: %s = %f\n', i, modNames{i}, params(i))
end
%%
clf;
modNames{16} = 'k_{PEVK,A} (low Ca)';
modNames{17} = 'k_{PEVK,D} (low Ca)';
mm = min(floor(log10(mod)));

n = length(mod);
polarplot(linspace(2*pi/n, 2*pi, n), log10(mod) - mm, 'x-', LineWidth=2);
% evalCombined(mod, mod, 1:n);
% set(gca, 'ThetaAxisUnits', 'radians', 'thetatick', linspace(2*pi/n, 2*pi, n));
set(gca, 'ThetaAxisUnits', 'radians', ...
    'thetatick', linspace(2*pi/n, 2*pi, n), 'thetaticklabel', modNames, ...
    'Rlim', [-3 1]-mm, 'RTick', [-2 -1 0 1]-mm, 'RTickLabel', [0.01 0.1 1 10]);
hold on;

%%
% mod = [0.0231    0.2275    0.4870    1.0234    0.7909   0.2929    0.1113    0.2652    0.0218    0.8988    0.7426    1.8401    0.7510    1.2811    1.6663];
% mod = ones(1, 15);
% modSel = [1 2 3 5 6 11 12 15];
% modSel = [7, 8, 9];
% modSel = [1 2 3 5 6 7 8 9 11 12 15];
% 
% modSel = [1 2 3 5 6 10 11 12 13 15];
% modSel = [1 2 3 5 6 8 10 11 12 13 15];
% mod = [    0.0140    0.3630    0.4241    1.0234    0.7953    0.2941    0.2514     0.8819    0.0135    0.8988    0.7571    4.3872    0.7510    1.2811    1.6652];
% cost = 301.7
options = optimset('Display','iter', 'TolFun', 1e-3, 'Algorithm','sqp', 'TolX', 0.01, 'PlotFcns', @optimplotfval, 'MaxIter', 500);
% x0 = x0([1:4 6:10 13]);
% x0 = [0.5795    1.6029    0.9047    1.0569    0.9236    0.6627    1.0497    0.9465    1.0641    0.7124];

% best in the worst so far optMods =     1.3705    1.4770    1.1905    1.1205    0.8670    0.3074    0.6812    0.9407    1.3393    0.8915    1.1628    1.0166
% modSel = [2 3 4 6 7 8 9 12 14 15];

% init = mod(modSel);
% init = ones(1, 10);
% modSel = [7, 9];
% mod =  [0.0299    0.8084    0.4798    1.0000    1.0862    1.6973    1.0000    1.0000    1.0000    1.2983    1.0000    1.0000    1.0000    0.1    0.50000];
% reducing the mu to nonsensitive value, 
% modSel = [1 2 3 5 6 10]; mod =  [0.0293    0.8574    0.4871    1.0000    0.9528    1.7817    1.0000    1.0000    1.0000    1.3369    1.0000    1.0000    1.0000    0.1000    0.5000];
% optimizing for the ramp-ups in logspace
% modSel = [1 2 3 4 5 6 7 8 9 10];mod = [0.4804    0.1794    0.9727    1.6729 0.7547   1.1479e+03 0.0245    0.1301    0.6788    1.4584 3.6372    0.2425    0.0030    0.1000 0.5000];
% optimizing for high pCa
modSel = [2 3 4 5 6 7 8 9 10 11 12 13]; mod = [0.4804    0.1794    0.9727    1.6729 0.7547   1.1479e+03 0.0245    0.1301    0.6788    1.4584 3.6372    0.2425    0.0030    0.1000 0.5000];


%%
modSel = [1:10 14:19];
modSel = 1:19;
modSel = [1:15 20];
% got as a 1e-3 cutoff on prelim run sensitivity
modSel = 0 < [1   0   0   0   0   1   1   1 1   1   0   1   1   1   1   0 0   0   0   0   1];
mod = m;
mod(20) = NaN;
%%
% mod(12) = 1;mod(14) = 10;
% mod(15) = 1;
% modSel = [1:6 10 15];
% modSel = [1:6 12 14 15 16]
% modSel = [1:6 7 8 9 22];
init = params(modSel);
% evalFunc = @(optMods) evalCombined(optMods, mod, modSel);

evalLin = @(optMods) evalCombined(optMods, params, modSel, [11])
x = fminsearch(evalLin, init, options);
params(modSel) = x;
%%
% mod(modSel) = x;
% mod
% save("mod.mat", "mod")
% modSel = [1     2     3     4     5     6     7     8     9    22 23];
init = params(modSel);
% evalFunc = @(optMods) evalCombined(optMods, mod, modSel);
evalLin = @(optParams) evalCombined(optParams, params, modSel, [4.4])
x = fminsearch(evalLin, init, options);
%% optim in log param space
% mod = ones(1, 17);
% mod = [mod ones(1, 21 - length(mod))];
% mod(21) = mod(10)
% modSel = [1:10 14:17];
% modSel = [11 12 13];
% modSel = 1:17;
% modSel = [1:17 21];

% Ca sensitive only
% modSel = [7 8 9 20];
% pca 5.75
% mod(modSel) = [0.1052    0.7715    0.2400  590.8530];
% pca 6
% mod(modSel) = [0.0274    1.0723    0.1016 2.1787e+03];
% 4.3e3
% mod = [0.1644    0.0822    0.9186    1.2568    0.4835    1.1405    0.1305    0.8418    0.2723    1.6966         0    1.4225    0.0027    0.5050    2.5843       NaN       NaN    1.0000 1.0000       NaN       NaN    0.9618];
% for all ramps
% mod = [0.1608    0.0820    0.9177    1.2529    0.4815    1.1445    0.1302    0.8398    0.2715    1.6819         0    1.4225    0.0027    0.5112    2.5843       NaN       NaN    1.0000 1.0000       NaN       NaN    0.9612];
% using prescribed Fss_noCa and FssCa
% mod = [0.3262    0.0269    0.8916    1.1183    0.4793   1.1058    6.1801    0.7867    0.1298    1.5060         0    1.4225    0.0027    0.5112    2.5843       NaN       NaN    1.0000    1.0000       NaN    4.1961    0.9254];
% updated the data to 60s decay
% mod = [0.1877    0.1925    0.9268    1.3007    0.5661    1.3564    0.1461    1.2518    0.2298    1.7830         0    1.4225    0.0027    0.7030    2.5843       NaN       NaN    1.0000    1.0000       NaN       NaN    0.8902];

% modSel = [1:10 12 14 15 22];
% fixed bug
% mod = [0.1460, 0.4296, 0.8464, 1.3103, 0.6893, 3.1370, 0.0388, 1.1946, 0.1355, 1.7664,      0, 1.4225, 0.0027, 0.6908, 2.5843,    NaN,    NaN, 1.0000, 1.0000,    NaN,    NaN, 0.6207];

% optim all ramps
% mod = [0.1360, 0.4811, 0.8499, 1.3379, 0.6881, 2.9849, 0.0849, 1.1904, 0.1302, 1.6995,      0, 1.4243, 0.0027, 0.7668, 2.4840,    NaN,    NaN, 1.0000, 1.0000,    NaN,    NaN, 0.4930]

% pCa 11 only
% mod = [1.3111    0.0572    1.0882    0.9839    0.5485    0.0672    0.0849    1.1904    0.1302    1.5060         0    3.1310    0.0027    0.3576    2.4840       NaN       NaN    1.0000    1.0000       NaN       NaN    0.4930];
% modSel = [1:6 11 12 14];

% optim for tail only: pCa11 0.1s ramp, 300s linear sampling
% mod = [1.8533    0.0821    1.1321    0.9230    0.6111    0.0301    0.0849    1.1904    0.1302    1.5060         0    2.8146    0.0027    0.8606    2.4840       NaN       NaN    1.0000    0.9400       NaN       NaN    0.4930];
init = max(-10, log10(params(modSel)));
evalLogCombined = @(logMod) evalCombined(10.^logMod, params, modSel, [4.4]);
x = fminsearch(evalLogCombined, init, options);
params(modSel) = 10.^x; 
% list params
a = [modSel; params(modSel)]; sprintf('%d: %1.3g\n', a(:))
% list params for save
disp(['params = [' sprintf('%1.3g, ', params(:)) '];'])


% params = [367, 2.39e+04, 2.3, 9, 2.33, 3.24e+06, 1.4, 17.8, 4.44e+03, 4.89, 1.01e-08, 12.8, 0.00389, 0.678, 0, NaN, NaN, 0.9, 0.175, NaN, NaN, 6.96e+03, 0, ];
% modSel = [7 8 9 22]; params([7 8 9 22]) = [10 50 params([1 2])];params = [367, 3e+04, 2.3, 9, 2.33, 3.24e+06, 7.77, 47.2, 1.4e+03, 4.89, 1.01e-08, 12.8, 0.00389, 0.678, 0, NaN, NaN, 0.9, 0.175, NaN, NaN, 1.42e+04, 0, ];
% params = [367, 3e+04, 2.3, 9, 2.33, 3.24e+06, 0.698, 7.38, 1.33e+03, 4.89, 1.01e-08, 12.8, 0.00389, 0.678, 0, NaN, NaN, 0.9, 0.175, NaN, NaN, 6.9e+03, 0, ];

%% optim with bounded params
options = optimset('Display','iter', 'TolFun', 1e-3, 'Algorithm','sqp', 'UseParallel', true, ...
    'TolX', 0.01, 'PlotFcns', @optimplotfval, 'MaxIter', 500);

% modSel = [11 12 13];
init = log10(mod(modSel));
evalLogCombined = @(logMod) evalCombined(10.^logMod, mod, modSel);
[x, FV, EF, OUT] = fmincon(evalLogCombined, init, [], [], [], [], log10(lb(modSel)), log(ub(modSel)), [], options);
%% test in GA
% parpool
ga_Opts = optimoptions('ga', ...
    'PopulationSize',64, ...            % 250
    'Display','iter', ...
    'MaxStallGenerations',8, ...  % 10
    'UseParallel',true);
modSel = 1:14;
% Ng = length(modSel);
ub = 100*ones(1, length(mod));
ub([3, 4, 5, 10, 14]) = [10 10 10, 10, 10];
lb = 0.001*ones(1, length(mod));
lb([3, 4, 5, 10, 14]) = [.1 .1 .1 .1 .1];
evalLogCombined = @(logMod) evalCombined(10.^logMod, mod, modSel);
init = log10(mod(modSel));

[p_OptimGA,Res_OptimGA,~,~,FinPopGA,FinScoreGA] = ...
    ga(evalLogCombined,length(modSel), ...
    [],[],[],[],...
    log10(lb(modSel)),log10(ub(modSel)),[],ga_Opts);

mod(modSel) = 10.^p_OptimGA;
% use fminserach afterwards
%% result
mod = [2.0398    0.9359    4.3424    2.3068    0.4784 0.6410    0.8054    0.8220    0.2923    0.7200 1.3634    1.0000   -2.8000    1.3022    0.8551];
mod = [2.0398    0.9359    4.3397    2.2756    0.4784 0.6412    0.8053    0.8333    0.2923    0.7200 1.3634    1.0468   -2.8000    1.3200    0.8759];

% proximal chain Ca dependency
mod = [1.1592    1.0379    0.9763    0.9779    1.1237    1.0935    0.9365    1.0882    0.9846    0.8931];
% Fix for nS
mod = [1.1697    1.0418    0.9774    0.9737    0.9858    1.0265    0.9403    1.0837    0.9889    0.8988 1 1 1];

% mods for pCa 6: modSel = [7, 9], x = [0.1657, 0.3895];
mod([7,9]) = [0.1657, 0.3895];
% thus, the profile is:
pCa = [11, 6, 4];
kC   = [10203,10203*4.78*0.3895,10203*4.78*0.9889];      % proximal chain force constant
kA   = [0, 16.44*0.1657,16.44*0.9403];
figure(44);clf;
% subplot(211);
plot(pCa, kA, 'x:', LineWidth=2);
ylabel('kA Value');
yyaxis right;
plot(pCa, kC, '+:', LineWidth=2);
ylabel('kC Value');

set ( gca, 'xdir', 'reverse' );
xlabel('pCa');
legend('kA (PEVK attachment)', 'kC (stiffening proximal chain)')
%% test pCa
figure(90);clf;hold on;
isolateRunCombinedModel(mod, 5.5, true)
figure(91);clf;hold on;
isolateRunCombinedModel(mod, 5.75, true)
figure(92);clf;hold on;
isolateRunCombinedModel(mod, 6, true)
figure(93);clf;hold on;
isolateRunCombinedModel(mod, 6.25, true)

%% Optimizing the reduced set to pCa
options = optimset('Display','iter', 'TolFun', 1e-6, 'Algorithm','sqp', 'UseParallel', false, ...
    'TolX', 1e-4, 'PlotFcns', @optimplotfval, 'MaxIter', 50);

% mod = [0.2880    0.0763    0.9277    1.7014    0.5508  632.7487    0.2108    0.5573    0.2802    1.5806    1.0005    1.3387    0.1604    0.5548    1.1785    NaN    NaN    1.0000    1.0000  415.4241];
mod5_5 = params;mod5_8 = params;mod6=params;
modSel = [7  9 22];
modNames(modSel);


% lower and upper bounds are a functions of previous mods
ub = @(m)[m(7) m(9) m(22)];
lb = @(m)[1e-6 m(1) m(2)];
m_ub = nan(size(mod));m_ub(modSel) = ub(params);
m_lb = nan(size(mod));m_lb(modSel) = lb(params);
clf;
plotParams(params, [-1 3]);hold on;
plotParams(m_ub, [-1 3]);hold on;
plotParams(m_lb, [-1 3]);
%%
% evalNewCombined = @(curmod) evalCombined(curmod, mod, 1:length(mod), [5.5]);
% modt = mod5_5;
% modt(8) = mod(8);
% evalNewCombined(modt)
%%
tic
% 2.3214
% evalLin([4.9816 1.3607e3 3.9386e4])
%
% mod5_5(modSel) = [1.1816 1.2607e3 3.9386e4];
mod5_5(modSel) = [2.4850 1.2607e+03 3.9386e+04];
evalLin([3 370.4209 2.028e+04])
%%
modSel = [7 8 9];
evalCombFixRat = @(params) evalCombined([params(1:21) params(2)*params(9)/params(1) params(23)], params, 1:23, [4.4]);
evalLin = @(curmod) evalCombFixRat([params(1:6) curmod params(10:end)]);
x = fminsearch(evalLin, params(modSel), options);
%%
tic
mod = mod5_5;
mod(22) = NaN; % auto settings
evalCombined(mod, mod, 1:23, [5.5]);
% evalLin(params(modSel))

toc
%%
% save('ModsNoAss.mat', 'mod', 'mod6','mod5_8','mod5_5')

% init = max(-10, log10(mod5_5(modSel)));
% evalLogCombined = @(logMod) evalCombined(10.^logMod, mod, modSel, [5.5]);
params(22) = NaN;
mod5_5 = params;mod5_8 = params;mod6=params;
modSel = [7 9];

% x = fminsearch(evalLogCombined, init, options);
evalLin = @(curmod) evalCombined([curmod params(2)*params(9)/params(1)], mod5_5, [modSel 22], [5.5]);
% x = fmincon(evalLin, mod5_5(modSel), [], [], [], [], lb(params), ub(params), [], options);
% x = fmincon(evalLin, mod5_5(modSel), [], [], [], [], [], [], [], options);
x = fminsearch(evalLin, mod5_5(modSel), options);
% mod5_5 = mod;
mod5_5(modSel) = x;
%%
% init = max(-10, log10(mod5_8(modSel)));
% evalLogCombined = @(logMod) evalCombined(10.^logMod, params, modSel, [5.75]);
% x = fminsearch(evalLogCombined, init, options);
% x = fmincon(evalLogCombined, init, [], [], [], [], log10(lb(mod5_5)), log10(ub(mod5_5)), [], options);
% mod5_8(modSel) = 10.^x;
evalLin = @(curmod) evalCombined(curmod, mod5_5, modSel, [5.75]);
x = fminsearch(evalLin, mod5_5(modSel), options);
mod5_8(modSel) = x;

%%
% init = max(-10, log10(mod6(modSel)));
% evalLogCombined = @(logMod) evalCombined(10.^logMod, params, modSel, [6]);
% x = fminsearch(evalLogCombined, init, options);
% x = fmincon(evalLogCombined, init, [], [], [], [], log10(lb(mod5_8)), log10(ub(mod5_8)), [], options);
% mod6(modSel) = 10.^x;

evalLin = @(curmod) evalCombined(curmod, mod5_8, modSel, [6]);
x = fminsearch(evalLin, mod5_8(modSel), options);
mod6(modSel) = x;
%%
% compare 
figure(39);
clf;
sclale = [-1 6];
plotParams(params, sclale)
hold on;
plotParams(mod5_5, sclale);
plotParams(mod5_8, sclale);
plotParams(mod6, sclale);
%
% evalCombined(mod6(modSel), mod, modSel, [5.5])
%%
mod = params;
mod11 = params;
mod11(modSel) = [1e-6 mod(1)]
modSet = [mod;mod5_5;mod5_8;mod6;mod11];
pcax = [4.4, 5.5, 5.75, 6, 11];
% modSel = 1:20;
% modSet(:, modSel)
figure(40);clf;
tiledlayout(1,4);
for i = 1:size(modSel, 2)
    nexttile();
    semilogy(-pcax, modSet(:, modSel(i)), 'x-');

    
    % pca11
    % semilogy(-pcax, modSet(:, modSel(i)), '-');

    legend(modNames{modSel(i)});
end

disp(['mod5_5 = [' sprintf('%1.3g, ', mod5_5(modSel)) '];'])
disp(['mod5_8 = [' sprintf('%1.3g, ', mod5_8(modSel)) '];'])
disp(['mod6 = [' sprintf('%1.3g, ', mod6(modSel)) '];'])
%% save last iteration
params(22) = NaN;
modSel = [7 9];
mod5_5 = params;mod5_8 = params;mod6 = params;
mod5_5(modSel) = [0.000165, 1.2e+03];
mod5_8(modSel) = [0.000188, 1.04e+03];
mod6(modSel) = [0.000238, 396];
%% Running all
modSel = 1:23;
evalCombined(mod5_5, params, modSel, [5.5]);
evalCombined(mod5_8, params, modSel, [5.75]);
evalCombined(mod6, params, modSel, [6]);

%%
function totalCost = evalCombined(optMods, mod, modSel, pCas)

    if nargin < 4
        % evaluate pCa 4.4 and 10 by default
        pCas = [4.4 10];
    end
    %normal - optimizing for all
    % modSel = 1:15;

    % modSel = [1 2 3 5 6 10];
    % optimizing only subset of mods
    % mod = [1.1697    1.0418    0.9774    0.9737    0.9858    1.0265    0.9403    1.0837    0.9889    0.8988 1 1 1];
    % mod = [1.16970000000000	0.928400000000000	0.977400000000000	1.02340000000000	1.01370000000000	1.10320000000000	0.937900000000000	1.19500000000000	0.909900000000000	0.898800000000000	1	1	1 1 1];
    % mod = [0.0165    0.7035    0.4403    1.0234    1.0077    0.5754    0.9379    1.1950    0.9099    0.8988    1.0000    1.1421    1.4792    1.1156    2.9834];
    % mod = [0.0185    0.8479    0.4307    1.0234    1.0326 0.5971    0.9379    1.1950    0.9099    0.8988 1.0000    1.4450    0.7510    1.2811    2.7365];
    % for absolute average
    % mod = [0.0343    0.7101    0.4731    1.0234    1.0916    1.9353    0.9379 1.1950    0.9099    0.8988    0.5952    2.0416    0.7510    1.2811 4.1891];
    % optim for -log10 weighting
    % mod = [0.928351480405723	0.928351480405723	1.01367539052550	1.02336567490158	1.01367539052550	1.10319213342611	0.937882365838957	1.19500150970587	0.909890571615859 1 1 1 1];
    % reset the search
    % mod = ones(1, 15);
    % pCa 11 decay with reduced mu
    % mod =  [0.0299    0.8084    0.4798    1.0000    1.0862    1.6973    1.0000    1.0000    1.0000    1.2983    1.0000    1.0000    1.0000    0.1    0.50000];
    % modSel = [1 2 3 5 6 10];
    

    % % mod([1:4 6:10 13]) = optMods;
    % modSel = [11, 12, 13];
    % modSel = [1 2 3 5 6 11 12 15];
    % modSel = [7, 8, 9];
    % modSel = [1 2 3 5 6 10];
    
    % store the init
    baseMods = mod;
    mod(modSel) = optMods;

    drawPlots = true;
    % drawPlots = false;
    totalCost = 0;

    if drawPlots && 1
        figInd = 100;
        try 
            set(groot,'CurrentFigure',figInd); % replace figure(indFig) without stealing the focus
        catch 
            f = figure(figInd); % if indFig is not an exsting figure it creates it (and steal the focus)
            f.Name = 'Parameter modifiers';
        end      
        clf;
        % plotParams(baseMods, [-1, 4]);hold on;
        plotParams(mod./baseMods, [-1, 1]);hold off;
    end


    %% pCa 4
    if ismember(4.4, pCas)
        pCa = 4.4;
        plotOnBackground(drawPlots, pCa);
        % RunCombinedModel;
        cost = isolateRunCombinedModel(mod, pCa, drawPlots);
        totalCost = totalCost + cost;
    end
    %% pCa 5.5
    if ismember(5.5, pCas)
        pCa = 5.5;
        plotOnBackground(drawPlots, pCa);
        % RunCombinedModel;
        cost = isolateRunCombinedModel(mod, pCa, drawPlots);
        totalCost = totalCost + cost;
    end    
    %% pCa 5.75
    if ismember(5.75, pCas)
        pCa = 5.75;
        plotOnBackground(drawPlots, pCa);
        % RunCombinedModel;
        cost = isolateRunCombinedModel(mod, pCa, drawPlots);
        totalCost = totalCost + cost;
    end

    %% pCa 6
    if ismember(6, pCas)
        pCa = 6;
        plotOnBackground(drawPlots, pCa);
        % RunCombinedModel;
        cost = isolateRunCombinedModel(mod, pCa, drawPlots);
        totalCost = totalCost + cost;
    end

    %% no Ca, but evaluate PEVK binding
    if ismember(10, pCas)    
        pCa = 10;
        plotOnBackground(drawPlots, pCa);
        % RunCombinedModel;
        cost = isolateRunCombinedModel(mod, pCa, drawPlots);
        totalCost = totalCost + cost*10;
    end
    %% no Ca, PEVK binding not evaluated
    if ismember(11, pCas)    
        pCa = 11;
        plotOnBackground(drawPlots, pCa);
        % RunCombinedModel;
        cost = isolateRunCombinedModel(mod, pCa, drawPlots);
        totalCost = totalCost + cost*10;
    end    
% return
end 

function plotOnBackground(drawPlots, pCa, cols)
    if nargin < 3
        cols = 2;
    end
    figInd = 100 + round(pCa*10);
    if drawPlots
        try 
            set(groot,'CurrentFigure',figInd); % replace figure(indFig) without stealing the focus
        catch 
            f = figure(figInd); % if indFig is not an exsting figure it creates it (and steal the focus)
            f.Name = ['pCa ' num2str(pCa)];
            if cols == 1
                aspect = 2;
                % normal size of 2-col figure on page is 7.2 inches
                % matlab's pixel is 1/96 of an inch
                f.Position = [300 200 3.5*96 3.5*96/aspect];
            else
                aspect = 1.5;
                f.Position = [300 200 7.2*96 7.2*96/aspect];
            end
        end
    end
end


function cost = isolateRunCombinedModel(mod, pCa, drawPlots)
% just to isolate the script, so the variables can't intervene
    % drawPlots = true;
    RunCombinedModel;
end

function plotParams(mod, mm, resetGca)

    % modNames = {'k_p(NoCa)', 'k_d', 'n_p', 'n_U', 'n_d', 'alphaU', 'k_{PEVK,A}', 'k_{PEVK,D}', 'k_p(highCa)', 'Fss', 'b', 'c', 'd', 'mu', 'alphaF_0','k_{PEVK,A} (low Ca)', 'k_{PEVK,D} (low Ca)', 'Lref', 'delU'};
    modNames = {'k_p(NoCa)', 'k_d', 'n_p', 'n_U', 'n_d', 'alphaU', 'k_{PEVK,A}', 'k_{PEVK,D}', 'k_p(highCa)', 'Fss', 'b', 'c', 'd', 'mu', 'alphaF_0','k_{PEVK,A} (low Ca)', 'k_{PEVK,D} (low Ca)', 'Lref', 'delU', ...
        'AlphaU_pCa', 'Fss_pCa', 'kd_pCa', 'kDf'};

    indxs = strcat(string(1:length(modNames)), ':');
    modNames = strcat(indxs, modNames);

    if nargin < 2 || isempty(mm)
        % try
        %     mima = get(gca, 'RLim');
        %     mi = mima(1);ma=mima(2);
        % catch
            mi = min(floor(log10(mod)));
            ma = mi + 3;
        % end
    else
        mi = mm(1);
        ma = mm(2);
    end
    n = length(mod);
    
    polarplot(linspace(2*pi/n, 2*pi, n), max(1e-3, log10(mod) - mi), 'x-', LineWidth=2);
    
    if ~strcmp('degrees', get(gca, 'ThetaAxisUnits')) || (nargin >= 3 && resetGca)
        % degrees means it is not adjusted yet 
        % or the param has not been provided 
        % or we want to reset
        return
    end
    rng = mi:ma;
    set(gca, 'ThetaAxisUnits', 'radians', ...
        'thetatick', linspace(2*pi/n, 2*pi, n), 'thetaticklabel', modNames, ...
        'Rlim', [min(rng) max(rng)]-mi, 'RTick', rng - mi, 'RTickLabel', 10.^rng);
end

function [c0, cost_sap, cost_sam] = runSa(pCa, mod, saSet)
    drawPlots = false;
    cost_sap = []; % SA plus
    cost_sam = []; % SA minus

    % pCa = 4;
    if drawPlots
        figure(100);
    end
    cost = isolateRunCombinedModel(mod, pCa, drawPlots);
    c0 = cost;
    %
    % saSet = [1:9 12:15];
    % saSet = [1:8 13];
    % saSet = [14 15];
    
    SAFact = 1.05;
    for i_m = saSet
        mod(i_m) = mod(i_m)*SAFact;
        fprintf('Mod %g is up to %g..', i_m, mod(i_m));
        % figure(i_m)
        cost = isolateRunCombinedModel(mod, pCa, drawPlots);
        cost_sap(i_m) = cost;
        
        mod(i_m) = mod(i_m)/SAFact/SAFact;
        fprintf('costing %1.4e€ and down to %g...', cost, mod(i_m));
        cost = isolateRunCombinedModel(mod, pCa, drawPlots);
        cost_sam(i_m) = cost;
        mod(i_m) = mod(i_m)*SAFact;
        fprintf('costing %1.4e€. \n', cost);
    end
end