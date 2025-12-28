% --- START OF FILE main_olg_v4.m ---

% Huggett (1996) 复制，包含动态人口结构
% 修改版：采用“混合时间单位”，包含个人养老金计划 (PPS),
% PAYG工资税率外生 (theta), PAYG福利内生 (b_payg)
% V4: Incorporates Heer (2020) style government finance (tau_k, tau_c, tau_l, G, B)
%     and exogenous PAYG tax rate theta.

clc; clear; close all;
fprintf('=== Huggett 模型: 动态人口, PPS, 外生 PAYG 税率, Heer GBC ===\n');

%% 1. 固定模型参数与人口设置
fprintf('\n--- 1. 初始化参数 ---\n');
cS = main_olg_v4_utils_pps.ParameterValues_HuggettStyle(); % cS: 包含所有常数和参数的结构体
paramS = struct(); % paramS: 包含派生参数的结构体

fprintf('参数已加载。nk=%d (非PPS资产网格点数), nw=%d (劳动效率状态数)。\n', cS.nk, cS.nw);
fprintf('年度年龄范围: %d-%d (%d 年)。年龄组数: %d (%d 个工作年龄组)。\n', ...
        cS.age1_orig, cS.ageLast_orig, cS.aD_orig, cS.aD_new, cS.aR_new);
fprintf('VFI 中使用年度 beta = %.4f。外生 PAYG 税率 (theta): %.3f\n', cS.beta, cS.theta);
fprintf('税率: tau_k=%.2f, tau_c=%.2f, tau_l=%.2f\n', cS.tau_k, cS.tau_c, cS.tau_l);
fprintf('政府财政: G/Y=%.2f, B/Y=%.2f\n', cS.gov_exp_frac_Y, cS.gov_debt_frac_Y);
if cS.pps_active
    fprintf('个人养老金计划 (PPS) 已激活。\n');
    fprintf('  最大缴费比例: %.4f, 领取期税率: %.2f, 回报率溢价: %.3f, 领取率: %.2f\n', ...
        cS.pps_max_contrib_frac, cS.pps_tax_rate_withdrawal, cS.pps_return_rate_premium, cS.pps_withdrawal_rate);
    fprintf('  PPS 计入总资本 K: %d, PPS 可遗赠: %d\n', cS.pps_in_K, cS.pps_bequeathable);
else
    fprintf('个人养老金计划 (PPS) 未激活。\n');
end

%% 2. 模拟人口动态至稳态 (年龄组层面)
fprintf('\n--- 2. 模拟人口动态 (年龄组层面) ---\n');
popS = main_olg_v4_utils_pps.initPopulation(cS);
popS = main_olg_v4_utils_pps.populationDynamics(popS, cS);
[Z_ss, dep_ratio_ss, bgp_reached, bgp_period] = main_olg_v4_utils_pps.detectSteadyStatePopulation(popS, cS);

fprintf('\n--- 人口模拟总结 ---\n');
fprintf('实际模拟期数: %d\n', length(popS.totalPop)-1);
if bgp_reached, fprintf('人口稳态 (年龄组) 在第 %d 期达到。\n', bgp_period);
else, fprintf('人口稳态 (年龄组) 未达到。使用第 %d 期的最终数据。\n', bgp_period); end
fprintf('稳态抚养比 (退休人口数/工作人口数): %.4f\n', dep_ratio_ss);

paramS.Z_ss_counts = Z_ss; % Absolute counts for each age GROUP in steady state

Z_ss_total = sum(Z_ss);
Z_ss_norm_group = zeros(cS.aD_new, 1);
if Z_ss_total > 1e-9, Z_ss_norm_group = Z_ss / Z_ss_total; end
paramS.ageMassV = Z_ss_norm_group(:); % Normalized steady-state age GROUP distribution (sums to 1)

Z_ss_norm_annual = zeros(cS.aD_orig, 1);
if Z_ss_total > 1e-9
    for a_new_map_idx = 1:cS.aD_new % Corrected loop variable name
        annual_indices = cS.physAgeMap{a_new_map_idx};
        group_mass = Z_ss_norm_group(a_new_map_idx);
        num_years_in_group = length(annual_indices);
        if num_years_in_group > 0
            mass_per_year = group_mass / num_years_in_group;
            Z_ss_norm_annual(annual_indices) = mass_per_year;
        end
    end
    if sum(Z_ss_norm_annual) > 1e-9 && abs(sum(Z_ss_norm_annual) - 1.0) > 1e-6 % Check if it sums to 1
        Z_ss_norm_annual = Z_ss_norm_annual / sum(Z_ss_norm_annual); % Re-normalize if not
    elseif sum(Z_ss_norm_annual) < 1e-9
        Z_ss_norm_annual(:) = 1/cS.aD_orig; % Fallback if sum is zero
    end
    fprintf('已推导出近似的年度稳态人口分布 (归一化)。\n');
else
    warning('稳态人口为零，无法推导年度分布。');
    Z_ss_norm_annual(:) = 1/cS.aD_orig;
end
paramS.Z_ss_norm_annual = Z_ss_norm_annual; % Normalized steady-state ANNUAL distribution (sums to 1)

if Z_ss_total > 1e-9 && length(popS.totalPop) > 1 && popS.totalPop(end-1) > 1e-9
    % popS.totalPop is history of Z_ss_total (sum of absolute counts per group)
    % popS.totalPop(end) / popS.totalPop(end-1) is growth factor over cS.yearStep years
    pop_growth_factor_per_group_period = popS.totalPop(end)/popS.totalPop(end-1);
    paramS.popGrowth_ss_group_period = pop_growth_factor_per_group_period - 1;
    paramS.popGrowthForDebt = pop_growth_factor_per_group_period^(1/cS.yearStep) - 1; % Annualized
    fprintf('稳态人口增长率 (来自模拟, %d年期): %.4f (年化用于债务: %.4f)\n', cS.yearStep, paramS.popGrowth_ss_group_period, paramS.popGrowthForDebt);
else
    warning('无法从人口模拟计算稳态增长率，使用 cS.popGrowth_orig = %.4f for debt dynamics', cS.popGrowth_orig);
    paramS.popGrowthForDebt = cS.popGrowth_orig; % Annual BGP growth rate for debt
end


%% 3. 预计算劳动供给和禀赋过程
fprintf('\n--- 3. 预计算劳动 ---\n');
[paramS.leLogGridV, paramS.leTrProbM, paramS.leProb1V] = main_olg_v4_utils_pps.EarningProcess_olgm(cS);
paramS.leGridV = exp(paramS.leLogGridV(:));
fprintf('劳动禀赋过程已校准 (nw=%d 个状态)。\n', cS.nw);
paramS.ageEffV_orig = cS.ageEffV_orig; paramS.ageEffV_new = cS.ageEffV_new;
fprintf('为 %d 个个体模拟 %d 年的年度劳动禀赋...\n', cS.nSim, cS.aD_orig);
eIdxM = main_olg_v4_utils_pps.LaborEndowSimulation_olgm(cS, paramS);
% L is total efficiency labor supply per capita of total population
[~, L] = main_olg_v4_utils_pps.LaborSupply_Huggett(eIdxM, cS, paramS, paramS.ageMassV); % Use paramS.ageMassV (normalized group dist)
fprintf('总劳动供给 (L_ss, 效率单位, 总体人均): %.4f\n', L);
if L <= 0 && sum(paramS.Z_ss_counts(1:cS.aR_new)) > 1e-9
    error('尽管存在正的工作人口，总劳动供给仍为零。');
elseif L <= 0, warning('总劳动供给为零或负。'); L = 1e-6; end

%% 4. 通过迭代求解一般均衡 (K*, TR_total*)
fprintf('\n--- 4. 求解一般均衡 (外生 PAYG 税率, PPS, Heer GBC) ---\n');

% 初始猜测值
KGuess = 49 ;     % 总生产性资本的猜测值
TRTotalGuess =0.4411; % 总净转移支付的猜测值 (政府净转移 + 意外遗赠)

% 迭代参数
maxIter = 200;
tolLevel = 1e-4;
dampK = 0.1;
dampTR = 0.1;
iter = 0; devNorm = inf;

fprintf('开始均衡迭代...\n');
fprintf('Iter |   K Guess  | TRTotGuess | RnetHHfac | WnetHH  | PAYG_b_iter | K Mod N-P | K Mod PPS | TRTot Mod |   K Dev    | TRTot Dev  |   Norm     | Time\n');
fprintf('----------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n');

eq_converged = false;
K_hist = zeros(maxIter+1,1); TR_hist = zeros(maxIter+1,1); Dev_hist = zeros(maxIter+1,1);
R_k_net_factor_hist = zeros(maxIter+1,1); B_payg_hist_iter_log = zeros(maxIter+1,1); % Renamed to avoid conflict
K_hist(1) = KGuess; TR_hist(1) = TRTotalGuess; Dev_hist(1) = devNorm;
R_k_net_factor_hist(1)=NaN; B_payg_hist_iter_log(1)=NaN;
KModel_nonpps = NaN; KModel_pps = NaN; TRTotalModel = NaN;

for iter = 1:maxIter
    iter_start_time = tic;

    % --- 第 4a 步: 计算税前价格, PAYG福利, 和家庭净价 ---
    [R_mkt_gross_factor_iter, MPL_gross_iter] = main_olg_v4_utils_pps.HHPrices_Huggett(KGuess, L, cS);
    r_mkt_gross_iter = R_mkt_gross_factor_iter - 1;

    theta_payg_exog_iter = cS.theta;
    total_gross_wage_bill_iter = MPL_gross_iter * L; % L is per capita, so this is per capita wage bill
    total_payg_revenue_iter = theta_payg_exog_iter * total_gross_wage_bill_iter; % Per capita revenue

    % paramS.Z_ss_counts are absolute counts. We need mass of retirees.
    mass_retirees_group_iter = sum(paramS.ageMassV(cS.aR_new + 1 : cS.aD_new)); % Share of retirees in pop
    b_payg_iter = 0;
    if mass_retirees_group_iter > 1e-9 && total_payg_revenue_iter >= 0
        b_payg_iter = total_payg_revenue_iter / mass_retirees_group_iter; % Benefit per retiree
    elseif mass_retirees_group_iter <= 1e-9 && total_payg_revenue_iter > 0
         warning('Iter %d: Positive PAYG revenue (%.2e) but no retirees! b_payg set to 0.', iter, total_payg_revenue_iter);
    end
    b_payg_iter = max(0, b_payg_iter);

    r_k_net_hh_iter = r_mkt_gross_iter * (1 - cS.tau_k);
    R_k_net_factor_hh_iter = 1 + r_k_net_hh_iter;
    w_net_hh_iter = MPL_gross_iter * (1 - theta_payg_exog_iter - cS.tau_l);
    w_net_hh_iter = max(0, w_net_hh_iter);

    bV_payg_for_vfi_iter = zeros(1, cS.aD_new);
    if cS.aR_new < cS.aD_new
        bV_payg_for_vfi_iter(cS.aR_new + 1 : cS.aD_new) = b_payg_iter;
    end

    % --- 第 4b 步: 求解家庭问题 (VFI) ---
    [cPolM, kPolM, cPpsPolM, ~] = main_olg_v4_utils_pps.HHSolution_VFI_Huggett(R_k_net_factor_hh_iter, w_net_hh_iter, TRTotalGuess, bV_payg_for_vfi_iter, paramS, cS);

    % --- 第 4c 步: 模拟家庭决策 ---
    [kHistM_non_pps, kPpsHistM, cHistM] = main_olg_v4_utils_pps.HHSimulation_olgm(kPolM, cPpsPolM, cPolM, eIdxM, R_k_net_factor_hh_iter, w_net_hh_iter, TRTotalGuess, bV_payg_for_vfi_iter, paramS, cS);

    % --- 第 4d 步: 计算模型的总生产性资本和总消费 ---
    % Aggregation uses Z_ss_norm_annual (sums to 1), results are per capita of total economy
    KModel_nonpps = mean(kHistM_non_pps, 1) * paramS.Z_ss_norm_annual;
    KModel_nonpps = max(1e-6, KModel_nonpps);
    KModel_pps = 0;
    if cS.pps_active && cS.pps_in_K && cS.pps_max_contrib_frac > 1e-9
        KModel_pps = mean(kPpsHistM, 1) * paramS.Z_ss_norm_annual;
        KModel_pps = max(0, KModel_pps);
    end
    KModel = KModel_nonpps + KModel_pps; % Per capita K

    CModel_iter = mean(cHistM,1) * paramS.Z_ss_norm_annual; % Per capita C

    % --- 第 4e 步: 计算模型的年度总意外遗赠 (T_bequest_Model) ---
    kprimeHistM_for_bequest = zeros(cS.nSim, cS.aD_orig);
    if cS.aD_orig > 1
        kprimeHistM_for_bequest(:, 1:cS.aD_orig-1) = kHistM_non_pps(:, 2:cS.aD_orig);
    end
    ageDeathMassV_annual = paramS.Z_ss_norm_annual(:) .* cS.d_orig(:); % Prob of dying at this age * mass at this age
    mean_bequest_value_iter = mean(kprimeHistM_for_bequest * R_k_net_factor_hh_iter, 1); % Avg bequest value if death occurs
    TotalBequests_iter_val_per_capita = sum(mean_bequest_value_iter(:) .* ageDeathMassV_annual(:)); % Total bequests per capita of economy
    T_bequest_Model = TotalBequests_iter_val_per_capita / (1 + paramS.popGrowthForDebt); % Distributed per capita next period
    T_bequest_Model = max(0, T_bequest_Model);

    % --- 第 4f 步: 政府预算平衡与模型隐含的总转移支付 (TRTotalModel) ---
    % All terms should be per capita of the economy
    Y_iter = cS.A * (KGuess^cS.alpha) * (L^(1-cS.alpha)); % KGuess and L are per capita
    G_fixed_val_iter = cS.gov_exp_frac_Y * Y_iter;
    B_fixed_val_iter = cS.gov_debt_frac_Y * Y_iter;

    LaborTaxRevenue_iter = (theta_payg_exog_iter + cS.tau_l) * MPL_gross_iter * L; % L is per capita labor
    CapitalTaxRevenue_iter = r_mkt_gross_iter * KGuess * cS.tau_k; % KGuess is per capita capital
    ConsumptionTaxRevenue_iter = CModel_iter * cS.tau_c; % CModel_iter is per capita consumption
    TotalTaxRevenue_iter = LaborTaxRevenue_iter + CapitalTaxRevenue_iter + ConsumptionTaxRevenue_iter; % All per capita

    r_b_iter = r_k_net_hh_iter;

    payg_outlay_iter = b_payg_iter * mass_retirees_group_iter; % This is (benefit per retiree) * (share of retirees) = per capita outlay

    TR_gov_Model = TotalTaxRevenue_iter - G_fixed_val_iter - payg_outlay_iter ...
                   - (r_b_iter - paramS.popGrowthForDebt) * B_fixed_val_iter;
    TR_gov_Model = max(0, TR_gov_Model); % Per capita gov transfers

    TRTotalModel = TR_gov_Model + T_bequest_Model; % Per capita total transfers

    % --- 第 4g 步: 计算偏差并检查收敛 ---
    KDev = KGuess - KModel;
    TRTotalDev = TRTotalGuess - TRTotalModel;
    devNorm = sqrt(KDev^2 + TRTotalDev^2);

    iter_time = toc(iter_start_time);
    fprintf('%4d | %10.4f | %10.4f | %9.4f | %9.4f | %11.4f | %10.4f | %10.4f | %10.4f | %10.4f | %10.4f | %10.4e | %.2fs\n', ...
             iter, KGuess, TRTotalGuess, R_k_net_factor_hh_iter, w_net_hh_iter, b_payg_iter, KModel_nonpps, KModel_pps, TRTotalModel, KDev, TRTotalDev, devNorm, iter_time);

     K_hist(iter+1) = KModel; TR_hist(iter+1) = TRTotalModel; Dev_hist(iter+1) = devNorm;
     R_k_net_factor_hist(iter+1) = R_k_net_factor_hh_iter; B_payg_hist_iter_log(iter+1) = b_payg_iter;

    if devNorm < tolLevel
        fprintf('均衡已收敛!\n');
        eq_converged = true;
        KGuess = KModel; TRTotalGuess = TRTotalModel;
        break;
    end

    KGuess = KGuess - dampK * KDev;
    TRTotalGuess = TRTotalGuess - dampTR * TRTotalDev;
    KGuess = max(1e-6, KGuess); TRTotalGuess = max(0, TRTotalGuess);
end

if ~eq_converged
    fprintf('\n警告: %d 次迭代后均衡未收敛。\n', maxIter);
    fprintf('最终偏差范数: %.4e\n', devNorm);
    if iter == maxIter
        KGuess = KModel; TRTotalGuess = TRTotalModel;
    end
end

K_eq = KGuess; TR_total_eq = TRTotalGuess;
R_k_net_factor_hh_eq = R_k_net_factor_hist(iter+1);

% --- 重新计算最终均衡价格与组成部分 ---
[R_mkt_gross_factor_eq, MPL_gross_eq] = main_olg_v4_utils_pps.HHPrices_Huggett(K_eq, L, cS);
r_mkt_gross_eq = R_mkt_gross_factor_eq - 1;

theta_payg_eq = cS.theta;
total_gross_wage_bill_eq_per_capita = MPL_gross_eq * L;
total_payg_revenue_eq_per_capita = theta_payg_eq * total_gross_wage_bill_eq_per_capita;
mass_retirees_group_eq = sum(paramS.ageMassV(cS.aR_new + 1 : cS.aD_new));
b_payg_eq = 0;
if mass_retirees_group_eq > 1e-9 && total_payg_revenue_eq_per_capita >= 0
    b_payg_eq = total_payg_revenue_eq_per_capita / mass_retirees_group_eq;
end
b_payg_eq = max(0, b_payg_eq);

w_net_hh_eq = MPL_gross_eq * (1 - theta_payg_eq - cS.tau_l); w_net_hh_eq = max(0,w_net_hh_eq);
bV_eq_new = zeros(1, cS.aD_new);
if cS.aR_new < cS.aD_new, bV_eq_new(cS.aR_new + 1 : cS.aD_new) = b_payg_eq; end

[cPolM_eq, kPolM_eq, cPpsPolM_eq, valueM_eq] = main_olg_v4_utils_pps.HHSolution_VFI_Huggett(R_k_net_factor_hh_eq, w_net_hh_eq, TR_total_eq, bV_eq_new, paramS, cS);
fprintf('模拟均衡状态下的最终分布...\n')
[kHistM_eq, kPpsHistM_eq, cHistM_eq] = main_olg_v4_utils_pps.HHSimulation_olgm(kPolM_eq, cPpsPolM_eq, cPolM_eq, eIdxM, R_k_net_factor_hh_eq, w_net_hh_eq, TR_total_eq, bV_eq_new, paramS, cS);
fprintf('最终模拟完成。\n')

C_eq = mean(cHistM_eq,1) * paramS.Z_ss_norm_annual;
Y_eq = cS.A * (K_eq^cS.alpha) * (L^(1-cS.alpha));
G_eq = cS.gov_exp_frac_Y * Y_eq;
B_eq = cS.gov_debt_frac_Y * Y_eq;

kprimeHistM_for_bequest_eq = zeros(cS.nSim, cS.aD_orig);
if cS.aD_orig > 1
    kprimeHistM_for_bequest_eq(:, 1:cS.aD_orig-1) = kHistM_eq(:, 2:cS.aD_orig);
end
ageDeathMassV_annual_eq = paramS.Z_ss_norm_annual(:) .* cS.d_orig(:);
mean_bequest_value_eq = mean(kprimeHistM_for_bequest_eq * R_k_net_factor_hh_eq, 1);
TotalBequests_eq_val_per_capita = sum(mean_bequest_value_eq(:) .* ageDeathMassV_annual_eq(:));
T_bequest_eq = TotalBequests_eq_val_per_capita / (1 + paramS.popGrowthForDebt);
T_bequest_eq = max(0, T_bequest_eq);

TR_gov_eq = TR_total_eq - T_bequest_eq;

LaborTax_eq = (theta_payg_eq + cS.tau_l) * MPL_gross_eq * L;
CapitalTax_eq = r_mkt_gross_eq * K_eq * cS.tau_k;
ConsumptionTax_eq = C_eq * cS.tau_c;
TotalTax_eq = LaborTax_eq + CapitalTax_eq + ConsumptionTax_eq;
PAYG_outlay_eq_per_capita = b_payg_eq * mass_retirees_group_eq;
r_b_eq = R_k_net_factor_hh_eq - 1;

GovBudgetResidual = TotalTax_eq - G_eq - TR_gov_eq - PAYG_outlay_eq_per_capita - (r_b_eq - paramS.popGrowthForDebt)*B_eq;
fprintf('GBC Check at Equilibrium: Residual = %.4e (should be close to zero)\n', GovBudgetResidual);
if abs(GovBudgetResidual) > 1e-3 * abs(Y_eq) && abs(Y_eq) > 1e-9
    warning('GBC residual is large at equilibrium! |Res| = %.2e, Y = %.2e', abs(GovBudgetResidual), Y_eq);
elseif abs(GovBudgetResidual) > 1e-3 && abs(Y_eq) <= 1e-9
    warning('GBC residual is large at equilibrium and Y is near zero! |Res| = %.2e', abs(GovBudgetResidual));
end


%% 5. 分析和绘制结果
fprintf('\n--- 5. 均衡结果与绘图 (外生 PAYG 税率 + PPS + Heer GBC) ---\n');
K_nonpps_eq = mean(kHistM_eq, 1) * paramS.Z_ss_norm_annual;
K_pps_eq = 0; if cS.pps_active && cS.pps_max_contrib_frac > 1e-9, K_pps_eq = mean(kPpsHistM_eq, 1) * paramS.Z_ss_norm_annual; end
TotalAssets_eq = K_nonpps_eq + K_pps_eq; % Per capita total assets

fprintf('均衡总生产性资本 (K*): %.4f (总体人均)\n', K_eq);
fprintf('  分解: 非PPS资本 K: %.4f, PPS资本 K: %.4f (均为总体人均)\n', K_nonpps_eq, K_pps_eq);
fprintf('均衡家庭总资产 (非PPS+PPS): %.4f (总体人均)\n', TotalAssets_eq);
fprintf('均衡总劳动 (L, 效率单位): %.4f (总体人均)\n', L);
fprintf('均衡总产出 (Y* 年度):   %.4f (总体人均)\n', Y_eq);
fprintf('均衡市场利率因子 (R_mkt_gross_factor*): %.4f (r_mkt_gross*=%.4f)\n', R_mkt_gross_factor_eq, r_mkt_gross_eq);
fprintf('  家庭净回报率因子 (R_k_net_factor_hh*): %.4f (r_k_net_hh*=%.4f)\n', R_k_net_factor_hh_eq, R_k_net_factor_hh_eq-1);
fprintf('均衡总工资率 (MPL_gross*): %.4f\n', MPL_gross_eq);
fprintf('均衡外生PAYG工资税率 (theta_payg): %.4f\n', theta_payg_eq);
fprintf('  其他劳动税率 (tau_l fixed): %.4f\n', cS.tau_l);
fprintf('  资本所得税率 (tau_k fixed): %.4f\n', cS.tau_k);
fprintf('  消费税率 (tau_c fixed): %.4f\n', cS.tau_c);
fprintf('均衡家庭净工资率 (w_net_hh*): %.4f\n', w_net_hh_eq);
fprintf('均衡内生PAYG福利 (b_payg*): %.4f (每位退休人员)\n', b_payg_eq);
fprintf('均衡总净转移支付 (TR_total*): %.4f (总体人均)\n', TR_total_eq);
fprintf('  其中意外遗赠 (T_bequest*): %.4f (总体人均)\n', T_bequest_eq);
fprintf('  其中政府净转移 (TR_gov*): %.4f (总体人均)\n', TR_gov_eq);
fprintf('均衡政府消费 (G*): %.4f (G/Y = %.2f) (总体人均)\n', G_eq, G_eq/Y_eq);
fprintf('均衡政府债务 (B*): %.4f (B/Y = %.2f) (总体人均)\n', B_eq, B_eq/Y_eq);
fprintf('均衡 K/Y 比率 (年度): %.4f\n', K_eq / Y_eq);
fprintf('均衡 C/Y 比率 (年度): %.4f\n', C_eq / Y_eq);

% --- 计算替代率 ---
mass_workers_group_eq = sum(paramS.ageMassV(1:cS.aR_new)); % Share of workers in pop
avg_gross_wage_per_worker_eq = 0;
if mass_workers_group_eq > 1e-9 && L > 0 && MPL_gross_eq > 0
    avg_gross_wage_per_worker_eq = (MPL_gross_eq * L) / mass_workers_group_eq;
end

payg_replacement_rate = 0;
if avg_gross_wage_per_worker_eq > 1e-9
    payg_replacement_rate = b_payg_eq / avg_gross_wage_per_worker_eq;
end
fprintf('PAYG 替代率 (b_payg / 平均工人总工资): %.4f (%.2f%%)\n', payg_replacement_rate, payg_replacement_rate*100);

avg_pps_benefit_for_payg_retiree = 0;
if cS.pps_active && cS.pps_max_contrib_frac > 1e-9 && mass_retirees_group_eq > 1e-9
    avg_kPps_by_annual_age = mean(kPpsHistM_eq, 1); % 1 x aD_orig
    pps_withdrawal_ages_indices = cS.pps_withdrawal_age_min : cS.aD_orig; % Annual indices

    if ~isempty(pps_withdrawal_ages_indices) && max(pps_withdrawal_ages_indices) <= cS.aD_orig && min(pps_withdrawal_ages_indices) >=1
        avg_kPps_at_withdrawal_annual_ages = avg_kPps_by_annual_age(pps_withdrawal_ages_indices);
        avg_pps_withdrawal_pretax_by_annual_age = avg_kPps_at_withdrawal_annual_ages * cS.pps_withdrawal_rate;
        avg_pps_withdrawal_posttax_by_annual_age = avg_pps_withdrawal_pretax_by_annual_age * (1 - cS.pps_tax_rate_withdrawal);

        Z_annual_mass_at_withdrawal_ages = paramS.Z_ss_norm_annual(pps_withdrawal_ages_indices);
        Total_PPS_PostTax_Withdrawals_PerCapitaEconomy = sum(avg_pps_withdrawal_posttax_by_annual_age(:) .* Z_annual_mass_at_withdrawal_ages(:));
        
        avg_pps_benefit_for_payg_retiree = Total_PPS_PostTax_Withdrawals_PerCapitaEconomy / mass_retirees_group_eq;
    else
        fprintf('警告: PPS 领取年龄指数超出范围，PPS平均福利未计算。\n');
    end
end
avg_pps_benefit_for_payg_retiree = max(0, avg_pps_benefit_for_payg_retiree);

total_replacement_rate = 0;
if avg_gross_wage_per_worker_eq > 1e-9
    total_replacement_rate = (b_payg_eq + avg_pps_benefit_for_payg_retiree) / avg_gross_wage_per_worker_eq;
end
fprintf('  其中平均PPS福利 (按PAYG退休人员计): %.4f\n', avg_pps_benefit_for_payg_retiree);
fprintf('PAYG+PPS 替代率 ((b_payg + pps_benefit) / 平均工人总工资): %.4f (%.2f%%)\n', total_replacement_rate, total_replacement_rate*100);


% --- 绘图 ---
figure('Name', '均衡收敛过程 (Heer GBC, Exo Theta)');
subplot(4,1,1); plot(0:iter, K_hist(1:iter+1), '-o'); title('K 收敛'); ylabel('K (总体人均)'); grid on; xlim([0 iter]);
subplot(4,1,2); plot(0:iter, TR_hist(1:iter+1), '-o'); title('TR_{total} 收敛'); ylabel('TR_{total} (总体人均)'); grid on; xlim([0 iter]);
subplot(4,1,3); plot(0:iter, B_payg_hist_iter_log(1:iter+1), '-o'); title('b_{PAYG} 收敛'); ylabel('b_{PAYG} (每位退休人员)'); grid on; xlim([0 iter]);
subplot(4,1,4); semilogy(1:iter, Dev_hist(2:iter+1), '-o'); title('偏差范数收敛'); ylabel('Norm'); xlabel('迭代次数'); grid on; xlim([1 iter]);
sgtitle('均衡迭代收敛过程 (Heer GBC, Exo Theta)');

plot_a_new = min(4, cS.aR_new);
if plot_a_new <= cS.aD_new && plot_a_new > 0 && ~isempty(cS.physAgeMap{plot_a_new}) && cS.nk > 1
    physAgeStart = cS.physAgeV_new(plot_a_new);
    physAgeEnd = cS.physAgeV_orig(cS.physAgeMap{plot_a_new}(end));
    plot_title_suffix = sprintf('年龄组 %d (年龄 %d-%d)', plot_a_new, physAgeStart, physAgeEnd);
    figure('Name', ['政策函数 (Heer GBC, Exo Theta): ' plot_title_suffix]);
    subplot(1, 3, 1);
    plot(cS.kGridV, cPolM_eq(:, 1, plot_a_new), 'b-', 'LineWidth', 1, 'DisplayName', sprintf('e=%d (低)',1)); hold on;
    plot(cS.kGridV, cPolM_eq(:, round(cS.nw/2), plot_a_new), 'k--', 'LineWidth', 1, 'DisplayName', sprintf('e=%d (中)',round(cS.nw/2)));
    plot(cS.kGridV, cPolM_eq(:, cS.nw, plot_a_new), 'r:', 'LineWidth', 1, 'DisplayName', sprintf('e=%d (高)',cS.nw)); hold off;
    xlabel('期初非PPS资本 (k)'); ylabel('最优消费 (c_{quantity})'); title(['消费政策 ' plot_title_suffix]); legend('Location', 'best'); grid on;
    subplot(1, 3, 2);
    plot(cS.kGridV, kPolM_eq(:, 1, plot_a_new), 'b-', 'LineWidth', 1, 'DisplayName', sprintf('e=%d (低)',1)); hold on;
    plot(cS.kGridV, kPolM_eq(:, round(cS.nw/2), plot_a_new), 'k--', 'LineWidth', 1, 'DisplayName', sprintf('e=%d (中)',round(cS.nw/2)));
    plot(cS.kGridV, kPolM_eq(:, cS.nw, plot_a_new), 'r:', 'LineWidth', 1, 'DisplayName', sprintf('e=%d (高)',cS.nw));
    plot(cS.kGridV, cS.kGridV, 'g--', 'DisplayName', 'k'' = k'); hold off;
    xlabel('期初非PPS资本 (k)'); ylabel('下一期非PPS资本 (k'')'); title(['非PPS储蓄政策 ' plot_title_suffix]); legend('Location', 'best'); grid on;
    subplot(1, 3, 3);
    if cS.pps_active && cS.pps_max_contrib_frac > 1e-9
        plot(cS.kGridV, cPpsPolM_eq(:, 1, plot_a_new), 'b-', 'LineWidth', 1, 'DisplayName', sprintf('e=%d (低)',1)); hold on;
        plot(cS.kGridV, cPpsPolM_eq(:, round(cS.nw/2), plot_a_new), 'k--', 'LineWidth', 1, 'DisplayName', sprintf('e=%d (中)',round(cS.nw/2)));
        plot(cS.kGridV, cPpsPolM_eq(:, cS.nw, plot_a_new), 'r:', 'LineWidth', 1, 'DisplayName', sprintf('e=%d (高)',cS.nw)); hold off;
        ylabel('PPS缴费 (c_{pps})'); title(['PPS缴费政策 ' plot_title_suffix]); legend('Location', 'best');
    else, title('PPS 未缴费或未激活'); axis off; end
     xlabel('期初非PPS资本 (k)'); grid on;
else, fprintf('跳过政策函数绘图 - plot_a_new 无效, physAgeMap 为空, 或 nk<=1。\n'); end

avgK_byAge_annual = mean(kHistM_eq, 1);
avgKpps_byAge_annual = zeros(1, cS.aD_orig);
if cS.pps_active && cS.pps_max_contrib_frac > 1e-9, avgKpps_byAge_annual = mean(kPpsHistM_eq, 1); end
avgTotalAssets_byAge_annual = avgK_byAge_annual + avgKpps_byAge_annual;

figure('Name', '生命周期资本剖面 (Heer GBC, Exo Theta)');
plot(cS.physAgeV_orig, avgK_byAge_annual, 'bo-', 'LineWidth', 1.5, 'DisplayName', '平均非PPS资本 (k)'); hold on;
if cS.pps_active && cS.pps_max_contrib_frac > 1e-9
    plot(cS.physAgeV_orig, avgKpps_byAge_annual, 'ro-', 'LineWidth', 1.5, 'DisplayName', '平均PPS资本 (k_{pps})');
    plot(cS.physAgeV_orig, avgTotalAssets_byAge_annual, 'k--', 'LineWidth', 1.5, 'DisplayName', '平均总资产 (k+k_{pps})');
end
xline(cS.ageRetire_orig, 'g--', 'LineWidth', 1, 'DisplayName', '退休开始'); hold off;
title('按年度年龄划分的平均期初资本持有量'); xlabel('真实年龄'); ylabel('平均资本'); legend('Location', 'best'); grid on; xlim([cS.age1_orig, cS.ageLast_orig]);

fprintf('\n--- 分析完成 (外生 PAYG 税率 + PPS + Heer GBC) ---\n');

% --- END OF FILE main_olg_v4.m ---