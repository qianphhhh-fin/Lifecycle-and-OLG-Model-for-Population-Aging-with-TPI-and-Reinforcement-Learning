% --- START OF FILE debug_hh_budget.m ---
% =========================================================================
% == SCRIPT: Household Budget Discrepancy Debugger
% == 目的：专门用于诊断和解决家庭部门理论资源与实际运用不匹配的问题。
% =========================================================================
clear; close all;
addpath(pwd);
fprintf('=== 家庭预算闭环专项调试脚本 ===\n\n');

%% 1. 初始化一个完全相同的环境
fprintf('--- 1. 初始化OLG环境 (零税率调试模式) ---\n');
cS = main_olg_v14_utils.ParameterValues_HuggettStyle();
paramS = struct();
cS.sigma = 3.0; cS.beta = 0.97; cS.cFloor = 0.05; cS.nSim = 1000;
cS.phi_bequest = 3; cS.sigma_bequest = cS.sigma;
cS.start_year = 1997; cS.end_year = 2102; cS.time_step = 5;
cS.alpha = 0.5; cS.ddk = 1 - (1 - 0.05)^cS.time_step;
% --- 零税率设定 ---
cS.tau_k = 0.0; cS.tau_l = 0.0; cS.tau_c = 0.0;
cS.pps_tax_rate_withdrawal = 0.0;
% --- 其他参数 ---
cS.gov_exp_frac_Y = 0.15; cS.rho_prime_payg = 0.5; cS.pps_activation_year = 2022;
cS.pps_return_rate_premium = 0.01; cS.pps_withdrawal_rate = 0.15;
cS.pps_contrib_limit = 9999; cS.pps_max_contrib_frac = 0.1;
ngrid = 20; cS.nk = ngrid; cS.nkpps = ngrid; cS.nkprime = ngrid; cS.npps = ngrid;
cS.A = 1.0; cS = main_olg_v14_utils.generateGrids(cS);
cS = main_olg_v14_utils.calcaulte_theta_payg_path(cS, false);
[Z_path, A_path, T] = main_olg_v14_utils.load_exogenous_paths(cS);
sim_years = cS.start_year:cS.time_Step:(cS.start_year + (T-1)*cS.time_Step);
cS.sim_years = sim_years;
[paramS.leLogGridV, paramS.leTrProbM, paramS.leProb1V] = main_olg_v14_utils.EarningProcess_olgm(cS);
paramS.leGridV = exp(paramS.leLogGridV(:));
paramS.ageEffV_new = cS.ageEffV_new;
eIdxM = main_olg_v14_utils.LaborEndowSimulation_olgm_AgeGroup(cS, paramS);
fprintf('✅ OLG环境初始化完成。\n\n');

%% 2. [宏观层面审计] 重现问题
t = 2;
fprintf('--- 2. 重现宏观层面的预算不匹配问题 (t=%d) ---\n', t);
K_pvt_t = 1.5; K_pps_t = 0.1; B_p_t = 0.2; B_g_t = 0.3;
TR_total_t = 0.01; Z_t = Z_path(:, t); A_t = A_path(t);
cS.pps_active = (sim_years(t) >= cS.pps_activation_year);
M_t = main_olg_v14_utils.get_prices_and_policy_at_t(t, K_pvt_t, B_p_t, K_pps_t, B_g_t, Z_t, A_t, cS, paramS, eIdxM);

% --- 理论总资源 (宏观加总) ---
TotalWealth_t = K_pvt_t + K_pps_t;
r_mkt = M_t.r_mkt_period;
r_net = M_t.r_net_period;
TotalReturn_k_pvt = K_pvt_t * r_net;
TotalReturn_k_pps = K_pps_t * (r_mkt + cS.pps_return_rate_premium);
GrossLaborIncome = M_t.w_t * M_t.L_t;
PensionBenefit_agg = M_t.b_t * sum(Z_t(cS.aR_new+1:end));
BequestReceived_agg = TR_total_t * sum(Z_t .* cS.s_1yr_transitionV);
Theoretical_TotalSources_HH = TotalWealth_t + TotalReturn_k_pvt + TotalReturn_k_pps + GrossLaborIncome + PensionBenefit_agg + BequestReceived_agg;

% --- 运行模拟 ---
[K_pvt_next, K_pps_next, C_t, Total_Cpps_t, Total_PpsTax_t, total_accidental_bequest] = ...
    main_olg_v14_utils.simulate_private_capital_forward(M_t, Z_t, cS, paramS, eIdxM, TR_total_t);

% --- 实际总运用 (宏观加总) ---
TotalWealth_next = K_pvt_next + K_pps_next;
ConsumptionExpenditure_agg = C_t * (1 + cS.tau_c);
PaygTax_agg = cS.theta_path(t) * GrossLaborIncome;
LaborTax_agg = cS.tau_l * max(0, GrossLaborIncome - Total_Cpps_t);
PpsContrib_agg = Total_Cpps_t;
PpsWithdrawalTax_agg = Total_PpsTax_t;
BequestLeft_agg = total_accidental_bequest;
Actual_TotalUses_HH = TotalWealth_next + ConsumptionExpenditure_agg + PaygTax_agg + LaborTax_agg + PpsContrib_agg + PpsWithdrawalTax_agg + BequestLeft_agg;

% --- 打印宏观结果 ---
HouseholdBudgetResidual = Theoretical_TotalSources_HH - Actual_TotalUses_HH;
fprintf('%-35s | %15.4f\n', '理论家庭总资源', Theoretical_TotalSources_HH);
fprintf('%-35s | %15.4f\n', '实际家庭总运用', Actual_TotalUses_HH);
fprintf('%-35s | %15.4e\n', '宏观预算残差', HouseholdBudgetResidual);
fprintf('--------------------------------------------------------\n\n');

%% 3. [微观层面诊断] 进入VFI内部，检查单个家庭的资源计算
fprintf('--- 3. 微观诊断：检查VFI内部单个状态点的资源计算 ---\n');

% --- a. 选择一个有代表性的状态点进行检查 ---
a_idx_check = 5;      % 一个典型的中年工作年龄组
ik_check = 10;        % 中等水平的非PPS资产
ikpps_check = 5;      % 中等水平的PPS资产
ie_check = 3;         % 平均的劳动效率
icpps_check_idx = 2;  % 一个非零的PPS缴费选择

% --- b. 获取该状态点对应的变量 ---
k_state = cS.kGridV(ik_check);
k_pps_state = cS.kppsGridV(ikpps_check);
epsilon_state = paramS.leGridV(ie_check);
b_age_val = 0; % 因为是工作年龄
if a_idx_check > cS.aR_new, b_age_val = M_t.b_t; end

% --- c. 手动计算该【单个】家庭的理论资源 ---
% 注意：这里的TR_total_t是给所有幸存者的总额，单个幸存者收到的是 TR_total_t / sum(Z_t .* s)
% 为简化，我们直接使用总额，因为VFI中也是这样处理的。
[pretax_non_capital_income_ind, c_pps_choice_ind] = main_olg_v14_utils.HHIncome_Huggett(M_t.w_t, TR_total_t, b_age_val, 0, a_idx_check, paramS, cS, epsilon_state);
% 假设一个缴费选择
gross_labor_income_ind = M_t.w_t * cS.ageEffV_new(a_idx_check) * epsilon_state;
cpps_grid = linspace(0, max(0, gross_labor_income_ind * cS.pps_max_contrib_frac), cS.npps);
c_pps_choice_ind = cpps_grid(icpps_check_idx);
[pretax_non_capital_income_ind, ~] = main_olg_v14_utils.HHIncome_Huggett(M_t.w_t, TR_total_t, b_age_val, c_pps_choice_ind, a_idx_check, paramS, cS, epsilon_state);

k_end_value_ind = k_state * (1 + r_net);
pps_withdrawal_gross_ind = 0; if a_idx_check >= cS.aR_new && cS.pps_active, pps_withdrawal_gross_ind = k_pps_state * cS.pps_withdrawal_rate; end
pps_withdrawal_net_ind = pps_withdrawal_gross_ind * (1 - cS.pps_tax_rate_withdrawal);
% 理论税前总资源
Theoretical_TotalResources_pretax_ind = k_end_value_ind + pretax_non_capital_income_ind + pps_withdrawal_net_ind;
% 理论税负
payg_tax_ind = cS.theta_path(t) * gross_labor_income_ind;
labor_tax_ind = cS.tau_l * max(0, gross_labor_income_ind - c_pps_choice_ind);
% 理论可支配资源
Theoretical_DisposableResources_ind = Theoretical_TotalResources_pretax_ind - payg_tax_ind - labor_tax_ind - c_pps_choice_ind;

fprintf('--- 检查点: a=%d, k=%.2f, k_pps=%.2f, e=%.2f, c_pps=%.2f ---\n', ...
    a_idx_check, k_state, k_pps_state, epsilon_state, c_pps_choice_ind);
fprintf('%-35s | %15.4f\n', '理论可支配资源 (手动计算)', Theoretical_DisposableResources_ind);

% --- d. 在VFI函数内部相同位置，计算模型感知的资源 ---
% (这部分代码直接从 HHSolutionByAge_VFI_GridSearch 复制而来)
[pretax_non_capital_income_vfi, ~] = main_olg_v14_utils.HHIncome_Huggett(M_t.w_t, TR_total_t, b_age_val, c_pps_choice_ind, a_idx_check, paramS, cS, epsilon_state);
k_end_value_vfi = k_state * (1 + M_t.r_net_period);
pps_withdrawal_gross_vfi = 0; if a_idx_check >= cS.aR_new && cS.pps_active, pps_withdrawal_gross_vfi = k_pps_state * cS.pps_withdrawal_rate; end
pps_withdrawal_net_vfi = pps_withdrawal_gross_vfi * (1 - cS.pps_tax_rate_withdrawal);
total_resources_pretax_vfi = k_end_value_vfi + pretax_non_capital_income_vfi + pps_withdrawal_net_vfi;
labor_income_gross_vfi = 0; if a_idx_check <= cS.aR_new, labor_income_gross_vfi = M_t.w_t * cS.ageEffV_new(a_idx_check) * epsilon_state; end
payg_tax_vfi = cS.theta_path(t) * labor_income_gross_vfi;
labor_tax_vfi = cS.tau_l * max(0, labor_income_gross_vfi - c_pps_choice_ind);
VFI_DisposableResources = total_resources_pretax_vfi - payg_tax_vfi - labor_tax_vfi - c_pps_choice_ind;

fprintf('%-35s | %15.4f\n', 'VFI内部感知的可支配资源', VFI_DisposableResources);

% --- e. 比较差异 ---
Micro_Residual = Theoretical_DisposableResources_ind - VFI_DisposableResources;
fprintf('%-35s | %15.4e\n', '微观资源计算残差', Micro_Residual);

if abs(Micro_Residual) < 1e-9
    fprintf('\n✅ 微观诊断通过！VFI内部资源计算与理论一致。\n');
    fprintf('   问题可能出在从VFI策略到宏观聚合的 simulate_private_capital_forward 或 HHSimulation_olgm 环节。\n');
else
    fprintf('\n❌ 微观诊断失败！VFI内部资源计算与理论【不一致】！\n');
    fprintf('   请仔细比对 debug_hh_budget.m 中手动计算的每一步与 HHSolutionByAge_VFI_GridSearch 中的对应行。\n');
end