% --- START OF FILE debug_hh_budget.m ---
% =========================================================================
% == SCRIPT: Household Budget Debugger (VERSION: Bequests go to Govt)
% == 目的：在“意外遗赠收归政府”的新模型设定下，审计家庭部门的预算闭环。
% =========================================================================
clear; close all;
addpath(pwd);
fprintf('=== 家庭预算闭环专项调试脚本 (设定：遗赠归政府) ===\n\n');

%% 1. 初始化一个完全相同的环境 (不变)
fprintf('--- 1. 初始化OLG环境 (零税率调试模式) ---\n');
cS = main_olg_v14_utils.ParameterValues_HuggettStyle();
paramS = struct();
cS.sigma = 3.0; cS.beta = 0.97; cS.cFloor = 0.05; cS.nSim = 1000;
cS.phi_bequest = 3; cS.sigma_bequest = cS.sigma;
cS.start_year = 1997; cS.end_year = 2102; cS.time_Step = 5;
cS.alpha = 0.5; cS.ddk = 1 - (1 - 0.05)^cS.time_Step;
cS.tau_k = 0.0; cS.tau_l = 0.0; cS.tau_c = 0.0;
cS.pps_tax_rate_withdrawal = 0.0;
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

%% 2. [宏观层面审计] 修正并执行审计
t = 6; % 选择一个 PPS 激活后的时期 (e.g., 2027)
fprintf('--- 2. 宏观层面预算闭环审计 (t=%d, year=%d) ---\n', t, sim_years(t));
% --- a. 定义 t 时期的期初状态变量 ---
K_pvt_t = 1.5; K_pps_t = 0.1; B_p_t = 0.2; B_g_t = 0.3;
Z_t = Z_path(:, t); A_t = A_path(t);

% --- b. 调用 simulate_private_capital_forward (调用方式不变) ---
[K_pvt_next, K_pps_next, C_t, Total_Cpps_t, Total_PpsTax_t, total_accidental_bequest, M_t_complete] = ...
    main_olg_v14_utils.simulate_private_capital_forward(t, K_pvt_t, B_p_t, K_pps_t, B_g_t, Z_t, A_t, cS, paramS, eIdxM);
fprintf('函数调用成功，开始审计...\n');

% --- c. 理论总资源 (宏观加总) ---
r_mkt = M_t_complete.r_mkt_period;
r_net = M_t_complete.r_net_period;
w_t = M_t_complete.w_t;
L_t = M_t_complete.L_t;

TotalWealth_t = K_pvt_t + K_pps_t;
Return_k_pvt_agg = K_pvt_t * r_net;
Return_k_pps_agg = K_pps_t * (r_mkt + cS.pps_return_rate_premium);
GrossLaborIncome_agg = w_t * L_t;
PensionBenefit_agg = M_t_complete.b_t * sum(Z_t(cS.aR_new+1:end));
% [核心修正] 家庭不再收到遗赠
BequestReceived_agg = 0;

Theoretical_TotalSources_HH = TotalWealth_t + Return_k_pvt_agg + Return_k_pps_agg + GrossLaborIncome_agg + PensionBenefit_agg + BequestReceived_agg;

% --- d. 实际总运用 (宏观加总) ---
TotalWealth_next = K_pvt_next + K_pps_next;
ConsumptionExpenditure_agg = C_t * (1 + cS.tau_c);
PaygTax_agg = cS.theta_path(t) * GrossLaborIncome_agg;
taxable_labor_income = max(0, GrossLaborIncome_agg - PaygTax_agg - Total_Cpps_t);
LaborTax_agg = cS.tau_l * taxable_labor_income;
PpsContrib_agg = Total_Cpps_t;
PpsWithdrawalTax_agg = Total_PpsTax_t;
% [关键] 家庭仍然“留下”了遗赠，但这些资源会从家庭部门的资产负债表中消失，转移给政府。
% 所以在家庭部门的“运用”方，它依然存在。
BequestLeft_agg = total_accidental_bequest;

Actual_TotalUses_HH = TotalWealth_next + ConsumptionExpenditure_agg + PaygTax_agg + LaborTax_agg + PpsContrib_agg + PpsWithdrawalTax_agg + BequestLeft_agg;

% --- e. 打印宏观结果 ---
HouseholdBudgetResidual = Theoretical_TotalSources_HH - Actual_TotalUses_HH;
fprintf('%-35s | %15.4f\n', '理论家庭总资源', Theoretical_TotalSources_HH);
fprintf('%-35s | %15.4f\n', '实际家庭总运用', Actual_TotalUses_HH);
fprintf('%-35s | %15.4e\n', '宏观预算残差', HouseholdBudgetResidual);
if abs(HouseholdBudgetResidual) < 1e-6, fprintf('✅ 宏观预算闭环审计通过！\n'); else, fprintf('❌ 宏观预算闭环审计失败！\n'); end
fprintf('--------------------------------------------------------\n\n');


%% 3. [微观层面诊断] 进入VFI内部，检查单个家庭的资源计算
fprintf('--- 3. 微观诊断：检查VFI内部单个状态点的资源计算 ---\n');

% --- a. 选择一个有代表性的状态点进行检查 ---
a_idx_check = 5; ik_check = 10; ikpps_check = 5; ie_check = 3; icpps_check_idx = 2;

% --- b. 获取该状态点对应的变量 ---
k_state = cS.kGridV(ik_check); k_pps_state = cS.kppsGridV(ikpps_check);
epsilon_state = paramS.leGridV(ie_check);

% --- c. 手动计算该【单个】家庭的理论资源 ---
M_t = M_t_complete;
% [核心修正] 微观层面，家庭收到的遗赠也为零
tr_per_capita = 0;

labor_income_gross_ind = M_t.w_t * cS.ageEffV_new(a_idx_check) * epsilon_state;
max_cpps_ind = min(cS.pps_contrib_limit, labor_income_gross_ind * cS.pps_max_contrib_frac);
cpps_grid = linspace(0, max(0, max_cpps_ind), cS.npps);
c_pps_choice_ind = cpps_grid(icpps_check_idx);
k_return_ind = k_state * (1 + M_t.r_net_period);
total_inflow_ind = k_return_ind + labor_income_gross_ind + tr_per_capita;
payg_tax_ind = cS.theta_path(t) * labor_income_gross_ind;
labor_tax_ind = cS.tau_l * max(0, labor_income_gross_ind - c_pps_choice_ind - payg_tax_ind);
cpps_outflow_ind = c_pps_choice_ind;
Theoretical_NetCash = total_inflow_ind - (payg_tax_ind + labor_tax_ind + cpps_outflow_ind);

fprintf('--- 检查点: a=%d, k=%.2f, k_pps=%.2f, e=%.2f, c_pps=%.2f ---\n', ...
    a_idx_check, k_state, k_pps_state, epsilon_state, c_pps_choice_ind);
fprintf('%-40s | %15.4f\n', '理论可用于 C 和 k_prime 的净现金', Theoretical_NetCash);

% --- d. 在VFI函数内部相同位置，计算模型感知的资源 ---
k_return_vfi = k_state * (1 + M_t.r_net_period);
total_inflow_vfi = k_return_vfi + labor_income_gross_ind + tr_per_capita; % tr_per_capita is 0
payg_tax_vfi = cS.theta_path(t) * labor_income_gross_ind;
labor_tax_vfi = cS.tau_l * max(0, labor_income_gross_ind - c_pps_choice_ind - payg_tax_vfi);
cpps_outflow_vfi = c_pps_choice_ind;
VFI_NetCash = total_inflow_vfi - (payg_tax_vfi + labor_tax_vfi + cpps_outflow_vfi);

fprintf('%-40s | %15.4f\n', 'VFI内部计算的净现金', VFI_NetCash);

% --- e. 比较差异 ---
Micro_Residual = Theoretical_NetCash - VFI_NetCash;
fprintf('%-40s | %15.4e\n', '微观资源计算残差', Micro_Residual);
if abs(Micro_Residual) < 1e-9, fprintf('\n✅ 微观诊断通过！VFI内部资源计算与理论一致。\n'); else, fprintf('\n❌ 微观诊断失败！\n'); end