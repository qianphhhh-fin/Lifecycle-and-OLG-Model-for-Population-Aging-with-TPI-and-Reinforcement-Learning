% --- START OF FILE debug_noDeath.m (Corrected v5 - Final Audit Passed) ---

% =========================================================================
% == SCRIPT: Aiyagari/Huggett Model Single-Period Accounting Audit
% == (Simplified from OLG model with no death)
% =========================================================================
clear; close all;
addpath(pwd);
fprintf('=== Aiyagari/Huggett 模型单期会计审计 (无死亡风险) ===\n\n');

%% 1. 初始化环境和参数
fprintf('--- 1. 初始化 Aiyagari/Huggett 环境...\n');
cS = main_olg_v14_utils.ParameterValues_HuggettStyle();
paramS = struct();

% --- 基本参数 ---
cS.sigma = 3.0; cS.beta = 0.97; cS.cFloor = 0.05; cS.nSim = 5000;
cS.alpha = 0.5; cS.ddk = 0.08; 
cS.A = 1.0;

% --- 税收政策 ---
cS.tau_k = 0.25; cS.tau_l = 0.08; cS.tau_c = 0.04;
cS.gov_exp_frac_Y = 0.15;

% --- 简化设置 ---
cS.aD_new = 2; 
cS.s_1yr_transitionV = 1; cS.phi_bequest = 0; cS.theta_path = 0; cS.pps_active = false;

% --- 网格和冲击 ---
ngrid = 10; cS.nk = ngrid; cS.nkprime = ngrid; cS.nkpps = 1;
cS = main_olg_v14_utils.generateGrids(cS);
[paramS.leLogGridV, paramS.leTrProbM, paramS.leProb1V] = main_olg_v14_utils.EarningProcess_olgm(cS);
paramS.leGridV = exp(paramS.leLogGridV(:));
eIdxM = main_olg_v14_utils.MarkovChainSimulation_AgeGroup(cS.nSim, cS, paramS.leProb1V, paramS.leTrProbM);

fprintf('✅ Aiyagari/Huggett 环境初始化完成。\n\n');

%% 2. [核心审计] 对单个时期进行会计闭环检验
fprintf('--- 2. 开始对模型进行单期隔离审计 ---\n');
fprintf('--------------------------------------------------------\n\n');

% --- 步骤 A: 设定 t 期的初始状态 ---
K_t = 1.5; % 假设期初总资本

% --- 使用模拟样本的平均效率作为总劳动供给 ---
e_idx_now = eIdxM(:, 1);
epsilon_sim_dist = paramS.leGridV(e_idx_now);
L_t = mean(epsilon_sim_dist);
fprintf('   纠正后的总劳动供给 L_t (基于模拟样本均值): %.4f\n', L_t);

% --- 步骤 B: 获取当期价格 ---
[r_mkt_t, w_t, Y_t] = main_olg_v14_utils.HHPrices_Huggett(K_t, L_t, cS.A, cS);
r_net_t = r_mkt_t * (1 - cS.tau_k);
M_t = struct('w_t', w_t, 'r_mkt_period', r_mkt_t, 'r_net_period', r_net_t);

% --- 步骤 C: 【模拟前】计算家庭部门的理论总资源 (Total Sources) ---
fprintf('--- 审计点 1: 家庭部门总资源核算 ---\n');
TotalWealth_t = K_t;
TotalReturn_k = K_t * r_net_t;
GrossLaborIncome = w_t * L_t; 
Theoretical_TotalSources_HH = TotalWealth_t + TotalReturn_k + GrossLaborIncome;
fprintf('%-35s | %15.4f\n', '理论家庭总资源', Theoretical_TotalSources_HH);
fprintf('\n');

% --- 步骤 D: 【运行模拟】 ---
fprintf('--- 正在调用 main_olg_v14_utils_noDeath.solve_and_simulate_onestep... ---\n');
k_initial_dist = ones(cS.nSim, 1) * K_t; 
[K_next, C_t] = main_olg_v14_utils_noDeath.solve_and_simulate_onestep(M_t, cS, paramS, eIdxM, k_initial_dist);
fprintf('--- 调用完成 ---\n\n');

% --- 步骤 E: 【模拟后】计算家庭部门的实际总运用 (Total Uses) ---
fprintf('--- 审计点 2: 家庭部门总运用核算 ---\n');
TotalWealth_next = K_next;
ConsumptionExpenditure_agg = C_t * (1 + cS.tau_c);
LaborTax_agg = cS.tau_l * GrossLaborIncome; 
Actual_TotalUses_HH = TotalWealth_next + ConsumptionExpenditure_agg + LaborTax_agg;
fprintf('%-35s | %15.4f\n', '实际家庭总运用', Actual_TotalUses_HH);
fprintf('\n');

% --- 家庭部门预算检验 ---
HouseholdBudgetResidual = Theoretical_TotalSources_HH - Actual_TotalUses_HH;
fprintf('%-35s | %15.4e | %s\n', 'A. 家庭部门预算检验', HouseholdBudgetResidual, '理论总资源 - 实际总运用');
fprintf('\n');

% --- 步骤 F: 继续进行完整的宏观审计 ---
fprintf('--- 审计点 3: 完整的宏观会计恒等式 ---\n');

% --- 供给-收入检验 ---
GrossCapitalReturn = (r_mkt_t + cS.ddk) * K_t;
Supply_Income_Residual = Y_t - (GrossLaborIncome + GrossCapitalReturn);
fprintf('%-35s | %15.4e | %s\n', 'B. 供给-收入检验', Supply_Income_Residual, 'Y - wL - MPK*K');

% =========================================================================
% == [最终审计修正] 正确核算包含政府储蓄的总投资 I
% =========================================================================
% 1. 计算政府总税收收入 (Total Tax Revenue)
CapitalTax_agg = cS.tau_k * r_mkt_t * K_t;
ConsumptionTax_agg = cS.tau_c * C_t;
TotalTaxRevenue = LaborTax_agg + CapitalTax_agg + ConsumptionTax_agg;

% 2. 计算政府支出和政府储蓄 (S_gov)
G_t = cS.gov_exp_frac_Y * Y_t;
S_gov = TotalTaxRevenue - G_t;

% 3. 计算私人净储蓄 (S_pvt)
S_pvt = K_next - K_t;

% 4. 计算总的净投资 (Net Investment = Total Savings)
I_net = S_pvt + S_gov;

% 5. 计算总的粗投资 (Gross Investment)
I_gross = I_net + cS.ddk * K_t;
% =========================================================================

% --- 商品市场出清检验 ---
Market_Clearing_Residual = Y_t - (C_t + I_gross + G_t);
fprintf('%-35s | %15.4e | %s\n', 'C. 市场出清残差', Market_Clearing_Residual, 'Y - C - I - G');
fprintf('\n');

% --- 最终结论 ---
if abs(Market_Clearing_Residual) < 1e-9 && abs(Supply_Income_Residual) < 1e-9 && abs(HouseholdBudgetResidual) < 1e-9
    fprintf('\n✅ 审计通过！Aiyagari/Huggett模型所有核心会计恒等式精确闭合！\n');
else
    fprintf('\n❌ 审计失败！请检查仍然存在的差异。\n');
    fprintf('   家庭预算残差: %.4e\n', HouseholdBudgetResidual);
    fprintf('   市场出清残差: %.4e\n', Market_Clearing_Residual);
end