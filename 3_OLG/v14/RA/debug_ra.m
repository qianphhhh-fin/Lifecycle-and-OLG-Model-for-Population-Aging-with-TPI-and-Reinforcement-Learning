% --- START OF FILE debug_ra.m ---

% =========================================================================
% == SCRIPT: debug_ra.m (RA Model Audit - Final Corrected Version v3)
% == 目的：使用一个永续存在的代表性代理(RA)模型来精确检验宏观会计恒等式。
% =========================================================================
clear; close all;
fprintf('=== RA 模型调试脚本 (Accounting Audit) ===\n\n');

%% 1. 加载环境 (适用于RA模型)
fprintf('--- 1. 正在为RA模型重建环境...\n');
addpath(pwd);
cS_base = main_olg_v14_utils.ParameterValues_HuggettStyle();
cS_base.sigma = 3.0; cS_base.beta = 0.97; cS_base.cFloor = 0.05;
cS_base.start_year = 1997; cS_base.end_year = 2102; cS_base.time_step = 5;
cS_base.alpha = 0.5; cS_base.ddk = 1 - (1 - 0.05)^cS_base.time_step;
cS_base.tau_k = 0.25; cS_base.tau_l = 0.08; cS_base.tau_c = 0.04; cS_base.gov_exp_frac_Y = 0.15;
cS_base.A = 1.0;
cS = main_olg_v14_utils_ra.ParameterValues_RA(cS_base);
[~, A_path, T] = main_olg_v14_utils_ra.load_exogenous_paths_RA(cS);
sim_years = cS.start_year:cS.time_step:(cS.start_year + (T-1)*cS.time_step);
fprintf('✅ RA环境重建完成。\n\n');

%% 2. 设置 t=1 的初始状态并演化到 t=2
fprintf('--- 2. 正在设置 t=1 状态并演化到 t=2 ...\n');
paths = struct();
paths.K_pvt = zeros(1, T); paths.B_g = zeros(1, T); paths.K_physical = zeros(1, T);
Y_annual_1_guess = 0.3; KY_ratio_start_year = 1.994; Bg_Y_ratio_start_year = 0.08;
paths.K_physical(1) = KY_ratio_start_year * Y_annual_1_guess;
Y_annual_1 = cS.A * paths.K_physical(1)^cS.alpha * cS.L_fixed^(1-cS.alpha);
paths.B_g(1) = Bg_Y_ratio_start_year * Y_annual_1;
paths.K_pvt(1) = paths.K_physical(1) + paths.B_g(1);

% --- 从 t=1 演化到 t=2 ---
t=1;
M_1 = main_olg_v14_utils_ra.get_prices_and_policy_at_t_RA(t, paths.K_pvt(1), paths.B_g(1), A_path(1), cS);
[K_pvt_2, C_1] = main_olg_v14_utils_ra.simulate_RA_capital_forward(M_1, cS);
paths.K_pvt(2) = K_pvt_2;

% =========================================================================
% == [核心修正] 使用正确的离散时间GBC来更新债务, 确保 t=2 的状态正确 ==
% =========================================================================
r_mkt1 = M_1.r_mkt_annual;
Tax_l1 = cS.tau_l * M_1.w_t * M_1.L_t;
Tax_k1 = cS.tau_k * r_mkt1 * M_1.K_physical_t;
Tax_c1 = cS.tau_c * C_1;
TotalTax1 = Tax_l1 + Tax_k1 + Tax_c1;
G1 = cS.gov_exp_frac_Y * M_1.Y_t;
paths.B_g(2) = paths.B_g(1) * (1 + r_mkt1) + G1 - TotalTax1;

fprintf('✅ 已获得 t=2 的初始状态。\n\n');

% =========================================================================
% == 代表性代理(RA)模型的最终审计 ==
% =========================================================================
%% 3. [最终审计] 剖析 t=2 的会计闭环 (RA版)
t = 2;
fprintf('--- 3. 开始对 t = %d (年份 %d) 进行RA模型审计 ---\n', t, sim_years(t));
fprintf('--------------------------------------------------------\n\n');

% --- 步骤 A: 获取所有状态、价格和决策 ---
K_pvt_t = paths.K_pvt(t);
B_g_t = paths.B_g(t);
M_t = main_olg_v14_utils_ra.get_prices_and_policy_at_t_RA(t, K_pvt_t, B_g_t, A_path(t), cS);
[K_pvt_next, C_t] = main_olg_v14_utils_ra.simulate_RA_capital_forward(M_t, cS);

% --- 步骤 B: [带来源标注] 构造国民经济核算总账 ---
fprintf('--- 国民经济核算总账 (t=%d, RA模型) ---\n', t);
fprintf('%-35s | %15s | %s\n', '项目', '金额 (每期总量)', '计算来源/函数');
fprintf('%s\n', repmat('-', 75, 1));

% 1. 供给与收入检验 (Y_gross vs. Factor Payments + Depreciation)
Y_t = M_t.Y_t;
r_mkt = M_t.r_mkt_annual;
Depreciation = cS.ddk * M_t.K_physical_t;
% [核心] 总收入 = 工资 + 资本毛回报 (MPK * K)
GrossFactorIncome = M_t.w_t * M_t.L_t + (r_mkt + cS.ddk) * M_t.K_physical_t;
Supply_Income_Residual = Y_t - GrossFactorIncome;
fprintf('%-35s | %15.4e | %s\n', 'A. 供给-收入检验 (Y-wL-MPK*K)', Supply_Income_Residual, '欧拉定理');
fprintf('   - 总供给 (GDP)             | %15.4f | M_t.Y_t\n', Y_t);
fprintf('   - 总要素收入 (wL+MPK*K)    | %15.4f | wL+(r+d)K\n', GrossFactorIncome);
fprintf('\n');

% 2. 家庭预算检验 (来源与运用)
NetLaborIncome_flow = M_t.w_t * M_t.L_t * (1 - cS.tau_l);
NetCapitalIncome_flow = r_mkt * M_t.K_physical_t * (1 - cS.tau_k);
BondIncome_flow = r_mkt * M_t.B_g_t;
DisposableIncome_flow = NetLaborIncome_flow + NetCapitalIncome_flow + BondIncome_flow;
TotalSources = M_t.K_pvt_t + DisposableIncome_flow;
TotalUses = K_pvt_next + C_t * (1 + cS.tau_c);
HouseholdBudgetResidual = TotalSources - TotalUses;
fprintf('%-35s | %15.4e | %s\n', 'B. 家庭预算检验 (来源-运用)', HouseholdBudgetResidual, '手动审计计算');
fprintf('   - 资金总来源 (K_t + Yd)    | %15.4f | K_pvt_t + Disp. Income\n', TotalSources);
fprintf('   - 资金总运用 (K_{t+1} + C_exp) | %15.4f | K_pvt_next + Cons. Exp.\n', TotalUses);
fprintf('\n');

% 3. 最终检验：市场出清 (Y_gross = C + I_gross + G)
G_t = cS.gov_exp_frac_Y * Y_t;
Tax_l = cS.tau_l * M_t.w_t * M_t.L_t;
Tax_k = cS.tau_k * r_mkt * M_t.K_physical_t;
Tax_c = cS.tau_c * C_t;
TotalTax = Tax_l + Tax_k + Tax_c;
B_g_next = B_g_t * (1 + r_mkt) + G_t - TotalTax;
K_physical_next = K_pvt_next - B_g_next;
I_gross = K_physical_next - M_t.K_physical_t + Depreciation;
Market_Clearing_Residual = Y_t - (C_t + I_gross + G_t);
fprintf('%-35s | %15.4e | %s\n', 'C. 市场出清残差 (Y-C-I_gross-G)', Market_Clearing_Residual, '最终检验');

if abs(Market_Clearing_Residual) < 1e-9 && abs(Supply_Income_Residual) < 1e-9 && abs(HouseholdBudgetResidual) < 1e-9
    fprintf('\n✅ 审计通过！RA模型所有会计恒等式精确闭合！\n');
else
    fprintf('\n❌ 审计失败！请检查仍然存在的差异。\n');
end