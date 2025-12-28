% --- START OF FILE debug_final_identity.m ---
% =========================================================================
% == SCRIPT: Final Identity Decomposition Debugger
% == 目的：通过分解国民收入恒等式，最终定位会计不一致的来源。
% =========================================================================
clear; close all;
addpath(pwd);
fprintf('=== 国民收入恒等式分解调试脚本 ===\n\n');

%% 1. 初始化一个完全相同的环境
fprintf('--- 1. 初始化OLG环境 (零税率调试模式) ---\n');
cS = main_olg_v14_utils.ParameterValues_HuggettStyle();
paramS = struct();
cS.sigma = 3.0; cS.beta = 0.97; cS.cFloor = 0.05; cS.nSim = 1000;
cS.phi_bequest = 3; cS.sigma_bequest = cS.sigma;
cS.start_year = 1997; cS.end_year = 2102; cS.time_step = 5;
cS.alpha = 0.5; cS.ddk = 1 - (1 - 0.05)^cS.time_step;
cS.tau_k = 0.0; cS.tau_l = 0.0; cS.tau_c = 0.0; cS.pps_tax_rate_withdrawal = 0.0;
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

%% 2. [核心诊断] 分解国民收入恒等式
t = 2;
fprintf('--- 2. 剖析 t=%d 的国民收入恒等式 ---\n', t);
K_pvt_t = 1.5; K_pps_t = 0.1; B_p_t = 0.2; B_g_t = 0.3;
TR_total_t = 0.01; Z_t = Z_path(:, t); A_t = A_path(t);
cS.pps_active = (sim_years(t) >= cS.pps_activation_year);

% --- 步骤 A: 运行模型一个周期 ---
M_t = main_olg_v14_utils.get_prices_and_policy_at_t(t, K_pvt_t, B_p_t, K_pps_t, B_g_t, Z_t, A_t, cS, paramS, eIdxM);
[K_pvt_next, K_pps_next, C_t, Total_Cpps_t, Total_PpsTax_t, ~] = ...
    main_olg_v14_utils.simulate_private_capital_forward(M_t, Z_t, cS, paramS, eIdxM, TR_total_t);
[B_p_next, PensionSurplus_t] = main_olg_v14_utils.update_pension_fund(B_p_t, M_t, Z_t, cS);
[B_g_next, G_t, TotalTax_t] = main_olg_v14_utils.update_gov_debt(B_g_t, C_t, M_t, Total_Cpps_t, K_pps_t, Total_PpsTax_t, cS, PensionSurplus_t);

% --- 步骤 B: [恒等式分解] ---
fprintf('\n--- 恒等式: S_pvt + S_g + S_p = I_gross ---\n');
fprintf('%-35s | %15s \n', '项目', '金额');
fprintf('%s\n', repmat('-', 55, 1));

% --- a. 计算恒等式左侧 (储蓄来源) ---
S_pvt = (K_pvt_next + K_pps_next) - (K_pvt_t + K_pps_t);
S_g = TotalTax_t - G_t - (B_g_next - B_g_t - M_t.r_mkt_period * B_g_t); % 政府储蓄 = T - G - 新增净借款
S_p = PensionSurplus_t - (B_p_next - B_p_t - M_t.r_mkt_period * B_p_t); % 养老金储蓄 = 盈余 - 新增净借款
LHS = S_pvt + S_g + S_p;
fprintf('%-35s | %15.4f \n', '私人储蓄 (S_pvt)', S_pvt);
fprintf('%-35s | %15.4f \n', '政府储蓄 (S_g)', S_g);
fprintf('%-35s | %15.4f \n', '养老金储蓄 (S_p)', S_p);
fprintf('%-35s | %15.4f \n', '国民总储蓄 (LHS)', LHS);
fprintf('%s\n', repmat('-', 55, 1));

% --- b. 计算恒等式右侧 (投资去向) ---
K_physical_t = M_t.K_physical_t;
K_physical_next = K_pvt_next + K_pps_next + B_p_next - B_g_next;
I_gross = K_physical_next - K_physical_t + cS.ddk * K_physical_t;
RHS = I_gross;
fprintf('%-35s | %15.4f \n', '国民总投资 (RHS)', RHS);
fprintf('%s\n', repmat('-', 55, 1));

% --- c. 计算最终残差 ---
Final_Residual = LHS - RHS;
fprintf('%-35s | %15.4e \n', '最终残差 (LHS - RHS)', Final_Residual);

if abs(Final_Residual) < 1e-6
    fprintf('\n✅ 恒等式检验通过！模型会计闭环。\n');
else
    fprintf('\n❌ 恒等式检验失败！请仔细检查上面每一项的计算。\n');
end