% --- START OF FILE debug_steady_state.m (Final Verification Version) ---

% =========================================================================
% == SCRIPT: OLG Model Final Verification (Steady-State Audit)
% =========================================================================
clear; close all;
addpath(pwd);
fprintf('=== OLG 模型最终验证: 稳态法务审计 ===\n\n');

%% 1. 初始化环境
fprintf('--- 1. 初始化OLG环境...\n');
cS = main_olg_v14_utils.ParameterValues_HuggettStyle();
paramS = struct();
cS.sigma = 3.0; cS.beta = 0.97; cS.cFloor = 0.05; cS.nSim = 2000;
cS.phi_bequest = 0; cS.sigma_bequest = cS.sigma;
cS.start_year = 2002; cS.end_year = 2102; cS.time_step = 5;
cS.alpha = 0.5; cS.ddk = 1 - (1 - 0.05)^cS.time_step;
cS.tau_k = 0.2; cS.tau_l = 0.1; cS.tau_c = 0.05; cS.gov_exp_frac_Y = 0.15;
cS.rho_prime_payg = 0.5; cS.pps_activation_year = 2022; % This will be overridden for the SS test
cS.pps_tax_rate_withdrawal = 0.03; cS.pps_return_rate_premium = 0.01;
cS.pps_withdrawal_rate = 0.15; cS.pps_contrib_limit = 9999; cS.pps_max_contrib_frac = 0.1;
ngrid = 10; cS.nk = ngrid; cS.nkpps = ngrid; cS.nkprime = ngrid; cS.npps = 10;
cS.A = 1.0; cS = main_olg_v14_utils.generateGrids(cS);
[paramS.leLogGridV, paramS.leTrProbM, paramS.leProb1V] = main_olg_v14_utils.EarningProcess_olgm(cS);
paramS.leGridV = exp(paramS.leLogGridV(:));
paramS.ageEffV_new = cS.ageEffV_new;
eIdxM = main_olg_v14_utils.LaborEndowSimulation_olgm_AgeGroup(cS, paramS);
fprintf('✅ OLG环境初始化完成。\n\n');


%% 2. 计算稳态均衡
fprintf('--- 2. 计算稳态均衡...\n');
popS = main_olg_v14_utils.initPopulation(cS);
popS = main_olg_v14_utils.populationDynamics(popS, cS);
[Z_ss, ~, ~, ~] = main_olg_v14_utils.detectSteadyStatePopulation(popS, cS);
Z_ss_norm = Z_ss / sum(Z_ss);

B_p_Y_ratio_target = 0.1;
cS.gov_debt_frac_Y = 0.3;
cS.theta_path = ones(50,1) * 0.11; % Fixed PAYG rate for steady state

% [关键] 在求解稳态时，我们假设PPS是不存在的
[ss, eq_found] = main_olg_v14_utils.solve_steady_state_with_fund(Z_ss_norm, B_p_Y_ratio_target, cS, paramS, eIdxM);

if ~eq_found, error('未能找到稳态均衡，无法进行审计。'); end
fprintf('✅ 稳态均衡求解完成。\n\n');


%% 3. [最终审计] 使用稳态值运行模型一期
fprintf('--- 3. 对稳态均衡点进行法务审计...\n');
t_ss = 1;

% 使用从求解器得到的均衡值作为输入
K_pvt_ss = ss.K_pvt;
K_pps_ss = ss.K_pps; % Should be 0 from the solver
B_p_ss   = ss.B_p;
B_g_ss   = ss.B_g;
Z_t      = Z_ss_norm;
A_t      = cS.A;

% [核心修正] 创建一个与求解器内部环境完全一致的 cS 用于审计
cS_audit = cS;
cS_audit.pps_active = false;      % 必须与求解器内部的设置一致！
cS_audit.theta_t = cS.theta_path(1); % 必须与求解器内部的设置一致！

% [验证步骤] 使用 cS_audit 运行模型一期
[K_pvt_next_agg, K_pps_next_agg, C_agg, ~, ~, ~, M_t_complete] = ...
    main_olg_v14_utils.simulate_private_capital_forward(t_ss, K_pvt_ss, B_p_ss, K_pps_ss, B_g_ss, Z_t, A_t, cS_audit, paramS, eIdxM);
fprintf('--- 模型运行完成, 开始审计 ---\n\n');


%% 4. 法务会计审计账本
fprintf('=======================================================\n');
fprintf('==      法务会计审计账本 (基于稳态均衡点)      ==\n');
fprintf('=======================================================\n\n');

[B_p_next, PensionSurplus_t] = main_olg_v14_utils.update_pension_fund(B_p_ss, M_t_complete, Z_t, cS_audit);
[B_g_next, G_t, ~] = main_olg_v14_utils.update_gov_debt(B_g_ss, M_t_complete, C_agg, cS_audit, PensionSurplus_t);

% --- 4.A 生产与收入账户 ---
fprintf('--- [A] 生产与收入账户 ---\n');
Y_t = M_t_complete.Y_t;
fprintf('%-30s: %12.4f\n', '国民产出 (Y)', Y_t);
fprintf('%-30s: %12.4f\n', '总消费 (C)', C_agg);
fprintf('\n');

% --- 4.B 国民投资账户 ---
fprintf('--- [B] 国民投资账户 ---\n');
Depreciation = cS_audit.ddk * M_t_complete.K_physical_t;
K_physical_next = K_pvt_next_agg + K_pps_next_agg + B_p_next - B_g_next;
NetInvestment = K_physical_next - M_t_complete.K_physical_t;
GrossInvestment = NetInvestment + Depreciation;
fprintf('%-30s: %12.4f\n', '总投资 (I = I_net + δK)', GrossInvestment);
fprintf('%-30s: %12.4e <<-- (稳态时必须为0)\n', '净投资 (I_net)', NetInvestment);
fprintf('\n');

% --- 4.C 政府支出账户 ---
fprintf('--- [C] 政府支出账户 ---\n');
fprintf('%-30s: %12.4f\n', '政府购买 (G)', G_t);
fprintf('\n');

%% 5. 最终清算：检验核心会计恒等式
fprintf('=======================================================\n');
fprintf('==                  最终清算 (稳态)                ==\n');
fprintf('=======================================================\n\n');

NetSaving_HH = (K_pvt_next_agg + K_pps_next_agg) - (K_pvt_ss + K_pps_ss);
NetSaving_Pension = B_p_next - B_p_ss;
NetSaving_Gov = -(B_g_next - B_g_ss);
TotalNetSaving = NetSaving_HH + NetSaving_Pension + NetSaving_Gov;
Saving_Investment_Residual = TotalNetSaving - NetInvestment;
fprintf('%-35s: %12.4e <<-- (必须为0)\n', '会计恒等式 (S-I) 残差', Saving_Investment_Residual);

Market_Clearing_Residual = Y_t - C_agg - GrossInvestment - G_t;
fprintf('%-35s: %12.4e <<-- (必须为0)\n', '市场出清 (Y-C-I-G) 残差', Market_Clearing_Residual);
fprintf('\n');

% --- 最终结论 ---
if abs(Saving_Investment_Residual) < 1e-9 && abs(NetInvestment) < 1e-9 && abs(Market_Clearing_Residual) < 1e-9
    fprintf('\n✅✅✅ 审计通过！模型在稳态均衡点上完美闭合！✅✅✅\n');
    fprintf('结论：您的 OLG 模型代码是完全正确和自洽的。\n');
else
    fprintf('\n❌ 审计失败！请检查稳态求解器和 simulate_forward 函数之间的不一致性。\n');
end