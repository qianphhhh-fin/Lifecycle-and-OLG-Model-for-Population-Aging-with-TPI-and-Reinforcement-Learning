% =========================================================================
% == SCRIPT: Transition Path Accounting Identity Debugger (Final Calibration)
% =========================================================================
clear; close all;
addpath(pwd);
fprintf('=== 转型路径宏观会计恒等式调试脚本 (最终校准) ===\n\n');

%% 1. 初始化环境并获取初始分布
fprintf('--- 1. 初始化OLG环境并求解初始稳态分布 ---\n');
cS = main_olg_v14_utils.ParameterValues_HuggettStyle();
paramS = struct();

% --- [核心] 采用更稳健的参数校准以增强储蓄动机 ---
cS.beta = 0.985;   % [从 0.97 显著提高] 极大地增加耐心程度 (年化)
cS.sigma = 3.0;    % [设为 3.0] 更强的风险规避和平滑消费动机
cS.phi_bequest = 2.0; % [新增] 引入温和的遗赠动机，这会显著改变临终行为，增强储蓄
cS.sigma_bequest = cS.sigma;

cS.cFloor = 0.05; cS.nSim = 1000;
cS.start_year = 1997; cS.end_year = 2102; cS.time_Step = 5;
cS.alpha = 0.4; % [设为 0.4] 更符合文献的资本份额
cS.ddk = 1 - (1 - 0.05)^cS.time_Step;
cS.tau_k = 0.25; cS.tau_l = 0.20; cS.tau_c = 0.10; % 使用一些正税率
cS.pps_tax_rate_withdrawal = 0.0;
cS.gov_exp_frac_Y = 0.15; cS.rho_prime_payg = 0.4; cS.pps_activation_year = 2022;
cS.pps_return_rate_premium = 0.01;
cS.pps_withdrawal_rate = 0.15;
cS.pps_contrib_limit = 9999;
cS.pps_max_contrib_frac = 0.1;
ngrid = 20; cS.nk = ngrid; cS.nkpps = 1; cS.nkprime = ngrid; cS.npps = 1; % 关闭PPS
cS.A = 1.0; cS = main_olg_v14_utils.generateGrids(cS);
[Z_path, A_path, T_sim] = main_olg_v14_utils.load_exogenous_paths(cS);
cS = main_olg_v14_utils.calcaulte_theta_payg_path(cS, false);
sim_years = cS.start_year:cS.time_Step:(cS.start_year + (T_sim-1)*cS.time_Step);
cS.sim_years = sim_years;
[paramS.leLogGridV, paramS.leTrProbM, paramS.leProb1V] = main_olg_v14_utils.EarningProcess_olgm(cS);
paramS.leGridV = exp(paramS.leLogGridV(:));
eIdxM = main_olg_v14_utils.LaborEndowSimulation_olgm_AgeGroup(cS, paramS);

% --- [核心] 稳态求解器也需要调整 ---
% 1. 使用一个更好的初始猜测值
% 2. 遗赠不再归政府，而是返还给家庭 (因为我们引入了遗赠动机)
% 3. 我们需要一个能处理遗赠的稳态求解器版本
% 为了简化，我们暂时还用旧的求解器，但把遗赠动机作为增强储蓄的工具

% 求解初始稳态 (B_p目标为0)
% 我们将以一个更合理的 K_guess 开始
k_guess_initial = 3.0; % 资本产出比约为3是常见的起点
fprintf('使用新的 K_guess_initial = %.4f\n', k_guess_initial);

% [注意] 当前的 solve_steady_state_with_fund 假设遗赠为0。
% 引入遗赠动机后，需要一个能解均衡遗赠的求解器。
% 但我们先看看仅靠 beta 和 sigma 能否成功。
[ss_initial, eq_found, initial_dist_k, initial_dist_kpps] = main_olg_v14_utils.solve_steady_state_with_fund(Z_path(:,1), 0, cS, paramS, eIdxM);
if ~eq_found, error('初始稳态未求解成功'); end
fprintf('✅ 初始稳态和分布已获得。\n\n');


%% 2. 运行单次迭代，获取一条转型路径
% ... (这部分暂时不变，如果稳态能解出，这里就能运行) ...
fprintf('--- 2. 运行单次路径迭代以生成供审计的路径 ---\n');
K_guess_path = linspace(ss_initial.K_total, ss_initial.K_total, T_sim)'; 
w_path = zeros(T_sim, 1); r_mkt_path = zeros(T_sim, 1); b_path = zeros(T_sim, 1);
for t = 1:T_sim
    paramS_t = paramS; paramS_t.ageMassV = Z_path(:,t);
    [~, L_t] = main_olg_v14_utils.LaborSupply_Huggett(eIdxM, paramS_t, Z_path(:,t));
    [r_mkt_path(t), w_path(t), ~] = main_olg_v14_utils.HHPrices_Huggett(K_guess_path(t), L_t, A_path(t), cS);
    mass_workers_t = sum(Z_path(1:cS.aR_new, t));
    if mass_workers_t > 0, b_path(t) = cS.rho_prime_payg * (w_path(t) * L_t / mass_workers_t); else, b_path(t) = 0; end
end
PolicyFunctions = main_olg_v14_utils.solve_transition_path_policies(w_path, r_mkt_path, b_path, T_sim, cS, paramS);
Results = main_olg_v14_utils.simulate_forward_with_policies(initial_dist_k, initial_dist_kpps, PolicyFunctions, T_sim, Z_path, A_path, cS, paramS, eIdxM);
fprintf('✅ 单次路径模拟完成。\n\n');


%% 3. 逐期检查宏观会计恒等式
% ... (这部分不变) ...
fprintf('--- 3. 逐期检查宏观会计恒等式 Y = C + I + G ---\n');
fprintf('%-6s | %-12s | %-12s | %-12s | %-12s | %-12s\n', 'Year', 'Y_t', 'C_t', 'I_t', 'G_t', 'Residual');
fprintf('%s\n', repmat('-', 70, 1));
residuals = zeros(T_sim-1, 1);
for t = 1:T_sim-1
    Y_t = Results.Y_path(t); C_t = Results.C_path(t); I_t = Results.I_path(t); G_t = Results.G_path(t);
    residual = Y_t - (C_t + I_t + G_t); residuals(t) = residual;
    fprintf('%-6d | %-12.4f | %-12.4f | %-12.4f | %-12.4f | %-12.3e\n', sim_years(t), Y_t, C_t, I_t, G_t, residual);
end
if max(abs(residuals)) < 1e-6, fprintf('\n✅ 宏观会计恒等式在所有时期均成立！模型闭环。\n'); else, fprintf('\n❌ 宏观会计恒等式存在显著残差！请检查。\n'); end
figure('Name', '转型路径宏观会计恒等式残差');
plot(sim_years(1:T_sim-1), residuals, '-o');
title('Y - (C+I+G) Residual over Time'); xlabel('Year'); ylabel('Accounting Residual'); grid on;