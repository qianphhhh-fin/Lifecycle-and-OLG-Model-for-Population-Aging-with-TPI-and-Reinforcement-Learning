% =========================================================================
% == SCRIPT: Transition Path Accounting Identity Debugger (FIXED)
% == 目的：在单次迭代产生的转型路径上，逐期检查宏观会计恒等式
% == Y = C + I + G 是否成立。
% =========================================================================
clear; close all;
addpath(pwd);
fprintf('=== 转型路径宏观会计恒等式调试脚本 ===\n\n');

%% 1. 初始化环境并获取初始分布
fprintf('--- 1. 初始化OLG环境并求解初始稳态分布 ---\n');
cS = main_olg_v14_utils.ParameterValues_HuggettStyle();
paramS = struct();
 
ngrid = 30; cS.nk = ngrid; cS.nkpps = ngrid; cS.nkprime = ngrid; cS.npps = 5;
cS = main_olg_v14_utils.generateGrids(cS);
[Z_path, A_path, T_sim] = main_olg_v14_utils.load_exogenous_paths(cS);
cS = main_olg_v14_utils.calcaulte_theta_payg_path(cS, false);
sim_years = cS.start_year:cS.time_Step:(cS.start_year + (T_sim-1)*cS.time_Step);
cS.sim_years = sim_years;
[paramS.leLogGridV, paramS.leTrProbM, paramS.leProb1V] = main_olg_v14_utils.EarningProcess_olgm(cS);
paramS.leGridV = exp(paramS.leLogGridV(:));
eIdxM = main_olg_v14_utils.LaborEndowSimulation_olgm_AgeGroup(cS, paramS);

% 求解初始稳态 (假设遗赠归政府，B_p目标为0)
[ss_initial, eq_found, initial_dist_k, initial_dist_kpps] = main_olg_v14_utils.solve_steady_state_with_fund(Z_path(:,1), 0, cS, paramS, eIdxM);
if ~eq_found, error('初始稳态未求解成功'); end
fprintf('✅ 初始稳态和分布已获得。\n\n');


%% 2. 运行单次迭代，获取一条转型路径
fprintf('--- 2. 运行单次路径迭代以生成供审计的路径 ---\n');

% a. 猜测一个资本路径 (简单线性插值)
K_guess_path = linspace(ss_initial.K_total, ss_initial.K_total * 1.2, T_sim)'; % 猜测一个变化的路径更有趣

% b. 计算价格路径
w_path = zeros(T_sim, 1);
r_mkt_path = zeros(T_sim, 1);
b_path = zeros(T_sim, 1);
L_path = zeros(T_sim, 1); % <--- 【新增】存储劳动供给路径
for t = 1:T_sim
    % ---【核心修正】---
    % 修正了 LaborSupply_Huggett 的调用参数，并删除了多余的 paramS_t
    [~, L_t] = main_olg_v14_utils.LaborSupply_Huggett(eIdxM, cS, paramS, Z_path(:,t));
    L_path(t) = L_t; % <--- 【新增】存储
    
    M_prices_t = main_olg_v14_utils.get_prices_at_t(K_guess_path(t), L_t, A_path(t), cS);
    r_mkt_path(t) = M_prices_t.r_mkt_t;
    w_path(t) = M_prices_t.w_t;

    mass_workers_t = sum(Z_path(1:cS.aR_new, t));
    if mass_workers_t > 0
        b_path(t) = cS.rho_prime_payg * (w_path(t) * L_t / mass_workers_t);
    else
        b_path(t) = 0;
    end
end

% c. 反向求解所有策略函数
PolicyFunctions = main_olg_v14_utils.solve_transition_path_policies(w_path, r_mkt_path, b_path, T_sim, cS, paramS);
save('PolicyFunctions.mat', 'PolicyFunctions');

% load('PolicyFunctions.mat')
% d. 前向模拟，得到宏观路径结果
% <--- 【核心修改】将价格路径传入模拟函数
Results = main_olg_v14_utils.simulate_forward_with_policies(initial_dist_k, initial_dist_kpps, PolicyFunctions, T_sim, Z_path, A_path, w_path, r_mkt_path, b_path, L_path, cS, paramS, eIdxM);
fprintf('✅ 单次路径模拟完成。\n\n');


%% 3. 逐期检查宏观会计恒等式
fprintf('--- 3. 逐期检查宏观会计恒等式 Y = C + I + G ---\n');
fprintf('%-6s | %-12s | %-12s | %-12s | %-12s | %-12s\n', ...
    'Year', 'Y_t', 'C_t', 'I_t', 'G_t', 'Residual');
fprintf('%s\n', repmat('-', 70, 1));

residuals = zeros(T_sim-1, 1);
for t = 1:T_sim-1
    Y_t = Results.Y_path(t);
    C_t = Results.C_path(t);
    I_t = Results.I_path(t);
    G_t = Results.G_path(t);
    
    residual = Y_t - (C_t + I_t + G_t);
    residuals(t) = residual;
    
    fprintf('%-6d | %-12.4f | %-12.4f | %-12.4f | %-12.4f | %-12.3e\n', ...
        sim_years(t), Y_t, C_t, I_t, G_t, residual);
end

% 将容忍度放宽一点，因为数值插值会带来微小误差
if max(abs(residuals)) < 1e-9 
    fprintf('\n✅ 宏观会计恒等式在所有时期均成立！模型闭环。\n');
else
    fprintf('\n❌ 宏观会计恒等式存在显著残差！请检查。\n');
    fprintf('   最大残差为: %.3e\n', max(abs(residuals)));
end

figure('Name', '转型路径宏观会计恒等式残差');
plot(sim_years(1:T_sim-1), residuals, '-o');
title('Y - (C+I+G) Residual over Time');
xlabel('Year');
ylabel('Accounting Residual');
grid on;