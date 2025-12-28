% =========================================================================
% == SCRIPT: test_diffN_ss.m
% == 版本: [v2.0 - 会计逻辑对齐最终版]
% ==
% == 目的:
% ==   - [!!! 核心修正 !!!] 删除了脚本末尾所有手动、重复且不精确的总量
% ==     聚合代码。现在完全依赖由 SS.solve_steady_state 返回的、内部
% ==     会计一致的 ss0 和 ssF 结构体。
% ==   - 这确保了保存到 .mat 文件中的 ss0 和 ssF 结构体自身就是完备的，
% ==     包含了转轨求解器所需的所有精确会计变量（如 Bequest_gen_hat_raw_ss）。
% ==   - 修正了人口路径生成函数 local_generate_pop_path... 中的一个
% ==     逻辑错误，以提高其在各种冲击情景下的稳健性。
% =========================================================================
clear; close all; clear classes;
addpath(pwd);
fprintf('=== 受控人口动态稳态求解脚本 (v2.0 - 会计逻辑对齐最终版) ===\n\n');

%% --- 1. 实验设置与参数初始化 ---
fprintf('--- 1. 实验设置与参数初始化 ---\n');
output_filename = 'SS/data_for_diff_transition.mat';
max_sim_periods = 200; % 模拟期数
convergence_tol = 1e-8; % 收敛容忍度

cS = utils.ParameterValues();
cS.pps_active = false;
cS.nk = 40; cS.nkpps = 1; cS.nkprime = 40;
cS = utils.generateGrids(cS);
% 此处预先计算 nw_expanded 是为了方便，后面会为ss0和ssF分别重新计算
[~,~,~, cS.nw_expanded] = utils.EarningProcess_AgeDependent(cS); 

n_ss0_annual = 0.005; % 初始稳态的人口年增长率
fprintf('   初始稳态人口年增长率 n_ss0 = %.4f\n', n_ss0_annual);


%% --- 2. [核心] 构建由出生率冲击驱动的人口转轨路径 ---
fprintf('\n--- 2. 构建由【出生率冲击】驱动的人口转轨路径 ---\n');

% --- 冲击设置 ---
shock_start_period = 3;           % 冲击在第10个模型期开始
shock_end_period = 4;             % 冲击在第11个模型期结束
birth_rate_multiplier = 0.99;      % 冲击期间，出生率变为原来的99%

% 构造出生率冲击乘子路径
birth_rate_multiplier_path = ones(1, max_sim_periods);
birth_rate_multiplier_path(shock_start_period : shock_end_period) = birth_rate_multiplier;
fprintf('   [冲击设定] 出生率将在第 %d-%d 期变为原水平的 %.2f 倍。\n', shock_start_period, shock_end_period, birth_rate_multiplier);

% 调用本地函数生成人口路径
[Z_path_raw, Z_path_norm, n_ssF_annual, T_sim] = local_generate_pop_path_from_birth_shock(...
    n_ss0_annual, cS.s_pathV(:, end), max_sim_periods, convergence_tol, cS, ...
    birth_rate_multiplier_path);

fprintf('   ✅ 人工人口路径生成完毕。\n');
fprintf('   内生决定的收敛期数 T_sim = %d\n', T_sim);
fprintf('   终期稳态人口年增长率 n_ssF = %.4f (预期值应回归到 %.4f)\n', n_ssF_annual, n_ss0_annual);
if abs(n_ssF_annual - n_ss0_annual) > 1e-7
    warning('检验提示：终期增长率(%.6f)与初始增长率(%.6f)未完全一致。', n_ssF_annual, n_ss0_annual);
end

% 将生成的路径存入 cS 结构体
cS.T_sim = T_sim;
cS.Z_path_raw = Z_path_raw;
cS.s_pathV = repmat(cS.s_pathV(:,end), 1, T_sim);
cS.Z_path = Z_path_norm;
cS.A_path = ones(1, T_sim);
cS.theta_path = repmat(cS.theta_path(end), 1, T_sim);
cS.g_A_ss = 0.0; % 假设技术无增长以隔离人口动态效应


%% --- 3. 求解初始稳态 (ss0) ---
fprintf('\n--- 3. 求解初始稳态 (ss0) ---\n');
cS0 = cS;
cS0.n_ss = n_ss0_annual;
cS0.s_pathV = cS.s_pathV(:, 1); % 使用第一期的生存率
Z0_norm = Z_path_norm(:, 1);

paramS0 = struct();
[paramS0.leGridV, paramS0.TrProbM_by_age, paramS0.leProb1V, cS0.nw_expanded] = utils.EarningProcess_AgeDependent(cS0);

params_for_ss0 = struct(...
    'Z', Z0_norm, 'A', 1.0, 'g_A_ss', cS0.g_A_ss, 'n_ss', cS0.n_ss);

fprintf('   ⚙️  启动稳态求解器 (求解 ss0)...\n');
[ss0, dist0_norm, ~, ~] = SS.solve_steady_state(cS0, paramS0, params_for_ss0, true, 'lsqnonlin');
if isempty(ss0), error('初始稳态(ss0)求解失败！'); else fprintf('✅ 初始稳态(ss0)求解成功！\n'); end


%% --- 4. 求解终期稳态 (ssF) ---
fprintf('\n--- 4. 求解终期稳态 (ssF) ---\n');
cSF = cS;
cSF.n_ss = n_ssF_annual;
cSF.s_pathV = cS.s_pathV(:, end); % 使用最后一期的生存率
ZF_norm = Z_path_norm(:, end);

paramSF = struct();
[paramSF.leGridV, paramSF.TrProbM_by_age, paramSF.leProb1V, cSF.nw_expanded] = utils.EarningProcess_AgeDependent(cSF);

params_for_ssF = struct(...
    'Z', ZF_norm, 'A', 1.0, 'g_A_ss', cSF.g_A_ss, 'n_ss', cSF.n_ss);

fprintf('   ⚙️  启动稳态求解器 (求解 ssF)...\n');
[ssF, ~, polF, valF] = SS.solve_steady_state(cSF, paramSF, params_for_ssF, true, 'lsqnonlin');
if isempty(ssF), error('终期稳态(ssF)求解失败！'); else fprintf('✅ 终期稳态(ssF)求解成功！\n'); end


%% --- 5. 保存测试所需数据 ---
fprintf('\n--- 5. 保存转轨测试所需的数据 ---\n');
if ~exist('SS', 'dir'), mkdir('SS'); end

% 从归一化分布和第一期的绝对总人口，构造初始的绝对分布
dist0_abs = dist0_norm .* reshape(Z_path_raw(:,1), [1,1,1,cS.aD_new]);

% [!!! 核心修正: 不再手动聚合，直接使用求解器返回的完备结构体 !!!]
% ss0 和 ssF 中已经包含了转轨所需的所有会计变量，如 Bequest_gen_hat_raw_ss
data_for_diff_transition = struct(...
    'ss0', ss0, 'ssF', ssF, ...
    'dist0', dist0_abs, ...
    'polF', polF, 'valF', valF, ...
    'cS', cS, ...
    'paramS0', paramS0, 'paramSF', paramSF);

save(output_filename, 'data_for_diff_transition', '-v7.3');
fprintf('✅ 所有受控实验数据已成功保存至:\n   %s\n', output_filename);

fprintf('\n--- 受控人口动态稳态求解脚本执行完毕 ---\n');


function [Z_path_raw_final, Z_path_norm_final, n_ssF_annual, T_final] = local_generate_pop_path_from_birth_shock(n_ss0_annual, s_pathV_ss, max_sim_periods, tol, cS, birth_rate_multiplier_path)
    % =========================================================================
    % == 本地函数: local_generate_pop_path_from_birth_shock (v1.5 - 稳健收敛版)
    % == 核心修正:
    % ==   - 修正了基准出生率的计算方法，确保与动态模拟逻辑精确一致。
    % ==   - [!!!] 修正了收敛判断逻辑，确保只有在冲击完全结束后才开始
    % ==     检查分布是否收敛到新的稳态，提高了在各种情景下的稳健性。
    % =========================================================================
    
    fprintf('   本地函数 local_generate_pop_path_from_birth_shock (v1.5) 已启动...\n');
    
    % 1. 计算初始稳态的理论分布
    [Z_theory_0, ~] = utils.compute_theoretical_ss_dist(s_pathV_ss, n_ss0_annual, cS.time_Step, cS.aD_new);

    % 2. 根据稳态增长恒等式反解出与动态模拟逻辑一致的基准“动态出生率”b
    n_period_0 = (1 + n_ss0_annual)^cS.time_Step - 1;
    survivors_share_of_prev_pop_ss = sum(Z_theory_0(1:end-1) .* s_pathV_ss(1:end-1));
    baseline_birth_rate_period = (1 + n_period_0) - survivors_share_of_prev_pop_ss;
    fprintf('      根据稳态恒等式反算出精确的基准时期出生率: %.8f\n', baseline_birth_rate_period);

    % 3. 在“份额”空间进行前向模拟
    Z_norm_path = zeros(cS.aD_new, max_sim_periods);
    Z_norm_path(:, 1) = Z_theory_0;
    
    n_period_path = zeros(1, max_sim_periods - 1);
    converged = false;
    T_final = max_sim_periods;
    
    last_shock_period = find(birth_rate_multiplier_path ~= 1, 1, 'last');
    if isempty(last_shock_period), last_shock_period = 0; end

    for t = 1:(max_sim_periods - 1)
        survivors_share_of_prev_pop = Z_norm_path(1:end-1, t) .* s_pathV_ss(1:end-1);
        birth_rate_t = baseline_birth_rate_period * birth_rate_multiplier_path(t);
        growth_factor_t = birth_rate_t + sum(survivors_share_of_prev_pop);
        n_period_path(t) = growth_factor_t - 1;
        
        Z_norm_path(1, t+1) = birth_rate_t / growth_factor_t;
        Z_norm_path(2:end, t+1) = survivors_share_of_prev_pop / growth_factor_t;
        
        % [!!! 核心修正: 只有在冲击结束后才开始检查收敛 !!!]
        if ~converged && t > last_shock_period
            diff = max(abs(Z_norm_path(:, t+1) - Z_norm_path(:, t)));
            if diff < tol
                % 收敛后，再多模拟20期以确保平滑，然后裁剪
                T_final = t + 20; 
                T_final = min(T_final, max_sim_periods);
                converged = true;
                fprintf('      人口分布在第 %d 期达到新稳态，延长模拟至 %d 期。\n', t+1, T_final);
            end
        end
    end
    
    if ~converged, warning('在达到最大期数 %d 后，人口分布仍未收敛！', max_sim_periods); end

    % 4. 裁剪最终路径
    Z_path_norm_final = Z_norm_path(:, 1:T_final);
    n_period_path_final = n_period_path(1:T_final-1);

    % 5. 一次性构造绝对人口路径 Z_path_raw
    Z_path_raw_final = zeros(cS.aD_new, T_final);
    pop_total_path = ones(1, T_final); % 以1为初始总人口
    for t = 1:(T_final - 1)
        pop_total_path(t+1) = pop_total_path(t) * (1 + n_period_path_final(t));
    end
    Z_path_raw_final = Z_path_norm_final .* pop_total_path;
    
    % 6. 返回最终的、精确的稳态增长率
    n_ssF_annual = (1 + n_period_path_final(end))^(1/cS.time_Step) - 1;
end