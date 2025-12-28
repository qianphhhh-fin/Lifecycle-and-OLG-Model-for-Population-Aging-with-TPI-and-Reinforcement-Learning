% =========================================================================
% == SCRIPT: test_diffN_ss.m
% == 版本: [v1.5 - BGP模拟逻辑修正版]
% ==
% == 目的:
% ==   - 修正 local_generate_pop_path 的内部模拟逻辑，使其与BGP的
% ==     理论定义完全一致，确保在无冲击时能精确维持初始增长率。
% ==   - 关键修改：在计算新生儿数量时，显式地引入(1+n)增长因子。
% =========================================================================
clear; close all; clear classes;
addpath(pwd);
fprintf('=== 受控人口动态稳态求解脚本 (test_diffN_ss.m) ===\n\n');

%% --- 1. 实验设置与参数初始化 ---
fprintf('--- 1. 实验设置与参数初始化 ---\n');
output_filename = 'SS/data_for_diff_transition.mat';
max_sim_periods = 1000; 
convergence_tol = 1e-7; 

cS = utils.ParameterValues();
cS.pps_active = false;
cS.nk = 40; cS.nkpps = 1; cS.nkprime = 40;
cS = utils.generateGrids(cS);
[~,~,~, cS.nw_expanded] = utils.EarningProcess_AgeDependent(cS);



n_ss0_annual = 0.005;
fprintf('   初始稳态人口年增长率 n_ss0 = %.4f\n', n_ss0_annual);


%% --- 2. [核心] 构建人工控制的人口转轨路径 ---
fprintf('\n--- 2. 构建人工控制的人口转轨路径 ---\n');


% 冲击设置为1，进行无冲击的基准测试
shock_start_period = 1;
birth_rate_multiplier_shock = 0.99; 

[Z_path_raw, Z_path_norm, n_ssF_annual, T_sim] = local_generate_pop_path(...
    n_ss0_annual, cS.s_pathV(:, end), max_sim_periods, convergence_tol, cS, ...
    shock_start_period, birth_rate_multiplier_shock);

fprintf('   ✅ 人工人口路径生成完毕。\n');
fprintf('   内生决定的收敛期数 T_sim = %d\n', T_sim);
fprintf('   终期稳态人口年增长率 n_ssF = %.4f\n', n_ssF_annual);
% if abs(n_ssF_annual - n_ss0_annual) < 1e-6
%     fprintf('   ✅ 无冲击检验通过：终期增长率与初始增长率一致。\n');
% else
%     fprintf('无冲击检验失败！终期增长率(%.4f)与初始增长率(%.4f)不一致。', n_ssF_annual, n_ss0_annual);
% end


% 将生成的路径存入 cS 结构体
cS.T_sim = T_sim;
cS.Z_path_raw = Z_path_raw;
cS.s_pathV = repmat(cS.s_pathV(:,end), 1, T_sim);
cS.Z_path = Z_path_norm;
cS.A_path = ones(1, T_sim); 
cS.theta_path = repmat(cS.theta_path(end), 1, T_sim);
cS.g_A_ss = 0.0; 


%% --- 3. 求解初始稳态 (ss0) ---
fprintf('\n--- 3. 求解初始稳态 (ss0) ---\n');
cS0 = cS;
cS0.n_ss = n_ss0_annual;
cS0.s_pathV = cS.s_pathV(:, end);

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
cSF.s_pathV = cS.s_pathV(:, end);

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

dist0_abs = dist0_norm * sum(Z_path_raw(:,1));

data_for_diff_transition = struct(...
    'ss0', ss0, 'ssF', ssF, ...
    'dist0', dist0_abs, ...
    'polF', polF, 'valF', valF, ...
    'cS', cS, ... 
    'paramS0', paramS0, 'paramSF', paramSF);

save(output_filename, 'data_for_diff_transition', '-v7.3');
fprintf('✅ 所有受控实验数据已成功保存至:\n   %s\n', output_filename);

fprintf('\n--- 受控人口动态稳态求解脚本执行完毕 ---\n');


function [Z_path_raw_final, Z_path_norm_final, n_ssF_annual, T_final] = local_generate_pop_path(n_ss0_annual, s_pathV, max_sim_periods, tol, cS, shock_start, shock_multiplier)
    % 本地函数: 生成一个由出生率冲击驱动的、受控的人口转轨路径
    % [修改] v1.7 - 与 simulate_dist_forward 逻辑完全对齐
    
    fprintf('   本地函数 local_generate_pop_path (v1.7) 已启动...\n');
    
    % 1. 计算初始稳态的理论分布
    Z_theory_0 = utils.compute_theoretical_ss_dist(s_pathV, n_ss0_annual, cS.time_Step, cS.aD_new);
    
    % 2. 构建时期增长率路径。冲击直接作用于增长率本身。
    n0_period = (1 + n_ss0_annual)^cS.time_Step - 1;
    n_ssF_annual = n_ss0_annual * shock_multiplier;
    nF_period = (1 + n_ssF_annual)^cS.time_Step - 1;
    
    n_period_path = ones(1, max_sim_periods) * n0_period;
    n_period_path(shock_start : end) = nF_period;
    fprintf('      人口时期增长率将从 %.5f 在第 %d 期永久变为 %.5f。\n', n0_period, shock_start, nF_period);

    % 3. 使用与 simulate_dist_forward 相同的 BGP 会计恒等式进行模拟
    Z_path_raw = zeros(cS.aD_new, max_sim_periods);
    Z_path_raw(:, 1) = Z_theory_0; 
    
    Z_norm_prev = Z_path_raw(:, 1) / sum(Z_path_raw(:, 1));
    converged = false;
    T_final = max_sim_periods;

    for t = 1:(max_sim_periods - 1)
        % a. 计算下一期的幸存者 (逻辑与 simulate_dist_forward 一致)
        % 注意: 我们模拟的是所有年龄组，所以用 s_pathV
        survivors_tp1 = Z_path_raw(1:end-1, t) .* s_pathV(1:end-1);
        
        % b. 根据BGP会计恒等式，精确计算新生儿数量
        % Total(t+1) = Total(t) * (1 + n_period(t))
        % Total(t+1) = Newborns(t+1) + Sum(Survivors_from_t)
        total_pop_t = sum(Z_path_raw(:, t));
        total_pop_tp1 = total_pop_t * (1 + n_period_path(t));
        new_entrants_tp1 = total_pop_tp1 - sum(survivors_tp1);
        
        % 确保新生儿数量非负
        if new_entrants_tp1 < 0
            warning('在 t=%d 时计算出负的新生儿数量，检查参数。', t);
            new_entrants_tp1 = 0;
        end
        
        % c. 合并得到下一期的总人口分布
        Z_path_raw(1, t+1) = new_entrants_tp1;
        Z_path_raw(2:end, t+1) = survivors_tp1;
        
        % d. 检查收敛
        Z_norm_current = Z_path_raw(:, t+1) / sum(Z_path_raw(:, t+1));
        diff = max(abs(Z_norm_current - Z_norm_prev));
        
        if ~converged && diff < tol && t >= shock_start 
            T_final = t + 1; 
            converged = true;
            fprintf('      人口分布在第 %d 期达到新稳态。\n', T_final);
        end
        Z_norm_prev = Z_norm_current;
    end
    
    if ~converged, warning('在达到最大期数 %d 后，人口分布仍未收敛！', max_sim_periods); end

    % 4. 裁剪最终的路径矩阵
    Z_path_raw_final = Z_path_raw(:, 1:T_final);
    Z_path_norm_final = Z_path_raw_final ./ sum(Z_path_raw_final, 1);
    
    % 5. 内部一致性检验
    fprintf('      正在进行内部一致性检验...\n');
    Z_simulated_F = Z_path_norm_final(:, end);
    % 理论分布应该用我们设定的目标增长率 n_ssF_annual 来计算
    Z_theoretical_F = utils.compute_theoretical_ss_dist(s_pathV, n_ssF_annual, cS.time_Step, cS.aD_new);
    validation_error = max(abs(Z_simulated_F - Z_theoretical_F));
    fprintf('      模拟出的最终分布与理论分布的最大差值: %.3e\n', validation_error);
    if validation_error > tol % 使用与收敛相同的容忍度
        fprintf('人口模拟验证失败！模拟结果与理论不符。');
    else
        fprintf('      ✅ 检验通过：人口模拟过程内部自洽。\n');
    end
end