% =========================================================================
% == SCRIPT: test_diffN_ss.m
% == 版本: [v2.6 - A_path 修正与添加版]
% ==
% == 目的:
% ==   - [!!! 关键修正 !!!] 添加了技术路径 (A_path) 的定义和保存。
% ==     TPI求解器需要一个明确的 A_path 来正确计算家庭部门的跨期决策，
% ==     缺少此路径是导致 g>0 时转轨失败的关键原因之一。
% ==   - 修正了 dist0 的构造逻辑，确保其在会计上绝对正确。
% =========================================================================
clear; close all; clear classes;
addpath(pwd);
fprintf('=== 受控人口动态稳态求解脚本 (v2.6 - A_path 修正与添加版) ===\n\n');

%% --- 1. 实验设置与参数初始化 ---
fprintf('--- 1. 实验设置与参数初始化 ---\n');
output_filename = 'SS/data_for_diff_transition.mat';
max_sim_periods = 500; 
convergence_tol = 1e-10; 

cS = utils.ParameterValues();
cS.pps_active = false;
cS.nk = 60; cS.nkpps = 1;cS.n_pps_rate_grid = 1;
cS = utils.generateGrids(cS);
[~,~,~, cS.nw_expanded] = utils.EarningProcess_AgeDependent(cS); 

n_ss0_annual = 0.005; % 初始稳态的人口年增长率
cS.g_A_ss = 0.01; % [!] 设定一个正的技术进步率进行测试
fprintf('   初始稳态人口年增长率 n_ss0 = %.4f\n', n_ss0_annual);
fprintf('   技术进步年增长率 g_A_ss = %.4f\n', cS.g_A_ss);


%% --- 2. [核心] 构建由出生率冲击驱动的人口转轨路径 ---
fprintf('\n--- 2. 构建由【出生率冲击】驱动的人口转轨路径 ---\n');

shock_start_period = 1;
shock_end_period = 2;
birth_rate_multiplier = 0.8;

birth_rate_multiplier_path = ones(1, max_sim_periods);
birth_rate_multiplier_path(shock_start_period : shock_end_period) = birth_rate_multiplier;
fprintf('   [冲击设定] 出生率将在第 %d-%d 期变为原水平的 %.2f 倍。\n', shock_start_period, shock_end_period, birth_rate_multiplier);

[Z_path_raw, Z_path_norm, n_ssF_annual, T_sim] = local_generate_pop_path_from_birth_shock(...
    n_ss0_annual, cS.s_pathV(:, end), max_sim_periods, convergence_tol, cS, ...
    birth_rate_multiplier_path);

fprintf('   ✅ 人工人口路径生成完毕。\n');
fprintf('   内生决定的收敛期数 T_sim = %d\n', T_sim);
fprintf('   终期稳态人口年增长率 n_ssF = %.4f (预期值应回归到 %.4f)\n', n_ssF_annual, n_ss0_annual);
if abs(n_ssF_annual - n_ss0_annual) > 1e-7
    warning('检验提示：终期增长率(%.6f)与初始增长率(%.6f)未完全一致。', n_ssF_annual, n_ss0_annual);
end

cS.T_sim = T_sim;
cS.Z_path_raw = Z_path_raw;
cS.s_pathV = repmat(cS.s_pathV(:,end), 1, T_sim);
cS.Z_path = Z_path_norm;
cS.theta_path = repmat(cS.theta_path(end), 1, T_sim);

% [!!! 核心修正: 定义并添加A_path !!!]
g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
cS.A_path = (1 + g_A_period).^(0:(T_sim-1));
fprintf('   ✅ 已构建基于 g_A_ss 的技术水平路径 A_path (T=%d)。\n', T_sim);


%% --- 3. 求解初始稳态 (ss0) ---
fprintf('\n--- 3. 求解初始稳态 (ss0) ---\n');
cS0 = cS;
cS0.n_ss = n_ss0_annual;
cS0.s_pathV = cS.s_pathV(:, 1);
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
if ~exist('TRANS', 'dir'), mkdir('TRANS'); end
if ~exist('SS', 'dir'), mkdir('SS'); end

Z0_norm_reshaped = reshape(Z0_norm, [1,1,1,cS.aD_new]);
Z0_norm_safe = Z0_norm_reshaped;
Z0_norm_safe(Z0_norm_safe < 1e-12) = 1; 
dist_cond_ss0 = dist0_norm ./ Z0_norm_safe;

Z_t1_abs_reshaped = reshape(Z_path_raw(:,1), [1,1,1,cS.aD_new]);
dist0_abs = dist_cond_ss0 .* Z_t1_abs_reshaped;
fprintf('   [信息] 已使用正确的会计方法重新构造 t=1 的初始分布(dist0)。\n');
fprintf('   检查: 初始稳态总人口 (sum(dist0_norm)): %.4f\n', sum(dist0_norm, 'all'));
fprintf('   检查: 转轨期t=1总人口 (sum(dist0_abs)):  %.4f\n', sum(dist0_abs, 'all'));


data_for_diff_transition = struct(...
    'ss0', ss0, 'ssF', ssF, ...
    'dist0', dist0_abs, ...
    'polF', polF, 'valF', valF, ...
    'cS', cS, ...
    'paramS0', paramS0, 'paramSF', paramSF);

save(output_filename, 'data_for_diff_transition', '-v7.3');
fprintf('✅ 所有受控实验数据已成功保存至:\n   %s\n', output_filename);


%% --- 6. [新增] 比较初始与终期稳态 ---
fprintf('\n--- 6. 比较初始稳态(ss0)与终期稳态(ssF) ---\n');
vars_to_compare = {'r_mkt', 'w_hat', 'K_private_hat', 'L_hat', 'C_agg', 'b_hat', 'Bequest_generated_agg', 'TR_distributed_agg'};
fprintf('%-25s | %15s | %15s | %15s | %12s\n', '变量名', '初始稳态 (ss0)', '终期稳态 (ssF)', '绝对差异', '相对差异 (%)');
fprintf(repmat('-', 1, 80));
fprintf('\n');
for i = 1:length(vars_to_compare)
    var_name = vars_to_compare{i};
    if isfield(ss0, var_name) && isfield(ssF, var_name)
        val0 = ss0.(var_name);
        valF = ssF.(var_name);
        abs_diff = valF - val0;
        if abs(val0) > 1e-9, rel_diff_pct = (abs_diff / val0) * 100; else, rel_diff_pct = NaN; end
        fprintf('%-25s | %15.6f | %15.6f | %15.6e | %12.4e\n', var_name, val0, valF, abs_diff, rel_diff_pct);
    else
        fprintf('%-25s | <变量不存在>\n', var_name);
    end
end
fprintf(repmat('-', 1, 80));
fprintf('\n');

fprintf('\n--- 受控人口动态稳态求解脚本执行完毕 ---\n');


function [Z_path_raw_final, Z_path_norm_final, n_ssF_annual, T_final] = local_generate_pop_path_from_birth_shock(n_ss0_annual, s_pathV_ss, max_sim_periods, tol, cS, birth_rate_multiplier_path)
    % local_generate_pop_path_from_birth_shock (v1.5)
    fprintf('   本地函数 local_generate_pop_path_from_birth_shock (v1.5) 已启动...\n');
    [Z_theory_0, ~] = population.compute_theoretical_ss_dist(s_pathV_ss, n_ss0_annual, cS.time_Step, cS.aD_new);
    n_period_0 = (1 + n_ss0_annual)^cS.time_Step - 1;
    survivors_share_of_prev_pop_ss = sum(Z_theory_0(1:end-1) .* s_pathV_ss(1:end-1));
    baseline_birth_rate_period = (1 + n_period_0) - survivors_share_of_prev_pop_ss;
    fprintf('      根据稳态恒等式反算出精确的基准时期出生率: %.8f\n', baseline_birth_rate_period);
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
        if ~converged && t > last_shock_period
            diff = max(abs(Z_norm_path(:, t+1) - Z_norm_path(:, t)));
            if diff < tol
                T_final = t + 10; 
                T_final = min(T_final, max_sim_periods);
                converged = true;
                fprintf('      人口分布在第 %d 期达到新稳态，延长模拟至 %d 期。\n', t+1, T_final);
            end
        end
    end
    if ~converged, warning('在达到最大期数 %d 后，人口分布仍未收敛！', max_sim_periods); end
    Z_path_norm_final = Z_norm_path(:, 1:T_final);
    n_period_path_final = n_period_path(1:T_final-1);
    Z_path_raw_final = zeros(cS.aD_new, T_final);
    pop_total_path = ones(1, T_final);
    for t = 1:(T_final - 1)
        pop_total_path(t+1) = pop_total_path(t) * (1 + n_period_path_final(t));
    end
    Z_path_raw_final = Z_path_norm_final .* pop_total_path;
    n_ssF_annual = (1 + n_period_path_final(end))^(1/cS.time_Step) - 1;
end