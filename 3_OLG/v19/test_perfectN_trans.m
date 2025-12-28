% =========================================================================
% == SCRIPT: test_perfectN_trans.m
% == 版本: [v2.5 - 最终会计对齐版]
% ==
% == 目的:
% ==   - [!!!] 使用 `ss0.Bequest_gen_hat_raw_ss` (未经BGP调整的原始总量)
% ==     来构造初始的人均遗赠猜测路径 (`bequest_hat_pc_path_guess`)。
% ==   - 这确保了初始猜测的会计口径与 TPI 引擎 (`calculate_paths_and_errors`)
% ==     内部处理 t=1 期遗赠的方式完全一致。
% ==   - 删除了所有之前版本中关于遗赠路径的多余和错误的构造逻辑。
% =========================================================================
clear; close all; clear classes;
addpath(pwd);
fprintf('=== TPI求解器零冲击诊断脚本 (v2.5 - 最终会计对齐版) ===\n\n');

%% --- 1. 加载“黄金标准”稳态数据 ---
fprintf('--- 1. 加载黄金标准稳态数据 ---\n');
data_filename = 'SS/data_for_perfect_transition.mat';
if ~exist(data_filename, 'file'), error('黄金标准稳待数据文件 "%s" 不存在。', data_filename); end
load(data_filename, 'data_for_perfect_transition');
fprintf('   ✅ 已加载数据: %s\n', data_filename);

ssF = data_for_perfect_transition.ssF;
distF = data_for_perfect_transition.distF;
valF = data_for_perfect_transition.valF;
polF = data_for_perfect_transition.polF;
cS = data_for_perfect_transition.cS;
paramSF = data_for_perfect_transition.paramSF;

cS.T_sim = 10;
T = cS.T_sim;
fprintf('   [快速诊断模式] T已被临时设为: %d\n', T);

%% --- 2. 手动构建一个完美的、自洽的、平稳增长的人口环境 ---
fprintf('\n--- 2. 手动构建完美的、平稳增长的BGP环境 ---\n');
ss0 = ssF;
dist0 = distF;
fprintf('   ✅ 初始稳态(ss0)和分布(dist0)已设为与终期稳态(ssF)相同。\n');

Z_ss_abs_dist = sum(distF, [1,2,3]);
Z_ss_norm = Z_ss_abs_dist(:) / sum(Z_ss_abs_dist(:));
total_pop_ss0 = sum(dist0, 'all');

g_n_period = (1 + cS.n_ss)^cS.time_Step - 1;
pop_growth_factors = (1 + g_n_period).^(0:T-1);
cS.Z_path_raw = (total_pop_ss0 * Z_ss_norm) * pop_growth_factors;
cS.Z_path = cS.Z_path_raw ./ sum(cS.Z_path_raw, 1);
fprintf('   ✅ 总人口路径 (Z_path_raw) 已被构建为以 n_ss=%.4f 恒定增长。\n', cS.n_ss);

s_pathV_ss = cS.s_pathV;
if size(s_pathV_ss, 2) > 1, s_pathV_ss = s_pathV_ss(:, end); end
cS.s_pathV = repmat(s_pathV_ss, 1, T);

g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
A_path_flat = (1 + g_A_period).^(0:T-1);
cS.A_path = A_path_flat;

if isfield(cS, 'theta_path'), cS.theta_path = repmat(cS.theta_path(end), 1, T);
else, cS.theta_path = zeros(1, T); end

%% --- 3. 构造初始猜测：平坦的价格路径和BGP一致的流量路径 ---
fprintf('\n--- 3. 构造BGP一致的猜测路径 ---\n');
r_path_guess = ones(1, T) * ss0.r_mkt;
w_hat_path_guess = ones(1, T) * ss0.w_hat;
b_hat_path_guess = ones(1, T) * ss0.b_hat;

% TR路径: 人均转移支付是平坦的
tr0_hat_pc = ss0.TR_distributed_agg / total_pop_ss0;
tr_hat_pc_path_guess = ones(1, T) * tr0_hat_pc;

% [!!! 核心修正: Bequest路径 !!!]
% TPI循环的迭代变量是【人均遗赠】。
% 在BGP下，这个量应该是常数。
% 我们用 ss0 中未经BGP调整的原始总量，除以总人口，得到正确的人均量。
bequest0_hat_pc = ss0.Bequest_gen_hat_raw_ss / total_pop_ss0;
bequest_hat_pc_path_guess = ones(1, T) * bequest0_hat_pc;

fprintf('   ✅ 所有价格和流量路径已构造完毕。\n');


%% --- 4. 运行TPI循环 ---
fprintf('\n--- 4. 启动TPI循环进行诊断 ---\n');
max_iter = 150;
tolerance = 1e-7;
damping.r = 0.2; damping.w = 0.2; damping.beq = 0.2; damping.tr = 0.2; damping.b_hat = 0.2;

r_path_iter = r_path_guess;
w_hat_path_iter = w_hat_path_guess;
bequest_hat_pc_path_iter = bequest_hat_pc_path_guess;
tr_hat_pc_path_iter = tr_hat_pc_path_guess;
b_hat_path_iter = b_hat_path_guess;

for iter = 1:max_iter
    [target_paths, errors, Pol_path, ~, aggr_supply] = TPI.calculate_paths_and_errors(...
        r_path_iter, w_hat_path_iter, bequest_hat_pc_path_iter, ...
        tr_hat_pc_path_iter, b_hat_path_iter, ...
        ss0, ssF, cS, paramSF, valF, polF, dist0);

    policy_error = 0;
    % 我们只关心第一次迭代的策略误差，之后会因宏观误差而自然偏离
    if iter == 1
        pol_t1 = Pol_path{1};
        for ia = 1:cS.aD_new
            err_c = max(abs(pol_t1(ia).c - polF(ia).c), [], 'all');
            err_k_prime = max(abs(pol_t1(ia).k_prime - polF(ia).k_prime), [], 'all');
            policy_error = max([policy_error, err_c, err_k_prime]);
            if isfield(polF(ia), 'kpps_prime') && ~isempty(polF(ia).kpps_prime)
                 err_kpps_prime = max(abs(pol_t1(ia).kpps_prime - polF(ia).kpps_prime), [], 'all');
                 policy_error = max(policy_error, err_kpps_prime);
            end
        end
    end

    fprintf('Iter [%2d/%2d] | Error: %.2e (r:%.1e, w:%.1e, bq:%.1e, tr:%.1e, b:%.1e) | Pol_Err(t=1): %.2e\n', ...
        iter, max_iter, errors.total, errors.r, errors.w, errors.beq, errors.tr, errors.b_hat, policy_error);

    if errors.total < tolerance
        fprintf('\n✅ 诊断成功！宏观路径收敛。\n');
        break;
    end
    if iter == max_iter, fprintf('\n❌ 诊断失败！未在最大迭代次数内收敛。\n'); end

    % 更新猜测路径，注意遗赠的特殊性
    r_path_iter(2:T) = (1 - damping.r) * r_path_iter(2:T) + damping.r * target_paths.r_path(2:T);
    w_hat_path_iter(2:T) = (1 - damping.w) * w_hat_path_iter(2:T) + damping.w * target_paths.w_hat_path(2:T);
    tr_hat_pc_path_iter(2:T) = (1 - damping.tr) * tr_hat_pc_path_iter(2:T) + damping.tr * target_paths.TR_hat_pc_path(2:T);
    b_hat_path_iter(2:T) = (1 - damping.b_hat) * b_hat_path_iter(2:T) + damping.b_hat * target_paths.b_hat_path(2:T);
    
    % [!!!] 更新遗赠路径，根据修正后的 target_paths 定义
    bequest_hat_pc_path_iter(2:T) = (1 - damping.beq) * bequest_hat_pc_path_iter(2:T) + damping.beq * target_paths.bequest_hat_pc_path(2:T);
end

%% --- 5. 结果分析 ---
if errors.total < tolerance
    fprintf('\n--- 5. 最终路径检验 ---\n');
    max_dev_r = max(abs(target_paths.r_path - ssF.r_mkt));
    max_dev_w = max(abs(target_paths.w_hat_path - ssF.w_hat));
    fprintf('   利率路径与稳态的最大偏差: %.3e\n', max_dev_r);
    fprintf('   工资路径与稳态的最大偏差: %.3e\n', max_dev_w);
    
    % 检查最终的人均遗赠路径是否等于我们构造的初始值
    max_dev_beq_pc = max(abs(target_paths. bequest_gen_raw_pc_path - bequest0_hat_pc));
    fprintf('   人均遗赠路径与BGP人均值的最大偏差: %.3e\n', max_dev_beq_pc);

    if max_dev_r < 1e-5 && max_dev_w < 1e-5 && max_dev_beq_pc < 1e-5
        fprintf('   ✅ 结论：求解出的路径确实是平坦的稳态路径。TPI核心逻辑通过检验。\n');
    else
        fprintf('   ⚠️ 警告：路径收敛了，但并非平坦路径。可能存在微小的不一致性。\n');
    end
end

%% --- 6. [新增] 最终国民账户核算检验 ---
fprintf('\n--- 6. 对最终收敛路径进行国民账户核算检验 ---\n');
if errors.total < tolerance
    
    % --- 步骤 6a: 使用收敛的路径，重新计算一次最终的宏观量 ---
    % 这一步确保我们使用的是最精确的聚合结果进行检验
    [final_paths, ~, ~, ~, final_aggr_supply] = TPI.calculate_paths_and_errors(...
        r_path_iter, w_hat_path_iter, bequest_hat_pc_path_iter, ...
        tr_hat_pc_path_iter, b_hat_path_iter, ...
        ss0, ssF, cS, paramSF, valF, polF, dist0);

    % --- 步骤 6b: 整理所有宏观路径到 results 结构体中 ---
    results = struct();
    total_pop_path = sum(cS.Z_path_raw, 1);
    
    % 价格与水平量路径
    results.r_path = final_paths.r_path;
    results.w_path = final_paths.w_hat_path .* cS.A_path;
    K_p_path_level = final_aggr_supply.K_p_hat_total_path .* cS.A_path;
    results.K_p_path = K_p_path_level;
    results.C_path = final_aggr_supply.C_hat_pc_path .* total_pop_path .* cS.A_path;
    
    % 产出路径
    K_g_hat_total_path = ssF.K_public_hat * total_pop_path;
    Y_hat_total_path = zeros(1, T);
    for t = 1:T
        prices_t = firm.get_prices_at_t(final_aggr_supply.K_p_hat_total_path(t), K_g_hat_total_path(t), final_aggr_supply.L_path(t), cS);
        Y_hat_total_path(t) = prices_t.Y_hat_t;
    end
    results.Y_path = Y_hat_total_path .* cS.A_path;
    
    % 政府支出路径
    results.I_g_path = cS.I_g_to_Y_ratio_ss * results.Y_path;
    G_c_path = cS.G_c_to_Y_ratio_ss * results.Y_path;
    
    % [核心诊断] 计算两种口径的投资
    % 方法1: 根据资本积累公式 (反映家庭储蓄行为)
    I_p_from_K_law = [K_p_path_level(2:end) - (1-cS.ddk) * K_p_path_level(1:end-1), NaN];
    g_n_period_F = (1 + cS.n_ss)^cS.time_Step - 1;
    g_A_period_F = (1 + cS.g_A_ss)^cS.time_Step - 1;
    g_total_F_period = (1 + g_A_period_F) * (1 + g_n_period_F) - 1;
    I_p_from_K_law(end) = (g_total_F_period + cS.ddk) * K_p_path_level(end);
    results.I_p_from_K_law = I_p_from_K_law;
    
    % 方法2: 根据NIPA恒等式 (作为残差，确保账户平衡)
    I_p_from_NIPA = results.Y_path - results.C_path - results.I_g_path - G_c_path;
    results.I_p_path = I_p_from_NIPA;

    % --- 步骤 6c: 调用NIPA检验函数 ---
    if ismethod('TPI', 'check_transition_NIPA')
        TPI.check_transition_NIPA(results, cS, ssF);
    else
        fprintf('   (NIPA检验函数 TPI.check_transition_NIPA 未找到，跳过检验。)\n');
    end
    
else
    fprintf('   ❌ TPI求解器未能收敛，无法进行NIPA检验。\n');
end

fprintf('\n--- 黄金标准诊断脚本NIPA检验部分执行完毕 ---\n');