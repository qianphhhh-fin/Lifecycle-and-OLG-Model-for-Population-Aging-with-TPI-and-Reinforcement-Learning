% ========================================================================
% == SCRIPT: test_diffN_trans.m
% == 版本: [v3.0 - Kg内生迭代整合版]
% ==
% == 目的:
% ==   - [!!! 核心整合 !!!] 将内生公共资本(K_g)的迭代逻辑完全整合进
% ==     受控人口冲击的转轨求解框架中。
% ==   - 所有迭代变量(r, w, b_hat, TR, BQ, K_g)的单位全部统一为
% ==     【有效人均】，以确保数值稳定性和理论一致性。
% ==   - 初始猜测路径基于稳态的有效人均值构造。
% ==   - TPI循环现在同时求解价格和公共资本的均衡路径。
% ==   - 结果后处理部分将有效人均K_g路径重构为总量，用于NIPA核算。
% =========================================================================
clear; close all;
addpath(pwd);
fprintf('=== 受控人口动态转轨求解脚本 (v3.0 - Kg内生迭代整合版) ===\n\n');

%% --- 1. 加载受控实验的稳态数据 ---
fprintf('--- 1. 加载受控实验的稳态数据 ---\n');
data_filename = 'SS/data_for_diff_transition.mat';
if ~exist(data_filename, 'file'), error('受控实验数据文件 "%s" 不存在，请先运行 test_diffN_ss.m。', data_filename); end
load(data_filename, 'data_for_diff_transition');
fprintf('   ✅ 已加载受控实验数据: %s\n', data_filename);

ss0 = data_for_diff_transition.ss0;
ssF = data_for_diff_transition.ssF;
dist0 = data_for_diff_transition.dist0;
valF = data_for_diff_transition.valF;
polF = data_for_diff_transition.polF;
cS = data_for_diff_transition.cS;
paramSF = data_for_diff_transition.paramSF;

if isempty(dist0), error('加载的 dist0 为空。'); end
if ~isfield(cS, 'start_year'), cS.start_year = 2024; end

max_iter = 200;
tolerance = 1e-6;
T = cS.T_sim;
damp_f = 0.5; % 增加了K_g迭代，使用更小的阻尼因子可能更稳健
damping = struct('r', damp_f, 'w', damp_f, 'beq', damp_f, 'tr', damp_f, 'b_hat', damp_f, 'k_g', damp_f);

fprintf('   模拟期数 T = %d\n', T);
fprintf('   阻尼因子: r=%.2f, w=%.2f, bq=%.2f, tr=%.2f, b_hat=%.2f, k_g=%.2f\n', ...
    damping.r, damping.w, damping.beq, damping.tr, damping.b_hat, damping.k_g);

%% --- 2. 构造【有效人均】口径的宏观路径初始猜测 ---
fprintf('\n--- 2. 构造【有效人均】口径的宏观路径初始猜测 ---\n');

% --- 2a. 计算稳态的【有效人均】值 ---
r_pe_ss0 = ss0.r_mkt;
w_hat_pe_ss0 = ss0.w_hat;
b_hat_pe_ss0 = ss0.b_hat;
k_g_hat_pe_ss0 = ss0.K_public_hat; % ss0求解时已是人均/有效人均
total_pop_ss0 = sum(dist0(:));
bequest_gen_raw_pe_ss0 = ss0.Bequest_gen_hat_raw_ss / total_pop_ss0;
tr_pe_ss0 = ss0.TR_distributed_agg / total_pop_ss0;

r_pe_ssF = ssF.r_mkt;
w_hat_pe_ssF = ssF.w_hat;
b_hat_pe_ssF = ssF.b_hat;
k_g_hat_pe_ssF = ssF.K_public_hat; % ssF求解时已是人均/有效人均
bequest_gen_raw_pe_ssF = ssF.Bequest_gen_hat_raw_ss;
tr_pe_ssF = ssF.TR_distributed_agg;

% --- 2b. 构造平滑的【有效人均】路径 ---
r_path_guess = linspace(r_pe_ss0, r_pe_ssF, T);
w_hat_path_guess = linspace(w_hat_pe_ss0, w_hat_pe_ssF, T);
b_hat_path_guess = linspace(b_hat_pe_ss0, b_hat_pe_ssF, T);
bequest_gen_raw_pc_path_guess = linspace(bequest_gen_raw_pe_ss0, bequest_gen_raw_pe_ssF, T);
TR_hat_pc_path_guess = linspace(tr_pe_ss0, tr_pe_ssF, T);
K_g_hat_path_guess = linspace(k_g_hat_pe_ss0, k_g_hat_pe_ssF, T);

fprintf('   ✅ 已生成包含K_g在内的【有效人均】口径初始猜测路径。\n');


%% --- 3. [核心] 时序路径迭代 (TPI) 循环 ---
fprintf('\n--- 3. 启动TPI循环 (最大迭代: %d, 容忍度: %.1e) ---\n', max_iter, tolerance);
r_path_iter = r_path_guess;
w_hat_path_iter = w_hat_path_guess;
bequest_gen_raw_pc_path_iter = bequest_gen_raw_pc_path_guess;
TR_hat_pc_path_iter = TR_hat_pc_path_guess;
b_hat_path_iter = b_hat_path_guess;
K_g_hat_path_iter = K_g_hat_path_guess; % K_g加入迭代

total_time_start = tic;
converged = false;
for iter = 1:max_iter
    iter_time_start = tic;

    % [!!!] 调用包含K_g迭代的TPI核心引擎
    [target_paths, errors, ~, ~, ~] = TPI.calculate_paths_and_errors(...
        r_path_iter, w_hat_path_iter, bequest_gen_raw_pc_path_iter, ...
        TR_hat_pc_path_iter, b_hat_path_iter, K_g_hat_path_iter, ...
        ss0, ssF, cS, paramSF, valF, polF, dist0);

    iter_time_elapsed = toc(iter_time_start);
    fprintf('Iter [%3d/%3d] | Error: %.3e (r:%.2e, w:%.2e, bq:%.2e, tr:%.2e, b:%.2e, kg:%.2e) | Time: %.2f s\n', ...
        iter, max_iter, errors.total, errors.r, errors.w, errors.beq, errors.tr, errors.b_hat, errors.k_g, iter_time_elapsed);

    if errors.total < tolerance
        fprintf('\n✅ TPI算法在第 %d 次迭代后成功收敛!\n', iter);
        converged = true;
        break;
    end

    % 价格和流量变量的更新范围是 1 到 T-1
    update_range_flow = 1:(T-1);
    r_path_iter(update_range_flow) = (1 - damping.r) * r_path_iter(update_range_flow) + damping.r * target_paths.r_path(update_range_flow);
    w_hat_path_iter(update_range_flow) = (1 - damping.w) * w_hat_path_iter(update_range_flow) + damping.w * target_paths.w_hat_path(update_range_flow);
    TR_hat_pc_path_iter(update_range_flow) = (1 - damping.tr) * TR_hat_pc_path_iter(update_range_flow) + damping.tr * target_paths.TR_hat_pc_path(update_range_flow);
    b_hat_path_iter(update_range_flow) = (1 - damping.b_hat) * b_hat_path_iter(update_range_flow) + damping.b_hat * target_paths.b_hat_path(update_range_flow);
    bequest_gen_raw_pc_path_iter(update_range_flow) = (1 - damping.beq) * bequest_gen_raw_pc_path_iter(update_range_flow) + damping.beq * target_paths.bequest_gen_raw_pc_path(update_range_flow);
    
    % 状态变量 K_g 的更新范围应该是 2 到 T
    update_range_state = 2:T;
    K_g_hat_path_iter(update_range_state) = (1 - damping.k_g) * K_g_hat_path_iter(update_range_state) + damping.k_g * target_paths.K_g_hat_path(update_range_state);
    
    if iter == max_iter, warning('TPI在达到最大迭代次数后仍未收敛！'); end
end
total_time_elapsed = toc(total_time_start);
fprintf('--- TPI循环结束，总用时: %.2f 秒 ---\n', total_time_elapsed);

%% --- 4. 结果整理、核算与可视化 ---
if converged
    fprintf('\n--- 4. 结果整理、核算与可视化 ---\n');
    
    [final_paths, ~, ~, ~, final_aggr_supply] = TPI.calculate_paths_and_errors(...
        r_path_iter, w_hat_path_iter, bequest_gen_raw_pc_path_iter, ...
        TR_hat_pc_path_iter, b_hat_path_iter, K_g_hat_path_iter, ...
        ss0, ssF, cS, paramSF, valF, polF, dist0);
    
    results = struct();
    total_pop_path = sum(cS.Z_path_raw, 1);
    results.r_path = final_paths.r_path;
    results.w_path = final_paths.w_hat_path .* cS.A_path;
    
    K_p_hat_total_path = final_aggr_supply.K_p_hat_total_path;
    K_p_path_level = K_p_hat_total_path .* cS.A_path;
    results.K_p_path = K_p_path_level;
    results.C_path = final_aggr_supply.C_hat_pc_path .* total_pop_path .* cS.A_path;

    % [!!!] 从有效人均K_g路径重构总量路径用于核算
    K_g_hat_pe_path = final_paths.K_g_hat_path;
    K_g_hat_total_path = K_g_hat_pe_path .* total_pop_path;
    
    Y_hat_total_path = zeros(1, T);
    for t = 1:T
        prices_t = firm.get_prices_at_t(K_p_hat_total_path(t), K_g_hat_total_path(t), final_aggr_supply.L_path(t), cS);
        Y_hat_total_path(t) = prices_t.Y_hat_t;
    end
    results.Y_path = Y_hat_total_path .* cS.A_path;
    results.I_g_path = cS.I_g_to_Y_ratio_ss * results.Y_path;

    I_p_from_K_law = [K_p_path_level(2:end) - (1-cS.ddk) * K_p_path_level(1:end-1), NaN];
    n_period_F = sum(cS.Z_path_raw(:, T)) / sum(cS.Z_path_raw(:, T-1)) - 1;
    g_A_period_F = (1+cS.g_A_ss)^cS.time_Step - 1;
    g_total_F_period = (1 + g_A_period_F) * (1 + n_period_F) - 1;
    I_p_from_K_law(end) = (g_total_F_period + cS.ddk) * K_p_path_level(end);
    results.I_p_from_K_law = I_p_from_K_law;

    G_c_path = cS.G_c_to_Y_ratio_ss * results.Y_path;
    I_p_from_NIPA = results.Y_path - results.C_path - results.I_g_path - G_c_path;
    results.I_p_path = I_p_from_NIPA;
    
    if ismethod('TPI', 'check_transition_NIPA')
        TPI.check_transition_NIPA(results, cS, ss0, ssF);
    end

    if ismethod('TPI', 'plot_population_shock_transition')
        TPI.plot_population_shock_transition(results, cS, ssF, '受控人口冲击实验');
    else
        fprintf('   (可视化函数 TPI.plot_population_shock_transition 未找到，跳过绘图。)\n');
    end
    
    fprintf('   ✅ 结果整理完毕。\n');
else
    fprintf('\n--- 4. 求解失败 ---\n');
    fprintf('   ❌ TPI求解器未能收敛，无法生成可靠结果。\n');
end