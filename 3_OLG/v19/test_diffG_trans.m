% ========================================================================
% == SCRIPT: test_diffG_trans.m
% == 版本: [v3.1 - 迭代路径保存与加载功能]
% ==
% == 目的:
% ==   - [!!! 新增功能 !!!] 增加了在每次运行TPI后，自动保存最后一次
% ==     迭代的路径变量。
% ==   - [!!! 新增功能 !!!] 在运行开始时，会检查是否存在已保存的迭代
% ==     路径文件。如果存在且长度匹配，则从该文件加载路径作为初始猜测，
% ==     实现“断点续算”或微调。如果不存在或长度不匹配，则重新生成
% ==     线性猜测路径。
% ==   - 核心迭代逻辑与v3.0保持一致。
% =========================================================================
clear; close all;
addpath(pwd);
fprintf('=== 受控技术动态(g_A Shock)转轨求解脚本 (v3.1) ===\n\n');

%% --- 1. 加载受控实验的稳态与路径数据 ---
fprintf('--- 1. 加载受控技术冲击实验数据 ---\n');
data_filename = 'SS/data_for_diff_G_transition.mat';
iter_results_filename = 'TRANS/iter_results_diff_G.mat'; % 定义保存迭代结果的文件

if ~exist(data_filename, 'file'), error('技术冲击实验数据文件 "%s" 不存在，请先运行 test_diffG_ss.m。', data_filename); end
load(data_filename, 'data_for_diff_G_transition');
fprintf('   ✅ 已加载技术冲击实验数据: %s\n', data_filename);

ss0 = data_for_diff_G_transition.ss0;
ssF = data_for_diff_G_transition.ssF;
dist0 = data_for_diff_G_transition.dist0;
valF = data_for_diff_G_transition.valF;
polF = data_for_diff_G_transition.polF;
cS = data_for_diff_G_transition.cS;
paramSF = data_for_diff_G_transition.paramSF;

if isempty(dist0), error('加载的 dist0 为空。'); end
if ~isfield(cS, 'start_year'), cS.start_year = 2024; end

max_iter = 200;
tolerance = 1e-5; 
T = cS.T_sim;
damp_f = 0.5; 
damping = struct('r', damp_f, 'w', damp_f, 'beq', damp_f, 'tr', damp_f, 'b_hat', damp_f, 'k_g', damp_f);

fprintf('   模拟期数 T = %d\n', T);
fprintf('   阻尼因子: r=%.2f, w=%.2f, bq=%.2f, tr=%.2f, b_hat=%.2f, k_g=%.2f\n', ...
    damping.r, damping.w, damping.beq, damping.tr, damping.b_hat, damping.k_g);

% --- 2b. 构造或加载初始猜测路径 (修正后版本) ---
fprintf('\n--- 2. 构造或加载【有效人均】口径的初始猜测路径 ---\n');

guess_paths_loaded = false;
if exist(iter_results_filename, 'file')
    fprintf('   发现已保存的迭代文件: %s\n', iter_results_filename);
    saved_data = load(iter_results_filename);
    if isfield(saved_data, 'T_saved') && saved_data.T_saved == T
        fprintf('   路径长度匹配 (T=%d)，将加载已保存的路径作为初始猜测。\n', T);
        r_path_guess = saved_data.r_path_iter;
        w_hat_path_guess = saved_data.w_hat_path_iter;
        b_hat_path_guess = saved_data.b_hat_path_iter;
        bequest_gen_raw_pc_path_guess = saved_data.bequest_gen_raw_pc_path_iter;
        TR_hat_pc_path_guess = saved_data.TR_hat_pc_path_iter;
        K_g_hat_path_guess = saved_data.K_g_hat_path_iter;
        guess_paths_loaded = true;
    else
        fprintf('   路径长度不匹配 (当前 T=%d, 已存 T=%d)，将重新创建猜测路径。\n', T, saved_data.T_saved);
    end
end

if ~guess_paths_loaded
    fprintf('   未加载任何文件，将从头创建线性插值的猜测路径。\n');
    % --- 2a. 计算稳态的【有效人均】值 (修正)---
    r_pe_ss0 = ss0.r_mkt;
    w_hat_pe_ss0 = ss0.w_hat;
    b_hat_pe_ss0 = ss0.b_hat;
    k_g_hat_pe_ss0 = ss0.K_public_hat;
    total_pop_ss0 = sum(dist0(:));
    bequest_gen_raw_pe_ss0 = ss0.Bequest_gen_hat_raw_ss / total_pop_ss0;
    tr_pe_ss0 = ss0.TR_distributed_agg / total_pop_ss0;

    r_pe_ssF = ssF.r_mkt;
    w_hat_pe_ssF = ssF.w_hat;
    b_hat_pe_ssF = ssF.b_hat;
    k_g_hat_pe_ssF = ssF.K_public_hat;
    total_pop_ssF = sum(cS.Z_path_raw(:, end)); % 使用终期总人口
    % [!!! 修正 !!!] 确保终期稳态值也转换为人均口径
    bequest_gen_raw_pe_ssF = ssF.Bequest_gen_hat_raw_ss / total_pop_ssF; 
    tr_pe_ssF = ssF.TR_distributed_agg / total_pop_ssF;

    % --- 2b. 构造平滑的【有效人均】路径 ---
    r_path_guess = linspace(r_pe_ss0, r_pe_ssF, T);
    w_hat_path_guess = linspace(w_hat_pe_ss0, w_hat_pe_ssF, T);
    b_hat_path_guess = linspace(b_hat_pe_ss0, b_hat_pe_ssF, T);
    bequest_gen_raw_pc_path_guess = linspace(bequest_gen_raw_pe_ss0, bequest_gen_raw_pe_ssF, T);
    TR_hat_pc_path_guess = linspace(tr_pe_ss0, tr_pe_ssF, T);
    K_g_hat_path_guess = linspace(k_g_hat_pe_ss0, k_g_hat_pe_ssF, T);
    
    fprintf('   ✅ 已生成【有效人均】口径的初始猜测路径。\n');
end


%% --- 3. [核心] 时序路径迭代 (TPI) 循环 ---
fprintf('\n--- 3. 启动TPI循环 (最大迭代: %d, 容忍度: %.1e) ---\n', max_iter, tolerance);
r_path_iter = r_path_guess;
w_hat_path_iter = w_hat_path_guess;
bequest_gen_raw_pc_path_iter = bequest_gen_raw_pc_path_guess;
TR_hat_pc_path_iter = TR_hat_pc_path_guess;
b_hat_path_iter = b_hat_path_guess;
K_g_hat_path_iter = K_g_hat_path_guess;

total_time_start = tic;
converged = false;
for iter = 1:max_iter
    iter_time_start = tic;

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
    
    % 状态变量 K_g 的更新范围是 2 到 T
    update_range_state = 2:T;
    K_g_hat_path_iter(update_range_state) = (1 - damping.k_g) * K_g_hat_path_iter(update_range_state) + damping.k_g * target_paths.K_g_hat_path(update_range_state);
    
    if iter == max_iter, warning('TPI在达到最大迭代次数后仍未收敛！'); end
end
total_time_elapsed = toc(total_time_start);
fprintf('--- TPI循环结束，总用时: %.2f 秒 ---\n', total_time_elapsed);

% --- [!!! 新增功能: 保存最后一次的迭代结果 !!!] ---
if ~exist('TRANS', 'dir'), mkdir('TRANS'); end
T_saved = T;
save(iter_results_filename, ...
    'r_path_iter', 'w_hat_path_iter', 'bequest_gen_raw_pc_path_iter', ...
    'TR_hat_pc_path_iter', 'b_hat_path_iter', 'K_g_hat_path_iter', 'T_saved');
fprintf('   ✅ 已将最后一次迭代的路径保存至: %s\n', iter_results_filename);


%% --- 4. 结果整理、核算与可视化 ---
if converged
    fprintf('\n--- 4. 结果整理、核算与可视化 ---\n');
    
    [final_paths, ~, ~, ~, final_aggr_supply] = TPI.calculate_paths_and_errors(...
        r_path_iter, w_hat_path_iter, bequest_gen_raw_pc_path_iter, ...
        TR_hat_pc_path_iter, b_hat_path_iter, K_g_hat_path_iter, ...
        ss0, ssF, cS, paramSF, valF, polF, dist0);
    
    results = struct();
    results.r_path = final_paths.r_path;
    results.w_path = final_paths.w_hat_path .* cS.A_path;
    
    % 从有效人均量重构总量
    total_pop_path = sum(cS.Z_path_raw, 1);
    K_p_hat_total_path = final_aggr_supply.K_p_hat_total_path; 
    K_p_path_level = K_p_hat_total_path .* cS.A_path;
    results.K_p_path = K_p_path_level;
    results.C_path = final_aggr_supply.C_hat_pc_path .* total_pop_path .* cS.A_path;
    
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
    g_A_period_F = cS.A_path(end) / cS.A_path(end-1) - 1;
    g_total_F_period = (1 + g_A_period_F) * (1 + n_period_F) - 1;
    I_p_from_K_law(end) = (g_total_F_period + cS.ddk) * K_p_path_level(end);
    results.I_p_from_K_law = I_p_from_K_law;

    G_c_path = cS.G_c_to_Y_ratio_ss * results.Y_path;
    results.I_p_path = results.Y_path - results.C_path - results.I_g_path - G_c_path;
    
    if ismethod('TPI', 'check_transition_NIPA')
        TPI.check_transition_NIPA(results, cS, ss0, ssF);
    end

    if ismethod('TPI', 'plot_transition_results')
        TPI.plot_transition_results(results, cS, ssF, '受控技术冲击(g_A)实验');
    end
    
    fprintf('   ✅ 结果整理完毕。\n');
else
    fprintf('\n--- 4. 求解失败 ---\n');
    fprintf('   ❌ TPI求解器未能收敛，无法生成可靠结果。\n');
end

fprintf('\n--- 受控技术动态(g_A Shock)转轨求解脚本执行完毕 ---\n');