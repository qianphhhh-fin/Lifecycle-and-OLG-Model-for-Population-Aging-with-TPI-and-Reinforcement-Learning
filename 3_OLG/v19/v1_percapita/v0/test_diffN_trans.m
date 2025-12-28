% 脚本: test_diffN_trans.m

% ========================================================================
% == SCRIPT: test_diffN_trans.m
% == 版本: [v2.5 - 遗赠变量名对齐最终版]
% ==
% == 目的:
% ==   - [!!!] 将所有与人均遗赠相关的迭代变量名从 `bequest_gen_hat_pc_path`
% ==     修改为 `bequest_gen_raw_pc_path`。
% ==   - 这包括初始猜测路径、TPI循环中的迭代变量、缓存文件的加载/保存，
% ==     以及对TPI核心引擎的调用。
% ==   - 此修改确保了整个脚本在变量命名上与核心函数的接口完全一致。
% =========================================================================
% clear; close all; clc; % Commented out for function conversion
addpath(pwd);
fprintf('=== 受控人口动态转轨求解脚本 (v2.5 - 遗赠变量名对齐最终版) ===\n\n');

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
tolerance = 1e-7;
T = cS.T_sim;
damp_f = 0.5;
damping = struct('r', damp_f, 'w', damp_f, 'beq', damp_f, 'tr', damp_f, 'b_hat', damp_f);

fprintf('   模拟期数 T = %d\n', T);
fprintf('   阻尼因子: r=%.2f, w=%.2f, bq=%.2f, tr=%.2f, b_hat=%.2f\n', ...
    damping.r, damping.w, damping.beq, damping.tr, damping.b_hat);

%% --- 2. 构造或加载宏观路径的初始猜测 ---
path_cache_filename = 'TRANS/converged_path_cache.mat';
guess_loaded = false;

if exist(path_cache_filename, 'file')
    fprintf('\n--- 2a. 尝试从缓存加载初始猜测路径 ---\n');
    try
        loaded_data = load(path_cache_filename);
        % [!!! 核心修正: 检查并加载正确的遗赠变量名 !!!]
        if isfield(loaded_data, 'r_path_iter') && isfield(loaded_data, 'bequest_gen_raw_pc_path_iter') && length(loaded_data.r_path_iter) == T
            r_path_guess = loaded_data.r_path_iter;
            w_hat_path_guess = loaded_data.w_hat_path_iter;
            bequest_gen_raw_pc_path_guess = loaded_data.bequest_gen_raw_pc_path_iter; 
            TR_hat_pc_path_guess = loaded_data.TR_hat_pc_path_iter;
            b_hat_path_guess = loaded_data.b_hat_path_iter;
            guess_loaded = true;
            fprintf('   ✅ 成功从 "%s" 加载初始猜测路径。\n', path_cache_filename);
        else
            fprintf('   ⚠️ 缓存文件中的路径长度或变量名与当前设置不匹配，将重新构造。\n');
        end
    catch ME
        fprintf('   ❌ 加载缓存文件时出错: %s\n', ME.message);
    end
end

if ~guess_loaded
    fprintf('\n--- 2b. 构造【人均量插值】的宏观路径初始猜测 ---\n');
    r_path_guess = linspace(ss0.r_mkt, ssF.r_mkt, T);
    w_hat_path_guess = linspace(ss0.w_hat, ssF.w_hat, T);
    b_hat_path_guess = linspace(ss0.b_hat, ssF.b_hat, T);
    
    total_pop_ss0 = sum(dist0, 'all');
    if abs(total_pop_ss0) < 1e-9, error('初始总人口为零。'); end

    % 使用未经BGP调整的原始总量，除以初始总人口，得到 t=1 的真实【人均原始遗赠产出】
    bequest_gen_raw_pc_ss0 = ss0.Bequest_gen_hat_raw_ss / total_pop_ss0;
    tr_pc_ss0 = ss0.TR_distributed_agg / total_pop_ss0;
    
    % ssF的人均量可以直接使用，因为其总人口为1
    bequest_gen_raw_pc_ssF = ssF.Bequest_gen_hat_raw_ss; 
    tr_pc_ssF = ssF.TR_distributed_agg;
    
    % [!!! 核心修正: 统一变量名 !!!]
    bequest_gen_raw_pc_path_guess = linspace(bequest_gen_raw_pc_ss0, bequest_gen_raw_pc_ssF, T);
    TR_hat_pc_path_guess = linspace(tr_pc_ss0, tr_pc_ssF, T);
    fprintf('   ✅ 已通过正确的人均量计算逻辑，生成了会计一致的初始猜测路径。\n');
end

%% --- 3. [核心] 时序路径迭代 (TPI) 循环 ---
fprintf('\n--- 3. 启动TPI循环 (最大迭代: %d, 容忍度: %.1e) ---\n', max_iter, tolerance);
r_path_iter = r_path_guess;
w_hat_path_iter = w_hat_path_guess;
% [!!! 核心修正: 统一变量名 !!!]
bequest_gen_raw_pc_path_iter = bequest_gen_raw_pc_path_guess;
TR_hat_pc_path_iter = TR_hat_pc_path_guess;
b_hat_path_iter = b_hat_path_guess;

total_time_start = tic;
converged = false;
for iter = 1:max_iter
    iter_time_start = tic;

    % [!!! 核心修正: 调用时使用正确的变量名 !!!]
    [target_paths, errors, ~, ~, ~] = TPI.calculate_paths_and_errors(...
        r_path_iter, w_hat_path_iter, bequest_gen_raw_pc_path_iter, ...
        TR_hat_pc_path_iter, b_hat_path_iter, ...
        ss0, ssF, cS, paramSF, valF, polF, dist0);

    iter_time_elapsed = toc(iter_time_start);
    fprintf('Iter [%3d/%3d] | Error: %.3e (r:%.2e, w:%.2e, bq:%.2e, tr:%.2e, b:%.2e) | Time: %.2f s\n', ...
        iter, max_iter, errors.total, errors.r, errors.w, errors.beq, errors.tr, errors.b_hat, iter_time_elapsed);

    if errors.total < tolerance
        fprintf('\n✅ TPI算法在第 %d 次迭代后成功收敛!\n', iter);
        converged = true;
        break;
    end

    update_range = 1:(T-2);
    r_path_iter(update_range) = (1 - damping.r) * r_path_iter(update_range) + damping.r * target_paths.r_path(update_range);
    w_hat_path_iter(update_range) = (1 - damping.w) * w_hat_path_iter(update_range) + damping.w * target_paths.w_hat_path(update_range);
    TR_hat_pc_path_iter(update_range) = (1 - damping.tr) * TR_hat_pc_path_iter(update_range) + damping.tr * target_paths.TR_hat_pc_path(update_range);
    b_hat_path_iter(update_range) = (1 - damping.b_hat) * b_hat_path_iter(update_range) + damping.b_hat * target_paths.b_hat_path(update_range);
    
    % [!!! 核心修正: 更新时使用正确的变量名 !!!]
    bequest_gen_raw_pc_path_iter(update_range) = (1 - damping.beq) * bequest_gen_raw_pc_path_iter(update_range) + damping.beq * target_paths.bequest_gen_raw_pc_path(update_range);

    if iter == max_iter, warning('TPI在达到最大迭代次数后仍未收敛！'); end
end
total_time_elapsed = toc(total_time_start);
fprintf('--- TPI循环结束，总用时: %.2f 秒 ---\n', total_time_elapsed);

%% --- 4. 结果整理、保存与可视化 ---
if converged
    fprintf('\n--- 4. 结果整理、保存与可视化 ---\n');
    
    fprintf('   正在保存收敛路径至缓存文件...\n');
    try
        if ~exist('TRANS', 'dir'), mkdir('TRANS'); end
        save(path_cache_filename, 'r_path_iter', 'w_hat_path_iter', 'bequest_gen_raw_pc_path_iter', 'TR_hat_pc_path_iter', 'b_hat_path_iter');
        fprintf('   ✅ 路径已成功保存至: %s\n', path_cache_filename);
    catch ME
        warning('保存路径缓存文件时出错: %s', ME.message);
    end

    [final_paths, ~, ~, ~, final_aggr_supply] = TPI.calculate_paths_and_errors(...
        r_path_iter, w_hat_path_iter, bequest_gen_raw_pc_path_iter, ...
        TR_hat_pc_path_iter, b_hat_path_iter, ...
        ss0, ssF, cS, paramSF, valF, polF, dist0);
    
    results = struct();
    total_pop_path = sum(cS.Z_path_raw, 1);
    results.r_path = final_paths.r_path;
    results.w_path = final_paths.w_hat_path .* cS.A_path;
    
    K_p_path_level = final_aggr_supply.K_p_hat_total_path .* cS.A_path;
    results.K_p_path = K_p_path_level;
    results.C_path = final_aggr_supply.C_hat_pc_path .* total_pop_path .* cS.A_path;

    K_g_hat_total_path = (ssF.K_public_hat / sum(cS.Z_path(:,end))) * total_pop_path;
    Y_hat_total_path = zeros(1, T);
    for t = 1:T
        prices_t = firm.get_prices_at_t(final_aggr_supply.K_p_hat_total_path(t), K_g_hat_total_path(t), final_aggr_supply.L_path(t), cS);
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

    % [!!! 核心修改: 调用新的专用绘图函数 !!!]
    if ismethod('TPI', 'plot_population_shock_transition')
        cS = rmfield( cS ,'start_year');
        TPI.plot_population_shock_transition(results, cS, ssF, '受控人口冲击实验');
    else
        fprintf('   (可视化函数 TPI.plot_population_shock_transition 未找到，跳过绘图。)\n');
    end
    
    fprintf('   ✅ 结果整理完毕。\n');
else
    fprintf('\n--- 4. 求解失败 ---\n');
    fprintf('   ❌ TPI求解器未能收敛，无法生成可靠结果。\n');
end

