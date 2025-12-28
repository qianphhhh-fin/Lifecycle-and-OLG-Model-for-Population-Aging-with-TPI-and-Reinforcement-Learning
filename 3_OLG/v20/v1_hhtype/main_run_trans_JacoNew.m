% =========================================================================
% == SCRIPT: main_run_transition.m
% == 版本: [v9.0 - K_g内生化与断点续算重构版]
% ==
% == 核心修改:
% ==   - 完全采纳 test_diffG_trans.m 的高级迭代框架。
% ==   - [!!!] 公共资本路径(K_g_path)现在是TPI迭代求解的内生变量之一。
% ==   - [!!!] 新增“断点续算”功能：
% ==     - 迭代开始前，检查是否存在 'TRANS/iter_results_main.mat'。
% ==     - 若存在且T匹配，则加载上次的迭代结果作为本次的初始猜测。
% ==     - 否则，重新生成基于稳态的线性猜测路径。
% ==   - [!!!] 每次迭代完成后，自动保存最新的路径变量，以备下次使用。
% ==   - 更新逻辑与变量口径与 test_diffG_trans.m 完全对齐。
% =========================================================================
clear; close all; 
addpath(pwd);
fprintf('=== OLG模型转轨动态求解脚本 (TPI)  ===\n\n');

%% --- 1. 加载稳态数据与初始化 ---
fprintf('--- 1. 加载稳态数据与初始化 ---\n');
data_filename = 'SS/data_for_transition.mat';
iter_results_filename = 'TRANS/iter_results_main.mat'; % 定义保存迭代结果的文件

if ~exist(data_filename, 'file'), error('转轨数据文件 "%s" 不存在。', data_filename); end
load(data_filename, 'data_for_transition');
fprintf('   ✅ 已加载稳态和参数: %s\n', data_filename);
fprintf('   [检测] 转轨期将以 pps_active = %s 模式运行。\n', string(data_for_transition.cS.pps_active));

ss0 = data_for_transition.ss0;
ssF = data_for_transition.ssF;
dist0 = data_for_transition.dist0;
valF = data_for_transition.valF;
polF = data_for_transition.polF;
cS = data_for_transition.cS;
paramSF = data_for_transition.paramSF;

if isempty(dist0), error('加载的 dist0 为空。'); end

max_iter = 200;
tolerance = 1e-5;
T = cS.T_sim;
damp_f = 0.5;
damping = struct('r', damp_f, 'w', damp_f, 'beq', damp_f, 'tr', damp_f, 'b_hat', damp_f, 'k_g', damp_f);

fprintf('   模拟期数 T = %d\n', T);
fprintf('   阻尼因子统一为: %.2f\n', damp_f);

%% --- 2. 构造或加载【有效人均】口径的初始猜测路径 ---
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
    % --- 计算稳态的【有效人均】值 ---
    r_pe_ss0 = ss0.r_mkt;
    w_hat_pe_ss0 = ss0.w_hat;
    b_hat_pe_ss0 = ss0.b_hat;
    k_g_hat_pe_ss0 = ss0.K_public_hat;
    total_pop_ss0 = 1; % sum(dist0(:));
    bequest_gen_raw_pe_ss0 = ss0.Bequest_gen_hat_raw_ss / total_pop_ss0;
    tr_pe_ss0 = ss0.TR_distributed_agg / total_pop_ss0;

    r_pe_ssF = ssF.r_mkt;
    w_hat_pe_ssF = ssF.w_hat;
    b_hat_pe_ssF = ssF.b_hat;
    k_g_hat_pe_ssF = ssF.K_public_hat;
    total_pop_ssF = sum(cS.Z_path_raw(:, end));
    bequest_gen_raw_pe_ssF = ssF.Bequest_gen_hat_raw_ss; 
    tr_pe_ssF = ssF.TR_distributed_agg;

    % --- 构造平滑的【有效人均】路径 ---
    r_path_guess = linspace(r_pe_ss0, r_pe_ssF, T);
    w_hat_path_guess = linspace(w_hat_pe_ss0, w_hat_pe_ssF, T);
    b_hat_path_guess = linspace(b_hat_pe_ss0, b_hat_pe_ssF, T);
    bequest_gen_raw_pc_path_guess = linspace(bequest_gen_raw_pe_ss0, bequest_gen_raw_pe_ssF, T);
    TR_hat_pc_path_guess = linspace(tr_pe_ss0, tr_pe_ssF, T);
    K_g_hat_path_guess = linspace(k_g_hat_pe_ss0, k_g_hat_pe_ssF, T);
    
    fprintf('   ✅ 已生成【有效人均】口径的线性初始猜测路径。\n');
end


%% --- 3. [核心] 高效混合求解循环 (TPI预热 + Bad Broyden法) ---
fprintf('\n--- 3. 启动高效混合求解循环 (最大迭代: %d, 容忍度: %.1e) ---\n', max_iter, tolerance);

% --- 3.1. 混合法参数设置 ---
max_iter_tpi    = 20;     % 初始TPI迭代次数 (用于预热)
lambda          = 1.0;    % Broyden法步长的阻尼因子 (通常设为1.0)
broyden_start_iter = max_iter_tpi + 1;
fprintf('   策略: %d次标准TPI预热, 然后切换到最多%d次Bad Broyden法。\n\n', ...
    max_iter_tpi, max_iter - max_iter_tpi);

% --- 3.2. 初始化迭代路径 ---
r_path_iter = r_path_guess;
w_hat_path_iter = w_hat_path_guess;
bequest_gen_raw_pc_path_iter = bequest_gen_raw_pc_path_guess;
TR_hat_pc_path_iter = TR_hat_pc_path_guess;
b_hat_path_iter = b_hat_path_guess;
K_g_hat_path_iter = K_g_hat_path_guess;

total_time_start = tic;
converged = false;
H_inv = []; % 初始化逆雅可比矩阵为空

% --- 3.3. 主循环 ---
for iter = 1:max_iter
    iter_time_start = tic;
    
    % --- 3.3.1. 将当前路径打包成一个长向量 X_k ---
    update_len_flow = T - 1;
    update_len_state = T - 2; % K_g: t=2...T-1
    
    X_k = [r_path_iter(1:update_len_flow), ...
           w_hat_path_iter(1:update_len_flow), ...
           bequest_gen_raw_pc_path_iter(1:update_len_flow), ...
           TR_hat_pc_path_iter(1:update_len_flow), ...
           b_hat_path_iter(1:update_len_flow), ...
           K_g_hat_path_iter(2:T-1)]';
           
    % --- 3.3.2. 定义变量在X_k中的位置映射 (map) ---
    len_vec = [update_len_flow, update_len_flow, update_len_flow, update_len_flow, update_len_flow, update_len_state];
    end_indices = cumsum(len_vec);
    start_indices = [1, end_indices(1:end-1)+1];
    iter_vars_map = struct(...
        'r',     start_indices(1):end_indices(1), ...
        'w_hat', start_indices(2):end_indices(2), ...
        'beq',   start_indices(3):end_indices(3), ...
        'tr',    start_indices(4):end_indices(4), ...
        'b_hat', start_indices(5):end_indices(5), ...
        'k_g',   start_indices(6):end_indices(6));
        
    % --- 3.3.3. 定义误差函数句柄 ---
    const_paths = struct('K_g_hat_path_t1', K_g_hat_path_guess(1));
    error_func_handle = @(X) TPI.tpi_error_wrapper(X, T, iter_vars_map, const_paths, ...
                                                    ss0, ssF, cS, paramSF, valF, polF, dist0);

    % --- 3.3.4. 计算当前误差向量 F_k ---
    F_k = error_func_handle(X_k);
    current_error_norm = norm(F_k, Inf);

    % --- 3.3.5. 检查收敛 ---
    if current_error_norm < tolerance
        fprintf('\n✅ 求解器在第 %d 次迭代后成功收敛!\n', iter);
        converged = true;
        break;
    end
    
    if iter < broyden_start_iter
        % ==========================
        %  阶段一: 标准TPI迭代 (预热)
        % ==========================
        [target_paths, ~] = TPI.calculate_paths_and_errors(...
            r_path_iter, w_hat_path_iter, bequest_gen_raw_pc_path_iter, ...
            TR_hat_pc_path_iter, b_hat_path_iter, K_g_hat_path_iter, ...
            ss0, ssF, cS, paramSF, valF, polF, dist0);

        iter_time_elapsed = toc(iter_time_start);
        fprintf('[TPI Iter %2d/%2d] | Error Norm (Inf): %.3e | Time: %.2f s\n', ...
            iter, max_iter_tpi, current_error_norm, iter_time_elapsed);

        % 更新路径 (使用简单的阻尼迭代)
        r_path_iter(1:update_len_flow)     = (1 - damp_f) * r_path_iter(1:update_len_flow)     + damp_f * target_paths.r_path(1:update_len_flow);
        w_hat_path_iter(1:update_len_flow) = (1 - damp_f) * w_hat_path_iter(1:update_len_flow) + damp_f * target_paths.w_hat_path(1:update_len_flow);
        TR_hat_pc_path_iter(1:update_len_flow) = (1 - damp_f) * TR_hat_pc_path_iter(1:update_len_flow) + damp_f * target_paths.TR_hat_pc_path(1:update_len_flow);
        b_hat_path_iter(1:update_len_flow) = (1 - damp_f) * b_hat_path_iter(1:update_len_flow) + damp_f * target_paths.b_hat_path(1:update_len_flow);
        bequest_gen_raw_pc_path_iter(1:update_len_flow) = (1 - damp_f) * bequest_gen_raw_pc_path_iter(1:update_len_flow) + damp_f * target_paths.bequest_gen_raw_pc_path(1:update_len_flow);
        K_g_hat_path_iter(2:T-1)             = (1 - damp_f) * K_g_hat_path_iter(2:T-1)             + damp_f * target_paths.K_g_hat_path(2:T-1);
        
        if iter == max_iter_tpi
             fprintf('\n--- TPI预热结束，计算雅可比矩阵并切换至Bad Broyden法 ---\n');
        end

    else
        % ===================================
        %  阶段二: "Bad Broyden" 法 (主求解)
        % ===================================
        if isempty(H_inv)
            % [昂贵操作] 仅在第一次进入Broyden阶段时计算一次 J 并求逆
            J_0 = TPI.calculate_jacobian_numerical(error_func_handle, X_k);
            fprintf('   [Broyden] 正在求雅可比矩阵的逆... ');
            H_inv = inv(J_0);
            fprintf('完成。\n');
        end
        
        % [廉价操作] 使用固定的 H_inv 和最新的 F_k 计算更新步长
        delta_X = -H_inv * F_k;
        
        % 更新 X_{k+1}
        X_k_new = X_k + lambda * delta_X;
        
        % 解包 X_k_new 更新迭代路径
        r_path_iter(1:update_len_flow)         = X_k_new(iter_vars_map.r)';
        w_hat_path_iter(1:update_len_flow)     = X_k_new(iter_vars_map.w_hat)';
        bequest_gen_raw_pc_path_iter(1:update_len_flow) = X_k_new(iter_vars_map.beq)';
        TR_hat_pc_path_iter(1:update_len_flow) = X_k_new(iter_vars_map.tr)';
        b_hat_path_iter(1:update_len_flow)     = X_k_new(iter_vars_map.b_hat)';
        K_g_hat_path_iter(2:T-1)               = X_k_new(iter_vars_map.k_g)';
        
        iter_time_elapsed = toc(iter_time_start);
        fprintf('[Broyden Iter %2d] | Error Norm (Inf): %.3e | Time: %.2f s\n', ...
            iter-max_iter_tpi, current_error_norm, iter_time_elapsed);
    end
    
    if iter == max_iter && ~converged, warning('求解器在达到最大迭代次数后仍未收敛！'); end
end

total_time_elapsed = toc(total_time_start);
fprintf('--- 混合求解循环结束，总用时: %.2f 秒 ---\n', total_time_elapsed);

% --- 收尾工作：确保路径的最后一期是稳态值并保存 ---
if converged
    r_path_iter(T) = ssF.r_mkt;
    w_hat_path_iter(T) = ssF.w_hat;
    total_pop_ssF = sum(cS.Z_path_raw(:, end));
    bequest_gen_raw_pc_path_iter(T) = ssF.Bequest_gen_hat_raw_ss / total_pop_ssF;
    TR_hat_pc_path_iter(T) = ssF.TR_distributed_agg / total_pop_ssF;
    b_hat_path_iter(T) = ssF.b_hat;
    K_g_hat_path_iter(T) = ssF.K_public_hat;
end

if ~exist('TRANS', 'dir'), mkdir('TRANS'); end
T_saved = T;
save(iter_results_filename, ...
    'r_path_iter', 'w_hat_path_iter', 'bequest_gen_raw_pc_path_iter', ...
    'TR_hat_pc_path_iter', 'b_hat_path_iter', 'K_g_hat_path_iter', 'T_saved');
fprintf('   ✅ 已将最后一次迭代的路径保存至: %s\n', iter_results_filename);
%% --- 4. 结果整理、核算与可视化 ---
if converged
    fprintf('\n--- 4. 结果整理、核算与可视化 ---\n');
    
    [final_paths, errors, final_Dist_path, final_Pol_path, final_aggr_supply] = TPI.calculate_paths_and_errors(...
        r_path_iter, w_hat_path_iter, bequest_gen_raw_pc_path_iter, ...
        TR_hat_pc_path_iter, b_hat_path_iter, K_g_hat_path_iter, ...
        ss0, ssF, cS, paramSF, valF, polF, dist0);
    
    results = struct();
    results.r_path = final_paths.r_path;
    results.w_path = final_paths.w_hat_path .* cS.A_path;
    
    total_pop_path = sum(cS.Z_path_raw, 1);
    K_p_hat_total_path = final_aggr_supply.K_p_hat_total_path; 
    K_p_path_level = K_p_hat_total_path .* cS.A_path;
    results.K_p_path = K_p_path_level;
    results.C_path = final_aggr_supply.C_hat_pc_path .* total_pop_path .* cS.A_path;
    
    K_g_hat_pe_path = final_paths.K_g_hat_path; 
    K_g_hat_total_path = K_g_hat_pe_path .* total_pop_path; 
    results.K_g_path = K_g_hat_total_path .* cS.A_path;
    
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
    
    % 保存最终结果
    if ~exist('TRANS', 'dir'), mkdir('TRANS'); end
    results_filename = 'TRANS/TPI_results_main_PPS.mat';
    save(results_filename, 'results', 'cS', 'final_Pol_path', 'final_Dist_path', 'ss0', 'ssF', '-v7.3');
    fprintf('   ✅ 完整转轨结果已保存至: %s\n', results_filename);
    
  % --- main_run_trans.m (末尾部分) ---

    % 可视化与检验
    TPI.check_transition_NIPA(results, cS, ss0, ssF);
    if cS.pps_active, suffix = '带PPS情景'; else, suffix = '无PPS情景'; end
    
    % [!!! 核心修改: 调用新的、功能更强的绘图函数 !!!]
    TPI.plot_transition_with_demographics(results, cS, ssF, suffix);
    
    fprintf('   ✅ 结果整理与可视化完毕。\n');
else
    fprintf('\n--- 4. 求解失败 ---\n');
    fprintf('   ❌ TPI求解器未能收敛，无法生成可靠结果。\n');
end

fprintf('\n--- 转轨动态求解脚本执行完毕 ---\n');