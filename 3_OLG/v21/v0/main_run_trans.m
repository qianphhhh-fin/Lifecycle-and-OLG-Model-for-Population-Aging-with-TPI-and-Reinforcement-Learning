% =========================================================================
% == SCRIPT: main_run_transition.m
% == 版本: [v9.4 - 滞后福利公式+保守阻尼]
% ==
% == 核心修改:
% ==   - 采用极其保守的阻尼系数来应对增强的系统耦合性。
% ==   - TPI引擎的调用接口保持不变，因为所有滞后逻辑都在引擎内部处理。
% =========================================================================
clear; close all;
addpath(pwd);
fprintf('=== OLG异质性模型转轨动态求解脚本 (TPI) v9.4 ===\n\n');

%% --- 1. 加载稳态数据与初始化 ---
fprintf('--- 1. 加载稳态数据与初始化 ---\n');
data_filename = 'SS/data_for_het_transition_nopps.mat';

if ~exist(data_filename, 'file'), error('转轨数据文件 "%s" 不存在。', data_filename); end
loaded_data_struct = load(data_filename, 'data_for_transition');
data = loaded_data_struct.data_for_transition;

fprintf('   ✅ 已加载异质性稳态和参数: %s\n', data_filename);

ss0 = data.ss0;
ssF = data.ssF;
dist0_h = data.dist0_h;
valF_h = data.valF_h;
polF_h = data.polF_h;
cS = data.cS;
paramSF = data.paramSF;

if isempty(dist0_h), error('加载的 dist0_h 为空。'); end
if ~iscell(valF_h), error('加载的 valF_h 不是元胞数组 (cell array)，请检查 data_for_het_transition.mat 的生成脚本。'); end
nH = cS.nTypes;

max_iter = 300; % 增加迭代次数以适应更慢的收敛
tolerance = 1e-5;
T = cS.T_sim;

% [核心修改] 使用极度保守的阻尼系数
fprintf('   采用更为保守的阻尼系数以提高收敛精度...\n');
damp_very_slow = 0.01; % 极度保守
damp_slow = 0.1;     % 保守
damp_fast = 0.1;      % 相对快速

damping = struct(...
    'adj',   damp_very_slow, ... % 新的强耦合变量，必须极度保守
    'w',     damp_slow,      ... % 工资直接影响adj，也需保守
    'r',     damp_slow,      ... % 利率间接受影响
    'beq',   damp_fast,      ... 
    'tr',    damp_fast,      ... 
    'k_g',   damp_fast       ...  
    );
%% --- 2. 构造或加载猜测路径 ---
fprintf('\n--- 2. 构造或加载猜测路径 ---\n');

if cS.pps_active
    iter_results_filename = 'TRANS/iter_results_het_pps.mat';
else
    iter_results_filename = 'TRANS/iter_results_het_nopps.mat';
end

guess_paths_loaded = false;
if exist(iter_results_filename, 'file')
    fprintf('   发现已保存的迭代文件: %s\n', iter_results_filename);
    saved_data = load(iter_results_filename);
    if isfield(saved_data, 'T_saved') && saved_data.T_saved == T && isfield(saved_data, 'adj_path_iter')
        fprintf('   路径长度与变量匹配 (T=%d, 存在adj_path)，加载已存路径。\n', T);
        r_path_guess = saved_data.r_path_iter;
        w_hat_path_guess = saved_data.w_hat_path_iter;
        adj_path_guess = saved_data.adj_path_iter;
        bequest_gen_raw_pc_path_guess = saved_data.bequest_gen_raw_pc_path_iter; 
        TR_hat_pc_path_guess = saved_data.TR_hat_pc_path_iter;
        K_g_hat_path_guess = saved_data.K_g_hat_path_iter;
        guess_paths_loaded = true;
    else
        fprintf('   路径或变量不匹配，将重新创建猜测路径。\n');
    end
end

if ~guess_paths_loaded
    fprintf('   未加载任何文件，将从头创建线性插值的猜测路径。\n');
    r_path_guess = linspace(ss0.r_mkt, ssF.r_mkt, T);
    w_hat_path_guess = linspace(ss0.w_hat, ssF.w_hat, T);

    TR_hat_pc_path_guess = linspace(ss0.TR_distributed_agg , ssF.TR_distributed_agg, T);
    K_g_hat_path_guess = linspace(ss0.K_public_hat, ssF.K_public_hat, T);

    bequest_gen_raw_pc_path_guess = linspace(ss0.Bequest_gen_hat_raw_ss, ssF.Bequest_gen_hat_raw_ss, T);
    
    adj_path_guess = linspace(ss0.adj, ssF.adj, T);
    adj_path_guess(10:end) = adj_path_guess(end) ;
    
    fprintf('   ✅ 已生成基于调整因子(adj)的线性初始猜测路径。\n');
end

%% --- 3. [核心] 时序路径迭代 (TPI) 循环 ---
plot_config = struct('T', T, 'cS', cS, 'ss0', ss0, 'ssF', ssF, 'status', 'initialize');
iter_data_init = struct('r_path_guess', r_path_guess);
h_fig_r = TPI.plot_convergence_realtime([], iter_data_init, plot_config);

fprintf('\n--- 3. 启动TPI循环 (最大迭代: %d, 容忍度: %.1e) ---\n', max_iter, tolerance);
r_path_iter = r_path_guess;
w_hat_path_iter = w_hat_path_guess;
adj_path_iter = adj_path_guess; 
bequest_gen_raw_pc_path_iter = bequest_gen_raw_pc_path_guess;
TR_hat_pc_path_iter = TR_hat_pc_path_guess;
K_g_hat_path_iter = K_g_hat_path_guess;

total_time_start = tic;
converged = false;
for iter = 1:max_iter
    iter_time_start = tic;

    [target_paths, errors] = TPI.calculate_paths_and_errors(...
        r_path_iter, w_hat_path_iter, bequest_gen_raw_pc_path_iter, ...
        TR_hat_pc_path_iter, adj_path_iter, K_g_hat_path_iter, ...
        ss0, ssF, cS, paramSF, valF_h, polF_h, dist0_h);

    iter_time_elapsed = toc(iter_time_start);
    fprintf('Iter [%3d/%3d] | Error: %.3e (r:%.2e, w:%.2e, bq:%.2e, tr:%.2e, adj:%.2e, kg:%.2e) | Time: %.2f s\n', ...
        iter, max_iter, errors.total, errors.r, errors.w, errors.beq, errors.tr, errors.adj, errors.k_g, iter_time_elapsed);

    if errors.total < tolerance
        fprintf('\n✅ TPI算法在第 %d 次迭代后成功收敛!\n', iter);
        converged = true;
        break;
    end

    update_range_flow = 1:(T-1);
    r_path_iter(update_range_flow) = (1 - damping.r) * r_path_iter(update_range_flow) + damping.r * target_paths.r_path(update_range_flow);
    w_hat_path_iter(update_range_flow) = (1 - damping.w) * w_hat_path_iter(update_range_flow) + damping.w * target_paths.w_hat_path(update_range_flow);
    TR_hat_pc_path_iter(update_range_flow) = (1 - damping.tr) * TR_hat_pc_path_iter(update_range_flow) + damping.tr * target_paths.TR_hat_pc_path(update_range_flow);
    
    adj_path_iter(update_range_flow) = (1 - damping.adj) * adj_path_iter(update_range_flow) + damping.adj * target_paths.adj_path(update_range_flow);

    bequest_gen_raw_pc_path_iter(update_range_flow) = (1 - damping.beq) * bequest_gen_raw_pc_path_iter(update_range_flow) + damping.beq * target_paths.bequest_gen_raw_pc_path(update_range_flow);

    update_range_state = 2:T;
    K_g_hat_path_iter(update_range_state) = (1 - damping.k_g) * K_g_hat_path_iter(update_range_state) + damping.k_g * target_paths.K_g_hat_path(update_range_state);

    if iter == max_iter, warning('TPI在达到最大迭代次数后仍未收敛！'); end

    if exist('h_fig_r', 'var') && ishandle(h_fig_r)
        plot_config.status = 'update';
        iter_data_update = struct('iter', iter, 'r_path_iter', r_path_iter, 'error_total', errors.total);
        TPI.plot_convergence_realtime(h_fig_r, iter_data_update, plot_config);
    end    
end
total_time_elapsed = toc(total_time_start);

if exist('h_fig_r', 'var') && ishandle(h_fig_r)
    plot_config.status = 'finalize';
    iter_data_final = struct('iter', iter, 'r_path_iter', r_path_iter);
    TPI.plot_convergence_realtime(h_fig_r, iter_data_final, plot_config);
end

if ~exist('TRANS', 'dir'), mkdir('TRANS'); end
T_saved = T;
save(iter_results_filename, ...
    'r_path_iter', 'w_hat_path_iter', 'bequest_gen_raw_pc_path_iter', ...
    'TR_hat_pc_path_iter', 'adj_path_iter', 'K_g_hat_path_iter', 'T_saved');
fprintf('   ✅ 已将最后一次迭代的路径保存至: %s\n', iter_results_filename);

%% --- 4. 结果整理、核算与可视化 ---
if converged
    fprintf('\n--- 4. 结果整理、核算与可视化 ---\n');

    [final_paths, ~, final_Pol_path_h, final_Dist_path_h, final_aggr_supply] = TPI.calculate_paths_and_errors(...
        r_path_iter, w_hat_path_iter, bequest_gen_raw_pc_path_iter, ...
        TR_hat_pc_path_iter, adj_path_iter, K_g_hat_path_iter, ...
        ss0, ssF, cS, paramSF, valF_h, polF_h, dist0_h);

    results = struct();
    results.r_path = final_paths.r_path;
    results.w_path = final_paths.w_hat_path .* cS.A_path;
    results.adj_path = final_paths.adj_path; 

    K_p_hat_total_path = final_aggr_supply.K_p_hat_total_path;
    results.K_p_path = K_p_hat_total_path .* cS.A_path;
    results.C_path = final_aggr_supply.C_hat_total_path .* cS.A_path;

    total_pop_path = sum(cS.Z_path_raw, 1);
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

    I_p_from_K_law = [results.K_p_path(2:end) - (1-cS.ddk) * results.K_p_path(1:end-1), NaN];
    n_period_F = sum(cS.Z_path_raw(:, T)) / sum(cS.Z_path_raw(:, T-1)) - 1;
    g_A_period_F = cS.A_path(end) / cS.A_path(end-1) - 1;
    g_total_F_period = (1 + g_A_period_F) * (1 + n_period_F) - 1;
    I_p_from_K_law(end) = (g_total_F_period + cS.ddk) * results.K_p_path(end);
    results.I_p_from_K_law = I_p_from_K_law;

    G_c_path = cS.G_c_to_Y_ratio_ss * results.Y_path;
    results.I_p_path = results.Y_path - results.C_path - results.I_g_path - G_c_path;

    if ~exist('TRANS', 'dir'), mkdir('TRANS'); end
    if cS.pps_active
        results_filename = 'TRANS/TPI_results_het_pps.mat';
    else
        results_filename = 'TRANS/TPI_results_het_nopps.mat';
    end
    save(results_filename, 'results', 'cS', 'final_Pol_path_h', 'final_Dist_path_h', 'ss0', 'ssF', '-v7.3');
    fprintf('   ✅ 完整转轨结果已保存至: %s\n', results_filename);

    TPI.check_transition_NIPA(results, cS, ss0, ssF);
    TPI.plot_transition_with_demographics(results, cS, ssF, '异质性家庭情景');
else
    fprintf('\n--- 4. 求解失败 ---\n');
end

fprintf('\n--- 转轨动态求解脚本执行完毕 ---\n');