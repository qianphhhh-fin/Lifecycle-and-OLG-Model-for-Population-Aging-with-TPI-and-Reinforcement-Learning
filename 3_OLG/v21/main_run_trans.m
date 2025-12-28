% =========================================================================
% == SCRIPT: main_run_transition.m
% == 版本: [v12.0 - 最终均衡求解版]
% ==
% == 核心修改:
% ==   - [!!!] TPI算法被重构，以解决会计不一致问题。
% ==   - 不再猜测 K_payg 路径，而是猜测 K_payg/Y 的比率路径 (kappa_path)。
% ==   - 这模仿了稳态求解器的成功逻辑，在每次迭代内部强制实现宏观均衡。
% ==   - 阻尼系数被调整以适应新的、耦合性更强的算法。
% =========================================================================
clear; close all;
addpath(pwd);
fprintf('=== OLG异质性模型转轨动态求解脚本 (TPI) v12.0 ===\n\n');

%% --- 1. 加载稳态数据与初始化 ---
fprintf('--- 1. 加载稳态数据与初始化 ---\n');
data_filename = 'SS/data_for_het_transition_nopps.mat'; 

if ~exist(data_filename, 'file'), error('转轨数据文件 "%s" 不存在。', data_filename); end
loaded_data_struct = load(data_filename);
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
if ~iscell(valF_h), error('加载的 valF_h 不是元胞数组，请检查 data_for_het_transition.mat 的生成脚本。'); end
nH = cS.nTypes;

max_iter = 200;
tolerance = 1e-5;
T = cS.T_sim;

fprintf('   采用保守的阻尼系数以确保收敛...\n');
damp_state = 0.3;
damp_flow = 0.5;
damping = struct(...
    'r',     damp_flow, ...
    'w',     damp_flow, ...
    'beq',   damp_flow, ...
    'tr',    damp_flow, ...
    'k_g',   damp_state, ...
    'kappa', damp_state ...  % [!!!] 新的猜测变量
    );

%% --- 2. 构造或加载猜测路径 ---
fprintf('\n--- 2. 构造或加载猜测路径 ---\n');
iter_results_filename = 'TRANS/iter_results_het_payg_nopps.mat';

guess_paths_loaded = false;
if exist(iter_results_filename, 'file')
    fprintf('   发现已保存的迭代文件: %s\n', iter_results_filename);
    saved_data = load(iter_results_filename);
    if isfield(saved_data, 'T_saved') && saved_data.T_saved == T && isfield(saved_data, 'kappa_path_iter')
        fprintf('   路径长度与变量匹配 (T=%d, 存在kappa_path)，加载已存路径。\n', T);
        r_path_guess = saved_data.r_path_iter;
        w_hat_path_guess = saved_data.w_hat_path_iter;
        bequest_gen_raw_pc_path_guess = saved_data.bequest_gen_raw_pc_path_iter;
        TR_hat_pc_path_guess = saved_data.TR_hat_pc_path_iter;
        K_g_hat_pc_path_guess = saved_data.K_g_hat_pc_path_iter;
        kappa_path_guess = saved_data.kappa_path_iter; % [!!!]
        guess_paths_loaded = true;
    else
        fprintf('   路径或变量不匹配，将重新创建猜测路径。\n');
    end
end

if ~guess_paths_loaded
    fprintf('   未加载任何文件，将从头创建线性插值的猜测路径。\n');
    r_path_guess = linspace(ss0.r_mkt, ssF.r_mkt, T);
    w_hat_path_guess = linspace(ss0.w_hat, ssF.w_hat, T);
    bequest_gen_raw_pc_path_guess = linspace(ss0.Bequest_gen_hat_raw_ss, ssF.Bequest_gen_hat_raw_ss, T);
    TR_hat_pc_path_guess = linspace(ss0.TR_distributed_agg, ssF.TR_distributed_agg, T);
    K_g_hat_pc_path_guess = linspace(ss0.K_public_hat, ssF.K_public_hat, T);
    
    % [!!!] 创建 kappa_path 的初始猜测
    kappa0 = ss0.K_payg_hat / ss0.Y_from_production_hat;
    kappaF = 0; % 终点基金为0
    kappa_path_guess = linspace(kappa0, kappaF, T);
        
    fprintf('   ✅ 已生成所有猜测路径。\n');
end

%% --- 3. [核心] 时序路径迭代 (TPI) 循环 ---
fprintf('   预计算外生劳动供给路径...\n');
L_path_h_exog = zeros(nH, T);
for h = 1:nH
    cS_h = cS; cS_h.ageEffV_new = cS.ageEffV_new_h(:, h);
    for t = 1:T
        pop_by_age_h_t = cS.Z_path_raw(:, t) * cS.type_weights(h);
        L_t_h = 0;
        for ia = 1:cS.aR_new
            L_t_h = L_t_h + pop_by_age_h_t(ia) * cS_h.ageEffV_new(ia) * paramSF.leGridV(1);
        end
        L_path_h_exog(h, t) = L_t_h;
    end
end
fprintf('   ✅ 外生劳动供给路径计算完毕。\n');

plot_config = struct('T', T, 'cS', cS, 'ss0', ss0, 'ssF', ssF, 'status', 'initialize');
iter_data_init = struct('r_path_guess', r_path_guess);
h_fig_r = TPI.plot_convergence_realtime([], iter_data_init, plot_config);

fprintf('\n--- 3. 启动TPI循环 (最大迭代: %d, 容忍度: %.1e) ---\n', max_iter, tolerance);
r_path_iter = r_path_guess;
w_hat_path_iter = w_hat_path_guess;
bequest_gen_raw_pc_path_iter = bequest_gen_raw_pc_path_guess;
TR_hat_pc_path_iter = TR_hat_pc_path_guess;
K_g_hat_pc_path_iter = K_g_hat_pc_path_guess;
kappa_path_iter = kappa_path_guess; % [!!!]

errors_format = 'Err: %.2e (r:%.1e,w:%.1e,bq:%.1e,tr:%.1e,kg:%.1e,kap:%.1e)';

total_time_start = tic;
converged = false;
for iter = 1:max_iter
    iter_time_start = tic;

    [target_paths, errors] = TPI.calculate_paths_and_errors(...
        r_path_iter, w_hat_path_iter, bequest_gen_raw_pc_path_iter, ...
        TR_hat_pc_path_iter, K_g_hat_pc_path_iter, kappa_path_iter, ... % [!!!]
        L_path_h_exog, ...
        ss0, ssF, cS, paramSF, valF_h, polF_h, dist0_h);

    iter_time_elapsed = toc(iter_time_start);
    fprintf(['Iter [%3d/%3d] | ', errors_format, ' | Time: %.2f s\n'], ...
        iter, max_iter, errors.total, errors.r, errors.w, errors.beq, errors.tr, errors.k_g, errors.kappa, iter_time_elapsed);

    if errors.total < tolerance
        fprintf('\n✅ TPI算法在第 %d 次迭代后成功收敛!\n', iter);
        converged = true;
        break;
    end

    % --- 路径更新 ---
    update_range_flow = 1:(T-1);
    r_path_iter(update_range_flow) = (1 - damping.r) * r_path_iter(update_range_flow) + damping.r * target_paths.r_path(update_range_flow);
    w_hat_path_iter(update_range_flow) = (1 - damping.w) * w_hat_path_iter(update_range_flow) + damping.w * target_paths.w_hat_path(update_range_flow);
    TR_hat_pc_path_iter(update_range_flow) = (1 - damping.tr) * TR_hat_pc_path_iter(update_range_flow) + damping.tr * target_paths.TR_hat_pc_path(update_range_flow);
    bequest_gen_raw_pc_path_iter(update_range_flow) = (1 - damping.beq) * bequest_gen_raw_pc_path_iter(update_range_flow) + damping.beq * target_paths.bequest_gen_raw_pc_path(update_range_flow);
    
    update_range_state = 2:T;
    K_g_hat_pc_path_iter(update_range_state) = (1 - damping.k_g) * K_g_hat_pc_path_iter(update_range_state) + damping.k_g * target_paths.K_g_hat_pc_path(update_range_state);
    kappa_path_iter(update_range_state) = (1 - damping.kappa) * kappa_path_iter(update_range_state) + damping.kappa * target_paths.kappa_path(update_range_state); % [!!!]

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
    'TR_hat_pc_path_iter',  'K_g_hat_pc_path_iter', 'kappa_path_iter','T_saved'); % [!!!]
fprintf('   ✅ 已将最后一次迭代的路径保存至: %s\n', iter_results_filename);

%% --- 4. 结果整理、核算与可视化 ---
if converged
    fprintf('\n--- 4. 结果整理、核算与可视化 ---\n');

    [final_paths, ~, final_Pol_path_h, final_Dist_path_h, ~] = TPI.calculate_paths_and_errors(...
        r_path_iter, w_hat_path_iter, bequest_gen_raw_pc_path_iter, ...
        TR_hat_pc_path_iter, K_g_hat_pc_path_iter, kappa_path_iter, ... % [!!!]
        L_path_h_exog, ...
        ss0, ssF, cS, paramSF, valF_h, polF_h, dist0_h);

    % --- 从 final_paths 中提取最终的、完全一致的宏观路径 ---
    results = final_paths.final_consistent_paths;
    
    % --- 保存与可视化 ---
    if ~exist('TRANS', 'dir'), mkdir('TRANS'); end
    if cS.pps_active
        results_filename = 'TRANS/TPI_results_het_pps.mat';
    else
        results_filename = 'TRANS/TPI_results_het_nopps.mat';
    end
    save(results_filename, 'results', 'cS', 'final_Pol_path_h', 'final_Dist_path_h', 'ss0', 'ssF', '-v7.3');
    fprintf('   ✅ 完整转轨结果已保存至: %s\n', results_filename);

    TPI.check_transition_NIPA(results, cS, ss0, ssF);
    TPI.plot_transition_with_demographics(results, cS, ssF, '异质性家庭情景 (含PAYG基金)');
else
    fprintf('\n--- 4. 求解失败 ---\n');
end

fprintf('\n--- 转轨动态求解脚本执行完毕 ---\n');