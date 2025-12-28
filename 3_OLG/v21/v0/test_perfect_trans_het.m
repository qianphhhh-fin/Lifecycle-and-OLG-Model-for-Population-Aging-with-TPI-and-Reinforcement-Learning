% =========================================================================
% == SCRIPT: test_perfect_trans_het.m
% == 版本: [v1.0 - 异质性黄金标准转轨测试]
% ==
% == 目的:
% ==   - 加载由 test_perfect_ss_het.m 生成的两个完美稳态。
% ==   - 在这两个稳态之间运行TPI，检验算法在理想环境下的收敛性。
% ==   - 使用之前确定的超保守阻尼策略。
% =========================================================================
clear; close all;
addpath(pwd);
fprintf('=== 异质性模型理想化环境转轨测试脚本 (v1.0) ===\n\n');

%% --- 1. 加载完美的稳态数据 ---
fprintf('--- 1. 加载理想化稳态数据与初始化 ---\n');
data_filename = 'SS/data_for_perfect_het_transition.mat';

if ~exist(data_filename, 'file'), error('完美转轨数据文件 "%s" 不存在。', data_filename); end
loaded_data_struct = load(data_filename, 'data_for_perfect_transition');
data = loaded_data_struct.data_for_perfect_transition;

fprintf('   ✅ 已加载理想化异质性稳态和参数: %s\n', data_filename);

ss0 = data.ss0;
ssF = data.ssF;
dist0_h = data.dist0_h;
valF_h = data.valF_h;
polF_h = data.polF_h;
cS = data.cS;
paramSF = data.paramSF;

nH = cS.nTypes;
max_iter = 500; 
tolerance = 1e-5;
T = cS.T_sim;

% 使用之前确定的超保守阻尼系数
fprintf('   采用【极其保守】的阻尼系数以应对强耦合系统...\n');
damp_ultra_slow = 0.02; 
damp_very_slow = 0.1;   
damp_slow = 0.25;       
damp_fast = 0.5;        

damping = struct(...
    'adj',   damp_ultra_slow, ... 
    'w',     damp_very_slow,  ... 
    'r',     damp_slow,       ... 
    'beq',   damp_fast,       ... 
    'tr',    damp_fast,       ... 
    'k_g',   damp_fast        ...  
    );

%% --- 2. 构造初始猜测路径 (在完美稳态间线性插值) ---
fprintf('\n--- 2. 构造线性猜测路径 ---\n');
r_path_guess = linspace(ss0.r_mkt, ssF.r_mkt, T);
w_hat_path_guess = linspace(ss0.w_hat, ssF.w_hat, T);
adj_path_guess = linspace(ss0.adj, ssF.adj, T);

% [注意] 这里的dist0_h是绝对分布，而ss0的量是基于归一化人口算的
% 为了构造一致的人均量，我们需要使用理论人口总质量
total_pop_ss0 = sum(sum(cS.Z_path_raw(:,1)) .* cS.type_weights');
total_pop_ssF = sum(sum(cS.Z_path_raw(:,end)) .* cS.type_weights');

tr_pc_ss0 = ss0.TR_distributed_agg;
tr_pc_ssF = ssF.TR_distributed_agg;
TR_hat_pc_path_guess = linspace(tr_pc_ss0, tr_pc_ssF, T);

kg_pc_ss0 = ss0.K_public_hat;
kg_pc_ssF = ssF.K_public_hat;
K_g_hat_path_guess = linspace(kg_pc_ss0, kg_pc_ssF, T);

bq_pc_ss0 = ss0.Bequest_gen_hat_raw_ss;
bq_pc_ssF = ssF.Bequest_gen_hat_raw_ss;
bequest_gen_raw_pc_path_guess = linspace(bq_pc_ss0, bq_pc_ssF, T);

fprintf('   ✅ 已生成基于正确【人均】单位的线性初始猜测路径。\n');


%% --- 3. [核心] 运行TPI循环 ---
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

%% --- 4. 结果分析 ---
if converged
    fprintf('\n--- 4. 结果分析与可视化 ---\n');

    [final_paths, ~, ~, ~, ~] = TPI.calculate_paths_and_errors(...
        r_path_iter, w_hat_path_iter, bequest_gen_raw_pc_path_iter, ...
        TR_hat_pc_path_iter, adj_path_iter, K_g_hat_path_iter, ...
        ss0, ssF, cS, paramSF, valF_h, polF_h, dist0_h);

    % (此处可以添加更详细的国民账户核算和绘图，但当前主要关注收敛性)
    fprintf('✅ 模型在理想化环境下成功收敛，证明TPI算法在新框架下是可行的。\n');
    TPI.plot_convergence_realtime(h_fig_r, struct('iter', iter, 'r_path_iter', r_path_iter), struct('status', 'finalize'));

else
    fprintf('\n--- 4. 求解失败 ---\n');
    fprintf('⚠️ 模型在理想化环境下依然无法收敛。这表明问题根植于模型的核心反馈机制或TPI算法本身。\n');
end

fprintf('\n--- 理想化环境转轨测试脚本执行完毕 ---\n');