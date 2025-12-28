% =========================================================================
% == SCRIPT: main_run_transition.m
% == 版本: [v8.1 - 端点锚定稳定版]
% ==
% == 核心修改:
% ==   - 采纳建议，在TPI迭代中严格固定路径的期末值(t=T)为ssF的对应值。
% ==   - 路径更新的范围从(2:T)调整为(2:T-1)，确保终点值不被更新。
% ==   - 这将显著提升算法的稳定性和收敛速度。
% =========================================================================
clear; close all; clc;
addpath(pwd);
fprintf('=== OLG模型转轨动态求解脚本 (TPI) - v8.1 端点锚定稳定版 ===\n\n');

%% --- 1. 加载稳态数据与初始化 ---
fprintf('--- 1. 加载稳态数据与初始化 ---\n');
data_filename = 'SS/data_for_transition.mat';
if ~exist(data_filename, 'file'), error('转轨数据文件 "%s" 不存在。', data_filename); end
load(data_filename, 'data_for_transition');
fprintf('   ✅ 已加载稳态和参数: %s\n', data_filename);
fprintf('   [检测] 转轨期将以 pps_active = %s 模式运行。\n', string(data_for_transition.cS.pps_active));

ss0 = data_for_transition.ss0;
ssF = data_for_transition.ssF;
dist0 = data_for_transition.dist0; % dist0 已经是与 Z_path_raw(:,1) 对齐的绝对分布
valF = data_for_transition.valF;
polF = data_for_transition.polF;
cS = data_for_transition.cS;
paramSF = data_for_transition.paramSF;

if isempty(dist0), error('加载的 dist0 为空。'); end

max_iter = 200;
tolerance = 1e-5;

% [稳定化策略] 保守的阻尼因子
damping.r = 0.1;
damping.w = 0.1;
damping.beq = 0.1;
damping.tr = 0.1;
damping.b_hat = 0.1; % 为新增的 b_hat 路径设置阻尼

T = cS.T_sim;
fprintf('   模拟期数 T = %d\n', T);
fprintf('   [稳定化策略] 使用保守的阻尼因子: r=%.2f, w=%.2f, bq=%.2f, tr=%.2f, b_hat=%.2f\n', ...
    damping.r, damping.w, damping.beq, damping.tr, damping.b_hat);

%% --- 2. 构造宏观路径的初始猜测 ---
fprintf('\n--- 2. 构造宏观路径的初始猜测 ---\n');
r_path_guess = linspace(ss0.r_mkt, ssF.r_mkt, T);
w_hat_path_guess = linspace(ss0.w_hat, ssF.w_hat, T);
b_hat_path_guess = linspace(ss0.b_hat, ssF.b_hat, T); % 新增b_hat路径猜测

mass_newborns0_norm = cS.Z_path(1,1);
mass_newbornsF_norm = cS.Z_path(1,end);

beq0_per_newborn_hat = ss0.Bequest_generated_agg / mass_newborns0_norm;
beqF_per_newborn_hat = ssF.Bequest_generated_agg / mass_newbornsF_norm;
Bequest_per_newborn_path_guess = linspace(beq0_per_newborn_hat, beqF_per_newborn_hat, T);

tr0_per_hh_hat = ss0.TR_distributed_agg;
trF_per_hh_hat = ssF.TR_distributed_agg;
TR_per_hh_path_guess = linspace(tr0_per_hh_hat, trF_per_hh_hat, T);

fprintf('   ✅ 所有价格和人均路径的初始猜测已生成 (已修正为与归一化人口兼容)。\n');

%% --- 3. [核心] 时序路径迭代 (TPI) 循环 ---
fprintf('\n--- 3. 启动TPI循环 (最大迭代: %d, 容忍度: %.1e) ---\n', max_iter, tolerance);
r_path_iter = r_path_guess;
w_hat_path_iter = w_hat_path_guess;
Bequest_per_newborn_path_iter = Bequest_per_newborn_path_guess;
TR_per_hh_path_iter = TR_per_hh_path_guess;
b_hat_path_iter = b_hat_path_guess; % 初始化b_hat迭代路径

total_time_start = tic;
for iter = 1:max_iter
    iter_time_start = tic;

    [target_paths, errors] = TPI.calculate_paths_and_errors(...
        r_path_iter, w_hat_path_iter, Bequest_per_newborn_path_iter, ...
        TR_per_hh_path_iter, b_hat_path_iter, ...
        ss0, ssF, cS, paramSF, valF, polF, dist0);

    iter_time_elapsed = toc(iter_time_start);
    fprintf('Iter [%3d/%3d] | Error: %.3e (r:%.2e, w:%.2e, bq:%.2e, tr:%.2e, b:%.2e) | Time: %.2f s\n', ...
        iter, max_iter, errors.total, errors.r, errors.w, errors.beq, errors.tr, errors.b_hat, iter_time_elapsed);

    if errors.total < tolerance, fprintf('\n✅ TPI算法在第 %d 次迭代后成功收敛!\n', iter); break; end

    % [!!! 核心修正: 更新范围调整 !!!]
    % 对 t=2...T-1 的中间路径进行阻尼更新
    update_range = 2:(T-1);
    r_path_iter(update_range) = (1 - damping.r) * r_path_iter(update_range) + damping.r * target_paths.r_path(update_range);
    w_hat_path_iter(update_range) = (1 - damping.w) * w_hat_path_iter(update_range) + damping.w * target_paths.w_hat_path(update_range);
    Bequest_per_newborn_path_iter(update_range) = (1 - damping.beq) * Bequest_per_newborn_path_iter(update_range) + damping.beq * target_paths.Bequest_per_newborn_path(update_range);
    TR_per_hh_path_iter(update_range) = (1 - damping.tr) * TR_per_hh_path_iter(update_range) + damping.tr * target_paths.TR_per_hh_path(update_range);
    b_hat_path_iter(update_range) = (1 - damping.b_hat) * b_hat_path_iter(update_range) + damping.b_hat * target_paths.b_hat_path(update_range);

    % [稳定化策略] 严格锚定 t=1 的初始条件
    r_path_iter(1) = ss0.r_mkt;
    w_hat_path_iter(1) = ss0.w_hat;
    Bequest_per_newborn_path_iter(1) = beq0_per_newborn_hat;
    TR_per_hh_path_iter(1) = tr0_per_hh_hat;
    b_hat_path_iter(1) = ss0.b_hat;
    
    % [!!! 核心修正: 新增终点锚定 !!!]
    % 严格锚定 t=T 的终点条件
    r_path_iter(T) = ssF.r_mkt;
    w_hat_path_iter(T) = ssF.w_hat;
    Bequest_per_newborn_path_iter(T) = beqF_per_newborn_hat;
    TR_per_hh_path_iter(T) = trF_per_hh_hat;
    b_hat_path_iter(T) = ssF.b_hat;


    if iter == max_iter, warning('TPI在达到最大迭代次数后仍未收敛！'); end
end
total_time_elapsed = toc(total_time_start);
fprintf('--- TPI求解完成，总耗时: %.2f 秒 ---\n', total_time_elapsed);

%% --- 4 & 5. 保存、可视化、检验 ---
if errors.total < tolerance
    fprintf('\n--- 4. 生成并保存最终转轨路径结果 ---\n');
    if ~exist('TRANS', 'dir'), mkdir('TRANS'); end
    results_filename = 'TRANS/TPI_results.mat';
    
    results = struct();
    results.r_path = r_path_iter;
    results.w_hat_path = w_hat_path_iter;
    results.Bequest_per_newborn_path = Bequest_per_newborn_path_iter;
    results.TR_per_hh_path = TR_per_hh_path_iter;
    results.b_hat_path = b_hat_path_iter;
    
    final_pathS = TPI.fill_auxiliary_paths(results, cS);
    [final_Pol_path, ~] = household.backward_hh(final_pathS, cS, paramSF, valF, polF);
    final_Dist_path = distribution.simulate_dist_forward(dist0, final_Pol_path, cS, paramSF, results.Bequest_per_newborn_path);
    final_aggr = aggregates.get_path_aggregates(final_Dist_path, final_Pol_path, cS, ss0, paramSF);
    
    % 从聚合结果中提取路径
    results.L_path = final_aggr.L_path;
    results.K_p_hat_path = final_aggr.K_p_hat_path;
    results.C_hat_path = final_aggr.C_hat_path;
    
    % 计算水平量 (level) 路径
    results.w_path = results.w_hat_path .* cS.A_path;
    results.K_p_path = results.K_p_hat_path .* cS.A_path;
    results.C_path = results.C_hat_path .* cS.A_path;
    
    % 计算总量路径
    total_pop_path = sum(cS.Z_path_raw, 1);
    results.Bequest_path_total = results.Bequest_per_newborn_path .* cS.Z_path_raw(1,:);
    results.TR_path_total = results.TR_per_hh_path .* total_pop_path .* cS.A_path;
    
    % 计算其余宏观量
    results.K_g_path = ssF.K_public_hat .* (total_pop_path / total_pop_path(1)) .* cS.A_path;
    Y_hat_path_final = zeros(1, T);
    for t = 1:T
        prices_t = firm.get_prices_at_t(results.K_p_hat_path(t), ssF.K_public_hat * (total_pop_path(t)/total_pop_path(1)), results.L_path(t), cS);
        Y_hat_path_final(t) = prices_t.Y_hat_t;
    end
    results.Y_path = Y_hat_path_final .* cS.A_path;

    K_p_end_path_level = [results.K_p_path(2:T), ssF.K_private_hat * cS.A_path(end) * (1+cS.g_A_ss)^cS.time_Step * (1+cS.n_ss)^cS.time_Step];
    results.I_p_path = K_p_end_path_level - (1 - cS.ddk) .* results.K_p_path;
    K_g_end_path_level = [results.K_g_path(2:T), ssF.K_public_hat * cS.A_path(end) * (1+cS.g_A_ss)^cS.time_Step * (1+cS.n_ss)^cS.time_Step];
    results.I_g_path = K_g_end_path_level - (1 - cS.ddk_g) .* results.K_g_path;
    
    save(results_filename, 'results', 'cS', 'final_Pol_path', 'final_Dist_path', 'ss0', 'ssF', '-v7.3');
    fprintf('   ✅ 完整转轨结果已保存至: %s\n', results_filename);
    
    fprintf('\n--- 5. 可视化与检验 ---\n');
    if cS.pps_active, suffix = '带PPS情景'; else, suffix = '无PPS情景'; end
    TPI.plot_transition_results(results, cS, ssF, suffix);
    TPI.check_transition_NIPA(results, cS, ssF);
else
    fprintf('\n--- TPI未收敛，跳过最终结果处理 ---\n');
end

fprintf('\n--- 转轨动态求解脚本执行完毕 ---\n');