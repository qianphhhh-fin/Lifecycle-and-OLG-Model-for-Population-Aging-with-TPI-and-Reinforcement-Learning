% =========================================================================
% == SCRIPT: main_run_transition.m
% == 版本: [v7.0 - 数值稳定最终版]
% ==
% == 核心修改:
% ==   - [稳定化策略1] 采用更保守、更激进的阻尼因子，特别是针对
% ==     反馈回路最强的工资(w)和其衍生的遗赠(beq)与转移(tr)路径。
% ==   - [稳定化策略2] 严格锚定初始条件。在每次迭代的阻尼更新之后，
% ==     强制将所有路径在 t=1 时的值重置为初始稳态ss0的对应值。
% ==     这可以防止来自迭代的误差"污染"固定的初始状态，是稳定TPI的
% ==     关键技巧。
% ==   - [代码清晰化] 更新了注释，明确解释了为何需要这些数值策略。
% =========================================================================
clear; close all; clc;
addpath(pwd);
fprintf('=== OLG模型转轨动态求解脚本 (TPI) - v7.0 数值稳定最终版 ===\n\n');

%% --- 1. 加载稳态数据与初始化 ---
fprintf('--- 1. 加载稳态数据与初始化 ---\n');
data_filename = 'SS/data_for_transition_nopps.mat';
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
paramS0 = data_for_transition.paramS0;
paramSF = data_for_transition.paramSF;

if isempty(dist0), error('加载的 dist0 为空。'); end

max_iter = 200; 
tolerance = 1e-6; 

% ======================== [核心修正 1: 阻尼因子] ========================
% 对于内生劳动供给模型，工资的反馈回路非常敏感，需要较小的阻尼因子
damping.r = 0.1;      % 资本市场通常较稳定，可以使用稍大的阻尼
damping.w = 0.1;       % [!!!] 工资路径需要非常保守的更新
damping.beq = 0.1;     % 遗赠和TR都与工资高度相关，同样需要保守
damping.tr = 0.1;
% =======================================================================

T = cS.T_sim;
fprintf('   模拟期数 T = %d\n', T);
fprintf('   [稳定化策略] 使用保守的阻尼因子: r=%.2f, w=%.2f, beq=%.2f, tr=%.2f\n', ...
    damping.r, damping.w, damping.beq, damping.tr);

%% --- 2. 构造宏观路径的初始猜测 ---
fprintf('\n--- 2. 构造宏观路径的初始猜测 ---\n');
r_path_guess = linspace(ss0.r_mkt, ssF.r_mkt, T);
w_hat_path_guess = linspace(ss0.w_hat, ssF.w_hat, T);

total_pop0 = sum(cS.Z_path_raw(:,1));
total_popF = sum(cS.Z_path_raw(:,end));
beq0_per_capita_hat = ss0.Bequest_distributed_agg / total_pop0;
beqF_per_capita_hat = ssF.Bequest_distributed_agg / total_popF;
Bequest_hat_path_guess = linspace(beq0_per_capita_hat, beqF_per_capita_hat, T);
tr0_hat_per_hh = (ss0.TR_distributed_agg / total_pop0) / cS.A_path(1);
trF_hat_per_hh = (ssF.TR_distributed_agg / total_popF) / cS.A_path(end);
TR_hat_path_guess = linspace(tr0_hat_per_hh, trF_hat_per_hh, T);
fprintf('   初始猜测已生成。\n');

%% --- 3. [核心] 时序路径迭代 (TPI) 循环 ---
fprintf('\n--- 3. 启动TPI循环 (最大迭代: %d, 容忍度: %.1e) ---\n', max_iter, tolerance);
r_path_iter = r_path_guess; 
w_hat_path_iter = w_hat_path_guess;
Bequest_hat_path_iter = Bequest_hat_path_guess;
TR_hat_path_iter = TR_hat_path_guess;

total_time_start = tic;
for iter = 1:max_iter
    iter_time_start = tic;
    
    [target_paths, errors] = TPI.calculate_paths_and_errors(...
        r_path_iter, w_hat_path_iter, Bequest_hat_path_iter, TR_hat_path_iter, ...
        ss0, ssF, cS, paramSF, valF, polF, dist0);

    iter_time_elapsed = toc(iter_time_start);
    fprintf('Iter [%3d/%3d] | Error: %.3e (r:%.2e, w:%.2e, bq:%.2e, tr:%.2e) | Time: %.2f s\n', ...
        iter, max_iter, errors.total, errors.r, errors.w, errors.beq, errors.tr, iter_time_elapsed);

    if errors.total < tolerance, fprintf('\n✅ TPI算法在第 %d 次迭代后成功收敛!\n', iter); break; end
    
    % --- [核心修正 2: 阻尼更新与锚定] ---
    % a. 对 t=2...T 的路径进行阻尼更新
    r_path_iter(2:T) = (1 - damping.r) * r_path_iter(2:T) + damping.r * target_paths.r_path(2:T);
    w_hat_path_iter(2:T) = (1 - damping.w) * w_hat_path_iter(2:T) + damping.w * target_paths.w_hat_path(2:T);
    Bequest_hat_path_iter(2:T) = (1 - damping.beq) * Bequest_hat_path_iter(2:T) + damping.beq * target_paths.Bequest_hat_path(2:T);
    TR_hat_path_iter(2:T) = (1 - damping.tr) * TR_hat_path_iter(2:T) + damping.tr * target_paths.TR_hat_path(2:T);

    % b. [!!!] 严格锚定 t=1 的初始条件，防止误差污染
    % 这一步至关重要，它确保了每次迭代都从一个固定的、正确的起点出发。
    r_path_iter(1) = ss0.r_mkt;
    w_hat_path_iter(1) = ss0.w_hat;
    Bequest_hat_path_iter(1) = beq0_per_capita_hat;
    TR_hat_path_iter(1) = tr0_hat_per_hh;
    % --- [修正结束] ---

    if iter == max_iter, warning('TPI在达到最大迭代次数后仍未收敛！'); end
end
total_time_elapsed = toc(total_time_start);
fprintf('--- TPI求解完成，总耗时: %.2f 秒 ---\n', total_time_elapsed);

%% --- 4 & 5. 保存、可视化、检验 ---
% (此部分代码无需修改)
if errors.total < tolerance
    fprintf('\n--- 4. 生成并保存最终转轨路径结果 ---\n');
    if ~exist('TRANS', 'dir'), mkdir('TRANS'); end
    if cS.pps_active
        results_filename = 'TRANS/TPI_results_pps.mat';
    else
        results_filename = 'TRANS/TPI_results_nopps.mat';
    end
    
    results = struct();
    results.r_path = r_path_iter;
    results.w_hat_path = w_hat_path_iter;
    results.Bequest_hat_path = Bequest_hat_path_iter;
    results.TR_hat_path = TR_hat_path_iter;
    
    final_pathS = TPI.fill_auxiliary_paths(results, cS); 
    [final_Pol_path, ~] = household.backward_hh(final_pathS, cS, paramSF, valF, polF);
    final_Dist_path = distribution.simulate_dist_forward(dist0, final_Pol_path, cS, paramSF, final_pathS.beq_transfer_hat_path);
    final_aggr = aggregates.get_path_aggregates(final_Dist_path, final_Pol_path, cS, ss0, paramSF);
    
    results.L_path = final_aggr.L_path;
    results.K_p_hat_path = final_aggr.K_p_hat_path;
    results.C_hat_path = final_aggr.C_hat_path;
    results.w_path = results.w_hat_path .* cS.A_path;
    results.K_p_path = results.K_p_hat_path .* cS.A_path;
    results.C_path = results.C_hat_path .* cS.A_path;
    
    total_pop_path = sum(cS.Z_path_raw, 1);
    results.Bequest_path_total = results.Bequest_hat_path .* total_pop_path;
    results.TR_path_total = results.TR_hat_path .* total_pop_path .* cS.A_path;
    
    results.K_g_path = ssF.K_public_hat .* cS.A_path;
    results.Y_hat_path = (results.K_p_hat_path.^cS.alpha) .* (ssF.K_public_hat.^cS.gamma) .* (results.L_path.^(1-cS.alpha-cS.gamma));
    results.Y_path = results.Y_hat_path .* cS.A_path;

    K_p_end_path = [results.K_p_path(2:T), ssF.K_private_hat * cS.A_path(end) * (1+cS.g_A_ss)^cS.time_Step * (1+cS.n_ss)^cS.time_Step];
    results.I_p_path = K_p_end_path - (1 - cS.ddk) .* results.K_p_path;
    K_g_end_path = [results.K_g_path(2:T), ssF.K_public_hat * cS.A_path(end) * (1+cS.g_A_ss)^cS.time_Step * (1+cS.n_ss)^cS.time_Step];
    results.I_g_path = K_g_end_path - (1-cS.ddk_g) .* results.K_g_path;
    
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