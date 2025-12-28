% ========================================================================
% == SCRIPT: test_diffN_trans.m
% == 版本: [v2.1 - 人均量插值修正版]
% ==
% == 目的:
% ==   - [!!! 关键逻辑修正 !!!] 纠正了构造初始猜测路径的严重错误。
% ==   - 由于 ss0 和 ssF 是在总人口为1的归一化环境下求解的，因此
% ==     ss0.TR_distributed_agg 和 ss0.Bequest_gen_hat_raw_ss 等
% ==     变量实际上是【人均量】。
% ==   - 新的逻辑直接对这些【人均量】进行线性插值，以构造正确的
% ==     【人均量猜测路径】，不再进行任何错误的总量转换。
% ==   - 这确保了TPI循环的初始输入在量级和单位上都是正确的。
% =========================================================================
clear; close all; clc;
addpath(pwd);
fprintf('=== 受控人口动态转轨求解脚本 (v2.1 - 人均量插值修正版) ===\n\n');

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

max_iter = 200;
tolerance = 1e-6;
T = cS.T_sim;

damping = struct('r', 0.2, 'w', 0.2, 'beq', 0.2, 'tr', 0.2, 'b_hat', 0.2);

fprintf('   模拟期数 T = %d\n', T);
fprintf('   阻尼因子: r=%.2f, w=%.2f, bq=%.2f, tr=%.2f, b_hat=%.2f\n', ...
    damping.r, damping.w, damping.beq, damping.tr, damping.b_hat);

%% --- 2. 构造宏观路径的初始猜测 (人均量插值修正版) ---
fprintf('\n--- 2. 构造【人均量插值】的宏观路径初始猜测 ---\n');

% 2a. 对价格和养老金福利(均为 интенсивные величины)进行插值
r_path_guess = linspace(ss0.r_mkt, ssF.r_mkt, T);
w_hat_path_guess = linspace(ss0.w_hat, ssF.w_hat, T);
b_hat_path_guess = linspace(ss0.b_hat, ssF.b_hat, T);
fprintf('   ✅ 价格和养老金路径已构造。\n');

% [!!! 关键逻辑修正 !!!]
% ss0和ssF是在总人口为1的环境下解出的，所以带_agg或_raw后缀的变量均为人均量。
% 我们应该直接对这些人均量进行插值。

% 2b. 构造【人均遗赠】猜测路径 (bequest_hat_pc_path)
bequest_pc_ss0 = ss0.Bequest_gen_hat_raw_ss; % 这是人均原始遗赠
bequest_pc_ssF = ssF.Bequest_gen_hat_raw_ss; % 这是人均原始遗赠
bequest_hat_pc_path_guess = linspace(bequest_pc_ss0, bequest_pc_ssF, T);
fprintf('   ✅ 人均遗赠猜测路径已构造。\n');

% 2c. 构造【人均转移支付】猜测路径 (TR_hat_pc_path_guess)
tr_pc_ss0 = ss0.TR_distributed_agg; % 这是人均转移支付
tr_pc_ssF = ssF.TR_distributed_agg; % 这是人均转移支付
TR_hat_pc_path_guess = linspace(tr_pc_ss0, tr_pc_ssF, T);
fprintf('   ✅ 人均转移支付猜测路径已构造。\n');


%% --- 3. [核心] 时序路径迭代 (TPI) 循环 ---
fprintf('\n--- 3. 启动TPI循环 (最大迭代: %d, 容忍度: %.1e) ---\n', max_iter, tolerance);
r_path_iter = r_path_guess;
w_hat_path_iter = w_hat_path_guess;
bequest_hat_pc_path_iter = bequest_hat_pc_path_guess;
TR_hat_pc_path_iter = TR_hat_pc_path_guess;
b_hat_path_iter = b_hat_path_guess;

total_time_start = tic;
for iter = 1:max_iter
    iter_time_start = tic;

    [target_paths, errors, ~, ~, aggr_supply] = TPI.calculate_paths_and_errors(...
        r_path_iter, w_hat_path_iter, bequest_hat_pc_path_iter, ...
        TR_hat_pc_path_iter, b_hat_path_iter, ...
        ss0, ssF, cS, paramSF, valF, polF, dist0);

    iter_time_elapsed = toc(iter_time_start);
    fprintf('Iter [%3d/%3d] | Error: %.3e (r:%.2e, w:%.2e, bq:%.2e, tr:%.2e, b:%.2e) | Time: %.2f s\n', ...
        iter, max_iter, errors.total, errors.r, errors.w, errors.beq, errors.tr, errors.b_hat, iter_time_elapsed);

    if errors.total < tolerance, fprintf('\n✅ TPI算法在第 %d 次迭代后成功收敛!\n', iter); break; end

    update_range = 2:(T);
    r_path_iter(update_range) = (1 - damping.r) * r_path_iter(update_range) + damping.r * target_paths.r_path(update_range);
    w_hat_path_iter(update_range) = (1 - damping.w) * w_hat_path_iter(update_range) + damping.w * target_paths.w_hat_path(update_range);
    TR_hat_pc_path_iter(update_range) = (1 - damping.tr) * TR_hat_pc_path_iter(update_range) + damping.tr * target_paths.TR_hat_pc_path(update_range);
    b_hat_path_iter(update_range) = (1 - damping.b_hat) * b_hat_path_iter(update_range) + damping.b_hat * target_paths.b_hat_path(update_range);
    bequest_hat_pc_path_iter(update_range) = (1 - damping.beq) * bequest_hat_pc_path_iter(update_range) + damping.beq * target_paths.bequest_hat_pc_path(update_range);

    if iter == max_iter, warning('TPI在达到最大迭代次数后仍未收敛！'); end
end
total_time_elapsed = toc(total_time_start);
fprintf('--- TPI求解完成，总耗时: %.2f 秒 ---\n', total_time_elapsed);

%% --- 4. 结果整理与可视化 ---
if errors.total < tolerance
    fprintf('\n--- 4. 结果整理与可视化 ---\n');
    
    % 在收敛后，再完整运行一次以获取最终的、一致的路径和聚合量
    [final_paths, ~, ~, ~, final_aggr_supply] = TPI.calculate_paths_and_errors(...
        r_path_iter, w_hat_path_iter, bequest_hat_pc_path_iter, ...
        TR_hat_pc_path_iter, b_hat_path_iter, ...
        ss0, ssF, cS, paramSF, valF, polF, dist0);
    
    results = struct();
    total_pop_path = sum(cS.Z_path_raw, 1);
    results.r_path = final_paths.r_path;
    results.w_path = final_paths.w_hat_path .* cS.A_path;
    results.K_p_path = final_aggr_supply.K_p_hat_total_path .* cS.A_path;
    results.C_path = final_aggr_supply.C_hat_pc_path .* total_pop_path .* cS.A_path;

    % 构造国民账户所需的其他变量
    Y_hat_total_path_supply = (final_paths.r_path + cS.ddk) .* final_aggr_supply.K_p_hat_total_path + final_paths.w_hat_path .* final_aggr_supply.L_path; % 这是一个简化的Y，仅用于绘图和检查
    results.Y_path = Y_hat_total_path_supply .* cS.A_path;
    
    if ismethod('TPI', 'check_transition_NIPA')
        % 为 NIPA 检查构造更精确的变量
        K_g_hat_path = ssF.K_public_hat * ones(1, T); % 假设公共资本固定在终期水平
        Y_from_prod_hat_path = zeros(1, T);
        for t=1:T
            prices_t = firm.get_prices_at_t(final_aggr_supply.K_p_hat_total_path(t), K_g_hat_path(t), final_aggr_supply.L_path(t), cS);
            Y_from_prod_hat_path(t) = prices_t.Y_hat_t;
        end
        results.Y_path = Y_from_prod_hat_path .* cS.A_path;
        
        K_p_path_level = final_aggr_supply.K_p_hat_total_path .* cS.A_path;
        results.I_p_path = [K_p_path_level(2:end) - (1-cS.ddk) * K_p_path_level(1:end-1), nan];
        results.I_p_path(end) = results.I_p_path(end-1); % 简单外插
        results.I_g_path = cS.I_g_to_Y_ratio_ss * results.Y_path;
        
        TPI.check_transition_NIPA(results, cS, ssF);
    end

    if ismethod('TPI', 'plot_transition_results')
        TPI.plot_transition_results(results, cS, ssF, '受控人口冲击实验');
    else
        fprintf('   (可视化函数 TPI.plot_transition_results 未找到，跳过绘图。)\n');
    end
    
    fprintf('   ✅ 结果整理完毕。\n');
else
    fprintf('\n--- 4. 求解失败 ---\n');
    fprintf('   ❌ TPI求解器未能收敛，无法生成可靠结果。\n');
end

fprintf('\n--- 受控人口动态转轨求解脚本执行完毕 ---\n');