% =========================================================================
% == SCRIPT: test_diffN_trans.m
% == 版本: [v1.4 - 最终人均化接口版]
% ==
% == 目的:
% ==   - 更新对TPI核心引擎的调用，使用新的、严格人均化的变量
% ==     (Bequest_hat_pc_path, TR_hat_pc_path) 进行迭代。
% ==   - 确保初始猜测的构造与新接口完全匹配。
% =========================================================================
clear; close all; clc;
addpath(pwd);
fprintf('=== 受控人口动态转轨求解脚本 (v1.4 - 最终人均化接口版) ===\n\n');

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

damping.r = 0.2; damping.w = 0.2; damping.beq = 0.2; damping.tr = 0.2; damping.b_hat = 0.2;

T = cS.T_sim;
fprintf('   模拟期数 T = %d\n', T);
fprintf('   阻尼因子: r=%.2f, w=%.2f, bq=%.2f, tr=%.2f, b_hat=%.2f\n', ...
    damping.r, damping.w, damping.beq, damping.tr, damping.b_hat);

%% --- 2. 构造宏观路径的初始猜测 (完全人均化) ---
fprintf('\n--- 2. 构造【完全人均化】的宏观路径初始猜测 ---\n');
r_path_guess = linspace(ss0.r_mkt, ssF.r_mkt, T);
w_hat_path_guess = linspace(ss0.w_hat, ssF.w_hat, T);
b_hat_path_guess = linspace(ss0.b_hat, ssF.b_hat, T);

% [!!! 核心修正: 所有计算都基于总人口进行人均化 !!!]
total_pop_ss0 = sum(dist0, 'all'); % 初始绝对总人口
total_pop_ssF = sum(cS.Z_path_raw(:,end)); % 终期绝对总人口

if total_pop_ss0 < 1e-9, total_pop_ss0 = 1; end
if total_pop_ssF < 1e-9, total_pop_ssF = 1; end

% 人均标准化遗赠 (总遗赠 / 总人口)
beq0_hat_pc = ss0.Bequest_generated_agg / total_pop_ss0;
beqF_hat_pc = ssF.Bequest_generated_agg / total_pop_ssF;
Bequest_hat_pc_path_guess = linspace(beq0_hat_pc, beqF_hat_pc, T);

% 人均标准化转移支付 (总TR / 总人口)
tr0_hat_pc = ss0.TR_distributed_agg / total_pop_ss0;
trF_hat_pc = ssF.TR_distributed_agg / total_pop_ssF;
TR_hat_pc_path_guess = linspace(tr0_hat_pc, trF_hat_pc, T);

fprintf('   ✅ 路径初始猜测已生成 (人均标准化)。\n');

%% --- 3. [核心] 时序路径迭代 (TPI) 循环 ---
fprintf('\n--- 3. 启动TPI循环 (最大迭代: %d, 容忍度: %.1e) ---\n', max_iter, tolerance);
r_path_iter = r_path_guess;
w_hat_path_iter = w_hat_path_guess;
Bequest_hat_pc_path_iter = Bequest_hat_pc_path_guess;
TR_hat_pc_path_iter = TR_hat_pc_path_guess;
b_hat_path_iter = b_hat_path_guess;

total_time_start = tic;
for iter = 1:max_iter
    iter_time_start = tic;

    % [!!! 核心修正: 调用新的人均化接口 !!!]
    [target_paths, errors] = TPI.calculate_paths_and_errors(...
        r_path_iter, w_hat_path_iter, Bequest_hat_pc_path_iter, ...
        TR_hat_pc_path_iter, b_hat_path_iter, ...
        ss0, ssF, cS, paramSF, valF, polF, dist0);

    iter_time_elapsed = toc(iter_time_start);
    fprintf('Iter [%3d/%3d] | Error: %.3e (r:%.2e, w:%.2e, bq:%.2e, tr:%.2e, b:%.2e) | Time: %.2f s\n', ...
        iter, max_iter, errors.total, errors.r, errors.w, errors.beq, errors.tr, errors.b_hat, iter_time_elapsed);

    if errors.total < tolerance, fprintf('\n✅ TPI算法在第 %d 次迭代后成功收敛!\n', iter); break; end

    update_range = 2:(T-1);
    r_path_iter(update_range) = (1 - damping.r) * r_path_iter(update_range) + damping.r * target_paths.r_path(update_range);
    w_hat_path_iter(update_range) = (1 - damping.w) * w_hat_path_iter(update_range) + damping.w * target_paths.w_hat_path(update_range);
    Bequest_hat_pc_path_iter(update_range) = (1 - damping.beq) * Bequest_hat_pc_path_iter(update_range) + damping.beq * target_paths.Bequest_hat_pc_path(update_range);
    TR_hat_pc_path_iter(update_range) = (1 - damping.tr) * TR_hat_pc_path_iter(update_range) + damping.tr * target_paths.TR_hat_pc_path(update_range);
    b_hat_path_iter(update_range) = (1 - damping.b_hat) * b_hat_path_iter(update_range) + damping.b_hat * target_paths.b_hat_path(update_range);

    r_path_iter(1) = ss0.r_mkt; w_hat_path_iter(1) = ss0.w_hat;
    Bequest_hat_pc_path_iter(1) = beq0_hat_pc;
    TR_hat_pc_path_iter(1) = tr0_hat_pc; b_hat_path_iter(1) = ss0.b_hat;
    
    r_path_iter(T) = ssF.r_mkt; w_hat_path_iter(T) = ssF.w_hat;
    Bequest_hat_pc_path_iter(T) = beqF_hat_pc;
    TR_hat_pc_path_iter(T) = trF_hat_pc; b_hat_path_iter(T) = ssF.b_hat;

    if iter == max_iter, warning('TPI在达到最大迭代次数后仍未收敛！'); end
end
total_time_elapsed = toc(total_time_start);
fprintf('--- TPI求解完成，总耗时: %.2f 秒 ---\n', total_time_elapsed);

%% --- 4. 结果分析与检验 ---
if errors.total < tolerance
    fprintf('\n--- 4. 检验结论 ---\n');
    fprintf('   ✅ 结论：在完全人均化的框架下，TPI求解器成功收敛。\n');
else
    fprintf('\n--- 4. 检验结论 ---\n');
    fprintf('   ❌ 警告：在完全人均化的框架下，TPI求解器未能收敛。\n');
end

fprintf('\n--- 受控人口动态转轨求解脚本执行完毕 ---\n');