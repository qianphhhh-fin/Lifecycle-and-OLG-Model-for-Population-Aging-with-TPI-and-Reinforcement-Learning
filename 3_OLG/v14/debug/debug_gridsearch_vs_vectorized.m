% =========================================================================
% == SCRIPT: debug_gridsearch_vs_vectorized.m
% == 目的：比较循环版(GridSearch)和向量化版(Vectorized)VFI求解器的
% ==       准确性、一致性和运行速度。
% =========================================================================
clear; clc; close all;
addpath(pwd);

fprintf('--- VFI求解器对比调试脚本 (GridSearch vs. Vectorized) ---\n');

%% 1. 构建一个统一的模拟环境
fprintf('\n--- 1. 构建统一的模拟环境 ---\n');
cS = main_olg_v14_utils.ParameterValues_HuggettStyle();
paramS = struct();

% 使用较大的网格来测试性能差异
ngrid = 15; nprime = 20; npps_choice = 10;
cS.nk = ngrid; cS.nkpps = ngrid; cS.nkprime = nprime; cS.npps = npps_choice;
fprintf('测试网格大小: nk=%d, nkpps=%d, nkprime=%d, npps=%d\n', ngrid, ngrid, nprime, npps_choice);
cS = main_olg_v14_utils.generateGrids(cS);

[paramS.leLogGridV, paramS.leTrProbM, paramS.leProb1V] = main_olg_v14_utils.EarningProcess_olgm(cS);
paramS.leGridV = exp(paramS.leLogGridV(:));

% 设定一些合理的市场价格和参数
M_age = struct('w_t', 1.5, 'r_mkt_period', 0.04, 'r_net_period', 0.04, 'b_t', 0.5);
cS.pps_active = true;
cS.theta_t = 0.11;
tr_per_capita_age = 0.05;
b_age_val = M_age.b_t;

% 伪造一个下一期的价值函数 V'
fprintf('伪造下一期的价值函数 V''...\n');
vPrime_fake = randn(cS.nk, cS.nkpps, cS.nw);

% 选择一个特定的年龄和效率状态进行比较
a_idx_test = 5;  % 测试工作期
% a_idx_test = 12; % 测试退休期
ie_idx_test = 3; % 测试中间效率状态

fprintf('✅ 环境构建完毕。将在年龄 a=%d, 效率状态 ie=%d 上进行对比。\n', a_idx_test, ie_idx_test);

%% 2. 运行两个版本的求解器并计时
fprintf('\n--- 2. 分别运行两个版本的求解器 ---\n');

% a. 运行循环版 (GridSearch)
tic;
[cPol_gs, kPol_gs, cPpsPol_gs, val_gs] = main_olg_v14_utils.HHSolutionByAge_VFI_GridSearch(a_idx_test, vPrime_fake, M_age, tr_per_capita_age, b_age_val, paramS, cS);
time_gs = toc;
fprintf('循环版 (GridSearch) 运行时间: %.4f 秒\n', time_gs);

% b. 运行向量化版 (Vectorized)
tic;
[cPol_vec, kPol_vec, cPpsPol_vec, val_vec] = main_olg_v14_utils.HHSolutionByAge_VFI_Vectorized(a_idx_test, vPrime_fake, M_age, tr_per_capita_age, b_age_val, paramS, cS);
time_vec = toc;
fprintf('向量化版 (Vectorized) 运行时间: %.4f 秒\n', time_vec);

% 计算性能提升
speedup = time_gs / time_vec;
fprintf('速度提升倍数: %.2f x\n', speedup);
if speedup > 1.2
    fprintf('✅ 向量化版本实现了显著加速！\n');
elseif speedup > 1.05
    fprintf('✅ 向量化版本实现了适度加速。\n');
else
    fprintf('⚠️  向量化版本加速效果有限，可能需要更大的网格才能体现优势。\n');
end

% 计算总的计算量作为参考
total_states = cS.nk * cS.nkpps * cS.nw;
total_decisions = cS.nkprime * cS.npps;
theoretical_operations = total_states * total_decisions;
fprintf('\n--- 性能分析 ---\n');
fprintf('总状态点数: %d\n', total_states);
fprintf('每个状态的平均决策组合数: ~%d\n', total_decisions);
fprintf('理论计算量级: ~%.1e 操作\n', theoretical_operations);

%% 3. 深度对比结果的准确性和一致性
fprintf('\n--- 3. 深度对比结果的一致性 ---\n');

% 我们只比较在 ie_idx_test 这个效率状态下的结果
val_gs_slice = val_gs(:,:,ie_idx_test);
val_vec_slice = val_vec(:,:,ie_idx_test);

kPol_gs_slice = kPol_gs(:,:,ie_idx_test);
kPol_vec_slice = kPol_vec(:,:,ie_idx_test);

cPpsPol_gs_slice = cPpsPol_gs(:,:,ie_idx_test);
cPpsPol_vec_slice = cPpsPol_vec(:,:,ie_idx_test);

cPol_gs_slice = cPol_gs(:,:,ie_idx_test);
cPol_vec_slice = cPol_vec(:,:,ie_idx_test);

% 计算详细差异统计
diff_val_matrix = abs(val_gs_slice - val_vec_slice);
diff_kPol_matrix = abs(kPol_gs_slice - kPol_vec_slice);
diff_cPpsPol_matrix = abs(cPpsPol_gs_slice - cPpsPol_vec_slice);
diff_cPol_matrix = abs(cPol_gs_slice - cPol_vec_slice);

fprintf('=== 价值函数差异分析 ===\n');
fprintf('最大绝对差: %e\n', max(diff_val_matrix, [], 'all'));
fprintf('平均绝对差: %e\n', mean(diff_val_matrix, 'all'));
fprintf('差异>1e-6的点数: %d / %d (%.1f%%)\n', sum(diff_val_matrix(:) > 1e-6), numel(diff_val_matrix), sum(diff_val_matrix(:) > 1e-6)/numel(diff_val_matrix)*100);

fprintf('\n=== 储蓄策略 k'' 差异分析 ===\n');
fprintf('最大绝对差: %e\n', max(diff_kPol_matrix, [], 'all'));
fprintf('平均绝对差: %e\n', mean(diff_kPol_matrix, 'all'));
fprintf('差异>1e-6的点数: %d / %d (%.1f%%)\n', sum(diff_kPol_matrix(:) > 1e-6), numel(diff_kPol_matrix), sum(diff_kPol_matrix(:) > 1e-6)/numel(diff_kPol_matrix)*100);

fprintf('\n=== PPS缴费策略差异分析 ===\n');
fprintf('最大绝对差: %e\n', max(diff_cPpsPol_matrix, [], 'all'));
fprintf('平均绝对差: %e\n', mean(diff_cPpsPol_matrix, 'all'));
fprintf('差异>1e-6的点数: %d / %d (%.1f%%)\n', sum(diff_cPpsPol_matrix(:) > 1e-6), numel(diff_cPpsPol_matrix), sum(diff_cPpsPol_matrix(:) > 1e-6)/numel(diff_cPpsPol_matrix)*100);

% 找出差异最大的几个点进行详细分析
[max_kPol_diff, max_idx] = max(diff_kPol_matrix(:));
[max_i, max_j] = ind2sub(size(diff_kPol_matrix), max_idx);

fprintf('\n=== 最大差异点详细分析 ===\n');
fprintf('位置: (k=%d, k_pps=%d) -> (%.4f, %.4f)\n', max_i, max_j, cS.kGridV(max_i), cS.kppsGridV(max_j));
fprintf('GridSearch: k''=%.6f, c_pps=%.6f, val=%.6f\n', kPol_gs_slice(max_i, max_j), cPpsPol_gs_slice(max_i, max_j), val_gs_slice(max_i, max_j));
fprintf('Vectorized: k''=%.6f, c_pps=%.6f, val=%.6f\n', kPol_vec_slice(max_i, max_j), cPpsPol_vec_slice(max_i, max_j), val_vec_slice(max_i, max_j));

% 检查是否有明显的系统性偏差
kPol_bias = mean(kPol_gs_slice - kPol_vec_slice, 'all');
cPpsPol_bias = mean(cPpsPol_gs_slice - cPpsPol_vec_slice, 'all');
val_bias = mean(val_gs_slice - val_vec_slice, 'all');

fprintf('\n=== 系统性偏差检查 ===\n');
fprintf('储蓄策略平均偏差: %e\n', kPol_bias);
fprintf('PPS缴费策略平均偏差: %e\n', cPpsPol_bias);
fprintf('价值函数平均偏差: %e\n', val_bias);

% 判断结果
tolerance_high = 1e-4;
tolerance_low = 1e-6;

if max(diff_kPol_matrix, [], 'all') < tolerance_high && max(diff_cPpsPol_matrix, [], 'all') < tolerance_high
    fprintf('\n✅ 结论：两个版本的策略函数结果高度一致！向量化版本正确。\n');
elseif mean(diff_kPol_matrix, 'all') < tolerance_low && mean(diff_cPpsPol_matrix, 'all') < tolerance_low
    fprintf('\n⚠️  结论：两个版本大体一致，但存在少数异常点。可能是数值精度或边界处理差异。\n');
else
    fprintf('\n❌ 警告：两个版本的策略函数结果存在显著系统性差异！需要检查向量化版本的逻辑。\n');
end

%% 4. 差异热力图可视化
figure('Name', '差异分析热力图', 'Position', [50, 50, 1200, 400]);

subplot(1, 3, 1);
imagesc(cS.kGridV, cS.kppsGridV, diff_val_matrix');
colorbar; title('价值函数绝对差异'); xlabel('k'); ylabel('k_{pps}');
set(gca, 'YDir', 'normal');

subplot(1, 3, 2);
imagesc(cS.kGridV, cS.kppsGridV, diff_kPol_matrix');
colorbar; title("储蓄策略 k' 绝对差异"); xlabel('k'); ylabel('k_{pps}');
set(gca, 'YDir', 'normal');

subplot(1, 3, 3);
imagesc(cS.kGridV, cS.kppsGridV, diff_cPpsPol_matrix');
colorbar; title("PPS缴费策略 c_{pps} 绝对差异"); xlabel('k'); ylabel('k_{pps}');
set(gca, 'YDir', 'normal');

sgtitle('差异分布热力图 (颜色越深差异越大)', 'FontSize', 12);

%% 5. 策略函数对比可视化
figure('Name', 'GridSearch vs. Vectorized 对比', 'Position', [100, 100, 1600, 800]);

% 价值函数对比
subplot(2, 4, 1);
surf(cS.kGridV, cS.kppsGridV, val_gs_slice'); view(30,30);
title('价值函数 V (GridSearch)'); xlabel('k'); ylabel('k_{pps}');
subplot(2, 4, 5);
surf(cS.kGridV, cS.kppsGridV, val_vec_slice'); view(30,30);
title('价值函数 V (Vectorized)'); xlabel('k'); ylabel('k_{pps}');

% 储蓄策略 k' 对比
subplot(2, 4, 2);
surf(cS.kGridV, cS.kppsGridV, kPol_gs_slice'); view(30,30);
title("储蓄策略 k' (GridSearch)"); xlabel('k'); ylabel('k_{pps}');
subplot(2, 4, 6);
surf(cS.kGridV, cS.kppsGridV, kPol_vec_slice'); view(30,30);
title("储蓄策略 k' (Vectorized)"); xlabel('k'); ylabel('k_{pps}');

% PPS缴费策略 c_pps 对比
subplot(2, 4, 3);
surf(cS.kGridV, cS.kppsGridV, cPpsPol_gs_slice'); view(30,30);
title("PPS缴费策略 c_{pps} (GridSearch)"); xlabel('k'); ylabel('k_{pps}');
subplot(2, 4, 7);
surf(cS.kGridV, cS.kppsGridV, cPpsPol_vec_slice'); view(30,30);
title("PPS缴费策略 c_{pps} (Vectorized)"); xlabel('k'); ylabel('k_{pps}');

% 消费策略 c 对比
subplot(2, 4, 4);
surf(cS.kGridV, cS.kppsGridV, cPol_gs_slice'); view(30,30);
title("消费策略 c (GridSearch)"); xlabel('k'); ylabel('k_{pps}');
subplot(2, 4, 8);
surf(cS.kGridV, cS.kppsGridV, cPol_vec_slice'); view(30,30);
title("消费策略 c (Vectorized)"); xlabel('k'); ylabel('k_{pps}');

sgtitle(sprintf('求解器对比: a=%d, ie=%d', a_idx_test, ie_idx_test), 'FontSize', 14, 'FontWeight', 'bold');

%% 6. 多网格大小性能测试
fprintf('\n--- 6. 多网格大小性能对比 ---\n');
grid_sizes = [8, 12, 16, 20];
speedups = zeros(size(grid_sizes));

for i = 1:length(grid_sizes)
    ngrid_test = grid_sizes(i);
    fprintf('\n测试网格大小: %dx%d\n', ngrid_test, ngrid_test);
    
    % 重新设置参数
    cS_test = cS;
    cS_test.nk = ngrid_test; cS_test.nkpps = ngrid_test;
    cS_test = main_olg_v14_utils.generateGrids(cS_test);
    
    % 创建对应大小的价值函数
    vPrime_test = randn(ngrid_test, ngrid_test, cS.nw);
    
    % 测试GridSearch
    tic;
    [~, ~, ~, ~] = main_olg_v14_utils.HHSolutionByAge_VFI_GridSearch(a_idx_test, vPrime_test, M_age, tr_per_capita_age, b_age_val, paramS, cS_test);
    time_gs_test = toc;
    
    % 测试Vectorized
    tic;
    [~, ~, ~, ~] = main_olg_v14_utils.HHSolutionByAge_VFI_Vectorized(a_idx_test, vPrime_test, M_age, tr_per_capita_age, b_age_val, paramS, cS_test);
    time_vec_test = toc;
    
    speedups(i) = time_gs_test / time_vec_test;
    fprintf('GridSearch: %.4fs, Vectorized: %.4fs, 加速倍数: %.2fx\n', time_gs_test, time_vec_test, speedups(i));
end

% 绘制性能曲线
figure('Name', '性能对比曲线', 'Position', [100, 600, 800, 400]);
subplot(1, 2, 1);
plot(grid_sizes.^2, speedups, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('网格点数 (nk × nkpps)');
ylabel('速度提升倍数');
title('向量化加速效果 vs 网格大小');
grid on;
ylim([0.8, max(speedups)*1.1]);
yline(1, 'Color', 'r', 'LineStyle', '--', 'Alpha', 0.7);

subplot(1, 2, 2);
theoretical_ops = (grid_sizes.^2 * cS.nw) .* (cS.nkprime * cS.npps);
loglog(theoretical_ops, speedups, 'ro-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('理论计算量');
ylabel('速度提升倍数');
title('加速效果 vs 计算复杂度');
grid on;

% 总结
fprintf('\n--- 性能总结 ---\n');
[max_speedup, max_idx] = max(speedups);
fprintf('最大加速倍数: %.2fx (网格大小: %dx%d)\n', max_speedup, grid_sizes(max_idx), grid_sizes(max_idx));
if max_speedup > 1.5
    fprintf('✅ 向量化版本在大网格下实现了显著加速！\n');
elseif max_speedup > 1.2
    fprintf('✅ 向量化版本实现了适度加速。\n');
else
    fprintf('⚠️  向量化版本加速效果有限，可能需要进一步优化。\n');
end