% =========================================================================
% == 脚本: test_matrix_for.m (V2.0 - 修正网格生成逻辑)
% == 功能: 
% ==   1. 比较 value_function_iteration (循环版) 和 
% ==      value_function_iteration_matrix (矩阵版) 的输出结果是否一致。
% ==   2. 量化矩阵化带来的性能提升。
% == 核心修正:
% ==   在修改 nW, nF 等网格点数后，必须手动重新生成 wGridV, fGridV 等
% ==   网格向量，以确保维度匹配。
% =========================================================================

clear all;
close all;
clc;

fprintf('正在初始化模型...\n');
addpath(pwd); % 确保utils.m在路径中

% --- 1. 初始化并设置一个"微型"模型用于快速测试 ---
cS = utils.set_parameters();

% 覆盖默认参数以进行快速测试
cS.nW = 7;       % 财富网格点
cS.nF = 4;       % 养老金网格点
cS.nC = 7;       % 消费/储蓄选择网格点
cS.nAlpha = 5;   % 风险资产配置网格点
cS.nQ = 4;       % 养老金缴费率网格点
cS.nShocks = 3;  % 随机冲击离散点数

% --- [!!! 核心修正: 重新生成网格向量 !!!] ---
cS.wGridV = exp(linspace(log(cS.w_min), log(cS.w_max), cS.nW))';
temp_fGridV = exp(linspace(log(cS.f_min+1e-6), log(cS.f_max), cS.nF-1))';
cS.fGridV = [0; temp_fGridV];
cS.alphaGridV = linspace(0, 1, cS.nAlpha)';
cS.qGridV = linspace(0, cS.Q_max, cS.nQ)';
cS.savingsFracGridV = linspace(0.001, 0.999, cS.nC)';

m_std = 3;
[zNodes, zProbMat] = utils.tauchen(cS.nShocks, 0, cS.sigma_z, 0, m_std);
[uNodes, ~] = utils.tauchen(cS.nShocks, 0, cS.sigma_u, 0, m_std);
[epsNodes, ~] = utils.tauchen(cS.nShocks, 0, cS.sigma_eps, 0, m_std);
cS.zNodes = zNodes;
cS.uNodes = uNodes;
cS.epsNodes = epsNodes;
cS.shockProbs = zProbMat(1, :)';
cS.R_shock_V = cS.R_f + cS.mu + cS.epsNodes;
% --- [!!! 修正结束 !!!] ---


fprintf('测试将在以下"微型"网格上进行:\n');
fprintf('  nW=%d, nF=%d, nC=%d, nAlpha=%d, nQ=%d, nShocks=%d\n', ...
        cS.nW, cS.nF, cS.nC, cS.nAlpha, cS.nQ, cS.nShocks);

% --- 2. 运行并计时：value_function_iteration (循环版) ---
fprintf('\n============================================================\n');
fprintf('开始运行: value_function_iteration (循环版)...\n');
fprintf('这是一个基准测试，可能需要几分钟，请耐心等待...\n');
tic;
[polS_for, valS_for] = utils.value_function_iteration(cS);
time_for = toc;
fprintf('循环版 VFI 完成。耗时: %.4f 秒\n', time_for);


% --- 3. 运行并计时：value_function_iteration_matrix (矩阵版) ---
fprintf('\n============================================================\n');
fprintf('开始运行: value_function_iteration_matrix (矩阵版)...\n');
tic;
% 确保您已经在 utils.m 中将新函数命名为 value_function_iteration_matrix
[polS_matrix, valS_matrix] = utils.value_function_iteration_matrix(cS);
time_matrix = toc;
fprintf('矩阵版 VFI 完成。耗时: %.4f 秒\n', time_matrix);


% --- 4. 比较结果 ---
fprintf('\n============================================================\n');
fprintf('正在比较两个函数的输出结果...\n');

tol = 1e-8; % 设置一个很小的容差来处理浮点数误差
results_are_identical = true;

% 比较 valS
max_diff_valS = max(abs(valS_for - valS_matrix), [], 'all', 'omitnan');
fprintf('值函数 (valS) 的最大绝对差值: %e\n', max_diff_valS);
if max_diff_valS > tol
    results_are_identical = false;
    fprintf('  [!!] 警告: valS 差异超过容差!\n');
end

% 比较 polS
for t = 1 : cS.tn
    % 比较消费策略 c
    if ~isempty(polS_for.c{t}) && ~isempty(polS_matrix.c{t})
        max_diff_c = max(abs(polS_for.c{t} - polS_matrix.c{t}), [], 'all', 'omitnan');
        if max_diff_c > tol
            results_are_identical = false;
            fprintf('  [!!] 警告: 消费策略(c)在 t=%d 时差异超过容差 (MaxDiff: %e)\n', t, max_diff_c);
        end
    end
    
    % 比较风险资产配置策略 alpha
    if ~isempty(polS_for.alpha{t}) && ~isempty(polS_matrix.alpha{t})
        max_diff_alpha = max(abs(polS_for.alpha{t} - polS_matrix.alpha{t}), [], 'all', 'omitnan');
        if max_diff_alpha > tol
            results_are_identical = false;
            fprintf('  [!!] 警告: 投资策略(alpha)在 t=%d 时差异超过容差 (MaxDiff: %e)\n', t, max_diff_alpha);
        end
    end
    
    % 比较养老金缴费策略 q
    if ~isempty(polS_for.q{t}) && ~isempty(polS_matrix.q{t})
        max_diff_q = max(abs(polS_for.q{t} - polS_matrix.q{t}), [], 'all', 'omitnan');
        if max_diff_q > tol
            results_are_identical = false;
            fprintf('  [!!] 警告: 缴费策略(q)在 t=%d 时差异超过容差 (MaxDiff: %e)\n', t, max_diff_q);
        end
    end
end

if results_are_identical
    fprintf('\n[SUCCESS] 所有值函数和策略函数的差异均在容差 (%e) 范围内。\n', tol);
else
    fprintf('\n[FAILURE] 函数输出存在显著差异，请检查代码逻辑。\n');
end


% --- 5. 报告性能 ---
fprintf('\n============================================================\n');
fprintf('性能报告:\n');
fprintf('  - 循环版 (for-loop) 耗时: %.4f 秒\n', time_for);
fprintf('  - 矩阵版 (matrix)  耗时: %.4f 秒\n', time_matrix);
if time_matrix > 0.001
    speedup_factor = time_for / time_matrix;
    fprintf('\n[结论] 矩阵化版本是循环版本的 %.2f 倍快。\n', speedup_factor);
else
    fprintf('\n[结论] 矩阵版执行时间过短，无法计算有效的提速比。\n');
end
fprintf('============================================================\n');