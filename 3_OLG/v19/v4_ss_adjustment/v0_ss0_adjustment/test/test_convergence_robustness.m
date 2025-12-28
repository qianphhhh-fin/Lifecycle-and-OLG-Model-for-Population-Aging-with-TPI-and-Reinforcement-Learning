% =========================================================================
% == SCRIPT: test_convergence_robustness.m
% == 版本: [v1.0]
% ==
% == 目的:
% ==   1. 对稳态求解器的收敛性进行系统的稳健性检验。
% ==   2. 遍历多个关键参数(beta, sigma, phi_bequest, theta, g_A, n_ss)
% ==      的不同组合。
% ==   3. 使用并行计算(parfor)加速测试过程。
% ==   4. 记录每种参数组合下的最终均衡误差，并报告求解失败的情况。
% ==   5. 使用 parfor_progress.m 提供进度条监控。
% ==
% == 如何使用:
% ==   1. 确保你已经成功运行了 `main_run_SS.m` 并生成了
% ==      `SS/data_for_transition.mat` 文件。
% ==   2. 确保 `parfor_progress.m` 文件在你的MATLAB路径中。
% ==      (可从 MATLAB File Exchange 下载)。
% ==   3. 直接运行此脚本。
% =========================================================================

clear; close all; clc;
addpath(pwd);
fprintf('=== 稳态求解器收敛性稳健性测试脚本 (v1.0) ===\n\n');

%% --- 1. 全局设置与初始化 ---

cS = utils.ParameterValues();

cS.nk = 30;

% 设定模拟的起始年份
cS.ss0_year = 2023;
cS.start_year = 2023;
cS.pop_data_last_year = 2030;


%% --- 2. [核心] 生成完整的外生路径 ---
[Z_path, Z_path_raw, cS] = population.generate_Z_path(cS, false); % 可视化检查
% Z_path_raw = Z_path_raw/sum(Z_path_raw(:,1)); % 将Z_path_raw的第一期总和归一化为1
T_sim = size(Z_path,2);
% A_path = utils.generate_tfp_path(cS,true);
g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
cS.A_path = (1 + g_A_period).^(0:(T_sim-1));

cS = utils.calcaulte_theta_payg_path(cS, false); % 可视化检查
theta_path = cS.theta_path;

[~, ~, ~, cS.nw_expanded] = utils.EarningProcess_AgeDependent(cS);

% 将核心路径存入cS，以便传递
cS.Z_path = Z_path; % 归一化分布
cS.Z_path_raw = Z_path_raw; % 绝对人口

%% --- 3. 求解稳态 (ss0, t=0) ---
fprintf('\n\n--- 3. 求解初始伪稳态 (ss0, t=0) ---\n');
fprintf('   特征: %d年人口结构, 无PPS制度\n', cS.ss0_year);

cS_base = cS;
cS_base.pps_active = false; % 初始稳态不含PPS
cS_base.nkpps = 1; cS_base.n_pps_rate_grid = 1;
cS_base = utils.generateGrids(cS_base);

cS_base.n_ss = ((sum(Z_path_raw(:,2))/sum(Z_path_raw(:,1)))^(1/cS.time_Step))-1;
cS_base.s_pathV = cS.s_pathV(:,1); % % 使用第一期的生存率,这个不设置的话就会出现“警告: 下界的长度小于 length(x)；将用 -Inf 填充缺失的下界。 ”
cS_base.theta_path = theta_path(1);
% cS_base.phi_bequest = 0.5;


% [关键] 使用与转轨完全一致的 t=1 人口分布 (绝对量)
Z0_abs = Z_path_raw(:,1);
Z_base = Z0_abs / sum(Z0_abs); % 求解器需要归一化分布

paramS_base = struct();
[paramS_base.leGridV, paramS_base.TrProbM_by_age, paramS_base.leProb1V, ~] = utils.EarningProcess_AgeDependent(cS_base);


% --- 定义要测试的参数网格 ---
beta_vec        = [0.98];         % 折现因子 (年化)
sigma_vec       = [2.0, 3.0, 4.0, 5.0,6.0,7.0];            % 风险厌恶系数
phi_bequest_vec = [0.1];            % 遗赠动机
theta_vec       = [0.06, 0.11, 0.15];         % PAYG 缴费率
g_A_ss_vec      = [0.02];               % 技术增长率 (年化)
n_ss_vec        = [cS_base.n_ss];              % 人口增长率 (年化)

% 使用 ndgrid 创建所有参数组合
[beta_grid, sigma_grid, phi_grid, theta_grid, gA_grid, n_grid] = ...
    ndgrid(beta_vec, sigma_vec, phi_bequest_vec, theta_vec, g_A_ss_vec, n_ss_vec);

% 将网格展平为向量，便于 parfor 迭代
num_sims = numel(beta_grid);
beta_flat        = beta_grid(:);
sigma_flat       = sigma_grid(:);
phi_flat         = phi_grid(:);
theta_flat       = theta_grid(:);
gA_flat          = gA_grid(:);
n_flat           = n_grid(:);

fprintf('   将对 %d 种不同的参数组合进行测试。\n', num_sims);
fprintf('   测试可能需要较长时间，请耐心等待。\n');


%% --- 2. 并行计算所有参数组合的稳态解 ---
fprintf('\n--- 2. 启动并行计算 ---\n');

% 初始化结果存储向量
final_errors_flat = zeros(num_sims, 1);

% 启动并行池 (如果尚未启动)
% 初始化进度条
try
    parfor_progress(num_sims);
catch
    warning('parfor_progress.m 未找到。将不会显示进度条。请从MATLAB File Exchange下载。');
end


tic;
for i = 1:num_sims
    % a. 为当前迭代创建一个独立的 cS 结构体副本
    cS_iter = cS_base;

    % b. 将当前组合的参数赋值给 cS_iter
    cS_iter.beta = beta_flat(i);
    cS_iter.sigma = sigma_flat(i);
    cS_iter.phi_bequest = phi_flat(i);
    cS_iter.theta_path = theta_flat(i); % 注意 theta_path 被设为标量
    cS_iter.g_A_ss = gA_flat(i);
    cS_iter.n_ss = n_flat(i);
    
    % c. 准备传递给求解器的外部参数
    params_ext_iter = struct(...
        'Z', Z_base, ...
        'A', 1.0, ...
        'g_A_ss', cS_iter.g_A_ss, ...
        'n_ss', cS_iter.n_ss);

    max_error = NaN; % 默认值为 NaN，表示失败

    try
        % d. 调用稳态求解器 (关闭详细输出以保持命令行清洁)
        [ss_result, ~, ~, ~] = SS.solve_steady_state(cS_iter, paramS_base, params_ext_iter, false);
        
        % e. 如果求解成功，则重新计算最终误差以确保一致性
        if ~isempty(ss_result)
            mass_total = sum(Z_base);
            bq_per_capita = ss_result.Bequest_distributed_agg / mass_total;
            tr_per_capita = ss_result.TR_distributed_agg / mass_total;

            x_eq_final = [ss_result.r_mkt, ss_result.w_hat, bq_per_capita, tr_per_capita, ss_result.b_hat];
            
            % 获取最终的误差向量
            errors_vec = SS.system_of_equations(x_eq_final, Z_base, cS_iter, paramS_base, params_ext_iter);
            
            % 记录所有市场中的最大绝对误差
            max_error = max(abs(errors_vec));
        end

    catch ME
        % f. 如果求解器本身抛出异常，也记录为失败
        % (在并行循环中，不建议打印大量错误信息，所以我们只记录结果)
        max_error = NaN;
    end
    
    % g. 存储当前参数组合下的最大误差
    final_errors_flat(i) = max_error;
    
    % h. 更新进度条
    try
        parfor_progress;
    catch
        % 如果没有进度条函数，则不执行任何操作
    end
end
elapsed_time = toc;

% 关闭进度条
try
    parfor_progress(0);
catch
end


%% --- 3. 分析并报告结果 ---
fprintf('\n\n--- 3. 测试结果分析 ---\n');
fprintf('   并行计算完成，总耗时: %.2f 秒 (%.2f 分钟).\n', elapsed_time, elapsed_time/60);

% --- 用户设定 ---
convergence_tol = 1e-7; % <--- 在这里设定你的收敛容忍度
max_results_to_display = 100; % 最多显示多少条成功的结果

% --- 统计与分析 ---
% 统计求解失败的次数
num_failures = sum(isnan(final_errors_flat));
num_sims_attempted = num_sims;

fprintf('\n   总共尝试了 %d 种参数组合。\n', num_sims_attempted);
fprintf('   其中，求解器未能成功返回结果 (失败) 的有 %d 个 (%.2f %%)\n', ...
    num_failures, (num_failures/num_sims_attempted)*100);

% 找出所有误差低于容忍度的 "成功" 组合
successful_indices = find(final_errors_flat < convergence_tol);
num_successful = length(successful_indices);

fprintf('   在求解成功的案例中，有 %d 个组合的最终误差小于设定的容忍度 (%.1e)。\n', ...
    num_successful, convergence_tol);

fprintf('\n--- 收敛质量达标的参数组合列表 (误差 < %.1e) ---\n', convergence_tol);

if num_successful == 0
    fprintf('   没有任何参数组合的误差低于设定的容忍度。\n');
else
    % 准备表格数据
    results_table = table(...
        beta_flat(successful_indices), ...
        sigma_flat(successful_indices), ...
        phi_flat(successful_indices), ...
        theta_flat(successful_indices), ...
        gA_flat(successful_indices), ...
        n_flat(successful_indices), ...
        final_errors_flat(successful_indices), ...
        'VariableNames', {'beta', 'sigma', 'phi_bequest', 'theta', 'g_A_ss', 'n_ss', 'MaxError'});
        
    % 按误差从小到大排序
    results_table = sortrows(results_table, 'MaxError');
    
    % 显示表格
    num_to_display = min(num_successful, max_results_to_display);
    fprintf('   (显示前 %d 个误差最小的结果)\n\n', num_to_display);
    disp(results_table(1:num_to_display, :));

    if num_successful > max_results_to_display
        fprintf('\n   ... 以及另外 %d 个成功的组合未显示。\n', num_successful - max_results_to_display);
    end
    
    % (可选) 保存成功的表格结果到CSV文件
    try
        writetable(results_table, 'SS/successful_convergence_combinations.csv');
        fprintf('\n   所有 %d 个成功组合的详细列表已保存至: SS/successful_convergence_combinations.csv\n', num_successful);
    catch ME
        warning('无法保存CSV文件，可能SS文件夹不存在或无写入权限。错误: %s', ME.message);
    end
end


% (可选) 保存完整的原始误差矩阵供后续详细分析
save('SS/convergence_test_results.mat', 'final_errors_flat', 'beta_vec', 'sigma_vec', 'phi_bequest_vec', 'theta_vec', 'g_A_ss_vec', 'n_ss_vec');
fprintf('\n   完整的原始误差数据已保存至 SS/convergence_test_results.mat\n');

fprintf('\n--- 稳健性测试脚本执行完毕 ---\n');