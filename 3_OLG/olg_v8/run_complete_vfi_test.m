% --- 网格搜索 vs fmincon 连续优化 VFI 比较测试脚本 ---

clc; % 清除命令行窗口
clear; % 清除工作区变量
close all; % 关闭所有图形窗口

fprintf('=== 网格搜索 vs fmincon 连续优化 VFI 多参数比较测试 ===\n');

% 1. 获取模型参数
fprintf('\n--- 1. 加载模型参数 ---\n');
try
    cS = main_olg_v8_utils.ParameterValues_HuggettStyle();
    fprintf('模型参数已成功加载。\n');
    fprintf('    基础状态变量参数：nk=%d, nkpps=%d, nw=%d\n', cS.nk, cS.nkpps, cS.nw);
    fprintf('    基础PPS离散选择点数：%d\n', cS.n_pps_choice_grid_points);
    
catch ME
    fprintf('错误：加载模型参数失败。\n');
    fprintf('错误信息: %s\n', ME.message);
    return; % 参数加载失败则退出
end

% 2. 启动并行环境
fprintf('\n--- 2. 启动并行计算环境 ---\n');
try
    pool = gcp('nocreate');
    if isempty(pool)
        fprintf('启动并行池...\n');
        parpool();
        pool = gcp('nocreate');
    end
    if ~isempty(pool)
        fprintf('并行池已激活，工作进程数: %d\n', pool.NumWorkers);
    else
        fprintf('警告：并行池未能启动，将使用串行计算。\n');
    end
catch ME_pool
    fprintf('警告：并行池操作失败: %s\n', ME_pool.message);
end

% 3. 定义测试参数组合
fprintf('\n--- 3. 定义测试参数组合 ---\n');

% 测试案例定义
test_cases = struct();

% 案例1: 年轻工作者（低PPS缴费）
test_cases(1).name = '年轻工作者';
test_cases(1).age_idx = 3;
test_cases(1).nk = 15;
test_cases(1).nkpps = 15;
test_cases(1).nw = 3;
test_cases(1).n_pps_choice_grid_points = 8;

% 案例2: 中年工作者（中等PPS缴费）
test_cases(2).name = '中年工作者';
test_cases(2).age_idx = 7;
test_cases(2).nk = 20;
test_cases(2).nkpps = 20;
test_cases(2).nw = 5;
test_cases(2).n_pps_choice_grid_points = 12;

% 案例3: 临近退休者（高PPS缴费）
test_cases(3).name = '临近退休者';
test_cases(3).age_idx = 11;
test_cases(3).nk = 25;
test_cases(3).nkpps = 25;
test_cases(3).nw = 5;
test_cases(3).n_pps_choice_grid_points = 15;

% 案例4: 小规模测试（快速验证）
test_cases(4).name = '小规模测试';
test_cases(4).age_idx = 5;
test_cases(4).nk = 10;
test_cases(4).nkpps = 10;
test_cases(4).nw = 3;
test_cases(4).n_pps_choice_grid_points = 6;

% 案例5: 大规模测试（性能压力测试）
test_cases(5).name = '大规模测试';
test_cases(5).age_idx = 8;
test_cases(5).nk = 30;
test_cases(5).nkpps = 30;
test_cases(5).nw = 7;
test_cases(5).n_pps_choice_grid_points = 20;

fprintf('定义了 %d 个测试案例\n', length(test_cases));

% 4. 执行多案例比较测试
fprintf('\n--- 4. 执行多案例比较测试 ---\n');

% 存储所有结果 - 使用cell数组方式初始化
all_results = [];

for case_idx = 1:length(test_cases)
    tc = test_cases(case_idx);
    fprintf('\n=== 案例 %d: %s ===\n', case_idx, tc.name);
    fprintf('参数: 年龄组=%d, nk=%d, nkpps=%d, nw=%d, PPS格点=%d\n', ...
        tc.age_idx, tc.nk, tc.nkpps, tc.nw, tc.n_pps_choice_grid_points);
    
    % 设置测试参数
    cS_test = cS;
    cS_test.nk = tc.nk;
    cS_test.nkpps = tc.nkpps;
    cS_test.nw = tc.nw;
    cS_test.n_pps_choice_grid_points = tc.n_pps_choice_grid_points;
    
    % 重新生成网格
    cS_test = main_olg_v8_utils.generateGrids(cS_test);
    
    % 计算实际年龄
    actual_age = cS_test.age1_orig + tc.age_idx - 1;
    fprintf('实际年龄: %d岁\n', actual_age);
    
    try
        % 运行单案例比较
        case_result = run_single_case_comparison(tc, cS_test);
        
        % 第一次赋值时初始化结构体数组
        if isempty(all_results)
            all_results = case_result;
        else
            all_results(end+1) = case_result;
        end
        
        fprintf('案例 %d 完成: 网格=%.2fs, fmincon=%.2fs, 速度比=%.1fx\n', ...
            case_idx, case_result.grid_time, case_result.fmincon_time, ...
            case_result.speed_ratio);
        
    catch ME_case
        fprintf('案例 %d 失败: %s\n', case_idx, ME_case.message);
        % 创建失败案例的结构体
        failed_case.success = false;
        failed_case.case_name = tc.name;
        failed_case.error = ME_case.message;
        failed_case.grid_time = NaN;
        failed_case.fmincon_time = NaN;
        failed_case.speed_ratio = NaN;
        failed_case.state_space_size = tc.nk * tc.nkpps * tc.nw;
        failed_case.pps_grid_points = tc.n_pps_choice_grid_points;
        failed_case.c_diff_max = NaN;
        failed_case.c_diff_mean = NaN;
        failed_case.k_diff_max = NaN;
        failed_case.k_diff_mean = NaN;
        failed_case.cpps_diff_max = NaN;
        failed_case.cpps_diff_mean = NaN;
        failed_case.val_diff_max = NaN;
        failed_case.val_diff_mean = NaN;
        failed_case.grid_smoothness = NaN;
        failed_case.fmincon_smoothness = NaN;
        failed_case.grid_stats = struct();
        failed_case.fmincon_stats = struct();
        
        % k决策变量相关字段
        failed_case.k_grid_smoothness = NaN;
        failed_case.k_fmincon_smoothness = NaN;
        failed_case.k_grid_stats = struct();
        failed_case.k_fmincon_stats = struct();
        
        % 值函数优劣比较字段
        failed_case.val_grid_mean = NaN;
        failed_case.val_fmincon_mean = NaN;
        failed_case.val_superiority = NaN;
        failed_case.fmincon_better_ratio = NaN;
        failed_case.max_val_improvement = NaN;
        failed_case.min_val_improvement = NaN;
        failed_case.significant_improvement_count = NaN;
        failed_case.significant_worsening_count = NaN;
        failed_case.total_states = NaN;
        
        % 第一次赋值时初始化结构体数组
        if isempty(all_results)
            all_results = failed_case;
        else
            all_results(end+1) = failed_case;
        end
    end
end

% 5. 综合结果分析
fprintf('\n--- 5. 综合结果分析 ---\n');
analyze_all_results(all_results, test_cases);

fprintf('\n=== 测试完成 ===\n');

% ===== 本地函数定义 =====

function case_result = run_single_case_comparison(tc, cS_test)
    % 运行单个案例的比较测试
    
    % 构造测试用的输入参数
    vPrime_dummy = zeros(cS_test.nk, cS_test.nkpps, cS_test.nw);
    R_k_net_factor = 1.05;
    w_gross = 1.0;
    TR_total = 0;
    b_age_val = 0;
    
    % 构造家庭决策参数
    paramS_test.leTrProbM = eye(cS_test.nw);
    paramS_test.tau_l = cS_test.tau_l_init_guess;
    paramS_test.theta_payg_actual_for_hh = 0.1;
    paramS_test.pps_tax_deferral_active = false;
    
    % 构造劳动效率网格
    epsilon_grid = linspace(0.5, 1.5, cS_test.nw);
    
    % 测试网格搜索方法
    fprintf('  运行网格搜索方法...\n');
    tic;
    [cPol_grid, kPol_grid, cPpsPol_grid, val_grid] = ...
        main_olg_v8_utils.HHSolutionByAge_VFI_Huggett_v8(...
        tc.age_idx, vPrime_dummy, R_k_net_factor, w_gross, ...
        TR_total, b_age_val, paramS_test, cS_test, epsilon_grid);
    grid_time = toc;
    
    % 测试fmincon连续优化方法
    fprintf('  运行fmincon连续优化方法...\n');
    tic;
    [cPol_fmincon, kPol_fmincon, cPpsPol_fmincon, val_fmincon] = ...
        main_olg_v8_utils.HHSolutionByAge_VFI_Huggett_v8_fmincon(...
        tc.age_idx, vPrime_dummy, R_k_net_factor, w_gross, ...
        TR_total, b_age_val, paramS_test, cS_test, epsilon_grid);
    fmincon_time = toc;
    
    % 计算差异
    c_diff_max = max(abs(cPol_grid(:) - cPol_fmincon(:)));
    c_diff_mean = mean(abs(cPol_grid(:) - cPol_fmincon(:)));
    k_diff_max = max(abs(kPol_grid(:) - kPol_fmincon(:)));
    k_diff_mean = mean(abs(kPol_grid(:) - kPol_fmincon(:)));
    cpps_diff_max = max(abs(cPpsPol_grid(:) - cPpsPol_fmincon(:)));
    cpps_diff_mean = mean(abs(cPpsPol_grid(:) - cPpsPol_fmincon(:)));
    val_diff_max = max(abs(val_grid(:) - val_fmincon(:)));
    val_diff_mean = mean(abs(val_grid(:) - val_fmincon(:)));
    
    % 值函数优劣比较分析
    val_grid_mean = mean(val_grid(:));
    val_fmincon_mean = mean(val_fmincon(:));
    val_superiority = val_fmincon_mean - val_grid_mean;  % 正值表示fmincon更优
    
    % 计算在多少个状态点上fmincon更优
    fmincon_better_count = sum(val_fmincon(:) > val_grid(:));
    total_states = numel(val_grid);
    fmincon_better_ratio = fmincon_better_count / total_states;
    
    % 计算值函数提升幅度的统计
    val_improvement = val_fmincon(:) - val_grid(:);
    max_val_improvement = max(val_improvement);
    min_val_improvement = min(val_improvement);
    significant_improvement_count = sum(val_improvement > 0.001);  % 显著改善的状态数
    significant_worsening_count = sum(val_improvement < -0.001);   % 显著恶化的状态数
    
    % 分析平滑性
    [grid_smoothness, grid_stats] = analyze_smoothness_metrics(cPpsPol_grid);
    [fmincon_smoothness, fmincon_stats] = analyze_smoothness_metrics(cPpsPol_fmincon);
    
    % 分析k决策变量的平滑性
    [k_grid_smoothness, k_grid_stats] = analyze_smoothness_metrics(kPol_grid);
    [k_fmincon_smoothness, k_fmincon_stats] = analyze_smoothness_metrics(kPol_fmincon);
    
    % 存储结果
    case_result.success = true;
    case_result.case_name = tc.name;
    case_result.grid_time = grid_time;
    case_result.fmincon_time = fmincon_time;
    case_result.speed_ratio = fmincon_time / grid_time;
    case_result.state_space_size = tc.nk * tc.nkpps * tc.nw;
    case_result.pps_grid_points = tc.n_pps_choice_grid_points;
    
    % 差异指标
    case_result.c_diff_max = c_diff_max;
    case_result.c_diff_mean = c_diff_mean;
    case_result.k_diff_max = k_diff_max;
    case_result.k_diff_mean = k_diff_mean;
    case_result.cpps_diff_max = cpps_diff_max;
    case_result.cpps_diff_mean = cpps_diff_mean;
    case_result.val_diff_max = val_diff_max;
    case_result.val_diff_mean = val_diff_mean;
    
    % 值函数优劣比较指标
    case_result.val_grid_mean = val_grid_mean;
    case_result.val_fmincon_mean = val_fmincon_mean;
    case_result.val_superiority = val_superiority;
    case_result.fmincon_better_ratio = fmincon_better_ratio;
    case_result.max_val_improvement = max_val_improvement;
    case_result.min_val_improvement = min_val_improvement;
    case_result.significant_improvement_count = significant_improvement_count;
    case_result.significant_worsening_count = significant_worsening_count;
    case_result.total_states = total_states;
    
    % 平滑性指标
    case_result.grid_smoothness = grid_smoothness;
    case_result.fmincon_smoothness = fmincon_smoothness;
    case_result.grid_stats = grid_stats;
    case_result.fmincon_stats = fmincon_stats;
    
    % k决策变量平滑性指标
    case_result.k_grid_smoothness = k_grid_smoothness;
    case_result.k_fmincon_smoothness = k_fmincon_smoothness;
    case_result.k_grid_stats = k_grid_stats;
    case_result.k_fmincon_stats = k_fmincon_stats;
    
    fprintf('  完成: 网格=%.2fs, fmincon=%.2fs, K最大差异=%.6f, CPPS最大差异=%.6f\n', ...
        grid_time, fmincon_time, k_diff_max, cpps_diff_max);
    fprintf('  值函数比较: fmincon平均值=%.6f, 网格平均值=%.6f, 优势=%.6f\n', ...
        val_fmincon_mean, val_grid_mean, val_superiority);
    fprintf('  fmincon更优的状态比例: %.1f%% (%d/%d)\n', ...
        fmincon_better_ratio*100, fmincon_better_count, total_states);
end

function [smoothness_score, stats] = analyze_smoothness_metrics(cpps_policy)
    % 计算平滑性指标
    
    mean_cpps_by_kpps = squeeze(mean(cpps_policy, [1,3]));
    
    stats.min_val = min(mean_cpps_by_kpps);
    stats.max_val = max(mean_cpps_by_kpps);
    stats.mean_val = mean(mean_cpps_by_kpps);
    
    if length(mean_cpps_by_kpps) > 1
        differences = diff(mean_cpps_by_kpps);
        stats.max_jump = max(abs(differences));
        stats.std_diff = std(differences);
        
        % 计算平滑性得分（越小越平滑）
        smoothness_score = stats.max_jump;
    else
        stats.max_jump = 0;
        stats.std_diff = 0;
        smoothness_score = 0;
    end
end

function analyze_all_results(all_results, test_cases)
    % 分析所有案例的综合结果
    
    successful_cases = [all_results.success];
    if sum(successful_cases) == 0
        fprintf('没有成功的测试案例。\n');
        return;
    end
    
    results = all_results(successful_cases);
    cases = test_cases(successful_cases);
    
    fprintf('\n=== 性能总结 ===\n');
    fprintf('%-15s | %-8s | %-8s | %-8s | %-10s | %-10s | %-10s | %-8s\n', ...
        '案例名称', '网格(s)', 'fmincon(s)', '速度比', 'K最大差异', 'PPS最大差异', '值函数优势', 'PPS格点');
    fprintf('%s\n', repmat('-', 110, 1));
    
    for i = 1:length(results)
        fprintf('%-15s | %8.2f | %8.2f | %8.1fx | %10.6f | %10.6f | %10.6f | %8d\n', ...
            results(i).case_name, results(i).grid_time, results(i).fmincon_time, ...
            results(i).speed_ratio, results(i).k_diff_max, results(i).cpps_diff_max, ...
            results(i).val_superiority, results(i).pps_grid_points);
    end
    
    % 总体统计
    avg_speed_ratio = mean([results.speed_ratio]);
    min_speed_ratio = min([results.speed_ratio]);
    max_speed_ratio = max([results.speed_ratio]);
    
    fprintf('\n总体性能统计:\n');
    fprintf('  平均速度比 (fmincon/网格): %.2fx\n', avg_speed_ratio);
    fprintf('  最快情况: %.2fx\n', min_speed_ratio);
    fprintf('  最慢情况: %.2fx\n', max_speed_ratio);
    
    % 差异分析
    fprintf('\n=== 解的差异分析 ===\n');
    
    % K决策变量（储蓄）差异分析
    avg_k_diff = mean([results.k_diff_max]);
    max_k_diff = max([results.k_diff_max]);
    avg_k_diff_mean = mean([results.k_diff_mean]);
    
    fprintf('储蓄决策(K)差异:\n');
    fprintf('  平均最大差异: %.6f\n', avg_k_diff);
    fprintf('  所有案例最大差异: %.6f\n', max_k_diff);
    fprintf('  平均均值差异: %.6f\n', avg_k_diff_mean);
    
    if max_k_diff > 0.1
        fprintf('  → K决策存在显著差异，两种方法储蓄策略明显不同\n');
    elseif max_k_diff > 0.01
        fprintf('  → K决策存在中等差异，储蓄策略有一定差别\n');
    else
        fprintf('  → K决策差异很小，储蓄策略基本一致\n');
    end
    
    % PPS缴费策略差异分析
    avg_cpps_diff = mean([results.cpps_diff_max]);
    max_cpps_diff = max([results.cpps_diff_max]);
    avg_cpps_diff_mean = mean([results.cpps_diff_mean]);
    
    fprintf('\nPPS缴费策略差异:\n');
    fprintf('  平均最大差异: %.6f\n', avg_cpps_diff);
    fprintf('  所有案例最大差异: %.6f\n', max_cpps_diff);
    fprintf('  平均均值差异: %.6f\n', avg_cpps_diff_mean);
    
    if max_cpps_diff > 0.01
        fprintf('  → PPS缴费存在显著差异，fmincon找到了不同的解\n');
    elseif max_cpps_diff > 0.001
        fprintf('  → PPS缴费存在中等差异，可能由数值精度导致\n');
    else
        fprintf('  → PPS缴费差异很小，两种方法基本一致\n');
    end
    
    % 消费决策差异分析
    avg_c_diff = mean([results.c_diff_max]);
    max_c_diff = max([results.c_diff_max]);
    avg_c_diff_mean = mean([results.c_diff_mean]);
    
    fprintf('\n消费决策(C)差异:\n');
    fprintf('  平均最大差异: %.6f\n', avg_c_diff);
    fprintf('  所有案例最大差异: %.6f\n', max_c_diff);
    fprintf('  平均均值差异: %.6f\n', avg_c_diff_mean);
    
    % 价值函数差异分析
    avg_val_diff = mean([results.val_diff_max]);
    max_val_diff = max([results.val_diff_max]);
    
    fprintf('\n价值函数差异:\n');
    fprintf('  平均最大差异: %.6f\n', avg_val_diff);
    fprintf('  所有案例最大差异: %.6f\n', max_val_diff);
    
    % 值函数优劣分析
    fprintf('\n=== 值函数优劣分析 ===\n');
    
    val_superiorities = [results.val_superiority];
    fmincon_better_ratios = [results.fmincon_better_ratio];
    
    avg_val_superiority = mean(val_superiorities);
    avg_fmincon_better_ratio = mean(fmincon_better_ratios);
    
    fprintf('整体优劣比较:\n');
    fprintf('  fmincon平均优势: %.6f\n', avg_val_superiority);
    fprintf('  fmincon更优状态平均比例: %.1f%%\n', avg_fmincon_better_ratio*100);
    
    % 案例级别的优劣统计
    fmincon_wins = sum(val_superiorities > 0);
    grid_wins = sum(val_superiorities < 0);
    ties = sum(abs(val_superiorities) < 1e-6);
    
    fprintf('\n案例级别胜负统计:\n');
    fprintf('  fmincon更优案例数: %d/%d\n', fmincon_wins, length(results));
    fprintf('  网格搜索更优案例数: %d/%d\n', grid_wins, length(results));
    fprintf('  基本相等案例数: %d/%d\n', ties, length(results));
    
    % 详细的改善统计
    total_improvements = sum([results.significant_improvement_count]);
    total_worsenings = sum([results.significant_worsening_count]);
    total_all_states = sum([results.total_states]);
    
    fprintf('\n详细改善统计 (阈值=0.001):\n');
    fprintf('  显著改善的状态总数: %d/%d (%.1f%%)\n', ...
        total_improvements, total_all_states, total_improvements/total_all_states*100);
    fprintf('  显著恶化的状态总数: %d/%d (%.1f%%)\n', ...
        total_worsenings, total_all_states, total_worsenings/total_all_states*100);
    
    % 最大改善和恶化
    max_improvement = max([results.max_val_improvement]);
    min_improvement = min([results.min_val_improvement]);
    
    fprintf('\n极值改善情况:\n');
    fprintf('  最大值函数提升: %.6f\n', max_improvement);
    fprintf('  最大值函数下降: %.6f\n', abs(min_improvement));
    
    % 优劣结论
    fprintf('\n值函数优劣结论:\n');
    if avg_val_superiority > 0.001
        fprintf('  → fmincon方法显著优于网格搜索，找到了更优的解\n');
    elseif avg_val_superiority > 0
        fprintf('  → fmincon方法略优于网格搜索\n');
    elseif avg_val_superiority > -0.001
        fprintf('  → 两种方法基本相当\n');
    else
        fprintf('  → 网格搜索方法优于fmincon\n');
    end
    
    if avg_fmincon_better_ratio > 0.6
        fprintf('  → fmincon在大多数状态下表现更好\n');
    elseif avg_fmincon_better_ratio > 0.4
        fprintf('  → 两种方法在不同状态下各有优势\n');
    else
        fprintf('  → 网格搜索在大多数状态下表现更好\n');
    end
    
    % 平滑性分析
    fprintf('\n=== 平滑性比较 ===\n');
    
    % PPS缴费策略平滑性
    grid_smoothness = [results.grid_smoothness];
    fmincon_smoothness = [results.fmincon_smoothness];
    
    fprintf('PPS缴费策略平滑性指标 (越小越平滑):\n');
    fprintf('  网格搜索平均: %.6f\n', mean(grid_smoothness));
    fprintf('  fmincon平均: %.6f\n', mean(fmincon_smoothness));
    
    if mean(fmincon_smoothness) < mean(grid_smoothness)
        fprintf('  → fmincon的PPS缴费策略更平滑\n');
    else
        fprintf('  → 网格搜索的PPS缴费策略更平滑\n');
    end
    
    % K决策变量（储蓄）平滑性
    k_grid_smoothness = [results.k_grid_smoothness];
    k_fmincon_smoothness = [results.k_fmincon_smoothness];
    
    fprintf('\n储蓄决策(K)平滑性指标 (越小越平滑):\n');
    fprintf('  网格搜索平均: %.6f\n', mean(k_grid_smoothness));
    fprintf('  fmincon平均: %.6f\n', mean(k_fmincon_smoothness));
    
    if mean(k_fmincon_smoothness) < mean(k_grid_smoothness)
        fprintf('  → fmincon的储蓄策略更平滑\n');
    else
        fprintf('  → 网格搜索的储蓄策略更平滑\n');
    end
    
    % 整体平滑性比较
    fprintf('\n整体平滑性评估:\n');
    pps_smooth_advantage = mean(grid_smoothness) - mean(fmincon_smoothness);
    k_smooth_advantage = mean(k_grid_smoothness) - mean(k_fmincon_smoothness);
    
    fprintf('  PPS策略：fmincon优势 = %.6f\n', pps_smooth_advantage);
    fprintf('  储蓄策略：fmincon优势 = %.6f\n', k_smooth_advantage);
    
    if pps_smooth_advantage > 0 && k_smooth_advantage > 0
        fprintf('  → fmincon在两个决策变量上都更平滑\n');
    elseif pps_smooth_advantage > 0 || k_smooth_advantage > 0
        fprintf('  → fmincon在某些决策变量上更平滑\n');
    else
        fprintf('  → 网格搜索方法整体上更平滑\n');
    end
    
    % 规模效应分析
    fprintf('\n=== 规模效应分析 ===\n');
    state_sizes = [results.state_space_size];
    speed_ratios = [results.speed_ratio];
    
    [sorted_sizes, sort_idx] = sort(state_sizes);
    sorted_ratios = speed_ratios(sort_idx);
    
    fprintf('按状态空间大小排序的速度比:\n');
    for i = 1:length(sorted_sizes)
        fprintf('  状态空间 %d: 速度比 %.2fx\n', sorted_sizes(i), sorted_ratios(i));
    end
    
    % 建议
    fprintf('\n=== 使用建议 ===\n');
    
    % 获取需要的变量
    max_k_diff = max([results.k_diff_max]);
    
    % 基于计算效率的建议
    if avg_speed_ratio < 1.5
        fprintf('• 计算效率: fmincon方法计算时间可接受，建议优先使用\n');
    elseif avg_speed_ratio < 3
        fprintf('• 计算效率: fmincon方法计算时间适中，可在需要高精度时使用\n');
    else
        fprintf('• 计算效率: fmincon方法计算时间较长，建议仅在对精度要求极高时使用\n');
    end
    
    % 基于解的质量的建议
    if avg_val_superiority > 0.001
        fprintf('• 解的质量: fmincon显著优于网格搜索，强烈建议使用fmincon\n');
    elseif avg_val_superiority > 0
        fprintf('• 解的质量: fmincon略优于网格搜索，建议使用fmincon\n');
    elseif avg_val_superiority > -0.001
        fprintf('• 解的质量: 两种方法质量相当，可根据其他因素选择\n');
    else
        fprintf('• 解的质量: 网格搜索优于fmincon，建议使用网格搜索\n');
    end
    
    % 基于差异程度的建议
    if max_cpps_diff > 0.01 || max_k_diff > 0.1
        fprintf('• 解的差异: 两种方法存在显著差异，建议根据具体需求选择方法\n');
    else
        fprintf('• 解的差异: 两种方法结果基本一致，可优先考虑计算效率\n');
    end
    
    % 综合建议
    fprintf('\n综合建议:\n');
    if avg_val_superiority > 0.001 && avg_speed_ratio < 3
        fprintf('  → 推荐使用fmincon方法：解的质量更优且计算时间可接受\n');
    elseif avg_val_superiority > 0 && avg_speed_ratio < 5
        fprintf('  → 建议使用fmincon方法：解的质量略优，可承受额外计算成本\n');
    elseif avg_val_superiority < -0.001
        fprintf('  → 推荐使用网格搜索方法：解的质量更优\n');
    elseif avg_speed_ratio > 5
        fprintf('  → 推荐使用网格搜索方法：计算效率更高且解的质量相当\n');
    else
        fprintf('  → 两种方法各有优劣，建议根据具体应用场景选择\n');
    end
end 