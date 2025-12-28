% =====================================================================
% == 比较VFI求解器 (v3: 聚焦 Hybrid vs Vectorized) ==
% =====================================================================
%
% 目标:
% 在固定的宏观环境下，精确对比两种核心VFI求解器的性能：
% 1. 'hybrid' (20x20): 中等精度，结合了网格搜索和连续优化，精度较高。
% 2. 'vectorized_grid' (80x80): 高精度，全矩阵化，速度快。
%
% 评估维度:
% - 策略质量 (平均终身效用)
% - 计算效率 (求解时间)
% - 生成的策略函数 (平均消费路径)
%
% =====================================================================

%% 1. 初始化环境和参数
clear; clc; close all;
addpath(pwd); % 确保 utils 文件在 MATLAB 路径中

fprintf('--- 1. 设置测试参数 ---\n');

% [核心修改] 定义要测试的求解器和网格配置
% 这个结构与Python的vfi_grid_configs列表一一对应
% configurations_to_test = { ...
%     struct('label', 'VFI_Med_hybrid (20x20)', 'nk', 20, 'nkpps', 20, 'nkprime', 20, 'npps', 20, 'solver', 'hybrid'), ...
%     struct('label', 'VFI_vectorized (80x80)', 'nk', 80, 'nkpps', 80, 'nkprime', 80, 'npps', 80, 'solver', 'vectorized_grid') ...
% };
configurations_to_test = { ...
    struct('label', 'VFI_hybrid (20x20)', 'nk', 20, 'nkpps', 20, 'nkprime', 20, 'npps', 20, 'solver', 'hybrid'), ...
    struct('label', 'VFI_vectorized (50x50)', 'nk', 50, 'nkpps', 50, 'nkprime', 50, 'npps', 50, 'solver', 'vectorized_grid'), ...
};
% configurations_to_test = { ...
%     struct('label', 'VFI_grid (5x5)', 'nk',5, 'nkpps', 5, 'nkprime',5, 'npps', 5, 'solver', 'grid'), ...
%     struct('label', 'VFI_hybrid (5x5)', 'nk',5, 'nkpps', 5, 'nkprime',5, 'npps', 5, 'solver', 'hybrid'), ...
%     struct('label', 'VFI_vectorized (50x50)', 'nk',50, 'nkpps', 50, 'nkprime',50, 'npps', 50, 'solver', 'vectorized_grid'), ...
% 
% };

% 模拟参数
n_sim = 500; % 与Python脚本对齐
random_seed = 42;
rng(random_seed); % 设置随机种子

%% 2. 准备固定的宏观环境和效率冲击路径
fprintf('--- 2. 准备宏观环境和效率路径 ---\n');

% 从Python脚本生成的.mat文件中加载统一的效率冲击路径
shock_path_filename = 'eIdxM_group_global.mat';
if exist(shock_path_filename, 'file')
    load(shock_path_filename, 'eIdxM_group_global');
    fprintf('✅ 从 %s 加载了全局效率路径 (size: %d x %d)。\n', ...
        shock_path_filename, size(eIdxM_group_global, 1), size(eIdxM_group_global, 2));
else
    % 如果没有找到文件，生成默认的效率冲击路径
    fprintf('⚠️ 找不到效率冲击路径文件 %s，将生成默认路径。\n', shock_path_filename);
    % 使用一个临时的参数结构来生成效率路径
    temp_cS = struct();
    temp_cS.nSim = n_sim;
    temp_cS.aD_new = 16; % 默认年龄组数
    temp_cS.nw = 5; % 默认效率状态数
    [temp_leLogGridV, temp_leTrProbM, temp_leProb1V] = main_olg_v10_utils.EarningProcess_olgm(temp_cS);
    eIdxM_group_global = main_olg_v10_utils.LaborEndowSimulation_olgm_AgeGroup(temp_cS, ...
        struct('leProb1V', temp_leProb1V, 'leTrProbM', temp_leTrProbM));
end

% 定义固定的宏观环境 (确保与Python的M_fixed完全一致)
M_FIXED = struct(...
    'R_k_net_factor', 1.03, ...
    'w_gross', 1.8, ...
    'TR_total', 0.0, ...
    'b_payg_avg_retiree', 0.4, ...
    'tau_l', 0.15, ...
    'theta_payg_actual', 0.10 ... 
);
fprintf('✅ 宏观环境已设置 (w=%.2f, TR=%.2f, theta=%.2f)。\n', ...
    M_FIXED.w_gross, M_FIXED.TR_total, M_FIXED.theta_payg_actual);

%% 3. 运行分析循环
all_results = {};

for i_config = 1:length(configurations_to_test)
    config = configurations_to_test{i_config};
    
    fprintf('\n============================================================\n');
    fprintf('运行配置: %s\n', config.label);
    fprintf('nk=%d, nkpps=%d, nkprime=%d, npps=%d, solver=''%s''\n', ...
        config.nk, config.nkpps, config.nkprime, config.npps, config.solver);
    fprintf('============================================================\n');
    
    % a. 为当前配置生成参数
    cS = main_olg_v10_utils.ParameterValues_HuggettStyle();
    cS.nk = config.nk;
    cS.nkpps = config.nkpps;
    cS.nkprime = config.nkprime;
    cS.npps = config.npps;
    cS.nSim = n_sim;
    
    % [新增] 设置v10新增的必要参数
    % --- V10 固定的政策参数 ---
    cS.rho_prime_payg = 0.2; 
    cS.theta_payg = 0.2; 
    cS.tau_l = 0.1;
    B_p_Y_ratio_target = 0.05; % [新] 目标：初始稳态养老金基金占GDP的5%
    
    cS = main_olg_v10_utils.generateGrids(cS);
    
    paramS = struct();
    % [对齐] 确保随机过程参数与Python一致
    cS.lePersistence = 0.90; 
    cS.leShockStd = 0.15;
    [paramS.leLogGridV, paramS.leTrProbM, paramS.leProb1V] = main_olg_v10_utils.EarningProcess_olgm(cS);
    paramS.leGridV = exp(paramS.leLogGridV);
    
    % [对齐] 确保税收参数与Python一致
    paramS.tau_l = M_FIXED.tau_l;
    paramS.theta_payg_actual_for_hh = M_FIXED.theta_payg_actual;
    paramS.pps_tax_deferral_active = cS.pps_active;

    bV_payg_vfi = zeros(1, cS.aD_new);
    if cS.aR_new < cS.aD_new
        bV_payg_vfi((cS.aR_new+1):end) = M_FIXED.b_payg_avg_retiree;
    end

    % b. 求解VFI (使用已修复的函数)
    tic;
    [cPolM, kPolM, cPpsPolM, ~] = main_olg_v10_utils.HHSolution_VFI_Huggett(...
        M_FIXED.R_k_net_factor, M_FIXED.w_gross, M_FIXED.TR_total, ...
        bV_payg_vfi, paramS, cS, config.solver);
    time_solve = toc;
    fprintf('VFI求解耗时: %.2f 秒。\n', time_solve);

    % c. 模拟生命周期路径
    tic;
    [k_path, ~, c_path, ~] = main_olg_v10_utils.HHSimulation_olgm(...
        kPolM, cPpsPolM, cPolM, eIdxM_group_global, ...
        M_FIXED.R_k_net_factor, M_FIXED.w_gross, M_FIXED.TR_total, ...
        bV_payg_vfi, paramS, cS);
    time_sim = toc;
    fprintf('模拟耗时: %.2f 秒。\n', time_sim);

    % d. 计算终身效用
    utility_vfi = calculate_lifetime_utility_local(c_path, cS, true);
    
    % e. 存储结果
    result_data = config; % 复制配置信息
    result_data.mean_utility = mean(utility_vfi);
    result_data.std_utility = std(utility_vfi);
    result_data.c_path = c_path;
    result_data.k_path = k_path;
    result_data.solve_time = time_solve;
    result_data.sim_time = time_sim;
    all_results{end+1} = result_data;
    
    fprintf('📈 结果: 平均效用 = %.4f (标准差 = %.4f)\n', result_data.mean_utility, result_data.std_utility);
end

%% 4. 可视化结果
fprintf('\n--- 4. 可视化结果 ---\n');
plot_focused_comparison(all_results);


%% 辅助函数 (Local Functions)

function lifetime_utility = calculate_lifetime_utility_local(c_path, cS, use_survival_prob)
    % 计算终身效用，适配v10版本
    [n_sim, aD] = size(c_path);
    lifetime_utility = zeros(n_sim, 1);
    beta = cS.beta;
    s_transitionV = cS.s_1yr_transitionV;

    for i_sim = 1:n_sim
        utility_sum = 0.0;
        cumulative_discount = 1.0;
        for a_group = 1:aD
            [~, u] = main_olg_v10_utils.CES_utility(c_path(i_sim, a_group), cS.sigma, cS);
            utility_sum = utility_sum + cumulative_discount * u;
            if a_group < aD
                survival_factor = 1.0;
                if use_survival_prob
                    survival_factor = s_transitionV(a_group);
                end
                cumulative_discount = cumulative_discount * (beta * survival_factor);
            end
        end
        lifetime_utility(i_sim) = utility_sum;
    end
end

function plot_focused_comparison(results_cell)
    % [核心修改] 增强的绘图函数，用于聚焦对比
    
    figure('Name', 'VFI Solver Focused Comparison', 'Position', [100, 100, 1200, 500]);
    
    labels = cellfun(@(x) x.label, results_cell, 'UniformOutput', false);
    utilities = cellfun(@(x) x.mean_utility, results_cell);
    solve_times = cellfun(@(x) x.solve_time, results_cell);
    
    % --- 子图1: 平均效用对比 ---
    subplot(1, 2, 1);
    b1 = bar(utilities);
    set(gca, 'xticklabel', labels, 'XTickLabelRotation', 20);
    ylabel('平均终身效用');
    title('策略质量 (平均效用)');
    grid on;
    % 添加数值标签
    xtips1 = 1:length(labels);
    ytips1 = b1.YData;
    text(xtips1, ytips1, sprintfc('%.4f', ytips1), 'HorizontalAlignment','center', 'VerticalAlignment','bottom');
    ylim([min(utilities)-0.5, max(utilities)+0.5]); % 调整Y轴范围

    % --- 子图2: 平均消费路径 & 求解时间对比 ---
    ax2 = subplot(1, 2, 2);
    hold(ax2, 'on');
    
    colors = lines(length(results_cell));
    bar_data = [];
    bar_labels = {};
    
    for i = 1:length(results_cell)
        res = results_cell{i};
        mean_c_path = mean(res.c_path, 1);
        plot(ax2, mean_c_path, 'o-', 'LineWidth', 1.5, 'Color', colors(i,:), 'DisplayName', res.label);
        
        bar_data(i) = res.solve_time;
        bar_labels{i} = sprintf('%s\n(%.2f s)', res.label, res.solve_time);
    end
    
    hold(ax2, 'off');
    xlabel(ax2, '年龄组');
    ylabel(ax2, '平均消费');
    title(ax2, '平均消费路径');
    legend(ax2, 'show', 'Location', 'best', 'Interpreter', 'none');
    grid(ax2, 'on');
    
    % 在右侧使用第二个Y轴显示求解时间
    yyaxis(ax2, 'right');
    bar(ax2, bar_data, 0.4, 'FaceAlpha', 0.2);
    ylabel(ax2, '求解时间 (秒)');
    set(ax2, 'YColor', 'k');
    set(ax2, 'xtick', 1:length(labels), 'xticklabel', []); % 隐藏右侧x轴标签
    ylim(ax2, [0, max(bar_data) * 1.2]);

    sgtitle('聚焦对比: 不同网格精度的Vectorized求解器', 'FontSize', 16, 'FontWeight', 'bold');
end