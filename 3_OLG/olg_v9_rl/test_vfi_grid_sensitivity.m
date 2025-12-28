% --- test_vfi_grid_sensitivity.m ---
% VFI 网格密度敏感性测试脚本 (纯 MATLAB 版本)
%
% 目标：
% 1. 在纯MATLAB环境中，系统性地研究 nk 和 nkpps 的密度对VFI策略性能的影响。
% 2. 对比 'linear' 和 'spline' 两种插值方法下的表现。
% 3. 可视化结果，以诊断性能反常的根本原因。

clear;
close all;
clc;

% --- 1. 设置测试参数 ---
fprintf('--- 1. 设置测试参数 ---\n');

% 定义要测试的网格配置 {nk, nkpps}
grid_configurations = [
    5, 5;
];

% 定义要测试的插值方法
interpolation_methods = { 'linear'};

% 模拟参数
n_sim = 500;
random_seed = 42;
rng(random_seed); % 设置随机种子以保证效率冲击路径一致

% --- 2. 准备固定的宏观环境和效率冲击路径 ---
fprintf('--- 2. 准备宏观环境和效率路径 ---\n');

% 创建一个临时的 cS 来生成效率冲击路径
cS_temp = main_olg_v8_utils.ParameterValues_HuggettStyle();
paramS_temp = struct();
[paramS_temp.leLogGridV, paramS_temp.leTrProbM, paramS_temp.leProb1V] = main_olg_v8_utils.EarningProcess_olgm(cS_temp);
paramS_temp.leGridV = exp(paramS_temp.leLogGridV);

% 生成一个统一的效率冲击路径，用于所有测试
eIdxM_group_global = main_olg_v8_utils.MarkovChainSimulation_AgeGroup(n_sim, cS_temp, paramS_temp.leProb1V, paramS_temp.leTrProbM);
fprintf('✅ 全局效率路径已生成 (size: %d x %d)。\n', size(eIdxM_group_global, 1), size(eIdxM_group_global, 2));

% 定义固定的宏观环境
M_FIXED = struct(...
    'R_k_net_factor', 1.03, ...
    'w_gross', 2.0, ...
    'TR_total', 0.1, ...
    'b_payg_avg_retiree', 0.4, ...
    'tau_l', 0.15, ...
    'theta_payg_actual', 0.12 ...
);

% --- 3. 运行敏感性分析循环 ---
all_results = {};
result_idx = 1;

for i_grid = 1:size(grid_configurations, 1)
    nk = grid_configurations(i_grid, 1);
    nkpps = grid_configurations(i_grid, 2);

    for i_interp = 1:length(interpolation_methods)
        interp_method = interpolation_methods{i_interp};
        
        fprintf('\n============================================================\n');
        fprintf('运行配置: nk=%d, nkpps=%d, interp=''%s''\n', nk, nkpps, interp_method);
        fprintf('============================================================\n');
        
        % a. 为当前配置生成参数
        cS = main_olg_v8_utils.ParameterValues_HuggettStyle();
        cS.nk = nk;
        cS.nkpps = nkpps;
        % [核心] 设定插值方法
        cS.interpolation_method = interp_method; 
        
        % 重新生成依赖于网格的参数
        cS = main_olg_v8_utils.generateGrids(cS);
        
        paramS = struct();
        [paramS.leLogGridV, paramS.leTrProbM, paramS.leProb1V] = main_olg_v8_utils.EarningProcess_olgm(cS);
        paramS.leGridV = exp(paramS.leLogGridV);
        paramS.tau_l = M_FIXED.tau_l;
        paramS.theta_payg_actual_for_hh = M_FIXED.theta_payg_actual;
        paramS.pps_tax_deferral_active = cS.pps_active;
        
        bV_payg_vfi = zeros(cS.aD_new, 1);
        if cS.aR_new < cS.aD_new
            bV_payg_vfi(cS.aR_new:end) = M_FIXED.b_payg_avg_retiree;
        end

        % b. 求解VFI
        tic;
        [cPolM_q, kPolM, cPpsPolM_choice, ~] = main_olg_v8_utils.HHSolution_VFI_Huggett(...
            M_FIXED.R_k_net_factor, M_FIXED.w_gross, M_FIXED.TR_total, ...
            bV_payg_vfi, paramS, cS,'grid');
        fprintf('VFI求解耗时: %.2f 秒。\n', toc);

        % c. 模拟生命周期路径
        tic;
        [k_path, kpps_path, c_path, cpps_path] = main_olg_v8_utils.HHSimulation_olgm(...
            kPolM, cPpsPolM_choice, cPolM_q, eIdxM_group_global, ...
            M_FIXED.R_k_net_factor, M_FIXED.w_gross, M_FIXED.TR_total, ...
            bV_payg_vfi, paramS, cS);
        fprintf('模拟耗时: %.2f 秒。\n', toc);

        % d. 计算终身效用
        utility_vfi = calculate_lifetime_utility(c_path, cS, true);
        
        % e. 存储结果
        all_results{result_idx} = struct(...
            'nk', nk, ...
            'nkpps', nkpps, ...
            'interp_method', interp_method, ...
            'mean_utility', mean(utility_vfi), ...
            'std_utility', std(utility_vfi), ...
            'c_path', c_path, ...
            'k_path', k_path ...
        );
        result_idx = result_idx + 1;
        
        fprintf('📈 结果: 平均效用 = %.4f\n', mean(utility_vfi));
    end
end

% --- 4. 可视化结果 ---
fprintf('\n--- 4. 可视化结果 ---\n');
plot_sensitivity_results(all_results);


% --- 辅助函数 ---
function lifetime_utility = calculate_lifetime_utility(c_path, cS, use_survival_prob)
    [n_sim, aD] = size(c_path);
    lifetime_utility = zeros(n_sim, 1);
    beta = cS.beta;
    s_transitionV = cS.s_1yr_transitionV;

    for i_sim = 1:n_sim
        utility_sum = 0.0;
        cumulative_discount = 1.0;
        for a_group = 1:aD
            [~, u] = main_olg_v8_utils.CES_utility(c_path(i_sim, a_group), cS.sigma, cS);
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

function plot_sensitivity_results(results)
    figure('Name', 'VFI Grid Sensitivity Analysis', 'Position', [100, 100, 1400, 900]);
    
    % [核心修正] 使用 cellfun 来创建逻辑索引，以安全地筛选元胞数组
    
    % 1. 找到所有 interp_method 为 'linear' 的元素的索引
    is_linear_idx = cellfun(@(x) strcmp(x.interp_method, 'linear'), results);
    
    % 2. 找到所有 interp_method 为 'spline' 的元素的索引
    is_spline_idx = cellfun(@(x) strcmp(x.interp_method, 'spline'), results);
    
    % 3. 使用逻辑索引来筛选出对应的元胞
    linear_results_cell = results(is_linear_idx);
    spline_results_cell = results(is_spline_idx);
    
    % [可选但推荐] 将筛选后的元胞数组转换为结构体数组，便于后续访问
    if ~isempty(linear_results_cell)
        linear_results = [linear_results_cell{:}];
    else
        linear_results = [];
    end
    if ~isempty(spline_results_cell)
        spline_results = [spline_results_cell{:}];
    else
        spline_results = [];
    end

    % --- 后续的绘图代码现在可以正常工作了 ---
    
    % 1. 效用 vs 网格点总数
    subplot(2, 2, 1);
    hold on;
    if ~isempty(linear_results)
        % 现在 linear_results 是一个结构体数组，可以直接用点索引
        total_points_lin = [linear_results.nk] .* [linear_results.nkpps];
        mean_utils_lin = [linear_results.mean_utility];
        plot(total_points_lin, mean_utils_lin, 'o-', 'LineWidth', 1.5, 'DisplayName', 'linear');
    end
    if ~isempty(spline_results)
        total_points_spl = [spline_results.nk] .* [spline_results.nkpps];
        mean_utils_spl = [spline_results.mean_utility];
        plot(total_points_spl, mean_utils_spl, 's--', 'LineWidth', 1.5, 'DisplayName', 'spline');
    end
    hold off;
    xlabel('总状态点数 (nk * nkpps)');
    ylabel('平均终身效用');
    title('VFI性能 vs. 网格总点数');
    legend('show', 'Location', 'best');
    grid on;

    % ... 其他子图的逻辑也需要做类似调整 ...
    % 例如，在绘制消费路径时：
    
    % 3. 平均消费路径对比 (linear)
    subplot(2, 2, 3);
    hold on;
    if ~isempty(linear_results)
        colors = parula(length(linear_results));
        for i = 1:length(linear_results)
            r = linear_results(i); % 直接索引结构体数组
            mean_c_path = mean(r.c_path, 1);
            plot(mean_c_path, 'Color', colors(i,:), 'LineWidth', 1.5, 'DisplayName', sprintf('(%d,%d)', r.nk, r.nkpps));
        end
    end
    hold off;
    title("平均消费路径 (interp='linear')");
    % ...
    
    % ... (对所有使用 results 的地方进行类似修改)
end