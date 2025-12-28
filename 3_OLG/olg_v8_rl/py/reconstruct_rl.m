% MATLAB Script to Reconstruct RL Network and Compare with VFI

clear; clc; close all;

%% 1. 加载RL模型参数
% -------------------------------------------------------------------------
disp('--- 1. 加载RL模型参数 ---');
paramsFile = './py/best_model_sb3/best_model_params.json';
if ~isfile(paramsFile)
    error('RL参数文件未找到: %s', paramsFile);
end
params = jsondecode(fileread(paramsFile));
disp('RL参数加载成功!');

%% 2. 加载OLG模型基础参数
% -------------------------------------------------------------------------
disp('--- 2. 加载OLG模型基础参数 ---');
% 确保MATLAB路径包含 main_olg_v8_utils
if ~(exist('main_olg_v8_utils', 'class') == 8)
    error('未找到 main_olg_v8_utils 类。请确保已将其路径添加到MATLAB中。');
end
cS = main_olg_v8_utils.ParameterValues_HuggettStyle();
disp('OLG基础参数加载成功!');

%% 3. 重建RL网络
% -------------------------------------------------------------------------
disp('--- 3. 重建RL网络 ---');
% 从ONNX文件加载网络
onnxFile = './py/best_model_sb3/best_model.onnx';
if ~isfile(onnxFile)
    error('MATLAB: ONNX模型文件未找到: %s', onnxFile);
end
try
    disp(['MATLAB: 正在从 ' onnxFile ' 导入ONNX网络...']);
    % 使用推荐的 importNetworkFromONNX 函数, 并移除无效的 'TargetNetwork' 参数
    rl_net = importNetworkFromONNX(onnxFile, "InputDataFormats", "BC", "OutputDataFormats", "BC");
    disp('MATLAB: ONNX网络导入成功!');
catch ME
    disp('MATLAB: ONNX网络导入失败。');
    disp(ME.message);
    rethrow(ME);
end

%% 4. 运行VFI vs RL 对比分析
% -------------------------------------------------------------------------
disp('--- 4. 运行VFI vs RL 对比分析 ---');
n_sim = 500; % 模拟个体数量
compare_rl_vs_vfi(rl_net, params, cS, n_sim);

disp('====================================================');
disp('脚本执行完成!');
disp('====================================================');

%% COMPARISON FUNCTION
% =========================================================================
function results = compare_rl_vs_vfi(rl_net, rl_params, cS, n_sim)
    % 主函数，用于比较RL和VFI的生命周期表现
    
    % --- Part 1: 运行VFI求解器 ---
    disp(' ');
    disp('--- 步骤 1: 运行MATLAB VFI求解器 ---');
    vfi_results = run_vfi_solver(cS);
    if ~vfi_results.success
        error('VFI求解失败，无法继续比较。');
    end

    % --- Part 2: 模拟生命周期 ---
    disp(' ');
    disp('--- 步骤 2: 模拟生命周期轨迹 ---');
    sim_results = simulate_lifecycles(rl_net, rl_params, vfi_results, cS, n_sim);

    % --- Part 3: 计算统计数据 ---
    disp(' ');
    disp('--- 步骤 3: 计算和显示统计结果 ---');
    results = analyze_simulation_results(sim_results);
    
    % --- Part 4: 绘制图表 ---
    disp(' ');
    disp('--- 步骤 4: 绘制比较图表 ---');
    plot_comparison(results);
end


%% HELPER FUNCTIONS
% =========================================================================

function vfi_results = run_vfi_solver(cS)
    % 运行VFI求解器，获取策略函数
    
    M_test = struct(...
        'R_k_net_factor', 1.03, ...
        'w_gross', 2.0, ...
        'TR_total', 0.1, ...
        'b_payg_avg_retiree', 0.4, ...
        'tau_l', 0.15, ...
        'theta_payg_actual', 0.12 ...
    );
    
    [leLogGridV, leTrProbM, leProb1V] = main_olg_v8_utils.EarningProcess_olgm(cS);
    paramS_vfi = struct(...
        'leLogGridV', leLogGridV, ...
        'leTrProbM', leTrProbM, ...
        'leProb1V', leProb1V, ...
        'leGridV', exp(leLogGridV), ...
        'ageEffV_new', cS.ageEffV_new, ...
        'tau_l', M_test.tau_l, ...
        'theta_payg_actual_for_hh', M_test.theta_payg_actual, ...
        'pps_tax_deferral_active', cS.pps_active ...
    );

    bV_payg_vfi = zeros(1, cS.aD_new);
    if cS.aR_new < cS.aD_new
        bV_payg_vfi(cS.aR_new:end) = M_test.b_payg_avg_retiree;
    end

    tic;
    [cPolM, kPolM, cPpsPolM, VPolM] = main_olg_v8_utils.HHSolution_VFI_Huggett(...
        M_test.R_k_net_factor, M_test.w_gross, M_test.TR_total, bV_payg_vfi, paramS_vfi, cS);
    vfi_time = toc;
    
    fprintf('VFI求解完成，耗时: %.2f 秒\n', vfi_time);

    vfi_results = struct(...
        'success', true, ...
        'cPolM', cPolM, 'kPolM', kPolM, 'cPpsPolM', cPpsPolM, 'VPolM', VPolM, ...
        'paramS_vfi', paramS_vfi, 'vfi_time', vfi_time, 'M_test', M_test...
    );
end

function sim_results = simulate_lifecycles(rl_net, rl_params, vfi_results, cS, n_sim)
    % 模拟RL和VFI两种策略下的生命周期
    
    % 设置随机数种子以保证结果可复现
    rng(42); 
    
    % 提取所需参数
    aD_new = cS.aD_new;
    beta = cS.beta;
    kMin = cS.kMin;
    kppsMin = cS.kppsMin;
    leProb1V = vfi_results.paramS_vfi.leProb1V;
    leTrProbM = vfi_results.paramS_vfi.leTrProbM;
    kGridV = cS.kGridV(:);
    kppsGridV = cS.kppsGridV(:);
    nw = cS.nw;

    % 获取与Python脚本对齐的PPS回报率
    pps_return_premium = 0;
    if isfield(cS, 'pps_return_rate_premium')
        pps_return_premium = cS.pps_return_rate_premium;
    end
    R_pps = 1 + (vfi_results.M_test.R_k_net_factor - 1) + pps_return_premium;
    fprintf('信息: 使用与Python对齐的PPS回报率 (R_pps = %.4f)\n', R_pps);

    % 提取动作空间边界
    action_low = rl_params.action_space_low;
    action_high = rl_params.action_space_high;
    
    % 初始化结果矩阵
    lifetime_utility_rl = zeros(n_sim, 1);
    k_path_rl = zeros(n_sim, aD_new);
    c_path_rl = zeros(n_sim, aD_new);
    cpps_path_rl = zeros(n_sim, aD_new);

    lifetime_utility_vfi = zeros(n_sim, 1);
    k_path_vfi = zeros(n_sim, aD_new);
    c_path_vfi = zeros(n_sim, aD_new);
    cpps_path_vfi = zeros(n_sim, aD_new);
    
    fprintf('开始模拟 %d 个个体的生命周期...\n', n_sim);
    
    for i_sim = 1:n_sim
        k_rl = kMin; kpps_rl = kppsMin;
        k_vfi = kMin; kpps_vfi = kppsMin;
        
        % 初始效率冲击
        eps_idx = find(rand <= cumsum(leProb1V), 1, 'first');
        
        utility_sum_rl = 0;
        utility_sum_vfi = 0;
        
        for age_idx = 1:aD_new
            % --- RL 模拟 ---
            obs_rl = [k_rl; kpps_rl; age_idx; eps_idx; ...
                      vfi_results.M_test.R_k_net_factor-1; ...
                      vfi_results.M_test.w_gross; ...
                      vfi_results.M_test.TR_total; ...
                      vfi_results.M_test.b_payg_avg_retiree; ...
                      cS.pps_tax_rate_withdrawal; ...
                      age_idx >= cS.aR_new];
            
            action_rl = predict_for_olg(rl_net, obs_rl, action_low, action_high);
            % 根据错误信息，网络只输出2个值 (c, k_next)。
            c_rl = action_rl(1); 
            k_next_rl = action_rl(2); 
            % 假设 cpps_rl不由RL网络直接决定。在没有明确规则的情况下，暂时设为0。
            cpps_rl = 0;
            
            k_path_rl(i_sim, age_idx) = k_rl;
            c_path_rl(i_sim, age_idx) = c_rl;
            cpps_path_rl(i_sim, age_idx) = cpps_rl;
            [~, u_rl] = main_olg_v8_utils.CES_utility(c_rl, cS.sigma, cS);
            utility_sum_rl = utility_sum_rl + (beta^(age_idx-1)) * u_rl;
            
            % --- VFI 模拟 ---
            k_idx = find(abs(kGridV - k_vfi) == min(abs(kGridV - k_vfi)), 1);
            kpps_idx = find(abs(kppsGridV - kpps_vfi) == min(abs(kppsGridV - kpps_vfi)), 1);
            
            c_vfi = vfi_results.cPolM(k_idx, kpps_idx, eps_idx, age_idx);
            k_next_vfi = vfi_results.kPolM(k_idx, kpps_idx, eps_idx, age_idx);
            cpps_vfi = vfi_results.cPpsPolM(k_idx, kpps_idx, eps_idx, age_idx);

            k_path_vfi(i_sim, age_idx) = k_vfi;
            c_path_vfi(i_sim, age_idx) = c_vfi;
            cpps_path_vfi(i_sim, age_idx) = cpps_vfi;
            [~, u_vfi] = main_olg_v8_utils.CES_utility(c_vfi, cS.sigma, cS);
            utility_sum_vfi = utility_sum_vfi + (beta^(age_idx-1)) * u_vfi;
            
            % --- 状态转移 (为下一期做准备) ---
            % 更新下一期的资产状态 (RL)
            k_rl = k_next_rl;
            % 修正 kpps 的累积方式，使其与 Python 版本对齐 (先加后乘，并限制范围)
            kpps_rl = (kpps_rl + cpps_rl) * R_pps;
            kpps_rl = max(cS.kppsMin, min(cS.kppsMax, kpps_rl));

            % 更新下一期的资产状态 (VFI)
            k_vfi = k_next_vfi;
            % 修正 kpps 的累积方式，使其与 Python 版本对齐 (先加后乘，并限制范围)
            kpps_vfi = (kpps_vfi + cpps_vfi) * R_pps;
            kpps_vfi = max(cS.kppsMin, min(cS.kppsMax, kpps_vfi));
            
            % 更新下一期的效率冲击
            if age_idx < aD_new
                eps_idx = find(rand <= cumsum(leTrProbM(eps_idx,:)), 1, 'first');
            end
        end
        lifetime_utility_rl(i_sim) = utility_sum_rl;
        lifetime_utility_vfi(i_sim) = utility_sum_vfi;
    end
    
    sim_results = pack_results(lifetime_utility_rl, lifetime_utility_vfi, ...
        k_path_rl, k_path_vfi, c_path_rl, c_path_vfi, cpps_path_rl, cpps_path_vfi);
end

function results_pack = pack_results(varargin)
    % 将多个结果打包成一个结构体
    names = {'lifetime_utility_rl', 'lifetime_utility_vfi', ...
        'k_path_rl', 'k_path_vfi', 'c_path_rl', 'c_path_vfi', 'cpps_path_rl', 'cpps_path_vfi'};
    results_pack = struct();
    for i = 1:length(varargin)
        results_pack.(names{i}) = varargin{i};
    end
end

function stats = analyze_simulation_results(sim_results)
    % 计算统计数据
    stats = sim_results;
    stats.mean_utility_rl = mean(sim_results.lifetime_utility_rl);
    stats.std_utility_rl = std(sim_results.lifetime_utility_rl);
    stats.mean_utility_vfi = mean(sim_results.lifetime_utility_vfi);
    stats.std_utility_vfi = std(sim_results.lifetime_utility_vfi);
    stats.utility_diff = stats.mean_utility_rl - stats.mean_utility_vfi;
    stats.utility_improvement_pct = (stats.utility_diff / abs(stats.mean_utility_vfi)) * 100;
    [~, p_value] = ttest(sim_results.lifetime_utility_rl, sim_results.lifetime_utility_vfi);
    stats.p_value = p_value;

    fprintf('\n--- 结果摘要 ---\n');
    fprintf('VFI  平均效用: %.4f ± %.4f\n', stats.mean_utility_vfi, stats.std_utility_vfi);
    fprintf('RL   平均效用: %.4f ± %.4f\n', stats.mean_utility_rl, stats.std_utility_rl);
    fprintf('差异 (RL-VFI): %.4f (%.2f%%)\n', stats.utility_diff, stats.utility_improvement_pct);
    % fprintf('p-value: %.4f (%s)\n', p_value, if(p_value < 0.05, '显著', '不显著'));
end

function plot_comparison(results)
    % 绘制比较图表
    fig = figure('Name', 'VFI vs RL Comparison', 'Position', [100, 100, 1200, 800]);
    sgtitle(sprintf('VFI vs RL 生命周期比较 (n=%d)', size(results.k_path_rl,1)), 'FontSize', 16);

    % 子图1: 效用分布
    subplot(2,3,1);
    hold on;
    histogram(results.lifetime_utility_vfi, 'FaceColor', 'r', 'FaceAlpha', 0.7);
    histogram(results.lifetime_utility_rl, 'FaceColor', 'b', 'FaceAlpha', 0.7);
    xline(results.mean_utility_vfi, 'r--', 'LineWidth', 2);
    xline(results.mean_utility_rl, 'b--', 'LineWidth', 2);
    title('效用分布比较'); xlabel('生涯总效用'); ylabel('频数');
    legend('VFI', 'RL', 'VFI Mean', 'RL Mean');
    grid on; hold off;

    % 子图2: 个体效用对比
    subplot(2,3,2);
    scatter(results.lifetime_utility_vfi, results.lifetime_utility_rl, 15, 'filled', 'MarkerFaceAlpha', 0.6);
    hold on;
    ref_min = min(min(results.lifetime_utility_vfi), min(results.lifetime_utility_rl));
    ref_max = max(max(results.lifetime_utility_vfi), max(results.lifetime_utility_rl));
    plot([ref_min, ref_max], [ref_min, ref_max], 'k--', 'LineWidth', 1.5);
    title('个体效用对比'); xlabel('VFI 效用'); ylabel('RL 效用');
    grid on; hold off;
    
    % 子图3: 效用差异分布
    subplot(2,3,3);
    utility_diffs = results.lifetime_utility_rl - results.lifetime_utility_vfi;
    histogram(utility_diffs, 'FaceColor', 'g', 'FaceAlpha', 0.7);
    hold on;
    xline(mean(utility_diffs), 'g-', 'LineWidth', 2, 'Label', sprintf('Mean Diff: %.4f', mean(utility_diffs)));
    xline(0, 'k--', 'LineWidth', 1);
    title('效用差异 (RL - VFI)'); xlabel('效用差异'); ylabel('频数');
    grid on; hold off;

    % 子图4: 平均资产路径
    subplot(2,3,4);
    plot(mean(results.k_path_vfi, 1), 'r-', 'LineWidth', 2);
    hold on;
    plot(mean(results.k_path_rl, 1), 'b--', 'LineWidth', 2);
    title('平均资产路径'); xlabel('年龄组'); ylabel('资产');
    legend('VFI', 'RL');
    grid on; hold off;

    % 子图5: 平均消费路径
    subplot(2,3,5);
    plot(mean(results.c_path_vfi, 1), 'r-', 'LineWidth', 2);
    hold on;
    plot(mean(results.c_path_rl, 1), 'b--', 'LineWidth', 2);
    title('平均消费路径'); xlabel('年龄组'); ylabel('消费');
    legend('VFI', 'RL');
    grid on; hold off;

    % 子图6: 平均PPS缴费路径
    subplot(2,3,6);
    plot(mean(results.cpps_path_vfi, 1), 'r-', 'LineWidth', 2);
    hold on;
    plot(mean(results.cpps_path_rl, 1), 'b--', 'LineWidth', 2);
    title('平均PPS缴费路径'); xlabel('年龄组'); ylabel('PPS缴费');
    legend('VFI', 'RL');
    grid on; hold off;
end

function action = predict_for_olg(net, observation, low, high)
    % 封装的预测函数，用于OLG模拟 (与predict_from_reconstructed_net一致)
    
    % 为ONNX导入的网络，输入需要是 (Batch, Features) 格式
    observation_row = single(observation(:)'); % 确保是行向量 (1, num_features)
    dlX = dlarray(observation_row, 'BC');      % 'BC' for (Batch, Channel/Features)
    
    % 进行预测
    dlY = predict(net, dlX);
    raw_action = extractdata(dlY);
    
    % 将动作从 [-1, 1] 缩放到 [low, high]
    % 注意：要确保 low 和 high 是列向量以便正确广播
    action = low(:) + (high(:) - low(:)) .* (raw_action(:) + 1) / 2;
end