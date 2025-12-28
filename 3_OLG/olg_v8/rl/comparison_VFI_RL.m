% --- comparison_VFI_RL.m ---
% VFI与RL结果比较脚本
% 读取RL训练结果并与VFI方法进行比较

clc;
clear;
close all;

fprintf('=== VFI vs RL Comparison Analysis ===\n');

%% 1. 读取RL结果
fprintf('\n--- 1. Loading RL Results ---\n');

result_folder = 'rl_result';
if ~exist(result_folder, 'dir')
    error('RL result folder "%s" not found. Please run main_olg_v8_rl_demo.m first.', result_folder);
end

% 读取训练好的智能体
agent_file = fullfile(result_folder, 'trained_agent.mat');
if exist(agent_file, 'file')
    load(agent_file, 'agent');
    fprintf('  ✓ Loaded trained agent from: %s\n', agent_file);
else
    error('Trained agent file not found: %s', agent_file);
end

% 读取环境参数
env_params_file = fullfile(result_folder, 'environment_params.mat');
if exist(env_params_file, 'file')
    load(env_params_file, 'cS_v8', 'fixedMacro', 'envConstantParams');
    fprintf('  ✓ Loaded environment parameters from: %s\n', env_params_file);
else
    error('Environment parameters file not found: %s', env_params_file);
end

% 读取策略函数数据
policy_file = fullfile(result_folder, 'policy_data.mat');
if exist(policy_file, 'file')
    load(policy_file, 'policy_data');
    fprintf('  ✓ Loaded policy data from: %s\n', policy_file);
else
    fprintf('  Warning: Policy data file not found: %s\n', policy_file);
    policy_data = struct();
end

% 读取训练统计（可选）
training_file = fullfile(result_folder, 'training_stats.mat');
if exist(training_file, 'file')
    load(training_file, 'trainingStats');
    fprintf('  ✓ Loaded training statistics from: %s\n', training_file);
else
    fprintf('  Note: Training statistics file not found: %s\n', training_file);
end

% 读取时间戳
timestamp_file = fullfile(result_folder, 'timestamp.mat');
if exist(timestamp_file, 'file')
    load(timestamp_file, 'timestamp');
    fprintf('  ✓ RL results timestamp: %s\n', timestamp.datestr);
else
    fprintf('  Note: Timestamp file not found\n');
end

%% 2. 重新创建环境（用于RL策略评估）
fprintf('\n--- 2. Recreating Environment ---\n');
env = OLGV8EnvDemo(envConstantParams);
fprintf('  ✓ Environment recreated successfully\n');

%% 3. 运行VFI求解
fprintf('\n--- 3. Running VFI Solution ---\n');

% 构建VFI所需的参数结构
paramS_vfi = struct();
paramS_vfi.tau_l = fixedMacro.tau_l;
paramS_vfi.theta_payg_actual_for_hh = fixedMacro.theta_payg_actual_for_hh;
paramS_vfi.pps_tax_deferral_active = cS_v8.pps_active;

% 生成劳动禀赋过程参数
try
    [leLogGridV_temp, leTrProbM_temp, leProb1V_temp] = main_olg_v8_utils.EarningProcess_olgm(cS_v8);
    paramS_vfi.leGridV = exp(leLogGridV_temp(:));
    paramS_vfi.leTrProbM = leTrProbM_temp;
    paramS_vfi.leProb1V = leProb1V_temp;
    fprintf('  ✓ Generated earning process parameters\n');
    fprintf('    leGridV length: %d, cS.nw: %d\n', length(paramS_vfi.leGridV), cS_v8.nw);
catch ME
    fprintf('  Warning: Failed to generate earning process: %s\n', ME.message);
    fprintf('  Using simplified parameters\n');
    % 使用简化的参数，确保长度与cS_v8.nw匹配
    if cS_v8.nw == 3
        paramS_vfi.leGridV = [0.8, 1.0, 1.2];
        paramS_vfi.leTrProbM = [0.7 0.2 0.1; 0.2 0.6 0.2; 0.1 0.2 0.7];
        paramS_vfi.leProb1V = [0.33, 0.34, 0.33];
    else
        % 为其他nw值生成参数
        paramS_vfi.leGridV = linspace(0.8, 1.2, cS_v8.nw);
        paramS_vfi.leTrProbM = eye(cS_v8.nw) * 0.8 + ones(cS_v8.nw) * 0.2 / cS_v8.nw;
        paramS_vfi.leProb1V = ones(1, cS_v8.nw) / cS_v8.nw;
    end
    fprintf('    Generated simplified parameters: leGridV length: %d, cS.nw: %d\n', length(paramS_vfi.leGridV), cS_v8.nw);
end

% 确保leGridV长度与cS.nw匹配
if length(paramS_vfi.leGridV) ~= cS_v8.nw
    fprintf('  Warning: leGridV length (%d) != cS.nw (%d), adjusting...\n', length(paramS_vfi.leGridV), cS_v8.nw);
    if length(paramS_vfi.leGridV) > cS_v8.nw
        paramS_vfi.leGridV = paramS_vfi.leGridV(1:cS_v8.nw);
    else
        % 扩展到所需长度
        paramS_vfi.leGridV = [paramS_vfi.leGridV(:)', ones(1, cS_v8.nw - length(paramS_vfi.leGridV))];
    end
    fprintf('    Adjusted leGridV length: %d\n', length(paramS_vfi.leGridV));
end

% 确保bV_payg向量长度正确
if length(fixedMacro.bV_payg) < cS_v8.aD_new
    fprintf('  Extending bV_payg vector from length %d to %d\n', length(fixedMacro.bV_payg), cS_v8.aD_new);
    bV_payg_extended = zeros(1, cS_v8.aD_new);
    bV_payg_extended(1:length(fixedMacro.bV_payg)) = fixedMacro.bV_payg;
    % 为剩余的年龄组设置默认值（退休年龄组获得福利）
    retirement_start_group = ceil((cS_v8.ageRetire_orig - cS_v8.age1_orig + 1) / cS_v8.yearStep) + 1;
    if retirement_start_group <= cS_v8.aD_new
        bV_payg_extended(retirement_start_group:cS_v8.aD_new) = 0.2;
    end
else
    bV_payg_extended = fixedMacro.bV_payg(1:cS_v8.aD_new);
end

% 验证网格参数一致性
fprintf('  Grid verification before VFI: nk=%d, length(kGridV)=%d, nkpps=%d, length(kppsGridV)=%d, nw=%d\n', ...
    cS_v8.nk, length(cS_v8.kGridV), cS_v8.nkpps, length(cS_v8.kppsGridV), cS_v8.nw);

% 强制确保网格长度匹配
if length(cS_v8.kGridV) ~= cS_v8.nk
    fprintf('    Warning: Regenerating kGridV to match nk\n');
    cS_v8.kGridV = linspace(cS_v8.kMin, cS_v8.kMax, cS_v8.nk)';
end
if length(cS_v8.kppsGridV) ~= cS_v8.nkpps
    fprintf('    Warning: Regenerating kppsGridV to match nkpps\n');
    cS_v8.kppsGridV = linspace(cS_v8.kppsMin, cS_v8.kppsMax, cS_v8.nkpps)';
end

% 添加VFI所需的其他参数
cS_v8.fminbnd_TolX = 1e-6;
cS_v8.fminbnd_Display = 'off';
cS_v8.n_pps_choice_grid_points = 5;
cS_v8.pps_tax_rate_withdrawal = 0.1;
cS_v8.aR_new = ceil(2/3 * cS_v8.aD_new); % 工作年龄组数量

% 调用VFI求解器
try
    fprintf('  Running VFI solver...\n');
    tic;
    [cPolM_vfi, kPolM_vfi, cPpsPolM_vfi, valM_vfi] = main_olg_v8_utils.HHSolution_VFI_Huggett(...
        fixedMacro.R_k_net_factor_hh, fixedMacro.MPL_gross, fixedMacro.TR_total, ...
        bV_payg_extended, paramS_vfi, cS_v8);
    vfi_time = toc;
    fprintf('  ✓ VFI solution completed successfully in %.2f seconds\n', vfi_time);
    vfi_available = true;
catch ME
    fprintf('  ✗ VFI solution failed: %s\n', ME.message);
    fprintf('    Detailed error info:\n');
    fprintf('      Message: %s\n', ME.message);
    if ~isempty(ME.stack)
        fprintf('      Location: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
    end
    fprintf('    Will skip VFI comparison\n');
    vfi_available = false;
end

%% 4. 策略函数比较
if vfi_available
    fprintf('\n--- 4. Policy Function Comparison ---\n');
    
    % 选择比较的状态点
    compare_a_idx = min(round(cS_v8.aD_new / 2), cS_v8.aD_new);
    if compare_a_idx == 0, compare_a_idx = 1; end
    compare_ie_idx = round(cS_v8.nw / 2);
    compare_ikpps_idx = round(cS_v8.nkpps / 2);
    if compare_ikpps_idx == 0, compare_ikpps_idx = 1; end
    
    fprintf('  Comparing at age_idx=%d, eps_idx=%d, k_pps_idx=%d\n', ...
        compare_a_idx, compare_ie_idx, compare_ikpps_idx);
    
    % 创建资产网格用于比较
    k_compare_grid = linspace(cS_v8.kMin, cS_v8.kMax, 15);
    k_pps_compare_val = cS_v8.kppsGridV(compare_ikpps_idx);
    
    % 提取VFI策略
    k_prime_vfi = zeros(size(k_compare_grid));
    c_pps_vfi = zeros(size(k_compare_grid));
    
    for i = 1:length(k_compare_grid)
        % 找到最近的网格点
        [~, k_idx] = min(abs(cS_v8.kGridV - k_compare_grid(i)));
        k_prime_vfi(i) = kPolM_vfi(k_idx, compare_ikpps_idx, compare_ie_idx, compare_a_idx);
        c_pps_vfi(i) = cPpsPolM_vfi(k_idx, compare_ikpps_idx, compare_ie_idx, compare_a_idx);
    end
    
    % 获取RL策略
    k_prime_rl = zeros(size(k_compare_grid));
    c_pps_rl = zeros(size(k_compare_grid));
    
    for i = 1:length(k_compare_grid)
        % 构建状态
        state_unnorm = [k_compare_grid(i); k_pps_compare_val; compare_ie_idx; compare_a_idx];
        
        % 归一化状态
        try
            obs_norm = env.normalizeObservation(state_unnorm);
        catch
            obs_norm = normalizeState(state_unnorm, cS_v8);
        end
        
        % 获取RL动作
        action_norm = getAction(agent, obs_norm);
        if iscell(action_norm)
            action_vec = action_norm{1};
        else
            action_vec = action_norm;
        end
        
        % 反归一化动作得到实际策略值
        % 对于k_prime: 从[-1,1]映射到[kMin, kMax]
        k_prime_rl(i) = cS_v8.kMin + (action_vec(1) + 1)/2 * (cS_v8.kMax - cS_v8.kMin);
        
        % 对于c_pps: 需要考虑约束
        % 简化处理：假设在允许范围内线性映射
        epsilon_val_for_cpps = 1.0; % 默认值
        if isfield(paramS_vfi, 'leGridV') && length(paramS_vfi.leGridV) >= compare_ie_idx
            epsilon_val_for_cpps = paramS_vfi.leGridV(compare_ie_idx);
        end
        max_cpps_approx = fixedMacro.MPL_gross * cS_v8.ageEffV_new(compare_a_idx) * ...
            epsilon_val_for_cpps * cS_v8.pps_max_contrib_frac;
        max_cpps_approx = min(max_cpps_approx, cS_v8.pps_annual_contrib_limit);
        c_pps_rl(i) = (action_vec(2) + 1)/2 * max_cpps_approx;
    end
    
    % 绘制策略函数比较
    figure('Name', 'VFI vs RL Policy Comparison', 'Position', [100, 100, 1200, 500]);
    
    subplot(1,2,1);
    plot(k_compare_grid, k_prime_vfi, 'b-o', 'LineWidth', 2, 'DisplayName', 'VFI');
    hold on;
    plot(k_compare_grid, k_prime_rl, 'r--s', 'LineWidth', 2, 'DisplayName', 'RL');
    plot(k_compare_grid, k_compare_grid, 'k:', 'LineWidth', 1, 'DisplayName', '45° line');
    hold off;
    xlabel('Current Assets k');
    ylabel('Next Period Assets k''');
    title(sprintf('Asset Policy: k''(k) at age=%d, ε=%d', compare_a_idx, compare_ie_idx));
    legend('Location', 'best');
    grid on;
    
    subplot(1,2,2);
    plot(k_compare_grid, c_pps_vfi, 'b-o', 'LineWidth', 2, 'DisplayName', 'VFI');
    hold on;
    plot(k_compare_grid, c_pps_rl, 'r--s', 'LineWidth', 2, 'DisplayName', 'RL');
    hold off;
    xlabel('Current Assets k');
    ylabel('PPS Contribution c_{pps}');
    title(sprintf('PPS Policy: c_{pps}(k) at age=%d, ε=%d', compare_a_idx, compare_ie_idx));
    legend('Location', 'best');
    grid on;
    
    % 保存策略比较图
    saveas(gcf, fullfile(result_folder, 'policy_comparison.png'));
    fprintf('  ✓ Policy comparison plot saved\n');
    
    %% 5. 计算策略差异统计
    fprintf('\n--- 5. Policy Difference Statistics ---\n');
    
    k_prime_diff = abs(k_prime_vfi - k_prime_rl);
    c_pps_diff = abs(c_pps_vfi - c_pps_rl);
    
    fprintf('  Asset policy k'' differences:\n');
    fprintf('    Mean absolute difference: %.4f\n', mean(k_prime_diff));
    fprintf('    Max absolute difference:  %.4f\n', max(k_prime_diff));
    fprintf('    Relative difference (%%):  %.2f%%\n', 100 * mean(k_prime_diff) / mean(k_prime_vfi));
    
    fprintf('  PPS policy c_pps differences:\n');
    fprintf('    Mean absolute difference: %.4f\n', mean(c_pps_diff));
    fprintf('    Max absolute difference:  %.4f\n', max(c_pps_diff));
    if mean(c_pps_vfi) > 1e-6
        fprintf('    Relative difference (%%):  %.2f%%\n', 100 * mean(c_pps_diff) / mean(c_pps_vfi));
    end
    
    %% 6. 值函数比较
    fprintf('\n--- 6. Value Function Comparison ---\n');
    
    % 计算RL隐含的值函数（通过模拟剩余生命周期）
    rl_values = zeros(size(k_compare_grid));
    for i = 1:length(k_compare_grid)
        % 修正的值函数估计：从当前年龄开始模拟到生命结束
        env_temp = OLGV8EnvDemo(envConstantParams);
        env_temp.State = [k_compare_grid(i); k_pps_compare_val; compare_ie_idx; compare_a_idx];
        
        cumulative_reward = 0;
        discount_factor = min(0.99, cS_v8.beta^cS_v8.yearStep);
        discount = 1;
        
        % 计算剩余生命周期步数
        remaining_steps = cS_v8.aD_new - compare_a_idx + 1;
        
        for step_idx = 1:remaining_steps
            obs = env_temp.normalizeObservation(env_temp.State);
            action = getAction(agent, obs);
            if iscell(action)
                action = action{1};
            end
            
            [next_obs, reward, isDone, ~] = env_temp.step(action);
            cumulative_reward = cumulative_reward + discount * reward;
            discount = discount * discount_factor;
            
            if isDone
                break;
            end
        end
        rl_values(i) = cumulative_reward;
    end
    
    % 提取VFI值函数
    vfi_values = zeros(size(k_compare_grid));
    for i = 1:length(k_compare_grid)
        [~, k_idx] = min(abs(cS_v8.kGridV - k_compare_grid(i)));
        vfi_values(i) = valM_vfi(k_idx, compare_ikpps_idx, compare_ie_idx, compare_a_idx);
    end
    
    % 绘制值函数比较
    figure('Name', 'Value Function Comparison', 'Position', [200, 200, 800, 400]);
    plot(k_compare_grid, vfi_values, 'b-o', 'LineWidth', 2, 'DisplayName', 'VFI Value');
    hold on;
    plot(k_compare_grid, rl_values, 'r--s', 'LineWidth', 2, 'DisplayName', 'RL Value (simulated)');
    hold off;
    xlabel('Current Assets k');
    ylabel('Value Function V(k)');
    title(sprintf('Value Function Comparison at age=%d, ε=%d (remaining steps=%d)', ...
        compare_a_idx, compare_ie_idx, cS_v8.aD_new - compare_a_idx + 1));
    legend('Location', 'best');
    grid on;
    
    % 保存值函数比较图
    saveas(gcf, fullfile(result_folder, 'value_comparison.png'));
    fprintf('  ✓ Value function comparison plot saved\n');
    
    % 值函数差异统计和调试信息
    value_diff = abs(vfi_values - rl_values);
    fprintf('  Value function differences:\n');
    fprintf('    VFI value range: [%.4f, %.4f]\n', min(vfi_values), max(vfi_values));
    fprintf('    RL value range:  [%.4f, %.4f]\n', min(rl_values), max(rl_values));
    fprintf('    Mean absolute difference: %.4f\n', mean(value_diff));
    fprintf('    Max absolute difference:  %.4f\n', max(value_diff));
    if mean(abs(vfi_values)) > 1e-6
        fprintf('    Relative difference (%%):  %.2f%%\n', 100 * mean(value_diff) / mean(abs(vfi_values)));
    end
    fprintf('    Discount factor used: %.4f\n', discount_factor);
    fprintf('    Remaining life steps: %d\n', cS_v8.aD_new - compare_a_idx + 1);
    
    %% 7. 生成详细的比较报告
    fprintf('\n--- 7. Detailed Comparison Report ---\n');
    generateComparisonReport(k_compare_grid, k_prime_vfi, k_prime_rl, c_pps_vfi, c_pps_rl, ...
        vfi_values, rl_values, compare_a_idx, compare_ie_idx);
    
    %% 8. 保存比较结果
    fprintf('\n--- 8. Saving Comparison Results ---\n');
    
    comparison_results = struct();
    comparison_results.k_compare_grid = k_compare_grid;
    comparison_results.k_prime_vfi = k_prime_vfi;
    comparison_results.k_prime_rl = k_prime_rl;
    comparison_results.c_pps_vfi = c_pps_vfi;
    comparison_results.c_pps_rl = c_pps_rl;
    comparison_results.vfi_values = vfi_values;
    comparison_results.rl_values = rl_values;
    comparison_results.compare_a_idx = compare_a_idx;
    comparison_results.compare_ie_idx = compare_ie_idx;
    comparison_results.compare_ikpps_idx = compare_ikpps_idx;
    comparison_results.vfi_time = vfi_time;
    comparison_results.timestamp = datetime('now');
    
    % 计算统计指标
    comparison_results.stats = struct();
    comparison_results.stats.k_prime_mean_abs_diff = mean(k_prime_diff);
    comparison_results.stats.k_prime_max_abs_diff = max(k_prime_diff);
    comparison_results.stats.k_prime_rel_diff_pct = 100 * mean(k_prime_diff) / mean(k_prime_vfi);
    comparison_results.stats.c_pps_mean_abs_diff = mean(c_pps_diff);
    comparison_results.stats.c_pps_max_abs_diff = max(c_pps_diff);
    if mean(c_pps_vfi) > 1e-6
        comparison_results.stats.c_pps_rel_diff_pct = 100 * mean(c_pps_diff) / mean(c_pps_vfi);
    else
        comparison_results.stats.c_pps_rel_diff_pct = NaN;
    end
    comparison_results.stats.value_mean_abs_diff = mean(value_diff);
    comparison_results.stats.value_max_abs_diff = max(value_diff);
    comparison_results.stats.value_rel_diff_pct = 100 * mean(value_diff) / mean(abs(vfi_values));
    
    comparison_file = fullfile(result_folder, 'comparison_results.mat');
    save(comparison_file, 'comparison_results');
    fprintf('  ✓ Comparison results saved to: %s\n', comparison_file);
    
else
    fprintf('\n--- VFI Comparison Skipped ---\n');
    fprintf('  VFI solution was not available\n');
end

%% 9. 计算效率指标
fprintf('\n--- 9. Computational Efficiency Comparison ---\n');
if exist('trainingStats', 'var')
    fprintf('  RL training episodes: %d\n', length(trainingStats.EpisodeReward));
    fprintf('  RL training time: Available in training monitor\n');
end
if vfi_available
    fprintf('  VFI solution time: %.2f seconds\n', vfi_time);
end
fprintf('  VFI advantage: Guaranteed global optimum (given discretization)\n');
fprintf('  RL advantage: Can handle continuous state/action spaces naturally\n');
fprintf('  RL advantage: No need to re-solve when parameters change slightly\n');

fprintf('\n=== VFI vs RL Comparison Completed ===\n');

%% 本地辅助函数
function obs_norm = normalizeState(state_unnorm, cS)
    % 本地归一化函数
    obs_norm = zeros(4,1);
    
    % 安全地获取边界值
    if isfield(cS, 'kMin')
        kMin_val = cS.kMin;
    else
        kMin_val = 0;
    end
    
    if isfield(cS, 'kMax')
        kMax_val = cS.kMax;
    else
        kMax_val = 10;
    end
    
    if isfield(cS, 'kppsMin')
        kppsMin_val = cS.kppsMin;
    else
        kppsMin_val = 0;
    end
    
    if isfield(cS, 'kppsMax')
        kppsMax_val = cS.kppsMax;
    else
        kppsMax_val = 5;
    end
    
    if isfield(cS, 'nw')
        nw_val = cS.nw;
    else
        nw_val = 3;
    end
    
    if isfield(cS, 'aD_new')
        aD_new_val = cS.aD_new;
    else
        aD_new_val = 10;
    end
    
    obs_norm(1) = (state_unnorm(1) - kMin_val) / (kMax_val - kMin_val + 1e-6);
    obs_norm(2) = (state_unnorm(2) - kppsMin_val) / (kppsMax_val - kppsMin_val + 1e-6);
    obs_norm(3) = (state_unnorm(3)-1) / (nw_val -1 + 1e-6);
    if nw_val == 1, obs_norm(3) = 0.5; end
    obs_norm(4) = (state_unnorm(4)-1) / (aD_new_val -1 + 1e-6);
    obs_norm = max(0, min(1, obs_norm));
end

function generateComparisonReport(k_grid, k_prime_vfi, k_prime_rl, c_pps_vfi, c_pps_rl, ...
    vfi_values, rl_values, age_idx, eps_idx)
    % 生成详细的VFI与RL比较报告
    
    fprintf('  === Detailed Comparison Report ===\n');
    fprintf('    Analysis at age_idx=%d, eps_idx=%d\n', age_idx, eps_idx);
    fprintf('    Asset range: [%.2f, %.2f]\n', min(k_grid), max(k_grid));
    
    % 策略函数相关性分析
    if length(k_prime_vfi) > 2 && length(k_prime_rl) > 2
        corr_k_prime = corrcoef(k_prime_vfi, k_prime_rl);
        fprintf('    Asset policy correlation: %.4f\n', corr_k_prime(1,2));
        
        if any(c_pps_vfi > 1e-6) && any(c_pps_rl > 1e-6)
            corr_c_pps = corrcoef(c_pps_vfi, c_pps_rl);
            fprintf('    PPS policy correlation: %.4f\n', corr_c_pps(1,2));
        else
            fprintf('    PPS policy correlation: N/A (zero contributions)\n');
        end
        
        if length(vfi_values) > 2 && length(rl_values) > 2
            corr_values = corrcoef(vfi_values, rl_values);
            fprintf('    Value function correlation: %.4f\n', corr_values(1,2));
        end
    end
    
    % 策略函数斜率比较
    if length(k_grid) > 1
        dk = k_grid(2) - k_grid(1);
        
        % 资产策略斜率
        slope_vfi_k = diff(k_prime_vfi) / dk;
        slope_rl_k = diff(k_prime_rl) / dk;
        
        fprintf('    Asset policy slopes:\n');
        fprintf('      VFI mean slope: %.4f (std: %.4f)\n', mean(slope_vfi_k), std(slope_vfi_k));
        fprintf('      RL mean slope:  %.4f (std: %.4f)\n', mean(slope_rl_k), std(slope_rl_k));
        
        % PPS策略斜率
        if any(c_pps_vfi > 1e-6) || any(c_pps_rl > 1e-6)
            slope_vfi_cpps = diff(c_pps_vfi) / dk;
            slope_rl_cpps = diff(c_pps_rl) / dk;
            
            fprintf('    PPS policy slopes:\n');
            fprintf('      VFI mean slope: %.4f (std: %.4f)\n', mean(slope_vfi_cpps), std(slope_vfi_cpps));
            fprintf('      RL mean slope:  %.4f (std: %.4f)\n', mean(slope_rl_cpps), std(slope_rl_cpps));
        end
    end
    
    % 策略函数单调性检查
    vfi_k_monotonic = all(diff(k_prime_vfi) >= -1e-6);
    rl_k_monotonic = all(diff(k_prime_rl) >= -1e-6);
    
    fprintf('    Policy monotonicity:\n');
    fprintf('      VFI asset policy monotonic: %s\n', char(string(vfi_k_monotonic)));
    fprintf('      RL asset policy monotonic:  %s\n', char(string(rl_k_monotonic)));
    
    % 边界行为分析
    fprintf('    Boundary behavior:\n');
    fprintf('      At k_min (%.2f): VFI k''=%.3f, RL k''=%.3f\n', ...
        k_grid(1), k_prime_vfi(1), k_prime_rl(1));
    fprintf('      At k_max (%.2f): VFI k''=%.3f, RL k''=%.3f\n', ...
        k_grid(end), k_prime_vfi(end), k_prime_rl(end));
    
    % 储蓄率比较
    savings_rate_vfi = (k_prime_vfi - k_grid') ./ (k_grid' + 1e-6);
    savings_rate_rl = (k_prime_rl - k_grid') ./ (k_grid' + 1e-6);
    
    fprintf('    Savings rates:\n');
    fprintf('      VFI mean: %.3f (std: %.3f)\n', mean(savings_rate_vfi), std(savings_rate_vfi));
    fprintf('      RL mean:  %.3f (std: %.3f)\n', mean(savings_rate_rl), std(savings_rate_rl));
    
    fprintf('  === End of Detailed Report ===\n');
end 