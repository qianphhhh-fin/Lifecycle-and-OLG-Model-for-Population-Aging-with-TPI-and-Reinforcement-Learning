% 简化版本：比较VFI和RL（SAC Agent）的优化结果
% 使用相同的宏观和微观参数比较两种方法

clc;
clear;
close all;
fprintf('=== 简化版本：比较 VFI 和 RL（SAC Agent）===\n');

%% 1. 初始化参数
fprintf('\n--- 1. 初始化参数 ---\n');
cS = main_olg_v8_utils.ParameterValues_HuggettStyle();

% 确保网格已生成
if ~isfield(cS, 'kppsGridV') || isempty(cS.kppsGridV)
    cS = main_olg_v8_utils.UpdateGrids(cS);
end

% 设置固定宏观参数用于公平比较
M_fixed.R_k_net_factor = 1.03;
M_fixed.w_gross = 2.0;
M_fixed.TR_total = 0.1;
M_fixed.b_payg_avg_retiree = 0.4;
M_fixed.tau_l = 0.15;
M_fixed.theta_payg_actual = 0.12;
M_fixed.b_payg_avg_for_obs = M_fixed.b_payg_avg_retiree;

fprintf('固定宏观参数:\n');
fprintf('  资本净回报率: %.3f\n', M_fixed.R_k_net_factor);
fprintf('  总工资: %.3f\n', M_fixed.w_gross);
fprintf('  劳动税率: %.3f\n', M_fixed.tau_l);

%% 2. 准备共同参数
[leLogGridV, leTrProbM, leProb1V] = main_olg_v8_utils.EarningProcess_olgm(cS);
leGridV = exp(leLogGridV(:));

paramS_common = struct();
paramS_common.leLogGridV = leLogGridV;
paramS_common.leTrProbM = leTrProbM;
paramS_common.leProb1V = leProb1V;
paramS_common.leGridV = leGridV;
paramS_common.ageEffV_new = cS.ageEffV_new;
paramS_common.tau_l = M_fixed.tau_l;
paramS_common.theta_payg_actual_for_hh = M_fixed.theta_payg_actual;
paramS_common.pps_tax_deferral_active = cS.pps_active;

% 构建PAYG福利向量
bV_payg = zeros(1, cS.aD_new);
if cS.aR_new < cS.aD_new
    bV_payg(cS.aR_new + 1 : cS.aD_new) = M_fixed.b_payg_avg_retiree;
end

%% 3. 运行VFI方法
fprintf('\n--- 2. 运行VFI方法 ---\n');
tic;
[cPolM_vfi, kPolM_vfi, cPpsPolM_vfi, VPolM_vfi] = ...
    main_olg_v8_utils.HHSolution_VFI_Huggett(...
    M_fixed.R_k_net_factor, M_fixed.w_gross, M_fixed.TR_total, bV_payg, paramS_common, cS);
vfi_time = toc;

fprintf('VFI完成: 耗时=%.2f秒, 策略矩阵尺寸=%s\n', vfi_time, mat2str(size(cPolM_vfi)));

%% 4. 尝试加载RL方法
fprintf('\n--- 3. 尝试加载RL方法 ---\n');
rl_available = false;
sac_file = 'final_SAC_Agent_OLG_R2024b.mat';

if exist(sac_file, 'file')
    try
        load(sac_file, 'finalAgent', 'obsInfo', 'actInfo');
        if exist('finalAgent', 'var') && isa(finalAgent, 'rl.agent.AbstractAgent')
            fprintf('SAC agent加载成功\n');
            
            % 转换策略
            tic;
            [cPolM_rl, kPolM_rl, cPpsPolM_rl, ~] = convert_SAC_to_PolicyFunctions(...
                finalAgent, cS, paramS_common, M_fixed);
            rl_time = toc;
            
            fprintf('RL策略转换完成，耗时=%.2f秒\n', rl_time);
            rl_available = true;
        end
    catch ME
        fprintf('RL加载失败: %s\n', ME.message);
    end
else
    fprintf('未找到SAC agent文件: %s\n', sac_file);
end

%% 5. 比较分析（如果RL可用）
if rl_available
    fprintf('\n--- 4. 比较分析 ---\n');
    
    % 选择一些代表性状态点进行比较
    nTest = 20;
    utility_diff_points = zeros(nTest, 1);
    
    fprintf('在%d个代表性状态点比较策略...\n', nTest);
    
    for i = 1:nTest
        % 随机选择状态
        k_idx = randi(cS.nk);
        kpps_idx = randi(cS.nkpps);
        eps_idx = randi(cS.nw);
        age_idx = randi(cS.aD_new);
        
        % 获取策略值
        c_vfi = cPolM_vfi(k_idx, kpps_idx, eps_idx, age_idx);
        c_rl = cPolM_rl(k_idx, kpps_idx, eps_idx, age_idx);
        
        % 计算当期效用差异
        [~, u_vfi] = main_olg_v8_utils.CES_utility(c_vfi, cS.sigma, cS);
        [~, u_rl] = main_olg_v8_utils.CES_utility(c_rl, cS.sigma, cS);
        
        utility_diff_points(i) = u_rl - u_vfi;
    end
    
    % 统计结果
    mean_diff = mean(utility_diff_points);
    std_diff = std(utility_diff_points);
    
    fprintf('\n当期效用差异统计 (RL - VFI):\n');
    fprintf('  均值: %.6f\n', mean_diff);
    fprintf('  标准差: %.6f\n', std_diff);
    fprintf('  改进的点数: %d/%d (%.1f%%)\n', ...
        sum(utility_diff_points > 0), nTest, 100*sum(utility_diff_points > 0)/nTest);
    
    if abs(mean_diff) < 1e-4
        fprintf('  >>> 两种方法表现相当\n');
    elseif mean_diff > 0
        fprintf('  >>> RL方法平均表现更好\n');
    else
        fprintf('  >>> VFI方法平均表现更好\n');
    end
    
    % 简单可视化
    figure('Name', '策略比较', 'Position', [100, 100, 1000, 400]);
    
    subplot(1, 3, 1);
    histogram(utility_diff_points, 8);
    xlabel('效用差异 (RL - VFI)');
    ylabel('频数');
    title('效用差异分布');
    xline(0, 'r--', 'LineWidth', 2);
    grid on;
    
    subplot(1, 3, 2);
    plot(1:nTest, utility_diff_points, 'bo-');
    xlabel('测试点');
    ylabel('效用差异');
    title('各测试点差异');
    yline(0, 'r--');
    yline(mean_diff, 'g-', 'LineWidth', 2);
    grid on;
    
    subplot(1, 3, 3);
    bar([sum(utility_diff_points < 0), sum(utility_diff_points > 0)]);
    set(gca, 'XTickLabel', {'VFI更好', 'RL更好'});
    ylabel('点数');
    title('优劣分布');
    
    % 比较计算效率
    fprintf('\n计算效率比较:\n');
    fprintf('  VFI时间: %.2f秒\n', vfi_time);
    fprintf('  RL转换时间: %.2f秒\n', rl_time);
    fprintf('  速度比 (VFI/RL): %.2fx\n', vfi_time/rl_time);
    
    % 保存简单结果
    simple_results = struct();
    simple_results.mean_utility_diff = mean_diff;
    simple_results.std_utility_diff = std_diff;
    simple_results.rl_better_count = sum(utility_diff_points > 0);
    simple_results.total_tests = nTest;
    simple_results.vfi_time = vfi_time;
    simple_results.rl_time = rl_time;
    simple_results.M_fixed = M_fixed;
    
    save('simple_comparison_results.mat', 'simple_results');
    fprintf('\n结果已保存到 simple_comparison_results.mat\n');
    
else
    fprintf('\nRL方法不可用，仅显示VFI结果\n');
    fprintf('VFI求解时间: %.2f秒\n', vfi_time);
    fprintf('请先运行 main_olg_v8_SAC.m 训练SAC agent\n');
end

fprintf('\n=== 简化比较完成 ===\n'); 