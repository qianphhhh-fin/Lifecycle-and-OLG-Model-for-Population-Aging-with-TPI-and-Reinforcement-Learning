% --- 分析网格搜索VFI策略函数平滑性和稳定性脚本 ---

clc; clear; close all;

fprintf('=== 网格搜索VFI策略函数平滑性和稳定性分析 ===\n');

%% 1. 初始化参数和环境
fprintf('\n--- 1. 初始化模型参数 ---\n');
% try
    cS = main_olg_v8_utils.ParameterValues_HuggettStyle();
    
    % 手动设置PAYG替代率参数（与main_olg_v8.m保持一致）
    cS.rho_prime_payg_fixed = 0.20; % 20%的PAYG替代率
    cS.use_continuous_optimization = false; % 强制使用网格搜索
    
    fprintf('模型参数加载成功\n');
    fprintf('格点设置: nk=%d, nkpps=%d, nw=%d, PPS选择点=%d\n', ...
        cS.nk, cS.nkpps, cS.nw, cS.n_pps_choice_grid_points);
    fprintf('PAYG替代率设置: %.3f\n', cS.rho_prime_payg_fixed);
    fprintf('VFI方法: 网格搜索\n');
% catch ME
%     fprintf('错误：参数加载失败 - %s\n', ME.message);
%     return;
% end

% 设置分析参数
analysis_config = struct();
analysis_config.test_ages = [3, 7, 11]; % 测试年龄组
analysis_config.test_epsilon_idx = 2;   % 测试效率状态
analysis_config.n_stability_runs = 5;   % 稳定性测试运行次数
analysis_config.perturbation_scale = 1e-6; % 扰动规模

%% 2. 计算辅助参数
fprintf('\n--- 2. 计算辅助参数 ---\n');
% try
    % 计算人口结构 - 使用与main_olg_v8.m相同的调用序列
    popS = main_olg_v8_utils.initPopulation(cS);
    popS = main_olg_v8_utils.populationDynamics(popS, cS);
    [Z_ss, ~, ~, ~] = main_olg_v8_utils.detectSteadyStatePopulation(popS, cS);
    
    % 计算人口参数
    Z_ss_total = sum(Z_ss);
    Z_ss_norm_group = zeros(cS.aD_new,1);
    if Z_ss_total > 1e-9
        Z_ss_norm_group = Z_ss / Z_ss_total;
    end
    paramS.ageMassV = Z_ss_norm_group(:);
    paramS.mass_workers_group = sum(paramS.ageMassV(1:cS.aR_new));
    if paramS.mass_workers_group < 1e-9
        error('模型中稳态工作人口占比为零或过小。');
    end
    
    % 计算年化人口增长率
    if Z_ss_total > 1e-9 && length(popS.totalPop) > 1 && popS.totalPop(end-1) > 1e-9
        pop_growth_factor_per_group_period = popS.totalPop(end) / popS.totalPop(end-1);
        paramS.popGrowthForDebt = pop_growth_factor_per_group_period^(1/cS.yearStep) - 1;
    else
        paramS.popGrowthForDebt = cS.popGrowth_orig;
    end
    
    % 计算劳动供给和禀赋过程
    [paramS.leLogGridV, paramS.leTrProbM, paramS.leProb1V] = main_olg_v8_utils.EarningProcess_olgm(cS);
    paramS.leGridV = exp(paramS.leLogGridV(:));
    paramS.ageEffV_new = cS.ageEffV_new;
    
    eIdxM = main_olg_v8_utils.LaborEndowSimulation_olgm(cS, paramS);
    [~, L_per_capita] = main_olg_v8_utils.LaborSupply_Huggett(eIdxM, cS, paramS, paramS.ageMassV);
    
    if L_per_capita <= 0
        L_per_capita = 1e-6;
        warning('L_per_capita异常，已重置');
    end
    paramS.L_per_capita = L_per_capita;
    
    fprintf('人口和劳动效率参数计算完成\n');
    fprintf('总劳动供给: %.4f\n', L_per_capita);
    
% catch ME
%     fprintf('错误：辅助参数计算失败 - %s\n', ME.message);
%     return;
% end

%% 3. 设置测试用的外生变量
fprintf('\n--- 3. 设置测试环境 ---\n');

% 使用合理的外生价格和政策参数进行测试
test_params = struct();
test_params.K_test = 20.0;  % 测试用资本存量
test_params.tau_l_test = 0.15;  % 测试用所得税率
test_params.theta_payg_test = 0.12;  % 测试用PAYG税率
test_params.TR_total_test = 0.5;     % 测试用转移支付

% 计算对应的价格
[R_k_net_factor_test, w_gross_test] = main_olg_v8_utils.HHPrices_Huggett(test_params.K_test, paramS.L_per_capita, cS);
r_k_net_test = R_k_net_factor_test - 1;

% 计算PAYG福利
avg_worker_wage = (w_gross_test * paramS.L_per_capita) / paramS.mass_workers_group;
b_payg_test = cS.rho_prime_payg_fixed * avg_worker_wage;
bV_test = zeros(1, cS.aD_new);
if cS.aR_new < cS.aD_new
    bV_test(cS.aR_new + 1 : cS.aD_new) = b_payg_test;
end

% 设置家庭参数
paramS_test = paramS;
paramS_test.tau_l = test_params.tau_l_test;
paramS_test.theta_payg_actual_for_hh = test_params.theta_payg_test;
paramS_test.pps_tax_deferral_active = cS.pps_active;

fprintf('测试环境设置完成\n');
fprintf('  K=%.2f, 税后资本回报率=%.4f, 毛工资率=%.4f\n', test_params.K_test, r_k_net_test, w_gross_test);
fprintf('  tau_l=%.3f, theta_payg=%.3f, b_payg=%.4f\n', test_params.tau_l_test, test_params.theta_payg_test, b_payg_test);

%% 4. 运行单次完整VFI（使用网格搜索）
fprintf('\n--- 4. 运行网格搜索VFI ---\n');

fprintf('强制使用网格搜索进行VFI...\n');
vfi_start_time = tic;

% 调用VFI主函数（确保使用网格搜索版本）
[cPolM_grid, kPolM_grid, cPpsPolM_grid, valM_grid] = main_olg_v8_utils.HHSolution_VFI_Huggett(...
    R_k_net_factor_test, w_gross_test, test_params.TR_total_test, bV_test, paramS_test, cS);

vfi_time = toc(vfi_start_time);
fprintf('网格搜索VFI完成，耗时: %.2f秒\n', vfi_time);

%% 5. 策略函数平滑性分析
fprintf('\n--- 5. 策略函数平滑性分析 ---\n');

smoothness_results = []; % 初始化为空数组

for age_idx = 1:length(analysis_config.test_ages)
    a_test = analysis_config.test_ages(age_idx);
    ie_test = analysis_config.test_epsilon_idx;
    
    fprintf('\n分析年龄组 %d (约%d岁):\n', a_test, cS.physAgeV_new(a_test));
    
    % 提取该年龄组的策略函数
    c_policy_slice = squeeze(cPolM_grid(:, :, ie_test, a_test));
    k_policy_slice = squeeze(kPolM_grid(:, :, ie_test, a_test));
    cpps_policy_slice = squeeze(cPpsPolM_grid(:, :, ie_test, a_test));
    val_slice = squeeze(valM_grid(:, :, ie_test, a_test));
    
    % 计算平滑性指标
    age_smoothness = analyze_policy_smoothness(c_policy_slice, k_policy_slice, cpps_policy_slice, val_slice, cS);
    age_smoothness.age_idx = a_test;
    age_smoothness.real_age = cS.physAgeV_new(a_test);
    
    % 使用动态赋值方式
    if isempty(smoothness_results)
        smoothness_results = age_smoothness;
    else
        smoothness_results(end+1) = age_smoothness;
    end
    
    % 打印结果
    fprintf('  消费策略平滑性:\n');
    fprintf('    最大梯度: %.6f, 平均梯度: %.6f, 方差: %.6f\n', ...
        age_smoothness.c_max_gradient, age_smoothness.c_mean_gradient, age_smoothness.c_gradient_var);
    fprintf('  储蓄策略平滑性:\n');
    fprintf('    最大梯度: %.6f, 平均梯度: %.6f, 方差: %.6f\n', ...
        age_smoothness.k_max_gradient, age_smoothness.k_mean_gradient, age_smoothness.k_gradient_var);
    fprintf('  PPS策略平滑性:\n');
    fprintf('    最大梯度: %.6f, 平均梯度: %.6f, 方差: %.6f\n', ...
        age_smoothness.cpps_max_gradient, age_smoothness.cpps_mean_gradient, age_smoothness.cpps_gradient_var);
    fprintf('  值函数平滑性:\n');
    fprintf('    最大梯度: %.6f, 平均梯度: %.6f, 方差: %.6f\n', ...
        age_smoothness.val_max_gradient, age_smoothness.val_mean_gradient, age_smoothness.val_gradient_var);
end

%% 6. 数值稳定性测试
fprintf('\n--- 6. 数值稳定性测试 ---\n');

stability_results = []; % 初始化为空数组

for run_idx = 1:analysis_config.n_stability_runs
    fprintf('稳定性测试运行 %d/%d...', run_idx, analysis_config.n_stability_runs);
    
    % 对价格施加微小扰动
    R_perturbed = R_k_net_factor_test * (1 + analysis_config.perturbation_scale * randn());
    w_perturbed = w_gross_test * (1 + analysis_config.perturbation_scale * randn());
    TR_perturbed = test_params.TR_total_test * (1 + analysis_config.perturbation_scale * randn());
    
    % 运行扰动后的VFI
    [cPolM_pert, kPolM_pert, cPpsPolM_pert, valM_pert] = main_olg_v8_utils.HHSolution_VFI_Huggett(...
        R_perturbed, w_perturbed, TR_perturbed, bV_test, paramS_test, cS);
    
    % 计算差异
    c_diff = abs(cPolM_grid - cPolM_pert);
    k_diff = abs(kPolM_grid - kPolM_pert);
    cpps_diff = abs(cPpsPolM_grid - cPpsPolM_pert);
    val_diff = abs(valM_grid - valM_pert);
    
    % 创建当前运行的结果结构体
    run_result = struct();
    run_result.c_max_diff = max(c_diff(:));
    run_result.c_mean_diff = mean(c_diff(:));
    run_result.k_max_diff = max(k_diff(:));
    run_result.k_mean_diff = mean(k_diff(:));
    run_result.cpps_max_diff = max(cpps_diff(:));
    run_result.cpps_mean_diff = mean(cpps_diff(:));
    run_result.val_max_diff = max(val_diff(:));
    run_result.val_mean_diff = mean(val_diff(:));
    
    % 动态赋值
    if isempty(stability_results)
        stability_results = run_result;
    else
        stability_results(end+1) = run_result;
    end
    
    fprintf(' 完成\n');
end

% 计算稳定性统计
c_stability_ratio = mean([stability_results.c_max_diff]) / analysis_config.perturbation_scale;
k_stability_ratio = mean([stability_results.k_max_diff]) / analysis_config.perturbation_scale;
cpps_stability_ratio = mean([stability_results.cpps_max_diff]) / analysis_config.perturbation_scale;
val_stability_ratio = mean([stability_results.val_max_diff]) / analysis_config.perturbation_scale;

fprintf('\n稳定性分析结果:\n');
fprintf('  消费策略稳定性比率: %.2f (理想值接近1)\n', c_stability_ratio);
fprintf('  储蓄策略稳定性比率: %.2f\n', k_stability_ratio);
fprintf('  PPS策略稳定性比率: %.2f\n', cpps_stability_ratio);
fprintf('  值函数稳定性比率: %.2f\n', val_stability_ratio);

if max([c_stability_ratio, k_stability_ratio, cpps_stability_ratio]) > 100
    fprintf('  警告: 检测到数值不稳定性！策略函数对微小扰动过度敏感。\n');
elseif max([c_stability_ratio, k_stability_ratio, cpps_stability_ratio]) > 10
    fprintf('  注意: 策略函数对扰动较为敏感，可能影响外层均衡收敛。\n');
else
    fprintf('  良好: 策略函数数值稳定。\n');
end

%% 7. 绘制策略函数
fprintf('\n--- 7. 绘制策略函数图形 ---\n');

create_policy_plots(cPolM_grid, kPolM_grid, cPpsPolM_grid, valM_grid, cS, analysis_config.test_ages);

%% 8. 边界行为分析
fprintf('\n--- 8. 边界行为分析 ---\n');

boundary_analysis = analyze_boundary_behavior(cPolM_grid, kPolM_grid, cPpsPolM_grid, cS);

fprintf('边界行为分析:\n');
fprintf('  消费策略触及下界比例: %.2f%%\n', boundary_analysis.c_floor_fraction * 100);
fprintf('  储蓄策略触及下界比例: %.2f%%\n', boundary_analysis.k_min_fraction * 100);
fprintf('  储蓄策略触及上界比例: %.2f%%\n', boundary_analysis.k_max_fraction * 100);
fprintf('  PPS缴费为零比例: %.2f%%\n', boundary_analysis.cpps_zero_fraction * 100);

%% 9. 收敛性诊断
fprintf('\n--- 9. 收敛性诊断建议 ---\n');

% 基于分析结果给出诊断
provide_convergence_diagnostics(smoothness_results, stability_results, boundary_analysis, c_stability_ratio, k_stability_ratio);

fprintf('\n=== 分析完成 ===\n');

%% ===== 辅助函数定义 =====

function smoothness = analyze_policy_smoothness(c_pol, k_pol, cpps_pol, val_pol, cS)
    % 分析策略函数的平滑性
    
    smoothness = struct();
    
    % 计算梯度（使用有限差分）
    if size(c_pol, 1) > 1
        % 沿k方向的梯度
        c_grad_k = diff(c_pol, 1, 1);
        k_grad_k = diff(k_pol, 1, 1);
        cpps_grad_k = diff(cpps_pol, 1, 1);
        val_grad_k = diff(val_pol, 1, 1);
        
        smoothness.c_max_gradient = max(abs(c_grad_k(:)));
        smoothness.c_mean_gradient = mean(abs(c_grad_k(:)));
        smoothness.c_gradient_var = var(c_grad_k(:));
        
        smoothness.k_max_gradient = max(abs(k_grad_k(:)));
        smoothness.k_mean_gradient = mean(abs(k_grad_k(:)));
        smoothness.k_gradient_var = var(k_grad_k(:));
        
        smoothness.cpps_max_gradient = max(abs(cpps_grad_k(:)));
        smoothness.cpps_mean_gradient = mean(abs(cpps_grad_k(:)));
        smoothness.cpps_gradient_var = var(cpps_grad_k(:));
        
        smoothness.val_max_gradient = max(abs(val_grad_k(:)));
        smoothness.val_mean_gradient = mean(abs(val_grad_k(:)));
        smoothness.val_gradient_var = var(val_grad_k(:));
    else
        smoothness.c_max_gradient = 0;
        smoothness.c_mean_gradient = 0;
        smoothness.c_gradient_var = 0;
        smoothness.k_max_gradient = 0;
        smoothness.k_mean_gradient = 0;
        smoothness.k_gradient_var = 0;
        smoothness.cpps_max_gradient = 0;
        smoothness.cpps_mean_gradient = 0;
        smoothness.cpps_gradient_var = 0;
        smoothness.val_max_gradient = 0;
        smoothness.val_mean_gradient = 0;
        smoothness.val_gradient_var = 0;
    end
end

function create_policy_plots(cPolM, kPolM, cPpsPolM, valM, cS, test_ages)
    % 创建策略函数图形
    
    figure('Position', [100, 100, 1200, 800]);
    
    for age_idx = 1:length(test_ages)
        a_test = test_ages(age_idx);
        ie_test = 2; % 中等效率状态
        
        % 选择中间的kpps状态进行绘图
        ikpps_mid = max(1, round(cS.nkpps / 2));
        
        % 提取策略函数切片
        c_slice = squeeze(cPolM(:, ikpps_mid, ie_test, a_test));
        k_slice = squeeze(kPolM(:, ikpps_mid, ie_test, a_test));
        cpps_slice = squeeze(cPpsPolM(:, ikpps_mid, ie_test, a_test));
        val_slice = squeeze(valM(:, ikpps_mid, ie_test, a_test));
        
        % 消费策略
        subplot(4, length(test_ages), age_idx);
        plot(cS.kGridV, c_slice, 'b-', 'LineWidth', 1.5);
        title(sprintf('消费策略 (年龄%d)', cS.physAgeV_new(a_test)));
        xlabel('当前资产k'); ylabel('消费c');
        grid on;
        
        % 储蓄策略
        subplot(4, length(test_ages), length(test_ages) + age_idx);
        plot(cS.kGridV, k_slice, 'r-', 'LineWidth', 1.5);
        hold on;
        plot(cS.kGridV, cS.kGridV, 'k--', 'LineWidth', 0.8);
        title(sprintf('储蓄策略 (年龄%d)', cS.physAgeV_new(a_test)));
        xlabel('当前资产k'); ylabel('下期资产k''');
        legend('最优储蓄', '45度线', 'Location', 'best');
        grid on;
        
        % PPS缴费策略
        subplot(4, length(test_ages), 2*length(test_ages) + age_idx);
        plot(cS.kGridV, cpps_slice, 'g-', 'LineWidth', 1.5);
        title(sprintf('PPS缴费策略 (年龄%d)', cS.physAgeV_new(a_test)));
        xlabel('当前资产k'); ylabel('PPS缴费c_{pps}');
        grid on;
        
        % 值函数
        subplot(4, length(test_ages), 3*length(test_ages) + age_idx);
        plot(cS.kGridV, val_slice, 'm-', 'LineWidth', 1.5);
        title(sprintf('值函数 (年龄%d)', cS.physAgeV_new(a_test)));
        xlabel('当前资产k'); ylabel('值函数V');
        grid on;
    end
    
    sgtitle('网格搜索VFI策略函数 (效率状态2, 中等PPS资产水平)');
end

function boundary_analysis = analyze_boundary_behavior(cPolM, kPolM, cPpsPolM, cS)
    % 分析边界行为
    
    total_states = numel(cPolM);
    
    boundary_analysis = struct();
    boundary_analysis.c_floor_fraction = sum(cPolM(:) <= cS.cFloor + 1e-6) / total_states;
    boundary_analysis.k_min_fraction = sum(kPolM(:) <= cS.kMin + 1e-6) / total_states;
    boundary_analysis.k_max_fraction = sum(kPolM(:) >= cS.kMax - 1e-6) / total_states;
    boundary_analysis.cpps_zero_fraction = sum(cPpsPolM(:) <= 1e-6) / total_states;
end

function provide_convergence_diagnostics(smoothness_results, stability_results, boundary_analysis, c_stability_ratio, k_stability_ratio)
    % 提供收敛性诊断建议
    
    fprintf('基于分析结果的诊断建议:\n\n');
    
    % 平滑性诊断
    avg_c_gradient_var = mean([smoothness_results.c_gradient_var]);
    avg_k_gradient_var = mean([smoothness_results.k_gradient_var]);
    
    if avg_c_gradient_var > 1e-3 || avg_k_gradient_var > 1e-3
        fprintf('1. 平滑性问题:\n');
        fprintf('   - 策略函数梯度方差较大，可能导致聚合行为不稳定\n');
        fprintf('   - 建议: 增加阻尼系数 (damp_K_v5, damp_tau_l_v5 → 0.1-0.2)\n\n');
    else
        fprintf('1. 平滑性良好:\n');
        fprintf('   - 策略函数足够平滑，有利于外层收敛\n\n');
    end
    
    % 稳定性诊断
    if c_stability_ratio > 10 || k_stability_ratio > 10
        fprintf('2. 数值稳定性问题:\n');
        fprintf('   - 策略函数对微小扰动过度敏感\n');
        fprintf('   - 建议: 放宽外层容忍度 (tol_K_tau_l → 5e-4, gbc_tol → 3e-3)\n\n');
    else
        fprintf('2. 数值稳定性良好:\n');
        fprintf('   - 策略函数对扰动的反应适中\n\n');
    end
    
    % 边界行为诊断
    if boundary_analysis.k_min_fraction > 0.1 || boundary_analysis.k_max_fraction > 0.05
        fprintf('3. 边界行为警告:\n');
        fprintf('   - 过多状态触及资产边界 (%.1f%% 下界, %.1f%% 上界)\n', ...
            boundary_analysis.k_min_fraction*100, boundary_analysis.k_max_fraction*100);
        fprintf('   - 建议: 扩大资产网格范围或调整网格密度\n\n');
    else
        fprintf('3. 边界行为正常:\n');
        fprintf('   - 资产网格设置合理\n\n');
    end
    
    % 综合建议
    fprintf('4. 综合建议:\n');
    if c_stability_ratio > 10
        fprintf('   - 网格搜索法在当前设置下表现出一定的数值不稳定性。');
        fprintf('   - 建议: 检查网格密度，或在VFI中使用更高精度的插值。');
    else
        fprintf('   - 网格搜索法总体表现良好。');
    end
end 