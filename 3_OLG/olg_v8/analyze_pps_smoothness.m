%% PPS缴费策略平滑性分析脚本
% 该脚本专门检查VFI求解后的PPS缴费策略是否平滑
% 基于fmincon_objective_helper_proportional（无梯度版本）

clear; clc; close all;

fprintf('=== PPS缴费策略平滑性分析 ===\n\n');

%% 1. 参数设置和初始化
fprintf('--- 1. 加载参数和设置 ---\n');
try
    cS = main_olg_v8_utils.ParameterValues_HuggettStyle();
    fprintf('✓ 参数加载成功\n');
    fprintf('  网格设置: nk=%d, nkpps=%d, nw=%d\n', cS.nk, cS.nkpps, cS.nw);
catch ME
    fprintf('✗ 参数加载失败: %s\n', ME.message);
    return;
end

% 初始化完整的paramS结构体（仿照main_olg_v8.m）
fprintf('--- 1.1 初始化人口参数 ---\n');
paramS = struct();

% 模拟人口动态
popS = main_olg_v8_utils.initPopulation(cS);
popS = main_olg_v8_utils.populationDynamics(popS, cS);
[Z_ss, ~, ~, ~] = main_olg_v8_utils.detectSteadyStatePopulation(popS, cS);
paramS.Z_ss_counts = Z_ss;
Z_ss_total = sum(Z_ss);
Z_ss_norm_group = zeros(cS.aD_new,1);
if Z_ss_total > 1e-9
    Z_ss_norm_group = Z_ss / Z_ss_total;
end
paramS.ageMassV = Z_ss_norm_group(:);
paramS.mass_workers_group = sum(paramS.ageMassV(1:cS.aR_new));

% 年度人口分布
Z_ss_norm_annual = zeros(cS.aD_orig,1);
if Z_ss_total > 1e-9
    for a_new_map_idx = 1:cS.aD_new
        annual_indices_in_group = cS.physAgeMap{a_new_map_idx};
        group_mass_fraction = Z_ss_norm_group(a_new_map_idx);
        num_years_in_this_group = length(annual_indices_in_group);
        if num_years_in_this_group > 0
            mass_per_year_in_group = group_mass_fraction / num_years_in_this_group;
            Z_ss_norm_annual(annual_indices_in_group) = mass_per_year_in_group;
        end
    end
    if sum(Z_ss_norm_annual) > 1e-9 && abs(sum(Z_ss_norm_annual) - 1.0) > 1e-6
        Z_ss_norm_annual = Z_ss_norm_annual / sum(Z_ss_norm_annual);
    elseif sum(Z_ss_norm_annual) < 1e-9
        Z_ss_norm_annual(:) = 1/cS.aD_orig;
    end
else
    Z_ss_norm_annual(:) = 1/cS.aD_orig;
end
paramS.Z_ss_norm_annual = Z_ss_norm_annual;

if Z_ss_total > 1e-9 && length(popS.totalPop) > 1 && popS.totalPop(end-1) > 1e-9
    pop_growth_factor_per_group_period = popS.totalPop(end) / popS.totalPop(end-1);
    paramS.popGrowthForDebt = pop_growth_factor_per_group_period^(1/cS.yearStep) - 1;
else
    paramS.popGrowthForDebt = cS.popGrowth_orig;
end

fprintf('--- 1.2 初始化劳动禀赋参数 ---\n');
% 预计算劳动供给和禀赋过程
[paramS.leLogGridV, paramS.leTrProbM, paramS.leProb1V] = main_olg_v8_utils.EarningProcess_olgm(cS);
paramS.leGridV = exp(paramS.leLogGridV(:));
paramS.ageEffV_new = cS.ageEffV_new;

% 模拟劳动禀赋过程
eIdxM = main_olg_v8_utils.LaborEndowSimulation_olgm(cS, paramS);
[~, L_per_capita] = main_olg_v8_utils.LaborSupply_Huggett(eIdxM, cS, paramS, paramS.ageMassV);
if L_per_capita <= 0
    L_per_capita = 1e-6;
    warning('L_per_capita 为零或负，已重置为1e-6。');
end
paramS.L_per_capita = L_per_capita;

fprintf('✓ 劳动禀赋参数初始化完成\n');
fprintf('  总劳动供给 (L, 效率单位): %.4f\n', L_per_capita);
fprintf('  稳态工人占比: %.4f\n', paramS.mass_workers_group);

% 设置必要的宏观变量（简化版）
cS.rho_prime_payg_fixed = 0.20;
fixedMacro = struct();
fixedMacro.R_k_net_factor_hh = 1.035;
fixedMacro.MPL_gross = 1.0;
fixedMacro.TR_total = 0.1;

% 确保使用无梯度版本
cS.use_continuous_optimization = true;
cS.enable_fmincon_gradients = false;  % 强制使用无梯度版本
fprintf('✓ 设置为无梯度版本优化\n');

%% 2. 运行单次VFI以获取策略函数
fprintf('\n--- 2. 运行VFI求解 ---\n');
tic;

% 选择一个中间年龄组进行分析
test_age = round(cS.aD_new * 0.1);  % 约60%生命周期位置
fprintf('分析年龄组: %d (共%d个年龄组)\n', test_age, cS.aD_new);

% 创建简化的期望价值函数（用于测试）
vPrime_test = ones(cS.nk, cS.nkpps, cS.nw) * 50;  % 简化的下期价值函数
for ik = 1:cS.nk
    for ikpps = 1:cS.nkpps
        for iw = 1:cS.nw
            vPrime_test(ik, ikpps, iw) = 50 + log(cS.kGridV(ik) + 1) + ...
                sqrt(cS.kppsGridV(ikpps) + 1) + 0.1 * iw;
        end
    end
end

% 设置宏观变量
R_k_test = fixedMacro.R_k_net_factor_hh;
w_gross_test = fixedMacro.MPL_gross;
TR_total_test = fixedMacro.TR_total;
b_age_test = 0.05;

% 设置年龄相关参数（现在使用完整的paramS）
paramS_test = paramS;  % 使用完整的paramS
paramS_test.tau_l = 0.10;
paramS_test.theta_payg_actual_for_hh = 0.08;
paramS_test.pps_tax_deferral_active = true;
epsilon_test = linspace(0.5, 1.5, cS.nw);

% 运行VFI
try
    [cPol, kPol, cPpsPol, val] = main_olg_v8_utils.HHSolutionByAge_VFI_Huggett_v8(...
        test_age, vPrime_test, R_k_test, w_gross_test, TR_total_test, b_age_test, ...
        paramS_test, cS, epsilon_test);
    
    vfi_time = toc;
    fprintf('✓ VFI求解完成，耗时: %.2f秒\n', vfi_time);
catch ME
    fprintf('✗ VFI求解失败: %s\n', ME.message);
    return;
end

%% 3. PPS策略平滑性分析
fprintf('\n--- 3. PPS策略平滑性分析 ---\n');

% 基本统计
pps_values = cPpsPol(:);
pps_nonzero = pps_values(pps_values > 1e-10);

fprintf('PPS缴费策略基本统计:\n');
fprintf('  最小值: %.6f\n', min(pps_values));
fprintf('  最大值: %.6f\n', max(pps_values));
fprintf('  平均值: %.6f\n', mean(pps_values));
fprintf('  非零比例: %.1f%%\n', length(pps_nonzero)/length(pps_values)*100);

%% 4. 跨状态变量的平滑性检查

%% 4.1 固定收入冲击，分析跨资产的平滑性
fprintf('\n--- 4.1 跨资产状态的平滑性分析 ---\n');

% 选择中位收入冲击
ie_mid = round(cS.nw/2);
epsilon_mid = epsilon_test(ie_mid);

% 固定PPS资产水平，分析跨普通资产的平滑性
ikpps_levels = [1, round(cS.nkpps/3), round(2*cS.nkpps/3), cS.nkpps];
smoothness_metrics = struct();

for i = 1:length(ikpps_levels)
    ikpps = ikpps_levels(i);
    kpps_level = cS.kppsGridV(ikpps);
    
    % 提取跨k的PPS策略
    pps_across_k = cPpsPol(:, ikpps, ie_mid);
    k_grid = cS.kGridV;
    
    % 计算平滑性指标
    % 1. 一阶差分（近似导数）
    if length(k_grid) > 1
        dk = diff(k_grid);
        dpps_dk = diff(pps_across_k) ./ dk;
        
        % 2. 二阶差分（近似二阶导数）
        if length(dpps_dk) > 1
            d2pps_dk2 = diff(dpps_dk) ./ dk(1:end-1);
            
            % 平滑性指标
            smoothness_metrics(i).kpps_level = kpps_level;
            smoothness_metrics(i).max_first_deriv = max(abs(dpps_dk));
            smoothness_metrics(i).std_first_deriv = std(dpps_dk);
            smoothness_metrics(i).max_second_deriv = max(abs(d2pps_dk2));
            smoothness_metrics(i).std_second_deriv = std(d2pps_dk2);
            
            % 检查跳跃点（一阶导数的大变化）
            large_jumps = find(abs(dpps_dk) > 3*std(dpps_dk));
            smoothness_metrics(i).num_jumps = length(large_jumps);
            smoothness_metrics(i).jump_locations = k_grid(large_jumps+1);
            
            fprintf('PPS资产水平 %.2f:\n', kpps_level);
            fprintf('  最大一阶导数: %.4f\n', smoothness_metrics(i).max_first_deriv);
            fprintf('  一阶导数标准差: %.4f\n', smoothness_metrics(i).std_first_deriv);
            fprintf('  最大二阶导数: %.4f\n', smoothness_metrics(i).max_second_deriv);
            fprintf('  检测到跳跃点数量: %d\n', smoothness_metrics(i).num_jumps);
            if smoothness_metrics(i).num_jumps > 0
                fprintf('  跳跃位置: [%.2f]\n', smoothness_metrics(i).jump_locations);
            end
        end
    end
end

%% 4.2 固定普通资产，分析跨PPS资产的平滑性
fprintf('\n--- 4.2 跨PPS资产状态的平滑性分析 ---\n');

ik_levels = [1, round(cS.nk/3), round(2*cS.nk/3), cS.nk];

for i = 1:length(ik_levels)
    ik = ik_levels(i);
    k_level = cS.kGridV(ik);
    
    % 提取跨k_pps的PPS策略
    pps_across_kpps = squeeze(cPpsPol(ik, :, ie_mid));
    kpps_grid = cS.kppsGridV;
    
    % 计算平滑性指标
    if length(kpps_grid) > 1
        dkpps = diff(kpps_grid);
        dpps_dkpps = diff(pps_across_kpps) ./ dkpps;
        
        % 检查单调性
        is_monotonic = all(dpps_dkpps(:) >= -1e-10) || all(dpps_dkpps(:) <= 1e-10);
        
        fprintf('普通资产水平 %.2f:\n', k_level);
        fprintf('  最大PPS-资产导数: %.4f\n', max(abs(dpps_dkpps)));
        fprintf('  PPS-资产导数标准差: %.4f\n', std(dpps_dkpps));
        fprintf('  是否单调: %s\n', mat2str(is_monotonic));
        
        % 检查非单调性点
        if ~is_monotonic
            sign_changes = find(diff(sign(dpps_dkpps)) ~= 0);
            fprintf('  非单调点数量: %d\n', length(sign_changes));
        end
    end
end

%% 5. 约束激活分析
fprintf('\n--- 5. 约束激活分析 ---\n');

% 检查PPS缴费约束的激活情况
max_pps_grid = zeros(cS.nk, cS.nkpps, cS.nw);
for ik = 1:cS.nk
    for ikpps = 1:cS.nkpps
        for ie = 1:cS.nw
            k_state = cS.kGridV(ik);
            k_pps_state = cS.kppsGridV(ikpps);
            epsilon_state = epsilon_test(ie);
            
            % 计算该状态下的最大可缴费额
            if test_age <= cS.aR_new  % 工作期
                [resources_after_0pps, labor_income_gross, ~] = main_olg_v8_utils.HHIncome_Huggett(...
                    k_state, R_k_test, w_gross_test, TR_total_test, b_age_test, ...
                    0, test_age, paramS_test, cS, epsilon_state);
                
                if isfield(cS, 'pps_contrib_rate_for_hh') && cS.pps_contrib_rate_for_hh > 0
                    max_pps_grid(ik, ikpps, ie) = labor_income_gross * cS.pps_contrib_rate_for_hh;
                else
                    max_pps_grid(ik, ikpps, ie) = labor_income_gross * 0.15;  % 假设15%上限
                end
            else
                max_pps_grid(ik, ikpps, ie) = 0;  % 退休期不缴费
            end
        end
    end
end

% 计算约束激活情况
constraint_active = abs(cPpsPol - max_pps_grid) < 1e-8 & max_pps_grid > 1e-8;
zero_contrib = cPpsPol < 1e-8;

fprintf('约束激活统计:\n');
fprintf('  上限约束激活比例: %.1f%%\n', sum(constraint_active(:))/numel(constraint_active)*100);
fprintf('  零缴费比例: %.1f%%\n', sum(zero_contrib(:))/numel(zero_contrib)*100);
fprintf('  内点解比例: %.1f%%\n', sum(~constraint_active(:) & ~zero_contrib(:))/numel(constraint_active)*100);

%% 6. 可视化分析
fprintf('\n--- 6. 生成可视化图表 ---\n');

% 创建图形窗口
figure('Position', [100, 100, 1200, 800]);

% 子图1: 跨普通资产的PPS策略（固定PPS资产和收入冲击）
subplot(2, 3, 1);
ikpps_mid = round(cS.nkpps/2);
plot(cS.kGridV, cPpsPol(:, ikpps_mid, ie_mid), 'b-', 'LineWidth', 2);
xlabel('普通资产 k');
ylabel('PPS缴费');
title(sprintf('PPS策略 vs 普通资产\n(k_{pps}=%.2f, ε=%.2f)', cS.kppsGridV(ikpps_mid), epsilon_mid));
grid on;

% 子图2: 跨PPS资产的PPS策略（固定普通资产和收入冲击）
subplot(2, 3, 2);
ik_mid = round(cS.nk/2);
plot(cS.kppsGridV, squeeze(cPpsPol(ik_mid, :, ie_mid)), 'r-', 'LineWidth', 2);
xlabel('PPS资产 k_{pps}');
ylabel('PPS缴费');
title(sprintf('PPS策略 vs PPS资产\n(k=%.2f, ε=%.2f)', cS.kGridV(ik_mid), epsilon_mid));
grid on;

% 子图3: 跨收入冲击的PPS策略
subplot(2, 3, 3);
pps_across_epsilon = squeeze(cPpsPol(ik_mid, ikpps_mid, :));
plot(epsilon_test, pps_across_epsilon, 'g-', 'LineWidth', 2);
xlabel('收入冲击 ε');
ylabel('PPS缴费');
title(sprintf('PPS策略 vs 收入冲击\n(k=%.2f, k_{pps}=%.2f)', cS.kGridV(ik_mid), cS.kppsGridV(ikpps_mid)));
grid on;

% 子图4: PPS策略热图（k vs k_pps，固定收入冲击）
subplot(2, 3, 4);
imagesc(cS.kppsGridV, cS.kGridV, cPpsPol(:, :, ie_mid));
colorbar;
xlabel('PPS资产 k_{pps}');
ylabel('普通资产 k');
title('PPS策略热图');
axis xy;

% 子图5: 约束激活区域
subplot(2, 3, 5);
constraint_viz = double(constraint_active(:, :, ie_mid)) + 2*double(zero_contrib(:, :, ie_mid));
imagesc(cS.kppsGridV, cS.kGridV, constraint_viz);
colorbar;
cmap = [1 1 1; 0 0 1; 1 0 0];  % 白色=内点，蓝色=上限约束，红色=零缴费
colormap(gca, cmap);
xlabel('PPS资产 k_{pps}');
ylabel('普通资产 k');
title('约束激活区域');
axis xy;

% 子图6: 一阶导数（平滑性指标）
subplot(2, 3, 6);
if exist('dpps_dk', 'var')
    plot(k_grid(2:end), dpps_dk, 'k-', 'LineWidth', 1.5);
    xlabel('普通资产 k');
    ylabel('dPPS/dk');
    title('PPS策略一阶导数');
    grid on;
    
    % 标记跳跃点
    if exist('large_jumps', 'var') && ~isempty(large_jumps)
        hold on;
        scatter(k_grid(large_jumps+1), dpps_dk(large_jumps), 50, 'r', 'filled');
        legend('一阶导数', '跳跃点');
        hold off;
    end
end

sgtitle('PPS缴费策略平滑性分析', 'FontSize', 14, 'FontWeight', 'bold');

%% 7. 不平滑原因诊断
fprintf('\n--- 7. 不平滑原因诊断 ---\n');

% 检查可能的不平滑原因
reasons = {};

% 1. 约束激活导致的不平滑
if sum(constraint_active(:)) > 0.1 * numel(constraint_active)
    reasons{end+1} = '约束激活: 上限约束激活比例过高';
end

% 2. 网格效应
if cS.nk < 50 || cS.nkpps < 30
    reasons{end+1} = '网格密度不足: 可能导致策略函数不平滑';
end

% 3. 优化精度问题
total_smoothness_issues = 0;
for i = 1:length(smoothness_metrics)
    if isfield(smoothness_metrics(i), 'num_jumps') && smoothness_metrics(i).num_jumps > 0
        total_smoothness_issues = total_smoothness_issues + smoothness_metrics(i).num_jumps;
    end
end

if total_smoothness_issues > 5
    reasons{end+1} = '优化精度: 检测到多个跳跃点，可能是数值优化精度不足';
end

% 4. 参数设置问题
if cS.sigma > 5
    reasons{end+1} = '风险厌恶系数过高: 可能导致策略函数剧烈变化';
end

fprintf('可能的不平滑原因:\n');
if isempty(reasons)
    fprintf('  未检测到明显的不平滑问题\n');
else
    for i = 1:length(reasons)
        fprintf('  %d. %s\n', i, reasons{i});
    end
end

%% 8. 改进建议
fprintf('\n--- 8. 改进建议 ---\n');

suggestions = {};

if sum(constraint_active(:)) > 0.2 * numel(constraint_active)
    suggestions{end+1} = '增加网格密度，特别是在约束边界附近';
    suggestions{end+1} = '考虑使用自适应网格';
end

if total_smoothness_issues > 3
    suggestions{end+1} = '提高fmincon的收敛容忍度';
    suggestions{end+1} = '使用更多的初始猜测点';
    suggestions{end+1} = '考虑使用有梯度的目标函数版本';
end

if cS.nk < 100
    suggestions{end+1} = sprintf('增加普通资产网格点数 (当前: %d, 建议: >100)', cS.nk);
end

if cS.nkpps < 50
    suggestions{end+1} = sprintf('增加PPS资产网格点数 (当前: %d, 建议: >50)', cS.nkpps);
end

fprintf('改进建议:\n');
if isempty(suggestions)
    fprintf('  当前设置看起来合理\n');
else
    for i = 1:length(suggestions)
        fprintf('  %d. %s\n', i, suggestions{i});
    end
end

%% 9. 保存结果
% fprintf('\n--- 9. 保存分析结果 ---\n');
% 
% % 保存数据和图表
% save_filename = sprintf('pps_smoothness_analysis_%s.mat', datestr(now, 'yyyymmdd_HHMMSS'));
% save(save_filename, 'cPpsPol', 'smoothness_metrics', 'constraint_active', 'zero_contrib', ...
%      'cS', 'test_age', 'reasons', 'suggestions');
% 
% % 保存图表
% fig_filename = sprintf('pps_smoothness_plots_%s.png', datestr(now, 'yyyymmdd_HHMMSS'));
% saveas(gcf, fig_filename);
% 
% fprintf('✓ 分析结果已保存:\n');
% fprintf('  数据文件: %s\n', save_filename);
% fprintf('  图表文件: %s\n', fig_filename);

fprintf('\n=== PPS缴费策略平滑性分析完成 ===\n');