% --- main_olg_v8.m (VFI最终对齐版) ---
%
% 目的:
%   - 作为OLG模型的纯VFI求解器。
%   - 在所有参数、宏观设定、均衡求解逻辑上与Python RL版本完全平行。
%   - 其计算出的均衡结果可直接与RL版本的均衡结果进行比较。
%
% 核心设定:
%   - 固定PAYG替代率 (rho_prime_payg_fixed)。
%   - 内生求解劳动所得税率 (tau_l) 以平衡政府一般预算。
%   - 使用离散网格搜索进行值函数迭代(VFI)。
%   - 完全基于年龄组进行模拟和计算。

clc; clear; close all; % 清理工作区
fprintf('=== OLG 模型 V8 (内生PPS缴费, 固定 Rho_prime_payg) [VFI最终对齐版] ===\n');
fprintf('    (与Python RL版本平行, 用于直接比较)\n');

% 关闭不必要的插值警告
warning('off', 'MATLAB:griddedInterpolant:MeshgridFailPointWarnId');
warning('off', 'MATLAB:griddedInterpolant:InterpEmptyGridId');

%% 1. 固定模型参数与人口设置
fprintf('\n--- 1. 初始化参数 ---\n');
cS = main_olg_v8_utils.ParameterValues_HuggettStyle(); % 调用对齐后的参数设置函数
paramS = struct(); % 初始化派生参数结构体

% 关键政策参数：固定PAYG替代率
cS.rho_prime_payg_fixed = 0.2;
fprintf('>>> V8: 固定 PAYG 替代率 (rho_prime_payg_fixed): %.3f\n', cS.rho_prime_payg_fixed);

%% 2. 模拟人口动态至稳态
fprintf('\n--- 2. 模拟人口动态 ---\n');
popS = main_olg_v8_utils.initPopulation(cS); % 初始化人口
popS = main_olg_v8_utils.populationDynamics(popS, cS); % 模拟动态演进
[Z_ss, ~, ~, ~] = main_olg_v8_utils.detectSteadyStatePopulation(popS, cS); % 获取稳态分布

paramS.Z_ss_counts = Z_ss;
Z_ss_total = sum(Z_ss);
if Z_ss_total < 1e-9, error('稳态总人口为零。'); end

paramS.ageMassV = Z_ss / Z_ss_total; % 年龄组人口占比
paramS.mass_workers_group = sum(paramS.ageMassV(1:cS.aR_new)); % 工作人口占比
if paramS.mass_workers_group < 1e-9, error('稳态工作人口占比过小。'); end

% 计算年化稳态人口增长率 (用于债务动态)
if length(popS.totalPop) > 1 && popS.totalPop(end-1) > 0
    paramS.popGrowthForDebt = (popS.totalPop(end) / popS.totalPop(end-1))^(1/cS.yearStep) - 1;
else
    paramS.popGrowthForDebt = 0.01; % 备用值
end
fprintf('人口参数计算完毕。年化稳态人口增长率: %.4f\n', paramS.popGrowthForDebt);

%% 3. 预计算劳动供给和禀赋过程
fprintf('\n--- 3. 预计算劳动 ---\n');
[paramS.leLogGridV, paramS.leTrProbM, paramS.leProb1V] = main_olg_v8_utils.EarningProcess_olgm(cS);
paramS.leGridV = exp(paramS.leLogGridV(:));
paramS.ageEffV_new = cS.ageEffV_new;

% 使用年龄组劳动禀赋模拟
eIdxM_group = main_olg_v8_utils.LaborEndowSimulation_olgm_AgeGroup(cS, paramS);

% 计算人均有效劳动供给
[~, L_per_capita] = main_olg_v8_utils.LaborSupply_Huggett(eIdxM_group, cS, paramS, paramS.ageMassV);
paramS.L_per_capita = max(L_per_capita, 1e-6); % 确保不为零

%% 4. 求解一般均衡
fprintf('\n--- 4. 求解一般均衡 (固定 rho_prime_payg_fixed=%.3f) ---\n', cS.rho_prime_payg_fixed);
K_global_guess = 2.0; % 总资本存量的初始猜测值

solve_start_time = tic;
[K_eq, tau_l_eq, gbc_residual_eq, eq_found, final_eq_solution_details] = ...
    main_olg_v8_utils.solve_K_tau_l_for_rho_prime(cS.rho_prime_payg_fixed, K_global_guess, cS, paramS, eIdxM_group);
solve_time = toc(solve_start_time);

fprintf('均衡求解完成。耗时: %.2f 秒。\n', solve_time);
fprintf('  均衡求解器返回状态: eq_found = %d\n', eq_found);
fprintf('  均衡结果: K_eq = %.4f, tau_l_eq = %.4f, GBC 残差 = %.3e\n', K_eq, tau_l_eq, gbc_residual_eq);
if ~eq_found, error('未能找到均衡解。'); end

%% 5. 分析和绘制最终均衡结果
fprintf('\n--- 5. 最终均衡结果与绘图 (rho_prime_payg_fixed=%.4f, tau_l_eq=%.4f) ---\n', cS.rho_prime_payg_fixed, tau_l_eq);

% 准备最终模拟所需的参数
paramS_eq = paramS;
paramS_eq.tau_l = tau_l_eq;
paramS_eq.theta_payg_actual_for_hh = final_eq_solution_details.theta_payg;
paramS_eq.pps_tax_deferral_active = cS.pps_active;

[R_mkt_gross_factor_eq, MPL_gross_eq] = main_olg_v8_utils.HHPrices_Huggett(K_eq, paramS.L_per_capita, cS);
R_k_net_factor_hh_eq = 1 + (R_mkt_gross_factor_eq - 1) * (1 - cS.tau_k);
bV_eq = zeros(1, cS.aD_new);
bV_eq(cS.aR_new+1:end) = final_eq_solution_details.b_payg;
TR_total_eq = 0.0; % 与RL版本对齐，无意外遗赠

% 使用均衡价格和税率，进行最终的VFI求解和模拟
[cPolM_eq, kPolM_eq, cPpsPolM_eq, ~] = main_olg_v8_utils.HHSolution_VFI_Huggett(R_k_net_factor_hh_eq, MPL_gross_eq, TR_total_eq, bV_eq, paramS_eq, cS);
[kHistM_eq, kPpsHistM_eq, cHistM_eq, ~] = main_olg_v8_utils.HHSimulation_olgm(kPolM_eq, cPpsPolM_eq, cPolM_eq, eIdxM_group, R_k_net_factor_hh_eq, MPL_gross_eq, TR_total_eq, bV_eq, paramS_eq, cS);

% 计算最终均衡的宏观总量
K_nonpps_agg = mean(kHistM_eq, 1) * paramS.ageMassV;
K_pps_agg = mean(kPpsHistM_eq, 1) * paramS.ageMassV;
Actual_K_eq = K_nonpps_agg + K_pps_agg;

% 打印最终均衡结果
fprintf('\n--- VFI 最终均衡汇总 ---\n');
fprintf('均衡总生产性资本 (K*): %.4f (非PPS: %.4f, PPS: %.4f)\n', Actual_K_eq, K_nonpps_agg, K_pps_agg);
fprintf('均衡总劳动 (L*): %.4f\n', paramS.L_per_capita);
fprintf('均衡总产出 (Y*): %.4f\n', main_olg_v8_utils.HHPrices_Huggett(Actual_K_eq, paramS.L_per_capita, cS));
fprintf('均衡工资率 (w*): %.4f, 净利率 (r*): %.4f\n', MPL_gross_eq, R_mkt_gross_factor_eq - 1 - cS.ddk);
fprintf('均衡劳动税率 (tau_l*): %.4f\n', tau_l_eq);
fprintf('均衡PAYG税率 (theta*): %.4f\n', final_eq_solution_details.theta_payg);
fprintf('目标替代率 (rho''): %.4f, 实际替代率: %.4f\n', cS.rho_prime_payg_fixed, final_eq_solution_details.b_payg / (MPL_gross_eq * paramS.L_per_capita / paramS.mass_workers_group));

%% 6. 绘图：均衡策略函数
fprintf('\n绘制最终均衡的策略函数...\n');
plot_a_idx = min(round(cS.aR_new / 2), cS.aD_new);
if plot_a_idx == 0, plot_a_idx = 1; end
plot_ie_idx = round(cS.nw / 2);

plot_nkpps_to_show = min(3, cS.nkpps);
if cS.nkpps > 0
    plot_ikpps_indices = round(linspace(1, cS.nkpps, plot_nkpps_to_show));
else
    plot_ikpps_indices = [];
end

figure_title_suffix = sprintf('年龄组 %d (约%d岁), 效率状态 %d', plot_a_idx, round(cS.age1_orig + (plot_a_idx-1)*cS.yearStep), plot_ie_idx);

if cS.nk > 1 && ~isempty(plot_ikpps_indices)
    % 非PPS储蓄策略
    figure('Name', ['VFI: 非PPS储蓄策略 k''(k | k_pps): ' figure_title_suffix]);
    hold on;
    colors = lines(plot_nkpps_to_show);
    for i = 1:plot_nkpps_to_show
        ikpps = plot_ikpps_indices(i);
        plot(cS.kGridV, squeeze(kPolM_eq(:, ikpps, plot_ie_idx, plot_a_idx)), 'LineWidth', 1.5, 'DisplayName', sprintf('k_{pps}=%.2f', cS.kppsGridV(ikpps)), 'Color', colors(i,:));
    end
    plot(cS.kGridV, cS.kGridV, 'k--', 'DisplayName', 'k''=k (45度线)');
    hold off;
    xlabel('当前非PPS资产 k'); ylabel('下一期非PPS资产 k''');
    title({'VFI: 非PPS储蓄策略 k''(k | k_{pps})'; figure_title_suffix});
    legend('Location', 'best'); grid on;

    % PPS缴费策略
    if cS.pps_active
        plot_a_idx_pps = min(5, cS.aR_new);
        figure_title_suffix_pps = sprintf('年龄组 %d (约%d岁), 效率状态 %d', plot_a_idx_pps, round(cS.age1_orig + (plot_a_idx_pps-1)*cS.yearStep), plot_ie_idx);
        figure('Name', ['VFI: PPS缴费策略 c_{pps}(k | k_pps): ' figure_title_suffix_pps]);
        hold on;
        for i = 1:plot_nkpps_to_show
            ikpps = plot_ikpps_indices(i);
            plot(cS.kGridV, squeeze(cPpsPolM_eq(:, ikpps, plot_ie_idx, plot_a_idx_pps)), 'LineWidth', 1.5, 'DisplayName', sprintf('k_{pps}=%.2f', cS.kppsGridV(ikpps)), 'Color', colors(i,:));
        end
        hold off;
        xlabel('当前非PPS资产 k'); ylabel('PPS缴费 c_{pps}');
        title({'VFI: PPS缴费策略 c_{pps}(k | k_{pps})'; figure_title_suffix_pps});
        legend('Location', 'best'); grid on;
    end
end

fprintf('\n--- VFI版本分析完成 ---\n');