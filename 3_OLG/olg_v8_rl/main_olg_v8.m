% --- START OF FILE main_olg_v8.m ---

% OLG 模型 V8 (内生PPS缴费决策, PPS所得税递延, VFI w k_pps状态):
% 目标: 求解给定 PAYG 替代率 (rho_prime_payg_fixed) 下的均衡
% PPS缴费: 个体优化选择PPS缴费额，但受收入比例上限和年度绝对上限约束。
% 其他特性同Baseline:
%   - PPS 缴费可从所得税前扣除 (所得税率为 tau_l)。
%   - tau_l 内生调整以平衡政府一般预算 (TR_gov = 0)。
%   - PAYG 税率 (theta_payg) 内生决定，但有上限 cS.theta_payg_max。
%   - VFI 状态变量仍然包含 k_pps (PPS资产)。

clc; % 清除命令行窗口
clear; % 清除工作区变量
close all; % 关闭所有图形窗口
fprintf('=== OLG 模型 V8 (内生PPS缴费, 固定 Rho_prime_payg, VFI w k_pps) ===\n');
fprintf('    (Rho_prime_payg 固定, TR_gov=0, tau_l 内生, theta_payg 有上限)\n');
fprintf('    (VFI 状态: k, k_pps, eps; PPS缴费: 内生选择，有比例和绝对上限)\n');
% 关闭插值相关警告
warning('off', 'MATLAB:griddedInterpolant:MeshgridFailPointWarnId');
warning('off', 'MATLAB:griddedInterpolant:InterpEmptyGridId');
warning('off', 'MATLAB:griddedInterpolant:DegenerateGridId');


%% 1. 固定模型参数与人口设置
fprintf('\n--- 1. 初始化参数 ---\n');
cS = main_olg_v8_utils.ParameterValues_HuggettStyle(); % 调用v8的参数设置函数
paramS = struct(); % 初始化一个结构体用于存储派生参数

% *** Fixed PAYG Replacement Rate (as in Baseline) ***
cS.rho_prime_payg_fixed = 0.4; % 例如: 设定一个固定的PAYG替代率 (40%)
fprintf('>>> V8: 固定 PAYG 替代率 (rho_prime_payg_fixed): %.3f\n', cS.rho_prime_payg_fixed);
% ***************************************************************

fprintf('参数已加载。nk=%d, nkpps=%d, nw=%d, nPpsChoiceGrid=%d。\n', cS.nk, cS.nkpps, cS.nw, cS.n_pps_choice_grid_points); % 打印网格点数
fprintf('年度年龄范围: %d-%d。模型年龄组数: %d。\n', cS.age1_orig, cS.ageLast_orig, cS.aD_new); % 打印年龄信息
fprintf('固定税率: tau_k=%.2f, tau_c=%.2f。G/Y=%.2f, B/Y=%.2f。\n', cS.tau_k, cS.tau_c, cS.gov_exp_frac_Y, cS.gov_debt_frac_Y); % 打印固定税率和财政目标
fprintf('PAYG 税率上限 (theta_payg_max): %.3f\n', cS.theta_payg_max); % 打印PAYG税率上限
fprintf('所得税率 tau_l 范围: [%.3f, %.3f], 总劳动税上限: %.3f\n', cS.tau_l_min, cS.tau_l_max, cS.max_total_labor_tax); % 打印所得税率范围和总税负上限
fprintf('PPS 年度缴费上限 (绝对值): %.2f, 比例上限 (法定): %.2f\n', cS.pps_annual_contrib_limit, cS.pps_max_contrib_frac); % 打印PPS缴费上限
% fprintf('V7 PPS 缴费规则 (部分年龄组的比例): \n'); % V8中此规则不再直接用于决定缴费
% if cS.aD_new > 0
%     for i_disp_pps = 1:min(5, cS.aD_new)
%         fprintf('  年龄组 %d: %.3f\n', i_disp_pps, cS.pps_fixed_contrib_schedule_frac(i_disp_pps));
%     end
%     if cS.aD_new > 5, fprintf('  ...\n'); end
% end


%% 2. 模拟人口动态至稳态
fprintf('\n--- 2. 模拟人口动态 ---\n');
popS = main_olg_v8_utils.initPopulation(cS); % 初始化人口结构
popS = main_olg_v8_utils.populationDynamics(popS, cS); % 模拟人口动态演进
[Z_ss, ~, ~, ~] = main_olg_v8_utils.detectSteadyStatePopulation(popS, cS); % 检测稳态人口分布
paramS.Z_ss_counts = Z_ss; % 存储稳态人口计数 (按年龄组)
Z_ss_total = sum(Z_ss); % 计算稳态总人口
Z_ss_norm_group = zeros(cS.aD_new,1); % 初始化归一化的年龄组人口分布
if Z_ss_total > 1e-9
    Z_ss_norm_group = Z_ss / Z_ss_total; % 计算归一化分布
end
paramS.ageMassV = Z_ss_norm_group(:); % 存储为列向量 (年龄组人口占比)
paramS.mass_workers_group = sum(paramS.ageMassV(1:cS.aR_new)); % 计算工作年龄人口占比 (基于年龄组)
if paramS.mass_workers_group < 1e-9
    error('模型中稳态工作人口占比为零或过小。');
end

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
fprintf('人口参数计算完毕。年化稳态人口增长率 (用于债务): %.4f\n', paramS.popGrowthForDebt);
fprintf('稳态工人占比 (基于年龄组人口): %.4f\n', paramS.mass_workers_group);

%% 3. 预计算劳动供给和禀赋过程
fprintf('\n--- 3. 预计算劳动 ---\n');
[paramS.leLogGridV, paramS.leTrProbM, paramS.leProb1V] = main_olg_v8_utils.EarningProcess_olgm(cS);
paramS.leGridV = exp(paramS.leLogGridV(:));
paramS.ageEffV_new = cS.ageEffV_new;

eIdxM = main_olg_v8_utils.LaborEndowSimulation_olgm(cS, paramS);

[~, L_per_capita] = main_olg_v8_utils.LaborSupply_Huggett(eIdxM, cS, paramS, paramS.ageMassV);
fprintf('总劳动供给 (L, 效率单位, 总体人均): %.4f\n', L_per_capita);
if L_per_capita <= 0
    L_per_capita = 1e-6;
    warning('L_per_capita 为零或负，已重置为1e-6。');
end
paramS.L_per_capita = L_per_capita;

% 快速检查劳动效率网格
[leLogGridV, leTrProbM, leProb1V] = main_olg_v8_utils.EarningProcess_olgm(cS);
leGridV = exp(leLogGridV(:));

fprintf('劳动效率状态网格值：\n');
for i = 1:length(leGridV)
    fprintf('状态 %d: %.4f\n', i, leGridV(i));
end
fprintf('\n初始分布概率：\n');
for i = 1:length(leProb1V)
    fprintf('状态 %d: %.4f\n', i, leProb1V(i));
end

%% 4. 求解一般均衡 (给定固定的 rho_prime_payg_fixed)
fprintf('\n--- 4. 求解一般均衡 (固定 rho_prime_payg_fixed=%.3f) ---\n', cS.rho_prime_payg_fixed);
fprintf('  当前格点参数: n_k=%d, n_kpps=%d, n_w=%d, n_PpsChoiceGrid=%d\n', cS.nk, cS.nkpps, cS.nw, cS.n_pps_choice_grid_points);

K_global_guess = 2 ; % 总资本存量的初始猜测值
paramS_for_solver = paramS;
paramS_for_solver.suppress_inner_print_header = false;
paramS_for_solver.suppress_initial_theta_print = false;

fprintf('调用均衡求解器 solve_K_tau_l_for_rho_prime (V8) with fixed rho_prime_payg_fixed...\n');
solve_start_time = tic;
[K_eq, tau_l_eq, gbc_residual_eq, eq_found, final_eq_solution_details] = ...
    main_olg_v8_utils.solve_K_tau_l_for_rho_prime(cS.rho_prime_payg_fixed, K_global_guess, cS, paramS_for_solver, eIdxM);
solve_time = toc(solve_start_time);

fprintf('均衡求解完成。耗时: %.2f 秒。\n', solve_time);
fprintf('  均衡求解器返回状态: eq_found = %d\n', eq_found);
fprintf('  均衡结果: K_eq = %.4f, tau_l_eq = %.4f, GBC 残差 = %.3e\n', K_eq, tau_l_eq, gbc_residual_eq);

if ~eq_found || isnan(K_eq) || isnan(tau_l_eq)
    error('未能为固定的 rho_prime_payg_fixed = %.3f 找到均衡解。请检查模型参数或初始猜测值。', cS.rho_prime_payg_fixed);
end
if abs(gbc_residual_eq) > cS.gbc_tol_for_internal_loop * 10
    warning('最终均衡的GBC残差 (%.2e) 较大。可能均衡解的质量不高。', gbc_residual_eq);
end

theta_payg_optimal_calc = NaN;
if isfield(final_eq_solution_details, 'theta_payg_required_before_cap')
    theta_payg_optimal_calc = final_eq_solution_details.theta_payg_required_before_cap;
end
fprintf('  (理论所需PAYG税率，未考虑上限前: %.4f)\n', theta_payg_optimal_calc);

%% 5. 分析和绘制最终均衡结果
fprintf('\n--- 5. 最终均衡结果与绘图 (rho_prime_payg_fixed=%.4f, tau_l_eq=%.4f, TR_gov=0) ---\n', cS.rho_prime_payg_fixed, tau_l_eq);

if isempty(fields(final_eq_solution_details)) || isnan(K_eq) || ...
   ~isfield(final_eq_solution_details, 'MPL_gross') || isnan(final_eq_solution_details.MPL_gross)
    error('最终均衡的详细信息未能获取或无效，无法进行分析。');
end

paramS_eq = paramS;
paramS_eq.tau_l = tau_l_eq;

if isfield(final_eq_solution_details, 'theta_payg') && isfinite(final_eq_solution_details.theta_payg)
    paramS_eq.theta_payg_actual_for_hh = final_eq_solution_details.theta_payg;
else
    warning('final_eq_solution_details.theta_payg 未找到或无效，将基于均衡rho和约束重新计算。');
    [~, temp_MPL_gross_for_theta_calc] = main_olg_v8_utils.HHPrices_Huggett(K_eq, paramS.L_per_capita, cS);
    temp_mass_retirees = sum(paramS.ageMassV(cS.aR_new+1:cS.aD_new));
    temp_theta_req = cS.rho_prime_payg_fixed * (temp_mass_retirees / paramS.mass_workers_group);
    temp_theta_req = max(0, temp_theta_req);
    temp_theta_act = min(temp_theta_req, cS.theta_payg_max);
    if (temp_theta_act + tau_l_eq) > cS.max_total_labor_tax
        temp_theta_act = max(0, cS.max_total_labor_tax - tau_l_eq);
    end
    paramS_eq.theta_payg_actual_for_hh = max(0, temp_theta_act);
    fprintf('  重新计算的实际PAYG税率 (用于最终VFI): %.4f (理论需求: %.4f)\n', paramS_eq.theta_payg_actual_for_hh, temp_theta_req);
end
paramS_eq.pps_tax_deferral_active = cS.pps_active;

[R_mkt_gross_factor_eq_final, MPL_gross_eq_final] = main_olg_v8_utils.HHPrices_Huggett(K_eq, paramS.L_per_capita, cS);
r_mkt_gross_eq_final = R_mkt_gross_factor_eq_final - 1;
r_k_net_hh_eq_final = r_mkt_gross_eq_final * (1 - cS.tau_k);
R_k_net_factor_hh_eq_final = 1 + r_k_net_hh_eq_final;

avg_worker_gross_wage_eq_final = (MPL_gross_eq_final * paramS.L_per_capita) / paramS.mass_workers_group;
b_payg_eq_final = cS.rho_prime_payg_fixed * avg_worker_gross_wage_eq_final;
bV_eq_new_final = zeros(1, cS.aD_new);
if cS.aR_new < cS.aD_new
    bV_eq_new_final(cS.aR_new + 1 : cS.aD_new) = b_payg_eq_final;
end

T_bequest_eq_final = 0;
if isfield(final_eq_solution_details, 'T_bequest_Model') && isfinite(final_eq_solution_details.T_bequest_Model)
    T_bequest_eq_final = final_eq_solution_details.T_bequest_Model;
else
    warning('T_bequest_Model 未在 final_eq_solution_details 中找到或无效。将使用0。');
end
TR_total_eq_final_for_vfi = T_bequest_eq_final;

fprintf('最终 VFI 调用参数: MPL_gross=%.4f, tau_l=%.4f, theta_payg_actual=%.4f, TR_total=%.4f (T_bequest)\n', ...
    MPL_gross_eq_final, paramS_eq.tau_l, paramS_eq.theta_payg_actual_for_hh, TR_total_eq_final_for_vfi);

fprintf('调用最终的 HHSolution_VFI_Huggett (V8)...\n');
[cPolM_eq, kPolM_eq, cPpsPolM_choice_eq, valM_equilibrium] = main_olg_v8_utils.HHSolution_VFI_Huggett(...
    R_k_net_factor_hh_eq_final, MPL_gross_eq_final, TR_total_eq_final_for_vfi, bV_eq_new_final, paramS_eq, cS);

fprintf('模拟最终均衡的分布 (HHSimulation_olgm V8)...\n');
[kHistM_eq, kPpsHistM_eq, cHistM_eq] = main_olg_v8_utils.HHSimulation_olgm(...
    kPolM_eq, cPpsPolM_choice_eq, cPolM_eq, eIdxM, ... % Note: cPpsPolM_choice_eq used here
    R_k_net_factor_hh_eq_final, MPL_gross_eq_final, ...
    TR_total_eq_final_for_vfi, bV_eq_new_final, paramS_eq, cS);

K_nonpps_eq_agg = mean(kHistM_eq, 1) * paramS.Z_ss_norm_annual;
K_pps_eq_agg = 0;
if cS.pps_active && cS.pps_in_K && (cS.pps_max_contrib_frac > 0 || cS.pps_annual_contrib_limit > 0) && ~isempty(kPpsHistM_eq)
    K_pps_eq_agg = mean(kPpsHistM_eq, 1) * paramS.Z_ss_norm_annual;
end
Actual_K_eq_final = K_nonpps_eq_agg + K_pps_eq_agg;

C_eq_final = mean(cHistM_eq,1) * paramS.Z_ss_norm_annual;
Y_eq_final = cS.A * (Actual_K_eq_final^cS.alpha) * (paramS.L_per_capita^(1-cS.alpha));
G_eq_final = cS.gov_exp_frac_Y * Y_eq_final;
B_eq_final = cS.gov_debt_frac_Y * Y_eq_final;

fprintf('\n--- V8 最终均衡汇总 ---\n');
fprintf('K_eq (来自均衡求解器): %.4f, K_eq (来自最终模拟): %.4f\n', K_eq, Actual_K_eq_final);
if abs(K_eq - Actual_K_eq_final) > 2e-2 && K_eq > 1e-9
    warning('K_eq from solver and K from final simulation differ significantly by %.3e. Y, G, B将使用最终模拟的K值计算。', abs(K_eq - Actual_K_eq_final));
end

fprintf('均衡总生产性资本 (K*): %.4f (总体人均)\n', Actual_K_eq_final);
fprintf('  其中: 非PPS资本 K_non_pps: %.4f, PPS资本 K_pps: %.4f\n', K_nonpps_eq_agg, K_pps_eq_agg);
fprintf('均衡总劳动 (L, 效率单位, 总体人均): %.4f\n', paramS.L_per_capita);
fprintf('均衡总产出 (Y*): %.4f\n', Y_eq_final);
fprintf('均衡市场毛回报率因子 (R_mkt_gross*): %.4f (对应 r_mkt_gross*=%.4f)\n', R_mkt_gross_factor_eq_final, r_mkt_gross_eq_final);
fprintf('  家庭税后资本净回报率因子 (R_k_net_hh*): %.4f (对应 r_k_net_hh*=%.4f)\n', R_k_net_factor_hh_eq_final, r_k_net_hh_eq_final);
fprintf('均衡市场总工资率 (MPL_gross*): %.4f\n', MPL_gross_eq_final);
fprintf('目标PAYG替代率 (rho_prime_payg_fixed*): %.4f (b_payg / avg_worker_gross_wage)\n', cS.rho_prime_payg_fixed);
fprintf('均衡内生实际PAYG税率 (theta_payg_eq*, 上限 %.3f): %.4f\n', cS.theta_payg_max, paramS_eq.theta_payg_actual_for_hh);
if isfield(final_eq_solution_details, 'theta_payg_required_before_cap')
    fprintf('  (理论所需PAYG税率，未考虑上限前: %.4f)\n', final_eq_solution_details.theta_payg_required_before_cap);
end
fprintf('均衡内生"所得"税率 (tau_l_eq*): %.4f\n', tau_l_eq);
fprintf('  固定资本所得税率 (tau_k): %.2f, 固定消费税率 (tau_c): %.2f\n', cS.tau_k, cS.tau_c);

w_net_hh_approx_display = MPL_gross_eq_final * (1 - paramS_eq.theta_payg_actual_for_hh - paramS_eq.tau_l );
fprintf('均衡近似家庭平均净工资率 (w_gross * (1-theta_act-tau_l)): %.4f (仅为说明，未精确考虑PPS对tau_l基数的普遍影响)\n', w_net_hh_approx_display);

fprintf('均衡PAYG福利 (b_payg*, 每位退休者): %.4f\n', b_payg_eq_final);
fprintf('均衡总净转移支付 (TR_total*, 总体人均, TR_gov=0): %.4f\n', TR_total_eq_final_for_vfi);
fprintf('  其中意外遗赠 (T_bequest*): %.4f\n', T_bequest_eq_final);
fprintf('  其中政府直接转移 (TR_gov*): 0.0000 (按设定)\n');
fprintf('均衡政府消费 (G*): %.4f (G/Y* = %.3f)\n', G_eq_final, G_eq_final/Y_eq_final);
fprintf('均衡政府债务 (B*): %.4f (B/Y* = %.3f)\n', B_eq_final, B_eq_final/Y_eq_final);
fprintf('均衡 K*/Y* 比率: %.4f\n', Actual_K_eq_final / Y_eq_final );
fprintf('均衡 C*/Y* 比率: %.4f\n', C_eq_final / Y_eq_final);

achieved_replacement_rate_final = 0;
if avg_worker_gross_wage_eq_final > 1e-9
    achieved_replacement_rate_final = b_payg_eq_final / avg_worker_gross_wage_eq_final;
end
fprintf('实际达成替代率 (b_payg / avg_worker_gross_wage): %.4f (应接近 rho_prime_payg_fixed*)\n', achieved_replacement_rate_final);
if abs(achieved_replacement_rate_final - cS.rho_prime_payg_fixed) > 1e-3 && cS.rho_prime_payg_fixed > 1e-9
    warning('最终达成的替代率与目标替代率差异较大 (差异: %.3e)。请检查计算一致性。', abs(achieved_replacement_rate_final - cS.rho_prime_payg_fixed));
end

final_gbc_residual = main_olg_v8_utils.check_gbc_residual(Actual_K_eq_final, C_eq_final, Y_eq_final, ...
    G_eq_final, B_eq_final, MPL_gross_eq_final, r_mkt_gross_eq_final, ...
    paramS_eq.theta_payg_actual_for_hh, tau_l_eq, ...
    b_payg_eq_final, T_bequest_eq_final, 0, cS, paramS_eq);
fprintf('最终GBC(一般预算)检查 @ 均衡状态: Residual = %.4e\n', final_gbc_residual);
if abs(final_gbc_residual) > 1e-2
    warning('最终GBC残差 (%.3e) 较大。可能需要调整内层循环的 gbc_tol_for_internal_loop 或检查一致性。', final_gbc_residual);
end

fprintf('\n绘制最终均衡的策略函数...\n');
plot_a_idx = min(round(cS.aR_new / 2), cS.aD_new);
% 
if plot_a_idx == 0, plot_a_idx = 1; end
plot_ie_idx = round(cS.nw / 2);

plot_nkpps_to_show = min(3, cS.nkpps);
plot_ikpps_indices = [];
if cS.nkpps > 0
    plot_ikpps_indices = round(linspace(1, cS.nkpps, plot_nkpps_to_show));
else
    warning('nkpps为0，无法绘制k_pps相关的策略函数切片。');
end

figure_title_suffix_base = sprintf('年龄组 %d (真实年龄约 %d岁), 效率状态 %d', ...
    plot_a_idx, cS.physAgeV_new(plot_a_idx), plot_ie_idx);

if cS.nk > 1 && ~isempty(plot_ikpps_indices)
    figure('Name', ['V8: 非PPS储蓄策略 k''(k | k_pps): ' figure_title_suffix_base]);
    hold on;
    colors = lines(plot_nkpps_to_show);
    legend_entries_k_prime = cell(plot_nkpps_to_show + 1, 1);
    for i_plot = 1:plot_nkpps_to_show
        ikpps = plot_ikpps_indices(i_plot);
        k_prime_slice = squeeze(kPolM_eq(:, ikpps, plot_ie_idx, plot_a_idx));
        plot(cS.kGridV, k_prime_slice, 'LineWidth', 1.5, 'DisplayName', sprintf('k_{pps}=%.2f', cS.kppsGridV(ikpps)), 'Color', colors(i_plot,:));
        legend_entries_k_prime{i_plot} = sprintf('k_{pps}=%.2f', cS.kppsGridV(ikpps));
    end
    plot(cS.kGridV, cS.kGridV, 'k--', 'DisplayName', 'k''=k (45度线)');
    legend_entries_k_prime{plot_nkpps_to_show + 1} = 'k''=k';
    hold off;
    xlabel('当前非PPS资产 k');
    ylabel('下一期非PPS资产 k''');
    title({'V8: 非PPS储蓄策略 k''(k | k_{pps})'; figure_title_suffix_base});
    legend(legend_entries_k_prime, 'Location', 'best'); grid on;

    if cS.pps_active
        plot_a_idx_pps_contrib = min(5, cS.aR_new); % Choose a working age for PPS contribution plot
        if plot_a_idx_pps_contrib == 0, plot_a_idx_pps_contrib = 1; end
        figure_title_suffix_pps_contrib = sprintf('年龄组 %d (真实年龄约 %d岁), 效率状态 %d', ...
            plot_a_idx_pps_contrib, cS.physAgeV_new(plot_a_idx_pps_contrib), plot_ie_idx);

        figure('Name', ['V8: PPS缴费策略(选择) c_{pps}(k | k_pps): ' figure_title_suffix_pps_contrib]);
        hold on;
        legend_entries_cpps = cell(plot_nkpps_to_show, 1);
        for i_plot = 1:plot_nkpps_to_show
            ikpps = plot_ikpps_indices(i_plot);
            cpps_slice_choice = squeeze(cPpsPolM_choice_eq(:, ikpps, plot_ie_idx, plot_a_idx_pps_contrib)); % Use choice policy and correct age
            plot(cS.kGridV, cpps_slice_choice, 'LineWidth', 1.5, 'DisplayName', sprintf('k_{pps}=%.2f', cS.kppsGridV(ikpps)), 'Color', colors(i_plot,:));
            legend_entries_cpps{i_plot} = sprintf('k_{pps}=%.2f', cS.kppsGridV(ikpps));
        end
        hold off;
        xlabel('当前非PPS资产 k');
        ylabel('PPS缴费 c_{pps} (内生选择)');
        title({'V8: PPS缴费策略 (内生选择) c_{pps}(k | k_{pps})'; figure_title_suffix_pps_contrib});
        legend(legend_entries_cpps, 'Location', 'best'); grid on;
    end

    figure('Name', ['V8: 消费策略 c(k | k_pps): ' figure_title_suffix_base]);
    hold on;
    legend_entries_c = cell(plot_nkpps_to_show, 1);
    for i_plot = 1:plot_nkpps_to_show
        ikpps = plot_ikpps_indices(i_plot);
        c_slice = squeeze(cPolM_eq(:, ikpps, plot_ie_idx, plot_a_idx));
        plot(cS.kGridV, c_slice, 'LineWidth', 1.5, 'DisplayName', sprintf('k_{pps}=%.2f', cS.kppsGridV(ikpps)), 'Color', colors(i_plot,:));
        legend_entries_c{i_plot} = sprintf('k_{pps}=%.2f', cS.kppsGridV(ikpps));
    end
    hold off;
    xlabel('当前非PPS资产 k');
    ylabel('消费量 c');
    title({'V8: 消费策略 c(k | k_{pps})'; figure_title_suffix_base});
    legend(legend_entries_c, 'Location', 'best'); grid on;

elseif cS.nk > 1 && cS.nkpps == 1 % Only one k_pps point
    figure_title_suffix_k_only = sprintf('年龄组 %d (真实年龄约 %d岁), 效率状态 %d, k_{pps}=%.2f (固定)', ...
                                        plot_a_idx, cS.physAgeV_new(plot_a_idx), plot_ie_idx, cS.kppsGridV(1));
    figure('Name', ['V8: 策略函数 (k变化, k_pps固定): ' figure_title_suffix_k_only]);
    
    subplot(1,3,1);
    plot(cS.kGridV, squeeze(kPolM_eq(:, 1, plot_ie_idx, plot_a_idx)), 'b-o', 'LineWidth', 1.5);
    hold on; plot(cS.kGridV, cS.kGridV, 'k--'); hold off;
    title('非PPS储蓄 k''(k)'); xlabel('k'); grid on; legend('k''','k''=k', 'Location','best');

    subplot(1,3,2);
    plot_a_idx_pps_contrib_single_kpps = min(5, cS.aR_new);
    if plot_a_idx_pps_contrib_single_kpps == 0, plot_a_idx_pps_contrib_single_kpps = 1; end
    if cS.pps_active
        plot(cS.kGridV, squeeze(cPpsPolM_choice_eq(:, 1, plot_ie_idx, plot_a_idx_pps_contrib_single_kpps)), 'r-o', 'LineWidth', 1.5); % Use choice policy
        title('PPS缴费 c_{pps}(k) (选择)'); xlabel('k'); grid on;
    else
        plot(cS.kGridV, zeros(cS.nk,1), 'r-o', 'LineWidth', 1.5);
        title('PPS缴费 c_{pps}(k) (未激活)'); xlabel('k'); grid on;
    end

    subplot(1,3,3);
    plot(cS.kGridV, squeeze(cPolM_eq(:, 1, plot_ie_idx, plot_a_idx)), 'g-o', 'LineWidth', 1.5);
    title('消费 c(k)'); xlabel('k'); grid on;
    
    sgtitle(['V8: 策略函数切片: ' figure_title_suffix_k_only]);
else
    fprintf('无法绘制策略函数：nk或nkpps维度不足，或plot_ikpps_indices为空。\n');
    if cS.nk <= 1, fprintf('  原因: nk = %d (需要 > 1)\n', cS.nk); end
    if cS.nkpps == 0 && cS.pps_active, fprintf('  原因: nkpps = 0 (需要 > 0 才能绘制多条k_pps线, 或pps未激活)\n'); end
    if isempty(plot_ikpps_indices) && cS.nkpps > 1, fprintf('  原因: plot_ikpps_indices 为空 (内部逻辑错误)\n'); end
end

fprintf('\n--- V8 OLG 模型 (内生PPS缴费, 固定 Rho_prime_payg) 分析完成 ---\n');