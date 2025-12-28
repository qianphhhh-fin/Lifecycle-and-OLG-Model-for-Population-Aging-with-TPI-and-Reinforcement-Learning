% --- START OF FILE main_olg_v6.m ---

% OLG Model V6 (PPS Income Tax Deferral, VFI with k_pps state):
% Find max sustainable PAYG replacement rate (rho_prime_payg)
% where rho_prime_payg = b_payg / avg_worker_gross_wage.
% PPS contributions are deductible from income before income tax (represented by tau_l).
% Achieved by endogenously adjusting tau_l to balance GBC (with TR_gov = 0),
% PAYG tax rate (theta_payg) is endogenous but capped at cS.theta_payg_max.
% VFI state now includes k_pps.

clc; clear; close all;
fprintf('=== OLG Model V6 (PPS Tax Deferral, VFI w k_pps): Max Sustainable Rho_prime_payg ===\n');
fprintf('    (Rho_prime_payg = b_payg / avg_worker_gross_wage)\n');
fprintf('    (TR_gov=0, tau_l endogenous, theta_payg capped, VFI state: k, k_pps, eps)\n');

%% 1. 固定模型参数与人口设置
fprintf('\n--- 1. 初始化参数 ---\n');
cS = main_olg_v6_utils.ParameterValues_HuggettStyle();
paramS = struct();

fprintf('参数已加载。nk=%d, nkpps=%d, nw=%d。\n', cS.nk, cS.nkpps, cS.nw); % Added nkpps
fprintf('年度年龄范围: %d-%d。年龄组数: %d。\n', cS.age1_orig, cS.ageLast_orig, cS.aD_new);
fprintf('固定税率: tau_k=%.2f, tau_c=%.2f。G/Y=%.2f, B/Y=%.2f。\n', cS.tau_k, cS.tau_c, cS.gov_exp_frac_Y, cS.gov_debt_frac_Y);
fprintf('PAYG 税率上限 (theta_payg_max): %.3f\n', cS.theta_payg_max);
fprintf('所得税率tau_l范围: [%.3f, %.3f], 总劳动税上限: %.3f\n', cS.tau_l_min, cS.tau_l_max, cS.max_total_labor_tax);
fprintf('PPS 年度缴费上限 (绝对值): %.2f, 比例上限: %.2f\n', cS.pps_annual_contrib_limit, cS.pps_max_contrib_frac);


%% 2. 模拟人口动态至稳态
fprintf('\n--- 2. 模拟人口动态 ---\n');
popS = main_olg_v6_utils.initPopulation(cS);
popS = main_olg_v6_utils.populationDynamics(popS, cS);
[Z_ss, ~, bgp_reached, bgp_period] = main_olg_v6_utils.detectSteadyStatePopulation(popS, cS);
paramS.Z_ss_counts = Z_ss;
Z_ss_total = sum(Z_ss);
Z_ss_norm_group = zeros(cS.aD_new,1); if Z_ss_total > 1e-9, Z_ss_norm_group = Z_ss / Z_ss_total; end
paramS.ageMassV = Z_ss_norm_group(:);
paramS.mass_workers_group = sum(paramS.ageMassV(1:cS.aR_new));
if paramS.mass_workers_group < 1e-9, error('Mass of workers is zero.'); end
Z_ss_norm_annual = zeros(cS.aD_orig,1);
if Z_ss_total > 1e-9
    for a_new_map_idx = 1:cS.aD_new, annual_indices = cS.physAgeMap{a_new_map_idx}; group_mass = Z_ss_norm_group(a_new_map_idx);
        num_years_in_group = length(annual_indices); if num_years_in_group > 0, mass_per_year = group_mass / num_years_in_group; Z_ss_norm_annual(annual_indices) = mass_per_year; end; end
    if sum(Z_ss_norm_annual) > 1e-9 && abs(sum(Z_ss_norm_annual) - 1.0) > 1e-6, Z_ss_norm_annual = Z_ss_norm_annual / sum(Z_ss_norm_annual);
    elseif sum(Z_ss_norm_annual) < 1e-9, Z_ss_norm_annual(:) = 1/cS.aD_orig; end
else, Z_ss_norm_annual(:) = 1/cS.aD_orig; end
paramS.Z_ss_norm_annual = Z_ss_norm_annual;
if Z_ss_total > 1e-9 && length(popS.totalPop) > 1 && popS.totalPop(end-1) > 1e-9
    pop_growth_factor_per_group_period = popS.totalPop(end)/popS.totalPop(end-1);
    paramS.popGrowthForDebt = pop_growth_factor_per_group_period^(1/cS.yearStep) - 1;
else, paramS.popGrowthForDebt = cS.popGrowth_orig; end
fprintf('人口参数计算完毕。年化稳态人口增长率: %.4f\n', paramS.popGrowthForDebt);
fprintf('稳态工人占比 (基于年龄组): %.4f\n', paramS.mass_workers_group);

%% 3. 预计算劳动供给和禀赋过程
fprintf('\n--- 3. 预计算劳动 ---\n');
[paramS.leLogGridV, paramS.leTrProbM, paramS.leProb1V] = main_olg_v6_utils.EarningProcess_olgm(cS);
paramS.leGridV = exp(paramS.leLogGridV(:)); paramS.ageEffV_new = cS.ageEffV_new;
eIdxM = main_olg_v6_utils.LaborEndowSimulation_olgm(cS, paramS);
[~, L_per_capita] = main_olg_v6_utils.LaborSupply_Huggett(eIdxM, cS, paramS, paramS.ageMassV);
fprintf('总劳动供给 (L, 效率单位, 总体人均): %.4f\n', L_per_capita);
if L_per_capita <= 0, L_per_capita = 1e-6; warning('L_per_capita was zero.'); end
paramS.L_per_capita = L_per_capita;

%% 4. 寻找最大可持续替代率 rho_prime_payg_max
fprintf('\n--- 4. 寻找最大可持续替代率 rho_prime_payg_max (VFI w k_pps state) ---\n');

rho_payg_min = 0.01;
rho_payg_max_search_upper_bound = 0.7; 
max_iter_rho_search = 25; 
tol_rho_search = 1e-3;    

K_global_guess = 15.0; 

rho_low = rho_payg_min;
rho_high = rho_payg_max_search_upper_bound;
rho_prime_payg_optimal = rho_payg_min; 
K_optimal = K_global_guess; 
tau_l_optimal = cS.tau_l_init_guess; 
theta_payg_optimal_calc = 0; 
final_eq_solution_details = struct(); 
final_eq_solution_found = false;

fprintf('开始搜索最大可持续替代率 rho_prime_payg ...\n');
fprintf('IterRho | Rho_low  | Rho_high | Rho_try  | K_tot_sol| K_pps_sol| Tau_l_sol| Theta_g_req| GBC_Res  | Feasible | Time\n');
fprintf('-----------------------------------------------------------------------------------------------------------------------------------\n');
paramS_for_inner_loop = paramS; 

for iter_rho = 1:max_iter_rho_search
    rho_prime_try = (rho_low + rho_high) / 2;
    iter_rho_start_time = tic;

    if iter_rho > 1 
        paramS_for_inner_loop.suppress_inner_print_header = true;
        paramS_for_inner_loop.suppress_initial_theta_print = true;
    else
        paramS_for_inner_loop.suppress_inner_print_header = false;
        paramS_for_inner_loop.suppress_initial_theta_print = false;
    end

    [K_solution, tau_l_solution, gbc_residual, eq_found_for_rho_try, solution_details_iter] = ...
        main_olg_v6_utils.solve_K_tau_l_for_rho_prime(rho_prime_try, K_global_guess, cS, paramS_for_inner_loop, eIdxM);

    iter_rho_time = toc(iter_rho_start_time);
    
    current_theta_payg_calc = NaN; K_pps_model_from_inner = NaN;
    if isfield(solution_details_iter, 'theta_payg_required_before_cap'), current_theta_payg_calc = solution_details_iter.theta_payg_required_before_cap; end
    if isfield(solution_details_iter, 'K_model_pps'), K_pps_model_from_inner = solution_details_iter.K_model_pps; end

    fprintf('%7d | %8.4f | %8.4f | %8.4f | %8.4f | %8.4f | %8.4f | %11.4f | %8.2e | %8d | %.2fs\n', ...
        iter_rho, rho_low, rho_high, rho_prime_try, ...
        K_solution, K_pps_model_from_inner, tau_l_solution, ... 
        current_theta_payg_calc, gbc_residual, eq_found_for_rho_try, iter_rho_time);

    if eq_found_for_rho_try
        rho_prime_payg_optimal = rho_prime_try; K_optimal = K_solution; tau_l_optimal = tau_l_solution;
        theta_payg_optimal_calc = current_theta_payg_calc; final_eq_solution_details = solution_details_iter; 
        rho_low = rho_prime_try; final_eq_solution_found = true; K_global_guess = K_solution; 
    else, rho_high = rho_prime_try; end

    if (rho_high - rho_low) < tol_rho_search, fprintf('替代率 (rho_prime_payg) 搜索收敛。\n'); break; end
end

if ~final_eq_solution_found
    warning('未能找到满足所有条件的最大可持续替代率 rho_prime_payg。将使用最后记录的可行值或最小值。');
    if isempty(fields(final_eq_solution_details)) && rho_prime_payg_optimal > rho_payg_min + tol_rho_search / 2 
         fprintf('最后尝试的可行 rho_prime_payg_optimal = %.4f，重新计算其均衡...\n', rho_prime_payg_optimal);
        paramS_for_inner_loop.suppress_inner_print_header = false; paramS_for_inner_loop.suppress_initial_theta_print = false;
        [K_eq_fallback, tau_l_eq_fallback, ~, eq_found_fallback, final_eq_solution_details_fallback] = ... 
             main_olg_v6_utils.solve_K_tau_l_for_rho_prime(rho_prime_payg_optimal, K_global_guess, cS, paramS_for_inner_loop, eIdxM);
        if eq_found_fallback && ~isempty(fields(final_eq_solution_details_fallback)) && isfield(final_eq_solution_details_fallback,'MPL_gross') && ~isnan(final_eq_solution_details_fallback.MPL_gross)
            K_optimal = K_eq_fallback; tau_l_optimal = tau_l_eq_fallback; final_eq_solution_details = final_eq_solution_details_fallback;
            if isfield(final_eq_solution_details, 'theta_payg_required_before_cap'), theta_payg_optimal_calc = final_eq_solution_details.theta_payg_required_before_cap; else, theta_payg_optimal_calc = NaN; end
        else, error('重新计算最后的rho_prime_payg_optimal失败或未返回有效结果。'); end
    elseif isempty(fields(final_eq_solution_details)), error('在指定的PAYG税率上限和财政约束下，即使最低替代率也无法持续。'); end
    rho_prime_payg_eq = rho_prime_payg_optimal; K_eq = K_optimal; tau_l_eq = tau_l_optimal; 
else
    rho_prime_payg_eq = rho_prime_payg_optimal; K_eq = K_optimal; tau_l_eq = tau_l_optimal;
    fprintf('找到的最大可持续替代率 rho_prime_payg_eq = %.4f (理论theta_payg_req=%.4f)\n', rho_prime_payg_eq, theta_payg_optimal_calc);
end


%% 5. 分析和绘制最终均衡结果
fprintf('\n--- 5. 最终均衡结果与绘图 (rho_prime_payg_eq=%.4f, tau_l_eq=%.4f, TR_gov=0) ---\n', rho_prime_payg_eq, tau_l_eq);

if isempty(fields(final_eq_solution_details)) || isnan(K_eq) || ~isfield(final_eq_solution_details, 'MPL_gross') || isnan(final_eq_solution_details.MPL_gross)
    error('最终均衡的详细信息未能获取或无效（例如MPL_gross缺失），无法进行分析。');
end

paramS_eq = paramS; 
paramS_eq.tau_l = tau_l_eq; 
if isfield(final_eq_solution_details, 'theta_payg')
    paramS_eq.theta_payg_actual_for_hh = final_eq_solution_details.theta_payg;
else 
    temp_MPL_gross_for_theta_calc = main_olg_v6_utils.HHPrices_Huggett(K_eq, paramS.L_per_capita, cS); % Use K_eq
    temp_avg_worker_wage = (temp_MPL_gross_for_theta_calc * paramS.L_per_capita) / paramS.mass_workers_group;
    temp_b_payg = rho_prime_payg_eq * temp_avg_worker_wage;
    temp_theta_req = rho_prime_payg_eq * (sum(paramS.ageMassV(cS.aR_new+1:cS.aD_new)) / paramS.mass_workers_group);
    temp_theta_act = min(temp_theta_req, cS.theta_payg_max);
    if (temp_theta_act + tau_l_eq) > cS.max_total_labor_tax
        temp_theta_act = max(0, cS.max_total_labor_tax - tau_l_eq);
    end
    paramS_eq.theta_payg_actual_for_hh = max(0, temp_theta_act);
    warning('final_eq_solution_details.theta_payg (actual) not found, recomputed for final VFI.');
end

[R_mkt_gross_factor_eq_final, MPL_gross_eq_final] = main_olg_v6_utils.HHPrices_Huggett(K_eq, paramS.L_per_capita, cS);
r_mkt_gross_eq_final = R_mkt_gross_factor_eq_final - 1;
r_k_net_hh_eq_final = r_mkt_gross_eq_final * (1 - cS.tau_k);
R_k_net_factor_hh_eq_final = 1 + r_k_net_hh_eq_final;

avg_worker_gross_wage_eq_final = (MPL_gross_eq_final * paramS.L_per_capita) / paramS.mass_workers_group;
b_payg_eq_final = rho_prime_payg_eq * avg_worker_gross_wage_eq_final;
bV_eq_new_final = zeros(1, cS.aD_new);
if cS.aR_new < cS.aD_new, bV_eq_new_final(cS.aR_new + 1 : cS.aD_new) = b_payg_eq_final; end

T_bequest_eq_final = 0; 
if isfield(final_eq_solution_details, 'T_bequest_Model') && ~isnan(final_eq_solution_details.T_bequest_Model)
    T_bequest_eq_final = final_eq_solution_details.T_bequest_Model;
else, warning('T_bequest_Model not found in final_eq_solution_details for final VFI. Using 0.'); end
TR_total_eq_final_for_vfi = T_bequest_eq_final;

fprintf('Final VFI call with: MPL_gross=%.4f, tau_l=%.4f, theta_payg_actual=%.4f, TR_total=%.4f\n', ...
    MPL_gross_eq_final, paramS_eq.tau_l, paramS_eq.theta_payg_actual_for_hh, TR_total_eq_final_for_vfi);

% --- 调用最终的VFI ---
[cPolM_eq, kPolM_eq, cPpsPolM_eq, valueM_eq] = main_olg_v6_utils.HHSolution_VFI_Huggett(...
    R_k_net_factor_hh_eq_final, MPL_gross_eq_final, TR_total_eq_final_for_vfi, bV_eq_new_final, paramS_eq, cS);
% cPolM_eq, kPolM_eq, cPpsPolM_eq 是 4D 矩阵 (nk, nkpps, nw, naD_new)

fprintf('Simulating final equilibrium distribution...\n');
[kHistM_eq, kPpsHistM_eq, cHistM_eq] = main_olg_v6_utils.HHSimulation_olgm(...
    kPolM_eq, cPpsPolM_eq, cPolM_eq, eIdxM, ...
    R_k_net_factor_hh_eq_final, MPL_gross_eq_final, TR_total_eq_final_for_vfi, bV_eq_new_final, paramS_eq, cS);

% ... (之前的宏观量计算和打印保持不变) ...
K_nonpps_eq_agg = mean(kHistM_eq, 1) * paramS.Z_ss_norm_annual;
K_pps_eq_agg = 0; 
if cS.pps_active && cS.pps_in_K && cS.pps_max_contrib_frac > 0 && ~isempty(kPpsHistM_eq) 
    K_pps_eq_agg = mean(kPpsHistM_eq, 1) * paramS.Z_ss_norm_annual; 
end
Actual_K_eq_final = K_nonpps_eq_agg + K_pps_eq_agg;
C_eq_final = mean(cHistM_eq,1) * paramS.Z_ss_norm_annual;
Y_eq_final = cS.A * (Actual_K_eq_final^cS.alpha) * (paramS.L_per_capita^(1-cS.alpha));
G_eq_final = cS.gov_exp_frac_Y * Y_eq_final;
B_eq_final = cS.gov_debt_frac_Y * Y_eq_final;

fprintf('K_eq from outer loop: %.4f, K_eq from final sim: %.4f\n', K_eq, Actual_K_eq_final);
if abs(K_eq - Actual_K_eq_final) > 2e-2 && K_eq > 1e-9 
    warning('K_eq from rho search and K from final simulation differ significantly by %.3e. Y, G, B will use final sim K.', abs(K_eq - Actual_K_eq_final));
    Y_eq_final = cS.A * (Actual_K_eq_final^cS.alpha) * (paramS.L_per_capita^(1-cS.alpha));
    G_eq_final = cS.gov_exp_frac_Y * Y_eq_final; B_eq_final = cS.gov_debt_frac_Y * Y_eq_final;
end

fprintf('均衡总生产性资本 (K*): %.4f (总体人均)\n', Actual_K_eq_final);
fprintf('  其中: 非PPS资本 K: %.4f, PPS资本 K: %.4f\n', K_nonpps_eq_agg, K_pps_eq_agg);
fprintf('均衡总劳动 (L): %.4f\n', paramS.L_per_capita);
fprintf('均衡总产出 (Y*): %.4f\n', Y_eq_final);
fprintf('均衡市场利率因子 (R_mkt*): %.4f (r_mkt*=%.4f)\n', R_mkt_gross_factor_eq_final, r_mkt_gross_eq_final);
fprintf('  家庭净回报率因子 (R_k_net_hh*): %.4f (r_k_net_hh*=%.4f)\n', R_k_net_factor_hh_eq_final, r_k_net_hh_eq_final);
fprintf('均衡总工资率 (MPL_gross*): %.4f\n', MPL_gross_eq_final);
fprintf('目标PAYG替代率 (rho_prime_payg_eq*): %.4f (b_payg / avg_worker_gross_wage)\n', rho_prime_payg_eq);
fprintf('均衡内生PAYG税率 (theta_payg_eq*, capped at %.3f): %.4f\n', cS.theta_payg_max, paramS_eq.theta_payg_actual_for_hh);
if isfield(final_eq_solution_details, 'theta_payg_required_before_cap')
    fprintf('  (理论所需 PAYG 税率，未设上限前: %.4f)\n', final_eq_solution_details.theta_payg_required_before_cap);
end
fprintf('均衡内生"所得"税率 (tau_l_eq*): %.4f\n', tau_l_eq);
fprintf('  固定资本所得税率 (tau_k): %.4f, 固定消费税率 (tau_c): %.4f\n', cS.tau_k, cS.tau_c);
w_net_hh_display = MPL_gross_eq_final * (1 - paramS_eq.theta_payg_actual_for_hh - tau_l_eq);
fprintf('均衡家庭净工资率 (w_net_hh*): %.4f (基于最终税率)\n', w_net_hh_display);
fprintf('均衡PAYG福利 (b_payg*): %.4f (每位退休人员)\n', b_payg_eq_final);
fprintf('均衡总净转移支付 (TR_total*): %.4f (总体人均, TR_gov=0)\n', TR_total_eq_final_for_vfi);
fprintf('  其中意外遗赠 (T_bequest*): %.4f\n', T_bequest_eq_final);
fprintf('  其中政府净转移 (TR_gov*): 0.0000 (按设定)\n');
fprintf('均衡政府消费 (G*): %.4f (G/Y=%.2f)\n', G_eq_final, G_eq_final/Y_eq_final);
fprintf('均衡政府债务 (B*): %.4f (B/Y=%.2f)\n', B_eq_final, B_eq_final/Y_eq_final);
fprintf('均衡 K/Y 比率: %.4f\n', Actual_K_eq_final / Y_eq_final );
fprintf('均衡 C/Y 比率: %.4f\n', C_eq_final / Y_eq_final);

achieved_replacement_rate_final = 0;
if avg_worker_gross_wage_eq_final > 1e-9, achieved_replacement_rate_final = b_payg_eq_final / avg_worker_gross_wage_eq_final; end
fprintf('实际达成替代率 (b_payg / avg_worker_gross_wage): %.4f (应接近 rho_prime_payg_eq)\n', achieved_replacement_rate_final);
if abs(achieved_replacement_rate_final - rho_prime_payg_eq) > 1e-2 && rho_prime_payg_eq > 1e-9, warning('最终达成的替代率与目标替代率差异较大。'); end

final_gbc_residual = main_olg_v6_utils.check_gbc_residual(Actual_K_eq_final, C_eq_final, Y_eq_final, G_eq_final, B_eq_final, MPL_gross_eq_final, r_mkt_gross_eq_final, paramS_eq.theta_payg_actual_for_hh, tau_l_eq, b_payg_eq_final, T_bequest_eq_final, 0, cS, paramS_eq);
fprintf('最终GBC Check at Equilibrium: Residual = %.4e\n', final_gbc_residual);

% --- 新增：绘制策略函数 ---
fprintf('\n绘制最终均衡的策略函数...\n');

% 选择要绘图的参数
plot_a_idx = min(round(cS.aR_new / 2), cS.aD_new); % 例如，工作期中间的一个年龄组
if plot_a_idx == 0, plot_a_idx = 1; end
plot_ie_idx = round(cS.nw / 2);                 % 例如，中间的劳动效率状态
% 选择几个k_pps的网格点进行绘图，或者绘制一个二维图
plot_nkpps_to_show = min(3, cS.nkpps); % 最多显示3条k_pps的线
plot_ikpps_indices = round(linspace(1, cS.nkpps, plot_nkpps_to_show));
if cS.nkpps == 0, plot_ikpps_indices = []; end % Handle case with no kpps grid

figure_title_suffix_base = sprintf('年龄组 %d (真实年龄约 %d), 效率状态 %d', ...
    plot_a_idx, cS.physAgeV_new(plot_a_idx), plot_ie_idx);

if cS.nk > 1 && ~isempty(plot_ikpps_indices)
    % 绘制 k_prime(k) for different k_pps levels
    figure('Name', ['非PPS储蓄策略 k''(k | k_pps): ' figure_title_suffix_base]);
    hold on;
    colors = lines(plot_nkpps_to_show);
    for i_plot = 1:plot_nkpps_to_show
        ikpps = plot_ikpps_indices(i_plot);
        k_prime_slice = squeeze(kPolM_eq(:, ikpps, plot_ie_idx, plot_a_idx));
        plot(cS.kGridV, k_prime_slice, 'LineWidth', 1.5, 'DisplayName', sprintf('k_{pps}=%.2f', cS.kppsGridV(ikpps)), 'Color', colors(i_plot,:));
    end
    plot(cS.kGridV, cS.kGridV, 'k--', 'DisplayName', 'k''=k');
    hold off;
    xlabel('当前非PPS资产 k'); ylabel('下一期非PPS资产 k''');
    title({'非PPS储蓄策略 k''(k | k_{pps})'; figure_title_suffix_base});
    legend('show', 'Location', 'best'); grid on;

    % 绘制 c_pps(k) for different k_pps levels
    if cS.pps_active && cS.pps_max_contrib_frac > 0
        figure('Name', ['PPS缴费策略 c_{pps}(k | k_pps): ' figure_title_suffix_base]);
        hold on;
        for i_plot = 1:plot_nkpps_to_show
            ikpps = plot_ikpps_indices(i_plot);
            cpps_slice = squeeze(cPpsPolM_eq(:, ikpps, plot_ie_idx, plot_a_idx));
            plot(cS.kGridV, cpps_slice, 'LineWidth', 1.5, 'DisplayName', sprintf('k_{pps}=%.2f', cS.kppsGridV(ikpps)), 'Color', colors(i_plot,:));
        end
        hold off;
        xlabel('当前非PPS资产 k'); ylabel('PPS缴费 c_{pps}');
        title({'PPS缴费策略 c_{pps}(k | k_{pps})'; figure_title_suffix_base});
        legend('show', 'Location', 'best'); grid on;
    end

    % 绘制 c(k) for different k_pps levels
    figure('Name', ['消费策略 c(k | k_pps): ' figure_title_suffix_base]);
    hold on;
    for i_plot = 1:plot_nkpps_to_show
        ikpps = plot_ikpps_indices(i_plot);
        c_slice = squeeze(cPolM_eq(:, ikpps, plot_ie_idx, plot_a_idx));
        plot(cS.kGridV, c_slice, 'LineWidth', 1.5, 'DisplayName', sprintf('k_{pps}=%.2f', cS.kppsGridV(ikpps)), 'Color', colors(i_plot,:));
    end
    hold off;
    xlabel('当前非PPS资产 k'); ylabel('消费 c');
    title({'消费策略 c(k | k_{pps})'; figure_title_suffix_base});
    legend('show', 'Location', 'best'); grid on;

elseif cS.nk > 1 && cS.nkpps == 1 % Only k varies, k_pps is effectively scalar (or not a state)
    figure('Name', ['策略函数: ' figure_title_suffix_base ' (k_{pps} fixed/scalar)']);
    subplot(1,3,1);
    plot(cS.kGridV, squeeze(kPolM_eq(:, 1, plot_ie_idx, plot_a_idx)), 'b-o'); title('k''(k)'); xlabel('k'); grid on;
    hold on; plot(cS.kGridV, cS.kGridV, 'k--'); hold off;
    subplot(1,3,2);
    plot(cS.kGridV, squeeze(cPpsPolM_eq(:, 1, plot_ie_idx, plot_a_idx)), 'r-o'); title('c_{pps}(k)'); xlabel('k'); grid on;
    subplot(1,3,3);
    plot(cS.kGridV, squeeze(cPolM_eq(:, 1, plot_ie_idx, plot_a_idx)), 'g-o'); title('c(k)'); xlabel('k'); grid on;
    sgtitle(['策略函数切片: ' figure_title_suffix_base]);
else
    fprintf('无法绘制策略函数：nk或nkpps维度不足或plot_ikpps_indices为空。\n');
end

fprintf('\n--- V6 (PPS Income Tax Deferral with k_pps state in VFI) 分析完成 ---\n');