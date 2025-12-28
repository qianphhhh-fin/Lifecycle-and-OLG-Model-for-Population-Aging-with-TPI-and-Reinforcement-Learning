% --- START OF FILE main_olg_v5.m ---

% OLG Model V5 (Revised with theta_payg_max):
% Find max sustainable PAYG replacement rate (rho_prime_payg)
% where rho_prime_payg = b_payg / avg_worker_gross_wage.
% Achieved by endogenously adjusting tau_l to balance GBC (with TR_gov = 0),
% AND PAYG tax rate (theta_payg) is endogenous but capped at cS.theta_payg_max.
% If required theta_payg > cS.theta_payg_max, the rho_prime_payg is considered unsustainable.

clc; clear; close all;
fprintf('=== OLG Model V5 (theta_payg_max): Max Sustainable Rho_prime_payg ===\n');
fprintf('    (Rho_prime_payg = b_payg / avg_worker_gross_wage)\n');
fprintf('    (TR_gov=0, tau_l endogenous, theta_payg endogenous with upper cap)\n');


%% 1. 固定模型参数与人口设置
fprintf('\n--- 1. 初始化参数 ---\n');
cS = main_olg_v5_utils.ParameterValues_HuggettStyle();
paramS = struct();

fprintf('参数已加载。nk=%d, nw=%d。\n', cS.nk, cS.nw);
fprintf('年度年龄范围: %d-%d。年龄组数: %d。\n', cS.age1_orig, cS.ageLast_orig, cS.aD_new);
fprintf('税率: tau_k=%.2f, tau_c=%.2f (固定)。G/Y=%.2f, B/Y=%.2f (固定)。\n', cS.tau_k, cS.tau_c, cS.gov_exp_frac_Y, cS.gov_debt_frac_Y);
fprintf('PAYG 税率上限 (theta_payg_max): %.3f\n', cS.theta_payg_max);

%% 2. 模拟人口动态至稳态
fprintf('\n--- 2. 模拟人口动态 ---\n');
popS = main_olg_v5_utils.initPopulation(cS);
popS = main_olg_v5_utils.populationDynamics(popS, cS);
[Z_ss, ~, bgp_reached, bgp_period] = main_olg_v5_utils.detectSteadyStatePopulation(popS, cS);
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
[paramS.leLogGridV, paramS.leTrProbM, paramS.leProb1V] = main_olg_v5_utils.EarningProcess_olgm(cS);
paramS.leGridV = exp(paramS.leLogGridV(:)); paramS.ageEffV_new = cS.ageEffV_new;
eIdxM = main_olg_v5_utils.LaborEndowSimulation_olgm(cS, paramS);
[~, L_per_capita] = main_olg_v5_utils.LaborSupply_Huggett(eIdxM, cS, paramS, paramS.ageMassV);
fprintf('总劳动供给 (L, 效率单位, 总体人均): %.4f\n', L_per_capita);
if L_per_capita <= 0, L_per_capita = 1e-6; warning('L_per_capita was zero.'); end
paramS.L_per_capita = L_per_capita;

%% 4. 寻找最大可持续替代率 rho_prime_payg_max
fprintf('\n--- 4. 寻找最大可持续替代率 rho_prime_payg_max (内生 tau_l, TR_gov=0, theta_payg capped) ---\n');

rho_payg_min = 0.01;
rho_payg_max_search_upper_bound = 0.7;
max_iter_rho_search = 25; % Increased iterations for outer search
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
fprintf('IterRho | Rho_low  | Rho_high | Rho_try  | K_sol    | Tau_l_sol| Theta_g_calc| GBC_Res  | Feasible | Time\n');
fprintf('-----------------------------------------------------------------------------------------------------------------\n');

for iter_rho = 1:max_iter_rho_search
    rho_prime_try = (rho_low + rho_high) / 2;
    iter_rho_start_time = tic;

    [K_solution, tau_l_solution, gbc_residual, eq_found_for_rho_try, solution_details_iter] = ...
        main_olg_v5_utils.solve_K_tau_l_for_rho_prime(rho_prime_try, K_global_guess, cS, paramS, eIdxM);

    iter_rho_time = toc(iter_rho_start_time);
    
    % Feasibility now also depends on whether theta_payg hit its cap AND was still insufficient
    % The 'eq_found_for_rho_try' from the inner function already considers if theta_payg_required > cS.theta_payg_max
    % and would return eq_found_for_rho_try = false in that case.
    
    % We retrieve the calculated theta_payg (before capping, if needed for display) or the capped one
    current_theta_payg_calc = NaN;
    if isfield(solution_details_iter, 'theta_payg_required_before_cap')
        current_theta_payg_calc = solution_details_iter.theta_payg_required_before_cap;
    elseif isfield(solution_details_iter, 'theta_payg') % if not capped or detail not saved
        current_theta_payg_calc = solution_details_iter.theta_payg;
    end


    fprintf('%7d | %8.4f | %8.4f | %8.4f | %8.4f | %8.4f | %11.4f | %8.2e | %8d | %.2fs\n', ...
        iter_rho, rho_low, rho_high, rho_prime_try, K_solution, tau_l_solution, current_theta_payg_calc, gbc_residual, eq_found_for_rho_try, iter_rho_time);

    if eq_found_for_rho_try % Inner loop converged AND theta_payg constraint met AND GBC is fine
        rho_prime_payg_optimal = rho_prime_try;
        K_optimal = K_solution;
        tau_l_optimal = tau_l_solution;
        theta_payg_optimal_calc = current_theta_payg_calc;
        final_eq_solution_details = solution_details_iter;
        rho_low = rho_prime_try;
        final_eq_solution_found = true;
        K_global_guess = K_solution;
    else
        rho_high = rho_prime_try;
    end

    if (rho_high - rho_low) < tol_rho_search
        fprintf('替代率 (rho_prime_payg) 搜索收敛。\n');
        break;
    end
end

if ~final_eq_solution_found
    warning('未能找到满足所有条件的最大可持续替代率 rho_prime_payg。将使用最后记录的可行值（如果有）或最小值。');
    % If rho_prime_payg_optimal is still at its initial rho_payg_min and no solution was ever found,
    % then final_eq_solution_details might be empty.
    if isempty(fields(final_eq_solution_details))
        if rho_prime_payg_optimal > rho_payg_min + tol_rho_search % Check if it moved from initial min
             fprintf('最后尝试的可行 rho_prime_payg_optimal = %.4f，重新计算其均衡...\n', rho_prime_payg_optimal);
            [K_eq, tau_l_eq, ~, ~, final_eq_solution_details] = ...
                 main_olg_v5_utils.solve_K_tau_l_for_rho_prime(rho_prime_payg_optimal, K_global_guess, cS, paramS, eIdxM);
             theta_payg_optimal_calc = final_eq_solution_details.theta_payg_required_before_cap; % Or .theta_payg
        else
            error('在指定的PAYG税率上限下，即使最低替代率也无法持续。请检查参数或theta_payg_max。');
        end
    end
    % Use the stored optimal values
    rho_prime_payg_eq = rho_prime_payg_optimal;
    K_eq = K_optimal;
    tau_l_eq = tau_l_optimal;
else
    rho_prime_payg_eq = rho_prime_payg_optimal;
    K_eq = K_optimal;
    tau_l_eq = tau_l_optimal;
    fprintf('找到的最大可持续替代率 rho_prime_payg_eq = %.4f (实际theta_payg_calc=%.4f)\n', rho_prime_payg_eq, theta_payg_optimal_calc);
end


%% 5. 分析和绘制最终均衡结果
fprintf('\n--- 5. 最终均衡结果与绘图 (rho_prime_payg_eq=%.4f, tau_l_eq=%.4f, TR_gov=0) ---\n', rho_prime_payg_eq, tau_l_eq);

if isempty(fields(final_eq_solution_details)) || isnan(K_eq)
    error('最终均衡的详细信息未能获取或无效，无法进行分析。');
end

paramS_eq = paramS;
paramS_eq.tau_l = tau_l_eq;

% These are from the successful run of solve_K_tau_l_for_rho_prime for rho_prime_payg_eq
MPL_gross_eq = final_eq_solution_details.MPL_gross;
theta_payg_eq = final_eq_solution_details.theta_payg; % This is the capped theta_payg
b_payg_eq = final_eq_solution_details.b_payg;
T_bequest_eq = final_eq_solution_details.T_bequest_Model;
TR_total_eq_final = T_bequest_eq;

r_mkt_gross_eq = final_eq_solution_details.R_mkt_gross -1;
r_k_net_hh_eq = r_mkt_gross_eq * (1 - cS.tau_k);
R_k_net_factor_hh_eq = 1 + r_k_net_hh_eq;
w_net_hh_eq = MPL_gross_eq * (1 - theta_payg_eq - tau_l_eq); % Use the actual theta_payg_eq applied
w_net_hh_eq = max(0, w_net_hh_eq);

bV_eq_new = zeros(1, cS.aD_new);
if cS.aR_new < cS.aD_new, bV_eq_new(cS.aR_new + 1 : cS.aD_new) = b_payg_eq; end

[cPolM_eq, kPolM_eq, cPpsPolM_eq, ~] = main_olg_v5_utils.HHSolution_VFI_Huggett(R_k_net_factor_hh_eq, w_net_hh_eq, TR_total_eq_final, bV_eq_new, paramS_eq, cS);
[kHistM_eq, kPpsHistM_eq, cHistM_eq] = main_olg_v5_utils.HHSimulation_olgm(kPolM_eq, cPpsPolM_eq, cPolM_eq, eIdxM, R_k_net_factor_hh_eq, w_net_hh_eq, TR_total_eq_final, bV_eq_new, paramS_eq, cS);

K_nonpps_eq = mean(kHistM_eq, 1) * paramS.Z_ss_norm_annual;
K_pps_eq = 0; if cS.pps_active && cS.pps_in_K && cS.pps_max_contrib_frac > 1e-9, K_pps_eq = mean(kPpsHistM_eq, 1) * paramS.Z_ss_norm_annual; end
Actual_K_eq_from_sim = K_nonpps_eq + K_pps_eq;
C_eq = mean(cHistM_eq,1) * paramS.Z_ss_norm_annual;
Y_eq = cS.A * (Actual_K_eq_from_sim^cS.alpha) * (paramS.L_per_capita^(1-cS.alpha));
G_eq = cS.gov_exp_frac_Y * Y_eq;
B_eq = cS.gov_debt_frac_Y * Y_eq;

fprintf('K_eq from outer loop: %.4f, K_eq from final sim: %.4f\n', K_eq, Actual_K_eq_from_sim);
if abs(K_eq - Actual_K_eq_from_sim) > 1e-2 && K_eq > 1e-9 % Relaxed tolerance for this check
    warning('K_eq from rho search and final sim differ by %.3e.', abs(K_eq - Actual_K_eq_from_sim));
end

fprintf('均衡总生产性资本 (K*): %.4f (总体人均)\n', Actual_K_eq_from_sim);
fprintf('  其中: 非PPS资本 K: %.4f, PPS资本 K: %.4f\n', K_nonpps_eq, K_pps_eq);
fprintf('均衡总劳动 (L): %.4f\n', paramS.L_per_capita);
fprintf('均衡总产出 (Y*): %.4f\n', Y_eq);
fprintf('均衡市场利率因子 (R_mkt*): %.4f (r_mkt*=%.4f)\n', final_eq_solution_details.R_mkt_gross, r_mkt_gross_eq);
fprintf('  家庭净回报率因子 (R_k_net_hh*): %.4f (r_k_net_hh*=%.4f)\n', R_k_net_factor_hh_eq, r_k_net_hh_eq-1);
fprintf('均衡总工资率 (MPL_gross*): %.4f\n', MPL_gross_eq);
fprintf('目标PAYG替代率 (rho_prime_payg_eq*): %.4f (b_payg / avg_worker_gross_wage)\n', rho_prime_payg_eq);
fprintf('均衡内生PAYG税率 (theta_payg_eq*, capped at %.3f): %.4f\n', cS.theta_payg_max, theta_payg_eq);
if isfield(final_eq_solution_details, 'theta_payg_required_before_cap')
    fprintf('  (理论所需 PAYG 税率，未设上限前: %.4f)\n', final_eq_solution_details.theta_payg_required_before_cap);
end
fprintf('均衡内生额外劳动税率 (tau_l_eq*): %.4f\n', tau_l_eq);
fprintf('  固定资本所得税率 (tau_k): %.4f, 固定消费税率 (tau_c): %.4f\n', cS.tau_k, cS.tau_c);
fprintf('均衡家庭净工资率 (w_net_hh*): %.4f\n', w_net_hh_eq);
fprintf('均衡PAYG福利 (b_payg*): %.4f (每位退休人员)\n', b_payg_eq);
fprintf('均衡总净转移支付 (TR_total*): %.4f (总体人均, TR_gov=0)\n', TR_total_eq_final);
fprintf('  其中意外遗赠 (T_bequest*): %.4f\n', T_bequest_eq);
fprintf('  其中政府净转移 (TR_gov*): 0.0000 (按设定)\n');
fprintf('均衡政府消费 (G*): %.4f (G/Y=%.2f)\n', G_eq, G_eq/Y_eq);
fprintf('均衡政府债务 (B*): %.4f (B/Y=%.2f)\n', B_eq, B_eq/Y_eq);
fprintf('均衡 K/Y 比率: %.4f\n', Actual_K_eq_from_sim / Y_eq );
fprintf('均衡 C/Y 比率: %.4f\n', C_eq / Y_eq);

avg_worker_gross_wage_eq_final = 0;
if paramS.mass_workers_group > 1e-9 && paramS.L_per_capita > 0 && MPL_gross_eq > 0
    avg_worker_gross_wage_eq_final = (MPL_gross_eq * paramS.L_per_capita) / paramS.mass_workers_group;
end
achieved_replacement_rate_final = 0;
if avg_worker_gross_wage_eq_final > 1e-9
    achieved_replacement_rate_final = b_payg_eq / avg_worker_gross_wage_eq_final;
end
fprintf('实际达成替代率 (b_payg / avg_worker_gross_wage): %.4f (应接近 rho_prime_payg_eq)\n', achieved_replacement_rate_final);
if abs(achieved_replacement_rate_final - rho_prime_payg_eq) > 1e-3 && rho_prime_payg_eq > 0
    warning('最终达成的替代率与目标替代率差异较大。这可能因为theta_payg被上限约束。');
end

final_gbc_residual = main_olg_v5_utils.check_gbc_residual(Actual_K_eq_from_sim, C_eq, Y_eq, G_eq, B_eq, ...
    MPL_gross_eq, r_mkt_gross_eq, theta_payg_eq, tau_l_eq, b_payg_eq, ...
    T_bequest_eq, 0, cS, paramS_eq);
fprintf('最终GBC Check at Equilibrium: Residual = %.4e\n', final_gbc_residual);

fprintf('\n--- V5 (theta_payg_max) 分析完成 ---\n');