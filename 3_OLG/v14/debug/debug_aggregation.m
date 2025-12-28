% =========================================================================
% == 调试脚本: debug_aggregation.m [v3 - 终极版]
% == 目的: 隔离并测试宏观总量（特别是劳动供给L）的聚合逻辑。
% =========================================================================
clear; close all;
addpath(pwd);
fprintf('=== [聚合逻辑专项调试脚本] ===\n\n');

%% 1. 初始化环境 (不变)
fprintf('--- 1. 初始化环境 ---\n');
cS = main_olg_v14_utils.ParameterValues_HuggettStyle();
paramS = struct();
ngrid = 50; cS.nk = ngrid; cS.nkpps = ngrid; cS.nkprime = ngrid; cS.npps = 5;
cS = main_olg_v14_utils.generateGrids(cS);
[Z_path, ~, ~] = main_olg_v14_utils.load_exogenous_paths(cS);
Z_ss_norm = Z_path(:,1);
cS = main_olg_v14_utils.calcaulte_theta_payg_path(cS, false);
[paramS.leLogGridV, paramS.leTrProbM, paramS.leProb1V] = main_olg_v14_utils.EarningProcess_olgm(cS);
paramS.leGridV = exp(paramS.leLogGridV(:));
fprintf('✅ 环境初始化完成。\n\n');


%% 2. 【核心修正】更精确地计算理论劳动供给 L_theory
fprintf('--- 2. 定义一个固定的宏观价格测试点 ---\n');
A_ss = cS.A;
theta_ss = cS.theta_path(1);

% =====================================================================
% == 【终极修正】通过模拟 cohort 生命周期来计算正确的 L_theory
% =====================================================================
% a. 初始化一个 cohort 在 20 岁时的禀赋分布
e_dist_by_age = zeros(cS.aD_new, cS.nw);
e_dist_by_age(1, :) = paramS.leProb1V'; % 年龄组1的分布是初始分布

% b. 模拟这个 cohort 的禀赋分布如何随年龄演进
for ia = 1:(cS.aD_new - 1)
    % e_dist_by_age(ia, :) 是 ia 年龄的分布
    % P' * dist_a(:) 得到 a+1 年龄的分布
    e_dist_by_age(ia+1, :) = e_dist_by_age(ia, :) * paramS.leTrProbM;
end

% c. 计算每个工作年龄组的平均禀赋冲击 E[epsilon | age]
mean_e_by_age = e_dist_by_age * paramS.leGridV(:);

% d. 计算总劳动供给 L_theory
% L = sum_{a=1 to a_R} ( Z_a * age_eff_a * E[epsilon | age=a] )
L_theory = 0;
for ia = 1:cS.aR_new
    L_theory = L_theory + Z_ss_norm(ia) * cS.ageEffV_new(ia) * mean_e_by_age(ia);
end
% =====================================================================

K_test = 1.2; % 选取一个合理的资本水平作为测试基准

% 基于这个【新计算】的 L_theory，生成一套完全一致的测试价格
M_test_prices = main_olg_v14_utils.get_prices_at_t(K_test, L_theory, A_ss, cS);
M_test_for_hh = M_test_prices;
mass_retirees_ss = sum(Z_ss_norm(cS.aR_new+1:end));
total_wage_bill_theory = M_test_prices.w_t * L_theory;

if mass_retirees_ss > 1e-9
    M_test_for_hh.b_t = (total_wage_bill_theory) / mass_retirees_ss * theta_ss;
else
    M_test_for_hh.b_t = 0;
end

fprintf('   测试价格点已生成:\n');
fprintf('   K = %.4f, L_theory (精确版) = %.4f, w = %.4f, r_mkt = %.4f, b = %.4f\n\n', ...
        K_test, L_theory, M_test_for_hh.w_t, M_test_for_hh.r_mkt_t, M_test_for_hh.b_t);


%% 3. 基于测试价格，求解微观策略和稳态分布 (不变)
fprintf('--- 3. 求解微观策略和稳态分布 ---\n');
cS_ss = cS; cS_ss.pps_active = false; cS_ss.theta_t = theta_ss;
if ~cS_ss.pps_active, cS_ss.nkpps = 1; cS_ss.npps = 1; cS_ss = main_olg_v14_utils.generateGrids(cS_ss); end
[~, kPolM, ~, ~] = main_olg_v14_utils.HHSolution_VFI_Huggett(M_test_for_hh, 0, paramS, cS_ss);
k_prime_idx = get_policy_index_matrix(kPolM, cS_ss);
Dist = solve_steady_state_distribution(k_prime_idx, paramS, cS_ss);
fprintf('✅ 微观数据生成完毕。\n\n');


%% 4. 【核心测试】从微观分布聚合宏观总量 (不变)
fprintf('--- 4. 【核心测试】从微观分布聚合宏观总量 ---\n');
[K_agg, C_agg, Tax_agg, Bequest_agg, L_agg] = ...
    aggregate_from_dist_nonstochastic(Dist, k_prime_idx, M_test_for_hh, Z_ss_norm, paramS, cS_ss);
fprintf('✅ 聚合完成。\n\n');


%% 5. 对比与诊断 (不变)
fprintf('--- 5. 对比与诊断 ---\n');
fprintf('%-25s | %12s | %12s | %12s\n', '变量', '理论值', '聚合值', '差异');
fprintf('%s\n', repmat('-', 65, 1));
fprintf('%-25s | %12.6f | %12.6f | %12.3e\n', '劳动供给 (L)', L_theory, L_agg, L_theory - L_agg);
TotalWageBill_agg = M_test_prices.w_t * L_agg;
fprintf('%-25s | %12.6f | %12.6f | %12.3e\n', '工资总额 (w*L)', total_wage_bill_theory, TotalWageBill_agg, total_wage_bill_theory - TotalWageBill_agg);
Y_theory = M_test_prices.Y_t;
I_agg = K_agg * cS.ddk;
G_agg = Tax_agg + Bequest_agg;
Y_agg_exp = C_agg + I_agg + G_agg;
fprintf('%-25s | %12.6f | %12.6f | %12.3e\n', '总产出 (Y)', Y_theory, Y_agg_exp, Y_theory - Y_agg_exp);
fprintf('%s\n', repmat('-', 65, 1));

if abs(L_theory - L_agg) < 1e-9
    fprintf('\n✅✅✅ 最终胜利！理论劳动供给与聚合劳动供给完全一致！✅✅✅\n');
    fprintf('   问题的根源在于，必须通过模拟生命周期来计算理论劳动供给。\n');
    fprintf('   现在，请将本脚本第2部分计算 L_theory 的代码块，\n');
    fprintf('   原封不动地替换你原始求解器中计算 L_ss 的部分。\n');
else
    fprintf('\n❌❌❌ 仍有差异！这表明 L_agg 的计算方式与我们对 L_theory 的新理解仍有不符。❌❌❌\n');
    fprintf('   请仔细检查 aggregate_from_dist_nonstochastic 中 L_agg 的聚合逻辑。\n');
end

%% === 辅助函数 (保持不变) ===
function k_prime_idx = get_policy_index_matrix(kPolM, cS), k_prime_idx = zeros(cS.nk, cS.nw, cS.aD_new, 'uint16'); for ia = 1:cS.aD_new, for ie = 1:cS.nw, for ik = 1:cS.nk, [~, idx] = min(abs(cS.kGridV - kPolM(ik, 1, ie, ia))); k_prime_idx(ik, ie, ia) = idx; end, end, end, end
function Dist = solve_steady_state_distribution(k_prime_idx, paramS, cS), Dist = zeros(cS.nk, cS.nw, cS.aD_new); Dist(1, :, 1) = paramS.leProb1V'; for ia = 1:(cS.aD_new - 1), Dist_next_age = zeros(cS.nk, cS.nw); for ik = 1:cS.nk, for ie = 1:cS.nw, cond_prob_at_state = Dist(ik, ie, ia); if cond_prob_at_state < 1e-20, continue; end, ik_prime = k_prime_idx(ik, ie, ia); transition_probs_e = paramS.leTrProbM(ie, :); Dist_next_age(ik_prime, :) = Dist_next_age(ik_prime, :) + cond_prob_at_state * transition_probs_e; end, end, if abs(sum(Dist_next_age, 'all') - 1.0) > 1e-6, warning('年龄组 %d 的分布总和 (%.4f) 不为1！', ia + 1, sum(Dist_next_age, 'all')); end, Dist(:, :, ia+1) = Dist_next_age; end, end
function [K_agg, C_agg, Tax_agg, Bequest_agg, L_agg] = aggregate_from_dist_nonstochastic(Dist, k_prime_idx, M_sim, Z_ss_norm, paramS, cS), K_agg = 0; C_agg = 0; Tax_agg = 0; Bequest_agg = 0; L_agg = 0; for ia = 1:cS.aD_new, age_pop_share = Z_ss_norm(ia); if age_pop_share < 1e-20, continue; end, for ie = 1:cS.nw, for ik = 1:cS.nk, cond_prob = Dist(ik, ie, ia); if cond_prob < 1e-20, continue; end, mass = cond_prob * age_pop_share; k_now = cS.kGridV(ik); epsilon_val = paramS.leGridV(ie); idx_k_prime = k_prime_idx(ik, ie, ia); k_prime = cS.kGridV(idx_k_prime); [c_val, tax_val, ~] = backout_c_tax_nonstochastic(k_now, k_prime, ia, epsilon_val, M_sim, cS, paramS); C_agg = C_agg + c_val * mass; Tax_agg = Tax_agg + tax_val * mass; K_agg = K_agg + k_prime * mass * cS.s_pathV(ia); Bequest_agg = Bequest_agg + k_prime * mass * (1 - cS.s_pathV(ia)); if ia <= cS.aR_new, individual_labor_supply = cS.ageEffV_new(ia) * epsilon_val; L_agg = L_agg + individual_labor_supply * mass; end, end, end, end, end
function [c_val, tax_val, labor_income_gross] = backout_c_tax_nonstochastic(k_now, k_prime, ia, epsilon_val, M_sim, cS, paramS), tr_per_capita = 0; cpps_decision = 0; cash_inflows_gross = k_now * (1 + M_sim.r_mkt_t) + tr_per_capita; labor_income_gross = 0; if ia <= cS.aR_new, labor_income_gross = M_sim.w_t * cS.ageEffV_new(ia) * epsilon_val; cash_inflows_gross = cash_inflows_gross + labor_income_gross; else, cash_inflows_gross = cash_inflows_gross + M_sim.b_t; end, payg_tax = cS.theta_t * labor_income_gross; labor_tax = cS.tau_l * max(0, labor_income_gross - cpps_decision - payg_tax); capital_tax = cS.tau_k * M_sim.r_mkt_t * k_now; general_taxes = labor_tax + capital_tax; cash_outflows_non_c = payg_tax + general_taxes + k_prime + cpps_decision; c_expend = cash_inflows_gross - cash_outflows_non_c; c_val = max(cS.cFloor, c_expend / (1 + cS.tau_c)); consumption_tax = c_val * cS.tau_c; tax_val = general_taxes + consumption_tax; end