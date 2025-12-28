% =========================================================================
% == 调试脚本: debug_asset_motion.m
% == 目的: 通过“总资产演进恒等式”侧面验证国民账户的内部一致性。
% == 方法:
% == 1. 求解模型稳态均衡。
% == 2. 分别计算恒等式 δ*K = Y - C - G 的两端：
% ==    - 左端 (LHS): 总投资 (由聚合资本 K_agg 决定)
% ==    - 右端 (RHS): 总储蓄 (由生产函数 Y_prod 和聚合的 C, G 决定)
% == 3. 检查两端的残差，验证存量-流量的会计闭环。
% =========================================================================
clear; close all;
addpath(pwd);
fprintf('=== [资产演进恒等式检验脚本] ===\n\n');

%% 1. 初始化环境 (与主脚本一致)
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


%% 2. 求解稳态均衡 (使用最终修正版的求解器)
fprintf('--- 2. 求解稳态均衡 ---\n');
[ss, eq_found, ~, ~] = solve_steady_state_nonstochastic(Z_ss_norm, cS, paramS); 
if ~eq_found, error('稳态求解失败'); end
fprintf('✅ 稳态求解完成。\n\n');


%% 3. 【核心检验】总资产演进恒等式：δ*K = Y - C - G
fprintf('--- 3. 【核心检验】总资产演进恒等式：δ*K = Y - C - G ---\n');

% --- 3.1 计算左端 (LHS): 总投资 I = δ*K ---
% 这个I完全由最终聚合出的资本存量决定。
LHS_Investment = ss.I; % 直接使用ss中已经算好的 I = ddk * K_physical

fprintf('   --- 左端 (LHS): 总投资 --- \n');
fprintf('   聚合资本存量 (K_agg) : %.6f\n', ss.K_physical);
fprintf('   折旧率 (δ)           : %.4f\n', cS.ddk);
fprintf('   => 总投资 (I = δ*K)    : %.6f\n\n', LHS_Investment);


% --- 3.2 计算右端 (RHS): 总储蓄 S = Y - C - G ---
% 这里的Y，我们刻意使用来自“生产函数”的那个，以进行交叉检验。
% C 和 G 则使用聚合值。
RHS_Saving = ss.Y_from_production - ss.C - ss.G;

fprintf('   --- 右端 (RHS): 总储蓄 --- \n');
fprintf('   生产函数产出 (Y_prod)  : %.6f\n', ss.Y_from_production);
fprintf('   聚合消费 (C_agg)       : %.6f\n', ss.C);
fprintf('   聚合政府购买 (G_agg)   : %.6f\n', ss.G);
fprintf('   => 总储蓄 (S = Y-C-G)   : %.6f\n\n', RHS_Saving);


% --- 3.3 对比与诊断 ---
residual = LHS_Investment - RHS_Saving;

fprintf('   --- 最终对比 --- \n');
fprintf('   LHS (总投资 I)       : %.6f\n', LHS_Investment);
fprintf('   RHS (总储蓄 S)       : %.6f\n', RHS_Saving);
fprintf('   检验残差 (I - S)     : %.3e\n\n', residual);


if abs(residual) < 1e-6
    fprintf('✅✅✅ 完美！资产演进恒等式精确成立 (I=S)！✅✅✅\n');
    fprintf('   这强有力地证明了模型的生产侧 (Y_prod) 和支出/储蓄侧 (C_agg, I_agg, G_agg) 是完全自洽的。\n');
else
    fprintf('❌❌❌ 失败！资产演进恒等式不成立！❌❌❌\n');
    fprintf('   残差 (I-S) 不为零，意味着生产出的价值在分配、消费、储蓄后，并没有精确地等于新投资。\n');
    fprintf('   这表明 Y_prod, C_agg, G_agg, I_agg 中至少有一项的计算与其他项不一致。\n');
    fprintf('   由于我们已验证 L, 这通常指向 K, C, G 的聚合或其依赖的税收计算存在问题。\n');
end


%% ========================================================================
%  == 所有需要的辅助函数 (均为最终修正版)
%  ========================================================================

function [ss, eq_found, Dist, k_prime_idx] = solve_steady_state_nonstochastic(Z_ss_norm, cS, paramS)
    K_physical_guess = 1.2;
    max_iter_K = 100; tol_K = 1e-8; damp_K = 0.4;
    eq_found = false; A_ss = cS.A; theta_ss = cS.theta_path(1);
    
    e_dist_by_age = zeros(cS.aD_new, cS.nw);
    e_dist_by_age(1, :) = paramS.leProb1V';
    for ia = 1:(cS.aD_new - 1)
        e_dist_by_age(ia+1, :) = e_dist_by_age(ia, :) * paramS.leTrProbM;
    end
    mean_e_by_age = e_dist_by_age * paramS.leGridV(:);
    L_ss = 0;
    for ia = 1:cS.aR_new
        L_ss = L_ss + Z_ss_norm(ia) * cS.ageEffV_new(ia) * mean_e_by_age(ia);
    end
    fprintf('   理论总劳动供给 L_ss: %.6f\n', L_ss);

    for iter_K = 1:max_iter_K
        M_ss_prices = main_olg_v14_utils.get_prices_at_t(K_physical_guess, L_ss, A_ss, cS);
        M_ss_for_hh = M_ss_prices;
        mass_retirees_ss = sum(Z_ss_norm(cS.aR_new+1:end));
        total_wage_bill = M_ss_prices.w_t * L_ss; 
        if mass_retirees_ss > 1e-9, M_ss_for_hh.b_t = (theta_ss * total_wage_bill) / mass_retirees_ss; else, M_ss_for_hh.b_t = 0; end
        
        cS_ss = cS; cS_ss.pps_active = false; cS_ss.theta_t = theta_ss;
        if ~cS_ss.pps_active, cS_ss.nkpps = 1; cS_ss.npps = 1; cS_ss = main_olg_v14_utils.generateGrids(cS_ss); end
        [~, kPolM, ~, ~] = main_olg_v14_utils.HHSolution_VFI_Huggett(M_ss_for_hh, 0, paramS, cS_ss);
        
        k_prime_idx = get_policy_index_matrix(kPolM, cS_ss);
        Dist = solve_steady_state_distribution(k_prime_idx, paramS, cS_ss);
        
        [K_pvt_model, ~, ~, ~, ~] = aggregate_from_dist_nonstochastic(Dist, k_prime_idx, M_ss_for_hh, Z_ss_norm, paramS, cS_ss);
        K_physical_model = K_pvt_model;
        
        K_error = K_physical_guess - K_physical_model;
        if abs(K_error) < tol_K, eq_found = true; break; end
        K_physical_guess = (1 - damp_K) * K_physical_guess + damp_K * K_physical_model;
    end

    if ~eq_found, warning('主循环(K)未收敛'); ss=struct(); return; end

    K_final = K_physical_model;
    M_final_prices = main_olg_v14_utils.get_prices_at_t(K_final, L_ss, A_ss, cS);
    M_final_for_hh = M_final_prices;
    mass_retirees_ss = sum(Z_ss_norm(cS.aR_new+1:end));
    total_wage_bill = M_final_prices.w_t * L_ss;
    if mass_retirees_ss > 1e-9, M_final_for_hh.b_t = (theta_ss * total_wage_bill) / mass_retirees_ss; else, M_final_for_hh.b_t = 0; end
    
    cS_ss = cS; cS_ss.pps_active = false; cS_ss.theta_t = theta_ss;
    if ~cS_ss.pps_active, cS_ss.nkpps = 1; cS_ss.npps = 1; cS_ss = main_olg_v14_utils.generateGrids(cS_ss); end
    [~, kPolM_final, ~, ~] = main_olg_v14_utils.HHSolution_VFI_Huggett(M_final_for_hh, 0, paramS, cS_ss);
    
    k_prime_idx_final = get_policy_index_matrix(kPolM_final, cS_ss);
    Dist_final = solve_steady_state_distribution(k_prime_idx_final, paramS, cS_ss);
    
    [K_pvt_final, C_final, Tax_final, Bequest_final, L_agg_final] = aggregate_from_dist_nonstochastic(Dist_final, k_prime_idx_final, M_final_for_hh, Z_ss_norm, paramS, cS_ss);
    
    I_final = cS.ddk * K_pvt_final;
    G_final = Tax_final + Bequest_final;
    Y_final_from_expenditure = C_final + I_final + G_final;

    ss = struct();
    ss.K_physical = K_pvt_final;
    ss.L = L_ss;
    ss.L_aggregated = L_agg_final;
    ss.Y = Y_final_from_expenditure;
    ss.w = M_final_prices.w_t;
    ss.r_mkt = M_final_prices.r_mkt_t;
    ss.C = C_final;
    ss.I = I_final;
    ss.G = G_final;
    ss.Y_from_production = M_final_prices.Y_t;
    
    Dist = Dist_final;
    k_prime_idx = k_prime_idx_final;
end

function k_prime_idx = get_policy_index_matrix(kPolM, cS), k_prime_idx = zeros(cS.nk, cS.nw, cS.aD_new, 'uint16'); for ia = 1:cS.aD_new, for ie = 1:cS.nw, for ik = 1:cS.nk, [~, idx] = min(abs(cS.kGridV - kPolM(ik, 1, ie, ia))); k_prime_idx(ik, ie, ia) = idx; end, end, end, end
function Dist = solve_steady_state_distribution(k_prime_idx, paramS, cS), Dist = zeros(cS.nk, cS.nw, cS.aD_new); Dist(1, :, 1) = paramS.leProb1V'; for ia = 1:(cS.aD_new - 1), Dist_next_age = zeros(cS.nk, cS.nw); for ik = 1:cS.nk, for ie = 1:cS.nw, cond_prob_at_state = Dist(ik, ie, ia); if cond_prob_at_state < 1e-20, continue; end, ik_prime = k_prime_idx(ik, ie, ia); transition_probs_e = paramS.leTrProbM(ie, :); Dist_next_age(ik_prime, :) = Dist_next_age(ik_prime, :) + cond_prob_at_state * transition_probs_e; end, end, Dist(:, :, ia+1) = Dist_next_age; end, end
function [K_agg, C_agg, Tax_agg, Bequest_agg, L_agg] = aggregate_from_dist_nonstochastic(Dist, k_prime_idx, M_sim, Z_ss_norm, paramS, cS), K_agg = 0; C_agg = 0; Tax_agg = 0; Bequest_agg = 0; L_agg = 0; for ia = 1:cS.aD_new, age_pop_share = Z_ss_norm(ia); if age_pop_share < 1e-20, continue; end, for ie = 1:cS.nw, for ik = 1:cS.nk, cond_prob = Dist(ik, ie, ia); if cond_prob < 1e-20, continue; end, mass = cond_prob * age_pop_share; k_now = cS.kGridV(ik); epsilon_val = paramS.leGridV(ie); idx_k_prime = k_prime_idx(ik, ie, ia); k_prime = cS.kGridV(idx_k_prime); [c_val, tax_val, ~] = backout_c_tax_nonstochastic(k_now, k_prime, ia, epsilon_val, M_sim, cS, paramS); C_agg = C_agg + c_val * mass; Tax_agg = Tax_agg + tax_val * mass; K_agg = K_agg + k_prime * mass * cS.s_pathV(ia); Bequest_agg = Bequest_agg + k_prime * mass * (1 - cS.s_pathV(ia)); if ia <= cS.aR_new, individual_labor_supply = cS.ageEffV_new(ia) * epsilon_val; L_agg = L_agg + individual_labor_supply * mass; end, end, end, end, end
function [c_val, tax_val, labor_income_gross] = backout_c_tax_nonstochastic(k_now, k_prime, ia, epsilon_val, M_sim, cS, paramS), tr_per_capita = 0; cpps_decision = 0; cash_inflows_gross = k_now * (1 + M_sim.r_mkt_t) + tr_per_capita; labor_income_gross = 0; if ia <= cS.aR_new, labor_income_gross = M_sim.w_t * cS.ageEffV_new(ia) * epsilon_val; cash_inflows_gross = cash_inflows_gross + labor_income_gross; else, cash_inflows_gross = cash_inflows_gross + M_sim.b_t; end, payg_tax = cS.theta_t * labor_income_gross; labor_tax = cS.tau_l * max(0, labor_income_gross - cpps_decision - payg_tax); capital_tax = cS.tau_k * M_sim.r_mkt_t * k_now; general_taxes = labor_tax + capital_tax; cash_outflows_non_c = payg_tax + general_taxes + k_prime + cpps_decision; c_expend = cash_inflows_gross - cash_outflows_non_c; c_val = max(cS.cFloor, c_expend / (1 + cS.tau_c)); consumption_tax = c_val * cS.tau_c; tax_val = general_taxes + consumption_tax; end