% =========================================================================
% == 调试脚本: debug_full_accounting.m [v3 - 最终修正版]
% == 目的: 修正政府预算闭环，实现完全的会计自洽。
% =========================================================================
clear; close all;
addpath(pwd);
fprintf('=== [完整会计矩阵调试脚本] ===\n\n');

%% 1. 初始化环境 (不变)
fprintf('--- 1. 初始化环境 ---\n');
cS = main_olg_v14_utils.ParameterValues_HuggettStyle();
paramS = struct();
ngrid = 30; cS.nk = ngrid; cS.nkpps = ngrid; cS.nkprime = ngrid; cS.npps = 5;
cS = main_olg_v14_utils.generateGrids(cS);
[Z_path, ~, ~] = main_olg_v14_utils.load_exogenous_paths(cS);
Z_ss_norm = Z_path(:,1);
cS = main_olg_v14_utils.calcaulte_theta_payg_path(cS, false);
[paramS.leLogGridV, paramS.leTrProbM, paramS.leProb1V] = main_olg_v14_utils.EarningProcess_olgm(cS);
paramS.leGridV = exp(paramS.leLogGridV(:));
fprintf('✅ 环境初始化完成。\n\n');

%% 2. 求解稳态均衡 (使用最终修正版的求解器)
fprintf('--- 2. 求解稳态均衡 ---\n');
[ss, eq_found, Dist, k_prime_idx] = solve_steady_state_nonstochastic_final_v2(Z_ss_norm, cS, paramS); 
if ~eq_found, error('稳态求解失败'); end
fprintf('✅ 稳态求解完成。\n\n');

%% 3. 构建详细的会计矩阵
fprintf('--- 3. 构建详细的会计矩阵 ---\n');
% ... (这部分与上一版完全一致，因为ss结构体已经包含了所有正确的值)
Y_gdp = ss.Y_from_production;
Depreciation = ss.I;
Y_ndp = Y_gdp - Depreciation;
NetLaborIncome = ss.w * ss.L;
NetCapitalIncome = ss.r_mkt * ss.K_physical;
FactorIncome_Net = NetLaborIncome + NetCapitalIncome;

fprintf('--- A. 从生产到收入的核算 ---\n');
fprintf('   国内生产总值 (GDP)          : %.6f\n', Y_gdp);
fprintf('   资本折旧 (δK)               : %.6f\n', Depreciation);
fprintf('   国内生产净值 (NDP)          : %.6f\n', Y_ndp);
fprintf('   NDP vs 净要素收入(wL+rK) 残差: %.3e\n\n', Y_ndp - FactorIncome_Net);

% --- B. 国民可支配总收入 (GNDI) 的使用 ---
Consumption_Total = ss.C;
GovPurchase_Total = ss.G;
NationalSaving_from_NDI = Y_ndp - Consumption_Total - GovPurchase_Total;
fprintf('--- B. 国民可支配总收入 (GNDI) 的使用 ---\n');
fprintf('   国民可支配净收入 (GNDI)     : %.6f\n', Y_ndp);
fprintf('   最终消费 (C)                  : %.6f\n', Consumption_Total);
fprintf('   政府购买 (G)                  : %.6f\n', GovPurchase_Total);
fprintf('   => 国民净储蓄 (S = GNDI-C-G)   : %.6f\n\n', NationalSaving_from_NDI);

% --- C. 分部门储蓄检验 ---
cS_agg = cS; cS_agg.theta_t = cS.theta_path(1);
M_final_for_hh = struct('w_t', ss.w, 'r_mkt_t', ss.r_mkt, 'b_t', ss.b);
[~, ~, Tax_Agg, Bequest_Agg, ~, PensionIn_Agg, PensionOut_Agg, ~, ~, ~] = ...
    aggregate_full_accounting(Dist, k_prime_idx, M_final_for_hh, cS_agg, Z_ss_norm, paramS);

PrivateDisposableIncome = NetLaborIncome + NetCapitalIncome + PensionOut_Agg - (Tax_Agg + PensionIn_Agg);
PrivateNetSaving = PrivateDisposableIncome - ss.C;
% 【修正】政府净储蓄 = (税收-政府购买) + 资本转移收入(遗赠)
GovNetSaving = (Tax_Agg - ss.G) + Bequest_Agg;
PensionNetSaving = PensionIn_Agg - PensionOut_Agg;
fprintf('--- C. 分部门储蓄检验 ---\n');
fprintf('   私人净储蓄 (S_pvt)            : %.6f\n', PrivateNetSaving);
fprintf('   政府净储蓄 (S_gov)            : %.6f\n', GovNetSaving);
fprintf('   养老金净储蓄 (S_pen)          : %.6f\n', PensionNetSaving);
fprintf('   加总得到的国民净储蓄         : %.6f\n', PrivateNetSaving + GovNetSaving + PensionNetSaving);
fprintf('   国民净储蓄 vs 分部门加总 残差 : %.3e\n\n', NationalSaving_from_NDI - (PrivateNetSaving + GovNetSaving + PensionNetSaving));

% --- D. 最终的储蓄-投资 恒等式 (净额概念) ---
NetInvestment = ss.I - Depreciation; % 稳态下为0
fprintf('--- D. 最终的储蓄-投资恒等式检验 (净额概念) ---\n');
fprintf('   国民净储蓄 (S_net)            : %.6f\n', NationalSaving_from_NDI);
fprintf('   国民净投资 (I_net)            : %.6f\n', NetInvestment);
fprintf('   S_net vs I_net 残差           : %.3e\n', NationalSaving_from_NDI - NetInvestment);

if abs(NationalSaving_from_NDI - NetInvestment) < 1e-6
    fprintf('\n✅✅✅ 最终胜利！国民净储蓄精确等于国民净投资 (0)！✅✅✅\n');
    fprintf('   这证明了整个经济的会计系统是完全闭环和自洽的。\n');
else
    fprintf('\n❌❌❌ 最后的失败！净储蓄不等于净投资！请仔细检查上面的每一个账户！❌❌❌\n');
end

%% ========================================================================
%  == 所有需要的辅助函数
%  ========================================================================
function [ss, eq_found, Dist, k_prime_idx] = solve_steady_state_nonstochastic_final(Z_ss_norm, cS, paramS)
    K_physical_guess = 1.2; max_iter_K = 100; tol_K = 1e-8; damp_K = 0.4; eq_found = false; A_ss = cS.A; theta_ss = cS.theta_path(1);
    e_dist_by_age = zeros(cS.aD_new, cS.nw); e_dist_by_age(1, :) = paramS.leProb1V';
    for ia = 1:(cS.aD_new - 1), e_dist_by_age(ia+1, :) = e_dist_by_age(ia, :) * paramS.leTrProbM; end
    mean_e_by_age = e_dist_by_age * paramS.leGridV(:); L_ss = 0;
    for ia = 1:cS.aR_new, L_ss = L_ss + Z_ss_norm(ia) * cS.ageEffV_new(ia) * mean_e_by_age(ia); end
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
        [K_pvt_model, ~, ~, ~, ~, ~, ~, ~, ~, ~] = aggregate_full_accounting(Dist, k_prime_idx, M_ss_for_hh, cS_ss, Z_ss_norm, paramS);
        K_physical_model = K_pvt_model; K_error = K_physical_guess - K_physical_model;
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
    [~, kPolM_final, ~, ~] = main_olg_v14_utils.HHSolution_VFI_Huggett(M_final_for_hh, 0, paramS, cS_ss);
    k_prime_idx_final = get_policy_index_matrix(kPolM_final, cS_ss);
    Dist_final = solve_steady_state_distribution(k_prime_idx_final, paramS, cS_ss);
    [K_pvt_final, C_final, Tax_final, ~, ~, ~, ~, ~, ~, ~] = aggregate_full_accounting(Dist_final, k_prime_idx_final, M_final_for_hh, cS_ss, Z_ss_norm, paramS);
    
    % 【最终修正】G只由当期税收决定
    G_final = Tax_final;
    I_final = cS.ddk * K_pvt_final;
    Y_final_from_expenditure = C_final + I_final + G_final;

    ss = struct(); ss.K_physical = K_pvt_final; ss.L = L_ss; ss.Y = Y_final_from_expenditure;
    ss.w = M_final_prices.w_t; ss.r_mkt = M_final_prices.r_mkt_t; ss.b = M_final_for_hh.b_t;
    ss.C = C_final; ss.I = I_final; ss.G = G_final; ss.Y_from_production = M_final_prices.Y_t;
    Dist = Dist_final; k_prime_idx = k_prime_idx_final;
end

% =========================================================================
% == 函数: solve_steady_state_nonstochastic_final_v2
% == 目的: 使用正确的经济均衡条件 (S_net = 0) 来求解模型。
% == 方法: 不再使用简单的资本不动点迭代，而是使用求根算法(fzero)
% ==       来寻找使国民净储蓄为零的均衡资本水平。
% =========================================================================
function [ss, eq_found, Dist, k_prime_idx] = solve_steady_state_nonstochastic_final_v2(Z_ss_norm, cS, paramS)
    fprintf('\n--- 开始求解稳态 (使用 S_net = 0 作为均衡条件) ---\n');

    % 定义一个函数，该函数输入K，输出国民净储蓄S_net
    % fzero 将寻找使这个函数输出为0的K
    excess_saving_function = @(K_guess) calculate_net_saving(K_guess, Z_ss_norm, cS, paramS);

    % 使用fzero求根器来寻找均衡资本 K_eq
    % 我们需要一个包含解的区间，例如 [0.1, 5.0]
    k_bracket = [0.1, 5.0]; 
    options = optimset('TolX', 1e-9, 'Display', 'iter');
    
    try
        [K_eq, fval, exitflag] = fzero(excess_saving_function, k_bracket, options);
        eq_found = (exitflag > 0);
    catch ME
        warning('fzero 求解失败: %s', ME.message);
        K_eq = NaN;
        eq_found = false;
    end

    if ~eq_found
        warning('主循环(fzero)未能找到均衡解');
        ss=struct(); Dist=[]; k_prime_idx=[]; return;
    end
    
    fprintf('✅ 均衡收敛！ 均衡资本 K_eq = %.6f, 最终净储蓄 S_net = %.3e\n', K_eq, fval);

    % --- 收敛后，使用最终的均衡资本K_eq重新计算所有宏观量以确保自洽 ---
    [~, ss, Dist, k_prime_idx] = calculate_net_saving(K_eq, Z_ss_norm, cS, paramS);
end

% =========================================================================
% == 辅助函数: calculate_net_saving
% == 输入: 资本K_guess 和其他参数
% == 输出: 对应的国民净储蓄 S_net, 以及包含所有宏观量的结构体ss (可选)
% =========================================================================
function [S_net, ss, Dist, k_prime_idx] = calculate_net_saving(K_guess, Z_ss_norm, cS, paramS)
    % --- 1. 基于K_guess，计算所有价格和微观决策 ---
    A_ss = cS.A; theta_ss = cS.theta_path(1);
    
    % 计算L_ss (这段代码已经是正确的)
    e_dist_by_age = zeros(cS.aD_new, cS.nw); e_dist_by_age(1, :) = paramS.leProb1V';
    for ia = 1:(cS.aD_new - 1), e_dist_by_age(ia+1, :) = e_dist_by_age(ia, :) * paramS.leTrProbM; end
    mean_e_by_age = e_dist_by_age * paramS.leGridV(:); L_ss = 0;
    for ia = 1:cS.aR_new, L_ss = L_ss + Z_ss_norm(ia) * cS.ageEffV_new(ia) * mean_e_by_age(ia); end

    M_prices = main_olg_v14_utils.get_prices_at_t(K_guess, L_ss, A_ss, cS);
    M_for_hh = M_prices;
    mass_retirees_ss = sum(Z_ss_norm(cS.aR_new+1:end));
    total_wage_bill = M_prices.w_t * L_ss;
    if mass_retirees_ss > 1e-9, M_for_hh.b_t = (theta_ss * total_wage_bill) / mass_retirees_ss; else, M_for_hh.b_t = 0; end
    
    cS_ss = cS; cS_ss.pps_active = false; cS_ss.theta_t = theta_ss;
    if ~cS_ss.pps_active, cS_ss.nkpps = 1; cS_ss.npps = 1; cS_ss = main_olg_v14_utils.generateGrids(cS_ss); end
    
    [~, kPolM, ~, ~] = main_olg_v14_utils.HHSolution_VFI_Huggett(M_for_hh, 0, paramS, cS_ss);
    k_prime_idx = get_policy_index_matrix(kPolM, cS_ss);
    Dist = solve_steady_state_distribution(k_prime_idx, paramS, cS_ss);
    
    % --- 2. 聚合得到 C 和 G ---
    [~, C_final, Tax_final, ~, ~, ~, ~, ~, ~, ~] = aggregate_full_accounting(Dist, k_prime_idx, M_for_hh, cS_ss, Z_ss_norm, paramS);
    
    % 使用正确的政府预算规则
    G_final = Tax_final;
    
    % --- 3. 计算国民净储蓄 S_net ---
    GDP = M_prices.Y_t;
    Depreciation = cS.ddk * K_guess;
    NDP = GDP - Depreciation;
    S_net = NDP - C_final - G_final;
    
    % --- 4. (可选) 打包所有宏观量，供最终收敛后使用 ---
    if nargout > 1
        I_final = Depreciation; % 稳态下总投资=折旧
        Y_exp = C_final + I_final + G_final;
        ss = struct();
        ss.K_physical = K_guess;
        ss.L = L_ss;
        ss.Y = Y_exp;
        ss.Y_from_production = GDP;
        ss.C = C_final;
        ss.I = I_final;
        ss.G = G_final;
        ss.w = M_prices.w_t;
        ss.r_mkt = M_prices.r_mkt_t;
        ss.b = M_for_hh.b_t;
    end
end

function k_prime_idx = get_policy_index_matrix(kPolM, cS), k_prime_idx = zeros(cS.nk, cS.nw, cS.aD_new, 'uint16'); for ia = 1:cS.aD_new, for ie = 1:cS.nw, for ik = 1:cS.nk, [~, idx] = min(abs(cS.kGridV - kPolM(ik, 1, ie, ia))); k_prime_idx(ik, ie, ia) = idx; end, end, end, end
function Dist = solve_steady_state_distribution(k_prime_idx, paramS, cS), Dist = zeros(cS.nk, cS.nw, cS.aD_new); Dist(1, :, 1) = paramS.leProb1V'; for ia = 1:(cS.aD_new - 1), Dist_next_age = zeros(cS.nk, cS.nw); for ik = 1:cS.nk, for ie = 1:cS.nw, cond_prob_at_state = Dist(ik, ie, ia); if cond_prob_at_state < 1e-20, continue; end, ik_prime = k_prime_idx(ik, ie, ia); transition_probs_e = paramS.leTrProbM(ie, :); Dist_next_age(ik_prime, :) = Dist_next_age(ik_prime, :) + cond_prob_at_state * transition_probs_e; end, end, Dist(:, :, ia+1) = Dist_next_age; end, end
function [K_agg, C_agg, Tax_agg, Bequest_agg, L_agg, PensionIn_agg, PensionOut_agg, LaborTax_agg, CapitalTax_agg, ConsumpTax_agg] = aggregate_full_accounting(Dist, k_prime_idx, M_sim, cS, Z_ss_norm, paramS), K_agg=0; C_agg=0; Tax_agg=0; Bequest_agg=0; L_agg=0; PensionIn_agg=0; PensionOut_agg=0; LaborTax_agg=0; CapitalTax_agg=0; ConsumpTax_agg=0; for ia=1:cS.aD_new, age_pop_share = Z_ss_norm(ia); if age_pop_share<1e-20, continue; end; for ie=1:cS.nw, for ik=1:cS.nk, cond_prob=Dist(ik,ie,ia); if cond_prob<1e-20, continue; end; mass=cond_prob*age_pop_share; k_now=cS.kGridV(ik); epsilon_val=paramS.leGridV(ie); idx_k_prime=k_prime_idx(ik,ie,ia); k_prime=cS.kGridV(idx_k_prime); [c_val,tax_val,~,payg_tax,labor_tax,capital_tax,consump_tax,pension_benefit]=backout_full_accounting(k_now,k_prime,ia,epsilon_val,M_sim,cS,paramS); C_agg=C_agg+c_val*mass; Tax_agg=Tax_agg+tax_val*mass; K_agg=K_agg+k_prime*mass*cS.s_pathV(ia); Bequest_agg=Bequest_agg+k_prime*mass*(1-cS.s_pathV(ia)); if ia<=cS.aR_new, L_agg=L_agg+(cS.ageEffV_new(ia)*epsilon_val)*mass; end; PensionIn_agg=PensionIn_agg+payg_tax*mass; PensionOut_agg=PensionOut_agg+pension_benefit*mass; LaborTax_agg=LaborTax_agg+labor_tax*mass; CapitalTax_agg=CapitalTax_agg+capital_tax*mass; ConsumpTax_agg=ConsumpTax_agg+consump_tax*mass; end, end, end, end
function [c_val, tax_val, labor_income_gross, payg_tax, labor_tax, capital_tax, consumption_tax, pension_benefit] = backout_full_accounting(k_now, k_prime, ia, epsilon_val, M_sim, cS, paramS), tr_per_capita=0; cpps_decision=0; cash_inflows_gross=k_now*(1+M_sim.r_mkt_t)+tr_per_capita; labor_income_gross=0; pension_benefit=0; if ia<=cS.aR_new, labor_income_gross=M_sim.w_t*cS.ageEffV_new(ia)*epsilon_val; cash_inflows_gross=cash_inflows_gross+labor_income_gross; else, pension_benefit=M_sim.b_t; cash_inflows_gross=cash_inflows_gross+pension_benefit; end; payg_tax=cS.theta_t*labor_income_gross; labor_tax=cS.tau_l*max(0,labor_income_gross-cpps_decision-payg_tax); capital_tax=cS.tau_k*M_sim.r_mkt_t*k_now; cash_outflows_non_c=payg_tax+labor_tax+capital_tax+k_prime+cpps_decision; c_expend=cash_inflows_gross-cash_outflows_non_c; c_val=max(cS.cFloor, c_expend/(1+cS.tau_c)); consumption_tax=c_val*cS.tau_c; tax_val=labor_tax+capital_tax+consumption_tax; end

