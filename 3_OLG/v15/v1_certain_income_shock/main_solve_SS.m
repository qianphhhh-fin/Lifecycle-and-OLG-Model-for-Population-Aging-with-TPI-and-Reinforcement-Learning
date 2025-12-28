% =========================================================================
% == SCRIPT: main_solve_SS.m [vPaper.2 - 重大冲击模型]
% == 目的: 在遗赠税模型的基础上，加入确定性的生命周期重大支出冲击，
% ==         以增强家庭的预防性储蓄动机，并观察对资本积累(K/Y)的影响。
% == 核心修正:
% == 1. 调用新的VFI求解器和核算函数，这些函数内部处理了重大冲击支出。
% == 2. 在国民账户中，总消费(C)被明确分解为有效用消费和冲击性支出。
% == 3. 在最终报告中，展示冲击支出的规模，并计算最终的资本产出比(K/Y)。
% =========================================================================
clear; close all;                           % 清空工作空间和关闭图形窗口
addpath(pwd);                               % 将当前目录添加到MATLAB路径
fprintf('=== [重大冲击版] 稳态模型求解与国民核算报告 (遗赠税模型) ===\n\n');

%% 1. 初始化环境
fprintf('--- 1. 初始化环境并生成模拟所需输入 ---\n');

% [核心修正] 调用包含冲击参数的新版参数设置函数
cS = main_olg_v15_utils.ParameterValues_HuggettStyle(); 
paramS = struct();
fprintf('   稳态财政设定：政府预算平衡 (G=T), 意外遗赠被政府100%%征收。\n');
cS.gov_debt_frac_Y = 0;

ngrid = 50; cS.nk = ngrid; cS.nkpps = ngrid; cS.nkprime = ngrid; cS.npps = 5;
cS = main_olg_v15_utils.generateGrids(cS);
[~, ~, ~] = main_olg_v15_utils.load_exogenous_paths(cS);
cS = main_olg_v15_utils.calcaulte_theta_payg_path(cS, false);
[paramS.leLogGridV, paramS.leTrProbM, paramS.leProb1V] = main_olg_v15_utils.EarningProcess_olgm(cS);
paramS.leGridV = exp(paramS.leLogGridV(:));
eIdxM = main_olg_v15_utils.LaborEndowSimulation_olgm_AgeGroup(cS, paramS);

% 1. 基于存活率计算稳态年龄人口质量
age_mass = ones(cS.aD_new, 1);
for ia = 1:(cS.aD_new - 1)
    age_mass(ia+1) = age_mass(ia) * cS.s_pathV(ia);
end

% 2. 归一化，得到每个年龄组占总人口的比例
Z_ss_model_implied = age_mass / sum(age_mass);

% 3. 用这个内生的、自洽的分布，替换掉从外部加载的 Z_path(:,1)
Z_ss_norm = Z_ss_model_implied;
fprintf('   已生成模型内生稳态人口分布，将用此分布求解。\n');

%% 2. 求解稳态 (资本市场出清法)
fprintf('\n--- 2. 调用稳态求解器 (以 资产供给 = 资产需求 为目标) ---\n');

[ss, eq_found, Dist, k_prime_idx] = solve_steady_state_iter_K(Z_ss_norm, cS, paramS);
if ~eq_found, error('资本市场出清法求解失败'); end

%% 3. 结果验证与展示
fprintf('\n\n--- 3. 结果验证与展示 ---\n');
fprintf('✅✅✅ 求解成功！均衡资本 K = %.8f ✅✅✅\n', ss.K_physical);

% [核心修正] 调用新版展示函数
display_national_accounts_bequest_tax(ss, Dist, k_prime_idx, cS, paramS);


% =========================================================================
% ==                        函数定义部分                                 ==
% =========================================================================

function [K_agg,C_utility_agg,Tax_agg,Bequest_tax_agg,L_agg,PensionIn_agg,PensionOut_agg,ShockExp_agg]=aggregate_full_accounting_BequestTax(Dist,k_prime_idx,M_sim,cS,paramS)
    % [vPaper.2 - 重大冲击版]
    % 核心变更: 分别聚合有效用消费(C_utility_agg)和无用的冲击支出(ShockExp_agg)。
    if ~isfield(cS, 'theta_t'), cS.theta_t = cS.theta_path(1); end
    K_agg=0; C_utility_agg=0; Tax_agg=0; Bequest_tax_agg=0; L_agg=0; PensionIn_agg=0; PensionOut_agg=0; ShockExp_agg=0;
    
    for ia=1:cS.aD_new
        for ie=1:cS.nw
            for ik=1:cS.nk
        mass = Dist(ik,ie,ia);
        if mass < 1e-20, continue; end
        k_now=cS.kGridV(ik); epsilon_val=paramS.leGridV(ie);
        idx_k_prime=k_prime_idx(ik,ie,ia); k_prime=cS.kGridV(idx_k_prime);
        
        % [核心修正] 调用新版核算函数，获取冲击支出
        [c_val,tax_val,labor_income_gross,payg_tax,~,~,~,pension_benefit,shock_exp_val] = ...
            backout_full_accounting_BequestTax(k_now,k_prime,ia,epsilon_val,M_sim,cS);
        
        C_utility_agg = C_utility_agg + c_val * mass;
        ShockExp_agg = ShockExp_agg + shock_exp_val * mass; % 聚合冲击支出
        Tax_agg = Tax_agg + tax_val * mass;
        
        prob_survive = cS.s_pathV(ia);
        K_agg = K_agg + k_prime * mass * prob_survive;
        prob_death = 1 - prob_survive;
        Bequest_tax_agg = Bequest_tax_agg + k_prime * mass * prob_death;
        
        if ia<=cS.aR_new, L_agg = L_agg + (cS.ageEffV_new(ia)*epsilon_val) * mass; end
        PensionIn_agg = PensionIn_agg + payg_tax * mass;
        PensionOut_agg = PensionOut_agg + pension_benefit * mass;
    end, end, end
end

function [c_val,tax_val,labor_income_gross,payg_tax,labor_tax,capital_tax,consumption_tax,pension_benefit,shock_expenditure_val]=backout_full_accounting_BequestTax(k_now,k_prime,ia,epsilon_val,M_sim,cS)
    % [vPaper.2 - 重大冲击版]
    % 核心变更: 精确复制VFI求解器中的逻辑，计算并返回重大冲击支出。
    cpps_decision=0;
    if ~isfield(cS, 'theta_t'), cS.theta_t = cS.theta_path(1); end
    
    market_return_factor = 1 + M_sim.r_mkt_t;
    capital_income = k_now * M_sim.r_mkt_t;
    
    % 1. 计算总流入
    cash_inflows_gross = k_now * market_return_factor;
    labor_income_gross = 0; pension_benefit = 0;
    if ia <= cS.aR_new
        labor_income_gross = M_sim.w_t * cS.ageEffV_new(ia) * epsilon_val;
        cash_inflows_gross = cash_inflows_gross + labor_income_gross;
    else
        pension_benefit = M_sim.b_t;
        cash_inflows_gross = cash_inflows_gross + pension_benefit;
    end
    
    % 2. 计算各项税收
    payg_tax = cS.theta_t * labor_income_gross;
    capital_tax = cS.tau_k * capital_income;
    labor_tax = cS.tau_l * max(0, labor_income_gross - cpps_decision - payg_tax);
    
    % 3. 计算税后净现金流
    net_cash_before_shock = cash_inflows_gross - (payg_tax + labor_tax + capital_tax);

    % 4. [核心修正] 计算重大冲击支出
    shock_expenditure_val = 0;
    is_shock_age = ismember(ia, [cS.shock_age_young_idx, cS.shock_age_old_idx]);
    if cS.shock_active && is_shock_age
        shock_expenditure_val = cS.shock_frac * net_cash_before_shock;
    end

    % 5. 计算常规消费支出 (基于剩余现金)
    net_cash_after_shock = net_cash_before_shock - shock_expenditure_val;
    c_expend = net_cash_after_shock - k_prime;
    
    c_val = max(cS.cFloor, c_expend / (1 + cS.tau_c));
    consumption_tax = c_val * cS.tau_c;
    tax_val = labor_tax + capital_tax + consumption_tax;
end

function [cPol_age_q, kPol_age, cPpsPol_age_choice, val_age] = HHSolutionByAge_VFI_Vectorized_BequestTax(a_idx, vPrime_kkppse_next, M_age, b_age_val, paramS_age, cS)
    % =========================================================================
    % == 函数: HHSolutionByAge_VFI_Vectorized_BequestTax [vPaper.2 - 重大冲击版]
    % == 核心变更:
    % == 1. 在计算完可用于消费和储蓄的净现金后，首先扣除确定性的重大支出。
    % == 2. 家庭的最优决策是基于【扣除冲击后】的剩余资源进行的。
    % == 3. 此支出不计入效用函数，也不纳消费税。
    % =========================================================================

    % --- 0. 初始化 ---
    val_age = -1e20 * ones(cS.nk, cS.nkpps, cS.nw);
    cPol_age_q = zeros(cS.nk, cS.nkpps, cS.nw);
    kPol_age = zeros(cS.nk, cS.nkpps, cS.nw);
    cPpsPol_age_choice = zeros(cS.nk, cS.nkpps, cS.nw);

    % --- 1. 最后生命周期 (逻辑不变) ---
    if a_idx == cS.aD_new
         [K_grid, ~] = ndgrid(cS.kGridV, cS.kppsGridV);
         pretax_non_capital_income = b_age_val;
         capital_tax = cS.tau_k .* M_age.r_mkt_t .* K_grid;
         total_resources = K_grid .* (1 + M_age.r_mkt_t) + pretax_non_capital_income - capital_tax;
         
         final_c_expenditure = total_resources;
         final_bequest = zeros(size(total_resources));
         final_c = max(cS.cFloor, final_c_expenditure / (1 + cS.tau_c));
         [~, util_c] = main_olg_v15_utils.CES_utility(final_c, cS.sigma, cS);
         final_v = util_c;
         
         for ie = 1:cS.nw
             cPol_age_q(:,:,ie) = final_c;
             val_age(:,:,ie) = final_v;
             kPol_age(:,:,ie) = final_bequest;
             cPpsPol_age_choice(:,:,ie) = 0;
         end
        return;
    end
    
    % --- 2. 准备期望价值插值器 (逻辑不变) ---
    effective_discount_factor = (cS.beta ^ cS.time_Step) * cS.s_pathV(a_idx); 
    EV_matrix = zeros(cS.nk, cS.nkpps, cS.nw);
    for ie_current = 1:cS.nw
        transition_probs = paramS_age.leTrProbM(ie_current, :);
        if isempty(vPrime_kkppse_next), EV_matrix = zeros(cS.nk, cS.nkpps, cS.nw); break; end
        vPrime_reshaped = reshape(vPrime_kkppse_next, [cS.nk * cS.nkpps, cS.nw]);
        EV_slice = vPrime_reshaped * transition_probs';
        EV_matrix(:, :, ie_current) = reshape(EV_slice, [cS.nk, cS.nkpps]);
    end
    EV_interpolants = cell(cS.nw, 1);
    is_pps_disabled = (cS.nkpps == 1);
    for ie_current=1:cS.nw
        if is_pps_disabled, EV_interpolants{ie_current}=griddedInterpolant(cS.kGridV,squeeze(EV_matrix(:,1,ie_current)),'pchip','none');
        else, EV_interpolants{ie_current}=griddedInterpolant({cS.kGridV, cS.kppsGridV},EV_matrix(:,:,ie_current),'pchip','none'); end
    end
    
    % --- 3. 准备当期价格和回报率 (逻辑不变) ---
    market_return_factor = 1 + M_age.r_mkt_t;

    % --- 4. 逐状态求解主循环 ---
    for ie = 1:cS.nw, epsilon_state = paramS_age.leGridV(ie); ev_interpolant = EV_interpolants{ie};
        for ik = 1:cS.nk, for ikpps = 1:cS.nkpps
            k_state = cS.kGridV(ik);
            best_val = -1e20; best_k_prime = cS.kMin; best_c = cS.cFloor;
            
            % 4a. 计算各项收入和税收 (逻辑不变)
            labor_income_gross = 0; pension_benefit = 0;
            capital_income = k_state * M_age.r_mkt_t; capital_tax = cS.tau_k * capital_income;
            if a_idx <= cS.aR_new
                labor_income_gross = M_age.w_t*cS.ageEffV_new(a_idx)*epsilon_state;
                payg_tax = cS.theta_t*labor_income_gross;
                labor_tax = cS.tau_l*max(0, labor_income_gross - payg_tax);
            else
                pension_benefit = b_age_val; payg_tax = 0; labor_tax = 0;
            end
            
            % 4b. 构建【扣税后】的净现金流
            cash_inflow_gross = k_state * market_return_factor + labor_income_gross + pension_benefit;
            net_cash_for_c_k_prime_before_shock = cash_inflow_gross - (payg_tax + labor_tax + capital_tax);
            
            % 4c. [核心修正] 扣除确定性重大支出
            shock_expenditure = 0;
            is_shock_age = ismember(a_idx, [cS.shock_age_young_idx, cS.shock_age_old_idx]);
            if cS.shock_active && is_shock_age
                shock_expenditure = cS.shock_frac * net_cash_for_c_k_prime_before_shock;
            end
            net_cash_for_c_k_prime_after_shock = net_cash_for_c_k_prime_before_shock - shock_expenditure;
            
            % 4d. 基于【剩余现金】进行最优储蓄决策
            k_prime_max = net_cash_for_c_k_prime_after_shock - cS.cFloor * (1 + cS.tau_c);
            if k_prime_max < cS.kMin, k_prime_grid = [cS.kMin]; else, k_prime_grid = unique([linspace(cS.kMin, k_prime_max, cS.nkprime), cS.kGridV(cS.kGridV <= k_prime_max)']); end
            
            for k_prime_choice = k_prime_grid
                c_expend = net_cash_for_c_k_prime_after_shock - k_prime_choice;
                if c_expend < cS.cFloor*(1+cS.tau_c), continue; end
                c_choice = c_expend / (1 + cS.tau_c);
                [~,util] = main_olg_v15_utils.CES_utility(c_choice, cS.sigma, cS);
                
                k_pps_prime = 0; % 在此简化模型中，pps为0
                if is_pps_disabled, ev = ev_interpolant(k_prime_choice); else, ev = ev_interpolant(k_prime_choice, k_pps_prime); end
                if isnan(ev), ev = -1e10; end
                
                current_val = util + effective_discount_factor * ev;
                if current_val > best_val, best_val = current_val; best_c = c_choice; best_k_prime = k_prime_choice; end
            end
            
            % 4e. 存储最优决策
            if isinf(best_val)||isnan(best_val), val_age(ik,ikpps,ie)=-1e20; kPol_age(ik,ikpps,ie)=cS.kMin; cPpsPol_age_choice(ik,ikpps,ie)=0; cPol_age_q(ik,ikpps,ie)=cS.cFloor;
            else, val_age(ik,ikpps,ie)=best_val; kPol_age(ik,ikpps,ie)=best_k_prime; cPpsPol_age_choice(ik,ikpps,ie)=0; cPol_age_q(ik,ikpps,ie)=best_c; end
    end,end,end
end


function [ss, eq_found, Dist, k_prime_idx] = solve_steady_state_iter_K(Z_ss_norm, cS, paramS)
    % [v2.1 - 逻辑修正] fzero 调用逻辑修正
    
    % fzero_wrapper_K 直接计算并返回资本市场的误差
    capital_market_wrapper = @(K_guess) fzero_wrapper_K(K_guess, Z_ss_norm, cS, paramS);
    
    k_bracket=[0.05, 10]; 
    options=optimset('TolX',1e-15,'Display','iter','MaxIter',500);
    [K_eq, fval, exitflag] = fzero(capital_market_wrapper, k_bracket, options); 
    
    eq_found = (exitflag > 0);
    if ~eq_found
        warning('iter_K 未能找到均衡解'); 
        ss=struct(); Dist=[]; k_prime_idx=[]; 
        return; 
    end
    
    fprintf('   iter_K 收敛！均衡资本 K_eq = %.8f, 最终K误差 = %.3e\n', K_eq, fval);
    
    % [核心修正] 使用修正后的输出顺序来获取所有最终结果
    % 我们现在请求所有5个输出，并正确地将它们分配给变量
    [~, ~, ss, Dist, k_prime_idx] = calculate_net_saving(K_eq, Z_ss_norm, cS, paramS);
end

function K_error_out = fzero_wrapper_K(K_guess, Z_ss_norm, cS, paramS)
    % [v2.1 - 逻辑修正]
    % 此函数现在直接调用 calculate_net_saving，并计算误差。
    % 它只关心第一个输出参数 K_model。
    K_model = calculate_net_saving(K_guess, Z_ss_norm, cS, paramS);
    K_error_out = K_guess - K_model;
end

function [K_model_out, S_net, ss, Dist, k_prime_idx] = calculate_net_saving(K_guess, Z_ss_norm, cS, paramS)
    % [v2.1 - 逻辑修正]
    % 核心修正: 调整了输出参数的顺序。
    % 1. K_model_out (模型计算出的资本需求) 现在是第一个输出，专门服务于fzero求解器。
    % 2. S_net 和完整的 ss 结构体等作为后续输出，供最终计算使用。
    % 3. 增加了 nargout > 1 的判断，避免在求解循环中进行不必要的详细计算。

    if K_guess <= 0, K_guess = 1e-8; end
    
    % --- 价格计算 (不变) ---
    A_ss=cS.A; theta_ss=cS.theta_path(1); e_dist_by_age=zeros(cS.aD_new,cS.nw);
    e_dist_by_age(1,:)=paramS.leProb1V'; for ia=1:(cS.aD_new-1),e_dist_by_age(ia+1,:)=e_dist_by_age(ia,:)*paramS.leTrProbM;end
    mean_e_by_age=e_dist_by_age*paramS.leGridV(:); L_ss=0;
    for ia=1:cS.aR_new,L_ss=L_ss+Z_ss_norm(ia)*cS.ageEffV_new(ia)*mean_e_by_age(ia);end
    M_prices=main_olg_v15_utils.get_prices_at_t(K_guess,L_ss,A_ss,cS);
    M_for_hh=M_prices;
    mass_retirees_ss=sum(Z_ss_norm(cS.aR_new+1:end));
    total_wage_bill=M_prices.w_t*L_ss;
    if mass_retirees_ss>1e-9,M_for_hh.b_t=(theta_ss*total_wage_bill)/mass_retirees_ss;else,M_for_hh.b_t=0;end
    cS_ss=cS;cS_ss.pps_active=false;cS_ss.theta_t=theta_ss;
    if ~cS_ss.pps_active,cS_ss.nkpps=1;cS_ss.npps=1;cS_ss=main_olg_v15_utils.generateGrids(cS_ss);end
    
    % --- 家庭问题求解 (不变) ---
    [~, kPolM, ~, ~]=HHSolution_VFI_Huggett_BequestTax(M_for_hh, paramS, cS_ss);
    k_prime_idx=get_policy_index_matrix(kPolM, cS_ss);
    Dist=solve_steady_state_distribution(k_prime_idx, paramS, cS_ss, Z_ss_norm);
    
    % --- 聚合 (不变) ---
    [K_model_out, C_utility_final, Tax_final, Bequest_tax_final, L_agg_check, ~, ~, ShockExp_final] = ...
        aggregate_full_accounting_BequestTax(Dist, k_prime_idx, M_for_hh, cS_ss, paramS);
    
    % 默认输出初始化
    S_net = [];
    ss = struct();

    % [核心修正] 只有在需要详细结果时 (即非fzero求解循环中)，才计算和打包ss结构体
    if nargout > 1
        if abs(L_agg_check - L_ss) > 1e-8, warning('聚合劳动供给(%.4f)与初始计算(%.4f)不一致', L_agg_check, L_ss); end
        
        C_total_final = C_utility_final + ShockExp_final;
        G_final = Tax_final + Bequest_tax_final; 
        
        GDP=M_prices.Y_t; 
        Depreciation=cS.ddk*K_guess; 
        NDP=GDP-Depreciation;
        S_net=NDP - C_total_final - G_final;
        
        I_final_gross = Depreciation;
        Y_exp = C_total_final + I_final_gross + G_final;
        
        ss.K_physical=K_guess;
        ss.L=L_ss;
        ss.Y=Y_exp;
        ss.Y_from_production=GDP;
        ss.C_total=C_total_final; 
        ss.C_utility=C_utility_final; 
        ss.C_shock=ShockExp_final;
        ss.I=I_final_gross;
        ss.G=G_final;
        ss.w=M_prices.w_t;
        ss.r_mkt=M_prices.r_mkt_t;
        ss.b=M_for_hh.b_t;
        ss.Bequest_tax = Bequest_tax_final; 
        ss.Regular_tax = Tax_final;
    end
end
function [cPolM, kPolM, cPpsPolM, valM] = HHSolution_VFI_Huggett_BequestTax(M_vfi, paramS_vfi, cS_vfi)
    valM = -Inf(cS_vfi.nk, cS_vfi.nkpps, cS_vfi.nw, cS_vfi.aD_new);
    cPolM = zeros(cS_vfi.nk, cS_vfi.nkpps, cS_vfi.nw, cS_vfi.aD_new);
    kPolM = zeros(cS_vfi.nk, cS_vfi.nkpps, cS_vfi.nw, cS_vfi.aD_new);
    cPpsPolM = zeros(cS_vfi.nk, cS_vfi.nkpps, cS_vfi.nw, cS_vfi.aD_new);
    bV_payg_vfi = zeros(1, cS_vfi.aD_new);
    if cS_vfi.aR_new < cS_vfi.aD_new
        retirement_indices = (cS_vfi.aR_new+1):cS_vfi.aD_new;
        bV_payg_vfi(retirement_indices) = M_vfi.b_t;
    end
    for a_idx = cS_vfi.aD_new : -1 : 1
        vPrime_kkppse_next = [];
        if a_idx < cS_vfi.aD_new
            vPrime_kkppse_next = valM(:,:,:,a_idx+1);
        end
        % [核心修正] 调用包含冲击逻辑的新版VFI求解器
        [cPolM(:,:,:,a_idx), kPolM(:,:,:,a_idx), cPpsPolM(:,:,:,a_idx), valM(:,:,:,a_idx)] = ...
            HHSolutionByAge_VFI_Vectorized_BequestTax(a_idx, vPrime_kkppse_next, M_vfi, bV_payg_vfi(a_idx), paramS_vfi, cS_vfi);
    end
end

function k_prime_idx=get_policy_index_matrix(kPolM,cS)
    k_prime_idx=zeros(cS.nk,cS.nw,cS.aD_new,'uint16');
    for ia=1:cS.aD_new,for ie=1:cS.nw,for ik=1:cS.nk
        [~,idx]=min(abs(cS.kGridV-kPolM(ik,1,ie,ia)));k_prime_idx(ik,ie,ia)=idx;
    end,end,end
end

function Dist=solve_steady_state_distribution(k_prime_idx,paramS,cS,Z_ss_norm)
    Dist = zeros(cS.nk, cS.nw, cS.aD_new);
    dist_newborn_ke = zeros(cS.nk, cS.nw);
    dist_newborn_ke(1, :) = paramS.leProb1V';
    Dist(:, :, 1) = dist_newborn_ke * Z_ss_norm(1);
    for ia = 1:(cS.aD_new - 1)
        dist_ia_ke = Dist(:,:,ia);
        dist_ia_plus_1_ke = zeros(cS.nk, cS.nw);
        for ik = 1:cS.nk, for ie = 1:cS.nw
            mass_at_state = dist_ia_ke(ik, ie);
            if mass_at_state < 1e-20, continue; end
            ik_prime = k_prime_idx(ik, ie, ia);
            transition_probs_e = paramS.leTrProbM(ie, :);
            mass_surviving = mass_at_state * cS.s_pathV(ia);
            dist_ia_plus_1_ke(ik_prime, :) = dist_ia_plus_1_ke(ik_prime, :) + mass_surviving * transition_probs_e;
        end, end
        Dist(:,:,ia+1) = dist_ia_plus_1_ke;
    end
    if abs(sum(Dist, 'all') - 1.0) > 1e-6
         warning('稳态联合分布总和(%.8f)不为1，可能存在会计不一致。', sum(Dist, 'all'));
    end
end

function display_national_accounts_bequest_tax(ss, Dist, k_prime_idx, cS, paramS)
    fprintf('\n\n========================================================================\n');
    fprintf('===     国民经济核算详细报告 (重大冲击模型版)     ===\n');
    fprintf('========================================================================\n');
    
    % --- 0. 重新聚合所有流量，确保一致性 ---
    M_sim = struct('r_mkt_t', ss.r_mkt, 'w_t', ss.w, 'b_t', ss.b);
    % [核心修正] 调用新版聚合函数
    [~, C_utility_agg, Tax_agg, Bequest_tax_agg, ~, PensionIn_agg, PensionOut_agg, ShockExp_agg] = ...
        aggregate_full_accounting_BequestTax(Dist, k_prime_idx, M_sim, cS, paramS);
    
    % --- 1. 计算宏观总量 ---
    Y_prod = ss.Y_from_production;
    K_supply = ss.K_physical;
    Depreciation = cS.ddk * K_supply;
    G_agg = Tax_agg + Bequest_tax_agg;
    I_agg_gross = Depreciation;
    C_total_agg = C_utility_agg + ShockExp_agg; % 总消费
    Y_exp_actual = C_total_agg + I_agg_gross + G_agg;
    
    % --- A. 宏观产出与支出核算 ---
    fprintf('--- A. 宏观产出与支出核算 ---\n');
    fprintf('   生产法 GDP (Y_prod):         %.6f\n', Y_prod);
    fprintf('   支出法 GDP (C_total+I+G):    %.6f\n', Y_exp_actual);
    fprintf('   ------------------------------------\n');
    fprintf('   核算误差 (Y_exp - Y_prod):     %.3e (此值应为0)\n', Y_exp_actual - Y_prod);
    fprintf('   总消费 (C_total):            %.6f (占生产法GDP: %.2f%%)\n', C_total_agg, C_total_agg/Y_prod*100);
    fprintf('     - 常规消费 (有作用):       %.6f\n', C_utility_agg);
    fprintf('     - 重大冲击支出 (无作用):   %.6f\n', ShockExp_agg);
    fprintf('   总投资 (I):                  %.6f (占生产法GDP: %.2f%%)\n', I_agg_gross, I_agg_gross/Y_prod*100);
    fprintf('   政府支出 (G):                %.6f (占生产法GDP: %.2f%%)\n', G_agg, G_agg/Y_prod*100);
    
    % --- B. 政府与养老金核算 ---
    fprintf('\n--- B. 政府与养老金体系核算 ---\n');
    Total_Gov_Revenue = Tax_agg + Bequest_tax_agg;
    fprintf('   政府总收入 (常规税+遗赠税):  %.6f\n', Total_Gov_Revenue);
    fprintf('   政府支出 (G):                %.6f\n', G_agg);
    fprintf('   政府预算平衡 (收入-支出):      %.3e\n', Total_Gov_Revenue - G_agg);
    fprintf('   ------------------------------------\n');
    fprintf('   养老金体系收入 (缴款):       %.6f\n', PensionIn_agg);
    fprintf('   养老金体系支出 (发放):       %.6f\n', PensionOut_agg);
    fprintf('   养老金体系平衡 (收入-支出):  %.3e\n', PensionIn_agg - PensionOut_agg);

    % --- C. 宏观一致性最终检验 ---
    fprintf('\n--- C. 宏观一致性最终检验 ---\n');
    NetInvestment = I_agg_gross - Depreciation; % 稳态下净投资为0
    NetSaving_National = Y_prod - C_total_agg - G_agg - Depreciation;
    fprintf('   经济净储蓄 (S_net = Y-C-G-D): %.6f\n', NetSaving_National);
    fprintf('   经济净投资 (I_net = I-D):     %.6f\n', NetInvestment);
    fprintf('   ------------------------------------\n');
    fprintf('   储蓄-投资缺口 (S_net - I_net): %.3e (此值在稳态下应为0)\n', NetSaving_National - NetInvestment);
    
    % --- D. [最终结果] 资本产出比 ---
    fprintf('\n--- D. 关键宏观比率 ---\n');
    fprintf('   资本产出比 (K/Y):            %.4f\n', K_supply / Y_prod);

    fprintf('\n========================================================================\n');
end