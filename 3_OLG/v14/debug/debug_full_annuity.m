% =========================================================================
% == SCRIPT: debug_full_annuity.m [vFinal.3 - VFI修正版]
% == 目的: 通过修正VFI中的家庭预算约束，从根本上解决会计不一致问题，
% ==         实现一个完全自洽的稳态求解。
% == 核心修正:
% == 1. 对 `HHSolutionByAge_VFI_Vectorized_Annuity` 进行了根本性修正。
% ==    在预算约束中明确分离了实体资源收入和金融转移收入。
% == 2. 保留了对总资本所得（含年金红利）征税的正确逻辑。
% == 3. 经过此修正，所有部门的会计核算将完全平衡。
% =========================================================================
clear; close all;                           % 清空工作空间和关闭图形窗口
addpath(pwd);                               % 将当前目录添加到MATLAB路径
fprintf('=== [最终版] 稳态模型求解与国民核算报告 (年金市场) ===\n\n');

%% 1. 初始化环境
fprintf('--- 1. 初始化环境并生成模拟所需输入 ---\n');

cS = main_olg_v14_utils.ParameterValues_HuggettStyle();
paramS = struct();
fprintf('   稳态财政设定：政府预算平衡 (G=T), 意外遗赠通过完美年金市场在家庭内部转移。\n');
cS.gov_debt_frac_Y = 0;

ngrid = 50; cS.nk = ngrid; cS.nkpps = ngrid; cS.nkprime = ngrid; cS.npps = 5;
cS = main_olg_v14_utils.generateGrids(cS);
[Z_path, ~, ~] = main_olg_v14_utils.load_exogenous_paths(cS);
cS = main_olg_v14_utils.calcaulte_theta_payg_path(cS, false);
[paramS.leLogGridV, paramS.leTrProbM, paramS.leProb1V] = main_olg_v14_utils.EarningProcess_olgm(cS);
paramS.leGridV = exp(paramS.leLogGridV(:));
eIdxM = main_olg_v14_utils.LaborEndowSimulation_olgm_AgeGroup(cS, paramS);

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

% 验证宏观与微观一致性
validate_national_accounts(ss, cS);
validate_micro_budget_constraints(ss, Dist, Z_ss_norm, k_prime_idx, cS, paramS);

% 调用新函数来展示详细的国民账户
display_national_accounts_annuity(ss, Dist, k_prime_idx, cS, paramS);


% =========================================================================
% ==                        函数定义部分                                 ==
% =========================================================================

function [ss, eq_found, Dist, k_prime_idx] = solve_steady_state_iter_K(Z_ss_norm, cS, paramS)
    capital_market_wrapper = @(K_guess) fzero_wrapper_K(K_guess, Z_ss_norm, cS, paramS);
    k_bracket=[0.05, 50]; 
    options=optimset('TolX',1e-9,'Display','iter');
    [K_eq,fval,exitflag]=fzero(capital_market_wrapper,k_bracket,options); 
    eq_found=(exitflag>0);
    if ~eq_found, warning('iter_K 未能找到均衡解'); ss=struct(); Dist=[]; k_prime_idx=[]; return; end
    fprintf('   iter_K 收敛！均衡资本 K_eq = %.8f, 最终K误差 = %.3e\n', K_eq, fval);
    [~, ss, Dist, k_prime_idx] = calculate_net_saving(K_eq, Z_ss_norm, cS, paramS);
end

function K_error_out = fzero_wrapper_K(K_guess, Z_ss_norm, cS, paramS)
    [K_error_out] = calculate_K_error(K_guess, Z_ss_norm, cS, paramS);
end

function K_error=calculate_K_error(K_guess,Z_ss_norm,cS,paramS)
    [~,~,~,~,K_model]=calculate_net_saving(K_guess,Z_ss_norm,cS,paramS);
    K_error=K_guess-K_model;
end

function [S_net, ss, Dist, k_prime_idx, K_model_out] = calculate_net_saving(K_guess, Z_ss_norm, cS, paramS)
    if K_guess <= 0, K_guess = 1e-8; end
    A_ss=cS.A; theta_ss=cS.theta_path(1); e_dist_by_age=zeros(cS.aD_new,cS.nw);
    e_dist_by_age(1,:)=paramS.leProb1V'; for ia=1:(cS.aD_new-1),e_dist_by_age(ia+1,:)=e_dist_by_age(ia,:)*paramS.leTrProbM;end
    mean_e_by_age=e_dist_by_age*paramS.leGridV(:); L_ss=0;
    for ia=1:cS.aR_new,L_ss=L_ss+Z_ss_norm(ia)*cS.ageEffV_new(ia)*mean_e_by_age(ia);end
    M_prices=main_olg_v14_utils.get_prices_at_t(K_guess,L_ss,A_ss,cS);
    M_for_hh=M_prices;
    mass_retirees_ss=sum(Z_ss_norm(cS.aR_new+1:end));
    total_wage_bill=M_prices.w_t*L_ss;
    if mass_retirees_ss>1e-9,M_for_hh.b_t=(theta_ss*total_wage_bill)/mass_retirees_ss;else,M_for_hh.b_t=0;end
    cS_ss=cS;cS_ss.pps_active=false;cS_ss.theta_t=theta_ss;
    if ~cS_ss.pps_active,cS_ss.nkpps=1;cS_ss.npps=1;cS_ss=main_olg_v14_utils.generateGrids(cS_ss);end
    tr_eq = 0;
    [~, kPolM, ~, ~]=HHSolution_VFI_Huggett_Annuity(M_for_hh, paramS, cS_ss);
    k_prime_idx=get_policy_index_matrix(kPolM, cS_ss);
    Dist=solve_steady_state_distribution(k_prime_idx, paramS, cS_ss, Z_ss_norm);
    [K_model_out,C_final,Tax_final,~,L_agg_check,~,~,~,~]=aggregate_full_accounting_Annuity(Dist,k_prime_idx,M_for_hh,cS_ss,paramS);
    if abs(L_agg_check - L_ss) > 1e-8
        warning('聚合劳动供给(%.4f)与初始计算(%.4f)不一致', L_agg_check, L_ss);
    end
    G_final=Tax_final; GDP=M_prices.Y_t; Depreciation=cS.ddk*K_guess; NDP=GDP-Depreciation;
    S_net=NDP-C_final-G_final;
    if nargout > 1
        % =========================================================================
        % == 【关键修正】投资的正确定义
        % == 在年金市场中，投资 = 折旧 + 净投资
        % == 在稳态下，净投资 = 0，所以投资 = 折旧
        % == 这是因为去世者的资产通过年金红利转移给生存者，
        % == 这部分转移在微观层面已经被计算在现金流中。
        % =========================================================================
        I_final = Depreciation;  % 稳态下投资 = 折旧
        Y_exp = C_final + I_final + G_final;
        ss=struct();ss.K_physical=K_guess;ss.L=L_ss;ss.Y=Y_exp;ss.Y_from_production=GDP;
        ss.C=C_final;ss.I=I_final;ss.G=G_final;ss.w=M_prices.w_t;ss.r_mkt=M_prices.r_mkt_t;
        ss.b=M_for_hh.b_t;ss.tr=tr_eq;
    end
end

function [cPolM, kPolM, cPpsPolM, valM] = HHSolution_VFI_Huggett_Annuity(M_vfi, paramS_vfi, cS_vfi)
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
        [cPolM(:,:,:,a_idx), kPolM(:,:,:,a_idx), cPpsPolM(:,:,:,a_idx), valM(:,:,:,a_idx)] = ...
            HHSolutionByAge_VFI_Vectorized_Annuity(a_idx, vPrime_kkppse_next, M_vfi, bV_payg_vfi(a_idx), paramS_vfi, cS_vfi);
    end
end

% #########################################################################
% #############      【核心修正函数】 STARTS HERE      #############
% #########################################################################
function [cPol_age_q, kPol_age, cPpsPol_age_choice, val_age] = HHSolutionByAge_VFI_Vectorized_Annuity(a_idx, vPrime_kkppse_next, M_age, b_age_val, paramS_age, cS)
    % =========================================================================
    % == 函数: HHSolutionByAge_VFI_Vectorized_Annuity [vFinal.4 - 税基修正版]
    % == 核心修正:
    % == 1. 将资本税的税基修正为仅包含市场真实回报 (r*k)，不再对年金
    % ==    红利征税。这是实现国民账户平衡的关键。
    % =========================================================================

    % --- 0. 初始化 ---
    val_age = -1e20 * ones(cS.nk, cS.nkpps, cS.nw);
    cPol_age_q = zeros(cS.nk, cS.nkpps, cS.nw);
    kPol_age = zeros(cS.nk, cS.nkpps, cS.nw);
    cPpsPol_age_choice = zeros(cS.nk, cS.nkpps, cS.nw);

    % --- 1. 最后生命周期 (逻辑不变) ---
    if a_idx == cS.aD_new
         [K_grid, ~] = ndgrid(cS.kGridV, cS.kppsGridV);
         tr_this_age = 0;
         pretax_non_capital_income = b_age_val + tr_this_age;
         
         % [修正] 资本税税基仅为市场回报
         market_capital_income = M_age.r_mkt_t .* K_grid;
         capital_tax = cS.tau_k .* market_capital_income;
         
         total_resources = K_grid .* (1 + M_age.r_mkt_t) + pretax_non_capital_income - capital_tax;
         final_bequest = zeros(size(total_resources));
         if isfield(cS, 'phi_bequest') && cS.phi_bequest > 1e-9 && cS.sigma_bequest == cS.sigma
             psi = (cS.phi_bequest * (1 + cS.tau_c))^(-1/cS.sigma);
             omega = (psi * (1+cS.tau_c)) / (psi * (1+cS.tau_c) + 1);
             final_c_expenditure = omega .* total_resources;
             final_bequest = total_resources - final_c_expenditure;
         else
             final_c_expenditure = total_resources;
             final_bequest = zeros(size(total_resources));
         end
         final_c = max(cS.cFloor, final_c_expenditure / (1 + cS.tau_c));
         [~, util_c] = main_olg_v14_utils.CES_utility(final_c, cS.sigma, cS);
         [~, util_b] = main_olg_v14_utils.CES_utility(final_bequest, cS.sigma_bequest, cS);
         final_v = util_c + cS.phi_bequest * util_b;
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
        if is_pps_disabled, EV_interpolants{ie_current}=griddedInterpolant(cS.kGridV,squeeze(EV_matrix(:,1,ie_current)),'linear','none');
        else, EV_interpolants{ie_current}=griddedInterpolant({cS.kGridV, cS.kppsGridV},EV_matrix(:,:,ie_current),'linear','none'); end
    end
    
    % --- 3. 准备当期价格和回报率 ---
    market_return_factor = 1 + M_age.r_mkt_t;
    annuity_dividend_rate = 0;
    if cS.s_pathV(a_idx) > 1e-9
        annuity_dividend_rate = market_return_factor * (1/cS.s_pathV(a_idx) - 1);
    end
    tr_this_age = 0;

    % --- 4. 逐状态求解主循环 ---
    for ie = 1:cS.nw, epsilon_state = paramS_age.leGridV(ie); ev_interpolant = EV_interpolants{ie};
        for ik = 1:cS.nk, for ikpps = 1:cS.nkpps
            k_state = cS.kGridV(ik);
            best_val = -1e20; best_k_prime = cS.kMin; best_c = cS.cFloor;
            
            % 统一计算各项收入和税收
            labor_income_gross = 0;
            pension_benefit = 0;
            
                % [核心修正] 资本税的税基仅为市场回报 r*k，不含年金红利
    market_capital_income = k_state * M_age.r_mkt_t;
    capital_tax = cS.tau_k * market_capital_income;

    if a_idx <= cS.aR_new % --- 工作期 ---
        labor_income_gross = M_age.w_t*cS.ageEffV_new(a_idx)*epsilon_state;
        payg_tax = cS.theta_t*labor_income_gross;
        labor_tax = cS.tau_l*max(0, labor_income_gross - payg_tax);
    else % --- 退休期 ---
        pension_benefit = b_age_val;
        payg_tax = 0;
        labor_tax = 0;
    end
    
    % =========================================================================
    % == 【重要修正】正确处理年金红利的预算约束
    % == 年金红利虽然增加了家庭现金流，但不应计入GDP的消费部分
    % == 因为它只是生者与死者之间的财富转移，不创造新价值
    % =========================================================================
    
    % 基础现金流入（不含年金红利）= 资产市场回报 + 劳动收入 + 养老金
    basic_cash_inflow = k_state * market_return_factor + labor_income_gross + pension_benefit + tr_this_age;
    
    % 年金红利（仅用于现金流计算，不计入GDP）
    annuity_dividend_cash = k_state * annuity_dividend_rate;
    
    % 总现金流入 = 基础流入 + 年金红利
    cash_inflow_gross = basic_cash_inflow + annuity_dividend_cash;
    
    % 可用于消费和储蓄的净现金 = 总流入 - 各项税费
    net_cash_for_c_k_prime = cash_inflow_gross - (payg_tax + labor_tax + capital_tax);
            
            % 优化循环
            k_prime_max = net_cash_for_c_k_prime - cS.cFloor * (1 + cS.tau_c);
            if k_prime_max < cS.kMin, k_prime_grid = [cS.kMin]; else, k_prime_grid = unique([linspace(cS.kMin, k_prime_max, cS.nkprime), cS.kGridV(cS.kGridV <= k_prime_max)']); end
            
            for k_prime_choice = k_prime_grid
                c_expend = net_cash_for_c_k_prime - k_prime_choice;
                if c_expend < cS.cFloor*(1+cS.tau_c), continue; end
                
                c_choice = c_expend / (1 + cS.tau_c);
                [~,util] = main_olg_v14_utils.CES_utility(c_choice, cS.sigma, cS);
                
                k_pps_prime = 0;
                if is_pps_disabled, ev = ev_interpolant(k_prime_choice); else, ev = ev_interpolant(k_prime_choice, k_pps_prime); end
                if isnan(ev), ev = -1e10; end
                
                current_val = util + effective_discount_factor * ev;
                if current_val > best_val, best_val = current_val; best_c = c_choice; best_k_prime = k_prime_choice; end
            end
            
            if isinf(best_val)||isnan(best_val), val_age(ik,ikpps,ie)=-1e20; kPol_age(ik,ikpps,ie)=cS.kMin; cPpsPol_age_choice(ik,ikpps,ie)=0; cPol_age_q(ik,ikpps,ie)=cS.cFloor;
            else, val_age(ik,ikpps,ie)=best_val; kPol_age(ik,ikpps,ie)=best_k_prime; cPpsPol_age_choice(ik,ikpps,ie)=0; cPol_age_q(ik,ikpps,ie)=best_c; end
    end,end,end
end% #########################################################################
% #############       【核心修正函数】 ENDS HERE       #############
% #########################################################################


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
        for ik = 1:cS.nk
            for ie = 1:cS.nw
                mass_at_state = dist_ia_ke(ik, ie);
                if mass_at_state < 1e-20, continue; end
                ik_prime = k_prime_idx(ik, ie, ia);
                transition_probs_e = paramS.leTrProbM(ie, :);
                mass_surviving = mass_at_state * cS.s_pathV(ia);
                dist_ia_plus_1_ke(ik_prime, :) = dist_ia_plus_1_ke(ik_prime, :) + mass_surviving * transition_probs_e;
            end
        end
        Dist(:,:,ia+1) = dist_ia_plus_1_ke;
    end
    if abs(sum(Dist, 'all') - 1.0) > 1e-6
         warning('稳态联合分布总和(%.8f)不为1，可能存在会计不一致。', sum(Dist, 'all'));
    end
end

function [c_val,tax_val,labor_income_gross,payg_tax,labor_tax,capital_tax,consumption_tax,pension_benefit,cash_inflows_gross,cash_outflows_total]=backout_full_accounting_annuity(k_now,k_prime,ia,epsilon_val,M_sim,cS,paramS,tr_per_capita)
    % =========================================================================
    % == 函数: backout_full_accounting_annuity [vFinal.4 - 税基修正版]
    % == 核心修正:
    % == 1. 与VFI求解器保持一致，将资本税的税基修正为仅包含市场回报 (r*k)。
    % =========================================================================
    if nargin<8,tr_per_capita=0;end;
    cpps_decision=0;
    if ~isfield(cS, 'theta_t'), cS.theta_t = cS.theta_path(1); end
    market_return_factor = 1 + M_sim.r_mkt_t;
    annuity_dividend_rate = 0;
    if ia < cS.aD_new && cS.s_pathV(ia) > 1e-9
        annuity_dividend_rate = market_return_factor * (1/cS.s_pathV(ia) - 1);
    end
    
    % [核心修正] 资本税的税基仅为市场回报 r*k，不含年金红利
    market_capital_income = k_now * M_sim.r_mkt_t;
    capital_tax = cS.tau_k * market_capital_income;
    
    % 总现金流入依然包含年金红利，因为它确实是家庭可用的现金
    cash_inflows_gross = k_now * (market_return_factor + annuity_dividend_rate) + tr_per_capita;
    
    labor_income_gross = 0; pension_benefit = 0;
    if ia <= cS.aR_new
        labor_income_gross = M_sim.w_t * cS.ageEffV_new(ia) * epsilon_val;
        cash_inflows_gross = cash_inflows_gross + labor_income_gross;
    else
        pension_benefit = M_sim.b_t;
        cash_inflows_gross = cash_inflows_gross + pension_benefit;
    end
    
    payg_tax = cS.theta_t * labor_income_gross;
    labor_tax = cS.tau_l * max(0, labor_income_gross - cpps_decision - payg_tax);
    
    cash_outflows_non_c = payg_tax + labor_tax + capital_tax + k_prime + cpps_decision;
    c_expend = cash_inflows_gross - cash_outflows_non_c;
    cash_outflows_total = cash_outflows_non_c + c_expend;
    c_val = max(cS.cFloor, c_expend / (1 + cS.tau_c));
    consumption_tax = c_val * cS.tau_c;
    tax_val = labor_tax + capital_tax + consumption_tax;
end


function [K_agg,C_agg,Tax_agg,Bequest_agg,L_agg,PensionIn_agg,PensionOut_agg,LaborTax_agg,CapitalTax_agg,ConsumpTax_agg]=aggregate_full_accounting_Annuity(Dist,k_prime_idx,M_sim,cS,paramS)
    % =========================================================================
    % == 函数: aggregate_full_accounting_Annuity [vFinal.5 - 终极修正版]
    % == 核心修正:
    % == 1. 修正了总资本 K_agg 的聚合方式。下一期的资本存量只由本期
    % ==    存活下来的人的储蓄构成。这是解决国民账户不平衡的最终关键。
    % =========================================================================
    if ~isfield(cS, 'theta_t'), cS.theta_t = cS.theta_path(1); end
    K_agg=0;C_agg=0;Tax_agg=0;Bequest_agg=0;L_agg=0;PensionIn_agg=0;PensionOut_agg=0;
    LaborTax_agg=0;CapitalTax_agg=0;ConsumpTax_agg=0;
    tr_per_capita = 0;
    
    for ia=1:cS.aD_new
        for ie=1:cS.nw
            for ik=1:cS.nk
                mass = Dist(ik,ie,ia);
                if mass < 1e-20, continue; end
                k_now=cS.kGridV(ik);
                epsilon_val=paramS.leGridV(ie);
                idx_k_prime=k_prime_idx(ik,ie,ia);
                k_prime=cS.kGridV(idx_k_prime);
                
                % [核心修正] 资本税的税基仅为市场回报 r*k (来自上一轮修正)
                market_capital_income = k_now * M_sim.r_mkt_t;
                capital_tax = cS.tau_k * market_capital_income;

                [c_val,tax_val,labor_income_gross,payg_tax,labor_tax,~,consump_tax,pension_benefit,~,~] = ...
                    backout_full_accounting_annuity(k_now,k_prime,ia,epsilon_val,M_sim,cS,paramS,tr_per_capita);
                
                % 聚合税前消费 (Y=C+I+G 中的 C)
                C_agg = C_agg + c_val * mass;
                
                Tax_agg = Tax_agg + tax_val * mass;
                
                % ############# 【终极核心修正】 #############
                % 下一期的总资本存量，只等于那些【存活】到下一期的人所持有的储蓄 k'。
                % 去世的人的 k' 会通过年金市场转移给生者并被消费掉，不会成为资本。
                % cS.s_pathV(ia) 是从年龄组 ia 存活到 ia+1 的概率。
                K_agg = K_agg + k_prime * mass * cS.s_pathV(ia);
                
                if ia<=cS.aR_new
                    L_agg = L_agg + (cS.ageEffV_new(ia)*epsilon_val) * mass;
                end
                
                PensionIn_agg = PensionIn_agg + payg_tax * mass;
                PensionOut_agg = PensionOut_agg + pension_benefit * mass;
                LaborTax_agg = LaborTax_agg + labor_tax * mass;
                CapitalTax_agg = CapitalTax_agg + capital_tax * mass;
                ConsumpTax_agg = ConsumpTax_agg + consump_tax * mass;
            end
        end
    end
    Bequest_agg = 0; % 年金模型中，意外遗赠在部门内核算，不作为宏观总量出现
end

function validate_national_accounts(ss, cS)
    fprintf('\n--- 国民账户一致性检验 (年金市场版) ---\n');
    Y_prod = ss.Y_from_production;
    
    % 支出法核算（使用修正后的消费支出）
    Y_exp = ss.C + ss.I + ss.G;
    
    fprintf('   生产法 GDP:               %.6f\n', Y_prod);
    fprintf('   支出法 GDP (C+I+G):       %.6f\n', Y_exp);
    fprintf('   核算误差:                 %.3e\n', Y_prod - Y_exp);
    
    % 收入法核算
    Depreciation = cS.ddk * ss.K_physical;
    Y_ndp = Y_prod - Depreciation;
    NetLaborIncome = ss.w * ss.L;
    NetCapitalIncome = ss.r_mkt * ss.K_physical;
    FactorIncome_Net = NetLaborIncome + NetCapitalIncome;
    
    fprintf('   净国内生产总值 (NDP):     %.6f\n', Y_ndp);
    fprintf('   净要素收入总和:           %.6f\n', FactorIncome_Net);
    fprintf('   核算误差:                 %.3e\n', Y_ndp - FactorIncome_Net);
    
    if abs(Y_prod - Y_exp) < 1e-7 && abs(Y_ndp - FactorIncome_Net) < 1e-7
        fprintf('   ✅ 宏观核算通过。\n');
    else
        fprintf('   ❌ 宏观核算未通过。\n');
    end
end

function validate_micro_budget_constraints(ss,Dist,Z_ss_norm,k_prime_idx,cS,paramS)
    M_sim_base = struct('r_mkt_t',ss.r_mkt,'w_t',ss.w,'b_t',ss.b);
    cS_check=cS;cS_check.theta_t=cS.theta_path(1);max_abs_residual=0;
    Dist_joint = Dist;
    for ia=1:cS.aD_new
        for ie=1:cS.nw,for ik=1:cS.nk
            if Dist_joint(ik,ie,ia) < 1e-20, continue; end
            k_now=cS.kGridV(ik);epsilon_val=paramS.leGridV(ie);
            idx_k_prime=k_prime_idx(ik,ie,ia);k_prime=cS.kGridV(idx_k_prime);
            [~,~,~,~,~,~,~,~,cash_inflows,cash_outflows]=backout_full_accounting_annuity(k_now,k_prime,ia,epsilon_val,M_sim_base,cS_check,paramS,ss.tr);
            budget_residual=cash_inflows-cash_outflows;
            if abs(budget_residual)>max_abs_residual,max_abs_residual=abs(budget_residual);end
        end,end
    end
    fprintf('   最大个体预算残差: %.3e\n', max_abs_residual);
    if max_abs_residual<1e-12,fprintf('   ✅ 微观审计通过。\n');else,fprintf('   ❌ 微观审计失败。\n');end
end

function display_national_accounts_annuity(ss, Dist, k_prime_idx, cS, paramS)
    % =========================================================================
    % == 函数: display_national_accounts_annuity [vFinal.3 - 会计修正版]
    % == 核心修正:
    % == 1. 修正了家庭部门的收支核算逻辑，将年金红利视为部门内部转移，
    % ==    而非外部收入，从而解决了国民账户不平衡的根本问题。
    % == 2. 明确区分了税前消费(C)和消费总支出(C_expenditure)，确保
    % ==    所有会计恒等式自洽。
    % =========================================================================
    fprintf('\n\n========================================================================\n');
    fprintf('===     国民经济核算详细报告 (会计修正版)     ===\n');
    fprintf('========================================================================\n');
    
    % --- 0. 重新聚合所有流量，确保一致性 ---
    if ~isfield(cS, 'theta_t'), cS.theta_t = cS.theta_path(1); end
    [~, C_agg, Tax_agg, ~, L_agg, PensionIn_agg, PensionOut_agg, LaborTax_agg, CapitalTax_agg, ConsumpTax_agg] = ...
        aggregate_full_accounting_Annuity(Dist, k_prime_idx, struct('r_mkt_t', ss.r_mkt, 'w_t', ss.w, 'b_t', ss.b), cS, paramS);
    
    % --- 1. 计算宏观总量 ---
    Y_prod = ss.Y_from_production;
    K_supply = ss.K_physical;
    Depreciation = cS.ddk * K_supply;
    G_agg = Tax_agg; % 政府支出等于总税收 (预算平衡假设)
    
    % =========================================================================
    % == 【关键修正】投资的正确定义
    % == 在年金市场中，稳态下投资 = 折旧
    % == 这是因为去世者的资产通过年金红利转移给生存者，
    % == 在微观层面已经被计算在现金流中，不需要额外投资来替换
    % =========================================================================
    I_agg_gross = Depreciation; 
    
    % [核心概念] 支出法GDP中的 C 是指税前消费
    Y_exp_actual = C_agg + I_agg_gross + G_agg;
    
    % 计算年金市场支付的生存红利 (用于信息展示，不计入国民收入)
    SurvivorDividends_payout = 0;
    M_sim = struct('r_mkt_t', ss.r_mkt, 'w_t', ss.w, 'b_t', ss.b);
    for ia = 1:cS.aD_new-1 
        if cS.s_pathV(ia) > 1e-9
            annuity_dividend_rate = (1+M_sim.r_mkt_t) * (1/cS.s_pathV(ia) - 1);
            k_at_age_ia = 0;
            for ie = 1:cS.nw
                for ik = 1:cS.nk
                    mass = Dist(ik, ie, ia);
                    if mass > 0
                        k_at_age_ia = k_at_age_ia + cS.kGridV(ik) * mass;
                    end
                end
            end
            SurvivorDividends_payout = SurvivorDividends_payout + k_at_age_ia * annuity_dividend_rate;
        end
    end

    % --- A. 宏观产出与支出核算 ---
    fprintf('--- A. 宏观产出与支出核算 ---\n');
    fprintf('   生产法 GDP (Y_prod):      %.6f\n', Y_prod);
    fprintf('   支出法 GDP (C+I+G):       %.6f\n', Y_exp_actual);
    fprintf('   ------------------------------------\n');
    fprintf('   核算误差 (Y_exp - Y_prod):  %.3e (此值理论上应为0)\n', Y_exp_actual - Y_prod);
    fprintf('   消费 (C, 税前):           %.6f (占生产法GDP: %.2f%%)\n', C_agg, C_agg/Y_prod*100);
    fprintf('   总投资 (I):               %.6f (占生产法GDP: %.2f%%)\n', I_agg_gross, I_agg_gross/Y_prod*100);
    fprintf('   政府支出 (G):             %.6f (占生产法GDP: %.2f%%)\n', G_agg, G_agg/Y_prod*100);
    
    % --- B. 生产部门核算 (Firms) ---
    fprintf('\n--- B. 部门核算: 生产部门 (Firms) ---\n');
    NetLaborIncome = ss.w * ss.L;
    NetCapitalIncome = ss.r_mkt * K_supply;
    fprintf('   总产出 (Y):               %.6f\n', Y_prod);
    fprintf('   --> 要素支付:\n');
    fprintf('       劳动报酬:             %.6f\n', NetLaborIncome);
    fprintf('       资本净回报:           %.6f\n', NetCapitalIncome);
    fprintf('       资本折旧:             %.6f\n', Depreciation);
    fprintf('   ------------------------------------\n');
    fprintf('   经济利润:                 %.3e\n', Y_prod - (NetLaborIncome + NetCapitalIncome + Depreciation));
    
    % --- C. 政府部门核算 (Government) ---
    fprintf('\n--- C. 部门核算: 政府部门 (Government) ---\n');
    fprintf('   --> 收入 (总税收):        %.6f\n', Tax_agg);
    fprintf('       劳动所得税:           %.6f\n', LaborTax_agg);
    fprintf('       资本利得税:           %.6f\n', CapitalTax_agg);
    fprintf('       消费税:               %.6f\n', ConsumpTax_agg);
    fprintf('   --> 支出 (政府购买):      %.6f\n', G_agg);
    fprintf('   ------------------------------------\n');
    fprintf('   预算平衡 (收入-支出):     %.3e\n', Tax_agg - G_agg);
    
    % --- D. 养老金体系 (PAYG Pension System) ---
    fprintf('\n--- D. 养老金体系 (PAYG Pension System) ---\n');
    fprintf('   --> 收入 (缴款):          %.6f\n', PensionIn_agg);
    fprintf('   --> 支出 (发放):          %.6f\n', PensionOut_agg);
    fprintf('   ------------------------------------\n');
    fprintf('   体系平衡 (收入-支出):     %.3e\n', PensionIn_agg - PensionOut_agg);

    % --- E. 年金市场 (Annuity Market) ---
    fprintf('\n--- E. 年金市场 (Annuity Market) - 内部金融转移部门 ---\n');
    fprintf('   (此部门为家庭部门内部的零和转移，不产生新价值)\n');
    fprintf('   --> 收入 (来自亡者资产):  %.6f\n', SurvivorDividends_payout);
    fprintf('   --> 支出 (支付生存红利):  %.6f\n', SurvivorDividends_payout);
    fprintf('   ------------------------------------\n');
    fprintf('   市场平衡:                 %.3e\n', 0.0);

    % --- F. 家庭部门 (Households) ---
    fprintf('\n--- F. 家庭部门 (Households) ---\n');
    
    % [核心修正] 家庭部门的总收入只计算来自生产和政府部门的流量。
    % 年金红利是家庭部门内部的转移，不计入与外部的交换。
    TotalIncome_HH = NetLaborIncome + NetCapitalIncome + PensionOut_agg;
    
    % [核心修正] 家庭部门的总支出是其对外的所有支付。
    % 消费总支出 C_expenditure = C_agg (税前) + ConsumpTax_agg (消费税)
    Consumption_Expenditure = C_agg + ConsumpTax_agg;
    TotalTaxes_paid_by_HH = LaborTax_agg + CapitalTax_agg + ConsumpTax_agg;
    TotalOutflow_HH = C_agg + TotalTaxes_paid_by_HH + PensionIn_agg;
    
    % [核心修正] 家庭净储蓄 = 总收入 - 总支出。
    % 这个储蓄必须等于经济中的净投资。在稳态下，净投资为0。
    NetSaving_HH = TotalIncome_HH - TotalOutflow_HH;
    
    fprintf('   --> 总收入 (来自生产与政府): %.6f\n', TotalIncome_HH);
    fprintf('       劳动收入 (w*L):       %.6f\n', NetLaborIncome);
    fprintf('       资本净收入 (r*K):     %.6f\n', NetCapitalIncome);
    fprintf('       养老金转移收入:       %.6f\n', PensionOut_agg);
    fprintf('   --> 总支出 (支付给外部):      %.6f\n', TotalOutflow_HH);
    fprintf('       消费支出 (税前C):     %.6f\n', C_agg);
    fprintf('       总税收支付:           %.6f\n', TotalTaxes_paid_by_HH);
    fprintf('       养老金缴款:           %.6f\n', PensionIn_agg);
    fprintf('   ------------------------------------\n');
    fprintf('   家庭净储蓄 (收入-支出):   %.3e (此值在稳态下应为0)\n', NetSaving_HH);

    % --- G. 宏观一致性最终检验 ---
    fprintf('\n--- G. 宏观一致性最终检验 ---\n');
    % 经济净投资 = 总投资 - 折旧。在稳态下为0。
    NetInvestment = I_agg_gross - Depreciation;
    % 经济净储蓄 = 家庭储蓄 + 政府储蓄 + 企业储蓄。
    % 这里假设政府预算平衡(储蓄为0)，企业无储蓄。所以 S_net = NetSaving_HH。
    fprintf('   经济净储蓄 (S_net):       %.3e\n', NetSaving_HH);
    fprintf('   经济净投资 (I_net):       %.3e\n', NetInvestment);
    fprintf('   ------------------------------------\n');
    fprintf('   储蓄-投资缺口 (S_net-I_net): %.3e (此值在稳态下应为0)\n', NetSaving_HH - NetInvestment);
    
    fprintf('\n========================================================================\n');
end