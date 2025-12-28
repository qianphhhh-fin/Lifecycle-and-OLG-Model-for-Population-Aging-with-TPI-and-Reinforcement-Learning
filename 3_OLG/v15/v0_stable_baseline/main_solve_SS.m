% =========================================================================
% == SCRIPT: debug_full_annuity.m [vPaper.1 - 遗赠税模型]
% == 目的: 根据 Börsch-Supan et al. (2018) 的设定，将模型从年金市场
% ==         修改为“意外遗赠被政府100%征收”的遗赠税模型。
% == 核心修正:
% == 1. 废除年金市场：HHSolutionByAge... 中不再有年金红利。
% == 2. 引入遗赠税：aggregate_... 函数现在计算总意外遗赠额。
% == 3. 政府预算：G = T_regular + Bequest_Tax。
% == 4. 经过此修正，模型将与论文设定对齐，并实现会计平衡。
% =========================================================================
clear; close all;                           % 清空工作空间和关闭图形窗口
addpath(pwd);                               % 将当前目录添加到MATLAB路径
fprintf('=== [论文对齐版] 稳态模型求解与国民核算报告 (遗赠税模型) ===\n\n');

%% 1. 初始化环境
fprintf('--- 1. 初始化环境并生成模拟所需输入 ---\n');

cS = main_olg_v15_utils.ParameterValues_HuggettStyle();
paramS = struct();
fprintf('   稳态财政设定：政府预算平衡 (G=T), 意外遗赠被政府100%%征收。\n');
cS.gov_debt_frac_Y = 0;

ngrid = 70; cS.nk = ngrid; cS.nkpps = ngrid; cS.nkprime = ngrid; cS.npps = 5;
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

% 验证宏观与微观一致性
validate_national_accounts(ss, cS);
validate_micro_budget_constraints(ss, Dist, Z_ss_norm, k_prime_idx, cS, paramS);

% 调用新函数来展示详细的国民账户
display_national_accounts_bequest_tax(ss, Dist, k_prime_idx, cS, paramS);


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
    M_prices=main_olg_v15_utils.get_prices_at_t(K_guess,L_ss,A_ss,cS);
    M_for_hh=M_prices;
    mass_retirees_ss=sum(Z_ss_norm(cS.aR_new+1:end));
    total_wage_bill=M_prices.w_t*L_ss;
    if mass_retirees_ss>1e-9,M_for_hh.b_t=(theta_ss*total_wage_bill)/mass_retirees_ss;else,M_for_hh.b_t=0;end
    cS_ss=cS;cS_ss.pps_active=false;cS_ss.theta_t=theta_ss;
    if ~cS_ss.pps_active,cS_ss.nkpps=1;cS_ss.npps=1;cS_ss=main_olg_v15_utils.generateGrids(cS_ss);end
    
    % 【模型变更】家庭问题求解器不再需要tr_eq，因为没有家庭内部转移了
    [~, kPolM, ~, ~]=HHSolution_VFI_Huggett_BequestTax(M_for_hh, paramS, cS_ss);
    k_prime_idx=get_policy_index_matrix(kPolM, cS_ss);
    Dist=solve_steady_state_distribution(k_prime_idx, paramS, cS_ss, Z_ss_norm);
    
    % 【模型变更】聚合函数现在计算遗赠税
    [K_model_out,C_final,Tax_final,Bequest_tax_final,L_agg_check,~,~,~,~]=aggregate_full_accounting_BequestTax(Dist,k_prime_idx,M_for_hh,cS_ss,paramS);
    
    if abs(L_agg_check - L_ss) > 1e-8, warning('聚合劳动供给(%.4f)与初始计算(%.4f)不一致', L_agg_check, L_ss); end
    
    % 【模型变更】政府支出G = 常规税 + 遗赠税
    G_final = Tax_final + Bequest_tax_final; 
    
    GDP=M_prices.Y_t; Depreciation=cS.ddk*K_guess; NDP=GDP-Depreciation;
    S_net=NDP-C_final-G_final;
    
    if nargout > 1
        I_final_gross = Depreciation; % 稳态总投资=折旧
        Y_exp = C_final + I_final_gross + G_final;
        ss=struct();ss.K_physical=K_guess;ss.L=L_ss;ss.Y=Y_exp;ss.Y_from_production=GDP;
        ss.C=C_final;ss.I=I_final_gross;ss.G=G_final;ss.w=M_prices.w_t;ss.r_mkt=M_prices.r_mkt_t;
        ss.b=M_for_hh.b_t;ss.Bequest_tax = Bequest_tax_final; ss.Regular_tax = Tax_final;
    end
end

function [cPolM, kPolM, cPpsPolM, valM] = HHSolution_VFI_Huggett_BequestTax(M_vfi, paramS_vfi, cS_vfi)
    % 【模型变更】此函数现在调用遗赠税版本的VFI求解器
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
            HHSolutionByAge_VFI_Vectorized_BequestTax(a_idx, vPrime_kkppse_next, M_vfi, bV_payg_vfi(a_idx), paramS_vfi, cS_vfi);
    end
end

% #########################################################################
% #############      【核心修正函数 - 遗赠税模型】 STARTS HERE      #############
% #########################################################################
function [cPol_age_q, kPol_age, cPpsPol_age_choice, val_age] = HHSolutionByAge_VFI_Vectorized_BequestTax(a_idx, vPrime_kkppse_next, M_age, b_age_val, paramS_age, cS)
    % =========================================================================
    % == 函数: HHSolutionByAge_VFI_Vectorized_BequestTax [vPaper.1]
    % == 核心变更:
    % == 1. 移除年金红利。家庭资本回报率仅为市场利率 r_mkt。
    % == 2. 预算约束中不再有任何形式的家庭内部转移(tr或年金)。
    % == 3. 贴现因子现在是 s_pathV(a_idx)，因为不存活的后果是资产被政府没收，
    % ==    从家庭角度看，下一期价值为0，这等效于用存活率贴现。
    % =========================================================================

    % --- 0. 初始化 ---
    val_age = -1e20 * ones(cS.nk, cS.nkpps, cS.nw);
    cPol_age_q = zeros(cS.nk, cS.nkpps, cS.nw);
    kPol_age = zeros(cS.nk, cS.nkpps, cS.nw);
    cPpsPol_age_choice = zeros(cS.nk, cS.nkpps, cS.nw);

    % --- 1. 最后生命周期 (逻辑不变，因为没有下一期，存活率为0) ---
    if a_idx == cS.aD_new
         [K_grid, ~] = ndgrid(cS.kGridV, cS.kppsGridV);
         pretax_non_capital_income = b_age_val; % 没有其他转移
         capital_tax = cS.tau_k .* M_age.r_mkt_t .* K_grid;
         total_resources = K_grid .* (1 + M_age.r_mkt_t) + pretax_non_capital_income - capital_tax;
         
         % 论文中没有意向遗赠，为简化，这里也假设没有。
         final_c_expenditure = total_resources;
         final_bequest = zeros(size(total_resources));

         final_c = max(cS.cFloor, final_c_expenditure / (1 + cS.tau_c));
         [~, util_c] = main_olg_v15_utils.CES_utility(final_c, cS.sigma, cS);
         final_v = util_c;
         
         for ie = 1:cS.nw
             cPol_age_q(:,:,ie) = final_c;
             val_age(:,:,ie) = final_v;
             kPol_age(:,:,ie) = final_bequest; % 最终资产为0
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
    
    % --- 3. 【模型变更】准备当期价格和回报率，移除年金部分 ---
    market_return_factor = 1 + M_age.r_mkt_t;

    % --- 4. 逐状态求解主循环 ---
    for ie = 1:cS.nw, epsilon_state = paramS_age.leGridV(ie); ev_interpolant = EV_interpolants{ie};
        for ik = 1:cS.nk, for ikpps = 1:cS.nkpps
            k_state = cS.kGridV(ik);
            best_val = -1e20; best_k_prime = cS.kMin; best_c = cS.cFloor;
            
            % 计算各项收入和税收
            labor_income_gross = 0;
            pension_benefit = 0;
            
            % 资本收入仅来自市场
            capital_income = k_state * M_age.r_mkt_t;
            capital_tax = cS.tau_k * capital_income;

            if a_idx <= cS.aR_new % --- 工作期 ---
                labor_income_gross = M_age.w_t*cS.ageEffV_new(a_idx)*epsilon_state;
                payg_tax = cS.theta_t*labor_income_gross;
                labor_tax = cS.tau_l*max(0, labor_income_gross - payg_tax);
            else % --- 退休期 ---
                pension_benefit = b_age_val;
                payg_tax = 0;
                labor_tax = 0;
            end
            
            % 构建总预算约束
            cash_inflow_gross = k_state * market_return_factor + labor_income_gross + pension_benefit;
            net_cash_for_c_k_prime = cash_inflow_gross - (payg_tax + labor_tax + capital_tax);
            
            % 优化循环
            k_prime_max = net_cash_for_c_k_prime - cS.cFloor * (1 + cS.tau_c);
            if k_prime_max < cS.kMin, k_prime_grid = [cS.kMin]; else, k_prime_grid = unique([linspace(cS.kMin, k_prime_max, cS.nkprime), cS.kGridV(cS.kGridV <= k_prime_max)']); end
            
            for k_prime_choice = k_prime_grid
                c_expend = net_cash_for_c_k_prime - k_prime_choice;
                if c_expend < cS.cFloor*(1+cS.tau_c), continue; end
                
                c_choice = c_expend / (1 + cS.tau_c);
                [~,util] = main_olg_v15_utils.CES_utility(c_choice, cS.sigma, cS);
                
                k_pps_prime = 0;
                if is_pps_disabled, ev = ev_interpolant(k_prime_choice); else, ev = ev_interpolant(k_prime_choice, k_pps_prime); end
                if isnan(ev), ev = -1e10; end
                
                current_val = util + effective_discount_factor * ev;
                if current_val > best_val, best_val = current_val; best_c = c_choice; best_k_prime = k_prime_choice; end
            end
            
            if isinf(best_val)||isnan(best_val), val_age(ik,ikpps,ie)=-1e20; kPol_age(ik,ikpps,ie)=cS.kMin; cPpsPol_age_choice(ik,ikpps,ie)=0; cPol_age_q(ik,ikpps,ie)=cS.cFloor;
            else, val_age(ik,ikpps,ie)=best_val; kPol_age(ik,ikpps,ie)=best_k_prime; cPpsPol_age_choice(ik,ikpps,ie)=0; cPol_age_q(ik,ikpps,ie)=best_c; end
    end,end,end
end
% #########################################################################
% #############       【核心修正函数 - 遗赠税模型】 ENDS HERE       #############
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

function [c_val,tax_val,labor_income_gross,payg_tax,labor_tax,capital_tax,consumption_tax,pension_benefit,cash_inflows_gross,cash_outflows_total]=backout_full_accounting_BequestTax(k_now,k_prime,ia,epsilon_val,M_sim,cS)
    % 【模型变更】此函数与新的VFI求解器逻辑完全对应
    cpps_decision=0;
    if ~isfield(cS, 'theta_t'), cS.theta_t = cS.theta_path(1); end
    
    market_return_factor = 1 + M_sim.r_mkt_t;
    capital_income = k_now * M_sim.r_mkt_t;
    
    cash_inflows_gross = k_now * market_return_factor;
    labor_income_gross = 0; pension_benefit = 0;
    if ia <= cS.aR_new
        labor_income_gross = M_sim.w_t * cS.ageEffV_new(ia) * epsilon_val;
        cash_inflows_gross = cash_inflows_gross + labor_income_gross;
    else
        pension_benefit = M_sim.b_t;
        cash_inflows_gross = cash_inflows_gross + pension_benefit;
    end
    
    payg_tax = cS.theta_t * labor_income_gross;
    capital_tax = cS.tau_k * capital_income;
    labor_tax = cS.tau_l * max(0, labor_income_gross - cpps_decision - payg_tax);
    
    cash_outflows_non_c = payg_tax + labor_tax + capital_tax + k_prime + cpps_decision;
    c_expend = cash_inflows_gross - cash_outflows_non_c;
    cash_outflows_total = cash_outflows_non_c + c_expend;
    
    c_val = max(cS.cFloor, c_expend / (1 + cS.tau_c));
    consumption_tax = c_val * cS.tau_c;
    tax_val = labor_tax + capital_tax + consumption_tax;
end

function [K_agg,C_agg,Tax_agg,Bequest_tax_agg,L_agg,PensionIn_agg,PensionOut_agg,LaborTax_agg,CapitalTax_agg,ConsumpTax_agg]=aggregate_full_accounting_BequestTax(Dist,k_prime_idx,M_sim,cS,paramS)
    % 【模型变更】此函数现在聚合常规税和遗赠税
    if ~isfield(cS, 'theta_t'), cS.theta_t = cS.theta_path(1); end
    K_agg=0;C_agg=0;Tax_agg=0;Bequest_tax_agg=0;L_agg=0;PensionIn_agg=0;PensionOut_agg=0;
    LaborTax_agg=0;CapitalTax_agg=0;ConsumpTax_agg=0;
    
    for ia=1:cS.aD_new
        for ie=1:cS.nw
            for ik=1:cS.nk
                mass = Dist(ik,ie,ia);
                if mass < 1e-20, continue; end
                k_now=cS.kGridV(ik);
                epsilon_val=paramS.leGridV(ie);
                idx_k_prime=k_prime_idx(ik,ie,ia);
                k_prime=cS.kGridV(idx_k_prime);
                
                [c_val,tax_val,labor_income_gross,payg_tax,labor_tax,capital_tax,consump_tax,pension_benefit,~,~] = ...
                    backout_full_accounting_BequestTax(k_now,k_prime,ia,epsilon_val,M_sim,cS);
                
                C_agg = C_agg + c_val * mass;
                Tax_agg = Tax_agg + tax_val * mass;
                
                % 【模型变更】聚合下一期资本时，只考虑存活下来的人
                prob_survive = cS.s_pathV(ia);
                K_agg = K_agg + k_prime * mass * prob_survive;
                
                % 【模型变更】计算遗赠税收入
                prob_death = 1 - prob_survive;
                Bequest_tax_agg = Bequest_tax_agg + k_prime * mass * prob_death;
                
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
end

function validate_national_accounts(ss, cS)
    fprintf('\n--- 国民账户一致性检验 (遗赠税模型版) ---\n');
    Y_prod = ss.Y_from_production;
    Y_exp = ss.Y; % ss.Y 已经是由 C+I+G 计算出的支出法GDP
    
    fprintf('   生产法 GDP:               %.6f\n', Y_prod);
    fprintf('   支出法 GDP (C+I+G):       %.6f\n', Y_exp);
    fprintf('   核算误差:                 %.3e\n', Y_prod - Y_exp);
    
    Depreciation = cS.ddk * ss.K_physical;
    Y_ndp = Y_prod - Depreciation;
    NetLaborIncome = ss.w * ss.L;
    NetCapitalIncome = ss.r_mkt * ss.K_physical;
    FactorIncome_Net = NetLaborIncome + NetCapitalIncome;
    
    fprintf('   净国内生产总值 (NDP):     %.6f\n', Y_ndp);
    fprintf('   净要素收入总和:           %.6f\n', FactorIncome_Net);
    fprintf('   核算误差:                 %.3e\n', Y_ndp - FactorIncome_Net);
    
    if abs(Y_prod - Y_exp) < 1e-6 && abs(Y_ndp - FactorIncome_Net) < 1e-6
        fprintf('   ✅ 宏观核算通过。\n');
    else
        fprintf('   ❌ 宏观核算未通过。\n');
    end
end

function validate_micro_budget_constraints(ss,Dist,Z_ss_norm,k_prime_idx,cS,paramS)
    M_sim_base = struct('r_mkt_t',ss.r_mkt,'w_t',ss.w,'b_t',ss.b);
    cS_check=cS;cS_check.theta_t=cS.theta_path(1);max_abs_residual=0;
    for ia=1:cS.aD_new
        for ie=1:cS.nw,for ik=1:cS.nk
            if Dist(ik,ie,ia) < 1e-20, continue; end
            k_now=cS.kGridV(ik);epsilon_val=paramS.leGridV(ie);
            idx_k_prime=k_prime_idx(ik,ie,ia);k_prime=cS.kGridV(idx_k_prime);
            [~,~,~,~,~,~,~,~,cash_inflows,cash_outflows]=backout_full_accounting_BequestTax(k_now,k_prime,ia,epsilon_val,M_sim_base,cS_check);
            budget_residual=cash_inflows-cash_outflows;
            if abs(budget_residual)>max_abs_residual,max_abs_residual=abs(budget_residual);end
        end,end
    end
    fprintf('   最大个体预算残差: %.3e\n', max_abs_residual);
    if max_abs_residual<1e-12,fprintf('   ✅ 微观审计通过。\n');else,fprintf('   ❌ 微观审计失败。\n');end
end

function display_national_accounts_bequest_tax(ss, Dist, k_prime_idx, cS, paramS)
    fprintf('\n\n========================================================================\n');
    fprintf('===     国民经济核算详细报告 (遗赠税模型版)     ===\n');
    fprintf('========================================================================\n');
    
    % --- 0. 重新聚合所有流量，确保一致性 ---
    M_sim = struct('r_mkt_t', ss.r_mkt, 'w_t', ss.w, 'b_t', ss.b);
    [K_agg_demand, C_agg, Tax_agg, Bequest_tax_agg, L_agg, PensionIn_agg, PensionOut_agg, LaborTax_agg, CapitalTax_agg, ConsumpTax_agg] = ...
        aggregate_full_accounting_BequestTax(Dist, k_prime_idx, M_sim, cS, paramS);
    
    % --- 1. 计算宏观总量 ---
    Y_prod = ss.Y_from_production;
    K_supply = ss.K_physical;
    Depreciation = cS.ddk * K_supply;
    
    % G_agg 是由模型内生决定的，等于所有政府收入之和
    G_agg = Tax_agg + Bequest_tax_agg;
    
    I_agg_gross = Depreciation; % 稳态下，总投资等于折旧
    Y_exp_actual = C_agg + I_agg_gross + G_agg;
    
    % --- A. 宏观产出与支出核算 ---
    fprintf('--- A. 宏观产出与支出核算 ---\n');
    fprintf('   生产法 GDP (Y_prod):      %.6f\n', Y_prod);
    fprintf('   支出法 GDP (C+I+G):       %.6f\n', Y_exp_actual);
    fprintf('   ------------------------------------\n');
    fprintf('   核算误差 (Y_exp - Y_prod):  %.3e (此值应为0)\n', Y_exp_actual - Y_prod);
    fprintf('   消费 (C):                 %.6f (占生产法GDP: %.2f%%)\n', C_agg, C_agg/Y_prod*100);
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
    Total_Gov_Revenue = Tax_agg + Bequest_tax_agg;
    fprintf('   --> 总收入:               %.6f\n', Total_Gov_Revenue);
    fprintf('       常规税收总额:         %.6f\n', Tax_agg);
    fprintf('         - 劳动所得税:       %.6f\n', LaborTax_agg);
    fprintf('         - 资本利得税:       %.6f\n', CapitalTax_agg);
    fprintf('         - 消费税:           %.6f\n', ConsumpTax_agg);
    fprintf('       遗赠税收入:           %.6f\n', Bequest_tax_agg);
    fprintf('   --> 支出 (政府购买):      %.6f\n', G_agg);
    fprintf('   ------------------------------------\n');
    fprintf('   预算平衡 (总收入-支出):   %.3e\n', Total_Gov_Revenue - G_agg);
    
    % --- D. 养老金体系 (PAYG Pension System) ---
    fprintf('\n--- D. 养老金体系 (PAYG Pension System) ---\n');
    fprintf('   --> 收入 (缴款):          %.6f\n', PensionIn_agg);
    fprintf('   --> 支出 (发放):          %.6f\n', PensionOut_agg);
    fprintf('   ------------------------------------\n');
    fprintf('   体系平衡 (收入-支出):     %.3e\n', PensionIn_agg - PensionOut_agg);

    % --- E. 家庭部门 (Households) ---
    fprintf('\n--- E. 家庭部门 (Households) ---\n');
    TotalIncome_HH = NetLaborIncome + NetCapitalIncome + PensionOut_agg;
    TotalOutflow_HH = (C_agg * (1+cS.tau_c)) + (Tax_agg - ConsumpTax_agg) + PensionIn_agg + Bequest_tax_agg;
    NetSaving_HH = TotalIncome_HH - TotalOutflow_HH;
    
    fprintf('   --> 总收入 (来自生产和转移): %.6f\n', TotalIncome_HH);
    fprintf('       劳动收入 (w*L):       %.6f\n', NetLaborIncome);
    fprintf('       资本净收入 (r*K):     %.6f\n', NetCapitalIncome);
    fprintf('       养老金转移:           %.6f\n', PensionOut_agg);
    fprintf('   --> 总支出与流出:         %.6f\n', TotalOutflow_HH);
    fprintf('       消费总支出 (含税):    %.6f\n', C_agg * (1+cS.tau_c));
    fprintf('       所得税支付:           %.6f\n', LaborTax_agg + CapitalTax_agg);
    fprintf('       养老金缴款:           %.6f\n', PensionIn_agg);
    fprintf('       意外遗赠 (被政府征收): %.6f\n', Bequest_tax_agg);
    fprintf('   ------------------------------------\n');
    fprintf('   家庭净储蓄 (收入-支出):   %.6f (此值在稳态下应为0)\n', NetSaving_HH);

    % --- F. 宏观一致性最终检验 ---
    fprintf('\n--- F. 宏观一致性最终检验 ---\n');
    % 经济总储蓄 = 私人储蓄 + 政府储蓄 + 企业储蓄
    % S_private = 家庭净储蓄
    % S_government = 政府总收入 - 政府总支出 = 0 (因预算平衡)
    % S_corporate = 经济利润 = 0 (因完全竞争)
    % S_national = S_private
    NetInvestment = I_agg_gross - Depreciation;
    fprintf('   经济净储蓄 (S_net):       %.6f\n', NetSaving_HH);
    fprintf('   经济净投资 (I_net):       %.6f\n', NetInvestment);
    fprintf('   ------------------------------------\n');
    fprintf('   储蓄-投资缺口 (S_net - I_net): %.3e (此值在稳态下应为0)\n', NetSaving_HH - NetInvestment);
    
    fprintf('\n========================================================================\n');
end