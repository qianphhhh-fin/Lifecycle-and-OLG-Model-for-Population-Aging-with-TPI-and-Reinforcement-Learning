% --- household.m ---
classdef household
    methods (Static)


        function [polS, valS] = VFI_solver(M_vfi, paramS_vfi, cS_vfi)
            % =========================================================================
            % == 函数: VFI_solver (版本 v2.2 - 养老金(b_hat)兼容版)
            % == 目的: 使稳态VFI求解器能够接收和处理养老金收入。
            % =========================================================================

            if isfield(cS_vfi, 'pps_active') && cS_vfi.pps_active
                valS = -Inf(cS_vfi.nk, cS_vfi.nkpps, cS_vfi.nw_expanded, cS_vfi.aD_new);
            else
                valS = -Inf(cS_vfi.nk, 1, cS_vfi.nw_expanded, cS_vfi.aD_new);
            end

            polS_cell = cell(cS_vfi.aD_new, 1);

            % [核心修正] 从 M_vfi 中获取养老金福利 b_hat，并构建年龄相关的养老金向量
            b_payg_hat_age_vec = zeros(1, cS_vfi.aD_new);
            if isfield(M_vfi, 'b_hat_t') && M_vfi.b_hat_t > 0
                % 仅为退休后年龄组分配养老金
                b_payg_hat_age_vec((cS_vfi.aR_new+1):end) = M_vfi.b_hat_t;
            end

            for a_idx = cS_vfi.aD_new : -1 : 1
                vPrime_kkppse_next = [];
                if a_idx < cS_vfi.aD_new
                    vPrime_kkppse_next = valS(:,:,:,a_idx+1);
                end

                % [核心修正] 将特定年龄的养老金 b_age_val 传递给引擎
                b_age_val = b_payg_hat_age_vec(a_idx);

                if isfield(cS_vfi, 'pps_active') && cS_vfi.pps_active
                    [val_age, pol_age] = household.VFI_by_age_PPS(a_idx, vPrime_kkppse_next, M_vfi, b_age_val, paramS_vfi, cS_vfi);
                else
                    [val_age, pol_age] = household.VFI_by_age(a_idx, vPrime_kkppse_next, M_vfi, b_age_val, paramS_vfi, cS_vfi);
                end

                valS(:,:,:,a_idx) = val_age;
                polS_cell{a_idx} = pol_age;
            end
            polS = [polS_cell{:}];
        end


       function [val_age, pol_age] = VFI_by_age_PPS(a_idx, vPrime_kkppse_next, M_age, b_age_val, paramS_age, cS)
            % =========================================================================
            % == 函数: VFI_by_age_PPS (版本 v1.6 - BGP逻辑完全重构版)
            % ==
            % == 核心修正:
            % ==   - 彻底重写了函数核心，确保所有计算都在`hat`(去技术趋势)空间内进行。
            % ==   - [!!!] PPS提取额(pps_withdrawal_hat)现在被正确计算为
            % ==     当期`hat`资产赚取利息后的一个比例。
            % ==   - [!!!] PPS资产的积累方程(kpps_prime_hat)现在遵循正确的BGP动态：
            % ==     kpps_prime_hat = (期末hat值) / (1+g_A)。
            % ==   - [!!!] 家庭预算约束(cash_on_hand_hat_vec)现在是单位完全一致的。
            % ==   - 此修正旨在从根本上解决PPS模块在 n>0 或 g_A>0 时的BGP不一致问题。
            % =========================================================================
            nk_search_grid = 150;
            val_age = -1e20 * ones(cS.nk, cS.nkpps, cS.nw_expanded);
            pol_age = struct(...
                'c', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
                'k_prime', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
                'kpps_prime', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
                'tax_regular', zeros(cS.nk, cS.nkpps, cS.nw_expanded));

            if isfield(M_age, 'g_A_t_plus_1')
                g_A_period = M_age.g_A_t_plus_1;
            else
                g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
            end

            growth_factor_bgp = (1 + g_A_period);
            discount_factor_V_prime = cS.beta * growth_factor_bgp^(1 - cS.sigma);

            k_grid_vec = cS.kGridV;
            kpps_grid_vec = cS.kppsGridV;
            vPrime_interpolants = cell(cS.nw_expanded, 1);
            if a_idx < cS.aD_new
                for ie_next = 1:cS.nw_expanded
                    vPrime_interpolants{ie_next} = griddedInterpolant({k_grid_vec, kpps_grid_vec}, vPrime_kkppse_next(:, :, ie_next), 'linear', 'linear');
                end
            end
            
            num_tasks = cS.nkpps * cS.nw_expanded;
            val_results_cell = cell(num_tasks, 1);
            pol_results_cell = cell(num_tasks, 1);
            
            parfor task_idx = 1:num_tasks
                [ikpps, ie] = ind2sub([cS.nkpps, cS.nw_expanded], task_idx);
                kpps_hat_current = kpps_grid_vec(ikpps);

                % --- [!!! 核心修正 !!!] 所有流量都在 hat 空间内计算 ---
                capital_return_hat_vec = k_grid_vec * (1 + M_age.r_mkt_t);
                labor_income_hat = 0;
                pps_contrib_hat = 0;
                pps_withdrawal_hat = 0;
                pps_tax_hat = 0;
                pension_income_hat = 0;

                if a_idx <= cS.aR_new
                    labor_income_hat = M_age.w_t * cS.ageEffV_new(a_idx) * paramS_age.leGridV(ie);
                    pps_contrib_hat = cS.pps_contrib_rate * labor_income_hat;
                else
                    pension_income_hat = b_age_val;
                    % 提取额是【期初hat资产赚取利息后】的一个比例
                    pps_withdrawal_hat = cS.pps_withdrawal_rate * (kpps_hat_current * (1 + M_age.r_mkt_t));
                    pps_tax_hat = cS.pps_tax_rate_withdrawal * pps_withdrawal_hat;
                end

                capital_tax_hat_vec = cS.tau_k .* (k_grid_vec * M_age.r_mkt_t);
                labor_tax_hat = cS.tau_l * labor_income_hat;
                payg_contribution_hat = M_age.theta_t * labor_income_hat;

                % [!!! 核心修正 !!!] 完全基于 hat 变量的预算约束
                cash_on_hand_hat_vec = capital_return_hat_vec + labor_income_hat + M_age.tr_per_hh + pension_income_hat + pps_withdrawal_hat ...
                    - (capital_tax_hat_vec + labor_tax_hat + pps_contrib_hat + payg_contribution_hat + pps_tax_hat);

                % [!!! 核心修正 !!!] 正确的 PPS 资产积累方程
                kpps_hat_end_of_period = kpps_hat_current * (1 + M_age.r_mkt_t) + pps_contrib_hat - pps_withdrawal_hat;
                kpps_prime_hat = max(cS.kppsMin, kpps_hat_end_of_period / growth_factor_bgp);

                % --- VFI 核心逻辑 (与 VFI_by_age 保持一致) ---
                k_prime_max_vec = (cash_on_hand_hat_vec - cS.cFloor * (1 + cS.tau_c)) / growth_factor_bgp;
                k_prime_max_vec(k_prime_max_vec < cS.kMin) = cS.kMin;
                k_prime_choices_mat = cS.kMin + (k_prime_max_vec - cS.kMin) .* linspace(0, 1, nk_search_grid);
                c_expend_choices_mat = cash_on_hand_hat_vec - k_prime_choices_mat * growth_factor_bgp;

                c_choices_mat = c_expend_choices_mat / (1 + cS.tau_c);
                invalid_choice_mask = (c_choices_mat < cS.cFloor);
                c_choices_mat(invalid_choice_mask) = cS.cFloor;
                [~, util_c_mat] = utils.CES_utility(c_choices_mat, cS.sigma, cS);
                util_c_mat(invalid_choice_mask) = -1e20;
                future_value_mat = zeros(cS.nk, nk_search_grid);
                
                if a_idx < cS.aD_new
                    trans_mat_row = paramS_age.TrProbM_by_age{a_idx + 1}(ie, :);
                    k_prime_choices_flat = k_prime_choices_mat(:);
                    kpps_prime_vec = ones(length(k_prime_choices_flat), 1) * kpps_prime_hat;
                    v_prime_interp_mat = zeros(length(k_prime_choices_flat), cS.nw_expanded);
                    for ie_next = 1:cS.nw_expanded
                        v_prime_interp_mat(:, ie_next) = vPrime_interpolants{ie_next}(k_prime_choices_flat, kpps_prime_vec);
                    end
                    ev_flat_vec = v_prime_interp_mat * trans_mat_row';
                    ev_on_choices_mat = reshape(ev_flat_vec, cS.nk, nk_search_grid);
                    
                    total_wealth_prime_mat = k_prime_choices_mat + kpps_prime_hat;
                    util_bequest_mat = utils.bequest_utility(total_wealth_prime_mat, cS);
                    
                    survival_rate = cS.s_pathV(a_idx);
                    if size(survival_rate,2) > 1, survival_rate = survival_rate(1); end
                    
                    future_value_mat = discount_factor_V_prime * (survival_rate * ev_on_choices_mat + (1 - survival_rate) * util_bequest_mat);
                else
                    total_wealth_prime_mat = k_prime_choices_mat + kpps_prime_hat;
                    future_value_mat = discount_factor_V_prime * utils.bequest_utility(total_wealth_prime_mat, cS);
                end
                
                total_value_mat = util_c_mat + future_value_mat;

                [best_val_vec, best_idx_vec] = max(total_value_mat, [], 2);
                linear_indices = sub2ind(size(k_prime_choices_mat), (1:cS.nk)', best_idx_vec);
                
                % 保存结果
                c_final_vec = c_choices_mat(linear_indices);
                tax_final_vec = capital_tax_hat_vec + labor_tax_hat + payg_contribution_hat + pps_tax_hat + c_final_vec * cS.tau_c;

                val_results_cell{task_idx} = best_val_vec;
                pol_results_cell{task_idx} = struct(...
                    'c', c_final_vec, ...
                    'k_prime', k_prime_choices_mat(linear_indices), ...
                    'kpps_prime', repmat(kpps_prime_hat, cS.nk, 1), ...
                    'tax_regular', tax_final_vec);
            end
            
            % 将并行计算结果重组回矩阵形式
            fields = fieldnames(pol_age);
            for task_idx = 1:num_tasks
                [ikpps, ie] = ind2sub([cS.nkpps, cS.nw_expanded], task_idx);
                val_age(:, ikpps, ie) = val_results_cell{task_idx};
                current_pol_slice = pol_results_cell{task_idx};
                for i = 1:length(fields)
                    field = fields{i};
                    pol_age.(field)(:, ikpps, ie) = current_pol_slice.(field);
                end
            end
       end
       
       
       function [val_age, pol_age] = VFI_by_age(a_idx, vPrime_kkppse_next, M_age, b_age_val, paramS_age, cS)
            % =========================================================================
            % == 函数: VFI_by_age (版本 v1.9 - 统一贴现修正版)
            % == 核心修正:
            % ==   - 将所有来自未来的效用（期望价值函数 EV' 和遗赠效用 U_bequest）
            % ==     都统一用一个贴现因子 discount_factor_V_prime 来处理。
            % ==   - discount_factor_V_prime = cS.beta * (1 + g_A_period)^(1 - cS.sigma)
            % ==   - 这个因子正确地将下一期的、去趋势化的价值，贴现回当期的可比单位。
            % ==     这从根本上统一了所有未来效用的处理方式，确保了BGP的动态一致性。
            % =========================================================================
            nk_search_grid = 150;
            val_age = -1e20 * ones(cS.nk, 1, cS.nw_expanded);
            pol_age = struct('c', zeros(cS.nk, 1, cS.nw_expanded), ...
                'k_prime', zeros(cS.nk, 1, cS.nw_expanded), ...
                'tax_regular', zeros(cS.nk, 1, cS.nw_expanded));

            if isfield(M_age, 'g_A_t_plus_1') 
                g_A_period = M_age.g_A_t_plus_1;
            else
                g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
            end
            
            growth_factor_bgp = (1 + g_A_period);
            % [!!! 核心修正 !!!] 定义一个统一的、用于所有未来效用项的贴现因子
            discount_factor_V_prime = cS.beta * growth_factor_bgp^(1 - cS.sigma);


            k_grid_vec = cS.kGridV;
            vPrime_interpolants = cell(cS.nw_expanded, 1);
            if a_idx < cS.aD_new
                vPrime_slice = squeeze(vPrime_kkppse_next(:, 1, :));
                for ie_next = 1:cS.nw_expanded
                    vPrime_interpolants{ie_next} = griddedInterpolant(k_grid_vec, vPrime_slice(:, ie_next), 'pchip', 'pchip');
                end
            end
            for ie = 1:cS.nw_expanded
                capital_return = k_grid_vec * (1 + M_age.r_mkt_t);
                labor_income_gross = 0;
                pension_income = 0;
                if a_idx <= cS.aR_new
                    labor_income_gross = M_age.w_t * cS.ageEffV_new(a_idx) * paramS_age.leGridV(ie);
                else
                    pension_income = b_age_val;
                end

                capital_tax_vec = cS.tau_k .* (k_grid_vec * M_age.r_mkt_t);
                labor_tax = cS.tau_l * labor_income_gross;
                payg_contribution = M_age.theta_t * labor_income_gross;
                cash_on_hand_vec = capital_return + labor_income_gross + M_age.tr_per_hh + pension_income - (capital_tax_vec + labor_tax + payg_contribution);
                k_prime_max_vec = (cash_on_hand_vec - cS.cFloor * (1 + cS.tau_c)) / growth_factor_bgp;
                k_prime_max_vec(k_prime_max_vec < cS.kMin) = cS.kMin;
                k_prime_choices_mat = cS.kMin + (k_prime_max_vec - cS.kMin) .* linspace(0, 1, nk_search_grid);
                c_expend_choices_mat = cash_on_hand_vec - k_prime_choices_mat * growth_factor_bgp;

                c_choices_mat = c_expend_choices_mat / (1 + cS.tau_c);
                invalid_choice_mask = (c_choices_mat < cS.cFloor);
                c_choices_mat(invalid_choice_mask) = cS.cFloor;
                [~, util_c_mat] = utils.CES_utility(c_choices_mat, cS.sigma, cS);
                util_c_mat(invalid_choice_mask) = -1e20;
                future_value_mat = zeros(cS.nk, nk_search_grid);
                if a_idx < cS.aD_new
                    trans_mat_row = paramS_age.TrProbM_by_age{a_idx + 1}(ie, :);
                    Vprime_interp_mat = cell2mat(cellfun(@(f) f(k_prime_choices_mat(:)), vPrime_interpolants, 'UniformOutput', false)');
                    ev_flat_vec = Vprime_interp_mat * trans_mat_row';
                    ev_on_choices_mat = reshape(ev_flat_vec, cS.nk, nk_search_grid);
                    survival_rate = cS.s_pathV(a_idx);
                    if size(survival_rate,2) > 1, survival_rate = survival_rate(1); end
                    util_bequest = utils.bequest_utility(k_prime_choices_mat, cS);

                    % [!!! 核心修正 !!!] 将来自未来的期望价值和遗赠效用统一用 discount_factor_V_prime 贴现
                    future_value_mat = discount_factor_V_prime * (survival_rate * ev_on_choices_mat + (1 - survival_rate) * util_bequest);
                else
                    % [!!! 核心修正 !!!] 将来自未来的遗赠效用统一用 discount_factor_V_prime 贴现
                    future_value_mat = discount_factor_V_prime * utils.bequest_utility(k_prime_choices_mat, cS);
                end

                % [!!! 核心修正 !!!] 当期价值 = 当期效用 + 已被正确贴现的未来价值
                total_value_mat = util_c_mat + future_value_mat;
                
                [best_val_vec, best_idx_vec] = max(total_value_mat, [], 2);
                linear_indices = sub2ind(size(k_prime_choices_mat), (1:cS.nk)', best_idx_vec);
                val_age(:, 1, ie) = best_val_vec;
                pol_age.c(:, 1, ie) = c_choices_mat(linear_indices);
                pol_age.k_prime(:, 1, ie) = k_prime_choices_mat(linear_indices);

                pol_age.tax_regular(:, 1, ie) = capital_tax_vec + labor_tax + payg_contribution + c_choices_mat(linear_indices) .* cS.tau_c;
            end
        end
        
        
        function [Pol_path, Val_path] = backward_hh(pathS, cS, paramSF, valF, polF)
            % =========================================================================
            % == 函数: backward_hh (版本 v1.7 - TR人均量接口修正版)
            % ==
            % == 核心修正:
            % ==   - 从 pathS 中读取 tr_per_hh_hat_path (由TPI.m提前计算好的人均量)。
            % ==   - 将 M_vfi_t.tr_per_hh 赋值为对应时期的人均转移支付。
            % ==   - 删除了不再需要的 g_n_path 计算和传递。
            % =========================================================================
            T = cS.T_sim;
            Pol_path = cell(1, T);
            Val_path = zeros(cS.nk, cS.nkpps, cS.nw_expanded, cS.aD_new, T);
            Val_path(:,:,:,:,T) = valF;
            Pol_path{T} = polF;
            A_path_extended = [cS.A_path, cS.A_path(end) * (1 + ((1+cS.g_A_ss)^cS.time_Step-1))];

            for t = (T-1) : -1 : 1
                g_A_onestep_ahead = (A_path_extended(t+1) / A_path_extended(t)) - 1;
                M_vfi_t = struct('r_mkt_t', pathS.r_path(t), ...
                    'w_t', pathS.w_hat_path(t), ...
                    'g_A_t_plus_1', g_A_onestep_ahead);

                % [!!! 核心修正 !!!]
                % 直接使用由TPI.m计算好的人均转移支付路径
                M_vfi_t.tr_per_hh = pathS.tr_per_hh_hat_path(t);
                M_vfi_t.theta_t = pathS.theta_path(t);

                cS_t = cS;
                if size(cS.s_pathV, 2) > 1
                    cS_t.s_pathV = cS.s_pathV(:, t);
                end
                Vprime_t_plus_1 = Val_path(:,:,:,:,t+1);

                [valS_t, polS_t] = household.VFI_transition_engine_by_age(M_vfi_t, pathS.b_hat_path(t), paramSF, cS_t, Vprime_t_plus_1);

                Val_path(:,:,:,:,t) = valS_t;
                Pol_path{t} = polS_t;
            end
        end
        
        
        
        
        function [valS_t, polS_t] = VFI_transition_engine_by_age(M_vfi_t, b_payg_hat_t, paramS_t, cS_t, Vprime_t_plus_1)
            % =========================================================================
            % == 函数: VFI_transition_engine_by_age (版本 v1.1)
            % == 目的:
            % ==   - 作为转轨路径VFI的顶层调度器。
            % ==   - 按年龄从大到小(aD -> 1)循环，为每个年龄组调用核心VFI引擎。
            % ==   - 负责正确分配特定年龄的养老金福利。
            % == 注意:
            % ==   - 此函数本身不包含经济学逻辑，仅做调度和数据传递。
            % ==   - 核心的经济学计算发生在 VFI_by_age 和 VFI_by_age_PPS 中。
            % =========================================================================
            valS_t = -Inf(cS_t.nk, cS_t.nkpps, cS_t.nw_expanded, cS_t.aD_new);
            polS_cell = cell(cS_t.aD_new, 1);

            bV_payg_vfi_hat = zeros(1, cS_t.aD_new);
            if cS_t.aR_new < cS_t.aD_new
                bV_payg_vfi_hat((cS_t.aR_new+1):cS_t.aD_new) = b_payg_hat_t;
            end

            for a_idx = cS_t.aD_new : -1 : 1
                vPrime_next_age = [];
                if a_idx < cS_t.aD_new
                    vPrime_next_age = Vprime_t_plus_1(:,:,:,a_idx+1);
                end
                b_age_val_hat = bV_payg_vfi_hat(a_idx);

                if isfield(cS_t, 'pps_active') && cS_t.pps_active
                    [val_age, pol_age] = household.VFI_by_age_PPS(a_idx, vPrime_next_age, M_vfi_t, b_age_val_hat, paramS_t, cS_t);
                else
                    [val_age, pol_age] = household.VFI_by_age(a_idx, vPrime_next_age, M_vfi_t, b_age_val_hat, paramS_t, cS_t);
                end

                val_reshaped = reshape(val_age, [cS_t.nk, cS_t.nkpps, cS_t.nw_expanded]);
                valS_t(:,:,:,a_idx) = val_reshaped;
                polS_cell{a_idx} = pol_age;
            end
            polS_t = [polS_cell{:}];
        end


    end
end
