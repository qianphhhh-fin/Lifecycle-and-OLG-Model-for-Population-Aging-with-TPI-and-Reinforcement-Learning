% --- household.m ---

classdef household
    methods (Static)


        %         function [val_age, pol_age] = VFI_by_age_PPS(a_idx, vPrime_kkppse_next, M_age, b_age_val, paramS_age, cS)
        %     % = =========================================================================
        %     % == 函数: VFI_by_age_PPS_endo (版本 v1.4 - 插值维度匹配修正版)
        %     % 内生PPS缴费, cS.pps_fixed为缴费上限.
        %     % == 核心修改:
        %     % ==   - [!!! 关键BUG修复 !!!] 修正了调用 griddedInterpolant 时，
        %     % ==     坐标输入数组 (k_prime 和 kpps_prime) 维度不匹配的错误。
        %     % ==   - 问题: griddedInterpolant 不支持隐式广播。k_prime_choices_tens
        %     % ==     的维度包含 nk_search，而 kpps_prime_hat_tens 不包含，导致
        %     % ==     二者 size 不同。
        %     % ==   - 解决方案: 在计算完 k_prime_choices_tens 之后，立即创建一个
        %     % ==     kpps_prime_hat_tens 的显式广播版本 kpps_prime_hat_tens_matched。
        %     % ==     这通过 kpps_prime_hat_tens + zeros(size(k_prime_choices_tens))
        %     % ==     实现，确保了 kpps_prime_hat_tens_matched 的维度与
        %     % ==     k_prime_choices_tens 完全一致。随后所有需要用到 kpps_prime
        %     % ==     的地方都使用这个匹配版本。
        %     % =========================================================================
        %     nk_search_grid = 150;
        % 
        %     % 1. 获取时期增长因子和贴现因子
        %     if isfield(M_age, 'g_A_t_plus_1')
        %         g_A_period = M_age.g_A_t_plus_1;
        %     else
        %         g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
        %     end
        %     growth_factor_bgp = (1 + g_A_period);
        %     discount_factor_V_prime = cS.beta * growth_factor_bgp^(1 - cS.sigma);
        % 
        %     % 2. 为下一期价值函数创建插值器
        %     k_grid_vec = cS.kGridV;
        %     kpps_grid_vec = cS.kppsGridV;
        %     vPrime_interpolants = cell(cS.nw_expanded, 1);
        %     if a_idx < cS.aD_new
        %         for ie_next = 1:cS.nw_expanded
        %             if cS.nkpps > 1 && cS.pps_active
        %                 vPrime_interpolants{ie_next} = griddedInterpolant({k_grid_vec, kpps_grid_vec}, vPrime_kkppse_next(:, :, ie_next), 'makima','makima');
        %             else
        %                 vPrime_interpolants{ie_next} = griddedInterpolant(k_grid_vec, squeeze(vPrime_kkppse_next(:, :, ie_next)), 'pchip','pchip');
        %             end
        %         end
        %     end
        % 
        %     % 3. 构建状态变量网格 (k, kpps, e)
        %     [k_mat, kpps_mat, le_mat] = ndgrid(cS.kGridV, cS.kppsGridV, paramS_age.leGridV);
        % 
        %     % 4. 计算与决策无关的部分
        %     labor_supply_mat = 0;
        %     is_working_age = (a_idx <= cS.aR_new);
        %     if is_working_age
        %         labor_supply_mat = cS.ageEffV_new(a_idx) .* le_mat;
        %     end
        %     labor_income_hat_mat = M_age.w_t * labor_supply_mat;
        % 
        %     capital_return_hat_mat = k_mat * (1 + M_age.r_mkt_t);
        %     capital_tax_hat_mat = cS.tau_k .* (k_mat * M_age.r_mkt_t);
        %     labor_tax_hat_mat = cS.tau_l * labor_income_hat_mat;
        %     payg_contribution_hat_mat = M_age.theta_t * labor_income_hat_mat;
        % 
        %     % 退休期收入与PPS提取 (与PPS缴费决策无关)
        %     pension_income_hat_mat = zeros(size(k_mat));
        %     pps_withdrawal_hat_mat = zeros(size(k_mat));
        %     pps_tax_hat_mat = zeros(size(k_mat));
        %     if ~is_working_age
        %         pension_income_hat_mat(:) = b_age_val;
        %         if cS.pps_active
        %             pps_withdrawal_hat_mat = cS.pps_withdrawal_rate .* (kpps_mat * (1 + M_age.r_mkt_t));
        %             pps_tax_hat_mat = cS.pps_tax_rate_withdrawal .* pps_withdrawal_hat_mat;
        %         end
        %     end
        % 
        %     % 5. 构建选择网格 (k' 和 pps_rate) 并进行向量化计算
        %     % 5.1 构建 pps_rate 选择网格
        %     if is_working_age && cS.pps_active && cS.pps_fixed > 0 && cS.n_pps_rate_grid > 1 % 这里cS.pps_fixed相当于cS.pps_max！！！！！！！！！！！！
        %         n_pps_rate_grid_eff = cS.n_pps_rate_grid;
        %         pps_rate_choices_vec = linspace(0, cS.pps_fixed, n_pps_rate_grid_eff);
        %     else
        %         n_pps_rate_grid_eff = 1;
        %         pps_rate_choices_vec = 0;
        %     end
        %     pps_rate_choices_reshaped = reshape(pps_rate_choices_vec, [1, 1, 1, 1, n_pps_rate_grid_eff]);
        % 
        %     % 5.2 计算依赖 pps_rate 选择的量
        %     pps_contrib_hat_tens = labor_income_hat_mat .* pps_rate_choices_reshaped;
        %     kpps_hat_end_of_period_tens = kpps_mat * (1 + M_age.r_mkt_t) + pps_contrib_hat_tens - pps_withdrawal_hat_mat;
        %     kpps_prime_hat_tens = max(cS.kppsMin, kpps_hat_end_of_period_tens / growth_factor_bgp);
        % 
        %     base_coh_hat_mat = capital_return_hat_mat + labor_income_hat_mat + M_age.tr_per_hh + pension_income_hat_mat + pps_withdrawal_hat_mat ...
        %         - (capital_tax_hat_mat + labor_tax_hat_mat + payg_contribution_hat_mat + pps_tax_hat_mat);
        % 
        %     cash_on_hand_hat_tens = base_coh_hat_mat - pps_contrib_hat_tens;
        % 
        %     % 5.3 构建 k' 选择网格
        %     k_prime_max_tens = (cash_on_hand_hat_tens - cS.cFloor * (1 + cS.tau_c)) / growth_factor_bgp;
        %     k_prime_max_tens(k_prime_max_tens < cS.kMin) = cS.kMin;
        % 
        %     k_prime_search_grid = linspace(0, 1, nk_search_grid);
        %     k_prime_search_grid_reshaped = reshape(k_prime_search_grid, [1, 1, 1, nk_search_grid, 1]);
        % 
        %     k_prime_choices_tens = cS.kMin + (k_prime_max_tens - cS.kMin) .* k_prime_search_grid_reshaped;
        % 
        %     % [!!! 核心修正 !!!] 创建一个与 k_prime_choices_tens 维度完全匹配的 kpps_prime 张量
        %     kpps_prime_hat_tens_matched = kpps_prime_hat_tens + zeros(size(k_prime_choices_tens));
        % 
        %     % 5.4 计算消费和当期效用
        %     c_expend_choices_tens = cash_on_hand_hat_tens - k_prime_choices_tens * growth_factor_bgp;
        %     c_choices_tens = c_expend_choices_tens / (1 + cS.tau_c);
        %     invalid_choice_mask = (c_choices_tens < cS.cFloor);
        %     c_choices_tens(invalid_choice_mask) = cS.cFloor;
        % 
        %     util_c_tens = (c_choices_tens.^(1 - cS.sigma))./(1 - cS.sigma);
        %     util_c_tens(invalid_choice_mask) = -1e20;
        % 
        %     % 5.5 计算遗赠效用
        %     total_wealth_prime_tens = k_prime_choices_tens + kpps_prime_hat_tens_matched;
        %     util_bequest_tens = utils.bequest_utility(total_wealth_prime_tens, cS);
        % 
        %     % 5.6 计算未来价值的期望
        %     future_value_tens = zeros(size(k_prime_choices_tens));
        %     if a_idx < cS.aD_new
        %         ev_on_choices_tens = zeros(size(k_prime_choices_tens));
        %         trans_mat = paramS_age.TrProbM_by_age{a_idx};
        % 
        %         for ie_next = 1:cS.nw_expanded
        %             trans_prob_col = trans_mat(:, ie_next);
        %             if any(trans_prob_col > 0)
        %                 if cS.nkpps > 1 && cS.pps_active
        %                     % [!!! 核心修正 !!!] 使用维度匹配的 kpps_prime 张量
        %                     v_prime_interp_tens = vPrime_interpolants{ie_next}(k_prime_choices_tens, kpps_prime_hat_tens_matched);
        %                 else
        %                     v_prime_interp_tens = vPrime_interpolants{ie_next}(k_prime_choices_tens);
        %                 end
        %                 ev_on_choices_tens = ev_on_choices_tens + v_prime_interp_tens .* reshape(trans_prob_col, [1, 1, cS.nw_expanded, 1, 1]);
        %             end
        %         end
        % 
        %         survival_rate = cS.s_pathV(a_idx);
        %         if size(survival_rate,2) > 1, survival_rate = survival_rate(1); end
        %         future_value_tens = discount_factor_V_prime * (survival_rate * ev_on_choices_tens + (1 - survival_rate) * util_bequest_tens);
        %     else
        %         future_value_tens = discount_factor_V_prime * util_bequest_tens;
        %     end
        % 
        %     total_value_tens = util_c_tens + future_value_tens;
        % 
        %     % 6. 寻找联合最优决策
        %     if n_pps_rate_grid_eff > 1
        %         S_5d = size(total_value_tens);
        %         total_value_reshaped = reshape(total_value_tens, [S_5d(1), S_5d(2), S_5d(3), S_5d(4)*S_5d(5)]);
        %         [val_age, best_joint_idx_mat] = max(total_value_reshaped, [], 4, 'omitnan');
        %         [best_k_prime_idx_mat, best_pps_rate_idx_mat] = ind2sub([nk_search_grid, n_pps_rate_grid_eff], best_joint_idx_mat);
        %     else
        %         [val_age, best_k_prime_idx_mat] = max(total_value_tens, [], 4, 'omitnan');
        %         best_pps_rate_idx_mat = ones(size(val_age));
        %     end
        % 
        %     % 7. 提取最优策略
        %     S = size(c_choices_tens);
        %     [I, J, K] = ndgrid(1:S(1), 1:S(2), 1:S(3));
        % 
        %     if n_pps_rate_grid_eff > 1
        %         linear_indices_final = sub2ind(S, I, J, K, best_k_prime_idx_mat, best_pps_rate_idx_mat);
        %     else
        %         linear_indices_final = sub2ind(S, I, J, K, best_k_prime_idx_mat);
        %     end
        % 
        %     c_final_mat = c_choices_tens(linear_indices_final);
        %     k_prime_final_mat = k_prime_choices_tens(linear_indices_final);
        %     % [!!! 核心修正 !!!] 使用维度匹配的张量进行最终索引
        %     kpps_prime_final_mat = kpps_prime_hat_tens_matched(linear_indices_final);
        % 
        %     tax_consumption_mat = c_final_mat .* cS.tau_c;
        %     tax_final_mat = capital_tax_hat_mat + labor_tax_hat_mat + pps_tax_hat_mat + payg_contribution_hat_mat + tax_consumption_mat;
        % 
        %     % 8. 保存结果到 policy struct
        %     pol_age = struct(...
        %         'c', c_final_mat, ...
        %         'l', labor_supply_mat, ...
        %         'k_prime', k_prime_final_mat, ...
        %         'kpps_prime', kpps_prime_final_mat, ...
        %         'tax_regular', tax_final_mat, ...
        %         'pps_contrib_rate', pps_rate_choices_vec(best_pps_rate_idx_mat)); % 记录最优缴费率
        % end

        function [val_age, pol_age] = VFI_by_age_PPS(a_idx, vPrime_kkppse_next, M_age, b_age_val, paramS_age, cS)
            % = =========================================================================
            % == 函数: VFI_by_age_PPS (版本 v1.5 - EET税收递延精确修正版)
            % == 核心修改:
            % ==   - [!!! 关键理论修正 !!!] 修正了原函数未能正确实现个人养老金
            % ==     缴费税收递延 (EET) 的根本性错误。
            % ==   - 原逻辑: 劳动税基于总劳动收入计算 (Tax = tau_l * Y)，PPS缴费
            % ==     被错误地当作一项税后支出。
            % ==   - 新逻辑: 劳动税的计算被移至PPS缴费决策之后。对于每一个可能
            % ==     的PPS缴费选择Q, 劳动税都基于正确的税基 (Y - Q) 进行计算,
            % ==     即 Tax = tau_l * (Y - Q)。这确保了模型的预算约束与EET
            % ==     理论完全匹配，正确地反映了税收激励。
            % ==   - 保持了对无PPS情况 (pps_active=false 或 pps_fixed=0) 的兼容性。
            % =========================================================================
            nk_search_grid = 150;

            % 1. 获取时期增长因子和贴现因子
            if isfield(M_age, 'g_A_t_plus_1')
                g_A_period = M_age.g_A_t_plus_1;
            else
                g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
            end
            growth_factor_bgp = (1 + g_A_period);
            discount_factor_V_prime = cS.beta * growth_factor_bgp^(1 - cS.sigma);

            % 2. 为下一期价值函数创建插值器
            k_grid_vec = cS.kGridV;
            kpps_grid_vec = cS.kppsGridV;
            vPrime_interpolants = cell(cS.nw_expanded, 1);
            if a_idx < cS.aD_new
                for ie_next = 1:cS.nw_expanded
                    if cS.nkpps > 1 && cS.pps_active
                        vPrime_interpolants{ie_next} = griddedInterpolant({k_grid_vec, kpps_grid_vec}, vPrime_kkppse_next(:, :, ie_next), 'makima','makima');
                    else
                        vPrime_interpolants{ie_next} = griddedInterpolant(k_grid_vec, squeeze(vPrime_kkppse_next(:, :, ie_next)), 'pchip','pchip');
                    end
                end
            end

            % 3. 构建状态变量网格 (k, kpps, e)
            [k_mat, kpps_mat, le_mat] = ndgrid(cS.kGridV, cS.kppsGridV, paramS_age.leGridV);

            % 4. 计算不直接依赖于当期PPS缴费决策的收入与支出项
            is_working_age = (a_idx <= cS.aR_new);
            labor_supply_mat = 0;
            if is_working_age
                labor_supply_mat = cS.ageEffV_new(a_idx) .* le_mat;
            end
            labor_income_hat_mat = M_age.w_t * labor_supply_mat;
            capital_return_hat_mat = k_mat * (1 + M_age.r_mkt_t);
            capital_tax_hat_mat = cS.tau_k .* (k_mat * M_age.r_mkt_t);
            payg_contribution_hat_mat = M_age.theta_t * labor_income_hat_mat;

            pension_income_hat_mat = zeros(size(k_mat));
            pps_withdrawal_hat_mat = zeros(size(k_mat));
            pps_tax_hat_mat = zeros(size(k_mat));
            if ~is_working_age
                pension_income_hat_mat(:) = b_age_val;
                if cS.pps_active
                    pps_withdrawal_hat_mat = cS.pps_withdrawal_rate .* (kpps_mat * (1 + M_age.r_mkt_t));
                    pps_tax_hat_mat = cS.pps_tax_rate_withdrawal .* pps_withdrawal_hat_mat;
                end
            end

            % 5. 构建选择网格 (k' 和 pps_rate)
            % 5.1 构建 pps_rate 选择网格
            if is_working_age && cS.pps_active && cS.pps_fixed > 0
                n_pps_rate_grid_eff = cS.n_pps_rate_grid;
                pps_rate_choices_vec = linspace(0, cS.pps_fixed, n_pps_rate_grid_eff);
            else
                n_pps_rate_grid_eff = 1;
                pps_rate_choices_vec = 0;
            end
            pps_rate_choices_reshaped = reshape(pps_rate_choices_vec, [1, 1, 1, 1, n_pps_rate_grid_eff]);

            % 6. [!!! 核心修正: 实现税收递延 !!!] 计算依赖于PPS缴费决策的现金流
            % 6.1 PPS缴费额 (level)
            pps_contrib_hat_tens = labor_income_hat_mat .* pps_rate_choices_reshaped;
            
            % 6.2 计算应税劳动收入和劳动收入税 (现在是张量)
            taxable_labor_income_tens = labor_income_hat_mat - pps_contrib_hat_tens;
            labor_tax_hat_tens = cS.tau_l .* taxable_labor_income_tens; % 税基是扣除PPS缴费后的收入

            % 6.3 计算总的手头现金 (Cash-on-Hand)
            cash_on_hand_hat_tens = (capital_return_hat_mat + labor_income_hat_mat + M_age.tr_per_hh + pension_income_hat_mat + pps_withdrawal_hat_mat) ...
                                  - (capital_tax_hat_mat + payg_contribution_hat_mat + pps_tax_hat_mat) ...
                                  - (labor_tax_hat_tens + pps_contrib_hat_tens);

            % 7. 计算下一期资产状态 k' 和 kpps'
            kpps_hat_end_of_period_tens = kpps_mat * (1 + M_age.r_mkt_t) + pps_contrib_hat_tens - pps_withdrawal_hat_mat;
            kpps_prime_hat_tens = max(cS.kppsMin, kpps_hat_end_of_period_tens / growth_factor_bgp);

            k_prime_max_tens = (cash_on_hand_hat_tens - cS.cFloor * (1 + cS.tau_c)) / growth_factor_bgp;
            k_prime_max_tens(k_prime_max_tens < cS.kMin) = cS.kMin;
            k_prime_search_grid = linspace(0, 1, nk_search_grid);
            k_prime_search_grid_reshaped = reshape(k_prime_search_grid, [1, 1, 1, nk_search_grid, 1]);
            k_prime_choices_tens = cS.kMin + (k_prime_max_tens - cS.kMin) .* k_prime_search_grid_reshaped;
            kpps_prime_hat_tens_matched = kpps_prime_hat_tens + zeros(size(k_prime_choices_tens));

            % 8. 计算当期效用和未来价值的期望
            c_expend_choices_tens = cash_on_hand_hat_tens - k_prime_choices_tens * growth_factor_bgp;
            c_choices_tens = c_expend_choices_tens / (1 + cS.tau_c);
            invalid_choice_mask = (c_choices_tens < cS.cFloor);
            c_choices_tens(invalid_choice_mask) = cS.cFloor;
            util_c_tens = (c_choices_tens.^(1 - cS.sigma))./(1 - cS.sigma);
            util_c_tens(invalid_choice_mask) = -1e20;

            total_wealth_prime_tens = k_prime_choices_tens + kpps_prime_hat_tens_matched;
            util_bequest_tens = utils.bequest_utility(total_wealth_prime_tens, cS);

            future_value_tens = zeros(size(k_prime_choices_tens));
            if a_idx < cS.aD_new
                ev_on_choices_tens = zeros(size(k_prime_choices_tens));
                trans_mat = paramS_age.TrProbM_by_age{a_idx};

                for ie_next = 1:cS.nw_expanded
                    trans_prob_col = trans_mat(:, ie_next);
                    if any(trans_prob_col > 0)
                        if cS.nkpps > 1 && cS.pps_active
                            v_prime_interp_tens = vPrime_interpolants{ie_next}(k_prime_choices_tens, kpps_prime_hat_tens_matched);
                        else
                            v_prime_interp_tens = vPrime_interpolants{ie_next}(k_prime_choices_tens);
                        end
                        ev_on_choices_tens = ev_on_choices_tens + v_prime_interp_tens .* reshape(trans_prob_col, [1, 1, cS.nw_expanded, 1, 1]);
                    end
                end

                survival_rate = cS.s_pathV(a_idx);
                if size(survival_rate,2) > 1, survival_rate = survival_rate(1); end
                future_value_tens = discount_factor_V_prime * (survival_rate * ev_on_choices_tens + (1 - survival_rate) * util_bequest_tens);
            else
                future_value_tens = discount_factor_V_prime * util_bequest_tens;
            end

            total_value_tens = util_c_tens + future_value_tens;

            % 9. 寻找联合最优决策
            if n_pps_rate_grid_eff > 1
                S_5d = size(total_value_tens);
                total_value_reshaped = reshape(total_value_tens, [S_5d(1), S_5d(2), S_5d(3), S_5d(4)*S_5d(5)]);
                [val_age, best_joint_idx_mat] = max(total_value_reshaped, [], 4, 'omitnan');
                [best_k_prime_idx_mat, best_pps_rate_idx_mat] = ind2sub([nk_search_grid, n_pps_rate_grid_eff], best_joint_idx_mat);
            else
                [val_age, best_k_prime_idx_mat] = max(total_value_tens, [], 4, 'omitnan');
                best_pps_rate_idx_mat = ones(size(val_age));
            end

            % 10. 提取最优策略
            S = size(c_choices_tens);
            [I, J, K] = ndgrid(1:S(1), 1:S(2), 1:S(3));
            if n_pps_rate_grid_eff > 1
                linear_indices_final = sub2ind(S, I, J, K, best_k_prime_idx_mat, best_pps_rate_idx_mat);
            else
                linear_indices_final = sub2ind(S, I, J, K, best_k_prime_idx_mat);
            end

            c_final_mat = c_choices_tens(linear_indices_final);
            k_prime_final_mat = k_prime_choices_tens(linear_indices_final);
            kpps_prime_final_mat = kpps_prime_hat_tens_matched(linear_indices_final);
            pps_rate_final_mat = pps_rate_choices_vec(best_pps_rate_idx_mat);

            % [!!! 核心修正 !!!] 根据最优决策重新计算总税收
            pps_contrib_final_mat = labor_income_hat_mat .* pps_rate_final_mat;
            taxable_labor_income_final_mat = labor_income_hat_mat - pps_contrib_final_mat;
            labor_tax_final_mat = cS.tau_l .* taxable_labor_income_final_mat;
            tax_consumption_mat = c_final_mat .* cS.tau_c;
            tax_final_mat = capital_tax_hat_mat + labor_tax_final_mat + pps_tax_hat_mat + payg_contribution_hat_mat + tax_consumption_mat;

            % 11. 保存结果到 policy struct
            pol_age = struct(...
                'c', c_final_mat, ...
                'l', labor_supply_mat, ...
                'k_prime', k_prime_final_mat, ...
                'kpps_prime', kpps_prime_final_mat, ...
                'tax_regular', tax_final_mat, ...
                'pps_contrib_rate', pps_rate_final_mat);
        end

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

                % if isfield(cS_vfi, 'pps_active') && cS_vfi.pps_active
                %     [val_age, pol_age] = household.VFI_by_age_PPS(a_idx, vPrime_kkppse_next, M_vfi, b_age_val, paramS_vfi, cS_vfi);
                % else
                %     [val_age, pol_age] = household.VFI_by_age(a_idx, vPrime_kkppse_next, M_vfi, b_age_val, paramS_vfi, cS_vfi);
                % end

                [val_age, pol_age] = household.VFI_by_age_PPS(a_idx, vPrime_kkppse_next, M_vfi, b_age_val, paramS_vfi, cS_vfi);

                valS(:,:,:,a_idx) = val_age;
                polS_cell{a_idx} = pol_age;
            end
            polS = [polS_cell{:}];
        end

        function [val_age, pol_age] = VFI_by_age(a_idx, vPrime_kkppse_next, M_age, b_age_val, paramS_age, cS)
            % =========================================================================
            % == 函数: VFI_by_age (版本 v2.3 - 跨期预算约束精确修正版)
            % == 核心修正:
            % ==   - [!!! 关键会计修正 !!!] 明确了 hat 变量的跨期关系。
            % ==   - 当期手头现金 Cash_on_Hand_hat(t) 用于当期消费 c_hat(t) 和
            % ==     【为下一期准备的】储蓄 K_prime_level(t)。
            % ==   - K_prime_level(t) 转换到下一期的 hat 单位是
            % ==     k_prime_hat(t) = K_prime_level(t) / A(t+1)
            % ==                     = (K_prime_level(t) / A(t)) / (A(t+1)/A(t))
            % ==                     = k_prime_hat_intermediate / (1+g_A)
            % ==     此修正确保了家庭面对的预算约束在增长环境下是正确的。
            % =========================================================================
            nk_search_grid = 150;
            val_age = -1e20 * ones(cS.nk, 1, cS.nw_expanded);
            pol_age = struct('c', zeros(cS.nk, 1, cS.nw_expanded), ...
                'l', zeros(cS.nk, 1, cS.nw_expanded), ...
                'k_prime', zeros(cS.nk, 1, cS.nw_expanded), ...
                'tax_regular', zeros(cS.nk, 1, cS.nw_expanded));

            if isfield(M_age, 'g_A_t_plus_1')
                g_A_period = M_age.g_A_t_plus_1;
            else
                g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
            end

            growth_factor_bgp = (1 + g_A_period);
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
                labor_supply_ie = 0;
                if a_idx <= cS.aR_new
                    labor_supply_ie = cS.ageEffV_new(a_idx) * paramS_age.leGridV(ie);
                end
                labor_income_gross = M_age.w_t * labor_supply_ie;

                pension_income = 0;
                if a_idx > cS.aR_new
                    pension_income = b_age_val;
                end

                capital_return = k_grid_vec * (1 + M_age.r_mkt_t);
                capital_tax_vec = cS.tau_k .* (k_grid_vec * M_age.r_mkt_t);
                labor_tax = cS.tau_l * labor_income_gross;
                payg_contribution = M_age.theta_t * labor_income_gross;

                % Cash on hand 是 t 期的 hat 量
                cash_on_hand_vec = capital_return + labor_income_gross + M_age.tr_per_hh + pension_income - (capital_tax_vec + labor_tax + payg_contribution);

                % [!!! 核心会计修正 !!!]
                % 从 CoH(t) 中拿出的用于储蓄的部分，是相对于 A(t) 的水平。
                % 我们要选择的 k_prime_hat(t) 是相对于 A(t+1) 的水平。
                % 所以 CoH_hat(t) = c_hat(t) + k_prime_level(t)/A(t)
                %                = c_hat(t) + (k_prime_hat(t)*A(t+1))/A(t)
                %                = c_hat(t) + k_prime_hat(t)*(1+g_A)
                k_prime_max_vec = (cash_on_hand_vec - cS.cFloor * (1 + cS.tau_c)) / growth_factor_bgp;
                k_prime_max_vec(k_prime_max_vec < cS.kMin) = cS.kMin;

                k_prime_choices_mat = cS.kMin + (k_prime_max_vec - cS.kMin) .* linspace(0, 1, nk_search_grid);

                % 从预算约束反解出消费
                c_expend_choices_mat = cash_on_hand_vec - k_prime_choices_mat * growth_factor_bgp;

                c_choices_mat = c_expend_choices_mat / (1 + cS.tau_c);
                invalid_choice_mask = (c_choices_mat < cS.cFloor);
                c_choices_mat(invalid_choice_mask) = cS.cFloor;

                [~, util_c_mat] = utils.CES_utility(c_choices_mat, cS.sigma, cS);
                util_c_mat(invalid_choice_mask) = -1e20;

                future_value_mat = zeros(cS.nk, nk_search_grid);
                if a_idx < cS.aD_new
                    trans_mat_row = paramS_age.TrProbM_by_age{a_idx}(ie, :);
                    Vprime_interp_mat = cell2mat(cellfun(@(f) f(k_prime_choices_mat(:)), vPrime_interpolants, 'UniformOutput', false)');
                    ev_flat_vec = Vprime_interp_mat * trans_mat_row';
                    ev_on_choices_mat = reshape(ev_flat_vec, cS.nk, nk_search_grid);
                    survival_rate = cS.s_pathV(a_idx);
                    if size(survival_rate,2) > 1, survival_rate = survival_rate(1); end
                    util_bequest = utils.bequest_utility(k_prime_choices_mat, cS);

                    future_value_mat = discount_factor_V_prime * (survival_rate * ev_on_choices_mat + (1 - survival_rate) * util_bequest);
                else
                    future_value_mat = discount_factor_V_prime * utils.bequest_utility(k_prime_choices_mat, cS);
                end

                total_value_mat = util_c_mat + future_value_mat;

                [best_val_vec, best_idx_vec] = max(total_value_mat, [], 2);
                linear_indices = sub2ind(size(k_prime_choices_mat), (1:cS.nk)', best_idx_vec);

                val_age(:, 1, ie) = best_val_vec;
                pol_age.c(:, 1, ie) = c_choices_mat(linear_indices);
                pol_age.k_prime(:, 1, ie) = k_prime_choices_mat(linear_indices);
                pol_age.l(:, 1, ie) = labor_supply_ie;
                pol_age.tax_regular(:, 1, ie) = capital_tax_vec + labor_tax + payg_contribution + c_choices_mat(linear_indices) .* cS.tau_c;
            end
        end


        function [Pol_path, Val_path] = backward_hh(pathS, cS, paramSF, valF, polF)
            % =========================================================================
            % == 函数: backward_hh (v1.2 - 时变退休年龄版)
            % == 核心修改:
            % ==   - [!!!] 新增了对时变退休年龄路径 cS.aR_new_path 的处理。
            % ==   - 在每个时期 t 的循环内部，创建一个临时的 cS_t 结构体。
            % ==   - 将该时期对应的存活率 s_pathV(:,t) 和退休年龄索引 aR_new_path(t)
            % ==     赋值给 cS_t，然后将 cS_t 传递给下一层求解器。
            % ==   - 这确保了VFI引擎在每个时期都使用正确的退休年龄进行决策。
            % =========================================================================
            T = cS.T_sim;
            Pol_path = cell(1, T);
            Val_path = zeros(cS.nk, cS.nkpps, cS.nw_expanded, cS.aD_new, T);
            Val_path(:,:,:,:,T) = valF;
            Pol_path{T} = polF;

            % 保存原始的、可能是多期的生存率路径
            s_pathV_full_path = cS.s_pathV;
            aR_new_full_path = cS.aR_new_path; % 获取退休年龄路径
            if cS.pps_active
            pps_fixed_full_path = cS.pps_fixed_path; % [!!! 新增 !!!] 获取PPS缴费率路径
            end

            for t = (T-1) : -1 : 1
                g_A_onestep_ahead = pathS.g_A_path(t);

                M_vfi_t = struct('r_mkt_t', pathS.r_path(t), ...
                    'w_t', pathS.w_hat_path(t), ...
                    'g_A_t_plus_1', g_A_onestep_ahead);

                M_vfi_t.tr_per_hh = pathS.tr_per_hh_hat_path(t);
                M_vfi_t.theta_t = pathS.theta_path(t);

                % === 核心修改: 创建一个包含时期特定参数的 cS_t ===
                cS_t = cS; % 复制基础参数
                if size(s_pathV_full_path, 2) > 1
                    cS_t.s_pathV = s_pathV_full_path(:, t);
                end
                cS_t.aR_new = aR_new_full_path(t); % 设置当期的退休年龄索引
                if cS.pps_active
                cS_t.pps_fixed = pps_fixed_full_path(t); % [!!! 新增 !!!] 设置当期的PPS缴费率
                end

                Vprime_t_plus_1 = Val_path(:,:,:,:,t+1);

                % 将包含正确生存率和退休年龄的 cS_t 传递下去
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

                % if isfield(cS_t, 'pps_active') && cS_t.pps_active
                [val_age, pol_age] = household.VFI_by_age_PPS(a_idx, vPrime_next_age, M_vfi_t, b_age_val_hat, paramS_t, cS_t);
                % else
                %     [val_age, pol_age] = household.VFI_by_age(a_idx, vPrime_next_age, M_vfi_t, b_age_val_hat, paramS_t, cS_t);
                % end

                val_reshaped = reshape(val_age, [cS_t.nk, cS_t.nkpps, cS_t.nw_expanded]);
                valS_t(:,:,:,a_idx) = val_reshaped;
                polS_cell{a_idx} = pol_age;
            end
            polS_t = [polS_cell{:}];
        end

        function [muM, utilM] = CES_utility(cM, sigma, cS)
            % =========================================================================
            % == 函数: CES_utility (v2.0 - 向量化安全修正版)
            % == 核心修正:
            % ==   - [!!! 关键稳健性修正 !!!] 替换了非向量化的 `if cM<=0`
            % ==     判断。旧代码在处理向量输入时，仅判断第一个元素，
            % ==     当向量中包含非正值时会导致 NaN 或 Inf。
            % ==   - 新方法: 使用 `max(cM, cS.cFloor)` 对消费输入进行
            % ==     向量化处理，确保所有元素都大于等于消费下限，然后再
            % ==     代入效用函数计算。这保证了函数对任意向量输入的健壮性。
            % =========================================================================

            % 使用 max 函数进行向量化安全处理，确保消费不低于下限
            c_adj = max(cM, cS.cFloor);

            if abs(sigma - 1) < 1e-6
                utilM = log(c_adj);
                muM = 1./c_adj;
            else
                utilM = (c_adj.^(1-sigma))./(1-sigma);
                muM = c_adj.^(-sigma);
            end
        end

        function util_beq = bequest_utility(k_prime, cS)
            % =========================================================================
            % == 函数: bequest_utility (v2.0 - 向量化安全修正版)
            % == 核心修正:
            % ==   - [!!! 关键稳健性修正 !!!] 与CES_utility的修正类似，
            % ==     使用 `max(k_prime, cS.cFloor)` 对遗赠财富输入进行
            % ==     向量化处理，以防止非正值输入导致计算错误。
            % ==   - 这确保了函数在任何求解器中（VFI, EGM）的行为都
            % ==     是稳健和可预测的。
            % =========================================================================
            if cS.phi_bequest <= 0
                util_beq = zeros(size(k_prime)); % 确保输出维度正确
                return;
            end

            % 使用 max 函数进行向量化安全处理
            k_adj = max(k_prime, cS.cFloor);

            if abs(cS.sigma - 1) < 1e-6
                util_beq = cS.phi_bequest * log(k_adj);
            else
                util_beq = cS.phi_bequest * (k_adj.^(1-cS.sigma))./(1-cS.sigma);
            end
        end

        function mu_beq = bequest_utility_deriv(k_prime, cS)
            % =========================================================================
            % == 函数: bequest_utility_deriv (版本 v1.0 - 内联实现版)
            % == 目的: 计算遗赠的边际效用 U_beq'(k')。
            % ==       d/dk [ φ * k^(1-σ) / (1-σ) ] = φ * k^(-σ)
            % ==       此函数被内联到 household 类中以消除对 utils 的依赖。
            % =========================================================================
            if cS.phi_bequest <= 0
                mu_beq = zeros(size(k_prime)); % 确保输出维度正确
                return;
            end

            % 使用 max 函数进行向量化安全处理，避免 k_prime <= 0 导致计算错误
            k_adj = max(k_prime, cS.cFloor);

            if abs(cS.sigma - 1) < 1e-6
                % 对数效用形式的导数: φ / k
                mu_beq = cS.phi_bequest ./ k_adj;
            else
                % CRRA效用形式的导数: φ * k^(-σ)
                mu_beq = cS.phi_bequest .* (k_adj.^(-cS.sigma));
            end
        end



        function run_pps_compare()
            % =========================================================================
            % == 测试函数: run_pps_compare (v2.2 - 调用签名修正版)
            % == 核心修改:
            % ==   - [!!! 关键修正 !!!] 根据用户指正，恢复了对
            % ==     utils.EarningProcess_AgeDependent 函数的正确调用签名。
            % ==     该函数返回四个独立的输出，而非两个结构体。
            % ==   - 保持了 v2.1 版本的模块化结构，对 nw=1 和 nw=5 的情况
            % ==     分别执行完整的测试流程。
            % =========================================================================

            fprintf('=======================================================\n');
            fprintf('==  开启 PPS VFI 求解器对比测试程序  ==\n');
            fprintf('==  (外生固定缴费 vs 内生选择缴费)  ==\n');
            fprintf('=======================================================\n\n');

            % --- 1. 初始化通用参数 ---
            cS_base = utils.ParameterValues();
            cS_base.nk = 40;
            cS_base.nkpps = 30;
            cS_base.g_A_ss = 0.015;
            cS_base.ageEffV_new = cS_base.ageEffV_new_h(:,1);
            cS_base.s_pathV = cS_base.s_pathV(:, 1);
            cS_base.pps_active = true;
            M1 = struct('r_mkt_t', 0.04, 'w_t', 1.5, 'tr_per_hh', 0.01, 'b_hat_t', 0.5, 'theta_t', 0.11);

            % --- 2. 运行 nw = 1 的测试案例 ---
            cS1 = cS_base;
            cS1.nw = 1;
            household.run_test_case(cS1, M1);

            % --- 3. 运行 nw = 5 的测试案例 ---
            % cS5 = cS_base;
            % cS5.nw = 5;
            % household.run_test_case(cS5, M1);

            fprintf('\n所有测试完成。\n');
        end

        % --- 内部辅助函数，执行一套完整的对比测试 ---
        function run_test_case(cS, M1)
            nw_val = cS.nw;
            fprintf('\n\n#############################################\n');
            fprintf('###   开始执行测试案例 (nw = %d)   ###\n', nw_val);
            fprintf('#############################################\n\n');

            % --- 1. 生成依赖于 nw 的参数 ---
            paramS = struct(); % 初始化 paramS
            % [!!! 关键修正 !!!] 恢复正确的四输出调用方式
            [paramS.leGridV, paramS.TrProbM_by_age, paramS.leProb1V, cS.nw_expanded] = utils.EarningProcess_AgeDependent(cS);

            b_payg_hat_age_vec = zeros(1, cS.aD_new);
            b_payg_hat_age_vec((cS.aR_new+1):end) = M1.b_hat_t;

            % --- 2. 测试模块 1: 零缴费等价性检验 ---
            fprintf('--- 测试 1 (nw=%d): 零缴费等价性检验 ---\n', nw_val);
            fprintf('   - 模型 A: 外生固定缴费, cS.pps_fixed = 0\n');
            fprintf('   - 模型 B: 内生选择缴费, cS.pps_max = 0\n');

            % --- 模型 A ---
            cS_A = cS;
            cS_A.pps_fixed = 0;
            cS_A = utils.generateGrids(cS_A);
            valS_A = -Inf(cS_A.nk, cS_A.nkpps, cS_A.nw_expanded, cS_A.aD_new);
            polS_A_cell = cell(cS_A.aD_new, 1);
            tic;
            for a_idx = cS_A.aD_new : -1 : 1
                vPrime_kkppse_next = [];
                if a_idx < cS_A.aD_new, vPrime_kkppse_next = valS_A(:,:,:,a_idx+1); end
                [val_age, pol_age] = household.VFI_by_age_PPS(a_idx, vPrime_kkppse_next, M1, b_payg_hat_age_vec(a_idx), paramS, cS_A);
                valS_A(:,:,:,a_idx) = val_age;
                polS_A_cell{a_idx} = pol_age;
            end
            polS_A = [polS_A_cell{:}];
            time_A = toc;
            fprintf('模型 A 求解完成，耗时: %.2f 秒\n', time_A);

            % --- 模型 B ---
            cS_B = cS;
            cS_B.pps_fixed = 0;
            cS_B.n_pps_rate_grid = 1;
            cS_B = utils.generateGrids(cS_B);
            valS_B = -Inf(cS_B.nk, cS_B.nkpps, cS_B.nw_expanded, cS_B.aD_new);
            polS_B_cell = cell(cS_B.aD_new, 1);
            tic;
            for a_idx = cS_B.aD_new : -1 : 1
                vPrime_kkppse_next = [];
                if a_idx < cS_B.aD_new, vPrime_kkppse_next = valS_B(:,:,:,a_idx+1); end
                [val_age, pol_age] = household.VFI_by_age_PPS_endo(a_idx, vPrime_kkppse_next, M1, b_payg_hat_age_vec(a_idx), paramS, cS_B);
                valS_B(:,:,:,a_idx) = val_age;
                polS_B_cell{a_idx} = pol_age;
            end
            polS_B = [polS_B_cell{:}];
            time_B = toc;
            fprintf('模型 B 求解完成，耗时: %.2f 秒\n\n', time_B);

            diff_val = max(abs(valS_A - valS_B), [], 'all');
            fprintf('结果对比 (Max Abs Diff): 价值函数差异: %.4e\n', diff_val);
            if diff_val < 1e-9, fprintf('结论: 等价性验证通过。\n\n'); else fprintf('结论: 等价性验证失败。\n\n'); end

            % --- 3. 测试模块 2: 正缴费对比 ---
            fprintf('--- 测试 2 (nw=%d): 正缴费对比 ---\n', nw_val);
            fprintf('   - 模型 C: 外生固定缴费, cS.pps_fixed = 0.1\n');
            fprintf('   - 模型 D: 内生选择缴费, cS.pps_max = 0.1, n_grid = 5\n');

            % --- 模型 C ---
            cS_C = cS;
            cS_C.pps_fixed = 0.1;
            cS_C = utils.generateGrids(cS_C);
            valS_C = -Inf(cS_C.nk, cS_C.nkpps, cS_C.nw_expanded, cS_C.aD_new);
            polS_C_cell = cell(cS_C.aD_new, 1);
            tic;
            for a_idx = cS_C.aD_new : -1 : 1
                vPrime_kkppse_next = [];
                if a_idx < cS_C.aD_new, vPrime_kkppse_next = valS_C(:,:,:,a_idx+1); end
                [val_age, pol_age] = household.VFI_by_age_PPS(a_idx, vPrime_kkppse_next, M1, b_payg_hat_age_vec(a_idx), paramS, cS_C);
                valS_C(:,:,:,a_idx) = val_age;
                polS_C_cell{a_idx} = pol_age;
            end
            polS_C = [polS_C_cell{:}];
            time_C = toc;
            fprintf('模型 C 求解完成，耗时: %.2f 秒\n', time_C);

            % --- 模型 D ---
            cS_D = cS;
            cS_D.pps_fixed = 0.1;
            cS_D.n_pps_rate_grid = 5;
            cS_D = utils.generateGrids(cS_D);
            valS_D = -Inf(cS_D.nk, cS_D.nkpps, cS_D.nw_expanded, cS_D.aD_new);
            polS_D_cell = cell(cS_D.aD_new, 1);
            tic;
            for a_idx = cS_D.aD_new : -1 : 1
                vPrime_kkppse_next = [];
                if a_idx < cS_D.aD_new, vPrime_kkppse_next = valS_D(:,:,:,a_idx+1); end
                [val_age, pol_age] = household.VFI_by_age_PPS_endo(a_idx, vPrime_kkppse_next, M1, b_payg_hat_age_vec(a_idx), paramS, cS_D);
                valS_D(:,:,:,a_idx) = val_age;
                polS_D_cell{a_idx} = pol_age;
            end
            polS_D = [polS_D_cell{:}];
            time_D = toc;
            fprintf('模型 D 求解完成，耗时: %.2f 秒\n\n', time_D);

            fprintf('结果对比 (全状态空间均值 for nw=%d):\n', nw_val);
            fprintf('%-18s | %-15s | %-15s\n', '变量', '模型 C (外生)', '模型 D (内生)');
            disp(repmat('-', 1, 55));
            fprintf('%-18s | %-15.4f | %-15.4f\n', '价值函数', mean(valS_C,'all'), mean(valS_D,'all'));
            fprintf('%-18s | %-15.4f | %-15.4f\n', '消费', mean(cell2mat(arrayfun(@(x) x.c(:), polS_C, 'UniformOutput', false))), mean(cell2mat(arrayfun(@(x) x.c(:), polS_D, 'UniformOutput', false))));
            fprintf('%-18s | %-15.4f | %-15.4f\n', '常规储蓄 (k_prime)', mean(cell2mat(arrayfun(@(x) x.k_prime(:), polS_C, 'UniformOutput', false))), mean(cell2mat(arrayfun(@(x) x.k_prime(:), polS_D, 'UniformOutput', false))));
            fprintf('%-18s | %-15s | %-15.4f\n', 'PPS缴费率', '0.1 (固定)', mean(cell2mat(arrayfun(@(x) x.pps_contrib_rate(:), polS_D, 'UniformOutput', false))));
        end


        function run_vfi_test()
            % =========================================================================
            % == 测试函数: run_vfi_test (V2.2 - 测试逻辑修正版)
            % == 核心变更:
            % ==   - [!!! 关键测试逻辑修正 !!!] 在“有PPS”的测试案例中，当
            % ==     目标是验证 pps_max = 0 的情况时，必须将 cS.nkpps 也设为1。
            % ==   - 这确保了 kpps 状态维度被坍缩至一个点(0)，使得该测试与
            % ==     “无PPS”的案例在经济上是等价和可比的。
            % ==   - 此前的测试将 nkpps 设为20，导致即使在0缴费率下，家庭依然
            % ==     拥有外生的、增长的kpps财富，从而造成结果不一致。
            % =========================================================================

            fprintf('============================================\n');
            fprintf('==  开启 household.VFI_solver 测试程序  ==\n');
            fprintf('============================================\n\n');

            % --- 1. 初始化和参数设置 ---
            cS = utils.ParameterValues();

            % 为了加速测试，可以适当减小网格密度
            cS.nk = 30;
            cS.nkpps = 1; % 默认值
            cS.nw = 1;
            cS.n_pps_rate_grid = 20;
            cS.pps_max = 0.1; % 关键测试参数

            [paramS.leGridV, paramS.TrProbM_by_age, paramS.leProb1V, cS.nw_expanded] = utils.EarningProcess_AgeDependent(cS);
            cS.s_pathV = cS.s_pathV(:, 1);

            M1 = struct();
            M1.r_mkt_t = 0.04;
            M1.w_t = 1.5;
            M1.tr_per_hh = 0.01;
            M1.b_hat_t = 0.5;
            M1.theta_t = 0.11;

            % --- 2. 测试案例 1: 无 PPS 系统 ---
            fprintf('--- 测试 1: 模型无 PPS (pps_active = false) ---\n');
            cS.pps_active = false;
            cS.nkpps = 1;
            cS_no_pps = utils.generateGrids(cS);

            tic;
            [polS1, valS1] = household.VFI_solver(M1, paramS, cS_no_pps);
            elapsed_time = toc;
            fprintf('求解完成，耗时: %.2f 秒\n', elapsed_time);

            ages_to_check = [3, cS.aR_new, cS.aD_new-1];
            fprintf('%-10s %-15s %-15s %-15s\n', '年龄(a)', 'Mean(Value)', 'Mean(c)', 'Mean(k_prime)');
            for age = ages_to_check
                val_mean = mean(valS1(:,:,:,age), 'all');
                c_mean = mean(polS1(age).c, 'all');
                k_prime_mean = mean(polS1(age).k_prime, 'all');
                fprintf('%-10d %-15.4f %-15.4f %-15.4f\n', age, val_mean, c_mean, k_prime_mean);
            end

            % --- 3. 测试案例 2: 有 PPS 系统 (但配置为与无PPS等价) ---
            fprintf('\n--- 测试 2: 模型有 PPS (pps_active = true, 但 pps_max=0, nkpps=1) ---\n');
            cS.pps_active = true;
            % [!!! 关键测试逻辑修正 !!!]
            % 为使测试2与测试1可比，当 pps_max=0 时，nkpps 必须为1
            cS.nkpps = 30;
            cS_with_pps = utils.generateGrids(cS);

            tic;
            [polS2, valS2] = household.VFI_solver(M1, paramS, cS_with_pps);
            elapsed_time = toc;
            fprintf('求解完成，耗时: %.2f 秒\n', elapsed_time);

            fprintf('%-10s %-15s %-15s %-15s %-15s\n', '年龄(a)', 'Mean(Value)', 'Mean(c)', 'Mean(k_prime)', 'Mean(kpps_prm)');
            for age = ages_to_check
                val_mean = mean(valS2(:,:,:,age), 'all');
                c_mean = mean(polS2(age).c, 'all');
                k_prime_mean = mean(polS2(age).k_prime, 'all');
                kpps_prime_mean = mean(polS2(age).kpps_prime, 'all');
                fprintf('%-10d %-15.4f %-15.4f %-15.4f %-15.4f\n', age, val_mean, c_mean, k_prime_mean, kpps_prime_mean);
            end

            fprintf('\n请验证测试1和测试2的结果是否一致。\n');
            fprintf('如果一致, 说明在 pps_max=0 时, 两个模型实现了等价。\n\n');
        end


        function run_vfi_egm_test()
            % =========================================================================
            % == 测试函数: run_vfi_egm_test (版本 v3.1 - 深度诊断增强版)
            % == 核心修正:
            % ==   - 增加了 "深度诊断" 环节。在完成整体对比后，如果发现显著差异，
            % ==     此环节将聚焦于误差首次出现的年龄组 (通常是最后一个年龄组 aD_new)。
            % ==   - 对该特定年龄的特定状态点(例如 k 网格的中间点)，分别调用
            % ==     VFI 和 EGM 的核心计算逻辑，并并排打印出价值函数的所有构成
            % ==     部分：c, k', kpps', u(c), u_bequest, discount, V_total。
            % ==   - 这种精细化的输出能够清晰地揭示两种方法在计算价值时的第一个
            % ==     分歧点，从而为最终解决问题提供了决定性的线索。
            % =========================================================================

            fprintf('======================================================\n');
            fprintf('==  开启 household VFI vs. EGM 全面对比测试程序  ==\n');
            fprintf('======================================================\n\n');

            % --- 1. 初始化和参数设置 ---
            cS = utils.ParameterValues();

            cS.nk = 40;
            cS.nkpps = 10;
            cS.nw = 1;
            cS.n_pps_rate_grid = 10;
            cS.pps_max = 0.15;
            cS.pps_active = true;

            [paramS.leGridV, paramS.TrProbM_by_age, paramS.leProb1V, cS.nw_expanded] = utils.EarningProcess_AgeDependent(cS);
            cS.s_pathV = cS.s_pathV(:, 1);
            cS_test = utils.generateGrids(cS);

            M1 = struct('r_mkt_t', 0.04, 'w_t', 1.5, 'tr_per_hh', 0.01, 'b_hat_t', 0.5, 'theta_t', 0.11);

            % --- 2. 使用 VFI (蛮力法) 求解 ---
            fprintf('--- (1/2) 正在使用 VFI (蛮力法) 求解模型 ---\n');
            tic;
            [polS_vfi, valS_vfi] = household.VFI_solver(M1, paramS, cS_test, false);
            vfi_time = toc;
            fprintf('VFI 求解完成。耗时: %.2f 秒\n\n', vfi_time);

            % --- 3. 使用 EGM 求解 ---
            fprintf('--- (2/2) 正在使用 EGM (终点网格法) 求解模型 ---\n');
            tic;
            [polS_egm, valS_egm] = household.VFI_solver(M1, paramS, cS_test, true);
            egm_time = toc;
            fprintf('EGM 求解完成。耗时: %.2f 秒\n\n', egm_time);

            % --- 4. 结果对比分析 ---
            fprintf('--- 开始对比 VFI 和 EGM 的求解结果 ---\n');
            fprintf('求解器性能: EGM 速度 / VFI 速度 = %.2f\n\n', vfi_time / egm_time);

            fprintf('%-5s | %-20s | %-20s | %-20s | %-20s\n', ...
                'Age', 'Diff (Value)', 'Diff (Consumption)', 'Diff (k prime)', 'Diff (kpps prime)');
            fprintf('%-5s | %-10s %-10s | %-10s %-10s | %-10s %-10s | %-10s %-10s\n', ...
                '', 'Max Abs', 'Mean Abs', 'Max Abs', 'Mean Abs', 'Max Abs', 'Mean Abs', 'Max Abs', 'Mean Abs');
            disp(repmat('-', 1, 105));

            max_diff_summary = struct('val', 0, 'c', 0, 'k_prime', 0, 'kpps_prime', 0);

            for a_idx = 1:cS_test.aD_new
                val_vfi_age = valS_vfi(:,:,:,a_idx);
                val_egm_age = valS_egm(:,:,:,a_idx);
                pol_vfi_age = polS_vfi(a_idx);
                pol_egm_age = polS_egm(a_idx);
                diff_val = abs(val_vfi_age - val_egm_age);
                diff_c = abs(pol_vfi_age.c - pol_egm_age.c);
                diff_k_prime = abs(pol_vfi_age.k_prime - pol_egm_age.k_prime);
                diff_kpps_prime = abs(pol_vfi_age.kpps_prime - pol_egm_age.kpps_prime);

                if max(diff_val,[],'all') > max_diff_summary.val, max_diff_summary.val = max(diff_val,[],'all'); end
                if max(diff_c,[],'all') > max_diff_summary.c, max_diff_summary.c = max(diff_c,[],'all'); end
                if max(diff_k_prime,[],'all') > max_diff_summary.k_prime, max_diff_summary.k_prime = max(diff_k_prime,[],'all'); end
                if max(diff_kpps_prime,[],'all') > max_diff_summary.kpps_prime, max_diff_summary.kpps_prime = max(diff_kpps_prime,[],'all'); end

                fprintf('%-5d | %-10.4e %-10.4e | %-10.4e %-10.4e | %-10.4e %-10.4e | %-10.4e %-10.4e\n', ...
                    a_idx, max(diff_val, [], 'all'), mean(diff_val, 'all'), max(diff_c, [], 'all'), mean(diff_c, 'all'), ...
                    max(diff_k_prime, [], 'all'), mean(diff_k_prime, 'all'), max(diff_kpps_prime, [], 'all'), mean(diff_kpps_prime, 'all'));
            end

            disp(repmat('-', 1, 105));
            fprintf('\n测试完成。全局最大差异总结：\n');
            fprintf('  - Value Func (Max|V_vfi - V_egm|):       %.4e\n', max_diff_summary.val);
            fprintf('  - Consumption (Max|c_vfi - c_egm|):     %.4e\n', max_diff_summary.c);
            fprintf('  - Capital Sav. (Max|k_vfi - k_egm|):      %.4e\n', max_diff_summary.k_prime);
            fprintf('  - PPS Sav. (Max|kpps_vfi - kpps_egm|): %.4e\n', max_diff_summary.kpps_prime);

            % --- 5. 深度诊断环节 ---
            if max_diff_summary.val > 1e-3 % 如果存在显著差异
                fprintf('\n--- 深度诊断: 检测到显著差异，聚焦最后一个年龄组 (a=%d) ---\n', cS_test.aD_new);

                a_idx = cS_test.aD_new;
                % 选择一个状态点进行诊断 (例如：无劳动收入冲击，k和kpps网格的中间点)
                ie = 1;
                ik = round(cS_test.nk / 2);
                ikpps = round(cS_test.nkpps / 2);

                fprintf('诊断状态点: a_idx=%d, k_idx=%d (k=%.2f), kpps_idx=%d (kpps=%.2f), ie=%d\n\n',...
                    a_idx, ik, cS_test.kGridV(ik), ikpps, cS_test.kppsGridV(ikpps), ie);

                fprintf('%-10s | %-12s | %-12s | %-12s | %-12s | %-12s | %-12s | %-12s\n', ...
                    'Method', 'c', 'k_prime', 'kpps_prime', 'util(c)', 'util_bequest', 'Discount', 'Total Value');
                disp(repmat('-', 1, 110));

                % VFI 分解
                pol_vfi_pt = polS_vfi(a_idx);
                c_vfi = pol_vfi_pt.c(ik, ikpps, ie);
                k_p_vfi = pol_vfi_pt.k_prime(ik, ikpps, ie);
                kpps_p_vfi = pol_vfi_pt.kpps_prime(ik, ikpps, ie);
                util_c_vfi = (c_vfi.^(1-cS_test.sigma))./(1-cS_test.sigma);
                total_wealth_p_vfi = k_p_vfi + kpps_p_vfi;
                util_b_vfi = cS_test.phi_bequest * (max(cS_test.cFloor, total_wealth_p_vfi).^(1-cS_test.sigma))./(1-cS_test.sigma);
                g_A_period = (1 + cS_test.g_A_ss)^cS_test.time_Step - 1;
                discount = cS_test.beta * (1 + g_A_period)^(1-cS_test.sigma);
                val_vfi_calc = util_c_vfi + discount * util_b_vfi;
                fprintf('%-10s | %-12.4f | %-12.4f | %-12.4f | %-12.4f | %-12.4f | %-12.4f | %-12.4f\n', ...
                    'VFI', c_vfi, k_p_vfi, kpps_p_vfi, util_c_vfi, util_b_vfi, discount, val_vfi_calc);

                % EGM 分解
                pol_egm_pt = polS_egm(a_idx);
                c_egm = pol_egm_pt.c(ik, ikpps, ie);
                k_p_egm = pol_egm_pt.k_prime(ik, ikpps, ie);
                kpps_p_egm = pol_egm_pt.kpps_prime(ik, ikpps, ie);
                util_c_egm = (c_egm.^(1-cS_test.sigma))./(1-cS_test.sigma);
                total_wealth_p_egm = k_p_egm + kpps_p_egm;
                util_b_egm = cS_test.phi_bequest * (max(cS_test.cFloor, total_wealth_p_egm).^(1-cS_test.sigma))./(1-cS_test.sigma);
                val_egm_calc = util_c_egm + discount * util_b_egm;
                fprintf('%-10s | %-12.4f | %-12.4f | %-12.4f | %-12.4f | %-12.4f | %-12.4f | %-12.4f\n', ...
                    'EGM', c_egm, k_p_egm, kpps_p_egm, util_c_egm, util_b_egm, discount, val_egm_calc);

                fprintf('\n诊断说明: 请对比两行数据。Total Value应等于 u(c) + u_bequest * Discount。\n');
                fprintf('如果策略(c, k_prime)相似，但 u(c) 或 u_bequest 不同，则说明效用函数实现有误。\n');
                fprintf('如果所有组件都相似，但从主程序中得到的 valS_vfi 和 valS_egm 不同，则说明数据传递/存储有误。\n');
            end
        end

    end


end



