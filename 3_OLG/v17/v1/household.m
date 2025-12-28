% --- household.m ---
classdef household
    methods (Static)
      function [polS, valS] = VFI_solver(M_vfi, paramS_vfi, cS_vfi)
            % =========================================================================
            % == 函数: VFI_solver
            % == 版本: [v2.0 - PPS 路由版]
            % ==
            % == 核心修改:
            % ==   - 根据全局开关 cS_vfi.pps_active 决定调用哪个版本的年龄求解器。
            % ==   - 初始化 valS 时考虑 kpps 维度。
            % ==   - [注意] 简化的养老金逻辑 bV_payg_vfi 保持不变，实际模型中可能需要
            % ==     与PPS更复杂的互动。
            % =========================================================================
            
            if isfield(cS_vfi, 'pps_active') && cS_vfi.pps_active
                % fprintf('   [家庭求解] 检测到 pps_active = true。使用 VFI_by_age_PPS 求解器。\n');
                % 初始化包含 kpps 维度的价值函数矩阵
                valS = -Inf(cS_vfi.nk, cS_vfi.nkpps, cS_vfi.nw_expanded, cS_vfi.aD_new);
            else
                % fprintf('   [家庭求解] 检测到 pps_active = false。使用 VFI_by_age (无PPS) 求解器。\n');
                % 原有的初始化
                valS = -Inf(cS_vfi.nk, 1, cS_vfi.nw_expanded, cS_vfi.aD_new);
            end

            polS_cell = cell(cS_vfi.aD_new, 1);
            bV_payg_vfi = zeros(1, cS_vfi.aD_new); % 简化养老金逻辑
            
            % --- 反向迭代 ---
            for a_idx = cS_vfi.aD_new : -1 : 1
                vPrime_kkppse_next = [];
                if a_idx < cS_vfi.aD_new
                    vPrime_kkppse_next = valS(:,:,:,a_idx+1);
                end
                
                % --- [核心路由逻辑] ---
                if isfield(cS_vfi, 'pps_active') && cS_vfi.pps_active
                    [val_age, pol_age] = household.VFI_by_age_PPS(a_idx, vPrime_kkppse_next, M_vfi, 0, paramS_vfi, cS_vfi);
                else
                    [val_age, pol_age] = household.VFI_by_age(a_idx, vPrime_kkppse_next, M_vfi, 0, paramS_vfi, cS_vfi);
                end
                
                valS(:,:,:,a_idx) = val_age;
                polS_cell{a_idx} = pol_age;
            end
            polS = [polS_cell{:}];
        end
        
        function [val_age, pol_age] = VFI_by_age_PPS(a_idx, vPrime_kkppse_next, M_age, b_age_val, paramS_age, cS)
            % =========================================================================
            % == 函数: VFI_by_age_PPS
            % == 版本: [v1.0 - 全新实现]
            % ==
            % == 目的:
            % ==   求解带个人养老金账户(PPS)的家庭问题，保留并行计算。
            % =========================================================================

            % --- 1. 初始化 ---
            nk_search_grid = 150; % 储蓄选择的搜索密度

            % 初始化输出变量 (包含 kpps 维度)
            val_age = -1e20 * ones(cS.nk, cS.nkpps, cS.nw_expanded);
            pol_age = struct(...
                'c', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
                'k_prime', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
                'kpps_prime', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
                'tax_regular', zeros(cS.nk, cS.nkpps, cS.nw_expanded));

            % 核心经济常量
            g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
            discount_factor_utility = cS.beta * (1 + g_A_period)^(1 - cS.sigma);
            k_grid_vec = cS.kGridV;
            kpps_grid_vec = cS.kppsGridV; % PPS 资产网格

            % --- 2. 准备下一期值函数的插值对象 ---
            vPrime_interpolants = cell(cS.nw_expanded, 1);
            if a_idx < cS.aD_new
                % 为每个 e' 创建一个 V(k, kpps) 的二维插值
                for ie_next = 1:cS.nw_expanded
                    vPrime_interpolants{ie_next} = griddedInterpolant({k_grid_vec, kpps_grid_vec}, vPrime_kkppse_next(:, :, ie_next), 'linear', 'linear');
                end
            end
            
            % --- 3. [并行化] 遍历所有 (kpps, e) 状态组合 ---
            num_tasks = cS.nkpps * cS.nw_expanded;
            % 预分配cell以收集parfor的结果
            val_results_cell = cell(num_tasks, 1);
            pol_results_cell = cell(num_tasks, 1);

            parfor task_idx = 1:num_tasks
                % 从一维任务索引解码为 (ikpps, ie)
                [ikpps, ie] = ind2sub([cS.nkpps, cS.nw_expanded], task_idx);
                kpps_current = kpps_grid_vec(ikpps);

                % --- 4. 计算预算约束 (对每个 k 状态) ---
                capital_return = k_grid_vec * (1 + M_age.r_mkt_t);
                labor_income_gross = 0;
                pps_contrib = 0;
                pps_withdrawal = 0;
                pps_tax = 0;

                if a_idx <= cS.aR_new % 工作期
                    labor_income_gross = M_age.w_t * cS.ageEffV_new(a_idx) * paramS_age.leGridV(ie);
                    pps_contrib = cS.pps_contrib_rate * labor_income_gross;
                else % 退休期
                    % [注意] 这里的 kpps_current 是去趋势化的值，需要先“复原”再计算提款
                    kpps_total_value_level = kpps_current * (1 + M_age.r_mkt_t) * (1 + g_A_period);
                    pps_withdrawal = cS.pps_withdrawal_rate * kpps_total_value_level;
                    pps_tax = cS.pps_tax_rate_withdrawal * pps_withdrawal;
                end
                
                capital_tax_vec = cS.tau_k .* (k_grid_vec * M_age.r_mkt_t);
                labor_tax = cS.tau_l * labor_income_gross;
                
                % [BGP对齐] cash_on_hand是水平(level)值
                cash_on_hand_vec = capital_return + labor_income_gross + M_age.tr_per_hh ...
                    - (capital_tax_vec + labor_tax + pps_contrib - (pps_withdrawal - pps_tax));
                
                % --- 5. 计算下一期PPS资产 (确定性) ---
                % [BGP对齐] 下一期PPS资产必须去趋势化
                kpps_prime_level = (kpps_current * (1 + M_age.r_mkt_t)) + (pps_contrib - pps_withdrawal);
                kpps_prime_detrended = max(cS.kppsMin, kpps_prime_level / (1 + g_A_period));
                
                % --- 6. 求解最优 k' 和 c ---
                k_prime_max_vec = (cash_on_hand_vec - cS.cFloor * (1 + cS.tau_c)) / (1 + g_A_period);
                k_prime_max_vec(k_prime_max_vec < cS.kMin) = cS.kMin;
                
                k_prime_choices_mat = cS.kMin + (k_prime_max_vec - cS.kMin) .* linspace(0, 1, nk_search_grid);
                
                c_expend_choices_mat = cash_on_hand_vec - k_prime_choices_mat * (1 + g_A_period);
                c_choices_mat = c_expend_choices_mat / (1 + cS.tau_c);

                invalid_choice_mask = (c_choices_mat < cS.cFloor);
                c_choices_mat(invalid_choice_mask) = cS.cFloor;
                [~, util_c_mat] = utils.CES_utility(c_choices_mat, cS.sigma, cS);
                util_c_mat(invalid_choice_mask) = -1e20;

                future_value_mat = zeros(cS.nk, nk_search_grid);
                if a_idx < cS.aD_new
                    trans_mat_row = paramS_age.TrProbM_by_age{a_idx + 1}(ie, :);
                    
                    % 准备插值点 (k', kpps')
                    k_prime_choices_flat = k_prime_choices_mat(:);
                    kpps_prime_vec = ones(length(k_prime_choices_flat), 1) * kpps_prime_detrended;
                    
                    % 使用二维插值获得期望值 V(k', kpps' | e)
                    v_prime_interp_mat = zeros(length(k_prime_choices_flat), cS.nw_expanded);
                    for ie_next = 1:cS.nw_expanded
                        v_prime_interp_mat(:, ie_next) = vPrime_interpolants{ie_next}(k_prime_choices_flat, kpps_prime_vec);
                    end
                    ev_flat_vec = v_prime_interp_mat * trans_mat_row';
                    ev_on_choices_mat = reshape(ev_flat_vec, cS.nk, nk_search_grid);
                    
                    % [注意] 遗赠动机作用于总财富 (k' + kpps')
                    total_wealth_prime_mat = k_prime_choices_mat + kpps_prime_detrended;
                    util_bequest_mat = utils.bequest_utility(total_wealth_prime_mat, cS);

                    future_value_mat = cS.s_pathV(a_idx) * ev_on_choices_mat + (1 - cS.s_pathV(a_idx)) * util_bequest_mat;
                else % 最后一期
                    total_wealth_prime_mat = k_prime_choices_mat + kpps_prime_detrended;
                    future_value_mat = utils.bequest_utility(total_wealth_prime_mat, cS);
                end

                total_value_mat = util_c_mat + discount_factor_utility * future_value_mat;

                % --- 7. 存储此 (ikpps, ie) 状态下的最优决策 ---
                [best_val_vec, best_idx_vec] = max(total_value_mat, [], 2);
                linear_indices = sub2ind(size(k_prime_choices_mat), (1:cS.nk)', best_idx_vec);
                
                % 将结果存入cell
                val_results_cell{task_idx} = best_val_vec;
                pol_results_cell{task_idx} = struct(...
                    'c', c_choices_mat(linear_indices), ...
                    'k_prime', k_prime_choices_mat(linear_indices), ...
                    'kpps_prime', repmat(kpps_prime_detrended, cS.nk, 1), ... % kpps'不依赖于k
                    'tax_regular', capital_tax_vec + labor_tax + c_choices_mat(linear_indices) * cS.tau_c + pps_tax);
            end

            % --- 8. 将并行计算的结果安全地组装回最终矩阵 ---
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
            % 单个年龄的VFI求解
            nk_search_grid = 150;
            val_age = -1e20 * ones(cS.nk, 1, cS.nw_expanded);
            pol_age = struct('c', zeros(cS.nk, 1, cS.nw_expanded), ...
                'k_prime', zeros(cS.nk, 1, cS.nw_expanded), ...
                'tax_regular', zeros(cS.nk, 1, cS.nw_expanded)); 

            g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
            discount_factor_utility = cS.beta * (1 + g_A_period)^(1 - cS.sigma);
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
                if a_idx <= cS.aR_new
                    labor_income_gross = M_age.w_t * cS.ageEffV_new(a_idx) * paramS_age.leGridV(ie);
                end
                
                % 简化税收：只考虑资本税和劳动税
                capital_tax_vec = cS.tau_k .* (k_grid_vec * M_age.r_mkt_t);
                labor_tax = cS.tau_l * labor_income_gross;
                
                % 预算约束
                cash_on_hand_vec = capital_return + labor_income_gross + M_age.tr_per_hh - (capital_tax_vec + labor_tax);
                available_for_c_and_s_vec = cash_on_hand_vec; % 假设无冲击支出

                k_prime_max_vec = (available_for_c_and_s_vec - cS.cFloor * (1 + cS.tau_c)) / (1 + g_A_period);
                k_prime_max_vec(k_prime_max_vec < cS.kMin) = cS.kMin;
                
                k_prime_choices_mat = cS.kMin + (k_prime_max_vec - cS.kMin) .* linspace(0, 1, nk_search_grid);
                c_expend_choices_mat = available_for_c_and_s_vec - k_prime_choices_mat * (1 + g_A_period);
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
                    future_value_mat = cS.s_pathV(a_idx) * ev_on_choices_mat + (1 - cS.s_pathV(a_idx)) * utils.bequest_utility(k_prime_choices_mat, cS);
                else
                    future_value_mat = utils.bequest_utility(k_prime_choices_mat, cS);
                end
                
                total_value_mat = util_c_mat + discount_factor_utility * future_value_mat;
                [best_val_vec, best_idx_vec] = max(total_value_mat, [], 2);
                linear_indices = sub2ind(size(k_prime_choices_mat), (1:cS.nk)', best_idx_vec);
                
                val_age(:, 1, ie) = best_val_vec;
                pol_age.c(:, 1, ie) = c_choices_mat(linear_indices);
                pol_age.k_prime(:, 1, ie) = k_prime_choices_mat(linear_indices);
                pol_age.tax_regular(:, 1, ie) = capital_tax_vec + labor_tax + c_choices_mat(linear_indices) .* cS.tau_c;
            end
        end
    end
end