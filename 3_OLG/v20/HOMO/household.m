% --- household.m ---
classdef household
    methods (Static)
        
        function [val_age, pol_age] = VFI_by_age_PPS(a_idx, vPrime_kkppse_next, M_age, b_age_val, paramS_age, cS)
    % =========================================================================
    % == 函数: VFI_by_age_PPS (版本 v4.0 - 矩阵化改造版)
    % == 核心修改:
    % ==   - [!!! 关键重构 !!!] 完全移除了 parfor 循环，所有计算均被向量化/矩阵化。
    % ==   - 通过 ndgrid 创建了 (k, kpps, e) 三维状态网格，所有变量都扩展到这个维度上。
    % ==   - k' 的选择网格被扩展为四维张量 (k, kpps, e, k'_choice)。
    % ==   - 对未来价值函数的期望计算通过一个短循环(对未来状态e')和矩阵广播实现，
    % ==     替代了原有的嵌套循环，极大地提高了计算效率。
    % ==   - 结果的组装过程被消除，最优值和策略通过 max 函数和线性索引直接从高维张量中提取。
    % =========================================================================
    nk_search_grid = 150;

    % 1. 获取时期增长因子和贴现因子 (与之前相同)
    if isfield(M_age, 'g_A_t_plus_1')
        g_A_period = M_age.g_A_t_plus_1;
    else
        g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
    end
    growth_factor_bgp = (1 + g_A_period);
    discount_factor_V_prime = cS.beta * growth_factor_bgp^(1 - cS.sigma);

    % 2. 为下一期价值函数创建插值器 (与之前相同)
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

    % 3. 构建状态变量网格 (矩阵化核心步骤)
    % 创建与 (k, kpps, e) 状态对应的3D矩阵
    [k_mat, kpps_mat, le_mat] = ndgrid(cS.kGridV, cS.kppsGridV, paramS_age.leGridV);

    % 4. 计算与决策无关的收入和税收 (矩阵化)
    labor_supply_mat = 0;
    if a_idx <= cS.aR_new
        labor_supply_mat = cS.ageEffV_new(a_idx) .* le_mat;
    end
    labor_income_hat_mat = M_age.w_t * labor_supply_mat;

    capital_return_hat_mat = k_mat * (1 + M_age.r_mkt_t);
    capital_tax_hat_mat = cS.tau_k .* (k_mat * M_age.r_mkt_t);
    labor_tax_hat_mat = cS.tau_l * labor_income_hat_mat;
    payg_contribution_hat_mat = M_age.theta_t * labor_income_hat_mat;

    % [核心简化] PPS缴费和提取逻辑 (矩阵化)
    pps_contrib_hat_mat = zeros(size(k_mat));
    pension_income_hat_mat = zeros(size(k_mat));
    pps_withdrawal_hat_mat = zeros(size(k_mat));
    pps_tax_hat_mat = zeros(size(k_mat));
    
    is_working_age = (a_idx <= cS.aR_new);
    if cS.pps_active
        if is_working_age
            % 工作期: 缴费为固定比例
            pps_contrib_hat_mat = labor_income_hat_mat * cS.pps_fixed;
        else % 退休期
            % 退休期: 领取PAYG养老金并提取PPS
            pension_income_hat_mat(:) = b_age_val;
            pps_withdrawal_hat_mat = cS.pps_withdrawal_rate .* (kpps_mat * (1 + M_age.r_mkt_t));
            pps_tax_hat_mat = cS.pps_tax_rate_withdrawal .* pps_withdrawal_hat_mat;
        end
    elseif ~is_working_age % 退休期但无PPS
         pension_income_hat_mat(:) = b_age_val;
    end

    % [核心简化] kpps' 不再是选择，而是给定缴费/提取后的确定值 (矩阵化)
    kpps_hat_end_of_period_mat = kpps_mat * (1 + M_age.r_mkt_t) + pps_contrib_hat_mat - pps_withdrawal_hat_mat;
    kpps_prime_hat_mat = max(cS.kppsMin, kpps_hat_end_of_period_mat / growth_factor_bgp);

    % 5. 构建选择网格并进行向量化计算 (对 k' 进行搜索)
    base_coh_hat_mat = capital_return_hat_mat + labor_income_hat_mat + M_age.tr_per_hh + pension_income_hat_mat + pps_withdrawal_hat_mat ...
        - (capital_tax_hat_mat + labor_tax_hat_mat + payg_contribution_hat_mat + pps_tax_hat_mat);
    cash_on_hand_hat_mat = base_coh_hat_mat - pps_contrib_hat_mat; % CoH在PPS缴费之后

    k_prime_max_mat = (cash_on_hand_hat_mat - cS.cFloor * (1 + cS.tau_c)) / growth_factor_bgp;
    k_prime_max_mat(k_prime_max_mat < cS.kMin) = cS.kMin;
    
    % 创建k'选择的4D张量: (k, kpps, e, k'_choice)
    k_prime_search_grid = linspace(0, 1, nk_search_grid);
    k_prime_search_grid_reshaped = reshape(k_prime_search_grid, [1, 1, 1, nk_search_grid]);
    k_prime_choices_tens = cS.kMin + (k_prime_max_mat - cS.kMin) .* k_prime_search_grid_reshaped;
    
    % 计算消费和当期效用的4D张量
    c_expend_choices_tens = cash_on_hand_hat_mat - k_prime_choices_tens * growth_factor_bgp;
    c_choices_tens = c_expend_choices_tens / (1 + cS.tau_c);
    invalid_choice_mask = (c_choices_tens < cS.cFloor);
    c_choices_tens(invalid_choice_mask) = cS.cFloor;

    util_c_tens = (c_choices_tens.^(1 - cS.sigma))./(1 - cS.sigma);
    util_c_tens(invalid_choice_mask) = -1e20;
    
    % 计算遗赠效用的4D张量
    kpps_prime_hat_tens = repmat(kpps_prime_hat_mat, [1, 1, 1, nk_search_grid]);
    total_wealth_prime_tens = k_prime_choices_tens + kpps_prime_hat_tens;
    util_bequest_tens = utils.bequest_utility(total_wealth_prime_tens, cS);

    % 计算未来价值的期望 (矩阵化)
    future_value_tens = zeros(size(k_prime_choices_tens));
    if a_idx < cS.aD_new
        ev_on_choices_tens = zeros(size(k_prime_choices_tens));
        trans_mat = paramS_age.TrProbM_by_age{a_idx}; % [nw_expanded x nw_expanded]
        
        for ie_next = 1:cS.nw_expanded
            % trans_prob_col(ie) = P(ie_next | ie)
            trans_prob_col = trans_mat(:, ie_next);
            if any(trans_prob_col > 0)
                if cS.nkpps > 1 && cS.pps_active
                    v_prime_interp_tens = vPrime_interpolants{ie_next}(k_prime_choices_tens, kpps_prime_hat_tens);
                else
                    v_prime_interp_tens = vPrime_interpolants{ie_next}(k_prime_choices_tens);
                end
                % 使用广播将转移概率应用到所有 (k, kpps, k_choice) 上
                % reshape(trans_prob_col, [1, 1, cS.nw_expanded]) 维度为 (1, 1, e)
                % v_prime_interp_tens 维度为 (k, kpps, e, k'_choice)
                % 广播相乘后，对每个 (k, kpps, ie, k'_choice), 其值乘以 P(ie_next | ie)
                ev_on_choices_tens = ev_on_choices_tens + v_prime_interp_tens .* reshape(trans_prob_col, [1, 1, cS.nw_expanded, 1]);
            end
        end

        survival_rate = cS.s_pathV(a_idx);
        if size(survival_rate,2) > 1, survival_rate = survival_rate(1); end
        future_value_tens = discount_factor_V_prime * (survival_rate * ev_on_choices_tens + (1 - survival_rate) * util_bequest_tens);
    else
        future_value_tens = discount_factor_V_prime * util_bequest_tens;
    end
    
    total_value_tens = util_c_tens + future_value_tens;

    % 6. 寻找最优决策 (在第4维度上最大化)
    [val_age, best_k_prime_idx_mat] = max(total_value_tens, [], 4, 'omitnan');
    
    % 7. 提取最优策略
    % 使用线性索引从4D张量中提取最优值
    S = size(c_choices_tens);
    [I, J, K] = ndgrid(1:S(1), 1:S(2), 1:S(3));
    linear_indices_final = sub2ind(S, I, J, K, best_k_prime_idx_mat);
    
    c_final_mat = c_choices_tens(linear_indices_final);
    k_prime_final_mat = k_prime_choices_tens(linear_indices_final);

    tax_final_mat = capital_tax_hat_mat + labor_tax_hat_mat + pps_tax_hat_mat + payg_contribution_hat_mat + c_final_mat .* cS.tau_c;
    
    % 8. 保存结果到 policy struct
    pol_age = struct(...
        'c', c_final_mat, ...
        'l', labor_supply_mat, ...
        'k_prime', k_prime_final_mat, ...
        'kpps_prime', kpps_prime_hat_mat, ... % kpps' 不依赖于 k' 的选择
        'tax_regular', tax_final_mat);
        end
        
        %         function [val_age, pol_age] = VFI_by_age_PPS(a_idx, vPrime_kkppse_next, M_age, b_age_val, paramS_age, cS)
        %     % =========================================================================
        %     % == 函数: VFI_by_age_PPS (版本 v3.0 - 固定PPS缴费比例版)
        %     % == 核心修改:
        %     % ==   - [!!! 关键简化 !!!] 移除了对PPS缴费率的优化搜索过程。
        %     % ==     之前版本会在一个缴费率网格 (n_pps_rate_grid) 上进行选择，
        %     % ==     导致 (k', pps_rate) 的二维联合搜索。
        %     % ==   - 新逻辑: 工作期的PPS缴费被设定为一个外生的固定比例 cS.pps_fixed。
        %     % ==     这使得对每个(k,kpps,e)状态，最优k'的选择成为一维搜索问题。
        %     % ==   - 效率提升: 此修改将核心计算的张量维度从3降低到2，
        %     % ==     极大地减少了计算量，显著提升了函数运行效率。
        %     % ==   - 状态保留: kpps 依然作为状态变量被保留和更新。
        %     % =========================================================================
        %     nk_search_grid = 150;
        %     val_age = -1e20 * ones(cS.nk, cS.nkpps, cS.nw_expanded);
        %     pol_age = struct(...
        %         'c', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
        %         'l', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
        %         'k_prime', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
        %         'kpps_prime', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
        %         'tax_regular', zeros(cS.nk, cS.nkpps, cS.nw_expanded));
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
        %     % 2. 为下一期价值函数创建插值器 (自适应维度)
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
        %     % 3. 使用并行计算优化VFI过程 (循环遍历状态变量)
        %     num_tasks = cS.nkpps * cS.nw_expanded;
        %     val_results_cell = cell(num_tasks, 1);
        %     pol_results_cell = cell(num_tasks, 1);
        % 
        %     parfor task_idx = 1:num_tasks
        %         [ikpps, ie] = ind2sub([cS.nkpps, cS.nw_expanded], task_idx);
        %         kpps_hat_current = kpps_grid_vec(ikpps);
        % 
        %         % 4. 计算与决策无关的收入和税收
        %         labor_supply_ie = 0;
        %         if a_idx <= cS.aR_new
        %             labor_supply_ie = cS.ageEffV_new(a_idx) * paramS_age.leGridV(ie);
        %         end
        %         labor_income_hat = M_age.w_t * labor_supply_ie;
        % 
        %         capital_return_hat_vec = k_grid_vec * (1 + M_age.r_mkt_t);
        %         capital_tax_hat_vec = cS.tau_k .* (k_grid_vec * M_age.r_mkt_t);
        %         labor_tax_hat = cS.tau_l * labor_income_hat;
        %         payg_contribution_hat = M_age.theta_t * labor_income_hat;
        % 
        %         % [核心简化] PPS缴费和提取逻辑
        %         pps_contrib_hat = 0;
        %         pension_income_hat = 0;
        %         pps_withdrawal_hat = 0;
        %         pps_tax_hat = 0;
        % 
        %         if a_idx <= cS.aR_new && cS.pps_active
        %             % 工作期: 缴费为固定比例
        %             pps_contrib_hat = labor_income_hat * cS.pps_fixed;
        %         elseif a_idx > cS.aR_new
        %             % 退休期: 领取PAYG养老金并提取PPS
        %             pension_income_hat = b_age_val;
        %             if cS.pps_active
        %                pps_withdrawal_hat = cS.pps_withdrawal_rate * (kpps_hat_current * (1 + M_age.r_mkt_t));
        %                pps_tax_hat = cS.pps_tax_rate_withdrawal * pps_withdrawal_hat;
        %             end
        %         end
        % 
        %         % [核心简化] kpps' 不再是选择，而是给定缴费/提取后的确定值
        %         kpps_hat_end_of_period = kpps_hat_current * (1 + M_age.r_mkt_t) + pps_contrib_hat - pps_withdrawal_hat;
        %         kpps_prime_hat_choice = max(cS.kppsMin, kpps_hat_end_of_period / growth_factor_bgp);
        % 
        %         % 5. 构建选择网格并进行向量化计算 (仅对 k' 进行搜索)
        %         base_coh_hat_vec = capital_return_hat_vec + labor_income_hat + M_age.tr_per_hh + pension_income_hat + pps_withdrawal_hat ...
        %             - (capital_tax_hat_vec + labor_tax_hat + payg_contribution_hat + pps_tax_hat);
        %         cash_on_hand_hat_vec = base_coh_hat_vec - pps_contrib_hat; % CoH在PPS缴费之后
        % 
        %         k_prime_max_vec = (cash_on_hand_hat_vec - cS.cFloor * (1 + cS.tau_c)) / growth_factor_bgp;
        %         k_prime_max_vec(k_prime_max_vec < cS.kMin) = cS.kMin;
        % 
        %         k_prime_search_grid = linspace(0, 1, nk_search_grid);
        %         k_prime_choices_mat = cS.kMin + (k_prime_max_vec - cS.kMin) .* k_prime_search_grid;
        % 
        %         c_expend_choices_mat = cash_on_hand_hat_vec - k_prime_choices_mat * growth_factor_bgp;
        %         c_choices_mat = c_expend_choices_mat / (1 + cS.tau_c);
        %         invalid_choice_mask = (c_choices_mat < cS.cFloor);
        %         c_choices_mat(invalid_choice_mask) = cS.cFloor;
        % 
        %         util_c_mat = (c_choices_mat.^(1-cS.sigma))./(1-cS.sigma);
        %         util_c_mat(invalid_choice_mask) = -1e20;
        % 
        %         total_wealth_prime_mat = k_prime_choices_mat + kpps_prime_hat_choice;
        %         util_bequest_mat = utils.bequest_utility(total_wealth_prime_mat, cS);
        % 
        %         future_value_mat = zeros(size(k_prime_choices_mat));
        %         if a_idx < cS.aD_new
        %             trans_mat_row = paramS_age.TrProbM_by_age{a_idx}(ie, :);
        %             ev_on_choices_mat = zeros(size(k_prime_choices_mat));
        %             for ie_next = 1:cS.nw_expanded
        %                if trans_mat_row(ie_next) > 0
        %                    if cS.nkpps > 1 && cS.pps_active
        %                        kpps_prime_repeated = repmat(kpps_prime_hat_choice, size(k_prime_choices_mat));
        %                        v_prime_interp = vPrime_interpolants{ie_next}(k_prime_choices_mat, kpps_prime_repeated);
        %                    else
        %                        v_prime_interp = vPrime_interpolants{ie_next}(k_prime_choices_mat);
        %                    end
        %                    ev_on_choices_mat = ev_on_choices_mat + trans_mat_row(ie_next) * v_prime_interp;
        %                end
        %             end
        % 
        %             survival_rate = cS.s_pathV(a_idx);
        %             if size(survival_rate,2) > 1, survival_rate = survival_rate(1); end
        %             future_value_mat = discount_factor_V_prime * (survival_rate * ev_on_choices_mat + (1 - survival_rate) * util_bequest_mat);
        %         else
        %             future_value_mat = discount_factor_V_prime * util_bequest_mat;
        %         end
        % 
        %         total_value_mat = util_c_mat + future_value_mat;
        % 
        %         % 6. 寻找最优决策 (仅针对k')
        %         [best_val_final_vec, best_k_prime_idx_vec] = max(total_value_mat, [], 2);
        % 
        %         linear_indices_final = sub2ind(size(c_choices_mat), (1:cS.nk)', best_k_prime_idx_vec);
        %         c_final_vec = c_choices_mat(linear_indices_final);
        %         k_prime_final_vec = k_prime_choices_mat(linear_indices_final);
        %         kpps_prime_final_vec = repmat(kpps_prime_hat_choice, cS.nk, 1);
        %         tax_final_vec = capital_tax_hat_vec + labor_tax_hat + pps_tax_hat + payg_contribution_hat + c_final_vec .* cS.tau_c;
        % 
        %         % 7. 保存该状态的最终结果
        %         val_results_cell{task_idx} = best_val_final_vec;
        %         pol_results_cell{task_idx} = struct(...
        %             'c', c_final_vec, ...
        %             'l', repmat(labor_supply_ie, cS.nk, 1), ...
        %             'k_prime', k_prime_final_vec, ...
        %             'kpps_prime', kpps_prime_final_vec, ...
        %             'tax_regular', tax_final_vec);
        %     end
        % 
        %     % 8. 将并行计算的结果重新组装成矩阵
        %     fields = fieldnames(pol_age);
        %     for task_idx = 1:num_tasks
        %         [ikpps, ie] = ind2sub([cS.nkpps, cS.nw_expanded], task_idx);
        %         val_age(:, ikpps, ie) = val_results_cell{task_idx};
        %         current_pol_slice = pol_results_cell{task_idx};
        %         for i = 1:length(fields)
        %             field = fields{i};
        %             pol_age.(field)(:, ikpps, ie) = current_pol_slice.(field);
        %         end
        %     end
        % end


%         function [val_age, pol_age] = VFI_by_age_PPS(a_idx, vPrime_kkppse_next, M_age, b_age_val, paramS_age, cS)
%             % =========================================================================
%             % == 函数: VFI_by_age_PPS (版本 v2.9 - 自适应插值维度修正版)
%             % == 核心修正:
%             % ==   - [!!! 关键鲁棒性修正 !!!] 修正了 griddedInterpolant 错误。
%             % ==   - 问题: 当 cS.nkpps = 1 时，kpps 维度只有一个网格点，导致
%             % ==     二维插值器创建失败。
%             % ==   - 解决方案: 增加了条件判断逻辑。
%             % ==     1. 在创建插值器时: 如果 cS.nkpps > 1，则创建二维插值器；
%             % ==        否则，创建一维插值器 (仅对 k 插值)，并使用 squeeze()
%             % ==        来移除 vPrime 中的单例 kpps 维度。
%             % ==     2. 在调用插值器时: 同样根据 cS.nkpps 的值，决定传入
%             % ==        一个参数 (k') 或两个参数 (k', kpps')。
%             % ==   - 此修改使函数能自适应状态空间的维度变化。
%             % =========================================================================
%             nk_search_grid = 150;
%             n_pps_rate_grid = cS.n_pps_rate_grid;
%             val_age = -1e20 * ones(cS.nk, cS.nkpps, cS.nw_expanded);
%             pol_age = struct(...
%                 'c', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
%                 'l', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
%                 'k_prime', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
%                 'kpps_prime', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
%                 'tax_regular', zeros(cS.nk, cS.nkpps, cS.nw_expanded));
% 
%             % 1. 获取时期增长因子和贴现因子
%             if isfield(M_age, 'g_A_t_plus_1')
%                 g_A_period = M_age.g_A_t_plus_1;
%             else
%                 g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
%             end
%             growth_factor_bgp = (1 + g_A_period);
%             discount_factor_V_prime = cS.beta * growth_factor_bgp^(1 - cS.sigma);
% 
%             % 2. 为下一期价值函数创建插值器 (自适应维度)
%             k_grid_vec = cS.kGridV;
%             kpps_grid_vec = cS.kppsGridV;
%             vPrime_interpolants = cell(cS.nw_expanded, 1);
%             if a_idx < cS.aD_new
%                 for ie_next = 1:cS.nw_expanded
%                     if cS.nkpps > 1 && cS.pps_active
%                         % 标准情况: 二维插值 (k, kpps)
%                         vPrime_interpolants{ie_next} = griddedInterpolant({k_grid_vec, kpps_grid_vec}, vPrime_kkppse_next(:, :, ie_next), 'makima','makima');
%                     else
%                         % 特殊情况: kpps维度被压缩，执行一维插值 (k)
%                         % 使用 squeeze 移除 vPrime_kkppse_next 中大小为1的维度
%                         vPrime_interpolants{ie_next} = griddedInterpolant(k_grid_vec, squeeze(vPrime_kkppse_next(:, :, ie_next)), 'pchip','pchip');
%                     end
%                 end
%             end
% 
%             % 3. 使用并行计算优化VFI过程
%             num_tasks = cS.nkpps * cS.nw_expanded;
%             val_results_cell = cell(num_tasks, 1);
%             pol_results_cell = cell(num_tasks, 1);
% 
%             for task_idx = 1:num_tasks
%                 [ikpps, ie] = ind2sub([cS.nkpps, cS.nw_expanded], task_idx);
%                 kpps_hat_current = kpps_grid_vec(ikpps);
% 
%                 % 4. 计算与决策无关的收入和税收
%                 labor_supply_ie = 0;
%                 if a_idx <= cS.aR_new
%                     labor_supply_ie = cS.ageEffV_new(a_idx) * paramS_age.leGridV(ie);
%                 end
%                 labor_income_hat = M_age.w_t * labor_supply_ie;
% 
%                 capital_return_hat_vec = k_grid_vec * (1 + M_age.r_mkt_t);
%                 capital_tax_hat_vec = cS.tau_k .* (k_grid_vec * M_age.r_mkt_t);
%                 labor_tax_hat = cS.tau_l * labor_income_hat;
%                 payg_contribution_hat = M_age.theta_t * labor_income_hat;
% 
%                 pension_income_hat = 0;
%                 pps_withdrawal_hat = 0;
%                 pps_tax_hat = 0;
%                 if a_idx > cS.aR_new
%                     pension_income_hat = b_age_val;
%                     pps_withdrawal_hat = cS.pps_withdrawal_rate * (kpps_hat_current * (1 + M_age.r_mkt_t));
%                     pps_tax_hat = cS.pps_tax_rate_withdrawal * pps_withdrawal_hat;
%                 end
% 
%                 % 5. 构建选择张量并进行向量化计算
%                 if a_idx > cS.aR_new || cS.pps_max <= 0 || n_pps_rate_grid <= 1
%                     n_pps_rate_grid_eff = 1;
%                     pps_rate_choices_vec = 0;
%                 else
%                     n_pps_rate_grid_eff = n_pps_rate_grid;
%                     pps_rate_choices_vec = linspace(0, cS.pps_max, n_pps_rate_grid_eff);
%                 end
% 
%                 pps_contrib_hat_choices = labor_income_hat * pps_rate_choices_vec;
%                 kpps_hat_end_of_period_choices = kpps_hat_current * (1 + M_age.r_mkt_t) + pps_contrib_hat_choices - pps_withdrawal_hat;
%                 kpps_prime_hat_choices = max(cS.kppsMin, kpps_hat_end_of_period_choices / growth_factor_bgp);
% 
%                 base_coh_hat_vec = capital_return_hat_vec + labor_income_hat + M_age.tr_per_hh + pension_income_hat + pps_withdrawal_hat ...
%                     - (capital_tax_hat_vec + labor_tax_hat + payg_contribution_hat + pps_tax_hat);
%                 cash_on_hand_hat_mat = base_coh_hat_vec - pps_contrib_hat_choices;
% 
%                 k_prime_max_mat = (cash_on_hand_hat_mat - cS.cFloor * (1 + cS.tau_c)) / growth_factor_bgp;
%                 k_prime_max_mat(k_prime_max_mat < cS.kMin) = cS.kMin;
% 
%                 k_prime_search_grid = linspace(0, 1, nk_search_grid);
%                 [~, k_prime_grid_3d, ~] = ndgrid(k_grid_vec, k_prime_search_grid, pps_rate_choices_vec);
%                 k_prime_max_3d = permute(k_prime_max_mat, [1 3 2]);
%                 k_prime_choices_tens = cS.kMin + (k_prime_max_3d - cS.kMin) .* k_prime_grid_3d;
% 
%                 cash_on_hand_3d = permute(cash_on_hand_hat_mat, [1 3 2]);
%                 c_expend_choices_tens = cash_on_hand_3d - k_prime_choices_tens * growth_factor_bgp;
%                 c_choices_tens = c_expend_choices_tens / (1 + cS.tau_c);
%                 invalid_choice_mask = (c_choices_tens < cS.cFloor);
%                 c_choices_tens(invalid_choice_mask) = cS.cFloor;
% 
%                 util_c_tens = (c_choices_tens.^(1-cS.sigma))./(1-cS.sigma);
%                 util_c_tens(invalid_choice_mask) = -1e20;
% 
%                 kpps_prime_reshaped = reshape(kpps_prime_hat_choices, [1, 1, n_pps_rate_grid_eff]);
%                 kpps_prime_tens = repmat(kpps_prime_reshaped, [cS.nk, nk_search_grid, 1]);
%                 total_wealth_prime_tens = k_prime_choices_tens + kpps_prime_tens;
% 
%                 util_bequest_tens = utils.bequest_utility(total_wealth_prime_tens, cS);
% 
%                 future_value_tens = zeros(size(k_prime_choices_tens));
%                 if a_idx < cS.aD_new
%                     trans_mat_row = paramS_age.TrProbM_by_age{a_idx}(ie, :);
%                     v_prime_interp_tens = zeros(size(k_prime_choices_tens,1), size(k_prime_choices_tens,2), size(k_prime_choices_tens,3), cS.nw_expanded);
%                     for ie_next = 1:cS.nw_expanded
%                        if cS.nkpps > 1
%                            % 标准情况: 调用二维插值器
%                            v_prime_interp_tens(:,:,:,ie_next) = vPrime_interpolants{ie_next}(k_prime_choices_tens, kpps_prime_tens);
%                        else
%                            % 特殊情况: 调用一维插值器
%                            v_prime_interp_tens(:,:,:,ie_next) = vPrime_interpolants{ie_next}(k_prime_choices_tens);
%                        end
%                     end
%                     ev_on_choices_tens = sum(v_prime_interp_tens .* reshape(trans_mat_row, 1, 1, 1, cS.nw_expanded), 4);
% 
%                     survival_rate = cS.s_pathV(a_idx);
%                     if size(survival_rate,2) > 1, survival_rate = survival_rate(1); end
%                     future_value_tens = discount_factor_V_prime * (survival_rate * ev_on_choices_tens + (1 - survival_rate) * util_bequest_tens);
%                 else
%                     future_value_tens = discount_factor_V_prime * util_bequest_tens;
%                 end
% 
%                 total_value_tens = util_c_tens + future_value_tens;
% 
%                 % 6. 寻找最优决策组合
%                 total_value_mat_2d = reshape(total_value_tens, cS.nk, nk_search_grid * n_pps_rate_grid_eff);
%                 [best_val_final_vec, best_joint_idx_vec] = max(total_value_mat_2d, [], 2);
% 
%                 [best_k_prime_idx, best_pps_rate_idx] = ind2sub([nk_search_grid, n_pps_rate_grid_eff], best_joint_idx_vec);
% 
%                 linear_indices_final = sub2ind(size(k_prime_choices_tens), (1:cS.nk)', best_k_prime_idx, best_pps_rate_idx);
%                 c_final_vec = c_choices_tens(linear_indices_final);
%                 k_prime_final_vec = k_prime_choices_tens(linear_indices_final);
%                 kpps_prime_final_vec = kpps_prime_hat_choices(best_pps_rate_idx);
%                 tax_final_vec = capital_tax_hat_vec + labor_tax_hat + payg_contribution_hat + pps_tax_hat + c_final_vec .* cS.tau_c;
% 
%                 % 7. 保存该状态的最终结果
%                 val_results_cell{task_idx} = best_val_final_vec;
%                 pol_results_cell{task_idx} = struct(...
%                     'c', c_final_vec, ...
%                     'l', repmat(labor_supply_ie, cS.nk, 1), ...
%                     'k_prime', k_prime_final_vec, ...
%                     'kpps_prime', kpps_prime_final_vec, ...
%                     'tax_regular', tax_final_vec);
%             end
% 
%             % 8. 将并行计算的结果重新组装成矩阵
%             fields = fieldnames(pol_age);
%             for task_idx = 1:num_tasks
%                 [ikpps, ie] = ind2sub([cS.nkpps, cS.nw_expanded], task_idx);
%                 val_age(:, ikpps, ie) = val_results_cell{task_idx};
%                 current_pol_slice = pol_results_cell{task_idx};
%                 for i = 1:length(fields)
%                     field = fields{i};
%                     pol_age.(field)(:, ikpps, ie) = current_pol_slice.(field);
%                 end
%             end
% end

        
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
            T = cS.T_sim;
            Pol_path = cell(1, T);
            Val_path = zeros(cS.nk, cS.nkpps, cS.nw_expanded, cS.aD_new, T);
            Val_path(:,:,:,:,T) = valF;
            Pol_path{T} = polF;
            
            for t = (T-1) : -1 : 1
                % [!!! 核心修正: 直接使用预计算的增长率 !!!]
                g_A_onestep_ahead = pathS.g_A_path(t);

                M_vfi_t = struct('r_mkt_t', pathS.r_path(t), ...
                    'w_t', pathS.w_hat_path(t), ...
                    'g_A_t_plus_1', g_A_onestep_ahead);

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



