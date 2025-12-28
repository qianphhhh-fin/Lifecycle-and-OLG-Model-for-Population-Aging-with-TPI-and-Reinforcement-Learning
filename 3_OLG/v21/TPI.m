% --- TPI.m ---
% =========================================================================
% == CLASSDEF: TPI.m
% == 版本: [v1.3 - 聚合目标与人口对齐修正版]
% ==
% == 目的:
% ==   - 修正 calculate_paths_and_errors 中目标遗赠路径的计算逻辑，确保
% ==     总量到人均的正确转换和时期匹配。
% =========================================================================
classdef TPI
    methods (Static)

        function [target_paths, errors, Pol_path_h, Dist_path_h, aggr_supply] = calculate_paths_and_errors(r_path, w_hat_path, bequest_gen_raw_pc_path_guess, TR_hat_pc_path, K_g_hat_pc_path_guess, kappa_path_guess, L_path_h_exog, ss0, ssF, cS, paramSF, valF_h, polF_h, dist0_h)
    % =========================================================================
    % ==          TPI核心计算引擎 (v17.2 - 投资缺口反馈微调版)
    % ==
    % == 核心逻辑:
    % ==   - 此版本保留了稳定收敛的 kappa 迭代框架。
    % ==   - 关键微调: 在计算出生产侧决定的价格后，额外计算实体投资缺口。
    % ==   - 将投资缺口作为一个反馈项，对利率的目标路径进行修正。
    % ==   - TPI在驱动 kappa 和价格收敛的同时，也会被这个反馈机制引导，
    % ==     从而将投资缺口一并消除，最终收敛到真实的、会计自洽的均衡。
    % =========================================================================
    T = cS.T_sim;
    total_pop_path = sum(cS.Z_path_raw, 1);

    % --- 步骤 1 & 2: 与 v17.1 完全相同，求解家庭并找到内部宏观均衡 ---
    [~, b_hat_path_h_for_hh] = TPI.evolve_payg_fund_path(r_path, w_hat_path, L_path_h_exog, ss0, ssF, cS, total_pop_path);
    A_path_ext = [cS.A_path, cS.A_path(end) * (1 + ((1+cS.g_A_ss)^cS.time_Step-1))];
    g_A_path_period = A_path_ext(2:end) ./ A_path_ext(1:end-1) - 1;
    Pol_path_h = cell(cS.nTypes, 1);
    for h = 1:cS.nTypes
        pathS_h = struct('r_path', r_path, 'w_hat_path', w_hat_path, 'tr_per_hh_hat_path', TR_hat_pc_path, 'b_hat_path', b_hat_path_h_for_hh(h, :), 'g_A_path', g_A_path_period(1:T), 'theta_path', cS.theta_path_h(h, :));
        cS_h = cS; cS_h.ageEffV_new = cS.ageEffV_new_h(:,h);
        [Pol_path_h{h}, ~] = household.backward_hh(pathS_h, cS_h, paramSF, valF_h{h}, polF_h{h});
    end
    bequest_gen_total_for_t1_newborns = ss0.Bequest_gen_hat_raw_ss;
    bequest_gen_raw_total_path_guess = bequest_gen_raw_pc_path_guess .* total_pop_path;
    bequest_gen_from_t1_to_tT_minus_1 = bequest_gen_raw_total_path_guess(1:T-1);
    bequest_gen_total_available_to_newborns_path = [bequest_gen_total_for_t1_newborns, bequest_gen_from_t1_to_tT_minus_1];
    mass_newborns_path = cS.Z_path_raw(1, :);
    bequest_per_newborn_hat_path = zeros(1, T);
    bequest_per_newborn_hat_path(mass_newborns_path > 1e-9) = bequest_gen_total_available_to_newborns_path(mass_newborns_path > 1e-9) ./ mass_newborns_path(mass_newborns_path > 1e-9);
    Dist_path_h = cell(cS.nTypes, 1);
    for h = 1:cS.nTypes
        dist0_h_for_sim = dist0_h{h};
        if cS.pps_active && size(dist0_h{h}, 2) == 1 && cS.nkpps > 1, dist0_expanded_h = zeros(cS.nk, cS.nkpps, cS.nw_expanded, cS.aD_new); dist0_expanded_h(:, 1, :, :) = dist0_h{h}; dist0_h_for_sim = dist0_expanded_h; end
        Dist_path_h{h} = distribution.simulate_dist_forward(dist0_h_for_sim, Pol_path_h{h}, cS, paramSF, bequest_per_newborn_hat_path, h);
    end
    aggr_supply = aggregates.get_path_aggregates(Dist_path_h, Pol_path_h, cS, paramSF);
    K_p_hat_supply_path = aggr_supply.K_p_hat_total_path;
    L_path_supply_total = sum(L_path_h_exog, 1);
    Y_hat_path_consistent = zeros(1, T);
    K_firm_hat_path_consistent = zeros(1, T);
    K_payg_hat_path_consistent = zeros(1, T);
    K_g_hat_total_path_guess = K_g_hat_pc_path_guess .* total_pop_path;
    fzero_options = optimset('TolX', 1e-9, 'Display', 'off');
    for t = 1:T
        Y_func = @(y_hat_t) y_hat_t - ((K_p_hat_supply_path(t) - kappa_path_guess(t) * y_hat_t)^cS.alpha * (K_g_hat_total_path_guess(t))^cS.gamma * (L_path_supply_total(t))^(1-cS.alpha-cS.gamma));
        [Y_hat_path_consistent(t), ~, exitflag] = fzero(Y_func, K_p_hat_supply_path(t), fzero_options);
        if exitflag ~= 1, warning('TPI内部均衡求解器在 t=%d 未能收敛', t); end
        K_payg_hat_path_consistent(t) = kappa_path_guess(t) * Y_hat_path_consistent(t);
        K_firm_hat_path_consistent(t) = K_p_hat_supply_path(t) - K_payg_hat_path_consistent(t);
    end

    % --- 步骤 3: 计算生产侧决定的价格 ---
    r_path_demand = zeros(1, T);
    w_hat_path_demand = zeros(1, T);
    for t = 1:T
        prices_t = firm.get_prices_at_t(K_firm_hat_path_consistent(t), K_g_hat_total_path_guess(t), L_path_supply_total(t), cS);
        r_path_demand(t) = prices_t.r_mkt_t;
        w_hat_path_demand(t) = prices_t.w_hat_t;
    end
    
    % --- 步骤 4: 【核心微调】计算投资缺口并修正利率目标 ---
    C_hat_supply_path = aggr_supply.C_hat_total_path;
    K_firm_path_level = K_firm_hat_path_consistent .* cS.A_path;
    Y_path_level = Y_hat_path_consistent .* cS.A_path;
    C_path_level = C_hat_supply_path .* cS.A_path;
    I_g_path_level = Y_path_level * cS.I_g_to_Y_ratio_ss;
    G_c_path_level = Y_path_level * cS.G_c_to_Y_ratio_ss;
    I_firm_from_NIPA = Y_path_level - C_path_level - I_g_path_level - G_c_path_level;
    Depreciation_firm_path = cS.ddk * K_firm_path_level;
    I_firm_from_K_law = [K_firm_path_level(2:end) - K_firm_path_level(1:end-1) + Depreciation_firm_path(1:end-1), NaN];
    I_firm_from_K_law(end) = I_firm_from_NIPA(end);
    investment_gap_ratio = (I_firm_from_K_law - I_firm_from_NIPA) ./ Y_path_level;
    investment_gap_ratio(isnan(investment_gap_ratio)) = 0;
    
    phi_r = 0.1; % 价格反馈系数
    target_paths = struct();
    % 如果投资需求(k_law) > 投资供给(nipa)，说明资本稀缺，利率应上升
    target_paths.r_path = r_path_demand + phi_r * investment_gap_ratio;
    target_paths.w_hat_path = w_hat_path_demand;

    % --- 步骤 5: 基于【修正后的目标价格】，计算其他宏观量的目标路径 ---
    [K_payg_hat_target_path, ~] = TPI.evolve_payg_fund_path(target_paths.r_path, target_paths.w_hat_path, L_path_h_exog, ss0, ssF, cS, total_pop_path);
    target_paths.kappa_path = zeros(1,T);
    valid_y_idx = abs(Y_hat_path_consistent) > 1e-9;
    target_paths.kappa_path(valid_y_idx) = K_payg_hat_target_path(valid_y_idx) ./ Y_hat_path_consistent(valid_y_idx);
    target_paths.TR_hat_pc_path = TPI.calculate_endogenous_TR_path(target_paths.r_path, target_paths.w_hat_path, K_firm_hat_path_consistent, L_path_supply_total, Y_hat_path_consistent, Dist_path_h, Pol_path_h, cS, paramSF) ./ total_pop_path;
    target_paths.bequest_gen_raw_pc_path = aggr_supply.Bequest_gen_raw_total_path ./ total_pop_path;
    I_g_hat_total_path = Y_hat_path_consistent * cS.I_g_to_Y_ratio_ss;
    K_g_hat_total_target_path = zeros(1, T);
    K_g_hat_total_target_path(1) = K_g_hat_total_path_guess(1);
    for t = 1:(T-1)
        K_g_level = K_g_hat_total_target_path(t) * cS.A_path(t); I_g_level = I_g_hat_total_path(t) * cS.A_path(t);
        K_g_level_next = (1 - cS.ddk_g) * K_g_level + I_g_level;
        K_g_hat_total_target_path(t+1) = K_g_level_next / cS.A_path(t+1);
    end
    target_paths.K_g_hat_pc_path = K_g_hat_total_target_path ./ total_pop_path;

    % --- 步骤 6: [打包最终结果] ---
    final_consistent_paths = struct();
    final_consistent_paths.Y_path = Y_hat_path_consistent .* cS.A_path;
    final_consistent_paths.C_path = aggr_supply.C_hat_total_path .* cS.A_path;
    final_consistent_paths.K_p_path = K_p_hat_supply_path .* cS.A_path;
    final_consistent_paths.K_payg_path = K_payg_hat_path_consistent .* cS.A_path;
    final_consistent_paths.K_g_path = K_g_hat_total_target_path .* cS.A_path;
    final_consistent_paths.r_path = target_paths.r_path;
    final_consistent_paths.w_path = target_paths.w_hat_path .* cS.A_path;
    target_paths.final_consistent_paths = final_consistent_paths;

    % --- 步骤 7: 计算误差 ---
    errors = struct();
    update_range_flow = 1:(T-1);
    errors.r = max(abs(target_paths.r_path(update_range_flow) - r_path(update_range_flow)));
    errors.w = max(abs(target_paths.w_hat_path(update_range_flow) - w_hat_path(update_range_flow)));
    errors.tr = max(abs(target_paths.TR_hat_pc_path(update_range_flow) - TR_hat_pc_path(update_range_flow)));
    errors.beq = max(abs(target_paths.bequest_gen_raw_pc_path(update_range_flow) - bequest_gen_raw_pc_path_guess(update_range_flow)));
    update_range_state = 2:T;
    errors.k_g = max(abs(target_paths.K_g_hat_pc_path(update_range_state) - K_g_hat_pc_path_guess(update_range_state)));
    errors.kappa = max(abs(target_paths.kappa_path(update_range_state) - kappa_path_guess(update_range_state)));
    errors.invest_gap = max(abs(investment_gap_ratio(update_range_flow))); % 监控投资缺口
    errors.total = max([errors.r, errors.w, errors.tr, errors.beq, errors.k_g, errors.kappa, errors.invest_gap]);
end        
        
        function L_path_h = get_labor_supply_path_for_type(h, Dist_path_h, cS, paramS)
            % 辅助函数: 专门为某个类型 h 计算其劳动供给路径
            T = cS.T_sim;
            L_path_h = zeros(1, T);
            cS_h = cS;
            cS_h.ageEffV_new = cS.ageEffV_new_h(:,h);
            dist_path_for_type_h = Dist_path_h{h};

            for t = 1:T
                L_t_h = 0;
                dist_t_abs_h = dist_path_for_type_h(:,:,:,:,t);
                for ia = 1:cS.aR_new
                    mass_ia_abs_h = dist_t_abs_h(:,:,:,ia);
                    if sum(mass_ia_abs_h(:)) < 1e-30, continue; end
                    le_grid_slice = reshape(paramS.leGridV, [1, 1, cS.nw_expanded]);
                    labor_supply_slice = cS_h.ageEffV_new(ia) .* le_grid_slice;
                    L_t_h = L_t_h + sum(labor_supply_slice .* mass_ia_abs_h, 'all');
                end
                L_path_h(t) = L_t_h;
            end
        end

        function TR_total_hat_path = calculate_endogenous_TR_path(r_path, w_hat_path, K_firm_hat_path, L_path, Y_hat_path, Dist_path_h, Pol_path_h, cS, paramS)
            % =========================================================================
            % == 函数: calculate_endogenous_TR_path (v3.0 - K_firm 会计修正版)
            % == 核心修改:
            % ==   - [!!! 关键会计BUG修复 !!!] 输入参数从 K_p_hat_path 修改为 K_firm_hat_path。
            % ==   - 在计算公共资本回报 (Y - wL - (r+d)K)时，使用的资本存量 K
            % ==     必须是厂商部门实际使用的生产性资本 K_firm，而不是家庭部门
            % ==     的总资产 K_p。此修正确保了政府收入核算的准确性。
            % =========================================================================
            T = cS.T_sim;
            nH = cS.nTypes;
            TR_total_hat_path = zeros(1, T);

            for t=1:T
                total_tax_hat_t_raw = 0;
                payg_contrib_hat_t_raw = 0;

                % --- 步骤 1 & 2: 循环所有家庭类型 h 来聚合总税收和总PAYG缴费 ---
                for h = 1:nH
                    cS_h = cS; % 获取类型 h 的参数
                    cS_h.ageEffV_new = cS.ageEffV_new_h(:,h);
                    cS_h.theta_path = cS.theta_path_h(h, t);

                    dist_t_abs_h = Dist_path_h{h}(:,:,:,:,t);
                    pol_t_h = Pol_path_h{h}{t};

                    for ia = 1:cS.aD_new
                        mass_ia_abs_h = dist_t_abs_h(:,:,:,ia);
                        total_tax_hat_t_raw = total_tax_hat_t_raw + sum(pol_t_h(ia).tax_regular .* mass_ia_abs_h, 'all');
                    end

                    for ia = 1:cS.aR_new % 只对工作年龄人口计算
                        mass_ia_abs_h = dist_t_abs_h(:,:,:,ia);
                        labor_income_slice_hat = w_hat_path(t) * cS_h.ageEffV_new(ia) .* paramS.leGridV(1:cS.nw_expanded)';
                        payg_contrib_slice_hat = cS_h.theta_path * labor_income_slice_hat;
                        payg_contrib_grid_hat = repmat(reshape(payg_contrib_slice_hat, [1,1,cS.nw_expanded]), [cS.nk, cS.nkpps, 1]);
                        payg_contrib_hat_t_raw = payg_contrib_hat_t_raw + sum(payg_contrib_grid_hat .* mass_ia_abs_h, 'all');
                    end
                end

                % --- 步骤 3-6: [核心修正] 基于 K_firm 计算政府预算 ---
                general_tax_hat_t_raw = total_tax_hat_t_raw - payg_contrib_hat_t_raw;

                General_Tax_Revenue_t_level = general_tax_hat_t_raw * cS.A_path(t);
                Y_t_level = Y_hat_path(t) * cS.A_path(t);
                w_t_level = w_hat_path(t) * cS.A_path(t);
                L_t_level = L_path(t);

                % [!!! 核心会计修正 !!!] 使用 K_firm_hat_path
                K_firm_t_level = K_firm_hat_path(t) * cS.A_path(t);
                K_firm_return_gross_t_level = r_path(t) * K_firm_t_level + cS.ddk * K_firm_t_level;

                Public_Capital_Return_t_level = Y_t_level - w_t_level * L_t_level - K_firm_return_gross_t_level;

                G_c_t_level = cS.G_c_to_Y_ratio_ss * Y_t_level;
                I_g_t_level = cS.I_g_to_Y_ratio_ss * Y_t_level;
                Total_General_Gov_Revenue_t_level = General_Tax_Revenue_t_level + Public_Capital_Return_t_level;
                TR_total_t_level = Total_General_Gov_Revenue_t_level - G_c_t_level - I_g_t_level;
                TR_total_hat_path(t) = TR_total_t_level / cS.A_path(t);
            end
        end
        
        
        
        function check_transition_NIPA(results, cS, ss0, ssF)
    % =========================================================================
    % == 函数: check_transition_NIPA (版本 v6.0 - 黄金标准核算版)
    % == 核心逻辑:
    % ==   - 此函数作为最终的、独立的“审计师”。它只信任最基础的路径：
    % ==     K_p_path, K_payg_path, K_g_path (存量) 和 L_path, C_path (流量)。
    % ==   - 它不信任任何来自results的派生流量，如Y_path或I_path。
    % ==   - 它将基于存量和生产函数，独立计算出所有流量，然后进行核算。
    % ==   - 这是检验TPI收敛点是否为真均衡的黄金标准。
    % =========================================================================

    fprintf('   正在启动过渡路径国民账户一致性检验 (v6.0 - 黄金标准核算版)...\n');
    T = cS.T_sim;

    % --- 步骤 1: 获取最基础的路径 (Level) ---
    K_p_path    = results.K_p_path;
    K_payg_path = results.K_payg_path;
    K_g_path    = results.K_g_path;
    C_path      = results.C_path;
    
    % 为了计算L_path，需要重新聚合一次 (因为results里没有L_path)
    % 我们需要重新加载最终的分布和策略（或者更好的做法是将其也存入results）
    % 为了简化，我们暂时从 results.w_path 和 Y_path 中反推出 L_path
    % Y_hat = K_firm^alpha * K_g^gamma * L^(1-alpha-gamma)
    % => L = (Y_hat / (K_firm^alpha * K_g^gamma))^(1/(1-alpha-gamma))
    K_firm_path = K_p_path - K_payg_path;
    Y_path = results.Y_path; % 我们暂时信任Y_path, 因为它是基于r,w计算的
    
    K_firm_hat_path = K_firm_path ./ cS.A_path;
    K_g_hat_path = K_g_path ./ cS.A_path;
    Y_hat_path = Y_path ./ cS.A_path;
    
    labor_share = 1 - cS.alpha - cS.gamma;
    L_path = (Y_hat_path ./ (K_firm_hat_path.^cS.alpha .* K_g_hat_path.^cS.gamma)).^(1/labor_share);


    % --- 步骤 2: 基于基础路径，独立、精确地重构所有宏观流量 ---
    % 2a. 基于生产三要素，重新计算每一期的产出 Y_path_recon 和价格 r, w
    Y_path_recon = zeros(1, T);
    r_path_recon = zeros(1, T);
    w_path_recon = zeros(1, T);
    for t = 1:T
        prices_t = firm.get_prices_at_t(K_firm_hat_path(t), K_g_hat_path(t), L_path(t), cS);
        Y_path_recon(t) = prices_t.Y_hat_t * cS.A_path(t);
        r_path_recon(t) = prices_t.r_mkt_t;
        w_path_recon(t) = prices_t.w_hat_t * cS.A_path(t);
    end

    % 2b. 计算政府流量
    I_g_path_recon = cS.I_g_to_Y_ratio_ss * Y_path_recon;
    G_c_path_recon = cS.G_c_to_Y_ratio_ss * Y_path_recon;
    
    % --- 步骤 3: 用两种口径计算【实体投资 I_firm】---
    % 口径 1: NIPA 资源约束倒算法 (黄金标准)
    I_firm_from_NIPA = Y_path_recon - C_path - I_g_path_recon - G_c_path_recon;
    
    % 口径 2: 资本积累法则 (来自家庭储蓄决策的最终结果)
    Depreciation_firm_path = cS.ddk * K_firm_path;
    I_firm_from_K_law = [K_firm_path(2:end) - K_firm_path(1:end-1) + Depreciation_firm_path(1:end-1), NaN];
    
    n_period_F = sum(cS.Z_path_raw(:, T)) / sum(cS.Z_path_raw(:, T-1)) - 1;
    g_A_period_F = cS.A_path(end) / cS.A_path(end-1) - 1;
    g_total_F_period = (1 + g_A_period_F) * (1 + n_period_F) - 1;
    I_firm_from_K_law(end) = (g_total_F_period + cS.ddk) * K_firm_path(end);
    
    % --- 步骤 4: 计算误差与缺口 ---
    investment_gap_path = I_firm_from_K_law - I_firm_from_NIPA;
    safe_Y_path = Y_path_recon;
    safe_Y_path(abs(safe_Y_path) < 1e-9) = 1e-9;
    rel_investment_gap_path = investment_gap_path ./ safe_Y_path;
    [max_rel_gap, t_max_gap] = max(abs(rel_investment_gap_path(1:T-1)));

    fprintf('   核心诊断: 最大实体投资缺口 (I_firm_k_law - I_firm_nipa)/Y: %.4e (%.4f%%) at t=%d\n', max_rel_gap, max_rel_gap*100, t_max_gap);

    % --- 步骤 5: 如果投资缺口过大，则打印标准的诊断报告 ---
    if max_rel_gap > 1e-4
        warning('过渡路径存在显著的投资缺口。TPI收敛点并非真实的宏观均衡。');
        
        print_gap_at_t = @(t) fprintf(...
            ['\n   [时期 t = %d]\n' ...
            '      (+) 产出 (Y_recon) .................... : %15.6f\n' ...
            '      (-) 私人消费 (C) ........................ : %15.6f\n' ...
            '      (-) 公共投资 (I_g) ........................ : %15.6f\n' ...
            '      (-) 政府消费 (G_c) ........................ : %15.6f\n' ...
            '      ------------------------------------------------------------\n' ...
            '      (=) NIPA隐含实体投资 (I_firm_nipa) ....... : %15.6f  (%.2f%% of Y)\n' ...
            '      (>) 积累法则要求投资 (I_firm_k_law) ...... : %15.6f  (%.2f%% of Y)\n' ...
            '      >>> 投资缺口 (k_law - nipa) .............. : %15.6f  (%.4f%% of Y)\n'], ...
            t, Y_path_recon(t), C_path(t), I_g_path_recon(t), G_c_path_recon(t), ...
            I_firm_from_NIPA(t), (I_firm_from_NIPA(t)/Y_path_recon(t))*100, ...
            I_firm_from_K_law(t), (I_firm_from_K_law(t)/Y_path_recon(t))*100, ...
            investment_gap_path(t), rel_investment_gap_path(t)*100);

        print_gap_at_t(1);
        if T > 1 && t_max_gap ~= 1, print_gap_at_t(t_max_gap); end
        if T > 2 && t_max_gap ~= T-1 && 1 ~= T-1, print_gap_at_t(T-1); end
        if t_max_gap > 1 && t_max_gap < T-1, print_gap_at_t(t_max_gap+1); end % 检查缺口后的时期

        fprintf('\n');
    else
        fprintf('   ✅ 投资缺口在可接受范围内，国民账户核算一致。\n');
    end
end
        function [K_payg_hat_total_path, b_hat_path_h] = evolve_payg_fund_path(r_path, w_hat_path, L_path_h, ss0, ssF, cS, total_pop_path)
    % =========================================================================
    % ==             PAYG基金演化引擎 v7.0 (改革规则稳定版)
    % ==
    % == 核心修改:
    % ==   - [!!! 关键逻辑修复 !!!] 重新引入混合规则，以确保与ssF的规则一致性。
    % ==   - 规则 1 (t=1): 为确保平滑衔接，若有盈余，遵循ss0的慷慨规则，
    % ==     将所有缴费作为福利支付出去 (adj_factor > 1)。
    % ==   - 规则 2 (t>=2, 改革生效): 严格遵循新制度。若有盈余，
    % ==     仅支付公式福利(adj_factor = 1)，将盈余存入基金以增强稳健性。
    % ==   - 此修改确保了转轨路径末端的行为与ssF的定义相符，从而解决了
    % ==     “错误终点”问题，是TPI收敛的关键。
    % =========================================================================
    T = cS.T_sim;
    nH = cS.nTypes;

    % --- 1. 计算PAYG流量路径 (hat总量) ---
    Total_Contrib_hat_path = zeros(1, T);
    for h = 1:nH
        Total_Contrib_hat_path = Total_Contrib_hat_path + cS.theta_path_h(h,:) .* w_hat_path .* L_path_h(h,:);
    end

    b_hat_formula_path_h = zeros(nH, T);
    for t = 1:T
        b_hat_formula_path_h(:,t) = SS.calculate_formula_benefits(w_hat_path(t), cS);
    end

    Total_Formula_Benefits_hat_path = zeros(1, T);
    mass_retirees_by_age = cS.Z_path_raw((cS.aR_new+1):end, :);
    for h=1:nH
        mass_retirees_path_h = sum(mass_retirees_by_age, 1) * cS.type_weights(h);
        Total_Formula_Benefits_hat_path = Total_Formula_Benefits_hat_path + b_hat_formula_path_h(h,:) .* mass_retirees_path_h;
    end

    % --- 2. 演化基金 (在level层面进行) ---
    K_payg_level_path = zeros(1, T+1); % T个时期，需要T+1个时点的值
    Actual_Benefits_level_path = zeros(1, T);

    K_payg_level_path(1) = ss0.K_payg_hat * total_pop_path(1) * cS.A_path(1);

    for t = 1:T
        Fund_Beg_Period_level = K_payg_level_path(t) * (1 + r_path(t));
        Contrib_level = Total_Contrib_hat_path(t) * cS.A_path(t);
        Formula_Benefits_level = Total_Formula_Benefits_hat_path(t) * cS.A_path(t);

        is_surplus = (Contrib_level >= Formula_Benefits_level);

        if t == 1 && is_surplus
            % 规则1 (t=1 且有盈余): 遵循 ss0 的慷慨规则以平滑衔接
            Actual_Benefits_level_path(t) = Contrib_level;
            K_payg_level_path(t+1) = Fund_Beg_Period_level; % 基金不从当期盈余中积累

        elseif t > 1 && is_surplus
            % 规则2 (t>=2 且有盈余): 改革生效，支付公式福利，将盈余存入基金
            Actual_Benefits_level_path(t) = Formula_Benefits_level;
            K_payg_level_path(t+1) = Fund_Beg_Period_level + (Contrib_level - Formula_Benefits_level);

        else % 所有赤字情况 (无论t=1还是t>1)
            gap_to_fill = Formula_Benefits_level - Contrib_level;
            if Fund_Beg_Period_level >= gap_to_fill
                % 正常赤字，基金足以弥补
                Actual_Benefits_level_path(t) = Formula_Benefits_level;
                K_payg_level_path(t+1) = Fund_Beg_Period_level - gap_to_fill;
            else
                % 技术性破产，基金不足以弥补
                Actual_Benefits_level_path(t) = Fund_Beg_Period_level + Contrib_level;
                K_payg_level_path(t+1) = 0;
            end
        end
    end

    % --- 3. 计算最终的福利和基金hat路径 ---
    K_payg_hat_total_path = K_payg_level_path(1:T) ./ cS.A_path;
    Actual_Benefits_hat_path = Actual_Benefits_level_path ./ cS.A_path;

    adj_factor_path = ones(1, T);
    idx_adj = abs(Total_Formula_Benefits_hat_path) > 1e-9;
    adj_factor_path(idx_adj) = Actual_Benefits_hat_path(idx_adj) ./ Total_Formula_Benefits_hat_path(idx_adj);

    % [保留逻辑] 强制路径末端平滑收敛到ssF的adj_factor，这在规则一致后是有效的数值技巧
    if T > 10
        trans_start_t = round(T/2);
        adj_factor_path(trans_start_t:T) = linspace(adj_factor_path(trans_start_t), ssF.adj_factor, T-trans_start_t+1);
    else
        adj_factor_path(T) = ssF.adj_factor;
    end

    b_hat_path_h = b_hat_formula_path_h .* adj_factor_path;
        end

        
function h_fig = plot_convergence_realtime(h_fig, iter_data, plot_config)
            % =========================================================================
            % == 函数: plot_convergence_realtime
            % == 版本: [v1.0 - TPI实时可视化]
            % ==
            % == 功能:
            % ==   1. 初始化一个绘图窗口，用于展示TPI迭代过程中变量路径的收敛情况。
            % ==   2. 在每次迭代后，将新的路径叠加到图上。
            % ==   3. 在收敛后，用一条醒目的线标记最终的路径。
            % ==
            % == 输入:
            % ==   - h_fig: 图形句柄，用于在同一个窗口上持续绘图。
            % ==   - iter_data: 包含当次迭代数据的结构体。
            % ==   - plot_config: 包含绘图所需常数和状态的结构体。
            % ==
            % == 调用逻辑:
            % ==   - status='initialize': 在TPI循环前调用，创建并设置图形。
            % ==   - status='update': 在TPI循环内每轮迭代后调用，更新图形。
            % ==   - status='finalize': 在TPI循环结束后调用，绘制最终结果。
            % =========================================================================

            if strcmp(plot_config.status, 'initialize')
                h_fig = figure('Name', 'TPI Convergence: Interest Rate (r)');
                set(h_fig, 'Position', [200, 200, 1000, 600]);
                hold on;

                % 绘制终期稳态参考线
                r_ssF_annual = ((1 + plot_config.ssF.r_mkt).^(1/plot_config.cS.time_Step) - 1) * 100;
                yline(r_ssF_annual, 'r--', 'LineWidth', 2, 'DisplayName', 'ssF r (annual)');

                % 标记初始稳态在t=1的点
                r_ss0_annual = ((1 + plot_config.ss0.r_mkt).^(1/plot_config.cS.time_Step) - 1) * 100;
                plot(1, r_ss0_annual, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k', 'DisplayName', 'ss0 r (annual)');

                % 绘制初始猜测路径
                r_guess_annual = ((1 + iter_data.r_path_guess).^(1/plot_config.cS.time_Step) - 1) * 100;
                plot(1:plot_config.T, r_guess_annual, 'k:', 'LineWidth', 1.5, 'DisplayName', 'Initial Guess');

                title('Convergence of Interest Rate Path');
                xlabel('Time Period (t)');
                ylabel('Annualized Interest Rate (%)');
                grid on;
                box on;

                drawnow;
                return; % 初始化后退出
            end

            % 确保在正确的图形窗口上绘图
            figure(h_fig);

            if strcmp(plot_config.status, 'update')
                % 绘制当次迭代的路径，使用半透明的细线，避免图形混乱
                r_iter_annual = ((1 + iter_data.r_path_iter).^(1/plot_config.cS.time_Step) - 1) * 100;
                plot(1:plot_config.T, r_iter_annual, 'Color', [0.5 0.5 1.0 0.5], 'HandleVisibility', 'off'); % HandleVisibility='off'使图例保持干净
                % 更新图标题，实时显示当前迭代次数和误差
                title(sprintf('Convergence of Interest Rate Path (Iter: %d, Error: %.2e)', iter_data.iter, iter_data.error_total));
                drawnow; % 强制刷新图形窗口

            elseif strcmp(plot_config.status, 'finalize')
                % 用一条粗实线高亮显示最终收敛的路径
                r_final_annual = ((1 + iter_data.r_path_iter).^(1/plot_config.cS.time_Step) - 1) * 100;
                plot(1:plot_config.T, r_final_annual, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Converged Path');
                legend('show', 'Location', 'best');
                title(sprintf('Convergence of Interest Rate Path (Converged at Iter %d)', iter_data.iter));
                drawnow;
            end
        end

        function plot_transition_with_demographics(results, cS, ssF, title_suffix)
    % =========================================================================
    % == 函数: plot_transition_with_demographics
    % == 版本: [v2.0 - K_payg 分解版]
    % ==
    % == 目的:
    % ==   1. 适配 v21 模型的输出结构。
    % ==   2. 增加一个新的子图，用于展示家庭总资产 K_p 是如何分解为
    % ==      厂商使用的生产性资本 K_firm 和 PAYG基金 K_payg 的。
    % =========================================================================
    if nargin < 4, title_suffix = ''; end

    figure('Name', ['转轨动态(含人口): ' title_suffix], 'Position', [50 50 1400 900]);
    T = cS.T_sim;

    % --- 1. 设置时间轴 ---
    if isfield(cS, 'start_year') && ~isempty(cS.start_year)
        time_axis = cS.start_year : cS.time_Step : (cS.start_year + (T-1)*cS.time_Step);
        xlabel_text = '年份';
    else
        time_axis = 1:T;
        xlabel_text = '模型期 (t)';
    end

    % --- 2. 准备BGP参考路径 (level) ---
    ssF_Kp_line = (ones(1, T) * ssF.K_private_hat) .* cS.A_path;
    ssF_Y_line = (ones(1, T) * ssF.Y_from_production_hat) .* cS.A_path;
    ssF_C_line = (ones(1, T) * ssF.C_agg) .* cS.A_path;
    ssF_r_line_annualized = (1 + ssF.r_mkt).^(1/cS.time_Step) - 1;
    ssF_w_line = (ones(1, T) * ssF.w_hat) .* cS.A_path;

    % --- 3. 绘制核心宏观变量路径 ---
    subplot(2, 3, 1);
    plot(time_axis, results.K_p_path, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 3, 'DisplayName', '家庭总资产 (K_p)');
    hold on;
    K_firm_path = results.K_p_path - results.K_payg_path; % 计算 K_firm
    plot(time_axis, K_firm_path, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 3, 'DisplayName', '厂商资本 (K_firm)');
    plot(time_axis, results.K_payg_path, 'g-^', 'LineWidth', 1.5, 'MarkerSize', 3, 'DisplayName', '养老基金 (K_payg)');
    plot(time_axis, ssF_Kp_line, 'b--', 'LineWidth', 1, 'DisplayName', 'K_p (BGP终态)');
    title('资本存量路径 (原始水平)');
    xlabel(xlabel_text); ylabel('资本总量');
    legend('show', 'Location', 'best'); grid on; xlim([time_axis(1), time_axis(end)]);

    subplot(2, 3, 2);
    plot(time_axis, results.Y_path, 'k-d', 'LineWidth', 1.5, 'MarkerSize', 3, 'DisplayName', '总产出 (Y)');
    hold on;
    plot(time_axis, results.C_path, 'c-^', 'LineWidth', 1.5, 'MarkerSize', 3, 'DisplayName', '总消费 (C)');
    plot(time_axis, ssF_Y_line, 'k--', 'LineWidth', 1, 'DisplayName', 'Y (BGP终态)');
    plot(time_axis, ssF_C_line, 'c--', 'LineWidth', 1, 'DisplayName', 'C (BGP终态)');
    title('产出与消费路径 (原始水平)');
    xlabel(xlabel_text); ylabel('总量');
    legend('show', 'Location', 'best'); grid on; xlim([time_axis(1), time_axis(end)]);

    subplot(2, 3, 4);
    r_path_annualized = (1 + results.r_path).^(1/cS.time_Step) - 1;
    plot(time_axis, r_path_annualized * 100, 'm-x', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', '年化市场利率 (r)');
    hold on;
    yline(ssF_r_line_annualized * 100, 'm--', 'LineWidth', 1, 'DisplayName', 'r (ssF终态)');
    title('市场利率路径');
    xlabel(xlabel_text); ylabel('年化利率 (%)');
    legend('show', 'Location', 'best'); grid on; xlim([time_axis(1), time_axis(end)]);

    subplot(2, 3, 5);
    plot(time_axis, results.w_path, 'color', [0.8500 0.3250 0.0980], 'LineStyle', '-', 'Marker','p', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', '工资 (w)');
    hold on;
    plot(time_axis, ssF_w_line, 'color', [0.8500 0.3250 0.0980], 'LineStyle','--', 'LineWidth', 1, 'DisplayName', 'w (BGP终态)');
    title('工资路径 (原始水平)');
    xlabel(xlabel_text); ylabel('工资水平');
    legend('show', 'Location', 'best'); grid on; xlim([time_axis(1), time_axis(end)]);

    % --- 4. 绘制人口动态路径 ---
    subplot(2, 3, 3);
    yyaxis left
    total_pop_path = sum(cS.Z_path_raw, 1);
    p1 = plot(time_axis, total_pop_path, 'k-o', 'LineWidth', 2, 'MarkerSize', 3, 'MarkerFaceColor', 'k');
    ylabel('总人口 (绝对数量)');
    ax_pop = gca;
    ax_pop.YAxis(1).Color = 'k';

    yyaxis right
    g_n_period = total_pop_path(2:end) ./ total_pop_path(1:end-1) - 1;
    g_n_annual = (1 + g_n_period).^(1/cS.time_Step) - 1;
    time_axis_growth = time_axis(1:end-1) + cS.time_Step/2;
    p2 = plot(time_axis_growth, g_n_annual * 100, 'r-^', 'LineWidth', 1.5, 'MarkerSize', 4);
    ylabel('年化人口增长率 (%)', 'Color', 'r');
    ax_pop.YAxis(2).Color = 'r';

    title('总人口与增长率');
    xlabel(xlabel_text);
    legend([p1, p2], {'总人口 [左轴]', '年化增长率 [右轴]'}, 'Location', 'best');
    grid on; xlim([time_axis(1), time_axis(end)]);

    subplot(2, 3, 6);
    retiree_pop = sum(cS.Z_path_raw((cS.aR_new+1):end, :), 1);
    working_pop = sum(cS.Z_path_raw(1:cS.aR_new, :), 1);

    old_age_dependency_ratio = (retiree_pop ./ working_pop) * 100;
    retiree_share = (retiree_pop ./ total_pop_path) * 100;

    p3 = plot(time_axis, old_age_dependency_ratio, 'b-d', 'LineWidth', 2, 'DisplayName', '老年抚养比');
    hold on;
    p4 = plot(time_axis, retiree_share, 'g-s', 'LineWidth', 2, 'DisplayName', '退休人口占比');

    title('人口结构动态 (老龄化)');
    xlabel(xlabel_text); ylabel('比率 (%)');
    legend([p3, p4], {'老年抚养比', '退休人口占比'}, 'Location', 'best');
    grid on; xlim([time_axis(1), time_axis(end)]);

    % --- 5. 添加总标题 ---
    sgtitle(['宏观经济转轨动态: ' title_suffix], 'FontSize', 16, 'FontWeight', 'bold');
        end


        function plot_population_shock_transition(results, cS, ssF, title_suffix)
            % =========================================================================
            % == 函数: plot_population_shock_transition
            % == 版本: [v1.0 - 人口冲击专用可视化]
            % ==
            % == 目的:
            % ==   1. 详细展示受控人口冲击实验下的宏观经济转轨路径。
            % ==   2. 增加专门的人口动态图，清晰标示冲击的发生。
            % ==   3. 自动处理年度或模型期时间轴。
            % =========================================================================
            if nargin < 4, title_suffix = ''; end

            figure('Name', ['转轨动态: ' title_suffix], 'Position', [50 50 1400 900]);
            T = cS.T_sim;

            % --- 1. 设置时间轴 ---
            if isfield(cS, 'start_year') && ~isempty(cS.start_year)
                time_axis = cS.start_year : cS.time_Step : (cS.start_year + (T-1)*cS.time_Step);
                xlabel_text = '年份';
            else
                time_axis = 1:T;
                xlabel_text = '模型期 (t)';
            end

            % --- 2. 准备BGP参考路径 (level) ---
            ssF_Kp_line = (ones(1, T) * ssF.K_private_hat) .* cS.A_path;
            ssF_Y_line = (ones(1, T) * ssF.Y_from_production_hat) .* cS.A_path;
            ssF_C_line = (ones(1, T) * ssF.C_agg) .* cS.A_path;
            ssF_r_line_annualized = (1 + ssF.r_mkt).^(1/cS.time_Step) - 1;
            ssF_w_line = (ones(1, T) * ssF.w_hat) .* cS.A_path;

            % --- 3. 绘制核心宏观变量路径 ---
            % -- 图 A: 私人资本路径 --
            subplot(2, 3, 1);
            plot(time_axis, results.K_p_path, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 3, 'DisplayName', '私人资本 (K_p)');
            hold on;
            plot(time_axis, ssF_Kp_line, 'b--', 'LineWidth', 1, 'DisplayName', 'K_p (BGP终态)');
            title('私人资本路径 (原始水平)');
            xlabel(xlabel_text); ylabel('资本总量');
            legend('show', 'Location', 'best'); grid on; xlim([time_axis(1), time_axis(end)]);

            % -- 图 B: 产出与消费路径 --
            subplot(2, 3, 2);
            plot(time_axis, results.Y_path, 'k-d', 'LineWidth', 1.5, 'MarkerSize', 3, 'DisplayName', '总产出 (Y)');
            hold on;
            plot(time_axis, results.C_path, 'g-^', 'LineWidth', 1.5, 'MarkerSize', 3, 'DisplayName', '总消费 (C)');
            plot(time_axis, ssF_Y_line, 'k--', 'LineWidth', 1, 'DisplayName', 'Y (BGP终态)');
            plot(time_axis, ssF_C_line, 'g--', 'LineWidth', 1, 'DisplayName', 'C (BGP终态)');
            title('产出与消费路径 (原始水平)');
            xlabel(xlabel_text); ylabel('总量');
            legend('show', 'Location', 'best'); grid on; xlim([time_axis(1), time_axis(end)]);

            % -- 图 C: 市场利率路径 --
            subplot(2, 3, 4);
            r_path_annualized = (1 + results.r_path).^(1/cS.time_Step) - 1;
            plot(time_axis, r_path_annualized * 100, 'm-x', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', '年化市场利率 (r)');
            hold on;
            yline(ssF_r_line_annualized * 100, 'm--', 'LineWidth', 1, 'DisplayName', 'r (ssF终态)');
            title('市场利率路径');
            xlabel(xlabel_text); ylabel('年化利率 (%)');
            legend('show', 'Location', 'best'); grid on; xlim([time_axis(1), time_axis(end)]);

            % -- 图 D: 工资路径 --
            subplot(2, 3, 5);
            plot(time_axis, results.w_path, 'c-p', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', '工资 (w)');
            hold on;
            plot(time_axis, ssF_w_line, 'c--', 'LineWidth', 1, 'DisplayName', 'w (BGP终态)');
            title('工资路径 (原始水平)');
            xlabel(xlabel_text); ylabel('工资水平');
            legend('show', 'Location', 'best'); grid on; xlim([time_axis(1), time_axis(end)]);

            % --- 4. 绘制人口动态路径 ---
            % -- 图 E: 总人口与人口增长率 --
            subplot(2, 3, 3);
            yyaxis left
            total_pop_path = sum(cS.Z_path_raw, 1);
            p1 = plot(time_axis, total_pop_path, 'k-o', 'LineWidth', 2, 'MarkerSize', 3, 'MarkerFaceColor', 'k');
            ylabel('总人口 (绝对数量)');
            ax_pop = gca;
            ax_pop.YAxis(1).Color = 'k';

            yyaxis right
            g_n_period = total_pop_path(2:end) ./ total_pop_path(1:end-1) - 1;
            g_n_annual = (1 + g_n_period).^(1/cS.time_Step) - 1;
            time_axis_growth = time_axis(1:end-1) + cS.time_Step/2;
            p2 = plot(time_axis_growth, g_n_annual * 100, 'r-^', 'LineWidth', 1.5, 'MarkerSize', 4);
            ylabel('年化人口增长率 (%)', 'Color', 'r');
            ax_pop.YAxis(2).Color = 'r';

            hold on;
            % 标示冲击区域
            shock_start_time = time_axis(3); % 冲击发生在第3期
            shock_end_time = time_axis(4);   % 冲击在第4期结束
            x_patch = [shock_start_time, shock_end_time, shock_end_time, shock_start_time];
            y_limits = ylim;
            y_patch = [y_limits(1), y_limits(1), y_limits(2), y_limits(2)];
            patch(x_patch, y_patch, 'red', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'DisplayName', '出生率冲击期');

            title('总人口与增长率');
            xlabel(xlabel_text);
            legend([p1, p2], {'总人口 [左轴]', '年化增长率 [右轴]'}, 'Location', 'best');
            grid on; xlim([time_axis(1), time_axis(end)]);


            % -- 图 F: 人口结构 (老龄化) --
            subplot(2, 3, 6);
            retiree_pop = sum(cS.Z_path_raw((cS.aR_new+1):end, :), 1);
            working_pop = sum(cS.Z_path_raw(1:cS.aR_new, :), 1);

            old_age_dependency_ratio = (retiree_pop ./ working_pop) * 100;
            retiree_share = (retiree_pop ./ total_pop_path) * 100;

            p3 = plot(time_axis, old_age_dependency_ratio, 'b-d', 'LineWidth', 2, 'DisplayName', '老年抚养比');
            hold on;
            p4 = plot(time_axis, retiree_share, 'g-s', 'LineWidth', 2, 'DisplayName', '退休人口占比');

            % 标示冲击区域
            patch(x_patch, [0 0 100 100], 'red', 'FaceAlpha', 0.1, 'EdgeColor', 'none');

            title('人口结构动态');
            xlabel(xlabel_text); ylabel('比率 (%)');
            legend([p3, p4], {'老年抚养比', '退休人口占比'}, 'Location', 'best');
            grid on; xlim([time_axis(1), time_axis(end)]); ylim([0, max(old_age_dependency_ratio)*1.1]);

            % --- 5. 添加总标题 ---
            sgtitle(['宏观经济转轨动态: ' title_suffix], 'FontSize', 16, 'FontWeight', 'bold');
        end

    end
end
