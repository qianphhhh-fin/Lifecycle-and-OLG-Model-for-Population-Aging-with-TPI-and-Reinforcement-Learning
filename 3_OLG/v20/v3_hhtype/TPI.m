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

        function [target_paths, errors, Pol_path_h, Dist_path_h, aggr_supply] = calculate_paths_and_errors(r_path, w_hat_path, bequest_gen_raw_pc_path_guess, TR_hat_pc_path, b_hat_path_h_guess, K_g_hat_path_pe_guess, ss0, ssF, cS, paramSF, valF_h, polF_h, dist0_h)
    % =========================================================================
    % ==          TPI核心计算引擎 (v9.10 - 遗赠会计精确修正版)
    % ==
    % == 核心修正:
    % ==   - [!!! BUG修复 !!!] 修复了因访问 t=0 人口数据导致的索引错误。
    % ==   - [!!! 逻辑简化与精确化 !!!]
    % ==     1. 直接使用 ss0.Bequest_generated_agg 作为 t=1 新生儿可获得的
    % ==        【调整后】总遗赠，因为它在稳态计算时已包含人口增长调整。
    % ==     2. 对转轨期 t=1..T-1 产生的遗赠猜测值，逐期除以正确的跨期
    % ==        人口增长因子 pop(t+1)/pop(t)，得到调整后的遗赠路径。
    % ==     3. 在计算新的目标遗赠路径时，也应用了同样的人口增长调整。
    % =========================================================================
    T = cS.T_sim;
    nH = cS.nTypes;

    % --- 1. [基于旧猜测] 求解家庭策略和分布，聚合出宏观供给 ---
    A_path_ext = [cS.A_path, cS.A_path(end) * (1 + ((1+cS.g_A_ss)^cS.time_Step-1))];
    g_A_path_period = A_path_ext(2:end) ./ A_path_ext(1:end-1) - 1;

    Pol_path_h = cell(nH, 1);
    for h = 1:nH
        pathS_h = struct();
        pathS_h.r_path = r_path;
        pathS_h.w_hat_path = w_hat_path;
        pathS_h.tr_per_hh_hat_path = TR_hat_pc_path;
        pathS_h.b_hat_path = b_hat_path_h_guess(h, :);
        pathS_h.g_A_path = g_A_path_period(1:T);
        pathS_h.theta_path = cS.theta_path_h(h, :);
        cS_h = cS; cS_h.ageEffV_new = cS.ageEffV_new_h(:,h);
        [Pol_path_h{h}, ~] = household.backward_hh(pathS_h, cS_h, paramSF, valF_h{h}, polF_h{h});
    end
    
    total_pop_path = sum(cS.Z_path_raw, 1);
    
    % [!!! 核心修正: 精确构造新生儿可获得的遗赠路径 !!!]
    
    % 步骤 a: 获取 t=1 新生儿可获得的【调整后】总遗赠 (直接来自ss0)
    bequest_adj_total_for_t1 = ss0.Bequest_generated_agg .* total_pop_path(1);  % ss0.Bequest_generated_agg是总人口为1的总量所以也是人均量
    
    % 步骤 b: 获取猜测的 t=1...T-1 期产生的【原始】总遗赠
    bequest_gen_raw_total_path_guess = bequest_gen_raw_pc_path_guess .* total_pop_path;
    
    % 步骤 c: 计算 t=2...T 新生儿可获得的【调整后】总遗赠
    %         t期产生的原始遗赠，需除以(1+g_n(t+1)) = pop(t+1)/pop(t)
    pop_growth_factors = 1 ;% total_pop_path(2:T) ./ total_pop_path(1:T-1);
    bequest_adj_total_for_t2_to_T = bequest_gen_raw_total_path_guess(1:T-1) ./ pop_growth_factors;

    % 步骤 d: 组合成完整的、t=1...T 期【新生儿可获得的、已调整】总遗赠路径
    bequest_gen_total_available_to_newborns_path = [bequest_adj_total_for_t1, bequest_adj_total_for_t2_to_T];
    
    % 步骤 e: 计算每个新生儿收到的遗赠
    mass_newborns_path = cS.Z_path_raw(1, :);
    bequest_per_newborn_hat_path = zeros(1, T);
    bequest_per_newborn_hat_path(mass_newborns_path > 1e-9) = ...
        bequest_gen_total_available_to_newborns_path(mass_newborns_path > 1e-9) ./ mass_newborns_path(mass_newborns_path > 1e-9);

    Dist_path_h = cell(nH, 1);
    for h = 1:nH
        dist0_h_for_sim = dist0_h{h};
        if cS.pps_active && size(dist0_h{h}, 2) == 1 && cS.nkpps > 1
            dist0_expanded_h = zeros(cS.nk, cS.nkpps, cS.nw_expanded, cS.aD_new);
            dist0_expanded_h(:, 1, :, :) = dist0_h{h};
            dist0_h_for_sim = dist0_expanded_h;
        end
        cS_h = cS; cS_h.s_pathV = cS.s_pathV;
        Dist_path_h{h} = distribution.simulate_dist_forward(dist0_h_for_sim, Pol_path_h{h}, cS_h, paramSF, bequest_per_newborn_hat_path, h);
    end

    aggr_supply = aggregates.get_path_aggregates(Dist_path_h, Pol_path_h, cS, paramSF);
    L_path_supply = aggr_supply.L_path;
    K_p_hat_total_path_supply = aggr_supply.K_p_hat_total_path;

    % --- 2. [计算最终目标] 基于宏观供给，计算本轮最终的价格和宏观流量 ---
    target_paths = struct();
    target_paths.r_path = zeros(1, T);
    target_paths.w_hat_path = zeros(1, T);
    Y_hat_total_path_calc = zeros(1, T);
    K_g_hat_total_path_guess = K_g_hat_path_pe_guess .* total_pop_path;
    
    for t = 1:T
        prices_t = firm.get_prices_at_t(K_p_hat_total_path_supply(t), K_g_hat_total_path_guess(t), L_path_supply(t), cS);
        target_paths.r_path(t) = prices_t.r_mkt_t;
        target_paths.w_hat_path(t) = prices_t.w_hat_t;
        Y_hat_total_path_calc(t) = prices_t.Y_hat_t;
    end
    
    % --- 3. [强化自洽性] 基于本轮最终价格和分布，重新计算所有政府和遗赠路径 ---
    target_paths.b_hat_path_h = zeros(nH, T);
    for h = 1:nH
        L_path_supply_h = TPI.get_labor_supply_path_for_type(h, Dist_path_h, cS, paramSF);
        total_pension_pot_hat_path_h = cS.theta_path_h(h,:) .* target_paths.w_hat_path .* L_path_supply_h;
        mass_retirees_path_h = sum(cS.Z_path_raw((cS.aR_new+1):end, :), 1) * cS.type_weights(h);
        target_paths.b_hat_path_h(h, mass_retirees_path_h > 1e-9) = total_pension_pot_hat_path_h(mass_retirees_path_h > 1e-9) ./ mass_retirees_path_h(mass_retirees_path_h > 1e-9);
    end

    TR_total_hat_target_path = TPI.calculate_endogenous_TR_path(target_paths.r_path, target_paths.w_hat_path, K_p_hat_total_path_supply, L_path_supply, Y_hat_total_path_calc, Dist_path_h, Pol_path_h, cS, paramSF);
    target_paths.TR_hat_pc_path = TR_total_hat_target_path ./ total_pop_path;
    
    % [!!!] 遗赠路径的目标值也需要进行同样的人口调整
    bequest_gen_raw_new_path = aggr_supply.Bequest_gen_raw_total_path;
    bequest_adj_new_path = bequest_gen_raw_new_path(1:T-1) ./ pop_growth_factors;
    % 将调整后的总量转回下一轮的人均猜测值
    bequest_pc_new_guess = zeros(1, T);
    bequest_pc_new_guess(1:T-1) = bequest_gen_raw_new_path(1:T-1) ./ total_pop_path(1:T-1); % 目标值用未调整的原始遗赠
    bequest_pc_new_guess(T) = bequest_pc_new_guess(T-1); % 假设最后一期不变
    target_paths.bequest_gen_raw_pc_path = bequest_pc_new_guess;

    target_paths.K_g_hat_path = zeros(1, T);
    target_paths_K_g_hat_total = zeros(1, T);
    I_g_hat_total_path_calc = cS.I_g_to_Y_ratio_ss * Y_hat_total_path_calc; 
    target_paths_K_g_hat_total(1) = K_g_hat_total_path_guess(1);
    for t = 1:(T-1)
        target_paths_K_g_hat_total(t+1) = ((1 - cS.ddk_g) * target_paths_K_g_hat_total(t) + I_g_hat_total_path_calc(t)) / (1 + g_A_path_period(t));
    end
    target_paths.K_g_hat_path = target_paths_K_g_hat_total ./ total_pop_path;

    % --- 4. 计算最终误差 ---
    errors = struct();
    update_range_flow = 1:(T-1);
    errors.r = max(abs(target_paths.r_path(update_range_flow) - r_path(update_range_flow)));
    errors.w = max(abs(target_paths.w_hat_path(update_range_flow) - w_hat_path(update_range_flow)));
    errors.b_hat = max(abs(target_paths.b_hat_path_h(:, update_range_flow) - b_hat_path_h_guess(:, update_range_flow)), [], 'all');
    errors.tr = max(abs(target_paths.TR_hat_pc_path(update_range_flow) - TR_hat_pc_path(update_range_flow)));
    errors.beq = max(abs(target_paths.bequest_gen_raw_pc_path(update_range_flow) - bequest_gen_raw_pc_path_guess(update_range_flow)));

    update_range_state = 2:T;
    errors.k_g = max(abs(target_paths.K_g_hat_path(update_range_state) - K_g_hat_path_pe_guess(update_range_state)));
    errors.total = max([errors.r, errors.w, errors.b_hat, errors.beq, errors.tr, errors.k_g]);
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

        function TR_total_hat_path = calculate_endogenous_TR_path(r_path, w_hat_path, K_p_hat_path, L_path, Y_hat_path, Dist_path_h, Pol_path_h, cS, paramS)
            % =========================================================================
            % == 函数: calculate_endogenous_TR_path (v2.0 - 异质性版)
            % == 核心修改:
            % ==   - 函数接口接收异质性家庭的分布和策略单元数组 (Dist_path_h, Pol_path_h)。
            % ==   - 内部聚合税收和PAYG缴费时，通过循环遍历所有家庭类型 h 来完成。
            % =========================================================================
            T = cS.T_sim;
            nH = cS.nTypes;
            TR_total_hat_path = zeros(1, T);

            for t=1:T
                total_tax_hat_t_raw = 0;
                payg_contrib_hat_t_raw = 0;

                % --- 步骤 1 & 2: [核心修改] 循环所有家庭类型 h 来聚合总税收和总PAYG缴费 ---
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

                % --- 步骤 3-6: 逻辑与之前版本相同 ---
                general_tax_hat_t_raw = total_tax_hat_t_raw - payg_contrib_hat_t_raw;

                General_Tax_Revenue_t_level = general_tax_hat_t_raw * cS.A_path(t);
                Y_t_level = Y_hat_path(t) * cS.A_path(t);
                w_t_level = w_hat_path(t) * cS.A_path(t);
                L_t_level = L_path(t);
                K_p_t_level = K_p_hat_path(t) * cS.A_path(t);
                K_p_return_gross_t_level = r_path(t) * K_p_t_level + cS.ddk * K_p_t_level;
                Public_Capital_Return_t_level = Y_t_level - w_t_level * L_t_level - K_p_return_gross_t_level;
                G_c_t_level = cS.G_c_to_Y_ratio_ss * Y_t_level;
                I_g_t_level = cS.I_g_to_Y_ratio_ss * Y_t_level;
                Total_General_Gov_Revenue_t_level = General_Tax_Revenue_t_level + Public_Capital_Return_t_level;
                TR_total_t_level = Total_General_Gov_Revenue_t_level - G_c_t_level - I_g_t_level;
                TR_total_hat_path(t) = TR_total_t_level / cS.A_path(t);
            end
        end

        function check_transition_NIPA(results, cS, ss0, ssF)
            % =========================================================================
            % == 函数: check_transition_NIPA (版本 v3.1 - 深度调试版)
            % == 核心修改:
            % ==   - 增加了详细的 t=1 时刻变量诊断信息，用于追踪投资缺口的根源。
            % ==   - 打印 K(1), K(2), hat(K(1)), hat(K(2)) 等关键存量。
            % ==   - 报告 hat(Y(1)), hat(C(1)) 等关键流量。
            % =========================================================================

            fprintf('   正在启动过渡路径国民账户一致性检验 (v3.1 - 深度调试版)...\n');
            T = cS.T_sim;

            % --- 1. 获取所有流量和存量路径 (Level) ---
            Y_path = results.Y_path;
            C_path = results.C_path;
            I_g_path = results.I_g_path;
            G_path = cS.G_c_to_Y_ratio_ss * Y_path;
            K_p_path = results.K_p_path;

            % --- 2. 获取两种口径的投资 (Level) ---
            I_p_from_NIPA = results.I_p_path;
            I_p_from_K_law = results.I_p_from_K_law;

            % --- 3. 计算NIPA闭合情况和投资缺口 ---
            Total_Uses_path = C_path + I_p_from_NIPA + I_g_path + G_path;
            nipa_error_path = Y_path - Total_Uses_path;
            investment_gap_path = I_p_from_K_law - I_p_from_NIPA;

            safe_Y_path = Y_path;
            safe_Y_path(abs(safe_Y_path) < 1e-9) = 1e-9;
            rel_investment_gap_path = investment_gap_path ./ safe_Y_path;

            [max_rel_gap, t_max_gap] = max(abs(rel_investment_gap_path(1:T-1)));

            fprintf('   国民账户资源约束误差 (Y - [C+I_nipa+Ig+Gc]): %.4e (应接近于0)\n', max(abs(nipa_error_path(1:T-1))));
            fprintf('   核心诊断: 最大投资缺口 (I_k_law - I_nipa)/Y: %.4e (%.4f%%) at t=%d\n', max_rel_gap, max_rel_gap*100, t_max_gap);

            % --- 4. [!!! 深度调试信息 !!!] ---
            fprintf('\n   --- [深度调试] t=1 时刻详细诊断 ---\n');
            % 获取去趋势量 (hat量)
            total_pop_path = sum(cS.Z_path_raw, 1);
            K_p_hat_path = K_p_path ./ (cS.A_path .* total_pop_path);
            Y_hat_path = Y_path ./ (cS.A_path .* total_pop_path);
            C_hat_path = C_path ./ (cS.A_path .* total_pop_path);

            % 打印存量信息
            fprintf('      [存量 Stock]\n');
            fprintf('         初始稳态 ss0.K_private_hat .......... : %15.6f (来自稳态求解器)\n', ss0.K_private_hat);
            fprintf('         t=1 期初资本 K_p_hat(1) ............... : %15.6f (来自TPI聚合dist0)\n', K_p_hat_path(1));
            fprintf('         t=2 期初资本 K_p_hat(2) ............... : %15.6f (来自TPI聚合dist(2))\n', K_p_hat_path(2));
            fprintf('         K_p(1)_level .......................... : %15.6f\n', K_p_path(1));
            fprintf('         K_p(2)_level .......................... : %15.6f\n', K_p_path(2));

            % 打印流量信息
            fprintf('      [流量 Flow]\n');
            fprintf('         t=1 产出 Y_hat(1) ...................... : %15.6f\n', Y_hat_path(1));
            fprintf('         t=1 消费 C_hat(1) ...................... : %15.6f\n', C_hat_path(1));
            fprintf('         I_p_from_K_law(1) = K(2)-(1-d)K(1) ..... : %15.6f\n', I_p_from_K_law(1));
            fprintf('         I_p_from_NIPA(1) = Y(1)-C(1)-... ...... : %15.6f\n', I_p_from_NIPA(1));

            % 打印关键比率
            fprintf('      [关键比率]\n');
            fprintf('         ss0 K/Y ............................... : %15.6f\n', ss0.K_private_hat / ss0.Y_from_production_hat);
            fprintf('         t=1 K/Y ............................... : %15.6f\n', K_p_hat_path(1) / Y_hat_path(1));

            % 打印诊断结论
            k_diff = K_p_hat_path(1) - ss0.K_private_hat;
            fprintf('      [诊断]\n');
            fprintf('         K_p_hat(1) 与 ss0.K_private_hat 的差异 : %.4e\n', k_diff);
            if abs(k_diff / ss0.K_private_hat) > 1e-5
                fprintf('         >>> 警告: t=1的期初资本总量与初始稳态严重不符！这很可能是问题的根源。\n');
                fprintf('             请检查 test_diffN_ss.m 中 dist0 的构造 和 aggregates.m 中两种聚合方式是否完全一致。\n');
            else
                fprintf('         >>> 信息: t=1的期初资本总量与初始稳态基本一致。问题可能出在 K_p_hat(2) 的计算上。\n');
            end


            % --- 5. 如果投资缺口过大，则打印标准的诊断报告 ---
            if max_rel_gap > 1e-4
                warning('过渡路径存在显著的投资缺口，表明储蓄决策与宏观资源约束不完全一致。');
                fprintf('\n   --- 标准诊断报告 ---\n');

                print_gap_at_t = @(t) fprintf(...
                    ['\n   [时期 t = %d]\n' ...
                    '      (+) 产出 (Y) .................... : %15.6f\n' ...
                    '      (-) 私人消费 (C) ................ : %15.6f\n' ...
                    '      (-) 公共投资 (I_g) ................ : %15.6f\n' ...
                    '      (-) 政府消费 (G_c) ................ : %15.6f\n' ...
                    '      ------------------------------------------------------------\n' ...
                    '      (=) NIPA隐含投资 (I_nipa) ....... : %15.6f  (%.2f%% of Y)\n' ...
                    '      (>) 资本积累要求投资 (I_k_law) .. : %15.6f  (%.2f%% of Y)\n' ...
                    '      >>> 投资缺口 (I_k_law - I_nipa) . : %15.6f  (%.4f%% of Y)\n'], ...
                    t, Y_path(t), C_path(t), I_g_path(t), G_path(t), ...
                    I_p_from_NIPA(t), (I_p_from_NIPA(t)/Y_path(t))*100, ...
                    I_p_from_K_law(t), (I_p_from_K_law(t)/Y_path(t))*100, ...
                    investment_gap_path(t), rel_investment_gap_path(t)*100);

                print_gap_at_t(1);
                if t_max_gap ~= 1 && T > 1, print_gap_at_t(t_max_gap); end
                if T > 2 && t_max_gap ~= T-1 && 1 ~= T-1, print_gap_at_t(T-1); end
                fprintf('\n');
            else
                fprintf('   ✅ 投资缺口在可接受范围内，储蓄决策与宏观资源基本一致。\n');
            end
        end


        % --- TPI.m ---
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
            % == 版本: [v1.0 - 人口动态增强版]
            % ==
            % == 目的:
            % ==   1. 在一张图上综合展示宏观经济转轨路径与背后的人口动态。
            % ==   2. 增加专门的人口总量、增长率和结构（老龄化）图。
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
            plot(time_axis, results.K_p_path, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 3, 'DisplayName', '私人资本 (K_p)');
            hold on;
            plot(time_axis, ssF_Kp_line, 'b--', 'LineWidth', 1, 'DisplayName', 'K_p (BGP终态)');
            title('私人资本路径 (原始水平)');
            xlabel(xlabel_text); ylabel('资本总量');
            legend('show', 'Location', 'best'); grid on; xlim([time_axis(1), time_axis(end)]);

            subplot(2, 3, 2);
            plot(time_axis, results.Y_path, 'k-d', 'LineWidth', 1.5, 'MarkerSize', 3, 'DisplayName', '总产出 (Y)');
            hold on;
            plot(time_axis, results.C_path, 'g-^', 'LineWidth', 1.5, 'MarkerSize', 3, 'DisplayName', '总消费 (C)');
            plot(time_axis, ssF_Y_line, 'k--', 'LineWidth', 1, 'DisplayName', 'Y (BGP终态)');
            plot(time_axis, ssF_C_line, 'g--', 'LineWidth', 1, 'DisplayName', 'C (BGP终态)');
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
            plot(time_axis, results.w_path, 'c-p', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', '工资 (w)');
            hold on;
            plot(time_axis, ssF_w_line, 'c--', 'LineWidth', 1, 'DisplayName', 'w (BGP终态)');
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
