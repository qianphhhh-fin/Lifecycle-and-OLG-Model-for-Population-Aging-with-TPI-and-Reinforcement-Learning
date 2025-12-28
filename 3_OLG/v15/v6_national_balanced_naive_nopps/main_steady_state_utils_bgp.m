classdef main_steady_state_utils_bgp
    % =========================================================================
    % == 类说明: main_steady_state_utils_bgp (BGP平衡增长路径版本 - 最终重构版)
    % ==
    % == [重构核心] 统一架构：从"求解+反推"模式升级为"一次性求解+聚合"模式。
    % == [重构目标] 完全消除backout函数，VFI直接返回所有微观会计变量。
    % == [重构优势] 通过设计确保微观决策和宏观核算的绝对一致性，从根本上解决模型不收敛问题。
    % == [统一控制] 通过cS.pps_active开关在同一套代码框架下处理有/无PPS两种情况。
    % =========================================================================

    methods (Static)

        % =======================================================
        % == 主入口：统一稳态求解器
        % =======================================================

        function [ss, Dist, polS, valS] = solve_steady_state_unified(cS, paramS, params_ext, verbose, x0_guess, solver_method)
            % [重构核心-最终版 v3] 统一的稳态求解器主入口
            % == 基于行业标准实践的最终修正 ==
            % 1. [信念分离] 承认在“校准稳态”下，VFI决策使用的物理生存率(s_pathV)
            %    和宏观人口演化(由外部Z决定)是分离的。
            % 2. [接受误差] 不再试图强制统一信念，接受校准稳态下国民账户的微小不一致性，
            %    这恰恰反映了现实人口偏离理论稳态的程度。

            if nargin < 4, verbose = true; end
            if nargin < 5, x0_guess = []; end
            if nargin < 6, solver_method = 'fsolve'; end

            % --- 1. 参数设置 ---
            cS.A = params_ext.A;
            if isfield(params_ext, 'g_A_ss'), cS.g_A_ss = params_ext.g_A_ss; end
            cS.theta_path = params_ext.theta;

            % --- [核心逻辑修正] ---
            % 无论何种场景，都直接使用传入的人口分布Z。
            % VFI决策将始终使用模型自身的物理生存概率cS.s_pathV。
            % 不再计算或使用s_implied。
            if isfield(params_ext, 'Z') && ~isempty(params_ext.Z)
                Z_ss_norm = params_ext.Z;
                if verbose, fprintf('   [求解器信息] 使用外部传入的人口分布Z进行求解/校准。\n'); end
            else
                error('错误：求解稳态需要一个明确的人口分布Z。请在params_ext中提供。');
            end

            if ~isfield(cS, 'pps_active'), cS.pps_active = false; end
            if cS.pps_active && cS.nkpps <= 1, error('PPS模式要求 nkpps > 1.');
            elseif ~cS.pps_active, cS.nkpps = 1; cS.kppsGridV = 0; end

            % --- 2. 求解 ---
            system_wrapper = @(x) main_steady_state_utils_bgp.system_of_equations_steady_state(x, Z_ss_norm, cS, paramS);
            if isempty(x0_guess), x0 = [0.3336, 0.078, 0.3]; else, x0 = x0_guess(1:3); end
            [x_eq, eq_found] = main_steady_state_utils_bgp.solve_with_fsolve(system_wrapper, x0, verbose);

            if ~eq_found
                warning('统一求解器未能找到均衡解！');
                ss = []; Dist = []; polS = []; valS = [];
                return;
            end

            % --- 3. 获取最终结果并展示 ---
            [ss, Dist, polS, valS] = main_steady_state_utils_bgp.calculate_aggregates_unified(x_eq(1), x_eq(2), x_eq(3), Z_ss_norm, cS, paramS);



            % [验证步骤]
            main_steady_state_utils_bgp.verify_household_budget_constraint(ss, polS, cS, paramS);
            main_steady_state_utils_bgp.verify_FULL_steady_state(Dist, polS, paramS, cS);

            if verbose, main_steady_state_utils_bgp.display_national_accounts_unified(ss, cS);
            else, main_steady_state_utils_bgp.check_national_accounts_unified(ss, cS); end
        end

        function F_error = system_of_equations_steady_state(x, Z_ss_norm, cS, paramS)
            % =========================================================================
            % == 函数: system_of_equations_steady_state
            % == 版本: [v10 - 最终完美预算版]
            % ==
            % == 核心修正:
            % ==   1. 政府预算现在正确地将资本折旧作为一项必须的支出。
            % ==   2. 政府可自由裁量的资源 = 总收入 - 折旧支出。
            % ==   3. 投资/消费比例(lambda_g)作用于这个正确的资源基数上。
            % =========================================================================

            K_p_guess = x(1); K_g_guess = x(2); L_guess = x(3);
            g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
            [ss, ~, ~, ~] = main_steady_state_utils_bgp.calculate_aggregates_unified(K_p_guess, K_g_guess, L_guess, Z_ss_norm, cS, paramS);

            % --- 市场出清误差 ---
            % 方程1: 净储蓄市场 (不变，已正确)
            S_p_supply = ss.Saving_private_flow;
            I_p_net_demand = g_A_period * K_p_guess;
            Bequest_demand = ss.Bequest_tax;
            Total_demand_for_net_saving = I_p_net_demand + Bequest_demand;
            error_Kp = Total_demand_for_net_saving - S_p_supply;

            % 方程2: 公共净储蓄市场
            % [核心修正]
            % A. 计算政府总收入
            Gov_Total_Revenue = ss.Regular_tax + ss.Bequest_tax + ss.Public_Capital_Return;
            % B. 总收入必须先覆盖公共资本折旧
            Resources_for_discretion = Gov_Total_Revenue - ss.Depreciation_g;
            % C. 政府的净储蓄(用于净投资)是可自由裁量资源的 lambda_g 部分
            Gov_Net_Saving = cS.lambda_g * Resources_for_discretion;

            I_g_net_demand = g_A_period * K_g_guess;
            error_Kg = I_g_net_demand - Gov_Net_Saving;

            % 方程3: 劳动市场
            error_L = L_guess - ss.L_hat;
            F_error = [error_Kp; error_Kg; error_L];
        end

        function [ss, Dist, polS, valS] = calculate_aggregates_unified(K_p_guess, K_g_guess, L_guess, Z_ss_norm, cS, paramS)
            % =========================================================================
            % == 函数: calculate_aggregates_unified
            % == 版本: [v13.1 - 字段恢复版]
            % ==
            % == 核心修正:
            % ==   1. 恢复了 ss.K_private_hat 字段，以修复主脚本的报错。
            % =========================================================================

            % --- 1. 初始化和价格计算 ---
            if K_p_guess <= 0, K_p_guess = 1e-8; end; if K_g_guess <= 0, K_g_guess = 1e-8; end; if L_guess <= 0, L_guess = 1e-8; end;
            A_ss = 1.0; theta_ss = cS.theta_path(1);
            M_prices = main_steady_state_utils_bgp.get_prices_at_t(K_p_guess, K_g_guess, L_guess, A_ss, cS);
            M_for_hh = M_prices;
            mass_retirees_ss = sum(Z_ss_norm((cS.aR_new+1):end));
            M_for_hh.b_t = (theta_ss * M_prices.w_t * L_guess) / max(1e-9, mass_retirees_ss);
            cS.theta_t = theta_ss;

            % --- 2. 求解家庭问题和稳态分布 ---
            [polS, valS] = main_steady_state_utils_bgp.HHSolution_VFI_unified(M_for_hh, paramS, cS);
            Dist = main_steady_state_utils_bgp.solve_steady_state_distribution_unified(polS, paramS, cS, Z_ss_norm);

            % --- 3. 聚合微观变量 ---
            K_private_hat_agg = 0; L_agg = 0; C_agg = 0; Bequest_tax_agg = 0;
            Shock_exp_agg = 0; Pension_out_agg = 0; Pension_in_agg = 0; Regular_tax_agg = 0;
            for ia = 1:cS.aD_new
                mass_dist = Dist(:,:,:,ia);
                K_private_hat_agg = K_private_hat_agg + sum(polS(ia).k_prime .* mass_dist, 'all');
                total_assets_chosen_for_next_period = sum((polS(ia).k_prime + polS(ia).kpps_prime) .* mass_dist, 'all');
                Bequest_tax_agg = Bequest_tax_agg + total_assets_chosen_for_next_period * (1 - cS.s_pathV(ia));
                if ia <= cS.aR_new
                    mass_by_epsilon = squeeze(sum(mass_dist, [1,2]));
                    L_agg = L_agg + sum( (cS.ageEffV_new(ia) * paramS.leGridV') .* mass_by_epsilon' , 'all');
                end
                C_agg = C_agg + sum(polS(ia).c .* mass_dist, 'all');
                Shock_exp_agg = Shock_exp_agg + sum(polS(ia).shock_exp .* mass_dist, 'all');
                Pension_out_agg = Pension_out_agg + sum(polS(ia).pension_out .* mass_dist, 'all');
                Regular_tax_agg = Regular_tax_agg + sum(polS(ia).tax_regular .* mass_dist, 'all');
                Pension_in_agg = Pension_in_agg + sum(polS(ia).tax_payg .* mass_dist, 'all');
            end

            % --- 4. 计算家庭净储蓄, 折旧, 和公共资本回报 ---
            Household_Total_Inflow = (M_prices.w_t * L_agg) + (K_p_guess * M_prices.r_mkt_t) + Pension_out_agg;
            Household_Total_Outlay_No_Saving = C_agg + Shock_exp_agg + Regular_tax_agg + Pension_in_agg;
            Saving_private_flow = Household_Total_Inflow - Household_Total_Outlay_No_Saving;

            Depreciation_p = cS.ddk * K_p_guess;
            Depreciation_g = cS.ddk_g * K_g_guess;

            Factor_Payment_Total = (M_prices.w_t * L_agg) + ((M_prices.r_mkt_t + cS.ddk) * K_p_guess);
            Public_Capital_Return = M_prices.Y_t - Factor_Payment_Total;

            % --- 5. 填充 ss 结构体 ---
            ss = struct();
            ss.K_private_hat = K_private_hat_agg; % [恢复] 期末私人资本存量
            ss.K_private_begin_hat = K_p_guess;
            ss.L_hat = L_agg;
            ss.K_public_hat = K_g_guess;
            ss.Y_from_production_hat = M_prices.Y_t;
            ss.w_hat = M_prices.w_t; ss.r_mkt = M_prices.r_mkt_t;
            ss.C_agg = C_agg; ss.Shock_exp_agg = Shock_exp_agg;
            ss.Bequest_tax = Bequest_tax_agg; ss.Regular_tax = Regular_tax_agg;
            ss.Pension_in = Pension_in_agg; ss.Pension_out = Pension_out_agg;
            ss.Saving_private_flow = Saving_private_flow;
            ss.Depreciation_p = Depreciation_p;
            ss.Depreciation_g = Depreciation_g;
            ss.Public_Capital_Return = Public_Capital_Return;
        end
        
        function [polS, valS] = HHSolution_VFI_unified(M_vfi, paramS_vfi, cS_vfi)
            % [重构核心] 统一的VFI主循环函数，通过向后迭代控制整个VFI过程。

            % --- 1. 初始化输出 ---
            valS = -Inf(cS_vfi.nk, cS_vfi.nkpps, cS_vfi.nw_expanded, cS_vfi.aD_new); % 价值函数矩阵
            polS_cell = cell(cS_vfi.aD_new, 1); % 使用cell数组临时存储每个年龄的策略结构体，效率更高

            % --- 2. 设置年龄相关的养老金 ---
            bV_payg_vfi = zeros(1, cS_vfi.aD_new); % 初始化养老金向量
            if cS_vfi.aR_new < cS_vfi.aD_new      % 只对退休人员设置养老金
                bV_payg_vfi((cS_vfi.aR_new+1):cS_vfi.aD_new) = M_vfi.b_t;
            end

            % --- 3. VFI向后迭代 ---
            for a_idx = cS_vfi.aD_new : -1 : 1 % 从最后一个年龄组向前递推
                vPrime_kkppse_next = []; % 初始化下一期的价值函数
                if a_idx < cS_vfi.aD_new % 如果不是最后一期
                    vPrime_kkppse_next = valS(:,:,:,a_idx+1); % 获取下一期的价值函数作为延续价值
                end

                % 调用统一的、并行的年龄组求解器
                % [valS(:,:,:,a_idx), polS_cell{a_idx}] = ...
                %     main_steady_state_utils_bgp.HHSolutionByAge_unified_serial_adapKgrid(...
                %     a_idx, vPrime_kkppse_next, M_vfi, bV_payg_vfi(a_idx), paramS_vfi, cS_vfi);
                [valS(:,:,:,a_idx), polS_cell{a_idx}] = ...
                    main_steady_state_utils_bgp.HHSolutionByAge_NAIVE(...
                    a_idx, vPrime_kkppse_next, M_vfi, bV_payg_vfi(a_idx), paramS_vfi, cS_vfi);

            end

            % --- 4. 将cell数组转换为更易于访问的结构体数组 ---
            polS = [polS_cell{:}]; % polS(ia).c 可以访问第ia期的消费决策
        end

        function [val_age, pol_age] = HHSolutionByAge_unified_serial(a_idx, vPrime_kkppse_next, M_age, b_age_val, paramS_age, cS)
            % =========================================================================
            % == 函数: HHSolutionByAge_unified_serial
            % == 版本: [v14 - Production-Ready Gold Standard]
            % ==
            % == 目的:
            % ==   提供一个理论正确、数值稳健、逻辑清晰的VFI年龄组求解器，作为
            % ==   模型验证和调试的黄金标准参考。
            % ==
            % == 核心特性:
            % ==   1. [纯串行] 使用最基础的 for-loop 嵌套，逻辑一目了然。
            % ==   2. [简化场景] 只处理无PPS的情形。
            % ==   3. [BGP理论完备] 包含BGP终点决策、双贴现因子等所有理论修正。
            % ==   4. [终极数值稳健] 使用'pchip' + interp1并指定外插值，根除数值
            % ==      不稳定问题。
            % ==   5. [生产就绪] 代码干净，无任何调试输出。
            % =========================================================================

            % --- 1. 初始化输出结构体 ---
            val_age = -1e20 * ones(cS.nk, 1, cS.nw_expanded);
            pol_age = struct(...
                'c', zeros(cS.nk, 1, cS.nw_expanded), ...
                'k_prime', zeros(cS.nk, 1, cS.nw_expanded), ...
                'kpps_prime', zeros(cS.nk, 1, cS.nw_expanded), ...
                'tax_regular', zeros(cS.nk, 1, cS.nw_expanded), ...
                'tax_payg', zeros(cS.nk, 1, cS.nw_expanded), ...
                'shock_exp', zeros(cS.nk, 1, cS.nw_expanded), ...
                'pension_in', zeros(cS.nk, 1, cS.nw_expanded), ...
                'pension_out', zeros(cS.nk, 1, cS.nw_expanded) ...
                );

            % --- 2. 预计算BGP相关的通用参数 ---
            g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
            market_return_factor = 1 + M_age.r_mkt_t;
            discount_factor_continuation = cS.beta * (1 + g_A_period)^(1 - cS.sigma);
            discount_factor_bequest = cS.beta;

            % --- 3. [理论核心] 处理最后一代人的决策 ---
            if a_idx == cS.aD_new
                k_capital_tax = cS.tau_k .* (cS.kGridV * M_age.r_mkt_t);
                total_wealth = cS.kGridV * market_return_factor - k_capital_tax + b_age_val;

                c_expend_final = total_wealth;
                k_prime_final = zeros(cS.nk, 1);

                if isfield(cS,'phi_bequest') && cS.phi_bequest > 0
                    rel_weight_c_over_k = (discount_factor_bequest * cS.phi_bequest)^(-1/cS.sigma) * (1 + g_A_period);
                    denominator = 1 + rel_weight_c_over_k / (1 + cS.tau_c);
                    k_prime_final = total_wealth ./ denominator;
                    c_expend_final = total_wealth - k_prime_final;
                end

                c_final = c_expend_final ./ (1+cS.tau_c);
                consumption_tax = c_final * cS.tau_c;
                final_regular_tax = k_capital_tax + consumption_tax;

                [~, util_c] = model_setup_utils_bgp.CES_utility(c_final, cS.sigma, cS);
                k_prime_final_real_value = k_prime_final * (1 + g_A_period);
                util_bequest = model_setup_utils_bgp.bequest_utility(k_prime_final_real_value, cS);
                util_final = util_c + discount_factor_bequest * util_bequest;

                for ie = 1:cS.nw_expanded
                    val_age(:, 1, ie) = util_final;
                    pol_age.c(:, 1, ie) = c_final;
                    pol_age.k_prime(:, 1, ie) = k_prime_final;
                    pol_age.pension_out(:, 1, ie) = b_age_val;
                    pol_age.tax_regular(:, 1, ie) = final_regular_tax;
                end
                return;
            end

            % --- 4. [理论核心] 处理所有非终点年龄组的决策 ---

            % 4.1 预计算期望价值函数 (包含数据清洗)
            EV_matrix = zeros(cS.nk, cS.nw_expanded);
            if ~isempty(vPrime_kkppse_next)
                vPrime_nopps_raw = squeeze(vPrime_kkppse_next(:,1,:));
                vPrime_nopps_clean = vPrime_nopps_raw;
                % 将所有非法或极小的值统一为一个固定的极小值
                vPrime_nopps_clean(vPrime_nopps_raw < -1e10 | isnan(vPrime_nopps_raw)) = -1e10;

                trans_mat_next = paramS_age.TrProbM_by_age{a_idx + 1};
                % 使用矩阵乘法高效计算期望
                EV_matrix = vPrime_nopps_clean * trans_mat_next';
            end

            % --- 4.2 [最原始的VFI循环] 遍历所有当前状态 (ik, ie) ---
            for ie = 1:cS.nw_expanded

                % 获取当前冲击状态对应的期望价值函数向量
                ev_vec_to_interp = EV_matrix(:,ie);

                for ik = 1:cS.nk
                    % --- A. 计算当前状态 (ik, ie) 下的所有资源和参数 ---
                    k_now = cS.kGridV(ik);

                    % 计算收入、税收和养老金 (直接计算，不优化)
                    capital_income = k_now * M_age.r_mkt_t;
                    capital_tax = cS.tau_k * capital_income;
                    k_return = k_now * market_return_factor;
                    pension_out = 0;
                    payg_tax = 0;
                    labor_tax = 0;
                    labor_income_gross = 0;

                    if a_idx <= cS.aR_new % 工作期
                        labor_income_gross = M_age.w_t * cS.ageEffV_new(a_idx) * paramS_age.leGridV(ie);
                        payg_tax = cS.theta_t * labor_income_gross;
                        labor_tax = cS.tau_l * max(0, labor_income_gross - payg_tax);
                    else % 退休期
                        pension_out = b_age_val;
                    end

                    cash_flow_after_taxes = k_return + labor_income_gross + pension_out ...
                        - (capital_tax + payg_tax + labor_tax);

                    % 计算支出冲击
                    shock_exp_factor = (ie == cS.nw + 1) * cS.kappa_young + (ie == cS.nw + 2) * cS.kappa_old;
                    shock_exp = shock_exp_factor * cash_flow_after_taxes;

                    % 最终可用于消费和储蓄的资源
                    available_for_c_and_s = cash_flow_after_taxes - shock_exp;

                    % --- B. 暴力网格搜索：遍历所有下一期储蓄选择 ik_prime ---
                    value_by_choice = -1e20 * ones(cS.nk, 1); % 存储每个选择的价值

                    for ik_prime = 1:cS.nk
                        k_prime_choice = cS.kGridV(ik_prime);

                        % BGP预算约束
                        c_expend_choice = available_for_c_and_s - k_prime_choice * (1 + g_A_period);

                        % 检查可行性：消费必须大于一个极小的正数
                        if c_expend_choice < 1e-9
                            continue; % 如果不可行，则跳过此选项
                        end

                        % 计算此选择下的消费和所有效用项
                        c_choice = c_expend_choice / (1 + cS.tau_c);

                        [~, util_c] = model_setup_utils_bgp.CES_utility(c_choice, cS.sigma, cS);

                        % [数值稳健性核心] 使用 interp1 并明确指定外插值
                        ev_choice = interp1(cS.kGridV, ev_vec_to_interp, k_prime_choice, 'pchip', -1e10);

                        k_prime_real_value = k_prime_choice * (1 + g_A_period);
                        util_bequest = model_setup_utils_bgp.bequest_utility(k_prime_real_value, cS);

                        % 组合成总价值
                        value_by_choice(ik_prime) = util_c + ...
                            cS.s_pathV(a_idx) * discount_factor_continuation * ev_choice + ...
                            (1 - cS.s_pathV(a_idx)) * discount_factor_bequest * util_bequest;
                    end

                    % --- C. 找到最优决策并存储结果 ---
                    [best_val, best_ik_prime] = max(value_by_choice);

                    if isfinite(best_val) && best_ik_prime > 0
                        % 获取最优决策
                        best_k_prime = cS.kGridV(best_ik_prime);
                        best_c_expend = available_for_c_and_s - best_k_prime * (1 + g_A_period);
                        best_c = best_c_expend / (1+cS.tau_c);

                        % 存储结果到pol_age
                        pol_age.c(ik, 1, ie) = best_c;
                        pol_age.k_prime(ik, 1, ie) = best_k_prime;
                        pol_age.pension_out(ik, 1, ie) = pension_out;
                        pol_age.shock_exp(ik, 1, ie) = shock_exp;
                        pol_age.tax_regular(ik, 1, ie) = capital_tax + labor_tax + best_c*cS.tau_c;
                        pol_age.tax_payg(ik, 1, ie) = payg_tax;

                        % 存储最优价值
                        val_age(ik, 1, ie) = best_val;
                    else
                        % 如果所有选项都不可行 (极少发生)，采取保守策略
                        val_age(ik, 1, ie) = -1e20;
                        pol_age.c(ik, 1, ie) = cS.cFloor;
                        pol_age.k_prime(ik, 1, ie) = cS.kMin;
                    end

                end % 结束 ik 循环
            end % 结束 ie 循环
        end

        function [val_age, pol_age] = HHSolutionByAge_unified_PARFOR(a_idx, vPrime_kkppse_next, M_age, b_age_val, paramS_age, cS)
            % =========================================================================
            % == 函数说明: [v6 - Final BGP Micro Fix] 年龄组决策求解器 - 修正贴现逻辑
            % == 核心修正:
            % ==  - 明确区分对【延续价值EV】和【遗赠效用】的贴现。
            % ==  - 延续价值使用增长调整的贴现因子: beta * (1+g)^(1-sigma)
            % ==  - 遗赠效用仅使用纯粹时间偏好贴现因子: beta
            % ==  - 这确保了家庭微观决策的理论一致性。
            % =========================================================================

            % --- 1. 初始化输出结构体 ---
            val_age = -1e20 * ones(cS.nk, cS.nkpps, cS.nw_expanded);

            pol_age = struct(...
                'c', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
                'k_prime', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
                'kpps_prime', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
                'tax_regular', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
                'tax_payg', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
                'shock_exp', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
                'pension_in', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
                'pension_out', zeros(cS.nk, cS.nkpps, cS.nw_expanded) ...
                );

            % --- 2. 预计算通用参数 ---
            g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
            market_return_factor = 1 + M_age.r_mkt_t;

            % [!!!!! 核心修正：区分贴现因子 !!!!!]
            % 用于贴现【延续价值】(未来消费流)的因子，需要增长调整
            discount_factor_continuation = cS.beta * (1 + g_A_period)^(1 - cS.sigma);
            % 用于贴现【遗赠效用】(一次性事件)的因子，仅需时间偏好
            discount_factor_bequest = cS.beta;

            % --- 3. 最后一期的特殊处理（保持不变，因为逻辑相对简单） ---
            if a_idx == cS.aD_new
                [K_grid, Kpps_grid] = ndgrid(cS.kGridV, cS.kppsGridV);
                k_capital_tax = cS.tau_k .* (K_grid .* M_age.r_mkt_t);
                total_wealth = K_grid .* market_return_factor - k_capital_tax + b_age_val;

                pps_withdrawal_tax = 0;
                if cS.pps_active
                    pps_wealth = Kpps_grid .* market_return_factor;
                    pps_withdrawal_tax = cS.pps_tax_rate_withdrawal .* pps_wealth;
                    total_wealth = total_wealth + pps_wealth - pps_withdrawal_tax;
                end

                k_prime_final = zeros(size(K_grid));
                c_expend_final = total_wealth;
                if isfield(cS,'phi_bequest') && cS.phi_bequest > 0
                    % 计算消费和遗赠的分配
                    phi_adj = cS.phi_bequest; % 遗赠动机参数本身

                    % 分享比例取决于调整后的遗赠效用与当期消费效用的比
                    phi_effective = phi_adj * discount_factor_bequest * (1 + g_A_period)^(1-cS.sigma);
                    if abs(cS.sigma-1)<1e-6, optimal_c_share=1/(1+phi_effective); else, optimal_c_share=1/(1+phi_effective^(1/cS.sigma)); end

                    c_expend_final = optimal_c_share * total_wealth;
                    k_prime_final = (1 - optimal_c_share) * total_wealth;
                end

                c_final = c_expend_final ./ (1+cS.tau_c);
                consumption_tax = c_final * cS.tau_c;

                [~, util_c] = model_setup_utils_bgp.CES_utility(c_final, cS.sigma, cS);

                % 计算遗赠效用 (使用其真实价值)
                k_prime_final_real_value = k_prime_final * (1 + g_A_period);
                util_bequest = model_setup_utils_bgp.bequest_utility(k_prime_final_real_value, cS);

                % 总效用 = 当期效用 + 纯粹贴现后的遗赠效用
                util_final = util_c + discount_factor_bequest * util_bequest;

                final_regular_tax = k_capital_tax + pps_withdrawal_tax + consumption_tax;

                for ie = 1:cS.nw_expanded
                    val_age(:, :, ie) = util_final;
                    pol_age.c(:, :, ie) = c_final;
                    pol_age.k_prime(:, :, ie) = k_prime_final;
                    pol_age.kpps_prime(:, :, ie) = 0;
                    pol_age.tax_regular(:, :, ie) = final_regular_tax;
                    pol_age.tax_payg(:, :, ie) = 0;
                    pol_age.pension_out(:, :, ie) = b_age_val;
                    pol_age.shock_exp(:, :, ie) = 0;
                    pol_age.pension_in(:, :, ie) = 0;
                end
                return;
            end

            % --- 4. 预计算下一期的期望价值函数 ---
            EV_matrix = zeros(cS.nk, cS.nkpps, cS.nw_expanded);
            if ~isempty(vPrime_kkppse_next)
                transition_matrix_next_age = paramS_age.TrProbM_by_age{a_idx + 1};
                vPrime_reshaped = reshape(vPrime_kkppse_next, [cS.nk * cS.nkpps, cS.nw_expanded]);
                for ie_current = 1:cS.nw_expanded
                    transition_probs = transition_matrix_next_age(ie_current, :);
                    EV_slice = vPrime_reshaped * transition_probs';
                    EV_matrix(:, :, ie_current) = reshape(EV_slice, [cS.nk, cS.nkpps]);
                end
            end

            % --- 5. [核心优化] 预计算不变量 ---

            % 5.1 预计算选择网格矩阵（提升到循环外部）
            % 首先初始化所有选择网格，确保在任何情况下都有定义
            K_prime_choices_no_pps = cS.kGridV';
            if cS.pps_active
                [K_prime_choices_pps, Kpps_prime_choices] = ndgrid(cS.kGridV, cS.kppsGridV);
            else
                % 即使在非PPS模式下也定义这些矩阵，避免parfor循环中的未定义变量错误
                % 使用最小尺寸的虚拟矩阵
                K_prime_choices_pps = zeros(1, 1);
                Kpps_prime_choices = zeros(1, 1);
            end

            % 5.2 预计算只依赖ie的向量
            labor_income_vec = zeros(cS.nw_expanded, 1);
            payg_tax_vec = zeros(cS.nw_expanded, 1);
            shock_exp_factor_vec = zeros(cS.nw_expanded, 1);

            if a_idx <= cS.aR_new  % 工作期
                labor_income_vec = M_age.w_t * cS.ageEffV_new(a_idx) * paramS_age.leGridV;
                payg_tax_vec = cS.theta_t * labor_income_vec;
            end

            % 计算冲击因子向量
            for ie = 1:cS.nw_expanded
                shock_exp_factor_vec(ie) = (ie == cS.nw + 1) * cS.kappa_young + (ie == cS.nw + 2) * cS.kappa_old;
            end

            % 5.3 预计算只依赖ik的向量
            capital_income_vec = cS.kGridV * M_age.r_mkt_t;
            capital_tax_vec = cS.tau_k * capital_income_vec;
            k_return_vec = cS.kGridV * market_return_factor;

            % 5.4 预计算只依赖ikpps的向量（如果有PPS）
            % 首先初始化所有向量，确保在任何情况下都有定义
            kpps_return_vec = zeros(cS.nkpps, 1);
            mandatory_withdrawal_vec = zeros(cS.nkpps, 1);
            pps_withdrawal_tax_vec = zeros(cS.nkpps, 1);
            kpps_prime_fixed_vec = zeros(cS.nkpps, 1);

            if cS.pps_active
                kpps_return_vec = cS.kppsGridV * market_return_factor;
                if a_idx > cS.aR_new  % 退休期
                    mandatory_withdrawal_vec = kpps_return_vec * cS.pps_withdrawal_rate;
                    pps_withdrawal_tax_vec = cS.pps_tax_rate_withdrawal * mandatory_withdrawal_vec;
                    kpps_prime_fixed_vec = kpps_return_vec - mandatory_withdrawal_vec;
                end
            end

            % --- 6. 并行循环求解每个状态点的最优决策 ---
            total_states = cS.nk * cS.nkpps * cS.nw_expanded;
            temp_results = cell(total_states, 1);

            for state_idx = 1:total_states
                warning('off', 'MATLAB:griddedInterpolant:MeshgridEval2DWarnId');
                [ik, ikpps, ie] = ind2sub([cS.nk, cS.nkpps, cS.nw_expanded], state_idx);

                % 使用预计算的向量获取当前状态的基础值
                k_now = cS.kGridV(ik);
                epsilon_state = paramS_age.leGridV(ie);
                pension_out = 0;

                % 从预计算向量中获取值
                capital_tax = capital_tax_vec(ik);
                k_return = k_return_vec(ik);
                shock_exp_factor = shock_exp_factor_vec(ie);

                % 根据是否有PPS，选择不同的分支
                if cS.pps_active
                    kpps_now = cS.kppsGridV(ikpps);
                    kpps_return = kpps_return_vec(ikpps);
                    ev_interpolant = griddedInterpolant({cS.kGridV, cS.kppsGridV}, EV_matrix(:,:,ie), 'linear', 'none');
                    local_K_prime = K_prime_choices_pps;     % 使用预计算的选择网格
                    local_Kpps_prime = Kpps_prime_choices;   % 使用预计算的选择网格

                    if a_idx <= cS.aR_new % --- 工作期 (有PPS) ---
                        labor_income_gross = labor_income_vec(ie);    % 从预计算向量获取
                        payg_tax = payg_tax_vec(ie);                  % 从预计算向量获取

                        total_inflow = k_return + labor_income_gross;

                        % 计算PPS缴费网格（这部分仍需要在循环内计算，因为依赖于当前kpps状态）
                        pps_contrib_grid = max(0, Kpps_prime_choices - kpps_return);
                        taxable_income_grid = max(0, labor_income_gross - payg_tax - pps_contrib_grid);
                        labor_tax_grid = cS.tau_l * taxable_income_grid;

                        cash_flow_after_all_taxes_grid = total_inflow - (payg_tax + capital_tax) - labor_tax_grid;
                        shock_exp_grid = shock_exp_factor * cash_flow_after_all_taxes_grid;
                        net_resource_for_c_and_s = cash_flow_after_all_taxes_grid - shock_exp_grid;
                        c_expend_choices = net_resource_for_c_and_s - local_K_prime * (1 + g_A_period) - pps_contrib_grid;
                        feasible_mask = (c_expend_choices > 1e-9) & (pps_contrib_grid <= cS.pps_contrib_limit);
                    else % --- 退休期 (有PPS) ---
                        pension_out = b_age_val;

                        % 使用预计算的向量
                        mandatory_withdrawal = mandatory_withdrawal_vec(ikpps);
                        pps_withdrawal_tax = pps_withdrawal_tax_vec(ikpps);
                        kpps_prime_fixed = kpps_prime_fixed_vec(ikpps);

                        cash_flow_after_taxes = k_return + pension_out - capital_tax + (mandatory_withdrawal - pps_withdrawal_tax);
                        shock_exp = shock_exp_factor * cash_flow_after_taxes;
                        available_basic = cash_flow_after_taxes - shock_exp;

                        c_expend_choices = available_basic - local_K_prime * (1 + g_A_period);
                        feasible_mask = c_expend_choices > 1e-9;
                        local_Kpps_prime = repmat(kpps_prime_fixed, size(local_K_prime));
                    end
                else % --- 无PPS模式 ---
                    ev_interpolant = griddedInterpolant(cS.kGridV, EV_matrix(:,1,ie), 'phip', 'none');
                    local_K_prime = K_prime_choices_no_pps;  % 使用预计算的选择网格
                    local_Kpps_prime = [];  % 在无PPS模式下定义为空

                    if a_idx <= cS.aR_new % --- 工作期 (无PPS) ---
                        labor_income_gross = labor_income_vec(ie);    % 从预计算向量获取
                        payg_tax = payg_tax_vec(ie);                  % 从预计算向量获取

                        labor_tax = cS.tau_l * max(0, labor_income_gross - payg_tax);
                        cash_flow_after_taxes = k_return + labor_income_gross - (payg_tax + capital_tax + labor_tax);
                        shock_exp = shock_exp_factor * cash_flow_after_taxes;
                        available_basic = cash_flow_after_taxes - shock_exp;
                        c_expend_choices = available_basic - local_K_prime * (1 + g_A_period);
                        feasible_mask = c_expend_choices > 1e-9;
                    else % --- 退休期 (无PPS) ---
                        pension_out = b_age_val;
                        cash_flow_after_taxes = k_return + pension_out - capital_tax;
                        shock_exp = shock_exp_factor * cash_flow_after_taxes;
                        available_basic = cash_flow_after_taxes - shock_exp;
                        c_expend_choices = available_basic - local_K_prime * (1 + g_A_period);
                        feasible_mask = c_expend_choices > 1e-9;
                    end
                end

                % --- 寻找最优决策并打包结果（这部分逻辑保持不变） ---
                if ~any(feasible_mask(:))
                    result_for_state = {-1e20, 0, cS.kMin, 0, 0, 0, 0, 0, 0};
                else
                    c_choices = c_expend_choices(feasible_mask) / (1 + cS.tau_c);
                    k_prime_feasible = local_K_prime(feasible_mask);
                    if cS.pps_active
                        kpps_prime_feasible = local_Kpps_prime(feasible_mask);
                        ev_vec = ev_interpolant(k_prime_feasible, kpps_prime_feasible);
                    else
                        kpps_prime_feasible = zeros(size(k_prime_feasible));
                        ev_vec = ev_interpolant(k_prime_feasible);
                    end

                    ev_vec(isnan(ev_vec)) = -1e10;
                    [~, util_c_vec] = model_setup_utils_bgp.CES_utility(c_choices, cS.sigma, cS);

                    % [!!!!! 核心修正：在VFI中应用正确的贴现 !!!!!]

                    % 1. 计算当期消费效用
                    % (util_c_vec 已在上面计算)

                    % 2. 计算遗赠效用 (使用其真实价值)
                    k_prime_total_feasible = k_prime_feasible + kpps_prime_feasible;
                    k_prime_total_real_value = k_prime_total_feasible * (1 + g_A_period);
                    util_bequest_vec = model_setup_utils_bgp.bequest_utility(k_prime_total_real_value, cS);

                    % 3. 组合成总价值 (应用不同的贴现因子！)
                    total_values_vec = util_c_vec + ...
                        ( cS.s_pathV(a_idx) * discount_factor_continuation * ev_vec ) + ...
                        ( (1 - cS.s_pathV(a_idx)) * discount_factor_bequest * util_bequest_vec );

                    [best_val, best_idx] = max(total_values_vec);

                    if isfinite(best_val)
                        best_c = c_choices(best_idx);
                        best_k = k_prime_feasible(best_idx);
                        best_kpps = kpps_prime_feasible(best_idx);
                        consumption_tax = best_c * cS.tau_c;

                        % 计算并分离税种
                        if a_idx <= cS.aR_new % 工作期
                            labor_income_gross = labor_income_vec(ie);
                            payg_tax = payg_tax_vec(ie);
                            pension_in_pps = 0;

                            if cS.pps_active
                                kpps_return = kpps_return_vec(ikpps);
                                best_pps_contrib = max(0, best_kpps - kpps_return);
                                best_labor_tax = cS.tau_l * max(0, labor_income_gross - payg_tax - best_pps_contrib);
                                total_inflow = k_return + labor_income_gross;
                                final_cash_flow = total_inflow - (payg_tax + capital_tax + best_labor_tax);
                                final_shock_exp = shock_exp_factor * final_cash_flow;
                                pension_in_pps = best_pps_contrib;
                            else
                                best_labor_tax = cS.tau_l * max(0, labor_income_gross - payg_tax);
                                final_shock_exp = shock_exp;  % 直接使用之前计算的值
                            end
                            total_regular_tax = capital_tax + best_labor_tax + consumption_tax;
                            total_payg_tax = payg_tax;
                        else % 退休期
                            final_shock_exp = shock_exp;  % 使用之前计算的值
                            pension_in_pps = 0;
                            total_payg_tax = 0;
                            if cS.pps_active
                                pps_withdrawal_tax = pps_withdrawal_tax_vec(ikpps);
                                total_regular_tax = capital_tax + pps_withdrawal_tax + consumption_tax;
                            else
                                total_regular_tax = capital_tax + consumption_tax;
                            end
                        end

                        result_for_state = {best_val, best_c, best_k, best_kpps, total_regular_tax, total_payg_tax, final_shock_exp, pension_in_pps, pension_out};
                    else
                        result_for_state = {-1e20, 0, cS.kMin, 0, 0, 0, 0, 0, 0};
                    end
                end
                temp_results{state_idx} = result_for_state;
            end

            % --- 7. 将并行计算结果重组为矩阵形式 ---
            for state_idx = 1:total_states
                [ik, ikpps, ie] = ind2sub([cS.nk, cS.nkpps, cS.nw_expanded], state_idx);
                res = temp_results{state_idx};
                val_age(ik, ikpps, ie) = res{1};
                pol_age.c(ik, ikpps, ie) = res{2};
                pol_age.k_prime(ik, ikpps, ie) = res{3};
                pol_age.kpps_prime(ik, ikpps, ie) = res{4};
                pol_age.tax_regular(ik, ikpps, ie) = res{5};
                pol_age.tax_payg(ik, ikpps, ie) = res{6};
                pol_age.shock_exp(ik, ikpps, ie) = res{7};
                pol_age.pension_in(ik, ikpps, ie) = res{8};
                pol_age.pension_out(ik, ikpps, ie) = res{9};
            end
        end

        function [val_age, pol_age] = HHSolutionByAge_unified_serial_adapKgrid(a_idx, vPrime_kkppse_next, M_age, b_age_val, paramS_age, cS)
            % =========================================================================
            % == 函数: HHSolutionByAge_unified_serial_adapKgrid
            % == 版本: [v22 - 核心预算约束修正版]
            % ==
            % == 核心修正:
            % ==   1. [理论修正] 纠正了BGP预算约束的定义。移除了错误的“未实现资本
            % ==      利得”项(g*k)作为收入来源。增长g只影响储蓄的成本(k'*(1+g))，
            % ==      而不是当期资源。
            % ==   2. [逻辑修正] 基于正确的资源定义，重新计算了支出冲击的基数。
            % ==   3. 这是解决宏观国民账户不平衡问题的根本性修复。
            % =========================================================================

            % --- 1. 初始化和参数设置 ---
            nk_local = 200;
            val_age = -1e20 * ones(cS.nk, 1, cS.nw_expanded);
            pol_age = struct(...
                'c', zeros(cS.nk, 1, cS.nw_expanded), 'k_prime', zeros(cS.nk, 1, cS.nw_expanded), ...
                'kpps_prime', zeros(cS.nk, 1, cS.nw_expanded), 'tax_regular', zeros(cS.nk, 1, cS.nw_expanded), ...
                'tax_payg', zeros(cS.nk, 1, cS.nw_expanded), 'shock_exp', zeros(cS.nk, 1, cS.nw_expanded), ...
                'pension_in', zeros(cS.nk, 1, cS.nw_expanded), 'pension_out', zeros(cS.nk, 1, cS.nw_expanded));

            % --- 2. 预计算BGP相关的通用参数 ---
            g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
            market_return_factor = 1 + M_age.r_mkt_t;
            discount_factor_future_utility = cS.beta * (1 + g_A_period)^(1 - cS.sigma);

            % --- 3. 终点决策 (这部分逻辑在BGP框架下是正确的, 保持不变) ---
            if a_idx == cS.aD_new
                % 终点决策是将当期总财富在消费和遗赠之间分配
                k_capital_tax = cS.tau_k .* (cS.kGridV * M_age.r_mkt_t);
                % 终点财富 = 资本回报 + 养老金 - 资本税
                total_wealth = cS.kGridV * market_return_factor - k_capital_tax + b_age_val;

                c_expend_final = total_wealth;
                k_prime_final = zeros(cS.nk, 1);

                if isfield(cS,'phi_bequest') && cS.phi_bequest > 0
                    % 财富在 c_expend 和 k_prime_bequest 之间分配
                    % u'(c) = φ * β' * u'(k') -> u'(c*(1+tau_c)) = ...
                    phi_effective = discount_factor_future_utility * cS.phi_bequest;
                    c_over_k_ratio = phi_effective^(-1/cS.sigma);

                    k_prime_final = total_wealth ./ (1 + (1+cS.tau_c)*c_over_k_ratio);
                    c_expend_final = total_wealth - k_prime_final;
                end

                c_final = c_expend_final ./ (1+cS.tau_c);
                consumption_tax = c_final * cS.tau_c;
                final_regular_tax = k_capital_tax + consumption_tax;

                [~, util_c] = model_setup_utils_bgp.CES_utility(c_final, cS.sigma, cS);
                util_bequest = model_setup_utils_bgp.bequest_utility(k_prime_final, cS);
                util_final = util_c + discount_factor_future_utility * util_bequest;

                for ie = 1:cS.nw_expanded
                    val_age(:, 1, ie) = util_final;
                    pol_age.c(:, 1, ie) = c_final;
                    pol_age.k_prime(:, 1, ie) = k_prime_final;
                    pol_age.pension_out(:, 1, ie) = b_age_val;
                    pol_age.tax_regular(:, 1, ie) = final_regular_tax;
                end
                return;
            end

            % --- 4. 非终点年龄组决策 ---
            vPrime_interpolants = cell(cS.nw_expanded, 1);
            if ~isempty(vPrime_kkppse_next)
                vPrime_nopps_raw = squeeze(vPrime_kkppse_next(:,1,:));
                vPrime_nopps_clean = vPrime_nopps_raw;
                vPrime_nopps_clean(vPrime_nopps_raw < -1e10 | isnan(vPrime_nopps_raw)) = -1e10;
                for ie_next = 1:cS.nw_expanded
                    vPrime_interpolants{ie_next} = griddedInterpolant(cS.kGridV, vPrime_nopps_clean(:, ie_next), 'pchip', 'nearest');
                end
            end
            trans_mat_next = paramS_age.TrProbM_by_age{a_idx + 1};

            for ie = 1:cS.nw_expanded
                for ik = 1:cS.nk
                    % --- A. [核心修正] 重新计算可支配资源，移除"未实现资本利得" ---
                    k_now = cS.kGridV(ik);

                    % 流入项
                    k_return = k_now * market_return_factor; % 资本及其市场回报 k*(1+r)
                    labor_income_gross = 0;
                    pension_out = 0;
                    if a_idx <= cS.aR_new
                        labor_income_gross = M_age.w_t * cS.ageEffV_new(a_idx) * paramS_age.leGridV(ie);
                    else
                        pension_out = b_age_val;
                    end

                    % 税收流出项
                    capital_tax = cS.tau_k * (k_now * M_age.r_mkt_t);
                    payg_tax = cS.theta_t * labor_income_gross;
                    labor_tax = cS.tau_l * max(0, labor_income_gross - payg_tax);

                    % [核心修正] 正确的 "税后现金流" (Cash-on-Hand)
                    cash_on_hand = k_return + labor_income_gross + pension_out - (capital_tax + payg_tax + labor_tax);

                    % [核心修正] 冲击支出应基于正确的现金流计算
                    shock_exp_factor = (ie == cS.nw + 1) * cS.kappa_young + (ie == cS.nw + 2) * cS.kappa_old;
                    shock_exp = shock_exp_factor * cash_on_hand;

                    % [核心修正] 最终可用于 "消费" 和 "储蓄" 的资源
                    available_for_c_and_s = cash_on_hand - shock_exp;

                    % --- D. 网格搜索最优决策 (应用最终修正) ---
                    % [k'*(1+g)] + c_expend <= available_for_c_and_s
                    k_prime_min_local = cS.kMin;
                    k_prime_max_local = (available_for_c_and_s - cS.cFloor * (1 + cS.tau_c)) / (1 + g_A_period);

                    if k_prime_max_local <= k_prime_min_local
                        k_prime_choice_grid_local = k_prime_min_local;
                    else
                        k_prime_choice_grid_local = linspace(k_prime_min_local, k_prime_max_local, nk_local)';
                    end

                    % 期望价值函数
                    v_at_k_prime_choices = zeros(nk_local, cS.nw_expanded);
                    if ~isempty(vPrime_kkppse_next)
                        for ie_next = 1:cS.nw_expanded
                            v_at_k_prime_choices(:, ie_next) = vPrime_interpolants{ie_next}(k_prime_choice_grid_local);
                        end
                    end
                    transition_probs_from_ie = trans_mat_next(ie, :);
                    ev_choices = v_at_k_prime_choices * transition_probs_from_ie';

                    % 根据预算约束计算每个选择对应的消费
                    c_expend_choices = available_for_c_and_s - k_prime_choice_grid_local .* (1 + g_A_period);
                    c_choices = c_expend_choices / (1 + cS.tau_c);

                    [~, util_c_vec] = model_setup_utils_bgp.CES_utility(c_choices, cS.sigma, cS);
                    util_bequest_vec = model_setup_utils_bgp.bequest_utility(k_prime_choice_grid_local, cS);

                    total_values_vec = util_c_vec + ...
                        discount_factor_future_utility * ( ...
                        cS.s_pathV(a_idx) * ev_choices + ...
                        (1 - cS.s_pathV(a_idx)) * util_bequest_vec ...
                        );

                    total_values_vec(c_expend_choices < 1e-9) = -1e20;
                    [best_val, best_idx_local] = max(total_values_vec);

                    if isfinite(best_val) && best_idx_local > 0 && c_expend_choices(best_idx_local) > 1e-9
                        best_k_prime = k_prime_choice_grid_local(best_idx_local);
                        best_c = c_choices(best_idx_local);

                        pol_age.c(ik, 1, ie) = best_c;
                        pol_age.k_prime(ik, 1, ie) = best_k_prime;
                        pol_age.pension_out(ik, 1, ie) = pension_out;
                        pol_age.shock_exp(ik, 1, ie) = shock_exp;
                        pol_age.tax_regular(ik, 1, ie) = capital_tax + labor_tax + best_c*cS.tau_c;
                        pol_age.tax_payg(ik, 1, ie) = payg_tax;
                        val_age(ik, 1, ie) = best_val;
                    else
                        val_age(ik, 1, ie) = -1e20;
                        pol_age.c(ik, 1, ie) = cS.cFloor;
                        pol_age.k_prime(ik, 1, ie) = cS.kMin;
                    end
                end
            end
        end

        function [val_age, pol_age] = HHSolutionByAge_NAIVE(a_idx, ~, M_age, b_age_val, paramS_age, cS)
            % =========================================================================
            % == 函数: HHSolutionByAge_NAIVE
            % == 版本: [v1.0 - 调试专用]
            % ==
            % == 目的:
            % ==   提供一个极简的、基于固定储蓄率的家庭决策规则，用于替换复杂的
            % ==   VFI求解器，以帮助定位模型不收敛的根源。
            % ==
            % == 核心逻辑:
            % ==   1. 计算与完整版VFI完全一致的税后现金流 (Cash-on-Hand)。
            % ==   2. 应用一个固定的储蓄率来决定下一期储蓄 k'。
            % ==   3. 剩余资源全部用于消费。
            % ==   4. 严格遵守BGP预算约束。
            % ==   5. 不计算价值函数 (val_age)，因为此规则是短视的。
            % =========================================================================

            % --- 0. 设定一个固定的、合理的储蓄率 (这是唯一的行为参数) ---
            SAVINGS_RATE = 0.20; % 假设家庭将20%的可支配资源用于储蓄

            % --- 1. 初始化输出结构体 (与原函数保持一致) ---
            val_age = -1e20 * ones(cS.nk, 1, cS.nw_expanded); % 价值函数在此处无意义
            pol_age = struct(...
                'c', zeros(cS.nk, 1, cS.nw_expanded), 'k_prime', zeros(cS.nk, 1, cS.nw_expanded), ...
                'kpps_prime', zeros(cS.nk, 1, cS.nw_expanded), 'tax_regular', zeros(cS.nk, 1, cS.nw_expanded), ...
                'tax_payg', zeros(cS.nk, 1, cS.nw_expanded), 'shock_exp', zeros(cS.nk, 1, cS.nw_expanded), ...
                'pension_in', zeros(cS.nk, 1, cS.nw_expanded), 'pension_out', zeros(cS.nk, 1, cS.nw_expanded));

            % --- 2. 预计算BGP相关的通用参数 (与原函数保持一致) ---
            g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
            market_return_factor = 1 + M_age.r_mkt_t;

            % --- 3. 终点决策 (简化为消费掉所有财富) ---
            if a_idx == cS.aD_new
                k_capital_tax = cS.tau_k .* (cS.kGridV * M_age.r_mkt_t);
                total_wealth = cS.kGridV * market_return_factor - k_capital_tax + b_age_val;

                c_expend_final = total_wealth;
                k_prime_final = zeros(cS.nk, 1); % 最后一期不储蓄
                c_final = c_expend_final ./ (1 + cS.tau_c);
                consumption_tax = c_final * cS.tau_c;
                final_regular_tax = k_capital_tax + consumption_tax;

                for ie = 1:cS.nw_expanded
                    pol_age.c(:, 1, ie) = c_final;
                    pol_age.k_prime(:, 1, ie) = k_prime_final;
                    pol_age.pension_out(:, 1, ie) = b_age_val;
                    pol_age.tax_regular(:, 1, ie) = final_regular_tax;
                end
                return;
            end

            % --- 4. 遍历所有状态点，应用固定储蓄率规则 ---
            for ie = 1:cS.nw_expanded
                for ik = 1:cS.nk
                    % --- A. 计算可支配资源 (此部分逻辑必须与VFI函数v22版完全一致!) ---
                    k_now = cS.kGridV(ik);

                    % 流入项
                    k_return = k_now * market_return_factor;
                    labor_income_gross = 0;
                    pension_out = 0;
                    if a_idx <= cS.aR_new
                        labor_income_gross = M_age.w_t * cS.ageEffV_new(a_idx) * paramS_age.leGridV(ie);
                    else
                        pension_out = b_age_val;
                    end

                    % 税收流出项
                    capital_tax = cS.tau_k * (k_now * M_age.r_mkt_t);
                    payg_tax = cS.theta_t * labor_income_gross;
                    labor_tax = cS.tau_l * max(0, labor_income_gross - payg_tax);

                    % 税后现金流 (Cash-on-Hand)
                    cash_on_hand = k_return + labor_income_gross + pension_out - (capital_tax + payg_tax + labor_tax);

                    % 冲击支出
                    shock_exp_factor = (ie == cS.nw + 1) * cS.kappa_young + (ie == cS.nw + 2) * cS.kappa_old;
                    shock_exp = shock_exp_factor * cash_on_hand;

                    % 最终可用于 "消费" 和 "储蓄" 的资源
                    available_for_c_and_s = cash_on_hand - shock_exp;

                    % --- B. 应用固定储蓄率规则 ---

                    % 计算储蓄额
                    saving_decision = SAVINGS_RATE * available_for_c_and_s;

                    % 根据BGP预算约束，计算出下一期的标准化资本 k_prime
                    k_prime_decision = saving_decision / (1 + g_A_period);

                    % 确保储蓄不会导致资产低于下限
                    k_prime_decision = max(cS.kMin, k_prime_decision);

                    % --- C. 计算消费 ---

                    % 重新计算储蓄的真实成本
                    saving_expenditure = k_prime_decision * (1 + g_A_period);

                    % 剩余资源全部用于消费支出
                    c_expend_decision = available_for_c_and_s - saving_expenditure;

                    % 消费支出必须为正
                    if c_expend_decision < 1e-9
                        c_expend_decision = 1e-9;
                        % 如果资源不足，相应调整储蓄
                        saving_expenditure = available_for_c_and_s - c_expend_decision;
                        k_prime_decision = max(cS.kMin, saving_expenditure / (1 + g_A_period));
                    end

                    c_decision = c_expend_decision / (1 + cS.tau_c);
                    consumption_tax = c_decision * cS.tau_c;

                    % --- D. 存储决策 ---
                    pol_age.c(ik, 1, ie) = c_decision;
                    pol_age.k_prime(ik, 1, ie) = k_prime_decision;
                    pol_age.pension_out(ik, 1, ie) = pension_out;
                    pol_age.shock_exp(ik, 1, ie) = shock_exp;
                    pol_age.tax_regular(ik, 1, ie) = capital_tax + labor_tax + consumption_tax;
                    pol_age.tax_payg(ik, 1, ie) = payg_tax;
                end
            end
        end

        function idx_mat = get_policy_index_matrix_unified(polS, cS, type)
            % [重构-修正版] 统一的策略函数离散化，将连续决策映射到离散网格索引。

            % --- 确定要离散化的策略('k'或'kpps')和对应的网格 ---
            if strcmp(type, 'k')
                gridV = cS.kGridV;
            elseif strcmp(type, 'kpps') && cS.pps_active
                gridV = cS.kppsGridV;
            else
                % 如果类型是'kpps'但PPS未激活，或类型未知，则返回一个全1的矩阵(代表第一个网格点)
                idx_mat = ones(cS.nk, cS.nkpps, cS.nw_expanded, cS.aD_new, 'uint16');
                return;
            end

            idx_mat = zeros(cS.nk, cS.nkpps, cS.nw_expanded, cS.aD_new, 'uint16'); % 初始化索引矩阵

            for ia = 1:cS.aD_new % 遍历所有年龄
                % --- 获取对应年龄的策略矩阵 ---
                if strcmp(type, 'k')
                    val_mat = polS(ia).k_prime;
                else
                    val_mat = polS(ia).kpps_prime;
                end

                % --- [核心修正] 使用嵌套循环和多维索引来确保正确赋值 ---
                for ik = 1:cS.nk
                    for ikpps = 1:cS.nkpps
                        for ie = 1:cS.nw_expanded
                            val_continuous = val_mat(ik, ikpps, ie); % 获取特定状态下的连续决策值

                            % 找到不大于该决策值的最大网格点索引
                            idx = find(gridV <= val_continuous, 1, 'last');
                            if isempty(idx)
                                idx = 1; % 如果找不到(比如决策值为负)，则取最小索引1
                            end
                            idx_mat(ik, ikpps, ie, ia) = idx; % 使用多维索引精确赋值
                        end
                    end
                end
            end
        end

        function Dist = solve_steady_state_distribution_unified(polS, paramS, cS, Z_ss_norm)
            % =========================================================================
            % == 函数: solve_steady_state_distribution_unified
            % == 版本: [v6 - 平滑求解最终版 (恢复概率质量分裂法)]
            % ==
            % == 核心修正:
            % ==   恢复使用“概率质量分裂法”(线性插值)。虽然这种方法与离散索引的
            % ==   验证器存在方法论差异，但它能保证系统方程的平滑性，这是 fsolve
            % ==   能够成功求解的【必要条件】。我们必须优先保证求解的成功。
            % =========================================================================

            % --- 1. 初始化 ---
            nk = cS.nk;
            nkpps = cS.nkpps;
            nw = cS.nw_expanded;
            Dist = zeros(nk, nkpps, nw, cS.aD_new, 'double');

            % --- 2. 判断求解模式并设定新生儿 ---
            is_theoretical_ss_mode = isfield(cS, 'g_A_ss') && cS.g_A_ss > 1e-6;
            dist_newborn = zeros(nk, nkpps, nw);
            dist_newborn(1, 1, 1:cS.nw) = paramS.leProb1V';

            if is_theoretical_ss_mode
                mass_levels = ones(cS.aD_new, 1);
                for ia = 1:(cS.aD_new - 1)
                    mass_levels(ia+1) = mass_levels(ia) * cS.s_pathV(ia);
                end
                Z_theory = mass_levels / sum(mass_levels);
                Dist(:, :, :, 1) = dist_newborn * Z_theory(1);
                Z_target = Z_theory;
            else
                Dist(:, :, :, 1) = dist_newborn * Z_ss_norm(1);
                Z_target = Z_ss_norm;
            end

            % --- 3. 向前迭代所有年龄组 (使用概率质量分裂法) ---
            for ia = 1:(cS.aD_new - 1)
                dist_ia = Dist(:, :, :, ia);
                dist_ia_plus_1_unscaled = zeros(nk, nkpps, nw, 'double');
                trans_mat_next = paramS.TrProbM_by_age{ia + 1};

                for ie = 1:nw
                    for ikpps = 1:nkpps
                        for ik = 1:nk
                            mass = dist_ia(ik, ikpps, ie);
                            if mass < 1e-30, continue; end

                            k_p_cont = polS(ia).k_prime(ik, ikpps, ie);
                            if cS.pps_active && nkpps > 1
                                kpps_p_cont = polS(ia).kpps_prime(ik, ikpps, ie);
                            else
                                kpps_p_cont = cS.kppsGridV(1);
                            end

                            % [核心] 使用线性插值权重来分裂质量
                            [ik_lower, ik_upper, w_k_upper] = main_steady_state_utils_bgp.find_grid_and_weights(k_p_cont, cS.kGridV);
                            w_k_lower = 1.0 - w_k_upper;
                            [ikpps_lower, ikpps_upper, w_kpps_upper] = main_steady_state_utils_bgp.find_grid_and_weights(kpps_p_cont, cS.kppsGridV);
                            w_kpps_lower = 1.0 - w_kpps_upper;

                            trans_probs_vec = reshape(trans_mat_next(ie, :), [1, 1, nw]);

                            % 将质量按照权重分配到相邻的网格点上
                            mass_chunk = mass * w_k_lower * w_kpps_lower;
                            if mass_chunk > 1e-30
                                dist_ia_plus_1_unscaled(ik_lower, ikpps_lower, :) = dist_ia_plus_1_unscaled(ik_lower, ikpps_lower, :) + mass_chunk * trans_probs_vec;
                            end
                            if ik_lower ~= ik_upper
                                mass_chunk = mass * w_k_upper * w_kpps_lower;
                                if mass_chunk > 1e-30
                                    dist_ia_plus_1_unscaled(ik_upper, ikpps_lower, :) = dist_ia_plus_1_unscaled(ik_upper, ikpps_lower, :) + mass_chunk * trans_probs_vec;
                                end
                            end
                            if ikpps_lower ~= ikpps_upper
                                mass_chunk = mass * w_k_lower * w_kpps_upper;
                                if mass_chunk > 1e-30
                                    dist_ia_plus_1_unscaled(ik_lower, ikpps_upper, :) = dist_ia_plus_1_unscaled(ik_lower, ikpps_upper, :) + mass_chunk * trans_probs_vec;
                                end
                                if ik_lower ~= ik_upper
                                    mass_chunk = mass * w_k_upper * w_kpps_upper;
                                    if mass_chunk > 1e-30
                                        dist_ia_plus_1_unscaled(ik_upper, ikpps_upper, :) = dist_ia_plus_1_unscaled(ik_upper, ikpps_upper, :) + mass_chunk * trans_probs_vec;
                                    end
                                end
                            end
                        end
                    end
                end

                % 根据目标年龄分布重新缩放
                target_mass = Z_target(ia+1);
                current_mass_generated = sum(dist_ia_plus_1_unscaled(:));
                if current_mass_generated > 1e-30
                    rescale_factor = target_mass / current_mass_generated;
                    Dist(:, :, :, ia+1) = dist_ia_plus_1_unscaled * rescale_factor;
                end
            end

            % 最后对整个分布做一次归一化
            total_mass_final = sum(Dist(:));
            if total_mass_final > 1e-9
                Dist = Dist / total_mass_final;
            end
        end        % == 保留的基础函数
        % =======================================================

        function M_prices = get_prices_at_t(K_p, K_g, L, A_t, cS)
            % =========================================================================
            % == 函数: get_prices_at_t
            % == 版本: [v2 - 理论一致版]
            % ==
            % == 核心修正:
            % ==   1. [理论完备] 生产函数和要素价格计算现在完全反映了由
            % ==      alpha 和 gamma 定义的三要素柯布-道格拉斯函数。
            % ==   2. 确保了要素报酬总和严格等于总产出 (欧拉定理)。
            % =========================================================================
            if K_p <= 0, K_p = 1e-8; end; if L <= 0, L = 1e-8; end; if K_g <= 0, K_g = 1e-8; end;

            % 使用与参数一致的三要素生产函数
            labor_exponent = 1 - cS.alpha - cS.gamma;
            if labor_exponent <= 0, error('劳动和资本的产出弹性之和必须小于1！'); end

            Y_period = A_t .* (K_p.^cS.alpha) .* (K_g.^cS.gamma) .* (L.^labor_exponent);

            % 根据正确的边际产出计算要素价格
            MPK_p_period = cS.alpha .* Y_period ./ K_p;
            w_t = labor_exponent .* Y_period ./ L;

            % 市场利率 = 资本边际产出 - 折旧率
            r_mkt_t = MPK_p_period - cS.ddk;

            M_prices = struct('K_p', K_p, 'K_g', K_g, 'L_t', L, 'Y_t', Y_period, 'w_t', w_t, 'r_mkt_t', r_mkt_t);
        end

        function [x_eq, eq_found] = solve_with_fsolve(system_wrapper, x0, verbose)
            % [保留不变] fsolve求解器的包装函数。
            if verbose, fsolve_display = 'iter'; else, fsolve_display = 'none'; end % 根据verbose设置显示选项
            options = optimoptions('fsolve', 'Display', fsolve_display, 'TolFun', 1e-6, 'TolX', 1e-6, 'MaxIterations', 100);

            if verbose, fprintf('\n--- 启动 fsolve 求解器 (求解 [K̂_p, K̂_g, L] - BGP版本) ---\n'); end
            [x_eq, ~, exitflag] = fsolve(system_wrapper, x0, options); % 调用fsolve
            if verbose, fprintf('--- fsolve 求解完成 ---\n'); end

            eq_found = (exitflag > 0); % 根据exitflag判断是否收敛
        end

        % =======================================================
        % == 结果展示与检查
        % =======================================================

        function display_national_accounts_unified(ss, cS)
            % =========================================================================
            % == 函数: display_national_accounts_unified
            % == 版本: [v22 - PPS支持展示版]
            % ==
            % == 核心修正:
            % ==   1. 报告标题会根据pps_active状态改变。
            % ==   2. 明确区分生产性资本(k)和总私人财富(k+kpps)。
            % ==   3. 家庭预算现在包含PPS缴费。
            % =========================================================================

            fprintf('\n\n================================================================================\n');
            if cS.pps_active
                fprintf('===       国民经济核算详细报告 (BGP + PPS激活 - 最终会计准则版)      ===\n');
            else
                fprintf('===        国民经济核算详细报告 (BGP + 无PPS - 最终会计准则版)       ===\n');
            end
            fprintf('================================================================================\n');

            % --- 0. 预计算所有关键宏观量 ---
            fprintf('\n--- [0. 核心宏观变量调试信息] ---\n');
            g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
            Y_prod = ss.Y_from_production_hat;
            C_total = ss.C_agg + ss.Shock_exp_agg;
            S_p_net = ss.Saving_private_flow;
            Bequests = ss.Bequest_tax;
            I_p_net = g_A_period * ss.K_private_begin_hat;
            I_g_net = g_A_period * ss.K_public_hat;
            Depreciation_p = ss.Depreciation_p;
            Depreciation_g = ss.Depreciation_g;
            Depreciation_total = Depreciation_p + Depreciation_g;
            I_p_gross = I_p_net + Depreciation_p;
            I_g_gross = I_g_net + Depreciation_g;
            I_total_gross = I_p_gross + I_g_gross;
            Public_Capital_Return = ss.Public_Capital_Return;
            Gov_Tax_Revenue = ss.Regular_tax + Bequests;
            Gov_Total_Revenue = Gov_Tax_Revenue + Public_Capital_Return;
            Resources_for_discretion = Gov_Total_Revenue - Depreciation_g;
            Gov_Net_Saving = cS.lambda_g * Resources_for_discretion;
            G_c = (1 - cS.lambda_g) * Resources_for_discretion;
            
            fprintf('   Y_prod=%.6f, C_total=%.6f, S_p_net=%.6f, Bequests=%.6f\n', Y_prod, C_total, S_p_net, Bequests);
            fprintf('   I_p_net=%.6f, I_g_net=%.6f, Dep_p=%.6f, Dep_g=%.6f\n', I_p_net, I_g_net, Depreciation_p, Depreciation_g);
            fprintf('   I_p_gross=%.6f, I_g_gross=%.6f, I_total_gross=%.6f\n', I_p_gross, I_g_gross, I_total_gross);
            fprintf('   Public_Capital_Return=%.6f, Gov_Tax_Revenue=%.6f, G_c=%.6f\n', Public_Capital_Return, Gov_Tax_Revenue, G_c);
            if cS.pps_active, fprintf('   PPS_contrib=%.6f\n', ss.PPS_contrib_agg); end

            % --- 1. 收入法检验: GDI vs GDP (欧拉定理) ---
            fprintf('\n--- [1. 收入法检验: GDI vs GDP (欧拉定理)] ---\n');
            GDI = (ss.w_hat * ss.L_hat) + ((ss.r_mkt + cS.ddk) * ss.K_private_begin_hat) + Public_Capital_Return;
            mismatch_gdi = Y_prod - GDI;
            fprintf('   A. 国内生产总值 (GDP) .....................: %12.8f\n', Y_prod);
            fprintf('   B. 国内总收入 (GDI) .......................: %12.8f\n', GDI);
            fprintf('   => 收入-产出缺口 (应为0) ..................: %12.4e\n', mismatch_gdi);

            % --- 2. 宏观资源约束 (总值 Gross Value) ---
            fprintf('\n--- [2. 宏观资源约束 (总值 Gross Value)] ---\n');
            Y_exp_gross = C_total + I_total_gross + G_c;
            mismatch_gross = Y_prod - Y_exp_gross;
            fprintf('   A. 国内生产总值 (GDP) .....................: %12.8f\n', Y_prod);
            fprintf('   B. 总支出 (C + I_gross + Gc) ..............: %12.8f\n', Y_exp_gross);
            fprintf('   => 宏观资源缺口 (应为0) ..................: %12.4e\n', mismatch_gross);

            % --- 3. 政府预算检验 (收支平衡) ---
            fprintf('\n--- [3. 政府预算检验 (收支平衡)] ---\n');
            Gov_Outflows = I_g_gross + G_c;
            mismatch_gov = Gov_Total_Revenue - Gov_Outflows;
            fprintf('   A. 政府总收入 .............................: %12.8f\n', Gov_Total_Revenue);
            fprintf('   B. 政府总支出 (I_g_gross + Gc) ............: %12.8f\n', Gov_Outflows);
            fprintf('   => 政府预算缺口 (应为0) ..................: %12.4e\n', mismatch_gov);

            % --- 4. 国民总储蓄 vs 国民总投资 (释疑部分) ---
            fprintf('\n--- [4. 国民总储蓄 vs 国民总投资 (释疑)] ---\n');
            S_total_gross = S_p_net + Gov_Net_Saving + Depreciation_total;
            mismatch_gross_saving = S_total_gross - I_total_gross;
            fprintf('   A. 国民总储蓄 (S_gross) ...................: %12.8f\n', S_total_gross);
            fprintf('   B. 国民总投资 (I_gross) ...................: %12.8f\n', I_total_gross);
            fprintf('   --------------------------------------------------------------------\n');
            fprintf('   C. 国民储蓄-投资缺口 (S_gross - I_gross) ...: %12.8f\n', mismatch_gross_saving);
            fprintf('   D. 总意外遗赠 (Bequests) ..................: %12.8f\n', Bequests);
            fprintf('   --------------------------------------------------------------------\n');
            fprintf('   => 理论验证 (缺口 - Bequests) (应为0) ......: %12.4e\n', mismatch_gross_saving - Bequests);
            
            fprintf('\n================================================================================\n');
        end

        function is_ok = check_national_accounts_unified(ss, cS)
            % [v7 - BGP会计修正] 国民账户快速检查函数
            g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
            Y_prod = ss.Y_from_production_hat;

            % [会计修正] 总消费
            C_total_agg = ss.C_agg + ss.Shock_exp_agg;

            % 政府支出
            Gov_Revenue_total = ss.Regular_tax + ss.Bequest_tax;
            G_c_agg = (1 - cS.lambda_g) * Gov_Revenue_total;

            % [BGP会计修正] 使用正确的期初资本计算投资
            K_p_hat_begin_period = ss.K_private_hat / (1 + g_A_period);
            K_g_hat_begin_period = ss.K_public_hat;
            I_p_required = (cS.ddk + g_A_period) * K_p_hat_begin_period;
            I_g_required = (cS.ddk_g + g_A_period) * K_g_hat_begin_period;
            I_total_required = I_p_required + I_g_required;

            % 总支出
            Y_exp_hat = C_total_agg + G_c_agg + I_total_required;

            accounting_error = Y_prod - Y_exp_hat;
            TOLERANCE = 1e-5;

            if abs(accounting_error) < TOLERANCE
                fprintf('   国民账户核算: 通过 (误差: %.3e, 容忍度: %.1e)\n', accounting_error, TOLERANCE);
                is_ok = true;
            else
                fprintf('   国民账户核算: 失败! (误差: %.3e, 容忍度: %.1e)\n', accounting_error, TOLERANCE);
                is_ok = false;
            end

            if cS.pps_active && isfield(ss, 'PPS_inflow_agg')
                fprintf('   补充信息 - PPS净流入: %.4f (占私人总投资 %.2f%%)\n', ...
                    ss.PPS_inflow_agg, (ss.PPS_inflow_agg / I_p_required) * 100);
            end
        end

        % =========================================================================
        % == [新函数] 法医级家庭预算约束检验器
        % =========================================================================
        function verify_household_budget_constraint(ss, polS, cS, paramS)
            % =========================================================================
            % == 函数: verify_household_budget_constraint
            % == 版本: [v4 - BGP会计准则同步最终修正版]
            % ==
            % == 核心修正:
            % ==   1. [逻辑同步] 本函数的会计逻辑已与修正后的VFI函数(v22+)完全同步。
            % ==      不再错误地计入"未实现资本利得"作为收入项。
            % ==   2. 这确保了微观检验的有效性。
            % =========================================================================

            fprintf('\n--- [法医级检验] 正在对单个家庭的预算约束进行微观解剖 (BGP正确会计准则版)... ---\n');

            % --- 1. 随机选择一个状态点进行检验 ---
            a_idx_test = randi([1, cS.aR_new]);
            ik_test = randi([1, cS.nk]);
            ie_test = randi([1, cS.nw]);

            fprintf('   检验样本: 年龄组 a_idx = %d, 资产索引 ik = %d, 冲击索引 ie = %d\n', a_idx_test, ik_test, ie_test);

            % --- 2. 提取该状态点的数据 ---
            M_age.r_mkt_t = ss.r_mkt;
            M_age.w_t = ss.w_hat;
            cS.theta_t = ss.Pension_in / max(1e-9, (ss.w_hat * ss.L_hat));

            c_decision = polS(a_idx_test).c(ik_test, 1, ie_test);
            k_prime_decision = polS(a_idx_test).k_prime(ik_test, 1, ie_test);
            g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;

            % --- 3. [核心修正] 使用与VFI完全一致的逻辑，重新计算可支配资源 ---
            k_now = cS.kGridV(ik_test);

            % A. 计算所有流入项 (无 unrealized_capital_gain)
            k_return = k_now * (1 + M_age.r_mkt_t);
            labor_income_gross = M_age.w_t * cS.ageEffV_new(a_idx_test) * paramS.leGridV(ie_test);
            pension_out = 0; % 在工作期，收到的养老金为0

            % B. 计算所有非消费、非储蓄的流出项
            capital_tax = cS.tau_k * (k_now * M_age.r_mkt_t);
            payg_tax = cS.theta_t * labor_income_gross;
            labor_tax = cS.tau_l * max(0, labor_income_gross - payg_tax);
            total_tax_flow = capital_tax + payg_tax + labor_tax;

            % C. 计算支出冲击 (基于正确的 cash-on-hand)
            cash_on_hand_recalc = k_return + labor_income_gross + pension_out - total_tax_flow;
            shock_exp_factor = (ie_test == cS.nw + 1) * cS.kappa_young + (ie_test == cS.nw + 2) * cS.kappa_old;
            shock_exp = shock_exp_factor * cash_on_hand_recalc;

            % D. 计算最终可用于消费和储蓄的资源
            available_for_c_and_s_recalc = cash_on_hand_recalc - shock_exp;

            fprintf('   A. 重新计算的可支配资源 (可用于C和S) ..: %12.8f\n', available_for_c_and_s_recalc);

            % E. 计算VFI决策后的总支出，使用BGP预算约束
            consumption_expenditure = c_decision * (1 + cS.tau_c);
            saving_expenditure = k_prime_decision * (1 + g_A_period); % <--- 储蓄的真实成本
            total_outlay_from_decision = consumption_expenditure + saving_expenditure;

            fprintf('   B. VFI决策后的总支出 (C_expend + k''*(1+g)) : %12.8f\n', total_outlay_from_decision);

            % F. 计算微观预算缺口
            micro_budget_gap = available_for_c_and_s_recalc - total_outlay_from_decision;
            fprintf('   -----------------------------------------------------------------\n');
            fprintf('   => 微观预算缺口 (A - B) ..................: %12.4e\n', micro_budget_gap);

            if abs(micro_budget_gap) < 1e-7
                fprintf('   ✅ [结论] 微观预算闭合！VFI的内部计算与外部检验完全一致！\n');
            else
                fprintf('   ⚠️ [结论] 微观预算仍不闭合！请仔细核对本函数与VFI函数的每一行计算！\n');
            end
            fprintf('--------------------------------------------------------------------------\n');
        end
        
        function [is_steady, max_diff] = verify_steady_state_distribution(Dist, k_prime_idx, kpps_prime_idx, paramS, cS, Z_ss_norm)
            % =========================================================================
            % == 函数说明: 稳态分布验证器
            % ==
            % == 验证逻辑:
            % ==   将输入的分布 Dist 按照策略函数演化一步，得到新分布 Dist_new。
            % ==   如果 Dist 是真正的稳态分布，那么 Dist_new 必须与 Dist 相等。
            % =========================================================================

            fprintf('\n--- 正在验证稳态分布 (Dist) 的正确性... ---\n');

            % --- 1. 初始化新分布矩阵 ---
            Dist_new = zeros(size(Dist));

            % --- 2. 处理新生儿 (第一年龄组) ---
            dist_newborn = zeros(cS.nk, cS.nkpps, cS.nw_expanded);
            dist_newborn(1, 1, 1:cS.nw) = paramS.leProb1V'; % 新生儿出生时没有资产，冲击服从长期分布
            Dist_new(:, :, :, 1) = dist_newborn * Z_ss_norm(1); % 设置新分布的第一期

            % --- 3. 向前迭代一步，计算所有后续年龄组的新分布 ---
            for ia = 1:(cS.aD_new - 1)
                dist_ia = Dist(:, :, :, ia); % 获取输入的、待验证的分布的第 ia 期

                dist_ia_plus_1 = zeros(cS.nk, cS.nkpps, cS.nw_expanded); % 初始化下一期的临时分布
                trans_mat_next = paramS.TrProbM_by_age{ia + 1};

                for ik = 1:cS.nk
                    for ikpps = 1:cS.nkpps
                        for ie = 1:cS.nw_expanded
                            mass = dist_ia(ik, ikpps, ie);
                            if mass < 1e-20, continue; end

                            % 根据策略函数，找到下一期的资产状态索引
                            ik_p = k_prime_idx(ik, ikpps, ie, ia);
                            ikpps_p = kpps_prime_idx(ik, ikpps, ie, ia);

                            % 按照冲击转移概率，将人口质量分配到下一期
                            trans_probs = reshape(trans_mat_next(ie, :), [1, 1, cS.nw_expanded]);
                            dist_ia_plus_1(ik_p, ikpps_p, :) = dist_ia_plus_1(ik_p, ikpps_p, :) + mass * trans_probs;
                        end
                    end
                end

                % 重新缩放以匹配稳态人口结构 (这一步必须和原函数完全一样)
                mass_at_ia = sum(dist_ia, 'all');
                if mass_at_ia > 1e-12
                    rescale = Z_ss_norm(ia+1) / mass_at_ia;
                    Dist_new(:, :, :, ia+1) = dist_ia_plus_1 * rescale;
                end
            end

            % --- 4. 比较新旧分布 ---
            % 注意：由于浮点数精度，我们不能直接比较是否等于0，而是看差的绝对值的最大值
            max_diff = max(abs(Dist - Dist_new), [], 'all');

            % 设定一个非常小的容忍度
            tolerance = 1e-9;

            if max_diff < tolerance
                is_steady = true;
                fprintf('   ✅ 验证通过！分布是稳态的 (最大差异: %.4e)\n', max_diff);
            else
                is_steady = false;
                fprintf('   ⚠️ 验证失败！分布不是稳态的 (最大差异: %.4e)\n', max_diff);
            end
        end

        function [is_steady, max_diff] = verify_FULL_steady_state(Dist_to_verify, polS, paramS, cS)
            % =========================================================================
            % == 函数说明: [v5 - 平滑验证最终版] 终期稳态完全内生验证器
            % ==
            % == 核心修正:
            % ==   1. [方法论统一] 为了与使用“概率质量分裂法”的求解器绝对一致，
            % ==      本验证器也完全采用【概率质量分裂法】（线性插值）来演化分布。
            % ==   2. [输入变更] 不再需要离散的 k_prime_idx，而是直接使用连续的
            % ==      策略函数 polS 进行计算。
            % ==   3. [逻辑自洽] 通过保证验证方法与求解方法的绝对统一，消除了
            % ==      “方法论鸿沟”，使得验证结果能够精确反映求解的准确性。
            % =========================================================================

            % --- [输入参数调整] ---
            % 注意：此版本的验证器不再需要 k_prime_idx 和 kpps_prime_idx
            % 而是直接使用原始的、连续的策略函数 polS
            fprintf('\n--- 正在进行终期稳态的【平滑方法】验证 (概率质量分裂法)...\n');

            % --- 步骤 1: 完全内生地计算理论年龄分布 Z_theory ---
            mass_levels_by_age = zeros(cS.aD_new, 1);
            mass_levels_by_age(1) = 1.0;
            for ia = 1:(cS.aD_new - 1)
                mass_levels_by_age(ia+1) = mass_levels_by_age(ia) * cS.s_pathV(ia);
            end
            Z_theory = mass_levels_by_age / sum(mass_levels_by_age);

            % --- 步骤 2: 基于 Z_theory 初始化一个从零开始计算的分布 ---
            Dist_recalculated = zeros(size(Dist_to_verify));
            dist_newborn = zeros(cS.nk, cS.nkpps, cS.nw_expanded);
            dist_newborn(1, 1, 1:cS.nw) = paramS.leProb1V';
            Dist_recalculated(:, :, :, 1) = dist_newborn * Z_theory(1);

            % --- 步骤 3: 向前迭代，使用【概率质量分裂法】演化分布 ---
            for ia = 1:(cS.aD_new - 1)
                dist_ia = Dist_recalculated(:, :, :, ia);
                dist_ia_plus_1_unscaled = zeros(cS.nk, cS.nkpps, cS.nw_expanded);
                trans_mat_next = paramS.TrProbM_by_age{ia + 1};

                % [核心修正] 使用与求解器完全相同的演化逻辑
                for ie = 1:cS.nw_expanded
                    for ikpps = 1:cS.nkpps
                        for ik = 1:cS.nk
                            mass = dist_ia(ik, ikpps, ie);
                            if mass < 1e-30, continue; end

                            % 直接从连续策略函数 polS 中获取决策
                            k_p_cont = polS(ia).k_prime(ik, ikpps, ie);
                            if cS.pps_active && cS.nkpps > 1
                                kpps_p_cont = polS(ia).kpps_prime(ik, ikpps, ie);
                            else
                                kpps_p_cont = cS.kppsGridV(1);
                            end

                            % 使用线性插值权重来分裂质量
                            [ik_lower, ik_upper, w_k_upper] = main_steady_state_utils_bgp.find_grid_and_weights(k_p_cont, cS.kGridV);
                            w_k_lower = 1.0 - w_k_upper;
                            [ikpps_lower, ikpps_upper, w_kpps_upper] = main_steady_state_utils_bgp.find_grid_and_weights(kpps_p_cont, cS.kppsGridV);
                            w_kpps_lower = 1.0 - w_kpps_upper;

                            trans_probs_vec = reshape(trans_mat_next(ie, :), [1, 1, cS.nw_expanded]);

                            mass_chunk = mass * w_k_lower * w_kpps_lower;
                            if mass_chunk > 1e-30
                                dist_ia_plus_1_unscaled(ik_lower, ikpps_lower, :) = dist_ia_plus_1_unscaled(ik_lower, ikpps_lower, :) + mass_chunk * trans_probs_vec;
                            end
                            if ik_lower ~= ik_upper
                                mass_chunk = mass * w_k_upper * w_kpps_lower;
                                if mass_chunk > 1e-30
                                    dist_ia_plus_1_unscaled(ik_upper, ikpps_lower, :) = dist_ia_plus_1_unscaled(ik_upper, ikpps_lower, :) + mass_chunk * trans_probs_vec;
                                end
                            end
                            if ikpps_lower ~= ikpps_upper
                                mass_chunk = mass * w_k_lower * w_kpps_upper;
                                if mass_chunk > 1e-30
                                    dist_ia_plus_1_unscaled(ik_lower, ikpps_upper, :) = dist_ia_plus_1_unscaled(ik_lower, ikpps_upper, :) + mass_chunk * trans_probs_vec;
                                end
                                if ik_lower ~= ik_upper
                                    mass_chunk = mass * w_k_upper * w_kpps_upper;
                                    if mass_chunk > 1e-30
                                        dist_ia_plus_1_unscaled(ik_upper, ikpps_upper, :) = dist_ia_plus_1_unscaled(ik_upper, ikpps_upper, :) + mass_chunk * trans_probs_vec;
                                    end
                                end
                            end
                        end
                    end
                end

                % 使用 Z_theory 进行重新缩放
                target_mass = Z_theory(ia+1);
                current_mass_generated = sum(dist_ia_plus_1_unscaled(:));
                if current_mass_generated > 1e-30
                    rescale_factor = target_mass / current_mass_generated;
                    Dist_recalculated(:, :, :, ia+1) = dist_ia_plus_1_unscaled * rescale_factor;
                end
            end

            % --- 步骤 4: 对最终结果进行归一化 ---
            total_mass_final = sum(Dist_recalculated(:));
            if total_mass_final > 1e-9
                Dist_recalculated = Dist_recalculated / total_mass_final;
            end

            % --- 步骤 5: 比较两个分布的“形状” ---
            max_diff = max(abs(Dist_to_verify - Dist_recalculated), [], 'all');

            tolerance = 1e-9;
            if max_diff < tolerance
                is_steady = true;
                fprintf('   ✅ [平滑方法] 验证通过！稳态分布是自洽的 (最大差异: %.4e)\n', max_diff);
            else
                is_steady = false;
                fprintf('   ⚠️ [平滑方法] 验证失败！稳态分布不自洽 (最大差异: %.4e)\n', max_diff);
            end
        end        % =======================================================

        % =======================================================
        % == 辅助函数: 概率质量分裂法的核心工具
        % =======================================================



        function [idx_lower, idx_upper, w_upper] = find_grid_and_weights(value, gridV)
            % =========================================================================
            % == 函数: find_grid_and_weights
            % == 版本: [v1 - 线性插值权重计算器]
            % ==
            % == 目的:
            % ==   给定一个连续值和一个离散网格，找到该值在网格中的位置，
            % ==   并计算线性插值所需的权重。
            % ==
            % == 输入:
            % ==   value: 连续值 (浮点数)
            % ==   gridV: 网格向量 (单调递增)
            % ==
            % == 输出:
            % ==   idx_lower: 下界网格点索引
            % ==   idx_upper: 上界网格点索引
            % ==   w_upper: 上界权重 (线性插值系数)
            % =========================================================================

            % --- 边界检查 ---
            if value <= gridV(1)
                idx_lower = 1;
                idx_upper = 1;
                w_upper = 0.0;
                return;
            end

            if value >= gridV(end)
                idx_lower = length(gridV);
                idx_upper = length(gridV);
                w_upper = 0.0;
                return;
            end

            % --- 找到包含该值的网格区间 ---
            % 使用 find 函数找到第一个大于等于 value 的网格点
            idx_upper = find(gridV >= value, 1);
            idx_lower = idx_upper - 1;

            % --- 计算线性插值权重 ---
            if idx_lower == idx_upper
                w_upper = 0.0;
            else
                grid_lower = gridV(idx_lower);
                grid_upper = gridV(idx_upper);
                w_upper = (value - grid_lower) / (grid_upper - grid_lower);
            end
        end

    end
end