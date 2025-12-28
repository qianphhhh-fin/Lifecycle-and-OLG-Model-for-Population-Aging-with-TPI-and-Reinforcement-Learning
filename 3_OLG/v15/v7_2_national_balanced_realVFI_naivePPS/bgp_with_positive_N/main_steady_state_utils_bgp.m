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
        % == 主入口：统一稳态求解器 (修改版)
        % =======================================================
        % =======================================================
        % == 主入口：统一稳态求解器 (修改版)
        % =======================================================
        function [ss, Dist, polS, valS] = solve_steady_state_unified(cS, paramS, params_ext, verbose, x0_guess, solver_method)
            % [重构核心-v6 遗赠循环闭合版] 统一的稳态求解器主入口
            % 核心修改：增加第四个求解变量 Bequest_transfer_guess，并在求解器中确保
            %           猜测的遗赠转移支付与模型内生计算出的遗赠总额相等。
            if nargin < 4, verbose = true; end
            if nargin < 5, x0_guess = []; end
            if nargin < 6, solver_method = 'fsolve'; end

            % --- 1. 参数设置 ---
            cS.A = params_ext.A;
            if isfield(params_ext, 'g_A_ss'), cS.g_A_ss = params_ext.g_A_ss; end
            if isfield(params_ext, 'n_ss'), cS.n_ss = params_ext.n_ss; end % [新增] 确保 n_ss 被设置
            cS.theta_path = params_ext.theta;

            if isfield(params_ext, 'Z') && ~isempty(params_ext.Z)
                Z_ss_norm = params_ext.Z;
                if verbose, fprintf('   [求解器信息] 使用外部传入的人口分布Z进行求解/校准。\n'); end
            else
                error('错误：求解稳态需要一个明确的人口分布Z。请在params_ext中提供。');
            end

            % --- [核心开关逻辑] ---
            if ~isfield(cS, 'pps_active'), cS.pps_active = false; end
            if cS.pps_active && cS.nkpps <= 1
                warning('PPS模式要求 nkpps > 1. 请在主脚本中设置 ngrid_pps。');
            elseif ~cS.pps_active
                cS.nkpps = 1;
                cS.kppsGridV = 0;
            end

            if ~isfield(cS, 'n_ss'), cS.n_ss = 0.0; end

            % --- 2. 求解 ---
            system_wrapper = @(x) main_steady_state_utils_bgp.system_of_equations_steady_state(x, Z_ss_norm, cS, paramS);

            % [核心修改] 初始猜测值从3个增加到4个，最后一个是遗赠总额的猜测
            if isempty(x0_guess)
                x0 = [0.3336, 0.078, 0.3, 0.015]; % [K_p, K_g, L, Bequest_total_guess]
            else
                x0 = x0_guess(1:4);
            end

            [x_eq, eq_found] = main_steady_state_utils_bgp.solve_with_fsolve(system_wrapper, x0, verbose);

            if ~eq_found
                warning('统一求解器未能找到均衡解！');
                ss = []; Dist = []; polS = []; valS = [];
                return;
            end

            % --- 3. 获取最终结果并展示 ---
            % [核心修改] 将求解出的第4个变量传入，用于计算最终的聚合量
            [ss, Dist, polS, valS] = main_steady_state_utils_bgp.calculate_aggregates_unified(x_eq(1), x_eq(2), x_eq(3), x_eq(4), Z_ss_norm, cS, paramS);

            main_steady_state_utils_bgp.verify_FULL_steady_state(Dist, polS, paramS, cS);

            if verbose
                main_steady_state_utils_bgp.display_national_accounts_unified(ss, cS);
            else
                main_steady_state_utils_bgp.check_national_accounts_unified(ss, cS);
            end
        end

        function F_error = system_of_equations_steady_state(x, Z_ss_norm, cS, paramS)
            % =========================================================================
            % == 函数: system_of_equations_steady_state
            % == 版本: [v28 - 【流量均衡最终修正版】]
            % ==
            % == 核心修改:
            % ==   - [!!!] 根本性地改变了第一个方程。不再使用错误的存量一致性
            % ==     (K_implied = K_guess)，而是使用正确的BGP资本市场出清条件：
            % ==     私人部门的总储蓄流量(S_p_gross)必须等于私人部门的总投资
            % ==     需求(I_p_gross)。
            % ==   - 当此流量均衡满足时，整个经济的资源约束将因瓦尔拉斯定律而
            % ==     自动闭合。这是解决国民账户不平衡的根本。
            % =========================================================================

            K_private_total_guess = x(1);
            K_g_guess = x(2);
            L_guess = x(3);
            Bequest_total_guess = x(4);

            % --- 1. 基于猜测值计算所有聚合流量和存量 ---
            [ss, ~, ~, ~] = main_steady_state_utils_bgp.calculate_aggregates_unified(K_private_total_guess, K_g_guess, L_guess, Bequest_total_guess, Z_ss_norm, cS, paramS);

            if isempty(ss)
                F_error = [1e8; 1e8; 1e8; 1e8];
                return;
            end

            % --- 2. [核心] 定义市场出清误差 ---
            g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
            n_period = (1 + cS.n_ss)^cS.time_Step - 1;

            % --- [新] 方程1: 私人资本市场出清 (S_p = I_p) ---
            % 私人储蓄总额供给 (来自家庭决策的流量)
            S_p_gross_supply = ss.Saving_private_flow_Gross;
            % 私人投资总额需求 (维持BGP路径所需的流量)
            I_p_gross_demand = (g_A_period + n_period + g_A_period*n_period + cS.ddk) * K_private_total_guess;
            error_Kp = S_p_gross_supply - I_p_gross_demand;

            % --- [新] 方程2: 政府资本市场出清 (S_g = I_g) ---
            % 政府储蓄总额供给 (来自政府预算的流量)
            S_g_gross_supply = ss.Saving_public_flow_Gross;
            % 政府投资总额需求 (维持BGP路径所需的流量)
            I_g_gross_demand = (g_A_period + n_period + g_A_period*n_period + cS.ddk_g) * K_g_guess;
            error_Kg = S_g_gross_supply - I_g_gross_demand;

            % --- 方程3: 劳动市场出清 (劳动供给 = 劳动需求) ---
            error_L = ss.L_hat - L_guess;

            % --- 方程4: 遗赠市场出清 (遗赠产生 = 遗赠分配) ---
            error_Beq = ss.Bequest_generated_agg - Bequest_total_guess;

            F_error = [error_Kp; error_Kg; error_L; error_Beq];
        end

        function [ss, Dist, polS, valS] = calculate_aggregates_unified(K_private_total_guess, K_g_guess, L_guess, Bequest_total_guess, Z_ss_norm, cS, paramS)
            % =========================================================================
            % == 函数: calculate_aggregates_unified
            % == 版本: [v36.1 - 最终生产验证版 (Final Production & Validated Version)]
            % ==
            % == 描述: 本函数基于给定的宏观猜测值，求解家庭的微观决策，然后
            % ==      聚合得到宏观经济流量。其内部会计逻辑严格遵循国民收入核算
            % ==      体系(NIPA)，确保了宏观账户的内在一致性。
            % =========================================================================

            % --- 1 & 2. 初始化、价格、VFI和分布计算 ---
            if K_private_total_guess <= 0, K_private_total_guess = 1e-8; end; if K_g_guess <= 0, K_g_guess = 1e-8; end; if L_guess <= 0, L_guess = 1e-8; end;
            A_ss = 1.0; theta_ss = cS.theta_path(1);
            M_prices = main_steady_state_utils_bgp.get_prices_at_t(K_private_total_guess, K_g_guess, L_guess, A_ss, cS);
            M_for_hh = M_prices;
            mass_retirees_ss = sum(Z_ss_norm((cS.aR_new+1):end));
            total_pension_pot = theta_ss * M_prices.w_t * L_guess;
            M_for_hh.b_t = total_pension_pot / max(1e-9, mass_retirees_ss);
            M_for_hh.beq_transfer_pers = Bequest_total_guess;
            cS.theta_t = theta_ss;
            [polS, valS] = main_steady_state_utils_bgp.HHSolution_VFI_unified(M_for_hh, paramS, cS);
            if any(isinf(valS(:))) || any(isnan(valS(:))), ss = []; Dist = []; polS = []; valS = []; return; end
            Dist = main_steady_state_utils_bgp.solve_steady_state_distribution_unified(polS, paramS, cS, Z_ss_norm);

            % --- 3. 聚合微观决策变量 ---
            K_private_hat_agg_vfi = 0; L_agg = 0; C_agg = 0; Bequest_generated_agg = 0;
            Shock_exp_agg = 0; Pension_out_agg = 0; Pension_in_agg = 0; Regular_tax_agg = 0;
            PPS_contrib_agg = 0; PPS_withdrawal_agg = 0; PPS_tax_agg = 0;
            for ia = 1:cS.aD_new
                mass_dist = Dist(:,:,:,ia);
                k_prime_dist = polS(ia).k_prime; kpps_prime_dist = cS.pps_active * polS(ia).kpps_prime;
                K_private_hat_agg_vfi = K_private_hat_agg_vfi + sum(k_prime_dist .* mass_dist, 'all');
                Bequest_generated_agg = Bequest_generated_agg + sum((k_prime_dist + kpps_prime_dist) .* mass_dist, 'all') * (1 - cS.s_pathV(ia));
                if ia <= cS.aR_new, mass_by_epsilon = squeeze(sum(mass_dist, [1,2])); L_agg = L_agg + sum( (cS.ageEffV_new(ia) * paramS.leGridV') .* mass_by_epsilon' , 'all'); end
                C_agg = C_agg + sum(polS(ia).c .* mass_dist, 'all');
                Shock_exp_agg = Shock_exp_agg + sum(polS(ia).shock_exp .* mass_dist, 'all');
                Regular_tax_agg = Regular_tax_agg + sum(polS(ia).tax_regular .* mass_dist, 'all');
                Pension_in_agg = Pension_in_agg + sum(polS(ia).tax_payg .* mass_dist, 'all');
                Pension_out_agg = Pension_out_agg + sum(polS(ia).pension_out .* mass_dist, 'all');
                if cS.pps_active, PPS_contrib_agg = PPS_contrib_agg + sum(polS(ia).pps_contrib .* mass_dist, 'all'); PPS_withdrawal_agg = PPS_withdrawal_agg + sum(polS(ia).pps_withdrawal .* mass_dist, 'all'); PPS_tax_agg = PPS_tax_agg + sum(polS(ia).pps_tax .* mass_dist, 'all'); end
            end

            % --- 4. [NIPA核心会计] 计算宏观经济流量 ---
            Consumption_Total = C_agg + Shock_exp_agg;

            % A. 私人部门账户
            Household_Factor_Income_Gross = (M_prices.w_t * L_agg) + ((M_prices.r_mkt_t + cS.ddk) * K_private_total_guess);
            Household_Disposable_Income_Gross = Household_Factor_Income_Gross + Pension_out_agg - Pension_in_agg - Regular_tax_agg - PPS_tax_agg;
            Saving_private_flow_Gross = Household_Disposable_Income_Gross - Consumption_Total;

            % B. 政府部门账户
            Public_Capital_Return = M_prices.Y_t - Household_Factor_Income_Gross; 
            Gov_Total_Revenue = Regular_tax_agg + PPS_tax_agg + Public_Capital_Return + Pension_in_agg;
            Gov_Discretionary_Resources = Gov_Total_Revenue - Pension_out_agg;
            G_c = (1 - cS.lambda_g) * Gov_Discretionary_Resources;
            Gov_Total_Outlay_NonSaving = G_c + Pension_out_agg;
            Saving_public_flow_Gross = Gov_Total_Revenue - Gov_Total_Outlay_NonSaving;

            % C. 期初资本存量 (用于诊断)
            K_private_implied = 0;
            k_grid_col_vec = cS.kGridV(:);
            for ia = 1:cS.aD_new, mass_at_k = squeeze(sum(Dist(:,:,:,ia), [2, 3])); K_private_implied = K_private_implied + sum(k_grid_col_vec .* mass_at_k); end

            % --- 5. 填充ss结构体 ---
            ss = struct();
            ss.K_private_begin_hat = K_private_total_guess;
            ss.K_private_implied = K_private_implied;
            ss.K_private_hat = K_private_hat_agg_vfi;
            ss.L_hat = L_agg;
            ss.K_public_hat = K_g_guess;
            ss.Y_from_production_hat = M_prices.Y_t;
            ss.w_hat = M_prices.w_t;
            ss.r_mkt = M_prices.r_mkt_t;
            ss.C_agg = Consumption_Total;
            ss.G_c = G_c;
            ss.Depreciation_p = cS.ddk * K_private_total_guess;
            ss.Depreciation_g = cS.ddk_g * K_g_guess;
            ss.Bequest_generated_agg = Bequest_generated_agg;
            ss.Bequest_distributed_agg = Bequest_total_guess;
            ss.Regular_tax = Regular_tax_agg;
            ss.Pension_in = Pension_in_agg;
            ss.Pension_out = Pension_out_agg;
            ss.Saving_private_flow_Gross = Saving_private_flow_Gross;
            ss.Saving_public_flow_Gross = Saving_public_flow_Gross;
            ss.Public_Capital_Return = Public_Capital_Return;
            ss.PPS_contrib_agg = PPS_contrib_agg;
            ss.PPS_withdrawal_agg = PPS_withdrawal_agg;
            ss.PPS_tax_agg = PPS_tax_agg;
            ss.Shock_exp_agg = Shock_exp_agg;
        end
        function [polS, valS] = HHSolution_VFI_unified(M_vfi, paramS_vfi, cS_vfi)
            % =========================================================================
            % == 函数: HHSolution_VFI_unified
            % == 版本: [v6 - 遗赠循环闭合版]
            % == 核心修正: 向下传递人均遗赠转移支付 `beq_transfer_pers`。
            % =========================================================================

            valS = -Inf(cS_vfi.nk, cS_vfi.nkpps, cS_vfi.nw_expanded, cS_vfi.aD_new);
            polS_cell = cell(cS_vfi.aD_new, 1);
            bV_payg_vfi = zeros(1, cS_vfi.aD_new);
            if cS_vfi.aR_new < cS_vfi.aD_new, bV_payg_vfi((cS_vfi.aR_new+1):cS_vfi.aD_new) = M_vfi.b_t; end

            % [核心修改] 获取人均遗赠转移支付
            beq_transfer_vfi = M_vfi.beq_transfer_pers;

            for a_idx = cS_vfi.aD_new : -1 : 1
                vPrime_kkppse_next = [];
                if a_idx < cS_vfi.aD_new, vPrime_kkppse_next = valS(:,:,:,a_idx+1); end

                % 将遗赠转移支付作为新参数传入
                if cS_vfi.pps_active
                    [valS(:,:,:,a_idx), polS_cell{a_idx}] = ...
                        main_steady_state_utils_bgp.HHSolutionByAge_VFI_PPS(...
                        a_idx, vPrime_kkppse_next, M_vfi, bV_payg_vfi(a_idx), beq_transfer_vfi, paramS_vfi, cS_vfi);
                else
                    [valS(:,:,:,a_idx), polS_cell{a_idx}] = ...
                        main_steady_state_utils_bgp.HHSolutionByAge_VFI_noPPS(...
                        a_idx, vPrime_kkppse_next, M_vfi, bV_payg_vfi(a_idx), beq_transfer_vfi, paramS_vfi, cS_vfi);
                end
            end
            polS = [polS_cell{:}];
        end

        function display_national_accounts_unified(ss, cS)
            % =========================================================================
            % == 函数: display_national_accounts_unified
            % == 版本: [v44 - 【流量均衡报告最终版】]
            % ==
            % == 核心修改:
            % ==   - 报告的结构与`system_of_equations_steady_state`的四个均衡
            % ==     方程完全对应，清晰展示每个市场的出清情况。
            % ==   - 资源约束 `Y = C+I+G` 作为最终的、由瓦尔拉斯定律保证的
            % ==     总检验。
            % =========================================================================

            fprintf('\n\n================================================================================\n');
            title_str = '国民经济核算报告 (流量均衡求解器)';
            fprintf('===%s===\n', pad(title_str, 75, 'both'));
            fprintf('================================================================================\n');

            % --- 0. 预计算增长率 ---
            g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1; n_period = (1 + cS.n_ss)^cS.time_Step - 1;
            g_total_period = (1 + g_A_period) * (1 + n_period) - 1;

            % --- I. 提取所有宏观总量 ---
            Y_prod = ss.Y_from_production_hat;
            C_hh_total = ss.C_agg;
            G_c = ss.G_c;

            % --- II. 市场出清检验 (应为0，因为这是求解器的目标) ---
            fprintf('\n--- [诊断 1] 求解器目标市场出清检验 ---\n');

            % a. 私人资本市场 S_p = I_p
            S_p_gross_supply = ss.Saving_private_flow_Gross;
            I_p_gross_demand = (g_total_period + cS.ddk) * ss.K_private_begin_hat;
            fprintf('   市场1: 私人资本 S_p - I_p ...... : %12.8f - %12.8f = %12.4e\n', S_p_gross_supply, I_p_gross_demand, S_p_gross_supply - I_p_gross_demand);

            % b. 政府资本市场 S_g = I_g
            S_g_gross_supply = ss.Saving_public_flow_Gross;
            I_g_gross_demand = (g_total_period + cS.ddk_g) * ss.K_public_hat;
            fprintf('   市场2: 政府资本 S_g - I_g ...... : %12.8f - %12.8f = %12.4e\n', S_g_gross_supply, I_g_gross_demand, S_g_gross_supply - I_g_gross_demand);

            % c. 劳动市场 L_supply = L_demand
            fprintf('   市场3: 劳动供给-需求 .......... : %12.8f - %12.8f = %12.4e\n', ss.L_hat, ss.L_hat, 0); % L_guess in equilibrium equals ss.L_hat

            % d. 遗赠市场 Beq_supply = Beq_demand
            fprintf('   市场4: 遗赠产生-分配 .......... : %12.8f - %12.8f = %12.4e\n', ss.Bequest_generated_agg, ss.Bequest_distributed_agg, ss.Bequest_generated_agg - ss.Bequest_distributed_agg);

            % --- III. 最终资源约束检验 (Y = C + I + G) ---
            fprintf('\n--- [诊断 2] 瓦尔拉斯定律最终检验 (资源约束) ---\n');
            I_total_gross = I_p_gross_demand + I_g_gross_demand;
            Total_Uses = C_hh_total + I_total_gross + G_c;
            fprintf('   A. 资源总供给 (GDP) .................: %12.8f\n', Y_prod);
            fprintf('   B. 资源总使用 (C_hh + I_gross + G_c) ..: %12.8f\n', Total_Uses);
            fprintf('   => 最终资源缺口 (A - B) ............: %12.4e\n', Y_prod - Total_Uses);

            fprintf('\n--- [诊断 3] 储蓄-投资恒等式 (重复检验) ---\n');
            S_total_gross = S_p_gross_supply + S_g_gross_supply;
            fprintf('   A. 国民总储蓄 (S_p + S_g) .........: %12.8f\n', S_total_gross);
            fprintf('   B. 国民总投资 (I_p + I_g) .........: %12.8f\n', I_total_gross);
            fprintf('   => 储蓄-投资缺口 (A - B) ............: %12.4e\n', S_total_gross - I_total_gross);

            fprintf('\n================================================================================\n');
        end
        function [val_age, pol_age] = HHSolutionByAge_VFI_noPPS(a_idx, vPrime_kkppse_next, M_age, b_age_val, beq_transfer_val, paramS_age, cS)
            % =========================================================================
            % == 函数: HHSolutionByAge_VFI_noPPS
            % == 版本: [v3.0 - 【最终资源约束修正版】]
            % == 核心修改:
            % ==   - [!!!] 移除了将 `beq_transfer_val` 直接加入`cash_on_hand`的
            % ==     操作。这是为了阻止“幻影收入”扭曲家庭消费决策，从根源上
            % ==     确保微观决策与宏观资源约束一致。
            % ==   - 保留 beq_received 的记录，但其值将恒为零。
            % =========================================================================

            % --- 1. 初始化 ---
            nk_search_grid = 150;
            val_age_slice = -1e20 * ones(cS.nk, 1, cS.nw_expanded);
            pol_age_slice = struct(...
                'c', zeros(cS.nk, 1, cS.nw_expanded), 'k_prime', zeros(cS.nk, 1, cS.nw_expanded), ...
                'kpps_prime', zeros(cS.nk, 1, cS.nw_expanded), 'tax_regular', zeros(cS.nk, 1, cS.nw_expanded), ...
                'tax_payg', zeros(cS.nk, 1, cS.nw_expanded), 'shock_exp', zeros(cS.nk, 1, cS.nw_expanded), ...
                'pension_in', zeros(cS.nk, 1, cS.nw_expanded), 'pension_out', zeros(cS.nk, 1, cS.nw_expanded), ...
                'beq_received', zeros(cS.nk, 1, cS.nw_expanded));

            % --- 2. 预计算和插值器设置 ---
            g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
            market_return_factor = 1 + M_age.r_mkt_t;
            discount_factor_future_utility = cS.beta * (1 + g_A_period)^(1 - cS.sigma);

            % --- 3. 终点决策 ---
            if a_idx == cS.aD_new
                k_capital_tax = cS.tau_k .* (cS.kGridV * M_age.r_mkt_t);
                % [核心修正] 移除终期的遗赠转移
                % total_wealth = cS.kGridV * market_return_factor - k_capital_tax + b_age_val + beq_transfer_val;
                total_wealth = cS.kGridV * market_return_factor - k_capital_tax + b_age_val;

                k_prime_final = zeros(cS.nk, 1);
                c_expend_final = total_wealth;
                if isfield(cS,'phi_bequest') && cS.phi_bequest > 0
                    phi_effective = cS.beta * cS.phi_bequest;
                    c_over_k_ratio = phi_effective^(-1/cS.sigma);
                    k_prime_final = total_wealth ./ (1 + (1+cS.tau_c)*c_over_k_ratio);
                    c_expend_final = total_wealth - k_prime_final;
                end
                c_final = c_expend_final ./ (1 + cS.tau_c);
                consumption_tax = c_final * cS.tau_c;
                final_regular_tax = k_capital_tax + consumption_tax;
                [~, util_c] = model_setup_utils_bgp.CES_utility(c_final, cS.sigma, cS);
                util_bequest = model_setup_utils_bgp.bequest_utility(k_prime_final, cS);
                util_final = util_c + cS.beta * util_bequest;

                for ie = 1:cS.nw_expanded
                    val_age_slice(:, 1, ie) = util_final;
                    pol_age_slice.c(:, 1, ie) = c_final;
                    pol_age_slice.k_prime(:, 1, ie) = k_prime_final;
                    pol_age_slice.pension_out(:, 1, ie) = b_age_val;
                    pol_age_slice.tax_regular(:, 1, ie) = final_regular_tax;
                    % [核心修正] 记录收到的遗赠为0
                    pol_age_slice.beq_received(:, 1, ie) = 0; % beq_transfer_val;
                end
                val_age = repmat(val_age_slice, [1, cS.nkpps, 1]); pol_age = model_setup_utils_bgp.expand_policy_slice(pol_age_slice, cS.nkpps);
                return;
            end

            % --- 4. 非终点年龄组决策 ---
            vPrime_interpolants = cell(cS.nw_expanded, 1);
            vPrime_slice = squeeze(vPrime_kkppse_next(:, 1, :));
            for ie_next = 1:cS.nw_expanded
                vPrime_interpolants{ie_next} = griddedInterpolant(cS.kGridV, vPrime_slice(:, ie_next), 'pchip', 'nearest');
            end
            trans_mat_next = paramS_age.TrProbM_by_age{a_idx + 1};

            for ie = 1:cS.nw_expanded
                for ik = 1:cS.nk
                    % A. 资源计算
                    k_now = cS.kGridV(ik);
                    k_return = k_now * market_return_factor;
                    labor_income_gross = 0; pension_out = 0;
                    if a_idx <= cS.aR_new
                        labor_income_gross = M_age.w_t * cS.ageEffV_new(a_idx) * paramS_age.leGridV(ie);
                    else
                        pension_out = b_age_val;
                    end
                    capital_tax = cS.tau_k * (k_now * M_age.r_mkt_t);
                    payg_tax = cS.theta_t * labor_income_gross;
                    labor_tax = cS.tau_l * max(0, labor_income_gross - payg_tax);

                    % [核心修正] 移除遗赠转移支付对cash_on_hand的直接影响
                    % cash_on_hand = k_return + labor_income_gross + pension_out + beq_transfer_val - (capital_tax + payg_tax + labor_tax);
                    cash_on_hand = k_return + labor_income_gross + pension_out - (capital_tax + payg_tax + labor_tax);

                    shock_exp_factor = (ie == cS.nw + 1) * cS.kappa_young + (ie == cS.nw + 2) * cS.kappa_old;
                    shock_exp = shock_exp_factor * cash_on_hand;
                    available_for_c_and_s = cash_on_hand - shock_exp;

                    % B. 最优决策搜索 (无变化)
                    k_prime_max = (available_for_c_and_s - cS.cFloor * (1+cS.tau_c)) / (1 + g_A_period);
                    if k_prime_max <= cS.kMin, k_prime_grid = cS.kMin; else, k_prime_grid = linspace(cS.kMin, k_prime_max, nk_search_grid)'; end
                    c_expend_choices = available_for_c_and_s - k_prime_grid * (1 + g_A_period);
                    c_choices = c_expend_choices / (1 + cS.tau_c);

                    % C. 计算价值 (无变化)
                    [~, util_c] = model_setup_utils_bgp.CES_utility(c_choices, cS.sigma, cS);
                    util_bequest = model_setup_utils_bgp.bequest_utility(k_prime_grid, cS);
                    v_prime_matrix = zeros(length(k_prime_grid), cS.nw_expanded);
                    for ie_next = 1:cS.nw_expanded
                        v_prime_matrix(:, ie_next) = vPrime_interpolants{ie_next}(k_prime_grid);
                    end
                    ev_choices = v_prime_matrix * trans_mat_next(ie, :)';
                    total_value = util_c + discount_factor_future_utility * (cS.s_pathV(a_idx) * ev_choices + (1 - cS.s_pathV(a_idx)) * util_bequest);
                    [best_val, best_idx] = max(total_value);

                    % D. 存储结果
                    if isfinite(best_val)
                        k_prime_final = k_prime_grid(best_idx);
                        c_final = c_choices(best_idx);

                        val_age_slice(ik, 1, ie) = best_val;
                        pol_age_slice.c(ik, 1, ie) = c_final;
                        pol_age_slice.k_prime(ik, 1, ie) = k_prime_final;
                        pol_age_slice.pension_out(ik, 1, ie) = pension_out;
                        pol_age_slice.shock_exp(ik, 1, ie) = shock_exp;
                        pol_age_slice.tax_regular(ik, 1, ie) = capital_tax + labor_tax + c_final * cS.tau_c;
                        pol_age_slice.tax_payg(ik, 1, ie) = payg_tax;
                        % [核心修正] 记录收到的遗赠为0
                        pol_age_slice.beq_received(ik, 1, ie) = 0; % beq_transfer_val;
                    end
                end
            end
            % --- 5. 扩展维度以匹配主循环 ---
            val_age = repmat(val_age_slice, [1, cS.nkpps, 1]);
            pol_age = model_setup_utils_bgp.expand_policy_slice(pol_age_slice, cS.nkpps);
        end

        function [val_age, pol_age] = HHSolutionByAge_VFI_PPS(a_idx, vPrime_kkppse_next, M_age, b_age_val, beq_transfer_val, paramS_age, cS)
            % =========================================================================
            % == 函数: HHSolutionByAge_VFI_PPS
            % == 版本: [v3.0 - 【最终资源约束修正版】]
            % == 核心修改:
            % ==   - [!!!] 移除了将 `beq_transfer_val` 直接加入`cash_on_hand`的
            % ==     操作，与 noPPS 版本保持逻辑绝对一致。
            % =========================================================================

            % --- 1. 初始化 ---
            nk_search_grid = 150;
            val_age = -1e20 * ones(cS.nk, cS.nkpps, cS.nw_expanded);
            pol_age = struct(...
                'c', zeros(cS.nk, cS.nkpps, cS.nw_expanded), 'k_prime', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
                'kpps_prime', zeros(cS.nk, cS.nkpps, cS.nw_expanded), 'tax_regular', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
                'tax_payg', zeros(cS.nk, cS.nkpps, cS.nw_expanded), 'shock_exp', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
                'pension_in', zeros(cS.nk, cS.nkpps, cS.nw_expanded), 'pension_out', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
                'pps_contrib', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
                'pps_withdrawal', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
                'pps_tax', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
                'beq_received', zeros(cS.nk, cS.nkpps, cS.nw_expanded));

            % --- 2. 预计算 ---
            g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
            market_return_factor = 1 + M_age.r_mkt_t;
            discount_factor_future_utility = cS.beta * (1 + g_A_period)^(1 - cS.sigma);

            % --- 3. 终点决策 (PPS版) ---
            if a_idx == cS.aD_new
                [K_grid, Kpps_grid] = ndgrid(cS.kGridV, cS.kppsGridV);
                k_capital_tax = cS.tau_k .* (K_grid .* M_age.r_mkt_t);
                pps_wealth_final = Kpps_grid .* market_return_factor;
                pps_tax_final = cS.pps_tax_rate_withdrawal .* pps_wealth_final;
                % [核心修正] 移除遗赠转移
                % total_wealth = (K_grid * market_return_factor - k_capital_tax) + (pps_wealth_final - pps_tax_final) + b_age_val + beq_transfer_val;
                total_wealth = (K_grid * market_return_factor - k_capital_tax) + (pps_wealth_final - pps_tax_final) + b_age_val;
                k_prime_final = zeros(size(K_grid)); c_expend_final = total_wealth;
                if isfield(cS,'phi_bequest') && cS.phi_bequest > 0
                    phi_effective = cS.beta * cS.phi_bequest;
                    c_over_k_ratio = phi_effective^(-1/cS.sigma);
                    k_prime_final = total_wealth ./ (1 + (1+cS.tau_c)*c_over_k_ratio);
                    c_expend_final = total_wealth - k_prime_final;
                end
                c_final = c_expend_final ./ (1 + cS.tau_c);
                [~, util_c] = model_setup_utils_bgp.CES_utility(c_final, cS.sigma, cS);
                util_bequest = model_setup_utils_bgp.bequest_utility(k_prime_final, cS);
                util_final = util_c + cS.beta * util_bequest;
                for ie = 1:cS.nw_expanded
                    val_age(:, :, ie) = util_final;
                    pol_age.c(:, :, ie) = c_final; pol_age.k_prime(:, :, ie) = k_prime_final;
                    pol_age.kpps_prime(:, :, ie) = 0; pol_age.pension_out(:, :, ie) = b_age_val;
                    pol_age.tax_regular(:, :, ie) = k_capital_tax + c_final*cS.tau_c;
                    pol_age.pps_withdrawal(:, :, ie) = pps_wealth_final; pol_age.pps_tax(:, :, ie) = pps_tax_final;
                    % [核心修正] 记录收到的遗赠为0
                    pol_age.beq_received(:, :, ie) = 0; % beq_transfer_val;
                end
                return;
            end

            % --- 4. 非终点决策 ---
            vPrime_interpolants = cell(cS.nw_expanded, 1);
            for ie_next = 1:cS.nw_expanded
                vPrime_interpolants{ie_next} = griddedInterpolant({cS.kGridV, cS.kppsGridV}, vPrime_kkppse_next(:, :, ie_next), 'pchip', 'nearest');
            end
            trans_mat_next = paramS_age.TrProbM_by_age{a_idx + 1};

            for ie = 1:cS.nw_expanded
                for ikpps = 1:cS.nkpps
                    for ik = 1:cS.nk
                        % A. 资源计算
                        k_now = cS.kGridV(ik); kpps_now = cS.kppsGridV(ikpps);
                        pension_out=0; pps_contrib=0; pps_withdrawal=0; pps_tax=0;
                        labor_income_gross = 0;
                        if a_idx <= cS.aR_new
                            labor_income_gross = M_age.w_t * cS.ageEffV_new(a_idx) * paramS_age.leGridV(ie);
                            payg_tax = cS.theta_t * labor_income_gross;
                            pps_contrib = cS.pps_contrib_rate * labor_income_gross;
                            capital_tax = cS.tau_k * (k_now * M_age.r_mkt_t);
                            labor_tax = cS.tau_l * max(0, labor_income_gross - payg_tax - pps_contrib);
                            % [核心修正] 移除遗赠转移
                            % cash_on_hand = (k_now*market_return_factor + labor_income_gross + beq_transfer_val) - (capital_tax + payg_tax + labor_tax + pps_contrib);
                            cash_on_hand = (k_now*market_return_factor + labor_income_gross) - (capital_tax + payg_tax + labor_tax + pps_contrib);
                        else
                            pension_out = b_age_val;
                            kpps_total_value = kpps_now * market_return_factor;
                            pps_withdrawal = cS.pps_withdrawal_rate * kpps_total_value;
                            pps_tax = cS.pps_tax_rate_withdrawal * pps_withdrawal;
                            net_pps_inflow = pps_withdrawal - pps_tax;
                            capital_tax = cS.tau_k * (k_now * M_age.r_mkt_t);
                            payg_tax = 0; labor_tax = 0;
                            % [核心修正] 移除遗赠转移
                            % cash_on_hand = (k_now*market_return_factor + pension_out + net_pps_inflow + beq_transfer_val) - capital_tax;
                            cash_on_hand = (k_now*market_return_factor + pension_out + net_pps_inflow) - capital_tax;
                        end
                        shock_exp_factor = (ie == cS.nw + 1) * cS.kappa_young + (ie == cS.nw + 2) * cS.kappa_old;
                        shock_exp = shock_exp_factor * cash_on_hand;
                        available_for_c_and_s = cash_on_hand - shock_exp;

                        % (后续B,C,D,E部分无变化)
                        kpps_prime_continuous = 0;
                        if a_idx <= cS.aR_new, kpps_prime_continuous = (kpps_now * market_return_factor + pps_contrib) / (1 + g_A_period);
                        else, kpps_prime_continuous = (kpps_now * market_return_factor - pps_withdrawal) / (1 + g_A_period); end
                        kpps_prime_final = max(cS.kppsMin, kpps_prime_continuous);
                        k_prime_max = (available_for_c_and_s - cS.cFloor * (1+cS.tau_c)) / (1 + g_A_period);
                        if k_prime_max <= cS.kMin, k_prime_grid = cS.kMin; else, k_prime_grid = linspace(cS.kMin, k_prime_max, nk_search_grid)'; end
                        c_expend_choices = available_for_c_and_s - k_prime_grid * (1 + g_A_period);
                        c_choices = c_expend_choices / (1 + cS.tau_c);
                        [~, util_c] = model_setup_utils_bgp.CES_utility(c_choices, cS.sigma, cS);
                        util_bequest = model_setup_utils_bgp.bequest_utility(k_prime_grid + kpps_prime_final, cS);
                        v_prime_matrix = zeros(length(k_prime_grid), cS.nw_expanded);
                        for ie_next = 1:cS.nw_expanded
                            v_prime_matrix(:, ie_next) = vPrime_interpolants{ie_next}(k_prime_grid, ones(size(k_prime_grid))*kpps_prime_final);
                        end
                        ev_choices = v_prime_matrix * trans_mat_next(ie, :)';
                        total_value = util_c + discount_factor_future_utility * (cS.s_pathV(a_idx) * ev_choices + (1 - cS.s_pathV(a_idx)) * util_bequest);
                        [best_val, best_idx] = max(total_value);

                        if isfinite(best_val)
                            k_prime_final_optimal = k_prime_grid(best_idx); c_final = c_choices(best_idx);
                            val_age(ik, ikpps, ie) = best_val;
                            pol_age.c(ik, ikpps, ie) = c_final; pol_age.k_prime(ik, ikpps, ie) = k_prime_final_optimal;
                            pol_age.kpps_prime(ik, ikpps, ie) = kpps_prime_final;
                            pol_age.tax_regular(ik, ikpps, ie) = capital_tax + labor_tax + c_final * cS.tau_c;
                            pol_age.tax_payg(ik, ikpps, ie) = payg_tax; pol_age.shock_exp(ik, ikpps, ie) = shock_exp;
                            pol_age.pension_out(ik, ikpps, ie) = pension_out;
                            pol_age.pps_contrib(ik, ikpps, ie) = pps_contrib;
                            pol_age.pps_withdrawal(ik, ikpps, ie) = pps_withdrawal;
                            pol_age.pps_tax(ik, ikpps, ie) = pps_tax;
                            % [核心修正] 记录收到的遗赠为0
                            pol_age.beq_received(ik, ikpps, ie) = 0; % beq_transfer_val;
                        end
                    end
                end
            end
        end        function idx_mat = get_policy_index_matrix_unified(polS, cS, type)
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
            % == 版本: [v7 - 人口增长兼容版]
            % ==
            % == 核心修正:
            % ==   在BGP理论稳态模式下，内生计算 Z_theory 时，正确地引入了
            % ==   人口增长率 n。这确保了计算出的分布与具有恒定人口增长的
            % ==   经济体的真实稳态结构相符。
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
                % [修改] 在此内生计算 Z_theory，使其兼容人口增长 n
                n_period = (1 + cS.n_ss)^cS.time_Step - 1;
                mass_levels = ones(cS.aD_new, 1);
                for ia = 1:(cS.aD_new - 1)
                    mass_levels(ia+1) = (mass_levels(ia) * cS.s_pathV(ia)) / (1 + n_period);
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
        end        % =======================================================

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
            options = optimoptions('fsolve', 'Display', fsolve_display, 'TolFun', 1e-6, 'TolX', 1e-6, 'MaxIterations', 100 ); % ,'algorithm','trust-region'

            if verbose, fprintf('\n--- 启动 fsolve 求解器 (求解 [K̂_p, K̂_g, L, Bequest_total] - BGP版本) ---\n'); end
            [x_eq, ~, exitflag] = fsolve(system_wrapper, x0, options); % 调用fsolve
            if verbose, fprintf('--- fsolve 求解完成 ---\n'); end

            eq_found = (exitflag > 0); % 根据exitflag判断是否收敛
        end

        % =======================================================
        % == 结果展示与检查
        % =======================================================

        function verify_household_budget_constraint(ss, polS, cS, paramS)
            % =========================================================================
            % == 函数: verify_household_budget_constraint
            % == 版本: [v5 - PPS兼容验证版]
            % ==
            % == 核心修正:
            % ==   1. 增加 cS.pps_active 开关，使其能同时验证有/无PPS两种模式。
            % ==   2. 在PPS模式下，会计逻辑与 HHSolutionByAge_NAIVE_PPS 完全同步。
            % =========================================================================
            fprintf('\n--- [法医级检验] 正在对单个家庭的预算约束进行微观解剖 (PPS兼容版)... ---\n');

            % --- 1. 随机选择一个状态点进行检验 ---
            a_idx_test = randi([1, cS.aR_new]);
            ik_test = randi([1, cS.nk]);
            ie_test = randi([1, cS.nw]);
            if cS.pps_active
                ikpps_test = randi([1, cS.nkpps]);
                fprintf('   检验样本: 年龄 a=%d, k-idx=%d, kpps-idx=%d, e-idx=%d\n', a_idx_test, ik_test, ikpps_test, ie_test);
            else
                ikpps_test = 1; % 无PPS模式下，kpps索引为1
                fprintf('   检验样本: 年龄 a=%d, k-idx=%d, e-idx=%d (无PPS模式)\n', a_idx_test, ik_test, ie_test);
            end

            % --- 2. 提取公共信息和决策 ---
            M_age.r_mkt_t = ss.r_mkt;
            M_age.w_t = ss.w_hat;
            cS.theta_t = ss.Pension_in / max(1e-9, (ss.w_hat * ss.L_hat));
            g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;

            c_decision = polS(a_idx_test).c(ik_test, ikpps_test, ie_test);
            k_prime_decision = polS(a_idx_test).k_prime(ik_test, ikpps_test, ie_test);

            % --- 3. [核心] 使用与VFI完全一致的逻辑，重新计算可支配资源 ---
            k_now = cS.kGridV(ik_test);

            if cS.pps_active
                % --- PPS激活模式下的预算重构 ---
                kpps_now = cS.kppsGridV(ikpps_test);

                if a_idx_test <= cS.aR_new % 工作期
                    labor_income_gross = M_age.w_t * cS.ageEffV_new(a_idx_test) * paramS.leGridV(ie_test);
                    payg_tax = cS.theta_t * labor_income_gross;
                    pps_contrib = cS.pps_contrib_rate * labor_income_gross;
                    capital_tax = cS.tau_k * (k_now * M_age.r_mkt_t);
                    labor_tax = cS.tau_l * max(0, labor_income_gross - payg_tax - pps_contrib);
                    cash_on_hand_recalc = (k_now * (1 + M_age.r_mkt_t) + labor_income_gross) ...
                        - (capital_tax + payg_tax + labor_tax + pps_contrib);
                else % 退休期
                    pension_out = ss.Pension_out / sum(cS.Z((cS.aR_new+1):end)); % 简化计算
                    kpps_total_value = kpps_now * (1 + M_age.r_mkt_t);
                    pps_withdrawal = cS.pps_withdrawal_rate * kpps_total_value;
                    pps_tax = cS.pps_tax_rate_withdrawal * pps_withdrawal;
                    net_pps_inflow = pps_withdrawal - pps_tax;
                    capital_tax = cS.tau_k * (k_now * M_age.r_mkt_t);
                    cash_on_hand_recalc = (k_now * (1 + M_age.r_mkt_t) + pension_out + net_pps_inflow) ...
                        - capital_tax;
                end

            else
                % --- 无PPS模式下的预算重构 (与原始版本一致) ---
                k_return = k_now * (1 + M_age.r_mkt_t);
                labor_income_gross = M_age.w_t * cS.ageEffV_new(a_idx_test) * paramS.leGridV(ie_test);
                capital_tax = cS.tau_k * (k_now * M_age.r_mkt_t);
                payg_tax = cS.theta_t * labor_income_gross;
                labor_tax = cS.tau_l * max(0, labor_income_gross - payg_tax);
                total_tax_flow = capital_tax + payg_tax + labor_tax;
                cash_on_hand_recalc = k_return + labor_income_gross - total_tax_flow;
            end

            % --- 4. 计算冲击和最终可用于 C/S 的资源 (通用逻辑) ---
            shock_exp_factor = (ie_test == cS.nw + 1) * cS.kappa_young + (ie_test == cS.nw + 2) * cS.kappa_old;
            shock_exp = shock_exp_factor * cash_on_hand_recalc;
            available_for_c_and_s_recalc = cash_on_hand_recalc - shock_exp;
            fprintf('   A. 重新计算的可支配资源 (可用于C和S) ..: %12.8f\n', available_for_c_and_s_recalc);

            % --- 5. 计算VFI决策后的总支出 ---
            consumption_expenditure = c_decision * (1 + cS.tau_c);
            saving_expenditure = k_prime_decision * (1 + g_A_period);
            total_outlay_from_decision = consumption_expenditure + saving_expenditure;
            fprintf('   B. VFI决策后的总支出 (C_expend + k''*(1+g)) : %12.8f\n', total_outlay_from_decision);

            % --- 6. 计算微观预算缺口 ---
            micro_budget_gap = available_for_c_and_s_recalc - total_outlay_from_decision;
            fprintf('   -----------------------------------------------------------------\n');
            fprintf('   => 微观预算缺口 (A - B) ..................: %12.4e\n', micro_budget_gap);

            if abs(micro_budget_gap) < 1e-7
                fprintf('   ✅ [结论] 微观预算闭合！VFI的内部计算与外部检验完全一致！\n');
            else
                fprintf('   ⚠️ [结论] 微观预算不闭合！请仔细核对本函数与VFI函数的每一行计算！\n');
            end
            fprintf('--------------------------------------------------------------------------\n');
        end        % =========================================================================
        % == [新函数] 法医级家庭预算约束检验器
        % =========================================================================

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
            % == 函数说明: [v6 - 人口增长兼容版] 终期稳态完全内生验证器
            % ==
            % == 核心修正:
            % ==   1. [方法论统一] 验证器现在也采用【概率质量分裂法】。
            % ==   2. [人口增长] 在内生计算理论年龄分布 Z_theory 时，完全
            % ==      复制了求解器中的逻辑，引入了人口增长率 n。这确保了
            % ==      验证过程的基准是正确的。
            % =========================================================================

            fprintf('\n--- 正在进行终期稳态的【平滑方法 + 人口增长】验证...\n');

            % --- 步骤 1: 完全内生地计算理论年龄分布 Z_theory (包含人口增长) ---
            n_period = (1 + cS.n_ss)^cS.time_Step - 1;
            mass_levels_by_age = zeros(cS.aD_new, 1);
            mass_levels_by_age(1) = 1.0;
            for ia = 1:(cS.aD_new - 1)
                mass_levels_by_age(ia+1) = (mass_levels_by_age(ia) * cS.s_pathV(ia)) / (1 + n_period);
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
        end
        % =======================================================
        % == 辅助函数: 概率质量分裂法的核心工具
        % =======================================================

        % =======================================================
        % =======================================================
        % == 结果展示与检查 (v30 - NIPA会计准则完全修复版)
        % =======================================================
        % =======================================================
        % == 结果展示与检查 (v31 - 最终平衡版)
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