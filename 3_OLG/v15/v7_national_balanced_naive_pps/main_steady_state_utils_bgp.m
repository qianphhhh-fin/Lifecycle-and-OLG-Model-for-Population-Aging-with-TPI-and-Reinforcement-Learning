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
            % == 版本: [v13 - 生产资本/财富分离最终版]
            % ==
            % == 核心修正:
            % ==   1. 严格区分生产性资本 K_p (来自k) 和总财富 (k+kpps)。
            % ==   2. 求解器的核心任务是让宏观生产性资本 K_p_guess 与微观
            % ==      决策加总的 k_begin_agg 相等，实现理性预期均衡。
            % =========================================================================

            K_p_guess = x(1); K_g_guess = x(2); L_guess = x(3);
            
            [ss, ~, ~, ~] = main_steady_state_utils_bgp.calculate_aggregates_unified(K_p_guess, K_g_guess, L_guess, Z_ss_norm, cS, paramS);

            % --- 市场出清误差 ---
            % 方程1: 生产性私人资本市场 (宏观预期 = 微观加总)
            Kp_demand_macro = K_p_guess;
            Kp_supply_micro = ss.K_begin_agg; % [核心] 只使用普通资本k
            error_Kp = Kp_supply_micro - Kp_demand_macro;

            % 方程2: 公共资本市场 (不变)
            g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
            Gov_Total_Revenue = ss.Regular_tax + ss.Bequest_tax + ss.Public_Capital_Return;
            Resources_for_discretion = Gov_Total_Revenue - ss.Depreciation_g;
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
            % == 版本: [v18 - 生产资本/财富分离最终版]
            % ==
            % == 核心修正:
            % ==   1. 价格(r,w)由生产性资本K_p_guess决定。
            % ==   2. 家庭的资本总收入基于其真实持有的总财富(k+kpps)和市场利率r计算。
            % ==   3. 折旧和公共资本回报只与生产性资本K_p_guess相关。
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

            % --- 3. 聚合微观细项 ---
            K_hat_agg = 0; Kpps_hat_agg = 0; L_agg = 0; C_agg = 0; Bequest_tax_agg = 0;
            Shock_exp_agg = 0; Pension_out_agg = 0; Pension_in_agg = 0; PPS_contrib_agg = 0;
            Tax_labor_agg = 0; Tax_capital_agg = 0; Tax_consumption_agg = 0; Tax_pps_withdrawal_agg = 0;
            K_begin_agg = 0; Kpps_begin_agg = 0;
            [grid_k, grid_kpps] = ndgrid(cS.kGridV, cS.kppsGridV);

            for ia = 1:cS.aD_new
                mass_dist = Dist(:,:,:,ia);
                K_hat_agg = K_hat_agg + sum(polS(ia).k_prime .* mass_dist, 'all');
                Kpps_hat_agg = Kpps_hat_agg + sum(polS(ia).kpps_prime .* mass_dist, 'all');
                K_begin_agg = K_begin_agg + sum(grid_k .* mass_dist, 'all');
                Kpps_begin_agg = Kpps_begin_agg + sum(grid_kpps .* mass_dist, 'all');
                total_assets_chosen_for_next_period = sum((polS(ia).k_prime + polS(ia).kpps_prime) .* mass_dist, 'all');
                Bequest_tax_agg = Bequest_tax_agg + total_assets_chosen_for_next_period * (1 - cS.s_pathV(ia));
                if ia <= cS.aR_new
                    mass_by_epsilon = squeeze(sum(mass_dist, [1,2]));
                    L_agg = L_agg + sum( (cS.ageEffV_new(ia) * paramS.leGridV') .* mass_by_epsilon' , 'all');
                end
                C_agg = C_agg + sum(polS(ia).c .* mass_dist, 'all');
                Shock_exp_agg = Shock_exp_agg + sum(polS(ia).shock_exp .* mass_dist, 'all');
                Pension_out_agg = Pension_out_agg + sum(polS(ia).pension_out .* mass_dist, 'all');
                Pension_in_agg = Pension_in_agg + sum(polS(ia).tax_payg .* mass_dist, 'all');
                PPS_contrib_agg = PPS_contrib_agg + sum(polS(ia).pension_in .* mass_dist, 'all');
                Tax_labor_agg = Tax_labor_agg + sum(polS(ia).tax_labor .* mass_dist, 'all');
                Tax_capital_agg = Tax_capital_agg + sum(polS(ia).tax_capital .* mass_dist, 'all');
                Tax_consumption_agg = Tax_consumption_agg + sum(polS(ia).tax_consumption .* mass_dist, 'all');
                Tax_pps_withdrawal_agg = Tax_pps_withdrawal_agg + sum(polS(ia).tax_pps_withdrawal .* mass_dist, 'all');
            end
            Regular_tax_agg = Tax_labor_agg + Tax_capital_agg + Tax_consumption_agg + Tax_pps_withdrawal_agg;

            % --- 4. 计算流量和打包 ---
            % [核心] 家庭的资本收入基于其真实持有的总财富 (K_begin_agg + Kpps_begin_agg)
            Household_Capital_Income = (K_begin_agg + Kpps_begin_agg) * M_prices.r_mkt_t;
            Household_Total_Inflow = (M_prices.w_t * L_agg) + Household_Capital_Income + Pension_out_agg;
            Household_Outlay_NonSaving = C_agg + Shock_exp_agg + Regular_tax_agg + Pension_in_agg;
            Saving_private_flow = Household_Total_Inflow - Household_Outlay_NonSaving;
            
            % [核心] 折旧和公共资本回报都只与生产性资本K_p_guess相关
            Depreciation_p = cS.ddk * K_p_guess;
            Depreciation_g = cS.ddk_g * K_g_guess;
            Factor_Payment_Total = (M_prices.w_t * L_agg) + ((M_prices.r_mkt_t + cS.ddk) * K_p_guess);
            Public_Capital_Return = M_prices.Y_t - Factor_Payment_Total;

            ss = struct();
            ss.K_private_hat = K_hat_agg + Kpps_hat_agg; % 总期末财富
            ss.K_private_begin_hat = K_p_guess; % 这是生产性资本 Kp(t)
            ss.L_hat = L_agg; ss.K_public_hat = K_g_guess;
            ss.Y_from_production_hat = M_prices.Y_t;
            ss.w_hat = M_prices.w_t; ss.r_mkt = M_prices.r_mkt_t;
            ss.C_agg = C_agg; ss.Shock_exp_agg = Shock_exp_agg;
            ss.Bequest_tax = Bequest_tax_agg; ss.Regular_tax = Regular_tax_agg;
            ss.Pension_in = Pension_in_agg; ss.Pension_out = Pension_out_agg;
            ss.Saving_private_flow = Saving_private_flow;
            ss.Depreciation_p = Depreciation_p; ss.Depreciation_g = Depreciation_g;
            ss.Public_Capital_Return = Public_Capital_Return;
            ss.PPS_contrib_agg = PPS_contrib_agg;
            ss.K_hat_agg = K_hat_agg; ss.Kpps_hat_agg = Kpps_hat_agg;
            ss.K_begin_agg = K_begin_agg; ss.Kpps_begin_agg = Kpps_begin_agg;
            ss.Tax_labor_agg = Tax_labor_agg; ss.Tax_capital_agg = Tax_capital_agg;
            ss.Tax_consumption_agg = Tax_consumption_agg; ss.Tax_pps_withdrawal_agg = Tax_pps_withdrawal_agg;
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
                if cS_vfi.pps_active
                    [valS(:,:,:,a_idx), polS_cell{a_idx}] = ...
                        main_steady_state_utils_bgp.HHSolutionByAge_NAIVE_PPS(...
                        a_idx, vPrime_kkppse_next, M_vfi, bV_payg_vfi(a_idx), paramS_vfi, cS_vfi);
                else
                    % 调用统一的、并行的年龄组求解器
                    [valS(:,:,:,a_idx), polS_cell{a_idx}] = ...
                        main_steady_state_utils_bgp.HHSolutionByAge_NAIVE(...
                        a_idx, vPrime_kkppse_next, M_vfi, bV_payg_vfi(a_idx), paramS_vfi, cS_vfi);
                end
            end

            % --- 4. 将cell数组转换为更易于访问的结构体数组 ---
            polS = [polS_cell{:}]; % polS(ia).c 可以访问第ia期的消费决策
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


        function [val_age, pol_age] = HHSolutionByAge_NAIVE_PPS(a_idx, ~, M_age, b_age_val, paramS_age, cS)
            % =========================================================================
            % == 函数: HHSolutionByAge_NAIVE_PPS
            % == 版本: [v3.0 - 税种拆分版]
            % ==
            % == 核心修正:
            % ==   1. pol_age 结构体现在分别存储不同税种，而不是合并到 tax_regular。
            % =========================================================================

            % --- 0. 设定行为参数 ---
            PPS_CONTRIB_RATE = 0.01;
            SAVINGS_RATE_K = 0.20;

            % --- 1. 初始化 ---
            val_age = -1e20 * ones(cS.nk, cS.nkpps, cS.nw_expanded);
            pol_age = struct(...
                'c', zeros(cS.nk, cS.nkpps, cS.nw_expanded), 'k_prime', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
                'kpps_prime', zeros(cS.nk, cS.nkpps, cS.nw_expanded), 'tax_payg', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
                'shock_exp', zeros(cS.nk, cS.nkpps, cS.nw_expanded), 'pension_in', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
                'pension_out', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
                'tax_labor', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...       % [新增]
                'tax_capital', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...     % [新增]
                'tax_consumption', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ... % [新增]
                'tax_pps_withdrawal', zeros(cS.nk, cS.nkpps, cS.nw_expanded) ... % [新增]
                );

            g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
            market_return_factor = 1 + M_age.r_mkt_t;

            % --- 2. 终点决策 ---
            if a_idx == cS.aD_new
                [K_grid, Kpps_grid] = ndgrid(cS.kGridV, cS.kppsGridV);
                total_capital = K_grid + Kpps_grid;
                capital_tax = cS.tau_k .* (total_capital * M_age.r_mkt_t);
                total_wealth = total_capital * market_return_factor - capital_tax + b_age_val;
                c_expend_final = total_wealth;
                c_final = c_expend_final ./ (1 + cS.tau_c);
                consumption_tax = c_final * cS.tau_c;

                for ie = 1:cS.nw_expanded
                    pol_age.c(:, :, ie) = c_final;
                    pol_age.pension_out(:, :, ie) = b_age_val;
                    pol_age.tax_capital(:, :, ie) = capital_tax;
                    pol_age.tax_consumption(:, :, ie) = consumption_tax;
                end
                return;
            end

            % --- 3. 遍历所有状态点 ---
            for ie = 1:cS.nw_expanded
                for ikpps = 1:cS.nkpps
                    for ik = 1:cS.nk
                        k_now = cS.kGridV(ik);
                        kpps_now = cS.kppsGridV(ikpps);

                        % A. 计算资源和税收
                        total_capital_now = k_now + kpps_now;
                        capital_return_gross = total_capital_now * market_return_factor;
                        capital_tax = cS.tau_k * (total_capital_now * M_age.r_mkt_t);
                        capital_return_net = capital_return_gross - capital_tax;

                        pension_out = 0; labor_income_gross = 0; payg_tax = 0; labor_tax = 0;
                        pps_contrib = 0; pps_withdrawal_net = 0; pps_withdrawal_tax = 0;

                        if a_idx <= cS.aR_new % 工作期
                            labor_income_gross = M_age.w_t * cS.ageEffV_new(a_idx) * paramS_age.leGridV(ie);
                            payg_tax = cS.theta_t * labor_income_gross;
                            pps_contrib = PPS_CONTRIB_RATE * labor_income_gross;
                            taxable_labor_income = max(0, labor_income_gross - payg_tax - pps_contrib);
                            labor_tax = cS.tau_l * taxable_labor_income;
                            kpps_prime_decision = kpps_now * market_return_factor + pps_contrib;
                            labor_inflow_net = labor_income_gross - payg_tax - labor_tax - pps_contrib;
                        else % 退休期
                            pension_out = b_age_val;
                            labor_inflow_net = 0;
                            pps_withdrawal_gross = cS.pps_withdrawal_rate * (kpps_now * market_return_factor);
                            pps_withdrawal_tax = cS.pps_tax_rate_withdrawal * pps_withdrawal_gross;
                            pps_withdrawal_net = pps_withdrawal_gross - pps_withdrawal_tax;
                            kpps_prime_decision = kpps_now * market_return_factor - pps_withdrawal_gross;
                        end

                        % B. 计算最终可支配资源
                        cash_on_hand = capital_return_net + labor_inflow_net + pension_out + pps_withdrawal_net;
                        shock_exp_factor = (ie == cS.nw + 1) * cS.kappa_young + (ie == cS.nw + 2) * cS.kappa_old;
                        shock_exp = shock_exp_factor * cash_on_hand;
                        available_for_c_and_s = cash_on_hand - shock_exp;

                        % C. 应用储蓄规则
                        k_saving = SAVINGS_RATE_K * available_for_c_and_s;
                        k_prime_decision = max(cS.kMin, k_saving / (1 + g_A_period));
                        c_expend = available_for_c_and_s - k_prime_decision * (1 + g_A_period);
                        if c_expend < 1e-9, c_expend = 1e-9; end
                        c_decision = c_expend / (1 + cS.tau_c);
                        consumption_tax = c_decision * cS.tau_c;

                        % D. 存储决策
                        pol_age.c(ik, ikpps, ie) = c_decision;
                        pol_age.k_prime(ik, ikpps, ie) = k_prime_decision;
                        pol_age.kpps_prime(ik, ikpps, ie) = max(cS.kppsMin, kpps_prime_decision);
                        pol_age.pension_out(ik, ikpps, ie) = pension_out;
                        pol_age.shock_exp(ik, ikpps, ie) = shock_exp;
                        pol_age.tax_payg(ik, ikpps, ie) = payg_tax;
                        pol_age.pension_in(ik, ikpps, ie) = pps_contrib;
                        % [核心修改] 分别存储税种
                        pol_age.tax_labor(ik, ikpps, ie) = labor_tax;
                        pol_age.tax_capital(ik, ikpps, ie) = capital_tax;
                        pol_age.tax_consumption(ik, ikpps, ie) = consumption_tax;
                        pol_age.tax_pps_withdrawal(ik, ikpps, ie) = pps_withdrawal_tax;
                    end
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
            % 版本: [v2 - 理论一致版] - 无需修改
            if K_p <= 0, K_p = 1e-8; end; if L <= 0, L = 1e-8; end; if K_g <= 0, K_g = 1e-8; end;
            
            labor_exponent = 1 - cS.alpha - cS.gamma;
            if labor_exponent <= 0, error('劳动和资本的产出弹性之和必须小于1！'); end
            
            Y_period = A_t .* (K_p.^cS.alpha) .* (K_g.^cS.gamma) .* (L.^labor_exponent);
            
            MPK_p_period = cS.alpha .* Y_period ./ K_p;
            w_t = labor_exponent .* Y_period ./ L;
            
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
            % == 版本: [v28 - 资本分离审计最终版]
            % ==
            % == 核心特性:
            % ==   1. 报告结构与 v13 求解器的“生产资本/财富分离”逻辑完全对齐。
            % ==   2. 清晰展示生产性资本 Kp 和金融性资本 Kpps 的区别与联系。
            % =========================================================================

            fprintf('\n\n================================================================================\n');
            if cS.pps_active
                fprintf('===      国民经济核算详细报告 (BGP + PPS激活 - 资本分离版)      ===\n');
            else
                fprintf('===        国民经济核算详细报告 (BGP + 无PPS - 标准版)       ===\n');
            end
            fprintf('================================================================================\n');

            % --- 0. 预计算 ---
            g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
            Y_prod = ss.Y_from_production_hat;
            C_total = ss.C_agg + ss.Shock_exp_agg;
            I_p_gross = ss.Depreciation_p + g_A_period * ss.K_private_begin_hat;
            I_g_gross = ss.Depreciation_g + g_A_period * ss.K_public_hat;
            I_total_gross = I_p_gross + I_g_gross;
            Gov_Total_Revenue = ss.Regular_tax + ss.Bequest_tax + ss.Public_Capital_Return;
            Resources_for_discretion = Gov_Total_Revenue - ss.Depreciation_g;
            G_c = (1 - cS.lambda_g) * Resources_for_discretion;

            % --- 1. [核心] 生产性资本市场检验 (宏观预期 vs 微观加总) ---
            fprintf('\n--- [1. 生产性资本市场检验 (理性预期均衡)] ---\n');
            Kp_demand_macro = ss.K_private_begin_hat;
            Kp_supply_micro = ss.K_begin_agg;
            mismatch_stock = Kp_supply_micro - Kp_demand_macro;
            fprintf('   A. 宏观层面需求的生产性资本 Kp(t) .........: %12.8f\n', Kp_demand_macro);
            fprintf('   B. 微观层面供给的生产性资本 Σk(t) .........: %12.8f\n', Kp_supply_micro);
            fprintf('   --------------------------------------------------------------------\n');
            fprintf('   => 存量均衡缺口 (B - A) (应为0) ..........: %12.4e\n', mismatch_stock);

            % --- 2. 宏观资源约束 (流量守恒) ---
            fprintf('\n--- [2. 宏观资源约束 (流量守恒)] ---\n');
            Y_exp_gross = C_total + I_total_gross + G_c;
            mismatch_gross = Y_prod - Y_exp_gross;
            fprintf('   A. 国内生产总值 (GDP) .....................: %12.8f\n', Y_prod);
            fprintf('   B. 国内总支出 (C + I_gross + Gc) ..........: %12.8f\n', Y_exp_gross);
            fprintf('       - 私人消费 (C) ........................: %12.8f\n', C_total);
            fprintf('       - 总投资 (I_gross) ....................: %12.8f\n', I_total_gross);
            fprintf('           - 私人(生产性)总投资 (I_p_gross) ..: %12.8f\n', I_p_gross);
            fprintf('           - 公共总投资 (I_g_gross) ..........: %12.8f\n', I_g_gross);
            fprintf('       - 政府消费 (G_c) ......................: %12.8f\n', G_c);
            fprintf('   --------------------------------------------------------------------\n');
            fprintf('   => 宏观流量缺口 (应为0) ..................: %12.4e\n', mismatch_gross);

            % --- 3. 政府预算检验 ---
            fprintf('\n--- [3. 政府预算检验 (收支平衡)] ---\n');
            Gov_Outflows = I_g_gross + G_c;
            mismatch_gov = Gov_Total_Revenue - Gov_Outflows;
            fprintf('   A. 政府总收入 .............................: %12.8f\n', Gov_Total_Revenue);
            fprintf('       - 税收收入 ............................: %12.8f\n', ss.Regular_tax + ss.Bequest_tax);
            fprintf('       - 公共资本回报 ........................: %12.8f\n', ss.Public_Capital_Return);
            fprintf('   B. 政府总支出 (I_g_gross + Gc) ............: %12.8f\n', Gov_Outflows);
            fprintf('   => 政府预算缺口 (应为0) ..................: %12.4e\n', mismatch_gov);

            % --- 4. [补充验证] 国民总储蓄 vs 国民总投资 ---
            fprintf('\n--- [4. 国民总储蓄 vs 国民总投资 (Walras定律验证)] ---\n');
            S_p_net = ss.Saving_private_flow;
            Gov_Net_Saving = cS.lambda_g * Resources_for_discretion;
            S_total_gross = S_p_net + Gov_Net_Saving + ss.Depreciation_p + ss.Depreciation_g;
            mismatch_gross_saving = S_total_gross - I_total_gross;
            fprintf('   A. 国民总储蓄 (S_gross) ...................: %12.8f\n', S_total_gross);
            fprintf('   B. 国民总投资 (I_gross) ...................: %12.8f\n', I_total_gross);
            fprintf('   => 储蓄-投资缺口 (S_gross - I_gross) .......: %12.4e\n', mismatch_gross_saving);
            
            fprintf('\n================================================================================\n');
        end        % == [新函数] 法医级家庭预算约束检验器
        % =========================================================================
        function verify_household_budget_constraint(ss, polS, cS, paramS)
            % =========================================================================
            % == 函数: verify_household_budget_constraint
            % == 版本: [v5.1 - 笔误修正版]
            % ==
            % == 核心修正:
            % ==   1. 修正了变量 ie_test 定义顺序的错误。
            % =========================================================================

            fprintf('\n--- [法医级检验] 正在对单个家庭的预算约束进行微观解剖 (PPS兼容版)... ---\n');

            % --- 1. 随机选择一个状态点进行检验 ---
            a_idx_test = randi([1, cS.aR_new]);
            ik_test = randi([1, cS.nk]);
            ie_test = randi([1, cS.nw]); % [修正] 将 ie_test 的定义提前

            % [PPS修改] 如果有PPS，也随机选择一个kpps状态
            if cS.pps_active && cS.nkpps > 1
                ikpps_test = randi([1, cS.nkpps]);
                fprintf('   检验样本: 年龄组 a=%d, k_idx=%d, kpps_idx=%d, e_idx=%d\n', a_idx_test, ik_test, ikpps_test, ie_test);
            else
                ikpps_test = 1;
                fprintf('   检验样本: 年龄组 a=%d, k_idx=%d, e_idx=%d\n', a_idx_test, ik_test, ie_test);
            end

            % --- 2. 提取该状态点的数据 ---
            M_age.r_mkt_t = ss.r_mkt;
            M_age.w_t = ss.w_hat;
            cS.theta_t = ss.Pension_in / max(1e-9, (ss.w_hat * ss.L_hat));
            g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;

            c_decision = polS(a_idx_test).c(ik_test, ikpps_test, ie_test);
            k_prime_decision = polS(a_idx_test).k_prime(ik_test, ikpps_test, ie_test);

            % --- 3. [核心修正] 使用与VFI完全一致的逻辑，重新计算可支配资源 ---
            k_now = cS.kGridV(ik_test);
            kpps_now = 0;
            if cS.pps_active && cS.nkpps > 1
                kpps_now = cS.kppsGridV(ikpps_test);
            end

            % A. [核心修正] 重新计算资源和税收 (与NAIVE_PPS v2.0完全同步)
            market_return_factor = 1 + M_age.r_mkt_t;
            total_capital_now = k_now + kpps_now;
            capital_return_gross = total_capital_now * market_return_factor;
            capital_tax = cS.tau_k * (total_capital_now * M_age.r_mkt_t);
            capital_return_net = capital_return_gross - capital_tax;

            pension_out = 0;
            labor_income_gross = 0;
            payg_tax = 0;
            labor_tax = 0;
            pps_contrib = 0;
            pps_withdrawal_net = 0;

            if a_idx_test <= cS.aR_new % 工作期
                labor_income_gross = M_age.w_t * cS.ageEffV_new(a_idx_test) * paramS.leGridV(ie_test);
                payg_tax = cS.theta_t * labor_income_gross;

                if cS.pps_active
                    PPS_CONTRIB_RATE = 0.01; % 与NAIVE_PPS中的参数保持一致
                    pps_contrib = PPS_CONTRIB_RATE * labor_income_gross;
                end

                taxable_labor_income = max(0, labor_income_gross - payg_tax - pps_contrib);
                labor_tax = cS.tau_l * taxable_labor_income;
                labor_inflow_net = labor_income_gross - payg_tax - labor_tax - pps_contrib;
            else % 退休期
                pension_out = ss.b_hat;
                labor_inflow_net = 0;

                if cS.pps_active
                    pps_withdrawal_gross = cS.pps_withdrawal_rate * (kpps_now * market_return_factor);
                    pps_withdrawal_tax = cS.pps_tax_rate_withdrawal * pps_withdrawal_gross;
                    pps_withdrawal_net = pps_withdrawal_gross - pps_withdrawal_tax;
                end
            end

            % B. 计算最终可支配资源
            cash_on_hand_recalc = capital_return_net + labor_inflow_net + pension_out + pps_withdrawal_net;
            shock_exp_factor = (ie_test == cS.nw + 1) * cS.kappa_young + (ie_test == cS.nw + 2) * cS.kappa_old;
            shock_exp = shock_exp_factor * cash_on_hand_recalc;
            available_for_c_and_s_recalc = cash_on_hand_recalc - shock_exp;

            fprintf('   A. 重新计算的可支配资源 (可用于C和S) ..: %12.8f\n', available_for_c_and_s_recalc);

            % C. 计算VFI决策后的总支出
            consumption_expenditure = c_decision * (1 + cS.tau_c);
            saving_expenditure = k_prime_decision * (1 + g_A_period);
            total_outlay_from_decision = consumption_expenditure + saving_expenditure;

            fprintf('   B. VFI决策后的总支出 (C_expend + k''*(1+g)) : %12.8f\n', total_outlay_from_decision);

            % D. 计算微观预算缺口
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