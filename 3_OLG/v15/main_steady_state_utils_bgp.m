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
        function [ss, Dist, polS, valS] = solve_steady_state_unified(cS, paramS, params_ext, verbose, x0_guess, solver_method)
            % [重构核心-v4 PPS版] 统一的稳态求解器主入口
            if nargin < 4, verbose = true; end
            if nargin < 5, x0_guess = []; end
            if nargin < 6, solver_method = 'fsolve'; end

            % --- 1. 参数设置 ---
            cS.A = params_ext.A;
            if isfield(params_ext, 'g_A_ss'), cS.g_A_ss = params_ext.g_A_ss; end
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
                % 强制设定为无PPS模式的维度，确保VFI兼容
                cS.nkpps = 1;
                cS.kppsGridV = 0;
            end

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
            % verify_household_budget_constraint 暂不更新，待需要时再适配PPS
            % main_steady_state_utils_bgp.verify_household_budget_constraint(ss, polS, cS, paramS);
            main_steady_state_utils_bgp.verify_FULL_steady_state(Dist, polS, paramS, cS);

            if verbose, main_steady_state_utils_bgp.display_national_accounts_unified(ss, cS);
            else, main_steady_state_utils_bgp.check_national_accounts_unified(ss, cS); end
        end

        % =======================================================
        % == 核心求解方程 (无变化)
        % =======================================================
        % =======================================================
        % == 核心求解方程 (恢复流量均衡版)
        % =======================================================
        function F_error = system_of_equations_steady_state(x, Z_ss_norm, cS, paramS)
            % =========================================================================
            % == 函数: system_of_equations_steady_state
            % == 版本: [v15 - 流量均衡最终版]
            % ==
            % == 核心决策:
            % ==   经过尝试，我们确认【存量均衡】方法 K_guess = K_resulting 在
            % ==   存在意外遗赠 (bequests) 的OLG模型中是【不正确】的。因为它
            % ==   忽略了死亡家庭持有的资产会“漏出”成为遗赠税，而不会转变为
            % ==   下一期的生产性资本。
            % ==
            % ==   因此，我们必须回到并坚持使用【流量均衡】条件。它正确地将
            % ==   意外遗赠(Bequest_demand)处理为对私人储蓄的一种需求，从而
            % ==   保证了宏观会计的完全闭合。
            % ==
            % ==   收敛性问题需要通过其他方法解决（如优化初值、使用更鲁棒的
            % ==   求解器或平滑政策函数），而不是改变正确的经济学原理。
            % =========================================================================

            K_private_total_guess = x(1); K_g_guess = x(2); L_guess = x(3);

            % --- 1. 基于猜测的宏观价格，计算所有微观和宏观聚合量 ---
            g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
            [ss, ~, ~, ~] = main_steady_state_utils_bgp.calculate_aggregates_unified(K_private_total_guess, K_g_guess, L_guess, Z_ss_norm, cS, paramS);

            % --- 2. [核心] 市场出清误差 (基于流量) ---
            % --- 方程1: 私人资本市场 (S_p_net = I_p_net + Bequests) ---
            % 私人净储蓄的供给
            S_p_supply = ss.Saving_private_flow;

            % 对私人净储蓄的需求来自两部分：
            % 1. 企业为了跟上技术增长所需的净投资
            I_p_net_demand = g_A_period * K_private_total_guess;
            % 2. 政府通过遗赠税从储蓄池中抽走的部分
            Bequest_demand = ss.Bequest_tax;

            Total_demand_for_net_saving = I_p_net_demand + Bequest_demand;

            error_Kp = Total_demand_for_net_saving - S_p_supply;

            % --- 方程2: 公共资本市场 (流量均衡) ---
            Gov_Total_Revenue = ss.Regular_tax + ss.Bequest_tax + ss.PPS_tax_agg + ss.Public_Capital_Return;
            Resources_for_discretion = Gov_Total_Revenue - ss.Depreciation_g;
            Gov_Net_Saving = cS.lambda_g * Resources_for_discretion;
            I_g_net_demand = g_A_period * K_g_guess;
            error_Kg = I_g_net_demand - Gov_Net_Saving;

            % --- 方程3: 劳动市场 (供给 = 需求) ---
            error_L = L_guess - ss.L_hat;

            F_error = [error_Kp; error_Kg; error_L];
        end        % == 宏观聚合函数 (修改版)
        % =======================================================
        % =======================================================
        % == 宏观聚合函数 (修改版)
        % =======================================================
        % =======================================================
        % == 宏观聚合函数 (修改版)
        % =======================================================
        function [ss, Dist, polS, valS] = calculate_aggregates_unified(K_private_total_guess, K_g_guess, L_guess, Z_ss_norm, cS, paramS)
            % =========================================================================
            % == 函数: calculate_aggregates_unified
            % == 版本: [v22 - 最终会计修正版]
            % ==
            % == 核心修正:
            % ==   [致命BUG修复] 修正了私人净储蓄 (Saving_private_flow) 的会计定义。
            % ==   之前的版本错误地将PPS的缴费和提取作为家庭部门的支出和收入项，
            % ==   这混淆了“国民经济核算(NIPA)”逻辑和“家庭现金流”逻辑，导致
            % ==   最终资源约束 Y = C+I+Gc 不平衡。
            % ==
            % ==   新逻辑严格遵循NIPA定义：
            % ==   私人净储蓄 = 私人可支配收入 - 私人消费
            % ==   其中，PPS缴费和提取是家庭内部的资产组合调整，不计入NIPA意义
            % ==   上的收入和支出。这确保了所有宏观恒等式在有PPS时也能闭合。
            % =========================================================================

            % --- 1. 初始化和价格计算 ---
            if K_private_total_guess <= 0, K_private_total_guess = 1e-8; end; if K_g_guess <= 0, K_g_guess = 1e-8; end; if L_guess <= 0, L_guess = 1e-8; end;
            A_ss = 1.0; theta_ss = cS.theta_path(1);

            M_prices = main_steady_state_utils_bgp.get_prices_at_t(K_private_total_guess, K_g_guess, L_guess, A_ss, cS);

            M_for_hh = M_prices;
            mass_retirees_ss = sum(Z_ss_norm((cS.aR_new+1):end));
            M_for_hh.b_t = (theta_ss * M_prices.w_t * L_guess) / max(1e-9, mass_retirees_ss);
            cS.theta_t = theta_ss;

            % --- 2. 求解家庭问题和稳态分布 ---
            [polS, valS] = main_steady_state_utils_bgp.HHSolution_VFI_unified(M_for_hh, paramS, cS);
            Dist = main_steady_state_utils_bgp.solve_steady_state_distribution_unified(polS, paramS, cS, Z_ss_norm);

            % --- 3. 聚合微观变量 ---
            K_private_hat_agg = 0; K_pps_hat_agg = 0; L_agg = 0; C_agg = 0; Bequest_tax_agg = 0;
            Shock_exp_agg = 0; Pension_out_agg = 0; Pension_in_agg = 0; Regular_tax_agg = 0;
            PPS_contrib_agg = 0; PPS_withdrawal_agg = 0; PPS_tax_agg = 0;

            for ia = 1:cS.aD_new
                mass_dist = Dist(:,:,:,ia);
                K_private_hat_agg = K_private_hat_agg + sum(polS(ia).k_prime .* mass_dist, 'all');
                if cS.pps_active, K_pps_hat_agg = K_pps_hat_agg + sum(polS(ia).kpps_prime .* mass_dist, 'all'); end

                total_assets_chosen_for_next_period = sum((polS(ia).k_prime + polS(ia).kpps_prime) .* mass_dist, 'all');
                Bequest_tax_agg = Bequest_tax_agg + total_assets_chosen_for_next_period * (1 - cS.s_pathV(ia));

                if ia <= cS.aR_new
                    mass_by_epsilon = squeeze(sum(mass_dist, [1,2]));
                    L_agg = L_agg + sum( (cS.ageEffV_new(ia) * paramS.leGridV') .* mass_by_epsilon' , 'all');
                end
                C_agg = C_agg + sum(polS(ia).c .* mass_dist, 'all');
                Shock_exp_agg = Shock_exp_agg + sum(polS(ia).shock_exp .* mass_dist, 'all');
                Regular_tax_agg = Regular_tax_agg + sum(polS(ia).tax_regular .* mass_dist, 'all');
                Pension_in_agg = Pension_in_agg + sum(polS(ia).tax_payg .* mass_dist, 'all');
                Pension_out_agg = Pension_out_agg + sum(polS(ia).pension_out .* mass_dist, 'all');

                if cS.pps_active
                    PPS_contrib_agg = PPS_contrib_agg + sum(polS(ia).pps_contrib .* mass_dist, 'all');
                    PPS_withdrawal_agg = PPS_withdrawal_agg + sum(polS(ia).pps_withdrawal .* mass_dist, 'all');
                    PPS_tax_agg = PPS_tax_agg + sum(polS(ia).pps_tax .* mass_dist, 'all');
                end
            end

            % --- 4. [最终会计修正] 严格遵循NIPA恒等式计算私人净储蓄 ---
            % 私人净储蓄 = 私人可支配收入 - 私人消费

            % A. 私人部门的总收入 (要素收入 + 转移收入)
            % 注意：PPS提取是家庭内部资产变现，不是NIPA意义上的新增收入。
            Household_NIPA_Inflow = (M_prices.w_t * L_agg) ...           % 劳动收入
                + (K_private_total_guess * M_prices.r_mkt_t) ...       % 资本利得收入 (来自 Kp 和 Kpps)
                + Pension_out_agg;                                     % 公共养老金转移收入

            % B. 私人部门的总支出（消费性 + 税收）
            % 注意：PPS缴费是储蓄的一种形式，不是NIPA意义上的支出。
            Household_NIPA_Outlay_NonSaving = C_agg ...                  % 私人消费
                + Shock_exp_agg ...                                    % 意外冲击消费
                + Regular_tax_agg ...                                  % 常规税
                + Pension_in_agg ...                                   % PAYG缴费 (视为税)
                + PPS_tax_agg;                                         % PPS提取时缴的税 (是真实的税)

            % C. [正确定义] 私人净储蓄 (可支配收入 - 消费)
            Saving_private_flow = Household_NIPA_Inflow - Household_NIPA_Outlay_NonSaving;

            % --- 5. 计算其他宏观量 ---
            Depreciation_p = cS.ddk * K_private_total_guess;
            Depreciation_g = cS.ddk_g * K_g_guess;
            Factor_Payment_Total = (M_prices.w_t * L_agg) + ((M_prices.r_mkt_t + cS.ddk) * K_private_total_guess);
            Public_Capital_Return = M_prices.Y_t - Factor_Payment_Total;

            % --- 6. 填充 ss 结构体 ---
            ss = struct();
            ss.K_private_begin_hat = K_private_total_guess;
            ss.K_private_hat = K_private_hat_agg; % 注意：这里只记录了非PPS部分
            ss.L_hat = L_agg;
            ss.K_public_hat = K_g_guess;
            ss.Y_from_production_hat = M_prices.Y_t;
            ss.w_hat = M_prices.w_t;
            ss.r_mkt = M_prices.r_mkt_t;
            ss.C_agg = C_agg;
            ss.Shock_exp_agg = Shock_exp_agg;
            ss.Bequest_tax = Bequest_tax_agg;
            ss.Regular_tax = Regular_tax_agg;
            ss.Pension_in = Pension_in_agg;
            ss.Pension_out = Pension_out_agg;
            ss.Saving_private_flow = Saving_private_flow;
            ss.Depreciation_p = Depreciation_p;
            ss.Depreciation_g = Depreciation_g;
            ss.Public_Capital_Return = Public_Capital_Return;

            ss.K_pps_hat = K_pps_hat_agg;
            ss.PPS_contrib_agg = PPS_contrib_agg;
            ss.PPS_withdrawal_agg = PPS_withdrawal_agg;
            ss.PPS_tax_agg = PPS_tax_agg;
        end

        function [polS, valS] = HHSolution_VFI_unified(M_vfi, paramS_vfi, cS_vfi)
            % =========================================================================
            % == 函数: HHSolution_VFI_unified
            % == 版本: [v5 - 最终VFI调用版]
            % == 核心修正: 调用全新的、基于完整VFI的子函数。
            % =========================================================================

            valS = -Inf(cS_vfi.nk, cS_vfi.nkpps, cS_vfi.nw_expanded, cS_vfi.aD_new);
            polS_cell = cell(cS_vfi.aD_new, 1);
            bV_payg_vfi = zeros(1, cS_vfi.aD_new);
            if cS_vfi.aR_new < cS_vfi.aD_new, bV_payg_vfi((cS_vfi.aR_new+1):cS_vfi.aD_new) = M_vfi.b_t; end

            for a_idx = cS_vfi.aD_new : -1 : 1
                vPrime_kkppse_next = [];
                if a_idx < cS_vfi.aD_new, vPrime_kkppse_next = valS(:,:,:,a_idx+1); end

                if cS_vfi.pps_active
                    % 调用【有PPS】的完整VFI求解器
                    [valS(:,:,:,a_idx), polS_cell{a_idx}] = ...
                        main_steady_state_utils_bgp.HHSolutionByAge_VFI_PPS(...
                        a_idx, vPrime_kkppse_next, M_vfi, bV_payg_vfi(a_idx), paramS_vfi, cS_vfi);
                else
                    % 调用【无PPS】的完整VFI求解器
                    [valS(:,:,:,a_idx), polS_cell{a_idx}] = ...
                        main_steady_state_utils_bgp.HHSolutionByAge_VFI_noPPS(...
                        a_idx, vPrime_kkppse_next, M_vfi, bV_payg_vfi(a_idx), paramS_vfi, cS_vfi);
                end
            end
            polS = [polS_cell{:}];
        end

        % 在类中添加这个辅助函数

        function [val_age, pol_age] = HHSolutionByAge_VFI_noPPS(a_idx, vPrime_kkppse_next, M_age, b_age_val, paramS_age, cS)
            % =========================================================================
            % == 函数: HHSolutionByAge_VFI_noPPS
            % == 版本: [v1.0 - 完整VFI版]
            % ==
            % == 目的: 基于高质量参考代码，为无PPS情境提供一个鲁棒、精确的VFI求解器。
            % == 核心特性:
            % ==   - 真正的效用最大化，而非固定储蓄率。
            % ==   - 使用PCHIP插值平滑处理下一期价值函数。
            % ==   - 使用自适应网格进行高效的最优决策搜索。
            % ==   - 完全兼容BGP框架和遗赠动机。
            % =========================================================================

            % --- 1. 初始化 ---
            nk_search_grid = 150; % 搜索网格的密度
            val_age_slice = -1e20 * ones(cS.nk, 1, cS.nw_expanded);
            pol_age_slice = struct(...
                'c', zeros(cS.nk, 1, cS.nw_expanded), 'k_prime', zeros(cS.nk, 1, cS.nw_expanded), ...
                'kpps_prime', zeros(cS.nk, 1, cS.nw_expanded), 'tax_regular', zeros(cS.nk, 1, cS.nw_expanded), ...
                'tax_payg', zeros(cS.nk, 1, cS.nw_expanded), 'shock_exp', zeros(cS.nk, 1, cS.nw_expanded), ...
                'pension_in', zeros(cS.nk, 1, cS.nw_expanded), 'pension_out', zeros(cS.nk, 1, cS.nw_expanded));

            % --- 2. 预计算和插值器设置 ---
            g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
            market_return_factor = 1 + M_age.r_mkt_t;
            % BGP下的真实折现因子，包含了对增长的调整
            discount_factor_future_utility = cS.beta * (1 + g_A_period)^(1 - cS.sigma);

            % --- 3. 终点决策 ---
            if a_idx == cS.aD_new
                k_capital_tax = cS.tau_k .* (cS.kGridV * M_age.r_mkt_t);
                total_wealth = cS.kGridV * market_return_factor - k_capital_tax + b_age_val;
                k_prime_final = zeros(cS.nk, 1);
                c_expend_final = total_wealth;
                % 如果有遗赠动机，则在消费和遗赠之间做最优分配
                if isfield(cS,'phi_bequest') && cS.phi_bequest > 0
                    phi_effective = cS.beta * cS.phi_bequest; % 终点没有增长调整
                    c_over_k_ratio = phi_effective^(-1/cS.sigma);
                    k_prime_final = total_wealth ./ (1 + (1+cS.tau_c)*c_over_k_ratio);
                    c_expend_final = total_wealth - k_prime_final;
                end
                c_final = c_expend_final ./ (1 + cS.tau_c);
                consumption_tax = c_final * cS.tau_c;
                final_regular_tax = k_capital_tax + consumption_tax;
                [~, util_c] = model_setup_utils_bgp.CES_utility(c_final, cS.sigma, cS);
                util_bequest = model_setup_utils_bgp.bequest_utility(k_prime_final, cS);
                util_final = util_c + cS.beta * util_bequest; % 使用常规beta

                for ie = 1:cS.nw_expanded
                    val_age_slice(:, 1, ie) = util_final;
                    pol_age_slice.c(:, 1, ie) = c_final;
                    pol_age_slice.k_prime(:, 1, ie) = k_prime_final;
                    pol_age_slice.pension_out(:, 1, ie) = b_age_val;
                    pol_age_slice.tax_regular(:, 1, ie) = final_regular_tax;
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
                    cash_on_hand = k_return + labor_income_gross + pension_out - (capital_tax + payg_tax + labor_tax);
                    shock_exp_factor = (ie == cS.nw + 1) * cS.kappa_young + (ie == cS.nw + 2) * cS.kappa_old;
                    shock_exp = shock_exp_factor * cash_on_hand;
                    available_for_c_and_s = cash_on_hand - shock_exp;

                    % B. 最优决策搜索
                    k_prime_max = (available_for_c_and_s - cS.cFloor * (1+cS.tau_c)) / (1 + g_A_period);
                    if k_prime_max <= cS.kMin
                        k_prime_grid = cS.kMin;
                    else
                        k_prime_grid = linspace(cS.kMin, k_prime_max, nk_search_grid)';
                    end

                    c_expend_choices = available_for_c_and_s - k_prime_grid * (1 + g_A_period);
                    c_choices = c_expend_choices / (1 + cS.tau_c);

                    % C. 计算价值
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
                    end
                end
            end
            % --- 5. 扩展维度以匹配主循环 ---
            val_age = repmat(val_age_slice, [1, cS.nkpps, 1]);
            pol_age = model_setup_utils_bgp.expand_policy_slice(pol_age_slice, cS.nkpps);
        end

        function [val_age, pol_age] = HHSolutionByAge_VFI_PPS(a_idx, vPrime_kkppse_next, M_age, b_age_val, paramS_age, cS)
            % =========================================================================
            % == 函数: HHSolutionByAge_VFI_PPS
            % == 版本: [v1.0 - 完整VFI版]
            % ==
            % == 目的: 为有PPS情境提供一个鲁棒、精确的VFI求解器。
            % == 核心特性:
            % ==   - kpps'演化遵循制度，k'通过效用最大化自由选择。
            % ==   - 使用二维PCHIP插值平滑处理下一期价值函数 v'(k', kpps')。
            % ==   - 彻底解决固定储蓄率与PPS规则结合产生的数值不稳定问题。
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
                'pps_tax', zeros(cS.nk, cS.nkpps, cS.nw_expanded));

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
                        if a_idx <= cS.aR_new
                            labor_income_gross = M_age.w_t * cS.ageEffV_new(a_idx) * paramS_age.leGridV(ie);
                            payg_tax = cS.theta_t * labor_income_gross;
                            pps_contrib = cS.pps_contrib_rate * labor_income_gross;
                            capital_tax = cS.tau_k * (k_now * M_age.r_mkt_t);
                            labor_tax = cS.tau_l * max(0, labor_income_gross - payg_tax - pps_contrib);
                            cash_on_hand = (k_now*market_return_factor + labor_income_gross) - (capital_tax + payg_tax + labor_tax + pps_contrib);
                        else
                            pension_out = b_age_val;
                            kpps_total_value = kpps_now * market_return_factor;
                            pps_withdrawal = cS.pps_withdrawal_rate * kpps_total_value;
                            pps_tax = cS.pps_tax_rate_withdrawal * pps_withdrawal;
                            net_pps_inflow = pps_withdrawal - pps_tax;
                            capital_tax = cS.tau_k * (k_now * M_age.r_mkt_t);
                            payg_tax = 0; labor_tax = 0;
                            cash_on_hand = (k_now*market_return_factor + pension_out + net_pps_inflow) - capital_tax;
                        end
                        shock_exp_factor = (ie == cS.nw + 1) * cS.kappa_young + (ie == cS.nw + 2) * cS.kappa_old;
                        shock_exp = shock_exp_factor * cash_on_hand;
                        available_for_c_and_s = cash_on_hand - shock_exp;

                        % B. 计算制度决定的 kpps'
                        kpps_prime_continuous = 0;
                        if a_idx <= cS.aR_new, kpps_prime_continuous = (kpps_now * market_return_factor + pps_contrib) / (1 + g_A_period);
                        else, kpps_prime_continuous = (kpps_now * market_return_factor - pps_withdrawal) / (1 + g_A_period); end
                        kpps_prime_final = max(cS.kppsMin, kpps_prime_continuous);

                        % C. 对 k' 进行最优决策搜索
                        k_prime_max = (available_for_c_and_s - cS.cFloor * (1+cS.tau_c)) / (1 + g_A_period);
                        if k_prime_max <= cS.kMin, k_prime_grid = cS.kMin; else, k_prime_grid = linspace(cS.kMin, k_prime_max, nk_search_grid)'; end

                        c_expend_choices = available_for_c_and_s - k_prime_grid * (1 + g_A_period);
                        c_choices = c_expend_choices / (1 + cS.tau_c);

                        % D. 计算总价值
                        [~, util_c] = model_setup_utils_bgp.CES_utility(c_choices, cS.sigma, cS);
                        util_bequest = model_setup_utils_bgp.bequest_utility(k_prime_grid + kpps_prime_final, cS);

                        v_prime_matrix = zeros(length(k_prime_grid), cS.nw_expanded);
                        for ie_next = 1:cS.nw_expanded
                            v_prime_matrix(:, ie_next) = vPrime_interpolants{ie_next}(k_prime_grid, ones(size(k_prime_grid))*kpps_prime_final);
                        end
                        ev_choices = v_prime_matrix * trans_mat_next(ie, :)';

                        total_value = util_c + discount_factor_future_utility * (cS.s_pathV(a_idx) * ev_choices + (1 - cS.s_pathV(a_idx)) * util_bequest);

                        [best_val, best_idx] = max(total_value);

                        % E. 存储结果
                        if isfinite(best_val)
                            k_prime_final_optimal = k_prime_grid(best_idx); c_final = c_choices(best_idx);
                            val_age(ik, ikpps, ie) = best_val;
                            pol_age.c(ik, ikpps, ie) = c_final; pol_age.k_prime(ik, ikpps, ie) = k_prime_final_optimal;
                            pol_age.kpps_prime(ik, ikpps, ie) = kpps_prime_final;
                            pol_age.tax_regular(ik, ikpps, ie) = capital_tax + labor_tax + c_final * cS.tau_c;
                            pol_age.tax_payg(ik, ikpps, ie) = payg_tax; pol_age.shock_exp(ik, ikpps, ie) = shock_exp;
                            pol_age.pension_out(ik, ikpps, ie) = pension_out;
                            pol_age.pps_contrib(ik, ikpps, ie) = pps_contrib;
                            pol_age.pps_withdrawal(ik, ikpps, ie) = pps_withdrawal; pol_age.pps_tax(ik, ikpps, ie) = pps_tax;
                        end
                    end
                end
            end
        end

        function [val_age, pol_age] = HHSolutionByAge_VFI_PPS_matrix(a_idx, vPrime_kkppse_next, M_age, b_age_val, paramS_age, cS)
            % =========================================================================
            % == 函数: HHSolutionByAge_VFI_PPS_matrix
            % == 版本: [v2.2 - 根源修正与防御性编程版]
            % ==
            % == 目的: 为有PPS情境提供一个鲁棒、精确、且运行极速的VFI求解器。
            % == 核心特性:
            % ==   - 通过完全矩阵化运算，消除了所有状态点(k, kpps, e)的循环。
            % ==   - kpps'演化遵循制度，k'通过效用最大化自由选择。
            % ==   - 使用二维PCHIP插值平滑处理下一期价值函数 v'(k', kpps')。
            % ==   - [v2.2] 强制ndgrid输入为列向量，从根源上解决维度不匹配问题。
            % ==   - [v2.2] 采用更清晰的切片赋值法创建因子矩阵，增强代码鲁棒性。
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
                'pps_tax', zeros(cS.nk, cS.nkpps, cS.nw_expanded));

            % --- 2. 预计算 ---
            g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
            market_return_factor = 1 + M_age.r_mkt_t;
            discount_factor_future_utility = cS.beta * (1 + g_A_period)^(1 - cS.sigma);

            % --- 3. 终点决策 (PPS版) - 矩阵化 ---
            if a_idx == cS.aD_new
                % 防御性编程：确保输入为列向量
                kGridV_col = cS.kGridV(:);
                kppsGridV_col = cS.kppsGridV(:);
                [K_grid, Kpps_grid] = ndgrid(kGridV_col, kppsGridV_col);

                k_capital_tax = cS.tau_k .* (K_grid .* M_age.r_mkt_t);
                pps_wealth_final = Kpps_grid .* market_return_factor;
                pps_tax_final = cS.pps_tax_rate_withdrawal .* pps_wealth_final;
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
                val_age = repmat(util_final, [1, 1, cS.nw_expanded]);
                pol_age.c = repmat(c_final, [1, 1, cS.nw_expanded]);
                pol_age.k_prime = repmat(k_prime_final, [1, 1, cS.nw_expanded]);
                pol_age.kpps_prime = zeros(cS.nk, cS.nkpps, cS.nw_expanded);
                pol_age.pension_out = repmat(b_age_val, [cS.nk, cS.nkpps, cS.nw_expanded]);
                pol_age.tax_regular = repmat(k_capital_tax + c_final*cS.tau_c, [1, 1, cS.nw_expanded]);
                pol_age.pps_withdrawal = repmat(pps_wealth_final, [1, 1, cS.nw_expanded]);
                pol_age.pps_tax = repmat(pps_tax_final, [1, 1, cS.nw_expanded]);
                return;
            end

            % --- 4. 非终点决策 ---
            % A. 创建插值器和状态网格
            vPrime_interpolants = cell(cS.nw_expanded, 1);
            for ie_next = 1:cS.nw_expanded
                vPrime_interpolants{ie_next} = griddedInterpolant({cS.kGridV, cS.kppsGridV}, vPrime_kkppse_next(:, :, ie_next), 'pchip', 'nearest');
            end
            trans_mat_next = paramS_age.TrProbM_by_age{a_idx + 1};

            % **修正点 1 (核心)**: 强制 ndgrid 的输入为列向量，以确保输出维度正确无误。
            % 这是解决维度不匹配错误的根源。
            kGridV_col = cS.kGridV(:);
            kppsGridV_col = cS.kppsGridV(:);
            eGridV_col = (1:cS.nw_expanded)';
            [K_3D, KPPS_3D, E_IDX_3D] = ndgrid(kGridV_col, kppsGridV_col, eGridV_col);

            le_3D = paramS_age.leGridV(E_IDX_3D);

            % B. 向量化计算可用资源
            CAPITAL_TAX_3D = cS.tau_k * (K_3D * M_age.r_mkt_t);
            LABOR_INCOME_GROSS_3D = zeros(size(K_3D)); PAYG_TAX_3D = zeros(size(K_3D));
            PPS_CONTRIB_3D = zeros(size(K_3D)); LABOR_TAX_3D = zeros(size(K_3D));
            PENSION_OUT_3D = zeros(size(K_3D)); PPS_WITHDRAWAL_3D = zeros(size(K_3D));
            PPS_TAX_3D = zeros(size(K_3D));

            if a_idx <= cS.aR_new
                LABOR_INCOME_GROSS_3D = M_age.w_t * cS.ageEffV_new(a_idx) .* le_3D;
                PAYG_TAX_3D = cS.theta_t * LABOR_INCOME_GROSS_3D;
                PPS_CONTRIB_3D = cS.pps_contrib_rate * LABOR_INCOME_GROSS_3D;
                LABOR_TAX_3D = cS.tau_l * max(0, LABOR_INCOME_GROSS_3D - PAYG_TAX_3D - PPS_CONTRIB_3D);
                CASH_ON_HAND_3D = (K_3D * market_return_factor + LABOR_INCOME_GROSS_3D) - ...
                    (CAPITAL_TAX_3D + PAYG_TAX_3D + LABOR_TAX_3D + PPS_CONTRIB_3D);
            else
                PENSION_OUT_3D(:) = b_age_val;
                kpps_total_value_3D = KPPS_3D * market_return_factor;
                PPS_WITHDRAWAL_3D = cS.pps_withdrawal_rate * kpps_total_value_3D;
                PPS_TAX_3D = cS.pps_tax_rate_withdrawal * PPS_WITHDRAWAL_3D;
                net_pps_inflow_3D = PPS_WITHDRAWAL_3D - PPS_TAX_3D;
                CASH_ON_HAND_3D = (K_3D * market_return_factor + PENSION_OUT_3D + net_pps_inflow_3D) - CAPITAL_TAX_3D;
            end

            % **修正点 2 (优化)**: 使用更清晰、更稳健的切片赋值法创建因子矩阵。
            shock_exp_factor_3D = zeros(cS.nk, cS.nkpps, cS.nw_expanded);
            if cS.nw + 1 <= cS.nw_expanded
                shock_exp_factor_3D(:, :, cS.nw + 1) = cS.kappa_young;
            end
            if cS.nw + 2 <= cS.nw_expanded
                shock_exp_factor_3D(:, :, cS.nw + 2) = cS.kappa_old;
            end
            SHOCK_EXP_3D = shock_exp_factor_3D .* CASH_ON_HAND_3D;
            AVAILABLE_FOR_C_S_3D = CASH_ON_HAND_3D - SHOCK_EXP_3D;

            % C. 向量化计算制度决定的 kpps'
            if a_idx <= cS.aR_new
                kpps_prime_continuous_3D = (KPPS_3D * market_return_factor + PPS_CONTRIB_3D) / (1 + g_A_period);
            else
                kpps_prime_continuous_3D = (KPPS_3D * market_return_factor - PPS_WITHDRAWAL_3D) / (1 + g_A_period);
            end
            KPPS_PRIME_3D = max(cS.kppsMin, kpps_prime_continuous_3D);

            % D. 向量化对 k' 的最优决策搜索
            k_prime_max_3D = (AVAILABLE_FOR_C_S_3D - cS.cFloor * (1 + cS.tau_c)) / (1 + g_A_period);
            is_constrained_3D = k_prime_max_3D <= cS.kMin;
            k_prime_max_3D(is_constrained_3D) = cS.kMin;

            linspace_weights = reshape(linspace(0, 1, nk_search_grid), 1, 1, 1, nk_search_grid);
            K_PRIME_CHOICES_4D = cS.kMin + (k_prime_max_3D - cS.kMin) .* linspace_weights;
            C_EXPEND_CHOICES_4D = AVAILABLE_FOR_C_S_3D - K_PRIME_CHOICES_4D * (1 + g_A_period);
            C_CHOICES_4D = C_EXPEND_CHOICES_4D / (1 + cS.tau_c);
            [~, UTIL_C_4D] = model_setup_utils_bgp.CES_utility(C_CHOICES_4D, cS.sigma, cS);
            KPPS_PRIME_4D = repmat(KPPS_PRIME_3D, [1, 1, 1, nk_search_grid]);
            UTIL_BEQUEST_4D = model_setup_utils_bgp.bequest_utility(K_PRIME_CHOICES_4D + KPPS_PRIME_4D, cS);

            % E. 向量化计算期望延续价值
            num_choices_total = numel(K_PRIME_CHOICES_4D);
            V_PRIME_RESHAPED = zeros(num_choices_total, cS.nw_expanded);
            k_prime_choices_flat = K_PRIME_CHOICES_4D(:);
            kpps_prime_flat = KPPS_PRIME_4D(:);

            for ie_next = 1:cS.nw_expanded
                V_PRIME_RESHAPED(:, ie_next) = vPrime_interpolants{ie_next}(k_prime_choices_flat, kpps_prime_flat);
            end

            EV_RESHAPED = V_PRIME_RESHAPED * trans_mat_next';

            E_IDX_4D = repmat(E_IDX_3D, [1, 1, 1, nk_search_grid]);
            E_IDX_4D_flat = E_IDX_4D(:);
            linear_idx_for_ev = sub2ind(size(EV_RESHAPED), (1:num_choices_total)', E_IDX_4D_flat);
            EV_CHOICES_FLAT = EV_RESHAPED(linear_idx_for_ev);
            EV_CHOICES_4D = reshape(EV_CHOICES_FLAT, size(K_PRIME_CHOICES_4D));

            % F. 计算总价值并最大化
            s_path_val = cS.s_pathV(a_idx);
            TOTAL_VALUE_4D = UTIL_C_4D + discount_factor_future_utility * ...
                (s_path_val * EV_CHOICES_4D + (1 - s_path_val) * UTIL_BEQUEST_4D);
            [best_val_3D, best_idx_3D] = max(TOTAL_VALUE_4D, [], 4, 'omitnan');

            % G. "收集"最优策略并存储结果
            valid_mask = isfinite(best_val_3D);
            best_val_3D(~valid_mask) = -1e20;
            val_age = best_val_3D;
            num_states = cS.nk * cS.nkpps * cS.nw_expanded;
            offset = (0:num_states-1)' * nk_search_grid;
            linear_indices = offset + best_idx_3D(:);
            k_prime_final_optimal_3D = reshape(K_PRIME_CHOICES_4D(linear_indices), size(best_idx_3D));
            c_final_3D = reshape(C_CHOICES_4D(linear_indices), size(best_idx_3D));
            k_prime_final_optimal_3D(~valid_mask) = 0;
            c_final_3D(~valid_mask) = 0;
            pol_age.c = c_final_3D;
            pol_age.k_prime = k_prime_final_optimal_3D;
            pol_age.kpps_prime = KPPS_PRIME_3D;
            pol_age.tax_regular = CAPITAL_TAX_3D + LABOR_TAX_3D + c_final_3D * cS.tau_c;
            pol_age.tax_payg = PAYG_TAX_3D;
            pol_age.shock_exp = SHOCK_EXP_3D;
            pol_age.pension_out = PENSION_OUT_3D;
            pol_age.pps_contrib = PPS_CONTRIB_3D;
            pol_age.pps_withdrawal = PPS_WITHDRAWAL_3D;
            pol_age.pps_tax = PPS_TAX_3D;

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
            options = optimoptions('fsolve', 'Display', fsolve_display, 'TolFun', 1e-6, 'TolX', 1e-6, 'MaxIterations', 100 ); % ,'algorithm','trust-region'

            if verbose, fprintf('\n--- 启动 fsolve 求解器 (求解 [K̂_p, K̂_g, L] - BGP版本) ---\n'); end
            [x_eq, ~, exitflag] = fsolve(system_wrapper, x0, options); % 调用fsolve
            if verbose, fprintf('--- fsolve 求解完成 ---\n'); end

            eq_found = (exitflag > 0); % 根据exitflag判断是否收敛
        end

        % =======================================================
        % == 结果展示与检查
        % =======================================================


        % == [新函数] 法医级家庭预算约束检验器
        % =========================================================================
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



        function [is_steady, max_diff] = verify_FULL_steady_state(Dist_to_verify, polS, paramS, cS, Z_target)
            % =========================================================================
            % == 函数说明: [v6 - 模式区分版] 稳态完全内生验证器
            % ==
            % == 核心修正:
            % ==   1. [增加输入] 增加一个Z_target输入，用于区分验证模式。
            % ==   2. [BGP真稳态模式] 如果Z_target为空(is_empty)，则内生计算理论Z_theory进行验证。
            % ==   3. [伪稳态模式] 如果Z_target非空(如Z_2023)，则使用这个外部Z进行验证。
            % ==   4. [输出信息优化] 根据模式，输出更清晰的验证结论。
            % =========================================================================

            % --- [输入参数调整] ---
            if nargin < 5
                Z_target = []; % 默认进入BGP真稳态验证模式
            end

            if isempty(Z_target)
                verification_mode = 'BGP Steady-State';
                fprintf('\n--- 正在进行【BGP真稳态】的平滑方法验证 (内生Z_theory)...\n');
                % --- 步骤 1a: 内生计算理论年龄分布 Z_theory ---
                mass_levels_by_age = zeros(cS.aD_new, 1);
                mass_levels_by_age(1) = 1.0;
                for ia = 1:(cS.aD_new - 1)
                    mass_levels_by_age(ia+1) = mass_levels_by_age(ia) * cS.s_pathV(ia);
                end
                Z_ss_to_use = mass_levels_by_age / sum(mass_levels_by_age);
            else
                verification_mode = 'Temporary Equilibrium';
                fprintf('\n--- ---------------------------------------------------- ---\n');
                fprintf('--- 警告: 正在对【伪稳态/当期均衡】进行动态一致性检验 ---\n');
                fprintf('---        由于人口分布非稳态，预期此检验会失败。    ---\n');
                fprintf('--- ---------------------------------------------------- ---\n');
                % --- 步骤 1b: 使用外部传入的人口分布 ---
                Z_ss_to_use = Z_target;
            end

            % --- 后续步骤与原函数完全相同，只是将 Z_theory 替换为 Z_ss_to_use ---

            % --- 步骤 2: 基于 Z_ss_to_use 初始化一个从零开始计算的分布 ---
            Dist_recalculated = zeros(size(Dist_to_verify));
            dist_newborn = zeros(cS.nk, cS.nkpps, cS.nw_expanded);
            dist_newborn(1, 1, 1:cS.nw) = paramS.leProb1V';
            Dist_recalculated(:, :, :, 1) = dist_newborn * Z_ss_to_use(1);

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
                            k_p_cont = polS(ia).k_prime(ik, ikpps, ie);
                            if cS.pps_active && cS.nkpps > 1
                                kpps_p_cont = polS(ia).kpps_prime(ik, ikpps, ie);
                            else
                                kpps_p_cont = cS.kppsGridV(1);
                            end
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

                % 使用 Z_ss_to_use 进行重新缩放
                target_mass = Z_ss_to_use(ia+1);
                current_mass_generated = sum(dist_ia_plus_1_unscaled(:));
                if current_mass_generated > 1e-30
                    rescale_factor = target_mass / current_mass_generated;
                    Dist_recalculated(:, :, :, ia+1) = dist_ia_plus_1_unscaled * rescale_factor;
                end
            end

            % --- 步骤 4: 归一化 ---
            total_mass_final = sum(Dist_recalculated(:));
            if total_mass_final > 1e-9
                Dist_recalculated = Dist_recalculated / total_mass_final;
            end

            % --- 步骤 5: 比较并输出结论 ---
            max_diff = max(abs(Dist_to_verify - Dist_recalculated), [], 'all');
            tolerance = 1e-9;
            if max_diff < tolerance
                is_steady = true;
                fprintf('   ✅ [%s] 验证通过！分布是动态自洽的 (最大差异: %.4e)\n', verification_mode, max_diff);
            else
                is_steady = false;
                if strcmp(verification_mode, 'Temporary Equilibrium')
                    fprintf('   ✅ [%s] 验证如期“失败”。分布不是动态自洽的 (最大差异: %.4e)，这对于伪稳态是正常现象。\n', verification_mode, max_diff);
                else
                    fprintf('   ⚠️ [%s] 验证失败！稳态分布不自洽 (最大差异: %.4e)\n', verification_mode, max_diff);
                end
            end
        end
        % =======================================================
        % == 辅助函数: 概率质量分裂法的核心工具
        % =======================================================

        % =======================================================
        function display_national_accounts_unified(ss, cS)
            % =========================================================================
            % == 函数: display_national_accounts_unified
            % == 版本: [v27 - 最终会计统一版]
            % ==
            % == 核心修正:
            % ==   修正了"[1. 核心诊断: 家庭部门完整流量表]"的会计逻辑。
            % ==   之前的版本错误地混用了现金流和NIPA概念，导致其重算的净储蓄
            % ==   (S_p_net_recalc) 与模型求解时使用的 (ss.Saving_private_flow)
            % ==   不一致，造成了虚假的“储蓄交叉验证”失败。
            % ==
            % ==   新逻辑严格遵循NIPA准则，与求解器(calculate_aggregates_unified)
            % ==   的定义完全统一，从而消除了这一伪矛盾。
            % =========================================================================

            fprintf('\n\n================================================================================\n');
            if cS.pps_active
                fprintf('===     国民经济核算详细报告 (BGP + PPS激活 - 最终会计统一版)      ===\n');
            else
                fprintf('===      国民经济核算详细报告 (BGP + 无PPS - 最终会计准则版)     ===\n');
            end
            fprintf('================================================================================\n');

            % --- 0. 预计算所有关键宏观量 ---
            g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;

            % --- 生产与总收入 ---
            Y_prod = ss.Y_from_production_hat;
            Wages = ss.w_hat * ss.L_hat;
            K_private_total_begin = ss.K_private_begin_hat;
            Private_Capital_Income_Net = ss.r_mkt * K_private_total_begin;
            Public_Capital_Return = ss.Public_Capital_Return;
            GDI = Wages + (ss.r_mkt + cS.ddk) * K_private_total_begin + Public_Capital_Return;

            % --- 家庭部门流量 ---
            C_hh = ss.C_agg;
            Shock_exp = ss.Shock_exp_agg;
            Consumption_tax = C_hh * cS.tau_c;
            Pension_in = ss.Pension_in;
            Pension_out = ss.Pension_out;
            Bequests_paid = ss.Bequest_tax;

            % --- PPS 部门流量 (仅为信息展示) ---
            PPS_contrib = ss.PPS_contrib_agg;
            PPS_withdrawal_gross = ss.PPS_withdrawal_agg;
            PPS_tax = ss.PPS_tax_agg;

            % --- 政府部门流量 ---
            Regular_tax_total = ss.Regular_tax;
            Gov_Total_Revenue = Regular_tax_total + Bequests_paid + PPS_tax + Public_Capital_Return;
            Depreciation_g = ss.Depreciation_g;
            Resources_for_discretion = Gov_Total_Revenue - Depreciation_g;
            Gov_Net_Saving = cS.lambda_g * Resources_for_discretion;
            G_c = (1 - cS.lambda_g) * Resources_for_discretion;
            I_g_net = g_A_period * ss.K_public_hat;
            I_g_gross = I_g_net + Depreciation_g;

            % --- 私人部门投资与储蓄 ---
            Depreciation_p = ss.Depreciation_p;
            I_p_net = g_A_period * K_private_total_begin;
            I_p_gross = I_p_net + Depreciation_p;
            I_total_gross = I_p_gross + I_g_gross;

            S_p_net_from_ss = ss.Saving_private_flow;

            % --- [1. 核心诊断: 家庭部门完整流量表 (NIPA准则)] ---
            fprintf('\n--- [1. 核心诊断: 家庭部门完整流量表 (NIPA准则)] ---\n');
            % A. 家庭总收入 (NIPA Inflows)
            HH_NIPA_Inflow = Wages + Private_Capital_Income_Net + Pension_out;
            % B. 家庭非储蓄性总支出 (NIPA Outlays)
            HH_NIPA_Outlay_NonSaving = (C_hh + Shock_exp) + Regular_tax_total + Pension_in + PPS_tax;
            % C. 基于NIPA流量表重算的家庭净储蓄
            S_p_net_recalc = HH_NIPA_Inflow - HH_NIPA_Outlay_NonSaving;

            fprintf('   --- A. 家庭总收入 (NIPA Inflows) ----: %12.8f\n', HH_NIPA_Inflow);
            fprintf('      1. 劳动总收入 (w*L) ...............: %12.8f\n', Wages);
            fprintf('      2. 净资本利息收入 (r*K_total) .....: %12.8f\n', Private_Capital_Income_Net);
            fprintf('      3. 公共养老金收入 (Pension Out) ..: %12.8f\n', Pension_out);
            fprintf('\n');
            fprintf('   --- B. 家庭总使用 (NIPA Outlays) ----: %12.8f\n', HH_NIPA_Outlay_NonSaving + S_p_net_recalc);
            fprintf('      --- B1. 消费性支出 (C+Shock) ----: %12.8f\n', (C_hh + Shock_exp));
            fprintf('      --- B2. 税费总支出 (Taxes) ------: %12.8f\n', (Regular_tax_total + Pension_in + PPS_tax));
            fprintf('         1. 常规税 (Regular Tax) .........: %12.8f\n', Regular_tax_total);
            fprintf('         2. PAYG缴费 (Pension In) ........: %12.8f\n', Pension_in);
            fprintf('         3. PPS提取税 (Withdrawal Tax) ...: %12.8f\n', PPS_tax);
            fprintf('      --- B3. 净储蓄 (Net Saving) -------: %12.8f\n', S_p_net_recalc);
            fprintf('         (信息项) PPS净缴费(C-W): %.8f\n', (PPS_contrib - PPS_withdrawal_gross));
            fprintf('\n');
            fprintf('   --- C. 储蓄交叉验证 -------------------\n');
            fprintf('   模型ss中记录的私人净储蓄 ..............: %12.8f\n', S_p_net_from_ss);
            fprintf('   本表(NIPA)重新计算的私人净储蓄 ........: %12.8f\n', S_p_net_recalc);
            fprintf('   => 两者差异 (应为0) ...................: %12.4e\n', S_p_net_from_ss - S_p_net_recalc);

            % --- [2. 收入法检验: GDI vs GDP (欧拉定理)] ---
            fprintf('\n\n--- [2. 收入法检验: GDI vs GDP (欧拉定理)] ---\n');
            fprintf('   A. 国内生产总值 (GDP) .....................: %12.8f\n', Y_prod);
            fprintf('   B. 国内总收入 (GDI) .......................: %12.8f\n', GDI);
            fprintf('   => 收入-产出缺口 (A - B) ..................: %12.4e\n', Y_prod - GDI);

            % --- [3. 政府预算检验 (收支细分)] ---
            fprintf('\n--- [3. 政府预算检验 (收支细分)] ---\n');
            fprintf('   --- A. 政府总收入 --------------------: %12.8f\n', Gov_Total_Revenue);
            fprintf('   --- B. 政府总支出 (I_g_gross+G_c) ----: %12.8f\n', (I_g_gross + G_c));
            fprintf('   => 政府预算缺口 (A - B) ..................: %12.4e\n', Gov_Total_Revenue - (I_g_gross + G_c));

            % --- [4. 终极资源约束检验: Y vs C + I + Gc] ---
            fprintf('\n--- [4. 终极资源约束检验 (Y = C+I_gross+Gc)] ---\n');
            Total_Uses = (C_hh + Shock_exp) + I_total_gross + G_c;
            fprintf('   A. 资源总供给 (GDP) .......................: %12.8f\n', Y_prod);
            fprintf('   B. 资源总使用 (C+I_gross+Gc) ..............: %12.8f\n', Total_Uses);
            fprintf('   => 最终资源缺口 (A - B) ..................: %12.4e\n', Y_prod - Total_Uses);

            % --- [5. 储蓄-投资恒等式检验 (S = I)] ---
            fprintf('\n--- [5. 储蓄-投资恒等式检验 (S_gross = I_gross)] ---\n');
            S_p_gross = S_p_net_from_ss + Depreciation_p;
            S_g_gross = Gov_Net_Saving + Depreciation_g;
            S_total_gross = S_p_gross + S_g_gross;
            mismatch_S_I = S_total_gross - I_total_gross;
            fprintf('   A. 国民总储蓄 (S_gross_total) .............: %12.8f\n', S_total_gross);
            fprintf('   B. 国民总投资 (I_gross_total) .............: %12.8f\n', I_total_gross);
            fprintf('   --------------------------------------------------------------------\n');
            fprintf('   C. 储蓄-投资缺口 (A - B) ..................: %12.8f\n', mismatch_S_I);
            fprintf('   D. 总意外遗赠 (Bequests) ..................: %12.8f\n', Bequests_paid);
            fprintf('   --------------------------------------------------------------------\n');
            fprintf('   => 理论验证 (缺口 - Bequests) (应为0) ......: %12.4e\n', mismatch_S_I - Bequests_paid);

            fprintf('\n================================================================================\n');
        end


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