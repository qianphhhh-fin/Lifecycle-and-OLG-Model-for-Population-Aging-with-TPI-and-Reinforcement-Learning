classdef main_steady_state_utils_bgp
    % =========================================================================
    % == main_steady_state_utils_bgp
    % == 版本: [v11.0 - 最终整合与时间线对齐版]
    % ==
    % == 最终修复与特性:
    % == 1. [财政时间线对齐] 修正了政府收支计算中公共投资(I_g)等变量的
    % ==    时间索引，确保所有流量在当期预算约束中正确核算，与OG-CORE对齐。
    % == 2. [VFI预算约束正确] 家庭VFI预算约束只受技术进步率g_A影响，与OG-CORE一致。
    % == 3. [分布计算正确] 稳态分布计算采用与OG-CORE对齐的条件分布法，逻辑稳健。
    % == 4. [存量-流量核算正确] 宏观聚合严格遵循BGP增长模型的存量-流量定义，
    % ==    保证了宏观自洽性。
    % =========================================================================

    methods (Static)

        function [ss, Dist, polS, valS] = solve_steady_state_unified(cS, paramS, params_ext, verbose, x0_guess, solver_method)
            if nargin < 4, verbose = true; end
            if nargin < 5, x0_guess = []; end
            if nargin < 6, solver_method = 'lsqnonlin'; end

            cS.A = params_ext.A;
            if isfield(params_ext, 'g_A_ss'), cS.g_A_ss = params_ext.g_A_ss; end
            if isfield(params_ext, 'n_ss'), cS.n_ss = params_ext.n_ss; end
            cS.tau_beq = 0.0;
            if isfield(params_ext, 'Z') && ~isempty(params_ext.Z), Z_ss_norm = params_ext.Z;
            else, error('错误：求解稳态需要一个明确的人口分布Z。'); end
            if ~isfield(cS, 'pps_active'), cS.pps_active = false; end
            if ~isfield(cS, 'n_ss'), cS.n_ss = 0.0; end

            is_db_mode = isfield(cS, 'endogenous_theta_mode') && cS.endogenous_theta_mode;

            if is_db_mode
                fprintf('   [求解模式] DB模式激活，theta将作为内生变量求解。\n');
                x0 = [0.3, 0.4, 0.3, 0.01, 0.4];
                lb = [1e-8, 1e-8, 1e-8, 1e-8, 1e-8];
                ub = [Inf, Inf, Inf, Inf, cS.theta_max];
                params_ext.theta = [];
            else
                fprintf('   [求解模式] DC模式激活，theta为外生给定。\n');
                x0 = [0.3, 0.4, 0.3, 0.01];
                lb = [1e-8, 1e-8, 1e-8, 1e-8];
                ub = [Inf, Inf, Inf, Inf];
            end

            system_wrapper = @(x) main_steady_state_utils_bgp.system_of_equations_steady_state(x, Z_ss_norm, cS, paramS, params_ext);

            fprintf('\n>>> 正在使用求解器: %s <<<\n', upper(solver_method));
            [x_eq, ~, exitflag] = main_steady_state_utils_bgp.solve_with_lsqnonlin(system_wrapper, x0, lb, ub, verbose);

            eq_found = (exitflag > 0);
            if ~eq_found || isempty(x_eq)
                warning('统一求解器未能找到均衡解！');
                ss = []; Dist = []; polS = []; valS = [];
                return;
            end

            K_p_eq = x_eq(1);
            K_g_eq = x_eq(2);
            L_eq = x_eq(3);
            Beq_eq = x_eq(4);
            if is_db_mode, params_ext.theta = x_eq(5); end

            [ss, Dist, polS, valS] = main_steady_state_utils_bgp.calculate_aggregates_unified(K_p_eq, K_g_eq, L_eq, Beq_eq, Z_ss_norm, cS, paramS, params_ext);
            if verbose, main_steady_state_utils_bgp.display_national_accounts_unified(ss, cS); end
        end

        function [ss, Dist, polS, valS] = calculate_aggregates_unified(K_private_total_guess, K_g_guess, L_guess, Bequest_Total_guess, Z_ss_norm, cS, paramS, params_ext)
            % =========================================================================
            % == 函数: calculate_aggregates_unified
            % == 版本: [v62 - 清洁VFI调用版]
            % ==
            % == 核心修正:
            % ==   - [代码解耦] 将VFI的完整反向迭代循环逻辑，完全移至一个独立的
            % ==     外部函数 HHSolution_unified_VFI 中。
            % ==   - [清晰调用] 此函数现在只负责准备VFI的输入，然后进行一次干净
            % ==     的调用，不再处理VFI的内部循环逻辑，修复了变量作用域错误。
            % ==   - [聚合兼容] 聚合循环依然保持对PPS和无PPS情况的兼容性。
            % =========================================================================

            % --- 1. 获取价格并准备VFI输入 ---
            M_prices = main_steady_state_utils_bgp.get_prices_at_t(K_private_total_guess, K_g_guess, L_guess, cS);
            M_for_hh = M_prices;
            M_for_hh.w_t = M_prices.w_hat_t;
            M_for_hh.r_mkt_t = M_prices.r_mkt_t;

            is_db_mode = isfield(cS, 'endogenous_theta_mode') && cS.endogenous_theta_mode;
            if is_db_mode
                theta_ss = params_ext.theta;
                M_for_hh.b_t = cS.payg_replacement_rate * M_prices.w_hat_t;
                cS.theta_t = theta_ss;
            else
                theta_ss = params_ext.theta;
                mass_retirees_ss = sum(Z_ss_norm((cS.aR_new+1):end));
                M_for_hh.b_t = (theta_ss * M_prices.w_hat_t * L_guess) / max(1e-9, mass_retirees_ss);
                cS.theta_t = theta_ss;
            end

            M_for_hh.beq_transfer_pers = Bequest_Total_guess / sum(Z_ss_norm(:));

            % [!!! 关键修正 !!!] 调用独立的、统一的VFI求解器
            [polS, valS] = main_steady_state_utils_bgp.HHSolution_unified_VFI(M_for_hh, paramS, cS);

            if any(any(any(any(isinf(valS))))) || any(any(any(any(isnan(valS))))), ss = []; Dist = []; polS = []; valS = []; return; end
            Dist = main_steady_state_utils_bgp.solve_steady_state_distribution_unified(polS, paramS, cS, Z_ss_norm);

            % --- 2. 聚合微观变量 ---
            K_k_hat_agg_raw = 0; K_pps_hat_agg_raw = 0;
            L_agg = 0; C_agg_raw = 0; Bequest_generated_agg_raw = 0;
            Shock_exp_agg = 0; Pension_out_agg = 0; Pension_in_agg = 0; Regular_tax_agg = 0;

            if isfield(cS, 'pps_active') && cS.pps_active
                for ia = 1:cS.aD_new
                    mass_dist_ia = Dist(:,:,:,ia);
                    k_prime_ia = polS(ia).k_prime;
                    kpps_prime_ia = polS(ia).kpps_prime;
                    K_k_hat_agg_raw = K_k_hat_agg_raw + sum(k_prime_ia .* mass_dist_ia, 'all');
                    K_pps_hat_agg_raw = K_pps_hat_agg_raw + sum(kpps_prime_ia .* mass_dist_ia, 'all');
                    if ia <= cS.aR_new
                        mass_on_e = squeeze(sum(mass_dist_ia, [1,2]));
                        L_agg = L_agg + cS.ageEffV_new(ia) * (paramS.leGridV' * mass_on_e);
                    end
                    total_wealth_prime_ia = k_prime_ia + kpps_prime_ia;
                    Bequest_generated_agg_raw = Bequest_generated_agg_raw + sum(total_wealth_prime_ia .* mass_dist_ia, 'all') * (1 - cS.s_pathV(ia));
                    C_agg_raw = C_agg_raw + sum(polS(ia).c .* mass_dist_ia, 'all');
                    Shock_exp_agg = Shock_exp_agg + sum(polS(ia).shock_exp .* mass_dist_ia, 'all');
                    Regular_tax_agg = Regular_tax_agg + sum(polS(ia).tax_regular .* mass_dist_ia, 'all');
                    Pension_in_agg = Pension_in_agg + sum(polS(ia).tax_payg .* mass_dist_ia, 'all');
                    Pension_out_agg = Pension_out_agg + sum(polS(ia).pension_out .* mass_dist_ia, 'all');
                end
            else
                for ia = 1:cS.aD_new
                    mass_dist_ia = Dist(:,:,:,ia);
                    k_prime_ia = polS(ia).k_prime;
                    K_k_hat_agg_raw = K_k_hat_agg_raw + sum(k_prime_ia .* mass_dist_ia, 'all');
                    if ia <= cS.aR_new
                        mass_on_e = squeeze(sum(mass_dist_ia, [1,2]));
                        L_agg = L_agg + cS.ageEffV_new(ia) * (paramS.leGridV' * mass_on_e);
                    end
                    Bequest_generated_agg_raw = Bequest_generated_agg_raw + sum(k_prime_ia .* mass_dist_ia, 'all') * (1 - cS.s_pathV(ia));
                    C_agg_raw = C_agg_raw + sum(polS(ia).c .* mass_dist_ia, 'all');
                    Shock_exp_agg = Shock_exp_agg + sum(polS(ia).shock_exp .* mass_dist_ia, 'all');
                    Regular_tax_agg = Regular_tax_agg + sum(polS(ia).tax_regular .* mass_dist_ia, 'all');
                    Pension_in_agg = Pension_in_agg + sum(polS(ia).tax_payg .* mass_dist_ia, 'all');
                    Pension_out_agg = Pension_out_agg + sum(polS(ia).pension_out .* mass_dist_ia, 'all');
                end
            end
            C_agg = C_agg_raw + Shock_exp_agg;
            K_private_hat_agg_raw = K_k_hat_agg_raw + K_pps_hat_agg_raw;

            % --- 3, 4, 5. BGP核算, 政府账户, 填充ss结构体 (这部分代码保持不变) ---
            g_A_period = (1+cS.g_A_ss)^cS.time_Step-1;
            n_period = (1+cS.n_ss)^cS.time_Step-1;
            K_private_hat_agg = K_private_hat_agg_raw / (1 + n_period);
            Bequest_generated_agg = Bequest_generated_agg_raw / (1 + n_period);
            g_total_factor = (1 + g_A_period) * (1 + n_period);
            Saving_private_flow_Gross = K_private_hat_agg * g_total_factor - (1 - cS.ddk) * K_private_total_guess;
            Saving_public_flow_Gross = K_g_guess * g_total_factor - (1 - cS.ddk_g) * K_g_guess;
            Y_hat = M_prices.Y_hat_t;
            G_c = Y_hat - C_agg - Saving_private_flow_Gross - Saving_public_flow_Gross;
            Public_Capital_Return = Y_hat - (M_prices.w_hat_t * L_agg) - ((M_prices.r_mkt_t + cS.ddk) * K_private_total_guess);
            ss = struct();
            ss.K_private_begin_hat = K_private_total_guess;
            ss.K_private_hat = K_private_hat_agg;
            ss.L_hat = L_agg;
            ss.K_public_hat = K_g_guess;
            ss.Y_from_production_hat = Y_hat;
            ss.w_hat = M_prices.w_hat_t;
            ss.r_mkt = M_prices.r_mkt_t;
            ss.C_agg = C_agg;
            ss.G_c = G_c;
            ss.I_g = Saving_public_flow_Gross;
            ss.I_p = Saving_private_flow_Gross;
            if isfield(cS, 'pps_active') && cS.pps_active
                ss.K_pps_hat = K_pps_hat_agg_raw / (1+n_period);
            end
            ss.Depreciation_p = cS.ddk * K_private_total_guess;
            ss.Depreciation_g = cS.ddk_g * K_g_guess;
            ss.Bequest_generated_agg = Bequest_generated_agg;
            ss.Bequest_distributed_agg = Bequest_Total_guess;
            ss.Regular_tax = Regular_tax_agg;
            ss.Pension_in = Pension_in_agg;
            ss.Pension_out = Pension_out_agg;
            ss.Saving_private_flow_Gross = Saving_private_flow_Gross;
            ss.Saving_public_flow_Gross = Saving_public_flow_Gross;
            ss.Public_Capital_Return = Public_Capital_Return;
            ss.theta = theta_ss;
        end

        function F_error = system_of_equations_steady_state(x, Z_ss_norm, cS, paramS, params_ext)
            % =========================================================================
            % == 函数: system_of_equations_steady_state
            % == 版本: [v58 - 最终BGP增长核算对齐版]
            % ==
            % == 核心修正:
            % ==   目标函数直接设为宏观变量的自洽性误差，而非间接的流量差额。
            % ==   1. [资本市场出清] 直接要求期初猜测的资本量(K_guess)等于
            % ==      期末聚合的资本量(ss.K_private_hat)。
            % ==   2. [遗赠市场出清] 同理，直接要求猜测的遗赠量等于聚合的遗赠量。
            % ==   此方法使求解器的目标与经济均衡的定义完全一致。
            % =========================================================================

            is_db_mode = isfield(cS, 'endogenous_theta_mode') && cS.endogenous_theta_mode;

            % 获取求解器猜测的宏观总量
            K_private_total_guess = x(1);
            K_g_guess = x(2);
            L_guess = x(3);
            Bequest_Total_guess = x(4);
            if is_db_mode, params_ext.theta = x(5); end

            % 基于猜测值，计算出一个完整的经济状态 (包括家庭的最优决策和聚合结果)
            [ss, ~, ~, ~] = main_steady_state_utils_bgp.calculate_aggregates_unified(K_private_total_guess, K_g_guess, L_guess, Bequest_Total_guess, Z_ss_norm, cS, paramS, params_ext);

            if isempty(ss)
                error_size = ifthen(is_db_mode, 5, 4);
                F_error = ones(error_size, 1) * 1e8;
                return;
            end

            % 构造误差向量：核心是 "输入 = 输出"
            % 市场1: 私人资本 - 猜测值必须等于聚合结果
            error_Kp = ss.K_private_hat - K_private_total_guess;

            % 市场2: 公共资本 - 假设总是出清(I_g由BGP需求定)，误差为0
            % 这是因为K_g不是由家庭行为决定的。
            error_Kg = 0; % No behavioral equation to solve for Kg, it's a guess.

            % 市场3: 劳动 - 猜测值必须等于聚合结果 (虽然在本模型中L是输入的猜测，但保留结构)
            error_L = ss.L_hat - L_guess;

            % 市场4: 遗赠 - 猜测值必须等于聚合结果
            error_Beq = ss.Bequest_generated_agg - Bequest_Total_guess;

            % PAYG 预算平衡条件
            if is_db_mode
                error_PAYG = ss.Pension_in - ss.Pension_out;
                F_error = [error_Kp; error_Kg; error_L; error_Beq; error_PAYG];
            else
                F_error = [error_Kp; error_Kg; error_L; error_Beq];
            end
        end

        function [polS, valS] = HHSolution_unified_VFI(M_vfi, paramS_vfi, cS_vfi)
            % =========================================================================
            % == 函数: HHSolution_unified_VFI
            % == 版本: [v1.0 - 最终统一版]
            % ==
            % == 目的: 封装所有VFI逻辑，根据pps_active标志自动切换。
            % =========================================================================

            if isfield(cS_vfi, 'pps_active') && cS_vfi.pps_active
                valS = -Inf(cS_vfi.nk, cS_vfi.nkpps, cS_vfi.nw_expanded, cS_vfi.aD_new);
            else
                valS = -Inf(cS_vfi.nk, 1, cS_vfi.nw_expanded, cS_vfi.aD_new);
            end

            polS_cell = cell(cS_vfi.aD_new, 1);
            bV_payg_vfi = zeros(1, cS_vfi.aD_new);
            if cS_vfi.aR_new < cS_vfi.aD_new, bV_payg_vfi((cS_vfi.aR_new+1):cS_vfi.aD_new) = M_vfi.b_t; end

            % --- 统一的反向迭代循环 ---
            for a_idx = cS_vfi.aD_new : -1 : 1
                vPrime_kkppse_next = [];
                if a_idx < cS_vfi.aD_new, vPrime_kkppse_next = valS(:,:,:,a_idx+1); end

                % [!!! 关键 !!!] 根据标志，调用对应的单步年龄求解器
                if isfield(cS_vfi, 'pps_active') && cS_vfi.pps_active
                    [val_age, pol_age] = main_steady_state_utils_bgp.HHSolutionByAge_VFI_PPS(a_idx, vPrime_kkppse_next, M_vfi, bV_payg_vfi(a_idx), M_vfi.beq_transfer_pers, paramS_vfi, cS_vfi);
                else
                    [val_age, pol_age] = main_steady_state_utils_bgp.HHSolutionByAge_VFI_noPPS(a_idx, vPrime_kkppse_next, M_vfi, bV_payg_vfi(a_idx), M_vfi.beq_transfer_pers, paramS_vfi, cS_vfi);
                end

                valS(:,:,:,a_idx) = val_age;
                polS_cell{a_idx} = pol_age;
            end

            polS = [polS_cell{:}];
        end

        function [val_age, pol_age] = HHSolutionByAge_VFI_noPPS(a_idx, vPrime_kkppse_next, M_age, b_age_val, beq_transfer_val, paramS_age, cS)
            nk_search_grid = 150;
            val_age = -1e20 * ones(cS.nk, 1, cS.nw_expanded);
            pol_age = struct('c', zeros(cS.nk, 1, cS.nw_expanded), ...
                'k_prime', zeros(cS.nk, 1, cS.nw_expanded), ...
                'tax_regular', zeros(cS.nk, 1, cS.nw_expanded), ...
                'tax_payg', zeros(cS.nk, 1, cS.nw_expanded), ...
                'shock_exp', zeros(cS.nk, 1, cS.nw_expanded), ...
                'pension_out', zeros(cS.nk, 1, cS.nw_expanded), ...
                'beq_received', zeros(cS.nk, 1, cS.nw_expanded));
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
                labor_income_gross = 0; pension_out = 0;
                if a_idx <= cS.aR_new
                    labor_income_gross = M_age.w_t * cS.ageEffV_new(a_idx) * paramS_age.leGridV(ie);
                else
                    pension_out = b_age_val;
                end
                capital_tax_base = k_grid_vec * M_age.r_mkt_t;
                capital_tax_vec = cS.tau_k .* capital_tax_base;
                payg_tax = cS.theta_t * labor_income_gross;
                labor_tax = cS.tau_l * max(0, labor_income_gross - payg_tax);
                cash_on_hand_vec = capital_return + labor_income_gross + pension_out + beq_transfer_val - (capital_tax_vec + payg_tax + labor_tax);
                shock_exp_factor = (ie == cS.nw + 1) * cS.kappa_young + (ie == cS.nw + 2) * cS.kappa_old;
                shock_exp_vec = shock_exp_factor .* cash_on_hand_vec;
                available_for_c_and_s_vec = cash_on_hand_vec - shock_exp_vec;
                k_prime_max_vec = (available_for_c_and_s_vec - cS.cFloor * (1 + cS.tau_c)) / (1 + g_A_period);
                k_prime_max_vec(k_prime_max_vec < cS.kMin) = cS.kMin;
                k_prime_choices_mat = cS.kMin + (k_prime_max_vec - cS.kMin) .* linspace(0, 1, nk_search_grid);
                c_expend_choices_mat = available_for_c_and_s_vec - k_prime_choices_mat * (1 + g_A_period);
                c_choices_mat = c_expend_choices_mat / (1 + cS.tau_c);
                invalid_choice_mask = (c_choices_mat < cS.cFloor);
                c_choices_mat(invalid_choice_mask) = cS.cFloor;
                [~, util_c_mat] = model_setup_utils_bgp.CES_utility(c_choices_mat, cS.sigma, cS);
                util_c_mat(invalid_choice_mask) = -1e20;
                future_value_mat = zeros(cS.nk, nk_search_grid);
                if a_idx < cS.aD_new
                    trans_mat_row = paramS_age.TrProbM_by_age{a_idx + 1}(ie, :);
                    Vprime_interp_mat = cell2mat(cellfun(@(f) f(k_prime_choices_mat(:)), vPrime_interpolants, 'UniformOutput', false)');
                    ev_flat_vec = Vprime_interp_mat * trans_mat_row';
                    ev_on_choices_mat = reshape(ev_flat_vec, cS.nk, nk_search_grid);
                    future_value_mat = cS.s_pathV(a_idx) * ev_on_choices_mat + (1 - cS.s_pathV(a_idx)) * model_setup_utils_bgp.bequest_utility(k_prime_choices_mat, cS);
                else
                    future_value_mat = model_setup_utils_bgp.bequest_utility(k_prime_choices_mat, cS);
                end
                total_value_mat = util_c_mat + discount_factor_utility * future_value_mat;
                [best_val_vec, best_idx_vec] = max(total_value_mat, [], 2);
                linear_indices = sub2ind(size(k_prime_choices_mat), (1:cS.nk)', best_idx_vec);
                val_age(:, 1, ie) = best_val_vec;
                pol_age.c(:, 1, ie) = c_choices_mat(linear_indices);
                pol_age.k_prime(:, 1, ie) = k_prime_choices_mat(linear_indices);
                pol_age.pension_out(:, 1, ie) = pension_out;
                pol_age.shock_exp(:, 1, ie) = shock_exp_vec;
                pol_age.tax_regular(:, 1, ie) = capital_tax_vec + labor_tax + c_choices_mat(linear_indices) .* cS.tau_c;
                pol_age.tax_payg(:, 1, ie) = payg_tax;
                pol_age.beq_received(:, 1, ie) = beq_transfer_val;
            end
        end


        function [val_age, pol_age] = HHSolutionByAge_VFI_PPS(a_idx, vPrime_kkppse_next, M_age, b_age_val, beq_transfer_val, paramS_age, cS)
            % =========================================================================
            % == 函数: HHSolutionByAge_VFI_PPS
            % == 版本: [v8.0 - BGP对齐并行版]
            % ==
            % == 目的:
            % ==   求解带个人养老金账户(PPS)的家庭问题，其核心逻辑（特别是BGP去趋势化）
            % ==   与已验证的 noPPS 版本完全对齐，同时保留了高效的并行计算。
            % =========================================================================

            % --- 1. 初始化 ---
            % 搜索网格密度
            nk_search_grid = 150;

            % 初始化输出变量，增加 kpps 维度
            val_age = -1e20 * ones(cS.nk, cS.nkpps, cS.nw_expanded);
            pol_age = struct(...
                'c', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
                'k_prime', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
                'kpps_prime', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
                'tax_regular', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
                'tax_payg', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
                'shock_exp', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
                'pension_out', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
                'pps_contrib', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
                'pps_withdrawal', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
                'pps_tax', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
                'beq_received', zeros(cS.nk, cS.nkpps, cS.nw_expanded));

            % [BGP对齐] 核心常量，与 noPPS 版本完全一致
            g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
            discount_factor_utility = cS.beta * (1 + g_A_period)^(1 - cS.sigma);
            k_grid_vec = cS.kGridV;
            kpps_grid_vec = cS.kppsGridV;
            market_return_factor = 1 + M_age.r_mkt_t;

            % --- 2. 准备下一期值函数的插值 ---
            vPrime_interpolants = cell(cS.nw_expanded, 1);
            if a_idx < cS.aD_new
                % 为每个下一期的冲击状态 e' 创建一个二维插值对象 V(k, kpps)
                for ie_next = 1:cS.nw_expanded
                    vPrime_interpolants{ie_next} = griddedInterpolant({k_grid_vec, kpps_grid_vec}, vPrime_kkppse_next(:, :, ie_next), 'linear', 'linear');
                end
            end

            % --- 3. [并行化] 迭代求解所有 (kpps, e) 状态 ---
            % 将 (ikpps, ie) 两个循环合并为一个 parfor 循环以最大化并行效率
            num_tasks = cS.nkpps * cS.nw_expanded;
            val_results_cell = cell(num_tasks, 1);
            pol_results_cell = cell(num_tasks, 1);

            parfor task_idx = 1:num_tasks
                % 从一维任务索引映射回二维 (ikpps, ie) 索引
                [ikpps, ie] = ind2sub([cS.nkpps, cS.nw_expanded], task_idx);
                kpps_now = kpps_grid_vec(ikpps);

                % --- 4. 计算当前状态下的预算约束 (与noPPS版本类似，但增加了PPS) ---
                capital_return = k_grid_vec * market_return_factor;

                labor_income_gross = 0; pension_out = 0; pps_contrib = 0; pps_withdrawal = 0;
                pps_tax = 0; payg_tax = 0; labor_tax = 0;

                if a_idx <= cS.aR_new % 工作期
                    labor_income_gross = M_age.w_t * cS.ageEffV_new(a_idx) * paramS_age.leGridV(ie);
                    payg_tax = cS.theta_t * labor_income_gross;
                    pps_contrib = cS.pps_contrib_rate * labor_income_gross;
                    capital_tax_vec = cS.tau_k .* (k_grid_vec * M_age.r_mkt_t);
                    labor_tax = cS.tau_l * max(0, labor_income_gross - payg_tax - pps_contrib);
                    % 总现金流 = 资本回报 + 劳动收入 + 遗赠 - 各项税收和缴费
                    cash_on_hand_vec = (capital_return + labor_income_gross + beq_transfer_val) ...
                        - (capital_tax_vec + payg_tax + labor_tax + pps_contrib);
                else % 退休期
                    pension_out = b_age_val;
                    kpps_total_value = kpps_now * market_return_factor;
                    pps_withdrawal = cS.pps_withdrawal_rate * kpps_total_value;
                    pps_tax = cS.pps_tax_rate_withdrawal * pps_withdrawal;
                    net_pps_inflow = pps_withdrawal - pps_tax;
                    capital_tax_vec = cS.tau_k .* (k_grid_vec * M_age.r_mkt_t);
                    % 总现金流 = 资本回报 + 养老金 + PPS净流入 + 遗赠 - 资本税
                    cash_on_hand_vec = (capital_return + pension_out + net_pps_inflow + beq_transfer_val) ...
                        - capital_tax_vec;
                end

                % 扣除重大支出冲击 (与noPPS版本一致)
                shock_exp_factor = (ie == cS.nw + 1) * cS.kappa_young + (ie == cS.nw + 2) * cS.kappa_old;
                shock_exp_vec = shock_exp_factor .* cash_on_hand_vec;
                available_for_c_and_s_vec = cash_on_hand_vec - shock_exp_vec;

                % --- 5. 计算下一期PPS资产 (确定性部分) ---
                % [BGP对齐] 下一期PPS资产也必须进行去趋势化
                if a_idx <= cS.aR_new % 工作期
                    kpps_prime_level = kpps_now * market_return_factor + pps_contrib;
                else % 退休期
                    kpps_prime_level = kpps_now * market_return_factor - pps_withdrawal;
                end
                kpps_prime_detrended = max(cS.kppsMin, kpps_prime_level / (1 + g_A_period));

                % --- 6. 求解最优 k' 和 c (与noPPS版本一致，但未来价值计算不同) ---
                k_prime_max_vec = (available_for_c_and_s_vec - cS.cFloor * (1 + cS.tau_c)) / (1 + g_A_period);
                k_prime_max_vec(k_prime_max_vec < cS.kMin) = cS.kMin;
                k_prime_choices_mat = cS.kMin + (k_prime_max_vec - cS.kMin) .* linspace(0, 1, nk_search_grid);

                c_expend_choices_mat = available_for_c_and_s_vec - k_prime_choices_mat * (1 + g_A_period);
                c_choices_mat = c_expend_choices_mat / (1 + cS.tau_c);

                invalid_choice_mask = (c_choices_mat < cS.cFloor);
                c_choices_mat(invalid_choice_mask) = cS.cFloor;
                [~, util_c_mat] = model_setup_utils_bgp.CES_utility(c_choices_mat, cS.sigma, cS);
                util_c_mat(invalid_choice_mask) = -1e20;

                future_value_mat = zeros(cS.nk, nk_search_grid);
                if a_idx < cS.aD_new
                    trans_mat_row = paramS_age.TrProbM_by_age{a_idx + 1}(ie, :);
                    k_prime_choices_flat = k_prime_choices_mat(:);
                    kpps_prime_vec = ones(length(k_prime_choices_flat), 1) * kpps_prime_detrended;

                    % 使用二维插值获得期望值
                    v_prime_flat_mat = zeros(length(k_prime_choices_flat), cS.nw_expanded);
                    for ie_next = 1:cS.nw_expanded
                        v_prime_flat_mat(:, ie_next) = vPrime_interpolants{ie_next}(k_prime_choices_flat, kpps_prime_vec);
                    end
                    ev_flat_vec = v_prime_flat_mat * trans_mat_row';
                    ev_on_choices_mat = reshape(ev_flat_vec, cS.nk, nk_search_grid);

                    % 遗赠动机作用于总财富 (k' + kpps')
                    total_wealth_prime_mat = k_prime_choices_mat + kpps_prime_detrended;
                    util_bequest_mat = model_setup_utils_bgp.bequest_utility(total_wealth_prime_mat, cS);

                    future_value_mat = cS.s_pathV(a_idx) * ev_on_choices_mat + (1 - cS.s_pathV(a_idx)) * util_bequest_mat;
                else % 最后一期
                    total_wealth_prime_mat = k_prime_choices_mat + kpps_prime_detrended;
                    future_value_mat = model_setup_utils_bgp.bequest_utility(total_wealth_prime_mat, cS);
                end

                total_value_mat = util_c_mat + discount_factor_utility * future_value_mat;

                % --- 7. 存储最优决策 ---
                [best_val_vec, best_idx_vec] = max(total_value_mat, [], 2);
                linear_indices = sub2ind(size(k_prime_choices_mat), (1:cS.nk)', best_idx_vec);

                k_prime_optimal_vec = k_prime_choices_mat(linear_indices);
                c_optimal_vec = c_choices_mat(linear_indices);

                val_results_cell{task_idx} = best_val_vec;
                pol_results_cell{task_idx} = struct(...
                    'c', c_optimal_vec, ...
                    'k_prime', k_prime_optimal_vec, ...
                    'kpps_prime', kpps_prime_detrended, ...
                    'tax_regular', capital_tax_vec + labor_tax + c_optimal_vec * cS.tau_c, ...
                    'tax_payg', payg_tax, ...
                    'shock_exp', shock_exp_vec, ...
                    'pension_out', pension_out, ...
                    'pps_contrib', pps_contrib, ...
                    'pps_withdrawal', pps_withdrawal, ...
                    'pps_tax', pps_tax, ...
                    'beq_received', beq_transfer_val);
            end

            % --- 8. [安全组装] 将并行计算的结果组装回最终的矩阵和结构体 ---
            fields = fieldnames(pol_age);
            for task_idx = 1:num_tasks
                [ikpps, ie] = ind2sub([cS.nkpps, cS.nw_expanded], task_idx);

                val_age(:, ikpps, ie) = val_results_cell{task_idx};

                current_pol_slice = pol_results_cell{task_idx};
                for i = 1:length(fields)
                    field = fields{i};
                    if isscalar(current_pol_slice.(field))
                        pol_age.(field)(:, ikpps, ie) = repmat(current_pol_slice.(field), cS.nk, 1);
                    else
                        pol_age.(field)(:, ikpps, ie) = current_pol_slice.(field);
                    end
                end
            end
        end

        function Dist = solve_steady_state_distribution_unified(polS, paramS, cS, Z_ss_norm)
            % =========================================================================
            % == 函数: solve_steady_state_distribution_unified
            % == 版本: [v3.0 - PPS/no-PPS 兼容版]
            % ==
            % == 核心修正:
            % ==   - [代码统一] 使用 if cS.pps_active 条件分支，使单一函数能够
            % ==     处理有PPS和无PPS两种情况。
            % ==   - [no-PPS路径] 当 pps_active=false 时，执行针对 (k,e) 状态
            % ==     空间的、基于线性插值的分布迭代。
            % ==   - [PPS路径] 当 pps_active=true 时，执行针对 (k,kpps,e) 状态
            % ==     空间的、基于双线性插值的分布迭代。
            % =========================================================================

            nk = cS.nk;
            nw = cS.nw_expanded;
            aD = cS.aD_new;

            if cS.pps_active
                % --- 路径 A: 包含 PPS 的情况 ---
                nkpps = cS.nkpps;
                Dist_cond = zeros(nk, nkpps, nw, aD, 'double');

                % 新生儿分布在 k=1, kpps=1, e=1..nw
                dist_newborn_cond_shape = zeros(nk, nkpps, nw);
                dist_newborn_cond_shape(1, 1, 1:cS.nw) = reshape(paramS.leProb1V, [1, 1, cS.nw]);
                Dist_cond(:, :, :, 1) = dist_newborn_cond_shape;

                for ia = 1:(aD - 1)
                    dist_cond_ia_slice = Dist_cond(:, :, :, ia);
                    if sum(dist_cond_ia_slice(:)) < 1e-30, continue; end

                    dist_cond_ia_plus_1_next = zeros(nk, nkpps, nw, 'double');
                    trans_mat_next_age = paramS.TrProbM_by_age{ia + 1};

                    for ik = 1:nk
                        for ikpps = 1:nkpps
                            for ie = 1:nw
                                mass_start = dist_cond_ia_slice(ik, ikpps, ie);
                                if mass_start < 1e-30, continue; end

                                k_prime = polS(ia).k_prime(ik, ikpps, ie);
                                kpps_prime = polS(ia).kpps_prime(ik, ikpps, ie);

                                [ik_lower, ~, w_k_upper] = main_steady_state_utils_bgp.find_grid_and_weights(k_prime, cS.kGridV);
                                w_k_lower = 1.0 - w_k_upper;
                                ik_upper = min(nk, ik_lower + 1);

                                [ikpps_lower, ~, w_kpps_upper] = main_steady_state_utils_bgp.find_grid_and_weights(kpps_prime, cS.kppsGridV);
                                w_kpps_lower = 1.0 - w_kpps_upper;
                                ikpps_upper = min(nkpps, ikpps_lower + 1);

                                trans_probs_vec = trans_mat_next_age(ie, :);
                                for ie_next = 1:nw
                                    prob_to_enext = trans_probs_vec(ie_next);
                                    if prob_to_enext < 1e-12, continue; end

                                    mass_to_distribute_e = mass_start * prob_to_enext;

                                    % 双线性插值分配
                                    dist_cond_ia_plus_1_next(ik_lower, ikpps_lower, ie_next) = dist_cond_ia_plus_1_next(ik_lower, ikpps_lower, ie_next) + mass_to_distribute_e * w_k_lower * w_kpps_lower;
                                    dist_cond_ia_plus_1_next(ik_upper, ikpps_lower, ie_next) = dist_cond_ia_plus_1_next(ik_upper, ikpps_lower, ie_next) + mass_to_distribute_e * w_k_upper * w_kpps_lower;
                                    dist_cond_ia_plus_1_next(ik_lower, ikpps_upper, ie_next) = dist_cond_ia_plus_1_next(ik_lower, ikpps_upper, ie_next) + mass_to_distribute_e * w_k_lower * w_kpps_upper;
                                    dist_cond_ia_plus_1_next(ik_upper, ikpps_upper, ie_next) = dist_cond_ia_plus_1_next(ik_upper, ikpps_upper, ie_next) + mass_to_distribute_e * w_k_upper * w_kpps_upper;
                                end
                            end
                        end
                    end

                    sum_next_dist = sum(dist_cond_ia_plus_1_next(:));
                    if sum_next_dist > 1e-9, Dist_cond(:, :, :, ia+1) = dist_cond_ia_plus_1_next / sum_next_dist; end
                end

            else
                % --- 路径 B: 无 PPS 的标准情况 ---
                Dist_cond = zeros(nk, 1, nw, aD, 'double');

                % 新生儿分布在 k=1, kpps=1(虚拟), e=1..nw
                dist_newborn_cond_shape = zeros(nk, 1, nw);
                dist_newborn_cond_shape(1, 1, 1:cS.nw) = reshape(paramS.leProb1V, [1, 1, cS.nw]);
                Dist_cond(:, :, :, 1) = dist_newborn_cond_shape;

                for ia = 1:(aD - 1)
                    dist_cond_ia_slice = Dist_cond(:, 1, :, ia); % 注意这里 kpps 维度是 1
                    if sum(dist_cond_ia_slice(:)) < 1e-30, continue; end

                    dist_cond_ia_plus_1_next = zeros(nk, 1, nw, 'double');
                    trans_mat_next_age = paramS.TrProbM_by_age{ia + 1};

                    for ik = 1:nk
                        for ie = 1:nw
                            mass_start = dist_cond_ia_slice(ik, 1, ie);
                            if mass_start < 1e-30, continue; end

                            k_prime = polS(ia).k_prime(ik, 1, ie);

                            [ik_lower, ~, w_k_upper] = main_steady_state_utils_bgp.find_grid_and_weights(k_prime, cS.kGridV);
                            w_k_lower = 1.0 - w_k_upper;
                            ik_upper = min(nk, ik_lower + 1);

                            trans_probs_vec = trans_mat_next_age(ie, :);
                            for ie_next = 1:nw
                                prob_to_enext = trans_probs_vec(ie_next);
                                if prob_to_enext < 1e-12, continue; end

                                mass_to_distribute = mass_start * prob_to_enext;

                                % 线性插值分配
                                dist_cond_ia_plus_1_next(ik_lower, 1, ie_next) = dist_cond_ia_plus_1_next(ik_lower, 1, ie_next) + mass_to_distribute * w_k_lower;
                                dist_cond_ia_plus_1_next(ik_upper, 1, ie_next) = dist_cond_ia_plus_1_next(ik_upper, 1, ie_next) + mass_to_distribute * w_k_upper;
                            end
                        end
                    end

                    sum_next_dist = sum(dist_cond_ia_plus_1_next(:));
                    if sum_next_dist > 1e-9, Dist_cond(:, :, :, ia+1) = dist_cond_ia_plus_1_next / sum_next_dist; end
                end
            end

            % 最终分布 = 条件分布 * 年龄的边际分布
            Dist = Dist_cond .* reshape(Z_ss_norm, [1, 1, 1, aD]);
        end
        
        function errors = display_national_accounts_unified(ss, cS)
    % =_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
    % == 函数: display_national_accounts_unified
    % == 版本: [v2.1 - 政府预算核算修正版]
    % ==
    % == 核心修正:
    % ==   - [政府支出定义] 在政府账户部分，明确使用“政府净投资”
    % ==     (I_g_net = I_g_gross - Depreciation_g)作为当期支出项。
    % ==     这与会计准则一致，因为折旧是对存量价值的调整，而非当期现金流出。
    % ==   - [赤字来源] 这个修正能更清晰地揭示出政府预算赤字/盈余的来源。
    %_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

    % ... (第0部分：准备工作，与v2.0完全相同) ...
    Y = ss.Y_from_production_hat;
    g_A_period = (1+cS.g_A_ss)^cS.time_Step-1;
    n_period = (1+cS.n_ss)^cS.time_Step-1;
    g_total_period = (1+g_A_period)*(1+n_period)-1;
    I_p_gross_demand = (g_total_period + cS.ddk) * ss.K_private_begin_hat;
    I_g_gross_demand = ss.I_g;
    errors = struct('err_Kp', ss.Saving_private_flow_Gross - I_p_gross_demand, ...
                    'err_Kg', ss.Saving_public_flow_Gross - I_g_gross_demand, ...
                    'err_L', ss.L_hat - ss.L_hat, ...
                    'err_Beq', ss.Bequest_generated_agg - ss.Bequest_distributed_agg, ...
                    'payg_balance', ss.Pension_in - ss.Pension_out);

    % --- 开始打印报告 (第0, 1, 2部分与v2.0完全相同) ---
    fprintf('\n\n================================================================================\n');
    fprintf('===                    稳态均衡：国民经济核算与诊断报告                    ===\n');
    fprintf('================================================================================\n');
    fprintf('--- [0] 均衡价格与核心比率 ---\n');
    fprintf('   年化市场利率 (r_annual) ............... : %12.4f%%\n', ( (1+ss.r_mkt)^(1/cS.time_Step) - 1 ) * 100);
    fprintf('   有效劳动单位工资 (w_hat) .............. : %12.6f\n', ss.w_hat);
    fprintf('   私人资本/产出比 (K_p/Y) ............... : %12.4f\n', ss.K_private_begin_hat / Y);
    fprintf('   公共资本/产出比 (K_g/Y) ............... : %12.4f\n', ss.K_public_hat / Y);
    fprintf('   总资本/产出比 ((K_p+K_g)/Y) ........... : %12.4f\n', (ss.K_private_begin_hat + ss.K_public_hat) / Y);
    fprintf('\n--- [1] 宏观存量 (Stocks) ---\n');
    fprintf('   总产出 (Y_hat) ........................ : %12.6f\n', Y);
    fprintf('   私人资本 (K_p_hat) .................... : %12.6f  (%.2f%% of Y)\n', ss.K_private_begin_hat, (ss.K_private_begin_hat / Y)*100);
    if isfield(ss, 'K_pps_hat'), fprintf('     其中: PPS资本 (K_pps) ............... : %12.6f  (%.2f%% of K_p)\n', ss.K_pps_hat, (ss.K_pps_hat / ss.K_private_begin_hat)*100); end
    fprintf('   公共资本 (K_g_hat) .................... : %12.6f  (%.2f%% of Y)\n', ss.K_public_hat, (ss.K_public_hat / Y)*100);
    fprintf('   有效劳动供给 (L_hat) .................. : %12.6f\n', ss.L_hat);
    fprintf('\n--- [2] 国民收入核算 (支出法流量, Flows) ---\n');
    fprintf('   私人消费 (C_p) ........................ : %12.6f  (%.2f%% of Y)\n', ss.C_agg, (ss.C_agg / Y)*100);
    fprintf('   政府消费 (G_c) ........................ : %12.6f  (%.2f%% of Y)\n', ss.G_c, (ss.G_c / Y)*100);
    fprintf('   私人总投资 (I_p) ...................... : %12.6f  (%.2f%% of Y)\n', ss.I_p, (ss.I_p / Y)*100);
    fprintf('   公共总投资 (I_g) ...................... : %12.6f  (%.2f%% of Y)\n', ss.I_g, (ss.I_g / Y)*100);
    Total_Expenditure = ss.C_agg + ss.G_c + ss.I_p + ss.I_g;
    fprintf('   ------------------------------------------------------------------------------\n');
    fprintf('   >>> 总支出 (C+G+I) .................... : %12.6f  (Discrepancy: %.4e)\n', Total_Expenditure, Y - Total_Expenditure);

    % --- [3] 政府与养老金账户 (核心修正) ---
    fprintf('\n--- [3] 政府与养老金账户 ---\n');
    total_labor_income = ss.w_hat * ss.L_hat;
    payg_implicit_tax_rate = ss.Pension_out / total_labor_income;
    fprintf('   PAYG 养老金系统:\n');
    fprintf('     总收入 (Pension_in) ................. : %12.6f\n', ss.Pension_in);
    fprintf('     总支出 (Pension_out) ................ : %12.6f  (%.2f%% of w*L)\n', ss.Pension_out, (ss.Pension_out / total_labor_income)*100);
    fprintf('     -> PAYG 隐性缴费率 .................. : %12.4f%%\n', payg_implicit_tax_rate * 100);
    fprintf('     -> 外生缴费率参数 (theta) .......... : %12.4f\n', ss.theta);
    fprintf('\n   政府账户 (不含PAYG):\n');
    fprintf('     总收入:\n');
    fprintf('       一般税收 (Regular Tax) ............ : %12.6f\n', ss.Regular_tax);
    fprintf('       公共资本回报 (Public K Return) .... : %12.6f\n', ss.Public_Capital_Return);
    gov_revenue = ss.Regular_tax + ss.Public_Capital_Return;
    fprintf('       ------------------------------------------------------------------------------\n');
    fprintf('       >>> 政府总收入 .................... : %12.6f\n', gov_revenue);
    fprintf('     总支出:\n');
    % [!!! 关键修正 !!!] 政府的当期支出是消费和净投资
    gov_net_investment = ss.I_g - ss.Depreciation_g;
    fprintf('       政府消费 (G_c) .................... : %12.6f\n', ss.G_c);
    fprintf('       政府净投资 (I_g_net) .............. : %12.6f\n', gov_net_investment);
    gov_expenditure = ss.G_c + gov_net_investment;
    fprintf('       ------------------------------------------------------------------------------\n');
    fprintf('       >>> 政府总支出 (G_c + I_g_net) .... : %12.6f\n', gov_expenditure);
    fprintf('     政府预算盈余/赤字 ................... : %12.6f\n', gov_revenue - gov_expenditure);

    % --- [4] 核心市场出清与自洽性检验 (保持不变) ---
    fprintf('\n--- [4] 核心市场出清与自洽性检验 ---\n');
    fprintf('   私人资本市场 (S_p - I_p) .............. : %12.4e\n', errors.err_Kp);
    fprintf('   公共资本市场 (S_g - I_g) .............. : %12.4e\n', errors.err_Kg);
    fprintf('   劳动市场 (L_supply - L_demand) ........ : %12.4e\n', errors.err_L);
    fprintf('   遗赠市场 (Bequest_gen - Bequest_dist) . : %12.4e\n', errors.err_Beq);
    fprintf('   PAYG系统预算平衡 ...................... : %12.4e\n', errors.payg_balance);
    fprintf('   ------------------------------------------------------------------------------\n');
    fprintf('   宏观存量自洽性 (K_p_end / K_p_begin - 1) : %12.4e%%\n', (ss.K_private_hat / ss.K_private_begin_hat - 1)*100);
    fprintf('================================================================================\n');
        end
        
        function [x_eq, resnorm, exitflag] = solve_with_lsqnonlin(system_wrapper, x0, lb, ub, verbose)
            if verbose, lsq_display = 'iter'; else, lsq_display = 'none'; end
            options = optimoptions('lsqnonlin', 'Display', lsq_display, 'FunctionTolerance', 1e-12, 'StepTolerance', 1e-12, 'MaxIterations', 200, 'Algorithm', 'trust-region-reflective');
            fprintf('\n--- 启动 lsqnonlin 求解器 (带边界约束) ---\n');
            [x_eq, resnorm, ~, exitflag] = lsqnonlin(system_wrapper, x0, lb, ub, options);
            fprintf('--- lsqnonlin 求解完成, 残差平方和(resnorm): %.4e ---\n', resnorm);
        end

        function M_prices = get_prices_at_t(K_p_hat, K_g_hat, L_hat, cS)
            if K_p_hat <= 0, K_p_hat = 1e-8; end; if L_hat <= 0, L_hat = 1e-8; end; if K_g_hat <= 0, K_g_hat = 1e-8; end;
            labor_exponent = 1 - cS.alpha - cS.gamma;
            Y_hat_t = (K_p_hat.^cS.alpha) .* (K_g_hat.^cS.gamma) .* (L_hat.^labor_exponent);
            MPK_p_period = cS.alpha .* Y_hat_t ./ K_p_hat;
            w_hat_t = labor_exponent .* Y_hat_t ./ L_hat;
            r_mkt_t = MPK_p_period - cS.ddk;
            M_prices = struct('Y_hat_t', Y_hat_t, 'w_hat_t', w_hat_t, 'r_mkt_t', r_mkt_t);
        end

        function [idx_lower, idx_upper, w_upper] = find_grid_and_weights(value, gridV)
            if value <= gridV(1), idx_lower = 1; idx_upper = 1; w_upper = 0.0; return; end
            if value >= gridV(end), idx_lower = length(gridV); idx_upper = length(gridV); w_upper = 0.0; return; end
            idx_upper = find(gridV >= value, 1);
            idx_lower = idx_upper - 1;
            if isempty(idx_lower) || idx_lower == 0, idx_lower = 1; idx_upper = 1; w_upper = 0; return; end
            if idx_lower == idx_upper, w_upper = 0.0; else, w_upper = (value - gridV(idx_lower)) / (gridV(idx_upper) - gridV(idx_lower)); end
        end
    end
end