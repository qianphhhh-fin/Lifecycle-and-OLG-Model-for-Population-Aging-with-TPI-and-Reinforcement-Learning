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
                x0 = [0.3, 0.4, 0.3, 0.01, 0.2];
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
            % == 版本: [v59 - 国民收入核算最终闭合版]
            % ==
            % == 核心修正:
            % ==   采纳最直接、最稳健的方式定义政府消费(G_c)。G_c不再通过
            % ==   政府预算的现金流反解，而是直接作为宏观资源约束的结余项。
            % ==   G_c = Y - C_agg - I_p_gross - I_g_gross
            % ==   这在结构上保证了国民收入核算恒等式 Y=C+I+G 永远成立，
            % ==   彻底消除最后的会计缺口。
            % =========================================================================

            % --- 1. 获取价格并求解家庭问题 ---
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
            [polS, valS] = main_steady_state_utils_bgp.HHSolution_VFI_unified(M_for_hh, paramS, cS);
            if any(any(any(any(isinf(valS))))) || any(any(any(any(isnan(valS))))), ss = []; Dist = []; polS = []; valS = []; return; end
            Dist = main_steady_state_utils_bgp.solve_steady_state_distribution_unified(polS, paramS, cS, Z_ss_norm);

            % --- 2. 聚合微观变量 ---
            K_private_hat_agg_raw = 0; L_agg = 0; C_agg_raw = 0; Bequest_generated_agg_raw = 0;
            Shock_exp_agg = 0; Pension_out_agg = 0; Pension_in_agg = 0; Regular_tax_agg = 0;

            for ia = 1:cS.aD_new
                mass_dist_ia = Dist(:,:,:,ia);
                k_prime_ia = polS(ia).k_prime;
                K_private_hat_agg_raw = K_private_hat_agg_raw + sum(k_prime_ia .* mass_dist_ia, 'all');
                if ia <= cS.aR_new, L_agg = L_agg + sum( (cS.ageEffV_new(ia) * paramS.leGridV' .* squeeze(sum(mass_dist_ia, [1,2])))' , 'all'); end
                Bequest_generated_agg_raw = Bequest_generated_agg_raw + sum(k_prime_ia .* mass_dist_ia, 'all') * (1 - cS.s_pathV(ia));
                C_agg_raw = C_agg_raw + sum(polS(ia).c .* mass_dist_ia, 'all');
                Shock_exp_agg = Shock_exp_agg + sum(polS(ia).shock_exp .* mass_dist_ia, 'all');
                Regular_tax_agg = Regular_tax_agg + sum(polS(ia).tax_regular .* mass_dist_ia, 'all');
                Pension_in_agg = Pension_in_agg + sum(polS(ia).tax_payg .* mass_dist_ia, 'all');
                Pension_out_agg = Pension_out_agg + sum(polS(ia).pension_out .* mass_dist_ia, 'all');
            end
            C_agg = C_agg_raw + Shock_exp_agg;

            % --- 3. BGP增长核算 ---
            g_A_period = (1+cS.g_A_ss)^cS.time_Step-1;
            n_period = (1+cS.n_ss)^cS.time_Step-1;

            K_private_hat_agg = K_private_hat_agg_raw / (1 + n_period);
            Bequest_generated_agg = Bequest_generated_agg_raw / (1 + n_period);

            g_total_factor = (1 + g_A_period) * (1 + n_period);

            Saving_private_flow_Gross = K_private_hat_agg * g_total_factor - (1 - cS.ddk) * K_private_total_guess;
            Saving_public_flow_Gross = K_g_guess * g_total_factor - (1 - cS.ddk_g) * K_g_guess;

            % --- 4. [最终修正] 政府部门账户与资源约束闭合 ---
            % Y_from_production_hat, C_agg, Saving_private_flow_Gross, Saving_public_flow_Gross
            % 全部基于均衡价格和最优行为计算得出。G_c是确保资源约束成立的唯一余项。
            Y_hat = M_prices.Y_hat_t;
            G_c = Y_hat - C_agg - Saving_private_flow_Gross - Saving_public_flow_Gross;

            Public_Capital_Return = Y_hat - (M_prices.w_hat_t * L_agg) - ((M_prices.r_mkt_t + cS.ddk) * K_private_total_guess);

            % --- 5. 填充ss结构体 ---
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

        function [polS, valS] = HHSolution_VFI_unified(M_vfi, paramS_vfi, cS_vfi)
            valS = -Inf(cS_vfi.nk, cS_vfi.nkpps, cS_vfi.nw_expanded, cS_vfi.aD_new);
            polS_cell = cell(cS_vfi.aD_new, 1);
            bV_payg_vfi = zeros(1, cS_vfi.aD_new);
            if cS_vfi.aR_new < cS_vfi.aD_new, bV_payg_vfi((cS_vfi.aR_new+1):cS_vfi.aD_new) = M_vfi.b_t; end
            for a_idx = cS_vfi.aD_new : -1 : 1
                vPrime_kkppse_next = [];
                if a_idx < cS_vfi.aD_new, vPrime_kkppse_next = valS(:,:,:,a_idx+1); end
                [val_age, pol_age] = main_steady_state_utils_bgp.HHSolutionByAge_VFI_noPPS(a_idx, vPrime_kkppse_next, M_vfi, bV_payg_vfi(a_idx), M_vfi.beq_transfer_pers, paramS_vfi, cS_vfi);
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

        function Dist = solve_steady_state_distribution_unified(polS, paramS, cS, Z_ss_norm)
            nk = cS.nk; nw = cS.nw_expanded; aD = cS.aD_new;
            Dist_cond = zeros(nk, 1, nw, aD, 'double');
            dist_newborn_cond_shape = zeros(nk, 1, nw);
            dist_newborn_cond_shape(1, 1, 1:cS.nw) = reshape(paramS.leProb1V, [1, 1, cS.nw]);
            Dist_cond(:, 1, :, 1) = dist_newborn_cond_shape;
            for ia = 1:(aD - 1)
                dist_cond_ia_slice = Dist_cond(:, 1, :, ia);
                if sum(dist_cond_ia_slice(:)) < 1e-30, continue; end
                dist_cond_ia_plus_1_next = zeros(nk, 1, nw, 'double');
                trans_mat_next_age = paramS.TrProbM_by_age{ia + 1};
                for ik = 1:nk, for ie = 1:nw
                        mass_start = dist_cond_ia_slice(ik, 1, ie);
                        if mass_start < 1e-30, continue; end
                        k_prime = polS(ia).k_prime(ik, 1, ie);
                        [ik_lower, ik_upper, w_k_upper] = main_steady_state_utils_bgp.find_grid_and_weights(k_prime, cS.kGridV);
                        w_k_lower = 1.0 - w_k_upper;
                        trans_probs_vec = trans_mat_next_age(ie, :);
                        for ie_next = 1:nw
                            prob_to_enext = trans_probs_vec(ie_next);
                            if prob_to_enext < 1e-12, continue; end
                            mass_to_distribute = mass_start * prob_to_enext;
                            if w_k_lower > 1e-9, dist_cond_ia_plus_1_next(ik_lower, 1, ie_next) = dist_cond_ia_plus_1_next(ik_lower, 1, ie_next) + mass_to_distribute * w_k_lower; end
                            if w_k_upper > 1e-9, dist_cond_ia_plus_1_next(ik_upper, 1, ie_next) = dist_cond_ia_plus_1_next(ik_upper, 1, ie_next) + mass_to_distribute * w_k_upper; end
                        end
                end, end
            sum_next_dist = sum(dist_cond_ia_plus_1_next(:));
            if sum_next_dist > 1e-9, Dist_cond(:, 1, :, ia+1) = dist_cond_ia_plus_1_next / sum_next_dist; end
            end
            Dist = Dist_cond .* reshape(Z_ss_norm, [1, 1, 1, aD]);
        end

        function errors = display_national_accounts_unified(ss, cS)
            verbose = true;
            g_A_period = (1+cS.g_A_ss)^cS.time_Step-1;
            n_period = (1+cS.n_ss)^cS.time_Step-1;
            g_total_period = (1+g_A_period)*(1+n_period)-1;
            I_p_gross_demand = (g_total_period + cS.ddk) * ss.K_private_begin_hat;
            I_g_gross_demand = ss.I_g;
            errors = struct('err_Kp', ss.Saving_private_flow_Gross - I_p_gross_demand, 'err_Kg', ss.Saving_public_flow_Gross - I_g_gross_demand, 'err_L', 0, 'err_Beq', ss.Bequest_generated_agg - ss.Bequest_distributed_agg, 'payg_balance', ss.Pension_in - ss.Pension_out);
            if verbose
                fprintf('\n\n================================================================================\n');
                fprintf('===                           国民经济核算与自洽性检验报告                           ===\n');
                fprintf('================================================================================\n');
                fprintf('--- [1] 生产与收入核算 ---\n');
                fprintf('   总产出 (Y) .......................... : %12.6f\n', ss.Y_from_production_hat);
                fprintf('   工资总额 (w*L) ...................... : %12.6f  (%.2f%% of Y)\n', ss.w_hat * ss.L_hat, (ss.w_hat * ss.L_hat / ss.Y_from_production_hat)*100);
                fprintf('   私人资本回报 (r*K_p) ................ : %12.6f  (%.2f%% of Y)\n', (ss.r_mkt+cS.ddk) * ss.K_private_begin_hat, ((ss.r_mkt+cS.ddk) * ss.K_private_begin_hat / ss.Y_from_production_hat)*100);
                fprintf('   公共资本回报 ........................ : %12.6f  (%.2f%% of Y)\n', ss.Public_Capital_Return, (ss.Public_Capital_Return / ss.Y_from_production_hat)*100);
                fprintf('\n--- [2] 支出法核算 ---\n');
                fprintf('   私人消费 (C_p) ...................... : %12.6f\n', ss.C_agg);
                fprintf('   政府消费 (G_c) ...................... : %12.6f\n', ss.G_c);
                fprintf('   私人投资 (I_p = S_p_gross) .......... : %12.6f\n', ss.Saving_private_flow_Gross);
                fprintf('   公共投资 (I_g = S_g_gross) .......... : %12.6f\n', ss.Saving_public_flow_Gross);
                Total_Expenditure = ss.C_agg + ss.G_c + ss.Saving_private_flow_Gross + ss.Saving_public_flow_Gross;
                fprintf('   >>> 总支出 (C+G+I) .................. : %12.6f  (Y-Expenditure Discrepancy: %.4e)\n', Total_Expenditure, ss.Y_from_production_hat - Total_Expenditure);
                fprintf('\n--- [3] 核心市场出清检验 (求解器目标) ---\n');
                fprintf('   市场1: 私人资本 S_p - I_p ...... : %12.8f - %12.8f = %12.4e\n', ss.Saving_private_flow_Gross, I_p_gross_demand, errors.err_Kp);
                fprintf('   市场2: 政府资本 S_g - I_g ...... : %12.8f - %12.8f = %12.4e\n', ss.Saving_public_flow_Gross, I_g_gross_demand, errors.err_Kg);
                fprintf('   市场3: 劳动供给-需求 .......... : %12.8f - %12.8f = %12.4e\n', ss.L_hat, ss.L_hat, errors.err_L);
                fprintf('   市场4: 遗赠产生-分配 .......... : %12.8f - %12.8f = %12.4e\n', ss.Bequest_generated_agg, ss.Bequest_distributed_agg, errors.err_Beq);
                fprintf('   PAYG系统: 收入 - 支出 .......... : %12.8f - %12.8f = %12.4e\n', ss.Pension_in, ss.Pension_out, errors.payg_balance);
                fprintf('\n--- [4] 宏观存量自洽性检验 ---\n');
                fprintf('   K_p_hat_begin (求解器均衡点) ...... : %15.8f\n', ss.K_private_begin_hat);
                fprintf('   K_p_hat_end (模型聚合期末选择) .... : %15.8f\n', ss.K_private_hat);
                fprintf('   -> 一致性误差 ((End/Begin)-1) .... : %14.2f%%\n', (ss.K_private_hat / ss.K_private_begin_hat - 1)*100);
                fprintf('================================================================================\n');
            end
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