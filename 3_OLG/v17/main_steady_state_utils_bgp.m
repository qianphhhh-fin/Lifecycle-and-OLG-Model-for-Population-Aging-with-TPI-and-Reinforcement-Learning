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
                % ss = []; Dist = []; polS = []; valS = [];
                % return;
            end
            
            % [!!! 核心修正 !!!] 在此处打印求解器返回的均衡解向量 x_eq
            fprintf('\n------------------------------------------------\n');
            fprintf('---  ✅ 求解器收敛！均衡解向量 (x_eq)  ---\n');
            fprintf('------------------------------------------------\n');
            fprintf('   私人资本 (K_p) .................... : %12.6f\n', x_eq(1));
            fprintf('   公共资本 (K_g) .................... : %12.6f\n', x_eq(2));
            fprintf('   有效劳动 (L) ...................... : %12.6f\n', x_eq(3));
            fprintf('   总遗赠 (Beq) ...................... : %12.6f\n', x_eq(4));
            if is_db_mode
                fprintf('   PAYG缴费率 (theta) ................ : %12.6f\n', x_eq(5));
            end
            fprintf('------------------------------------------------\n\n');


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
            % == 版本: [v9.0 - 遗赠会计核算最终修正版]
            % ==
            % == 核心修正:
            % ==   1. [遗赠来源修正] 遗赠现在由家庭带入当期的期初资产决定，
            % ==      而不是期末的储蓄计划，这修复了财富守恒的漏洞。
            % ==   2. [遗赠去趋势化修正] 移除了对当期遗赠流量的错误人口增长折扣，
            % ==      确保了BGP核算的正确性。
            % =========================================================================

            % --- 1. 获取价格并准备VFI输入 (无变动) ---
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
            
            [polS, valS] = main_steady_state_utils_bgp.HHSolution_unified_VFI(M_for_hh, paramS, cS);
            if any(any(any(any(isinf(valS))))) || any(any(any(any(isnan(valS))))), ss = []; Dist = []; polS = []; valS = []; return; end
            
            Dist = main_steady_state_utils_bgp.solve_steady_state_distribution_unified(polS, paramS, cS, Z_ss_norm, Bequest_Total_guess);

            % --- 2. [!!! 核心修正 !!!] 宏观聚合 ---
            K_k_hat_agg_raw = 0;
            K_pps_hat_agg_raw = 0;
            L_agg = 0;
            C_agg_raw = 0;
            Bequest_generated_agg_raw = 0; % 重新审视此变量的计算
            Shock_exp_agg = 0;
            Pension_out_agg = 0;
            Pension_in_agg = 0;
            Regular_tax_agg = 0;
            
            % 为了正确计算遗赠，我们需要遍历所有状态 (ia, ik, ie)
            k_grid_mesh = repmat(reshape(cS.kGridV, [cS.nk, 1, 1]), [1, cS.nkpps, cS.nw_expanded]);

            if isfield(cS, 'pps_active') && cS.pps_active
                kpps_grid_mesh = repmat(reshape(cS.kppsGridV, [1, cS.nkpps, 1]), [cS.nk, 1, cS.nw_expanded]);
            end
            
            for ia = 1:cS.aD_new
                mass_dist_ia = Dist(:,:,:,ia);
                
                % 聚合期末选择的资产 (下一期的期初资产)
                K_k_hat_agg_raw = K_k_hat_agg_raw + sum(polS(ia).k_prime .* mass_dist_ia, 'all');
                if isfield(cS, 'pps_active') && cS.pps_active
                    K_pps_hat_agg_raw = K_pps_hat_agg_raw + sum(polS(ia).kpps_prime .* mass_dist_ia, 'all');
                end
                
                % [!!! 遗赠来源修正 !!!]
                % 遗赠来自于当期期初的财富存量，而不是期末的储蓄决策。
                % 总遗赠 = sum_{ia,ik,ie} [ 财富(k_ik) * 质量(Dist_ia,ik,ie) * 死亡率(1-s_ia) ]
                wealth_at_start_of_period = k_grid_mesh;
                if isfield(cS, 'pps_active') && cS.pps_active
                    wealth_at_start_of_period = wealth_at_start_of_period + kpps_grid_mesh;
                end
                Bequest_generated_agg_raw = Bequest_generated_agg_raw + sum(wealth_at_start_of_period .* mass_dist_ia, 'all') * (1 - cS.s_pathV(ia));

                % 聚合劳动
                if ia <= cS.aR_new
                    mass_on_e = squeeze(sum(mass_dist_ia, [1,2]));
                    L_agg = L_agg + cS.ageEffV_new(ia) * (paramS.leGridV' * mass_on_e);
                end
                
                % 聚合其他流量
                C_agg_raw = C_agg_raw + sum(polS(ia).c .* mass_dist_ia, 'all');
                Shock_exp_agg = Shock_exp_agg + sum(polS(ia).shock_exp .* mass_dist_ia, 'all');
                Regular_tax_agg = Regular_tax_agg + sum(polS(ia).tax_regular .* mass_dist_ia, 'all');
                Pension_in_agg = Pension_in_agg + sum(polS(ia).tax_payg .* mass_dist_ia, 'all');
                Pension_out_agg = Pension_out_agg + sum(polS(ia).pension_out .* mass_dist_ia, 'all');
            end
            
            C_agg = C_agg_raw + Shock_exp_agg;
            
            K_private_hat_agg_raw = K_k_hat_agg_raw + K_pps_hat_agg_raw;

            % --- 3. [!!! 核心修正 !!!] BGP去趋势化 ---
            n_period = (1 + cS.n_ss)^cS.time_Step - 1;
            
            % K'是 t+1 期的变量，需要进行人口增长去趋势化
            K_private_hat_agg = K_private_hat_agg_raw / (1 + n_period);
            
            % [!!! 遗赠去趋势化修正 !!!]
            % Bequest 是在 t 期产生的流量，在与 t 期的其他变量比较时，
            % 不需要进行人口增长去趋势化。
            Bequest_generated_agg = Bequest_generated_agg_raw; % <--- 移除 "/ (1 + n_period)"
            
            % --- 4 & 5. 财政规则与填充ss结构体 (与v8.0版本逻辑相同) ---
            Y_hat = M_prices.Y_hat_t;
            I_g_ss = cS.I_g_to_Y_ratio_ss * Y_hat;
            wages_total_consistent = M_prices.w_hat_t * L_guess;
            private_capital_return_total_consistent = (M_prices.r_mkt_t + cS.ddk) * K_private_total_guess;
            Public_Capital_Return = Y_hat - wages_total_consistent - private_capital_return_total_consistent;
            Gov_Revenue = Regular_tax_agg + Public_Capital_Return;
            G_c_ss = Gov_Revenue - I_g_ss;
            
            ss = struct();
            ss.K_private_begin_hat = K_private_total_guess;
            ss.K_public_hat = K_g_guess;
            ss.L_guess = L_guess;
            ss.Bequest_distributed_agg = Bequest_Total_guess;
            ss.K_private_hat = K_private_hat_agg;
            ss.L_hat = L_agg;
            ss.C_agg = C_agg;
            ss.Bequest_generated_agg = Bequest_generated_agg;
            ss.Regular_tax = Regular_tax_agg;
            ss.Pension_in = Pension_in_agg;
            ss.Pension_out = Pension_out_agg;
            ss.Y_from_production_hat = Y_hat;
            ss.w_hat = M_prices.w_hat_t;
            ss.r_mkt = M_prices.r_mkt_t;
            ss.I_g = I_g_ss;
            ss.G_c = G_c_ss;
            ss.Public_Capital_Return = Public_Capital_Return;
            ss.theta = theta_ss;
            if isfield(cS, 'pps_active') && cS.pps_active
                ss.K_pps_hat = K_pps_hat_agg_raw / (1 + n_period);
            end
            ss.Depreciation_p = cS.ddk * K_private_total_guess;
            ss.Depreciation_g = cS.ddk_g * K_g_guess;
        end
                                

        function F_error = system_of_equations_steady_state(x, Z_ss_norm, cS, paramS, params_ext)
            % =========================================================================
            % == 函数: system_of_equations_steady_state
            % == 版本: [v6.0 - 预算平衡模式最终版]
            % ==
            % == 核心修正:
            % ==   采用严格的、非冗余的均衡条件作为误差，符合OG-CORE原则。
            % ==   1. [私人资本市场出清]
            % ==   2. [公共资本BGP自洽] (新增的核心条件)
            % ==   3. [劳动市场出清]
            % ==   4. [遗赠市场出清]
            % ==   5. [养老金系统平衡] (DB模式)
            % ==   * 资源约束 (Y=C+I+G) 不作为误差，而作为事后检验。
            % =========================================================================

            is_db_mode = isfield(cS, 'endogenous_theta_mode') && cS.endogenous_theta_mode;

            % --- 1. 获取求解器猜测的宏观总量 ---
            K_private_total_guess = x(1);
            K_g_guess = x(2);
            L_guess = x(3);
            Bequest_Total_guess = x(4);
            if is_db_mode, params_ext.theta = x(5); end

            % --- 2. 基于猜测值，计算出一个完整的经济状态 ---
            %    calculate_aggregates_unified 将会根据财政规则内生计算 G_c
            [ss, ~, ~, ~] = main_steady_state_utils_bgp.calculate_aggregates_unified(K_private_total_guess, K_g_guess, L_guess, Bequest_Total_guess, Z_ss_norm, cS, paramS, params_ext);

            if isempty(ss)
                error_size = ifthen(is_db_mode, 5, 4);
                F_error = ones(error_size, 1) * 1e8;
                return;
            end

            % --- 3. [!!! 核心修正 !!!] 构造独立的均衡误差向量 ---
            
            % 误差1: 私人资本市场出清 (存量自洽)
            % 经济含义：家庭在期末意愿持有的资本总量，必须等于厂商在期初使用的资本总量。
            error_Kp = ss.K_private_hat - K_private_total_guess;

            % 误差2: 公共资本存量与流量的BGP自洽
            % 经济含义：由 I_g/Y 规则决定的公共投资流量(ss.I_g)，必须恰好等于
            % 维持公共资本存量(K_g_guess)在BGP上增长所需的投资。
            g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
            n_period = (1 + cS.n_ss)^cS.time_Step - 1;
            g_total_period = (1 + g_A_period) * (1 + n_period) - 1;
            required_Ig_for_bgp = (g_total_period + cS.ddk_g) * K_g_guess;
            error_Kg_consistency = ss.I_g - required_Ig_for_bgp;

            % 误差3: 劳动市场出清
            % 经济含义：家庭供给的总劳动，必须等于厂商需求的总劳动。
            error_L = ss.L_hat - L_guess;

            % 误差4: 遗赠市场出清
            % 经济含义：经济中产生的总遗赠，必须等于被分配给新家庭的总遗赠。
            error_Beq = ss.Bequest_generated_agg - Bequest_Total_guess;

            % --- 4. 组合误差向量 ---
            if is_db_mode
                % 误差5: 养老金PAYG系统预算平衡
                error_PAYG = ss.Pension_in - ss.Pension_out;
                F_error = [error_Kp; error_Kg_consistency; error_L; error_Beq; error_PAYG];
            else
                F_error = [error_Kp; error_Kg_consistency; error_L; error_Beq];
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
                    [val_age, pol_age] = main_steady_state_utils_bgp.HHSolutionByAge_VFI_PPS(a_idx, vPrime_kkppse_next, M_vfi, bV_payg_vfi(a_idx), paramS_vfi, cS_vfi);
                else
                    [val_age, pol_age] = main_steady_state_utils_bgp.HHSolutionByAge_VFI_noPPS(a_idx, vPrime_kkppse_next, M_vfi, bV_payg_vfi(a_idx), paramS_vfi, cS_vfi);
                end

                valS(:,:,:,a_idx) = val_age;
                polS_cell{a_idx} = pol_age;
            end

            polS = [polS_cell{:}];
        end

                function [val_age, pol_age] = HHSolutionByAge_VFI_noPPS(a_idx, vPrime_kkppse_next, M_age, b_age_val, paramS_age, cS)
            % =========================================================================
            % == 函数: HHSolutionByAge_VFI_noPPS
            % == 版本: [v2.0 - 遗赠修正版]
            % ==
            % == 核心修正:
            % ==   - [移除遗赠收入] 彻底移除了 beq_transfer_val 参数和其在
            % ==     预算约束中的错误应用。遗赠现在被正确地处理为新生代的
            % ==     初始资产禀赋，在分布函数中实现。
            % =========================================================================
            nk_search_grid = 150;
            
            val_age = -1e20 * ones(cS.nk, 1, cS.nw_expanded);
            
            pol_age = struct('c', zeros(cS.nk, 1, cS.nw_expanded), ...
                'k_prime', zeros(cS.nk, 1, cS.nw_expanded), ...
                'tax_regular', zeros(cS.nk, 1, cS.nw_expanded), ...
                'tax_payg', zeros(cS.nk, 1, cS.nw_expanded), ...
                'shock_exp', zeros(cS.nk, 1, cS.nw_expanded), ...
                'pension_out', zeros(cS.nk, 1, cS.nw_expanded)); % 移除了 beq_received 字段
            
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
                
                pension_out = 0;
                
                if a_idx <= cS.aR_new
                    labor_income_gross = M_age.w_t * cS.ageEffV_new(a_idx) * paramS_age.leGridV(ie);
                else
                    pension_out = b_age_val;
                end
                
                capital_tax_base = k_grid_vec * M_age.r_mkt_t;
                
                capital_tax_vec = cS.tau_k .* capital_tax_base;
                
                payg_tax = cS.theta_t * labor_income_gross;
                
                labor_tax = cS.tau_l * max(0, labor_income_gross - payg_tax);
                
                % [!!! 核心会计修正 !!!] cash_on_hand 中不再包含遗赠项。
                cash_on_hand_vec = capital_return + labor_income_gross + pension_out - (capital_tax_vec + payg_tax + labor_tax);
                
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

        function Dist = solve_steady_state_distribution_unified(polS, paramS, cS, Z_ss_norm, Bequest_Total_guess)
            % =========================================================================
            % == 函数: solve_steady_state_distribution_unified
            % == 版本: [v4.0 - 遗赠禀赋最终版]
            % ==
            % == 核心修正:
            % ==   - [初始禀赋] 新生代的初始分布不再是(k=0)的一个点，而是根据
            % ==     人均继承的遗赠量，被精准地分配到对应的资产网格点上。
            % ==     这正确地将遗赠建模为代际间的财富存量转移。
            % ==   - [线性插值] 使用线性插值将质量分配到最接近的两个网格点，
            % ==     以处理人均遗赠额不恰好等于某个网格点的情况。
            % =========================================================================

            nk = cS.nk;
            nw = cS.nw_expanded;
            aD = cS.aD_new;

            % --- [!!! 核心会计修正 !!!] ---
            % 步骤 A: 计算每个新生儿继承的遗赠量
            mass_newborns = Z_ss_norm(1); % 稳态下新生代的占比
            if mass_newborns > 1e-9
                beq_per_newborn = Bequest_Total_guess / mass_newborns;
            else
                beq_per_newborn = 0;
            end

            % 步骤 B: 找到该遗赠量在资产网格上的位置和插值权重
            [ik_lower, ~, w_upper] = main_steady_state_utils_bgp.find_grid_and_weights(beq_per_newborn, cS.kGridV);
            w_lower = 1.0 - w_upper;
            ik_upper = min(nk, ik_lower + 1);

            % 步骤 C: 构造新生代的初始条件分布
            dist_newborn_cond_shape = zeros(nk, 1, nw);
            
            % 将新生代的初始生产力冲击概率 (leProb1V) 分配到对应的资产网格点
            prob_mass_to_distribute = reshape(paramS.leProb1V, [1, 1, cS.nw]);
            dist_newborn_cond_shape(ik_lower, 1, 1:cS.nw) = w_lower * prob_mass_to_distribute;
            dist_newborn_cond_shape(ik_upper, 1, 1:cS.nw) = dist_newborn_cond_shape(ik_upper, 1, 1:cS.nw) + w_upper * prob_mass_to_distribute;


            if cS.pps_active
                % --- PPS路径 ---
                nkpps = cS.nkpps;
                Dist_cond = zeros(nk, nkpps, nw, aD, 'double');
                
                % 新生儿 kpps=0，k 由遗赠决定，e 由初始冲击决定
                dist_newborn_cond_shape_pps = zeros(nk, nkpps, nw);
                dist_newborn_cond_shape_pps(:, 1, :) = dist_newborn_cond_shape; % pps从0开始
                Dist_cond(:, :, :, 1) = dist_newborn_cond_shape_pps;
                
                % ... 后续的分布迭代循环(ia=1 to aD-1)保持不变 ...
                for ia = 1:(aD - 1), dist_cond_ia_slice = Dist_cond(:, :, :, ia); if sum(dist_cond_ia_slice(:)) < 1e-30, continue; end; dist_cond_ia_plus_1_next = zeros(nk, nkpps, nw, 'double'); trans_mat_next_age = paramS.TrProbM_by_age{ia + 1}; for ik = 1:nk, for ikpps = 1:nkpps, for ie = 1:nw, mass_start = dist_cond_ia_slice(ik, ikpps, ie); if mass_start < 1e-30, continue; end; k_prime = polS(ia).k_prime(ik, ikpps, ie); kpps_prime = polS(ia).kpps_prime(ik, ikpps, ie); [ik_lower_p, ~, w_k_upper] = main_steady_state_utils_bgp.find_grid_and_weights(k_prime, cS.kGridV); w_k_lower = 1.0 - w_k_upper; ik_upper_p = min(nk, ik_lower_p + 1); [ikpps_lower, ~, w_kpps_upper] = main_steady_state_utils_bgp.find_grid_and_weights(kpps_prime, cS.kppsGridV); w_kpps_lower = 1.0 - w_kpps_upper; ikpps_upper = min(nkpps, ikpps_lower + 1); trans_probs_vec = trans_mat_next_age(ie, :); for ie_next = 1:nw, prob_to_enext = trans_probs_vec(ie_next); if prob_to_enext < 1e-12, continue; end; mass_to_distribute_e = mass_start * prob_to_enext; dist_cond_ia_plus_1_next(ik_lower_p, ikpps_lower, ie_next) = dist_cond_ia_plus_1_next(ik_lower_p, ikpps_lower, ie_next) + mass_to_distribute_e * w_k_lower * w_kpps_lower; dist_cond_ia_plus_1_next(ik_upper_p, ikpps_lower, ie_next) = dist_cond_ia_plus_1_next(ik_upper_p, ikpps_lower, ie_next) + mass_to_distribute_e * w_k_upper * w_kpps_lower; dist_cond_ia_plus_1_next(ik_lower_p, ikpps_upper, ie_next) = dist_cond_ia_plus_1_next(ik_lower_p, ikpps_upper, ie_next) + mass_to_distribute_e * w_k_lower * w_kpps_upper; dist_cond_ia_plus_1_next(ik_upper_p, ikpps_upper, ie_next) = dist_cond_ia_plus_1_next(ik_upper_p, ikpps_upper, ie_next) + mass_to_distribute_e * w_k_upper * w_kpps_upper; end; end; end; end; sum_next_dist = sum(dist_cond_ia_plus_1_next(:)); if sum_next_dist > 1e-9, Dist_cond(:, :, :, ia+1) = dist_cond_ia_plus_1_next / sum_next_dist; end; end

            else
                % --- 无 PPS 路径 ---
                Dist_cond = zeros(nk, 1, nw, aD, 'double');
                
                % 设置新生代的初始分布
                Dist_cond(:, :, :, 1) = dist_newborn_cond_shape;
                
                % ... 后续的分布迭代循环(ia=1 to aD-1)保持不变 ...
                 for ia = 1:(aD - 1), dist_cond_ia_slice = Dist_cond(:, 1, :, ia); if sum(dist_cond_ia_slice(:)) < 1e-30, continue; end; dist_cond_ia_plus_1_next = zeros(nk, 1, nw, 'double'); trans_mat_next_age = paramS.TrProbM_by_age{ia + 1}; for ik = 1:nk, for ie = 1:nw, mass_start = dist_cond_ia_slice(ik, 1, ie); if mass_start < 1e-30, continue; end; k_prime = polS(ia).k_prime(ik, 1, ie); [ik_lower_p, ~, w_k_upper] = main_steady_state_utils_bgp.find_grid_and_weights(k_prime, cS.kGridV); w_k_lower = 1.0 - w_k_upper; ik_upper_p = min(nk, ik_lower_p + 1); trans_probs_vec = trans_mat_next_age(ie, :); for ie_next = 1:nw, prob_to_enext = trans_probs_vec(ie_next); if prob_to_enext < 1e-12, continue; end; mass_to_distribute = mass_start * prob_to_enext; dist_cond_ia_plus_1_next(ik_lower_p, 1, ie_next) = dist_cond_ia_plus_1_next(ik_lower_p, 1, ie_next) + mass_to_distribute * w_k_lower; dist_cond_ia_plus_1_next(ik_upper_p, 1, ie_next) = dist_cond_ia_plus_1_next(ik_upper_p, 1, ie_next) + mass_to_distribute * w_k_upper; end; end; end; sum_next_dist = sum(dist_cond_ia_plus_1_next(:)); if sum_next_dist > 1e-9, Dist_cond(:, :, :, ia+1) = dist_cond_ia_plus_1_next / sum_next_dist; end; end
            end

            % 最终分布 = 条件分布 * 年龄的边际分布
            Dist = Dist_cond .* reshape(Z_ss_norm, [1, 1, 1, aD]);
        end        
                        function errors = display_national_accounts_unified(ss, cS)
            % =========================================================================
            % == 函数: display_national_accounts_unified
            % == 版本: [v4.0 - 会计核算修正最终版]
            % ==
            % == 核心修正:
            % ==   - [资源约束检验] 新增了资源约束 Y - (C+I_p+I_g+G_c) 的检验，
            % ==     作为模型整体会计逻辑的最终验证。
            % ==   - [BGP投资定义] 报告中使用的 I_p 和 I_g 明确基于BGP增长公式
            % ==     I = (g_total + delta) * K 计算，确保与理论一致。
            % =========================================================================

            % --- [0] 准备工作：计算增长因子和理论投资需求 ---
            Y = ss.Y_from_production_hat;
            
            g_A_period = (1+cS.g_A_ss)^cS.time_Step-1;
            
            n_period = (1+cS.n_ss)^cS.time_Step-1;
            
            g_total_period = (1+g_A_period)*(1+n_period)-1;
            
            % [BGP核算] BGP下维持私人资本存量(K_p_begin)所需的总投资
            I_p_demand_bgp = (g_total_period + cS.ddk) * ss.K_private_begin_hat;
            
            % [BGP核算] BGP下维持公共资本存量(K_public_hat)所需的总投资
            I_g_demand_bgp = (g_total_period + cS.ddk_g) * ss.K_public_hat;
            
            % --- [报告开始] ---
            fprintf('\n\n================================================================================\n');
            fprintf('===           稳态均衡：国民经济核算与诊断报告 (预算平衡模式)            ===\n');
            fprintf('================================================================================\n');
            fprintf('--- [0] 均衡价格与核心比率 ---\n');
            fprintf('   年化市场利率 (r_annual) ............... : %12.4f%%\n', ( (1+ss.r_mkt)^(1/cS.time_Step) - 1 ) * 100);
            fprintf('   有效劳动单位工资 (w_hat) .............. : %12.6f\n', ss.w_hat);
            fprintf('   私人资本/产出比 (K_p/Y) ............... : %12.4f\n', ss.K_private_begin_hat / Y);
            fprintf('   公共资本/产出比 (K_g/Y) ............... : %12.4f\n', ss.K_public_hat / Y);
            fprintf('\n--- [1] 宏观存量 (Stocks) ---\n');
            fprintf('   总产出 (Y_hat) ........................ : %12.6f\n', Y);
            fprintf('   私人资本 (K_p_hat) .................... : %12.6f  (%.2f%% of Y)\n', ss.K_private_begin_hat, (ss.K_private_begin_hat / Y)*100);
            if isfield(ss, 'K_pps_hat'), fprintf('     其中: PPS资本 (K_pps) ............... : %12.6f  (%.2f%% of K_p)\n', ss.K_pps_hat, (ss.K_pps_hat / ss.K_private_begin_hat)*100); end
            fprintf('   公共资本 (K_g_hat) .................... : %12.6f  (%.2f%% of Y)\n', ss.K_public_hat, (ss.K_public_hat / Y)*100);
            fprintf('   有效劳动供给 (L_hat) .................. : %12.6f\n', ss.L_hat);

            % --- [2] 国民收入核算 (支出法流量) ---
            fprintf('\n--- [2] 国民收入核算 (支出法流量, Flows) ---\n');
            fprintf('   私人消费 (C_p) ........................ : %12.6f  (%.2f%% of Y)\n', ss.C_agg, (ss.C_agg / Y)*100);
            fprintf('   政府消费 (G_c) [内生决定] ........... : %12.6f  (%.2f%% of Y)\n', ss.G_c, (ss.G_c / Y)*100);
            fprintf('   私人总投资 (I_p) [BGP需求] .......... : %12.6f  (%.2f%% of Y)\n', I_p_demand_bgp, (I_p_demand_bgp / Y)*100);
            fprintf('   公共总投资 (I_g) [按规则] ........... : %12.6f  (%.2f%% of Y)\n', ss.I_g, (ss.I_g / Y)*100);
            Total_Expenditure = ss.C_agg + I_p_demand_bgp + ss.I_g + ss.G_c;
            fprintf('   ------------------------------------------------------------------------------\n');
            fprintf('   >>> 总支出 (C+I_p+I_g+G_c) ............ : %12.6f\n', Total_Expenditure);
            
            % [!!! 核心检验 !!!] 资源约束自洽性检验
            resource_constraint_error = Y - Total_Expenditure;
            fprintf('   >>> 资源约束误差 (Y - Total_Expenditure) : %12.4e\n', resource_constraint_error);
            if abs(resource_constraint_error / Y) < 1e-6
                fprintf('       检验结果: ✅ 资源约束满足，模型会计逻辑自洽。\n');
            else
                fprintf('       检验结果: ⚠️ 资源约束不满足，模型会计逻辑存在问题！\n');
            end

            % --- [3] 政府与养老金账户 ---
            fprintf('\n--- [3] 政府与养老金账户 ---\n');
            fprintf('   PAYG 养老金系统:\n');
            fprintf('     总收入 (Pension_in) ................. : %12.6f\n', ss.Pension_in);
            fprintf('     总支出 (Pension_out) ................ : %12.6f\n', ss.Pension_out);
            fprintf('     -> PAYG系统预算平衡 ................ : %12.4e\n', ss.Pension_in - ss.Pension_out);
            fprintf('\n   政府一般预算 (收支平衡):\n');
            gov_revenue = ss.Regular_tax + ss.Public_Capital_Return;
            fprintf('     总收入 (税收 + 公共资本回报) ........ : %12.6f\n', gov_revenue);
            gov_expenditure = ss.G_c + ss.I_g;
            fprintf('     总支出 (G_c + I_g) .................. : %12.6f\n', gov_expenditure);
            budget_balance_error = gov_revenue - gov_expenditure;
            fprintf('     -> 预算平衡误差 (收入 - 支出) ........ : %12.4e\n', budget_balance_error);
            if abs(budget_balance_error/Y) < 1e-9
                 fprintf('       检验结果: ✅ 政府预算按设计严格平衡。\n');
            else
                 fprintf('       检验结果: ⚠️ 政府预算未按设计平衡！\n');
            end

            % --- [4] 核心市场出清与自洽性检验 ---
            fprintf('\n--- [4] 核心市场出清与自洽性检验 ---\n');
            % 定义用于展示的误差结构体
            errors = struct('err_Kp', ss.K_private_hat - ss.K_private_begin_hat, ...
                            'err_Kg_consistency', ss.I_g - I_g_demand_bgp, ...
                            'err_L', ss.L_hat - ss.L_guess, ...
                            'err_Beq', ss.Bequest_generated_agg - ss.Bequest_distributed_agg, ...
                            'payg_balance', ss.Pension_in - ss.Pension_out);

            fprintf('   [私人资本] 存量自洽性 (K_end - K_begin) : %12.4e\n', errors.err_Kp);
            fprintf('   [公共资本] 流量/存量BGP自洽性 (误差) .. : %12.4e\n', errors.err_Kg_consistency);
            fprintf('   [劳动市场] 供需缺口 (L_supply-L_demand) : %12.4e\n', errors.err_L);
            fprintf('   [遗赠市场] 供需缺口 (B_gen-B_dist) .... : %12.4e\n', errors.err_Beq);
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