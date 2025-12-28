classdef main_steady_state_utils
    % =========================================================================
    % == 类说明: main_steady_state_utils
    % ==
    % == 包含两套稳态求解器：
    % == 1. 原版本（不含PPS）: solve_steady_state_complete [使用固定网格]
    % == 2. PPS扩展版本: solve_steady_state_complete_with_pps [使用自适应网格]
    % ==
    % == 核心特性：
    % == - 智能自适应网格系统：外部基于智能猜测值预设网格，避免求解器内反复调整
    % == - 多种求解器支持：fsolve, surrogateopt, hybrid, robust
    % == - 智能边界设置：求解器边界自动调整以匹配网格范围
    % == - 双轨制设计：初始稳态用固定网格，新稳态用预设自适应网格
    % ==
    % == 使用示例：
    % == % 求解不含PPS的稳态（用于初始稳态）
    % == [ss_old, ~, ~, ~] = main_steady_state_utils.solve_steady_state_complete(cS, paramS, params_ext, true);
    % ==
    % == % 求解含PPS的稳态（用于终点稳态，带智能自适应网格）
    % == cS_new.pps_active = true;
    % == % 可选：基于智能猜测值预设自适应网格
    % == k_max_adaptive = 10 * initial_guess_K_p;  % 网格上限为猜测值的10倍
    % == cS_new = model_setup_utils.generateGrids(cS_new, 'k_max', k_max_adaptive);
    % == [ss_new, ~, V_new, k_pol_new, cPps_pol_new] = main_steady_state_utils.solve_steady_state_complete_with_pps(cS_new, paramS, params_ext, true);
    % ==
    % == % 使用不同求解器的示例
    % == [ss_fsolve, ~, ~, ~, ~] = main_steady_state_utils.solve_steady_state_complete_with_pps(cS, paramS, params_ext, true, [], 'fsolve');
    % == [ss_surrogate, ~, ~, ~, ~] = main_steady_state_utils.solve_steady_state_complete_with_pps(cS, paramS, params_ext, true, [], 'surrogateopt');
    % == [ss_hybrid, ~, ~, ~, ~] = main_steady_state_utils.solve_steady_state_complete_with_pps(cS, paramS, params_ext, true, [], 'hybrid');
    % == [ss_robust, ~, ~, ~, ~] = main_steady_state_utils.solve_steady_state_complete_with_pps(cS, paramS, params_ext, true, [], 'robust');
    % ==
    % == % 测试功能
    % == main_steady_state_utils.test_adaptive_grid_system();    % 测试智能自适应网格系统
    % == main_steady_state_utils.test_pps_steady_state_solver(); % 测试PPS稳态求解器
    % == main_steady_state_utils.test_solver_comparison();       % 测试不同求解器的性能比较
    % =========================================================================

    methods (Static)


        % =======================================================
        % == 阶段一：稳态求解器 (作为过渡态的起点)
        % =======================================================


        function [ss, Dist, V, kPolM] = solve_steady_state_complete(cS_ss, paramS, params_ext, verbose)
            % [修改] 增加 verbose 输入参数
            if nargin < 4, verbose = true; end % 默认详细输出

            % 1. 使用模型内生的人口分布
            age_mass = ones(cS_ss.aD_new, 1);
            for ia = 1:(cS_ss.aD_new - 1)
                age_mass(ia+1) = age_mass(ia) * cS_ss.s_pathV(ia);
            end
            % 如果外部提供了特定年份的人口分布(例如校准时)，则使用外部的
            if isfield(params_ext, 'Z') && ~isempty(params_ext.Z)
                Z_ss_norm = params_ext.Z;
            else
                Z_ss_norm = age_mass / sum(age_mass);
            end

            % 2. 覆盖cS中的外生变量
            cS_ss.A = params_ext.A;
            if isfield(params_ext, 'g_A_ss'), cS_ss.g_A_ss = params_ext.g_A_ss; end
            cS_ss.theta_path = params_ext.theta; % 确保是标量

            % [新增] 计算稳态的隐含存活率向量
            prob_survive_implied_ss0 = zeros(cS_ss.aD_new, 1);
            for ia = 1:(cS_ss.aD_new - 1)
                if Z_ss_norm(ia) > 1e-12
                    prob_survive_implied_ss0(ia) = Z_ss_norm(ia+1) / Z_ss_norm(ia);
                else
                    prob_survive_implied_ss0(ia) = 0;
                end
            end
            prob_survive_implied_ss0(cS_ss.aD_new) = 0;
            cS_ss.prob_survive_implied_ss0 = prob_survive_implied_ss0;

            % 3. [核心] 调用稳态求解器，并传递 verbose 参数
            % [修改] 确保函数能返回 V 和 kPolM
            [ss, eq_found, Dist, ~, V, kPolM] = ...
                main_steady_state_utils.solve_steady_state_iter_Kg(Z_ss_norm, cS_ss, paramS, verbose);

            if ~eq_found
                warning('稳态求解失败！');
                ss = []; Dist = []; V = []; kPolM = []; % 求解失败返回空
            end
        end

        function [ss, eq_found, Dist, k_prime_idx, V, kPolM] = solve_steady_state_iter_Kg(Z_ss_norm, cS, paramS, verbose)
            % [修改] 增加 verbose 输入参数并修改返回值
            if nargin < 4, verbose = true; end

            system_wrapper = @(x) main_steady_state_utils.system_of_equations_Kg(x, Z_ss_norm, cS, paramS);

            k_p_guess_initial = 3.5;
            k_g_guess_initial = 1.0;
            x0 = [k_p_guess_initial, k_g_guess_initial];

            if verbose
                fsolve_display = 'iter';
            else
                fsolve_display = 'none'; % 安静模式
            end
            options = optimoptions('fsolve', 'Display', fsolve_display, 'TolFun', 1e-9, 'TolX', 1e-9, 'MaxIterations', 500);

            if verbose, fprintf('\n--- 启动 fsolve 求解器 (求解 [K_p, K_g]) ---\n'); end
            [x_eq, ~, exitflag] = fsolve(system_wrapper, x0, options);
            if verbose, fprintf('--- fsolve 求解完成 ---\n'); end

            eq_found = (exitflag > 0);
            if ~eq_found
                if verbose, warning('fsolve 未能找到均衡解 (exitflag = %d)', exitflag); end
                ss=[]; Dist=[]; k_prime_idx=[]; V=[]; kPolM=[];
                return;
            end

            K_p_eq = x_eq(1);
            K_g_eq = x_eq(2);

            % --- 求解成功后，进行最后一次计算以获得所有变量 ---
            % [修改] 调用已更新的函数，获取 V 和 kPolM
            [~, ss, Dist, k_prime_idx, V, kPolM] = main_steady_state_utils.calculate_aggregates_for_Kg(K_p_eq, K_g_eq, Z_ss_norm, cS, paramS);

            % --- [修改] 根据 verbose 决定输出详细程度 ---
            if verbose
                main_steady_state_utils.display_national_accounts_gov_investment(ss, Dist, k_prime_idx, cS, paramS);
            else
                % 在校准循环中，只做一个快速检查
                main_steady_state_utils.check_national_accounts(ss, Dist, k_prime_idx, cS, paramS);
            end
        end

        function F_error = system_of_equations_Kg(x, Z_ss_norm, cS, paramS)
            % [vPaper.5 - 政府投资版]
            % 目的: 为 fsolve 提供一个误差向量 F_error = [error_Kp; error_Kg]

            K_p_guess = x(1);
            K_g_guess = x(2);

            % 调用核心计算函数，得到模型隐含的宏观聚合量
            [K_p_model, ss] = main_steady_state_utils.calculate_aggregates_for_Kg(K_p_guess, K_g_guess, Z_ss_norm, cS, paramS);

            % 方程1的误差: 私人资本供给 - 私人资本需求
            % (fsolve 需要 error=0，所以我们定义为 guess - model_outcome)
            error_Kp = K_p_guess - K_p_model;

            % 方程2的误差: 公共投资 - 公共资本折旧
            Gov_Revenue_total = ss.Regular_tax + ss.Bequest_tax;
            I_g_model = cS.lambda_g * Gov_Revenue_total;
            Depreciation_g_model = cS.ddk_g * K_g_guess;
            error_Kg = I_g_model - Depreciation_g_model;

            F_error = [error_Kp; error_Kg];
        end

                function [K_p_model_out, ss, Dist, k_prime_idx, V, kPolM] = calculate_aggregates_for_Kg(K_p_guess, K_g_guess, Z_ss_norm, cS, paramS)
            % [vPaper.5 - 政府投资版]
            % [修改] 此函数现在是核心，计算并返回所有稳态结果，包括 V 和 kPolM
            % [注意] 初始稳态使用固定网格，不使用自适应网格系统

            % ... (函数前半部分，计算劳动供给 L_guess 的循环保持不变) ...
            if K_p_guess <= 0, K_p_guess = 1e-8; end
            if K_g_guess <= 0, K_g_guess = 1e-8; end
            A_ss = cS.A; theta_ss = cS.theta_path(1);

            L_iter_tol = 1e-7; L_iter_max = 50; L_damping = 0.5;

            mean_e_by_age = zeros(cS.aD_new,1);
            e_dist_by_age = zeros(cS.aD_new, cS.nw_expanded);
            e_dist_by_age(1, 1:cS.nw) = paramS.leProb1V';
            mean_e_by_age(1) = e_dist_by_age(1, :) * paramS.leGridV(:);
            for ia = 1:(cS.aD_new - 1)
                e_dist_by_age(ia+1, :) = e_dist_by_age(ia, :) * paramS.TrProbM_by_age{ia+1};
                mean_e_by_age(ia+1) = e_dist_by_age(ia+1, :) * paramS.leGridV(:);
            end
            L_guess = 0;
            for ia = 1:cS.aR_new
                L_guess = L_guess + Z_ss_norm(ia) * cS.ageEffV_new(ia) * mean_e_by_age(ia);
            end

            for iter = 1:L_iter_max
                M_prices = main_steady_state_utils.get_prices_at_t(K_p_guess, K_g_guess, L_guess, A_ss, cS);
                M_for_hh = M_prices;

                total_wage_bill = M_prices.w_t * L_guess;
                mass_retirees_ss = sum(Z_ss_norm((cS.aR_new+1):end));
                if mass_retirees_ss > 1e-9
                    M_for_hh.b_t = (theta_ss * total_wage_bill) / mass_retirees_ss;
                else, M_for_hh.b_t = 0; end

                cS_ss = cS; cS_ss.theta_t = theta_ss;

                % [修改] 在这里只计算一次家庭问题，以找到劳动供给
                [~, temp_kPolM, ~, ~] = main_steady_state_utils.HHSolution_VFI(M_for_hh, paramS, cS_ss);
                temp_k_prime_idx = main_steady_state_utils.get_policy_index_matrix(temp_kPolM, cS_ss);
                temp_Dist = main_steady_state_utils.solve_steady_state_distribution(temp_k_prime_idx, paramS, cS_ss, Z_ss_norm);

                L_model = 0;
                for ia = 1:cS.aR_new, for ie = 1:cS.nw_expanded, for ik = 1:cS.nk
                            mass = temp_Dist(ik, ie, ia);
                            if mass > 0
                                epsilon_val = paramS.leGridV(ie);
                                L_model = L_model + (cS.ageEffV_new(ia) * epsilon_val) * mass;
                            end
                end,end,end

        L_error = abs(L_model - L_guess);
        if L_error < L_iter_tol, break; end
        L_guess = L_damping * L_guess + (1 - L_damping) * L_model;
            end

            % --- [核心修改] 劳动供给收敛后，计算最终的 V, kPol, Dist 和宏观量 ---
            M_prices = main_steady_state_utils.get_prices_at_t(K_p_guess, K_g_guess, L_guess, A_ss, cS);
            M_for_hh = M_prices;
            total_wage_bill = M_prices.w_t * L_guess;
            mass_retirees_ss = sum(Z_ss_norm((cS.aR_new+1):end));
            if mass_retirees_ss > 1e-9, M_for_hh.b_t = (theta_ss * total_wage_bill) / mass_retirees_ss;
            else, M_for_hh.b_t = 0; end

            % 调用VFI得到最终的价值和政策函数
            [~, kPolM, ~, V] = main_steady_state_utils.HHSolution_VFI(M_for_hh, paramS, cS_ss);
            k_prime_idx = main_steady_state_utils.get_policy_index_matrix(kPolM, cS_ss);
            Dist = main_steady_state_utils.solve_steady_state_distribution(k_prime_idx, paramS, cS_ss, Z_ss_norm);



            % [修改] 使用最终的分布和决策规则来聚合，并确保使用稳态隐含存活率
            [K_p_model_out, C_utility_final, Tax_final, Bequest_tax_final, ~, ~, ~, ~] = ...
                main_steady_state_utils.aggregate_expenditure_shock(Dist, k_prime_idx, M_for_hh, cS_ss, paramS);

            % 填充完整的 ss 结构体
            ss = struct();
            ss.K_private = K_p_guess;
            ss.K_public = K_g_guess;
            ss.K_total = K_p_guess + K_g_guess;
            ss.L = L_guess;
            ss.Y_from_production = M_prices.Y_t;
            ss.w = M_prices.w_t;
            ss.r_mkt = M_prices.r_mkt_t;
            ss.b = M_for_hh.b_t;
            ss.Bequest_tax = Bequest_tax_final;
            ss.Regular_tax = Tax_final;
        end

        function [cPolM, kPolM, cPpsPolM, valM] = HHSolution_VFI(M_vfi, paramS_vfi, cS_vfi)
            valM = -Inf(cS_vfi.nk, cS_vfi.nw_expanded, cS_vfi.aD_new);
            cPolM = zeros(cS_vfi.nk, cS_vfi.nw_expanded, cS_vfi.aD_new);
            kPolM = zeros(cS_vfi.nk, cS_vfi.nw_expanded, cS_vfi.aD_new);
            cPpsPolM = zeros(cS_vfi.nk, cS_vfi.nw_expanded, cS_vfi.aD_new);
            bV_payg_vfi = zeros(1, cS_vfi.aD_new);
            if cS_vfi.aR_new < cS_vfi.aD_new
                bV_payg_vfi((cS_vfi.aR_new+1):cS_vfi.aD_new) = M_vfi.b_t;
            end
            for a_idx = cS_vfi.aD_new : -1 : 1
                vPrime_kkppse_next = [];
                if a_idx < cS_vfi.aD_new
                    vPrime_kkppse_next = valM(:,:,a_idx+1);
                end
                [cPolM(:,:,a_idx), kPolM(:,:,a_idx), cPpsPolM(:,:,a_idx), valM(:,:,a_idx)] = ...
                    main_steady_state_utils.HHSolutionByAge_VFI_ExpenditureShock_Vectorized(a_idx, vPrime_kkppse_next, M_vfi, bV_payg_vfi(a_idx), paramS_vfi, cS_vfi);
            end
        end

        function [cPol_age_q, kPol_age, cPpsPol_age_choice, val_age] = HHSolutionByAge_VFI_ExpenditureShock(a_idx, vPrime_kkppse_next, M_age, b_age_val, paramS_age, cS)
            % [vPaper.4.9 - 暖光遗赠动机版]
            % 初始稳态版本：不含PPS决策

            val_age = -1e20 * ones(cS.nk, cS.nw_expanded);
            cPol_age_q = zeros(cS.nk, cS.nw_expanded);
            kPol_age = zeros(cS.nk, cS.nw_expanded);
            cPpsPol_age_choice = zeros(cS.nk, cS.nw_expanded);

            % 最后一期处理
            if a_idx == cS.aD_new
                [K_grid, ~] = ndgrid(cS.kGridV, ones(1, 1));  % 只使用普通资本
                pretax_non_capital_income = b_age_val;
                capital_tax = cS.tau_k .* M_age.r_mkt_t .* K_grid;
                total_resources = K_grid .* (1 + M_age.r_mkt_t) + pretax_non_capital_income - capital_tax;

                if cS.phi_bequest > 0
                    % 有遗赠动机：在消费和遗赠之间进行最优权衡（与向量化版本一致）
                    % 最优化问题：max U(c) + φ_bequest × U(bequest)  s.t. c + bequest = total_wealth

                    % 使用解析解求解最优分配（基于相同的CRRA效用函数）
                    if abs(cS.sigma - 1) < 1e-6
                        % 对数效用的情况
                        optimal_c_share = 1 / (1 + cS.phi_bequest);
                    else
                        % CRRA效用的情况
                        optimal_c_share = 1 / (1 + cS.phi_bequest^(1/cS.sigma));
                    end

                    % 计算最优消费和遗赠
                    c_final = optimal_c_share * total_resources / (1 + cS.tau_c);
                    k_prime_final = (1 - optimal_c_share) * total_resources;

                    % 计算效用
                    [~, util_c] = model_setup_utils.CES_utility(c_final, cS.sigma, cS);
                    util_bequest = model_setup_utils.bequest_utility(k_prime_final, cS);
                    util_final = util_c + util_bequest;

                    for ie = 1:cS.nw_expanded
                        cPol_age_q(:,ie) = c_final;
                        kPol_age(:,ie) = k_prime_final; % 记录遗赠额
                        val_age(:,ie) = util_final;
                    end
                else
                    % 无遗赠动机：消费掉所有资源
                    c_expend_final = total_resources;
                    final_c = c_expend_final / (1 + cS.tau_c);
                    [~, util_c] = model_setup_utils.CES_utility(final_c, cS.sigma, cS);
                    for ie = 1:cS.nw_expanded
                        cPol_age_q(:,ie) = final_c;
                        val_age(:,ie) = util_c;
                        kPol_age(:,ie) = 0;
                    end
                end
                cPpsPol_age_choice(:,:) = 0; % 最后一期没有PPS选择
                return;
            end

            effective_discount_factor = (cS.beta ^ cS.time_Step) * cS.s_pathV(a_idx);
            bequest_discount_factor = (cS.beta ^ cS.time_Step) * (1 - cS.s_pathV(a_idx));

            EV_matrix = zeros(cS.nk, cS.nw_expanded);
            if ~isempty(vPrime_kkppse_next)
                transition_matrix_next_age = paramS_age.TrProbM_by_age{a_idx + 1};
                for ie_current = 1:cS.nw_expanded
                    transition_probs = transition_matrix_next_age(ie_current, :);
                    EV_slice = vPrime_kkppse_next * transition_probs';
                    EV_matrix(:, ie_current) = EV_slice;
                end
            end
            EV_interpolants = cell(cS.nw_expanded, 1);
            for ie_current = 1:cS.nw_expanded
                EV_interpolants{ie_current} = griddedInterpolant(cS.kGridV, EV_matrix(:,ie_current), 'pchip', 'none');
            end

            market_return_factor = 1 + M_age.r_mkt_t;
            for ie = 1:cS.nw_expanded
                epsilon_state = paramS_age.leGridV(ie);
                ev_interpolant = EV_interpolants{ie};
                for ik = 1:cS.nk
                    k_state = cS.kGridV(ik);
                    best_val = -1e20; best_k_prime = cS.kMin; best_c = 0;

                    capital_income = k_state * M_age.r_mkt_t;
                    labor_income_gross = 0; pension_benefit = 0;
                    if a_idx <= cS.aR_new
                        labor_income_gross = M_age.w_t*cS.ageEffV_new(a_idx)*epsilon_state;
                    else
                        pension_benefit = b_age_val;
                    end
                    payg_tax = cS.theta_t*labor_income_gross;
                    labor_tax = cS.tau_l*max(0, labor_income_gross - payg_tax);
                    capital_tax = cS.tau_k * capital_income;

                    net_cash_before_shock = k_state * market_return_factor + labor_income_gross + pension_benefit - (payg_tax + labor_tax + capital_tax);

                    shock_expenditure = 0;
                    if ie == cS.nw + 1
                        shock_expenditure = cS.kappa_young * net_cash_before_shock;
                    elseif ie == cS.nw + 2
                        shock_expenditure = cS.kappa_old * net_cash_before_shock;
                    end

                    net_cash_for_c_k_prime = net_cash_before_shock - shock_expenditure;

                    for k_prime_choice = cS.kGridV'
                        c_expend = net_cash_for_c_k_prime - k_prime_choice;

                        if c_expend <= 0
                            break;
                        end

                        c_choice = c_expend / (1 + cS.tau_c);
                        [~,util] = model_setup_utils.CES_utility(c_choice, cS.sigma, cS);

                        ev = ev_interpolant(k_prime_choice);
                        if isnan(ev), ev = -1e10; end

                        % 在当期价值计算中加入遗赠效用
                        util_bequest = model_setup_utils.bequest_utility(k_prime_choice, cS);

                        current_val = util + effective_discount_factor * ev + bequest_discount_factor * util_bequest;

                        if current_val > best_val, best_val = current_val; best_c = c_choice; best_k_prime = k_prime_choice; end
                    end

                    if isinf(best_val)||isnan(best_val)
                        val_age(ik,ie)=-1e20; kPol_age(ik,ie)=cS.kMin; cPpsPol_age_choice(ik,ie)=0; cPol_age_q(ik,ie)=0;
                    else
                        val_age(ik,ie)=best_val; kPol_age(ik,ie)=best_k_prime; cPpsPol_age_choice(ik,ie)=0; cPol_age_q(ik,ie)=best_c;
                    end
                end
            end
        end

        function [cPol_age_q, kPol_age, cPpsPol_age_choice, val_age] = HHSolutionByAge_VFI_ExpenditureShock_Vectorized(a_idx, vPrime_kkppse_next, M_age, b_age_val, paramS_age, cS)
            % [向量化版本 - 初始稳态]
            % [已修正] 统一了最后一期的最优消费/遗赠权衡逻辑

            val_age = -1e20 * ones(cS.nk, cS.nw_expanded);
            cPol_age_q = zeros(cS.nk, cS.nw_expanded);
            kPol_age = zeros(cS.nk, cS.nw_expanded);
            cPpsPol_age_choice = zeros(cS.nk, cS.nw_expanded);

            % --- 最后一期的处理逻辑 [已修正] ---
            if a_idx == cS.aD_new
                [K_grid, ~] = ndgrid(cS.kGridV, ones(1, 1));
                pretax_non_capital_income = b_age_val;
                capital_tax = cS.tau_k .* M_age.r_mkt_t .* K_grid;
                total_resources = K_grid .* (1 + M_age.r_mkt_t) + pretax_non_capital_income - capital_tax;

                if cS.phi_bequest > 0
                    % 有遗赠动机：在消费和遗赠之间进行最优权衡
                    % 最优化问题：max U(c) + φ_bequest × U(bequest)  s.t. c(1+tau_c) + bequest = total_resources

                    % 使用解析解求解最优分配（基于相同的CRRA效用函数）
                    if abs(cS.sigma - 1) < 1e-6
                        % 对数效用的情况
                        optimal_c_expend_share = 1 / (1 + cS.phi_bequest * (1+cS.tau_c));
                    else
                        % CRRA效用的情况
                        optimal_c_expend_share = 1 / (1 + cS.phi_bequest^(1/cS.sigma) * (1+cS.tau_c)^((1-cS.sigma)/cS.sigma) );
                    end

                    % 计算最优消费支出和遗赠
                    c_expend_final = optimal_c_expend_share * total_resources;
                    c_final = c_expend_final / (1 + cS.tau_c);
                    k_prime_final = (1 - optimal_c_expend_share) * total_resources;

                    % 计算效用
                    [~, util_c] = model_setup_utils.CES_utility(c_final, cS.sigma, cS);
                    util_bequest = model_setup_utils.bequest_utility(k_prime_final, cS);
                    util_final = util_c + util_bequest;

                    for ie = 1:cS.nw_expanded
                        cPol_age_q(:,ie) = c_final;
                        kPol_age(:,ie) = k_prime_final;
                        val_age(:,ie) = util_final;
                    end
                else
                    % 无遗赠动机：消费掉所有资源
                    c_expend_final = total_resources;
                    final_c = c_expend_final / (1 + cS.tau_c);
                    [~, util_c] = model_setup_utils.CES_utility(final_c, cS.sigma, cS);
                    for ie = 1:cS.nw_expanded
                        cPol_age_q(:,ie) = final_c;
                        val_age(:,ie) = util_c;
                        kPol_age(:,ie) = 0;
                    end
                end
                cPpsPol_age_choice(:,:) = 0;
                return;
            end

            % --- 向量化计算开始 ---
            effective_discount_factor = (cS.beta ^ cS.time_Step) * cS.s_pathV(a_idx);
            bequest_discount_factor = (cS.beta ^ cS.time_Step) * (1 - cS.s_pathV(a_idx));

            % 计算期望价值矩阵
            EV_matrix = zeros(cS.nk, cS.nw_expanded);
            if ~isempty(vPrime_kkppse_next)
                transition_matrix_next_age = paramS_age.TrProbM_by_age{a_idx + 1};
                for ie_current = 1:cS.nw_expanded
                    transition_probs = transition_matrix_next_age(ie_current, :);
                    EV_slice = vPrime_kkppse_next * transition_probs';
                    EV_matrix(:, ie_current) = EV_slice;
                end
            end

            % 准备期望价值插值器
            EV_interpolants = cell(cS.nw_expanded, 1);
            for ie_current = 1:cS.nw_expanded
                EV_interpolants{ie_current} = griddedInterpolant(cS.kGridV, EV_matrix(:,ie_current), 'pchip', 'none');
            end

            market_return_factor = 1 + M_age.r_mkt_t;

            % --- 核心向量化部分：同时计算所有状态组合 ---
            for ie = 1:cS.nw_expanded
                epsilon_state = paramS_age.leGridV(ie);
                ev_interpolant = EV_interpolants{ie};

                % 为当前收入冲击状态创建向量化的状态网格
                K_states = cS.kGridV; % [nk x 1]

                % 计算收入和税收（向量化）
                capital_income = K_states * M_age.r_mkt_t; % [nk x 1]

                if a_idx <= cS.aR_new
                    labor_income_gross = M_age.w_t * cS.ageEffV_new(a_idx) * epsilon_state;
                    pension_benefit = 0;
                else
                    labor_income_gross = 0;
                    pension_benefit = b_age_val;
                end

                payg_tax = cS.theta_t * labor_income_gross;
                labor_tax = cS.tau_l * max(0, labor_income_gross - payg_tax);
                capital_tax = cS.tau_k * capital_income; % [nk x 1]

                net_cash_before_shock = K_states * market_return_factor + labor_income_gross + pension_benefit - (payg_tax + labor_tax + capital_tax); % [nk x 1]

                % 计算冲击支出
                shock_expenditure = 0;
                if ie == cS.nw + 1
                    shock_expenditure = cS.kappa_young * net_cash_before_shock;
                elseif ie == cS.nw + 2
                    shock_expenditure = cS.kappa_old * net_cash_before_shock;
                end

                net_cash_for_c_k_prime = net_cash_before_shock - shock_expenditure; % [nk x 1]

                % 对于每个资本状态，计算最优选择
                for ik = 1:cS.nk
                    available_cash = net_cash_for_c_k_prime(ik);

                    % 创建所有可能的储蓄选择向量
                    k_prime_choices = cS.kGridV'; % [1 x nk_choices]

                    % 计算对应的消费支出
                    c_expend_choices = available_cash - k_prime_choices; % [1 x nk_choices]

                    % 找到可行的选择（消费支出 > 0）
                    feasible_mask = c_expend_choices > 0; % [1 x nk_choices]

                    if ~any(feasible_mask)
                        % 没有可行选择，设置为最低值
                        val_age(ik, ie) = -1e20;
                        kPol_age(ik, ie) = cS.kMin;
                        cPpsPol_age_choice(ik, ie) = 0;
                        cPol_age_q(ik, ie) = 0;
                        continue;
                    end

                    % 只保留可行的选择
                    feasible_k_primes = k_prime_choices(feasible_mask); % [1 x n_feasible]
                    feasible_c_expends = c_expend_choices(feasible_mask); % [1 x n_feasible]

                    % 向量化计算所有可行选择的消费和效用
                    feasible_c_choices = feasible_c_expends / (1 + cS.tau_c); % [1 x n_feasible]

                    % 计算消费效用（向量化）
                    [~, util_c_vec] = model_setup_utils.CES_utility(feasible_c_choices, cS.sigma, cS);

                    % 计算期望价值（向量化插值）
                    ev_vec = ev_interpolant(feasible_k_primes'); % [n_feasible x 1]
                    ev_vec = ev_vec'; % [1 x n_feasible]

                    % 处理NaN值
                    ev_vec(isnan(ev_vec)) = -1e10;

                    % 计算遗赠效用（向量化）
                    util_bequest_vec = model_setup_utils.bequest_utility(feasible_k_primes, cS); % [1 x n_feasible]

                    % 计算总效用
                    total_values = util_c_vec + effective_discount_factor * ev_vec + bequest_discount_factor * util_bequest_vec; % [1 x n_feasible]

                    % 找到最优选择
                    [best_val, best_idx] = max(total_values);

                    if isfinite(best_val)
                        val_age(ik, ie) = best_val;
                        kPol_age(ik, ie) = feasible_k_primes(best_idx);
                        cPol_age_q(ik, ie) = feasible_c_choices(best_idx);
                        cPpsPol_age_choice(ik, ie) = 0;
                    else
                        val_age(ik, ie) = -1e20;
                        kPol_age(ik, ie) = cS.kMin;
                        cPpsPol_age_choice(ik, ie) = 0;
                        cPol_age_q(ik, ie) = 0;
                    end
                end
            end
        end

        function k_prime_idx=get_policy_index_matrix(kPolM,cS)
            % [vPaper.4.1 - 已修正离散化规则]
            % 修正: 更改了离散化方法。不再寻找"最近"的网格点，
            %       而是寻找"可负担的最高"网格点，以确保预算约束在
            %       离散化后永不被违反，从而避免 c_expend < 0。
            k_prime_idx=zeros(cS.nk,cS.nw_expanded,cS.aD_new,'uint16');
            for ia=1:cS.aD_new
                for ie=1:cS.nw_expanded
                    for ik=1:cS.nk
                        % 获取连续空间下的最优储蓄决策
                        k_prime_continuous = kPolM(ik,ie,ia);

                        % 寻找所有小于等于该最优决策的网格点的索引
                        affordable_indices = find(cS.kGridV <= k_prime_continuous);

                        if isempty(affordable_indices)
                            % 如果没有任何网格点可负担（例如k_prime_continuous为负）
                            % 则强制选择最低的资产水平（通常是k=0，索引为1）
                            idx = 1;
                        else
                            % 在所有可负担的网格点中，选择最高的那个
                            idx = affordable_indices(end);
                        end
                        k_prime_idx(ik,ie,ia) = idx;
                    end
                end
            end
        end

        % --- 校准模块的目标函数 (嵌套函数) ---
        function error_sq = calibration_objective(x, cS_base, paramS_calib, params_ext_calib, target_ky, lb, ub)
            % [核心修改] 增加了 lb 和 ub 输入，并加入了边界惩罚逻辑
            % 目标函数: 返回 (模型K/Y - 目标K/Y)^2 + 边界惩罚
            % 输入:
            %   x: 待校准参数向量 [beta, gamma, lambda_g]
            %   cS_base: 基础参数结构体
            %   paramS_calib: 冲击过程参数
            %   params_ext_calib: 校准年份的外生变量 (Z, A, theta)
            %   target_ky: 目标K/Y比率
            %   lb: 参数下界向量
            %   ub: 参数上界向量

            % 1. 检查参数是否越界
            if any(x < lb) || any(x > ub)
                % 如果任何一个参数越界，返回一个巨大的惩罚值
                violation = max(0, lb - x) + max(0, x - ub); % 计算越界距离
                penalty = 1e8 + 1e6 * sum(violation.^2); % 与越界距离成正比的惩罚
                error_sq = penalty;
                fprintf('   参数 [beta=%.4f, gamma=%.3f, lambda_g=%.3f] 超出边界，惩罚值=%.3e\n', x(1), x(2), x(3), penalty);
                return;
            end

            % 2. 额外的参数合理性检查
            if any(x <= 0) || x(1) >= 2 || x(2) >= 1 || x(3) >= 1
                error_sq = 1e10;
                fprintf('   参数 [beta=%.4f, gamma=%.3f, lambda_g=%.3f] 超出合理范围，返回巨大误差\n', x(1), x(2), x(3));
                return;
            end

            % 3. 创建一个临时的cS副本并更新参数
            cS_temp = cS_base;
            cS_temp.beta = x(1);
            cS_temp.gamma = x(2);
            cS_temp.lambda_g = x(3);

            % 4. 使用 try-catch 包装稳态求解器
            try
                [ss_calib, ~, ~, ~] = main_steady_state_utils.solve_steady_state_complete(cS_temp, paramS_calib, params_ext_calib, false);

                % 5. 检查求解是否成功以及结果是否合理
                if isempty(ss_calib) || ~isfield(ss_calib, 'K_private') || ss_calib.Y_from_production < 1e-6
                    error_sq = 1e9; % 求解失败，返回大误差但小于边界惩罚
                    fprintf('   尝试参数 [beta=%.4f, gamma=%.3f, lambda_g=%.3f] -> 求解失败\n', x(1), x(2), x(3));
                    return;
                end

                % 6. 额外的合理性检查
                if ss_calib.K_private < 0 || ss_calib.K_public < 0 || ss_calib.L < 0
                    error_sq = 1e9;
                    fprintf('   尝试参数 [beta=%.4f, gamma=%.3f, lambda_g=%.3f] -> 负值结果\n', x(1), x(2), x(3));
                    return;
                end

                % 7. 计算模型生成的K/Y比率
                model_ky = (ss_calib.K_private + ss_calib.K_public) / ss_calib.Y_from_production;

                % 8. 检查K/Y是否合理
                if model_ky < 0 || model_ky > 20 || isnan(model_ky) || isinf(model_ky)
                    error_sq = 1e9;
                    fprintf('   尝试参数 [beta=%.4f, gamma=%.3f, lambda_g=%.3f] -> K/Y异常=%.4f\n', x(1), x(2), x(3), model_ky);
                    return;
                end

                % 9. 计算平方误差
                error_sq = (model_ky - target_ky)^2;

                fprintf('   尝试参数 [beta=%.4f, gamma=%.3f, lambda_g=%.3f] -> K/Y=%.4f, 误差^2=%.3e\n', ...
                    x(1), x(2), x(3), model_ky, error_sq);

            catch ME
                fprintf('   尝试参数 [beta=%.4f, gamma=%.3f, lambda_g=%.3f] -> 异常: %s\n', x(1), x(2), x(3), ME.message);
                error_sq = 1e9; % 发生错误返回大误差但小于边界惩罚
            end
        end

        function dist_next = evolve_distribution_one_step(dist_now, k_prime_idx_now, Z_path_t, Z_path_t_plus_1, cS, paramS)
            % [v_trans - 修正版]
            % 核心修改:
            % 1. 明确函数输入为当前和下一期的人口路径 Z_path，使逻辑更清晰。
            % 2. 使用 Z_path(ia+1, t+1) / Z_path(ia, t) 作为隐含存活率来演化分布，
            %    确保与外部人口数据完全一致。这与稳态求解的逻辑是一致的。

            dist_next = zeros(cS.nk, cS.nw_expanded, cS.aD_new);

            % 1. 新生儿分布 (在 t+1 期)
            dist_newborn_ke = zeros(cS.nk, cS.nw_expanded);
            dist_newborn_ke(1, 1:cS.nw) = paramS.leProb1V';
            dist_next(:, :, 1) = dist_newborn_ke * Z_path_t_plus_1(1);

            % 2. 幸存者演化 (从 t 到 t+1)
            for ia = 1:(cS.aD_new - 1)
                dist_ia_ke = dist_now(:,:,ia); % t 期的分布
                dist_ia_plus_1_ke_unscaled = zeros(cS.nk, cS.nw_expanded);
                transition_matrix_next_age = paramS.TrProbM_by_age{ia + 1};

                % 计算未缩放的下一期分布 (假设存活率为1)
                for ik = 1:cS.nk
                    for ie = 1:cS.nw_expanded
                        mass_at_state = dist_ia_ke(ik, ie);
                        if mass_at_state < 1e-20, continue; end

                        ik_prime = k_prime_idx_now(ik, 1, ie, ia); % 假设 nkpps=1
                        transition_probs_e = transition_matrix_next_age(ie, :);

                        dist_ia_plus_1_ke_unscaled(ik_prime, :) = dist_ia_plus_1_ke_unscaled(ik_prime, :) + mass_at_state * transition_probs_e;
                    end
                end

                % 使用 Z_path 来重新缩放，确保与外部数据吻合
                mass_at_ia_t = Z_path_t(ia); % t 期的总人口
                if mass_at_ia_t > 1e-12
                    rescale_factor = Z_path_t_plus_1(ia+1) / mass_at_ia_t;
                    dist_next(:,:,ia+1) = dist_ia_plus_1_ke_unscaled * rescale_factor;
                else
                    dist_next(:,:,ia+1) = zeros(cS.nk, cS.nw_expanded);
                end
            end

            % 归一化总人口 (作为安全检查)
            total_mass_next = sum(dist_next, 'all');
            if abs(total_mass_next - 1.0) > 1e-4
                dist_next = dist_next / total_mass_next;
            end
        end

        function K_agg = aggregate_k_prime(dist_now, k_prime_idx_now, prob_survive_implied, cS)
            % [修改] 函数现在显式接收隐含存活率向量作为参数
            K_agg = 0;
            total_mass = sum(dist_now, 'all');
            if total_mass < 1e-9, return; end
            dist_norm = dist_now / total_mass;

            for ia = 1:cS.aD_new
                % [核心修改] 使用传入的、与宏观演化一致的隐含存活率
                prob_survive = prob_survive_implied(ia);

                for ie = 1:cS.nw_expanded
                    for ik = 1:cS.nk
                        mass = dist_norm(ik, ie, ia);
                        if mass < 1e-20, continue; end
                        idx_k_prime = k_prime_idx_now(ik, ie, ia);
                        k_prime = cS.kGridV(idx_k_prime);

                        % 下一期的资本存量是当期储蓄 k' 乘以存活下来的人口
                        K_agg = K_agg + k_prime * mass * prob_survive;
                    end
                end
            end
        end

        function [L_agg, Tax_agg, Bequest_tax_agg] = aggregate_flow_vars_at_t(dist_t, k_prime_idx_t, M_sim_t, prob_survive_implied, cS, paramS)
            % [修改] 函数现在显式接收隐含存活率向量作为参数
            L_agg=0; Tax_agg=0; Bequest_tax_agg=0;
            total_mass = sum(dist_t, 'all');
            if total_mass < 1e-9, return; end
            dist_norm = dist_t / total_mass;

            for ia = 1:cS.aD_new
                for ie = 1:cS.nw_expanded
                    for ik = 1:cS.nk
                        mass = dist_norm(ik,ie,ia);
                        if mass < 1e-20, continue; end

                        k_now = cS.kGridV(ik);
                        epsilon_val = paramS.leGridV(ie);
                        idx_k_prime = k_prime_idx_t(ik, ie, ia);
                        k_prime = cS.kGridV(idx_k_prime);

                        % 使用已有的 backout 函数
                        [~,tax_val,~,~,~,~,~,~] = main_steady_state_utils.backout_accounting_expenditure_shock(k_now, k_prime, ia, ie, epsilon_val, M_sim_t, cS);

                        Tax_agg = Tax_agg + tax_val * mass;

                        % [核心修改] 使用传入的隐含死亡率
                        prob_death = 1 - prob_survive_implied(ia);
                        Bequest_tax_agg = Bequest_tax_agg + k_prime * mass * prob_death;

                        if ia <= cS.aR_new
                            L_agg = L_agg + (cS.ageEffV_new(ia)*epsilon_val) * mass;
                        end
                    end
                end
            end
        end

        function Dist=solve_steady_state_distribution(k_prime_idx,paramS,cS,Z_ss_norm)
            % [v15 - 修正版]
            % 核心修正: 修复了当使用外部人口数据(Z_ss_norm)时，
            %           由于同时使用内部存活率(s_pathV)导致的会计不一致问题。
            % 解决方案: 不再使用 s_pathV。而是通过 Z_ss_norm(ia+1)/Z_ss_norm(ia)
            %           来反推隐含的存活率，以保证最终分布与外部数据完全吻合。

            Dist = zeros(cS.nk, cS.nw_expanded, cS.aD_new);

            % 1. 初始化新生儿分布
            % 新生儿只分布在k=0，e为初始冲击分布的状态
            dist_newborn_ke = zeros(cS.nk, cS.nw_expanded);
            dist_newborn_ke(1, 1:cS.nw) = paramS.leProb1V';

            % 新生儿的总人口由 Z_ss_norm(1) 决定
            Dist(:, :, 1) = dist_newborn_ke * Z_ss_norm(1);

            % 2. 迭代计算幸存者的分布
            for ia = 1:(cS.aD_new - 1)
                dist_ia_ke = Dist(:,:,ia); % 当前年龄的分布
                dist_ia_plus_1_ke_unscaled = zeros(cS.nk, cS.nw_expanded);
                transition_matrix_next_age = paramS.TrProbM_by_age{ia + 1};

                % 计算未缩放的下一期分布 (假设存活率为1)
                for ik = 1:cS.nk
                    for ie = 1:cS.nw_expanded
                        mass_at_state = dist_ia_ke(ik, ie);
                        if mass_at_state < 1e-20, continue; end

                        ik_prime = k_prime_idx(ik, ie, ia);
                        transition_probs_e = transition_matrix_next_age(ie, :);

                        % [核心修改] 不再乘以 s_pathV(ia)
                        % 只是将当前质量按照决策和冲击转移到下一期
                        dist_ia_plus_1_ke_unscaled(ik_prime, :) = dist_ia_plus_1_ke_unscaled(ik_prime, :) + mass_at_state * transition_probs_e;
                    end
                end

                % [核心修改] 使用 Z_ss_norm 来重新缩放下一期的分布
                mass_at_ia = sum(dist_ia_ke, 'all'); % 这应该等于 Z_ss_norm(ia)

                if mass_at_ia > 1e-12
                    % 计算隐含的存活率/缩放因子
                    rescale_factor = Z_ss_norm(ia+1) / mass_at_ia;
                    % 应用缩放因子，确保下一期总人口等于 Z_ss_norm(ia+1)
                    Dist(:,:,ia+1) = dist_ia_plus_1_ke_unscaled * rescale_factor;
                else
                    % 如果当前年龄组没人，下一年龄组也没人
                    Dist(:,:,ia+1) = zeros(cS.nk, cS.nw_expanded);
                end
            end

            % 3. 最终检查
            final_sum = sum(Dist, 'all');
            if abs(final_sum - 1.0) > 1e-6
                % 这个警告现在理论上不应该再被触发
                warning('稳态联合分布总和(%.8f)不为1，可能存在会计不一致。', final_sum);
            end
        end

        function [K_agg,C_utility_agg,Tax_agg,Bequest_tax_agg,L_agg,PensionIn_agg,PensionOut_agg,ShockExp_agg] = aggregate_expenditure_shock(Dist, k_prime_idx, M_sim, cS, paramS)
            % [修改] 函数根据上下文使用正确的存活率：稳态时使用prob_survive_implied_ss0

            if ~isfield(cS, 'theta_t'), cS.theta_t = cS.theta_path(1); end
            K_agg=0; C_utility_agg=0; Tax_agg=0; Bequest_tax_agg=0; L_agg=0; PensionIn_agg=0; PensionOut_agg=0; ShockExp_agg=0;

            % [新增] 根据上下文确定使用哪个存活率向量
            if isfield(cS, 'prob_survive_implied_ss0') && ~isempty(cS.prob_survive_implied_ss0)
                % 稳态求解：使用基于稳态人口分布的隐含存活率
                prob_survive_implied = cS.prob_survive_implied_ss0;
            elseif isfield(cS, 'prob_survive_implied_trans') && ~isempty(cS.prob_survive_implied_trans)
                % 过渡态求解：使用基于当期人口分布的隐含存活率
                prob_survive_implied = cS.prob_survive_implied_trans;
            else
                % 备用方案：使用外生存活率（但这可能导致会计不一致）
                warning('未找到隐含存活率，使用外生存活率 s_pathV');
                prob_survive_implied = cS.s_pathV;
            end

            for ia=1:cS.aD_new
                for ie=1:cS.nw_expanded
                    for ik=1:cS.nk
                        mass = Dist(ik,ie,ia);
                        if mass < 1e-20, continue; end
                        k_now=cS.kGridV(ik); epsilon_val=paramS.leGridV(ie);
                        idx_k_prime=k_prime_idx(ik,ie,ia); k_prime=cS.kGridV(idx_k_prime);

                        [c_val,tax_val,shock_exp,payg_tax,~,~,~,pension_benefit] = ...
                            main_steady_state_utils.backout_accounting_expenditure_shock(k_now,k_prime,ia,ie,epsilon_val,M_sim,cS);

                        C_utility_agg = C_utility_agg + c_val * mass;
                        Tax_agg = Tax_agg + tax_val * mass;
                        ShockExp_agg = ShockExp_agg + shock_exp * mass;

                        % [核心修改] 使用上下文相关的隐含存活/死亡率
                        prob_survive = prob_survive_implied(ia);
                        K_agg = K_agg + k_prime * mass * prob_survive;
                        prob_death = 1 - prob_survive;
                        Bequest_tax_agg = Bequest_tax_agg + k_prime * mass * prob_death;

                        if ia<=cS.aR_new, L_agg = L_agg + (cS.ageEffV_new(ia)*epsilon_val) * mass; end
                        PensionIn_agg = PensionIn_agg + payg_tax * mass;
                        PensionOut_agg = PensionOut_agg + pension_benefit * mass;
                    end
                end
            end
        end

        function [c_val,tax_val,shock_expenditure,payg_tax,labor_tax,capital_tax,consumption_tax,pension_benefit]=backout_accounting_expenditure_shock(k_now,k_prime,ia,ie,epsilon_val,M_sim,cS)
            % [vPaper.4.8 - 弹性消费下限版]
            % 核心变更: 移除了 c_val = max(cFloor, ...) 的强制修正。
            %           现在函数直接报告由预算约束反解出的实际消费值，
            %           即使它可能低于cFloor。

            if ~isfield(cS, 'theta_t'), cS.theta_t = cS.theta_path(1); end

            market_return_factor = 1 + M_sim.r_mkt_t;
            capital_income = k_now * M_sim.r_mkt_t;
            labor_income_gross = 0; pension_benefit = 0;
            if ia <= cS.aR_new
                labor_income_gross = M_sim.w_t * cS.ageEffV_new(ia) * epsilon_val;
            else
                pension_benefit = M_sim.b_t;
            end
            total_inflow = k_now * market_return_factor + labor_income_gross + pension_benefit;

            payg_tax = cS.theta_t * labor_income_gross;
            capital_tax = cS.tau_k * capital_income;
            labor_tax = cS.tau_l * max(0, labor_income_gross - payg_tax);
            net_cash_before_shock = total_inflow - (payg_tax + capital_tax + labor_tax);

            shock_expenditure = 0;
            if ie == cS.nw + 1
                shock_expenditure = cS.kappa_young * net_cash_before_shock;
            elseif ie == cS.nw + 2
                shock_expenditure = cS.kappa_old * net_cash_before_shock;
            end

            net_cash_for_c_k_prime = net_cash_before_shock - shock_expenditure;

            % [!!!!! 核心修改 !!!!!]
            c_expend_available = net_cash_for_c_k_prime - k_prime;

            % 我们不再用max(cFloor,...)来强制修正消费。
            % 我们相信上游的VFI和离散化步骤会保证c_expend_available是正的。
            c_val = c_expend_available / (1 + cS.tau_c);

            % 如果c_val最终为负或零，说明上游有更严重的问题，让它在后续计算中暴露出来。

            consumption_tax = c_val * cS.tau_c;
            tax_val = labor_tax + capital_tax + consumption_tax;

            % --- 预算平衡检验依然保留，作为最终防线 ---
            total_consumption_expenditure = c_val * (1 + cS.tau_c);
            total_tax_outflow = payg_tax + capital_tax + labor_tax;
            total_outflow = k_prime + total_consumption_expenditure + shock_expenditure + total_tax_outflow;
            budget_gap = total_inflow - total_outflow;

            if abs(budget_gap) > 1e-9
                % 现在，预算缺口应该非常接近于0，因为我们没有"丢弃"任何现金
                error('微观预算约束被违反！状态: (ia=%d, ie=%d, ik=%d), 预算缺口: %.3e', ...
                    ia, ie, find(abs(cS.kGridV - k_now) < 1e-10, 1), budget_gap);
            end
        end

        function M_prices = get_prices_at_t(K_p, K_g, L, A_t, cS)
            % [vPaper.5.1 - 政府投资修正版]
            % 核心修改: 将公共资本 K_g 作为外部性，修正生产函数，确保国民账户平衡。
            %           Y = (A * K_g^gamma) * K_p^alpha * L^(1-alpha)

            if K_p <= 0, K_p = 1e-8; end; if L <= 0, L = 1e-8; end; if K_g <= 0, K_g = 1e-8; end;

            % 有效的总要素生产率，受公共资本影响
            A_effective = A_t .* (K_g.^cS.gamma);

            % 生产函数，现在对私人要素是常数规模报酬
            Y_period = A_effective .* (K_p.^cS.alpha) .* (L.^(1-cS.alpha));

            % 私人资本的边际产出
            MPK_p_period = cS.alpha .* Y_period ./ K_p;

            % 工资 (现在占总产出的 1-alpha 份额)
            w_t = (1-cS.alpha) .* Y_period ./ L;

            % 私人资本的市场回报率
            r_mkt_t = MPK_p_period - cS.ddk;

            M_prices = struct('K_p', K_p, 'K_g', K_g, 'L_t', L, 'Y_t', Y_period, 'w_t', w_t, 'r_mkt_t', r_mkt_t);
        end


        % =======================================================
        % == 阶段二：稳态求解器 (作为过渡态的终点，包含PPS决策和额外的PPS状态变量)
        % =======================================================


        % =======================================================
        % == 结果分析与展示
        % =======================================================

        function is_ok = check_national_accounts(ss, Dist, k_prime_idx, cS, paramS)
            % [新增] 用于校准循环的简化版国民账户检查
            M_sim = struct('r_mkt_t', ss.r_mkt, 'w_t', ss.w, 'b_t', ss.b);

            [~, C_utility_agg, Tax_agg, Bequest_tax_agg, ~, ~, ~, ShockExp_agg] = ...
                main_steady_state_utils.aggregate_expenditure_shock(Dist, k_prime_idx, M_sim, cS, paramS);

            Y_prod = ss.Y_from_production;
            Gov_Revenue_total = Tax_agg + Bequest_tax_agg;
            I_g_agg = cS.lambda_g * Gov_Revenue_total;
            G_c_agg = (1 - cS.lambda_g) * Gov_Revenue_total;
            I_p_agg_gross = cS.ddk * ss.K_private;
            I_total_agg = I_p_agg_gross + I_g_agg;
            C_total_agg = C_utility_agg + ShockExp_agg;
            Y_exp_actual = C_total_agg + I_total_agg + G_c_agg;

            error_gdp = abs(Y_exp_actual - Y_prod);
            if error_gdp < 1e-4
                fprintf('   国民账户核算: 通过 (误差: %.3e)\n', error_gdp);
                is_ok = true;
            else
                fprintf('   国民账户核算: 失败! (误差: %.3e)\n', error_gdp);
                is_ok = false;
            end
        end

        function display_national_accounts_gov_investment(ss, Dist, k_prime_idx, cS, paramS)
            fprintf('\n\n========================================================================\n');
            fprintf('===     国民经济核算详细报告 (政府投资模型版)     ===\n');
            fprintf('========================================================================\n');

            M_sim = struct('r_mkt_t', ss.r_mkt, 'w_t', ss.w, 'b_t', ss.b);


            % [修改] aggregate_expenditure_shock函数会自动选择正确的存活率向量
            [~, C_utility_agg, Tax_agg, Bequest_tax_agg, L_agg_check, PensionIn_agg, PensionOut_agg, ShockExp_agg] = ...
                main_steady_state_utils.aggregate_expenditure_shock(Dist, k_prime_idx, M_sim, cS, paramS);

            % --- 宏观变量定义 ---
            Y_prod = ss.Y_from_production;
            K_p = ss.K_private;
            K_g = ss.K_public;

            Gov_Revenue_total = Tax_agg + Bequest_tax_agg;
            I_g_agg = cS.lambda_g * Gov_Revenue_total;         % 政府投资
            G_c_agg = (1 - cS.lambda_g) * Gov_Revenue_total; % 政府消费

            I_p_agg_gross = cS.ddk * K_p; % 私人部门总投资 = 折旧
            I_total_agg = I_p_agg_gross + I_g_agg; % 经济总投资

            C_total_agg = C_utility_agg + ShockExp_agg; % 总消费

            % 支出法GDP
            Y_exp_actual = C_total_agg + I_total_agg + G_c_agg;

            fprintf('--- A. 宏观产出与支出核算 ---\n');
            fprintf('   生产法 GDP (Y_prod):         %.6f\n', Y_prod);
            fprintf('   支出法 GDP (C+I_total+G_c):  %.6f\n', Y_exp_actual);
            fprintf('   ------------------------------------\n');
            fprintf('   核算误差 (Y_exp - Y_prod):     %.3e (此值应接近0)\n', Y_exp_actual - Y_prod);
            fprintf('   总消费 (C):                  %.6f (占GDP: %.2f%%)\n', C_total_agg, C_total_agg/Y_prod*100);
            fprintf('     - 常规消费 (有作用):       %.6f\n', C_utility_agg);
            fprintf('     - 重大冲击支出 (无作用):   %.6f\n', ShockExp_agg);
            fprintf('   总投资 (I_total = I_p+I_g):  %.6f (占GDP: %.2f%%)\n', I_total_agg, I_total_agg/Y_prod*100);
            fprintf('     - 私人投资 (I_p):          %.6f\n', I_p_agg_gross);
            fprintf('     - 政府投资 (I_g):          %.6f\n', I_g_agg);
            fprintf('   政府消费 (G_c):              %.6f (占GDP: %.2f%%)\n', G_c_agg, G_c_agg/Y_prod*100);

            fprintf('\n--- B. 政府与养老金体系核算 ---\n');
            fprintf('   政府总收入 (T):              %.6f\n', Gov_Revenue_total);
            fprintf('   政府总支出 (G_c+I_g):        %.6f\n', G_c_agg + I_g_agg);
            fprintf('   政府预算平衡 (T - G_c - I_g):%.3e (此值应接近0)\n', Gov_Revenue_total - (G_c_agg + I_g_agg));
            fprintf('   ------------------------------------\n');
            fprintf('   养老金体系平衡 (收入-支出):  %.3e (此值应接近0)\n', PensionIn_agg - PensionOut_agg);
            fprintf('   劳动供给核算误差 (L_agg-L_ss): %.3e\n', L_agg_check - ss.L);

            fprintf('\n--- C. 关键宏观比率 ---\n');
            fprintf('   私人资本产出比 (K_p/Y):      %.4f\n', K_p / Y_prod);
            fprintf('   公共资本产出比 (K_g/Y):      %.4f\n', K_g / Y_prod);
            fprintf('   总资本产出比 (K_total/Y):    %.4f\n', (K_p + K_g) / Y_prod);
            fprintf('\n========================================================================\n');
        end

        function plot_saving_rate_by_age(ss, Dist, k_prime_idx, cS, paramS)
            % =========================================================================
            % == FUNCTION: plot_saving_rate_by_age
            % == 目的: 计算并绘制出模型中各个年龄组的平均储蓄率。
            % == 储蓄率定义: S_rate = (Aggregate k_prime) / (Aggregate Net Cash for C/k')
            % =========================================================================

            fprintf('\n   正在计算各年龄组储蓄率...\n');

            % 初始化用于聚合的向量
            age_saving_agg = zeros(cS.aD_new, 1);
            age_disposable_resources_agg = zeros(cS.aD_new, 1);

            M_sim = struct('r_mkt_t', ss.r_mkt, 'w_t', ss.w, 'b_t', ss.b);

            % 遍历所有状态来计算每个年龄组的总储蓄和总可支配资源
            for ia = 1:cS.aD_new
                for ik = 1:cS.nk
                    for ie = 1:cS.nw_expanded
                        mass = Dist(ik, ie, ia);
                        if mass < 1e-20
                            continue;
                        end

                        % 获取当前状态和决策
                        k_now = cS.kGridV(ik);
                        epsilon_val = paramS.leGridV(ie);
                        idx_k_prime = k_prime_idx(ik, ie, ia);
                        k_prime = cS.kGridV(idx_k_prime);

                        % --- 为该代理人计算可用于消费和储蓄的净现金流 ---
                        % 这部分逻辑与 backout_accounting_... 函数中的预算约束计算一致
                        capital_income = k_now * M_sim.r_mkt_t;
                        labor_income_gross = 0;
                        pension_benefit = 0;
                        if ia <= cS.aR_new
                            labor_income_gross = M_sim.w_t * cS.ageEffV_new(ia) * epsilon_val;
                        else
                            pension_benefit = M_sim.b_t;
                        end

                        payg_tax = cS.theta_path(1) * labor_income_gross;
                        labor_tax = cS.tau_l * max(0, labor_income_gross - payg_tax);
                        capital_tax = cS.tau_k * capital_income;

                        net_cash_before_shock = k_now * (1 + M_sim.r_mkt_t) + labor_income_gross + pension_benefit - (payg_tax + labor_tax + capital_tax);

                        shock_expenditure = 0;
                        if ie == cS.nw + 1
                            shock_expenditure = cS.kappa_young * net_cash_before_shock;
                        elseif ie == cS.nw + 2
                            shock_expenditure = cS.kappa_old * net_cash_before_shock;
                        end

                        net_cash_for_c_k_prime = net_cash_before_shock - shock_expenditure;

                        % 聚合该代理人的储蓄和可支配资源
                        age_saving_agg(ia) = age_saving_agg(ia) + mass * k_prime;
                        age_disposable_resources_agg(ia) = age_disposable_resources_agg(ia) + mass * net_cash_for_c_k_prime;
                    end
                end
            end

            % 计算最终的储蓄率
            saving_rate_by_age = age_saving_agg ./ age_disposable_resources_agg;

            % 处理可能因分母为0导致的NaN（例如最后一期）
            saving_rate_by_age(isnan(saving_rate_by_age)) = 0;

            % --- 绘制图形 ---
            % 创建用于绘图的年龄向量（取年龄区间的中间点）
            age_vector_plot = cS.age1_orig + ((1:cS.aD_new) - 1) * cS.time_Step + floor(cS.time_Step / 2);

            figure('Name', '年龄-储蓄率剖面');
            plot(age_vector_plot, saving_rate_by_age * 100, 'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
            hold on;

            % 绘制一条垂线以标记退休年龄
            y_limits = ylim;
            plot([cS.ageRetire_orig, cS.ageRetire_orig], y_limits, 'k--', 'LineWidth', 1.5, 'DisplayName', sprintf('退休年龄 (%d岁)', cS.ageRetire_orig));

            hold off;
            grid on;
            title('各年龄组的平均储蓄率');
            xlabel('年龄 (岁)');
            ylabel('储蓄率 (%)');
            legend('储蓄率', '退休年龄', 'Location', 'best');
            xlim([cS.age1_orig, cS.ageLast_orig]);
        end



        % =======================================================
        % == 阶段三：终点稳态求解器 (包含PPS决策的扩展版本)
        % =======================================================
        function [cPol_age_q, kPol_age, cPpsPol_age_choice, val_age] = HHSolutionByAge_VFI_ExpenditureShock_Vectorized_with_pps_ParFor(a_idx, vPrime_kkppse_next, M_age, b_age_val, paramS_age, cS)
            % [parfor并行版本] 在关键循环中使用并行计算的PPS家庭决策函数
            % 核心改进：将 (ik, ikpps) 双重循环合并为单个parfor循环

            % 抑制MESHGRID插值性能警告
            warning('off', 'MATLAB:griddedInterpolant:MeshgridEval2DWarnId');
            warning('off', 'MATLAB:griddedInterpolant:MeshgridEval2DWarn');
            warning('off', 'MATLAB:griddedInterpolant:MeshgridEval2D');

            val_age = -1e20 * ones(cS.nk, cS.nkpps, cS.nw_expanded);
            cPol_age_q = zeros(cS.nk, cS.nkpps, cS.nw_expanded);
            kPol_age = zeros(cS.nk, cS.nkpps, cS.nw_expanded);
            cPpsPol_age_choice = zeros(cS.nk, cS.nkpps, cS.nw_expanded);

            % --- 最后一期的处理逻辑 ---
            if a_idx == cS.aD_new
                [K_grid, Kpps_grid] = ndgrid(cS.kGridV, cS.kppsGridV);

                % 计算所有资产的税后总财富
                k_after_return = K_grid .* (1 + M_age.r_mkt_t);
                k_capital_tax = cS.tau_k .* (K_grid .* M_age.r_mkt_t);
                k_after_tax_value = k_after_return - k_capital_tax;

                kpps_after_return = Kpps_grid .* (1 + M_age.r_mkt_t);
                kpps_withdrawal_tax = cS.pps_tax_rate_withdrawal .* kpps_after_return;
                kpps_after_tax_value = kpps_after_return - kpps_withdrawal_tax;

                pension_after_tax = b_age_val;
                total_after_tax_wealth = k_after_tax_value + kpps_after_tax_value + pension_after_tax;

                if cS.phi_bequest > 0
                    % 有遗赠动机
                    if abs(cS.sigma - 1) < 1e-6
                        optimal_c_share = 1 / (1 + cS.phi_bequest);
                    else
                        optimal_c_share = 1 / (1 + cS.phi_bequest^(1/cS.sigma));
                    end

                    c_expend_final = optimal_c_share * total_after_tax_wealth;
                    c_final = c_expend_final / (1 + cS.tau_c);
                    optimal_bequest = (1 - optimal_c_share) * total_after_tax_wealth;

                    [~, util_c] = model_setup_utils.CES_utility(c_final, cS.sigma, cS);
                    util_bequest = model_setup_utils.bequest_utility(optimal_bequest, cS);
                    util_final = util_c + util_bequest;

                    k_prime_final = optimal_bequest;
                    kpps_prime_final = zeros(size(Kpps_grid));

                    for ie = 1:cS.nw_expanded
                        cPol_age_q(:,:,ie) = c_final;
                        kPol_age(:,:,ie) = k_prime_final;
                        cPpsPol_age_choice(:,:,ie) = kpps_prime_final;
                        val_age(:,:,ie) = util_final;
                    end
                else
                    % 无遗赠动机
                    c_expend_final = total_after_tax_wealth;
                    final_c = c_expend_final / (1 + cS.tau_c);
                    [~, util_c] = model_setup_utils.CES_utility(final_c, cS.sigma, cS);

                    for ie = 1:cS.nw_expanded
                        cPol_age_q(:,:,ie) = final_c;
                        val_age(:,:,ie) = util_c;
                        kPol_age(:,:,ie) = 0;
                        cPpsPol_age_choice(:,:,ie) = 0;
                    end
                end
                return;
            end

            % --- 并行化计算开始 ---
            effective_discount_factor = (cS.beta ^ cS.time_Step) * cS.s_pathV(a_idx);
            bequest_discount_factor = (cS.beta ^ cS.time_Step) * (1 - cS.s_pathV(a_idx));

            % 计算期望价值矩阵
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

            % 准备期望价值插值器
            EV_interpolants = cell(cS.nw_expanded, 1);
            for ie_current = 1:cS.nw_expanded
                EV_interpolants{ie_current} = griddedInterpolant({cS.kGridV, cS.kppsGridV}, EV_matrix(:,:,ie_current), 'linear', 'none');
            end

            market_return_factor = 1 + M_age.r_mkt_t;

            % --- 核心parfor并行化部分 ---
            for ie = 1:cS.nw_expanded
                epsilon_state = paramS_age.leGridV(ie);
                ev_interpolant = EV_interpolants{ie};

                % 为当前收入冲击状态创建向量化的状态网格
                [K_states, Kpps_states] = ndgrid(cS.kGridV, cS.kppsGridV);

                % 计算收益
                k_income = K_states * M_age.r_mkt_t;
                kpps_income = Kpps_states * M_age.r_mkt_t;
                capital_tax = cS.tau_k * k_income;

                if a_idx <= cS.aR_new
                    % 工作期
                    labor_income_gross = M_age.w_t * cS.ageEffV_new(a_idx) * epsilon_state;
                    pension_benefit = 0;
                    payg_tax = cS.theta_t * labor_income_gross;
                    basic_resources = K_states * market_return_factor + labor_income_gross - (payg_tax + capital_tax);
                else
                    % 退休期
                    labor_income_gross = 0;
                    pension_benefit = b_age_val;
                    payg_tax = 0;
                    basic_resources = K_states * market_return_factor + pension_benefit - capital_tax;
                end

                % 计算冲击支出
                shock_expenditure = 0;
                if ie == cS.nw + 1
                    shock_expenditure = cS.kappa_young * basic_resources;
                elseif ie == cS.nw + 2
                    shock_expenditure = cS.kappa_old * basic_resources;
                end

                net_basic_resources = basic_resources - shock_expenditure;

                % 合并(ik, ikpps)双重循环为单个索引
                total_state_combinations = cS.nk * cS.nkpps;

                % 预分配输出数组
                val_results = zeros(total_state_combinations, 1);
                cPol_results = zeros(total_state_combinations, 1);
                kPol_results = zeros(total_state_combinations, 1);
                cPpsPol_results = zeros(total_state_combinations, 1);

                % 使用parfor并行化状态组合
                parfor state_idx = 1:total_state_combinations
                    % 将线性索引转换为(ik, ikpps)
                    [ik, ikpps] = ind2sub([cS.nk, cS.nkpps], state_idx);

                    available_basic = net_basic_resources(ik, ikpps);
                    current_kpps = cS.kppsGridV(ikpps);

                    % 初始化所有局部变量
                    kpps_prime_fixed_local = 0;
                    actual_kpps_contribution = [];

                    if a_idx <= cS.aR_new
                        % 工作期决策逻辑
                        [K_prime_choices, Kpps_prime_choices] = ndgrid(cS.kGridV, cS.kppsGridV);

                        kpps_contribution = Kpps_prime_choices - current_kpps * (1 + M_age.r_mkt_t);
                        kpps_contribution = max(0, kpps_contribution);

                        if labor_income_gross > 0
                            taxable_labor_income = max(0, labor_income_gross - payg_tax - kpps_contribution);
                            labor_tax_grid = cS.tau_l * taxable_labor_income;
                        else
                            labor_tax_grid = zeros(size(kpps_contribution));
                        end

                        actual_kpps_contribution = Kpps_prime_choices - current_kpps * (1 + M_age.r_mkt_t);
                        actual_kpps_contribution = max(0, actual_kpps_contribution);

                        c_expend_choices = available_basic - K_prime_choices - actual_kpps_contribution - labor_tax_grid;

                    else
                        % 退休期决策逻辑
                        current_kpps_total = current_kpps * (1 + M_age.r_mkt_t);
                        period_withdrawal_rate = cS.pps_withdrawal_rate;
                        mandatory_withdrawal = current_kpps_total * period_withdrawal_rate;
                        kpps_prime_fixed_local = current_kpps_total - mandatory_withdrawal;
                        pps_withdrawal_tax = cS.pps_tax_rate_withdrawal * mandatory_withdrawal;
                        net_pps_withdrawal = mandatory_withdrawal - pps_withdrawal_tax;

                        K_prime_choices = cS.kGridV';
                        c_expend_choices = available_basic + net_pps_withdrawal - K_prime_choices;
                        Kpps_prime_choices = repmat(kpps_prime_fixed_local, size(K_prime_choices));

                        % 为退休期设置虚拟的actual_kpps_contribution
                        actual_kpps_contribution = zeros(size(K_prime_choices));
                    end

                    % 可行性检查
                    feasible_mask = c_expend_choices > 0;

                    if a_idx <= cS.aR_new
                        max_contribution_mask = actual_kpps_contribution <= (labor_income_gross * cS.pps_max_contrib_frac);
                        contrib_limit_mask = actual_kpps_contribution <= cS.pps_contrib_limit;
                        feasible_mask = feasible_mask & max_contribution_mask & contrib_limit_mask;
                    end

                    if ~any(feasible_mask(:))
                        % 没有可行选择
                        val_results(state_idx) = -1e20;
                        kPol_results(state_idx) = cS.kMin;
                        cPpsPol_results(state_idx) = cS.kppsGridV(1);
                        cPol_results(state_idx) = 0;
                        continue;
                    end

                    % 只保留可行的选择
                    feasible_k_primes = K_prime_choices(feasible_mask);
                    feasible_kpps_primes = Kpps_prime_choices(feasible_mask);
                    feasible_c_expends = c_expend_choices(feasible_mask);

                    % 计算消费和效用
                    feasible_c_choices = feasible_c_expends / (1 + cS.tau_c);
                    [~, util_c_vec] = model_setup_utils.CES_utility(feasible_c_choices, cS.sigma, cS);

                    % 计算期望价值 (确保数据格式正确)
                    feasible_k_primes_col = feasible_k_primes(:);
                    feasible_kpps_primes_col = feasible_kpps_primes(:);
                    ev_vec = ev_interpolant(feasible_k_primes_col, feasible_kpps_primes_col);
                    ev_vec(isnan(ev_vec)) = -1e10;

                    % 计算遗赠效用
                    total_bequest = feasible_k_primes + feasible_kpps_primes;
                    util_bequest_vec = model_setup_utils.bequest_utility(total_bequest, cS);

                    % 确保所有向量都是列向量，避免维度不匹配
                    util_c_vec = util_c_vec(:);
                    ev_vec = ev_vec(:);
                    util_bequest_vec = util_bequest_vec(:);

                    % 计算总效用
                    total_values = util_c_vec + effective_discount_factor * ev_vec + bequest_discount_factor * util_bequest_vec;

                    % 找到最优选择 (确保返回标量)
                    [best_val, best_idx] = max(total_values(:));

                    if isfinite(best_val)
                        val_results(state_idx) = best_val;
                        kPol_results(state_idx) = feasible_k_primes(best_idx);
                        cPol_results(state_idx) = feasible_c_choices(best_idx);
                        cPpsPol_results(state_idx) = feasible_kpps_primes(best_idx);
                    else
                        val_results(state_idx) = -1e20;
                        kPol_results(state_idx) = cS.kMin;
                        if a_idx <= cS.aR_new
                            cPpsPol_results(state_idx) = cS.kppsGridV(1);
                        else
                            cPpsPol_results(state_idx) = kpps_prime_fixed_local;
                        end
                        cPol_results(state_idx) = 0;
                    end
                end

                % 将线性结果转换为矩阵形式
                val_age(:, :, ie) = reshape(val_results, [cS.nk, cS.nkpps]);
                cPol_age_q(:, :, ie) = reshape(cPol_results, [cS.nk, cS.nkpps]);
                kPol_age(:, :, ie) = reshape(kPol_results, [cS.nk, cS.nkpps]);
                cPpsPol_age_choice(:, :, ie) = reshape(cPpsPol_results, [cS.nk, cS.nkpps]);
            end
        end

        function [cPol_age_q, kPol_age, cPpsPol_age_choice, val_age] = HHSolutionByAge_VFI_ExpenditureShock_Vectorized_with_pps(a_idx, vPrime_kkppse_next, M_age, b_age_val, paramS_age, cS)
            % [PPS扩展版本 - 修正的税收递延制度 v2]
            % 核心修改：修正最后一期遗赠逻辑，使其与预算约束一致

            % 抑制MESHGRID插值性能警告
            warning('off', 'MATLAB:griddedInterpolant:MeshgridEval2DWarnId');

            val_age = -1e20 * ones(cS.nk, cS.nkpps, cS.nw_expanded);
            cPol_age_q = zeros(cS.nk, cS.nkpps, cS.nw_expanded);
            kPol_age = zeros(cS.nk, cS.nkpps, cS.nw_expanded);
            cPpsPol_age_choice = zeros(cS.nk, cS.nkpps, cS.nw_expanded);

            % --- 最后一期的处理逻辑 [已修正] ---
            if a_idx == cS.aD_new
                [K_grid, Kpps_grid] = ndgrid(cS.kGridV, cS.kppsGridV);

                % 1. 计算所有资产的税后总财富
                k_after_return = K_grid .* (1 + M_age.r_mkt_t);
                k_capital_tax = cS.tau_k .* (K_grid .* M_age.r_mkt_t);
                k_after_tax_value = k_after_return - k_capital_tax;

                kpps_after_return = Kpps_grid .* (1 + M_age.r_mkt_t);
                kpps_withdrawal_tax = cS.pps_tax_rate_withdrawal .* kpps_after_return;
                kpps_after_tax_value = kpps_after_return - kpps_withdrawal_tax;

                pension_after_tax = b_age_val;
                total_after_tax_wealth = k_after_tax_value + kpps_after_tax_value + pension_after_tax;

                if cS.phi_bequest > 0
                    % 2. 在消费和遗赠之间进行最优权衡
                    if abs(cS.sigma - 1) < 1e-6
                        optimal_c_share = 1 / (1 + cS.phi_bequest);
                    else
                        optimal_c_share = 1 / (1 + cS.phi_bequest^(1/cS.sigma));
                    end

                    c_expend_final = optimal_c_share * total_after_tax_wealth;
                    c_final = c_expend_final / (1 + cS.tau_c);
                    optimal_bequest = (1 - optimal_c_share) * total_after_tax_wealth;

                    % 3. 计算效用
                    [~, util_c] = model_setup_utils.CES_utility(c_final, cS.sigma, cS);
                    util_bequest = model_setup_utils.bequest_utility(optimal_bequest, cS);
                    util_final = util_c + util_bequest;

                    % 4. [核心修正] 将最优遗赠分配到策略函数中
                    %    一个合理的假设是，所有财富被整合到流动性最强的常规账户中进行遗赠
                    k_prime_final = optimal_bequest;
                    kpps_prime_final = zeros(size(Kpps_grid)); % PPS账户被清空

                    for ie = 1:cS.nw_expanded
                        cPol_age_q(:,:,ie) = c_final;
                        kPol_age(:,:,ie) = k_prime_final;
                        cPpsPol_age_choice(:,:,ie) = kpps_prime_final;
                        val_age(:,:,ie) = util_final;
                    end
                else
                    % 无遗赠动机：消费掉所有资源
                    c_expend_final = total_after_tax_wealth;
                    final_c = c_expend_final / (1 + cS.tau_c);
                    [~, util_c] = model_setup_utils.CES_utility(final_c, cS.sigma, cS);

                    for ie = 1:cS.nw_expanded
                        cPol_age_q(:,:,ie) = final_c;
                        val_age(:,:,ie) = util_c;
                        kPol_age(:,:,ie) = 0;
                        cPpsPol_age_choice(:,:,ie) = 0;
                    end
                end
                return;
            end

            % --- 向量化计算开始 ---
            effective_discount_factor = (cS.beta ^ cS.time_Step) * cS.s_pathV(a_idx);
            bequest_discount_factor = (cS.beta ^ cS.time_Step) * (1 - cS.s_pathV(a_idx));

            % 计算期望价值矩阵（现在是四维的）
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

            % 准备二维期望价值插值器
            EV_interpolants = cell(cS.nw_expanded, 1);
            for ie_current = 1:cS.nw_expanded
                EV_interpolants{ie_current} = griddedInterpolant({cS.kGridV, cS.kppsGridV}, EV_matrix(:,:,ie_current), 'linear', 'none');
            end

            market_return_factor = 1 + M_age.r_mkt_t;

            % --- 核心向量化部分：工作期与退休期分别处理 ---
            for ie = 1:cS.nw_expanded
                epsilon_state = paramS_age.leGridV(ie);
                ev_interpolant = EV_interpolants{ie};

                % 为当前收入冲击状态创建向量化的状态网格
                [K_states, Kpps_states] = ndgrid(cS.kGridV, cS.kppsGridV); % [nk x nkpps]

                % [修正] 分别计算普通资产和PPS资产的收益
                k_income = K_states * M_age.r_mkt_t;         % [nk x nkpps] 普通资产收益
                kpps_income = Kpps_states * M_age.r_mkt_t;   % [nk x nkpps] PPS资产收益

                % [修正] 只对普通资产收益征收资本利得税，PPS投资收益免税
                capital_tax = cS.tau_k * k_income; % [nk x nkpps]

                if a_idx <= cS.aR_new
                    % === 工作期：PPS资金锁定，不能取出 ===
                    labor_income_gross = M_age.w_t * cS.ageEffV_new(a_idx) * epsilon_state;
                    pension_benefit = 0;

                    payg_tax = cS.theta_t * labor_income_gross;

                    % [修正] 工作期可支配总资源不包括PPS资产（锁定状态）
                    basic_resources = K_states * market_return_factor + labor_income_gross - (payg_tax + capital_tax); % [nk x nkpps]

                else
                    % === 退休期：可以取出PPS资金 ===
                    labor_income_gross = 0;
                    pension_benefit = b_age_val;
                    payg_tax = 0;

                    % [修正] 退休期可支配总资源包括PPS资产，但取出时需征税
                    % 这里我们假设可以灵活取出PPS资金，具体取出多少在后面的优化中决定
                    basic_resources = K_states * market_return_factor + pension_benefit - capital_tax; % [nk x nkpps]
                end

                % 计算冲击支出（基于基础资源）
                shock_expenditure = 0;
                if ie == cS.nw + 1
                    shock_expenditure = cS.kappa_young * basic_resources;
                elseif ie == cS.nw + 2
                    shock_expenditure = cS.kappa_old * basic_resources;
                end

                net_basic_resources = basic_resources - shock_expenditure; % [nk x nkpps]

                % 对于每个 (k, kpps) 状态组合，计算最优选择
                for ik = 1:cS.nk
                    for ikpps = 1:cS.nkpps
                        available_basic = net_basic_resources(ik, ikpps);
                        current_kpps = cS.kppsGridV(ikpps);

                        if a_idx <= cS.aR_new
                            % === 工作期决策逻辑 ===
                            % 决策变量：k' (下一期普通资产) 和 pps_contribution (新增PPS缴费)
                            % 预算约束：c + k' + pps_contribution = available_basic - labor_tax_adjusted

                            % 创建决策网格：k' 和 下一期PPS总额 kpps'
                            [K_prime_choices, Kpps_prime_choices] = ndgrid(cS.kGridV, cS.kppsGridV); % [nk x nkpps]

                            % 计算新增PPS缴费
                            kpps_contribution = Kpps_prime_choices - current_kpps * (1 + M_age.r_mkt_t); % [nk x nkpps]

                            % 确保PPS缴费非负
                            kpps_contribution = max(0, kpps_contribution);

                            % [修正] 计算税前扣除后的劳动税
                            if labor_income_gross > 0
                                taxable_labor_income = max(0, labor_income_gross - payg_tax - kpps_contribution);
                                labor_tax_grid = cS.tau_l * taxable_labor_income; % [nk x nkpps]
                            else
                                labor_tax_grid = zeros(size(kpps_contribution));
                            end

                            % 重新计算实际PPS缴费（基于减税后的能力）
                            % 这里我们重新计算kpps_contribution，确保它能从税前收入中支付
                            actual_kpps_contribution = Kpps_prime_choices - current_kpps * (1 + M_age.r_mkt_t);
                            actual_kpps_contribution = max(0, actual_kpps_contribution);

                            % 计算消费支出
                            c_expend_choices = available_basic - K_prime_choices - actual_kpps_contribution - labor_tax_grid; % [nk x nkpps]

                        else
                            % === 退休期决策逻辑（外生固定取出率） ===
                            % 决策变量：仅k' (下一期普通资产)
                            % PPS取出金额由外生固定取出率决定

                            % [修正] 外生固定PPS取出率
                            current_kpps_total = current_kpps * (1 + M_age.r_mkt_t);

                            % 计算5年期强制取出率 (基于年度取出率)
                            period_withdrawal_rate = cS.pps_withdrawal_rate;
                            mandatory_withdrawal = current_kpps_total * period_withdrawal_rate;

                            % 下期PPS资产水平（外生确定）
                            kpps_prime_fixed = current_kpps_total - mandatory_withdrawal;

                            % [修正] 对PPS取出征税
                            pps_withdrawal_tax = cS.pps_tax_rate_withdrawal * mandatory_withdrawal;

                            % 计算净PPS取出金额
                            net_pps_withdrawal = mandatory_withdrawal - pps_withdrawal_tax;

                            % 创建决策网格：仅对k'进行优化
                            K_prime_choices = cS.kGridV'; % [1 x nk]

                            % 计算消费支出
                            c_expend_choices = available_basic + net_pps_withdrawal - K_prime_choices; % [1 x nk]

                            % 为所有K_prime_choices复制相同的kpps_prime_fixed
                            Kpps_prime_choices = repmat(kpps_prime_fixed, size(K_prime_choices)); % [1 x nk]
                        end

                        % 找到可行的选择（消费支出 > 0）
                        feasible_mask = c_expend_choices > 0; % [nk x nkpps] 工作期 或 [1 x nk] 退休期

                        % [修正] 在工作期添加额外的可行性检查
                        if a_idx <= cS.aR_new
                            % 工作期：确保PPS缴费不超过合理限制
                            max_contribution_mask = actual_kpps_contribution <= (labor_income_gross * cS.pps_max_contrib_frac); % 最多缴费指定比例的工资
                            contrib_limit_mask = actual_kpps_contribution <= cS.pps_contrib_limit; % 不超过缴费限额
                            feasible_mask = feasible_mask & max_contribution_mask & contrib_limit_mask;
                        else
                            % 退休期：PPS取出金额外生确定，只需检查消费可行性
                            % 不需要额外的可行性检查，因为取出金额已经外生固定
                        end

                        if ~any(feasible_mask(:))
                            % 没有可行选择，设置为最低值
                            val_age(ik, ikpps, ie) = -1e20;
                            kPol_age(ik, ikpps, ie) = cS.kMin;
                            cPpsPol_age_choice(ik, ikpps, ie) = cS.kppsGridV(1);  % 使用PPS网格的最小值
                            cPol_age_q(ik, ikpps, ie) = 0;
                            continue;
                        end

                        % 只保留可行的选择
                        feasible_k_primes = K_prime_choices(feasible_mask);
                        feasible_kpps_primes = Kpps_prime_choices(feasible_mask);
                        feasible_c_expends = c_expend_choices(feasible_mask);

                        % 向量化计算所有可行选择的消费和效用
                        feasible_c_choices = feasible_c_expends / (1 + cS.tau_c);

                        % 计算消费效用（向量化）
                        [~, util_c_vec] = model_setup_utils.CES_utility(feasible_c_choices, cS.sigma, cS);

                        % 计算期望价值（向量化插值）
                        ev_vec = ev_interpolant(feasible_k_primes, feasible_kpps_primes);

                        % 处理NaN值
                        ev_vec(isnan(ev_vec)) = -1e10;

                        % 计算遗赠效用（向量化）- 基于总遗赠额
                        total_bequest = feasible_k_primes + feasible_kpps_primes;
                        util_bequest_vec = model_setup_utils.bequest_utility(total_bequest, cS);

                        % 计算总效用
                        total_values = util_c_vec + effective_discount_factor * ev_vec + bequest_discount_factor * util_bequest_vec;

                        % 找到最优选择
                        [best_val, best_idx] = max(total_values);

                        if isfinite(best_val)
                            val_age(ik, ikpps, ie) = best_val;
                            kPol_age(ik, ikpps, ie) = feasible_k_primes(best_idx);
                            cPol_age_q(ik, ikpps, ie) = feasible_c_choices(best_idx);
                            cPpsPol_age_choice(ik, ikpps, ie) = feasible_kpps_primes(best_idx);
                        else
                            val_age(ik, ikpps, ie) = -1e20;
                            kPol_age(ik, ikpps, ie) = cS.kMin;
                            % [修正] 对于退休期，即使没有可行解，也要记录外生的kpps_prime
                            if a_idx <= cS.aR_new
                                cPpsPol_age_choice(ik, ikpps, ie) = cS.kppsGridV(1);  % 工作期使用PPS网格的最小值
                            else
                                cPpsPol_age_choice(ik, ikpps, ie) = kpps_prime_fixed;  % 退休期使用外生固定值
                            end
                            cPol_age_q(ik, ikpps, ie) = 0;
                        end
                    end
                end
            end
        end

        function [cPolM, kPolM, cPpsPolM, valM] = HHSolution_VFI_with_pps(M_vfi, paramS_vfi, cS_vfi)
            % [PPS扩展版本] 求解包含PPS决策的家庭问题
            valM = -Inf(cS_vfi.nk, cS_vfi.nkpps, cS_vfi.nw_expanded, cS_vfi.aD_new);
            cPolM = zeros(cS_vfi.nk, cS_vfi.nkpps, cS_vfi.nw_expanded, cS_vfi.aD_new);
            kPolM = zeros(cS_vfi.nk, cS_vfi.nkpps, cS_vfi.nw_expanded, cS_vfi.aD_new);
            cPpsPolM = zeros(cS_vfi.nk, cS_vfi.nkpps, cS_vfi.nw_expanded, cS_vfi.aD_new);

            bV_payg_vfi = zeros(1, cS_vfi.aD_new);
            if cS_vfi.aR_new < cS_vfi.aD_new
                bV_payg_vfi((cS_vfi.aR_new+1):cS_vfi.aD_new) = M_vfi.b_t;
            end

            for a_idx = cS_vfi.aD_new : -1 : 1
                vPrime_kkppse_next = [];
                if a_idx < cS_vfi.aD_new
                    vPrime_kkppse_next = valM(:,:,:,a_idx+1);
                end
                [cPolM(:,:,:,a_idx), kPolM(:,:,:,a_idx), cPpsPolM(:,:,:,a_idx), valM(:,:,:,a_idx)] = ...
                    main_steady_state_utils.HHSolutionByAge_VFI_ExpenditureShock_Vectorized_with_pps_ParFor(a_idx, vPrime_kkppse_next, M_vfi, bV_payg_vfi(a_idx), paramS_vfi, cS_vfi);
            end
        end

        function [k_prime_idx, kpps_prime_idx] = get_policy_index_matrix_with_pps(kPolM, cPpsPolM, cS)
            % [PPS扩展版本] 将连续的二维储蓄决策离散化为网格索引
            k_prime_idx = zeros(cS.nk, cS.nkpps, cS.nw_expanded, cS.aD_new, 'uint16');
            kpps_prime_idx = zeros(cS.nk, cS.nkpps, cS.nw_expanded, cS.aD_new, 'uint16');

            for ia = 1:cS.aD_new
                for ie = 1:cS.nw_expanded
                    for ik = 1:cS.nk
                        for ikpps = 1:cS.nkpps
                            % 获取连续空间下的最优储蓄决策
                            k_prime_continuous = kPolM(ik, ikpps, ie, ia);
                            kpps_prime_continuous = cPpsPolM(ik, ikpps, ie, ia);

                            % 离散化普通储蓄决策
                            affordable_k_indices = find(cS.kGridV <= k_prime_continuous);
                            if isempty(affordable_k_indices)
                                idx_k = 1;
                            else
                                idx_k = affordable_k_indices(end);
                            end

                            % 离散化PPS储蓄决策
                            affordable_kpps_indices = find(cS.kppsGridV <= kpps_prime_continuous);
                            if isempty(affordable_kpps_indices)
                                idx_kpps = 1;
                            else
                                idx_kpps = affordable_kpps_indices(end);
                            end

                            k_prime_idx(ik, ikpps, ie, ia) = idx_k;
                            kpps_prime_idx(ik, ikpps, ie, ia) = idx_kpps;
                        end
                    end
                end
            end
        end

        function Dist = solve_steady_state_distribution_with_pps(k_prime_idx, kpps_prime_idx, paramS, cS, Z_ss_norm)
            % [PPS扩展版本] 求解包含PPS状态的稳态分布
            Dist = zeros(cS.nk, cS.nkpps, cS.nw_expanded, cS.aD_new);

            % 1. 初始化新生儿分布
            % 新生儿只分布在 (k=0, kpps=0) 的状态
            dist_newborn_kkppse = zeros(cS.nk, cS.nkpps, cS.nw_expanded);
            % [修正] 处理维度匹配问题：确保正确设置新生儿分布
            newborn_probs = reshape(paramS.leProb1V, [1, 1, length(paramS.leProb1V)]);
            dist_newborn_kkppse(1, 1, 1:cS.nw) = newborn_probs;

            % 新生儿的总人口由 Z_ss_norm(1) 决定
            Dist(:, :, :, 1) = dist_newborn_kkppse * Z_ss_norm(1);

            % 2. 迭代计算幸存者的分布
            for ia = 1:(cS.aD_new - 1)
                dist_ia_kkppse = Dist(:, :, :, ia); % 当前年龄的四维分布
                dist_ia_plus_1_kkppse_unscaled = zeros(cS.nk, cS.nkpps, cS.nw_expanded);
                transition_matrix_next_age = paramS.TrProbM_by_age{ia + 1};

                % 计算未缩放的下一期分布
                for ik = 1:cS.nk
                    for ikpps = 1:cS.nkpps
                        for ie = 1:cS.nw_expanded
                            mass_at_state = dist_ia_kkppse(ik, ikpps, ie);
                            if mass_at_state < 1e-20, continue; end

                            ik_prime = k_prime_idx(ik, ikpps, ie, ia);
                            ikpps_prime = kpps_prime_idx(ik, ikpps, ie, ia);
                            transition_probs_e = transition_matrix_next_age(ie, :);

                            % 将质量转移到新的状态
                            % [修正] 处理维度匹配问题：确保右侧维度为 [1, 1, nw_expanded]
                            transition_probs_reshaped = reshape(transition_probs_e, [1, 1, length(transition_probs_e)]);
                            dist_ia_plus_1_kkppse_unscaled(ik_prime, ikpps_prime, :) = ...
                                dist_ia_plus_1_kkppse_unscaled(ik_prime, ikpps_prime, :) + mass_at_state * transition_probs_reshaped;
                        end
                    end
                end

                % 使用 Z_ss_norm 来重新缩放下一期的分布
                mass_at_ia = sum(dist_ia_kkppse, 'all');

                if mass_at_ia > 1e-12
                    rescale_factor = Z_ss_norm(ia+1) / mass_at_ia;
                    Dist(:, :, :, ia+1) = dist_ia_plus_1_kkppse_unscaled * rescale_factor;
                else
                    Dist(:, :, :, ia+1) = zeros(cS.nk, cS.nkpps, cS.nw_expanded);
                end
            end

            % 3. 最终检查
            final_sum = sum(Dist, 'all');
            if abs(final_sum - 1.0) > 1e-6
                warning('稳态联合分布总和(%.8f)不为1，可能存在会计不一致。', final_sum);
            end
        end

        function [c_val, tax_val, shock_expenditure, payg_tax, labor_tax, capital_tax, consumption_tax, pension_benefit] = backout_accounting_expenditure_shock_with_pps(k_now, kpps_now, k_prime, kpps_prime, ia, ie, epsilon_val, M_sim, cS)
            % [PPS扩展版本 - 完全重构的会计镜像]
            % [已修正预算检验逻辑 - v3]
            % 严格镜像HHSolution函数中的预算约束逻辑，确保完美的会计一致性
            % 统一预算恒等式：[期初k + 期初k_pps] + 收入 = 消费 + [期末k' + 期末k_pps'] + 税收 + 冲击支出

            if ~isfield(cS, 'theta_t'), cS.theta_t = cS.theta_path(1); end

            % === 第1步：计算总流入 (Total Inflow) ===
            market_return_factor = 1 + M_sim.r_mkt_t;
            inflow_from_k = k_now * market_return_factor;
            inflow_from_kpps = kpps_now * market_return_factor;

            labor_income_gross = 0; pension_benefit = 0;
            if ia <= cS.aR_new
                labor_income_gross = M_sim.w_t * cS.ageEffV_new(ia) * epsilon_val;
            else
                pension_benefit = M_sim.b_t;
            end
            inflow_from_income = labor_income_gross + pension_benefit;

            total_inflow = inflow_from_k + inflow_from_kpps + inflow_from_income;

            % === 第2步：计算税收（精确镜像HHSolution的逻辑） ===
            payg_tax = cS.theta_t * labor_income_gross;

            if ia <= cS.aR_new
                % === 工作期税收计算 ===
                capital_tax = cS.tau_k * (k_now * M_sim.r_mkt_t);
                pps_contribution = max(0, kpps_prime - inflow_from_kpps);
                taxable_labor_income = max(0, labor_income_gross - payg_tax - pps_contribution);
                labor_tax = cS.tau_l * taxable_labor_income;
            else
                % === 退休期税收计算（外生固定取出率） [已修正] ===
                capital_tax = cS.tau_k * (k_now * M_sim.r_mkt_t);

                % [核心修正] 完美镜像VFI中的外生提取规则
                current_kpps_total_accounting = inflow_from_kpps;
                period_withdrawal_rate = cS.pps_withdrawal_rate;

                % 强制提取额完全由规则确定，与kpps_prime无关
                kpps_withdrawal = current_kpps_total_accounting * period_withdrawal_rate;

                % 对这部分强制提取额征税
                labor_tax = cS.pps_tax_rate_withdrawal * kpps_withdrawal;
            end

            total_tax_paid_before_consumption_tax = payg_tax + capital_tax + labor_tax;

            % === 第3步：计算冲击支出（精确镜像决策函数逻辑） ===
            if ia <= cS.aR_new
                basic_resources_for_shock = inflow_from_k + labor_income_gross - (payg_tax + capital_tax);
            else
                basic_resources_for_shock = inflow_from_k + pension_benefit - capital_tax;
            end

            shock_expenditure = 0;
            if ie == cS.nw + 1
                shock_expenditure = cS.kappa_young * basic_resources_for_shock;
            elseif ie == cS.nw + 2
                shock_expenditure = cS.kappa_old * basic_resources_for_shock;
            end

            % === 第4步：反算含税消费支出 ===
            outflow_non_consumption = k_prime + kpps_prime + total_tax_paid_before_consumption_tax + shock_expenditure;
            c_expend_available = total_inflow - outflow_non_consumption;

            % === 第5步：计算消费和消费税 ===
            c_val = c_expend_available / (1 + cS.tau_c);
            consumption_tax = c_val * cS.tau_c;
            tax_val = total_tax_paid_before_consumption_tax + consumption_tax;

            % === 第6步：精确的预算平衡检验 ===
            total_consumption_expenditure = c_val * (1 + cS.tau_c);
            actual_total_outflow = k_prime + kpps_prime + total_tax_paid_before_consumption_tax + shock_expenditure + total_consumption_expenditure;
            budget_gap = total_inflow - actual_total_outflow;

            if abs(budget_gap) > 1e-8
                error('微观预算约束被违反！状态: (ia=%d, ie=%d, ik=%d, ikpps=%d), 预算缺口: %.3e\n总流入: %.6f, 总流出: %.6f', ...
                    ia, ie, find(abs(cS.kGridV - k_now) < 1e-10, 1), find(abs(cS.kppsGridV - kpps_now) < 1e-10, 1), budget_gap, total_inflow, actual_total_outflow);
            end
        end
        
        function [K_agg, C_utility_agg, Tax_agg, Bequest_tax_agg, L_agg, PensionIn_agg, PensionOut_agg, ShockExp_agg] = aggregate_expenditure_shock_with_pps(Dist, k_prime_idx, kpps_prime_idx, M_sim, cS, paramS)
            % [PPS扩展版本] 聚合包含PPS的宏观变量
            if ~isfield(cS, 'theta_t'), cS.theta_t = cS.theta_path(1); end
            K_agg = 0; C_utility_agg = 0; Tax_agg = 0; Bequest_tax_agg = 0; L_agg = 0; PensionIn_agg = 0; PensionOut_agg = 0; ShockExp_agg = 0;

            % 根据上下文确定使用哪个存活率向量
            if isfield(cS, 'prob_survive_implied_ss0') && ~isempty(cS.prob_survive_implied_ss0)
                prob_survive_implied = cS.prob_survive_implied_ss0;
            elseif isfield(cS, 'prob_survive_implied_trans') && ~isempty(cS.prob_survive_implied_trans)
                prob_survive_implied = cS.prob_survive_implied_trans;
            else
                warning('未找到隐含存活率，使用外生存活率 s_pathV');
                prob_survive_implied = cS.s_pathV;
            end

            for ia = 1:cS.aD_new
                for ie = 1:cS.nw_expanded
                    for ik = 1:cS.nk
                        for ikpps = 1:cS.nkpps
                            mass = Dist(ik, ikpps, ie, ia);
                            if mass < 1e-20, continue; end

                            k_now = cS.kGridV(ik);
                            kpps_now = cS.kppsGridV(ikpps);
                            epsilon_val = paramS.leGridV(ie);

                            idx_k_prime = k_prime_idx(ik, ikpps, ie, ia);
                            idx_kpps_prime = kpps_prime_idx(ik, ikpps, ie, ia);
                            k_prime = cS.kGridV(idx_k_prime);
                            kpps_prime = cS.kppsGridV(idx_kpps_prime);

                            [c_val, tax_val, shock_exp, payg_tax, ~, ~, ~, pension_benefit] = ...
                                main_steady_state_utils.backout_accounting_expenditure_shock_with_pps(k_now, kpps_now, k_prime, kpps_prime, ia, ie, epsilon_val, M_sim, cS);

                            C_utility_agg = C_utility_agg + c_val * mass;
                            Tax_agg = Tax_agg + tax_val * mass;
                            ShockExp_agg = ShockExp_agg + shock_exp * mass;

                            % 核心修改：总资本包括两种资产
                            prob_survive = prob_survive_implied(ia);
                            K_agg = K_agg + (k_prime + kpps_prime) * mass * prob_survive;
                            prob_death = 1 - prob_survive;
                            Bequest_tax_agg = Bequest_tax_agg + (k_prime + kpps_prime) * mass * prob_death;

                            if ia <= cS.aR_new, L_agg = L_agg + (cS.ageEffV_new(ia) * epsilon_val) * mass; end
                            PensionIn_agg = PensionIn_agg + payg_tax * mass;
                            PensionOut_agg = PensionOut_agg + pension_benefit * mass;
                        end
                    end
                end
            end
        end

                function [K_p_model_out, ss, Dist, k_prime_idx, kpps_prime_idx, V, kPolM, cPpsPolM] = calculate_aggregates_for_Kg_with_pps(K_p_guess, K_g_guess, Z_ss_norm, cS, paramS)
            % [PPS扩展版本] 计算包含PPS的稳态聚合量
            % [修改] 动态网格设置移至 main_run_transition.m 中，此函数专注于聚合计算
            if K_p_guess <= 0, K_p_guess = 1e-8; end
            if K_g_guess <= 0, K_g_guess = 1e-8; end
            A_ss = cS.A; theta_ss = cS.theta_path(1);

            L_iter_tol = 1e-7; L_iter_max = 50; L_damping = 0.5;

            % 计算年龄-冲击联合分布
            mean_e_by_age = zeros(cS.aD_new, 1);
            e_dist_by_age = zeros(cS.aD_new, cS.nw_expanded);
            e_dist_by_age(1, 1:cS.nw) = paramS.leProb1V';
            mean_e_by_age(1) = e_dist_by_age(1, :) * paramS.leGridV(:);
            for ia = 1:(cS.aD_new - 1)
                e_dist_by_age(ia+1, :) = e_dist_by_age(ia, :) * paramS.TrProbM_by_age{ia+1};
                mean_e_by_age(ia+1) = e_dist_by_age(ia+1, :) * paramS.leGridV(:);
            end

            L_guess = 0;
            for ia = 1:cS.aR_new
                L_guess = L_guess + Z_ss_norm(ia) * cS.ageEffV_new(ia) * mean_e_by_age(ia);
            end

            % 劳动供给迭代
            for iter = 1:L_iter_max
                M_prices = main_steady_state_utils.get_prices_at_t(K_p_guess, K_g_guess, L_guess, A_ss, cS);
                M_for_hh = M_prices;

                total_wage_bill = M_prices.w_t * L_guess;
                mass_retirees_ss = sum(Z_ss_norm((cS.aR_new+1):end));
                if mass_retirees_ss > 1e-9
                    M_for_hh.b_t = (theta_ss * total_wage_bill) / mass_retirees_ss;
                else
                    M_for_hh.b_t = 0;
                end

                % 确保PPS激活
                cS_ss = cS;
                cS_ss.pps_active = true;
                cS_ss.theta_t = theta_ss;

                % [调试] 确保网格设置正确
                if cS_ss.nkpps == 1
                    error('PPS网格未正确设置！cS.nkpps = %d，应该 > 1', cS_ss.nkpps);
                end

                % 计算临时劳动供给
                [~, temp_kPolM, temp_cPpsPolM, ~] = main_steady_state_utils.HHSolution_VFI_with_pps(M_for_hh, paramS, cS_ss);
                [temp_k_prime_idx, temp_kpps_prime_idx] = main_steady_state_utils.get_policy_index_matrix_with_pps(temp_kPolM, temp_cPpsPolM, cS_ss);
                temp_Dist = main_steady_state_utils.solve_steady_state_distribution_with_pps(temp_k_prime_idx, temp_kpps_prime_idx, paramS, cS_ss, Z_ss_norm);

                L_model = 0;
                for ia = 1:cS.aR_new
                    for ie = 1:cS.nw_expanded
                        for ik = 1:cS.nk
                            for ikpps = 1:cS.nkpps
                                mass = temp_Dist(ik, ikpps, ie, ia);
                                if mass > 0
                                    epsilon_val = paramS.leGridV(ie);
                                    L_model = L_model + (cS.ageEffV_new(ia) * epsilon_val) * mass;
                                end
                            end
                        end
                    end
                end

                L_error = abs(L_model - L_guess);
                if L_error < L_iter_tol, break; end
                L_guess = L_damping * L_guess + (1 - L_damping) * L_model;
            end

            % 计算最终的 V, kPol, cPpsPol, Dist 和宏观量
            M_prices = main_steady_state_utils.get_prices_at_t(K_p_guess, K_g_guess, L_guess, A_ss, cS);
            M_for_hh = M_prices;
            total_wage_bill = M_prices.w_t * L_guess;
            mass_retirees_ss = sum(Z_ss_norm((cS.aR_new+1):end));
            if mass_retirees_ss > 1e-9
                M_for_hh.b_t = (theta_ss * total_wage_bill) / mass_retirees_ss;
            else
                M_for_hh.b_t = 0;
            end

            % 调用PPS版本的VFI得到最终的价值和政策函数
            [~, kPolM, cPpsPolM, V] = main_steady_state_utils.HHSolution_VFI_with_pps(M_for_hh, paramS, cS_ss);
            [k_prime_idx, kpps_prime_idx] = main_steady_state_utils.get_policy_index_matrix_with_pps(kPolM, cPpsPolM, cS_ss);
            Dist = main_steady_state_utils.solve_steady_state_distribution_with_pps(k_prime_idx, kpps_prime_idx, paramS, cS_ss, Z_ss_norm);

            % 使用PPS版本的聚合函数
            [K_p_model_out, C_utility_final, Tax_final, Bequest_tax_final, ~, ~, ~, ~] = ...
                main_steady_state_utils.aggregate_expenditure_shock_with_pps(Dist, k_prime_idx, kpps_prime_idx, M_for_hh, cS_ss, paramS);

            % 填充完整的 ss 结构体
            ss = struct();
            ss.K_private = K_p_guess;
            ss.K_public = K_g_guess;
            ss.K_total = K_p_guess + K_g_guess;
            ss.L = L_guess;
            ss.Y_from_production = M_prices.Y_t;
            ss.w = M_prices.w_t;
            ss.r_mkt = M_prices.r_mkt_t;
            ss.b = M_for_hh.b_t;
            ss.Bequest_tax = Bequest_tax_final;
            ss.Regular_tax = Tax_final;
        end

        function F_error = system_of_equations_Kg_with_pps(x, Z_ss_norm, cS, paramS)
            % [PPS扩展版本] 为fsolve提供包含PPS的方程组
            K_p_guess = x(1);
            K_g_guess = x(2);

            % 调用PPS版本的核心计算函数
            [K_p_model, ss] = main_steady_state_utils.calculate_aggregates_for_Kg_with_pps(K_p_guess, K_g_guess, Z_ss_norm, cS, paramS);

            % 方程1的误差: 私人资本供给 - 私人资本需求
            error_Kp = K_p_guess - K_p_model;

            % 方程2的误差: 公共投资 - 公共资本折旧
            Gov_Revenue_total = ss.Regular_tax + ss.Bequest_tax;
            I_g_model = cS.lambda_g * Gov_Revenue_total;
            Depreciation_g_model = cS.ddk_g * K_g_guess;
            error_Kg = I_g_model - Depreciation_g_model;

            F_error = [error_Kp; error_Kg];
        end

        function [ss, eq_found, Dist, k_prime_idx, kpps_prime_idx, V, kPolM, cPpsPolM] = solve_steady_state_iter_Kg_with_pps(Z_ss_norm, cS, paramS, verbose, x0_guess, solver_method)
            % [PPS扩展版本] 迭代求解包含PPS的稳态
            % [新增] 支持智能初始猜测值，提高大TFP变化时的求解效率
            % [新增] 支持多种求解器：'fsolve' (默认), 'surrogateopt', 'hybrid'
            if nargin < 4, verbose = true; end
            if nargin < 6, solver_method = 'fsolve'; end  % 默认使用fsolve

            system_wrapper = @(x) main_steady_state_utils.system_of_equations_Kg_with_pps(x, Z_ss_norm, cS, paramS);

            % [修改] 使用智能初始猜测值，如果没有提供则使用默认值
            if nargin < 5 || isempty(x0_guess)
                % 默认初始猜测值（适用于小TFP的情况）
                k_p_guess_initial = 3.5;
                k_g_guess_initial = 1.0;
                x0 = [k_p_guess_initial, k_g_guess_initial];
                if verbose
                    fprintf('   使用默认初始猜测值: Kp=%.2f, Kg=%.2f\n', x0(1), x0(2));
                end
            else
                % 使用提供的智能初始猜测值
                x0 = x0_guess;
                if verbose
                    fprintf('   使用智能初始猜测值: Kp=%.2f, Kg=%.2f\n', x0(1), x0(2));
                end
            end

            % [新增] 根据求解器方法选择不同的求解策略
            switch lower(solver_method)
                case 'fsolve'
                    [x_eq, eq_found] = main_steady_state_utils.solve_with_fsolve(system_wrapper, x0, verbose);

                case 'surrogateopt'
                    [x_eq, eq_found] = main_steady_state_utils.solve_with_surrogateopt(system_wrapper, x0, verbose);

                case 'hybrid'
                    % 简单混合策略：先尝试fsolve，失败则用surrogateopt
                    [x_eq, eq_found] = main_steady_state_utils.solve_with_fsolve(system_wrapper, x0, verbose);
                    if ~eq_found
                        if verbose, fprintf('   fsolve失败，切换到surrogateopt...\n'); end
                        [x_eq, eq_found] = main_steady_state_utils.solve_with_surrogateopt(system_wrapper, x0, verbose);
                    end

                case 'robust'
                    % 鲁棒混合策略：多阶段求解，包括随机初始化
                    [x_eq, eq_found] = main_steady_state_utils.solve_with_hybrid_robust(system_wrapper, x0, verbose);

                otherwise
                    error('未知的求解器方法: %s。支持的方法: fsolve, surrogateopt, hybrid, robust', solver_method);
            end

            if ~eq_found
                if verbose, warning('所有求解器都未能找到均衡解'); end
                ss = []; Dist = []; k_prime_idx = []; kpps_prime_idx = []; V = []; kPolM = []; cPpsPolM = [];
                return;
            end

            K_p_eq = x_eq(1);
            K_g_eq = x_eq(2);

            % 求解成功后，进行最后一次计算以获得所有变量
            [~, ss, Dist, k_prime_idx, kpps_prime_idx, V, kPolM, cPpsPolM] = main_steady_state_utils.calculate_aggregates_for_Kg_with_pps(K_p_eq, K_g_eq, Z_ss_norm, cS, paramS);

            % 根据 verbose 决定输出详细程度
            if verbose
                % 需要适配四维分布的显示函数
                main_steady_state_utils.display_national_accounts_gov_investment_with_pps(ss, Dist, k_prime_idx, kpps_prime_idx, cS, paramS);
            else
                main_steady_state_utils.check_national_accounts_with_pps(ss, Dist, k_prime_idx, kpps_prime_idx, cS, paramS);
            end
        end

        function [ss, Dist, V, kPolM, cPpsPolM] = solve_steady_state_complete_with_pps(cS_ss, paramS, params_ext, verbose, x0_guess, solver_method)
            % [PPS扩展版本 - 修正的税收递延制度] 完整的包含PPS的稳态求解器
            % [新增] 支持智能初始猜测值，提高大TFP变化时的求解效率
            % [新增] 支持多种求解器：'fsolve' (默认), 'surrogateopt', 'hybrid'
            if nargin < 4, verbose = true; end
            if nargin < 5, x0_guess = []; end  % 默认为空，让下层函数使用默认值
            if nargin < 6, solver_method = 'fsolve'; end  % 默认使用fsolve

            % 1. 使用模型内生的人口分布
            age_mass = ones(cS_ss.aD_new, 1);
            for ia = 1:(cS_ss.aD_new - 1)
                age_mass(ia+1) = age_mass(ia) * cS_ss.s_pathV(ia);
            end

            % 如果外部提供了特定年份的人口分布，则使用外部的
            if isfield(params_ext, 'Z') && ~isempty(params_ext.Z)
                Z_ss_norm = params_ext.Z;
            else
                Z_ss_norm = age_mass / sum(age_mass);
            end

            % 2. 覆盖cS中的外生变量
            cS_ss.A = params_ext.A;
            if isfield(params_ext, 'g_A_ss'), cS_ss.g_A_ss = params_ext.g_A_ss; end
            cS_ss.theta_path = params_ext.theta;

            % [修正] 确保PPS激活并验证相关参数
            cS_ss.pps_active = true;

            % [新增] 验证PPS相关参数设置
            if ~isfield(cS_ss, 'nkpps') || cS_ss.nkpps <= 1
                error('PPS稳态求解器要求 cS.nkpps > 1，当前 nkpps = %d', cS_ss.nkpps);
            end
            if ~isfield(cS_ss, 'kppsGridV') || length(cS_ss.kppsGridV) ~= cS_ss.nkpps
                error('PPS网格 kppsGridV 长度(%d) 与 nkpps(%d) 不匹配', length(cS_ss.kppsGridV), cS_ss.nkpps);
            end

            % if verbose
            %     fprintf('PPS制度参数检查:\n');
            %     fprintf('  - PPS网格点数: %d\n', cS_ss.nkpps);
            %     fprintf('  - PPS资产范围: [%.3f, %.3f]\n', min(cS_ss.kppsGridV), max(cS_ss.kppsGridV));
            %     fprintf('  - 资本利得税率: %.2f%% (仅适用于普通资产)\n', cS_ss.tau_k * 100);
            %     fprintf('  - 劳动税率: %.2f%% (适用于工资和PPS取出)\n', cS_ss.tau_l * 100);
            % end

            % 计算稳态的隐含存活率向量
            prob_survive_implied_ss0 = zeros(cS_ss.aD_new, 1);
            for ia = 1:(cS_ss.aD_new - 1)
                if Z_ss_norm(ia) > 1e-12
                    prob_survive_implied_ss0(ia) = Z_ss_norm(ia+1) / Z_ss_norm(ia);
                else
                    prob_survive_implied_ss0(ia) = 0;
                end
            end
            prob_survive_implied_ss0(cS_ss.aD_new) = 0;
            cS_ss.prob_survive_implied_ss0 = prob_survive_implied_ss0;

            % 3. 调用PPS版本的稳态求解器
            [ss, eq_found, Dist, k_prime_idx, kpps_prime_idx, V, kPolM, cPpsPolM] = ...
                main_steady_state_utils.solve_steady_state_iter_Kg_with_pps(Z_ss_norm, cS_ss, paramS, verbose, x0_guess, solver_method);

            if ~eq_found
                warning('稳态求解失败！');
                ss = []; Dist = []; V = []; kPolM = []; cPpsPolM = [];
            end
        end

        % 辅助函数：适配四维分布的国民账户检查
        function is_ok = check_national_accounts_with_pps(ss, Dist, k_prime_idx, kpps_prime_idx, cS, paramS)
            % [PPS扩展版本] 用于校准循环的简化版国民账户检查
            M_sim = struct('r_mkt_t', ss.r_mkt, 'w_t', ss.w, 'b_t', ss.b);

            [~, C_utility_agg, Tax_agg, Bequest_tax_agg, ~, ~, ~, ShockExp_agg] = ...
                main_steady_state_utils.aggregate_expenditure_shock_with_pps(Dist, k_prime_idx, kpps_prime_idx, M_sim, cS, paramS);

            Y_prod = ss.Y_from_production;
            Gov_Revenue_total = Tax_agg + Bequest_tax_agg;
            I_g_agg = cS.lambda_g * Gov_Revenue_total;
            G_c_agg = (1 - cS.lambda_g) * Gov_Revenue_total;
            I_p_agg_gross = cS.ddk * ss.K_private;
            I_total_agg = I_p_agg_gross + I_g_agg;
            C_total_agg = C_utility_agg + ShockExp_agg;
            Y_exp_actual = C_total_agg + I_total_agg + G_c_agg;

            error_gdp = abs(Y_exp_actual - Y_prod);
            if error_gdp < 1e-4
                fprintf('   国民账户核算 (PPS版本): 通过 (误差: %.3e)\n', error_gdp);
                is_ok = true;
            else
                fprintf('   国民账户核算 (PPS版本): 失败! (误差: %.3e)\n', error_gdp);
                is_ok = false;
            end
        end

        function display_national_accounts_gov_investment_with_pps(ss, Dist, k_prime_idx, kpps_prime_idx, cS, paramS)
            % [PPS扩展版本] 详细的国民账户报告
            fprintf('\n\n========================================================================\n');
            fprintf('===     国民经济核算详细报告 (政府投资模型版 - 包含PPS)     ===\n');
            fprintf('========================================================================\n');

            M_sim = struct('r_mkt_t', ss.r_mkt, 'w_t', ss.w, 'b_t', ss.b);

            [~, C_utility_agg, Tax_agg, Bequest_tax_agg, L_agg_check, PensionIn_agg, PensionOut_agg, ShockExp_agg] = ...
                main_steady_state_utils.aggregate_expenditure_shock_with_pps(Dist, k_prime_idx, kpps_prime_idx, M_sim, cS, paramS);

            % 宏观变量定义
            Y_prod = ss.Y_from_production;
            K_p = ss.K_private;
            K_g = ss.K_public;

            Gov_Revenue_total = Tax_agg + Bequest_tax_agg;
            I_g_agg = cS.lambda_g * Gov_Revenue_total;
            G_c_agg = (1 - cS.lambda_g) * Gov_Revenue_total;
            I_p_agg_gross = cS.ddk * K_p;
            I_total_agg = I_p_agg_gross + I_g_agg;
            C_total_agg = C_utility_agg + ShockExp_agg;
            Y_exp_actual = C_total_agg + I_total_agg + G_c_agg;

            fprintf('--- A. 宏观产出与支出核算 (包含PPS) ---\n');
            fprintf('   生产法 GDP (Y_prod):         %.6f\n', Y_prod);
            fprintf('   支出法 GDP (C+I_total+G_c):  %.6f\n', Y_exp_actual);
            fprintf('   ------------------------------------\n');
            fprintf('   核算误差 (Y_exp - Y_prod):     %.3e (此值应接近0)\n', Y_exp_actual - Y_prod);
            fprintf('   总消费 (C):                  %.6f (占GDP: %.2f%%)\n', C_total_agg, C_total_agg/Y_prod*100);
            fprintf('   总投资 (I_total = I_p+I_g):  %.6f (占GDP: %.2f%%)\n', I_total_agg, I_total_agg/Y_prod*100);
            fprintf('   政府消费 (G_c):              %.6f (占GDP: %.2f%%)\n', G_c_agg, G_c_agg/Y_prod*100);

            fprintf('\n--- B. 政府与养老金体系核算 ---\n');
            fprintf('   政府总收入 (T):              %.6f\n', Gov_Revenue_total);
            fprintf('   政府总支出 (G_c+I_g):        %.6f\n', G_c_agg + I_g_agg);
            fprintf('   政府预算平衡 (T - G_c - I_g):%.3e (此值应接近0)\n', Gov_Revenue_total - (G_c_agg + I_g_agg));
            fprintf('   养老金体系平衡 (收入-支出):  %.3e (此值应接近0)\n', PensionIn_agg - PensionOut_agg);
            fprintf('   劳动供给核算误差 (L_agg-L_ss): %.3e\n', L_agg_check - ss.L);

            fprintf('\n--- C. 关键宏观比率 ---\n');
            fprintf('   私人资本产出比 (K_p/Y):      %.4f\n', K_p / Y_prod);
            fprintf('   公共资本产出比 (K_g/Y):      %.4f\n', K_g / Y_prod);
            fprintf('   总资本产出比 (K_total/Y):    %.4f\n', (K_p + K_g) / Y_prod);
            fprintf('   [注意] 总资本现在包括PPS资产\n');
            fprintf('\n========================================================================\n');
        end

        % 测试函数：验证修正后的PPS稳态求解器
        function test_pps_steady_state_solver()
            % [修正版本] 验证修正后的税收递延PPS制度是否正确工作
            fprintf('\n=== 测试修正后的PPS稳态求解器 ===\n');

            try
                % 创建简化的测试参数
                cS = model_setup_utils.ParameterValues();
                cS.time_Step = 5;
                cS.T_sim = 10;
                cS.nk = 15;     % 增加网格点以提高精度
                cS.nkpps = 15;  % 增加PPS网格点
                cS.npps = 3;
                cS.pps_active = true;

                % 设置合理的税率参数
                cS.tau_k = 0.20;  % 资本利得税率
                cS.tau_l = 0.25;  % 劳动税率（也用于PPS取出征税）
                cS.tau_c = 0.05;  % 消费税率

                cS = model_setup_utils.generateGrids(cS);

                % 创建简化的参数
                paramS = struct();
                [paramS.leGridV, paramS.TrProbM_by_age, paramS.leProb1V, cS.nw_expanded] = ...
                    model_setup_utils.EarningProcess_AgeDependent(cS);

                % 创建测试的外部参数
                params_ext = struct();
                params_ext.A = 1.0;
                params_ext.theta = 0.15;
                params_ext.Z = ones(cS.aD_new, 1) / cS.aD_new; % 均匀分布

                % 测试修正后的PPS版本
                fprintf('   正在测试修正后的PPS稳态求解器...\n');
                fprintf('   PPS制度特征: 工作期锁定+缴费税前扣除, 退休期取出征税, 投资收益免税\n');
                [ss_pps, ~, ~, ~, ~] = main_steady_state_utils.solve_steady_state_complete_with_pps(cS, paramS, params_ext, false);

                if ~isempty(ss_pps)
                    fprintf('   ✅ 修正后PPS稳态求解器测试通过\n');
                    fprintf('   结果: K_p=%.4f, K_g=%.4f, Y=%.4f, K/Y=%.4f\n', ...
                        ss_pps.K_private, ss_pps.K_public, ss_pps.Y_from_production, ...
                        (ss_pps.K_private + ss_pps.K_public)/ss_pps.Y_from_production);
                else
                    fprintf('   ❌ 修正后PPS稳态求解器测试失败\n');
                end

                % 比较与原版本的差异
                fprintf('   正在测试原版本的稳态求解器作为对比...\n');
                cS_old = cS;
                cS_old.pps_active = false;
                cS_old.nkpps = 1;
                cS_old = model_setup_utils.generateGrids(cS_old);

                [ss_old, ~, ~, ~] = main_steady_state_utils.solve_steady_state_complete(cS_old, paramS, params_ext, false);

                if ~isempty(ss_old)
                    fprintf('   ✅ 原版本稳态求解器测试通过\n');
                    fprintf('   结果: K_p=%.4f, K_g=%.4f, Y=%.4f, K/Y=%.4f\n', ...
                        ss_old.K_private, ss_old.K_public, ss_old.Y_from_production, ...
                        (ss_old.K_private + ss_old.K_public)/ss_old.Y_from_production);

                    % [修正] 更详细的差异分析
                    fprintf('   经济影响分析 (PPS制度 vs 无PPS):\n');
                    k_total_change = (ss_pps.K_private + ss_pps.K_public) - (ss_old.K_private + ss_old.K_public);
                    k_total_pct = k_total_change / (ss_old.K_private + ss_old.K_public) * 100;
                    tax_change = (ss_pps.Regular_tax + ss_pps.Bequest_tax) - (ss_old.Regular_tax + ss_old.Bequest_tax);
                    tax_pct = tax_change / (ss_old.Regular_tax + ss_old.Bequest_tax) * 100;

                    if k_total_change > 0
                        effect_desc = '储蓄激励效果';
                    else
                        effect_desc = '储蓄抑制效果';
                    end
                    fprintf('     总资本存量变化: %.4f (%.2f%%) - %s\n', ...
                        k_total_change, k_total_pct, effect_desc);

                    if tax_change > 0
                        tax_desc = '税收增加';
                    else
                        tax_desc = '税收减少';
                    end
                    fprintf('     政府总税收变化: %.4f (%.2f%%) - %s\n', ...
                        tax_change, tax_pct, tax_desc);
                    fprintf('     产出变化: %.4f (%.2f%%)\n', ...
                        ss_pps.Y_from_production - ss_old.Y_from_production, ...
                        (ss_pps.Y_from_production - ss_old.Y_from_production) / ss_old.Y_from_production * 100);

                    % [新增] 经济学直觉检验
                    fprintf('   经济学直觉检验:\n');
                    if k_total_change > 0
                        fprintf('     ✅ PPS税收优惠提高了总储蓄，符合理论预期\n');
                    else
                        fprintf('     ⚠️ PPS制度导致总储蓄下降，可能需要检查参数设置\n');
                    end

                    if abs(tax_pct) < 5
                        fprintf('     ✅ 政府税收变化适中，PPS制度具有财政可持续性\n');
                    else
                        fprintf('     ⚠️ 政府税收变化较大(%.1f%%)，需要关注财政影响\n', tax_pct);
                    end

                else
                    fprintf('   ❌ 原版本稳态求解器测试失败\n');
                end

            catch ME
                fprintf('   ❌ 测试过程中发生错误: %s\n', ME.message);
                if ~isempty(ME.stack)
                    fprintf('   错误位置: %s (第%d行)\n', ME.stack(1).name, ME.stack(1).line);
                end
            end

            fprintf('=== 修正后PPS稳态求解器测试完成 ===\n\n');
        end

        function test_adaptive_grid_system()
            % [新增] 测试智能自适应网格系统 (外部预设版本)
            fprintf('\n=== 智能自适应网格系统测试 (外部预设版本) ===\n');

            try
                % 创建测试参数
                cS = model_setup_utils.ParameterValues();
                cS.time_Step = 5;
                cS.nk = 10;
                cS.nkpps = 10;
                cS.pps_active = true;

                % 初始网格（默认情况）
                cS_default = model_setup_utils.generateGrids(cS);
                default_k_max = max(cS_default.kGridV);
                default_kpps_max = max(cS_default.kppsGridV);

                fprintf('   默认网格范围:\n');
                fprintf('     常规资产: [%.2f, %.2f]\n', min(cS_default.kGridV), default_k_max);
                fprintf('     PPS资产: [%.2f, %.2f]\n', min(cS_default.kppsGridV), default_kpps_max);

                % 测试不同的K_p猜测值
                test_cases = [
                    struct('K_p', 5, 'scenario', '常规情况'),
                    struct('K_p', 50, 'scenario', '中等增长'),
                    struct('K_p', 200, 'scenario', '极端增长')
                    ];

                fprintf('\n   自适应网格测试:\n');
                for i = 1:length(test_cases)
                    K_p_test = test_cases(i).K_p;
                    scenario = test_cases(i).scenario;

                    % 计算自适应网格上限
                    GRID_SCALING_FACTOR = 10;
                    k_max_adaptive = GRID_SCALING_FACTOR * K_p_test;
                    kpps_max_adaptive = 0.5 * k_max_adaptive;

                    % 生成自适应网格
                    cS_adaptive = model_setup_utils.generateGrids(cS, 'k_max', k_max_adaptive, 'kpps_max', kpps_max_adaptive);

                    fprintf('   %s (K_p=%.0f):\n', scenario, K_p_test);
                    fprintf('     自适应常规资产: [%.2f, %.2f] (扩展%.1f倍)\n', ...
                        min(cS_adaptive.kGridV), max(cS_adaptive.kGridV), max(cS_adaptive.kGridV)/default_k_max);
                    fprintf('     自适应PPS资产: [%.2f, %.2f] (扩展%.1f倍)\n', ...
                        min(cS_adaptive.kppsGridV), max(cS_adaptive.kppsGridV), max(cS_adaptive.kppsGridV)/default_kpps_max);

                    % 验证网格覆盖度
                    coverage_ratio = K_p_test / max(cS_adaptive.kGridV);
                    if coverage_ratio < 0.8
                        fprintf('     ✅ 网格覆盖充足 (猜测值占网格上限%.1f%%)\n', coverage_ratio * 100);
                    else
                        fprintf('     ⚠️ 网格覆盖不足 (猜测值占网格上限%.1f%%)\n', coverage_ratio * 100);
                    end
                end

                fprintf('\n   ✅ 智能自适应网格系统测试完成\n');
                fprintf('   核心优势: 基于智能猜测值预设网格范围，一次性优化配置\n');
                fprintf('   适用范围: 外部调用generateGrids预设，求解器内部专注计算\n');
                fprintf('   性能特点: 避免求解过程中反复网格调整，提高计算效率\n');

            catch ME
                fprintf('   ❌ 测试过程中发生错误: %s\n', ME.message);
            end

            fprintf('=== 智能自适应网格系统测试完成 (外部预设+双轨制) ===\n\n');
        end

        function test_solver_comparison()
            % [新增] 测试不同求解器的性能比较
            fprintf('\n=== 求解器性能比较测试 ===\n');

            try
                % 创建测试参数
                cS = model_setup_utils.ParameterValues();
                cS.time_Step = 5;
                cS.nk = 10;     % 使用较小的网格以加快测试
                cS.nkpps = 10;
                cS.pps_active = true;
                cS = model_setup_utils.generateGrids(cS);

                paramS = struct();
                [paramS.leGridV, paramS.TrProbM_by_age, paramS.leProb1V, cS.nw_expanded] = ...
                    model_setup_utils.EarningProcess_AgeDependent(cS);

                params_ext = struct('A', 2.0, 'theta', 0.15, 'Z', ones(cS.aD_new, 1) / cS.aD_new);

                % 测试不同求解器
                solvers = {'fsolve', 'surrogateopt', 'hybrid'};
                results = struct();

                for i = 1:length(solvers)
                    solver_name = solvers{i};
                    fprintf('\n   测试求解器: %s\n', solver_name);

                    tic;
                    [ss_result, ~, ~, ~, ~] = main_steady_state_utils.solve_steady_state_complete_with_pps(...
                        cS, paramS, params_ext, false, [], solver_name);
                    elapsed_time = toc;

                    if ~isempty(ss_result)
                        results.(solver_name).success = true;
                        results.(solver_name).time = elapsed_time;
                        results.(solver_name).K_total = ss_result.K_private + ss_result.K_public;
                        results.(solver_name).Y = ss_result.Y_from_production;
                        results.(solver_name).KY_ratio = results.(solver_name).K_total / results.(solver_name).Y;

                        fprintf('     ✅ 成功 - 耗时: %.2f秒, K/Y=%.4f\n', elapsed_time, results.(solver_name).KY_ratio);
                    else
                        results.(solver_name).success = false;
                        results.(solver_name).time = elapsed_time;
                        fprintf('     ❌ 失败 - 耗时: %.2f秒\n', elapsed_time);
                    end
                end

                % 总结比较
                fprintf('\n   === 求解器比较总结 ===\n');
                fprintf('   求解器\t\t成功\t耗时(秒)\tK/Y比率\n');
                fprintf('   %s\n', repmat('-', 1, 50));

                for i = 1:length(solvers)
                    solver_name = solvers{i};
                    if results.(solver_name).success
                        fprintf('   %-12s\t✅\t%.2f\t\t%.4f\n', solver_name, ...
                            results.(solver_name).time, results.(solver_name).KY_ratio);
                    else
                        fprintf('   %-12s\t❌\t%.2f\t\t--\n', solver_name, results.(solver_name).time);
                    end
                end

                fprintf('\n   建议:\n');
                fprintf('   - 对于常规问题，优先使用 fsolve (速度最快)\n');
                fprintf('   - 对于复杂参数，使用 surrogateopt (最鲁棒)\n');
                fprintf('   - 不确定时，使用 hybrid (平衡速度和鲁棒性)\n');

            catch ME
                fprintf('   ❌ 测试过程中发生错误: %s\n', ME.message);
            end

            fprintf('=== 求解器性能比较测试完成 ===\n\n');
        end

        % =======================================================
        % == 辅助函数：多种求解器实现
        % =======================================================

        function [x_eq, eq_found] = solve_with_fsolve(system_wrapper, x0, verbose)
            % [fsolve求解器] 使用传统的fsolve方法求解方程组
            if verbose
                fsolve_display = 'iter';
            else
                fsolve_display = 'none';
            end
            options = optimoptions('fsolve', 'Display', fsolve_display, 'TolFun', 1e-9, 'TolX', 1e-9, 'MaxIterations', 500);

            if verbose, fprintf('\n--- 启动 fsolve 求解器 (求解 [K_p, K_g] - PPS版本) ---\n'); end
            [x_eq, ~, exitflag] = fsolve(system_wrapper, x0, options);
            if verbose, fprintf('--- fsolve 求解完成 ---\n'); end

            eq_found = (exitflag > 0);
        end

        function [x_eq, eq_found] = solve_with_surrogateopt(system_wrapper, x0, verbose)
            % [surrogateopt求解器] 使用全局优化方法求解方程组
            % 核心思想：将F(x)=0的方程组问题转换为min ||F(x)||^2的优化问题

            % 定义目标函数：最小化方程组残差的平方和
            objective_function = @(x) norm(system_wrapper(x), 2)^2;

            % [智能边界] 基于初始猜测值动态调整搜索边界
            % 这与我们的自适应网格系统协同工作
            lb = [0.1, 0.01];    % K_p, K_g的下界

            % 根据初始猜测值设置上界，支持极端TFP增长情况
            if max(x0) > 20
                % 对于极端情况，扩展搜索空间
                scale_factor = max(2, max(x0) / 20);
                ub = [200 * scale_factor, 200 * scale_factor];
            else
                ub = [100, 100];     % 常规情况的上界
            end

            % 确保初始猜测值在界内
            x0 = max(lb, min(ub, x0));

            % 设置surrogateopt选项
            if verbose
                display_option = 'iter';
            else
                display_option = 'none';
            end

            options = optimoptions('surrogateopt', ...
                'Display', display_option, ...
                'MaxFunctionEvaluations', 1000, ...
                'ObjectiveLimit', 1e-18, ...  % 当目标函数值小于此值时停止
                'UseParallel', false, ...
                'InitialPoints', x0);  % 指定初始点

            if verbose, fprintf('\n--- 启动 surrogateopt 求解器 (求解 [K_p, K_g] - PPS版本) ---\n'); end

            % [修正] surrogateopt 不需要初始点作为第四个参数，使用 InitialPoints 选项代替
            [x_eq, fval, exitflag] = surrogateopt(objective_function, lb, ub, options);

            if verbose, fprintf('--- surrogateopt 求解完成 ---\n'); end

            % 判断是否收敛：目标函数值足够小
            tolerance = 1e-6;
            eq_found = (exitflag > 0) && (fval < tolerance);

            if verbose
                if eq_found
                    fprintf('   surrogateopt收敛成功: 目标函数值 = %.3e\n', fval);
                else
                    fprintf('   surrogateopt收敛失败: 目标函数值 = %.3e (exitflag = %d)\n', fval, exitflag);
                end
            end
        end

        function [x_eq, eq_found] = solve_with_hybrid_robust(system_wrapper, x0, verbose)
            % [鲁棒混合求解器] 增强的混合求解策略
            % 1. 首先尝试fsolve（快速）
            % 2. 如果失败，尝试surrogateopt（鲁棒）
            % 3. 如果还失败，尝试多个随机初始点的surrogateopt

            if verbose, fprintf('\n--- 启动混合鲁棒求解器 ---\n'); end

            % 第1阶段：尝试fsolve
            [x_eq, eq_found] = main_steady_state_utils.solve_with_fsolve(system_wrapper, x0, verbose);
            if eq_found
                if verbose, fprintf('   混合求解器: fsolve成功\n'); end
                return;
            end

            % 第2阶段：尝试surrogateopt（使用原初始点）
            if verbose, fprintf('   混合求解器: fsolve失败，尝试surrogateopt...\n'); end
            [x_eq, eq_found] = main_steady_state_utils.solve_with_surrogateopt(system_wrapper, x0, verbose);
            if eq_found
                if verbose, fprintf('   混合求解器: surrogateopt成功\n'); end
                return;
            end

            % 第3阶段：多点随机初始化的surrogateopt
            if verbose, fprintf('   混合求解器: 单点surrogateopt失败，尝试多点随机初始化...\n'); end

            % 第3阶段：使用多个随机初始点的surrogateopt
            num_random_starts = 5;
            lb = [0.1, 0.01];

            % [智能边界] 与其他求解器保持一致的边界设置
            if max(x0) > 20
                scale_factor = max(2, max(x0) / 20);
                ub = [200 * scale_factor, 200 * scale_factor];
            else
                ub = [100, 100];
            end

            % 生成多个随机初始点矩阵
            x0_matrix = zeros(num_random_starts, 2);
            for i = 1:num_random_starts
                x0_matrix(i, :) = lb + (ub - lb) .* rand(1, 2);
                if verbose, fprintf('   随机初始点 %d: [%.2f, %.2f]\n', i, x0_matrix(i, 1), x0_matrix(i, 2)); end
            end

            % 使用多点初始化的surrogateopt
            objective_function = @(x) norm(system_wrapper(x), 2)^2;
            options = optimoptions('surrogateopt', ...
                'Display', 'none', ...
                'MaxFunctionEvaluations', 1500, ...
                'ObjectiveLimit', 1e-18, ...
                'UseParallel', false, ...
                'InitialPoints', x0_matrix);

            if verbose, fprintf('   启动多点随机初始化的surrogateopt (共%d个初始点)...\n', num_random_starts); end
            [x_eq, fval, exitflag] = surrogateopt(objective_function, lb, ub, options);

            tolerance = 1e-6;
            eq_found = (exitflag > 0) && (fval < tolerance);

            if eq_found && verbose
                fprintf('   混合求解器: 多点surrogateopt成功，目标函数值 = %.3e\n', fval);
            end

            % 如果所有方法都失败了
            if verbose, fprintf('   混合求解器: 所有方法都失败了\n'); end
            x_eq = x0;  % 返回初始猜测值
            eq_found = false;
        end

    end
end