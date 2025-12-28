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
            % =========================================================================
            % == 函数: solve_steady_state_unified
            % == 版本: [v10 - 【求解器可切换最终版】]
            % ==
            % == 核心功能:
            % ==   - 增加 `solver_method` 参数，允许用户在 'fsolve' 和 'lsqnonlin' 之间切换。
            % ==   - 'fsolve': 默认无约束求解。
            % ==   - 'fsolve_transform': 使用变量变换法施加非负约束。
            % ==   - 'lsqnonlin': 使用边界约束进行稳健求解（推荐）。
            % ==   - 整合所有逻辑，使其成为一个统一、灵活的求解器入口。
            % =========================================================================
            if nargin < 4, verbose = true; end
            if nargin < 5, x0_guess = []; end
            if nargin < 6, solver_method = 'lsqnonlin'; end % 默认使用推荐的 lsqnonlin

            % --- 1. 参数设置 (不变) ---
            cS.A = params_ext.A;
            if isfield(params_ext, 'g_A_ss'), cS.g_A_ss = params_ext.g_A_ss; end
            if isfield(params_ext, 'n_ss'), cS.n_ss = params_ext.n_ss; end
            cS.theta_path = params_ext.theta; cS.tau_beq = 0.0;
            if isfield(params_ext, 'Z') && ~isempty(params_ext.Z), Z_ss_norm = params_ext.Z;
            else, error('错误：求解稳态需要一个明确的人口分布Z。'); end
            if ~isfield(cS, 'pps_active'), cS.pps_active = false; end
            if cS.pps_active && cS.nkpps <= 1, warning('PPS模式要求 nkpps > 1.');
            elseif ~cS.pps_active, cS.nkpps = 1; cS.kppsGridV = 0; end
            if ~isfield(cS, 'n_ss'), cS.n_ss = 0.0; end

            % --- 2. 准备求解 ---
            if isempty(x0_guess)
                x0 = [0.3336, 0.078, 0.3, 0.015];
            else
                x0 = x0_guess(1:4);
            end

            % --- 3. [核心] 根据选择的求解器进行调度 ---
            fprintf('\n>>> 正在使用求解器: %s <<<\n', upper(solver_method));

            switch lower(solver_method)
                case 'lsqnonlin'
                    system_wrapper = @(x) main_steady_state_utils_bgp.system_of_equations_steady_state(x, Z_ss_norm, cS, paramS, 'none');
                    lb = [1e-8, 1e-8, 1e-8, 1e-8];
                    ub = [Inf, Inf, Inf, Inf];
                    [x_eq, eq_found] = main_steady_state_utils_bgp.solve_with_lsqnonlin(system_wrapper, x0, lb, ub, verbose);
                    Bequest_final = x_eq(4);

                case 'fsolve_transform'
                    system_wrapper = @(x) main_steady_state_utils_bgp.system_of_equations_steady_state(x, Z_ss_norm, cS, paramS, 'transform');
                    x0_transformed = x0;
                    x0_transformed(4) = sqrt(max(0, x0(4)));
                    fprintf('   [约束求解] 使用fsolve+变量变换。原始遗赠猜测=%.4f, 变换后猜测=%.4f\n', x0(4), x0_transformed(4));
                    [x_eq_transformed, eq_found] = main_steady_state_utils_bgp.solve_with_fsolve(system_wrapper, x0_transformed, verbose);
                    if eq_found
                        x_eq = x_eq_transformed;
                        x_eq(4) = x_eq_transformed(4)^2; % 结果反变换
                        Bequest_final = x_eq(4);
                    else
                        x_eq = []; % 如果未找到解，返回空
                    end

                case 'fsolve'
                    system_wrapper = @(x) main_steady_state_utils_bgp.system_of_equations_steady_state(x, Z_ss_norm, cS, paramS, 'none');
                    fprintf('   [无约束求解] 使用标准fsolve。警告：可能尝试负值。\n');
                    [x_eq, eq_found] = main_steady_state_utils_bgp.solve_with_fsolve(system_wrapper, x0, verbose);
                    if eq_found, Bequest_final = x_eq(4); end

                otherwise
                    error('未知的求解器方法: %s。请选择 ''lsqnonlin'', ''fsolve_transform'', 或 ''fsolve''。', solver_method);
            end

            if ~eq_found || isempty(x_eq)
                warning('统一求解器未能找到均衡解！');
                ss = []; Dist = []; polS = []; valS = [];
                return;
            end

            % --- 4. 获取最终结果 ---
            [ss, Dist, polS, valS] = main_steady_state_utils_bgp.calculate_aggregates_unified(x_eq(1), x_eq(2), x_eq(3), Bequest_final, Z_ss_norm, cS, paramS);

            % --- 5. 验证与展示 ---
            if verbose, main_steady_state_utils_bgp.display_national_accounts_unified(ss, cS, Dist, polS, paramS); end
        end

        function F_error = system_of_equations_steady_state(x, Z_ss_norm, cS, paramS, method)
            % =========================================================================
            % == 函数: system_of_equations_steady_state
            % == 版本: [v33 - 带实时误差监控版]
            % =========================================================================
            if nargin < 5, method = 'none'; end

            K_private_total_guess = x(1);
            K_g_guess = x(2);
            L_guess = x(3);

            % [核心修改] 根据方法标志进行反变换
            if strcmp(method, 'transform')
                sqrt_Beq = x(4);
                Bequest_Total_guess = sqrt_Beq^2;
            else % 默认(none)或lsqnonlin的情况
                Bequest_Total_guess = x(4);
            end

            % 提前检查猜测值，防止无效输入导致VFI崩溃
            if K_private_total_guess <= 0 || K_g_guess <= 0 || L_guess <= 0
                F_error = [1e8; 1e8; 1e8; 1e8];
                return;
            end

            [ss, ~, ~, ~] = main_steady_state_utils_bgp.calculate_aggregates_unified(K_private_total_guess, K_g_guess, L_guess, Bequest_Total_guess, Z_ss_norm, cS, paramS);

            if isempty(ss)
                F_error = [1e8; 1e8; 1e8; 1e8];
                return;
            end

            % 市场出清方程部分保持不变
            g_A_period = (1+cS.g_A_ss)^cS.time_Step-1; n_period = (1+cS.n_ss)^cS.time_Step-1;
            g_total_period = (1+g_A_period)*(1+n_period)-1; % 使用BGP总增长率

            I_p_gross_demand = (g_total_period + cS.ddk) * K_private_total_guess;
            error_Kp = ss.Saving_private_flow_Gross - I_p_gross_demand;

            I_g_gross_demand = (g_total_period + cS.ddk_g) * K_g_guess;
            error_Kg = ss.Saving_public_flow_Gross - I_g_gross_demand;

            error_L = ss.L_hat - L_guess;

            error_Beq = ss.Bequest_generated_agg - Bequest_Total_guess;

            % --- [ 新增的调试代码行 ] ---
            % 实时打印四个市场的误差，用于监控求解过程
            % fprintf('   [Solver Debug] Errors -> Kp: %+.4e | Kg: %+.4e | L: %+.4e | Beq: %+.4e\n', error_Kp, error_Kg, error_L, error_Beq);
            % --- [ 新增代码结束 ] ---

            F_error = [error_Kp; error_Kg; error_L; error_Beq];
        end

        % =========================================================================
% == 完整函数代码: main_steady_state_utils_bgp.m -> calculate_aggregates_unified
% == 版本: [v44 - 价格函数接口校正版]
% ==
% == 核心修改:
% ==   - [一致性] 本函数现在完全依赖于以“去趋势化(hat)”单位作为输入和输出
% ==     的 get_prices_at_t (v4) 函数。
% ==   - [字段名修正] 所有从价格函数 M_prices 中获取的变量（如工资、产出），
% ==     都明确使用 w_hat_t 和 Y_hat_t，以消除歧义。
% ==   - [逻辑正确性] 所有后续的计算，包括养老金、家庭收入、政府收入等，
% ==     都基于这些正确的“去趋势化”价格进行，确保了整个核算框架的内在一致性。
% =========================================================================
function [ss, Dist, polS, valS] = calculate_aggregates_unified(K_private_total_guess, K_g_guess, L_guess, Bequest_Total_guess, Z_ss_norm, cS, paramS)
    % --- 1. 初始化、价格计算 ---
    if K_private_total_guess <= 0, K_private_total_guess = 1e-8; end; if K_g_guess <= 0, K_g_guess = 1e-8; end; if L_guess <= 0, L_guess = 1e-8; end;
    
    % [调用修正] 在稳态求解中，所有输入给价格函数的都是hat单位, A_t的角色由A_hat=1.0扮演
    M_prices = main_steady_state_utils_bgp.get_prices_at_t(K_private_total_guess, K_g_guess, L_guess, 1.0, cS);
    M_for_hh = M_prices;
    % [字段名修正] 明确传递给VFI的是去趋势化的工资
    M_for_hh.w_t = M_prices.w_hat_t; 
    M_for_hh.r_mkt_t = M_prices.r_mkt_t;

    % --- 2. 养老金参数设定 (引入上限) ---
    mass_retirees_ss = sum(Z_ss_norm((cS.aR_new+1):end));

    if isfield(cS, 'endogenous_theta_mode') && cS.endogenous_theta_mode
        % DB 模式
        % [字段名修正] 养老金目标基于去趋势化的工资w_hat
        b_t_target = cS.payg_replacement_rate * M_prices.w_hat_t;
        total_pension_outlay_target = b_t_target * mass_retirees_ss;
        total_wage_base = M_prices.w_hat_t * L_guess;
        if total_wage_base > 1e-12
            theta_theoretical = total_pension_outlay_target / total_wage_base;
        else, theta_theoretical = 0; end
        if ~isfield(cS, 'theta_max'), cS.theta_max = 1.0; end
        theta_ss = min(theta_theoretical, cS.theta_max);
        total_pension_revenue_actual = theta_ss * total_wage_base;
        b_t_actual = total_pension_revenue_actual / max(1e-9, mass_retirees_ss);
        M_for_hh.b_t = b_t_actual;
        cS.theta_t = theta_ss;
    else
        % DC 模式
        theta_ss = cS.theta_path(1);
        % [字段名修正] 养老金总额基于去趋势化的工资w_hat
        total_pension_pot = theta_ss * M_prices.w_hat_t * L_guess;
        M_for_hh.b_t = total_pension_pot / max(1e-9, mass_retirees_ss);
        cS.theta_t = theta_ss;
    end

    % --- 3. VFI和分布计算 ---
    total_population_mass = sum(Z_ss_norm(:));
    M_for_hh.beq_transfer_pers = Bequest_Total_guess / max(1e-12, total_population_mass);
    [polS, valS] = main_steady_state_utils_bgp.HHSolution_VFI_unified(M_for_hh, paramS, cS);
    if any(isinf(valS(:))) || any(isnan(valS(:))), ss = []; Dist = []; polS = []; valS = []; return; end
    Dist = main_steady_state_utils_bgp.solve_steady_state_distribution_unified(polS, paramS, cS, Z_ss_norm);

    % --- 4. 聚合微观决策变量 ---
    K_private_hat_agg_vfi = 0; K_pps_hat_agg_vfi = 0;
    L_agg = 0; C_agg = 0; Bequest_generated_agg = 0;
    Shock_exp_agg = 0; Pension_out_agg = 0; Pension_in_agg = 0; Regular_tax_agg = 0;
    PPS_contrib_agg = 0; PPS_withdrawal_agg = 0; PPS_tax_agg = 0;
    for ia = 1:cS.aD_new
        mass_dist = Dist(:,:,:,ia);
        k_prime_dist = polS(ia).k_prime;
        K_private_hat_agg_vfi = K_private_hat_agg_vfi + sum(k_prime_dist .* mass_dist, 'all');
        kpps_prime_dist = 0;
        if cS.pps_active, kpps_prime_dist = polS(ia).kpps_prime;
            K_pps_hat_agg_vfi = K_pps_hat_agg_vfi + sum(kpps_prime_dist .* mass_dist, 'all');
            PPS_contrib_agg = PPS_contrib_agg + sum(polS(ia).pps_contrib .* mass_dist, 'all');
            PPS_withdrawal_agg = PPS_withdrawal_agg + sum(polS(ia).pps_withdrawal .* mass_dist, 'all');
            PPS_tax_agg = PPS_tax_agg + sum(polS(ia).pps_tax .* mass_dist, 'all'); end
        Bequest_generated_agg = Bequest_generated_agg + sum((k_prime_dist + kpps_prime_dist) .* mass_dist, 'all') * (1 - cS.s_pathV(ia));
        if ia <= cS.aR_new, mass_by_epsilon = squeeze(sum(mass_dist, [1,2])); L_agg = L_agg + sum( (cS.ageEffV_new(ia) * paramS.leGridV') .* mass_by_epsilon' , 'all'); end
        C_agg = C_agg + sum(polS(ia).c .* mass_dist, 'all');
        Shock_exp_agg = Shock_exp_agg + sum(polS(ia).shock_exp .* mass_dist, 'all');
        Regular_tax_agg = Regular_tax_agg + sum(polS(ia).tax_regular .* mass_dist, 'all');
        Pension_in_agg = Pension_in_agg + sum(polS(ia).tax_payg .* mass_dist, 'all');
        Pension_out_agg = Pension_out_agg + sum(polS(ia).pension_out .* mass_dist, 'all');
    end

    % --- [关键修正] 强制PAYG预算平衡 ---
    Pension_in_agg = Pension_out_agg;

    % --- 5. 宏观经济流量核算 ---
    Consumption_Total = C_agg + Shock_exp_agg;
    % [字段名修正] 家庭要素收入基于去趋势化的工资 w_hat
    Household_Factor_Income_Gross = (M_prices.w_hat_t * L_agg) + ((M_prices.r_mkt_t + cS.ddk) * K_private_total_guess);
    Household_Disposable_Income_Gross = Household_Factor_Income_Gross + Pension_out_agg - Pension_in_agg - Regular_tax_agg - PPS_tax_agg;
    Saving_private_flow_Gross = Household_Disposable_Income_Gross - Consumption_Total;
    % [字段名修正] 公共资本回报基于去趋势化的产出 Y_hat
    Public_Capital_Return = M_prices.Y_hat_t - Household_Factor_Income_Gross;
    Gov_Total_Revenue = Regular_tax_agg + PPS_tax_agg + Public_Capital_Return + Pension_in_agg;
    Gov_Discretionary_Resources = Gov_Total_Revenue - Pension_out_agg;
    G_c = (1 - cS.lambda_g) * Gov_Discretionary_Resources;
    Gov_Total_Outlay_NonSaving = G_c + Pension_out_agg;
    Saving_public_flow_Gross = Gov_Total_Revenue - Gov_Total_Outlay_NonSaving;

    % --- 6. 存量核算 (不变) ---
    K_private_implied = 0; k_grid_col_vec = cS.kGridV(:);
    for ia = 1:cS.aD_new, mass_at_k = squeeze(sum(Dist(:,:,:,ia), [2, 3])); K_private_implied = K_private_implied + sum(k_grid_col_vec .* mass_at_k); end

    % --- 7. 填充ss结构体 ---
    ss = struct();
    ss.K_private_begin_hat = K_private_total_guess; ss.K_private_implied = K_private_implied;
    ss.K_private_hat = K_private_hat_agg_vfi; ss.K_pps_hat = K_pps_hat_agg_vfi;
    ss.L_hat = L_agg; ss.K_public_hat = K_g_guess;
    % [字段名修正] 存储去趋势化的产出和工资
    ss.Y_from_production_hat = M_prices.Y_hat_t; ss.w_hat = M_prices.w_hat_t; ss.r_mkt = M_prices.r_mkt_t;
    ss.C_agg = Consumption_Total; ss.G_c = G_c;
    ss.Depreciation_p = cS.ddk * K_private_total_guess; ss.Depreciation_g = cS.ddk_g * K_g_guess;
    ss.Bequest_generated_agg = Bequest_generated_agg; ss.Bequest_distributed_agg = Bequest_Total_guess;
    ss.Bequest_tax_revenue = 0; ss.Regular_tax = Regular_tax_agg;
    ss.Pension_in = Pension_in_agg; ss.Pension_out = Pension_out_agg;
    ss.Saving_private_flow_Gross = Saving_private_flow_Gross; ss.Saving_public_flow_Gross = Saving_public_flow_Gross;
    ss.Public_Capital_Return = Public_Capital_Return;
    ss.PPS_contrib_agg = PPS_contrib_agg; ss.PPS_withdrawal_agg = PPS_withdrawal_agg; ss.PPS_tax_agg = PPS_tax_agg;
    ss.Shock_exp_agg = Shock_exp_agg;
    ss.theta = theta_ss;
end


        % =========================================================================
% == 完整函数代码: main_steady_state_utils_bgp.m -> HHSolution_VFI_unified
% == 版本: [v2.0 - VFI调用路径强制一致性修复版]
% ==
% == 核心修正:
% ==   - [!!!根本性修复!!!] 解决了稳态VFI和转轨VFI使用不同代码分支的
% ==     核心不一致问题。
% ==   - [正确实现] 在此函数中，我们显式计算稳态的BGP增长因子 `g_A_period_ss`，
% ==     并将其以 `g_A_t_plus_1` 的字段名添加到传递给底层VFI的结构体中。
% ==   - [一致性保证] 此修改强制稳态求解器与转轨求解器在调用核心VFI时，
% ==     走过完全相同的 `if isfield(M_age, 'g_A_t_plus_1')` 分支。这从
% ==     根源上消除了两者之间因代码路径不同而产生的任何数值行为差异，
% ==     确保了转轨路径在期末能精确收敛到稳态解。
% =========================================================================
        function [polS, valS] = HHSolution_VFI_unified(M_vfi, paramS_vfi, cS_vfi)
            valS = -Inf(cS_vfi.nk, cS_vfi.nkpps, cS_vfi.nw_expanded, cS_vfi.aD_new);
            polS_cell = cell(cS_vfi.aD_new, 1);
            bV_payg_vfi = zeros(1, cS_vfi.aD_new);
            if cS_vfi.aR_new < cS_vfi.aD_new, bV_payg_vfi((cS_vfi.aR_new+1):cS_vfi.aD_new) = M_vfi.b_t; end

            beq_transfer_vfi = M_vfi.beq_transfer_pers;

            % --- [核心修复：强制VFI调用路径一致] ---
            % 为了确保与转轨路径求解器100%的行为一致性，我们在此显式计算稳态
            % BGP增长因子，并使用与转轨求解器相同的字段名 `g_A_t_plus_1` 
            % 将其传递给底层的VFI函数。
            g_A_period_ss = (1 + cS_vfi.g_A_ss)^cS_vfi.time_Step - 1;
            M_vfi.g_A_t_plus_1 = g_A_period_ss;
            % 这一修改确保了无论是稳态还是转轨，HHSolutionByAge_VFI_* 函数
            % 内部都将执行相同的 `if isfield(...)` 逻辑分支，从而根除了
            % 任何潜在的、因代码路径分岔导致的行为不一致。
            % --- [修复结束] ---

            for a_idx = cS_vfi.aD_new : -1 : 1
                vPrime_kkppse_next = [];
                if a_idx < cS_vfi.aD_new, vPrime_kkppse_next = valS(:,:,:,a_idx+1); end

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

        % --- [v6.0 最终版] 本地矩阵化 NoPPS 函数 (定制化搜索空间) ---
        function [val_age, pol_age] = HHSolutionByAge_VFI_noPPS(a_idx, vPrime_kkppse_next, M_age, b_age_val, beq_transfer_val, paramS_age, cS)

            nk_search_grid = 150;
            val_age_slice = -1e20 * ones(cS.nk, 1, cS.nw_expanded);
            pol_age_slice = struct(...
                'c', zeros(cS.nk, 1, cS.nw_expanded), 'k_prime', zeros(cS.nk, 1, cS.nw_expanded), ...
                'kpps_prime', zeros(cS.nk, 1, cS.nw_expanded), 'tax_regular', zeros(cS.nk, 1, cS.nw_expanded), ...
                'tax_payg', zeros(cS.nk, 1, cS.nw_expanded), 'shock_exp', zeros(cS.nk, 1, cS.nw_expanded), ...
                'pension_out', zeros(cS.nk, 1, cS.nw_expanded), 'beq_received', zeros(cS.nk, 1, cS.nw_expanded));

            % 根据调用上下文（转轨或稳态）来确定增长率 g_A_period
            if isfield(M_age, 'g_A_t_plus_1')
                % 转轨调用: 使用传入的、随时间变化的下一期增长率
                g_A_period = M_age.g_A_t_plus_1;
            else
                % 稳态调用: 使用恒定的稳态增长率
                g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
            end
            market_return_factor = 1 + M_age.r_mkt_t;
            discount_factor_future_utility = cS.beta * (1 + g_A_period)^(1 - cS.sigma);
            k_grid_vec = cS.kGridV;

            if a_idx == cS.aD_new % 终期决策不变
                k_capital_tax_vec = cS.tau_k .* (k_grid_vec * M_age.r_mkt_t);
                total_wealth_vec = k_grid_vec * market_return_factor - k_capital_tax_vec + b_age_val + beq_transfer_val;
                k_prime_final_vec = zeros(cS.nk, 1);
                c_expend_final_vec = total_wealth_vec;
                if isfield(cS,'phi_bequest') && cS.phi_bequest > 0
                    phi_effective = cS.beta * cS.phi_bequest;
                    c_over_k_ratio = phi_effective^(-1/cS.sigma);
                    k_prime_final_vec = total_wealth_vec ./ (1 + (1+cS.tau_c)*c_over_k_ratio);
                    c_expend_final_vec = total_wealth_vec - k_prime_final_vec;
                end
                c_final_vec = c_expend_final_vec ./ (1 + cS.tau_c);
                consumption_tax_vec = c_final_vec * cS.tau_c;
                final_regular_tax_vec = k_capital_tax_vec + consumption_tax_vec;
                [~, util_c] = model_setup_utils_bgp.CES_utility(c_final_vec, cS.sigma, cS);
                util_bequest = model_setup_utils_bgp.bequest_utility(k_prime_final_vec, cS);
                util_final_vec = util_c + cS.beta * util_bequest;
                util_final_vec(total_wealth_vec < 0) = -1e20;

                for ie = 1:cS.nw_expanded
                    val_age_slice(:, 1, ie) = util_final_vec;
                    pol_age_slice.c(:, 1, ie) = c_final_vec; pol_age_slice.k_prime(:, 1, ie) = k_prime_final_vec;
                    pol_age_slice.pension_out(:, 1, ie) = b_age_val;
                    pol_age_slice.tax_regular(:, 1, ie) = final_regular_tax_vec;
                    pol_age_slice.beq_received(:, 1, ie) = beq_transfer_val;
                end
                val_age = repmat(val_age_slice, [1, cS.nkpps, 1]);
                pol_age = model_setup_utils_bgp.expand_policy_slice(pol_age_slice, cS.nkpps);
                return;
            end

            vPrime_interpolants = cell(cS.nw_expanded, 1);
            vPrime_slice = squeeze(vPrime_kkppse_next(:, 1, :));
            for ie_next = 1:cS.nw_expanded
                vPrime_interpolants{ie_next} = griddedInterpolant(k_grid_vec, vPrime_slice(:, ie_next), 'pchip', 'nearest');
            end
            trans_mat_next = paramS_age.TrProbM_by_age{a_idx + 1};

            for ie = 1:cS.nw_expanded
                k_return_vec = k_grid_vec * market_return_factor;
                labor_income_gross = 0; pension_out = 0;
                if a_idx <= cS.aR_new, labor_income_gross = M_age.w_t * cS.ageEffV_new(a_idx) * paramS_age.leGridV(ie);
                else, pension_out = b_age_val; end
                capital_tax_vec = cS.tau_k .* (k_grid_vec * M_age.r_mkt_t);
                payg_tax = cS.theta_t * labor_income_gross;
                labor_tax = cS.tau_l * max(0, labor_income_gross - payg_tax);
                cash_on_hand_vec = k_return_vec + labor_income_gross + pension_out + beq_transfer_val - (capital_tax_vec + payg_tax + labor_tax);
                shock_exp_factor = (ie == cS.nw + 1) * cS.kappa_young + (ie == cS.nw + 2) * cS.kappa_old;
                shock_exp_vec = shock_exp_factor .* cash_on_hand_vec;
                available_for_c_and_s_vec = cash_on_hand_vec - shock_exp_vec;

                % A. 创建定制化的 k' 搜索矩阵 [nk x nk_search_grid]
                k_prime_max_vec = (available_for_c_and_s_vec - cS.cFloor * (1+cS.tau_c)) / (1 + g_A_period);
                k_prime_max_vec(k_prime_max_vec < cS.kMin) = cS.kMin;

                linspace_weights = linspace(0, 1, nk_search_grid);
                k_prime_choices_mat = cS.kMin + (k_prime_max_vec - cS.kMin) .* linspace_weights;

                % B. 计算消费选择矩阵
                c_expend_choices_mat = available_for_c_and_s_vec - k_prime_choices_mat .* (1 + g_A_period);
                c_choices_mat = c_expend_choices_mat / (1 + cS.tau_c);

                % C. 计算价值矩阵
                invalid_choice_mask = (c_choices_mat < cS.cFloor);
                c_choices_mat(invalid_choice_mask) = cS.cFloor;
                [~, util_c_mat] = model_setup_utils_bgp.CES_utility(c_choices_mat, cS.sigma, cS);
                util_c_mat(invalid_choice_mask) = -1e20;

                % D. 在定制化网格上插值未来价值
                % Reshape for efficient interpolation
                k_prime_choices_flat = k_prime_choices_mat(:);
                v_prime_flat_mat = zeros(length(k_prime_choices_flat), cS.nw_expanded);
                for ie_next = 1:cS.nw_expanded
                    v_prime_flat_mat(:, ie_next) = vPrime_interpolants{ie_next}(k_prime_choices_flat);
                end
                ev_flat_vec = v_prime_flat_mat * trans_mat_next(ie, :)';
                ev_on_choices_mat = reshape(ev_flat_vec, cS.nk, nk_search_grid);

                util_bequest_mat = model_setup_utils_bgp.bequest_utility(k_prime_choices_mat, cS);

                future_value_mat = cS.s_pathV(a_idx) * ev_on_choices_mat + (1 - cS.s_pathV(a_idx)) * util_bequest_mat;

                total_value_mat = util_c_mat + discount_factor_future_utility * future_value_mat;

                [best_val_vec, best_idx_vec] = max(total_value_mat, [], 2);

                % E. 存储结果
                linear_indices = sub2ind(size(k_prime_choices_mat), (1:cS.nk)', best_idx_vec);
                k_prime_final_vec = k_prime_choices_mat(linear_indices);
                c_final_vec = c_choices_mat(linear_indices);

                val_age_slice(:, 1, ie) = best_val_vec;
                pol_age_slice.c(:, 1, ie) = c_final_vec;
                pol_age_slice.k_prime(:, 1, ie) = k_prime_final_vec;
                pol_age_slice.pension_out(:, 1, ie) = pension_out;
                pol_age_slice.shock_exp(:, 1, ie) = shock_exp_vec;
                pol_age_slice.tax_regular(:, 1, ie) = capital_tax_vec + labor_tax + c_final_vec .* cS.tau_c;
                pol_age_slice.tax_payg(:, 1, ie) = payg_tax;
                pol_age_slice.beq_received(:, 1, ie) = beq_transfer_val;
            end

            val_age = repmat(val_age_slice, [1, cS.nkpps, 1]);
            pol_age = model_setup_utils_bgp.expand_policy_slice(pol_age_slice, cS.nkpps);
        end

        % --- [v6.3 二维并行终极优化版] 本地矩阵化 PPS 函数 ---
        function [val_age, pol_age] = HHSolutionByAge_VFI_PPS(a_idx, vPrime_kkppse_next, M_age, b_age_val, beq_transfer_val, paramS_age, cS)
            % =========================================================================
            % == 函数: local_HHSolutionByAge_VFI_PPS
            % == 版本: [v6.3 - 【二维并行终极优化版】]
            % ==
            % == 核心修改:
            % ==   - [终极性能优化] 采用“扁平化”策略，将 (ie, ikpps) 两个维度的
            % ==     循环重构为一个单一的线性 parfor 循环。
            % ==   - 这使得总共 (nw_expanded * nkpps) 个独立任务都可以被并行调度，
            % ==     从而最大化地利用多核CPU资源，进一步提升计算速度。
            % =========================================================================

            nk_search_grid = 150;
            val_age = -1e20 * ones(cS.nk, cS.nkpps, cS.nw_expanded);
            pol_age = struct(...
                'c', zeros(cS.nk, cS.nkpps, cS.nw_expanded), 'k_prime', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
                'kpps_prime', zeros(cS.nk, cS.nkpps, cS.nw_expanded), 'tax_regular', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
                'tax_payg', zeros(cS.nk, cS.nkpps, cS.nw_expanded), 'shock_exp', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
                'pension_out', zeros(cS.nk, cS.nkpps, cS.nw_expanded), 'pps_contrib', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
                'pps_withdrawal', zeros(cS.nk, cS.nkpps, cS.nw_expanded), 'pps_tax', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
                'beq_received', zeros(cS.nk, cS.nkpps, cS.nw_expanded));

            % [!!! 关键修正 !!!]
            % 根据调用上下文（转轨或稳态）来确定增长率 g_A_period
            if isfield(M_age, 'g_A_t_plus_1')
                % 转轨调用: 使用传入的、随时间变化的下一期增长率
                g_A_period = M_age.g_A_t_plus_1;
            else
                % 稳态调用: 使用恒定的稳态增长率
                g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
            end
            market_return_factor = 1 + M_age.r_mkt_t;
            discount_factor_future_utility = cS.beta * (1 + g_A_period)^(1 - cS.sigma);
            k_grid_vec = cS.kGridV;

            if a_idx == cS.aD_new % 终期决策不变
                [K_grid_mat, Kpps_grid_mat] = ndgrid(k_grid_vec, cS.kppsGridV);
                k_capital_tax_mat = cS.tau_k .* (K_grid_mat .* M_age.r_mkt_t);
                pps_wealth_final_mat = Kpps_grid_mat .* market_return_factor;
                pps_tax_final_mat = cS.pps_tax_rate_withdrawal .* pps_wealth_final_mat;
                total_wealth_mat = (K_grid_mat * market_return_factor - k_capital_tax_mat) + (pps_wealth_final_mat - pps_tax_final_mat) + b_age_val + beq_transfer_val;
                k_prime_final_mat = zeros(size(K_grid_mat)); c_expend_final_mat = total_wealth_mat;
                if isfield(cS,'phi_bequest') && cS.phi_bequest > 0, phi_effective = cS.beta * cS.phi_bequest; c_over_k_ratio = phi_effective^(-1/cS.sigma); k_prime_final_mat = total_wealth_mat ./ (1 + (1+cS.tau_c)*c_over_k_ratio); c_expend_final_mat = total_wealth_mat - k_prime_final_mat; end
                c_final_mat = c_expend_final_mat ./ (1 + cS.tau_c);
                consumption_tax_mat = c_final_mat * cS.tau_c;
                [~, util_c] = model_setup_utils_bgp.CES_utility(c_final_mat, cS.sigma, cS);
                util_bequest = model_setup_utils_bgp.bequest_utility(k_prime_final_mat, cS);
                util_final_mat = util_c + cS.beta * util_bequest;
                util_final_mat(total_wealth_mat < 0) = -1e20;

                for ie = 1:cS.nw_expanded
                    val_age(:, :, ie) = util_final_mat;
                    pol_age.c(:, :, ie) = c_final_mat; pol_age.k_prime(:, :, ie) = k_prime_final_mat;
                    pol_age.kpps_prime(:, :, ie) = 0; pol_age.pension_out(:, :, ie) = b_age_val;
                    pol_age.tax_regular(:, :, ie) = k_capital_tax_mat + consumption_tax_mat;
                    pol_age.pps_withdrawal(:, :, ie) = pps_wealth_final_mat; pol_age.pps_tax(:, :, ie) = pps_tax_final_mat;
                    pol_age.beq_received(:, :, ie) = beq_transfer_val;
                end
                return;
            end

            vPrime_interpolants = cell(cS.nw_expanded, 1);
            for ie_next = 1:cS.nw_expanded
                vPrime_interpolants{ie_next} = griddedInterpolant({k_grid_vec, cS.kppsGridV}, vPrime_kkppse_next(:, :, ie_next), 'pchip', 'nearest');
            end
            trans_mat_next = paramS_age.TrProbM_by_age{a_idx + 1};

            % [核心修改] 创建一个一维的任务列表
            num_tasks = cS.nw_expanded * cS.nkpps;
            % 创建Cell数组来收集所有任务的结果
            val_results_cell = cell(num_tasks, 1);
            pol_results_cell = cell(num_tasks, 1);

            parfor task_idx = 1:num_tasks
                % 从一维任务索引映射回二维 (ie, ikpps) 索引
                [ikpps, ie] = ind2sub([cS.nkpps, cS.nw_expanded], task_idx);

                % --- 这部分代码是单次迭代的核心计算，逻辑基本不变 ---
                % --- 只是现在它不再在一个ikpps循环内部，而是直接执行 ---

                kpps_now = cS.kppsGridV(ikpps);
                pension_out=0; pps_contrib=0; pps_withdrawal=0; pps_tax=0; labor_income_gross = 0;
                if a_idx <= cS.aR_new
                    labor_income_gross = M_age.w_t * cS.ageEffV_new(a_idx) * paramS_age.leGridV(ie);
                    payg_tax = cS.theta_t * labor_income_gross;
                    pps_contrib = cS.pps_contrib_rate * labor_income_gross;
                    capital_tax_vec = cS.tau_k .* (k_grid_vec * M_age.r_mkt_t);
                    labor_tax = cS.tau_l * max(0, labor_income_gross - payg_tax - pps_contrib);
                    cash_on_hand_vec = (k_grid_vec * market_return_factor + labor_income_gross + beq_transfer_val) ...
                        - (capital_tax_vec + payg_tax + labor_tax + pps_contrib);
                else
                    pension_out = b_age_val;
                    kpps_total_value = kpps_now * market_return_factor;
                    pps_withdrawal = cS.pps_withdrawal_rate * kpps_total_value;
                    pps_tax = cS.pps_tax_rate_withdrawal * pps_withdrawal;
                    net_pps_inflow = pps_withdrawal - pps_tax;
                    capital_tax_vec = cS.tau_k .* (k_grid_vec * M_age.r_mkt_t);
                    payg_tax = 0; labor_tax = 0;
                    cash_on_hand_vec = (k_grid_vec * market_return_factor + pension_out + net_pps_inflow + beq_transfer_val) ...
                        - capital_tax_vec;
                end

                shock_exp_factor = (ie == cS.nw + 1) * cS.kappa_young + (ie == cS.nw + 2) * cS.kappa_old;
                shock_exp_vec = shock_exp_factor .* cash_on_hand_vec;
                available_for_c_and_s_vec = cash_on_hand_vec - shock_exp_vec;

                if a_idx <= cS.aR_new, kpps_prime_continuous = (kpps_now * market_return_factor + pps_contrib) / (1 + g_A_period);
                else, kpps_prime_continuous = (kpps_now * market_return_factor - pps_withdrawal) / (1 + g_A_period); end
                kpps_prime_final = max(cS.kppsMin, kpps_prime_continuous);

                k_prime_max_vec = (available_for_c_and_s_vec - cS.cFloor * (1+cS.tau_c)) / (1 + g_A_period);
                k_prime_max_vec(k_prime_max_vec < cS.kMin) = cS.kMin;
                linspace_weights = linspace(0, 1, nk_search_grid);
                k_prime_choices_mat = cS.kMin + (k_prime_max_vec - cS.kMin) .* linspace_weights;

                c_expend_choices_mat = available_for_c_and_s_vec - k_prime_choices_mat .* (1 + g_A_period);
                c_choices_mat = c_expend_choices_mat / (1 + cS.tau_c);

                invalid_choice_mask = (c_choices_mat < cS.cFloor);
                c_choices_mat(invalid_choice_mask) = cS.cFloor;
                [~, util_c_mat] = model_setup_utils_bgp.CES_utility(c_choices_mat, cS.sigma, cS);
                util_c_mat(invalid_choice_mask) = -1e20;

                k_prime_choices_flat = k_prime_choices_mat(:);
                kpps_prime_vec = ones(length(k_prime_choices_flat), 1) * kpps_prime_final;
                v_prime_flat_mat = zeros(length(k_prime_choices_flat), cS.nw_expanded);
                for ie_next = 1:cS.nw_expanded
                    v_prime_flat_mat(:, ie_next) = vPrime_interpolants{ie_next}(k_prime_choices_flat, kpps_prime_vec);
                end
                ev_flat_vec = v_prime_flat_mat * trans_mat_next(ie, :)';
                ev_on_choices_mat = reshape(ev_flat_vec, cS.nk, nk_search_grid);

                util_bequest_mat = model_setup_utils_bgp.bequest_utility(k_prime_choices_mat + kpps_prime_final, cS);

                future_value_mat = cS.s_pathV(a_idx) * ev_on_choices_mat + (1 - cS.s_pathV(a_idx)) * util_bequest_mat;
                total_value_mat = util_c_mat + discount_factor_future_utility * future_value_mat;

                [best_val_vec, best_idx_vec] = max(total_value_mat, [], 2);

                linear_indices = sub2ind(size(k_prime_choices_mat), (1:cS.nk)', best_idx_vec);
                k_prime_final_optimal_vec = k_prime_choices_mat(linear_indices);
                c_final_vec = c_choices_mat(linear_indices);

                % 将单次迭代的结果存入临时变量
                val_result_slice = best_val_vec;
                pol_result_slice = struct(...
                    'c', c_final_vec, 'k_prime', k_prime_final_optimal_vec, 'kpps_prime', kpps_prime_final, ...
                    'tax_regular', capital_tax_vec + labor_tax + c_final_vec * cS.tau_c, 'tax_payg', payg_tax, ...
                    'shock_exp', shock_exp_vec, 'pension_out', pension_out, 'pps_contrib', pps_contrib, ...
                    'pps_withdrawal', pps_withdrawal, 'pps_tax', pps_tax, 'beq_received', beq_transfer_val);

                % 将结果存入 Cell 数组中
                val_results_cell{task_idx} = val_result_slice;
                pol_results_cell{task_idx} = pol_result_slice;
            end

            % [核心修改] 在并行循环结束后，安全地将结果组装回最终的矩阵和结构体
            fields = fieldnames(pol_age);
            for task_idx = 1:num_tasks
                [ikpps, ie] = ind2sub([cS.nkpps, cS.nw_expanded], task_idx);

                % 组装价值函数
                val_age(:, ikpps, ie) = val_results_cell{task_idx};

                % 组装策略函数
                current_pol = pol_results_cell{task_idx};
                for i = 1:length(fields)
                    field = fields{i};
                    % 对于标量场，需要用 repmat 扩展维度
                    if isscalar(current_pol.(field))
                        pol_age.(field)(:, ikpps, ie) = repmat(current_pol.(field), cS.nk, 1);
                    else
                        pol_age.(field)(:, ikpps, ie) = current_pol.(field);
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

        % =======================================================
        % == 辅助函数:
        % =======================================================
        function errors = display_national_accounts_unified(ss, cS, Dist_to_verify, polS, paramS, verbose)
            % =========================================================================
            % == 函数: display_national_accounts_unified
            % == 版本: [v55 - 增加'已实现替代率'汇报]
            % ==
            % == 核心修改:
            % ==   1. [保留原貌] 函数打印到命令行的报告格式和内容与 v52 版兼容。
            % ==   2. [增加开关] 增加布尔型输入参数 `verbose` (默认为 true)。
            % ==   3. [增加输出] 函数总是返回一个包含核心误差的 `errors` 结构体。
            % ==   4. [新增指标] 在DB模式的报告中，增加汇报“已实现(毛)替代率”，
            % ==      该指标反映了在缴费率上限等约束下，实际达成的替代率水平，
            % ==      与政策目标形成直接对比。
            % =========================================================================

            if nargin < 6, verbose = true; end % 默认行为是打印报告

            % --- [第1步] 计算所有需要展示和返回的指标 ---
            g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
            n_period = (1 + cS.n_ss)^cS.time_Step - 1;
            g_total_period = (1 + g_A_period) * (1 + n_period) - 1;
            I_p_gross_demand = (g_total_period + cS.ddk) * ss.K_private_begin_hat;
            I_g_gross_demand = (g_total_period + cS.ddk_g) * ss.K_public_hat;
            I_total_gross = I_p_gross_demand + I_g_gross_demand;
            Total_Uses = ss.C_agg + I_total_gross + ss.G_c;

            err_Kp = ss.Saving_private_flow_Gross - I_p_gross_demand;
            err_Kg = ss.Saving_public_flow_Gross - I_g_gross_demand;
            err_L = 0; % by construction
            err_Beq = ss.Bequest_generated_agg - ss.Bequest_distributed_agg;
            err_NIPA = ss.Y_from_production_hat - Total_Uses;
            payg_balance = ss.Pension_in - ss.Pension_out;

            mass_retirees = sum(Dist_to_verify(:,:,:,(cS.aR_new+1):end), 'all');

            % a) 计算替代率指标
            avg_gross_wage = ss.w_hat; % 平均毛工资用于计算毛替代率
            net_worker_wage = ss.w_hat * (1 - ss.theta) * (1 - cS.tau_l); % 平均净工资用于计算净替代率
            avg_pension_actual = 0;
            if mass_retirees > 1e-9
                avg_pension_actual = ss.Pension_out / mass_retirees;
            end

            realized_gross_rr = 0;
            if avg_gross_wage > 1e-9
                realized_gross_rr = avg_pension_actual / avg_gross_wage;
            end

            net_payg_rr = 0;
            if net_worker_wage > 1e-9
                net_payg_rr = avg_pension_actual / net_worker_wage;
            end

            net_pps_rr = 0;
            if cS.pps_active && net_worker_wage > 1e-9 && mass_retirees > 1e-9
                avg_pps_net_withdrawal = (ss.PPS_withdrawal_agg - ss.PPS_tax_agg) / mass_retirees;
                net_pps_rr = avg_pps_net_withdrawal / net_worker_wage;
            end
            total_net_rr = net_payg_rr + net_pps_rr;


            % --- [第2步] 将需要返回的误差打包 ---
            errors = struct(...
                'err_Kp', err_Kp, ...
                'err_Kg', err_Kg, ...
                'err_L', err_L, ...
                'err_Beq', err_Beq, ...
                'payg_balance', payg_balance, ...
                'err_NIPA', err_NIPA ...
                );

            % --- [第3步] 如果 verbose 为 true, 则打印完整报告 ---
            if verbose
                fprintf('\n\n================================================================================\n');
                title_str = '国民经济核算与自洽性检验报告 (全诊断版)';
                fprintf('===%s===\n', pad(title_str, 75, 'both'));

                if isfield(cS, 'endogenous_theta_mode') && cS.endogenous_theta_mode
                    fprintf('--- 养老金模式: [内生缴费率] (Defined Benefit / DB)\n');
                    fprintf('   政策目标: (毛)替代率 ................... = %.2f %%\n', cS.payg_replacement_rate * 100);
                    % --- [新增汇报行] ---
                    fprintf('   已实现: (毛)替代率 ................... = %.2f %%\n', realized_gross_rr * 100);
                    % --- [新增结束] ---
                    if isfield(cS, 'theta_max'), fprintf('   政策约束: 最大缴费率上限 ............. = %.2f %%\n', cS.theta_max * 100); end
                else
                    fprintf('--- 养老金模式: [外生缴费率] (Defined Contribution / DC)\n');
                end
                fprintf('   政策设定: 最终缴费率 (Theta) ......... = %.4f (%.2f %%)\n', ss.theta, ss.theta * 100);

                fprintf('   ----------------------------------------------------------------------------\n');
                fprintf('   --- [核心指标] 净替代率 (基于税后工资 w_net=%.4f) ---\n', net_worker_wage);
                fprintf('   公共养老金(PAYG) 贡献的净替代率 ....... = %.2f %%\n', net_payg_rr * 100);
                if cS.pps_active
                    fprintf('   私人养老金(PPS) 贡献的净替代率 ....... = %.2f %%\n', net_pps_rr * 100);
                    fprintf('   ----------------------------------------------------------------------------\n');
                end
                fprintf('   => 总计净替代率 (衡量生活水平) ......... = %.2f %%\n', total_net_rr * 100);

                fprintf('================================================================================\n');

                % --- 诊断 1: 市场出清检验 ---
                fprintf('\n--- [诊断 1] 求解器目标市场出清检验 ---\n');
                fprintf('   市场1: 私人资本 S_p - I_p ...... : %12.8f - %12.8f = %12.4e\n', ss.Saving_private_flow_Gross, I_p_gross_demand, err_Kp);
                fprintf('   市场2: 政府资本 S_g - I_g ...... : %12.8f - %12.8f = %12.4e\n', ss.Saving_public_flow_Gross, I_g_gross_demand, err_Kg);
                fprintf('   市场3: 劳动供给-需求 .......... : %12.8f - %12.8f = %12.4e\n', ss.L_hat, ss.L_hat, err_L);
                fprintf('   市场4: 遗赠产生-分配 .......... : %12.8f - %12.8f = %12.4e\n', ss.Bequest_generated_agg, ss.Bequest_distributed_agg, err_Beq);
                fprintf('   PAYG系统: 收入 - 支出 .......... : %12.8f - %12.8f = %12.4e\n', ss.Pension_in, ss.Pension_out, payg_balance);

                if payg_balance < -1e-9 && payg_balance / ss.Y_from_production_hat < -1e-5 % 增加一个相对重要的阈值
                    payg_shortfall = -payg_balance;
                    non_pension_revenue = ss.Regular_tax + ss.PPS_tax_agg + ss.Public_Capital_Return;
                    if non_pension_revenue > 1e-9
                        crowding_out_ratio = payg_shortfall / non_pension_revenue;
                        fprintf('      => 该缺口相当于政府非养老金总收入的 %.2f %%，造成等额挤占。\n', crowding_out_ratio * 100);
                    else
                        fprintf('      => 存在缺口(%.4f)，但政府非养老金收入为零，无法计算挤占比例。\n', payg_shortfall);
                    end
                end
            end
        end

        function [x_eq, eq_found] = solve_with_fsolve(system_wrapper, x0, verbose)
            % [保留不变] fsolve求解器的包装函数
            if verbose, fsolve_display = 'iter'; else, fsolve_display = 'none'; end
            options = optimoptions('fsolve', 'Display', fsolve_display, 'TolFun', 1e-8, 'TolX', 1e-8, 'MaxIterations', 100);

            if verbose, fprintf('\n--- 启动 fsolve 求解器 ---\n'); end
            [x_eq, ~, exitflag] = fsolve(system_wrapper, x0, options);
            if verbose, fprintf('--- fsolve 求解完成 ---\n'); end

            eq_found = (exitflag > 0);
        end

        function [x_eq, eq_found] = solve_with_lsqnonlin(system_wrapper, x0, lb, ub, verbose)
            % [保留不变] lsqnonlin 求解器的包装函数
            if verbose, lsq_display = 'iter'; else, lsq_display = 'none'; end
            options = optimoptions('lsqnonlin', 'Display', lsq_display, 'FunctionTolerance', 1e-20, 'StepTolerance', 1e-15, 'MaxIterations', 200, 'Algorithm', 'trust-region-reflective');

            if verbose, fprintf('\n--- 启动 lsqnonlin 求解器 (带边界约束) ---\n'); end
            [x_eq, resnorm, ~, exitflag] = lsqnonlin(system_wrapper, x0, lb, ub, options);
            if verbose, fprintf('--- lsqnonlin 求解完成, 残差平方和(resnorm): %.4e ---\n', resnorm); end

            eq_found = (exitflag > 0);
        end

        % 这是 main_steady_state_utils_bgp.m 中需要修改的单个函数
        % 其他函数无需改动

        function M_prices = get_prices_at_t(K_p_hat, K_g_hat, L_hat, A_t, cS)
            % =========================================================================
            % == 函数: get_prices_at_t
            % == 版本: [v4 - 全去趋势化输入/输出校正版]
            % ==
            % == 核心修正:
            % ==   1. [输入统一] 输入变量被明确为去趋势化的 "hat" 单位
            % ==      (K_p_hat, K_g_hat, L_hat)。这与转轨函数完全一致。
            % ==   2. [输出统一] 函数返回去趋势化的 Y_hat 和 w_hat。
            % ==   3. [逻辑简化] A_t 参数实际上不再需要，因为所有计算都在
            % ==      去趋势化的空间内进行。保留它是为了接口兼容，但在函数内部
            % ==      它被忽略，因为 A_hat=1。
            % ==   4. [命名修正] 输出结构体字段名被修正为 Y_hat 和 w_hat，
            % ==      消除歧义。
            % =========================================================================
            if K_p_hat <= 0, K_p_hat = 1e-8; end; if L_hat <= 0, L_hat = 1e-8; end; if K_g_hat <= 0, K_g_hat = 1e-8; end;

            labor_exponent = 1 - cS.alpha - cS.gamma;
            if labor_exponent <= 0, error('劳动和资本的产出弹性之和必须小于1！'); end

            % 在去趋势化的世界里，A_hat 总是 1.0
            A_hat = 1.0;
            
            % 计算去趋势化的产出 Y_hat
            Y_hat_t = A_hat .* (K_p_hat.^cS.alpha) .* (K_g_hat.^cS.gamma) .* (L_hat.^labor_exponent);

            % 根据边际产出计算要素价格 (所有价格都是 "hat" 或稳态单位)
            MPK_p_period = cS.alpha .* Y_hat_t ./ K_p_hat;
            w_hat_t = labor_exponent .* Y_hat_t ./ L_hat; 

            % 市场利率 r_mkt_t 是平稳的，不受 A_t 直接影响
            r_mkt_t = MPK_p_period - cS.ddk;

            % 返回的结构体字段名也明确为 hat
            M_prices = struct('K_p_hat', K_p_hat, 'K_g_hat', K_g_hat, 'L_t_hat', L_hat, 'Y_hat_t', Y_hat_t, 'w_hat_t', w_hat_t, 'r_mkt_t', r_mkt_t);
        end

        


        function [idx_lower, idx_upper, w_upper] = find_grid_and_weights(value, gridV)
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

