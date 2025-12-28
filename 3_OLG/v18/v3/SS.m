% --- SS.m ---
classdef SS
    methods (Static)
        function [ss, Dist, polS, valS] = solve_steady_state(cS, paramS, params_ext, verbose, solver_method)
            % =========================================================================
            % == 函数: solve_steady_state
            % == 版本: [v1.2 - 顶层计算最终版]
            % ==
            % == 核心修改:
            % ==   - 将所有宏观计算(Y, K_g, 价格, 国民账户)都提升到这个顶层函数中。
            % ==   - 确保在找到均衡解后，使用完全相同的计算逻辑生成最终的ss结构体。
            % =========================================================================

            if nargin < 4, verbose = true; end
            if nargin < 5, solver_method = 'lsqnonlin'; end

            cS.A = params_ext.A;
            if isfield(params_ext, 'g_A_ss'), cS.g_A_ss = params_ext.g_A_ss; end
            if isfield(params_ext, 'n_ss'), cS.n_ss = params_ext.n_ss; end
            if isfield(params_ext, 'Z') && ~isempty(params_ext.Z), Z_ss_norm = params_ext.Z;
            else, error('错误：求解稳态需要一个明确的人口分布Z。'); end

            % 预计算L_supply_exog (这部分正确，无需修改)
            expected_le = sum(paramS.leGridV .* paramS.leProb1V);
            L_supply_exog = 0;
            for ia = 1:cS.aR_new
                L_supply_exog = L_supply_exog + (cS.ageEffV_new(ia) * expected_le * Z_ss_norm(ia));
            end
            fprintf('   [预计算] 外生总有效劳动供给 (L_supply_exog): %.6f\n', L_supply_exog);

            % 设置初始猜测值和求解器
            r_guess_init = 0.04;
            w_guess_init = 1.5;
            K_guess_for_BQ = 3.0 * (cS.nw / 1.5);
            BQ_guess_init = mean(1-cS.s_pathV) * K_guess_for_BQ / sum(Z_ss_norm(1:cS.aD_new));
            TR_guess_init = 0.0;
            x0 = [r_guess_init, w_guess_init, BQ_guess_init, TR_guess_init];
            lb = [1e-8, 1e-8, 1e-8, -Inf];
            ub = [Inf, Inf, Inf, Inf];
            system_wrapper = @(x) SS.system_of_equations(x, Z_ss_norm, L_supply_exog, cS, paramS, params_ext);
            options = optimoptions('lsqnonlin', 'Display', 'iter', 'FunctionTolerance', 1e-12, 'StepTolerance', 1e-12, 'MaxIterations', 200);
            [x_eq, ~, ~, exitflag] = lsqnonlin(system_wrapper, x0, lb, ub, options);

            if (exitflag <= 0), error('求解器未能找到均衡解！'); end
            fprintf('\n--- ✅ 求解器收敛！均衡解 [r, w, BQ, TR]: [%.6f, %.6f, %.6f, %.6f] ---\n', x_eq(1), x_eq(2), x_eq(3), x_eq(4));

            % --- [修改] 使用均衡解，并用与求解器内部完全相同的逻辑生成最终ss ---
            [~, ss, Dist, polS, valS] = SS.system_of_equations(x_eq, Z_ss_norm, L_supply_exog, cS, paramS, params_ext);
        end

        function [F_error, ss, Dist, polS, valS] = system_of_equations(x_guess, Z_ss_norm, L_supply_exog, cS, paramS, params_ext)
            % =========================================================================
            % == 函数: system_of_equations
            % == 版本: [v1.3 - PAYG预算分离版]
            % ==
            % == 核心修改:
            % ==   - 在求解TR时，明确将PAYG缴费从政府一般收入中剔除，
            % ==     确保求解器和报告函数的政府预算定义完全一致。
            % =========================================================================
            r_guess = x_guess(1);
            w_guess = x_guess(2);
            BQ_guess = x_guess(3);
            TR_guess = x_guess(4);

            % --- 流程 1: 微观聚合 ---
            [macro_state, Dist, polS, valS] = SS.run_micro_to_macro(r_guess, w_guess, BQ_guess, TR_guess, Z_ss_norm, L_supply_exog, cS, paramS, params_ext);

            % --- 流程 2: 宏观计算 ---
            g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
            n_period = (1 + cS.n_ss)^cS.time_Step - 1;
            g_total_period = (1 + g_A_period) * (1 + n_period) - 1;
            Kg_Y_ratio_factor = cS.I_g_to_Y_ratio_ss / (g_total_period + cS.ddk_g);
            term1 = macro_state.K_private_hat^cS.alpha;
            term2 = macro_state.L_hat^(1 - cS.alpha - cS.gamma);
            term3 = Kg_Y_ratio_factor^cS.gamma;
            Y_hat_calc = (term1 * term2 * term3)^(1 / (1 - cS.gamma));
            K_g_hat_calc = Kg_Y_ratio_factor * Y_hat_calc;

            prices_new = firm.get_prices_at_t(macro_state.K_private_hat, K_g_hat_calc, macro_state.L_hat, cS);
            r_new = prices_new.r_mkt_t;
            w_new = prices_new.w_hat_t;

            BQ_new = macro_state.Bequest_generated_agg;
            
            % [核心修正] 从政府预算约束计算引致的转移支付
            % a. 计算PAYG缴费总额
            PAYG_contributions_agg = cS.theta_path(end) * w_new * macro_state.L_hat;
            
            % b. 计算一般税收 (从总聚合税收中剔除PAYG缴费)
            General_Tax_Revenue_calc = macro_state.Total_Tax_Revenue - PAYG_contributions_agg;
            
            % c. 计算公共资本回报
            Public_Capital_Return_calc = prices_new.Y_hat_t - (w_new * macro_state.L_hat) - ((r_new + cS.ddk) * macro_state.K_private_hat);
            
            % d. 计算政府一般收入
            Total_General_Gov_Revenue_calc = General_Tax_Revenue_calc + Public_Capital_Return_calc;
            
            % e. 计算政府支出
            G_consumption_spending = cS.G_c_to_Y_ratio_ss * prices_new.Y_hat_t;
            G_investment_spending = cS.I_g_to_Y_ratio_ss * prices_new.Y_hat_t;
            
            % f. 求解使得一般预算平衡的 TR
            TR_new = Total_General_Gov_Revenue_calc - (G_consumption_spending + G_investment_spending);

            % --- 流程 3: 计算误差 ---
            error_r = r_guess - r_new;
            error_w = w_guess - w_new;
            error_BQ = BQ_guess - BQ_new;
            error_TR = TR_guess - TR_new;
            F_error = [error_r; error_w; error_BQ; error_TR];

            % --- 流程 4: 打包最终的ss结构体 ---
            if nargout > 1
                ss = struct();
                ss.Y_from_production_hat = prices_new.Y_hat_t;
                ss.K_private_hat = macro_state.K_private_hat;
                ss.K_public_hat = K_g_hat_calc;
                ss.L_hat = macro_state.L_hat;
                ss.r_mkt = r_new;
                ss.w_hat = w_new;
                ss.C_agg = macro_state.C_agg;
                ss.I_g = G_investment_spending;
                ss.Total_Tax_Revenue = macro_state.Total_Tax_Revenue; % 仍然保存包含PAYG的总税收，以便报告函数分解
                ss.Bequest_generated_agg = BQ_new;
                ss.Bequest_distributed_agg = BQ_new;
                ss.TR_distributed_agg = TR_new;
                ss.Public_Capital_Return = Public_Capital_Return_calc;
                ss.G_c = G_consumption_spending;
            end
        end

        function [macro_state, Dist, polS, valS] = run_micro_to_macro(r_mkt, w_hat, BQ, TR, Z_ss_norm, L_supply_exog, cS, paramS, params_ext)
            % =========================================================================
            % == 函数: run_micro_to_macro
            % == 版本: [v1.3 - 养老金缴费(theta)兼容版]
            % ==
            % == 核心修改:
            % ==   - 将稳态的养老金缴费率 theta_ss 添加到传递给VFI的结构体中。
            % =========================================================================

            % --- 1. 计算稳态的标准化养老金(b_hat) ---
            total_pension_pot_hat = cS.theta_path(end) * w_hat * L_supply_exog;
            mass_retirees = sum(Z_ss_norm((cS.aR_new+1):end));
            b_hat = total_pension_pot_hat / max(1e-9, mass_retirees);

            % --- 2. 给定价格、养老金福利和缴费率，解家庭问题 ---
            M_for_hh.r_mkt_t = r_mkt;
            M_for_hh.w_t = w_hat;
            M_for_hh.tr_per_hh = TR / sum(Z_ss_norm);
            M_for_hh.b_hat_t = b_hat;
            M_for_hh.theta_t = cS.theta_path(end); % [核心修正] 传递稳态缴费率

            [polS, valS] = household.VFI_solver(M_for_hh, paramS, cS);

            % --- 3. 给定策略，解稳态分布 ---
            Dist = distribution.get_ss_dist(polS, paramS, cS, Z_ss_norm, BQ);

            % --- 4. 聚合出宏观状态 ---
            [K_private_hat_agg, C_agg, Bequest_generated_agg, Total_Tax_Revenue] = aggregates.get_aggregates(Dist, polS, cS, paramS);

            % --- 5. 将聚合结果打包返回 ---
            macro_state = struct();
            macro_state.K_private_hat = K_private_hat_agg;
            macro_state.C_agg = C_agg;
            macro_state.Bequest_generated_agg = Bequest_generated_agg;
            macro_state.Total_Tax_Revenue = Total_Tax_Revenue;
            macro_state.L_hat = L_supply_exog;
        end


        function display_national_accounts(ss, cS, paramS, Z_ss_norm, report_filename)
            % =========================================================================
            % == 函数: display_national_accounts
            % == 版本: [v7.7 - PAYG 核算分离版]
            % ==
            % == 核心功能:
            % ==   - 在报告中明确分离 PAYG 系统的收入（缴费）和支出（福利）。
            % ==   - 政府的一般收入和支出不再包含 PAYG 流量。
            % ==   - 这使得政府预算平衡的检验更加清晰。
            % =========================================================================

            % --- 0. 准备工作：文件句柄和BGP增长因子 ---
            try
                fileID = fopen(report_filename, 'w', 'n', 'UTF-8');
                if fileID == -1, error('无法创建或打开报告文件: %s', report_filename); end
            catch ME
                warning('创建报告文件时出错，将输出到命令行。错误');
                fileID = 1;
            end

            g_A_period = (1+cS.g_A_ss)^cS.time_Step-1;
            n_period = (1+cS.n_ss)^cS.time_Step-1;
            g_total_period = (1+g_A_period)*(1+n_period)-1;
            Y = ss.Y_from_production_hat;

            % --- 1. [核心修正] 重新计算 PAYG 流量和真实税收 ---
            Dist = SS.get_final_dist(ss, cS, paramS, Z_ss_norm);
            polS = SS.get_final_polS(ss, cS, paramS, Z_ss_norm);

            PAYG_contributions_agg = 0;
            PAYG_benefits_agg = 0;
            General_Tax_Revenue = 0; % 税收，不含PAYG缴费

            % 重新聚合PPS相关的流量 (如果激活)
            pps_active = isfield(cS, 'pps_active') && cS.pps_active;
            K_k_hat_agg = 0; K_pps_hat_agg = 0;
            PPS_contrib_agg = 0; PPS_withdrawal_agg = 0; PPS_tax_agg = 0;

            for ia = 1:cS.aD_new
                mass_dist_ia = Dist(:,:,:,ia);
                
                % PAYG 流量计算
                labor_income_slice_full = ss.w_hat * cS.ageEffV_new(ia) .* paramS.leGridV(1:cS.nw_expanded)';
                payg_contrib_slice_full = cS.theta_path(end) * labor_income_slice_full;
                payg_contrib_grid = repmat(reshape(payg_contrib_slice_full, [1,1,cS.nw_expanded]), [cS.nk, cS.nkpps, 1]);

                if ia <= cS.aR_new
                    PAYG_contributions_agg = PAYG_contributions_agg + sum(payg_contrib_grid .* mass_dist_ia, 'all');
                else
                    total_pension_pot_hat = cS.theta_path(end) * ss.w_hat * ss.L_hat;
                    mass_retirees = sum(Z_ss_norm((cS.aR_new+1):end));
                    b_hat_ss = total_pension_pot_hat / max(1e-9, mass_retirees);
                    PAYG_benefits_agg = PAYG_benefits_agg + sum(b_hat_ss .* mass_dist_ia, 'all');
                end
                
                % 从聚合的 tax_regular 中减去 PAYG 缴费，得到纯粹的税收
                tax_regular_slice = polS(ia).tax_regular;
                general_tax_slice = tax_regular_slice - payg_contrib_grid;
                General_Tax_Revenue = General_Tax_Revenue + sum(general_tax_slice .* mass_dist_ia, 'all');
                
                % PPS 流量和存量计算 (如果激活)
                if pps_active
                    kpps_prime_slice = polS(ia).kpps_prime;
                    K_pps_hat_agg = K_pps_hat_agg + sum(kpps_prime_slice .* mass_dist_ia, 'all');

                    if ia <= cS.aR_new
                        pps_contrib_slice_full = cS.pps_contrib_rate * labor_income_slice_full;
                        pps_contrib_grid = repmat(reshape(pps_contrib_slice_full, [1,1,cS.nw_expanded]), [cS.nk, cS.nkpps, 1]);
                        PPS_contrib_agg = PPS_contrib_agg + sum(pps_contrib_grid .* mass_dist_ia, 'all');
                    else
                        kpps_current_grid = cS.kppsGridV;
                        kpps_current_full = repmat(reshape(kpps_current_grid, [1, cS.nkpps, 1]), [cS.nk, 1, cS.nw_expanded]);
                        kpps_level = kpps_current_full * (1 + ss.r_mkt) * (1 + g_A_period);
                        pps_withdrawal_slice = cS.pps_withdrawal_rate * kpps_level;
                        pps_tax_slice = cS.pps_tax_rate_withdrawal * pps_withdrawal_slice;
                        PPS_withdrawal_agg = PPS_withdrawal_agg + sum(pps_withdrawal_slice .* mass_dist_ia, 'all');
                        PPS_tax_agg = PPS_tax_agg + sum(pps_tax_slice .* mass_dist_ia, 'all');
                    end
                end
                K_k_hat_agg = K_k_hat_agg + sum(polS(ia).k_prime .* mass_dist_ia, 'all');
            end
            K_k_hat_agg = K_k_hat_agg / (1 + n_period);
            K_pps_hat_agg = K_pps_hat_agg / (1 + n_period);
            General_Tax_Revenue = General_Tax_Revenue - PPS_tax_agg; % 确保一般税收不含PPS税

            % --- 2. 准备国民账户组件 ---
            I_p_demand_bgp = (g_total_period + cS.ddk) * ss.K_private_hat;
            I_g_display = cS.I_g_to_Y_ratio_ss * Y;
            G_c_display = cS.G_c_to_Y_ratio_ss * Y;
            I_p_accounting = Y - ss.C_agg - I_g_display - G_c_display;
            Total_Expenditure = ss.C_agg + I_p_accounting + I_g_display + G_c_display;
            resource_constraint_error = Y - Total_Expenditure;

            Depreciation_p = cS.ddk * ss.K_private_hat;
            Depreciation_g = cS.ddk_g * ss.K_public_hat;
            Total_Labor_Income = ss.w_hat * ss.L_hat;
            Total_Capital_Income_Gross = (ss.r_mkt + cS.ddk) * ss.K_private_hat;
            Public_Capital_Return_display = Y - Total_Labor_Income - Total_Capital_Income_Gross;
            
            Total_General_Gov_Revenue = General_Tax_Revenue + PPS_tax_agg + Public_Capital_Return_display;
            Total_General_Gov_Expenditure = G_c_display + I_g_display + ss.TR_distributed_agg;
            general_budget_balance_error = Total_General_Gov_Revenue - Total_General_Gov_Expenditure;
            payg_balance_error = PAYG_contributions_agg - PAYG_benefits_agg;

            % --- 报告开始 ---
            fprintfAndLog(fileID, '\n\n=========================================================================================================\n');
            if pps_active
                fprintfAndLog(fileID, '###                 稳 态 国 民 经 济 核 算 报 告 (v7.7 - 带PPS, PAYG分离版)              ###\n');
            else
                fprintfAndLog(fileID, '###                       稳 态 国 民 经 济 核 算 报 告 (v7.7 - PAYG分离版)                      ###\n');
            end
            fprintfAndLog(fileID, '=========================================================================================================\n');

            fprintfAndLog(fileID, '\n--- [ I. 宏观总量与国民支出 ] ---\n');
            fprintfAndLog(fileID, '   A. 生产 (GDP) 与核心比率:\n');
            fprintfAndLog(fileID, '      国内生产总值 (Y) .......................... : %15.6f\n', Y);
            fprintfAndLog(fileID, '      私人资本/产出比 (K_p/Y) ................... : %15.6f\n', ss.K_private_hat / Y);
            fprintfAndLog(fileID, '      公共资本/产出比 (K_g/Y) ................... : %15.6f\n', ss.K_public_hat / Y);
            fprintfAndLog(fileID, '   B. 均衡价格:\n');
            fprintfAndLog(fileID, '      真实利率 (r, 模型期) ...................... : %15.6f (年化: %.4f %%)\n', ss.r_mkt, ((1+ss.r_mkt)^(1/cS.time_Step) - 1)*100);
            fprintfAndLog(fileID, '      真实工资率 (w, 单位有效劳动) .............. : %15.6f\n', ss.w_hat);

            fprintfAndLog(fileID, '   C. 国民支出 (Y = C + I_p + I_g + G_c) 与资源约束检验:\n');
            fprintfAndLog(fileID, '      (+) 私人消费 (C) .......................... : %15.6f  (%6.2f %% of Y)\n', ss.C_agg, (ss.C_agg / Y)*100);
            fprintfAndLog(fileID, '      (+) 私人总投资 (I_p, 会计值) ............ : %15.6f  (%6.2f %% of Y)\n', I_p_accounting, (I_p_accounting / Y)*100);
            fprintfAndLog(fileID, '      (+) 公共总投资 (I_g) ...................... : %15.6f  (%6.2f %% of Y)\n', I_g_display, (I_g_display / Y)*100);
            fprintfAndLog(fileID, '      (+) 政府消费 (G_c) ........................ : %15.6f  (%6.2f %% of Y)\n', G_c_display, (G_c_display / Y)*100);
            fprintfAndLog(fileID, '      --------------------------------------------------------------------------------------------------\n');
            fprintfAndLog(fileID, '      (=) 总支出 (C + I_p_acct + I_g + G_c) ..... : %15.6f\n', Total_Expenditure);
            fprintfAndLog(fileID, '      --------------------------------------------------------------------------------------------------\n');
            fprintfAndLog(fileID, '      >>> 资源约束误差 (Y - 总支出) ............. : %15.4e (Y的 %.4e %%)\n', resource_constraint_error, (resource_constraint_error / Y)*100);
            if abs(resource_constraint_error/Y) < 1e-6
                fprintfAndLog(fileID, '          核验结果: ✅ 资源约束在会计上严格满足。\n');
            else
                fprintfAndLog(fileID, '          核验结果: ⚠️ 资源约束在会计上存在误差！\n');
            end

            fprintfAndLog(fileID, '\n   D. 投资一致性检验 (仅适用于真稳态):\n');
            fprintfAndLog(fileID, '      BGP理论投资需求 (g+d)K .................. : %15.6f\n', I_p_demand_bgp);
            fprintfAndLog(fileID, '      会计投资 (Y-C-I_g-G_c) .................... : %15.6f\n', I_p_accounting);
            investment_gap = I_p_accounting - I_p_demand_bgp;
            fprintfAndLog(fileID, '      >>> 投资缺口 (会计值 - 理论值) .......... : %15.4e (Y的 %.4f %%)\n', investment_gap, (investment_gap/Y)*100);
            if abs(investment_gap/Y) < 1e-4
                fprintfAndLog(fileID, '          检验结果: ✅ 经济体处于或接近真实的BGP稳态。\n');
            else
                fprintfAndLog(fileID, '          检验结果: ⚠️  经济体偏离BGP稳态 (预期中的伪稳态特征)。\n');
            end

            % --- [ II. 厂商部门 ] ---
            fprintfAndLog(fileID, '\n\n--- [ II. 厂商部门 ] ---\n');
            fprintfAndLog(fileID, '   A. 生产要素投入 (存量):\n');
            fprintfAndLog(fileID, '      总私人资本存量 (K_p) ...................... : %15.6f\n', ss.K_private_hat);
            if pps_active
                fprintfAndLog(fileID, '         其中: 常规资本 (K_k) ................. : %15.6f  (%6.2f %% of K_p)\n', K_k_hat_agg, (K_k_hat_agg / ss.K_private_hat) * 100);
                fprintfAndLog(fileID, '         其中: PPS资本 (K_pps) ................ : %15.6f  (%6.2f %% of K_p)\n', K_pps_hat_agg, (K_pps_hat_agg / ss.K_private_hat) * 100);
            end
            fprintfAndLog(fileID, '      公共资本存量 (K_g) ........................ : %15.6f\n', ss.K_public_hat);
            fprintfAndLog(fileID, '      有效劳动需求 (L) .......................... : %15.6f\n', ss.L_hat);
            fprintfAndLog(fileID, '   B. 要素报酬 (流量):\n');
            fprintfAndLog(fileID, '      支付工资总额 (w*L) ........................ : %15.6f  (Y的 %6.2f %%)\n', Total_Labor_Income, (Total_Labor_Income/Y)*100);
            fprintfAndLog(fileID, '      私人资本总回报 ((r+d)*K_p) ................ : %15.6f  (Y的 %6.2f %%)\n', Total_Capital_Income_Gross, (Total_Capital_Income_Gross/Y)*100);
            fprintfAndLog(fileID, '      公共资本回报 (Y - wL - (r+d)K_p) ........ : %15.6f  (Y的 %6.2f %%)\n', Public_Capital_Return_display, (Public_Capital_Return_display/Y)*100);
            fprintfAndLog(fileID, '   C. 投资与折旧 (流量):\n');
            fprintfAndLog(fileID, '      私人总投资 (I_p, BGP理论值) ............... : %15.6f\n', I_p_demand_bgp);
            fprintfAndLog(fileID, '      私人资本折旧 (d*K_p) ...................... : %15.6f\n', Depreciation_p);
            fprintfAndLog(fileID, '      私人净投资 (I_p_bgp - d*K_p) .............. : %15.6f\n', I_p_demand_bgp - Depreciation_p);

            % --- [ III. 家庭部门 ] ---
            fprintfAndLog(fileID, '\n\n--- [ III. 家庭部门 ] ---\n');
            fprintfAndLog(fileID, '   A. 收入来源 (流量):\n');
            fprintfAndLog(fileID, '      劳动收入 (来自厂商) ....................... : %15.6f\n', Total_Labor_Income);
            Net_Capital_Income = ss.r_mkt * ss.K_private_hat;
            fprintfAndLog(fileID, '      净资本收入 (r*K_p) ........................ : %15.6f\n', Net_Capital_Income);
            fprintfAndLog(fileID, '      养老金福利 ................................ : %15.6f\n', PAYG_benefits_agg);
            if pps_active
                fprintfAndLog(fileID, '      PPS账户提取 ............................... : %15.6f\n', PPS_withdrawal_agg);
            end
            fprintfAndLog(fileID, '      收到的遗赠 ................................ : %15.6f\n', ss.Bequest_distributed_agg);
            fprintfAndLog(fileID, '      收到的转移支付 ............................ : %15.6f\n', ss.TR_distributed_agg);
            Total_Household_Income = Total_Labor_Income + Net_Capital_Income + PAYG_benefits_agg + PPS_withdrawal_agg + ss.Bequest_distributed_agg + ss.TR_distributed_agg;
            fprintfAndLog(fileID, '      --------------------------------------------------------------------------------------------------\n');
            fprintfAndLog(fileID, '      >>> 家庭总收入 ............................ : %15.6f\n', Total_Household_Income);
            fprintfAndLog(fileID, '   B. 收入用途 (流量):\n');
            fprintfAndLog(fileID, '      消费支出 (C) .............................. : %15.6f\n', ss.C_agg);
            fprintfAndLog(fileID, '      支付的税收 (资本税/劳动税/消费税) ......... : %15.6f\n', General_Tax_Revenue);
            fprintfAndLog(fileID, '      养老金缴费 ................................ : %15.6f\n', PAYG_contributions_agg);
            if pps_active
                fprintfAndLog(fileID, '      PPS账户缴费 (储蓄) ...................... : %15.6f\n', PPS_contrib_agg);
                fprintfAndLog(fileID, '      PPS提取税 ................................. : %15.6f\n', PPS_tax_agg);
            end
            Household_Savings = Total_Household_Income - ss.C_agg - General_Tax_Revenue - PAYG_contributions_agg - PPS_contrib_agg - PPS_tax_agg;
            fprintfAndLog(fileID, '      净储蓄 .................................... : %15.6f\n', Household_Savings);
            fprintfAndLog(fileID, '      --------------------------------------------------------------------------------------------------\n');
            fprintfAndLog(fileID, '      >>> 家庭总支出与储蓄 ...................... : %15.6f\n', Total_Household_Income);
            fprintfAndLog(fileID, '   C. 财富与遗赠 (存量与流量):\n');
            fprintfAndLog(fileID, '      持有总资产 (K_p, 下一期初) ................ : %15.6f\n', ss.K_private_hat);
            fprintfAndLog(fileID, '      产生的总遗赠 (留给下一代) ................. : %15.6f\n', ss.Bequest_generated_agg);

            % --- [ IV. 政府部门 ] ---
            fprintfAndLog(fileID, '\n\n--- [ IV. 政府与养老金体系 ] ---\n');
            fprintfAndLog(fileID, '   A. 一般政府预算:\n');
            fprintfAndLog(fileID, '      (+) 一般税收收入 .......................... : %15.6f\n', General_Tax_Revenue);
            if pps_active
                fprintfAndLog(fileID, '      (+) PPS提取税 ............................. : %15.6f\n', PPS_tax_agg);
            end
            fprintfAndLog(fileID, '      (+) 公共资本的隐性回报 .................... : %15.6f\n', Public_Capital_Return_display);
            fprintfAndLog(fileID, '      --------------------------------------------------------------------------------------------------\n');
            fprintfAndLog(fileID, '      (=) 政府总收入 ............................ : %15.6f\n', Total_General_Gov_Revenue);
            fprintfAndLog(fileID, '      (-) 政府消费 (G_c) ........................ : %15.6f\n', G_c_display);
            fprintfAndLog(fileID, '      (-) 政府总投资 (I_g) ...................... : %15.6f\n', I_g_display);
            fprintfAndLog(fileID, '         (公共资本折旧 d_g*K_g) ................ : (%14.6f)\n', Depreciation_g);
            fprintfAndLog(fileID, '         (公共净投资 I_g - d_g*K_g) ............ : (%14.6f)\n', I_g_display - Depreciation_g);
            fprintfAndLog(fileID, '      (-) 对家庭的转移支付 (TR) ................. : %15.6f\n', ss.TR_distributed_agg);
            fprintfAndLog(fileID, '      --------------------------------------------------------------------------------------------------\n');
            fprintfAndLog(fileID, '      (=) 政府总支出 ............................ : %15.6f\n', Total_General_Gov_Expenditure);
            fprintfAndLog(fileID, '      >>> 预算平衡误差 (收入 - 支出) ............ : %15.4e\n', general_budget_balance_error);
            if abs(general_budget_balance_error/Y) < 1e-9
                fprintfAndLog(fileID, '          核验结果: ✅ 政府一般预算按设计严格平衡。\n');
            else
                fprintfAndLog(fileID, '          核验结果: ⚠️ 政府一般预算未按设计平衡！\n');
            end
            fprintfAndLog(fileID, '   B. 现收现付(PAYG)养老金体系:\n');
            fprintfAndLog(fileID, '      (+) 养老金总缴费 .......................... : %15.6f\n', PAYG_contributions_agg);
            fprintfAndLog(fileID, '      (-) 养老金总福利 .......................... : %15.6f\n', PAYG_benefits_agg);
            fprintfAndLog(fileID, '      >>> PAYG体系余额 (缴费 - 福利) ............ : %15.4e\n', payg_balance_error);
            if abs(payg_balance_error/Y) < 1e-9
                fprintfAndLog(fileID, '          核验结果: ✅ PAYG体系按设计严格平衡。\n');
            else
                fprintfAndLog(fileID, '          核验结果: ⚠️ PAYG体系未按设计平衡！\n');
            end
            
            % --- [ V. 核心均衡条件检验 ] ---
            fprintfAndLog(fileID, '\n\n--- [ V. 核心均衡条件检验 (来自求解器的最终误差) ] ---\n');
            final_errors = SS.system_of_equations([ss.r_mkt, ss.w_hat, ss.Bequest_distributed_agg, ss.TR_distributed_agg], Z_ss_norm, ss.L_hat, cS, paramS, struct('A',1.0));
            fprintfAndLog(fileID, '   1. 资本市场出清 (r_guess - r_new) ......... : %15.4e\n', final_errors(1));
            fprintfAndLog(fileID, '   2. 劳动市场出清 (w_guess - w_new) ....... : %15.4e\n', final_errors(2));
            fprintfAndLog(fileID, '   3. 遗赠市场出清 (BQ_guess - BQ_new) ... : %15.4e\n', final_errors(3));
            fprintfAndLog(fileID, '   4. 政府预算平衡 (TR_guess - TR_new) ..... : %15.4e\n', final_errors(4));

            fprintfAndLog(fileID, '\n=========================================================================================================\n');
            fprintfAndLog(fileID, '###                                      报 告 结 束                                      ###\n');
            fprintfAndLog(fileID, '=========================================================================================================\n\n');

            if fileID ~= 1
                fclose(fileID);
                fprintf('\n报告已成功保存至: %s\n', report_filename);
            end

            function fprintfAndLog(fileID, formatSpec, varargin)
                fprintf(1, formatSpec, varargin{:});
                if fileID ~= 1, fprintf(fileID, formatSpec, varargin{:}); end
            end
        end
        function Dist = get_final_dist(ss, cS, paramS, Z_ss_norm)
            [~, ~, Dist] = SS.system_of_equations([ss.r_mkt, ss.w_hat, ss.Bequest_distributed_agg, ss.TR_distributed_agg], Z_ss_norm, ss.L_hat, cS, paramS, struct('A',1.0));
        end

        function polS = get_final_polS(ss, cS, paramS, Z_ss_norm)
            % =========================================================================
            % == 函数: get_final_polS (版本 v1.1 - theta_t 兼容版)
            % ==
            % == 核心修正:
            % ==   - 在构建 M_for_hh 时，添加 theta_t 字段，以匹配 VFI_solver 的新接口。
            % =========================================================================
            M_for_hh.r_mkt_t = ss.r_mkt;
            M_for_hh.w_t = ss.w_hat;
            M_for_hh.tr_per_hh = ss.TR_distributed_agg / sum(Z_ss_norm);
            
            % [核心修正] 计算并添加 b_hat_t 和 theta_t 字段
            total_pension_pot_hat = cS.theta_path(end) * ss.w_hat * ss.L_hat;
            mass_retirees = sum(Z_ss_norm((cS.aR_new+1):end));
            M_for_hh.b_hat_t = total_pension_pot_hat / max(1e-9, mass_retirees);
            M_for_hh.theta_t = cS.theta_path(end);

            [polS, ~] = household.VFI_solver(M_for_hh, paramS, cS);
        end

    end
end