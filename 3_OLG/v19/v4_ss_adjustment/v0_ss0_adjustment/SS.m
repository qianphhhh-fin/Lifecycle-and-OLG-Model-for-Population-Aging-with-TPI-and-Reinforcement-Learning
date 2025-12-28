% --- SS.m ---
classdef SS
    methods (Static)
        function [ss, Dist, polS, valS] = solve_steady_state(cS, paramS, params_ext, verbose, is_bgp_ss)
            % =========================================================================
            % == 函数: solve_steady_state
            % == 版本: [v1.5 - 伪稳态/BGP稳态区分版]
            % ==
            % == 核心修改:
            % ==   - 增加布尔输入参数 'is_bgp_ss'。
            % ==   - 将此标志传递给 system_of_equations，以控制是否强制执行BGP资本积累约束。
            % =========================================================================
            if verbose
                verbose='iter';
            else
                verbose='off';
            end
            cS.A = params_ext.A;
            if isfield(params_ext, 'g_A_ss'), cS.g_A_ss = params_ext.g_A_ss; end
            if isfield(params_ext, 'n_ss'), cS.n_ss = params_ext.n_ss; end
            if isfield(params_ext, 'Z') && ~isempty(params_ext.Z), Z_ss_norm = params_ext.Z;
            else, error('错误：求解稳态需要一个明确的人口分布Z。'); end

            if ~cS.pps_active
                cS.nkpps = 1; cS.n_pps_rate_grid = 1; % 没开启的话一定要设置为1！！！！！
            end

            % 设置初始猜测值
            r_guess_init = 0.04;
            w_guess_init = 1.5;
            K_guess_for_BQ = 2.5;
            BQ_guess_init = mean(1-cS.s_pathV) * K_guess_for_BQ / sum(Z_ss_norm(1:cS.aD_new));
            TR_guess_init = 0.0;
            b_hat_guess_init = 0.05;

            x0 = [r_guess_init, w_guess_init, BQ_guess_init, TR_guess_init, b_hat_guess_init];
            lb = [1e-8, 1e-8, 1e-8, -Inf, 0];
            ub = [Inf, Inf, Inf, Inf, Inf];

            % 将 is_bgp_ss 标志传递给核心方程组
            system_wrapper = @(x) SS.system_of_equations(x, Z_ss_norm, cS, paramS, params_ext, is_bgp_ss);
  
    options = optimoptions('lsqnonlin', 'Display', verbose, ...
    'FunctionTolerance', 1e-14, 'StepTolerance', 1e-14, 'OptimalityTolerance', 1e-14, ...
    'MaxIterations', 500, 'Algorithm', 'levenberg-marquardt'); % 

            [x_eq, ~, ~, exitflag] = lsqnonlin(system_wrapper, x0, lb, ub, options);

            % 使用均衡解生成最终结果
            [~, ss, Dist, polS, valS] = SS.system_of_equations(x_eq, Z_ss_norm, cS, paramS, params_ext, is_bgp_ss);
        end

        % --- SS.m ---
        function [F_error, ss, Dist, polS, valS] = system_of_equations(x_guess, Z_ss_norm, cS, paramS, params_ext, is_bgp_ss)
            % =========================================================================
            % == 函数: system_of_equations
            % == 版本: [v3.2 - 伪稳态/BGP稳态区分版]
            % ==
            % == 核心修改:
            % ==   - 接收 'is_bgp_ss' 标志。
            % ==   - 如果 is_bgp_ss 为 false (即求解伪稳态 ss0),
            % ==     则在计算后将 error_K 强制设为0，从而将其从求解目标中移除。
            % ==     这允许伪稳态存在一个投资缺口，使得其他流量市场可以精确出清。
            % =========================================================================
            r_guess = x_guess(1);
            w_guess = x_guess(2);
            bq_guess = x_guess(3);
            tr_guess = x_guess(4);
            b_hat_guess = x_guess(5);

            % --- 流程 1: 微观聚合 (家庭部门) ---
            [aggr_out, Dist, polS, valS] = SS.run_micro_to_macro(r_guess, w_guess, bq_guess, tr_guess, b_hat_guess, Z_ss_norm, cS, paramS, params_ext);
            L_supply_agg = aggr_out.L_hat;

            % --- 流程 2: [核心逻辑] 计算资本供需和误差 ---
            K_supply_from_dist = aggregates.get_aggregates_from_dist(Dist, cS);
            K_supply_agg_t = K_supply_from_dist.K_p_hat_total;
            K_supply_agg_tplus1_raw = aggr_out.K_p_hat_tplus1_raw;
            n_period = (1 + cS.n_ss)^cS.time_Step - 1;
            error_K = K_supply_agg_t - (K_supply_agg_tplus1_raw / (1 + n_period));

            % --- [!!! 关键修改 !!!] ---
            % 如果不是求解BGP稳态，则忽略资本存量积累误差
            if ~is_bgp_ss
                error_K = 0;
            end

            % --- 流程 3: 宏观计算 (厂商部门) 与要素价格误差 ---
            g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
            g_total_period = (1 + g_A_period) * (1 + n_period) - 1;

            Kg_Y_ratio_factor = cS.I_g_to_Y_ratio_ss / (g_total_period + cS.ddk_g);
            Y_hat_calc = (K_supply_agg_t^cS.alpha * L_supply_agg^(1 - cS.alpha - cS.gamma) * Kg_Y_ratio_factor^cS.gamma)^(1 / (1 - cS.gamma));
            K_g_hat_calc = Kg_Y_ratio_factor * Y_hat_calc;

            prices_new = firm.get_prices_at_t(K_supply_agg_t, K_g_hat_calc, L_supply_agg, cS);
            r_new = prices_new.r_mkt_t;
            w_new = prices_new.w_hat_t;

            error_r = r_guess - r_new;
            error_w = w_guess - w_new;

            % --- 流程 4: 计算其他市场的出清误差 ---
            mass_total = sum(Z_ss_norm);
            bequest_generated_adj = aggr_out.Bequest_generated_agg / (1 + n_period);
            bq_new = bequest_generated_adj / mass_total;
            error_bq = bq_guess - bq_new;

            PAYG_contributions_agg = cS.theta_path(end) * w_new * L_supply_agg;
            General_Tax_Revenue_calc = aggr_out.Total_Tax_Revenue - PAYG_contributions_agg;
            Public_Capital_Return_calc = prices_new.Y_hat_t - (w_new * L_supply_agg) - ((r_new + cS.ddk) * K_supply_agg_t);
            Total_General_Gov_Revenue_calc = General_Tax_Revenue_calc + Public_Capital_Return_calc;
            G_consumption_spending = cS.G_c_to_Y_ratio_ss * prices_new.Y_hat_t;
            G_investment_spending = cS.I_g_to_Y_ratio_ss * prices_new.Y_hat_t;
            tr_new = (Total_General_Gov_Revenue_calc - (G_consumption_spending + G_investment_spending)) / mass_total;
            error_tr = tr_guess - tr_new;

            mass_retirees = sum(Z_ss_norm((cS.aR_new+1):end));
            total_pension_paid_hat = b_hat_guess * mass_retirees;
            total_pension_contrib_hat = cS.theta_path(end) * w_guess * L_supply_agg;
            error_b_hat = total_pension_contrib_hat - total_pension_paid_hat;

            F_error = [error_r; error_w; error_bq; error_tr; error_b_hat; error_K];

            % --- 流程 5: 打包最终的ss结构体 ---
            if nargout > 1
                ss = struct();
                ss.Y_from_production_hat = prices_new.Y_hat_t;
                ss.K_private_hat = K_supply_agg_t;
                ss.K_public_hat = K_g_hat_calc;
                ss.L_hat = L_supply_agg;
                ss.r_mkt = r_new;
                ss.w_hat = w_new;
                ss.C_agg = aggr_out.C_agg;
                ss.I_g = G_investment_spending;
                ss.Total_Tax_Revenue = aggr_out.Total_Tax_Revenue;
                ss.Bequest_generated_agg = bequest_generated_adj;
                ss.Bequest_distributed_agg = bequest_generated_adj;
                ss.Bequest_gen_hat_raw_ss = aggr_out.Bequest_generated_agg;
                ss.TR_distributed_agg = tr_new * mass_total;
                ss.Public_Capital_Return = Public_Capital_Return_calc;
                ss.G_c = G_consumption_spending;
                ss.b_hat = b_hat_guess;
            end
        end
        % --- SS.m ---
        function [macro_state, Dist, polS, valS] = run_micro_to_macro(r_mkt, w_hat, bq, tr, b_hat_guess, Z_ss_norm, cS, paramS, params_ext)
            % =========================================================================
            % == 函数: run_micro_to_macro
            % == 版本: [v1.7 - 接口对齐修复版]
            % ==
            % == 核心修改:
            % ==   - 修正了对 aggregates.get_aggregates 的调用方式。
            % ==   - 现在正确地接收一个结构体 aggrS，然后从该结构体中解包
            % ==     出所需的聚合变量。
            % =========================================================================

            b_hat_for_vfi = b_hat_guess;
            M_for_hh.r_mkt_t = r_mkt;
            M_for_hh.w_t = w_hat;
            M_for_hh.tr_per_hh = tr;
            M_for_hh.b_hat_t = b_hat_for_vfi;
            M_for_hh.theta_t = cS.theta_path(end);

            [polS, valS] = household.VFI_solver(M_for_hh, paramS, cS);
            Dist = distribution.get_ss_dist(polS, paramS, cS, Z_ss_norm, bq);

            % --- [核心修正] 调用新的聚合函数，并接收其返回的结构体 ---
            aggrS = aggregates.get_aggregates(Dist, polS, cS, paramS);

            % --- 将聚合结果从结构体中解包并打包返回 ---
            macro_state = struct();
            macro_state.K_p_hat_tplus1_raw = aggrS.K_p_hat_tplus1_raw;
            macro_state.C_agg = aggrS.C_agg;
            macro_state.Bequest_generated_agg = aggrS.Bequest_generated_agg;
            macro_state.Total_Tax_Revenue = aggrS.Total_Tax_Revenue;
            macro_state.L_hat = aggrS.L_hat;
        end

        function report_struct = display_national_accounts(ss, cS, paramS, Z_ss_norm, report_filename, print_flag, is_bgp_ss)
            % =========================================================================
            % == 函数: display_national_accounts
            % == 版本: [v8.4 - is_bgp_ss 标志传递版]
            % ==
            % == 核心修改:
            % ==   - 增加输入参数 is_bgp_ss。
            % ==   - 在调用 get_final_dist 和 get_final_polS 时传递此标志，
            % ==     以确保在重新计算分布时，使用正确的稳态类型假设。
            % =========================================================================

            if nargin < 6
                print_flag = true; % 默认打印报告
            end
            if nargin < 5
                report_filename = 'ss_national_accounts_report.txt';
            end


            % --- 0. 准备工作：文件句柄和BGP增长因子 ---
            if print_flag
                try
                    fileID = fopen(report_filename, 'w', 'n', 'UTF-8');
                    if fileID == -1, error('无法创建或打开报告文件: %s', report_filename); end
                catch ME
                    warning('创建报告文件时出错，将输出到命令行。错误: %s', ME.message);
                    fileID = 1;
                end
            else
                fileID = -1; % -1 表示不写入文件
            end

            g_A_period = (1+cS.g_A_ss)^cS.time_Step-1;
            n_period = (1+cS.n_ss)^cS.time_Step-1;
            g_total_period = (1+g_A_period)*(1+n_period)-1;
            Y = ss.Y_from_production_hat;
            mass_total = sum(Z_ss_norm);

            % --- 1. 重新计算 PAYG 流量和真实税收 ---
            % [!!! 核心修改 !!!] 将 is_bgp_ss 标志传递下去
            Dist = SS.get_final_dist(ss, cS, paramS, Z_ss_norm, is_bgp_ss);
            polS = SS.get_final_polS(ss, cS, paramS, Z_ss_norm);

            PAYG_contributions_agg = 0;
            PAYG_benefits_agg = 0;
            General_Tax_Revenue = 0;

            pps_active = isfield(cS, 'pps_active') && cS.pps_active;
            K_k_hat_agg = 0; K_pps_hat_agg = 0;
            PPS_contrib_agg = 0; PPS_withdrawal_agg = 0; PPS_tax_agg = 0;

            for ia = 1:cS.aD_new
                mass_dist_ia = Dist(:,:,:,ia);

                labor_income_slice_full = ss.w_hat * cS.ageEffV_new(ia) .* paramS.leGridV(1:cS.nw_expanded)';
                payg_contrib_slice_full = cS.theta_path(end) * labor_income_slice_full;
                payg_contrib_grid = repmat(reshape(payg_contrib_slice_full, [1,1,cS.nw_expanded]), [cS.nk, cS.nkpps, 1]);

                if ia <= cS.aR_new
                    PAYG_contributions_agg = PAYG_contributions_agg + sum(payg_contrib_grid .* mass_dist_ia, 'all');
                else
                    PAYG_benefits_agg = PAYG_benefits_agg + sum(ss.b_hat .* mass_dist_ia, 'all');
                end

                tax_regular_slice = polS(ia).tax_regular;
                general_tax_slice = tax_regular_slice - payg_contrib_grid;
                General_Tax_Revenue = General_Tax_Revenue + sum(general_tax_slice .* mass_dist_ia, 'all');

                if pps_active
                    kpps_prime_slice = polS(ia).kpps_prime;
                    K_pps_hat_agg = K_pps_hat_agg + sum(kpps_prime_slice .* mass_dist_ia, 'all');
                    if ia <= cS.aR_new
                        pps_contrib_slice_full = cS.pps_contrib_rate * labor_income_slice_full;
                        pps_contrib_grid_sim = repmat(reshape(pps_contrib_slice_full, [1,1,cS.nw_expanded]), [cS.nk, cS.nkpps, 1]);
                        PPS_contrib_agg = PPS_contrib_agg + sum(pps_contrib_grid_sim .* mass_dist_ia, 'all');

                    else
                        kpps_current_grid = cS.kppsGridV;
                        kpps_current_full = repmat(reshape(kpps_current_grid, [1, cS.nkpps, 1]), [cS.nk, 1, cS.nw_expanded]);
                        kpps_level = kpps_current_full * (1 + ss.r_mkt);
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
            General_Tax_Revenue = General_Tax_Revenue - PPS_tax_agg;

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
            Net_Capital_Income = ss.r_mkt * ss.K_private_hat;
            Total_Household_Income = Total_Labor_Income + Net_Capital_Income + PAYG_benefits_agg + PPS_withdrawal_agg + ss.Bequest_distributed_agg + ss.TR_distributed_agg;
            Household_Savings = Total_Household_Income - ss.C_agg - General_Tax_Revenue - PAYG_contributions_agg - PPS_contrib_agg - PPS_tax_agg;
            Total_General_Gov_Revenue = General_Tax_Revenue + PPS_tax_agg + Public_Capital_Return_display;
            Total_General_Gov_Expenditure = G_c_display + I_g_display + ss.TR_distributed_agg;
            general_budget_balance_error = Total_General_Gov_Revenue - Total_General_Gov_Expenditure;
            payg_balance_error = PAYG_contributions_agg - PAYG_benefits_agg;

            % --- 3. [核心] 创建并填充 report_struct ---
            report_struct = struct();
            report_struct.Y = Y;
            report_struct.K_p = ss.K_private_hat;
            report_struct.K_g = ss.K_public_hat;
            report_struct.L = ss.L_hat;
            report_struct.r_period = ss.r_mkt;
            report_struct.r_annual = ((1+ss.r_mkt)^(1/cS.time_Step) - 1)*100;
            report_struct.w = ss.w_hat;
            report_struct.C = ss.C_agg;
            report_struct.I_p_acct = I_p_accounting;
            report_struct.I_p_bgp = I_p_demand_bgp;
            report_struct.I_g = I_g_display;
            report_struct.G_c = G_c_display;
            report_struct.Depreciation_p = Depreciation_p;
            report_struct.Depreciation_g = Depreciation_g;
            report_struct.Kp_Y_ratio = ss.K_private_hat / Y;
            report_struct.C_Y_ratio = (ss.C_agg / Y)*100;
            report_struct.Ip_Y_ratio = (I_p_accounting / Y)*100;
            report_struct.Total_Labor_Income = Total_Labor_Income;
            report_struct.Net_Capital_Income = Net_Capital_Income;
            report_struct.Public_Capital_Return = Public_Capital_Return_display;
            report_struct.Household_Savings = Household_Savings;
            report_struct.RC_Error = resource_constraint_error;
            report_struct.Invest_Gap = I_p_accounting - I_p_demand_bgp;
            report_struct.RC_Error_pct = (resource_constraint_error / Y)*100;
            report_struct.Invest_Gap_pct = (report_struct.Invest_Gap / Y)*100;
            report_struct.TR = ss.TR_distributed_agg;
            report_struct.Tax_General = General_Tax_Revenue;
            report_struct.Gov_Budget_Balance_Error = general_budget_balance_error;
            report_struct.PAYG_contrib = PAYG_contributions_agg;
            report_struct.PAYG_benefit = PAYG_benefits_agg;
            report_struct.PAYG_Balance_Error = payg_balance_error;

            if pps_active
                report_struct.K_k = K_k_hat_agg;
                report_struct.K_pps = K_pps_hat_agg;
                report_struct.PPS_contrib = PPS_contrib_agg;
                report_struct.PPS_withdrawal = PPS_withdrawal_agg;
                report_struct.PPS_tax = PPS_tax_agg;
            end

            % --- 报告开始 ---
            if print_flag
                fprintfAndLog(fileID, '\n\n=========================================================================================================\n');
                if pps_active
                    fprintfAndLog(fileID, '###                 稳 态 国 民 经 济 核 算 报 告 (v8.3 - PPS变量输出增强版)                     ###\n');
                else
                    fprintfAndLog(fileID, '###                       稳 态 国 民 经 济 核 算 报 告 (v8.2 - 最终完整版)                     ###\n');
                end
                fprintfAndLog(fileID, '=========================================================================================================\n');

                fprintfAndLog(fileID, '\n--- [ I. 宏观总量与国民支出 ] ---\n');
                fprintfAndLog(fileID, '   A. 生产 (GDP) 与核心比率:\n');
                fprintfAndLog(fileID, '      国内生产总值 (Y) .......................... : %15.6f\n', report_struct.Y);
                fprintfAndLog(fileID, '      私人资本/产出比 (K_p/Y) ................... : %15.6f\n', report_struct.Kp_Y_ratio);
                fprintfAndLog(fileID, '      公共资本/产出比 (K_g/Y) ................... : %15.6f\n', report_struct.K_g / report_struct.Y);
                fprintfAndLog(fileID, '   B. 均衡价格:\n');
                fprintfAndLog(fileID, '      真实利率 (r, 模型期) ...................... : %15.6f (年化: %.4f %%)\n', report_struct.r_period, report_struct.r_annual);
                fprintfAndLog(fileID, '      真实工资率 (w, 单位有效劳动) .............. : %15.6f\n', report_struct.w);

                fprintfAndLog(fileID, '   C. 国民支出 (Y = C + I_p + I_g + G_c) 与资源约束检验:\n');
                fprintfAndLog(fileID, '      (+) 私人消费 (C) .......................... : %15.6f  (%6.2f %% of Y)\n', report_struct.C, report_struct.C_Y_ratio);
                fprintfAndLog(fileID, '      (+) 私人总投资 (I_p, 会计值) ............ : %15.6f  (%6.2f %% of Y)\n', report_struct.I_p_acct, report_struct.Ip_Y_ratio);
                fprintfAndLog(fileID, '      (+) 公共总投资 (I_g) ...................... : %15.6f  (%6.2f %% of Y)\n', report_struct.I_g, (report_struct.I_g / Y)*100);
                fprintfAndLog(fileID, '      (+) 政府消费 (G_c) ........................ : %15.6f  (%6.2f %% of Y)\n', report_struct.G_c, (report_struct.G_c / Y)*100);
                fprintfAndLog(fileID, '      --------------------------------------------------------------------------------------------------\n');
                fprintfAndLog(fileID, '      (=) 总支出 (C + I_p_acct + I_g + G_c) ..... : %15.6f\n', Total_Expenditure);
                fprintfAndLog(fileID, '      --------------------------------------------------------------------------------------------------\n');
                fprintfAndLog(fileID, '      >>> 资源约束误差 (Y - 总支出) ............. : %15.4e (Y的 %.4e %%)\n', report_struct.RC_Error, report_struct.RC_Error_pct);
                if abs(report_struct.RC_Error_pct) < 1e-6, fprintfAndLog(fileID, '          核验结果: ✅ 资源约束在会计上严格满足。\n');
                else, fprintfAndLog(fileID, '          核验结果: ⚠️ 资源约束在会计上存在误差！\n'); end

                fprintfAndLog(fileID, '\n   D. 投资一致性检验 (仅适用于真稳态):\n');
                fprintfAndLog(fileID, '      BGP理论投资需求 (g+d)K .................. : %15.6f\n', report_struct.I_p_bgp);
                fprintfAndLog(fileID, '      会计投资 (Y-C-I_g-G_c) .................... : %15.6f\n', report_struct.I_p_acct);
                fprintfAndLog(fileID, '      >>> 投资缺口 (会计值 - 理论值) .......... : %15.4e (Y的 %.4f %%)\n', report_struct.Invest_Gap, report_struct.Invest_Gap_pct);
                if abs(report_struct.Invest_Gap_pct) < 1e-4, fprintfAndLog(fileID, '          检验结果: ✅ 经济体处于或接近真实的BGP稳态。\n');
                else, fprintfAndLog(fileID, '          检验结果: ⚠️  经济体偏离BGP稳态 (预期中的伪稳态特征)。\n'); end

                fprintfAndLog(fileID, '\n\n--- [ II. 厂商部门 ] ---\n');
                fprintfAndLog(fileID, '   A. 生产要素投入 (存量):\n');
                fprintfAndLog(fileID, '      总私人资本存量 (K_p) ...................... : %15.6f\n', report_struct.K_p);
                if pps_active, fprintfAndLog(fileID, '         其中: 常规资本 (K_k) ................. : %15.6f  (%6.2f %% of K_p)\n', K_k_hat_agg, (K_k_hat_agg / report_struct.K_p) * 100); fprintfAndLog(fileID, '         其中: PPS资本 (K_pps) ................ : %15.6f  (%6.2f %% of K_p)\n', K_pps_hat_agg, (K_pps_hat_agg / report_struct.K_p) * 100); end
                fprintfAndLog(fileID, '      公共资本存量 (K_g) ........................ : %15.6f\n', report_struct.K_g);
                fprintfAndLog(fileID, '      有效劳动需求 (L) .......................... : %15.6f\n', report_struct.L);
                fprintfAndLog(fileID, '   B. 要素报酬 (流量):\n');
                fprintfAndLog(fileID, '      支付工资总额 (w*L) ........................ : %15.6f  (Y的 %6.2f %%)\n', report_struct.Total_Labor_Income, (report_struct.Total_Labor_Income/Y)*100);
                fprintfAndLog(fileID, '      私人资本总回报 ((r+d)*K_p) ................ : %15.6f  (Y的 %6.2f %%)\n', Total_Capital_Income_Gross, (Total_Capital_Income_Gross/Y)*100);
                fprintfAndLog(fileID, '      公共资本回报 (Y - wL - (r+d)K_p) ........ : %15.6f  (Y的 %6.2f %%)\n', report_struct.Public_Capital_Return, (report_struct.Public_Capital_Return/Y)*100);
                fprintfAndLog(fileID, '   C. 投资与折旧 (流量):\n');
                fprintfAndLog(fileID, '      私人总投资 (I_p, BGP理论值) ............... : %15.6f\n', report_struct.I_p_bgp);
                fprintfAndLog(fileID, '      私人资本折旧 (d*K_p) ...................... : %15.6f\n', report_struct.Depreciation_p);
                fprintfAndLog(fileID, '      私人净投资 (I_p_bgp - d*K_p) .............. : %15.6f\n', report_struct.I_p_bgp - report_struct.Depreciation_p);

                fprintfAndLog(fileID, '\n\n--- [ III. 家庭部门 ] ---\n');
                fprintfAndLog(fileID, '   A. 收入来源 (流量):\n');
                fprintfAndLog(fileID, '      劳动收入 (来自厂商) ....................... : %15.6f\n', report_struct.Total_Labor_Income);
                fprintfAndLog(fileID, '      净资本收入 (r*K_p) ........................ : %15.6f\n', report_struct.Net_Capital_Income);
                fprintfAndLog(fileID, '      养老金福利 ................................ : %15.6f\n', report_struct.PAYG_benefit);
                if pps_active, fprintfAndLog(fileID, '      PPS账户提取 ............................... : %15.6f\n', PPS_withdrawal_agg); end
                fprintfAndLog(fileID, '      收到的遗赠 ................................ : %15.6f\n', ss.Bequest_distributed_agg);
                fprintfAndLog(fileID, '      收到的转移支付 ............................ : %15.6f\n', report_struct.TR);
                fprintfAndLog(fileID, '      --------------------------------------------------------------------------------------------------\n');
                fprintfAndLog(fileID, '      >>> 家庭总收入 ............................ : %15.6f\n', Total_Household_Income);
                fprintfAndLog(fileID, '   B. 收入用途 (流量):\n');
                fprintfAndLog(fileID, '      消费支出 (C) .............................. : %15.6f\n', report_struct.C);
                fprintfAndLog(fileID, '      支付的税收 (资本税/劳动税/消费税) ......... : %15.6f\n', report_struct.Tax_General);
                fprintfAndLog(fileID, '      养老金缴费 ................................ : %15.6f\n', report_struct.PAYG_contrib);
                if pps_active, fprintfAndLog(fileID, '      PPS账户缴费 (储蓄) ...................... : %15.6f\n', PPS_contrib_agg); fprintfAndLog(fileID, '      PPS提取税 ................................. : %15.6f\n', PPS_tax_agg); end
                fprintfAndLog(fileID, '      净储蓄 .................................... : %15.6f\n', report_struct.Household_Savings);
                fprintfAndLog(fileID, '      --------------------------------------------------------------------------------------------------\n');
                fprintfAndLog(fileID, '      >>> 家庭总支出与储蓄 ...................... : %15.6f\n', Total_Household_Income);
                fprintfAndLog(fileID, '   C. 财富与遗赠 (存量与流量):\n');
                fprintfAndLog(fileID, '      持有总资产 (K_p, 下一期初) ................ : %15.6f\n', report_struct.K_p);
                fprintfAndLog(fileID, '      产生的总遗赠 (留给下一代) ................. : %15.6f\n', ss.Bequest_generated_agg);

                fprintfAndLog(fileID, '\n\n--- [ IV. 政府与养老金体系 ] ---\n');
                fprintfAndLog(fileID, '   A. 一般政府预算:\n');
                fprintfAndLog(fileID, '      (+) 一般税收收入 .......................... : %15.6f\n', report_struct.Tax_General);
                if pps_active, fprintfAndLog(fileID, '      (+) PPS提取税 ............................. : %15.6f\n', PPS_tax_agg); end
                fprintfAndLog(fileID, '      (+) 公共资本的隐性回报 .................... : %15.6f\n', report_struct.Public_Capital_Return);
                fprintfAndLog(fileID, '      --------------------------------------------------------------------------------------------------\n');
                fprintfAndLog(fileID, '      (=) 政府总收入 ............................ : %15.6f\n', Total_General_Gov_Revenue);
                fprintfAndLog(fileID, '      (-) 政府消费 (G_c) ........................ : %15.6f\n', report_struct.G_c);
                fprintfAndLog(fileID, '      (-) 政府总投资 (I_g) ...................... : %15.6f\n', report_struct.I_g);
                fprintfAndLog(fileID, '         (公共资本折旧 d_g*K_g) ................ : (%14.6f)\n', Depreciation_g);
                fprintfAndLog(fileID, '         (公共净投资 I_g - d_g*K_g) ............ : (%14.6f)\n', report_struct.I_g - Depreciation_g);
                fprintfAndLog(fileID, '      (-) 对家庭的转移支付 (TR) ................. : %15.6f\n', report_struct.TR);
                fprintfAndLog(fileID, '      --------------------------------------------------------------------------------------------------\n');
                fprintfAndLog(fileID, '      (=) 政府总支出 ............................ : %15.6f\n', Total_General_Gov_Expenditure);
                fprintfAndLog(fileID, '      >>> 预算平衡误差 (收入 - 支出) ............ : %15.4e\n', report_struct.Gov_Budget_Balance_Error);
                if abs(report_struct.Gov_Budget_Balance_Error/Y) < 1e-9, fprintfAndLog(fileID, '          核验结果: ✅ 政府一般预算按设计严格平衡。\n');
                else, fprintfAndLog(fileID, '          核验结果: ⚠️ 政府一般预算未按设计平衡！\n'); end
                fprintfAndLog(fileID, '   B. 现收现付(PAYG)养老金体系:\n');
                fprintfAndLog(fileID, '      (+) 养老金总缴费 .......................... : %15.6f\n', report_struct.PAYG_contrib);
                fprintfAndLog(fileID, '      (-) 养老金总福利 .......................... : %15.6f\n', report_struct.PAYG_benefit);
                fprintfAndLog(fileID, '      >>> PAYG体系余额 (缴费 - 福利) ............ : %15.4e\n', report_struct.PAYG_Balance_Error);
                if abs(report_struct.PAYG_Balance_Error/Y) < 1e-9, fprintfAndLog(fileID, '          核验结果: ✅ PAYG体系按设计严格平衡。\n');
                else, fprintfAndLog(fileID, '          核验结果: ⚠️ PAYG体系未按设计平衡！\n'); end

                fprintfAndLog(fileID, '\n\n--- [ V. 核心均衡条件检验 (来自求解器的最终误差) ] ---\n');
                bq_per_capita_final = ss.Bequest_distributed_agg / mass_total;
                tr_per_capita_final = ss.TR_distributed_agg / mass_total;
                x_final_pc = [ss.r_mkt, ss.w_hat, bq_per_capita_final, tr_per_capita_final, ss.b_hat];
                final_errors = SS.system_of_equations(x_final_pc, Z_ss_norm, cS, paramS, struct('A',1.0), is_bgp_ss);
                fprintfAndLog(fileID, '   1. 资本市场出清 (r_guess - r_new) ......... : %15.4e\n', final_errors(1));
                fprintfAndLog(fileID, '   2. 劳动市场出清 (w_guess - w_new) ....... : %15.4e\n', final_errors(2));
                fprintfAndLog(fileID, '   3. 遗赠市场出清 (BQ_guess - BQ_new) ... : %15.4e\n', final_errors(3));
                fprintfAndLog(fileID, '   4. 政府预算平衡 (TR_guess - TR_new) ..... : %15.4e\n', final_errors(4));
                fprintfAndLog(fileID, '   5. PAYG体系平衡 (Contrib - Benefit) .... : %15.4e\n', final_errors(5));
                if (length(final_errors) > 5)
                    fprintfAndLog(fileID, '   6. BGP资本积累 (K_t - K_t+1_adj) ........ : %15.4e\n', final_errors(6));
                end

                fprintfAndLog(fileID, '\n=========================================================================================================\n');
                fprintfAndLog(fileID, '###                                      报 告 结 束                                      ###\n');
                fprintfAndLog(fileID, '=========================================================================================================\n\n');

                if fileID > 1
                    fclose(fileID);
                    fprintf('\n报告已成功保存至: %s\n', report_filename);
                end
            end % end of print_flag condition

            function fprintfAndLog(fileID, formatSpec, varargin)
                if print_flag
                    fprintf(1, formatSpec, varargin{:});
                    if fileID > 1, fprintf(fileID, formatSpec, varargin{:}); end
                end
            end
        end

        function Dist = get_final_dist(ss, cS, paramS, Z_ss_norm, is_bgp_ss)
            % =========================================================================
            % == 函数: get_final_dist
            % == 版本: [v1.4 - is_bgp_ss 标志传递版]
            % ==
            % == 核心修改:
            % ==   - 增加输入参数 is_bgp_ss。
            % ==   - 在调用 system_of_equations 时，将此标志传递下去。
            % =========================================================================
            mass_total = sum(Z_ss_norm);
            bq_per_capita = ss.Bequest_distributed_agg / mass_total;
            tr_per_capita = ss.TR_distributed_agg / mass_total;

            x_eq_final = [ss.r_mkt, ss.w_hat, bq_per_capita, tr_per_capita, ss.b_hat];
            % [!!! 核心修改 !!!] 将 is_bgp_ss 标志传递给核心方程组
            [~, ~, Dist] = SS.system_of_equations(x_eq_final, Z_ss_norm, cS, paramS, struct('A',1.0), is_bgp_ss);
        end



        function polS = get_final_polS(ss, cS, paramS, Z_ss_norm)
            % =========================================================================
            % == 函数: get_final_polS
            % == 版本: [v1.3 - 养老金福利(b_hat)求解版]
            % ==
            % == 核心修改:
            % ==   - 在构建 M_for_hh 时，直接使用 ss 中存储的均衡养老金福利 ss.b_hat。
            % =========================================================================
            M_for_hh.r_mkt_t = ss.r_mkt;
            M_for_hh.w_t = ss.w_hat;
            M_for_hh.tr_per_hh = ss.TR_distributed_agg / sum(Z_ss_norm);

            % [核心修正] 直接使用最终的、均衡的 b_hat
            M_for_hh.b_hat_t = ss.b_hat;
            M_for_hh.theta_t = cS.theta_path(end);

            [polS, ~] = household.VFI_solver(M_for_hh, paramS, cS);
        end

    end
end