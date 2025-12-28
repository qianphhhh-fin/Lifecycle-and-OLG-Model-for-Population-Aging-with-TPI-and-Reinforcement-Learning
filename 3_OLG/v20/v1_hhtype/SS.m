% --- SS.m ---
classdef SS
    methods (Static)


       function [ss, Dist, polS, valS] = solve_steady_state_ss0(cS, paramS, params_ext)
            % =========================================================================
            % == 函数: solve_steady_state_ss0 (新增)
            % == 版本: [v1.0 - 初始时期专用求解器]
            % ==
            % == 目的:
            % ==   - 求解转轨路径的初始时期(t=1)，此时资本市场【不】要求出清。
            % ==   - 利率和工资由厂商的需求侧（即给定的K和L）决定。
            % ==   - 求解器寻找一组(r, w, TR, b, BQ)，使得除了资本市场外的
            % ==     所有其他市场（劳动、商品、政府预算、养老金、遗赠）出清。
            % =========================================================================
            verbose = 'iter'; % 在调试初始时期时，建议始终显示迭代过程
            cS.A = params_ext.A;
            if isfield(params_ext, 'g_A_ss'), cS.g_A_ss = params_ext.g_A_ss; end
            if isfield(params_ext, 'n_ss'), cS.n_ss = params_ext.n_ss; end
            if isfield(params_ext, 'Z') && ~isempty(params_ext.Z), Z_ss_norm = params_ext.Z;
            else, error('错误：求解初始稳态需要一个明确的人口分布Z。'); end
            
            % --- 初始猜测值 ---
            r_guess_init = 0.04;
            w_guess_init = 1.5;
            tr_guess_init = 0.0;
            b_hat_by_type_guess = [0.1; 0.08; 0.03; 0.02];
            
            % [核心差异] K/Y比率是外生目标，不再内生求解资本存量本身
            % 但在求解器内部，我们仍需迭代遗赠和分布，所以需要一个BQ的初始猜测
            K_guess_for_BQ_total = 3.5; % 使用一个合理的资本水平来猜测遗赠规模
            BQ_guess_init_by_type = mean(1-cS.s_pathV) * K_guess_for_BQ_total * cS.type_weights;

            % 构造求解向量
            x0 = [r_guess_init; w_guess_init; tr_guess_init; b_hat_by_type_guess; BQ_guess_init_by_type];
            lb = [0; 1e-8; -Inf; repmat(0, cS.nTypes, 1); repmat(1e-8, cS.nTypes, 1)];
            ub = [Inf; Inf; Inf; repmat(Inf, cS.nTypes, 1); repmat(Inf, cS.nTypes, 1)];

            % 创建匿名函数，将额外参数传入误差函数
            system_wrapper = @(x) SS.system_of_equations_ss0(x, Z_ss_norm, cS, paramS, params_ext);

            options = optimoptions('lsqnonlin', 'Display', verbose, ...
                'FunctionTolerance', 1e-14, 'StepTolerance', 1e-14, 'OptimalityTolerance', 1e-14, ...
                'MaxIterations', 500, 'Algorithm', 'trust-region-reflective');

            [x_eq, ~, ~, exitflag] = lsqnonlin(system_wrapper, x0, lb, ub, options);

            % [事后处理] 使用求解出的均衡价格，计算最终的宏观量和分布
            % is_final_run 标志位告诉内层函数，这是最后一次运行，需要返回所有结果
            is_final_run = true;
            [~, ss, Dist, polS, valS] = SS.system_of_equations_ss0(x_eq, Z_ss_norm, cS, paramS, params_ext, is_final_run);
        end

        function [F_error, ss, Dist, polS, valS] = system_of_equations_ss0(x_guess, Z_ss_norm, cS, paramS, params_ext, is_final_run)
            % =========================================================================
            % == 函数: system_of_equations_ss0 (新增)
            % == 版本: [v1.0]
            % ==
            % == 目的:
            % ==   - 为 solve_steady_state_ss0 提供误差向量。
            % ==   - 核心逻辑:
            % ==     1. 给定猜测的价格，通过VFI和分布模拟，得到一个【内生的】资本和劳动供给。
            % ==     2. 用这个内生的K和L，通过【厂商需求侧】反算出理论价格。
            % ==     3. 计算猜测价格与理论价格的误差。
            % ==     4. 同时计算其他所有市场的预算平衡误差。
            % ==     5. 【不】计算资本市场BGP积累误差。
            % =========================================================================
            
            if nargin < 6, is_final_run = false; end
            ss = []; Dist = []; polS = []; valS = [];

            % --- 1. 解码猜测值 ---
            r_guess = x_guess(1);
            w_guess = x_guess(2);
            tr_guess = x_guess(3);
            b_hat_by_type_guess = x_guess(4 : 3+cS.nTypes);
            bq_by_type_guess = x_guess(4+cS.nTypes : 3+2*cS.nTypes);

            % --- 2. 微观聚合，得到一个内生的K和L ---
            [aggr_out, Dist, polS, valS] = SS.run_micro_to_macro(r_guess, w_guess, bq_by_type_guess, tr_guess, b_hat_by_type_guess, Z_ss_norm, cS, paramS, params_ext);
            L_supply_agg_total = aggr_out.L_hat;
            K_supply_from_dist = aggregates.get_aggregates_from_dist(Dist, cS);
            K_supply_agg_t = K_supply_from_dist.K_p_hat_total;

            % --- 3. 价格由需求侧决定 ---
            prices_demand_side = firm.get_prices_from_KL(K_supply_agg_t, L_supply_agg_total, cS, params_ext);
            r_demand = prices_demand_side.r_mkt_t;
            w_demand = prices_demand_side.w_hat_t;

            % --- 4. 计算误差向量 ---
            % a. 要素价格误差
            error_r = r_guess - r_demand;
            error_w = w_guess - w_demand;
            
            % b. 政府、养老金、遗赠市场误差 (这部分逻辑与标准稳态求解器完全相同)
            Y_hat_demand = prices_demand_side.Y_hat_t;
            n_period = (1 + cS.n_ss)^cS.time_Step - 1;

            % b.1 遗赠市场
            bequest_generated_by_type_adj = aggr_out.Bequest_generated_by_type_agg / (1 + n_period);
            error_bq_by_type = bq_by_type_guess - bequest_generated_by_type_adj;

            % b.2 PAYG预算
            PAYG_contrib_by_type = zeros(cS.nTypes, 1);
            mass_retirees_by_type = zeros(cS.nTypes, 1);
            L_hat_by_type = aggr_out.L_hat_by_type;
            labor_income_by_type = w_demand * L_hat_by_type; % 使用需求侧工资
            for i_type = 1:cS.nTypes
                Dist_type = Dist(:,:,:,:,i_type);
                if i_type <= 2, theta_ss_type = cS.theta_path_urban(end);
                else, theta_ss_type = cS.theta_path_resident(end); end
                PAYG_contrib_by_type(i_type) = theta_ss_type * labor_income_by_type(i_type);
                mass_retirees_this_type = 0;
                for ia = (cS.aR_new+1):cS.aD_new, mass_retirees_this_type = mass_retirees_this_type + sum(Dist_type(:,:,:,ia), 'all'); end
                mass_retirees_by_type(i_type) = mass_retirees_this_type;
            end
            total_pension_paid_by_type = b_hat_by_type_guess .* mass_retirees_by_type;
            error_b_by_type = PAYG_contrib_by_type - total_pension_paid_by_type;

            % b.3 政府一般预算
            PAYG_contributions_total_agg = sum(PAYG_contrib_by_type);
            General_Tax_Revenue_calc = aggr_out.Total_Tax_Revenue - PAYG_contributions_total_agg;
            Public_Capital_Return_calc = Y_hat_demand - (w_demand * L_supply_agg_total) - ((r_demand + cS.ddk) * K_supply_agg_t);
            Total_General_Gov_Revenue_calc = General_Tax_Revenue_calc + Public_Capital_Return_calc;
            G_consumption_spending = cS.G_c_to_Y_ratio_ss * Y_hat_demand;
            G_investment_spending = cS.I_g_to_Y_ratio_ss * Y_hat_demand;
            mass_total = sum(Z_ss_norm);
            tr_new = (Total_General_Gov_Revenue_calc - (G_consumption_spending + G_investment_spending)) / mass_total;
            error_tr = tr_guess - tr_new;

            % 最终误差向量 (注意：没有资本积累误差)
            F_error = [error_r; error_w; error_tr; error_b_by_type; error_bq_by_type];

            % --- 5. 如果是最后一次运行，打包ss结构体 ---
            if is_final_run
                ss = struct();
                ss.Y_from_production_hat = Y_hat_demand;
                ss.K_private_hat = K_supply_agg_t;
                ss.K_public_hat = prices_demand_side.Y_hat_t * (cS.I_g_to_Y_ratio_ss / ((1 + (1+params_ext.g_A_ss)^cS.time_Step-1)*(1+(1+params_ext.n_ss)^cS.time_Step-1)-1 + cS.ddk_g));
                ss.L_hat = L_supply_agg_total;
                ss.r_mkt = r_demand;
                ss.w_hat = w_demand;
                ss.C_agg = aggr_out.C_agg;
                ss.I_g = G_investment_spending;
                ss.Total_Tax_Revenue = aggr_out.Total_Tax_Revenue;
                ss.Bequest_generated_by_type_agg = bequest_generated_by_type_adj;
                ss.Bequest_generated_agg = sum(bequest_generated_by_type_adj);
                ss.Bequest_distributed_agg = ss.Bequest_generated_agg;
                ss.Bequest_gen_hat_raw_ss = aggr_out.Bequest_generated_agg;
                ss.TR_distributed_agg = tr_new * mass_total;
                ss.Public_Capital_Return = Public_Capital_Return_calc;
                ss.G_c = G_consumption_spending;
                ss.b_hat_by_type = b_hat_by_type_guess;
                
                % [重要] 记录此时的投资缺口，它反映了t=1的非均衡性质
                K_supply_agg_tplus1_raw = aggr_out.K_p_hat_tplus1_raw;
                ss.K_p_hat_tplus1_raw = K_supply_agg_tplus1_raw;
                ss.investment_gap_ss0 = K_supply_agg_tplus1_raw / (1 + n_period) - K_supply_agg_t;
            end
        end

        function [ss, Dist, polS, valS] = solve_steady_state(cS, paramS, params_ext, verbose, is_bgp_ss)
            % =========================================================================
            % == 函数: solve_steady_state
            % == 版本: [v1.9]
            % ==
            % == 核心修改:
            % ==   - 将 ss0 的特殊求解逻辑移至新的专用函数 solve_steady_state_ss0。
            % ==   - 此函数现在只负责求解真正的BGP稳态。
            % ==   - 增加了对 is_bgp_ss 输入的检查，确保不会被误用。
            % =========================================================================
            if ~is_bgp_ss
                error('solve_steady_state 函数现在只用于求解BGP稳态 (is_bgp_ss=true)。请使用 solve_steady_state_ss0 求解初始时期。');
            end

            if verbose, verbose='iter'; else, verbose='off'; end
            cS.A = params_ext.A;
            if isfield(params_ext, 'g_A_ss'), cS.g_A_ss = params_ext.g_A_ss; end
            if isfield(params_ext, 'n_ss'), cS.n_ss = params_ext.n_ss; end
            if isfield(params_ext, 'Z') && ~isempty(params_ext.Z), Z_ss_norm = params_ext.Z;
            else, error('错误：求解稳态需要一个明确的人口分布Z。'); end

            if ~cS.pps_active, cS.nkpps = 1; cS.n_pps_rate_grid = 1; end

            r_guess_init = 0.04;
            w_guess_init = 1.5;
            tr_guess_init = 0.0;
            b_hat_by_type_guess = 0.1 * ones(cS.type_weights,1);
            K_guess_for_BQ_total = 3.5;
            BQ_guess_init_by_type = mean(1-cS.s_pathV) * K_guess_for_BQ_total * cS.type_weights;

            x0 = [r_guess_init; w_guess_init; tr_guess_init; b_hat_by_type_guess; BQ_guess_init_by_type];
            lb = [0; 1e-8; -Inf; repmat(0, cS.nTypes, 1); repmat(1e-8, cS.nTypes, 1)];
            ub = [Inf; Inf; Inf; repmat(Inf, cS.nTypes, 1); repmat(Inf, cS.nTypes, 1)];

            % 注意，这里调用的是标准的 system_of_equations
            system_wrapper = @(x) SS.system_of_equations(x, Z_ss_norm, cS, paramS, params_ext, true);

            options = optimoptions('lsqnonlin', 'Display', verbose, ...
                'FunctionTolerance', 1e-12, 'StepTolerance', 1e-12, 'OptimalityTolerance', 1e-12, ...
                'MaxIterations', 500); % , 'Algorithm', 'trust-region-reflective'

            [x_eq, ~, ~, exitflag] = lsqnonlin(system_wrapper, x0, lb, ub, options);

            is_final_run = true;
            [~, ss, Dist, polS, valS] = SS.system_of_equations(x_eq, Z_ss_norm, cS, paramS, params_ext, true, is_final_run);
        end


        % ... system_of_equations 函数保持原样，但我会把is_final_run的逻辑加进去 ...
        function [F_error, ss, Dist, polS, valS] = system_of_equations(x_guess, Z_ss_norm, cS, paramS, params_ext, is_bgp_ss, is_final_run)
            if nargin < 7, is_final_run = false; end
            if nargin < 6, is_bgp_ss = true; end
            ss = []; Dist = []; polS = []; valS = [];

            r_guess = x_guess(1);
            w_guess = x_guess(2);
            tr_guess = x_guess(3);
            b_hat_by_type_guess = x_guess(4 : 3+cS.nTypes);
            bq_by_type_guess = x_guess(4+cS.nTypes : 3+2*cS.nTypes);

            [aggr_out, Dist, polS, valS] = SS.run_micro_to_macro(r_guess, w_guess, bq_by_type_guess, tr_guess, b_hat_by_type_guess, Z_ss_norm, cS, paramS, params_ext);
            L_supply_agg_total = aggr_out.L_hat;
            
            K_supply_from_dist = aggregates.get_aggregates_from_dist(Dist, cS);
            K_supply_agg_t = K_supply_from_dist.K_p_hat_total;
            K_supply_agg_tplus1_raw = aggr_out.K_p_hat_tplus1_raw;
            n_period = (1 + cS.n_ss)^cS.time_Step - 1;
            error_K = K_supply_agg_t - (K_supply_agg_tplus1_raw / (1 + n_period));
            
            if ~is_bgp_ss, error_K = 0; end

            prices_new = firm.get_prices_from_KL(K_supply_agg_t, L_supply_agg_total, cS, params_ext);
            r_new = prices_new.r_mkt_t;
            w_new = prices_new.w_hat_t;
            error_r = r_guess - r_new;
            error_w = w_guess - w_new;

            PAYG_contrib_by_type = zeros(cS.nTypes, 1);
            mass_retirees_by_type = zeros(cS.nTypes, 1);
            L_hat_by_type = aggr_out.L_hat_by_type;
            labor_income_by_type = w_new * L_hat_by_type;
            for i_type = 1:cS.nTypes
                Dist_type = Dist(:,:,:,:,i_type);
                if i_type <= 2, theta_ss_type = cS.theta_path_urban(end);
                else, theta_ss_type = cS.theta_path_resident(end); end
                PAYG_contrib_by_type(i_type) = theta_ss_type * labor_income_by_type(i_type);
                mass_retirees_this_type = 0;
                for ia = (cS.aR_new+1):cS.aD_new, mass_retirees_this_type = mass_retirees_this_type + sum(Dist_type(:,:,:,ia), 'all'); end
                mass_retirees_by_type(i_type) = mass_retirees_this_type;
            end
            bequest_generated_by_type_adj = aggr_out.Bequest_generated_by_type_agg / (1 + n_period);
            error_bq_by_type = bq_by_type_guess - bequest_generated_by_type_adj;
            total_pension_paid_by_type = b_hat_by_type_guess .* mass_retirees_by_type;
            error_b_by_type = PAYG_contrib_by_type - total_pension_paid_by_type;
            PAYG_contributions_total_agg = sum(PAYG_contrib_by_type);
            General_Tax_Revenue_calc = aggr_out.Total_Tax_Revenue - PAYG_contributions_total_agg;
            Public_Capital_Return_calc = prices_new.Y_hat_t - (w_new * L_supply_agg_total) - ((r_new + cS.ddk) * K_supply_agg_t);
            Total_General_Gov_Revenue_calc = General_Tax_Revenue_calc + Public_Capital_Return_calc;
            G_consumption_spending = cS.G_c_to_Y_ratio_ss * prices_new.Y_hat_t;
            G_investment_spending = cS.I_g_to_Y_ratio_ss * prices_new.Y_hat_t;
            mass_total = sum(Z_ss_norm);
            tr_new = (Total_General_Gov_Revenue_calc - (G_consumption_spending + G_investment_spending)) / mass_total;
            error_tr = tr_guess - tr_new;

            F_error = [error_r; error_w; error_tr; error_b_by_type; error_bq_by_type; error_K];

            if is_final_run
                ss = struct();
                ss.Y_from_production_hat = prices_new.Y_hat_t;
                ss.K_private_hat = K_supply_agg_t;
                ss.K_public_hat = prices_new.Y_hat_t * (cS.I_g_to_Y_ratio_ss / ((1 + (1+params_ext.g_A_ss)^cS.time_Step-1)*(1+(1+params_ext.n_ss)^cS.time_Step-1)-1 + cS.ddk_g));
                ss.L_hat = L_supply_agg_total;
                ss.r_mkt = r_new;
                ss.w_hat = w_new;
                ss.C_agg = aggr_out.C_agg;
                ss.I_g = G_investment_spending;
                ss.Total_Tax_Revenue = aggr_out.Total_Tax_Revenue;
                ss.Bequest_generated_by_type_agg = bequest_generated_by_type_adj;
                ss.Bequest_generated_agg = sum(bequest_generated_by_type_adj);
                ss.Bequest_distributed_agg = ss.Bequest_generated_agg;
                ss.Bequest_gen_hat_raw_ss = aggr_out.Bequest_generated_agg;
                ss.TR_distributed_agg = tr_new * mass_total;
                ss.Public_Capital_Return = Public_Capital_Return_calc;
                ss.G_c = G_consumption_spending;
                ss.b_hat_by_type = b_hat_by_type_guess;
            end
        end



        function [macro_state, Dist, polS, valS] = run_micro_to_macro(r_mkt, w_hat, bq_by_type, tr, b_hat_by_type_guess, Z_ss_norm, cS, paramS, params_ext)
            % =========================================================================
            % == 函数: run_micro_to_macro
            % == 版本: [v2.4 - 聚合逻辑修正版]
            % ==
            % == 核心修改:
            % ==   - 新增透传 `L_hat_by_type`，确保上游函数可以获得
            % ==     唯一的、一致的按类型划分的劳动供给数据。
            % =========================================================================

            M_for_hh.r_mkt_t = r_mkt;
            M_for_hh.w_t = w_hat;
            M_for_hh.tr_per_hh = tr;
            M_for_hh.b_hat_by_type = b_hat_by_type_guess;

            [polS, valS] = household.VFI_solver(M_for_hh, paramS, cS);

            Dist = distribution.get_ss_dist(polS, paramS, cS, Z_ss_norm, bq_by_type);

            aggrS = aggregates.get_aggregates(Dist, polS, cS, paramS);

            macro_state = struct();
            macro_state.K_p_hat_tplus1_raw = aggrS.K_p_hat_tplus1_raw;
            macro_state.C_agg = aggrS.C_agg;
            macro_state.Bequest_generated_by_type_agg = aggrS.Bequest_generated_by_type_agg;
            macro_state.Bequest_generated_agg = aggrS.Bequest_generated_agg;
            macro_state.Total_Tax_Revenue = aggrS.Total_Tax_Revenue;
            macro_state.L_hat = aggrS.L_hat;
            macro_state.L_hat_by_type = aggrS.L_hat_by_type; % [新增]
        end





        function report_struct = display_national_accounts(ss, cS, paramS, Z_ss_norm, report_filename, print_flag, is_bgp_ss)
            % =========================================================================
            % == 函数: display_national_accounts
            % == 版本: [v9.5 - 核算变量与报告修正版]
            % ==
            % == 核心修改:
            % ==   - [!!! BUG修复 !!!] 重新定义了在重构中丢失的家庭部门总收入
            % ==     (Total_Household_Income) 和总储蓄 (Household_Savings) 变量。
            % ==   - [!!! BUG修复 !!!] 修正了新的分类型报告部分对未定义变量
            % ==     (如 Household_Savings_by_type) 的调用错误，改为从
            % ==     report_struct.by_type 中获取正确数据。
            % ==   - 整合并清理了报告逻辑，确保所有打印的变量都已预先计算并存入
            % ==     report_struct，增强了代码的稳健性和可读性。
            % =========================================================================

            if nargin < 6, print_flag = true; end
            if nargin < 5, report_filename = 'ss_national_accounts_report.txt'; end

            % --- 0. 准备工作 ---
            if print_flag
                try
                    fileID = fopen(report_filename, 'w', 'n', 'UTF-8');
                    if fileID == -1, error('无法创建或打开报告文件: %s', report_filename); end
                catch ME
                    warning('创建报告文件时出错，将输出到命令行。错误: %s', ME.message);
                    fileID = 1;
                end
            else, fileID = -1; end

            g_A_period = (1+cS.g_A_ss)^cS.time_Step-1;
            n_period = (1+cS.n_ss)^cS.time_Step-1;
            g_total_period = (1+g_A_period)*(1+n_period)-1;
            Y = ss.Y_from_production_hat;
            mass_total = sum(Z_ss_norm);

            % --- 1. 重新计算支持异质性类型的分布和策略 ---
            Dist_5D = SS.get_final_dist(ss, cS, paramS, Z_ss_norm, is_bgp_ss);
            polS_by_type = SS.get_final_polS(ss, cS, paramS, Z_ss_norm);

            % --- 2. 循环聚合所有类型的宏观流量 (总量 + 分类型) ---
            C_by_type = zeros(cS.nTypes, 1);
            L_by_type = zeros(cS.nTypes, 1);
            K_p_by_type = zeros(cS.nTypes, 1);
            PAYG_contributions_by_type = zeros(cS.nTypes, 1);
            PAYG_benefits_by_type = zeros(cS.nTypes, 1);
            Tax_regular_by_type = zeros(cS.nTypes, 1);
            pps_active = isfield(cS, 'pps_active') && cS.pps_active;
            PPS_contrib_by_type = zeros(cS.nTypes, 1);
            PPS_withdrawal_by_type = zeros(cS.nTypes, 1);
            PPS_tax_by_type = zeros(cS.nTypes, 1);

            for i_type = 1:cS.nTypes
                Dist_type = Dist_5D(:,:,:,:,i_type);
                polS_type = polS_by_type{i_type};
                b_hat_this_type = ss.b_hat_by_type(i_type);
                if i_type <= 2, theta_ss_type = cS.theta_path_urban(end);
                else, theta_ss_type = cS.theta_path_resident(end); end

                k_grid_full = repmat(cS.kGridV, [1, cS.nkpps, cS.nw_expanded]);
                kpps_grid_full = repmat(reshape(cS.kppsGridV, [1, cS.nkpps]), [cS.nk, 1, cS.nw_expanded]);

                for ia = 1:cS.aD_new
                    mass_dist_ia = Dist_type(:,:,:,ia);
                    if sum(mass_dist_ia, 'all') < 1e-30, continue; end

                    K_p_by_type(i_type) = K_p_by_type(i_type) + sum(k_grid_full .* mass_dist_ia, 'all');
                    if pps_active, K_p_by_type(i_type) = K_p_by_type(i_type) + sum(kpps_grid_full .* mass_dist_ia, 'all'); end

                    C_by_type(i_type) = C_by_type(i_type) + sum(polS_type(ia).c .* mass_dist_ia, 'all');
                    Tax_regular_by_type(i_type) = Tax_regular_by_type(i_type) + sum(polS_type(ia).tax_regular .* mass_dist_ia, 'all');

                    if ia <= cS.aR_new
                        L_by_type(i_type) = L_by_type(i_type) + sum(polS_type(ia).l .* mass_dist_ia, 'all');
                        labor_income_slice_full = ss.w_hat * cS.ageEff_by_type(i_type, ia) .* paramS.leGridV(1:cS.nw_expanded)';
                        PAYG_contributions_by_type(i_type) = PAYG_contributions_by_type(i_type) + sum(theta_ss_type .* labor_income_slice_full .* mass_dist_ia, 'all');
                        if pps_active
                            pps_contrib_slice_full = cS.pps_fixed * labor_income_slice_full;
                            PPS_contrib_by_type(i_type) = PPS_contrib_by_type(i_type) + sum(pps_contrib_slice_full .* mass_dist_ia, 'all');
                        end
                    else
                        PAYG_benefits_by_type(i_type) = PAYG_benefits_by_type(i_type) + sum(b_hat_this_type .* mass_dist_ia, 'all');
                        if pps_active
                            kpps_level = kpps_grid_full * (1 + ss.r_mkt);
                            pps_withdrawal_slice = cS.pps_withdrawal_rate * kpps_level;
                            pps_tax_slice = cS.pps_tax_rate_withdrawal * pps_withdrawal_slice;
                            PPS_withdrawal_by_type(i_type) = PPS_withdrawal_by_type(i_type) + sum(pps_withdrawal_slice .* mass_dist_ia, 'all');
                            PPS_tax_by_type(i_type) = PPS_tax_by_type(i_type) + sum(pps_tax_slice .* mass_dist_ia, 'all');
                        end
                    end
                end
            end

            % --- 3. 准备国民账户组件 ---
            PAYG_contributions_agg = sum(PAYG_contributions_by_type);
            PAYG_benefits_agg = sum(PAYG_benefits_by_type);
            PPS_contrib_agg = sum(PPS_contrib_by_type);
            PPS_withdrawal_agg = sum(PPS_withdrawal_by_type);
            PPS_tax_agg = sum(PPS_tax_by_type);
            General_Tax_Revenue = sum(Tax_regular_by_type) - PAYG_contributions_agg - PPS_tax_agg;

            I_p_demand_bgp = (g_total_period + cS.ddk) * ss.K_private_hat;
            I_g_display = cS.I_g_to_Y_ratio_ss * Y;
            G_c_display = cS.G_c_to_Y_ratio_ss * Y;
            I_p_accounting = Y - sum(C_by_type) - I_g_display - G_c_display;
            Total_Expenditure = sum(C_by_type) + I_p_accounting + I_g_display + G_c_display;
            resource_constraint_error = Y - Total_Expenditure;
            Depreciation_p = cS.ddk * ss.K_private_hat;
            Depreciation_g = cS.ddk_g * ss.K_public_hat;
            Total_Labor_Income = ss.w_hat * ss.L_hat;
            Total_Capital_Income_Gross = (ss.r_mkt + cS.ddk) * ss.K_private_hat;
            Public_Capital_Return_display = Y - Total_Labor_Income - Total_Capital_Income_Gross;

            % [!!! 核心修正: 重新定义家庭账户变量 !!!]
            Net_Capital_Income = ss.r_mkt * ss.K_private_hat;
            Total_Household_Income = Total_Labor_Income + Net_Capital_Income + PAYG_benefits_agg + PPS_withdrawal_agg + ss.Bequest_distributed_agg + ss.TR_distributed_agg;
            Household_Savings = Total_Household_Income - sum(C_by_type) - General_Tax_Revenue - PAYG_contributions_agg - PPS_contrib_agg - PPS_tax_agg;

            Total_General_Gov_Revenue = General_Tax_Revenue + PPS_tax_agg + Public_Capital_Return_display;
            Total_General_Gov_Expenditure = G_c_display + I_g_display + ss.TR_distributed_agg;
            general_budget_balance_error = Total_General_Gov_Revenue - Total_General_Gov_Expenditure;
            payg_balance_error = PAYG_contributions_agg - PAYG_benefits_agg;
            payg_balance_by_type_error = PAYG_contributions_by_type - PAYG_benefits_by_type;

            % --- 4. 创建并填充 report_struct ---
            report_struct = ss;
            report_struct.Y = Y; 
            report_struct.K_p = ss.K_private_hat; 
            report_struct.K_g = ss.K_public_hat;
            report_struct.L = ss.L_hat;
            report_struct.r_period = ss.r_mkt; 
            report_struct.w = ss.w_hat;
            report_struct.r_annual = ((1+ss.r_mkt)^(1/cS.time_Step) - 1)*100;
            report_struct.C = sum(C_by_type); report_struct.I_p_acct = I_p_accounting;
            report_struct.I_p_bgp = I_p_demand_bgp; report_struct.I_g = I_g_display; report_struct.G_c = G_c_display;
            report_struct.Depreciation_p = Depreciation_p; report_struct.Depreciation_g = Depreciation_g;
            report_struct.Kp_Y_ratio = ss.K_private_hat / Y; report_struct.C_Y_ratio = (report_struct.C / Y)*100;
            report_struct.Ip_Y_ratio = (I_p_accounting / Y)*100;
            report_struct.Total_Labor_Income = Total_Labor_Income;
            report_struct.Net_Capital_Income = Net_Capital_Income;
            report_struct.Household_Savings = Household_Savings; % [修正] 添加总储蓄
            report_struct.Public_Capital_Return = Public_Capital_Return_display;
            report_struct.RC_Error = resource_constraint_error; report_struct.Invest_Gap = I_p_accounting - I_p_demand_bgp;
            report_struct.RC_Error_pct = (resource_constraint_error / Y)*100; report_struct.Invest_Gap_pct = (report_struct.Invest_Gap / Y)*100;
            report_struct.TR = ss.TR_distributed_agg; report_struct.Tax_General = General_Tax_Revenue;
            report_struct.Gov_Budget_Balance_Error = general_budget_balance_error;
            report_struct.PAYG_contrib = PAYG_contributions_agg; report_struct.PAYG_benefit = PAYG_benefits_agg;
            report_struct.PAYG_Balance_Error = payg_balance_error;
            report_struct.PAYG_contrib_by_type = PAYG_contributions_by_type;
            report_struct.PAYG_benefit_by_type = PAYG_benefits_by_type;
            report_struct.PAYG_balance_by_type_error = payg_balance_by_type_error;
            if pps_active
                report_struct.PPS_contrib = PPS_contrib_agg;
                report_struct.PPS_withdrawal = PPS_withdrawal_agg;
                report_struct.PPS_tax = PPS_tax_agg;
            end
            by_type = struct();
            by_type.Pop_Share = cS.type_weights * 100;
            by_type.C = C_by_type; by_type.L = L_by_type; by_type.K_p = K_p_by_type;
            by_type.LaborIncome = L_by_type * ss.w_hat;
            by_type.CapitalIncome = K_p_by_type * ss.r_mkt;
            by_type.Bequest_received = ss.Bequest_generated_by_type_agg;
            by_type.TR_received = ss.TR_distributed_agg * cS.type_weights;
            by_type.PAYG_benefit = PAYG_benefits_by_type;
            by_type.PAYG_contrib = PAYG_contributions_by_type;
            by_type.Tax_General = Tax_regular_by_type - PAYG_contributions_by_type - PPS_tax_by_type;
            by_type.PPS_contrib = PPS_contrib_by_type; by_type.PPS_withdrawal = PPS_withdrawal_by_type; by_type.PPS_tax = PPS_tax_by_type;
            by_type.TotalIncome = by_type.LaborIncome + by_type.CapitalIncome + by_type.Bequest_received + by_type.TR_received + by_type.PAYG_benefit + by_type.PPS_withdrawal;
            by_type.Savings = by_type.TotalIncome - by_type.C - by_type.Tax_General - by_type.PAYG_contrib - by_type.PPS_contrib - by_type.PPS_tax;
            report_struct.by_type = by_type;

            % --- 5. 报告打印与最终检验 ---
            if print_flag
                type_names = {'高收城镇', '中低收城镇', '高收居民', '中低收居民'};
                fprintfAndLog(fileID, '\n\n=========================================================================================================\n');
                if pps_active
                    fprintfAndLog(fileID, '###            稳 态 国 民 经 济 核 算 报 告 (v9.5 - 核算修正版)                    ###\n');
                else
                    fprintfAndLog(fileID, '###                 稳 态 国 民 经 济 核 算 报 告 (v9.5 - 核算修正版)                      ###\n');
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
                fprintfAndLog(fileID, '      公共资本存量 (K_g) ........................ : %15.6f\n', report_struct.K_g);
                fprintfAndLog(fileID, '      有效劳动需求 (L) .......................... : %15.6f\n', report_struct.L);
                fprintfAndLog(fileID, '   B. 要素报酬 (流量):\n');
                fprintfAndLog(fileID, '      支付工资总额 (w*L) ........................ : %15.6f  (Y的 %6.2f %%)\n', report_struct.Total_Labor_Income, (report_struct.Total_Labor_Income/Y)*100);
                fprintfAndLog(fileID, '      私人资本总回报 ((r+d)*K_p) ................ : %15.6f  (Y的 %6.2f %%)\n', Total_Capital_Income_Gross, (Total_Capital_Income_Gross/Y)*100);
                fprintfAndLog(fileID, '      公共资本回报 (Y - wL - (r+d)K_p) ........ : %15.6f  (Y的 %6.2f %%)\n', report_struct.Public_Capital_Return, (report_struct.Public_Capital_Return/Y)*100);

                fprintfAndLog(fileID, '\n\n--- [ III. 家庭部门 ] ---\n');
                fprintfAndLog(fileID, '   A. 收入来源 (流量):\n');
                fprintfAndLog(fileID, '      劳动收入 (来自厂商) ....................... : %15.6f\n', report_struct.Total_Labor_Income);
                fprintfAndLog(fileID, '      净资本收入 (r*K_p) ........................ : %15.6f\n', report_struct.Net_Capital_Income);
                fprintfAndLog(fileID, '      养老金福利 ................................ : %15.6f\n', report_struct.PAYG_benefit);
                if pps_active, fprintfAndLog(fileID, '      PPS账户提取 ............................... : %15.6f\n', report_struct.PPS_withdrawal); end
                fprintfAndLog(fileID, '      收到的遗赠 ................................ : %15.6f\n', ss.Bequest_distributed_agg);
                fprintfAndLog(fileID, '      收到的转移支付 ............................ : %15.6f\n', report_struct.TR);
                fprintfAndLog(fileID, '      --------------------------------------------------------------------------------------------------\n');
                fprintfAndLog(fileID, '      >>> 家庭总收入 ............................ : %15.6f\n', Total_Household_Income);
                fprintfAndLog(fileID, '   B. 收入用途 (流量):\n');
                fprintfAndLog(fileID, '      消费支出 (C) .............................. : %15.6f\n', report_struct.C);
                fprintfAndLog(fileID, '      支付的一般税收 ............................ : %15.6f\n', report_struct.Tax_General);
                fprintfAndLog(fileID, '      养老金缴费 ................................ : %15.6f\n', report_struct.PAYG_contrib);
                if pps_active, fprintfAndLog(fileID, '      PPS账户缴费 (储蓄) ...................... : %15.6f\n', report_struct.PPS_contrib); fprintfAndLog(fileID, '      PPS提取税 ................................. : %15.6f\n', report_struct.PPS_tax); end
                fprintfAndLog(fileID, '      净储蓄 .................................... : %15.6f\n', report_struct.Household_Savings);
                fprintfAndLog(fileID, '      --------------------------------------------------------------------------------------------------\n');
                fprintfAndLog(fileID, '      >>> 家庭总支出与储蓄 ...................... : %15.6f\n', Total_Household_Income);

                fprintfAndLog(fileID, '\n\n--- [ IV. 政府与养老金体系 ] ---\n');
                fprintfAndLog(fileID, '   A. 一般政府预算:\n');
                fprintfAndLog(fileID, '      (+) 一般税收收入 .......................... : %15.6f\n', report_struct.Tax_General);
                if pps_active, fprintfAndLog(fileID, '      (+) PPS提取税 ............................. : %15.6f\n', report_struct.PPS_tax); end
                fprintfAndLog(fileID, '      (+) 公共资本的隐性回报 .................... : %15.6f\n', report_struct.Public_Capital_Return);
                fprintfAndLog(fileID, '      --------------------------------------------------------------------------------------------------\n');
                fprintfAndLog(fileID, '      (=) 政府总收入 ............................ : %15.6f\n', Total_General_Gov_Revenue);
                fprintfAndLog(fileID, '      (-) 政府消费 (G_c) ........................ : %15.6f\n', report_struct.G_c);
                fprintfAndLog(fileID, '      (-) 政府总投资 (I_g) ...................... : %15.6f\n', report_struct.I_g);
                fprintfAndLog(fileID, '      (-) 对家庭的转移支付 (TR) ................. : %15.6f\n', report_struct.TR);
                fprintfAndLog(fileID, '      --------------------------------------------------------------------------------------------------\n');
                fprintfAndLog(fileID, '      (=) 政府总支出 ............................ : %15.6f\n', Total_General_Gov_Expenditure);
                fprintfAndLog(fileID, '      >>> 预算平衡误差 (收入 - 支出) ............ : %15.4e\n', report_struct.Gov_Budget_Balance_Error);
                if abs(report_struct.Gov_Budget_Balance_Error/Y) < 1e-9, fprintfAndLog(fileID, '          核验结果: ✅ 政府一般预算按设计严格平衡。\n');
                else, fprintfAndLog(fileID, '          核验结果: ⚠️ 政府一般预算未按设计平衡！\n'); end
                fprintfAndLog(fileID, '   B. 现收现付(PAYG)养老金体系 (按类型细分):\n');
                fprintfAndLog(fileID, '      %-15s | %15s | %15s | %15s | %15s\n', '类型', '福利水平(b_hat)', '总缴费', '总福利', '余额(误差)');
                fprintfAndLog(fileID, '      %s\n', repmat('-', 1, 72));
                for i_type = 1:cS.nTypes
                    fprintfAndLog(fileID, '      %-15s | %15.6f | %15.6f | %15.6f | %15.4e\n', ...
                        type_names{i_type}, ss.b_hat_by_type(i_type), ...
                        PAYG_contributions_by_type(i_type), PAYG_benefits_by_type(i_type), payg_balance_by_type_error(i_type));
                end
                fprintfAndLog(fileID, '      %s\n', repmat('-', 1, 72));
                fprintfAndLog(fileID, '      (+) 养老金总缴费 .......................... : %15.6f\n', PAYG_contributions_agg);
                fprintfAndLog(fileID, '      (-) 养老金总福利 .......................... : %15.6f\n', PAYG_benefits_agg);
                fprintfAndLog(fileID, '      >>> PAYG体系总余额 (缴费 - 福利) .......... : %15.4e\n', payg_balance_error);
                if abs(payg_balance_error/Y) < 1e-9, fprintfAndLog(fileID, '          核验结果: ✅ PAYG体系按设计严格平衡。\n');
                else, fprintfAndLog(fileID, '          核验结果: ⚠️ PAYG体系未按设计平衡！\n'); end

                % [!!! 修正部分 !!!]
                fprintfAndLog(fileID, '\n\n--- [ V. 按家庭类型细分的宏观账户 ] ---\n');
                fprintfAndLog(fileID, '   A. 各类型绝对量贡献:\n');
                fprintfAndLog(fileID, '      %-20s | %15s | %15s | %15s | %15s | %15s\n', '变量', type_names{1}, type_names{2}, type_names{3}, type_names{4}, '经济体总量');
                fprintfAndLog(fileID, '      %s\n', repmat('-', 1, 110));
                fprintfAndLog(fileID, '      %-20s | %15.4f | %15.4f | %15.4f | %15.4f | %15.4f\n', '人口权重 (%)', report_struct.by_type.Pop_Share, 100.0);
                fprintfAndLog(fileID, '      %-20s | %15.4f | %15.4f | %15.4f | %15.4f | %15.4f\n', '私人资本 (K_p)', report_struct.by_type.K_p, sum(report_struct.by_type.K_p));
                fprintfAndLog(fileID, '      %-20s | %15.4f | %15.4f | %15.4f | %15.4f | %15.4f\n', '有效劳动 (L)', report_struct.by_type.L, sum(report_struct.by_type.L));
                fprintfAndLog(fileID, '      %-20s | %15.4f | %15.4f | %15.4f | %15.4f | %15.4f\n', '总消费 (C)', report_struct.by_type.C, sum(report_struct.by_type.C));
                fprintfAndLog(fileID, '      %-20s | %15.4f | %15.4f | %15.4f | %15.4f | %15.4f\n', '总储蓄 (S)', report_struct.by_type.Savings, sum(report_struct.by_type.Savings));
                fprintfAndLog(fileID, '      %-20s | %15.4f | %15.4f | %15.4f | %15.4f | %15.4f\n', '总收入 (Income)', report_struct.by_type.TotalIncome, sum(report_struct.by_type.TotalIncome));
                fprintfAndLog(fileID, '      %s\n', repmat('-', 1, 110));
                fprintfAndLog(fileID, '\n   B. 各类型相对份额 (%% of Total):\n');
                fprintfAndLog(fileID, '      %-20s | %15s | %15s | %15s | %15s | %15s\n', '变量', type_names{1}, type_names{2}, type_names{3}, type_names{4}, '合计');
                fprintfAndLog(fileID, '      %s\n', repmat('-', 1, 110));
                fprintfAndLog(fileID, '      %-20s | %15.2f%% | %15.2f%% | %15.2f%% | %15.2f%% | %15.2f%%\n', '私人资本 (K_p)', report_struct.by_type.K_p/sum(report_struct.by_type.K_p)*100, 100.0);
                fprintfAndLog(fileID, '      %-20s | %15.2f%% | %15.2f%% | %15.2f%% | %15.2f%% | %15.2f%%\n', '有效劳动 (L)', report_struct.by_type.L/sum(report_struct.by_type.L)*100, 100.0);
                fprintfAndLog(fileID, '      %-20s | %15.2f%% | %15.2f%% | %15.2f%% | %15.2f%% | %15.2f%%\n', '总消费 (C)', report_struct.by_type.C/sum(report_struct.by_type.C)*100, 100.0);
                fprintfAndLog(fileID, '      %-20s | %15.2f%% | %15.2f%% | %15.2f%% | %15.2f%% | %15.2f%%\n', '总储蓄 (S)', report_struct.by_type.Savings/sum(report_struct.by_type.Savings)*100, 100.0);
                fprintfAndLog(fileID, '      %-20s | %15.2f%% | %15.2f%% | %15.2f%% | %15.2f%% | %15.2f%%\n', '总收入 (Income)', report_struct.by_type.TotalIncome/sum(report_struct.by_type.TotalIncome)*100, 100.0);
                fprintfAndLog(fileID, '      %s\n', repmat('-', 1, 110));

                fprintfAndLog(fileID, '\n\n--- [ VI. 核心均衡条件检验 (来自求解器的最终误差) ] ---\n');
                tr_per_capita_final = ss.TR_distributed_agg / mass_total;
                b_hat_by_type_final = ss.b_hat_by_type;
                bq_by_type_final = ss.Bequest_generated_by_type_agg;
                x_final_vec = [ss.r_mkt; ss.w_hat; tr_per_capita_final; b_hat_by_type_final; bq_by_type_final];
                final_errors = SS.system_of_equations(x_final_vec, Z_ss_norm, cS, paramS, struct('A',1.0), is_bgp_ss);
                fprintfAndLog(fileID, '   1. 资本市场出清 (r_guess - r_new) ......... : %15.4e\n', final_errors(1));
                fprintfAndLog(fileID, '   2. 劳动市场出清 (w_guess - w_new) ....... : %15.4e\n', final_errors(2));
                fprintfAndLog(fileID, '   3. 政府预算平衡 (TR_guess - TR_new) ..... : %15.4e\n', final_errors(3));
                base_idx = 4;
                for i_type = 1:cS.nTypes
                    fprintfAndLog(fileID, '   %d.%d. PAYG平衡 (类型 %d) .................. : %15.4e\n', base_idx, i_type, i_type, final_errors(base_idx-1+i_type));
                end
                base_idx = base_idx + cS.nTypes;
                for i_type = 1:cS.nTypes
                    fprintfAndLog(fileID, '   %d.%d. 遗赠市场出清 (类型 %d) ............ : %15.4e\n', base_idx, i_type, i_type, final_errors(base_idx-1+i_type));
                end
                base_idx = base_idx + cS.nTypes;
                if (length(final_errors) >= base_idx)
                    fprintfAndLog(fileID, '   %d. BGP资本积累 (K_t - K_t+1_adj) ........ : %15.4e\n', base_idx, final_errors(base_idx));
                end

                fprintfAndLog(fileID, '\n\n=========================================================================================================\n');
                fprintfAndLog(fileID, '###                                      报 告 结 束                                      ###\n');
                fprintfAndLog(fileID, '=========================================================================================================\n\n');

                if fileID > 1, fclose(fileID); fprintf('\n报告已成功保存至: %s\n', report_filename); end
            end

            function fprintfAndLog(fileID, formatSpec, varargin)
                if print_flag
                    fprintf(1, formatSpec, varargin{:});
                    if fileID > 1, fprintf(fileID, formatSpec, varargin{:}); end
                end
            end
        end
        
        function report_struct = display_national_accounts_ss0(ss, cS, paramS, Z_ss_norm, report_filename, print_flag)
            % =========================================================================
            % == 函数: display_national_accounts_ss0
            % == 版本: [v1.1 - 增加资本结构汇报]
            % ==
            % == 核心修改:
            % ==   - 明确计算公共资本存量 K_g 及其与产出的比率 K_g/Y。
            % ==   - 在报告中新增一个“资本结构”部分，分别汇报 K_p/Y, K_g/Y,
            % ==     以及加总后的总资本产出比 Total K/Y。
            % =========================================================================

            if nargin < 6, print_flag = true; end
            if nargin < 5, report_filename = 'ss0_national_accounts_report.txt'; end

            % --- 0. 准备工作 ---
            if print_flag
                try, fileID = fopen(report_filename, 'w', 'n', 'UTF-8');
                catch ME, warning('无法创建报告文件 %s。错误: %s', report_filename, ME.message); fileID = 1; end
            else, fileID = -1; end

            Y = ss.Y_from_production_hat;
            n_period = (1 + cS.n_ss)^cS.time_Step - 1;

            % --- 1. 重新聚合关键流量 (与之前相同) ---
            Dist_5D = SS.get_final_dist_ss0(ss, cS, paramS, Z_ss_norm);
            polS_by_type = SS.get_final_polS_ss0(ss, cS, paramS, Z_ss_norm); % 使用ss0的辅助函数
            
            C_agg = 0;
            for i_type = 1:cS.nTypes
                Dist_type = Dist_5D(:,:,:,:,i_type);
                polS_type = polS_by_type{i_type};
                for ia = 1:cS.aD_new
                    mass_dist_ia = Dist_type(:,:,:,ia);
                    if sum(mass_dist_ia,'all') > 1e-30
                       C_agg = C_agg + sum(polS_type(ia).c .* mass_dist_ia, 'all');
                    end
                end
            end

            % --- 2. 准备国民账户组件 (与之前相同) ---
            I_g_display = ss.I_g;
            G_c_display = ss.G_c;
            I_p_accounting = Y - C_agg - I_g_display - G_c_display;
            Total_Expenditure = C_agg + I_p_accounting + I_g_display + G_c_display;
            resource_constraint_error = Y - Total_Expenditure;
            K_p_tplus1_raw = ss.K_p_hat_tplus1_raw;
            K_p_t = ss.K_private_hat;
            I_p_behavioral = K_p_tplus1_raw - (1 - cS.ddk) * K_p_t;

            % --- 3. [核心修改] 创建并填充 report_struct (增加资本结构) ---
            report_struct = ss;
            report_struct.Y = Y; 
            report_struct.K_p = ss.K_private_hat; 
            report_struct.K_g = ss.K_public_hat; % 从ss结构体中直接获取
            report_struct.L = ss.L_hat;
            report_struct.r_annual = ((1+ss.r_mkt)^(1/cS.time_Step) - 1)*100;
            report_struct.C = C_agg; 
            report_struct.I_p_acct = I_p_accounting;
            
            % 计算并添加资本结构比率
            report_struct.Kp_Y_ratio = ss.K_private_hat / Y;
            report_struct.Kg_Y_ratio = ss.K_public_hat / Y;
            report_struct.TotalK_Y_ratio = report_struct.Kp_Y_ratio + report_struct.Kg_Y_ratio;

            % --- 4. [核心修改] 报告打印 (增加资本结构部分) ---
            if print_flag
                fprintfAndLog(fileID, '\n\n=================================================================================\n');
                fprintfAndLog(fileID, '###      初始时期 (t=1) 宏观经济快照 (资本市场非均衡版)      ###\n');
                fprintfAndLog(fileID, '=================================================================================\n');

                fprintfAndLog(fileID, '\n--- [ I. 宏观总量与价格 ] ---\n');
                fprintfAndLog(fileID, '   国内生产总值 (Y) .......................... : %15.6f\n', report_struct.Y);
                fprintfAndLog(fileID, '   有效劳动 (L) .............................. : %15.6f\n', report_struct.L);
                fprintfAndLog(fileID, '   真实利率 (r, 年化) ........................ : %15.4f %%\n', report_struct.r_annual);
                fprintfAndLog(fileID, '   真实工资率 (w) ............................ : %15.6f\n', report_struct.w_hat);
                
                fprintfAndLog(fileID, '\n--- [ II. 资本结构 (t=1) ] ---\n');
                fprintfAndLog(fileID, '   当期私人资本存量 (K_p, 内生) ............ : %15.6f\n', report_struct.K_p);
                fprintfAndLog(fileID, '   当期公共资本存量 (K_g, 政策决定) ........ : %15.6f\n', report_struct.K_g);
                fprintfAndLog(fileID, '   --------------------------------------------------------------------------\n');
                fprintfAndLog(fileID, '   私人资本/产出比 (K_p/Y) ................... : %15.6f\n', report_struct.Kp_Y_ratio);
                fprintfAndLog(fileID, '   公共资本/产出比 (K_g/Y) ................... : %15.6f\n', report_struct.Kg_Y_ratio);
                fprintfAndLog(fileID, '   >>> 总资本/产出比 ((K_p+K_g)/Y) ........... : %15.6f\n', report_struct.TotalK_Y_ratio);

                fprintfAndLog(fileID, '\n--- [ III. 商品市场流量 (t=1) ] ---\n');
                fprintfAndLog(fileID, '   A. 国民支出 (Y = C + I_p + I_g + G_c):\n');
                fprintfAndLog(fileID, '      (+) 私人消费 (C) .......................... : %15.6f  (%6.2f %% of Y)\n', report_struct.C, (report_struct.C/Y)*100);
                fprintfAndLog(fileID, '      (+) 私人总投资 (I_p, 会计值) ............ : %15.6f  (%6.2f %% of Y)\n', report_struct.I_p_acct, (report_struct.I_p_acct/Y)*100);
                fprintfAndLog(fileID, '      (+) 公共总投资 (I_g) ...................... : %15.6f  (%6.2f %% of Y)\n', I_g_display, (I_g_display/Y)*100);
                fprintfAndLog(fileID, '      (+) 政府消费 (G_c) ........................ : %15.6f  (%6.2f %% of Y)\n', G_c_display, (G_c_display/Y)*100);
                fprintfAndLog(fileID, '      --------------------------------------------------------------------------\n');
                fprintfAndLog(fileID, '      (=) 总支出 ................................ : %15.6f\n', Total_Expenditure);
                fprintfAndLog(fileID, '      >>> 资源约束误差 (Y - 总支出) ............. : %15.4e\n', resource_constraint_error);
                if abs(resource_constraint_error/Y) < 1e-9, fprintfAndLog(fileID, '          核验结果: ✅ 商品市场流量平衡。\n'); end

                fprintfAndLog(fileID, '\n   B. 资本积累动态 (t=1 -> t=2):\n');
                fprintfAndLog(fileID, '      家庭意愿投资 (I_p, 行为值) .............. : %15.6f\n', I_p_behavioral);
                fprintfAndLog(fileID, '          (定义: K_t+1 - (1-d)K_t)\n');
                fprintfAndLog(fileID, '      会计核算投资 (I_p, 剩余值) .............. : %15.6f\n', I_p_accounting);
                fprintfAndLog(fileID, '          (定义: Y - C - G)\n');
                fprintfAndLog(fileID, '      >>> 投资不一致性 .......................... : %15.4e\n', I_p_accounting - I_p_behavioral);
                fprintfAndLog(fileID, '          说明: 此差异是正常的，它反映了t=1的非均衡性质。该差异将在\n');
                fprintfAndLog(fileID, '                转轨路径的后续时期(t>=2)中被动态吸收和消除。\n');
                
                fprintfAndLog(fileID, '\n\n=================================================================================\n');
                fprintfAndLog(fileID, '###                            报 告 结 束                            ###\n');
                fprintfAndLog(fileID, '=================================================================================\n\n');

                if fileID > 1, fclose(fileID); fprintf('\n初始时期报告已成功保存至: %s\n', report_filename); end
            end

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
            % == 版本: [v1.6 - b_hat 字段修正版]
            % ==
            % == 核心修改:
            % ==   - [!!! BUG修复 !!!] 修正了 x_eq_final 向量的构造。
            % ==   - 不再访问已不存在的 ss.b_hat，而是读取新的向量字段 ss.b_hat_by_type。
            % ==   - 确保构造的 x_eq_final 向量与 system_of_equations 的输入格式完全匹配。
            % =========================================================================
            mass_total = sum(Z_ss_norm);
            tr_per_capita = ss.TR_distributed_agg / mass_total;
            bq_by_type_final = ss.Bequest_generated_by_type_agg;

            % [核心修正] 从 ss 中获取按类型划分的 b_hat 向量
            b_hat_by_type_final = ss.b_hat_by_type;

            % [核心修正] 构造与求解器一致的、正确长度的 x_eq_final 向量
            x_eq_final = [ss.r_mkt; ss.w_hat; tr_per_capita; b_hat_by_type_final; bq_by_type_final];

            [~, ~, Dist] = SS.system_of_equations(x_eq_final, Z_ss_norm, cS, paramS, struct('A',1.0), is_bgp_ss);
            
        end

        function polS = get_final_polS(ss, cS, paramS, Z_ss_norm)
            % =========================================================================
            % == 函数: get_final_polS
            % == 版本: [v1.5 - b_hat 字段修正版]
            % ==
            % == 核心修改:
            % ==   - [!!! BUG修复 !!!] 在构建 M_for_hh 时，不再访问已不存在的 ss.b_hat。
            % ==   - 现在正确地将 ss.b_hat_by_type 向量赋给 M_for_hh.b_hat_by_type，
            % ==     以匹配 VFI_solver 的输入要求。
            % =========================================================================
            M_for_hh.r_mkt_t = ss.r_mkt;
            M_for_hh.w_t = ss.w_hat;
            M_for_hh.tr_per_hh = ss.TR_distributed_agg / sum(Z_ss_norm);

            % [核心修正] 直接使用最终的、按类型划分的 b_hat 向量
            M_for_hh.b_hat_by_type = ss.b_hat_by_type;

            [polS, ~] = household.VFI_solver(M_for_hh, paramS, cS);
        end    
    
    
            function Dist = get_final_dist_ss0(ss, cS, paramS, Z_ss_norm)
            mass_total = sum(Z_ss_norm);
            tr_per_capita = ss.TR_distributed_agg / mass_total;
            bq_by_type_final = ss.Bequest_generated_by_type_agg;
            b_hat_by_type_final = ss.b_hat_by_type;
            x_eq_final = [ss.r_mkt; ss.w_hat; tr_per_capita; b_hat_by_type_final; bq_by_type_final];
            
            params_ext = struct('A', 1.0, 'g_A_ss', cS.g_A_ss, 'n_ss', cS.n_ss, 'Z', Z_ss_norm);
            [~, ~, Dist] = SS.system_of_equations_ss0(x_eq_final, Z_ss_norm, cS, paramS, params_ext, true);
        end

        function polS = get_final_polS_ss0(ss, cS, paramS, Z_ss_norm) % 虽然内容一样，但为了清晰分开
            M_for_hh.r_mkt_t = ss.r_mkt;
            M_for_hh.w_t = ss.w_hat;
            M_for_hh.tr_per_hh = ss.TR_distributed_agg / sum(Z_ss_norm);
            M_for_hh.b_hat_by_type = ss.b_hat_by_type;

            [polS, ~] = household.VFI_solver(M_for_hh, paramS, cS);
        end
    
    
    end



end