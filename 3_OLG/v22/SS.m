classdef SS
    methods (Static)
        
        function [ss, Dist, polS, valS] = solve_steady_state(cS, paramS, params_ext, verbose, is_bgp_ss)
            % =========================================================================
            % == 函数: solve_steady_state
            % == 版本: [v2.2 - 统一PAYG调整因子版]
            % ==
            % == 核心修改:
            % ==   - 将求解变量从 [r, w, tr, b_hat_h(1..n), bq_pc]
            % ==     修改为 [r, w, tr, adj, bq_pc]。
            % ==   - adj 是一个统一的调整因子，用于确保整个PAYG体系收支平衡。
            % ==   - b_hat_h 不再是求解变量，而是基于 w_hat 和 adj 计算得出。
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

            if ~isfield(cS, 'num_hh_types') || cS.num_hh_types <= 0
                error('未在cS中定义家庭类型数量(num_hh_types)或其值无效。');
            end
            nH = cS.num_hh_types;

            if ~cS.pps_active
                cS.nkpps = 1; cS.n_pps_rate_grid = 1; 
            end

            % --- 设置初始猜测值 (变量数: 3 + 1 + 1 = 5) ---
            r_guess_init = 0.04;
            w_guess_init = 1.5;
            TR_guess_init = 0.0;
            adj_guess_init = 1.0; % 调整因子初始猜测为1
            K_guess_for_BQ = 2.5;
            BQ_total_pc_guess_init = mean(1-cS.s_pathV) * K_guess_for_BQ / sum(Z_ss_norm,'all');

            % x_guess = [r, w, tr, adj, bq_total_pc]
            x0 = [r_guess_init; w_guess_init; TR_guess_init; adj_guess_init; BQ_total_pc_guess_init];
            
            % --- 设置边界 ---
            lb = [-0.1; 1e-8; -Inf; 0; 1e-8]; % adj 必须为非负
            ub = [Inf; Inf; Inf; Inf; Inf];

            system_wrapper = @(x) SS.system_of_equations(x, Z_ss_norm, cS, paramS, params_ext, is_bgp_ss);
  
            options = optimoptions('lsqnonlin', 'Display', verbose, ...
            'FunctionTolerance', 1e-12, 'StepTolerance', 1e-12, 'OptimalityTolerance', 1e-12, ...
            'MaxIterations', 500); % ,'Algorithm','levenberg-marquardt'

            [x_eq, ~, ~, exitflag] = lsqnonlin(system_wrapper, x0, lb, ub, options);

            [~, ss, Dist, polS, valS] = SS.system_of_equations(x_eq, Z_ss_norm, cS, paramS, params_ext, is_bgp_ss);
        end

        function [F_error, ss, Dist, polS, valS] = system_of_equations(x_guess, Z_ss_norm, cS, paramS, params_ext, is_bgp_ss)
            % =========================================================================
            % == 函数: system_of_equations
            % == 版本: [v5.1 - 福利公式滞后工资版]
            % ==
            % == 核心修改:
            % ==   - 在稳态下，w_hat(t-1) = w_hat(t)，所以直接用 w_guess
            % ==     同时作为当期工资和计算福利基数的工资，保持逻辑一致性。
            % =========================================================================
            nH = cS.num_hh_types;

            % --- 流程 0: 解包猜测向量 ---
            r_guess = x_guess(1);
            w_guess = x_guess(2);
            tr_guess = x_guess(3);
            adj_guess = x_guess(4); 
            bq_guess = x_guess(5); 

            % --- 流程 1: 初始化聚合变量与存储容器 ---
            L_supply_agg_total = 0;
            K_supply_agg_t_total = 0;
            K_supply_agg_tplus1_raw_total = 0;
            Total_Tax_Revenue_total = 0;
            Bequest_generated_raw_total = 0;

            Dist = cell(nH, 1);
            polS = cell(nH, 1);
            valS = cell(nH, 1);
            ss_h = cell(nH, 1);
            aggr_out_h =  cell(nH, 1);

            % [核心修改] 在稳态下，福利基数也由当期(即稳态)工资决定
            b_hat_formula_h = SS.calculate_formula_benefits(w_guess, cS);

            % --- 流程 2: [核心循环] 遍历各类家庭，求解微观并聚合 ---
            for h = 1:nH
                cS_h = cS;
                paramS_h = paramS;
                if isfield(cS_h, 'theta_path_h'), cS_h.theta_path = cS.theta_path_h(h, end); end
                if isfield(cS_h, 'ageEffV_new_h'), cS_h.ageEffV_new = cS.ageEffV_new_h(:,h); end

                b_hat_for_vfi = adj_guess * b_hat_formula_h(h);

                [aggr_out_h{h}, Dist{h}, polS{h}, valS{h}] = SS.run_micro_to_macro(r_guess, w_guess, bq_guess, tr_guess, b_hat_for_vfi, Z_ss_norm(:, h), cS_h, paramS_h, params_ext);

                L_supply_agg_total = L_supply_agg_total + aggr_out_h{h}.L_hat;
                K_supply_from_dist_h = aggregates.get_aggregates_from_dist(Dist{h}, cS_h);
                K_supply_agg_t_total = K_supply_agg_t_total + K_supply_from_dist_h.K_p_hat_total;
                K_supply_agg_tplus1_raw_total = K_supply_agg_tplus1_raw_total + aggr_out_h{h}.K_p_hat_tplus1_raw;
                Bequest_generated_raw_total = Bequest_generated_raw_total + aggr_out_h{h}.Bequest_generated_agg;
                Total_Tax_Revenue_total = Total_Tax_Revenue_total + aggr_out_h{h}.Total_Tax_Revenue;

                ss_h{h}.L_hat = aggr_out_h{h}.L_hat;
                ss_h{h}.K_hat = K_supply_from_dist_h.K_p_hat_total;
                ss_h{h}.C_agg = aggr_out_h{h}.C_agg;
            end

            % --- 流程 3: [宏观层面] 市场出清 ---
            n_period = (1 + cS.n_ss)^cS.time_Step - 1;

            g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
            g_total_period = (1 + g_A_period) * (1 + n_period) - 1;
            
            Kg_Y_ratio_factor = cS.I_g_to_Y_ratio_ss / ( g_total_period + cS.ddk_g); % 
            % if ~is_bgp_ss
            %     Kg_Y_ratio_factor = cS.I_g_to_Y_ratio_ss / ( cS.ddk_g); %
            % end
            Y_hat_calc = (K_supply_agg_t_total^cS.alpha * L_supply_agg_total^(1 - cS.alpha - cS.gamma) * Kg_Y_ratio_factor^cS.gamma)^(1 / (1 - cS.gamma));
            K_g_hat_calc = Kg_Y_ratio_factor * Y_hat_calc;
            prices_new = firm.get_prices_at_t(K_supply_agg_t_total, K_g_hat_calc, L_supply_agg_total, cS);
            r_new = prices_new.r_mkt_t;
            w_new = prices_new.w_hat_t;

            error_r = r_guess - r_new;
            error_w = w_guess - w_new;

            PAYG_contributions_total = 0;
            PAYG_formula_benefits_total = 0;
            for h = 1:nH
                payg_contrib_h = cS.theta_path_h(h, end) * w_guess * aggr_out_h{h}.L_hat;
                PAYG_contributions_total = PAYG_contributions_total + payg_contrib_h;

                mass_retirees_h = sum(Z_ss_norm((cS.aR_new+1):end, h));
                formula_benefit_h = b_hat_formula_h(h) * mass_retirees_h;
                PAYG_formula_benefits_total = PAYG_formula_benefits_total + formula_benefit_h;
            end

            total_benefits_paid_out = adj_guess * PAYG_formula_benefits_total;
            error_payg = PAYG_contributions_total - total_benefits_paid_out;

            General_Tax_Revenue_calc = Total_Tax_Revenue_total - PAYG_contributions_total;
            Public_Capital_Return_calc = prices_new.Y_hat_t - (w_new * L_supply_agg_total) - ((r_new + cS.ddk) * K_supply_agg_t_total);
            Total_General_Gov_Revenue_calc = General_Tax_Revenue_calc + Public_Capital_Return_calc;
            G_consumption_spending = cS.G_c_to_Y_ratio_ss * prices_new.Y_hat_t;
            G_investment_spending = cS.I_g_to_Y_ratio_ss * prices_new.Y_hat_t;
            mass_total = sum(Z_ss_norm(:));
            tr_new = (Total_General_Gov_Revenue_calc - (G_consumption_spending + G_investment_spending)) / mass_total;
            error_tr = tr_guess - tr_new;

                        % [!!! 核心修改 !!!] 为遗赠资产计入利息
            % bequest_generated_adj = Bequest_generated_raw_total * (1 + r_guess) / (1 + n_period);
            
            bequest_generated_adj = Bequest_generated_raw_total  / (1 + n_period);
            % if ~is_bgp_ss
            %     bequest_generated_adj = Bequest_generated_raw_total;
            % end
            bq_new = 0;
            if mass_total > 1e-12, bq_new = bequest_generated_adj / mass_total; end
            error_bq = bq_guess - bq_new;

            error_K = K_supply_agg_t_total - (K_supply_agg_tplus1_raw_total / (1 + n_period));
            if ~is_bgp_ss
                error_K = 0; 
                % error_bq = 0;
            end

            F_error = [error_r; error_w; error_tr; error_payg; error_bq; error_K]; %

            if nargout > 1
                ss = struct();
                ss.Y_from_production_hat = prices_new.Y_hat_t;
                ss.K_private_hat = K_supply_agg_t_total;
                ss.K_public_hat = K_g_hat_calc;
                ss.L_hat = L_supply_agg_total;
                ss.r_mkt = r_new;
                ss.w_hat = w_new;
                ss.C_agg = sum(cellfun(@(x) x.C_agg, ss_h));
                ss.I_g = G_investment_spending;
                ss.Total_Tax_Revenue = Total_Tax_Revenue_total;
                ss.Bequest_generated_agg = bequest_generated_adj; 
                ss.Bequest_distributed_agg = bequest_generated_adj;
                ss.Bequest_gen_hat_raw_ss = Bequest_generated_raw_total;
                ss.TR_distributed_agg = tr_new * mass_total;
                ss.Public_Capital_Return = Public_Capital_Return_calc;
                ss.G_c = G_consumption_spending;
                ss.adj = adj_guess; 
                ss.b_hat_formula_h = b_hat_formula_h;
                ss.b_hat_h = adj_guess * b_hat_formula_h; 
                ss.ss_by_type = ss_h;
            end
        end        
        
        
        function [macro_state, Dist, polS, valS] = run_micro_to_macro(r_mkt, w_hat, bq, tr, b_hat_guess, Z_ss_norm, cS, paramS, params_ext)
            % =========================================================================
            % == 函数: run_micro_to_macro
            % == 版本: [v1.9 - 变量名修正版]
            % ==
            % == 核心修改:
            % ==   - 修正了调用 household.VFI_solver 时使用了未定义的变量 M_vfi 的错误。
            % ==   - 正确的变量名应为 M_for_hh。
            % =========================================================================

            b_hat_for_vfi = b_hat_guess;
            M_for_hh.r_mkt_t = r_mkt;
            M_for_hh.w_t = w_hat;
            M_for_hh.tr_per_hh = tr;
            M_for_hh.b_hat_t = b_hat_for_vfi;
            M_for_hh.theta_t = cS.theta_path(end); % cS中应已包含类型h的theta

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
        
        function b_hat_formula_h = calculate_formula_benefits(w_hat_for_formula, cS)
            % =========================================================================
            % ==             养老金福利公式计算器 v1.1 (接口标准化)
            % ==
            % == 核心修改:
            % ==   - 将输入参数名改为 w_hat_for_formula，以明确其作用。
            % ==     在稳态下，它就是当期工资；在转轨中，它将是滞后工资。
            % =========================================================================
            nH = cS.nTypes;
            b_hat_formula_h = zeros(nH, 1);

            eff_at_retirement_age = cS.ageEffV_new_h(cS.aR_new, :); 

            % [核心修改] 使用传入的特定工资来计算福利基数
            social_avg_wage = w_hat_for_formula;

            for h = 1:nH
                individual_indexed_wage = w_hat_for_formula * eff_at_retirement_age(h);

                if h <= 2 
                    b_basic = (social_avg_wage + individual_indexed_wage) / 2 * cS.n_contrib_years * 0.01;
                    b_personal = individual_indexed_wage * cS.personal_account_factor_urban;
                    b_hat_formula_h(h) = b_basic + b_personal;

                else 
                    b_basic = social_avg_wage * cS.resident_basic_benefit_ratio;
                    b_personal = individual_indexed_wage * cS.resident_personal_factor;
                    b_hat_formula_h(h) = b_basic + b_personal;
                end
            end
        end

        function report_struct = display_national_accounts(ss, cS, paramS, Z_ss_norm, report_filename, print_flag, is_bgp_ss)
            % =========================================================================
            % == 函数: display_national_accounts
            % == 版本: [v9.6 - PAYG报告加总逻辑修正版]
            % ==
            % == 核心修正:
            % ==   - [!!! 关键BUG修复 !!!] 修正了[V. 分类型PAYG体系分解报告]中
            % ==     “总量”列的计算错误。
            % ==   - 旧逻辑: 直接对各类型的人均福利向量(ss.b_hat_h)求和，这在
            % ==     经济上是无意义的。
            % ==   - 新逻辑:
            % ==     1. 在主循环中，正确地聚合了“公式化福利”和“最终福利”的总量。
            % ==     2. 在调用 print_row 时，传入了正确计算出的聚合总量，
            % ==        确保了报告的内部一致性。
            % =========================================================================

            if nargin < 6, print_flag = true; end
            if nargin < 5, report_filename = 'ss_national_accounts_report.txt'; end
            
            nH = cS.num_hh_types;

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
                fileID = -1; 
            end

            g_A_period = (1+cS.g_A_ss)^cS.time_Step-1;
            n_period = (1+cS.n_ss)^cS.time_Step-1;
            g_total_period = (1+g_A_period)*(1+n_period)-1;
            Y = ss.Y_from_production_hat;
            mass_total = sum(Z_ss_norm(:));

            % --- 1. 重新计算各类流量 (核心循环) ---
            Dist_by_type = SS.get_final_dist(ss, cS, paramS, Z_ss_norm, is_bgp_ss);
            polS_by_type = SS.get_final_polS(ss, cS, paramS, Z_ss_norm);

            % 初始化分类型和总量聚合变量
            PAYG_contributions_agg_h = zeros(nH, 1);
            PAYG_benefits_agg_h = zeros(nH, 1);
            PAYG_formula_benefits_agg_h = zeros(nH, 1); % [!!! 核心修正 !!!] 新增
            General_Tax_Revenue_total = 0;
            C_agg_total = 0;
            L_hat_total = 0;
            K_p_hat_tplus1_raw_total = 0; 
            Bequest_generated_raw_total = 0;
            
            pps_active = isfield(cS, 'pps_active') && cS.pps_active;

            for h = 1:nH
                Dist_h = Dist_by_type{h};
                polS_h = polS_by_type{h};
                cS_h = cS; 
                cS_h.theta_path = cS.theta_path_h(h, end);
                cS_h.ageEffV_new = cS.ageEffV_new_h(:,h);

                for ia = 1:cS.aD_new
                    mass_dist_ia_h = Dist_h(:,:,:,ia);
                    if sum(mass_dist_ia_h, 'all') < 1e-30, continue; end
                    
                    pol_ia_h = polS_h(ia);
                    
                    C_agg_total = C_agg_total + sum(pol_ia_h.c .* mass_dist_ia_h, 'all');
                    
                    if ia <= cS.aR_new
                        labor_supply_slice = cS_h.ageEffV_new(ia) .* paramS.leGridV'; 
                        labor_supply_grid = repmat(reshape(labor_supply_slice, [1,1,cS.nw_expanded]), [cS.nk, cS.nkpps, 1]);
                        L_hat_total = L_hat_total + sum(labor_supply_grid .* mass_dist_ia_h, 'all');
                    end

                    labor_income_slice_full = ss.w_hat * cS_h.ageEffV_new(ia) .* paramS.leGridV(1:cS.nw_expanded)';
                    payg_contrib_slice_full = cS_h.theta_path * labor_income_slice_full;
                    payg_contrib_grid = repmat(reshape(payg_contrib_slice_full, [1,1,cS.nw_expanded]), [cS.nk, cS.nkpps, 1]);

                    if ia <= cS.aR_new
                        PAYG_contributions_agg_h(h) = PAYG_contributions_agg_h(h) + sum(payg_contrib_grid .* mass_dist_ia_h, 'all');
                    else
                        % [!!! 核心修正 !!!] 同时聚合最终福利和公式化福利
                        PAYG_benefits_agg_h(h) = PAYG_benefits_agg_h(h) + sum(ss.b_hat_h(h) .* mass_dist_ia_h, 'all');
                        PAYG_formula_benefits_agg_h(h) = PAYG_formula_benefits_agg_h(h) + sum(ss.b_hat_formula_h(h) .* mass_dist_ia_h, 'all');
                    end

                    tax_regular_slice = pol_ia_h.tax_regular;
                    general_tax_slice = tax_regular_slice - payg_contrib_grid;
                    General_Tax_Revenue_total = General_Tax_Revenue_total + sum(general_tax_slice .* mass_dist_ia_h, 'all');
                    
                    K_p_hat_tplus1_raw_total = K_p_hat_tplus1_raw_total + sum(pol_ia_h.k_prime .* mass_dist_ia_h, 'all');
                    if pps_active
                         k_prime_for_bequest = pol_ia_h.k_prime + pol_ia_h.kpps_prime;
                         K_p_hat_tplus1_raw_total = K_p_hat_tplus1_raw_total + sum(pol_ia_h.kpps_prime .* mass_dist_ia_h, 'all');
                    else
                         k_prime_for_bequest = pol_ia_h.k_prime;
                    end
                    prob_death = (1 - cS_h.s_pathV(ia));
                    Bequest_generated_raw_total = Bequest_generated_raw_total + sum(k_prime_for_bequest .* mass_dist_ia_h, 'all') * prob_death;
                end
            end
            
            PAYG_contributions_agg_total = sum(PAYG_contributions_agg_h);
            PAYG_benefits_agg_total = sum(PAYG_benefits_agg_h);
            PAYG_formula_benefits_agg_total = sum(PAYG_formula_benefits_agg_h); % [!!! 核心修正 !!!]
            
            % --- 2. 准备国民账户组件 (总量) ---
            I_p_demand_bgp = (g_total_period + cS.ddk) * ss.K_private_hat;
            I_g_display = cS.I_g_to_Y_ratio_ss * Y;
            G_c_display = cS.G_c_to_Y_ratio_ss * Y;
            I_p_accounting = Y - C_agg_total - I_g_display - G_c_display;
            Total_Expenditure = C_agg_total + I_p_accounting + I_g_display + G_c_display;
            resource_constraint_error = Y - Total_Expenditure;
            Depreciation_p = cS.ddk * ss.K_private_hat;
            Depreciation_g = cS.ddk_g * ss.K_public_hat;
            Total_Labor_Income = ss.w_hat * ss.L_hat;
            Total_Capital_Income_Gross = (ss.r_mkt + cS.ddk) * ss.K_private_hat;
            Public_Capital_Return_display = Y - Total_Labor_Income - Total_Capital_Income_Gross;
            Net_Capital_Income = ss.r_mkt * ss.K_private_hat;
            Total_Household_Income = Total_Labor_Income + Net_Capital_Income + PAYG_benefits_agg_total + ss.Bequest_distributed_agg + ss.TR_distributed_agg;
            Household_Savings = Total_Household_Income - C_agg_total - General_Tax_Revenue_total - PAYG_contributions_agg_total;
            
            Total_General_Gov_Revenue = General_Tax_Revenue_total + Public_Capital_Return_display;
            Total_General_Gov_Expenditure = G_c_display + I_g_display + ss.TR_distributed_agg;
            
            % --- 3. 创建并填充 report_struct ---
            report_struct = struct();
            report_struct.Y = Y;
            report_struct.K_p = ss.K_private_hat;
            report_struct.K_g = ss.K_public_hat;
            report_struct.L = ss.L_hat;
            report_struct.r_period = ss.r_mkt;
            report_struct.r_annual = ((1+ss.r_mkt)^(1/cS.time_Step) - 1)*100;
            report_struct.w = ss.w_hat;
            report_struct.C = C_agg_total;
            report_struct.I_p_acct = I_p_accounting;
            report_struct.I_p_bgp = I_p_demand_bgp;
            report_struct.I_g = I_g_display;
            report_struct.G_c = G_c_display;
            report_struct.Depreciation_p = Depreciation_p;
            report_struct.Depreciation_g = Depreciation_g;
            report_struct.Kp_Y_ratio = ss.K_private_hat / Y;
            report_struct.C_Y_ratio = (C_agg_total / Y)*100;
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
            report_struct.Tax_General = General_Tax_Revenue_total;
            report_struct.PAYG_contrib = PAYG_contributions_agg_total;
            report_struct.PAYG_benefit = PAYG_benefits_agg_total;


            % --- 报告开始 ---
            if print_flag
                fprintfAndLog(fileID, '\n\n=========================================================================================================\n');
                fprintfAndLog(fileID, '###            稳 态 国 民 经 济 核 算 报 告 (v9.6 - PAYG报告加总逻辑修正版)           ###\n');
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
                fprintfAndLog(fileID, '   C. 投资与折旧 (流量):\n');
                fprintfAndLog(fileID, '      私人总投资 (I_p, BGP理论值) ............... : %15.6f\n', report_struct.I_p_bgp);
                fprintfAndLog(fileID, '      私人资本折旧 (d*K_p) ...................... : %15.6f\n', report_struct.Depreciation_p);
                fprintfAndLog(fileID, '      私人净投资 (I_p_bgp - d*K_p) .............. : %15.6f\n', report_struct.I_p_bgp - report_struct.Depreciation_p);

                fprintfAndLog(fileID, '\n\n--- [ III. 家庭部门 ] ---\n');
                fprintfAndLog(fileID, '   A. 收入来源 (流量):\n');
                fprintfAndLog(fileID, '      劳动收入 (来自厂商) ....................... : %15.6f\n', report_struct.Total_Labor_Income);
                fprintfAndLog(fileID, '      净资本收入 (r*K_p) ........................ : %15.6f\n', report_struct.Net_Capital_Income);
                fprintfAndLog(fileID, '      养老金福利 ................................ : %15.6f\n', report_struct.PAYG_benefit);
                fprintfAndLog(fileID, '      收到的遗赠 ................................ : %15.6f\n', ss.Bequest_distributed_agg);
                fprintfAndLog(fileID, '      收到的转移支付 ............................ : %15.6f\n', report_struct.TR);
                fprintfAndLog(fileID, '      --------------------------------------------------------------------------------------------------\n');
                fprintfAndLog(fileID, '      >>> 家庭总收入 ............................ : %15.6f\n', Total_Household_Income);
                fprintfAndLog(fileID, '   B. 收入用途 (流量):\n');
                fprintfAndLog(fileID, '      消费支出 (C) .............................. : %15.6f\n', report_struct.C);
                fprintfAndLog(fileID, '      支付的税收 (资本税/劳动税/消费税) ......... : %15.6f\n', report_struct.Tax_General);
                fprintfAndLog(fileID, '      养老金缴费 ................................ : %15.6f\n', report_struct.PAYG_contrib);
                fprintfAndLog(fileID, '      净储蓄 .................................... : %15.6f\n', report_struct.Household_Savings);
                fprintfAndLog(fileID, '      --------------------------------------------------------------------------------------------------\n');
                fprintfAndLog(fileID, '      >>> 家庭总支出与储蓄 ...................... : %15.6f\n', Total_Household_Income);
                fprintfAndLog(fileID, '   C. 财富与遗赠 (存量与流量):\n');
                fprintfAndLog(fileID, '      持有总资产 (K_p, 下一期初) ................ : %15.6f\n', report_struct.K_p);
                fprintfAndLog(fileID, '      产生的总遗赠 (留给下一代) ................. : %15.6f\n', ss.Bequest_generated_agg);

                fprintfAndLog(fileID, '\n\n--- [ IV. 政府与养老金体系 ] ---\n');
                fprintfAndLog(fileID, '   A. 一般政府预算:\n');
                fprintfAndLog(fileID, '      (+) 一般税收收入 .......................... : %15.6f\n', report_struct.Tax_General);
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
                fprintfAndLog(fileID, '      >>> 预算平衡误差 (收入 - 支出) ............ : %15.4e\n', Total_General_Gov_Revenue - Total_General_Gov_Expenditure);
                
                fprintfAndLog(fileID, '\n   B. 现收现付(PAYG)养老金体系 (总量):\n');
                fprintfAndLog(fileID, '      均衡调整因子 (adj) ........................ : %15.6f\n', ss.adj);
                fprintfAndLog(fileID, '      (+) 养老金总缴费 .......................... : %15.6f\n', report_struct.PAYG_contrib);
                fprintfAndLog(fileID, '      (-) 养老金总福利 .......................... : %15.6f\n', report_struct.PAYG_benefit);
                fprintfAndLog(fileID, '      >>> PAYG体系总余额 (缴费 - 福利) .......... : %15.4e\n', PAYG_contributions_agg_total - PAYG_benefits_agg_total);
                                
                fprintfAndLog(fileID, '\n\n--- [ V. 分类型PAYG体系分解报告 ] ---\n');
                header_format = '%-25s |';
                cell_format_val = ' %12.4f |';
                cell_format_pct = ' (%5.1f%%) |';
                for h=1:nH, header_format = [header_format, '     类型 %-15d |']; end
                header_format = [header_format, ' %12s |\n'];
                fprintfAndLog(fileID, header_format, '变量', 1:nH, '总量');
                fprintfAndLog(fileID, '%s\n', repmat('-', 1, 30 + nH*25 + 15));
                
                print_row = @(label, values_h, total_val) ...
                    fprintfAndLog(fileID, ['%-25s |' repmat([cell_format_val, cell_format_pct], 1, nH) cell_format_val '\n'], ...
                    label, compose_print_array(values_h, total_val), total_val);
                
                % [!!! 核心修正 !!!] 使用正确聚合的总量
                print_row('PAYG 缴费', PAYG_contributions_agg_h, report_struct.PAYG_contrib);
                print_row('PAYG 公式化福利', PAYG_formula_benefits_agg_h, PAYG_formula_benefits_agg_total);
                print_row('PAYG 最终福利(含adj)', PAYG_benefits_agg_h, report_struct.PAYG_benefit);

                fprintfAndLog(fileID, '\n\n--- [ VI. 核心均衡条件检验 (来自最终均衡解) ] ---\n');
                
                final_prices = firm.get_prices_at_t(ss.K_private_hat, ss.K_public_hat, ss.L_hat, cS);
                error_r_final = ss.r_mkt - final_prices.r_mkt_t;
                error_w_final = ss.w_hat - final_prices.w_hat_t;
                
                general_tax_final = ss.Total_Tax_Revenue - report_struct.PAYG_contrib;
                gov_revenue_final = general_tax_final + ss.Public_Capital_Return;
                gov_expenditure_final = ss.G_c + ss.I_g + ss.TR_distributed_agg;
                error_tr_final = gov_revenue_final - gov_expenditure_final;
                
                error_payg_final = PAYG_contributions_agg_total - PAYG_benefits_agg_total;

                bequest_adj_final = Bequest_generated_raw_total / (1 + n_period);
                error_bq_final = bequest_adj_final - ss.Bequest_distributed_agg;

                error_K_final = 0;
                if is_bgp_ss
                    error_K_final = ss.K_private_hat - (K_p_hat_tplus1_raw_total / (1 + n_period));
                end
                
                fprintfAndLog(fileID, '   1. 资本市场出清 (r_supply - r_demand) ..... : %15.4e\n', error_r_final);
                fprintfAndLog(fileID, '   2. 劳动市场出清 (w_supply - w_demand) ..... : %15.4e\n', error_w_final);
                fprintfAndLog(fileID, '   3. 政府预算平衡 (Gov_Rev - Gov_Exp) ..... : %15.4e\n', error_tr_final);
                fprintfAndLog(fileID, '   4. PAYG体系平衡 (Total Contrib - Total Ben) : %15.4e\n', error_payg_final);
                fprintfAndLog(fileID, '   5. 遗赠市场出清 (Bq_gen - Bq_dist) ....... : %15.4e\n', error_bq_final);
                fprintfAndLog(fileID, '   6. BGP资本积累 (K_supply - K_demand) ..... : %15.4e\n', error_K_final);


                fprintfAndLog(fileID, '\n=========================================================================================================\n');
                fprintfAndLog(fileID, '###                                      报 告 结 束                                      ###\n');
                fprintfAndLog(fileID, '=========================================================================================================\n\n');

                if fileID > 1
                    fclose(fileID);
                    fprintf('\n报告已成功保存至: %s\n', report_filename);
                end
            end 

            function fprintfAndLog(fileID, formatSpec, varargin)
                if print_flag
                    fprintf(1, formatSpec, varargin{:});
                    if fileID > 1, fprintf(fileID, formatSpec, varargin{:}); end
                end
            end
            
            function print_array = compose_print_array(values_h, total_val)
                print_array = [];
                for val_i = values_h(:)'
                    pct = 0;
                    if abs(total_val) > 1e-9, pct = (val_i / total_val) * 100; end
                    print_array = [print_array, val_i, pct];
                end
            end
        end
        
        function Dist_by_type = get_final_dist(ss, cS, paramS, Z_ss_norm, is_bgp_ss)
            % =========================================================================
            % == 函数: get_final_dist
            % == 版本: [v2.2 - 统一PAYG调整因子版]
            % ==
            % == 核心修改:
            % ==   - 构造x_eq_final时，用ss.adj替换了ss.b_hat_h向量。
            % =========================================================================
            mass_total = sum(Z_ss_norm(:));
            nH = cS.num_hh_types;
            
            tr_per_capita = ss.TR_distributed_agg / mass_total;
            bq_per_capita = ss.Bequest_distributed_agg / mass_total;

            % [核心修改] 使用 adj 构造最终的x向量
            x_eq_final = [ss.r_mkt; ss.w_hat; tr_per_capita; ss.adj; bq_per_capita];
            
            [~, ~, Dist_by_type] = SS.system_of_equations(x_eq_final, Z_ss_norm, cS, paramS, struct('A',1.0), is_bgp_ss);
        end
        
        function polS_by_type = get_final_polS(ss, cS, paramS, Z_ss_norm)
            % =========================================================================
            % == 函数: get_final_polS
            % == 版本: [v2.2 - 统一PAYG调整因子版]
            % ==
            % == 核心修改:
            % ==   - 传递给VFI的 b_hat_t 现在是最终的、经过adj调整后的福利。
            % =========================================================================
            nH = cS.num_hh_types;
            polS_by_type = cell(nH, 1);
            mass_total = sum(Z_ss_norm(:));

            for h = 1:nH
                cS_h = cS; 
                cS_h.theta_path = cS.theta_path_h(h);
                cS_h.ageEffV_new = cS.ageEffV_new_h(:,h);
                
                M_for_hh.r_mkt_t = ss.r_mkt;
                M_for_hh.w_t = ss.w_hat;
                M_for_hh.tr_per_hh = ss.TR_distributed_agg / mass_total;
                M_for_hh.b_hat_t = ss.b_hat_h(h); % [核心修改] 使用最终福利
                M_for_hh.theta_t = cS_h.theta_path(end);

                [polS_by_type{h}, ~] = household.VFI_solver(M_for_hh, paramS, cS_h);
            end
        end

    end
    
end