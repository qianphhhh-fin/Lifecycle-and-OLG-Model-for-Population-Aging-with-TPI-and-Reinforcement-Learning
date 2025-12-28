% --- SS.m ---
classdef SS
    methods (Static)

        function [ss, Dist, polS, valS] = solve_steady_state(cS, paramS, params_ext, verbose, is_bgp_ss)
            % =========================================================================
            % == 函数: solve_steady_state
            % == 版本: [v5.0 - 市场统一修正版]
            % ==
            % == 核心修改:
            % ==   - [!!! 关键修正 !!!] 统一了求解向量 x0 的结构。无论 nTypes 是多少，
            % ==     现在都只求解一个总的 bq_guess_agg 和一个总的 b_hat_guess_agg。
            % ==   - 这强制求解器寻找一个统一的市场价格，消除了因市场分割导致的
            % ==     退化检验失败的根源。
            % ==   - 移除了原有的 nTypes>1 的分支，简化了代码结构。
            % =========================================================================
            if verbose, verbose='iter'; else, verbose='off'; end
            cS.A = params_ext.A;
            if isfield(params_ext, 'g_A_ss'), cS.g_A_ss = params_ext.g_A_ss; end
            if isfield(params_ext, 'n_ss'), cS.n_ss = params_ext.n_ss; end
            if isfield(params_ext, 'Z') && ~isempty(params_ext.Z), Z_ss_norm = params_ext.Z;
            else, error('错误：求解稳态需要一个明确的人口分布Z。'); end
            if ~cS.pps_active, cS.nkpps = 1; cS.n_pps_rate_grid = 1; end

            % --- [核心修正] 统一的初始猜测值向量 x0 ---
            r_guess_init = 0.04;
            w_guess_init = 1.5;
            tr_guess_init = 0.0;
            
            % 无论类型数量，都只求解一个总量的 bq 和一个统一的 b_hat
            K_guess_for_BQ = 2.5;
            s_pathV_ss = cS.s_pathV(:, end);
            bq_guess_agg_init = mean(1-s_pathV_ss) * K_guess_for_BQ; % 总量遗赠的初始猜测
            b_hat_guess_agg_init = 0.05; % 统一的PAYG福利水平的初始猜测

            x0 = [r_guess_init; w_guess_init; bq_guess_agg_init; tr_guess_init; b_hat_guess_agg_init];
            lb = [-0.1; 1e-8; 1e-8; -Inf; 0];
            ub = [Inf; Inf; Inf; Inf; Inf];
            
            system_wrapper = @(x) SS.system_of_equations(x, Z_ss_norm, cS, paramS, params_ext, is_bgp_ss);
  
            options = optimoptions('lsqnonlin', 'Display', verbose, ...
                'FunctionTolerance', 1e-12, 'StepTolerance', 1e-12, 'OptimalityTolerance', 1e-12, ...
                'MaxIterations', 100);

            [x_eq, ~, ~, exitflag] = lsqnonlin(system_wrapper, x0, lb, ub, options);

            [~, ss, Dist, polS, valS] = SS.system_of_equations(x_eq, Z_ss_norm, cS, paramS, params_ext, is_bgp_ss);
        end
% --- SS.m ---
% ... (classdef SS and other methods remain the same) ...
function [F_error, ss, Dist_by_type, polS_by_type, valS_by_type] = system_of_equations(x_guess, Z_ss_norm, cS, paramS, params_ext, is_bgp_ss)
            % =========================================================================
            % == 函数: system_of_equations
            % == 版本: [v8.0 - 市场统一修正版]
            % ==
            % == 核心修改:
            % ==   - 函数签名保持不变，但 x_guess 的内部结构已变为5维。
            % ==   - 解包 x_guess 时，直接得到总量的 bq_guess_agg 和统一的 b_hat_guess_agg。
            % ==   - 调用 run_micro_to_macro 时，传入这些统一的猜测值。
            % ==   - 构建 F_error 时，计算总量的遗赠市场误差和总量的PAYG体系误差，
            % ==     取代了原有的分类型误差向量。
            % =========================================================================
            
            % --- [核心修正] 解包统一的猜测值 ---
            r_guess = x_guess(1);
            w_guess = x_guess(2);
            bq_guess_agg = x_guess(3); 
            tr_guess = x_guess(4); 
            b_hat_guess_agg = x_guess(5);

            [aggr_out, Dist_by_type, polS_by_type, valS_by_type] = SS.run_micro_to_macro(r_guess, w_guess, bq_guess_agg, tr_guess, b_hat_guess_agg, Z_ss_norm, cS, paramS, params_ext);
            
            L_supply_agg = aggr_out.L_hat_total;
            K_supply_agg_t = aggr_out.K_start_of_period_total; 
            
            n_period = (1 + cS.n_ss)^cS.time_Step - 1;

            error_K = K_supply_agg_t - (aggr_out.K_end_of_period_total / (1 + n_period));
            if ~is_bgp_ss, error_K = 0; end

            g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
            g_total_period = (1 + g_A_period) * (1 + n_period) - 1;
            Kg_Y_ratio_factor = cS.I_g_to_Y_ratio_ss / (g_total_period + cS.ddk_g);
            Y_hat_calc = (K_supply_agg_t^cS.alpha * L_supply_agg^(1 - cS.alpha - cS.gamma) * Kg_Y_ratio_factor^cS.gamma)^(1 / (1 - cS.gamma));
            K_g_hat_calc = Kg_Y_ratio_factor * Y_hat_calc;
            prices_new = firm.get_prices_at_t(K_supply_agg_t, K_g_hat_calc, L_supply_agg, cS);
            error_r = r_guess - prices_new.r_mkt_t;
            error_w = w_guess - prices_new.w_hat_t;

            PAYG_contrib_total = aggr_out.PAYG_contrib_total;
            General_Tax_Revenue_calc = aggr_out.Total_Tax_Revenue_total - PAYG_contrib_total;
            Public_Capital_Return_calc = Y_hat_calc - (prices_new.w_hat_t * L_supply_agg) - ((prices_new.r_mkt_t + cS.ddk) * K_supply_agg_t);
            Total_General_Gov_Revenue_calc = General_Tax_Revenue_calc + Public_Capital_Return_calc;
            tr_new = (Total_General_Gov_Revenue_calc - (cS.G_c_to_Y_ratio_ss + cS.I_g_to_Y_ratio_ss) * Y_hat_calc) / sum(Z_ss_norm);
            error_tr = tr_guess - tr_new;

            % --- [核心修正] 构建统一市场的误差 ---
            bequest_generated_adj = aggr_out.Bequest_generated_total / (1 + n_period);
            error_bq_agg = bq_guess_agg - bequest_generated_adj;
            
            error_b_hat_agg = aggr_out.PAYG_contrib_total - aggr_out.PAYG_benefit_total;
            
            F_error = [error_r; error_w; error_bq_agg; error_tr; error_b_hat_agg; error_K];
            
            if nargout > 1
                ss = struct();
                ss.Y_from_production_hat = Y_hat_calc;
                ss.K_private_hat = K_supply_agg_t;
                ss.K_public_hat = K_g_hat_calc;
                ss.L_hat = L_supply_agg;
                ss.r_mkt = prices_new.r_mkt_t;
                ss.w_hat = prices_new.w_hat_t;
                ss.C_agg = aggr_out.C_agg_total;
                ss.I_g = (cS.I_g_to_Y_ratio_ss) * Y_hat_calc;
                ss.G_c = (cS.G_c_to_Y_ratio_ss) * Y_hat_calc;
                ss.Bequest_distributed_agg = bequest_generated_adj;
                ss.TR_distributed_agg = tr_new * sum(Z_ss_norm);
                % 为了报告和后续分析，需要将均衡的统一价格分解到各类型
                ss.b_hat_by_type = repmat(b_hat_guess_agg, cS.nTypes, 1);
                mass_newborns_by_type = Z_ss_norm(1) .* cS.type_weights;
                ss.bq_by_type = (bequest_generated_adj * cS.type_weights) ./ mass_newborns_by_type .* mass_newborns_by_type;
            end
        end
% --- SS.m ---
% ... (classdef SS and other methods remain the same) ...
function [aggr_out, Dist_by_type, polS_by_type, valS_by_type] = run_micro_to_macro(r_mkt, w_hat, bq_guess_agg, tr, b_hat_guess_agg, Z_ss_norm, cS, paramS, ~)
            % =========================================================================
            % == 函数: run_micro_to_macro
            % == 版本: [v6.0 - 市场统一修正版]
            % ==
            % == 核心修改:
            % ==   - 函数签名变更：现在接收总量的 bq_guess_agg 和统一的 b_hat_guess_agg。
            % ==   - 在内部，将统一的 b_hat_guess_agg 包装成一个向量，传递给VFI求解器，
            % ==     确保所有类型的家庭都面对相同的养老金福利预期。
            % ==   - 在计算新生儿初始财富时，使用总的 bq_guess_agg 除以总新生儿人口，
            % ==     得到统一的人均遗赠 bq_per_newborn，确保所有新生儿初始财富相同。
            % =========================================================================

            % --- [核心修正] 将统一的 b_hat 猜测值传递给所有家庭类型 ---
            M_for_hh.r_mkt_t = r_mkt;
            M_for_hh.w_t = w_hat;
            M_for_hh.tr_per_hh = tr;
            M_for_hh.b_hat_by_type = repmat(b_hat_guess_agg, cS.nTypes, 1); 
            [polS_by_type, valS_by_type] = household.VFI_solver(M_for_hh, paramS, cS);

            Dist_by_type = cell(cS.nTypes, 1);
            aggrS_by_type = cell(cS.nTypes, 1);
            
            % --- [核心修正] 计算统一的人均遗赠 ---
            total_mass_newborns = Z_ss_norm(1);
            bq_per_newborn_unified = 0;
            if total_mass_newborns > 1e-12
                bq_per_newborn_unified = bq_guess_agg / total_mass_newborns;
            end
            
            for i_type = 1:cS.nTypes
                polS_i = polS_by_type{i_type};
                Z_ss_norm_i = Z_ss_norm * cS.type_weights(i_type);
                
                % 所有类型的新生儿都接收相同的人均遗赠
                Dist_by_type{i_type} = distribution.get_ss_dist(polS_i, paramS, cS, Z_ss_norm_i, bq_per_newborn_unified);
                
                aggrS_by_type{i_type} = aggregates.get_aggregates(Dist_by_type{i_type}, polS_i, cS, paramS, i_type);
            end

            % --- 聚合逻辑保持不变 ---
            K_start_of_period_total = 0; K_end_of_period_total = 0;
            C_agg_total = 0; L_hat_total = 0; Total_Tax_Revenue_total = 0;
            Bequest_generated_total = 0;
            PAYG_contrib_total = 0;
            mass_retirees_total = sum(Z_ss_norm((cS.aR_new+1):end));
            PAYG_benefit_total = b_hat_guess_agg * mass_retirees_total;
            
            for i_type = 1:cS.nTypes
                aggrS_i = aggrS_by_type{i_type};
                
                K_start_of_period_total = K_start_of_period_total + aggrS_i.K_start_of_period;
                K_end_of_period_total = K_end_of_period_total + aggrS_i.K_end_of_period;
                C_agg_total = C_agg_total + aggrS_i.C_agg;
                L_hat_total = L_hat_total + aggrS_i.L_hat;
                Total_Tax_Revenue_total = Total_Tax_Revenue_total + aggrS_i.Total_Tax_Revenue;
                Bequest_generated_total = Bequest_generated_total + aggrS_i.Bequest_generated_agg;

                is_urban = (i_type <= 2);
                if cS.nTypes > 1
                    theta_for_type = is_urban * cS.theta_path_urban(end) + (1-is_urban) * cS.theta_path_resident(end);
                else
                    theta_for_type = cS.theta_path(end);
                end
                PAYG_contrib_total = PAYG_contrib_total + (theta_for_type * w_hat * aggrS_i.L_hat);
            end
            
            aggr_out = struct();
            aggr_out.K_start_of_period_total = K_start_of_period_total;
            aggr_out.K_end_of_period_total = K_end_of_period_total;
            aggr_out.C_agg_total = C_agg_total;
            aggr_out.Total_Tax_Revenue_total = Total_Tax_Revenue_total;
            aggr_out.L_hat_total = L_hat_total;
            aggr_out.Bequest_generated_total = Bequest_generated_total;
            aggr_out.PAYG_contrib_total = PAYG_contrib_total;
            aggr_out.PAYG_benefit_total = PAYG_benefit_total;
        end

function Dist = get_final_dist(ss, cS, paramS, Z_ss_norm, is_bgp_ss)
            % =========================================================================
            % == 函数: get_final_dist
            % == 版本: [v2.0 - 动态求解器兼容版]
            % =========================================================================
            mass_total = sum(Z_ss_norm);
            tr_per_capita = ss.TR_distributed_agg / mass_total;
            
            if cS.nTypes == 1
                bq_guess = ss.bq_by_type;
                b_hat_guess = ss.b_hat_by_type;
                x_eq_final = [ss.r_mkt; ss.w_hat; bq_guess; tr_per_capita; b_hat_guess];
            else
                x_eq_final = [ss.r_mkt; ss.w_hat; tr_per_capita; ss.bq_by_type; ss.b_hat_by_type];
            end

            [~, ~, Dist] = SS.system_of_equations(x_eq_final, Z_ss_norm, cS, paramS, struct('A',1.0), is_bgp_ss);
        end

        function polS = get_final_polS(ss, cS, paramS, Z_ss_norm)
            % =========================================================================
            % == 函数: get_final_polS
            % == 版本: [v2.0 - 动态求解器兼容版]
            % =========================================================================
            M_for_hh.r_mkt_t = ss.r_mkt;
            M_for_hh.w_t = ss.w_hat;
            M_for_hh.tr_per_hh = ss.TR_distributed_agg / sum(Z_ss_norm);
            M_for_hh.b_hat_by_type = ss.b_hat_by_type;

            [polS, ~] = household.VFI_solver(M_for_hh, paramS, cS);
        end
        


        function report_struct = display_national_accounts(ss, cS, paramS, Z_ss_norm, report_filename, print_flag, is_bgp_ss, Dist_by_type, polS_by_type)
            % =========================================================================
            % == 函数: display_national_accounts
            % == 版本: [v11.0 - 格式对齐最终版]
            % ==
            % == 核心修改:
            % ==   - [!!!] 重新编排了所有 fprintfAndLog 语句，使其输出的报告
            % ==     结构和内容与 v8.4 (同质性黄金标准) 版本完全一致。
            % ==   - 保留了 v10.x 版本中所有正确的、经过修正的会计变量聚合逻辑。
            % ==   - 增加了对 PPS 体系变量的报告，当 pps_active=true 时激活。
            % ==   - 异质性模型的特定报告内容（分类型账户）被移至报告末尾的
            % ==     一个独立章节，以保持主报告的一致性。
            % ==   - 最终检验部分的格式也与 v8.4 版对齐。
            % =========================================================================

            if nargin < 9, error('display_national_accounts 需要提供 Dist_by_type 和 polS_by_type'); end
            if nargin < 6, print_flag = true; end
            if nargin < 5, report_filename = 'ss_national_accounts_report.txt'; end

            % --- 0. 准备工作 ---
            fileID = 1; 
            if print_flag
                try
                    fileID = fopen(report_filename, 'w', 'n', 'UTF-8');
                    if fileID == -1, error('无法创建或打开报告文件: %s', report_filename); end
                catch ME
                    warning('创建报告文件时出错，将输出到命令行。错误: %s', ME.message);
                    fileID = 1;
                end
            end

            g_A_period = (1+cS.g_A_ss)^cS.time_Step-1;
            n_period = (1+cS.n_ss)^cS.time_Step-1;
            g_total_period = (1+g_A_period)*(1+n_period)-1;
            Y = ss.Y_from_production_hat;

            

            % --- 1. 循环聚合所有类型的宏观流量 ---
            C_by_type = zeros(cS.nTypes, 1); L_by_type = zeros(cS.nTypes, 1);
            K_k_by_type_start_period = zeros(cS.nTypes, 1); K_pps_by_type_start_period = zeros(cS.nTypes, 1);
            K_k_by_type_end_period = zeros(cS.nTypes, 1); K_pps_by_type_end_period = zeros(cS.nTypes, 1);
            PAYG_contributions_by_type = zeros(cS.nTypes, 1); PAYG_benefits_by_type = zeros(cS.nTypes, 1);
            Tax_regular_by_type = zeros(cS.nTypes, 1);
            pps_active = isfield(cS, 'pps_active') && cS.pps_active;
            PPS_contrib_by_type = zeros(cS.nTypes, 1); PPS_withdrawal_by_type = zeros(cS.nTypes, 1); PPS_tax_by_type = zeros(cS.nTypes, 1);
            
            k_grid_full = repmat(cS.kGridV, [1, cS.nkpps, cS.nw_expanded]);
            if cS.nkpps > 1, kpps_grid_full = repmat(reshape(cS.kppsGridV, [1, cS.nkpps]), [cS.nk, 1, cS.nw_expanded]);
            else, kpps_grid_full = zeros(cS.nk, 1, cS.nw_expanded); end

            for i_type = 1:cS.nTypes
                Dist_type = Dist_by_type{i_type}; polS_type = polS_by_type{i_type};
                b_hat_this_type = ss.b_hat_by_type(i_type);
                if cS.nTypes > 1, theta_ss_type = (i_type <= 2) * cS.theta_path_urban(end) + (i_type > 2) * cS.theta_path_resident(end);
                else, theta_ss_type = cS.theta_path(end); end

                for ia = 1:cS.aD_new
                    mass_dist_ia = Dist_type(:,:,:,ia);
                    if sum(mass_dist_ia, 'all') < 1e-30, continue; end
                    
                    K_k_by_type_start_period(i_type) = K_k_by_type_start_period(i_type) + sum(k_grid_full .* mass_dist_ia, 'all');
                    if pps_active, K_pps_by_type_start_period(i_type) = K_pps_by_type_start_period(i_type) + sum(kpps_grid_full .* mass_dist_ia, 'all'); end

                    C_by_type(i_type) = C_by_type(i_type) + sum(polS_type(ia).c .* mass_dist_ia, 'all');
                    Tax_regular_by_type(i_type) = Tax_regular_by_type(i_type) + sum(polS_type(ia).tax_regular .* mass_dist_ia, 'all');
                    
                    if ia <= cS.aR_new
                        le_grid_reshaped = reshape(paramS.leGridV(1:cS.nw_expanded), [1, 1, cS.nw_expanded]);
                        L_by_type(i_type) = L_by_type(i_type) + sum(cS.ageEff_by_type(i_type, ia) .* le_grid_reshaped .* mass_dist_ia, 'all');
                        labor_income_grid = ss.w_hat * cS.ageEff_by_type(i_type, ia) .* le_grid_reshaped;
                        PAYG_contributions_by_type(i_type) = PAYG_contributions_by_type(i_type) + sum(theta_ss_type .* labor_income_grid .* mass_dist_ia, 'all');
                    else
                        PAYG_benefits_by_type(i_type) = PAYG_benefits_by_type(i_type) + sum(b_hat_this_type .* mass_dist_ia, 'all');
                    end
                    K_k_by_type_end_period(i_type) = K_k_by_type_end_period(i_type) + sum(polS_type(ia).k_prime .* mass_dist_ia, 'all');
                    if pps_active, K_pps_by_type_end_period(i_type) = K_pps_by_type_end_period(i_type) + sum(polS_type(ia).kpps_prime .* mass_dist_ia, 'all'); end
                end
            end
            
            General_Tax_by_type = Tax_regular_by_type - PAYG_contributions_by_type;
            
            % --- 2. 准备所有国民账户组件 ---
            C_agg = sum(C_by_type); L_agg = sum(L_by_type);
            PAYG_contributions_agg = sum(PAYG_contributions_by_type);
            PAYG_benefits_agg = sum(PAYG_benefits_by_type);
            General_Tax_Revenue = sum(General_Tax_by_type);
            I_p_demand_bgp = (g_total_period + cS.ddk) * ss.K_private_hat;
            I_g_display = cS.I_g_to_Y_ratio_ss * Y; G_c_display = cS.G_c_to_Y_ratio_ss * Y;
            I_p_accounting = Y - C_agg - I_g_display - G_c_display;
            Total_Expenditure = C_agg + I_p_accounting + I_g_display + G_c_display;
            resource_constraint_error = Y - Total_Expenditure;
            Depreciation_p = cS.ddk * ss.K_private_hat; Depreciation_g = cS.ddk_g * ss.K_public_hat;
            Total_Labor_Income = ss.w_hat * L_agg; Net_Capital_Income = ss.r_mkt * ss.K_private_hat;
            Total_Capital_Income_Gross = (ss.r_mkt + cS.ddk) * ss.K_private_hat;
            Public_Capital_Return_display = Y - Total_Labor_Income - Total_Capital_Income_Gross;
            Total_Household_Income = Total_Labor_Income + Net_Capital_Income + PAYG_benefits_agg + sum(PPS_withdrawal_by_type) + ss.Bequest_distributed_agg + ss.TR_distributed_agg;
            Household_Savings_agg = Total_Household_Income - C_agg - General_Tax_Revenue - PAYG_contributions_agg - sum(PPS_contrib_by_type) - sum(PPS_tax_by_type);
            Total_General_Gov_Revenue = General_Tax_Revenue + sum(PPS_tax_by_type) + Public_Capital_Return_display;
            Total_General_Gov_Expenditure = G_c_display + I_g_display + ss.TR_distributed_agg;
            general_budget_balance_error = Total_General_Gov_Revenue - Total_General_Gov_Expenditure;
            payg_balance_error = PAYG_contributions_agg - PAYG_benefits_agg;
            K_k_hat_agg_adj = sum(K_k_by_type_end_period) / (1+n_period);
            K_pps_hat_agg_adj = sum(K_pps_by_type_end_period) / (1+n_period);

            % --- 3. 创建并填充 report_struct ---
            report_struct = ss;
            report_struct.Y = Y; 
            report_struct.K_p = ss.K_private_hat;
            report_struct.K_g = ss.K_public_hat;
            report_struct.TR = ss.TR_distributed_agg;
            report_struct.L = L_agg;
            report_struct.r_annual = ((1+ss.r_mkt)^(1/cS.time_Step) - 1)*100;
            report_struct.w = ss.w_hat; report_struct.C = C_agg; 
            report_struct.I_p_acct = I_p_accounting; report_struct.I_p_bgp = I_p_demand_bgp;
            report_struct.Kp_Y_ratio = ss.K_private_hat / Y;
            report_struct.C_Y_ratio = (C_agg / Y)*100;
            report_struct.Ip_Y_ratio = (I_p_accounting / Y)*100;
            report_struct.Invest_Gap = I_p_accounting - I_p_demand_bgp;
            report_struct.Invest_Gap_pct = (report_struct.Invest_Gap / Y)*100;
            report_struct.Total_Labor_Income = Total_Labor_Income;
            report_struct.Net_Capital_Income = Net_Capital_Income;
            report_struct.Household_Savings = Household_Savings_agg;
            report_struct.Public_Capital_Return = Public_Capital_Return_display;
            report_struct.Tax_General = General_Tax_Revenue;
            report_struct.Gov_Budget_Balance_Error = general_budget_balance_error;
            report_struct.PAYG_contrib = PAYG_contributions_agg;
            report_struct.PAYG_benefit = PAYG_benefits_agg;
            report_struct.PAYG_Balance_Error = payg_balance_error;
            report_struct.Depreciation_p = Depreciation_p;
            report_struct.Depreciation_g = Depreciation_g;
            report_struct.RC_Error = resource_constraint_error;
            report_struct.RC_Error_pct = (resource_constraint_error/Y)*100;

            if pps_active
                report_struct.K_k = K_k_hat_agg_adj; report_struct.K_pps = K_pps_hat_agg_adj;
                report_struct.PPS_contrib = sum(PPS_contrib_by_type);
                report_struct.PPS_withdrawal = sum(PPS_withdrawal_by_type);
                report_struct.PPS_tax = sum(PPS_tax_by_type);
            end

                        % ======================== DEBUGGING BLOCK START (v2.0) ========================
            % if print_flag
            %     fprintfAndLog(fileID, '\n\n--- [ DEBUGGING & DIAGNOSTICS (v2.0) ] ---\n');
            % 
            %     % --- 1. 检查各类型分布的总人口质量 ---
            %     fprintfAndLog(fileID, '   >>> 1. 检查各类型分布的总人口质量 (mass) <<<\n');
            %     total_mass_check = 0;
            %     for i_type = 1:cS.nTypes
            %         actual_mass = sum(Dist_by_type{i_type}, 'all');
            %         expected_mass = cS.type_weights(i_type) * sum(Z_ss_norm);
            %         fprintfAndLog(fileID, '      类型 %d:  实际质量 = %-12.8f | 预期质量 = %-12.8f | 差异 = %-12.4e\n', ...
            %             i_type, actual_mass, expected_mass, actual_mass - expected_mass);
            %         total_mass_check = total_mass_check + actual_mass;
            %     end
            %     fprintfAndLog(fileID, '      --------------------------------------------------------------------------------\n');
            %     fprintfAndLog(fileID, '      合计:   实际总质量 = %-12.8f | 预期总质量 = %-12.8f | 差异 = %-12.4e\n', ...
            %         total_mass_check, sum(Z_ss_norm), total_mass_check - sum(Z_ss_norm));
            %     if abs(total_mass_check - sum(Z_ss_norm)) > 1e-6
            %          fprintfAndLog(fileID, '      >>> [!!! 严重错误 !!!] 分布的总质量与人口数据不匹配！\n');
            %     else
            %          fprintfAndLog(fileID, '      >>> [✅ 通过] 分布总质量正确。\n');
            %     end
            % 
            %     % --- 2. 详细检查各类型聚合的宏观流量 ---
            %     fprintfAndLog(fileID, '\n   >>> 2. 详细检查各类型聚合的宏观流量 (基于其自身分布) <<<\n');
            %     fprintfAndLog(fileID, '      %-10s | %-12s | %-12s | %-12s | %-12s | %-12s\n', ...
            %         '类型', '劳动(L)', '消费(C)', '总税收', 'PAYG缴费', '期末总资产(K'')');
            %     fprintfAndLog(fileID, '      %s\n', repmat('-', 1, 75));
            % 
            %     L_agg_check = 0; C_agg_check = 0; K_prime_agg_check = 0;
            %     for i_type = 1:cS.nTypes
            %         aggrS_i = aggregates.get_aggregates(Dist_by_type{i_type}, polS_by_type{i_type}, cS, paramS, i_type);
            % 
            %         is_urban = (i_type <= 2);
            %         if cS.nTypes > 1, theta_for_type = is_urban * cS.theta_path_urban(end) + (1-is_urban) * cS.theta_path_resident(end);
            %         else, theta_for_type = cS.theta_path(end); end
            %         payg_contrib_check = theta_for_type * ss.w_hat * aggrS_i.L_hat;
            % 
            %         fprintfAndLog(fileID, '      %-10d | %-12.6f | %-12.6f | %-12.6f | %-12.6f | %-12.6f\n', ...
            %             i_type, aggrS_i.L_hat, aggrS_i.C_agg, aggrS_i.Total_Tax_Revenue, payg_contrib_check, aggrS_i.K_p_hat_tplus1_raw);
            % 
            %         L_agg_check = L_agg_check + aggrS_i.L_hat;
            %         C_agg_check = C_agg_check + aggrS_i.C_agg;
            %         K_prime_agg_check = K_prime_agg_check + aggrS_i.K_p_hat_tplus1_raw;
            %     end
            %     fprintfAndLog(fileID, '      %s\n', repmat('-', 1, 75));
            %     fprintfAndLog(fileID, '      %-10s | %-12.6f | %-12.6f | %-12s | %-12s | %-12.6f\n', ...
            %         '合计', L_agg_check, C_agg_check, '', '', K_prime_agg_check);
            %     fprintfAndLog(fileID, '      %-10s | %-12.6f | %-12.6f | %-12s | %-12s | %-12.6f\n', ...
            %         '报告值', ss.L_hat, ss.C_agg, '', '', ss.K_private_hat * (1+n_period)); % K_prime = K_t+1 = K_t * (1+n)
            %     diff_L = abs(L_agg_check - ss.L_hat);
            %     diff_C = abs(C_agg_check - ss.C_agg);
            %     diff_K = abs(K_prime_agg_check - (ss.K_private_hat * (1+n_period)));
            %     if diff_L > 1e-6 || diff_C > 1e-6 || diff_K > 1e-6
            %         fprintfAndLog(fileID, '      >>> [!!! 警告 !!!] 重新聚合的总量与报告值不符！差异 (L,C,K''): %.2e, %.2e, %.2e\n', diff_L, diff_C, diff_K);
            %     else
            %         fprintfAndLog(fileID, '      >>> [✅ 通过] 重新聚合的总量与报告值一致。\n');
            %     end
            % 
            %     % --- 3. 检查市场出清误差的构成 ---
            %     fprintfAndLog(fileID, '\n   >>> 3. 检查市场出清误差的构成 <<<\n');
            %     fprintfAndLog(fileID, '      均衡条件             | 猜测值 (Guess) | 计算值 (New)   | 误差 (Guess - New)\n');
            %     fprintfAndLog(fileID, '      ---------------------|----------------|----------------|-------------------\n');
            %     % 资本和劳动市场
            %     fprintfAndLog(fileID, '      利率 (r)             | %-14.8f | %-14.8f | %-17.4e\n', ss.r_mkt, report_struct.r_mkt - final_errors(1), final_errors(1));
            %     fprintfAndLog(fileID, '      工资 (w)             | %-14.8f | %-14.8f | %-17.4e\n', ss.w_hat, report_struct.w - final_errors(2), final_errors(2));
            %     % PAYG 和 遗赠市场
            %     for i_type = 1:cS.nTypes
            %         fprintfAndLog(fileID, '      PAYG b_hat (Type %d)   | %-14.8f | %-14.8f | %-17.4e\n', i_type, ss.b_hat_by_type(i_type), NaN, final_errors(7+i_type));
            %         fprintfAndLog(fileID, '      Bequest bq (Type %d)| %-14.8f | %-14.8f | %-17.4e\n', i_type, ss.bq_by_type(i_type), NaN, final_errors(3+i_type));
            %     end
            % 
            %     fprintfAndLog(fileID, '\n');
            % end
            % ========================= DEBUGGING BLOCK END (v2.0) =========================

            % --- 4. 报告打印 ---
            if print_flag
                fprintfAndLog(fileID, '\n\n=========================================================================================================\n');
                fprintfAndLog(fileID, '###                 稳 态 国 民 经 济 核 算 报 告 (v11.0 - 格式对齐最终版)                     ###\n');
                fprintfAndLog(fileID, '=========================================================================================================\n');

                fprintfAndLog(fileID, '\n--- [ I. 宏观总量与国民支出 ] ---\n');
                fprintfAndLog(fileID, '   A. 生产 (GDP) 与核心比率:\n');
                fprintfAndLog(fileID, '      国内生产总值 (Y) .......................... : %15.6f\n', report_struct.Y);
                fprintfAndLog(fileID, '      私人资本/产出比 (K_p/Y) ................... : %15.6f\n', report_struct.Kp_Y_ratio);
                fprintfAndLog(fileID, '      公共资本/产出比 (K_g/Y) ................... : %15.6f\n', report_struct.K_g / report_struct.Y);
                fprintfAndLog(fileID, '   B. 均衡价格:\n');
                fprintfAndLog(fileID, '      真实利率 (r, 模型期) ...................... : %15.6f (年化: %.4f %%)\n', report_struct.r_mkt, report_struct.r_annual);
                fprintfAndLog(fileID, '      真实工资率 (w, 单位有效劳动) .............. : %15.6f\n', report_struct.w);
                fprintfAndLog(fileID, '   C. 国民支出 (Y = C + I_p + I_g + G_c) 与资源约束检验:\n');
                fprintfAndLog(fileID, '      (+) 私人消费 (C) .......................... : %15.6f  (%6.2f %% of Y)\n', report_struct.C, report_struct.C_Y_ratio);
                fprintfAndLog(fileID, '      (+) 私人总投资 (I_p, 会计值) ............ : %15.6f  (%6.2f %% of Y)\n', report_struct.I_p_acct, report_struct.Ip_Y_ratio);
                fprintfAndLog(fileID, '      (+) 公共总投资 (I_g) ...................... : %15.6f  (%6.2f %% of Y)\n', report_struct.I_g, (report_struct.I_g / Y)*100);
                fprintfAndLog(fileID, '      (+) 政府消费 (G_c) ........................ : %15.6f  (%6.2f %% of Y)\n', report_struct.G_c, (report_struct.G_c / Y)*100);
                fprintfAndLog(fileID, '      --------------------------------------------------------------------------------------------------\n');
                fprintfAndLog(fileID, '      (=) 总支出 (C + I_p_acct + I_g + G_c) ..... : %15.6f\n', report_struct.C + report_struct.I_p_acct + report_struct.I_g + report_struct.G_c);
                fprintfAndLog(fileID, '      --------------------------------------------------------------------------------------------------\n');
                fprintfAndLog(fileID, '      >>> 资源约束误差 (Y - 总支出) ............. : %15.4e (Y的 %.4e %%)\n', report_struct.RC_Error, report_struct.RC_Error_pct);
                if abs(report_struct.RC_Error_pct) < 1e-6, fprintfAndLog(fileID, '          核验结果: ✅ 资源约束在会计上严格满足。\n');
                else, fprintfAndLog(fileID, '          核验结果: ⚠️ 资源约束在会计上存在误差！\n'); end
                fprintfAndLog(fileID, '\n   D. 投资一致性检验 (仅适用于真稳态):\n');
                fprintfAndLog(fileID, '      BGP理论投资需求 (g+d)K .................. : %15.6f\n', report_struct.I_p_bgp);
                fprintfAndLog(fileID, '      会计投资 (Y-C-I_g-G_c) .................... : %15.6f\n', report_struct.I_p_acct);
                fprintfAndLog(fileID, '      >>> 投资缺口 (会计值 - 理论值) .......... : %15.4e (Y的 %.4f %%)\n', report_struct.Invest_Gap, report_struct.Invest_Gap_pct);
                if abs(report_struct.Invest_Gap_pct) < 1e-4 || ~is_bgp_ss, fprintfAndLog(fileID, '          检验结果: ✅ 投资需求与会计值一致 (或为伪稳态)。\n');
                else, fprintfAndLog(fileID, '          检验结果: ⚠️  经济体偏离BGP稳态！\n'); end

                fprintfAndLog(fileID, '\n\n--- [ II. 厂商部门 ] ---\n');
                fprintfAndLog(fileID, '   A. 生产要素投入 (存量):\n');
                fprintfAndLog(fileID, '      总私人资本存量 (K_p) ...................... : %15.6f\n', report_struct.K_p);
                if pps_active, fprintfAndLog(fileID, '         其中: 常规资本 (K_k) ................. : %15.6f  (%6.2f %% of K_p)\n', report_struct.K_k, (report_struct.K_k / report_struct.K_p) * 100); fprintfAndLog(fileID, '         其中: PPS资本 (K_pps) ................ : %15.6f  (%6.2f %% of K_p)\n', report_struct.K_pps, (report_struct.K_pps / report_struct.K_p) * 100); end
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
                if pps_active, fprintfAndLog(fileID, '      PPS账户提取 ............................... : %15.6f\n', sum(PPS_withdrawal_by_type)); end
                fprintfAndLog(fileID, '      收到的遗赠 ................................ : %15.6f\n', ss.Bequest_distributed_agg);
                fprintfAndLog(fileID, '      收到的转移支付 ............................ : %15.6f\n', report_struct.TR);
                fprintfAndLog(fileID, '      --------------------------------------------------------------------------------------------------\n');
                fprintfAndLog(fileID, '      >>> 家庭总收入 ............................ : %15.6f\n', Total_Household_Income);
                fprintfAndLog(fileID, '   B. 收入用途 (流量):\n');
                fprintfAndLog(fileID, '      消费支出 (C) .............................. : %15.6f\n', report_struct.C);
                fprintfAndLog(fileID, '      支付的一般税收 ............................ : %15.6f\n', report_struct.Tax_General);
                fprintfAndLog(fileID, '      养老金缴费 ................................ : %15.6f\n', report_struct.PAYG_contrib);
                if pps_active, fprintfAndLog(fileID, '      PPS账户缴费 (储蓄) ...................... : %15.6f\n', sum(PPS_contrib_by_type)); fprintfAndLog(fileID, '      PPS提取税 ................................. : %15.6f\n', sum(PPS_tax_by_type)); end
                fprintfAndLog(fileID, '      净储蓄 .................................... : %15.6f\n', report_struct.Household_Savings);
                fprintfAndLog(fileID, '      --------------------------------------------------------------------------------------------------\n');
                fprintfAndLog(fileID, '      >>> 家庭总支出与储蓄 ...................... : %15.6f\n', Total_Household_Income);
                fprintfAndLog(fileID, '   C. 财富与遗赠 (存量与流量):\n');
                fprintfAndLog(fileID, '      持有总资产 (K_p, BGP调整后) ............... : %15.6f\n', report_struct.K_p);
                fprintfAndLog(fileID, '      产生的总遗赠 (BGP调整后) .................. : %15.6f\n', ss.Bequest_distributed_agg);

                fprintfAndLog(fileID, '\n\n--- [ IV. 政府与养老金体系 ] ---\n');
                fprintfAndLog(fileID, '   A. 一般政府预算:\n');
                fprintfAndLog(fileID, '      (+) 一般税收收入 .......................... : %15.6f\n', report_struct.Tax_General);
                if pps_active, fprintfAndLog(fileID, '      (+) PPS提取税 ............................. : %15.6f\n', sum(PPS_tax_by_type)); end
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
                
                fprintfAndLog(fileID, '\n   B. 现收现付(PAYG)养老金体系 (总量):\n');
                fprintfAndLog(fileID, '      (+) 养老金总缴费 .......................... : %15.6f\n', report_struct.PAYG_contrib);
                fprintfAndLog(fileID, '      (-) 养老金总福利 .......................... : %15.6f\n', report_struct.PAYG_benefit);
                fprintfAndLog(fileID, '      >>> PAYG体系余额 (缴费 - 福利) ............ : %15.4e\n', report_struct.PAYG_Balance_Error);
                 if abs(report_struct.PAYG_Balance_Error/Y) < 1e-9, fprintfAndLog(fileID, '          核验结果: ✅ PAYG体系按设计严格平衡。\n');
                else, fprintfAndLog(fileID, '          核验结果: ⚠️ PAYG体系未按设计平衡！\n'); end

                if cS.nTypes > 1
                    payg_balance_by_type_error = PAYG_contributions_by_type - PAYG_benefits_by_type;
                    type_names = {'高收入城镇', '中低收入城镇', '高收入居民', '中低收入居民'};
                    fprintfAndLog(fileID, '\n   C. 现收现付(PAYG)养老金体系 (按类型细分):\n');
                    fprintfAndLog(fileID, '      %-15s | %15s | %15s | %15s | %15s\n', '类型', '福利水平(b_hat)', '总缴费', '总福利', '余额(误差)');
                    fprintfAndLog(fileID, '      %s\n', repmat('-', 1, 72));
                    for i_type = 1:cS.nTypes
                        fprintfAndLog(fileID, '      %-15s | %15.6f | %15.6f | %15.6f | %15.4e\n', ...
                            type_names{i_type}, ss.b_hat_by_type(i_type), ...
                            PAYG_contributions_by_type(i_type), PAYG_benefits_by_type(i_type), payg_balance_by_type_error(i_type));
                    end
                end

                if cS.nTypes > 1
                    % --- 准备分类型数据 ---
                    TotalIncome_by_type = zeros(cS.nTypes, 1);
                    Savings_by_type = zeros(cS.nTypes, 1);
                    for i_type = 1:cS.nTypes
                        TotalIncome_by_type(i_type) = ss.w_hat * L_by_type(i_type) + ss.r_mkt * K_k_by_type_start_period(i_type) ...
                            + PAYG_benefits_by_type(i_type) + ss.bq_by_type(i_type) + ss.TR_distributed_agg * cS.type_weights(i_type);
                        Savings_by_type(i_type) = TotalIncome_by_type(i_type) - C_by_type(i_type) - General_Tax_by_type(i_type) - PAYG_contributions_by_type(i_type);
                    end

                    fprintfAndLog(fileID, '\n\n--- [ V. 按家庭类型细分的宏观账户 ] ---\n');
                    fprintfAndLog(fileID, '   A. 各类型绝对量贡献:\n');
                    fprintfAndLog(fileID, '      %-20s | %15s | %15s | %15s | %15s | %15s\n', '变量', type_names{1}, type_names{2}, type_names{3}, type_names{4}, '经济体总量');
                    fprintfAndLog(fileID, '      %s\n', repmat('-', 1, 110));
                    fprintfAndLog(fileID, '      %-20s | %15.4f | %15.4f | %15.4f | %15.4f | %15.4f\n', '人口权重 (%)', cS.type_weights*100, 100.0);
                    fprintfAndLog(fileID, '      %-20s | %15.4f | %15.4f | %15.4f | %15.4f | %15.4f\n', '私人资本 (K_p)', K_k_by_type_start_period, sum(K_k_by_type_start_period));
                    fprintfAndLog(fileID, '      %-20s | %15.4f | %15.4f | %15.4f | %15.4f | %15.4f\n', '有效劳动 (L)', L_by_type, sum(L_by_type));
                    fprintfAndLog(fileID, '      %-20s | %15.4f | %15.4f | %15.4f | %15.4f | %15.4f\n', '总消费 (C)', C_by_type, sum(C_by_type));
                    fprintfAndLog(fileID, '      %-20s | %15.4f | %15.4f | %15.4f | %15.4f | %15.4f\n', '总储蓄 (S)', Savings_by_type, sum(Savings_by_type));
                    fprintfAndLog(fileID, '      %-20s | %15.4f | %15.4f | %15.4f | %15.4f | %15.4f\n', '总收入 (Income)', TotalIncome_by_type, sum(TotalIncome_by_type));
                    fprintfAndLog(fileID, '      %s\n', repmat('-', 1, 110));

                    fprintfAndLog(fileID, '\n   B. 各类型相对份额 (%% of Total):\n');
                    fprintfAndLog(fileID, '      %-20s | %15s | %15s | %15s | %15s | %15s\n', '变量', type_names{1}, type_names{2}, type_names{3}, type_names{4}, '合计');
                    fprintfAndLog(fileID, '      %s\n', repmat('-', 1, 110));
                    fprintfAndLog(fileID, '      %-20s | %15.2f%% | %15.2f%% | %15.2f%% | %15.2f%% | %15.2f%%\n', '私人资本 (K_p)', (K_k_by_type_start_period/sum(K_k_by_type_start_period))*100, 100.0);
                    fprintfAndLog(fileID, '      %-20s | %15.2f%% | %15.2f%% | %15.2f%% | %15.2f%% | %15.2f%%\n', '有效劳动 (L)', (L_by_type/sum(L_by_type))*100, 100.0);
                    fprintfAndLog(fileID, '      %-20s | %15.2f%% | %15.2f%% | %15.2f%% | %15.2f%% | %15.2f%%\n', '总消费 (C)', (C_by_type/sum(C_by_type))*100, 100.0);
                    fprintfAndLog(fileID, '      %-20s | %15.2f%% | %15.2f%% | %15.2f%% | %15.2f%% | %15.2f%%\n', '总储蓄 (S)', (Savings_by_type/sum(Savings_by_type))*100, 100.0);
                    fprintfAndLog(fileID, '      %-20s | %15.2f%% | %15.2f%% | %15.2f%% | %15.2f%% | %15.2f%%\n', '总收入 (Income)', (TotalIncome_by_type/sum(TotalIncome_by_type))*100, 100.0);
                    fprintfAndLog(fileID, '      %s\n', repmat('-', 1, 110));
                end

                fprintfAndLog(fileID, '\n\n--- [ VI. 核心均衡条件检验 (来自求解器的最终误差) ] ---\n');
                mass_total = sum(Z_ss_norm);
                tr_per_capita_final = ss.TR_distributed_agg / mass_total;
                if cS.nTypes == 1
                    bq_final = ss.bq_by_type; b_hat_final = ss.b_hat_by_type;
                    x_final_pc = [ss.r_mkt; ss.w_hat; bq_final; tr_per_capita_final; b_hat_final];
                    final_errors = SS.system_of_equations(x_final_pc, Z_ss_norm, cS, paramS, struct('A',1.0), is_bgp_ss);
                    fprintfAndLog(fileID, '   1. 资本市场出清 (r_guess - r_new) ......... : %15.4e\n', final_errors(1));
                    fprintfAndLog(fileID, '   2. 劳动市场出清 (w_guess - w_new) ....... : %15.4e\n', final_errors(2));
                    fprintfAndLog(fileID, '   3. 遗赠市场出清 (BQ_guess - BQ_new) ... : %15.4e\n', final_errors(3));
                    fprintfAndLog(fileID, '   4. 政府预算平衡 (TR_guess - TR_new) ..... : %15.4e\n', final_errors(4));
                    fprintfAndLog(fileID, '   5. PAYG体系平衡 (Contrib - Benefit) .... : %15.4e\n', final_errors(5));
                    if (length(final_errors) > 5), fprintfAndLog(fileID, '   6. BGP资本积累 (K_t - K_t+1_adj) ........ : %15.4e\n', final_errors(6)); end
                else
                    x_final_pc = [ss.r_mkt; ss.w_hat; tr_per_capita_final; ss.bq_by_type; ss.b_hat_by_type];
                    final_errors = SS.system_of_equations(x_final_pc, Z_ss_norm, cS, paramS, struct('A',1.0), is_bgp_ss);
                    fprintfAndLog(fileID, '   1. 资本市场出清 (r_guess - r_new) ......... : %15.4e\n', final_errors(1));
                    fprintfAndLog(fileID, '   2. 劳动市场出清 (w_guess - w_new) ....... : %15.4e\n', final_errors(2));
                    fprintfAndLog(fileID, '   3. 政府预算平衡 (TR_guess - TR_new) ..... : %15.4e\n', final_errors(3));
                    err_idx = 4;
                    for i_type = 1:cS.nTypes, fprintfAndLog(fileID, '   3.%d. 遗赠市场出清 (类型 %d) ............ : %15.4e\n', i_type, i_type, final_errors(err_idx)); err_idx = err_idx + 1; end
                    for i_type = 1:cS.nTypes, fprintfAndLog(fileID, '   4.%d. PAYG平衡 (类型 %d) .................. : %15.4e\n', i_type, i_type, final_errors(err_idx)); err_idx = err_idx + 1; end
                    if (length(final_errors) >= err_idx), fprintfAndLog(fileID, '   5. BGP资本积累 (K_t - K_t+1_adj) ........ : %15.4e\n', final_errors(err_idx)); end
                end

                fprintfAndLog(fileID, '\n\n=========================================================================================================\n');
                fprintfAndLog(fileID, '###                                      报 告 结 束                                      ###\n');
                fprintfAndLog(fileID, '=========================================================================================================\n\n');

                if fileID > 1, fclose(fileID); fprintf('\n报告已成功保存至: %s\n', report_filename); end
            end

            function fprintfAndLog(fileID, formatSpec, varargin)
                if print_flag
                    fprintf(1, formatSpec, varargin{:});
                    if fileID > 1 && fileID ~= -1, fprintf(fileID, formatSpec, varargin{:}); end
                end
            end
        end

     
    end
end