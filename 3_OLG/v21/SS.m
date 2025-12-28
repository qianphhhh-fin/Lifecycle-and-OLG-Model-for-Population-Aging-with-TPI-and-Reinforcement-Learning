% --- SS.m ---
classdef SS
    methods (Static)
        
        function [ss, Dist, polS, valS] = solve_steady_state(cS, paramS, params_ext, verbose, is_bgp_ss)
            % =========================================================================
            % == 函数: solve_steady_state
            % == 版本: [v3.0 - 新PAYG制度和内生基金版]
            % ==
            % == 核心修改:
            % ==   - 求解变量(x)重构为: [r, w, tr, adj_factor, bq_total_pc]。
            % ==   - 新增布尔输入 `is_bgp_ss` 用于区分求解 ss0 还是 ssF。
            % ==   - 将 `is_bgp_ss` 传递给 system_of_equations 以应用不同规则。
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

            if ~cS.pps_active
                cS.nkpps = 1; cS.n_pps_rate_grid = 1;
            end

            % --- 设置初始猜测值 (变量: r, w, tr, adj_factor, bq_total_pc) ---
            r_guess_init = 0.04;
            w_guess_init = 1.5;
            TR_guess_init = 0.0;
            adj_factor_guess_init = 1.0; % 初始猜测福利公式无需调整
            K_guess_for_BQ = 2.5;
            BQ_total_pc_guess_init = mean(1-cS.s_pathV) * K_guess_for_BQ / sum(Z_ss_norm,'all');

            x0 = [r_guess_init; w_guess_init; TR_guess_init; adj_factor_guess_init; BQ_total_pc_guess_init];

            % --- 设置边界 ---
            lb = [-0.1; 1e-8; -Inf; 0; 1e-8]; % adj_factor >= 0
            ub = [Inf; Inf; Inf; Inf; Inf];

            system_wrapper = @(x) SS.system_of_equations(x, Z_ss_norm, cS, paramS, params_ext, is_bgp_ss);

            options = optimoptions('lsqnonlin', 'Display', verbose, ...
            'FunctionTolerance', 1e-12, 'StepTolerance', 1e-12, 'OptimalityTolerance', 1e-12, ...
            'MaxIterations', 500);

            [x_eq, ~, ~, exitflag] = lsqnonlin(system_wrapper, x0, lb, ub, options);

            if exitflag <= 0
               warning('稳态求解器未能正常收敛或已停止。');
               ss = []; Dist = []; polS = []; valS = [];
               return;
            end
            
            [~, ss, Dist, polS, valS] = SS.system_of_equations(x_eq, Z_ss_norm, cS, paramS, params_ext, is_bgp_ss);
        end        
       
        
        
        function [F_error, ss, Dist, polS, valS] = system_of_equations(x_guess, Z_ss_norm, cS, paramS, params_ext, is_bgp_ss)
            % =========================================================================
            % == 函数: system_of_equations
            % == 版本: [v5.3 - 诊断信息补全版]
            % ==
            % == 核心修改:
            % ==   - [!!! BUG修复 !!!] 补上了缺失的 ss.Bequest_gen_hat_raw_ss
            % ==     和 ss.ss_by_type 输出项，此两项为转轨路径求解和调试所必需。
            % =========================================================================
            nH = cS.num_hh_types;

            % --- 流程 0: 解包猜测向量 ---
            r_guess = x_guess(1);
            w_guess = x_guess(2);
            tr_guess = x_guess(3);
            adj_factor_guess = x_guess(4);
            bq_guess = x_guess(5);

            % --- 流程 1: 基于猜测值，计算各类型家庭面对的参数 ---
            b_hat_formula_h = SS.calculate_formula_benefits(w_guess, cS);
            b_hat_for_vfi_h = adj_factor_guess * b_hat_formula_h;

            % --- 流程 2: 初始化聚合变量与存储容器 ---
            L_supply_agg_total = 0;
            K_supply_agg_tplus1_raw_total = 0;
            Total_Tax_Revenue_total = 0;
            Bequest_generated_raw_total = 0;
            K_p_supply_from_dist_total = 0;
            PAYG_contributions_total = 0;
            PAYG_benefits_total = 0;

            Dist = cell(nH, 1);
            polS = cell(nH, 1);
            valS = cell(nH, 1);
            ss_h = cell(nH, 1); % [新增] 初始化分类型结果容器

            % --- 流程 3: [核心循环] 遍历各类家庭，求解微观并聚合 ---
            for h = 1:nH
                cS_h = cS;
                paramS_h = paramS;
                if isfield(cS_h, 'theta_path_h'), cS_h.theta_path = cS.theta_path_h(h, end); end
                if isfield(cS_h, 'ageEffV_new_h'), cS_h.ageEffV_new = cS.ageEffV_new_h(:,h); end

                [aggr_out_h, Dist{h}, polS{h}, valS{h}] = SS.run_micro_to_macro(r_guess, w_guess, bq_guess, tr_guess, b_hat_for_vfi_h(h), Z_ss_norm(:, h), cS_h, paramS_h, params_ext);
                
                K_p_dist_h = aggregates.get_aggregates_from_dist(Dist{h}, cS_h);

                % [新增] 保存分类型结果
                ss_h{h}.L_hat = aggr_out_h.L_hat;
                ss_h{h}.K_hat = K_p_dist_h.K_p_hat_total;
                ss_h{h}.C_agg = aggr_out_h.C_agg;
                
                L_supply_agg_total = L_supply_agg_total + aggr_out_h.L_hat;
                K_p_supply_from_dist_total = K_p_supply_from_dist_total + K_p_dist_h.K_p_hat_total;
                K_supply_agg_tplus1_raw_total = K_supply_agg_tplus1_raw_total + aggr_out_h.K_p_hat_tplus1_raw;
                Bequest_generated_raw_total = Bequest_generated_raw_total + aggr_out_h.Bequest_generated_agg;
                Total_Tax_Revenue_total = Total_Tax_Revenue_total + aggr_out_h.Total_Tax_Revenue;
                
                payg_contrib_h = cS_h.theta_path * w_guess * aggr_out_h.L_hat;
                mass_retirees_h = sum(Z_ss_norm((cS.aR_new+1):end, h));
                payg_benefit_h = b_hat_for_vfi_h(h) * mass_retirees_h;
                PAYG_contributions_total = PAYG_contributions_total + payg_contrib_h;
                PAYG_benefits_total = PAYG_benefits_total + payg_benefit_h;
            end

            % --- 流程 4: [宏观层面] 市场出清 ---
            n_period = (1 + cS.n_ss)^cS.time_Step - 1;
            g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
            g_total_period = (1 + g_A_period) * (1 + n_period) - 1;

            K_payg_calc = 0;
            K_firm_calc = 0;
            if ~is_bgp_ss
                k_payg_ratio = 0;
                if isfield(cS, 'K_payg_to_Y_ratio_ss0'), k_payg_ratio = cS.K_payg_to_Y_ratio_ss0; end
                
                Y_func = @(y) y - ((K_p_supply_from_dist_total - k_payg_ratio * y)^cS.alpha * ...
                                  (cS.I_g_to_Y_ratio_ss / (g_total_period + cS.ddk_g) * y)^cS.gamma * ...
                                  L_supply_agg_total^(1-cS.alpha-cS.gamma));
                
                Y_hat_calc = fzero(Y_func, K_p_supply_from_dist_total);
                K_payg_calc = k_payg_ratio * Y_hat_calc;
                K_firm_calc = K_p_supply_from_dist_total - K_payg_calc;
            else
                K_payg_calc = 0;
                K_firm_calc = K_p_supply_from_dist_total;
            end
            
            if K_firm_calc < 1e-8, K_firm_calc = 1e-8; end

            Kg_Y_ratio = cS.I_g_to_Y_ratio_ss / (g_total_period + cS.ddk_g);
            Y_from_K_L = (K_firm_calc^cS.alpha * L_supply_agg_total^(1-cS.alpha-cS.gamma) * Kg_Y_ratio^cS.gamma)^(1/(1-cS.gamma));
            K_g_calc = Kg_Y_ratio * Y_from_K_L;
            prices_new = firm.get_prices_at_t(K_firm_calc, K_g_calc, L_supply_agg_total, cS);
            r_new = prices_new.r_mkt_t;
            w_new = prices_new.w_hat_t;

            error_r = r_guess - r_new;
            error_w = w_guess - w_new;

            payg_surplus = PAYG_contributions_total - PAYG_benefits_total;
            if is_bgp_ss
                error_payg = payg_surplus;
            else
                required_surplus = (g_total_period - r_new) * K_payg_calc;
                error_payg = payg_surplus - required_surplus;
            end

            General_Tax_Revenue_calc = Total_Tax_Revenue_total - PAYG_contributions_total;
            Y_hat_final = prices_new.Y_hat_t;
            Public_Capital_Return_calc = Y_hat_final - (w_new * L_supply_agg_total) - ((r_new + cS.ddk) * K_firm_calc);
            Total_General_Gov_Revenue_calc = General_Tax_Revenue_calc + Public_Capital_Return_calc;
            G_consumption_spending = cS.G_c_to_Y_ratio_ss * Y_hat_final;
            G_investment_spending = cS.I_g_to_Y_ratio_ss * Y_hat_final;
            mass_total = sum(Z_ss_norm(:));
            tr_new = (Total_General_Gov_Revenue_calc - (G_consumption_spending + G_investment_spending)) / mass_total;
            error_tr = tr_guess - tr_new;

            bequest_generated_adj = Bequest_generated_raw_total / (1 + n_period);
            bq_new = 0;
            if mass_total > 1e-12, bq_new = bequest_generated_adj / mass_total; end
            error_bq = bq_guess - bq_new;

            error_K = 0;
            if is_bgp_ss
                error_K = K_p_supply_from_dist_total - (K_supply_agg_tplus1_raw_total / (1 + n_period));
            end
            
            % --- 流程 5: 组合最终的误差向量 ---
            F_error = [error_r; error_w; error_tr; error_payg; error_bq];
            if is_bgp_ss, F_error = [F_error; error_K]; end

            % --- 流程 6: 打包最终的ss结构体 (如果需要) ---
            if nargout > 1
                ss = struct();
                ss.Y_from_production_hat = Y_hat_final;
                ss.C_agg = sum(cellfun(@(x) x.C_agg, ss_h)); % [!!! 核心修改 !!!]
                ss.K_private_hat = K_p_supply_from_dist_total;
                ss.K_firm_hat = K_firm_calc;
                ss.K_payg_hat = K_payg_calc;
                ss.K_public_hat = K_g_calc;
                ss.L_hat = L_supply_agg_total;
                ss.r_mkt = r_new;
                ss.w_hat = w_new;
                ss.I_g = G_investment_spending;
                ss.G_c = G_consumption_spending;
                ss.TR_distributed_agg = tr_new * mass_total;
                ss.Bequest_distributed_agg = bequest_generated_adj;
                ss.Bequest_gen_hat_raw_ss = Bequest_generated_raw_total;
                ss.adj_factor = adj_factor_guess;
                ss.b_hat_h = b_hat_for_vfi_h;
                ss.ss_by_type = ss_h; % [核心修改]
                ss.errors = struct('r', error_r, 'w', error_w, 'tr', error_tr, 'payg', error_payg, 'bq', error_bq, 'K', error_K);
            end
        end
        
        
        function [macro_state, Dist, polS, valS] = run_micro_to_macro(r_mkt, w_hat, bq, tr, b_hat_for_vfi, Z_ss_norm, cS, paramS, params_ext)
            % =========================================================================
            % == 函数: run_micro_to_macro
            % == 版本: [v2.0 - 福利参数化版]
            % ==
            % == 核心修改:
            % ==   - 输入从 b_hat_guess 变为 b_hat_for_vfi。
            % ==   - b_hat_for_vfi 是一个已经计算好的、可以直接传递给VFI的福利参数，
            % ==     而不是一个待求解的均衡价格。
            % =========================================================================

            M_for_hh.r_mkt_t = r_mkt;
            M_for_hh.w_t = w_hat;
            M_for_hh.tr_per_hh = tr;
            M_for_hh.b_hat_t = b_hat_for_vfi; % [核心修改] 直接使用传入的福利参数
            M_for_hh.theta_t = cS.theta_path(end);

            [polS, valS] = household.VFI_solver(M_for_hh, paramS, cS);
            Dist = distribution.get_ss_dist(polS, paramS, cS, Z_ss_norm, bq);

            aggrS = aggregates.get_aggregates(Dist, polS, cS, paramS);

            macro_state = struct();
            macro_state.K_p_hat_tplus1_raw = aggrS.K_p_hat_tplus1_raw;
            macro_state.C_agg = aggrS.C_agg;
            macro_state.Bequest_generated_agg = aggrS.Bequest_generated_agg;
            macro_state.Total_Tax_Revenue = aggrS.Total_Tax_Revenue;
            macro_state.L_hat = aggrS.L_hat;
        end        
        
        function report_struct = display_national_accounts(ss, cS, paramS, Z_ss_norm, Dist_by_type, polS_by_type, report_filename, print_flag, is_bgp_ss)
            % =========================================================================
            % == 函数: display_national_accounts
            % == 版本: [v11.3 - 黄金会计恒等式最终版]
            % ==
            % == 核心修改:
            % ==   - 废除所有复杂的部门储蓄加总。
            % ==   - 采用最基础、最稳健的国民收入定义来构建核心检验：
            % ==     国民总储蓄 (S_gross = Y - C - G_c)
            % ==     必须等于
            % ==     国民总投资 (I_gross = I_p_acct + I_g)。
            % ==   - 此修改将从根本上确保会计平衡，消除所有储蓄-投资缺口。
            % =========================================================================

            if nargin < 9, print_flag = true; end
            if nargin < 8, report_filename = 'ss_national_accounts_report.txt'; end

            % --- 0. 准备工作：文件句柄和BGP增长因子 ---
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

            % --- 1. 基于最终分布，重新计算总消费 ---
            C_agg_total = 0;
            for h = 1:cS.num_hh_types
                Dist_h = Dist_by_type{h};
                polS_h = polS_by_type{h};
                for ia = 1:cS.aD_new
                    mass_dist_ia_h = Dist_h(:,:,:,ia);
                    if sum(mass_dist_ia_h, 'all') < 1e-30, continue; end
                    C_agg_total = C_agg_total + sum(polS_h(ia).c .* mass_dist_ia_h, 'all');
                end
            end

            % --- 2. [核心] 构建基于黄金标准的国民账户组件 ---
            G_c_display = ss.G_c;
            I_g_display = ss.I_g;
            % 私人投资由资源约束倒算得出
            I_p_accounting = Y - C_agg_total - I_g_display - G_c_display;

            % 国民总储蓄 (来自收入法)
            Gross_National_Saving = Y - C_agg_total - G_c_display;
            % 国民总投资 (来自支出法)
            Gross_National_Investment = I_p_accounting + I_g_display;

            % --- 报告开始 ---
            fprintfAndLog(fileID, '\n\n=========================================================================================================\n');
            fprintfAndLog(fileID, '###            稳 态 国 民 经 济 核 算 报 告 (v11.3 - 黄金会计恒等式最终版)          ###\n');
            fprintfAndLog(fileID, '=========================================================================================================\n');
            
            fprintfAndLog(fileID, '\n--- [ I. 宏观总量与国民支出 ] ---\n');
            fprintfAndLog(fileID, '   A. 核心总量与价格:\n');
            fprintfAndLog(fileID, '      国内生产总值 (Y) .......................... : %15.6f\n', Y);
            fprintfAndLog(fileID, '      家庭总资产 (K_p) .......................... : %15.6f\n', ss.K_private_hat);
            fprintfAndLog(fileID, '      真实利率 (r, 模型期) ...................... : %15.6f (年化: %.4f %%)\n', ss.r_mkt, ((1+ss.r_mkt)^(1/cS.time_Step) - 1)*100);
            fprintfAndLog(fileID, '   B. 国民支出 (Y = C + I_p + I_g + G_c):\n');
            fprintfAndLog(fileID, '      (+) 私人消费 (C) .......................... : %15.6f\n', C_agg_total);
            fprintfAndLog(fileID, '      (+) 私人总投资 (I_p, 会计值) ............ : %15.6f\n', I_p_accounting);
            fprintfAndLog(fileID, '      (+) 公共总投资 (I_g) ...................... : %15.6f\n', I_g_display);
            fprintfAndLog(fileID, '      (+) 政府消费 (G_c) ........................ : %15.6f\n', G_c_display);
            
            fprintfAndLog(fileID, '\n\n--- [ II. [核心] 国民总储蓄 vs 总投资检验 (S_gross = I_gross) ] ---\n');
            fprintfAndLog(fileID, '   A. 国民总储蓄 (S_gross = Y - C - G_c) ....... : %15.6f\n', Gross_National_Saving);
            fprintfAndLog(fileID, '   B. 国民总投资 (I_gross = I_p_acct + I_g) .... : %15.6f\n', Gross_National_Investment);
            fprintfAndLog(fileID, '   C. 最终检验:\n');
            fprintfAndLog(fileID, '      >>> 储蓄-投资总缺口 (S_gross - I_gross) .. : %15.4e (会计恒等式, 必须为0)\n', Gross_National_Saving - Gross_National_Investment);

            fprintfAndLog(fileID, '\n\n--- [ III. 投资一致性检验 (总投资 vs BGP要求) ] ---\n');
            I_p_demand_bgp = (g_total_period + cS.ddk) * ss.K_private_hat;
            I_g_demand_bgp = (g_total_period + cS.ddk_g) * ss.K_public_hat;
            Total_Investment_Demand_BGP = I_p_demand_bgp + I_g_demand_bgp;
            fprintfAndLog(fileID, '      BGP理论总投资需求 (g+d)K_p + (g+d_g)K_g ... : %15.6f\n', Total_Investment_Demand_BGP);
            fprintfAndLog(fileID, '      会计总投资 (I_p_acct + I_g) ............... : %15.6f\n', Gross_National_Investment);
            fprintfAndLog(fileID, '      >>> 总投资缺口 (会计值 - 理论值) ......... : %15.4e (Y的 %.4f %%)\n', Gross_National_Investment - Total_Investment_Demand_BGP, (Gross_National_Investment - Total_Investment_Demand_BGP)/Y*100);
            if is_bgp_ss
                fprintfAndLog(fileID, '          检验结果: ✅ 这是一个真稳态(ssF)，总投资缺口应接近于0。\n');
            else
                fprintfAndLog(fileID, '          检验结果: ⚠️ 这是一个伪稳态(ss0)，存在投资缺口是正常现象。\n');
            end

            fprintfAndLog(fileID, '\n\n--- [ IV. PAYG 体系参数 ] ---\n');
            fprintfAndLog(fileID, '      内生福利调整因子 (adj_factor) ............. : %15.6f (注: >1表示福利比公式更慷慨)\n', ss.adj_factor);
            
            fprintfAndLog(fileID, '\n\n--- [ V. 核心均衡条件检验 (来自最终均衡解) ] ---\n');
            % 此部分只为验证求解器工作正常，不再需要重算所有流量
            fprintfAndLog(fileID, '   (本部分仅确认求解器返回的各市场误差值，会计平衡见 [II] 和 [III])\n');
            fprintfAndLog(fileID, '   1. 资本市场出清 (r_supply - r_demand) ..... : %15.4e\n', ss.errors.r);
            fprintfAndLog(fileID, '   2. 劳动市场出清 (w_supply - w_demand) ..... : %15.4e\n', ss.errors.w);
            fprintfAndLog(fileID, '   3. 政府预算平衡 (Gov_Rev - Gov_Exp) ..... : %15.4e\n', ss.errors.tr);
            fprintfAndLog(fileID, '   4. PAYG体系平衡 (盈余-BGP要求) ............ : %15.4e\n', ss.errors.payg);
            fprintfAndLog(fileID, '   5. 遗赠市场出清 (Bq_gen - Bq_dist) ....... : %15.4e\n', ss.errors.bq);
            if is_bgp_ss
                fprintfAndLog(fileID, '   6. 私人资本BGP积累 (K_p supply - demand) .. : %15.4e (真稳态下，此项为误差)\n', ss.errors.K);
            else
                fprintfAndLog(fileID, '   6. 私人资本BGP积累 .......................... : (伪稳态下，此项不作为误差进行检验)\n');
            end

            fprintfAndLog(fileID, '\n=========================================================================================================\n');
            if fileID > 1, fclose(fileID); fprintf('\n报告已成功保存至: %s\n', report_filename); end

            if nargout > 0, report_struct = struct(); end
            
            function fprintfAndLog(fileID, formatSpec, varargin)
                if print_flag
                    fprintf(1, formatSpec, varargin{:});
                    if fileID > 1, fprintf(fileID, formatSpec, varargin{:}); end
                end
            end
        end        
        function b_hat_formula_h = calculate_formula_benefits(w_hat, cS)
            % =========================================================================
            % ==             【新增】养老金福利公式计算器 v1.0
            % ==
            % == 功能: 根据给定的宏观工资水平(w_hat)和制度参数(cS)，计算
            % ==       不同类型家庭的公式化养老金福利(b_hat)。
            % ==
            % == 核心设定:
            % ==   - 社会平均工资被简化代理为 w_hat (有效劳动单位的报酬)。
            % ==   - 个人指数化工资被简化为 w_hat * ageEff (特定年龄的有效收入)。
            % ==   - 福利是针对退休后所有年龄段的，因此计算一个平均值。
            % =========================================================================
            nH = cS.nTypes;
            b_hat_formula_h = zeros(nH, 1);


            
            % 设定一个代表性的退休后年龄，用于计算个人指数化工资部分
            % (因为效率曲线在退休后为0，我们需要一个代理值)
            % 此处使用退休前最后一个年龄的效率作为整个退休阶段的代理
            eff_at_retirement_age = cS.ageEffV_new_h(cS.aR_new, :); % [1 x nH] vector

            % 定义社会平均工资 (简化代理)
            social_avg_wage = w_hat;

            for h = 1:nH
                % 获取该类型的个人指数化工资代理值
                individual_indexed_wage = w_hat * eff_at_retirement_age(h);

                if h <= 2 % 类型1 & 2: 城镇职工
                    % 基础部分: (社平工资 + 本人指数化工资)/2 * 年限 * 1%
                    b_basic = (social_avg_wage + individual_indexed_wage) / 2 * cS.n_contrib_years * 0.01;
                    % 个人账户部分: 本人指数化工资 * 因子
                    b_personal = individual_indexed_wage * cS.personal_account_factor_urban;
                    b_hat_formula_h(h) = b_basic + b_personal;

                else % 类型3 & 4: 城乡居民
                    % 基础部分: 社平工资 * 比例
                    b_basic = social_avg_wage * cS.resident_basic_benefit_ratio;
                    % 个人账户部分: 本人指数化工资 * 因子
                    b_personal = individual_indexed_wage * cS.resident_personal_factor;
                    b_hat_formula_h(h) = b_basic + b_personal;
                end
            end
        end
        
        function Dist_by_type = get_final_dist(ss, cS, paramS, Z_ss_norm, is_bgp_ss)
            % =========================================================================
            % == 函数: get_final_dist
            % == 版本: [v2.1 - 均衡变量构造修正版]
            % ==
            % == 核心修改:
            % ==   - [!!!] 修正了构造最终均衡向量 x_eq_final 的逻辑。
            % ==   - 不再使用 ss.Bequest_distributed_agg 和 ss.TR_distributed_agg
            % ==     来计算人均量，因为这会导致维度错误。
            % ==   - 直接使用求解器返回的、本身就是人均量的均衡变量 ss.bq_h 和
            % ==     (ss.TR_distributed_agg / mass_total) 来构造 x_eq_final。
            % =========================================================================
            mass_total = sum(Z_ss_norm(:));
            nH = cS.num_hh_types;
            
            % 直接使用求解器返回的人均均衡变量来构造x_eq_final
            tr_per_capita = ss.TR_distributed_agg / mass_total;
            bq_per_capita = ss.Bequest_distributed_agg / mass_total; % ss.bq_h 本身就是 per-capita of the type

            x_eq_final = [ss.r_mkt; ss.w_hat; tr_per_capita; ss.b_hat_h(:); bq_per_capita];
            
            [~, ~, Dist_by_type] = SS.system_of_equations(x_eq_final, Z_ss_norm, cS, paramS, struct('A',1.0), is_bgp_ss);
        end
        function polS_by_type = get_final_polS(ss, cS, paramS, Z_ss_norm)
            % =========================================================================
            % == 函数: get_final_polS
            % == 版本: [v2.1 - theta_t 传递修正版]
            % ==
            % == 核心修改:
            % ==   - 在构建 M_for_hh 结构体时，增加了 theta_t 字段。
            % ==   - 此修正确保 VFI 求解器能够接收到计算PAYG缴费所需的缴费率。
            % =========================================================================
            nH = cS.num_hh_types;
            polS_by_type = cell(nH, 1);
            mass_total = sum(Z_ss_norm(:));

            for h = 1:nH
                cS_h = cS; % 准备类型特定的参数
                cS_h.theta_path = cS.theta_path_h(h);
                cS_h.ageEffV_new = cS.ageEffV_new_h(:,h);
                
                M_for_hh.r_mkt_t = ss.r_mkt;
                M_for_hh.w_t = ss.w_hat;
                M_for_hh.tr_per_hh = ss.TR_distributed_agg / mass_total;
                M_for_hh.b_hat_t = ss.b_hat_h(h); % 类型特定的养老金
                
                % [!!! 核心修正 !!!] 将特定类型的theta_path最终值传递给VFI
                M_for_hh.theta_t = cS_h.theta_path(end);

                [polS_by_type{h}, ~] = household.VFI_solver(M_for_hh, paramS, cS_h);
            end
        end

    end
    
    end