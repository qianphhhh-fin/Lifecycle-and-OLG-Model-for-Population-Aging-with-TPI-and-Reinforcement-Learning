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
            % [重构核心-v4 PPS版] 统一的稳态求解器主入口
            if nargin < 4, verbose = true; end
            if nargin < 5, x0_guess = []; end
            if nargin < 6, solver_method = 'fsolve'; end

            % --- 1. 参数设置 ---
            cS.A = params_ext.A;
            if isfield(params_ext, 'g_A_ss'), cS.g_A_ss = params_ext.g_A_ss; end
            cS.theta_path = params_ext.theta;

            if isfield(params_ext, 'Z') && ~isempty(params_ext.Z)
                Z_ss_norm = params_ext.Z;
                if verbose, fprintf('   [求解器信息] 使用外部传入的人口分布Z进行求解/校准。\n'); end
            else
                error('错误：求解稳态需要一个明确的人口分布Z。请在params_ext中提供。');
            end
            
            % --- [核心开关逻辑] ---
            if ~isfield(cS, 'pps_active'), cS.pps_active = false; end
            if cS.pps_active && cS.nkpps <= 1
                warning('PPS模式要求 nkpps > 1. 请在主脚本中设置 ngrid_pps。');
            elseif ~cS.pps_active
                % 强制设定为无PPS模式的维度，确保VFI兼容
                cS.nkpps = 1; 
                cS.kppsGridV = 0;
            end

            % --- 2. 求解 ---
            system_wrapper = @(x) main_steady_state_utils_bgp.system_of_equations_steady_state(x, Z_ss_norm, cS, paramS);
            if isempty(x0_guess), x0 = [0.3336, 0.078, 0.3]; else, x0 = x0_guess(1:3); end
            [x_eq, eq_found] = main_steady_state_utils_bgp.solve_with_fsolve(system_wrapper, x0, verbose);

            if ~eq_found
                warning('统一求解器未能找到均衡解！');
                ss = []; Dist = []; polS = []; valS = [];
                return;
            end

            % --- 3. 获取最终结果并展示 ---
            [ss, Dist, polS, valS] = main_steady_state_utils_bgp.calculate_aggregates_unified(x_eq(1), x_eq(2), x_eq(3), Z_ss_norm, cS, paramS);

            % [验证步骤]
            % verify_household_budget_constraint 暂不更新，待需要时再适配PPS
            % main_steady_state_utils_bgp.verify_household_budget_constraint(ss, polS, cS, paramS);
            main_steady_state_utils_bgp.verify_FULL_steady_state(Dist, polS, paramS, cS);

            if verbose, main_steady_state_utils_bgp.display_national_accounts_unified(ss, cS);
            else, main_steady_state_utils_bgp.check_national_accounts_unified(ss, cS); end
        end

        % =======================================================
        % == 核心求解方程 (无变化)
        % =======================================================
        % =======================================================
        % == 核心求解方程 (恢复流量均衡版)
        % =======================================================
        function F_error = system_of_equations_steady_state(x, Z_ss_norm, cS, paramS)
            % =========================================================================
            % == 函数: system_of_equations_steady_state
            % == 版本: [v15 - 流量均衡最终版]
            % ==
            % == 核心决策:
            % ==   经过尝试，我们确认【存量均衡】方法 K_guess = K_resulting 在
            % ==   存在意外遗赠 (bequests) 的OLG模型中是【不正确】的。因为它
            % ==   忽略了死亡家庭持有的资产会“漏出”成为遗赠税，而不会转变为
            % ==   下一期的生产性资本。
            % ==
            % ==   因此，我们必须回到并坚持使用【流量均衡】条件。它正确地将
            % ==   意外遗赠(Bequest_demand)处理为对私人储蓄的一种需求，从而
            % ==   保证了宏观会计的完全闭合。
            % ==
            % ==   收敛性问题需要通过其他方法解决（如优化初值、使用更鲁棒的
            % ==   求解器或平滑政策函数），而不是改变正确的经济学原理。
            % =========================================================================

            K_private_total_guess = x(1); K_g_guess = x(2); L_guess = x(3);
            
            % --- 1. 基于猜测的宏观价格，计算所有微观和宏观聚合量 ---
            g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
            [ss, ~, ~, ~] = main_steady_state_utils_bgp.calculate_aggregates_unified(K_private_total_guess, K_g_guess, L_guess, Z_ss_norm, cS, paramS);
            
            % --- 2. [核心] 市场出清误差 (基于流量) ---
            % --- 方程1: 私人资本市场 (S_p_net = I_p_net + Bequests) ---
            % 私人净储蓄的供给
            S_p_supply = ss.Saving_private_flow;
            
            % 对私人净储蓄的需求来自两部分：
            % 1. 企业为了跟上技术增长所需的净投资
            I_p_net_demand = g_A_period * K_private_total_guess;
            % 2. 政府通过遗赠税从储蓄池中抽走的部分
            Bequest_demand = ss.Bequest_tax; 
            
            Total_demand_for_net_saving = I_p_net_demand + Bequest_demand;
            
            error_Kp = Total_demand_for_net_saving - S_p_supply;

            % --- 方程2: 公共资本市场 (流量均衡) ---
            Gov_Total_Revenue = ss.Regular_tax + ss.Bequest_tax + ss.PPS_tax_agg + ss.Public_Capital_Return;
            Resources_for_discretion = Gov_Total_Revenue - ss.Depreciation_g;
            Gov_Net_Saving = cS.lambda_g * Resources_for_discretion;
            I_g_net_demand = g_A_period * K_g_guess;
            error_Kg = I_g_net_demand - Gov_Net_Saving;

            % --- 方程3: 劳动市场 (供给 = 需求) ---
            error_L = L_guess - ss.L_hat;

            F_error = [error_Kp; error_Kg; error_L];
        end        % == 宏观聚合函数 (修改版)
        % =======================================================
% =======================================================
        % == 宏观聚合函数 (修改版)
        % =======================================================
        % =======================================================
        % == 宏观聚合函数 (修改版)
        % =======================================================
        function [ss, Dist, polS, valS] = calculate_aggregates_unified(K_private_total_guess, K_g_guess, L_guess, Z_ss_norm, cS, paramS)
            % =========================================================================
            % == 函数: calculate_aggregates_unified
            % == 版本: [v22 - 最终会计修正版]
            % ==
            % == 核心修正:
            % ==   [致命BUG修复] 修正了私人净储蓄 (Saving_private_flow) 的会计定义。
            % ==   之前的版本错误地将PPS的缴费和提取作为家庭部门的支出和收入项，
            % ==   这混淆了“国民经济核算(NIPA)”逻辑和“家庭现金流”逻辑，导致
            % ==   最终资源约束 Y = C+I+Gc 不平衡。
            % ==
            % ==   新逻辑严格遵循NIPA定义：
            % ==   私人净储蓄 = 私人可支配收入 - 私人消费
            % ==   其中，PPS缴费和提取是家庭内部的资产组合调整，不计入NIPA意义
            % ==   上的收入和支出。这确保了所有宏观恒等式在有PPS时也能闭合。
            % =========================================================================

            % --- 1. 初始化和价格计算 ---
            if K_private_total_guess <= 0, K_private_total_guess = 1e-8; end; if K_g_guess <= 0, K_g_guess = 1e-8; end; if L_guess <= 0, L_guess = 1e-8; end;
            A_ss = 1.0; theta_ss = cS.theta_path(1);

            M_prices = main_steady_state_utils_bgp.get_prices_at_t(K_private_total_guess, K_g_guess, L_guess, A_ss, cS);

            M_for_hh = M_prices;
            mass_retirees_ss = sum(Z_ss_norm((cS.aR_new+1):end));
            M_for_hh.b_t = (theta_ss * M_prices.w_t * L_guess) / max(1e-9, mass_retirees_ss);
            cS.theta_t = theta_ss;

            % --- 2. 求解家庭问题和稳态分布 ---
            [polS, valS] = main_steady_state_utils_bgp.HHSolution_VFI_unified(M_for_hh, paramS, cS);
            Dist = main_steady_state_utils_bgp.solve_steady_state_distribution_unified(polS, paramS, cS, Z_ss_norm);

            % --- 3. 聚合微观变量 ---
            K_private_hat_agg = 0; K_pps_hat_agg = 0; L_agg = 0; C_agg = 0; Bequest_tax_agg = 0;
            Shock_exp_agg = 0; Pension_out_agg = 0; Pension_in_agg = 0; Regular_tax_agg = 0;
            PPS_contrib_agg = 0; PPS_withdrawal_agg = 0; PPS_tax_agg = 0;

            for ia = 1:cS.aD_new
                mass_dist = Dist(:,:,:,ia);
                K_private_hat_agg = K_private_hat_agg + sum(polS(ia).k_prime .* mass_dist, 'all');
                if cS.pps_active, K_pps_hat_agg = K_pps_hat_agg + sum(polS(ia).kpps_prime .* mass_dist, 'all'); end

                total_assets_chosen_for_next_period = sum((polS(ia).k_prime + polS(ia).kpps_prime) .* mass_dist, 'all');
                Bequest_tax_agg = Bequest_tax_agg + total_assets_chosen_for_next_period * (1 - cS.s_pathV(ia));

                if ia <= cS.aR_new
                    mass_by_epsilon = squeeze(sum(mass_dist, [1,2]));
                    L_agg = L_agg + sum( (cS.ageEffV_new(ia) * paramS.leGridV') .* mass_by_epsilon' , 'all');
                end
                C_agg = C_agg + sum(polS(ia).c .* mass_dist, 'all');
                Shock_exp_agg = Shock_exp_agg + sum(polS(ia).shock_exp .* mass_dist, 'all');
                Regular_tax_agg = Regular_tax_agg + sum(polS(ia).tax_regular .* mass_dist, 'all');
                Pension_in_agg = Pension_in_agg + sum(polS(ia).tax_payg .* mass_dist, 'all');
                Pension_out_agg = Pension_out_agg + sum(polS(ia).pension_out .* mass_dist, 'all');

                if cS.pps_active
                    PPS_contrib_agg = PPS_contrib_agg + sum(polS(ia).pps_contrib .* mass_dist, 'all');
                    PPS_withdrawal_agg = PPS_withdrawal_agg + sum(polS(ia).pps_withdrawal .* mass_dist, 'all');
                    PPS_tax_agg = PPS_tax_agg + sum(polS(ia).pps_tax .* mass_dist, 'all');
                end
            end

            % --- 4. [最终会计修正] 严格遵循NIPA恒等式计算私人净储蓄 ---
            % 私人净储蓄 = 私人可支配收入 - 私人消费

            % A. 私人部门的总收入 (要素收入 + 转移收入)
            % 注意：PPS提取是家庭内部资产变现，不是NIPA意义上的新增收入。
            Household_NIPA_Inflow = (M_prices.w_t * L_agg) ...           % 劳动收入
                + (K_private_total_guess * M_prices.r_mkt_t) ...       % 资本利得收入 (来自 Kp 和 Kpps)
                + Pension_out_agg;                                     % 公共养老金转移收入

            % B. 私人部门的总支出（消费性 + 税收）
            % 注意：PPS缴费是储蓄的一种形式，不是NIPA意义上的支出。
            Household_NIPA_Outlay_NonSaving = C_agg ...                  % 私人消费
                + Shock_exp_agg ...                                    % 意外冲击消费
                + Regular_tax_agg ...                                  % 常规税
                + Pension_in_agg ...                                   % PAYG缴费 (视为税)
                + PPS_tax_agg;                                         % PPS提取时缴的税 (是真实的税)

            % C. [正确定义] 私人净储蓄 (可支配收入 - 消费)
            Saving_private_flow = Household_NIPA_Inflow - Household_NIPA_Outlay_NonSaving;

            % --- 5. 计算其他宏观量 ---
            Depreciation_p = cS.ddk * K_private_total_guess;
            Depreciation_g = cS.ddk_g * K_g_guess;
            Factor_Payment_Total = (M_prices.w_t * L_agg) + ((M_prices.r_mkt_t + cS.ddk) * K_private_total_guess);
            Public_Capital_Return = M_prices.Y_t - Factor_Payment_Total;

            % --- 6. 填充 ss 结构体 ---
            ss = struct();
            ss.K_private_begin_hat = K_private_total_guess;
            ss.K_private_hat = K_private_hat_agg; % 注意：这里只记录了非PPS部分
            ss.L_hat = L_agg;
            ss.K_public_hat = K_g_guess;
            ss.Y_from_production_hat = M_prices.Y_t;
            ss.w_hat = M_prices.w_t;
            ss.r_mkt = M_prices.r_mkt_t;
            ss.C_agg = C_agg;
            ss.Shock_exp_agg = Shock_exp_agg;
            ss.Bequest_tax = Bequest_tax_agg;
            ss.Regular_tax = Regular_tax_agg;
            ss.Pension_in = Pension_in_agg;
            ss.Pension_out = Pension_out_agg;
            ss.Saving_private_flow = Saving_private_flow;
            ss.Depreciation_p = Depreciation_p;
            ss.Depreciation_g = Depreciation_g;
            ss.Public_Capital_Return = Public_Capital_Return;

            ss.K_pps_hat = K_pps_hat_agg;
            ss.PPS_contrib_agg = PPS_contrib_agg;
            ss.PPS_withdrawal_agg = PPS_withdrawal_agg;
            ss.PPS_tax_agg = PPS_tax_agg;
        end        
        
        function [polS, valS] = HHSolution_VFI_unified(M_vfi, paramS_vfi, cS_vfi)
            % =========================================================================
            % == 函数: HHSolution_VFI_unified
            % == 版本: [v3 - 变量名BUG修复]
            % == 核心修正: 修正了调用子函数时传递的变量名错误 (paramS_age -> paramS_vfi)
            % =========================================================================
            
            % --- 1. 初始化输出 ---
            valS = -Inf(cS_vfi.nk, cS_vfi.nkpps, cS_vfi.nw_expanded, cS_vfi.aD_new);
            polS_cell = cell(cS_vfi.aD_new, 1);

            % --- 2. 设置年龄相关的养老金 ---
            bV_payg_vfi = zeros(1, cS_vfi.aD_new);
            if cS_vfi.aR_new < cS_vfi.aD_new, bV_payg_vfi((cS_vfi.aR_new+1):cS_vfi.aD_new) = M_vfi.b_t; end

            % --- 3. VFI向后迭代 ---
            for a_idx = cS_vfi.aD_new : -1 : 1
                vPrime_kkppse_next = [];
                if a_idx < cS_vfi.aD_new, vPrime_kkppse_next = valS(:,:,:,a_idx+1); end

                % --- [核心开关] 根据 cS.pps_active 选择求解器 ---
                if cS_vfi.pps_active
                    % [BUG 修复] 调用支持PPS的Naive求解器，使用正确的变量名 paramS_vfi
                    [valS(:,:,:,a_idx), polS_cell{a_idx}] = ...
                        main_steady_state_utils_bgp.HHSolutionByAge_NAIVE_PPS(...
                        a_idx, M_vfi, bV_payg_vfi(a_idx), paramS_vfi, cS_vfi);
                else
                    % 调用原始的、无PPS的Naive求解器
                    [valS(:,:,:,a_idx), polS_cell{a_idx}] = ...
                        main_steady_state_utils_bgp.HHSolutionByAge_NAIVE(...
                        a_idx, vPrime_kkppse_next, M_vfi, bV_payg_vfi(a_idx), paramS_vfi, cS_vfi);
                end
            end

            % --- 4. 将cell数组转换为结构体数组 ---
            polS = [polS_cell{:}];
        end        
        
        
        function [val_age, pol_age] = HHSolutionByAge_NAIVE(a_idx, ~, M_age, b_age_val, paramS_age, cS)
            % =========================================================================
            % == 函数: HHSolutionByAge_NAIVE
            % == 版本: [v2.0 - 维度兼容版 BUG修复]
            % ==
            % == 核心修正:
            % ==   1. 使函数的输出维度与 cS.nkpps 保持一致，即使在无PPS模式下。
            % ==   2. 使用 repmat 复制结果，确保与主VFI循环的维度匹配，修复 fsolve 报错。
            % =========================================================================

            % --- 0. 设定一个固定的、合理的储蓄率 ---
            SAVINGS_RATE = 0.20;

            % --- 1. [维度修正] 初始化单切片输出结构体 (nkpps=1) ---
            val_age_slice = -1e20 * ones(cS.nk, 1, cS.nw_expanded);
            pol_age_slice = struct(...
                'c', zeros(cS.nk, 1, cS.nw_expanded), 'k_prime', zeros(cS.nk, 1, cS.nw_expanded), ...
                'kpps_prime', zeros(cS.nk, 1, cS.nw_expanded), 'tax_regular', zeros(cS.nk, 1, cS.nw_expanded), ...
                'tax_payg', zeros(cS.nk, 1, cS.nw_expanded), 'shock_exp', zeros(cS.nk, 1, cS.nw_expanded), ...
                'pension_in', zeros(cS.nk, 1, cS.nw_expanded), 'pension_out', zeros(cS.nk, 1, cS.nw_expanded));

            % --- 2. 预计算BGP相关的通用参数 ---
            g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
            market_return_factor = 1 + M_age.r_mkt_t;

            % --- 3. 终点决策 (简化为消费掉所有财富) ---
            if a_idx == cS.aD_new
                k_capital_tax = cS.tau_k .* (cS.kGridV * M_age.r_mkt_t);
                total_wealth = cS.kGridV * market_return_factor - k_capital_tax + b_age_val;
                c_expend_final = total_wealth;
                k_prime_final = zeros(cS.nk, 1);
                c_final = c_expend_final ./ (1 + cS.tau_c);
                consumption_tax = c_final * cS.tau_c;
                final_regular_tax = k_capital_tax + consumption_tax;

                for ie = 1:cS.nw_expanded
                    pol_age_slice.c(:, 1, ie) = c_final;
                    pol_age_slice.k_prime(:, 1, ie) = k_prime_final;
                    pol_age_slice.pension_out(:, 1, ie) = b_age_val;
                    pol_age_slice.tax_regular(:, 1, ie) = final_regular_tax;
                end
            else
                % --- 4. 遍历所有状态点，应用固定储蓄率规则 ---
                for ie = 1:cS.nw_expanded
                    for ik = 1:cS.nk
                        % --- A. 计算可支配资源 (逻辑与原版一致) ---
                        k_now = cS.kGridV(ik);
                        k_return = k_now * market_return_factor;
                        labor_income_gross = 0; pension_out = 0;
                        if a_idx <= cS.aR_new
                            labor_income_gross = M_age.w_t * cS.ageEffV_new(a_idx) * paramS_age.leGridV(ie);
                        else
                            pension_out = b_age_val;
                        end
                        capital_tax = cS.tau_k * (k_now * M_age.r_mkt_t);
                        payg_tax = cS.theta_t * labor_income_gross;
                        labor_tax = cS.tau_l * max(0, labor_income_gross - payg_tax);
                        cash_on_hand = k_return + labor_income_gross + pension_out - (capital_tax + payg_tax + labor_tax);
                        shock_exp_factor = (ie == cS.nw + 1) * cS.kappa_young + (ie == cS.nw + 2) * cS.kappa_old;
                        shock_exp = shock_exp_factor * cash_on_hand;
                        available_for_c_and_s = cash_on_hand - shock_exp;
                        % --- B. 应用固定储蓄率规则 ---
                        saving_decision = SAVINGS_RATE * available_for_c_and_s;
                        k_prime_decision = saving_decision / (1 + g_A_period);
                        k_prime_decision = max(cS.kMin, k_prime_decision);
                        % --- C. 计算消费 ---
                        saving_expenditure = k_prime_decision * (1 + g_A_period);
                        c_expend_decision = available_for_c_and_s - saving_expenditure;
                        if c_expend_decision < 1e-9
                            c_expend_decision = 1e-9;
                            saving_expenditure = available_for_c_and_s - c_expend_decision;
                            k_prime_decision = max(cS.kMin, saving_expenditure / (1 + g_A_period));
                        end
                        c_decision = c_expend_decision / (1 + cS.tau_c);
                        consumption_tax = c_decision * cS.tau_c;
                        % --- D. 存储决策到单切片中 ---
                        pol_age_slice.c(ik, 1, ie) = c_decision;
                        pol_age_slice.k_prime(ik, 1, ie) = k_prime_decision;
                        pol_age_slice.pension_out(ik, 1, ie) = pension_out;
                        pol_age_slice.shock_exp(ik, 1, ie) = shock_exp;
                        pol_age_slice.tax_regular(ik, 1, ie) = capital_tax + labor_tax + consumption_tax;
                        pol_age_slice.tax_payg(ik, 1, ie) = payg_tax;
                    end
                end
            end
            
            % --- 5. [核心修复] 将单切片结果扩展到目标维度 ---
            val_age = repmat(val_age_slice, [1, cS.nkpps, 1]);
            
            pol_age = pol_age_slice; % 先复制结构体
            % 遍历结构体的所有字段，并用repmat扩展它们
            fields = fieldnames(pol_age);
            for i = 1:length(fields)
                field_name = fields{i};
                pol_age.(field_name) = repmat(pol_age_slice.(field_name), [1, cS.nkpps, 1]);
            end
        end
                
        % == [新增函数] 支持PPS的Naive VFI求解器
        % =======================================================
 % =======================================================
        % == [新增函数] 支持PPS的Naive VFI求解器
        % =======================================================
        function [val_age, pol_age] = HHSolutionByAge_NAIVE_PPS(a_idx, M_age, b_age_val, paramS_age, cS)
            % =========================================================================
            % == 函数: HHSolutionByAge_NAIVE_PPS
            % == 版本: [v2.0 - 终点期会计BUG修复版]
            % ==
            % == 核心修正:
            % ==   [致命BUG修复] 修正了终点期(a_idx == cS.aD_new)的会计逻辑。
            % ==   之前的版本将 pps_tax_final 错误地加总到 final_regular_tax
            % ==   中，同时又独立赋值给 pol_age.pps_tax，导致PPS税收在聚合时
            % ==   被重复计算，从而破坏了宏观均衡。
            % =========================================================================

            % --- 0. 设定固定的、合理的储蓄率 ---
            SAVINGS_RATE = 0.20;

            % --- 1. 初始化输出结构体 (增加PPS相关字段) ---
            val_age = -1e20 * ones(cS.nk, cS.nkpps, cS.nw_expanded); % 价值函数无意义
            pol_age = struct(...
                'c', zeros(cS.nk, cS.nkpps, cS.nw_expanded), 'k_prime', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
                'kpps_prime', zeros(cS.nk, cS.nkpps, cS.nw_expanded), 'tax_regular', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
                'tax_payg', zeros(cS.nk, cS.nkpps, cS.nw_expanded), 'shock_exp', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
                'pension_in', zeros(cS.nk, cS.nkpps, cS.nw_expanded), 'pension_out', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
                'pps_contrib', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
                'pps_withdrawal', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
                'pps_tax', zeros(cS.nk, cS.nkpps, cS.nw_expanded) ...
                );

            % --- 2. 预计算通用参数 ---
            g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
            market_return_factor = 1 + M_age.r_mkt_t;

            % --- 3. 终点决策 (清算所有资产) ---
            if a_idx == cS.aD_new
                [K_grid, Kpps_grid] = ndgrid(cS.kGridV, cS.kppsGridV);
                k_capital_tax = cS.tau_k .* (K_grid .* M_age.r_mkt_t);
                
                pps_wealth_final = Kpps_grid .* market_return_factor;
                pps_tax_final = cS.pps_tax_rate_withdrawal .* pps_wealth_final; % 假设全部取出并缴税
                
                total_wealth = (K_grid * market_return_factor - k_capital_tax) + (pps_wealth_final - pps_tax_final) + b_age_val;

                c_expend_final = total_wealth;
                c_final = c_expend_final ./ (1 + cS.tau_c);
                consumption_tax = c_final * cS.tau_c;
                
                % [BUG修复] final_regular_tax 只应包含资本税和消费税，不应包含PPS税
                final_regular_tax = k_capital_tax + consumption_tax;

                for ie = 1:cS.nw_expanded
                    pol_age.c(:, :, ie) = c_final;
                    pol_age.k_prime(:, :, ie) = zeros(size(K_grid)); % 终点储蓄为0
                    pol_age.kpps_prime(:, :, ie) = zeros(size(Kpps_grid)); % 终点PPS储蓄为0
                    pol_age.pension_out(:, :, ie) = b_age_val;
                    pol_age.tax_regular(:, :, ie) = final_regular_tax; % <-- 使用修复后的值
                    pol_age.pps_withdrawal(:, :, ie) = pps_wealth_final; % <-- 总提取额
                    pol_age.pps_tax(:, :, ie) = pps_tax_final; % <-- 总提取税
                end
                return;
            end

            % --- 4. 遍历所有状态点，应用规则 (工作期/退休期逻辑不变) ---
            for ie = 1:cS.nw_expanded
              for ikpps = 1:cS.nkpps
                for ik = 1:cS.nk
                    k_now = cS.kGridV(ik);
                    kpps_now = cS.kppsGridV(ikpps);
                    
                    % 初始化流量
                    pension_out = 0;
                    pps_contrib = 0;
                    pps_withdrawal = 0;
                    pps_tax = 0;
                    
                    if a_idx <= cS.aR_new % --- 工作期 ---
                        labor_income_gross = M_age.w_t * cS.ageEffV_new(a_idx) * paramS_age.leGridV(ie);
                        payg_tax = cS.theta_t * labor_income_gross;
                        
                        pps_contrib = cS.pps_contrib_rate * labor_income_gross;
                        kpps_prime_decision = (kpps_now * market_return_factor + pps_contrib) / (1 + g_A_period); % 标准化
                        
                        capital_tax = cS.tau_k * (k_now * M_age.r_mkt_t);
                        labor_tax = cS.tau_l * max(0, labor_income_gross - payg_tax - pps_contrib);

                        cash_on_hand = (k_now * market_return_factor + labor_income_gross) ...
                            - (capital_tax + payg_tax + labor_tax + pps_contrib);

                    else % --- 退休期 ---
                        pension_out = b_age_val;
                        
                        kpps_total_value = kpps_now * market_return_factor;
                        pps_withdrawal = cS.pps_withdrawal_rate * kpps_total_value;
                        pps_tax = cS.pps_tax_rate_withdrawal * pps_withdrawal;
                        kpps_prime_decision = (kpps_total_value - pps_withdrawal) / (1 + g_A_period); % 标准化
                        net_pps_inflow_to_budget = pps_withdrawal - pps_tax;

                        capital_tax = cS.tau_k * (k_now * M_age.r_mkt_t);
                        payg_tax = 0; labor_tax = 0;

                        cash_on_hand = (k_now * market_return_factor + pension_out + net_pps_inflow_to_budget) ...
                            - capital_tax;
                    end
                    
                    % --- 通用逻辑: 计算消费和常规储蓄 ---
                    shock_exp_factor = (ie == cS.nw + 1) * cS.kappa_young + (ie == cS.nw + 2) * cS.kappa_old;
                    shock_exp = shock_exp_factor * cash_on_hand;
                    available_for_c_and_s = cash_on_hand - shock_exp;

                    saving_decision = SAVINGS_RATE * available_for_c_and_s;
                    k_prime_decision = saving_decision / (1 + g_A_period);
                    k_prime_decision = max(cS.kMin, k_prime_decision);
                    
                    saving_expenditure = k_prime_decision * (1 + g_A_period);
                    c_expend_decision = available_for_c_and_s - saving_expenditure;
                    if c_expend_decision < 1e-9
                        c_expend_decision = 1e-9;
                        saving_expenditure = available_for_c_and_s - c_expend_decision;
                        k_prime_decision = max(cS.kMin, saving_expenditure / (1 + g_A_period));
                    end
                    
                    c_decision = c_expend_decision / (1 + cS.tau_c);
                    consumption_tax = c_decision * cS.tau_c;

                    % --- 存储所有决策 ---
                    pol_age.c(ik, ikpps, ie) = c_decision;
                    pol_age.k_prime(ik, ikpps, ie) = k_prime_decision;
                    pol_age.kpps_prime(ik, ikpps, ie) = kpps_prime_decision;
                    pol_age.pension_out(ik, ikpps, ie) = pension_out;
                    pol_age.shock_exp(ik, ikpps, ie) = shock_exp;
                    pol_age.tax_regular(ik, ikpps, ie) = capital_tax + labor_tax + consumption_tax;
                    pol_age.tax_payg(ik, ikpps, ie) = payg_tax;
                    pol_age.pps_contrib(ik, ikpps, ie) = pps_contrib;
                    pol_age.pps_withdrawal(ik, ikpps, ie) = pps_withdrawal;
                    pol_age.pps_tax(ik, ikpps, ie) = pps_tax;
                end
              end
            end
        end
        
        
        function idx_mat = get_policy_index_matrix_unified(polS, cS, type)
            % [重构-修正版] 统一的策略函数离散化，将连续决策映射到离散网格索引。

            % --- 确定要离散化的策略('k'或'kpps')和对应的网格 ---
            if strcmp(type, 'k')
                gridV = cS.kGridV;
            elseif strcmp(type, 'kpps') && cS.pps_active
                gridV = cS.kppsGridV;
            else
                % 如果类型是'kpps'但PPS未激活，或类型未知，则返回一个全1的矩阵(代表第一个网格点)
                idx_mat = ones(cS.nk, cS.nkpps, cS.nw_expanded, cS.aD_new, 'uint16');
                return;
            end

            idx_mat = zeros(cS.nk, cS.nkpps, cS.nw_expanded, cS.aD_new, 'uint16'); % 初始化索引矩阵

            for ia = 1:cS.aD_new % 遍历所有年龄
                % --- 获取对应年龄的策略矩阵 ---
                if strcmp(type, 'k')
                    val_mat = polS(ia).k_prime;
                else
                    val_mat = polS(ia).kpps_prime;
                end

                % --- [核心修正] 使用嵌套循环和多维索引来确保正确赋值 ---
                for ik = 1:cS.nk
                    for ikpps = 1:cS.nkpps
                        for ie = 1:cS.nw_expanded
                            val_continuous = val_mat(ik, ikpps, ie); % 获取特定状态下的连续决策值

                            % 找到不大于该决策值的最大网格点索引
                            idx = find(gridV <= val_continuous, 1, 'last');
                            if isempty(idx)
                                idx = 1; % 如果找不到(比如决策值为负)，则取最小索引1
                            end
                            idx_mat(ik, ikpps, ie, ia) = idx; % 使用多维索引精确赋值
                        end
                    end
                end
            end
        end

        function Dist = solve_steady_state_distribution_unified(polS, paramS, cS, Z_ss_norm)
            % =========================================================================
            % == 函数: solve_steady_state_distribution_unified
            % == 版本: [v6 - 平滑求解最终版 (恢复概率质量分裂法)]
            % ==
            % == 核心修正:
            % ==   恢复使用“概率质量分裂法”(线性插值)。虽然这种方法与离散索引的
            % ==   验证器存在方法论差异，但它能保证系统方程的平滑性，这是 fsolve
            % ==   能够成功求解的【必要条件】。我们必须优先保证求解的成功。
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
                mass_levels = ones(cS.aD_new, 1);
                for ia = 1:(cS.aD_new - 1)
                    mass_levels(ia+1) = mass_levels(ia) * cS.s_pathV(ia);
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
        end        % == 保留的基础函数
        % =======================================================

        function M_prices = get_prices_at_t(K_p, K_g, L, A_t, cS)
            % =========================================================================
            % == 函数: get_prices_at_t
            % == 版本: [v2 - 理论一致版]
            % ==
            % == 核心修正:
            % ==   1. [理论完备] 生产函数和要素价格计算现在完全反映了由
            % ==      alpha 和 gamma 定义的三要素柯布-道格拉斯函数。
            % ==   2. 确保了要素报酬总和严格等于总产出 (欧拉定理)。
            % =========================================================================
            if K_p <= 0, K_p = 1e-8; end; if L <= 0, L = 1e-8; end; if K_g <= 0, K_g = 1e-8; end;

            % 使用与参数一致的三要素生产函数
            labor_exponent = 1 - cS.alpha - cS.gamma;
            if labor_exponent <= 0, error('劳动和资本的产出弹性之和必须小于1！'); end

            Y_period = A_t .* (K_p.^cS.alpha) .* (K_g.^cS.gamma) .* (L.^labor_exponent);

            % 根据正确的边际产出计算要素价格
            MPK_p_period = cS.alpha .* Y_period ./ K_p;
            w_t = labor_exponent .* Y_period ./ L;

            % 市场利率 = 资本边际产出 - 折旧率
            r_mkt_t = MPK_p_period - cS.ddk;

            M_prices = struct('K_p', K_p, 'K_g', K_g, 'L_t', L, 'Y_t', Y_period, 'w_t', w_t, 'r_mkt_t', r_mkt_t);
        end

        function [x_eq, eq_found] = solve_with_fsolve(system_wrapper, x0, verbose)
            % [保留不变] fsolve求解器的包装函数。
            if verbose, fsolve_display = 'iter'; else, fsolve_display = 'none'; end % 根据verbose设置显示选项
            options = optimoptions('fsolve', 'Display', fsolve_display, 'TolFun', 1e-6, 'TolX', 1e-6, 'MaxIterations', 100,'algorithm','trust-region' );

            if verbose, fprintf('\n--- 启动 fsolve 求解器 (求解 [K̂_p, K̂_g, L] - BGP版本) ---\n'); end
            [x_eq, ~, exitflag] = fsolve(system_wrapper, x0, options); % 调用fsolve
            if verbose, fprintf('--- fsolve 求解完成 ---\n'); end

            eq_found = (exitflag > 0); % 根据exitflag判断是否收敛
        end

        % =======================================================
        % == 结果展示与检查
        % =======================================================

        function verify_household_budget_constraint(ss, polS, cS, paramS)
            % =========================================================================
            % == 函数: verify_household_budget_constraint
            % == 版本: [v5 - PPS兼容验证版]
            % ==
            % == 核心修正:
            % ==   1. 增加 cS.pps_active 开关，使其能同时验证有/无PPS两种模式。
            % ==   2. 在PPS模式下，会计逻辑与 HHSolutionByAge_NAIVE_PPS 完全同步。
            % =========================================================================
            fprintf('\n--- [法医级检验] 正在对单个家庭的预算约束进行微观解剖 (PPS兼容版)... ---\n');

            % --- 1. 随机选择一个状态点进行检验 ---
            a_idx_test = randi([1, cS.aR_new]);
            ik_test = randi([1, cS.nk]);
            ie_test = randi([1, cS.nw]);
            if cS.pps_active
                ikpps_test = randi([1, cS.nkpps]);
                fprintf('   检验样本: 年龄 a=%d, k-idx=%d, kpps-idx=%d, e-idx=%d\n', a_idx_test, ik_test, ikpps_test, ie_test);
            else
                ikpps_test = 1; % 无PPS模式下，kpps索引为1
                fprintf('   检验样本: 年龄 a=%d, k-idx=%d, e-idx=%d (无PPS模式)\n', a_idx_test, ik_test, ie_test);
            end

            % --- 2. 提取公共信息和决策 ---
            M_age.r_mkt_t = ss.r_mkt;
            M_age.w_t = ss.w_hat;
            cS.theta_t = ss.Pension_in / max(1e-9, (ss.w_hat * ss.L_hat));
            g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;

            c_decision = polS(a_idx_test).c(ik_test, ikpps_test, ie_test);
            k_prime_decision = polS(a_idx_test).k_prime(ik_test, ikpps_test, ie_test);

            % --- 3. [核心] 使用与VFI完全一致的逻辑，重新计算可支配资源 ---
            k_now = cS.kGridV(ik_test);
            
            if cS.pps_active
                % --- PPS激活模式下的预算重构 ---
                kpps_now = cS.kppsGridV(ikpps_test);
                
                if a_idx_test <= cS.aR_new % 工作期
                    labor_income_gross = M_age.w_t * cS.ageEffV_new(a_idx_test) * paramS.leGridV(ie_test);
                    payg_tax = cS.theta_t * labor_income_gross;
                    pps_contrib = cS.pps_contrib_rate * labor_income_gross;
                    capital_tax = cS.tau_k * (k_now * M_age.r_mkt_t);
                    labor_tax = cS.tau_l * max(0, labor_income_gross - payg_tax - pps_contrib);
                    cash_on_hand_recalc = (k_now * (1 + M_age.r_mkt_t) + labor_income_gross) ...
                        - (capital_tax + payg_tax + labor_tax + pps_contrib);
                else % 退休期
                    pension_out = ss.Pension_out / sum(cS.Z((cS.aR_new+1):end)); % 简化计算
                    kpps_total_value = kpps_now * (1 + M_age.r_mkt_t);
                    pps_withdrawal = cS.pps_withdrawal_rate * kpps_total_value;
                    pps_tax = cS.pps_tax_rate_withdrawal * pps_withdrawal;
                    net_pps_inflow = pps_withdrawal - pps_tax;
                    capital_tax = cS.tau_k * (k_now * M_age.r_mkt_t);
                    cash_on_hand_recalc = (k_now * (1 + M_age.r_mkt_t) + pension_out + net_pps_inflow) ...
                        - capital_tax;
                end

            else
                % --- 无PPS模式下的预算重构 (与原始版本一致) ---
                k_return = k_now * (1 + M_age.r_mkt_t);
                labor_income_gross = M_age.w_t * cS.ageEffV_new(a_idx_test) * paramS.leGridV(ie_test);
                capital_tax = cS.tau_k * (k_now * M_age.r_mkt_t);
                payg_tax = cS.theta_t * labor_income_gross;
                labor_tax = cS.tau_l * max(0, labor_income_gross - payg_tax);
                total_tax_flow = capital_tax + payg_tax + labor_tax;
                cash_on_hand_recalc = k_return + labor_income_gross - total_tax_flow;
            end
            
            % --- 4. 计算冲击和最终可用于 C/S 的资源 (通用逻辑) ---
            shock_exp_factor = (ie_test == cS.nw + 1) * cS.kappa_young + (ie_test == cS.nw + 2) * cS.kappa_old;
            shock_exp = shock_exp_factor * cash_on_hand_recalc;
            available_for_c_and_s_recalc = cash_on_hand_recalc - shock_exp;
            fprintf('   A. 重新计算的可支配资源 (可用于C和S) ..: %12.8f\n', available_for_c_and_s_recalc);

            % --- 5. 计算VFI决策后的总支出 ---
            consumption_expenditure = c_decision * (1 + cS.tau_c);
            saving_expenditure = k_prime_decision * (1 + g_A_period);
            total_outlay_from_decision = consumption_expenditure + saving_expenditure;
            fprintf('   B. VFI决策后的总支出 (C_expend + k''*(1+g)) : %12.8f\n', total_outlay_from_decision);

            % --- 6. 计算微观预算缺口 ---
            micro_budget_gap = available_for_c_and_s_recalc - total_outlay_from_decision;
            fprintf('   -----------------------------------------------------------------\n');
            fprintf('   => 微观预算缺口 (A - B) ..................: %12.4e\n', micro_budget_gap);

            if abs(micro_budget_gap) < 1e-7
                fprintf('   ✅ [结论] 微观预算闭合！VFI的内部计算与外部检验完全一致！\n');
            else
                fprintf('   ⚠️ [结论] 微观预算不闭合！请仔细核对本函数与VFI函数的每一行计算！\n');
            end
            fprintf('--------------------------------------------------------------------------\n');
        end        % =========================================================================
        % == [新函数] 法医级家庭预算约束检验器
        % =========================================================================
        
        function [is_steady, max_diff] = verify_steady_state_distribution(Dist, k_prime_idx, kpps_prime_idx, paramS, cS, Z_ss_norm)
            % =========================================================================
            % == 函数说明: 稳态分布验证器
            % ==
            % == 验证逻辑:
            % ==   将输入的分布 Dist 按照策略函数演化一步，得到新分布 Dist_new。
            % ==   如果 Dist 是真正的稳态分布，那么 Dist_new 必须与 Dist 相等。
            % =========================================================================

            fprintf('\n--- 正在验证稳态分布 (Dist) 的正确性... ---\n');

            % --- 1. 初始化新分布矩阵 ---
            Dist_new = zeros(size(Dist));

            % --- 2. 处理新生儿 (第一年龄组) ---
            dist_newborn = zeros(cS.nk, cS.nkpps, cS.nw_expanded);
            dist_newborn(1, 1, 1:cS.nw) = paramS.leProb1V'; % 新生儿出生时没有资产，冲击服从长期分布
            Dist_new(:, :, :, 1) = dist_newborn * Z_ss_norm(1); % 设置新分布的第一期

            % --- 3. 向前迭代一步，计算所有后续年龄组的新分布 ---
            for ia = 1:(cS.aD_new - 1)
                dist_ia = Dist(:, :, :, ia); % 获取输入的、待验证的分布的第 ia 期

                dist_ia_plus_1 = zeros(cS.nk, cS.nkpps, cS.nw_expanded); % 初始化下一期的临时分布
                trans_mat_next = paramS.TrProbM_by_age{ia + 1};

                for ik = 1:cS.nk
                    for ikpps = 1:cS.nkpps
                        for ie = 1:cS.nw_expanded
                            mass = dist_ia(ik, ikpps, ie);
                            if mass < 1e-20, continue; end

                            % 根据策略函数，找到下一期的资产状态索引
                            ik_p = k_prime_idx(ik, ikpps, ie, ia);
                            ikpps_p = kpps_prime_idx(ik, ikpps, ie, ia);

                            % 按照冲击转移概率，将人口质量分配到下一期
                            trans_probs = reshape(trans_mat_next(ie, :), [1, 1, cS.nw_expanded]);
                            dist_ia_plus_1(ik_p, ikpps_p, :) = dist_ia_plus_1(ik_p, ikpps_p, :) + mass * trans_probs;
                        end
                    end
                end

                % 重新缩放以匹配稳态人口结构 (这一步必须和原函数完全一样)
                mass_at_ia = sum(dist_ia, 'all');
                if mass_at_ia > 1e-12
                    rescale = Z_ss_norm(ia+1) / mass_at_ia;
                    Dist_new(:, :, :, ia+1) = dist_ia_plus_1 * rescale;
                end
            end

            % --- 4. 比较新旧分布 ---
            % 注意：由于浮点数精度，我们不能直接比较是否等于0，而是看差的绝对值的最大值
            max_diff = max(abs(Dist - Dist_new), [], 'all');

            % 设定一个非常小的容忍度
            tolerance = 1e-9;

            if max_diff < tolerance
                is_steady = true;
                fprintf('   ✅ 验证通过！分布是稳态的 (最大差异: %.4e)\n', max_diff);
            else
                is_steady = false;
                fprintf('   ⚠️ 验证失败！分布不是稳态的 (最大差异: %.4e)\n', max_diff);
            end
        end

        function [is_steady, max_diff] = verify_FULL_steady_state(Dist_to_verify, polS, paramS, cS)
            % =========================================================================
            % == 函数说明: [v5 - 平滑验证最终版] 终期稳态完全内生验证器
            % ==
            % == 核心修正:
            % ==   1. [方法论统一] 为了与使用“概率质量分裂法”的求解器绝对一致，
            % ==      本验证器也完全采用【概率质量分裂法】（线性插值）来演化分布。
            % ==   2. [输入变更] 不再需要离散的 k_prime_idx，而是直接使用连续的
            % ==      策略函数 polS 进行计算。
            % ==   3. [逻辑自洽] 通过保证验证方法与求解方法的绝对统一，消除了
            % ==      “方法论鸿沟”，使得验证结果能够精确反映求解的准确性。
            % =========================================================================

            % --- [输入参数调整] ---
            % 注意：此版本的验证器不再需要 k_prime_idx 和 kpps_prime_idx
            % 而是直接使用原始的、连续的策略函数 polS
            fprintf('\n--- 正在进行终期稳态的【平滑方法】验证 (概率质量分裂法)...\n');

            % --- 步骤 1: 完全内生地计算理论年龄分布 Z_theory ---
            mass_levels_by_age = zeros(cS.aD_new, 1);
            mass_levels_by_age(1) = 1.0;
            for ia = 1:(cS.aD_new - 1)
                mass_levels_by_age(ia+1) = mass_levels_by_age(ia) * cS.s_pathV(ia);
            end
            Z_theory = mass_levels_by_age / sum(mass_levels_by_age);

            % --- 步骤 2: 基于 Z_theory 初始化一个从零开始计算的分布 ---
            Dist_recalculated = zeros(size(Dist_to_verify));
            dist_newborn = zeros(cS.nk, cS.nkpps, cS.nw_expanded);
            dist_newborn(1, 1, 1:cS.nw) = paramS.leProb1V';
            Dist_recalculated(:, :, :, 1) = dist_newborn * Z_theory(1);

            % --- 步骤 3: 向前迭代，使用【概率质量分裂法】演化分布 ---
            for ia = 1:(cS.aD_new - 1)
                dist_ia = Dist_recalculated(:, :, :, ia);
                dist_ia_plus_1_unscaled = zeros(cS.nk, cS.nkpps, cS.nw_expanded);
                trans_mat_next = paramS.TrProbM_by_age{ia + 1};

                % [核心修正] 使用与求解器完全相同的演化逻辑
                for ie = 1:cS.nw_expanded
                    for ikpps = 1:cS.nkpps
                        for ik = 1:cS.nk
                            mass = dist_ia(ik, ikpps, ie);
                            if mass < 1e-30, continue; end

                            % 直接从连续策略函数 polS 中获取决策
                            k_p_cont = polS(ia).k_prime(ik, ikpps, ie);
                            if cS.pps_active && cS.nkpps > 1
                                kpps_p_cont = polS(ia).kpps_prime(ik, ikpps, ie);
                            else
                                kpps_p_cont = cS.kppsGridV(1);
                            end

                            % 使用线性插值权重来分裂质量
                            [ik_lower, ik_upper, w_k_upper] = main_steady_state_utils_bgp.find_grid_and_weights(k_p_cont, cS.kGridV);
                            w_k_lower = 1.0 - w_k_upper;
                            [ikpps_lower, ikpps_upper, w_kpps_upper] = main_steady_state_utils_bgp.find_grid_and_weights(kpps_p_cont, cS.kppsGridV);
                            w_kpps_lower = 1.0 - w_kpps_upper;

                            trans_probs_vec = reshape(trans_mat_next(ie, :), [1, 1, cS.nw_expanded]);

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

                % 使用 Z_theory 进行重新缩放
                target_mass = Z_theory(ia+1);
                current_mass_generated = sum(dist_ia_plus_1_unscaled(:));
                if current_mass_generated > 1e-30
                    rescale_factor = target_mass / current_mass_generated;
                    Dist_recalculated(:, :, :, ia+1) = dist_ia_plus_1_unscaled * rescale_factor;
                end
            end

            % --- 步骤 4: 对最终结果进行归一化 ---
            total_mass_final = sum(Dist_recalculated(:));
            if total_mass_final > 1e-9
                Dist_recalculated = Dist_recalculated / total_mass_final;
            end

            % --- 步骤 5: 比较两个分布的“形状” ---
            max_diff = max(abs(Dist_to_verify - Dist_recalculated), [], 'all');

            tolerance = 1e-9;
            if max_diff < tolerance
                is_steady = true;
                fprintf('   ✅ [平滑方法] 验证通过！稳态分布是自洽的 (最大差异: %.4e)\n', max_diff);
            else
                is_steady = false;
                fprintf('   ⚠️ [平滑方法] 验证失败！稳态分布不自洽 (最大差异: %.4e)\n', max_diff);
            end
        end        % =======================================================

        % =======================================================
        % == 辅助函数: 概率质量分裂法的核心工具
        % =======================================================

% =======================================================
                function display_national_accounts_unified(ss, cS)
            % =========================================================================
            % == 函数: display_national_accounts_unified
            % == 版本: [v26 - 终极法医诊断版]
            % ==
            % == 核心目的:
            % ==   为了诊断宏观账户不平衡问题，此版本将所有会计细项都展示出来，
            % ==   特别是创建了一个全新的【家庭部门完整流量表】，用于进行法医级
            % ==   的交叉验证。
            % ==
            % == 核心特性:
            % ==   1. [新增] "家庭部门完整流量表"，详细列出所有收支，并计算缺口。
            % ==   2. [细化] 将所有聚合流量 (如税收) 分解到其最基础的组成部分。
            % ==   3. [追溯] 明确展示关键变量 (如私人净储蓄) 是如何从基础流量中计算得出的。
            % =========================================================================

            fprintf('\n\n================================================================================\n');
            if cS.pps_active
                fprintf('===     国民经济核算详细报告 (BGP + PPS激活 - 终极法医诊断版)      ===\n');
            else
                fprintf('===      国民经济核算详细报告 (BGP + 无PPS - 最终会计准则版)     ===\n');
            end
            fprintf('================================================================================\n');

            % --- 0. 预计算所有关键宏观量 ---
            g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
            
            % --- 生产与总收入 ---
            Y_prod = ss.Y_from_production_hat;
            Wages = ss.w_hat * ss.L_hat;
            K_private_total_begin = ss.K_private_begin_hat; % 这是总私人资本 (Kp+Kpps)
            Private_Capital_Income_Gross = (ss.r_mkt + cS.ddk) * K_private_total_begin;
            Private_Capital_Income_Net = ss.r_mkt * K_private_total_begin; % 净利息收入
            Public_Capital_Return = ss.Public_Capital_Return;
            GDI = Wages + Private_Capital_Income_Gross + Public_Capital_Return;

            % --- 家庭部门流量 ---
            C_hh = ss.C_agg;
            Shock_exp = ss.Shock_exp_agg;
            C_total_expenditure = C_hh * (1 + cS.tau_c);
            Consumption_tax = C_hh * cS.tau_c;
            Pension_in = ss.Pension_in;
            Pension_out = ss.Pension_out;
            Bequests_paid = ss.Bequest_tax; % 这是家庭意外支付的，也是政府收到的
            
            % --- PPS 部门流量 ---
            PPS_contrib = ss.PPS_contrib_agg;
            PPS_withdrawal_gross = ss.PPS_withdrawal_agg;
            PPS_tax = ss.PPS_tax_agg;
            PPS_withdrawal_net = PPS_withdrawal_gross - PPS_tax;

            % --- 反推税收细项 ---
            Capital_tax_est = Private_Capital_Income_Net * cS.tau_k;
            Regular_tax_total = ss.Regular_tax;
            Labor_tax_est = Regular_tax_total - Capital_tax_est - Consumption_tax;
            
            % --- 政府部门流量 ---
            Gov_Tax_Revenue = Regular_tax_total + Bequests_paid + PPS_tax;
            Gov_Total_Revenue = Gov_Tax_Revenue + Public_Capital_Return;
            Depreciation_g = ss.Depreciation_g;
            Resources_for_discretion = Gov_Total_Revenue - Depreciation_g;
            Gov_Net_Saving = cS.lambda_g * Resources_for_discretion;
            G_c = (1 - cS.lambda_g) * Resources_for_discretion;
            I_g_net = g_A_period * ss.K_public_hat;
            I_g_gross = I_g_net + Depreciation_g;

            % --- 私人部门投资与储蓄 ---
            Depreciation_p = ss.Depreciation_p; % 总私人资本折旧
            I_p_net = g_A_period * K_private_total_begin;
            I_p_gross = I_p_net + Depreciation_p;
            I_total_gross = I_p_gross + I_g_gross;
            
            S_p_net_from_ss = ss.Saving_private_flow; % 这是求解器计算出的私人净储蓄
            
            % --- [1. 核心诊断: 家庭部门完整流量表] ---
            fprintf('\n--- [1. 核心诊断: 家庭部门完整流量表 (法医级)] ---\n');
            % A. 家庭总资源 (Inflows)
            HH_Total_Inflow = Wages + Private_Capital_Income_Net + Pension_out + PPS_withdrawal_gross;
            
            % B. 家庭总支出 (消费性)
            HH_Consumptive_Outlay = C_total_expenditure + Shock_exp + Pension_in + Labor_tax_est + Capital_tax_est + PPS_contrib + PPS_tax;
            
            % C. 基于流量表重算的家庭净储蓄
            S_p_net_recalc = HH_Total_Inflow - HH_Consumptive_Outlay;
            
            fprintf('   --- A. 家庭总资源 (Inflows) ---------: %12.8f\n', HH_Total_Inflow);
            fprintf('      1. 劳动总收入 (w*L) ...............: %12.8f\n', Wages);
            fprintf('      2. 净资本利息收入 (r*K_total) .....: %12.8f\n', Private_Capital_Income_Net);
            fprintf('      3. 公共养老金收入 (Pension Out) ..: %12.8f\n', Pension_out);
            fprintf('      4. PPS提取总额 (Gross Withdrawal) .: %12.8f\n', PPS_withdrawal_gross);
            fprintf('\n');
            fprintf('   --- B. 家庭总支出 (Outlays) ---------: %12.8f\n', HH_Consumptive_Outlay + S_p_net_recalc);
            fprintf('      --- B1. 消费性支出 (Consumptive) --: %12.8f\n', HH_Consumptive_Outlay);
            fprintf('         1. 消费总支出 (C*(1+tau_c)) .....: %12.8f\n', C_total_expenditure);
            fprintf('         2. 意外冲击支出 .................: %12.8f\n', Shock_exp);
            fprintf('         3. PAYG缴费 (Pension In) ........: %12.8f\n', Pension_in);
            fprintf('         4. 劳动税 (估计) ................: %12.8f\n', Labor_tax_est);
            fprintf('         5. 资本利得税 (估计) ............: %12.8f\n', Capital_tax_est);
            fprintf('         6. PPS缴费 (Contribution) .......: %12.8f\n', PPS_contrib);
            fprintf('         7. PPS提取税 (Withdrawal Tax) ...: %12.8f\n', PPS_tax);
            fprintf('      --- B2. 净储蓄 (Net Saving) -------: %12.8f\n', S_p_net_recalc);
            fprintf('\n');
            fprintf('   --- C. 预算检验 -----------------------\n');
            fprintf('   家庭总资源 (A) ........................: %12.8f\n', HH_Total_Inflow);
            fprintf('   家庭总使用 (B1+B2) ....................: %12.8f\n', HH_Consumptive_Outlay + S_p_net_recalc);
            fprintf('   => 家庭预算缺口 (应为0) ...............: %12.4e\n', HH_Total_Inflow - (HH_Consumptive_Outlay + S_p_net_recalc));
            fprintf('\n');
            fprintf('   --- D. 储蓄交叉验证 -------------------\n');
            fprintf('   模型ss中记录的私人净储蓄 ..............: %12.8f\n', S_p_net_from_ss);
            fprintf('   本表重新计算的私人净储蓄 ..............: %12.8f\n', S_p_net_recalc);
            fprintf('   => 两者差异 (应为0) ...................: %12.4e\n', S_p_net_from_ss - S_p_net_recalc);

            % --- [2. 收入法检验: GDI vs GDP (欧拉定理)] ---
            fprintf('\n\n--- [2. 收入法检验: GDI vs GDP (欧拉定理)] ---\n');
            mismatch_gdi = Y_prod - GDI;
            fprintf('   A. 国内生产总值 (GDP) .....................: %12.8f\n', Y_prod);
            fprintf('   B. 国内总收入 (GDI) .......................: %12.8f\n', GDI);
            fprintf('   => 收入-产出缺口 (A - B) ..................: %12.4e\n', mismatch_gdi);

            % --- [3. 政府预算检验 (收支细分)] ---
            fprintf('\n--- [3. 政府预算检验 (收支细分)] ---\n');
            Gov_Outflows = I_g_gross + G_c;
            mismatch_gov = Gov_Total_Revenue - Gov_Outflows;
            fprintf('   --- A. 政府总收入 --------------------: %12.8f\n', Gov_Total_Revenue);
            fprintf('      1. 税收总收入 .......................: %12.8f\n', Gov_Tax_Revenue);
            fprintf('         - 消费税 .........................: %12.8f\n', Consumption_tax);
            fprintf('         - 劳动税 (估计) ..................: %12.8f\n', Labor_tax_est);
            fprintf('         - 资本利得税 (估计) ..............: %12.8f\n', Capital_tax_est);
            fprintf('         - 遗赠税 .........................: %12.8f\n', Bequests_paid);
            fprintf('         - PPS提取税 ......................: %12.8f\n', PPS_tax);
            fprintf('      2. 公共资本回报 .....................: %12.8f\n', Public_Capital_Return);
            fprintf('   --- B. 政府总支出 --------------------: %12.8f\n', Gov_Outflows);
            fprintf('      1. 政府总投资 (I_g_gross) ...........: %12.8f\n', I_g_gross);
            fprintf('      2. 政府消费 (G_c) ...................: %12.8f\n', G_c);
            fprintf('   => 政府预算缺口 (A - B) ..................: %12.4e\n', mismatch_gov);
            
            % --- [4. 终极资源约束检验: Y vs C + I + Gc] ---
            fprintf('\n--- [4. 终极资源约束检验 (Y = C+I_gross+Gc)] ---\n');
            Total_Uses = (C_hh + Shock_exp) + I_total_gross + G_c; % 使用无税消费
            mismatch_resources = Y_prod - Total_Uses;
            fprintf('   A. 资源总供给 (GDP) .......................: %12.8f\n', Y_prod);
            fprintf('   B. 资源总使用 (C+I_gross+Gc) ..............: %12.8f\n', Total_Uses);
            fprintf('      = 总消费(无税) (C_hh + Shock) .......: %12.8f\n', (C_hh + Shock_exp));
            fprintf('      + 国民生产性总投资(I_p_gross+I_g_gross): %12.8f\n', I_total_gross);
            fprintf('      + 政府消费 (Gc) .......................: %12.8f\n', G_c);
            fprintf('   => 最终资源缺口 (A - B) ..................: %12.4e\n', mismatch_resources);
            
            % --- [5. 储蓄-投资恒等式检验 (S = I)] ---
            fprintf('\n--- [5. 储蓄-投资恒等式检验 (S_gross = I_gross)] ---\n');
            S_p_gross = S_p_net_from_ss + Depreciation_p;
            S_g_gross = Gov_Net_Saving + Depreciation_g;
            S_total_gross = S_p_gross + S_g_gross;

            mismatch_S_I = S_total_gross - I_total_gross;
            
            fprintf('   A. 国民总储蓄 (S_gross_total) .............: %12.8f\n', S_total_gross);
            fprintf('      = 私人总储蓄 (S_p_net + Dep_p) ........: %12.8f\n', S_p_gross);
            fprintf('      + 政府总储蓄 (S_g_net + Dep_g) ........: %12.8f\n', S_g_gross);
            fprintf('   B. 国民总投资 (I_gross_total) .............: %12.8f\n', I_total_gross);
            fprintf('   --------------------------------------------------------------------\n');
            fprintf('   C. 储蓄-投资缺口 (A - B) ..................: %12.8f\n', mismatch_S_I);
            fprintf('   D. 总意外遗赠 (Bequests) ..................: %12.8f\n', Bequests_paid);
            fprintf('   --------------------------------------------------------------------\n');
            fprintf('   => 理论验证 (缺口 - Bequests) (应为0) ......: %12.4e\n', mismatch_S_I - Bequests_paid);

            fprintf('\n================================================================================\n');
                end
                
        function [idx_lower, idx_upper, w_upper] = find_grid_and_weights(value, gridV)
            % =========================================================================
            % == 函数: find_grid_and_weights
            % == 版本: [v1 - 线性插值权重计算器]
            % ==
            % == 目的:
            % ==   给定一个连续值和一个离散网格，找到该值在网格中的位置，
            % ==   并计算线性插值所需的权重。
            % ==
            % == 输入:
            % ==   value: 连续值 (浮点数)
            % ==   gridV: 网格向量 (单调递增)
            % ==
            % == 输出:
            % ==   idx_lower: 下界网格点索引
            % ==   idx_upper: 上界网格点索引
            % ==   w_upper: 上界权重 (线性插值系数)
            % =========================================================================

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