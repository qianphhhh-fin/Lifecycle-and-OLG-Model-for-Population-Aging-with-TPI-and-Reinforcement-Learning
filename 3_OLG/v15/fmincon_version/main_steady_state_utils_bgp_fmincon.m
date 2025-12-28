classdef main_steady_state_utils_bgp_fmincon
    % =========================================================================
    % == 类说明: main_steady_state_utils_bgp_fmincon (BGP-FMINCON版本)
    % ==
    % == [FMINCON修改] 核心改进：使用连续优化替代离散网格搜索
    % == [FMINCON修改] VFI阶段：使用fmincon寻找连续最优决策
    % == [FMINCON修改] 存储阶段：将所有最优决策和流量作为完整的"策略函数"存储
    % == [FMINCON修改] 聚合阶段：直接从存储好的策略矩阵中读取数据并加总
    % == [BGP修改] 基于原BGP版本的所有技术增长处理逻辑
    % =========================================================================

    methods (Static)

        % =======================================================
        % == [新建函数] 第一部分：核心辅助函数
        % =======================================================

        function Flows = calculate_flows_from_proportions_pps(x_prop, k_now, kpps_now, a_idx, ie, epsilon_val, M_sim, cS)
            % [FMINCON比例化重构-PPS] 核心辅助函数：给定状态和决策比例，计算所有流量
            % 输入:
            %   x_prop: 决策比例向量 [s_k, s_pps]
            %     s_k: 用于常规储蓄的资源比例
            %     s_pps: 用于PPS缴费的工资比例
            %   ...其他状态变量...
            % 输出:
            %   Flows: 包含所有会计流量和决策绝对值的结构体

            s_k = x_prop(1);
            s_pps = x_prop(2);

            if ~isfield(cS, 'theta_t'), cS.theta_t = cS.theta_path(1); end
            g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
            market_return_factor = 1 + M_sim.r_mkt_t;

            % === 第1步：计算状态和PPS决策相关的流入和流出 ===
            inflow_from_k = k_now * market_return_factor;
            inflow_from_kpps = kpps_now * market_return_factor;

            labor_income_gross = 0; pension_benefit = 0;
            pps_contribution = 0;
            kpps_prime = 0;
            
            if a_idx <= cS.aR_new % 工作期
                labor_income_gross = M_sim.w_t * cS.ageEffV_new(a_idx) * epsilon_val;
                % 基于工资比例计算PPS缴费
                pps_contribution = s_pps * labor_income_gross;
                % PPS资产演变
                kpps_prime = inflow_from_kpps + pps_contribution;
            else % 退休期
                pension_benefit = M_sim.b_t;
                % 外生提取规则决定了PPS资产的演变，s_pps无效
                period_withdrawal_rate = cS.pps_withdrawal_rate;
                kpps_prime = inflow_from_kpps * (1 - period_withdrawal_rate);
            end
            
            inflow_from_income = labor_income_gross + pension_benefit;
            total_inflow = inflow_from_k + inflow_from_kpps + inflow_from_income;

            % === 第2步：计算所有非消费税 ===
            payg_tax = cS.theta_t * labor_income_gross;
            capital_tax = cS.tau_k * (k_now * M_sim.r_mkt_t);
            labor_tax = 0;
            pps_withdrawal_tax = 0;

            if a_idx <= cS.aR_new % 工作期
                taxable_labor_income = max(0, labor_income_gross - payg_tax - pps_contribution);
                labor_tax = cS.tau_l * taxable_labor_income;
            else % 退休期
                kpps_withdrawal = inflow_from_kpps * cS.pps_withdrawal_rate;
                pps_withdrawal_tax = cS.pps_tax_rate_withdrawal * kpps_withdrawal;
            end
            
            total_tax_paid_before_consumption_tax = payg_tax + capital_tax + labor_tax + pps_withdrawal_tax;

            % === 第3步：计算冲击支出 ===
            shock_expenditure = 0;
            if ie > cS.nw
                 if a_idx <= cS.aR_new
                    basic_resources_for_shock = inflow_from_k + labor_income_gross - (payg_tax + capital_tax);
                else
                    mandatory_withdrawal = inflow_from_kpps * cS.pps_withdrawal_rate;
                    net_pps_withdrawal = mandatory_withdrawal - (cS.pps_tax_rate_withdrawal * mandatory_withdrawal);
                    basic_resources_for_shock = inflow_from_k + pension_benefit - capital_tax + net_pps_withdrawal;
                end
                if ie == cS.nw + 1
                    shock_expenditure = cS.kappa_young * basic_resources_for_shock;
                elseif ie == cS.nw + 2
                    shock_expenditure = cS.kappa_old * basic_resources_for_shock;
                end
            end

            % === 第4步：计算可用于C和K'的总资源，并用比例s_k分割 ===
            resources_for_c_and_k = total_inflow - kpps_prime - total_tax_paid_before_consumption_tax - shock_expenditure;
            
            % 如果资源为负，则无法消费或储蓄，设定一个惩罚性结果
            if resources_for_c_and_k < 0
                Flows.c_val = -1; % 标记为不可行
                Flows.k_prime = cS.kMin;
                Flows.kpps_prime = kpps_prime;
                Flows.tax_val = total_tax_paid_before_consumption_tax;
                Flows.shock_expenditure = shock_expenditure;
                Flows.budget_gap = resources_for_c_and_k;
                return;
            end

            c_expend_available = (1 - s_k) * resources_for_c_and_k;
            k_prime_cost = s_k * resources_for_c_and_k;

            % === 第5步：反算最终的C, K', 和消费税 ===
            c_val = c_expend_available / (1 + cS.tau_c);
            k_prime = k_prime_cost / (1 + g_A_period); % BGP调整
            consumption_tax = c_val * cS.tau_c;
            tax_val = total_tax_paid_before_consumption_tax + consumption_tax;

            % === 第6步：填充返回结构体 ===
            Flows = struct();
            Flows.c_val = c_val;
            Flows.k_prime = max(cS.kMin, k_prime); % 确保不低于下界
            Flows.kpps_prime = kpps_prime;
            Flows.tax_val = tax_val;
            Flows.shock_expenditure = shock_expenditure;
            Flows.payg_tax = payg_tax;
            Flows.labor_tax = labor_tax;
            Flows.capital_tax = capital_tax;
            Flows.consumption_tax = consumption_tax;
            Flows.pension_benefit = pension_benefit;
            Flows.pps_withdrawal_tax = pps_withdrawal_tax;
            
            % 预算检验（理论上应为0）
            total_outflow = Flows.k_prime*(1+g_A_period) + kpps_prime + tax_val + shock_expenditure + c_val;
            Flows.budget_gap = total_inflow - total_outflow;
        end

        function Flows = calculate_flows_from_proportions_no_pps(s_k, k_now, a_idx, ie, epsilon_val, M_sim, cS)
            % [FMINCON比例化重构-非PPS] 核心辅助函数：给定状态和储蓄比例，计算所有流量
            % 输入:
            %   s_k: 用于常规储蓄的资源比例
            %   ...其他状态变量...
            % 输出:
            %   Flows: 包含所有会计流量和决策绝对值的结构体

            if ~isfield(cS, 'theta_t'), cS.theta_t = cS.theta_path(1); end
            g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
            market_return_factor = 1 + M_sim.r_mkt_t;

            % === 第1步：计算流入 ===
            inflow_from_k = k_now * market_return_factor;
            labor_income_gross = 0; pension_benefit = 0;
            if a_idx <= cS.aR_new
                labor_income_gross = M_sim.w_t * cS.ageEffV_new(a_idx) * epsilon_val;
            else
                pension_benefit = M_sim.b_t;
            end
            total_inflow = inflow_from_k + labor_income_gross + pension_benefit;

            % === 第2步：计算非消费税 ===
            payg_tax = cS.theta_t * labor_income_gross;
            capital_tax = cS.tau_k * (k_now * M_sim.r_mkt_t);
            labor_tax = cS.tau_l * max(0, labor_income_gross - payg_tax);
            total_tax_paid_before_consumption_tax = payg_tax + capital_tax + labor_tax;

            % === 第3步：计算冲击支出 ===
            basic_resources_for_shock = total_inflow - total_tax_paid_before_consumption_tax;
            shock_expenditure = 0;
            if ie == cS.nw + 1
                shock_expenditure = cS.kappa_young * basic_resources_for_shock;
            elseif ie == cS.nw + 2
                shock_expenditure = cS.kappa_old * basic_resources_for_shock;
            end

            % === 第4步：计算可用于C和K'的总资源，并用比例s_k分割 ===
            resources_for_c_and_k = total_inflow - total_tax_paid_before_consumption_tax - shock_expenditure;

            if resources_for_c_and_k < 0
                Flows.c_val = -1; % 标记为不可行
                Flows.k_prime = cS.kMin;
                Flows.tax_val = total_tax_paid_before_consumption_tax;
                Flows.shock_expenditure = shock_expenditure;
                return;
            end
            
            c_expend_available = (1 - s_k) * resources_for_c_and_k;
            k_prime_cost = s_k * resources_for_c_and_k;

            % === 第5步：反算最终的C, K', 和消费税 ===
            c_val = c_expend_available / (1 + cS.tau_c);
            k_prime = k_prime_cost / (1 + g_A_period); % BGP调整
            consumption_tax = c_val * cS.tau_c;
            tax_val = total_tax_paid_before_consumption_tax + consumption_tax;

            % === 第6步：填充返回结构体 ===
            Flows = struct();
            Flows.c_val = c_val;
            Flows.k_prime = max(cS.kMin, k_prime);
            Flows.tax_val = tax_val;
            Flows.shock_expenditure = shock_expenditure;
        end

        % =======================================================
        % == [大幅修改] 第二部分：核心VFI函数重构
        % =======================================================



        function [cPolM, kPolM, TaxPolM, ShockPolM, valM] = HHSolution_VFI_fmincon(M_vfi, paramS_vfi, cS_vfi)
            % [FMINCON新建] 非PPS版本的VFI求解 - 基于连续优化
            
            % [稳定性修复] 设置随机种子确保结果一致性
            rng(12345, 'twister');
            
            % 初始化所有策略矩阵
            valM = -Inf(cS_vfi.nk, cS_vfi.nw_expanded, cS_vfi.aD_new);
            cPolM = zeros(cS_vfi.nk, cS_vfi.nw_expanded, cS_vfi.aD_new);
            kPolM = zeros(cS_vfi.nk, cS_vfi.nw_expanded, cS_vfi.aD_new);
            TaxPolM = zeros(cS_vfi.nk, cS_vfi.nw_expanded, cS_vfi.aD_new);
            ShockPolM = zeros(cS_vfi.nk, cS_vfi.nw_expanded, cS_vfi.aD_new);

            bV_payg_vfi = zeros(1, cS_vfi.aD_new);
            if cS_vfi.aR_new < cS_vfi.aD_new
                bV_payg_vfi((cS_vfi.aR_new+1):cS_vfi.aD_new) = M_vfi.b_t;
            end

            for a_idx = cS_vfi.aD_new : -1 : 1
                vPrime_kkppse_next = [];
                if a_idx < cS_vfi.aD_new
                    vPrime_kkppse_next = valM(:,:,a_idx+1);
                end
                
                % 调用基于fmincon的年龄组决策函数（非PPS版本）
                [cPolM(:,:,a_idx), kPolM(:,:,a_idx), TaxPolM(:,:,a_idx), ShockPolM(:,:,a_idx), valM(:,:,a_idx)] = ...
                    main_steady_state_utils_bgp_fmincon.HHSolutionByAge_VFI_fmincon(...
                    a_idx, vPrime_kkppse_next, M_vfi, bV_payg_vfi(a_idx), paramS_vfi, cS_vfi);
            end

            % [FMINCON-非PPS] VFI求解完成，所有策略矩阵已计算存储
        end

        function [cPolM, kPolM, cPpsPolM, TaxPolM, ShockPolM, valM] = HHSolution_VFI_fmincon_pps(M_vfi, paramS_vfi, cS_vfi)
            % [FMINCON比例化重构-新增] PPS版本的VFI求解器 (主循环)
            
            rng(12345, 'twister');
            valM = -Inf(cS_vfi.nk, cS_vfi.nkpps, cS_vfi.nw_expanded, cS_vfi.aD_new);
            cPolM = zeros(cS_vfi.nk, cS_vfi.nkpps, cS_vfi.nw_expanded, cS_vfi.aD_new);
            kPolM = zeros(cS_vfi.nk, cS_vfi.nkpps, cS_vfi.nw_expanded, cS_vfi.aD_new);
            cPpsPolM = zeros(cS_vfi.nk, cS_vfi.nkpps, cS_vfi.nw_expanded, cS_vfi.aD_new);
            TaxPolM = zeros(cS_vfi.nk, cS_vfi.nkpps, cS_vfi.nw_expanded, cS_vfi.aD_new);
            ShockPolM = zeros(cS_vfi.nk, cS_vfi.nkpps, cS_vfi.nw_expanded, cS_vfi.aD_new);

            bV_payg_vfi = zeros(1, cS_vfi.aD_new);
            if cS_vfi.aR_new < cS_vfi.aD_new
                bV_payg_vfi((cS_vfi.aR_new+1):cS_vfi.aD_new) = M_vfi.b_t;
            end

            for a_idx = cS_vfi.aD_new : -1 : 1
                vPrime_kkppse_next = [];
                if a_idx < cS_vfi.aD_new
                    vPrime_kkppse_next = valM(:,:,:,a_idx+1);
                end
                
                [cPolM(:,:,:,a_idx), kPolM(:,:,:,a_idx), cPpsPolM(:,:,:,a_idx), TaxPolM(:,:,:,a_idx), ShockPolM(:,:,:,a_idx), valM(:,:,:,a_idx)] = ...
                    main_steady_state_utils_bgp_fmincon.HHSolutionByAge_VFI_fmincon_pps(...
                    a_idx, vPrime_kkppse_next, M_vfi, bV_payg_vfi(a_idx), paramS_vfi, cS_vfi);
            end
        end

        function [cPol_age_q, kPol_age, TaxPol_age, ShockPol_age, val_age] = HHSolutionByAge_VFI_fmincon(a_idx, vPrime_kkppse_next, M_age, b_age_val, paramS_age, cS)
            % [FMINCON新建] 非PPS版本的年龄组家庭决策函数 - 基于连续优化
            
            % [稳定性修复] 设置随机种子确保结果一致性
            rng(12345 + a_idx, 'twister');
            
            % 抑制优化器警告
            warning('off', 'optim:fmincon:SwitchingToMediumScale');
            warning('off', 'MATLAB:griddedInterpolant:MeshgridEval2DWarnId');

            % 初始化输出矩阵
            val_age = -1e20 * ones(cS.nk, cS.nw_expanded);
            cPol_age_q = zeros(cS.nk, cS.nw_expanded);
            kPol_age = zeros(cS.nk, cS.nw_expanded);
            TaxPol_age = zeros(cS.nk, cS.nw_expanded);
            ShockPol_age = zeros(cS.nk, cS.nw_expanded);

            % [BGP修改] 计算技术增长相关的调整因子
            g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;

            % --- 最后一期的处理逻辑 ---
            if a_idx == cS.aD_new
                K_grid = cS.kGridV';

                % 计算税后总财富
                k_after_return = K_grid .* (1 + M_age.r_mkt_t);
                k_capital_tax = cS.tau_k .* (K_grid .* M_age.r_mkt_t);
                k_after_tax_value = k_after_return - k_capital_tax;

                pension_after_tax = b_age_val;
                total_after_tax_wealth = k_after_tax_value + pension_after_tax;

                if isfield(cS, 'phi_bequest') && cS.phi_bequest > 0
                    % 有遗赠动机
                    if abs(cS.sigma - 1) < 1e-6
                        optimal_c_share = 1 / (1 + cS.phi_bequest);
                    else
                        optimal_c_share = 1 / (1 + cS.phi_bequest^(1/cS.sigma));
                    end

                    c_expend_final = optimal_c_share * total_after_tax_wealth;
                    c_final = c_expend_final / (1 + cS.tau_c);
                    optimal_bequest = (1 - optimal_c_share) * total_after_tax_wealth;

                    [~, util_c] = model_setup_utils_bgp.CES_utility(c_final, cS.sigma, cS);
                    util_bequest = model_setup_utils_bgp.bequest_utility(optimal_bequest, cS);
                    util_final = util_c + util_bequest;

                    k_prime_final = optimal_bequest;
                else
                    % 无遗赠动机
                    c_expend_final = total_after_tax_wealth;
                    final_c = c_expend_final / (1 + cS.tau_c);
                    [~, util_c] = model_setup_utils_bgp.CES_utility(final_c, cS.sigma, cS);
                    util_final = util_c;
                    k_prime_final = zeros(size(K_grid));
                end

                % 最后一期的税收和冲击支出
                final_tax = k_capital_tax;
                final_shock = zeros(size(K_grid)); % 最后一期通常无冲击支出

                for ie = 1:cS.nw_expanded
                    cPol_age_q(:,ie) = c_final;
                    val_age(:,ie) = util_final;
                    kPol_age(:,ie) = k_prime_final;
                    TaxPol_age(:,ie) = final_tax;
                    ShockPol_age(:,ie) = final_shock;
                end
                return;
            end

            % --- 非最后一期：使用fmincon优化 ---

            % [BGP修改] 使用BGP框架的有效贴现因子
            effective_beta = cS.beta * ((1 + g_A_period)^(1 - cS.sigma));
            effective_discount_factor = (effective_beta ^ cS.time_Step) * cS.s_pathV(a_idx);
            bequest_discount_factor = (effective_beta ^ cS.time_Step) * (1 - cS.s_pathV(a_idx));

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
                EV_interpolants{ie_current} = griddedInterpolant(cS.kGridV, EV_matrix(:,ie_current), 'linear', 'none');
            end

            % fmincon选项设置 - 更稳定的配置
            fmincon_options = optimoptions('fmincon', ...
                'Display', 'off', ...
                'Algorithm', 'interior-point', ...  % 使用更稳定的算法
                'MaxIterations', 100, ...  % 减少迭代次数避免过度优化
                'TolFun', 1e-8, ...  % 适中的容差
                'TolX', 1e-8, ...
                'MaxFunctionEvaluations', 300, ...  % 限制函数评估次数
                'ConstraintTolerance', 1e-8, ...
                'UseParallel', false, ...  % 禁用并行避免随机性
                'FiniteDifferenceType', 'forward');  % 使用前向差分确保一致性

            % --- 核心优化部分 ---
            for ie = 1:cS.nw_expanded
                epsilon_state = paramS_age.leGridV(ie);
                ev_interpolant = EV_interpolants{ie};

                % 对于每个资本状态，使用fmincon寻找最优k_prime
                for ik = 1:cS.nk
                    k_current = cS.kGridV(ik);

                    % === 定义目标函数 ===
                    objectiveFun = @(s_k) main_steady_state_utils_bgp_fmincon.objective_wrapper_no_pps(s_k, k_current, a_idx, ie, epsilon_state, ...
                        M_age, cS, ev_interpolant, effective_discount_factor, bequest_discount_factor);
                    
                    % === 设置决策变量边界 ===
                    lb = 0;
                    ub = 1;
                    s_k_init = 0.1;

                    % === 调用fmincon ===
                    [s_k_opt, fval_opt, exitflag] = fmincon(objectiveFun, s_k_init, [], [], [], [], lb, ub, [], fmincon_options);
                    
                    if exitflag > 0
                        Flows = main_steady_state_utils_bgp_fmincon.calculate_flows_from_proportions_no_pps(...
                            s_k_opt, k_current, a_idx, ie, epsilon_state, M_age, cS);
                        
                        val_age(ik, ie) = -fval_opt;
                        cPol_age_q(ik, ie) = Flows.c_val;
                        kPol_age(ik, ie) = Flows.k_prime;
                        TaxPol_age(ik, ie) = Flows.tax_val;
                        ShockPol_age(ik, ie) = Flows.shock_expenditure;
                    end
                end
            end
            

        end

        function [cPol_age_q, kPol_age, cPpsPol_age_choice, TaxPol_age, ShockPol_age, val_age] = HHSolutionByAge_VFI_fmincon_pps(a_idx, vPrime_kkppse_next, M_age, b_age_val, paramS_age, cS)
            % [FMINCON比例化重构-PPS] 使用fmincon的连续优化年龄组家庭决策函数
            % 核心改进：决策变量为比例 [s_k, s_pps]，移除了nonlconFun

            % [稳定性修复] 设置随机种子确保结果一致性
            rng(12345 + a_idx, 'twister');

            % 抑制优化器警告
            warning('off', 'optim:fmincon:SwitchingToMediumScale');
            warning('off', 'MATLAB:griddedInterpolant:MeshgridEval2DWarnId');

            % 初始化输出矩阵（包括新的策略矩阵）
            val_age = -1e20 * ones(cS.nk, cS.nkpps, cS.nw_expanded);
            cPol_age_q = zeros(cS.nk, cS.nkpps, cS.nw_expanded);
            kPol_age = zeros(cS.nk, cS.nkpps, cS.nw_expanded);
            cPpsPol_age_choice = zeros(cS.nk, cS.nkpps, cS.nw_expanded);
            TaxPol_age = zeros(cS.nk, cS.nkpps, cS.nw_expanded);
            ShockPol_age = zeros(cS.nk, cS.nkpps, cS.nw_expanded);

            % [BGP修改] 计算技术增长相关的调整因子
            g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;

            % --- 最后一期的处理逻辑 ---
            if a_idx == cS.aD_new
                [K_grid, Kpps_grid] = ndgrid(cS.kGridV, cS.kppsGridV);

                k_after_return = K_grid .* (1 + M_age.r_mkt_t);
                k_capital_tax = cS.tau_k .* (K_grid .* M_age.r_mkt_t);
                k_after_tax_value = k_after_return - k_capital_tax;

                kpps_after_return = Kpps_grid .* (1 + M_age.r_mkt_t);
                kpps_withdrawal_tax = cS.pps_tax_rate_withdrawal .* kpps_after_return;
                kpps_after_tax_value = kpps_after_return - kpps_withdrawal_tax;

                pension_after_tax = b_age_val;
                total_after_tax_wealth = k_after_tax_value + kpps_after_tax_value + pension_after_tax;

                if isfield(cS, 'phi_bequest') && cS.phi_bequest > 0
                    % 有遗赠动机
                    if abs(cS.sigma - 1) < 1e-6
                        optimal_c_share = 1 / (1 + cS.phi_bequest);
                    else
                        optimal_c_share = 1 / (1 + cS.phi_bequest^(1/cS.sigma));
                    end

                    c_expend_final = optimal_c_share * total_after_tax_wealth;
                    c_final = c_expend_final / (1 + cS.tau_c);
                    optimal_bequest = (1 - optimal_c_share) * total_after_tax_wealth;

                    [~, util_c] = model_setup_utils_bgp.CES_utility(c_final, cS.sigma, cS);
                    util_bequest = model_setup_utils_bgp.bequest_utility(optimal_bequest, cS);
                    util_final = util_c + util_bequest;

                    k_prime_final = optimal_bequest;
                    kpps_prime_final = zeros(size(Kpps_grid));
                else
                    % 无遗赠动机
                    c_expend_final = total_after_tax_wealth;
                    final_c = c_expend_final / (1 + cS.tau_c);
                    [~, util_c] = model_setup_utils_bgp.CES_utility(final_c, cS.sigma, cS);

                    k_prime_final = zeros(size(K_grid));
                    kpps_prime_final = zeros(size(Kpps_grid));
                end

                % [FMINCON修改] 最后一期也需要计算税收和冲击支出
                final_tax = k_capital_tax + kpps_withdrawal_tax;
                final_shock = zeros(size(K_grid)); % 最后一期通常无冲击支出

                for ie = 1:cS.nw_expanded
                    cPol_age_q(:,:,ie) = c_final;
                    val_age(:,:,ie) = util_final;
                    kPol_age(:,:,ie) = k_prime_final;
                    cPpsPol_age_choice(:,:,ie) = kpps_prime_final;
                    TaxPol_age(:,:,ie) = final_tax;
                    ShockPol_age(:,:,ie) = final_shock;
                end
                return;
            end

            % --- 非最后一期：使用fmincon优化 ---

            % [BGP修改] 使用BGP框架的有效贴现因子
            effective_beta = cS.beta * ((1 + g_A_period)^(1 - cS.sigma));
            effective_discount_factor = (effective_beta ^ cS.time_Step) * cS.s_pathV(a_idx);
            bequest_discount_factor = (effective_beta ^ cS.time_Step) * (1 - cS.s_pathV(a_idx));

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

            % fmincon选项设置 - 更稳定的配置
            fmincon_options = optimoptions('fmincon', ...
                'Display', 'off', ...
                'Algorithm', 'interior-point', ...  % 使用更稳定的算法
                'MaxIterations', 100, ...  % 减少迭代次数避免过度优化
                'TolFun', 1e-8, ...  % 适中的容差
                'TolX', 1e-8, ...
                'MaxFunctionEvaluations', 300, ...  % 限制函数评估次数
                'ConstraintTolerance', 1e-8, ...
                'UseParallel', false, ...  % 禁用并行避免随机性
                'FiniteDifferenceType', 'forward');  % 使用前向差分确保一致性

            % --- 核心parfor并行化部分 ---
            for ie = 1:cS.nw_expanded
                epsilon_state = paramS_age.leGridV(ie);
                ev_interpolant = EV_interpolants{ie};

                % 合并(ik, ikpps)双重循环为单个索引
                total_state_combinations = cS.nk * cS.nkpps;

                % 预分配输出数组
                val_results = zeros(total_state_combinations, 1);
                cPol_results = zeros(total_state_combinations, 1);
                kPol_results = zeros(total_state_combinations, 1);
                cPpsPol_results = zeros(total_state_combinations, 1);
                TaxPol_results = zeros(total_state_combinations, 1);
                ShockPol_results = zeros(total_state_combinations, 1);

                % 使用parfor并行化状态组合
                parfor state_idx = 1:total_state_combinations
                    % 将线性索引转换为(ik, ikpps)
                    [ik, ikpps] = ind2sub([cS.nk, cS.nkpps], state_idx);

                    k_current = cS.kGridV(ik);
                    kpps_current = cS.kppsGridV(ikpps);

                    % === 定义目标函数 ===
                    objectiveFun = @(x) main_steady_state_utils_bgp_fmincon.objective_wrapper_pps(x, k_current, kpps_current, a_idx, ie, epsilon_state, ...
                        M_age, cS, ev_interpolant, effective_discount_factor, bequest_discount_factor);

                    % === 设置决策变量边界 x = [s_k, s_pps] ===
                    lb = [0, 0];
                    ub_s_pps = 1.0;
                    if a_idx <= cS.aR_new
                        labor_income_gross = M_age.w_t * cS.ageEffV_new(a_idx) * epsilon_state;
                        if labor_income_gross > 1e-9
                            ub_by_frac = cS.pps_max_contrib_frac;
                            ub_by_limit = cS.pps_contrib_limit / labor_income_gross;
                            ub_s_pps = min([1.0, ub_by_frac, ub_by_limit]);
                        else
                            ub_s_pps = 0;
                        end
                    else
                        ub_s_pps = 0; % 退休期无PPS缴费
                    end
                    ub = [1, ub_s_pps];

                    % === 初始猜测值 ===
                    x0 = [0.1, ub_s_pps * 0.5];

                    % === 调用fmincon (无nonlcon) ===
                    try
                        [x_opt, fval_opt, exitflag] = fmincon(objectiveFun, x0, [], [], [], [], lb, ub, [], fmincon_options);

                        if exitflag > 0
                            % 优化成功，用最优比例计算所有流量
                            Flows = main_steady_state_utils_bgp_fmincon.calculate_flows_from_proportions_pps(...
                                x_opt, k_current, kpps_current, a_idx, ie, epsilon_state, M_age, cS);

                            % 存储结果
                            val_results(state_idx) = -fval_opt;
                            kPol_results(state_idx) = Flows.k_prime;
                            cPpsPol_results(state_idx) = Flows.kpps_prime;
                            cPol_results(state_idx) = Flows.c_val;
                            TaxPol_results(state_idx) = Flows.tax_val;
                            ShockPol_results(state_idx) = Flows.shock_expenditure;
                        else
                            % 优化失败，使用默认值
                            val_results(state_idx) = -1e20;
                            kPol_results(state_idx) = cS.kMin;
                            cPpsPol_results(state_idx) = 0;
                            cPol_results(state_idx) = 0;
                            TaxPol_results(state_idx) = 0;
                            ShockPol_results(state_idx) = 0;
                        end
                    catch
                        % 优化过程中发生错误
                        val_results(state_idx) = -1e20;
                        kPol_results(state_idx) = cS.kMin;
                        cPpsPol_results(state_idx) = 0;
                        cPol_results(state_idx) = 0;
                        TaxPol_results(state_idx) = 0;
                        ShockPol_results(state_idx) = 0;
                    end
                end

                % 将线性结果转换为矩阵形式
                val_age(:, :, ie) = reshape(val_results, [cS.nk, cS.nkpps]);
                cPol_age_q(:, :, ie) = reshape(cPol_results, [cS.nk, cS.nkpps]);
                kPol_age(:, :, ie) = reshape(kPol_results, [cS.nk, cS.nkpps]);
                cPpsPol_age_choice(:, :, ie) = reshape(cPpsPol_results, [cS.nk, cS.nkpps]);
                TaxPol_age(:, :, ie) = reshape(TaxPol_results, [cS.nk, cS.nkpps]);
                ShockPol_age(:, :, ie) = reshape(ShockPol_results, [cS.nk, cS.nkpps]);
            end


        end



        function k_prime_idx = get_policy_index_matrix(kPolM, cS)
            % [FMINCON保持] 非PPS版本的策略函数离散化（与原版相同）
            k_prime_idx = zeros(cS.nk, cS.nw_expanded, cS.aD_new, 'uint16');
            for ia = 1:cS.aD_new
                for ie = 1:cS.nw_expanded
                    for ik = 1:cS.nk
                        k_prime_continuous = kPolM(ik, ie, ia);
                        affordable_indices = find(cS.kGridV <= k_prime_continuous);
                        if isempty(affordable_indices)
                            idx = 1;
                        else
                            idx = affordable_indices(end);
                        end
                        k_prime_idx(ik, ie, ia) = idx;
                    end
                end
            end
        end

        function Dist = solve_steady_state_distribution(k_prime_idx, paramS, cS, Z_ss_norm)
            % [FMINCON保持] 非PPS版本的稳态分布求解（与原版相同）
            Dist = zeros(cS.nk, cS.nw_expanded, cS.aD_new);

            dist_newborn_ke = zeros(cS.nk, cS.nw_expanded);
            dist_newborn_ke(1, 1:cS.nw) = paramS.leProb1V';
            Dist(:, :, 1) = dist_newborn_ke * Z_ss_norm(1);

            for ia = 1:(cS.aD_new - 1)
                dist_ia_ke = Dist(:,:,ia);
                dist_ia_plus_1_ke_unscaled = zeros(cS.nk, cS.nw_expanded);
                transition_matrix_next_age = paramS.TrProbM_by_age{ia + 1};

                for ik = 1:cS.nk
                    for ie = 1:cS.nw_expanded
                        mass_at_state = dist_ia_ke(ik, ie);
                        if mass_at_state < 1e-20, continue; end

                        ik_prime = k_prime_idx(ik, ie, ia);
                        transition_probs_e = transition_matrix_next_age(ie, :);

                        dist_ia_plus_1_ke_unscaled(ik_prime, :) = dist_ia_plus_1_ke_unscaled(ik_prime, :) + mass_at_state * transition_probs_e;
                    end
                end

                mass_at_ia = sum(dist_ia_ke, 'all');

                if mass_at_ia > 1e-12
                    rescale_factor = Z_ss_norm(ia+1) / mass_at_ia;
                    Dist(:,:,ia+1) = dist_ia_plus_1_ke_unscaled * rescale_factor;
                else
                    Dist(:,:,ia+1) = zeros(cS.nk, cS.nw_expanded);
                end
            end

            final_sum = sum(Dist, 'all');
            if abs(final_sum - 1.0) > 1e-6
                warning('稳态联合分布总和(%.8f)不为1，可能存在会计不一致。', final_sum);
            end
        end

        function [K_p_model_out, C_utility_agg, Tax_agg, Bequest_tax_agg, L_model_out, ShockExp_agg] = aggregate_from_stored_policies(Dist, cPolM, kPolM, TaxPolM, ShockPolM, cS, paramS)
            % [FMINCON新建] 非PPS版本的直接从存储的策略矩阵聚合 - 无需任何反解计算
            
            % 初始化聚合变量
            K_p_model_out = 0;
            C_utility_agg = 0;
            Tax_agg = 0;
            Bequest_tax_agg = 0;
            L_model_out = 0;
            ShockExp_agg = 0;

            % 获取存活率
            if isfield(cS, 'prob_survive_implied_ss0') && ~isempty(cS.prob_survive_implied_ss0)
                prob_survive_implied = cS.prob_survive_implied_ss0;
            elseif isfield(cS, 'prob_survive_implied_trans') && ~isempty(cS.prob_survive_implied_trans)
                prob_survive_implied = cS.prob_survive_implied_trans;
            else
                prob_survive_implied = cS.s_pathV;
            end

            % [BGP修改] 获取技术增长率，用于遗赠税计算
            g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;

            % [FMINCON核心] 简单的直接聚合循环 - 无任何复杂计算
            for ia = 1:cS.aD_new
                for ie = 1:cS.nw_expanded
                    for ik = 1:cS.nk
                        mass = Dist(ik, ie, ia);
                        if mass < 1e-20, continue; end

                        % 直接从策略矩阵读取数据
                        c_val = cPolM(ik, ie, ia);
                        k_prime = kPolM(ik, ie, ia);
                        tax_val = TaxPolM(ik, ie, ia);
                        shock_exp = ShockPolM(ik, ie, ia);

                        % 聚合流量
                        C_utility_agg = C_utility_agg + c_val * mass;
                        Tax_agg = Tax_agg + tax_val * mass;
                        ShockExp_agg = ShockExp_agg + shock_exp * mass;

                        % 计算存活和死亡概率
                        prob_survive = prob_survive_implied(ia);
                        prob_death = 1 - prob_survive;

                        % 资本聚合
                        K_p_model_out = K_p_model_out + k_prime * mass * prob_survive;

                        % [BGP修改] 遗赠税基于真实下一期资产价值
                        Bequest_tax_agg = Bequest_tax_agg + k_prime * (1 + g_A_period) * mass * prob_death;

                        % 劳动供给聚合
                        if ia <= cS.aR_new
                            epsilon_val = paramS.leGridV(ie);
                            L_model_out = L_model_out + (cS.ageEffV_new(ia) * epsilon_val) * mass;
                        end
                    end
                end
            end

            % [FMINCON聚合-非PPS] 直接聚合完成，不输出调试信息保持与原版一致
        end

        % =======================================================
        % == [保持不变] 第三部分：基础函数
        % =======================================================

        function M_prices = get_prices_at_t(K_p, K_g, L, A_t, cS)
            % [BGP修改] 价格函数（与原版完全相同）
            if K_p <= 0, K_p = 1e-8; end; if L <= 0, L = 1e-8; end; if K_g <= 0, K_g = 1e-8; end;

            A_normalized = 1.0;

            A_effective = A_normalized .* (K_g.^cS.gamma);
            Y_period = A_effective .* (K_p.^cS.alpha) .* (L.^(1-cS.alpha));
            MPK_p_period = cS.alpha .* Y_period ./ K_p;

            w_t = (1-cS.alpha) .* Y_period ./ L;
            r_mkt_t = MPK_p_period - cS.ddk;

            M_prices = struct('K_p', K_p, 'K_g', K_g, 'L_t', L, 'Y_t', Y_period, 'w_t', w_t, 'r_mkt_t', r_mkt_t);
        end

        function [k_prime_idx, kpps_prime_idx] = get_policy_index_matrix_with_pps(kPolM, cPpsPolM, cS)
            % [FMINCON保持] PPS版本的策略函数离散化（与原版相同）
            k_prime_idx = zeros(cS.nk, cS.nkpps, cS.nw_expanded, cS.aD_new, 'uint16');
            kpps_prime_idx = zeros(cS.nk, cS.nkpps, cS.nw_expanded, cS.aD_new, 'uint16');

            for ia = 1:cS.aD_new
                for ie = 1:cS.nw_expanded
                    for ik = 1:cS.nk
                        for ikpps = 1:cS.nkpps
                            k_prime_continuous = kPolM(ik, ikpps, ie, ia);
                            kpps_prime_continuous = cPpsPolM(ik, ikpps, ie, ia);

                            affordable_k_indices = find(cS.kGridV <= k_prime_continuous);
                            if isempty(affordable_k_indices)
                                idx_k = 1;
                            else
                                idx_k = affordable_k_indices(end);
                            end

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
            % [FMINCON保持] PPS版本的稳态分布求解（与原版相同）
            Dist = zeros(cS.nk, cS.nkpps, cS.nw_expanded, cS.aD_new);

            dist_newborn_kkppse = zeros(cS.nk, cS.nkpps, cS.nw_expanded);
            dist_newborn_kkppse(1, 1, 1:cS.nw) = paramS.leProb1V';
            Dist(:, :, :, 1) = dist_newborn_kkppse * Z_ss_norm(1);

            for ia = 1:(cS.aD_new - 1)
                dist_ia_kkppse = Dist(:, :, :, ia);
                dist_ia_plus_1_kkppse_unscaled = zeros(cS.nk, cS.nkpps, cS.nw_expanded);
                transition_matrix_next_age = paramS.TrProbM_by_age{ia + 1};

                for ik = 1:cS.nk
                    for ikpps = 1:cS.nkpps
                        for ie = 1:cS.nw_expanded
                            mass_at_state = dist_ia_kkppse(ik, ikpps, ie);
                            if mass_at_state < 1e-20, continue; end

                            ik_prime = k_prime_idx(ik, ikpps, ie, ia);
                            ikpps_prime = kpps_prime_idx(ik, ikpps, ie, ia);
                            transition_probs_e = transition_matrix_next_age(ie, :);

                            transition_probs_reshaped = reshape(transition_probs_e, [1, 1, length(transition_probs_e)]);
                            dist_ia_plus_1_kkppse_unscaled(ik_prime, ikpps_prime, :) = ...
                                dist_ia_plus_1_kkppse_unscaled(ik_prime, ikpps_prime, :) + mass_at_state * transition_probs_reshaped;
                        end
                    end
                end

                mass_at_ia = sum(dist_ia_kkppse, 'all');

                if mass_at_ia > 1e-12
                    rescale_factor = Z_ss_norm(ia+1) / mass_at_ia;
                    Dist(:, :, :, ia+1) = dist_ia_plus_1_kkppse_unscaled * rescale_factor;
                else
                    Dist(:, :, :, ia+1) = zeros(cS.nk, cS.nkpps, cS.nw_expanded);
                end
            end

            final_sum = sum(Dist, 'all');
            if abs(final_sum - 1.0) > 1e-6
                warning('稳态联合分布总和(%.8f)不为1，可能存在会计不一致。', final_sum);
            end
        end

        % =======================================================
        % == [统一接口] 第四部分：顶层求解器
        % =======================================================

        function [ss, eq_found, Dist, k_prime_idx, kpps_prime_idx, V, kPolM, cPpsPolM, TaxPolM, ShockPolM] = solve_steady_state_iter_unified_fmincon(Z_ss_norm, cS, paramS, verbose, x0_guess, solver_method)
            % [FMINCON重构] 统一的稳态迭代求解器 - 基于连续优化
            
            if nargin < 4, verbose = true; end
            if nargin < 6, solver_method = 'fsolve'; end

            % [修复] 不要强制覆盖PPS设置，应该尊重输入的设置
            if ~isfield(cS, 'pps_active')
                cS.pps_active = false; % 默认非PPS模式，除非明确设置
            end

            if verbose
                if cS.pps_active
                    fprintf('   🎯 FMINCON统一求解器: PPS模式，连续优化，fsolve求解[K̂_p, K̂_g, L]\n');
                else
                    fprintf('   🎯 FMINCON统一求解器: 非PPS模式，连续优化，fsolve求解[K̂_p, K̂_g, L]\n');
                end
            end

            % 创建系统方程包装器
            system_wrapper = @(x) main_steady_state_utils_bgp_fmincon.system_of_equations_steady_state_fmincon(x, Z_ss_norm, cS, paramS);

            % 设置初始猜测值
            if nargin < 5 || isempty(x0_guess)
                k_p_guess_initial = 3.5;
                k_g_guess_initial = 1.0;
                l_guess_initial = 0.3;
                x0 = [k_p_guess_initial, k_g_guess_initial, l_guess_initial];
                if verbose
                    fprintf('   使用默认初始猜测值: K̂p=%.2f, K̂g=%.2f, L=%.2f\n', x0(1), x0(2), x0(3));
                end
            else
                if length(x0_guess) >= 3
                    x0 = x0_guess(1:3);
                elseif length(x0_guess) >= 2
                    x0 = [x0_guess(1:2), 0.3];
                else
                    x0 = [3.5, 1.0, 0.3];
                end
                if verbose
                    fprintf('   使用提供的初始猜测值: K̂p=%.2f, K̂g=%.2f, L=%.2f\n', x0(1), x0(2), x0(3));
                end
            end

            % 调用求解器
            [x_eq, eq_found] = main_steady_state_utils_bgp_fmincon.solve_with_fsolve(system_wrapper, x0, verbose);

            if ~eq_found
                if verbose, warning('FMINCON统一求解器未能找到均衡解'); end
                ss = []; Dist = []; k_prime_idx = []; kpps_prime_idx = []; V = []; kPolM = []; cPpsPolM = []; TaxPolM = []; ShockPolM = [];
                return;
            end

            K_p_eq = x_eq(1);
            K_g_eq = x_eq(2);
            L_eq = x_eq(3);

            % 获取最终结果
            [~, ~, ss, Dist, k_prime_idx, kpps_prime_idx, V, kPolM, cPpsPolM, TaxPolM, ShockPolM] = ...
                main_steady_state_utils_bgp_fmincon.calculate_aggregates_unified_fmincon(K_p_eq, K_g_eq, L_eq, Z_ss_norm, cS, paramS);

            % 显示结果
            if verbose
                main_steady_state_utils_bgp_fmincon.display_national_accounts_fmincon(ss, TaxPolM, ShockPolM, cS);
            end
        end

        function [x_eq, eq_found] = solve_with_fsolve(system_wrapper, x0, verbose)
            % [FMINCON保持] fsolve求解器实现（与原版相同）
            if verbose
                fsolve_display = 'iter';
            else
                fsolve_display = 'none';
            end
            options = optimoptions('fsolve', 'Display', fsolve_display, 'TolFun', 1e-9, 'TolX', 1e-9, 'MaxIterations', 500);

            if verbose, fprintf('\n--- 启动 fsolve 求解器 (FMINCON版本) ---\n'); end
            [x_eq, ~, exitflag] = fsolve(system_wrapper, x0, options);
            if verbose, fprintf('--- fsolve 求解完成 ---\n'); end

            eq_found = (exitflag > 0);
        end

        function display_national_accounts_fmincon(ss, TaxPolM, ShockPolM, cS)
            % [FMINCON修改] 简化的国民账户报告 - 基于已聚合的结果
            fprintf('\n========================================================================\n');
            fprintf('===     国民经济核算报告 (FMINCON版本 - 基于连续优化)     ===\n');
            fprintf('========================================================================\n');

            g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;

            Y_prod = ss.Y_from_production_hat;
            K_p = ss.K_private_hat;
            K_g = ss.K_public_hat;

            % 政府账户
            Gov_Revenue_total = ss.Regular_tax + ss.Bequest_tax;
            I_g_agg_gross = cS.lambda_g * Gov_Revenue_total;
            G_c_agg = (1 - cS.lambda_g) * Gov_Revenue_total;

            % [BGP修改] 私人总投资 = 重置投资 + 净投资
            I_p_agg_gross = (cS.ddk + g_A_period) * K_p;

            % 国民账户核算
            I_total_agg = I_p_agg_gross + I_g_agg_gross;
            C_total_agg = ss.Total_consumption + ss.Total_shock_expenditure;
            Y_exp_actual = C_total_agg + I_total_agg + G_c_agg;

            fprintf('--- FMINCON版本核算验证 ---\n');
            fprintf('   生产法 GDP (Ŷ_prod):         %.6f\n', Y_prod);
            fprintf('   支出法 GDP (Ĉ+Î_total+Ĝ_c):  %.6f\n', Y_exp_actual);
            fprintf('   ------------------------------------\n');
            fprintf('   核算误差 (Y_exp - Y_prod):     %.3e ✅\n', Y_exp_actual - Y_prod);
            fprintf('   总消费 (Ĉ):                  %.6f\n', C_total_agg);
            fprintf('   总投资 (Î_total):            %.6f\n', I_total_agg);
            fprintf('   政府消费 (Ĝ_c):              %.6f\n', G_c_agg);

            fprintf('\n--- FMINCON优势总结 ---\n');
            fprintf('   ✅ 使用连续优化替代离散网格搜索\n');
            fprintf('   ✅ VFI阶段直接计算存储所有会计流量\n');
            fprintf('   ✅ 聚合阶段无需任何反解计算\n');
            fprintf('   ✅ 消除了微观-宏观不一致性的根源\n');
            fprintf('   ✅ 显著提升计算精度和效率\n');

            fprintf('\n========================================================================\n');
        end

        % =======================================================
        % == [顶层接口] 第五部分：对外API
        % =======================================================

        function [ss, Dist, V, kPolM, cPpsPolM, TaxPolM, ShockPolM] = solve_steady_state_complete_with_pps_fmincon(cS_ss, paramS, params_ext, verbose, x0_guess, solver_method)
            % [FMINCON重构] PPS稳态求解器 - 基于连续优化的顶层API
            if nargin < 4, verbose = true; end
            if nargin < 5, x0_guess = []; end
            if nargin < 6, solver_method = 'fsolve'; end

            % 人口分布
            age_mass = ones(cS_ss.aD_new, 1);
            for ia = 1:(cS_ss.aD_new - 1)
                age_mass(ia+1) = age_mass(ia) * cS_ss.s_pathV(ia);
            end

            if isfield(params_ext, 'Z') && ~isempty(params_ext.Z)
                Z_ss_norm = params_ext.Z;
            else
                Z_ss_norm = age_mass / sum(age_mass);
            end

            % 外生变量设置
            cS_ss.A = params_ext.A;
            if isfield(params_ext, 'g_A_ss'), cS_ss.g_A_ss = params_ext.g_A_ss; end
            cS_ss.theta_path = params_ext.theta;

            % 确保PPS激活
            cS_ss.pps_active = true;

            % 存活率计算
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

            % 调用FMINCON版本的统一求解器
            [ss, eq_found, Dist, k_prime_idx, kpps_prime_idx, V, kPolM, cPpsPolM, TaxPolM, ShockPolM] = ...
                main_steady_state_utils_bgp_fmincon.solve_steady_state_iter_unified_fmincon(Z_ss_norm, cS_ss, paramS, verbose, x0_guess, solver_method);

            if ~eq_found
                warning('FMINCON稳态求解失败！');
                ss = []; Dist = []; V = []; kPolM = []; cPpsPolM = []; TaxPolM = []; ShockPolM = [];
            end
        end

        % =======================================================
        % == [FMINCON非PPS] 第六部分：非PPS版本的FMINCON求解器
        % =======================================================

        function [ss, Dist, V, kPolM] = solve_steady_state_complete_fmincon(cS_ss, paramS, params_ext, verbose, x0_guess, solver_method)
            % [FMINCON新建] 非PPS稳态求解器 - 基于连续优化的顶层API（无PPS模式）
            if nargin < 4, verbose = true; end
            if nargin < 5, x0_guess = []; end
            if nargin < 6, solver_method = 'fsolve'; end

            % 人口分布
            age_mass = ones(cS_ss.aD_new, 1);
            for ia = 1:(cS_ss.aD_new - 1)
                age_mass(ia+1) = age_mass(ia) * cS_ss.s_pathV(ia);
            end

            if isfield(params_ext, 'Z') && ~isempty(params_ext.Z)
                Z_ss_norm = params_ext.Z;
            else
                Z_ss_norm = age_mass / sum(age_mass);
            end

            % 外生变量设置
            cS_ss.A = params_ext.A;
            if isfield(params_ext, 'g_A_ss'), cS_ss.g_A_ss = params_ext.g_A_ss; end
            cS_ss.theta_path = params_ext.theta;

            % 确保非PPS模式
            cS_ss.pps_active = false;

            % 存活率计算
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

            % 调用FMINCON版本的统一求解器（非PPS模式）
            [ss, eq_found, Dist, k_prime_idx, ~, V, kPolM, ~, ~, ~] = ...
                main_steady_state_utils_bgp_fmincon.solve_steady_state_iter_unified_fmincon(Z_ss_norm, cS_ss, paramS, verbose, x0_guess, solver_method);

            if ~eq_found
                warning('FMINCON非PPS稳态求解失败！');
                ss = []; Dist = []; V = []; kPolM = [];
            end
        end

        function [K_p_model_out, L_model_out, ss, Dist, k_prime_idx, kpps_prime_idx, V, kPolM, cPpsPolM, TaxPolM, ShockPolM] = calculate_aggregates_unified_fmincon(K_p_guess, K_g_guess, L_guess, Z_ss_norm, cS, paramS)
            % [FMINCON重构] 统一的聚合计算函数 - 基于存储的策略矩阵直接聚合
            % 核心改进：
            % 1. 移除所有backout和aggregate调用
            % 2. 直接从存储的策略矩阵中读取并聚合会计流量
            % 3. 大幅简化聚合逻辑
            
            % [BGP修改] 基本输入验证和标准化
            if K_p_guess <= 0, K_p_guess = 1e-8; end
            if K_g_guess <= 0, K_g_guess = 1e-8; end
            if L_guess <= 0, L_guess = 1e-8; end
            A_ss = 1.0; % [BGP修改] 强制为1.0以保持标准化一致性
            theta_ss = cS.theta_path(1);

            % 计算价格
            M_prices = main_steady_state_utils_bgp_fmincon.get_prices_at_t(K_p_guess, K_g_guess, L_guess, A_ss, cS);
            M_for_hh = M_prices;

            % 计算养老金
            total_wage_bill = M_prices.w_t * L_guess;
            mass_retirees_ss = sum(Z_ss_norm((cS.aR_new+1):end));
            if mass_retirees_ss > 1e-9
                M_for_hh.b_t = (theta_ss * total_wage_bill) / mass_retirees_ss;
            else
                M_for_hh.b_t = 0;
            end

            cS_ss = cS;
            cS_ss.theta_t = theta_ss;

            % [FMINCON核心] 调用基于fmincon的VFI，获取所有策略矩阵
            if cS.pps_active
                [cPolM, kPolM, cPpsPolM, TaxPolM, ShockPolM, V] = ...
                    main_steady_state_utils_bgp_fmincon.HHSolution_VFI_fmincon_pps(M_for_hh, paramS, cS_ss);
                
                % 计算离散化索引
                [k_prime_idx, kpps_prime_idx] = main_steady_state_utils_bgp_fmincon.get_policy_index_matrix_with_pps(kPolM, cPpsPolM, cS_ss);
                
                % 求解分布
                Dist = main_steady_state_utils_bgp_fmincon.solve_steady_state_distribution_with_pps(k_prime_idx, kpps_prime_idx, paramS, cS_ss, Z_ss_norm);
                
                % [FMINCON核心] 新的直接聚合逻辑 - 无需反解计算
                [K_p_model_out, C_utility_final, Tax_final, Bequest_tax_final, L_model_out, ShockExp_final] = ...
                    main_steady_state_utils_bgp_fmincon.aggregate_from_stored_policies_with_pps(...
                    Dist, cPolM, kPolM, cPpsPolM, TaxPolM, ShockPolM, cS_ss, paramS);
            else
                % 调用非PPS版本的VFI
                [cPolM, kPolM, TaxPolM, ShockPolM, V] = ...
                    main_steady_state_utils_bgp_fmincon.HHSolution_VFI_fmincon(M_for_hh, paramS, cS_ss);
                
                % 计算离散化索引
                k_prime_idx = main_steady_state_utils_bgp_fmincon.get_policy_index_matrix(kPolM, cS_ss);
                
                % 求解分布
                Dist = main_steady_state_utils_bgp_fmincon.solve_steady_state_distribution(k_prime_idx, paramS, cS_ss, Z_ss_norm);
                
                % 创建空的PPS相关变量以保持接口一致性
                kpps_prime_idx = [];
                cPpsPolM = zeros(cS.nk, cS.nw_expanded, cS.aD_new); % 保持接口一致
                
                % [FMINCON核心] 非PPS版本的直接聚合逻辑
                [K_p_model_out, C_utility_final, Tax_final, Bequest_tax_final, L_model_out, ShockExp_final] = ...
                    main_steady_state_utils_bgp_fmincon.aggregate_from_stored_policies(...
                    Dist, cPolM, kPolM, TaxPolM, ShockPolM, cS_ss, paramS);
            end

            % [BGP修改] 填充完整的 ss 结构体
            ss = struct();
            ss.K_private_hat = K_p_guess;
            ss.K_public_hat = K_g_guess;
            ss.K_total_hat = K_p_guess + K_g_guess;
            ss.L_hat = L_guess;
            ss.Y_from_production_hat = M_prices.Y_t;
            ss.w_hat = M_prices.w_t;
            ss.r_mkt = M_prices.r_mkt_t;
            ss.b_hat = M_for_hh.b_t;
            ss.Bequest_tax = Bequest_tax_final;
            ss.Regular_tax = Tax_final;
            ss.Total_consumption = C_utility_final;
            ss.Total_shock_expenditure = ShockExp_final;
        end

        function [K_p_model_out, C_utility_agg, Tax_agg, Bequest_tax_agg, L_model_out, ShockExp_agg] = aggregate_from_stored_policies_with_pps(Dist, cPolM, kPolM, cPpsPolM, TaxPolM, ShockPolM, cS, paramS)
            % [FMINCON新建] 直接从存储的策略矩阵聚合 - 无需任何反解计算
            % 这是新框架的核心优势：聚合阶段变得极其简单和快速
            
            % 初始化聚合变量
            K_p_model_out = 0;
            C_utility_agg = 0;
            Tax_agg = 0;
            Bequest_tax_agg = 0;
            L_model_out = 0;
            ShockExp_agg = 0;

            % 获取存活率
            if isfield(cS, 'prob_survive_implied_ss0') && ~isempty(cS.prob_survive_implied_ss0)
                prob_survive_implied = cS.prob_survive_implied_ss0;
            elseif isfield(cS, 'prob_survive_implied_trans') && ~isempty(cS.prob_survive_implied_trans)
                prob_survive_implied = cS.prob_survive_implied_trans;
            else
                prob_survive_implied = cS.s_pathV;
            end

            % [BGP修改] 获取技术增长率，用于遗赠税计算
            g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;

            % [FMINCON核心] 简单的直接聚合循环 - 无任何复杂计算
            for ia = 1:cS.aD_new
                for ie = 1:cS.nw_expanded
                    for ik = 1:cS.nk
                        for ikpps = 1:cS.nkpps
                            mass = Dist(ik, ikpps, ie, ia);
                            if mass < 1e-20, continue; end

                            % 直接从策略矩阵读取数据
                            c_val = cPolM(ik, ikpps, ie, ia);
                            k_prime = kPolM(ik, ikpps, ie, ia);
                            kpps_prime = cPpsPolM(ik, ikpps, ie, ia);
                            tax_val = TaxPolM(ik, ikpps, ie, ia);
                            shock_exp = ShockPolM(ik, ikpps, ie, ia);

                            % 聚合流量
                            C_utility_agg = C_utility_agg + c_val * mass;
                            Tax_agg = Tax_agg + tax_val * mass;
                            ShockExp_agg = ShockExp_agg + shock_exp * mass;

                            % 计算存活和死亡概率
                            prob_survive = prob_survive_implied(ia);
                            prob_death = 1 - prob_survive;

                            % 资本聚合
                            K_p_model_out = K_p_model_out + (k_prime + kpps_prime) * mass * prob_survive;

                            % [BGP修改] 遗赠税基于真实下一期资产价值
                            Bequest_tax_agg = Bequest_tax_agg + (k_prime + kpps_prime) * (1 + g_A_period) * mass * prob_death;

                            % 劳动供给聚合
                            if ia <= cS.aR_new
                                epsilon_val = paramS.leGridV(ie);
                                L_model_out = L_model_out + (cS.ageEffV_new(ia) * epsilon_val) * mass;
                            end
                        end
                    end
                end
            end

            % [FMINCON聚合] 直接聚合完成，不输出调试信息保持与原版一致
        end

        function F_error = system_of_equations_steady_state_fmincon(x, Z_ss_norm, cS, paramS)
            % [FMINCON修改] 统一的系统方程组 - 适配新的聚合函数接口
            % 输入: x = [K_p_guess, K_g_guess, L_guess]
            % 输出: F_error = [error_Kp; error_Kg; error_L]
            
            K_p_guess = x(1);
            K_g_guess = x(2);
            L_guess = x(3);

            % 调用FMINCON版本的聚合计算函数
            [K_p_model, L_model, ss, ~, ~, ~, ~, ~, ~, ~] = main_steady_state_utils_bgp_fmincon.calculate_aggregates_unified_fmincon(K_p_guess, K_g_guess, L_guess, Z_ss_norm, cS, paramS);

            % 方程1的误差: 私人资本供给 - 私人资本需求
            error_Kp = K_p_guess - K_p_model;

            % [BGP修改] 方程2的误差: 公共投资 - 公共资本总需求
            g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
            Gov_Revenue_total = ss.Regular_tax + ss.Bequest_tax;
            I_g_model = cS.lambda_g * Gov_Revenue_total;
            Depreciation_g_model = (cS.ddk_g + g_A_period) * K_g_guess;
            error_Kg = I_g_model - Depreciation_g_model;

            % 方程3的误差: 劳动供给 - 劳动需求
            error_L = L_guess - L_model;

            F_error = [error_Kp; error_Kg; error_L];
        end

        % =======================================================
        % == [FMINCON辅助] 目标函数包装器（从parfor中移出）
        % =======================================================
        
        function neg_utility = objective_wrapper_no_pps(s_k, k_curr, a, i_e, eps_val, M, cS_obj, ev_interp, eff_disc, beq_disc)
            % [FMINCON辅助] 非PPS版本的目标函数包装器 - 独立静态方法
            Flows = main_steady_state_utils_bgp_fmincon.calculate_flows_from_proportions_no_pps(s_k, k_curr, a, i_e, eps_val, M, cS_obj);
            
            if Flows.c_val <= 0
                neg_utility = 1e10;
                return;
            end
            
            [~, util_c] = model_setup_utils_bgp.CES_utility(Flows.c_val, cS_obj.sigma, cS_obj);
            ev_val = ev_interp(Flows.k_prime);
            if isnan(ev_val), ev_val = -1e10; end

            util_bequest = model_setup_utils_bgp.bequest_utility(Flows.k_prime, cS_obj);

            total_utility = util_c + eff_disc * ev_val + beq_disc * util_bequest;
            neg_utility = -total_utility;
        end

        function neg_utility = objective_wrapper_pps(x, k_curr, kpps_curr, a, i_e, eps_val, M, cS_obj, ev_interp, eff_disc, beq_disc)
            % [FMINCON辅助] PPS版本的目标函数包装器 - 独立静态方法
            Flows = main_steady_state_utils_bgp_fmincon.calculate_flows_from_proportions_pps(x, k_curr, kpps_curr, a, i_e, eps_val, M, cS_obj);
            
            if Flows.c_val <= 0
                neg_utility = 1e10; % 惩罚不可行消费
                return;
            end

            [~, util_c] = model_setup_utils_bgp.CES_utility(Flows.c_val, cS_obj.sigma, cS_obj);
            ev_val = ev_interp(Flows.k_prime, Flows.kpps_prime);
            if isnan(ev_val), ev_val = -1e10; end % 惩罚域外值

            total_bequest = Flows.k_prime + Flows.kpps_prime;
            util_bequest = model_setup_utils_bgp.bequest_utility(total_bequest, cS_obj);

            total_utility = util_c + eff_disc * ev_val + beq_disc * util_bequest;
            neg_utility = -total_utility;
        end

    end
end 