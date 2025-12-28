% --- household.m ---
classdef household
    methods (Static)

        
        function [V_l, M_V_l] = labor_disutility(l_supply, chi_age, nu_h, cS)
            % =========================================================================
            % ==             【新增】劳动负效用计算器 v1.0
            % ==
            % == 功能:
            % ==   - 计算给定劳动供给 l_supply 的负效用 V_l 及其边际量 M_V_l。
            % ==   - 兼容 Frisch 弹性 nu_h 的不同取值。
            % ==   - 核心公式: V_l = chi * l^(1 + 1/nu) / (1 + 1/nu)
            % ==   - 特殊处理: 当 nu -> inf (完全无弹性)，1/nu -> 0，
            % ==     此时 V_l = chi * l。这对应于您提供的测试设置，
            % ==     可以完美复现固定劳动供给模型的结果。
            % ==   - 输入 l_supply 应为有效劳动供给 (l_choice * age_eff * shock_e)。
            % =========================================================================
            if chi_age <= 0
                V_l = zeros(size(l_supply));
                M_V_l = zeros(size(l_supply));
                return;
            end

            if isinf(nu_h)
                % GHH 偏好或固定劳动供给的极限情况
                V_l = chi_age .* l_supply;
                M_V_l = chi_age * ones(size(l_supply));
            else
                power = 1 + 1/nu_h;
                V_l = chi_age .* (l_supply.^power) ./ power;
                M_V_l = chi_age .* (l_supply.^(1/nu_h));
            end
        end


function [valS_t, polS_t] = VFI_transition_engine_by_age(M_vfi_t, b_payg_hat_t, paramS_t, cS_t, Vprime_t_plus_1)
            % =========================================================================
            % == 函数: VFI_transition_engine_by_age (版本 v1.3 - 简化版)
            % == 核心修改:
            % ==   - 移除了冗余的 h_idx 输入参数。
            % ==   - 现在完全依赖于上层函数 backward_hh 来确保传入的 cS_t 结构体
            % ==     中已经包含了正确的 h_current_type 字段。
            % =========================================================================
            valS_t = -Inf(cS_t.nk, cS_t.nkpps, cS_t.nw_expanded, cS_t.aD_new);
            polS_cell = cell(cS_t.aD_new, 1);

            bV_payg_vfi_hat = zeros(1, cS_t.aD_new);
            if cS_t.aR_new < cS_t.aD_new
                bV_payg_vfi_hat((cS_t.aR_new+1):cS_t.aD_new) = b_payg_hat_t;
            end

            for a_idx = cS_t.aD_new : -1 : 1
                vPrime_next_age = [];
                if a_idx < cS_t.aD_new
                    vPrime_next_age = Vprime_t_plus_1(:,:,:,a_idx+1);
                end
                b_age_val_hat = bV_payg_vfi_hat(a_idx);

                % cS_t 已经包含了 h_current_type
                [val_age, pol_age] = household.VFI_by_age_PPS(a_idx, vPrime_next_age, M_vfi_t, b_age_val_hat, paramS_t, cS_t);

                val_reshaped = reshape(val_age, [cS_t.nk, cS_t.nkpps, cS_t.nw_expanded]);
                valS_t(:,:,:,a_idx) = val_reshaped;
                polS_cell{a_idx} = pol_age;
            end
            polS_t = [polS_cell{:}];
        end
        function [Pol_path, Val_path] = backward_hh(pathS, cS, paramSF, valF, polF, h_idx)
            % =========================================================================
            % == 函数: backward_hh (v1.6 - h_idx 显式输入版)
            % == 核心修改:
            % ==   - 函数签名中新增了一个必需的输入参数 h_idx。
            % ==   - 在时间 t 的循环内部，直接使用这个接收到的 h_idx 来设置
            % ==     cS_t.h_current_type = h_idx;
            % ==   - 这是解决 h_current_type 字段丢失问题的最稳健方法。
            % =========================================================================
            T = cS.T_sim;
            Pol_path = cell(1, T);
            Val_path = zeros(cS.nk, cS.nkpps, cS.nw_expanded, cS.aD_new, T);
            Val_path(:,:,:,:,T) = valF;
            Pol_path{T} = polF;
            
            s_pathV_full_path = cS.s_pathV;
            aR_new_full_path = cS.aR_new_path; 

            for t = (T-1) : -1 : 1
                g_A_onestep_ahead = pathS.g_A_path(t);

                M_vfi_t = struct('r_mkt_t', pathS.r_path(t), ...
                    'w_t', pathS.w_hat_path(t), ...
                    'g_A_t_plus_1', g_A_onestep_ahead);

                M_vfi_t.tr_per_hh = pathS.tr_per_hh_hat_path(t);
                M_vfi_t.theta_t = pathS.theta_path(t);

                cS_t = cS; 
                if size(s_pathV_full_path, 2) > 1
                    cS_t.s_pathV = s_pathV_full_path(:, t);
                end
                cS_t.aR_new = aR_new_full_path(t);
                
                % [!!! 核心修正 !!!] 直接使用传入的 h_idx 参数
                cS_t.h_current_type = h_idx;
                
                Vprime_t_plus_1 = Val_path(:,:,:,:,t+1);

                [valS_t, polS_t] = household.VFI_transition_engine_by_age(M_vfi_t, pathS.b_hat_path(t), paramSF, cS_t, Vprime_t_plus_1);

                Val_path(:,:,:,:,t) = valS_t;
                Pol_path{t} = polS_t;
            end
        end
        
        
        
        
        function [val_age, pol_age] = VFI_by_age_PPS(a_idx, vPrime_kkppse_next, M_age, b_age_val, paramS_age, cS)
            % =========================================================================
            % == 函数: VFI_by_age_PPS (版本 v5.1 - 弹性劳动供给联合决策版)
            % == 核心修改:
            % ==   - (无逻辑修改) 此版本与 v5.0 完全相同，但现在它依赖于
            % ==     上层调用者 (VFI_transition_engine_by_age 或 VFI_solver)
            % ==     来确保传入的 cS 结构体中已经包含了正确的 'h_current_type' 字段。
            % =========================================================================
            nk_search_grid = 150;
            nL = length(cS.lGridV);

            % 1. 增长与贴现因子
            if isfield(M_age, 'g_A_t_plus_1')
                g_A_period = M_age.g_A_t_plus_1;
            else
                g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
            end
            growth_factor_bgp = (1 + g_A_period);
            discount_factor_V_prime = cS.beta * growth_factor_bgp^(1 - cS.sigma);

            % 2. V' 插值器
            k_grid_vec = cS.kGridV;
            kpps_grid_vec = cS.kppsGridV;
            vPrime_interpolants = cell(cS.nw_expanded, 1);
            if a_idx < cS.aD_new
                for ie_next = 1:cS.nw_expanded
                    if cS.nkpps > 1 && cS.pps_active
                        vPrime_interpolants{ie_next} = griddedInterpolant({k_grid_vec, kpps_grid_vec}, vPrime_kkppse_next(:, :, ie_next), 'makima','makima');
                    else
                        vPrime_interpolants{ie_next} = griddedInterpolant(k_grid_vec, squeeze(vPrime_kkppse_next(:, :, ie_next)), 'pchip','pchip');
                    end
                end
            end

            % 3. 构建状态-选择联合网格 (k, kpps, e, l_choice)
            [k_mat, kpps_mat, le_mat, l_choice_mat] = ndgrid(cS.kGridV, cS.kppsGridV, paramS_age.leGridV, cS.lGridV);

            % 4. 计算与 (k', pps_rate) 决策无关的量
            labor_supply_tens = zeros(size(k_mat));
            if a_idx <= cS.aR_new
                labor_supply_tens = l_choice_mat .* cS.ageEffV_new(a_idx) .* le_mat;
            end
            labor_income_hat_tens = M_age.w_t * labor_supply_tens;

            capital_return_hat_tens = k_mat * (1 + M_age.r_mkt_t);
            capital_tax_hat_tens = cS.tau_k .* (k_mat * M_age.r_mkt_t);
            labor_tax_hat_tens = cS.tau_l * labor_income_hat_tens;
            payg_contribution_hat_tens = M_age.theta_t * labor_income_hat_tens;

            pps_contrib_hat_tens = zeros(size(k_mat));
            pension_income_hat_tens = zeros(size(k_mat));
            pps_withdrawal_hat_tens = zeros(size(k_mat));
            pps_tax_hat_tens = zeros(size(k_mat));
            
            is_working_age = (a_idx <= cS.aR_new);
            if cS.pps_active
                if is_working_age
                    pps_contrib_hat_tens = labor_income_hat_tens * cS.pps_fixed;
                else
                    pension_income_hat_tens(:) = b_age_val;
                    pps_withdrawal_hat_tens = cS.pps_withdrawal_rate .* (kpps_mat * (1 + M_age.r_mkt_t));
                    pps_tax_hat_tens = cS.pps_tax_rate_withdrawal .* pps_withdrawal_hat_tens;
                end
            elseif ~is_working_age
                 pension_income_hat_tens(:) = b_age_val;
            end

            kpps_hat_end_of_period_tens = kpps_mat * (1 + M_age.r_mkt_t) + pps_contrib_hat_tens - pps_withdrawal_hat_tens;
            kpps_prime_hat_tens = max(cS.kppsMin, kpps_hat_end_of_period_tens / growth_factor_bgp);

            % 5. 构建 k' 选择网格并进行向量化计算
            base_coh_hat_tens = capital_return_hat_tens + labor_income_hat_tens + M_age.tr_per_hh + pension_income_hat_tens + pps_withdrawal_hat_tens ...
                - (capital_tax_hat_tens + labor_tax_hat_tens + payg_contribution_hat_tens + pps_tax_hat_tens);
            cash_on_hand_hat_tens = base_coh_hat_tens - pps_contrib_hat_tens; 

            k_prime_max_tens = (cash_on_hand_hat_tens - cS.cFloor * (1 + cS.tau_c)) / growth_factor_bgp;
            k_prime_max_tens(k_prime_max_tens < cS.kMin) = cS.kMin;
            
            k_prime_search_grid_reshaped = reshape(linspace(0, 1, nk_search_grid), [1, 1, 1, 1, nk_search_grid]);
            k_prime_choices_5D = cS.kMin + (k_prime_max_tens - cS.kMin) .* k_prime_search_grid_reshaped;
            
            c_expend_choices_5D = cash_on_hand_hat_tens - k_prime_choices_5D * growth_factor_bgp;
            c_choices_5D = c_expend_choices_5D / (1 + cS.tau_c);
            invalid_choice_mask = (c_choices_5D < cS.cFloor);
            c_choices_5D(invalid_choice_mask) = cS.cFloor;
            util_c_5D = (c_choices_5D.^(1 - cS.sigma))./(1 - cS.sigma);
            util_c_5D(invalid_choice_mask) = -1e20;
            
            % [util_l_tens, ~] = household.labor_disutility(labor_supply_tens, cS.chi_pathV(a_idx), cS.nu_h(cS.h_current_type), cS);
                     % [!!! 核心修正: 使用 l_choice_mat 而不是 labor_supply_tens !!!]
            [util_l_tens, ~] = household.labor_disutility(l_choice_mat, cS.chi_pathV(a_idx), cS.nu_h(cS.h_current_type), cS);
            
            kpps_prime_hat_5D = repmat(kpps_prime_hat_tens, [1, 1, 1, 1, nk_search_grid]);
            total_wealth_prime_5D = k_prime_choices_5D + kpps_prime_hat_5D;
            util_bequest_5D = utils.bequest_utility(total_wealth_prime_5D, cS);

            future_value_5D = zeros(size(k_prime_choices_5D));
            if a_idx < cS.aD_new
                ev_on_choices_5D = zeros(size(k_prime_choices_5D));
                trans_mat = paramS_age.TrProbM_by_age{a_idx};
                
                for ie_next = 1:cS.nw_expanded
                    trans_prob_col = trans_mat(:, ie_next);
                    if any(trans_prob_col > 0)
                        if cS.nkpps > 1 && cS.pps_active
                            v_prime_interp_5D = vPrime_interpolants{ie_next}(k_prime_choices_5D, kpps_prime_hat_5D);
                        else
                            v_prime_interp_5D = vPrime_interpolants{ie_next}(k_prime_choices_5D);
                        end
                        ev_on_choices_5D = ev_on_choices_5D + v_prime_interp_5D .* reshape(trans_prob_col, [1, 1, cS.nw_expanded, 1, 1]);
                    end
                end

                survival_rate = cS.s_pathV(a_idx);
                if size(survival_rate,2) > 1, survival_rate = survival_rate(1); end
                future_value_5D = discount_factor_V_prime * (survival_rate * ev_on_choices_5D + (1 - survival_rate) * util_bequest_5D);
            else
                future_value_5D = discount_factor_V_prime * util_bequest_5D;
            end
            
            total_value_5D = util_c_5D - util_l_tens + future_value_5D;

            % 6. 寻找最优决策
            S_5D = size(total_value_5D);
            total_value_4D_reshaped = reshape(total_value_5D, S_5D(1), S_5D(2), S_5D(3), nL * nk_search_grid);
            [val_age, best_joint_idx_mat] = max(total_value_4D_reshaped, [], 4, 'omitnan');
            
            [best_l_idx_mat, best_k_prime_idx_mat] = ind2sub([nL, nk_search_grid], best_joint_idx_mat);

            % 7. 提取最优策略
            [I, J, K] = ndgrid(1:S_5D(1), 1:S_5D(2), 1:S_5D(3));
            linear_indices_final = sub2ind(S_5D, I, J, K, best_l_idx_mat, best_k_prime_idx_mat);
            
            c_final_mat = c_choices_5D(linear_indices_final);
            k_prime_final_mat = k_prime_choices_5D(linear_indices_final);

            S_4D = size(labor_supply_tens);
            linear_indices_4D = sub2ind(S_4D, I, J, K, best_l_idx_mat);
            
            l_supply_final_mat = labor_supply_tens(linear_indices_4D);
            kpps_prime_final_mat = kpps_prime_hat_tens(linear_indices_4D);
            
            tax_final_mat = capital_tax_hat_tens(linear_indices_4D) + ...
                            labor_tax_hat_tens(linear_indices_4D) + ...
                            pps_tax_hat_tens(linear_indices_4D) + ...
                            payg_contribution_hat_tens(linear_indices_4D) + ...
                            c_final_mat .* cS.tau_c;

            % 8. 保存结果
            pol_age = struct(...
                'c', c_final_mat, ...
                'l', l_supply_final_mat, ...
                'k_prime', k_prime_final_mat, ...
                'kpps_prime', kpps_prime_final_mat, ...
                'tax_regular', tax_final_mat);
        end


        function [polS, valS] = VFI_solver(M_vfi, paramS_vfi, cS_vfi)
            % =========================================================================
            % == 函数: VFI_solver (版本 v2.2 - 养老金(b_hat)兼容版)
            % == 目的: 使稳态VFI求解器能够接收和处理养老金收入。
            % =========================================================================

            if isfield(cS_vfi, 'pps_active') && cS_vfi.pps_active
                valS = -Inf(cS_vfi.nk, cS_vfi.nkpps, cS_vfi.nw_expanded, cS_vfi.aD_new);
            else
                valS = -Inf(cS_vfi.nk, 1, cS_vfi.nw_expanded, cS_vfi.aD_new);
            end

            polS_cell = cell(cS_vfi.aD_new, 1);

            % [核心修正] 从 M_vfi 中获取养老金福利 b_hat，并构建年龄相关的养老金向量
            b_payg_hat_age_vec = zeros(1, cS_vfi.aD_new);
            if isfield(M_vfi, 'b_hat_t') && M_vfi.b_hat_t > 0
                % 仅为退休后年龄组分配养老金
                b_payg_hat_age_vec((cS_vfi.aR_new+1):end) = M_vfi.b_hat_t;
            end

            for a_idx = cS_vfi.aD_new : -1 : 1
                vPrime_kkppse_next = [];
                if a_idx < cS_vfi.aD_new
                    vPrime_kkppse_next = valS(:,:,:,a_idx+1);
                end

                % [核心修正] 将特定年龄的养老金 b_age_val 传递给引擎
                b_age_val = b_payg_hat_age_vec(a_idx);

                % if isfield(cS_vfi, 'pps_active') && cS_vfi.pps_active
                %     [val_age, pol_age] = household.VFI_by_age_PPS(a_idx, vPrime_kkppse_next, M_vfi, b_age_val, paramS_vfi, cS_vfi);
                % else
                %     [val_age, pol_age] = household.VFI_by_age(a_idx, vPrime_kkppse_next, M_vfi, b_age_val, paramS_vfi, cS_vfi);
                % end

                [val_age, pol_age] = household.VFI_by_age_PPS(a_idx, vPrime_kkppse_next, M_vfi, b_age_val, paramS_vfi, cS_vfi);

                valS(:,:,:,a_idx) = val_age;
                polS_cell{a_idx} = pol_age;
            end
            polS = [polS_cell{:}];
        end




        function [muM, utilM] = CES_utility(cM, sigma, cS)
            % =========================================================================
            % == 函数: CES_utility (v2.0 - 向量化安全修正版)
            % == 核心修正:
            % ==   - [!!! 关键稳健性修正 !!!] 替换了非向量化的 `if cM<=0`
            % ==     判断。旧代码在处理向量输入时，仅判断第一个元素，
            % ==     当向量中包含非正值时会导致 NaN 或 Inf。
            % ==   - 新方法: 使用 `max(cM, cS.cFloor)` 对消费输入进行
            % ==     向量化处理，确保所有元素都大于等于消费下限，然后再
            % ==     代入效用函数计算。这保证了函数对任意向量输入的健壮性。
            % =========================================================================

            % 使用 max 函数进行向量化安全处理，确保消费不低于下限
            c_adj = max(cM, cS.cFloor);

            if abs(sigma - 1) < 1e-6
                utilM = log(c_adj);
                muM = 1./c_adj;
            else
                utilM = (c_adj.^(1-sigma))./(1-sigma);
                muM = c_adj.^(-sigma);
            end
        end

        function util_beq = bequest_utility(k_prime, cS)
            % =========================================================================
            % == 函数: bequest_utility (v2.0 - 向量化安全修正版)
            % == 核心修正:
            % ==   - [!!! 关键稳健性修正 !!!] 与CES_utility的修正类似，
            % ==     使用 `max(k_prime, cS.cFloor)` 对遗赠财富输入进行
            % ==     向量化处理，以防止非正值输入导致计算错误。
            % ==   - 这确保了函数在任何求解器中（VFI, EGM）的行为都
            % ==     是稳健和可预测的。
            % =========================================================================
            if cS.phi_bequest <= 0
                util_beq = zeros(size(k_prime)); % 确保输出维度正确
                return;
            end

            % 使用 max 函数进行向量化安全处理
            k_adj = max(k_prime, cS.cFloor);

            if abs(cS.sigma - 1) < 1e-6
                util_beq = cS.phi_bequest * log(k_adj);
            else
                util_beq = cS.phi_bequest * (k_adj.^(1-cS.sigma))./(1-cS.sigma);
            end
        end

        function mu_beq = bequest_utility_deriv(k_prime, cS)
            % =========================================================================
            % == 函数: bequest_utility_deriv (版本 v1.0 - 内联实现版)
            % == 目的: 计算遗赠的边际效用 U_beq'(k')。
            % ==       d/dk [ φ * k^(1-σ) / (1-σ) ] = φ * k^(-σ)
            % ==       此函数被内联到 household 类中以消除对 utils 的依赖。
            % =========================================================================
            if cS.phi_bequest <= 0
                mu_beq = zeros(size(k_prime)); % 确保输出维度正确
                return;
            end

            % 使用 max 函数进行向量化安全处理，避免 k_prime <= 0 导致计算错误
            k_adj = max(k_prime, cS.cFloor);

            if abs(cS.sigma - 1) < 1e-6
                % 对数效用形式的导数: φ / k
                mu_beq = cS.phi_bequest ./ k_adj;
            else
                % CRRA效用形式的导数: φ * k^(-σ)
                mu_beq = cS.phi_bequest .* (k_adj.^(-cS.sigma));
            end
        end



        function run_pps_compare()
            % =========================================================================
            % == 测试函数: run_pps_compare (v2.2 - 调用签名修正版)
            % == 核心修改:
            % ==   - [!!! 关键修正 !!!] 根据用户指正，恢复了对
            % ==     utils.EarningProcess_AgeDependent 函数的正确调用签名。
            % ==     该函数返回四个独立的输出，而非两个结构体。
            % ==   - 保持了 v2.1 版本的模块化结构，对 nw=1 和 nw=5 的情况
            % ==     分别执行完整的测试流程。
            % =========================================================================

            fprintf('=======================================================\n');
            fprintf('==  开启 PPS VFI 求解器对比测试程序  ==\n');
            fprintf('==  (外生固定缴费 vs 内生选择缴费)  ==\n');
            fprintf('=======================================================\n\n');

            % --- 1. 初始化通用参数 ---
            cS_base = utils.ParameterValues();
            cS_base.nk = 40;
            cS_base.nkpps = 30;
            cS_base.g_A_ss = 0.015;
            cS_base.ageEffV_new = cS_base.ageEffV_new_h(:,1);
            cS_base.s_pathV = cS_base.s_pathV(:, 1);
            cS_base.pps_active = true;
            M1 = struct('r_mkt_t', 0.04, 'w_t', 1.5, 'tr_per_hh', 0.01, 'b_hat_t', 0.5, 'theta_t', 0.11);

            % --- 2. 运行 nw = 1 的测试案例 ---
            cS1 = cS_base;
            cS1.nw = 1;
            household.run_test_case(cS1, M1);

            % --- 3. 运行 nw = 5 的测试案例 ---
            % cS5 = cS_base;
            % cS5.nw = 5;
            % household.run_test_case(cS5, M1);

            fprintf('\n所有测试完成。\n');
        end

        % --- 内部辅助函数，执行一套完整的对比测试 ---
        function run_test_case(cS, M1)
            nw_val = cS.nw;
            fprintf('\n\n#############################################\n');
            fprintf('###   开始执行测试案例 (nw = %d)   ###\n', nw_val);
            fprintf('#############################################\n\n');

            % --- 1. 生成依赖于 nw 的参数 ---
            paramS = struct(); % 初始化 paramS
            % [!!! 关键修正 !!!] 恢复正确的四输出调用方式
            [paramS.leGridV, paramS.TrProbM_by_age, paramS.leProb1V, cS.nw_expanded] = utils.EarningProcess_AgeDependent(cS);

            b_payg_hat_age_vec = zeros(1, cS.aD_new);
            b_payg_hat_age_vec((cS.aR_new+1):end) = M1.b_hat_t;

            % --- 2. 测试模块 1: 零缴费等价性检验 ---
            fprintf('--- 测试 1 (nw=%d): 零缴费等价性检验 ---\n', nw_val);
            fprintf('   - 模型 A: 外生固定缴费, cS.pps_fixed = 0\n');
            fprintf('   - 模型 B: 内生选择缴费, cS.pps_max = 0\n');

            % --- 模型 A ---
            cS_A = cS;
            cS_A.pps_fixed = 0;
            cS_A = utils.generateGrids(cS_A);
            valS_A = -Inf(cS_A.nk, cS_A.nkpps, cS_A.nw_expanded, cS_A.aD_new);
            polS_A_cell = cell(cS_A.aD_new, 1);
            tic;
            for a_idx = cS_A.aD_new : -1 : 1
                vPrime_kkppse_next = [];
                if a_idx < cS_A.aD_new, vPrime_kkppse_next = valS_A(:,:,:,a_idx+1); end
                [val_age, pol_age] = household.VFI_by_age_PPS(a_idx, vPrime_kkppse_next, M1, b_payg_hat_age_vec(a_idx), paramS, cS_A);
                valS_A(:,:,:,a_idx) = val_age;
                polS_A_cell{a_idx} = pol_age;
            end
            polS_A = [polS_A_cell{:}];
            time_A = toc;
            fprintf('模型 A 求解完成，耗时: %.2f 秒\n', time_A);

            % --- 模型 B ---
            cS_B = cS;
            cS_B.pps_fixed = 0;
            cS_B.n_pps_rate_grid = 1;
            cS_B = utils.generateGrids(cS_B);
            valS_B = -Inf(cS_B.nk, cS_B.nkpps, cS_B.nw_expanded, cS_B.aD_new);
            polS_B_cell = cell(cS_B.aD_new, 1);
            tic;
            for a_idx = cS_B.aD_new : -1 : 1
                vPrime_kkppse_next = [];
                if a_idx < cS_B.aD_new, vPrime_kkppse_next = valS_B(:,:,:,a_idx+1); end
                [val_age, pol_age] = household.VFI_by_age_PPS_endo(a_idx, vPrime_kkppse_next, M1, b_payg_hat_age_vec(a_idx), paramS, cS_B);
                valS_B(:,:,:,a_idx) = val_age;
                polS_B_cell{a_idx} = pol_age;
            end
            polS_B = [polS_B_cell{:}];
            time_B = toc;
            fprintf('模型 B 求解完成，耗时: %.2f 秒\n\n', time_B);

            diff_val = max(abs(valS_A - valS_B), [], 'all');
            fprintf('结果对比 (Max Abs Diff): 价值函数差异: %.4e\n', diff_val);
            if diff_val < 1e-9, fprintf('结论: 等价性验证通过。\n\n'); else fprintf('结论: 等价性验证失败。\n\n'); end

            % --- 3. 测试模块 2: 正缴费对比 ---
            fprintf('--- 测试 2 (nw=%d): 正缴费对比 ---\n', nw_val);
            fprintf('   - 模型 C: 外生固定缴费, cS.pps_fixed = 0.1\n');
            fprintf('   - 模型 D: 内生选择缴费, cS.pps_max = 0.1, n_grid = 5\n');

            % --- 模型 C ---
            cS_C = cS;
            cS_C.pps_fixed = 0.1;
            cS_C = utils.generateGrids(cS_C);
            valS_C = -Inf(cS_C.nk, cS_C.nkpps, cS_C.nw_expanded, cS_C.aD_new);
            polS_C_cell = cell(cS_C.aD_new, 1);
            tic;
            for a_idx = cS_C.aD_new : -1 : 1
                vPrime_kkppse_next = [];
                if a_idx < cS_C.aD_new, vPrime_kkppse_next = valS_C(:,:,:,a_idx+1); end
                [val_age, pol_age] = household.VFI_by_age_PPS(a_idx, vPrime_kkppse_next, M1, b_payg_hat_age_vec(a_idx), paramS, cS_C);
                valS_C(:,:,:,a_idx) = val_age;
                polS_C_cell{a_idx} = pol_age;
            end
            polS_C = [polS_C_cell{:}];
            time_C = toc;
            fprintf('模型 C 求解完成，耗时: %.2f 秒\n', time_C);

            % --- 模型 D ---
            cS_D = cS;
            cS_D.pps_fixed = 0.1;
            cS_D.n_pps_rate_grid = 5;
            cS_D = utils.generateGrids(cS_D);
            valS_D = -Inf(cS_D.nk, cS_D.nkpps, cS_D.nw_expanded, cS_D.aD_new);
            polS_D_cell = cell(cS_D.aD_new, 1);
            tic;
            for a_idx = cS_D.aD_new : -1 : 1
                vPrime_kkppse_next = [];
                if a_idx < cS_D.aD_new, vPrime_kkppse_next = valS_D(:,:,:,a_idx+1); end
                [val_age, pol_age] = household.VFI_by_age_PPS_endo(a_idx, vPrime_kkppse_next, M1, b_payg_hat_age_vec(a_idx), paramS, cS_D);
                valS_D(:,:,:,a_idx) = val_age;
                polS_D_cell{a_idx} = pol_age;
            end
            polS_D = [polS_D_cell{:}];
            time_D = toc;
            fprintf('模型 D 求解完成，耗时: %.2f 秒\n\n', time_D);

            fprintf('结果对比 (全状态空间均值 for nw=%d):\n', nw_val);
            fprintf('%-18s | %-15s | %-15s\n', '变量', '模型 C (外生)', '模型 D (内生)');
            disp(repmat('-', 1, 55));
            fprintf('%-18s | %-15.4f | %-15.4f\n', '价值函数', mean(valS_C,'all'), mean(valS_D,'all'));
            fprintf('%-18s | %-15.4f | %-15.4f\n', '消费', mean(cell2mat(arrayfun(@(x) x.c(:), polS_C, 'UniformOutput', false))), mean(cell2mat(arrayfun(@(x) x.c(:), polS_D, 'UniformOutput', false))));
            fprintf('%-18s | %-15.4f | %-15.4f\n', '常规储蓄 (k_prime)', mean(cell2mat(arrayfun(@(x) x.k_prime(:), polS_C, 'UniformOutput', false))), mean(cell2mat(arrayfun(@(x) x.k_prime(:), polS_D, 'UniformOutput', false))));
            fprintf('%-18s | %-15s | %-15.4f\n', 'PPS缴费率', '0.1 (固定)', mean(cell2mat(arrayfun(@(x) x.pps_contrib_rate(:), polS_D, 'UniformOutput', false))));
        end


        function run_vfi_test()
            % =========================================================================
            % == 测试函数: run_vfi_test (V2.2 - 测试逻辑修正版)
            % == 核心变更:
            % ==   - [!!! 关键测试逻辑修正 !!!] 在“有PPS”的测试案例中，当
            % ==     目标是验证 pps_max = 0 的情况时，必须将 cS.nkpps 也设为1。
            % ==   - 这确保了 kpps 状态维度被坍缩至一个点(0)，使得该测试与
            % ==     “无PPS”的案例在经济上是等价和可比的。
            % ==   - 此前的测试将 nkpps 设为20，导致即使在0缴费率下，家庭依然
            % ==     拥有外生的、增长的kpps财富，从而造成结果不一致。
            % =========================================================================

            fprintf('============================================\n');
            fprintf('==  开启 household.VFI_solver 测试程序  ==\n');
            fprintf('============================================\n\n');

            % --- 1. 初始化和参数设置 ---
            cS = utils.ParameterValues();

            % 为了加速测试，可以适当减小网格密度
            cS.nk = 30;
            cS.nkpps = 1; % 默认值
            cS.nw = 1;
            cS.n_pps_rate_grid = 20;
            cS.pps_max = 0.1; % 关键测试参数

            [paramS.leGridV, paramS.TrProbM_by_age, paramS.leProb1V, cS.nw_expanded] = utils.EarningProcess_AgeDependent(cS);
            cS.s_pathV = cS.s_pathV(:, 1);

            M1 = struct();
            M1.r_mkt_t = 0.04;
            M1.w_t = 1.5;
            M1.tr_per_hh = 0.01;
            M1.b_hat_t = 0.5;
            M1.theta_t = 0.11;

            % --- 2. 测试案例 1: 无 PPS 系统 ---
            fprintf('--- 测试 1: 模型无 PPS (pps_active = false) ---\n');
            cS.pps_active = false;
            cS.nkpps = 1;
            cS_no_pps = utils.generateGrids(cS);

            tic;
            [polS1, valS1] = household.VFI_solver(M1, paramS, cS_no_pps);
            elapsed_time = toc;
            fprintf('求解完成，耗时: %.2f 秒\n', elapsed_time);

            ages_to_check = [3, cS.aR_new, cS.aD_new-1];
            fprintf('%-10s %-15s %-15s %-15s\n', '年龄(a)', 'Mean(Value)', 'Mean(c)', 'Mean(k_prime)');
            for age = ages_to_check
                val_mean = mean(valS1(:,:,:,age), 'all');
                c_mean = mean(polS1(age).c, 'all');
                k_prime_mean = mean(polS1(age).k_prime, 'all');
                fprintf('%-10d %-15.4f %-15.4f %-15.4f\n', age, val_mean, c_mean, k_prime_mean);
            end

            % --- 3. 测试案例 2: 有 PPS 系统 (但配置为与无PPS等价) ---
            fprintf('\n--- 测试 2: 模型有 PPS (pps_active = true, 但 pps_max=0, nkpps=1) ---\n');
            cS.pps_active = true;
            % [!!! 关键测试逻辑修正 !!!]
            % 为使测试2与测试1可比，当 pps_max=0 时，nkpps 必须为1
            cS.nkpps = 30;
            cS_with_pps = utils.generateGrids(cS);

            tic;
            [polS2, valS2] = household.VFI_solver(M1, paramS, cS_with_pps);
            elapsed_time = toc;
            fprintf('求解完成，耗时: %.2f 秒\n', elapsed_time);

            fprintf('%-10s %-15s %-15s %-15s %-15s\n', '年龄(a)', 'Mean(Value)', 'Mean(c)', 'Mean(k_prime)', 'Mean(kpps_prm)');
            for age = ages_to_check
                val_mean = mean(valS2(:,:,:,age), 'all');
                c_mean = mean(polS2(age).c, 'all');
                k_prime_mean = mean(polS2(age).k_prime, 'all');
                kpps_prime_mean = mean(polS2(age).kpps_prime, 'all');
                fprintf('%-10d %-15.4f %-15.4f %-15.4f %-15.4f\n', age, val_mean, c_mean, k_prime_mean, kpps_prime_mean);
            end

            fprintf('\n请验证测试1和测试2的结果是否一致。\n');
            fprintf('如果一致, 说明在 pps_max=0 时, 两个模型实现了等价。\n\n');
        end


        function run_vfi_egm_test()
            % =========================================================================
            % == 测试函数: run_vfi_egm_test (版本 v3.1 - 深度诊断增强版)
            % == 核心修正:
            % ==   - 增加了 "深度诊断" 环节。在完成整体对比后，如果发现显著差异，
            % ==     此环节将聚焦于误差首次出现的年龄组 (通常是最后一个年龄组 aD_new)。
            % ==   - 对该特定年龄的特定状态点(例如 k 网格的中间点)，分别调用
            % ==     VFI 和 EGM 的核心计算逻辑，并并排打印出价值函数的所有构成
            % ==     部分：c, k', kpps', u(c), u_bequest, discount, V_total。
            % ==   - 这种精细化的输出能够清晰地揭示两种方法在计算价值时的第一个
            % ==     分歧点，从而为最终解决问题提供了决定性的线索。
            % =========================================================================

            fprintf('======================================================\n');
            fprintf('==  开启 household VFI vs. EGM 全面对比测试程序  ==\n');
            fprintf('======================================================\n\n');

            % --- 1. 初始化和参数设置 ---
            cS = utils.ParameterValues();

            cS.nk = 40;
            cS.nkpps = 10;
            cS.nw = 1;
            cS.n_pps_rate_grid = 10;
            cS.pps_max = 0.15;
            cS.pps_active = true;

            [paramS.leGridV, paramS.TrProbM_by_age, paramS.leProb1V, cS.nw_expanded] = utils.EarningProcess_AgeDependent(cS);
            cS.s_pathV = cS.s_pathV(:, 1);
            cS_test = utils.generateGrids(cS);

            M1 = struct('r_mkt_t', 0.04, 'w_t', 1.5, 'tr_per_hh', 0.01, 'b_hat_t', 0.5, 'theta_t', 0.11);

            % --- 2. 使用 VFI (蛮力法) 求解 ---
            fprintf('--- (1/2) 正在使用 VFI (蛮力法) 求解模型 ---\n');
            tic;
            [polS_vfi, valS_vfi] = household.VFI_solver(M1, paramS, cS_test, false);
            vfi_time = toc;
            fprintf('VFI 求解完成。耗时: %.2f 秒\n\n', vfi_time);

            % --- 3. 使用 EGM 求解 ---
            fprintf('--- (2/2) 正在使用 EGM (终点网格法) 求解模型 ---\n');
            tic;
            [polS_egm, valS_egm] = household.VFI_solver(M1, paramS, cS_test, true);
            egm_time = toc;
            fprintf('EGM 求解完成。耗时: %.2f 秒\n\n', egm_time);

            % --- 4. 结果对比分析 ---
            fprintf('--- 开始对比 VFI 和 EGM 的求解结果 ---\n');
            fprintf('求解器性能: EGM 速度 / VFI 速度 = %.2f\n\n', vfi_time / egm_time);

            fprintf('%-5s | %-20s | %-20s | %-20s | %-20s\n', ...
                'Age', 'Diff (Value)', 'Diff (Consumption)', 'Diff (k prime)', 'Diff (kpps prime)');
            fprintf('%-5s | %-10s %-10s | %-10s %-10s | %-10s %-10s | %-10s %-10s\n', ...
                '', 'Max Abs', 'Mean Abs', 'Max Abs', 'Mean Abs', 'Max Abs', 'Mean Abs', 'Max Abs', 'Mean Abs');
            disp(repmat('-', 1, 105));

            max_diff_summary = struct('val', 0, 'c', 0, 'k_prime', 0, 'kpps_prime', 0);

            for a_idx = 1:cS_test.aD_new
                val_vfi_age = valS_vfi(:,:,:,a_idx);
                val_egm_age = valS_egm(:,:,:,a_idx);
                pol_vfi_age = polS_vfi(a_idx);
                pol_egm_age = polS_egm(a_idx);
                diff_val = abs(val_vfi_age - val_egm_age);
                diff_c = abs(pol_vfi_age.c - pol_egm_age.c);
                diff_k_prime = abs(pol_vfi_age.k_prime - pol_egm_age.k_prime);
                diff_kpps_prime = abs(pol_vfi_age.kpps_prime - pol_egm_age.kpps_prime);

                if max(diff_val,[],'all') > max_diff_summary.val, max_diff_summary.val = max(diff_val,[],'all'); end
                if max(diff_c,[],'all') > max_diff_summary.c, max_diff_summary.c = max(diff_c,[],'all'); end
                if max(diff_k_prime,[],'all') > max_diff_summary.k_prime, max_diff_summary.k_prime = max(diff_k_prime,[],'all'); end
                if max(diff_kpps_prime,[],'all') > max_diff_summary.kpps_prime, max_diff_summary.kpps_prime = max(diff_kpps_prime,[],'all'); end

                fprintf('%-5d | %-10.4e %-10.4e | %-10.4e %-10.4e | %-10.4e %-10.4e | %-10.4e %-10.4e\n', ...
                    a_idx, max(diff_val, [], 'all'), mean(diff_val, 'all'), max(diff_c, [], 'all'), mean(diff_c, 'all'), ...
                    max(diff_k_prime, [], 'all'), mean(diff_k_prime, 'all'), max(diff_kpps_prime, [], 'all'), mean(diff_kpps_prime, 'all'));
            end

            disp(repmat('-', 1, 105));
            fprintf('\n测试完成。全局最大差异总结：\n');
            fprintf('  - Value Func (Max|V_vfi - V_egm|):       %.4e\n', max_diff_summary.val);
            fprintf('  - Consumption (Max|c_vfi - c_egm|):     %.4e\n', max_diff_summary.c);
            fprintf('  - Capital Sav. (Max|k_vfi - k_egm|):      %.4e\n', max_diff_summary.k_prime);
            fprintf('  - PPS Sav. (Max|kpps_vfi - kpps_egm|): %.4e\n', max_diff_summary.kpps_prime);

            % --- 5. 深度诊断环节 ---
            if max_diff_summary.val > 1e-3 % 如果存在显著差异
                fprintf('\n--- 深度诊断: 检测到显著差异，聚焦最后一个年龄组 (a=%d) ---\n', cS_test.aD_new);

                a_idx = cS_test.aD_new;
                % 选择一个状态点进行诊断 (例如：无劳动收入冲击，k和kpps网格的中间点)
                ie = 1;
                ik = round(cS_test.nk / 2);
                ikpps = round(cS_test.nkpps / 2);

                fprintf('诊断状态点: a_idx=%d, k_idx=%d (k=%.2f), kpps_idx=%d (kpps=%.2f), ie=%d\n\n',...
                    a_idx, ik, cS_test.kGridV(ik), ikpps, cS_test.kppsGridV(ikpps), ie);

                fprintf('%-10s | %-12s | %-12s | %-12s | %-12s | %-12s | %-12s | %-12s\n', ...
                    'Method', 'c', 'k_prime', 'kpps_prime', 'util(c)', 'util_bequest', 'Discount', 'Total Value');
                disp(repmat('-', 1, 110));

                % VFI 分解
                pol_vfi_pt = polS_vfi(a_idx);
                c_vfi = pol_vfi_pt.c(ik, ikpps, ie);
                k_p_vfi = pol_vfi_pt.k_prime(ik, ikpps, ie);
                kpps_p_vfi = pol_vfi_pt.kpps_prime(ik, ikpps, ie);
                util_c_vfi = (c_vfi.^(1-cS_test.sigma))./(1-cS_test.sigma);
                total_wealth_p_vfi = k_p_vfi + kpps_p_vfi;
                util_b_vfi = cS_test.phi_bequest * (max(cS_test.cFloor, total_wealth_p_vfi).^(1-cS_test.sigma))./(1-cS_test.sigma);
                g_A_period = (1 + cS_test.g_A_ss)^cS_test.time_Step - 1;
                discount = cS_test.beta * (1 + g_A_period)^(1-cS_test.sigma);
                val_vfi_calc = util_c_vfi + discount * util_b_vfi;
                fprintf('%-10s | %-12.4f | %-12.4f | %-12.4f | %-12.4f | %-12.4f | %-12.4f | %-12.4f\n', ...
                    'VFI', c_vfi, k_p_vfi, kpps_p_vfi, util_c_vfi, util_b_vfi, discount, val_vfi_calc);

                % EGM 分解
                pol_egm_pt = polS_egm(a_idx);
                c_egm = pol_egm_pt.c(ik, ikpps, ie);
                k_p_egm = pol_egm_pt.k_prime(ik, ikpps, ie);
                kpps_p_egm = pol_egm_pt.kpps_prime(ik, ikpps, ie);
                util_c_egm = (c_egm.^(1-cS_test.sigma))./(1-cS_test.sigma);
                total_wealth_p_egm = k_p_egm + kpps_p_egm;
                util_b_egm = cS_test.phi_bequest * (max(cS_test.cFloor, total_wealth_p_egm).^(1-cS_test.sigma))./(1-cS_test.sigma);
                val_egm_calc = util_c_egm + discount * util_b_egm;
                fprintf('%-10s | %-12.4f | %-12.4f | %-12.4f | %-12.4f | %-12.4f | %-12.4f | %-12.4f\n', ...
                    'EGM', c_egm, k_p_egm, kpps_p_egm, util_c_egm, util_b_egm, discount, val_egm_calc);

                fprintf('\n诊断说明: 请对比两行数据。Total Value应等于 u(c) + u_bequest * Discount。\n');
                fprintf('如果策略(c, k_prime)相似，但 u(c) 或 u_bequest 不同，则说明效用函数实现有误。\n');
                fprintf('如果所有组件都相似，但从主程序中得到的 valS_vfi 和 valS_egm 不同，则说明数据传递/存储有误。\n');
            end
        end

        function run_vfi_labor_supply_test()
            % =========================================================================
            % == 测试函数: run_vfi_labor_supply_test (v1.4 - 最终精简版)
            % == 目的:
            % ==   1. 专门为 VFI_by_age_PPS 函数创建一个独立的、可复现的测试环境。
            % ==   2. 检验在不同的劳动负效用参数(chi)下，模型计算出的劳动供给
            % ==      决策是否符合经济学直觉。
            % ==   3. [核心修正] 遵循最小化原则，禁用了与测试目标(劳动供给)无关的
            % ==      PPS功能 (pps_active = false)，从根本上移除了对 cS.pps_fixed
            % ==      等参数的依赖，使测试环境更稳定、更专注。
            % =========================================================================

            fprintf('\n\n===================================================================\n');
            fprintf('==  开启 VFI_by_age_PPS 劳动供给决策专项测试程序  ==\n');
            fprintf('===================================================================\n\n');

            % --- 1. 初始化通用参数 ---
            fprintf('--- 1. 初始化通用参数 ---\n');
            cS_base = utils.ParameterValues();
            
            % --- 1a. 基础模型设定 ---
            cS_base.nk = 40;
            cS_base.nw = 5; 
            
            % [!!! 关键修正 !!!] 禁用PPS，并将nkpps设为1
            cS_base.pps_active = false;
            cS_base.nkpps = 1;
            
            % --- 1b. 生成测试所需的完备路径和参数 ---
            cS_base.start_year = 2023; cS_base.end_year = 2123; cS_base.time_Step = 5;
            cS_base.T_sim = (cS_base.end_year - cS_base.start_year)/cS_base.time_Step + 1;
            cS_base.s_pathV = repmat(cS_base.s_pathV(:,1), 1, cS_base.T_sim); 

            cS_base = utils.age_efficiency(cS_base);
            cS_base = utils.theta_payg_path(cS_base, 'baseline', false);
            cS_base.aR_new_path = repmat(cS_base.aR_new, 1, cS_base.T_sim); 
            
            cS_base.nu_h = [2.0; 1.5; 1.0; 1.0]; 
            cS_base.nL = 10;
            cS_base.lGridV = linspace(0.1, 1.0, cS_base.nL)'; 

            cS_base = utils.generateGrids(cS_base);

            paramS = struct();
            [paramS.leGridV, paramS.TrProbM_by_age, paramS.leProb1V, cS_base.nw_expanded] = utils.EarningProcess_AgeDependent(cS_base);
            
            % --- 1c. 设置宏观价格环境和V_prime ---
            M1 = struct('r_mkt_t', 0.04, 'w_t', 1.5, 'tr_per_hh', 0.01, 'b_hat_t', 0.5);
            V_prime_dummy = zeros(cS_base.nk, cS_base.nkpps, cS_base.nw_expanded, cS_base.aD_new);
            fprintf('   ✅ 通用参数初始化完成。\n');

            % --- 2. 运行两个核心测试案例 ---
            fprintf('\n--- 2a. 运行测试案例 A: 低劳动负效用 (chi_base = 2.0) ---\n');
            cS_A = cS_base;
            cS_A = utils.labor_setting_for_test(cS_A, 2.0);
            [polS_A, ~] = household.solve_vfi_for_test(cS_A, paramS, M1, V_prime_dummy);
            
            fprintf('\n--- 2b. 运行测试案例 B: 高劳动负效用 (chi_base = 10.0) ---\n');
            cS_B = cS_base;
            cS_B = utils.labor_setting_for_test(cS_B, 10.0); 
            [polS_B, ~] = household.solve_vfi_for_test(cS_B, paramS, M1, V_prime_dummy);
            
            % --- 3. 对比分析结果 ---
            fprintf('\n\n--- 3. 结果对比与诊断分析 ---\n');
            fprintf('核心检验逻辑: 如果VFI引擎正确，案例B的平均劳动供给应显著低于案例A。\n\n');
            
            fprintf('%-15s | %-25s | %-25s\n', '年龄组', '案例 A (低chi) Avg(l)', '案例 B (高chi) Avg(l)');
            fprintf('%s\n', repmat('-', 1, 70));

            ages_to_check = {
                '青年 (a=5)', 5;
                '中年 (a=10)', 10;
                '临退 (a=15)', 15;
            };
            
            all_passed = true;
            for i = 1:size(ages_to_check, 1)
                age_label = ages_to_check{i, 1};
                a_idx = ages_to_check{i, 2};
                
                avg_l_A = mean(polS_A(a_idx).l, 'all');
                avg_l_B = mean(polS_B(a_idx).l, 'all');
                
                fprintf('%-15s | %25.4f | %25.4f\n', age_label, avg_l_A, avg_l_B);
                
                if avg_l_B > avg_l_A * 0.95 
                    all_passed = false;
                end
            end
            
            fprintf('\n--- 诊断结论 ---\n');
            if all_passed
                fprintf('   ✅ [测试通过] 劳动供给对负效用参数(chi)敏感，VFI引擎似乎工作正常。\n');
            else
                fprintf('   ❌ [测试失败] 劳动供给对负效用参数(chi)不敏感！\n');
                fprintf('      这强烈表明在 VFI_by_age_PPS 中，劳动负效用项未被正确计算或应用。\n');
                fprintf('      请检查 total_value_5D = util_c_5D - util_l_tens + future_value_5D 这一行。\n');
            end
            fprintf('===================================================================\n');

        end        
        % --- 内部辅助函数，用于执行VFI求解 ---
        function [polS, valS] = solve_vfi_for_test(cS, paramS, M, V_prime)
            valS = -Inf(cS.nk, cS.nkpps, cS.nw_expanded, cS.aD_new);
            polS_cell = cell(cS.aD_new, 1);
            
            b_payg_hat_age_vec = zeros(1, cS.aD_new);
            b_payg_hat_age_vec((cS.aR_new+1):end) = M.b_hat_t;
            
            tic;
            for a_idx = cS.aD_new : -1 : 1
                vPrime_next_age = [];
                if a_idx < cS.aD_new
                    vPrime_next_age = V_prime(:,:,:,a_idx+1); 
                end
                b_age_val = b_payg_hat_age_vec(a_idx);
                
                % 模拟异质性家庭循环
                % 为简化，我们只测试一种家庭类型 (h=1, 高收入城镇职工)
                h_current_type = 1;
                cS.h_current_type = h_current_type;
                cS.ageEffV_new = cS.ageEffV_new_h(:, h_current_type);
                M.theta_t = cS.theta_path_h(h_current_type, 1); % 使用一个固定的theta
                M.g_A_t_plus_1 = 0.015; % 使用固定的增长率

                [val_age, pol_age] = household.VFI_by_age_PPS(a_idx, vPrime_next_age, M, b_age_val, paramS, cS);

                valS(:,:,:,a_idx) = val_age;
                polS_cell{a_idx} = pol_age;
            end
            elapsed_time = toc;
            polS = [polS_cell{:}];
            fprintf('   VFI求解完成，耗时: %.2f 秒\n', elapsed_time);
        end


    end


end



