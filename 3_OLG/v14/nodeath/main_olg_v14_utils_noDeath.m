% --- START OF FILE main_olg_v14_utils_noDeath.m (Corrected v3 - Final) ---

classdef main_olg_v14_utils_noDeath

    methods (Static)

        % =========================================================================
        % == 1. 核心求解与模拟入口 (新)
        % =========================================================================
        function [K_next_agg, C_agg] = solve_and_simulate_onestep(M_t, cS, paramS, eIdxM, k_initial_dist)
            fprintf('   (noDeath) 正在求解无限期VFI...\n');
            [cPolM, kPolM, ~] = main_olg_v14_utils_noDeath.HHSolution_VFI_Huggett_noDeath(M_t, paramS, cS);
            fprintf('   (noDeath) VFI求解完成。\n');
            
            fprintf('   (noDeath) 正在模拟家庭决策...\n');
            [K_next_agg, C_agg] = main_olg_v14_utils_noDeath.HHSimulation_noDeath(kPolM, cPolM, eIdxM, M_t, paramS, cS, k_initial_dist);
            fprintf('   (noDeath) 模拟完成。\n');
        end


        % =========================================================================
        % == 2. 家庭问题求解函数 (Aiyagari/Huggett 风格)
        % =========================================================================

        function [cPolM, kPolM, valM_converged] = HHSolution_VFI_Huggett_noDeath(M_vfi, paramS_vfi, cS_vfi)
            max_iter = 500; tol = 1e-6; dist = 1e6; iter = 0;
            valM_old = zeros(cS_vfi.nk, cS_vfi.nw);
            while dist > tol && iter < max_iter
                iter = iter + 1;
                [cPolM, kPolM, valM_new] = main_olg_v14_utils_noDeath.HHSolutionByPeriod_VFI_GridSearch_noDeath(valM_old, M_vfi, paramS_vfi, cS_vfi);
                dist = max(abs(valM_new(:) - valM_old(:)));
                valM_old = valM_new;
                if mod(iter, 25) == 0, fprintf('      VFI (noDeath) 迭代 %d, 距离 = %.4e\n', iter, dist); end
            end
            if iter == max_iter, warning('VFI (noDeath) 在 %d 次迭代内未收敛。', max_iter); end
            valM_converged = valM_new;
        end

        function [cPol_period, kPol_period, val_period] = HHSolutionByPeriod_VFI_GridSearch_noDeath(vPrime_ke_next, M_age, paramS_age, cS)
            % [审计修正版 v3]
            % 核心修正: 预算约束使用税前口径，与宏观审计保持一致。
            
            val_period    = -Inf(cS.nk, cS.nw);
            cPol_period   = zeros(cS.nk, cS.nw);
            kPol_period   = zeros(cS.nk, cS.nw);

            EV_matrix = zeros(cS.nk, cS.nw);
            for ie_current = 1:cS.nw
                transition_probs = paramS_age.leTrProbM(ie_current, :);
                EV_slice = vPrime_ke_next * transition_probs';
                EV_matrix(:, ie_current) = EV_slice;
            end
            
            ev_interpolants = cell(cS.nw, 1);
            for ie_current = 1:cS.nw, ev_interpolants{ie_current} = griddedInterpolant(cS.kGridV, EV_matrix(:, ie_current), 'linear', 'linear'); end

            for ik = 1:cS.nk
                for ie = 1:cS.nw
                    k_state = cS.kGridV(ik);
                    epsilon_state = paramS_age.leGridV(ie);
                    ev_interpolant = ev_interpolants{ie};
                    best_val = -Inf; best_c = cS.cFloor; best_k_prime = cS.kMin;
                    
                    % [核心修正] 1. 计算税前总资源
                    labor_income_gross = M_age.w_t * epsilon_state;
                    capital_income_net = k_state * M_age.r_net_period;
                    total_sources = k_state + capital_income_net + labor_income_gross; % 税前口径

                    % [核心修正] 2. 计算总支出
                    labor_tax = cS.tau_l * labor_income_gross;
                    
                    % 可用于消费和储蓄的净资源
                    disposable_resources = total_sources - labor_tax;

                    k_prime_max_budget = disposable_resources - cS.cFloor * (1 + cS.tau_c);
                    if k_prime_max_budget < cS.kMin, continue; end

                    k_prime_grid = linspace(cS.kMin, k_prime_max_budget, cS.nkprime);
                    for k_prime_choice = k_prime_grid
                        c_expend = disposable_resources - k_prime_choice;
                        c_choice = max(cS.cFloor, c_expend / (1 + cS.tau_c));
                        [~, util] = main_olg_v14_utils.CES_utility(c_choice, cS.sigma, cS);
                        k_prime_clamped = max(cS.kGridV(1), min(cS.kGridV(end), k_prime_choice));
                        ev = ev_interpolant(k_prime_clamped);
                        current_val = util + cS.beta * ev;
                        if current_val > best_val
                            best_val = current_val; best_c = c_choice; best_k_prime = k_prime_choice;
                        end
                    end
                    val_period(ik, ie) = best_val; cPol_period(ik, ie) = best_c; kPol_period(ik, ie) = best_k_prime;
                end
            end
        end

        % =========================================================================
        % == 3. 模拟与聚合函数 (Aiyagari/Huggett 风格)
        % =========================================================================
        function [K_next_agg, C_agg] = HHSimulation_noDeath(kPolM, cPolM, eIdxM, M_sim, paramS_sim, cS_sim, k_initial_dist)
            % [审计修正版 v3]
            % 核心修正: 预算约束使用税前口径，与宏观审计和VFI求解器保持一致。

            nSim = size(eIdxM, 1);
            kPolInterp = cell(cS_sim.nw, 1);
            for ie = 1:cS_sim.nw
                if cS_sim.nk > 1, kPolInterp{ie} = griddedInterpolant(cS_sim.kGridV, kPolM(:, ie), 'linear', 'linear');
                else, kPolInterp{ie} = @(k) kPolM(1, ie); end
            end

            k_now = k_initial_dist; e_idx_now = eIdxM(:, 1);
            k_next = zeros(nSim, 1);

            % 1. 模拟下一期资本 k_next
            for ie = 1:cS_sim.nw
                idx_sim = find(e_idx_now == ie);
                if isempty(idx_sim), continue; end
                k_next(idx_sim) = kPolInterp{ie}(k_now(idx_sim));
            end
            
            % 2. 使用预算约束反推当期消费 c_decision
            epsilon_state = paramS_sim.leGridV(e_idx_now);
            labor_income_gross = M_sim.w_t .* epsilon_state;
            capital_income_net = k_now .* M_sim.r_net_period;
            
            % [核心修正] 使用税前口径计算总资源
            total_sources = k_now + capital_income_net + labor_income_gross;
            
            % [核心修正] 将劳动税作为一项支出
            labor_tax = cS_sim.tau_l .* labor_income_gross;
            
            % 消费支出 = (总资源 - 劳动税) - 为下一期储蓄的资本
            c_expenditure = (total_sources - labor_tax) - k_next;
            
            c_decision = c_expenditure ./ (1 + cS_sim.tau_c);
            
            % --- 聚合 ---
            K_next_agg = mean(k_next);
            C_agg = mean(c_decision);
        end

    end % End of Static Methods
end