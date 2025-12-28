classdef main_olg_v6_utils
    % OLG 模型 v6 的工具函数 (VFI包含k_pps作为状态变量, PPS所得税递延)

    methods (Static)

        function cS = ParameterValues_HuggettStyle()
            cS = main_olg_v4_utils_pps.ParameterValues_HuggettStyle();

            cS.tau_l_init_guess = 0.05;
            cS.tau_l_min = 0.00;
            cS.tau_l_max = 0.3; % 调整一般劳动税上限
            cS.max_total_labor_tax = 0.6;

            cS.theta_payg_max = 0.35;

            cS.pps_active = true;
            cS.pps_annual_contrib_limit = 1.5;
            cS.pps_max_contrib_frac = 0.15;
            cS.pps_in_K = true;
            cS.pps_bequeathable = false;
            cS.pps_tax_rate_withdrawal = 0.03;
            cS.pps_return_rate_premium = 0.1;
            cS.pps_withdrawal_rate = 0.15;
            cS.pps_contribution_age_max_idx = cS.aR_idx_orig - 1;
            cS.pps_withdrawal_age_min_idx = cS.aR_idx_orig;

            cS.nkpps = 20;
            cS.kppsMin = 0;
            cS.kppsMax = cS.kMax / 2; % 示例：PPS资产上限为非PPS的一半
            if cS.nkpps > 0 % Ensure kppsMax is at least a small positive if nkpps > 0
                cS.kppsMax = max(cS.kppsMax, 1e-3);
            end
            power_kpps = 1.5;
            if cS.nkpps > 1
                kppsGridV_temp = cS.kppsMin + (cS.kppsMax - cS.kppsMin) * (linspace(0, 1, cS.nkpps).^power_kpps);
                kppsGridV_temp(1) = cS.kppsMin;
            elseif cS.nkpps == 1
                kppsGridV_temp = cS.kppsMin; % Or an average value if preferred for single point
            else % nkpps = 0
                kppsGridV_temp = [];
            end
            cS.kppsGridV = kppsGridV_temp(:);

            cS.max_iter_K_tau_l = 100;
            cS.tol_K_tau_l = 1e-4;
            cS.damp_K_v5 = 0.3;
            cS.damp_tau_l_v5 = 0.3;
            cS.gbc_tol_for_internal_loop = 1e-3;
            cS.gbc_tol_for_rho_search = 1e-2;

            cS.max_stagnation_iters = 10;
            cS.min_norm_improvement_frac = 1e-3;
            cS.max_tau_l_boundary_strikes = 5;
        end

        % --- 人口和基础经济函数 (不变) ---
        function popS = initPopulation(cS); popS = main_olg_v4_utils_pps.initPopulation(cS); end
        function popS = populationDynamics(popS, cS); popS = main_olg_v4_utils_pps.populationDynamics(popS, cS); end
        function [Z_ss, dep_ratio_ss, bgp_reached, bgp_period] = detectSteadyStatePopulation(popS, cS); [Z_ss, dep_ratio_ss, bgp_reached, bgp_period] = main_olg_v4_utils_pps.detectSteadyStatePopulation(popS, cS); end
        function [logGridV, trProbM, prob1V] = EarningProcess_olgm(cS); [logGridV, trProbM, prob1V] = main_olg_v4_utils_pps.EarningProcess_olgm(cS); end
        function eIdxM = LaborEndowSimulation_olgm(cS, paramS); eIdxM = main_olg_v4_utils_pps.LaborEndowSimulation_olgm(cS, paramS); end
        function [HHlaborM_group, L] = LaborSupply_Huggett(eIdxM, cS, paramS, Z_ss_norm_group); [HHlaborM_group, L] = main_olg_v4_utils_pps.LaborSupply_Huggett(eIdxM, cS, paramS, Z_ss_norm_group); end
        function [R_market_gross_factor, MPL_gross] = HHPrices_Huggett(K_productive, L_total_eff, cS); [R_market_gross_factor, MPL_gross] = main_olg_v4_utils_pps.HHPrices_Huggett(K_productive, L_total_eff, cS); end
        function [muM, utilM] = CES_utility(cM_quantity, sig, cS); [muM, utilM] = main_olg_v4_utils_pps.CES_utility(cM_quantity, sig, cS); end

        % --- V6 家庭收入计算 (PPS所得税递延) ---
        function [resources_for_c_and_k_prime, labor_income_gross_state, pps_deduction_actual_state] = HHIncome_Huggett(k_now_val, R_k_net, w_gross, ...
                TR_total, b_payg_val, c_pps_choice_val, ...
                a_idx, paramS_hh, cS, epsilon_val)
            labor_income_gross_state = 0; pps_deduction_actual_state = 0; non_capital_income = 0;
            if a_idx <= cS.aR_new
                age_efficiency = cS.ageEffV_new(a_idx);
                labor_income_gross_state = w_gross * age_efficiency * epsilon_val;
                c_pps_choice_val = max(0, c_pps_choice_val);
                max_pps_by_frac = labor_income_gross_state * cS.pps_max_contrib_frac;
                pps_deduction_for_limit_check = min(c_pps_choice_val, cS.pps_annual_contrib_limit);
                pps_deduction_for_limit_check = min(pps_deduction_for_limit_check, max_pps_by_frac);

                if isfield(paramS_hh, 'pps_tax_deferral_active') && paramS_hh.pps_tax_deferral_active
                    pps_deduction_actual_state = pps_deduction_for_limit_check;
                else
                    pps_deduction_actual_state = 0; % No tax deferral if flag is false or missing
                end

                labor_income_taxable_for_tau_l = labor_income_gross_state - pps_deduction_actual_state;
                labor_income_taxable_for_tau_l = max(0, labor_income_taxable_for_tau_l);

                income_tax_tau_l = labor_income_taxable_for_tau_l * paramS_hh.tau_l;
                payg_tax_theta = labor_income_gross_state * paramS_hh.theta_payg_actual_for_hh;

                labor_income_net_of_all_taxes = labor_income_gross_state - income_tax_tau_l - payg_tax_theta;
                non_capital_income = labor_income_net_of_all_taxes + TR_total + b_payg_val;
            else
                non_capital_income = TR_total + b_payg_val;
            end
            capital_income = R_k_net * k_now_val;
            resources_for_c_and_k_prime = capital_income + non_capital_income - c_pps_choice_val;
            if ~isfinite(resources_for_c_and_k_prime), resources_for_c_and_k_prime = -1e10; end
        end

        % --- V6 VFI 主函数 (状态空间扩展为 k, k_pps, epsilon) ---
        function [cPolM_q, kPolM, cPpsPolM, valM] = HHSolution_VFI_Huggett(R_k_net, w_gross, TR_tot, bV_payg, paramS_vfi, cS)
            cPolM_q  = zeros(cS.nk, cS.nkpps, cS.nw, cS.aD_new); kPolM  = zeros(cS.nk, cS.nkpps, cS.nw, cS.aD_new);
            cPpsPolM = zeros(cS.nk, cS.nkpps, cS.nw, cS.aD_new); valM = -Inf(cS.nk, cS.nkpps, cS.nw, cS.aD_new); % Initialize valM to -Inf

            % --- 确保输出变量在函数开始时就被定义和预分配 ---
            cPol_age_q = zeros(cS.nk, cS.nkpps, cS.nw);
            kPol_age   = zeros(cS.nk, cS.nkpps, cS.nw);
            cPpsPol_age= zeros(cS.nk, cS.nkpps, cS.nw);
            val_age    = -Inf(cS.nk, cS.nkpps, cS.nw); % Default to -Inf

            fprintf('  VFI (HHSolution_VFI_Huggett): Starting backwards induction...\n');
            total_age_groups_vfi = cS.aD_new;
            for a_idx = cS.aD_new : -1 : 1
                if mod(cS.aD_new - a_idx + 1, 1) == 0 || a_idx == cS.aD_new || a_idx == 1
                    progress_pct_vfi = (cS.aD_new - a_idx + 1) / total_age_groups_vfi * 100;
                    fprintf('    VFI for age group %2d / %2d (approx %.0f%%)...\n', a_idx, total_age_groups_vfi, progress_pct_vfi);
                end
                vPrime_kkppse_next = [];
                if a_idx < cS.aD_new, vPrime_kkppse_next = valM(:,:,:,a_idx+1); end
                eps_grid = paramS_vfi.leGridV;
                [cPolM_q(:,:,:,a_idx), kPolM(:,:,:,a_idx), cPpsPolM(:,:,:,a_idx), valM(:,:,:,a_idx)] = ...
                    main_olg_v6_utils.HHSolutionByAge_VFI_Huggett(a_idx, vPrime_kkppse_next, ...
                    R_k_net, w_gross, TR_tot, bV_payg(a_idx), paramS_vfi, cS, eps_grid);
            end
            fprintf('  VFI (HHSolution_VFI_Huggett): Finished.\n');
        end

        % --- V6 VFI 按年龄组求解 (状态 k, k_pps, epsilon; 优化 k', c_pps) ---
        function [cPol_age_q, kPol_age, cPpsPol_age, val_age] = HHSolutionByAge_VFI_Huggett(...
                a_idx, vPrime_kkppse_next, R_k_net_age, w_gross_age, TR_total_age, b_age_val, paramS_age, cS, epsilon_grid)

            % --- 确保输出变量在函数开始时就被定义和预分配 ---
            cPol_age_q = zeros(cS.nk, cS.nkpps, cS.nw);
            kPol_age   = zeros(cS.nk, cS.nkpps, cS.nw);
            cPpsPol_age= zeros(cS.nk, cS.nkpps, cS.nw);
            val_age    = -Inf(cS.nk, cS.nkpps, cS.nw); % Default to -Inf

            model_age_group_start_year_idx = cS.physAgeMap{a_idx}(1);
            is_pps_contrib_eligible = (model_age_group_start_year_idx <= cS.pps_contribution_age_max_idx && ...
                cS.pps_active && (cS.pps_max_contrib_frac > 0 || cS.pps_annual_contrib_limit > 0) );

            fmincon_opts = optimoptions('fmincon','Display','none','Algorithm','sqp', 'TolFun', 1e-7, 'TolX', 1e-7, 'MaxFunctionEvaluations', 3000, 'MaxIterations', 750, 'StepTolerance', 1e-7);
            fminbnd_opts = optimset('TolX', 1e-6, 'Display', 'none');

            if a_idx == cS.aD_new
                % Last period logic directly assigns to output variables
                for ik_temp = 1:cS.nk
                    for ikpps_temp = 1:cS.nkpps
                        k_now_val_last = cS.kGridV(ik_temp);
                        k_pps_now_val_last = cS.kppsGridV(ikpps_temp);
                        for ie_temp = 1:cS.nw
                            [resources_non_pps, ~, ~] = main_olg_v6_utils.HHIncome_Huggett(...
                                k_now_val_last, R_k_net_age, w_gross_age, TR_total_age, b_age_val, ...
                                0, a_idx, paramS_age, cS, epsilon_grid(ie_temp));

                            pps_final_withdrawal_pretax = k_pps_now_val_last;
                            pps_final_withdrawal_net = pps_final_withdrawal_pretax * (1 - cS.pps_tax_rate_withdrawal);
                            total_resources_for_consumption = resources_non_pps + pps_final_withdrawal_net;

                            cPol_age_q(ik_temp,ikpps_temp,ie_temp) = max(cS.cFloor, total_resources_for_consumption / (1 + cS.tau_c) );
                            kPol_age(ik_temp,ikpps_temp,ie_temp) = cS.kMin;
                            cPpsPol_age(ik_temp,ikpps_temp,ie_temp) = 0;
                            [~, val_age(ik_temp,ikpps_temp,ie_temp)] = main_olg_v6_utils.CES_utility(cPol_age_q(ik_temp,ikpps_temp,ie_temp), cS.sigma, cS);
                        end; end; end
            else % Not the last period
                EV_interpolants = cell(cS.nw, 1);
                for ie_current = 1:cS.nw
                    EV_for_interp = zeros(cS.nk, cS.nkpps);
                    for ik_next = 1:cS.nk
                        for ikpps_next = 1:cS.nkpps
                            expected_v_sum = 0;
                            for ie_next = 1:cS.nw
                                expected_v_sum = expected_v_sum + vPrime_kkppse_next(ik_next, ikpps_next, ie_next) * paramS_age.leTrProbM(ie_current, ie_next); 
                            end
                            EV_for_interp(ik_next, ikpps_next) = expected_v_sum; 
                        end
                    end

                    if cS.nk > 1 && cS.nkpps > 1
                        EV_interpolants{ie_current} = griddedInterpolant({cS.kGridV, cS.kppsGridV}, EV_for_interp, 'linear');
                    elseif cS.nk > 1, EV_interpolants{ie_current} = griddedInterpolant(cS.kGridV, EV_for_interp(:,1), 'linear');
                    elseif cS.nkpps > 1, EV_interpolants{ie_current} = griddedInterpolant(cS.kppsGridV, EV_for_interp(1,:)', 'linear');
                    else, EV_interpolants{ie_current} = @(k_s, kp_s) EV_for_interp(1,1); 
                    end
                end



                % Temporary variables for parfor loop slices
                c_results_cell = cell(cS.nk,1);
                k_results_cell = cell(cS.nk,1);
                cpps_results_cell = cell(cS.nk,1);
                v_results_cell = cell(cS.nk,1);

                for ik = 1:cS.nk % Re-enabled parfor
                    c_slice_par = zeros(cS.nkpps, cS.nw);
                    k_slice_par = zeros(cS.nkpps, cS.nw);
                    cpps_slice_par = zeros(cS.nkpps, cS.nw);
                    v_slice_par = -Inf(cS.nkpps, cS.nw);

                    % These options are broadcast to workers
                    fmincon_opts_pf = fmincon_opts;
                    fminbnd_opts_pf = fminbnd_opts;

                    for ikpps = 1:cS.nkpps
                        for ie = 1 : cS.nw
                            k_state = cS.kGridV(ik); k_pps_state = cS.kppsGridV(ikpps); epsilon_state = epsilon_grid(ie);
                            EV_interpolant_for_state_pf = EV_interpolants{ie};
                            obj_fun = @(x_choices) -main_olg_v6_utils.BellmanForKprimeKppsCpps(x_choices(1), x_choices(2), k_state, k_pps_state, epsilon_state, a_idx, R_k_net_age, w_gross_age, TR_total_age, b_age_val, EV_interpolant_for_state_pf, cS, paramS_age);
                            lab_inc_g_state = 0; 
                            if a_idx <= cS.aR_new
                                age_eff_state = cS.ageEffV_new(a_idx);
                                lab_inc_g_state = w_gross_age * age_eff_state * epsilon_state; 
                            end
                            max_cpps_this_state = 0; 
                            if is_pps_contrib_eligible && lab_inc_g_state > 0
                                max_cpps_by_frac_state = lab_inc_g_state * cS.pps_max_contrib_frac; 
                                max_cpps_this_state = min(max_cpps_by_frac_state, cS.pps_annual_contrib_limit); 
                                max_cpps_this_state = max(0, max_cpps_this_state); 
                            end
                            lb_optim = [0; cS.kMin]; ub_optim = [max_cpps_this_state; cS.kMax];
                            cpps_guess_init = 0; 
                            if is_pps_contrib_eligible && max_cpps_this_state > 0
                                cpps_guess_init = 0.05 * max_cpps_this_state; 
                            end
                            k_guess_init = k_state;
                            x0_optim = [cpps_guess_init; k_guess_init]; x0_optim(1) = max(lb_optim(1), min(ub_optim(1), x0_optim(1))); x0_optim(2) = max(lb_optim(2), min(ub_optim(2), x0_optim(2)));
                            [temp_res_for_k_opt_fallback, ~, ~] = main_olg_v6_utils.HHIncome_Huggett(k_state,R_k_net_age,w_gross_age,TR_total_age,b_age_val,0,a_idx,paramS_age,cS,epsilon_state);
                            if ~is_pps_contrib_eligible || max_cpps_this_state < 1e-7 
                                cpps_val = 0; 
                                ev_1d_for_fallback = @(kprime_val) main_olg_v6_utils.CallInterpolator(EV_interpolant_for_state_pf, kprime_val, cS.kppsGridV(1), cS);
                                if cS.nkpps == 1 && isa(EV_interpolant_for_state_pf, 'griddedInterpolant')
                                    ev_1d_for_fallback = EV_interpolant_for_state_pf; 
                                elseif cS.nkpps == 1 && isa(EV_interpolant_for_state_pf, 'function_handle')
                                    ev_1d_for_fallback = EV_interpolant_for_state_pf; 
                                end
                                [c_val, k_p_val, v_val] = main_olg_v6_utils.HHSolutionByOneState_OptK_Mod(a_idx, temp_res_for_k_opt_fallback, ev_1d_for_fallback, fminbnd_opts_pf, cS, paramS_age);
                            else
                                try
                                    nonlcon_fun = @(x_nlc) main_olg_v6_utils.fmincon_nonlcon_v6(x_nlc, k_state, k_pps_state, epsilon_state, a_idx, R_k_net_age, w_gross_age, TR_total_age, b_age_val, cS, paramS_age);
                                    [x_opt, neg_v_opt, exitflag, output_fmincon] = fmincon(obj_fun, x0_optim, [],[],[],[], lb_optim, ub_optim, nonlcon_fun, fmincon_opts_pf);
                                    
                                if exitflag <= 0 
                                    fprintf('DEBUG: fmincon failed for state (a=%d, k_idx=%d, kpps_idx=%d, e_idx=%d) with exitflag=%d\n', a_idx, ik, ikpps, ie, exitflag);
                                    disp(output_fmincon); % 打印fmincon的输出结构体
                                    error('fmincon_failed_exitflag_v6'); 
                                end
                                    cpps_val = x_opt(1);
                                    k_p_val = x_opt(2);
                                    v_val = -neg_v_opt;
                                    [res_final, ~, ~] = main_olg_v6_utils.HHIncome_Huggett(k_state,R_k_net_age,w_gross_age,TR_total_age,b_age_val,cpps_val,a_idx,paramS_age,cS,epsilon_state);
                                    c_val = (res_final - k_p_val) / (1+cS.tau_c);
                                    c_val = max(cS.cFloor, c_val);
                                catch
                                    warning('1-d interp')
                                    ME_fmincon_v6, cpps_val = 0;
                                    ev_1d_for_fallback = @(kprime_val) main_olg_v6_utils.CallInterpolator(EV_interpolant_for_state_pf, kprime_val, cS.kppsGridV(1), cS);
                                    if cS.nkpps == 1 && isa(EV_interpolant_for_state_pf, 'griddedInterpolant'), ev_1d_for_fallback = EV_interpolant_for_state_pf; elseif cS.nkpps == 1 && isa(EV_interpolant_for_state_pf, 'function_handle'), ev_1d_for_fallback = EV_interpolant_for_state_pf; end; [c_val, k_p_val, v_val] = main_olg_v6_utils.HHSolutionByOneState_OptK_Mod(a_idx, temp_res_for_k_opt_fallback, ev_1d_for_fallback, fminbnd_opts_pf, cS, paramS_age); end
                            end; c_slice_par(ikpps, ie)=c_val; k_slice_par(ikpps, ie)=k_p_val; cpps_slice_par(ikpps, ie)=cpps_val; v_slice_par(ikpps, ie)=v_val;
                        end; end % End ie, ikpps loops
                    c_results_cell{ik} = c_slice_par; k_results_cell{ik} = k_slice_par;
                    cpps_results_cell{ik} = cpps_slice_par; v_results_cell{ik} = v_slice_par;
                end % End ik parfor loop

                % Assemble results from cell arrays directly into output variables
                for ik_assemble = 1:cS.nk
                    if ~isempty(c_results_cell{ik_assemble})
                        val_age(ik_assemble,:,:) = v_results_cell{ik_assemble};
                        cPol_age_q(ik_assemble,:,:) = c_results_cell{ik_assemble};
                        kPol_age(ik_assemble,:,:) = k_results_cell{ik_assemble};
                        cPpsPol_age(ik_assemble,:,:) = cpps_results_cell{ik_assemble};
                    else
                        warning('HHSolutionByAge: Empty cell result for ik_assemble = %d, age_idx = %d. Value function may be incorrect.', ik_assemble, a_idx);
                        % val_age already initialized to -Inf, others to 0
                    end
                end
            end % End if a_idx == cS.aD_new / else
        end % End HHSolutionByAge_VFI_Huggett

        function V_out = BellmanForKprimeKppsCpps(c_pps_choice, k_prime_choice, k_now, k_pps_now, epsilon_now, a_idx_bellman, R_k_net, w_gross, TR_tot, b_payg_val, EV_next_interpolant_2D, cS_bellman, paramS_bellman)
            c_pps_choice = max(0, c_pps_choice); k_prime_choice = max(cS_bellman.kMin, min(cS_bellman.kMax, k_prime_choice));
            [resources_for_c_kprime, ~, ~] = main_olg_v6_utils.HHIncome_Huggett(k_now, R_k_net, w_gross, TR_tot, b_payg_val, c_pps_choice, a_idx_bellman, paramS_bellman, cS_bellman, epsilon_now);
            consump_budget_expenditure = resources_for_c_kprime - k_prime_choice; min_consump_expenditure = cS_bellman.cFloor * (1 + cS_bellman.tau_c);
            if consump_budget_expenditure < min_consump_expenditure - 1e-6, V_out = -1e20 - 1e10 * (min_consump_expenditure - consump_budget_expenditure); return; end
            consump_quantity = consump_budget_expenditure / (1 + cS_bellman.tau_c);
            [~, util_current] = main_olg_v6_utils.CES_utility(consump_quantity, cS_bellman.sigma, cS_bellman);
            if a_idx_bellman < cS_bellman.aD_new
                model_age_group_start_year_idx_bellman = cS_bellman.physAgeMap{a_idx_bellman}(1); pps_wd_pretax_this_period = 0; current_annual_model_age_idx = model_age_group_start_year_idx_bellman; is_retired_annual_idx_check = (current_annual_model_age_idx >= cS_bellman.aR_idx_orig);
                if is_retired_annual_idx_check && current_annual_model_age_idx >= cS_bellman.pps_withdrawal_age_min_idx && cS_bellman.pps_active, pps_wd_pretax_this_period = k_pps_now * cS_bellman.pps_withdrawal_rate; end
                pps_return_factor = 1 + ( (R_k_net - 1) + cS_bellman.pps_return_rate_premium ); k_pps_prime_val = (k_pps_now + c_pps_choice - pps_wd_pretax_this_period) * pps_return_factor; k_pps_prime_val = max(cS_bellman.kppsMin, min(cS_bellman.kppsMax, k_pps_prime_val));
                EV_next = main_olg_v6_utils.CallInterpolator(EV_next_interpolant_2D, k_prime_choice, k_pps_prime_val, cS_bellman); % Use CallInterpolator
                if ~isfinite(EV_next), EV_next = -1e12; end; s_transition = cS_bellman.s_1yr_transitionV(a_idx_bellman); V_out = util_current + cS_bellman.beta * s_transition * EV_next; else, V_out = util_current; end
            if ~isfinite(V_out), V_out = -1e12; end
        end

        function [c_ineq, ceq] = fmincon_nonlcon_v6(x_nlc, k_now_nlc, k_pps_now_nlc_unused, epsilon_nlc, a_idx_nlc, R_k_net_nlc, w_gross_nlc, TR_total_nlc, b_age_val_nlc, cS_nlc, paramS_nlc)
            ceq = []; c_pps_nlc = x_nlc(1); k_prime_nlc = x_nlc(2);
            [resources_for_c_kprime_nlc, ~, ~] = main_olg_v6_utils.HHIncome_Huggett(k_now_nlc, R_k_net_nlc, w_gross_nlc, TR_total_nlc, b_age_val_nlc, c_pps_nlc, a_idx_nlc, paramS_nlc, cS_nlc, epsilon_nlc);
            consumption_expenditure_nlc = resources_for_c_kprime_nlc - k_prime_nlc; min_consumption_expenditure_nlc = cS_nlc.cFloor * (1 + cS_nlc.tau_c);
            c_ineq = min_consumption_expenditure_nlc - consumption_expenditure_nlc;
        end

        function [c_quantity, kPrime, ValueFunc] = HHSolutionByOneState_OptK_Mod(a_new_osa_k, budget_for_c_expenditure_and_kprime, vPofK_int_osa_k_1D, fminbnd_opts_in, cS, paramS_osa_k)
            kPMin_osa_k = cS.kMin; kPMax_osa_k = budget_for_c_expenditure_and_kprime - (cS.cFloor * (1 + cS.tau_c)); ValueFunc = -Inf;
            function negV_nested = negBellmanObjective_nested(kP_nested)
                [negV_nested, ~] = BellmanInner_nested(kP_nested); end
            function [negVal_nested, Val_nested] = BellmanInner_nested(kP_inner_nested)
                cons_quantity_nested = max(cS.cFloor, (budget_for_c_expenditure_and_kprime - kP_inner_nested) / (1 + cS.tau_c) ); [~, util_val_nested] = main_olg_v6_utils.CES_utility(cons_quantity_nested, cS.sigma, cS); if ~isfinite(util_val_nested), util_val_nested = -1e12; end; Val_nested = util_val_nested; if a_new_osa_k < cS.aD_new, evF_val_nested = -Inf; try kP_eval_nested = max(cS.kGridV(1),min(cS.kGridV(end),kP_inner_nested)); evF_val_nested=vPofK_int_osa_k_1D(kP_eval_nested); catch; if kP_inner_nested<cS.kGridV(1),evF_val_nested=vPofK_int_osa_k_1D(cS.kGridV(1));else,evF_val_nested=vPofK_int_osa_k_1D(cS.kGridV(end));end;end; if ~isfinite(evF_val_nested),evF_val_nested=-1e12;end; s_transition_nested = cS.s_1yr_transitionV(a_new_osa_k); Val_nested = util_val_nested + cS.beta*s_transition_nested*evF_val_nested;end; if ~isfinite(Val_nested), negVal_nested=1e12; Val_nested=-1e12; else negVal_nested = -Val_nested; end; end
        if kPMax_osa_k <= kPMin_osa_k, kPrime = kPMin_osa_k; c_quantity = max(cS.cFloor, (budget_for_c_expenditure_and_kprime - kPrime) / (1 + cS.tau_c) ); [~, u_c_osa_k] = main_olg_v6_utils.CES_utility(c_quantity, cS.sigma, cS); if a_new_osa_k < cS.aD_new, EV_c_osa_k = -Inf; try kP_eval_osa_k = max(cS.kGridV(1),min(cS.kGridV(end),kPrime)); EV_c_osa_k=vPofK_int_osa_k_1D(kP_eval_osa_k); catch; if kPrime<cS.kGridV(1),EV_c_osa_k=vPofK_int_osa_k_1D(cS.kGridV(1));else,EV_c_osa_k=vPofK_int_osa_k_1D(cS.kGridV(end));end;end; if ~isfinite(EV_c_osa_k),EV_c_osa_k=-1e12;end; s_transition_osa_k = cS.s_1yr_transitionV(a_new_osa_k); ValueFunc = u_c_osa_k + cS.beta*s_transition_osa_k*EV_c_osa_k; else, ValueFunc = u_c_osa_k; end; if ~isfinite(ValueFunc), ValueFunc = -1e12; end
        else, kPMin_opt_osa_k = max(cS.kMin, kPMin_osa_k); kPMax_opt_osa_k = max(kPMin_opt_osa_k + 1e-9, min(cS.kMax, kPMax_osa_k)); if kPMin_opt_osa_k >= kPMax_opt_osa_k, kPrime_opt_osa_k = kPMin_opt_osa_k; [negV_osa_k, ~] = BellmanInner_nested(kPrime_opt_osa_k); ValueFunc = -negV_osa_k; else, obj_osa_k = @(kP_osa_k) negBellmanObjective_nested(kP_osa_k); [kPrime_opt_osa_k, negV_osa_k, eflag_osa_k] = fminbnd(obj_osa_k, kPMin_opt_osa_k, kPMax_opt_osa_k, fminbnd_opts_in); if eflag_osa_k <= 0 || abs(kPrime_opt_osa_k-kPMin_opt_osa_k)<1e-7 || abs(kPrime_opt_osa_k-kPMax_opt_osa_k)<1e-7, [nV_min_osa_k,~]=BellmanInner_nested(kPMin_opt_osa_k); [nV_max_osa_k,~]=BellmanInner_nested(kPMax_opt_osa_k); if nV_min_osa_k <= nV_max_osa_k + 1e-9 , kPrime_opt_osa_k=kPMin_opt_osa_k; negV_osa_k=nV_min_osa_k; else, kPrime_opt_osa_k=kPMax_opt_osa_k; negV_osa_k=nV_max_osa_k; end; end; ValueFunc = -negV_osa_k; end; kPrime = kPrime_opt_osa_k; c_quantity = max(cS.cFloor, (budget_for_c_expenditure_and_kprime - kPrime) / (1 + cS.tau_c) ); if ~isfinite(ValueFunc), ValueFunc = -1e12; end; end
        kPrime = max(cS.kMin, min(cS.kMax, kPrime)); c_quantity = max(cS.cFloor, (budget_for_c_expenditure_and_kprime - kPrime) / (1 + cS.tau_c) ); [~, u_final_c_osa_k] = main_olg_v6_utils.CES_utility(c_quantity, cS.sigma, cS); if a_new_osa_k < cS.aD_new, EV_A_final_osa_k = -Inf; try kP_eval_final_osa_k = max(cS.kGridV(1),min(cS.kGridV(end),kPrime)); EV_A_final_osa_k=vPofK_int_osa_k_1D(kP_eval_final_osa_k); catch; if kPrime<cS.kGridV(1),EV_A_final_osa_k=vPofK_int_osa_k_1D(cS.kGridV(1));else,EV_A_final_osa_k=vPofK_int_osa_k_1D(cS.kGridV(end));end;end; if ~isfinite(EV_A_final_osa_k),EV_A_final_osa_k=-1e12;end; s_transition_final_osa_k = cS.s_1yr_transitionV(a_new_osa_k); ValueFunc = u_final_c_osa_k + cS.beta*s_transition_final_osa_k*EV_A_final_osa_k; else, ValueFunc = u_final_c_osa_k; end; if ~isfinite(kPrime), kPrime = cS.kMin; end; if ~isfinite(c_quantity), c_quantity = cS.cFloor; end; if ~isfinite(ValueFunc), ValueFunc = -1e12; end
        end

        % --- HHSimulation_olgm (家庭模拟，现在需要处理k_pps作为状态) ---
        function [kHistM, kPpsHistM, cHistM] = HHSimulation_olgm(kPolM_4D, cPpsPolM_4D, cPolM_consump_q_4D, eIdxM, R_k_net_hh, w_gross_sim, TR_total_sim, bV_payg_sim, paramS_sim, cS)
            ageToGroupMap=zeros(cS.aD_orig,1); for a_map_idx=1:cS.aD_new, idx_map=cS.physAgeMap{a_map_idx}; if ~isempty(idx_map), ageToGroupMap(idx_map)=a_map_idx; end; end
            nSim=size(eIdxM,1); kHistM=zeros(nSim,cS.aD_orig+1); kPpsHistM=zeros(nSim,cS.aD_orig+1); cHistM=zeros(nSim,cS.aD_orig);
            leGridV_col=paramS_sim.leGridV(:);

            kPolInterp    = cell(cS.nw, cS.aD_new); cPpsPolInterp = cell(cS.nw, cS.aD_new); cPolqInterp   = cell(cS.nw, cS.aD_new);
            for ia = 1:cS.aD_new, for ie = 1:cS.nw
                    if cS.nk > 1 && cS.nkpps > 1
                        kPolInterp{ie,ia}    = griddedInterpolant({cS.kGridV, cS.kppsGridV}, squeeze(kPolM_4D(:,:,ie,ia)), 'linear', 'nearest');
                        cPpsPolInterp{ie,ia} = griddedInterpolant({cS.kGridV, cS.kppsGridV}, squeeze(cPpsPolM_4D(:,:,ie,ia)), 'linear', 'nearest');
                        cPolqInterp{ie,ia}   = griddedInterpolant({cS.kGridV, cS.kppsGridV}, squeeze(cPolM_consump_q_4D(:,:,ie,ia)), 'linear', 'nearest');
                    elseif cS.nk > 1 && cS.nkpps == 1 % k_pps is scalar
                        kPolInterp{ie,ia}    = griddedInterpolant(cS.kGridV, squeeze(kPolM_4D(:,1,ie,ia)), 'linear', 'nearest');
                        cPpsPolInterp{ie,ia} = griddedInterpolant(cS.kGridV, squeeze(cPpsPolM_4D(:,1,ie,ia)), 'linear', 'nearest');
                        cPolqInterp{ie,ia}   = griddedInterpolant(cS.kGridV, squeeze(cPolM_consump_q_4D(:,1,ie,ia)), 'linear', 'nearest');
                    elseif cS.nk == 1 && cS.nkpps > 1 % k is scalar
                        kPolInterp{ie,ia}    = griddedInterpolant(cS.kppsGridV, squeeze(kPolM_4D(1,:,ie,ia))', 'linear', 'nearest');
                        cPpsPolInterp{ie,ia} = griddedInterpolant(cS.kppsGridV, squeeze(cPpsPolM_4D(1,:,ie,ia))', 'linear', 'nearest');
                        cPolqInterp{ie,ia}   = griddedInterpolant(cS.kppsGridV, squeeze(cPolM_consump_q_4D(1,:,ie,ia))', 'linear', 'nearest');
                    elseif cS.nk == 1 && cS.nkpps == 1 % Both scalar
                        kPolInterp{ie,ia}    = @(x,y) kPolM_4D(1,1,ie,ia);
                        cPpsPolInterp{ie,ia} = @(x,y) cPpsPolM_4D(1,1,ie,ia);
                        cPolqInterp{ie,ia}   = @(x,y) cPolM_consump_q_4D(1,1,ie,ia);
                    else % nk=0 or nkpps=0 (empty grids) - should not happen if grids are well-defined
                        error('HHSimulation_olgm: nk or nkpps is zero, cannot create interpolants.');
                    end
            end; end

        pps_return_net_annual=(R_k_net_hh-1) + cS.pps_return_rate_premium;

        for a_orig_loop=1:cS.aD_orig
            a_new_group_idx = ageToGroupMap(a_orig_loop); kNowV_ann = kHistM(:,a_orig_loop); kPpsNowV_ann = kPpsHistM(:,a_orig_loop);
            kNextNonPpsV_ann_sim = zeros(nSim,1); cPpsDecisionV_ann_sim = zeros(nSim,1); cConsumpValV_q_ann_sim = zeros(nSim,1); kPpsNextV_ann_sim = zeros(nSim,1);
            model_age_current_idx_sim = a_orig_loop;
            is_working_age_sim = (model_age_current_idx_sim < cS.aR_idx_orig);
            is_pps_contrib_eligible_sim = (model_age_current_idx_sim <= cS.pps_contribution_age_max_idx && cS.pps_active);
            is_pps_withdrawal_eligible_sim = (~is_working_age_sim && cS.pps_active && model_age_current_idx_sim >= cS.pps_withdrawal_age_min_idx);
            pps_withdrawal_pretax_this_year = zeros(nSim,1); if is_pps_withdrawal_eligible_sim, pps_withdrawal_pretax_this_year = kPpsNowV_ann * cS.pps_withdrawal_rate; end

            for ie_sim_loop = 1 : cS.nw
                simIdxV = find(eIdxM(:, a_orig_loop) == ie_sim_loop); if isempty(simIdxV), continue; end
                kNow_cl = max(cS.kGridV(1), min(cS.kGridV(end), kNowV_ann(simIdxV)));
                kPpsNow_cl = max(cS.kppsGridV(1), min(cS.kppsGridV(end), kPpsNowV_ann(simIdxV)));

                if cS.nk > 1 && cS.nkpps > 1
                    kNextNonPpsV_ann_sim(simIdxV) = kPolInterp{ie_sim_loop, a_new_group_idx}(kNow_cl, kPpsNow_cl);
                    cPpsChosenByPol = cPpsPolInterp{ie_sim_loop, a_new_group_idx}(kNow_cl, kPpsNow_cl);
                    cConsumpValV_q_ann_sim(simIdxV) = cPolqInterp{ie_sim_loop, a_new_group_idx}(kNow_cl, kPpsNow_cl);
                elseif cS.nk > 1
                    kNextNonPpsV_ann_sim(simIdxV) = kPolInterp{ie_sim_loop, a_new_group_idx}(kNow_cl);
                    cPpsChosenByPol = cPpsPolInterp{ie_sim_loop, a_new_group_idx}(kNow_cl);
                    cConsumpValV_q_ann_sim(simIdxV) = cPolqInterp{ie_sim_loop, a_new_group_idx}(kNow_cl);
                elseif cS.nkpps > 1
                    kNextNonPpsV_ann_sim(simIdxV) = kPolInterp{ie_sim_loop, a_new_group_idx}(kPpsNow_cl); % Should be kPolInterp(kPpsNow_cl) if k is scalar
                    cPpsChosenByPol = cPpsPolInterp{ie_sim_loop, a_new_group_idx}(kPpsNow_cl);
                    cConsumpValV_q_ann_sim(simIdxV) = cPolqInterp{ie_sim_loop, a_new_group_idx}(kPpsNow_cl);
                else % Both scalar
                    kNextNonPpsV_ann_sim(simIdxV) = kPolInterp{ie_sim_loop, a_new_group_idx}(kNow_cl, kPpsNow_cl); % Call function handle
                    cPpsChosenByPol = cPpsPolInterp{ie_sim_loop, a_new_group_idx}(kNow_cl, kPpsNow_cl);
                    cConsumpValV_q_ann_sim(simIdxV) = cPolqInterp{ie_sim_loop, a_new_group_idx}(kNow_cl, kPpsNow_cl);
                end

                if is_pps_contrib_eligible_sim && cS.pps_max_contrib_frac > 0 % Check also if contrib frac > 0
                    current_gross_labor_income = w_gross_sim * cS.ageEffV_new(a_new_group_idx) * leGridV_col(ie_sim_loop);
                    max_cpps_by_frac = current_gross_labor_income * cS.pps_max_contrib_frac;
                    actual_cpps_contrib = min(cPpsChosenByPol, cS.pps_annual_contrib_limit);
                    actual_cpps_contrib = min(actual_cpps_contrib, max_cpps_by_frac);
                    cPpsDecisionV_ann_sim(simIdxV) = max(0, actual_cpps_contrib);
                else, cPpsDecisionV_ann_sim(simIdxV) = 0; end

                if cS.pps_active, kPpsNextV_ann_sim(simIdxV) = (kPpsNowV_ann(simIdxV) + cPpsDecisionV_ann_sim(simIdxV) - pps_withdrawal_pretax_this_year(simIdxV)) * (1 + pps_return_net_annual); kPpsNextV_ann_sim(simIdxV) = max(cS.kppsMin, min(cS.kppsMax, kPpsNextV_ann_sim(simIdxV))); else, kPpsNextV_ann_sim(simIdxV) = kPpsNowV_ann(simIdxV); end
            end
            kHistM(:,a_orig_loop+1)=max(cS.kMin,min(cS.kMax,kNextNonPpsV_ann_sim)); kPpsHistM(:,a_orig_loop+1)=max(cS.kppsMin,min(cS.kppsMax,kPpsNextV_ann_sim)); cHistM(:,a_orig_loop)=max(cS.cFloor,cConsumpValV_q_ann_sim);
        end
        kHistM = kHistM(:,1:cS.aD_orig); kPpsHistM = kPpsHistM(:,1:cS.aD_orig);
        end

        % --- 核心均衡求解器 (solve_K_tau_l_for_rho_prime) ---
        function [K_sol, tau_l_sol, gbc_res_final, converged_and_feasible, solution_details] = solve_K_tau_l_for_rho_prime(rho_prime_payg_target, K_init_guess, cS, paramS_in, eIdxM)
            K_guess = K_init_guess; tau_l_guess = cS.tau_l_init_guess; L_pc = paramS_in.L_per_capita; mass_workers = paramS_in.mass_workers_group;
            maxIter = cS.max_iter_K_tau_l; tol_norm = cS.tol_K_tau_l; dampK = cS.damp_K_v5; damp_tau_l = cS.damp_tau_l_v5;
            converged_and_feasible = false; K_sol = NaN; tau_l_sol = NaN; gbc_res_final = Inf; solution_details = struct();
            mass_retirees_calc = sum(paramS_in.ageMassV(cS.aR_new + 1 : cS.aD_new)); theta_payg_required = 0;
            if mass_workers > 1e-9, theta_payg_required = rho_prime_payg_target * (mass_retirees_calc / mass_workers); else, if rho_prime_payg_target > 1e-9, theta_payg_required = Inf; else, theta_payg_required = 0; end; end
            theta_payg_required = max(0, theta_payg_required); solution_details.theta_payg_required_before_cap = theta_payg_required;
            if theta_payg_required > cS.theta_payg_max + 1e-5, if ~isfield(paramS_in, 'suppress_initial_theta_print') || ~paramS_in.suppress_initial_theta_print, fprintf('  solve_K_tau_l_for_rho_prime: rho_prime_target=%.4f leads to theta_req=%.4f > theta_max=%.3f. Infeasible upfront.\n', rho_prime_payg_target, theta_payg_required, cS.theta_payg_max); end; converged_and_feasible = false; K_sol = K_init_guess; tau_l_sol = tau_l_guess; gbc_res_final = Inf; solution_details.theta_payg = min(theta_payg_required, cS.theta_payg_max); solution_details.MPL_gross = NaN; solution_details.R_mkt_gross = NaN; solution_details.b_payg = NaN; solution_details.T_bequest_Model = NaN; solution_details.C_model = NaN; solution_details.Y_model = NaN; solution_details.K_model_pps = NaN; solution_details.K_model_non_pps = NaN; return; end
            stagnation_counter = 0; prev_devNorm = Inf; tau_l_boundary_strike_count = 0;
            if ~isfield(paramS_in, 'suppress_inner_print_header') || ~paramS_in.suppress_inner_print_header, fprintf('  solve_K_tau_l_for_rho_prime: rho_prime_target=%.4f (theta_req=%.4f), K_init=%.2f, tau_l_init=%.3f\n', rho_prime_payg_target, theta_payg_required, K_guess, tau_l_guess); fprintf('  IterKTL | K_guess  | tau_l_gs | MPL_g    | theta_act| K_tot_mod| K_pps_mod| GBC_res  | K_dev    | tau_l_dev| Norm     | Improv   | Strikes\n'); fprintf('  -------------------------------------------------------------------------------------------------------------------------------------------\n'); end
            MPL_g = NaN; R_mkt_g = NaN; theta_payg_actual = NaN; b_payg = NaN; T_bequest_model_iter = NaN; C_model = NaN; Y_for_gbc = NaN; gbc_res = Inf; K_model = K_guess; K_dev = Inf; K_model_pps_iter = NaN; K_model_nonpps = NaN;
            for iter_ktl = 1:maxIter
                [R_mkt_g, MPL_g] = main_olg_v6_utils.HHPrices_Huggett(K_guess, L_pc, cS); r_mkt_g = R_mkt_g - 1;
                avg_worker_gross_wage = 0; if mass_workers > 1e-9 && L_pc > 0 && MPL_g > 0, avg_worker_gross_wage = (MPL_g * L_pc) / mass_workers; end; b_payg = rho_prime_payg_target * avg_worker_gross_wage; b_payg = max(0, b_payg);
                theta_payg_actual = theta_payg_required; if (theta_payg_actual + tau_l_guess) > cS.max_total_labor_tax, theta_payg_actual = max(0, cS.max_total_labor_tax - tau_l_guess); end; theta_payg_actual = max(0, theta_payg_actual);
                r_k_net_hh = r_mkt_g * (1 - cS.tau_k); R_k_net_hh = 1 + r_k_net_hh;
                bV_payg_vec = zeros(1, cS.aD_new); if cS.aR_new < cS.aD_new, bV_payg_vec(cS.aR_new+1:cS.aD_new) = b_payg; end
                paramS_iter = paramS_in; paramS_iter.tau_l = tau_l_guess; paramS_iter.theta_payg_actual_for_hh = theta_payg_actual;
                TR_total_for_vfi_guess = 0.05 * MPL_g; max_vfi_tr_iter = 5; tol_vfi_tr = 1e-3; T_bequest_model_iter = TR_total_for_vfi_guess;
                cPolM_4D_iter = zeros(cS.nk, cS.nkpps, cS.nw, cS.aD_new); % Initialize for this iteration
                kPolM_4D_iter = zeros(cS.nk, cS.nkpps, cS.nw, cS.aD_new);
                cPpsPolM_4D_iter = zeros(cS.nk, cS.nkpps, cS.nw, cS.aD_new);
                for i_vfi_tr = 1:max_vfi_tr_iter, [cPolM_vfi_temp, kPolM_vfi_temp, cPpsPolM_vfi_temp, ~] = main_olg_v6_utils.HHSolution_VFI_Huggett(R_k_net_hh, MPL_g, TR_total_for_vfi_guess, bV_payg_vec, paramS_iter, cS); [kHistM_sim_vfi_temp, kPpsHistM_sim_vfi_temp, ~] = main_olg_v6_utils.HHSimulation_olgm(kPolM_vfi_temp, cPpsPolM_vfi_temp, cPolM_vfi_temp, eIdxM, R_k_net_hh, MPL_g, TR_total_for_vfi_guess, bV_payg_vec, paramS_iter, cS); kprime_for_bequest_temp = zeros(cS.nSim, cS.aD_orig); if cS.aD_orig > 1, kprime_for_bequest_temp(:, 1:cS.aD_orig-1) = kHistM_sim_vfi_temp(:, 2:cS.aD_orig); end; ageDeathMass_annual_temp = paramS_in.Z_ss_norm_annual(:) .* cS.d_orig(:); mean_bequest_val_temp = mean(kprime_for_bequest_temp * R_k_net_hh, 1); TotalBequests_pc_temp = sum(mean_bequest_val_temp(:) .* ageDeathMass_annual_temp(:)); T_bequest_model_iter_new = TotalBequests_pc_temp / (1 + paramS_in.popGrowthForDebt); T_bequest_model_iter_new = max(0, T_bequest_model_iter_new); T_bequest_model_iter = T_bequest_model_iter_new; cPolM_4D_iter = cPolM_vfi_temp; kPolM_4D_iter = kPolM_vfi_temp; cPpsPolM_4D_iter = cPpsPolM_vfi_temp; if abs(T_bequest_model_iter_new - TR_total_for_vfi_guess) < tol_vfi_tr || i_vfi_tr == max_vfi_tr_iter, TR_total_for_vfi = T_bequest_model_iter_new; break; end; TR_total_for_vfi_guess = 0.5 * TR_total_for_vfi_guess + 0.5 * T_bequest_model_iter_new; end; T_bequest_model_iter = TR_total_for_vfi;
                [kHistM_non_pps_sim, kPpsHistM_sim, cHistM_sim] = main_olg_v6_utils.HHSimulation_olgm(kPolM_4D_iter, cPpsPolM_4D_iter, cPolM_4D_iter, eIdxM, R_k_net_hh, MPL_g, TR_total_for_vfi, bV_payg_vec, paramS_iter, cS);
                K_model_nonpps = mean(kHistM_non_pps_sim, 1) * paramS_in.Z_ss_norm_annual; K_model_pps_iter = 0; if cS.pps_active && cS.pps_in_K && cS.pps_max_contrib_frac > 0, K_model_pps_iter = mean(kPpsHistM_sim, 1) * paramS_in.Z_ss_norm_annual; K_model_pps_iter = max(0, K_model_pps_iter); end; K_model = K_model_nonpps + K_model_pps_iter; K_model = max(1e-6, K_model); C_model = mean(cHistM_sim,1) * paramS_in.Z_ss_norm_annual;
                Y_for_gbc = cS.A * (K_guess^cS.alpha) * (L_pc^(1-cS.alpha)); G_val = cS.gov_exp_frac_Y * Y_for_gbc; B_val = cS.gov_debt_frac_Y * Y_for_gbc;
                gbc_res = main_olg_v6_utils.check_gbc_residual(K_guess, C_model, Y_for_gbc, G_val, B_val, MPL_g, r_mkt_g, theta_payg_actual, tau_l_guess, b_payg, T_bequest_model_iter, 0, cS, paramS_in);
                K_dev = K_guess - K_model; tau_l_dev_raw = gbc_res / (MPL_g * L_pc + 1e-9); current_devNorm = sqrt(K_dev^2 + (gbc_res)^2 ); norm_improvement = prev_devNorm - current_devNorm;
                fprintf('  %7d | %8.4f | %8.4f | %8.4f | %8.4f | %8.4f | %8.4f | %8.2e | %8.4f | %8.4f | %8.2e | %.1e | %7d\n', iter_ktl, K_guess, tau_l_guess, MPL_g, theta_payg_actual, K_model, K_model_pps_iter, gbc_res, K_dev, tau_l_dev_raw, current_devNorm, norm_improvement, tau_l_boundary_strike_count);
                payg_fully_funded_by_actual_theta = (theta_payg_actual >= theta_payg_required - 1e-5);
                if current_devNorm < tol_norm && abs(gbc_res) < cS.gbc_tol_for_internal_loop && payg_fully_funded_by_actual_theta, converged_and_feasible = true; K_sol = K_model; tau_l_sol = tau_l_guess; gbc_res_final = gbc_res; solution_details.R_mkt_gross = R_mkt_g; solution_details.MPL_gross = MPL_g; solution_details.theta_payg = theta_payg_actual; solution_details.b_payg = b_payg; solution_details.T_bequest_Model = T_bequest_model_iter; solution_details.C_model = C_model; solution_details.Y_model = Y_for_gbc; solution_details.K_model_pps = K_model_pps_iter; solution_details.K_model_non_pps = K_model_nonpps; fprintf('  Convergence for K and tau_l achieved (rho_prime_target=%.4f, theta_act=%.4f).\n', rho_prime_payg_target, theta_payg_actual); break;
                elseif current_devNorm < tol_norm && abs(gbc_res) < cS.gbc_tol_for_internal_loop && ~payg_fully_funded_by_actual_theta, fprintf('  Converged K, tau_l, GBC for rho_prime=%.4f, BUT actual theta_payg (%.4f) (due to total tax cap) < required (%.4f). Marking as infeasible.\n', rho_prime_payg_target, theta_payg_actual, theta_payg_required); converged_and_feasible = false; K_sol = K_model; tau_l_sol = tau_l_guess; gbc_res_final = gbc_res; solution_details.R_mkt_gross = R_mkt_g; solution_details.MPL_gross = MPL_g; solution_details.theta_payg = theta_payg_actual; solution_details.b_payg = b_payg; solution_details.T_bequest_Model = T_bequest_model_iter; solution_details.C_model = C_model; solution_details.Y_model = Y_for_gbc; solution_details.K_model_pps = K_model_pps_iter; solution_details.K_model_non_pps = K_model_nonpps; break; end
                K_guess_next_iter = K_guess - dampK * K_dev; K_guess = max(1e-3, K_guess_next_iter);
                new_tau_l_unconstrained = tau_l_guess - damp_tau_l * tau_l_dev_raw; tau_l_next_iter_constrained = max(cS.tau_l_min, min(cS.tau_l_max, new_tau_l_unconstrained));
                is_tau_l_at_boundary = (abs(tau_l_next_iter_constrained - cS.tau_l_max) < 1e-7 && new_tau_l_unconstrained >= cS.tau_l_max) || (abs(tau_l_next_iter_constrained - cS.tau_l_min) < 1e-7 && new_tau_l_unconstrained <= cS.tau_l_min);
                if is_tau_l_at_boundary && abs(gbc_res) > cS.gbc_tol_for_internal_loop, tau_l_boundary_strike_count = tau_l_boundary_strike_count + 1; else, tau_l_boundary_strike_count = 0; end; tau_l_guess = tau_l_next_iter_constrained;
                if tau_l_boundary_strike_count >= cS.max_tau_l_boundary_strikes && abs(gbc_res) > cS.gbc_tol_for_internal_loop, fprintf('  WARNING: tau_l at boundary (%.4f) for %d strikes and GBC (%.2e) not balanced. Aborting for rho_prime=%.4f.\n', tau_l_guess, tau_l_boundary_strike_count, gbc_res, rho_prime_payg_target); converged_and_feasible = false; K_sol = K_model; tau_l_sol = tau_l_guess; gbc_res_final = gbc_res; solution_details.R_mkt_gross = R_mkt_g; solution_details.MPL_gross = MPL_g; solution_details.theta_payg = theta_payg_actual; solution_details.b_payg = b_payg; solution_details.T_bequest_Model = T_bequest_model_iter; solution_details.C_model = C_model; solution_details.Y_model = Y_for_gbc; solution_details.K_model_pps = K_model_pps_iter; solution_details.K_model_non_pps = K_model_nonpps; break; end
                if iter_ktl > 1, if norm_improvement < (cS.min_norm_improvement_frac * prev_devNorm) && current_devNorm > tol_norm, stagnation_counter = stagnation_counter + 1; else, stagnation_counter = 0; end; end; prev_devNorm = current_devNorm;
                if stagnation_counter >= cS.max_stagnation_iters && current_devNorm > tol_norm, fprintf('  WARNING: Norm stagnation detected after %d iterations. Norm: %.2e > Tol: %.1e. Aborting for rho_prime=%.4f.\n', iter_ktl, current_devNorm, tol_norm, rho_prime_payg_target); converged_and_feasible = false; K_sol = K_model; tau_l_sol = tau_l_guess; gbc_res_final = gbc_res; solution_details.R_mkt_gross = R_mkt_g; solution_details.MPL_gross = MPL_g; solution_details.theta_payg = theta_payg_actual; solution_details.b_payg = b_payg; solution_details.T_bequest_Model = T_bequest_model_iter; solution_details.C_model = C_model; solution_details.Y_model = Y_for_gbc; solution_details.K_model_pps = K_model_pps_iter; solution_details.K_model_non_pps = K_model_nonpps; break; end
            end
            if ~converged_and_feasible && iter_ktl == maxIter, fprintf('  WARNING: K and tau_l iteration Maxed out or was not feasible for rho_prime_target=%.4f.\n', rho_prime_payg_target); K_sol = K_model; tau_l_sol = tau_l_guess; gbc_res_final = gbc_res; if exist('MPL_g', 'var'), solution_details.R_mkt_gross = R_mkt_g; solution_details.MPL_gross = MPL_g; solution_details.theta_payg = theta_payg_actual; solution_details.b_payg = b_payg; solution_details.T_bequest_Model = T_bequest_model_iter; solution_details.C_model = C_model; solution_details.Y_model = Y_for_gbc; solution_details.K_model_pps = K_model_pps_iter; solution_details.K_model_non_pps = K_model_nonpps; else, solution_details.R_mkt_gross = NaN; solution_details.MPL_gross = NaN; solution_details.theta_payg = NaN; solution_details.b_payg = NaN; solution_details.T_bequest_Model = NaN; solution_details.C_model = NaN; solution_details.Y_model = NaN; solution_details.K_model_pps = NaN; solution_details.K_model_non_pps = NaN; end; end
            if ~isfield(solution_details, 'theta_payg_required_before_cap'), solution_details.theta_payg_required_before_cap = theta_payg_required; end
            if ~isfield(solution_details, 'K_model_pps') && exist('K_model_pps_iter','var'), solution_details.K_model_pps = K_model_pps_iter; elseif ~isfield(solution_details, 'K_model_pps'), solution_details.K_model_pps = NaN; end
            if ~isfield(solution_details, 'K_model_non_pps') && exist('K_model_nonpps','var'), solution_details.K_model_non_pps = K_model_nonpps; elseif ~isfield(solution_details, 'K_model_non_pps'), solution_details.K_model_non_pps = NaN; end
        end

        function gbc_residual = check_gbc_residual(K_val, C_val, Y_val, G_val, B_val, MPL_gross_val, r_mkt_gross_val, theta_payg_val_actual, tau_l_val, b_payg_val_target, T_bequest_val, TR_gov_val, cS, paramS_loc)
            L_pc_loc = paramS_loc.L_per_capita;
            LaborTaxRev_general_part = tau_l_val * MPL_gross_val * L_pc_loc;
            CapitalTaxRev = r_mkt_gross_val * K_val * cS.tau_k;
            ConsumptionTaxRev = C_val * cS.tau_c;
            GeneralRevenue = LaborTaxRev_general_part + CapitalTaxRev + ConsumptionTaxRev;
            GovConsumption = G_val;
            r_b_val = r_mkt_gross_val * (1 - cS.tau_k);
            DebtService = (r_b_val - paramS_loc.popGrowthForDebt) * B_val;
            GeneralOutlays = GovConsumption + DebtService + TR_gov_val;
            gbc_residual = GeneralRevenue - GeneralOutlays;
        end

        % --- NEW: Helper function to call interpolant safely ---
        function ev_val = CallInterpolator(interpolant_obj, k_val, k_pps_val, cS_local)
            ev_val = -Inf;
            try
                if isa(interpolant_obj, 'griddedInterpolant')
                    if cS_local.nk > 1 && cS_local.nkpps > 1
                        ev_val = interpolant_obj(k_val, k_pps_val);
                    elseif cS_local.nk > 1 % nkpps is 1 or 0
                        ev_val = interpolant_obj(k_val);
                    elseif cS_local.nkpps > 1 % nk is 1 or 0
                        ev_val = interpolant_obj(k_pps_val);
                    else % Both nk and nkpps are 1 or less (scalar value stored in function handle way)
                        if isscalar(interpolant_obj.Values)
                            ev_val = interpolant_obj.Values;
                        else % Fallback for a 1x1 grid but interpolant object still needs query
                            ev_val = interpolant_obj(k_val, k_pps_val); % This might error if not a function handle
                        end
                    end
                elseif isa(interpolant_obj, 'function_handle')
                    ev_val = interpolant_obj(k_val, k_pps_val); % Assumes fh takes 2 args
                else
                    warning('CallInterpolator: Unhandled interpolant type.'); ev_val = -1e11;
                end
            catch ME_call_interp
                warning('CallInterpolator: Error during interpolation: %s. Clamping and retrying.', ME_call_interp.message);
                k_cl = max(cS_local.kGridV(1), min(cS_local.kGridV(end), k_val));
                k_pps_cl = 0; if cS_local.nkpps > 0 && ~isempty(cS_local.kppsGridV), k_pps_cl = max(cS_local.kppsGridV(1), min(cS_local.kppsGridV(end), k_pps_val)); end
                try
                    if isa(interpolant_obj, 'griddedInterpolant')
                        if cS_local.nk > 1 && cS_local.nkpps > 1, ev_val = interpolant_obj(k_cl, k_pps_cl);
                        elseif cS_local.nk > 1, ev_val = interpolant_obj(k_cl);
                        elseif cS_local.nkpps > 1, ev_val = interpolant_obj(k_pps_cl);
                        else, if isa(interpolant_obj, 'function_handle'), ev_val = interpolant_obj(k_cl, k_pps_cl); else, if isscalar(interpolant_obj.Values), ev_val = interpolant_obj.Values; else, ev_val = -1e11; end; end; end
                    elseif isa(interpolant_obj, 'function_handle'), ev_val = interpolant_obj(k_cl, k_pps_cl);
                    end
                catch
                    ev_val = -1e11; % Final fallback
                end
            end
            if ~isfinite(ev_val), ev_val = -1e12; end
        end

    end % End methods (Static)
end % End classdef