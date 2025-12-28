% --- START OF FILE main_olg_v5_utils.m ---

classdef main_olg_v5_utils
    % OLG 模型 v5 的工具函数
    % 寻找最大可持续替代率 rho_prime_payg = b_payg / avg_worker_gross_wage
    % 内生调整 tau_l 使 TR_gov=0
    % PAYG 税率 theta_payg 内生但有上限 cS.theta_payg_max
    % 优化：1. 内层循环开始前检查theta_payg_required是否超限
    %       2. 内层循环中，若Norm停滞，提前终止
    %       3. 新增：内层循环中，若tau_l连续多次达到边界且GBC未平衡，提前终止

    methods (Static)

        function cS = ParameterValues_HuggettStyle()
            cS = main_olg_v4_utils_pps.ParameterValues_HuggettStyle();

            cS.tau_l_init_guess = 0.05;
            cS.tau_l_min = -0.05;
            cS.tau_l_max = 0.3;
            cS.max_total_labor_tax = 0.5;

            cS.theta_payg_max = 0.30;

            cS.max_iter_K_tau_l = 100;
            cS.tol_K_tau_l = 1e-4;
            cS.damp_K_v5 = 0.5;
            cS.damp_tau_l_v5 = 0.5;
            cS.gbc_tol_for_internal_loop = 1e-3;
            cS.gbc_tol_for_rho_search = 1e-2;

            cS.max_stagnation_iters = 10; % For Norm stagnation
            cS.min_norm_improvement_frac = 1e-3;
            
            % 新增：tau_l 连续碰壁检测参数
            cS.max_tau_l_boundary_strikes = 5; % tau_l连续碰壁的最大次数
        end

        function popS = initPopulation(cS); popS = main_olg_v4_utils_pps.initPopulation(cS); end
        function popS = populationDynamics(popS, cS); popS = main_olg_v4_utils_pps.populationDynamics(popS, cS); end
        function [Z_ss, dep_ratio_ss, bgp_reached, bgp_period] = detectSteadyStatePopulation(popS, cS); [Z_ss, dep_ratio_ss, bgp_reached, bgp_period] = main_olg_v4_utils_pps.detectSteadyStatePopulation(popS, cS); end
        function [logGridV, trProbM, prob1V] = EarningProcess_olgm(cS); [logGridV, trProbM, prob1V] = main_olg_v4_utils_pps.EarningProcess_olgm(cS); end
        function eIdxM = LaborEndowSimulation_olgm(cS, paramS); eIdxM = main_olg_v4_utils_pps.LaborEndowSimulation_olgm(cS, paramS); end
        function [HHlaborM_group, L] = LaborSupply_Huggett(eIdxM, cS, paramS, Z_ss_norm_group); [HHlaborM_group, L] = main_olg_v4_utils_pps.LaborSupply_Huggett(eIdxM, cS, paramS, Z_ss_norm_group); end
        function [R_market_gross_factor, MPL_gross] = HHPrices_Huggett(K_productive, L_total_eff, cS); [R_market_gross_factor, MPL_gross] = main_olg_v4_utils_pps.HHPrices_Huggett(K_productive, L_total_eff, cS); end
        function [cPolM_quantity, kPolM, cPpsPolM, valueM] = HHSolution_VFI_Huggett(R_k_net_factor_vfi, w_net_hh_vfi, TR_total_vfi, bV_payg_vfi, paramS_vfi, cS); [cPolM_quantity, kPolM, cPpsPolM, valueM] = main_olg_v4_utils_pps.HHSolution_VFI_Huggett(R_k_net_factor_vfi, w_net_hh_vfi, TR_total_vfi, bV_payg_vfi, paramS_vfi, cS); end
        function [kHistM, kPpsHistM, cHistM] = HHSimulation_olgm(kPolM, cPpsPolM, cPolM_consump, eIdxM, R_k_net_factor_hh, w_net_hh, TR_total, bV_payg, paramS, cS); [kHistM, kPpsHistM, cHistM] = main_olg_v4_utils_pps.HHSimulation_olgm(kPolM, cPpsPolM, cPolM_consump, eIdxM, R_k_net_factor_hh, w_net_hh, TR_total, bV_payg, paramS, cS);end
        function [muM, utilM] = CES_utility(cM_quantity, sig, cS); [muM, utilM] = main_olg_v4_utils_pps.CES_utility(cM_quantity, sig, cS); end
        function incomeM_resources = HHIncome_Huggett(k_non_pps_grid, R_k_net_factor_hh_price, w_net_hh_price, TR_total_price, b_payg_price_group, a_new_idx_income, paramS_income, cS); incomeM_resources = main_olg_v4_utils_pps.HHIncome_Huggett(k_non_pps_grid, R_k_net_factor_hh_price, w_net_hh_price, TR_total_price, b_payg_price_group, a_new_idx_income, paramS_income, cS); end

        function [K_sol, tau_l_sol, gbc_res_final, converged_and_feasible, solution_details] = solve_K_tau_l_for_rho_prime(rho_prime_payg_target, K_init_guess, cS, paramS_in, eIdxM)
            K_guess = K_init_guess;
            tau_l_guess = cS.tau_l_init_guess;
            L_pc = paramS_in.L_per_capita;
            mass_workers = paramS_in.mass_workers_group;

            maxIter = cS.max_iter_K_tau_l;
            tol_norm = cS.tol_K_tau_l;
            dampK = cS.damp_K_v5;
            damp_tau_l = cS.damp_tau_l_v5;

            converged_and_feasible = false;
            K_sol = NaN; tau_l_sol = NaN; gbc_res_final = Inf;
            solution_details = struct();
            
            mass_retirees_calc = sum(paramS_in.ageMassV(cS.aR_new + 1 : cS.aD_new));
            theta_payg_required = 0;
            if mass_workers > 1e-9
                theta_payg_required = rho_prime_payg_target * (mass_retirees_calc / mass_workers);
            else
                if rho_prime_payg_target > 1e-9, theta_payg_required = Inf; else, theta_payg_required = 0; end
            end
            theta_payg_required = max(0, theta_payg_required);
            solution_details.theta_payg_required_before_cap = theta_payg_required;

            if theta_payg_required > cS.theta_payg_max + 1e-5 
                if ~isfield(paramS_in, 'suppress_initial_theta_print') || ~paramS_in.suppress_initial_theta_print
                    fprintf('  solve_K_tau_l_for_rho_prime: rho_prime_target=%.4f leads to theta_req=%.4f > theta_max=%.3f. Infeasible upfront.\n', ...
                             rho_prime_payg_target, theta_payg_required, cS.theta_payg_max);
                end
                converged_and_feasible = false; K_sol = K_init_guess; tau_l_sol = tau_l_guess; gbc_res_final = Inf; 
                solution_details.theta_payg = min(theta_payg_required, cS.theta_payg_max);
                solution_details.MPL_gross = NaN; solution_details.R_mkt_gross = NaN; solution_details.b_payg = NaN; 
                solution_details.T_bequest_Model = NaN; solution_details.C_model = NaN; solution_details.Y_model = NaN;
                return; 
            end

            stagnation_counter = 0;
            prev_devNorm = Inf;
            tau_l_boundary_strike_count = 0; % NEW: Counter for tau_l hitting boundary

            if ~isfield(paramS_in, 'suppress_inner_print_header') || ~paramS_in.suppress_inner_print_header
                fprintf('  solve_K_tau_l_for_rho_prime: rho_prime_target=%.4f (theta_req=%.4f), K_init=%.2f, tau_l_init=%.3f\n', rho_prime_payg_target, theta_payg_required, K_guess, tau_l_guess);
                fprintf('  IterKTL | K_guess  | tau_l_gs | MPL_g    | theta_act| K_model  | GBC_res  | K_dev    | tau_l_dev| Norm     | Improv   | Strikes\n');
                fprintf('  -----------------------------------------------------------------------------------------------------------------------------\n');
            end
            
            MPL_g = NaN; R_mkt_g = NaN; theta_payg_actual = NaN; b_payg = NaN;
            T_bequest_model_iter = NaN; C_model = NaN; Y_for_gbc = NaN; gbc_res = Inf;
            K_model = K_guess; K_dev = Inf;

            for iter_ktl = 1:maxIter
                [R_mkt_g, MPL_g] = main_olg_v5_utils.HHPrices_Huggett(K_guess, L_pc, cS);
                r_mkt_g = R_mkt_g - 1;

                avg_worker_gross_wage = 0;
                if mass_workers > 1e-9 && L_pc > 0 && MPL_g > 0, avg_worker_gross_wage = (MPL_g * L_pc) / mass_workers; end
                b_payg = rho_prime_payg_target * avg_worker_gross_wage; b_payg = max(0, b_payg);

                theta_payg_actual = theta_payg_required; 
                if (theta_payg_actual + tau_l_guess) > cS.max_total_labor_tax
                    theta_payg_actual = max(0, cS.max_total_labor_tax - tau_l_guess);
                end
                theta_payg_actual = max(0, theta_payg_actual);

                r_k_net_hh = r_mkt_g * (1 - cS.tau_k); R_k_net_hh = 1 + r_k_net_hh;
                w_net_hh = MPL_g * (1 - theta_payg_actual - tau_l_guess); w_net_hh = max(0, w_net_hh);
                bV_payg = zeros(1, cS.aD_new); if cS.aR_new < cS.aD_new, bV_payg(cS.aR_new+1:cS.aD_new) = b_payg; end

                paramS_iter = paramS_in; paramS_iter.tau_l = tau_l_guess;
                TR_total_for_vfi_guess = 0.05 * MPL_g; max_vfi_tr_iter = 5; tol_vfi_tr = 1e-3;
                T_bequest_model_iter = TR_total_for_vfi_guess; cPolM = []; kPolM = []; cPpsPolM = [];

                for i_vfi_tr = 1:max_vfi_tr_iter
                    [cPolM_vfi_temp, kPolM_vfi_temp, cPpsPolM_vfi_temp, ~] = main_olg_v5_utils.HHSolution_VFI_Huggett(R_k_net_hh, w_net_hh, TR_total_for_vfi_guess, bV_payg, paramS_iter, cS);
                    [kHistM_sim_vfi_temp, ~, ~] = main_olg_v5_utils.HHSimulation_olgm(kPolM_vfi_temp, cPpsPolM_vfi_temp, cPolM_vfi_temp, eIdxM, R_k_net_hh, w_net_hh, TR_total_for_vfi_guess, bV_payg, paramS_iter, cS);
                    kprime_for_bequest_temp = zeros(cS.nSim, cS.aD_orig); if cS.aD_orig > 1, kprime_for_bequest_temp(:, 1:cS.aD_orig-1) = kHistM_sim_vfi_temp(:, 2:cS.aD_orig); end
                    ageDeathMass_annual_temp = paramS_in.Z_ss_norm_annual(:) .* cS.d_orig(:); mean_bequest_val_temp = mean(kprime_for_bequest_temp * R_k_net_hh, 1);
                    TotalBequests_pc_temp = sum(mean_bequest_val_temp(:) .* ageDeathMass_annual_temp(:));
                    T_bequest_model_iter_new = TotalBequests_pc_temp / (1 + paramS_in.popGrowthForDebt); T_bequest_model_iter_new = max(0, T_bequest_model_iter_new);
                    T_bequest_model_iter = T_bequest_model_iter_new; cPolM = cPolM_vfi_temp; kPolM = kPolM_vfi_temp; cPpsPolM = cPpsPolM_vfi_temp;
                    if abs(T_bequest_model_iter_new - TR_total_for_vfi_guess) < tol_vfi_tr || i_vfi_tr == max_vfi_tr_iter, TR_total_for_vfi = T_bequest_model_iter_new; break; end
                    TR_total_for_vfi_guess = 0.5 * TR_total_for_vfi_guess + 0.5 * T_bequest_model_iter_new;
                end; T_bequest_model_iter = TR_total_for_vfi;

                [kHistM_non_pps, kPpsHistM, cHistM] = main_olg_v5_utils.HHSimulation_olgm(kPolM, cPpsPolM, cPolM, eIdxM, R_k_net_hh, w_net_hh, TR_total_for_vfi, bV_payg, paramS_iter, cS);
                K_model_nonpps = mean(kHistM_non_pps, 1) * paramS_in.Z_ss_norm_annual;
                K_model_pps = 0; if cS.pps_active && cS.pps_in_K && cS.pps_max_contrib_frac > 1e-9, K_model_pps = mean(kPpsHistM, 1) * paramS_in.Z_ss_norm_annual; end
                K_model = K_model_nonpps + K_model_pps; K_model = max(1e-6, K_model);
                C_model = mean(cHistM,1) * paramS_in.Z_ss_norm_annual;

                Y_for_gbc = cS.A * (K_guess^cS.alpha) * (L_pc^(1-cS.alpha));
                G_val = cS.gov_exp_frac_Y * Y_for_gbc; B_val = cS.gov_debt_frac_Y * Y_for_gbc;
                gbc_res = main_olg_v5_utils.check_gbc_residual(K_guess, C_model, Y_for_gbc, G_val, B_val, MPL_g, r_mkt_g, theta_payg_actual, tau_l_guess, b_payg, T_bequest_model_iter, 0, cS, paramS_in);

                K_dev = K_guess - K_model;
                tau_l_dev_raw = gbc_res / (MPL_g * L_pc + 1e-9);
                current_devNorm = sqrt(K_dev^2 + (gbc_res)^2 ); 
                norm_improvement = prev_devNorm - current_devNorm;

                fprintf('  %7d | %8.4f | %8.4f | %8.4f | %8.4f | %8.4f | %8.2e | %8.4f | %8.4f | %8.2e | %.1e | %7d\n', ...
                         iter_ktl, K_guess, tau_l_guess, MPL_g, theta_payg_actual, K_model, gbc_res, K_dev, tau_l_dev_raw, current_devNorm, norm_improvement, tau_l_boundary_strike_count);
                
                payg_fully_funded_by_actual_theta = (theta_payg_actual >= theta_payg_required - 1e-5);

                if current_devNorm < tol_norm && abs(gbc_res) < cS.gbc_tol_for_internal_loop && payg_fully_funded_by_actual_theta
                    converged_and_feasible = true; K_sol = K_model; tau_l_sol = tau_l_guess; gbc_res_final = gbc_res;
                    solution_details.R_mkt_gross = R_mkt_g; solution_details.MPL_gross = MPL_g; solution_details.theta_payg = theta_payg_actual;
                    solution_details.b_payg = b_payg; solution_details.T_bequest_Model = T_bequest_model_iter;
                    solution_details.C_model = C_model; solution_details.Y_model = Y_for_gbc;
                    fprintf('  Convergence for K and tau_l achieved (rho_prime_target=%.4f, theta_act=%.4f).\n', rho_prime_payg_target, theta_payg_actual);
                    break;
                elseif current_devNorm < tol_norm && abs(gbc_res) < cS.gbc_tol_for_internal_loop && ~payg_fully_funded_by_actual_theta
                    fprintf('  Converged K, tau_l, GBC for rho_prime=%.4f, BUT actual theta_payg (%.4f) (due to total tax cap) < required (%.4f). Marking as infeasible.\n', rho_prime_payg_target, theta_payg_actual, theta_payg_required);
                    converged_and_feasible = false; K_sol = K_model; tau_l_sol = tau_l_guess; gbc_res_final = gbc_res;
                    solution_details.R_mkt_gross = R_mkt_g; solution_details.MPL_gross = MPL_g; solution_details.theta_payg = theta_payg_actual;
                    solution_details.b_payg = b_payg; solution_details.T_bequest_Model = T_bequest_model_iter;
                    solution_details.C_model = C_model; solution_details.Y_model = Y_for_gbc;
                    break; 
                end
                
                K_guess_next_iter = K_guess - dampK * K_dev; K_guess = max(1e-3, K_guess_next_iter);
                
                new_tau_l_unconstrained = tau_l_guess - damp_tau_l * tau_l_dev_raw; % if gbc_res > 0 (surplus), dev_raw > 0, tau_l decreases. Correct.
                
                tau_l_next_iter_constrained = max(cS.tau_l_min, min(cS.tau_l_max, new_tau_l_unconstrained));
                
                % Check if tau_l is hitting a boundary AND GBC is still significantly off
                is_tau_l_at_boundary = (abs(tau_l_next_iter_constrained - cS.tau_l_max) < 1e-7 && new_tau_l_unconstrained >= cS.tau_l_max) || ...
                                       (abs(tau_l_next_iter_constrained - cS.tau_l_min) < 1e-7 && new_tau_l_unconstrained <= cS.tau_l_min);
                
                if is_tau_l_at_boundary && abs(gbc_res) > cS.gbc_tol_for_internal_loop
                    tau_l_boundary_strike_count = tau_l_boundary_strike_count + 1;
                else
                    tau_l_boundary_strike_count = 0; % Reset if not at boundary or GBC is fine
                end
                
                tau_l_guess = tau_l_next_iter_constrained; % Apply boundary constrained tau_l

                % Total labor tax cap (may further adjust tau_l_guess if theta_payg_actual is high)
                % theta_payg_actual was already determined based on previous tau_l_guess and total cap.
                % If new tau_l_guess makes total tax too high, it means original theta_payg_actual
                % (which is theta_payg_required here) was too high for *this new* tau_l.
                % This implies an inconsistency if not handled carefully.
                % A simpler approach: the theta_payg_actual calc at start of loop already factors in current tau_l_guess.
                % The tau_l update then tries to fix GBC. If this tau_l hits its own boundary, that's the primary issue.

                if tau_l_boundary_strike_count >= cS.max_tau_l_boundary_strikes && abs(gbc_res) > cS.gbc_tol_for_internal_loop
                     fprintf('  WARNING: tau_l at boundary (%.4f) for %d strikes and GBC (%.2e) not balanced. Aborting for rho_prime=%.4f.\n', ...
                              tau_l_guess, tau_l_boundary_strike_count, gbc_res, rho_prime_payg_target);
                     converged_and_feasible = false; K_sol = K_model; tau_l_sol = tau_l_guess; gbc_res_final = gbc_res;
                     solution_details.R_mkt_gross = R_mkt_g; solution_details.MPL_gross = MPL_g; solution_details.theta_payg = theta_payg_actual;
                     solution_details.b_payg = b_payg; solution_details.T_bequest_Model = T_bequest_model_iter;
                     solution_details.C_model = C_model; solution_details.Y_model = Y_for_gbc;
                     break;
                end
                
                if iter_ktl > 1
                    if norm_improvement < (cS.min_norm_improvement_frac * prev_devNorm) && current_devNorm > tol_norm
                        stagnation_counter = stagnation_counter + 1;
                    else, stagnation_counter = 0; end
                end; prev_devNorm = current_devNorm;

                if stagnation_counter >= cS.max_stagnation_iters && current_devNorm > tol_norm
                    fprintf('  WARNING: Norm stagnation detected after %d iterations. Norm: %.2e > Tol: %.1e. Aborting for rho_prime=%.4f.\n', iter_ktl, current_devNorm, tol_norm, rho_prime_payg_target);
                    converged_and_feasible = false; K_sol = K_model; tau_l_sol = tau_l_guess; gbc_res_final = gbc_res;
                    solution_details.R_mkt_gross = R_mkt_g; solution_details.MPL_gross = MPL_g; solution_details.theta_payg = theta_payg_actual;
                    solution_details.b_payg = b_payg; solution_details.T_bequest_Model = T_bequest_model_iter;
                    solution_details.C_model = C_model; solution_details.Y_model = Y_for_gbc;
                    break;
                end
            end

            if ~converged_and_feasible && iter_ktl == maxIter
                fprintf('  WARNING: K and tau_l iteration Maxed out or was not feasible for rho_prime_target=%.4f.\n', rho_prime_payg_target);
                K_sol = K_model; tau_l_sol = tau_l_guess; gbc_res_final = gbc_res;
                if exist('MPL_g', 'var')
                    solution_details.R_mkt_gross = R_mkt_g; solution_details.MPL_gross = MPL_g; solution_details.theta_payg = theta_payg_actual;
                    solution_details.b_payg = b_payg; solution_details.T_bequest_Model = T_bequest_model_iter;
                    solution_details.C_model = C_model; solution_details.Y_model = Y_for_gbc;
                else
                    solution_details.R_mkt_gross = NaN; solution_details.MPL_gross = NaN; solution_details.theta_payg = NaN; 
                    solution_details.b_payg = NaN; solution_details.T_bequest_Model = NaN; solution_details.C_model = NaN; solution_details.Y_model = NaN;
                end
            end
            if ~isfield(solution_details, 'theta_payg_required_before_cap')
                solution_details.theta_payg_required_before_cap = theta_payg_required;
            end
        end

        function gbc_residual = check_gbc_residual(K_val, C_val, Y_val, G_val, B_val, ...
                                            MPL_gross_val, r_mkt_gross_val, ...
                                            theta_payg_val_actual, tau_l_val, b_payg_val_target, ...
                                            T_bequest_val, TR_gov_val, ...
                                            cS, paramS_loc)
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
    end % End methods (Static)
end % End classdef