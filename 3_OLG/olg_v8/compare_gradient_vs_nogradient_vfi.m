% --- test_gradient_calculation_enhanced.m ---
% 脚本：全面检测目标函数的解析梯度与数值梯度的差异 (增强版)

clc;
clear;
close all;
fprintf('=== 目标函数梯度检测脚本 (增强版) ===\n');

%% 1. 加载模型参数和必要的辅助数据 (与之前类似)
fprintf('\n--- 1. 加载参数 ---\n');
cS = main_olg_v8_utils.ParameterValues_HuggettStyle();
paramS_test = struct();
[paramS_test.leLogGridV, paramS_test.leTrProbM, paramS_test.leProb1V] = main_olg_v8_utils.EarningProcess_olgm(cS);
paramS_test.leGridV = exp(paramS_test.leLogGridV(:));
paramS_test.tau_l = 0.15; 
paramS_test.theta_payg_actual_for_hh = 0.10;
paramS_test.pps_tax_deferral_active = cS.pps_active;

fprintf('正在创建简化的下一期期望值函数插值器...\n');
vPrime_dummy = rand(cS.nk, cS.nkpps, cS.nw) * 50 - 100;
EV_matrix_dummy = zeros(cS.nk, cS.nkpps, cS.nw);
for ie_current_test = 1:cS.nw
    transition_probs_test = paramS_test.leTrProbM(ie_current_test, :);
    EV_slice_test = sum(vPrime_dummy .* reshape(transition_probs_test, 1, 1, cS.nw), 3);
    EV_matrix_dummy(:, :, ie_current_test) = EV_slice_test;
end
EV_interpolants_test = cell(cS.nw, 1);
for ie_current_test = 1:cS.nw
    if cS.nk > 1 && cS.nkpps > 1
        EV_interpolants_test{ie_current_test} = griddedInterpolant({cS.kGridV, cS.kppsGridV}, EV_matrix_dummy(:, :, ie_current_test), 'spline', 'spline');
    elseif cS.nk > 1 && cS.nkpps == 1, EV_interpolants_test{ie_current_test} = griddedInterpolant(cS.kGridV, EV_matrix_dummy(:, 1, ie_current_test), 'spline', 'spline');
    elseif cS.nk == 1 && cS.nkpps > 1, EV_interpolants_test{ie_current_test} = griddedInterpolant(cS.kppsGridV, squeeze(EV_matrix_dummy(1, :, ie_current_test)), 'spline', 'spline');
    else, EV_interpolants_test{ie_current_test} = @(k_s, kp_s) EV_matrix_dummy(1, 1, ie_current_test);
    end
end
fprintf('插值器创建完毕。\n');

R_k_net_factor_age_test = 1.02;
w_gross_age_test = 1.5;
TR_total_age_test = 0.05;
b_age_val_test = 0.1;

%% 2. 定义测试状态点和决策变量 (与之前类似)
fprintf('\n--- 2. 定义测试点和参数 ---\n');
test_ages_idx = [1, round(cS.aR_new/2), cS.aR_new-1];
test_ages_idx = unique(max(1, min(test_ages_idx, cS.aD_new-1)));
test_k_indices = round(linspace(1, cS.nk, min(3, cS.nk)));
test_kpps_indices = round(linspace(1, cS.nkpps, min(3, cS.nkpps)));
if cS.nkpps == 0, test_kpps_indices = 1; end
test_ie_indices = round(linspace(1, cS.nw, min(3, cS.nw)));
test_x_prop_values = [
    0.0, 0.0; 0.5, 0.5; 1.0, 1.0; 0.1, 0.8; 0.8, 0.1;
    0.0, 0.5; 0.5, 0.0; 0.0, 1.0; 1.0, 0.0; % More boundary cases for x_prop
    rand(1,2); rand(1,2)
];
test_x_prop_values(end-1:end,:) = max(0, min(1, test_x_prop_values(end-1:end,:)));

h_fd = 1e-7; 
gradient_diff_threshold = 1e-4;
func_val_diff_threshold = 1e-6; % Stricter for function value

total_tests = 0;
failed_tests_fval = 0;
failed_tests_grad = 0;

%% 3. 辅助函数：用于提取中间变量和梯度的修改版目标函数
% 我们需要在 fmincon_objective_helper_proportional_with_deri 内部能够访问中间变量
% 为了测试，最简单的方法是暂时修改该函数，使其返回一个包含所有中间计算结果的结构体
% 或者，在测试脚本中重新实现这些计算步骤，但这容易引入不一致。
% 暂时假设我们可以通过某种方式（例如，全局变量，或修改函数签名使其返回更多）
% 获取这些中间导数。为了这个测试脚本，我们将“假设”这些中间解析导数是
% 从 _with_deri 函数中正确计算并提取出来的。

% 对于数值梯度，我们将定义小函数来获取中间步骤的值
value_wrapper_actual_c_pps = @(pps_prop_in, max_cpps) pps_prop_in * max_cpps; % Simplified

value_wrapper_resources_after_pps = @(acpps_in, k_s, k_pps_s, eps_s, aidx, ie_s, Rk, wg, TR, b, pS, cSf, max_cpps_ignored_here) ...
    main_olg_v8_utils.HHIncome_Huggett(k_s, Rk, wg, TR, b, acpps_in, aidx, pS, cSf, eps_s);

value_wrapper_actual_k_prime = @(kp_prop_in, res_after_pps_in, cSf_local) ...
    calculate_actual_k_prime_for_test(kp_prop_in, res_after_pps_in, cSf_local);

value_wrapper_current_c = @(akp_in, res_after_pps_in, cSf_local) ...
    calculate_current_c_for_test(akp_in, res_after_pps_in, cSf_local);

value_wrapper_k_pps_prime = @(acpps_in, k_pps_s_in, aidx_in, Rk_in, cSf_local) ...
    calculate_k_pps_prime_for_test(acpps_in, k_pps_s_in, aidx_in, Rk_in, cSf_local);


%% 4. 循环测试
fprintf('\n--- 3. 开始梯度检测循环 ---\n');
% Header
fprintf('%-4s|%-4s|%-4s|%-4s|%-12s|%-9s|%-9s|%-9s|%-9s|%-9s|%-9s|%-7s| Status\n', ...
    'Age','kId','kpId','eId', 'x_prop', 'f_noG', 'f_wG', 'dFdP_A','dFdP_N','dFdK_A','dFdK_N','RelDiff');
disp(repmat('-', 1, 130));

for ia_idx_loop = 1:length(test_ages_idx)
    a_idx = test_ages_idx(ia_idx_loop);

    for ik_idx_loop = 1:length(test_k_indices)
        ik = test_k_indices(ik_idx_loop);
        k_state = cS.kGridV(ik);

        for ikpps_idx_loop = 1:length(test_kpps_indices)
            ikpps = test_kpps_indices(ikpps_idx_loop);
            if cS.nkpps > 0, k_pps_state = cS.kppsGridV(ikpps); else, k_pps_state = 0; end

            for ie_idx_loop = 1:length(test_ie_indices)
                ie = test_ie_indices(ie_idx_loop);
                epsilon_state = paramS_test.leGridV(ie);

                max_permissible_cpps_current_state = 0;
                model_age_group_start_year_idx_test = cS.physAgeMap{a_idx}(1);
                is_pps_eligible_test = (a_idx <= cS.aR_new && model_age_group_start_year_idx_test <= cS.pps_contribution_age_max_idx && cS.pps_active);
                if is_pps_eligible_test
                    age_efficiency_test = cS.ageEffV_new(a_idx);
                    current_gross_labor_income_test = w_gross_age_test * age_efficiency_test * epsilon_state;
                    if current_gross_labor_income_test > 1e-6
                        max_cpps_by_frac_test = current_gross_labor_income_test * cS.pps_max_contrib_frac;
                        max_permissible_cpps_current_state = min(cS.pps_annual_contrib_limit, max_cpps_by_frac_test);
                        max_permissible_cpps_current_state = max(0, max_permissible_cpps_current_state);
                    end
                end
                
                for i_xprop = 1:size(test_x_prop_values, 1)
                    x_prop_test_orig = test_x_prop_values(i_xprop, :);
                    x_prop_test = x_prop_test_orig; % Make a copy for modification
                    
                    if max_permissible_cpps_current_state < 1e-9
                        x_prop_test(1) = 0; % Force pps_proportion to 0 if no contribution possible
                    end
                    x_prop_test = max(0, min(1, x_prop_test)); % Ensure in [0,1]

                    total_tests = total_tests + 1;
                    current_status = "PASS";
                    error_details = "";

                    % === Test Function Value Consistency ===
                    f_val_no_grad = main_olg_v8_utils.fmincon_objective_helper_proportional(...
                        x_prop_test, k_state, k_pps_state, epsilon_state, a_idx, ie, ...
                        R_k_net_factor_age_test, w_gross_age_test, TR_total_age_test, b_age_val_test, ...
                        paramS_test, cS, EV_interpolants_test, max_permissible_cpps_current_state);

                    [f_val_with_grad, grad_analytic] = main_olg_v8_utils.fmincon_objective_helper_proportional_with_deri(...
                        x_prop_test, k_state, k_pps_state, epsilon_state, a_idx, ie, ...
                        R_k_net_factor_age_test, w_gross_age_test, TR_total_age_test, b_age_val_test, ...
                        paramS_test, cS, EV_interpolants_test, max_permissible_cpps_current_state);

                    func_val_diff = abs(f_val_no_grad - f_val_with_grad);
                    if func_val_diff > func_val_diff_threshold
                        current_status = "FAIL_FVAL";
                        failed_tests_fval = failed_tests_fval + 1;
                        error_details = sprintf('fVal diff: %.2e; ', func_val_diff);
                    end

                    % === Test Gradient Numerically (only if function values are consistent) ===
                    grad_numeric = zeros(size(x_prop_test));
                    if ~strcmp(current_status, "FAIL_FVAL")
                        for i_dim = 1:length(x_prop_test)
                            x_plus_h = x_prop_test;
                            x_minus_h = x_prop_test;
                            
                            % Central difference, careful at bounds [0,1]
                            if x_prop_test(i_dim) + h_fd <= 1
                                x_plus_h(i_dim) = x_prop_test(i_dim) + h_fd;
                            else % At upper bound, use backward difference effectively for this point
                                x_plus_h(i_dim) = x_prop_test(i_dim); % Will make f_plus = f_center
                            end
                            
                            if x_prop_test(i_dim) - h_fd >= 0
                                x_minus_h(i_dim) = x_prop_test(i_dim) - h_fd;
                            else % At lower bound
                                x_minus_h(i_dim) = x_prop_test(i_dim); % Will make f_minus = f_center
                            end

                            % Use the non-gradient version for f_plus and f_minus
                            f_plus = main_olg_v8_utils.fmincon_objective_helper_proportional(...
                                x_plus_h, k_state, k_pps_state, epsilon_state, a_idx, ie, ...
                                R_k_net_factor_age_test, w_gross_age_test, TR_total_age_test, b_age_val_test, ...
                                paramS_test, cS, EV_interpolants_test, max_permissible_cpps_current_state);
                            
                            f_minus = main_olg_v8_utils.fmincon_objective_helper_proportional(...
                                x_minus_h, k_state, k_pps_state, epsilon_state, a_idx, ie, ...
                                R_k_net_factor_age_test, w_gross_age_test, TR_total_age_test, b_age_val_test, ...
                                paramS_test, cS, EV_interpolants_test, max_permissible_cpps_current_state);
                            
                            denominator = (x_plus_h(i_dim) - x_minus_h(i_dim));
                            if abs(denominator) < 1e-12 % Avoid division by zero if no change in x
                                grad_numeric(i_dim) = 0;
                            else
                                grad_numeric(i_dim) = (f_plus - f_minus) / denominator;
                            end
                        end

                        grad_abs_diff = norm(grad_analytic - grad_numeric);
                        grad_rel_diff = grad_abs_diff / (0.5 * (norm(grad_analytic) + norm(grad_numeric)) + 1e-8);

                        if grad_rel_diff > gradient_diff_threshold
                            current_status = "FAIL_GRAD";
                            failed_tests_grad = failed_tests_grad + 1;
                            error_details = [error_details, sprintf('Grad RelDiff: %.2e; An=[%.1e,%.1e],Num=[%.1e,%.1e]', ...
                                grad_rel_diff, grad_analytic(1), grad_analytic(2), grad_numeric(1), grad_numeric(2))];
                        end
                    else
                         grad_rel_diff = NaN; % Cannot compute grad diff if fvals are different
                    end
                    
                    fprintf('%-4d|%-4d|%-4d|%-4d| [%.2f, %.2f]   |%-9.2e|%-9.2e|%-9.2e|%-9.2e|%-9.2e|%-9.2e|%-7.2e| %s\n', ...
                        a_idx, ik, ikpps, ie, x_prop_test(1), x_prop_test(2), ...
                        f_val_no_grad, f_val_with_grad, ...
                        grad_analytic(1), grad_numeric(1), grad_analytic(2), grad_numeric(2), ...
                        grad_rel_diff, current_status);
                    if ~strcmp(current_status, "PASS")
                        disp(['  Details: ', error_details]);
                    end
                end
            end
        end
    end
end

disp(repmat('-', 1, 130));
fprintf('\n--- 4. 检测完成 ---\n');
fprintf('总测试次数: %d\n', total_tests);
fprintf('函数值不一致次数 (FAIL_FVAL): %d\n', failed_tests_fval);
fprintf('梯度不一致次数 (FAIL_GRAD): %d\n', failed_tests_grad);
if failed_tests_fval == 0 && failed_tests_grad == 0
    fprintf('所有检测通过！\n');
else
    fprintf('检测存在失败项，请检查！\n');
end


%% Helper functions for testing intermediate value calculations (if needed)
% These should exactly mirror the logic within your fmincon_objective_helper_proportional
% (the one you trust for function values)

function actual_k_prime = calculate_actual_k_prime_for_test(k_prime_proportion, resources_after_pps, cS_local)
    consumption_floor_spending = cS_local.cFloor * (1 + cS_local.tau_c);
    resources_for_kprime_and_c_above_floor = resources_after_pps - consumption_floor_spending;
    
    k_prime_val_calc = 0;
    if resources_for_kprime_and_c_above_floor >= 0
        k_prime_val_calc = k_prime_proportion * resources_for_kprime_and_c_above_floor;
        k_prime_val_calc = max(cS_local.kMin, min(k_prime_val_calc, resources_for_kprime_and_c_above_floor));
    else
        k_prime_val_calc = cS_local.kMin;
    end
    actual_k_prime = max(cS_local.kMin, min(k_prime_val_calc, cS_local.kMax));
end

function current_c = calculate_current_c_for_test(actual_k_prime, resources_after_pps, cS_local)
    consumption_expenditure = resources_after_pps - actual_k_prime;
    current_c_unclamped = consumption_expenditure / (1 + cS_local.tau_c);
    current_c = max(cS_local.cFloor, current_c_unclamped);
end

function k_pps_prime = calculate_k_pps_prime_for_test(actual_c_pps, k_pps_state, a_idx, R_k_net_factor_age, cS_local)
    pps_withdrawal = 0;
    annual_age_check = cS_local.physAgeMap{a_idx}(1); % Assuming cS.physAgeMap is available
    is_retired = (a_idx > cS_local.aR_new);
    if is_retired && annual_age_check >= cS_local.pps_withdrawal_age_min_idx && cS_local.pps_active
        pps_withdrawal = k_pps_state * cS_local.pps_withdrawal_rate;
    end
    pps_return_factor = 1 + ((R_k_net_factor_age - 1) + cS_local.pps_return_rate_premium);
    k_pps_prime_unclamped = (k_pps_state + actual_c_pps - pps_withdrawal) * pps_return_factor;
    k_pps_prime = max(cS_local.kppsMin, min(cS_local.kppsMax, k_pps_prime_unclamped));
end