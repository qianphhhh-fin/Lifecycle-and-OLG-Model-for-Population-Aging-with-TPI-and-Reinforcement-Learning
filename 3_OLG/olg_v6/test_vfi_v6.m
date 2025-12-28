% --- test_vfi_v6.m ---
function test_vfi_v6() % Renamed main function
    clc; clear; 
    close all;
    fprintf('=== 测试 OLG 模型 V6 的 VFI - 检查价值函数对 k_pps 的敏感性 ===\n');
    fprintf('    (PPS所得税递延, VFI状态: k, k_pps, epsilon)\n');

    %% 1. 加载模型参数
    fprintf('\n--- 1. 加载参数 ---\n');
    cS_orig = main_olg_v6_utils.ParameterValues_HuggettStyle(); 

    % --- 选择一个激励最强的场景进行测试 ---
    scenario_to_debug = struct('name', 'Premium + Tax Deferral', 'pps_premium', 0.1, 'pps_tax_deferral', true);
    
    fprintf('\n\n========================================================================\n');
    fprintf('=== 排查场景: %s ===\n', scenario_to_debug.name);
    fprintf('========================================================================\n');

    cS = cS_orig; % Use a copy for scenario-specific changes
    cS.pps_return_rate_premium = scenario_to_debug.pps_premium;
        
    paramS = struct(); 
    [logGridV_temp, trProbM_temp, prob1V_temp] = main_olg_v6_utils.EarningProcess_olgm(cS);
    paramS.leLogGridV = logGridV_temp; 
    paramS.leGridV = exp(logGridV_temp); 
    paramS.leTrProbM = trProbM_temp; 
    paramS.leProb1V = prob1V_temp;       
    paramS.ageMassV = ones(cS.aD_new, 1) / cS.aD_new; 
    paramS.Z_ss_norm_annual = ones(cS.aD_orig, 1) / cS.aD_orig; 
    paramS.popGrowthForDebt = cS.popGrowth_orig;
    paramS.mass_workers_group = sum(paramS.ageMassV(1:cS.aR_new));
    if paramS.mass_workers_group < 1e-9 && cS.aR_new > 0, paramS.mass_workers_group = 0.5; end
    paramS.L_per_capita = 0.5; 
    paramS.ageEffV_new = cS.ageEffV_new;

    fprintf('参数: nk=%d, nkpps=%d, nw=%d, naD_new=%d\n', cS.nk, cS.nkpps, cS.nw, cS.aD_new);
    fprintf('PPS 回报溢价: %.3f, PPS 年度缴费上限: %.2f, PPS 缴费比例上限: %.2f\n', ...
        cS.pps_return_rate_premium, cS.pps_annual_contrib_limit, cS.pps_max_contrib_frac);
    fprintf('PPS 所得税递延在此场景中: %s\n', mat2str(scenario_to_debug.pps_tax_deferral));

    R_k_net_factor_test = 1.02; MPL_gross_test = 1.2; TR_total_test = 0.1; b_payg_per_retiree_test = 0.5; 
    bV_payg_test = zeros(1, cS.aD_new); if cS.aR_new < cS.aD_new, bV_payg_test(cS.aR_new + 1 : cS.aD_new) = b_payg_per_retiree_test; end
    
    paramS_for_vfi = paramS; % paramS for VFI call
    paramS_for_vfi.tau_l = 0.15; 
    paramS_for_vfi.theta_payg_actual_for_hh = 0.10; 
    paramS_for_vfi.pps_tax_deferral_active = scenario_to_debug.pps_tax_deferral; % Used by local HHIncome if that's called

    fprintf('测试参数: R_k_net=%.3f, MPL_gross=%.2f, TR_total=%.2f, b_payg=%.2f, tau_l=%.3f, theta_payg=%.3f\n', ...
            R_k_net_factor_test, MPL_gross_test, TR_total_test, b_payg_per_retiree_test, paramS_for_vfi.tau_l, paramS_for_vfi.theta_payg_actual_for_hh);

    %% 2. 执行完整的VFI求解 (调用 main_olg_v6_utils 中的版本)
    fprintf('\n--- 场景 [%s]: 执行完整的 main_olg_v6_utils.HHSolution_VFI_Huggett 以获取 V_{a+1} ---\n', scenario_to_debug.name);
    tic;
    % valM_full_scenario is V(k, k_pps, epsilon, age_group)
    [cPol_q_vfi_results, kPol_vfi_results, cPpsPol_vfi_results, valM_full_scenario] = ...
        main_olg_v6_utils.HHSolution_VFI_Huggett(R_k_net_factor_test, MPL_gross_test, ...
                                               TR_total_test, bV_payg_test, ...
                                               paramS_for_vfi, cS); % Pass current cS
    vfi_full_time = toc;
    fprintf('场景 [%s]: 完整VFI执行完毕。耗时: %.2f 秒。\n', scenario_to_debug.name, vfi_full_time);

    if any(isnan(valM_full_scenario(:))) || any(isinf(valM_full_scenario(:)))
        warning('场景 [%s]: 完整VFI结果中包含NaN或Inf值！后续排查可能不准确。', scenario_to_debug.name); 
        fprintf('DEBUG: valM_full_scenario for scenario "%s" contains NaN or Inf.\n', scenario_to_debug.name);
        num_nans = sum(isnan(valM_full_scenario(:))); num_infs = sum(isinf(valM_full_scenario(:)));
        fprintf('  Number of NaNs: %d\n', num_nans); fprintf('  Number of Infs: %d\n', num_infs);
        if num_nans > 0 && num_nans < 200 && numel(valM_full_scenario)>0, [r_nan,c_nan,p_nan,t_nan]=ind2sub(size(valM_full_scenario),find(isnan(valM_full_scenario))); fprintf('  NaNs at (k,kpps,e,age):\n'); disp([r_nan,c_nan,p_nan,t_nan]); end
        if num_infs > 0 && num_infs < 200 && numel(valM_full_scenario)>0, [r_inf,c_inf,p_inf,t_inf]=ind2sub(size(valM_full_scenario),find(isinf(valM_full_scenario))); fprintf('  Infs at (k,kpps,e,age):\n'); disp([r_inf,c_inf,p_inf,t_inf]); end
        return; 
    end

    %% 3. 检查价值函数 V(k, k_pps, eps, age) 对 k_pps 的敏感性
    fprintf('\n--- 场景 [%s]: 检查价值函数 V(k, k_pps, eps, age) 对 k_pps 的敏感性 ---\n', scenario_to_debug.name);

    ik_test_indices = round(linspace(1, cS.nk, min(3, cS.nk))); 
    if isempty(ik_test_indices) && cS.nk > 0, ik_test_indices = 1; elseif cS.nk == 0, ik_test_indices = []; end % Handle empty cS.nk
    ie_test_idx = round(cS.nw / 2); 
    if ie_test_idx == 0 && cS.nw > 0, ie_test_idx = 1; elseif cS.nw == 0, ie_test_idx = 0; end

    age_groups_to_plot = [1, min(round(cS.aR_new/2),cS.aD_new) , cS.aR_new, min(cS.aR_new+1, cS.aD_new), max(1, cS.aD_new-1), cS.aD_new];
    age_groups_to_plot = unique(age_groups_to_plot); 
    age_groups_to_plot(age_groups_to_plot > cS.aD_new | age_groups_to_plot < 1) = [];

    if ie_test_idx > 0 && ~isempty(ik_test_indices)
        fprintf('将为以下非PPS资产k的索引打印价值函数: '); fprintf('%d ', ik_test_indices); fprintf('\n');
        fprintf('将为劳动效率状态索引 %d (值 %.2f) 打印价值函数。\n', ie_test_idx, paramS_for_vfi.leGridV(ie_test_idx));

        for ia_idx_vf = age_groups_to_plot 
            fprintf('\n--- 年龄组 a = %d (真实年龄约 %d-%d) ---\n', ia_idx_vf, cS.physAgeV_new(ia_idx_vf), cS.physAgeV_new(ia_idx_vf)+cS.yearStep-1);
            fprintf('%-10s |', 'k_pps');
            for k_val_print_idx = 1:length(ik_test_indices)
                fprintf('%-12s |', sprintf('V(k=%.1f)', cS.kGridV(ik_test_indices(k_val_print_idx))));
            end
            fprintf('\n%s\n', repmat('-',1, 10 + length(ik_test_indices)*15 ));
            if cS.nkpps == 0 || isempty(cS.kppsGridV), fprintf('  没有定义PPS资产网格 (nkpps=%d)。跳过此年龄组的PPS敏感性检查。\n', cS.nkpps); continue; end

            for ikpps_idx_vf = 1:cS.nkpps 
                k_pps_val_vf = cS.kppsGridV(ikpps_idx_vf);
                fprintf('%-10.3f |', k_pps_val_vf);
                for ik_idx_vf_local_idx = 1:length(ik_test_indices)
                    actual_ik_idx = ik_test_indices(ik_idx_vf_local_idx);
                    value_point = valM_full_scenario(actual_ik_idx, ikpps_idx_vf, ie_test_idx, ia_idx_vf);
                    fprintf('%-12.3f |', value_point);
                end; fprintf('\n');
            end
        end
    else
        fprintf('无法进行价值函数敏感性检查，nw或nk可能为0。\n');
    end
    fprintf('\n--- 价值函数敏感性检查结束 ---\n');
    
    %% 4. （可选）手动测试特定 c_pps 选择 (使用此文件末尾定义的局部Bellman和HHIncome)
    % 此部分用于更细致的单点分析，可以根据需要取消注释和修改
    %{
    fprintf('\n--- 场景 [%s]: 排查特定状态下的PPS决策 (通过局部函数手动测试) ---\n', scenario_to_debug.name);
    a_idx_debug = min(round(cS.aR_new / 2) + 2, cS.aD_new); if a_idx_debug == 0, a_idx_debug = 1; end
    if a_idx_debug >= cS.aD_new && cS.aD_new > 1, a_idx_debug = cS.aD_new -1; elseif a_idx_debug >= cS.aD_new && cS.aD_new == 1, fprintf('只有一个年龄组，无法进行VFI排查的未来价值步骤。\n'); return; end
    k_idx_debug = round(cS.nk / 3); k_pps_idx_debug_start = 1; e_idx_debug = round(cS.nw / 2);                
    k_now_debug = cS.kGridV(k_idx_debug); k_pps_now_debug = 0; if cS.nkpps > 0 && ~isempty(cS.kppsGridV), k_pps_now_debug = cS.kppsGridV(k_pps_idx_debug_start); end
    epsilon_now_debug = paramS_for_vfi.leGridV(e_idx_debug); b_age_val_debug = bV_payg_test(a_idx_debug);
    fprintf('排查状态：AgeGroupIdx=%d (年:%d-%d), k_now=%.2f, k_pps_now=%.2f, epsilon_idx=%d (val=%.2f)\n', ...
        a_idx_debug, cS.physAgeV_new(a_idx_debug), cS.physAgeV_new(a_idx_debug)+cS.yearStep-1, k_now_debug, k_pps_now_debug, e_idx_debug, epsilon_now_debug);

    vPrime_kkppse_next_debug = valM_full_scenario(:,:,:,a_idx_debug+1);
    EV_for_interp_debug_val = zeros(cS.nk, cS.nkpps);
    for ik_next = 1:cS.nk, for ikpps_next = 1:cS.nkpps, expected_v_sum = 0; for ie_next = 1:cS.nw, expected_v_sum = expected_v_sum + vPrime_kkppse_next_debug(ik_next, ikpps_next, ie_next) * paramS_for_vfi.leTrProbM(e_idx_debug, ie_next); end; EV_for_interp_debug_val(ik_next, ikpps_next) = expected_v_sum; end; end
    EV_next_interpolant_debug_obj = [];
    if cS.nk > 1 && cS.nkpps > 1, EV_next_interpolant_debug_obj = griddedInterpolant({cS.kGridV, cS.kppsGridV}, EV_for_interp_debug_val, 'linear'); 
    elseif cS.nk > 1, EV_next_interpolant_debug_obj = griddedInterpolant(cS.kGridV, EV_for_interp_debug_val(:,1), 'linear');
    elseif cS.nkpps > 1, EV_next_interpolant_debug_obj = griddedInterpolant(cS.kppsGridV, EV_for_interp_debug_val(1,:)', 'linear');
    else, EV_next_interpolant_debug_obj = @(k_s, kp_s) EV_for_interp_debug_val(1,1); end
    
    % ... (之前的手动测试 c_pps 和 k_prime 的循环，确保调用的是此文件末尾的局部函数) ...
    % ... (例如 CallInterpolator_local_test, HHIncome_Huggett_local_test, BellmanForKprimeOnly_GivenCpps_local_test_verbose) ...
    %}

    fprintf('\n--- 场景 [%s] 排查结束 ---\n', scenario_to_debug.name);

fprintf('\n--- 所有场景排查结束 ---\n');

% end % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< END OF MAIN FUNCTION

% =========================================================================
% ==================== LOCAL HELPER FUNCTIONS FOR TESTING =================
% =========================================================================
% (HHIncome_Huggett_local_test, BellmanForKprimeOnly_GivenCpps_local_test_verbose, CallInterpolator_local_test, CES_utility_local_test)
% 请将上一轮回复中这四个局部函数的定义粘贴到这里。
% 为简洁起见，我在这里省略了它们，但它们对于运行 Section 4 的可选部分是必需的。
% 如果您主要关注 Section 3 的价值函数敏感性检查，则暂时不需要它们。

