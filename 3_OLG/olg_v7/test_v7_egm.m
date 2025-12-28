% --- START OF FILE test_v7_egm.m ---
% 目的: 比较基于fminbnd的VFI与EGM方法的求解速度和结果
% 修改: 在test脚本内部处理vPrime维度，不修改utils文件。
%       EGM部分为框架，需要填充核心逻辑。

clc;
clear;
close all;

fprintf('=== 测试和比较 VFI (fminbnd) 与 EGM ===\n');

%% 1. 初始化模型参数和基本设置
fprintf('\n--- 1. 初始化参数 ---\n');
cS = main_olg_v7_utils.ParameterValues_HuggettStyle();
paramS = struct(); % 用于存储派生参数

% --- 测试参数设置 ---
cS.nk = 30;     % 减小网格点以便快速测试
cS.nw = 5;
cS.nkpps = 1;   % 关键: 测试k_pps维度退化的情况
if cS.nkpps == 1 && (isempty(cS.kppsGridV) || length(cS.kppsGridV) ~=1)
    cS.kppsGridV = 0; % 如果nkpps=1, 设置一个标量kppsGridV
    fprintf('nkpps=1, cS.kppsGridV 已设置为 0。\n');
elseif cS.nkpps > 1 && (isempty(cS.kppsGridV) || length(cS.kppsGridV) ~= cS.nkpps)
    power_kpps_test = 1.5; % 与utils中一致
    kppsGridV_temp_test = cS.kppsMin + (cS.kppsMax - cS.kppsMin) * (linspace(0, 1, cS.nkpps).^power_kpps_test);
    if cS.nkpps > 0, kppsGridV_temp_test(1) = cS.kppsMin; end
    cS.kppsGridV = kppsGridV_temp_test(:);
    fprintf('nkpps>1, cS.kppsGridV 已重新生成。\n');
end
% cS.aD_new = 10; % 可以减少年龄组以加速测试，但会影响结果可比性

fprintf('参数已加载/调整。nk=%d, nkpps=%d, nw=%d, aD_new=%d\n', cS.nk, cS.nkpps, cS.nw, cS.aD_new);

% --- 确保 paramS 中包含必要参数 ---
if ~isfield(paramS, 'leGridV') || ~isfield(paramS, 'leTrProbM') || ~isfield(paramS, 'leProb1V') || ...
   isempty(paramS.leGridV) || isempty(paramS.leTrProbM) || isempty(paramS.leProb1V)
    fprintf('劳动禀赋参数缺失，正在计算...\n');
    [logGridV_temp, trProbM_temp, prob1V_temp] = main_olg_v7_utils.EarningProcess_olgm(cS);
    paramS.leLogGridV = logGridV_temp;
    paramS.leGridV = exp(logGridV_temp(:));
    paramS.leTrProbM = trProbM_temp;
    paramS.leProb1V = prob1V_temp;
    fprintf('劳动禀赋参数已计算并存入 paramS。\n');
end
if ~isfield(paramS, 'ageEffV_new') || isempty(paramS.ageEffV_new), paramS.ageEffV_new = cS.ageEffV_new; end

% --- 固定的经济环境 ---
R_k_net_factor_fixed = 1.02;
MPL_gross_fixed = 1.5;
TR_total_fixed = 0.05;
bV_payg_fixed = zeros(1, cS.aD_new);
if cS.aR_new < cS.aD_new, bV_payg_fixed(cS.aR_new+1:cS.aD_new) = 0.1; end

paramS_fixed_for_vfi = paramS;
paramS_fixed_for_vfi.tau_l = 0.1; % 假设所得税率
paramS_fixed_for_vfi.theta_payg_actual_for_hh = 0.1; % 假设PAYG税率
paramS_fixed_for_vfi.pps_tax_deferral_active = cS.pps_active; % 确保PPS递延状态正确

fprintf('固定经济环境: R_k_net=%.3f, MPL_gross=%.3f, TR=%.3f\n', ...
    R_k_net_factor_fixed, MPL_gross_fixed, TR_total_fixed);


%% 2. 使用现有的VFI方法 (基于fminbnd) 求解
fprintf('\n--- 2. 使用基于fminbnd的VFI求解 ---\n');
tic;

valM_fminbnd_full = -Inf(cS.nk, cS.nkpps, cS.nw, cS.aD_new);
cPolM_q_fminbnd_full  = zeros(cS.nk, cS.nkpps, cS.nw, cS.aD_new);
kPolM_fminbnd_full  = zeros(cS.nk, cS.nkpps, cS.nw, cS.aD_new);
cPpsPolM_rule_fminbnd_full = zeros(cS.nk, cS.nkpps, cS.nw, cS.aD_new);

for a_idx_fminbnd = cS.aD_new : -1 : 1
    vPrime_next_fminbnd_input = -Inf(cS.nk, cS.nkpps, cS.nw); % 总是创建3D输入
    if a_idx_fminbnd < cS.aD_new
        vPrime_next_fminbnd_input = valM_fminbnd_full(:,:,:,a_idx_fminbnd+1);
    end

    % 调用原始的HHSolutionByAge_VFI_Huggett_v7
    [c_age_temp, k_age_temp, cpps_age_temp, v_age_temp] = ...
        main_olg_v7_utils.HHSolutionByAge_VFI_Huggett_v7(...
        a_idx_fminbnd, vPrime_next_fminbnd_input, ...
        R_k_net_factor_fixed, MPL_gross_fixed, TR_total_fixed, bV_payg_fixed(a_idx_fminbnd), ...
        paramS_fixed_for_vfi, cS, paramS_fixed_for_vfi.leGridV);

    cPolM_q_fminbnd_full(:,:,:,a_idx_fminbnd) = c_age_temp;
    kPolM_fminbnd_full(:,:,:,a_idx_fminbnd) = k_age_temp;
    cPpsPolM_rule_fminbnd_full(:,:,:,a_idx_fminbnd) = cpps_age_temp;
    valM_fminbnd_full(:,:,:,a_idx_fminbnd) = v_age_temp;
end
time_fminbnd = toc;
fprintf('基于fminbnd的VFI求解完成，耗时: %.4f 秒\n', time_fminbnd);

%% 3. 使用EGM方法求解 (作为内嵌函数实现)
fprintf('\n--- 3. 使用EGM方法求解 (框架，假设k_pps已移出VFI状态) ---\n');
tic;

% --- 初始化EGM结果存储 (假设k_pps已移出状态，策略和价值是 nk x nw x aD_new) ---
cPolM_q_egm  = zeros(cS.nk, cS.nw, cS.aD_new);
kPolM_egm  = zeros(cS.nk, cS.nw, cS.aD_new);
valM_egm = -Inf(cS.nk, cS.nw, cS.aD_new);
cPpsPolM_rule_egm = zeros(cS.nk, cS.nw, cS.aD_new); % 规则化PPS缴费


% --- EGM 主循环 (逆向迭代) ---
for a_idx_egm_outer = cS.aD_new : -1 : 1
    fprintf('  EGM求解年龄组: %d\n', a_idx_egm_outer);
    val_func_next_input_egm = -Inf(cS.nk, cS.nw); % V_{a+1}(k',eps')
    if a_idx_egm_outer < cS.aD_new
        val_func_next_input_egm = valM_egm(:,:,a_idx_egm_outer+1);
    end

    % 调用内嵌的EGM求解函数
    [cPolM_q_egm(:,:,a_idx_egm_outer), ...
     kPolM_egm(:,:,a_idx_egm_outer), ...
     valM_egm(:,:,a_idx_egm_outer), ...
     cPpsPolM_rule_egm(:,:,a_idx_egm_outer)] = ...
        solve_age_EGM_nested(a_idx_egm_outer, val_func_next_input_egm, cS, paramS_fixed_for_vfi, ...
                             R_k_net_factor_fixed, MPL_gross_fixed, TR_total_fixed, bV_payg_fixed);
end
time_egm = toc;
fprintf('EGM求解完成，耗时: %.4f 秒\n', time_egm);


%% 4. 比较结果 (与之前类似)
fprintf('\n--- 4. 比较结果 ---\n');
compare_a_idx = round(cS.aD_new / 2);
compare_ie_idx = round(cS.nw / 2);
compare_ikpps_idx = 1; % 因为我们主要关注 nkpps=1 的情况来简化EGM

c_fminbnd = squeeze(cPolM_q_fminbnd_full(:, compare_ikpps_idx, compare_ie_idx, compare_a_idx));
c_egm     = squeeze(cPolM_q_egm(:, compare_ie_idx, compare_a_idx));
k_fminbnd = squeeze(kPolM_fminbnd_full(:, compare_ikpps_idx, compare_ie_idx, compare_a_idx));
k_egm     = squeeze(kPolM_egm(:, compare_ie_idx, compare_a_idx));
v_fminbnd = squeeze(valM_fminbnd_full(:, compare_ikpps_idx, compare_ie_idx, compare_a_idx));
v_egm     = squeeze(valM_egm(:, compare_ie_idx, compare_a_idx));

% ... (绘图和差异统计代码与上一个回复中的类似，这里省略以保持简洁) ...
% --- 比较消费策略 ---
figure('Name', 'EGM vs Fminbnd: Consumption');
plot(cS.kGridV, c_fminbnd, 'b-o', 'DisplayName', 'VFI (fminbnd)'); hold on;
plot(cS.kGridV, c_egm, 'r-x', 'DisplayName', 'EGM (框架)'); hold off;
xlabel('当前资产 k'); ylabel('消费 c');
title(sprintf('消费策略比较: Age %d, Eps %d', compare_a_idx, compare_ie_idx));
legend show; grid on;

% --- 比较储蓄策略 ---
figure('Name', 'EGM vs Fminbnd: Savings');
plot(cS.kGridV, k_fminbnd, 'b-o', 'DisplayName', 'VFI (fminbnd)'); hold on;
plot(cS.kGridV, k_egm, 'r-x', 'DisplayName', 'EGM (框架)');
plot(cS.kGridV, cS.kGridV, 'k--', 'DisplayName', 'k''=k'); hold off;
xlabel('当前资产 k'); ylabel('下一期资产 k''');
title(sprintf('储蓄策略比较: Age %d, Eps %d', compare_a_idx, compare_ie_idx));
legend show; grid on;

diff_c = abs(c_fminbnd - c_egm); diff_k = abs(k_fminbnd - k_egm);
fprintf('\n策略差异 (Age %d, Eps %d):\n', compare_a_idx, compare_ie_idx);
fprintf('  最大消费差异: %.4e, 平均消费差异: %.4e\n', max(diff_c), mean(diff_c));
fprintf('  最大储蓄差异: %.4e, 平均储蓄差异: %.4e\n', max(diff_k), mean(diff_k));


%% 5. 内嵌 EGM 求解函数和辅助函数

function [c_policy_age, k_policy_age, v_func_age, c_pps_rule_age_out] = ...
    solve_age_EGM_nested(a_idx_egm, val_func_next_ke_egm, cS_egm, paramS_egm, ...
                         R_k_net_egm, MPL_g_egm, TR_t_egm, bV_payg_egm)
    % 输入:
    %   a_idx_egm: 当前年龄组
    %   val_func_next_ke_egm: 下一期价值函数 V_{a+1}(k',eps'), 维度 (nk x nw)
    %   cS_egm, paramS_egm: 参数结构体
    %   R_k_net_egm, MPL_g_egm, TR_t_egm, bV_payg_egm(a_idx): 当前经济环境价格

    % --- 初始化本年龄组的结果 ---
    c_policy_age = zeros(cS_egm.nk, cS_egm.nw);
    k_policy_age = zeros(cS_egm.nk, cS_egm.nw);
    v_func_age   = -Inf(cS_egm.nk, cS_egm.nw);
    c_pps_rule_age_out = zeros(cS_egm.nk, cS_egm.nw);

    % --- 处理最后一个年龄组 (与fminbnd版本类似) ---
    if a_idx_egm == cS_egm.aD_new
        for ik_egm_last = 1:cS_egm.nk
            for ie_egm_last = 1:cS_egm.nw
                k_now_val_last = cS_egm.kGridV(ik_egm_last);
                epsilon_val_last = paramS_egm.leGridV(ie_egm_last);
                cpps_last_egm = 0; % 规则决定
                c_pps_rule_age_out(ik_egm_last, ie_egm_last) = cpps_last_egm;

                [resources_last, ~, ~] = main_olg_v7_utils.HHIncome_Huggett(...
                    k_now_val_last, R_k_net_egm, MPL_g_egm, TR_t_egm, bV_payg_egm(a_idx_egm), ...
                    cpps_last_egm, a_idx_egm, paramS_egm, cS_egm, epsilon_val_last);
                
                c_policy_age(ik_egm_last, ie_egm_last) = max(cS_egm.cFloor, resources_last / (1 + cS_egm.tau_c));
                k_policy_age(ik_egm_last, ie_egm_last) = cS_egm.kMin;
                [~, v_func_age(ik_egm_last, ie_egm_last)] = main_olg_v7_utils.CES_utility(...
                    c_policy_age(ik_egm_last, ie_egm_last), cS_egm.sigma, cS_egm);
            end
        end
        return;
    end

    % --- EGM 步骤 (适用于 a_idx_egm < cS_egm.aD_new) ---
    % 预计算下一期期望边际效用乘以回报因子等
    E_marg_util_next_times_R_at_k_prime_grid_egm = zeros(cS_egm.nk, cS_egm.nw); % (nk_prime_grid x nw_current_eps)
    for ie_now = 1:cS_egm.nw % 当前效率状态
        for ik_prime = 1:cS_egm.nk % 下一期资产网格点 k'_j
            k_prime_val = cS_egm.kGridV(ik_prime);
            expected_marg_util_next_sum = 0;
            for ie_next = 1:cS_egm.nw % 下一期效率状态
                % *** EGM核心步骤1: 获取下一期在(k_prime_val, ie_next)的消费和边际效用 ***
                % 这需要 c_policy_next(k',eps') 或从 V_next(k',eps')推导
                % 这里用占位符，需要您用真实逻辑替换
                if isempty(val_func_next_ke_egm) || any(isnan(val_func_next_ke_egm(:))) || any(isinf(val_func_next_ke_egm(:)))
                     % 如果下一期价值函数无效（例如第一个年龄组，或a_idx = cS.aD_new-1）
                     % 这里的逻辑需要根据您如何初始化V_{J}来决定，或者如何处理V_{J+1}（通常为0）
                     c_next_val_hypothetical = cS_egm.cFloor + 0.01; % 极简占位符
                else
                    % 假设有一个函数能从价值函数得到消费，或者上一轮EGM已算出c_policy_next
                    % c_next_val_hypothetical = get_c_from_v(k_prime_val, val_func_next_ke_egm(ik_prime, ie_next), ...);
                    % 为了能运行，使用一个非常粗略的近似：假设消费占资产的一部分
                    c_next_val_hypothetical = max(cS_egm.cFloor, k_prime_val * 0.1 + TR_t_egm + bV_payg_egm(a_idx_egm+1)); % 极度简化
                end
                marg_util_c_next = local_MarginalUtility_CES(c_next_val_hypothetical, cS_egm.sigma);
                expected_marg_util_next_sum = expected_marg_util_next_sum + ...
                    marg_util_c_next * paramS_egm.leTrProbM(ie_now, ie_next);
            end
            E_marg_util_next_times_R_at_k_prime_grid_egm(ik_prime, ie_now) = ...
                cS_egm.beta * cS_egm.s_1yr_transitionV(a_idx_egm) * R_k_net_egm * expected_marg_util_next_sum;
        end
    end

    % --- 对于每个当前效率状态 ie_now ---
    for ie_now = 1:cS_egm.nw
        % --- EGM核心步骤2: 从欧拉方程右端反解当期消费 ---
        rhs_euler = E_marg_util_next_times_R_at_k_prime_grid_egm(:, ie_now); % (nk x 1) for current ie_now
        c_now_on_k_prime_grid = local_inv_MarginalUtility_CES(rhs_euler, cS_egm.sigma); % (nk x 1)

        % --- 处理消费约束 ---
        c_now_on_k_prime_grid = max(cS_egm.cFloor, c_now_on_k_prime_grid);

        % --- EGM核心步骤3: 计算当期规则化PPS缴费和非资本收入 ---
        epsilon_val_egm_curr = paramS_egm.leGridV(ie_now);
        cpps_eff_contrib_curr = calculate_cpps_rule_egm(a_idx_egm, epsilon_val_egm_curr, MPL_g_egm, paramS_egm, cS_egm);
        non_capital_income_curr = get_non_capital_income_egm_local(a_idx_egm, epsilon_val_egm_curr, ...
            MPL_g_egm, TR_t_egm, bV_payg_egm(a_idx_egm), cpps_eff_contrib_curr, paramS_egm, cS_egm, R_k_net_egm);

        % --- EGM核心步骤4: 反解当期期末总资源 (end-of-period assets before k_next choice) ---
        %   end_of_period_total_assets_A = (1+tau_c)c_now + k_prime
        %   (k_prime_grid 是下一期资产选择的网格)
        end_of_period_total_assets_A = (1 + cS_egm.tau_c) .* c_now_on_k_prime_grid + cS_egm.kGridV; % (nk x 1)

        % --- EGM核心步骤5: 从期末总资源反解当期期初资产 k_now ---
        %   A = (1+r_net)k_now + non_capital_income_curr
        %   k_now = (A - non_capital_income_curr) / (1+r_net)
        k_now_endogenous_grid = (end_of_period_total_assets_A - non_capital_income_curr) / R_k_net_egm; % (nk x 1)

        % --- 清理和插值 (与之前框架类似，但要处理约束) ---
        %   *** EGM的关键步骤：处理借贷约束和插值 ***
        %   1. 找到 k_now_endogenous_grid < cS_egm.kMin 的点，这些点是借贷约束的。
        %      在这些点，c_now 需要重新计算，k_prime = cS_egm.kMin。
        %   2. 对于非约束的点，进行插值。
        %   此处为简化版，不完整包含严格的约束处理逻辑，仅做基本插值：
        valid_pts = ~isnan(k_now_endogenous_grid) & ~isinf(k_now_endogenous_grid) & (k_now_endogenous_grid >= cS_egm.kMin - 1e-6); % 添加kMin检查
        k_now_clean = k_now_endogenous_grid(valid_pts);
        c_now_clean = c_now_on_k_prime_grid(valid_pts);

        if length(k_now_clean) < 2
            warning('EGM test: Age %d, Eps %d, not enough valid endogenous grid points for interpolation.', a_idx_egm, ie_now);
            % Fallback: use a simple rule or previous fminbnd result for this slice
            temp_c_slice_fminbnd = squeeze(cPolM_q_fminbnd_full(:,cS_egm.nkpps==1,ie_now,a_idx_egm)); % 假设nkpps=1
            c_policy_age(:, ie_now) = temp_c_slice_fminbnd;
        else
            [k_now_sorted, sort_idx] = sort(k_now_clean);
            c_now_sorted = c_now_clean(sort_idx);
            [k_now_unique, unique_idx] = unique(k_now_sorted, 'stable');
            c_now_for_interp = c_now_sorted(unique_idx);

            if length(k_now_unique) < 2
                 warning('EGM test: Age %d, Eps %d, not enough unique endogenous grid points.', a_idx_egm, ie_now);
                 temp_c_slice_fminbnd = squeeze(cPolM_q_fminbnd_full(:,cS_egm.nkpps==1,ie_now,a_idx_egm));
                 c_policy_age(:, ie_now) = temp_c_slice_fminbnd;
            else
                c_policy_age(:, ie_now) = interp1(k_now_unique, c_now_for_interp, cS_egm.kGridV, 'linear'); % No 'extrap' for now
                % Handle extrapolation manually for points outside the k_now_unique range
                c_policy_age(cS_egm.kGridV < k_now_unique(1), ie_now) = c_now_for_interp(1);
                c_policy_age(cS_egm.kGridV > k_now_unique(end), ie_now) = c_now_for_interp(end);
                c_policy_age(isnan(c_policy_age(:,ie_now)), ie_now) = cS_egm.cFloor; % Fill any remaining NaNs
            end
        end
        c_policy_age(:, ie_now) = max(cS_egm.cFloor, c_policy_age(:, ie_now));

        % --- 计算储蓄策略和价值函数 (基于EGM得到的消费策略) ---
        for ik_final = 1:cS_egm.nk
            k_val_final = cS_egm.kGridV(ik_final);
            c_val_final = c_policy_age(ik_final, ie_now);
            cpps_final = calculate_cpps_rule_egm(a_idx_egm, paramS_egm.leGridV(ie_now), MPL_g_egm, paramS_egm, cS_egm);
            c_pps_rule_age_out(ik_final, ie_now) = cpps_final;
            
            non_cap_inc_final = get_non_capital_income_egm_local(a_idx_egm, paramS_egm.leGridV(ie_now), ...
                MPL_g_egm, TR_t_egm, bV_payg_egm(a_idx_egm), cpps_final, paramS_egm, cS_egm, R_k_net_egm);
            
            budget_final_egm = R_k_net_egm * k_val_final + non_cap_inc_final;
            k_policy_age(ik_final, ie_now) = budget_final_egm - (1 + cS_egm.tau_c) * c_val_final;
            k_policy_age(ik_final, ie_now) = max(cS_egm.kMin, min(cS_egm.kMax, k_policy_age(ik_final, ie_now)));

            [~, util_curr_egm] = main_olg_v7_utils.CES_utility(c_val_final, cS_egm.sigma, cS_egm);
            v_func_age(ik_final, ie_now) = util_curr_egm;
            if a_idx_egm < cS_egm.aD_new
                k_p_val_egm = k_policy_age(ik_final, ie_now);
                k_p_val_clamped_egm = max(cS_egm.kGridV(1), min(cS_egm.kGridV(end), k_p_val_egm));
                
                expected_v_next_sum_egm = 0;
                for ie_next_val = 1:cS_egm.nw
                    % 这里需要插值 val_func_next_ke_egm(k_p_val_clamped_egm, ie_next_val)
                    % 为简化，假设 k_p_val_clamped_egm 恰好是网格点或使用最近邻
                    idx_k_next = find(cS_egm.kGridV >= k_p_val_clamped_egm, 1, 'first');
                    if isempty(idx_k_next), idx_k_next = cS_egm.nk; end
                    
                    ev_next_val = -Inf;
                    if ~isempty(val_func_next_ke_egm) && idx_k_next <= size(val_func_next_ke_egm,1) && ie_next_val <= size(val_func_next_ke_egm,2)
                        ev_next_val = val_func_next_ke_egm(idx_k_next, ie_next_val);
                    end
                    if ~isfinite(ev_next_val), ev_next_val = -1e12; end % 增加isfinite检查

                    expected_v_next_sum_egm = expected_v_next_sum_egm + ...
                        ev_next_val * paramS_egm.leTrProbM(ie_now, ie_next_val);
                end
                 if ~isfinite(expected_v_next_sum_egm), expected_v_next_sum_egm = -1e12; end

                v_func_age(ik_final, ie_now) = util_curr_egm + ...
                    cS_egm.beta * cS_egm.s_1yr_transitionV(a_idx_egm) * expected_v_next_sum_egm;
            end
            if ~isfinite(v_func_age(ik_final, ie_now)), v_func_age(ik_final, ie_now) = -1e12; end
        end
    end % 结束 ie_now 循环
end % 结束 solve_age_EGM_nested 函数


% --- 辅助内嵌函数 for EGM ---
function margU = local_MarginalUtility_CES(cons_q, sigma_val_local)
    if abs(sigma_val_local - 1) < 1e-7
        margU = 1 ./ cons_q;
    else
        margU = cons_q .^ (-sigma_val_local);
    end
    margU(cons_q <= 1e-9) = 1e12; % 处理极小或零消费
    margU(~isfinite(margU)) = 1e12; % 处理NaN/Inf
end

function margU_inv = local_inv_MarginalUtility_CES(margU_val, sigma_val_local)
    if abs(sigma_val_local-0) < 1e-9, error('Sigma cannot be zero for inv_MarginalUtility_CES'); end
    margU_inv = margU_val.^(-1/sigma_val_local);
    margU_inv(margU_val <= 1e-12) = 1e10; % 对应无穷边际效用，消费极大
    margU_inv(~isfinite(margU_inv)) = 1e-6; % 处理NaN/Inf, 返回一个小的正消费
end

function cpps_contrib = calculate_cpps_rule_egm(a_idx_calc, epsilon_val_calc, MPL_g_calc, paramS_calc, cS_calc)
    cpps_contrib = 0;
    model_age_group_start_year_idx_calc = cS_calc.physAgeMap{a_idx_calc}(1);
    is_eligible_calc = (a_idx_calc <= cS_calc.aR_new && ...
                       model_age_group_start_year_idx_calc <= cS_calc.pps_contribution_age_max_idx && ...
                       cS_calc.pps_active && ...
                       (cS_calc.pps_max_contrib_frac > 0 || cS_calc.pps_annual_contrib_limit > 0) && ...
                       cS_calc.pps_fixed_contrib_schedule_frac(a_idx_calc) > 0);
    if is_eligible_calc
        age_eff_calc = cS_calc.ageEffV_new(a_idx_calc);
        lab_inc_calc = MPL_g_calc * age_eff_calc * epsilon_val_calc;
        contrib_rate_calc = cS_calc.pps_fixed_contrib_schedule_frac(a_idx_calc);
        cpps_desired_calc = lab_inc_calc * contrib_rate_calc;
        max_pps_income_calc = lab_inc_calc * cS_calc.pps_max_contrib_frac;
        cpps_contrib = min(cpps_desired_calc, cS_calc.pps_annual_contrib_limit);
        cpps_contrib = min(cpps_contrib, max_pps_income_calc);
        cpps_contrib = max(0, cpps_contrib);
    end
end

function non_cap_inc_val = get_non_capital_income_egm_local(a_idx_gni, epsilon_val_gni, ...
                                              MPL_g_gni, TR_t_gni, b_payg_val_gni, ...
                                              cpps_contrib_gni, paramS_gni, cS_gni, R_k_net_gni)
    % 这个函数需要返回的是预算约束中不依赖于当前k的部分:
    % (1+r_net)k_now + [ non_capital_income_at_current_period ]
    % non_capital_income_at_current_period = gross_labor_income - taxes - pps_contrib_as_expenditure + transfers + payg_benefits
    
    % 调用HHIncome_Huggett(k_now=0)得到的是:
    % resources_at_zero_k = (0) + (non_capital_income_gross_of_budget_item_cpps - cpps_contrib_as_budget_item)
    % 这正是我们需要的。
    [resources_at_zero_k_val, ~, ~] = main_olg_v7_utils.HHIncome_Huggett(...
        0, R_k_net_gni, MPL_g_gni, TR_t_gni, b_payg_val_gni, ...
        cpps_contrib_gni, a_idx_gni, paramS_gni, cS_gni, epsilon_val_gni);
    non_cap_inc_val = resources_at_zero_k_val;
end

% --- END OF FILE test_v7_egm.m ---