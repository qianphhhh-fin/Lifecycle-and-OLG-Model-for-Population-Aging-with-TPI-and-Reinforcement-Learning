% =========================================================================
% == SCRIPT: test_aggregation_logic.m
% == 版本:   [v10.0 - 劳动供给外生性终极检验]
% == 目的:   通过预计算理论总劳动，并与聚合结果对比，来最终确定
% ==         总劳动L在当前模型中究竟是外生参数还是内生变量。
% =========================================================================
clear; close all; clear classes;
addpath(pwd);
fprintf('=== 独立单元测试脚本 for calculate_aggregates_unified [v10.0] ===\n\n');

%% --- 1. 环境设定 (不变) ---
fprintf('--- 1. 正在构建求解终期稳态(ssF)所需的完整环境...\n');
ngrid = 100; ngrid_pps = 1;
cS = model_setup_utils_bgp.ParameterValues();
cS.nk = ngrid; cS.nkpps = ngrid_pps; cS.nkprime = ngrid; cS.npps = ngrid_pps;
cS = model_setup_utils_bgp.generateGrids(cS);
cS.endogenous_theta_mode = false;
cS.pps_active = false; cS.nw = 5; % <-- [重要] 激活劳动效率冲击以进行最严格的测试
cS.ss0_year = 2023; cS.start_year = 2023;
[Z_path, Z_path_raw, A_path, cS] = model_setup_utils_bgp.generate_exo_paths(cS, false);
cS = model_setup_utils_bgp.calcaulte_theta_payg_path(cS, false);
theta_path = cS.theta_path;
cSF = cS; cSF.pps_active = false; 
cSF.s_pathV = cS.s_pathV(:,end); 
total_pop_path = sum(Z_path_raw, 1);
pop_growth_rate_path_period = (total_pop_path(2:end) ./ total_pop_path(1:end-1)) - 1;
if isempty(pop_growth_rate_path_period); cSF.n_ss = 0; else; cSF.n_ss = (1+pop_growth_rate_path_period(end))^(1/cSF.time_Step)-1; end
cSF.g_A_ss = cS.g_A_ss; 
theta_for_ssF = theta_path(end);
Z_ss_norm_F = Z_path(:,end);
paramSF = struct();
fprintf('   [!] 正在激活劳动效率冲击 (nw=5) 以进行最严格测试...\n');
[paramSF.leGridV, paramSF.TrProbM_by_age, paramSF.leProb1V, cSF.nw_expanded] = model_setup_utils_bgp.EarningProcess_AgeDependent(cSF);
paramSF.leLogGridV = log(paramSF.leGridV(1:cSF.nw));
fprintf('   ✅ 环境构建完成。\n\n');

%% --- 2. 预计算理论总劳动 L_theoretical ---
fprintf('--- 2. 正在预计算理论上的总有效劳动 (L_theoretical)...\n');
% 理论总劳动 = sum over working ages { age_mass * age_efficiency * avg_efficiency_shock }
L_theoretical = 0;
% 对于新生儿，他们的效率冲击分布是leProb1V
% 对于其他年龄，其效率分布是上一代转移过来的。在稳态下，所有年龄组的效率分布都应收敛到长期稳态分布。
% 长期稳态效率分布是转移矩阵的特征向量
[leEigVec, leEigVal] = eig(paramSF.TrProbM_by_age{end}'); % 使用最后一个年龄段的转移矩阵（或者任何一个，因为它们应该相同）
[~, max_idx] = max(diag(leEigVal));
le_stationary_dist = abs(leEigVec(:, max_idx) / sum(leEigVec(:, max_idx)));

avg_efficiency_shock = sum(paramSF.leGridV(1:cS.nw) .* le_stationary_dist);
fprintf('   - 稳态下的平均效率冲击值: %.6f\n', avg_efficiency_shock);

for ia = 1:cSF.aR_new
    L_theoretical = L_theoretical + Z_ss_norm_F(ia) * cSF.ageEffV_new(ia) * avg_efficiency_shock;
end
fprintf('   - ✅ 预计算出的理论总劳动 L_theoretical: %.8f\n\n', L_theoretical);


%% --- 3. 设定与理论值不符的输入 ---
fprintf('--- 3. 设定与理论值不符的宏观输入...\n');
r_guess_F   = 0.04; 
L_guess_F   = L_theoretical + 0.1; % **[核心测试]** 故意使用一个错误的L_guess
Kg_guess_F  = 1.5;      
Beq_guess_F = 0.02;     
Kp_from_r = main_steady_state_utils_bgp.get_Kp_from_r(r_guess_F, Kg_guess_F, L_guess_F, cSF);
fprintf('   - 使用的L_guess (故意错误): %.8f\n', L_guess_F);
fprintf('   - 其他宏观输入: Kp=%.4f, Kg=%.4f, Beq=%.4f\n\n', Kp_from_r, Kg_guess_F, Beq_guess_F);

%% --- 4. [核心] 执行单次计算 ---
fprintf('\n--- 4. 启动单次聚合计算以检验 L_agg...\n');
[ssF, ~] = calculate_aggregates_unified_instrumented(...
    Kp_from_r, Kg_guess_F, L_guess_F, Beq_guess_F, Z_ss_norm_F, cSF, paramSF);

%% --- 5. 最终诊断分析 ---
if ~isempty(ssF)
    fprintf('\n--- 5. 最终诊断分析 ---\n');
    L_agg_result = ssF.L_hat;
    diff = L_agg_result - L_theoretical;
    
    fprintf('   - 预计算的理论总劳动 L_theoretical: %.8f\n', L_theoretical);
    fprintf('   - 最终聚合出的总劳动 L_agg:         %.8f\n', L_agg_result);
    fprintf('   - 两者差异 (L_agg - L_theoretical): %.4e\n', diff);

    if abs(diff) < 1e-6
         fprintf('\n   [结论]: ✅✅✅ 假设成立！L_agg 与 L_theoretical 完全一致！✅✅✅\n\n');
         fprintf('   这无可辩驳地证明了:\n');
         fprintf('   1. 总劳动L在您的模型中是外生给定的常数。\n');
         fprintf('   2. `L_agg` 的聚合计算是正确的。\n');
         fprintf('   3. 根本问题在于使用了 `5x5` 求解器去求解一个实际上是 `4x4` 的问题。\n');
         fprintf('   [最终行动方案]: 必须将求解框架重构为 `4x4` 系统。\n');
    else
         fprintf('\n   [结论]: 💥 假设被推翻！L_agg 与理论值不符。💥\n');
         fprintf('   这意味着 `solve_steady_state_distribution_unified` 函数的输出 `Dist`\n');
         fprintf('   内生地依赖于价格，导致了 `L_agg` 的变化。\n');
    end
else
    fprintf('--- 5. 函数执行失败 (返回为空) ---\n');
end

%% ================================================================
%% ==           本地化的聚合与核算函数 (与v8.0相同)             ==
%% ================================================================

function [ss, Dist] = calculate_aggregates_unified_instrumented(K_private_total_guess, K_g_guess, L_guess, Bequest_Total_guess, Z_ss_norm, cS, paramS)
    if K_private_total_guess <= 0, K_private_total_guess = 1e-8; end; if K_g_guess <= 0, K_g_guess = 1e-8; end; if L_guess <= 0, L_guess = 1e-8; end;
    M_prices = main_steady_state_utils_bgp.get_prices_at_t(K_private_total_guess, K_g_guess, L_guess, 1.0, cS);
    M_for_hh = M_prices; M_for_hh.w_t = M_prices.w_hat_t; M_for_hh.r_mkt_t = M_prices.r_mkt_t;
    theta_ss = cS.theta_path(1);
    Household_Factor_Income_Gross = (M_prices.w_hat_t * L_guess) + (M_prices.r_mkt_t * K_private_total_guess);
    total_pension_pot_guess = theta_ss * M_prices.w_hat_t * L_guess;
    mass_retirees_ss = sum(Z_ss_norm((cS.aR_new+1):end));
    M_for_hh.b_t = total_pension_pot_guess / max(1e-9, mass_retirees_ss);
    M_for_hh.theta_t = theta_ss;
    total_population_mass = sum(Z_ss_norm(:));
    M_for_hh.beq_transfer_pers = Bequest_Total_guess / max(1e-12, total_population_mass);
    
    [polS, ~] = main_steady_state_utils_bgp.HHSolution_VFI_unified(M_for_hh, paramS, cS);

    Dist = main_steady_state_utils_bgp.solve_steady_state_distribution_unified(polS, paramS, cS, Z_ss_norm);
    
    L_agg = 0;
    for ia = 1:cS.aD_new
        if ia <= cS.aR_new
            mass_by_epsilon = squeeze(sum(Dist(:,:,:,ia), [1,2]));
            % 这里假设paramS.leGridV的维度是(nw_expanded, 1)或(1, nw_expanded)
            L_agg = L_agg + sum( (cS.ageEffV_new(ia) .* paramS.leGridV') .* mass_by_epsilon' , 'all');
        end
    end
    
    ss = struct();
    ss.L_hat = L_agg;
end