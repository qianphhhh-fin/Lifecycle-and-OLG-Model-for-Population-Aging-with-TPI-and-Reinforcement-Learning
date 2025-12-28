% --- 文件开头部分不变 ---
clc; clear; close all;
fprintf('=== Huggett Model with Dynamic Population (V2 Simplified - Corrected b Calc) ===\n'); % 更新标题

%% 1. 固定模型参数与人口设置
fprintf('\n--- 1. Initializing Parameters (Mixed Annual/Group) ---\n');
cS = main_olg_v2_utils.ParameterValues_HuggettStyle();
paramS = struct(); % 初始化空的paramS
fprintf('Parameters loaded. nk=%d, nw=%d.\n', cS.nk, cS.nw);
fprintf('Annual ages: %d-%d (%d yrs). Groups: %d (%d working).\n', ...
        cS.age1_orig, cS.ageLast_orig, cS.aD_orig, cS.aD_new, cS.aR_new);
fprintf('Using ANNUAL beta = %.4f and 1-yr transition survival in VFI.\n', cS.beta);

%% 2. Simulate Population Dynamics to Steady State (Group Level)
fprintf('\n--- 2. Simulating Population Dynamics (Group Level) ---\n');
popS = main_olg_v2_utils.initPopulation(cS);
popS = main_olg_v2_utils.populationDynamics(popS, cS);
[Z_ss, dep_ratio_ss, bgp_reached, bgp_period] = main_olg_v2_utils.detectSteadyStatePopulation(popS, cS); % Z_ss 是稳态人口数量

fprintf('\n--- Population Simulation Summary ---\n');
fprintf('Actual simulation periods: %d\n', length(popS.totalPop)-1);
if bgp_reached
    fprintf('Population steady state (groups) reached at period: %d\n', bgp_period);
else
    fprintf('Population steady state (groups) NOT reached. Using final period %d.\n', bgp_period);
end
fprintf('Steady State Dependency Ratio (Counts Retire/Work): %.4f\n', dep_ratio_ss);

% 将稳态人口数量 Z_ss 存入 paramS 以便传递
paramS.Z_ss_counts = Z_ss; % <--- 新增，传递稳态人口数

% Normalize steady-state GROUP population distribution
Z_ss_total = sum(Z_ss);
Z_ss_norm_group = zeros(cS.aD_new,1);
if Z_ss_total > 1e-9
    Z_ss_norm_group = Z_ss / Z_ss_total;
end
paramS.ageMassV = Z_ss_norm_group(:); % Store GROUP mass distribution in paramS for L calculation (L的计算用占比)

% *** Approximate ANNUAL steady-state distribution from GROUP distribution ***
Z_ss_norm_annual = zeros(cS.aD_orig, 1);
if Z_ss_total > 1e-9 % 使用 Z_ss_total 判断
    for a_new = 1:cS.aD_new
        annual_indices = cS.physAgeMap{a_new};
        group_mass = Z_ss_norm_group(a_new);
        num_years_in_group = length(annual_indices);
        if num_years_in_group > 0
            mass_per_year = group_mass / num_years_in_group;
            Z_ss_norm_annual(annual_indices) = mass_per_year;
        end
    end
    if sum(Z_ss_norm_annual) > 1e-9
        Z_ss_norm_annual = Z_ss_norm_annual / sum(Z_ss_norm_annual);
    end
    fprintf('Derived approximate ANNUAL steady-state population distribution.\n');
else
    warning('Steady state population is zero. Cannot derive annual distribution.');
end


%% 3. Precompute Labor Supply and Endowment Process
fprintf('\n--- 3. Precomputing Labor ---\n');
[paramS.leLogGridV, paramS.leTrProbM, paramS.leProb1V] = main_olg_v2_utils.EarningProcess_olgm(cS);
paramS.leGridV = exp(paramS.leLogGridV(:));
fprintf('Labor endowment process calibrated (nw=%d states).\n', cS.nw);
paramS.ageEffV_orig = cS.ageEffV_orig;
paramS.ageEffV_new  = cS.ageEffV_new;
fprintf('Simulating ANNUAL labor endowments for %d individuals over %d years...\n', cS.nSim, cS.aD_orig);
eIdxM = main_olg_v2_utils.LaborEndowSimulation_olgm(cS, paramS);
[~, L] = main_olg_v2_utils.LaborSupply_Huggett(eIdxM, cS, paramS, Z_ss_norm_group); % L的计算使用归一化质量
fprintf('Aggregate Labor Supply (L_ss): %.4f\n', L);
if L <= 0 && sum(Z_ss(1:cS.aR_new)) > 1e-9 % 检查基于原始计数
    error('Aggregate Labor Supply is zero despite positive working population.');
elseif L <= 0
     warning('Aggregate Labor Supply is zero or negative. Setting to small positive value.');
     L = 1e-6;
end


%% 4. Solve for General Equilibrium (K*, T*) via Iteration
fprintf('\n--- 4. Solving for General Equilibrium (Mixed Time Units) ---\n');
KGuess =  29.8828; % 初始猜测值
TGuess = 0.5054; % 初始猜测值
maxIter = 300; % 增加迭代次数
tolLevel = 1e-5;
dampK = 0.1;  % 调整阻尼因子
dampT = 0.1;
iter = 0; devNorm = inf;
fprintf('Starting equilibrium iteration...\n');
fprintf('Iter |   K Guess  |  T Guess   |   K Model  |   T Model  |   K Dev    |   T Dev    |   Norm     | Time\n');
fprintf('------------------------------------------------------------------------------------------------------\n');
eq_converged = false; K_hist = zeros(maxIter+1,1); T_hist = zeros(maxIter+1,1); Dev_hist = zeros(maxIter+1,1);
K_hist(1) = KGuess; T_hist(1) = TGuess; Dev_hist(1) = devNorm;
KModel=NaN; TModel=NaN; % 初始化

for iter = 1:maxIter
    iter_start_time = tic;

    % Step 3a (Inside Loop): Compute ANNUAL prices given KGuess, L
    % 传递 paramS (其中包含 paramS.Z_ss_counts) 给价格函数
    [~, R_iter, w_iter, b_iter] = main_olg_v2_utils.HHPrices_Huggett(KGuess, L, cS, paramS); % <--- 修改调用

    bV_new = zeros(1, cS.aD_new);
    if cS.aR_new < cS.aD_new, bV_new(cS.aR_new + 1 : cS.aD_new) = b_iter; end
    [cPolM, kPolM, ~] = utils.HHSolution_VFI_Huggett(R_iter, w_iter, TGuess, bV_new, paramS, cS);
    % [kHistM, ~] = utils.HHSimulation_olgm(kPolM, cPolM, eIdxM, cS);
        % [cPolM, kPolM, cPpsPolM_dummy, ~] = main_olg_v2_utils.HHSolution_VFI_Huggett(R_iter, w_iter, TGuess, bV_new, paramS, cS); % HHSolution_VFI_Huggett现在也输出cPpsPolM
    % HHSimulation_olgm 现在需要 cPpsPolM 作为输入
    kPpsHistM_dummy = zeros(cS.nSim, cS.aD_orig); % 创建一个空的 kPpsHistM，因为这个版本没有它
    % 或者，如果 main_olg_v2_utils.HHSimulation_olgm 现在总是期望 cPpsPolM:
    cPpsPolM_for_sim = zeros(cS.nk, cS.nw, cS.aD_new); % 如果 pps_max_contrib_frac=0, 这个应该是全零
    if exist('cPpsPolM_dummy','var') && ~isempty(cPpsPolM_dummy) % 确保 cPpsPolM_dummy 被赋值
        cPpsPolM_for_sim = cPpsPolM_dummy;
    end
    [kHistM, ~, ~] = main_olg_v2_utils.HHSimulation_olgm(kPolM, cPpsPolM_for_sim, cPolM, eIdxM, R_iter, w_iter, TGuess, bV_new, paramS, cS); % 假设它返回 kHistM, kPpsHistM, cHistM
    KModel = mean(kHistM, 1) * Z_ss_norm_annual;
    KModel = max(1e-6, KModel);
    kprimeHistM = zeros(cS.nSim, cS.aD_orig);
    if cS.aD_orig > 1, kprimeHistM(:, 1:cS.aD_orig-1) = kHistM(:, 2:cS.aD_orig); end
    ageDeathMassV_annual = Z_ss_norm_annual(:) .* cS.d_orig(:);
    mean_bequest_value_by_age_annual = mean(kprimeHistM * R_iter, 1);
    TotalBequests = sum(mean_bequest_value_by_age_annual(:) .* ageDeathMassV_annual(:));
    TModel = TotalBequests / (1 + cS.popGrowth_orig);
    TModel = max(0, TModel);
    KDev = KGuess - KModel; TDev = TGuess - TModel;
    devNorm = sqrt(KDev^2 + TDev^2);
    iter_time = toc(iter_start_time);
    fprintf('%4d | %10.4f | %10.4f | %10.4f | %10.4f | %10.4f | %10.4f | %10.4e | %.2fs\n', ...
             iter, KGuess, TGuess, KModel, TModel, KDev, TDev, devNorm, iter_time);
     K_hist(iter+1) = KModel; T_hist(iter+1) = TModel; Dev_hist(iter+1) = devNorm;
    if devNorm < tolLevel
        fprintf('Equilibrium converged!\n'); eq_converged = true; KGuess = KModel; TGuess = TModel; break;
    end
    KGuess = KGuess - dampK * KDev; TGuess = TGuess - dampT * TDev;
    KGuess = max(1e-6, KGuess); TGuess = max(0, TGuess);
end
if ~eq_converged
    fprintf('\nWarning: Equilibrium did not converge after %d iterations.\n', maxIter);
    fprintf('Final Deviation Norm: %.4e\n', devNorm); KGuess = KModel; TGuess = TModel;
end
K_eq = KGuess; T_eq = TGuess;

% --- Step 5a: Recalculate final prices and b_eq using Z_ss_counts from paramS ---
[Y_eq, R_eq, w_eq, b_eq_final, MPL_eq_final] = main_olg_v2_utils.HHPrices_Huggett(K_eq, L, cS, paramS); % <--- 修改调用并获取MPL_eq_final

bV_eq_new = zeros(1, cS.aD_new);
if cS.aR_new < cS.aD_new, bV_eq_new(cS.aR_new + 1 : cS.aD_new) = b_eq_final; end
[cPolM_eq, kPolM_eq, valueM_eq] = utils.HHSolution_VFI_Huggett(R_eq, w_eq, T_eq, bV_eq_new, paramS, cS);
% [kHistM_eq, cHistM_eq] = utils.HHSimulation_olgm(kPolM_eq, cPolM_eq, eIdxM, cS); % Optional

%% 5. Analyze and Plot Results
fprintf('\n--- 5. Equilibrium Results and Plots ---\n');
fprintf('Equilibrium Aggregate Capital (K*): %.4f\n', K_eq);
fprintf('Equilibrium Aggregate Labor (L):    %.4f\n', L);
fprintf('Equilibrium Output (Y* annual):   %.4f\n', Y_eq);
fprintf('Equilibrium Interest Factor (R* annual): %.4f (Rate r*=%.4f)\n', R_eq, R_eq - 1);
fprintf('Equilibrium Wage Rate (w* annual):     %.4f\n', w_eq);
fprintf('Equilibrium SS Benefit (b* annual, using dynamic pop):    %.4f\n', b_eq_final); % 使用最终的b
fprintf('Equilibrium Bequest Transfer (T* annual): %.4f\n', T_eq);
fprintf('Equilibrium K/Y Ratio (annual): %.4f\n', K_eq / Y_eq);

% --- 计算并打印有效平均替代率 ---
total_retirees_from_Z_ss = sum(paramS.Z_ss_counts(cS.aR_new + 1 : cS.aD_new));
effective_replacement_rate_calc = 0;
if total_retirees_from_Z_ss > 1e-9 && L > 1e-9 % L是总有效劳动供给
    % 替代率 = (theta * MPL * L / 总退休人数) / MPL = theta * L / 总退休人数
    effective_replacement_rate_calc = (cS.theta * (L / total_retirees_from_Z_ss)) * 100;
end

fprintf('Effective Average PAYG Replacement Rate (theta * L_eff_total / Retiree_count_ss): %.2f%%\n', effective_replacement_rate_calc);
fprintf('  (Calculated using: cS.theta=%.2f, L_eff_total=%.4f, Retiree_count_ss=%.2f)\n', cS.theta, L, total_retirees_from_Z_ss);
fprintf('  (This means b_eq_final = %.2f * MPL_eq_final)\n', effective_replacement_rate_calc / 100);
fprintf('  (MPL_eq_final used for check: %.4f)\n', MPL_eq_final);


% --- 绘图部分不变 ---
% ... (Plot 1, 2, 3) ...

% --- Plots ---

% Plot 1: Convergence History
figure('Name', 'Equilibrium Convergence (V2 Modified)');
subplot(3,1,1); plot(0:iter, K_hist(1:iter+1), '-o'); title('K Convergence'); ylabel('Aggregate Capital'); grid on; xlim([0 iter]);
subplot(3,1,2); plot(0:iter, T_hist(1:iter+1), '-o'); title('T Convergence'); ylabel('Annual Bequest Transfer'); grid on; xlim([0 iter]);
subplot(3,1,3); semilogy(1:iter, Dev_hist(2:iter+1), '-o'); title('Deviation Norm Convergence'); ylabel('Norm'); xlabel('Iteration'); grid on; xlim([1 iter]);
sgtitle('Equilibrium Iteration Convergence (Mixed Time Units)');

% Plot 2: Policy Function Slices (Example: Group 4 ~ age 37-41)
plot_a_new = 4; % Example working age group index
if plot_a_new > cS.aD_new, plot_a_new = cS.aR_new; end % Ensure valid index
physAgeStart = cS.physAgeV_new(plot_a_new);
physAgeEnd = cS.physAgeV_orig(cS.physAgeMap{plot_a_new}(end));
plot_title_suffix = sprintf('Group %d (Age %d-%d)', plot_a_new, physAgeStart, physAgeEnd);

figure('Name', ['Policy Functions: ' plot_title_suffix]);
subplot(1, 2, 1); % Consumption policy (annual c)
plot(cS.kGridV, cPolM_eq(:, 1, plot_a_new), 'b-', 'LineWidth', 1, 'DisplayName', sprintf('e=%d (Low)',1)); hold on;
plot(cS.kGridV, cPolM_eq(:, round(cS.nw/2), plot_a_new), 'k--', 'LineWidth', 1, 'DisplayName', sprintf('e=%d (Mid)',round(cS.nw/2)));
plot(cS.kGridV, cPolM_eq(:, cS.nw, plot_a_new), 'r:', 'LineWidth', 1, 'DisplayName', sprintf('e=%d (High)',cS.nw)); hold off;
xlabel('Start-of-Year Capital (k)'); ylabel('Optimal Annual Consumption (c)');
title(['Consumption Policy c(k, e, a_{new}) ' plot_title_suffix]);
legend('Location', 'northwest'); grid on;

subplot(1, 2, 2); % Saving policy (capital for next year k')
plot(cS.kGridV, kPolM_eq(:, 1, plot_a_new), 'b-', 'LineWidth', 1, 'DisplayName', sprintf('e=%d (Low)',1)); hold on;
plot(cS.kGridV, kPolM_eq(:, round(cS.nw/2), plot_a_new), 'k--', 'LineWidth', 1, 'DisplayName', sprintf('e=%d (Mid)',round(cS.nw/2)));
plot(cS.kGridV, kPolM_eq(:, cS.nw, plot_a_new), 'r:', 'LineWidth', 1, 'DisplayName', sprintf('e=%d (High)',cS.nw));
plot(cS.kGridV, cS.kGridV, 'g--', 'DisplayName', 'k'' = k'); % 45-degree line
hold off;
xlabel('Start-of-Year Capital (k)'); ylabel('Next Year Capital (k'')');
title(['Saving Policy k''(k, e, a_{new}) ' plot_title_suffix]);
legend('Location', 'northwest'); grid on;

% Plot 3: Lifecycle Profile of Average Capital (Start-of-Year)
% Requires simulation results kHistM_eq
if exist('kHistM_eq', 'var')
    avgK_byAge_annual = mean(kHistM_eq, 1); % Average capital holdings at the start of each ANNUAL age
    figure('Name', 'Lifecycle Capital Profile (Start-of-Year)');
    plot(cS.physAgeV_orig, avgK_byAge_annual, 'bo-', 'LineWidth', 1.5);
    hold on;
    xline(cS.ageRetire_orig, 'r--', 'LineWidth', 1, 'DisplayName', 'Retirement Starts');
    hold off;
    title('Average Start-of-Year Capital Holdings by Annual Age');
    xlabel('Physical Age');
    ylabel('Average Capital (k)');
    legend('Avg Capital', 'Retirement', 'Location', 'northwest');
    grid on;
    xlim([cS.age1_orig, cS.ageLast_orig]);
else
    fprintf('\nSkipping lifecycle plot (kHistM_eq not generated).\n');
end



% --- END OF FILE main_olg_v2_modified.m ---