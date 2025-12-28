% --- START OF FILE main_olg_v1.m ---

% 扩展的Huggett(1996)模型
% 加入老龄化人口结构、公共养老金和个人养老金系统
% 基于Huggett_1996_Main_New.m和OLG模型设计
% 使用 parfor 加速家户问题求解
% 基于 v3 utils (对齐 LaTeX 理论部分)

%--------------------------------------------------------------------
%{
扩展模型特点:
(1). 老龄化年龄结构 (16组, 22-101岁)
(2). 公共养老金系统 (PAYG, 赤字累积)
(3). 个人养老金账户 (工作期缴费q, 退休期提取)
(4). 异质性: 劳动禀赋冲击 (AR(1))
(5). 意外遗产: 死亡代理人资产转移给新生代
%}

%% 环境设置
clc; clear; close all;
rng('default'); % For reproducibility

% --- Start Parallel Pool ---
numCores = 20; % Set desired number of cores
if isempty(gcp('nocreate'))
    fprintf('Starting parallel pool with %d workers...\n', numCores);
    try
        parpool(numCores);
    catch ME_pool
        warning('Could not start parallel pool with %d workers: %s. Using default.', numCores, ME_pool.message);
        try parpool(); catch, warning('Failed to start default parallel pool.'); end
    end
else
    disp('Parallel pool already exists.');
    poolobj = gcp('nocreate');
    fprintf('Existing pool has %d workers.\n', poolobj.NumWorkers);
end
% --------------------------

%% 固定模型参数设置
fprintf('=== 模型初始化 ===\n');
cS = main_olg_v1_utils.initParameters();
[paramS.leLogGridV, paramS.leTrProbM, paramS.leProb1V] = main_olg_v1_utils.EarningProcess_olgm(cS);
paramS.leGridV = exp(paramS.leLogGridV);
paramS.ageEffV = cS.h;
fprintf('Parameter Initialization Complete. kMax=%.1f, pMax=%.1f\n', cS.kMax, cS.pMax);

%% 人口结构动态
fprintf('\n=== 步骤1: 模拟人口结构动态过程 ===\n');
fprintf('最大模拟周期数: %d 期（寻找稳态）\n', cS.max_periods);
popS = main_olg_v1_utils.initPopulation(cS);
popS = main_olg_v1_utils.populationDynamics(popS, cS);
cS.nPeriods = size(popS.Z, 2); % Update actual periods simulated
[Z_ss, dep_ratio_ss, bgp_reached, bgp_period] = main_olg_v1_utils.detectSteadyStatePopulation(popS, cS);

fprintf('\n=== 人口动态模拟结果概要 ===\n');
fprintf('实际模拟期数: %d\n', cS.nPeriods);
if bgp_reached
    fprintf('人口稳态: 在第%d期达到\n', bgp_period);
    t_ss = bgp_period; % Use steady state period index
else
    fprintf('人口稳态: 未在%d期内检测到。\n', cS.max_periods);
    fprintf('警告: 均衡将在最后模拟期 (%d) 的人口结构下计算。\n', cS.nPeriods);
    t_ss = cS.nPeriods; % Use the last simulated period
end
fprintf('稳态人口年龄结构已绘制。\n');

% Ensure Z_ss corresponds to the target period t_ss
Z_ss = popS.Z(:, t_ss);
if isfield(popS, 'dependencyRatio') && length(popS.dependencyRatio) >= t_ss
    dep_ratio_ss = popS.dependencyRatio(t_ss); % Get correct dependency ratio
else
    % Recalculate if needed
    working_pop_ss = sum(Z_ss(1:cS.workingAgeMaxIdx));
    retired_pop_ss = sum(Z_ss(cS.retirementAgeIdx:end));
    if working_pop_ss > 0, dep_ratio_ss = retired_pop_ss / working_pop_ss; else dep_ratio_ss = inf; end
end


%% 求解模型均衡
fprintf('\n=== 步骤2: 求解经济均衡 (稳态) ===\n');

% 计算稳态总劳动供给
L_ss = main_olg_v1_utils.computeAggLabor(Z_ss, cS, paramS); % Z_ss is population counts for the target period
fprintf('稳态总劳动 L_ss: %.4f\n', L_ss);
if L_ss <= 0
    error('稳态总劳动供给 L_ss 为零或负，模型无法继续。请检查人口或效率参数。');
end

% 初始猜测 K (总私人资产 K+P), T (意外遗产), D (养老金累积赤字)
% Need realistic guesses, maybe based on K/Y or K/L ratios
Y_guess_rough = (cS.kMax/4)^cS.alpha * L_ss^(1-cS.alpha); % Rough guess assuming K is 1/4 of max grid
KGuess =  3.5 * Y_guess_rough; % Target K/Y around 3.5
TGuess = 0.05;  % Bequest per capita
DGuess = Y_guess_rough * 0.5;   % Initial debt guess as fraction of rough GDP
fprintf('Initial Guesses: K=%.2f, T=%.4f, D=%.2f\n', KGuess, TGuess, DGuess);


% 均衡迭代参数
maxIter_eq = cS.maxIter_eq;
tolLevel_eq = cS.tolLevel_eq;
dampK = cS.dampK;
dampT = cS.dampT;
dampD = cS.dampD;

iter = 0;
devNorm = inf;
K_history = zeros(1, maxIter_eq+1); K_history(1) = KGuess;
T_history = zeros(1, maxIter_eq+1); T_history(1) = TGuess;
D_history = zeros(1, maxIter_eq+1); D_history(1) = DGuess;

fprintf('\n迭代寻找均衡 (K* (总资产), T*, D*)\n');
fprintf('%4s | %12s | %12s | %12s | %12s | %8s\n', 'Iter', 'KGuess', 'TGuess', 'DGuess', 'Norm', 'Time (s)');
fprintf('----------------------------------------------------------------------\n');

% Store last successful results in case of non-convergence
last_K_private_Model = NaN; last_P_Model = NaN; last_TModel = NaN; last_DModel_next = NaN;
last_agg_C = NaN; last_agg_PP_contrib = NaN; last_R = NaN; last_w = NaN; last_b = NaN;
last_deficit = NaN; last_Y = NaN;
last_public_revenue = NaN; last_public_expenditure = NaN;

eq_converged = false;
iter_timer = tic;

for iter = 1:maxIter_eq
    iter_start_time = tic;
    fprintf('Iter %d Start: K=%.2f, T=%.4f, D=%.2f\n', iter, KGuess, TGuess, DGuess);

    % 1. 计算价格和养老金流 (给定 KGuess = K_firm)
    try
        [Y, R, w, b, deficit_flow, public_pension_revenue, public_pension_expenditure] = main_olg_v1_utils.computePrices(KGuess, L_ss, DGuess, Z_ss, cS);
        if ~all(isfinite([Y, R, w, b, deficit_flow, public_pension_revenue, public_pension_expenditure]))
            error('NaN or Inf encountered in computePrices.');
        end
        fprintf('  Prices: Y=%.2f, R=%.5f, w=%.4f, b=%.4f, deficit=%.4f\n', Y, R, w, b, deficit_flow);
    catch ME_prices
        fprintf('Error in computePrices at iter %d: %s\n', iter, ME_prices.message);
        fprintf('Stopping iteration.\n'); break;
    end

    % 2. 求解家户问题 (使用并行化 parfor)
    try
        [cPolM, kPolM, pPolM, qPolM, valueM] = main_olg_v1_utils.solveHouseholdProblem(R, w, TGuess, b, cS, paramS);
    catch ME_hh
        fprintf('Error in solveHouseholdProblem at iter %d: %s\n', iter, ME_hh.message);
        fprintf('Stopping iteration.\n'); break;
    end

    % 3. 计算模型宏观总量 (给定政策函数和价格)
    try
        [K_private_Model, P_Model, TModel, DModel_next, agg_C, agg_pension_contrib_private, ~, ~, ~] = main_olg_v1_utils.computeAggregates(kPolM, pPolM, qPolM, cPolM, Z_ss, TGuess, DGuess, deficit_flow, R, w, cS, paramS);
         if ~all(isfinite([K_private_Model, P_Model, TModel, DModel_next, agg_C, agg_pension_contrib_private]))
             % Check for NaNs that might come from computeAggregates if simulation failed badly
             error('NaN or Inf encountered in computeAggregates results.');
         end
    catch ME_agg
        fprintf('Error in computeAggregates at iter %d: %s\n', iter, ME_agg.message);
        fprintf('Stopping iteration.\n'); break;
    end

    % Store last valid results from this iteration
    last_K_private_Model = K_private_Model; last_P_Model = P_Model; last_TModel = TModel;
    last_DModel_next = DModel_next; last_agg_C = agg_C; last_agg_PP_contrib = agg_pension_contrib_private;
    last_R = R; last_w = w; last_b = b; last_deficit = deficit_flow; last_Y = Y;
    last_public_revenue = public_pension_revenue; last_public_expenditure = public_pension_expenditure;

    % 4. 计算偏差
    KDev = KGuess - K_private_Model;
    TDev = TGuess - TModel;
    DDev = DGuess - DModel_next;
    % Use relative deviation for K if KGuess is large, absolute otherwise?
    % Or just a weighted norm. Let's use absolute for now.
    devNorm = sqrt(KDev^2 + TDev^2 + DDev^2);

    iter_time = toc(iter_start_time);
    fprintf('%4d | %12.4f | %12.6f | %12.4f | %12.6e | %8.2f\n', iter, KGuess, TGuess, DGuess, devNorm, iter_time);
    fprintf('     | K_Mod=%10.4f P_Mod=%10.4f T_Mod=%10.6f D_Next=%10.4f\n', K_private_Model, P_Model, TModel, DModel_next);

    % 5. 检查收敛
    if devNorm < tolLevel_eq && isfinite(devNorm)
        fprintf('均衡已收敛! K*(k+p)=%.4f, T*=%.6f, D*=%.4f\n', K_private_Model, TModel, DModel_next); % Use model values at convergence
        eq_converged = true;
        % Store final converged values
        KGuess = K_private_Model; TGuess = TModel; DGuess = DModel_next;
        break;
    end

    % 6. 更新猜测值 (阻尼)
    KGuessNew = KGuess - dampK * KDev;
    TGuessNew = TGuess - dampT * TDev;
    DGuessNew = DGuess - dampD * DDev;

    % Add bounds/safeguards for guesses
    KGuess = max(cS.kMin * sum(Z_ss)*1.1, KGuessNew); % Ensure K > slightly above min level
    TGuess = max(0, TGuessNew); % Ensure T >= 0
    % Add a check for runaway D? Optional.
    % if abs(DGuessNew) > 10 * abs(Y), warning('Debt guess seems large relative to Y'); end
    DGuess = DGuessNew;

    % Store history for plotting
    K_history(iter+1) = KGuess;
    T_history(iter+1) = TGuess;
    D_history(iter+1) = DGuess;

    % Sanity check for oscillations or divergence
    if iter > 10 && devNorm > 1e6 % Arbitrary large number
        fprintf('Warning: Deviation norm large (%.2e), potential divergence.\n', devNorm);
        % Optional: break if diverging too much
        % break;
    end


end % End equilibrium iteration loop

total_eq_time = toc(iter_timer);
fprintf('Equilibrium search finished in %.2f seconds.\n', total_eq_time);

% Assign final results based on convergence status
if ~eq_converged
    fprintf('\n警告：达到最大迭代次数 %d 或遇到错误，但未收敛。\n', maxIter_eq);
    fprintf('最后偏差范数: %.6e\n', devNorm);
    if isnan(last_K_private_Model)
        error('Equilibrium iteration failed completely, no valid results stored.');
    else
        fprintf('使用最后一次成功迭代的结果。\n');
        K_private_Model = last_K_private_Model; P_Model = last_P_Model; TModel = last_TModel;
        DModel_next = last_DModel_next; agg_C = last_agg_C; agg_pension_contrib_private = last_agg_PP_contrib;
        R = last_R; w = last_w; b = last_b; deficit_flow = last_deficit; Y = last_Y;
        public_pension_revenue = last_public_revenue; public_pension_expenditure = last_public_expenditure;
        % Keep KGuess, TGuess, DGuess as they were at the start of the last successful iteration
        KGuess = K_history(iter); TGuess = T_history(iter); DGuess = D_history(iter);
    end
else
    % If converged, use the values from the end of the loop
    K_private_Model = KGuess; % Converged guess IS the model result
    P_Model = last_P_Model; % Use P_Model from the last iteration
    TModel = TGuess;
    DModel_next = DGuess;
    agg_C = last_agg_C; agg_pension_contrib_private = last_agg_PP_contrib;
    R = last_R; w = last_w; b = last_b; deficit_flow = last_deficit; Y = last_Y;
    public_pension_revenue = last_public_revenue; public_pension_expenditure = last_public_expenditure;
end

% 存储均衡结果
equilibriumS.K_private = K_private_Model; % Total private assets (k+p)
equilibriumS.K_firm = K_private_Model;    % Assume K_firm = K_private in eq.
equilibriumS.P_private = P_Model;         % Personal pension assets (p)
equilibriumS.K_nonPension = K_private_Model - P_Model; % k assets
equilibriumS.T = TModel;                  % Bequest per capita
equilibriumS.D = DModel_next;             % Equilibrium Debt level
equilibriumS.L = L_ss;
equilibriumS.R = R;
equilibriumS.w = w;
equilibriumS.b = b;
equilibriumS.deficit_flow = deficit_flow;
equilibriumS.public_pension_revenue = public_pension_revenue;
equilibriumS.public_pension_expenditure = public_pension_expenditure;
equilibriumS.Y = Y;
equilibriumS.Z_ss = Z_ss;
equilibriumS.dep_ratio_ss = dep_ratio_ss;
equilibriumS.agg_C = agg_C;

%% 计算稳态特征并生成图表
fprintf('\n=== 步骤3 & 4: 计算稳态特征并生成图表 ===\n');

% Re-solve HH problem and simulate *once* with equilibrium values to get final policies & distributions
fprintf('  Calculating final policies and distributions using equilibrium values...\n');
R_eq = equilibriumS.R; w_eq = equilibriumS.w; T_eq = equilibriumS.T;
b_eq = equilibriumS.b; D_eq = equilibriumS.D; K_eq = equilibriumS.K_firm;

if ~all(isfinite([R_eq, w_eq, T_eq, b_eq, D_eq, K_eq]))
    error('Final equilibrium values contain NaN or Inf. Cannot proceed.');
end

% Re-solve HH problem for final policies
[cPolM_eq, kPolM_eq, pPolM_eq, qPolM_eq, valueM_eq] = main_olg_v1_utils.solveHouseholdProblem(R_eq, w_eq, T_eq, b_eq, cS, paramS);

% Re-run simulation to get distributions corresponding to equilibrium policies
% Calculate final deficit flow using equilibrium values
[~, ~, ~, ~, deficit_flow_eq, ~, ~] = main_olg_v1_utils.computePrices(K_eq, L_ss, D_eq, Z_ss, cS);
% Run aggregation to get sim results
[~, ~, ~, ~, ~, ~, sim_k_eq, sim_p_eq, sim_e_idx_eq] = main_olg_v1_utils.computeAggregates(kPolM_eq, pPolM_eq, qPolM_eq, cPolM_eq, Z_ss, T_eq, D_eq, deficit_flow_eq, R_eq, w_eq, cS, paramS);

% --- Calculate average assets by age ---
avg_private_assets = zeros(cS.aD, 1); avg_k_assets = zeros(cS.aD, 1); avg_p_assets = zeros(cS.aD, 1);
fprintf('  Calculating average assets by age...\n');
for a_plot = 1:cS.aD % Use different loop var
    avg_assets_a = main_olg_v1_utils.computeAverageAssetsByAge(sim_k_eq, sim_p_eq, a_plot, cS);
    avg_private_assets(a_plot) = avg_assets_a;
    valid_k = isfinite(sim_k_eq(:,a_plot)); avg_k = 0; if any(valid_k), avg_k = mean(sim_k_eq(valid_k, a_plot)); end
    valid_p = isfinite(sim_p_eq(:,a_plot)); avg_p = 0; if any(valid_p), avg_p = mean(sim_p_eq(valid_p, a_plot)); end
    avg_k_assets(a_plot) = avg_k; avg_p_assets(a_plot) = avg_p;
end
fprintf('  Average assets by age calculated.\n');

%% 图表
fprintf('\n=== 图表生成 ===\n');

% 1. 迭代收敛图
figure('Name', 'Equilibrium Convergence');
plot_iter = 1:iter; % Use actual number of iterations run
if isempty(plot_iter), plot_iter=1; end % Handle case where loop didn't run
subplot(3,1,1); plot(plot_iter, K_history(2:iter+1), '-o'); title('K Guess Convergence'); ylabel('K (k+p)'); grid on; xlim([1 max(1,iter)]);
subplot(3,1,2); plot(plot_iter, T_history(2:iter+1), '-o'); title('T Guess Convergence'); ylabel('T (Bequest)'); grid on; xlim([1 max(1,iter)]);
subplot(3,1,3); plot(plot_iter, D_history(2:iter+1), '-o'); title('D Guess Convergence'); ylabel('D (Pension Debt)'); xlabel('Iteration'); grid on; xlim([1 max(1,iter)]);
sgtitle('Equilibrium Iteration Convergence');

% 2. 政策函数切片图
target_age_idx = 4; target_e_idx = floor(cS.nE/2) + 1; target_p_idx = floor(cS.nP/2) + 1;
target_e_idx = max(1, min(cS.nE, target_e_idx)); target_p_idx = max(1, min(cS.nP, target_p_idx));
fixed_p_val = cS.pGridV(target_p_idx); fixed_e_val = paramS.leGridV(target_e_idx);

figure('Name', sprintf('Policy Slices Age %d',target_age_idx));
subplot(2,2,1); plot(cS.kGridV, kPolM_eq(:, target_p_idx, target_e_idx, target_age_idx), 'b'); title('k'' vs k'); xlabel('Current k'); ylabel('Next k'''); grid on; xlim([cS.kGridV(1) cS.kGridV(end)]);
subplot(2,2,2); plot(cS.kGridV, pPolM_eq(:, target_p_idx, target_e_idx, target_age_idx), 'r'); title('p'' vs k'); xlabel('Current k'); ylabel('Next p'''); grid on; xlim([cS.kGridV(1) cS.kGridV(end)]);
subplot(2,2,3); plot(cS.kGridV, qPolM_eq(:, target_p_idx, target_e_idx, target_age_idx)*100, 'g'); title('q (%) vs k'); xlabel('Current k'); ylabel('q (%)'); grid on; xlim([cS.kGridV(1) cS.kGridV(end)]);
subplot(2,2,4); plot(cS.kGridV, cPolM_eq(:, target_p_idx, target_e_idx, target_age_idx), 'k'); title('c vs k'); xlabel('Current k'); ylabel('c'); grid on; xlim([cS.kGridV(1) cS.kGridV(end)]);
sgtitle(sprintf('Policy Function Slices (Age Group %d, p=%.2f, e=%.2f)', target_age_idx, fixed_p_val, fixed_e_val));

target_age_idx_retire = 12;
figure('Name', sprintf('Policy Slices Age %d',target_age_idx_retire));
subplot(1,2,1); plot(cS.kGridV, kPolM_eq(:, target_p_idx, target_e_idx, target_age_idx_retire), 'b'); title('k'' vs k'); xlabel('Current k'); ylabel('Next k'''); grid on; xlim([cS.kGridV(1) cS.kGridV(end)]);
subplot(1,2,2); plot(cS.kGridV, cPolM_eq(:, target_p_idx, target_e_idx, target_age_idx_retire), 'k'); title('c vs k'); xlabel('Current k'); ylabel('c'); grid on; xlim([cS.kGridV(1) cS.kGridV(end)]);
sgtitle(sprintf('Policy Function Slices (Age Group %d, p=%.2f, e=%.2f)', target_age_idx_retire, fixed_p_val, fixed_e_val));

% 3. 生命周期资产图
figure('Name', 'Lifecycle Asset Profiles');
age_midpoints = cS.age1 + 2.5 + (0:(cS.aD-1))*5;
plot(age_midpoints, avg_private_assets, 'b-o', 'LineWidth', 2, 'DisplayName', 'Total Private Assets (k+p)'); hold on;
plot(age_midpoints, avg_k_assets, 'k--x', 'LineWidth', 1.5, 'DisplayName', 'Capital Assets (k)');
plot(age_midpoints, avg_p_assets, 'm-.s', 'LineWidth', 1.5, 'DisplayName', 'Personal Pension Assets (p)');
xline(cS.age1 + (cS.retirementAgeIdx-1)*5, 'r--', 'DisplayName', 'Retirement Start'); hold off;
title('生命周期平均资产'); xlabel('年龄'); ylabel('平均资产'); legend('Location', 'best'); grid on;

%% 输出稳态结果
kl_ratio = 0; if equilibriumS.L > 0, kl_ratio = equilibriumS.K_firm / equilibriumS.L; end
kp_ratio_private = 0; if equilibriumS.L > 0, kp_ratio_private = equilibriumS.K_private / equilibriumS.L; end
k_to_p_ratio = Inf; if equilibriumS.P_private > 1e-9, k_to_p_ratio = equilibriumS.K_nonPension / equilibriumS.P_private; end
y_per_capita = 0; total_pop = sum(equilibriumS.Z_ss); if total_pop > 0, y_per_capita = equilibriumS.Y / total_pop; end
deficit_ratio_Y = 0; if abs(equilibriumS.Y) > 1e-9, deficit_ratio_Y = equilibriumS.deficit_flow / equilibriumS.Y; end
debt_ratio_Y = 0; if abs(equilibriumS.Y) > 1e-9, debt_ratio_Y = equilibriumS.D / equilibriumS.Y; end

fprintf('\n================模型稳态结果概要================\n');
fprintf('人口稳态期 (或模拟最终期): %d\n', t_ss);
fprintf('收敛状态: %s (Iterations: %d)\n', mat2str(eq_converged), iter);
fprintf('稳态人口抚养比: %.4f\n', equilibriumS.dep_ratio_ss);
fprintf('稳态总劳动供给 L: %.4f\n', equilibriumS.L);
fprintf('稳态总产出 Y: %.4f\n', equilibriumS.Y);
fprintf('稳态人均产出: %.4f\n', y_per_capita);
fprintf('稳态总私人资产 K*(k+p): %.4f\n', equilibriumS.K_private);
fprintf('  其中：资本资产 (k): %.4f\n', equilibriumS.K_nonPension);
fprintf('  其中：个人养老金资产 (p): %.4f\n', equilibriumS.P_private);
fprintf('  资本/个人养老金 资产比 (k/p): %.4f\n', k_to_p_ratio);
fprintf('稳态资本/劳动比 (K_firm/L): %.4f\n', kl_ratio);
fprintf('稳态总私人资产/劳动比 ((k+p)/L): %.4f\n', kp_ratio_private);
fprintf('稳态利率 R: %.5f (净利率 r: %.3f%%)\n', equilibriumS.R, (equilibriumS.R-1)*100);
fprintf('稳态工资率 w: %.4f\n', equilibriumS.w);
fprintf('稳态公共养老金福利 b: %.4f\n', equilibriumS.b);
fprintf('稳态意外遗产转移 T (人均): %.6f\n', equilibriumS.T);
fprintf('稳态公共养老金收入: %.4f\n', equilibriumS.public_pension_revenue);
fprintf('稳态公共养老金支出: %.4f\n', equilibriumS.public_pension_expenditure);
fprintf('稳态养老金赤字流/GDP: %.4f%%\n', deficit_ratio_Y*100);
fprintf('稳态累积养老金债务/GDP: %.4f%%\n', debt_ratio_Y*100);
fprintf('稳态总消费 C: %.4f\n', equilibriumS.agg_C);
I_res = equilibriumS.Y - equilibriumS.agg_C; % C+I = Y (assuming no G in model for now)
fprintf('稳态总投资 I (Y-C): %.4f\n', I_res);
% Sanity Check: Goods Market Clearing (Y approx C+I+Deficit?) - deficit not included in basic Y=C+I
fprintf('Sanity Check: Y - C - I = %.4f (should be close to 0 if no G)\n', equilibriumS.Y - equilibriumS.agg_C - I_res);
fprintf('==============================================\n');

% --- Clean up Parallel Pool ---
poolobj = gcp('nocreate');
if ~isempty(poolobj)
    disp('Shutting down parallel pool.');
    delete(poolobj);
end
% -----------------------------

fprintf('主程序运行完毕。\n');
% --- END OF FILE main_olg_v1.m ---