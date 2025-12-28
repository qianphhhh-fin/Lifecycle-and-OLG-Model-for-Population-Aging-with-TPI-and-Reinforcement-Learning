% --- START OF FILE main_olg_v2_pps.m ---

% Huggett (1996) 复制，包含动态人口结构
% 修改版：采用“混合时间单位”，包含个人养老金计划 (PPS),
% 并且PAYG工资税率内生，替代率外生

clc; clear; close all;
fprintf('=== Huggett 模型: 动态人口, PPS, 内生 PAYG 税率 ===\n');

%% 1. 固定模型参数与人口设置
fprintf('\n--- 1. 初始化参数 ---\n');
cS = main_olg_v2_utils_pps.ParameterValues_HuggettStyle(); % cS: 包含所有常数和参数的结构体
paramS = struct(); % paramS: 包含派生参数的结构体

fprintf('参数已加载。nk=%d (非PPS资产网格点数), nw=%d (劳动效率状态数)。\n', cS.nk, cS.nw);
fprintf('年度年龄范围: %d-%d (%d 年)。年龄组数: %d (%d 个工作年龄组)。\n', ...
        cS.age1_orig, cS.ageLast_orig, cS.aD_orig, cS.aD_new, cS.aR_new);
fprintf('VFI 中使用年度 beta = %.4f。外生 PAYG 替代率: %.2f\n', cS.beta, cS.pension_replacement_rate);
if cS.pps_active
    fprintf('个人养老金计划 (PPS) 已激活。\n');
    fprintf('  最大缴费比例: %.4f, 领取期税率: %.2f, 回报率溢价: %.3f, 领取率: %.2f\n', ...
        cS.pps_max_contrib_frac, cS.pps_tax_rate_withdrawal, cS.pps_return_rate_premium, cS.pps_withdrawal_rate);
    fprintf('  PPS 计入总资本 K: %d, PPS 可遗赠: %d\n', cS.pps_in_K, cS.pps_bequeathable);
else
    fprintf('个人养老金计划 (PPS) 未激活。\n');
end

%% 2. 模拟人口动态至稳态 (年龄组层面)
fprintf('\n--- 2. 模拟人口动态 (年龄组层面) ---\n');
popS = main_olg_v2_utils_pps.initPopulation(cS);
popS = main_olg_v2_utils_pps.populationDynamics(popS, cS);
[Z_ss, dep_ratio_ss, bgp_reached, bgp_period] = main_olg_v2_utils_pps.detectSteadyStatePopulation(popS, cS);

fprintf('\n--- 人口模拟总结 ---\n');
fprintf('实际模拟期数: %d\n', length(popS.totalPop)-1);
if bgp_reached, fprintf('人口稳态 (年龄组) 在第 %d 期达到。\n', bgp_period);
else, fprintf('人口稳态 (年龄组) 未达到。使用第 %d 期的最终数据。\n', bgp_period); end
fprintf('稳态抚养比 (退休人口数/工作人口数): %.4f\n', dep_ratio_ss);

paramS.Z_ss_counts = Z_ss; % 将稳态人口数量存入paramS

Z_ss_total = sum(Z_ss);
Z_ss_norm_group = zeros(cS.aD_new, 1);
if Z_ss_total > 1e-9, Z_ss_norm_group = Z_ss / Z_ss_total; end
paramS.ageMassV = Z_ss_norm_group(:);

Z_ss_norm_annual = zeros(cS.aD_orig, 1);
if Z_ss_total > 1e-9
    for a_new = 1:cS.aD_new, annual_indices = cS.physAgeMap{a_new}; group_mass = Z_ss_norm_group(a_new); num_years_in_group = length(annual_indices); if num_years_in_group > 0, mass_per_year = group_mass / num_years_in_group; Z_ss_norm_annual(annual_indices) = mass_per_year; end; end
    if sum(Z_ss_norm_annual) > 1e-9, Z_ss_norm_annual = Z_ss_norm_annual / sum(Z_ss_norm_annual); end
    fprintf('已推导出近似的年度稳态人口分布。\n');
else, warning('稳态人口为零，无法推导年度分布。'); Z_ss_norm_annual(:) = 1/cS.aD_orig; end % 如果人口为0，设为均匀分布以避免后续错误

%% 3. 预计算劳动供给和禀赋过程
fprintf('\n--- 3. 预计算劳动 ---\n');
[paramS.leLogGridV, paramS.leTrProbM, paramS.leProb1V] = main_olg_v2_utils_pps.EarningProcess_olgm(cS);
paramS.leGridV = exp(paramS.leLogGridV(:));
fprintf('劳动禀赋过程已校准 (nw=%d 个状态)。\n', cS.nw);
paramS.ageEffV_orig = cS.ageEffV_orig; paramS.ageEffV_new = cS.ageEffV_new;
fprintf('为 %d 个个体模拟 %d 年的年度劳动禀赋...\n', cS.nSim, cS.aD_orig);
eIdxM = main_olg_v2_utils_pps.LaborEndowSimulation_olgm(cS, paramS);
[~, L] = main_olg_v2_utils_pps.LaborSupply_Huggett(eIdxM, cS, paramS, Z_ss_norm_group);
fprintf('总劳动供给 (L_ss, 效率单位): %.4f\n', L);
if L <= 0 && sum(paramS.Z_ss_counts(1:cS.aR_new)) > 1e-9 % 使用paramS中的Z_ss_counts
    error('尽管存在正的工作人口，总劳动供给仍为零。');
elseif L <= 0, warning('总劳动供给为零或负。'); L = 1e-6; end

%% 4. 通过迭代求解一般均衡 (K*, T*, theta_payg*)
fprintf('\n--- 4. 求解一般均衡 (外生替代率, PPS) ---\n');

% 初始猜测值 (基于之前能收敛的简化模型NoPPS的结果，并为PPS留出一些空间)

% cS.pension_replacement_rate = 0.30;
KGuess = 27.3262 ;  % 总生产性资本的猜测值
TGuess = 0.1747;   % 年度 lump-sum 转移支付的猜测值

% 迭代参数
maxIter = 200; % 增加迭代次数
tolLevel = 1e-4;
dampK = 0.1;   % 使用较小的阻尼因子(如果cS.pension_replacement_rate太大了，工资需要用于PAYG的比例太高，需要较小的dampK才能收敛）
dampT = 0.1;
iter = 0; devNorm = inf;

fprintf('开始均衡迭代...\n');
fprintf('Iter |   K Guess  |  T Guess   | R_market  | ThetaPAYG | K Mod N-P | K Mod PPS  |  T Model   |   K Dev    |   T Dev    |   Norm     | Time\n');
fprintf('---------------------------------------------------------------------------------------------------------------------------------------------------\n');

eq_converged = false;
K_hist = zeros(maxIter+1,1); T_hist = zeros(maxIter+1,1); Dev_hist = zeros(maxIter+1,1);
R_mkt_hist = zeros(maxIter+1,1); Theta_payg_hist = zeros(maxIter+1,1); % 记录这些价格
K_hist(1) = KGuess; T_hist(1) = TGuess; Dev_hist(1) = devNorm; R_mkt_hist(1)=NaN; Theta_payg_hist(1)=NaN;
KModel_nonpps = NaN; KModel_pps = NaN; TModel = NaN;

for iter = 1:maxIter
    iter_start_time = tic;

    % --- 第 4a 步: 计算税前价格, PAYG福利目标, 和内生PAYG税率 ---
    [R_market_iter, MPL_gross_iter] = main_olg_v2_utils_pps.HHPrices_Huggett(KGuess, L, cS); % 只获取税前市场价格

    b_payg_target_iter = cS.pension_replacement_rate * MPL_gross_iter;
    b_payg_target_iter = max(0, b_payg_target_iter);

    retiree_count_iter = sum(paramS.Z_ss_counts(cS.aR_new + 1 : cS.aD_new)); % 使用paramS中存储的稳态人口数
    total_pension_outlay_target_iter = b_payg_target_iter * retiree_count_iter;
    total_gross_wage_bill_iter = MPL_gross_iter * L;

    theta_payg_needed_iter = 0;
    if total_gross_wage_bill_iter > 1e-9
        theta_payg_needed_iter = total_pension_outlay_target_iter / total_gross_wage_bill_iter;
    else
        if total_pension_outlay_target_iter > 1e-9, warning('Iter %d: Gross wage bill zero but pension outlay positive!', iter); theta_payg_needed_iter = cS.max_payg_payroll_tax_rate; end
    end
    theta_payg_actual_iter = max(0, min(theta_payg_needed_iter, cS.max_payg_payroll_tax_rate));

    w_net_iter = MPL_gross_iter * (1 - theta_payg_actual_iter); % 家庭面临的净工资
    w_net_iter = max(0, w_net_iter);

    bV_payg_for_vfi_iter = zeros(1, cS.aD_new);
    if cS.aR_new < cS.aD_new
        bV_payg_for_vfi_iter(cS.aR_new + 1 : cS.aD_new) = b_payg_target_iter;
    end

    % --- 第 4b 步: 求解家庭问题 (VFI) ---
    % VFI 的状态变量是 k (代表非PPS资产)
    [cPolM, kPolM, cPpsPolM, ~] = main_olg_v2_utils_pps.HHSolution_VFI_Huggett(R_market_iter, w_net_iter, TGuess, bV_payg_for_vfi_iter, paramS, cS);

    % --- 第 4c 步: 模拟家庭决策 ---
    [kHistM_non_pps, kPpsHistM, cHistM] = main_olg_v2_utils_pps.HHSimulation_olgm(kPolM, cPpsPolM, cPolM, eIdxM, R_market_iter, w_net_iter, TGuess, bV_payg_for_vfi_iter, paramS, cS);

    % --- 第 4d 步: 计算模型的总生产性资本 ---
    KModel_nonpps = mean(kHistM_non_pps, 1) * Z_ss_norm_annual;
    KModel_nonpps = max(1e-6, KModel_nonpps);
    KModel_pps = 0;
    if cS.pps_active && cS.pps_in_K && cS.pps_max_contrib_frac > 1e-9 % 仅当PPS实际贡献时
        KModel_pps = mean(kPpsHistM, 1) * Z_ss_norm_annual;
        KModel_pps = max(0, KModel_pps);
    end
    KModel = KModel_nonpps + KModel_pps;

    % --- 第 4e 步: 计算模型的年度总意外遗赠 (TModel) ---
    kprimeHistM_for_bequest = zeros(cS.nSim, cS.aD_orig);
    if cS.aD_orig > 1
        % 假设遗赠基于非PPS资产的期末价值
        kprimeHistM_for_bequest(:, 1:cS.aD_orig-1) = kHistM_non_pps(:, 2:cS.aD_orig);
    end
    ageDeathMassV_annual = Z_ss_norm_annual(:) .* cS.d_orig(:);
    mean_bequest_value = mean(kprimeHistM_for_bequest * R_market_iter, 1);
    TotalBequests = sum(mean_bequest_value(:) .* ageDeathMassV_annual(:));
    TModel = TotalBequests / (1 + cS.popGrowth_orig);
    TModel = max(0, TModel);

    % --- 第 4f 步: 计算偏差并检查收敛 ---
    KDev = KGuess - KModel;
    TDev = TGuess - TModel;
    devNorm = sqrt(KDev^2 + TDev^2);

    iter_time = toc(iter_start_time);
    fprintf('%4d | %10.4f | %10.4f | %9.4f | %9.4f | %10.4f | %10.4f | %10.4f | %10.4f | %10.4f | %10.4e | %.2fs\n', ...
             iter, KGuess, TGuess, R_market_iter, theta_payg_actual_iter, KModel_nonpps, KModel_pps, TModel, KDev, TDev, devNorm, iter_time);

     K_hist(iter+1) = KModel; T_hist(iter+1) = TModel; Dev_hist(iter+1) = devNorm;
     R_mkt_hist(iter+1) = R_market_iter; Theta_payg_hist(iter+1) = theta_payg_actual_iter;

    if devNorm < tolLevel
        fprintf('均衡已收敛!\n');
        eq_converged = true;
        KGuess = KModel; TGuess = TModel;
        break;
    end

    KGuess = KGuess - dampK * KDev;
    TGuess = TGuess - dampT * TDev;
    KGuess = max(1e-6, KGuess); TGuess = max(0, TGuess);
end

if ~eq_converged
    fprintf('\n警告: %d 次迭代后均衡未收敛。\n', maxIter);
    fprintf('最终偏差范数: %.4e\n', devNorm);
    if iter == maxIter % 确保使用最后一次迭代计算的值
        KGuess = KModel; TGuess = TModel;
    end
end

K_eq = KGuess; T_eq = TGuess;
theta_payg_eq = Theta_payg_hist(iter+1); % 获取最后一次迭代的theta

% --- 重新计算最终均衡价格 ---
[R_market_eq, MPL_gross_eq] = main_olg_v2_utils_pps.HHPrices_Huggett(K_eq, L, cS); % 注意：HHPrices_Huggett现在只返回这两个
Y_eq = cS.A * (K_eq^cS.alpha) * (L^(1-cS.alpha)); % 手动计算Y_eq
b_payg_eq = cS.pension_replacement_rate * MPL_gross_eq; b_payg_eq = max(0, b_payg_eq);
w_net_eq = MPL_gross_eq * (1 - theta_payg_eq); w_net_eq = max(0, w_net_eq);

bV_eq_new = zeros(1, cS.aD_new);
if cS.aR_new < cS.aD_new, bV_eq_new(cS.aR_new + 1 : cS.aD_new) = b_payg_eq; end
[cPolM_eq, kPolM_eq, cPpsPolM_eq, valueM_eq] = main_olg_v2_utils_pps.HHSolution_VFI_Huggett(R_market_eq, w_net_eq, T_eq, bV_eq_new, paramS, cS);

fprintf('模拟均衡状态下的最终分布...\n')
[kHistM_eq, kPpsHistM_eq, cHistM_eq] = main_olg_v2_utils_pps.HHSimulation_olgm(kPolM_eq, cPpsPolM_eq, cPolM_eq, eIdxM, R_market_eq, w_net_eq, T_eq, bV_eq_new, paramS, cS);
fprintf('最终模拟完成。\n')

%% 5. 分析和绘制结果
fprintf('\n--- 5. 均衡结果与绘图 (外生替代率 + PPS) ---\n');
K_nonpps_eq = mean(kHistM_eq, 1) * Z_ss_norm_annual;
K_pps_eq = 0; if cS.pps_active && cS.pps_max_contrib_frac > 1e-9, K_pps_eq = mean(kPpsHistM_eq, 1) * Z_ss_norm_annual; end
TotalAssets_eq = K_nonpps_eq + K_pps_eq;

fprintf('均衡总生产性资本 (K*): %.4f\n', K_eq);
fprintf('  分解: 非PPS资本 K: %.4f, PPS资本 K: %.4f\n', K_nonpps_eq, K_pps_eq);
fprintf('均衡家庭总资产 (非PPS+PPS): %.4f\n', TotalAssets_eq);
fprintf('均衡总劳动 (L):    %.4f\n', L);
fprintf('均衡总产出 (Y* 年度):   %.4f\n', Y_eq);
fprintf('均衡市场利率因子 (R_market* 年度): %.4f (利率 r_mkt*=%.4f)\n', R_market_eq, R_market_eq - 1);
fprintf('均衡总工资率 (MPL_gross*): %.4f\n', MPL_gross_eq);
fprintf('均衡内生PAYG工资税率 (theta_payg*): %.4f\n', theta_payg_eq);
fprintf('均衡净工资率 (w* PAYG税后): %.4f\n', w_net_eq);
fprintf('均衡PAYG福利 (b* 基于替代率): %.4f\n', b_payg_eq);
fprintf('均衡意外遗赠转移支付 (T* 年度): %.4f\n', T_eq);
fprintf('均衡 K/Y 比率 (年度): %.4f\n', K_eq / Y_eq);

% --- 绘图 ---
figure('Name', '均衡收敛过程 (外生替代率 + PPS)');
subplot(4,1,1); plot(0:iter, K_hist(1:iter+1), '-o'); title('K 收敛'); ylabel('K'); grid on; xlim([0 iter]);
subplot(4,1,2); plot(0:iter, T_hist(1:iter+1), '-o'); title('T 收敛'); ylabel('T'); grid on; xlim([0 iter]);
subplot(4,1,3); plot(0:iter, Theta_payg_hist(1:iter+1), '-o'); title('Theta PAYG 收敛'); ylabel('\theta_{PAYG}'); grid on; xlim([0 iter]);
subplot(4,1,4); semilogy(1:iter, Dev_hist(2:iter+1), '-o'); title('偏差范数收敛'); ylabel('Norm'); xlabel('迭代次数'); grid on; xlim([1 iter]);
sgtitle('均衡迭代收敛过程 (外生替代率 + PPS)');

plot_a_new = min(4, cS.aR_new);
if plot_a_new <= cS.aD_new && plot_a_new > 0 && ~isempty(cS.physAgeMap{plot_a_new})
    physAgeStart = cS.physAgeV_new(plot_a_new);
    physAgeEnd = cS.physAgeV_orig(cS.physAgeMap{plot_a_new}(end));
    plot_title_suffix = sprintf('年龄组 %d (年龄 %d-%d)', plot_a_new, physAgeStart, physAgeEnd);
    figure('Name', ['政策函数 (外生替代率 + PPS): ' plot_title_suffix]);
    subplot(1, 3, 1);
    plot(cS.kGridV, cPolM_eq(:, 1, plot_a_new), 'b-', 'LineWidth', 1, 'DisplayName', sprintf('e=%d (低)',1)); hold on;
    plot(cS.kGridV, cPolM_eq(:, round(cS.nw/2), plot_a_new), 'k--', 'LineWidth', 1, 'DisplayName', sprintf('e=%d (中)',round(cS.nw/2)));
    plot(cS.kGridV, cPolM_eq(:, cS.nw, plot_a_new), 'r:', 'LineWidth', 1, 'DisplayName', sprintf('e=%d (高)',cS.nw)); hold off;
    xlabel('期初非PPS资本 (k)'); ylabel('最优消费 (c)'); title(['消费政策 ' plot_title_suffix]); legend('Location', 'best'); grid on;
    subplot(1, 3, 2);
    plot(cS.kGridV, kPolM_eq(:, 1, plot_a_new), 'b-', 'LineWidth', 1, 'DisplayName', sprintf('e=%d (低)',1)); hold on;
    plot(cS.kGridV, kPolM_eq(:, round(cS.nw/2), plot_a_new), 'k--', 'LineWidth', 1, 'DisplayName', sprintf('e=%d (中)',round(cS.nw/2)));
    plot(cS.kGridV, kPolM_eq(:, cS.nw, plot_a_new), 'r:', 'LineWidth', 1, 'DisplayName', sprintf('e=%d (高)',cS.nw));
    plot(cS.kGridV, cS.kGridV, 'g--', 'DisplayName', 'k'' = k'); hold off;
    xlabel('期初非PPS资本 (k)'); ylabel('下一期非PPS资本 (k'')'); title(['非PPS储蓄政策 ' plot_title_suffix]); legend('Location', 'best'); grid on;
    subplot(1, 3, 3);
    if cS.pps_active && cS.pps_max_contrib_frac > 1e-9
        plot(cS.kGridV, cPpsPolM_eq(:, 1, plot_a_new), 'b-', 'LineWidth', 1, 'DisplayName', sprintf('e=%d (低)',1)); hold on;
        plot(cS.kGridV, cPpsPolM_eq(:, round(cS.nw/2), plot_a_new), 'k--', 'LineWidth', 1, 'DisplayName', sprintf('e=%d (中)',round(cS.nw/2)));
        plot(cS.kGridV, cPpsPolM_eq(:, cS.nw, plot_a_new), 'r:', 'LineWidth', 1, 'DisplayName', sprintf('e=%d (高)',cS.nw)); hold off;
        ylabel('PPS缴费 (c_{pps})'); title(['PPS缴费政策 ' plot_title_suffix]); legend('Location', 'best');
    else, title('PPS 未缴费或未激活'); axis off; end % 修改标题
     xlabel('期初非PPS资本 (k)'); grid on;
else, fprintf('跳过政策函数绘图 - 无效的 plot_a_new 索引或空的 physAgeMap。\n'); end

avgK_byAge_annual = mean(kHistM_eq, 1);
avgKpps_byAge_annual = zeros(1, cS.aD_orig);
if cS.pps_active && cS.pps_max_contrib_frac > 1e-9, avgKpps_byAge_annual = mean(kPpsHistM_eq, 1); end
avgTotalAssets_byAge_annual = avgK_byAge_annual + avgKpps_byAge_annual; % 如果PPS有资产才加

figure('Name', '生命周期资本剖面 (外生替代率 + PPS)');
plot(cS.physAgeV_orig, avgK_byAge_annual, 'bo-', 'LineWidth', 1.5, 'DisplayName', '平均非PPS资本 (k)'); hold on;
if cS.pps_active && cS.pps_max_contrib_frac > 1e-9
    plot(cS.physAgeV_orig, avgKpps_byAge_annual, 'ro-', 'LineWidth', 1.5, 'DisplayName', '平均PPS资本 (k_{pps})');
    plot(cS.physAgeV_orig, avgTotalAssets_byAge_annual, 'k--', 'LineWidth', 1.5, 'DisplayName', '平均总资产 (k+k_{pps})');
end
xline(cS.ageRetire_orig, 'g--', 'LineWidth', 1, 'DisplayName', '退休开始'); hold off;
title('按年度年龄划分的平均期初资本持有量'); xlabel('真实年龄'); ylabel('平均资本'); legend('Location', 'best'); grid on; xlim([cS.age1_orig, cS.ageLast_orig]);

fprintf('\n--- 分析完成 (外生替代率 + PPS) ---\n');