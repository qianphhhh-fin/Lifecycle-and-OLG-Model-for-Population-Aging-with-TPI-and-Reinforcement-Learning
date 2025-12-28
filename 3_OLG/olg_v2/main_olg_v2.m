% --- START OF FILE main_olg_v2.m ---

% OLG 模型: 简化模型PAYG (外生theta, 内生b) + PPS决策
% - 动态人口结构
% - "混合时间单位" VFI 和模拟

clc; clear; close all;
fprintf('=== OLG Model: Simplified PAYG (Exo Theta) + PPS ===\n');

%% 1. 固定模型参数与人口设置
fprintf('\n--- 1. 初始化参数 ---\n');
cS = main_olg_v2_utils.ParameterValues_HuggettStyle(); % 使用新的utils文件名 (main_olg_v2_utils.m)
paramS = struct(); % 初始化空的paramS，用于存储派生参数

fprintf('参数已加载。nk=%d (非PPS资产网格点数), nw=%d (劳动效率状态数)。\n', cS.nk, cS.nw);
fprintf('年度年龄范围: %d-%d (%d 年)。年龄组数: %d (%d 个工作年龄组)。\n', ...
        cS.age1_orig, cS.ageLast_orig, cS.aD_orig, cS.aD_new, cS.aR_new);
fprintf('VFI 中使用年度 beta = %.4f。外生 PAYG 贡献率 (theta) = %.2f\n', cS.beta, cS.theta);
if cS.pps_active
    fprintf('个人养老金计划 (PPS) 已激活。\n');
    fprintf('  最大缴费比例: %.2f, 领取期税率: %.2f, 回报率溢价: %.3f, 领取率: %.2f\n', ...
        cS.pps_max_contrib_frac, cS.pps_tax_rate_withdrawal, cS.pps_return_rate_premium, cS.pps_withdrawal_rate);
    fprintf('  PPS 计入总资本 K: %d, PPS 可遗赠: %d\n', cS.pps_in_K, cS.pps_bequeathable);
else
    fprintf('个人养老金计划 (PPS) 未激活。\n');
end

%% 2. 模拟人口动态至稳态 (年龄组层面)
fprintf('\n--- 2. 模拟人口动态 (年龄组层面) ---\n');
popS = main_olg_v2_utils.initPopulation(cS);
popS = main_olg_v2_utils.populationDynamics(popS, cS);
[Z_ss, dep_ratio_ss, bgp_reached, bgp_period] = main_olg_v2_utils.detectSteadyStatePopulation(popS, cS); % Z_ss 是稳态人口数量

fprintf('\n--- 人口模拟总结 ---\n');
fprintf('实际模拟期数: %d\n', length(popS.totalPop)-1);
if bgp_reached, fprintf('人口稳态 (年龄组) 在第 %d 期达到。\n', bgp_period);
else, fprintf('人口稳态 (年龄组) 未达到。使用第 %d 期的最终数据。\n', bgp_period); end
fprintf('稳态抚养比 (退休人口数/工作人口数): %.4f\n', dep_ratio_ss);

% 将稳态人口数量 Z_ss 存入 paramS 以便传递给价格函数
paramS.Z_ss_counts = Z_ss; % 用于在HHPrices_Huggett中计算b

% 归一化稳态年龄组人口分布
Z_ss_total = sum(Z_ss);
Z_ss_norm_group = zeros(cS.aD_new,1);
if Z_ss_total > 1e-9, Z_ss_norm_group = Z_ss / Z_ss_total; end
paramS.ageMassV = Z_ss_norm_group(:); % 用于计算总劳动L (L的计算用占比)

% 从年龄组分布推导近似的年度稳态分布 (用于K和T的加总)
Z_ss_norm_annual = zeros(cS.aD_orig, 1);
if Z_ss_total > 1e-9
    for a_new = 1:cS.aD_new
        annual_indices = cS.physAgeMap{a_new};
        group_mass = Z_ss_norm_group(a_new);
        num_years_in_group = length(annual_indices);
        if num_years_in_group > 0
            mass_per_year = group_mass / num_years_in_group;
            Z_ss_norm_annual(annual_indices) = mass_per_year;
        end
    end
    if sum(Z_ss_norm_annual) > 1e-9, Z_ss_norm_annual = Z_ss_norm_annual / sum(Z_ss_norm_annual); end
    fprintf('已推导出近似的年度稳态人口分布。\n');
else, warning('稳态人口为零，无法推导年度分布。'); end

%% 3. 预计算劳动供给和禀赋过程
fprintf('\n--- 3. 预计算劳动 ---\n');
[paramS.leLogGridV, paramS.leTrProbM, paramS.leProb1V] = main_olg_v2_utils.EarningProcess_olgm(cS);
paramS.leGridV = exp(paramS.leLogGridV(:));
fprintf('劳动禀赋过程已校准 (nw=%d 个状态)。\n', cS.nw);
paramS.ageEffV_orig = cS.ageEffV_orig; paramS.ageEffV_new = cS.ageEffV_new;
fprintf('为 %d 个个体模拟 %d 年的年度劳动禀赋...\n', cS.nSim, cS.aD_orig);
eIdxM = main_olg_v2_utils.LaborEndowSimulation_olgm(cS, paramS);
[~, L] = main_olg_v2_utils.LaborSupply_Huggett(eIdxM, cS, paramS, Z_ss_norm_group);
fprintf('总劳动供给 (L_ss, 效率单位): %.4f\n', L);
if L <= 0 && sum(Z_ss(1:cS.aR_new)) > 1e-9
    error('尽管存在正的工作人口，总劳动供给仍为零。');
elseif L <= 0, warning('总劳动供给为零或负。'); L = 1e-6; end

%% 4. 通过迭代求解一般均衡 (K*, T*)
fprintf('\n--- 4. 求解一般均衡 (简化PAYG + PPS) ---\n');

% 初始猜测值 (可以基于之前的收敛结果调整)
KGuess =  52.8894;  % 总生产性资本的猜测值
TGuess =  0.4030;   % 年度 lump-sum 转移支付的猜测值

% 迭代参数
maxIter = 500; % 增加迭代次数
tolLevel = 1e-5;
dampK = 0.1;   % 尝试较小的阻尼因子
dampT = 0.1;
iter = 0; devNorm = inf;

fprintf('开始均衡迭代...\n');
fprintf('Iter |   K Guess  |  T Guess   | R_market  | PAYG_b | K Mod N-P | K Mod PPS  |  T Model   |   K Dev    |   T Dev    |   Norm     | Time\n');
fprintf('---------------------------------------------------------------------------------------------------------------------------------------------------\n');

eq_converged = false;
K_hist = zeros(maxIter+1,1); T_hist = zeros(maxIter+1,1); Dev_hist = zeros(maxIter+1,1);
R_hist = zeros(maxIter+1,1); B_payg_hist = zeros(maxIter+1,1); % 用于记录R和b
K_hist(1) = KGuess; T_hist(1) = TGuess; Dev_hist(1) = devNorm; R_hist(1)=NaN; B_payg_hist(1)=NaN;
KModel_nonpps = NaN; KModel_pps = NaN; TModel = NaN; % 初始化用于打印

for iter = 1:maxIter
    iter_start_time = tic;

    % --- 第 4a 步: 计算价格 (R, w, b) ---
    % HHPrices_Huggett 现在使用 paramS.Z_ss_counts 计算 b
    [~, R_iter, w_iter, b_iter_payg, ~] = main_olg_v2_utils.HHPrices_Huggett(KGuess, L, cS, paramS);
    % R_iter: 市场回报率因子 (简化模型中假设无资本税)
    % w_iter: 净工资 (已扣除固定的cS.theta)
    % b_iter_payg: 内生计算的人均PAYG福利

    % 为VFI创建PAYG福利向量
    bV_payg_for_vfi = zeros(1, cS.aD_new);
    if cS.aR_new < cS.aD_new
        bV_payg_for_vfi(cS.aR_new + 1 : cS.aD_new) = b_iter_payg;
    end

    % --- 第 4b 步: 求解家庭问题 (VFI，包含PPS选择逻辑) ---
    % VFI 的状态变量是 k_non_pps
    [cPolM, kPolM, cPpsPolM, ~] = main_olg_v2_utils.HHSolution_VFI_Huggett(R_iter, w_iter, TGuess, bV_payg_for_vfi, paramS, cS);

    % --- 第 4c 步: 模拟家庭决策 (年度，追踪 k_non_pps 和 k_pps) ---
    [kHistM_non_pps, kPpsHistM, cHistM] = main_olg_v2_utils.HHSimulation_olgm(kPolM, cPpsPolM, cPolM, eIdxM, R_iter, w_iter, TGuess, bV_payg_for_vfi, paramS, cS);

    % --- 第 4d 步: 计算模型的总生产性资本 ---
    KModel_nonpps = mean(kHistM_non_pps, 1) * Z_ss_norm_annual;
    KModel_nonpps = max(1e-6, KModel_nonpps);
    KModel_pps = 0;
    if cS.pps_active && cS.pps_in_K
        KModel_pps = mean(kPpsHistM, 1) * Z_ss_norm_annual;
        KModel_pps = max(0, KModel_pps);
    end
    KModel = KModel_nonpps + KModel_pps; % 总生产性资本

    % --- 第 4e 步: 计算模型的年度总意外遗赠 (TModel) ---
    % 假设只有非PPS资产 kHistM_non_pps 可以遗赠 (如果 cS.pps_bequeathable = false)
    kprimeHistM_for_bequest = zeros(cS.nSim, cS.aD_orig);
    if cS.aD_orig > 1
        kprimeHistM_for_bequest(:, 1:cS.aD_orig-1) = kHistM_non_pps(:, 2:cS.aD_orig);
    end
    ageDeathMassV_annual = Z_ss_norm_annual(:) .* cS.d_orig(:);
    mean_bequest_value = mean(kprimeHistM_for_bequest * R_iter, 1); % 遗赠资产按市场回报率估值
    TotalBequests = sum(mean_bequest_value(:) .* ageDeathMassV_annual(:));
    TModel = TotalBequests / (1 + cS.popGrowth_orig);
    TModel = max(0, TModel);

    % --- 第 4f 步: 计算偏差并检查收敛 ---
    KDev = KGuess - KModel;
    TDev = TGuess - TModel;
    devNorm = sqrt(KDev^2 + TDev^2);

    iter_time = toc(iter_start_time);
    fprintf('%4d | %10.4f | %10.4f | %9.4f | %6.4f | %10.4f | %10.4f | %10.4f | %10.4f | %10.4f | %10.4e | %.2fs\n', ...
             iter, KGuess, TGuess, R_iter, b_iter_payg, KModel_nonpps, KModel_pps, TModel, KDev, TDev, devNorm, iter_time);

     K_hist(iter+1) = KModel; T_hist(iter+1) = TModel; Dev_hist(iter+1) = devNorm;
     R_hist(iter+1) = R_iter; B_payg_hist(iter+1) = b_iter_payg;

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
    KGuess = KModel; TGuess = TModel;
end

K_eq = KGuess; T_eq = TGuess;

% --- 重新计算最终均衡价格 ---
[Y_eq, R_eq, w_eq, b_eq_final, MPL_eq_final] = main_olg_v2_utils.HHPrices_Huggett(K_eq, L, cS, paramS);

bV_eq_new = zeros(1, cS.aD_new);
if cS.aR_new < cS.aD_new, bV_eq_new(cS.aR_new + 1 : cS.aD_new) = b_eq_final; end
[cPolM_eq, kPolM_eq, cPpsPolM_eq, valueM_eq] = main_olg_v2_utils.HHSolution_VFI_Huggett(R_eq, w_eq, T_eq, bV_eq_new, paramS, cS);

fprintf('模拟均衡状态下的最终分布...\n')
[kHistM_eq, kPpsHistM_eq, cHistM_eq] = main_olg_v2_utils.HHSimulation_olgm(kPolM_eq, cPpsPolM_eq, cPolM_eq, eIdxM, R_eq, w_eq, T_eq, bV_eq_new, paramS, cS);
fprintf('最终模拟完成。\n')

%% 5. 分析和绘制结果
fprintf('\n--- 5. 均衡结果与绘图 (简化PAYG + PPS) ---\n');
K_nonpps_eq = mean(kHistM_eq, 1) * Z_ss_norm_annual;
K_pps_eq = 0; if cS.pps_active, K_pps_eq = mean(kPpsHistM_eq, 1) * Z_ss_norm_annual; end
TotalAssets_eq = K_nonpps_eq; if cS.pps_active, TotalAssets_eq = K_nonpps_eq + K_pps_eq; end % 家庭总资产

fprintf('均衡总生产性资本 (K*): %.4f\n', K_eq);
fprintf('  分解: 非PPS资本 K: %.4f, PPS资本 K: %.4f\n', K_nonpps_eq, K_pps_eq);
fprintf('均衡家庭总资产 (非PPS+PPS): %.4f\n', TotalAssets_eq);
fprintf('均衡总劳动 (L):    %.4f\n', L);
fprintf('均衡总产出 (Y* 年度):   %.4f\n', Y_eq);
fprintf('均衡市场利率因子 (R* 年度): %.4f (利率 r*=%.4f)\n', R_eq, R_eq - 1);
fprintf('均衡总工资率 (MPL_gross*): %.4f\n', MPL_eq_final);
fprintf('均衡净工资率 (w* PAYG税后): %.4f\n', w_eq);
fprintf('均衡PAYG福利 (b* 内生): %.4f\n', b_eq_final);
fprintf('均衡意外遗赠转移支付 (T* 年度): %.4f\n', T_eq);
fprintf('均衡 K/Y 比率 (年度): %.4f\n', K_eq / Y_eq);

% 计算并打印有效平均替代率
total_retirees_from_Z_ss_final = sum(paramS.Z_ss_counts(cS.aR_new + 1 : cS.aD_new));
effective_replacement_rate_final = 0;
if total_retirees_from_Z_ss_final > 1e-9 && L > 1e-9
    effective_replacement_rate_final = (cS.theta * (L / total_retirees_from_Z_ss_final)) * 100;
end
fprintf('有效平均PAYG替代率 (theta * L_eff / Retiree_count): %.2f%%\n', effective_replacement_rate_final);
fprintf('  (基于: cS.theta=%.2f, L=%.4f, RetireeCount=%.2f, MPL_eq=%.4f, b_eq=%.4f)\n', ...
    cS.theta, L, total_retirees_from_Z_ss_final, MPL_eq_final, b_eq_final);

% --- 绘图 ---
figure('Name', '均衡收敛过程 (简化PAYG + PPS)');
subplot(3,1,1); plot(0:iter, K_hist(1:iter+1), '-o'); title('K 收敛'); ylabel('K'); grid on; xlim([0 iter]);
subplot(3,1,2); plot(0:iter, T_hist(1:iter+1), '-o'); title('T 收敛'); ylabel('T'); grid on; xlim([0 iter]);
subplot(3,1,3); semilogy(1:iter, Dev_hist(2:iter+1), '-o'); title('偏差范数收敛'); ylabel('Norm'); xlabel('迭代次数'); grid on; xlim([1 iter]);
sgtitle('均衡迭代收敛过程 (简化PAYG + PPS)');

plot_a_new = min(4, cS.aR_new);
if plot_a_new <= cS.aD_new && plot_a_new > 0 && ~isempty(cS.physAgeMap{plot_a_new})
    physAgeStart = cS.physAgeV_new(plot_a_new);
    physAgeEnd = cS.physAgeV_orig(cS.physAgeMap{plot_a_new}(end));
    plot_title_suffix = sprintf('年龄组 %d (年龄 %d-%d)', plot_a_new, physAgeStart, physAgeEnd);
    figure('Name', ['政策函数 (简化PAYG + PPS): ' plot_title_suffix]);
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
    if cS.pps_active
        plot(cS.kGridV, cPpsPolM_eq(:, 1, plot_a_new), 'b-', 'LineWidth', 1, 'DisplayName', sprintf('e=%d (低)',1)); hold on;
        plot(cS.kGridV, cPpsPolM_eq(:, round(cS.nw/2), plot_a_new), 'k--', 'LineWidth', 1, 'DisplayName', sprintf('e=%d (中)',round(cS.nw/2)));
        plot(cS.kGridV, cPpsPolM_eq(:, cS.nw, plot_a_new), 'r:', 'LineWidth', 1, 'DisplayName', sprintf('e=%d (高)',cS.nw)); hold off;
        ylabel('PPS缴费 (c_{pps})'); title(['PPS缴费政策 ' plot_title_suffix]); legend('Location', 'best');
    else, title('PPS 未激活'); axis off; end
     xlabel('期初非PPS资本 (k)'); grid on;
else, fprintf('跳过政策函数绘图 - 无效的 plot_a_new 索引或空的 physAgeMap。\n'); end

avgK_byAge_annual = mean(kHistM_eq, 1);
avgKpps_byAge_annual = zeros(1, cS.aD_orig);
if cS.pps_active, avgKpps_byAge_annual = mean(kPpsHistM_eq, 1); end
avgTotalAssets_byAge_annual = avgK_byAge_annual;
if cS.pps_active, avgTotalAssets_byAge_annual = avgK_byAge_annual + avgKpps_byAge_annual; end

figure('Name', '生命周期资本剖面 (简化PAYG + PPS)');
plot(cS.physAgeV_orig, avgK_byAge_annual, 'bo-', 'LineWidth', 1.5, 'DisplayName', '平均非PPS资本 (k)'); hold on;
if cS.pps_active
    plot(cS.physAgeV_orig, avgKpps_byAge_annual, 'ro-', 'LineWidth', 1.5, 'DisplayName', '平均PPS资本 (k_{pps})');
    plot(cS.physAgeV_orig, avgTotalAssets_byAge_annual, 'k--', 'LineWidth', 1.5, 'DisplayName', '平均总资产 (k+k_{pps})');
end
xline(cS.ageRetire_orig, 'g--', 'LineWidth', 1, 'DisplayName', '退休开始'); hold off;
title('按年度年龄划分的平均期初资本持有量'); xlabel('真实年龄'); ylabel('平均资本'); legend('Location', 'best'); grid on; xlim([cS.age1_orig, cS.ageLast_orig]);

fprintf('\n--- 分析完成 (简化PAYG + PPS) ---\n');