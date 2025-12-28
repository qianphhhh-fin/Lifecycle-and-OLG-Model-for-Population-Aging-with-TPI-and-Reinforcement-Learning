% --- START OF FILE main_olg_v3_pps.m ---

% Huggett (1996) 复制，包含动态人口结构
% 修改版：采用“混合时间单位”，包含个人养老金计划 (PPS),
% 并且PAYG工资税率内生，替代率外生
% V3: Incorporates Heer (2020) style government finance (tau_k, tau_c, tau_l, G, B)

clc; clear; close all;
fprintf('=== Huggett 模型: 动态人口, PPS, 内生 PAYG 税率, Heer GBC ===\n');

%% 1. 固定模型参数与人口设置
fprintf('\n--- 1. 初始化参数 ---\n');
cS = main_olg_v3_utils_pps.ParameterValues_HuggettStyle(); % cS: 包含所有常数和参数的结构体
paramS = struct(); % paramS: 包含派生参数的结构体

fprintf('参数已加载。nk=%d (非PPS资产网格点数), nw=%d (劳动效率状态数)。\n', cS.nk, cS.nw);
fprintf('年度年龄范围: %d-%d (%d 年)。年龄组数: %d (%d 个工作年龄组)。\n', ...
        cS.age1_orig, cS.ageLast_orig, cS.aD_orig, cS.aD_new, cS.aR_new);
fprintf('VFI 中使用年度 beta = %.4f。外生 PAYG 替代率: %.2f\n', cS.beta, cS.pension_replacement_rate);
fprintf('税率: tau_k=%.2f, tau_c=%.2f, tau_l=%.2f\n', cS.tau_k, cS.tau_c, cS.tau_l);
fprintf('政府财政: G/Y=%.2f, B/Y=%.2f\n', cS.gov_exp_frac_Y, cS.gov_debt_frac_Y);
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
popS = main_olg_v3_utils_pps.initPopulation(cS);
popS = main_olg_v3_utils_pps.populationDynamics(popS, cS);
[Z_ss, dep_ratio_ss, bgp_reached, bgp_period] = main_olg_v3_utils_pps.detectSteadyStatePopulation(popS, cS);

fprintf('\n--- 人口模拟总结 ---\n');
fprintf('实际模拟期数: %d\n', length(popS.totalPop)-1);
if bgp_reached, fprintf('人口稳态 (年龄组) 在第 %d 期达到。\n', bgp_period);
else, fprintf('人口稳态 (年龄组) 未达到。使用第 %d 期的最终数据。\n', bgp_period); end
fprintf('稳态抚养比 (退休人口数/工作人口数): %.4f\n', dep_ratio_ss);

paramS.Z_ss_counts = Z_ss; 

Z_ss_total = sum(Z_ss);
Z_ss_norm_group = zeros(cS.aD_new, 1);
if Z_ss_total > 1e-9, Z_ss_norm_group = Z_ss / Z_ss_total; end
paramS.ageMassV = Z_ss_norm_group(:); % Normalized steady-state age group distribution

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
    if sum(Z_ss_norm_annual) > 1e-9
        Z_ss_norm_annual = Z_ss_norm_annual / sum(Z_ss_norm_annual);
    end
    fprintf('已推导出近似的年度稳态人口分布。\n');
else
    warning('稳态人口为零，无法推导年度分布。');
    Z_ss_norm_annual(:) = 1/cS.aD_orig; 
end
paramS.Z_ss_norm_annual = Z_ss_norm_annual; % Store for potential use

% Calculate steady-state population growth rate from Z_ss for GBC consistency
if Z_ss_total > 1e-9 && length(popS.totalPop) > 1 && popS.totalPop(end-1) > 1e-9
    paramS.popGrowth_ss = popS.totalPop(end)/popS.totalPop(end-1) - 1;
    fprintf('稳态人口增长率 (来自模拟): %.4f (年度化近似: %.4f)\n', paramS.popGrowth_ss, (1+paramS.popGrowth_ss)^(1/cS.yearStep) -1 );
    % For GBC, we often use the annual BGP growth rate.
    % If popS.totalPop is group-level, then paramS.popGrowth_ss is per-group-period.
    % Let's stick to cS.popGrowth_orig for BGP debt dynamics for now, or refine this.
    paramS.popGrowthForDebt = cS.popGrowth_orig; % Annual BGP growth rate for debt
else
    warning('无法从人口模拟计算稳态增长率，使用 cS.popGrowth_orig = %.4f', cS.popGrowth_orig);
    paramS.popGrowthForDebt = cS.popGrowth_orig;
end


%% 3. 预计算劳动供给和禀赋过程
fprintf('\n--- 3. 预计算劳动 ---\n');
[paramS.leLogGridV, paramS.leTrProbM, paramS.leProb1V] = main_olg_v3_utils_pps.EarningProcess_olgm(cS);
paramS.leGridV = exp(paramS.leLogGridV(:));
fprintf('劳动禀赋过程已校准 (nw=%d 个状态)。\n', cS.nw);
paramS.ageEffV_orig = cS.ageEffV_orig; paramS.ageEffV_new = cS.ageEffV_new;
fprintf('为 %d 个个体模拟 %d 年的年度劳动禀赋...\n', cS.nSim, cS.aD_orig);
eIdxM = main_olg_v3_utils_pps.LaborEndowSimulation_olgm(cS, paramS);
[~, L] = main_olg_v3_utils_pps.LaborSupply_Huggett(eIdxM, cS, paramS, Z_ss_norm_group);
fprintf('总劳动供给 (L_ss, 效率单位): %.4f\n', L);
if L <= 0 && sum(paramS.Z_ss_counts(1:cS.aR_new)) > 1e-9 
    error('尽管存在正的工作人口，总劳动供给仍为零。');
elseif L <= 0, warning('总劳动供给为零或负。'); L = 1e-6; end

%% 4. 通过迭代求解一般均衡 (K*, TR_total*, theta_payg*)
fprintf('\n--- 4. 求解一般均衡 (外生替代率, PPS, Heer GBC) ---\n');

% 初始猜测值
% （1）设置为与V2 no pps相同的均衡
% KGuess = 9.2364   ; % 总生产性资本的猜测值 (Adjusted from previous runs)
% TRTotalGuess =  1.0266;   % 总净转移支付的猜测值 (政府净转移 + 意外遗赠) (Adjusted)

% (2) cS.pension_replacement_rate = 0.30; 
KGuess = 21.7207;
TRTotalGuess = 0.1347;

% (3) cS.pension_replacement_rate = 0.30;  cS.pps_max_contrib_frac = 0; 
KGuess = 11.6947;
TRTotalGuess =  0.2626;



% 迭代参数
maxIter = 200;
tolLevel = 1e-4;
dampK = 0.1;   % Adjusted damping
dampTR = 0.1;
iter = 0; devNorm = inf;

fprintf('开始均衡迭代...\n');
fprintf('Iter |   K Guess  | TRTotGuess | RnetHHfac | WnetHH  | ThetaPAYG | K Mod N-P | K Mod PPS | TRTot Mod |   K Dev    | TRTot Dev  |   Norm     | Time\n');
fprintf('--------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n');

eq_converged = false;
K_hist = zeros(maxIter+1,1); TR_hist = zeros(maxIter+1,1); Dev_hist = zeros(maxIter+1,1); 
R_k_net_factor_hist = zeros(maxIter+1,1); Theta_payg_hist = zeros(maxIter+1,1); 
K_hist(1) = KGuess; TR_hist(1) = TRTotalGuess; Dev_hist(1) = devNorm; R_k_net_factor_hist(1)=NaN; Theta_payg_hist(1)=NaN;
KModel_nonpps = NaN; KModel_pps = NaN; TRTotalModel = NaN; 

for iter = 1:maxIter
    iter_start_time = tic;

    % --- 第 4a 步: 计算税前价格, PAYG福利目标, 和内生PAYG税率 ---
    [R_mkt_gross_factor_iter, MPL_gross_iter] = main_olg_v3_utils_pps.HHPrices_Huggett(KGuess, L, cS);
    r_mkt_gross_iter = R_mkt_gross_factor_iter - 1; % This is MPK_gross - ddk

    b_payg_target_iter = cS.pension_replacement_rate * MPL_gross_iter;
    b_payg_target_iter = max(0, b_payg_target_iter);
    retiree_count_iter = sum(paramS.Z_ss_counts(cS.aR_new + 1 : cS.aD_new)); 
    total_pension_outlay_target_iter = b_payg_target_iter * retiree_count_iter;
    total_gross_wage_bill_iter = MPL_gross_iter * L;
    theta_payg_needed_iter = 0;
    if total_gross_wage_bill_iter > 1e-9
        theta_payg_needed_iter = total_pension_outlay_target_iter / total_gross_wage_bill_iter;
    else
        if total_pension_outlay_target_iter > 1e-9, warning('Iter %d: Gross wage bill zero but pension outlay positive!', iter); theta_payg_needed_iter = cS.max_payg_payroll_tax_rate; end
    end
    theta_payg_actual_iter = max(0, min(theta_payg_needed_iter, cS.max_payg_payroll_tax_rate));

    % Net prices for households
    r_k_net_hh_iter = r_mkt_gross_iter * (1 - cS.tau_k); % Net rate of return on capital for HH
    R_k_net_factor_hh_iter = 1 + r_k_net_hh_iter; 
    
    w_net_hh_iter = MPL_gross_iter * (1 - theta_payg_actual_iter - cS.tau_l); 
    w_net_hh_iter = max(0, w_net_hh_iter);

    bV_payg_for_vfi_iter = zeros(1, cS.aD_new);
    if cS.aR_new < cS.aD_new
        bV_payg_for_vfi_iter(cS.aR_new + 1 : cS.aD_new) = b_payg_target_iter;
    end

    % --- 第 4b 步: 求解家庭问题 (VFI) ---
    [cPolM, kPolM, cPpsPolM, ~] = main_olg_v3_utils_pps.HHSolution_VFI_Huggett(R_k_net_factor_hh_iter, w_net_hh_iter, TRTotalGuess, bV_payg_for_vfi_iter, paramS, cS);

    % --- 第 4c 步: 模拟家庭决策 ---
    [kHistM_non_pps, kPpsHistM, cHistM] = main_olg_v3_utils_pps.HHSimulation_olgm(kPolM, cPpsPolM, cPolM, eIdxM, R_k_net_factor_hh_iter, w_net_hh_iter, TRTotalGuess, bV_payg_for_vfi_iter, paramS, cS);

    % --- 第 4d 步: 计算模型的总生产性资本和总消费 ---
    KModel_nonpps = mean(kHistM_non_pps, 1) * Z_ss_norm_annual; % Aggregate non-PPS K
    KModel_nonpps = max(1e-6, KModel_nonpps);
    KModel_pps = 0;
    if cS.pps_active && cS.pps_in_K && cS.pps_max_contrib_frac > 1e-9
        KModel_pps = mean(kPpsHistM, 1) * Z_ss_norm_annual; % Aggregate PPS K
        KModel_pps = max(0, KModel_pps);
    end
    KModel = KModel_nonpps + KModel_pps;
    
    CModel_iter = mean(cHistM,1) * Z_ss_norm_annual; % Aggregate consumption quantity

    % --- 第 4e 步: 计算模型的年度总意外遗赠 (T_bequest_Model) ---
    kprimeHistM_for_bequest = zeros(cS.nSim, cS.aD_orig);
    if cS.aD_orig > 1
        kprimeHistM_for_bequest(:, 1:cS.aD_orig-1) = kHistM_non_pps(:, 2:cS.aD_orig); % k' for bequest
    end
    ageDeathMassV_annual = Z_ss_norm_annual(:) .* cS.d_orig(:);
    % Bequests are valued at the net return factor households would have received
    mean_bequest_value_iter = mean(kprimeHistM_for_bequest * R_k_net_factor_hh_iter, 1); 
    TotalBequests_iter_val = sum(mean_bequest_value_iter(:) .* ageDeathMassV_annual(:));
    T_bequest_Model = TotalBequests_iter_val / (1 + paramS.popGrowthForDebt); % Distributed to next gen growing at popGrowthForDebt
    T_bequest_Model = max(0, T_bequest_Model);

    % --- 第 4f 步: 政府预算平衡与模型隐含的总转移支付 (TRTotalModel) ---
    Y_iter = cS.A * (KGuess^cS.alpha) * (L^(1-cS.alpha)); 
    G_fixed_val_iter = cS.gov_exp_frac_Y * Y_iter;       
    B_fixed_val_iter = cS.gov_debt_frac_Y * Y_iter;       
    
    LaborTaxRevenue_iter = (theta_payg_actual_iter + cS.tau_l) * MPL_gross_iter * L;
    CapitalTaxRevenue_iter = (r_mkt_gross_iter) * KGuess * cS.tau_k; % Tax on (MPK_gross - ddk) * K_guess
    ConsumptionTaxRevenue_iter = CModel_iter * cS.tau_c; % Tax on consumption quantity
    TotalTaxRevenue_iter = LaborTaxRevenue_iter + CapitalTaxRevenue_iter + ConsumptionTaxRevenue_iter;

    % Interest rate on government debt. Assume it's the net rate households get on capital.
    r_b_iter = r_k_net_hh_iter; 
    
    % Government budget constraint to find TR_gov_Model:
    % G + TR_gov + PAYG_pensions + r_b*B = TotalTaxRev + n*B (where n*B is net new debt for SS)
    % TR_gov = TotalTaxRev - G - PAYG_pensions - (r_b - n)*B
    TR_gov_Model = TotalTaxRevenue_iter - G_fixed_val_iter - total_pension_outlay_target_iter ...
                   - (r_b_iter - paramS.popGrowthForDebt) * B_fixed_val_iter; 
    TR_gov_Model = max(0, TR_gov_Model); % Assume government transfers are non-negative

    TRTotalModel = TR_gov_Model + T_bequest_Model;

    % --- 第 4g 步: 计算偏差并检查收敛 ---
    KDev = KGuess - KModel;
    TRTotalDev = TRTotalGuess - TRTotalModel;
    devNorm = sqrt(KDev^2 + TRTotalDev^2);

    iter_time = toc(iter_start_time);
    fprintf('%4d | %10.4f | %10.4f | %9.4f | %9.4f | %9.4f | %10.4f | %10.4f | %10.4f | %10.4f | %10.4f | %10.4e | %.2fs\n', ...
             iter, KGuess, TRTotalGuess, R_k_net_factor_hh_iter, w_net_hh_iter, theta_payg_actual_iter, KModel_nonpps, KModel_pps, TRTotalModel, KDev, TRTotalDev, devNorm, iter_time);

     K_hist(iter+1) = KModel; TR_hist(iter+1) = TRTotalModel; Dev_hist(iter+1) = devNorm;
     R_k_net_factor_hist(iter+1) = R_k_net_factor_hh_iter; Theta_payg_hist(iter+1) = theta_payg_actual_iter;

    if devNorm < tolLevel
        fprintf('均衡已收敛!\n');
        eq_converged = true;
        KGuess = KModel; TRTotalGuess = TRTotalModel;
        break;
    end

    KGuess = KGuess - dampK * KDev;
    TRTotalGuess = TRTotalGuess - dampTR * TRTotalDev;
    KGuess = max(1e-6, KGuess); TRTotalGuess = max(0, TRTotalGuess); 
end

if ~eq_converged
    fprintf('\n警告: %d 次迭代后均衡未收敛。\n', maxIter);
    fprintf('最终偏差范数: %.4e\n', devNorm);
    % Use the model-implied values from the last iteration if not converged
    if iter == maxIter 
        KGuess = KModel; TRTotalGuess = TRTotalModel;
    end
end

K_eq = KGuess; TR_total_eq = TRTotalGuess; 
theta_payg_eq = Theta_payg_hist(iter+1); 
R_k_net_factor_hh_eq = R_k_net_factor_hist(iter+1); % This is the net factor for HH

% --- 重新计算最终均衡价格与组成部分 ---
[R_mkt_gross_factor_eq, MPL_gross_eq] = main_olg_v3_utils_pps.HHPrices_Huggett(K_eq, L, cS);
r_mkt_gross_eq = R_mkt_gross_factor_eq - 1; % MPK_gross - ddk

w_net_hh_eq = MPL_gross_eq * (1 - theta_payg_eq - cS.tau_l); w_net_hh_eq = max(0,w_net_hh_eq);
b_payg_eq = cS.pension_replacement_rate * MPL_gross_eq; b_payg_eq = max(0, b_payg_eq);
bV_eq_new = zeros(1, cS.aD_new);
if cS.aR_new < cS.aD_new, bV_eq_new(cS.aR_new + 1 : cS.aD_new) = b_payg_eq; end

% Re-solve VFI and simulate HH with equilibrium prices/transfers
[cPolM_eq, kPolM_eq, cPpsPolM_eq, valueM_eq] = main_olg_v3_utils_pps.HHSolution_VFI_Huggett(R_k_net_factor_hh_eq, w_net_hh_eq, TR_total_eq, bV_eq_new, paramS, cS);
fprintf('模拟均衡状态下的最终分布...\n')
[kHistM_eq, kPpsHistM_eq, cHistM_eq] = main_olg_v3_utils_pps.HHSimulation_olgm(kPolM_eq, cPpsPolM_eq, cPolM_eq, eIdxM, R_k_net_factor_hh_eq, w_net_hh_eq, TR_total_eq, bV_eq_new, paramS, cS);
fprintf('最终模拟完成。\n')

C_eq = mean(cHistM_eq,1) * Z_ss_norm_annual; % Equilibrium aggregate consumption quantity
Y_eq = cS.A * (K_eq^cS.alpha) * (L^(1-cS.alpha));
G_eq = cS.gov_exp_frac_Y * Y_eq;
B_eq = cS.gov_debt_frac_Y * Y_eq;

kprimeHistM_for_bequest_eq = zeros(cS.nSim, cS.aD_orig);
if cS.aD_orig > 1
    kprimeHistM_for_bequest_eq(:, 1:cS.aD_orig-1) = kHistM_eq(:, 2:cS.aD_orig);
end
mean_bequest_value_eq = mean(kprimeHistM_for_bequest_eq * R_k_net_factor_hh_eq, 1); % Valued at net factor
TotalBequests_eq_val = sum(mean_bequest_value_eq(:) .* (Z_ss_norm_annual(:) .* cS.d_orig(:)));
T_bequest_eq = TotalBequests_eq_val / (1 + paramS.popGrowthForDebt);
T_bequest_eq = max(0, T_bequest_eq);

TR_gov_eq = TR_total_eq - T_bequest_eq; 

% GBC check at equilibrium
LaborTax_eq = (theta_payg_eq + cS.tau_l) * MPL_gross_eq * L;
CapitalTax_eq = r_mkt_gross_eq * K_eq * cS.tau_k;
ConsumptionTax_eq = C_eq * cS.tau_c;
TotalTax_eq = LaborTax_eq + CapitalTax_eq + ConsumptionTax_eq;
PAYG_outlay_eq = b_payg_eq * sum(paramS.Z_ss_counts(cS.aR_new + 1 : cS.aD_new));
r_b_eq = R_k_net_factor_hh_eq - 1; % Interest rate on debt assumed to be HH net capital return rate

GovBudgetResidual = TotalTax_eq - G_eq - TR_gov_eq - PAYG_outlay_eq - (r_b_eq - paramS.popGrowthForDebt)*B_eq;
fprintf('GBC Check at Equilibrium: Residual = %.4e (should be close to zero)\n', GovBudgetResidual);
if abs(GovBudgetResidual) > 1e-3, warning('GBC residual is large at equilibrium!'); end


%% 5. 分析和绘制结果
fprintf('\n--- 5. 均衡结果与绘图 (外生替代率 + PPS + Heer GBC) ---\n');
K_nonpps_eq = mean(kHistM_eq, 1) * Z_ss_norm_annual;
K_pps_eq = 0; if cS.pps_active && cS.pps_max_contrib_frac > 1e-9, K_pps_eq = mean(kPpsHistM_eq, 1) * Z_ss_norm_annual; end
TotalAssets_eq = K_nonpps_eq + K_pps_eq;

fprintf('均衡总生产性资本 (K*): %.4f\n', K_eq);
fprintf('  分解: 非PPS资本 K: %.4f, PPS资本 K: %.4f\n', K_nonpps_eq, K_pps_eq);
fprintf('均衡家庭总资产 (非PPS+PPS): %.4f\n', TotalAssets_eq);
fprintf('均衡总劳动 (L):    %.4f\n', L);
fprintf('均衡总产出 (Y* 年度):   %.4f\n', Y_eq);
fprintf('均衡市场利率因子 (R_mkt_gross_factor*): %.4f (r_mkt_gross*=%.4f)\n', R_mkt_gross_factor_eq, r_mkt_gross_eq);
fprintf('  家庭净回报率因子 (R_k_net_factor_hh*): %.4f (r_k_net_hh*=%.4f)\n', R_k_net_factor_hh_eq, R_k_net_factor_hh_eq-1);
fprintf('均衡总工资率 (MPL_gross*): %.4f\n', MPL_gross_eq);
fprintf('均衡内生PAYG工资税率 (theta_payg*): %.4f\n', theta_payg_eq);
fprintf('  其他劳动税率 (tau_l fixed): %.4f\n', cS.tau_l);
fprintf('  资本所得税率 (tau_k fixed): %.4f\n', cS.tau_k);
fprintf('  消费税率 (tau_c fixed): %.4f\n', cS.tau_c);
fprintf('均衡家庭净工资率 (w_net_hh*): %.4f\n', w_net_hh_eq);
fprintf('均衡PAYG福利 (b_payg*): %.4f\n', b_payg_eq);
fprintf('均衡总净转移支付 (TR_total*): %.4f\n', TR_total_eq);
fprintf('  其中意外遗赠 (T_bequest*): %.4f\n', T_bequest_eq);
fprintf('  其中政府净转移 (TR_gov*): %.4f\n', TR_gov_eq);
fprintf('均衡政府消费 (G*): %.4f (G/Y = %.2f)\n', G_eq, G_eq/Y_eq);
fprintf('均衡政府债务 (B*): %.4f (B/Y = %.2f)\n', B_eq, B_eq/Y_eq);
fprintf('均衡 K/Y 比率 (年度): %.4f\n', K_eq / Y_eq);
fprintf('均衡 C/Y 比率 (年度): %.4f\n', C_eq / Y_eq);

% --- 绘图 ---
figure('Name', '均衡收敛过程 (Heer GBC)');
subplot(4,1,1); plot(0:iter, K_hist(1:iter+1), '-o'); title('K 收敛'); ylabel('K'); grid on; xlim([0 iter]);
subplot(4,1,2); plot(0:iter, TR_hist(1:iter+1), '-o'); title('TR_{total} 收敛'); ylabel('TR_{total}'); grid on; xlim([0 iter]);
subplot(4,1,3); plot(0:iter, Theta_payg_hist(1:iter+1), '-o'); title('Theta PAYG 收敛'); ylabel('\theta_{PAYG}'); grid on; xlim([0 iter]);
subplot(4,1,4); semilogy(1:iter, Dev_hist(2:iter+1), '-o'); title('偏差范数收敛'); ylabel('Norm'); xlabel('迭代次数'); grid on; xlim([1 iter]);
sgtitle('均衡迭代收敛过程 (Heer GBC)');

plot_a_new = min(4, cS.aR_new);
if plot_a_new <= cS.aD_new && plot_a_new > 0 && ~isempty(cS.physAgeMap{plot_a_new}) && cS.nk > 1 % Added cS.nk condition
    physAgeStart = cS.physAgeV_new(plot_a_new);
    physAgeEnd = cS.physAgeV_orig(cS.physAgeMap{plot_a_new}(end));
    plot_title_suffix = sprintf('年龄组 %d (年龄 %d-%d)', plot_a_new, physAgeStart, physAgeEnd);
    figure('Name', ['政策函数 (Heer GBC): ' plot_title_suffix]);
    subplot(1, 3, 1);
    plot(cS.kGridV, cPolM_eq(:, 1, plot_a_new), 'b-', 'LineWidth', 1, 'DisplayName', sprintf('e=%d (低)',1)); hold on;
    plot(cS.kGridV, cPolM_eq(:, round(cS.nw/2), plot_a_new), 'k--', 'LineWidth', 1, 'DisplayName', sprintf('e=%d (中)',round(cS.nw/2)));
    plot(cS.kGridV, cPolM_eq(:, cS.nw, plot_a_new), 'r:', 'LineWidth', 1, 'DisplayName', sprintf('e=%d (高)',cS.nw)); hold off;
    xlabel('期初非PPS资本 (k)'); ylabel('最优消费 (c_{quantity})'); title(['消费政策 ' plot_title_suffix]); legend('Location', 'best'); grid on;
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
    else, title('PPS 未缴费或未激活'); axis off; end 
     xlabel('期初非PPS资本 (k)'); grid on;
else, fprintf('跳过政策函数绘图 - plot_a_new 无效, physAgeMap 为空, 或 nk<=1。\n'); end

avgK_byAge_annual = mean(kHistM_eq, 1);
avgKpps_byAge_annual = zeros(1, cS.aD_orig);
if cS.pps_active && cS.pps_max_contrib_frac > 1e-9, avgKpps_byAge_annual = mean(kPpsHistM_eq, 1); end
avgTotalAssets_byAge_annual = avgK_byAge_annual + avgKpps_byAge_annual; 

figure('Name', '生命周期资本剖面 (Heer GBC)');
plot(cS.physAgeV_orig, avgK_byAge_annual, 'bo-', 'LineWidth', 1.5, 'DisplayName', '平均非PPS资本 (k)'); hold on;
if cS.pps_active && cS.pps_max_contrib_frac > 1e-9
    plot(cS.physAgeV_orig, avgKpps_byAge_annual, 'ro-', 'LineWidth', 1.5, 'DisplayName', '平均PPS资本 (k_{pps})');
    plot(cS.physAgeV_orig, avgTotalAssets_byAge_annual, 'k--', 'LineWidth', 1.5, 'DisplayName', '平均总资产 (k+k_{pps})');
end
xline(cS.ageRetire_orig, 'g--', 'LineWidth', 1, 'DisplayName', '退休开始'); hold off;
title('按年度年龄划分的平均期初资本持有量'); xlabel('真实年龄'); ylabel('平均资本'); legend('Location', 'best'); grid on; xlim([cS.age1_orig, cS.ageLast_orig]);

fprintf('\n--- 分析完成 (外生替代率 + PPS + Heer GBC) ---\n');

% --- END OF FILE main_olg_v3_pps.m ---