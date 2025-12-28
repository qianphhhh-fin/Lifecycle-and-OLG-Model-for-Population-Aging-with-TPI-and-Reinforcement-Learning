% --- START OF FILE main_olg_v14.m ---
%
% OLG 模型 v14: 基于真实数据校准和外生路径的过渡动态分析
%
% 核心变化:
% 1. 不再求解初始稳态，而是从start_year年的校准数据点出发。
% 2. 引入了外生的、随时间变化的TFP路径(A_path)和人口路径(Z_path)。
% 3. 社保基金(B_p)和个人养老金(K_pps)都从零开始，内生演化。
% 4. PPS系统在模拟的特定年份(2022年)被激活。
%
% ---

clear; close all;
addpath(pwd);
fprintf('=== OLG 模型 v14: 基于数据校准的过渡动态分析 ===\n');

%% 1. 初始化参数和环境
fprintf('\n--- 1. 初始化参数 ---\n');
cS = main_olg_v14_utils.ParameterValues_HuggettStyle();
paramS = struct();

% --- 家庭偏好参数 ---
cS.sigma = 3.0;
cS.beta = 0.97; % 5年期
cS.cFloor = 0.05;
cS.nSim = 1000;
cS.phi_bequest = 3; 
cS.sigma_bequest = cS.sigma;

cS.start_year = 1997;
cS.end_year = 2102;
cS.time_step = 5;

% --- 生产技术参数 (基于start_year年校准) ---
cS.alpha = 0.5; % [强烈建议] 大幅提高资本份额以解决负利率问题
cS.ddk = 1 - (1 - 0.05)^cS.time_step; % 5年期折旧率 (年化8%)

% --- 政府财政参数 (基于start_year年校准) ---
cS.tau_k = 0.25;
cS.tau_l = 0.08;
cS.tau_c = 0.04;
cS.gov_exp_frac_Y = 0.15;

% --- 网格参数 ---
ngrid = 20;
cS.nk = ngrid;
cS.nkpps = ngrid;
cS.nkprime = ngrid;
cS.npps = ngrid;
cS.A = 1.0; % [v14修正] 为网格生成提供一个基准A值
cS = main_olg_v14_utils.generateGrids(cS);
    
% --- 养老金和PPS制度参数 ---
cS.rho_prime_payg = 0.5; 
% cS.theta_payg = 0.05; 
cS.pps_active = false; % 初始为关闭状态
cS.pps_activation_year = 2022; % PPS激活年份
cS.pps_tax_rate_withdrawal = 0.03;
cS.pps_return_rate_premium = 0.01;
cS.pps_withdrawal_rate = 0.15;
cS.pps_contrib_limit = 9999;
cS.pps_max_contrib_frac = 0.1;

% [新的调用方式]
fprintf('正在构建基于现实制度差异的有效养老金缴费率路径 (theta_path)...\n');
% 将时间轴参数和绘图开关传递给新函数
cS = main_olg_v14_utils.calcaulte_theta_payg_path(cS, true);
fprintf('✅ 复合有效缴费率路径构建完成。\n');

%% 2. [v14核心] 加载外生数据路径
[Z_path, A_path, T] = main_olg_v14_utils.load_exogenous_paths(cS);
sim_years = cS.start_year:cS.time_Step:(cS.start_year + (T-1)*cS.time_Step);
cS.sim_years = sim_years; % 将年份向量存入cS


% --- 生成劳动效率冲击路径 (与之前相同) ---
[paramS.leLogGridV, paramS.leTrProbM, paramS.leProb1V] = main_olg_v14_utils.EarningProcess_olgm(cS);
paramS.leGridV = exp(paramS.leLogGridV(:));
paramS.ageEffV_new = cS.ageEffV_new;
eIdxM = main_olg_v14_utils.LaborEndowSimulation_olgm_AgeGroup(cS, paramS);
fprintf('✅ 效率冲击路径生成完成。\n');

% --- 在 main_olg_v14.m 的第2节 "初始化模型参数" 中 ---



% =========================================================================
% == 在 main_olg_v14.m 中，替换整个第3节 ==
% =========================================================================
% =========================================================================
% == 替换 main 脚本中的第3节 ==
% =========================================================================
%% 3. [v14 最终修正版] 基于数据校准start_year年的初始状态 (t=1)
fprintf('\n--- 3. 校准start_year年 (t=1) 的初始状态 ---\n');

% 初始化路径
paths = struct();
paths.K_pvt = zeros(1, T); paths.K_pps = zeros(1, T); paths.B_p = zeros(1, T);
paths.B_g = zeros(1, T); paths.w = zeros(1, T); paths.r_net = zeros(1, T);
paths.L = zeros(1, T); paths.Y = zeros(1, T); paths.G = zeros(1, T);
paths.C = zeros(1, T); paths.Investment = zeros(1, T-1);
paths.market_clearing_residual = zeros(1, T-1);
paths.Bp_Y_ratio = zeros(1, T);
paths.K_physical = zeros(1, T);

% 校准参数
KY_ratio_start_year = 1.994;
Bg_Y_ratio_start_year = 0.08;
A_1 = A_path(1);
Z_1 = Z_path(:, 1);

% [核心修正] 使用迭代法寻找自洽的初始宏观状态
fprintf('   通过迭代寻找自洽的 Y_annual(1) 和 K_physical(1)...\n');
Y_annual_guess = 0.3; % 从一个合理的年化产出猜测开始
tol = 1e-7; max_iter = 100; damp = 0.7; 

for iter = 1:max_iter
    K_physical_target = KY_ratio_start_year * Y_annual_guess;
    B_g_target = Bg_Y_ratio_start_year * Y_annual_guess;
    K_pvt_1 = K_physical_target + B_g_target;
    
    [~, L_1_main] = main_olg_v14_utils.LaborSupply_Huggett(eIdxM, cS, paramS, Z_1);
    [~, ~, Y_annual_new] = main_olg_v14_utils.HHPrices_Huggett(K_physical_target, L_1_main, A_1, cS);
    
    if abs(Y_annual_new - Y_annual_guess) < tol
        fprintf('   ✅ 自洽的初始状态已找到 (Y_annual=%.4f, K_phys=%.4f)。\n', Y_annual_new, K_physical_target);
        break;
    end
    Y_annual_guess = (1 - damp) * Y_annual_guess + damp * Y_annual_new;
end
if iter == max_iter, warning('初始状态迭代在%d次内未收敛。', max_iter); end

% 将收敛后的值赋给路径
Y_annual_1 = Y_annual_new;
paths.Y(1) = Y_annual_1 * cS.time_Step;
paths.K_physical(1) = KY_ratio_start_year * Y_annual_1;
paths.B_g(1) = Bg_Y_ratio_start_year * Y_annual_1;
paths.B_p(1) = 0;
paths.K_pps(1) = 0;
paths.K_pvt(1) = paths.K_physical(1) + paths.B_g(1);

% 获取最终的价格和宏观变量
M_1 = main_olg_v14_utils.get_prices_and_policy_at_t(1, paths.K_pvt(1), paths.B_p(1), paths.K_pps(1), paths.B_g(1), Z_1, A_1, cS, paramS, eIdxM);
paths.L(1) = M_1.L_t;
paths.w(1) = M_1.w_t;
paths.r_net(1) = M_1.r_net_annual;
paths.G(1) = cS.gov_exp_frac_Y * paths.Y(1);
paths.Bp_Y_ratio(1) = paths.B_p(1) / paths.Y(1);

fprintf('✅ start_year年初始状态校准完成: Y_total_5y=%.3f, K_physical=%.3f, K_pvt=%.3f, B_g=%.3f\n', ...
    paths.Y(1), paths.K_physical(1), paths.K_pvt(1), paths.B_g(1));
% =========================================================================
% == 替换 main 脚本中的第4节和第5节 ==
% =========================================================================
%% 4. [v14 最终修正版] 过渡动态主循环 (含正确的遗赠机制)
fprintf('\n--- 4. 开始向前模拟过渡动态 (start_year - %d) ---\n', sim_years(end));

% 打印表头
fprintf('\n%6s | %8s | %12s | %12s | %12s | %10s | %12s | %12s | %12s\n', ...
    'Time', 'Year', 'Physical K', 'Pension Fund', 'Govt Debt', 'Output', 'Net Rate', 'Bp/Y Ratio', 'Residual');
fprintf('%s\n', repmat('-', 1, 130));

% 初始化 t=1 的意外遗赠 (来自 t=0，假设为0)
accidental_bequest_t_minus_1 = 0;

% --- 主循环 (现在从 t=1 到 T-1) ---
for t = 1:(T-1)
    % a. 激活PPS系统
    if sim_years(t) >= cS.pps_activation_year && ~cS.pps_active
        cS.pps_active = true;
        fprintf('--- 激活个人养老金(PPS)系统于 %d (t=%d) ---\n', sim_years(t), t);
    end

    % b. 获取 t 期的宏观环境和供给侧变量 Y(t)
    M_t = main_olg_v14_utils.get_prices_and_policy_at_t(t, paths.K_pvt(t), paths.B_p(t), paths.K_pps(t), paths.B_g(t), Z_path(:,t), A_path(t), cS, paramS, eIdxM);
    paths.Y(t) = M_t.Y_t;
    paths.w(t) = M_t.w_t;
    paths.r_net(t) = M_t.r_net_annual;
    paths.L(t) = M_t.L_t;
    paths.Bp_Y_ratio(t) = paths.B_p(t) / paths.Y(t);

    % c. [核心] 模拟家庭在 t 期的决策
    mass_survived_t = sum(Z_path(:,t) .* cS.s_1yr_transitionV);
    if mass_survived_t > 0
        TR_total_t = accidental_bequest_t_minus_1 / mass_survived_t;
    else
        TR_total_t = 0;
    end
    
    [K_pvt_next, K_pps_next, C_t, Total_Cpps_t, Total_PpsTax_t, accidental_bequest_t] = ...
        main_olg_v14_utils.simulate_private_capital_forward(M_t, Z_path(:,t), cS, paramS, eIdxM, TR_total_t);
    
    paths.C(t) = C_t;
    
    % d. 状态变量演化到 t+1
    paths.K_pvt(t+1) = K_pvt_next;
    paths.K_pps(t+1) = K_pps_next;
    [paths.B_p(t+1), PensionSurplus_t] = main_olg_v14_utils.update_pension_fund(paths.B_p(t), M_t, Z_path(:,t), cS);
    paths.B_g(t+1) = main_olg_v14_utils.update_gov_debt(paths.B_g(t), paths.C(t), M_t, Total_Cpps_t, paths.K_pps(t), Total_PpsTax_t, cS, PensionSurplus_t);
    paths.K_physical(t+1) = paths.K_pvt(t+1) + paths.K_pps(t+1) + paths.B_p(t+1) - paths.B_g(t+1);
    
    % e. 计算投资、政府支出和市场出清残差
    paths.Investment(t) = paths.K_physical(t+1) - (1 - cS.ddk) * paths.K_physical(t);
    paths.G(t) = cS.gov_exp_frac_Y * paths.Y(t);
    paths.market_clearing_residual(t) = paths.Y(t) - (paths.C(t) + paths.Investment(t) + paths.G(t));
    
    % f. 为下一次循环准备遗赠
    accidental_bequest_t_minus_1 = accidental_bequest_t;
    
    % 实时监控
    fprintf('%6d | %8d | %12.3f | %12.3f | %12.3f | %10.3f | %11.2f%% | %11.2f%% | %12.3e\n', ...
        t, sim_years(t), paths.K_physical(t), paths.B_p(t), paths.B_g(t), paths.Y(t), paths.r_net(t)*100, paths.Bp_Y_ratio(t)*100, paths.market_clearing_residual(t));
end

% 处理最后一期 T 的变量
t=T;
M_T = main_olg_v14_utils.get_prices_and_policy_at_t(T, paths.K_pvt(T), paths.B_p(T), paths.K_pps(T), paths.B_g(T), Z_path(:,T), A_path(T), cS, paramS, eIdxM);
paths.Y(T) = M_T.Y_t;
paths.w(T) = M_T.w_t;
paths.r_net(T) = M_T.r_net_annual;
paths.L(T) = M_T.L_t;
paths.Bp_Y_ratio(T) = paths.B_p(T) / paths.Y(T);

% 商品市场出清统计报告
fprintf('\n--- 商品市场出清检验结果 ---\n');
valid_residuals = paths.market_clearing_residual(1:T-1);
valid_Y = paths.Y(1:T-1);
max_abs_residual = max(abs(valid_residuals));
mean_abs_residual = mean(abs(valid_residuals));
if mean(valid_Y) ~= 0
    max_rel_residual_pct = (max_abs_residual / mean(valid_Y)) * 100;
    mean_rel_residual_pct = (mean_abs_residual / mean(valid_Y)) * 100;
    fprintf('最大绝对残差: %.6f (占平均GDP的 %.4f%%)\n', max_abs_residual, max_rel_residual_pct);
    fprintf('平均绝对残差: %.6f (占平均GDP的 %.4f%%)\n', mean_abs_residual, mean_rel_residual_pct);
else
    fprintf('最大绝对残差: %.6f\n', max_abs_residual);
    fprintf('平均绝对残差: %.6f\n', mean_abs_residual);
end
if max_abs_residual < 1e-5
    fprintf('✅ 商品市场近似出清，残差在可接受范围内。\n');
else
    fprintf('⚠️ 商品市场残差较大，请检查代码逻辑。\n');
end

%% 5. 可视化结果
fprintf('\n--- 5. 绘制过渡动态路径 ---\n');
figure('Name', 'OLG v14: 过渡动态与市场出清', 'Position', [100, 100, 1600, 1000]);
time_axis_years = sim_years;
time_axis_residual = sim_years(1:T-1);

subplot(3,3,1); plot(time_axis_years, paths.Y, 'k-', 'LineWidth', 2); title('总产出 Y (5年总量)'); xlabel('年份'); grid on;
subplot(3,3,2); plot(time_axis_years, paths.K_physical, 'b-', 'LineWidth', 2); title('物理总资本 K_{physical}'); xlabel('年份'); grid on;
subplot(3,3,3); plot(time_axis_years, paths.B_p, 'r-', 'LineWidth', 2); title('养老金基金 B_{p}'); xlabel('年份'); grid on;
subplot(3,3,4); plot(time_axis_years, paths.Bp_Y_ratio*100, 'm-', 'LineWidth', 2); title('养老金基金占GDP比例 B_{p}/Y (%)'); xlabel('年份'); ylabel('%'); grid on;
subplot(3,3,5); plot(time_axis_years, paths.K_pps, 'c-', 'LineWidth', 2); title('个人养老金 K_{pps}'); xlabel('年份'); grid on;
subplot(3,3,6); plot(time_axis_years, paths.B_g, 'g-', 'LineWidth', 2); title('政府债务 B_{g}'); xlabel('年份'); grid on;
subplot(3,3,7); plot(time_axis_residual, paths.Investment(1:T-1), 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980]); title('总投资 I'); xlabel('年份'); grid on;
subplot(3,3,8); bar(time_axis_residual, paths.market_clearing_residual(1:T-1), 'FaceColor', [0.5 0.5 0.5]); title('商品市场残差 Y-C-I-G'); xlabel('年份'); grid on; max_res_display = max(abs(paths.market_clearing_residual(1:T-1))) * 1.2; if max_res_display > 0, ylim([-max_res_display, max_res_display]); end
subplot(3,3,9); yyaxis left; plot(time_axis_years, paths.r_net*100, 'b-'); ylabel('年化净利率 (%)'); yyaxis right; plot(time_axis_years, paths.w, 'r-'); ylabel('工资'); title('价格路径'); xlabel('年份'); grid on; legend('r_{net}', 'w', 'Location', 'best');
sgtitle('OLG 模型 v14 过渡动态分析 (1997-2102) - 最终修正版', 'FontSize', 16, 'FontWeight', 'bold');