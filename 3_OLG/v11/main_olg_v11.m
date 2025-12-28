% --- START OF FILE main_olg_v11.m ---
%
% OLG 模型 V11: 基于真实数据校准和外生路径的过渡动态分析
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
fprintf('=== OLG 模型 V11: 基于数据校准的过渡动态分析 ===\n');

%% 1. 初始化参数和环境
fprintf('\n--- 1. 初始化参数 ---\n');
cS = main_olg_v11_utils.ParameterValues_HuggettStyle();
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
cS.alpha = 0.45; % [强烈建议] 大幅提高资本份额以解决负利率问题
cS.ddk = 1 - (1 - 0.08)^cS.time_step; % 5年期折旧率 (年化8%)

% --- 政府财政参数 (基于start_year年校准) ---
cS.tau_k = 0.25;
cS.tau_l = 0.08;
cS.tau_c = 0.04;
cS.gov_exp_frac_Y = 0.15;

% --- 网格参数 ---
ngrid = 50;
cS.nk = ngrid;
cS.nkpps = ngrid;
cS.nkprime = ngrid;
cS.npps = ngrid;
cS.A = 1.0; % [V11修正] 为网格生成提供一个基准A值
cS = main_olg_v11_utils.generateGrids(cS);
    
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
cS = main_olg_v11_utils.calcaulte_theta_payg_path(cS, true);
fprintf('✅ 复合有效缴费率路径构建完成。\n');

%% 2. [V11核心] 加载外生数据路径
[Z_path, A_path, T] = main_olg_v11_utils.load_exogenous_paths(cS);
sim_years = cS.start_year:cS.time_Step:(cS.start_year + (T-1)*cS.time_Step);

% --- 生成劳动效率冲击路径 (与之前相同) ---
[paramS.leLogGridV, paramS.leTrProbM, paramS.leProb1V] = main_olg_v11_utils.EarningProcess_olgm(cS);
paramS.leGridV = exp(paramS.leLogGridV(:));
paramS.ageEffV_new = cS.ageEffV_new;
eIdxM = main_olg_v11_utils.LaborEndowSimulation_olgm_AgeGroup(cS, paramS);
fprintf('✅ 效率冲击路径生成完成。\n');

% --- 在 main_olg_v11.m 的第2节 "初始化模型参数" 中 ---



%% 3. [V11核心 - 已修正] 基于数据校准start_year年的初始状态 (t=1)
fprintf('\n--- 3. 校准start_year年 (t=1) 的初始状态 ---\n');

% 初始化所有状态变量的路径
paths = struct();
paths.K_pvt = zeros(1, T);
paths.K_pps = zeros(1, T);
paths.B_p = zeros(1, T);
paths.B_g = zeros(1, T);
paths.w = zeros(1, T);
paths.r_net = zeros(1, T);
paths.L = zeros(1, T);
paths.Y = zeros(1, T);
paths.G = zeros(1, T);
paths.C = zeros(1, T);
paths.Investment = zeros(1, T-1);
paths.market_clearing_residual = zeros(1, T-1);
paths.Bp_Y_ratio = zeros(1, T);

% --- 校准参数 ---
KY_ratio_start_year = 1.994;
Bg_Y_ratio_start_year = 0.08;
A_1 = A_path(1);
Z_1 = Z_path(:, 1);

% --- [核心修正] 使用迭代法寻找自洽的初始宏观状态 ---
fprintf('   通过迭代寻找自洽的 Y(1) 和 K(1)...\n');
Y_guess = 1.0; % 从一个任意的Y猜测开始
for iter = 1:50
    K_guess = KY_ratio_start_year * Y_guess;
    
    % 分解资本
    B_p_1 = 0;
    K_pps_1 = 0;
    B_g_1 = Bg_Y_ratio_start_year * Y_guess; % 债务也与Y的猜测值挂钩
    K_pvt_1 = K_guess - B_p_1 - K_pps_1 - B_g_1;
    
    % 基于当前的K, 计算价格和新的Y
    M_1_iter = main_olg_v11_utils.get_prices_and_policy_at_t(1, K_pvt_1, B_p_1, K_pps_1, Z_1, A_1, cS, paramS, eIdxM);
    Y_new = M_1_iter.Y_t;
    
    % 检查收敛
    if abs(Y_new - Y_guess) < 1e-7
        fprintf('   ✅ 自洽的初始状态已找到 (Y=%.4f)。\n', Y_new);
        break;
    end
    
    % 更新猜测值
    Y_guess = 0.5 * Y_guess + 0.5 * Y_new;
end
if iter == 20
    warning('初始状态迭代未收敛。');
end

% --- 将收敛后的值赋给路径 ---
paths.Y(1) = Y_new;
K_total_1 = KY_ratio_start_year * paths.Y(1);
paths.B_p(1) = 0;
paths.K_pps(1) = 0;
paths.B_g(1) = Bg_Y_ratio_start_year * paths.Y(1);
paths.K_pvt(1) = K_total_1 - paths.B_p(1) - paths.K_pps(1) - paths.B_g(1);

% 获取最终的价格和宏观变量
M_1 = main_olg_v11_utils.get_prices_and_policy_at_t(1, paths.K_pvt(1), paths.B_p(1), paths.K_pps(1), Z_1, A_1, cS, paramS, eIdxM);
paths.L(1) = M_1.L_t;
paths.w(1) = M_1.w_t;
paths.r_net(1) = M_1.r_net_t;
paths.G(1) = cS.gov_exp_frac_Y * paths.Y(1);
paths.Bp_Y_ratio(1) = paths.B_p(1) / paths.Y(1);

fprintf('✅ start_year年初始状态校准完成: Y=%.3f, K_total=%.3f, K_pvt=%.3f, B_g=%.3f\n', ...
    paths.Y(1), K_total_1, paths.K_pvt(1), paths.B_g(1));

%% 4. [V11核心] 过渡动态主循环 (t=1 to T-1)
fprintf('\n--- 4. 开始向前模拟过渡动态 (start_year - %d) ---\n', sim_years(end));

% 打印表头
fprintf('\n%6s | %8s | %12s | %12s | %12s | %12s | %10s | %12s | %12s\n', ...
    'Time', 'Year', 'Private K', 'PPS K', 'Pension Fund', 'Govt Debt', 'Output', 'Net Rate', 'Bp/Y Ratio');
fprintf('%s\n', repmat('-', 1, 120));

% t=1 的状态需要单独处理，因为它决定了 t=2 的状态
% 我们需要计算出 C(1) 和下一期的资本存量
t = 1;
M_t = main_olg_v11_utils.get_prices_and_policy_at_t(t, paths.K_pvt(t), paths.B_p(t), paths.K_pps(t), Z_path(:,t), A_path(t), cS, paramS, eIdxM);

[K_pvt_next, K_pps_next, C_t, Total_Cpps_t, Total_PpsWithdrawalTax_t] = main_olg_v11_utils.simulate_private_capital_forward(M_t, Z_path(:,t), cS, paramS, eIdxM);
paths.C(1) = C_t; % 记录当期消费

% 打印 t=1 (start_year年) 的状态
fprintf('%6d | %8d | %12.3f | %12.3f | %12.3f | %12.3f | %10.3f | %11.2f%% | %11.2f%%\n', ...
    t, sim_years(t), paths.K_pvt(t), paths.K_pps(t), paths.B_p(t), paths.B_g(t), paths.Y(t), paths.r_net(t)*100, paths.Bp_Y_ratio(t)*100);

% 向前模拟主循环 (现在从 t=2 到 T)
for t = 2:T
    % a. 激活PPS系统（如果到达指定年份）
    if sim_years(t) >= cS.pps_activation_year && ~cS.pps_active
        cS.pps_active = true;
        fprintf('--- 激活个人养老金(PPS)系统于 %d (t=%d) ---\n', sim_years(t), t);
    end

    % b. 获取当前状态 (t-1 的决策决定了 t 的状态)
    paths.K_pvt(t) = K_pvt_next;
    paths.K_pps(t) = K_pps_next;
    paths.B_p(t) = main_olg_v11_utils.update_pension_fund(paths.B_p(t-1), M_t, Z_path(:,t-1), cS);
    paths.B_g(t) = main_olg_v11_utils.update_gov_debt(paths.B_g(t-1), paths.C(t-1), M_t, Total_Cpps_t, paths.K_pps(t-1), Total_PpsWithdrawalTax_t, cS);

    % c. 计算当期价格和宏观变量
    M_t = main_olg_v11_utils.get_prices_and_policy_at_t(t, paths.K_pvt(t), paths.B_p(t), paths.K_pps(t), Z_path(:,t), A_path(t), cS, paramS, eIdxM);
    
    % d. 模拟家庭决策，得到下一期资本和当期消费
    [K_pvt_next, K_pps_next, C_t, Total_Cpps_t, Total_PpsWithdrawalTax_t] = main_olg_v11_utils.simulate_private_capital_forward(M_t, Z_path(:,t), cS, paramS, eIdxM);
    
    % e. 记录当期宏观变量
    paths.w(t) = M_t.w_t;
    paths.r_net(t) = M_t.r_net_t;
    paths.L(t) = M_t.L_t;
    paths.Y(t) = M_t.Y_t;
    paths.C(t) = C_t;
    paths.Bp_Y_ratio(t) = paths.B_p(t) / paths.Y(t);
    
    % f. 计算投资和检验商品市场出清 (针对 t-1 时期)
    K_total_t_minus_1 = paths.K_pvt(t-1) + paths.B_p(t-1) + paths.K_pps(t-1);
    K_total_t = paths.K_pvt(t) + paths.B_p(t) + paths.K_pps(t);
    paths.Investment(t-1) = K_total_t - (1 - cS.ddk) * K_total_t_minus_1;
    paths.G(t-1) = cS.gov_exp_frac_Y * paths.Y(t-1);
    paths.market_clearing_residual(t-1) = paths.Y(t-1) - (paths.C(t-1) + paths.Investment(t-1) + paths.G(t-1));
    
    % 实时监控
    fprintf('%6d | %8d | %12.3f | %12.3f | %12.3f | %12.3f | %10.3f | %11.2f%% | %11.2f%%\n', ...
        t, sim_years(t), paths.K_pvt(t), paths.K_pps(t), paths.B_p(t), paths.B_g(t), paths.Y(t), paths.r_net(t)*100, paths.Bp_Y_ratio(t)*100);
end

fprintf('%s\n', repmat('-', 1, 120));
fprintf('✅ 过渡动态模拟完成。\n');

% 商品市场出清统计
max_residual = max(abs(paths.market_clearing_residual));
mean_residual = mean(abs(paths.market_clearing_residual));
fprintf('\n--- 商品市场出清检验结果 ---\n');
fprintf('最大残差: %.6f (%.3f%% of 平均GDP)\n', max_residual, max_residual/mean(paths.Y(1:T-1))*100);
fprintf('平均残差: %.6f (%.3f%% of 平均GDP)\n', mean_residual, mean_residual/mean(paths.Y(1:T-1))*100);
if max_residual < 0.001 * mean(paths.Y(1:T-1))
    fprintf('✅ 商品市场近似出清，残差在可接受范围内。\n');
else
    fprintf('⚠️ 商品市场残差较大，可能存在数值问题。\n');
end

%% 5. 可视化结果
fprintf('\n--- 5. 绘制过渡动态路径 ---\n');
figure('Name', 'OLG V11: 基于数据校准的过渡动态', 'Position', [100, 100, 1600, 1000]);
time_axis_years = sim_years;

subplot(3,3,1);
plot(time_axis_years, paths.Y, 'k-', 'LineWidth', 2);
title('总产出 Y (start_year=1)'); xlabel('年份'); grid on;

subplot(3,3,2);
plot(time_axis_years, paths.K_total, 'b-', 'LineWidth', 2);
title('总资本 K_{total}'); xlabel('年份'); grid on;

subplot(3,3,3);
plot(time_axis_years, paths.B_p, 'r-', 'LineWidth', 2);
title('养老金基金 B_{p}'); xlabel('年份'); grid on;

subplot(3,3,4);
plot(time_axis_years, paths.Bp_Y_ratio*100, 'm-', 'LineWidth', 2);
title('养老金基金占GDP比例 B_{p}/Y (%)'); xlabel('年份'); ylabel('%'); grid on;

subplot(3,3,5);
plot(time_axis_years, paths.K_pps, 'c-', 'LineWidth', 2);
title('个人养老金 K_{pps}'); xlabel('年份'); grid on;

subplot(3,3,6);
plot(time_axis_years, paths.B_g, 'g-', 'LineWidth', 2);
title('政府债务 B_{g}'); xlabel('年份'); grid on;

subplot(3,3,7);
plot(time_axis_years(1:T-1), paths.Investment, 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980]);
title('总投资 I'); xlabel('年份'); grid on;

subplot(3,3,8);
plot(time_axis_years(1:T-1), paths.market_clearing_residual, 'LineWidth', 2, 'Color', [0.5 0.5 0.5]);
title('商品市场残差 Y-C-I-G'); xlabel('年份'); grid on;

subplot(3,3,9);
yyaxis left
plot(time_axis_years, paths.r_net*100, 'b-');
ylabel('净利率 (%)');
yyaxis right
plot(time_axis_years, paths.w, 'r-');
ylabel('工资');
title('价格路径'); xlabel('年份'); grid on;
legend('r_{net}', 'w', 'Location', 'best');

sgtitle('OLG 模型 V11 过渡动态分析 (1997-2102)', 'FontSize', 16, 'FontWeight', 'bold');