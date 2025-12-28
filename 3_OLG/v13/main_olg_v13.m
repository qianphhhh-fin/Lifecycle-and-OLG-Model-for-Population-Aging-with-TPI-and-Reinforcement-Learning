% --- START OF FILE main_olg_v13.m ---
%
% OLG 模型 v13: 基于真实数据校准和外生路径的过渡动态分析
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
fprintf('=== OLG 模型 v13: 基于数据校准的过渡动态分析 ===\n');

%% 1. 初始化参数和环境
fprintf('\n--- 1. 初始化参数 ---\n');
cS = main_olg_v13_utils.ParameterValues_HuggettStyle();
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
cS.alpha = 0.55; % [强烈建议] 大幅提高资本份额以解决负利率问题
cS.ddk = 1 - (1 - 0.05)^cS.time_step; % 5年期折旧率 (年化8%)

% --- 政府财政参数 (基于start_year年校准) ---
cS.tau_k = 0.25;
cS.tau_l = 0.08;
cS.tau_c = 0.04;
cS.gov_exp_frac_Y = 0.15;

% --- 网格参数 ---
ngrid = 30;
cS.nk = ngrid;
cS.nkpps = ngrid;
cS.nkprime = ngrid;
cS.npps = ngrid;
cS.A = 1.0; % [v13修正] 为网格生成提供一个基准A值
cS = main_olg_v13_utils.generateGrids(cS);
    
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

% [新增] 调试参数
cS.debug_period = 2; % 设置为 > 0 的整数来激活对该期的详细调试
                     % 设置为 0 或-1 来关闭调试

% [新的调用方式]
fprintf('正在构建基于现实制度差异的有效养老金缴费率路径 (theta_path)...\n');
% 将时间轴参数和绘图开关传递给新函数
cS = main_olg_v13_utils.calcaulte_theta_payg_path(cS, true);
fprintf('✅ 复合有效缴费率路径构建完成。\n');

%% 2. [v13核心] 加载外生数据路径
[Z_path, A_path, T] = main_olg_v13_utils.load_exogenous_paths(cS);
sim_years = cS.start_year:cS.time_Step:(cS.start_year + (T-1)*cS.time_Step);

% --- 生成劳动效率冲击路径 (与之前相同) ---
[paramS.leLogGridV, paramS.leTrProbM, paramS.leProb1V] = main_olg_v13_utils.EarningProcess_olgm(cS);
paramS.leGridV = exp(paramS.leLogGridV(:));
paramS.ageEffV_new = cS.ageEffV_new;
eIdxM = main_olg_v13_utils.LaborEndowSimulation_olgm_AgeGroup(cS, paramS);
fprintf('✅ 效率冲击路径生成完成。\n');

% --- 在 main_olg_v13.m 的第2节 "初始化模型参数" 中 ---



%% 3. [v13核心 - 已修正单位统一的校准] 基于数据校准start_year年的初始状态 (t=1)
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
paths.K_physical = zeros(1, T); % 新增物理资本路径

% --- 校准参数 ---
KY_ratio_start_year = 1.994;
Bg_Y_ratio_start_year = 0.08;
A_1 = A_path(1);
Z_1 = Z_path(:, 1);

% --- [核心修正] 校准循环在【年化】单位下进行，以避免放大效应 ---
fprintf('   通过迭代寻找自洽的【年化】Y_annual(1) 和 K_physical(1)...\n');
Y_annual_guess = 1.5; % 从一个合理的【年化】Y猜测开始
tol = 1e-7;
max_iter = 100;
damp = 0.7; % 阻尼系数

for iter = 1:max_iter
    % a. 基于 Y_annual_guess 计算目标宏观存量
    %    K/Y 和 B_g/Y 比率通常是基于年化GDP，所以这里计算正确
    K_physical_target = KY_ratio_start_year * Y_annual_guess;
    B_g_target = Bg_Y_ratio_start_year * Y_annual_guess;
    
    % b. 模型设定：初始时刻社保基金和PPS为零
    B_p_1 = 0;
    K_pps_1 = 0;
    
    % c. [正确核算] 倒推出 K_pvt_1
    K_pvt_1 = K_physical_target + B_g_target - B_p_1 - K_pps_1;
    
    if K_pvt_1 < 0
        error('初始校准失败：计算出的初始私人资本 K_pvt(1) 为负，请检查校准目标参数(K/Y, B_g/Y)。');
    end
    
    % d. 基于当前的资本存量，计算价格和新的【年化】Y
    M_1_iter = main_olg_v13_utils.get_prices_and_policy_at_t(1, K_pvt_1, B_p_1, K_pps_1, B_g_target, Z_1, A_1, cS, paramS, eIdxM);
    
    % [核心] get_prices...返回的 M_t.Y_t 是5年期总量，我们需要用年化产出进行比较
    Y_annual_new = M_1_iter.Y_t / cS.time_Step; 
    
    % e. 检查收敛
    if abs(Y_annual_new - Y_annual_guess) < tol
        fprintf('   ✅ 自洽的初始【年化】状态已找到 (Y_annual=%.4f, K_phys=%.4f)。\n', Y_annual_new, K_physical_target);
        break;
    end
    
    % f. 更新猜测值
    Y_annual_guess = (1 - damp) * Y_annual_guess + damp * Y_annual_new;
end
if iter == max_iter
    warning('初始状态迭代在%d次内未收敛。', max_iter);
end

% --- 将收敛后的值赋给路径 ---
Y_annual_final = Y_annual_new;
paths.Y(1) = Y_annual_final * cS.time_Step; % 存储为5年期总量
paths.K_physical(1) = KY_ratio_start_year * Y_annual_final;
paths.B_g(1) = Bg_Y_ratio_start_year * Y_annual_final;
paths.B_p(1) = 0;
paths.K_pps(1) = 0;
paths.K_pvt(1) = paths.K_physical(1) + paths.B_g(1) - paths.B_p(1) - paths.K_pps(1);

% 获取最终的价格和宏观变量
M_1 = main_olg_v13_utils.get_prices_and_policy_at_t(1, paths.K_pvt(1), paths.B_p(1), paths.K_pps(1), paths.B_g(1), Z_1, A_1, cS, paramS, eIdxM);
paths.L(1) = M_1.L_t;
paths.w(1) = M_1.w_t;
paths.r_net(1) = M_1.r_net_t;
paths.G(1) = cS.gov_exp_frac_Y * paths.Y(1);
paths.Bp_Y_ratio(1) = paths.B_p(1) / paths.Y(1);

fprintf('✅ start_year年初始状态校准完成: Y(5-yr)=%.3f, K_physical=%.3f, K_pvt=%.3f, B_g=%.3f\n', ...
    paths.Y(1), paths.K_physical(1), paths.K_pvt(1), paths.B_g(1));
%% 4. [v13核心 - 最终修正版平滑启动] 过渡动态主循环 (t=1 to T-1)
fprintf('\n--- 4. 开始向前模拟过渡动态 (start_year - %d) ---\n', sim_years(end));

% 打印表头
fprintf('\n%6s | %8s | %12s | %12s | %12s | %10s | %12s | %12s | %12s\n', ...
    'Time', 'Year', 'Physical K', 'Pension Fund', 'Govt Debt', 'Output', 'Net Rate', 'Bp/Y Ratio', 'Residual');
fprintf('%s\n', repmat('-', 1, 130));


% --- 步骤 1: 单独处理 t=1，实现平滑启动 ---
t = 1;
% a. 获取 t=1 的宏观环境 (已在第3节校准)
M_1 = main_olg_v13_utils.get_prices_and_policy_at_t(t, paths.K_pvt(1), paths.B_p(1), paths.K_pps(1), paths.B_g(1), Z_path(:,t), A_path(t), cS, paramS, eIdxM);

% b. 模拟 t=1 的家庭决策，得到 t+1 的私人资本 和 t 的模拟消费
[K_pvt_2, K_pps_2, C_t1_sim, Total_Cpps_t1, Total_PpsTax_t1] = main_olg_v13_utils.simulate_private_capital_forward(M_1, Z_path(:,t), cS, paramS, eIdxM);

% c. [反推法]
%    首先，使用模拟出的消费 C_t1_sim 来计算一个临时的 t=2 状态，目的是得到一个一致的投资 I(1)
B_p_2_temp = main_olg_v13_utils.update_pension_fund(paths.B_p(1), M_1, Z_path(:,t), cS);
B_g_2_temp = main_olg_v13_utils.update_gov_debt(paths.B_g(1), C_t1_sim, M_1, Total_Cpps_t1, paths.K_pps(1), Total_PpsTax_t1, cS);
K_physical_2_temp = K_pvt_2 + K_pps_2 + B_p_2_temp - B_g_2_temp;

% d. 计算 t=1 的投资 I(1) 和 政府支出 G(1)
paths.Investment(1) = K_physical_2_temp - (1 - cS.ddk) * paths.K_physical(1);
paths.G(1) = cS.gov_exp_frac_Y * paths.Y(1);

% e. 根据商品市场出清，计算 t=1 的【最终】消费 C(1)
paths.C(1) = paths.Y(1) - paths.Investment(1) - paths.G(1);

% f. [最终确定 t=2 状态]
%    使用这个最终的 C(1) 来计算 t=2 的所有状态变量，确保主循环的起点是完全一致的。
paths.K_pvt(2) = K_pvt_2;
paths.K_pps(2) = K_pps_2;
paths.B_p(2) = main_olg_v13_utils.update_pension_fund(paths.B_p(1), M_1, Z_path(:,t), cS);
paths.B_g(2) = main_olg_v13_utils.update_gov_debt(paths.B_g(1), paths.C(1), M_1, Total_Cpps_t1, paths.K_pps(1), Total_PpsTax_t1, cS);
paths.K_physical(2) = paths.K_pvt(2) + paths.K_pps(2) + paths.B_p(2) - paths.B_g(2);

% g. [最终检验] 重新计算一次I(1)和残差，确保它们为零
paths.Investment(1) = paths.K_physical(2) - (1-cS.ddk)*paths.K_physical(1);
paths.C(1) = paths.Y(1) - paths.Investment(1) - paths.G(1);
paths.market_clearing_residual(1) = paths.Y(1) - paths.C(1) - paths.Investment(1) - paths.G(1);

% 打印 t=1 的最终状态
fprintf('%6d | %8d | %12.3f | %12.3f | %12.3f | %10.3f | %11.2f%% | %11.2f%% | %12.3e\n', ...
    t, sim_years(t), paths.K_physical(1), paths.B_p(1), paths.B_g(1), paths.Y(1), paths.r_net(1)*100, paths.Bp_Y_ratio(1)*100, paths.market_clearing_residual(1));

% --- 步骤 2: 向前模拟主循环 (从 t=2 到 T-1) ---
for t = 2:(T-1)
    % a. 激活PPS系统
    if sim_years(t) >= cS.pps_activation_year && ~cS.pps_active
        cS.pps_active = true;
        fprintf('--- 激活个人养老金(PPS)系统于 %d (t=%d) ---\n', sim_years(t), t);
    end

    % b. 获取 t 期的宏观环境 Y(t) 等
    M_t = main_olg_v13_utils.get_prices_and_policy_at_t(t, paths.K_pvt(t), paths.B_p(t), paths.K_pps(t), paths.B_g(t), Z_path(:,t), A_path(t), cS, paramS, eIdxM);
    paths.Y(t) = M_t.Y_t;
    paths.w(t) = M_t.w_t;
    paths.r_net(t) = M_t.r_net_t;
    paths.L(t) = M_t.L_t;
    paths.Bp_Y_ratio(t) = paths.B_p(t) / paths.Y(t);

    % c. 模拟家庭决策，得到 C(t) 和 t+1 的私人资本
    [K_pvt_next, K_pps_next, C_t, Total_Cpps_t, Total_PpsTax_t] = main_olg_v13_utils.simulate_private_capital_forward(M_t, Z_path(:,t), cS, paramS, eIdxM);
    paths.C(t) = C_t;
    
    % d. 演化到 t+1
    paths.K_pvt(t+1) = K_pvt_next;
    paths.K_pps(t+1) = K_pps_next;
    paths.B_p(t+1) = main_olg_v13_utils.update_pension_fund(paths.B_p(t), M_t, Z_path(:,t), cS);
    paths.B_g(t+1) = main_olg_v13_utils.update_gov_debt(paths.B_g(t), paths.C(t), M_t, Total_Cpps_t, paths.K_pps(t), Total_PpsTax_t, cS);
    paths.K_physical(t+1) = paths.K_pvt(t+1) + paths.K_pps(t+1) + paths.B_p(t+1) - paths.B_g(t+1);
    
    % e. 计算 I(t), G(t)
    paths.Investment(t) = paths.K_physical(t+1) - (1 - cS.ddk) * paths.K_physical(t);
    paths.G(t) = cS.gov_exp_frac_Y * paths.Y(t);
    
    % f. 计算残差
    paths.market_clearing_residual(t) = paths.Y(t) - (paths.C(t) + paths.Investment(t) + paths.G(t));
    
    % [调用调试函数]
    if t == cS.debug_period
        main_olg_v13_utils.debug_national_accounts(t, paths, M_t, cS);
    end
    
    % 实时监控
    fprintf('%6d | %8d | %12.3f | %12.3f | %12.3f | %10.3f | %11.2f%% | %11.2f%% | %12.3e\n', ...
        t, sim_years(t), paths.K_physical(t), paths.B_p(t), paths.B_g(t), paths.Y(t), paths.r_net(t)*100, paths.Bp_Y_ratio(t)*100, paths.market_clearing_residual(t));
end


fprintf('%s\n', repmat('-', 1, 130));
fprintf('✅ 过渡动态模拟完成。\n');

% 处理最后一期 T 的变量 (没有投资和残差)
t=T;
M_T = main_olg_v13_utils.get_prices_and_policy_at_t(T, paths.K_pvt(T), paths.B_p(T), paths.K_pps(T), paths.B_g(T), Z_path(:,T), A_path(T), cS, paramS, eIdxM);
paths.Y(T) = M_T.Y_t;
paths.w(T) = M_T.w_t;
paths.r_net(T) = M_T.r_net_t;
paths.L(T) = M_T.L_t;
paths.Bp_Y_ratio(T) = paths.B_p(T) / paths.Y(T);

% [新增] 商品市场出清统计报告
fprintf('\n--- 商品市场出清检验结果 ---\n');
valid_residuals = paths.market_clearing_residual(1:T-1);
valid_Y = paths.Y(1:T-1);
max_abs_residual = max(abs(valid_residuals));
mean_abs_residual = mean(abs(valid_residuals));
max_rel_residual_pct = (max_abs_residual / mean(valid_Y)) * 100;
mean_rel_residual_pct = (mean_abs_residual / mean(valid_Y)) * 100;

fprintf('最大绝对残差: %.6f (占平均GDP的 %.4f%%)\n', max_abs_residual, max_rel_residual_pct);
fprintf('平均绝对残差: %.6f (占平均GDP的 %.4f%%)\n', mean_abs_residual, mean_rel_residual_pct);
if max_abs_residual < 1e-5
    fprintf('✅ 商品市场近似出清，残差在可接受范围内。\n');
else
    fprintf('⚠️ 商品市场残差较大，请检查代码逻辑，特别是资产演化方程。\n');
end

%% 5. 可视化结果
fprintf('\n--- 5. 绘制过渡动态路径 ---\n');
figure('Name', 'OLG v13: 过渡动态与市场出清', 'Position', [100, 100, 1600, 1000]);
time_axis_years = sim_years;
time_axis_residual = sim_years(1:T-1); % 残差和投资只到 T-1

subplot(3,3,1);
plot(time_axis_years, paths.Y, 'k-', 'LineWidth', 2);
title('总产出 Y'); xlabel('年份'); grid on;

subplot(3,3,2);
plot(time_axis_years, paths.K_physical, 'b-', 'LineWidth', 2);
title('物理总资本 K_{physical}'); xlabel('年份'); grid on;

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
plot(time_axis_residual, paths.Investment(1:T-1), 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980]);
title('总投资 I'); xlabel('年份'); grid on;

% [修改] 绘制市场出清残差图
subplot(3,3,8);
bar(time_axis_residual, paths.market_clearing_residual(1:T-1), 'FaceColor', [0.5 0.5 0.5]);
title('商品市场残差 Y-C-I-G'); xlabel('年份'); grid on;
% 自动调整Y轴以更好地显示小的残差
max_res_display = max(abs(paths.market_clearing_residual(1:T-1))) * 1.2;
if max_res_display > 0, ylim([-max_res_display, max_res_display]); end


subplot(3,3,9);
yyaxis left
plot(time_axis_years, paths.r_net*100, 'b-');
ylabel('净利率 (%)');
yyaxis right
plot(time_axis_years, paths.w, 'r-');
ylabel('工资');
title('价格路径'); xlabel('年份'); grid on;
legend('r_{net}', 'w', 'Location', 'best');

sgtitle('OLG 模型 v13 过渡动态分析 (1997-2102) - 含市场出清检验', 'FontSize', 16, 'FontWeight', 'bold');