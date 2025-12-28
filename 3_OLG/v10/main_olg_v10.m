% --- START OF FILE main_olg_v10.m (v10.3 - Final with SS Fund Calibration) ---

clear; close all;
addpath(pwd);
fprintf('=== OLG 模型 V10: 过渡动态分析 (基于含基金的初始稳态) ===\n');

%% 1. 初始化参数和环境
fprintf('\n--- 1. 初始化 ---\n');
T = 150; 
cS = main_olg_v10_utils.ParameterValues_HuggettStyle();
cS.nSim = 5000;
paramS = struct();

% --- V10 固定的政策参数 ---
cS.rho_prime_payg = 0.4; 
cS.theta_payg = 0.15; 
B_p_Y_ratio_target = 0.05; % [新] 目标：初始稳态养老金基金占GDP的5%

% --- [对齐] 家庭偏好参数 ---
cS.sigma = 3.0;     % 一个标准的、温和的风险厌恶系数。IES=1/3，行为稳健。
cS.beta = 0.97;     % 非常有耐心的家庭（5年期贴现因子）。
cS.cFloor = 0.05;
cS.nSim = 1000;
% 如果您已加入遗赠动机
cS.phi_bequest = 3; 
cS.sigma_bequest = cS.sigma;

% --- [对齐] 生产技术参数 ---
cS.A = 2;
cS.alpha = 0.36;
% 关键修改：设定一个与5年期匹配的、基于现实的折旧率
% 基于年化折旧率 delta_annual = 8% (0.08)
cS.ddk = 0.3; %1 - (1 - 0.08)^5;  % 约等于 0.341

% --- [对齐] 政府财政参数 ---
cS.tau_k = 0.20;    % 基于中国现实和主流文献的资本税率
cS.tau_c = 0.02;
cS.tau_l = 0.10;
cS.gov_exp_frac_Y = 0.15;
cS.gov_debt_frac_Y = 0.60;

ngrid = 50;
cS.nk = ngrid;
cS.nkpps = ngrid;
cS.nkprime = ngrid;
cS.npps = ngrid;
    
% --- 财政反馈规则参数 ---
cS.debt_ratio_target = 0.60;
cS.debt_feedback_strength = 0.5;


% --- [对齐] PPS制度参数 ---
cS.pps_active = true;
cS.pps_tax_rate_withdrawal = 0.03;
cS.pps_return_rate_premium = 0.01;
cS.pps_withdrawal_rate = 0.15;
cS.pps_in_K = true;
cS.pps_bequeathable = true;
cS.pps_contrib_limit = 9999;
cS.pps_max_contrib_frac = 0.1;

% 重新生成网格
cS = main_olg_v10_utils.generateGrids(cS);
fprintf('风险厌恶系数 sigma=%.1f, 目标初始养老金基金/GDP=%.1f%%\n', cS.sigma, B_p_Y_ratio_target*100);

%% 2. 计算外生路径
% ... (这部分完全不变) ...
fprintf('\n--- 2. 计算外生路径 ---\n');
popS = main_olg_v10_utils.initPopulation(cS);
popS = main_olg_v10_utils.populationDynamics(popS, cS);
Z_path = zeros(cS.aD_new, T + 1);
num_pop_periods = size(popS.Z, 2);
for t = 1:(T + 1)
    if t < num_pop_periods, Z_path(:, t) = popS.Z(:, t);
    else, Z_path(:, t) = popS.Z(:, end); end
end
Z_path_norm = Z_path ./ sum(Z_path, 1);
fprintf('✅ 人口过渡路径计算完成。\n');
[paramS.leLogGridV, paramS.leTrProbM, paramS.leProb1V] = main_olg_v10_utils.EarningProcess_olgm(cS);
paramS.leGridV = exp(paramS.leLogGridV(:));
paramS.ageEffV_new = cS.ageEffV_new;
eIdxM = main_olg_v10_utils.LaborEndowSimulation_olgm_AgeGroup(cS, paramS);
fprintf('✅ 效率冲击路径生成完成。\n');


%% 3. [核心] 求解初始稳态 (t=0)
fprintf('\n--- 3. 求解初始稳态 (有目标基金, 无初始PPS) ---\n');
Z_0 = Z_path_norm(:, 1);
[ss0, eq_found] = main_olg_v10_utils.solve_steady_state_with_fund(Z_0, B_p_Y_ratio_target, cS, paramS, eIdxM);
if ~eq_found, error('无法求解初始稳态，模拟无法继续。'); end

fprintf('✅ 初始稳态求解完成: K_total=%.3f, K_pvt=%.3f, B_p=%.3f, K_pps=%.3f, r_net=%.2f%%\n', ...
    ss0.K_total, ss0.K_pvt, ss0.B_p, ss0.K_pps, ss0.r_net*100);

%% 4. 过渡动态主循环
fprintf('\n--- 4. 开始向前模拟过渡动态 (PPS系统已激活) ---\n');

% [重要] 确保在过渡动态中，PPS系统是激活的
cS.pps_active = true;

% 初始化所有状态变量的路径
paths = struct();
paths.K_pvt = zeros(1, T + 1);
paths.K_pps = zeros(1, T + 1);
paths.B_p = zeros(1, T + 1);
paths.B_g = zeros(1, T + 1);
paths.w = zeros(1, T + 1);
paths.r_net = zeros(1, T + 1);
paths.L = zeros(1, T + 1);
paths.Y = zeros(1, T + 1);
paths.G = zeros(1, T + 1);
paths.C = zeros(1, T + 1);
paths.Investment = zeros(1, T);
paths.market_clearing_residual = zeros(1, T);

% [核心] 使用求解出的稳态值设定 t=0 的初始状态
paths.K_pvt(1) = ss0.K_pvt;
paths.K_pps(1) = ss0.K_pps;
paths.B_p(1) = ss0.B_p;
paths.B_g(1) = ss0.B_g;
% 记录 t=0 时的其他宏观变量
paths.w(1) = ss0.w;
paths.r_net(1) = ss0.r_net;
paths.L(1) = ss0.L;
paths.Y(1) = ss0.Y;
paths.G(1) = ss0.G;
paths.C(1) = ss0.C;

% 打印表头
fprintf('\n%s\n', repmat('=', 1, 130));
fprintf('%6s | %12s | %12s | %12s | %12s | %10s | %10s | %12s | %10s\n', ...
    'Time', 'Private K', 'PPS K', 'Pension Fund', 'Govt Debt', 'Output', 'Wage', 'Net Rate', 'Gov Spend');
fprintf('%6s | %12s | %12s | %12s | %12s | %10s | %10s | %12s | %10s\n', ...
    '(t)', '(K_pvt)', '(K_pps)', '(B_p)', '(B_g)', '(Y)', '(w)', '(r_net)', '(G)');
fprintf('%s\n', repmat('-', 1, 130));
% 打印 t=0 的状态
fprintf('%6d | %12.3f | %12.3f | %12.3f | %12.3f | %10.3f | %10.3f | %11.2f%% | %10.3f\n', ...
    0, paths.K_pvt(1), paths.K_pps(1), paths.B_p(1), paths.B_g(1), paths.Y(1), paths.w(1), paths.r_net(1)*100, paths.G(1));

% 向前模拟主循环 (现在从 t=1 到 T)
for t = 1:T
    % a. 获取当前状态 (t-1 的决策决定了 t 的状态)
    K_pvt_t = paths.K_pvt(t);
    K_pps_t = paths.K_pps(t);
    B_p_t = paths.B_p(t);
    B_g_t = paths.B_g(t);
    Z_t_norm = Z_path_norm(:, t);
    
    % b. 计算当期价格和宏观变量
    M_t = main_olg_v10_utils.get_prices_and_policy_at_t(K_pvt_t, B_p_t, K_pps_t, Z_t_norm, cS, paramS, eIdxM);
    M_t.current_t = t;
    
    % c. 模拟家庭决策，得到下一期私人资本和当期消费
    [paths.K_pvt(t+1), paths.K_pps(t+1), paths.C(t), Total_Cpps_t, Total_PpsWithdrawalTax_t] = main_olg_v10_utils.simulate_private_capital_forward(M_t, Z_t_norm, cS, paramS, eIdxM);
    
    % d. 更新其他状态变量
    paths.B_p(t+1) = main_olg_v10_utils.update_pension_fund(B_p_t, M_t, Z_t_norm, cS);
    [paths.B_g(t+1), paths.G(t)] = main_olg_v10_utils.update_gov_debt(B_g_t, paths.C(t),  M_t,Total_Cpps_t, K_pps_t,Total_PpsWithdrawalTax_t, cS, paramS);

    % e. 记录当期宏观变量
    paths.w(t) = M_t.w_t;
    paths.r_net(t) = M_t.r_net_t;
    paths.L(t) = M_t.L_t;
    paths.Y(t) = M_t.Y_t;
    
    % [新增] f. 计算投资和检验商品市场出清
    K_total_t = K_pvt_t + B_p_t + K_pps_t;
    K_total_t_plus_1 = paths.K_pvt(t+1) + paths.B_p(t+1) + paths.K_pps(t+1);
    paths.Investment(t) = K_total_t_plus_1 - (1 - cS.ddk) * K_total_t;
    paths.market_clearing_residual(t) = paths.Y(t) - (paths.C(t) + paths.Investment(t) + paths.G(t));
    
    % 商品市场出清检验: Y_t = C_t + I_t + G_t
    market_clearing_residual = paths.Y(t) - (paths.C(t) + paths.Investment(t) + paths.G(t));
    
    % 如果残差过大，发出警告
    if abs(market_clearing_residual) > 0.01 * paths.Y(t)  % 1%的容忍度
        fprintf('  ⚠️  t=%d: 商品市场残差 = %.4f (%.2f%% of GDP)\n', ...
            t, market_clearing_residual, market_clearing_residual/paths.Y(t)*100);
    end
    
    % g. 实时监控
    % if mod(t, 10) == 0 || t == T
        fprintf('%6d | %12.3f | %12.3f | %12.3f | %12.3f | %10.3f | %10.3f | %11.2f%% | %10.3f\n', ...
            t, paths.K_pvt(t), paths.K_pps(t), paths.B_p(t), paths.B_g(t), paths.Y(t), paths.w(t), paths.r_net(t)*100, paths.G(t));
    % end
end

fprintf('%s\n', repmat('-', 1, 130));
fprintf('✅ 过渡动态模拟完成。\n');

% 商品市场出清统计
max_residual = max(abs(paths.market_clearing_residual));
mean_residual = mean(abs(paths.market_clearing_residual));
fprintf('\n--- 商品市场出清检验结果 ---\n');
fprintf('最大残差: %.6f (%.3f%% of 平均GDP)\n', max_residual, max_residual/mean(paths.Y(1:T))*100);
fprintf('平均残差: %.6f (%.3f%% of 平均GDP)\n', mean_residual, mean_residual/mean(paths.Y(1:T))*100);
if max_residual < 0.001 * mean(paths.Y(1:T))
    fprintf('✅ 商品市场近似出清，残差在可接受范围内。\n');
else
    fprintf('⚠️ 商品市场残差较大，可能存在数值问题。\n');
end

%% 5. 可视化结果
fprintf('\n--- 5. 绘制过渡动态路径 ---\n');
figure('Name', 'OLG V10: 过渡动态分析', 'Position', [100, 100, 1600, 1000]);
time_axis = 0:T;

subplot(3,3,1);
plot(time_axis, paths.K_pvt, 'b-', 'LineWidth', 2);
title('私人资本 K_{pvt}'); xlabel('时期 (t)'); grid on;

subplot(3,3,2);
plot(time_axis, paths.K_pps, 'm-', 'LineWidth', 2);
title('个人养老金 K_{pps}'); xlabel('时期 (t)'); grid on;

subplot(3,3,3);
plot(time_axis, paths.B_p, 'r-', 'LineWidth', 2);
title('养老金基金 B_{p}'); xlabel('时期 (t)'); grid on;

subplot(3,3,4);
plot(time_axis, paths.B_g, 'g-', 'LineWidth', 2);
title('政府债务 B_{g}'); xlabel('时期 (t)'); grid on;

subplot(3,3,5);
plot(time_axis(1:T+1), [paths.K_pvt + paths.K_pps + paths.B_p], 'k-', 'LineWidth', 2);
title('总资本 K_{total}'); xlabel('时期 (t)'); grid on;

subplot(3,3,6);
plot(time_axis(1:T), paths.Investment, 'c-', 'LineWidth', 2);
title('总投资 I'); xlabel('时期 (t)'); grid on;

subplot(3,3,7);
plot(time_axis(1:T), paths.Y(1:T), 'y-', 'LineWidth', 2);
title('总产出 Y'); xlabel('时期 (t)'); grid on;

subplot(3,3,8);
plot(time_axis(1:T), paths.market_clearing_residual, 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980]);
title('商品市场残差 Y-C-I-G'); xlabel('时期 (t)'); grid on;
ylabel('绝对残差');

subplot(3,3,9);
plot(time_axis, paths.G, 'LineWidth', 2, 'Color', [0.5 0.5 0.5]);
title('政府支出 G'); xlabel('时期 (t)'); grid on;

sgtitle('OLG 模型 V10 过渡动态分析', 'FontSize', 16, 'FontWeight', 'bold');