% --- START OF FILE main_olg_v12.m ---
clear; close all;
addpath(pwd);
fprintf('=== OLG 模型 V12: 最终会计修正版 ===\n');

%% 1. 初始化参数和环境
fprintf('\n--- 1. 初始化参数 ---\n');
cS = main_olg_v12_utils.ParameterValues_HuggettStyle();
paramS = struct();
% ... (所有参数设置保持不变) ...
cS.sigma = 3.0; cS.beta = 0.97; cS.cFloor = 0.05; cS.nSim = 1000;
cS.phi_bequest = 3; cS.sigma_bequest = cS.sigma;
cS.start_year = 1997; cS.end_year = 2102; cS.time_step = 5;
cS.alpha = 0.45; cS.ddk = 1 - (1 - 0.08)^cS.time_step;
cS.tau_k = 0.25; cS.tau_l = 0.08; cS.tau_c = 0.04;
% cS.gov_exp_frac_Y = 0.15; % G将内生等于T
ngrid = 40; cS.nk = ngrid; cS.nkpps = ngrid; cS.nkprime = ngrid; cS.npps = ngrid;
cS.A = 1.0; cS = main_olg_v12_utils.generateGrids(cS);
cS.rho_urban_employee = 0.60; cS.rho_resident = 0.15;
cS.pps_active = false; cS.pps_activation_year = 2022;
cS.pps_tax_rate_withdrawal = 0.03; cS.pps_return_rate_premium = 0.01;
cS.pps_withdrawal_rate = 0.15; cS.pps_contrib_limit = 9999; cS.pps_max_contrib_frac = 0.1;
cS = main_olg_v12_utils.calcaulte_theta_payg_path(cS, true);
fprintf('✅ 复合有效缴费率与覆盖率路径构建完成。\n');

%% 2. 加载外生数据路径
[Z_path, A_path, T] = main_olg_v12_utils.load_exogenous_paths(cS);
sim_years = cS.start_year:cS.time_step:(cS.start_year + (T-1)*cS.time_step);
[paramS.leLogGridV, paramS.leTrProbM, paramS.leProb1V] = main_olg_v12_utils.EarningProcess_olgm(cS);
paramS.leGridV = exp(paramS.leLogGridV(:));
paramS.ageEffV_new = cS.ageEffV_new;
eIdxM = main_olg_v12_utils.LaborEndowSimulation_olgm_AgeGroup(cS, paramS);
fprintf('✅ 效率冲击路径生成完成。\n');

%% 3. 求解start_year年 (t=1) 的初始稳态
fprintf('\n--- 3. 求解start_year年 (t=1) 的初始稳态 ---\n');
Z_1 = Z_path(:, 1); A_1 = A_path(1);
[ss_1997, eq_found] = main_olg_v12_utils.solve_initial_steady_state(Z_1, A_1, cS, paramS, eIdxM);
if ~eq_found, error('无法找到初始稳态，程序终止。'); end

% 初始化所有状态变量的路径
paths = struct();
paths.K_pvt = zeros(1, T); paths.K_pps = zeros(1, T); paths.B_p = zeros(1, T);
paths.w = zeros(1, T); paths.r_net = zeros(1, T); paths.L = zeros(1, T);
paths.Y = zeros(1, T); paths.G = zeros(1, T); paths.C = zeros(1, T);
paths.Investment = zeros(1, T); paths.market_clearing_residual = zeros(1, T);
paths.K_prod = zeros(1, T);

% --- 将求解出的稳态值赋给路径的 t=1 时刻 ---
paths.Y(1) = ss_1997.Y; paths.K_pvt(1) = ss_1997.K_pvt; paths.K_pps(1) = ss_1997.K_pps;
paths.B_p(1) = ss_1997.B_p; paths.L(1) = ss_1997.L; paths.w(1) = ss_1997.w;
paths.r_net(1) = ss_1997.r_net; paths.C(1) = ss_1997.C; paths.G(1) = ss_1997.G;
paths.K_prod(1) = ss_1997.K_pvt + ss_1997.K_pps + ss_1997.B_p;
paths.Investment(1) = cS.ddk * paths.K_prod(1); % 稳态下 I=d*K
paths.market_clearing_residual(1) = paths.Y(1) - (paths.C(1) + paths.Investment(1) + paths.G(1));

fprintf('✅ 初始稳态求解完成: Y=%.3f, K_prod=%.3f, C=%.3f, I=%.3f, G=%.3f\n', ...
    paths.Y(1), paths.K_prod(1), paths.C(1), paths.Investment(1), paths.G(1));
fprintf('   稳态残差: %.3e\n', paths.market_clearing_residual(1));

%% 4. 过渡动态主循环
fprintf('\n--- 4. 开始向前模拟过渡动态 (start_year - %d) ---\n', sim_years(end));
fprintf('\n%6s | %8s | %12s | %12s | %12s | %10s | %10s | %12s\n', ...
    'Time', 'Year', 'Private K', 'PPS K', 'Pension Fund', 'Output', 'Net Rate', 'Residual');
fprintf('%s\n', repmat('-', 1, 100));
fprintf('%6d | %8d | %12.3f | %12.3f | %12.3f | %10.3f | %9.2f%% | %12.3e\n', ...
    1, sim_years(1), paths.K_pvt(1), paths.K_pps(1), paths.B_p(1), paths.Y(1), paths.r_net(1)*100, paths.market_clearing_residual(1));

% 使用审计版函数来获得 t=1 时刻的精确总量，用于演化到 t=2
M_t1 = main_olg_v12_utils.get_prices_and_policy_at_t(1, paths.K_pvt(1), paths.B_p(1), paths.K_pps(1), Z_path(:,1), A_path(1), cS, paramS, eIdxM);
[K_pvt_next, K_pps_next, ~, ~] = main_olg_v12_utils.simulate_private_capital_forward_audited(M_t1, Z_path(:,1), cS, paramS, eIdxM);

for t = 2:T
    if sim_years(t) >= cS.pps_activation_year && ~cS.pps_active
        cS.pps_active = true;
        fprintf('--- 激活个人养老金(PPS)系统于 %d (t=%d) ---\n', sim_years(t), t);
    end

    paths.K_pvt(t) = K_pvt_next;
    paths.K_pps(t) = K_pps_next;
    
    % [注意] 养老金和政府债务的更新依赖于 t-1 的价格
    M_t_minus_1 = main_olg_v12_utils.get_prices_and_policy_at_t(t-1, paths.K_pvt(t-1), paths.B_p(t-1), paths.K_pps(t-1), Z_path(:,t-1), A_path(t-1), cS, paramS, eIdxM);
    paths.B_p(t) = main_olg_v12_utils.update_pension_fund(paths.B_p(t-1), M_t_minus_1, Z_path(:,t-1), cS);
    
    paths.K_prod(t) = paths.K_pvt(t) + paths.K_pps(t) + paths.B_p(t);
    
    M_t = main_olg_v12_utils.get_prices_and_policy_at_t(t, paths.K_pvt(t), paths.B_p(t), paths.K_pps(t), Z_path(:,t), A_path(t), cS, paramS, eIdxM);
    
    % [核心修正] 每次都调用审计版来获得精确的C和T
    [K_pvt_next, K_pps_next, C_t, Aggregates_t] = main_olg_v12_utils.simulate_private_capital_forward_audited(M_t, Z_path(:,t), cS, paramS, eIdxM);

    paths.w(t) = M_t.w_t; paths.r_net(t) = M_t.r_net_t;
    paths.L(t) = M_t.L_t; paths.Y(t) = M_t.Y_t; paths.C(t) = C_t;
    
    % [核心修正] G = T, I = K' - (1-d)K
    T_t = Aggregates_t.TotalLaborTax + Aggregates_t.TotalCapitalTax + Aggregates_t.TotalConsumptionTax;
    paths.G(t) = T_t;
    paths.Investment(t) = paths.K_prod(t) - (1 - cS.ddk) * paths.K_prod(t-1);
    
    paths.market_clearing_residual(t) = paths.Y(t) - (paths.C(t) + paths.Investment(t) + paths.G(t));

    fprintf('%6d | %8d | %12.3f | %12.3f | %12.3f | %10.3f | %9.2f%% | %12.3e\n', ...
        t, sim_years(t), paths.K_pvt(t), paths.K_pps(t), paths.B_p(t), paths.Y(t), paths.r_net(t)*100, paths.market_clearing_residual(t));
end

fprintf('%s\n', repmat('-', 1, 100));
fprintf('✅ 过渡动态模拟完成。\n');

max_abs_residual = max(abs(paths.market_clearing_residual(2:end)));
fprintf('\n--- 商品市场出清检验最终总结 (t>1) ---\n');
fprintf('最大绝对残差: %.3e\n', max_abs_residual);
if max_abs_residual < 1e-4, fprintf('✅ 商品市场近似出清，模型数值稳定性良好。\n');
else, fprintf('⚠️ 商品市场最终残差较大，请检查模型实现。\n'); end

%% 5. 可视化结果
% ... (绘图部分保持不变, 但现在G是内生的) ...
figure('Name', 'OLG V12: 最终会计修正版', 'Position', [100, 100, 1600, 1000]);
time_axis_years = sim_years;
subplot(3,3,1); plot(time_axis_years, paths.Y, 'k-', 'LineWidth', 2); title('总产出 Y'); xlabel('年份'); grid on;
subplot(3,3,2); plot(time_axis_years, paths.K_prod, 'b-', 'LineWidth', 2); title('生产性资本 K_{prod}'); xlabel('年份'); grid on;
subplot(3,3,3); plot(time_axis_years, paths.B_p, 'r-', 'LineWidth', 2); title('养老金基金 B_{p}'); xlabel('年份'); grid on;
subplot(3,3,4); plot(time_axis_years, paths.G, 'g-', 'LineWidth', 2); title('政府购买 G (内生=T)'); xlabel('年份'); grid on;
subplot(3,3,5); plot(time_axis_years, paths.K_pps, 'c-', 'LineWidth', 2); title('个人养老金 K_{pps}'); xlabel('年份'); grid on;
subplot(3,3,6); plot(time_axis_years, paths.C, 'm-', 'LineWidth', 2); title('总消费 C'); xlabel('年份'); grid on;
subplot(3,3,7); plot(time_axis_years, paths.Investment, 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980]); title('总投资 I'); xlabel('年份'); grid on;
subplot(3,3,8); plot(time_axis_years, paths.market_clearing_residual, 'LineWidth', 2, 'Color', [0.5 0.5 0.5]); title('商品市场残差 Y-C-I-G'); xlabel('年份'); grid on;
subplot(3,3,9); yyaxis left; plot(time_axis_years, paths.r_net*100, 'b-'); ylabel('净利率 (%)');
yyaxis right; plot(time_axis_years, paths.w, 'r-'); ylabel('工资');
title('价格路径'); xlabel('年份'); grid on; legend('r_{net}', 'w', 'Location', 'best');
sgtitle('OLG 模型 V12 过渡动态分析 (最终会计修正版)', 'FontSize', 16, 'FontWeight', 'bold');
% --- END OF FILE main_olg_v12.m ---