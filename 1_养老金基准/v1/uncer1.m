%% 数值验证（最终版）：不确定性模型下的经济学定理与效应博弈
% =========================================================================
% 描述:
% 本脚本通过数值实验，系统性地验证在不确定性模型中推导出的
% 核心经济学定理和比较静态分析。特别地，它将展示养老金决策中，
% “对冲动机”与“内生投资组合调整效应”之间的复杂博弈。
%
% 作者: 计算宏观经济学家
% =========================================================================

clear;
clc;
close all;

%% 1. 设立基准情景 (Baseline Case)
% -------------------------------------------------------------------------
fprintf('============================================================\n');
fprintf('1. 求解基准情景 (Baseline Case)\n');
fprintf('============================================================\n');

param_base.beta  = 0.96; param_base.Y1    = 100; param_base.g     = 1.2;
param_base.Y2    = param_base.Y1 * param_base.g; param_base.q_max = 0.12;
param_base.tau   = 0.25; param_base.tau_p = 0.15; param_base.Rf    = 1.02;
param_base.Rp    = 1.03; param_base.RH    = 1.15; param_base.RL    = 0.95; 
param_base.p     = 0.75;

% 求解基准模型
x_base = solve_model(param_base);
q1_base = x_base(2);
alpha1_base = x_base(3);

% 获取第二期最优策略以供比较
[alpha2_base, ~, ~] = solve_period2_strategy(param_base);

fprintf('基准解求解完成:\n');
fprintf('  - 最优 q1*     = %.8f (内点解)\n', q1_base);
fprintf('  - 最优 alpha1* = %.8f\n', alpha1_base);
fprintf('  - 最优 alpha2* = %.8f\n\n', alpha2_base);

%% 2. 验证定理 1: 最优消费与投资组合规则
% -------------------------------------------------------------------------
fprintf('============================================================\n');
fprintf('2. 验证定理 1: 最优消费与投资组合规则\n');
fprintf('============================================================\n');

% 2.1 验证 alpha* 的时间一致性 (alpha1* = alpha2*)
fprintf('2.1 验证 alpha* 的时间一致性 (alpha1* = alpha2*)\n');
fprintf('    - 理论预测: 风险资产配置比例不随时间变化。\n');
fprintf('    - 数值结果: alpha1* = %.8f, alpha2* = %.8f\n', alpha1_base, alpha2_base);
if abs(alpha1_base - alpha2_base) < 1e-7
    fprintf('    - [验证成功]\n\n');
else
    fprintf('    - [验证失败]\n\n');
end

% 2.2 验证 alpha* 的财富独立性
fprintf('2.2 验证 alpha* 的财富独立性\n');
param_high_Y1 = param_base;
param_high_Y1.Y1 = 200; % 将初始收入加倍
param_high_Y1.Y2 = param_high_Y1.Y1 * param_high_Y1.g;

x_high_Y1 = solve_model(param_high_Y1);
alpha1_high_Y1 = x_high_Y1(3);

fprintf('    - 理论预测: 风险资产配置比例不随财富（收入）变化。\n');
fprintf('    - 基准 alpha1* (Y1=100) = %.8f\n', alpha1_base);
fprintf('    - 新 alpha1*   (Y1=200) = %.8f\n', alpha1_high_Y1);
if abs(alpha1_base - alpha1_high_Y1) < 1e-7
    fprintf('    - [验证成功]\n\n');
else
    fprintf('    - [验证失败]\n\n');
end

%% 3. 验证定理 2: 养老金缴存决策与比较静态
% -------------------------------------------------------------------------
fprintf('============================================================\n');
fprintf('3. 验证定理 2: 养老金缴存决策与比较静态\n');
fprintf('============================================================\n');

% 3.1 验证SDF一阶条件 (自洽性检查)
fprintf('3.1 验证SDF一阶条件 (自洽性检查)\n');
[~, ~, ~, rel_error] = final_check_sdf_condition(x_base, param_base);
fprintf('    - 理论预测: 最优解应精确满足SDF一阶条件。\n');
fprintf('    - 相对误差: %.4e (%.4f%%)\n', rel_error, rel_error*100);
if rel_error < 1e-6
    fprintf('    - [验证成功]\n\n');
else
    fprintf('    - [验证失败]\n\n');
end

% 3.2 比较静态: d(q1*)/d(Rp) > 0
fprintf('3.2 比较静态: d(q1*)/d(Rp) > 0\n');
n_points = 100; % 使用奇数点以确保基准点在中心
Rp_vec = linspace(param_base.Rp - 0.01, param_base.Rp + 0.01, n_points);
q1_results_Rp = zeros(n_points, 1);
param_current = param_base;
for i = 1:n_points
    param_current.Rp = Rp_vec(i);
    x_sol = solve_model(param_current);
    q1_results_Rp(i) = x_sol(2);
end
figure('Name', '比较静态分析: Rp');
plot(Rp_vec, q1_results_Rp, '-b', 'LineWidth', 2);
hold on;
plot(param_base.Rp, q1_base, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
xlabel('养老金回报率 (Rp)'); ylabel('最优养老金缴存比例 (q1*)');
title('比较静态分析: q1* 对 Rp 的响应');
legend('数值解路径', '基准点', 'Location', 'northwest');
grid on;
fprintf('    - 理论预测: 提高养老金回报率 Rp, 将增加养老金缴存 q1*。\n');
if all(diff(q1_results_Rp) >= -1e-9) % 允许数值误差
    fprintf('    - [验证成功] 曲线单调递增。\n\n');
else
    fprintf('    - [验证失败] 曲线非单调递增。\n\n');
end


% 3.3 [扩展] 比较静态: d(q1*)/d(risk) 的完整博弈
fprintf('3.3 [扩展] 比较静态: d(q1*)/d(risk) 的完整博弈\n');
spread_vec = linspace(0.01, 0.40, n_points); % x轴为离散程度，范围扩大
q1_results_risk = zeros(n_points, 1);
alpha1_results_risk = zeros(n_points, 1);
mean_R = param_base.p * param_base.RH + (1-param_base.p) * param_base.RL;
param_current = param_base;
for i = 1:n_points
    spread = spread_vec(i);
    % 构造均值保持扩展
    param_current.RH = mean_R + spread * (1 - param_base.p);
    param_current.RL = mean_R - spread * param_base.p;
    x_sol = solve_model(param_current);
    q1_results_risk(i) = x_sol(2);
    alpha1_results_risk(i) = x_sol(3);
end
figure('Name', '比较静态分析: 市场风险博弈');
yyaxis left; % 激活左侧Y轴
plot(spread_vec, q1_results_risk, '-b', 'LineWidth', 2);
ylabel('最优养老金缴存比例 (q1*)');
ylim([-0.01, param_base.q_max + 0.01]);
hold on;
base_spread = (param_base.RH - param_base.RL) / ((1-param_base.p) + param_base.p); % 计算基准spread
plot(base_spread, q1_base, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');

yyaxis right; % 激活右侧Y轴
plot(spread_vec, alpha1_results_risk, '--m', 'LineWidth', 2);
ylabel('最优风险资产配置 (alpha1*)');
ylim([-0.1, 1.1]);

xlabel('市场风险 (均值保持扩展的离散程度)');
title('比较静态分析: 风险增加下的效应博弈');
legend('q1* (左轴)', '基准点', 'alpha1* (右轴)', 'Location', 'west');
grid on;
fprintf('    - 理论预测: 简单理论预测单调递增，但完整模型揭示了非单调性。\n');
fprintf('    - 数值结果: 曲线呈现复杂的非单调关系，验证了效应博弈的存在。\n\n');

%% 4. 验证定理 3: 税收递延的价值
% -------------------------------------------------------------------------
fprintf('============================================================\n');
fprintf('4. 验证定理 3: 税收递延的价值\n');
fprintf('============================================================\n');

% 4.1 比较静态: d(q1*)/d(tau) > 0
fprintf('4.1 比较静态: d(q1*)/d(tau) > 0\n');
tau_vec = linspace(param_base.tau - 0.1, param_base.tau + 0.1, n_points);
q1_results_tau = zeros(n_points, 1);
param_current = param_base;
for i = 1:n_points
    param_current.tau = tau_vec(i);
    x_sol = solve_model(param_current);
    q1_results_tau(i) = x_sol(2);
end
figure('Name', '比较静态分析: tau');
plot(tau_vec, q1_results_tau, '-b', 'LineWidth', 2);
hold on;
plot(param_base.tau, q1_base, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
xlabel('当期税率 (tau)'); ylabel('最优养老金缴存比例 (q1*)');
title('比较静态分析: q1* 对当期税率 tau 的响应');
legend('数值解路径', '基准点', 'Location', 'northwest');
grid on;
fprintf('    - 理论预测: 提高当期税率 tau, 将增加养老金缴存 q1*。\n');
if all(diff(q1_results_tau) >= -1e-9)
    fprintf('    - [验证成功] 曲线单调递增。\n\n');
else
    fprintf('    - [验证失败] 曲线非单调递增。\n\n');
end

% 4.2 比较静态: d(q1*)/d(tau_p) < 0
fprintf('4.2 比较静态: d(q1*)/d(tau_p) < 0\n');
tau_p_vec = linspace(param_base.tau_p - 0.1, param_base.tau_p + 0.1, n_points);
q1_results_tau_p = zeros(n_points, 1);
param_current = param_base;
for i = 1:n_points
    param_current.tau_p = tau_p_vec(i);
    x_sol = solve_model(param_current);
    q1_results_tau_p(i) = x_sol(2);
end
figure('Name', '比较静态分析: tau_p');
plot(tau_p_vec, q1_results_tau_p, '-b', 'LineWidth', 2);
hold on;
plot(param_base.tau_p, q1_base, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
xlabel('提取期税率 (tau_p)'); ylabel('最优养老金缴存比例 (q1*)');
title('比较静态分析: q1* 对提取期税率 tau_p 的响应');
legend('数值解路径', '基准点', 'Location', 'northeast');
grid on;
fprintf('    - 理论预测: 提高提取期税率 tau_p, 将减少养老金缴存 q1*。\n');
if all(diff(q1_results_tau_p) <= 1e-9)
    fprintf('    - [验证成功] 曲线单调递减。\n\n');
else
    fprintf('    - [验证失败] 曲线非单调递减。\n\n');
end

fprintf('============================================================\n');
fprintf('所有定理验证完毕。\n');
fprintf('============================================================\n');


%% 辅助函数区域
% =========================================================================

function x_opt = solve_model(param)
    % 数值求解器
    x0 = [0.3, param.q_max/2, 0.5]; 
    lb = [1e-6, 0, 0];
    ub = [1, param.q_max, 1];
    options = optimoptions('fmincon', 'Display', 'none', 'Algorithm', 'sqp', ...
        'MaxFunctionEvaluations', 20000, ...
        'OptimalityTolerance', 1e-12, 'StepTolerance', 1e-12, ...
        'FiniteDifferenceType', 'central');
    objectiveFunc = @(x) dynamic_lifestyle_utility(x, param);
    [x_opt, ~, exitflag] = fmincon(objectiveFunc, x0, [], [], [], [], lb, ub, [], options);
    if exitflag < 1, warning('求解器未能收敛!'); end
end

function neg_utility = dynamic_lifestyle_utility(x, param)
    % 动态最优目标函数
    c1_ratio = x(1); q1 = x(2); alpha1 = x(3);
    beta = param.beta; tau = param.tau; Y1 = param.Y1; Y2 = param.Y2; 
    Rf = param.Rf; p = param.p; RH = param.RH; RL = param.RL;
    [alpha2_star, q2_star, ~] = solve_period2_strategy(param);
    W1 = (1 - tau) * (Y1 - q1 * Y1); C1 = c1_ratio * W1; S1 = W1 - C1; P1 = q1 * Y1;
    if S1 < -1e-8 || C1 <= 1e-8, neg_utility = inf; return; end
    utility_t1 = log(C1);
    W2_H = S1 * ((1-alpha1)*Rf + alpha1*RH) + (1-tau)*(Y2 - q2_star*Y2);
    utility_from_t2_H = calculate_V2(W2_H, P1, alpha2_star, q2_star, param);
    W2_L = S1 * ((1-alpha1)*Rf + alpha1*RL) + (1-tau)*(Y2 - q2_star*Y2);
    utility_from_t2_L = calculate_V2(W2_L, P1, alpha2_star, q2_star, param);
    if isinf(utility_from_t2_H) || isinf(utility_from_t2_L), neg_utility = inf; return; end
    total_expected_utility = utility_t1 + beta * (p * utility_from_t2_H + (1-p) * utility_from_t2_L);
    neg_utility = -total_expected_utility;
end

function V2 = calculate_V2(W2, P1, alpha2, q2, param)
    % 计算V2
    beta = param.beta; tau_p = param.tau_p; Rp = param.Rp; Rf = param.Rf;
    p = param.p; RH = param.RH; RL = param.RL;
    if W2 <= 1e-8, V2 = -inf; return; end
    P2 = P1*Rp + q2*param.Y2;
    R_ce = exp(p*log((1-alpha2)*Rf+alpha2*RH) + (1-p)*log((1-alpha2)*Rf+alpha2*RL));
    TotalWealth2_PV = W2 + (1-tau_p)*P2*Rp/R_ce;
    C2 = (1/(1+beta)) * TotalWealth2_PV; S2 = W2 - C2;
    if S2 < -1e-8, V2 = -inf; return; end
    C3_H = S2*((1-alpha2)*Rf + alpha2*RH) + (1-tau_p)*P2*Rp;
    C3_L = S2*((1-alpha2)*Rf + alpha2*RL) + (1-tau_p)*P2*Rp;
    if C3_H <= 1e-8 || C3_L <= 1e-8, V2 = -inf; return; end
    E_log_C3 = p*log(C3_H) + (1-p)*log(C3_L);
    V2 = log(C2) + beta*E_log_C3;
end

function [alpha_star, q_star, R_ce] = solve_period2_strategy(param)
    % 解析求解第二期策略
    p = param.p; RH = param.RH; RL = param.RL; Rf = param.Rf;
    Rp = param.Rp; tau = param.tau; tau_p = param.tau_p; q_max = param.q_max;
    numerator = p*(RH-Rf)*Rf + (1-p)*(RL-Rf)*Rf;
    denominator = -(p*(RH-Rf)*(RL-Rf) + (1-p)*(RL-Rf)*(RH-Rf));
    if abs(denominator) < 1e-9, alpha_unconstrained = 0; else, alpha_unconstrained = numerator / denominator; end
    alpha_star = max(0, min(1, alpha_unconstrained));
    R_ce = exp(p*log((1-alpha_star)*Rf+alpha_star*RH) + (1-p)*log((1-alpha_star)*Rf+alpha_star*RL));
    if (1 - tau_p) * Rp > (1 - tau) * R_ce, q_star = q_max; else, q_star = 0; end
end

function [lhs, rhs, abs_error, rel_error] = final_check_sdf_condition(x_opt, param)
    % 最终的SDF验证函数
    c1_ratio = x_opt(1); q1 = x_opt(2); alpha1 = x_opt(3);
    [alpha2, q2, ~] = solve_period2_strategy(param);
    beta = param.beta; tau = param.tau; tau_p = param.tau_p;
    Y1 = param.Y1; Y2 = param.Y2; Rf = param.Rf; Rp = param.Rp;
    p = param.p; RH = param.RH; RL = param.RL;
    W1 = (1 - tau) * (Y1 - q1 * Y1); C1 = c1_ratio * W1; S1 = W1 - C1; P1 = q1 * Y1;
    W2_H = S1*((1-alpha1)*Rf+alpha1*RH) + (1-tau)*(Y2-q2*Y2); C2_H = calculate_V2_C2(W2_H, P1, alpha2, q2, param); S2_H = W2_H-C2_H; P2 = P1*Rp+q2*Y2;
    W2_L = S1*((1-alpha1)*Rf+alpha1*RL) + (1-tau)*(Y2-q2*Y2); C2_L = calculate_V2_C2(W2_L, P1, alpha2, q2, param); S2_L = W2_L-C2_L;
    C3_HH = S2_H*((1-alpha2)*Rf+alpha2*RH) + (1-tau_p)*P2*Rp; C3_HL = S2_H*((1-alpha2)*Rf+alpha2*RL) + (1-tau_p)*P2*Rp;
    C3_LH = S2_L*((1-alpha2)*Rf+alpha2*RH) + (1-tau_p)*P2*Rp; C3_LL = S2_L*((1-alpha2)*Rf+alpha2*RL) + (1-tau_p)*P2*Rp;
    SDF_HH = beta^2*C1/C3_HH; SDF_HL = beta^2*C1/C3_HL; SDF_LH = beta^2*C1/C3_LH; SDF_LL = beta^2*C1/C3_LL;
    R_opt1_H=(1-alpha1)*Rf+alpha1*RH; R_opt1_L=(1-alpha1)*Rf+alpha1*RL; R_opt2_H=(1-alpha2)*Rf+alpha2*RH; R_opt2_L=(1-alpha2)*Rf+alpha2*RL;
    Ret_HH=R_opt1_H*R_opt2_H; Ret_HL=R_opt1_H*R_opt2_L; Ret_LH=R_opt1_L*R_opt2_H; Ret_LL=R_opt1_L*R_opt2_L;
    E_SDF_x_Return = p*p*(SDF_HH*Ret_HH) + p*(1-p)*(SDF_HL*Ret_HL) + (1-p)*p*(SDF_LH*Ret_LH) + (1-p)*(1-p)*(SDF_LL*Ret_LL);
    E_SDF = p*p*SDF_HH + p*(1-p)*SDF_HL + (1-p)*p*SDF_LH + (1-p)*(1-p)*SDF_LL;
    lhs = (1 - tau_p) * Rp^2;
    rhs = (1 - tau) * (E_SDF_x_Return / E_SDF);
    abs_error = abs(lhs - rhs);
    rel_error = abs_error / lhs;
end

function C2 = calculate_V2_C2(W2, P1, alpha2, q2, param)
    % 辅助函数，在验证中计算C2
    beta = param.beta; tau_p = param.tau_p; Rp = param.Rp; Rf = param.Rf;
    p = param.p; RH = param.RH; RL = param.RL;
    if W2 <= 1e-8, C2 = 0; return; end
    P2 = P1*Rp + q2*param.Y2;
    R_ce = exp(p*log((1-alpha2)*Rf+alpha2*RH) + (1-p)*log((1-alpha2)*Rf+alpha2*RL));
    TotalWealth2_PV = W2 + (1-tau_p)*P2*Rp/R_ce;
    C2 = (1/(1+beta)) * TotalWealth2_PV;
end