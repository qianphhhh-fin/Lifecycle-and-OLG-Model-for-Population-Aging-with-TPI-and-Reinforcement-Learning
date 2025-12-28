%% 验证三期生命周期模型(EET)的解析解与数值解的一致性
% =========================================================================
% 描述:
% 本脚本通过一个蒙特卡洛实验来验证EET税收模式下三期模型的解析解。
% 它会进行 n 次模拟，每次都随机生成一组经济参数。
% 在每次模拟中，它会分别用解析公式和数值优化(fmincon)来求解模型，
% 然后比对结果是否在容差范围内完全一致。
%
% 预期结果:
% 如果解析解是正确的，匹配率应为 100%。
% =========================================================================

clear;
clc;
close all;

%% 1. 模拟设置
n = 100; % 设置随机参数抽取的次数
match_count = 0; % 初始化匹配成功计数器
tolerance = 1e-5; % 设置判断相等的容差

fprintf('开始验证... 将进行 %d 次随机参数模拟。\n', n);
fprintf('============================================================\n');

%% 2. 主循环
for i = 1:n
    fprintf('\n--- 模拟运行: %d/%d ---\n', i, n);
    
    % --- 2.1 随机生成参数 ---
    param.beta  = 0.96;
    param.Y1    = 100;
    param.g     = 1.2;
    param.Y2    = param.Y1 * param.g;
    param.q_max = 0.12;
    
    % 在合理范围内随机抽取税率和收益率
    param.tau   = 0.10 + 0.30 * rand; % 税率在 [0.10, 0.40] 之间
    param.tau_p = 0.05 + 0.20 * rand; % 提取税率在 [0.05, 0.25] 之间
    param.Rf    = 1.01 + 0.05 * rand; % 无风险利率在 [1.01, 1.06] 之间
    param.Rp    = 1.01 + 0.07 * rand; % 养老金收益率在 [1.01, 1.08] 之间

    fprintf('当前参数: tau=%.3f, tau_p=%.3f, Rf=%.3f, Rp=%.3f\n', ...
            param.tau, param.tau_p, param.Rf, param.Rp);
            
    % --- 2.2 解析求解 ---
    sol_analytical = solve_analytical(param);
    
    % --- 2.3 数值求解 ---
    % 决策变量向量 x = [c1, q1, c2, q2]
    x0 = [0.3, param.q_max/2, 0.5, param.q_max/2]; 
    lb = [1e-6, 0, 1e-6, 0];
    ub = [1, param.q_max, 1, param.q_max];
    options = optimoptions('fmincon', 'Display', 'none', 'Algorithm', 'sqp');
    objectiveFunc = @(x) lifestyle_utility(x, param);
    
    [x_opt, ~, exitflag] = fmincon(objectiveFunc, x0, [], [], [], [], lb, ub, [], options);
    
    if exitflag < 1
        fprintf('!! 数值求解器未能成功收敛 (exitflag=%d)，跳过本次验证 !!\n', exitflag);
        continue;
    end
    
    sol_numerical.c1 = x_opt(1);
    sol_numerical.q1 = x_opt(2);
    sol_numerical.c2 = x_opt(3);
    sol_numerical.q2 = x_opt(4);
    
    % --- 2.4 结果比对 ---
    err_c1 = abs(sol_analytical.c1 - sol_numerical.c1);
    err_q1 = abs(sol_analytical.q1 - sol_numerical.q1);
    err_c2 = abs(sol_analytical.c2 - sol_numerical.c2);
    err_q2 = abs(sol_analytical.q2 - sol_numerical.q2);
    
    if err_c1 < tolerance && err_q1 < tolerance && err_c2 < tolerance && err_q2 < tolerance
        fprintf('结果: [成功] 解析解与数值解完全匹配。\n');
        match_count = match_count + 1;
    else
        fprintf('结果: [失败] 解析解与数值解不匹配！\n');
        fprintf('%15s | %15s | %15s\n', '变量', '解析解', '数值解');
        fprintf('--------------------------------------------------\n');
        fprintf('%15s | %15.8f | %15.8f\n', 'c1*', sol_analytical.c1, sol_numerical.c1);
        fprintf('%15s | %15.8f | %15.8f\n', 'q1*', sol_analytical.q1, sol_numerical.q1);
        fprintf('%15s | %15.8f | %15.8f\n', 'c2*', sol_analytical.c2, sol_numerical.c2);
        fprintf('%15s | %15.8f | %15.8f\n', 'q2*', sol_analytical.q2, sol_numerical.q2);
    end
end

%% 3. 最终总结
fprintf('\n============================================================\n');
fprintf('验证完成。\n');
fprintf('在 %d 次随机模拟中，有 %d 次解析解与数值解完全一致。\n', n, match_count);
fprintf('匹配率: %.1f%%\n', (match_count/n)*100);
fprintf('============================================================\n');


%% 辅助函数区域

function sol = solve_analytical(param)
    % =====================================================================
    % 使用解析公式求解模型
    % =====================================================================
    
    % 解包参数
    beta  = param.beta; tau = param.tau; tau_p = param.tau_p;
    Y1 = param.Y1; Y2 = param.Y2; q_max = param.q_max;
    Rf = param.Rf; Rp = param.Rp;

    % 1. 决定最优缴存比例 q1* 和 q2* (角点解)
    if (1 - tau_p) * Rp^2 > (1 - tau) * Rf^2
        q1_star = q_max;
    else
        q1_star = 0;
    end
    
    if (1 - tau_p) * Rp > (1 - tau) * Rf
        q2_star = q_max;
    else
        q2_star = 0;
    end
    
    % 2. 计算生命周期总财富现值 A1*
    term_Y1 = (1 - tau) * (Y1 - q1_star * Y1);
    term_Y2 = (1 - tau) * (Y2 - q2_star * Y2) / Rf;
    term_P  = (1 - tau_p) * (q1_star * Y1 * Rp + q2_star * Y2) * Rp / Rf^2;
    A1_star = term_Y1 + term_Y2 + term_P;
    
    % 3. 计算最优消费 C1* 和消费比例 c1*
    C1_star = (1 / (1 + beta + beta^2)) * A1_star;
    W1_star = (1 - tau) * (Y1 - q1_star * Y1);
    c1_star = C1_star / W1_star;
    
    % 4. 计算第二期的状态变量 W2* 和 P2*
    P1_star = q1_star * Y1;
    W2_star = (W1_star - C1_star) * Rf + (1 - tau) * (Y2 - q2_star * Y2);
    P2_star = P1_star * Rp + q2_star * Y2;
    
    % 5. 计算最优消费 C2* 和消费比例 c2*
    C2_star = (1 / (1 + beta)) * (W2_star + (1 - tau_p) * P2_star * Rp / Rf);
    c2_star = C2_star / W2_star;
    
    % 封装结果
    sol.c1 = c1_star;
    sol.q1 = q1_star;
    sol.c2 = c2_star;
    sol.q2 = q2_star;
end


function [neg_utility, C, W, P] = lifestyle_utility(x, param)
    % =====================================================================
    % 目标函数: 计算给定决策变量下的负总效用 (EET税收模式)
    % =====================================================================

    % 解包决策变量
    c1 = x(1); q1 = x(2); c2 = x(3); q2 = x(4);
    
    % 解包参数
    beta  = param.beta; tau   = param.tau; tau_p = param.tau_p;
    Y1    = param.Y1; Y2    = param.Y2; Rf    = param.Rf; Rp    = param.Rp;

    % 时期 1: 工作期 (EET 模式)
    W1 = (1 - tau) * (Y1 - q1 * Y1); 
    C1 = c1 * W1;
    P1 = q1 * Y1;
    if W1 - C1 < -1e-6, neg_utility = inf; return; end % 允许微小负储蓄以处理数值误差
    
    % 时期 2: 过渡期 (EET 模式)
    W2 = (W1 - C1) * Rf + (1 - tau) * (Y2 - q2 * Y2);
    C2 = c2 * W2;
    P2 = P1 * Rp + q2 * Y2;
    if W2 - C2 < -1e-6, neg_utility = inf; return; end
    
    % 时期 3: 退休期
    W3 = (W2 - C2) * Rf + (1 - tau_p) * P2 * Rp;
    C3 = W3;
    
    % 计算效用
    if C1 <= 0 || C2 <= 0 || C3 <= 0
        neg_utility = inf;
    else
        utility = log(C1) + beta * log(C2) + beta^2 * log(C3);
        neg_utility = -utility;
    end
end