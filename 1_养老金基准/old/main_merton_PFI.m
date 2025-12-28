clear
% 根据Samuelson论文"LIFETIME PORTFOLIO SELECTION BY DYNAMIC STOCHASTIC PROGRAMMING"
% 有限期Merton模型求解与数值模拟 (使用后向求解法backward induction)
% 验证CRRA效用函数下的消费和风险投资策略
% 模型特点：
% - 有限T期问题
% - 没有劳动收入
% - 没有养老金购买决策
% - 允许杠杆（风险资产投资比例允许为负或超过1）
% - 决策变量表示为比例：
%   C - 消费比例（占当前财富）
%   A - 风险资产投资比例（相对总财富，不是剩余投资）
%
% 根据Samuelson论文验证：
% - CRRA效用下，最优投资组合比例与财富水平无关
% - 最优风险资产配置满足公式：
%   ∫[(1-w)(1+r)+wZ]^(γ-1)(Z-1-r)dP(Z) = 0
% - 消费决策规则采用C*_T-i = c_i * W_T-i的形式

% 创建保存结果的目录
if ~exist('result_merton_PFI', 'dir')
    mkdir('result_merton_PFI');
end

% 设置随机种子以获得可重复的结果
rng(42);

%% 模型参数设置
% 效用函数参数
beta = 0.95;    % 折现因子 (对应论文中的1/(1+ρ))
rho = 1/beta - 1; % 时间偏好率
gamma = 0.5;    % CRRA效用函数参数 (γ<1)
% U(C) = C^γ/γ, γ<1，γ≠0
% U(C) = log(C), γ=0 (特殊情况)

% 收益率参数
Rf = 1.02;      % 无风险资产收益率
Rs_mean = 1.06;  % 风险资产平均收益率
sigma_s = 0.18;  % 风险资产收益率的标准差

% 状态变量网格设置
nW = 500;       % 财富网格点数
nR = 5;         % 风险资产收益率随机冲击的网格点数
T = 20;         % 有限期问题的期数

%% 初始化状态变量网格和随机冲击
% 财富网格 (线性均匀网格)
w_min = 0.1;    % 最小财富水平
w_max = 50;     % 最大财富水平
log_w_min = log(w_min);
log_w_max = log(w_max);
log_wealth_grid = linspace(log_w_min, log_w_max, nW)';
gW  = exp(log_wealth_grid);  % 对数等距网格转换为财富网格
% 设定分布类型：二项分布('binomial')或截断正态分布('truncated_normal')
dist_type = 'binomial'; % 可选: 'binomial', 'truncated_normal'

if strcmp(dist_type, 'binomial')
    % 使用Samuelson论文中的二项分布收益率模型
    % Prob{Z=λ} = 1/2 = Prob{Z=λ^(-1)}, λ>1
    
    % 计算满足条件的λ: λ > 1 + r + sqrt(2r + r^2)
    r_annual = Rf - 1;  % 年化无风险利率
    lambda_min = 1 + r_annual + sqrt(2*r_annual + r_annual^2);
    lambda = max(1.4, lambda_min * 1.01); % 确保λ足够大，论文提到的例子使用1.4
    
    % 定义二项分布的收益率和概率
    R_grid = [lambda, 1/lambda];
    R_prob = [0.5, 0.5];
    
    % 验证期望收益率满足条件
    expected_return = 0.5*lambda + 0.5/lambda;
    fprintf('使用二项分布收益率模型:\n');
    fprintf('λ值: %.4f (最小要求: %.4f)\n', lambda, lambda_min);
    fprintf('风险资产期望收益率: %.4f (无风险收益率: %.4f)\n', expected_return, Rf);
    fprintf('风险溢价: %.4f\n', expected_return - Rf);
    
    % 计算在这个设定下的最大安全杠杆水平
    w_max_safe = Rf / (Rf - min(R_grid));
    fprintf('最大安全杠杆水平(确保正收益): %.4f\n', w_max_safe);
    
    % 重新计算参数
    mu = expected_return - Rf;  % 风险溢价
    sigma_effective = sqrt((lambda - 1/lambda)^2 * 0.25); % 二项分布的标准差
    
else
    % 使用截断正态分布，确保收益率为非负
    % 1. 首先生成正态分布
    % 2. 然后截断，确保所有收益率为正
    mean_return = Rs_mean - 1;  % 平均超额收益率
    
    % 使用截断正态分布代替标准正态分布
    [R_grid, R_prob] = truncated_normal_grid(mean_return, sigma_s, nR, 0, Inf);
    R_grid = R_grid + 1;  % 转换回总收益率
    
    % 计算期望收益率
    expected_return = sum(R_grid .* R_prob);
    fprintf('使用截断正态分布收益率模型:\n');
    fprintf('风险资产期望收益率: %.4f (无风险收益率: %.4f)\n', expected_return, Rf);
    fprintf('风险溢价: %.4f\n', expected_return - Rf);
    
    % 计算标准差
    variance = sum(((R_grid - expected_return).^2) .* R_prob);
    sigma_effective = sqrt(variance);
    
    % 计算在这个设定下的最大安全杠杆水平
    w_max_safe = Rf / (Rf - min(R_grid));
    fprintf('最大安全杠杆水平(确保正收益): %.4f\n', w_max_safe);
    
    % 风险溢价
    mu = expected_return - Rf;
end

%% 计算Samuelson论文中的理论解
fprintf('计算Samuelson论文中的理论解...\n');

% 计算理论最优风险资产配置比例
% 根据Samuelson论文，我们需要解方程：
% ∫[(1-w)(1+r)+wZ]^(γ-1)(Z-1-r)dP(Z) = 0
% 使用fmincon求解
w0 = 0.5; % 初始猜测值
samuelson_opt = optimoptions('fmincon', 'Display', 'off', 'OptimalityTolerance', 1e-8,'Algorithm','interior-point');
[w_samuelson, ~] = fmincon(@(w) samuelson_condition(w, R_grid, R_prob, Rf, gamma), w0, [], [], [], [], -w_max_safe*0.9, w_max_safe*0.9, [], samuelson_opt);

% 计算Merton公式的理论最优风险资产配置比例 (用于比较)
w_merton = mu / (gamma * sigma_effective^2);  % 连续模型的理论值

fprintf('Samuelson公式求解的最优风险资产配置比例: %.4f\n', w_samuelson);
fprintf('Merton公式的最优风险资产配置比例: %.4f\n', w_merton);

% 计算风险调整后的收益率 E[(1 + r_p)^γ]
one_plus_r_star_gamma = 0;
for i = 1:length(R_grid)
    portfolio_return = (1-w_samuelson)*(Rf) + w_samuelson*R_grid(i);
    
    % 检查投资组合收益是否为负
    if portfolio_return <= 0
        fprintf('警告：在某些状态下投资组合收益为负(%.4f)，这超出了模型适用范围\n', portfolio_return);
        fprintf('风险资产收益率: %.4f, 风险资产配置比例: %.4f\n', R_grid(i), w_samuelson);
        fprintf('重新调整w_samuelson以确保所有状态下的收益为正...\n');
        
        % 使用安全的w值
        w_samuelson = w_max_safe * 0.9;  % 使用最大安全杠杆的90%
        fprintf('调整后的风险资产配置比例: %.4f\n', w_samuelson);
        
        % 重新计算所有投资组合收益
        one_plus_r_star_gamma = 0;
        for j = 1:length(R_grid)
            portfolio_return_j = (1-w_samuelson)*(Rf) + w_samuelson*R_grid(j);
            one_plus_r_star_gamma = one_plus_r_star_gamma + (portfolio_return_j^gamma) * R_prob(j);
        end
        break;
    end
    
    one_plus_r_star_gamma = one_plus_r_star_gamma + (portfolio_return^gamma) * R_prob(i);
end

% 根据(1+r*)^γ计算r*
r_star = (one_plus_r_star_gamma)^(1/gamma) - 1;
fprintf('风险调整收益率r* = %.4f\n', r_star);

% 计算Samuelson论文中的   
a_1 = (beta * one_plus_r_star_gamma)^(1/(gamma-1));
fprintf('增长率因子a_1 = %.4f\n', a_1);

% 计算每期的最优消费比例
c_optimal = zeros(T, 1);
for t = 1:T
    if t == T
        c_optimal(t) = 1.0;  % 终期全部消费
    else
        % 根据Samuelson论文中的公式:
        % c_i = a_1^i / (1 + a_1 + a_1^2 + ... + a_1^i)
        if abs(a_1 - 1) < 1e-10  % a_1接近1的情况
            c_optimal(t) = 1 / (1 + (T-t));
        else  % a_1不等于1的一般情况
            c_optimal(t) = a_1^(T-t) * (a_1-1) / (a_1^(T-t+1) - 1);
        end
    end
end

fprintf('Samuelson论文理论解:\n');
fprintf('最优风险资产配置比例 = %.4f\n', w_samuelson);
fprintf('风险调整收益率E[(1+r_p)^γ] = %.4f\n', one_plus_r_star_gamma);
fprintf('风险调整收益率r* = %.4f\n', r_star);
fprintf('a_1 = %.4f\n', a_1);
fprintf('最终期消费比例 = 1.0000\n');
fprintf('倒数第二期消费比例 = %.4f\n', c_optimal(T-1));
fprintf('倒数第三期消费比例 = %.4f\n', c_optimal(T-2));

%% 初始化后向求解的数据结构
% 策略变量 - 按时间、财富网格存储
C_star = zeros(nW, T);  % 最优消费策略
A_star = zeros(nW, T);  % 最优风险资产配置策略
V_star = zeros(nW, T+1);  % 最优值函数，包括终期

% 初始化终期值函数
% 终期最优策略是全部消费
for i = 1:nW
    if gamma == 0  % Bernoulli对数效用
        V_star(i, T+1) = log(gW(i));
    else
        V_star(i, T+1) = gW(i)^gamma / gamma;
    end
end

%% 后向求解 - 从终期往前求解每一期的最优策略（使用策略函数迭代）
fprintf('开始后向求解过程(使用策略函数迭代)...\n');
tic;

% 初始化优化选项
options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp', ...
    'MaxFunctionEvaluations', 1000, 'OptimalityTolerance', 1e-20);

% 从T-1期开始向前求解
for t = T:-1:1
    fprintf('求解第%d期...\n', t);
    
    % 创建临时数组以便在parfor中累积结果
    C_star_t = zeros(nW, 1);
    A_star_t = zeros(nW, 1);
    V_star_t = zeros(nW, 1);
    
    % 初始化策略 - 使用理论最优策略作为起点
    for i = 1:nW
        if t == T
            C_star_t(i) = 1.0;  % 终期全部消费
        else
            C_star_t(i) = c_optimal(t);  % 理论最优消费比例
        end
        A_star_t(i) = w_samuelson;  % 理论最优风险资产配置
    end
    
    % 策略函数迭代
    max_policy_iter = 20;   % 最大策略迭代次数
    policy_tol = 1e-8;     % 策略收敛判断阈值
    
    for iter = 1:max_policy_iter
        % 1. 策略评估 - 固定当前策略，计算值函数
        for i = 1:nW
            current_wealth = gW(i);
            c_ratio = C_star_t(i);
            a_ratio = A_star_t(i);
            
            % 计算当前策略下的值函数
            V_star_t(i) = policy_evaluation(c_ratio, a_ratio, current_wealth, ...
                V_star(:, t+1), gW, R_grid, R_prob, Rf, gamma, beta);
        end
        
        % 2. 策略改进 - 固定值函数，找到最优策略
        old_C = C_star_t;
        old_A = A_star_t;
        
        % 对每个财富网格点更新最优策略
        for i = 1:nW
            current_wealth = gW(i);
            
            % 终期处理：固定消费比例为1.0，无需优化
            if t == T
                C_star_t(i) = 1.0;  % 终期全部消费
                A_star_t(i) = w_samuelson;  % 风险资产配置理论值(实际上终期已无投资意义)
                continue;  % 跳过优化
            end
            
            % 非终期：初始猜测值为当前策略
            x0 = [C_star_t(i), A_star_t(i)];
            
            % 优化约束
            lb = [0, 0];    % 消费比例下限为0
            ub = [1, 1];     % 消费比例上限为1
            
            % 求解最优策略
            value_function = @(x) -period_value_function(x, current_wealth, V_star(:, t+1), gW, R_grid, R_prob, Rf, gamma, beta);
            [x_opt, ~] = fmincon(value_function, x0, [], [], [], [], lb, ub, [], options);
            
            % 更新策略
            C_star_t(i) = x_opt(1);
            A_star_t(i) = x_opt(2);
        end
        
        % 检查策略是否收敛
        C_change = max(abs(C_star_t - old_C));
        A_change = max(abs(A_star_t - old_A));
        
        fprintf('  策略迭代 #%d - 策略变化: C=%.6f, A=%.6f\n', iter, C_change, A_change);
        
        if max(C_change, A_change) < policy_tol
            fprintf('  策略已收敛，完成迭代\n');
            break;
        end
    end
    
    % 最终值函数评估
    for i = 1:nW
        current_wealth = gW(i);
        c_ratio = C_star_t(i);
        a_ratio = A_star_t(i);
        
        % 计算最终策略下的值函数
        V_star_t(i) = policy_evaluation(c_ratio, a_ratio, current_wealth, ...
            V_star(:, t+1), gW, R_grid, R_prob, Rf, gamma, beta);
    end
    
    % 将临时数组的结果复制到全局数组
    C_star(:, t) = C_star_t;
    A_star(:, t) = A_star_t;
    V_star(:, t) = V_star_t;
    
    % 验证结果是否与理论一致
    mean_c = mean(C_star(:, t));
    mean_a = mean(A_star(:, t));
    std_a = std(A_star(:, t));
    
    fprintf('第%d期平均消费比例: %.4f (理论值: %.4f), 差异: %.2f%%\n', ...
        t, mean_c, c_optimal(t), abs(mean_c-c_optimal(t))/c_optimal(t)*100);
    fprintf('第%d期平均风险资产配置比例: %.4f (理论值: %.4f), 标准差: %.6f\n', ...
        t, mean_a, w_samuelson, std_a);
end

solve_time = toc;
fprintf('后向求解完成，耗时%.2f秒\n', solve_time);

%% 验证Samuelson结论
% 计算每个时期风险资产配置比例的统计量
A_mean_by_period = median(A_star, 1);
A_std_by_period = std(A_star, 0, 1);
A_relative_std = A_std_by_period ./ A_mean_by_period;

% 计算每个时期消费比例的统计量
C_mean_by_period = median(C_star, 1);
C_std_by_period = std(C_star, 0, 1);

fprintf('\n验证Samuelson结论：\n');
fprintf('风险资产配置比例的时间平均值: %.4f\n', mean(A_mean_by_period));
fprintf('风险资产配置比例的平均标准差: %.6f\n', mean(A_std_by_period));
fprintf('风险资产配置比例的平均变异系数: %.4f%%\n', mean(A_relative_std)*100);
fprintf('理论值与数值平均值的差异: %.4f%%\n', abs(w_samuelson-mean(A_mean_by_period))/w_samuelson*100);

fprintf('消费比例的时间平均值: %.4f\n', mean(C_mean_by_period));
fprintf('消费比例的平均标准差: %.6f\n', mean(C_std_by_period));

%% 验证消费比例是否与理论一致
fprintf('\n消费比例与理论值比较：\n');
periods_to_check = [1, 2, 3, T-2, T-1, T];  % 检查特定期数
for t_idx = 1:length(periods_to_check)
    t = periods_to_check(t_idx);
    if t == T
        theory_c = 1.0;  % 终期全部消费
    else
        theory_c = c_optimal(t);
    end
    actual_c = mean(C_star(:, t));
    fprintf('第%d期消费比例: %.4f (理论值: %.4f), 差异: %.2f%%\n', ...
        t, actual_c, theory_c, abs(actual_c-theory_c)/theory_c*100);
end

%% 保存模型结果
save('result_merton_PFI/model_results_finite.mat', 'C_star', 'A_star', 'V_star', 'gW', 'R_grid', 'R_prob', ...
    'beta', 'gamma', 'Rf', 'Rs_mean', 'sigma_s', 'w_merton', 'w_samuelson', 'c_optimal', 'T');

% 保存特定时期的数据用于可视化
periods_to_save = [1, 2, 3, T-2, T-1, T];  % 保存特定期数的数据
for t_idx = 1:length(periods_to_save)
    t = periods_to_save(t_idx);
    writematrix(C_star(:, t), sprintf('result_merton_PFI/C_t%d.csv', t));
    writematrix(A_star(:, t), sprintf('result_merton_PFI/A_t%d.csv', t));
    writematrix(V_star(:, t), sprintf('result_merton_PFI/V_t%d.csv', t));
end
writematrix(gW, 'result_merton_PFI/gW.csv');

%% 绘制策略函数和值函数，验证与财富的关系
% 选择几个有代表性的时期进行绘制
periods_to_plot = [1, T/4, T/2, 3*T/4, T];
periods_to_plot = max(1, min(T, round(periods_to_plot)));

% 1. 消费策略
figure('Position', [100, 100, 1200, 400]);
subplot(1, 3, 1);
hold on;
for t_idx = 1:length(periods_to_plot)
    t = periods_to_plot(t_idx);
    plot(gW, C_star(:, t), 'LineWidth', 1.5);
end
hold off;
title('不同期数的最优消费比例');
xlabel('财富');
ylabel('消费比例');
legend(arrayfun(@(t) sprintf('第%d期', t), periods_to_plot, 'UniformOutput', false), 'Location', 'best');
grid on;

% 2. 风险资产配置策略
subplot(1, 3, 2);
hold on;
for t_idx = 1:length(periods_to_plot)
    t = periods_to_plot(t_idx);
    plot(gW, A_star(:, t), 'LineWidth', 1.5);
end
plot([gW(1), gW(end)], [w_samuelson, w_samuelson], 'k--', 'LineWidth', 1.5);
hold off;
title('不同期数的最优风险资产配置比例');
xlabel('财富');
ylabel('风险资产配置比例');
legend([arrayfun(@(t) sprintf('第%d期', t), periods_to_plot, 'UniformOutput', false), {'理论值'}], 'Location', 'best');
grid on;

% 3. 值函数
subplot(1, 3, 3);
hold on;
for t_idx = 1:length(periods_to_plot)
    t = periods_to_plot(t_idx);
    if gamma ~= 0  % 非Bernoulli效用
        plot(gW, gamma*V_star(:, t), 'LineWidth', 1.5);
    else  % Bernoulli对数效用
        plot(gW, V_star(:, t), 'LineWidth', 1.5);
    end
end
hold off;
if gamma ~= 0
    title('不同期数的值函数(V×γ)');
else
    title('不同期数的值函数');
end
xlabel('财富');
ylabel('值函数');
legend(arrayfun(@(t) sprintf('第%d期', t), periods_to_plot, 'UniformOutput', false), 'Location', 'best');
grid on;

saveas(gcf, 'result_merton_PFI/policy_functions_finite.png');

% 4. 随时间变化的策略
figure('Position', [100, 100, 800, 600]);

% 消费比例随时间变化
subplot(2, 1, 1);
plot(1:T, C_mean_by_period, 'b-', 'LineWidth', 2);
hold on;
plot(1:T, c_optimal, 'r--', 'LineWidth', 1.5);
hold off;
title('消费比例随时间变化');
xlabel('时间');
ylabel('消费比例');
legend('数值平均值', '理论值', 'Location', 'best');
grid on;

% 风险资产配置比例随时间变化
subplot(2, 1, 2);
plot(1:T, A_mean_by_period, 'g-', 'LineWidth', 2);
hold on;
plot([1, T], [w_samuelson, w_samuelson], 'k--', 'LineWidth', 1.5);
hold off;
title('风险资产配置比例随时间变化');
xlabel('时间');
ylabel('风险资产配置比例');
legend('数值平均值', '理论值', 'Location', 'best');
grid on;

saveas(gcf, 'result_merton_PFI/strategy_time_evolution.png');

fprintf('\n有限期Merton模型(Samuelson设定)求解完成。\n');
fprintf('理论最优风险资产配置比例: %.4f\n', w_samuelson);
fprintf('数值解平均风险资产配置比例: %.4f\n', mean(A_mean_by_period));
fprintf('第一期理论消费比例: %.4f\n', c_optimal(1));
fprintf('第一期数值解消费比例: %.4f\n', C_mean_by_period(1));

%% 辅助函数
function value = period_value_function(x, W, next_period_V, wealth_grid, R_grid, R_prob, Rf, gamma, beta)
    % 计算当期值函数
    % x(1): 消费比例
    % x(2): 风险资产配置比例
    
    c_ratio = x(1);  % 消费比例
    a_ratio = x(2);  % 风险资产配置比例（相对剩余财富）
    
    % 计算消费金额
    consumption = c_ratio * W;
    
    % 计算剩余投资金额
    remaining_wealth = W * (1 - c_ratio);
    
    % 计算当期效用
    if consumption <= 0
        % 惩罚负消费或零消费
        value = -1e10;
        return;
    elseif gamma == 0  % Bernoulli对数效用
        current_utility = log(consumption);
    else
        current_utility = (consumption^gamma) / gamma;
    end
    
    % 计算期望未来效用
    expected_future_utility = 0;
    
    for i = 1:length(R_grid)
        % 计算下一期财富
        next_wealth = remaining_wealth * ((1-a_ratio)*Rf + a_ratio*R_grid(i));
        
        if next_wealth <= 0
            % 惩罚负财富
            future_value = -1e10;
        else
            % 对下一期财富插值计算下一期值函数
            future_value = interp1(wealth_grid, next_period_V, next_wealth, 'spline');
        end
        
        expected_future_utility = expected_future_utility + R_prob(i) * future_value;
    end
    
    % 总值函数
    value = current_utility + beta * expected_future_utility;
end

function value = samuelson_condition(w, R_grid, R_prob, Rf, gamma)
    % 计算Samuelson论文中的积分条件: ∫[(1-w)(1+r)+wZ]^(γ-1)(Z-1-r)dP(Z) = 0
    n = length(R_grid);
    total = 0;
    
    for i = 1:n
        portfolio_return = (1-w)*(Rf) + w*R_grid(i);
        if gamma ~= 0
            term = (portfolio_return^(gamma-1)) * (R_grid(i)-Rf) * R_prob(i);
        else  % Bernoulli对数效用的特殊情况
            term = (1/portfolio_return) * (R_grid(i)-Rf) * R_prob(i);
        end
        total = total + term;
    end
    
    value = abs(total);  % 返回绝对值以便最小化
end

function [z_grid, z_prob] = truncated_normal_grid(mu, sigma, n, lower_bound, upper_bound)
    % 创建截断正态分布的离散近似
    % mu: 均值
    % sigma: 标准差
    % n: 网格点数
    % lower_bound: 下界 (通常为0确保非负)
    % upper_bound: 上界 (可以是Inf)
    
    % 计算截断正态分布的均值和方差
    if lower_bound == -Inf && upper_bound == Inf
        % 标准正态分布情况
        adjusted_mu = mu;
        adjusted_sigma = sigma;
    else
        % 计算截断范围
        if lower_bound == -Inf
            a = -Inf;
        else
            a = (lower_bound - mu) / sigma;
        end
        
        if upper_bound == Inf
            b = Inf;
        else
            b = (upper_bound - mu) / sigma;
        end
        
        % 计算标准正态PDF和CDF
        if a == -Inf
            pdf_a = 0;
        else
            pdf_a = normpdf(a);
        end
        
        if b == Inf
            pdf_b = 0;
        else
            pdf_b = normpdf(b);
        end
        
        cdf_a = normcdf(a);
        cdf_b = normcdf(b);
        
        % 调整均值和方差
        adjusted_mu = mu + sigma * (pdf_a - pdf_b) / (cdf_b - cdf_a);
        adjusted_sigma = sigma * sqrt(1 + (a*pdf_a - b*pdf_b)/(cdf_b - cdf_a) - ((pdf_a - pdf_b)/(cdf_b - cdf_a))^2);
    end
    
    % 创建截断正态分布的离散近似网格
    % 等概率方法 (quantile-based)
    probs = linspace(0, 1, n+1);
    probs = probs(2:end) - (probs(2)-probs(1))/2;  % 取中点
    
    % 计算截断正态分布的分位数
    if lower_bound == -Inf && upper_bound == Inf
        z_grid = norminv(probs, mu, sigma);
    else
        % 计算截断正态分布的CDF
        trunc_probs = cdf_a + probs * (cdf_b - cdf_a);
        z_grid = norminv(trunc_probs, mu, sigma);
    end
    
    % 限制在边界内
    z_grid = max(lower_bound, min(upper_bound, z_grid));
    
    % 均匀概率
    z_prob = ones(n, 1) / n;
    
    return;
end

% 添加策略评估函数
function value = policy_evaluation(c_ratio, a_ratio, W, next_period_V, wealth_grid, R_grid, R_prob, Rf, gamma, beta)
    % 计算给定策略下的值函数
    % c_ratio: 消费比例
    % a_ratio: 风险资产配置比例
    % W: 当前财富
    
    % 计算消费金额
    consumption = c_ratio * W;
    
    % 计算剩余投资金额
    remaining_wealth = W * (1 - c_ratio);
    
    % 计算当期效用
    if consumption <= 0
        value = -1e10;
        return;
    elseif gamma == 0  % Bernoulli对数效用
        current_utility = log(consumption);
    else
        current_utility = (consumption^gamma) / gamma;
    end
    
    % 计算期望未来效用
    expected_future_utility = 0;
    
    for i = 1:length(R_grid)
        % 计算下一期财富
        next_wealth = remaining_wealth * ((1-a_ratio)*Rf + a_ratio*R_grid(i));
        
        if next_wealth <= 0
            % 负财富情况下的值函数
            future_value = -1e10;
        else
            % 对下一期财富插值计算下一期值函数
            future_value = interp1(wealth_grid, next_period_V, next_wealth, 'spline');
        end
        
        expected_future_utility = expected_future_utility + R_prob(i) * future_value;
    end
    
    % 总值函数
    value = current_utility + beta * expected_future_utility;
end 