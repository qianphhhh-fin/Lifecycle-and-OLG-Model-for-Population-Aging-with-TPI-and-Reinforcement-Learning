%% main_3stages_VFI.m
% 三期模型求解与模拟（VFI方法和解析解比较）
% 参考自main_baseline_VFI.m

clear;
close all;

%% 参数设置
% 基础参数
gamma = 2;      % 风险规避系数
beta = 0.95;    % 时间折扣因子

% 收益率参数
Rf = 1.02;      % 无风险收益率
Rp_bar = 1.03;  % 养老金固定收益率
RH = 1.10;      % 风险资产高收益率
RL = 0.94;      % 风险资产低收益率
p = 0.6;        % 高收益率概率

% 计算期望收益率和方差
mu1 = p*log(RH) + (1-p)*log(RL);   % 对数期望收益率
sigma1 = sqrt(p*(log(RH)-mu1)^2 + (1-p)*(log(RL)-mu1)^2);  % 对数标准差

% 初始财富设置（代替随机收入）
W1_init = 1.0;  % 固定初始财富

fprintf('模型参数：\n');
fprintf('gamma = %.2f, beta = %.2f\n', gamma, beta);
fprintf('Rf = %.2f, Rp_bar = %.2f\n', Rf, Rp_bar);
fprintf('RH = %.2f, RL = %.2f, p = %.2f\n', RH, RL, p);
fprintf('初始财富 W1 = %.2f\n', W1_init);

%% VFI求解

% 定义网格
n_W = 50;                          % 财富网格点数量
W_max = 100;
W_min = 0.25;
l_maxcash = log(W_max);
l_mincash = log(W_min);
stepcash = (l_maxcash - l_mincash) / (n_W - 1);
for i1 = 1:n_W
    lgcash(i1) = l_mincash + (i1-1) * stepcash;
end
for i1 = 1:n_W
    W_grid(i1) = exp(lgcash(i1));
end

% 定义决策变量网格
% n_c = 50;                          % 消费比例网格
% c_grid = linspace(0.01, 0.99, n_c);  % 消费比例网格
% n_alpha = 50;                      % 风险资产比例网格
% alpha_grid = linspace(0, 1, n_alpha); % 风险资产比例网格
% n_q = 50;                          % 养老金比例网格
% q_grid = linspace(0, 0.99, n_q);     % 养老金比例网格

% 初始化值函数和策略函数
V3 = zeros(n_W, 1);                % 第三期值函数
V2 = zeros(n_W, 1);                % 第二期值函数
V1 = zeros(n_W, 1);                % 第一期值函数（对每个网格点）

% 策略函数
c2_policy = zeros(n_W, 1);          % 第二期消费比例
alpha2_policy = zeros(n_W, 1);      % 第二期风险资产比例
c1_policy = zeros(n_W, 1);          % 第一期消费比例
alpha1_policy = zeros(n_W, 1);      % 第一期风险资产比例
q1_policy = zeros(n_W, 1);          % 第一期养老金比例

% 设置优化选项
options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'active-set', ...
    'MaxFunctionEvaluations', 1000, 'MaxIterations', 400, ...
    'OptimalityTolerance', 1e-10, 'StepTolerance', 1e-10, ...
    'ConstraintTolerance', 1e-10);

fprintf('\n开始VFI求解...\n');
tic;

%% 第三期求解 (消费全部财富)
fprintf('求解第三期...\n');
for i = 1:n_W
    W3 = W_grid(i);
    V3(i) = beta^2 * (W3^(1-gamma))/(1-gamma);
end

%% 第二期求解
fprintf('求解第二期...\n');
for i = 1:n_W
    W2 = W_grid(i);
    
    % % 初始网格搜索，找到最优的消费比例和风险资产比例
    % best_value = -inf;
    % best_c = 0;
    % best_alpha = 0;
    % 
    % for ic = 1:n_c
    %     c2 = c_grid(ic);
    %     for ia = 1:n_alpha
    %         alpha2 = alpha_grid(ia);
    % 
    %         % 计算期望效用
    %         EV = 0;
    %         % 考虑两种风险资产收益情况
    %         Rp_H = alpha2 * RH + (1-alpha2) * Rf;  % 高收益组合收益率
    %         Rp_L = alpha2 * RL + (1-alpha2) * Rf;  % 低收益组合收益率
    % 
    %         W3_H = (W2 * (1-c2)) * Rp_H;  % 高收益下的第三期财富
    %         W3_L = (W2 * (1-c2)) * Rp_L;  % 低收益下的第三期财富
    % 
    %         % 插值计算第三期值函数
    %         if W3_H > W_max
    %             V3_H = beta^2 * (W3_H^(1-gamma))/(1-gamma);
    %         else
    %             V3_H = interp1(W_grid, V3, W3_H, 'spline', 'extrap');
    %         end
    % 
    %         if W3_L > W_max
    %             V3_L = beta^2 * (W3_L^(1-gamma))/(1-gamma);
    %         else
    %             V3_L = interp1(W_grid, V3, W3_L, 'spline', 'extrap');
    %         end
    % 
    %         % 计算期望效用
    %         EV = p * V3_H + (1-p) * V3_L;
    % 
    %         % 当期效用
    %         utility = beta * ((W2 * c2)^(1-gamma)) / (1-gamma);
    % 
    %         % 总效用
    %         value = utility + EV;
    % 
    %         % 更新最优决策
    %         if value > best_value
    %             best_value = value;
    %             best_c = c2;
    %             best_alpha = alpha2;
    %         end
    %     end
    % end
    
    % 精细化搜索 - 使用fmincon优化
    objective = @(x) -second_period_value(x, W2, beta, gamma, Rf, RH, RL, p, W_grid, V3, W_max);
    x0 = [0.5, 0.5];
    lb = [0.001, 0];
    ub = [0.999, 1];
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    nonlcon = [];
    
    [x_opt, fval] = fmincon(objective, x0, A, b, Aeq, beq, lb, ub, nonlcon, options);
    
    c2_policy(i) = x_opt(1);
    alpha2_policy(i) = x_opt(2);
    V2(i) = -fval;
end

% 计算解析K2值
% 解析解中的最优alpha2
alpha2_analytical = (p*(RH-Rf) - (1-p)*(Rf-RL)) / (gamma * (p*(RH-mean([RH, RL]))^2 + (1-p)*(RL-mean([RH, RL]))^2));
alpha2_analytical = min(max(alpha2_analytical, 0), 1);  % 确保在[0,1]范围内

% 计算最优c2
Rp_mean = p*(alpha2_analytical*RH + (1-alpha2_analytical)*Rf) + (1-p)*(alpha2_analytical*RL + (1-alpha2_analytical)*Rf);
Rp_gamma_mean = p*((alpha2_analytical*RH + (1-alpha2_analytical)*Rf)^(1-gamma)) + (1-p)*((alpha2_analytical*RL + (1-alpha2_analytical)*Rf)^(1-gamma));
c2_analytical = 1 / (1 + beta^(1/gamma) * (Rp_gamma_mean)^(1/gamma));

% 计算K2
K2 = beta * (c2_analytical^(1-gamma) + beta * (1-c2_analytical)^(1-gamma) * Rp_gamma_mean);

%% 第一期求解
fprintf('求解第一期...\n');
% 对每个网格点求解最优策略
for i = 1:n_W
    W1 = W_grid(i);
    
    % 精细化搜索 - 使用fmincon优化
    objective = @(x) -first_period_value(x, W1, beta, gamma, Rf, RH, RL, p, Rp_bar, W_grid, V2, W_min, W_max, K2);
    x0 = [0.8, 0.5, 0.5];
    lb = [0.001, 0, 0];
    ub = [0.999, 1, 0.999];
    A = [1 0 1];  % c + q <= 1的约束
    b = 0.999;
    Aeq = [];
    beq = [];
    nonlcon = [];
    
    [x_opt, fval] = fmincon(objective, x0, A, b, Aeq, beq, lb, ub, nonlcon, options);
    
    c1_policy(i) = x_opt(1);
    alpha1_policy(i) = x_opt(2);
    q1_policy(i) = x_opt(3);
    V1(i) = -fval;
    
    % 显示进度
    if mod(i, 10) == 0
        fprintf('  已完成第一期网格点 %d/%d (%.1f%%)\n', i, n_W, i/n_W*100);
    end
end

vfi_time = toc;
fprintf('VFI求解完成，耗时 %.2f 秒\n', vfi_time);

%% 计算解析解

fprintf('\n计算解析解...\n');
tic;

% 第二期解析解
fprintf('第二期解析解：\n');

% 计算E[R_{p2}^{1-gamma}]和期望收益率
Rp2_gamma_mean = 0;
Rp_mean = 0;

% 考虑所有可能的alpha2值
n_alpha_test = 1000;
alpha2_grid = linspace(0, 1, n_alpha_test);
best_value = -inf;
best_alpha2 = 0;

for ia = 1:n_alpha_test
    alpha2 = alpha2_grid(ia);
    % 计算组合收益率的期望值
    E_Rp2 = p * (alpha2 * RH + (1-alpha2) * Rf) + (1-p) * (alpha2 * RL + (1-alpha2) * Rf);
    
    % 计算E[R_{p2}^{1-gamma}]
    E_Rp2_gamma = p * ((alpha2 * RH + (1-alpha2) * Rf)^(1-gamma)) + (1-p) * ((alpha2 * RL + (1-alpha2) * Rf)^(1-gamma));
    
    % 计算c2
    c2 = 1 / (1 + beta^(1/gamma) * (E_Rp2_gamma)^(1/gamma));
    
    % 计算one-period值函数
    value = beta * (c2^(1-gamma))/(1-gamma) + beta^2 * ((1-c2)^(1-gamma)) * E_Rp2_gamma/(1-gamma);
    
    if value > best_value
        best_value = value;
        best_alpha2 = alpha2;
        best_c2 = c2;
        best_E_Rp2_gamma = E_Rp2_gamma;
    end
end

% 最优alpha2和c2
alpha2_analytical = best_alpha2;
c2_analytical = best_c2;
Rp2_gamma_mean = best_E_Rp2_gamma;

% 计算K2
K2 = beta * (c2_analytical^(1-gamma) + beta * (1-c2_analytical)^(1-gamma) * Rp2_gamma_mean);

fprintf('alpha2* = %.4f, c2* = %.4f, K2 = %.4f\n', alpha2_analytical, c2_analytical, K2);

% 第一期解析解
fprintf('计算第一期解析解（完整理论公式）...\n');

% 计算期望收益率和方差
R_mean = p*RH + (1-p)*RL;
var_R1 = p*(RH-R_mean)^2 + (1-p)*(RL-R_mean)^2;

% 创建投资组合收益率可能性
R1_values = [RH, RL];
R1_probs = [p, 1-p];

% 初步估计 alpha1, c1, q1 用于迭代计算
alpha1_initial = (R_mean - Rf) / (gamma * var_R1);
alpha1_initial = min(max(alpha1_initial, 0), 1);  % 确保在[0,1]范围内

c1_initial = 0.3;  % 初始猜测消费比例
q1_initial = 0.1;  % 初始猜测养老金比例

% 初始化迭代
alpha1_curr = alpha1_initial;
c1_curr = c1_initial;
q1_curr = q1_initial;

% 迭代参数
max_iter = 100;
tolerance = 1e-6;
converged = false;

fprintf('迭代计算完整解析解...\n');
for iter = 1:max_iter
    % 保存当前值
    alpha1_prev = alpha1_curr;
    c1_prev = c1_curr;
    q1_prev = q1_curr;
    
    % 为每个可能的收益率计算权重和期望值
    S_values = zeros(size(R1_values));
    weights = zeros(size(R1_values));
    
    for i = 1:length(R1_values)
        R1 = R1_values(i);
        Rp1 = alpha1_curr * R1 + (1-alpha1_curr) * Rf;
        
        % 计算财富状态 S
        S_values(i) = (W1_init - W1_init*c1_curr - W1_init*q1_curr) * Rp1 + W1_init*q1_curr * Rp_bar;
        
        % 计算边际效用权重 S^{-gamma}
        weights(i) = S_values(i)^(-gamma);
    end
    
    % 规范化权重
    normalized_weights = weights / sum(weights .* R1_probs);
    
    % 计算各种加权期望值
    E_S_gamma_Rp1 = 0;  % E[S^{-gamma} * Rp1]
    E_S_gamma = 0;      % E[S^{-gamma}]
    E_S_gamma_Rp1_Rp_diff = 0;  % E[S^{-gamma} * (Rp1 - Rp_bar)]
    E_S_gamma_Rp1_Rp_diff_squared = 0;  % E[S^{-gamma} * (Rp1 - Rp_bar)^2]
    E_S_gamma_R1_Rf_diff = 0;   % E[S^{-gamma} * (R1 - Rf)]
    E_S_gamma_R1_R1mean_diff_times_R1_Rf_diff = 0;  % E[S^{-gamma} * (R1 - E[R1]) * (R1 - Rf)]
    
    for i = 1:length(R1_values)
        R1 = R1_values(i);
        prob = R1_probs(i);
        Rp1 = alpha1_curr * R1 + (1-alpha1_curr) * Rf;
        
        E_S_gamma = E_S_gamma + prob * weights(i);
        E_S_gamma_Rp1 = E_S_gamma_Rp1 + prob * weights(i) * Rp1;
        E_S_gamma_Rp1_Rp_diff = E_S_gamma_Rp1_Rp_diff + prob * weights(i) * (Rp1 - Rp_bar);
        E_S_gamma_Rp1_Rp_diff_squared = E_S_gamma_Rp1_Rp_diff_squared + prob * weights(i) * (Rp1 - Rp_bar)^2;
        E_S_gamma_R1_Rf_diff = E_S_gamma_R1_Rf_diff + prob * weights(i) * (R1 - Rf);
        E_S_gamma_R1_R1mean_diff_times_R1_Rf_diff = E_S_gamma_R1_R1mean_diff_times_R1_Rf_diff + ...
            prob * weights(i) * (R1 - R_mean) * (R1 - Rf);
    end
    
    % 更新消费比例 c1* 根据理论公式
    consumption_adjustment = (E_S_gamma_Rp1 / E_S_gamma)^(1/gamma);
    c1_new = 1 / (1 + (beta * Rp2_gamma_mean^(1/gamma)) * consumption_adjustment);
    
    % 更新养老金比例 q1*
    if abs(E_S_gamma_Rp1_Rp_diff_squared) > 1e-10
        q1_new = (E_S_gamma_Rp1_Rp_diff / E_S_gamma_Rp1_Rp_diff_squared) * (1 - c1_new);
    else
        q1_new = 0;  % 如果分母接近零，设为0以避免数值问题
    end
    
    % 约束 q1 到合理范围
    q1_new = min(max(q1_new, 0), 1 - c1_new);
    
    % 更新风险资产比例 alpha1*
    if abs(E_S_gamma_R1_R1mean_diff_times_R1_Rf_diff) > 1e-10
        alpha1_new = E_S_gamma_R1_Rf_diff / E_S_gamma_R1_R1mean_diff_times_R1_Rf_diff;
    else
        alpha1_new = alpha1_initial;  % 如果分母接近零，使用初始估计
    end
    
    % 约束 alpha1 到合理范围
    alpha1_new = min(max(alpha1_new, 0), 1);
    
    % 应用松弛因子增加收敛稳定性
    relaxation = 0.7;
    c1_curr = relaxation * c1_new + (1 - relaxation) * c1_prev;
    q1_curr = relaxation * q1_new + (1 - relaxation) * q1_prev;
    alpha1_curr = relaxation * alpha1_new + (1 - relaxation) * alpha1_prev;
    
    % 检查收敛性
    error_c1 = abs(c1_curr - c1_prev);
    error_q1 = abs(q1_curr - q1_prev);
    error_alpha1 = abs(alpha1_curr - alpha1_prev);
    max_error = max([error_c1, error_q1, error_alpha1]);
    
    if mod(iter, 10) == 0
        fprintf('  迭代 %d: c1 = %.4f, q1 = %.4f, alpha1 = %.4f, 最大误差 = %.6f\n', ...
            iter, c1_curr, q1_curr, alpha1_curr, max_error);
    end
    
    if max_error < tolerance
        converged = true;
        break;
    end
end

% 最终解析解
c1_analytical = c1_curr;
q1_analytical = q1_curr;
alpha1_analytical = alpha1_curr;

% 计算未约束的养老金比例（理论原始公式）
q1_unconstrained = (E_S_gamma_Rp1_Rp_diff / E_S_gamma_Rp1_Rp_diff_squared) * (1 - c1_analytical);

% 理论组合收益率均值
Rp1_mean = alpha1_analytical*R_mean + (1-alpha1_analytical)*Rf;

% 检验理论条件
Rp_threshold_low = Rp1_mean - sqrt(gamma * var_R1 * alpha1_analytical^2);
Rp_threshold_high = Rp1_mean + sqrt(gamma * var_R1 * alpha1_analytical^2);
condition1 = R_mean > Rf;
condition2 = Rp_bar > Rp_threshold_low && Rp_bar < Rp_threshold_high;
condition3 = W1_init > W1_init*c1_analytical + W1_init*q1_analytical;

if converged
    fprintf('解析解迭代在 %d 次迭代后收敛, 最终误差 = %.6f\n', iter, max_error);
else
    fprintf('警告：解析解迭代在 %d 次迭代后仍未收敛, 最终误差 = %.6f\n', max_iter, max_error);
end

% 打印解析解结果
fprintf('\n完整理论解析解结果:\n');
fprintf('消费比例 c1* = %.4f\n', c1_analytical);
fprintf('风险资产比例 alpha1* = %.4f\n', alpha1_analytical);
fprintf('养老金比例 q1* = %.4f (未约束: %.4f)\n', q1_analytical, q1_unconstrained);

% 打印养老金收益率和组合收益率对比
fprintf('\n养老金收益率与投资组合收益率对比:\n');
fprintf('养老金固定收益率 Rp_bar = %.4f\n', Rp_bar);
fprintf('最优投资组合期望收益率 = %.4f (理论允许范围: %.4f-%.4f)\n', ...
    Rp1_mean, Rp_threshold_low, Rp_threshold_high);
fprintf('期望收益差值 = %.4f\n', Rp1_mean - Rp_bar);

% fprintf('\n解析解条件检验:\n');
% fprintf('条件1: E[R1] > Rf: %s (%.4f > %.4f)\n', condition1 ? '满足' : '不满足', R_mean, Rf);
% fprintf('条件2: Rp_bar在允许范围内: %s (%.4f < %.4f < %.4f)\n', ...
%     condition2 ? '满足' : '不满足', Rp_threshold_low, Rp_bar, Rp_threshold_high);
% fprintf('条件3: 可行性: %s (%.4f > %.4f + %.4f)\n', ...
%     condition3 ? '满足' : '不满足', W1_init, W1_init*c1_analytical, W1_init*q1_analytical);

analytical_time = toc;
fprintf('解析解计算完成，耗时 %.2f 秒\n', analytical_time);

% 定义中点索引（用于比较）
mid_point = ceil(n_W/2);

%% 数值模拟
fprintf('\n进行数值模拟...\n');
tic;

% 模拟参数
n_sim = 10000;        % 降低模拟次数以加快速度
n_periods = 3;        % 模拟期数 (三期模型)

% 初始化结果存储
sim_c1 = zeros(1, n_sim);   % 消费比例
sim_alpha1 = zeros(1, n_sim);
sim_q1 = zeros(1, n_sim);
sim_c2 = zeros(1, n_sim);   % 第二期消费
sim_alpha2 = zeros(1, n_sim);
sim_wealth = zeros(n_periods, n_sim);  % 三期的财富水平

% 设置随机数种子以保证可重复性
rng(12345);

% 为每个模拟生成随机收益序列
R_scenarios = zeros(n_periods-1, n_sim);  % 只需要两期的收益率
for t = 1:n_periods-1
    % 以概率p生成高收益率，以概率1-p生成低收益率
    high_return_idx = rand(1, n_sim) < p;
    R_scenarios(t, high_return_idx) = RH;
    R_scenarios(t, ~high_return_idx) = RL;
end

% 初始财富（使用网格中点）
W1_idx = ceil(n_W/2);  % 使用网格中点作为模拟的起始财富
W1 = W_grid(W1_idx);
sim_wealth(1, :) = W1;

% 第一期决策
% 使用数值解策略（从策略函数中查找）
c1_ratio = c1_policy(W1_idx);
alpha1_ratio = alpha1_policy(W1_idx);
q1_ratio = q1_policy(W1_idx);

% 记录第一期决策
sim_c1(:) = c1_ratio;
sim_alpha1(:) = alpha1_ratio;
sim_q1(:) = q1_ratio;

% 计算每个模拟的第二期财富
for i = 1:n_sim
    % 第一期投资组合收益率
    Rp1 = alpha1_ratio * R_scenarios(1, i) + (1-alpha1_ratio) * Rf;
    
    % 计算第二期财富
    W2 = (W1 - W1*c1_ratio - W1*q1_ratio) * Rp1 + W1*q1_ratio * Rp_bar;
    sim_wealth(2, i) = W2;
    
    % 第二期决策 - 使用插值获取策略
    if W2 <= W_min
        c2_ratio = c2_policy(1);
        alpha2_ratio = alpha2_policy(1);
    elseif W2 >= W_max
        c2_ratio = c2_policy(end);
        alpha2_ratio = alpha2_policy(end);
    else
        % 插值获取策略
        c2_ratio = interp1(W_grid, c2_policy, W2, 'linear');
        alpha2_ratio = interp1(W_grid, alpha2_policy, W2, 'linear');
    end
    
    % 记录第二期决策
    sim_c2(i) = c2_ratio;
    sim_alpha2(i) = alpha2_ratio;
    
    % 计算第三期财富
    Rp2 = alpha2_ratio * R_scenarios(2, i) + (1-alpha2_ratio) * Rf;
    W3 = (W2 * (1-c2_ratio)) * Rp2;
    sim_wealth(3, i) = W3;
end

% 计算所有模拟样本的平均值
% 第一期决策
overall_sim_c1 = mean(sim_c1);
overall_sim_alpha1 = mean(sim_alpha1);
overall_sim_q1 = mean(sim_q1);

% 第二期决策
overall_sim_c2 = mean(sim_c2);
overall_sim_alpha2 = mean(sim_alpha2);

% 财富
overall_sim_W2 = mean(sim_wealth(2,:));
overall_sim_W3 = mean(sim_wealth(3,:));

% 打印简化的结果对比
fprintf('\n=== 每期数值模拟解与解析解对比 ===\n');

fprintf('\n第一期决策对比:\n');
fprintf('                 消费比例(c1)  风险资产比例(alpha1)  养老金比例(q1)\n');
fprintf('数值模拟平均值:    %.4f          %.4f              %.4f\n', ...
    overall_sim_c1, overall_sim_alpha1, overall_sim_q1);
fprintf('数值解:           %.4f          %.4f              %.4f\n', ...
    c1_policy(W1_idx), alpha1_policy(W1_idx), q1_policy(W1_idx));
fprintf('解析解:           %.4f          %.4f              %.4f  (未约束q1: %.4f)\n', ...
    c1_analytical, alpha1_analytical, q1_analytical, q1_unconstrained);

fprintf('\n第二期决策对比:\n');
fprintf('                 消费比例(c2)  风险资产比例(alpha2)\n');
fprintf('数值模拟平均值:    %.4f          %.4f\n', ...
    overall_sim_c2, overall_sim_alpha2);
fprintf('数值解平均值:      %.4f          %.4f   (取中点值)\n', ...
    c2_policy(mid_point), alpha2_policy(mid_point));
fprintf('解析解值:          %.4f          %.4f\n', ...
    c2_analytical, alpha2_analytical);

fprintf('\n财富轨迹:\n');
fprintf('初始财富: %.4f\n', W1_init);
fprintf('平均第二期财富: %.4f\n', overall_sim_W2);
fprintf('平均第三期财富: %.4f\n', overall_sim_W3);

fprintf('\n相对误差分析 (|数值解-解析解|/解析解):\n');
fprintf('第一期: c1 = %.2f%%, alpha1 = %.2f%%, q1 = %.2f%%\n', ...
    abs(c1_policy(W1_idx)-c1_analytical)/c1_analytical*100, ...
    abs(alpha1_policy(W1_idx)-alpha1_analytical)/alpha1_analytical*100, ...
    abs(q1_policy(W1_idx)-q1_analytical)/max(0.0001, q1_analytical)*100);
fprintf('第二期: c2 = %.2f%%, alpha2 = %.2f%%\n', ...
    abs(c2_policy(mid_point)-c2_analytical)/c2_analytical*100, ...
    abs(alpha2_policy(mid_point)-alpha2_analytical)/alpha2_analytical*100);

sim_time = toc;
fprintf('\n数值模拟完成，耗时 %.2f 秒\n', sim_time);

%% 比较数值解和解析解
fprintf('\n比较数值解和解析解:\n');

% 定义要比较的财富点
W1_compare = W1_init;  % 使用初始设定的财富点进行比较
W1_idx = find(abs(W_grid - W1_compare) == min(abs(W_grid - W1_compare)), 1);  % 找到最接近的网格点

% 第二期比较
fprintf('第二期策略函数（取中点值）:\n');
fprintf('数值解: c2 = %.4f, alpha2 = %.4f\n', c2_policy(mid_point), alpha2_policy(mid_point));
fprintf('解析解: c2 = %.4f, alpha2 = %.4f\n', c2_analytical, alpha2_analytical);
fprintf('相对误差: c2 = %.2f%%, alpha2 = %.2f%%\n', ...
    abs(c2_policy(mid_point)-c2_analytical)/c2_analytical*100, ...
    abs(alpha2_policy(mid_point)-alpha2_analytical)/alpha2_analytical*100);

% 第一期比较
fprintf('\n第一期策略 (W1 = %.4f):\n', W_grid(W1_idx));
fprintf('数值解: c1 = %.4f, alpha1 = %.4f, q1 = %.4f\n', ...
    c1_policy(W1_idx), alpha1_policy(W1_idx), q1_policy(W1_idx));
fprintf('解析解: c1 = %.4f, alpha1 = %.4f, q1 = %.4f\n', ...
    c1_analytical, alpha1_analytical, q1_analytical);
fprintf('未约束q1值: %.4f\n', q1_unconstrained);

% 计算相对误差,避免除以0
c1_err = abs(c1_policy(W1_idx)-c1_analytical)/c1_analytical*100;
alpha1_err = abs(alpha1_policy(W1_idx)-alpha1_analytical)/alpha1_analytical*100;
if q1_analytical > 0
    q1_err = abs(q1_policy(W1_idx)-q1_analytical)/q1_analytical*100;
else
    q1_err = NaN; % 当解析解q1为0时,用NaN表示无法计算相对误差
end

fprintf('相对误差: c1 = %.2f%%, alpha1 = %.2f%%, q1 = ', c1_err, alpha1_err);
if isnan(q1_err)
    fprintf('无法计算(解析解为0)\n');
else
    fprintf('%.2f%%\n', q1_err);
end

% 绘制第一期策略函数
figure;
subplot(3,1,1);
plot(W_grid, c1_policy, 'b-', 'LineWidth', 2);
hold on;
plot([W_min, W_max], [c1_analytical, c1_analytical], 'r--', 'LineWidth', 2);
title('第一期消费策略');
xlabel('财富 (W_1)');
ylabel('消费比例 (c_1)');
legend('数值解', '解析解');
grid on;

subplot(3,1,2);
plot(W_grid, alpha1_policy, 'b-', 'LineWidth', 2);
hold on;
plot([W_min, W_max], [alpha1_analytical, alpha1_analytical], 'r--', 'LineWidth', 2);
title('第一期风险资产投资策略');
xlabel('财富 (W_1)');
ylabel('风险资产比例 (\alpha_1)');
legend('数值解', '解析解');
grid on;

subplot(3,1,3);
plot(W_grid, q1_policy, 'b-', 'LineWidth', 2);
hold on;
plot([W_min, W_max], [q1_analytical, q1_analytical], 'r--', 'LineWidth', 2);
title('第一期养老金策略');
xlabel('财富 (W_1)');
ylabel('养老金比例 (q_1)');
legend('数值解', '解析解');
grid on;

%% 辅助函数
function value = second_period_value(x, W2, beta, gamma, Rf, RH, RL, p, W_grid, V3, W_max)
    c2 = x(1);
    alpha2 = x(2);
    
    % 计算期望效用
    EV = 0;
    
    % 考虑两种风险资产收益情况
    Rp_H = alpha2 * RH + (1-alpha2) * Rf;  % 高收益组合收益率
    Rp_L = alpha2 * RL + (1-alpha2) * Rf;  % 低收益组合收益率
    
    W3_H = (W2 * (1-c2)) * Rp_H;  % 高收益下的第三期财富
    W3_L = (W2 * (1-c2)) * Rp_L;  % 低收益下的第三期财富
    
    % 插值计算第三期值函数
    if W3_H > W_max
        V3_H = beta^2 * (W3_H^(1-gamma))/(1-gamma);
    else
        V3_H = interp1(W_grid, V3, W3_H, 'spline', 'extrap');
    end
    
    if W3_L > W_max
        V3_L = beta^2 * (W3_L^(1-gamma))/(1-gamma);
    else
        V3_L = interp1(W_grid, V3, W3_L, 'spline', 'extrap');
    end
    
    % 计算期望效用
    EV = p * V3_H + (1-p) * V3_L;
    
    % 当期效用
    utility = beta * ((W2 * c2)^(1-gamma)) / (1-gamma);
    
    % 总效用
    value = utility + EV;
end

function value = first_period_value(x, W1, beta, gamma, Rf, RH, RL, p, Rp_bar, W_grid, V2, W_min, W_max, K2)
    c1 = x(1);
    alpha1 = x(2);
    q1 = x(3);
    
    % 计算期望效用
    EV = 0;
    
    % 考虑收入和风险资产收益的两种组合
    for i_r = 1:2
        if i_r == 1
            R1 = RH;
            prob_r = p;
        else
            R1 = RL;
            prob_r = 1-p;
        end
        
        % 计算投资组合收益率
        Rp1 = alpha1 * R1 + (1-alpha1) * Rf;
        
        % 计算第二期财富
        W2 = (W1 - W1*c1 - W1*q1) * Rp1 + W1*q1 * Rp_bar;
        
        % 插值计算第二期值函数
        if W2 > W_max
            V2_interp = K2 * (W2^(1-gamma))/(1-gamma);
        elseif W2 < W_min
            V2_interp = K2 * (W2^(1-gamma))/(1-gamma);
        else
            V2_interp = interp1(W_grid, V2, W2, 'spline', 'extrap');
        end
        
        % 更新期望效用
        EV = EV + prob_r * V2_interp;
    end
    
    % 当期效用
    utility = ((W1 * c1)^(1-gamma)) / (1-gamma);
    
    % 总效用
    value = utility + EV;
end 