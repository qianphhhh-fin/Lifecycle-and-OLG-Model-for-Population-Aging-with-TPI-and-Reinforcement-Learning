clear
% 三期生命周期模型求解（使用策略函数迭代方法）
% 前两期为工作阶段(t=1,2)，第三期为退休阶段(t=3)
% 所有决策变量都表示为比例：
% C - 消费比例（占当前财富）
% Q - 养老金购买比例（占当前财富）
% A - 风险资产投资比例（占剩余投资）
% 剩余部分 1-C-Q-A 为无风险资产投资比例
%
% 优化版本：
% - 使用interp2函数进行双线性插值
% - 使用矩阵运算代替循环，提高计算效率
% - 采用策略函数迭代(PFI)方法求解动态规划问题

% 创建保存结果的目录
if ~exist('result_3stages_PFI', 'dir')
    mkdir('result_3stages_PFI');
end

% 设置随机种子以获得可重复的结果
rng(42);

%% 模型参数设置
% 效用函数参数
beta = 0.95;    % 折现因子
gamma = 3.0;    % 相对风险厌恶系数

% 收入和收益率参数
sigma_eps = 0.3;  % 劳动收入的标准差
Rp = 1.03;      % 养老基金固定收益率
Rf = 1.02;      % 无风险资产收益率
Rs_mean = 1.06;  % 风险资产平均收益率(水平值，非对数)
sigma_s = 0.18;  % 风险资产对数收益率的标准差

% 状态变量网格设置
nW = 51;        % 现金财富网格点数
nF = 51;        % 养老金余额网格点数
nY = 3;         % 劳动收入随机冲击的网格点数
nR = 3;         % 风险资产收益率随机冲击的网格点数

% 策略函数迭代参数
max_iter = 200;  % 最大迭代次数
tol = 1e-7;      % 收敛容差
alpha = 0.8;     % 迭代步长

%% 初始化状态变量网格和随机冲击
% 现金财富网格 (对数均匀网格)
W_min = 0.1;
W_max = 50.0;
log_W_min = log(W_min);
log_W_max = log(W_max);
step_W = (log_W_max - log_W_min) / (nW - 1);
gW = zeros(nW, 1);
for i = 1:nW
    gW(i) = exp(log_W_min + (i-1) * step_W);
end

% 养老金余额网格 (对数均匀网格)
F_min = 0.0;  % 最小值为0，表示初期没有养老金
F_max = 50.0;
log_F_min = log(F_min + 0.1);  % 加一个小数以避免log(0)
log_F_max = log(F_max + 0.1);
step_F = (log_F_max - log_F_min) / (nF - 1);
gF = zeros(nF, 1);
for i = 1:nF
    gF(i) = exp(log_F_min + (i-1) * step_F) - 0.1;
    if gF(i) < 0
        gF(i) = 0;  % 确保没有负值
    end
end

% 计算对数正态分布的对数均值参数
mu_s = log(Rs_mean) - 0.5 * sigma_s^2;  % 对数正态分布的均值参数
r_f = log(Rf);  % 对数无风险收益率

% 定义离散化的收入和收益率随机冲击
% 劳动收入为对数正态分布
[Y_grid, Y_prob] = tauchen_method(0, 0, sigma_eps, nY, 3);
Y_grid = exp(Y_grid);  % 转换为水平值

% 风险资产收益率为对数正态分布
[logR_grid, R_prob] = tauchen_method(mu_s, 0, sigma_s, nR, 3);
R_grid = exp(logR_grid);  % 从对数空间转换为实际收益率，即对数正态分布

% 初始化策略函数和值函数
C = zeros(nW, nF, 3);  % 消费比例，三个时期
A = zeros(nW, nF, 3);  % 风险资产配置比例
Q = zeros(nW, nF, 3);  % 养老金购买比例
V = zeros(nW, nF, 3);  % 值函数

% 初始化联合概率权重 (对于工作期随机性)
joint_prob = zeros(nY, nR);
for i = 1:nY
    for j = 1:nR
        joint_prob(i, j) = Y_prob(i) * R_prob(j);
    end
end

%% 初始化优化选项
% 改进的SQP算法设置，更精确的求解参数
options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp', ...
    'MaxFunctionEvaluations', 1000, 'MaxIterations', 400, ...
    'OptimalityTolerance', 1e-8, 'StepTolerance', 1e-8, ...
    'ConstraintTolerance', 1e-8, 'FiniteDifferenceType', 'central');

% 多起点优化设置
n_starts = 3;  % 优化起点数量

%% 第3期求解 (退休期)
fprintf('求解第3期 (退休期) 问题...\n');

% 退休期的策略直接求解 (无需迭代，因为这是最后一期)
for i = 1:nW
    for j = 1:nF
        % 退休期收入就是养老金余额
        income = gF(j);
        
        % 总可用财富
        total_wealth = gW(i) + income;
        
        % 最优策略：消费所有可用财富
        % 注意：此处 C 表示消费比例，但第3期消费所有财富，所以设为1
        C(i, j, 3) = 1.0;
        A(i, j, 3) = 0;  % 无需投资配置
        Q(i, j, 3) = 0;  % 无需购买养老金
        
        % 计算值函数
        if total_wealth > 0
            V(i, j, 3) = (total_wealth^(1-gamma)) / (1-gamma);
        else
            V(i, j, 3) = -1e10;  % 一个极小值表示不可行的状态
        end
    end
end

fprintf('第3期求解完成\n');

%% 第2期求解 (工作期)
fprintf('求解第2期 (工作期) 问题...\n');
tic;

% 初始化第2期策略
for i = 1:nW
    for j = 1:nF
        C(i, j, 2) = 0.5;  % 初始消费为财富的50%
        A(i, j, 2) = 0.3;  % 初始风险资产配置为30%
        Q(i, j, 2) = 0.1;  % 初始养老金购买为财富的10%
    end
end

% 策略函数迭代 (第2期)
for iter = 1:max_iter
    % 创建新的策略函数
    C_new = zeros(nW, nF);
    A_new = zeros(nW, nF);
    Q_new = zeros(nW, nF);
    
    % 逐点优化
    parfor i = 1:nW
        for j = 1:nF
            % 优化约束
            lb = [0, 0, 0];  % 下界
            ub = [1, 1, 1];  % 上界
            
            % 非线性约束 (消费比例+养老金购买比例+风险资产投资比例<=1)
            nonlcon = @(x) deal(x(1) + x(2) + x(3) - 1, []);
            
            % 目标函数
            obj_fun = @(x) -fun_valuefunc_working(x, gW(i), gF(j), V(:,:,3), Y_grid, R_grid, joint_prob, Rf, Rp, gamma, beta, gW, gF, 2);
            
            % 多起点优化
            best_val = Inf;
            best_x = [0, 0, 0];
            
            for start = 1:n_starts
                if start == 1
                    % 第一个起点：使用当前策略
                    x0 = [C(i, j, 2), A(i, j, 2), Q(i, j, 2)];
                else
                    % 生成随机起点
                    % 确保满足约束 x1 + x2 + x3 <= 1
                    sum_x = 2;  % 初始化为大于1的值
                    while sum_x > 1
                        rand_x = rand(1, 3);
                        sum_x = sum(rand_x);
                        if sum_x > 1
                            rand_x = rand_x / (1.5 * sum_x);  % 缩小以提高接受率
                        end
                    end
                    x0 = rand_x;
                end
                
                % 执行优化
                [x_tmp, fval_tmp] = fmincon(obj_fun, x0, [], [], [], [], lb, ub, nonlcon, options);
                
                % 检查是否是更好的解
                if fval_tmp < best_val
                    best_val = fval_tmp;
                    best_x = x_tmp;
                end
            end
            
            % 更新策略为最优解
            C_new(i, j) = best_x(1);
            A_new(i, j) = best_x(2);
            Q_new(i, j) = best_x(3);
        end
    end
    
    % 计算策略变化
    policy_change = max(max(abs(C_new - C(:,:,2)))) + ...
                   max(max(abs(A_new - A(:,:,2)))) + ...
                   max(max(abs(Q_new - Q(:,:,2))));
    
    % 更新策略 (使用步长)
    C(:,:,2) = (1-alpha) * C(:,:,2) + alpha * C_new;
    A(:,:,2) = (1-alpha) * A(:,:,2) + alpha * A_new;
    Q(:,:,2) = (1-alpha) * Q(:,:,2) + alpha * Q_new;
    
    % 更新值函数
    for i = 1:nW
        for j = 1:nF
            V(i, j, 2) = fun_valuefunc_working([C(i,j,2), A(i,j,2), Q(i,j,2)], ...
                gW(i), gF(j), V(:,:,3), Y_grid, R_grid, joint_prob, Rf, Rp, gamma, beta, gW, gF, 2);
        end
    end
    
    % 检查收敛
    if policy_change < tol
        fprintf('第2期在%d次迭代后收敛\n', iter);
        break;
    end
    
    fprintf('第2期迭代%d次，策略变化量: %f\n', iter, policy_change);
    
    % 达到最大迭代次数
    if iter == max_iter
        fprintf('警告：第2期在%d次迭代后仍未收敛\n', max_iter);
    end
end

fprintf('第2期求解完成，耗时%.2f秒\n', toc);

%% 第1期求解 (工作期)
fprintf('求解第1期 (工作期) 问题...\n');
tic;

% 初始化第1期策略
for i = 1:nW
    for j = 1:nF
        C(i, j, 1) = 0.5;  % 初始消费为财富的50%
        A(i, j, 1) = 0.3;  % 初始风险资产配置为30%
        Q(i, j, 1) = 0.1;  % 初始养老金购买为财富的10%
    end
end

% 策略函数迭代 (第1期)
for iter = 1:max_iter
    % 创建新的策略函数
    C_new = zeros(nW, nF);
    A_new = zeros(nW, nF);
    Q_new = zeros(nW, nF);
    
    % 逐点优化
    parfor i = 1:nW
        for j = 1:nF
            % 优化约束
            lb = [0, 0, 0];  % 下界
            ub = [1, 1, 1];  % 上界
            
            % 非线性约束 (消费比例+养老金购买比例+风险资产投资比例<=1)
            nonlcon = @(x) deal(x(1) + x(2) + x(3) - 1, []);
            
            % 目标函数
            obj_fun = @(x) -fun_valuefunc_working(x, gW(i), gF(j), V(:,:,2), Y_grid, R_grid, joint_prob, Rf, Rp, gamma, beta, gW, gF, 1);
            
            % 多起点优化
            best_val = Inf;
            best_x = [0, 0, 0];
            
            for start = 1:n_starts
                if start == 1
                    % 第一个起点：使用当前策略
                    x0 = [C(i, j, 1), A(i, j, 1), Q(i, j, 1)];
                else
                    % 生成随机起点
                    % 确保满足约束 x1 + x2 + x3 <= 1
                    sum_x = 2;  % 初始化为大于1的值
                    while sum_x > 1
                        rand_x = rand(1, 3);
                        sum_x = sum(rand_x);
                        if sum_x > 1
                            rand_x = rand_x / (1.5 * sum_x);  % 缩小以提高接受率
                        end
                    end
                    x0 = rand_x;
                end
                
                % 执行优化
                [x_tmp, fval_tmp] = fmincon(obj_fun, x0, [], [], [], [], lb, ub, nonlcon, options);
                
                % 检查是否是更好的解
                if fval_tmp < best_val
                    best_val = fval_tmp;
                    best_x = x_tmp;
                end
            end
            
            % 更新策略为最优解
            C_new(i, j) = best_x(1);
            A_new(i, j) = best_x(2);
            Q_new(i, j) = best_x(3);
        end
    end
    
    % 计算策略变化
    policy_change = max(max(abs(C_new - C(:,:,1)))) + ...
                   max(max(abs(A_new - A(:,:,1)))) + ...
                   max(max(abs(Q_new - Q(:,:,1))));
    
    % 更新策略 (使用步长)
    C(:,:,1) = (1-alpha) * C(:,:,1) + alpha * C_new;
    A(:,:,1) = (1-alpha) * A(:,:,1) + alpha * A_new;
    Q(:,:,1) = (1-alpha) * Q(:,:,1) + alpha * Q_new;
    
    % 更新值函数
    for i = 1:nW
        for j = 1:nF
            V(i, j, 1) = fun_valuefunc_working([C(i,j,1), A(i,j,1), Q(i,j,1)], ...
                gW(i), gF(j), V(:,:,2), Y_grid, R_grid, joint_prob, Rf, Rp, gamma, beta, gW, gF, 1);
        end
    end
    
    % 检查收敛
    if policy_change < tol
        fprintf('第1期在%d次迭代后收敛\n', iter);
        break;
    end
    
    fprintf('第1期迭代%d次，策略变化量: %f\n', iter, policy_change);
    
    % 达到最大迭代次数
    if iter == max_iter
        fprintf('警告：第1期在%d次迭代后仍未收敛\n', max_iter);
    end
end

fprintf('第1期求解完成，耗时%.2f秒\n', toc);

%% 保存结果
save('result_3stages_PFI/model_results.mat', 'C', 'A', 'Q', 'V', ...
    'gW', 'gF', 'Y_grid', 'R_grid', 'beta', 'gamma', 'Rf', 'Rp', 'Rs_mean', ...
    'sigma_eps', 'sigma_s');

% 保存切片用于可视化
writematrix(C(:,1,:), 'result_3stages_PFI/C_slice.csv');
writematrix(A(:,1,:), 'result_3stages_PFI/A_slice.csv');
writematrix(Q(:,1,:), 'result_3stages_PFI/Q_slice.csv');
writematrix(V(:,1,:), 'result_3stages_PFI/V_slice.csv');
writematrix(gW, 'result_3stages_PFI/gW.csv');
writematrix(gF, 'result_3stages_PFI/gF.csv');

%% 解析解与数值解比较
fprintf('计算解析解并与数值解进行比较...\n');

% 使用上面已经计算的对数空间参数
% mu_s和r_f已经在上面初始化部分计算过了

%% 第2期解析解
% 根据公式(α₂ = (μₛ + 0.5σₛ²)/(γσₛ²))计算风险资产配置比例
alpha2_analytic = (mu_s + 0.5 * sigma_s^2 - r_f) / (gamma * sigma_s^2);

% 固定养老金购买比例（假设为外生给定）
eta2 = 0.1;  % 这里可以修改为其他合理的值

% 计算对数组合收益率的参数
log_portfolio_mean = alpha2_analytic * mu_s + (1-alpha2_analytic) * r_f;
log_portfolio_var = (alpha2_analytic * sigma_s)^2;

% 根据公式计算消费比例系数
kappa2_analytic = 1 / (1 + beta^(1/gamma) * exp((1-gamma)/gamma * (log_portfolio_mean + (1-gamma)/(2*gamma) * log_portfolio_var)));

fprintf('第2期解析解:\n');
fprintf('风险资产配置比例 (α₂): %.4f\n', alpha2_analytic);
fprintf('消费比例系数 (κ₂): %.4f\n', kappa2_analytic);

%% 第1期解析解
% 劳动收入与股票收益协方差（本例中假设为零，没有对冲需求）
sigma_sy = 0;  

% 根据公式(α₁ = α₂ - γσₛᵧ/σₛ²)计算风险资产配置比例
alpha1_analytic = alpha2_analytic - gamma * sigma_sy / sigma_s^2;

% 根据公式计算消费系数b
expected_portfolio_return = exp(log_portfolio_mean + 0.5 * log_portfolio_var);
b_coeff = 1 / (1 + beta^(1/gamma) * expected_portfolio_return^((1-gamma)/gamma));

% 固定养老金购买比例（假设为外生给定）
eta1 = 0.1;  % 这里可以修改为其他合理的值

% 简化消费比例计算（这是一个近似值）
kappa1_analytic = b_coeff / (1 - eta1);

fprintf('第1期解析解:\n');
fprintf('风险资产配置比例 (α₁): %.4f\n', alpha1_analytic);
fprintf('消费系数 (b): %.4f\n', b_coeff);
fprintf('消费比例系数 (κ₁): %.4f\n', kappa1_analytic);

%% 比较解析解与数值解
% 在现金网格上计算数值解的平均值（对养老金余额维度取第一个值）
A_numerical_p2 = mean(A(:,1,2));
C_numerical_p2 = mean(C(:,1,2));
A_numerical_p1 = mean(A(:,1,1));
C_numerical_p1 = mean(C(:,1,1));

fprintf('\n解析解与数值解比较:\n');
fprintf('第2期风险资产配置 - 解析解: %.4f, 数值解: %.4f, 差异: %.4f\n', ...
    alpha2_analytic, A_numerical_p2, abs(alpha2_analytic - A_numerical_p2));
fprintf('第2期消费比例 - 解析解: %.4f, 数值解: %.4f, 差异: %.4f\n', ...
    kappa2_analytic, C_numerical_p2, abs(kappa2_analytic - C_numerical_p2));
fprintf('第1期风险资产配置 - 解析解: %.4f, 数值解: %.4f, 差异: %.4f\n', ...
    alpha1_analytic, A_numerical_p1, abs(alpha1_analytic - A_numerical_p1));
fprintf('第1期消费比例 - 解析解: %.4f, 数值解: %.4f, 差异: %.4f\n', ...
    kappa1_analytic, C_numerical_p1, abs(kappa1_analytic - C_numerical_p1));

%% 可视化解析解与数值解的比较
% 创建图表
figure;
subplot(2,2,1);
bar([alpha2_analytic, A_numerical_p2]);
title('第2期风险资产配置比较');
set(gca, 'XTickLabel', {'解析解', '数值解'});
ylabel('风险资产比例');

subplot(2,2,2);
bar([kappa2_analytic, C_numerical_p2]);
title('第2期消费比例比较');
set(gca, 'XTickLabel', {'解析解', '数值解'});
ylabel('消费比例');

subplot(2,2,3);
bar([alpha1_analytic, A_numerical_p1]);
title('第1期风险资产配置比较');
set(gca, 'XTickLabel', {'解析解', '数值解'});
ylabel('风险资产比例');

subplot(2,2,4);
bar([kappa1_analytic, C_numerical_p1]);
title('第1期消费比例比较');
set(gca, 'XTickLabel', {'解析解', '数值解'});
ylabel('消费比例');

% 保存图表
saveas(gcf, 'result_3stages_PFI/analytic_vs_numerical.png');

% 绘制策略函数随财富变化的趋势图
figure;
% 第2期风险资产配置比例
subplot(2,2,1);
plot(gW, A(:,1,2), 'b-', 'LineWidth', 2);
hold on;
plot([gW(1), gW(end)], [alpha2_analytic, alpha2_analytic], 'r--', 'LineWidth', 2);
title('第2期风险资产配置');
xlabel('财富水平');
ylabel('风险资产比例');
legend('数值解', '解析解');
grid on;

% 第2期消费比例
subplot(2,2,2);
plot(gW, C(:,1,2), 'b-', 'LineWidth', 2);
hold on;
plot([gW(1), gW(end)], [kappa2_analytic, kappa2_analytic], 'r--', 'LineWidth', 2);
title('第2期消费比例');
xlabel('财富水平');
ylabel('消费比例');
legend('数值解', '解析解');
grid on;

% 第1期风险资产配置比例
subplot(2,2,3);
plot(gW, A(:,1,1), 'b-', 'LineWidth', 2);
hold on;
plot([gW(1), gW(end)], [alpha1_analytic, alpha1_analytic], 'r--', 'LineWidth', 2);
title('第1期风险资产配置');
xlabel('财富水平');
ylabel('风险资产比例');
legend('数值解', '解析解');
grid on;

% 第1期消费比例
subplot(2,2,4);
plot(gW, C(:,1,1), 'b-', 'LineWidth', 2);
hold on;
plot([gW(1), gW(end)], [kappa1_analytic, kappa1_analytic], 'r--', 'LineWidth', 2);
title('第1期消费比例');
xlabel('财富水平');
ylabel('消费比例');
legend('数值解', '解析解');
grid on;

% 保存图表
saveas(gcf, 'result_3stages_PFI/strategies_vs_wealth.png');

fprintf('解析解与数值解比较完成，结果已保存\n');

fprintf('三期生命周期模型求解完成，结果已保存到result_3stages_PFI文件夹\n');

%% 辅助函数
function [z_grid, z_prob] = tauchen_method(mu, rho, sigma, n, m)
    % Tauchen方法离散化AR(1)过程
    % mu: 均值
    % rho: 自相关系数
    % sigma: 创新项标准差
    % n: 网格点数
    % m: 状态空间范围倍数
    
    sigma_z = sigma / sqrt(1 - rho^2);  % 无条件标准差
    z_max = mu + m * sigma_z;
    z_min = mu - m * sigma_z;
    
    % 创建网格
    z_grid = linspace(z_min, z_max, n)';
    step = (z_max - z_min) / (n - 1);
    
    % 计算转移概率
    z_prob = zeros(n, 1);
    
    % 对于独立随机变量 (rho=0)，直接使用正态分布概率
    if rho == 0
        for i = 1:n
            if i == 1
                z_prob(i) = normcdf(z_grid(i) + step/2, mu, sigma);
            elseif i == n
                z_prob(i) = 1 - normcdf(z_grid(i) - step/2, mu, sigma);
            else
                z_prob(i) = normcdf(z_grid(i) + step/2, mu, sigma) - ...
                           normcdf(z_grid(i) - step/2, mu, sigma);
            end
        end
    end
    
    return;
end 