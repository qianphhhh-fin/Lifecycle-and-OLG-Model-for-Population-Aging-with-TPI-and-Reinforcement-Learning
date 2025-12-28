% =========================================================================
% solve_lifecycle_model.m
% 求解三期生命周期模型（无简化假设的完全数值解）
% =========================================================================

%% 1. 初始化
clear;
clc;
close all;

fprintf('开始求解三期生命周期模型...\n');

%% 2. 模型参数设定
% 将所有参数放入一个结构体中，方便传递
param.beta = 0.96;      % 时间偏好因子

% 收入
param.Y1 = 1.0;         % 第1期收入
param.Y2 = 1.2;         % 第2期收入

% 税收
param.tau = 0.2;        % 当期所得税
param.tau_p = 0.15;     % 退休期养老金税
param.q_max = 1;     % 养老金缴存比例上限

% 资产回报
param.Rf = 1.02;        % 无风险资产总回报
param.Rp = 1.03;        % 养老金账户总回报
param.RH = 1.15;        % 风险资产高回报
param.RL = 0.90;        % 风险资产低回报
param.p = 0.6;          % 高回报概率
% 检查风险溢价: E(Rs) > Rf
if (param.p * param.RH + (1-param.p) * param.RL <= param.Rf)
    warning('风险资产的期望回报不高于无风险回报，结果可能不符合直觉。');
end


%% 3. 数值方法设定
% 为第2期的状态变量构建网格
% 状态变量: (W2_begin, P1)
% W2_begin: 第2期初，由第1期储蓄S1带来的金融财富
% P1:       第1期末累积的养老金财富

% W2_begin 网格
nW2 = 80;               % W2网格点数量
max_S1 = (1-param.tau)*(1-0)*param.Y1; % S1的理论最大值
W2_grid_max = max_S1 * param.RH * 1.2; % 留一些余量
W2_grid = linspace(1e-5, W2_grid_max, nW2)';

% P1 网格
nP1 = 40;               % P1网格点数量
P1_grid_max = param.q_max * param.Y1 * 1.2;
P1_grid = linspace(0, P1_grid_max, nP1);

% 初始化结果矩阵
V2_mat = zeros(nW2, nP1);         % 价值函数
C2_opt_mat = zeros(nW2, nP1);     % 最优消费 C2*
q2_opt_mat = zeros(nW2, nP1);     % 最优缴存 q2*
alpha2_opt_mat = zeros(nW2, nP1); % 最优风险配置 alpha2*


%% 4. 逆向归纳：求解第2期
fprintf('正在求解第2期...\n');
tic;

% 使用并行计算加速（如果你的MATLAB有Parallel Computing Toolbox）
% parfor i = 1:nW2 
for i = 1:nW2
    % 内部循环不能并行，所以只并行外层
    for j = 1:nP1
        W2_begin = W2_grid(i);
        P1 = P1_grid(j);

        % 调用优化器求解该状态点的最优决策
        [solution, fval] = solve_period_2_problem(W2_begin, P1, param);
        
        % 存储结果
        V2_mat(i, j) = -fval; % fmincon求最小值，价值函数要取反
        C2_opt_mat(i, j) = solution(1);
        q2_opt_mat(i, j) = solution(2);
        alpha2_opt_mat(i, j) = solution(3);
    end
    if mod(i, 10) == 0
        fprintf('  已完成 %.0f%%...\n', (i/nW2)*100);
    end
end
toc;
fprintf('第2期求解完成。\n');


%% 5. 构建第2期价值函数的插值函数
% V2(W2_begin, P1)
V2_interp = griddedInterpolant({W2_grid, P1_grid}, V2_mat, "makima", "makima");


%% 6. 求解第1期
fprintf('正在求解第1期...\n');
tic;

% 目标函数
objective_t1 = @(x) objective_period_1(x, param, V2_interp);

% 决策变量 x = [C1, q1, alpha1]
% 初始猜测值
x0 = [0.6 * (1-param.tau)*param.Y1, 0.05, 0.5]; 
% 边界约束
% C1 > 0, 0 <= q1 <= q_max, 0 <= alpha1 <= 1
lb = [1e-6, 0, 0];
ub = [inf, param.q_max, 1];

% 非线性约束: 储蓄S1必须为正
% S1 = (1-tau)(1-q1)Y1 - C1 >= 0  => C1 - (1-tau)(1-q1)Y1 <= 0
nonlcon_t1 = @(x) savings_constraint_period_1(x, param);

% fmincon设置
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');

% 求解
[x1_opt, V1] = fmincon(objective_t1, x0, [], [], [], [], lb, ub, nonlcon_t1, options);

V1 = -V1;
C1_opt = x1_opt(1);
q1_opt = x1_opt(2);
alpha1_opt = x1_opt(3);

toc;
fprintf('第1期求解完成。\n\n');

%% 7. 结果展示
fprintf('==================== 最优决策结果 ====================\n');
fprintf('第1期最优决策:\n');
fprintf('  消费 C1*      = %.4f\n', C1_opt);
fprintf('  养老金缴存 q1* = %.4f (占收入 %.2f%%)\n', q1_opt, q1_opt*100);
fprintf('  风险资产配置 a1* = %.4f\n', alpha1_opt);
fprintf('======================================================\n');

%% 8. 绘图：可视化第2期的策略函数
figure('Name', '第2期策略函数 (Policy Functions)'); % , 'Position', [100, 100, 1200, 500]

% 创建网格用于绘图
[P1_mesh, W2_mesh] = meshgrid(P1_grid, W2_grid);

% 图1: 最优风险资产配置 alpha2*(W2_begin, P1)
subplot(1, 3, 1);
surf(P1_mesh, W2_mesh, alpha2_opt_mat);
shading interp;
xlabel('P_1');
ylabel('W_{2,begin}');
zlabel('\alpha_2^*');
title({'最优风险资产配置 \alpha_2^*', '(检验理论假设的关键图)'});
colorbar;
view(30, 25); % 调整视角

% 图2: 最优养老金缴存 q2*(W2_begin, P1)
subplot(1, 3, 2);
surf(P1_mesh, W2_mesh, q2_opt_mat);
shading interp;
xlabel('P_1');
ylabel('W_{2,begin}');
zlabel('q_2^*');
title('最优养老金缴存比例 q_2^*');
colorbar;
view(30, 25);

% 图3: 最优消费 C2*(W2_begin, P1)
subplot(1, 3, 3);
surf(P1_mesh, W2_mesh, C2_opt_mat);
shading interp;
xlabel('P_1');
ylabel('W_{2,begin}');
zlabel('C_2^*');
title('最优消费 C_2^*');
colorbar;
view(30, 25);

sgtitle('第2期最优决策作为状态变量 (W_{2,begin}, P_1) 的函数');

function [solution, fval] = solve_period_2_problem(W2_begin, P1, param)
% SOLVE_PERIOD_2_PROBLEM 使用fmincon求解第2期的优化问题
%
% 输入:
%   W2_begin - 第2期初的金融财富 (S1*Rs2)
%   P1       - 第1期末的养老金财富
%   param    - 参数结构体
%
% 输出:
%   solution - 最优决策 [C2*, q2*, alpha2*]
%   fval     - 目标函数的最小值 (-V2)

% 目标函数
objective = @(x) objective_period_2(x, W2_begin, P1, param);

% 决策变量 x = [C2, q2, alpha2]
% 初始猜测值
W2_available_guess = W2_begin + (1-param.tau)*(param.Y2 - 0.5*param.q_max*param.Y2);
x0 = [0.7 * W2_available_guess, 0.5*param.q_max, 0.5];

% 边界约束
% C2 > 0, 0 <= q2 <= q_max, 0 <= alpha2 <= 1
lb = [1e-6, 0, 0];
ub = [inf, param.q_max, 1];

% 非线性约束: 储蓄S2必须为正
% S2 = W2_begin + (1-tau)(Y2-q2*Y2) - C2 >= 0
% => C2 - (W2_begin + (1-tau)(Y2-q2*Y2)) <= 0
nonlcon = @(x) savings_constraint_period_2(x, W2_begin, param);

% fmincon 设置
options = optimoptions('fmincon', 'Display', 'none', 'Algorithm', 'sqp');

% 求解
[solution, fval] = fmincon(objective, x0, [], [], [], [], lb, ub, nonlcon, options);

end

% --- 内部辅助函数 ---

function neg_V2 = objective_period_2(x, W2_begin, P1, param)
    % 解包决策变量
    C2 = x(1);
    q2 = x(2);
    alpha2 = x(3);

    % 计算第2期可支配财富和储蓄
    W2_available = W2_begin + (1-param.tau)*(param.Y2 - q2*param.Y2);
    S2 = W2_available - C2;

    % 如果储蓄为负，返回一个巨大的惩罚值
    if S2 < 0
        neg_V2 = 1e10;
        return;
    end
    
    % 计算第2期末的养老金财富
    P2 = P1*param.Rp + q2*param.Y2;

    % 计算第3期消费的两种可能结果
    % C3 = S2 * R_opt3 + (1-tau_p)*P2*Rp
    R_opt3_H = (1-alpha2)*param.Rf + alpha2*param.RH;
    R_opt3_L = (1-alpha2)*param.Rf + alpha2*param.RL;
    
    C3_H = S2 * R_opt3_H + (1-param.tau_p)*P2*param.Rp;
    C3_L = S2 * R_opt3_L + (1-param.tau_p)*P2*param.Rp;
    
    % 如果C3可能为负（在高杠杆或极端回报下），也施加惩罚
    if C3_H <= 0 || C3_L <= 0
        neg_V2 = 1e10;
        return;
    end

    % 计算第2期开始的期望总效用
    E_log_C3 = param.p * log(C3_H) + (1-param.p) * log(C3_L);
    V2 = log(C2) + param.beta * E_log_C3;

    neg_V2 = -V2; % fmincon求最小值
end

function [c, ceq] = savings_constraint_period_2(x, W2_begin, param)
    C2 = x(1);
    q2 = x(2);
    
    % c(x) <= 0
    c = C2 - (W2_begin + (1-param.tau)*(param.Y2 - q2*param.Y2));
    % 无等式约束
    ceq = [];
end

function [c, ceq] = savings_constraint_period_1(x, param)
% SAVINGS_CONSTRAINT_PERIOD_1 第1期的储蓄非负约束

C1 = x(1);
q1 = x(2);

% c(x) <= 0
c = C1 - (1-param.tau)*(1-q1)*param.Y1;
% 无等式约束
ceq = [];
end

function neg_V1 = objective_period_1(x, param, V2_interp)
% OBJECTIVE_PERIOD_1 计算第1期的目标函数 (-V1)

% 解包决策变量
C1 = x(1);
q1 = x(2);
alpha1 = x(3);

% 计算第1期储蓄和养老金
S1 = (1-param.tau)*(1-q1)*param.Y1 - C1;
P1 = q1*param.Y1;

% 如果储蓄为负，返回巨大惩罚值
if S1 < 0
    neg_V1 = 1e10;
    return;
end

% 计算第2期初金融财富的两种可能结果
R_opt2_H = (1-alpha1)*param.Rf + alpha1*param.RH;
R_opt2_L = (1-alpha1)*param.Rf + alpha1*param.RL;

W2_begin_H = S1 * R_opt2_H;
W2_begin_L = S1 * R_opt2_L;

% 使用插值函数计算第2期价值函数的期望值
% 注意: V2_interp的输入是 {W2_grid, P1_grid}
V2_H = V2_interp(W2_begin_H, P1);
V2_L = V2_interp(W2_begin_L, P1);
E_V2 = param.p * V2_H + (1-param.p) * V2_L;

% 计算第1期的总价值函数
V1 = log(C1) + param.beta * E_V2;

neg_V1 = -V1; % fmincon求最小值
end