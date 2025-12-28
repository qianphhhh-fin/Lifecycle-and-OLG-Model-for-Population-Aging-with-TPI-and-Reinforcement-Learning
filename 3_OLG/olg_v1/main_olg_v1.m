%% main_olg.m
% 重叠世代模型(OLG)求解与模拟
% 使用值函数迭代(VFI)方法求解家庭最优决策问题
% 参考Python版本main_olg.py实现

% 清除工作区
clear;
close all;

% 创建保存图片和数据的目录
if ~exist('fig', 'dir')
    mkdir('fig');
end
if ~exist('data', 'dir')
    mkdir('data');
end

% 设置随机种子以保证结果可重复
rng(42);

% 运行模式选择
run_mode = 2;  % 1: 求解+模拟, 2: 仅模拟(从文件读取)

%% 模型参数初始化
% 模拟参数
periods = 100;  % 总模拟时期数
age_groups = 16;  % 16个五年期年龄组：22-101岁

% BGP检测参数
bgp_tolerance = 0.001;  % BGP稳定性容差
bgp_window = 5;  % 连续满足条件的窗口期

% 人口参数
x = -0.02;  % 新生代人口增长率调整因子

% 存活概率：从一个年龄组到下一个年龄组的存活概率
% 提高高龄组存活率以增加老龄人口
beta_surv = [0.995, 0.99, 0.985, 0.98, 0.975, 0.97, 0.965, 0.96, ...
    0.95, 0.94, 0.92, 0.89, 0.85, 0.80, 0.75];  % 显著提高高龄存活率

% 养老金系统参数
tau_g = 0.16;  % 公共养老金缴费率
lambda_pension = 0.45; % 养老金替代率

% 个人养老金参数
tau_p = 0.03;  % 个人养老金提取税率
phi = 0.005;  % 费率优惠因子
withdrawal_rate = 0.15;  % 退休期个人养老金提取率

% 生产函数参数
alpha = 0.36;  % 资本产出弹性
delta = 0.08;  % 资本折旧率
A = 1.0;  % 初始技术水平
growth_rate = 0.008;  % 技术进步率

% 劳动效率参数
% 0-7组(22-61岁)为工作人口，8-15组(62-101岁)为退休人口
h = [2.2, 2.5, 2.7, 2.9, 3.0, 3.1, 3.1, 3.0, ...  % 工作年龄组
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];   % 退休年龄组 (62-101岁)

% 效用函数参数
gamma = 2.0;  % 相对风险厌恶系数
theta = 0.04;  % 时间偏好率

% 市场参数
Rf = 1.02;  % 无风险收益率
Rp_mean = 1.04;  % 风险资产收益率期望
Rp_std = 0.18;  % 风险资产收益率标准差

% 政府参数
tau_c = 0.12;  % 消费税率
G_ratio = 0.18;  % 政府消费占GDP比例

% 退休参数
retirement_age = 8;  % 退休年龄组索引 (对应62岁开始)

% 值函数迭代参数
n_W = 50;  % 财富网格点数
n_P = 50;  % 个人养老金资产网格点数
min_W = 0.01;  % 最小财富水平
max_W = 10.0;  % 最大财富水平
min_P = 0.0;  % 最小个人养老金资产
max_P = 120.0;  % 最大个人养老金资产

% 风险资产收益率离散化
n_r = 5;  % 收益率离散化点数

% 优化问题参数
options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp', ...
    'MaxFunctionEvaluations', 1000, 'MaxIterations', 100, ...
    'OptimalityTolerance', 1e-6, 'StepTolerance', 1e-6);

%% 设置财富和养老金资产网格
% 对数均匀网格
W_grid = exp(linspace(log(min_W), log(max_W), n_W));
P_grid = linspace(min_P, max_P, n_P);

% 初始化值函数和策略函数
V = zeros(age_groups, n_W, n_P);  % 值函数
c_policy = zeros(age_groups, n_W, n_P);  % 消费率策略
alpha_policy = zeros(age_groups, n_W, n_P);  % 风险资产配置策略
q_policy = zeros(age_groups, n_W, n_P);  % 养老金缴费率策略

%% 离散化风险资产收益率分布
% 使用Gaussian quadrature方法
% 计算正态分布的分位数点
quantiles = linspace(0.1, 0.9, n_r);
r_grid = norminv(quantiles, Rp_mean, Rp_std);

% 计算各点的概率
r_prob = ones(n_r, 1) / n_r;

%% 主程序

if run_mode == 1
    % 运行值函数迭代
    fprintf('开始值函数迭代求解...\n');
    tic;
    
    % 终末期值函数 - 最后一期的最优策略是消费所有资产
    for i_w = 1:n_W
        for i_p = 1:n_P
            wealth = W_grid(i_w);
            pension = P_grid(i_p);
            
            % 最后一期提取所有养老金
            pension_benefit = pension;
            consumption = wealth + pension_benefit;
            
            % CRRA效用函数
            if gamma == 1
                V(age_groups, i_w, i_p) = log(consumption);
            else
                V(age_groups, i_w, i_p) = (consumption^(1-gamma) - 1) / (1-gamma);
            end
            
            c_policy(age_groups, i_w, i_p) = 1.0;  % 消费所有资产
            alpha_policy(age_groups, i_w, i_p) = 0.0;  % 不投资风险资产
            q_policy(age_groups, i_w, i_p) = 0.0;  % 不缴纳个人养老金
        end
    end
    
    % 从倒数第二期开始向前求解
    for age = (age_groups-1):-1:1
        fprintf('求解年龄组 %d 的值函数...\n', age);
        age_time = tic;
        
        % 是否为退休年龄
        is_retired = (age >= retirement_age);
        
        % 创建临时数组存储结果
        temp_c_policy = zeros(n_W, n_P);
        temp_alpha_policy = zeros(n_W, n_P);
        temp_q_policy = zeros(n_W, n_P);
        temp_V = zeros(n_W, n_P);
        
        if is_retired  % 退休期
            % 只在养老金维度上并行
            parfor i_p = 1:n_P
                % 为当前养老金水平创建临时数组
                local_c = zeros(n_W, 1);
                local_alpha = zeros(n_W, 1);
                local_q = zeros(n_W, 1);
                local_V = zeros(n_W, 1);
                
                pension = P_grid(i_p);
                % 养老金给付
                pension_benefit = pension * withdrawal_rate * (1 - tau_p);
                
                % 对每个财富水平求解
                for i_w = 1:n_W
                    wealth = W_grid(i_w);
                    
                    % 确定收入
                    income = lambda_pension + pension_benefit;
                    
                    % 可用资源
                    resources = wealth + income;
                    
                    % 优化问题
                    x0 = [0.5, 0.3];  % 初始猜测值 [消费率, 风险资产比例]
                    lb = [0.01, 0.0];  % 下界
                    ub = [0.99, 1.0];  % 上界
                    
                    % 优化函数
                    retired_obj = @(x) retired_value_function(x, age, wealth, pension, ...
                        squeeze(V(age+1,:,:)), r_grid, r_prob, Rf, lambda_pension, ...
                        tau_p, withdrawal_rate, gamma, theta, W_grid, P_grid, beta_surv(age));
                    
                    [policy, fval] = fmincon(retired_obj, x0, [], [], [], [], lb, ub, [], options);
                    
                    % 提取最优策略
                    c_rate = policy(1);
                    alpha = policy(2);
                    
                    % 存储到本地临时数组
                    local_c(i_w) = c_rate;
                    local_alpha(i_w) = alpha;
                    local_q(i_w) = 0.0;  % 退休期不缴纳养老金
                    local_V(i_w) = -fval;
                end
                
                % 将本地结果复制到临时结果数组
                temp_c_policy(:, i_p) = local_c;
                temp_alpha_policy(:, i_p) = local_alpha;
                temp_q_policy(:, i_p) = local_q;
                temp_V(:, i_p) = local_V;
            end
        else  % 工作期
            % 只在养老金维度上并行
            parfor i_p = 1:n_P
                % 为当前养老金水平创建临时数组
                local_c = zeros(n_W, 1);
                local_alpha = zeros(n_W, 1);
                local_q = zeros(n_W, 1);
                local_V = zeros(n_W, 1);
                
                pension = P_grid(i_p);
                
                % 对每个财富水平求解
                for i_w = 1:n_W
                    wealth = W_grid(i_w);
                    
                    % 劳动收入
                    labor_income = h(age);
                    
                    % 可用资源
                    resources = wealth + labor_income;
                    
                    % 优化问题
                    x0 = [0.5, 0.3, 0.1];  % 初始猜测值 [消费率, 风险资产比例, 养老金缴费率]
                    lb = [0.01, 0.0, 0.0];  % 下界
                    ub = [0.99, 1.0, 0.5];  % 上界
                    
                    % 优化函数
                    working_obj = @(x) working_value_function(x, age, wealth, pension, ...
                        squeeze(V(age+1,:,:)), r_grid, r_prob, Rf, Rp_mean, ...
                        labor_income, gamma, theta, W_grid, P_grid, beta_surv(age));
                    
                    [policy, fval] = fmincon(working_obj, x0, [], [], [], [], lb, ub, [], options);
                    
                    % 提取最优策略
                    c_rate = policy(1);
                    alpha = policy(2);
                    q = policy(3);
                    
                    % 存储到本地临时数组
                    local_c(i_w) = c_rate;
                    local_alpha(i_w) = alpha;
                    local_q(i_w) = q;
                    local_V(i_w) = -fval;
                end
                
                % 将本地结果复制到临时结果数组
                temp_c_policy(:, i_p) = local_c;
                temp_alpha_policy(:, i_p) = local_alpha;
                temp_q_policy(:, i_p) = local_q;
                temp_V(:, i_p) = local_V;
            end
        end
        
        % 将临时结果复制到主变量
        c_policy(age, :, :) = temp_c_policy;
        alpha_policy(age, :, :) = temp_alpha_policy;
        q_policy(age, :, :) = temp_q_policy;
        V(age, :, :) = temp_V;
        
        age_duration = toc(age_time);
        fprintf('  年龄组 %d 求解完成，耗时 %.2f 秒\n', age, age_duration);
    end
    
    vfi_time = toc;
    fprintf('值函数迭代求解完成，总耗时 %.2f 秒\n', vfi_time);
    
    % 保存求解结果
    fprintf('正在保存求解结果...\n');
    save_data = struct();
    save_data.V = V;
    save_data.c_policy = c_policy;
    save_data.alpha_policy = alpha_policy;
    save_data.q_policy = q_policy;
    save_data.W_grid = W_grid;
    save_data.P_grid = P_grid;
    save_data.r_grid = r_grid;
    save_data.r_prob = r_prob;
    save_data.params = struct();
    save_data.params.periods = periods;
    save_data.params.age_groups = age_groups;
    save_data.params.x = x;
    save_data.params.beta_surv = beta_surv;
    save_data.params.tau_g = tau_g;
    save_data.params.lambda_pension = lambda_pension;
    save_data.params.tau_p = tau_p;
    save_data.params.phi = phi;
    save_data.params.withdrawal_rate = withdrawal_rate;
    save_data.params.alpha = alpha;
    save_data.params.delta = delta;
    save_data.params.A = A;
    save_data.params.growth_rate = growth_rate;
    save_data.params.h = h;
    save_data.params.gamma = gamma;
    save_data.params.theta = theta;
    save_data.params.Rf = Rf;
    save_data.params.Rp_mean = Rp_mean;
    save_data.params.Rp_std = Rp_std;
    save_data.params.tau_c = tau_c;
    save_data.params.G_ratio = G_ratio;
    save_data.params.retirement_age = retirement_age;
    save_data.params.min_W = min_W;
    save_data.params.max_W = max_W;
    save_data.params.min_P = min_P;
    save_data.params.max_P = max_P;
    
    % 保存到文件
    save('data/vfi_results.mat', 'save_data');
    fprintf('求解结果已保存至 data/vfi_results.mat\n');

else
    % 从文件加载求解结果
    fprintf('从文件加载求解结果...\n');
    if ~exist('data/vfi_results.mat', 'file')
        error('求解结果文件不存在，请先运行求解模式 (run_mode=1)');
    end
    
    load('data/vfi_results.mat', 'save_data');
    
    % 恢复变量
    V = save_data.V;
    c_policy = save_data.c_policy;
    alpha_policy = save_data.alpha_policy;
    q_policy = save_data.q_policy;
    W_grid = save_data.W_grid;
    P_grid = save_data.P_grid;
    r_grid = save_data.r_grid;
    r_prob = save_data.r_prob;
    
    % 更新参数
    periods = save_data.params.periods;
    age_groups = save_data.params.age_groups;
    x = save_data.params.x;
    beta_surv = save_data.params.beta_surv;
    tau_g = save_data.params.tau_g;
    lambda_pension = save_data.params.lambda_pension;
    tau_p = save_data.params.tau_p;
    phi = save_data.params.phi;
    withdrawal_rate = save_data.params.withdrawal_rate;
    alpha = save_data.params.alpha;
    delta = save_data.params.delta;
    A = save_data.params.A;
    growth_rate = save_data.params.growth_rate;
    h = save_data.params.h;
    gamma = save_data.params.gamma;
    theta = save_data.params.theta;
    Rf = save_data.params.Rf;
    Rp_mean = save_data.params.Rp_mean;
    Rp_std = save_data.params.Rp_std;
    tau_c = save_data.params.tau_c;
    G_ratio = save_data.params.G_ratio;
    retirement_age = save_data.params.retirement_age;
    min_W = save_data.params.min_W;
    max_W = save_data.params.max_W;
    min_P = save_data.params.min_P;
    max_P = save_data.params.max_P;
    
    fprintf('求解结果加载完成\n');
end

% 运行蒙特卡洛模拟
fprintf('\n开始运行蒙特卡洛模拟 (100次)...\n');
sim_count = 10; % 模拟次数

% 存储所有模拟结果
% 初始化存储所有模拟结果的数组
all_r_history = zeros(periods, sim_count);           % 存储每次模拟的利率历史
all_w_history = zeros(periods, sim_count);           % 存储每次模拟的工资历史
all_Y_history = zeros(periods, sim_count);           % 存储每次模拟的产出历史
all_C_history = zeros(periods, sim_count);           % 存储每次模拟的消费历史
all_T_history = zeros(periods, sim_count);           % 存储每次模拟的税收历史
all_B_history = zeros(periods, sim_count);           % 存储每次模拟的养老金支出历史
all_deficit_history = zeros(periods, sim_count);     % 存储每次模拟的养老金赤字历史
all_replacement_rate_history = zeros(periods, sim_count); % 存储每次模拟的替代率历史
all_dependency_ratio_history = zeros(periods, sim_count); % 存储每次模拟的抚养比历史
all_K = zeros(periods+1, sim_count);                 % 存储每次模拟的资本存量历史（包括初始值）
all_P = zeros(periods+1, 16, sim_count);             % 存储每次模拟的个人养老金资产历史（按年龄组）
all_D = zeros(periods+1, sim_count);                 % 存储每次模拟的公共养老金累计赤字历史
all_S = zeros(periods+1, 16, sim_count);             % 存储每次模拟的私人金融资产历史（按年龄组）
all_alpha_history = zeros(periods, 16, sim_count);   % 存储每次模拟的风险投资比例历史（按年龄组）
all_c_rate_history = zeros(periods, 16, sim_count);  % 存储每次模拟的消费率历史（按年龄组）
all_q_rate_history = zeros(periods, 16, sim_count);  % 存储每次模拟的养老金缴费率历史（按年龄组）

% 执行多次模拟
for sim = 1:sim_count
    fprintf('运行模拟 %d/%d\n', sim, sim_count);
    
    % 使用run_simulation_with_random_returns而不是原来的run_simulation
    [r_history, w_history, Y_history, C_history, T_history, B_history, ...
     deficit_history, replacement_rate_history, dependency_ratio_history, ...
     Z, K, P, D, S, alpha_history, c_rate_history, q_rate_history] = run_simulation_with_random_returns(periods, c_policy, alpha_policy, q_policy, ...
     W_grid, P_grid, tau_g, lambda_pension, withdrawal_rate, tau_p, A, ...
     growth_rate, alpha, delta, Rf, Rp_mean, Rp_std, h, x, beta_surv, G_ratio, ...
     retirement_age, min_W, max_W, min_P, max_P);
    
    % 存储结果
    all_r_history(:, sim) = r_history;
    all_w_history(:, sim) = w_history;
    all_Y_history(:, sim) = Y_history;
    all_C_history(:, sim) = C_history;
    all_T_history(:, sim) = T_history;
    all_B_history(:, sim) = B_history;
    all_deficit_history(:, sim) = deficit_history;
    all_replacement_rate_history(:, sim) = replacement_rate_history;
    all_dependency_ratio_history(:, sim) = dependency_ratio_history;
    all_K(:, sim) = K;
    all_P(:,:,sim) = P;
    all_D(:, sim) = D;
    all_S(:,:,sim) = S;
    all_alpha_history(:,:,sim) = alpha_history;
    all_c_rate_history(:,:,sim) = c_rate_history;
    all_q_rate_history(:,:,sim) = q_rate_history;
end

% 计算平均结果
avg_r_history = mean(all_r_history, 2);
avg_w_history = mean(all_w_history, 2);
avg_Y_history = mean(all_Y_history, 2);
avg_C_history = mean(all_C_history, 2);
avg_T_history = mean(all_T_history, 2);
avg_B_history = mean(all_B_history, 2);
avg_deficit_history = mean(all_deficit_history, 2);
avg_replacement_rate_history = mean(all_replacement_rate_history, 2);
avg_dependency_ratio_history = mean(all_dependency_ratio_history, 2);
avg_K = mean(all_K, 2);
avg_P = mean(all_P, 3);
avg_D = mean(all_D, 2);
avg_S = mean(all_S, 3);
avg_alpha_history = mean(all_alpha_history, 3);
avg_c_rate_history = mean(all_c_rate_history, 3);
avg_q_rate_history = mean(all_q_rate_history, 3);

% 检测是否达到平衡增长路径(BGP)
fprintf('\n检测是否达到平衡增长路径(BGP)...\n');
[is_bgp, bgp_periods, bgp_metrics] = detect_balanced_growth_path(periods, avg_Y_history, avg_K(1:end-1), ...
    avg_C_history, avg_r_history, avg_w_history, Z, dependency_ratio_history, bgp_tolerance, bgp_window);

if is_bgp
    fprintf('系统已达到平衡增长路径！\n');
    fprintf('首次达到BGP的时期为: %d\n', bgp_periods(1));
    
    % 打印BGP特征
    fprintf('\nBGP特征:\n');
    fprintf('  - 最终经济增长率: %.4f%%\n', bgp_metrics.growth_rate * 100);
    fprintf('  - 资本产出比: %.4f\n', bgp_metrics.k_y_ratio);
    fprintf('  - 消费产出比: %.4f\n', bgp_metrics.c_y_ratio);
    fprintf('  - 老年抚养比: %.4f\n', bgp_metrics.dependency_ratio);
    
    % 绘制BGP指标
    plot_bgp_metrics(periods, is_bgp, bgp_periods, avg_Y_history, avg_K(1:end-1), ...
        avg_C_history, avg_r_history, avg_w_history, dependency_ratio_history);
else
    fprintf('系统在%d期内未达到平衡增长路径。\n', periods);
    fprintf('可能需要增加模拟期数或调整经济参数。\n');
end

% 绘制平均结果
% plot_monte_carlo_results(periods, avg_r_history, avg_w_history, avg_Y_history, avg_C_history, ...
%  avg_T_history, avg_B_history, avg_deficit_history, avg_replacement_rate_history, ...
%  avg_dependency_ratio_history, Z, avg_K, avg_P, avg_D, avg_S, c_policy, alpha_policy, ...
%  q_policy, W_grid, P_grid, retirement_age, all_deficit_history, all_Y_history);

% 绘制所有年龄组的风险投资比例
% plot_all_age_groups_alpha(periods, avg_alpha_history);

% 绘制所有年龄组的消费率
% plot_all_age_groups_consumption(periods, avg_c_rate_history);

% 绘制所有年龄组的养老金缴费率
% plot_all_age_groups_pension(periods, avg_q_rate_history, retirement_age);

% fprintf('蒙特卡洛模拟完成。结果已保存到fig文件夹\n');

%% 运行模拟主函数
function [r_history, w_history, Y_history, C_history, T_history, B_history, deficit_history,...
    replacement_rate_history, dependency_ratio_history, Z, K, P, D, S, alpha_history, c_rate_history, q_rate_history]...
    = run_simulation_with_random_returns(periods, c_policy, alpha_policy, q_policy, W_grid,...
    P_grid, tau_g, lambda_pension, withdrawal_rate, tau_p, A, growth_rate,...
    alpha, delta, Rf, Rp_mean, Rp_std, h, x, beta_surv, G_ratio, retirement_age, min_W, max_W, min_P, max_P)
    
    % 基本结构与run_simulation相同，但添加随机收益率
    [Z, K, P, D, S] = initialize_state();
    
    % 存储模拟结果
    r_history = zeros(periods, 1);
    w_history = zeros(periods, 1);
    Y_history = zeros(periods, 1);
    C_history = zeros(periods, 1);
    T_history = zeros(periods, 1);
    B_history = zeros(periods, 1);
    deficit_history = zeros(periods, 1);
    replacement_rate_history = zeros(periods, 1);
    dependency_ratio_history = zeros(periods, 1);
    alpha_history = zeros(periods, 16); % 记录每期每个年龄组的风险投资比例
    c_rate_history = zeros(periods, 16); % 记录每期每个年龄组的消费率
    q_rate_history = zeros(periods, 16); % 记录每期每个年龄组的养老金缴费率
    
    % 预分配完整向量
    Z = [Z; zeros(periods, 16)];
    K = [K; zeros(periods, 1)];
    P = [P; zeros(periods, 16)];
    D = [D; zeros(periods, 1)];
    S = [S; zeros(periods, 16)];
    
    % 运行模拟
    for t = 1:periods
        % 计算当期的生产和要素价格
        [Y, w, r, L] = calculate_factor_prices(K(t), Z(t,:), h, A, alpha, delta);
        
        % 更新技术水平
        A = A * (1 + growth_rate);
        
        % 公共养老金系统
        [T_g, B_g, deficit, pension_benefit] = pension_system(Z(t,:), w, tau_g, lambda_pension, h);
        
        % 使用带随机收益率的家庭优化决策函数
        [C, S_new, P_new, personal_pension_tax, alpha_values, c_rate_values, q_rate_values] = household_optimization_with_random_returns(t, w, r, S(t,:), P(t,:), Z(t,:),...
            c_policy, alpha_policy, q_policy, W_grid, P_grid, Rf, Rp_mean, Rp_std, lambda_pension, withdrawal_rate,...
            tau_p, h, retirement_age, min_W, max_W, min_P, max_P);
        
        % 记录每个年龄组的策略
        alpha_history(t, :) = alpha_values;
        c_rate_history(t, :) = c_rate_values;
        q_rate_history(t, :) = q_rate_values;
        
        % 资本市场出清
        [K_next, D_next] = capital_market_clearing(K(t), Y, C, deficit, personal_pension_tax, S_new, G_ratio, Rf, delta);
        
        % 更新状态变量
        K(t+1) = K_next;
        D(t+1) = D_next;
        S(t+1,:) = S_new;
        P(t+1,:) = P_new;
        
        % 更新人口
        if t < periods
            Z_new = update_population(Z, t, x, beta_surv);
            Z(t+1,:) = Z_new;
        end
        
        % 计算抚养比
        dependency_ratio = calculate_dependency_ratio(Z(t,:));
        
        % 记录历史数据
        r_history(t) = r;
        w_history(t) = w;
        Y_history(t) = Y;
        C_history(t) = C;
        T_history(t) = T_g;
        B_history(t) = B_g;
        deficit_history(t) = deficit;
        dependency_ratio_history(t) = dependency_ratio;
        
        % 计算总替代率（公共+个人）
        retired_pop = sum(Z(t,9:16));
        if retired_pop > 0
            avg_personal_pension = sum(P_new(9:16) .* Z(t,9:16)) / retired_pop;
        else
            avg_personal_pension = 0;
        end
        total_replacement_rate = (pension_benefit + avg_personal_pension) / w;
        replacement_rate_history(t) = total_replacement_rate;
    end
end

%% 带随机收益率的家庭优化函数
function [C_total, S_new, P_new, personal_pension_tax, alpha_values, c_rate_values, q_rate_values] = household_optimization_with_random_returns(t, w, r, S, P, Z,...
    c_policy, alpha_policy, q_policy, W_grid, P_grid, Rf, Rp_mean, Rp_std, lambda_pension, withdrawal_rate,...
    tau_p, h, retirement_age, min_W, max_W, min_P, max_P)
    
    % 大部分逻辑与原函数相同
    C_by_age = zeros(1, 16);
    S_new = zeros(1, 16);
    P_new = zeros(1, 16);
    personal_pension_taxes = zeros(1, 16);
    alpha_values = zeros(1, 16); % 存储每个年龄组的风险投资比例
    c_rate_values = zeros(1, 16); % 存储每个年龄组的消费率
    q_rate_values = zeros(1, 16); % 存储每个年龄组的养老金缴费率
    
    S = reshape(S, 1, []);
    P = reshape(P, 1, []);
    Z = reshape(Z, 1, []);
    
    % 对每个年龄组生成随机收益率
    R_risky = normrnd(Rp_mean, Rp_std, 1, 16);  % 每个年龄组的风险资产随机收益率
    
    % 对每个年龄组应用策略函数
    for a = 1:16
        % 获取当前年龄组的私人金融资产
        S_current = S(a);
        % 获取当前年龄组的个人养老金资产
        P_current = P(a);
        
        S_bounded = max(min(S_current, max_W), min_W);

        P_bounded = max(min(P_current, max_P), min_P);
        %  if P_current > 130
        %     disp([P_current])
        % end       
        c_rate = max(0, min(1, interpolate_policy(c_policy(a,:,:), S_bounded, P_bounded, W_grid, P_grid)));
        alpha = max(0, min(1, interpolate_policy(alpha_policy(a,:,:), S_bounded, P_bounded, W_grid, P_grid)));
        q_rate = max(0, min(1, interpolate_policy(q_policy(a,:,:), S_bounded, P_bounded, W_grid, P_grid)));

        
        
        % 记录策略值
        alpha_values(a) = alpha;
        c_rate_values(a) = c_rate;
        q_rate_values(a) = q_rate;
        
        if a <= retirement_age  % 工作人口
            labor_income = w * h(a);
            total_resources = S_current + labor_income;
            
            C_by_age(a) = c_rate * total_resources;
            savings = total_resources * (1 - c_rate);
            
            pension_contrib = savings * q_rate;
            remaining_savings = savings * (1 - q_rate);
            
            % 使用随机收益率进行资产更新
            R_portfolio = alpha * R_risky(a) + (1 - alpha) * Rf;
            S_new(a) = remaining_savings * R_portfolio;

            
            % 养老金资产也使用随机收益率
            P_new(a) = P_current * Rp_mean + pension_contrib;
            personal_pension_taxes(a) = 0;
            
        else  % 退休人口
            pension_withdrawal = P_current * withdrawal_rate;
            pension_benefit = pension_withdrawal * (1 - tau_p);
            personal_pension_taxes(a) = pension_withdrawal * tau_p;
            
            public_pension = lambda_pension * w;
            total_resources = S_current + public_pension + pension_benefit;
            
            C_by_age(a) = c_rate * total_resources;
            remaining_savings = total_resources * (1 - c_rate);
            
            % 使用随机收益率进行资产更新
            R_portfolio = alpha * R_risky(a) + (1 - alpha) * Rf;
            S_new(a) = remaining_savings * R_portfolio;       
            P_new(a) = P_current * (1 - withdrawal_rate);
        end
    end
    
    C_total = sum(C_by_age .* Z);
    personal_pension_tax = sum(personal_pension_taxes .* Z);
end

%% 绘制所有年龄组的风险投资比例
function plot_all_age_groups_alpha(periods, alpha_history)
    % 设置中文显示支持
    set(0, 'DefaultAxesFontName', 'SimHei');
    set(0, 'DefaultTextFontName', 'SimHei');
    
    % 创建4x4的子图
    figure('Position', [50, 50, 1500, 1000]);
    
    for age = 1:16
        subplot(4, 4, age);
        plot(1:periods, alpha_history(:, age), 'LineWidth', 2);
        title(sprintf('年龄组 %d (%d-%d岁)', age, 22+(age-1)*5, 26+(age-1)*5), 'FontSize', 12);
        xlabel('时期', 'FontSize', 10);
        ylabel('风险投资比例', 'FontSize', 10);
        ylim([0, 1]);
        grid on;
    end
    
    sgtitle('各年龄组风险投资比例随时间变化', 'FontSize', 16);
    
    % 调整图表布局
    set(gcf, 'PaperPositionMode', 'auto');
    tight_layout = get(gcf, 'Position');
    set(gcf, 'Position', [tight_layout(1), tight_layout(2), tight_layout(3), tight_layout(4)]);
    
    % 保存图表
    saveas(gcf, 'fig/all_age_groups_alpha.png');
end

%% 绘制所有年龄组的消费率
function plot_all_age_groups_consumption(periods, c_rate_history)
    % 设置中文显示支持
    set(0, 'DefaultAxesFontName', 'SimHei');
    set(0, 'DefaultTextFontName', 'SimHei');
    
    % 创建4x4的子图
    figure('Position', [50, 50, 1500, 1000]);
    
    for age = 1:16
        subplot(4, 4, age);
        plot(1:periods, c_rate_history(:, age), 'LineWidth', 2, 'Color', [0, 0.5, 0]);
        title(sprintf('年龄组 %d (%d-%d岁)', age, 22+(age-1)*5, 26+(age-1)*5), 'FontSize', 12);
        xlabel('时期', 'FontSize', 10);
        ylabel('消费率', 'FontSize', 10);
        ylim([0, 1]);
        grid on;
    end
    
    sgtitle('各年龄组消费率随时间变化', 'FontSize', 16);
    
    % 调整图表布局
    set(gcf, 'PaperPositionMode', 'auto');
    tight_layout = get(gcf, 'Position');
    set(gcf, 'Position', [tight_layout(1), tight_layout(2), tight_layout(3), tight_layout(4)]);
    
    % 保存图表
    saveas(gcf, 'fig/all_age_groups_consumption.png');
end

%% 绘制所有年龄组的养老金缴费率
function plot_all_age_groups_pension(periods, q_rate_history, retirement_age)
    % 设置中文显示支持
    set(0, 'DefaultAxesFontName', 'SimHei');
    set(0, 'DefaultTextFontName', 'SimHei');
    
    % 创建4x4的子图
    figure('Position', [50, 50, 1500, 1000]);
    
    for age = 1:16
        subplot(4, 4, age);
        plot(1:periods, q_rate_history(:, age), 'LineWidth', 2, 'Color', [0.8, 0.2, 0]);
        title(sprintf('年龄组 %d (%d-%d岁)', age, 22+(age-1)*5, 26+(age-1)*5), 'FontSize', 12);
        xlabel('时期', 'FontSize', 10);
        ylabel('养老金缴费率', 'FontSize', 10);
        
        if age <= retirement_age
            ylim([0, 0.6]);  % 工作年龄组
        else
            % ylim([0, 0.1]);  % 退休年龄组（应该接近0）
            text(periods/2, 0.05, '退休期', 'HorizontalAlignment', 'center', 'FontSize', 12, 'Color', [0.5 0.5 0.5]);
        end
        
        grid on;
    end
    
    sgtitle('各年龄组养老金缴费率随时间变化', 'FontSize', 16);
    
    % 调整图表布局
    set(gcf, 'PaperPositionMode', 'auto');
    tight_layout = get(gcf, 'Position');
    set(gcf, 'Position', [tight_layout(1), tight_layout(2), tight_layout(3), tight_layout(4)]);
    
    % 保存图表
    saveas(gcf, 'fig/all_age_groups_pension.png');
end

%% 绘制蒙特卡洛模拟结果，增加置信区间
function plot_monte_carlo_results(periods, r_history, w_history, Y_history, C_history, ...
 T_history, B_history, deficit_history, replacement_rate_history, ...
 dependency_ratio_history, Z, K, P, D, S, c_policy, alpha_policy, ...
 q_policy, W_grid, P_grid, retirement_age, all_deficit_history, all_Y_history)
    
    % 基本绘图逻辑与原函数相同，但增加置信区间
    % 设置中文显示支持
    set(0, 'DefaultAxesFontName', 'SimHei');
    set(0, 'DefaultTextFontName', 'SimHei');
    
    % 人口结构计算
    young_workers = sum(Z(1:periods, 1:4), 2);  % 22-41岁
    older_workers = sum(Z(1:periods, 5:8), 2);  % 42-61岁
    retirees = sum(Z(1:periods, 9:16), 2);  % 62-101岁
    total_population = young_workers + older_workers + retirees;
    
    young_ratio = young_workers ./ total_population * 100;
    older_ratio = older_workers ./ total_population * 100;
    retirees_ratio = retirees ./ total_population * 100;
    
    % 创建图表 - 添加95%置信区间
    figure('Position', [100, 100, 1200, 900]);
    
    % 养老金赤字主图 - 带置信区间
    subplot(2, 2, 1);
    
    % 计算赤字的95%置信区间
    deficit_lower = prctile(all_deficit_history, 2.5, 2);
    deficit_upper = prctile(all_deficit_history, 97.5, 2);
    
    % 绘制平均赤字和置信区间
    plot(1:periods, deficit_history, 'r-', 'LineWidth', 2.5);
    hold on;
    fill([1:periods, fliplr(1:periods)], [deficit_lower', fliplr(deficit_upper')], 'r', 'FaceAlpha', 0.2);
    
    title('养老金赤字演变 (带95%置信区间)', 'FontSize', 14);
    xlabel('时期', 'FontSize', 12);
    ylabel('赤字金额', 'FontSize', 12);
    legend('平均赤字', '95%置信区间', 'Location', 'best', 'FontSize', 12);
    grid on;
    
    % 其他图表与原函数类似，可以添加更多置信区间
    % ...

    saveas(gcf, 'fig/monte_carlo_results.png');
end

%% 退休期值函数计算函数
function obj_value = retired_value_function(x, age, wealth, pension, next_V, r_grid, r_prob, Rf, lambda_pension, tau_p, withdrawal_rate, gamma, theta, W_grid, P_grid, beta_surv)
    % 提取决策变量
    c_rate = x(1);
    alpha = x(2);
    
    % 养老金给付
    pension_benefit = pension * withdrawal_rate * (1 - tau_p);
    
    % 确定收入
    income = lambda_pension + pension_benefit;
    
    % 当期消费
    consumption = (wealth + income) * c_rate;
    
    % 当期效用
    if gamma == 1
        current_utility = log(consumption);
    else
        current_utility = (consumption^(1-gamma) - 1) / (1-gamma);
    end
    
    % 下期资产
    savings = (wealth + income) * (1 - c_rate);
    
    % 期望折现效用
    expected_value = 0;
    
    % 计算继续价值
    for i_r = 1:length(r_grid)
        % 下期风险资产收益率
        R_risky = r_grid(i_r);
        
        % 下期无风险利率
        R_safe = Rf;
        
        % 投资组合收益率
        R_portfolio = alpha * R_risky + (1 - alpha) * R_safe;
        
        % 下期财富
        next_wealth = savings * R_portfolio;
        
        % 限制在网格范围内
        next_wealth = max(min(next_wealth, W_grid(end)), W_grid(1));
        
        % 下期养老金资产（扣除提取部分）
        next_pension = pension * (1 - withdrawal_rate);
        next_pension = max(min(next_pension, P_grid(end)), P_grid(1));
        
        % 插值得到下期值函数
        next_value = interpolate_value(next_wealth, next_pension, next_V, W_grid, P_grid);
        
        % 累加期望值
        expected_value = expected_value + r_prob(i_r) * next_value;
    end
    
    % 乘以折现因子和存活概率
    expected_value = expected_value * (1 / (1 + theta)) * beta_surv;
    
    % 返回负的值函数（最小化问题）
    obj_value = -(current_utility + expected_value);
end

%% 工作期值函数计算函数
function obj_value = working_value_function(x, age, wealth, pension, next_V, r_grid, r_prob, Rf, Rp_mean, labor_income, gamma, theta, W_grid, P_grid, beta_surv)
    % 提取决策变量
    c_rate = x(1);
    alpha = x(2);
    q = x(3);
    
    % 资源约束
    resources = wealth + labor_income;
    
    % 当期消费
    consumption = resources * c_rate;
    
    % 储蓄
    savings = resources * (1 - c_rate);
    
    % 个人养老金缴费
    pension_contrib = savings * q;
    
    % 剩余资产投资
    remaining_savings = savings * (1 - q);
    
    % 当期效用
    if gamma == 1
        current_utility = log(consumption);
    else
        current_utility = (consumption^(1-gamma) - 1) / (1-gamma);
    end
    
    % 期望折现效用
    expected_value = 0;
    
    % 计算继续价值
    for i_r = 1:length(r_grid)
        % 下期风险资产收益率
        R_risky = r_grid(i_r);
        
        % 下期无风险利率
        R_safe = Rf;
        
        % 投资组合收益率
        R_portfolio = alpha * R_risky + (1 - alpha) * R_safe;
        
        % 下期财富
        next_wealth = remaining_savings * R_portfolio;
        
        % 限制在网格范围内
        next_wealth = max(min(next_wealth, W_grid(end)), W_grid(1));
        
        % 下期养老金资产（包括本期缴费和收益）
        next_pension = (pension + pension_contrib) * Rp_mean;
        next_pension = max(min(next_pension, P_grid(end)), P_grid(1));
        
        % 插值得到下期值函数
        next_value = interpolate_value(next_wealth, next_pension, next_V, W_grid, P_grid);
        
        % 累加期望值
        expected_value = expected_value + r_prob(i_r) * next_value;
    end
    
    % 乘以折现因子和存活概率
    expected_value = expected_value * (1 / (1 + theta)) * beta_surv;
    
    % 返回负的值函数（最小化问题）
    obj_value = -(current_utility + expected_value);
end

%% 插值函数
function v = interpolate_value(wealth, pension, V, W_grid, P_grid)
    % 使用二维插值计算值函数
    v = interp2(P_grid, W_grid, V, pension, wealth, 'spline');
end

%% 初始化状态变量
function [Z, K, P, D, S] = initialize_state()
    % 初始化状态变量 - 使用2023年中国人口结构数据
    initial_pop = [76.20905535, 86.45596319, 113.8702459, 98.60198303, 86.64117824, 102.7890433, 112.0217783, 99.04620047, ...
                  64.05142331, 66.93157492, 44.16815149, 25.40848066, 14.97325553, 6.872421945, 1.743059943, 0.216184341];
    
    Z = zeros(1, 16);  % 人口分布，行为时期，列为年龄组
    Z(1,:) = initial_pop;
    
    K = zeros(1, 1);
    K(1) = 3000.0;  % 初始资本存量
    
    % 个人养老金资产分布 - 16个年龄组，降低初始积累
    initial_pension = zeros(1, 16);
    for i = 1:8  % 工作年龄组1-8
        initial_pension(i) = 5.0 * i;  % 初始积累
    end
    
    P = zeros(1, 16);
    P(1,:) = initial_pension;  % 个人养老金资产
    
    D = zeros(1, 1);
    D(1) = 5.0;  % 初始公共养老金累计赤字
    
    % 初始化私人金融资产
    initial_saving = zeros(1, 16);
    for i = 1:16
        initial_saving(i) = 0.1 * i;  % 简单假设初始储蓄水平
    end
    
    S = zeros(1, 16);
    S(1,:) = initial_saving;  % 私人金融资产
end

%% 更新人口分布
function Z_new = update_population(Z, t, x, beta_surv)
    % 更新人口分布 - 16个年龄组
    Z_new = zeros(1, 16);
    
    % 最年轻组 (22-26岁) - 新生人口
    % 随时间变化的增长率：前期小负增长，后期大负增长
    growth_rate = x;
    if t < 6
        % 前6期较缓和的负增长
        growth_rate = -0.01 - 0.003 * t;  % 从-1%逐渐增大负增长率
    elseif t >= 6
        % 后期强化负增长
        growth_rate = -0.03 - 0.004 * min(t-6, 10);  % 最终达到-7%的负增长率
    end
    
    Z_new(1) = Z(t, 1) * (1 + growth_rate);
    
    % 其他年龄组 - 从上一组存活下来的人口
    for a = 2:16
        Z_new(a) = Z(t, a-1) * beta_surv(a-1);
    end
end

%% 生产函数
function Y = production_function(K, L, A, alpha)
    % Cobb-Douglas 生产函数
    Y = K^alpha * (A * L)^(1-alpha);
end

%% 计算要素价格
function [Y, w, r, L] = calculate_factor_prices(K, Z, h, A, alpha, delta)
    % 计算有效劳动力 - 仅考虑工作年龄组 (1-8)
    L = 0;
    for a = 1:8  % 1-8组为工作人口
        L = L + h(a) * Z(a);
    end
    
    % 计算总产出
    Y = production_function(K, L, A, alpha);
    
    % 要素价格：工资和利率
    w = (1 - alpha) * Y / L;
    r = alpha * Y / K - delta;
end

%% 插值策略函数
function policy = interpolate_policy(policy_grid, wealth, pension, W_grid, P_grid)
    % 使用二维插值计算策略函数
    policy_grid = squeeze(policy_grid);
    policy = interp2(P_grid, W_grid, policy_grid, pension, wealth, 'spline');
end

%% 养老金系统
function [T_g, B_g, deficit, pension_benefit] = pension_system(Z, w, tau_g, lambda_pension, h)
    % 公共养老金收入 - 来自工作人口 (1-8组)
    T_g = 0;
    for a = 1:8  % 工作年龄组
        T_g = T_g + tau_g * w * h(a) * Z(a);
    end
    
    % 退休时的工资率（简化为当期工资）
    retirement_wage = w;
    
    % 公共养老金支出 - 支付给退休人口 (9-16组)
    pension_benefit = lambda_pension * retirement_wage;
    B_g = 0;
    for a = 9:16  % 退休年龄组
        B_g = B_g + pension_benefit * Z(a);
    end
    
    % 养老金赤字
    deficit = max(0, B_g - T_g);
end

%% 资本市场出清
function [K_next, D_next] = capital_market_clearing(K, Y, C, deficit, personal_pension_tax, S_new, G_ratio, Rf, delta)
    % 政府消费
    G = G_ratio * Y;
    
    % 投资
    I = Y - C - G;
    
    % 更新资本存量
    K_next = (1 - delta) * K + I;
    
    % 更新养老金累计赤字
    D_next = K * (1 + Rf - 1) + deficit;
end

%% 计算抚养比
function dependency_ratio = calculate_dependency_ratio(Z)
    % 计算抚养比：退休人口/工作人口
    retired_population = sum(Z(9:16));
    working_population = sum(Z(1:8));
    
    if working_population > 0
        dependency_ratio = retired_population / working_population;
    else
        dependency_ratio = inf;  % 避免除零错误
    end
end

%% 检测平衡增长路径(BGP)函数
function [is_bgp, bgp_periods, bgp_metrics] = detect_balanced_growth_path(periods, Y_history, K_history, ...
    C_history, r_history, w_history, Z, dependency_ratio_history, tolerance, window)
    % 初始化
    is_bgp = false;
    bgp_periods = [];
    bgp_metrics = struct();
    
    % 计算各关键变量的增长率
    Y_growth = zeros(periods-1, 1);
    K_growth = zeros(periods-1, 1);
    C_growth = zeros(periods-1, 1);
    w_growth = zeros(periods-1, 1);
    
    for t = 2:periods
        Y_growth(t-1) = Y_history(t) / Y_history(t-1) - 1;
        K_growth(t-1) = K_history(t) / K_history(t-1) - 1;
        C_growth(t-1) = C_history(t) / C_history(t-1) - 1;
        w_growth(t-1) = w_history(t) / w_history(t-1) - 1;
    end
    
    % 计算人口年龄结构比例
    pop_structure = zeros(periods, 3);
    for t = 1:periods
        total_pop = sum(Z(t,:));
        % 年轻工人比例 (22-41岁)
        pop_structure(t,1) = sum(Z(t,1:4)) / total_pop;
        % 年长工人比例 (42-61岁)
        pop_structure(t,2) = sum(Z(t,5:8)) / total_pop;
        % 退休人口比例 (62+岁)
        pop_structure(t,3) = sum(Z(t,9:16)) / total_pop;
    end
    
    % 计算人口结构变化率
    pop_structure_change = zeros(periods-1, 3);
    for t = 2:periods
        for i = 1:3
            pop_structure_change(t-1,i) = abs(pop_structure(t,i) - pop_structure(t-1,i));
        end
    end
    
    % 连续稳定期检测
    stable_periods = 0;
    
    % 至少从第window+1期开始检测，以便有足够的数据计算增长率
    for t = window+1:periods
        % 确保所有索引都是正整数
        if t-window < 1
            continue;
        end
        
        % 检查过去window期内的稳定性
        
        % 1. 检查增长率稳定性
        Y_growth_stable = all(abs(diff(Y_growth(t-window:t-1))) < tolerance);
        K_growth_stable = all(abs(diff(K_growth(t-window:t-1))) < tolerance);
        C_growth_stable = all(abs(diff(C_growth(t-window:t-1))) < tolerance);
        
        % 2. 检查利率和工资稳定性
        r_stable = all(abs(diff(r_history(t-window:t-1))) < tolerance);
        % 使用预先计算的工资增长率，避免索引出错
        w_growth_stable = all(abs(diff(w_growth(t-window:t-1))) < tolerance);
        
        % 3. 检查人口结构稳定性
        pop_structure_stable = all(all(pop_structure_change(t-window:t-1,:) < tolerance));
        
        % 4. 检查老年抚养比稳定性
        dependency_ratio_stable = all(abs(diff(dependency_ratio_history(t-window:t-1))) < tolerance);
        
        % 所有条件都满足
        if Y_growth_stable && K_growth_stable && C_growth_stable && ...
           r_stable && w_growth_stable && pop_structure_stable && dependency_ratio_stable
            
            % 记录首次达到BGP的时期
            if isempty(bgp_periods)
                bgp_periods = t-window:t-1;
                
                % 记录BGP特征
                bgp_metrics.growth_rate = mean(Y_growth(t-window:t-1));
                bgp_metrics.k_y_ratio = mean(K_history(t-window:t-1) ./ Y_history(t-window:t-1));
                bgp_metrics.c_y_ratio = mean(C_history(t-window:t-1) ./ Y_history(t-window:t-1));
                bgp_metrics.dependency_ratio = mean(dependency_ratio_history(t-window:t-1));
            end
            
            stable_periods = stable_periods + 1;
        else
            % 重置稳定期计数
            stable_periods = 0;
        end
    end
    
    % 如果有任何稳定期，认为达到了BGP
    is_bgp = ~isempty(bgp_periods);
end

%% 绘制BGP指标图
function plot_bgp_metrics(periods, is_bgp, bgp_periods, Y_history, K_history, ...
    C_history, r_history, w_history, dependency_ratio_history)
    
    % 如果未达到BGP，不绘制
    if ~is_bgp
        return;
    end
    
    % 设置中文显示支持
    set(0, 'DefaultAxesFontName', 'SimHei');
    set(0, 'DefaultTextFontName', 'SimHei');
    
    % 计算增长率
    Y_growth = zeros(periods-1, 1);
    K_growth = zeros(periods-1, 1);
    C_growth = zeros(periods-1, 1);
    
    for t = 2:periods
        Y_growth(t-1) = Y_history(t) / Y_history(t-1) - 1;
        K_growth(t-1) = K_history(t) / K_history(t-1) - 1;
        C_growth(t-1) = C_history(t) / C_history(t-1) - 1;
    end
    
    % 创建图表
    figure('Position', [150, 150, 1200, 900]);
    
    % 创建六个子图
    % 1. 产出、资本和消费增长率
    subplot(2, 3, 1);
    plot(2:periods, Y_growth, 'b-', 'LineWidth', 2);
    hold on;
    plot(2:periods, K_growth, 'r-', 'LineWidth', 2);
    plot(2:periods, C_growth, 'g-', 'LineWidth', 2);
    if ~isempty(bgp_periods)
        % 标记BGP区域
        x_bgp = [min(bgp_periods), max(bgp_periods), max(bgp_periods), min(bgp_periods)];
        y_limits = ylim;
        patch([x_bgp(1), x_bgp(2), x_bgp(2), x_bgp(1)], [y_limits(1), y_limits(1), y_limits(2), y_limits(2)], ...
            [0.8, 0.8, 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        text(min(bgp_periods), y_limits(1) + 0.9 * (y_limits(2) - y_limits(1)), 'BGP', ...
            'FontSize', 12, 'HorizontalAlignment', 'center');
    end
    title('关键变量增长率', 'FontSize', 14);
    xlabel('时期', 'FontSize', 12);
    ylabel('增长率', 'FontSize', 12);
    legend('产出', '资本', '消费', 'Location', 'best', 'FontSize', 10);
    grid on;
    
    % 2. 利率
    subplot(2, 3, 2);
    plot(1:periods, r_history, 'b-', 'LineWidth', 2);
    if ~isempty(bgp_periods)
        % 标记BGP区域
        x_bgp = [min(bgp_periods), max(bgp_periods), max(bgp_periods), min(bgp_periods)];
        y_limits = ylim;
        patch([x_bgp(1), x_bgp(2), x_bgp(2), x_bgp(1)], [y_limits(1), y_limits(1), y_limits(2), y_limits(2)], ...
            [0.8, 0.8, 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    end
    title('利率', 'FontSize', 14);
    xlabel('时期', 'FontSize', 12);
    ylabel('利率', 'FontSize', 12);
    grid on;
    
    % 3. 工资
    subplot(2, 3, 3);
    plot(1:periods, w_history, 'r-', 'LineWidth', 2);
    if ~isempty(bgp_periods)
        % 标记BGP区域
        x_bgp = [min(bgp_periods), max(bgp_periods), max(bgp_periods), min(bgp_periods)];
        y_limits = ylim;
        patch([x_bgp(1), x_bgp(2), x_bgp(2), x_bgp(1)], [y_limits(1), y_limits(1), y_limits(2), y_limits(2)], ...
            [0.8, 0.8, 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    end
    title('工资率', 'FontSize', 14);
    xlabel('时期', 'FontSize', 12);
    ylabel('工资', 'FontSize', 12);
    grid on;
    
    % 4. 资本产出比
    subplot(2, 3, 4);
    K_Y_ratio = K_history ./ Y_history;
    plot(1:periods, K_Y_ratio, 'b-', 'LineWidth', 2);
    if ~isempty(bgp_periods)
        % 标记BGP区域
        x_bgp = [min(bgp_periods), max(bgp_periods), max(bgp_periods), min(bgp_periods)];
        y_limits = ylim;
        patch([x_bgp(1), x_bgp(2), x_bgp(2), x_bgp(1)], [y_limits(1), y_limits(1), y_limits(2), y_limits(2)], ...
            [0.8, 0.8, 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    end
    title('资本产出比', 'FontSize', 14);
    xlabel('时期', 'FontSize', 12);
    ylabel('K/Y', 'FontSize', 12);
    grid on;
    
    % 5. 消费产出比
    subplot(2, 3, 5);
    C_Y_ratio = C_history ./ Y_history;
    plot(1:periods, C_Y_ratio, 'g-', 'LineWidth', 2);
    if ~isempty(bgp_periods)
        % 标记BGP区域
        x_bgp = [min(bgp_periods), max(bgp_periods), max(bgp_periods), min(bgp_periods)];
        y_limits = ylim;
        patch([x_bgp(1), x_bgp(2), x_bgp(2), x_bgp(1)], [y_limits(1), y_limits(1), y_limits(2), y_limits(2)], ...
            [0.8, 0.8, 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    end
    title('消费产出比', 'FontSize', 14);
    xlabel('时期', 'FontSize', 12);
    ylabel('C/Y', 'FontSize', 12);
    grid on;
    
    % 6. 老年抚养比
    subplot(2, 3, 6);
    plot(1:periods, dependency_ratio_history, 'r-', 'LineWidth', 2);
    if ~isempty(bgp_periods)
        % 标记BGP区域
        x_bgp = [min(bgp_periods), max(bgp_periods), max(bgp_periods), min(bgp_periods)];
        y_limits = ylim;
        patch([x_bgp(1), x_bgp(2), x_bgp(2), x_bgp(1)], [y_limits(1), y_limits(1), y_limits(2), y_limits(2)], ...
            [0.8, 0.8, 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    end
    title('老年抚养比', 'FontSize', 14);
    xlabel('时期', 'FontSize', 12);
    ylabel('抚养比', 'FontSize', 12);
    grid on;
    
    sgtitle('平衡增长路径(BGP)特征', 'FontSize', 16);
    
    % 调整图表布局
    set(gcf, 'PaperPositionMode', 'auto');
    tight_layout = get(gcf, 'Position');
    set(gcf, 'Position', [tight_layout(1), tight_layout(2), tight_layout(3), tight_layout(4)]);
    
    % 保存图表
    saveas(gcf, 'fig/bgp_metrics.png');
end

