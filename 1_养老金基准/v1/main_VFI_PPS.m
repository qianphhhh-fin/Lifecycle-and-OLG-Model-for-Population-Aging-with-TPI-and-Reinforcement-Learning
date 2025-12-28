%% main_VFI_PPS.m
% 生命周期模型求解与模拟（值函数迭代方法）
% 严格遵循LaTeX理论文档描述
% 运行模式：求解+模拟

% 清除工作区
clear;
close all;
fclose('all'); % 关闭所有打开的文件

direct_simu = false;
save_path = 'result_pps/';
if ~exist(save_path, 'dir')
    mkdir(save_path);
end
diary(fullfile(save_path, 'log.txt')); % 将命令窗口输出保存到日志文件
diary on;

%% 模型求解和模拟
if ~direct_simu

% 设置随机种子以获得可重复的结果
rng(42);

% --- 1. 参数设定 ---
fprintf('1. 开始设置模型参数...\n');

% 人口结构与生命周期
tb = 25;    % 进入模型的年龄 (t=1 in theory)
tr = 65;    % 退休年龄 (K in theory)
td = 100;   % 死亡年龄 (T in theory)
tn = td - tb + 1;  % 总期数
work_period = tr - tb; % 工作期长度

% 偏好
beta = 0.96;    % 主观折现因子
gamma = 3.0;    % 相对风险规避系数 (1/psi)

% 收入过程 (遵循理论)
% ln Y_t = f(t) + M_t + u_t, M_t = M_{t-1} + z_t
% f(t) = b1*age + b2*age^2 (Cocco et al., 2005 使用的形式)
b1 = 0.046;
b2 = -0.00047;
sigma_u = sqrt(0.025); % 暂时性冲击标准差
sigma_z = sqrt(0.012); % 持久性冲击标准差

% 养老金体系与税收 (新参数)
Y_bar = 0.2;         % 社会基本养老金 (第一支柱)
tau_y = 0.25;        % 劳动所得税率
tau_q = 0.15;        % 养老金提取税率
Q_max = 0.18;        % 养老金缴费占收入的上限比例 (简化为收入比例)

% 资产回报率
R_f = 1.02;     % 无风险资产税后总回报率
R_p = 1.02;     % 个人养老金账户回报率 (基准设为无风险)
mu = 0.04;      % 风险资产超额回报
sigma_epsilon = 0.157; % 风险资产回报率标准差

% 状态变量网格数量
nW = 50;      % 流动性财富 W
nF = 40;      % 养老金财富 F
nM = 7;       % 持久性收入冲击 M
nEpsilon = 5; % 风险资产收益冲击 Epsilon
nU = 5;       % 暂时性收入冲击 U
nZ = 5;       % 持久性收入冲击 Z

% 并行计算
num_cores = feature('numcores');
fprintf('使用 %d 个CPU核心进行并行计算\n', num_cores);
if isempty(gcp('nocreate'))
    parpool('local', num_cores);
end

% --- 2. 离散化外生冲击 ---
fprintf('2. 开始离散化外生随机冲击...\n');
% 使用 Tauchen 方法离散化正态分布
[z_grid, z_prob] = discretizeAR1_Tauchen(0, 0, sigma_z, nZ, 2, struct('parallel',0)); % z_t ~ N(0, sigma_z^2)
[u_grid, u_prob] = discretizeAR1_Tauchen(0, 0, sigma_u, nU, 2, struct('parallel',0)); % u_t ~ N(0, sigma_u^2)
[epsilon_grid, epsilon_prob] = discretizeAR1_Tauchen(0, 0, sigma_epsilon, nEpsilon, 2, struct('parallel',0)); % epsilon_t ~ N(0, sigma_epsilon^2)
u_shocks = exp(u_grid);
z_shocks = z_grid;
R_risk_shocks = R_f + mu + epsilon_grid; % 风险资产总回报率的实现值

% 构建联合概率矩阵
% prob(z, u, epsilon) = prob(z) * prob(u) * prob(epsilon)
[P_Z, P_U, P_EPSILON] = ndgrid(z_prob, u_prob, epsilon_prob);
trans_prob_work = P_Z .* P_U .* P_EPSILON; % 工作期三冲击联合概率

% --- 3. 构建状态空间网格 ---
fprintf('3. 开始构建状态空间网格...\n');
% 确定性年龄收入 f(t)
age_grid = tb:td;
f_t = b1 * age_grid + b2 * age_grid.^2;

% 持久性收入冲击网格 M_t
m_max = 1.5;
m_min = -1.5;
gM = linspace(m_min, m_max, nM)';

% 流动性财富网格 W_t
w_max = 100;
w_min = 0.01;
gW = logspace(log10(w_min), log10(w_max), nW)';
gW(1) = 0; % 允许破产/零财富

% 养老金账户财富网格 F_t
f_max = 100;
f_min = 0.01;
gF = logspace(log10(f_min), log10(f_max), nF)';
gF(1) = 0;

% 退休期年金网格 P_bar (用于预计算退休期值函数)
p_max = 5;
p_min = 0;
nP = 30;
gP = linspace(p_min, p_max, nP)';

% 条件生存概率
survprob_data = [0.99975, 0.99974, 0.99973, 0.99972, 0.99971, 0.99969, 0.99968, 0.99966, 0.99963, 0.99961, ...
    0.99958, 0.99954, 0.99951, 0.99947, 0.99943, 0.99938, 0.99934, 0.99928, 0.99922, 0.99916, ...
    0.99908, 0.999, 0.99891, 0.99881, 0.9987, 0.99857, 0.99843, 0.99828, 0.99812, 0.99794, ...
    0.99776, 0.99756, 0.99735, 0.99713, 0.9969, 0.99666, 0.9964, 0.99613, 0.99583, 0.99551, ...
    0.99515, 0.99476, 0.99432, 0.99383, 0.9933, 0.9927, 0.99205, 0.99133, 0.99053, 0.98961, ...
    0.98852, 0.98718, 0.98553, 0.98346, 0.98089, 0.97772, 0.97391, 0.96943, 0.96429, 0.95854, ...
    0.95221, 0.94537, 0.93805, 0.93027, 0.92202, 0.91327, 0.90393, 0.89389, 0.88304, 0.87126, ...
    0.85846, 0.84452, 0.82935, 0.81282, 0.79485, 0.77543, 0.75458, 0.7324, 0.70893, 0.68424, ...
    0.68424, 0.68424, 0.68424, 0.68424, 0.68424, 0.68424, 0.68424, 0.68424, 0.68424, 0.68424, ...
    0.68424, 0.68424, 0.68424, 0.68424, 0.68424, 0.68424, 0.68424, 0.68424, 0.68424, 0.68424];
survprob = survprob_data(tb:td-1)'; % 年龄t到t+1的存活概率

% --- 4. 初始化策略和值函数数组 ---
fprintf('4. 开始初始化策略和值函数数组...\n');
% 工作期: 状态 (W, F, M, t)
V_work = zeros(nW, nF, nM, tn);
C_work = zeros(nW, nF, nM, tn);
Q_work = zeros(nW, nF, nM, tn);
A_work = zeros(nW, nF, nM, tn); % alpha_t

% 退休期: 状态 (W, P_bar, t) - 用于预计算
V_retired = zeros(nW, nP, tn);
C_retired = zeros(nW, nP, tn);
A_retired = zeros(nW, nP, tn);

% 优化选项
options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp', ...
    'MaxFunctionEvaluations', 1000, 'MaxIterations', 400, ...
    'OptimalityTolerance', 1e-8, 'StepTolerance', 1e-8, 'ConstraintTolerance', 1e-8);

% --- 5. 模型求解：逆向归纳法 ---
fprintf('5. 开始求解模型（逆向归纳法）...\n');

% 5.1 最后一期 (t=T)
fprintf('   - 求解最后一期 (t=%d)...\n', td);
% 在最后一期，没有遗赠动机，所有流动性财富都被消费
V_retired(:, :, tn) = (gW).^(1-gamma) / (1-gamma); % V(W,P) but P doesn't matter
V_retired(gW==0, :, tn) = -inf; % 0消费的效用为负无穷
C_retired(:, :, tn) = repmat(gW, 1, nP);
A_retired(:, :, tn) = 0;
% 对于工作期最后一年退休的情形，也适用
V_work(:, :, :, tn) = repmat(reshape((gW).^(1-gamma) / (1-gamma), [nW, 1, 1]), [1, nF, nM]);
V_work(gW==0, :, :, tn) = -inf;


% 5.2 退休期 (t = T-1 to K)
fprintf('   - 求解退休期 (t=%d to %d)...\n', td-1, tr);
retire_start_time = tic;
for t_idx = (tn-1):-1:(tr-tb+1)
    age = t_idx + tb - 1;
    fprintf('     正在求解年龄: %d (%d/%d)\n', age, (tn-1)-t_idx+1, (tn-1)-(tr-tb+1)+1);
    
    % 下一期的值函数，用于插值
    next_V = squeeze(V_retired(:, :, t_idx+1));
    
    % 对每个年金水平 P 并行计算
    temp_C = zeros(nW, nP);
    temp_A = zeros(nW, nP);
    temp_V = zeros(nW, nP);
    
    parfor iP = 1:nP
        P_bar = gP(iP);
        
        local_C = zeros(nW, 1);
        local_A = zeros(nW, 1);
        local_V = zeros(nW, 1);
        
        % 插值器
        interp_V_next = griddedInterpolant(gW, next_V(:, iP), 'spline', 'spline');

        for iW = 2:nW % 从iW=2开始，因为iW=1时W=0，无法消费
            W = gW(iW);
            
            % 决策变量 x = [C, alpha]
            % 约束: C > 0, C <= W + Y_bar + P_bar, 0 <= alpha <= 1
            available_resource = W + Y_bar + P_bar;
            x0 = [0.5 * available_resource, 0.5];
            lb = [1e-6, 0];
            ub = [available_resource, 1];
            
            obj_fun = @(x) fun_valuefunc_retired_pps(x, W, Y_bar, P_bar, ...
                interp_V_next, R_f, R_risk_shocks, epsilon_prob, ...
                gamma, beta, survprob(t_idx));
                
            [policy, fval] = fmincon(obj_fun, x0, [], [], [], [], lb, ub, [], options);

            local_C(iW) = policy(1);
            local_A(iW) = policy(2);
            local_V(iW) = -fval;
        end
        local_V(1) = -inf; % W=0时的值
        temp_C(:, iP) = local_C;
        temp_A(:, iP) = local_A;
        temp_V(:, iP) = local_V;
    end
    C_retired(:, :, t_idx) = temp_C;
    A_retired(:, :, t_idx) = temp_A;
    V_retired(:, :, t_idx) = temp_V;
end
retire_time = toc(retire_start_time);
fprintf('   退休期求解完成，耗时 %.2f 分钟\n', retire_time/60);

% 创建一个插值器，用于工作期最后一年到退休第一年的过渡
V_retired_interp = griddedInterpolant({gW, gP}, squeeze(V_retired(:,:,tr-tb+1)), 'spline', 'spline');

% 5.3 工作期 (t = K-1 to 1)
fprintf('   - 求解工作期 (t=%d to %d)...\n', tr-1, tb);
work_start_time = tic;
for t_idx = (tr-tb):-1:1
    age = t_idx + tb - 1;
    fprintf('     正在求解年龄: %d (%d/%d)\n', age, (tr-tb)-t_idx+1, tr-tb);

    temp_C = zeros(nW, nF, nM);
    temp_Q = zeros(nW, nF, nM);
    temp_A = zeros(nW, nF, nM);
    temp_V = zeros(nW, nF, nM);
    
    % 为并行循环准备插值器
    if t_idx == (tr-tb) % 如果是最后工作期
        is_last_work_period = true;
        next_V_interp = V_retired_interp; % 使用退休期插值器
    else
        is_last_work_period = false;
        % 创建下一期工作期值函数的插值器
        next_V = V_work(:, :, :, t_idx+1);
        next_V_interp = griddedInterpolant({gW, gF, gM}, next_V, 'spline', 'spline');
    end

    parfor iM = 1:nM
        M = gM(iM);
        
        local_C = zeros(nW, nF);
        local_Q = zeros(nW, nF);
        local_A = zeros(nW, nF);
        local_V = zeros(nW, nF);
        
        for iF = 1:nF
            F = gF(iF);
            for iW = 2:nW % W=0时无法决策
                W = gW(iW);

                % 决策变量 x = [C, Q, alpha]
                % 约束: C>0, 0<=Q<=Q_max*Y, 0<=alpha<=1, C+S = W+(1-tau_y)*(Y-Q)
                % 我们让优化器选择C,Q,alpha, 储蓄S由预算约束决定
                
                % 由于Y是随机的，约束依赖于冲击的实现值，这使得fmincon难以处理。
                % 解决方法：我们将决策变量定义为相对于某种确定性资源的比例。
                % 但为了严格遵循理论，我们直接优化C, Q, alpha。
                % 我们需要在目标函数内部处理约束。
                
                x0 = [0.3*W, 0.05, 0.5]; % 初始猜测 [C, Q, alpha]
                lb = [1e-6, 0, 0];
                ub = [inf, inf, 1]; % C, Q的上限在函数内部处理
                
                 obj_fun = @(x) fun_valuefunc_work_pps(x, W, F, M, ...
                    next_V_interp, is_last_work_period, ...
                    f_t(t_idx), u_shocks, z_shocks, R_risk_shocks, trans_prob_work, ...
                    R_f, R_p, tau_y, tau_q, Q_max, Y_bar, ...
                    gamma, beta, survprob(t_idx), gW, gF, gM, age, tr, td);

                [policy, fval] = fmincon(obj_fun, x0, [], [], [], [], lb, ub, [], options);

                local_C(iW, iF) = policy(1);
                local_Q(iW, iF) = policy(2);
                local_A(iW, iF) = policy(3);
                local_V(iW, iF) = -fval;
            end
        end
        local_V(1, :) = -inf;
        
        temp_C(:,:,iM) = local_C;
        temp_Q(:,:,iM) = local_Q;
        temp_A(:,:,iM) = local_A;
        temp_V(:,:,iM) = local_V;
    end
    C_work(:, :, :, t_idx) = temp_C;
    Q_work(:, :, :, t_idx) = temp_Q;
    A_work(:, :, :, t_idx) = temp_A;
    V_work(:, :, :, t_idx) = temp_V;
end
work_time = toc(work_start_time);
fprintf('   工作期求解完成，耗时 %.2f 分钟\n', work_time/60);

total_time = retire_time + work_time;
fprintf('模型求解完成，总耗时 %.2f 分钟\n', total_time/60);

% --- 6. 保存结果 ---
fprintf('6. 保存求解结果...\n');
model_results = struct();
model_results.C_work = C_work;
model_results.Q_work = Q_work;
model_results.A_work = A_work;
model_results.V_work = V_work;
model_results.C_retired = C_retired;
model_results.A_retired = A_retired;
model_results.V_retired = V_retired;
model_results.gW = gW;
model_results.gF = gF;
model_results.gM = gM;
model_results.gP = gP;
model_results.params = v2struct(tb, tr, td, tn, beta, gamma, b1, b2, ...
    sigma_u, sigma_z, Y_bar, tau_y, tau_q, Q_max, R_f, R_p, mu, ...
    sigma_epsilon, f_t, survprob);

save(fullfile(save_path, 'model_results_pps.mat'), 'model_results', '-v7.3');
fprintf('模型求解结果已保存到 %s\n', fullfile(save_path, 'model_results_pps.mat'));

%% 数值模拟
fprintf('开始进行数值模拟...\n');
simulate_model_pps(model_results, save_path);

else
%% 直接加载结果并模拟
fprintf('加载已有模型结果并进行模拟...\n');
load(fullfile(save_path, 'model_results_pps.mat'), 'model_results');
simulate_model_pps(model_results, save_path);
end
diary off;

