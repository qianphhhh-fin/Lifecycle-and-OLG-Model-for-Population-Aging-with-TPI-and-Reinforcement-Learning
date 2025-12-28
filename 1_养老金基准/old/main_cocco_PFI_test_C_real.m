%% main_baseline_PFI_test.m
% 生命周期模型求解与模拟（策略函数迭代方法）
% 参考Python版本main_cocco_PFI.py实现
% 运行模式：求解+模拟

% 清除工作区
clear;
close all;

%% 辅助函数定义

% 生成标准正态分布随机数
function r = f_randn(size)
    r = randn(size);
end

% 寻找值在网格中最接近的位置
function idx = f_ntoil_modified(value, grid, n)
    if value <= grid(1)
        idx = 1;
    elseif value >= grid(n)
        idx = n;
    else
        for i = 1:(n-1)
            if value >= grid(i) && value < grid(i+1)
                idx = i;
                return;
            end
        end
        idx = n;
    end
end

%% 模型求解和模拟

% 创建保存每期数据的目录
if ~exist('result_cocco_matlab_PFI', 'dir')
    mkdir('result_cocco_matlab_PFI');
end

% 设置随机种子以获得可重复的结果
rng(42);

% 变量定义
tb = 18;    % 初始开始的年纪
tr = 61;    % 退休年龄
td = 100;   % 死亡年龄
tn = td - tb + 1;  % 总期数

% 状态变量grid数量 - 匹配Python代码
ncash = 51;  % 手中现金
n = 5;       % 外生随机冲击的grid数量

% 外生参数
% 基础收入f的系数
aa = (-2.170042 + 2.700381);
b1 = 0.16818;
b2 = -0.0323371 / 10;
b3 = 0.0019704 / 100;

% 养老金相关 - 仅保留基本养老金
ret_fac = 0.6827;  % 退休后固定支付的工资（基本养老金）
pension_pct = 0;   % 工资扣除缴纳养老保险的比例，简化版设为0

% 随机冲击参数
smay = sqrt(0.169993);  % 白噪声shock的标准差
smav = sqrt(0.112572);  % 持续shock的标准差
corr_z_epsilon = 0.0;   % 工资收入白噪声与风险收益随机项的相关性
corr_u_epsilon = 0.0;   % 工资收入AR(1)随机项与风险收益随机项的相关性

% 效用函数参数
gamma = 3.84;   % Epstein-Zin的相对风险规避系数
beta = 0.95;    % Epstein-Zin的discount factor
psi = 0.15;     % Epstein-Zin的跨期替代弹性

% 资产收益率参数
rf = 1.02;    % 无风险总收入
mu = 0.04;    % 超额收益
sigr = 0.27;  % 风险资产收益率的标准差

% PFI迭代参数
max_iterations = 100;  % 最大迭代次数，与Python保持一致
tolerance = 1e-4;      % 收敛容差，与Python保持一致

% 初始化数组
survprob = zeros(tn - 1, 1);
beta2 = zeros(tn - 1, 1);
grid = zeros(n, 1);
weig = zeros(n, 1);
gret = zeros(n, 1);
ones_n_1 = ones(n, 1);
grid2 = zeros(n, 1);
yp = zeros(n, n);
yh = zeros(n, n);
nweig1 = zeros(n, n, n);
f_y = zeros(tr - tb + 1, 1);
gy = zeros(tr - tb, 1);
gyp = zeros(n, n, tn - 1);
gcash = zeros(ncash, 1);
lgcash = zeros(ncash, 1);
secd = zeros(ncash, 1);

% 简化为二维数组 - 去掉基金维度
C = zeros(ncash, tn);
V = zeros(ncash, tn);
A = ones(ncash, tn);

% 正态分布的离散化近似
gamma00 = 0;  % AR自相关系数
mew = 0;      % AR的常数项
sigma = 1;    % 白噪声的标准差
tauchenoptions.parallel=0;
[grid, weig_matrix] = discretizeAR1_Tauchen(mew, gamma00, sigma, n, 2, tauchenoptions);
weig = diag(weig_matrix);

% 条件生存概率
survprob_data = [0.99975, 0.99974, 0.99973, 0.99972, 0.99971, 0.99969, 0.99968, 0.99966, 0.99963, 0.99961, ...
    0.99958, 0.99954, 0.99951, 0.99947, 0.99943, 0.99938, 0.99934, 0.99928, 0.99922, 0.99916, ...
    0.99908, 0.999, 0.99891, 0.99881, 0.9987, 0.99857, 0.99843, 0.99828, 0.99812, 0.99794, ...
    0.99776, 0.99756, 0.99735, 0.99713, 0.9969, 0.99666, 0.9964, 0.99613, 0.99583, 0.99551, ...
    0.99515, 0.99476, 0.99432, 0.99383, 0.9933, 0.9927, 0.99205, 0.99133, 0.99053, 0.98961, ...
    0.98852, 0.98718, 0.98553, 0.98346, 0.98089, 0.97772, 0.97391, 0.96943, 0.96429, 0.95854, ...
    0.95221, 0.94537, 0.93805, 0.93027, 0.92202, 0.91327, 0.90393, 0.89389, 0.88304, 0.87126, ...
    0.85846, 0.84452, 0.82935, 0.81282, 0.79485, 0.77543, 0.75458, 0.7324, 0.70893, 0.68424, ...
    0.68424, 0.68424];
for i = 1:min(length(survprob_data), length(survprob))
    survprob(i) = survprob_data(i);
end

% 风险资产收益率
for i1 = 1:n
    gret(i1) = rf + mu + grid(i1) * sigr;
end

% 外生随机性的概率
for i6 = 1:n
    for i7 = 1:n
        for i8 = 1:n
            nweig1(i6, i7, i8) = weig(i6) * weig(i7) * weig(i8);
        end
    end
end

% Epstein-Zin效用函数参数
theta = (1.0 - gamma) / (1.0 - 1.0 / psi);
psi_1 = 1.0 - 1.0 / psi;
psi_2 = 1.0 / psi_1;

% cash-on-hand的grid
maxcash = 100;
mincash = 0.25;
l_maxcash = log(maxcash);
l_mincash = log(mincash);
stepcash = (l_maxcash - l_mincash) / (ncash - 1);

for i1 = 1:ncash
    lgcash(i1) = l_mincash + (i1-1) * stepcash;
end

for i1 = 1:ncash
    gcash(i1) = exp(lgcash(i1));
end

% 劳动收入（归一化）
for i1 = 1:n
    grid2 = grid(i1) .* corr_u_epsilon + grid .* ones(n, 1) .* sqrt(1 - corr_u_epsilon^2);
    yh(:, i1) = exp(grid2 * smay);  % 白噪声随机项u_t的grid
end

for i1 = 1:n
    grid2 = grid(i1) .* corr_z_epsilon + grid .* ones(n, 1) .* sqrt(1 - corr_z_epsilon^2);
    yp(:, i1) = grid2 * smav;  % permanent labor p_t 的随机游走项白噪声z_t的grid
end

% 随年龄变化的基础收入
for i1 = tb:tr
    f_y(i1 - tb + 1) = exp((aa + b1 * i1 + b2 * i1^2 + b3 * i1^3));
end

% 工作年龄
for i1 = tb:(tr-1)
    gy(i1 - tb + 1) = f_y(i1 - tb + 2) / f_y(i1 - tb + 1) - 1.0;  % 对基础工资收入进行normalization
    for i2 = 1:n
        gyp(:, i2, i1 - tb + 1) = exp(gy(i1 - tb + 1) * ones(n, 1) + yp(:, i2));  % normalized工资收入+持久收入shock
    end
end

% 退休年龄
for i1 = (tr - tb + 1):(tn - 1)
    for i2 = 1:n
        gyp(:, i2, i1) = exp(0.0 * ones(n, 1));
    end
end

% 终止期 - 简化为没有养老金维度
for i1 = 1:ncash
    C(i1, tn) = gcash(i1);    % 最后一期的现金全部用来消费
    A(i1, tn) = 0.0;          % 最后一期的投资全部为零
    V(i1, tn) = -(gcash(i1)^(1-gamma)); % 与Python代码一致的效用函数
end

% 设置优化选项 - 使优化参数接近Python的L-BFGS-B
options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'interior-point', ...
    'MaxFunctionEvaluations', 1000, 'MaxIterations', 400, ...
    'OptimalityTolerance', 1e-6, 'StepTolerance', 1e-6, ...
    'ConstraintTolerance', 1e-6);





%% 模型求解 - 使用策略函数迭代方法
fprintf('开始求解模型（策略函数迭代方法）...\n');

% 退休期
fprintf('开始求解退休期值函数...\n');
retire_start_time = tic; % 开始计时退休期求解

for i1 = 1:(td - tr + 1)
    t = tn - i1;
    period_start_time = tic; % 开始计时当前时期
    fprintf('退休期求解进度: %d/%d\n', i1, td - tr + 1);

    % 初始化策略 - 使用合理的猜测
    for i3 = 1:ncash
        % 初始策略猜测
        if i3 == 1
            C(i3, t) = 0.2;
            A(i3, t) = 0.2;
        else
            % 使用上一个状态的策略作为初始猜测
            C(i3, t) = C(i3-1, t);
            A(i3, t) = A(i3-1, t);
        end
    end


        % 策略改进 - 基于当前值函数更新策略
        C_new = zeros(ncash, 1);
        A_new = zeros(ncash, 1);
        V_new = zeros(ncash, 1);

        for i3 = 1:ncash
            x0 = [0.5, 0.8];
            lb = [0, 0];
            ub = [0.999 * gcash(i3), 1];

            if i3 ~= 1
                lb = [C_new(i3-1), 0];  % 使用当前迭代中前一个状态的最优策略
            end
            
            % 确保初始值在边界范围内
            x0(1) = max(min(x0(1), ub(1)), lb(1));
            x0(2) = max(min(x0(2), ub(2)), lb(2));

            % 优化函数，使用fmincon模拟Python的minimize函数
            [policy, fval, exitflag] = fmincon(@(x) evaluate_policy_retired(x, gcash(i3), V(:, t+1), gret, rf, ret_fac, gamma, beta, psi_1, psi_2, theta, gcash, survprob(t), weig), ...
                x0, [], [], [], [], lb, ub, [], options);
            
            C_new(i3) = policy(1);
            A_new(i3) = policy(2);
            V_new(i3) = -fval;
        end


        % 更新策略 
        C(:, t) = C_new;
        A(:, t) = A_new;
        V(:, t) = V_new;
                
    
    % 输出当前时期求解时间
    period_time = toc(period_start_time);
    fprintf('  退休期时期 %d 求解完成，耗时 %.2f 秒\n', t, period_time);
end

% 输出退休期总求解时间
retire_time = toc(retire_start_time);
fprintf('退休期求解完成，总耗时 %.2f 秒 (%.2f 分钟)\n', retire_time, retire_time/60);

% 工作期（退休前）
fprintf('开始求解工作期值函数...\n');
work_start_time = tic; % 开始计时工作期求解

for i1 = 1:(tr - tb)
    t = tr - tb - i1 + 1;
    period_start_time = tic; % 开始计时当前时期
    fprintf('工作期求解进度: %d/%d\n', i1, tr - tb);

    % 初始化策略 - 使用合理的猜测
    for i3 = 1:ncash
        % 初始策略猜测
        if i3 == 1
            C(i3, t) = 0.1;
            A(i3, t) = 0.1;
        else
            % 使用上一个状态的策略作为初始猜测
            C(i3, t) = C(i3-1, t);
            A(i3, t) = A(i3-1, t);
        end
    end
        % 策略改进 - 基于当前值函数更新策略
        C_new = zeros(ncash, 1);
        A_new = zeros(ncash, 1);
        V_new = zeros(ncash, 1);

        for i3 = 1:ncash
            x0 = [C(i3, t), A(i3, t)];
            lb = [0, 0];
            ub = [0.999 * gcash(i3), 1];

            if i3 ~= 1
                lb = [C_new(i3-1), 0];  % 使用当前迭代中前一个状态的最优策略
            end
            
            % 确保初始值在边界范围内
            x0(1) = max(min(x0(1), ub(1)), lb(1));
            x0(2) = max(min(x0(2), ub(2)), lb(2));

            % 优化函数，使用fmincon模拟Python的minimize函数
            [policy, fval, exitflag] = fmincon(@(x) evaluate_policy_work(x, gcash(i3), gyp(:, :, t), V(:, t+1), yh, gret, rf, gamma, beta, psi_1, psi_2, theta, gcash, survprob(t), n, weig), ...
                x0, [], [], [], [], lb, ub, [], options);
            
            C_new(i3) = policy(1);
            A_new(i3) = policy(2);
            V_new(i3) = -fval;
        end


        % 更新策略
        C(:, t) = C_new;
        A(:, t) = A_new;
        V(:, t) = V_new;
    
    % 输出当前时期求解时间
    period_time = toc(period_start_time);
    fprintf('  工作期时期 %d 求解完成，耗时 %.2f 秒\n', t, period_time);
end

% 输出工作期总求解时间
work_time = toc(work_start_time);
fprintf('工作期求解完成，总耗时 %.2f 秒 (%.2f 分钟)\n', work_time, work_time/60);

% 输出总求解时间
total_time = retire_time + work_time;
fprintf('模型求解完成，总耗时 %.2f 秒 (%.2f 分钟)\n', total_time, total_time/60);

% 保存结果
if ~exist('result_cocco_matlab_PFI', 'dir')
    mkdir('result_cocco_matlab_PFI');
end

% 保存完整的二维策略函数和值函数
writematrix(C, 'result_cocco_matlab_PFI/C.csv');
writematrix(A, 'result_cocco_matlab_PFI/A.csv');
writematrix(V, 'result_cocco_matlab_PFI/V.csv');
writematrix(gcash, 'result_cocco_matlab_PFI/gcash.csv');
save('result_cocco_matlab_PFI/model_results.mat', 'C', 'A', 'V', 'gcash', ...
    'tb', 'tr', 'td', 'tn', 'n', 'rf', 'mu', 'sigr', 'smay', 'smav', ...
    'ret_fac', 'gy');

fprintf('模型求解完成，结果已保存到result_cocco_matlab_PFI文件夹\n');

% 返回求解结果，供模拟使用
model_results = struct();
model_results.C = C;
model_results.A = A;
model_results.V = V;
model_results.gcash = gcash;
model_results.params = struct();
model_results.params.tb = tb;
model_results.params.tr = tr;
model_results.params.td = td;
model_results.params.tn = tn;
model_results.params.n = n;
model_results.params.rf = rf;
model_results.params.mu = mu;
model_results.params.sigr = sigr;
model_results.params.smay = smay;
model_results.params.smav = smav;
model_results.params.ret_fac = ret_fac;
model_results.params.gy = gy;

%% 数值模拟
% 调用模拟函数
fprintf('开始进行数值模拟...\n');
simulate_model(model_results);

%% 数值模拟函数
function simulate_model(model_results)
    % 数值模拟函数
    % 基于求解结果进行数值模拟
    %
    % 参数：
    %   model_results: 模型求解结果

    % 从模型结果中获取参数和数据
    C = model_results.C;
    A = model_results.A;
    gcash = model_results.gcash;
    
    % 获取参数
    params = model_results.params;
    tb = params.tb;
    tr = params.tr;
    td = params.td;
    tn = params.tn;
    n = params.n;
    rf = params.rf;
    mu = params.mu;
    sigr = params.sigr;
    smay = params.smay;
    smav = params.smav;
    ret_fac = params.ret_fac;
    gy = params.gy;
    
    fprintf('数值模拟开始...\n');
    
    % 设置随机种子以获得可重复的结果
    rng(42);
    
    % 数值模拟参数
    nsim = 10000;
    ones_nsim_1 = ones(nsim, 1);
    meanY = zeros(tn, 1);
    meanC = zeros(tn, 1);
    meanW = zeros(tn, 1);
    meanA = zeros(tn, 1);
    meanS = zeros(tn, 1);
    meanB = zeros(tn, 1);
    meanWY = zeros(tn, 1);
    meanalpha = zeros(tn, 1);
    meanGPY = zeros(tn, 1);
    cGPY = zeros(tn, 1);
    meanYs = zeros(tn, 1);
    meanCs = zeros(tn, 1);
    meanWs = zeros(tn, 1);
    simPY = zeros(tn, nsim);
    simGPY = zeros(tn, nsim);
    simY = zeros(tn, nsim);
    simC = zeros(tn, nsim);
    simW = zeros(tn, nsim);
    simA = zeros(tn, nsim);
    simS = zeros(tn, nsim);
    simB = zeros(tn, nsim);
    simW_Y = zeros(tn, nsim);
    simR = zeros(tn, nsim);
    eps_y = zeros(1, 1);
    simTY = zeros(1, 1);
    eps_r = zeros(1, 1);
    cash = zeros(tn, nsim);
    simC_pct = zeros(tn, nsim);
    
    % 1、模拟生成labor income
    fprintf('生成劳动收入模拟数据...\n');
    for i1 = 1:(nsim/2) % 另外一半模拟完全对称
        % 工作期第一期
        eps_y(1, 1) = randn(1); % N(0,1)
        simPY(1, i1) = eps_y(1, 1) * smav; % 初始的p
        simPY(1, nsim/2 + i1) = -eps_y(1, 1) * smav;
        simGPY(1, i1) = 1.0;
        simGPY(1, nsim/2 + i1) = 1.0;
        simTY(1, 1) = randn(1);
        simY(1, i1) = exp(simTY(1, 1) * smay);
        simY(1, nsim/2 + i1) = exp(-simTY(1, 1) * smay);
        
        % 工作期第2期~退休
        for i2 = 2:(tr-tb)
            w = i2 + tb - 1;
            eps_y(1, 1) = randn(1);
            simPY(i2, i1) = eps_y(1, 1) * smav + simPY(i2-1, i1);
            simPY(i2, nsim/2 + i1) = -eps_y(1, 1) * smav + simPY(i2-1, nsim/2 + i1);
            simGPY(i2, i1) = exp(gy(i2-1)) * exp(simPY(i2, i1)) / exp(simPY(i2-1, i1));
            simGPY(i2, nsim/2 + i1) = exp(gy(i2-1)) * exp(simPY(i2, nsim/2 + i1)) / exp(simPY(i2-1, nsim/2 + i1));
            simTY(1, 1) = randn(1);
            simY(i2, i1) = exp(simTY(1, 1) * smay);
            simY(i2, nsim/2 + i1) = exp(-simTY(1, 1) * smay);
        end
    end
    
    % 退休期
    for t = (tr-tb+1):tn
        simY(t, :) = ret_fac;
        simGPY(t, :) = 1.0;
    end
    
    % 2、模拟风险投资的收益率
    fprintf('生成风险投资收益率模拟数据...\n');
    for t = 1:tn
        for i1 = 1:(nsim/2)
            eps_r(1, 1) = randn(1);
            simR(t, i1) = mu + rf + sigr * eps_r(1, 1);
            simR(t, nsim/2 + i1) = mu + rf - sigr * eps_r(1, 1);
        end
    end
    
    % 3、从第一期开始迭代，得到各控制变量的值
    fprintf('模拟控制变量...\n');
    simW(:, :) = 0; % 初始财富设置为0
    
    for t = 1:tn
        fprintf('模拟进度: %d/%d\n', t, tn);
        
        % 向量化计算所有模拟的财富收入比和现金持有量
        if any(simY(t, :) == 0)
            warning('第%d期存在零收入，已调整为极小值', t);
            simY(t, simY(t, :) == 0) = 1e-6;
        end
        
        simW_Y(t, :) = simW(t, :) ./ simY(t, :); % 上期财富-本期工资收入比
        cash(t, :) = simW(t, :) + simY(t, :); % cash-on-hand
        
        % 处理cash超出范围的情况
        cash_t = cash(t, :);
        cash_t = max(min(cash_t, gcash(end)), gcash(1));
        
        % 应用插值获取策略 - 使用cubic插值，与Python保持一致
        simC(t, :) = interp1(gcash, C(:, t), cash_t, 'spline', 'extrap');
        simA(t, :) = interp1(gcash, A(:, t), cash_t, 'spline', 'extrap');
        
        % 确保约束条件满足
        simC(t, :) = max(min(simC(t, :), 0.9999 * cash(t, :)), 0);
        simA(t, :) = min(max(simA(t, :), 0), 1);
        
        % 计算各种模拟变量
        simC_pct(t, :) = simC(t, :) ./ cash(t, :);
        sav = cash(t, :) - simC(t, :); % 用于投资的金额
        simS(t, :) = simA(t, :) .* sav; % 风险投资额
        simS(t, :) = min(simS(t, :), sav); % 确保风险投资不超过总投资
        simB(t, :) = sav - simS(t, :); % 无风险投资额
        
        % 更新消费比例 - 与Python保持一致
        simC(t, :) = simC(t, :) ./ cash(t, :);
        
        % 计算下期财富
        if t < tn
            % 确保没有零除错误
            if any(simGPY(t+1, :) == 0)
                warning('第%d期存在零GPY，已调整为1', t+1);
                simGPY(t+1, simGPY(t+1, :) == 0) = 1;
            end
            
            simW(t+1, :) = (simB(t, :) * rf + simS(t, :) .* simR(t, :)) ./ simGPY(t+1, :);
        end
    end
    
    % 多次模拟path下变量平均值
    meanC = mean(simC, 2);
    meanC_pct = mean(simC_pct, 2);
    meanY = mean(simY, 2);
    meanW = mean(simW, 2);
    meanS = mean(simS, 2);
    meanB = mean(simB, 2);
    meanWY = mean(simW_Y, 2);
    meanalpha = mean(simA, 2);
    meanGPY = mean(simGPY, 2);
    
    % 保存模拟结果
    if ~exist('result_cocco_matlab_PFI/simulation', 'dir')
        mkdir('result_cocco_matlab_PFI/simulation');
    end
    
    % 创建表格保存模拟结果
    sim_results_table = table(meanC, meanC_pct, meanY, meanW, meanS, meanB, meanWY, meanalpha, meanGPY);
    writetable(sim_results_table, 'result_cocco_matlab_PFI/simulation/simulation_results.csv');
    
    % 保存原始模拟数据用于进一步分析
    save('result_cocco_matlab_PFI/simulation/raw_simulation.mat', 'simC', 'simA', 'simY', 'simW', 'simS', 'simB', 'simR', 'simGPY', '-v7.3');
    
    % 绘制结果
    figure;
    plot(meanC_pct, 'LineWidth', 2);
    hold on;
    plot(meanalpha, 'LineWidth', 2);
    hold off;
    legend('消费比例 (c)', '风险资产配置比例 (\alpha)', 'Location', 'best');
    xlabel('期数');
    ylabel('比例');
    title('消费比例和资产配置比例 (PFI方法)');
    grid on;
    
    % 保存图像
    saveas(gcf, 'result_cocco_matlab_PFI/simulation/results_plot.png');
    saveas(gcf, 'result_cocco_matlab_PFI/simulation/results_plot.fig');
    
    fprintf('数值模拟完成，结果已保存到result_cocco_matlab_PFI/simulation文件夹\n');
end


%% AR(1)过程离散化的Tauchen方法
function [z_grid, P] = discretizeAR1_Tauchen(mew, rho, sigma, znum, Tauchen_q, tauchenoptions)
% Tauchen方法离散化AR(1)过程为马尔可夫链
%
% 参数:
% mew - 常数项
% rho - 自相关系数
% sigma - 创新的标准差
% znum - 离散化后的状态数（必须是奇数）
% Tauchen_q - 定义网格点的超参数
% tauchenoptions - 额外选项
%
% 返回:
% z_grid - 离散的状态值
% P - 转移矩阵

if znum == 1
    z_grid = mew / (1 - rho);
    P = 1;
    return;
end

% 标准方法
zstar = mew / (1 - rho);  % z的期望值
sigmaz = sigma / sqrt(1 - rho^2);  % z的标准差

z_grid = zstar * ones(znum, 1) + linspace(-Tauchen_q * sigmaz, Tauchen_q * sigmaz, znum)';
omega = z_grid(2) - z_grid(1);  % 网格点等距

P = zeros(znum, znum);

for i = 1:znum
    for j = 1:znum
        if j == 1
            P(i, j) = normcdf((z_grid(j) + omega/2 - rho * z_grid(i) - mew) / sigma);
        elseif j == znum
            P(i, j) = 1 - normcdf((z_grid(j) - omega/2 - rho * z_grid(i) - mew) / sigma);
        else
            P(i, j) = normcdf((z_grid(j) + omega/2 - rho * z_grid(i) - mew) / sigma) - ...
                     normcdf((z_grid(j) - omega/2 - rho * z_grid(i) - mew) / sigma);
        end
    end
end
end


%% 退休期间进行策略评估的函数
function v = evaluate_policy_retired(policy, cash, nextV, gret, rf, ret_fac, gamma, beta, psi_1, psi_2, theta, gcash, survprob, weig)
% 对退休期间的给定策略进行评估，计算其值函数
%
% 参数:
% policy - 策略 [消费, 风险资产投资比例]
% cash - 手中现金
% nextV - 下期值函数
% gret - 随机收益率
% rf - 无风险收益率
% ret_fac - 固定养老金
% gamma - 相对风险规避系数
% beta - 贴现因子
% psi_1, psi_2, theta - Epstein-Zin效用函数参数
% gcash - 现金网格
% survprob - 存活概率
% weig - 随机状态权重
%
% 返回:
% v - 值函数值

auxVV = 0;
sav = cash - policy(1);

% 效用函数
u = -(policy(1))^(1-gamma);

% 下期的cash-on-hand
cash_1 = (rf * (1 - policy(2)) + gret .* policy(2)) .* sav + ret_fac;

% 限制cash_1的范围
cash_1 = max(min(cash_1, gcash(end)), gcash(1));

% 使用cubic插值，与Python保持一致
int_V = interp1(gcash, nextV, cash_1, 'spline', 'extrap');

% 计算期望值
auxVV = auxVV + weig' * (survprob * int_V);

% 和Python一致的返回值
if auxVV == 0
    v = -(u);
else
    v = -(u + beta * auxVV);
end
end

%% 工作期间的值函数计算
function v = evaluate_policy_work(x, cash, gyp, nextV, yh, gret, rf, gamma, beta, psi_1, psi_2, theta, gcash, survprob, n, weig)
% 工作期间的值函数计算
%
% 参数:
% x - 决策变量 [消费, 风险资产投资比例]
% cash - 手中现金
% gyp - 劳动收入增长率
% nextV - 下期值函数
% yh - 暂时收入冲击
% gret - 随机收益率
% rf - 无风险收益率
% gamma - 相对风险规避系数
% beta - 贴现因子
% psi_1, psi_2, theta - Epstein-Zin效用函数参数
% gcash - 现金网格
% survprob - 存活概率
% n - 随机状态数量
% weig - 随机状态权重
%
% 返回:
% v - 函数值

auxVV = 0;

% 与Python一致的效用函数
u = -(x(1))^(1-gamma);

% 创建索引网格
[I6, I7, I8] = ndgrid(1:n, 1:n, 1:n);

% 计算所有可能的储蓄值
sav_values = zeros(n, n, n);
cash_1_values = zeros(n, n, n);

% 循环计算各状态下的值
for i6 = 1:n
    for i7 = 1:n
        for i8 = 1:n
            % 计算储蓄值
            sav_values(i6, i7, i8) = (cash - x(1)) / gyp(i6, i8);
            
            % 计算投资组合收益率
            portfolio_return = rf * (1 - x(2)) + gret(i8) * x(2);
            
            % 计算下期现金值
            cash_1_values(i6, i7, i8) = portfolio_return * sav_values(i6, i7, i8) + yh(i7, i8);
        end
    end
end

% 限制cash_1的范围
cash_1_values = max(min(cash_1_values, gcash(end)), gcash(1));

% 对每个现金值应用插值
int_V_values = zeros(n, n, n);
for i6 = 1:n
    for i7 = 1:n
        for i8 = 1:n
            % 使用cubic插值，与Python保持一致
            int_V_values(i6, i7, i8) = interp1(gcash, nextV, cash_1_values(i6, i7, i8), 'spline', 'extrap');
        end
    end
end

% 计算加权和，类似于numpy.sum(weig * survprob * int_V_values * gyp)
weighted_sum = 0;
for i6 = 1:n
    for i7 = 1:n
        for i8 = 1:n
            weighted_sum = weighted_sum + weig(i6) * weig(i7) * weig(i8) * survprob * int_V_values(i6, i7, i8) * gyp(i6, i8);
        end
    end
end
auxVV = weighted_sum;

% 和Python一致的返回值
if auxVV == 0
    v = -(u);
else
    v = -(u + beta * auxVV);
end
end 