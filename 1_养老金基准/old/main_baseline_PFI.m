clear
% 使用策略函数迭代方法求解模型的最优消费和投资决策
% 使用parfor在fund和cash层面并行计算

% 创建保存每期数据的目录
if ~exist('result_baseline_matlab_PFI', 'dir')
    mkdir('result_baseline_matlab_PFI');
    mkdir('result_baseline_matlab_PFI/slices');
    mkdir('result_baseline_matlab_PFI/slices/retire');
    mkdir('result_baseline_matlab_PFI/slices/work');
end

% 设置随机种子以获得可重复的结果
rng(42);

% 变量定义
tb = 18;    % 初始开始的年纪
tr = 61;    % 退休年龄
td = 100;   % 死亡年龄
tn = td - tb + 1;  % 总期数

% 状态变量grid数量
ncash = 21;  % 手中现金
nfund = 21;  % 养老基金余额
n = 3;       % 外生随机冲击的grid数量

% 外生参数
% 基础收入f的系数
aa = (-2.170042 + 2.700381);
b1 = 0.16818;
b2 = -0.0323371 / 10;
b3 = 0.0019704 / 100;

% 养老金相关
ret_fac = 0.6827;  % 退休后固定支付的工资（基本养老金）
pension_pct = 0;   % 工资扣除缴纳基本养老保险的比例（个人账户部分）
rp = 1.03;    % 个人养老基金的无风险收益率

% 随机冲击参数
smay = sqrt(0.169993);  % 白噪声shock的标准差
smav = sqrt(0.112572);  % 持续shock的标准差
corr_z_epsilon = 0.0;      % 工资收入白噪声与风险收益随机项的相关性
corr_u_epsilon = 0.0;      % 工资收入AR(1)随机项与风险收益随机项的相关性

% 效用函数参数
gamma = 3.84;   % Epstein-Zin的相对风险规避系数
beta = 0.95;    % Epstein-Zin的discount factor
psi = 0.15;     % Epstein-Zin的跨期替代弹性

% 资产收益率参数
rf = 1.02;    % 无风险总收入
mu = 0.04;    % 超额收益
sigr = 0.27;  % 风险资产收益率的标准差

% PFI迭代参数
max_iterations = 400;  % 最大迭代次数
tolerance = 1e-6;      % 收敛容差
alpha = 0.8;  % 迭代步长

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
gfund = zeros(nfund, 1);
lgfund = zeros(nfund, 1);
secd = zeros(ncash, 1);
C = zeros(ncash, nfund, tn);
V = zeros(ncash, nfund, tn);
A = ones(ncash, nfund, tn);
Q = zeros(ncash, nfund, tn);  % 个人养老金购买决策

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
    gcash(i1) = exp(l_mincash + (i1-1) * stepcash);
end

% fund的grid
maxfund = 100;
minfund = 0.25;
l_maxfund = log(maxfund);
l_minfund = log(minfund);
stepfund = (l_maxfund - l_minfund) / (nfund - 1);

for i1 = 1:nfund
    gfund(i1) = exp(l_minfund + (i1-1) * stepfund);
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
    gy(i1 - tb + 1) = f_y(i1 - tb + 2) / f_y(i1 - tb + 1);  % 对基础工资收入进行normalization
    for i2 = 1:n
        gyp(:, i2, i1 - tb + 1) = gy(i1 - tb + 1) * ones(n, 1) .* exp(yp(:, i2));  % normalized工资收入+持久收入shock
    end
end

% 退休年龄
for i1 = (tr - tb + 1):(tn - 1)
    for i2 = 1:n
        gyp(:, i2, i1) = exp(0.0 * ones(n, 1));
    end
end

% 终止期
% 策略函数
for i1 = 1:ncash
    for i2 = 1:nfund
        C(i1, i2, tn) = 1; % 最后一期的现金全部用来消费
        A(i1, i2, tn) = 0.0;      % 最后一期的投资全部为零
        Q(i1, i2, tn) = 0.0;      % 最后一期不购买养老金
        V(i1, i2, tn) = (C(i1, i2, tn)*gcash(i1)) * (1 - beta)^psi_2;
    end
end

% 计算个人养老金给付
% 假设退休后按照剩余年份平均分配养老基金余额
pension_period = td - tr + 1;
pension_rate = 1 / pension_period;

% 确定可用CPU核心数
num_cores = 14; % feature('numcores');
fprintf('使用 %d 个CPU核心进行并行计算\n', num_cores);

% 设置本地优化选项
options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp', ...
    'MaxFunctionEvaluations', 1000, 'MaxIterations', 400, ...
    'OptimalityTolerance', 1e-8, 'StepTolerance', 1e-8, ...
    'ConstraintTolerance', 1e-8, 'FiniteDifferenceType', 'central');

% 添加全局优化选项
% pso_options = optimoptions('particleswarm', 'Display', 'off', ...
%     'SwarmSize', 50, 'MaxIterations', 100, 'FunctionTolerance', 1e-6, ...
%     'UseParallel',false);

% 策略函数迭代
% 退休期
fprintf('开始求解退休期值函数（使用策略函数迭代）...\n');
retire_start_time = tic; % 开始计时退休期求解

for i1 = 1:(td - tr + 1)
    t = tn - i1;
    period_start_time = tic; % 开始计时当前时期
    fprintf('退休期求解进度: %d/%d\n', i1, td - tr + 1);

    % 初始化策略
    for i2 = 1:nfund
        for i3 = 1:ncash
            C(i3, i2, t) = 0.5; % C(i3, i2, t+1);
            A(i3, i2, t) = 1; % A(i3, i2, t+1);
            Q(i3, i2, t) = 0.0;  % 退休期不需要购买个人养老金
        end
    end

    % 初始化迭代计数器
    iter_count = 0;
    
    % 策略迭代循环
    for iter = 1:max_iterations
        iter_count = iter;  % 更新迭代计数

 

        
        % 策略改进 - 基于当前值函数更新策略，使用并行计算
        C_new = zeros(ncash, nfund);
        A_new = zeros(ncash, nfund);

        % 创建所有(cash, fund)组合的任务列表
        tasks = zeros(ncash * nfund, 2);
        task_idx = 1;
        for i2 = 1:nfund
            for i3 = 1:ncash
                tasks(task_idx, :) = [i3, i2];
                task_idx = task_idx + 1;
            end
        end

        % 设置并行池
        if isempty(gcp('nocreate'))
            parpool('local', num_cores);
        end

        % 在每个(cash, fund)组合上并行
        results = zeros(size(tasks, 1), 3); % 用于存储结果的数组 [C, A, 优化状态]
        
        % 显示进度条初始化
        % fprintf('  PFI迭代 %d/%d: 0%%', iter, max_iterations);
        % progress_step = max(1, round(size(tasks, 1) / 50)); % 每处理50个任务更新一次进度条
        
        parfor task_idx = 1:size(tasks, 1)
            i3 = tasks(task_idx, 1);
            i2 = tasks(task_idx, 2);
            
            x0 = [C(i3, i2, t), A(i3, i2, t)];
            lb = [0, 0];
            ub = [1, 1];
            
            % 确保初始值在边界范围内
            x0 = max(min(x0, ub), lb);
            
            % 计算养老金给付
            pension_pay = pension_rate * gfund(i2);
            
            % 使用fmincon的SQP算法代替模式搜索
            policy = fmincon(@(x) baseline_fun_valuefunc_retired(x, gcash(i3), gfund(i2), ...
                V(:, i2, t+1), gret, rf, ret_fac, pension_pay, ...
                gamma, beta, psi_1, psi_2, theta, gcash, gfund, survprob(t), weig), ...
                x0, [], [], [], [], lb, ub, [], options);
            
            % 将结果存储到临时变量
            results(task_idx, :) = [policy(1), policy(2), 1]; % 1表示优化成功
        end

        % 更新进度条
        % fprintf('\r  PFI迭代 %d/%d: 100%%\n', iter, max_iterations);

        % 将并行结果整合回C_new和A_new
        for task_idx = 1:size(tasks, 1)
            i3 = tasks(task_idx, 1);
            i2 = tasks(task_idx, 2);
            C_new(i3, i2) = results(task_idx, 1);
            A_new(i3, i2) = results(task_idx, 2);
        end

        % 计算策略变化量
        policy_change = max(max(abs(C_new - C(:, :, t)))) + max(max(abs(A_new - A(:, :, t))));

        % 自适应步长计算
        adaptive_alpha = min(0.9, max(0.1, 1.0 / (1.0 + policy_change)));

        % 更新策略
        C(:, :, t) = C_new;
        A(:, :, t) = A_new;

         % 计算值函数 - 逐点计算而不是整个矩阵
        V_new = zeros(ncash, nfund);
        for i3 = 1:ncash
            for i2 = 1:nfund
                % 计算养老金给付
                pension_pay = pension_rate * gfund(i2);
                
                % 计算当前状态的值函数
                V_new(i3, i2) = -baseline_fun_valuefunc_retired([C(i3, i2, t), A(i3, i2, t)], ...
                    gcash(i3), gfund(i2), V(:, i2, t+1), gret, rf, ret_fac, pension_pay, ...
                    gamma, beta, psi_1, psi_2, theta, gcash, gfund, survprob(t), weig, t);
            end
        end
        V(:, :, t) = V_new;        

        % 报告当前policy change
        % fprintf('  退休期时期 %d, 迭代 %d: policy_change = %.6f\n', t, iter, policy_change);

        % 检查收敛
        if policy_change < tolerance
            fprintf('  退休期时期 %d 在 %d 次迭代后收敛, %.6f\n', t, iter_count, policy_change);
            break;
        end

        % 如果达到最大迭代次数仍未收敛
        if iter_count == max_iterations
            fprintf('  警告：退休期时期 %d 在 %d 次迭代后仍未收敛\n', t, max_iterations);
        end

        % % 保存当前时期的切片数据
        % slice_dir = sprintf('result_baseline_matlab_PFI/slices/retire/period_%d', t);
        % if ~exist(slice_dir, 'dir')
        %     mkdir(slice_dir);
        % end
        % 
        % writematrix(C(:, :, t), sprintf('%s/C_slice_t%d.csv', slice_dir, t));
        % writematrix(A(:, :, t), sprintf('%s/A_slice_t%d.csv', slice_dir, t));
        % writematrix(Q(:, :, t), sprintf('%s/Q_slice_t%d.csv', slice_dir, t));
        % writematrix(V(:, :, t), sprintf('%s/V_slice_t%d.csv', slice_dir, t));
        % 
        % % 保存为.mat文件
        % save(sprintf('%s/data_t%d.mat', slice_dir, t), 'C', 'A', 'Q', 'V');
    end
    
    % 输出当前时期求解时间
    period_time = toc(period_start_time);
    fprintf('  退休期时期 %d 求解完成，耗时 %.2f 秒\n', t, period_time);
end

% 输出退休期总求解时间
retire_time = toc(retire_start_time);
fprintf('退休期求解完成，总耗时 %.2f 秒 (%.2f 分钟)\n', retire_time, retire_time/60);

% 工作期（退休前）
fprintf('开始求解工作期值函数（使用策略函数迭代）...\n');
work_start_time = tic; % 开始计时工作期求解

for i1 = 1:(tr - tb)
    t = tr - tb - i1 + 1;
    period_start_time = tic; % 开始计时当前时期
    fprintf('工作期求解进度: %d/%d\n', i1, tr - tb);

    % 初始化策略
    for i2 = 1:nfund
        for i3 = 1:ncash
            C(i3, i2, t) = 0.5; % C(i3, i2, t+1);
            A(i3, i2, t) = 1; % A(i3, i2, t+1);
            Q(i3, i2, t) = 0.1;
        end
    end

    % 初始化迭代计数器
    iter_count = 0;
    
    % 策略迭代循环
    for iter = 1:max_iterations
        iter_count = iter;  % 更新迭代计数
        
        % 策略改进 - 基于当前值函数更新策略，使用并行计算
        C_new = zeros(ncash, nfund);
        A_new = zeros(ncash, nfund);
        Q_new = zeros(ncash, nfund);

        % 创建所有(cash, fund)组合的任务列表
        tasks = zeros(ncash * nfund, 5);  % [i3, i2, C, A, Q]
        task_idx = 1;
        for i2 = 1:nfund
            for i3 = 1:ncash
                tasks(task_idx, 1:2) = [i3, i2];
                task_idx = task_idx + 1;
            end
        end

        % 设置并行池
        if isempty(gcp('nocreate'))
            parpool('local', feature('numcores'));
        end

        % 在每个(cash, fund)组合上并行
        results = zeros(size(tasks, 1), 4); % 用于存储结果的数组 [C, A, Q, 优化状态]
        
        % 显示进度条初始化
        % fprintf('  PFI迭代 %d/%d: 0%%', iter, max_iterations);
        % progress_step = max(1, round(size(tasks, 1) / 50)); % 每处理50个任务更新一次进度条
        
        parfor task_idx = 1:size(tasks, 1)
            i3 = tasks(task_idx, 1);
            i2 = tasks(task_idx, 2);
            
            x0 = [C(i3, i2, t), A(i3, i2, t), Q(i3, i2, t)];
            lb = [0, 0, 0];
            ub = [1, 1, 1];
            
            % 确保初始值在边界范围内
            x0 = max(min(x0, ub), lb);
            
            % 使用fmincon的SQP算法代替模式搜索
            policy = fmincon(@(x) baseline_fun_valuefunc_work(x, gcash(i3), gfund(i2), ...
                gyp(:, :, t), V(:, :, t+1), yh, gret, rf, rp, pension_pct, ...
                gamma, beta, psi_1, psi_2, theta, gcash, gfund, survprob(t), n, nweig1), ...
                x0, [], [], [], [], lb, ub, [], options);
            
            % 将结果存储到临时变量
            results(task_idx, :) = [policy(1), policy(2), policy(3), 1]; % 1表示优化成功
        end

        % 更新进度条
        % fprintf('\r  PFI迭代 %d/%d: 100%%\n', iter, max_iterations);

        % 将并行结果整合回C_new, A_new和Q_new
        for task_idx = 1:size(tasks, 1)
            i3 = tasks(task_idx, 1);
            i2 = tasks(task_idx, 2);
            C_new(i3, i2) = results(task_idx, 1);
            A_new(i3, i2) = results(task_idx, 2);
            Q_new(i3, i2) = results(task_idx, 3);
        end

        % 计算策略变化量
        policy_change = max(max(abs(C_new - C(:, :, t)))) + ...
            max(max(abs(A_new - A(:, :, t)))) + ...
            max(max(abs(Q_new - Q(:, :, t))));

        % 自适应步长计算
        adaptive_alpha = min(0.9, max(0.1, 1.0 / (1.0 + policy_change)));

        % 更新策略
        C(:, :, t) = (1-alpha) * C(:, :, t) + alpha * C_new;
        A(:, :, t) = (1-alpha) * A(:, :, t) + alpha * A_new;
        Q(:, :, t) = (1-alpha) * Q(:, :, t) + alpha * Q_new;
        
        % 计算值函数 - 逐点计算而不是整个矩阵
        V_new = zeros(ncash, nfund);
        for i3 = 1:ncash
            for i2 = 1:nfund
                % 计算当前状态的值函数
                V_new(i3, i2) = -baseline_fun_valuefunc_work([C(i3, i2, t), A(i3, i2, t), Q(i3, i2, t)], ...
                    gcash(i3), gfund(i2), gyp(:, :, t), V(:, :, t+1), yh, gret, rf, rp, pension_pct, ...
                    gamma, beta, psi_1, psi_2, theta, gcash, gfund, survprob(t), n, nweig1, t);
            end
        end
        V(:, :, t) = V_new;     

        % 报告当前policy change
        % fprintf('  工作期时期 %d, 迭代 %d: policy_change = %.6f\n', t, iter, policy_change);

        % 检查收敛
        if policy_change < tolerance
            fprintf('  工作期时期 %d 在 %d 次迭代后收敛,  %.6f\n', t, iter_count, policy_change);
            break;
        end

        % 如果达到最大迭代次数仍未收敛
        if iter_count == max_iterations
            fprintf('  警告：工作期时期 %d 在 %d 次迭代后仍未收敛\n', t, max_iterations);
        end

        % % 保存当前时期的切片数据
        % slice_dir = sprintf('result_baseline_matlab_PFI/slices/work/period_%d', t);
        % if ~exist(slice_dir, 'dir')
        %     mkdir(slice_dir);
        % end
        % 
        % writematrix(C(:, :, t), sprintf('%s/C_slice_t%d.csv', slice_dir, t));
        % writematrix(A(:, :, t), sprintf('%s/A_slice_t%d.csv', slice_dir, t));
        % writematrix(Q(:, :, t), sprintf('%s/Q_slice_t%d.csv', slice_dir, t));
        % writematrix(V(:, :, t), sprintf('%s/V_slice_t%d.csv', slice_dir, t));
        % 
        % % 保存为.mat文件
        % save(sprintf('%s/data_t%d.mat', slice_dir, t), 'C', 'A', 'Q', 'V');
    end
    
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
if ~exist('result_baseline_matlab_PFI', 'dir')
    mkdir('result_baseline_matlab_PFI');
end

% 保存完整的三维策略函数和值函数
save('result_baseline_matlab_PFI/model_results.mat', 'C', 'A', 'Q', 'V', 'gcash', 'gfund', ...
    'tb', 'tr', 'td', 'tn', 'n', 'rf', 'rp', 'mu', 'sigr', 'smay', 'smav', ...
    'ret_fac', 'pension_pct', 'pension_rate', 'gy');

% 同时保存一个切片用于Excel查看
writematrix(C(:, 1, :), 'result_baseline_matlab_PFI/C_slice.csv');
writematrix(A(:, 1, :), 'result_baseline_matlab_PFI/A_slice.csv');
writematrix(Q(:, 1, :), 'result_baseline_matlab_PFI/Q_slice.csv');
writematrix(V(:, 1, :), 'result_baseline_matlab_PFI/V_slice.csv');
writematrix(gcash, 'result_baseline_matlab_PFI/gcash.csv');
writematrix(gfund, 'result_baseline_matlab_PFI/gfund.csv');

fprintf('模型求解完成，结果已保存到result_baseline_matlab_PFI文件夹\n');

% 返回求解结果，供模拟使用
model_results = struct();
model_results.C = C;
model_results.A = A;
model_results.Q = Q;
model_results.V = V;
model_results.gcash = gcash;
model_results.gfund = gfund;
model_results.params = struct();
model_results.params.tb = tb;
model_results.params.tr = tr;
model_results.params.td = td;
model_results.params.tn = tn;
model_results.params.n = n;
model_results.params.rf = rf;
model_results.params.rp = rp;
model_results.params.mu = mu;
model_results.params.sigr = sigr;
model_results.params.smay = smay;
model_results.params.smav = smav;
model_results.params.ret_fac = ret_fac;
model_results.params.pension_pct = pension_pct;
model_results.params.pension_rate = pension_rate;
model_results.params.gy = gy;


function v = baseline_fun_valuefunc_work(x, cash, fund, gyp, nextV, yh, gret, rf, rp, tau, gamma, beta, psi_1, psi_2, theta, gcash, gfund, survprob, n, weig, t)
% 工作期间的值函数计算
%
% 参数:
% x - 决策变量 [消费, 风险资产投资比例, 个人养老金购买比例]
% cash - 手中现金
% fund - 养老基金账户余额
% gyp - 劳动收入增长率
% nextV - 下期值函数
% yh - 暂时收入冲击
% gret - 随机收益率
% rf - 无风险收益率
% rp - 个人养老基金收益率
% tau - 用于缴纳基本养老金的工资比例
% gamma - 相对风险规避系数
% beta - 贴现因子
% psi_1, psi_2, theta - Epstein-Zin效用函数参数
% gcash - 现金网格
% gfund - 养老基金网格
% survprob - 存活概率
% n - 随机状态数量
% weig - 随机状态权重
% t - 当前时期（用于缓存，可选）
%
% 返回:
% v - 函数值

% Epstein-Zin效用函数
% if gamma == 1
%     u = log(cash * x(1));  % 对数效用
% else
u = (cash * x(1) + 1e-07)^psi_1;  % CRRA效用
% end

% 计算个人养老金购买
pension_contrib = x(3) * cash * (1 - x(1));  % 个人养老金购买

% 获取维度
n1 = min(n, size(gyp, 1));  % z_t维度
n2 = min(n, size(gyp, 2));  % eta_t维度
n3 = min(n, size(gyp, 3));  % epsilon_t维度

% 使用meshgrid创建网格
[I6, I7, I8] = meshgrid(1:n1, 1:n2, 1:n3);
I6 = permute(I6, [2, 1, 3]); % 调整维度顺序为(n1,n2,n3)
I7 = permute(I7, [2, 1, 3]);
I8 = permute(I8, [2, 1, 3]);

% 提取相应的gyp值
gyp_values = gyp(sub2ind(size(gyp), I6, I7, I8));

% 计算储蓄值 - 矩阵运算
sav_matrix = cash * (1 - x(1)) * (1 - x(3)) ./ gyp_values;

% 计算投资组合收益率
portfolio_return = rf * (1 - x(2)) + gret(I8) * x(2);

% 提取yh值
yh_values = yh(sub2ind(size(yh), I7, I8));

% 计算下期现金值
cash_1_matrix = portfolio_return .* sav_matrix + (1-tau) * yh_values;

% 限制cash_1的范围
cash_1_matrix = max(min(cash_1_matrix, gcash(end)), gcash(1));

% 计算下期养老基金账户余额
fund_1_matrix = (fund + pension_contrib) * rp ./ gyp_values + tau * yh_values;

% 限制fund_1的范围
fund_1_matrix = max(min(fund_1_matrix, gfund(end)), gfund(1));

% 使用meshgrid为插值准备网格
[X, Y] = meshgrid(gcash, gfund);

% 准备值函数网格
V_grid = permute(nextV, [2, 1]);

% 使用interp2进行插值
% 将3D矩阵重塑为2D以便批量插值
cash_1_flat = reshape(cash_1_matrix, [], 1);
fund_1_flat = reshape(fund_1_matrix, [], 1);

% 批量插值
int_V_flat = interp2(X, Y, V_grid, cash_1_flat, fund_1_flat, 'spline');

% 重塑回3D矩阵
int_V_matrix = reshape(int_V_flat, size(cash_1_matrix));
int_V_matrix = max(int_V_matrix, 1e-20); % 确保插值结果为正

% 提取权重
weig_values = weig(sub2ind(size(weig), I6, I7, I8));

% 计算期望值
expected_value = weig_values .* survprob .* ((int_V_matrix .* gyp_values).^(1-gamma));
auxVV = sum(expected_value(:));

% if auxVV == 0
%     v = -((1-beta)*u)^(psi_2);  % 没有未来值函数时
% else
% Epstein-Zin递归效用
v = -(((1-beta)*u + beta*((auxVV)^(1/theta)))^psi_2);
% end

end

function v = baseline_fun_valuefunc_retired(x, cash, fund, nextV, gret, rf, ret_fac, pension_pay, gamma, beta, psi_1, psi_2, theta, gcash, gfund, survprob, weig, t)
% 退休期间的值函数计算
%
% 参数:
% x - 决策变量 [消费, 风险资产投资比例]
% cash - 手中现金
% fund - 养老基金账户余额
% nextV - 下期值函数
% gret - 随机收益率
% rf - 无风险收益率
% ret_fac - 固定基本养老金
% pension_pay - 个人养老金给付
% gamma - 相对风险规避系数
% beta - 贴现因子
% psi_1, psi_2, theta - Epstein-Zin效用函数参数
% gcash - 现金网格
% gfund - 养老基金网格
% survprob - 存活概率
% weig - 随机状态权重
% t - 当前时期（用于缓存，可选）
%
% 返回:
% v - 函数值

auxVV = 0;
sav = cash * (1 - x(1));

% Epstein-Zin效用函数
u = (cash * x(1) +0.00001).^psi_1;  % CRRA效用


% 下期的cash-on-hand
cash_1 = (rf * (1 - x(2)) + gret .* x(2)) .* sav + pension_pay;

% 使用MATLAB的一维插值
int_V = interp1(gcash, nextV, cash_1,  'spline');
int_V = max(int_V, 1e-20); % 确保插值结果为正

% 计算期望值
auxVV = auxVV + weig' * (survprob * int_V.^(1-gamma));

% Epstein-Zin递归效用
v = -(((1-beta)*u + beta*((auxVV)^(1/theta)))^psi_2);

% 确保返回值是有效的双精度数
end

function [z_grid,P]=discretizeAR1_Tauchen(mew,rho,sigma,znum,Tauchen_q, tauchenoptions)


if exist('tauchenoptions','var')==0
    % Recommended choice for Parallel is 2 (on GPU). It is substantially faster (albeit only for very large grids; for small grids cpu is just as fast)
    tauchenoptions.parallel=1+(gpuDeviceCount>0);
else
    %Check tauchenoptions for missing fields, if there are some fill them with the defaults
    if isfield(tauchenoptions,'parallel')==0
        tauchenoptions.parallel=1+(gpuDeviceCount>0);
    end
end

% Check for a deterministic shifter
if exist('tauchenoptions.dshift','var')==0
    tauchenoptions.dshift=0;
end

if znum==1
    z_grid=mew/(1-rho); %expected value of z
    P=1;
    if tauchenoptions.parallel==2
        z_grid=gpuArray(z_grid);
        P=gpuArray(P);
    end
    return
end

% Note: tauchenoptions.dshift equals zero gives the Tauchen method. 
% For nonzero tauchenoptions.dshift this is actually implementing a non-standard Tauchen method.
if tauchenoptions.parallel==0 || tauchenoptions.parallel==1
    zstar=mew/(1-rho); %expected value of z
    sigmaz=sigma/sqrt(1-rho^2); %stddev of z
    
    z_grid=zstar*ones(znum,1) + linspace(-Tauchen_q*sigmaz,Tauchen_q*sigmaz,znum)';
    omega=z_grid(2)-z_grid(1); %Note that all the points are equidistant by construction.
    
    zi=z_grid*ones(1,znum);
%     zj=ones(znum,1)*z';
    zj=tauchenoptions.dshift*ones(znum,znum)+ones(znum,1)*z_grid';
    
    P_part1=normcdf(zj+omega/2-rho*zi,mew,sigma);
    P_part2=normcdf(zj-omega/2-rho*zi,mew,sigma);
    
    P=P_part1-P_part2;
    P(:,1)=P_part1(:,1);
    P(:,znum)=1-P_part2(:,znum);
    
elseif tauchenoptions.parallel==2 %Parallelize on GPU
    zstar=mew/(1-rho); %expected value of z
    sigmaz=sigma/sqrt(1-rho^2); %stddev of z
    
    z_grid=gpuArray(zstar*ones(znum,1) + linspace(-Tauchen_q*sigmaz,Tauchen_q*sigmaz,znum)');
    omega=z_grid(2)-z_grid(1); %Note that all the points are equidistant by construction.
    
    
    tauchenoptions.dshift=gpuArray(tauchenoptions.dshift*ones(1,znum));
    
    erfinput=arrayfun(@(zi,zj,omega,rho,mew,sigma) ((zj+omega/2-rho*zi)-mew)/sqrt(2*sigma^2), z_grid,tauchenoptions.dshift+z_grid',omega, rho,mew,sigma);
    P_part1=0.5*(1+erf(erfinput));
    
    erfinput=arrayfun(@(zi,zj,omega,rho,mew,sigma) ((zj-omega/2-rho*zi)-mew)/sqrt(2*sigma^2), z_grid,tauchenoptions.dshift+z_grid',omega, rho,mew,sigma);
    P_part2=0.5*(1+erf(erfinput));
    
    P=P_part1-P_part2;
    P(:,1)=P_part1(:,1);
    P(:,znum)=1-P_part2(:,znum);
    
end

end

