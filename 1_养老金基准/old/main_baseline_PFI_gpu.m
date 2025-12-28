clear
% 使用策略函数迭代方法求解模型的最优消费和投资决策
% 使用GPU并行计算和网格搜索替代fmincon

% 检查GPU可用性
if gpuDeviceCount > 0
    gpu_device = gpuDevice();
    fprintf('使用GPU: %s\n', gpu_device.Name);
    fprintf('GPU内存: %.2f GB\n', gpu_device.AvailableMemory/1024/1024/1024);
    use_gpu = true;
    
    % 重置GPU以清除内存
    reset(gpu_device);
else
    fprintf('未检测到GPU，使用CPU计算\n');
    use_gpu = false;
end


% 创建保存每期数据的目录
if ~exist('result_baseline_matlab_PFI_GPU', 'dir')
    mkdir('result_baseline_matlab_PFI_GPU');
    mkdir('result_baseline_matlab_PFI_GPU/slices');
    mkdir('result_baseline_matlab_PFI_GPU/slices/retire');
    mkdir('result_baseline_matlab_PFI_GPU/slices/work');
end

% 设置随机种子以获得可重复的结果
rng(42);

% 变量定义
tb = 18;    % 初始开始的年纪
tr = 61;    % 退休年龄
td = 100;   % 死亡年龄
tn = td - tb + 1;  % 总期数

% 状态变量grid数量
ncash = 5;  % 手中现金
nfund = 5;  % 养老基金余额
n = 3;       % 外生随机冲击的grid数量

% 网格搜索参数
grid_size_retired = [8, 8];  % 退休期网格搜索的每维网格点数 [消费, 投资]
grid_size_work = [6, 6, 6];  % 工作期网格搜索的每维网格点数 [消费, 投资, 养老金]

% 外生参数
% 基础收入f的系数
aa = (-2.170042 + 2.700381);
b1 = 0.16818;
b2 = -0.0323371 / 10;
b3 = 0.0019704 / 100;

% 养老金相关
ret_fac = 0.6827;  % 退休后固定支付的工资（基本养老金）
pension_pct = 0.08;   % 工资扣除缴纳基本养老保险的比例（个人账户部分）
rp = 1.04;    % 个人养老基金的无风险收益率

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
tolerance = 1e-4;      % 收敛容差

% 初始化数组
if use_gpu
    % 使用GPU数组
    survprob = gpuArray(zeros(tn - 1, 1));
    beta2 = gpuArray(zeros(tn - 1, 1));
    grid = gpuArray(zeros(n, 1));
    weig = gpuArray(zeros(n, 1));
    gret = gpuArray(zeros(n, 1));
    ones_n_1 = gpuArray(ones(n, 1));
    grid2 = gpuArray(zeros(n, 1));
    yp = gpuArray(zeros(n, n));
    yh = gpuArray(zeros(n, n));
    nweig1 = gpuArray(zeros(n, n, n));
    f_y = gpuArray(zeros(tr - tb + 1, 1));
    gy = gpuArray(zeros(tr - tb, 1));
    gyp = gpuArray(zeros(n, n, tn - 1));
    gcash = gpuArray(zeros(ncash, 1));
    lgcash = gpuArray(zeros(ncash, 1));
    gfund = gpuArray(zeros(nfund, 1));
    lgfund = gpuArray(zeros(nfund, 1));
    secd = gpuArray(zeros(ncash, 1));
    C = gpuArray(zeros(ncash, nfund, tn));
    V = gpuArray(zeros(ncash, nfund, tn));
    A = gpuArray(ones(ncash, nfund, tn));
    Q = gpuArray(zeros(ncash, nfund, tn));  % 个人养老金购买决策
else
    % 使用CPU数组
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
end

% 正态分布的离散化近似
gamma00 = 0;  % AR自相关系数
mew = 0;      % AR的常数项
sigma = 1;    % 白噪声的标准差
tauchenoptions.parallel=0;
[grid_cpu, weig_matrix_cpu] = discretizeAR1_Tauchen(mew, gamma00, sigma, n, 2, tauchenoptions);
if use_gpu
    grid = gpuArray(grid_cpu);
    weig_matrix = gpuArray(weig_matrix_cpu);
else
    grid = grid_cpu;
    weig_matrix = weig_matrix_cpu;
end
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
num_cores = 14; %feature('numcores');
fprintf('使用 %d 个PU核心进行并行计算\n', num_cores);

% 设置优化选项
options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp', ...
    'MaxFunctionEvaluations', 200, 'OptimalityTolerance', 1e-7);

% 策略函数迭代
% 退休期
fprintf('开始求解退休期值函数（使用GPU网格搜索）...\n');
retire_start_time = tic; % 开始计时退休期求解

for i1 = 1:(td - tr + 1)
    t = tn - i1;
    period_start_time = tic; % 开始计时当前时期
    fprintf('退休期求解进度: %d/%d\n', i1, td - tr + 1);

    % 初始化策略
    for i2 = 1:nfund
        for i3 = 1:ncash
            C(i3, i2, t) = C(i3, i2, t+1);
            A(i3, i2, t) = A(i3, i2, t+1);
            Q(i3, i2, t) = 0.0;  % 退休期不需要购买个人养老金
        end
    end

    % 初始化迭代计数器
    iter_count = 0;
    
    % 策略迭代循环
    for iter = 1:max_iterations
        iter_count = iter;  % 更新迭代计数
        
        % 策略改进 - 基于当前值函数更新策略，使用GPU网格搜索
        C_new = zeros(ncash, nfund, 'like', C);
        A_new = zeros(ncash, nfund, 'like', A);

        % 创建所有(cash, fund)组合的任务列表
        tasks = zeros(ncash * nfund, 2);
        task_idx = 1;
        for i2 = 1:nfund
            for i3 = 1:ncash
                tasks(task_idx, :) = [i3, i2];
                task_idx = task_idx + 1;
            end
        end

        % 在每个(cash, fund)组合上并行
        results = zeros(size(tasks, 1), 3); % 用于存储结果的数组 [C, A, 优化状态]
        
        parfor task_idx = 1:size(tasks, 1)
            i3 = tasks(task_idx, 1);
            i2 = tasks(task_idx, 2);
            
            % 初始猜测
            x0 = [C(i3, i2, t), A(i3, i2, t)];
            lb = [0, 0];
            ub = [1, 1];
            
            % 确保初始值在边界范围内
            x0 = max(min(x0, ub), lb);
            
            % 计算养老金给付
            pension_pay = pension_rate * gfund(i2);
            
            % 使用GPU网格搜索优化
            [policy, ~] = gpu_grid_search(@(x) gpu_baseline_fun_valuefunc_retired(x, gcash(i3), gfund(i2), ...
                V(:, i2, t+1), gret, rf, ret_fac, pension_pay, ...
                gamma, beta, psi_1, psi_2, theta, gcash, gfund, survprob(t), weig, t), ...
                lb, ub, grid_size_retired);
                
            % 将结果存储到临时变量
            results(task_idx, :) = [policy(1), policy(2), 1]; % 1表示优化成功
        end

        % 将并行结果整合回C_new和A_new
        for task_idx = 1:size(tasks, 1)
            i3 = tasks(task_idx, 1);
            i2 = tasks(task_idx, 2);
            C_new(i3, i2) = results(task_idx, 1);
            A_new(i3, i2) = results(task_idx, 2);
        end

        % 计算策略变化量
        policy_change = max(max(abs(C_new - C(:, :, t)))) + max(max(abs(A_new - A(:, :, t))));

        % 自适应步长
        alpha = 0.1;

        % 更新策略
        C(:, :, t) = (1-alpha) * C(:, :, t) + alpha * C_new;
        A(:, :, t) = (1-alpha) * A(:, :, t) + alpha * A_new;
        
        % 计算值函数 - 逐点计算而不是整个矩阵
        V_new = zeros(ncash, nfund, 'like', V);
        for i3 = 1:ncash
            for i2 = 1:nfund
                % 计算养老金给付
                pension_pay = pension_rate * gfund(i2);
                
                % 计算当前状态的值函数
                V_new(i3, i2) = -gpu_baseline_fun_valuefunc_retired([C(i3, i2, t), A(i3, i2, t)], ...
                    gcash(i3), gfund(i2), V(:, i2, t+1), gret, rf, ret_fac, pension_pay, ...
                    gamma, beta, psi_1, psi_2, theta, gcash, gfund, survprob(t), weig, t);
            end
        end
        V(:, :, t) = V_new;

        % 检查收敛
        if policy_change < tolerance
            fprintf('  退休期时期 %d 在 %d 次迭代后收敛\n', t, iter_count);
            break;
        end

        % 如果达到最大迭代次数仍未收敛
        if iter_count == max_iterations
            fprintf('  警告：退休期时期 %d 在 %d 次迭代后仍未收敛\n', t, max_iterations);
        end

        % 定期清理GPU内存
        if use_gpu && mod(iter, 10) == 0
            wait(gpuDevice);
        end
    end
    
    % 输出当前时期求解时间
    period_time = toc(period_start_time);
    if use_gpu
        fprintf('  退休期时期 %d 求解完成，耗时 %.2f 秒，GPU内存: %.2f GB\n', t, period_time, get_gpu_memory());
    else
        fprintf('  退休期时期 %d 求解完成，耗时 %.2f 秒\n', t, period_time);
    end
end

% 输出退休期总求解时间
retire_time = toc(retire_start_time);
fprintf('退休期求解完成，总耗时 %.2f 秒 (%.2f 分钟)\n', retire_time, retire_time/60);

% 工作期（退休前）
fprintf('开始求解工作期值函数（使用GPU网格搜索）...\n');
work_start_time = tic; % 开始计时工作期求解

for i1 = 1:(tr - tb)
    t = tr - tb - i1 + 1;
    period_start_time = tic; % 开始计时当前时期
    fprintf('工作期求解进度: %d/%d\n', i1, tr - tb);

    % 初始化策略
    for i2 = 1:nfund
        for i3 = 1:ncash
            C(i3, i2, t) = C(i3, i2, t+1);
            A(i3, i2, t) = A(i3, i2, t+1);
            Q(i3, i2, t) = 0.1;
        end
    end

    % 初始化迭代计数器
    iter_count = 0;
    
    % 策略迭代循环
    for iter = 1:max_iterations
        iter_count = iter;  % 更新迭代计数
        
        % 策略改进 - 基于当前值函数更新策略，使用GPU网格搜索
        C_new = zeros(ncash, nfund, 'like', C);
        A_new = zeros(ncash, nfund, 'like', A);
        Q_new = zeros(ncash, nfund, 'like', Q);

        % 创建所有(cash, fund)组合的任务列表
        tasks = zeros(ncash * nfund, 2);
        task_idx = 1;
        for i2 = 1:nfund
            for i3 = 1:ncash
                tasks(task_idx, :) = [i3, i2];
                task_idx = task_idx + 1;
            end
        end

        % 在每个(cash, fund)组合上并行
        results = zeros(size(tasks, 1), 4); % 用于存储结果的数组 [C, A, Q, 优化状态]
        
        for task_idx = 1:size(tasks, 1)
            i3 = tasks(task_idx, 1);
            i2 = tasks(task_idx, 2);
            
            % 初始猜测
            x0 = [C(i3, i2, t), A(i3, i2, t), Q(i3, i2, t)];
            lb = [0, 0, 0];
            ub = [1, 1, 1];
            
            % 确保初始值在边界范围内
            x0 = max(min(x0, ub), lb);
            
            % 使用GPU网格搜索优化
            [policy, ~] = gpu_grid_search(@(x) gpu_baseline_fun_valuefunc_work(x, gcash(i3), gfund(i2), ...
                gyp(:, :, t), V(:, :, t+1), yh, gret, rf, rp, pension_pct, ...
                gamma, beta, psi_1, psi_2, theta, gcash, gfund, survprob(t), n, nweig1, t), ...
                lb, ub, grid_size_work);
                
            % 将结果存储到临时变量
            results(task_idx, :) = [policy(1), policy(2), policy(3), 1]; % 1表示优化成功
        end

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

        % 自适应步长
        alpha = 0.1;

        % 更新策略
        C(:, :, t) = (1-alpha) * C(:, :, t) + alpha * C_new;
        A(:, :, t) = (1-alpha) * A(:, :, t) + alpha * A_new;
        Q(:, :, t) = (1-alpha) * Q(:, :, t) + alpha * Q_new;
        
        % 计算值函数 - 逐点计算而不是整个矩阵
        V_new = zeros(ncash, nfund, 'like', V);
        for i3 = 1:ncash
            for i2 = 1:nfund
                % 计算当前状态的值函数
                V_new(i3, i2) = -gpu_baseline_fun_valuefunc_work([C(i3, i2, t), A(i3, i2, t), Q(i3, i2, t)], ...
                    gcash(i3), gfund(i2), gyp(:, :, t), V(:, :, t+1), yh, gret, rf, rp, pension_pct, ...
                    gamma, beta, psi_1, psi_2, theta, gcash, gfund, survprob(t), n, nweig1, t);
            end
        end
        V(:, :, t) = V_new;

        % 检查收敛
        if policy_change < tolerance
            fprintf('  工作期时期 %d 在 %d 次迭代后收敛\n', t, iter_count);
            break;
        end

        % 如果达到最大迭代次数仍未收敛
        if iter_count == max_iterations
            fprintf('  警告：工作期时期 %d 在 %d 次迭代后仍未收敛\n', t, max_iterations);
        end

        % 定期清理GPU内存
        if use_gpu && mod(iter, 10) == 0
            wait(gpuDevice);
        end
    end
    
    % 输出当前时期求解时间
    period_time = toc(period_start_time);
    if use_gpu
        fprintf('  工作期时期 %d 求解完成，耗时 %.2f 秒，GPU内存: %.2f GB\n', t, period_time, get_gpu_memory());
    else
        fprintf('  工作期时期 %d 求解完成，耗时 %.2f 秒\n', t, period_time);
    end
end

% 输出工作期总求解时间
work_time = toc(work_start_time);
fprintf('工作期求解完成，总耗时 %.2f 秒 (%.2f 分钟)\n', work_time, work_time/60);

% 输出总求解时间
total_time = retire_time + work_time;
fprintf('模型求解完成，总耗时 %.2f 秒 (%.2f 分钟)\n', total_time, total_time/60);

% 保存结果
if ~exist('result_baseline_matlab_PFI_GPU', 'dir')
    mkdir('result_baseline_matlab_PFI_GPU');
end

% 将GPU数组转换为CPU数组以保存
C_cpu = gather(C);
A_cpu = gather(A);
Q_cpu = gather(Q);
V_cpu = gather(V);
gcash_cpu = gather(gcash);
gfund_cpu = gather(gfund);

% 保存完整的三维策略函数和值函数
save('result_baseline_matlab_PFI_GPU/model_results.mat', 'C_cpu', 'A_cpu', 'Q_cpu', 'V_cpu', 'gcash_cpu', 'gfund_cpu', ...
    'tb', 'tr', 'td', 'tn', 'n', 'rf', 'rp', 'mu', 'sigr', 'smay', 'smav', ...
    'ret_fac', 'pension_pct', 'pension_rate', 'gy');

% 同时保存一个切片用于Excel查看
writematrix(gather(C(:, 1, :)), 'result_baseline_matlab_PFI_GPU/C_slice.csv');
writematrix(gather(A(:, 1, :)), 'result_baseline_matlab_PFI_GPU/A_slice.csv');
writematrix(gather(Q(:, 1, :)), 'result_baseline_matlab_PFI_GPU/Q_slice.csv');
writematrix(gather(V(:, 1, :)), 'result_baseline_matlab_PFI_GPU/V_slice.csv');
writematrix(gather(gcash), 'result_baseline_matlab_PFI_GPU/gcash.csv');
writematrix(gather(gfund), 'result_baseline_matlab_PFI_GPU/gfund.csv');

fprintf('模型求解完成，结果已保存到result_baseline_matlab_PFI_GPU文件夹\n');

% 返回求解结果，供模拟使用
model_results = struct();
model_results.C = C_cpu;
model_results.A = A_cpu;
model_results.Q = Q_cpu;
model_results.V = V_cpu;
model_results.gcash = gcash_cpu;
model_results.gfund = gfund_cpu;
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

% 如果使用了GPU，释放GPU内存
if use_gpu
    gpuDevice(1);
end

% 添加GPU内存监控函数
function mem = get_gpu_memory()
    if gpuDeviceCount > 0
        d = gpuDevice();
        mem = d.AvailableMemory/1024/1024/1024; % GB
    else
        mem = 0;
    end
end




