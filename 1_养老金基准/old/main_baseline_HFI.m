clear
% 使用Howard策略迭代方法求解模型的最优消费和投资决策
% 使用parfor在fund和cash层面并行计算

% 创建保存每期数据的目录
if ~exist('result_baseline_matlab_HFI', 'dir')
    mkdir('result_baseline_matlab_HFI');
    mkdir('result_baseline_matlab_HFI/slices');
    mkdir('result_baseline_matlab_HFI/slices/retire');
    mkdir('result_baseline_matlab_HFI/slices/work');
end

% 设置随机种子以获得可重复的结果
rng(42);

% 变量定义
tb = 18;    % 初始开始的年纪
tr = 61;    % 退休年龄
td = 100;   % 死亡年龄
tn = td - tb + 1;  % 总期数

% 状态变量grid数量
ncash = 51;  % 手中现金
nfund = 51;  % 养老基金余额
n = 3;       % 外生随机冲击的grid数量

% 外生参数
% 基础收入f的系数
aa = (-2.170042 + 2.700381);
b1 = 0.16818;
b2 = -0.0323371 / 10;
b3 = 0.0019704 / 100;

% 养老金相关
ret_fac = 0.6827;  % 退休后固定支付的工资（基本养老金）
pension_pct = 0.08;   % 工资扣除缴纳基本养老保险的比例（个人账户部分）
rp = 1.01;    % 个人养老基金的无风险收益率

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

% Howard策略迭代参数
max_iterations = 200;  % 最大迭代次数
max_value_iterations = 100;  % 每次策略更新后的值函数迭代次数
tolerance = 1e-6;      % 收敛容差
value_tolerance = 1e-6; % 值函数迭代收敛容差

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
num_cores = 14; % feature('numcores');
fprintf('使用 %d 个CPU核心进行并行计算\n', num_cores);

% 设置优化选项
options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp', ...
    'MaxFunctionEvaluations', 200, 'OptimalityTolerance', 1e-7);

% Howard策略迭代
% 退休期
fprintf('开始求解退休期值函数（使用Howard策略迭代）...\n');
retire_start_time = tic; % 开始计时退休期求解

for i1 = 1:(td - tr + 1)
    t = tn - i1;
    period_start_time = tic; % 开始计时当前时期
    fprintf('退休期求解进度: %d/%d\n', i1, td - tr + 1);

    % 初始化策略
    for i2 = 1:nfund
        for i3 = 1:ncash
            % 合理猜测
            C(i3, i2, t) = C(i3, i2, t+1);
            A(i3, i2, t) = A(i3, i2, t+1);
            Q(i3, i2, t) = 0.0;  % 退休期不需要购买个人养老金

            % 随机猜测
            % c_perturb = C(i3, i2, t+1) + (rand() - 0.5) * 0.4;
            % a_perturb = A(i3, i2, t+1) + (rand() - 0.5) * 0.4;
            % q_perturb = 0.1 + (rand() - 0.5) * 0.4;
            % % 确保值在合理边界内 [0, 1]
            % C(i3, i2, t) = max(0, min(1, c_perturb));
            % A(i3, i2, t) = max(0, min(1, a_perturb));
            % Q(i3, i2, t) = max(0, min(1, q_perturb));
        end
    end

    % 初始化迭代计数器
    iter_count = 0;
    
    % Howard策略迭代循环
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
            
            % 使用fmincon优化
                [policy, ~] = fmincon(@(x) baseline_fun_valuefunc_retired(x, gcash(i3), gfund(i2), ...
                    V(:, i2, t+1), gret, rf, ret_fac, pension_pay, ...
                    gamma, beta, psi_1, psi_2, theta, gcash, gfund, survprob(t), weig), ...
                    x0, [], [], [], [], lb, ub, [], options);
                
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

        % 更新策略 - Howard策略迭代不使用步长，直接更新
        C(:, :, t) = C_new;
        A(:, :, t) = A_new;
        
        % Howard策略迭代的核心：多次值函数迭代
        V_old = zeros(ncash, nfund);
        
        % 值函数迭代循环
        for value_iter = 1:max_value_iterations
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
            
            % 计算值函数变化量
            value_change = max(max(abs(V_new - V_old)));
            V_old = V_new;
            fprintf('值函数误差%d \n',value_change)
            % 如果值函数收敛，提前退出值函数迭代
            if value_change < value_tolerance
                break;
            end
        end
        
        % 更新值函数
        V(:, :, t) = V_old;

        % 检查策略收敛
        if (policy_change < tolerance) %&& (policy_change>1e-12)
            
            fprintf('  退休期时期 %d 在 %d 次迭代后收敛, 误差为 %d \n', t, iter_count, policy_change);
            break;
        end

        % 如果达到最大迭代次数仍未收敛
        if iter_count == max_iterations
            fprintf('  警告：退休期时期 %d 在 %d 次迭代后仍未收敛\n', t, max_iterations);
        end
    end
    
    % 保存当前时期的切片数据
    slice_dir = sprintf('result_baseline_matlab_HFI/slices/retire/period_%d', t);
    if ~exist(slice_dir, 'dir')
        mkdir(slice_dir);
    end

    writematrix(C(:, :, t), sprintf('%s/C_slice_t%d.csv', slice_dir, t));
    writematrix(A(:, :, t), sprintf('%s/A_slice_t%d.csv', slice_dir, t));
    writematrix(Q(:, :, t), sprintf('%s/Q_slice_t%d.csv', slice_dir, t));
    writematrix(V(:, :, t), sprintf('%s/V_slice_t%d.csv', slice_dir, t));

    % 保存为.mat文件
    save(sprintf('%s/data_t%d.mat', slice_dir, t), 'C', 'A', 'Q', 'V');
    
    % 输出当前时期求解时间
    period_time = toc(period_start_time);
    fprintf('  退休期时期 %d 求解完成，耗时 %.2f 秒\n', t, period_time);
end

% 输出退休期总求解时间
retire_time = toc(retire_start_time);
fprintf('退休期求解完成，总耗时 %.2f 秒 (%.2f 分钟)\n', retire_time, retire_time/60);

% 工作期（退休前）
fprintf('开始求解工作期值函数（使用Howard策略迭代）...\n');
work_start_time = tic; % 开始计时工作期求解

for i1 = 1:(tr - tb)
    t = tr - tb - i1 + 1;
    period_start_time = tic; % 开始计时当前时期
    fprintf('工作期求解进度: %d/%d\n', i1, tr - tb);

    % 初始化策略
    for i2 = 1:nfund
        for i3 = 1:ncash
            % 合理猜测
            C(i3, i2, t) = C(i3, i2, t+1);
            A(i3, i2, t) = A(i3, i2, t+1);
            Q(i3, i2, t) = 0.1;

                        % 随机猜测
            % c_perturb = C(i3, i2, t+1) + (rand() - 0.5) * 0.4;
            % a_perturb = A(i3, i2, t+1) + (rand() - 0.5) * 0.4;
            % q_perturb = 0.1 + (rand() - 0.5) * 0.4;
            % % 确保值在合理边界内 [0, 1]
            % C(i3, i2, t) = max(0, min(1, c_perturb));
            % A(i3, i2, t) = max(0, min(1, a_perturb));
            % Q(i3, i2, t) = max(0, min(1, q_perturb));
        end
    end

    % 初始化迭代计数器
    iter_count = 0;
    
    % Howard策略迭代循环
    for iter = 1:max_iterations
        iter_count = iter;  % 更新迭代计数
        
        % 策略改进 - 基于当前值函数更新策略，使用并行计算
        C_new = zeros(ncash, nfund);
        A_new = zeros(ncash, nfund);
        Q_new = zeros(ncash, nfund);

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
        results = zeros(size(tasks, 1), 4); % 用于存储结果的数组 [C, A, Q, 优化状态]
        
        parfor task_idx = 1:size(tasks, 1)
            i3 = tasks(task_idx, 1);
            i2 = tasks(task_idx, 2);
            
            x0 = [C(i3, i2, t), A(i3, i2, t), Q(i3, i2, t)];
            lb = [0, 0, 0];
            ub = [1, 1, 1];
            
            % 确保初始值在边界范围内
            x0 = max(min(x0, ub), lb);
            
            % 使用fmincon优化
                [policy, ~] = fmincon(@(x) baseline_fun_valuefunc_work(x, gcash(i3), gfund(i2), ...
                    gyp(:, :, t), V(:, :, t+1), yh, gret, rf, rp, pension_pct, ...
                    gamma, beta, psi_1, psi_2, theta, gcash, gfund, survprob(t), n, nweig1), ...
                    x0, [], [], [], [], lb, ub, [], options);
                
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

        % 更新策略 - Howard策略迭代不使用步长，直接更新
        C(:, :, t) = C_new;
        A(:, :, t) = A_new;
        Q(:, :, t) = Q_new;
        
        % Howard策略迭代的核心：多次值函数迭代
        V_old = zeros(ncash, nfund);
        
        % 值函数迭代循环
        for value_iter = 1:max_value_iterations
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
            
            % 计算值函数变化量
            value_change = max(max(abs(V_new - V_old)));
            V_old = V_new;
            
            % 如果值函数收敛，提前退出值函数迭代
            if value_change < value_tolerance
                break;
            end
        end
        
        % 更新值函数
        V(:, :, t) = V_old;

        % 检查策略收敛
        if policy_change < tolerance
            fprintf('  工作期时期 %d 在 %d 次迭代后收敛\n, 误差为 %d \n', t, iter_count, policy_change);
            break;
        end

        % 如果达到最大迭代次数仍未收敛
        if iter_count == max_iterations
            fprintf('  警告：工作期时期 %d 在 %d 次迭代后仍未收敛\n', t, max_iterations);
        end
    end
    
    % 保存当前时期的切片数据
    slice_dir = sprintf('result_baseline_matlab_HFI/slices/work/period_%d', t);
    if ~exist(slice_dir, 'dir')
        mkdir(slice_dir);
    end

    writematrix(C(:, :, t), sprintf('%s/C_slice_t%d.csv', slice_dir, t));
    writematrix(A(:, :, t), sprintf('%s/A_slice_t%d.csv', slice_dir, t));
    writematrix(Q(:, :, t), sprintf('%s/Q_slice_t%d.csv', slice_dir, t));
    writematrix(V(:, :, t), sprintf('%s/V_slice_t%d.csv', slice_dir, t));

    % 保存为.mat文件
    save(sprintf('%s/data_t%d.mat', slice_dir, t), 'C', 'A', 'Q', 'V');
    
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
if ~exist('result_baseline_matlab_HFI', 'dir')
    mkdir('result_baseline_matlab_HFI');
end

% 保存完整的三维策略函数和值函数
save('result_baseline_matlab_HFI/model_results.mat', 'C', 'A', 'Q', 'V', 'gcash', 'gfund', ...
    'tb', 'tr', 'td', 'tn', 'n', 'rf', 'rp', 'mu', 'sigr', 'smay', 'smav', ...
    'ret_fac', 'pension_pct', 'pension_rate', 'gy');

% 同时保存一个切片用于Excel查看
writematrix(C(:, 1, :), 'result_baseline_matlab_HFI/C_slice.csv');
writematrix(A(:, 1, :), 'result_baseline_matlab_HFI/A_slice.csv');
writematrix(Q(:, 1, :), 'result_baseline_matlab_HFI/Q_slice.csv');
writematrix(V(:, 1, :), 'result_baseline_matlab_HFI/V_slice.csv');
writematrix(gcash, 'result_baseline_matlab_HFI/gcash.csv');
writematrix(gfund, 'result_baseline_matlab_HFI/gfund.csv');

fprintf('模型求解完成，结果已保存到result_baseline_matlab_HFI文件夹\n');

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



