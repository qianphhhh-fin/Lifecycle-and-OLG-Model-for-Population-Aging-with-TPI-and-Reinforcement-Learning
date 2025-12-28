function batch_baseline_VFI(varargin)
% batch_baseline_VFI - 求解养老金模型并可通过命令行参数修改关键参数
%
% 使用方法:
% matlab -nosplash -nodesktop -r "batch_baseline_VFI('tb',20,'tr',65,'rf',1.01);exit"
%
% 可设置参数:
% 'tb' - 初始开始年龄 (默认: 18)
% 'tr' - 退休年龄 (默认: 61)
% 'td' - 死亡年龄 (默认: 100)
% 'n' - 外生随机冲击的grid数量 (默认: 5)
% 'ncash' - 手中现金网格数 (默认: 51)
% 'nfund' - 养老基金账户余额网格数 (默认: 51)
% 'ret_fac' - 退休后固定支付的工资/基本养老金 (默认: 0)
% 'pension_pct' - 工资扣除缴纳养老保险的比例 (默认: 0)
% 'rp' - 个人养老基金的无风险收益率 (默认: 1.04)
% 'smay' - 白噪声shock的标准差 (默认: sqrt(0.169993))
% 'smav' - 持续shock的标准差 (默认: sqrt(0.112572))
% 'gamma' - Epstein-Zin的相对风险规避系数 (默认: 3.84)
% 'beta' - Epstein-Zin的discount factor (默认: 0.95)
% 'sigr' - 风险资产收益率的标准差 (默认: 0.27)
% 'rf' - 无风险总收入 (默认: 1.02)
% 'mu' - 超额收益 (默认: 0.04)
% 'save_path' - 保存结果的路径 (默认: 'result/baseline_no_ret_fac/')

% 命令行参数解析
p = inputParser;
addParameter(p, 'tb', 18, @isnumeric);          % 初始开始年龄
addParameter(p, 'tr', 61, @isnumeric);          % 退休年龄
addParameter(p, 'td', 100, @isnumeric);         % 死亡年龄
addParameter(p, 'n', 5, @isnumeric);            % 外生随机冲击的grid数量
addParameter(p, 'ncash', 51, @isnumeric);       % 手中现金网格数
addParameter(p, 'nfund', 51, @isnumeric);       % 养老基金账户余额网格数
addParameter(p, 'ret_fac', 0, @isnumeric);      % 退休后固定支付的工资（基本养老金）
addParameter(p, 'pension_pct', 0, @isnumeric);  % 工资扣除缴纳养老保险的比例
addParameter(p, 'rp', 1.04, @isnumeric);        % 个人养老基金的无风险收益率
addParameter(p, 'smay', sqrt(0.169993), @isnumeric); % 白噪声shock的标准差
addParameter(p, 'smav', sqrt(0.112572), @isnumeric); % 持续shock的标准差
addParameter(p, 'corr_z_epsilon', 0, @isnumeric);  % 工资收入白噪声与风险收益随机项的相关性
addParameter(p, 'corr_u_epsilon', 0, @isnumeric);  % 工资收入AR(1)随机项与风险收益随机项的相关性
addParameter(p, 'gamma', 3.84, @isnumeric);     % Epstein-Zin的相对风险规避系数
addParameter(p, 'beta', 0.95, @isnumeric);      % Epstein-Zin的discount factor
addParameter(p, 'sigr', 0.27, @isnumeric);      % 风险资产收益率的标准差
addParameter(p, 'rf', 1.02, @isnumeric);        % 无风险总收入
addParameter(p, 'mu', 0.04, @isnumeric);        % 超额收益
addParameter(p, 'save_path', 'result/sensitivity/', @ischar); % 保存路径
addParameter(p, 'filename', 'model_results.mat', @ischar); % 保存路径


parse(p, varargin{:});

% 只打印非默认参数，不输出换行符
% fprintf('批处理: ');
% defaults = {'tb', 18; 'tr', 61; 'td', 100; 'n', 5; 'ncash', 51; 'nfund', 51; ...
%     'ret_fac', 0; 'pension_pct', 0; 'rp', 1.04; 'smay', sqrt(0.169993); ...
%     'smav', sqrt(0.112572); 'gamma', 3.84; 'beta', 0.95; 'sigr', 0.27; ...
%     'rf', 1.02; 'mu', 0.04; 'save_path', 'result/sensitivity/'};
% 
% hasNonDefault = false;
% param_string = '';
% for i = 1:size(defaults, 1)
%     param = defaults{i, 1};
%     defaultVal = defaults{i, 2};
% 
%     % 处理数值参数
%     if isnumeric(defaultVal)
%         currentVal = p.Results.(param);
%         if abs(currentVal - defaultVal) > 1e-10
%             if hasNonDefault
%                 param_string = [param_string, ', '];
%             end
%             param_string = [param_string, sprintf('%s=%.4g', param, currentVal)];
%             hasNonDefault = true;
%         end
%         % 处理字符串参数
%     elseif ischar(defaultVal)
%         currentVal = p.Results.(param);
%         if ~strcmp(currentVal, defaultVal)
%             if hasNonDefault
%                 param_string = [param_string, ', '];
%             end
%             % 对路径特殊处理，避免打印过长
%             if strcmp(param, 'save_path')
%                 [~, folder, ~] = fileparts(currentVal);
%                 param_string = [param_string, sprintf('%s=%s', param, folder)];
%             else
%                 param_string = [param_string, sprintf('%s=%s', param, currentVal)];
%             end
%             hasNonDefault = true;
%         end
%     end
% end
% 
% if ~hasNonDefault
%     param_string = '默认参数';
% end
% 
% % 截断过长的参数字符串以确保显示在一行
% if length(param_string) > 40
%     param_string = [param_string(1:37), '...'];
% end
% 
% % 将参数字符串保存，用于进度条显示
% params_display = param_string;
% 
% % 设置保存路径
% save_path = p.Results.save_path;

% close all;

save_path = p.Results.save_path;
filename = p.Results.filename;
%% 模型求解和模拟
% 创建保存每期数据的目录
% if ~exist(save_path, 'dir')
%     mkdir(save_path);
% end

% 设置随机种子以获得可重复的结果
rng(42);

% 变量定义
tb = p.Results.tb;    % 初始开始的年纪
tr = p.Results.tr;    % 退休年龄
td = p.Results.td;    % 死亡年龄
tn = td - tb + 1;     % 总期数

% 状态变量grid数量 - 匹配Python代码
ncash = p.Results.ncash;  % 手中现金
nfund = p.Results.nfund;  % 养老基金账户余额 - 新增状态变量
n = p.Results.n;          % 外生随机冲击的grid数量

% 外生参数
% 基础收入f的系数
aa = (-2.170042 + 2.700381);
b1 = 0.16818;
b2 = -0.0323371 / 10;
b3 = 0.0019704 / 100;

% 养老金相关 - 仅保留基本养老金
ret_fac = p.Results.ret_fac;      % 退休后固定支付的工资（基本养老金）
pension_pct = p.Results.pension_pct;   % 工资扣除缴纳养老保险的比例，简化版设为0
rp = p.Results.rp;                % 个人养老基金的无风险收益率 - 新增参数

% 随机冲击参数
smay = p.Results.smay;  % 白噪声shock的标准差
smav = p.Results.smav;  % 持续shock的标准差
corr_z_epsilon = p.Results.corr_z_epsilon;     % 工资持续收入随机项与风险收益随机项的相关性
corr_u_epsilon = p.Results.corr_u_epsilon;     % 工资临时收入白噪声与风险收益随机项的相关性

% 效用函数参数
gamma = p.Results.gamma;   % Epstein-Zin的相对风险规避系数
beta = p.Results.beta;     % Epstein-Zin的discount factor
psi = 0.15;                % Epstein-Zin的跨期替代弹性

% 资产收益率参数
rf = p.Results.rf;     % 无风险总收入
mu = p.Results.mu;     % 超额收益
sigr = p.Results.sigr; % 风险资产收益率的标准差

% 确定可用CPU核心数
num_cores = 20; % feature('numcores');

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
gfund = zeros(nfund, 1);  % 养老基金网格 - 新增
lgfund = zeros(nfund, 1); % 养老基金网格的对数 - 新增
secd = zeros(ncash, 1);

% 扩展为三维数组 - 增加养老基金维度
C = zeros(ncash, nfund,  tn);
V = zeros(ncash, nfund,  tn);
A = ones(ncash, nfund, tn);
Q = zeros(ncash, nfund, tn);  % 个人养老金购买决策 - 新增

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

% 养老基金grid - 新增
maxfund = 100;
minfund = 0.01;
l_maxfund = log(maxfund);
l_minfund = log(minfund);
stepfund = (l_maxfund - l_minfund) / (nfund - 1);

for i1 = 1:nfund
    lgfund(i1) = l_minfund + (i1-1) * stepfund;
    gfund(i1) = exp(lgfund(i1));
end
gfund(1) = 0;

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

% 计算个人养老金给付 - 新增
% 假设退休后按照剩余年份平均分配养老基金余额
pension_period = td - tr + 1;
pension_rate = 1 / (pension_period -20);

% 终止期 - 扩展到三维
    for i1 = 1:ncash
        for i2 = 1:nfund
            C(i1, i2, tn) = 1;    % 最后一期的现金全部用来消费
            A(i1, i2, tn) = 0.0;          % 最后一期的投资全部为零
            Q(i1, i2, tn) = 0.0;          % 最后一期不购买养老金
            % V(i1, i2, tn) = (gcash(i1)) * (1 - beta)^psi_2; % 与Python代码一致的效用函数
            V(i1, i2,tn) = -(gcash(i1))^(1-gamma); % 与Python代码一致的效用函数
        end
    end

% 设置优化选项 - 使优化参数接近Python的L-BFGS-B
options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp', ...
    'MaxFunctionEvaluations', 1000, 'MaxIterations', 400, ...
    'OptimalityTolerance', 1e-10, 'StepTolerance', 1e-10, ...
    'ConstraintTolerance', 1e-10);

%% 模型求解 - 使用策略函数迭代方法
% 输出带有进度条的开始信息 - 参数和进度条在同一行
% fprintf('%s - 进度: [                              ]  0%%', params_display);

% 退休期
retire_start_time = tic; % 开始计时退休期求解
total_steps = (td - tr) + (tr - tb);
progress_percentage = 0;

for i1 = 1:(td - tr)
    t = tn - i1;
    period_start_time = tic; % 开始计时当前时期

        % 创建临时数组存储结果
        temp_C = zeros(ncash, nfund);
        temp_A = zeros(ncash, nfund);
        temp_V = zeros(ncash, nfund);
        
        % 只在养老基金维度上并行
        parfor i2 = 1:nfund
            % 为当前养老基金水平创建临时数组
            local_C = zeros(ncash, 1);
            local_A = zeros(ncash, 1);
            local_V = zeros(ncash, 1);
            
            % 计算养老金给付
            pension_pay = pension_rate * gfund(i2);
            
            % 对每个现金水平进行循环
            for i3 = 1:ncash
                x0 = [0.5, 0.5];
                lb = [0, 0];
                ub = [1, 1];

                % 确保初始值在边界范围内
                x0(1) = max(min(x0(1), ub(1)), lb(1));
                x0(2) = max(min(x0(2), ub(2)), lb(2));

                % 优化函数，使用fmincon模拟Python的minimize函数
                [policy, fval, exitflag] = fmincon(@(x) fun_valuefunc_retired(x, gcash(i3), gfund(i2), squeeze(V(:, i2,t+1))...
                    , gret, rf, ret_fac, pension_pay, gamma, beta, psi_1, psi_2, theta, gcash, gfund, survprob(t), weig), ...
                    x0, [], [], [], [], lb, ub, [], options);
                
                % 存储到本地临时数组
                local_C(i3) = policy(1);
                local_A(i3) = policy(2);
                local_V(i3) = -fval;
            end
            
            % 将本地结果复制到临时结果数组
            temp_C(:, i2) = local_C;
            temp_A(:, i2) = local_A;
            temp_V(:, i2) = local_V;
        end
        
        % 将临时结果复制到主变量
            C(:, :, t) = temp_C;
                    temp_A(temp_C>=0.99990)=1;
        temp_A(temp_C<=1-0.99990)=0;
            A(:, :,  t) = temp_A;
            V(:, :,  t) = temp_V;

    % 更新进度条 - 使用\r回到行首，整行重新显示参数和进度条
    % progress_percentage = progress_percentage + 1/total_steps * 100;
    % progress_bar = repmat(' ', 1, 30);
    % filled = round(progress_percentage * 30 / 100);
    % progress_bar(1:filled) = '=';
    % if filled > 0
    %     progress_bar(filled) = '>';
    % end
    % fprintf('\r%s - 进度: [%s] %3.0f%%', params_display, progress_bar, progress_percentage);
end

options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp', ...
    'MaxFunctionEvaluations', 1000, 'MaxIterations', 400, ...
    'OptimalityTolerance', 1e-10, 'StepTolerance', 1e-10, ...
    'ConstraintTolerance', 1e-10);

% 工作期（退休前）
work_start_time = tic; % 开始计时工作期求解

for i1 = 1:(tr - tb)
    t = tr - tb - i1 +1;
    period_start_time = tic; % 开始计时当前时期

        % 设置并行池

    

        % 在nfund层面使用parfor并行计算
        temp_C = zeros(ncash, nfund);
        temp_A = zeros(ncash, nfund);
        temp_Q = zeros(ncash, nfund);
        temp_V = zeros(ncash, nfund);  

        parfor i2 = 1:nfund
            % 创建临时数组存储当前fund层级的结果
            temp_C_i2 = zeros(ncash, 1);
            temp_A_i2 = zeros(ncash, 1);
            temp_Q_i2 = zeros(ncash, 1);
            temp_V_i2 = zeros(ncash, 1);
            
            for i3 = 1:ncash
                if i1==1
                     x0 = [ C(i3,i2,t+1), A(i3,i2,t+1), 0.5];
                else
                x0 = [ C(i3,i2,t+1), A(i3,i2,t+1), Q(i3,i2,t+1)];
                end
                lb = [0, 0, 0];
                ub = [1, 1, 1];

                % 确保初始值在边界范围内
                x0(1) = max(min(x0(1), ub(1)), lb(1));
                x0(2) = max(min(x0(2), ub(2)), lb(2));
                x0(3) = max(min(x0(3), ub(3)), lb(3));
                
                [policy, fval, exitflag] = fmincon(@(x) fun_valuefunc_work(x, gcash(i3), gfund(i2), gyp(:, :, t),...
                    squeeze(V(:, :, t+1)), yh, gret, rf, rp, pension_pct, gamma, beta, psi_1, psi_2, theta, gcash, gfund,...
                    survprob(t), n, nweig1), ...
                    x0, [], [], [], [], lb, ub, [], options);

                % 存储到临时数组
                temp_C_i2(i3) = policy(1);
                temp_A_i2(i3) = policy(2);
                temp_Q_i2(i3) = policy(3);
                temp_V_i2(i3) = -fval;
            end
            
            % 将当前fund层级的结果存储到主变量
            temp_C(:, i2) = temp_C_i2;
            temp_A(:, i2) = temp_A_i2;
            temp_Q(:, i2) = temp_Q_i2;
            temp_V(:, i2) = temp_V_i2;
        end

            % 将临时结果复制到主变量
        C(:, :, t) = temp_C;
        temp_A(temp_C>=0.99990 | temp_Q>=0.99990)=1;
        A(:, :,  t) = temp_A;
        temp_Q(temp_C>=0.99990)=0;
        Q(:, :,  t) = temp_Q;
        V(:, :, t) = temp_V;
    
    % 更新进度条 - 使用\r回到行首，整行重新显示参数和进度条
    % progress_percentage = progress_percentage + 1/total_steps * 100;
    % progress_bar = repmat(' ', 1, 30);
    % filled = round(progress_percentage * 30 / 100);
    % progress_bar(1:filled) = '=';
    % if filled > 0
    %     progress_bar(filled) = '>';
    % end
    % fprintf('\r%s - 进度: [%s] %3.0f%%', params_display, progress_bar, progress_percentage);
end

% 结束进度条，添加完成消息和保存路径
% fprintf('\r%s - 结果: [%s] 已保存到 %s\n', params_display, repmat('=', 1, 30), save_path);

% 保存一些二维切片，便于查看
% writematrix(C(:, 1, :), fullfile(save_path, 'C_slice.csv'));
% writematrix(A(:, 1, :), fullfile(save_path, 'A_slice.csv'));
% writematrix(Q(:, 1, :), fullfile(save_path, 'Q_slice.csv'));
% writematrix(V(:, 1, :), fullfile(save_path, 'V_slice.csv'));
% writematrix(gcash, fullfile(save_path, 'gcash.csv'));
% writematrix(gfund, fullfile(save_path, 'gfund.csv'));

% 返回求解结果，供模拟使用
model_results = struct();
model_results.C = C;
model_results.A = A;
model_results.Q = Q;  % 增加养老金购买决策
model_results.V = V;
model_results.gcash = gcash;
model_results.gfund = gfund;  % 增加养老基金网格
model_results.params = struct();
model_results.params.tb = tb;
model_results.params.tr = tr;
model_results.params.td = td;
model_results.params.tn = tn;
model_results.params.n = n;
model_results.params.rf = rf;
model_results.params.rp = rp;  % 增加养老基金收益率
model_results.params.mu = mu;
model_results.params.sigr = sigr;
model_results.params.smay = smay;
model_results.params.smav = smav;
model_results.params.ret_fac = ret_fac;
model_results.params.corr_z_epsilon = corr_z_epsilon;
model_results.params.corr_u_epsilon = corr_u_epsilon;
model_results.params.pension_pct = pension_pct;  % 增加养老金缴费比例
model_results.params.pension_rate = pension_rate;  % 增加养老金支付比例
model_results.params.gy = gy;

% 保存完整的三维策略函数和值函数
save(fullfile(save_path, filename), 'model_results');

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

%% 退休期间的值函数计算
% function v = cocco_fun_valuefunc_retired(x, cash, fund, nextV, gret, rf, ret_fac, pension_pay, gamma, beta, psi_1, psi_2, theta, gcash, gfund, survprob, weig)
% % 退休期间的值函数计算
% %
% % 参数:
% % x - 决策变量 [消费, 风险资产投资比例]
% % cash - 手中现金
% % fund - 养老基金账户余额
% % nextV - 下期值函数
% % gret - 随机收益率
% % rf - 无风险收益率
% % ret_fac - 固定养老金
% % pension_pay - 个人养老金给付
% % gamma - 相对风险规避系数
% % beta - 贴现因子
% % psi_1, psi_2, theta - Epstein-Zin效用函数参数
% % gcash - 现金网格
% % gfund - 养老基金网格
% % survprob - 存活概率
% % weig - 随机状态权重
% %
% % 返回:
% % v - 函数值
%
% % 获取网格大小
% ncash = length(gcash);
% nfund = length(gfund);
%
% auxVV = 0;
% sav = cash - x(1);
%
% % 与Python一致的效用函数
% u = -(x(1))^(1-gamma);
%
% % 下期的cash-on-hand
% cash_1 = (rf * (1 - x(2)) + gret .* x(2)) .* sav + ret_fac + pension_pay;
%
% % 限制cash_1的范围
% cash_1 = max(min(cash_1, gcash(end)), gcash(1));
%
% % 下期养老基金余额 - 保持不变，因为退休期间不再购买养老金
% fund_1 = fund;
%
% % 使用双线性插值计算下期值函数
% % 创建网格
% [X, Y] = meshgrid(gcash, gfund);
% % 将nextV转换为适合interp2的格式
% V_grid = reshape(nextV, [ncash, nfund])';
% % 对每个可能的随机收益计算下期值函数
% int_V = zeros(size(gret));
% for i = 1:length(gret)
%     int_V(i) = interp2(X, Y, V_grid, cash_1(i), fund_1, 'spline', 'extrap');
% end
%
% % 计算期望值
% auxVV = auxVV + weig' * (survprob * int_V);
%
% % 和Python一致的返回值
% if auxVV == 0
%     v = -(u);
% else
%     v = -(u + beta * auxVV);
% end
% end

    function v =  fun_valuefunc_retired(x, cash, fund, nextV, gret, rf, ret_fac, pension_pay, gamma, beta, psi_1, psi_2, theta, gcash, gfund, survprob, weig, t)
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
        sav = cash * (1  - x(1));

        % Epstein-Zin效用函数
        % u = (x(1) +0.00001).^psi_1;  % ez效用
        u = -(x(1) *cash)^(1-gamma);  % CRRA效用


        % 下期的cash-on-hand
        cash_1 = (rf * (1 - x(2)) + gret .* x(2)) .* sav + ret_fac + pension_pay;

        % 限制cash_1的范围
        cash_1 = max(min(cash_1, gcash(end)), gcash(1));

        % 使用MATLAB的一维插值
        int_V = interp1(gcash, nextV, cash_1,  'spline', 'extrap');
        % int_V = max(int_V, 1e-20); % 确保插值结果为正

        % 计算期望值
        % auxVV = auxVV + weig' * (survprob * int_V.^(1-gamma));
        auxVV = auxVV + weig' * survprob * int_V;

        % Epstein-Zin递归效用
        % v = -(((1-beta)*u + beta*((auxVV)^(1/theta)))^psi_2);
        v = -(u + beta*auxVV);

        % 确保返回值是有效的双精度数
    end

%% 工作期间的值函数计算

    function v = fun_valuefunc_work(x, cash, fund, gyp, nextV, yh, gret, rf, rp, tau, gamma, beta, psi_1, psi_2, theta, gcash, gfund, survprob, n, weig, t)
        % 工作期间的值函数计算
        %
        % 参数:
        % x - 决策变量 [消费比例, 风险资产投资比例, 个人养老金购买比例]
        % cash - 手中现金
        % fund - 养老基金账户余额
        % gyp - 劳动收入增长率
        % nextV - 下期值函数
        % yh - 暂时收入冲击
        % yh_realized - 当前已实现的临时收入冲击值
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

        % 计算消费和养老金购买
        consumption = x(1) * cash;  % 消费绝对值 = 消费比例 * 现金
        pension_contrib = x(3) * (1-x(1)) * cash;  % 个人养老金购买

        % 效用函数
        u = -(consumption)^(1-gamma);

        % 获取维度
        % n1 = min(n, size(gyp, 1));  % z_t维度
        % n2 = min(n, size(gyp, 2));  % eta_t维度
        % n3 = n;  % epsilon_t维度

        % 创建索引网格
        % [I6, I7, I8] = ndgrid(1:n1, 1:n2, 1:n3);
        % I6 = permute(I6, [1, 2, 3]); % 确保维度顺序正确
        % I7 = permute(I7, [1, 2, 3]);
        % I8 = permute(I8, [1, 2, 3]);

        % 预分配数组
        % sav_values = zeros(n1, n2, n3);
        % cash_1_values = zeros(n1, n2, n3);
        % fund_1_values = zeros(n1, n2, n3);
        % int_V_values = zeros(n1, n2, n3);
        %
        % % 循环计算各状态下的值
        % for i6 = 1:n1
        %     for i7 = 1:n2
        %         for i8 = 1:n3
        %             % 获取当前状态的收入增长率
        %             gyp_value = gyp(i6, i8);
        %
        %             % 计算储蓄值
        %             sav_values(i6, i7, i8) = cash * (1 - x(1)) * (1 - x(3)) / gyp_value;
        %
        %             % 计算投资组合收益率
        %             portfolio_return = rf * (1 - x(2)) + gret(i8) * x(2);
        %
        %             % 获取暂时收入冲击
        %             yh_value = yh(i7, i8);
        %
        %             % 计算下期现金值
        %             cash_1_values(i6, i7, i8) = portfolio_return * sav_values(i6, i7, i8) + (1-tau) * yh_value;
        %
        %             % 计算下期养老基金账户余额
        %             fund_1_values(i6, i7, i8) = (fund + pension_contrib) * rp / gyp_value + tau * yh_value;
        %         end
        %     end
        % end

        cash_value = repmat(cash, [n, n, n]);
        gyp_value = repmat(reshape(gyp, [n, 1, n]), [1, n, 1]); % 首先将gyp转换为三维数组[n1, 1, n3]，然后在第二维复制n2次
        gret_value = repmat(reshape(gret, [1,1,n]), [n, n, 1]); % 首先将gret转换为三维数组[1, 1, n3]，然后在第1和第2维复制n2次
        yh_value = repmat(reshape(yh, [1, n, n]), [n, 1, 1]); % 首先将转换为三维数组[1, n2, n3]，然后在第1维复制n1次


        sav_values = cash_value * (1 - x(1)) * (1 - x(3)) ./ gyp_value;
        portfolio_return = rf * (1 - x(2)) + gret_value * x(2);
        cash_1_values = portfolio_return .* sav_values + (1-tau) * yh_value;
        fund_1_values = (fund + pension_contrib) * rp ./ gyp_value + tau * yh_value;



        % 限制cash_1和fund_1的范围
        cash_1_values = max(min(cash_1_values, gcash(end)), gcash(1));
        fund_1_values = max(min(fund_1_values, gfund(end)), gfund(1));

        % % 使用meshgrid为插值准备网格
        % [X, Y] = meshgrid(gcash, gfund);
        %
        % % 对每个状态组合进行二维插值
        % int_V_values = min(interp2(X, Y, nextV', cash_1_values, fund_1_values, 'spline'),0);

        int_V_values = min(interp2(gfund, gcash, nextV, fund_1_values, cash_1_values,  'spline'),0);

        % 计算加权和
        auxVV = sum(weig .* (survprob * int_V_values .* gyp_value),'all');

        % weighted_sum = 0;
        % for i6 = 1:n1
        %     for i7 = 1:n2
        %         for i8 = 1:n3
        %             % 获取权重 - 修正索引错误
        %             weight = weig(i6, i7, i8);
        %             weighted_sum = weighted_sum + weight * survprob * int_V_values(i6, i7, i8) * gyp(i6, i8);
        %         end
        %     end
        % end
        % auxVV = weighted_sum;

        % 返回值函数
        % if auxVV >0
        %     'v>0!'
        % else
        v = -(u + beta * auxVV);
    end