%% 生命周期模型求解与模拟 (VFI方法, CRRA效用)
% 描述:
% 本脚本求解并模拟了一个包含内生储蓄和投资决策的生命周期模型。
% 模型特点:
% - CRRA 效用函数
% - 持久性和暂时性收入冲击
% - 工作期和退休期
% - 无风险和风险资产投资
% - 工资税 (tau_y)
%
% 结构:
% 1. setup_parameters: 初始化所有模型参数、网格和外生过程。
% 2. solve_model_vfi: 使用值函数迭代(VFI)逆向求解模型。
% 3. simulate_model: 基于求解的策略函数进行前向蒙特卡洛模拟。
%
% 作者: Gemini (根据用户需求修改)
% 日期: 2025-09-15

% --- 初始化 ---
clear;
close all;
clc;

% --- 主流程控制器 ---
fprintf('开始执行生命周期模型 (CRRA 效用)...\n');

% 1. 设置所有模型参数
cS = setup_parameters();
save_path = cS.save_path;
% ** 修改: 更改保存文件名以反映CRRA效用 **
file_path = fullfile(save_path, 'vfi_benchmark_results.mat');

% 2. 求解模型 (值函数迭代)
fprintf('\n===== 步骤 2: 求解模型 (VFI) =====\n');
% if ~exist(file_path,'file')
    vfi_results = solve_model_vfi(cS);
    if ~exist(save_path, 'dir'), mkdir(save_path); end
    save(file_path, 'vfi_results', 'cS');
    fprintf('模型求解结果已保存到 %s\n', file_path);
% else
%     load(file_path)
% end

% 3. 模拟模型
fprintf('\n===== 步骤 3: 模拟模型 =====\n');
simulate_model(vfi_results, cS);

fprintf('\n模型执行完毕。\n');



%% =====================================================================
%                       1. 参数设置函数
%  =====================================================================
function cS = setup_parameters()
fprintf('===== 步骤 1: 设置模型参数 =====\n');

cS = struct();

% --- A. 生命周期与人口结构 ---
cS.tb = 22;       % 初始工作年龄
cS.tr = 26;       % 退休年龄
cS.td = 30;      % 最高寿命
cS.tn = cS.td - cS.tb + 1; % 总期数

% --- B. 经济主体偏好 (CRRA) ---
cS.beta = 0.95;     % 折现因子
cS.gamma = 3.84;    % 相对风险规避系数
% ** CRRA 修改: 移除psi, 因为 psi = 1/gamma **
% cS.psi = 0.15;      % 跨期替代弹性 (Epstein-Zin specific)

% --- C. 收入过程 ---
% 确定性部分 g(t) = exp(aa + b1*t + b2*t^2 + b3*t^3)
cS.aa = -2.170042 + 2.700381;
cS.b1 = 0.16818;
cS.b2 = -0.0323371 / 10;
cS.b3 = 0.0019704 / 100;
% 随机部分
cS.smay =  0.1; %sqrt(0.169993); %; % 暂时冲击标准差 (u_t)
cS.smav =  0.1; % sqrt(0.112572); % ; % 持久冲击标准差 (z_t)
cS.corr_z_epsilon = 0; % 持久冲击与风险收益的相关性 (未使用)
cS.corr_u_epsilon = 0; % 暂时冲击与风险收益的相关性 (未使用)

% --- D. 资产与回报率 ---
cS.rf = 1.02;     % 无风险总回报率
cS.mu = 0.04;     % 风险资产超额回报率
cS.sigr = 0.27; %0.27;   % 风险资产回报率标准差

% --- E. 养老金与税收 ---
cS.ret_fac = 0.6827;      % 退休后基础养老金 (替代率)
cS.tau_y = 0;          % 工资税率

% --- F. 数值求解参数 ---
cS.n_shocks = 8;      % 离散化的随机冲击节点数
cS.ncash = 101;        % 现金持有量网格数
cS.maxcash = 5;     % 现金网格上限 (归一化)
cS.mincash = 0.25;    % 现金网格下限 (归一化)
cS.save_path = 'result/'; % 结果保存路径

% --- G. 派生参数与预计算 ---
% ** CRRA 修改: 移除psi相关的派生参数 **
% cS.psi_1 = 1.0 - 1.0 / cS.psi;
% cS.psi_2 = 1.0 / cS.psi_1;
% cS.theta = (1.0 - cS.gamma) / (1.0 - cS.psi_1);

% 生存概率
cS.survprob = zeros(cS.tn - 1, 1);
survprob_data = [0.99975, 0.99974, 0.99973, 0.99972, 0.99971, 0.99969, 0.99968, 0.99966, 0.99963, 0.99961, ...
    0.99958, 0.99954, 0.99951, 0.99947, 0.99943, 0.99938, 0.99934, 0.99928, 0.99922, 0.99916, ...
    0.99908, 0.999, 0.99891, 0.99881, 0.9987, 0.99857, 0.99843, 0.99828, 0.99812, 0.99794, ...
    0.99776, 0.99756, 0.99735, 0.99713, 0.9969, 0.99666, 0.9964, 0.99613, 0.99583, 0.99551, ...
    0.99515, 0.99476, 0.99432, 0.99383, 0.9933, 0.9927, 0.99205, 0.99133, 0.99053, 0.98961, ...
    0.98852, 0.98718, 0.98553, 0.98346, 0.98089, 0.97772, 0.97391, 0.96943, 0.96429, 0.95854, ...
    0.95221, 0.94537, 0.93805, 0.93027, 0.92202, 0.91327, 0.90393, 0.89389, 0.88304, 0.87126, ...
    0.85846, 0.84452, 0.82935, 0.81282, 0.79485, 0.77543, 0.75458, 0.7324, 0.70893, 0.68424, ...
    0.68424, 0.68424];
for i = 1:min(length(survprob_data), length(cS.survprob))
    cS.survprob(i) = 1; % survprob_data(i);
end

% 离散化正态分布 (Tauchen方法)
tauchenoptions.parallel=0;
[grid, weig_matrix] = discretizeAR1_Tauchen(0, 0, 1, cS.n_shocks, 2, tauchenoptions);
cS.shock_grid = grid;
cS.shock_weig = diag(weig_matrix);

% 风险资产回报率网格
cS.gret = cS.rf + cS.mu + cS.shock_grid * cS.sigr;

% 状态转移概率 (三维冲击)
cS.nweig1 = zeros(cS.n_shocks, cS.n_shocks, cS.n_shocks);
for i6 = 1:cS.n_shocks
    for i7 = 1:cS.n_shocks
        for i8 = 1:cS.n_shocks
            cS.nweig1(i6, i7, i8) = cS.shock_weig(i6) * cS.shock_weig(i7) * cS.shock_weig(i8);
        end
    end
end

% 状态变量网格
lgcash = linspace(log(cS.mincash), log(cS.maxcash), cS.ncash)';
cS.gcash = exp(lgcash);

% --- 收入过程计算 ---
% 确定性收入剖面
cS.f_y = zeros(cS.tr - cS.tb + 1, 1);
for i1 = cS.tb:cS.tr
    age = i1;
    cS.f_y(age - cS.tb + 1) = exp(cS.aa + cS.b1*age + cS.b2*age^2 + cS.b3*age^3-3);
end

% 收入过程网格
cS.yh = zeros(cS.n_shocks, cS.n_shocks); % 暂时冲击 exp(u_t)
cS.yp = zeros(cS.n_shocks, cS.n_shocks); % 持久冲击增长 exp(z_t)
cS.gyp = zeros(cS.n_shocks, cS.n_shocks, cS.tn - 1); % 总收入增长

for i1 = 1:cS.n_shocks
    grid2_u = cS.shock_grid(i1) .* cS.corr_u_epsilon + cS.shock_grid .* sqrt(1 - cS.corr_u_epsilon^2);
    cS.yh(:, i1) = exp(grid2_u * cS.smay);

    grid2_z = cS.shock_grid(i1) .* cS.corr_z_epsilon + cS.shock_grid .* sqrt(1 - cS.corr_z_epsilon^2);
    cS.yp(:, i1) = exp(grid2_z * cS.smav);
end

% 工作期收入增长
work_periods = cS.tr - cS.tb;
for t = 1:work_periods
    G_t = cS.f_y(t+1) / cS.f_y(t); % 确定性增长
    cS.gyp(:,:,t) = repmat(G_t, cS.n_shocks, cS.n_shocks) .* cS.yp;
end

% 退休期收入增长 (无增长)
cS.gyp(:,:,(work_periods+1):(cS.tn-1)) = 1.0;

fprintf('参数设置完成。\n');
end

%
%% =====================================================================
%                       2. VFI 求解器函数 (修正版)
%  =====================================================================
function results = solve_model_vfi(cS)

% 初始化策略和价值函数数组
V = zeros(cS.ncash,  cS.tn);
C_policy = zeros(cS.ncash, cS.tn);
A_policy = zeros(cS.ncash,  cS.tn);
Q_policy = zeros(cS.ncash, cS.tn); % This is unused but kept for structure

% --- 终端期 (t=T) ---
% ** CRRA 修改: 终端价值函数 V_T = C_T^(1-gamma)/(1-gamma) **
% 其中 C_T = W_T (即 cS.gcash)
if cS.gamma == 1
    V(:, cS.tn) = log(cS.gcash);
else
    V(:, cS.tn) = cS.gcash.^(1-cS.gamma) / (1-cS.gamma);
end
C_policy(:, cS.tn) = cS.gcash;
A_policy(:, cS.tn) = 0.0;
Q_policy(:, cS.tn) = 0.0;

% --- 逆向归纳求解 ---
fprintf('开始逆向归纳求解...\n');
options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'interior-point', ...
    'MaxFunctionEvaluations', 1000, 'MaxIterations', 400, ...
    'OptimalityTolerance', 1e-10, 'StepTolerance', 1e-10, ...
    'ConstraintTolerance', 1e-10);
constant_pension_payout = 0; % No PPS in this version
% -- 退休期 (t = T-1, ..., K) --
fprintf('求解退休期...\n');
retire_start_time = tic;
for t = (cS.tn-1):-1:(cS.tr-cS.tb+1)
    fprintf('  t = %d (年龄 %d)\n', t, t + cS.tb - 1);
    V_next = V(:,t+1); % 下一期值函数

    V_next_interp_obj = griddedInterpolant(cS.gcash, V_next,  'linear', 'linear');
    local_C = zeros(cS.ncash,1);
    local_A = zeros(cS.ncash,1);
    local_V = zeros(cS.ncash,1);
    for i_cash = 1:cS.ncash
        gcash_val = cS.gcash(i_cash);

        lb = [1e-6, 0];
        x0 = [0.5 * gcash_val, 1-1e-06];
        if i_cash > 1 && local_C(i_cash-1) > 0
            x0 = [local_C(i_cash-1), local_A(i_cash-1)];
        end

        ub = [gcash_val, 1];

        % 将插值对象传递给目标函数
        obj_fun = @(x) fun_valuefunc_retired(x, gcash_val, constant_pension_payout, V_next_interp_obj, cS, t);

        [policy, fval] = fmincon(obj_fun, x0, [], [], [], [], lb, ub, [], options);

        local_C(i_cash) = policy(1);
        local_A(i_cash) = policy(2);
        local_V(i_cash) = -fval;
    end
    C_policy(:, t) = local_C;
    A_policy(:, t) = local_A;
    V(:,  t) = local_V;
end
fprintf('退休期求解完成，耗时 %.2f 分钟\n', toc(retire_start_time)/60);

% -- 工作期 (t = K-1, ..., 1) --
fprintf('求解工作期...\n');
work_start_time = tic;
for t = (cS.tr - cS.tb):-1:1
    fprintf('  t = %d (年龄 %d)\n', t, t + cS.tb - 1);
    V_next = V(:,t+1);
    V_next_interp_obj = griddedInterpolant(cS.gcash, V_next,  'spline', 'linear');
    local_C = zeros(cS.ncash,1);
    local_A = zeros(cS.ncash,1);
    local_V = zeros(cS.ncash,1);
    for i_cash = 1:cS.ncash
        gcash_val = cS.gcash(i_cash);
        lb = [1e-3, 0];
        if i_cash > 1
            x0 = [local_C(i_cash-1), local_A(i_cash-1)];
        else
            x0 = [gcash_val * 0.5, 1-1e-06]; % C, alpha
        end

        ub = [gcash_val, 1];

        % 将插值对象传递给目标函数
        obj_fun = @(x) fun_valuefunc_work(x, gcash_val, V_next_interp_obj, cS, t);

        [policy, fval] = fmincon(obj_fun, x0, [], [], [], [], lb, ub, [], options);

        local_C(i_cash) = policy(1);
        local_A(i_cash) = policy(2);
        local_V(i_cash) = -fval;
    end
    C_policy(:,  t) = local_C;
    A_policy(:,  t) = local_A;
    V(:,  t) = local_V;
end
fprintf('工作期求解完成，耗时 %.2f 分钟\n', toc(work_start_time)/60);

% --- 封装结果 ---
results.V = V;
results.C_policy = C_policy;
results.A_policy = A_policy;
end
%% =====================================================================
%                VFI - 退休期值函数 (CRRA)
%  =====================================================================
function v = fun_valuefunc_retired(x, cash, pension_payout, V_next_interp_obj, cS, t)
% x = [消费 C, 风险资产比例 alpha]
alpha = x(2);
consumption = x(1);

if consumption <= 1e-8, v = 1e10; return; end

sav = cash - consumption;
if sav < 0, v = 1e10; return; end

% ** CRRA 修改: 即期效用函数 **
if cS.gamma == 1
    u_term = log(consumption);
else
    u_term = consumption^(1 - cS.gamma) / (1 - cS.gamma);
end

% 下一期现金持有量 (归一化)
ret_fac_norm = cS.ret_fac;
cash_1 = (cS.rf * (1 - alpha) + cS.gret * alpha) * sav + ret_fac_norm + pension_payout;

% 使用 griddedInterpolant 对象进行插值
int_V = V_next_interp_obj(cash_1);

% ** CRRA 修改: 计算期望 **
% V_t = u(C_t) + beta * E[V_{t+1}]
% 退休期收入增长为1, gyp=1
expected_term = cS.shock_weig' * int_V;

% 完整的 CRRA Bellman方程
value = u_term + cS.beta * cS.survprob(t) * expected_term;

v = -value; % fmincon 最小化
end

%% =====================================================================
%                VFI - 工作期值函数 (CRRA)
%  =====================================================================

function v = fun_valuefunc_work(x, cash, V_next_interp_obj, cS, t)
% x = [绝对消费 C, 风险资产比例 alpha]
consumption = x(1);
alpha       = x(2);

if consumption <= 1e-8, v = 1e10; return; end

liquid_sav = cash - consumption;
if liquid_sav < -1e-6, v = 1e10; return; end

% ** CRRA 修改: 即期效用函数 **
if cS.gamma == 1
    u_term = log(consumption);
else
    u_term = consumption^(1 - cS.gamma) / (1 - cS.gamma);
end

% 构造下一期状态的网格 (三维冲击)
n = cS.n_shocks;
gyp_3d = repmat(reshape(cS.gyp(:,:,t), [n, 1, n]), [1, n, 1]);
gret_3d = repmat(reshape(cS.gret, [1, 1, n]), [n, n, 1]);
yh_3d = repmat(reshape(cS.yh, [1, n, n]), [n, 1, 1]);

% 下一期归一化状态变量
portfolio_return = cS.rf * (1 - alpha) + gret_3d * alpha;
cash_1 = portfolio_return .* (liquid_sav ./ gyp_3d) + (1-cS.tau_y) .* yh_3d;

% 使用 griddedInterpolant 对象进行插值
int_V = V_next_interp_obj(cash_1);

% ** CRRA 修改: 构造期望项 **
% v_t = u(c_t) + beta * E[G_{t+1}^(1-gamma) * v_{t+1}]

    term_inside_exp = (gyp_3d.^(1 - cS.gamma)) .* int_V;

expected_term = sum(cS.nweig1 .* term_inside_exp, 'all');

% 完整的 CRRA Bellman方程
value = u_term + cS.beta * cS.survprob(t) * expected_term;

v = -value; % fmincon 最小化
end


%% =====================================================================
%                       3. 模拟器函数
%  =====================================================================
function simulate_model(results, cS)
% --- 初始化 ---
fprintf('模拟开始...\n');
rng(42); % 可复现性
nsim = 10000;

% 解包策略函数
C_policy = results.C_policy;
A_policy = results.A_policy;

% 初始化模拟数组
% 冲击和外生过程
simGPY = ones(cS.tn, nsim);       % 持久收入增长因子 G_t * exp(z_t)
simY_norm = zeros(cS.tn, nsim);    % 归一化暂时收入 exp(u_t)
simR = zeros(cS.tn, nsim);         % 风险回报率

% 归一化单位下的决策和状态变量
simW_norm = zeros(cS.tn, nsim);    % 归一化流动性财富 w_t = W_t/P_t
simC_norm = zeros(cS.tn, nsim);    % 归一化消费 c_t = C_t/P_t
simA = zeros(cS.tn, nsim);         % 风险资产配置 alpha_t (已经是比例)
simS_norm = zeros(cS.tn, nsim);    % 风险资产
simB_norm = zeros(cS.tn, nsim);    % 无风险资产

% --- 1. 模拟外生冲击和收入 (采用 Antithetic Variates) ---
fprintf('  生成收入和回报率路径...\n');
work_periods = cS.tr - cS.tb;

for i1 = 1:floor(nsim/2)
    z_shocks = randn(work_periods, 1) * cS.smav;
    u_shocks_log = randn(work_periods, 1) * cS.smay;
    r_shocks = randn(cS.tn, 1) * cS.sigr;

    % 正向冲击路径
    simGPY(2:work_periods+1, i1) = cS.f_y(2:work_periods+1) ./ cS.f_y(1:work_periods) .* exp(z_shocks);
    simY_norm(1:work_periods, i1) = exp(u_shocks_log);
    simR(:, i1) = cS.rf + cS.mu + r_shocks;

    % 反向冲击路径
    i2 = nsim/2 + i1;
    simGPY(2:work_periods+1, i2) = cS.f_y(2:work_periods+1) ./ cS.f_y(1:work_periods) .* exp(-z_shocks);
    simY_norm(1:work_periods, i2) = exp(-u_shocks_log);
    simR(:, i2) = cS.rf + cS.mu - r_shocks;
end

norm_base_pension = cS.ret_fac;
simY_norm(work_periods+1:cS.tn, :) = norm_base_pension;

% --- 2. 迭代模拟生命周期决策 (在归一化单位下) ---
fprintf('  模拟生命周期决策 (归一化单位)...\n');
for t = 1:cS.tn
    if t <= work_periods % 工作期
        norm_cash = simW_norm(t, :) + (1-cS.tau_y) * simY_norm(t, :);
    else % 退休期
        norm_cash = simW_norm(t, :) + simY_norm(t, :);
    end

    C_interp_obj = griddedInterpolant(cS.gcash, C_policy(:,t),  'spline', 'linear');
    A_interp_obj = griddedInterpolant(cS.gcash, A_policy(:,t),  'spline', 'linear');
    % 插值获取最优决策
    simC_norm(t, :) = C_interp_obj(norm_cash);
    simC_prop(t,:)= simC_norm(t, :)/norm_cash;
    simA(t, :)      = A_interp_obj(norm_cash);
    
    simC_norm(t, :) = max(min(simC_norm(t, :), norm_cash), 0);
    simA(t, :)      = max(min(simA(t, :), 1), 0);

    % 计算剩余储蓄
    if t <= work_periods % 工作期
        % 安全检查
        total_outflow = simC_norm(t, :);
        scale_factor = ones(1, nsim);
        idx_exceed = total_outflow > norm_cash;
        scale_factor(idx_exceed) = norm_cash(idx_exceed) ./ total_outflow(idx_exceed) * 0.9999;
        simC_norm(t, :) = simC_norm(t, :) .* scale_factor;

        liquid_sav_norm = norm_cash - simC_norm(t, :);
    else % 退休期
        simC_norm(t,:) = min(simC_norm(t,:), norm_cash * 0.9999);
        liquid_sav_norm = norm_cash - simC_norm(t, :);
    end

    simS_norm(t, :) = simA(t, :) .* liquid_sav_norm;
    simB_norm(t, :) = liquid_sav_norm - simS_norm(t, :);

    % 更新下一期状态 (归一化单位)
    if t < cS.tn
        portfolio_return_next = simB_norm(t, :) * cS.rf + simS_norm(t, :) .* simR(t, :);
        simW_norm(t+1, :) = portfolio_return_next ./ simGPY(t+1, :);
    end
end

% --- 3. 结果汇总与绘图 (在归一化单位下) ---
fprintf('  汇总归一化结果并绘图...\n');
meanC = mean(simC_norm, 2);
meanC_prop = mean(simC_prop,2);
meanY = mean(simY_norm, 2);
meanLiquidW = mean(simW_norm, 2);
meanA = mean(simA, 2);

ages = cS.tb:cS.td;

figure('Position', [100, 100, 1200, 800]);

subplot(2,2,1);
plot(ages, meanLiquidW, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Liquid Wealth (w)');
hold on;
plot(ages, meanC, 'r-', 'LineWidth', 2, 'DisplayName', 'Consumption (c)');
plot(ages, meanY, 'k:', 'LineWidth', 2, 'DisplayName', 'Income (y)');
xline(cS.tr, '--', 'Retirement', 'LineWidth', 1.5);
hold off;
title('Mean Lifecycle Profiles (Normalized by P_t)'); xlabel('Age');
legend('show', 'Location', 'northwest'); grid on; xlim([cS.tb, cS.td]);

subplot(2,2,2);
plot(ages, meanA, 'm-', 'LineWidth', 2, 'DisplayName', 'Risky Share (\alpha)');
hold on;
plot(ages,meanC_prop,'DisplayName', 'C-prop')
xline(cS.tr, '--', 'Retirement', 'LineWidth', 1.5);
hold off;
title('Mean Portfolio Choices'); xlabel('Age'); ylabel('Share');
legend('show', 'Location', 'best'); grid on; xlim([cS.tb, cS.td]); ylim([-0.05, 1.05]);

subplot(2,2,4);
% 计算 w/y 的均值 (mean of ratios)
simWY_norm = simW_norm ./ simY_norm;
simWY_norm(isinf(simWY_norm) | isnan(simWY_norm)) = NaN;
meanWY_norm = nanmean(simWY_norm, 2);
plot(ages, meanWY_norm, 'LineWidth', 2);
xline(cS.tr, '--', 'Retirement', 'LineWidth', 1.5);
title('Mean Wealth-to-Income Ratio (w/y)'); xlabel('Age');
grid on; xlim([cS.tb, cS.td]);

% ** 修改: 更新图表标题 **
sgtitle('Lifecycle Simulation Results (Normalized, CRRA)', 'FontSize', 16, 'FontWeight', 'bold');

% ** 修改: 更新图片保存名 **
output_filename = fullfile(cS.save_path, 'vfi_benchmark_profile.png');
print(gcf, output_filename, '-dpng', '-r300');
fprintf('模拟图形已保存到 %s\n', output_filename);
end

% =====================================================================
%                辅助函数
% =====================================================================
function [z_grid, P] = discretizeAR1_Tauchen(mew, rho, sigma, znum, Tauchen_q, ~)
% Tauchen方法离散化AR(1)过程
if znum == 1
    z_grid = mew / (1 - rho);
    P = 1;
    return;
end

zstar = mew / (1 - rho);
sigmaz = sigma / sqrt(1 - rho^2);

z_grid = zstar + linspace(-Tauchen_q * sigmaz, Tauchen_q * sigmaz, znum)';
omega = z_grid(2) - z_grid(1);

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