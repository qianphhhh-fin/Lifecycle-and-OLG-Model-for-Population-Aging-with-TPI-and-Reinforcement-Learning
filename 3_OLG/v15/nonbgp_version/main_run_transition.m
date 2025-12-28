% =========================================================================
% == SCRIPT: main_run_transition.m (REVISED FOR CALIBRATION & CLARITY)
% == 目的: OLG模型过渡路径求解器 - 顶层控制脚本
% == 核心修改:
% == 1. 新增初始稳态(SS0)的校准阶段，以匹配中国初始稳态年份宏观数据。
% == 2. 将所有模拟和求解器参数集中在主脚本顶部。
% == 3. 明确了外生路径加载函数的数据流。
% =========================================================================
clear; close all;
addpath(pwd);
fprintf('=== OLG模型过渡路径求解与分析 (v3 - 带校准流程) ===\n\n');

%% --- 1. 初始化环境、模拟范围与求解器参数 ---
fprintf('--- 1. 初始化环境与参数 ---\n');

% --- 步骤 1.1: 定义模拟范围与求解器参数 ---
fprintf('   定义模拟范围与求解器参数...\n');
% 基础设定
MODEL_START_YEAR = 2023;         % 数据的起始年份（稳态年份）
TIME_STEP = 5;                   % 模型每期代表的年数
T_SIM_MAX_PERIODS = 40;          % 模拟期数 (40期 * 5年/期 = 200年)

% 求解器参数
MAX_ITER_TRANS = 100;      % 过渡路径最大迭代次数
TOL_TRANS = 1e-4;          % 资本路径的收敛容忍度
LAMBDA_TRANS = 0.2;        % 松弛因子 (Damping factor)

% [新增] 稳态求解器选择 (适用于PPS版本)
% 可选方法: 'fsolve' (快速), 'surrogateopt' (鲁棒), 'hybrid' (平衡), 'robust' (最鲁棒)
% 推荐：复杂参数/大TFP变化时使用 'hybrid' 或 'robust'
STEADY_STATE_SOLVER = 'hybrid';  % 默认使用混合求解器

% --- 步骤 1.2: 加载模型物理参数 (初始猜测) ---
fprintf('   加载模型物理参数 (作为校准起点)...\n');
cS = model_setup_utils.ParameterValues();

% --- 步骤 1.3: 将模拟设定参数添加到 cS 结构体中 ---
fprintf('   整合所有参数到 cS 结构体...\n');
cS.time_Step = TIME_STEP;
cS.ss0_year = MODEL_START_YEAR;
cS.start_year = MODEL_START_YEAR+1;
cS.end_year = cS.start_year + (T_SIM_MAX_PERIODS - 1) * cS.time_Step;
cS.T_sim = T_SIM_MAX_PERIODS;

% 网格设定 (依赖于 cS.time_Step)
ngrid = 40; 
cS.nk = ngrid; 
cS.nkpps = ngrid; 
cS.nkprime = ngrid; 
cS.npps = ngrid;
% [注释] 这里的初始网格生成是安全的，仅用于设置网格点数等基础参数
% 在稳态求解过程中，内生自适应网格系统会根据宏观猜测值动态调整网格范围
cS = model_setup_utils.generateGrids(cS);

% --- 步骤 1.4: 生成冲击过程和加载外生路径 ---
fprintf('   生成冲击过程并加载外生路径...\n');
paramS = struct();
% 生成年龄依赖冲击“信号”过程
[paramS.leGridV, paramS.TrProbM_by_age, paramS.leProb1V, cS.nw_expanded] = ...
    model_setup_utils.EarningProcess_AgeDependent(cS);
paramS.leLogGridV = log(paramS.leGridV(1:cS.nw));

% [核心修正] 现在调用 load_exogenous_paths 无需额外参数
[Z_path, A_path] = model_setup_utils.load_exogenous_paths(cS, false);
cS = model_setup_utils.calcaulte_theta_payg_path(cS, false);

% 确保所有路径长度与 T_sim 一致
Z_path = Z_path(:, 1:cS.T_sim);
A_path = A_path(1:cS.T_sim);
cS.theta_path = cS.theta_path(1:cS.T_sim);
fprintf('   已加载所有外生路径，模拟 %d 期 (%d-%d年)。\n', cS.T_sim, cS.start_year, cS.end_year);


%% --- 2. [新增] 校准初始稳态 (t=0, 对应初始稳态年份) ---
fprintf('\n--- 2. 校准初始稳态以匹配 %d 年宏观目标 ---\n', cS.ss0_year);

% --- 步骤 2.1: 设定校准目标和待校准参数 ---
TARGET_KY_RATIO = 4.3; % 目标总资本产出比 (K_total / Y), 基于文献的合理设定
fprintf('   校准目标: K/Y = %.2f\n', TARGET_KY_RATIO);

% 定义待校准参数的初始猜测值和边界
% x = [beta, gamma, lambda_g]
x0 = [cS.beta, cS.gamma, cS.lambda_g];
lb = [0.985, 0.08, 0.3]; % 下界：提高gamma下界避免求解失败
ub = [1.03, 0.12, 0.35]; % 上界：适度提高beta上界以获得更高K/Y
fprintf('   待校准参数 (初始值): beta=%.4f, gamma=%.3f, lambda_g=%.3f\n', x0(1), x0(2), x0(3));

% --- 步骤 2.2: 获取校准所需的外部数据 ---
Z_ss0 = model_setup_utils.get_calibration_inputs(cS.ss0_year, cS);
params_for_calib = struct('Z', Z_ss0, 'A', 1.0, 'theta', 0.0); % A=1.0为归一化, theta稍后计算

% 计算初始稳态年份的theta值
temp_cS = cS;
temp_cS.start_year = cS.ss0_year;
temp_cS.end_year = cS.ss0_year;
temp_cS = model_setup_utils.calcaulte_theta_payg_path(temp_cS,false);
params_for_calib.theta = temp_cS.theta_path(1);
fprintf('   已加载 %d 年人口分布, 计算得到 theta = %.4f\n', cS.ss0_year, params_for_calib.theta);

% % --- 步骤 2.3: 运行优化器进行校准 ---
% objective_fun = @(x) main_transition_utils.calibration_objective(x, cS, paramS, params_for_calib, TARGET_KY_RATIO, lb, ub);
% 
% objective_fun([1.0480, 0.148, 0.319]);
% 
% % surrogateopt 专门为计算昂贵的函数设计，非常适合稳态求解问题
% options = optimoptions('surrogateopt', ...
%     'Display', 'iter', ...
%     'MaxFunctionEvaluations', 150, ...  % 相比fminsearch大大减少评估次数
%     'PlotFcn', 'surrogateoptplot', ...                  % 设为[]关闭绘图，或'surrogateoptplot'开启
%     'UseParallel', false, ...           % 如果有并行工具箱可设为true
%     'MinSampleDistance', 1e-3, ...      % 最小采样距离，避免重复采样
%     'MinSurrogatePoints', 15);          % 初始采样点数
% 
% fprintf('   启动 surrogateopt 进行校准 (全局优化，代理模型加速)...\n');
% fprintf('   搜索范围: beta[%.3f,%.3f], gamma[%.3f,%.3f], lambda_g[%.3f,%.3f]\n', ...
%     lb(1), ub(1), lb(2), ub(2), lb(3), ub(3));
% fprintf('   预计函数评估次数: ~150次 (相比网格搜索的数千次大幅减少)\n');
% 
% [x_calibrated, fval] = surrogateopt(objective_fun, lb, ub, options);
% fprintf('--- 校准完成! ---\n');
% 
% % --- 步骤 2.4: 更新 cS 结构体为校准后的参数 ---
% cS.beta = x_calibrated(1);
% cS.gamma = x_calibrated(2);
% cS.lambda_g = x_calibrated(3);
% fprintf('✅ 校准后参数: beta=%.4f, gamma=%.3f, lambda_g=%.3f\n', cS.beta, cS.gamma, cS.lambda_g);
% fprintf('   最终模型 K/Y 与目标的平方误差: %.3e\n', fval);


%% --- 3. 求解旧稳态 (起点, 使用校准后参数) ---
fprintf('\n--- 3. 求解改革前的旧稳态 (t=0, %d年) ---\n', cS.ss0_year);
cS_old = cS; % cS中已经是校准后的参数

% 注意：旧稳态的gamma和lambda_g应该与新稳态不同，这里为了演示设为一致
% 如果改革是关于gamma或lambda_g，则需要在这里设置 cS_old.gamma = 0 等
% 这里我们假设改革是关于其他参数，或者说旧稳态就是我们校准出的状态

% 使用校准时相同的外部参数，但这次要求详细输出
% [注意] 旧稳态使用固定网格系统，保持传统稳定性
[ss_old, Dist_old, ~, ~] = ...
    main_steady_state_utils.solve_steady_state_complete(cS_old, paramS, params_for_calib, false); % verbose = true

fprintf('✅ 旧稳态求解完成: Kp=%.4f, Kg=%.4f, K/Y=%.4f\n', ...
    ss_old.K_private, ss_old.K_public, (ss_old.K_private + ss_old.K_public)/ss_old.Y_from_production);


%% --- 4. 准备并求解新稳态 (终点) ---
fprintf('\n--- 4. 求解改革后的新稳态 (t=T) ---\n');

% --- 步骤 4.1: 加载完整的过渡路径外生数据 ---
fprintf('   加载完整的外生路径数据 (%d-%d)...\n', cS.start_year, cS.end_year);
[Z_path, A_path] = model_setup_utils.load_exogenous_paths(cS, false);
cS = model_setup_utils.calcaulte_theta_payg_path(cS, false);

% 确保所有路径长度与 T_sim 一致
Z_path = Z_path(:, 1:cS.T_sim);
A_path = A_path(1:cS.T_sim);
cS.theta_path = cS.theta_path(1:cS.T_sim);
fprintf('   已加载所有外生路径，模拟 %d 期。\n', cS.T_sim);

% --- 步骤 4.2: 求解新稳态 (包含PPS) ---
cS_new = cS; % 假设新稳态使用校准后的参数，但外部环境是T期末的

% 配置PPS参数
cS_new.pps_active = true;      % 激活PPS
if cS_new.nkpps == 1  % 如果当前网格是单维的，需要扩展
    cS_new.nkpps = 20;         % 设置PPS网格点数
    cS_new.npps = 5;           % 设置PPS分配份额选择数
    % [修改] 不在这里静态生成网格，让内生自适应网格系统处理
    % 这样网格范围会根据求解过程中的宏观猜测值动态调整
    fprintf('   PPS网格将由内生自适应系统根据宏观猜测值动态生成\n');
end

params_at_end = struct('Z', Z_path(:,end), 'A', A_path(end), 'theta', cS.theta_path(end));
params_at_end.g_A_ss = (A_path(end)/A_path(end-1))^(1/cS.time_Step) - 1;

% --- [新增] 基于旧稳态和TFP增长计算智能初始猜测值 ---
TFP_growth_ratio = A_path(end) / 1.0;  % A_path(end) / A_old, 其中A_old=1.0为归一化

k_p_guess_new = ss_old.K_private * TFP_growth_ratio * 3;
k_g_guess_new = ss_old.K_public * TFP_growth_ratio * 3;

x0_new_ss = [k_p_guess_new, k_g_guess_new];
fprintf('   智能初始猜测值: Kp=%.2f, Kg=%.2f (总计K=%.2f)\n', ...
    x0_new_ss(1), x0_new_ss(2), sum(x0_new_ss));

% ===== 基于智能初始猜测值设置自适应网格 =====
% 在求解器启动前一次性设置适当的网格范围，避免求解过程中反复调整
fprintf('   🔧 基于智能初始猜测值设置自适应网格...\n');

% 1. 定义网格上限的缩放因子
GRID_SCALING_FACTOR = 10;  % 假设最富有家庭资产约为平均水平的10倍

% 2. 基于智能初始猜测值计算网格上限
k_p_guess = x0_new_ss(1);  % 私人资本初始猜测值
k_max_adaptive = GRID_SCALING_FACTOR * k_p_guess;
kpps_max_adaptive = 0.5 * k_max_adaptive;  % PPS资产网格为常规资产的一半

% 3. 应用自适应网格到cS_new
cS_new = model_setup_utils.generateGrids(cS_new, 'k_max', k_max_adaptive, 'kpps_max', kpps_max_adaptive);

fprintf('   📊 网格范围已设置: k∈[%.1f, %.1f], kpps∈[%.1f, %.1f]\n', ...
    min(cS_new.kGridV), max(cS_new.kGridV), min(cS_new.kppsGridV), max(cS_new.kppsGridV));


fprintf('   ⚙️  启动稳态求解器...\n');
[ss_new, ~, V_new, k_pol_new, cPps_pol_new] = ...
    main_steady_state_utils.solve_steady_state_complete_with_pps(cS_new, paramS, params_at_end, true, x0_new_ss, 'fsolve');

fprintf('✅ 新稳态求解完成 (包含PPS): Kp=%.4f, Kg=%.4f, K/Y=%.4f\n', ...
    ss_new.K_private, ss_new.K_public, (ss_new.K_private + ss_new.K_public)/ss_new.Y_from_production);


%% --- 5. 求解过渡路径 ---
fprintf('\n--- 5. 启动过渡路径迭代求解器 ---\n');

% [注意] 当前的过渡路径求解器尚未更新为PPS版本
% 这里我们使用旧稳态(无PPS)作为起点，新稳态(含PPS)作为终点
% 完整的过渡路径求解器PPS版本需要进一步开发

K_p_path_guess = linspace(ss_old.K_private, ss_new.K_private, cS.T_sim)';
K_g_path_guess = linspace(ss_old.K_public, ss_new.K_public, cS.T_sim)';
L_path_guess = linspace(ss_old.L, ss_new.L, cS.T_sim)';

% [临时解决方案] 使用原始的过渡路径求解器
% 注意：这可能会导致终点不完全匹配，因为新稳态包含PPS决策
fprintf('   [警告] 使用简化的过渡路径求解器，终点稳态的PPS决策可能不完全匹配\n');

TransitionResults = main_transition_utils.solve_transition_path(...
    K_p_path_guess, K_g_path_guess, L_path_guess, ...
    Dist_old, V_new, k_pol_new, ss_new, ...
    Z_path, A_path, ...
    cS, paramS, ...
    MAX_ITER_TRANS, TOL_TRANS, LAMBDA_TRANS);

if ~TransitionResults.converged
    warning('过渡路径在达到最大迭代次数后未能收敛！');
else
    fprintf('✅✅✅ 过渡路径求解成功！ ✅✅✅\n');
end



%% --- 6. 保存所有结果 ---
fprintf('\n--- 6. 保存所有结果到 .mat 文件 ---\n');
% 整合所有需要保存的信息
Output.cS = cS;
Output.cS_new = cS_new;  % 包含PPS配置的新稳态参数
Output.paramS = paramS;
Output.ss_old = ss_old;
Output.ss_new = ss_new;
Output.TransitionResults = TransitionResults;
Output.ExogenousPaths = struct('Z_path', Z_path, 'A_path', A_path, 'theta_path', cS.theta_path);

% 保存PPS相关的政策函数
if exist('cPps_pol_new', 'var')
    Output.cPps_pol_new = cPps_pol_new;
    Output.V_new = V_new;
    Output.k_pol_new = k_pol_new;
    fprintf('   已包含PPS政策函数数据\n');
end

save('transition_results.mat', 'Output', '-v7.3');
fprintf('   结果已保存到: transition_results.mat\n');

%% --- 7. (可选) 基础可视化 ---
figure('Name', '资本过渡路径');
subplot(1,2,1)
plot(1:cS.T_sim, TransitionResults.K_p_path, 'b-o', 'LineWidth', 2);
hold on;
yline(ss_old.K_private, 'k--', 'LineWidth', 1.5);
yline(ss_new.K_private, 'r--', 'LineWidth', 1.5);
title('私人资本 (K_p) 过渡路径');
xlabel('时期 (t)'); ylabel('资本存量');
legend('过渡路径', '旧稳态', '新稳态 (含PPS)', 'Location', 'best');
grid on;

subplot(1,2,2)
plot(1:cS.T_sim, TransitionResults.K_g_path, 'b-s', 'LineWidth', 2);
hold on;
yline(ss_old.K_public, 'k--', 'LineWidth', 1.5);
yline(ss_new.K_public, 'r--', 'LineWidth', 1.5);
title('公共资本 (K_g) 过渡路径');
xlabel('时期 (t)'); ylabel('资本存量');
legend('过渡路径', '旧稳态', '新稳态 (含PPS)', 'Location', 'best');
grid on;

fprintf('\n--- 流程结束 ---\n');