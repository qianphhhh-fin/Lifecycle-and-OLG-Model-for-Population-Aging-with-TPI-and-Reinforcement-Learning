% =========================================================================
% == SCRIPT: 比较parfor优化版
% =========================================================================
clear; close all;
addpath(pwd);
fprintf('=== OLG模型过渡路径求解与分析 (BGP版本 - 稳态化模型) ===\n\n');

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

% 稳态求解器选择
STEADY_STATE_SOLVER = 'fsolve';  % 默认使用混合求解器

% --- 步骤 1.2: 加载模型物理参数 (初始猜测) ---
fprintf('   加载模型物理参数 (作为校准起点)...\n');
cS = model_setup_utils_bgp.ParameterValues();

% --- 步骤 1.3: 将模拟设定参数添加到 cS 结构体中 ---
fprintf('   整合所有参数到 cS 结构体...\n');
cS.time_Step = TIME_STEP;
cS.ss0_year = MODEL_START_YEAR;
cS.start_year = MODEL_START_YEAR+1;
cS.end_year = cS.start_year + (T_SIM_MAX_PERIODS - 1) * cS.time_Step;
cS.T_sim = T_SIM_MAX_PERIODS;

% 网格设定
ngrid = 40; 
cS.nk = ngrid; 
cS.nkpps = ngrid; 
cS.nkprime = ngrid; 
cS.npps = ngrid;
cS = model_setup_utils_bgp.generateGrids(cS);

% --- 步骤 1.4: 生成冲击过程和加载外生路径 ---
fprintf('   生成冲击过程...\n');
paramS = struct();
% 生成年龄依赖冲击"信号"过程
[paramS.leGridV, paramS.TrProbM_by_age, paramS.leProb1V, cS.nw_expanded] = ...
    model_setup_utils_bgp.EarningProcess_AgeDependent(cS);
paramS.leLogGridV = log(paramS.leGridV(1:cS.nw));

fprintf('   加载过渡路径所需的外生数据...\n');

% [核心修正] 调用重构后的 load_exogenous_paths 函数
% 明确指定 mode = 'transition' 来获取完整的 Z_path 和 A_path
% 第三个参数 plot_flag 设为 false，因为这里只是加载数据，不需要立即绘图
[Z_path, A_path] = model_setup_utils_bgp.load_exogenous_paths(cS);

% [逻辑修正] PAYG税率路径的计算应该在加载Z_path之后，因为它可能依赖于人口结构
cS = model_setup_utils_bgp.calcaulte_theta_payg_path(cS, false);

% 确保所有路径长度与 T_sim 一致
% 注意：load_exogenous_paths 内部已经保证了路径长度正确，但为保险起见可以保留
Z_path = Z_path(:, 1:cS.T_sim);
A_path = A_path(1:cS.T_sim);
cS.theta_path = cS.theta_path(1:cS.T_sim);
fprintf('   已加载所有外生路径，模拟 %d 期 (%d-%d年)。\n', cS.T_sim, cS.start_year, cS.end_year);

%% --- 2. [BGP修改] 校准与求解初始稳态 (t=0, 对应初始稳态年份) ---
fprintf('\n--- 2. 校准与求解初始“伪稳态”以匹配 %d 年宏观目标 ---\n', cS.ss0_year);

% --- 步骤 2.1: 设定校准目标和待校准参数 ---
TARGET_KY_RATIO = 4.3; % 目标总资本产出比 (K_total / Y)
fprintf('   校准目标: K/Y = %.2f\n', TARGET_KY_RATIO);

% 为初始稳态设定经济增长率（通常校准时设为0或一个较低的近期增长率）
cS.g_A_ss = 0.0; 
fprintf('   已设定初始稳态的长期技术年增长率: g_A_ss = %.3f\n', cS.g_A_ss);

% 定义待校准参数的初始猜测值和边界
x0_calib = [cS.beta, cS.gamma, cS.lambda_g];
lb_calib = [0.985, 0.08, 0.3];
ub_calib = [1.03, 0.12, 0.35];
fprintf('   待校准参数 (初始值): beta=%.4f, gamma=%.3f, lambda_g=%.3f\n', x0_calib(1), x0_calib(2), x0_calib(3));

% --- 步骤 2.2: 获取校准基准年(2023年)的外部数据 ---
fprintf('   正在为 %d 年加载外部数据...\n', cS.ss0_year);

% [核心] 获取2023年的现实人口分布，这将作为初始“伪稳态”的人口结构基础
Z_ss0 = model_setup_utils_bgp.get_calibration_inputs(cS.ss0_year, cS);

% 计算2023年对应的PAYG缴费率theta
temp_cS = cS;
temp_cS.start_year = cS.ss0_year;
temp_cS.end_year = cS.ss0_year;
temp_cS = model_setup_utils_bgp.calcaulte_theta_payg_path(temp_cS, false);
theta_ss0 = temp_cS.theta_path(1);
fprintf('   已加载 %d 年人口分布, 计算得到 theta = %.4f\n', cS.ss0_year, theta_ss0);

% [核心] 构建一个专门用于初始稳态求解的参数结构体
% 这个结构体将强制求解器使用2023年的现实人口数据 Z_ss0
params_for_ss0 = struct(...
    'Z', Z_ss0, ...            % <<-- 关键：传入2023年的人口分布
    'A', 1.0, ...              % 标准化初始技术水平为1
    'theta', theta_ss0, ...    % 2023年的theta
    'g_A_ss', cS.g_A_ss ...    % 初始稳态的技术增长率
);

% --- 步骤 2.3: 运行优化器进行校准 ---
% [此部分为演示目的暂时注释，实际使用时可取消注释]
% fprintf('   启动校准流程...\n');
% % 校准时，目标函数内部会调用稳态求解器，并传入 params_for_ss0
% objective_fun = @(x) main_steady_state_utils_bgp.calibration_objective(x, cS, paramS, params_for_ss0, TARGET_KY_RATIO, lb_calib, ub_calib);
% [~, calibrated_params] = fmincon(objective_fun, x0_calib, ...); 
% % 更新cS中的参数为校准后的结果
% cS.beta = calibrated_params(1); 
% cS.gamma = calibrated_params(2);
% cS.lambda_g = calibrated_params(3);
% fprintf('✅ 校准完成。\n');