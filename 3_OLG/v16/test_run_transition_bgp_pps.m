% --- START OF FILE main_run_transition_bgp_pps.m ---

% =========================================================================
% == SCRIPT: main_run_steady_state_test_PPS.m
% == 目的: [带PPS功能版] OLG模型稳态求解器 - 独立测试脚本
% == 核心:
% ==   1. [新功能] 引入并测试一个带特定规则的私人养老金(PPS)体系。
% ==   2. [可切换] 通过 cS.pps_active 开关控制是否激活PPS模式。
% ==   3. [兼容性] 当 cS.pps_active = false 时，结果应与原版完全一致。
% =========================================================================
clear; close all; clear classes;
addpath(pwd);
fprintf('=== OLG模型 BGP稳态求解器 (含PPS功能测试) ===\n\n');

%% --- 1. 初始化环境与参数 ---
fprintf('--- 1. 初始化环境与参数 ---\n');

% --- 步骤 1.1: 定义模拟和模型参数 ---
TIME_STEP = 5;
ngrid = 40;
ngrid_pps = 20; % 为PPS资产设定网格密度

% --- 步骤 1.2: 加载模型物理参数 ---
fprintf('   加载模型物理参数...\n');
cS = model_setup_utils_bgp.ParameterValues();

% --- 步骤 1.3: 设定BGP稳态的特定参数 ---
fprintf('   设定BGP稳态的经济参数...\n');
cS.time_Step = TIME_STEP;
cS.g_A_ss = 0.012; % 设定一个正的长期技术年增长率

% --- [核心开关] 激活PPS模式 ---
cS.pps_active = false;
fprintf('   [模式设定] PPS模式已激活: %s\n', mat2str(cS.pps_active));

if cS.pps_active
    % --- 设定PPS相关制度参数 ---
    cS.pps_contrib_rate = 0.03;         % PPS缴费率 (占劳动收入)
    cS.pps_withdrawal_rate = 0.10;      % 退休后每期强制提取比例
    cS.pps_tax_rate_withdrawal = 0.03;  % PPS提取时的税率
    fprintf('   设定PPS参数: 缴费率=%.2f%%, 提取率=%.2f%%, 提取税率=%.2f%%\n', ...
        cS.pps_contrib_rate*100, cS.pps_withdrawal_rate*100, cS.pps_tax_rate_withdrawal*100);
end

% --- [新增核心修改] 设定养老金模式 ---
cS.endogenous_theta_mode = true; % <--- 新增开关，true表示DB模式，false表示DC模式
if cS.endogenous_theta_mode
    cS.endogenous_theta_mode = true;
    cS.payg_replacement_rate = 0.4;
    cS.theta_max = 0.99; % <--- 在这里设置你的缴费率上限！
    fprintf('   [养老金模式] 内生缴费率(DB)模式激活, 目标替代率: %.2f%%, 缴费率上限: %.2f%%\n', ...
        cS.payg_replacement_rate * 100, cS.theta_max * 100);
else
    % 在DC模式下，使用ParameterValues()中设定的外生缴费率
    fprintf('   [养老金模式] 外生缴费率(DC)模式已激活, 缴费率 theta = %.4f\n', cS.theta_path);
end


% --- 步骤 1.4: 生成网格和收入过程 ---
if cS.pps_active
    cS.nk = ngrid; cS.nkpps = ngrid_pps; cS.nkprime = ngrid; cS.npps = ngrid_pps;
else
    cS.nk = ngrid; cS.nkpps = 1; cS.nkprime = ngrid; cS.npps = 1;
end
cS = model_setup_utils_bgp.generateGrids(cS);

fprintf('   生成收入过程...\n');
paramS = struct();
[paramS.leGridV, paramS.TrProbM_by_age, paramS.leProb1V, cS.nw_expanded] = ...
    model_setup_utils_bgp.EarningProcess_AgeDependent(cS);
paramS.leLogGridV = log(paramS.leGridV(1:cS.nw));

%% --- 2. [核心] 求解理论上的BGP稳态 ---
fprintf('\n--- 2. 求解一个独立的BGP理论稳态 (t=T) ---\n');

% --- 步骤 2.1: 为理论稳态设定外生环境 ---
cS.n_ss = -0.0184;     % [新增] 设定长期人口年增长率 (0.0 表示ZPG)
cS.s_pathV = cS.s_pathV(:,end);
fprintf('   设定长期增长率: 技术 g_A=%.3f, 人口 n=%.3f\n', cS.g_A_ss, cS.n_ss);
Z_steady_state = model_setup_utils_bgp.compute_theoretical_ss_dist(cS, cS.n_ss, cS.s_pathV);

% --- [逻辑修改] theta现在可能是内生的，所以不再需要在这里设定一个固定的theta值 ---
if cS.endogenous_theta_mode
    % 在DB模式下，我们传递一个占位符，因为它将在求解器内部被完全忽略和重新计算
    theta_for_solver = [];
    fprintf('   养老金缴费率 theta 将在求解器内部内生决定。\n');
else
    % 在DC模式下，我们使用一个预设的合理值
    theta_steady_state = 0.700; % 使用参数文件中的值
    fprintf('   设定终期稳态的养老金缴费率 theta = %.4f\n', theta_steady_state);
    theta_for_solver = theta_steady_state;
end

% 将所有外部参数打包
params_for_ss = struct(...
    'Z', Z_steady_state, ...
    'A', 1.0, ...
    'theta', theta_for_solver, ... % 传递占位符或实际值
    'g_A_ss', cS.g_A_ss,...
    'n_ss', cS.n_ss);

% --- 步骤 2.2: 调用统一稳态求解器 ---
tic;
fprintf('   ⚙️  启动统一稳态求解器 (BGP模式)...\n');
[ss, Dist, polS, valS] = ...
    main_steady_state_utils_bgp.solve_steady_state_unified(cS, paramS, params_for_ss, true, [], 'lsqnonlin');
toc;

%% --- 3. 结果分析 ---
if isempty(ss)
    error('BGP稳态求解失败！请检查代码逻辑或参数。');
else
    fprintf('✅ BGP稳态求解成功！\n');
    fprintf('\n--- 最终结果概览 ---\n');
    fprintf('   标准化结果: K̂p=%.4f, K̂g=%.4f, L̂=%.4f, Ŷ=%.4f\n', ...
        ss.K_private_hat, ss.K_public_hat, ss.L_hat, ss.Y_from_production_hat);
    if cS.pps_active
        % 原始代码
        fprintf('   PPS相关: K̂_pps=%.4f\n', ss.K_pps_hat/(1+(1+cS.g_A_ss)^cS.time_Step-1));
    end
    fprintf('   利率和工资: r=%.4f, w=%.4f\n', ss.r_mkt, ss.w_hat);
    % [新增] 报告内生计算出的theta
    if cS.endogenous_theta_mode
        fprintf('   内生养老金缴费率: theta=%.4f (%.2f%%)\n', ss.theta, ss.theta*100);
    end
end

% 无需保存或画图，只专注于验证。
fprintf('\n--- BGP稳态求解流程结束 (PPS功能测试) ---\n');