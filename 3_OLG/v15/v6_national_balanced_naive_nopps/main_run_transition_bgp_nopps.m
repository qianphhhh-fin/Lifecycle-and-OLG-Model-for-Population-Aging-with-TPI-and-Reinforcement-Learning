% =========================================================================
% == SCRIPT: main_run_steady_state_test_GOLDEN.m
% == 目的: [黄金标准版] OLG模型稳态求解器 - 独立测试脚本
% == 核心:
% ==   1. [独立性] 此脚本只专注于求解一个独立的、理论上自洽的BGP稳态，
% ==      以验证核心求解器和宏观账户的正确性。
% ==   2. [清晰性] 移除了所有与过渡路径相关的复杂设定，避免数据污染。
% ==   3. [可靠性] 这是验证模型能否正确求解的最终标准。
% =========================================================================
clear; close all;
addpath(pwd);
fprintf('=== OLG模型 BGP稳态求解器 黄金标准测试 ===\n\n');

%% --- 1. 初始化环境与参数 ---
fprintf('--- 1. 初始化环境与参数 ---\n');

% --- 步骤 1.1: 定义模拟和模型参数 ---
TIME_STEP = 5;
ngrid = 100;

% --- 步骤 1.2: 加载模型物理参数 ---
fprintf('   加载模型物理参数...\n');
cS = model_setup_utils_bgp.ParameterValues();

% --- 步骤 1.3: 设定BGP稳态的特定参数 ---
fprintf('   设定BGP稳态的经济参数...\n');
cS.time_Step = TIME_STEP;
cS.g_A_ss = 0.015; % 设定一个正的长期技术年增长率
cS.pps_active = false; % 设定为无PPS模式

% --- 步骤 1.4: 生成网格和收入过程 ---
cS.nk = ngrid; cS.nkpps = 1; cS.nkprime = ngrid; cS.npps = 1;
cS = model_setup_utils_bgp.generateGrids(cS);

fprintf('   生成收入过程...\n');
paramS = struct();
[paramS.leGridV, paramS.TrProbM_by_age, paramS.leProb1V, cS.nw_expanded] = ...
    model_setup_utils_bgp.EarningProcess_AgeDependent(cS);
paramS.leLogGridV = log(paramS.leGridV(1:cS.nw));

%% --- 2. [核心] 求解理论上的BGP稳态 ---
fprintf('\n--- 2. 求解一个独立的BGP理论稳态 (t=T) ---\n');

% --- 步骤 2.1: 为理论稳态设定外生环境 ---
fprintf('   计算理想化的ZPG理论稳态人口分布...\n');
Z_steady_state = model_setup_utils_bgp.compute_theoretical_ss_dist_zpg(cS);

% 终期稳态的theta值需要基于其自身的人口结构来计算
% 为了简化，我们暂时使用一个固定的、合理的值
theta_steady_state = 0.08;
fprintf('   设定终期稳态的养老金缴费率 theta = %.4f\n', theta_steady_state);

% 将所有外部参数打包
params_for_ss = struct(...
    'Z', Z_steady_state, ...
    'A', 1.0, ...
    'theta', theta_steady_state, ...
    'g_A_ss', cS.g_A_ss);

% --- 步骤 2.2: 调用统一稳态求解器 ---
tic;
fprintf('   ⚙️  启动统一稳态求解器 (BGP模式)...\n');
[ss, Dist, polS, valS] = ...
    main_steady_state_utils_bgp.solve_steady_state_unified(cS, paramS, params_for_ss, true);
toc;

%% --- 3. 结果分析 ---
if isempty(ss)
    error('BGP稳态求解失败！请检查代码逻辑或参数。');
else
    fprintf('✅ BGP稳态求解成功！\n');
    fprintf('\n--- 最终结果概览 ---\n');
    fprintf('   标准化结果: K̂p=%.4f, K̂g=%.4f, L̂=%.4f, Ŷ=%.4f\n', ...
        ss.K_private_hat/(1+(1+cS.g_A_ss)^cS.time_Step-1), ss.K_public_hat, ss.L_hat, ss.Y_from_production_hat);
    fprintf('   利率和工资: r=%.4f, w=%.4f\n', ss.r_mkt, ss.w_hat);
end

% 无需保存或画图，只专注于验证。
fprintf('\n--- BGP稳态求解流程结束 (黄金标准测试) ---\n');