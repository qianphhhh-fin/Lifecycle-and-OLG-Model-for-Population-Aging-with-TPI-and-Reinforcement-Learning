% =========================================================================
% == SCRIPT: main_run_transition_bgp.m
% == 版本: [v2.0 - 比较静态分析版]
% ==
% == 目的:
% ==   通过求解两个关键的稳态，为完整的转轨动态分析提供起点和终点。
% ==
% == 核心流程:
% ==   1. [路径生成] 调用内生路径生成器，确定总模拟期数(T_sim)并获取
% ==      人口(Z_path)、TFP(A_path)和养老金缴费率(theta_path)的完整路径。
% ==   2. [求解初始稳态 ss0]
% ==      - 特征: 2023年人口结构, 无长期增长(g=0), 无PPS制度。
% ==      - 作用: 代表经济在改革前的、现实的“伪稳态”起点。
% ==   3. [求解终期稳态 ssF]
% ==      - 特征: 内生决定的长期稳态人口结构, 有长期增长(g>0), 有PPS制度。
% ==      - 作用: 代表经济在所有冲击完成后，最终会收敛到的BGP稳态终点。
% =========================================================================
clear; close all; clear classes;
addpath(pwd);
fprintf('=== OLG模型 BGP比较静态分析脚本 (求解初始与终期稳态) ===\n\n');

%% --- 1. 全局设置与初始化 ---
fprintf('--- 1. 全局设置与初始化 ---\n');

TIME_STEP = 5;
ngrid = 10;       % 常规资产网格密度
ngrid_pps = 10;   % PPS资产网格密度

% 加载模型物理参数
fprintf('   加载模型物理参数...\n');
cS = model_setup_utils_bgp.ParameterValues();
cS.time_Step = TIME_STEP;

% 设定初始年份和模拟起点
cS.ss0_year = 2023;
cS.start_year = 2023;

%% --- 2. [核心] 生成完整的外生路径 ---
fprintf('\n--- 2. 生成完整的外生路径 (内生决定模拟长度) ---\n');

% 步骤 2.1: 调用内生路径生成器获取人口和TFP路径
% 注意：此函数现在会返回一个更新后的cS，其中包含了内生计算出的T_sim和end_year
[Z_path, A_path, cS] = model_setup_utils_bgp.generate_exo_paths(cS, true);

% 步骤 2.2: 调用函数获取PAYG缴费率路径
% 注意：calcaulte_theta_payg_path 需要 cS 中有 start_year 和 end_year
cS = model_setup_utils_bgp.calcaulte_theta_payg_path(cS, false); % false表示不画图
theta_path = cS.theta_path;


%% --- 3. 求解初始稳态 (ss0, t=0) ---
% fprintf('\n\n--- 3. 求解初始伪稳态 (ss0, t=0) ---\n');
% fprintf('   特征: %d年人口结构, 无长期增长(g=0), 无PPS制度\n', cS.ss0_year);
% 
% % --- 步骤 3.1: 为ss0设定特定的经济环境 ---
% cS0 = cS; % 创建一个独立的副本以防污染
% cS0.pps_active = false;
% cS0.g_A_ss = 0.0; % 初始稳态无增长
% 
% % 设定网格和收入过程 (无PPS)
% cS0.nk = ngrid; cS0.nkpps = 1; cS0.nkprime = ngrid; cS0.npps = 1;
% cS0 = model_setup_utils_bgp.generateGrids(cS0);
% paramS0 = struct();
% [paramS0.leGridV, paramS0.TrProbM_by_age, paramS0.leProb1V, cS0.nw_expanded] = ...
%     model_setup_utils_bgp.EarningProcess_AgeDependent(cS0);
% paramS0.leLogGridV = log(paramS0.leGridV(1:cS0.nw));
% 
% % 打包求解ss0所需的外生参数
% params_for_ss0 = struct(...
%     'Z', model_setup_utils_bgp.get_calibration_inputs(cS.ss0_year, cS0), ... % 2023年人口分布
%     'A', 1.0, ...
%     'theta', theta_path(1), ... % 使用路径的期初值
%     'g_A_ss', cS0.g_A_ss);     % g=0
% 
% % --- 步骤 3.2: 调用统一稳态求解器 ---
% tic;
% fprintf('   ⚙️  启动统一稳态求解器 (求解 ss0)...\n');
% [ss0, ~, ~, ~] = ...
%     main_steady_state_utils_bgp.solve_steady_state_unified(cS0, paramS0, params_for_ss0, true);
% toc;
% 
% if isempty(ss0)
%     error('初始稳态(ss0)求解失败！请检查代码逻辑或参数。');
% else
%     fprintf('✅ 初始稳态(ss0)求解成功！\n');
%     fprintf('   利率和工资: r=%.4f, w=%.4f\n', ss0.r_mkt, ss0.w_hat);
% end


%% --- 4. 求解终期稳态 (ssF, t=T) ---
fprintf('\n\n--- 4. 求解终期BGP稳态 (ssF, t=T) ---\n');
fprintf('   特征: 长期稳态人口结构, 有长期增长(g=%.4f), 有PPS制度\n', cS.g_A_ss);

% --- 步骤 4.1: 为ssF设定特定的经济环境 ---
cSF = cS; % 创建一个独立的副本
cSF.pps_active = true; % 激活PPS
cSF.g_A_ss = cS.g_A_ss; % 直接使用中心定义的参数

% 设定PPS相关制度参数 (这些可以放在ParameterValues中，或在这里明确设定)
cSF.pps_contrib_rate = 0.03;
cSF.pps_withdrawal_rate = 0.10;
cSF.pps_tax_rate_withdrawal = 0.03;
fprintf('   设定PPS参数: 缴费率=%.2f%%, 提取率=%.2f%%, 提取税率=%.2f%%\n', ...
    cSF.pps_contrib_rate*100, cSF.pps_withdrawal_rate*100, cSF.pps_tax_rate_withdrawal*100);

% 设定网格和收入过程 (有PPS)
cSF.nk = ngrid; cSF.nkpps = ngrid_pps; cSF.nkprime = ngrid; cSF.npps = ngrid_pps;
cSF = model_setup_utils_bgp.generateGrids(cSF);
paramSF = struct();
[paramSF.leGridV, paramSF.TrProbM_by_age, paramSF.leProb1V, cSF.nw_expanded] = ...
    model_setup_utils_bgp.EarningProcess_AgeDependent(cSF);
paramSF.leLogGridV = log(paramSF.leGridV(1:cSF.nw));

% 打包求解ssF所需的外生参数
params_for_ssF = struct(...
    'Z', Z_path(:, end), ...   % 使用路径的期末人口分布 (稳态分布)
    'A', 1.0, ...
    'theta', theta_path(end), ... % 使用路径的期末值
    'g_A_ss', cSF.g_A_ss);       % g > 0

% --- 步骤 4.2: 调用统一稳态求解器 ---
tic;
fprintf('   ⚙️  启动统一稳态求解器 (求解 ssF)...\n');
[ssF, ~, ~, ~] = ...
    main_steady_state_utils_bgp.solve_steady_state_unified(cSF, paramSF, params_for_ssF, true);
toc;

if isempty(ssF)
    error('终期稳态(ssF)求解失败！请检查代码逻辑或参数。');
else
    fprintf('✅ 终期稳态(ssF)求解成功！\n');
    fprintf('   利率和工资: r=%.4f, w=%.4f\n', ssF.r_mkt, ssF.w_hat);
end


%% --- 5. 结果对比与总结 ---
fprintf('\n\n===========================================================\n');
fprintf('===           比较静态分析: 关键宏观变量对比          ===\n');
fprintf('===========================================================\n');
fprintf('%-25s | %-15s | %-15s \n', '宏观变量', '初始稳态 (ss0)', '终期稳态 (ssF)');
fprintf('-----------------------------------------------------------\n');
% 为BGP下的资本做增长调整，使其与非增长稳态可比
% K_hat_F_comparable = ssF.K_private_hat / (1 + cSF.g_A_ss)^cSF.time_Step;
% K_pps_hat_F_comparable = ssF.K_pps_hat / (1 + cSF.g_A_ss)^cSF.time_Step;
% Y_hat_F_comparable = ssF.Y_from_production_hat / (1 + cSF.g_A_ss)^cSF.time_Step;
% 上述方法不完全对，因为利率和工资不同，直接比较标准化值更有意义。
fprintf('%-25s | %-15.4f | %-15.4f \n', '私人资本 (K_p_hat)', ss0.K_private_hat, ssF.K_private_hat);
fprintf('%-25s | %-15.4f | %-15.4f \n', 'PPS资本 (K_pps_hat)', ss0.K_pps_hat, ssF.K_pps_hat);
fprintf('%-25s | %-15.4f | %-15.4f \n', '公共资本 (K_g_hat)', ss0.K_public_hat, ssF.K_public_hat);
fprintf('%-25s | %-15.4f | %-15.4f \n', '总劳动 (L_hat)', ss0.L_hat, ssF.L_hat);
fprintf('%-25s | %-15.4f | %-15.4f \n', '总产出 (Y_hat)', ss0.Y_from_production_hat, ssF.Y_from_production_hat);
fprintf('-----------------------------------------------------------\n');
fprintf('%-25s | %-15.4f | %-15.4f \n', '市场利率 (r)', ss0.r_mkt, ssF.r_mkt);
fprintf('%-25s | %-15.4f | %-15.4f \n', '标准化工资 (w_hat)', ss0.w_hat, ssF.w_hat);
fprintf('===========================================================\n');

fprintf('\n--- 比较静态分析脚本执行完毕 ---\n');