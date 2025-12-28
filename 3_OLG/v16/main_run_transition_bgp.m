% =========================================================================
% == SCRIPT: main_run_transition_bgp.m
% == 版本: [v3.0 - 最终整合与模块化版]
% ==
% == 目的:
% ==   通过求解两个关键的稳态，为完整的转轨动态分析提供起点和终点。
% ==   此版本已完全兼容最新的utils功能，包括内生/外生theta和PPS模块。
% ==
% == 核心流程:
% ==   1. [全局设定] 在脚本顶部设定核心的政策开关。
% ==   2. [路径生成] 获取人口、TFP和养老金缴费率的完整路径。
% ==   3. [求解初始稳态 ss0] (try...end 模块)
% ==      - 强制设定为改革前环境 (无PPS, 可选外生theta)。
% ==   4. [求解终期稳态 ssF] (try...end 模块)
% ==      - 根据全局设定决定其环境 (可激活PPS, 可内生theta)。
% ==   5. [结果对比] 提供清晰的比较静态分析表格。
% =========================================================================
clear; close all; clear classes;
addpath(pwd);
fprintf('=== OLG模型 BGP比较静态分析脚本 (求解初始与终期稳态) ===\n\n');

%% --- 1. 全局设置与初始化 ---
fprintf('--- 1. 全局设置与初始化 ---\n');

ngrid = 50;       % 常规资产网格密度
ngrid_pps = 20;   % PPS资产网格密度

% 加载模型物理参数
fprintf('   加载模型物理参数...\n');
cS = model_setup_utils_bgp.ParameterValues();


% --- [核心] 设定全局政策情景开关 ---
cS.endogenous_theta_mode = true;  % true: DB模式(内生theta); false: DC模式(外生theta)
cS.pps_active = true;             % true: 激活PPS制度; false: 关闭PPS制度

% 设定初始年份和模拟起点
cS.ss0_year = 2023;
cS.start_year = 2023;

% --- 如果是DB模式，在此处设定全局的政策参数 ---
if cS.endogenous_theta_mode
    cS.payg_replacement_rate = 0.35; % 目标替代率
    cS.theta_max = 0.70;             % 缴费率上限
end


%% --- 2. [核心] 生成完整的外生路径 ---
fprintf('\n--- 2. 生成完整的外生路径 (内生决定模拟长度) ---\n');

% 步骤 2.1: 调用内生路径生成器获取人口和TFP路径
% 注意：此函数现在会返回一个更新后的cS，其中包含了内生计算出的T_sim和end_year
[Z_path, A_path, cS] = model_setup_utils_bgp.generate_exo_paths(cS, false); % false: 不画人口图

% 步骤 2.2: 调用函数获取PAYG缴费率路径 (在DC模式下会用到)
cS = model_setup_utils_bgp.calcaulte_theta_payg_path(cS, false);
theta_path = cS.theta_path;

ss0 = struct(); % 初始化以防求解失败
ssF = struct(); % 初始化以防求解失败

%% --- 3. 求解初始稳态 (ss0, t=0) ---
fprintf('\n\n--- 3. 求解初始伪稳态 (ss0, t=0) ---\n');
fprintf('   特征: %d年人口结构, 无长期增长(g=0), 无PPS制度\n', cS.ss0_year);

% --- 步骤 3.1: 为ss0设定特定的经济环境 ---
cS0 = cS; % 创建一个独立的副本以防污染
cS0.pps_active = false; % **强制**关闭PPS
cS0.g_A_ss = 0.0;     % 初始稳态无增长
cS0.n_ss = 0.0;       % 初始稳态人口零增长（注意：这只是为了计算方便，实际人口结构是给定的）
cS0.s_pathV = cS0.s_pathV(:,1); % 初始稳态使用初期存活率

% [养老金模式适配]
if cS0.endogenous_theta_mode
    fprintf('   [养老金模式] DB模式激活, 目标替代率=%.2f%%, 上限=%.2f%%\n', ...
        cS0.payg_replacement_rate*100, cS0.theta_max*100);
    theta_for_ss0 = []; % 在DB模式下，传递空值
else
    fprintf('   [养老金模式] DC模式激活, 外生缴费率=%.4f\n', theta_path(1));
    theta_for_ss0 = theta_path(1); % 使用路径的期初值
end

% 设定网格和收入过程 (无PPS)
cS0.nk = ngrid; cS0.nkpps = 1; cS0.nkprime = ngrid; cS0.npps = 1;
cS0 = model_setup_utils_bgp.generateGrids(cS0);
paramS0 = struct();
[paramS0.leGridV, paramS0.TrProbM_by_age, paramS0.leProb1V, cS0.nw_expanded] = ...
    model_setup_utils_bgp.EarningProcess_AgeDependent(cS0);
paramS0.leLogGridV = log(paramS0.leGridV(1:cS0.nw));

% 打包求解ss0所需的外生参数
params_for_ss0 = struct(...
    'Z', Z_path(:,1), ... % 2023年人口分布
    'A', 1.0, ...
    'theta', theta_for_ss0, ...
    'g_A_ss', cS0.g_A_ss,...
    'n_ss', cS0.n_ss);

% --- 步骤 3.2: 调用统一稳态求解器 ---
tic;
fprintf('   ⚙️  启动统一稳态求解器 (求解 ss0)...\n');
[ss0_result, dist0, pol0, val0] = ...
    main_steady_state_utils_bgp.solve_steady_state_unified(cS0, paramS0, params_for_ss0, true, [], 'lsqnonlin');
toc;

if isempty(ss0_result)
    warning('初始稳态(ss0)求解失败！');
else
    ss0 = ss0_result;
    fprintf('✅ 初始稳态(ss0)求解成功！\n');
    fprintf('   利率和工资: r=%.4f, w=%.4f\n', ss0.r_mkt, ss0.w_hat);

    % --- [新增代码开始] ---
    % 增加一个检查以防Y为零或负数导致计算错误
    if ss0.Y_from_production_hat > 1e-9
        % 1. 计算总资本 (包括私人、公共和PPS)
        total_capital_ss0 = ss0.K_private_hat + ss0.K_public_hat + ss0.K_pps_hat;
        ky_total_ratio_ss0 = total_capital_ss0 / ss0.Y_from_production_hat;
        
        % 2. 计算私人资本 (仅包括家庭和PPS账户)
        private_capital_ss0 = ss0.K_private_hat + ss0.K_pps_hat;
        ky_private_ratio_ss0 = private_capital_ss0 / ss0.Y_from_production_hat;
        
        % 3. 打印两个比率
        fprintf('   总资本产出比 (Total K/Y) .... = %.4f\n', ky_total_ratio_ss0);
        fprintf('   私人资本产出比 (Private Kp/Y) = %.4f\n', ky_private_ratio_ss0);
    else
        fprintf('   资本产出比: 无法计算 (Y_hat <= 0)\n');
    end
    % --- [新增代码结束] ---
end


%% --- 4. 求解终期稳态 (ssF, t=T) ---

fprintf('\n\n--- 4. 求解终期BGP稳态 (ssF, t=T) ---\n');

% --- 步骤 4.1: 为ssF设定特定的经济环境 ---
cSF = cS; % 创建一个独立的副本
cSF.pps_active = cS.pps_active; % **继承**全局PPS开关状态
cSF.s_pathV = cSF.s_pathV(:,end); % 初始稳态使用初期存活率

fprintf('   特征: 长期稳态人口(n=%.3f), 长期增长(g=%.3f), PPS激活=%s\n', ...
    cSF.n_ss, cSF.g_A_ss, mat2str(cSF.pps_active));

% [养老金模式适配]
if cSF.endogenous_theta_mode
    fprintf('   [养老金模式] DB模式激活, 目标替代率=%.2f%%, 上限=%.2f%%\n', ...
        cSF.payg_replacement_rate*100, cSF.theta_max*100);
    theta_for_ssF = [];
else
    fprintf('   [养老金模式] DC模式激活, 外生缴费率=%.4f\n', theta_path(end));
    theta_for_ssF = theta_path(end);
end

% theta_for_ssF = 0.4;
if cSF.pps_active
    % 设定PPS相关制度参数
    cSF.pps_contrib_rate = 0.03;
    cSF.pps_withdrawal_rate = 0.10;
    cSF.pps_tax_rate_withdrawal = 0.03;
    fprintf('   设定PPS参数: 缴费率=%.2f%%, 提取率=%.2f%%, 提取税率=%.2f%%\n', ...
        cSF.pps_contrib_rate*100, cSF.pps_withdrawal_rate*100, cSF.pps_tax_rate_withdrawal*100);
end

% 设定网格和收入过程
if cSF.pps_active
    cSF.nk = ngrid; cSF.nkpps = ngrid_pps; cSF.nkprime = ngrid; cSF.npps = ngrid_pps;
else
    cSF.nk = ngrid; cSF.nkpps = 1; cSF.nkprime = ngrid; cS.npps = 1;
end
cSF = model_setup_utils_bgp.generateGrids(cSF);
paramSF = struct();
[paramSF.leGridV, paramSF.TrProbM_by_age, paramSF.leProb1V, cSF.nw_expanded] = ...
    model_setup_utils_bgp.EarningProcess_AgeDependent(cSF);
paramSF.leLogGridV = log(paramSF.leGridV(1:cSF.nw));

% 打包求解ssF所需的外生参数
params_for_ssF = struct(...
    'Z', Z_path(:,end), ... % 内生计算稳态人口分布
    'A', 1.0, ...
    'theta', theta_for_ssF, ...
    'g_A_ss', cSF.g_A_ss,...
    'n_ss', cSF.n_ss);

% --- 步骤 4.2: 调用统一稳态求解器 ---
tic;
fprintf('   ⚙️  启动统一稳态求解器 (求解 ssF)...\n');
[ssF_result, distF, polF, valF] = ...
    main_steady_state_utils_bgp.solve_steady_state_unified(cSF, paramSF, params_for_ssF, true, [], 'lsqnonlin');
toc;

if isempty(ssF_result)
    warning('终期稳态(ssF)求解失败！');
else
    ssF = ssF_result;
    fprintf('✅ 终期稳态(ssF)求解成功！\n');
    fprintf('   利率和工资: r=%.4f, w=%.4f\n', ssF.r_mkt, ssF.w_hat);
end


%% --- 5. 结果对比与总结 ---
if ~isempty(fieldnames(ss0)) && ~isempty(fieldnames(ssF))
    fprintf('\n\n===========================================================\n');
    fprintf('===           比较静态分析: 关键宏观变量对比          ===\n');
    fprintf('===========================================================\n');
    fprintf('%-25s | %-15s | %-15s \n', '宏观变量', '初始稳态 (ss0)', '终期稳态 (ssF)');
    fprintf('-----------------------------------------------------------\n');
    fprintf('%-25s | %-15.4f | %-15.4f \n', '私人资本 (K_p_hat)', ss0.K_private_hat, ssF.K_private_hat);
    fprintf('%-25s | %-15.4f | %-15.4f \n', 'PPS资本 (K_pps_hat)', ss0.K_pps_hat, ssF.K_pps_hat);
    fprintf('%-25s | %-15.4f | %-15.4f \n', '公共资本 (K_g_hat)', ss0.K_public_hat, ssF.K_public_hat);
    fprintf('%-25s | %-15.4f | %-15.4f \n', '总劳动 (L_hat)', ss0.L_hat, ssF.L_hat);
    fprintf('%-25s | %-15.4f | %-15.4f \n', '总产出 (Y_hat)', ss0.Y_from_production_hat, ssF.Y_from_production_hat);
    fprintf('-----------------------------------------------------------\n');
    fprintf('%-25s | %-15.4f | %-15.4f \n', '市场利率 (r)', ss0.r_mkt, ssF.r_mkt);
    fprintf('%-25s | %-15.4f | %-15.4f \n', '标准化工资 (w_hat)', ss0.w_hat, ssF.w_hat);
    fprintf('%-25s | %-15.4f | %-15.4f \n', 'PAYG缴费率 (theta)', ss0.theta, ssF.theta);
    fprintf('===========================================================\n');
else
    fprintf('\n\n比较静态分析无法执行，因为至少一个稳态求解失败。\n');
end

fprintf('\n--- 比较静态分析脚本执行完毕 ---\n');