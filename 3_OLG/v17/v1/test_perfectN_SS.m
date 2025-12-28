% =========================================================================
% == SCRIPT: test_perfectN_SS_with_PPS.m (示例)
% == 版本: [v1.0 - PPS 激活版]
% =========================================================================
clear; close all; clear classes;
addpath(pwd); 
fprintf('=== OLG模型 (带PPS) - 模块化健康人口环境测试脚本 ===\n\n');

%% --- 1. 全局设置与初始化 ---
report_filename = '稳态国民经济核算报告_带PPS.txt';

ngrid_k = 40; % 常规资产网格点数
ngrid_kpps = 40; % PPS资产网格点数

fprintf('   加载模型物理参数...\n');
cS = utils.ParameterValues();

% --- [核心修改] 激活并配置PPS ---
cS.pps_active = true;
cS.pps_contrib_rate = 0.03;        % PPS的强制缴费率 (例如 8%)
cS.pps_withdrawal_rate = 0.1;     % 退休后每个模型期强制提取PPS资产的比例 (例如 5%)
cS.pps_tax_rate_withdrawal = 0.03; % 提取PPS资产时适用的税率 (例如 10%)
fprintf('   [配置] PPS模块已激活！缴费率=%.2f, 提款率=%.2f, 税率=%.2f\n', ...
    cS.pps_contrib_rate, cS.pps_withdrawal_rate, cS.pps_tax_rate_withdrawal);

fprintf('   设定并生成全局共享的资产网格...\n');
cS.nk = ngrid_k;
cS.nkpps = ngrid_kpps; % 设定PPS网格点数
cS.nkprime = ngrid_k;
cS = utils.generateGrids(cS); % generateGrids现在会生成两个网格

% 设定政策情景
cS.tau_k = 0.25; 
cS.tau_l = 0.20; 
cS.tau_c = 0.05; 

%% --- 2. 构建理想化稳态环境 (与之前相同) ---
fprintf('\n--- 2. 构建理想化稳态环境 ---\n');
cSF_test = cS; 
cSF_test.g_A_ss = 0.015;
cSF_test.I_g_to_Y_ratio_ss = 0.03; 
cSF_test.n_ss = 0.01;
cSF_test.s_pathV = cS.s_pathV(:, end);
Z_theory = utils.compute_theoretical_ss_dist(cSF_test.s_pathV, cSF_test.n_ss, cSF_test.time_Step, cSF_test.aD_new);

%% --- 3. 在理想环境下求解稳态 ---
fprintf('\n\n--- 3. 在理想环境下求解带PPS的稳态 (ssF_test) ---\n');

paramSF_test = struct();
[paramSF_test.leGridV, paramSF_test.TrProbM_by_age, paramSF_test.leProb1V, cSF_test.nw_expanded] = utils.EarningProcess_AgeDependent(cSF_test);

params_for_ssF_test = struct(...
    'Z', Z_theory, ...
    'A', 1.0, ...
    'g_A_ss', cSF_test.g_A_ss, ...
    'n_ss', cSF_test.n_ss);

tic;
fprintf('   ⚙️  启动统一稳态求解器 (调用 SS.solve_steady_state)...\n');
% --- 调用求解器，它现在会自动处理PPS逻辑 ---
[ssF_test, ~, ~, ~] = SS.solve_steady_state(cSF_test, paramSF_test, params_for_ssF_test, true, 'lsqnonlin');
toc;

if isempty(ssF_test)
    warning('带PPS的测试稳态(ssF_test)求解失败！脚本终止。');
    return;
else
    fprintf('✅ 带PPS的测试稳态(ssF_test)求解成功！\n');
end

%% --- 4. 结果分析与最终检验 (与之前相同) ---
if ~isempty(fieldnames(ssF_test))
    SS.display_national_accounts(ssF_test, cSF_test, paramSF_test, Z_theory, report_filename);
else
    fprintf('\n\n结果分析无法执行，因为测试稳态(ssF_test)求解失败。\n');
end

fprintf('\n--- 带PPS的测试脚本执行完毕 ---\n');