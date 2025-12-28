% =========================================================================
% == SCRIPT: test_perfectN_SS.m
% == 版本: [v1.2 - 黄金标准生成最终版]
% ==
% == 目的:
% ==   - 生成一个内部完全自洽的稳态，用作TPI诊断的“黄金标准”。
% ==   - 确保保存的数据 (cS, paramSF) 包含完全一致的人口动态参数。
% =========================================================================
clear; close all; clear classes;
addpath(pwd); 
fprintf('=== 理想化环境稳态求解与保存脚本 (v1.2) ===\n\n');

%% --- 1. 全局设置与初始化 ---
report_filename = '稳态国民经济核算报告_理论环境.txt';
output_data_filename = 'SS/data_for_perfect_transition.mat';

ngrid_k = 60;
ngrid_kpps = 20;

fprintf('   加载模型物理参数...\n');
cS = utils.ParameterValues();
cS.pps_active = true;
fprintf('   [配置] PPS模块已禁用。\n');

fprintf('   设定并生成全局共享的资产网格...\n');
cS.nk = ngrid_k; cS.nkpps = ngrid_kpps; cS.nkprime = ngrid_k;
cS = utils.generateGrids(cS);

cS.tau_k = 0.03; cS.tau_l = 0.06; cS.tau_c = 0.05;

%% --- 2. 构建理想化稳态环境 ---
fprintf('\n--- 2. 构建理想化稳态环境 ---\n');
cSF_test = cS; 
cSF_test.g_A_ss = 0.015;
cSF_test.I_g_to_Y_ratio_ss = 0.03;
cSF_test.n_ss = -0.01; % [核心设定] 强制人口零增长以简化测试
fprintf('   [核心设定] 理想环境的人口年增长率 n_ss 被强制设为: %.4f\n', cSF_test.n_ss);

% 使用最后一年（稳态）的生存率数据，并确保它是一个列向量
s_pathV_ss = cSF_test.s_pathV(:, end);
cSF_test.s_pathV = s_pathV_ss; % 强制 cS 结构体内部也使用这个唯一的生存率向量

% [核心逻辑] 根据生存率和零增长率，重新计算一个完全自洽的理论稳态人口分布
Z_theory = utils.compute_theoretical_ss_dist(cSF_test.s_pathV, cSF_test.n_ss, cSF_test.time_Step, cSF_test.aD_new);
fprintf('   ✅ 已根据 n_ss=%.4f 和 s_pathV 重新计算出完全自洽的理论人口分布 (Z_theory)。\n', cSF_test.n_ss);

%% --- 3. 在理想环境下求解稳态 ---
fprintf('\n\n--- 3. 在理想环境下求解稳态 (ssF_perfect) ---\n');

paramSF_test = struct();
[paramSF_test.leGridV, paramSF_test.TrProbM_by_age, paramSF_test.leProb1V, cSF_test.nw_expanded] = utils.EarningProcess_AgeDependent(cSF_test);

params_for_ss_solve = struct(...
    'Z', Z_theory, ...
    'A', 1.0, ...
    'g_A_ss', cSF_test.g_A_ss, ...
    'n_ss', cSF_test.n_ss);

tic;
fprintf('   ⚙️  启动统一稳态求解器 (调用 SS.solve_steady_state)...\n');
[ssF_perfect, distF_perfect, polF_perfect, valF_perfect] = SS.solve_steady_state(cSF_test, paramSF_test, params_for_ss_solve, true, 'lsqnonlin');
toc;

if isempty(ssF_perfect)
    error('理想化稳态(ssF_perfect)求解失败！脚本终止。');
else
    fprintf('✅ 理想化稳态(ssF_perfect)求解成功！\n');
end

%% --- 4. 结果分析与最终检验 ---
SS.display_national_accounts(ssF_perfect, cSF_test, paramSF_test, Z_theory, report_filename);

%% --- 5. 保存所有必要数据到 .mat 文件 ---
fprintf('\n\n--- 5. 保存结果用于TPI诊断 ---\n');
try
    data_for_perfect_transition = struct(...
        'ssF', ssF_perfect, ...
        'distF', distF_perfect, ...
        'polF', polF_perfect, ...
        'valF', valF_perfect, ...
        'cS', cSF_test, ...
        'paramSF', paramSF_test ...
    );
    
    save(output_data_filename, 'data_for_perfect_transition', '-v7.3');
    fprintf('✅ 所有TPI所需的数据已成功保存至: %s\n', output_data_filename);
catch ME
    warning('保存数据时发生错误！');
    fprintf('错误信息: %s\n', ME.message);
end

fprintf('\n--- 理想化环境稳态求解脚本执行完毕 ---\n');
