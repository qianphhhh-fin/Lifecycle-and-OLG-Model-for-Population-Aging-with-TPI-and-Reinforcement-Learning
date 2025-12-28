% =========================================================================
% == SCRIPT: test_perfectN_SS.m
% == 版本: [v2.0 - 异质性模型验证版]
% ==
% == 目的:
% ==   - 验证经过异质性改造后的模型，在 nTypes=1 的设定下，
% ==     能否精确复现原同质模型的“黄金标准”稳态结果。
% =========================================================================
clear; close all; clear classes;
addpath(pwd); 
fprintf('=== 异质性模型退化验证脚本 (v2.0) ===\n\n');

%% --- 1. 全局设置与初始化 ---
report_filename = '稳态国民经济核算报告_理论环境_HET_n1.txt';
output_data_filename = 'SS/data_for_perfect_transition_HET_n1.mat';

ngrid_k = 50;
ngrid_kpps = 1;

fprintf('   加载模型物理参数 (异质性版本)...\n');
cS = old_utils.ParameterValues();
% --- [核心验证设定] 将模型强制设为单类型同质家庭 ---
cS.nTypes = 1;
cS.type_weights = 1;
% 确保所有类型相关参数维度为1
cS.beta = cS.beta(1);
cS.phi_bequest = cS.phi_bequest(1);

% cS = utils.age_efficiency(cS);
cS.ageEff_by_type = cS.ageEffV_new';
cS.theta_path_urban = cS.theta_path;
% 默认使用城镇职工的养老金路径作为同质代表
fprintf('   [验证设定] 家庭类型 nTypes 已强制设为 1。\n');

cS.pps_active = false;
fprintf('   [配置] PPS模块已禁用。\n');

fprintf('   设定并生成全局共享的资产网格...\n');
cS.nk = ngrid_k; cS.nkpps = ngrid_kpps; cS.nkprime = ngrid_k;
cS = utils.generateGrids(cS);

cS.tau_k = 0.03; cS.tau_l = 0.06; cS.tau_c = 0.05;
%% --- 2. 构建理想化稳态环境 (与原脚本一致) ---
fprintf('\n--- 2. 构建理想化稳态环境 ---\n');
cSF_test = cS; 
cSF_test.g_A_ss = 0.015;
cSF_test.I_g_to_Y_ratio_ss = 0.03;
cSF_test.n_ss = 0.01; % [核心设定] 强制人口零增长以简化测试
fprintf('   [核心设定] 理想环境的人口年增长率 n_ss 被强制设为: %.4f\n', cSF_test.n_ss);

s_pathV_ss = cSF_test.s_pathV(:, end);
cSF_test.s_pathV = s_pathV_ss;

Z_theory = population.compute_theoretical_ss_dist(cSF_test.s_pathV(:,1), cSF_test.n_ss, cSF_test.time_Step, cSF_test.aD_new);
fprintf('   ✅ 已根据 n_ss=%.4f 和 s_pathV 重新计算出完全自洽的理论人口分布 (Z_theory)。\n', cSF_test.n_ss);

%% --- 3. 在理想环境下求解稳态 (使用异质性求解器) ---
fprintf('\n\n--- 3. 在理想环境下求解稳态 (ssF_perfect_het_n1) ---\n');

paramSF_test = struct();
[paramSF_test.leGridV, paramSF_test.TrProbM_by_age, paramSF_test.leProb1V, cSF_test.nw_expanded] = utils.EarningProcess_AgeDependent(cSF_test);

params_for_ss_solve = struct(...
    'Z', Z_theory, ...
    'A', 1.0, ...
    'g_A_ss', cSF_test.g_A_ss, ...
    'n_ss', cSF_test.n_ss);

tic;
fprintf('   ⚙️  启动统一稳态求解器 (调用异质性版本的 SS.solve_steady_state)...\n');
[ssF_perfect, distF_by_type, polF_by_type, valF_by_type] = SS.solve_steady_state(cSF_test, paramSF_test, params_for_ss_solve, true, true);
toc;

if isempty(ssF_perfect)
    error('理想化稳态(ssF_perfect_het_n1)求解失败！脚本终止。');
else
    fprintf('✅ 理想化稳态(ssF_perfect_het_n1)求解成功！\n');
end

%% --- 4. 结果分析与最终检验 ---
% [核心验证] 调用新的、能处理异质性结果的 display_national_accounts
SS.display_national_accounts(ssF_perfect, cSF_test, paramSF_test, Z_theory, report_filename, true, true, distF_by_type, polF_by_type);

%% --- 5. 保存所有必要数据到 .mat 文件 ---
fprintf('\n\n--- 5. 保存结果用于TPI诊断 ---\n');
try
    data_for_perfect_transition_het_n1 = struct(...
        'ssF', ssF_perfect, ...
        'distF_by_type', distF_by_type, ...
        'polF_by_type', polF_by_type, ...
        'valF_by_type', valF_by_type, ...
        'cS', cSF_test, ...
        'paramSF', paramSF_test ...
    );
    
    save(output_data_filename, 'data_for_perfect_transition_het_n1', '-v7.3');
    fprintf('✅ 所有TPI所需的数据已成功保存至: %s\n', output_data_filename);
    fprintf('\n>>> 请将此文件的报告与原同质模型的报告进行对比，以验证修改的正确性。<<<\n');
catch ME
    warning('保存数据时发生错误！');
    fprintf('错误信息: %s\n', ME.message);
end

fprintf('\n--- 异质性模型退化验证脚本执行完毕 ---\n');