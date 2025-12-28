% =========================================================================
% == SCRIPT: test_perfectN_SS_hhtype4.m
% == 版本: [v1.0 - 4类型家庭理想环境测试版]
% ==
% == 目的:
% ==   - 在一个理想化的、人口结构恒定的环境中，测试包含4类异质性
% ==     家庭的模型的稳态求解器能否成功收敛。
% ==   - 生成一个包含4类家庭的、内部完全自洽的稳态，作为后续分析
% ==     (如政策冲击比较) 的基准。
% ==   - 使用与同质性测试完全相同的人口动态参数，以隔离家庭异质性
% ==     带来的影响。
% =========================================================================
clear; close all; clear classes;
addpath(pwd); 
fprintf('=== 4类型异质性家庭理想环境稳态求解脚本 (v1.0) ===\n\n');

%% --- 1. 全局设置与初始化 ---
report_filename = '稳态国民经济核算报告_理论环境_HET_n4.txt';
output_data_filename = 'SS/data_for_perfect_transition_HET_n4.mat';

ngrid_k = 40;
ngrid_kpps = 1;

fprintf('   加载模型物理参数 (含4类异质性家庭)...\n');
cS = utils.ParameterValues();

% --- [核心设定] 确保加载了4种家庭类型 ---
if cS.nTypes ~= 4
    warning('ParameterValues 未能正确加载4种家庭类型，检测到 nTypes = %d', cS.nTypes);
end
fprintf('   [配置] 模型已加载 %d 类异质性家庭。\n', cS.nTypes);

cS.pps_active = false;
fprintf('   [配置] PPS模块已禁用。\n');

fprintf('   设定并生成全局共享的资产网格...\n');
cS.nk = ngrid_k; cS.nkpps = ngrid_kpps;
cS = utils.generateGrids(cS);
cS = utils.calcaulte_theta_payg_path(cS, false);
cS.tau_k = 0.03; cS.tau_l = 0.06; cS.tau_c = 0.05;
%% --- 2. 构建理想化稳态环境 (与同质性测试完全一致) ---
fprintf('\n--- 2. 构建理想化稳态环境 ---\n');
cSF_test = cS; 
cSF_test.g_A_ss = 0.015;
cSF_test.I_g_to_Y_ratio_ss = 0.03;
cSF_test.n_ss = 0.01;
fprintf('   [核心设定] 理想环境的人口年增长率 n_ss 被强制设为: %.4f\n', cSF_test.n_ss);

% 使用最后一年（稳态）的生存率数据
s_pathV_ss = cSF_test.s_pathV(:, end);
cSF_test.s_pathV = repmat(s_pathV_ss, 1, size(cS.s_pathV, 2)); % 确保维度一致

% 根据生存率和增长率，计算一个完全自洽的理论稳态人口分布
Z_theory = population.compute_theoretical_ss_dist(s_pathV_ss, cSF_test.n_ss, cSF_test.time_Step, cSF_test.aD_new);
fprintf('   ✅ 已根据 n_ss=%.4f 和 s_pathV 重新计算出完全自洽的理论人口分布 (Z_theory)。\n', cSF_test.n_ss);

%% --- 3. 在理想环境下求解稳态 (使用异质性求解器) ---
fprintf('\n\n--- 3. 在4类型家庭理想环境下求解稳态 (ssF_perfect_het_n4) ---\n');

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
    error('理想化稳态(ssF_perfect_het_n4)求解失败！脚本终止。');
else
    fprintf('✅ 理想化稳态(ssF_perfect_het_n4)求解成功！\n');
end

%% --- 4. 结果分析与最终检验 ---
% 调用能处理异质性结果的 display_national_accounts
SS.display_national_accounts(ssF_perfect, cSF_test, paramSF_test, Z_theory, report_filename, true, true, distF_by_type, polF_by_type);

%% --- 5. 保存所有必要数据到 .mat 文件 ---
fprintf('\n\n--- 5. 保存结果用于后续分析 ---\n');
try
    data_for_perfect_transition_het_n4 = struct(...
        'ssF', ssF_perfect, ...
        'distF_by_type', distF_by_type, ...
        'polF_by_type', polF_by_type, ...
        'valF_by_type', valF_by_type, ...
        'cS', cSF_test, ...
        'paramSF', paramSF_test ...
    );
    
    save(output_data_filename, 'data_for_perfect_transition_het_n4', '-v7.3');
    fprintf('✅ 所有4类型家庭稳态数据已成功保存至: %s\n', output_data_filename);
catch ME
    warning('保存数据时发生错误！');
    fprintf('错误信息: %s\n', ME.message);
end

fprintf('\n--- 4类型异质性家庭理想环境稳态求解脚本执行完毕 ---\n');