% =========================================================================
% == SCRIPT: test_perfectn_ss_hhsame.m
% == 版本: [v1.0 - 异质性框架验证版]
% ==
% == 目的:
% ==   - 在异质性家庭框架下，设定所有家庭类型具有完全相同的参数。
% ==   - 求解其稳态，并与同质代理人模型的“黄金标准”结果进行对比。
% ==   - 验证异质性框架在参数相同时能够正确退化为同质代理人情况。
% =========================================================================
clear; close all; clear classes;
addpath(pwd); 
fprintf('=== 异质性框架验证稳态求解脚本 (v1.0) ===\n\n');

%% --- 1. 全局设置与初始化 ---
report_filename = '稳态国民经济核算报告_HHSAME.txt';
output_data_filename = 'SS/data_for_hhsame_transition.mat';

ngrid_k = 50;
ngrid_kpps = 1;

fprintf('   加载模型物理参数...\n');
cS = utils.ParameterValues();
% 计算年龄效率剖面
ages_eff = (cS.age1_orig : cS.ageLast_orig)';              % 获取用于计算效率的物理年龄向量。
beta0=-13.215788; beta1=1.349514; beta2=-0.043363; beta3=0.000585; beta4=-0.000003; % 效率函数的多项式系数。
log_age_eff_orig = beta0 + beta1*ages_eff + beta2*(ages_eff.^2) + beta3*(ages_eff.^3) + beta4*(ages_eff.^4);
ageEffV_orig_unnormalized = exp(log_age_eff_orig);          % 计算未归一化的年龄效率。
ageEffV_orig_unnormalized((cS.aR_idx_orig):end) = 0;        % 退休后效率为0。
mean_efficiency_working = mean(ageEffV_orig_unnormalized(1:(cS.aR_idx_orig - 1))); % 计算工作期的平均效率用于归一化。
ageEffV_orig = ageEffV_orig_unnormalized / mean_efficiency_working; % 归一化，使得工作期平均效率为1。
cS.ageEffV_new = zeros(cS.aD_new, 1);                       % 初始化模型年龄组的平均效率向量。
for a = 1:cS.aD_new
    phys_ages_in_group = cS.physAgeMap{a};
    model_ages_indices = phys_ages_in_group - cS.age1_orig + 1;
    cS.ageEffV_new(a) = mean(ageEffV_orig(model_ages_indices));
end

cS.type_weights
cS.pps_active = false;
fprintf('   [配置] PPS模块已禁用。\n');

fprintf('   设定并生成全局共享的资产网格...\n');
cS.nk = ngrid_k; cS.nkpps = ngrid_kpps; cS.nkprime = ngrid_k;

cS = utils.generateGrids(cS);

cS.tau_k = 0.03; cS.tau_l = 0.06; cS.tau_c = 0.05;

%% --- 2. 构建理想化稳态环境 (核心修改) ---
fprintf('\n--- 2. 构建异质性但参数相同的理想化环境 ---\n');
cSF_test = cS; 
cSF_test.g_A_ss = 0.015;
cSF_test.I_g_to_Y_ratio_ss = 0.03;
cSF_test.n_ss = 0.01; 

fprintf('   [核心设定] 理想环境的人口年增长率 n_ss 被设为: %.4f\n', cSF_test.n_ss);

% --- [!!! 核心异质性设定 !!!] ---
cSF_test.num_hh_types = 4;
nH = cSF_test.num_hh_types;
% --- [!!! 人口分布处理 !!!] ---
s_pathV_ss = cSF_test.s_pathV(:, end);
cSF_test.s_pathV = s_pathV_ss; 
% 计算总体的理论人口分布
Z_theory_total = population.compute_theoretical_ss_dist(cSF_test.s_pathV, cSF_test.n_ss, cSF_test.time_Step, cSF_test.aD_new);
% 将总分布按比例分配给各类家庭，得到 [aD x nH] 矩阵
Z_theory_h = Z_theory_total * cSF_test.type_weights';
cSF_test.theta_path_h = repmat(cSF_test.theta_path, [nH, 1]);
cSF_test.ageEffV_new_h = repmat(cSF_test.ageEffV_new, [1, nH]);
fprintf('   ✅ 已根据 n_ss 和 s_pathV 计算出理论人口分布，并按比例分配给 %d 类家庭。\n', nH);

%% --- 3. 在理想环境下求解稳态 ---
fprintf('\n\n--- 3. 在理想环境下求解稳态 (ssF_hhsame) ---\n');

paramSF_test = struct();
[paramSF_test.leGridV, paramSF_test.TrProbM_by_age, paramSF_test.leProb1V, cSF_test.nw_expanded] = utils.EarningProcess_AgeDependent(cSF_test);

params_for_ss_solve = struct(...
    'Z', Z_theory_h, ... % 传入矩阵形式的分布
    'A', 1.0, ...
    'g_A_ss', cSF_test.g_A_ss, ...
    'n_ss', cSF_test.n_ss);

tic;
fprintf('   ⚙️  启动统一稳态求解器 (调用 SS.solve_steady_state)...\n');
[ssF_hhsame, distF_hhsame, polF_hhsame, valF_hhsame] = SS.solve_steady_state(cSF_test, paramSF_test, params_for_ss_solve, true, 'lsqnonlin');
toc;

if isempty(ssF_hhsame)
    error('理想化稳态(ssF_hhsame)求解失败！脚本终止。');
else
    fprintf('✅ 理想化稳态(ssF_hhsame)求解成功！\n');
end

%% --- 4. 结果分析与最终检验 ---
% 调用新的、能处理异质性家庭的报告函数
SS.display_national_accounts(ssF_hhsame, cSF_test, paramSF_test, Z_theory_h, report_filename, true, true);

fprintf('\n--- 验证脚本执行完毕 ---\n');
fprintf('下一步：请手动比较 %s 与 test_perfectN_SS.m 生成的报告。\n', report_filename);
fprintf('如果两个报告中的宏观总量完全一致，则证明异质性框架构建正确。\n');