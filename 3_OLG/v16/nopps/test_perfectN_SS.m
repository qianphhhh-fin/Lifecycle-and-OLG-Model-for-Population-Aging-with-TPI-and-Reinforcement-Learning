% =========================================================================
% == SCRIPT: main_run_SS_nopp_TEST.m
% == 版本: [v1.0 - 健康人口参数测试版]
% ==
% == 目的:
% ==   在一个由正出生率和低死亡率决定的“健康”稳态人口结构下，
% ==   测试最终版模型代码(v9.0)的收敛性和自洽性。
% ==   如果在此环境下模型能够完美求解，则证明代码本身是正确的，
% ==   而原脚本的收敛问题来自于极端老龄化的参数设定。
% =========================================================================
clear; close all; clear classes;
addpath(pwd);
fprintf('=== OLG模型(无PPS) - 健康人口环境测试脚本 ===\n\n');

%% --- 1. 全局设置与初始化 ---
fprintf('--- 1. 全局设置与初始化 ---\n');

ngrid = 50;
ngrid_pps = 40;

fprintf('   加载模型物理参数...\n');
cS = model_setup_utils_bgp.ParameterValues();
fprintf('   设定并生成全局共享的资产网格...\n');
cS.nk = ngrid; 
cS.nkpps = ngrid_pps; 
cS.nkprime = ngrid; 
cS.npps = ngrid_pps;
cS = model_setup_utils_bgp.generateGrids(cS);
fprintf('   ✅ 全局网格(kGridV等)已生成并存入cS。\n');

% --- 设定政策情景 ---
cS.endogenous_theta_mode = true;
cS.pps_active = false;
cS.payg_replacement_rate = 0.4; % 维持政策目标不变
cS.theta_max = 0.99;

%% --- 2. [核心] 构建一个“健康”的理想化稳态环境 ---
fprintf('\n--- 2. 构建理想化稳态环境 ---\n');

% --- 步骤 2.1: 设定长期经济增长率 ---
cSF_test = cS; % 创建一个用于测试的参数副本
cSF_test.g_A_ss = 0.015;  % 设定一个健康的长期年化技术进步率 (1.5%)

% --- 步骤 2.2: 设定一个能导致人口正增长的出生率 ---
% 这里的出生率是一个校准工具，我们需要它能抵消死亡并实现正增长
% 我们不直接使用出生率，而是设定一个目标长期人口年化增长率
cSF_test.n_ss = 0.01;   % 设定健康的长期年化人口增长率 (1%)

% --- 步骤 2.3: 使用模型最低的死亡率作为稳态死亡率 ---
% cS.s_pathV 的最后一列代表了2100年的存活率，通常是最低的死亡率
cSF_test.s_pathV = cS.s_pathV(:, end);

fprintf('   理想环境设定:\n');
fprintf('     - 年化技术增长率 (g_A): %.3f\n', cSF_test.g_A_ss);
fprintf('     - 年化人口增长率 (n):   %.3f\n', cSF_test.n_ss);

% --- 步骤 2.4: 计算理论稳态人口分布 Z_theory ---
% 我们需要一个新的辅助函数来完成这个计算。
fprintf('   正在计算对应的理论稳态人口分布 Z_theory ...\n');
Z_theory = local_compute_theoretical_ss_dist(cSF_test.s_pathV, cSF_test.n_ss, cSF_test.time_Step, cSF_test.aD_new);
fprintf('   ✅ 理论人口分布计算完成。\n');

% 可选：可视化检查新的人口金字塔
% figure('Name', '健康稳态人口结构');
% barh(Z_theory); title('理论稳态人口分布 (n > 0)'); ylabel('年龄组'); xlabel('占比');

%% --- 3. 在理想环境下求解稳态 ---
fprintf('\n\n--- 3. 在理想环境下求解稳态 (ssF_test) ---\n');

paramSF_test = struct();
[paramSF_test.leGridV, paramSF_test.TrProbM_by_age, paramSF_test.leProb1V, cSF_test.nw_expanded] = model_setup_utils_bgp.EarningProcess_AgeDependent(cSF_test);
paramSF_test.leLogGridV = log(paramSF_test.leGridV(1:cSF_test.nw));
params_for_ssF_test = struct('Z', Z_theory, 'A', 1.0, 'theta', [], 'g_A_ss', cSF_test.g_A_ss, 'n_ss', cSF_test.n_ss);

tic;
fprintf('   ⚙️  启动统一稳态求解器 (求解 ssF_test)...\n');
% [重要] 确保这里调用的是最终的、完全修正的求解器代码 (v9.0)
[ssF_test, distF_test, polF_test, ~] = main_steady_state_utils_bgp.solve_steady_state_unified(cSF_test, paramSF_test, params_for_ssF_test, true, [], 'lsqnonlin');
toc;

if isempty(ssF_test)
    warning('测试稳态(ssF_test)求解失败！');
    return;
else
    fprintf('✅ 测试稳态(ssF_test)求解成功！ r=%.4f, w=%.4f\n', ssF_test.r_mkt, ssF_test.w_hat);
end

%% --- 4. 结果分析与最终检验 ---
if ~isempty(fieldnames(ssF_test))
    fprintf('\n\n============================================================\n');
    fprintf('===         测试稳态 (ssF_test) 结果分析与检验         ===\n');
    fprintf('============================================================\n');
    
    fprintf('--- [核心] 资本存量稳态自洽性检验 ---\n');
    fprintf('   稳态存量是否自洽 (期初总量 == 期末总量)？\n');
    fprintf('   K_p_hat_begin (求解器均衡点)       : %15.8f\n', ssF_test.K_private_begin_hat);
    fprintf('   K_p_hat_end (模型聚合期末选择)     : %15.8f\n', ssF_test.K_private_hat);
    ratio = ssF_test.K_private_hat / ssF_test.K_private_begin_hat;
    fprintf('   -> 两者之比 (End / Begin)           : %15.8f\n', ratio);

    consistency_error = ratio - 1;
    fprintf('   -> 一致性误差 ((比率)-1)          : %14.2f%%\n', consistency_error * 100);

    if abs(consistency_error) < 1e-4 % 容忍度设为 0.01%
        fprintf('   \n   检验结果: ✅✅✅ 代码通过验证！模型在健康参数下完美收敛且自洽。\n');
    else
        fprintf('   \n   检验结果: ⚠️  代码仍存在深层问题。即使在理想环境下也无法实现存量自洽。\n');
    end

    fprintf('============================================================\n');
else
    fprintf('\n\n结果分析无法执行，因为测试稳态(ssF_test)求解失败。\n');
end

fprintf('\n--- 健康人口参数测试脚本执行完毕 ---\n');


%% --- 本地辅助函数 ---
function Z_theory = local_compute_theoretical_ss_dist(s_pathV_ss, n_ss_annual, time_Step, aD_new)
    % 功能: 计算一个具有恒定人口年增长率n的理论稳态人口分布。
    
    % 1. 将年化增长率转换为模型期的增长率
    n_period = (1 + n_ss_annual)^time_Step - 1;

    % 2. 以新生儿为基准1.0，计算各年龄组的相对人口规模
    mass_levels_by_age = zeros(aD_new, 1);
    mass_levels_by_age(1) = 1.0;
    for ia = 1:(aD_new - 1)
        % 核心逻辑：下一代的人口规模由上一代的存活者构成，但相对于
        % 新一代的新生儿，他们的相对规模因人口整体增长而被“稀释”。
        mass_levels_by_age(ia+1) = (mass_levels_by_age(ia) * s_pathV_ss(ia)) / (1 + n_period);
    end

    % 3. 将相对规模归一化，得到最终的概率分布
    Z_theory = mass_levels_by_age / sum(mass_levels_by_age);
end