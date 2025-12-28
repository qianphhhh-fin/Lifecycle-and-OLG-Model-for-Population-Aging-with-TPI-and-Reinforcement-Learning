% =========================================================================
% == SCRIPT: test_perfect_ss.m
% == 版本: [v1.0 - 理论完美稳态检验]
% ==
% == 目的:
% ==   - 检验核心求解器 SS.solve_steady_state 的准确性。
% ==   - 使用一个理论上完美的、具有恒定人口增长率和生存率的稳态人口
% ==     分布，来求解模型的BGP稳态。
% ==   - 在这种理想环境下，所有市场的出清误差，特别是【投资缺口】，
% ==     都应该收敛到接近机器精度的水平。
% ==   - 如果测试通过，则证明模型的VFI、聚合和市场出清逻辑是正确的。
% =========================================================================
clear; close all; clear classes;
addpath(pwd);
fprintf('=== OLG模型核心求解器检验脚本 (v1.0) ===\n\n');

%% --- 1. 设置一套用于测试的“完美”稳态参数 ---
fprintf('--- 1. 设置理论稳态的经济学参数 ---\n');

% --- 获取基础参数结构体 ---
cS = utils.ParameterValues();

% --- 覆盖关键参数以进行清晰的测试 ---
cS.pps_active = false; % 在完美测试中，先关闭PPS以简化问题
cS.nkpps = 1;
cS.nk = 50;
cS.nw = 1; % 关闭异质性冲击，让所有家庭同质
cS.nTypes = 1; % 只使用一类家庭进行测试
cS.type_weights = 1; % 类型权重为1

% --- 设定一个清晰的BGP环境 ---
n_ss_annual = 0.01;  % 年化人口增长率: 1%
g_A_ss_annual = 0.015; % 年化技术进步率: 1.5%

cS.n_ss = n_ss_annual;
cS.g_A_ss = g_A_ss_annual;

% --- 使用一组恒定的生存率 ---
% 我们从cS中提取第一年的生存率数据，并假设它在稳态中永远不变
cS.s_pathV = cS.s_pathV(:, 1); 

fprintf('   年化人口增长率 (n)  : %.2f%%\n', n_ss_annual*100);
fprintf('   年化技术进步率 (g_A) : %.2f%%\n', g_A_ss_annual*100);
fprintf('   测试家庭类型 (nTypes) : %d\n', cS.nTypes);

%% --- 2. [核心] 生成理论完美的稳态人口分布 ---
fprintf('\n--- 2. 生成理论完美的稳态人口分布 ---\n');

% 使用你的工具函数，根据 n 和 s 计算理论分布
[Z_theory_ss, ~] = population.compute_theoretical_ss_dist(cS.s_pathV, cS.n_ss, cS.time_Step, cS.aD_new);

fprintf('   ✅ 理论人口分布已生成。\n');



%% --- 3. 准备求解器所需的其他输入 ---
fprintf('\n--- 3. 准备求解器输入 ---\n');

cS.theta_path_urban = 0.2; % 使用一个缴费率
cS.theta_path_resident = 0.03; % 只是占位符

% --- 生成网格和劳动过程 ---
cS = utils.generateGrids(cS);
paramS = struct();
[paramS.leGridV, paramS.TrProbM_by_age, paramS.leProb1V, cS.nw_expanded] = utils.EarningProcess_AgeDependent(cS);

% --- 打包传递给求解器的外部参数 ---
params_for_ss = struct(...
    'Z', Z_theory_ss, ... % [核心] 使用理论人口分布
    'A', 1.0, ...
    'g_A_ss', cS.g_A_ss, ...
    'n_ss', cS.n_ss);

fprintf('   ✅ 所有输入准备就绪。\n');

%% --- 4. 调用求解器并进行检验 ---
fprintf('\n--- 4. 启动稳态求解器进行测试 ---\n');

tic;
% 调用标准BGP求解器
[ssF, ~, ~, ~] = SS.solve_steady_state(cS, paramS, params_for_ss, true, true);
toc;

if isempty(ssF)
    error('测试失败：稳态求解器未能收敛！');
else
    fprintf('✅ 稳态求解器成功返回结果。\n');
    
    fprintf('\n--- 5. 生成并检验最终国民账户 ---\n');
    % 生成最终报告
    report_file_test = 'SS/perfect_ss_test_report.txt';
    report_struct = SS.display_national_accounts(ssF, cS, paramS, Z_theory_ss, report_file_test, true, true);
    
    % [!!! 最终检验 !!!]
    fprintf('\n\n==================================================\n');
    fprintf('###              核 心 准 确 性 检 验             ###\n');
    fprintf('==================================================\n');
    
    investment_gap_pct = report_struct.Invest_Gap_pct;
    fprintf('   投资缺口 (I_acct - I_bgp) / Y : %.4e %%\n', investment_gap_pct);
    
    if abs(investment_gap_pct) < 1e-4
        fprintf('\n   [✅✅✅ 测试通过 ✅✅✅]\n');
        fprintf('   投资缺口极小，证明求解器在BGP设定下是精确的。\n');
    else
        fprintf('\n   [❌❌❌ 测试失败 ❌❌❌]\n');
        fprintf('   投资缺口过大，表明VFI、聚合或市场出清逻辑中可能存在与\n');
        fprintf('   BGP不一致的错误。请重点检查：\n');
        fprintf('     - household.m 中的贴现和增长调整因子。\n');
        fprintf('     - aggregates.m 中的资本聚合逻辑。\n');
        fprintf('     - SS.m/system_of_equations.m 中的BGP资本积累方程。\n');
    end
    fprintf('==================================================\n\n');
end