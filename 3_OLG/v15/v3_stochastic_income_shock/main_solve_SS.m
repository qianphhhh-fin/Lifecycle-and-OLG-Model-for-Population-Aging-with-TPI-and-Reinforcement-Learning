% =========================================================================
% == SCRIPT: main_solve_SS.m [vPaper.4 - 支出冲击模型]
% == 目的: 实现一个基于年龄依赖的、异质性的重大支出冲击模型，
% ==         以增强家庭的预防性储蓄动机，并观察对资本积累(K/Y)的影响。
% == 核心模型: 冲击体现为预算约束中的大额支出，而非收入下降。
% =========================================================================
clear; close all;                           % 清空工作空间和关闭图形窗口
addpath(pwd);                               % 将当前目录添加到MATLAB路径
fprintf('=== [支出冲击版] 稳态模型求解与国民核算报告 (遗赠税模型) ===\n\n');

%% 1. 初始化环境
fprintf('--- 1. 初始化环境并生成模拟所需输入 ---\n');

% [核心变更] 调用包含支出冲击参数的新版参数设置函数
cS = main_olg_v15_utils.ParameterValues_ExpenditureShock(); 
paramS = struct();
fprintf('   稳态财政设定：政府预算平衡 (G=T), 意外遗赠被政府100%%征收。\n');
cS.gov_debt_frac_Y = 0;

% 网格设定
ngrid = 100; cS.nk = ngrid; cS.nkpps = ngrid; cS.nkprime = ngrid; cS.npps = 5;
cS = main_olg_v15_utils.generateGrids(cS);

% [核心变更] 调用新的、生成年龄依赖冲击“信号”过程的函数
[paramS.leGridV, paramS.TrProbM_by_age, paramS.leProb1V, cS.nw_expanded] = ...
    main_olg_v15_utils.EarningProcess_AgeDependent(cS);
paramS.leLogGridV = log(paramS.leGridV(1:cS.nw)); % Log只对常规状态有意义
fprintf('   已生成年龄依赖的支出冲击“信号”过程。\n');

% 1. 基于存活率计算稳态年龄人口质量
age_mass = ones(cS.aD_new, 1);
for ia = 1:(cS.aD_new - 1)
    age_mass(ia+1) = age_mass(ia) * cS.s_pathV(ia);
end

% 2. 归一化，得到每个年龄组占总人口的比例
Z_ss_model_implied = age_mass / sum(age_mass);

% 3. 用这个内生的、自洽的分布，替换掉从外部加载的 Z_path(:,1)
Z_ss_norm = Z_ss_model_implied;
fprintf('   已生成模型内生稳态人口分布，将用此分布求解。\n');

%% 2. 求解稳态 (资本市场出清法)
fprintf('\n--- 2. 调用稳态求解器 (以 资产供给 = 资产需求 为目标) ---\n');

% [核心变更] 调用已适配支出冲击模型的静态求解器方法
[ss, eq_found, Dist, k_prime_idx] = main_olg_v15_utils.solve_steady_state_iter_K(Z_ss_norm, cS, paramS);
if ~eq_found, error('资本市场出清法求解失败'); end

%% 3. 结果验证与展示
fprintf('\n\n--- 3. 结果验证与展示 ---\n');
fprintf('✅✅✅ 求解成功！均衡资本 K = %.8f ✅✅✅\n', ss.K_physical);

% [核心变更] 调用适配支出冲击模型的新版展示函数
main_olg_v15_utils.display_national_accounts_expenditure_shock(ss, Dist, k_prime_idx, cS, paramS);