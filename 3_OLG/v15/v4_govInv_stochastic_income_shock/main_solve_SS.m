% =========================================================================
% == SCRIPT: main_solve_SS.m [vPaper.5 - 政府投资版]
% == 目的: 在支出冲击模型基础上，引入生产性政府投资，
% ==         分析其对资本积累(K/Y)的影响。
% == 核心模型: 政府支出分为消费性支出(Gc)和投资性支出(Ig)。
% =========================================================================
clear; close all;                           % 清空工作空间和关闭图形窗口
addpath(pwd);                               % 将当前目录添加到MATLAB路径
fprintf('=== [政府投资版] 稳态模型求解与国民核算报告 ===\n\n');

%% 1. 初始化环境
fprintf('--- 1. 初始化环境并生成模拟所需输入 ---\n');

% 调用包含政府投资新参数的参数设置函数
cS = main_olg_v15_utils.ParameterValues_ExpenditureShock(); 
paramS = struct();
fprintf('   稳态财政设定：政府投资开启，意外遗赠被政府100%%征收。\n');
cS.gov_debt_frac_Y = 0;

% 网格设定
ngrid = 100; cS.nk = ngrid; cS.nkpps = ngrid; cS.nkprime = ngrid; cS.npps = 5;
cS = main_olg_v15_utils.generateGrids(cS);

% 生成年龄依赖冲击“信号”过程
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

%% 2. 求解稳态 (资本市场出清法 - fsolve)
fprintf('\n--- 2. 调用稳态求解器 (求解私人与公共资本联合均衡) ---\n');

% [核心变更] 调用已适配政府投资模型的新版求解器
[ss, eq_found, Dist, k_prime_idx] = main_olg_v15_utils.solve_steady_state_iter_Kg(Z_ss_norm, cS, paramS);
% if ~eq_found, error('资本市场出清法(fsolve)求解失败'); end

%% 3. 结果验证与展示
fprintf('\n\n--- 3. 结果验证与展示 ---\n');
fprintf('✅✅✅ 求解成功！均衡私人资本 K_p = %.6f, 公共资本 K_g = %.6f ✅✅✅\n', ss.K_private, ss.K_public);
fprintf('✅✅✅ 总资本 K_total = %.6f ✅✅✅\n', ss.K_total);


% [核心变更] 调用适配政府投资模型的新版展示函数
main_olg_v15_utils.display_national_accounts_gov_investment(ss, Dist, k_prime_idx, cS, paramS);

%% 4. 绘制年龄剖面图
fprintf('\n--- 4. 绘制年龄剖面图 ---\n');
main_olg_v15_utils.plot_saving_rate_by_age(ss, Dist, k_prime_idx, cS, paramS);
fprintf('   年龄-储蓄率剖面图已生成。\n');