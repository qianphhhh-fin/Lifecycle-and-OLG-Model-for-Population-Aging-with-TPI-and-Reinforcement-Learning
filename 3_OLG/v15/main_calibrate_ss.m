% =========================================================================
% == SCRIPT: main_calibrate_ss.m [vPaper.4 - 支出冲击模型]
% == 目的: 自动校准支出冲击的概率，以匹配一个外部给定的
% ==         宏观目标(稳态资本产出比 K/Y)，然后求解并报告最终结果。
% =========================================================================
clear; close all;
addpath(pwd);
fprintf('=== [支出冲击版] 模型校准与稳态求解 ===\n\n');

%% 1. 初始化环境 (设定不依赖于校准的参数)
fprintf('--- 1. 初始化环境并生成模拟所需输入 ---\n');

cS = main_olg_v15_utils.ParameterValues_ExpenditureShock(); 
paramS = struct();
fprintf('   稳态财政设定：政府预算平衡 (G=T), 意外遗赠被政府100%%征收。\n');
cS.gov_debt_frac_Y = 0;

% 网格设定 (注意：由于pps_active=false，nkpps将由generateGrids自动处理)
ngrid = 100; cS.nk = ngrid; cS.nkpps = ngrid; cS.nkprime = ngrid; cS.npps = 5;
cS = main_olg_v15_utils.generateGrids(cS);

% 注意: 我们在这里暂时不调用 EarningProcess_AgeDependent,
% 因为它将在校准循环内部根据不同的冲击概率被反复调用。
% 我们只需要初始化一些基础值。
paramS.leGridV = []; % 将在循环内生成
paramS.TrProbM_by_age = {}; % 将在循环内生成
paramS.leProb1V = []; % 将在循环内生成
cS.nw_expanded = cS.nw + 2;

% 基于存活率计算稳态年龄人口质量 (这是模型内生的，与校准无关)
age_mass = ones(cS.aD_new, 1);
for ia = 1:(cS.aD_new - 1)
    age_mass(ia+1) = age_mass(ia) * cS.s_pathV(ia);
end
Z_ss_norm = age_mass / sum(age_mass);
fprintf('   环境初始化完成。\n');


%% 2. 校准冲击概率以匹配目标 K/Y
fprintf('\n--- 2. 启动校准程序 ---\n');

KY_TARGET = 3; % <--- 在这里设置您的目标K/Y比率 (例如，中国的实际数据约为3.5-4.0)
fprintf('   目标 K/Y 比率: %.4f\n', KY_TARGET);

% 将不变的参数打包，准备传递给目标函数
% @(p_common) 是一个匿名函数，它只接收一个参数p_common,
% 然后调用我们定义好的校准函数，并传入所有其他固定的参数。
calibration_handle = @(p_common) main_olg_v15_utils.calibrate_ky_target_objective(p_common, cS, paramS, Z_ss_norm, KY_TARGET);

% 设置冲击概率的搜索区间 [p_min, p_max]
% 我们知道p=0.02时K/Y约为2.03，我们需要更高的K/Y，所以解应该在更高的p值域。
% 我们需要找到一个区间，使得区间两端的误差符号相反。
% 例如，p=0.01的K/Y会低于3.5(负误差)，p=0.2的K/Y很可能会高于3.5(正误差)。
p_shock_bracket = [0.01, 0.25]; 

fprintf('   使用fzero在区间 [%.2f, %.2f] 内寻找最优冲击概率...\n', p_shock_bracket(1), p_shock_bracket(2));
options_fzero_outer = optimset('Display', 'iter', 'TolX', 1e-5);
calibration_handle(0.01)
% 调用fzero进行校准！
[p_shock_calibrated, fval, exitflag] = fzero(calibration_handle, p_shock_bracket, options_fzero_outer);

if exitflag <= 0
    error('校准失败！未能找到满足条件的冲击概率。请检查搜索区间或模型参数。');
end

fprintf('\n\n------------------- 校准完成！-------------------\n');
fprintf('✅ 目标 K/Y (%.4f) 已匹配 (最终误差: %.3e)。\n', KY_TARGET, fval);
fprintf('✅ 校准得到的年化冲击概率 p_common = %.6f\n', p_shock_calibrated);
fprintf('--------------------------------------------------\n');


%% 3. 使用校准后的参数，进行最终的、详细的求解与分析
fprintf('\n--- 3. 使用校准后的参数进行最终求解并展示国民账户 ---\n');

% 使用校准得到的最佳参数，更新最终的参数结构体
cS.p_shock_young_peak_annual = p_shock_calibrated;
cS.p_shock_old_peak_annual   = p_shock_calibrated;

% 使用最终参数，重新生成冲击过程
[paramS.leGridV, paramS.TrProbM_by_age, paramS.leProb1V, cS.nw_expanded] = ...
    main_olg_v15_utils.EarningProcess_AgeDependent(cS);
fprintf('   已生成最终的、基于校准后概率的冲击过程。\n');


% 调用您原有的稳态求解器，这次可以打开其内部的迭代显示
% 注意：这里的求解过程会比校准中的单次迭代要长，因为它会执行松弛法精炼
[ss, eq_found, Dist, k_prime_idx] = main_olg_v15_utils.solve_steady_state_iter_K(Z_ss_norm, cS, paramS);
if ~eq_found, error('使用校准参数进行最终求解时失败'); end


%% 4. 最终结果验证与展示
fprintf('\n\n--- 4. 最终结果验证与展示 ---\n');
fprintf('✅✅✅ 校准后模型求解成功！均衡资本 K = %.8f ✅✅✅\n', ss.K_physical);

% 展示最终的国民经济核算，并验证K/Y是否与目标一致
main_olg_v15_utils.display_national_accounts_expenditure_shock(ss, Dist, k_prime_idx, cS, paramS);