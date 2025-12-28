% =========================================================================
% == 脚本: solve_and_verify_ss.m
% == 目的: 作为一个完全自包含的单元，加载初始参数后，从头开始
% ==         求解并立即验证终期稳态。此方法确保求解过程和验证过程
% ==         使用完全相同的代码路径，杜绝任何因外部调用或状态污染
% ==         导致的不一致。
% == 版本 v1.0: 终极解决方案。
% =========================================================================
clear; clc;
fprintf('--- 启动稳态求解与验证一体化脚本 (v1.0) ---\n\n');

% --- 1. 加载求解器找到的均衡结果和参数 ---
try
    load('SS/data_for_transition_nopps.mat', 'data_for_transition');
    % 我们只使用加载的 cS 和 paramS 作为起点
    cS_init = data_for_transition.cS;
    paramS_init = data_for_transition.paramSF;
    Z_ss_norm = data_for_transition.cS.Z_path(:, end);
    fprintf('✅ 成功加载初始参数 cS 和 paramSF。\n');
catch ME
    error('无法加载 .mat 文件。请确保已运行最新的求解器并生成文件。\n错误信息: %s', ME.message);
end

% --- 2. 设置求解器 ---
fprintf('\n--- [第1步] 正在启动自包含的稳态求解器 ---\n');

% 定义匿名函数，将 cS_init, paramS_init, Z_ss_norm "锁定" 到误差函数中
% [重要] 使用 v49 版本的 main_steady_state_utils_bgp.m
system_wrapper = @(x) main_steady_state_utils_bgp.system_of_equations_price_based(x, Z_ss_norm, cS_init, paramS_init);

% 设定初始猜测值和边界
x0 = [0.04, 0.4, 1.5, 0.02]; 
lb = [1e-4, 1e-4, 1e-4, 1e-8];
ub = [10, 2.0, 10.0, 10.0];

options = optimoptions('lsqnonlin', 'Display', 'iter', 'FunctionTolerance', 1e-9, 'StepTolerance', 1e-9, 'MaxIterations', 200, 'Algorithm', 'trust-region-reflective');
[x_eq, resnorm, ~, exitflag] = lsqnonlin(system_wrapper, x0, lb, ub, options);

if exitflag <= 0 || resnorm > 1e-9
    error('求解器未能在此脚本中找到均衡解。模型可能存在根本性问题。');
end
fprintf('   ✅ 求解器在本脚本内部成功收敛！\n');

% --- 3. 用求出的解，独立计算最终的 ss 结构体 ---
fprintf('\n--- [第2步] 正在使用求出的均衡解计算最终宏观变量 ---\n');
r_eq = x_eq(1); L_eq = x_eq(2); Kg_eq = x_eq(3); Beq_eq = x_eq(4);
Kp_eq = main_steady_state_utils_bgp.get_Kp_from_r(r_eq, Kg_eq, L_eq, cS_init);

% 调用核心引擎，得到最终的、应是会计闭合的 ss 结构体
[ss_final, ~, ~, ~] = main_steady_state_utils_bgp.calculate_aggregates_unified(Kp_eq, Kg_eq, L_eq, Beq_eq, Z_ss_norm, cS_init, paramS_init);

if isempty(ss_final)
    error('在均衡解处计算最终聚合变量失败。');
end
fprintf('   ✅ 最终 ss 结构体已生成。\n');


% --- 4. 最终国民账户检验 ---
fprintf('\n\n--- [最终检验] 基于自包含求解结果的国民经济核算 ---\n');
g_A_period = (1 + cS_init.g_A_ss)^cS_init.time_Step - 1;
n_period = (1 + cS_init.n_ss)^cS_init.time_Step - 1;
g_total_period = (1 + g_A_period) * (1 + n_period) - 1;

Ip_demand = (g_total_period + cS_init.ddk) * ss_final.K_private_begin_hat;
Ig_demand = (g_total_period + cS_init.ddk_g) * ss_final.K_public_hat;

fprintf('[检验1] 私人资本流量市场 (S_p vs I_p)\n');
fprintf('  私人储蓄供给 (S_p) ......... : %.6f\n', ss_final.Saving_private_flow_Gross);
fprintf('  私人投资需求 (I_p) ................ : %.6f\n', Ip_demand);
fprintf('  --> 流量市场缺口 .................. : %.6e\n\n', ss_final.Saving_private_flow_Gross - Ip_demand);

fprintf('[检验2] 国民收入恒等式 (Y vs C+I+G)\n');
Total_Expenditure = ss_final.C_agg + Ip_demand + Ig_demand + ss_final.G_c;
fprintf('  总产出 (Y) .................. : %.6f\n', ss_final.Y_from_production_hat);
fprintf('  总支出 (C+I_p+I_g+G_c) ...... : %.6f\n', Total_Expenditure);
fprintf('  --> NIPA 缺口 ..................... : %.6e\n\n', ss_final.Y_from_production_hat - Total_Expenditure);

% 附加检验 PAYG 平衡
payg_gap = ss_final.Pension_in_agg - ss_final.Pension_out_agg;
fprintf('[检验3] 养老金系统 (In vs Out)\n');
fprintf('  --> PAYG 系统缺口 ................. : %.6e\n\n', payg_gap);


if abs(ss_final.Saving_private_flow_Gross - Ip_demand) < 1e-6 && abs(ss_final.Y_from_production_hat - Total_Expenditure) < 1e-6 && abs(payg_gap) < 1e-6
    fprintf('====================================================\n');
    fprintf('=== 结论: 模型会计一致性验证通过！所有流量均已闭合。 ===\n');
    fprintf('====================================================\n');
    fprintf('\n最终均衡解:\n');
    disp(ss_final);
else
    fprintf('=========================================================================\n');
    fprintf('=== 结论: 最终失败。即使在自包含脚本中，会计依然不闭合。      ===\n');
    fprintf('=== 这表明问题出在更深层次，很可能在 VFI 或分布计算的函数中。 ===\n');
    fprintf('=========================================================================\n');
end

fprintf('--- 诊断脚本执行完毕 ---\n');