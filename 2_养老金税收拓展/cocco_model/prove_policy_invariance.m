%% 证明最优策略在选择变量对数变换下的不变性
%
% 目标:
%   通过一个数值例子，证明对于一个标准的动态规划问题，
%   在原始空间 (C, alpha) 中求解，与在变换空间 (c=log(C), alpha) 中求解，
%   会得到完全相同的最优策略。
%   这验证了方法三(对数变换)的理论正确性。
%
% 作者: Gemini
% 日期: 2025-09-17

% --- 初始化 ---
clear;
close all;
clc;

fprintf('===== 证明选择变量变换下策略不变性的脚本开始 =====\n\n');

%% =====================================================================
%                       1. 设置参数和网格
%  =====================================================================
fprintf('--- 步骤 1: 设置参数和网格 ---\n');
cS.beta = 0.95; cS.gamma = 3.84; cS.rf = 1.02; cS.mu = 0.04; cS.sigr = 0.27;
cS.n_shocks = 5; cS.ncash = 51; cS.maxcash = 100; cS.mincash = 0.25; cS.C_min = 0.01;
cS.survprob_T_minus_1 = 0.68424;
lgcash = linspace(log(cS.mincash), log(cS.maxcash), cS.ncash)';
cS.gcash = exp(lgcash);
[grid] = discretizeAR1_Tauchen(0, 0, 1, cS.n_shocks, 2);
shock_pdf = normpdf(grid); cS.shock_weig = shock_pdf / sum(shock_pdf);
cS.gret = cS.rf + cS.mu + grid * cS.sigr;
fprintf('参数设置完成。\n\n');

%% =====================================================================
%                  2. 定义效用函数
%  =====================================================================
fprintf('--- 步骤 2: 定义原始效用函数 (U) ---\n');
U = @(C) (max(C, cS.C_min).^(1 - cS.gamma)) / (1 - cS.gamma);
fprintf('函数定义完成。\n\n');

%% =====================================================================
%                  3. 计算终端期价值函数
%  =====================================================================
fprintf('--- 步骤 3: 计算终端期 T 的价值函数 V_T ---\n');
C_T = cS.gcash;
V_T = U(C_T);
V_T_interp_obj = griddedInterpolant(cS.gcash, V_T, 'pchip', 'linear');
fprintf('终端期价值函数计算和插值完成。\n\n');

%% =====================================================================
%              4. 在 t=T-1 期进行反向迭代 (核心修改)
%  =====================================================================
fprintf('--- 步骤 4: 在 T-1 期求解贝尔曼方程 ---\n');

% 初始化策略数组
C_policy_A = zeros(cS.ncash, 1); A_policy_A = zeros(cS.ncash, 1);
C_policy_C = zeros(cS.ncash, 1); A_policy_C = zeros(cS.ncash, 1);

options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'interior-point', ...
    'OptimalityTolerance', 1e-18, 'StepTolerance', 1e-18);

fprintf('正在对状态空间上的每个点进行优化...\n');
tic;
for i_cash = 1:cS.ncash
    cash = cS.gcash(i_cash);
    
    % --- 为系统 A (原始空间: C, alpha) 求解 ---
    x0_A = [cash*0.5, 0.5]; 
    lb_A = [cS.C_min, 0]; 
    ub_A = [cash, 1];
    
    obj_fun_A = @(x) -fun_calculate_total_value(x, cash, U, V_T_interp_obj, cS);
    [policy_A, ~] = fmincon(obj_fun_A, x0_A, [], [], [], [], lb_A, ub_A, [], options);
    C_policy_A(i_cash) = policy_A(1);
    A_policy_A(i_cash) = policy_A(2);
    
    % --- 为系统 C (对数变换空间: c=log(C), alpha) 求解 ---
    % 定义新的选择变量 z = [c, alpha]
    % 目标函数：输入z，内部计算 C = exp(z(1))，然后计算价值
    obj_fun_C = @(z) -fun_calculate_total_value([exp(z(1)), z(2)], cash, U, V_T_interp_obj, cS);
    
    % 变换约束边界
    % C >= C_min  =>  exp(c) >= C_min  =>  c >= log(C_min)
    % C <= cash   =>  exp(c) <= cash   =>  c <= log(cash)
    z0_C = [log(cash*0.5), 0.5];
    lb_C = [log(cS.C_min), 0];
    ub_C = [log(cash), 1];
    
    [policy_C_raw, ~] = fmincon(obj_fun_C, z0_C, [], [], [], [], lb_C, ub_C, [], options);
    
    % **核心验证步骤**: 将最优的 c* 转换回 C*
    C_policy_C(i_cash) = exp(policy_C_raw(1));
    A_policy_C(i_cash) = policy_C_raw(2);
end
fprintf('优化完成，耗时 %.2f 秒。\n\n', toc);

%% =====================================================================
%                       5. 比较结果并可视化
%  =====================================================================
fprintf('--- 步骤 5: 比较两个系统的最优策略 ---\n');
max_diff_C = max(abs(C_policy_A - C_policy_C));
max_diff_A = max(abs(A_policy_A - A_policy_C));

fprintf('消费策略 C*(s) 的最大绝对差值: %e\n', max_diff_C);
fprintf('投资策略 a*(s) 的最大绝对差值: %e\n', max_diff_A);

if max_diff_C < 1e-6 && max_diff_A < 1e-6
    fprintf('\n结论: 两个系统计算出的最优策略在数值精度内完全相同。\n');
    fprintf('这正确地证明了，最优策略对于选择变量的对数变换是不变的。\n');
else
    fprintf('\n警告: 策略存在显著差异，请检查代码或优化器容差。\n');
end

figure('Position', [100, 100, 1200, 600], 'Name', '策略不变性数值证明 (对数变换)');
subplot(1, 2, 1);
plot(cS.gcash, C_policy_A, 'b-', 'LineWidth', 2.5, 'DisplayName', '系统 A: max U(C)'); hold on;
plot(cS.gcash, C_policy_C, 'r--','LineWidth', 1.5, 'DisplayName', '系统 C: max U(exp(c))');
grid on; box on; title('最优消费策略 C*(s) 对比'); xlabel('期初现金持有量 (s)'); ylabel('最优消费 C*');
legend('show', 'Location', 'northwest'); xlim([cS.mincash, cS.maxcash]);
subplot(1, 2, 2);
plot(cS.gcash, A_policy_A, 'b-', 'LineWidth', 2.5, 'DisplayName', '系统 A: max U(C)'); hold on;
plot(cS.gcash, A_policy_C, 'r--','LineWidth', 1.5, 'DisplayName', '系统 C: max U(exp(c))');
grid on; box on; title('最优投资组合策略 a*(s) 对比'); xlabel('期初现金持有量 (s)'); ylabel('最优风险资产比例 a*');
legend('show', 'Location', 'best'); xlim([cS.mincash, cS.maxcash]);
sgtitle('策略不变性证明: 对选择变量进行对数变换', 'FontSize', 16, 'FontWeight', 'bold');

fprintf('\n===== 脚本执行完毕 =====\n');

%% =====================================================================
%                       辅助函数
%  =====================================================================

function total_value = fun_calculate_total_value(x, cash, reward_func, V_next_interp_obj, cS)
    % 这个函数只计算贝尔曼方程的右侧（即Q值），不取负
    % 输入 x = [C, alpha]
    consumption = x(1);
    alpha = x(2);
    current_reward = reward_func(consumption);
    sav = cash - consumption;
    cash_next = (cS.rf * (1 - alpha) + cS.gret' * alpha) * sav;
    V_next = V_next_interp_obj(cash_next);
    expected_V_next = cS.shock_weig' * V_next';
    total_value = current_reward + cS.beta * cS.survprob_T_minus_1 * expected_V_next;
end

function [z_grid] = discretizeAR1_Tauchen(mew, rho, sigma, znum, Tauchen_q)
    if znum == 1, z_grid = mew / (1 - rho); return; end
    zstar = mew / (1 - rho);
    sigmaz = sigma / sqrt(1 - rho^2);
    z_grid = zstar + linspace(-Tauchen_q * sigmaz, Tauchen_q * sigmaz, znum)';
end