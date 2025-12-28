% ===================================================================
% test_solver_comparison.m - 比较不同利率下储蓄策略的测试脚本
% ===================================================================
clc; clear all; close all;
fprintf('--- 开始VFI求解器对比测试：利率 vs. 储蓄策略 ---\n');

%% 1. 初始化参数和环境
% (完全采用您提供的参数集)
fprintf('\n--- 1. 初始化参数 ---\n');
T = 150; 
cS = main_olg_v10_utils.ParameterValues_HuggettStyle();
paramS = struct();

% --- V10 固定的政策参数 ---
cS.rho_prime_payg = 0.4; 
cS.theta_payg = 0.15; 
B_p_Y_ratio_target = 0.05;

% --- 家庭偏好参数 ---
cS.sigma = 3.0;
cS.beta = 0.97;
cS.cFloor = 0.05;
cS.phi_bequest = 0; 
cS.sigma_bequest = cS.sigma;

% --- 生产技术参数 ---
cS.A = 2;
cS.alpha = 0.36;
cS.ddk = 1 - (1 - 0.08)^5;  % 约等于 0.341

% --- 政府财政参数 ---
cS.tau_k = 0.20;
cS.tau_c = 0.02;
cS.tau_l = 0.10;
cS.gov_exp_frac_Y = 0.15;
cS.gov_debt_frac_Y = 0.60;

% --- 财政反馈规则参数 ---
cS.debt_ratio_target = 0.60;
cS.debt_feedback_strength = 0.5;

% --- 网格和模拟参数 ---
ngrid = 50; % 使用一个更精细的网格以便绘图
cS.nk = ngrid;
cS.nkpps = ngrid;
cS.nkprime = ngrid;
cS.npps = ngrid;
cS.nSim = 1000;

% --- 完成参数初始化 ---
cS = main_olg_v10_utils.generateGrids(cS);
[paramS.leLogGridV, paramS.leTrProbM, paramS.leProb1V] = main_olg_v10_utils.EarningProcess_olgm(cS);
paramS.leGridV = exp(paramS.leLogGridV(:));
fprintf('参数加载完成。\n');


%% 2. 定义要测试的利率环境
% 定义一个包含负、零、正利率的列表
r_net_list = [-0.3,-0.1,-0.05, 0, 0.02, 0.04]; % 测试年化 -2%, 0%, 2%, 4%
R_k_net_factor_list = 1 + r_net_list;

% 存储结果
k_policy_results = cell(length(r_net_list), 1);
legend_labels = cell(length(r_net_list), 1);

%% 3. 循环求解并存储策略函数
fprintf('\n--- 2. 循环求解不同利率下的VFI ---\n');

% 定义一个固定的、简单的宏观环境（除了利率）
w_gross = 1.5;
TR_total = 0;
bV_payg = zeros(1, cS.aD_new);
bV_payg(cS.aR_new+1:end) = 0.3; 

for i = 1:length(r_net_list)
    R_k_net_factor_hh = (R_k_net_factor_list(i))^cS.yearStep; % <<<< 转换为5年期回报因子
    r_annual = r_net_list(i);
    
    fprintf('正在求解: 年化利率 r_net = %.2f%% (5年期回报因子 = %.4f)\n', r_annual*100, R_k_net_factor_hh);
    
    % 调用VFI求解器
    [~, kPolM, ~, ~] = main_olg_v10_utils.HHSolution_VFI_Huggett(...
        R_k_net_factor_hh, w_gross, TR_total, bV_payg, ...
        paramS, cS, 'vectorize_grid');
    
    % 存储我们关心的策略函数切片
    % (倒数第二期, 最低PPS资产, 最低效率冲击)
    a_idx_to_plot = cS.aD_new - 1;
    ikpps_to_plot = 1;
    ie_to_plot = 1;
    k_policy_results{i} = squeeze(kPolM(:, ikpps_to_plot, ie_to_plot, a_idx_to_plot));
    
    % 创建图例标签
    legend_labels{i} = sprintf('r = %.1f%%', r_annual * 100);
end

fprintf('所有VFI求解完成。\n');

%% 4. 可视化比较结果
fprintf('\n--- 3. 绘制比较图 ---\n');

figure('Name', '利率对储蓄策略的影响', 'Position', [100, 100, 800, 600]);
hold on;

colors = lines(length(r_net_list)); % 为不同利率生成不同的颜色

for i = 1:length(r_net_list)
    plot(cS.kGridV, k_policy_results{i}, 'LineWidth', 2, 'Color', colors(i,:));
end

hold off;

title(sprintf('不同利率下的储蓄策略函数 k''(k) \n (年龄 a=%d)', a_idx_to_plot), 'FontSize', 14);
xlabel('当前资产 k', 'FontSize', 12);
ylabel('下一期储蓄 k''', 'FontSize', 12);
legend(legend_labels, 'Location', 'northwest', 'FontSize', 10);
grid on;
axis tight;
% 添加一条 y=x 的45度线作为参考
line(xlim, ylim, 'Color', 'k', 'LineStyle', '--');
legend([legend_labels, {'k'' = k (45-degree line)'}]);

fprintf('绘图完成。请查看生成的图形。\n');