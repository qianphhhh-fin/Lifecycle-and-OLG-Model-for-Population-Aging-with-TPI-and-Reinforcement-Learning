% =========================================================================
% == SCRIPT: VFI Solver Debugger (Zero Asset Test)
% == 目的：专门测试当期初资产为零时，家庭的最优储蓄决策。
% ==       用于诊断模型是否存在“贫困陷阱”或无法从零资本启动的问题。
% =========================================================================
clear; close all;
addpath(pwd);
fprintf('=== VFI 求解器专项调试脚本 (零资产测试) ===\n\n');

%% 1. 初始化一个【固定】且【合理】的测试环境
fprintf('--- 1. 初始化一个固定的、经济合理的测试环境 ---\n');
cS = main_olg_v14_utils.ParameterValues_HuggettStyle();
paramS = struct();
% --- 使用与主程序一致的、可能导致问题的参数 ---
cS.beta = 0.97; % [关键] 使用可能导致问题的原始 beta 值
cS.sigma = 2.0; 
cS.cFloor = 0.05; 
cS.nSim = 100; cS.time_Step = 5;
cS.alpha = 0.5; cS.ddk = 1 - (1 - 0.05)^cS.time_Step;
cS.tau_k = 0.0; cS.tau_l = 0.0; cS.tau_c = 0.0; 
cS.gov_exp_frac_Y = 0.15; cS.rho_prime_payg = 0.5;
ngrid = 20; cS.nk = ngrid; cS.nkpps = 1; cS.nkprime = ngrid; cS.npps = 1;
cS.A = 1.0; cS = main_olg_v14_utils.generateGrids(cS);
[Z_path, ~, ~] = main_olg_v14_utils.load_exogenous_paths(cS);
cS.theta_path = zeros(size(Z_path, 2), 1) + 0.10;

[paramS.leLogGridV, paramS.leTrProbM, paramS.leProb1V] = main_olg_v14_utils.EarningProcess_olgm(cS);
paramS.leGridV = exp(paramS.leLogGridV(:));

M_age = struct();
M_age.w_t = 1.5;
M_age.r_mkt_period = 0.04 * cS.time_Step; % 约4%年化市场回报
M_age.r_net_period = M_age.r_mkt_period * (1 - cS.tau_k);
M_age.b_t = 0.5;
M_age.current_t = 1;
cS.theta_t = 0.10;
cS.pps_active = false;
fprintf('✅ 测试环境初始化完成。\n\n');


%% 2. 模拟VFI过程，聚焦一个年龄组
a_check = 5; % 检查一个中年工作年龄组 (a=5)
fprintf('--- 2. 检查年龄组 a = %d 的决策过程 ---\n', a_check);

fprintf('   假设下一期(a=%d)的价值函数 V'' 是 log(k)。\n', a_check + 1);
vPrime_kkppse_next = log(cS.kGridV + 1e-6)'; % 加一个小常数避免log(0)
vPrime_kkppse_next = repmat(vPrime_kkppse_next, [1, 1, cS.nw]);
vPrime_kkppse_next = reshape(vPrime_kkppse_next, [cS.nk, 1, cS.nw]);

fprintf('   调用 HHSolutionByAge_VFI_GridSearch...\n');
[~, kPol, ~, ~] = main_olg_v14_utils.HHSolutionByAge_VFI_GridSearch(...
    a_check, vPrime_kkppse_next, M_age, 0, 0, paramS, cS);


%% 3. [核心测试] 手动剖析【零资产】状态点的计算细节
fprintf('\n--- 3. 手动剖析【零资产】状态点的计算细节 ---\n');
ik_check = 1; % *** 固定在第一个资产点 (k=0) ***
ie_check = round(cS.nw / 2);

k_state = cS.kGridV(ik_check);
epsilon_state = paramS.leGridV(ie_check);
fprintf('   检查状态点: k = %.2f, ε = %.2f\n', k_state, epsilon_state);

% --- a. 计算预算约束 ---
labor_income_gross = M_age.w_t * cS.ageEffV_new(a_check) * epsilon_state;
payg_tax = cS.theta_t * labor_income_gross;
labor_tax = cS.tau_l * max(0, labor_income_gross - payg_tax);
k_return = k_state * (1 + M_age.r_net_period); % k=0时，此项为0
net_cash = k_return + labor_income_gross - payg_tax - labor_tax;

fprintf('   - 劳动收入 (毛): %.4f\n', labor_income_gross);
fprintf('   - 可用净现金 (C+k''): %.4f\n', net_cash);

% --- b. 遍历几个储蓄选项，手动计算其价值 ---
fprintf('\n   手动计算几个储蓄选项 k'' 的总价值:\n');
fprintf('%10s | %10s | %10s | %10s | %10s | %12s\n', 'k''', 'C', 'U(C)', 'EV', 'beta*S*EV', 'Total V');
fprintf('%s\n', repmat('-', 68, 1));

EV_matrix = zeros(cS.nk, 1, cS.nw);
for ie_current = 1:cS.nw
    transition_probs = paramS.leTrProbM(ie_current, :);
    vPrime_reshaped = reshape(vPrime_kkppse_next, [cS.nk, cS.nw]);
    EV_slice = vPrime_reshaped * transition_probs';
    EV_matrix(:, 1, ie_current) = EV_slice;
end
ev_interpolant = griddedInterpolant(cS.kGridV, squeeze(EV_matrix(:,1,ie_check)), 'linear');

beta_5yr = cS.beta ^ cS.time_Step;
s_5yr = cS.s_pathV(a_check);
effective_discount_factor = beta_5yr * s_5yr;

% 遍历从不储蓄到全部储蓄的选项
k_prime_options = linspace(cS.kMin, net_cash - cS.cFloor*(1+cS.tau_c), 5);
for k_prime_choice = k_prime_options
    if k_prime_choice < cS.kMin, continue; end
    c_expend = net_cash - k_prime_choice;
    c_choice = c_expend / (1 + cS.tau_c);
    if c_choice < cS.cFloor, continue; end
    
    [~, util] = main_olg_v14_utils.CES_utility(c_choice, cS.sigma, cS);
    k_prime_clamped = max(cS.kGridV(1), min(cS.kGridV(end), k_prime_choice));
    ev = ev_interpolant(k_prime_clamped);
    total_val = util + effective_discount_factor * ev;
    
    fprintf('%10.2f | %10.2f | %10.2f | %10.2f | %10.2f | %12.4f\n', ...
        k_prime_choice, c_choice, util, ev, effective_discount_factor * ev, total_val);
end

k_prime_from_solver = kPol(ik_check, 1, ie_check);
fprintf('\n   函数为 k=0 计算出的最优决策: k'' = %.4f\n', k_prime_from_solver);

if k_prime_from_solver > cS.kMin
    fprintf('✅ 诊断: 模型具有自启动能力。在k=0时，家庭依然选择正储蓄。\n');
    fprintf('   稳态崩溃问题更可能与高层聚合或迭代过程中的极端价格有关。\n');
else
    fprintf('❌ 诊断: 模型存在“贫困陷阱”。在k=0时，最优决策为零储蓄。\n');
    fprintf('   需要增强储蓄动机 (提高beta, 提高sigma) 才能让经济启动。\n');
end
fprintf('--- 调试结束 ---\n');