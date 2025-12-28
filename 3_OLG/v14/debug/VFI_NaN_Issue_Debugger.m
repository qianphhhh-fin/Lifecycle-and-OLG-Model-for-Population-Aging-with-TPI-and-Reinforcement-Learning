% =========================================================================
% == SCRIPT: VFI_NaN_Issue_Debugger_v2.m
% == 目的：修正v1版本，通过从网格上限开始，【必然】复现 VFI 中的 NaN 问题
% =========================================================================
clear; clc; close all;
addpath(pwd);

fprintf('--- VFI NaN 问题验证脚本 (v2 - 修正版) ---\n');

%% 1. 构建一个最小化的模拟环境 (与v1相同)
fprintf('\n--- 1. 构建最小化环境 ---\n');
cS = main_olg_v14_utils.ParameterValues_HuggettStyle();
paramS = struct();

ngrid = 20; cS.nk = ngrid; cS.nkpps = ngrid; cS.nkprime = ngrid; cS.npps = 3;
cS = main_olg_v14_utils.generateGrids(cS);

M_age = struct('w_t', 1.5, 'r_mkt_period', 0.04, 'r_net_period', 0.04);
cS.pps_active = true;
cS.theta_t = 0.11;
cS.pps_return_rate_premium = 0.01;

[paramS.leLogGridV, paramS.leTrProbM, paramS.leProb1V] = main_olg_v14_utils.EarningProcess_olgm(cS);
paramS.leGridV = exp(paramS.leLogGridV(:));

fprintf('伪造一个下一期的价值函数 V''...\n');
vPrime_fake = randn(cS.nk, cS.nkpps, cS.nw);
vPrime_interpolants = cell(cS.nw, 1);
for ie = 1:cS.nw
    vPrime_interpolants{ie} = griddedInterpolant({cS.kGridV, cS.kppsGridV}, squeeze(vPrime_fake(:,:,ie)), 'linear', 'none');
end
ev_interpolant = vPrime_interpolants{3};

fprintf('✅ 环境构建完毕。\n');

%% 2. 模拟一个家庭的决策，复现问题 (核心修正点)
fprintf('\n--- 2. 复现问题：k_pps_prime 超出网格范围 ---\n');

% --- 【核心修正】---
% 不再选择一个中间值，而是直接选择网格上的【最大值】作为当前状态。
% 这将保证 k_pps_prime 必然超出网格范围。
k_pps_state = cS.kppsGridV(end); 
c_pps_choice = 0.01; % 即使是一个非常小的缴费

pps_return_factor = 1 + M_age.r_mkt_period + cS.pps_return_rate_premium;
k_pps_prime = (k_pps_state + c_pps_choice) * pps_return_factor;

fprintf('家庭当前PPS资产 (k_pps_state): %.4f (选择了网格最大值)\n', k_pps_state);
fprintf('本期PPS缴费 (c_pps_choice):   %.4f\n', c_pps_choice);
fprintf('PPS账户回报因子:             %.4f\n', pps_return_factor);
fprintf('-----------------------------------------------------\n');
fprintf('计算出的下一期PPS资产 (k_pps_prime): %.4f\n', k_pps_prime);
fprintf('PPS资产网格上限 (cS.kppsMax):      %.4f\n', cS.kppsMax);

if k_pps_prime > cS.kppsMax
    fprintf('\n>>> 结论：计算出的 k_pps_prime (%.4f) 显著大于网格上限 (%.4f)！\n', k_pps_prime, cS.kppsMax);
else
    fprintf('\n>>> 未能复现问题。这不应该发生，请检查 generateGrids 和这里的逻辑。\n');
end

%% 3. 演示NaN的产生和传播 (与v1相同)
fprintf('\n--- 3. 演示NaN的产生：用超范围的k_pps_prime查询EV ---\n');

% 假设家庭选择了某个 k_prime_choice
k_prime_choice = cS.kGridV(end-1); % 同样选择一个较大的值

ev_problematic = ev_interpolant(k_prime_choice, k_pps_prime);

fprintf('用 (k_prime=%.2f, k_pps_prime=%.2f) 查询插值器...\n', k_prime_choice, k_pps_prime);
fprintf('>>> 得到的期望价值 EV: %f\n', ev_problematic);

if isnan(ev_problematic)
    fprintf('>>> 验证成功！由于查询点越界，EV返回了 NaN。\n');
    fprintf('    这个 NaN 会被您的代码转为-1e10，导致储蓄决策崩溃。\n');
else
    fprintf('>>> 未能复现NaN，请检查代码或参数。\n');
end

%% 4. 演示修复方案 (与v1相同)
fprintf('\n--- 4. 演示修复方案：对 k_pps_prime 进行边界约束 ---\n');

k_pps_prime_fixed = min(k_pps_prime, cS.kppsMax);

fprintf('原始 k_pps_prime: %.4f\n', k_pps_prime);
fprintf('修复后 (clamped) 的 k_pps_prime: %.4f\n', k_pps_prime_fixed);

ev_fixed = ev_interpolant(k_prime_choice, k_pps_prime_fixed);

fprintf('用修复后的 (k_prime=%.2f, k_pps_prime=%.2f) 查询插值器...\n', k_prime_choice, k_pps_prime_fixed);
fprintf('>>> 得到的新期望价值 EV: %f\n', ev_fixed);

if ~isnan(ev_fixed)
    fprintf('>>> 修复成功！EV现在是一个正常的数值，家庭可以进行有效的储蓄决策。\n');
else
    fprintf('>>> 修复失败，依然得到NaN，这不应该发生，请检查逻辑。\n');
end