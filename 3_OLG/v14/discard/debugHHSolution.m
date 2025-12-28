% --- START OF FILE debugHHSolution.m (Naive Policy Test Version) ---

% =========================================================================
% == SCRIPT: Unit Test for Naive Policy and HHSimulation
% =========================================================================
clear; close all;
addpath(pwd);
fprintf('=== 单元测试: Naive 策略下的会计平衡 ===\n\n');

%% 1. 设置环境 (与之前相同)
fprintf('--- 1. 初始化测试环境...\n');
cS = main_olg_v14_utils.ParameterValues_HuggettStyle();
cS.pps_active=true;
paramS = struct();
cS.sigma = 3.0; cS.beta = 0.97; cS.cFloor = 0.05; cS.nSim = 2000;
cS.phi_bequest = 0; cS.sigma_bequest = cS.sigma;
cS.start_year = 1997; cS.end_year = 2102; cS.time_step = 5;
cS.alpha = 0.5; cS.ddk = 1 - (1 - 0.05)^cS.time_step;
cS.tau_k = 0.2; cS.tau_l = 0.1; cS.tau_c = 0.05; cS.gov_exp_frac_Y = 0.15;
cS.rho_prime_payg = 0.5; cS.pps_activation_year = 2022;
cS.pps_tax_rate_withdrawal = 0.03; cS.pps_return_rate_premium = 0.01;
cS.pps_withdrawal_rate = 0.15; cS.pps_contrib_limit = 9999; cS.pps_max_contrib_frac = 0.1;
ngrid = 10; cS.nk = ngrid; cS.nkpps = ngrid; cS.nkprime = ngrid; cS.npps = 5;
cS.A = 1.0; cS = main_olg_v14_utils.generateGrids(cS);
cS.theta_t = 0.08; 
[paramS.leLogGridV, paramS.leTrProbM, paramS.leProb1V] = main_olg_v14_utils.EarningProcess_olgm(cS);
paramS.leGridV = exp(paramS.leLogGridV(:));
paramS.ageEffV_new = cS.ageEffV_new;
cS.aR_new = 9;
M_age = struct('w_t', 1.2, 'r_net_period', 0.04, 'r_mkt_period', 0.05, 'b_t', cS.rho_prime_payg * 1.2);
tr_per_capita_age = 0.02;
b_age_val_retiree = cS.rho_prime_payg * M_age.w_t;
% 为了调用模拟器，需要一个完整的 eIdxM
eIdxM_group = main_olg_v14_utils.LaborEndowSimulation_olgm_AgeGroup(cS, paramS);

fprintf('✅ 环境初始化完成。\n\n');

%% 2. [核心测试] 调用 Naive 策略并模拟
fprintf('--- 2. 调用 Naive 策略并使用 HHSimulation_olgm 模拟...\n');
% 确保 HHSolution_VFI_Huggett 中的开关设置为 true
[cPolM, kPolM, cPpsPolM, ~] = main_olg_v14_utils.HHSolution_VFI_Huggett(M_age, tr_per_capita_age, paramS, cS);

% 使用 HHSimulation_olgm 来反推微观历史
% 注意：这里的 cPolM 是 Naive 策略生成的，我们忽略它，只关心反推的结果
[kHistM, kPpsHistM, cHistM, cppsHistM] = main_olg_v14_utils.HHSimulation_olgm(kPolM, cPpsPolM, cPolM, eIdxM_group, M_age, paramS, cS, tr_per_capita_age);

fprintf('✅ 调用和模拟完成。\n\n');

%% 3. [审计] 手动、严格地核算现金流
fprintf('--- 3. 手动进行现金流审计 (基于 HHSimulation_olgm 的结果)...\n');
% 随机抽取一个体进行审计
iSim = 10; 

% --- 审计工人 (a=9) ---
a_worker = 9;
fprintf('--- 审计工作期 (a=%d), 个体 %d...\n', a_worker, iSim);
k_state_w = kHistM(iSim, a_worker);
epsilon_state_w = paramS.leGridV(eIdxM_group(iSim, a_worker));
% 来自模拟器的结果
c_worker_sim = cHistM(iSim, a_worker);
k_prime_worker_sim = kHistM(iSim, a_worker+1);
cpps_worker_sim = cppsHistM(iSim, a_worker);

% 手动计算流入
in_k_return_w = k_state_w * (1 + M_age.r_net_period);
in_labor_w = M_age.w_t * cS.ageEffV_new(a_worker) * epsilon_state_w;
in_bequest_w = tr_per_capita_age; % a=9 > 1
total_inflow_w = in_k_return_w + in_labor_w + in_bequest_w;

% 手动计算流出
out_consump_w = c_worker_sim * (1 + cS.tau_c);
out_k_prime_w = k_prime_worker_sim;
out_payg_tax_w = cS.theta_t * in_labor_w;
out_labor_tax_w = cS.tau_l * max(0, in_labor_w - cpps_worker_sim);
out_cpps_w = cpps_worker_sim;
total_outflow_w = out_consump_w + out_k_prime_w + out_payg_tax_w + out_labor_tax_w + out_cpps_w;
residual_worker = total_inflow_w - total_outflow_w;
fprintf('%-25s | %15.4e\n', '工作期现金流残差', residual_worker);

% --- 审计退休人员 (a=10) ---
a_retiree = 10;
fprintf('--- 审计退休期 (a=%d), 个体 %d...\n', a_retiree, iSim);
k_state_r = kHistM(iSim, a_retiree);
k_pps_state_r = kPpsHistM(iSim, a_retiree);
% 来自模拟器的结果
c_retiree_sim = cHistM(iSim, a_retiree);
k_prime_retiree_sim = kHistM(iSim, a_retiree+1);

% 手动计算流入
in_k_return_r = k_state_r * (1 + M_age.r_net_period);
in_pension_r = b_age_val_retiree;
in_bequest_r = tr_per_capita_age;
pps_w_gross_r = k_pps_state_r * cS.pps_withdrawal_rate;
in_pps_w_net_r = pps_w_gross_r * (1 - cS.pps_tax_rate_withdrawal);
total_inflow_r = in_k_return_r + in_pension_r + in_bequest_r + in_pps_w_net_r;

% 手动计算流出
out_consump_r = c_retiree_sim * (1 + cS.tau_c);
out_k_prime_r = k_prime_retiree_sim;
total_outflow_r = out_consump_r + out_k_prime_r;
residual_retiree = total_inflow_r - total_outflow_r;
fprintf('%-25s | %15.4e\n', '退休期现金流残差', residual_retiree);
fprintf('\n');

%% 4. 最终结论
if abs(residual_worker) < 1e-9 && abs(residual_retiree) < 1e-9
    fprintf('✅ 单元测试通过！HHSimulation_olgm 正确执行了 Naive 策略的预算约束！\n');
else
    fprintf('❌ 单元测试失败！HHSimulation_olgm 的会计逻辑仍然存在问题。\n');
end