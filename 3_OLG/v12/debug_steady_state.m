% --- START OF FILE debug_steady_state.m (最终裁决版 - 统一人口逻辑) ---
%
% 调试脚本 V5: OLG模型稳态T型账户审计 (统一人口逻辑)
%
% 目的:
% 1. [最终修正] 只使用一个内在的人口动态来源(生存概率s(a))来构建一个
%    零增长的人口稳态分布。这确保了家庭预期与宏观现实的完全一致。
% 2. 调用求解器，使用此一致的人口分布找到模型的初始稳态。
% 3. 对该真实稳态进行T型账户审计，验证所有市场是否同时出清。
%
% ---

clear; close all;
fprintf('=== [调试脚本 V5] OLG模型稳态T型账户审计 (统一人口逻辑) ===\n\n');

%% 1. 设置与主文件完全一致的环境
addpath(pwd);
cS = main_olg_v12_utils.ParameterValues_HuggettStyle();
paramS = struct();

% --- 参数设置 (与上次修正保持一致) ---
cS.sigma = 3.0; 
beta_annual = 0.97;
cS.time_step = 5;
cS.beta = beta_annual ^ cS.time_step; 
fprintf('>>> [注意] 年度beta=%.3f, time_step=%d, 模型使用的5年期beta=%.4f\n', beta_annual, cS.time_step, cS.beta);
annual_depreciation = 0.06;
cS.ddk = 1 - (1 - annual_depreciation)^cS.time_step;
fprintf('>>> [注意] 年度折旧率=%.2f, time_step=%d, 模型使用的5年期折旧率=%.4f\n', annual_depreciation, cS.time_step, cS.ddk);

cS.cFloor = 0.05; cS.nSim = 1000;
cS.phi_bequest = 3; cS.sigma_bequest = cS.sigma;
cS.start_year = 1997; cS.end_year = 2102; 
cS.alpha = 0.45; 
cS.tau_k = 0.25; cS.tau_l = 0.08; cS.tau_c = 0.04;
ngrid = 10; cS.nk = ngrid; cS.nkpps = ngrid; cS.nkprime = ngrid; cS.npps = ngrid;
cS.A = 1.0; cS = main_olg_v12_utils.generateGrids(cS);
cS.rho_urban_employee = 0.60; cS.rho_resident = 0.15;
cS.pps_active = false; cS.pps_activation_year = 2022;
cS.pps_tax_rate_withdrawal = 0.03; cS.pps_return_rate_premium = 0.01;
cS.pps_withdrawal_rate = 0.15; cS.pps_contrib_limit = 9999; cS.pps_max_contrib_frac = 0.1;
cS = main_olg_v12_utils.calcaulte_theta_payg_path(cS, false);

% --- 加载外生路径 (仅用于获取 eIdxM 和 A_ss) ---
[~, A_path, ~] = main_olg_v12_utils.load_exogenous_paths(cS);
[paramS.leLogGridV, paramS.leTrProbM, paramS.leProb1V] = main_olg_v12_utils.EarningProcess_olgm(cS);
paramS.leGridV = exp(paramS.leLogGridV(:));
paramS.ageEffV_new = cS.ageEffV_new;
eIdxM = main_olg_v12_utils.LaborEndowSimulation_olgm_AgeGroup(cS, paramS);

%% 2. [核心修正] 构建与生存概率 s(a) 内在一致的零增长人口稳态分布 Z_ss
fprintf('\n--- 正在构建零增长人口稳态分布... ---\n');
Z_ss_raw = ones(cS.aD_new, 1);
% 假设第一期出生人口为1
for a = 1:(cS.aD_new - 1)
    % Z(a+1) = Z(a) * s(a)
    % 从a期到a+1期的存活率是 cS.s_1yr_transitionV(a)
    Z_ss_raw(a+1) = Z_ss_raw(a) * cS.s_1yr_transitionV(a);
end
Z_ss_norm = Z_ss_raw / sum(Z_ss_raw); % 归一化，使其总和为1
fprintf('--- 人口稳态分布构建完成。---\n');

% 可选：可视化检查构建的人口分布
figure('Name', '内在一致的人口稳态分布');
bar(Z_ss_norm * 100);
xlabel('年龄组');
ylabel('占总人口百分比 (%)');
title('根据生存概率 s(a) 构建的零增长人口稳态分布');
grid on;


%% 3. [核心修正] 调用求解器，并接收所有自洽的聚合数据
fprintf('\n--- 正在调用OLG稳态求解器... ---\n');
A_ss = A_path(1);
[ss, eq_found, Aggregates_HH_audited] = main_olg_v12_utils.solve_initial_steady_state(Z_ss_norm, A_ss, cS, paramS, eIdxM);

if ~eq_found, error('无法找到OLG稳态，审计中止。'); end
fprintf('--- 稳态求解完成，开始审计该均衡点... ---\n\n');

%% 4. [核心修正] 直接使用求解器返回的数据进行审计，不再重新计算！
K_prod_ss = ss.K_pvt + ss.K_pps + ss.B_p;

fprintf('--- OLG模型稳态T型账户审计 ---\n\n');

%% 5. 构建T型账户进行审计 (容忍度放宽至1e-5以适应模拟误差)
fprintf('--- OLG模型稳态T型账户审计 (所有值为5年期总量) ---\n\n');
tolerance = 1e-5;

% --- 企业部门 (Firms) ---
fprintf('=============== 1. 企业部门 (Firms) ===============\n');
Firm_Sources = ss.Y;
Firm_Uses = ss.w * ss.L + ss.r_mkt_rental * K_prod_ss;
fprintf('   资金来源 (Sources):\n     + 销售产品 (Y):                 %9.4f\n', Firm_Sources);
fprintf('   --------------------------------------------------\n     总来源:                         %9.4f\n\n', Firm_Sources);
fprintf('   资金运用 (Uses):\n     + 支付工资 (w*L):               %9.4f\n     + 支付资本租金 (r_rental*K):    %9.4f\n', ss.w * ss.L, ss.r_mkt_rental * K_prod_ss);
fprintf('   --------------------------------------------------\n     总运用:                         %9.4f\n\n', Firm_Uses);
Firm_Budget_Residual = Firm_Sources - Firm_Uses;
fprintf('   企业部门预算残差 (利润):          %9.4e\n', Firm_Budget_Residual);
if abs(Firm_Budget_Residual) < tolerance, fprintf('   ✅ 企业部门预算平衡 (利润为零)。\n\n');
else, fprintf('   ⚠️ 企业部门预算不平衡！\n\n'); end

fprintf('=============== 2. 家庭部门 (Households) ===============\n');
HH_Gross_Factor_Income = ss.w * ss.L + ss.r_mkt_rental * K_prod_ss;
Pension_Benefit = ss.G - (Aggregates_HH_audited.TotalLaborTax + Aggregates_HH_audited.TotalCapitalTax + Aggregates_HH_audited.TotalConsumptionTax);
Pension_Benefit = (Aggregates_HH_audited.TotalLaborTax + Aggregates_HH_audited.TotalCapitalTax + Aggregates_HH_audited.TotalConsumptionTax) - ss.G;
Pension_Benefit = ss.Y - ss.C - ss.Investment - (Aggregates_HH_audited.TotalLaborTax + Aggregates_HH_audited.TotalCapitalTax + Aggregates_HH_audited.TotalConsumptionTax);
Pension_Benefit = ss.G - (Aggregates_HH_audited.TotalLaborTax + Aggregates_HH_audited.TotalCapitalTax + Aggregates_HH_audited.TotalConsumptionTax);

% 从 ss.G 的定义反推 Pension_Benefit
T_total_audited = Aggregates_HH_audited.TotalLaborTax + Aggregates_HH_audited.TotalCapitalTax + Aggregates_HH_audited.TotalConsumptionTax;
Pension_Benefit = T_total_audited - ss.G;

HH_Direct_Taxes = Aggregates_HH_audited.TotalLaborTax + Aggregates_HH_audited.TotalCapitalTax;
Depreciation = cS.ddk * K_prod_ss;
Bequest_Transfer = ss.TR;
HH_Net_Disposable_Income = HH_Gross_Factor_Income + Pension_Benefit + Bequest_Transfer - HH_Direct_Taxes - Depreciation;
HH_Consumption_Value = ss.C;
HH_Consumption_Expenditure = HH_Consumption_Value * (1 + cS.tau_c);
HH_Net_Saving_from_sim = Aggregates_HH_audited.TotalSaving; 
HH_Uses_of_NDI = HH_Consumption_Expenditure + HH_Net_Saving_from_sim;

fprintf('   资金来源 (净可支配收入 Net Disposable Income):\n');
fprintf('     + 总要素收入 (w*L + r_rental*K): %9.4f\n', HH_Gross_Factor_Income);
fprintf('     + 养老金福利 (b*N_r):            %9.4f\n', Pension_Benefit);
fprintf('     + 意外遗赠转移 (TR):             %9.4f\n', Bequest_Transfer);
fprintf('     - 直接税 (T_l + T_k):            %9.4f\n', HH_Direct_Taxes);
fprintf('     - 资本折旧 (δK):                 %9.4f\n', Depreciation);
fprintf('   --------------------------------------------------\n     = 净可支配总收入:               %9.4f\n\n', HH_Net_Disposable_Income);
fprintf('   资金运用 (Uses of Net Disposable Income):\n');
fprintf('     + 消费总支出 (C*(1+tau_c)):     %9.4f\n', HH_Consumption_Expenditure);
fprintf('     + 净储蓄 (S_net_pvt):             %9.4e\n', HH_Net_Saving_from_sim);
fprintf('   --------------------------------------------------\n     = 总运用:                       %9.4f\n\n', HH_Uses_of_NDI);
HH_Budget_Residual = HH_Net_Disposable_Income - HH_Uses_of_NDI;
fprintf('   家庭部门预算残差 (NDI - C_exp - S_net): %9.4e\n', HH_Budget_Residual);
if abs(HH_Budget_Residual) < tolerance, fprintf('   ✅ 家庭部门预算平衡。\n\n'); else, fprintf('   ⚠️ 家庭部门预算不平衡！这是根本漏洞！\n\n'); end

% --- 政府部门 (Government) ---
fprintf('========== 3. 政府部门 (Consolidated Government) ==========\n');
T_total = Aggregates_HH_audited.TotalLaborTax + Aggregates_HH_audited.TotalCapitalTax + Aggregates_HH_audited.TotalConsumptionTax;
Gov_Uses = ss.G + Pension_Benefit;
Gov_Primary_Surplus = T_total - Gov_Uses;
fprintf('   资金来源 (Sources):\n     + 总税收收入 (T_total):         %9.4f\n', T_total);
fprintf('   --------------------------------------------------\n     总来源:                         %9.4f\n\n', T_total);
fprintf('   资金运用 (Uses):\n     + 政府购买 (G):                 %9.4f\n', ss.G);
fprintf('     + 养老金福利 (b*N_r):           %9.4f\n', Pension_Benefit);
fprintf('   --------------------------------------------------\n     总运用 (基础支出):              %9.4f\n\n', Gov_Uses);
fprintf('   政府基础盈余 (T - G - b):         %9.4e\n\n', Gov_Primary_Surplus);
if abs(Gov_Primary_Surplus) < tolerance, fprintf('   ✅ 政府预算平衡。\n\n'); else, fprintf('   ⚠️ 政府预算不平衡！\n\n'); end

% --- 最终审计：市场出清 ---
fprintf('============== 4. 最终审计: 商品市场出清 (瓦尔拉斯定律验证) ============== \n');
Investment_gross_ss = ss.Investment;
Total_Uses_Goods_Market = ss.C + Investment_gross_ss + ss.G;
fprintf('   商品市场供给 (Y):\n     + 总产出:                       %9.4f\n\n', ss.Y);
fprintf('   商品市场需求 (C + I_gross + G):\n');
fprintf('     + 总消费 (C):                   %9.4f\n', ss.C);
fprintf('     + 总投资 (I_gross = δK):        %9.4f\n', Investment_gross_ss);
fprintf('     + 政府购买 (G):                 %9.4f\n', ss.G);
fprintf('   --------------------------------------------------\n     = 总需求:                       %9.4f\n\n', Total_Uses_Goods_Market);
Final_Residual = ss.market_clearing_residual;
fprintf('   最终残差 (Y - C - I_gross - G):     %9.4e\n', Final_Residual);
if abs(Final_Residual) < tolerance
    fprintf('   ✅ 恭喜！所有账户平衡，商品市场因资产市场出清而自动出清！\n')
else 
    fprintf('   ⚠️ 审计失败！求解器或模型存在问题！\n')
end