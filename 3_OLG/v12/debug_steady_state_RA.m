% --- START OF FILE debug_steady_state_RA.m ---
%
% 调试脚本 RA版: 代表性消费者模型 T型账户审计
%
% 目的:
% 1. 求解一个简化的代表性消费者模型的稳态。
% 2. 用与复杂模型完全相同的会计逻辑审计该模型的T型账户。
% 3. 验证核心预算约束和会计恒等式是否能在最简单的环境下成立。
%
% ---

clear; close all;
fprintf('=== [调试脚本 RA版] 代表性消费者模型 T型账户审计 ===\n\n');

%% 1. 设置参数 (与主文件保持一致)
addpath(pwd);
cS = struct();
% 标准的年化 beta
cS.beta = 0.96; 
cS.sigma = 3.0; 
cS.alpha = 0.45;
% 年化折旧率 8%
cS.ddk = 0.08; 
cS.tau_k = 0.25;
cS.tau_l = 0.08;
cS.tau_c = 0.04;
% cS.gov_exp_frac_Y = 0.15; % G将内生等于T
cS.A = 1.0;

%% 2. 求解RA模型的稳态
[ss, eq_found] = main_olg_v12_utils_simplified.solve_steady_state_RA(cS);
if ~eq_found
    error('无法求解RA稳态，审计中止');
end

fprintf('\n--- RA稳态结果 ---\n');
fprintf('Y=%.4f, K=%.4f, C=%.4f, I=%.4f, G=%.4f\n\n', ss.Y, ss.K, ss.C_value, ss.I_gross, ss.G);

%% 3. 构建T型账户进行审计
fprintf('--- RA模型 T型账户审计 ---\n\n');

% --- 企业部门 (Firms) ---
fprintf('=============== 1. 企业部门 (Firms) ===============\n');
Firm_Sources = ss.Y;
Firm_Uses = ss.w * ss.L + ss.r_rental * ss.K;
fprintf('   资金来源 (Sources):\n     + 销售产品 (Y):                 %9.4f\n', Firm_Sources);
fprintf('   --------------------------------------------------\n     总来源:                         %9.4f\n\n', Firm_Sources);
fprintf('   资金运用 (Uses):\n     + 支付工资 (w*L):               %9.4f\n     + 支付资本租金 (r_rental*K):    %9.4f\n', ss.w * ss.L, ss.r_rental * ss.K);
fprintf('   --------------------------------------------------\n     总运用:                         %9.4f\n\n', Firm_Uses);
Firm_Budget_Residual = Firm_Sources - Firm_Uses;
fprintf('   企业部门预算残差 (利润):          %9.4e\n', Firm_Budget_Residual);
if abs(Firm_Budget_Residual) < 1e-6, fprintf('   ✅ 企业部门预算平衡 (利润为零)。\n\n');
else, fprintf('   ⚠️ 企业部门预算不平衡！\n\n'); end

% --- 家庭部门 (Households) ---
fprintf('=============== 2. 家庭部门 (Households) ===============\n');
HH_Gross_Factor_Income = ss.w * ss.L + ss.r_rental * ss.K;
HH_Direct_Taxes = ss.T_labor + ss.T_capital;
Depreciation = ss.I_gross;
HH_Net_Disposable_Income = HH_Gross_Factor_Income - HH_Direct_Taxes - Depreciation;
HH_Consumption_Expenditure = ss.C_value * (1 + cS.tau_c);
HH_Net_Saving = 0; 
HH_Uses_of_NDI = HH_Consumption_Expenditure + HH_Net_Saving;

fprintf('   资金来源 (净可支配收入 Net Disposable Income):\n');
fprintf('     + 总要素收入 (w*L + r_rental*K): %9.4f\n', HH_Gross_Factor_Income);
fprintf('     - 直接税 (T_l + T_k):            %9.4f\n', HH_Direct_Taxes);
fprintf('     - 资本折旧 (δK):                 %9.4f\n', Depreciation);
fprintf('   --------------------------------------------------\n     = 净可支配总收入:               %9.4f\n\n', HH_Net_Disposable_Income);
fprintf('   资金运用 (Uses of Net Disposable Income):\n');
fprintf('     + 消费总支出 (C*(1+tau_c)):     %9.4f\n', HH_Consumption_Expenditure);
fprintf('     + 净储蓄 (S_net_pvt):             %9.4f\n', HH_Net_Saving);
fprintf('   --------------------------------------------------\n     = 总运用:                       %9.4f\n\n', HH_Uses_of_NDI);

HH_Budget_Residual = HH_Net_Disposable_Income - HH_Uses_of_NDI;
fprintf('   家庭部门预算残差 (NDI - C_exp - S_net): %9.4e\n', HH_Budget_Residual);
if abs(HH_Budget_Residual) < 1e-6, fprintf('   ✅ 家庭部门预算平衡。\n\n');
else, fprintf('   ⚠️ 家庭部门预算不平衡！\n\n'); end

% --- 政府部门 (Government) ---
fprintf('========== 3. 政府部门 (Consolidated Government) ==========\n');
Gov_Sources = ss.T_total;
Gov_Uses = ss.G; 
Gov_Primary_Surplus = Gov_Sources - Gov_Uses;
fprintf('   资金来源 (Sources):\n     + 总税收收入 (T_total):         %9.4f\n', Gov_Sources);
fprintf('   --------------------------------------------------\n     总来源:                         %9.4f\n\n', Gov_Sources);
fprintf('   资金运用 (Uses):\n     + 政府购买 (G):                 %9.4f\n', Gov_Uses);
fprintf('   --------------------------------------------------\n     总运用 (基础支出):              %9.4f\n\n', Gov_Uses);
fprintf('   政府基础盈余 (T - G):             %9.4f\n\n', Gov_Primary_Surplus);

% --- 最终审计：市场出清 ---
fprintf('============== 4. 最终审计: 商品市场出清 ============== \n');
Total_Uses_Goods_Market = ss.C_value + ss.I_gross + ss.G;

fprintf('   商品市场供给 (Y):\n     + 总产出:                       %9.4f\n\n', ss.Y);
fprintf('   商品市场需求 (C + I_gross + G):\n');
fprintf('     + 总消费 (C):                   %9.4f\n', ss.C_value);
fprintf('     + 总投资 (I_gross = δK):        %9.4f\n', ss.I_gross);
fprintf('     + 政府购买 (G):                 %9.4f\n', ss.G);
fprintf('   --------------------------------------------------\n     = 总需求:                       %9.4f\n\n', Total_Uses_Goods_Market);

Final_Residual = ss.Y - Total_Uses_Goods_Market;
fprintf('   最终残差 (Y - C - I_gross - G):     %9.4e\n', Final_Residual);
if abs(Final_Residual) < 1e-7, fprintf('   ✅ 恭喜！所有账户平衡，商品市场出清！\n');
else, fprintf('   ⚠️ 审计失败！核心会计逻辑存在问题！\n'); end

% --- END OF FILE debug_steady_state_RA.m ---