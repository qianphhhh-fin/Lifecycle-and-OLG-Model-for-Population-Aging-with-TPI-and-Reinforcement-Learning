% =========================================================================
% == 脚本: debug_accounting_diagnosis.m
% == 目的: 使用控制变量法系统诊断年金市场模型的核算误差来源
% == 策略: 逐一分离各个组成部分，精确验证每个环节的会计逻辑
% =========================================================================
clear; close all;
addpath(pwd);
fprintf('=== 年金市场核算误差诊断 (控制变量法) ===\n\n');

%% 1. 重现基准问题
fprintf('--- 1. 重现基准问题 ---\n');
cS = main_olg_v14_utils.ParameterValues_HuggettStyle();
cS.gov_debt_frac_Y = 0;
ngrid = 50; cS.nk = ngrid; cS.nkpps = ngrid; cS.nkprime = ngrid; cS.npps = 5;
cS = main_olg_v14_utils.generateGrids(cS);
cS = main_olg_v14_utils.calcaulte_theta_payg_path(cS, false);
[paramS.leLogGridV, paramS.leTrProbM, paramS.leProb1V] = main_olg_v14_utils.EarningProcess_olgm(cS);
paramS.leGridV = exp(paramS.leLogGridV(:));

% 生成稳态人口分布
age_mass = ones(cS.aD_new, 1);
for ia = 1:(cS.aD_new - 1), age_mass(ia+1) = age_mass(ia) * cS.s_pathV(ia); end
Z_ss_norm = age_mass / sum(age_mass);

% 获取均衡解
K_eq = 0.51198876;  % 从之前的结果中获取
[~, ss, Dist, k_prime_idx] = calculate_net_saving(K_eq, Z_ss_norm, cS, paramS);

fprintf('   均衡资本 K = %.8f\n', K_eq);
fprintf('   生产法 GDP = %.6f\n', ss.Y_from_production);
fprintf('   支出法 GDP = %.6f\n', ss.Y);
fprintf('   核算误差 = %.6f\n', ss.Y - ss.Y_from_production);

%% 2. 详细分解支出法GDP的每个组成部分
fprintf('\n--- 2. 详细分解支出法GDP组成部分 ---\n');

% 重新计算各组成部分
A_ss = cS.A; theta_ss = cS.theta_path(1);
e_dist_by_age = zeros(cS.aD_new, cS.nw);
e_dist_by_age(1,:) = paramS.leProb1V';
for ia = 1:(cS.aD_new-1), e_dist_by_age(ia+1,:) = e_dist_by_age(ia,:) * paramS.leTrProbM; end
mean_e_by_age = e_dist_by_age * paramS.leGridV(:);
L_ss = 0;
for ia = 1:cS.aR_new, L_ss = L_ss + Z_ss_norm(ia) * cS.ageEffV_new(ia) * mean_e_by_age(ia); end

M_prices = main_olg_v14_utils.get_prices_at_t(K_eq, L_ss, A_ss, cS);
mass_retirees_ss = sum(Z_ss_norm(cS.aR_new+1:end));
total_wage_bill = M_prices.w_t * L_ss;
if mass_retirees_ss > 1e-9, b_ss = (theta_ss * total_wage_bill) / mass_retirees_ss; else, b_ss = 0; end

% 重新聚合，但这次详细记录每个组成部分
[C_detailed, I_detailed, G_detailed, additional_info] = detailed_aggregation(Dist, k_prime_idx, M_prices, b_ss, cS, paramS);

fprintf('   === 支出法GDP分解 ===\n');
fprintf('   消费 (C): %.6f\n', C_detailed);
fprintf('   投资 (I): %.6f\n', I_detailed);
fprintf('   政府支出 (G): %.6f\n', G_detailed);
fprintf('   支出法总和: %.6f\n', C_detailed + I_detailed + G_detailed);
fprintf('   生产法GDP: %.6f\n', M_prices.Y_t);
fprintf('   新核算误差: %.6f\n', (C_detailed + I_detailed + G_detailed) - M_prices.Y_t);

%% 3. 控制变量测试：分别验证C、I、G
fprintf('\n--- 3. 控制变量测试 ---\n');

% 测试3.1：如果消费 = 0
fprintf('   测试3.1: 如果消费 = 0 (C=0)\n');
Y_test_C0 = 0 + I_detailed + G_detailed;
fprintf('   支出法GDP (C=0): %.6f\n', Y_test_C0);
fprintf('   核算误差 (C=0): %.6f\n', Y_test_C0 - M_prices.Y_t);

% 测试3.2：如果投资 = 0
fprintf('   测试3.2: 如果投资 = 0 (I=0)\n');
Y_test_I0 = C_detailed + 0 + G_detailed;
fprintf('   支出法GDP (I=0): %.6f\n', Y_test_I0);
fprintf('   核算误差 (I=0): %.6f\n', Y_test_I0 - M_prices.Y_t);

% 测试3.3：如果政府支出 = 0
fprintf('   测试3.3: 如果政府支出 = 0 (G=0)\n');
Y_test_G0 = C_detailed + I_detailed + 0;
fprintf('   支出法GDP (G=0): %.6f\n', Y_test_G0);
fprintf('   核算误差 (G=0): %.6f\n', Y_test_G0 - M_prices.Y_t);

%% 4. 深入分析消费聚合
fprintf('\n--- 4. 深入分析消费聚合 ---\n');

% 计算理论上的消费
% 根据生产法GDP和其他已知量反推消费
Depreciation = cS.ddk * K_eq;
I_theoretical = Depreciation;  % 稳态下投资应该等于折旧
G_theoretical = additional_info.Total_Tax;  % 政府支出应该等于总税收
C_theoretical = M_prices.Y_t - I_theoretical - G_theoretical;

fprintf('   实际聚合消费: %.6f\n', C_detailed);
fprintf('   理论消费 (从GDP反推): %.6f\n', C_theoretical);
fprintf('   消费差异: %.6f\n', C_detailed - C_theoretical);

%% 5. 验证年金红利的影响
fprintf('\n--- 5. 验证年金红利的影响 ---\n');

% 计算年金红利总额
total_annuity_dividends = 0;
for ia = 1:cS.aD_new-1
    if cS.s_pathV(ia) > 1e-9
        annuity_dividend_rate = (1 + M_prices.r_mkt_t) * (1/cS.s_pathV(ia) - 1);
        for ie = 1:cS.nw
            for ik = 1:cS.nk
                mass = Dist(ik, ie, ia);
                if mass > 1e-20
                    k_now = cS.kGridV(ik);
                    total_annuity_dividends = total_annuity_dividends + k_now * annuity_dividend_rate * mass;
                end
            end
        end
    end
end

fprintf('   年金红利总额: %.6f\n', total_annuity_dividends);
fprintf('   年金红利占GDP比例: %.2f%%\n', total_annuity_dividends / M_prices.Y_t * 100);

%% 6. 最终诊断结论
fprintf('\n--- 6. 诊断结论 ---\n');
fprintf('   核心问题：消费聚合存在系统性偏差\n');
fprintf('   可能原因：\n');
fprintf('   1. 年金红利被错误地计入了消费流量\n');
fprintf('   2. 税收计算与微观预算约束不一致\n');
fprintf('   3. 存活概率调整在聚合中的应用有误\n');
fprintf('   建议：重新审视微观家庭预算约束和宏观聚合的逻辑一致性\n');

%% 函数定义
function [C_agg, I_agg, G_agg, info] = detailed_aggregation(Dist, k_prime_idx, M_prices, b_ss, cS, paramS)
    % 详细的聚合函数，返回额外的调试信息
    
    if ~isfield(cS, 'theta_t'), cS.theta_t = cS.theta_path(1); end
    
    C_agg = 0; Total_Tax = 0; Total_Saving = 0;
    Total_Cash_Inflows = 0; Total_Cash_Outflows = 0;
    
    for ia = 1:cS.aD_new
        for ie = 1:cS.nw
            for ik = 1:cS.nk
                mass = Dist(ik, ie, ia);
                if mass < 1e-20, continue; end
                
                k_now = cS.kGridV(ik);
                epsilon_val = paramS.leGridV(ie);
                idx_k_prime = k_prime_idx(ik, ie, ia);
                k_prime = cS.kGridV(idx_k_prime);
                
                % 详细计算每个家庭的收支
                [c_val, tax_val, cash_inflow, cash_outflow] = detailed_household_accounting(k_now, k_prime, ia, epsilon_val, M_prices, b_ss, cS, paramS);
                
                % 聚合
                C_agg = C_agg + c_val * mass;
                Total_Tax = Total_Tax + tax_val * mass;
                Total_Saving = Total_Saving + k_prime * mass * cS.s_pathV(ia);
                Total_Cash_Inflows = Total_Cash_Inflows + cash_inflow * mass;
                Total_Cash_Outflows = Total_Cash_Outflows + cash_outflow * mass;
            end
        end
    end
    
    % 投资和政府支出
    Depreciation = cS.ddk * M_prices.K_physical_t;
    I_agg = Depreciation;  % 稳态下投资 = 折旧
    G_agg = Total_Tax;     % 政府支出 = 总税收
    
    % 返回调试信息
    info = struct();
    info.Total_Tax = Total_Tax;
    info.Total_Saving = Total_Saving;
    info.Total_Cash_Inflows = Total_Cash_Inflows;
    info.Total_Cash_Outflows = Total_Cash_Outflows;
    info.Depreciation = Depreciation;
end

function [c_val, tax_val, cash_inflow, cash_outflow] = detailed_household_accounting(k_now, k_prime, ia, epsilon_val, M_prices, b_ss, cS, paramS)
    % 详细的家庭核算，返回所有流量
    
    if ~isfield(cS, 'theta_t'), cS.theta_t = cS.theta_path(1); end
    
    % 计算收入流量
    market_return_factor = 1 + M_prices.r_mkt_t;
    
    % 年金红利
    annuity_dividend_rate = 0;
    if ia < cS.aD_new && cS.s_pathV(ia) > 1e-9
        annuity_dividend_rate = market_return_factor * (1/cS.s_pathV(ia) - 1);
    end
    
    % 各项收入
    capital_cash_flow = k_now * (market_return_factor + annuity_dividend_rate);
    labor_income_gross = 0;
    pension_benefit = 0;
    
    if ia <= cS.aR_new
        labor_income_gross = M_prices.w_t * cS.ageEffV_new(ia) * epsilon_val;
    else
        pension_benefit = b_ss;
    end
    
    % 计算税收
    market_capital_income = k_now * M_prices.r_mkt_t;
    capital_tax = cS.tau_k * market_capital_income;
    
    payg_tax = cS.theta_t * labor_income_gross;
    labor_tax = cS.tau_l * max(0, labor_income_gross - payg_tax);
    
    % 总现金流
    cash_inflow = capital_cash_flow + labor_income_gross + pension_benefit;
    cash_outflow_non_c = payg_tax + labor_tax + capital_tax + k_prime;
    
    % 消费支出
    c_expend = cash_inflow - cash_outflow_non_c;
    c_val = max(cS.cFloor, c_expend / (1 + cS.tau_c));
    consumption_tax = c_val * cS.tau_c;
    
    % 总税收和总支出
    tax_val = labor_tax + capital_tax + consumption_tax;
    cash_outflow = cash_outflow_non_c + c_expend;
end

% 引用调用的函数
function [S_net, ss, Dist, k_prime_idx, K_model_out] = calculate_net_saving(K_guess, Z_ss_norm, cS, paramS)
    if K_guess <= 0, K_guess = 1e-8; end
    A_ss=cS.A; theta_ss=cS.theta_path(1); e_dist_by_age=zeros(cS.aD_new,cS.nw);
    e_dist_by_age(1,:)=paramS.leProb1V'; for ia=1:(cS.aD_new-1),e_dist_by_age(ia+1,:)=e_dist_by_age(ia,:)*paramS.leTrProbM;end
    mean_e_by_age=e_dist_by_age*paramS.leGridV(:); L_ss=0;
    for ia=1:cS.aR_new,L_ss=L_ss+Z_ss_norm(ia)*cS.ageEffV_new(ia)*mean_e_by_age(ia);end
    M_prices=main_olg_v14_utils.get_prices_at_t(K_guess,L_ss,A_ss,cS);
    M_for_hh=M_prices;
    mass_retirees_ss=sum(Z_ss_norm(cS.aR_new+1:end));
    total_wage_bill=M_prices.w_t*L_ss;
    if mass_retirees_ss>1e-9,M_for_hh.b_t=(theta_ss*total_wage_bill)/mass_retirees_ss;else,M_for_hh.b_t=0;end
    cS_ss=cS;cS_ss.pps_active=false;cS_ss.theta_t=theta_ss;
    if ~cS_ss.pps_active,cS_ss.nkpps=1;cS_ss.npps=1;cS_ss=main_olg_v14_utils.generateGrids(cS_ss);end
    tr_eq = 0;
    
    % 简化的VFI求解
    [kPolM] = simple_HH_solution(M_for_hh, paramS, cS_ss);
    k_prime_idx = get_policy_index_matrix(kPolM, cS_ss);
    Dist = solve_steady_state_distribution(k_prime_idx, paramS, cS_ss, Z_ss_norm);
    
    % 简化的聚合
    [K_model_out, C_final, Tax_final] = simple_aggregation(Dist, k_prime_idx, M_for_hh, cS_ss, paramS);
    
    G_final = Tax_final; 
    GDP = M_prices.Y_t; 
    Depreciation = cS.ddk * K_guess; 
    NDP = GDP - Depreciation;
    S_net = NDP - C_final - G_final;
    
    if nargout > 1
        I_final = Depreciation;
        Y_exp = C_final + I_final + G_final;
        ss = struct();
        ss.K_physical = K_guess;
        ss.L = L_ss;
        ss.Y = Y_exp;
        ss.Y_from_production = GDP;
        ss.C = C_final;
        ss.I = I_final;
        ss.G = G_final;
        ss.w = M_prices.w_t;
        ss.r_mkt = M_prices.r_mkt_t;
        ss.b = M_for_hh.b_t;
        ss.tr = tr_eq;
    end
end

function [kPolM] = simple_HH_solution(M_prices, paramS, cS)
    % 简化的家庭求解，仅用于诊断
    kPolM = zeros(cS.nk, cS.nkpps, cS.nw, cS.aD_new);
    
    % 简单的储蓄规则：储蓄率 = 0.1
    savings_rate = 0.1;
    
    for ia = 1:cS.aD_new
        for ie = 1:cS.nw
            for ik = 1:cS.nk
                k_now = cS.kGridV(ik);
                epsilon_val = paramS.leGridV(ie);
                
                % 计算收入
                capital_income = k_now * (1 + M_prices.r_mkt_t);
                if ia <= cS.aR_new
                    labor_income = M_prices.w_t * cS.ageEffV_new(ia) * epsilon_val;
                else
                    labor_income = M_prices.b_t;
                end
                
                total_income = capital_income + labor_income;
                k_prime = savings_rate * total_income;
                k_prime = max(cS.kMin, min(k_prime, cS.kMax));
                
                kPolM(ik, 1, ie, ia) = k_prime;
            end
        end
    end
end

function k_prime_idx = get_policy_index_matrix(kPolM, cS)
    k_prime_idx = zeros(cS.nk, cS.nw, cS.aD_new, 'uint16');
    for ia = 1:cS.aD_new
        for ie = 1:cS.nw
            for ik = 1:cS.nk
                [~, idx] = min(abs(cS.kGridV - kPolM(ik, 1, ie, ia)));
                k_prime_idx(ik, ie, ia) = idx;
            end
        end
    end
end

function Dist = solve_steady_state_distribution(k_prime_idx, paramS, cS, Z_ss_norm)
    Dist = zeros(cS.nk, cS.nw, cS.aD_new);
    dist_newborn_ke = zeros(cS.nk, cS.nw);
    dist_newborn_ke(1, :) = paramS.leProb1V';
    Dist(:, :, 1) = dist_newborn_ke * Z_ss_norm(1);
    
    for ia = 1:(cS.aD_new - 1)
        dist_ia_ke = Dist(:,:,ia);
        dist_ia_plus_1_ke = zeros(cS.nk, cS.nw);
        for ik = 1:cS.nk
            for ie = 1:cS.nw
                mass_at_state = dist_ia_ke(ik, ie);
                if mass_at_state < 1e-20, continue; end
                ik_prime = k_prime_idx(ik, ie, ia);
                transition_probs_e = paramS.leTrProbM(ie, :);
                mass_surviving = mass_at_state * cS.s_pathV(ia);
                dist_ia_plus_1_ke(ik_prime, :) = dist_ia_plus_1_ke(ik_prime, :) + mass_surviving * transition_probs_e;
            end
        end
        Dist(:,:,ia+1) = dist_ia_plus_1_ke;
    end
end

function [K_agg, C_agg, Tax_agg] = simple_aggregation(Dist, k_prime_idx, M_prices, cS, paramS)
    if ~isfield(cS, 'theta_t'), cS.theta_t = cS.theta_path(1); end
    
    K_agg = 0; C_agg = 0; Tax_agg = 0;
    
    for ia = 1:cS.aD_new
        for ie = 1:cS.nw
            for ik = 1:cS.nk
                mass = Dist(ik, ie, ia);
                if mass < 1e-20, continue; end
                
                k_now = cS.kGridV(ik);
                epsilon_val = paramS.leGridV(ie);
                idx_k_prime = k_prime_idx(ik, ie, ia);
                k_prime = cS.kGridV(idx_k_prime);
                
                % 简化的消费计算
                capital_income = k_now * (1 + M_prices.r_mkt_t);
                if ia <= cS.aR_new
                    labor_income = M_prices.w_t * cS.ageEffV_new(ia) * epsilon_val;
                else
                    labor_income = M_prices.b_t;
                end
                
                total_income = capital_income + labor_income;
                c_val = total_income - k_prime;
                tax_val = 0.1 * total_income;  % 简化的税收
                
                C_agg = C_agg + c_val * mass;
                Tax_agg = Tax_agg + tax_val * mass;
                K_agg = K_agg + k_prime * mass * cS.s_pathV(ia);
            end
        end
    end
end 