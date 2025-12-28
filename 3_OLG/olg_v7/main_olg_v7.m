% --- START OF FILE main_olg_v7.m ---

% OLG 模型 V7 (PPS 所得税递延, VFI w k_pps 状态, PPS缴费规则化):
% 目标: 寻找最大可持续的 PAYG 替代率 (rho_prime_payg)
% 定义: rho_prime_payg = b_payg / avg_worker_gross_wage (PAYG福利 / 平均工人税前总工资)
% 约束: PPS 缴费可从所得税前扣除 (所得税率为 tau_l)。
%       tau_l 内生调整以平衡政府一般预算 (TR_gov = 0)。
%       PAYG 税率 (theta_payg) 内生决定，但有上限 cS.theta_payg_max。
%       VFI 状态变量仍然包含 k_pps (PPS资产)。
%       PPS 缴费行为规则化: 按年龄组特定比例缴存，但受法定上限约束。

clc; % 清除命令行窗口
clear; % 清除工作区变量
close all; % 关闭所有图形窗口
fprintf('=== OLG 模型 V7 (PPS 税收递延, VFI w k_pps, PPS缴费规则化): 最大可持续 Rho_prime_payg ===\n');
fprintf('    (Rho_prime_payg = b_payg / avg_worker_gross_wage)\n');
fprintf('    (TR_gov=0, tau_l 内生, theta_payg 有上限, VFI 状态: k, k_pps, eps)\n');
fprintf('    (PPS缴费: 基于年龄组的固定比例规则，有法定上限)\n');

%% 1. 固定模型参数与人口设置
fprintf('\n--- 1. 初始化参数 ---\n');
cS = main_olg_v7_utils.ParameterValues_HuggettStyle(); % 调用v7的参数设置函数
paramS = struct(); % 初始化一个结构体用于存储派生参数

fprintf('参数已加载。nk=%d, nkpps=%d, nw=%d。\n', cS.nk, cS.nkpps, cS.nw); % 打印网格点数
fprintf('年度年龄范围: %d-%d。模型年龄组数: %d。\n', cS.age1_orig, cS.ageLast_orig, cS.aD_new); % 打印年龄信息
fprintf('固定税率: tau_k=%.2f, tau_c=%.2f。G/Y=%.2f, B/Y=%.2f。\n', cS.tau_k, cS.tau_c, cS.gov_exp_frac_Y, cS.gov_debt_frac_Y); % 打印固定税率和财政目标
fprintf('PAYG 税率上限 (theta_payg_max): %.3f\n', cS.theta_payg_max); % 打印PAYG税率上限
fprintf('所得税率 tau_l 范围: [%.3f, %.3f], 总劳动税上限: %.3f\n', cS.tau_l_min, cS.tau_l_max, cS.max_total_labor_tax); % 打印所得税率范围和总税负上限
fprintf('PPS 年度缴费上限 (绝对值): %.2f, 比例上限 (法定): %.2f\n', cS.pps_annual_contrib_limit, cS.pps_max_contrib_frac); % 打印PPS缴费上限
fprintf('V7 PPS 缴费规则 (部分年龄组的比例): \n'); % 打印PPS规则化缴费的部分信息
if cS.aD_new > 0
    for i_disp_pps = 1:min(5, cS.aD_new) % 最多显示前5个年龄组的规则比例
        fprintf('  年龄组 %d: %.3f\n', i_disp_pps, cS.pps_fixed_contrib_schedule_frac(i_disp_pps));
    end
    if cS.aD_new > 5, fprintf('  ...\n'); end % 如果年龄组过多，用省略号表示
end


%% 2. 模拟人口动态至稳态
fprintf('\n--- 2. 模拟人口动态 ---\n');
popS = main_olg_v7_utils.initPopulation(cS); % 初始化人口结构
popS = main_olg_v7_utils.populationDynamics(popS, cS); % 模拟人口动态演进
[Z_ss, ~, ~, ~] = main_olg_v7_utils.detectSteadyStatePopulation(popS, cS); % 检测稳态人口分布
paramS.Z_ss_counts = Z_ss; % 存储稳态人口计数 (按年龄组)
Z_ss_total = sum(Z_ss); % 计算稳态总人口
Z_ss_norm_group = zeros(cS.aD_new,1); % 初始化归一化的年龄组人口分布
if Z_ss_total > 1e-9 % 确保总人口不为零
    Z_ss_norm_group = Z_ss / Z_ss_total; % 计算归一化分布
end
paramS.ageMassV = Z_ss_norm_group(:); % 存储为列向量 (年龄组人口占比)
paramS.mass_workers_group = sum(paramS.ageMassV(1:cS.aR_new)); % 计算工作年龄人口占比 (基于年龄组)
if paramS.mass_workers_group < 1e-9 % 检查工作人口占比是否过小
    error('模型中稳态工作人口占比为零或过小。'); 
end

% 将年龄组人口分布转换为年度人口分布 (用于计算遗赠等)
Z_ss_norm_annual = zeros(cS.aD_orig,1); % 初始化年度人口分布向量
if Z_ss_total > 1e-9 % 如果稳态总人口有效
    for a_new_map_idx = 1:cS.aD_new % 遍历每个年龄组
        annual_indices_in_group = cS.physAgeMap{a_new_map_idx}; % 获取该年龄组包含的年度年龄索引
        group_mass_fraction = Z_ss_norm_group(a_new_map_idx); % 该年龄组的人口占比
        num_years_in_this_group = length(annual_indices_in_group); % 该年龄组包含的年数
        if num_years_in_this_group > 0 % 如果年数有效
            mass_per_year_in_group = group_mass_fraction / num_years_in_this_group; % 该年龄组内每年的平均人口占比
            Z_ss_norm_annual(annual_indices_in_group) = mass_per_year_in_group; % 分配给对应的年度年龄
        end
    end
    if sum(Z_ss_norm_annual) > 1e-9 && abs(sum(Z_ss_norm_annual) - 1.0) > 1e-6 % 如果加总不为1 (有精度误差)
        Z_ss_norm_annual = Z_ss_norm_annual / sum(Z_ss_norm_annual); % 重新归一化
    elseif sum(Z_ss_norm_annual) < 1e-9 % 如果加总过小 (异常情况)
        Z_ss_norm_annual(:) = 1/cS.aD_orig; % 使用均匀分布作为备用
    end
else % 如果稳态总人口为零 (异常情况)
    Z_ss_norm_annual(:) = 1/cS.aD_orig; % 使用均匀分布作为备用
end
paramS.Z_ss_norm_annual = Z_ss_norm_annual; % 存储年度人口分布

% 计算稳态人口增长率 (用于政府债务动态)
if Z_ss_total > 1e-9 && length(popS.totalPop) > 1 && popS.totalPop(end-1) > 1e-9 % 如果历史总人口数据有效
    % 从人口动态模拟的最后两期计算每期(年龄组时长)的人口增长因子
    pop_growth_factor_per_group_period = popS.totalPop(end) / popS.totalPop(end-1);
    % 转换为年化人口增长率: (增长因子)^(1/每组年数) - 1
    paramS.popGrowthForDebt = pop_growth_factor_per_group_period^(1/cS.yearStep) - 1; 
else % 如果历史数据不足，使用参数中的原始设定值
    paramS.popGrowthForDebt = cS.popGrowth_orig; 
end
fprintf('人口参数计算完毕。年化稳态人口增长率 (用于债务): %.4f\n', paramS.popGrowthForDebt);
fprintf('稳态工人占比 (基于年龄组人口): %.4f\n', paramS.mass_workers_group);

%% 3. 预计算劳动供给和禀赋过程
fprintf('\n--- 3. 预计算劳动 ---\n');
% 生成劳动效率冲击的网格、转移矩阵和初始分布
[paramS.leLogGridV, paramS.leTrProbM, paramS.leProb1V] = main_olg_v7_utils.EarningProcess_olgm(cS);
paramS.leGridV = exp(paramS.leLogGridV(:)); % 效率水平 (取指数)
paramS.ageEffV_new = cS.ageEffV_new; % 年龄效率剖面 (按年龄组)

% 模拟个体的劳动效率冲击路径 (年度)
eIdxM = main_olg_v7_utils.LaborEndowSimulation_olgm(cS, paramS);

% 计算总体人均有效劳动供给 L_per_capita
[~, L_per_capita] = main_olg_v7_utils.LaborSupply_Huggett(eIdxM, cS, paramS, paramS.ageMassV);
fprintf('总劳动供给 (L, 效率单位, 总体人均): %.4f\n', L_per_capita);
if L_per_capita <= 0 % 如果劳动供给为零或负 (异常)
    L_per_capita = 1e-6; % 设置为一个非常小的正值
    warning('L_per_capita 为零或负，已重置为1e-6。');
end
paramS.L_per_capita = L_per_capita; % 存储计算得到的总劳动供给

%% 4. 寻找最大可持续替代率 rho_prime_payg_max (VFI w k_pps state, PPS缴费规则化)
fprintf('\n--- 4. 寻找最大可持续替代率 rho_prime_payg_max (VFI w k_pps state, PPS规则化) ---\n');

% --- 搜索参数设置 ---
rho_payg_min = 0.01; % PAYG替代率搜索下限
rho_payg_max_search_upper_bound = 0.7; % PAYG替代率搜索上限的初始猜测 (可调整)
max_iter_rho_search = 25; % 替代率二分法搜索的最大迭代次数
tol_rho_search = 1e-3;    % 替代率搜索的收敛容忍度

K_global_guess = 15.0; % 总资本存量的初始猜测值 (用于内层均衡求解)

% --- 初始化搜索循环变量 ---
rho_low = rho_payg_min; % 当前搜索区间的下限
rho_high = rho_payg_max_search_upper_bound; % 当前搜索区间的上限
rho_prime_payg_optimal = rho_payg_min; % 当前找到的最优可持续替代率 (初始化为最小值)
K_optimal = K_global_guess; % 对应最优替代率的资本存量
tau_l_optimal = cS.tau_l_init_guess; % 对应最优替代率的所得税率
theta_payg_optimal_calc = 0; % 对应最优替代率的理论PAYG税率 (未考虑上限前)
final_eq_solution_details = struct(); % 存储最优均衡的详细信息
final_eq_solution_found = false; % 标记是否找到了至少一个可行解

fprintf('开始搜索最大可持续替代率 rho_prime_payg ...\n');
fprintf('IterRho | Rho_low  | Rho_high | Rho_try  | K_tot_sol| K_pps_sol| Tau_l_sol| Theta_g_req| GBC_Res  | Feasible | Time(s)\n');
fprintf('---------------------------------------------------------------------------------------------------------------------------------------\n');
paramS_for_inner_loop = paramS; % 复制一份参数结构体，用于内层循环 (可能在循环中修改)

for iter_rho = 1:max_iter_rho_search % 开始二分法搜索循环
    rho_prime_try = (rho_low + rho_high) / 2; % 计算当前尝试的替代率 (区间中点)
    iter_rho_start_time = tic; % 记录当次迭代开始时间

    % --- 控制内层循环打印信息的标志 ---
    if iter_rho > 1 % 如果不是第一次替代率迭代
        paramS_for_inner_loop.suppress_inner_print_header = true; % 不打印内层K-tau_l循环的表头
        paramS_for_inner_loop.suppress_initial_theta_print = true; % 不打印内层循环初始theta_req的检查信息
    else % 如果是第一次替代率迭代
        paramS_for_inner_loop.suppress_inner_print_header = false; % 打印内层表头
        paramS_for_inner_loop.suppress_initial_theta_print = false; % 打印初始theta_req检查
    end

    % --- 调用内层均衡求解器 ---
    % 对于给定的 rho_prime_try，求解均衡的 K, tau_l，并检查可行性
    [K_solution_iter, tau_l_solution_iter, gbc_residual_iter, eq_found_for_rho_try_iter, solution_details_iter] = ...
        main_olg_v7_utils.solve_K_tau_l_for_rho_prime(rho_prime_try, K_global_guess, cS, paramS_for_inner_loop, eIdxM);

    iter_rho_time = toc(iter_rho_start_time); % 计算当次迭代耗时
    
    % 从解的细节中获取PAYG理论税率和PPS资本 (如果存在)
    current_theta_payg_calc_iter = NaN; K_pps_model_from_inner_iter = NaN;
    if isfield(solution_details_iter, 'theta_payg_required_before_cap')
        current_theta_payg_calc_iter = solution_details_iter.theta_payg_required_before_cap;
    end
    if isfield(solution_details_iter, 'K_model_pps')
        K_pps_model_from_inner_iter = solution_details_iter.K_model_pps;
    end

    % 打印当次替代率搜索的迭代结果
    fprintf('%7d | %8.4f | %8.4f | %8.4f | %8.4f | %8.4f | %8.4f | %11.4f | %8.2e | %8d | %.2f\n', ...
        iter_rho, rho_low, rho_high, rho_prime_try, ...
        K_solution_iter, K_pps_model_from_inner_iter, tau_l_solution_iter, ... 
        current_theta_payg_calc_iter, gbc_residual_iter, eq_found_for_rho_try_iter, iter_rho_time);

    % --- 更新搜索区间 ---
    if eq_found_for_rho_try_iter % 如果当前尝试的替代率是可持续的 (找到了均衡且满足约束)
        rho_prime_payg_optimal = rho_prime_try; % 更新最优替代率
        K_optimal = K_solution_iter; % 更新对应的K
        tau_l_optimal = tau_l_solution_iter; % 更新对应的tau_l
        theta_payg_optimal_calc = current_theta_payg_calc_iter; % 更新对应的理论theta
        final_eq_solution_details = solution_details_iter; % 存储均衡细节
        
        rho_low = rho_prime_try; % 提高搜索下限，尝试更高的替代率
        final_eq_solution_found = true; % 标记已找到至少一个可行解
        K_global_guess = K_solution_iter; % 使用本次求解的K作为下一次内层循环的K初值 (有助于加速收敛)
    else % 如果当前尝试的替代率不可持续
        rho_high = rho_prime_try; % 降低搜索上限
    end

    % 检查替代率搜索是否收敛
    if (rho_high - rho_low) < tol_rho_search
        fprintf('替代率 (rho_prime_payg) 搜索已收敛。\n');
        break; % 跳出替代率搜索循环
    end
end % 结束替代率搜索循环 (iter_rho)

% --- 处理搜索结果 ---
if ~final_eq_solution_found % 如果整个搜索过程没有找到任何一个可持续的替代率
    warning('未能找到满足所有条件的最大可持续替代率 rho_prime_payg。将使用最后记录的可行值或最小值。');
    % 检查是否曾经有过一个尝试成功的rho_prime_payg_optimal (即使最后一次二分迭代可能失败)
    % rho_prime_payg_optimal 初始化为 rho_payg_min
    if isempty(fields(final_eq_solution_details)) && rho_prime_payg_optimal > rho_payg_min + tol_rho_search / 2 
        % 如果 final_eq_solution_details 是空的，但 rho_prime_payg_optimal 曾被更新过 (说明至少有一个可行解被短暂记录过)
        % 这通常意味着 rho_low 曾被成功更新，但后续更高的 rho_try 失败了。
        % 此时 rho_prime_payg_optimal 记录的是最后一次成功的 rho_low。
        fprintf('最后记录的可行 rho_prime_payg_optimal = %.4f，尝试重新计算其均衡状态...\n', rho_prime_payg_optimal);
        paramS_for_inner_loop.suppress_inner_print_header = false; % 允许打印
        paramS_for_inner_loop.suppress_initial_theta_print = false; % 允许打印
        [K_eq_fallback, tau_l_eq_fallback, ~, eq_found_fallback, final_eq_solution_details_fallback] = ... 
             main_olg_v7_utils.solve_K_tau_l_for_rho_prime(rho_prime_payg_optimal, K_global_guess, cS, paramS_for_inner_loop, eIdxM);
        
        if eq_found_fallback && ~isempty(fields(final_eq_solution_details_fallback)) && ...
           isfield(final_eq_solution_details_fallback,'MPL_gross') && ~isnan(final_eq_solution_details_fallback.MPL_gross)
            % 如果重新计算成功且结果有效
            K_optimal = K_eq_fallback; 
            tau_l_optimal = tau_l_eq_fallback; 
            final_eq_solution_details = final_eq_solution_details_fallback; % 使用这个重新计算的细节
            if isfield(final_eq_solution_details, 'theta_payg_required_before_cap')
                theta_payg_optimal_calc = final_eq_solution_details.theta_payg_required_before_cap;
            else
                theta_payg_optimal_calc = NaN;
            end
        else
            % 如果重新计算失败，这是一个更严重的问题
            error('重新计算最后的 rho_prime_payg_optimal 失败或未返回有效结果。检查模型或参数。');
        end
    elseif isempty(fields(final_eq_solution_details)) % 如果 final_eq_solution_details 始终为空 (即 rho_low 从未被更新)
        % 这意味着即使是最低的 rho_payg_min 也无法持续
        error('在指定的PAYG税率上限和财政约束下，即使最低的替代率 (rho_payg_min=%.3f) 也无法持续。请检查参数或约束。', rho_payg_min);
    end
    % 最终确定的均衡替代率、资本和所得税率
    rho_prime_payg_eq = rho_prime_payg_optimal; 
    K_eq = K_optimal; 
    tau_l_eq = tau_l_optimal; 
    fprintf('由于未找到严格收敛的最大可持续替代率，将使用 rho_prime_payg_eq = %.4f (理论theta_payg_req=%.4f)\n', rho_prime_payg_eq, theta_payg_optimal_calc);
else % 如果搜索成功收敛到了一个解
    rho_prime_payg_eq = rho_prime_payg_optimal; 
    K_eq = K_optimal; 
    tau_l_eq = tau_l_optimal;
    fprintf('找到的最大可持续替代率 rho_prime_payg_eq = %.4f (理论theta_payg_req=%.4f)\n', rho_prime_payg_eq, theta_payg_optimal_calc);
end


%% 5. 分析和绘制最终均衡结果
fprintf('\n--- 5. 最终均衡结果与绘图 (rho_prime_payg_eq=%.4f, tau_l_eq=%.4f, TR_gov=0) ---\n', rho_prime_payg_eq, tau_l_eq);

% --- 检查最终均衡解的有效性 ---
if isempty(fields(final_eq_solution_details)) || isnan(K_eq) || ...
   ~isfield(final_eq_solution_details, 'MPL_gross') || isnan(final_eq_solution_details.MPL_gross)
    error('最终均衡的详细信息未能获取或无效 (例如MPL_gross缺失)，无法进行分析。');
end

% --- 设置用于最终VFI和模拟的参数 ---
paramS_eq = paramS; % 复制一份基础参数
paramS_eq.tau_l = tau_l_eq; % 设置均衡所得税率

% 从 final_eq_solution_details 中获取实际PAYG税率，如果不存在则重新计算
if isfield(final_eq_solution_details, 'theta_payg') && isfinite(final_eq_solution_details.theta_payg)
    paramS_eq.theta_payg_actual_for_hh = final_eq_solution_details.theta_payg; % 家庭面临的实际PAYG税率
else 
    % 如果缺失，则基于均衡结果和约束重新计算一个
    warning('final_eq_solution_details.theta_payg (均衡中的实际PAYG税率) 未找到或无效，将基于均衡rho和约束重新计算。');
    % 使用 K_eq (即K_optimal) 和 L_per_capita 计算市场价格
    [~, temp_MPL_gross_for_theta_calc] = main_olg_v7_utils.HHPrices_Huggett(K_eq, paramS.L_per_capita, cS);
    temp_avg_worker_wage = (temp_MPL_gross_for_theta_calc * paramS.L_per_capita) / paramS.mass_workers_group;
    % temp_b_payg = rho_prime_payg_eq * temp_avg_worker_wage; % 均衡PAYG福利 (这一步其实不需要计算theta)
    
    % 理论所需theta (与外层搜索的 final_eq_solution_details.theta_payg_required_before_cap 应一致)
    temp_mass_retirees = sum(paramS.ageMassV(cS.aR_new+1:cS.aD_new));
    temp_theta_req = rho_prime_payg_eq * (temp_mass_retirees / paramS.mass_workers_group);
    temp_theta_req = max(0, temp_theta_req);
    
    temp_theta_act = min(temp_theta_req, cS.theta_payg_max); % 应用PAYG税率上限
    if (temp_theta_act + tau_l_eq) > cS.max_total_labor_tax % 应用总劳动税负上限
        temp_theta_act = max(0, cS.max_total_labor_tax - tau_l_eq);
    end
    paramS_eq.theta_payg_actual_for_hh = max(0, temp_theta_act); % 确保非负
    fprintf('  重新计算的实际PAYG税率 (用于最终VFI): %.4f (理论需求: %.4f)\n', paramS_eq.theta_payg_actual_for_hh, temp_theta_req);
end
paramS_eq.pps_tax_deferral_active = cS.pps_active; % 确保PPS税收递延状态正确传递

% --- 计算最终均衡的市场价格和家庭回报率 ---
[R_mkt_gross_factor_eq_final, MPL_gross_eq_final] = main_olg_v7_utils.HHPrices_Huggett(K_eq, paramS.L_per_capita, cS); % 使用均衡K_eq
r_mkt_gross_eq_final = R_mkt_gross_factor_eq_final - 1; % 均衡市场净回报率 (MPK-delta)
r_k_net_hh_eq_final = r_mkt_gross_eq_final * (1 - cS.tau_k); % 家庭面临的税后资本净回报率
R_k_net_factor_hh_eq_final = 1 + r_k_net_hh_eq_final; % 家庭面临的税后资本净回报因子

% --- 计算最终均衡的PAYG福利和转移支付 ---
avg_worker_gross_wage_eq_final = (MPL_gross_eq_final * paramS.L_per_capita) / paramS.mass_workers_group; % 均衡状态下的平均工人工资
b_payg_eq_final = rho_prime_payg_eq * avg_worker_gross_wage_eq_final; % 均衡PAYG福利 (每位退休者)
bV_eq_new_final = zeros(1, cS.aD_new); % 初始化福利向量 (按年龄组)
if cS.aR_new < cS.aD_new % 如果存在退休年龄组
    bV_eq_new_final(cS.aR_new + 1 : cS.aD_new) = b_payg_eq_final; % 为退休组设置福利
end

T_bequest_eq_final = 0; % 初始化意外遗赠 (人均)
if isfield(final_eq_solution_details, 'T_bequest_Model') && isfinite(final_eq_solution_details.T_bequest_Model)
    T_bequest_eq_final = final_eq_solution_details.T_bequest_Model; % 从均衡细节中获取
else
    warning('T_bequest_Model 未在 final_eq_solution_details 中找到或无效 (用于最终VFI)。将使用0。');
    % 如果需要，可以进行一次小的迭代来估算 T_bequest，但这会增加复杂性。
    % 对于一致性，最好确保 solve_K_tau_l_for_rho_prime 总是返回一个T_bequest。
end
TR_total_eq_final_for_vfi = T_bequest_eq_final; % 总转移支付 (因 TR_gov=0)

fprintf('最终 VFI 调用参数: MPL_gross=%.4f, tau_l=%.4f, theta_payg_actual=%.4f, TR_total=%.4f (T_bequest)\n', ...
    MPL_gross_eq_final, paramS_eq.tau_l, paramS_eq.theta_payg_actual_for_hh, TR_total_eq_final_for_vfi);

% --- 调用最终的VFI求解家庭问题 ---
fprintf('调用最终的 HHSolution_VFI_Huggett (V7)...\n');
[cPolM_eq, kPolM_eq, cPpsPolM_eq_rule, ~] = main_olg_v7_utils.HHSolution_VFI_Huggett(...
    R_k_net_factor_hh_eq_final, MPL_gross_eq_final, TR_total_eq_final_for_vfi, bV_eq_new_final, paramS_eq, cS);
% 返回的策略 cPolM_eq (消费量), kPolM_eq (非PPS储蓄), cPpsPolM_eq_rule (规则化PPS缴费)
% 均为 4D 矩阵 (nk, nkpps, nw, naD_new)

% --- 模拟最终均衡的资产和消费分布 ---
fprintf('模拟最终均衡的分布 (HHSimulation_olgm V7)...\n');
[kHistM_eq, kPpsHistM_eq, cHistM_eq] = main_olg_v7_utils.HHSimulation_olgm(...
    kPolM_eq, cPpsPolM_eq_rule, cPolM_eq, eIdxM, ... % 策略函数和效率路径
    R_k_net_factor_hh_eq_final, MPL_gross_eq_final, ... % 价格
    TR_total_eq_final_for_vfi, bV_eq_new_final, paramS_eq, cS); % 转移支付和参数

% --- 计算最终均衡的宏观总量 ---
% (与v6逻辑类似，但使用从v7模拟得到的结果)
K_nonpps_eq_agg = mean(kHistM_eq, 1) * paramS.Z_ss_norm_annual; % 汇总非PPS资本 (总体人均)
K_pps_eq_agg = 0; % 初始化PPS资本
if cS.pps_active && cS.pps_in_K && (cS.pps_max_contrib_frac > 0 || cS.pps_annual_contrib_limit > 0) && ~isempty(kPpsHistM_eq) 
    K_pps_eq_agg = mean(kPpsHistM_eq, 1) * paramS.Z_ss_norm_annual; % 汇总PPS资本 (总体人均)
end
Actual_K_eq_final = K_nonpps_eq_agg + K_pps_eq_agg; % 最终模拟得到的总生产性资本

C_eq_final = mean(cHistM_eq,1) * paramS.Z_ss_norm_annual; % 最终模拟得到的总消费 (总体人均)
Y_eq_final = cS.A * (Actual_K_eq_final^cS.alpha) * (paramS.L_per_capita^(1-cS.alpha)); % 最终总产出
G_eq_final = cS.gov_exp_frac_Y * Y_eq_final; % 最终政府消费
B_eq_final = cS.gov_debt_frac_Y * Y_eq_final; % 最终政府债务

% --- 打印最终均衡结果 ---
fprintf('\n--- V7 最终均衡汇总 ---\n');
fprintf('K_eq (来自替代率搜索循环): %.4f, K_eq (来自最终模拟): %.4f\n', K_eq, Actual_K_eq_final);
if abs(K_eq - Actual_K_eq_final) > 2e-2 && K_eq > 1e-9 
    % 如果差异较大，提示用户，并可以考虑是否使用最终模拟的K来重新计算Y,G,B (当前已这样做)
    warning('K_eq from rho search and K from final simulation differ significantly by %.3e. Y, G, B将使用最终模拟的K值计算。', abs(K_eq - Actual_K_eq_final));
    % Y,G,B已经使用Actual_K_eq_final计算过了
end

fprintf('均衡总生产性资本 (K*): %.4f (总体人均)\n', Actual_K_eq_final);
fprintf('  其中: 非PPS资本 K_non_pps: %.4f, PPS资本 K_pps: %.4f\n', K_nonpps_eq_agg, K_pps_eq_agg);
fprintf('均衡总劳动 (L, 效率单位, 总体人均): %.4f\n', paramS.L_per_capita);
fprintf('均衡总产出 (Y*): %.4f\n', Y_eq_final);
fprintf('均衡市场毛回报率因子 (R_mkt_gross*): %.4f (对应 r_mkt_gross*=%.4f)\n', R_mkt_gross_factor_eq_final, r_mkt_gross_eq_final);
fprintf('  家庭税后资本净回报率因子 (R_k_net_hh*): %.4f (对应 r_k_net_hh*=%.4f)\n', R_k_net_factor_hh_eq_final, r_k_net_hh_eq_final);
fprintf('均衡市场总工资率 (MPL_gross*): %.4f\n', MPL_gross_eq_final);
fprintf('目标PAYG替代率 (rho_prime_payg_eq*): %.4f (b_payg / avg_worker_gross_wage)\n', rho_prime_payg_eq);
fprintf('均衡内生实际PAYG税率 (theta_payg_eq*, 上限 %.3f): %.4f\n', cS.theta_payg_max, paramS_eq.theta_payg_actual_for_hh);
if isfield(final_eq_solution_details, 'theta_payg_required_before_cap') % 如果有理论税率信息
    fprintf('  (理论所需PAYG税率，未考虑上限前: %.4f)\n', final_eq_solution_details.theta_payg_required_before_cap);
end
fprintf('均衡内生"所得"税率 (tau_l_eq*): %.4f\n', tau_l_eq);
fprintf('  固定资本所得税率 (tau_k): %.2f, 固定消费税率 (tau_c): %.2f\n', cS.tau_k, cS.tau_c);

% 计算家庭面临的净工资率 (扣除实际PAYG税和所得税，注意所得税基数已扣减PPS)
% 净工资的准确计算应该在家庭层面考虑PPS扣减。这里显示一个近似的总体平均概念。
% w_net_hh_display = MPL_gross_eq_final * (1 - paramS_eq.theta_payg_actual_for_hh - paramS_eq.tau_l_effective_avg);
% 其中 tau_l_effective_avg 是一个平均意义上的有效所得税率（考虑PPS抵扣后）。
% 简单显示：
w_net_hh_approx_display = MPL_gross_eq_final * (1 - paramS_eq.theta_payg_actual_for_hh - paramS_eq.tau_l );
fprintf('均衡近似家庭平均净工资率 (w_gross * (1-theta_act-tau_l)): %.4f (仅为说明，未精确考虑PPS对tau_l基数的普遍影响)\n', w_net_hh_approx_display);

fprintf('均衡PAYG福利 (b_payg*, 每位退休者): %.4f\n', b_payg_eq_final);
fprintf('均衡总净转移支付 (TR_total*, 总体人均, TR_gov=0): %.4f\n', TR_total_eq_final_for_vfi);
fprintf('  其中意外遗赠 (T_bequest*): %.4f\n', T_bequest_eq_final);
fprintf('  其中政府直接转移 (TR_gov*): 0.0000 (按设定)\n');
fprintf('均衡政府消费 (G*): %.4f (G/Y* = %.3f)\n', G_eq_final, G_eq_final/Y_eq_final);
fprintf('均衡政府债务 (B*): %.4f (B/Y* = %.3f)\n', B_eq_final, B_eq_final/Y_eq_final);
fprintf('均衡 K*/Y* 比率: %.4f\n', Actual_K_eq_final / Y_eq_final );
fprintf('均衡 C*/Y* 比率: %.4f\n', C_eq_final / Y_eq_final);

% 最终核对实际达成的替代率
achieved_replacement_rate_final = 0;
if avg_worker_gross_wage_eq_final > 1e-9 % 确保分母不为零
    achieved_replacement_rate_final = b_payg_eq_final / avg_worker_gross_wage_eq_final;
end
fprintf('实际达成替代率 (b_payg / avg_worker_gross_wage): %.4f (应接近 rho_prime_payg_eq*)\n', achieved_replacement_rate_final);
if abs(achieved_replacement_rate_final - rho_prime_payg_eq) > 1e-3 && rho_prime_payg_eq > 1e-9
    warning('最终达成的替代率与目标替代率差异较大 (差异: %.3e)。请检查计算一致性。', abs(achieved_replacement_rate_final - rho_prime_payg_eq));
end

% 最终政府预算平衡检查 (使用最终模拟和均衡值)
final_gbc_residual = main_olg_v7_utils.check_gbc_residual(Actual_K_eq_final, C_eq_final, Y_eq_final, ...
    G_eq_final, B_eq_final, MPL_gross_eq_final, r_mkt_gross_eq_final, ...
    paramS_eq.theta_payg_actual_for_hh, tau_l_eq, ...
    b_payg_eq_final, T_bequest_eq_final, 0, cS, paramS_eq);
fprintf('最终GBC(一般预算)检查 @ 均衡状态: Residual = %.4e\n', final_gbc_residual);
if abs(final_gbc_residual) > 1e-2 % 如果残差仍然较大 (相对于内层循环容忍度)
    warning('最终GBC残差 (%.3e) 较大。可能需要调整内层循环的 gbc_tol_for_internal_loop 或检查一致性。', final_gbc_residual);
end


% --- 新增：绘制策略函数 (与v6类似) ---
fprintf('\n绘制最终均衡的策略函数...\n');

% 选择要绘图的参数 (与v6保持一致)
plot_a_idx = min(round(cS.aR_new / 2), cS.aD_new); % 例如，工作期中间的一个年龄组
if plot_a_idx == 0, plot_a_idx = 1; end % 避免索引为0
plot_ie_idx = round(cS.nw / 2);                 % 例如，中间的劳动效率状态

% 选择几个k_pps的网格点进行绘图，或者绘制一个二维图
plot_nkpps_to_show = min(3, cS.nkpps); % 最多显示3条k_pps的线 (如果nkpps<3则显示实际数量)
plot_ikpps_indices = [];
if cS.nkpps > 0
    plot_ikpps_indices = round(linspace(1, cS.nkpps, plot_nkpps_to_show)); % 均匀选取几个k_pps索引
else % 如果nkpps为0 (不太可能，除非参数设置如此)
    warning('nkpps为0，无法绘制k_pps相关的策略函数切片。');
end

figure_title_suffix_base = sprintf('年龄组 %d (真实年龄约 %d岁), 效率状态 %d', ...
    plot_a_idx, cS.physAgeV_new(plot_a_idx), plot_ie_idx);

if cS.nk > 1 && ~isempty(plot_ikpps_indices) % 确保非PPS资产网格和选取的PPS资产索引有效
    % --- 绘制 k_prime(k | k_pps) ---
    figure('Name', ['V7: 非PPS储蓄策略 k''(k | k_pps): ' figure_title_suffix_base]);
    hold on;
    colors = lines(plot_nkpps_to_show); % 获取一组颜色用于绘图
    legend_entries_k_prime = cell(plot_nkpps_to_show + 1, 1);
    for i_plot = 1:plot_nkpps_to_show
        ikpps = plot_ikpps_indices(i_plot); % 当前要绘制的k_pps索引
        % 从4D策略矩阵中提取对应切片 (nk x 1)
        k_prime_slice = squeeze(kPolM_eq(:, ikpps, plot_ie_idx, plot_a_idx)); 
        plot(cS.kGridV, k_prime_slice, 'LineWidth', 1.5, 'DisplayName', sprintf('k_{pps}=%.2f', cS.kppsGridV(ikpps)), 'Color', colors(i_plot,:));
        legend_entries_k_prime{i_plot} = sprintf('k_{pps}=%.2f', cS.kppsGridV(ikpps));
    end
    plot(cS.kGridV, cS.kGridV, 'k--', 'DisplayName', 'k''=k (45度线)'); % 绘制45度线作为参考
    legend_entries_k_prime{plot_nkpps_to_show + 1} = 'k''=k';
    hold off;
    xlabel('当前非PPS资产 k'); 
    ylabel('下一期非PPS资产 k''');
    title({'V7: 非PPS储蓄策略 k''(k | k_{pps})'; figure_title_suffix_base});
    legend(legend_entries_k_prime, 'Location', 'best'); grid on;

    % --- 绘制 c_pps(k | k_pps) (基于规则的) ---
    if cS.pps_active % 只有当PPS激活时才绘制
        figure('Name', ['V7: PPS缴费策略(规则) c_{pps}(k | k_pps): ' figure_title_suffix_base]);
        hold on;
        legend_entries_cpps = cell(plot_nkpps_to_show, 1);
        for i_plot = 1:plot_nkpps_to_show
            ikpps = plot_ikpps_indices(i_plot);
             % 从4D策略矩阵中提取对应切片 (nk x 1)
            cpps_slice_rule = squeeze(cPpsPolM_eq_rule(:, ikpps, plot_ie_idx, plot_a_idx));
            plot(cS.kGridV, cpps_slice_rule, 'LineWidth', 1.5, 'DisplayName', sprintf('k_{pps}=%.2f', cS.kppsGridV(ikpps)), 'Color', colors(i_plot,:));
            legend_entries_cpps{i_plot} = sprintf('k_{pps}=%.2f', cS.kppsGridV(ikpps));
        end
        hold off;
        xlabel('当前非PPS资产 k'); 
        ylabel('PPS缴费 c_{pps} (基于规则)');
        title({'V7: PPS缴费策略 (基于规则) c_{pps}(k | k_{pps})'; figure_title_suffix_base});
        legend(legend_entries_cpps, 'Location', 'best'); grid on;
    end

    % --- 绘制 c(k | k_pps) (消费量) ---
    figure('Name', ['V7: 消费策略 c(k | k_pps): ' figure_title_suffix_base]);
    hold on;
    legend_entries_c = cell(plot_nkpps_to_show, 1);
    for i_plot = 1:plot_nkpps_to_show
        ikpps = plot_ikpps_indices(i_plot);
        % 从4D策略矩阵中提取对应切片 (nk x 1)
        c_slice = squeeze(cPolM_eq(:, ikpps, plot_ie_idx, plot_a_idx)); 
        plot(cS.kGridV, c_slice, 'LineWidth', 1.5, 'DisplayName', sprintf('k_{pps}=%.2f', cS.kppsGridV(ikpps)), 'Color', colors(i_plot,:));
        legend_entries_c{i_plot} = sprintf('k_{pps}=%.2f', cS.kppsGridV(ikpps));
    end
    hold off;
    xlabel('当前非PPS资产 k'); 
    ylabel('消费量 c');
    title({'V7: 消费策略 c(k | k_{pps})'; figure_title_suffix_base});
    legend(legend_entries_c, 'Location', 'best'); grid on;

elseif cS.nk > 1 && cS.nkpps == 1 % 如果只有k变化, k_pps是标量 (例如 nkpps=1, kppsGridV只有一个点)
    figure_title_suffix_k_only = sprintf('年龄组 %d (真实年龄约 %d岁), 效率状态 %d, k_{pps}=%.2f (固定)', ...
                                        plot_a_idx, cS.physAgeV_new(plot_a_idx), plot_ie_idx, cS.kppsGridV(1));
    figure('Name', ['V7: 策略函数 (k变化, k_pps固定): ' figure_title_suffix_k_only]);
    
    subplot(1,3,1); % 子图1: k'(k)
    plot(cS.kGridV, squeeze(kPolM_eq(:, 1, plot_ie_idx, plot_a_idx)), 'b-o', 'LineWidth', 1.5); 
    hold on; plot(cS.kGridV, cS.kGridV, 'k--'); hold off;
    title('非PPS储蓄 k''(k)'); xlabel('k'); grid on; legend('k''','k''=k', 'Location','best');

    subplot(1,3,2); % 子图2: c_pps(k) (规则化)
    if cS.pps_active
        plot(cS.kGridV, squeeze(cPpsPolM_eq_rule(:, 1, plot_ie_idx, plot_a_idx)), 'r-o', 'LineWidth', 1.5); 
        title('PPS缴费 c_{pps}(k) (规则)'); xlabel('k'); grid on;
    else
        plot(cS.kGridV, zeros(cS.nk,1), 'r-o', 'LineWidth', 1.5); % 如果PPS未激活，画一条零线
        title('PPS缴费 c_{pps}(k) (未激活)'); xlabel('k'); grid on;
    end

    subplot(1,3,3); % 子图3: c(k)
    plot(cS.kGridV, squeeze(cPolM_eq(:, 1, plot_ie_idx, plot_a_idx)), 'g-o', 'LineWidth', 1.5); 
    title('消费 c(k)'); xlabel('k'); grid on;
    
    sgtitle(['V7: 策略函数切片: ' figure_title_suffix_k_only]); % 整个图的标题
else
    fprintf('无法绘制策略函数：nk或nkpps维度不足，或plot_ikpps_indices为空。\n');
    if cS.nk <= 1, fprintf('  原因: nk = %d (需要 > 1)\n', cS.nk); end
    if cS.nkpps == 0, fprintf('  原因: nkpps = 0 (需要 > 0 才能绘制多条k_pps线)\n'); end
    if isempty(plot_ikpps_indices) && cS.nkpps > 1, fprintf('  原因: plot_ikpps_indices 为空 (内部逻辑错误)\n'); end
end

fprintf('\n--- V7 (PPS 税收递延, VFI w k_pps, PPS缴费规则化) 分析完成 ---\n');
% --- END OF FILE main_olg_v7.m ---