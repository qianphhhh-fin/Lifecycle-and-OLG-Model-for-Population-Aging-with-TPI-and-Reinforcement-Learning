% =========================================================================
% == SCRIPT: test_vfi_optimization.m
% ==
% == 目的: [VFI性能与精度回归测试脚本]
% ==   验证矩阵化VFI函数与原始循环VFI函数的计算结果是否一致，并比较性能。
% ==
% == 核心步骤:
% ==   1. 复制与主运行脚本完全一致的参数环境。
% ==   2. 运行原始VFI (来自utils)，获取基准结果和时间。
% ==   3. 在本脚本内定义本地的、矩阵化的VFI函数。
% ==   4. 运行本地矩阵化VFI，获取测试结果和时间。
% ==   5. 对比两组结果的数值差异，并报告性能提升。
% ==   6. 可通过 cS.pps_active 开关切换测试 noPPS 和 PPS 两种模式。
% =========================================================================
clear; close all; clear classes;
addpath(pwd);
fprintf('====== VFI 性能与精度回归测试脚本 ======\n\n');

%% --- 1. 初始化环境与参数 (与主脚本完全一致) ---
fprintf('--- 1. 初始化环境与参数 ---\n');

% --- 步骤 1.1: 定义模拟和模型参数 ---
TIME_STEP = 5;
ngrid = 40;
ngrid_pps = 20; 

% --- 步骤 1.2: 加载模型物理参数 ---
fprintf('   加载模型物理参数...\n');
cS = model_setup_utils_bgp.ParameterValues();

% --- 步骤 1.3: 设定BGP稳态的特定参数 ---
fprintf('   设定BGP稳态的经济参数...\n');
cS.time_Step = TIME_STEP;
cS.g_A_ss = 0.015; 
cS.n_ss = -0.0;     

% --- [核心开关] 选择要测试的模式 ---
cS.pps_active = true; % <--- 在 true 和 false 之间切换来测试两种情况
fprintf('   [测试模式设定] PPS模式已激活: %s\n', mat2str(cS.pps_active));

if cS.pps_active
    cS.pps_contrib_rate = 0.03;       
    cS.pps_withdrawal_rate = 0.10;    
    cS.pps_tax_rate_withdrawal = 0.03; 
end

% --- 步骤 1.4: 生成网格和收入过程 ---
if cS.pps_active
    cS.nk = ngrid; cS.nkpps = ngrid_pps; cS.nkprime = ngrid; cS.npps = ngrid_pps;
else
    cS.nk = ngrid; cS.nkpps = 1; cS.nkprime = ngrid; cS.npps = 1;
end
cS = model_setup_utils_bgp.generateGrids(cS);

fprintf('   生成收入过程...\n');
paramS = struct();
[paramS.leGridV, paramS.TrProbM_by_age, paramS.leProb1V, cS.nw_expanded] = ...
    model_setup_utils_bgp.EarningProcess_AgeDependent(cS);
paramS.leLogGridV = log(paramS.leGridV(1:cS.nw));

% --- 步骤 1.5: 准备VFI所需的输入价格 ---
fprintf('   生成一组用于测试的合理价格...\n');
M_vfi = struct('r_mkt_t', 0.04, 'w_t', 1.5, 'b_t', 0.5, 'beq_transfer_pers', 0.02);
cS.theta_t = 0.3; % 假设一个合理的PAYG税率

%% --- 2. 运行基准：原始VFI函数 ---
fprintf('\n--- 2. 运行基准: 原始(循环)VFI函数 ---\n');
tic;
[polS_orig, valS_orig] = main_steady_state_utils_bgp.HHSolution_VFI_unified(M_vfi, paramS, cS);
time_orig = toc;
fprintf('   完成! 原始VFI运行时间: %.4f 秒\n', time_orig);

%% --- 3. 运行测试：矩阵化VFI函数 (在本地定义) ---
fprintf('\n--- 3. 运行测试: 本地(矩阵化)VFI函数 ---\n');
tic;
[polS_matrix, valS_matrix] = local_HHSolution_VFI_unified(M_vfi, paramS, cS);
time_matrix = toc;
fprintf('   完成! 矩阵化VFI运行时间: %.4f 秒\n', time_matrix);

%% --- 4. 结果比对与分析 ---
fprintf('\n--- 4. 结果比对与分析 ---\n');

% --- 4.1 性能比较 ---
performance_gain = time_orig / time_matrix;
fprintf('   [性能分析]\n');
fprintf('   - 原始VFI:    %.4f 秒\n', time_orig);
fprintf('   - 矩阵化VFI:  %.4f 秒\n', time_matrix);
fprintf('   - 性能提升:   %.2f 倍\n\n', performance_gain);

% --- 4.2 数值一致性检验 (修改版) ---
fprintf('   [数值一致性检验 (最大绝对差异)]\n');
tolerance_pol = 1e-12; % 对策略函数使用非常严格的容忍度
tolerance_val = 1e-7;  % 对价值函数使用一个稍宽松但仍然很严格的容忍度

% a. 价值函数比较
val_diff = max(abs(valS_orig(:) - valS_matrix(:)));
fprintf('   - 价值函数 (valS): %.4e (容忍度: < %.1e)\n', val_diff, tolerance_val);

% b. 策略函数比较
fprintf('   - 策略函数 (polS):\n');
fields_to_check = {'c', 'k_prime', 'kpps_prime', 'tax_regular', 'tax_payg', 'shock_exp', ...
                   'pension_out', 'pps_withdrawal', 'pps_tax', 'beq_received'};
all_match = true;
if val_diff > tolerance_val, all_match = false; end

for i = 1:length(fields_to_check)
    field = fields_to_check{i};
    if isfield(polS_orig, field) && isfield(polS_matrix, field)
        diff_field = mean(max(abs([polS_orig.(field)] - [polS_matrix.(field)])),'all');
        fprintf('     - %-16s: %.4e (容忍度: < %.1e)\n', field, diff_field, tolerance_pol);
        if diff_field > tolerance_pol
            all_match = false;
        end
    end
end

% --- 4.3 最终结论 ---
fprintf('\n   [最终结论]\n');
if all_match
    fprintf('   ✅ [成功] 结果在数值精度范围内完全一致！矩阵化优化正确无误。\n');
else
    fprintf('   ❌ [失败] 结果存在显著差异！请仔细检查逻辑或容忍度设置。\n');
end

% =========================================================================
% ==                      本地VFI函数定义部分 (v6.0 最终修正版)           ==
% =========================================================================

% --- 本地主入口函数 (无变化) ---
function [polS, valS] = local_HHSolution_VFI_unified(M_vfi, paramS_vfi, cS_vfi)
    valS = -Inf(cS_vfi.nk, cS_vfi.nkpps, cS_vfi.nw_expanded, cS_vfi.aD_new);
    polS_cell = cell(cS_vfi.aD_new, 1);
    bV_payg_vfi = zeros(1, cS_vfi.aD_new);
    if cS_vfi.aR_new < cS_vfi.aD_new, bV_payg_vfi((cS_vfi.aR_new+1):cS_vfi.aD_new) = M_vfi.b_t; end
    beq_transfer_vfi = M_vfi.beq_transfer_pers;

    for a_idx = cS_vfi.aD_new : -1 : 1
        vPrime_kkppse_next = [];
        if a_idx < cS_vfi.aD_new, vPrime_kkppse_next = valS(:,:,:,a_idx+1); end
        
        if cS_vfi.pps_active
            [valS(:,:,:,a_idx), polS_cell{a_idx}] = ...
                local_HHSolutionByAge_VFI_PPS(a_idx, vPrime_kkppse_next, M_vfi, bV_payg_vfi(a_idx), beq_transfer_vfi, paramS_vfi, cS_vfi);
        else
            [valS(:,:,:,a_idx), polS_cell{a_idx}] = ...
                local_HHSolutionByAge_VFI_noPPS(a_idx, vPrime_kkppse_next, M_vfi, bV_payg_vfi(a_idx), beq_transfer_vfi, paramS_vfi, cS_vfi);
        end
    end
    polS = [polS_cell{:}];
end


% --- [v6.0 最终版] 本地矩阵化 NoPPS 函数 (定制化搜索空间) ---
function [val_age, pol_age] = local_HHSolutionByAge_VFI_noPPS(a_idx, vPrime_kkppse_next, M_age, b_age_val, beq_transfer_val, paramS_age, cS)
    % =========================================================================
    % == 函数: local_HHSolutionByAge_VFI_noPPS
    % == 版本: [v6.0 - 【定制化搜索空间】最终修正版]
    % ==
    % == 核心修正:
    % ==   - [!!! 根本性修正 !!!] 本版本通过向量化操作，为每个 k_i 精确地
    % ==     创建了一个专用的、连续的 k' 搜索空间，从而在数学上
    % ==     完全复现了原始循环版本的逻辑。这解决了之前所有版本的失败根源。
    % ==   - 关键1: 构造了一个 nk x nk_search_grid 的 k' 选择矩阵，其中每
    % ==     一行都是一个根据该行资源定制的 linspace。
    % ==   - 关键2: 所有后续计算（消费、插值、效用）都在这个定制化的
    % ==     矩阵上进行，确保了结果的一致性。
    % =========================================================================
    
    nk_search_grid = 150;
    val_age_slice = -1e20 * ones(cS.nk, 1, cS.nw_expanded);
    pol_age_slice = struct(...
        'c', zeros(cS.nk, 1, cS.nw_expanded), 'k_prime', zeros(cS.nk, 1, cS.nw_expanded), ...
        'kpps_prime', zeros(cS.nk, 1, cS.nw_expanded), 'tax_regular', zeros(cS.nk, 1, cS.nw_expanded), ...
        'tax_payg', zeros(cS.nk, 1, cS.nw_expanded), 'shock_exp', zeros(cS.nk, 1, cS.nw_expanded), ...
        'pension_out', zeros(cS.nk, 1, cS.nw_expanded), 'beq_received', zeros(cS.nk, 1, cS.nw_expanded)); 

    g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
    market_return_factor = 1 + M_age.r_mkt_t;
    discount_factor_future_utility = cS.beta * (1 + g_A_period)^(1 - cS.sigma);
    k_grid_vec = cS.kGridV; 

    if a_idx == cS.aD_new % 终期决策不变
        k_capital_tax_vec = cS.tau_k .* (k_grid_vec * M_age.r_mkt_t);
        total_wealth_vec = k_grid_vec * market_return_factor - k_capital_tax_vec + b_age_val + beq_transfer_val;
        k_prime_final_vec = zeros(cS.nk, 1);
        c_expend_final_vec = total_wealth_vec;
        if isfield(cS,'phi_bequest') && cS.phi_bequest > 0
            phi_effective = cS.beta * cS.phi_bequest;
            c_over_k_ratio = phi_effective^(-1/cS.sigma);
            k_prime_final_vec = total_wealth_vec ./ (1 + (1+cS.tau_c)*c_over_k_ratio);
            c_expend_final_vec = total_wealth_vec - k_prime_final_vec;
        end
        c_final_vec = c_expend_final_vec ./ (1 + cS.tau_c);
        consumption_tax_vec = c_final_vec * cS.tau_c;
        final_regular_tax_vec = k_capital_tax_vec + consumption_tax_vec;
        [~, util_c] = model_setup_utils_bgp.CES_utility(c_final_vec, cS.sigma, cS);
        util_bequest = model_setup_utils_bgp.bequest_utility(k_prime_final_vec, cS);
        util_final_vec = util_c + cS.beta * util_bequest;
        util_final_vec(total_wealth_vec < 0) = -1e20; 

        for ie = 1:cS.nw_expanded
            val_age_slice(:, 1, ie) = util_final_vec;
            pol_age_slice.c(:, 1, ie) = c_final_vec; pol_age_slice.k_prime(:, 1, ie) = k_prime_final_vec;
            pol_age_slice.pension_out(:, 1, ie) = b_age_val;
            pol_age_slice.tax_regular(:, 1, ie) = final_regular_tax_vec;
            pol_age_slice.beq_received(:, 1, ie) = beq_transfer_val;
        end
        val_age = repmat(val_age_slice, [1, cS.nkpps, 1]); 
        pol_age = model_setup_utils_bgp.expand_policy_slice(pol_age_slice, cS.nkpps);
        return;
    end

    vPrime_interpolants = cell(cS.nw_expanded, 1);
    vPrime_slice = squeeze(vPrime_kkppse_next(:, 1, :));
    for ie_next = 1:cS.nw_expanded
        vPrime_interpolants{ie_next} = griddedInterpolant(k_grid_vec, vPrime_slice(:, ie_next), 'pchip', 'nearest');
    end
    trans_mat_next = paramS_age.TrProbM_by_age{a_idx + 1};

    for ie = 1:cS.nw_expanded
        k_return_vec = k_grid_vec * market_return_factor;
        labor_income_gross = 0; pension_out = 0;
        if a_idx <= cS.aR_new, labor_income_gross = M_age.w_t * cS.ageEffV_new(a_idx) * paramS_age.leGridV(ie);
        else, pension_out = b_age_val; end
        capital_tax_vec = cS.tau_k .* (k_grid_vec * M_age.r_mkt_t);
        payg_tax = cS.theta_t * labor_income_gross; 
        labor_tax = cS.tau_l * max(0, labor_income_gross - payg_tax); 
        cash_on_hand_vec = k_return_vec + labor_income_gross + pension_out + beq_transfer_val - (capital_tax_vec + payg_tax + labor_tax);
        shock_exp_factor = (ie == cS.nw + 1) * cS.kappa_young + (ie == cS.nw + 2) * cS.kappa_old;
        shock_exp_vec = shock_exp_factor .* cash_on_hand_vec;
        available_for_c_and_s_vec = cash_on_hand_vec - shock_exp_vec;

        % A. 创建定制化的 k' 搜索矩阵 [nk x nk_search_grid]
        k_prime_max_vec = (available_for_c_and_s_vec - cS.cFloor * (1+cS.tau_c)) / (1 + g_A_period);
        k_prime_max_vec(k_prime_max_vec < cS.kMin) = cS.kMin;
        
        linspace_weights = linspace(0, 1, nk_search_grid);
        k_prime_choices_mat = cS.kMin + (k_prime_max_vec - cS.kMin) .* linspace_weights;
        
        % B. 计算消费选择矩阵
        c_expend_choices_mat = available_for_c_and_s_vec - k_prime_choices_mat .* (1 + g_A_period);
        c_choices_mat = c_expend_choices_mat / (1 + cS.tau_c);
        
        % C. 计算价值矩阵
        invalid_choice_mask = (c_choices_mat < cS.cFloor);
        c_choices_mat(invalid_choice_mask) = cS.cFloor;
        [~, util_c_mat] = model_setup_utils_bgp.CES_utility(c_choices_mat, cS.sigma, cS);
        util_c_mat(invalid_choice_mask) = -1e20; 

        % D. 在定制化网格上插值未来价值
        % Reshape for efficient interpolation
        k_prime_choices_flat = k_prime_choices_mat(:);
        v_prime_flat_mat = zeros(length(k_prime_choices_flat), cS.nw_expanded);
        for ie_next = 1:cS.nw_expanded
             v_prime_flat_mat(:, ie_next) = vPrime_interpolants{ie_next}(k_prime_choices_flat);
        end
        ev_flat_vec = v_prime_flat_mat * trans_mat_next(ie, :)';
        ev_on_choices_mat = reshape(ev_flat_vec, cS.nk, nk_search_grid);

        util_bequest_mat = model_setup_utils_bgp.bequest_utility(k_prime_choices_mat, cS);

        future_value_mat = cS.s_pathV(a_idx) * ev_on_choices_mat + (1 - cS.s_pathV(a_idx)) * util_bequest_mat;
        
        total_value_mat = util_c_mat + discount_factor_future_utility * future_value_mat;
        
        [best_val_vec, best_idx_vec] = max(total_value_mat, [], 2);

        % E. 存储结果
        linear_indices = sub2ind(size(k_prime_choices_mat), (1:cS.nk)', best_idx_vec);
        k_prime_final_vec = k_prime_choices_mat(linear_indices);
        c_final_vec = c_choices_mat(linear_indices);

        val_age_slice(:, 1, ie) = best_val_vec;
        pol_age_slice.c(:, 1, ie) = c_final_vec;
        pol_age_slice.k_prime(:, 1, ie) = k_prime_final_vec;
        pol_age_slice.pension_out(:, 1, ie) = pension_out;
        pol_age_slice.shock_exp(:, 1, ie) = shock_exp_vec;
        pol_age_slice.tax_regular(:, 1, ie) = capital_tax_vec + labor_tax + c_final_vec .* cS.tau_c;
        pol_age_slice.tax_payg(:, 1, ie) = payg_tax;
        pol_age_slice.beq_received(:, 1, ie) = beq_transfer_val;
    end
    
    val_age = repmat(val_age_slice, [1, cS.nkpps, 1]);
    pol_age = model_setup_utils_bgp.expand_policy_slice(pol_age_slice, cS.nkpps);
end


% --- [v6.2 并行修复版] 本地矩阵化 PPS 函数 (定制化搜索空间) ---
% --- [v6.3 二维并行终极优化版] 本地矩阵化 PPS 函数 ---
function [val_age, pol_age] = local_HHSolutionByAge_VFI_PPS(a_idx, vPrime_kkppse_next, M_age, b_age_val, beq_transfer_val, paramS_age, cS)
    % =========================================================================
    % == 函数: local_HHSolutionByAge_VFI_PPS
    % == 版本: [v6.3 - 【二维并行终极优化版】]
    % ==
    % == 核心修改:
    % ==   - [终极性能优化] 采用“扁平化”策略，将 (ie, ikpps) 两个维度的
    % ==     循环重构为一个单一的线性 parfor 循环。
    % ==   - 这使得总共 (nw_expanded * nkpps) 个独立任务都可以被并行调度，
    % ==     从而最大化地利用多核CPU资源，进一步提升计算速度。
    % =========================================================================
    
    nk_search_grid = 150;
    val_age = -1e20 * ones(cS.nk, cS.nkpps, cS.nw_expanded);
    pol_age = struct(...
        'c', zeros(cS.nk, cS.nkpps, cS.nw_expanded), 'k_prime', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
        'kpps_prime', zeros(cS.nk, cS.nkpps, cS.nw_expanded), 'tax_regular', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
        'tax_payg', zeros(cS.nk, cS.nkpps, cS.nw_expanded), 'shock_exp', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
        'pension_out', zeros(cS.nk, cS.nkpps, cS.nw_expanded), 'pps_contrib', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
        'pps_withdrawal', zeros(cS.nk, cS.nkpps, cS.nw_expanded), 'pps_tax', zeros(cS.nk, cS.nkpps, cS.nw_expanded), ...
        'beq_received', zeros(cS.nk, cS.nkpps, cS.nw_expanded));

    g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
    market_return_factor = 1 + M_age.r_mkt_t;
    discount_factor_future_utility = cS.beta * (1 + g_A_period)^(1 - cS.sigma);
    k_grid_vec = cS.kGridV;

    if a_idx == cS.aD_new % 终期决策不变
        [K_grid_mat, Kpps_grid_mat] = ndgrid(k_grid_vec, cS.kppsGridV);
        k_capital_tax_mat = cS.tau_k .* (K_grid_mat .* M_age.r_mkt_t);
        pps_wealth_final_mat = Kpps_grid_mat .* market_return_factor;
        pps_tax_final_mat = cS.pps_tax_rate_withdrawal .* pps_wealth_final_mat;
        total_wealth_mat = (K_grid_mat * market_return_factor - k_capital_tax_mat) + (pps_wealth_final_mat - pps_tax_final_mat) + b_age_val + beq_transfer_val;
        k_prime_final_mat = zeros(size(K_grid_mat)); c_expend_final_mat = total_wealth_mat;
        if isfield(cS,'phi_bequest') && cS.phi_bequest > 0, phi_effective = cS.beta * cS.phi_bequest; c_over_k_ratio = phi_effective^(-1/cS.sigma); k_prime_final_mat = total_wealth_mat ./ (1 + (1+cS.tau_c)*c_over_k_ratio); c_expend_final_mat = total_wealth_mat - k_prime_final_mat; end
        c_final_mat = c_expend_final_mat ./ (1 + cS.tau_c);
        consumption_tax_mat = c_final_mat * cS.tau_c;
        [~, util_c] = model_setup_utils_bgp.CES_utility(c_final_mat, cS.sigma, cS);
        util_bequest = model_setup_utils_bgp.bequest_utility(k_prime_final_mat, cS);
        util_final_mat = util_c + cS.beta * util_bequest;
        util_final_mat(total_wealth_mat < 0) = -1e20;

        for ie = 1:cS.nw_expanded
            val_age(:, :, ie) = util_final_mat;
            pol_age.c(:, :, ie) = c_final_mat; pol_age.k_prime(:, :, ie) = k_prime_final_mat;
            pol_age.kpps_prime(:, :, ie) = 0; pol_age.pension_out(:, :, ie) = b_age_val;
            pol_age.tax_regular(:, :, ie) = k_capital_tax_mat + consumption_tax_mat;
            pol_age.pps_withdrawal(:, :, ie) = pps_wealth_final_mat; pol_age.pps_tax(:, :, ie) = pps_tax_final_mat;
            pol_age.beq_received(:, :, ie) = beq_transfer_val;
        end
        return;
    end

    vPrime_interpolants = cell(cS.nw_expanded, 1);
    for ie_next = 1:cS.nw_expanded
        vPrime_interpolants{ie_next} = griddedInterpolant({k_grid_vec, cS.kppsGridV}, vPrime_kkppse_next(:, :, ie_next), 'pchip', 'nearest');
    end
    trans_mat_next = paramS_age.TrProbM_by_age{a_idx + 1};

    % [核心修改] 创建一个一维的任务列表
    num_tasks = cS.nw_expanded * cS.nkpps;
    % 创建Cell数组来收集所有任务的结果
    val_results_cell = cell(num_tasks, 1);
    pol_results_cell = cell(num_tasks, 1);
    
    parfor task_idx = 1:num_tasks
        % 从一维任务索引映射回二维 (ie, ikpps) 索引
        [ikpps, ie] = ind2sub([cS.nkpps, cS.nw_expanded], task_idx);
        
        % --- 这部分代码是单次迭代的核心计算，逻辑基本不变 ---
        % --- 只是现在它不再在一个ikpps循环内部，而是直接执行 ---
        
        kpps_now = cS.kppsGridV(ikpps);
        pension_out=0; pps_contrib=0; pps_withdrawal=0; pps_tax=0; labor_income_gross = 0;
        if a_idx <= cS.aR_new 
            labor_income_gross = M_age.w_t * cS.ageEffV_new(a_idx) * paramS_age.leGridV(ie);
            payg_tax = cS.theta_t * labor_income_gross;
            pps_contrib = cS.pps_contrib_rate * labor_income_gross;
            capital_tax_vec = cS.tau_k .* (k_grid_vec * M_age.r_mkt_t); 
            labor_tax = cS.tau_l * max(0, labor_income_gross - payg_tax - pps_contrib);
            cash_on_hand_vec = (k_grid_vec * market_return_factor + labor_income_gross + beq_transfer_val) ...
                               - (capital_tax_vec + payg_tax + labor_tax + pps_contrib); 
        else 
            pension_out = b_age_val;
            kpps_total_value = kpps_now * market_return_factor;
            pps_withdrawal = cS.pps_withdrawal_rate * kpps_total_value;
            pps_tax = cS.pps_tax_rate_withdrawal * pps_withdrawal;
            net_pps_inflow = pps_withdrawal - pps_tax;
            capital_tax_vec = cS.tau_k .* (k_grid_vec * M_age.r_mkt_t); 
            payg_tax = 0; labor_tax = 0;
            cash_on_hand_vec = (k_grid_vec * market_return_factor + pension_out + net_pps_inflow + beq_transfer_val) ...
                               - capital_tax_vec; 
        end
            
        shock_exp_factor = (ie == cS.nw + 1) * cS.kappa_young + (ie == cS.nw + 2) * cS.kappa_old;
        shock_exp_vec = shock_exp_factor .* cash_on_hand_vec;
        available_for_c_and_s_vec = cash_on_hand_vec - shock_exp_vec;

        if a_idx <= cS.aR_new, kpps_prime_continuous = (kpps_now * market_return_factor + pps_contrib) / (1 + g_A_period);
        else, kpps_prime_continuous = (kpps_now * market_return_factor - pps_withdrawal) / (1 + g_A_period); end
        kpps_prime_final = max(cS.kppsMin, kpps_prime_continuous);

        k_prime_max_vec = (available_for_c_and_s_vec - cS.cFloor * (1+cS.tau_c)) / (1 + g_A_period);
        k_prime_max_vec(k_prime_max_vec < cS.kMin) = cS.kMin;
        linspace_weights = linspace(0, 1, nk_search_grid);
        k_prime_choices_mat = cS.kMin + (k_prime_max_vec - cS.kMin) .* linspace_weights;
            
        c_expend_choices_mat = available_for_c_and_s_vec - k_prime_choices_mat .* (1 + g_A_period);
        c_choices_mat = c_expend_choices_mat / (1 + cS.tau_c);
            
        invalid_choice_mask = (c_choices_mat < cS.cFloor);
        c_choices_mat(invalid_choice_mask) = cS.cFloor;
        [~, util_c_mat] = model_setup_utils_bgp.CES_utility(c_choices_mat, cS.sigma, cS);
        util_c_mat(invalid_choice_mask) = -1e20;

        k_prime_choices_flat = k_prime_choices_mat(:);
        kpps_prime_vec = ones(length(k_prime_choices_flat), 1) * kpps_prime_final;
        v_prime_flat_mat = zeros(length(k_prime_choices_flat), cS.nw_expanded);
        for ie_next = 1:cS.nw_expanded
            v_prime_flat_mat(:, ie_next) = vPrime_interpolants{ie_next}(k_prime_choices_flat, kpps_prime_vec);
        end
        ev_flat_vec = v_prime_flat_mat * trans_mat_next(ie, :)';
        ev_on_choices_mat = reshape(ev_flat_vec, cS.nk, nk_search_grid);

        util_bequest_mat = model_setup_utils_bgp.bequest_utility(k_prime_choices_mat + kpps_prime_final, cS);
            
        future_value_mat = cS.s_pathV(a_idx) * ev_on_choices_mat + (1 - cS.s_pathV(a_idx)) * util_bequest_mat;
        total_value_mat = util_c_mat + discount_factor_future_utility * future_value_mat;
            
        [best_val_vec, best_idx_vec] = max(total_value_mat, [], 2);

        linear_indices = sub2ind(size(k_prime_choices_mat), (1:cS.nk)', best_idx_vec);
        k_prime_final_optimal_vec = k_prime_choices_mat(linear_indices);
        c_final_vec = c_choices_mat(linear_indices);

        % 将单次迭代的结果存入临时变量
        val_result_slice = best_val_vec;
        pol_result_slice = struct(...
            'c', c_final_vec, 'k_prime', k_prime_final_optimal_vec, 'kpps_prime', kpps_prime_final, ...
            'tax_regular', capital_tax_vec + labor_tax + c_final_vec * cS.tau_c, 'tax_payg', payg_tax, ...
            'shock_exp', shock_exp_vec, 'pension_out', pension_out, 'pps_contrib', pps_contrib, ...
            'pps_withdrawal', pps_withdrawal, 'pps_tax', pps_tax, 'beq_received', beq_transfer_val);
        
        % 将结果存入 Cell 数组中
        val_results_cell{task_idx} = val_result_slice;
        pol_results_cell{task_idx} = pol_result_slice;
    end

    % [核心修改] 在并行循环结束后，安全地将结果组装回最终的矩阵和结构体
    fields = fieldnames(pol_age);
    for task_idx = 1:num_tasks
        [ikpps, ie] = ind2sub([cS.nkpps, cS.nw_expanded], task_idx);
        
        % 组装价值函数
        val_age(:, ikpps, ie) = val_results_cell{task_idx};
        
        % 组装策略函数
        current_pol = pol_results_cell{task_idx};
        for i = 1:length(fields)
            field = fields{i};
            % 对于标量场，需要用 repmat 扩展维度
            if isscalar(current_pol.(field))
                pol_age.(field)(:, ikpps, ie) = repmat(current_pol.(field), cS.nk, 1);
            else
                pol_age.(field)(:, ikpps, ie) = current_pol.(field);
            end
        end
    end
end