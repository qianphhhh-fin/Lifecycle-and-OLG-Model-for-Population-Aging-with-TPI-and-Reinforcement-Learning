% --- START OF FILE main_olg_v8_nn.m ---

% OLG 模型 V8 (内生PPS缴费决策, PPS所得税递延, VFI w k_pps状态):
% 目标: 求解给定 PAYG 替代率 (rho_prime_payg_fixed) 下的均衡
% PPS缴费: 个体优化选择PPS缴费额，但受收入比例上限和年度绝对上限约束。
% 其他特性同Baseline:
%   - PPS 缴费可从所得税前扣除 (所得税率为 tau_l)。
%   - tau_l 内生调整以平衡政府一般预算 (TR_gov = 0)。
%   - PAYG 税率 (theta_payg) 内生决定，但有上限 cS.theta_payg_max。
%   - VFI 状态变量仍然包含 k_pps (PPS资产)。
%   - 新增: 可选择使用预训练的神经网络 (NN) 替代VFI内层循环 (分阶段模型)。

clc; % 清除命令行窗口
clear; % 清除工作区变量
close all; % 关闭所有图形窗口
fprintf('=== OLG 模型 V8 (内生PPS缴费, 固定 Rho_prime_payg, VFI w k_pps) - NN分阶段版 ===\n');
fprintf('    (Rho_prime_payg 固定, TR_gov=0, tau_l 内生, theta_payg 有上限)\n');
fprintf('    (VFI 状态: k, k_pps, eps; PPS缴费: 内生选择，有比例和绝对上限)\n');
fprintf('    (可选择使用NN分阶段模型进行VFI近似)\n');
% 关闭插值相关警告
warning('off', 'MATLAB:griddedInterpolant:MeshgridFailPointWarnId');
warning('off', 'MATLAB:griddedInterpolant:InterpEmptyGridId');
warning('off', 'MATLAB:griddedInterpolant:DegenerateGridId');


%% 1. 固定模型参数与人口设置
fprintf('\n--- 1. 初始化参数 ---\n');
cS = main_olg_v8_utils.ParameterValues_HuggettStyle(); % 调用v8的参数设置函数
paramS = struct(); % 初始化一个结构体用于存储派生参数

% *** 神经网络模型加载 (如果使用NN替代VFI) ***
nn_models_loaded = struct();
if cS.use_NN_for_VFI
    fprintf('\n--- 加载神经网络模型 (替代VFI) ---\n');
    models_loaded_successfully = true;
    
    % 首先尝试加载新的三网络文件
    if exist(cS.nn_three_networks_filename, 'file')
        fprintf('尝试加载三网络模型文件: %s\n', cS.nn_three_networks_filename);
        try
            data_three_nets = load(cS.nn_three_networks_filename);
            if isfield(data_three_nets, 'trainedNet_pps_classifier') && ...
               isfield(data_three_nets, 'trainedNet_k_prime_regressor') && ...
               isfield(data_three_nets, 'total_input_dim') && ...
               isfield(data_three_nets, 'use_b_payg_as_input')
                
                % 成功加载三网络模型
                nn_models_loaded.trainedNet_pps_classifier = data_three_nets.trainedNet_pps_classifier;
                nn_models_loaded.trainedNet_k_prime_regressor = data_three_nets.trainedNet_k_prime_regressor;
                nn_models_loaded.total_input_dim = data_three_nets.total_input_dim;
                nn_models_loaded.use_b_payg_as_input = data_three_nets.use_b_payg_as_input;
                
                % PPS回归器可能为空
                if isfield(data_three_nets, 'trainedNet_pps_regressor')
                    nn_models_loaded.trainedNet_pps_regressor = data_three_nets.trainedNet_pps_regressor;
                else
                    nn_models_loaded.trainedNet_pps_regressor = [];
                end
                
                % 加载其他可能的参数
                if isfield(data_three_nets, 'cS')
                    nn_models_loaded.cS = data_three_nets.cS;
                end
                
                fprintf('三网络模型加载成功！\n');
                fprintf('  - PPS分类器: 已加载\n');
                % fprintf('  - PPS回归器: %s\n', isempty(nn_models_loaded.trainedNet_pps_regressor) ? '未提供' : '已加载');
                fprintf('  - k_prime回归器: 已加载\n');
                fprintf('  - 输入维度: %d\n', nn_models_loaded.total_input_dim);
                % fprintf('  - 使用b_payg作为输入: %s\n', nn_models_loaded.use_b_payg_as_input ? 'true' : 'false');
                
                % 验证参数一致性
                if isfield(nn_models_loaded, 'cS')
                    nn_cs_check = nn_models_loaded.cS;
                    mismatch_nn = false;
                    if nn_cs_check.nk ~= cS.nk, fprintf('警告: NN nk (%d) 与当前 cS.nk (%d) 不匹配。\n', nn_cs_check.nk, cS.nk); mismatch_nn=true; end
                    if nn_cs_check.nkpps ~= cS.nkpps, fprintf('警告: NN nkpps (%d) 与当前 cS.nkpps (%d) 不匹配。\n', nn_cs_check.nkpps, cS.nkpps); mismatch_nn=true; end
                    if nn_cs_check.nw ~= cS.nw, fprintf('警告: NN nw (%d) 与当前 cS.nw (%d) 不匹配。\n', nn_cs_check.nw, cS.nw); mismatch_nn=true; end
                    if nn_cs_check.aD_new ~= cS.aD_new, fprintf('警告: NN aD_new (%d) 与当前 cS.aD_new (%d) 不匹配。\n', nn_cs_check.aD_new, cS.aD_new); mismatch_nn=true; end
                    if mismatch_nn
                        warning('加载的NN的参数对应的cS网格设置与当前cS不完全匹配。结果可能不准确或出错。');
                    end
                end
                
            else
                fprintf('三网络文件格式不正确，回退到分离文件加载...\n');
                models_loaded_successfully = false;
            end
        catch ME
            fprintf('三网络文件加载失败: %s\n回退到分离文件加载...\n', ME.message);
            models_loaded_successfully = false;
        end
    else
        fprintf('三网络文件不存在，尝试分离文件加载...\n');
        models_loaded_successfully = false;
    end
    
    % 如果三网络加载失败，回退到分离文件方式
    if ~models_loaded_successfully
        fprintf('使用分离文件方式加载神经网络...\n');
        models_loaded_successfully = true; % 重置标志
        
        % 加载 k_prime 模型
        if exist(cS.nn_model_k_prime_filename, 'file')
            fprintf('正在加载 k_prime 网络和参数从: %s\n', cS.nn_model_k_prime_filename);
            data_kp = load(cS.nn_model_k_prime_filename);
            if isfield(data_kp, 'trainedNet_k_prime') && isfield(data_kp, 'input_means_all')
                nn_models_loaded.k_prime_net = data_kp.trainedNet_k_prime;
                nn_models_loaded.k_prime_norm_params = data_kp;
            else
                warning('k_prime NN MAT文件 %s 缺少关键变量。', cS.nn_model_k_prime_filename);
                models_loaded_successfully = false;
            end
        else
            warning('找不到 k_prime NN MAT文件: %s', cS.nn_model_k_prime_filename);
            models_loaded_successfully = false;
        end

        % 加载 c_pps 分类器模型
        if models_loaded_successfully && exist(cS.nn_model_c_pps_classifier_filename, 'file')
            fprintf('正在加载 c_pps 分类器网络和参数从: %s\n', cS.nn_model_c_pps_classifier_filename);
            data_clf = load(cS.nn_model_c_pps_classifier_filename);
            if isfield(data_clf, 'trainedNet_classifier_pps') && isfield(data_clf, 'input_means_all')
                nn_models_loaded.c_pps_classifier_net = data_clf.trainedNet_classifier_pps;
                nn_models_loaded.c_pps_classifier_norm_params = data_clf; % 假设共享输入归一化参数 (来自k_prime或clf文件)
            else
                warning('c_pps 分类器 NN MAT文件 %s 缺少关键变量。', cS.nn_model_c_pps_classifier_filename);
                models_loaded_successfully = false;
            end
        else
            if models_loaded_successfully, warning('找不到 c_pps 分类器 NN MAT文件: %s', cS.nn_model_c_pps_classifier_filename); end
            models_loaded_successfully = false;
        end

        % 加载 c_pps 回归器模型
        if models_loaded_successfully && exist(cS.nn_model_c_pps_regressor_filename, 'file')
            fprintf('正在加载 c_pps 回归器网络和参数从: %s\n', cS.nn_model_c_pps_regressor_filename);
            data_reg = load(cS.nn_model_c_pps_regressor_filename);
            if isfield(data_reg, 'trainedNet_regressor_pps') && ~isempty(data_reg.trainedNet_regressor_pps) && ...
               isfield(data_reg, 'input_means_all') && isfield(data_reg, 'output_mean_c_pps_reg')
                nn_models_loaded.c_pps_regressor_net = data_reg.trainedNet_regressor_pps;
                nn_models_loaded.c_pps_regressor_norm_params = data_reg;
            else
                fprintf('信息: c_pps 回归器 NN MAT文件 %s 中模型为空或缺少参数。如果分类器预测缴费，c_pps金额将为0。\n', cS.nn_model_c_pps_regressor_filename);
                % nn_models_loaded.c_pps_regressor_net 保持为空
            end
        else
            if models_loaded_successfully, fprintf('信息: 未找到 c_pps 回归器 NN MAT文件: %s (可能是因为训练时c_pps>0样本不足)。\n', cS.nn_model_c_pps_regressor_filename); end
        end

        if models_loaded_successfully && ~isempty(nn_models_loaded.k_prime_net) && ~isempty(nn_models_loaded.c_pps_classifier_net)
            fprintf('神经网络模型加载完毕 (k_prime 和 c_pps分类器已加载；c_pps回归器可选)。\n');
            % 这里可以添加 cS 与 norm_params.cS 的验证逻辑
            if isfield(nn_models_loaded.k_prime_norm_params, 'cS')
                nn_cs_check = nn_models_loaded.k_prime_norm_params.cS;
                 mismatch_nn = false;
                if nn_cs_check.nk ~= cS.nk, fprintf('警告: NN nk (%d) 与当前 cS.nk (%d) 不匹配。\n', nn_cs_check.nk, cS.nk); mismatch_nn=true; end
                if nn_cs_check.nkpps ~= cS.nkpps, fprintf('警告: NN nkpps (%d) 与当前 cS.nkpps (%d) 不匹配。\n', nn_cs_check.nkpps, cS.nkpps); mismatch_nn=true; end
                if nn_cs_check.nw ~= cS.nw, fprintf('警告: NN nw (%d) 与当前 cS.nw (%d) 不匹配。\n', nn_cs_check.nw, cS.nw); mismatch_nn=true; end
                if nn_cs_check.aD_new ~= cS.aD_new, fprintf('警告: NN aD_new (%d) 与当前 cS.aD_new (%d) 不匹配。\n', nn_cs_check.aD_new, cS.aD_new); mismatch_nn=true; end
                if mismatch_nn
                    warning('加载的NN的参数对应的cS网格设置与当前cS不完全匹配。结果可能不准确或出错。');
                end
            end
        else
            warning('一个或多个必要的NN模型加载失败。将回退到完整VFI。');
            cS.use_NN_for_VFI = false;
        end
    end
    
    if ~models_loaded_successfully
        warning('所有神经网络加载方式都失败。将回退到完整VFI。');
        cS.use_NN_for_VFI = false;
    end
else
    fprintf('cS.use_NN_for_VFI 为 false，将不使用NN。\n');
end
% --------------------------------------

% *** Fixed PAYG Replacement Rate (as in Baseline) ***
cS.rho_prime_payg_fixed = 0.20; 
fprintf('>>> V8: 固定 PAYG 替代率 (rho_prime_payg_fixed): %.3f\n', cS.rho_prime_payg_fixed);

fprintf('参数已加载。nk=%d, nkpps=%d, nw=%d, nPpsChoiceGrid=%d。\n', cS.nk, cS.nkpps, cS.nw, cS.n_pps_choice_grid_points); 
fprintf('年度年龄范围: %d-%d。模型年龄组数: %d。\n', cS.age1_orig, cS.ageLast_orig, cS.aD_new); 
fprintf('固定税率: tau_k=%.2f, tau_c=%.2f。G/Y=%.2f, B/Y=%.2f。\n', cS.tau_k, cS.tau_c, cS.gov_exp_frac_Y, cS.gov_debt_frac_Y); 
fprintf('PAYG 税率上限 (theta_payg_max): %.3f\n', cS.theta_payg_max); 
fprintf('所得税率 tau_l 范围: [%.3f, %.3f], 总劳动税上限: %.3f\n', cS.tau_l_min, cS.tau_l_max, cS.max_total_labor_tax); 
fprintf('PPS 年度缴费上限 (绝对值): %.2f, 比例上限 (法定): %.2f\n', cS.pps_annual_contrib_limit, cS.pps_max_contrib_frac); 
fprintf('是否使用NN进行VFI近似: %s\n', mat2str(cS.use_NN_for_VFI)); 


%% 2. 模拟人口动态至稳态
% ... (此部分与你之前版本相同，此处省略以节省空间) ...
fprintf('\n--- 2. 模拟人口动态 ---\n');
popS = main_olg_v8_utils.initPopulation(cS); 
popS = main_olg_v8_utils.populationDynamics(popS, cS); 
[Z_ss, ~, ~, ~] = main_olg_v8_utils.detectSteadyStatePopulation(popS, cS); 
paramS.Z_ss_counts = Z_ss; 
Z_ss_total = sum(Z_ss); 
Z_ss_norm_group = zeros(cS.aD_new,1); 
if Z_ss_total > 1e-9, Z_ss_norm_group = Z_ss / Z_ss_total; end
paramS.ageMassV = Z_ss_norm_group(:); 
paramS.mass_workers_group = sum(paramS.ageMassV(1:cS.aR_new)); 
if paramS.mass_workers_group < 1e-9, error('模型中稳态工作人口占比为零或过小。'); end
Z_ss_norm_annual = zeros(cS.aD_orig,1);
if Z_ss_total > 1e-9
    for a_new_map_idx = 1:cS.aD_new
        annual_indices_in_group = cS.physAgeMap{a_new_map_idx};
        group_mass_fraction = Z_ss_norm_group(a_new_map_idx);
        num_years_in_this_group = length(annual_indices_in_group);
        if num_years_in_this_group > 0
            mass_per_year_in_group = group_mass_fraction / num_years_in_this_group;
            Z_ss_norm_annual(annual_indices_in_group) = mass_per_year_in_group;
        end
    end
    if sum(Z_ss_norm_annual) > 1e-9 && abs(sum(Z_ss_norm_annual) - 1.0) > 1e-6
        Z_ss_norm_annual = Z_ss_norm_annual / sum(Z_ss_norm_annual);
    elseif sum(Z_ss_norm_annual) < 1e-9, Z_ss_norm_annual(:) = 1/cS.aD_orig; end
else, Z_ss_norm_annual(:) = 1/cS.aD_orig; end
paramS.Z_ss_norm_annual = Z_ss_norm_annual;
if Z_ss_total > 1e-9 && length(popS.totalPop) > 1 && popS.totalPop(end-1) > 1e-9
    pop_growth_factor_per_group_period = popS.totalPop(end) / popS.totalPop(end-1);
    paramS.popGrowthForDebt = pop_growth_factor_per_group_period^(1/cS.yearStep) - 1;
else, paramS.popGrowthForDebt = cS.popGrowth_orig; end
fprintf('人口参数计算完毕。年化稳态人口增长率 (用于债务): %.4f\n', paramS.popGrowthForDebt);
fprintf('稳态工人占比 (基于年龄组人口): %.4f\n', paramS.mass_workers_group);

%% 3. 预计算劳动供给和禀赋过程
% ... (此部分与你之前版本相同，此处省略) ...
fprintf('\n--- 3. 预计算劳动 ---\n');
[paramS.leLogGridV, paramS.leTrProbM, paramS.leProb1V] = main_olg_v8_utils.EarningProcess_olgm(cS);
paramS.leGridV = exp(paramS.leLogGridV(:)); 
paramS.ageEffV_new = cS.ageEffV_new;
eIdxM = main_olg_v8_utils.LaborEndowSimulation_olgm(cS, paramS);
[~, L_per_capita] = main_olg_v8_utils.LaborSupply_Huggett(eIdxM, cS, paramS, paramS.ageMassV);
fprintf('总劳动供给 (L, 效率单位, 总体人均): %.4f\n', L_per_capita);
if L_per_capita <= 0, L_per_capita = 1e-6; warning('L_per_capita 为零或负，已重置为1e-6。'); end
paramS.L_per_capita = L_per_capita;
[leLogGridV_check, ~, leProb1V_check] = main_olg_v8_utils.EarningProcess_olgm(cS); % Renamed for clarity
leGridV_check_vals = exp(leLogGridV_check(:));
fprintf('劳动效率状态网格值：\n');
for i = 1:length(leGridV_check_vals), fprintf('状态 %d: %.4f\n', i, leGridV_check_vals(i)); end
fprintf('\n初始分布概率：\n');
for i = 1:length(leProb1V_check), fprintf('状态 %d: %.4f\n', i, leProb1V_check(i)); end


%% 4. 求解一般均衡 (给定固定的 rho_prime_payg_fixed)
fprintf('\n--- 4. 求解一般均衡 (固定 rho_prime_payg_fixed=%.3f) ---\n', cS.rho_prime_payg_fixed);
fprintf('  当前格点参数: n_k=%d, n_kpps=%d, n_w=%d, n_PpsChoiceGrid=%d\n', cS.nk, cS.nkpps, cS.nw, cS.n_pps_choice_grid_points);

K_global_guess =4.5 ; 
paramS_for_solver = paramS; % paramS 包含了 leGridV, Z_ss_norm_annual 等
paramS_for_solver.suppress_inner_print_header = false;
paramS_for_solver.suppress_initial_theta_print = false;

fprintf('调用均衡求解器 solve_K_tau_l_for_rho_prime_nn_twostage (V8) with fixed rho_prime_payg_fixed...\n');
vfi_mode_solver_str_ge = '完整VFI';

% 检查神经网络是否已正确加载（兼容两种加载方式）
nn_available = false;
if cS.use_NN_for_VFI
    % 检查三网络方式的字段
    if isfield(nn_models_loaded, 'trainedNet_pps_classifier') && ~isempty(nn_models_loaded.trainedNet_pps_classifier) && ...
       isfield(nn_models_loaded, 'trainedNet_k_prime_regressor') && ~isempty(nn_models_loaded.trainedNet_k_prime_regressor)
        nn_available = true;
        vfi_mode_solver_str_ge = '预训练NN (三网络)';
    % 检查分离文件方式的字段
    elseif isfield(nn_models_loaded, 'k_prime_net') && ~isempty(nn_models_loaded.k_prime_net) && ...
           isfield(nn_models_loaded, 'c_pps_classifier_net') && ~isempty(nn_models_loaded.c_pps_classifier_net)
        nn_available = true;
        vfi_mode_solver_str_ge = '预训练NN (分阶段)';
    end
end

if ~nn_available && cS.use_NN_for_VFI
    fprintf('    (警告: cS.use_NN_for_VFI为真但必要NN未加载，将回退到完整VFI)\n');
end
fprintf('    (模式: 使用 %s 进行家庭问题求解)\n', vfi_mode_solver_str_ge);

solve_start_time = tic;
[K_eq, tau_l_eq, gbc_residual_eq, eq_found, final_eq_solution_details] = ...
    main_olg_v8_utils.solve_K_tau_l_for_rho_prime_nn_twostage(cS.rho_prime_payg_fixed, K_global_guess, cS, paramS_for_solver, eIdxM, ...
    nn_models_loaded); % 传递包含所有NN模型和参数的结构体
solve_time = toc(solve_start_time);

% ... (后续的均衡结果打印和错误检查与你之前版本相同，此处省略) ...
fprintf('均衡求解完成。耗时: %.2f 秒。\n', solve_time);
fprintf('  均衡求解器返回状态: eq_found = %d\n', eq_found);
fprintf('  均衡结果: K_eq = %.4f, tau_l_eq = %.4f, GBC 残差 = %.3e\n', K_eq, tau_l_eq, gbc_residual_eq);
if ~eq_found || isnan(K_eq) || isnan(tau_l_eq)
    error('未能为固定的 rho_prime_payg_fixed = %.3f 找到均衡解。', cS.rho_prime_payg_fixed);
end
if abs(gbc_residual_eq) > cS.gbc_tol_for_internal_loop * 10
    warning('最终均衡的GBC残差 (%.2e) 较大。', gbc_residual_eq);
end
theta_payg_optimal_calc = NaN;
if isfield(final_eq_solution_details, 'theta_payg_required_before_cap')
    theta_payg_optimal_calc = final_eq_solution_details.theta_payg_required_before_cap;
end
fprintf('  (理论所需PAYG税率，未考虑上限前: %.4f)\n', theta_payg_optimal_calc);


%% 5. 分析和绘制最终均衡结果
fprintf('\n--- 5. 最终均衡结果与绘图 (rho_prime_payg_fixed=%.4f, tau_l_eq=%.4f, TR_gov=0) ---\n', cS.rho_prime_payg_fixed, tau_l_eq);
% ... (paramS_eq 的准备与之前版本相同) ...
if isempty(fields(final_eq_solution_details)) || isnan(K_eq) || ...
   ~isfield(final_eq_solution_details, 'MPL_gross') || isnan(final_eq_solution_details.MPL_gross)
    error('最终均衡的详细信息未能获取或无效，无法进行分析。');
end
paramS_eq = paramS; 
paramS_eq.tau_l = tau_l_eq;
if isfield(final_eq_solution_details, 'theta_payg') && isfinite(final_eq_solution_details.theta_payg)
    paramS_eq.theta_payg_actual_for_hh = final_eq_solution_details.theta_payg;
else
    % ... (重新计算 theta_payg_actual_for_hh 的逻辑) ...
end
paramS_eq.pps_tax_deferral_active = cS.pps_active;
[R_mkt_gross_factor_eq_final, MPL_gross_eq_final] = main_olg_v8_utils.HHPrices_Huggett(K_eq, paramS.L_per_capita, cS);
r_mkt_gross_eq_final = R_mkt_gross_factor_eq_final - 1;
r_k_net_hh_eq_final = r_mkt_gross_eq_final * (1 - cS.tau_k);
R_k_net_factor_hh_eq_final = 1 + r_k_net_hh_eq_final;
avg_worker_gross_wage_eq_final = (MPL_gross_eq_final * paramS.L_per_capita) / paramS.mass_workers_group;
b_payg_eq_final = cS.rho_prime_payg_fixed * avg_worker_gross_wage_eq_final;
bV_eq_new_final = zeros(1, cS.aD_new);
if cS.aR_new < cS.aD_new, bV_eq_new_final(cS.aR_new + 1 : cS.aD_new) = b_payg_eq_final; end
T_bequest_eq_final = 0;
if isfield(final_eq_solution_details, 'T_bequest_Model') && isfinite(final_eq_solution_details.T_bequest_Model)
    T_bequest_eq_final = final_eq_solution_details.T_bequest_Model;
end
TR_total_eq_final_for_vfi = T_bequest_eq_final;
fprintf('最终 家庭问题求解器 调用参数: MPL_gross=%.4f, tau_l=%.4f, theta_payg_actual=%.4f, TR_total=%.4f (T_bequest)\n', ...
    MPL_gross_eq_final, paramS_eq.tau_l, paramS_eq.theta_payg_actual_for_hh, TR_total_eq_final_for_vfi);


vfi_nn_mode_final_str_main = 'HHSolution_VFI_Huggett (V8)';

% 检查神经网络可用性（兼容两种加载方式）
nn_available_final = false;
three_networks_mode = false;

if cS.use_NN_for_VFI
    % 检查三网络方式的字段
    if isfield(nn_models_loaded, 'trainedNet_pps_classifier') && ~isempty(nn_models_loaded.trainedNet_pps_classifier) && ...
       isfield(nn_models_loaded, 'trainedNet_k_prime_regressor') && ~isempty(nn_models_loaded.trainedNet_k_prime_regressor)
        nn_available_final = true;
        three_networks_mode = true;
        vfi_nn_mode_final_str_main = 'HHSolution_NN_Huggett (V8)';
    % 检查分离文件方式的字段
    elseif isfield(nn_models_loaded, 'k_prime_net') && ~isempty(nn_models_loaded.k_prime_net) && ...
           isfield(nn_models_loaded, 'c_pps_classifier_net') && ~isempty(nn_models_loaded.c_pps_classifier_net)
        nn_available_final = true;
        three_networks_mode = false;
        vfi_nn_mode_final_str_main = 'HHSolution_NN_Huggett_two_stage (V8)';
    end
end

if nn_available_final
    fprintf('调用最终的 %s ...\n', vfi_nn_mode_final_str_main);
    
    if three_networks_mode
        % 使用新的三网络方式调用
        [kPolM_eq, cPpsPolM_choice_eq, cPolM_eq] = main_olg_v8_utils.HHSolution_NN_Huggett(...
            R_k_net_factor_hh_eq_final, MPL_gross_eq_final, TR_total_eq_final_for_vfi, bV_eq_new_final, ...
            paramS_eq, cS, ...
            nn_models_loaded.trainedNet_pps_classifier, ...
            nn_models_loaded.trainedNet_pps_regressor, ...
            nn_models_loaded.trainedNet_k_prime_regressor, ...
            nn_models_loaded.total_input_dim, ...
            nn_models_loaded.use_b_payg_as_input);
        valM_equilibrium = []; % 三网络方式不返回值函数
    else
        % 使用旧的分离文件方式调用
        [cPolM_eq, kPolM_eq, cPpsPolM_choice_eq, valM_equilibrium] = main_olg_v8_utils.HHSolution_NN_Huggett_two_stage(...
            R_k_net_factor_hh_eq_final, MPL_gross_eq_final, TR_total_eq_final_for_vfi, bV_eq_new_final, paramS_eq, cS, ...
            nn_models_loaded.k_prime_net, nn_models_loaded.k_prime_norm_params, ...
            nn_models_loaded.c_pps_classifier_net, nn_models_loaded.c_pps_classifier_norm_params, ...
            nn_models_loaded.c_pps_regressor_net, nn_models_loaded.c_pps_regressor_norm_params);
    end
else
    if cS.use_NN_for_VFI, fprintf('    (警告: cS.use_NN_for_VFI为真但必要NN未加载，将回退到完整VFI进行最终求解)\n'); end
    fprintf('调用最终的 %s ...\n', vfi_nn_mode_final_str_main);
    [cPolM_eq, kPolM_eq, cPpsPolM_choice_eq, valM_equilibrium] = main_olg_v8_utils.HHSolution_VFI_Huggett(...
        R_k_net_factor_hh_eq_final, MPL_gross_eq_final, TR_total_eq_final_for_vfi, bV_eq_new_final, paramS_eq, cS);
end

% ... (后续的模拟、汇总、GBC检查、绘图等与你之前版本相同，此处省略) ...
fprintf('模拟最终均衡的分布 (HHSimulation_olgm V8)...\n');
[kHistM_eq, kPpsHistM_eq, cHistM_eq] = main_olg_v8_utils.HHSimulation_olgm(...
    kPolM_eq, cPpsPolM_choice_eq, cPolM_eq, eIdxM, ... 
    R_k_net_factor_hh_eq_final, MPL_gross_eq_final, ...
    TR_total_eq_final_for_vfi, bV_eq_new_final, paramS_eq, cS);
K_nonpps_eq_agg = mean(kHistM_eq, 1) * paramS.Z_ss_norm_annual; K_pps_eq_agg = 0;
if cS.pps_active && cS.pps_in_K && (cS.pps_max_contrib_frac > 0 || cS.pps_annual_contrib_limit > 0) && ~isempty(kPpsHistM_eq)
    K_pps_eq_agg = mean(kPpsHistM_eq, 1) * paramS.Z_ss_norm_annual; end
Actual_K_eq_final = K_nonpps_eq_agg + K_pps_eq_agg;
C_eq_final = mean(cHistM_eq,1) * paramS.Z_ss_norm_annual;
Y_eq_final = cS.A * (Actual_K_eq_final^cS.alpha) * (paramS.L_per_capita^(1-cS.alpha));
G_eq_final = cS.gov_exp_frac_Y * Y_eq_final; B_eq_final = cS.gov_debt_frac_Y * Y_eq_final;
fprintf('\n--- V8 最终均衡汇总 ---\n');
% ... (所有汇总打印) ...
final_gbc_residual = main_olg_v8_utils.check_gbc_residual(Actual_K_eq_final, C_eq_final, Y_eq_final, ...
    G_eq_final, B_eq_final, MPL_gross_eq_final, r_mkt_gross_eq_final, ...
    paramS_eq.theta_payg_actual_for_hh, tau_l_eq, ...
    b_payg_eq_final, T_bequest_eq_final, 0, cS, paramS_eq);
fprintf('最终GBC(一般预算)检查 @ 均衡状态: Residual = %.4e\n', final_gbc_residual);
if abs(final_gbc_residual) > 1e-2 && ~(cS.use_NN_for_VFI && abs(final_gbc_residual) < 5e-2) 
    warning('最终GBC残差 (%.3e) 较大。', final_gbc_residual);
end
fprintf('\n绘制最终均衡的策略函数...\n');
% ... (绘图逻辑) ...
fprintf('\n--- V8 OLG 模型 (内生PPS缴费, 固定 Rho_prime_payg) 分析完成 ---\n');

% --- END OF FILE main_olg_v8_nn.m ---