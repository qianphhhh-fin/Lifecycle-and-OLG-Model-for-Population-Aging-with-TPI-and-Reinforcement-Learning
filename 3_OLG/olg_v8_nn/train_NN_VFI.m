% --- START OF MODIFIED FILE train_NN_VFI_three_networks.m ---

%% 0. 设置训练环境和参数
clc;
clear;
close all;
fprintf('=== 训练神经网络替代VFI (OLG V8) - PPS分类+回归+k_prime比例 ===\n');

% --- 用户选项：是否重新生成训练数据 ---
REGENERATE_DATA = false; % <<<< 设置为 true 则重新运行VFI生成数据
                        % <<<< 设置为 false 则尝试加载已保存的数据
DATA_FILENAME = 'generated_vfi_training_data_v8_proportion.mat'; % 保存/加载数据的文件名

% 加载基础模型参数
fprintf('--- 0.1. 加载基础模型参数 ---\n');
cS = main_olg_v8_utils.ParameterValues_HuggettStyle();
paramS_for_vfi = struct(); % 这个结构体主要用于VFI数据生成
fprintf('--- 0.2. 预计算劳动过程 (用于VFI数据生成) ---\n');
[paramS_for_vfi.leLogGridV, paramS_for_vfi.leTrProbM, paramS_for_vfi.leProb1V] = main_olg_v8_utils.EarningProcess_olgm(cS);
paramS_for_vfi.leGridV = exp(paramS_for_vfi.leLogGridV(:));
paramS_for_vfi.ageEffV_new = cS.ageEffV_new;

fprintf('--- 0.3. 设置神经网络训练参数 ---\n');
num_macro_samples = 100; 
points_per_vfi_sample = 5000; 
validation_split = 0.15;
test_split = 0.15;

% --- 宏观变量采样范围 ---
r_k_net_hh_min = 0.005; r_k_net_hh_max = 0.06;
R_k_net_factor_min = 1 + r_k_net_hh_min; R_k_net_factor_max = 1 + r_k_net_hh_max;
MPL_gross_min = 0.5 * cS.tgWage; MPL_gross_max = 1.5 * cS.tgWage;
tau_l_min_sample = cS.tau_l_min; tau_l_max_sample = cS.tau_l_max;
theta_payg_min = 0.01; theta_payg_max_train = cS.theta_payg_max * 0.95;
TR_total_min = 0; TR_total_max = 0.1 * MPL_gross_max;
b_payg_min_train = 0; b_payg_max_factor = 0.6;
use_b_payg_as_input = true;
cS.rho_prime_payg_fixed_for_training = 0.20;

% --- 确定输入输出维度 ---
num_k_inputs = 1; num_k_pps_inputs = 1;
num_eps_inputs = cS.nw; num_a_inputs = cS.aD_new;
num_macro_inputs = 5 + (use_b_payg_as_input * 1);
total_input_dim = num_k_inputs + num_k_pps_inputs + num_eps_inputs + num_a_inputs + num_macro_inputs;
fprintf('神经网络输入维度: %d\n', total_input_dim);
fprintf('将训练三个网络：\n');
fprintf('  1. PPS二分类器：是否参与PPS\n');
fprintf('  2. PPS比例回归器：缴费比例(0-1)\n');
fprintf('  3. k_prime比例回归器：储蓄比例(0-1)\n');

%% 1. 数据生成或加载
X_data_all = [];
Y_data_pps_binary = []; % PPS二分类标签
Y_data_pps_proportion = []; % PPS缴费比例 (仅对参与PPS的样本)
Y_data_k_prime_proportion = []; % k_prime比例

if REGENERATE_DATA || ~exist(DATA_FILENAME, 'file')
    if ~REGENERATE_DATA && ~exist(DATA_FILENAME, 'file')
        fprintf('未找到已保存的数据文件 "%s"，将重新生成数据。\n', DATA_FILENAME);
    end
    fprintf('\n--- 1. 开始生成训练数据 (可能非常耗时) ---\n');
    rng(123); % For reproducibility of sampling

    for i_sample = 1:num_macro_samples
        fprintf('  --- 生成宏观样本 %d / %d ---\n', i_sample, num_macro_samples);
        current_R_k_net_factor = R_k_net_factor_min + rand() * (R_k_net_factor_max - R_k_net_factor_min);
        current_MPL_gross = MPL_gross_min + rand() * (MPL_gross_max - MPL_gross_min);
        current_tau_l = tau_l_min_sample + rand() * (tau_l_max_sample - tau_l_min_sample);
        current_theta_payg_actual = theta_payg_min + rand() * (theta_payg_max_train - theta_payg_min);
        current_theta_payg_actual = min(current_theta_payg_actual, cS.max_total_labor_tax - current_tau_l);
        current_theta_payg_actual = max(0, current_theta_payg_actual);
        avg_worker_gross_wage_approx = current_MPL_gross; 
        current_b_payg = cS.rho_prime_payg_fixed_for_training * avg_worker_gross_wage_approx;
        current_b_payg = max(0, current_b_payg);
        if ~use_b_payg_as_input, current_b_payg = 0; end
        current_TR_total = TR_total_min + rand() * (TR_total_max - TR_total_min);
        
        fprintf('    当前宏观: R_k_net_fac=%.4f, MPL=%.3f, tau_l=%.3f, theta=%.3f, TR=%.3f, b_payg=%.3f\n', ...
                current_R_k_net_factor, current_MPL_gross, current_tau_l, current_theta_payg_actual, current_TR_total, current_b_payg);

        paramS_vfi_iter = paramS_for_vfi;
        paramS_vfi_iter.tau_l = current_tau_l;
        paramS_vfi_iter.theta_payg_actual_for_hh = current_theta_payg_actual;
        paramS_vfi_iter.pps_tax_deferral_active = cS.pps_active;
        bV_payg_vfi_iter = zeros(1, cS.aD_new);
        if cS.aR_new < cS.aD_new, bV_payg_vfi_iter(cS.aR_new + 1 : cS.aD_new) = current_b_payg; end

        vfi_timer = tic; fprintf('    调用 HHSolution_VFI_Huggett (V8) ... ');
        try
            [~, kPolM_vfi, cPpsPolM_choice_vfi, ~] = main_olg_v8_utils.HHSolution_VFI_Huggett(...
                current_R_k_net_factor, current_MPL_gross, current_TR_total, bV_payg_vfi_iter, paramS_vfi_iter, cS);
            fprintf('完成 (耗时 %.2f s)\n', toc(vfi_timer));
        catch ME_vfi
            fprintf('\n错误: VFI 求解失败 (样本 %d): %s\n', i_sample, ME_vfi.message); continue;
        end
        
        num_states_total = cS.nk * cS.nkpps * cS.nw * cS.aD_new;
        
        % 首先找到所有 c_pps > 0 的样本索引
        idx_c_pps_positive_flat = [];
        for idx_flat = 1:num_states_total
            [ik, ikpps, ie, ia] = ind2sub([cS.nk, cS.nkpps, cS.nw, cS.aD_new], idx_flat);
            if cPpsPolM_choice_vfi(ik, ikpps, ie, ia) > 1e-7
                idx_c_pps_positive_flat = [idx_c_pps_positive_flat, idx_flat];
            end
        end
        
        % 保留所有 c_pps > 0 的样本
        num_c_pps_positive_current = length(idx_c_pps_positive_flat);
        fprintf('    当前宏观样本中 c_pps > 0 的样本数: %d\n', num_c_pps_positive_current);
        
        % 从剩余样本中随机采样，补充到目标数量
        remaining_indices = setdiff(1:num_states_total, idx_c_pps_positive_flat);
        num_additional_needed = max(0, points_per_vfi_sample - num_c_pps_positive_current);
        num_additional_available = length(remaining_indices);
        num_additional_to_sample = min(num_additional_needed, num_additional_available);
        
        if num_additional_to_sample > 0
            additional_sampled_indices = remaining_indices(randperm(num_additional_available, num_additional_to_sample));
            sampled_indices = [idx_c_pps_positive_flat, additional_sampled_indices];
        else
            sampled_indices = idx_c_pps_positive_flat;
        end
        
        actual_points_to_sample = length(sampled_indices);
        fprintf('    保留 %d 个 c_pps > 0 样本，额外采样 %d 个样本，总计 %d 个样本\n', ...
                num_c_pps_positive_current, num_additional_to_sample, actual_points_to_sample);
        
        temp_X = zeros(actual_points_to_sample, total_input_dim);
        temp_Y_pps_binary = zeros(actual_points_to_sample, 1);
        temp_Y_pps_proportion = zeros(actual_points_to_sample, 1);
        temp_Y_k_prime_proportion = zeros(actual_points_to_sample, 1);
        
        count = 0;
        for idx_flat = sampled_indices
            count = count + 1;
            [ik, ikpps, ie, ia] = ind2sub([cS.nk, cS.nkpps, cS.nw, cS.aD_new], idx_flat);
            current_k = cS.kGridV(ik); current_k_pps = cS.kppsGridV(ikpps);
            eps_one_hot = zeros(1, cS.nw); eps_one_hot(ie) = 1;
            a_one_hot = zeros(1, cS.aD_new); a_one_hot(ia) = 1;
            macro_vars_vec = [current_R_k_net_factor, current_MPL_gross, current_tau_l, current_theta_payg_actual, current_TR_total];
            if use_b_payg_as_input, macro_vars_vec = [macro_vars_vec, current_b_payg]; end
            temp_X(count, :) = [current_k, current_k_pps, eps_one_hot, a_one_hot, macro_vars_vec];
            
            % 计算目标变量
            epsilon_val = paramS_vfi_iter.leGridV(ie);
            pps_binary = 0; pps_proportion = 0; k_prime_proportion = 0;
            
            % 计算PPS二分类和比例
            actual_c_pps = cPpsPolM_choice_vfi(ik, ikpps, ie, ia);
            if actual_c_pps > 1e-7
                pps_binary = 1; % 参与PPS
                
                % 计算PPS比例
                if ia <= cS.aR_new && cS.physAgeMap{ia}(1) <= cS.pps_contribution_age_max_idx && cS.pps_active
                    age_efficiency = cS.ageEffV_new(ia);
                    current_gross_labor_income = current_MPL_gross * age_efficiency * epsilon_val;
                    if current_gross_labor_income > 1e-6
                        max_cpps_by_frac = current_gross_labor_income * cS.pps_max_contrib_frac;
                        max_permissible_cpps = min(cS.pps_annual_contrib_limit, max_cpps_by_frac);
                        max_permissible_cpps = max(0, max_permissible_cpps);
                        if max_permissible_cpps > 1e-6
                            pps_proportion = min(1.0, actual_c_pps / max_permissible_cpps);
                        end
                    end
                end
            else
                pps_binary = 0; % 不参与PPS
                pps_proportion = 0; % 这个值对于不参与PPS的样本不重要
            end
            
            % 计算k_prime比例
            b_age_val = 0;
            if ia > cS.aR_new, b_age_val = current_b_payg; end
            
            [resources_for_c_k_prime, ~, ~] = main_olg_v8_utils.HHIncome_Huggett(...
                current_k, current_R_k_net_factor, current_MPL_gross, current_TR_total, b_age_val, ...
                actual_c_pps, ia, paramS_vfi_iter, cS, epsilon_val);
            
            consumption_floor_spending = cS.cFloor * (1 + cS.tau_c);
            resources_for_kprime_and_c_above_floor = resources_for_c_k_prime - consumption_floor_spending;
            
            if resources_for_kprime_and_c_above_floor > 1e-6
                k_prime_proportion = min(1.0, (kPolM_vfi(ik, ikpps, ie, ia) - cS.kMin) / resources_for_kprime_and_c_above_floor);
            end
            
            temp_Y_pps_binary(count, 1) = pps_binary;
            temp_Y_pps_proportion(count, 1) = pps_proportion;
            temp_Y_k_prime_proportion(count, 1) = k_prime_proportion;
        end
        
        X_data_all = [X_data_all; temp_X];
        Y_data_pps_binary = [Y_data_pps_binary; temp_Y_pps_binary];
        Y_data_pps_proportion = [Y_data_pps_proportion; temp_Y_pps_proportion];
        Y_data_k_prime_proportion = [Y_data_k_prime_proportion; temp_Y_k_prime_proportion];
        
        fprintf('    为样本 %d 添加了 %d 个数据点。总数据点数: %d\n', i_sample, actual_points_to_sample, size(X_data_all, 1));
    end

    if isempty(X_data_all)
        error('未能生成任何训练数据。请检查VFI过程或采样点数。');
    end
    fprintf('--- 数据生成完成。总共 %d 个数据点。---\n', size(X_data_all, 1));
    
    % 保存生成的数据
    fprintf('正在保存生成的训练数据到: %s\n', DATA_FILENAME);
    save(DATA_FILENAME, 'X_data_all', 'Y_data_pps_binary', 'Y_data_pps_proportion', 'Y_data_k_prime_proportion', 'cS', 'paramS_for_vfi', '-v7.3');
    fprintf('数据已保存。\n');
else
    fprintf('\n--- 1. 正在从 "%s" 加载已保存的训练数据 ---\n', DATA_FILENAME);
    
    % 首先检查文件是否存在
    if ~exist(DATA_FILENAME, 'file')
        error('数据文件 "%s" 不存在。请检查文件路径或设置 REGENERATE_DATA = true 重新生成数据。', DATA_FILENAME);
    end
    
    % 查看文件中的变量
    fprintf('检查数据文件中的变量...\n');
    file_info = whos('-file', DATA_FILENAME);
    fprintf('文件中包含的变量:\n');
    for i = 1:length(file_info)
        fprintf('  %s: [%s] %s\n', file_info(i).name, num2str(file_info(i).size), file_info(i).class);
    end
    
    % 尝试加载新格式的数据
    new_format_success = false;
    try
        load(DATA_FILENAME, 'X_data_all', 'Y_data_pps_binary', 'Y_data_pps_proportion', 'Y_data_k_prime_proportion', 'cS');
        if exist('X_data_all', 'var') && exist('Y_data_pps_binary', 'var') && exist('Y_data_pps_proportion', 'var') && exist('Y_data_k_prime_proportion', 'var')
            % 检查新格式数据的完整性
            if ~isempty(X_data_all) && ~isempty(Y_data_pps_binary) && ~isempty(Y_data_pps_proportion) && ~isempty(Y_data_k_prime_proportion)
                if size(Y_data_pps_binary, 1) == size(X_data_all, 1) && ...
                   size(Y_data_pps_proportion, 1) == size(X_data_all, 1) && ...
                   size(Y_data_k_prime_proportion, 1) == size(X_data_all, 1)
                    fprintf('成功加载三网络格式数据。\n');
                    fprintf('  X_data_all: [%d, %d]\n', size(X_data_all, 1), size(X_data_all, 2));
                    fprintf('  Y_data_pps_binary: [%d, %d]\n', size(Y_data_pps_binary, 1), size(Y_data_pps_binary, 2));
                    fprintf('  Y_data_pps_proportion: [%d, %d]\n', size(Y_data_pps_proportion, 1), size(Y_data_pps_proportion, 2));
                    fprintf('  Y_data_k_prime_proportion: [%d, %d]\n', size(Y_data_k_prime_proportion, 1), size(Y_data_k_prime_proportion, 2));
                    new_format_success = true;
                else
                    fprintf('新格式数据尺寸不匹配，将尝试旧格式。\n');
                end
            else
                fprintf('新格式数据为空，将尝试旧格式。\n');
            end
        else
            fprintf('新格式数据变量缺失，将尝试旧格式。\n');
        end
    catch ME_new
        fprintf('新格式加载失败: %s\n', ME_new.message);
    end
    
    if ~new_format_success
        % 尝试加载旧格式的比例数据
        fprintf('尝试加载旧版本的比例数据格式...\n');
        try
            % 先清除可能存在的变量
            clear Y_data_pps_binary Y_data_pps_proportion Y_data_k_prime_proportion Y_data_proportions
            
            load(DATA_FILENAME, 'X_data_all', 'Y_data_proportions', 'cS');
            
            % 检查加载的变量
            fprintf('旧格式数据加载检查:\n');
            if exist('X_data_all', 'var')
                fprintf('  X_data_all: [%d, %d]\n', size(X_data_all, 1), size(X_data_all, 2));
            else
                error('X_data_all 变量不存在');
            end
            
            if exist('Y_data_proportions', 'var')
                fprintf('  Y_data_proportions: [%d, %d]\n', size(Y_data_proportions, 1), size(Y_data_proportions, 2));
                if size(Y_data_proportions, 2) ~= 2
                    error('Y_data_proportions 应该有2列，实际有 %d 列', size(Y_data_proportions, 2));
                end
                if size(Y_data_proportions, 1) ~= size(X_data_all, 1)
                    error('Y_data_proportions 的行数 (%d) 与 X_data_all 的行数 (%d) 不匹配', size(Y_data_proportions, 1), size(X_data_all, 1));
                end
            else
                error('Y_data_proportions 变量不存在');
            end
            
            fprintf('检测到旧版本数据格式，正在转换...\n');
            
            % 从旧的比例数据生成三个新变量
            num_samples = size(Y_data_proportions, 1);
            
            % 创建新的Y变量
            Y_data_pps_proportion = Y_data_proportions(:, 1); % 第一列是PPS比例
            Y_data_k_prime_proportion = Y_data_proportions(:, 2); % 第二列是k_prime比例
            Y_data_pps_binary = zeros(num_samples, 1);
            
            % 生成PPS二分类标签：如果PPS比例>0则为1，否则为0
            Y_data_pps_binary(Y_data_pps_proportion > 1e-7) = 1;
            
            % 对于不参与PPS的样本，确保其PPS比例为0
            Y_data_pps_proportion(Y_data_pps_binary == 0) = 0;
            
            % 验证转换结果
            fprintf('转换后的数据验证:\n');
            fprintf('  Y_data_pps_binary: [%d, %d]\n', size(Y_data_pps_binary, 1), size(Y_data_pps_binary, 2));
            fprintf('  Y_data_pps_proportion: [%d, %d]\n', size(Y_data_pps_proportion, 1), size(Y_data_pps_proportion, 2));
            fprintf('  Y_data_k_prime_proportion: [%d, %d]\n', size(Y_data_k_prime_proportion, 1), size(Y_data_k_prime_proportion, 2));
            
            pps_positive_count = sum(Y_data_pps_binary);
            fprintf('从旧数据转换完成：\n');
            fprintf('  总样本数: %d\n', num_samples);
            fprintf('  PPS参与样本: %d (%.1f%%)\n', pps_positive_count, pps_positive_count/num_samples*100);
            fprintf('  PPS比例范围: [%.3f, %.3f]\n', min(Y_data_pps_proportion), max(Y_data_pps_proportion));
            fprintf('  k_prime比例范围: [%.3f, %.3f]\n', min(Y_data_k_prime_proportion), max(Y_data_k_prime_proportion));
            
        catch ME_old
            error('无法加载数据文件 "%s"。\n新格式加载失败: %s\n旧格式加载失败: %s\n请检查文件是否存在或重新生成数据。', DATA_FILENAME, ME_new.message, ME_old.message);
        end
    end
    
    fprintf('数据加载完成。总共 %d 个数据点。\n', size(X_data_all, 1));
end

%% 2. 数据预处理和统计
fprintf('\n--- 2. 数据预处理 (无归一化) ---\n');

% 检查数据一致性
fprintf('数据一致性检查:\n');
fprintf('  X_data_all 大小: [%d, %d]\n', size(X_data_all, 1), size(X_data_all, 2));
fprintf('  Y_data_pps_binary 大小: [%d, %d]\n', size(Y_data_pps_binary, 1), size(Y_data_pps_binary, 2));
fprintf('  Y_data_pps_proportion 大小: [%d, %d]\n', size(Y_data_pps_proportion, 1), size(Y_data_pps_proportion, 2));
fprintf('  Y_data_k_prime_proportion 大小: [%d, %d]\n', size(Y_data_k_prime_proportion, 1), size(Y_data_k_prime_proportion, 2));

% 确保所有Y数据的行数与X_data_all一致
num_samples_X = size(X_data_all, 1);
if size(Y_data_pps_binary, 1) ~= num_samples_X
    error('Y_data_pps_binary 的行数 (%d) 与 X_data_all 的行数 (%d) 不匹配', size(Y_data_pps_binary, 1), num_samples_X);
end
if size(Y_data_pps_proportion, 1) ~= num_samples_X
    error('Y_data_pps_proportion 的行数 (%d) 与 X_data_all 的行数 (%d) 不匹配', size(Y_data_pps_proportion, 1), num_samples_X);
end
if size(Y_data_k_prime_proportion, 1) ~= num_samples_X
    error('Y_data_k_prime_proportion 的行数 (%d) 与 X_data_all 的行数 (%d) 不匹配', size(Y_data_k_prime_proportion, 1), num_samples_X);
end

fprintf('数据一致性检查通过。\n');

% 检查数据范围
fprintf('输入数据统计:\n');
fprintf('  k 范围: [%.3f, %.3f]\n', min(X_data_all(:,1)), max(X_data_all(:,1)));
fprintf('  k_pps 范围: [%.3f, %.3f]\n', min(X_data_all(:,2)), max(X_data_all(:,2)));
macro_start_col = 2 + cS.nw + cS.aD_new + 1;
fprintf('  宏观变量范围:\n');
for i = 1:num_macro_inputs
    col_idx = macro_start_col + i - 1;
    fprintf('    宏观变量%d: [%.4f, %.4f]\n', i, min(X_data_all(:,col_idx)), max(X_data_all(:,col_idx)));
end

% 统计PPS参与情况
pps_positive_count = sum(Y_data_pps_binary);
pps_positive_rate = pps_positive_count / length(Y_data_pps_binary);
fprintf('输出数据统计:\n');
fprintf('  PPS参与: %d / %d (%.1f%%)\n', pps_positive_count, length(Y_data_pps_binary), pps_positive_rate * 100);

% 对于参与PPS的样本，统计比例分布
pps_positive_indices = find(Y_data_pps_binary > 0.5);
if ~isempty(pps_positive_indices)
    pps_proportions_positive = Y_data_pps_proportion(pps_positive_indices);
    fprintf('  PPS比例 (仅参与者): 范围 [%.3f, %.3f], 均值 %.3f\n', ...
        min(pps_proportions_positive), max(pps_proportions_positive), mean(pps_proportions_positive));
else
    fprintf('  PPS比例: 无参与者\n');
end

fprintf('  k_prime比例: 范围 [%.3f, %.3f], 均值 %.3f\n', ...
    min(Y_data_k_prime_proportion), max(Y_data_k_prime_proportion), mean(Y_data_k_prime_proportion));

% 划分训练/验证/测试集
num_total_points = size(X_data_all, 1);
fprintf('\n数据集划分准备:\n');
fprintf('  总样本数: %d\n', num_total_points);
fprintf('  验证集比例: %.1f%%, 测试集比例: %.1f%%\n', validation_split*100, test_split*100);

cv_all = cvpartition(num_total_points, 'HoldOut', validation_split + test_split);
idx_train = training(cv_all);
idx_val_test = test(cv_all);
num_val_test = sum(idx_val_test);

fprintf('初始划分: 训练集 %d, 验证+测试集 %d\n', sum(idx_train), num_val_test);

if num_val_test > 0
    cv_val = cvpartition(num_val_test, 'HoldOut', test_split / (validation_split + test_split));
    temp_indices = find(idx_val_test);
    idx_val = temp_indices(training(cv_val));
    idx_test = temp_indices(test(cv_val));
    fprintf('二次划分: 验证集 %d, 测试集 %d\n', length(idx_val), length(idx_test));
else
    idx_val = []; idx_test = [];
    warning('数据集太小，无法进行验证/测试集划分。');
end

% 准备训练集
fprintf('\n准备训练集数据...\n');
X_train = X_data_all(idx_train, :);
Y_train_pps_binary = Y_data_pps_binary(idx_train, :);
Y_train_pps_proportion = Y_data_pps_proportion(idx_train, :);
Y_train_k_prime_proportion = Y_data_k_prime_proportion(idx_train, :);

% 准备验证集
if ~isempty(idx_val)
    fprintf('准备验证集数据...\n');
    X_val = X_data_all(idx_val, :);
    Y_val_pps_binary = Y_data_pps_binary(idx_val, :);
    Y_val_pps_proportion = Y_data_pps_proportion(idx_val, :);
    Y_val_k_prime_proportion = Y_data_k_prime_proportion(idx_val, :);
else
    X_val = []; Y_val_pps_binary = []; Y_val_pps_proportion = []; Y_val_k_prime_proportion = [];
end

% 准备测试集
if ~isempty(idx_test)
    fprintf('准备测试集数据...\n');
    X_test = X_data_all(idx_test, :);
    Y_test_pps_binary = Y_data_pps_binary(idx_test, :);
    Y_test_pps_proportion = Y_data_pps_proportion(idx_test, :);
    Y_test_k_prime_proportion = Y_data_k_prime_proportion(idx_test, :);
else
    X_test = []; Y_test_pps_binary = []; Y_test_pps_proportion = []; Y_test_k_prime_proportion = [];
end

fprintf('数据集划分: 训练 %d, 验证 %d, 测试 %d\n', size(X_train,1), size(X_val,1), size(X_test,1));

%% 3. 训练三个神经网络
fprintf('\n--- 3. 训练三个神经网络 ---\n');

%% 3.1 训练PPS二分类器
fprintf('\n--- 3.1 训练PPS二分类器 ---\n');
Y_train_pps_binary_cat = categorical(Y_train_pps_binary, [0 1], {'false', 'true'});
if ~isempty(X_val)
    Y_val_pps_binary_cat = categorical(Y_val_pps_binary, [0 1], {'false', 'true'});
else
    Y_val_pps_binary_cat = [];
end

% 计算类别权重
class_counts = countcats(Y_train_pps_binary_cat);
fprintf('PPS分类器训练集类别分布: false=%d, true=%d\n', class_counts(1), class_counts(2));
if class_counts(2) > 0 && class_counts(1) > 0
    class_weights = [1, class_counts(1)/class_counts(2)]; % 给正样本更高权重
    fprintf('使用类别权重: false=%.2f, true=%.2f\n', class_weights(1), class_weights(2));
else
    class_weights = [1, 1];
    fprintf('警告: 类别不平衡严重，使用均匀权重\n');
end

layers_pps_classifier = [
    featureInputLayer(total_input_dim, 'Name', 'input_clf', 'Normalization', 'none')
    fullyConnectedLayer(128, 'Name', 'fc1_clf')
    reluLayer('Name', 'relu1_clf')
    dropoutLayer(0.2, 'Name', 'dropout1_clf')
    fullyConnectedLayer(64, 'Name', 'fc2_clf')
    reluLayer('Name', 'relu2_clf')
    dropoutLayer(0.2, 'Name', 'dropout2_clf')
    fullyConnectedLayer(2, 'Name', 'fc_output_clf') 
    softmaxLayer('Name', 'softmax_clf')
    classificationLayer('Name', 'output_clf', 'Classes', {'false', 'true'}, 'ClassWeights', class_weights)];

val_data_clf = {};
if ~isempty(X_val) && ~isempty(Y_val_pps_binary_cat)
    val_data_clf = {X_val, Y_val_pps_binary_cat};
end

options_pps_classifier = trainingOptions('adam', 'MaxEpochs', 100, 'MiniBatchSize', 128, ...
    'InitialLearnRate', 1e-3, 'LearnRateSchedule', 'piecewise', 'LearnRateDropFactor', 0.5, ...
    'LearnRateDropPeriod', 30, 'Shuffle', 'every-epoch', 'ValidationData', val_data_clf, ...
    'ValidationFrequency', 10, 'Plots', 'training-progress', 'Verbose', false, 'L2Regularization', 0.005);

[trainedNet_pps_classifier, trainInfo_pps_classifier] = trainNetwork(X_train, Y_train_pps_binary_cat, layers_pps_classifier, options_pps_classifier);
fprintf('PPS二分类器训练完成。\n');

%% 3.2 训练PPS比例回归器 (仅对参与PPS的样本)
fprintf('\n--- 3.2 训练PPS比例回归器 ---\n');
idx_train_pps_positive = find(Y_train_pps_binary > 0.5);
if length(idx_train_pps_positive) > 50 % 确保有足够的正样本
    X_train_pps_reg = X_train(idx_train_pps_positive, :);
    Y_train_pps_reg = Y_train_pps_proportion(idx_train_pps_positive, :);
    
    % 验证集中的正样本
    if ~isempty(X_val)
        idx_val_pps_positive = find(Y_val_pps_binary > 0.5);
        if ~isempty(idx_val_pps_positive)
            X_val_pps_reg = X_val(idx_val_pps_positive, :);
            Y_val_pps_reg = Y_val_pps_proportion(idx_val_pps_positive, :);
        else
            X_val_pps_reg = []; Y_val_pps_reg = [];
        end
    else
        X_val_pps_reg = []; Y_val_pps_reg = [];
    end
    
    fprintf('PPS回归器数据: 训练集 %d 样本, 验证集 %d 样本\n', size(X_train_pps_reg,1), size(X_val_pps_reg,1));
    
    layers_pps_regressor = [
        featureInputLayer(total_input_dim, 'Name', 'input_reg', 'Normalization', 'none')
        fullyConnectedLayer(128, 'Name', 'fc1_reg')
        reluLayer('Name', 'relu1_reg')
        fullyConnectedLayer(64, 'Name', 'fc2_reg')
        reluLayer('Name', 'relu2_reg')
        fullyConnectedLayer(1, 'Name', 'fc_output_reg')
        regressionLayer('Name', 'output_reg')];
    
    val_data_reg = {};
    if ~isempty(X_val_pps_reg)
        val_data_reg = {X_val_pps_reg, Y_val_pps_reg};
    end
    
    options_pps_regressor = trainingOptions('adam', 'MaxEpochs', 150, 'MiniBatchSize', 64, ...
        'InitialLearnRate', 1e-3, 'LearnRateSchedule', 'piecewise', 'LearnRateDropFactor', 0.5, ...
        'LearnRateDropPeriod', 40, 'Shuffle', 'every-epoch', 'ValidationData', val_data_reg, ...
        'ValidationFrequency', 10, 'Plots', 'training-progress', 'Verbose', false, 'L2Regularization', 0.001);
    
    [trainedNet_pps_regressor, trainInfo_pps_regressor] = trainNetwork(X_train_pps_reg, Y_train_pps_reg, layers_pps_regressor, options_pps_regressor);
    fprintf('PPS比例回归器训练完成。\n');
else
    fprintf('警告: PPS正样本太少 (%d)，跳过PPS比例回归器训练。\n', length(idx_train_pps_positive));
    trainedNet_pps_regressor = [];
    trainInfo_pps_regressor = [];
end

%% 3.3 训练k_prime比例回归器
fprintf('\n--- 3.3 训练k_prime比例回归器 ---\n');

layers_k_prime_regressor = [
    featureInputLayer(total_input_dim, 'Name', 'input_k', 'Normalization', 'none')
    fullyConnectedLayer(128, 'Name', 'fc1_k')
    reluLayer('Name', 'relu1_k')
    fullyConnectedLayer(64, 'Name', 'fc2_k')
    reluLayer('Name', 'relu2_k')
    fullyConnectedLayer(1, 'Name', 'fc_output_k')
    regressionLayer('Name', 'output_k')];

val_data_k = {};
if ~isempty(X_val)
    val_data_k = {X_val, Y_val_k_prime_proportion};
end

options_k_prime_regressor = trainingOptions('adam', 'MaxEpochs', 150, 'MiniBatchSize', 128, ...
    'InitialLearnRate', 1e-3, 'LearnRateSchedule', 'piecewise', 'LearnRateDropFactor', 0.5, ...
    'LearnRateDropPeriod', 40, 'Shuffle', 'every-epoch', 'ValidationData', val_data_k, ...
    'ValidationFrequency', 10, 'Plots', 'training-progress', 'Verbose', false, 'L2Regularization', 0.001);

[trainedNet_k_prime_regressor, trainInfo_k_prime_regressor] = trainNetwork(X_train, Y_train_k_prime_proportion, layers_k_prime_regressor, options_k_prime_regressor);
fprintf('k_prime比例回归器训练完成。\n');

%% 4. 保存网络和参数
fprintf('\n--- 4. 保存训练好的网络和参数 ---\n');
save('trained_VFI_NN_v8_three_networks.mat', ...
     'trainedNet_pps_classifier', 'trainInfo_pps_classifier', ...
     'trainedNet_pps_regressor', 'trainInfo_pps_regressor', ...
     'trainedNet_k_prime_regressor', 'trainInfo_k_prime_regressor', ...
     'cS', 'total_input_dim', 'use_b_payg_as_input');
fprintf('三个神经网络已保存到: trained_VFI_NN_v8_three_networks.mat\n');

%% 5. 测试网络预测
fprintf('\n--- 5. 测试集评估 ---\n');
if ~isempty(X_test)
    % 5.1 PPS分类器测试
    Y_test_pps_binary_cat = categorical(Y_test_pps_binary, [0 1], {'false', 'true'});
    Y_pred_pps_binary_cat = classify(trainedNet_pps_classifier, X_test);
    accuracy_pps_clf = sum(Y_pred_pps_binary_cat == Y_test_pps_binary_cat) / length(Y_test_pps_binary_cat);
    fprintf('PPS分类器准确率: %.3f\n', accuracy_pps_clf);
    
    % 5.2 PPS回归器测试 (仅对实际参与PPS的样本)
    if ~isempty(trainedNet_pps_regressor)
        idx_test_pps_actual_positive = find(Y_test_pps_binary > 0.5);
        if ~isempty(idx_test_pps_actual_positive)
            X_test_pps_actual = X_test(idx_test_pps_actual_positive, :);
            Y_test_pps_actual = Y_test_pps_proportion(idx_test_pps_actual_positive, :);
            Y_pred_pps_actual = predict(trainedNet_pps_regressor, X_test_pps_actual);
            Y_pred_pps_actual = max(0, min(1, Y_pred_pps_actual)); % 确保在[0,1]范围内
            rmse_pps_reg = sqrt(mean((Y_test_pps_actual - Y_pred_pps_actual).^2));
            fprintf('PPS比例回归器 RMSE (仅实际参与者): %.4f\n', rmse_pps_reg);
        else
            fprintf('测试集中无实际PPS参与者\n');
        end
    else
        fprintf('PPS比例回归器未训练\n');
    end
    
    % 5.3 k_prime回归器测试
    Y_pred_k_prime = predict(trainedNet_k_prime_regressor, X_test);
    Y_pred_k_prime = max(0, min(1, Y_pred_k_prime));
    rmse_k_prime = sqrt(mean((Y_test_k_prime_proportion - Y_pred_k_prime).^2));
    fprintf('k_prime比例回归器 RMSE: %.4f\n', rmse_k_prime);
    
    % 5.4 组合模型评估：模拟实际使用
    fprintf('\n--- 组合模型评估 ---\n');
    % 预测PPS参与概率
    Y_pred_pps_binary_prob = predict(trainedNet_pps_classifier, X_test);
    Y_pred_pps_binary_logical = (Y_pred_pps_binary_prob(:,2) > 0.5); % 第二列是'true'的概率
    
    % 对分类器预测为正的样本，用回归器预测比例
    Y_pred_pps_proportion_combined = zeros(size(Y_test_pps_proportion));
    if ~isempty(trainedNet_pps_regressor) && sum(Y_pred_pps_binary_logical) > 0
        idx_pred_positive = find(Y_pred_pps_binary_logical);
        X_pred_positive = X_test(idx_pred_positive, :);
        Y_pred_pps_subset = predict(trainedNet_pps_regressor, X_pred_positive);
        Y_pred_pps_subset = max(0, min(1, Y_pred_pps_subset));
        Y_pred_pps_proportion_combined(idx_pred_positive) = Y_pred_pps_subset;
    end
    
    % 评估组合模型的PPS预测
    combined_pps_rmse = sqrt(mean((Y_test_pps_proportion - Y_pred_pps_proportion_combined).^2));
    fprintf('组合模型PPS比例 RMSE: %.4f\n', combined_pps_rmse);
    
    % 可视化结果
    figure('Name', '三网络模型测试结果');
    
    subplot(2,2,1);
    confusionchart(Y_test_pps_binary_cat, Y_pred_pps_binary_cat);
    title('PPS分类器混淆矩阵');
    
    subplot(2,2,2);
    if ~isempty(trainedNet_pps_regressor) && ~isempty(idx_test_pps_actual_positive)
        scatter(Y_test_pps_actual(1:min(200,end)), Y_pred_pps_actual(1:min(200,end)), '.');
        hold on; plot([0,1], [0,1], 'r--'); hold off;
        xlabel('真实PPS比例'); ylabel('预测PPS比例');
        title('PPS比例回归 (实际参与者)'); grid on;
    else
        text(0.5, 0.5, '无PPS回归数据', 'HorizontalAlignment', 'center');
        title('PPS比例回归');
    end
    
    subplot(2,2,3);
    scatter(Y_test_k_prime_proportion(1:min(200,end)), Y_pred_k_prime(1:min(200,end)), '.');
    hold on; plot([0,1], [0,1], 'r--'); hold off;
    xlabel('真实k''比例'); ylabel('预测k''比例');
    title('k''比例回归'); grid on;
    
    subplot(2,2,4);
    scatter(Y_test_pps_proportion(1:min(200,end)), Y_pred_pps_proportion_combined(1:min(200,end)), '.');
    hold on; plot([0,1], [0,1], 'r--'); hold off;
    xlabel('真实PPS比例'); ylabel('组合模型预测');
    title('组合模型PPS预测'); grid on;
    
else
    fprintf('--- 5. 无测试集数据进行评估 ---\n');
end

fprintf('\n=== 三网络训练脚本执行完毕 ===\n');
% --- END OF MODIFIED FILE train_NN_VFI_three_networks.m ---