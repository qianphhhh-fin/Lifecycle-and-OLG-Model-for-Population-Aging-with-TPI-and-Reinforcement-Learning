% --- START OF FILE compare_rl_and_vfi.m ---
% 比较VFI和RL（SAC Agent）的优化结果 (期望化死亡率版本)
% 目标：在相同的宏观和微观参数下，比较两种方法的优化效果
% 特性：使用期望化死亡率，所有个体都完整生活到最大年龄

clc;
clear;
close all;
fprintf('=== 比较 VFI 和 RL（SAC Agent）的优化结果 ===\n');
fprintf('在相同的宏观和微观参数下比较两种方法的表现\n\n');

%% 1. 初始化共同参数
fprintf('--- 1. 初始化共同参数 ---\n');
cS = main_olg_v8_utils.ParameterValues_HuggettStyle();

% 确保网格已生成 (特别是 kppsGridV)
if ~isfield(cS, 'kppsGridV') || isempty(cS.kppsGridV)
    cS = main_olg_v8_utils.UpdateGrids(cS);
    fprintf('网格已生成\n');
end

fprintf('参数已加载：nk=%d, nkpps=%d, nw=%d\n', cS.nk, cS.nkpps, cS.nw);

% 设置固定的宏观经济环境（用于公平比较）
M_test.R_k_net_factor = 1.03;  % 净资本回报率
M_test.w_gross = 2.0;          % 总工资率
M_test.TR_total = 0.1;         % 总转移支付
M_test.b_payg_avg_retiree = 0.4; % 平均PAYG福利
M_test.tau_l = 0.15;           % 劳动所得税率
M_test.theta_payg_actual = 0.12; % 实际PAYG税率
M_test.b_payg_avg_for_obs = M_test.b_payg_avg_retiree;

fprintf('设定固定宏观参数用于比较：\n');
fprintf('  R_k_net_factor = %.3f\n', M_test.R_k_net_factor);
fprintf('  w_gross = %.3f\n', M_test.w_gross);
fprintf('  tau_l = %.3f\n', M_test.tau_l);
fprintf('  theta_payg_actual = %.3f\n', M_test.theta_payg_actual);

%% 2. 运行VFI方法
fprintf('\n--- 2. 运行VFI方法 ---\n');
vfi_start_time = tic;

% 初始化VFI需要的参数
paramS_vfi = struct();
[paramS_vfi.leLogGridV, paramS_vfi.leTrProbM, paramS_vfi.leProb1V] = main_olg_v8_utils.EarningProcess_olgm(cS);
paramS_vfi.leGridV = exp(paramS_vfi.leLogGridV(:));
paramS_vfi.ageEffV_new = cS.ageEffV_new;
paramS_vfi.tau_l = M_test.tau_l;
paramS_vfi.theta_payg_actual_for_hh = M_test.theta_payg_actual;
paramS_vfi.pps_tax_deferral_active = cS.pps_active;

% 构建PAYG福利向量
bV_payg_vfi = zeros(1, cS.aD_new);
if cS.aR_new < cS.aD_new
    bV_payg_vfi(cS.aR_new + 1 : cS.aD_new) = M_test.b_payg_avg_retiree;
end

% 运行VFI
fprintf('开始VFI求解...\n');
[cPolM_vfi, kPolM_vfi, cPpsPolM_vfi, VPolM_vfi] = ...
    main_olg_v8_utils.HHSolution_VFI_Huggett(...
    M_test.R_k_net_factor, M_test.w_gross, M_test.TR_total, bV_payg_vfi, paramS_vfi, cS);

vfi_time = toc(vfi_start_time);
fprintf('VFI求解完成。耗时: %.2f 秒\n', vfi_time);
fprintf('VFI策略矩阵尺寸: cPol=%s, kPol=%s, cPpsPol=%s\n', ...
    mat2str(size(cPolM_vfi)), mat2str(size(kPolM_vfi)), mat2str(size(cPpsPolM_vfi)));

%% 3. 加载RL（SAC）Agent并转换策略
fprintf('\n--- 3. 加载RL（SAC）Agent并转换策略 ---\n');

% 检查是否存在训练好的SAC Agent文件
sac_agent_file = 'final_SAC_Agent_OLG_R2024b.mat';
if ~exist(sac_agent_file, 'file')
    warning('未找到SAC Agent文件: %s', sac_agent_file);
    fprintf('请先运行 main_olg_v8_SAC.m 来训练SAC Agent。\n');
    return;
end

% 加载训练好的SAC Agent
rl_start_time = tic;
fprintf('加载训练好的SAC Agent...\n');

try
    load(sac_agent_file, 'finalAgent', 'cS', 'paramS_for_rl', 'rng_M');
    
    % 检查加载的Agent
    if ~exist('finalAgent', 'var') || ~isa(finalAgent, 'rl.agent.AbstractAgent')
        error('Agent文件中缺少有效的finalAgent');
    end
    
    fprintf('SAC Agent加载成功\n');
    
    % 使用当前测试的宏观参数M_test转换策略矩阵
    fprintf('使用M_test参数转换SAC策略矩阵...\n');
    fprintf('转换参数: R_k=%.3f, w=%.3f, tau_l=%.3f, theta=%.3f\n', ...
        M_test.R_k_net_factor, M_test.w_gross, M_test.tau_l, M_test.theta_payg_actual);
    
    % 调用转换函数 (4个参数版本)
    [cPolM_sac, kPolM_sac, cPpsPolM_choice_sac, ~] = convert_SAC_to_PolicyFunctions(...
        finalAgent, cS, paramS_for_rl, M_test);
    
    rl_time = toc(rl_start_time);
    fprintf('SAC策略矩阵转换完成。耗时: %.2f 秒\n', rl_time);
    fprintf('策略矩阵尺寸: cPol=%s, kPol=%s, cPpsPol=%s\n', ...
        mat2str(size(cPolM_sac)), mat2str(size(kPolM_sac)), mat2str(size(cPpsPolM_choice_sac)));
    
    rl_success = true;
catch ME
    rl_time = toc(rl_start_time);
    fprintf('SAC策略转换失败: %s\n', ME.message);
    fprintf('错误详情: %s\n', getReport(ME));
    rl_success = false;
end

% 初始化RL参数结构（用于状态索引查找）
if rl_success && exist('paramS_for_rl', 'var')
    % 使用从agent文件中加载的参数
    paramS_rl = paramS_for_rl;
    fprintf('使用从Agent文件加载的RL参数\n');
else
    % 重新生成RL参数
    paramS_rl = struct();
    [paramS_rl.leLogGridV, paramS_rl.leTrProbM, paramS_rl.leProb1V] = main_olg_v8_utils.EarningProcess_olgm(cS);
    paramS_rl.leGridV = exp(paramS_rl.leLogGridV(:));
    paramS_rl.ageEffV_new = cS.ageEffV_new;
    fprintf('重新生成RL参数\n');
end

%% 4. 生成随机策略（作为基准比较）
fprintf('\n--- 4. 生成随机策略（基准比较）---\n');
random_start_time = tic;

% 随机策略矩阵（与VFI和RL相同的维度）
cPolM_random = zeros(cS.nk, cS.nkpps, cS.nw, cS.aD_new);
kPolM_random = zeros(cS.nk, cS.nkpps, cS.nw, cS.aD_new);
cPpsPolM_random = zeros(cS.nk, cS.nkpps, cS.nw, cS.aD_new);

fprintf('生成随机策略矩阵...\n');

% 为每个状态生成随机但可行的策略
for age_idx = 1:cS.aD_new
    for i_k = 1:cS.nk
        for i_kpps = 1:cS.nkpps
            for i_eps = 1:cS.nw
                k_state = cS.kGridV(i_k);
                kpps_state = cS.kppsGridV(i_kpps);
                epsilon_val = paramS_rl.leGridV(i_eps);
                
                % 计算当前状态下的资源
                paramS_hh_random.tau_l = M_test.tau_l;
                paramS_hh_random.theta_payg_actual_for_hh = M_test.theta_payg_actual;
                paramS_hh_random.pps_tax_deferral_active = cS.pps_active;
                
                % 计算PAYG福利
                b_payg_this_age = 0;
                if age_idx > cS.aR_new
                    b_payg_this_age = M_test.b_payg_avg_retiree;
                end
                
                % 随机PPS缴费（在允许范围内）
                max_cpps_random = 0;
                if age_idx <= cS.aR_new && cS.pps_active
                    age_efficiency = cS.ageEffV_new(age_idx);
                    gross_labor_income = M_test.w_gross * age_efficiency * epsilon_val;
                    if gross_labor_income > 1e-6
                        max_cpps_by_frac = gross_labor_income * cS.pps_max_contrib_frac;
                        max_cpps_random = min(cS.pps_annual_contrib_limit, max_cpps_by_frac);
                        max_cpps_random = max(0, max_cpps_random);
                    end
                end
                cpps_random = rand() * max_cpps_random;  % 随机选择0到最大值之间
                
                % 计算扣除PPS缴费后的资源
                [resources_after_pps, ~, ~] = main_olg_v8_utils.HHIncome_Huggett(...
                    k_state, M_test.R_k_net_factor, M_test.w_gross, ...
                    M_test.TR_total, b_payg_this_age, cpps_random, ...
                    age_idx, paramS_hh_random, cS, epsilon_val);
                
                % 随机储蓄率（0到1之间）
                consumption_floor_spending = cS.cFloor * (1 + cS.tau_c);
                resources_for_kprime_c = resources_after_pps - consumption_floor_spending;
                
                if resources_for_kprime_c > 0
                    savings_rate = rand();  % 随机储蓄率
                    k_next_random = savings_rate * resources_for_kprime_c;
                    k_next_random = max(cS.kMin, min(k_next_random, cS.kMax));
                else
                    k_next_random = cS.kMin;
                end
                
                % 计算对应的消费
                consumption_expenditure = resources_after_pps - k_next_random;
                c_random = max(cS.cFloor, consumption_expenditure / (1 + cS.tau_c));
                
                % 存储随机策略
                cPolM_random(i_k, i_kpps, i_eps, age_idx) = c_random;
                kPolM_random(i_k, i_kpps, i_eps, age_idx) = k_next_random;
                cPpsPolM_random(i_k, i_kpps, i_eps, age_idx) = cpps_random;
            end
        end
    end
end

random_time = toc(random_start_time);
fprintf('随机策略生成完成。耗时: %.2f 秒\n', random_time);

%% 5. 模拟生命周期轨迹比较（三种策略，期望化死亡率）
if rl_success
    fprintf('\n--- 5. 模拟生命周期轨迹比较 (期望化死亡率) ---\n');
    
    nSim = 100; % 模拟个体数量
    fprintf('模拟 %d 个个体的完整生命周期轨迹 (固定%d期)...\n', nSim, cS.aD_new);
    
    % 初始化结果存储
    lifetime_utility_vfi = zeros(nSim, 1);
    lifetime_utility_rl = zeros(nSim, 1);
    lifetime_utility_random = zeros(nSim, 1);
    
    % 生命周期轨迹存储
    age_path = 1:cS.aD_new;
    k_path_vfi = zeros(nSim, cS.aD_new);
    k_path_rl = zeros(nSim, cS.aD_new);
    k_path_random = zeros(nSim, cS.aD_new);
    c_path_vfi = zeros(nSim, cS.aD_new);
    c_path_rl = zeros(nSim, cS.aD_new);
    c_path_random = zeros(nSim, cS.aD_new);
    cpps_path_vfi = zeros(nSim, cS.aD_new);
    cpps_path_rl = zeros(nSim, cS.aD_new);
    cpps_path_random = zeros(nSim, cS.aD_new);
    
    sim_start_time = tic;
    
    for i_sim = 1:nSim
        if mod(i_sim, 20) == 0
            fprintf('  进度: %d/%d\n', i_sim, nSim);
        end
        
        % 初始状态（三种方法使用相同的初始状态）
        k_current_vfi = cS.kMin;
        kpps_current_vfi = cS.kppsMin;
        k_current_rl = cS.kMin;
        kpps_current_rl = cS.kppsMin;
        k_current_random = cS.kMin;
        kpps_current_random = cS.kppsMin;
        eps_idx_current = find(rand() <= cumsum(paramS_vfi.leProb1V), 1, 'first');
        if isempty(eps_idx_current), eps_idx_current = 1; end
        
        utility_sum_vfi = 0;
        utility_sum_rl = 0;
        utility_sum_random = 0;
        
        for age_idx = 1:cS.aD_new
            % VFI策略 - 使用VFI的当前状态
            [k_idx_vfi, kpps_idx_vfi, eps_idx] = get_state_indices(k_current_vfi, kpps_current_vfi, eps_idx_current, cS, paramS_vfi);
            
            c_vfi = cPolM_vfi(k_idx_vfi, kpps_idx_vfi, eps_idx, age_idx);
            k_next_vfi = kPolM_vfi(k_idx_vfi, kpps_idx_vfi, eps_idx, age_idx);
            
            if size(cPpsPolM_vfi, 4) >= age_idx
                cpps_vfi = cPpsPolM_vfi(k_idx_vfi, kpps_idx_vfi, eps_idx, age_idx);
            else
                cpps_vfi = 0;
            end
            
            % RL策略 - 使用RL的当前状态
            [k_idx_rl, kpps_idx_rl, eps_idx] = get_state_indices(k_current_rl, kpps_current_rl, eps_idx_current, cS, paramS_rl);
            
            c_rl = cPolM_sac(k_idx_rl, kpps_idx_rl, eps_idx, age_idx);
            k_next_rl = kPolM_sac(k_idx_rl, kpps_idx_rl, eps_idx, age_idx);
            
            if size(cPpsPolM_choice_sac, 4) >= age_idx
                cpps_rl = cPpsPolM_choice_sac(k_idx_rl, kpps_idx_rl, eps_idx, age_idx);
            else
                cpps_rl = 0;
            end
            
            % 随机策略 - 使用随机策略的当前状态
            [k_idx_random, kpps_idx_random, eps_idx] = get_state_indices(k_current_random, kpps_current_random, eps_idx_current, cS, paramS_rl);
            
            c_random = cPolM_random(k_idx_random, kpps_idx_random, eps_idx, age_idx);
            k_next_random = kPolM_random(k_idx_random, kpps_idx_random, eps_idx, age_idx);
            
            if size(cPpsPolM_random, 4) >= age_idx
                cpps_random = cPpsPolM_random(k_idx_random, kpps_idx_random, eps_idx, age_idx);
            else
                cpps_random = 0;
            end
            
            % 存储轨迹
            k_path_vfi(i_sim, age_idx) = k_current_vfi;
            k_path_rl(i_sim, age_idx) = k_current_rl;
            k_path_random(i_sim, age_idx) = k_current_random;
            c_path_vfi(i_sim, age_idx) = c_vfi;
            c_path_rl(i_sim, age_idx) = c_rl;
            c_path_random(i_sim, age_idx) = c_random;
            cpps_path_vfi(i_sim, age_idx) = cpps_vfi;
            cpps_path_rl(i_sim, age_idx) = cpps_rl;
            cpps_path_random(i_sim, age_idx) = cpps_random;
            
            % 计算当期效用
            [~, u_vfi] = main_olg_v8_utils.CES_utility(c_vfi, cS.sigma, cS);
            [~, u_rl] = main_olg_v8_utils.CES_utility(c_rl, cS.sigma, cS);
            [~, u_random] = main_olg_v8_utils.CES_utility(c_random, cS.sigma, cS);
            
            % 累积折现效用
            discount_factor = cS.beta^(age_idx - 1);
            utility_sum_vfi = utility_sum_vfi + discount_factor * u_vfi;
            utility_sum_rl = utility_sum_rl + discount_factor * u_rl;
            utility_sum_random = utility_sum_random + discount_factor * u_random;
            
            % 独立更新状态 - 关键修改！
            k_current_vfi = k_next_vfi;  % VFI使用自己的状态演化
            k_current_rl = k_next_rl;    % RL使用自己的状态演化
            k_current_random = k_next_random; % 随机策略使用自己的状态演化
            
            % PPS资产演化 (简化处理，三种方法类似)
            % 这里可以进一步细化，但由于PPS演化规则复杂，暂时简化
            if cS.pps_active
                % 基础PPS演化（简化版本）
                pps_return_factor = 1 + ((M_test.R_k_net_factor - 1) + cS.pps_return_rate_premium);
                kpps_current_vfi = (kpps_current_vfi + cpps_vfi) * pps_return_factor;
                kpps_current_rl = (kpps_current_rl + cpps_rl) * pps_return_factor;
                kpps_current_random = (kpps_current_random + cpps_random) * pps_return_factor;
                
                % 限制在合理范围内
                kpps_current_vfi = max(cS.kppsMin, min(cS.kppsMax, kpps_current_vfi));
                kpps_current_rl = max(cS.kppsMin, min(cS.kppsMax, kpps_current_rl));
                kpps_current_random = max(cS.kppsMin, min(cS.kppsMax, kpps_current_random));
            end
            
            % 效率冲击演化（三种方法使用相同的冲击序列）
            if age_idx < cS.aD_new
                trans_probs = paramS_vfi.leTrProbM(eps_idx_current, :);
                eps_idx_current = find(rand() <= cumsum(trans_probs), 1, 'first');
                if isempty(eps_idx_current), eps_idx_current = eps_idx; end
            end
            
            % 期望化死亡率处理：移除随机死亡，确保与RL训练环境一致
            % 所有个体都完整生活到最大年龄 cS.aD_new
            % 死亡率通过VFI的Bellman方程期望值和RL的价值函数学习自动纳入
        end
        
        lifetime_utility_vfi(i_sim) = utility_sum_vfi;
        lifetime_utility_rl(i_sim) = utility_sum_rl;
        lifetime_utility_random(i_sim) = utility_sum_random;
    end
    
    sim_time = toc(sim_start_time);
    fprintf('生命周期模拟完成。耗时: %.2f 秒\n', sim_time);
    
    %% 6. 结果分析和比较
    fprintf('\n--- 6. 结果分析和比较 (期望化死亡率版本) ---\n');
    fprintf('注意：本次比较使用期望化死亡率，确保与RL训练环境一致\n');
    fprintf('     - 所有个体都完整生活%d期\n', cS.aD_new);
    fprintf('     - VFI和RL策略都基于相同的宏观参数M_test生成\n');
    fprintf('     - VFI和RL各自使用自己的最优决策路径\n');
    fprintf('     - 死亡率通过期望值纳入，无随机终止\n\n');
    
    % 统计摘要
    mean_utility_vfi = mean(lifetime_utility_vfi);
    mean_utility_rl = mean(lifetime_utility_rl);
    mean_utility_random = mean(lifetime_utility_random);
    std_utility_vfi = std(lifetime_utility_vfi);
    std_utility_rl = std(lifetime_utility_rl);
    std_utility_random = std(lifetime_utility_random);
    
    % 分析路径差异
    mean_k_diff = mean(mean(k_path_rl - k_path_vfi, 1));
    mean_c_diff = mean(mean(c_path_rl - c_path_vfi, 1));
    mean_cpps_diff = mean(mean(cpps_path_rl - cpps_path_vfi, 1));
    
    fprintf('生涯总折现效用比较 (基于 %d 个个体):\n', nSim);
    fprintf('  VFI方法:\n');
    fprintf('    均值: %.4f\n', mean_utility_vfi);
    fprintf('    标准差: %.4f\n', std_utility_vfi);
    fprintf('  RL方法 (SAC):\n');
    fprintf('    均值: %.4f\n', mean_utility_rl);
    fprintf('    标准差: %.4f\n', std_utility_rl);
    fprintf('  随机策略:\n');
    fprintf('    均值: %.4f\n', mean_utility_random);
    fprintf('    标准差: %.4f\n', std_utility_random);
    
    fprintf('\n路径差异分析 (RL - VFI 的平均值):\n');
    fprintf('  平均资产差异: %.4f\n', mean_k_diff);
    fprintf('  平均消费差异: %.4f\n', mean_c_diff);
    fprintf('  平均PPS缴费差异: %.4f\n', mean_cpps_diff);
    
    utility_diff_vfi_rl = mean_utility_rl - mean_utility_vfi;
    utility_improvement_pct_vfi_rl = (utility_diff_vfi_rl / abs(mean_utility_vfi)) * 100;
    
    utility_diff_vfi_random = mean_utility_random - mean_utility_vfi;
    utility_improvement_pct_vfi_random = (utility_diff_vfi_random / abs(mean_utility_vfi)) * 100;
    
    fprintf('\n效用差异分析:\n');
    fprintf('  RL - VFI = %.4f\n', utility_diff_vfi_rl);
    fprintf('  相对改进: %.2f%%\n', utility_improvement_pct_vfi_rl);
    
    fprintf('  随机策略 - VFI = %.4f\n', utility_diff_vfi_random);
    fprintf('  相对改进: %.2f%%\n', utility_improvement_pct_vfi_random);
    
    if utility_diff_vfi_rl > 0
        fprintf('  >>> RL方法表现更好！\n');
    elseif utility_diff_vfi_rl < 0
        fprintf('  >>> VFI方法表现更好！\n');
    else
        fprintf('  >>> 两种方法表现相当。\n');
    end
    
    % 统计显著性检验
    [h_rl_vfi, p_value_rl_vfi] = ttest2(lifetime_utility_rl, lifetime_utility_vfi);
    [h_random_vfi, p_value_random_vfi] = ttest2(lifetime_utility_random, lifetime_utility_vfi);
    fprintf('\n统计显著性检验 (双样本t检验):\n');
    fprintf('  RL - VFI: p值: %.4f\n', p_value_rl_vfi);
    fprintf('  随机策略 - VFI: p值: %.4f\n', p_value_random_vfi);
    % fprintf('  显著性水平0.05下的结果: %s\n', ...
    %     h ? '差异显著' : '差异不显著');
    
    %% 7. 可视化比较
    fprintf('\n--- 7. 生成比较图表 ---\n');
    
    % 图1：生涯效用分布比较
    figure('Name', 'VFI vs RL: 生涯效用分布比较', 'Position', [100, 100, 1200, 400]);
    
    subplot(1, 3, 1);
    histogram(lifetime_utility_vfi, 20, 'DisplayName', 'VFI');
    hold on;
    histogram(lifetime_utility_rl, 20, 'DisplayName', 'RL (SAC)');
    histogram(lifetime_utility_random, 20, 'DisplayName', '随机策略');
    xlabel('生涯总折现效用');
    ylabel('频数');
    title('效用分布比较');
    legend('show');
    grid on;
    
    subplot(1, 3, 2);
    boxplot([lifetime_utility_vfi, lifetime_utility_rl, lifetime_utility_random], 'Labels', {'VFI', 'RL', '随机策略'});
    ylabel('生涯总折现效用');
    title('效用分布箱线图');
    grid on;
    
    subplot(1, 3, 3);
    scatter(lifetime_utility_vfi, lifetime_utility_rl, 'filled', 'DisplayName', 'RL vs VFI');
    hold on;
    scatter(lifetime_utility_vfi, lifetime_utility_random, 'filled', 'DisplayName', '随机策略 vs VFI');
    
    % 绘制45度线
    min_val = min([lifetime_utility_vfi; lifetime_utility_rl; lifetime_utility_random]);
    max_val = max([lifetime_utility_vfi; lifetime_utility_rl; lifetime_utility_random]);
    plot([min_val, max_val], [min_val, max_val], 'r--', 'LineWidth', 2, 'DisplayName', '45度线');
    
    xlabel('VFI 生涯效用');
    ylabel('其他策略 生涯效用');
    title('个体效用对比');
    legend('RL', '随机策略', '45度线', 'Location', 'best');
    grid on;
    
    % 图2：平均生命周期轨迹比较
    figure('Name', 'VFI vs RL: 平均生命周期轨迹', 'Position', [150, 150, 1200, 800]);
    
    mean_k_vfi = mean(k_path_vfi, 1);
    mean_k_rl = mean(k_path_rl, 1);
    mean_k_random = mean(k_path_random, 1);
    mean_c_vfi = mean(c_path_vfi, 1);
    mean_c_rl = mean(c_path_rl, 1);
    mean_c_random = mean(c_path_random, 1);
    mean_cpps_vfi = mean(cpps_path_vfi, 1);
    mean_cpps_rl = mean(cpps_path_rl, 1);
    mean_cpps_random = mean(cpps_path_random, 1);
    
    subplot(2, 2, 1);
    plot(age_path, mean_k_vfi, 'b-', 'LineWidth', 2, 'DisplayName', 'VFI');
    hold on;
    plot(age_path, mean_k_rl, 'r--', 'LineWidth', 2, 'DisplayName', 'RL (SAC)');
    plot(age_path, mean_k_random, 'g--', 'LineWidth', 2, 'DisplayName', '随机策略');
    xlabel('年龄组');
    ylabel('平均非PPS资产');
    title('生命周期资产积累');
    legend('show');
    grid on;
    
    subplot(2, 2, 2);
    plot(age_path, mean_c_vfi, 'b-', 'LineWidth', 2, 'DisplayName', 'VFI');
    hold on;
    plot(age_path, mean_c_rl, 'r--', 'LineWidth', 2, 'DisplayName', 'RL (SAC)');
    plot(age_path, mean_c_random, 'g--', 'LineWidth', 2, 'DisplayName', '随机策略');
    xlabel('年龄组');
    ylabel('平均消费');
    title('生命周期消费');
    legend('show');
    grid on;
    
    subplot(2, 2, 3);
    plot(age_path, mean_cpps_vfi, 'b-', 'LineWidth', 2, 'DisplayName', 'VFI');
    hold on;
    plot(age_path, mean_cpps_rl, 'r--', 'LineWidth', 2, 'DisplayName', 'RL (SAC)');
    plot(age_path, mean_cpps_random, 'g--', 'LineWidth', 2, 'DisplayName', '随机策略');
    xlabel('年龄组');
    ylabel('平均PPS缴费');
    title('生命周期PPS缴费');
    legend('show');
    grid on;
    
    subplot(2, 2, 4);
    utility_diff_individual_rl = lifetime_utility_rl - lifetime_utility_vfi;
    utility_diff_individual_random = lifetime_utility_random - lifetime_utility_vfi;
    
    plot(1:nSim, utility_diff_individual_rl, 'ro', 'MarkerSize', 4, 'DisplayName', 'RL - VFI');
    hold on;
    plot(1:nSim, utility_diff_individual_random, 'go', 'MarkerSize', 4, 'DisplayName', '随机策略 - VFI');
    yline(0, 'k--', 'LineWidth', 2, 'DisplayName', '零线');
    yline(mean(utility_diff_individual_rl), 'r-', 'LineWidth', 2, 'DisplayName', 'RL平均差异');
    yline(mean(utility_diff_individual_random), 'g-', 'LineWidth', 2, 'DisplayName', '随机策略平均差异');
    xlabel('个体编号');
    ylabel('效用差异 (其他策略 - VFI)');
    title('个体效用差异');
    legend('show', 'Location', 'best');
    grid on;
    
    %% 8. 保存结果
    fprintf('\n--- 8. 保存比较结果 ---\n');
    
    comparison_results = struct();
    comparison_results.M_test = M_test;
    comparison_results.nSim = nSim;
    comparison_results.vfi_time = vfi_time;
    comparison_results.rl_time = rl_time;
    comparison_results.mean_utility_vfi = mean_utility_vfi;
    comparison_results.mean_utility_rl = mean_utility_rl;
    comparison_results.mean_utility_random = mean_utility_random;
    comparison_results.std_utility_vfi = std_utility_vfi;
    comparison_results.std_utility_rl = std_utility_rl;
    comparison_results.std_utility_random = std_utility_random;
    comparison_results.utility_diff_vfi_rl = utility_diff_vfi_rl;
    comparison_results.utility_improvement_pct_vfi_rl = utility_improvement_pct_vfi_rl;
    comparison_results.utility_diff_vfi_random = utility_diff_vfi_random;
    comparison_results.utility_improvement_pct_vfi_random = utility_improvement_pct_vfi_random;
    comparison_results.p_value_rl_vfi = p_value_rl_vfi;
    comparison_results.p_value_random_vfi = p_value_random_vfi;
    comparison_results.is_significant_rl_vfi = h_rl_vfi;
    comparison_results.is_significant_random_vfi = h_random_vfi;
    comparison_results.lifetime_utility_vfi = lifetime_utility_vfi;
    comparison_results.lifetime_utility_rl = lifetime_utility_rl;
    comparison_results.lifetime_utility_random = lifetime_utility_random;
    comparison_results.k_path_vfi = k_path_vfi;
    comparison_results.k_path_rl = k_path_rl;
    comparison_results.k_path_random = k_path_random;
    comparison_results.c_path_vfi = c_path_vfi;
    comparison_results.c_path_rl = c_path_rl;
    comparison_results.c_path_random = c_path_random;
    comparison_results.cpps_path_vfi = cpps_path_vfi;
    comparison_results.cpps_path_rl = cpps_path_rl;
    comparison_results.cpps_path_random = cpps_path_random;
    
    save('vfi_vs_rl_comparison_results.mat', 'comparison_results');
    fprintf('比较结果已保存到 vfi_vs_rl_comparison_results.mat\n');
    
else
    fprintf('\n由于SAC策略转换失败，无法进行详细比较。\n');
end

%% 9. 计算性能总结
fprintf('\n--- 9. 性能总结 ---\n');
fprintf('计算时间比较:\n');
fprintf('  VFI求解时间: %.2f 秒\n', vfi_time);
fprintf('  随机策略生成时间: %.2f 秒\n', random_time);
if rl_success
    fprintf('  RL策略加载时间: %.2f 秒\n', rl_time);
    % fprintf('  RL相对速度: %.2fx %s\n', vfi_time/rl_time, ...
    %     vfi_time > rl_time ? '(更慢)' : '(更快)');
end

fprintf('\n=== 比较分析完成 ===\n');

%%% 辅助函数
function [k_idx, kpps_idx, eps_idx] = get_state_indices(k_val, kpps_val, eps_idx_val, cS, paramS)
    % 获取状态变量在网格中的索引
    
    % 非PPS资产索引
    [~, k_idx] = min(abs(cS.kGridV - k_val));
    
    % PPS资产索引 - 需要先检查kppsGridV是否存在
    if isfield(cS, 'kppsGridV')
        [~, kpps_idx] = min(abs(cS.kppsGridV - kpps_val));
    else
        % 如果kppsGridV不存在，创建一个简单的网格
        kppsMax = max(10, cS.kMax/2); % 简单估计PPS资产上界
        kppsGridV_temp = linspace(cS.kppsMin, kppsMax, cS.nkpps);
        [~, kpps_idx] = min(abs(kppsGridV_temp - kpps_val));
    end
    
    % 效率索引
    eps_idx = eps_idx_val;
    
    % 确保索引在有效范围内
    k_idx = max(1, min(k_idx, cS.nk));
    kpps_idx = max(1, min(kpps_idx, cS.nkpps));
    eps_idx = max(1, min(eps_idx, cS.nw));
end
% --- END OF FILE compare_rl_and_vfi.m --- 