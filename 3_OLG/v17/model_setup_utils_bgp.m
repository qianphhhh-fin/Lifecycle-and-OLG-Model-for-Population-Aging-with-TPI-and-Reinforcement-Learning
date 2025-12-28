classdef model_setup_utils_bgp
    methods (Static)

        function cS = ParameterValues()
            % =========================================================================
            % == 函数: ParameterValues
            % == 版本: [v_dynamic_mortality - 动态死亡率版]
            % ==
            % == 核心修改:
            % ==   - 移除硬编码的`raw_mortality_data`。
            % ==   - 调用新的`load_mortality_data`辅助函数，从外部Excel文件加载时变的死亡率数据。
            % ==   - 将完整的死亡率路径存入cS，并为当前计算选择一个基准年份的死亡率，以保持向后兼容性。
            % =========================================================================

            % --- 核心时间与年龄参数 ---
            cS.time_Step = 5;                   % 模型中每个时期代表的年数 (例如，5年为一个年龄组)。
            cS.age1_orig = 20;                  % 个体进入模型的起始物理年龄 (例如，20岁)。
            cS.ageLast_orig = 98;               % 个体确定性死亡前的最大物理年龄 (例如，98岁)。
            cS.ageRetire_orig = 60;             % 个体的法定退休物理年龄 (例如，60岁)。

            % --- 行为参数 (家庭偏好) ---
            cS.beta = 0.98;                    % 主观折现因子 (年化)，反映了个体的耐心程度。模型中将根据time_Step进行调整。
            cS.sigma = 5;                       % 跨期替代弹性的倒数 (相对风险厌恶系数)，控制消费平滑的意愿。
            cS.phi_bequest = 5;                 % 遗赠效用权重，衡量“暖心式”(warm-glow)遗赠动机的强度。
            cS.cFloor = 0.005;                  % 最低消费水平，用于保证效用函数在消费趋近于零时有界。

            % --- 生产技术参数 ---
            cS.alpha = 0.35;                    % 产出中私人资本的份额 (资本产出弹性)。
            cS.gamma = 0.12;                    % 产出中公共资本的份额 (公共资本产出弹性)。
            cS.ddk = 1 - (1 - 0.015)^cS.time_Step; % 私人资本的折旧率 (模型期)，由年化率0.015转换而来。
            cS.ddk_g = cS.ddk;                  % 公共资本的折旧率 (模型期)，为简化设为与私人资本相同。
            cS.A = 1.0;                         % 全要素生产率(TFP)的初始基准水平，通常归一化为1。

            % --- 政府与政策参数 ---
            cS.tau_k = 0.02;                    % 资本利得税率。
            cS.tau_l = 0.06;                    % 劳动收入税率 (在PAYG缴费之外)。
            cS.tau_c = 0.03;                    % 消费税率。
            % --- 步骤 2.1b: 设定稳态下的目标政府债务/GDP比率 ---

            % --- 长期BGP(平衡增长路径)参数 ---
            cS.g_A_ss = 0.012;                  % 长期稳态下技术进步的年化增长率 (TFP growth rate)。
            cS.long_term_birth_rate_target_per_1000 = 3.5; % 长期稳态目标的年化粗出生率 (每千人)。
            cS.transition_period_years = 10;    % 外生路径从数据点过渡到长期稳态所需的年数。
            cS.curvature_param = 0.5;           % 过渡路径的曲率参数 (k<1为上凸收敛，k>1为下凸收敛)。
            cS.cbr_last_year = 2040;            % 联合国出生率数据使用的最后年份

            % --- 养老金体系参数 ---
            cS.pps_active = false;              % 控制是否激活私人养老金(PPS)模块的全局开关。
            cS.pps_contrib_rate = 0.03;         % PPS的强制缴费率 (占劳动收入的比例)。
            cS.pps_withdrawal_rate = 0.10;      % 退休后每个模型期强制提取PPS资产的比例。
            cS.pps_tax_rate_withdrawal = 0.03;  % 提取PPS资产时适用的税率。
            cS.theta_path = 0.06;               % 外生设定的PAYG养老金缴费率 (在DC模式下的基准值)。

            % --- 收入过程参数 ---
            cS.nw = 1;   % 5                       % 劳动效率冲击离散状态的数量。
            % 年龄节点
            cS.young_shock_start_age = 26;
            cS.young_shock_end_age = 35;
            cS.old_shock_start_age = 65;
            cS.old_shock_end_age = 88;
            % 冲击年化概率峰值
            cS.p_shock_young_peak_annual = 0; %0.01;
            cS.p_shock_old_peak_annual = 0; %0.01;
            % 冲击支出比例 (kappa)
            cS.kappa_young = 0.9;
            cS.kappa_old = 0.9;

            cS = model_setup_utils_bgp.generate_mortality_path(cS); % 生成死亡率数据

            % 计算年龄效率剖面
            ages_eff = (cS.age1_orig : cS.ageLast_orig)';              % 获取用于计算效率的物理年龄向量。
            beta0=-13.215788; beta1=1.349514; beta2=-0.043363; beta3=0.000585; beta4=-0.000003; % 效率函数的多项式系数。
            log_age_eff_orig = beta0 + beta1*ages_eff + beta2*(ages_eff.^2) + beta3*(ages_eff.^3) + beta4*(ages_eff.^4);
            ageEffV_orig_unnormalized = exp(log_age_eff_orig);          % 计算未归一化的年龄效率。
            ageEffV_orig_unnormalized((cS.aR_idx_orig):end) = 0;        % 退休后效率为0。
            mean_efficiency_working = mean(ageEffV_orig_unnormalized(1:(cS.aR_idx_orig - 1))); % 计算工作期的平均效率用于归一化。
            ageEffV_orig = ageEffV_orig_unnormalized / mean_efficiency_working; % 归一化，使得工作期平均效率为1。
            cS.ageEffV_new = zeros(cS.aD_new, 1);                       % 初始化模型年龄组的平均效率向量。
            for a = 1:cS.aD_new
                phys_ages_in_group = cS.physAgeMap{a};
                model_ages_indices = phys_ages_in_group - cS.age1_orig + 1;
                cS.ageEffV_new(a) = mean(ageEffV_orig(model_ages_indices));
            end
        end



        function [Z_out,Z_out_raw, A_path, cS] = generate_exo_paths(cS, graph_flag)
            % =========================================================================
            % == 函数: generate_exo_paths
            % == 版本: [v_fully_aligned_paths - 完全对齐路径版]
            % ==
            % == 核心修改:
            % ==   1. [路径完全对齐] 不仅Z_out和A_path，现在cS.s_path_all和cS.s_pathV的
            % ==      时间维度(列数)也与内生决定的模拟期数T_sim完全对齐。
            % ==   2. [即时生成] 在人口模拟循环中即时生成并存储每期使用的存活率向量。
            % ==   3. [最终裁剪] 在确定T_sim后，对所有路径矩阵进行统一裁剪，确保维度一致。
            % =========================================================================

            fprintf('\n--- 启动外生路径生成器 (v_fully_aligned_paths) ---\n');

            % --- 1. 初始化模拟参数 ---
            ss_convergence_tol = 1e-6;
            max_sim_periods = 200;

            % --- 2. 加载并处理【出生率】数据 ---
            try 
                cbr_data = readtable('data\人口\UN_PPP2024_CBR_birthper1000_China.xlsx');
            catch ME
                error('无法加载出生率数据文件: %s', ME.message);
            end
            var_names_cbr = cbr_data.Properties.VariableNames;
            year_cols_indices_cbr = find(startsWith(var_names_cbr, 'y'));
            un_years = str2double(cellfun(@(x) x(2:end), var_names_cbr(year_cols_indices_cbr), 'UniformOutput', false));
            un_cbr = table2array(cbr_data(1, year_cols_indices_cbr)) / 1000;

            % --- 3. 构建【年度出生率】路径并确定恒定规则 ---
            last_data_year = cS.cbr_last_year;
            fprintf('   人口相关路径将在 %d 年后保持恒定。\n', last_data_year);

            annual_years_full_range = cS.start_year:(cS.start_year + (max_sim_periods * cS.time_Step) -1);
            birth_rate_path_annual_full = zeros(size(annual_years_full_range));
            data_mask = annual_years_full_range <= last_data_year;
            birth_rate_path_annual_full(data_mask) = interp1(un_years, un_cbr, annual_years_full_range(data_mask), 'pchip', 'extrap');
            last_data_rate = birth_rate_path_annual_full(find(data_mask, 1, 'last'));
            ss_mask = annual_years_full_range > last_data_year;
            birth_rate_path_annual_full(ss_mask) = last_data_rate;
            
            % --- 4. 初始化人口分布 ---
            try
                pop_data_by_age = readtable('data\人口\china_population_by_age_headerFix.xlsx');
            catch ME
                error('无法加载分年龄人口数据文件: %s', ME.message); 
            end
            initial_pop_row = pop_data_by_age(pop_data_by_age.Year == cS.ss0_year, :);
            if isempty(initial_pop_row), error('在人口数据中找不到初始年份: %d', cS.ss0_year); end
            num_age_groups_all = size(cS.s_path_all, 1);
            initial_pop_aggregated = model_setup_utils_bgp.aggregate_population(initial_pop_row, num_age_groups_all, cS.time_Step);
            
            % --- 5. 执行人口动态前向模拟 ---
            fprintf('   正在启动人口动态模拟以确定 T_sim ...\n');
            Z_path_raw = zeros(num_age_groups_all, max_sim_periods);
            Z_path_raw(:, 1) = initial_pop_aggregated;
            Z_norm_prev = Z_path_raw(:, 1) / sum(Z_path_raw(:, 1));
            
            % [核心修改] 预分配用于存储模拟期间存活率路径的矩阵
            s_path_all_sim = zeros(num_age_groups_all, max_sim_periods);
            s_pathV_sim = zeros(cS.aD_new, max_sim_periods);
            
            converged = false;
            T_final = max_sim_periods;

            for t = 1:max_sim_periods
                % a. 确定当前模拟期 t 对应的日历年份
                current_year = cS.start_year + (t - 1) * cS.time_Step;
    
                % b. 获取当前年份的年度死亡率，并计算对应的年度存活率
                if current_year <= cS.mortality_years(end)
                    raw_survival_data_t = 1 - interp1(cS.mortality_years, cS.mortality_path_data, current_year);
                else
                    raw_survival_data_t = 1 - cS.mortality_path_data(end,:); % 超出数据年份，使用最后一年数据
                end
                raw_survival_data_t(end) = 0; % 确保100岁之后存活率为0

                % c. 计算并存储本期的 s_path_all 和 s_pathV
                % 计算 s_path_all_t
                s_path_all_t = zeros(num_age_groups_all, 1);
                for a = 1:(num_age_groups_all - 1)
                    s_period_prod = 1.0;
                    start_age = (a-1) * cS.time_Step;
                    for i = 1:cS.time_Step
                        age = start_age + i - 1;
                        if (age + 1) <= length(raw_survival_data_t)
                            s_period_prod = s_period_prod * raw_survival_data_t(age + 1);
                        else
                            s_period_prod = 0; break;
                        end
                    end
                    s_path_all_t(a) = s_period_prod;
                end
                s_path_all_sim(:, t) = s_path_all_t;

                % 计算 s_pathV_t
                s_pathV_t = zeros(cS.aD_new, 1);
                for a = 1:(cS.aD_new - 1)
                    s_period_prod = 1.0;
                    phys_ages_in_group = cS.physAgeMap{a};
                    for age_idx = 1:length(phys_ages_in_group)
                        age = phys_ages_in_group(age_idx);
                        if (age + 1) <= length(raw_survival_data_t)
                            s_period_prod = s_period_prod * raw_survival_data_t(age + 1);
                        else
                            s_period_prod = 0; break;
                        end
                    end
                    s_pathV_t(a) = s_period_prod;
                end
                s_pathV_sim(:, t) = s_pathV_t;

                % d. 模拟下一期人口 (如果还在最大模拟期内)
                if t < max_sim_periods
                    Z_path_raw(2:end, t+1) = Z_path_raw(1:end-1, t) .* s_path_all_sim(1:end-1, t);
    
                    total_pop_base = sum(Z_path_raw(:, t));
                    period_start_idx_in_annual = (t-1)*cS.time_Step + 1;
                    period_end_idx_in_annual = t*cS.time_Step;
                    avg_birth_rate_period = mean(birth_rate_path_annual_full(period_start_idx_in_annual:period_end_idx_in_annual));
                    new_entrants = total_pop_base * avg_birth_rate_period * cS.time_Step;
                    Z_path_raw(1, t+1) = new_entrants;
    
                    Z_norm_current = Z_path_raw(:, t+1) / sum(Z_path_raw(:, t+1));
                    diff = max(abs(Z_norm_current - Z_norm_prev));
                    if ~converged && diff < ss_convergence_tol
                        T_final = t + 1; converged = true;
                        fprintf('   ✅ 完整人口结构在第 %d 期 (年份: %d) 达到稳态。\n', T_final, cS.start_year + t*cS.time_Step);
                    end
                    Z_norm_prev = Z_norm_current;
                end
            end
            if ~converged, T_final = max_sim_periods; warning('人口模拟在达到最大期数 %d 后仍未收敛！', max_sim_periods); end

            % --- 6. [最终对齐] 裁剪所有路径，使其长度等于T_final ---
            cS.T_sim = T_final;
            cS.end_year = cS.start_year + (T_final - 1) * cS.time_Step;
            fprintf('   内生决定的模拟期数 T_sim = %d, 结束年份 = %d\n', cS.T_sim, cS.end_year);

            % a. 裁剪人口路径并格式化Z_out
            Z_path_raw_final = Z_path_raw(:, 1:cS.T_sim);
            model_age_start_index = floor(cS.age1_orig / cS.time_Step) + 1;
            model_age_end_index = model_age_start_index + cS.aD_new - 1;
            Z_out_raw = Z_path_raw_final(model_age_start_index : model_age_end_index, :);
            Z_out = Z_out_raw ./ sum(Z_out_raw, 1);
            if size(Z_out, 1) ~= cS.aD_new, warning('Z_out的行数(%d)与cS.aD_new(%d)不匹配！', size(Z_out, 1), cS.aD_new); end

            % b. [核心修改] 裁剪存活率路径矩阵
            cS.s_path_all = s_path_all_sim(:, 1:cS.T_sim);
            cS.s_pathV = s_pathV_sim(:, 1:cS.T_sim);
            fprintf('   ✅ 已将存活率路径 cS.s_path_all 和 cS.s_pathV 的长度裁剪为 %d 期。\n', cS.T_sim);

            % c. 计算人口增长率路径
            total_pop_path = sum(Z_path_raw_final, 1);
            pop_growth_rate_path_period = (total_pop_path(2:end) ./ total_pop_path(1:end-1)) - 1;
            pop_growth_rate_path_annual = (1 + pop_growth_rate_path_period).^(1/cS.time_Step) - 1;
            if converged && ~isempty(pop_growth_rate_path_annual), cS.n_ss = pop_growth_rate_path_annual(end); else, cS.n_ss = 0; end

            % --- 7. 生成TFP路径并进行可视化 ---
            model_sim_years = cS.start_year:cS.time_Step:cS.end_year;
            annual_years_vec = cS.start_year:cS.end_year;
            A_path = model_setup_utils_bgp.generate_tfp_path(cS, annual_years_vec, model_sim_years); % A_path的长度已由model_sim_years决定，自动对齐
            
            if graph_flag
                final_annual_mask = annual_years_full_range >= cS.start_year & annual_years_full_range <= cS.end_year;
                birth_rate_path_annual_for_plot = birth_rate_path_annual_full(final_annual_mask);
                if length(model_sim_years) >= 2
                    A_path_annual_for_plot = interp1(model_sim_years, A_path, annual_years_vec, 'pchip');
                else
                    A_path_annual_for_plot = A_path;
                end
                model_setup_utils_bgp.create_exogenous_dynamics_plots(cS, Z_out, A_path_annual_for_plot, ...
                    birth_rate_path_annual_for_plot, annual_years_vec, pop_growth_rate_path_annual);
            end
        end        
        
        function pop_aggregated = aggregate_population(pop_row, num_age_groups, time_step)
            pop_aggregated = zeros(num_age_groups, 1);
            for a = 1:num_age_groups
                age_start = (a-1) * time_step;
                age_end = a * time_step - 1;

                current_sum = 0;
                for age = age_start:age_end
                    col_name = sprintf('age_%d', age);
                    if ismember(col_name, pop_row.Properties.VariableNames)
                        current_sum = current_sum + pop_row.(col_name);
                    elseif age >= 100 % 特殊处理100+岁人口
                        col_name_100plus = 'age_100';
                        if ismember(col_name_100plus, pop_row.Properties.VariableNames)
                            current_sum = current_sum + pop_row.(col_name_100plus);
                        end
                        break; % 100+的人口已经全部加总，跳出内层循环
                    end
                end
                pop_aggregated(a) = current_sum;
            end
        end

        function cS = generateGrids(cS, varargin)
            % [BGP修改] 支持可选的网格上限参数，实现自适应网格
            % [BGP修改] 修改默认值计算逻辑，适应稳态化模型
            % 新签名: generateGrids(cS, 'k_max', value, 'kpps_max', value)
            % 这将网格生成从固定参数变为动态服务

            % 使用 inputParser 处理可选参数
            p = inputParser;

            % [BGP修改] 简化默认值计算逻辑
            % 在稳态化求解中，这些默认值几乎不会被用到，因为求解器会根据猜测的
            % 标准化资本 k̂_guess 来动态地传入新的 k_max
            default_kMax = 40; % 直接设为合理的常数
            default_kppsMax = 20; % 直接设为合理的常数

            % 添加可选参数
            addParameter(p, 'k_max', default_kMax, @isnumeric);
            addParameter(p, 'kpps_max', default_kppsMax, @isnumeric);
            parse(p, varargin{:});

            % 获取最终的网格上限值
            final_k_max = p.Results.k_max;
            final_kpps_max = p.Results.kpps_max;

            % 设置网格参数
            cS.kMin = 0;
            cS.kMax = final_k_max;
            cS.kppsMin = 0;
            cS.kppsMax = final_kpps_max;

            % 使用非线性网格以更好地捕捉财富分布
            power_k = 2;

            % 生成常规资产网格
            if cS.nk > 1
                kGridV_temp = cS.kMin + (cS.kMax - cS.kMin) * (linspace(0, 1, cS.nk).^power_k);
            else
                kGridV_temp = cS.kMin;
            end
            cS.kGridV = kGridV_temp(:);

            % 生成PPS资产网格
            if cS.nkpps > 1
                kppsGridV_temp = cS.kppsMin + (cS.kppsMax - cS.kppsMin) * (linspace(0, 1, cS.nkpps).^power_k);
            else
                kppsGridV_temp = cS.kppsMin;
            end
            cS.kppsGridV = kppsGridV_temp(:);
        end

% =========================================================================
% ==               最终修正版的 EarningProcess_AgeDependent
% =========================================================================
function [leGridV_expanded, TrProbM_by_age, leProb1V, nw_expanded] = EarningProcess_AgeDependent(cS)
    % =========================================================================
    % == 函数: EarningProcess_AgeDependent
    % == 版本: [v12 - 最终转置修正版]
    % ==
    % == 目的:
    % ==   纠正之前对 tauchen 函数输出的错误假设。
    % ==
    % == 核心修正:
    % ==   - [!!! 关键 !!!] 移除了对 tauchen 输出矩阵的转置操作。
    % ==     经过验证，tauchen 函数本身就生成了一个标准的“行随机”矩阵，
    % ==     额外的转置是错误的，并导致了行和不为1的警告。
    % ==     此修正确保了从源头开始就使用正确的转移概率矩阵。
    % =========================================================================

    % --- 1. 生成基础AR(1)过程 ---
    lePersistence = 0.9; leShockStd = 0.15^0.5; Tauchen_q = 2.0;
    [leLogGridV_normal, leTrProbM_from_tauchen] = model_setup_utils_bgp.tauchen(cS.nw, lePersistence, leShockStd, 0, Tauchen_q);

    % [!!! 根本性修正 !!!]
    % tauchen函数返回的是一个标准的行随机矩阵 (行和为1)。
    % 无需转置。
    leTrProbM_normal = leTrProbM_from_tauchen;

    % 检查修正是否有效（可选但推荐）
    row_sums = sum(leTrProbM_normal, 2);
    if any(abs(row_sums - 1.0) > 1e-6)
        % 如果这里仍然报警，说明 tauchen 函数的内部实现有问题。
        warning('EarningProcess_AgeDependent:RowSumCheck', ...
                '转移矩阵 leTrProbM_normal 的行和不为1，请检查 tauchen 函数的内部实现！');
    end
    
    % 现在基于正确的、行随机的 leTrProbM_normal 计算平稳分布
    % 注意：为了找到对应于特征值1的左特征向量（即平稳分布），我们需要转置矩阵
    [V, D] = eig(leTrProbM_normal'); 
    [~, c] = min(abs(diag(D)-1));
    % 确保 V(:,c) 是平稳分布向量
    leProb1V = abs(V(:,c)/sum(V(:,c))); % 新生儿的长期冲击分布

    % --- 2. [自动切换] 检测是否存在重大支出冲击 ---
    has_major_shocks = (isfield(cS, 'p_shock_young_peak_annual') && cS.p_shock_young_peak_annual > 0) || ...
        (isfield(cS, 'p_shock_old_peak_annual') && cS.p_shock_old_peak_annual > 0);

    if has_major_shocks
        % --- 路径 A: 生成包含重大支出冲击的扩展过程 ---
        % (此部分代码无需修改，因为它现在会接收到正确的行随机 leTrProbM_normal)
        fprintf('   正在生成【含】重大支出冲击的收入过程...\n');
        phys_ages = cS.age1_orig : cS.ageLast_orig;
        num_phys_ages = length(phys_ages);
        p_young_annual = zeros(num_phys_ages, 1);
        p_old_annual = zeros(num_phys_ages, 1);
        y_start_idx = find(phys_ages == cS.young_shock_start_age, 1);
        y_end_idx = find(phys_ages == cS.young_shock_end_age, 1);
        if ~isempty(y_start_idx) && ~isempty(y_end_idx)
            p_young_annual(y_start_idx:y_end_idx) = cS.p_shock_young_peak_annual;
        end
        o_start_idx = find(phys_ages == cS.old_shock_start_age, 1);
        o_end_idx = find(phys_ages == cS.old_shock_end_age, 1);
        if ~isempty(o_start_idx) && ~isempty(o_end_idx)
            p_old_annual(o_start_idx:o_end_idx) = cS.p_shock_old_peak_annual;
        end
        p_shock_path_model = zeros(cS.aD_new, 2);
        for a = 1:cS.aD_new
            phys_age_indices_in_group = cS.physAgeMap{a};
            p_no_shock_y = 1.0;
            p_no_shock_o = 1.0;
            for phys_idx = phys_age_indices_in_group
                if phys_idx > 0 && phys_idx <= num_phys_ages
                    p_no_shock_y = p_no_shock_y * (1 - p_young_annual(phys_idx));
                    p_no_shock_o = p_no_shock_o * (1 - p_old_annual(phys_idx));
                end
            end
            p_shock_path_model(a, 1) = 1 - p_no_shock_y;
            p_shock_path_model(a, 2) = 1 - p_no_shock_o;
        end
        leGridV_expanded = [exp(leLogGridV_normal(:)); 0; 0];
        nw_expanded = cS.nw + 2;
        TrProbM_by_age = cell(cS.aD_new, 1);
        for a = 1:cS.aD_new
            Tr_a = zeros(nw_expanded, nw_expanded);
            p_y = p_shock_path_model(a, 1);
            p_o = p_shock_path_model(a, 2);
            p_total_shock = min(1.0, p_y + p_o);
            p_y_norm = 0;
            p_o_norm = 0;
            if p_total_shock > 0
                p_y_norm = p_y / p_total_shock;
                p_o_norm = p_o / p_total_shock;
            end
            Tr_a(1:cS.nw, 1:cS.nw) = leTrProbM_normal * (1 - p_total_shock);
            Tr_a(1:cS.nw, cS.nw + 1) = p_total_shock * p_y_norm;
            Tr_a(1:cS.nw, cS.nw + 2) = p_total_shock * p_o_norm;
            Tr_a(cS.nw + 1, 1:cS.nw) = leProb1V';
            Tr_a(cS.nw + 2, 1:cS.nw) = leProb1V';
            row_sums_a = sum(Tr_a, 2);
            row_sums_a(row_sums_a == 0) = 1; % 防止除以零
            Tr_a = Tr_a ./ row_sums_a;
            TrProbM_by_age{a} = Tr_a;
        end
    else
        % --- 路径 B: 生成不含重大支出冲击的标准过程 ---
        fprintf('   检测到重大冲击参数为零，正在生成【无】重大支出冲击的标准收入过程...\n');
        nw_expanded = cS.nw;
        leGridV_expanded = exp(leLogGridV_normal(:));
        TrProbM_by_age = cell(cS.aD_new, 1);
        for a = 1:cS.aD_new
            % 现在 leTrProbM_normal 已经是正确的行随机矩阵
            TrProbM_by_age{a} = leTrProbM_normal;
        end
    end
end        
function [y_grid_out, trProbM_out] = tauchen(N, rho, sigma, mu, m)
                                    % --- [新增] 处理 N=1 (同质代理人) 的边缘情况 ---
            if N == 1
                y_grid_out = 0;      % 唯一的网格点是对数均值0
                trProbM_out = 1;     % 转移概率为100%
                return;              % 直接返回，跳过后续计算
            end
            
            % [BGP不变] 标准的、用于离散化AR(1)过程的工具函数
            std_y = sqrt(sigma^2 / (1-rho^2)); y_max = m*std_y; y_min = -y_max;
            y = linspace(y_min, y_max, N); d = y(2)-y(1);
            trProbM_out = zeros(N,N);
            for j=1:N, for k=1:N
                    m_k = rho*y(j) + mu;
                    if k==1, trProbM_out(j,k) = normcdf((y(1)-m_k+d/2)/sigma);
                    elseif k==N, trProbM_out(j,k) = 1 - normcdf((y(N)-m_k-d/2)/sigma);
                    else, trProbM_out(j,k) = normcdf((y(k)-m_k+d/2)/sigma) - normcdf((y(k)-m_k-d/2)/sigma); end
            end, end
        y_grid_out = y(:);
        end

        function [muM, utilM] = CES_utility(cM, sigma, cS)
            % [BGP修改] 此函数的作用是输入一个消费值，输出其效用
            % [BGP修改] 在标准化框架下，VFI求解器会向它传递标准化的消费 ĉ
            % [BGP修改] 这个函数会正确地计算 u(ĉ) = ĉ^(1-σ)/(1-σ)
            % [BGP修改] 无需知道绝对技术水平A，完全在标准化环境中操作
            if cM<=0
                c_adj = cS.cFloor;
                print('消费小于零')
            else
                c_adj = cM;
            end
            if abs(sigma - 1) < 1e-6
                utilM = log(c_adj);
                muM = 1./c_adj;
            else utilM = (c_adj.^(1-sigma))./(1-sigma);
                muM = c_adj.^(-sigma);
            end
        end

        function util_beq = bequest_utility(k_prime, cS)
            % [BGP修改] 此函数将被传递标准化的遗赠 k̂'，并正确计算其效用
            % [BGP修改] 在标准化框架下，遗赠效用函数处理标准化的遗赠值
            % [BGP修改] 无需知道绝对技术水平A，完全在标准化环境中操作
            if cS.phi_bequest <= 0
                util_beq = 0;
                return;
            end

            k_adj = max(cS.cFloor, k_prime);

            if abs(cS.sigma - 1) < 1e-6
                util_beq = cS.phi_bequest * log(k_adj);
            else
                util_beq = cS.phi_bequest * (k_adj.^(1-cS.sigma))./(1-cS.sigma);
            end
        end


        function A_path = generate_tfp_path(cS, annual_years_vec, model_sim_years)
            % =========================================================================
            % == 函数: generate_tfp_path (v_curvature_control - 曲率可控最终版)
            % ==
            % == 核心修正:
            % ==   调用新的smooth_transition函数，并传入一个曲率参数，
            % ==   确保TFP增长率路径是平滑收敛的。
            % =========================================================================

            fprintf('   正在构建TFP路径 (v_curvature_control - 曲率可控最终版)...\n');

            % --- 1. 加载数据和设定参数 ---
            tfp_data_pwt = readtable('data\PWT\china_pwt_data.xlsx');
            pwt_years = tfp_data_pwt.year(tfp_data_pwt.year >= 1979);
            pwt_tfp_level = tfp_data_pwt.rtfpna(tfp_data_pwt.year >= 1979);
            bai_projections = [2025, 0.0557; 2030, 0.0482; 2035, 0.0394; 2040, 0.0340; 2045, 0.0346; 2050, 0.0298];

            long_term_g_annual = cS.g_A_ss;
            fprintf('   TFP长期年化增长率将平滑收敛至: %.4f%% (曲率k=%.1f)\n', long_term_g_annual * 100, cS.curvature_param);

            data_end_year = bai_projections(end, 1);
            transition_end_year = data_end_year + cS.transition_period_years;

            % --- 2. 构建平滑的【年度增长率】路径 ---
            g_path_annual = zeros(size(annual_years_vec));
            g_anchor_years = [pwt_years(end); bai_projections(:,1)];
            g_anchor_rates = [0.058; bai_projections(:,2)];

            data_mask = annual_years_vec <= data_end_year;
            g_path_annual(data_mask) = interp1(g_anchor_years, g_anchor_rates, annual_years_vec(data_mask), 'pchip', 'extrap');
            last_data_g = g_path_annual(find(data_mask, 1, 'last'));

            transition_mask = annual_years_vec > data_end_year & annual_years_vec <= transition_end_year;
            if any(transition_mask)
                num_trans_years = sum(transition_mask);
                g_path_annual(transition_mask) = model_setup_utils_bgp.smooth_transition(last_data_g, long_term_g_annual, num_trans_years, cS.curvature_param);
            end

            ss_mask = annual_years_vec > transition_end_year;
            g_path_annual(ss_mask) = long_term_g_annual;

            % --- 3. 从增长率路径构建水平路径 ---
            A_path_annual_level = zeros(size(annual_years_vec));

            % 找到起始年份的TFP水平作为锚点
            start_year_tfp_level_anchor = interp1([pwt_years; 2100], [pwt_tfp_level; pwt_tfp_level(end)*(1.02)^50], cS.start_year, 'pchip', 'extrap');

            start_year_idx = find(annual_years_vec == cS.start_year, 1);
            A_path_annual_level(start_year_idx) = start_year_tfp_level_anchor;

            % 向前积分
            for t_idx = (start_year_idx + 1):length(annual_years_vec)
                A_path_annual_level(t_idx) = A_path_annual_level(t_idx - 1) * (1 + g_path_annual(t_idx));
            end
            % 向后积分
            for t_idx = (start_year_idx - 1):-1:1
                A_path_annual_level(t_idx) = A_path_annual_level(t_idx + 1) / (1 + g_path_annual(t_idx + 1));
            end

            % --- 4. 采样并归一化 ---
            A_path = interp1(annual_years_vec, A_path_annual_level, model_sim_years, 'linear');
            if A_path(1) ~= 0 && ~isnan(A_path(1))
                A_path = A_path / A_path(1);
            else
                warning('初始TFP水平为0或NaN，无法归一化。');
            end
            fprintf('✅ TFP路径构建完成。\n');
        end
        
        function create_plots(cS, plot_years, Z_path, A_path)
            fprintf('   正在生成外生路径可视化图表...\n');
            T_sim = size(Z_path, 2);
            figure('Name', '外生路径可视化', 'Position', [100, 100, 1200, 500]);
            subplot(1, 2, 1);
            selected_ages = [1, 4, 8, 12, 16];
            age_labels = {};
            for i = 1:length(selected_ages)
                age_start = cS.age1_orig + (selected_ages(i) - 1) * cS.time_Step;
                age_end = age_start + cS.time_Step - 1;
                age_labels{i} = sprintf('%d-%d岁', age_start, age_end);
            end
            colors = lines(length(selected_ages));
            for i = 1:length(selected_ages)
                plot(plot_years, Z_path(selected_ages(i), :) * 100, 'Color', colors(i,:), 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 4);
                hold on;
            end
            xlabel('年份'); ylabel('人口占比 (%)'); title('人口分布路径 (Z\_path)');
            legend(age_labels, 'Location', 'best', 'FontSize', 8); grid on;
            xlim([plot_years(1), plot_years(end)]);
            subplot(1, 2, 2);
            plot(plot_years, A_path, 'b-', 'LineWidth', 2.5, 'Marker', 's', 'MarkerSize', 5, 'MarkerFaceColor', 'b');
            xlabel('年份'); ylabel('TFP水平 (相对于起始年份)'); title('全要素生产率路径 (A\_path)');
            grid on; xlim([plot_years(1), plot_years(end)]);
            avg_growth_rate = (A_path(end)/A_path(1))^(1/(T_sim*cS.time_Step)) - 1;
            text(0.05, 0.95, sprintf('平均年增长率: %.2f%%', avg_growth_rate*100), 'Units', 'normalized', 'FontSize', 10, 'BackgroundColor', 'white');
            sgtitle(sprintf('外生路径概览 (%d-%d年, %d个时期)', cS.start_year, cS.end_year, T_sim), 'FontSize', 14);
            fprintf('✅ 外生路径可视化完成。\n');
        end

        function Z_ss0 = get_calibration_inputs(target_year, cS)
            % [BGP不变] 此函数的作用是为初始稳态提供校准年份的人口分布 Z_ss0
            % 这个功能和逻辑在新的框架下完全正确且必要
            fprintf('   为校准获取 %d 年的输入...\n', target_year);

            % 1. 加载人口数据
            pop_data = readtable('data\人口\population_by_age_group_all_years.xlsx', 'Sheet', 'pop_normalized');
            data_years = str2double(cellfun(@(x) x(2:end), pop_data.Properties.VariableNames(2:end), 'UniformOutput', false));

            % 2. 插值得到目标年份的人口分布
            Z_ss0_raw = interp1(data_years, pop_data{:, 2:end}', target_year, 'linear', 'extrap')';

            % 3. 归一化
            Z_ss0 = Z_ss0_raw / sum(Z_ss0_raw);
        end

        function cS = calcaulte_theta_payg_path(cS, graph_flag)
            % [BGP不变] 此函数定义了一条外生的政策路径（养老金缴费率）
            % 它不受模型内生增长机制的影响
            fprintf('正在构建基于【覆盖率插值】的有效养老金缴费率路径 (theta_path)...\n');

            % --- 1. 定义关键参数 ---
            theta_urban_employee_effective = 0.24;
            theta_resident_effective = 0.03;

            % [政策目标] 定义未来的目标覆盖率
            coverage_urban_final = 1;
            coverage_resident_final = 1;
            year_reach_final = 2050;

            % --- 2. 收集来自官方统计的真实数据点 (单位：万人) ---
            year_pax_urban = [1997, 2000, 2005, 2008, 2011, 2014, 2017, 2020, 2023];
            pax_urban      = [8671, 10448, 17487, 21890, 28391, 34124, 40293, 45638, 52112];

            year_pax_resident = [2008, 2011, 2014, 2017, 2020, 2023];
            pax_resident      = [0,    45586, 50075, 51255, 54227, 54475];

            year_laborpop_data = [1997, 2000, 2005, 2008, 2011, 2014, 2017, 2020, 2023];
            laborpop_data      = [78531, 82352, 87532, 90833, 94072, 92837, 91199, 89552, 86481];

            % --- 3. 计算历史数据点的【覆盖率】 ---
            coverage_urban_data = pax_urban ./ laborpop_data;

            laborpop_for_resident_years = interp1(year_laborpop_data, laborpop_data, year_pax_resident, 'linear', 'extrap');
            coverage_resident_data = pax_resident ./ laborpop_for_resident_years;

            % --- 4. 直接对【覆盖率】进行插值，生成完整的年度路径 ---
            annual_years_vector = cS.start_year:cS.end_year;

            % a. 插值生成年度【城镇职工覆盖率】路径
            interp_years_urban = [year_pax_urban, year_reach_final, cS.end_year];
            interp_coverage_urban = [coverage_urban_data, coverage_urban_final, coverage_urban_final];
            [unique_years_u, ia_u, ~] = unique(interp_years_urban);
            coverage_urban_annual = interp1(unique_years_u, interp_coverage_urban(ia_u), annual_years_vector, 'linear');

            % b. 插值生成年度【城乡居民覆盖率】路径
            interp_years_resident = [cS.start_year, year_pax_resident, year_reach_final, cS.end_year];
            interp_coverage_resident = [0, coverage_resident_data, coverage_resident_final, coverage_resident_final];
            [unique_years_r, ia_r, ~] = unique(interp_years_resident);
            coverage_resident_annual = interp1(unique_years_r, interp_coverage_resident(ia_r), annual_years_vector, 'linear');
            coverage_resident_annual(annual_years_vector < min(year_pax_resident)) = 0;

            % --- 5. 基于【覆盖率路径】计算最终的【年度有效缴费率】路径 ---
            theta_path_annual = (coverage_urban_annual * theta_urban_employee_effective) + ...
                (coverage_resident_annual * theta_resident_effective);

            % --- 6. 从年度路径中，提取出模型【5年期】所需的路径 ---
            model_year_indices = 1:(cS.end_year - cS.start_year + 1);
            model_year_indices_5yr = model_year_indices(1:cS.time_Step:end);

            T_sim = length(cS.start_year:cS.time_Step:cS.end_year);
            theta_path = theta_path_annual(model_year_indices_5yr);
            theta_path = theta_path(1:min(T_sim, length(theta_path)));

            % --- 7. 将最终路径存入 cS 结构体 ---
            cS.theta_path = theta_path;

            % --- 8. (可选) 可视化检查路径 ---
            if graph_flag
                T_plot = length(cS.theta_path);
                time_axis = cS.start_year:cS.time_Step:(cS.start_year + cS.time_Step*(T_plot-1));

                figure('Name', 'Effective PAYG Tax Rate Path (theta_path) - Corrected');
                plot(time_axis, cS.theta_path, 'k-s', 'LineWidth', 2, 'MarkerFaceColor', 'k', 'DisplayName', '总有效缴费率 (加权平均)');
                hold on;

                contribution_urban = coverage_urban_annual(model_year_indices_5yr(1:T_plot)) * theta_urban_employee_effective;
                contribution_resident = coverage_resident_annual(model_year_indices_5yr(1:T_plot)) * theta_resident_effective;

                plot(time_axis, contribution_urban, 'b--o', 'LineWidth', 1.5, 'DisplayName', '城镇职工体系贡献');
                plot(time_axis, contribution_resident, 'r--d', 'LineWidth', 1.5, 'DisplayName', '城乡居民体系贡献');

                title('模型使用的有效养老金缴费率路径 (基于覆盖率插值)');
                xlabel('年份');
                ylabel('有效缴费率 (θ_t)');
                legend('show', 'Location', 'best');
                grid on;
                ylim([0, max(cS.theta_path)*1.2]);
            end
            fprintf('✅ 基于覆盖率插值的缴费率路径构建完成。\n');
        end

        function plot_population_pyramid_simplified(Z_out)
            % 简化版人口金字塔可视化

            % 年龄组标签
            age_groups = arrayfun(@(a) sprintf('%d-%d', 20+(a-1)*5, 24+(a-1)*5), 1:16, 'UniformOutput', false);

            % 选择三个时间点
            t_initial = 1;
            t_peak_aging = 10; % 约50年后
            t_steady_state = 40;

            figure('Name', '人口结构金字塔演变', 'Position', [100, 100, 1400, 500]);

            % --- 子图1: 初始年份 ---
            subplot(1, 3, 1);
            barh(Z_out(:, t_initial) * 100, 'FaceColor', '#0072BD');
            set(gca, 'YDir', 'reverse');
            yticks(1:16);
            yticklabels(age_groups);
            xlabel('人口占比 (%)');
            title(sprintf('初始年 (t=%d)', t_initial));
            grid on;
            xlim([0, 12]);

            % --- 子图2: 老龄化高峰期 ---
            subplot(1, 3, 2);
            barh(Z_out(:, t_peak_aging) * 100, 'FaceColor', '#D95319');
            set(gca, 'YDir', 'reverse');
            yticks(1:16);
            yticklabels(age_groups);
            title(sprintf('老龄化高峰附近 (t=%d)', t_peak_aging));
            grid on;
            xlim([0, 12]);

            % --- 子图3: 长期稳态 ---
            subplot(1, 3, 3);
            barh(Z_out(:, t_steady_state) * 100, 'FaceColor', '#77AC30');
            set(gca, 'YDir', 'reverse');
            yticks(1:16);
            yticklabels(age_groups);
            title(sprintf('长期稳态 (t=%d)', t_steady_state));
            grid on;
            xlim([0, 12]);

            sgtitle('人口结构从“中年型”向“老年型”的动态演变', 'FontSize', 16, 'FontWeight', 'bold');
        end

        function Z_theory = compute_theoretical_ss_dist(cS, n_ss_annual, s_pathV)
            % =========================================================================
            % == 函数: compute_theoretical_ss_dist (原 compute_theoretical_ss_dist_zpg)
            % == 版本: [v2 - 人口增长兼容版]
            % ==
            % == 目的: 计算一个具有恒定人口年增长率n的理论稳态人口分布。
            % ==
            % == 核心修改:
            % ==   1. 增加一个输入参数 n_ss_annual (年化稳态人口增长率)。
            % ==   2. 在计算下一代人口的相对规模时，不仅考虑存活率 s_pathV，
            % ==      还用人口增长因子 (1+n_period) 进行折现。这反映了
            % ==      一个增长中的人口，年轻一代的规模相对更大。
            % ==   3. 当 n_ss_annual = 0 时，此函数的结果与原ZPG版本完全一致。
            % =========================================================================
            fprintf('   正在计算具有年增长率 n=%.4f 的理论稳态分布...\n', n_ss_annual);

            % 1. 将年化增长率转换为模型期的增长率
            n_period = (1 + n_ss_annual)^cS.time_Step - 1;

            % 2. 以新生儿为基准1.0，计算各年龄组的相对人口规模
            mass_levels_by_age = zeros(cS.aD_new, 1);
            mass_levels_by_age(1) = 1.0;
            for ia = 1:(cS.aD_new - 1)
                % 核心逻辑：下一代的人口规模由上一代的存活者构成，但相对于
                % 新一代的新生儿，他们的相对规模因人口整体增长而被“稀释”。
                % 因此，需要除以人口增长因子 (1 + n_period)。
                mass_levels_by_age(ia+1) = (mass_levels_by_age(ia) * s_pathV(ia)) / (1 + n_period);
            end

            % 3. 将相对规模归一化，得到最终的概率分布
            Z_theory = mass_levels_by_age / sum(mass_levels_by_age);

            fprintf('   ✅ 理论稳态分布计算完成。\n');
        end

        function pol_age_expanded = expand_policy_slice(pol_age_slice, nkpps)
            pol_age_expanded = pol_age_slice;
            fields = fieldnames(pol_age_expanded);
            for i = 1:length(fields)
                field_name = fields{i};
                pol_age_expanded.(field_name) = repmat(pol_age_slice.(field_name), [1, nkpps, 1]);
            end
        end

        function create_exogenous_dynamics_plots(cS, Z_out, A_path_annual, birth_rate_path_annual, annual_years_vec, pop_growth_rate_path_annual)
            % =========================================================================
            % == 函数: create_exogenous_dynamics_plots
            % == 版本: [v_viz_growth_rate - 可视化增长率最终版]
            % ==
            % == 核心功能:
            % ==   1. 接收所有关键外生路径（人口、TFP、出生率、人口增长率）。
            % ==   2. 使用双Y轴绘制功能强大的“人口动态”和“技术动态”仪表盘。
            % ==   3. 智能处理年度和模型期数据，确保曲线对齐。
            % =========================================================================

            fprintf('   正在生成外生动态可视化图表 (v_viz_growth_rate)...\n');

            % --- 创建图形窗口 ---
            figure('Name', '外生路径核心动态可视化', 'Position', [100, 100, 1400, 600]);

            % --- 图 1: 人口动态 (老龄化、增长率 vs 出生率) ---
            subplot(1, 2, 1);

            % -- 左轴: 结构与增长 --
            yyaxis left

            % 绘制退休人口占比 (基于模型期数据)
            retiree_indices = (cS.aR_new + 1):cS.aD_new;
            model_sim_years = cS.start_year:cS.time_Step:cS.end_year;
            retiree_proportion_path = sum(Z_out(retiree_indices, :), 1);
            p1 = plot(model_sim_years, retiree_proportion_path * 100, 'b-o', 'LineWidth', 2.5, 'MarkerSize', 4, 'MarkerFaceColor', 'b');
            hold on;

            % 绘制人口年增长率 (增长率对应的是时期之间的变化)
            years_for_growth_rate = model_sim_years(2:end) - cS.time_Step/2; % 将点绘制在时期中点
            if ~isempty(pop_growth_rate_path_annual)
                p2 = plot(years_for_growth_rate, pop_growth_rate_path_annual * 100, 'm-s', 'LineWidth', 2, 'MarkerSize', 4, 'MarkerFaceColor','m');
            end
            hold off;

            ylabel('占比或增长率 (%)');
            ax = gca;
            ax.YAxis(1).Color = 'k';

            % -- 右轴: 驱动因素 --
            yyaxis right

            % 绘制年化出生率 (年度数据)
            p3 = plot(annual_years_vec, birth_rate_path_annual * 1000, 'r--', 'LineWidth', 2);
            ylabel('年化出生率 (每千人)', 'Color', 'r');
            ax.YAxis(2).Color = 'r';

            % -- 图形格式 --
            title('人口动态：结构、增长与出生率');
            xlabel('年份');
            grid on;
            xlim([cS.start_year, cS.end_year]);
            legend_handles = [p1, p3];
            legend_texts = {'退休人口占比 (%) [左轴]', '年化出生率 (每千人) [右轴]'};
            if exist('p2', 'var'), legend_handles = [legend_handles, p2]; legend_texts{end+1} = '人口年增长率 (%) [左轴]'; end
            legend(legend_handles, legend_texts, 'Location', 'best');

            % --- 图 2: 技术动态 (TFP水平 vs 增长率) ---
            subplot(1, 2, 2);

            if ~isempty(A_path_annual)
                % 计算TFP年增长率 (年度数据)
                tfp_growth_rate_annual = (A_path_annual(2:end) ./ A_path_annual(1:end-1)) - 1;

                % 左轴: 绘制TFP水平 (归一化，年度数据)
                A_path_annual_normalized = A_path_annual / A_path_annual(1);
                yyaxis left
                plot(annual_years_vec, A_path_annual_normalized, 'k-s', 'LineWidth', 2.5, 'MarkerSize', 4, 'MarkerFaceColor', 'k');
                ylabel('TFP 水平 (相对于起始年份)', 'Color', 'k');
                ax_tech = gca;
                ax_tech.YAxis(1).Color = 'k';

                % 右轴: 绘制TFP年增长率 (年度数据)
                yyaxis right
                plot(annual_years_vec(2:end), tfp_growth_rate_annual * 100, 'g--', 'LineWidth', 2);
                hold on;
                yline(cS.g_A_ss * 100, 'g:', 'LineWidth', 1.5);
                hold off;
                ylabel('TFP 年增长率 (%)', 'Color', 'g');
                ax_tech.YAxis(2).Color = 'g';

                legend('TFP水平', 'TFP年增长率', sprintf('稳态增长率 (%.2f%%)', cS.g_A_ss*100), 'Location', 'best');
            else
                text(0.5, 0.5, 'TFP路径无法生成 (T_{sim} < 2)', 'HorizontalAlignment', 'center');
            end

            % 图形格式
            title('技术动态：TFP水平与增长率');
            xlabel('年份');
            grid on;
            xlim([cS.start_year, cS.end_year]);

            sgtitle(sprintf('外生路径概览 (%d-%d年)', cS.start_year, cS.end_year), 'FontSize', 16, 'FontWeight', 'bold');
            fprintf('✅ 外生动态可视化完成。\n');
        end

        function smooth_path = smooth_transition(start_val, end_val, num_steps, curvature)
            % =========================================================================
            % == 辅助函数: smooth_transition (v2.0 - 曲率可控版)
            % == 目的: 使用幂函数生成一个从start_val到end_val的、曲率可控的过渡路径。
            % ==
            % == 输入:
            % ==   start_val: 过渡路径的起始值
            % ==   end_val:   过渡路径的终点值
            % ==   num_steps: 过渡期的长度 (年份数)
            % ==   curvature: 曲率参数 'k'
            % ==              k > 1: 下凸式收敛 (二阶导>0, 变化率先慢后快)
            % ==              k = 1: 线性收敛 (等同于linspace)
            % ==              0 < k < 1: 上凸式收敛 (二阶导<0, 变化率先快后慢)
            % =========================================================================
            if nargin < 4
                curvature = 2; % 默认使用下凸式收敛
            end
            if num_steps <= 0
                smooth_path = [];
                return;
            end

            t = linspace(0, 1, num_steps);
            weights = t .^ curvature;
            smooth_path = (1 - weights) * start_val + weights * end_val;
        end

        function cS = generate_mortality_path(cS)
            % =========================================================================
            % == 函数: generate_mortality_path
            % == 版本: [v2.1 - 变量名修正版]
            % ==
            % == 目的:
            % ==   1. 从外部文件加载2023-2100年的年度死亡率数据。
            % ==   2. 计算并生成【时变】的存活率路径。
            % ==   3. 输出的 s_pathV 和 s_path_all 将是矩阵，其中每一列代表一个年份的存活率向量。
            % ==
            % == 核心修正:
            % ==   - 修正了由于变量名不匹配导致的“无法识别的函数或变量”错误。
            % ==   - 将数据加载和路径计算整合到一个函数中，并理顺了逻辑。
            % =========================================================================
            
            fprintf('--- 启动时变存活率路径生成器 (v2.1) ---\n');

            % --- 1. 从外部文件加载时变死亡率数据 ---
            filePath = 'data/人口/WPP_MORT_CHINA_2023_2100_matrix.xlsx';
            
            fprintf('   正在从 "%s" 加载年度死亡率数据...\n', filePath);
            try
                dataTable = readtable(filePath); % 使用readtable读取整个Excel文件
            catch ME
                error('无法加载死亡率数据文件: %s\n请检查文件路径和格式。错误详情: %s', filePath, ME.message); % 如果失败，抛出错误
            end
            
            % --- 2. 提取数据并验证 ---
            if ~ismember('Year', dataTable.Properties.VariableNames)
                error('死亡率数据文件中未找到 "Year" 列。'); % 验证年份列是否存在
            end
            yearsV = dataTable.Year; % 提取年份向量
            mortalityM = table2array(dataTable(:, 2:end)); % 提取死亡率数据矩阵 (从第二列到最后一列)
            
            if size(mortalityM, 2) ~= 101 % 验证数据维度是否为101 (0-100岁)
                warning('死亡率数据列数不为101。检测到 %d 列，请检查文件。', size(mortalityM, 2));
            end
            
            % --- 3. 将加载的原始数据存入cS结构体 ---
            % [BUG修复] 使用正确的变量名(yearsV, mortalityM)进行赋值
            cS.mortality_years = yearsV;
            cS.mortality_path_data = mortalityM;
            num_years = length(cS.mortality_years);
            fprintf('   ✅ 成功加载 %d 个年份 (%d-%d) 的死亡率数据。\n', num_years, min(yearsV), max(yearsV));

            % --- 4. 预计算与年龄相关的静态参数 ---
            cS.aD_orig = cS.ageLast_orig - cS.age1_orig + 1;
            cS.aR_idx_orig = cS.ageRetire_orig - cS.age1_orig + 1;
            cS.aW_orig = cS.aR_idx_orig - 1;
            cS.physAgeV_orig = (cS.age1_orig : cS.ageLast_orig)';

            cS.aD_new = ceil(cS.aD_orig / cS.time_Step);
            cS.aR_new = ceil((cS.ageRetire_orig - cS.age1_orig) / cS.time_Step);
            cS.physAgeMap = cell(cS.aD_new, 1);
            for a = 1:cS.aD_new
                startIdx = (a-1)*cS.time_Step + 1;
                endIdx = min(a*cS.time_Step, cS.aD_orig);
                cS.physAgeMap{a} = cS.age1_orig + (startIdx:endIdx) - 1;
            end
            num_age_groups_all = ceil((100+1) / cS.time_Step);

            % --- 5. 预分配输出矩阵 ---
            s_pathV_over_time = zeros(cS.aD_new, num_years);
            s_path_all_over_time = zeros(num_age_groups_all, num_years);
            
            fprintf('   正在为 %d 个年份计算模型期存活率路径...\n', num_years);

            % --- 6. 核心循环：为每一年计算存活率路径 ---
            for t = 1:num_years
                % a. 提取当前年份t的年度死亡率数据 (qx)
                raw_mortality_data_t = cS.mortality_path_data(t, :);
                raw_mortality_data_t(end) = 1.0; % 强制确保100岁确定性死亡

                % b. 计算当前年份t的年度存活率
                raw_survival_data_t = 1 - raw_mortality_data_t;
                
                % c. 计算【模型核心存活率 s_pathV】(用于VFI)
                s_pathV_t = zeros(cS.aD_new, 1);
                for a = 1:(cS.aD_new - 1)
                    s_period_prod = 1.0;
                    phys_ages_in_group = cS.physAgeMap{a};
                    for age_idx = 1:length(phys_ages_in_group)
                        age = phys_ages_in_group(age_idx);
                        % 年龄从0开始索引, 所以物理年龄age对应于数据中的第age+1个元素
                        if (age + 1) <= length(raw_survival_data_t)
                            s_period_prod = s_period_prod * raw_survival_data_t(age + 1);
                        else
                            s_period_prod = 0; % 超出数据范围则无法存活
                            break;
                        end
                    end
                    s_pathV_t(a) = s_period_prod;
                end
                s_pathV_t(cS.aD_new) = 0; % 最后一期确定性死亡
                
                % d. 计算【全年龄组存活率 s_path_all】(用于人口模拟)
                s_path_all_t = zeros(num_age_groups_all, 1);
                for a = 1:(num_age_groups_all - 1)
                    s_period_prod = 1.0;
                    start_age = (a-1) * cS.time_Step;
                    for i = 1:cS.time_Step
                        age = start_age + i - 1;
                        if (age + 1) <= length(raw_survival_data_t)
                            s_period_prod = s_period_prod * raw_survival_data_t(age + 1);
                        else
                            s_period_prod = 0;
                            break;
                        end
                    end
                    s_path_all_t(a) = s_period_prod;
                end
                s_path_all_t(end) = 0; % 最后一个年龄组确定性死亡

                % e. 将当年计算出的向量存入输出矩阵的对应列
                s_pathV_over_time(:, t) = s_pathV_t;
                s_path_all_over_time(:, t) = s_path_all_t;
            end

            % --- 7. 将时变路径矩阵赋值给cS结构体 ---
            cS.s_pathV = s_pathV_over_time;
            cS.s_path_all = s_path_all_over_time;
            
            fprintf('   ✅ 成功生成时变存活率路径矩阵。\n');
            fprintf('      - cS.s_pathV 维度: [%d, %d]\n', size(cS.s_pathV, 1), size(cS.s_pathV, 2));
            fprintf('      - cS.s_path_all 维度: [%d, %d]\n', size(cS.s_path_all, 1), size(cS.s_path_all, 2));
        end        
    end
end