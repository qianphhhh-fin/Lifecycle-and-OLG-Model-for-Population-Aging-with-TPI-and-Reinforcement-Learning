classdef population
    methods (Static)

        function [Z_out,Z_out_raw, cS] = generate_Z_path(cS, graph_flag)
            % =========================================================================
            % == 函数: generate_Z_path
            % == 版本: [v3.1 - 人口路径对齐版]
            % ==
            % == 核心修改:
            % ==   - 使用统一的 cS.pop_data_last_year 参数来控制出生率和死亡率
            % ==     路径的外推行为，确保两者同步进入稳态。
            % ==   - 在计算每期存活率时，增加了对 cS.pop_data_last_year 的判断，
            % ==     确保在预测期后使用恒定的死亡率数据。
            % =========================================================================

            fprintf('\n--- 启动外生路径生成器 (v3.1 - 人口路径对齐版) ---\n');

            % --- 1. 初始化模拟参数 ---
            ss_convergence_tol = 1e-6;
            max_sim_periods = 200;

            % --- 2. 加载并处理【出生率】数据 ---
            try
                cbr_data = readtable('..\data\人口\UN_PPP2024_CBR_birthper1000_China.xlsx');
            catch ME
                error('无法加载出生率数据文件: %s', ME.message);
            end
            var_names_cbr = cbr_data.Properties.VariableNames;
            year_cols_indices_cbr = find(startsWith(var_names_cbr, 'y'));
            un_years = str2double(cellfun(@(x) x(2:end), var_names_cbr(year_cols_indices_cbr), 'UniformOutput', false));
            un_cbr = table2array(cbr_data(1, year_cols_indices_cbr)) / 1000;

            % --- 3. 构建【年度出生率】路径并确定恒定规则 ---
            % [!!! 核心修改 !!!] 使用统一的截断年份
            last_data_year = cS.pop_data_last_year; 
            fprintf('   [对齐设定] 所有人口相关路径将在 %d 年后保持恒定。\n', last_data_year);

            annual_years_full_range = cS.start_year:(cS.start_year + (max_sim_periods * cS.time_Step) -1);
            birth_rate_path_annual_full = zeros(size(annual_years_full_range));
            data_mask = annual_years_full_range <= last_data_year;
            birth_rate_path_annual_full(data_mask) = interp1(un_years, un_cbr, annual_years_full_range(data_mask), 'pchip', 'extrap');
            last_data_rate = birth_rate_path_annual_full(find(data_mask, 1, 'last'));
            ss_mask = annual_years_full_range > last_data_year;
            birth_rate_path_annual_full(ss_mask) = last_data_rate;

            % --- 4. 初始化人口分布 ---
            try
                pop_data_by_age = readtable('..\data\人口\china_population_by_age_headerFix.xlsx');
            catch ME
                error('无法加载分年龄人口数据文件: %s', ME.message);
            end
            initial_pop_row = pop_data_by_age(pop_data_by_age.Year == cS.ss0_year, :);
            if isempty(initial_pop_row), error('在人口数据中找不到初始年份: %d', cS.ss0_year); end
            num_age_groups_all = size(cS.s_path_all, 1);
            initial_pop_aggregated = population.aggregate_population(initial_pop_row, num_age_groups_all, cS.time_Step);

            % --- 5. 执行人口动态前向模拟 ---
            fprintf('   正在启动人口动态模拟以确定 T_sim ...\n');
            Z_path_raw = zeros(num_age_groups_all, max_sim_periods);
            Z_path_raw(:, 1) = initial_pop_aggregated;
            Z_norm_prev = Z_path_raw(:, 1) / sum(Z_path_raw(:, 1));

            s_path_all_sim = zeros(num_age_groups_all, max_sim_periods);
            s_pathV_sim = zeros(cS.aD_new, max_sim_periods);

            converged = false;
            T_final = max_sim_periods;
            
            % 获取死亡率数据中最后一年的数据，用于外推
            mortality_data_ss_year = cS.mortality_path_data(cS.mortality_years == last_data_year, :);
            if isempty(mortality_data_ss_year)
                mortality_data_ss_year = cS.mortality_path_data(end, :); % 如果找不到精确年份，用最后一年
                warning('未在死亡率数据中找到指定的 pop_data_last_year (%d)，将使用最后一年 (%d) 的数据进行外推。', last_data_year, cS.mortality_years(end));
            end


            for t = 1:max_sim_periods
                current_year = cS.start_year + (t - 1) * cS.time_Step;

                % [!!! 核心修改 !!!] 根据截断年份决定使用哪套死亡率数据
                if current_year <= last_data_year
                    % 在预测期内，使用插值的时变死亡率
                    raw_mortality_data_t = interp1(cS.mortality_years, cS.mortality_path_data, current_year, 'linear', 'extrap');
                else
                    % 在预测期外，使用恒定的 `last_data_year` 的死亡率
                    raw_mortality_data_t = mortality_data_ss_year;
                end
                
                raw_survival_data_t = 1 - raw_mortality_data_t;
                raw_survival_data_t(end) = 0; % 确保100岁确定死亡

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

                if t < max_sim_periods
                    survivors_tp1_abs = Z_path_raw(1:end-1, t) .* s_path_all_sim(1:end-1, t);
                    total_survivors_tp1 = sum(survivors_tp1_abs);

                    period_start_year_tp1 = cS.start_year + t * cS.time_Step;
                    period_end_year_tp1 = period_start_year_tp1 + cS.time_Step - 1;

                    annual_indices = (annual_years_full_range >= period_start_year_tp1) & (annual_years_full_range <= period_end_year_tp1);
                    if ~any(annual_indices)
                        cbr_target_period = birth_rate_path_annual_full(end);
                    else
                        cbr_target_period = mean(birth_rate_path_annual_full(annual_indices));
                    end

                    cbr_target_model_period = cbr_target_period * cS.time_Step;

                    total_pop_tp1 = total_survivors_tp1 / (1 - cbr_target_model_period);
                    new_entrants_tp1 = total_pop_tp1 - total_survivors_tp1;

                    Z_path_raw(1, t+1) = new_entrants_tp1;
                    Z_path_raw(2:end, t+1) = survivors_tp1_abs;

                    Z_norm_current = Z_path_raw(:, t+1) / sum(Z_path_raw(:, t+1));
                    diff = max(abs(Z_norm_current - Z_norm_prev));
                    if ~converged && current_year > last_data_year && diff < ss_convergence_tol
                        T_final = t + 20; converged = true;
                        fprintf('   ✅ 完整人口结构在第 %d 期 (年份: %d) 达到稳态, 额外增加至 %d 期。\n', t, cS.start_year + t*cS.time_Step, T_final);
                    end
                    Z_norm_prev = Z_norm_current;
                end
            end
            if ~converged, T_final = max_sim_periods; warning('人口模拟在达到最大期数 %d 后仍未收敛！', max_sim_periods); end

            % --- 6. 裁剪所有路径，使其长度等于T_final ---
            cS.T_sim = T_final;
            cS.end_year = cS.start_year + (T_final - 1) * cS.time_Step;
            fprintf('   内生决定的模拟期数 T_sim = %d, 结束年份 = %d\n', cS.T_sim, cS.end_year);

            Z_path_raw_final = Z_path_raw(:, 1:cS.T_sim);
            cS.Z_path_raw_final = Z_path_raw_final;
            model_age_start_index = floor(cS.age1_orig / cS.time_Step) + 1;
            model_age_end_index = model_age_start_index + cS.aD_new - 1;
            Z_out_raw = Z_path_raw_final(model_age_start_index : model_age_end_index, :);
            Z_out = Z_out_raw ./ sum(Z_out_raw, 1);
            if size(Z_out, 1) ~= cS.aD_new, warning('Z_out的行数(%d)与cS.aD_new(%d)不匹配！', size(Z_out, 1), cS.aD_new); end

            cS.s_path_all = s_path_all_sim(:, 1:cS.T_sim);
            cS.s_pathV = s_pathV_sim(:, 1:cS.T_sim);
            fprintf('   ✅ 已将存活率路径 cS.s_path_all 和 cS.s_pathV 的长度裁剪为 %d 期。\n', cS.T_sim);

            total_pop_path = sum(Z_path_raw_final, 1);
            pop_growth_rate_path_period = (total_pop_path(2:end) ./ total_pop_path(1:end-1)) - 1;
            pop_growth_rate_path_annual = (1 + pop_growth_rate_path_period).^(1/cS.time_Step) - 1;
            if converged && ~isempty(pop_growth_rate_path_annual), cS.n_ss = pop_growth_rate_path_annual(end); else, cS.n_ss = 0; end

            

            % --- 8. [!!! 核心新增: 调用稳态检测器 !!!] ---
            if converged
                % 准备检测器所需的输入
                Z_sim_ss_final = Z_out(:, end); % 模拟出的最终稳态分布
                s_pathV_ss_final = cS.s_pathV(:, end); % 最终的稳态生存率
                n_ss_annual_final = cS.n_ss; % 模拟出的最终稳态年增长率

                % 调用检测器
                population.check_ss_population_consistency(Z_sim_ss_final, s_pathV_ss_final, n_ss_annual_final, cS);
            end

            if graph_flag
                annual_years_vec = cS.start_year:cS.end_year;
                final_annual_mask = annual_years_full_range >= cS.start_year & annual_years_full_range <= cS.end_year;
                if any(final_annual_mask)
                    birth_rate_path_annual_for_plot = birth_rate_path_annual_full(final_annual_mask);
                else
                    birth_rate_path_annual_for_plot = []; % Handle case where no data falls in the plot range
                end
                population.create_Z_dynamics_plots(cS, Z_out, ...
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
            filePath = '../data/人口/WPP_MORT_CHINA_2023_2100_matrix.xlsx';

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

                function create_Z_dynamics_plots(cS, Z_out, birth_rate_path_annual, annual_years_vec, pop_growth_rate_path_annual)
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


            sgtitle(sprintf('外生路径概览 (%d-%d年)', cS.start_year, cS.end_year), 'FontSize', 16, 'FontWeight', 'bold');
            fprintf('✅ 外生动态可视化完成。\n');
        end


        function [Z_theory, birth_rate_period] = compute_theoretical_ss_dist(s_pathV_ss, n_ss_annual, time_Step, aD_new)
            % =========================================================================
            % == 函数: compute_theoretical_ss_dist (v1.1 - 返回出生率版)
            % == 功能: 计算具有恒定人口年增长率n的理论稳态人口分布。
            % == 新增:
            % ==   - 返回第二个参数 `birth_rate_period`，即模型周期内的粗出生率
            % ==     (新生儿占总人口的比例)。
            % =========================================================================

            % 诊断：检查输入的生存率是否合理
            if any(s_pathV_ss(1:end-1) <= 0) || any(s_pathV_ss(1:end-1) > 1)
                warning('compute_theoretical_ss_dist:InvalidSurvivalRates', '输入的生存率s_pathV_ss包含无效值(<=0 或 >1)。');
            end
            if ~issorted(s_pathV_ss, 'descend')
                warning('compute_theoretical_ss_dist:NonMonotonic', '输入的生存率s_pathV_ss不是单调递减的，这可能导致不真实的年龄分布。');
            end

            n_period = (1 + n_ss_annual)^time_Step - 1;
            mass_levels_by_age = zeros(aD_new, 1);
            mass_levels_by_age(1) = 1.0; % 将新生儿设为基准 1
            for ia = 1:(aD_new - 1)
                mass_levels_by_age(ia+1) = (mass_levels_by_age(ia) * s_pathV_ss(ia)) / (1 + n_period);
            end

            % 归一化得到人口分布
            total_mass = sum(mass_levels_by_age);
            if total_mass > 1e-9
                Z_theory = mass_levels_by_age / total_mass;
            else
                Z_theory = zeros(aD_new, 1);
            end

            % 新生儿的占比就是模型周期的粗出生率
            birth_rate_period = Z_theory(1);
        end

        function check_ss_population_consistency(Z_sim_ss, s_pathV_ss, n_sim_ss_annual, cS, tol)
            % =========================================================================
            % == 函数: check_ss_population_consistency (v1.0)
            % == 目的: 交叉验证一个模拟出的稳态人口分布是否与理论稳态一致。
            % ==
            % == 流程:
            % ==   1. 接收一个模拟收敛后的人口分布 (Z_sim_ss)。
            % ==   2. 接收与之对应的生存率 (s_pathV_ss) 和人口年增长率 (n_sim_ss_annual)。
            % ==   3. 使用 utils.compute_theoretical_ss_dist 根据生存率和增长率，
            % ==      独立地、从理论上计算出应有的稳态分布 (Z_theory)。
            % ==   4. 比较 Z_sim_ss 和 Z_theory，如果差异过大则发出警告。
            % =========================================================================
            if nargin < 5
                tol = 1e-7; % 默认容忍度
            end

            fprintf('   [交叉验证] 正在启动稳态人口分布一致性检测器...\n');

            % --- 步骤 1: 根据模拟出的最终增长率和生存率，计算理论稳态分布 ---
            [Z_theory, ~] = population.compute_theoretical_ss_dist(s_pathV_ss, n_sim_ss_annual, cS.time_Step, cS.aD_new);

            % --- 步骤 2: 计算两个分布之间的最大绝对差异 ---
            max_abs_diff = max(abs(Z_sim_ss - Z_theory));

            fprintf('      模拟得到的最终年增长率 (n_ss) ....: %.8f\n', n_sim_ss_annual);
            fprintf('      模拟分布与理论分布的最大绝对差异 .: %.4e\n', max_abs_diff);

            % --- 步骤 3: 报告结果 ---
            if max_abs_diff < tol
                fprintf('   [✅ 验证成功] 模拟收敛的稳态人口分布与理论稳态完全一致。\n');
            else
                warning('generate_exo_paths:InconsistentSS', ...
                    ['[❌ 验证失败] 模拟收敛的稳态人口分布与理论稳态不一致！\n' ...
                    '差异 (%.4e) 已超出容忍度 (%.4e)。\n' ...
                    '这可能意味着前向模拟过程中的数值累积误差或逻辑不一致。\n' ...
                    '请检查 generate_exo_paths 中的人口迭代恒等式。'], max_abs_diff, tol);
            end
        end


    end
end