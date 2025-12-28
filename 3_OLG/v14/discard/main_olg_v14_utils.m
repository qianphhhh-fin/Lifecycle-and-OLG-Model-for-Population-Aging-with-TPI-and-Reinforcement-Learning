% =========================================================================
% == 在 main_olg_v14_utils.m 文件中，请使用以下代码块 ==
% == 替换或添加相应的函数。                            ==
% =========================================================================

classdef main_olg_v14_utils

    % ... 这里是您从v14复制过来的所有静态方法 ...

    methods (Static)

        % --- [v14 新增] 从Excel加载和处理外生路径 ---
        function [Z_path, A_path, T_sim] = load_exogenous_paths(cS)
            fprintf('--- [v14] 加载和处理外生数据路径 ---\n');

            % 1. 加载人口路径数据
            pop_data = readtable('..\data\人口\population_by_age_group_all_years.xlsx', 'Sheet', 'pop_normalized');

            % 模拟年份设置
            cS.start_year = 1997;
            cS.end_year = 2102;
            sim_years = cS.start_year:cS.time_Step:cS.end_year;
            T_sim = length(sim_years);

            Z_path_raw = zeros(cS.aD_new, T_sim);

            for t = 1:T_sim
                year_t = sim_years(t);
                % 找到数据中最接近的年份列
                [~, col_idx] = min(abs(str2double(pop_data.Properties.VariableNames(2:end)) - year_t));
                % +1 是因为第一列是年龄组标签
                Z_path_raw(:, t) = pop_data{:, col_idx + 1};
            end

            % 归一化人口分布
            Z_path = Z_path_raw ./ sum(Z_path_raw, 1);
            fprintf('✅ 人口路径加载完成 (1997-%d, %d个时期)。\n', cS.end_year, T_sim);

            % 2. 加载和处理TFP路径数据
            tfp_data_pwt = readtable('..\data\PWT\china_pwt_data.xlsx');

            % 白重恩(2017)对劳动生产率增长率的预测 (作为TFP增长率的代理)
            % 2021-25: 5.57%, 26-30: 4.82%, 31-35: 3.94%, 36-40: 3.40%, 41-45: 3.46%, 46-50: 2.98%
            % 之后收敛到长期增长率
            bai_projections = [
                2020, 0.0628; % 2016-2020年的均值
                2025, 0.0557;
                2030, 0.0482;
                2035, 0.0394;
                2040, 0.0340;
                2045, 0.0346;
                2050, 0.0298;
                2102, 0.0150; % 假设长期年均增长率为1.5%
                ];

            % 计算年度TFP增长率
            g_A_annual = zeros(cS.end_year - cS.start_year + 1, 1);

            % 使用PWT数据 (1997-2019)
            tfp_pwt_series = tfp_data_pwt.rtfpna(tfp_data_pwt.year >= cS.start_year & tfp_data_pwt.year <= 2019);
            for i = 1:(length(tfp_pwt_series)-1)
                g_A_annual(i) = tfp_pwt_series(i+1) / tfp_pwt_series(i) - 1;
            end

            % 使用预测数据 (2020-2102)
            current_year_idx = 2020 - cS.start_year; % g_A_annual(23) -> 2019年的增长率
            g_A_annual(current_year_idx) = (tfp_pwt_series(end) / tfp_pwt_series(end-1)) -1; % 2019年的增长率

            for i = 1:size(bai_projections, 1)-1
                proj_cS.start_year = bai_projections(i, 1);
                proj_cS.end_year = bai_projections(i+1, 1);
                proj_rate = bai_projections(i, 2);

                for year = (proj_cS.start_year + 1):proj_cS.end_year
                    if year <= cS.end_year
                        idx = year - cS.start_year + 1;
                        g_A_annual(idx) = proj_rate;
                    end
                end
            end

            % 将年度增长率转换为5年期增长因子
            g_A_5year_factor = zeros(T_sim, 1);
            for t = 1:T_sim
                year_t = sim_years(t);
                factor = 1.0;
                % 计算从 t 到 t+5 年的累积增长
                for i = 0:cS.time_Step-1
                    annual_idx = (year_t + i) - cS.start_year + 1;
                    if annual_idx <= length(g_A_annual)
                        factor = factor * (1 + g_A_annual(annual_idx));
                    else
                        factor = factor * (1 + bai_projections(end,2)); % 使用长期增长率
                    end
                end
                g_A_5year_factor(t) = factor;
            end

            % 构建TFP水平路径 (A_1997 = 1)
            A_path = ones(T_sim, 1);
            for t = 1:(T_sim - 1)
                A_path(t+1) = A_path(t) * g_A_5year_factor(t);
            end
            fprintf('✅ TFP路径构建完成 (1997-2102)。\n');
        end


        function [cS] = calcaulte_theta_payg_path(cS, graph_flag)
            % [v14.4 - 已修正插值逻辑]
            % 核心修正：从插值【绝对人数】改为直接插值【覆盖率】，以保证路径的平滑和目标的准确达成。

            fprintf('正在构建基于【覆盖率插值】的有效养老金缴费率路径 (theta_path)...\n');

            % --- 1. 定义关键参数 ---
            theta_urban_employee_effective = 0.20; % 城镇职工体系的名义缴费率
            theta_resident_effective = 0.03;       % 城乡居民体系的名义缴费率 (近似值)

            % [政策目标] 定义未来的目标覆盖率
            coverage_urban_final = 0.8;    % 目标：城镇职工覆盖率达到劳动人口的60%
            coverage_resident_final = 0.8; % 目标：城乡居民覆盖率达到劳动人口的35%
            year_reach_final = 2050;        % 假设在2050年达到并维持这个目标覆盖率

            % --- 2. 收集来自官方统计的真实数据点 (单位：万人) ---
            year_pax_urban = [1997, 2000, 2005, 2008, 2011, 2014, 2017, 2020, 2023];
            pax_urban      = [8671, 10448, 17487, 21890, 28391, 34124, 40293, 45638, 52112];

            year_pax_resident = [2008, 2011, 2014, 2017, 2020, 2023];
            pax_resident      = [0,    45586, 50075, 51255, 54227, 54475];

            year_laborpop_data = [1997, 2000, 2005, 2008, 2011, 2014, 2017, 2020, 2023];
            laborpop_data      = [78531, 82352, 87532, 90833, 94072, 92837, 91199, 89552, 86481];

            % --- 3. [核心修正] 计算历史数据点的【覆盖率】 ---
            % a. 计算城镇职工历史覆盖率
            coverage_urban_data = pax_urban ./ laborpop_data;

            % b. 计算城乡居民历史覆盖率
            %    需要先插值得到对应年份的劳动人口
            laborpop_for_resident_years = interp1(year_laborpop_data, laborpop_data, year_pax_resident, 'linear', 'extrap');
            coverage_resident_data = pax_resident ./ laborpop_for_resident_years;

            % --- 4. [核心修正] 直接对【覆盖率】进行插值，生成完整的年度路径 ---
            annual_years_vector = cS.start_year:cS.end_year;

            % a. 插值生成年度【城镇职工覆盖率】路径
            interp_years_urban = [year_pax_urban, year_reach_final, cS.end_year];
            interp_coverage_urban = [coverage_urban_data, coverage_urban_final, coverage_urban_final]; % 达到目标后保持不变
            [unique_years_u, ia_u, ~] = unique(interp_years_urban);
            coverage_urban_annual = interp1(unique_years_u, interp_coverage_urban(ia_u), annual_years_vector, 'linear');

            % b. 插值生成年度【城乡居民覆盖率】路径
            interp_years_resident = [cS.start_year, year_pax_resident, year_reach_final, cS.end_year];
            interp_coverage_resident = [0, coverage_resident_data, coverage_resident_final, coverage_resident_final]; % 达到目标后保持不变
            [unique_years_r, ia_r, ~] = unique(interp_years_resident);
            coverage_resident_annual = interp1(unique_years_r, interp_coverage_resident(ia_r), annual_years_vector, 'linear');
            coverage_resident_annual(annual_years_vector < min(year_pax_resident)) = 0; % 确保制度开始前为0

            % --- 5. 基于【覆盖率路径】计算最终的【年度有效缴费率】路径 ---
            theta_path_annual = (coverage_urban_annual * theta_urban_employee_effective) + ...
                (coverage_resident_annual * theta_resident_effective);

            % --- 6. 从年度路径中，提取出模型【5年期】所需的路径 ---
            model_year_indices = 1:(cS.end_year - cS.start_year + 1);
            model_year_indices_5yr = model_year_indices(1:cS.time_Step:end);

            T_sim = length(cS.start_year:cS.time_Step:cS.end_year);
            theta_path = theta_path_annual(model_year_indices_5yr);
            theta_path = theta_path(1:min(T_sim, length(theta_path))); % 确保长度正确

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
                ylim([0, max(cS.theta_path)*1.2]); % 调整Y轴范围
            end
            fprintf('✅ 基于覆盖率插值的缴费率路径构建完成。\n');
        end


        function cS = ParameterValues_HuggettStyle()
            % [VFI最终对齐版 v2] - 添加了均衡求解器参数
            cS.pps_active = true;
            % --- 人口结构基础参数 ---
            cS.age1_orig = 20;
            cS.ageLast_orig = 98;
            cS.ageRetire_orig = 65;
            cS.aD_orig = cS.ageLast_orig - cS.age1_orig + 1;
            cS.aR_idx_orig = cS.ageRetire_orig - cS.age1_orig + 1;
            cS.aW_orig = cS.aR_idx_orig - 1;
            cS.physAgeV_orig = (cS.age1_orig : cS.ageLast_orig)';

            % --- 年龄组聚合参数 ---
            cS.time_Step = 5;
            cS.aD_new = ceil(cS.aD_orig / cS.time_Step);
            cS.aR_new = ceil(cS.aW_orig / cS.time_Step);

            cS.physAgeMap = cell(cS.aD_new, 1);
            for a = 1:cS.aD_new
                startIdx = (a-1)*cS.time_Step + 1;
                endIdx = min(a*cS.time_Step, cS.aD_orig);
                cS.physAgeMap{a} = startIdx:endIdx;
            end

            d_orig_data = [0.00159,0.00169,0.00174,0.00172,0.00165,0.00156,0.00149,0.00145,0.00145,0.00149,0.00156,0.00163,0.00171,0.00181,0.00193,0.00207,0.00225,0.00246,0.00270,0.00299,0.00332,0.00368,0.00409,0.00454,0.00504,0.00558,0.00617,0.00686,0.00766,0.00865,0.00955,0.01058,0.01162,0.01264,0.01368,0.01475,0.01593,0.01730,0.01891,0.02074,0.02271,0.02476,0.02690,0.02912,0.03143,0.03389,0.03652,0.03930,0.04225,0.04538,0.04871,0.05230,0.05623,0.06060,0.06542,0.07066,0.07636,0.08271,0.08986,0.09788,0.10732,0.11799,0.12895,0.13920,0.14861,0.16039,0.17303,0.18665,0.20194,0.21877,0.23601,0.25289,0.26973,0.28612,0.30128,0.31416,0.32915,0.34450,0.36018];
            s_orig = 1 - d_orig_data;
            cS.s_1yr_transitionV = zeros(cS.aD_new, 1);
            for a = 1:(cS.aD_new - 1)
                lastYearIdxInGroup = cS.physAgeMap{a}(end);
                cS.s_1yr_transitionV(a) = s_orig(lastYearIdxInGroup);
            end

            % --- [对齐] 劳动效率冲击过程参数 ---
            cS.lePersistence = 0.9;
            cS.leShockStd = 0.1^0.5;
            cS.nw = 5;

            % --- [对齐] 年龄效率剖面 ---
            ageEffV_orig_temp = zeros(100, 1);
            ageEffV_orig_temp(20:72) = [linspace(0.3,1.5,36-20+1), 1.5.*ones(1,47-37+1), linspace(1.5,0.2,65-48+1), linspace(0.18,0,72-66+1)];
            ageEffV_orig = ageEffV_orig_temp(cS.age1_orig : cS.ageLast_orig);
            cS.ageEffV_new = zeros(cS.aD_new, 1);
            for a = 1:cS.aD_new
                cS.ageEffV_new(a) = mean(ageEffV_orig(cS.physAgeMap{a}));
            end


            % =================================================================

        end

        function popS = initPopulation(cS)
            % initPopulation - 初始化人口结构
            popS = struct();

            % 基于2023年中国人口结构的初始分布 (16个年龄组)
            initial_pop_dist = [76.2, 86.4, 113.8, 98.6, 86.6, 102.7, 112.0, 99.0, 64.0, 66.9, 44.1, 25.4, 14.9, 6.8, 1.7, 0.2];

            initial_total = sum(initial_pop_dist);

            if initial_total > 0 && length(initial_pop_dist) == cS.aD_new
                popS.Z = (initial_pop_dist / initial_total * 100)';
            else
                warning('初始人口数据不匹配或总和为零。将设置为均匀的初始年龄组人口分布。');
                popS.Z = ones(cS.aD_new, 1) * (100 / cS.aD_new);
            end

            popS.totalPop = sum(popS.Z(:, 1));

            if popS.totalPop(1) > 1e-9
                popS.ageDist = popS.Z(:, 1) / popS.totalPop(1);
            else
                popS.ageDist = zeros(cS.aD_new, 1);
            end

            popS.initialAgeDist = popS.ageDist;
            fprintf('初始年龄组人口已设置。总人口=%.2f (代表百分比基数)。\n', popS.totalPop(1));
        end

        function popS = populationDynamics(popS, cS)
            % populationDynamics - 模拟人口动态演进

            % 基于中国数据的年龄组间存活率
            % beta_surv_pop = [0.998, 0.996, 0.994, 0.992, 0.988, 0.984, 0.980, 0.976, ...
            %     0.970, 0.960, 0.945, 0.920, 0.880, 0.800, 0.680];
            % cS.survivalProbV_popdyn = [beta_surv_pop, 0]'; % 最后一个年龄组存活率为0

            max_periods = 50;
            bgp_tolerance = 0.001;
            bgp_window = 5;

            Z_history = zeros(cS.aD_new, max_periods + 1);
            totalPop_history = zeros(max_periods + 1, 1);
            ageDist_history = zeros(cS.aD_new, max_periods + 1);

            Z_history(:, 1) = popS.Z;
            totalPop_history(1) = popS.totalPop;
            ageDist_history(:, 1) = popS.ageDist;

            fprintf('人口动态模拟开始 (年龄组, 最大期数 = %d)...\n', max_periods);
            bgp_reached = false;

            for t = 1:max_periods
                Z_current = Z_history(:, t);
                Z_next = zeros(cS.aD_new, 1);

                % 时变的人口增长率
                if t < 5
                    growth_rate = -0.01 - 0.003 * t;
                else
                    growth_rate = -0.03 - 0.004 * min(t - 5, 10);
                end

                Z_next(1) = Z_current(1) * (1 + growth_rate);

                for a = 2:cS.aD_new
                    Z_next(a) = Z_current(a-1) * cS.s_1yr_transitionV(a-1);
                end

                Z_history(:, t+1) = Z_next;
                totalPop_history(t+1) = sum(Z_next);
                if totalPop_history(t+1) > 1e-9
                    ageDist_history(:, t+1) = Z_next / totalPop_history(t+1);
                end

                if t >= bgp_window
                    is_stable = true;
                    for w = 0:bgp_window-1
                        change = norm(ageDist_history(:, t+1-w) - ageDist_history(:, t-w));
                        if change >= bgp_tolerance
                            is_stable = false;
                            break;
                        end
                    end
                    if is_stable
                        fprintf('人口稳态在模拟期数 %d 达到。\n', t);
                        bgp_reached = true;
                        Z_history = Z_history(:, 1:t+1);
                        totalPop_history = totalPop_history(1:t+1);
                        ageDist_history = ageDist_history(:, 1:t+1);
                        break;
                    end
                end
            end

            popS.Z = Z_history;
            popS.totalPop = totalPop_history;
            popS.ageDist = ageDist_history;

            if ~bgp_reached
                warning('人口稳态未在 %d 期内达到。', max_periods);
            end
        end

        function [Z_ss, dependency_ratio_ss, bgp_reached, bgp_period] = detectSteadyStatePopulation(popS, cS)
            % detectSteadyStatePopulation - 从模拟的人口历史中检测并提取稳态分布

            n_periods = size(popS.Z, 2);
            bgp_window = 5;
            bgp_tolerance = 0.001;
            bgp_reached = false;
            bgp_period = n_periods;

            if n_periods < bgp_window + 1
                warning('人口模拟期数过短，无法进行稳态检查。');
                Z_ss = popS.Z(:, end);
            else
                for t = n_periods : -1 : bgp_window + 1
                    is_stable = true;
                    for w = 0:bgp_window-1
                        change = norm(popS.ageDist(:, t-w) - popS.ageDist(:, t-w-1));
                        if change >= bgp_tolerance
                            is_stable = false;
                            break;
                        end
                    end
                    if is_stable
                        bgp_reached = true;
                        bgp_period = t;
                        break;
                    end
                end
                Z_ss = popS.Z(:, bgp_period);
            end

            working_pop_ss = sum(Z_ss(1:cS.aR_new));
            retired_pop_ss = sum(Z_ss(cS.aR_new+1:end));
            dependency_ratio_ss = retired_pop_ss / working_pop_ss;

            % 绘图比较
            figure('Name', 'VFI: 初始 vs 稳态人口分布');
            bar_data = [popS.initialAgeDist * 100, (Z_ss / sum(Z_ss)) * 100];
            bar(bar_data);
            xlabel('年龄组');
            ylabel('占总人口百分比 (%)');
            title(sprintf('初始 vs 稳态人口分布 (稳态于第%d期)', bgp_period));
            legend('初始分布', '稳态分布', 'Location', 'best');
            grid on;
        end

        function [HHlaborM_group, L_total_eff_pc] = LaborSupply_Huggett(eIdxM_group, cS, paramS, Z_ss_norm_group)
            % LaborSupply_Huggett - 计算劳动供给（适配年龄组模拟）
            nSim = size(eIdxM_group, 1);
            HHlaborM_group = zeros(nSim, cS.aD_new);

            for a_group = 1:cS.aR_new % 只在工作年龄组计算
                eIdx_this_age = eIdxM_group(:, a_group);
                labor_eff = paramS.leGridV(eIdx_this_age);
                HHlaborM_group(:, a_group) = cS.ageEffV_new(a_group) * labor_eff;
            end

            mean_labor_per_working_group = mean(HHlaborM_group(:, 1:cS.aR_new), 1);
            L_total_eff_pc = sum(mean_labor_per_working_group' .* Z_ss_norm_group(1:cS.aR_new));
        end

        % =========================================================================
        % == 修正 get_prices_at_t 函数 ==
        % =========================================================================
        % =========================================================================
        % == 在 main_olg_v14_utils.m 中，做最后一次修改 ==
        % =========================================================================
        function M_prices = get_prices_at_t(K_total_t, L_t, A_t, cS)
            % [v14.5 最终简化版] 只计算价格，不处理人口和劳动供给。
            % 劳动供给 L_t 作为输入，由调用它的上层函数负责计算。

            % [核心修正] 确保 HHPrices_Huggett 调用正确
            [r_mkt_t, w_t] = main_olg_v14_utils.HHPrices_Huggett(K_total_t, L_t, A_t, cS);

            r_net_t = r_mkt_t * (1 - cS.tau_k);

            M_prices = struct();
            M_prices.K_total_t = K_total_t;
            M_prices.L_t = L_t;
            M_prices.Y_t = A_t * (K_total_t^cS.alpha) * (L_t^(1-cS.alpha));
            M_prices.w_t = w_t;
            M_prices.r_mkt_t = r_mkt_t;
            M_prices.r_net_t = r_net_t;
        end
        function cS = generateGrids(cS)
            % generateGrids - 根据当前的网格参数设置，重新生成资产网格。
            %
            % 输入：
            %   cS - 包含 nk, kMax, nkpps, kppsMax 等参数的结构体
            % 输出：
            %   cS - 更新了 kGridV 和 kppsGridV 的参数结构体

            % 重新生成非PPS资产网格 (kGridV)
            % --- [对齐] 资产网格参数 ---
            cS.tgKY = 3;
            cS.tgWage = (1-cS.alpha)*cS.A*((cS.tgKY/cS.A)^(cS.alpha/(1-cS.alpha)));

            cS.kMin = 0;
            cS.kMax = 15 * cS.tgWage;
            cS.kppsMin = 0;
            cS.kppsMax = cS.kMax / 2;

            power_k = 1.5; % 网格密度参数
            if cS.nk > 1
                kGridV_temp = cS.kMin + (cS.kMax - cS.kMin) * (linspace(0, 1, cS.nk).^power_k);
                kGridV_temp(1) = cS.kMin;
            elseif cS.nk == 1
                kGridV_temp = cS.kMin;
            else
                kGridV_temp = [];
            end
            cS.kGridV = kGridV_temp(:); % 确保是列向量

            % 重新生成PPS资产网格 (kppsGridV)
            power_kpps = 1.5; % PPS资产网格密度参数
            if cS.nkpps > 1
                kppsGridV_temp = cS.kppsMin + (cS.kppsMax - cS.kppsMin) * (linspace(0, 1, cS.nkpps).^power_kpps);
                kppsGridV_temp(1) = cS.kppsMin;
            elseif cS.nkpps == 1
                kppsGridV_temp = cS.kppsMin;
            else
                kppsGridV_temp = [];
            end
            cS.kppsGridV = kppsGridV_temp(:); % 确保是列向量

            % 输出确认信息
            % fprintf('网格已重新生成：nk=%d, nkpps=%d\n', cS.nk, cS.nkpps);
        end


        % =========================================================================
        % == 在 main_olg_v14_utils.m 中，请用这个版本替换 solve_steady_state_with_fund
        % =========================================================================
% =========================================================================
% == 在 main_olg_v14_utils.m 中，请用这个最终版本替换 solve_steady_state_with_fund
% =========================================================================
        function [ss, eq_found] = solve_steady_state_with_fund(Z_ss_norm, B_p_Y_ratio_target, cS, paramS, eIdxM)
    % [v14.6 最终收敛版] 在收敛后，基于最终的均衡K值重新计算所有量，以确保ss结构体完全自洽。

    fprintf('\n--- 开始求解初始稳态 (有目标基金 B_p, 无初始PPS) ---\n');

    K_guess = 1.8816;
    max_iter = 100;
    tol = 1e-4;
    damp = 0.5;
    eq_found = false;
    
    A_ss = cS.A;
    theta_ss = cS.theta_path(1);

    paramS_ss_for_L = paramS;
    paramS_ss_for_L.ageMassV = Z_ss_norm;
    [~, L_ss] = main_olg_v14_utils.LaborSupply_Huggett(eIdxM, cS, paramS_ss_for_L, Z_ss_norm);
    
    fprintf('%5s | %12s | %12s | %12s | %12s | %12s | %12s\n', ...
        'Iter', 'K_guess', 'K_model', 'K_pvt', 'B_p', 'r_net', 'K_error');
    fprintf('%s\n', repmat('-', 90, 1));

    K_total_model = 0; % Initialize

    for iter = 1:max_iter
        M_ss_prices = main_olg_v14_utils.get_prices_at_t(K_guess, L_ss, A_ss, cS);
        B_p_target = B_p_Y_ratio_target * M_ss_prices.Y_t;
        
        M_ss = M_ss_prices;
        M_ss.current_t = 1;
        
        mass_workers_ss = sum(Z_ss_norm(1:cS.aR_new));
        avg_wage_ss = M_ss.w_t * M_ss.L_t / mass_workers_ss;
        M_ss.b_t = cS.rho_prime_payg * avg_wage_ss;

        cS_ss = cS;
        cS_ss.pps_active = false;
        cS_ss.theta_t = theta_ss;
        
        B_g_t = cS.gov_debt_frac_Y * M_ss_prices.Y_t;
        
        [K_pvt_model, K_pps_model, ~, ~, ~, ~, ~] = ...
            main_olg_v14_utils.simulate_private_capital_forward(1, 0, B_p_target, 0, B_g_t, Z_ss_norm, A_ss, cS_ss, paramS, eIdxM);

        K_total_model = K_pvt_model + B_p_target + K_pps_model;
        K_error = K_guess - K_total_model;

        fprintf('%5d | %12.4f | %12.4f | %12.4f | %12.4f | %11.2f%% | %12.3e\n', ...
            iter, K_guess, K_total_model, K_pvt_model, B_p_target, M_ss_prices.r_net_t*100, K_error);

        if abs(K_error) < tol
            fprintf('✅ 初始稳态均衡收敛！\n');
            eq_found = true;
            break;
        end
        K_guess = (1 - damp) * K_guess + damp * K_total_model;
    end

    if ~eq_found, warning('初始稳态均衡未在最大迭代次数内收敛。'); end
    
    % --- [核心修正] 在循环结束后，使用最终收敛的 K_total_model 重新计算所有一切 ---
    fprintf('--- 基于最终均衡资本 K=%.4f 重新计算所有量以确保自洽性 ---\n', K_total_model);
    K_final = K_total_model;

    % 1. 重新计算价格和宏观总量
    M_final_prices = main_olg_v14_utils.get_prices_at_t(K_final, L_ss, A_ss, cS);
    B_p_final = B_p_Y_ratio_target * M_final_prices.Y_t;
    B_g_final = cS.gov_debt_frac_Y * M_final_prices.Y_t; % Use final Y

    % 2. 重新计算家庭决策和聚合量
    M_final = M_final_prices;
    M_final.current_t = 1;
    mass_workers_final = sum(Z_ss_norm(1:cS.aR_new));
    avg_wage_final = M_final.w_t * M_final.L_t / mass_workers_final;
    M_final.b_t = cS.rho_prime_payg * avg_wage_final;

    cS_final = cS;
    cS_final.pps_active = false;
    cS_final.theta_t = theta_ss;

    [K_pvt_final, K_pps_final, C_final, ~, ~, ~, ~] = ...
        main_olg_v14_utils.simulate_private_capital_forward(1, 0, B_p_final, 0, B_g_final, Z_ss_norm, A_ss, cS_final, paramS, eIdxM);
    
    % --- 打包最终的、完全自洽的稳态结果 ---
    ss = struct();
    ss.K_total = K_final;
    ss.K_pvt = K_pvt_final;
    ss.K_pps = K_pps_final;
    ss.B_p = B_p_final;
    ss.B_g = B_g_final;
    ss.L = M_final_prices.L_t;
    ss.Y = M_final_prices.Y_t;
    ss.w = M_final_prices.w_t;
    ss.r_mkt = M_final_prices.r_mkt_t;
    ss.r_net = M_final_prices.r_net_t;
    ss.C = C_final;
    ss.G = cS.gov_exp_frac_Y * ss.Y;

    % 最终自洽性检查 (可选，但推荐)
    final_check_K = ss.K_pvt + ss.K_pps + ss.B_p - ss.B_g;
    if abs(final_check_K - ss.K_total) > 1e-6
       warning('最终稳态打包存在不一致性！ K_assets(%.4f) != K_physical(%.4f)', final_check_K, ss.K_total);
    end
end
        
        function gbc_residual = check_gbc_residual(...
                K_t, C_t, Y_t, G_t, B_g_t, w_t, r_mkt_t, ...
                PensionRevenue_t, PensionOutlay_t, ... % [新] 传入养老金收支
                cS, paramS_t)
            % [v14.2 - 精确版GBC]
            % 明确区分一般预算和养老金预算，并处理可能的缺口。

            % 1. 计算一般财政收入 (General Revenue)
            %    a. 一般劳动税收入 (税基是总劳动收入)
            LaborTaxRev_general = cS.tau_l * w_t * paramS_t.L_per_capita;

            %    b. 资本税收入
            % [修正] 资本税基于税前市场回报率，不额外扣除折旧
            CapitalTaxRev = cS.tau_k * r_mkt_t * K_t;

            %    c. 消费税收入
            ConsumptionTaxRev = cS.tau_c * C_t;

            GeneralRevenue = LaborTaxRev_general + CapitalTaxRev + ConsumptionTaxRev;

            % 2. 计算一般财政支出 (General Outlays)
            %    a. 政府消费
            GovConsumption = G_t;

            %    b. 债务利息
            % [修正] 债务利息基于市场利率，不额外扣除折旧
            DebtService = r_mkt_t * B_g_t;

            %    c. [核心] 处理养老金缺口
            %       如果养老金有赤字 (支出 > 收入)，这个缺口必须由一般财政来弥补。
            pension_deficit_subsidy = max(0, PensionOutlay_t - PensionRevenue_t);

            GeneralOutlays = GovConsumption + DebtService + pension_deficit_subsidy;

            % 3. 计算预算缺口
            gbc_residual = GeneralRevenue - GeneralOutlays;
        end


        function [r_mkt_period, w_t, Y_period] = HHPrices_Huggett(K, L, A_t, cS)
            if K <= 0, K = 1e-8; end; if L <= 0, L = 1e-8; end
            Y_period = A_t * (K.^cS.alpha) .* (L.^(1-cS.alpha));
            MPK_period = cS.alpha * Y_period / K;
            w_t = (1-cS.alpha) * Y_period / L;
            r_mkt_period = MPK_period - cS.ddk;
        end

        function M_t = get_prices_and_policy_at_t(t, K_pvt_t, B_p_t, K_pps_t, B_g_t, Z_t_norm, A_t, cS, paramS, eIdxM)
            % [最终修正] 只计算利率r，不计算回报因子R
            K_physical_t = K_pvt_t + K_pps_t + B_p_t - B_g_t;
            if K_physical_t <= 0, K_physical_t = 1e-8; end
            paramS_t = paramS;
            paramS_t.ageMassV = Z_t_norm;
            [~, L_t] = main_olg_v14_utils.LaborSupply_Huggett(eIdxM, cS, paramS_t, Z_t_norm);
            [r_mkt_period, w_t, Y_period] = main_olg_v14_utils.HHPrices_Huggett(K_physical_t, L_t, A_t, cS);
            r_net_period = r_mkt_period * (1 - cS.tau_k);
            mass_workers_t = sum(Z_t_norm(1:cS.aR_new));
            if mass_workers_t > 0, avg_wage_period = w_t * L_t / mass_workers_t; else, avg_wage_period = 0; end
            b_period = cS.rho_prime_payg * avg_wage_period;

            M_t = struct();
            M_t.current_t = t; M_t.K_physical_t = K_physical_t; M_t.L_t = L_t;
            M_t.Y_t = Y_period; M_t.w_t = w_t; M_t.b_t = b_period;
            M_t.r_mkt_period = r_mkt_period; M_t.r_net_period = r_net_period;
            M_t.K_pvt_t = K_pvt_t; M_t.K_pps_t = K_pps_t; M_t.B_p_t = B_p_t; M_t.B_g_t = B_g_t;
        end

        function [B_p_next, PensionSurplus_period] = update_pension_fund(B_p_t, M_t, Z_t_norm, cS)
            theta_t = cS.theta_path(M_t.current_t);
            PensionRevenue_period = theta_t * M_t.w_t * M_t.L_t;
            mass_retirees_t = sum(Z_t_norm(cS.aR_new+1:end));
            PensionOutlay_period = M_t.b_t * mass_retirees_t;
            PensionSurplus_period = PensionRevenue_period - PensionOutlay_period;
            B_p_next = B_p_t * (1 + M_t.r_mkt_period) + PensionSurplus_period;
            B_p_next = max(0, B_p_next);
        end

 
        function [B_g_next, G_t, TotalTax_t] = update_gov_debt(B_g_t, M_t_complete, C_t, cS, PensionSurplus_period)
            % [最终审计通过版 v11 - 统一劳动税税基]
            % 核心修正: 劳动税税基 (w*L) 必须扣除所有养老金缴费 (强制PAYG 和 自愿PPS)，
            %           以与家庭微观决策的预算约束完全对齐。这是闭合Y=C+I+G的关键。

            % --- 1. 从 M_t 结构体中获取所有需要的当期变量 ---
            w_t              = M_t_complete.w_t;
            L_t              = M_t_complete.L_t;
            Y_t              = M_t_complete.Y_t;
            r_mkt_t          = M_t_complete.r_mkt_period;
            K_pvt_t          = M_t_complete.K_pvt_t;
            Total_Cpps_t     = M_t_complete.Total_Cpps_t; % Aggregate PPS contributions
            Total_PpsTax_t   = M_t_complete.Total_PpsTax_t;
            theta_t          = cS.theta_path(M_t_complete.current_t);

            % --- 2. 计算各项税收 (与微观层面完全对齐) ---

            % a. PAYG缴费总额 (这是给养老金体系的，不是一般税)
            payg_contributions_total = theta_t * w_t * L_t;

            % b. [核心修正] 计算劳动所得税 (Labor Income Tax)
            %    税基 = 总劳动收入 - PAYG缴费 - PPS缴费
            taxable_labor_income = max(0, w_t * L_t - payg_contributions_total - Total_Cpps_t);
            LaborTax_t = cS.tau_l * taxable_labor_income;

            % c. 资本利得税 (Capital Gains Tax)
            %    税基是家庭持有的、需纳税的私人资本 K_pvt_t (此逻辑保持不变)
            CapitalTax_t = cS.tau_k * r_mkt_t * K_pvt_t;

            % d. 消费税 和 PPS提现税 (Consumption and PPS Withdrawal Tax)
            ConsumptionTax_t = cS.tau_c * C_t;
            PpsWithdrawalTax_t = Total_PpsTax_t;

            % e. 总税收 (进入一般财政预算)
            TotalTax_t = LaborTax_t + CapitalTax_t + ConsumptionTax_t + PpsWithdrawalTax_t;

            % --- 3. 计算政府支出与主预算盈余 ---
            G_t = cS.gov_exp_frac_Y * Y_t;
            pension_deficit_subsidy = max(0, -PensionSurplus_period);
            PrimarySurplus_g = TotalTax_t - G_t - pension_deficit_subsidy;

            % --- 4. 更新政府债务 ---
            B_g_next = B_g_t * (1 + r_mkt_t) - PrimarySurplus_g;
        end
        
        function [pretax_non_capital_income, pps_contrib_expenditure] = HHIncome_Huggett(w_gross, tr_per_capita, b_payg_val, c_pps_chosen, a_idx, paramS_hh, cS, epsilon_val), pps_contrib_expenditure = c_pps_chosen; labor_income_gross = 0; if a_idx <= cS.aR_new, labor_income_gross = w_gross * cS.ageEffV_new(a_idx) * epsilon_val; end; pension_benefits = 0; if a_idx > cS.aR_new, pension_benefits = b_payg_val; end; pretax_non_capital_income = labor_income_gross + pension_benefits + tr_per_capita; end


        function [cPolM, kPolM, cPpsPolM, valM] = HHSolution_VFI_Huggett(M_vfi, tr_per_capita_vfi, paramS_vfi, cS_vfi)
            % 这个函数现在是标准的生命周期求解器，无需修改
            % 确保使用之前修正过的边界条件版本
            valM = -Inf(cS_vfi.nk, cS_vfi.nkpps, cS_vfi.nw, cS_vfi.aD_new);
            cPolM = zeros(cS_vfi.nk, cS_vfi.nkpps, cS_vfi.nw, cS_vfi.aD_new);
            kPolM = zeros(cS_vfi.nk, cS_vfi.nkpps, cS_vfi.nw, cS_vfi.aD_new);
            cPpsPolM = zeros(cS_vfi.nk, cS_vfi.nkpps, cS_vfi.nw, cS_vfi.aD_new);
            bV_payg_vfi = zeros(1, cS_vfi.aD_new);
            bV_payg_vfi(cS_vfi.aR_new+1:end) = M_vfi.b_t;

            for a_idx = cS_vfi.aD_new : -1 : 1
                vPrime_kkppse_next = [];
                if a_idx < cS_vfi.aD_new
                    vPrime_kkppse_next = valM(:,:,:,a_idx+1);
                end

                [cPolM(:,:,:,a_idx), kPolM(:,:,:,a_idx), cPpsPolM(:,:,:,a_idx), valM(:,:,:,a_idx)] = ...
                    main_olg_v14_utils.HHSolutionByAge_VFI_GridSearch(a_idx, vPrime_kkppse_next, M_vfi, tr_per_capita_vfi, bV_payg_vfi(a_idx), paramS_vfi, cS_vfi);
                % [cPolM(:,:,:,a_idx), kPolM(:,:,:,a_idx), cPpsPolM(:,:,:,a_idx), valM(:,:,:,a_idx)] = ...
                %     main_olg_v14_utils.HHSolutionByAge_NAIVE(a_idx,  M_vfi, tr_per_capita_vfi, bV_payg_vfi(a_idx), paramS_vfi, cS_vfi);
            end
        end
        % =========================================================================
        % == 在 main_olg_v14_utils.m 文件中，请在 HHSolutionByAge_VFI_GridSearch 旁边 ==
        % == 添加以下这个【新的】函数。                                        ==
        % =========================================================================

        function [cPol_age_q, kPol_age, cPpsPol_age_choice, val_age] = HHSolutionByAge_NAIVE(a_idx, M_age, tr_per_capita_age, b_age_val, paramS_age, cS)
            % [调试用] 一个简单的、透明的 "Naive" 决策策略
            % 工作期：税后非资本收入的一半消费，一半储蓄为k'。
            % 退休期：花光所有，k'=0。

            % 初始化输出矩阵
            val_age    = zeros(cS.nk, cS.nkpps, cS.nw); % val在这种策略下无意义
            cPol_age_q = zeros(cS.nk, cS.nkpps, cS.nw);
            kPol_age   = zeros(cS.nk, cS.nkpps, cS.nw);
            cPpsPol_age_choice = zeros(cS.nk, cS.nkpps, cS.nw); % Naive策略不使用PPS

            % 遍历所有状态点 (ik, ikpps, ie)
            for ik = 1:cS.nk
                for ikpps = 1:cS.nkpps
                    for ie = 1:cS.nw
                        k_state = cS.kGridV(ik);
                        k_pps_state = cS.kppsGridV(ikpps);
                        epsilon_state = paramS_age.leGridV(ie);

                        % --- [核心] 应用 Naive 决策规则 ---
                        if a_idx <= cS.aR_new % --- 工作期 ---
                            % 1. 计算税后非资本收入
                            labor_income_gross = M_age.w_t * cS.ageEffV_new(a_idx) * epsilon_state;
                            payg_tax = cS.theta_t * labor_income_gross;
                            labor_tax = cS.tau_l * labor_income_gross; % 简化：假设cpps=0

                            % 遗赠只给 a>1 的人
                            tr_this_age = 0; if a_idx > 1, tr_this_age = tr_per_capita_age; end

                            net_non_capital_income = labor_income_gross - payg_tax - labor_tax + tr_this_age;

                            % 2. 决策：一半储蓄，一半消费
                            k_prime_decision = 0.5 * net_non_capital_income;

                            % 3. 计算可用现金并反推消费
                            cash_in_hand = k_state * (1 + M_age.r_net_period) + net_non_capital_income;
                            c_expend = cash_in_hand - k_prime_decision;
                            c_decision = max(cS.cFloor, c_expend / (1 + cS.tau_c));

                            % 存储决策
                            kPol_age(ik, ikpps, ie) = k_prime_decision;
                            cPol_age_q(ik, ikpps, ie) = c_decision;

                        else % --- 退休期 ---
                            % 决策：花光所有，不留储蓄
                            kPol_age(ik, ikpps, ie) = cS.kMin; % k'=0

                            % 反推消费 (代码与 HHSolutionByAge_VFI_GridSearch 退休期部分一致)
                            pps_withdrawal_gross = k_pps_state * cS.pps_withdrawal_rate;
                            inflows_pps_withdrawal_net = pps_withdrawal_gross * (1 - cS.pps_tax_rate_withdrawal);
                            total_inflows = k_state * (1 + M_age.r_net_period) + b_age_val + tr_per_capita_age + inflows_pps_withdrawal_net;

                            net_cash = total_inflows;
                            k_prime_choice = kPol_age(ik, ikpps, ie); % 就是 0
                            c_expend = net_cash - k_prime_choice;
                            c_decision = max(cS.cFloor, c_expend / (1 + cS.tau_c));

                            cPol_age_q(ik, ikpps, ie) = c_decision;
                        end
                    end
                end
            end
        end

       
        function [cPol_age_q, kPol_age, cPpsPol_age_choice, val_age] = HHSolutionByAge_VFI_GridSearch(a_idx, vPrime_kkppse_next, M_age, tr_per_capita_age, b_age_val, paramS_age, cS)
            % [最终审计通过版 v16 - 统一劳动税税基]
            % 核心修正: 劳动税(labor_tax)的计算现在也扣除了强制PAYG缴费(payg_tax)，
            %           以确保家庭的预算约束与宏观GBC完全对齐。

            val_age    = -Inf(cS.nk, cS.nkpps, cS.nw);
            cPol_age_q = zeros(cS.nk, cS.nkpps, cS.nw);
            kPol_age   = zeros(cS.nk, cS.nkpps, cS.nw);
            cPpsPol_age_choice = zeros(cS.nk, cS.nkpps, cS.nw);

            % --- 最后一期逻辑 (不变) ---
            if a_idx == cS.aD_new, [K_grid, Kpps_grid] = ndgrid(cS.kGridV, cS.kppsGridV); [pretax_non_capital_income, ~] = main_olg_v14_utils.HHIncome_Huggett(0, tr_per_capita_age, b_age_val, 0, a_idx, paramS_age, cS, 0); k_end_value = K_grid .* (1 + M_age.r_net_period); pps_liquidation_net = Kpps_grid .* (1 + M_age.r_mkt_period + cS.pps_return_rate_premium) * (1 - cS.pps_tax_rate_withdrawal); total_resources = k_end_value + pps_liquidation_net + pretax_non_capital_income; final_bequest = zeros(size(total_resources)); util_b = 0; if isfield(cS, 'phi_bequest') && cS.phi_bequest > 1e-9, phi_adj = cS.phi_bequest * (1 + cS.tau_c)^(-1); omega = (1 + phi_adj^(-1/cS.sigma))^(-1); final_c_expenditure = omega .* total_resources; final_bequest = (1 - omega) .* total_resources; [~, util_b] = main_olg_v14_utils.CES_utility(final_bequest, cS.sigma_bequest, cS); else, final_c_expenditure = total_resources; end; final_c = max(cS.cFloor, final_c_expenditure / (1 + cS.tau_c)); [~, util_c] = main_olg_v14_utils.CES_utility(final_c, cS.sigma, cS); final_v = util_c; if isfield(cS, 'phi_bequest'), final_v = final_v + cS.phi_bequest * util_b; end; for ie = 1:cS.nw, cPol_age_q(:,:,ie) = final_c; val_age(:,:,ie) = final_v; kPol_age(:,:,ie) = final_bequest; cPpsPol_age_choice(:,:,ie) = 0; end; return; end

            % --- 非最后一期逻辑 (插值器部分不变) ---
            EV_matrix = zeros(cS.nk, cS.nkpps, cS.nw); for ie_current = 1:cS.nw, transition_probs = paramS_age.leTrProbM(ie_current, :); vPrime_reshaped = reshape(vPrime_kkppse_next, [cS.nk * cS.nkpps, cS.nw]); EV_slice = vPrime_reshaped * transition_probs'; EV_matrix(:, :, ie_current) = reshape(EV_slice, [cS.nk, cS.nkpps]); end; EV_interpolants = cell(cS.nw, 1); is_pps_disabled = (cS.nkpps == 1); for ie_current = 1:cS.nw, if is_pps_disabled, EV_interpolants{ie_current} = griddedInterpolant(cS.kGridV, squeeze(EV_matrix(:, 1, ie_current)), 'linear', 'linear'); else, EV_interpolants{ie_current} = griddedInterpolant({cS.kGridV, cS.kppsGridV}, EV_matrix(:, :, ie_current), 'linear', 'linear'); end, end;

            for ik = 1:cS.nk, for ikpps = 1:cS.nkpps, for ie = 1:cS.nw
                        k_state = cS.kGridV(ik); k_pps_state = cS.kppsGridV(ikpps); epsilon_state = paramS_age.leGridV(ie);
                        ev_interpolant = EV_interpolants{ie};
                        best_val = -Inf; best_c = cS.cFloor; best_k_prime = cS.kMin; best_c_pps = 0;

                        tr_this_age = 0; if a_idx > 1, tr_this_age = tr_per_capita_age; end

                        if a_idx <= cS.aR_new % --- 工作期 ---
                            labor_income_gross = M_age.w_t * cS.ageEffV_new(a_idx) * epsilon_state;
                            max_cpps = 0; if cS.pps_active, max_cpps = min(cS.pps_contrib_limit, labor_income_gross * cS.pps_max_contrib_frac); end
                            cpps_grid = linspace(0, max(0, max_cpps), cS.npps);
                            payg_tax = cS.theta_t * labor_income_gross; % PAYG tax is fixed for this state

                            for c_pps_choice = cpps_grid
                                k_return = k_state * (1 + M_age.r_net_period);
                                total_inflow = k_return + labor_income_gross + tr_this_age;

                                % [核心修正] 劳动税基现在也扣除了PAYG缴费
                                labor_tax = cS.tau_l * max(0, labor_income_gross - c_pps_choice - payg_tax);
                                cpps_outflow = c_pps_choice;

                                net_cash_for_c_and_k_prime = total_inflow - (payg_tax + labor_tax + cpps_outflow);

                                k_prime_max = net_cash_for_c_and_k_prime - cS.cFloor * (1 + cS.tau_c);
                                if k_prime_max < cS.kMin, k_prime_grid = [cS.kMin]; else, k_prime_grid = linspace(cS.kMin, k_prime_max, cS.nkprime); end

                                for k_prime_choice = k_prime_grid
                                    c_expend = net_cash_for_c_and_k_prime - k_prime_choice;
                                    if c_expend < cS.cFloor * (1 + cS.tau_c), continue; end
                                    c_choice = c_expend / (1 + cS.tau_c);
                                    [~, util] = main_olg_v14_utils.CES_utility(c_choice, cS.sigma, cS);

                                    pps_return_factor = 1 + M_age.r_mkt_period + cS.pps_return_rate_premium;
                                    k_pps_prime = (k_pps_state + c_pps_choice) * pps_return_factor;
                                    k_prime_clamped = max(cS.kGridV(1), min(cS.kGridV(end), k_prime_choice));
                                    k_pps_prime_clamped = max(cS.kppsGridV(1), min(cS.kppsGridV(end), k_pps_prime));

                                    if is_pps_disabled, ev = ev_interpolant(k_prime_clamped); else, ev = ev_interpolant(k_prime_clamped, k_pps_prime_clamped); end
                                    current_val = util + cS.beta * cS.s_1yr_transitionV(a_idx) * ev;
                                    if current_val > best_val, best_val = current_val; best_c = c_choice; best_k_prime = k_prime_choice; best_c_pps = c_pps_choice; end
                                end
                            end

                        else % --- 退休期 (逻辑不变) ---
                            c_pps_choice = 0;
                            k_return = k_state * (1 + M_age.r_net_period);
                            pps_withdrawal_gross = 0; if cS.pps_active, pps_withdrawal_gross = k_pps_state * cS.pps_withdrawal_rate; end
                            pps_withdrawal_net = pps_withdrawal_gross * (1 - cS.pps_tax_rate_withdrawal);
                            total_inflow = k_return + b_age_val + tr_this_age + pps_withdrawal_net;
                            net_cash = total_inflow;
                            k_prime_max = net_cash - cS.cFloor * (1 + cS.tau_c);

                            if k_prime_max < cS.kMin, k_prime_grid = [cS.kMin]; else, k_prime_grid = linspace(cS.kMin, k_prime_max, cS.nkprime); end

                            for k_prime_choice = k_prime_grid
                                c_expend = net_cash - k_prime_choice;
                                if c_expend < cS.cFloor * (1 + cS.tau_c), continue; end
                                c_choice = c_expend / (1 + cS.tau_c);
                                [~, util] = main_olg_v14_utils.CES_utility(c_choice, cS.sigma, cS);

                                pps_return_factor = 1 + M_age.r_mkt_period + cS.pps_return_rate_premium;
                                k_pps_prime = (k_pps_state - pps_withdrawal_gross) * pps_return_factor;
                                k_prime_clamped = max(cS.kGridV(1), min(cS.kGridV(end), k_prime_choice));
                                k_pps_prime_clamped = max(cS.kppsGridV(1), min(cS.kppsGridV(end), k_pps_prime));

                                if is_pps_disabled, ev = ev_interpolant(k_prime_clamped); else, ev = ev_interpolant(k_prime_clamped, k_pps_prime_clamped); end
                                current_val = util + cS.beta * cS.s_1yr_transitionV(a_idx) * ev;
                                if current_val > best_val, best_val = current_val; best_c = c_choice; best_k_prime = k_prime_choice; best_c_pps = c_pps_choice; end
                            end
                        end
                        val_age(ik, ikpps, ie) = best_val; cPol_age_q(ik, ikpps, ie) = best_c; kPol_age(ik, ikpps, ie) = best_k_prime; cPpsPol_age_choice(ik, ikpps, ie) = best_c_pps;
            end, end, end
        end
        % =========================================================================
        % == 在 main_olg_v14_utils.m 文件中，请使用以下代码块 ==
        % == 替换掉原有的 simulate_private_capital_forward 函数。==
        % =========================================================================

        % =========================================================================
        % == 在 main_olg_v14_utils.m 文件中，请使用以下代码块 ==
        % == 替换掉原有的 simulate_private_capital_forward 函数。==
        % =========================================================================
        % =========================================================================
        % == 在 main_olg_v14_utils.m 文件中，请使用以下代码块 ==
        % == 替换掉原有的 simulate_private_capital_forward 函数。==
        % =========================================================================
        % =========================================================================
        % == 在 main_olg_v14_utils.m 文件中，请使用以下代码块 ==
        % == 替换掉原有的 simulate_private_capital_forward 函数。==
        % == 这是最终的、将使模型平衡的正确版本。               ==
        % =========================================================================

        % =========================================================================
        % == 在 main_olg_v14_utils.m 文件中，请使用以下代码块 ==
        % == 替换掉原有的 simulate_private_capital_forward 函数。==
        % == 这是最终的、将使模型平衡的正确版本。               ==
        % =========================================================================
% =========================================================================
% == 在 main_olg_v14_utils.m 中，请用这个版本替换 simulate_private_capital_forward
% =========================================================================
function [K_pvt_next, K_pps_next, C_t, Total_Cpps_t, Total_PpsTax_t, total_accidental_bequest, M_t_complete] = simulate_private_capital_forward(t, K_pvt_t, B_p_t, K_pps_t, B_g_t, Z_t, A_t, cS, paramS, eIdxM)
    % [最终审计通过版 v18 - 稳态/转型路径兼容性修正]
    % 核心修正: 修复了在稳态求解中调用此函数时 cS.sim_years 字段不存在的问题。
    %           现在，函数会优先检查 cS.pps_active 是否已被外部设置(如在稳态求解器中)，
    %           如果未设置，才根据年份进行判断。

    % --- 1. 计算宏观价格和政策 ---
    K_physical_t = K_pvt_t + K_pps_t + B_p_t - B_g_t;
    if K_physical_t <= 0, K_physical_t = 1e-8; end
    
    paramS_t = paramS; 
    paramS_t.ageMassV = Z_t;
    [~, L_t] = main_olg_v14_utils.LaborSupply_Huggett(eIdxM, cS, paramS_t, Z_t);
    
    [r_mkt_period, w_t, Y_t] = main_olg_v14_utils.HHPrices_Huggett(K_physical_t, L_t, A_t, cS);
    r_net_period = r_mkt_period * (1 - cS.tau_k);
    
    mass_workers_t = sum(Z_t(1:cS.aR_new));
    if mass_workers_t > 0, avg_wage_period = w_t * L_t / mass_workers_t; else, avg_wage_period = 0; end
    b_period = cS.rho_prime_payg * avg_wage_period;

    M_t = struct();
    M_t.current_t = t; M_t.K_physical_t = K_physical_t; M_t.L_t = L_t; M_t.Y_t = Y_t;
    M_t.w_t = w_t; M_t.b_t = b_period; M_t.r_mkt_period = r_mkt_period; M_t.r_net_period = r_net_period;
    M_t.K_pvt_t = K_pvt_t; M_t.K_pps_t = K_pps_t; M_t.B_p_t = B_p_t; M_t.B_g_t = B_g_t;
    
    cS_vfi = cS; 
    
    % [核心修正] 检查 cS.pps_active 是否已由上层函数（如稳态求解器）预设
    if ~isfield(cS_vfi, 'pps_active')
        % 如果没有预设，则根据年份判断（用于转型路径）
        cS_vfi.pps_active = (cS.sim_years(t) >= cS_vfi.pps_activation_year);
    end
    % 如果 cS.pps_active 已经被设为 false (来自稳态求解器)，这里的逻辑将尊重该设置。
    
    if ~cS_vfi.pps_active
        cS_vfi.nkpps = 1; 
        cS_vfi.npps = 1; 
        cS_vfi = main_olg_v14_utils.generateGrids(cS_vfi); 
    end
    
    % 检查 cS.theta_t 是否已由上层函数（如稳态求解器）预设
    if ~isfield(cS_vfi, 'theta_t')
        cS_vfi.theta_t = cS.theta_path(t);
    end

    % --- 2. 完整的不动点迭代 (求解价值函数 V 和 遗赠 TR) (这部分不变) ---
    valM = zeros(cS_vfi.nk, cS_vfi.nkpps, cS_vfi.nw, cS_vfi.aD_new);
    TR_total_t_eq = 0;
    mass_survived = sum(Z_t .* cS_vfi.s_1yr_transitionV);
    
    for iter_main = 1:100 
        tr_per_survivor = 0;
        if mass_survived > 0, tr_per_survivor = TR_total_t_eq / mass_survived; end
        
        [~, kPolM, cPpsPolM, valM_new] = main_olg_v14_utils.HHSolution_VFI_Huggett(M_t, tr_per_survivor, paramS, cS_vfi);
        [kHistM, kPpsHistM, ~, ~] = main_olg_v14_utils.HHSimulation_olgm(kPolM, cPpsPolM, kPolM, eIdxM, M_t, paramS, cS_vfi, tr_per_survivor);

        mean_k_prime_by_age = mean(kHistM(:, 2:end), 1);
        mean_k_pps_prime_by_age = mean(kPpsHistM(:, 2:end), 1);
        prob_death_by_age = (1 - cS_vfi.s_1yr_transitionV);
        bequest_from_k = mean_k_prime_by_age' .* prob_death_by_age .* Z_t;
        bequest_from_k_pps = mean_k_pps_prime_by_age' .* prob_death_by_age .* Z_t;
        TR_total_t_new = sum(bequest_from_k) + sum(bequest_from_k_pps);
        
        dist_v = max(abs(valM_new(:) - valM(:)));
        dist_tr = abs(TR_total_t_new - TR_total_t_eq);
        if iter_main > 1 && dist_v < 1e-5 && dist_tr < 1e-5, break; end

        valM = valM_new;
        TR_total_t_eq = 0.5 * TR_total_t_eq + 0.5 * TR_total_t_new;
    end
    
    % --- 3. 使用最终收敛的策略和遗赠，进行最终模拟 ---
    tr_per_survivor_eq = 0;
    if mass_survived > 0, tr_per_survivor_eq = TR_total_t_eq / mass_survived; end
    
    [cPolM_final, kPolM_final, cPpsPolM_final, ~] = main_olg_v14_utils.HHSolution_VFI_Huggett(M_t, tr_per_survivor_eq, paramS, cS_vfi);
    [kHistM, kPpsHistM, cHistM, cppsHistM] = main_olg_v14_utils.HHSimulation_olgm(kPolM_final, cPpsPolM_final, cPolM_final, eIdxM, M_t, paramS, cS_vfi, tr_per_survivor_eq);

    % --- 4. 聚合 ---
    mean_k_prime_by_age = mean(kHistM(:, 2:end), 1);
    mean_k_pps_prime_by_age = mean(kPpsHistM(:, 2:end), 1);
    
    K_pvt_next = sum(mean_k_prime_by_age' .* Z_t .* cS.s_1yr_transitionV); 
    K_pps_next = sum(mean_k_pps_prime_by_age' .* Z_t .* cS.s_1yr_transitionV);
    
    C_t = sum(mean(cHistM, 1)' .* Z_t); 
    Total_Cpps_t = sum(mean(cppsHistM, 1)' .* Z_t);
    
    prob_death_by_age = (1 - cS.s_1yr_transitionV);
    total_accidental_bequest = sum(mean_k_prime_by_age' .* Z_t .* prob_death_by_age) + ...
                               sum(mean_k_pps_prime_by_age' .* Z_t .* prob_death_by_age);
    
    % --- 5. 精确核算PPS提现税 ---
    Total_PpsTax_t = 0;
    if cS_vfi.pps_active && cS.pps_withdrawal_rate > 0
        mean_k_pps_by_age = mean(kPpsHistM(:, 1:end-1), 1);
        retired_ages_idx = (cS.aR_new + 1):cS.aD_new; 
        if ~isempty(retired_ages_idx)
            k_pps_retirees_avg = mean_k_pps_by_age(retired_ages_idx); 
            mass_retirees_vec = Z_t(retired_ages_idx);
            total_pps_withdrawal_of_retirees = sum(k_pps_retirees_avg' .* mass_retirees_vec) * cS.pps_withdrawal_rate;
            Total_PpsTax_t = total_pps_withdrawal_of_retirees * cS.pps_tax_rate_withdrawal;
        end
    end
    
    % --- 6. 打包返回 ---
    M_t_complete = M_t; 
    M_t_complete.C_t = C_t;
    M_t_complete.Total_Cpps_t = Total_Cpps_t; 
    M_t_complete.Total_PpsTax_t = Total_PpsTax_t;
end        % =========================================================================
        % == 在 main_olg_v14_utils.m 文件中，请使用以下代码块 ==
        % == 替换掉原有的 HHSimulation_olgm 函数。这是最终的正确版本。==
        % =========================================================================

        % =========================================================================
        % == 在 main_olg_v14_utils.m 文件中，请使用以下代码块 ==
        % == 替换掉原有的 HHSimulation_olgm 函数。这是最终的正确版本。==
        % =========================================================================

        % =========================================================================
        % == 在 main_olg_v14_utils.m 文件中，请使用以下代码块 ==
        % == 替换掉原有的 HHSimulation_olgm 函数。这是最终的正确版本。==
        % =========================================================================

        % =========================================================================
        % == 在 main_olg_v14_utils.m 文件中，请使用以下代码块 ==
        % == 替换掉原有的 HHSimulation_olgm 函数。这是最终的正确版本。==
        % =========================================================================
        function [kHistM_out, kPpsHistM_out, cHistM_out, cppsHistM_out] = HHSimulation_olgm(kPolM, cPpsPolM_choice, cPolM_consump, eIdxM_group, M_sim, paramS_sim, cS_sim, tr_per_capita_sim)
            % [最终审计通过版 v10 - 统一劳动税税基]
            % 核心修正: 劳动税(labor_tax)的计算现在也扣除了强制PAYG缴费(payg_tax)，
            %           以确保在模拟中反推消费c时，使用的预算约束与VFI求解器完全对齐。

            nSim = size(eIdxM_group, 1); aD = cS_sim.aD_new;
            kHistM_out = zeros(nSim, aD + 1); kPpsHistM_out = zeros(nSim, aD + 1);
            cHistM_out = zeros(nSim, aD); cppsHistM_out = zeros(nSim, aD);

            % --- 1. 创建策略函数插值器 (不变) ---
            kPolInterp = cell(cS_sim.nw, aD); cPpsPolInterp = cell(cS_sim.nw, aD);
            for ia = 1:aD, for ie = 1:cS_sim.nw
                    if cS_sim.nk > 1 && cS_sim.nkpps > 1, kPolInterp{ie,ia} = griddedInterpolant({cS_sim.kGridV, cS_sim.kppsGridV}, squeeze(kPolM(:,:,ie,ia)), 'linear', 'linear'); cPpsPolInterp{ie,ia} = griddedInterpolant({cS_sim.kGridV, cS_sim.kppsGridV}, squeeze(cPpsPolM_choice(:,:,ie,ia)), 'linear', 'linear');
                    elseif cS_sim.nk > 1, kPolInterp{ie,ia} = griddedInterpolant(cS_sim.kGridV, squeeze(kPolM(:,1,ie,ia)), 'linear', 'linear'); cPpsPolInterp{ie,ia} = griddedInterpolant(cS_sim.kGridV, squeeze(cPpsPolM_choice(:,1,ie,ia)), 'linear', 'linear');
                    else, kPolInterp{ie,ia} = @(x,y) squeeze(kPolM(1,1,ie,ia)); cPpsPolInterp{ie,ia} = @(x,y) squeeze(cPpsPolM_choice(1,1,ie,ia)); end
            end, end

        % --- 2. 按年龄逐期模拟 ---
        for a_idx = 1:aD
            k_now = kHistM_out(:, a_idx); k_pps_now = kPpsHistM_out(:, a_idx);

            % --- 2a. 获取决策 (不变) ---
            k_next_decision = zeros(nSim, 1); cpps_decision = zeros(nSim, 1);
            for ie = 1:cS_sim.nw
                idx_sim = find(eIdxM_group(:, a_idx) == ie); if isempty(idx_sim), continue; end
                k_now_e = k_now(idx_sim); k_pps_now_e = k_pps_now(idx_sim);
                if cS_sim.nk > 1 && cS_sim.nkpps > 1, k_next_decision(idx_sim) = kPolInterp{ie, a_idx}(k_now_e, k_pps_now_e); cpps_decision(idx_sim) = cPpsPolInterp{ie, a_idx}(k_now_e, k_pps_now_e);
                else, k_next_decision(idx_sim) = kPolInterp{ie, a_idx}(k_now_e); cpps_decision(idx_sim) = cPpsPolInterp{ie, a_idx}(k_now_e); end
            end
            cppsHistM_out(:, a_idx) = max(0, cpps_decision);
            kHistM_out(:, a_idx + 1) = max(cS_sim.kMin, k_next_decision);

            % --- 2b. 使用严格对齐的预算约束反推消费 c ---
            tr_this_age = 0; if a_idx > 1, tr_this_age = tr_per_capita_sim; end

            if a_idx <= cS_sim.aR_new % --- 工作期 ---
                labor_income_gross = M_sim.w_t .* cS_sim.ageEffV_new(a_idx) .* paramS_sim.leGridV(eIdxM_group(:, a_idx));

                % 现金流入
                k_return = k_now .* (1 + M_sim.r_net_period);
                total_inflow = k_return + labor_income_gross + tr_this_age;

                % 现金流出 (不含消费)
                payg_tax = cS_sim.theta_t .* labor_income_gross;
                % [核心修正] 劳动税基现在也扣除了PAYG缴费
                labor_tax = cS_sim.tau_l .* max(0, labor_income_gross - cppsHistM_out(:, a_idx) - payg_tax);
                k_prime_outflow = kHistM_out(:, a_idx + 1);
                cpps_outflow = cppsHistM_out(:, a_idx);

                total_outflow_non_c = payg_tax + labor_tax + k_prime_outflow + cpps_outflow;
                c_expend = total_inflow - total_outflow_non_c;

            else % --- 退休期 ---
                % 现金流入
                k_return = k_now .* (1 + M_sim.r_net_period);
                b_age_val = M_sim.b_t;
                pps_withdrawal_gross = 0; if cS_sim.pps_active, pps_withdrawal_gross = k_pps_now .* cS_sim.pps_withdrawal_rate; end
                pps_withdrawal_net = pps_withdrawal_gross .* (1 - cS_sim.pps_tax_rate_withdrawal);
                total_inflow = k_return + b_age_val + tr_this_age + pps_withdrawal_net;

                % 现金流出 (不含消费)
                k_prime_outflow = kHistM_out(:, a_idx + 1);
                c_expend = total_inflow - k_prime_outflow;
            end

            cHistM_out(:, a_idx) = max(cS_sim.cFloor, c_expend ./ (1 + cS_sim.tau_c));

            % --- 2c. 更新下一期PPS资产 (逻辑不变) ---
            pps_return_factor = 1 + M_sim.r_mkt_period + cS_sim.pps_return_rate_premium;
            pps_withdrawal_gross_for_update = 0;
            if a_idx > cS_sim.aR_new && cS_sim.pps_active, pps_withdrawal_gross_for_update = k_pps_now .* cS_sim.pps_withdrawal_rate; end
            k_pps_after_withdrawal = k_pps_now - pps_withdrawal_gross_for_update;
            k_pps_next_unclamped = (k_pps_after_withdrawal + cppsHistM_out(:, a_idx)) * pps_return_factor;
            kPpsHistM_out(:, a_idx + 1) = max(cS_sim.kppsMin, min(cS_sim.kppsMax, k_pps_next_unclamped));
        end
        end

        function [muM, utilM] = CES_utility(cM, sigma, cS)
            c_adj = max(cS.cFloor, cM);
            if abs(sigma - 1) < 1e-6, utilM = log(c_adj); muM = 1./c_adj;
            else, utilM = (c_adj.^(1-sigma))./(1-sigma); muM = c_adj.^(-sigma); end
            utilM(cM < cS.cFloor) = -1e10 - (cS.cFloor - cM(cM < cS.cFloor))*1e10;
        end

        % =====================================================================
        % == 4. 模拟与分析函数 (保留全部) ==
        % =====================================================================

        % --- 劳动禀赋过程 ---
        function [leLogGridV, leTrProbM, leProb1V] = EarningProcess_olgm(cS)
            lePersistence = 0.90;     % 改为 0.90
            leShockStd = 0.15;        % 改为 0.15
            Tauchen_q = 2.0;
            [leLogGridV_raw, leTrProbM] = main_olg_v14_utils.tauchen(cS.nw, lePersistence, leShockStd, 0, Tauchen_q);
            leLogGridV = leLogGridV_raw - mean(leLogGridV_raw);
            [~, D] = eig(leTrProbM');
            [~, c] = min(abs(diag(D)-1));
            leProb1V = abs(D(:,c)/sum(D(:,c)));
        end

        function [y_grid_out, trProbM_out] = tauchen(N, rho, sigma, mu, m)
            % ... (与原版一致)
            std_y = sqrt(sigma^2 / (1-rho^2));
            y_max = m*std_y; y_min = -y_max;
            y = linspace(y_min, y_max, N);
            d = y(2)-y(1);
            trProbM_out = zeros(N,N);
            for j=1:N
                for k=1:N
                    m_k = rho*y(j) + mu;
                    if k==1, trProbM_out(j,k) = normcdf((y(1)-m_k+d/2)/sigma);
                    elseif k==N, trProbM_out(j,k) = 1 - normcdf((y(N)-m_k-d/2)/sigma);
                    else, trProbM_out(j,k) = normcdf((y(k)-m_k+d/2)/sigma) - normcdf((y(k)-m_k-d/2)/sigma); end
                end
            end
            y_grid_out = y(:);
        end

        function eIdxM_group = LaborEndowSimulation_olgm_AgeGroup(cS, paramS)
            eIdxM_group = main_olg_v14_utils.MarkovChainSimulation_AgeGroup(cS.nSim, cS, paramS.leProb1V, paramS.leTrProbM);
        end

        function eIdxM_group_out = MarkovChainSimulation_AgeGroup(num_simulations, cS, p0, P)
            rng(433);
            eIdxM_group_out = zeros(num_simulations, cS.aD_new, 'uint16');
            eIdxM_group_out(:,1) = 1 + sum(rand(num_simulations,1) > cumsum(p0(:)'), 2);
            for a=2:cS.aD_new
                eIdxM_group_out(:,a) = 1 + sum(rand(num_simulations,1) > cumsum(P(eIdxM_group_out(:,a-1),:), 2), 2);
            end
        end

        function ev = CallInterpolator(interpolant, k_prime, k_pps_prime, cS)
            % [新增] 安全的插值器调用函数，处理边界情况
            % 确保查询点在网格范围内，避免外插错误

            try
                % 钳位到网格边界
                k_prime_clamped = max(cS.kGridV(1), min(cS.kGridV(end), k_prime));
                k_pps_prime_clamped = max(cS.kppsGridV(1), min(cS.kppsGridV(end), k_pps_prime));

                % 调用插值器
                ev = interpolant(k_prime_clamped, k_pps_prime_clamped);

                % 处理可能的NaN或Inf
                if ~isfinite(ev)
                    ev = -1e10; % 设为极低值，避免选择这个点
                end
            catch
                % 插值失败时的后备方案
                ev = -1e10;
            end
        end

    end

    methods (Static, Access = private) % <<<<<< 新增一个私有静态方法块

        function neg_v = objective_for_k_prime_private(k_prime_choice, resources, k_pps_state, c_pps_choice, R_k_net_factor_age, a_idx, ev_interpolant, cS)
            % [移到这里] 这个函数现在是类的一个私有静态方法
            % 它的代码内容完全不变

            % 1. 计算消费
            c_expend = resources - k_prime_choice;
            c_choice = max(cS.cFloor, c_expend / (1 + cS.tau_c));

            % 2. 计算当期效用
            [~, util] = main_olg_v14_utils.CES_utility(c_choice, cS.sigma, cS);

            % 3. 计算下一期PPS资产
            pps_withdrawal = 0;
            if a_idx >= cS.aR_new, pps_withdrawal = k_pps_state * cS.pps_withdrawal_rate; end
            pps_return_factor = 1 + ((R_k_net_factor_age - 1) + cS.pps_return_rate_premium);
            k_pps_prime = (k_pps_state + c_pps_choice - pps_withdrawal) * pps_return_factor;

            % 4. [核心修复] 在调用插值器前，对所有查询点进行钳位
            %    fminbnd 优化的 k_prime_choice 已经由其边界约束，但为了安全再次钳位。
            k_prime_clamped = max(cS.kGridV(1), min(cS.kGridV(end), k_prime_choice));
            k_pps_prime_clamped = max(cS.kppsGridV(1), min(cS.kppsGridV(end), k_pps_prime));

            % 使用被钳位的点进行插值
            % 注意：CallInterpolator 函数内部已经有 try-catch 钳位逻辑，
            % 但在这里显式钳位是更稳健的做法，确保我们完全控制了输入。
            ev = main_olg_v14_utils.CallInterpolator(ev_interpolant, k_prime_clamped, k_pps_prime_clamped, cS);

            % 5. 计算总价值
            current_val = util + cS.beta * cS.s_1yr_transitionV(a_idx) * ev;

            % fminbnd是最小化器，所以返回负价值
            neg_v = -current_val;
        end
    end % End of Static Methods
end