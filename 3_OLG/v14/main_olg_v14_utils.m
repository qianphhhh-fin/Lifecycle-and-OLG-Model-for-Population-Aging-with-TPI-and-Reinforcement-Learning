% =========================================================================
% == 在 main_olg_v14_utils.m 文件中，请使用以下代码块 ==
% == 替换或添加相应的函数。                            ==
% =========================================================================

classdef main_olg_v14_utils

    % ... 这里是您从v14复制过来的所有静态方法 ...

    methods (Static)

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
                [~, col_idx] = min(abs(str2double(cellfun(@(x) x(2:end), pop_data.Properties.VariableNames(2:end), 'UniformOutput', false)) - year_t));
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
            % [VFI最终对齐版 v4] - 修正5年期存活率
            cS.pps_active = true;
            % --- 人口结构基础参数 ---
            cS.age1_orig = 20;
            cS.ageLast_orig = 98;
            cS.ageRetire_orig = 65;
            cS.aD_orig = cS.ageLast_orig - cS.age1_orig + 1;
            cS.aR_idx_orig = cS.ageRetire_orig - cS.age1_orig + 1;
            cS.aW_orig = cS.aR_idx_orig - 1;
            cS.physAgeV_orig = (cS.age1_orig : cS.ageLast_orig)';

            % 其他基本参数
            cS.beta = 0.995; % 年化 beta
            cS.sigma = 3; cS.cFloor = 0.05; cS.nSim = 5000;
            cS.sigma = 2.5;      % [非常高] 极度厌恶消费不平滑
            cS.phi_bequest = 0; % [关键] 引入温和但有效的遗赠动机
            cS.sigma_bequest = cS.sigma;
            cS.start_year = 1997; cS.end_year = 2102; cS.time_Step = 5;
            cS.alpha = 0.4; cS.ddk = 1 - (1 - 0.05)^cS.time_Step;
            cS.tau_k = 0.02; cS.tau_l = 0.06; cS.tau_c = 0.03; cS.pps_tax_rate_withdrawal = 0.03;
            cS.gov_exp_frac_Y = 0.15; cS.rho_prime_payg = 0.5; cS.pps_activation_year = 2022;
            cS.pps_return_rate_premium = 0.01;
            cS.pps_withdrawal_rate = 0.15;
            cS.pps_contrib_limit = 9999;
            cS.pps_max_contrib_frac = 0.1;
            cS.A = 1.0;

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

            % 原始年度死亡率数据
            d_orig_data = [0.00159,0.00169,0.00174,0.00172,0.00165,0.00156,0.00149,0.00145,0.00145,0.00149,0.00156,0.00163,0.00171,0.00181,0.00193,0.00207,0.00225,0.00246,0.00270,0.00299,0.00332,0.00368,0.00409,0.00454,0.00504,0.00558,0.00617,0.00686,0.00766,0.00865,0.00955,0.01058,0.01162,0.01264,0.01368,0.01475,0.01593,0.01730,0.01891,0.02074,0.02271,0.02476,0.02690,0.02912,0.03143,0.03389,0.03652,0.03930,0.04225,0.04538,0.04871,0.05230,0.05623,0.06060,0.06542,0.07066,0.07636,0.08271,0.08986,0.09788,0.10732,0.11799,0.12895,0.13920,0.14861,0.16039,0.17303,0.18665,0.20194,0.21877,0.23601,0.25289,0.26973,0.28612,0.30128,0.31416,0.32915,0.34450,0.36018];
            s_orig_1yr = 1 - d_orig_data; % 年度存活率

            % [核心修正] 计算5年期的存活概率
            cS.s_pathV = zeros(cS.aD_new, 1);
            for a = 1:(cS.aD_new - 1)
                % 从年龄组a存活到年龄组a+1的概率，是组内5年连续存活的概率
                age_indices = cS.physAgeMap{a};
                s_5yr = 1.0;
                for i = 1:length(age_indices)
                    s_5yr = s_5yr * s_orig_1yr(age_indices(i));
                end
                cS.s_pathV(a) = s_5yr;
            end
            cS.s_pathV(cS.aD_new) = 0; % 最后一个年龄组存活率为0

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

            % --- 政府债务参数 ---
            cS.gov_debt_frac_Y = 0.3; % 假设政府债务占GDP的30%
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

        function M_prices = get_prices_at_t(K, L, A_t, cS)

            % [v14.9 稳健版] 为利率增加上限以提高求解器稳定性
            if K <= 0
                K = 1e-8;
            end
            if L <= 0
                L = 1e-8;
            end

            % [核心修正] 使用按元素求幂 (.^) 和按元素除法 (./)
            Y_period = A_t .* (K.^cS.alpha) .* (L.^(1-cS.alpha));
            MPK_period = cS.alpha .* Y_period ./ K;
            w_t = (1-cS.alpha) .* Y_period ./ L;

            r_mkt_t = MPK_period - cS.ddk;
            r_net_t = r_mkt_t .* (1 - cS.tau_k);

            M_prices = struct();
            M_prices.K = K;
            M_prices.L_t = L;
            M_prices.Y_t = Y_period; % 使用已经计算好的 Y_period
            M_prices.w_t = w_t;
            M_prices.r_mkt_t = r_mkt_t;
            M_prices.r_net_period = r_net_t;
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

        function [ss, eq_found, initial_dist_k, initial_dist_kpps] = solve_steady_state_with_fund(Z_ss_norm, B_p_Y_ratio_target, cS, paramS, eIdxM)
            % [v14.9 - 逻辑重构版]
            fprintf('\n--- 开始求解初始稳态 (有目标基金 B_p, 无初始PPS) ---\n');
            initial_dist_k = []; initial_dist_kpps = [];
            K_guess = 3.0; % 从一个更合理的猜测开始，比如3.0
            max_iter = 500; tol = 1e-4; damp = 0.5; % 阻尼系数可以适当调高
            eq_found = false; A_ss = cS.A; theta_ss = cS.theta_path(1);
            paramS_ss_for_L = paramS; paramS_ss_for_L.ageMassV = Z_ss_norm;
            [~, L_ss] = main_olg_v14_utils.LaborSupply_Huggett(eIdxM, cS, paramS_ss_for_L, Z_ss_norm);

            fprintf('%5s | %12s | %12s | %12s | %12s | %12s | %12s | %12s\n', 'Iter', 'K_guess', 'K_model', 'K_pvt', 'K_pps', 'B_p', 'r_net', 'K_error');
            fprintf('%s\n', repmat('-', 105, 1));

            for iter = 1:max_iter
                % 1. 根据 K_guess 计算正确的价格
                M_ss_prices = main_olg_v14_utils.get_prices_at_t(K_guess, L_ss, A_ss, cS);
                B_p_target = B_p_Y_ratio_target * M_ss_prices.Y_t;

                % 2. 准备参数，传递给家庭求解器
                M_ss = M_ss_prices;
                M_ss.current_t = 1; % 设定时期
                mass_workers_ss = sum(Z_ss_norm(1:cS.aR_new));
                avg_wage_ss = M_ss.w_t * M_ss.L_t / mass_workers_ss;
                M_ss.b_t = cS.rho_prime_payg * avg_wage_ss;

                cS_ss = cS;
                cS_ss.pps_active = false; % 初始稳态无PPS
                cS_ss.theta_t = theta_ss;
                if ~cS_ss.pps_active
                    cS_ss.nkpps = 1; cS_ss.npps = 1; cS_ss = main_olg_v14_utils.generateGrids(cS_ss);
                end

                % 3. [核心修改] 调用新函数，在给定价格下聚合储蓄
                [K_pvt_model, K_pps_model, ~] = main_olg_v14_utils.aggregate_household_savings(M_ss, Z_ss_norm, cS_ss, paramS, eIdxM);

                % 4. 计算模型总资本并更新猜测
                K_total_model = K_pvt_model + K_pps_model + B_p_target; % 政府债务 B_g 是融资项，不计入资本存量
                K_error = K_guess - K_total_model;

                fprintf('%5d | %12.4f | %12.4f | %12.4f | %12.4f | %12.4f | %11.2f%% | %12.3e\n', iter, K_guess, K_total_model, K_pvt_model, K_pps_model, B_p_target, M_ss_prices.r_net_period*100, K_error);
                if abs(K_error) < tol, fprintf('✅ 初始稳态均衡收敛！\n'); eq_found = true; break; end
                K_guess = (1 - damp) * K_guess + damp * K_total_model;
                if K_guess < 0, K_guess = 0.1; end % 防止 K_guess 变成负数
            end

            if ~eq_found, warning('初始稳态均衡未在最大迭代次数内收敛。'); ss=struct(); return; end

            % --- 收敛后，重新计算所有量以确保最终一致性 ---
            fprintf('--- 基于最终均衡资本 K=%.4f 重新计算所有量以确保自洽性 ---\n', K_total_model);
            K_final = K_total_model;
            M_final_prices = main_olg_v14_utils.get_prices_at_t(K_final, L_ss, A_ss, cS);
            M_final = M_final_prices; M_final.current_t = 1;
            mass_workers_final = sum(Z_ss_norm(1:cS.aR_new));
            avg_wage_final = M_final.w_t * M_final.L_t / mass_workers_final;
            M_final.b_t = cS.rho_prime_payg * avg_wage_final;

            cS_final = cS; cS_final.pps_active = false; cS_final.theta_t = theta_ss;
            if ~cS_final.pps_active
                cS_final.nkpps = 1; cS_final.npps = 1; cS_final = main_olg_v14_utils.generateGrids(cS_final);
            end

            % 获取最终的策略和分布
            tr_per_survivor_eq = 0;
            [cPolM, kPolM, cPpsPolM, ~] = main_olg_v14_utils.HHSolution_VFI_Huggett(M_final, tr_per_survivor_eq, paramS, cS_final);
            [kHistM, kPpsHistM, cHistM, ~] = main_olg_v14_utils.HHSimulation_olgm(kPolM, cPpsPolM, cPolM, eIdxM, M_final, paramS, cS_final, tr_per_survivor_eq);
            initial_dist_k = kHistM;
            initial_dist_kpps = kPpsHistM;

            % 使用最终分布聚合
            [K_pvt_final, K_pps_final, C_final] = main_olg_v14_utils.aggregate_household_savings(M_final, Z_ss_norm, cS_final, paramS, eIdxM);
            B_p_final = B_p_Y_ratio_target * M_final_prices.Y_t;
            B_g_final = cS.gov_debt_frac_Y * M_final_prices.Y_t;

            ss = struct();
            ss.K_total = K_final; ss.K_pvt = K_pvt_final; ss.K_pps = K_pps_final;
            ss.B_p = B_p_final; ss.B_g = B_g_final; ss.L = M_final_prices.L_t;
            ss.Y = M_final_prices.Y_t; ss.w = M_final_prices.w_t; ss.r_mkt = M_final_prices.r_mkt_t;
            ss.r_net = M_final_prices.r_net_period; ss.C = C_final; ss.G = cS.gov_exp_frac_Y * ss.Y;
        end

        function [K_pvt_agg, K_pps_agg, C_agg] = aggregate_household_savings(M_prices_given, Z_dist, cS_agg, paramS, eIdxM)
            % [新函数]
            % 核心目的：给定一套宏观价格，求解家庭问题，并返回聚合后的储蓄和消费。
            % 这个函数【不】自己计算价格，而是严格使用给定的价格。

            % 1. 求解家庭最优策略
            %    假设意外遗赠为0，这是稳态下的标准做法。
            tr_per_survivor_eq = 0;
            [cPolM, kPolM, cPpsPolM, ~] = main_olg_v14_utils.HHSolution_VFI_Huggett(M_prices_given, tr_per_survivor_eq, paramS, cS_agg);

            % 2. 模拟家庭生命周期决策
            %    在稳态下，家庭的初始资产分布是内生的，所以我们从0开始模拟。
            %    注意：HHSimulation_olgm 返回的是一个历史矩阵，其中 kHistM(:,a) 是年龄a期初的资产。
            %    我们需要的是他们选择的下一期资产 k'，即 kHistM(:, a+1)。
            [kHistM, kPpsHistM, cHistM, ~] = main_olg_v14_utils.HHSimulation_olgm(kPolM, cPpsPolM, cPolM, eIdxM, M_prices_given, paramS, cS_agg, tr_per_survivor_eq);

            % 3. 聚合下一期的资本存量和当期消费
            mean_k_prime_by_age     = mean(kHistM(:, 2:end), 1); % k'
            mean_k_pps_prime_by_age = mean(kPpsHistM(:, 2:end), 1); % k_pps'

            % 稳态下的聚合：下一期的总资本 = 各年龄组选择的储蓄 * 各年龄组人口占比 * 存活到下一期的概率
            % 注意：在稳态下，下一期的分布 Z_{t+1} 和人口结构 Z_t 是一样的。
            % cS.s_pathV(a) 是从年龄组 a 存活到 a+1 的概率
            K_pvt_agg = sum(mean_k_prime_by_age' .* Z_dist .* cS_agg.s_pathV);
            K_pps_agg = sum(mean_k_pps_prime_by_age' .* Z_dist .* cS_agg.s_pathV);

            % 当期消费的聚合
            C_agg = sum(mean(cHistM, 1)' .* Z_dist);
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


        function [B_p_next, PensionSurplus_period] = update_pension_fund(B_p_t, M_t, Z_t_norm, cS)
            theta_t = cS.theta_path(M_t.current_t);
            PensionRevenue_period = theta_t * M_t.w_t * M_t.L_t;
            mass_retirees_t = sum(Z_t_norm(cS.aR_new+1:end));
            PensionOutlay_period = M_t.b_t * mass_retirees_t;
            PensionSurplus_period = PensionRevenue_period - PensionOutlay_period;
            B_p_next = B_p_t * (1 + M_t.r_mkt_t) + PensionSurplus_period;
            B_p_next = max(0, B_p_next);
        end

        function [B_g_next, G_t, TotalTax_t] = update_gov_debt(B_g_t, M_t_complete, C_t, cS, PensionSurplus_period, TotalBequestRevenue_t)
            % [版本：遗赠归政府]
            % 核心修正: 增加了一个新的输入参数 TotalBequestRevenue_t，
            %           并将其计入政府的基本财政盈余(Primary Surplus)。

            % --- 1. 从 M_t 结构体中获取所有需要的当期变量 (不变) ---
            w_t              = M_t_complete.w_t;
            L_t              = M_t_complete.L_t;
            Y_t              = M_t_complete.Y_t;
            r_mkt_t          = M_t_complete.r_mkt_t;
            K_pvt_t          = M_t_complete.K_pvt_t;
            Total_Cpps_t     = M_t_complete.Total_Cpps_t;
            Total_PpsTax_t   = M_t_complete.Total_PpsTax_t;
            theta_t          = cS.theta_path(M_t_complete.current_t);

            % --- 2. 计算各项税收 (不变) ---
            payg_contributions_total = theta_t * w_t * L_t;
            taxable_labor_income = max(0, w_t * L_t - payg_contributions_total - Total_Cpps_t);
            LaborTax_t = cS.tau_l * taxable_labor_income;
            CapitalTax_t = cS.tau_k * r_mkt_t * K_pvt_t;
            ConsumptionTax_t = cS.tau_c * C_t;
            PpsWithdrawalTax_t = Total_PpsTax_t;
            TotalTax_t = LaborTax_t + CapitalTax_t + ConsumptionTax_t + PpsWithdrawalTax_t;

            % --- 3. 计算政府支出与主预算盈余 ---
            G_t = cS.gov_exp_frac_Y * Y_t;
            pension_deficit_subsidy = max(0, -PensionSurplus_period);

            % [核心修正] 将意外遗赠加入政府基本盈余
            PrimarySurplus_g = (TotalTax_t + TotalBequestRevenue_t) - G_t - pension_deficit_subsidy;

            % --- 4. 更新政府债务 (不变) ---
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
            % 处理养老金收益的维度匹配问题
            if cS_vfi.aR_new < cS_vfi.aD_new
                retirement_indices = (cS_vfi.aR_new+1):cS_vfi.aD_new;
                num_retirement_periods = length(retirement_indices);

                if isscalar(M_vfi.b_t)
                    % 如果 b_t 是标量，则复制到所有退休期
                    bV_payg_vfi(retirement_indices) = M_vfi.b_t;
                else
                    % 如果 b_t 是向量，确保维度匹配
                    b_t_vec = M_vfi.b_t(:)';  % 转换为行向量
                    if length(b_t_vec) == num_retirement_periods
                        bV_payg_vfi(retirement_indices) = b_t_vec;
                    else
                        % 如果长度不匹配，使用第一个值填充所有退休期
                        bV_payg_vfi(retirement_indices) = b_t_vec(1);
                    end
                end
            end

            for a_idx = cS_vfi.aD_new : -1 : 1
                vPrime_kkppse_next = [];
                if a_idx < cS_vfi.aD_new
                    vPrime_kkppse_next = valM(:,:,:,a_idx+1);
                end

                [cPolM(:,:,:,a_idx), kPolM(:,:,:,a_idx), cPpsPolM(:,:,:,a_idx), valM(:,:,:,a_idx)] = ...
                    main_olg_v14_utils.HHSolutionByAge_VFI_Vectorized(a_idx, vPrime_kkppse_next, M_vfi, tr_per_capita_vfi, bV_payg_vfi(a_idx), paramS_vfi, cS_vfi);
                % [cPolM(:,:,:,a_idx), kPolM(:,:,:,a_idx), cPpsPolM(:,:,:,a_idx), valM(:,:,:,a_idx)] = ...
                %     main_olg_v14_utils.HHSolutionByAge_NAIVE(a_idx,  M_vfi, tr_per_capita_vfi, bV_payg_vfi(a_idx), paramS_vfi, cS_vfi);
            end
        end

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
            % [v19.1 - VFI求解器最终修复版]
            % 核心修复：在调用EV插值器之前，对 k_pps_prime 进行边界约束(clamping)，
            %           防止因查询点越界导致EV返回NaN，进而使储蓄决策崩溃至0。

            val_age = -1e20 * ones(cS.nk, cS.nkpps, cS.nw);
            cPol_age_q = zeros(cS.nk, cS.nkpps, cS.nw);
            kPol_age = zeros(cS.nk, cS.nkpps, cS.nw);
            cPpsPol_age_choice = zeros(cS.nk, cS.nkpps, cS.nw);

            % --- Last Age Period (Terminal Condition) ---
            if a_idx == cS.aD_new
                [K_grid, Kpps_grid] = ndgrid(cS.kGridV, cS.kppsGridV);
                [pretax_non_capital_income, ~] = main_olg_v14_utils.HHIncome_Huggett(0, tr_per_capita_age, b_age_val, 0, a_idx, paramS_age, cS, 0);
                k_end_value = K_grid .* (1 + M_age.r_net_period);
                pps_liquidation_net = Kpps_grid .* (1 + M_age.r_mkt_t + cS.pps_return_rate_premium) * (1 - cS.pps_tax_rate_withdrawal);
                total_resources = k_end_value + pps_liquidation_net + pretax_non_capital_income;
                final_bequest = zeros(size(total_resources));
                util_b = 0;
                if isfield(cS, 'phi_bequest') && cS.phi_bequest > 1e-9
                    phi_adj = cS.phi_bequest * (1 + cS.tau_c)^(-1);
                    omega = (1 + phi_adj^(-1/cS.sigma))^(-1);
                    final_c_expenditure = omega .* total_resources;
                    final_bequest = (1 - omega) .* total_resources;
                    [~, util_b] = main_olg_v14_utils.CES_utility(final_bequest, cS.sigma_bequest, cS);
                else
                    final_c_expenditure = total_resources;
                end
                final_c = max(cS.cFloor, final_c_expenditure / (1 + cS.tau_c));
                [~, util_c] = main_olg_v14_utils.CES_utility(final_c, cS.sigma, cS);
                final_v = util_c;
                if isfield(cS, 'phi_bequest'), final_v = final_v + cS.phi_bequest * util_b; end
                for ie = 1:cS.nw
                    cPol_age_q(:,:,ie) = final_c;
                    val_age(:,:,ie) = final_v;
                    kPol_age(:,:,ie) = final_bequest;
                    cPpsPol_age_choice(:,:,ie) = 0;
                end
                return;
            end

            % --- Standard Age Periods (Backward Induction) ---
            effective_discount_factor = (cS.beta ^ cS.time_Step) * cS.s_pathV(a_idx);

            % --- Prepare EV interpolants ---
            EV_matrix = zeros(cS.nk, cS.nkpps, cS.nw);
            for ie_current = 1:cS.nw
                transition_probs = paramS_age.leTrProbM(ie_current, :);
                if isempty(vPrime_kkppse_next), EV_matrix = zeros(cS.nk, cS.nkpps, cS.nw); break; end
                vPrime_reshaped = reshape(vPrime_kkppse_next, [cS.nk * cS.nkpps, cS.nw]);
                EV_slice = vPrime_reshaped * transition_probs';
                EV_matrix(:, :, ie_current) = reshape(EV_slice, [cS.nk, cS.nkpps]);
            end
            EV_interpolants = cell(cS.nw, 1);
            is_pps_disabled = (cS.nkpps == 1);
            for ie_current = 1:cS.nw
                if is_pps_disabled
                    EV_interpolants{ie_current} = griddedInterpolant(cS.kGridV, squeeze(EV_matrix(:, 1, ie_current)), 'linear', 'none');
                else
                    EV_interpolants{ie_current} = griddedInterpolant({cS.kGridV, cS.kppsGridV}, EV_matrix(:, :, ie_current), 'linear', 'none');
                end
            end

            % --- Loop over current states (k, k_pps, e) ---
            for ik = 1:cS.nk, for ikpps = 1:cS.nkpps, for ie = 1:cS.nw
                        k_state = cS.kGridV(ik); k_pps_state = cS.kppsGridV(ikpps); epsilon_state = paramS_age.leGridV(ie);
                        ev_interpolant = EV_interpolants{ie};
                        best_val = -1e20; best_c = cS.cFloor; best_k_prime = cS.kMin; best_c_pps = 0;
                        tr_this_age = 0; if a_idx > 1, tr_this_age = tr_per_capita_age; end

                        if a_idx <= cS.aR_new % --- Working Age ---
                            labor_income_gross = M_age.w_t * cS.ageEffV_new(a_idx) * epsilon_state;
                            max_cpps = 0; if cS.pps_active, max_cpps = min(cS.pps_contrib_limit, labor_income_gross * cS.pps_max_contrib_frac); end
                            cpps_grid = linspace(0, max(0, max_cpps), cS.npps);
                            payg_tax = cS.theta_t * labor_income_gross;

                            for c_pps_choice = cpps_grid
                                k_return = k_state * (1 + M_age.r_net_period);
                                total_inflow = k_return + labor_income_gross + tr_this_age;
                                labor_tax = cS.tau_l * max(0, labor_income_gross - c_pps_choice - payg_tax);
                                cpps_outflow = c_pps_choice;
                                net_cash_for_c_and_k_prime = total_inflow - (payg_tax + labor_tax + cpps_outflow);
                                k_prime_max = net_cash_for_c_and_k_prime - cS.cFloor * (1 + cS.tau_c);

                                if k_prime_max < cS.kMin, k_prime_grid = [cS.kMin]; else
                                    k_prime_grid_linear = linspace(cS.kMin, k_prime_max, cS.nkprime);
                                    k_prime_grid_global = cS.kGridV(cS.kGridV <= k_prime_max)';
                                    k_prime_grid = unique([k_prime_grid_linear, k_prime_grid_global]);
                                end

                                for k_prime_choice = k_prime_grid
                                    c_expend = net_cash_for_c_and_k_prime - k_prime_choice;
                                    if c_expend < cS.cFloor * (1 + cS.tau_c), continue; end
                                    c_choice = c_expend / (1 + cS.tau_c);
                                    [~, util] = main_olg_v14_utils.CES_utility(c_choice, cS.sigma, cS);

                                    pps_return_factor = 1 + M_age.r_mkt_t + cS.pps_return_rate_premium;
                                    k_pps_prime = (k_pps_state + c_pps_choice) * pps_return_factor;

                                    % !!!!! 核心修复 !!!!!
                                    % 在调用插值器前，确保 k_pps_prime 在网格范围内
                                    k_pps_prime = max(cS.kppsGridV(1), min(cS.kppsGridV(end), k_pps_prime));

                                    if is_pps_disabled, ev = ev_interpolant(k_prime_choice); else, ev = ev_interpolant(k_prime_choice, k_pps_prime); end
                                    if isnan(ev), ev = -1e10; end % 防御性代码

                                    current_val = util + effective_discount_factor * ev;
                                    if current_val > best_val
                                        best_val = current_val;
                                        best_c = c_choice;
                                        best_k_prime = k_prime_choice;
                                        best_c_pps = c_pps_choice;
                                    end
                                end
                            end
                        else % --- Retirement Age ---
                            c_pps_choice = 0;
                            k_return = k_state * (1 + M_age.r_net_period);
                            pps_withdrawal_gross = 0; if cS.pps_active, pps_withdrawal_gross = k_pps_state * cS.pps_withdrawal_rate; end
                            pps_withdrawal_net = pps_withdrawal_gross * (1 - cS.pps_tax_rate_withdrawal);
                            total_inflow = k_return + b_age_val + tr_this_age + pps_withdrawal_net;
                            net_cash = total_inflow;
                            k_prime_max = net_cash - cS.cFloor * (1 + cS.tau_c);

                            if k_prime_max < cS.kMin, k_prime_grid = [cS.kMin]; else
                                k_prime_grid_linear = linspace(cS.kMin, k_prime_max, cS.nkprime);
                                k_prime_grid_global = cS.kGridV(cS.kGridV <= k_prime_max)';
                                k_prime_grid = unique([k_prime_grid_linear, k_prime_grid_global]);
                            end

                            for k_prime_choice = k_prime_grid
                                if k_prime_choice < cS.kMin, continue; end
                                c_expend = net_cash - k_prime_choice;
                                if c_expend < cS.cFloor * (1 + cS.tau_c), continue; end
                                c_choice = c_expend / (1 + cS.tau_c);
                                [~, util] = main_olg_v14_utils.CES_utility(c_choice, cS.sigma, cS);

                                pps_return_factor = 1 + M_age.r_mkt_t + cS.pps_return_rate_premium;
                                k_pps_prime = (k_pps_state - pps_withdrawal_gross) * pps_return_factor;

                                % !!!!! 核心修复 !!!!!
                                % 在调用插值器前，确保 k_pps_prime 在网格范围内
                                k_pps_prime = max(cS.kppsGridV(1), min(cS.kppsGridV(end), k_pps_prime));

                                if is_pps_disabled, ev = ev_interpolant(k_prime_choice); else, ev = ev_interpolant(k_prime_choice, k_pps_prime); end
                                if isnan(ev), ev = -1e10; end % 防御性代码

                                current_val = util + effective_discount_factor * ev;
                                if current_val > best_val, best_val = current_val; best_c = c_choice; best_k_prime = k_prime_choice; best_c_pps = c_pps_choice; end
                            end
                        end

                        % --- Store optimal choices ---
                        if isinf(best_val) || isnan(best_val)
                            val_age(ik, ikpps, ie) = -1e20; cPol_age_q(ik, ikpps, ie) = cS.cFloor; kPol_age(ik, ikpps, ie) = cS.kMin; cPpsPol_age_choice(ik, ikpps, ie) = 0;
                        else
                            val_age(ik, ikpps, ie) = best_val; cPol_age_q(ik, ikpps, ie) = best_c; kPol_age(ik, ikpps, ie) = best_k_prime; cPpsPol_age_choice(ik, ikpps, ie) = best_c_pps;
                        end
            end, end, end
        end

        % =========================================================================
        % == 在 main_olg_v14_utils.m 中，请【替换】此函数
        % =========================================================================
        % =========================================================================
        % == 在 main_olg_v14_utils.m 中，请【替换】此函数
        % =========================================================================
        % =========================================================================
        % == 在 main_olg_v14_utils.m 中，请【替换】此函数
        % =========================================================================
        function [cPol_age_q, kPol_age, cPpsPol_age_choice, val_age] = HHSolutionByAge_VFI_Vectorized(a_idx, vPrime_kkppse_next, M_age, tr_per_capita_age, b_age_val, paramS_age, cS)
            % =========================================================================
            % == 函数: HHSolutionByAge_VFI_Vectorized [v6.3 - 遗赠动机最终版]
            % == 核心修正:
            % == 1. 彻底修正了最后生命周期的遗赠动机实现。使用了基于严格的
            % ==    一阶条件推导出的最优消费/遗赠分配公式。
            % == 2. 此修正解决了储蓄对利率反应不正常的根本问题，是 fzero
            % ==    得以正常工作的关键。
            % =========================================================================

            % --- 0. 初始化 (不变) ---
            val_age = -1e20 * ones(cS.nk, cS.nkpps, cS.nw);
            cPol_age_q = zeros(cS.nk, cS.nkpps, cS.nw);
            kPol_age = zeros(cS.nk, cS.nkpps, cS.nw);
            cPpsPol_age_choice = zeros(cS.nk, cS.nkpps, cS.nw);

            % --- 1. 最后生命周期 (Terminal Condition) ---
            if a_idx == cS.aD_new
                % [核心修正] 使用正确的消费/遗赠分配公式
                [K_grid, ~] = ndgrid(cS.kGridV, cS.kppsGridV);

                tr_this_age = tr_per_capita_age;
                pretax_non_capital_income = b_age_val + tr_this_age;

                capital_tax = cS.tau_k .* M_age.r_mkt_t .* K_grid;
                total_resources = K_grid .* (1 + M_age.r_mkt_t) + pretax_non_capital_income - capital_tax;

                final_bequest = zeros(size(total_resources));
                final_v = -1e20 * ones(size(total_resources));

                if isfield(cS, 'phi_bequest') && cS.phi_bequest > 1e-9 && cS.sigma_bequest == cS.sigma
                    % 基于 FOC: u'(c)/(1+tau_c) = phi * u'(b) => c/b = (phi*(1+tau_c))^(-1/sigma)
                    psi = (cS.phi_bequest * (1 + cS.tau_c))^(-1/cS.sigma);

                    % 预算约束: c*(1+tau_c) + b = R
                    % 联立求解得到消费支出占总资源的比例 omega
                    omega = (psi * (1+cS.tau_c)) / (psi * (1+cS.tau_c) + 1);

                    final_c_expenditure = omega .* total_resources;
                    final_bequest = total_resources - final_c_expenditure;

                    final_c = max(cS.cFloor, final_c_expenditure / (1 + cS.tau_c));

                    [~, util_c] = main_olg_v14_utils.CES_utility(final_c, cS.sigma, cS);
                    [~, util_b] = main_olg_v14_utils.CES_utility(final_bequest, cS.sigma_bequest, cS);

                    final_v = util_c + cS.phi_bequest * util_b;
                else
                    % 如果没有遗赠动机，或效用函数形式不匹配，则全部消费
                    final_c_expenditure = total_resources;
                    final_bequest = zeros(size(total_resources));
                    final_c = max(cS.cFloor, final_c_expenditure / (1 + cS.tau_c));
                    [~, util_c] = main_olg_v14_utils.CES_utility(final_c, cS.sigma, cS);
                    final_v = util_c;
                end

                for ie = 1:cS.nw
                    cPol_age_q(:,:,ie) = final_c;
                    val_age(:,:,ie) = final_v;
                    kPol_age(:,:,ie) = final_bequest;
                    cPpsPol_age_choice(:,:,ie) = 0;
                end
                return;
            end

            % --- 2. 准备期望价值插值器 (不变) ---
            effective_discount_factor = (cS.beta ^ cS.time_Step) * cS.s_pathV(a_idx);
            EV_matrix = zeros(cS.nk, cS.nkpps, cS.nw);
            for ie_current = 1:cS.nw
                transition_probs = paramS_age.leTrProbM(ie_current, :);
                if isempty(vPrime_kkppse_next), EV_matrix = zeros(cS.nk, cS.nkpps, cS.nw); break; end
                vPrime_reshaped = reshape(vPrime_kkppse_next, [cS.nk * cS.nkpps, cS.nw]);
                EV_slice = vPrime_reshaped * transition_probs';
                EV_matrix(:, :, ie_current) = reshape(EV_slice, [cS.nk, cS.nkpps]);
            end
            EV_interpolants = cell(cS.nw, 1);
            is_pps_disabled = (cS.nkpps == 1);
            for ie_current=1:cS.nw
                if is_pps_disabled, EV_interpolants{ie_current}=griddedInterpolant(cS.kGridV,squeeze(EV_matrix(:,1,ie_current)),'linear','none');
                else, EV_interpolants{ie_current}=griddedInterpolant({cS.kGridV, cS.kppsGridV},EV_matrix(:,:,ie_current),'linear','none'); end
            end

            % --- 3. 逐状态求解主循环 (不变) ---
            for ie = 1:cS.nw
                epsilon_state = paramS_age.leGridV(ie);
                ev_interpolant = EV_interpolants{ie};
                tr_this_age = tr_per_capita_age;

                for ik = 1:cS.nk, for ikpps = 1:cS.nkpps
                        k_state = cS.kGridV(ik); k_pps_state = cS.kppsGridV(ikpps);
                        best_val = -1e20; best_k_prime = cS.kMin; best_c_pps = 0; best_c = cS.cFloor;

                        if a_idx <= cS.aR_new % --- 工作期 ---
                            labor_income_gross = M_age.w_t*cS.ageEffV_new(a_idx)*epsilon_state;
                            max_cpps = 0; if cS.pps_active, max_cpps=min(cS.pps_contrib_limit,labor_income_gross*cS.pps_max_contrib_frac); end
                            cpps_grid_choices=linspace(0,max(0,max_cpps),cS.npps);
                            payg_tax=cS.theta_t*labor_income_gross;
                            cash_inflow_gross=k_state*(1+M_age.r_mkt_t)+labor_income_gross+tr_this_age;
                            capital_tax=cS.tau_k*M_age.r_mkt_t*k_state;
                            for c_pps_choice = cpps_grid_choices
                                labor_tax=cS.tau_l*max(0,labor_income_gross-c_pps_choice-payg_tax);
                                net_cash_for_c_k_prime=cash_inflow_gross-(payg_tax+labor_tax+capital_tax+c_pps_choice);
                                k_prime_max=net_cash_for_c_k_prime-cS.cFloor*(1+cS.tau_c);
                                if k_prime_max<cS.kMin, k_prime_grid=[cS.kMin]; else, k_prime_grid=unique([linspace(cS.kMin,k_prime_max,cS.nkprime), cS.kGridV(cS.kGridV<=k_prime_max)']); end
                                for k_prime_choice=k_prime_grid
                                    c_expend=net_cash_for_c_k_prime-k_prime_choice; if c_expend<cS.cFloor*(1+cS.tau_c), continue; end
                                    c_choice=c_expend/(1+cS.tau_c);
                                    [~,util]=main_olg_v14_utils.CES_utility(c_choice,cS.sigma,cS);
                                    k_pps_prime=(k_pps_state+c_pps_choice)*(1+M_age.r_mkt_t+cS.pps_return_rate_premium);
                                    k_pps_prime=max(cS.kppsGridV(1),min(cS.kppsGridV(end),k_pps_prime));
                                    if is_pps_disabled, ev=ev_interpolant(k_prime_choice); else, ev=ev_interpolant(k_prime_choice,k_pps_prime); end
                                    if isnan(ev), ev=-1e10; end
                                    current_val=util+effective_discount_factor*ev;
                                    if current_val>best_val, best_val=current_val; best_c=c_choice; best_k_prime=k_prime_choice; best_c_pps=c_pps_choice; end
                                end
                            end
                        else % --- 退休期 ---
                            c_pps_choice=0; b_val_scalar=b_age_val; if ~isscalar(b_age_val), b_val_scalar=b_age_val(1); end
                            cash_inflow_gross=k_state*(1+M_age.r_mkt_t)+b_val_scalar+tr_this_age;
                            capital_tax=cS.tau_k*M_age.r_mkt_t*k_state;
                            pps_withdrawal_gross=0; if cS.pps_active, pps_withdrawal_gross=k_pps_state*cS.pps_withdrawal_rate; end
                            pps_withdrawal_net=pps_withdrawal_gross*(1-cS.pps_tax_rate_withdrawal);
                            net_cash_for_c_k_prime=cash_inflow_gross-capital_tax+pps_withdrawal_net;
                            k_prime_max=net_cash_for_c_k_prime-cS.cFloor*(1+cS.tau_c);
                            if k_prime_max<cS.kMin, k_prime_grid=[cS.kMin]; else, k_prime_grid=unique([linspace(cS.kMin,k_prime_max,cS.nkprime), cS.kGridV(cS.kGridV<=k_prime_max)']); end
                            for k_prime_choice=k_prime_grid
                                c_expend=net_cash_for_c_k_prime-k_prime_choice; if c_expend<cS.cFloor*(1+cS.tau_c), continue; end
                                c_choice=c_expend/(1+cS.tau_c);
                                [~,util]=main_olg_v14_utils.CES_utility(c_choice,cS.sigma,cS);
                                k_pps_prime=(k_pps_state-pps_withdrawal_gross)*(1+M_age.r_mkt_t+cS.pps_return_rate_premium);
                                k_pps_prime=max(cS.kppsGridV(1),min(cS.kppsGridV(end),k_pps_prime));
                                if is_pps_disabled, ev=ev_interpolant(k_prime_choice); else, ev=ev_interpolant(k_prime_choice,k_pps_prime); end
                                if isnan(ev), ev=-1e10; end
                                current_val=util+effective_discount_factor*ev;
                                if current_val>best_val, best_val=current_val; best_c=c_choice; best_k_prime=k_prime_choice; end
                            end
                        end
                        if isinf(best_val)||isnan(best_val), val_age(ik,ikpps,ie)=-1e20; kPol_age(ik,ikpps,ie)=cS.kMin; cPpsPol_age_choice(ik,ikpps,ie)=0; cPol_age_q(ik,ikpps,ie)=cS.cFloor;
                        else, val_age(ik,ikpps,ie)=best_val; kPol_age(ik,ikpps,ie)=best_k_prime; cPpsPol_age_choice(ik,ikpps,ie)=best_c_pps; cPol_age_q(ik,ikpps,ie)=best_c; end
                end, end
            end
        end
        function PolicyFunctions = solve_transition_path_policies(w_path, r_mkt_path, b_path, T_sim, cS, paramS)
            % [v14.1 - 转型路径求解器修复版]
            % 核心功能：给定完整的价格路径，从期末T反向迭代，为每个时期t和每个年龄a
            %           计算出最优的决策函数(策略函数)和价值函数。

            fprintf('--- [转型路径] 开始反向求解所有策略函数 (t=%d to 1) ---\n', T_sim);

            PolicyFunctions = cell(T_sim, cS.aD_new);
            V_future_by_age = cell(cS.aD_new, 1);

            % --- 从最后一期 T 开始，反向时间循环 ---
            for t = T_sim:-1:1
                fprintf('   正在处理时期 t = %d\n', t)

                % --- 1. 准备当期(t)的环境 ---
                M_t = struct();
                M_t.w_t = w_path(t);
                M_t.r_mkt_t = r_mkt_path(t);
                M_t.r_net_period = r_mkt_path(t) * (1 - cS.tau_k);
                M_t.b_t = b_path(t);

                cS_vfi = cS;
                cS_vfi.theta_t = cS.theta_path(t);
                current_pps_active = (cS.sim_years(t) >= cS.pps_activation_year);
                cS_vfi.pps_active = current_pps_active;

                % 根据当期pps状态，设置正确的网格
                if ~cS_vfi.pps_active && cS_vfi.nkpps > 1
                    % 仅在需要时才修改和重新生成网格
                    cS_vfi.nkpps = 1; cS_vfi.npps = 1;
                    cS_vfi = main_olg_v14_utils.generateGrids(cS_vfi);
                end

                V_current_by_age = cell(cS.aD_new, 1);

                % --- 2. 年龄反向循环 ---
                for a = cS.aD_new:-1:1
                    vPrime_kkppse_next = [];
                    if a < cS.aD_new
                        % 从 t+1 期获取价值函数
                        vPrime_kkppse_next_full = V_future_by_age{a+1};

                        % !!!!! 核心修复 !!!!!
                        % 检查是否处于PPS激活的边界。如果是，则对价值函数进行降维。
                        future_pps_active = false;
                        if t < T_sim
                            future_pps_active = (cS.sim_years(t+1) >= cS.pps_activation_year);
                        end

                        if ~current_pps_active && future_pps_active
                            % 当前 t 期 pps 未激活, 但下一期 t+1 激活了。
                            % 家庭预期进入 t+1 期时的 k_pps 为 0。
                            % 我们只取 vPrime 在 k_pps=0 (即第1个网格点) 的那个切片。
                            vPrime_kkppse_next = vPrime_kkppse_next_full(:, 1, :);
                        else
                            % 在其他所有情况下，维度都是一致的，直接使用即可。
                            vPrime_kkppse_next = vPrime_kkppse_next_full;
                        end
                    end

                    % --- 3. 求解当前(t,a)的策略和价值 ---
                    [~, kPol, cPpsPol, val] = main_olg_v14_utils.HHSolutionByAge_VFI_Vectorized(a, vPrime_kkppse_next, M_t, 0, M_t.b_t, paramS, cS_vfi);

                    % --- 4. 存储策略函数的插值器 ---
                    kPolInterp_age = cell(cS_vfi.nw, 1);
                    cPpsPolInterp_age = cell(cS_vfi.nw, 1);
                    for ie = 1:cS_vfi.nw
                        kPol_slice = squeeze(kPol(:,:,ie));
                        cPpsPol_slice = squeeze(cPpsPol(:,:,ie));

                        if cS_vfi.nk > 1 && cS_vfi.nkpps > 1
                            kPolInterp_age{ie} = griddedInterpolant({cS_vfi.kGridV, cS_vfi.kppsGridV}, kPol_slice, 'linear', 'none');
                            cPpsPolInterp_age{ie} = griddedInterpolant({cS_vfi.kGridV, cS_vfi.kppsGridV}, cPpsPol_slice, 'linear', 'none');
                        elseif cS_vfi.nk > 1
                            kPolInterp_age{ie} = griddedInterpolant(cS_vfi.kGridV, kPol_slice, 'linear', 'none');
                            cPpsPolInterp_age{ie} = griddedInterpolant(cS_vfi.kGridV, cPpsPol_slice, 'linear', 'none');
                        else
                            % 处理 nk=1 的情况
                            kPolInterp_age{ie} = @(x,y) kPol_slice;
                            cPpsPolInterp_age{ie} = @(x,y) cPpsPol_slice;
                        end
                    end
                    PolicyFunctions{t, a} = {kPolInterp_age, cPpsPolInterp_age};

                    % 存储当前价值函数矩阵，供t-1期使用
                    V_current_by_age{a} = val;
                end

                % 整个t期的价值函数已算完，将其赋给 V_future，为 t-1 期做准备
                V_future_by_age = V_current_by_age;
            end
            fprintf('✅ 所有策略函数求解完毕。\n');
        end

        function [K_pvt_next, K_pps_next, C_t, Total_Cpps_t, Total_PpsTax_t, total_accidental_bequest, M_t_complete] = simulate_private_capital_forward(t, K_pvt_t, B_p_t, K_pps_t, B_g_t, Z_t, A_t, cS, paramS, eIdxM)
            % [版本：遗赠归政府 - 存活率字段修正]
            % 核心变化：将所有 cS.s_1yr_transitionV 的调用替换为 cS.s_pathV

            % --- 1. 计算宏观价格和政策 (不变) ---
            K_physical_t = K_pvt_t + K_pps_t + B_p_t - B_g_t;
            if K_physical_t <= 0, K_physical_t = 1e-8; end
            paramS_t = paramS; paramS_t.ageMassV = Z_t;
            [~, L_t] = main_olg_v14_utils.LaborSupply_Huggett(eIdxM, cS, paramS_t, Z_t);
            M_prices = main_olg_v14_utils.get_prices_at_t(K_physical_t, L_t, A_t, cS);
            r_mkt_t = M_prices.r_mkt_t;
            w_t = M_prices.w_t;
            Y_t = M_prices.Y_t;
            r_net_period = r_mkt_t * (1 - cS.tau_k);
            mass_workers_t = sum(Z_t(1:cS.aR_new));
            if mass_workers_t > 0, avg_wage_period = w_t * L_t / mass_workers_t; else, avg_wage_period = 0; end
            b_period = cS.rho_prime_payg * avg_wage_period;

            M_t = struct();
            M_t.current_t = t; M_t.K_physical_t = K_physical_t; M_t.L_t = L_t; M_t.Y_t = Y_t;
            M_t.w_t = w_t; M_t.b_t = b_period; M_t.r_mkt_t = r_mkt_t; M_t.r_net_period = r_net_period;
            M_t.K_pvt_t = K_pvt_t; M_t.K_pps_t = K_pps_t; M_t.B_p_t = B_p_t; M_t.B_g_t = B_g_t;

            cS_vfi = cS;
            if ~isfield(cS_vfi, 'pps_active'), cS_vfi.pps_active = (cS.sim_years(t) >= cS_vfi.pps_activation_year); end
            if ~cS_vfi.pps_active, cS_vfi.nkpps = 1; cS_vfi.npps = 1; cS_vfi = main_olg_v14_utils.generateGrids(cS_vfi); end
            if ~isfield(cS_vfi, 'theta_t'), cS_vfi.theta_t = cS.theta_path(t); end

            % --- 2. [核心简化] 单次求解和模拟 ---
            tr_per_survivor_eq = 0;
            [cPolM_final, kPolM_final, cPpsPolM_final, ~] = main_olg_v14_utils.HHSolution_VFI_Huggett(M_t, tr_per_survivor_eq, paramS, cS_vfi);
            [kHistM, kPpsHistM, cHistM, cppsHistM] = main_olg_v14_utils.HHSimulation_olgm(kPolM_final, cPpsPolM_final, cPolM_final, eIdxM, M_t, paramS, cS_vfi, tr_per_survivor_eq);

            % --- 3. 聚合 ---
            mean_k_prime_by_age = mean(kHistM(:, 2:end), 1);
            mean_k_pps_prime_by_age = mean(kPpsHistM(:, 2:end), 1);

            % [核心修正] 使用新的存活率字段 s_pathV
            K_pvt_next = sum(mean_k_prime_by_age' .* Z_t .* cS.s_pathV);
            K_pps_next = sum(mean_k_pps_prime_by_age' .* Z_t .* cS.s_pathV);

            C_t = sum(mean(cHistM, 1)' .* Z_t);
            Total_Cpps_t = sum(mean(cppsHistM, 1)' .* Z_t);

            % [核心修正] 使用新的存活率字段 s_pathV
            prob_death_by_age = (1 - cS.s_pathV);
            total_accidental_bequest = sum(mean_k_prime_by_age' .* Z_t .* prob_death_by_age) + ...
                sum(mean_k_pps_prime_by_age' .* Z_t .* prob_death_by_age);

            % --- 4. 精确核算PPS提现税 (不变) ---
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

            % --- 5. 打包返回 (不变) ---
            M_t_complete = M_t;
            M_t_complete.C_t = C_t;
            M_t_complete.Total_Cpps_t = Total_Cpps_t;
            M_t_complete.Total_PpsTax_t = Total_PpsTax_t;
        end

        function [kHistM_out, kPpsHistM_out, cHistM_out, cppsHistM_out] = HHSimulation_olgm(kPolM, cPpsPolM_choice, cPolM_consump, eIdxM_group, M_sim, paramS_sim, cS_sim, tr_per_capita_sim)
            % [v13 - 最终裁决版]
            % 这个版本在逻辑上与v12一致，但通过最清晰的实现方式，
            % 旨在终结任何关于变量覆盖或错误索引的怀疑。

            % --- 0. 初始化 ---
            nSim = size(eIdxM_group, 1);
            aD = cS_sim.aD_new;
            nw = cS_sim.nw;

            kHistM_out = zeros(nSim, aD + 1);
            kPpsHistM_out = zeros(nSim, aD + 1);
            cHistM_out = zeros(nSim, aD);
            cppsHistM_out = zeros(nSim, aD);

            % --- 1. 创建插值器库 ---
            % 为避免任何可能的混淆，我们使用一个不同的变量名
            PolicyInterp_k = cell(nw, aD);
            PolicyInterp_cpps = cell(nw, aD);
            for ia = 1:aD
                for ie = 1:nw
                    kPolM_slice = squeeze(kPolM(:, :, ie, ia));
                    cPpsPolM_slice = squeeze(cPpsPolM_choice(:, :, ie, ia));

                    if cS_sim.nk > 1 && cS_sim.nkpps > 1
                        PolicyInterp_k{ie,ia} = griddedInterpolant({cS_sim.kGridV, cS_sim.kppsGridV}, kPolM_slice, 'linear', 'none');
                        PolicyInterp_cpps{ie,ia} = griddedInterpolant({cS_sim.kGridV, cS_sim.kppsGridV}, cPpsPolM_slice, 'linear', 'none');
                    elseif cS_sim.nk > 1
                        PolicyInterp_k{ie,ia} = griddedInterpolant(cS_sim.kGridV, kPolM_slice, 'linear', 'none');
                        PolicyInterp_cpps{ie,ia} = griddedInterpolant(cS_sim.kGridV, cPpsPolM_slice, 'linear', 'none');
                    else
                        PolicyInterp_k{ie,ia} = @(x,y) kPolM_slice;
                        PolicyInterp_cpps{ie,ia} = @(x,y) cPpsPolM_slice;
                    end
                end
            end

            % --- 2. 按年龄逐期模拟 ---
            for a_idx = 1:aD
                k_now = kHistM_out(:, a_idx);
                k_pps_now = kPpsHistM_out(:, a_idx);

                k_next_decision = zeros(nSim, 1);
                cpps_decision = zeros(nSim, 1);

                for ie_idx = 1:nw
                    % 找到所有在当前年龄a_idx，效率状态为ie_idx的个体
                    current_hh_idx = (eIdxM_group(:, a_idx) == ie_idx);
                    if ~any(current_hh_idx), continue; end

                    % 从库中取出【唯一确定】的插值器
                    k_interpolator = PolicyInterp_k{ie_idx, a_idx};
                    cpps_interpolator = PolicyInterp_cpps{ie_idx, a_idx};

                    % 获取这些家庭的当前状态
                    k_now_e = k_now(current_hh_idx);
                    k_pps_now_e = k_pps_now(current_hh_idx);

                    % 应用插值器
                    if cS_sim.nk > 1 && cS_sim.nkpps > 1
                        k_next_decision(current_hh_idx) = k_interpolator(k_now_e, k_pps_now_e);
                        cpps_decision(current_hh_idx) = cpps_interpolator(k_now_e, k_pps_now_e);
                    else
                        k_next_decision(current_hh_idx) = k_interpolator(k_now_e);
                        cpps_decision(current_hh_idx) = cpps_interpolator(k_now_e);
                    end
                end

                k_next_decision(isnan(k_next_decision)) = cS_sim.kMin;
                cpps_decision(isnan(cpps_decision)) = 0;

                kHistM_out(:, a_idx + 1) = max(cS_sim.kMin, k_next_decision);
                cppsHistM_out(:, a_idx) = max(0, cpps_decision);

                % --- 3. 反推消费 (逻辑与之前完全一致) ---
                tr_this_age = 0; if a_idx > 1, tr_this_age = tr_per_capita_sim; end
                if a_idx <= cS_sim.aR_new
                    labor_income_gross = M_sim.w_t .* cS_sim.ageEffV_new(a_idx) .* paramS_sim.leGridV(eIdxM_group(:, a_idx));
                    k_return = k_now .* (1 + M_sim.r_net_period);
                    total_inflow = k_return + labor_income_gross + tr_this_age;
                    payg_tax = cS_sim.theta_t .* labor_income_gross;
                    labor_tax = cS_sim.tau_l .* max(0, labor_income_gross - cppsHistM_out(:, a_idx) - payg_tax);
                    k_prime_outflow = kHistM_out(:, a_idx + 1);
                    cpps_outflow = cppsHistM_out(:, a_idx);
                    total_outflow_non_c = payg_tax + labor_tax + k_prime_outflow + cpps_outflow;
                    c_expend = total_inflow - total_outflow_non_c;
                else
                    k_return = k_now .* (1 + M_sim.r_net_period);
                    b_age_val = M_sim.b_t;
                    pps_withdrawal_gross = 0; if cS_sim.pps_active, pps_withdrawal_gross = k_pps_now .* cS_sim.pps_withdrawal_rate; end
                    pps_withdrawal_net = pps_withdrawal_gross .* (1 - cS_sim.pps_tax_rate_withdrawal);
                    total_inflow = k_return + b_age_val + tr_this_age + pps_withdrawal_net;
                    k_prime_outflow = kHistM_out(:, a_idx + 1);
                    c_expend = total_inflow - k_prime_outflow;
                end
                cHistM_out(:, a_idx) = max(cS_sim.cFloor, c_expend ./ (1 + cS_sim.tau_c));

                % --- 4. 更新PPS资产 (逻辑与之前完全一致) ---
                pps_return_factor = 1 + M_sim.r_mkt_t + cS_sim.pps_return_rate_premium;
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

        % =========================================================================
        % == 在 main_olg_v14_utils.m 中，请【替换】此函数
        % =========================================================================

        function Results = simulate_forward_with_policies(initial_dist_k, initial_dist_kpps, PolicyFunctions, T_sim, Z_path, A_path, w_path, r_mkt_path, b_path, L_path, cS, paramS, eIdxM)
            % =========================================================================
            % == 【函数已最终修正】simulate_forward_with_policies
            % == 核心修正：
            % == 1. (来自上次修正) 严格使用传入的价格路径(w, r, b, L)进行模拟和核算。
            % == 2. (本次核心修正) 在每个时期结束后，正确处理代际更替（Aging）。
            % ==    将 t 期 a 年龄组选择的储蓄 k'，正确地传递给 t+1 期 a+1 年龄组作为期初资产。
            % =========================================================================
            fprintf('--- [转型路径] 开始基于策略函数进行前向模拟 (t=1 to %d) ---\n', T_sim);

            % 初始化资产历史面板
            k_panel = initial_dist_k(:, 1:cS.aD_new);
            kpps_panel = initial_dist_kpps(:, 1:cS.aD_new);

            % 初始化结果存储结构
            Results = struct();
            Results.K_pvt_path = zeros(T_sim, 1);
            Results.K_pps_path = zeros(T_sim, 1);
            Results.B_p_path = zeros(T_sim, 1);
            Results.B_g_path = zeros(T_sim, 1);
            Results.K_physical_path = zeros(T_sim, 1);
            Results.C_path = zeros(T_sim, 1);
            Results.Y_path = zeros(T_sim, 1);
            Results.I_path = zeros(T_sim, 1);
            Results.G_path = zeros(T_sim, 1);
            Results.L_path = L_path;
            Results.w_path = w_path;
            Results.r_mkt_path = r_mkt_path;

            % --- 时间前向循环 ---
            for t = 1:T_sim
                if mod(t,10) == 0, fprintf('   正在模拟时期 t = %d\n', t); end

                Z_t = Z_path(:, t);
                A_t = A_path(t);

                % --- 1. 聚合当期资本，使用传入的价格 ---
                K_pvt_t = sum(mean(k_panel, 1)' .* Z_t);
                K_pps_t = sum(mean(kpps_panel, 1)' .* Z_t);
                B_p_t = Results.B_p_path(t);
                B_g_t = Results.B_g_path(t);
                K_physical_t = K_pvt_t + K_pps_t + B_p_t - B_g_t;
                if K_physical_t <= 0, K_physical_t = 1e-8; end

                w_t = w_path(t);
                r_mkt_t = r_mkt_path(t);
                b_t = b_path(t);
                L_t = L_path(t);
                Y_t = A_t * (K_physical_t^cS.alpha) * (L_t^(1-cS.alpha));

                % --- 2. 应用策略，模拟当期决策 ---
                M_t = struct('w_t', w_t, 'r_mkt_t', r_mkt_t, 'r_net_period', r_mkt_t*(1-cS.tau_k), 'b_t', b_t);
                M_t.current_t = t;
                cS_sim = cS;
                cS_sim.theta_t = cS.theta_path(t);
                if isfield(cS, 'sim_years'), cS_sim.pps_active = (cS.sim_years(t) >= cS_sim.pps_activation_year); else, cS_sim.pps_active = true; end
                Policies_t = PolicyFunctions(t, :);

                % HHSimulation 返回的是一个历史矩阵，记录了从期初 k 到期末决策 k' 的过程
                [kHist_t, kppsHist_t, c_panel, ~] = main_olg_v14_utils.HHSimulation_olgm_transition(k_panel, kpps_panel, Policies_t, eIdxM, M_t, paramS, cS_sim);

                % --- 3. 聚合当期结果并存储 ---
                Results.C_path(t) = sum(mean(c_panel, 1)' .* Z_t);
                Results.Y_path(t) = Y_t;
                Results.G_path(t) = cS.gov_exp_frac_Y * Y_t;
                Results.K_pvt_path(t) = K_pvt_t;
                Results.K_pps_path(t) = K_pps_t;
                Results.K_physical_path(t) = K_physical_t;

                % --- 4. 【核心修正】准备下一期的期初资产，正确处理代际更替 ---

                % kHist_t 的维度是 nSim x (aD+1)。第 a+1 列存储的是 a 年龄组选择的 k'。
                % 我们先将所有人的 k' 决策提取出来，形成一个 nSim x aD 的决策矩阵。
                k_prime_choices = kHist_t(:, 2:end);
                kpps_prime_choices = kppsHist_t(:, 2:end);

                % 创建 t+1 期的空白资产面板
                k_panel_tp1 = zeros(size(k_panel));
                kpps_panel_tp1 = zeros(size(kpps_panel));

                % 代际更替：t期 a-1 年龄组的储蓄决策，成为 t+1期 a 年龄组的期初资产。
                % k_prime_choices 的第 a-1 列是 a-1 年龄组的决策。
                k_panel_tp1(:, 2:cS.aD_new) = k_prime_choices(:, 1:cS.aD_new-1);
                kpps_panel_tp1(:, 2:cS.aD_new) = kpps_prime_choices(:, 1:cS.aD_new-1);

                % 新生代（年龄组1）的期初资产为0
                % k_panel_tp1(:, 1) 已经是0，无需操作。

                % 更新面板状态，用于下一次循环
                k_panel = k_panel_tp1;
                kpps_panel = kpps_panel_tp1;
            end

            % --- 5. 计算投资路径 ---
            % 投资 = Gross Investment = K_physical(t+1) - (1-depreciation)*K_physical(t)
            for t = 1:T_sim-1
                K_physical_tp1 = Results.K_physical_path(t+1);
                Results.I_path(t) = K_physical_tp1 - (1-cS.ddk)*Results.K_physical_path(t);
            end
            Results.I_path(T_sim) = Results.I_path(T_sim-1); % 简单外推最后一期

            fprintf('✅ 前向模拟完成。\n');
        end


        function [kHistM_out, kPpsHistM_out, cHistM_out, cppsHistM_out] = HHSimulation_olgm_transition(k_panel_in, kpps_panel_in, Policies_t, eIdxM_group, M_sim, paramS_sim, cS_sim)
            % =========================================================================
            % == 【函数已最终修正】HHSimulation_olgm_transition
            % == 核心修正：
            % == 1. 此函数现在【严格使用】传入的 k_panel_in 和 kpps_panel_in 作为
            % ==    每个家庭在每个年龄组的期初资产。
            % == 2. 不再初始化 kHistM_out 和 kPpsHistM_out 为零，而是用传入的
            % ==    面板数据作为历史记录的起点。
            % == 3. 循环中直接使用 k_panel_in(:, a_idx) 作为当前状态，而不是
            % ==    依赖于前一列的计算结果。
            % =========================================================================

            nSim = size(eIdxM_group, 1);
            aD = cS_sim.aD_new;

            % --- 1. 【核心修正】初始化输出矩阵，并将传入的面板作为历史记录的“期初”部分 ---
            kHistM_out = zeros(nSim, aD + 1);
            kPpsHistM_out = zeros(nSim, aD + 1);
            cHistM_out = zeros(nSim, aD);
            cppsHistM_out = zeros(nSim, aD);

            % kHistM_out 的前 aD 列代表 nSim 个家庭在 aD 个年龄组的期初资产
            kHistM_out(:, 1:aD) = k_panel_in;
            kPpsHistM_out(:, 1:aD) = kpps_panel_in;

            tr_per_capita_sim = 0; % 在转型路径中，意外遗赠通过其他机制处理

            % --- 2. 按年龄逐期模拟决策 ---
            for a_idx = 1:aD
                % --- 【核心修正】直接从传入的面板中获取当前状态 ---
                k_now = k_panel_in(:, a_idx);
                k_pps_now = kpps_panel_in(:, a_idx);

                % 从策略函数库中获取当前年龄的插值器
                Policy_age = Policies_t{a_idx};
                kPolInterp_cell = Policy_age{1};
                cPpsPolInterp_cell = Policy_age{2};

                k_next_decision = zeros(nSim, 1);
                cpps_decision = zeros(nSim, 1);

                % 遍历所有效率冲击状态，为每个家庭找到其决策
                for ie = 1:cS_sim.nw
                    idx_sim = find(eIdxM_group(:, a_idx) == ie);
                    if isempty(idx_sim), continue; end

                    kPolInterp = kPolInterp_cell{ie};
                    cPpsPolInterp = cPpsPolInterp_cell{ie};

                    k_now_e = k_now(idx_sim);
                    k_pps_now_e = k_pps_now(idx_sim);

                    % 使用插值器得到决策
                    if numel(kPolInterp.GridVectors) > 1
                        k_next_decision(idx_sim) = kPolInterp(k_now_e, k_pps_now_e);
                        cpps_decision(idx_sim) = cPpsPolInterp(k_now_e, k_pps_now_e);
                    else
                        k_next_decision(idx_sim) = kPolInterp(k_now_e);
                        cpps_decision(idx_sim) = cPpsPolInterp(k_now_e);
                    end
                end

                % 处理插值可能产生的NaN
                k_next_decision(isnan(k_next_decision)) = cS_sim.kMin;
                cpps_decision(isnan(cpps_decision)) = 0;

                % --- 3. 存储当期决策 k' 和 cpps ---
                cppsHistM_out(:, a_idx) = max(0, cpps_decision);

                % k' 决策被存储在历史矩阵的下一列
                kHistM_out(:, a_idx + 1) = max(cS_sim.kMin, k_next_decision);

                % --- 4. 反推消费 (此逻辑保持不变，它依赖于正确的 k_now 和 k_prime) ---
                tr_this_age = 0; if a_idx > 1, tr_this_age = tr_per_capita_sim; end
                if a_idx <= cS_sim.aR_new
                    labor_income_gross = M_sim.w_t .* cS_sim.ageEffV_new(a_idx) .* paramS_sim.leGridV(eIdxM_group(:, a_idx));
                    k_return = k_now .* (1 + M_sim.r_net_period);
                    total_inflow = k_return + labor_income_gross + tr_this_age;
                    payg_tax = cS_sim.theta_t .* labor_income_gross;
                    labor_tax = cS_sim.tau_l .* max(0, labor_income_gross - cppsHistM_out(:, a_idx) - payg_tax);
                    k_prime_outflow = kHistM_out(:, a_idx + 1);
                    cpps_outflow = cppsHistM_out(:, a_idx);
                    total_outflow_non_c = payg_tax + labor_tax + k_prime_outflow + cpps_outflow;
                    c_expend = total_inflow - total_outflow_non_c;
                else
                    k_return = k_now .* (1 + M_sim.r_net_period);
                    b_age_val = M_sim.b_t;
                    pps_withdrawal_gross = 0; if cS_sim.pps_active, pps_withdrawal_gross = k_pps_now .* cS_sim.pps_withdrawal_rate; end
                    pps_withdrawal_net = pps_withdrawal_gross .* (1 - cS_sim.pps_tax_rate_withdrawal);
                    total_inflow = k_return + b_age_val + tr_this_age + pps_withdrawal_net;
                    k_prime_outflow = kHistM_out(:, a_idx + 1);
                    c_expend = total_inflow - k_prime_outflow;
                end
                cHistM_out(:, a_idx) = max(cS_sim.cFloor, c_expend ./ (1 + cS_sim.tau_c));

                % --- 5. 更新PPS资产 (此逻辑保持不变) ---
                pps_return_factor = 1 + M_sim.r_mkt_t + cS_sim.pps_return_rate_premium;
                pps_withdrawal_gross_for_update = 0;
                if a_idx > cS_sim.aR_new && cS_sim.pps_active, pps_withdrawal_gross_for_update = k_pps_now .* cS_sim.pps_withdrawal_rate; end
                k_pps_after_withdrawal = k_pps_now - pps_withdrawal_gross_for_update;
                k_pps_next_unclamped = (k_pps_after_withdrawal + cppsHistM_out(:, a_idx)) * pps_return_factor;
                kPpsHistM_out(:, a_idx + 1) = max(cS_sim.kppsMin, min(cS_sim.kppsMax, k_pps_next_unclamped));
            end
        end







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