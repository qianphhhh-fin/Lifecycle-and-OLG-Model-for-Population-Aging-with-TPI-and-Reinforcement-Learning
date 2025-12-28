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

        % =========================================================================
        % == 替换 utils 文件中的 HHPrices_Huggett 函数 ==
        % =========================================================================
        % =========================================================================
        % == 替换 utils 文件中的 HHPrices_Huggett 函数 (最终正确版) ==
        % =========================================================================
        function [r_mkt_period, w_t, Y_period] = HHPrices_Huggett(K, L, A_t, cS)
            % [最终审计通过版]
            % - 返回所有“每期”的量。
            % - 简化折旧逻辑，确保欧拉定理通过构造恒等。
            
            if K <= 0, K = 1e-8; end
            if L <= 0, L = 1e-8; end

            % 1. 计算每期总产出
            Y_period = A_t * (K.^cS.alpha) .* (L.^(1-cS.alpha));

            % 2. 计算资本的边际产出毛额 (MPK) 和工资
            MPK_period = cS.alpha * Y_period / K;
            w_t = (1-cS.alpha) * Y_period / L;

            % 3. 计算资本的净回报率 (r)
            %    r = MPK - delta, 其中 delta 是每期的折旧率 cS.ddk
            r_mkt_period = MPK_period - cS.ddk;
        end
        % =========================================================================
        % == 替换 utils 文件中的 get_prices_and_policy_at_t 函数 ==
        % =========================================================================
        % =========================================================================
        % == 在 utils 文件中，替换 get_prices_and_policy_at_t (最终版) ==
        % =========================================================================
% =========================================================================
% == 在 utils 文件中，替换 get_prices_and_policy_at_t (最终版 v2) ==
% =========================================================================
function M_t = get_prices_and_policy_at_t(t, K_pvt_t, B_p_t, K_pps_t, B_g_t, Z_t_norm, A_t, cS, paramS, eIdxM)
    % [最终审计修正版 v2]
    % - 修正了 Y_t 和 r_mkt_5y 的计算，确保供给-收入恒等式成立。
    
    K_physical_t = K_pvt_t + K_pps_t + B_p_t - B_g_t;
    if K_physical_t <= 0, K_physical_t = 1e-8; end
    
    paramS_t = paramS; paramS_t.ageMassV = Z_t_norm;
    [~, L_t] = main_olg_v14_utils.LaborSupply_Huggett(eIdxM, cS, paramS_t, Z_t_norm);

    [r_mkt_annual, w_t, Y_annual] = main_olg_v14_utils.HHPrices_Huggett(K_physical_t, L_t, A_t, cS);

    % [核心修正] 统一使用“单利”逻辑计算5年期总量
    % 1. 5年期总产出
    Y_t_5y = Y_annual * cS.time_Step;
    
    % 2. 5年期总劳动收入 = w_t * L_t * 5
    
    % 3. 5年期总资本净收入 (r_mkt * K * 5)
    r_mkt_5y_simple = r_mkt_annual * cS.time_Step;
    
    % 4. 5年期总折旧 = ddk_annual * K * 5 (等价于 ddk_5y * K)
    
    % 5. 家庭面对的税后回报因子 (这个仍然需要复利，因为它用于跨期决策)
    r_net_annual = r_mkt_annual * (1 - cS.tau_k);
    R_net_factor_5y = (1 + r_net_annual)^cS.time_Step;
    
    mass_workers_t = sum(Z_t_norm(1:cS.aR_new));
    if mass_workers_t > 0, avg_wage_t = w_t * L_t / mass_workers_t; else, avg_wage_t = 0; end
    b_t = cS.rho_prime_payg * avg_wage_t;

    M_t = struct();
    M_t.current_t = t;
    M_t.K_physical_t = K_physical_t;
    M_t.L_t = L_t;
    M_t.Y_t = Y_t_5y; % 存储5年总产出
    M_t.w_t = w_t;
    M_t.b_t = b_t;
    M_t.r_mkt_annual = r_mkt_annual;
    M_t.r_net_annual = r_net_annual;
    M_t.R_net_factor_5y = R_net_factor_5y; % 用于VFI的复利因子
    M_t.r_mkt_5y_simple = r_mkt_5y_simple; % [新增] 用于宏观会计的单利利率
    M_t.K_pvt_t = K_pvt_t; M_t.K_pps_t = K_pps_t; M_t.B_p_t = B_p_t; M_t.B_g_t = B_g_t;
end        % =========================================================================
        % == 替换 utils 文件中的 update_pension_fund 函数 ==
        % =========================================================================
        function [B_p_next, NetSurplus_annual] = update_pension_fund(B_p_t, M_t, Z_t_norm, cS)
            % [最终修正版] 返回当期盈余/赤字
            theta_t = cS.theta_path(M_t.current_t);
            PensionRevenue_annual = theta_t * M_t.w_t * M_t.L_t;
            mass_retirees_t = sum(Z_t_norm(cS.aR_new+1:end));
            PensionOutlay_annual = M_t.b_t * mass_retirees_t;
            NetSurplus_annual = PensionRevenue_annual - PensionOutlay_annual;

            r_a = M_t.r_mkt_annual;
            R_5y = (1 + r_a)^cS.time_Step;
            if abs(r_a) > 1e-6, FlowAcc_5y = (R_5y - 1) / r_a; else, FlowAcc_5y = cS.time_Step; end

            B_p_next = B_p_t * R_5y + NetSurplus_annual * FlowAcc_5y;
            B_p_next = max(0, B_p_next); % 基金本身不能为负
        end
        % --- START OF NEW FUNCTION in main_olg_v14_utils.m ---

        % --- START OF CORRECTED FUNCTION in main_olg_v14_utils.m ---

        % --- START OF CORRECTED FUNCTION in main_olg_v14_utils.m ---

        % =========================================================================
        % == 在 main_olg_v14_utils.m 中，替换整个 calcaulte_theta_payg_path 函数 ==
        % =========================================================================

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
        % --- END OF CORRECTED FUNCTION ---



        % =========================================================================
        % == 替换 utils 文件中的 update_gov_debt 函数 ==
        % =========================================================================
        function [B_g_next, G_t] = update_gov_debt(B_g_t, C_t, M_t, Total_Cpps_t, K_pps_t, Total_PpsTax_t, cS, PensionSurplus_annual)
            % [最终修正版] 将养老金盈余/赤字计入政府预算

            % 计算年化流量
            G_t = cS.gov_exp_frac_Y * M_t.Y_t; % G是5年总量
            G_annual = G_t / cS.time_Step;
            C_annual = C_t / cS.time_Step;
            Total_Cpps_annual = Total_Cpps_t / cS.time_Step;
            % Total_PpsTax_t 已经是年化流量，无需再除
            Total_PpsWithdrawalTax_annual = Total_PpsTax_t;

            GrossLaborIncome_annual = M_t.w_t * M_t.L_t;
            TaxableLaborIncome_annual = GrossLaborIncome_annual - Total_Cpps_annual;
            LaborTaxRevenue_annual = cS.tau_l * max(0, TaxableLaborIncome_annual);

            r_a = M_t.r_mkt_annual;
            TaxableCapitalStock_t = M_t.K_physical_t - K_pps_t;
            CapitalTaxRevenue_annual = cS.tau_k * r_a * TaxableCapitalStock_t;

            ConsumptionTaxRevenue_annual = cS.tau_c * C_annual;
            TaxRevenue_annual = LaborTaxRevenue_annual + CapitalTaxRevenue_annual + ConsumptionTaxRevenue_annual + Total_PpsWithdrawalTax_annual;

            % [核心] 政府的广义赤字 = 自身赤字 - 养老金盈余
            PrimaryDeficit_annual = G_annual - TaxRevenue_annual - PensionSurplus_annual;

            % 计算5年期回报因子和流量累积因子
            R_5y = (1 + r_a)^cS.time_Step;
            if abs(r_a) > 1e-6, FlowAcc_5y = (R_5y - 1) / r_a; else, FlowAcc_5y = cS.time_Step; end

            % 应用正确的多年期演化方程
            B_g_next = B_g_t * R_5y + PrimaryDeficit_annual * FlowAcc_5y;
        end

        function cS = ParameterValues_HuggettStyle()
            % [VFI最终对齐版 v2] - 添加了均衡求解器参数

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

        function M_prices = get_prices_at_t(K_total_t, Z_t_norm, cS, paramS_t, eIdxM)
            % [v14.3辅助] 仅根据总资本和人口，计算当期的要素价格和宏观总量。

            paramS_t.ageMassV = Z_t_norm;
            [~, L_t] = main_olg_v14_utils.LaborSupply_Huggett(eIdxM, cS, paramS_t, Z_t_norm);

            [R_mkt_factor, w_t] = main_olg_v14_utils.HHPrices_Huggett(K_total_t, L_t, cS);
            r_mkt_t = R_mkt_factor - 1;
            % [修正] 移除双重折旧扣除：家庭获得的净回报率应该是税后的市场利率
            r_net_t = r_mkt_t * (1 - cS.tau_k);

            M_prices = struct();
            M_prices.K_total_t = K_total_t;
            M_prices.L_t = L_t;
            M_prices.Y_t = cS.A * (K_total_t^cS.alpha) * (L_t^(1-cS.alpha));
            M_prices.w_t = w_t;
            M_prices.r_mkt_t = r_mkt_t;
            M_prices.r_net_t = r_net_t;
        end
        % --- 在 main_olg_v14_utils.m 中，替换 check_gbc_residual 函数 ---

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

        % --- 在 main_olg_v14_utils.m 的 methods (Static) 块中新增此函数 ---

        % --- 在 main_olg_v14_utils.m 的 methods (Static) 块中新增此函数 ---
        % --- 您可以删除或注释掉旧的 solve_steady_state 函数 ---

        % --- 在 main_olg_v14_utils.m 的 methods (Static) 块中新增此函数 ---
        % --- 您可以删除或注释掉旧的 solve_steady_state_endo_rho 函数 ---

        % --- 在 main_olg_v14_utils.m 的 methods (Static) 块中新增此函数 ---
        % --- 您可以删除或注释掉所有旧的 solve_... 函数 ---

        function [ss, eq_found] = solve_steady_state_with_fund(Z_ss_norm, B_p_Y_ratio_target, cS, paramS, eIdxM)
            % [v14.3核心, 已修正K_pps聚合] 求解包含目标养老金基金的稳态均衡。
            % - 迭代变量: 总资本 K_total
            % - 目标: 找到一个 K，使得 K_guess = K_pvt_model(K) + B_p_target(K)

            fprintf('\n--- 开始求解初始稳态 (有目标基金 B_p, 无初始PPS) ---\n');

            K_guess = 1.8816; % 初始猜测总资本
            max_iter = 100;
            tol = 1e-4;
            damp = 0.5;
            eq_found = false;

            fprintf('%5s | %12s | %12s | %12s | %12s | %12s | %12s\n', ...
                'Iter', 'K_guess', 'K_model', 'K_pvt', 'B_p', 'r_net', 'K_error');
            fprintf('%s\n', repmat('-', 90, 1));

            for iter = 1:max_iter
                % --- 步骤 1: 基于 K_guess 计算当期价格和宏观量 ---
                paramS_ss = paramS;
                paramS_ss.ageMassV = Z_ss_norm;

                % 在这个函数调用中，B_p_t 只是一个占位符，因为 K_total_t 直接由 K_guess 决定
                M_ss_prices = main_olg_v14_utils.get_prices_at_t(K_guess, Z_ss_norm, cS, paramS_ss, eIdxM);



                % --- 步骤 3: 求解家庭问题 (VFI) 并聚合私人资本 ---
                % 家庭面临的环境 M_ss 包含了刚刚计算出的价格
                M_ss = M_ss_prices;
                M_ss.current_t = -1; % 标记为稳态求解

                % --- [核心修正] 步骤 1.5: 计算当期政策变量 b_t 并添加到 M_ss ---
                mass_workers_ss = sum(Z_ss_norm(1:cS.aR_new));
                avg_wage_ss = M_ss.w_t * M_ss.L_t / mass_workers_ss;
                M_ss.b_t = cS.rho_prime_payg * avg_wage_ss; % 使用固定的替代率


                % --- 步骤 2: 计算目标基金规模 ---
                B_p_target = B_p_Y_ratio_target * M_ss.Y_t;

                % --- 修正结束 ---

                % 标记这是初始稳态求解，以禁用PPS
                paramS_ss.is_initial_steady_state = true;

                [K_pvt_model, K_pps_model, C_ss_model, Total_PpsTax_t1] = main_olg_v14_utils.simulate_private_capital_forward(M_ss, Z_ss_norm, cS, paramS_ss, eIdxM);

                % --- 步骤 4: 计算模型的总资本需求 ---
                K_total_model = K_pvt_model + B_p_target + K_pps_model;

                % --- 步骤 5: 计算误差并检查收敛 ---
                K_error = K_guess - K_total_model;

                fprintf('%5d | %12.4f | %12.4f | %12.4f | %12.4f | %11.2f%% | %12.3e\n', ...
                    iter, K_guess, K_total_model, K_pvt_model, B_p_target, M_ss_prices.r_net_t*100, K_error);

                if abs(K_error) < tol
                    fprintf('✅ 初始稳态均衡收敛！\n');
                    eq_found = true;
                    break;
                end

                % --- 步骤 6: 更新猜测 ---
                K_guess = (1 - damp) * K_guess + damp * K_total_model;
            end

            if ~eq_found, warning('初始稳态均衡未在最大迭代次数内收敛。'); end

            % --- 打包最终稳态结果 ---
            ss = struct();
            ss.K_total = K_total_model;
            ss.K_pvt = K_pvt_model;
            ss.B_p = B_p_target;
            ss.K_pps = K_pps_model;
            ss.L = M_ss_prices.L_t;
            ss.Y = M_ss_prices.Y_t;
            ss.w = M_ss_prices.w_t;
            ss.r_net = M_ss_prices.r_net_t;
            ss.C = C_ss_model;
            ss.G = cS.gov_exp_frac_Y * ss.Y;
            ss.B_g = cS.gov_debt_frac_Y * ss.Y;
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

        % --- 在 main_olg_v14_utils.m 中，替换整个 HHIncome_Huggett 函数 ---

        % =========================================================================
        % == 在 utils 文件中，替换 HHIncome_Huggett 函数 (最终修正版) ==
        % =========================================================================
        % =========================================================================
        % == 在 utils 文件中，替换 HHIncome_Huggett (最终版) ==
        % =========================================================================
        function [period_net_income, pps_contrib_expenditure] = HHIncome_Huggett(...
                R_k_mkt_factor, w_gross, TR_total, b_payg_val, ...
                c_pps_chosen, a_idx, paramS_hh, cS, epsilon_val, k_now_val)
            % [最终审计修正版]
            % - 接收【税前】市场回报因子 R_k_mkt_factor。
            % - 内部明确计算资本总收入、资本税、资本净收入。

            pps_contrib_expenditure = c_pps_chosen;

            % 1. 劳动净收入 (5年总量)
            labor_income_net_5y = 0;
            if a_idx <= cS.aR_new
                labor_income_gross_annual = w_gross * cS.ageEffV_new(a_idx) * epsilon_val;
                labor_income_gross_5y = labor_income_gross_annual * cS.time_Step;
                payg_tax_5y = cS.theta_t * labor_income_gross_5y;
                labor_tax_5y = cS.tau_l * max(0, labor_income_gross_5y - pps_contrib_expenditure);
                labor_income_net_5y = labor_income_gross_5y - payg_tax_5y - labor_tax_5y;
            end

            % 2. 资本净收入 (5年总量)
            capital_income_gross_5y = k_now_val .* (R_k_mkt_factor - 1);
            capital_tax_5y = cS.tau_k * capital_income_gross_5y;
            capital_income_net_5y = capital_income_gross_5y - capital_tax_5y;

            % 3. 其他净收入 (5年总量)
            pension_benefits_5y = 0;
            if a_idx > cS.aR_new
                pension_benefits_5y = b_payg_val * cS.time_Step;
            end
            transfer_income_5y = TR_total;

            % 4. 总的净收入流
            period_net_income = labor_income_net_5y + capital_income_net_5y + pension_benefits_5y + transfer_income_5y;
        end

        % =========================================================================
        % == 在 utils 文件中，替换 HHSolution_VFI_Huggett (最终版) ==
        % =========================================================================
        % =========================================================================
        % == 在 utils 文件中，修改 HHSolution_VFI_Huggett ==
        % =========================================================================
        function [cPolM_q, kPolM, cPpsPolM_choice, valM] = HHSolution_VFI_Huggett(...
                R_k_mkt_factor_vfi, w_gross_vfi, TR_total_vfi, bV_payg_vfi, r_mkt_annual_vfi, ... % <--- 新增参数
                paramS_vfi, cS, solverMethod)

            if nargin < 8, solverMethod = 'vectorized_grid'; end % nargin 增加 1

            valM = -Inf(cS.nk, cS.nkpps, cS.nw, cS.aD_new);
            cPolM_q = zeros(cS.nk, cS.nkpps, cS.nw, cS.aD_new);
            kPolM = zeros(cS.nk, cS.nkpps, cS.nw, cS.aD_new);
            cPpsPolM_choice = zeros(cS.nk, cS.nkpps, cS.nw, cS.aD_new);

            for a_idx = cS.aD_new : -1 : 1
                vPrime_kkppse_next = [];
                if a_idx < cS.aD_new
                    vPrime_kkppse_next = valM(:,:,:,a_idx+1);
                end

                % [核心修正] 将 r_mkt_annual_vfi 作为新参数传递
                [cPolM_q(:,:,:,a_idx), kPolM(:,:,:,a_idx), cPpsPolM_choice(:,:,:,a_idx), valM(:,:,:,a_idx)] = ...
                    main_olg_v14_utils.HHSolutionByAge_VFI_Huggett_VectorizedGrid(a_idx, vPrime_kkppse_next, ...
                    R_k_mkt_factor_vfi, w_gross_vfi, TR_total_vfi, bV_payg_vfi(a_idx), r_mkt_annual_vfi, ... % <--- 新增参数
                    paramS_vfi, cS);
            end
        end

        function [cPol_age_q, kPol_age, cPpsPol_age_choice, val_age] = HHSolutionByAge_VFI_Huggett_v9_GridSearch(...
                a_idx, vPrime_kkppse_next, R_k_net_factor_age, w_gross_age, TR_total_age, b_age_val, ...
                paramS_age, cS)
            % [最终修正版] 并行化的离散网格搜索VFI求解器

            % 初始化输出矩阵
            val_age    = -Inf(cS.nk, cS.nkpps, cS.nw);
            cPol_age_q = zeros(cS.nk, cS.nkpps, cS.nw);
            kPol_age   = zeros(cS.nk, cS.nkpps, cS.nw);
            cPpsPol_age_choice = zeros(cS.nk, cS.nkpps, cS.nw);

            % --- [核心修正] 最后一期逻辑 (带暖光遗赠动机) ---
            % --- [经过特殊情况处理的、真正终局的] 最后一期逻辑 ---
            if a_idx == cS.aD_new
                % 1. 计算当期的非财富收入 (只有养老金)
                [pension_income, ~, ~] = main_olg_v14_utils.HHIncome_Huggett(0, 0, 0, TR_total_age, b_age_val, 0, a_idx, paramS_age, cS, 0);

                % 2. 计算两种资产在赚取了各自的回报后的期末价值
                [K_grid, Kpps_grid] = ndgrid(cS.kGridV, cS.kppsGridV);
                final_value_k = K_grid .* R_k_net_factor_age;
                pps_return_factor = R_k_net_factor_age + cS.pps_return_rate_premium;
                final_value_kpps_gross = Kpps_grid .* pps_return_factor;
                final_value_kpps_net = final_value_kpps_gross * (1 - cS.pps_tax_rate_withdrawal);

                % 3. 计算可用于分配的总资源
                total_resources = final_value_k + final_value_kpps_net + pension_income;

                % ======================= [最终BUG修复的关键] =======================
                % 4. 根据是否有遗赠动机，决定如何分配最终资源

                if cS.phi_bequest > 1e-9 % 使用一个小的容差来判断是否有遗赠动机
                    % --- 情况A: 有遗赠动机 ---
                    % 求解消费和遗赠的最优分配
                    phi_adj = cS.phi_bequest * (1 + cS.tau_c)^(-1);
                    omega = (1 + phi_adj^(-1/cS.sigma))^(-1);

                    final_c_expenditure = omega .* total_resources;
                    final_bequest = (1 - omega) .* total_resources;

                    final_c = max(cS.cFloor, final_c_expenditure / (1 + cS.tau_c));

                    % 计算最终的价值函数
                    [~, util_c] = main_olg_v14_utils.CES_utility(final_c, cS.sigma, cS);
                    [~, util_b] = main_olg_v14_utils.CES_utility(final_bequest, cS.sigma_bequest, cS);
                    final_v = util_c + cS.phi_bequest * util_b;

                else
                    % --- 情况B: 没有遗赠动机 (phi_bequest = 0) ---
                    % 将所有资源用于消费，没有遗赠
                    final_bequest = zeros(size(total_resources));
                    final_c_expenditure = total_resources;
                    final_c = max(cS.cFloor, final_c_expenditure / (1 + cS.tau_c));

                    % 最终价值只来源于消费
                    [~, final_v] = main_olg_v14_utils.CES_utility(final_c, cS.sigma, cS);
                end
                % =================================================================

                % 5. 将结果赋给所有epsilon状态
                for ie = 1:cS.nw
                    cPol_age_q(:,:,ie) = final_c;
                    val_age(:,:,ie) = final_v;
                    kPol_age(:,:,ie) = final_bequest;
                    cPpsPol_age_choice(:,:,ie) = 0;
                end
                return;
            end

            % --- 非最后一期逻辑 ---

            % a. 计算期望未来价值矩阵 E[V']
            % (原始代码正确)
            EV_matrix = zeros(cS.nk, cS.nkpps, cS.nw);
            for ie_current = 1:cS.nw
                transition_probs = paramS_age.leTrProbM(ie_current, :);
                % 使用pagetimes进行高效的张量乘法
                EV_slice = pagemtimes(reshape(vPrime_kkppse_next, [cS.nk * cS.nkpps, cS.nw]), reshape(transition_probs, [cS.nw, 1]));
                EV_matrix(:, :, ie_current) = reshape(EV_slice, [cS.nk, cS.nkpps]);
            end

            % b. [修正] 创建插值器，并指定外插行为
            EV_interpolants = cell(cS.nw, 1);
            for ie_current = 1:cS.nw
                % 指定线性外插，避免NaN
                EV_interpolants{ie_current} = griddedInterpolant({cS.kGridV, cS.kppsGridV}, EV_matrix(:, :, ie_current), 'linear', 'linear');
            end

            % c. 并行循环遍历所有状态点
            for ik = 1:cS.nk
                % 为parfor创建临时变量，避免修改sliced variable
                val_slice = -Inf(cS.nkpps, cS.nw);
                c_slice = zeros(cS.nkpps, cS.nw);
                k_slice = zeros(cS.nkpps, cS.nw);
                cpps_slice = zeros(cS.nkpps, cS.nw);

                for ikpps = 1:cS.nkpps
                    for ie = 1:cS.nw
                        % 当前状态
                        k_state = cS.kGridV(ik);
                        k_pps_state = cS.kppsGridV(ikpps);
                        epsilon_state = paramS_age.leGridV(ie);

                        % 初始化最优值
                        best_val = -Inf;
                        best_c = cS.cFloor;
                        best_k_prime = cS.kMin;
                        best_c_pps = 0;

                        % 确定c_pps的搜索网格
                        max_cpps = 0;
                        if a_idx < cS.aR_new % 仅工作期可缴费
                            age_eff = cS.ageEffV_new(a_idx);
                            gross_labor_income = w_gross_age * age_eff * epsilon_state;
                            max_cpps = min(cS.pps_contrib_limit, gross_labor_income * cS.pps_max_contrib_frac);
                        end
                        cpps_grid = linspace(0, max(0, max_cpps), cS.npps);

                        % 遍历c_pps决策
                        for c_pps_choice = cpps_grid
                            [resources, ~, ~] = main_olg_v14_utils.HHIncome_Huggett(k_state, R_k_net_factor_age, w_gross_age, TR_total_age, b_age_val, c_pps_choice, a_idx, paramS_age, cS, epsilon_state);

                            % [核心修复] 将PPS提取的财富加入总资源
                            if a_idx >= cS.aR_new && cS.pps_active
                                pps_withdrawal_gross = k_pps_state * cS.pps_withdrawal_rate;
                                resources = resources + pps_withdrawal_gross * (1 - cS.pps_tax_rate_withdrawal);
                            end


                            k_prime_max_budget = resources - cS.cFloor * (1 + cS.tau_c);
                            if k_prime_max_budget < cS.kMin, continue; end % 如果资源不足以支付最低消费，跳过

                            k_prime_grid = linspace(cS.kMin, k_prime_max_budget, cS.nkprime);

                            % 遍历k'决策
                            temp = [];
                            for k_prime_choice = k_prime_grid
                                c_expend = resources - k_prime_choice;
                                c_choice = max(cS.cFloor, c_expend / (1 + cS.tau_c));

                                [~, util] = main_olg_v14_utils.CES_utility(c_choice, cS.sigma, cS);

                                pps_withdrawal = 0;
                                if a_idx >= cS.aR_new, pps_withdrawal = k_pps_state * cS.pps_withdrawal_rate; end
                                pps_return_factor = 1 + ((R_k_net_factor_age - 1) + cS.pps_return_rate_premium);
                                k_pps_prime = (k_pps_state + c_pps_choice - pps_withdrawal) * pps_return_factor;

                                % [核心修复] 钳位 k_prime_choice 和 k_pps_prime
                                k_prime_clamped = max(cS.kGridV(1), min(cS.kGridV(end), k_prime_choice));
                                k_pps_prime_clamped = max(cS.kppsGridV(1), min(cS.kppsGridV(end), k_pps_prime));

                                % 调用插值器
                                ev = EV_interpolants{ie}(k_prime_clamped, k_pps_prime_clamped);

                                current_val = util + cS.beta * cS.s_1yr_transitionV(a_idx) * ev;
                                % fprintf('kprime=%f, value=%f\n',k_prime_choice,current_val)
                                temp = [temp;[k_prime_choice,current_val]];

                                if current_val > best_val
                                    best_val = current_val;
                                    best_c = c_choice;
                                    best_k_prime = k_prime_choice;
                                    best_c_pps = c_pps_choice;
                                end
                            end
                        end

                        % 记录该状态点的最优解
                        val_slice(ikpps, ie) = best_val;
                        c_slice(ikpps, ie) = best_c;
                        k_slice(ikpps, ie) = best_k_prime;
                        cpps_slice(ikpps, ie) = best_c_pps;
                    end
                end

                % 将计算好的切片结果赋给主矩阵
                val_age(ik,:,:) = val_slice;
                cPol_age_q(ik,:,:) = c_slice;
                kPol_age(ik,:,:) = k_slice;
                cPpsPol_age_choice(ik,:,:) = cpps_slice;
            end
        end


        % =====================================================================
        % == [新] 2.B VFI 按年龄求解器 (高精度混合优化版) ==
        % =====================================================================
        function [cPol_age_q, kPol_age, cPpsPol_age_choice, val_age] = HHSolutionByAge_VFI_Huggett_HybridOptimizer(...
                a_idx, vPrime_kkppse_next, R_k_net_factor_age, w_gross_age, TR_total_age, b_age_val, ...
                paramS_age, cS)
            % [高精度版] 混合网格搜索与连续优化VFI求解器
            % - 对 c_pps 进行网格搜索。
            % - 对 k' (下一期资本) 使用 fminbnd 进行连续优化，以获得更高精度。
            % - 使用 'spline' 插值以更准确地估计期望未来价值。

            % 初始化输出矩阵
            val_age    = -Inf(cS.nk, cS.nkpps, cS.nw);
            cPol_age_q = zeros(cS.nk, cS.nkpps, cS.nw);
            kPol_age   = zeros(cS.nk, cS.nkpps, cS.nw);
            cPpsPol_age_choice = zeros(cS.nk, cS.nkpps, cS.nw);

            % --- [核心修正] 最后一期逻辑 (带暖光遗赠动机) ---
            if a_idx == cS.aD_new
                % 在最后一期，家庭需要在“当期消费c”和“留下遗赠b”之间分配其总财富。
                % 目标: max U(c) + phi * U(b)  s.t. c(1+tau_c) + b = total_resources

                % 1. 计算可用于分配的总资源
                [K_grid, Kpps_grid] = ndgrid(cS.kGridV, cS.kppsGridV);

                % 收入部分只包含非资本、非劳动收入（即养老金福利）
                [pension_income, ~, ~] = main_olg_v14_utils.HHIncome_Huggett(0, 0, 0, TR_total_age, b_age_val, 0, a_idx, paramS_age, cS, 0);

                % 财富部分是所有资产在赚取了当期利息后的期末价值
                wealth_at_start = K_grid + Kpps_grid * (1 - cS.pps_tax_rate_withdrawal);

                % ======================= [BUG修复的关键所在] =======================
                % 总资源 = (期初财富 * 当期回报因子) + 当期非资本收入
                total_resources = wealth_at_start * R_k_net_factor_age + pension_income;
                % =================================================================

                % 2. 求解消费和遗赠的最优分配
                % 对于CRRA效用，消费支出 c_spend = c*(1+tau_c) 和遗赠 b 之间有一个最优比例。
                phi_adj = cS.phi_bequest * (1 + cS.tau_c)^(-1);
                omega = (1 + phi_adj^(-1/cS.sigma))^(-1);

                % 计算最优消费和遗赠
                final_c_expenditure = omega .* total_resources;
                final_bequest = (1 - omega) .* total_resources;

                final_c = max(cS.cFloor, final_c_expenditure / (1 + cS.tau_c));

                % 3. 计算最终的价值函数
                [~, util_c] = main_olg_v14_utils.CES_utility(final_c, cS.sigma, cS);
                [~, util_b] = main_olg_v14_utils.CES_utility(final_bequest, cS.sigma_bequest, cS);

                final_v = util_c + cS.phi_bequest * util_b;

                % 4. 将结果赋给所有epsilon状态
                for ie = 1:cS.nw
                    cPol_age_q(:,:,ie) = final_c;
                    val_age(:,:,ie) = final_v;
                    % 在这个最优分配下，家庭不会持有任何资产到下一期，
                    % 但kPolM记录的是储蓄决策，这里可以理解为“留下的财富”，即遗赠。
                    kPol_age(:,:,ie) = final_bequest;
                    cPpsPol_age_choice(:,:,ie) = 0;
                end
                return;
            end

            % --- 非最后一期逻辑 ---

            % a. 计算期望未来价值矩阵 E[V'] (与原版相同)
            EV_matrix = zeros(cS.nk, cS.nkpps, cS.nw);
            for ie_current = 1:cS.nw
                transition_probs = paramS_age.leTrProbM(ie_current, :);
                EV_slice = pagemtimes(reshape(vPrime_kkppse_next, [cS.nk * cS.nkpps, cS.nw]), reshape(transition_probs, [cS.nw, 1]));
                EV_matrix(:, :, ie_current) = reshape(EV_slice, [cS.nk, cS.nkpps]);
            end

            % b. [改进] 创建插值器，使用 'spline' 方法提高精度，并指定外插
            EV_interpolants = cell(cS.nw, 1);
            for ie_current = 1:cS.nw
                EV_interpolants{ie_current} = griddedInterpolant({cS.kGridV, cS.kppsGridV}, EV_matrix(:, :, ie_current), 'spline', 'linear');
            end

            % c. 并行循环遍历所有状态点
            optim_options = optimset('Display', 'off', 'TolX', 1e-5); % fminbnd 优化选项

            parfor ik = 1:cS.nk
                % 为parfor创建临时变量
                val_slice = -Inf(cS.nkpps, cS.nw);
                c_slice = zeros(cS.nkpps, cS.nw);
                k_slice = zeros(cS.nkpps, cS.nw);
                cpps_slice = zeros(cS.nkpps, cS.nw);

                for ikpps = 1:cS.nkpps
                    for ie = 1:cS.nw
                        % 当前状态
                        k_state = cS.kGridV(ik);
                        k_pps_state = cS.kppsGridV(ikpps);
                        epsilon_state = paramS_age.leGridV(ie);
                        ev_interpolant = EV_interpolants{ie}; % 获取当前epsilon状态对应的插值器

                        % 初始化最优值
                        best_val_for_cpps_grid = -Inf;
                        best_c_for_cpps_grid = cS.cFloor;
                        best_k_prime_for_cpps_grid = cS.kMin;
                        best_c_pps_for_cpps_grid = 0;

                        % 确定c_pps的搜索网格
                        max_cpps = 0;
                        if a_idx < cS.aR_new
                            age_eff = cS.ageEffV_new(a_idx);
                            gross_labor_income = w_gross_age * age_eff * epsilon_state;
                            max_cpps = min(cS.pps_contrib_limit, gross_labor_income * cS.pps_max_contrib_frac);
                        end
                        cpps_grid = linspace(0, max(0, max_cpps), cS.npps);

                        % [改进] 遍历c_pps决策网格
                        for c_pps_choice = cpps_grid
                            % 1. 计算给定c_pps下的总资源
                            [resources, ~, ~] = main_olg_v14_utils.HHIncome_Huggett(k_state, R_k_net_factor_age, w_gross_age, TR_total_age, b_age_val, c_pps_choice, a_idx, paramS_age, cS, epsilon_state);

                            % [核心修复] 将PPS提取的财富加入总资源
                            if a_idx >= cS.aR_new && cS.pps_active
                                pps_withdrawal_gross = k_pps_state * cS.pps_withdrawal_rate;
                                resources = resources + pps_withdrawal_gross * (1 - cS.pps_tax_rate_withdrawal);
                            end

                            k_prime_max_budget = resources - cS.cFloor * (1 + cS.tau_c);
                            if k_prime_max_budget < cS.kMin, continue; end % 资源不足

                            % 2. 定义目标函数，用于对 k' 进行优化
                            % 2. 定义目标函数
                            % [重要修改] 现在调用一个独立的静态方法
                            objective_func = @(k_prime) main_olg_v14_utils.objective_for_k_prime_private(...
                                k_prime, resources, k_pps_state, c_pps_choice, ...
                                R_k_net_factor_age, a_idx, ev_interpolant, cS);

                            % 3. [核心改进] 使用fminbnd
                            [k_prime_opt, neg_val_opt] = fminbnd(objective_func, cS.kMin, k_prime_max_budget, optim_options);

                            current_max_val = -neg_val_opt;

                            % 4. 更新最优解
                            if current_max_val > best_val_for_cpps_grid
                                best_val_for_cpps_grid = current_max_val;
                                best_k_prime_for_cpps_grid = k_prime_opt;
                                best_c_pps_for_cpps_grid = c_pps_choice;

                                % 根据最优 k' 和 c_pps 计算对应的 c
                                c_expend = resources - best_k_prime_for_cpps_grid;
                                best_c_for_cpps_grid = max(cS.cFloor, c_expend / (1 + cS.tau_c));
                            end
                        end % 结束对 c_pps_grid 的循环

                        % 记录该状态点的最优解
                        val_slice(ikpps, ie) = best_val_for_cpps_grid;
                        c_slice(ikpps, ie) = best_c_for_cpps_grid;
                        k_slice(ikpps, ie) = best_k_prime_for_cpps_grid;
                        cpps_slice(ikpps, ie) = best_c_pps_for_cpps_grid;
                    end
                end

                % 将计算好的切片结果赋给主矩阵
                val_age(ik,:,:) = val_slice;
                cPol_age_q(ik,:,:) = c_slice;
                kPol_age(ik,:,:) = k_slice;
                cPpsPol_age_choice(ik,:,:) = cpps_slice;
            end
        end


        % =========================================================================
        % == 在 utils 文件中，替换整个 HHSolutionByAge_VFI_Huggett_VectorizedGrid 函数 ==
        % =========================================================================

        % =========================================================================
        % == 在 utils 文件中，替换整个 HHSolutionByAge_VFI_Huggett_VectorizedGrid 函数 ==
        % =========================================================================
        % =========================================================================
        % == 在 utils 文件中，替换 HHSolutionByAge..._VectorizedGrid (最终版) ==
        % =========================================================================
        % =========================================================================
        % == 在 utils 文件中，修改 HHSolutionByAge_VFI_Huggett_VectorizedGrid ==
        % =========================================================================
        function [cPol_age_q, kPol_age, cPpsPol_age_choice, val_age] = HHSolutionByAge_VFI_Huggett_VectorizedGrid(...
                a_idx, vPrime_kkppse_next, R_k_mkt_factor_age, w_gross_age, TR_total_age, b_age_val, r_mkt_annual_age, ... % <--- 新增参数
                paramS_age, cS)

            % ... (函数开头的初始化和最后一期逻辑保持不变) ...
            val_age = -Inf(cS.nk, cS.nkpps, cS.nw); cPol_age_q = zeros(cS.nk, cS.nkpps, cS.nw); kPol_age = zeros(cS.nk, cS.nkpps, cS.nw); cPpsPol_age_choice = zeros(cS.nk, cS.nkpps, cS.nw);
            if a_idx == cS.aD_new, [K_grid, Kpps_grid] = ndgrid(cS.kGridV, cS.kppsGridV); [period_net_income, ~] = main_olg_v14_utils.HHIncome_Huggett(R_k_mkt_factor_age, 0, TR_total_age, b_age_val, 0, a_idx, paramS_age, cS, 0, K_grid); total_resources = K_grid + Kpps_grid + period_net_income; if cS.phi_bequest > 1e-9, phi_adj = cS.phi_bequest * (1 + cS.tau_c)^(-1); omega = (1 + phi_adj^(-1/cS.sigma))^(-1); final_c_expenditure = omega .* total_resources; final_bequest = (1 - omega) .* total_resources; else, final_c_expenditure = total_resources; final_bequest = zeros(size(total_resources)); end, final_c = max(cS.cFloor, final_c_expenditure / (1 + cS.tau_c)); [~, util_c] = main_olg_v14_utils.CES_utility(final_c, cS.sigma, cS); [~, util_b] = main_olg_v14_utils.CES_utility(final_bequest, cS.sigma_bequest, cS); final_v = util_c + cS.phi_bequest * util_b; for ie = 1:cS.nw, cPol_age_q(:,:,ie) = final_c; val_age(:,:,ie) = final_v; kPol_age(:,:,ie) = final_bequest; cPpsPol_age_choice(:,:,ie) = 0; end, return; end
            EV_matrix = zeros(cS.nk, cS.nkpps, cS.nw); for ie_current = 1:cS.nw, transition_probs = paramS_age.leTrProbM(ie_current, :); vPrime_reshaped = reshape(vPrime_kkppse_next, [cS.nk * cS.nkpps, cS.nw]); EV_slice = vPrime_reshaped * transition_probs'; EV_matrix(:, :, ie_current) = reshape(EV_slice, [cS.nk, cS.nkpps]); end, EV_interpolants = cell(cS.nw, 1); is_pps_disabled = (cS.nkpps == 1); for ie_current = 1:cS.nw, if is_pps_disabled, EV_vector_k = squeeze(EV_matrix(:, 1, ie_current)); EV_interpolants{ie_current} = griddedInterpolant(cS.kGridV, EV_vector_k, 'linear', 'linear'); else, EV_interpolants{ie_current} = griddedInterpolant({cS.kGridV, cS.kppsGridV}, EV_matrix(:, :, ie_current), 'linear', 'linear'); end, end
            prop_k_prime_grid = linspace(0, 1, cS.nkprime)'; if cS.npps == 1, prop_cpps_grid = 0; else, prop_cpps_grid = linspace(0, cS.pps_max_contrib_frac, cS.npps)'; end

            for ie = 1:cS.nw
                % ... (大部分循环内部代码保持不变) ...
                epsilon_state = paramS_age.leGridV(ie); ev_interpolant = EV_interpolants{ie};
                k_state_4D = repmat(reshape(cS.kGridV, [cS.nk, 1, 1, 1]), [1, cS.nkpps, cS.nkprime, cS.npps]); kpps_state_4D = repmat(reshape(cS.kppsGridV, [1, cS.nkpps, 1, 1]), [cS.nk, 1, cS.nkprime, cS.npps]); prop_k_prime_4D = repmat(reshape(prop_k_prime_grid, [1, 1, cS.nkprime, 1]), [cS.nk, cS.nkpps, 1, cS.npps]); prop_cpps_4D = repmat(reshape(prop_cpps_grid, [1, 1, 1, cS.npps]), [cS.nk, cS.nkpps, cS.nkprime, 1]);
                actual_cpps_4D_input = zeros(size(prop_cpps_4D)); if a_idx <= cS.aR_new && cS.pps_active, age_efficiency = cS.ageEffV_new(a_idx); gross_labor_income_4D = w_gross_age * age_efficiency * epsilon_state; actual_cpps_4D_input = min(cS.pps_contrib_limit, gross_labor_income_4D .* prop_cpps_4D); end
                [period_net_income_4D, actual_cpps_4D] = main_olg_v14_utils.HHIncome_Huggett(R_k_mkt_factor_age, w_gross_age, TR_total_age, b_age_val, actual_cpps_4D_input, a_idx, paramS_age, cS, epsilon_state, k_state_4D);
                pps_withdrawal_val = 0; if a_idx >= cS.aR_new && cS.pps_active, pps_withdrawal_val = kpps_state_4D .* cS.pps_withdrawal_rate * (1 - cS.pps_tax_rate_withdrawal); end
                resources_4D = k_state_4D + kpps_state_4D + period_net_income_4D - actual_cpps_4D + pps_withdrawal_val;
                c_floor_spending = cS.cFloor * (1 + cS.tau_c); resources_above_floor_4D = max(0, resources_4D - c_floor_spending); actual_k_prime_4D = resources_above_floor_4D .* prop_k_prime_4D; c_expend_4D = resources_4D - actual_k_prime_4D; c_choice_4D = max(cS.cFloor, c_expend_4D / (1 + cS.tau_c));
                [~, util_4D] = main_olg_v14_utils.CES_utility(c_choice_4D, cS.sigma, cS);
                pps_withdrawal_for_evolution = 0; if a_idx >= cS.aR_new && cS.pps_active, pps_withdrawal_for_evolution = kpps_state_4D .* cS.pps_withdrawal_rate; end

                % [核心修正] 使用传递进来的 r_mkt_annual_age
                pps_return_factor_5y = (1 + r_mkt_annual_age + cS.pps_return_rate_premium)^cS.time_Step;

                k_pps_prime_4D = (kpps_state_4D + actual_cpps_4D - pps_withdrawal_for_evolution) .* pps_return_factor_5y;

                % ... (后续代码保持不变) ...
                k_prime_clamped = max(cS.kGridV(1), min(cS.kGridV(end), actual_k_prime_4D)); k_pps_prime_clamped = max(cS.kppsGridV(1), min(cS.kppsGridV(end), k_pps_prime_4D));
                if is_pps_disabled, ev_mat = ev_interpolant(k_prime_clamped); else, ev_mat = ev_interpolant(k_prime_clamped, k_pps_prime_clamped); end
                val_grid = util_4D + cS.beta * cS.s_1yr_transitionV(a_idx) * ev_mat; val_grid(c_expend_4D < 0) = -Inf;
                [val_max_k, idx_k_prime] = max(val_grid, [], 3); [val_max_kc, idx_cpps] = max(val_max_k, [], 4); val_age(:,:,ie) = squeeze(val_max_kc);
                [I, J] = ndgrid(1:cS.nk, 1:cS.nkpps); temp = squeeze(idx_cpps);
                if size(idx_k_prime, 4) == 1, linear_idx_k = sub2ind(size(squeeze(idx_k_prime(:,:,:,1))), I(:), J(:),temp(:)); else, linear_idx_k = sub2ind(size(idx_k_prime), I(:), J(:), ones(numel(I),1), temp(:)); end
                final_k_prime_prop_idx = reshape(idx_k_prime(linear_idx_k), [cS.nk, cS.nkpps]); best_prop_k_prime = prop_k_prime_grid(final_k_prime_prop_idx);
                final_cpps_prop_idx = squeeze(idx_cpps); best_prop_cpps = prop_cpps_grid(final_cpps_prop_idx);
                [k_state_2D, kpps_state_2D] = ndgrid(cS.kGridV, cS.kppsGridV);
                best_actual_cpps = zeros(size(best_prop_cpps));
                if a_idx <= cS.aR_new && cS.pps_active, age_efficiency = cS.ageEffV_new(a_idx); gross_labor_income_2D = w_gross_age * age_efficiency * epsilon_state; best_actual_cpps = min(cS.pps_contrib_limit, gross_labor_income_2D .* best_prop_cpps); end
                [period_net_income_final, ~] = main_olg_v14_utils.HHIncome_Huggett(R_k_mkt_factor_age, w_gross_age, TR_total_age, b_age_val, best_actual_cpps, a_idx, paramS_age, cS, epsilon_state, k_state_2D);
                pps_withdrawal_final = 0; if a_idx >= cS.aR_new && cS.pps_active, pps_withdrawal_final = kpps_state_2D .* cS.pps_withdrawal_rate * (1 - cS.pps_tax_rate_withdrawal); end
                resources_final = k_state_2D + kpps_state_2D + period_net_income_final - best_actual_cpps + pps_withdrawal_final;
                resources_above_floor_final = max(0, resources_final - c_floor_spending); best_k_prime = resources_above_floor_final .* best_prop_k_prime; c_expend_final = resources_final - best_k_prime;
                cPol_age_q(:,:,ie) = max(cS.cFloor, c_expend_final / (1 + cS.tau_c)); kPol_age(:,:,ie) = best_k_prime; cPpsPol_age_choice(:,:,ie) = best_actual_cpps;
            end
        end


        function [kHistM_out, kPpsHistM_out, cHistM_out, cppsHistM_out] = HHSimulation_olgm(...
                kPolM, cPpsPolM_choice, cPolM_consump, eIdxM_group, ...
                R_k_net, w, TR, bV_payg, paramS_sim, cS_sim)
            % [最终对齐版] - 与Python的HHSimulation_olgm_rl和VFI模拟器在物理过程上完全一致。
            % - 使用 griddedInterpolant 模拟家庭基于VFI策略的生命周期路径。
            % - 核心是按年龄组进行模拟，并确保所有状态演化逻辑正确。

            nSim = size(eIdxM_group, 1);
            aD = cS_sim.aD_new;

            % 初始化历史记录矩阵
            kHistM_out = zeros(nSim, aD);
            kPpsHistM_out = zeros(nSim, aD);
            cHistM_out = zeros(nSim, aD);
            cppsHistM_out = zeros(nSim, aD);

            % --- 1. 创建策略函数的插值器 ---
            kPolInterp = cell(cS_sim.nw, aD);
            cPpsPolInterp = cell(cS_sim.nw, aD);
            cPolInterp = cell(cS_sim.nw, aD);

            for ia = 1:aD
                for ie = 1:cS_sim.nw
                    % 根据资产网格维度选择合适的插值方法
                    if cS_sim.nk > 1 && cS_sim.nkpps > 1
                        kPolInterp{ie,ia} = griddedInterpolant({cS_sim.kGridV, cS_sim.kppsGridV}, squeeze(kPolM(:,:,ie,ia)), 'linear', 'linear');
                        cPpsPolInterp{ie,ia} = griddedInterpolant({cS_sim.kGridV, cS_sim.kppsGridV}, squeeze(cPpsPolM_choice(:,:,ie,ia)), 'linear', 'linear');
                        cPolInterp{ie,ia} = griddedInterpolant({cS_sim.kGridV, cS_sim.kppsGridV}, squeeze(cPolM_consump(:,:,ie,ia)), 'linear', 'linear');
                    elseif cS_sim.nk > 1
                        kPolInterp{ie,ia} = griddedInterpolant(cS_sim.kGridV, squeeze(kPolM(:,1,ie,ia)), 'linear', 'linear');
                        cPpsPolInterp{ie,ia} = griddedInterpolant(cS_sim.kGridV, squeeze(cPpsPolM_choice(:,1,ie,ia)), 'linear', 'linear');
                        cPolInterp{ie,ia} = griddedInterpolant(cS_sim.kGridV, squeeze(cPolM_consump(:,1,ie,ia)), 'linear', 'linear');
                    elseif cS_sim.nkpps > 1 % 只有kpps网格
                        kPolInterp{ie,ia} = griddedInterpolant(cS_sim.kppsGridV, squeeze(kPolM(1,:,ie,ia))', 'linear', 'linear');
                        cPpsPolInterp{ie,ia} = griddedInterpolant(cS_sim.kppsGridV, squeeze(cPpsPolM_choice(1,:,ie,ia))', 'linear', 'linear');
                        cPolInterp{ie,ia} = griddedInterpolant(cS_sim.kppsGridV, squeeze(cPolM_consump(1,:,ie,ia))', 'linear', 'linear');
                    else % nk=1, nkpps=1 的标量情况
                        kPolInterp{ie,ia} = @(x,y) squeeze(kPolM(1,1,ie,ia));
                        cPpsPolInterp{ie,ia} = @(x,y) squeeze(cPpsPolM_choice(1,1,ie,ia));
                        cPolInterp{ie,ia} = @(x,y) squeeze(cPolM_consump(1,1,ie,ia));
                    end
                end
            end

            % --- 2. 初始化状态和参数 ---
            pps_return_factor = 1 + ((R_k_net - 1) + cS_sim.pps_return_rate_premium);
            k_next = zeros(nSim, 1);
            k_pps_next = zeros(nSim, 1);

            % --- 3. 按年龄组进行前向模拟 ---
            for a_idx = 1:aD
                % 获取当前状态
                k_now = k_next;
                k_pps_now = k_pps_next;
                kHistM_out(:, a_idx) = k_now;
                kPpsHistM_out(:, a_idx) = k_pps_now;

                % 初始化当期决策向量
                k_prime_decision = zeros(nSim, 1);
                cpps_decision = zeros(nSim, 1);
                c_decision = zeros(nSim, 1);

                % --- 4. 根据效率冲击状态，查询策略并获取决策 ---
                for ie = 1:cS_sim.nw
                    % 找到当前效率状态对应的所有个体
                    idx_sim = find(eIdxM_group(:, a_idx) == ie);
                    if isempty(idx_sim), continue; end

                    k_now_e = k_now(idx_sim);
                    k_pps_now_e = k_pps_now(idx_sim);

                    % 验证代码
                    if any(k_now_e > cS_sim.kMax) || any(k_pps_now_e > cS_sim.kppsMax)
                        fprintf('警告: a_idx=%d, ie=%d, 有 %d 个体超出网格范围！\n', ...
                            a_idx, ie, sum(k_now_e > cS_sim.kMax | k_pps_now_e > cS_sim.kppsMax));
                    end

                    % 使用插值器获取决策
                    if cS_sim.nk > 1 && cS_sim.nkpps > 1
                        k_prime_decision(idx_sim) = kPolInterp{ie, a_idx}(k_now_e, k_pps_now_e);
                        cpps_decision(idx_sim) = cPpsPolInterp{ie, a_idx}(k_now_e, k_pps_now_e);
                        c_decision(idx_sim) = cPolInterp{ie, a_idx}(k_now_e, k_pps_now_e);
                    elseif cS_sim.nk > 1
                        k_prime_decision(idx_sim) = kPolInterp{ie, a_idx}(k_now_e);
                        cpps_decision(idx_sim) = cPpsPolInterp{ie, a_idx}(k_now_e);
                        c_decision(idx_sim) = cPolInterp{ie, a_idx}(k_now_e);
                    elseif cS_sim.nkpps > 1
                        k_prime_decision(idx_sim) = kPolInterp{ie, a_idx}(k_pps_now_e);
                        cpps_decision(idx_sim) = cPpsPolInterp{ie, a_idx}(k_pps_now_e);
                        c_decision(idx_sim) = cPolInterp{ie, a_idx}(k_pps_now_e);
                    else % 标量情况
                        k_prime_decision(idx_sim) = kPolInterp{ie, a_idx}(k_now_e, k_pps_now_e);
                        cpps_decision(idx_sim) = cPpsPolInterp{ie, a_idx}(k_now_e, k_pps_now_e);
                        c_decision(idx_sim) = cPolInterp{ie, a_idx}(k_now_e, k_pps_now_e);
                    end
                end

                % --- 5. 记录当期决策并演化到下一期状态 ---
                cHistM_out(:, a_idx) = max(cS_sim.cFloor, c_decision);
                cppsHistM_out(:, a_idx) = max(0, cpps_decision);

                if a_idx < aD
                    % 演化非PPS资产
                    k_next = max(cS_sim.kMin, min(cS_sim.kMax, k_prime_decision));

                    % 演化PPS资产
                    pps_withdrawal = 0;
                    if a_idx >= cS_sim.aR_new && cS_sim.pps_active
                        pps_withdrawal = k_pps_now * cS_sim.pps_withdrawal_rate;
                    end
                    k_pps_next_unclamped = (k_pps_now + cppsHistM_out(:, a_idx) - pps_withdrawal) * pps_return_factor;
                    k_pps_next = max(cS_sim.kppsMin, min(cS_sim.kppsMax, k_pps_next_unclamped));
                end
            end
        end


        % =========================================================================
        % =========================================================================
        % == 在 utils 文件中，替换 simulate_private_capital_forward (最终版) ==
        % =========================================================================
        % =========================================================================
        % == 在 utils 文件中，修改 simulate_private_capital_forward ==
        % =========================================================================
        % =========================================================================
        % == 在 utils 文件中，替换 simulate_private_capital_forward (最终、完整、正确版) ==
        % =========================================================================
        function [K_pvt_next, K_pps_next, C_t, Total_Cpps_t, Total_PpsWithdrawalTax_t, total_accidental_bequest] = simulate_private_capital_forward(M_t, Z_t_norm, cS, paramS, eIdxM, TR_total)
            % [v14 最终审计修正版 v3]
            % - 重新加入了输出变量的默认值初始化，以修复 "not assigned" 错误。
            % - 确保所有接口和参数传递都正确无误。

            % --- [关键修正] 初始化所有输出变量，特别是那些在条件块中赋值的 ---
            Total_Cpps_t = 0;
            Total_PpsWithdrawalTax_t = 0;
            total_accidental_bequest = 0;

            % 1. 准备VFI求解所需的参数
            R_k_mkt_factor_hh = M_t.R_mkt_factor_5y;
            w_gross = M_t.w_t;
            bV_payg = zeros(1, cS.aD_new);
            bV_payg(cS.aR_new+1:end) = M_t.b_t;
            paramS_vfi = paramS;

            % --- 严格的PPS开关 ---
            cS_vfi = cS;
            if cS.sim_years(M_t.current_t) < cS.pps_activation_year
                cS_vfi.pps_active = false;
            else
                cS_vfi.pps_active = true;
            end
            if ~cS_vfi.pps_active
                cS_vfi.nkpps = 1; cS_vfi.npps = 1;
                cS_vfi = main_olg_v14_utils.generateGrids(cS_vfi);
            end
            cS_vfi.theta_t = cS.theta_path(M_t.current_t);

            % 2. 调用VFI求解器 (传递所有需要的参数)
            [cPolM, kPolM, cPpsPolM, ~] = main_olg_v14_utils.HHSolution_VFI_Huggett(...
                R_k_mkt_factor_hh, w_gross, TR_total, bV_payg, M_t.r_mkt_annual, ...
                paramS_vfi, cS_vfi, 'vectorized_grid');

            % 3. 使用策略，模拟所有家庭在 t 期的决策
            [kHistM, kPpsHistM, cHistM, cppsHistM] = main_olg_v14_utils.HHSimulation_olgm(...
                kPolM, cPpsPolM, cPolM, eIdxM, ...
                M_t.R_net_factor_5y, w_gross, TR_total, bV_payg, paramS_vfi, cS_vfi);

            % 4. 聚合得到下一期的私人资本 K_pvt_{t+1} 和 K_pps_{t+1}
            k_prime_paths = kHistM(:, 2:end);
            mean_k_prime_by_age = mean(k_prime_paths, 1);
            weights = Z_t_norm(1:cS.aD_new-1);
            K_pvt_next = mean_k_prime_by_age * weights;

            k_pps_prime_paths = kPpsHistM(:, 2:end);
            mean_k_pps_prime_by_age = mean(k_pps_prime_paths, 1);
            K_pps_next = mean_k_pps_prime_by_age * weights;

            % 5. 聚合得到当期总消费 C_t
            mean_c_by_age = mean(cHistM, 1);
            C_t = mean_c_by_age * Z_t_norm;

            % 6. 计算本期产生的、供下一期使用的意外遗赠总额
            mass_died = Z_t_norm .* (1 - cS.s_1yr_transitionV);
            k_with_return = kHistM .* M_t.R_net_factor_5y;
            pps_return_factor_5y = (1 + M_t.r_mkt_annual + cS.pps_return_rate_premium)^cS.time_Step;
            k_pps_with_return = kPpsHistM .* pps_return_factor_5y;
            avg_wealth_by_age = mean(k_with_return + k_pps_with_return * (1-cS.pps_tax_rate_withdrawal), 1);
            total_accidental_bequest = sum(mass_died' .* avg_wealth_by_age);

            % 7. 聚合其他流量
            Total_Cpps_t = mean(cppsHistM, 1) * Z_t_norm;

            % [关键] 只有在PPS激活时才重新计算，否则使用默认值0
            if cS_vfi.pps_active && cS.pps_withdrawal_rate > 0
                mean_k_pps_by_age = mean(kPpsHistM, 1);
                retired_ages_idx = (cS.aR_new + 1):cS.aD_new;
                k_pps_retirees_avg = mean_k_pps_by_age(retired_ages_idx);
                mass_retirees_vec = Z_t_norm(retired_ages_idx);
                total_pps_capital_of_retirees = k_pps_retirees_avg * mass_retirees_vec;
                Total_PpsWithdrawalTax_t = total_pps_capital_of_retirees * cS.pps_withdrawal_rate * cS.pps_tax_rate_withdrawal;
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