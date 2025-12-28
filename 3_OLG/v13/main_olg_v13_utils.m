% =========================================================================
% == 在 main_olg_v13_utils.m 文件中，请使用以下代码块 ==get_p
% == 替换或添加相应的函数。                            ==
% =========================================================================

classdef main_olg_v13_utils

    % ... 这里是您从v13复制过来的所有静态方法 ...

    methods (Static)

        % --- [v13 新增] 从Excel加载和处理外生路径 ---
        function [Z_path, A_path, T_sim] = load_exogenous_paths(cS)
            fprintf('--- [v13] 加载和处理外生数据路径 ---\n');

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

        % --- [v13 修改] 价格计算函数，现在接收 A_t ---
        function [R_market_gross_factor, MPL_gross, Y_t_annual] = HHPrices_Huggett(K, L, A_t, cS) % 输出年化Y
            if K <= 0, K = 1e-6; end
            if L <= 0, L = 1e-6; end

            % [核心修正] Y的计算应该基于年化TFP
            % A_path是5年期累积增长的结果，我们需要反算出年化增长率
            g_5y = A_t / 1.0; % 假设 A(0)=1
            g_annual = (g_5y)^(1/cS.time_Step) - 1;
            A_t_annual = (1 + g_annual); % 这是一个近似，更准确的应该是A(t)的几何平均
            % 一个更简单、更稳健的方法是直接将A_t视为一个抽象的生产力水平

            % [更稳健的修正] Y_t_annual 计算年化产出。
            % A_t 是一个水平，L是年化流量，所以 K 也应被视为一个“服务流量”
            Y_t_annual = A_t * (K^cS.alpha) * (L^(1-cS.alpha));

            MPK_gross = cS.alpha * Y_t_annual / K;
            MPL_gross = (1-cS.alpha) * Y_t_annual / L;
            % [核心修正] R_market_gross_factor 是年化回报因子
            r_annual_gross = MPK_gross - (1-(1-cS.ddk)^(1/cS.time_Step)); % 使用年化折旧率
            R_market_gross_factor = 1 + r_annual_gross;
        end
        % --- [v13 修改] 价格和政策计算的封装函数 ---
        % --- [v13 修改] 价格和政策计算的封装函数 ---
        function M_t = get_prices_and_policy_at_t(t, K_pvt_t, B_p_t, K_pps_t, B_g_t, Z_t_norm, A_t, cS, paramS, eIdxM) % 新增 B_g_t 输入
            % a. [核心修正] 计算用于生产的物理总资本 K_physical
            K_physical_t = K_pvt_t + K_pps_t + B_p_t - B_g_t;

            if K_physical_t <= 0
                % 如果物理资本为负或零，这是一个严重问题，通常由不合理的参数导致
                % 比如政府债务 B_g 过高，超过了所有储蓄。
                % 返回一个错误或极差的结果，让上层迭代调整。
                warning('计算出的物理资本 K_physical_t (%.2f) 为负或零，迭代可能失败。', K_physical_t);
                K_physical_t = 1e-8; % 设置一个很小的正数以避免计算崩溃
            end

            % b. 计算劳动供给
            paramS_t = paramS;
            paramS_t.ageMassV = Z_t_norm;
            [~, L_t] = main_olg_v13_utils.LaborSupply_Huggett(eIdxM, cS, paramS_t, Z_t_norm);

            % c. 获取价格和【年化】总产出
            [R_mkt_factor_annual, w_t, Y_t_annual] = main_olg_v13_utils.HHPrices_Huggett(K_physical_t, L_t, A_t, cS);
            r_mkt_t_annual = R_mkt_factor_annual - 1;
            r_net_t_annual = r_mkt_t_annual * (1 - cS.tau_k);

            % d. 计算养老金福利 (基于固定的替代率)
            mass_workers_t = sum(Z_t_norm(1:cS.aR_new));
            if mass_workers_t > 0
                avg_wage_t = w_t * L_t / mass_workers_t;
            else
                avg_wage_t = 0; % 避免除以零
            end
            b_t = cS.rho_prime_payg * avg_wage_t;

            % e. 打包所有当期宏观变量
            M_t = struct();
            M_t.current_t = t;
            M_t.K_physical_t = K_physical_t;
            M_t.L_t = L_t;
            M_t.Y_t = Y_t_annual * cS.time_Step; % [核心] Y_t 存储为5年期总量，以匹配C,I,G
            M_t.w_t = w_t;                       % w 是年化工资率
            M_t.r_mkt_t = r_mkt_t_annual;        % r 是年化利率
            M_t.r_net_t = r_net_t_annual;        % r_net 是年化净利率
            M_t.b_t = b_t;                       % b 是年化养老金

            % [新] 将原始资产存量也存入，便于后续演化函数调用
            M_t.K_pvt_t = K_pvt_t;
            M_t.K_pps_t = K_pps_t;
            M_t.B_p_t = B_p_t;
            M_t.B_g_t = B_g_t;
            M_t.Z_t_norm = Z_t_norm;
            M_t.paramS = paramS_t; % paramS_t 是函数内部使用的变量
            M_t.eIdxM = eIdxM;
        end

        function B_p_next = update_pension_fund(B_p_t, M_t, Z_t_norm, cS)
            % [v13 结构性修正] 更新养老金基金，正确处理5年期流量与存量

            theta_t = cS.theta_path(M_t.current_t);

            % a. 计算年化流量 (Annual Flow)
            PensionRevenue_annual = theta_t * M_t.w_t * M_t.L_t;
            mass_retirees_t = sum(Z_t_norm(cS.aR_new+1:end));
            PensionOutlay_annual = M_t.b_t * mass_retirees_t;
            NetSurplus_annual = PensionRevenue_annual - PensionOutlay_annual;

            % b. 计算5年期回报因子和流量累积因子
            r_a = M_t.r_mkt_t; % 年化市场利率
            R_5y = (1 + r_a)^cS.time_Step;

            % 流量累积因子: ( (1+r)^5 - 1 ) / r
            if abs(r_a) > 1e-6
                FlowAcc_5y = (R_5y - 1) / r_a;
            else % 利率接近0时的特殊处理，避免除以0
                FlowAcc_5y = cS.time_Step;
            end

            % c. 应用正确的多年期演化方程
            B_p_next = B_p_t * R_5y + NetSurplus_annual * FlowAcc_5y;
            B_p_next = max(0, B_p_next);
        end

        % --- START OF NEW FUNCTION in main_olg_v13_utils.m ---

        % --- START OF NEW FUNCTION in main_olg_v13_utils.m ---

        % --- START OF CORRECTED FUNCTION in main_olg_v13_utils.m ---

        % --- START OF CORRECTED FUNCTION in main_olg_v13_utils.m ---

        % =========================================================================
        % == 在 main_olg_v13_utils.m 中，替换整个 calcaulte_theta_payg_path 函数 ==
        % =========================================================================

        function [cS] = calcaulte_theta_payg_path(cS, graph_flag)
            % [v13.4 - 已修正插值逻辑]
            % 核心修正：从插值【绝对人数】改为直接插值【覆盖率】，以保证路径的平滑和目标的准确达成。

            fprintf('正在构建基于【覆盖率插值】的有效养老金缴费率路径 (theta_path)...\n');

            % --- 1. 定义关键参数 ---
            theta_urban_employee_effective = 0.20; % 城镇职工体系的名义缴费率
            theta_resident_effective = 0.03;       % 城乡居民体系的名义缴费率 (近似值)

            % [政策目标] 定义未来的目标覆盖率
            coverage_urban_final = 0.80;    % 目标：城镇职工覆盖率达到劳动人口的60%
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

        % --- [v13 修改] 政府债务演化，正确处理5年期利率和流量 ---
        % =========================================================================
        % == 在您的 utils 文件中，替换整个 update_gov_debt 函数 ==
        % =========================================================================

        function [B_g_next, G_t] = update_gov_debt(B_g_t, C_t, M_t, Total_Cpps_t, K_pps_t, Total_PpsWithdrawalTax_t, cS)
            % [v13 结构性修正 - 字段名修复] 更新政府债务，正确处理5年期流量与存量

            % a. 计算年化流量 (Annual Flow)
            %    注意：M_t中的Y,C等是5年期总量，需要先转换为年均量
            %    [注] 您的代码可能已将Y,C等处理为年化，这里假设w*L, r*K等是年化流量
            Y_annual = M_t.Y_t / cS.time_Step;
            C_annual = C_t / cS.time_Step;
            Total_Cpps_annual = Total_Cpps_t / cS.time_Step;
            Total_PpsWithdrawalTax_annual = Total_PpsWithdrawalTax_t / cS.time_Step;

            G_annual = cS.gov_exp_frac_Y * Y_annual;
            G_t = G_annual * cS.time_Step; % 输出5年期总支出

            GrossLaborIncome_annual = (M_t.w_t * M_t.L_t);
            TaxableLaborIncome_annual = GrossLaborIncome_annual - Total_Cpps_annual;
            LaborTaxRevenue_annual = cS.tau_l * max(0, TaxableLaborIncome_annual);

            % --- [核心修正] ---
            % 使用新的字段名 K_physical_t 替代旧的 K_total_t
            % TaxableCapitalStock_t 是产生应税资本收入的资本存量
            TaxableCapitalStock_t = M_t.K_physical_t - K_pps_t;
            CapitalTaxRevenue_annual = cS.tau_k * M_t.r_mkt_t * TaxableCapitalStock_t;
            % --------------------

            ConsumptionTaxRevenue_annual = cS.tau_c * C_annual;

            TaxRevenue_annual = LaborTaxRevenue_annual + CapitalTaxRevenue_annual + ConsumptionTaxRevenue_annual + Total_PpsWithdrawalTax_annual;

            PrimaryDeficit_annual = G_annual - TaxRevenue_annual;

            % b. 计算5年期回报因子和流量累积因子
            r_a = M_t.r_mkt_t; % 年化市场利率
            R_5y = (1 + r_a)^cS.time_Step;
            if abs(r_a) > 1e-6
                FlowAcc_5y = (R_5y - 1) / r_a;
            else
                FlowAcc_5y = cS.time_Step;
            end

            % c. 应用正确的多年期演化方程
            B_g_next = B_g_t * R_5y + PrimaryDeficit_annual * FlowAcc_5y;
        end
        
        function cS = ParameterValues_HuggettStyle()
            % [VFI最终对齐版 v3 - 修正生存概率]

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
            
            % --- [核心修正] 计算正确的5年期生存概率 ---
            % 将旧变量重命名，以示区别和清晰
            cS.s_5yr_transitionV = zeros(cS.aD_new, 1);
            for a = 1:(cS.aD_new - 1)
                % 获取当前5岁年龄组中包含的所有原始年度年龄的索引
                age_indices_in_group = cS.physAgeMap{a};
                
                % 5年期生存率是该组内所有年度生存率的乘积
                cS.s_5yr_transitionV(a) = prod(s_orig(age_indices_in_group));
            end
            % 最后一个年龄组的生存率为0，因为模型在那里结束
            cS.s_5yr_transitionV(cS.aD_new) = 0;
            
            % 我们需要将旧的变量名 s_1yr_transitionV 替换为 s_5yr_transitionV
            % 在所有调用它的地方。最主要的是VFI求解器。

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
                    Z_next(a) = Z_current(a-1) * cS.s_5yr_transitionV(a-1);
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
            % [v13.3辅助] 仅根据总资本和人口，计算当期的要素价格和宏观总量。

            paramS_t.ageMassV = Z_t_norm;
            [~, L_t] = main_olg_v13_utils.LaborSupply_Huggett(eIdxM, cS, paramS_t, Z_t_norm);

            [R_mkt_factor, w_t] = main_olg_v13_utils.HHPrices_Huggett(K_total_t, L_t, cS);
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
        % --- 在 main_olg_v13_utils.m 中，替换 check_gbc_residual 函数 ---

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

        % --- 在 main_olg_v13_utils.m 的 methods (Static) 块中新增此函数 ---

        % --- 在 main_olg_v13_utils.m 的 methods (Static) 块中新增此函数 ---
        % --- 您可以删除或注释掉旧的 solve_steady_state 函数 ---

        % --- 在 main_olg_v13_utils.m 的 methods (Static) 块中新增此函数 ---
        % --- 您可以删除或注释掉旧的 solve_steady_state_endo_rho 函数 ---

        % --- 在 main_olg_v13_utils.m 的 methods (Static) 块中新增此函数 ---
        % --- 您可以删除或注释掉所有旧的 solve_... 函数 ---

        function [ss, eq_found] = solve_steady_state_with_fund(Z_ss_norm, B_p_Y_ratio_target, cS, paramS, eIdxM)
            % [v13.3核心, 已修正K_pps聚合] 求解包含目标养老金基金的稳态均衡。
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
                M_ss_prices = main_olg_v13_utils.get_prices_at_t(K_guess, Z_ss_norm, cS, paramS_ss, eIdxM);



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

                [K_pvt_model, K_pps_model, C_ss_model, ~] = main_olg_v13_utils.simulate_private_capital_forward(M_ss, Z_ss_norm, cS, paramS_ss, eIdxM);

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
            % [v13.2 - 精确版GBC]
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

        % --- 在 main_olg_v13_utils.m 中，替换整个 HHIncome_Huggett 函数 ---

        function [resources_for_c_and_k_prime, actual_pps_contribution_expenditure, non_capital_income] = HHIncome_Huggett(...
                k_now_val, R_k_net_factor, w_gross, TR_total, b_payg_val, ...
                c_pps_chosen, a_idx, paramS_hh, cS, epsilon_val)
            % [v13.2 - 精确税收递延版]
            % 精确模拟了 PAYG税 和 一般劳动税 的不同税基，
            % 以及个人养老金 (PPS) 对一般劳动税的税收递延效应。

            actual_pps_contribution_expenditure = 0;
            non_capital_income = 0;

            % --- 1. 计算非资本收入 (工作期) ---
            if a_idx <= cS.aR_new
                % a. 基础变量
                age_efficiency = cS.ageEffV_new(a_idx);
                labor_income_gross = w_gross * age_efficiency * epsilon_val;

                % b. [核心] 分步计算税收和净收入

                % (1) 计算并扣除 PAYG 税 (θ)
                %     PAYG税的税基是全部的税前劳动收入。
                payg_tax = cS.theta_t * labor_income_gross;

                % (2) 计算一般劳动税 (τ_l) 的应税收入
                %     这里体现了 PPS 的税收递延：c_pps_chosen 从应税收入中扣除。
                pps_deduction = 0;
                if isfield(paramS_hh, 'pps_tax_deferral_active') && paramS_hh.pps_tax_deferral_active
                    pps_deduction = c_pps_chosen;
                end
                labor_income_taxable_for_tau_l = labor_income_gross - pps_deduction;

                % (3) 计算并扣除一般劳动税 (τ_l)
                general_labor_tax = cS.tau_l * max(0, labor_income_taxable_for_tau_l);

                % (4) 计算税后劳动净收入
                %     从总收入中减去所有税收
                labor_income_net = labor_income_gross - payg_tax - general_labor_tax;

                % c. 计算当期总的非资本收入
                actual_pps_contribution_expenditure = c_pps_chosen;
                non_capital_income = labor_income_net + TR_total + b_payg_val;

            else % 退休年龄组 (逻辑不变)
                actual_pps_contribution_expenditure = 0;
                pps_withdrawal_net = 0; % 假设PPS提取在退休期收入计算中处理
                non_capital_income = TR_total + b_payg_val + pps_withdrawal_net;
            end

            % --- 2. 计算资本收入 (逻辑不变) ---
            capital_income_net_of_tax = k_now_val * (R_k_net_factor - 1);

            % --- 3. 计算可用于消费和新储蓄的总资源 (逻辑不变) ---
            resources_for_c_and_k_prime = k_now_val + capital_income_net_of_tax + non_capital_income - actual_pps_contribution_expenditure;
        end


        function [cPolM_q, kPolM, cPpsPolM_choice, valM] = HHSolution_VFI_Huggett(...
                R_k_net_factor_vfi, w_gross_vfi, TR_total_vfi, bV_payg_vfi, paramS_vfi, cS, solverMethod)
            % [带求解器切换功能的版本]
            % solverMethod: 'grid' 或 'hybrid'

            % --- 默认使用网格搜索法 ---
            if nargin < 7
                solverMethod = 'hybrid';
            end

            valM = -Inf(cS.nk, cS.nkpps, cS.nw, cS.aD_new);
            cPolM_q = zeros(cS.nk, cS.nkpps, cS.nw, cS.aD_new);
            kPolM = zeros(cS.nk, cS.nkpps, cS.nw, cS.aD_new);
            cPpsPolM_choice = zeros(cS.nk, cS.nkpps, cS.nw, cS.aD_new);
            %
            % fprintf('  🔄 VFI V9-compatible (%s Solver): 开始逆向迭代... [nk=%d, nkpps=%d, nkprime=%d, npps=%d]\n', ...
            %     solverMethod, cS.nk, cS.nkpps, cS.nkprime, cS.npps);
            % vfi_start_time = tic;
            % [进度条] 1. 初始化进度条
            progress_msg = ''; % 用于存储上一条消息的长度
            for a_idx = cS.aD_new : -1 : 1
                % [进度条] 2. 更新并打印进度
                % 首先，删除上一条消息
                % fprintf(repmat('\b', 1, length(progress_msg)));
                % 计算进度
                % percent_done = (cS.aD_new - a_idx + 1) / cS.aD_new * 100;
                % 创建新消息
                % progress_msg = sprintf('    正在处理年龄组 %2d/%d (%.0f%%)...', cS.aD_new - a_idx + 1, cS.aD_new, percent_done);
                % 打印新消息
                % fprintf('%s', progress_msg);


                vPrime_kkppse_next = [];
                if a_idx < cS.aD_new
                    vPrime_kkppse_next = valM(:,:,:,a_idx+1);
                end

                % --- [核心] 根据 solverMethod 选择求解器 ---
                if strcmpi(solverMethod, 'hybrid')
                    [cPolM_q(:,:,:,a_idx), kPolM(:,:,:,a_idx), cPpsPolM_choice(:,:,:,a_idx), valM(:,:,:,a_idx)] = ...
                        main_olg_v13_utils.HHSolutionByAge_VFI_Huggett_HybridOptimizer(a_idx, vPrime_kkppse_next, ...
                        R_k_net_factor_vfi, w_gross_vfi, TR_total_vfi, bV_payg_vfi(a_idx), paramS_vfi, cS);
                elseif strcmpi(solverMethod, 'vectorized_grid') % <<<<< 新增CASE
                    [cPolM_q(:,:,:,a_idx), kPolM(:,:,:,a_idx), cPpsPolM_choice(:,:,:,a_idx), valM(:,:,:,a_idx)] = ...
                        main_olg_v13_utils.HHSolutionByAge_VFI_Huggett_VectorizedGrid(a_idx, vPrime_kkppse_next, ...
                        R_k_net_factor_vfi, w_gross_vfi, TR_total_vfi, bV_payg_vfi(a_idx), paramS_vfi, cS);
                else % 默认或指定 'grid'
                    [cPolM_q(:,:,:,a_idx), kPolM(:,:,:,a_idx), cPpsPolM_choice(:,:,:,a_idx), valM(:,:,:,a_idx)] = ...
                        main_olg_v13_utils.HHSolutionByAge_VFI_Huggett_v9_GridSearch(a_idx, vPrime_kkppse_next, ...
                        R_k_net_factor_vfi, w_gross_vfi, TR_total_vfi, bV_payg_vfi(a_idx), paramS_vfi, cS);
                end

            end
            % [进度条] 3. 循环结束后换行
            % fprintf('\n');
            % total_time = toc(vfi_start_time);
            % fprintf('  ✅ VFI (%s Solver): 完成! 总耗时: %.2f秒\n', solverMethod, total_time);
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
                [pension_income, ~, ~] = main_olg_v13_utils.HHIncome_Huggett(0, 0, 0, TR_total_age, b_age_val, 0, a_idx, paramS_age, cS, 0);

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
                    [~, util_c] = main_olg_v13_utils.CES_utility(final_c, cS.sigma, cS);
                    [~, util_b] = main_olg_v13_utils.CES_utility(final_bequest, cS.sigma_bequest, cS);
                    final_v = util_c + cS.phi_bequest * util_b;

                else
                    % --- 情况B: 没有遗赠动机 (phi_bequest = 0) ---
                    % 将所有资源用于消费，没有遗赠
                    final_bequest = zeros(size(total_resources));
                    final_c_expenditure = total_resources;
                    final_c = max(cS.cFloor, final_c_expenditure / (1 + cS.tau_c));

                    % 最终价值只来源于消费
                    [~, final_v] = main_olg_v13_utils.CES_utility(final_c, cS.sigma, cS);
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
                            [resources, ~, ~] = main_olg_v13_utils.HHIncome_Huggett(k_state, R_k_net_factor_age, w_gross_age, TR_total_age, b_age_val, c_pps_choice, a_idx, paramS_age, cS, epsilon_state);

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

                                [~, util] = main_olg_v13_utils.CES_utility(c_choice, cS.sigma, cS);

                                pps_withdrawal = 0;
                                if a_idx >= cS.aR_new, pps_withdrawal = k_pps_state * cS.pps_withdrawal_rate; end
                                pps_return_factor = 1 + ((R_k_net_factor_age - 1) + cS.pps_return_rate_premium);
                                k_pps_prime = (k_pps_state + c_pps_choice - pps_withdrawal) * pps_return_factor;

                                % [核心修复] 钳位 k_prime_choice 和 k_pps_prime
                                k_prime_clamped = max(cS.kGridV(1), min(cS.kGridV(end), k_prime_choice));
                                k_pps_prime_clamped = max(cS.kppsGridV(1), min(cS.kppsGridV(end), k_pps_prime));

                                % 调用插值器
                                ev = EV_interpolants{ie}(k_prime_clamped, k_pps_prime_clamped);

                                current_val = util + cS.beta * cS.s_5yr_transitionV(a_idx) * ev;
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
                [pension_income, ~, ~] = main_olg_v13_utils.HHIncome_Huggett(0, 0, 0, TR_total_age, b_age_val, 0, a_idx, paramS_age, cS, 0);

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
                [~, util_c] = main_olg_v13_utils.CES_utility(final_c, cS.sigma, cS);
                [~, util_b] = main_olg_v13_utils.CES_utility(final_bequest, cS.sigma_bequest, cS);

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
                            [resources, ~, ~] = main_olg_v13_utils.HHIncome_Huggett(k_state, R_k_net_factor_age, w_gross_age, TR_total_age, b_age_val, c_pps_choice, a_idx, paramS_age, cS, epsilon_state);

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
                            objective_func = @(k_prime) main_olg_v13_utils.objective_for_k_prime_private(...
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


        function [cPol_age_q, kPol_age, cPpsPol_age_choice, val_age] = HHSolutionByAge_VFI_Huggett_VectorizedGrid(...
                a_idx, vPrime_kkppse_next, R_k_net_factor_5year, w_gross_age, TR_total_age, b_age_val, ...
                paramS_age, cS)

            % 初始化输出
            val_age = -Inf(cS.nk, cS.nkpps, cS.nw);
            cPol_age_q = zeros(cS.nk, cS.nkpps, cS.nw);
            kPol_age = zeros(cS.nk, cS.nkpps, cS.nw);
            cPpsPol_age_choice = zeros(cS.nk, cS.nkpps, cS.nw);

            % 准备通用的回报和累积因子
            r_net_annual_vfi = R_k_net_factor_5year^(1/cS.time_Step) - 1;
            r_pps_annual_vfi = r_net_annual_vfi + cS.pps_return_rate_premium;
            pps_return_factor_5year_vfi = (1 + r_pps_annual_vfi)^cS.time_Step;
            if abs(r_net_annual_vfi) > 1e-7
                FlowAcc_5y = ((1+r_net_annual_vfi)^cS.time_Step - 1) / r_net_annual_vfi;
            else
                FlowAcc_5y = cS.time_Step;
            end

            % --- 最终期逻辑 (保持不变) ---
            if a_idx == cS.aD_new
                [K_grid, Kpps_grid] = ndgrid(cS.kGridV, cS.kppsGridV);
                final_value_k = K_grid .* R_k_net_factor_5year;
                final_value_kpps_gross = Kpps_grid .* pps_return_factor_5year_vfi;
                final_value_kpps_net = final_value_kpps_gross * (1 - cS.pps_tax_rate_withdrawal);
                non_capital_income_annual = TR_total_age + b_age_val;
                non_capital_income_accumulated = non_capital_income_annual * FlowAcc_5y;
                total_resources = final_value_k + final_value_kpps_net + non_capital_income_accumulated;
                if cS.phi_bequest > 1e-9
                    phi_adj = cS.phi_bequest * (1 + cS.tau_c)^(-1);
                    omega = (1 + phi_adj^(-1/cS.sigma))^(-1);
                    final_c_expenditure = omega .* total_resources;
                    final_bequest = (1 - omega) .* total_resources;
                else
                    final_c_expenditure = total_resources;
                    final_bequest = zeros(size(total_resources));
                end
                c_annual = max(cS.cFloor, (final_c_expenditure / FlowAcc_5y) / (1 + cS.tau_c));
                [~, util_c] = main_olg_v13_utils.CES_utility(c_annual, cS.sigma, cS);
                [~, util_b] = main_olg_v13_utils.CES_utility(final_bequest, cS.sigma_bequest, cS);
                final_v = util_c * FlowAcc_5y + util_b;
                for ie = 1:cS.nw
                    cPol_age_q(:,:,ie) = c_annual;
                    val_age(:,:,ie) = final_v;
                    kPol_age(:,:,ie) = final_bequest;
                    cPpsPol_age_choice(:,:,ie) = 0;
                end
                return;
            end

            % --- 非最终期逻辑 ---
            EV_matrix = zeros(cS.nk, cS.nkpps, cS.nw);
            for ie_current = 1:cS.nw
                transition_probs = paramS_age.leTrProbM(ie_current, :);
                EV_slice = pagemtimes(reshape(vPrime_kkppse_next, [cS.nk * cS.nkpps, cS.nw]), reshape(transition_probs, [cS.nw, 1]));
                EV_matrix(:, :, ie_current) = reshape(EV_slice, [cS.nk, cS.nkpps]);
            end
            EV_interpolants = cell(cS.nw, 1);
            for ie_current = 1:cS.nw
                EV_interpolants{ie_current} = griddedInterpolant({cS.kGridV, cS.kppsGridV}, EV_matrix(:, :, ie_current), 'linear', 'linear');
            end

            if cS.pps_active
                prop_cpps_grid = linspace(0, cS.pps_max_contrib_frac, cS.npps)';
            else
                prop_cpps_grid = [0];
            end
            prop_k_prime_grid = linspace(0, 1, cS.nkprime)';
            
            for ie = 1:cS.nw
                epsilon_state = paramS_age.leGridV(ie);
                ev_interpolant = EV_interpolants{ie};
                
                npps_current = length(prop_cpps_grid);
                k_state_4D = repmat(reshape(cS.kGridV, [cS.nk, 1, 1, 1]), [1, cS.nkpps, cS.nkprime, npps_current]);
                kpps_state_4D = repmat(reshape(cS.kppsGridV, [1, cS.nkpps, 1, 1]), [cS.nk, 1, cS.nkprime, npps_current]);
                prop_k_prime_4D = repmat(reshape(prop_k_prime_grid, [1, 1, cS.nkprime, 1]), [cS.nk, cS.nkpps, 1, npps_current]);
                prop_cpps_4D = repmat(reshape(prop_cpps_grid, [1, 1, 1, npps_current]), [cS.nk, cS.nkpps, cS.nkprime, 1]);

                % --- [最终核心修正] 简化资源和决策变量 ---
                
                % 1. 计算【年化】非资本税后收入 和 【年化】PPS缴费
                non_capital_income_annual = 0;
                actual_cpps_annual_4D = zeros(size(prop_cpps_4D));
                if a_idx <= cS.aR_new
                    age_efficiency = cS.ageEffV_new(a_idx);
                    labor_income_gross = w_gross_age * age_efficiency * epsilon_state;
                    actual_cpps_annual_4D = min(cS.pps_contrib_limit, labor_income_gross .* prop_cpps_4D);
                    payg_tax = cS.theta_t * labor_income_gross;
                    labor_income_taxable_for_tau_l = labor_income_gross - actual_cpps_annual_4D;
                    general_labor_tax = cS.tau_l * max(0, labor_income_taxable_for_tau_l);
                    labor_income_net = labor_income_gross - payg_tax - general_labor_tax;
                    non_capital_income_annual = labor_income_net + TR_total_age + b_age_val;
                else
                    non_capital_income_annual = TR_total_age + b_age_val;
                end
                
                % 2. 计算【5年期末】的总资源
                wealth_at_end = k_state_4D .* R_k_net_factor_5year;
                non_capital_income_accumulated = non_capital_income_annual .* FlowAcc_5y;
                resources_at_end_4D = wealth_at_end + non_capital_income_accumulated;
                
                if a_idx >= cS.aR_new && cS.pps_active
                    pps_withdrawal_gross_annual = kpps_state_4D * cS.pps_withdrawal_rate;
                    pps_withdrawal_net_annual = pps_withdrawal_gross_annual * (1 - cS.pps_tax_rate_withdrawal);
                    resources_at_end_4D = resources_at_end_4D + pps_withdrawal_net_annual * FlowAcc_5y;
                end
                
                % 3. 决策：直接分配期末总资源
                % 最低消费要求也要转换为期末价值
                c_annual_floor_exp = cS.cFloor * (1 + cS.tau_c);
                c_total_floor_exp = c_annual_floor_exp * FlowAcc_5y;
                
                resources_above_floor_4D = max(0, resources_at_end_4D - c_total_floor_exp);
                actual_k_prime_4D = resources_above_floor_4D .* prop_k_prime_4D;
                c_expend_total_4D = resources_at_end_4D - actual_k_prime_4D;
                
                % 4. 计算效用
                % 为了计算效用，我们需要一个等价的年化消费
                c_annual_equiv_4D = max(cS.cFloor, (c_expend_total_4D ./ FlowAcc_5y) ./ (1 + cS.tau_c));
                [~, util_4D] = main_olg_v13_utils.CES_utility(c_annual_equiv_4D, cS.sigma, cS);
                
                % 5. 演化 PPS 资产
                pps_withdrawal_for_evolution_annual = 0;
                if a_idx >= cS.aR_new && cS.pps_active, pps_withdrawal_for_evolution_annual = kpps_state_4D .* cS.pps_withdrawal_rate; end
                
                k_pps_prime_4D = (kpps_state_4D .* pps_return_factor_5year_vfi) + (actual_cpps_annual_4D .* FlowAcc_5y) - (pps_withdrawal_for_evolution_annual .* FlowAcc_5y);

                % 6. 计算总价值
                k_prime_clamped = max(cS.kGridV(1), min(cS.kGridV(end), actual_k_prime_4D));
                k_pps_prime_clamped = max(cS.kppsGridV(1), min(cS.kppsGridV(end), k_pps_prime_4D));
                ev_mat = ev_interpolant(k_prime_clamped, k_pps_prime_clamped);
                val_grid = (util_4D .* FlowAcc_5y) + cS.beta * cS.s_5yr_transitionV(a_idx) * ev_mat;
                val_grid(c_expend_total_4D < 0) = -Inf;
                
                % 7. 寻找最优并提取结果 (代码不变)
                [val_max_k, idx_k_prime] = max(val_grid, [], 3);
                [val_max_kc, idx_cpps] = max(val_max_k, [], 4);
                val_age(:,:,ie) = squeeze(val_max_kc);

                [I, J] = ndgrid(1:cS.nk, 1:cS.nkpps);
                linear_idx_cpps = sub2ind(size(squeeze(idx_k_prime)), I(:), J(:), squeeze(idx_cpps(:)));
                best_k_prime_idx = idx_k_prime(linear_idx_cpps);
                kPol_age(:,:,ie) = reshape(actual_k_prime_4D(best_k_prime_idx), [cS.nk, cS.nkpps]);
                cPol_age_q(:,:,ie) = reshape(c_annual_equiv_4D(best_k_prime_idx), [cS.nk, cS.nkpps]);
                best_cpps_idx = squeeze(idx_cpps);
                cPpsPol_age_choice(:,:,ie) = reshape(actual_cpps_annual_4D(sub2ind(size(actual_cpps_annual_4D), I(:), J(:), ones(numel(I),1), best_cpps_idx(:))), [cS.nk, cS.nkpps]);
            end
        end        
        
        function [kHistM_out, kPpsHistM_out, cHistM_out, cppsHistM_out] = HHSimulation_olgm(...
                kPolM, cPpsPolM_choice, cPolM_consump, eIdxM_group, ...
                R_k_net_5year, w, TR, bV_payg, paramS_sim, cS_sim) % 重命名输入以明确是5年期因子

            nSim = size(eIdxM_group, 1);
            aD = cS_sim.aD_new;

            % 初始化历史记录矩阵
            kHistM_out = zeros(nSim, aD);
            kPpsHistM_out = zeros(nSim, aD);
            cHistM_out = zeros(nSim, aD);
            cppsHistM_out = zeros(nSim, aD);

            % --- 1. 创建策略函数的插值器 (代码不变) ---
            kPolInterp = cell(cS_sim.nw, aD);
            cPpsPolInterp = cell(cS_sim.nw, aD);
            cPolInterp = cell(cS_sim.nw, aD);

            for ia = 1:aD
                for ie = 1:cS_sim.nw
                    if cS_sim.nk > 1 && cS_sim.nkpps > 1
                        kPolInterp{ie,ia} = griddedInterpolant({cS_sim.kGridV, cS_sim.kppsGridV}, squeeze(kPolM(:,:,ie,ia)), 'linear', 'linear');
                        cPpsPolInterp{ie,ia} = griddedInterpolant({cS_sim.kGridV, cS_sim.kppsGridV}, squeeze(cPpsPolM_choice(:,:,ie,ia)), 'linear', 'linear');
                        cPolInterp{ie,ia} = griddedInterpolant({cS_sim.kGridV, cS_sim.kppsGridV}, squeeze(cPolM_consump(:,:,ie,ia)), 'linear', 'linear');
                    elseif cS_sim.nk > 1
                        kPolInterp{ie,ia} = griddedInterpolant(cS_sim.kGridV, squeeze(kPolM(:,1,ie,ia)), 'linear', 'linear');
                        cPpsPolInterp{ie,ia} = griddedInterpolant(cS_sim.kGridV, squeeze(cPpsPolM_choice(:,1,ie,ia)), 'linear', 'linear');
                        cPolInterp{ie,ia} = griddedInterpolant(cS_sim.kGridV, squeeze(cPolM_consump(:,1,ie,ia)), 'linear', 'linear');
                    else % 假定至少有一个网格点
                        kPolInterp{ie,ia} = @(x,y) squeeze(kPolM(1,1,ie,ia));
                        cPpsPolInterp{ie,ia} = @(x,y) squeeze(cPpsPolM_choice(1,1,ie,ia));
                        cPolInterp{ie,ia} = @(x,y) squeeze(cPolM_consump(1,1,ie,ia));
                    end
                end
            end

            % --- 2. 初始化状态和参数 ---
            % [核心修正] 正确计算PPS的5年期回报因子
            % 假设 cS_sim.pps_return_rate_premium 是一个年化的超额回报率
            r_net_annual_sim = R_k_net_5year^(1/cS_sim.time_Step) - 1;
            r_pps_annual_sim = r_net_annual_sim + cS_sim.pps_return_rate_premium;
            pps_return_factor_5year = (1 + r_pps_annual_sim)^cS_sim.time_Step;

            k_next = zeros(nSim, 1);
            k_pps_next = zeros(nSim, 1);

            % --- 3. 按年龄组进行前向模拟 ---
            for a_idx = 1:aD
                k_now = k_next;
                k_pps_now = k_pps_next;
                kHistM_out(:, a_idx) = k_now;
                kPpsHistM_out(:, a_idx) = k_pps_now;

                k_prime_decision = zeros(nSim, 1);
                cpps_decision = zeros(nSim, 1);
                c_decision = zeros(nSim, 1);

                % --- 4. 根据效率冲击状态，查询策略并获取决策 (代码不变) ---
                for ie = 1:cS_sim.nw
                    idx_sim = find(eIdxM_group(:, a_idx) == ie);
                    if isempty(idx_sim), continue; end

                    k_now_e = k_now(idx_sim);
                    k_pps_now_e = k_pps_now(idx_sim);

                    if cS_sim.nk > 1 && cS_sim.nkpps > 1
                        k_prime_decision(idx_sim) = kPolInterp{ie, a_idx}(k_now_e, k_pps_now_e);
                        cpps_decision(idx_sim) = cPpsPolInterp{ie, a_idx}(k_now_e, k_pps_now_e);
                        c_decision(idx_sim) = cPolInterp{ie, a_idx}(k_now_e, k_pps_now_e);
                    elseif cS_sim.nk > 1
                        k_prime_decision(idx_sim) = kPolInterp{ie, a_idx}(k_now_e);
                        cpps_decision(idx_sim) = cPpsPolInterp{ie, a_idx}(k_now_e);
                        c_decision(idx_sim) = cPolInterp{ie, a_idx}(k_now_e);
                    else
                        k_prime_decision(idx_sim) = kPolInterp{ie, a_idx}(k_now_e, k_pps_now_e);
                        cpps_decision(idx_sim) = cPpsPolInterp{ie, a_idx}(k_now_e, k_pps_now_e);
                        c_decision(idx_sim) = cPolInterp{ie, a_idx}(k_now_e, k_pps_now_e);
                    end
                end

                % --- 5. 记录当期决策并演化到下一期状态 ---
                cHistM_out(:, a_idx) = max(cS_sim.cFloor, c_decision);
                cppsHistM_out(:, a_idx) = max(0, cpps_decision);

                if a_idx < aD
                    k_next = max(cS_sim.kMin, min(cS_sim.kMax, k_prime_decision));

                    pps_withdrawal = 0;
                    if a_idx >= cS_sim.aR_new && cS_sim.pps_active
                        pps_withdrawal = k_pps_now * cS_sim.pps_withdrawal_rate;
                    end

                    % [核心修正] 使用正确计算的5年期PPS回报因子
                    k_pps_next_unclamped = (k_pps_now + cppsHistM_out(:, a_idx) - pps_withdrawal) * pps_return_factor_5year;
                    k_pps_next = max(cS_sim.kppsMin, min(cS_sim.kppsMax, k_pps_next_unclamped));
                end
            end
        end

        function [K_pvt_next, K_pps_next, C_t, Total_Cpps_t, Total_PpsWithdrawalTax_t] = simulate_private_capital_forward(M_t, Z_t_norm, cS, paramS, eIdxM)
            % [v13核心 - VFI决策版, 已修正K_pps聚合]
            % 在给定的宏观环境 M_t下，通过VFI求解家庭的最优策略，
            % 然后模拟所有家庭的行为，最终聚合得到下一期的私人总资本 K_pvt_{t+1}。

            % 1. 准备VFI求解所需的参数
            %    家庭基于静态预期，认为当前的宏观环境 M_t 会永远持续。

            % --- [核心修正] ---
            % VFI求解器模拟的是一个5年期的决策，因此必须使用5年期的回报因子！
            % M_t.r_net_t 存储的是【年化】净利率。
            r_net_annual = M_t.r_net_t;
            R_k_net_factor_5year = (1 + r_net_annual)^cS.time_Step;
            % --------------------

            w_gross = M_t.w_t; % 工资是年化流量, 这在VFI的收入计算中是正确的
            TR_total = 0; % 假设无意外遗赠

            bV_payg = zeros(1, cS.aD_new);
            bV_payg(cS.aR_new+1:end) = M_t.b_t; % 养老金是年化流量, 正确

            paramS_vfi = paramS;
            paramS_vfi.tau_l = cS.tau_l; % 使用固定的tau_l

            % --- 根据标志调整 cS 副本 ---
            cS_vfi = cS; % 创建一个副本，避免修改原始cS
            if isfield(paramS_vfi, 'is_initial_steady_state') && paramS_vfi.is_initial_steady_state
                cS_vfi.pps_active = false; % 暂时禁用PPS系统
            end
            % 从路径中获取当期的有效缴费率
            theta_t = cS.theta_path(M_t.current_t);

            % 将当期的缴费率临时写入副本结构体中
            cS_vfi.theta_t = theta_t;

            % 2. [核心] 调用VFI求解器，传入【5年期】的回报因子
            [cPolM, kPolM, cPpsPolM, ~] = main_olg_v13_utils.HHSolution_VFI_Huggett(...
                R_k_net_factor_5year, w_gross, TR_total, bV_payg, ...
                paramS_vfi, cS_vfi, 'vectorized_grid');

            % 3. [核心] 使用刚算出的策略，模拟所有家庭在 t 期的决策
            %    HHSimulation_olgm 也需要使用【5年期】的回报因子
            [kHistM, kPpsHistM, cHistM, cppsHistM] = main_olg_v13_utils.HHSimulation_olgm(...
                kPolM, cPpsPolM, cPolM, eIdxM, ...
                R_k_net_factor_5year, w_gross, TR_total, bV_payg, paramS_vfi, cS_vfi);

            % 4. [核心] 聚合得到下一期的私人资本 K_pvt_{t+1}
            % 提取所有储蓄决策 (k' 的路径)
            k_prime_paths = kHistM(:, 2:end); % 尺寸为 [nSim, aD_new-1]

            % 计算每个年龄组的平均储蓄决策
            mean_k_prime_by_age = mean(k_prime_paths, 1); % 尺寸为 [1, aD_new-1]

            % 使用 t 期的人口分布 Z_t_norm 来加权
            weights = Z_t_norm(1:cS.aD_new-1); % 尺寸为 [aD_new-1, 1]

            K_pvt_next = mean_k_prime_by_age * weights;

            % [新增] 聚合得到下一期的个人养老金总资产 K_pps_{t+1}
            k_pps_prime_paths = kPpsHistM(:, 2:end);
            mean_k_pps_prime_by_age = mean(k_pps_prime_paths, 1);
            K_pps_next = mean_k_pps_prime_by_age * weights;

            % 检查: 在初始稳态下，PPS资产和缴费应为零
            if isfield(paramS_vfi, 'is_initial_steady_state') && paramS_vfi.is_initial_steady_state
                if any(kPpsHistM(:) > 1e-6) || any(cppsHistM(:) > 1e-6)
                    warning('在初始稳态求解中，PPS资产或缴费不为零，请检查逻辑。');
                end
            end

            % 5. [新增] 聚合得到当期总消费 C_t
            mean_c_by_age = mean(cHistM, 1); % 每个年龄组的平均消费
            weights_c = Z_t_norm; % 当期所有年龄组的人口权重
            C_t = mean_c_by_age * weights_c;

            % 在 simulate_private_capital_forward 的聚合部分
            mean_cpps_by_age = mean(cppsHistM, 1);
            Total_Cpps_t = mean_cpps_by_age * Z_t_norm;

            % 6. 聚合当期因PPS提取而产生的总税收收入
            Total_PpsWithdrawalTax_t = 0;
            if cS.pps_active && cS.pps_withdrawal_rate > 0
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
            [leLogGridV_raw, leTrProbM] = main_olg_v13_utils.tauchen(cS.nw, lePersistence, leShockStd, 0, Tauchen_q);
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
            eIdxM_group = main_olg_v13_utils.MarkovChainSimulation_AgeGroup(cS.nSim, cS, paramS.leProb1V, paramS.leTrProbM);
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
        % =========================================================================
        % == 在 main_olg_v13_utils.m 中，在文件末尾的 end 之前，添加这个新函数 ==
        % =========================================================================

        function debug_national_accounts(t, paths, M_t, cS)
            % 这是一个强大的诊断工具，用于检查t期国民账户的平衡情况。

            fprintf('\n\n=============================================================\n');
            fprintf('--- National Accounts Debug for t = %d (Year: %d) ---\n', t, cS.start_year + (t-1)*cS.time_Step);
            fprintf('=============================================================\n');

            % --- 1. 商品市场 (Y = C + I + G) ---
            Y = paths.Y(t);
            C = paths.C(t);
            I = paths.Investment(t);
            G = paths.G(t);
            Residual = paths.market_clearing_residual(t);

            fprintf('\n--- 1. Goods Market (Y = C + I + G) ---\n');
            fprintf('Supply (Y)                     : %12.4f\n', Y);
            fprintf('-----------------------------------------------------\n');
            fprintf('Demand Components:\n');
            fprintf('  Consumption (C)              : %12.4f\n', C);
            fprintf('  Investment (I)               : %12.4f\n', I);
            fprintf('  Government Spending (G)      : %12.4f\n', G);
            fprintf('Total Demand (C+I+G)           : %12.4f\n', C+I+G);
            fprintf('-----------------------------------------------------\n');
            fprintf('Residual (Y - [C+I+G])         : %12.4f (%.4f%% of Y)\n', Residual, (Residual/Y)*100);

            % --- 2. 收入法GDP (Y = wL + rK) ---
            wL = M_t.w_t * M_t.L_t * cS.time_Step; % 5年期总劳动收入
            rK_gross = M_t.r_mkt_t * M_t.K_physical_t * cS.time_Step; % 5年期总资本收入

            fprintf('\n--- 2. Income-Side GDP (Y = wL + r_gross*K) ---\n');
            fprintf('Gross Labor Income (w*L)       : %12.4f\n', wL);
            fprintf('Gross Capital Income (r*K)     : %12.4f\n', rK_gross);
            fprintf('Total Factor Income (wL+rK)    : %12.4f\n', wL + rK_gross);
            fprintf('Production-Side Y              : %12.4f\n', Y);
            fprintf('-----------------------------------------------------\n');
            fprintf('Residual (Y - [wL+rK])         : %12.4f\n', Y - (wL + rK_gross));

            % --- 3. 储蓄-投资恒等式 (S_private + S_public = I) ---
            % 这需要计算总税收和总转移支付

            % 计算总税收
            Y_annual = Y / cS.time_Step;
            C_annual = C / cS.time_Step;
            % 为了获取详细的税收，我们可能需要重新计算一下
            [~, Total_Cpps_t, Total_PpsTax_t] = main_olg_v13_utils.simulate_private_capital_forward(M_t, M_t.Z_t_norm, cS, M_t.paramS, M_t.eIdxM); % 重新模拟以获取流量

            Total_Cpps_annual = Total_Cpps_t / cS.time_Step;
            Total_PpsTax_annual = Total_PpsTax_t / cS.time_Step;

            GrossLaborIncome_annual = M_t.w_t * M_t.L_t;
            TaxableLaborIncome_annual = GrossLaborIncome_annual - Total_Cpps_annual;
            LaborTaxRev_annual = cS.tau_l * max(0, TaxableLaborIncome_annual);

            TaxableCapitalStock = M_t.K_physical_t - M_t.K_pps_t;
            CapitalTaxRev_annual = cS.tau_k * M_t.r_mkt_t * TaxableCapitalStock;

            ConsumptionTaxRev_annual = cS.tau_c * C_annual;

            TotalTax_annual = LaborTaxRev_annual + CapitalTaxRev_annual + ConsumptionTaxRev_annual + Total_PpsTax_annual;

            % 计算总转移支付 (养老金)
            mass_retirees = sum(M_t.Z_t_norm(cS.aR_new+1:end));
            PensionOutlay_annual = M_t.b_t * mass_retirees;

            % 计算储蓄
            S_private_annual = (Y_annual - TotalTax_annual + PensionOutlay_annual) - C_annual;
            S_public_annual = TotalTax_annual - (G/cS.time_Step) - PensionOutlay_annual;
            S_national_annual = S_private_annual + S_public_annual;

            fprintf('\n--- 3. National Savings (S = Y - C - G) ---\n');
            fprintf('National Savings (Y-C-G)/step : %12.4f\n', (Y - C - G) / cS.time_Step);
            fprintf('Calculated S_national_annual   : %12.4f\n', S_national_annual);
            fprintf('Net Investment (I/step)        : %12.4f\n', I / cS.time_Step);
            fprintf('-----------------------------------------------------\n');
            fprintf('Residual (S_nat - I_net)       : %12.4f\n', S_national_annual - (I / cS.time_Step));

            fprintf('\n================== End of Debug Report ==================\n\n');
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
            [~, util] = main_olg_v13_utils.CES_utility(c_choice, cS.sigma, cS);

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
            ev = main_olg_v13_utils.CallInterpolator(ev_interpolant, k_prime_clamped, k_pps_prime_clamped, cS);

            % 5. 计算总价值
            current_val = util + cS.beta * cS.s_5yr_transitionV(a_idx) * ev;

            % fminbnd是最小化器，所以返回负价值
            neg_v = -current_val;
        end
    end % End of Static Methods
end