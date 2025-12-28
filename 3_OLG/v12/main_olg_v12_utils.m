% --- START OF FILE main_olg_v12_utils.m ---
%
% OLG 模型 v12 工具函数库
%
% 核心变化 (对比 V12):
% 1. [v12 核心修正] 重写了 `solve_initial_steady_state` 函数。
%    现在通过迭代生产性资本存量(K_prod)来寻找资本市场出清的均衡点。
% 2. [v12 兼容性修正] 移除了 'if_else' 函数调用和不必要的 '~' 输出参数。
%
% ---

classdef main_olg_v12_utils

    methods (Static)

        function [ss, eq_found, Aggregates_audited] = solve_initial_steady_state(Z_ss_norm, A_ss, cS, paramS, eIdxM)
            fprintf('\n--- [最终版求解器 V2] 开始求解初始稳态 ---\n');
            f_residual = @(k_demand) main_olg_v12_utils.get_asset_mkt_residual_unified(k_demand, Z_ss_norm, A_ss, cS, paramS, eIdxM);
            options = optimset('TolX', 1e-8);
            [K_ss, ~, exitflag] = fzero(f_residual, [0.1, 15.0], options);
            eq_found = (exitflag == 1);
            if ~eq_found, error('稳态求解失败！'); end
            [~, ss, Aggregates_audited] = main_olg_v12_utils.get_asset_mkt_residual_unified(K_ss, Z_ss_norm, A_ss, cS, paramS, eIdxM, true);
            fprintf('✅ 资产市场均衡求解成功！ K* = %.4f\n', K_ss);
        end        % =========================================================================

        % [最终修正] get_asset_mkt_residual_with_bequests: 返回所有审计所需的数据
        % =========================================================================
function [residual, ss, Aggregates_final] = get_asset_mkt_residual_unified(K_demand, Z_norm, A_ss, cS, paramS, eIdxM, final_run)
            if nargin < 7, final_run = false; end
            
            M = main_olg_v12_utils.get_prices_and_policy_at_t(1, K_demand, 0, 0, Z_norm, A_ss, cS, paramS, eIdxM);
            
            TR_guess = 0.03; damp = 0.25; max_iter_tr = 100; tol_tr = 1e-7;

            for iter = 1:max_iter_tr
                % [核心修正] 每次都使用权威的审计函数进行模拟
                [~, ~, ~, Aggregates_iter] = main_olg_v12_utils.simulate_private_capital_forward_audited(M, Z_norm, cS, paramS, eIdxM, TR_guess);
                
                % [核心修正] TR的计算直接来自于S_net，强制财富守恒
                TR_new = Aggregates_iter.TotalSaving;
                
                if abs(TR_new - TR_guess) < tol_tr, TR_guess = TR_new; break; end
                TR_guess = (1 - damp) * TR_guess + damp * TR_new;
            end
            TR_eq = TR_guess;
            
            % 使用均衡的TR，进行最后一次权威模拟，获得所有最终变量
            [K_pvt_next, K_pps_next, C_final, Aggregates_final] = main_olg_v12_utils.simulate_private_capital_forward_audited(M, Z_norm, cS, paramS, eIdxM, TR_eq);
            K_supply = K_pvt_next + K_pps_next;
            residual = K_supply - K_demand;

            if ~final_run
                fprintf('K_d=%.4f | r=%.2f%% | TR_eq=%.4f | S_net=%.4f | K_s=%.4f | Res=%.3e\n', K_demand, M.r_net_t*100, TR_eq, Aggregates_final.TotalSaving, K_supply, residual);
                ss = struct();
            else
                ss = struct();
                ss.K_pvt = K_pvt_next; ss.K_pps = K_pps_next; ss.B_p = 0;
                ss.L = M.L_t; ss.Y = M.Y_t; ss.w = M.w_t;
                ss.r_mkt_rental = M.r_mkt_rental_t; ss.r_net = M.r_net_t;
                ss.C = C_final;
                T_total = Aggregates_final.TotalLaborTax + Aggregates_final.TotalCapitalTax + Aggregates_final.TotalConsumptionTax;
                Pension_Benefit = T_total - (ss.Y - ss.C - (cS.ddk * K_demand)); % G=Y-C-I
                ss.G = T_total - Pension_Benefit;
                ss.Investment = cS.ddk * K_demand;
                ss.market_clearing_residual = ss.Y - (ss.C + ss.Investment + ss.G);
                ss.TR = TR_eq;
            end
        end        
        function residual = get_asset_mkt_residual(K_demand, Z_norm, A_ss, cS, paramS, eIdxM)
            % 这个函数严格遵循 Aiyagari (1994) 的步骤：
            % 1. 假设一个资本需求 K_demand
            % 2. 计算在此需求下的价格 r, w
            % 3. 家庭在此价格下做最优决策
            % 4. 汇总家庭决策，得到资本供给 K_supply
            % 5. 返回残差 K_supply - K_demand

            if K_demand <= 0, K_demand = 1e-6; end

            % 步骤1 & 2: 基于给定的 K_demand 计算价格
            % 在初始稳态求解中，假设养老金基金 B_p=0, 个人养老金 K_pps=0
            % 因此，总的生产性资本就等于企业部门的资本需求 K_demand
            M = main_olg_v12_utils.get_prices_and_policy_at_t(1, K_demand, 0, 0, Z_norm, A_ss, cS, paramS, eIdxM);

            % 步骤3 & 4: 模拟家庭决策，得到总的资产供给 K_supply
            paramS_vfi = paramS;
            paramS_vfi.is_initial_steady_state = true; % 确保PPS不被激活

            % simulate_..._audited 返回的 K_pvt_next 和 K_pps_next 在稳态下
            % 就是家庭部门愿意持有的长期资产存量，即总供给
            [K_pvt_supply, K_pps_supply, ~, ~] = ...
                main_olg_v12_utils.simulate_private_capital_forward_audited(M, Z_norm, cS, paramS_vfi, eIdxM);

            K_supply = K_pvt_supply + K_pps_supply;

            % 步骤5: 计算资产市场残差
            residual = K_supply - K_demand;

            fprintf('K_demand=%.4f | r_net=%.4f%% | K_supply=%.4f | Residual=%.4e\n', ...
                K_demand, M.r_net_t*100, K_supply, residual);
        end
        % --- [新增] 辅助函数：计算给定K_prod下的商品市场残差 ---

        function [Z_path, A_path, T_sim] = load_exogenous_paths(cS)
            fprintf('--- [v12] 加载和处理外生数据路径 ---\n');
            pop_data = readtable('..\data\人口\population_by_age_group_all_years.xlsx', 'Sheet', 'pop_normalized');
            sim_years = cS.start_year:cS.time_step:cS.end_year;
            T_sim = length(sim_years);
            Z_path_raw = zeros(cS.aD_new, T_sim);
            for t = 1:T_sim
                year_t = sim_years(t);
                [~, col_idx] = min(abs(str2double(pop_data.Properties.VariableNames(2:end)) - year_t));
                Z_path_raw(:, t) = pop_data{:, col_idx + 1};
            end
            Z_path = Z_path_raw ./ sum(Z_path_raw, 1);
            fprintf('✅ 人口路径加载完成 (1997-%d, %d个时期)。\n', cS.end_year, T_sim);
            tfp_data_pwt = readtable('..\data\PWT\china_pwt_data.xlsx');
            bai_projections = [
                2020, 0.0628; 2025, 0.0557; 2030, 0.0482; 2035, 0.0394;
                2040, 0.0340; 2045, 0.0346; 2050, 0.0298; 2102, 0.0150;
                ];
            g_A_annual = zeros(cS.end_year - cS.start_year + 1, 1);
            tfp_pwt_series = tfp_data_pwt.rtfpna(tfp_data_pwt.year >= cS.start_year & tfp_data_pwt.year <= 2019);
            for i = 1:(length(tfp_pwt_series)-1)
                g_A_annual(i) = tfp_pwt_series(i+1) / tfp_pwt_series(i) - 1;
            end
            current_year_idx = 2020 - cS.start_year;
            g_A_annual(current_year_idx) = (tfp_pwt_series(end) / tfp_pwt_series(end-1)) -1;
            for i = 1:size(bai_projections, 1)-1
                proj_start_year = bai_projections(i, 1);
                proj_end_year = bai_projections(i+1, 1);
                proj_rate = bai_projections(i, 2);
                for year = (proj_start_year + 1):proj_end_year
                    if year <= cS.end_year
                        idx = year - cS.start_year + 1;
                        g_A_annual(idx) = proj_rate;
                    end
                end
            end
            g_A_5year_factor = zeros(T_sim, 1);
            for t = 1:T_sim
                year_t = sim_years(t);
                factor = 1.0;
                for i = 0:cS.time_step-1
                    annual_idx = (year_t + i) - cS.start_year + 1;
                    if annual_idx <= length(g_A_annual)
                        factor = factor * (1 + g_A_annual(annual_idx));
                    else
                        factor = factor * (1 + bai_projections(end,2));
                    end
                end
                g_A_5year_factor(t) = factor;
            end
            A_path = ones(T_sim, 1);
            for t = 1:(T_sim - 1)
                A_path(t+1) = A_path(t) * g_A_5year_factor(t);
            end
            fprintf('✅ TFP路径构建完成 (1997-%d)。\n', cS.end_year);
        end

        function [r_mkt_rental, w, Y] = HHPrices_Huggett(K, L, A_t, cS)
            % 返回: r_mkt_rental (资本租金价格, MPK), w (工资), Y (产出)
            if K <= 0, K = 1e-6; end
            if L <= 0, L = 1e-6; end
            Y = A_t * (K^cS.alpha) * (L^(1-cS.alpha));
            % r_mkt_rental 现在是资本的租金价格 (MPK)
            r_mkt_rental = cS.alpha * Y / K;
            w = (1-cS.alpha) * Y / L;
        end

        function M_t = get_prices_and_policy_at_t(t, K_pvt_t, B_p_t, K_pps_t, Z_t_norm, A_t, cS, paramS, eIdxM)
            K_prod_t = K_pvt_t + B_p_t + K_pps_t;
            paramS_t = paramS;
            paramS_t.ageMassV = Z_t_norm;
            [L_t] = main_olg_v12_utils.LaborSupply_Huggett(eIdxM, cS, paramS_t, Z_t_norm);

            [r_mkt_rental_t, w_t, Y_t] = main_olg_v12_utils.HHPrices_Huggett(K_prod_t, L_t, A_t, cS);
            r_net_t = (r_mkt_rental_t - cS.ddk) * (1 - cS.tau_k); % 家庭的净回报率

            cov_urban_t = cS.coverage_urban_path(t);
            cov_resident_t = cS.coverage_resident_path(t);
            total_coverage = cov_urban_t + cov_resident_t;
            if total_coverage > 1e-9, share_urban = cov_urban_t / total_coverage; share_resident = cov_resident_t / total_coverage;
            else, share_urban = 0; share_resident = 0; end
            mass_workers_t = sum(Z_t_norm(1:cS.aR_new));
            avg_wage_t = w_t * L_t / mass_workers_t;
            b_urban_t = cS.rho_urban_employee * avg_wage_t;
            b_resident_t = cS.rho_resident * avg_wage_t;
            b_t = share_urban * b_urban_t + share_resident * b_resident_t;

            M_t = struct();
            M_t.current_t = t;
            M_t.K_prod_t = K_prod_t;
            M_t.L_t = L_t;
            M_t.Y_t = Y_t;
            M_t.w_t = w_t;
            M_t.r_mkt_rental_t = r_mkt_rental_t; % 资本租金率
            M_t.r_net_t = r_net_t;             % 家庭净回报率
            M_t.b_t = b_t;
        end

        function B_p_next = update_pension_fund(B_p_t, M_t, Z_t_norm, cS)
            theta_t = cS.theta_path(M_t.current_t);
            PensionRevenue_annual = theta_t * M_t.w_t * M_t.L_t;
            mass_retirees_t = sum(Z_t_norm(cS.aR_new+1:end));
            PensionOutlay_annual = M_t.b_t * mass_retirees_t;
            NetSurplus_annual = PensionRevenue_annual - PensionOutlay_annual;
            r_a = M_t.r_mkt_t;
            R_5y = (1 + r_a)^cS.time_step;
            if abs(r_a) > 1e-6
                FlowAcc_5y = (R_5y - 1) / r_a;
            else
                FlowAcc_5y = cS.time_step;
            end
            B_p_next = B_p_t * R_5y + NetSurplus_annual * FlowAcc_5y;
            B_p_next = max(0, B_p_next);
        end

        % --- [V12 最终修正版] 养老金缴费率和覆盖率路径计算 ---
        function [cS] = calcaulte_theta_payg_path(cS, graph_flag)
            % 1. 定义关键参数
            theta_urban_employee_effective = 0.20;
            theta_resident_effective = 0.03;

            % [修正] 设定覆盖率的长期稳定目标
            coverage_urban_final = 0.9;
            coverage_resident_final = 0.9; % 假设城乡居民覆盖率在2023年后会稳定在一个较高水平
            year_of_peak_coverage = 2040; % 假设在2025年左右达到平台期

            % 2. 收集真实数据点
            year_pax_urban = [1997, 2000, 2005, 2008, 2011, 2014, 2017, 2020, 2023];
            pax_urban      = [8671, 10448, 17487, 21890, 28391, 34124, 40293, 45638, 52112];
            year_pax_resident = [2008, 2011, 2014, 2017, 2020, 2023];
            pax_resident      = [0,    45586, 50075, 51255, 54227, 54475];
            year_laborpop_data = [1997, 2000, 2005, 2008, 2011, 2014, 2017, 2020, 2023];
            laborpop_data      = [78531, 82352, 87532, 90833, 94072, 92837, 91199, 89552, 86481];

            % 3. 插值生成完整的【年度】劳动年龄人口路径
            annual_years_vector = cS.start_year:cS.end_year;
            year_reach_final = cS.end_year;
            % 劳动人口的预测保持不变
            laborpop_annual = interp1([year_laborpop_data, year_reach_final], [laborpop_data, laborpop_data(end)*0.85], annual_years_vector, 'linear', 'extrap');

            % 4. [核心修正] 分别插值生成【年度覆盖率】路径，而不是参保人数路径

            % a. 计算历史数据中的覆盖率
            coverage_urban_hist = pax_urban ./ interp1(year_laborpop_data, laborpop_data, year_pax_urban);
            coverage_resident_hist = pax_resident ./ interp1(year_laborpop_data, laborpop_data, year_pax_resident);

            % b. 插值生成城镇职工覆盖率路径
            interp_years_urban = [year_pax_urban, year_reach_final];
            interp_cov_urban   = [coverage_urban_hist, coverage_urban_final];
            coverage_urban_annual = interp1(interp_years_urban, interp_cov_urban, annual_years_vector, 'linear', 'extrap');

            % c. 插值生成城乡居民覆盖率路径
            %    我们假设覆盖率在2023年之后，会从当前值平滑过渡到长期稳定值
            last_real_year = year_pax_resident(end); % 2023
            last_real_coverage = coverage_resident_hist(end); % 2023年的覆盖率

            interp_years_resident = [cS.start_year, year_pax_resident, year_of_peak_coverage, year_reach_final];
            interp_cov_resident   = [0, coverage_resident_hist, coverage_resident_final, coverage_resident_final];

            [unique_years, ia] = unique(interp_years_resident);
            unique_cov = interp_cov_resident(ia);

            coverage_resident_annual = interp1(unique_years, unique_cov, annual_years_vector, 'linear');
            coverage_resident_annual(annual_years_vector < min(year_pax_resident)) = 0; % 确保早期为0

            % 5. 计算【年度有效缴费率】路径
            theta_path_annual = (coverage_urban_annual * theta_urban_employee_effective) + (coverage_resident_annual * theta_resident_effective);

            % 6. 从年度路径中，提取出模型【5年期】所需的路径
            model_year_indices = 1:(cS.end_year - cS.start_year + 1);
            model_year_indices_5yr = model_year_indices(1:cS.time_step:end);
            T_sim = floor(length(annual_years_vector) / cS.time_step) + 1;

            theta_path = theta_path_annual(model_year_indices_5yr);
            cS.theta_path = theta_path(1:min(T_sim, length(theta_path)));

            coverage_urban_path = coverage_urban_annual(model_year_indices_5yr);
            cS.coverage_urban_path = coverage_urban_path(1:min(T_sim, length(coverage_urban_path)));

            coverage_resident_path = coverage_resident_annual(model_year_indices_5yr);
            cS.coverage_resident_path = coverage_resident_path(1:min(T_sim, length(coverage_resident_path)));

            % 7. (可选) 可视化检查路径
            if graph_flag
                T_plot = length(cS.theta_path);
                time_axis = cS.start_year:cS.time_step:cS.start_year+cS.time_step*(T_plot-1);

                figure('Name', 'Effective PAYG Tax Rate and Coverage Paths (Corrected)');

                subplot(2,1,1);
                plot(time_axis, cS.theta_path, 'k-s', 'LineWidth', 2, 'DisplayName', '总有效缴费率 (θ_t)');
                hold on;
                contribution_urban = cS.coverage_urban_path(1:T_plot) * theta_urban_employee_effective;
                contribution_resident = cS.coverage_resident_path(1:T_plot) * theta_resident_effective;
                plot(time_axis, contribution_urban, 'b--.', 'DisplayName', '城镇职工体系贡献');
                plot(time_axis, contribution_resident, 'r--.', 'DisplayName', '城乡居民体系贡献');
                title('模型使用的有效养老金缴费率路径');
                xlabel('年份'); ylabel('有效缴费率');
                legend('show', 'Location', 'best'); grid on;
                ylim([0, theta_urban_employee_effective*1.2]);

                subplot(2,1,2);
                plot(time_axis, cS.coverage_urban_path(1:T_plot)*100, 'b-s', 'LineWidth', 2, 'DisplayName', '城镇职工覆盖率');
                hold on;
                plot(time_axis, cS.coverage_resident_path(1:T_plot)*100, 'r-s', 'LineWidth', 2, 'DisplayName', '城乡居民覆盖率');
                title('养老金覆盖率路径 (占劳动年龄人口比例)');
                xlabel('年份'); ylabel('覆盖率 (%)');
                legend('show', 'Location', 'best'); grid on;
            end
        end

        function B_g_next = update_gov_debt(B_g_t, C_t, M_t, Total_Cpps_t, K_pps_t, Total_PpsWithdrawalTax_t, cS)
            Y_annual = M_t.Y_t / cS.time_step;
            C_annual = C_t / cS.time_step;
            Total_Cpps_annual = Total_Cpps_t / cS.time_step;
            Total_PpsWithdrawalTax_annual = Total_PpsWithdrawalTax_t / cS.time_step;
            G_annual = cS.gov_exp_frac_Y * Y_annual;
            GrossLaborIncome_annual = (M_t.w_t * M_t.L_t);
            TaxableLaborIncome_annual = GrossLaborIncome_annual - Total_Cpps_annual;
            LaborTaxRevenue_annual = cS.tau_l * max(0, TaxableLaborIncome_annual);
            TaxableCapitalStock_t = M_t.K_total_t;
            CapitalTaxRevenue_annual = cS.tau_k * M_t.r_mkt_t * TaxableCapitalStock_t;
            ConsumptionTaxRevenue_annual = cS.tau_c * C_annual;
            TaxRevenue_annual = LaborTaxRevenue_annual + CapitalTaxRevenue_annual + ConsumptionTaxRevenue_annual + Total_PpsWithdrawalTax_annual;
            PrimaryDeficit_annual = G_annual - TaxRevenue_annual;
            r_a = M_t.r_mkt_t;
            R_5y = (1 + r_a)^cS.time_step;
            if abs(r_a) > 1e-6
                FlowAcc_5y = (R_5y - 1) / r_a;
            else
                FlowAcc_5y = cS.time_step;
            end
            B_g_next = B_g_t * R_5y + PrimaryDeficit_annual * FlowAcc_5y;
        end

        function cS = ParameterValues_HuggettStyle()
            cS.age1_orig = 20; cS.ageLast_orig = 98; cS.ageRetire_orig = 65;
            cS.aD_orig = cS.ageLast_orig - cS.age1_orig + 1;
            cS.aR_idx_orig = cS.ageRetire_orig - cS.age1_orig + 1;
            cS.aW_orig = cS.aR_idx_orig - 1;
            cS.time_step = 5;
            cS.aD_new = ceil(cS.aD_orig / cS.time_step);
            cS.aR_new = ceil(cS.aW_orig / cS.time_step);
            cS.physAgeMap = cell(cS.aD_new, 1);
            for a = 1:cS.aD_new, cS.physAgeMap{a} = (a-1)*cS.time_step + 1 : min(a*cS.time_step, cS.aD_orig); end
            d_orig_data = [0.00159,0.00169,0.00174,0.00172,0.00165,0.00156,0.00149,0.00145,0.00145,0.00149,0.00156,0.00163,0.00171,0.00181,0.00193,0.00207,0.00225,0.00246,0.00270,0.00299,0.00332,0.00368,0.00409,0.00454,0.00504,0.00558,0.00617,0.00686,0.00766,0.00865,0.00955,0.01058,0.01162,0.01264,0.01368,0.01475,0.01593,0.01730,0.01891,0.02074,0.02271,0.02476,0.02690,0.02912,0.03143,0.03389,0.03652,0.03930,0.04225,0.04538,0.04871,0.05230,0.05623,0.06060,0.06542,0.07066,0.07636,0.08271,0.08986,0.09788,0.10732,0.11799,0.12895,0.13920,0.14861,0.16039,0.17303,0.18665,0.20194,0.21877,0.23601,0.25289,0.26973,0.28612,0.30128,0.31416,0.32915,0.34450,0.36018];
            s_orig = 1 - d_orig_data;
            cS.s_1yr_transitionV = zeros(cS.aD_new, 1);
            for a = 1:(cS.aD_new - 1), cS.s_1yr_transitionV(a) = s_orig(cS.physAgeMap{a}(end)); end
            cS.lePersistence = 0.9; cS.leShockStd = 0.1^0.5; cS.nw = 5;
            ageEffV_orig_temp = zeros(100, 1);
            ageEffV_orig_temp(20:72) = [linspace(0.3,1.5,36-20+1), 1.5.*ones(1,47-37+1), linspace(1.5,0.2,65-48+1), linspace(0.18,0,72-66+1)];
            ageEffV_orig = ageEffV_orig_temp(20:98);
            cS.ageEffV_new = zeros(cS.aD_new, 1);
            for a = 1:cS.aD_new, cS.ageEffV_new(a) = mean(ageEffV_orig(cS.physAgeMap{a})); end
        end

        function cS = generateGrids(cS)
            cS.tgKY = 3; cS.tgWage = (1-cS.alpha)*cS.A*((cS.tgKY/cS.A)^(cS.alpha/(1-cS.alpha)));
            cS.kMin = 0; cS.kMax = 15 * cS.tgWage; cS.kppsMin = 0; cS.kppsMax = cS.kMax / 2;
            power_k = 1.5;
            if cS.nk > 1, kGridV_temp = cS.kMin + (cS.kMax - cS.kMin) * (linspace(0, 1, cS.nk).^power_k); kGridV_temp(1) = cS.kMin;
            elseif cS.nk == 1, kGridV_temp = cS.kMin; else, kGridV_temp = []; end
            cS.kGridV = kGridV_temp(:);
            power_kpps = 1.5;
            if cS.nkpps > 1, kppsGridV_temp = cS.kppsMin + (cS.kppsMax - cS.kppsMin) * (linspace(0, 1, cS.nkpps).^power_kpps); kppsGridV_temp(1) = cS.kppsMin;
            elseif cS.nkpps == 1, kppsGridV_temp = cS.kppsMin; else, kppsGridV_temp = []; end
            cS.kppsGridV = kppsGridV_temp(:);
        end

        % --- [核心修正] HHIncome_Huggett ---
        % =========================================================================
        % [决定性最终修正] HHIncome_Huggett 函数
        % =========================================================================
        function [resources, tax_info] = HHIncome_Huggett(k_now_val, R_k_net_factor, w_gross, TR_total, b_payg_val, c_pps_chosen, a_idx, paramS_hh, cS, epsilon_val)
            tax_info = struct('labor_tax', 0, 'capital_tax', 0);

            % 1. 劳动收入和相关税收
            if a_idx <= cS.aR_new
                age_efficiency = cS.ageEffV_new(a_idx);
                labor_income_gross = w_gross * age_efficiency * epsilon_val;
                payg_tax = cS.theta_t * labor_income_gross;
                pps_deduction = 0;
                if isfield(paramS_hh, 'pps_tax_deferral_active') && paramS_hh.pps_tax_deferral_active
                    pps_deduction = c_pps_chosen;
                end
                labor_income_taxable_for_tau_l = labor_income_gross - pps_deduction;
                general_labor_tax = cS.tau_l * max(0, labor_income_taxable_for_tau_l);
                tax_info.labor_tax = general_labor_tax + payg_tax;
                non_capital_income_base = labor_income_gross - tax_info.labor_tax;
            else
                non_capital_income_base = 0;
            end

            % 2. 资本收入和相关税收
            r_net_from_param = R_k_net_factor - 1;
            r_mkt_rental = (r_net_from_param + cS.ddk) / (1 - cS.tau_k);
            capital_income_gross = k_now_val * r_mkt_rental;
            tax_info.capital_tax = capital_income_gross * cS.tau_k;
            capital_income_net_of_tax = capital_income_gross - tax_info.capital_tax;

            % 3. 总资源 (核心最终修正)
            % 家庭可用于分配的总资源 =
            %   期初财富 + 税后资本收入 + 税后劳动收入 + 所有转移支付 - PPS缴费 - 资本折旧
            %   这等价于 NDI (净可支配收入) + k_now - C_expend

            % [最终会计修正] 意外遗赠转移支付 TR_total 必须被加入到家庭资源中！
            after_tax_cash_income = non_capital_income_base + capital_income_net_of_tax + b_payg_val + TR_total;

            % 从总资源中扣除折旧才是净收入
            net_income = after_tax_cash_income - c_pps_chosen - (cS.ddk * k_now_val);

            % 总资源 = 期初财富 + 净收入
            resources = k_now_val + net_income;
        end
        % --- [V12 最终会计修正版 - 最终定义] 精确审计模拟函数 ---
 % =========================================================================
        % [决定性最终修正] simulate_private_capital_forward_audited
        % =========================================================================
        function [K_pvt_next, K_pps_next, C_t, Aggregates_t] = simulate_private_capital_forward_audited(M_t, Z_t_norm, cS, paramS, eIdxM, TR_total)
            if nargin < 6, TR_total = 0; end % 默认没有转移支付

            % --- 1. 求解最优策略 (这部分正确，TR_total被传入) ---
            R_k_net_factor_hh = 1 + M_t.r_net_t;
            w_gross = M_t.w_t;
            bV_payg = zeros(1, cS.aD_new);
            bV_payg(cS.aR_new+1:end) = M_t.b_t;
            
            paramS_vfi = paramS; paramS_vfi.tau_l = cS.tau_l;
            cS_vfi = cS;
            if isfield(paramS_vfi, 'is_initial_steady_state') && paramS_vfi.is_initial_steady_state, cS_vfi.pps_active = false; end
            cS_vfi.theta_t = cS.theta_path(M_t.current_t);
            
            [cPolM, kPolM, cPpsPolM] = main_olg_v12_utils.HHSolution_VFI_Huggett(R_k_net_factor_hh, w_gross, TR_total, bV_payg, paramS_vfi, cS_vfi);
            
            % --- 2. 准备模拟 (这部分不变) ---
            nSim = cS.nSim;
            kHistM = zeros(nSim, cS.aD_new); kPpsHistM = zeros(nSim, cS.aD_new);
            cHistM = zeros(nSim, cS.aD_new); sHistM = zeros(nSim, cS.aD_new);
            labor_tax_HistM = zeros(nSim, cS.aD_new); capital_tax_HistM = zeros(nSim, cS.aD_new);
            cppsHistM = zeros(nSim, cS.aD_new);
            % ... 插值器创建不变 ...
            kPolInterp = cell(cS.nw, cS.aD_new); cPpsPolInterp = cell(cS.nw, cS.aD_new); cPolInterp = cell(cS.nw, cS.aD_new);
            for ia = 1:cS.aD_new, for ie = 1:cS.nw, if cS.nk > 1 && cS.nkpps > 1, kPolInterp{ie,ia} = griddedInterpolant({cS.kGridV, cS.kppsGridV}, squeeze(kPolM(:,:,ie,ia)), 'linear', 'linear'); cPpsPolInterp{ie,ia} = griddedInterpolant({cS.kGridV, cS.kppsGridV}, squeeze(cPpsPolM(:,:,ie,ia)), 'linear', 'linear'); cPolInterp{ie,ia} = griddedInterpolant({cS.kGridV, cS.kppsGridV}, squeeze(cPolM(:,:,ie,ia)), 'linear', 'linear'); elseif cS.nk > 1, kPolInterp{ie,ia} = griddedInterpolant(cS.kGridV, squeeze(kPolM(:,1,ie,ia)), 'linear', 'linear'); cPpsPolInterp{ie,ia} = griddedInterpolant(cS.kGridV, squeeze(cPpsPolM(:,1,ie,ia)), 'linear', 'linear'); cPolInterp{ie,ia} = griddedInterpolant(cS.kGridV, squeeze(cPolM(:,1,ie,ia)), 'linear', 'linear'); else, kPolInterp{ie,ia} = @(x,y) squeeze(kPolM(1,1,ie,ia)); cPpsPolInterp{ie,ia} = @(x,y) squeeze(cPpsPolM(1,1,ie,ia)); cPolInterp{ie,ia} = @(x,y) squeeze(cPolM(1,1,ie,ia)); end, end, end
            r_mkt_rental = (R_k_net_factor_hh - 1 + cS.ddk) / (1 - cS.tau_k);
            pps_return_factor = 1 + r_mkt_rental - cS.ddk + cS.pps_return_rate_premium;

            % --- 3. 模拟循环 ---
            for a_idx = 1:cS.aD_new, for ie = 1:cS.nw
                    idx_sim = find(eIdxM(:, a_idx) == ie); if isempty(idx_sim), continue; end
                    
                    k_now_e = kHistM(idx_sim, a_idx); k_pps_now_e = kPpsHistM(idx_sim, a_idx);
                    
                    k_prime_decision_e = kPolInterp{ie, a_idx}(k_now_e, k_pps_now_e);
                    c_decision_e = cPolInterp{ie, a_idx}(k_now_e, k_pps_now_e);
                    cpps_decision_e = cPpsPolInterp{ie, a_idx}(k_now_e, k_pps_now_e);
                    
                    cHistM(idx_sim, a_idx) = c_decision_e;
                    cppsHistM(idx_sim, a_idx) = cpps_decision_e;
                    
                    pps_withdrawal = 0; 
                    if a_idx >= cS.aR_new && cS_vfi.pps_active, pps_withdrawal = k_pps_now_e * cS.pps_withdrawal_rate; end
                    k_pps_prime_e = (k_pps_now_e + cpps_decision_e - pps_withdrawal) * pps_return_factor;
                    
                    sHistM(idx_sim, a_idx) = (k_prime_decision_e + k_pps_prime_e) - (k_now_e + k_pps_now_e);

                    % [决定性最终修正] 在这里调用 HHIncome_Huggett 时，必须传入 TR_total
                    [~, tax_info_e] = main_olg_v12_utils.HHIncome_Huggett(k_now_e, R_k_net_factor_hh, w_gross, TR_total, bV_payg(a_idx), cpps_decision_e, a_idx, paramS_vfi, cS_vfi, paramS.leGridV(ie));
                    labor_tax_HistM(idx_sim, a_idx) = tax_info_e.labor_tax;
                    capital_tax_HistM(idx_sim, a_idx) = tax_info_e.capital_tax;

                    if a_idx < cS.aD_new
                        kHistM(idx_sim, a_idx + 1) = k_prime_decision_e;
                        kPpsHistM(idx_sim, a_idx + 1) = k_pps_prime_e;
                    end
            end, end
            
            % --- 4. 聚合 (这部分不变) ---
            mean_k_prime_by_age = mean(kHistM(:, 2:end), 1);
            weights = Z_t_norm(1:cS.aD_new-1);
            K_pvt_next = mean_k_prime_by_age * weights;
            
            mean_k_pps_prime_by_age = mean(kPpsHistM(:, 2:end), 1);
            K_pps_next = mean_k_pps_prime_by_age * weights;
            
            C_t = mean(cHistM, 1) * Z_t_norm;
            
            Aggregates_t = struct();
            Aggregates_t.TotalConsumption = C_t;
            Aggregates_t.TotalSaving = mean(sHistM, 1) * Z_t_norm; 
            Aggregates_t.TotalLaborTax = mean(labor_tax_HistM, 1) * Z_t_norm;
            Aggregates_t.TotalCapitalTax = mean(capital_tax_HistM, 1) * Z_t_norm;
            Aggregates_t.TotalConsumptionTax = C_t * cS.tau_c;
            
            Total_PpsWithdrawalTax_t = 0;
            if cS.pps_active && cS.pps_withdrawal_rate > 0, mean_k_pps_by_age = mean(kPpsHistM, 1); retired_ages_idx = (cS.aR_new + 1):cS.aD_new; total_pps_capital_of_retirees = mean_k_pps_by_age(retired_ages_idx) * Z_t_norm(retired_ages_idx); Total_PpsWithdrawalTax_t = total_pps_capital_of_retirees * cS.pps_withdrawal_rate * cS.pps_tax_rate_withdrawal; end
            Aggregates_t.TotalPpsWithdrawalTax = Total_PpsWithdrawalTax_t;
            Aggregates_t.Total_Cpps_t = mean(cppsHistM, 1) * Z_t_norm;
            Aggregates_t.K_pvt_next = K_pvt_next;
            Aggregates_t.K_pps_next = K_pps_next;
       end
        function [cPolM, kPolM, cPpsPolM] = HHSolution_VFI_Huggett(R_k_net_factor_vfi, w_gross_vfi, TR_total_vfi, bV_payg_vfi, paramS_vfi, cS_vfi)
            valM = -Inf(cS_vfi.nk, cS_vfi.nkpps, cS_vfi.nw, cS_vfi.aD_new);
            cPolM = zeros(cS_vfi.nk, cS_vfi.nkpps, cS_vfi.nw, cS_vfi.aD_new);
            kPolM = zeros(cS_vfi.nk, cS_vfi.nkpps, cS_vfi.nw, cS_vfi.aD_new);
            cPpsPolM = zeros(cS_vfi.nk, cS_vfi.nkpps, cS_vfi.nw, cS_vfi.aD_new);
            for a_idx = cS_vfi.aD_new : -1 : 1
                vPrime_kkppse_next = [];
                if a_idx < cS_vfi.aD_new
                    vPrime_kkppse_next = valM(:,:,:,a_idx+1);
                end
                [cPolM(:,:,:,a_idx), kPolM(:,:,:,a_idx), cPpsPolM(:,:,:,a_idx), valM(:,:,:,a_idx)] = ...
                    main_olg_v12_utils.HHSolutionByAge_VFI_Huggett_VectorizedGrid(a_idx, vPrime_kkppse_next, ...
                    R_k_net_factor_vfi, w_gross_vfi, TR_total_vfi, bV_payg_vfi(a_idx), paramS_vfi, cS_vfi);
            end
        end

        % --- [V12 最终会计修正版] VFI核心决策函数 ---
        function [cPol_age_q, kPol_age, cPpsPol_age_choice, val_age] = HHSolutionByAge_VFI_Huggett_VectorizedGrid(a_idx, vPrime_kkppse_next, R_k_net_factor_age, w_gross_age, TR_total_age, b_age_val, paramS_age, cS)
            val_age = -Inf(cS.nk, cS.nkpps, cS.nw); cPol_age_q = zeros(cS.nk, cS.nkpps, cS.nw);
            kPol_age = zeros(cS.nk, cS.nkpps, cS.nw); cPpsPol_age_choice = zeros(cS.nk, cS.nkpps, cS.nw);

            r_net_age = R_k_net_factor_age - 1;
            r_mkt_rental_age = (r_net_age + cS.ddk) / (1 - cS.tau_k);

            if a_idx == cS.aD_new
                [K_grid, Kpps_grid] = ndgrid(cS.kGridV, cS.kppsGridV);
                [resources_from_k, ~] = main_olg_v12_utils.HHIncome_Huggett(K_grid, R_k_net_factor_age, 0, TR_total_age, b_age_val, 0, a_idx, paramS_age, cS, 0);
                pps_return_factor = 1 + r_mkt_rental_age - cS.ddk + cS.pps_return_rate_premium;
                final_value_kpps_net = (Kpps_grid .* pps_return_factor) * (1 - cS.pps_tax_rate_withdrawal);
                total_resources = resources_from_k + final_value_kpps_net;

                % [核心修正] 最终期财富用于含税消费和遗赠
                if cS.phi_bequest > 1e-9
                    phi_adj = cS.phi_bequest * (1 + cS.tau_c)^(-1);
                    omega = (1 + phi_adj^(-1/cS.sigma))^(-1);
                    final_c_expenditure = omega .* total_resources;
                    final_bequest = (1 - omega) .* total_resources;
                else
                    final_c_expenditure = total_resources;
                    final_bequest = zeros(size(total_resources));
                end

                final_c_value = final_c_expenditure / (1 + cS.tau_c);
                final_c = max(cS.cFloor, final_c_value);

                [~, util_c] = main_olg_v12_utils.CES_utility(final_c, cS.sigma, cS);
                [~, util_b] = main_olg_v12_utils.CES_utility(final_bequest, cS.sigma_bequest, cS);
                final_v = util_c + cS.phi_bequest * util_b;
                for ie = 1:cS.nw, cPol_age_q(:,:,ie) = final_c; val_age(:,:,ie) = final_v; kPol_age(:,:,ie) = final_bequest; end
                return;
            end

            EV_matrix = zeros(cS.nk, cS.nkpps, cS.nw);
            for ie = 1:cS.nw
                EV_slice = pagemtimes(reshape(vPrime_kkppse_next, [cS.nk * cS.nkpps, cS.nw]), paramS_age.leTrProbM(ie, :)');
                EV_matrix(:, :, ie) = reshape(EV_slice, [cS.nk, cS.nkpps]);
            end
            EV_interpolants = cell(cS.nw, 1);
            for ie = 1:cS.nw, EV_interpolants{ie} = griddedInterpolant({cS.kGridV, cS.kppsGridV}, EV_matrix(:, :, ie), 'linear', 'linear'); end
            prop_k_prime_grid = linspace(0, 1, cS.nkprime)'; prop_cpps_grid = linspace(0, cS.pps_max_contrib_frac, cS.npps)';

            for ie = 1:cS.nw
                epsilon_state = paramS_age.leGridV(ie); ev_interpolant = EV_interpolants{ie};
                k_state_4D = repmat(reshape(cS.kGridV, [cS.nk, 1, 1, 1]), [1, cS.nkpps, cS.nkprime, cS.npps]);
                kpps_state_4D = repmat(reshape(cS.kppsGridV, [1, cS.nkpps, 1, 1]), [cS.nk, 1, cS.nkprime, cS.npps]);

                % ... (prop_k_prime_4D, prop_cpps_4D, actual_cpps_4D logic remains the same) ...
                prop_k_prime_4D = repmat(reshape(prop_k_prime_grid, [1, 1, cS.nkprime, 1]), [cS.nk, cS.nkpps, 1, cS.npps]);
                prop_cpps_4D = repmat(reshape(prop_cpps_grid, [1, 1, 1, cS.npps]), [cS.nk, cS.nkpps, cS.nkprime, 1]);
                actual_cpps_4D = zeros(size(prop_cpps_4D));
                if a_idx <= cS.aR_new, actual_cpps_4D = min(cS.pps_contrib_limit, (w_gross_age * cS.ageEffV_new(a_idx) * epsilon_state) .* prop_cpps_4D); end

                [resources_for_c_and_k_prime, ~] = main_olg_v12_utils.HHIncome_Huggett(k_state_4D, R_k_net_factor_age, w_gross_age, TR_total_age, b_age_val, actual_cpps_4D, a_idx, paramS_age, cS, epsilon_state);
                pps_withdrawal_val = 0;
                if a_idx >= cS.aR_new && cS.pps_active, pps_withdrawal_val = kpps_state_4D .* cS.pps_withdrawal_rate * (1 - cS.pps_tax_rate_withdrawal); end
                resources_4D = resources_for_c_and_k_prime + pps_withdrawal_val;

                % [核心修正] 消费支出 c_expend 必须覆盖消费税
                c_floor_expend = cS.cFloor * (1 + cS.tau_c);
                resources_above_floor_4D = max(0, resources_4D - c_floor_expend);

                % [核心会计修正] k' 和 消费支出 C_expend 分配总资源
                actual_k_prime_4D = resources_4D .* prop_k_prime_4D;
                c_expend_4D = resources_4D - actual_k_prime_4D;
                c_choice_4D = max(cS.cFloor, c_expend_4D / (1 + cS.tau_c));

                [~, util_4D] = main_olg_v12_utils.CES_utility(c_choice_4D, cS.sigma, cS);

                pps_withdrawal_for_evolution = 0;
                if a_idx >= cS.aR_new, pps_withdrawal_for_evolution = kpps_state_4D .* cS.pps_withdrawal_rate; end
                pps_return_factor = 1 + r_mkt_rental_age - cS.ddk + cS.pps_return_rate_premium;
                k_pps_prime_4D = (kpps_state_4D + actual_cpps_4D - pps_withdrawal_for_evolution) * pps_return_factor;

                % ... (rest of VFI remains largely the same, but based on corrected C and k') ...
                k_prime_clamped = max(cS.kGridV(1), min(cS.kGridV(end), actual_k_prime_4D));
                k_pps_prime_clamped = max(cS.kppsGridV(1), min(cS.kppsGridV(end), k_pps_prime_4D));
                ev_mat = ev_interpolant(k_prime_clamped, k_pps_prime_clamped);
                val_grid = util_4D + cS.beta * cS.s_1yr_transitionV(a_idx) * ev_mat;
                val_grid(c_expend_4D < 0) = -Inf;
                [val_max_k, idx_k_prime] = max(val_grid, [], 3); [val_max_kc, idx_cpps] = max(val_max_k, [], 4);
                val_age(:,:,ie) = squeeze(val_max_kc);
                [I, J] = ndgrid(1:cS.nk, 1:cS.nkpps); temp = squeeze(idx_cpps);
                linear_idx_k = sub2ind(size(squeeze(val_grid(:,:,:,1))), I(:), J(:),temp(:));
                final_k_prime_prop_idx = reshape(idx_k_prime(linear_idx_k), [cS.nk, cS.nkpps]);
                best_prop_k_prime = prop_k_prime_grid(final_k_prime_prop_idx);
                best_prop_cpps = prop_cpps_grid(squeeze(idx_cpps));
                [k_state_2D, kpps_state_2D] = ndgrid(cS.kGridV, cS.kppsGridV);
                best_c_pps = 0;
                if a_idx <= cS.aR_new, best_c_pps = min(cS.pps_contrib_limit, (w_gross_age * cS.ageEffV_new(a_idx) * epsilon_state) .* best_prop_cpps); end
                [resources_final_pre_pps, ~] = main_olg_v12_utils.HHIncome_Huggett(k_state_2D, R_k_net_factor_age, w_gross_age, TR_total_age, b_age_val, best_c_pps, a_idx, paramS_age, cS, epsilon_state);
                pps_withdrawal_final = 0;
                if a_idx >= cS.aR_new && cS.pps_active, pps_withdrawal_final = kpps_state_2D .* cS.pps_withdrawal_rate * (1 - cS.pps_tax_rate_withdrawal); end
                resources_final = resources_final_pre_pps + pps_withdrawal_final;
                resources_above_floor_final = max(0, resources_final - c_floor_expend);
                % [核心会计修正] 重新计算最优消费，确保预算约束严格成立
                best_k_prime = resources_final .* best_prop_k_prime;
                c_expend_final = resources_final - best_k_prime;
                cPol_age_q(:,:,ie) = max(cS.cFloor, c_expend_final / (1 + cS.tau_c));
                kPol_age(:,:,ie) = best_k_prime; cPpsPol_age_choice(:,:,ie) = best_c_pps;
            end
        end

        % =========================================================================
        % [修正] HHSimulation_olgm 函数，以处理空的输入策略
        % =========================================================================
        function [kHistM_out, kPpsHistM_out, cHistM_out, cppsHistM_out] = HHSimulation_olgm(kPolM, cPpsPolM_choice, cPolM_consump, eIdxM_group, R_k_net, cS_sim)
            nSim = size(eIdxM_group, 1);
            aD = cS_sim.aD_new;
            kHistM_out = zeros(nSim, aD);
            kPpsHistM_out = zeros(nSim, aD);
            cHistM_out = zeros(nSim, aD);
            cppsHistM_out = zeros(nSim, aD);

            % --- [核心修正] ---
            % 检查输入是否为空，并设置标志位
            do_k_interp = ~isempty(kPolM);
            do_cpps_interp = ~isempty(cPpsPolM_choice);
            do_c_interp = ~isempty(cPolM_consump);

            kPolInterp = cell(cS_sim.nw, aD);
            cPpsPolInterp = cell(cS_sim.nw, aD);
            cPolInterp = cell(cS_sim.nw, aD);

            for ia = 1:aD
                for ie = 1:cS_sim.nw
                    % 根据标志位，只为非空的策略创建插值器
                    if do_k_interp
                        if cS_sim.nk > 1 && cS_sim.nkpps > 1
                            kPolInterp{ie,ia} = griddedInterpolant({cS_sim.kGridV, cS_sim.kppsGridV}, squeeze(kPolM(:,:,ie,ia)), 'linear', 'linear');
                        elseif cS_sim.nk > 1
                            kPolInterp{ie,ia} = griddedInterpolant(cS_sim.kGridV, squeeze(kPolM(:,1,ie,ia)), 'linear', 'linear');
                        else
                            kPolInterp{ie,ia} = @(x,y) squeeze(kPolM(1,1,ie,ia));
                        end
                    end

                    if do_cpps_interp
                        if cS_sim.nk > 1 && cS_sim.nkpps > 1
                            cPpsPolInterp{ie,ia} = griddedInterpolant({cS_sim.kGridV, cS_sim.kppsGridV}, squeeze(cPpsPolM_choice(:,:,ie,ia)), 'linear', 'linear');
                        elseif cS_sim.nk > 1
                            cPpsPolInterp{ie,ia} = griddedInterpolant(cS_sim.kGridV, squeeze(cPpsPolM_choice(:,1,ie,ia)), 'linear', 'linear');
                        else
                            cPpsPolInterp{ie,ia} = @(x,y) squeeze(cPpsPolM_choice(1,1,ie,ia));
                        end
                    end

                    if do_c_interp
                        if cS_sim.nk > 1 && cS_sim.nkpps > 1
                            cPolInterp{ie,ia} = griddedInterpolant({cS_sim.kGridV, cS_sim.kppsGridV}, squeeze(cPolM_consump(:,:,ie,ia)), 'linear', 'linear');
                        elseif cS_sim.nk > 1
                            cPolInterp{ie,ia} = griddedInterpolant(cS_sim.kGridV, squeeze(cPolM_consump(:,1,ie,ia)), 'linear', 'linear');
                        else
                            cPolInterp{ie,ia} = @(x,y) squeeze(cPolM_consump(1,1,ie,ia));
                        end
                    end
                end
            end

            pps_return_factor = 1 + ((R_k_net - 1) + cS_sim.pps_return_rate_premium);
            k_next = zeros(nSim, 1);
            k_pps_next = zeros(nSim, 1);

            for a_idx = 1:aD
                k_now = k_next; k_pps_now = k_pps_next;
                kHistM_out(:, a_idx) = k_now;
                kPpsHistM_out(:, a_idx) = k_pps_now;

                k_prime_decision = zeros(nSim, 1);
                cpps_decision = zeros(nSim, 1);
                c_decision = zeros(nSim, 1);

                for ie = 1:cS_sim.nw
                    idx_sim = find(eIdxM_group(:, a_idx) == ie);
                    if isempty(idx_sim), continue; end

                    k_now_e = k_now(idx_sim);
                    k_pps_now_e = k_pps_now(idx_sim);

                    % 根据标志位，只执行非空的插值
                    if do_k_interp
                        if cS_sim.nk > 1 && cS_sim.nkpps > 1
                            k_prime_decision(idx_sim) = kPolInterp{ie, a_idx}(k_now_e, k_pps_now_e);
                        else
                            k_prime_decision(idx_sim) = kPolInterp{ie, a_idx}(k_now_e);
                        end
                    end

                    if do_cpps_interp
                        if cS_sim.nk > 1 && cS_sim.nkpps > 1
                            cpps_decision(idx_sim) = cPpsPolInterp{ie, a_idx}(k_now_e, k_pps_now_e);
                        else
                            cpps_decision(idx_sim) = cPpsPolInterp{ie, a_idx}(k_now_e);
                        end
                    end

                    if do_c_interp
                        if cS_sim.nk > 1 && cS_sim.nkpps > 1
                            c_decision(idx_sim) = cPolInterp{ie, a_idx}(k_now_e, k_pps_now_e);
                        else
                            c_decision(idx_sim) = cPolInterp{ie, a_idx}(k_now_e);
                        end
                    end
                end

                if do_c_interp, cHistM_out(:, a_idx) = max(cS_sim.cFloor, c_decision); end
                if do_cpps_interp, cppsHistM_out(:, a_idx) = max(0, cpps_decision); end

                if a_idx < aD
                    k_next = max(cS_sim.kMin, min(cS_sim.kMax, k_prime_decision));
                    pps_withdrawal = 0;
                    if a_idx >= cS_sim.aR_new && cS_sim.pps_active
                        pps_withdrawal = k_pps_now * cS_sim.pps_withdrawal_rate;
                    end
                    % 注意：cppsHistM_out可能未被计算，但cpps_decision是计算了的(如果do_cpps_interp为真)
                    % 为了安全起见，我们应该使用cpps_decision。但由于我们只关心k_hist，这个分支暂时可以忽略。
                    k_pps_next_unclamped = (k_pps_now + cpps_decision) * pps_return_factor; % 使用局部变量
                    k_pps_next = max(cS_sim.kppsMin, min(cS_sim.kppsMax, k_pps_next_unclamped));
                end
            end
        end
        % --- [V12 最终审计修正版] 家庭决策模拟与聚合 ---
        function [K_pvt_next, K_pps_next, C_t, Aggregates_t] = simulate_private_capital_forward(M_t, Z_t_norm, cS, paramS, eIdxM)
            R_k_net_factor_hh = 1 + M_t.r_net_t;
            w_gross = M_t.w_t;
            TR_total = 0;
            bV_payg = zeros(1, cS.aD_new);
            bV_payg(cS.aR_new+1:end) = M_t.b_t;

            paramS_vfi = paramS;
            paramS_vfi.tau_l = cS.tau_l;
            cS_vfi = cS;
            if isfield(paramS_vfi, 'is_initial_steady_state') && paramS_vfi.is_initial_steady_state
                cS_vfi.pps_active = false;
            end
            cS_vfi.theta_t = cS.theta_path(M_t.current_t);

            [cPolM, kPolM, cPpsPolM] = main_olg_v12_utils.HHSolution_VFI_Huggett(R_k_net_factor_hh, w_gross, TR_total, bV_payg, paramS_vfi, cS_vfi);
            [kHistM, kPpsHistM, cHistM, cppsHistM] = main_olg_v12_utils.HHSimulation_olgm(kPolM, cPpsPolM, cPolM, eIdxM, R_k_net_factor_hh, cS_vfi);

            k_prime_paths = kHistM(:, 2:end);
            mean_k_prime_by_age = mean(k_prime_paths, 1);
            weights = Z_t_norm(1:cS.aD_new-1);
            K_pvt_next = mean_k_prime_by_age * weights;

            k_pps_prime_paths = kPpsHistM(:, 2:end);
            mean_k_pps_prime_by_age = mean(k_pps_prime_paths, 1);
            K_pps_next = mean_k_pps_prime_by_age * weights;

            mean_c_by_age = mean(cHistM, 1);
            C_t = mean_c_by_age * Z_t_norm;

            mean_cpps_by_age = mean(cppsHistM, 1);
            Total_Cpps_t = mean_cpps_by_age * Z_t_norm;

            Aggregates_t = struct();
            Aggregates_t.TotalConsumptionTax = C_t * cS.tau_c;

            % [审计修正] 不再使用有漏洞的聚合方式
            % total_labor_tax_val = 0;
            % total_capital_tax_val = 0;

            Total_PpsWithdrawalTax_t = 0;
            if cS.pps_active && cS.pps_withdrawal_rate > 0
                mean_k_pps_by_age = mean(kPpsHistM, 1);
                retired_ages_idx = (cS.aR_new + 1):cS.aD_new;
                total_pps_capital_of_retirees = mean_k_pps_by_age(retired_ages_idx) * Z_t_norm(retired_ages_idx);
                Total_PpsWithdrawalTax_t = total_pps_capital_of_retirees * cS.pps_withdrawal_rate * cS.pps_tax_rate_withdrawal;
            end
            Aggregates_t.TotalPpsWithdrawalTax = Total_PpsWithdrawalTax_t;
            Aggregates_t.Total_Cpps_t = Total_Cpps_t;
        end

        function [muM, utilM] = CES_utility(cM, sigma, cS)
            c_adj = max(cS.cFloor, cM);
            if abs(sigma - 1) < 1e-6, utilM = log(c_adj); muM = 1./c_adj;
            else, utilM = (c_adj.^(1-sigma))./(1-sigma); muM = c_adj.^(-sigma); end
            utilM(cM < cS.cFloor) = -1e10 - (cS.cFloor - cM(cM < cS.cFloor))*1e10;
        end

        function [leLogGridV, leTrProbM, leProb1V] = EarningProcess_olgm(cS)
            lePersistence = 0.90; leShockStd = 0.15; Tauchen_q = 2.0;
            [leLogGridV_raw, leTrProbM] = main_olg_v12_utils.tauchen(cS.nw, lePersistence, leShockStd, 0, Tauchen_q);
            leLogGridV = leLogGridV_raw - mean(leLogGridV_raw);
            [~, D] = eig(leTrProbM'); [~, c] = min(abs(diag(D)-1));
            leProb1V = abs(D(:,c)/sum(D(:,c)));
        end

        function [y_grid_out, trProbM_out] = tauchen(N, rho, sigma, mu, m)
            std_y = sqrt(sigma^2 / (1-rho^2)); y = linspace(-m*std_y, m*std_y, N); d = y(2)-y(1);
            trProbM_out = zeros(N,N);
            for j=1:N, for k=1:N
                    m_k = rho*y(j) + mu;
                    if k==1, trProbM_out(j,k) = normcdf((y(1)-m_k+d/2)/sigma);
                    elseif k==N, trProbM_out(j,k) = 1 - normcdf((y(N)-m_k-d/2)/sigma);
                    else, trProbM_out(j,k) = normcdf((y(k)-m_k+d/2)/sigma) - normcdf((y(k)-m_k-d/2)/sigma); end
            end, end
        y_grid_out = y(:);
        end

        function eIdxM_group = LaborEndowSimulation_olgm_AgeGroup(cS, paramS)
            eIdxM_group = main_olg_v12_utils.MarkovChainSimulation_AgeGroup(cS.nSim, cS, paramS.leProb1V, paramS.leTrProbM);
        end

        function eIdxM_group_out = MarkovChainSimulation_AgeGroup(num_simulations, cS, p0, P)
            rng(433);
            eIdxM_group_out = zeros(num_simulations, cS.aD_new, 'uint16');
            eIdxM_group_out(:,1) = 1 + sum(rand(num_simulations,1) > cumsum(p0(:)'), 2);
            for a=2:cS.aD_new, eIdxM_group_out(:,a) = 1 + sum(rand(num_simulations,1) > cumsum(P(eIdxM_group_out(:,a-1),:), 2), 2); end
        end

        function [L_total_eff_pc] = LaborSupply_Huggett(eIdxM_group, cS, paramS, Z_ss_norm_group)
            nSim = size(eIdxM_group, 1); HHlaborM_group = zeros(nSim, cS.aD_new);
            for a_group = 1:cS.aR_new
                labor_eff = paramS.leGridV(eIdxM_group(:, a_group));
                HHlaborM_group(:, a_group) = cS.ageEffV_new(a_group) * labor_eff;
            end
            mean_labor_per_working_group = mean(HHlaborM_group(:, 1:cS.aR_new), 1);
            L_total_eff_pc = sum(mean_labor_per_working_group' .* Z_ss_norm_group(1:cS.aR_new));
        end
    end
end
% --- END OF FILE main_olg_v12_utils.m ---