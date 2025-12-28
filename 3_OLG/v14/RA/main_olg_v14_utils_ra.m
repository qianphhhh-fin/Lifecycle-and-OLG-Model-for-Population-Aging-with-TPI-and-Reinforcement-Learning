% --- START OF FILE main_olg_v14_utils_ra.m ---
% =========================================================================
% == UTILS for Representative Agent (RA) Model Audit (Final Corrected Version v2)
% == 目的：提供一个极度简化的RA模型环境，用于精确检验宏观会计恒等式。
% =========================================================================

classdef main_olg_v14_utils_ra

    methods (Static)

        function cS_ra = ParameterValues_RA(cS_base)
            cS_ra = cS_base;
            cS_ra.time_step = 1;
            cS_ra.L_fixed = 0.4;
            fields_to_remove = {
                'age1_orig', 'ageLast_orig', 'ageRetire_orig', 'aD_orig', 'aR_idx_orig', 'aW_orig', 'physAgeV_orig', 'aD_new', 'aR_new', 'physAgeMap', 's_1yr_transitionV', 'ageEffV_new', ...
                'lePersistence', 'leShockStd', 'nw', 'phi_bequest', 'sigma_bequest', ...
                'rho_prime_payg', 'pps_active', 'pps_activation_year', 'pps_tax_rate_withdrawal', 'pps_return_rate_premium', 'pps_withdrawal_rate', 'pps_contrib_limit', 'pps_max_contrib_frac', 'theta_path'
            };
            for i = 1:length(fields_to_remove)
                if isfield(cS_ra, fields_to_remove{i})
                    cS_ra = rmfield(cS_ra, fields_to_remove{i});
                end
            end
            ddk_annual_equiv = 1 - (1-cS_base.ddk)^(1/cS_base.time_step);
            cS_ra.ddk = 1 - (1 - ddk_annual_equiv)^cS_ra.time_step;
        end

        function [Z_path, A_path, T_sim] = load_exogenous_paths_RA(cS)
            fprintf('--- [RA] 正在独立加载和处理外生TFP路径 ---\n');
            Z_path = [];
            sim_years = cS.start_year:cS.time_Step:cS.end_year;
            T_sim = length(sim_years);
            tfp_data_pwt = readtable('..\data\PWT\china_pwt_data.xlsx');
            bai_projections = [2020, 0.0628; 2025, 0.0557; 2030, 0.0482; 2035, 0.0394; 2040, 0.0340; 2045, 0.0346; 2050, 0.0298; 2102, 0.0150];
            g_A_annual = zeros(cS.end_year - cS.start_year + 1, 1);
            tfp_pwt_series = tfp_data_pwt.rtfpna(tfp_data_pwt.year >= cS.start_year & tfp_data_pwt.year <= 2019);
            for i = 1:(length(tfp_pwt_series)-1), g_A_annual(i) = tfp_pwt_series(i+1) / tfp_pwt_series(i) - 1; end
            current_year_idx = 2020 - cS.start_year;
            g_A_annual(current_year_idx) = (tfp_pwt_series(end) / tfp_pwt_series(end-1)) -1;
            for i = 1:size(bai_projections, 1)-1
                proj_start_year = bai_projections(i, 1); proj_end_year = bai_projections(i+1, 1); proj_rate = bai_projections(i, 2);
                for year = (proj_start_year + 1):proj_end_year
                    if year <= cS.end_year, idx = year - cS.start_year + 1; g_A_annual(idx) = proj_rate; end
                end
            end
            g_A_period_factor = ones(T_sim, 1);
            for t = 1:T_sim
                annual_idx = sim_years(t) - cS.start_year + 1;
                if annual_idx <= length(g_A_annual), g_A_period_factor(t) = 1 + g_A_annual(annual_idx);
                else, g_A_period_factor(t) = 1 + bai_projections(end,2); end
            end
            A_path = cumprod([1; g_A_period_factor(1:end-1)]);
            fprintf('✅ TFP路径已独立加载 (%d个时期)。\n', T_sim);
        end

        function M_t = get_prices_and_policy_at_t_RA(t, K_pvt_t, B_g_t, A_t, cS)
            K_physical_t = K_pvt_t - B_g_t;
            if K_physical_t <= 0, K_physical_t = 1e-8; end
            L_t = cS.L_fixed;
            [r_mkt_period, w_t, Y_period] = main_olg_v14_utils.HHPrices_Huggett(K_physical_t, L_t, A_t, cS);
            M_t = struct();
            M_t.current_t = t;
            M_t.K_physical_t = K_physical_t; M_t.L_t = L_t;
            M_t.Y_t = Y_period; M_t.w_t = w_t;
            M_t.r_mkt_annual = r_mkt_period; % Use the same name for compatibility
            M_t.K_pvt_t = K_pvt_t; M_t.B_g_t = B_g_t;
        end

        % =========================================================================
        % == [核心修正] 确保家庭预算约束通过构造恒等 ==
        % =========================================================================
        function [K_pvt_next, C_t] = simulate_RA_capital_forward(M_t, cS)
            % 1. 计算当期可支配收入流 (Yd_flow)
            r_mkt = M_t.r_mkt_annual;
            NetLaborIncome_flow = M_t.w_t * M_t.L_t * (1 - cS.tau_l);
            NetCapitalIncome_flow = r_mkt * M_t.K_physical_t * (1 - cS.tau_k);
            BondIncome_flow = r_mkt * M_t.B_g_t;
            DisposableIncome_flow = NetLaborIncome_flow + NetCapitalIncome_flow + BondIncome_flow;

            % 2. 计算期末可用于分配的总资源
            TotalResources = M_t.K_pvt_t + DisposableIncome_flow;

            % 3. 确定消费和储蓄 (考虑cFloor)
            % 意愿消费支出
            Desired_ConsumptionExpenditure = (1 - cS.beta) * TotalResources;
            Desired_C = Desired_ConsumptionExpenditure / (1 + cS.tau_c);
            
            % 实际消费量
            C_t = max(cS.cFloor, Desired_C);
            
            % 实际消费支出
            Actual_ConsumptionExpenditure = C_t * (1 + cS.tau_c);
            
            % 实际储蓄 (由预算约束倒推)
            K_pvt_next = TotalResources - Actual_ConsumptionExpenditure;
        end
    end
end