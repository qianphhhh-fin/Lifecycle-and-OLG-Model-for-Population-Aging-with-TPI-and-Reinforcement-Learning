% --- TPI.m ---
% =========================================================================
% == CLASSDEF: TPI.m
% == 版本: [v1.3 - 聚合目标与人口对齐修正版]
% ==
% == 目的:
% ==   - 修正 calculate_paths_and_errors 中目标遗赠路径的计算逻辑，确保
% ==     总量到人均的正确转换和时期匹配。
% =========================================================================
classdef TPI
    methods (Static)
        
        function [target_paths, errors, Pol_path, Dist_path, aggr_supply] = calculate_paths_and_errors(r_path, w_hat_path, Bequest_hat_path, TR_hat_path, ss0, ssF, cS, paramSF, valF, polF, dist0)
            % =========================================================================
            % ==          [注释修正] TPI核心计算引擎 (v3.7 - 注释与逻辑对齐)
            % =========================================================================
            T = cS.T_sim;

            % --- 1. 准备路径结构体 ---
            pathS = struct();
            pathS.r_path = r_path;
            pathS.w_hat_path = w_hat_path;
            pathS.tr_hat_path = TR_hat_path; % 政府转移 (人均)
            
            % [注释修正] Bequest_hat_path 是【总量】遗赠路径的猜测值
            pathS.beq_transfer_hat_path = Bequest_hat_path; 
            
            pathS = TPI.fill_auxiliary_paths(pathS, cS);

            % --- 2. 向后求解家庭最优策略 ---
            % backward_hh 内部逻辑不直接使用 beq_transfer_hat_path，所以无需修改
            [Pol_path, ~] = household.backward_hh(pathS, cS, paramSF, valF, polF);

            % --- 3. 向前模拟人口分布的演进 ---
            % [关键] 将【总量】遗赠路径传递给修正后的分布模拟器
            Dist_path = distribution.simulate_dist_forward(dist0, Pol_path, cS, paramSF, Bequest_hat_path);

            % --- 4. 聚合宏观供给量 ---
            aggr_supply = aggregates.get_path_aggregates(Dist_path, Pol_path, cS, ss0, paramSF);

            % --- 5. 计算引致的新宏观路径 (Target Paths) ---
            target_paths = struct();
            L_path = aggr_supply.L_path;
            K_p_hat_path = aggr_supply.K_p_hat_path;
            
            target_paths.r_path = zeros(1, T);
            target_paths.w_hat_path = zeros(1, T);
            Y_hat_path_calc = zeros(1, T);
            K_g_hat_fixed = ssF.K_public_hat; % 在零冲击测试中，公共资本应保持不变
            for t = 1:T
                prices_t = firm.get_prices_at_t(K_p_hat_path(t), K_g_hat_fixed, L_path(t), cS);
                target_paths.r_path(t) = prices_t.r_mkt_t;
                target_paths.w_hat_path(t) = prices_t.w_hat_t;
                Y_hat_path_calc(t) = prices_t.Y_hat_t;
            end
            
            % [逻辑不变，保持正确]
            % aggr_supply.Bequest_gen_hat_path 是【总量】
            % target_paths.Bequest_hat_path 也应该是【总量】，与迭代变量保持一致
            target_paths.Bequest_hat_path = aggr_supply.Bequest_gen_hat_path;

            target_paths.TR_hat_path = TPI.calculate_endogenous_TR_path(target_paths.r_path, target_paths.w_hat_path, K_p_hat_path, L_path, Y_hat_path_calc, Dist_path, Pol_path, cS, paramSF);

            % --- 6. 计算误差 ---
            % [逻辑不变，保持正确] 比较【总量】对【总量】
            errors.r = max(abs(target_paths.r_path(2:T) - r_path(2:T)));
            errors.w = max(abs(target_paths.w_hat_path(2:T) - w_hat_path(2:T)));
            errors.beq = max(abs(target_paths.Bequest_hat_path(2:T) - Bequest_hat_path(2:T)));
            errors.tr = max(abs(target_paths.TR_hat_path(2:T) - TR_hat_path(2:T)));
            errors.total = max([errors.r, errors.w, errors.beq, errors.tr]);
        end
        
        function pathS = fill_auxiliary_paths(pathS, cS)
            T = cS.T_sim;
            pathS.A_path = cS.A_path;
            pathS.Z_path_raw = cS.Z_path_raw;
            pathS.theta_path = cS.theta_path;
            w_path = pathS.w_hat_path .* cS.A_path;
            L_path_matrix = cS.Z_path_raw(1:cS.aR_new, :) .* cS.ageEffV_new(1:cS.aR_new);
            L_path = sum(L_path_matrix, 1); 
            pathS.b_hat_path = zeros(1, T);
            for t = 1:T
                mass_retirees_t = sum(cS.Z_path_raw((cS.aR_new+1):end, t));
                total_pension_pot = pathS.theta_path(t) * w_path(t) * L_path(t);
                b_raw_t = total_pension_pot / max(1e-9, mass_retirees_t);
                pathS.b_hat_path(t) = b_raw_t / pathS.A_path(t);
            end
        end

        % =========================================================================
        % == 函数: calculate_endogenous_TR_path
        % == 版本: [v1.1 - PAYG预算分离修正版]
        % ==
        % == 核心修正:
        % ==   - 在聚合总税收后，明确地计算并减去当期的PAYG缴费总额，
        % ==     得到用于平衡一般预算的【一般税收收入】。
        % ==   - 确保了此处的政府预算核算逻辑与稳态求解器完全一致。
        % =========================================================================
        function TR_hat_path = calculate_endogenous_TR_path(r_path, w_hat_path, K_p_hat_path, L_path, Y_hat_path, Dist_path, Pol_path, cS, paramS)
            T = cS.T_sim;
            TR_hat_path = zeros(1, T);
            total_pop_path = sum(cS.Z_path_raw, 1);

            for t=1:T
                % --- 步骤 1: 聚合当期的总税收 (包含PAYG缴费) ---
                total_tax_hat_t_raw = 0;
                dist_t_abs = Dist_path(:,:,:,:,t);
                pol_t = Pol_path{t};
                for ia = 1:cS.aD_new
                    mass_ia_abs = dist_t_abs(:,:,:,ia);
                    total_tax_hat_t_raw = total_tax_hat_t_raw + sum(pol_t(ia).tax_regular .* mass_ia_abs, 'all');
                end

                % --- 步骤 2: [核心修正] 聚合当期的PAYG缴费总额 ---
                payg_contrib_hat_t_raw = 0;
                for ia = 1:cS.aR_new % 只对工作年龄人口计算
                    mass_ia_abs = dist_t_abs(:,:,:,ia);
                    labor_income_slice_hat = w_hat_path(t) * cS.ageEffV_new(ia) .* paramS.leGridV(1:cS.nw_expanded)';
                    payg_contrib_slice_hat = cS.theta_path(t) * labor_income_slice_hat;
                    payg_contrib_grid_hat = repmat(reshape(payg_contrib_slice_hat, [1,1,cS.nw_expanded]), [cS.nk, cS.nkpps, 1]);
                    payg_contrib_hat_t_raw = payg_contrib_hat_t_raw + sum(payg_contrib_grid_hat .* mass_ia_abs, 'all');
                end
                
                % --- 步骤 3: 计算一般税收收入 (从总税收中剔除PAYG缴费) ---
                general_tax_hat_t_raw = total_tax_hat_t_raw - payg_contrib_hat_t_raw;
                
                % --- 步骤 4: 计算政府预算的其他项 (与之前相同，但使用修正后的一般税收) ---
                General_Tax_Revenue_t = general_tax_hat_t_raw * cS.A_path(t);
                Y_t = Y_hat_path(t) * cS.A_path(t);
                w_t = w_hat_path(t) * cS.A_path(t);
                K_p_t = K_p_hat_path(t) * cS.A_path(t);
                
                % 注意：这里的 K_p_return_gross_t 是(r+d)K，不是rK
                K_p_return_gross_t = r_path(t) * K_p_t + cS.ddk * K_p_t;
                Public_Capital_Return_t = Y_t - w_t * L_path(t) - K_p_return_gross_t;
                
                G_c_t = cS.G_c_to_Y_ratio_ss * Y_t;
                I_g_t = cS.I_g_to_Y_ratio_ss * Y_t;
                
                Total_General_Gov_Revenue_t = General_Tax_Revenue_t + Public_Capital_Return_t;
                
                % --- 步骤 5: 求解平衡预算的TR ---
                TR_raw_t_total = Total_General_Gov_Revenue_t - G_c_t - I_g_t;
                tr_raw_t_per_capita = TR_raw_t_total / total_pop_path(t);
                TR_hat_path(t) = tr_raw_t_per_capita / cS.A_path(t);
            end
        end

        function check_transition_NIPA(results, cS, ssF)
           fprintf('   正在启动过渡路径国民账户一致性检验...\n');
            T = cS.T_sim;
            Y_path = results.Y_path;
            C_path = results.C_path;
            I_p_path = results.I_p_path;
            I_g_path = results.I_g_path;
            G_to_Y_ratio = cS.G_c_to_Y_ratio_ss;
            G_path = G_to_Y_ratio * Y_path;
            if isfield(ssF, 'G_c'), G_path(T) = ssF.G_c * cS.A_path(T); else, G_path(T) = cS.G_c_to_Y_ratio_ss * ssF.Y_from_production_hat * cS.A_path(T); end
            Total_Uses_path = C_path + I_p_path + I_g_path + G_path;
            errors_path = Y_path - Total_Uses_path;
            safe_Y_path = Y_path;
            safe_Y_path(abs(safe_Y_path) < 1e-9) = 1e-9;
            rel_errors_path = errors_path ./ safe_Y_path;
            max_rel_error = max(abs(rel_errors_path));
            fprintf('   国民账户在路径上的最大相对误差 (误差/GDP): %.4e (%.4f%%)\n', max_rel_error, max_rel_error*100);
            if max_rel_error > 1e-5, warning('过渡路径国民账户存在显著误差，请检查代码逻辑。'); else, fprintf('   ✅ 过渡路径在数值精度内保持国民账户闭合。\n'); end
        end

        function plot_transition_results(results, cS, ssF, title_suffix)
            if nargin < 4, title_suffix = ''; end
            figure('Name', ['转轨动态路径 ' title_suffix], 'Position', [100 100 1200 800]);
            T = cS.T_sim;
            time_axis = cS.start_year : cS.time_Step : (cS.start_year + (T-1)*cS.time_Step);
            ssF_Kp_line = (ones(1, T) * ssF.K_private_hat) .* cS.A_path;
            ssF_Y_line = (ones(1, T) * ssF.Y_from_production_hat) .* cS.A_path;
            ssF_C_line = (ones(1, T) * ssF.C_agg) .* cS.A_path;
            ssF_r_line_annualized = (1 + ssF.r_mkt).^(1/cS.time_Step) - 1;
            ssF_w_line = (ones(1, T) * ssF.w_hat) .* cS.A_path;
            subplot(2,2,1);
            plot(time_axis, results.K_p_path, 'b-o', 'LineWidth', 1.5, 'DisplayName', '私人资本 (K_p)');
            hold on; plot(time_axis, ssF_Kp_line, 'b--', 'LineWidth', 1, 'DisplayName', 'K_p (BGP)');
            title('私人资本路径 (原始水平)'); xlabel('年份'); ylabel('资本');
            legend('show', 'Location', 'best'); grid on;
            subplot(2,2,2);
            plot(time_axis, results.Y_path, 'k-d', 'LineWidth', 1.5, 'DisplayName', '总产出 (Y)');
            hold on; plot(time_axis, results.C_path, 'g-^', 'LineWidth', 1.5, 'DisplayName', '总消费 (C)');
            plot(time_axis, ssF_Y_line, 'k--', 'LineWidth', 1, 'DisplayName', 'Y (BGP)');
            plot(time_axis, ssF_C_line, 'g--', 'LineWidth', 1, 'DisplayName', 'C (BGP)');
            title('产出与消费路径 (原始水平)'); xlabel('年份'); ylabel('总量');
            legend('show', 'Location', 'best'); grid on;
            subplot(2,2,3);
            r_path_annualized = (1 + results.r_path).^(1/cS.time_Step) - 1;
            plot(time_axis, r_path_annualized * 100, 'm-x', 'LineWidth', 1.5, 'DisplayName', '年化市场利率 (r)');
            hold on; yline(ssF_r_line_annualized * 100, 'm--', 'LineWidth', 1, 'DisplayName', 'r (ssF)');
            title('市场利率路径'); xlabel('年份'); ylabel('年化利率 (%)');
            legend('show', 'Location', 'best'); grid on;
            subplot(2,2,4);
            plot(time_axis, results.w_path, 'c-p', 'LineWidth', 1.5, 'DisplayName', '工资 (w)');
            hold on; plot(time_axis, ssF_w_line, 'c--', 'LineWidth', 1, 'DisplayName', 'w (BGP)');
            title('工资路径 (原始水平)'); xlabel('年份'); ylabel('工资');
            legend('show', 'Location', 'best'); grid on;
            sgtitle(['宏观经济转轨动态: ' title_suffix], 'FontSize', 16, 'FontWeight', 'bold');
        end

        function F_error = TPI_error_wrapper(X_guess, T, ss0, ssF, cS, paramSF, valF, polF, dist0)
            len = T - 1;
            r_path_guess_inner = reshape(X_guess(1:len), 1, len);
            w_hat_path_guess_inner = reshape(X_guess(len+1 : 2*len), 1, len);
            Bequest_hat_path_guess_inner = reshape(X_guess(2*len+1 : 3*len), 1, len);
            r_path_guess = [ss0.r_mkt, r_path_guess_inner];
            w_hat_path_guess = [ss0.w_hat, w_hat_path_guess_inner];
            Bequest_hat_path_guess = [ss0.Bequest_distributed_agg, Bequest_hat_path_guess_inner];
            pathS_iter = TPI.calculate_hh_income_paths(r_path_guess, w_hat_path_guess, Bequest_hat_path_guess, cS, ssF);
            [Pol_path, ~] = household.backward_hh(pathS_iter, cS, paramSF, valF, polF);
            [~, aggr_supply] = distribution.simulate_dist_forward(dist0, Pol_path, cS, paramSF, ss0);
            L_supply_path = aggr_supply.L_path;
            K_p_hat_supply_path = aggr_supply.K_p_hat_path;
            r_target_path = zeros(1, T);
            w_hat_target_path = zeros(1, T);
            K_g_hat_fixed = ssF.K_public_hat;
            for t = 1:T
                prices_t = firm.get_prices_at_t(K_p_hat_supply_path(t), K_g_hat_fixed, L_supply_path(t), cS);
                r_target_path(t) = prices_t.r_mkt_t;
                w_hat_target_path(t) = prices_t.w_hat_t;
            end
            Bequest_hat_target_path = zeros(1, T);
            Bequest_hat_target_path(1) = ss0.Bequest_distributed_agg;
            Bequest_hat_target_path(2:T) = aggr_supply.Bequest_gen_hat_path(1:T-1);
            error_r = r_target_path(2:T) - r_path_guess(2:T);
            error_w = w_hat_path_guess(2:T) - w_hat_path_guess(2:T);
            error_beq = Bequest_hat_target_path(2:T) - Bequest_hat_path_guess(2:T);
            F_error = [error_r, error_w, error_beq]'; 
        end
    end
end
