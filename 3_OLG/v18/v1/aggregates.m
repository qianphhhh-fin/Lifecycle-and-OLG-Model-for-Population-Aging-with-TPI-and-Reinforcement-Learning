% --- aggregates.m ---
classdef aggregates
    methods (Static)


        function [K_private_hat_agg, C_agg, Bequest_generated_agg, Total_Tax_Revenue] = get_aggregates(Dist, polS, cS, paramS)
            % =========================================================================
            % == 函数: get_aggregates
            % == 版本: [v3.1 - PPS 聚合最终版]
            % ==
            % == 核心修改:
            % ==   - 确保在无PPS时，kpps相关聚合为0，代码健壮。
            % =========================================================================

            K_k_hat_agg_raw = 0;
            K_kpps_hat_agg_raw = 0;
            C_agg_raw = 0;
            Bequest_k_generated_raw = 0;
            Bequest_kpps_generated_raw = 0;
            Total_Tax_Revenue = 0;

            for ia = 1:cS.aD_new
                mass_dist_ia = Dist(:,:,:,ia);
                C_agg_raw = C_agg_raw + sum(polS(ia).c .* mass_dist_ia, 'all');
                Total_Tax_Revenue = Total_Tax_Revenue + sum(polS(ia).tax_regular .* mass_dist_ia, 'all');

                K_k_hat_agg_raw = K_k_hat_agg_raw + sum(polS(ia).k_prime .* mass_dist_ia, 'all');
                prob_death = (1 - cS.s_pathV(ia));
                Bequest_k_generated_raw = Bequest_k_generated_raw + sum(polS(ia).k_prime .* mass_dist_ia, 'all') * prob_death;

                % [PPS自适应聚合]
                if isfield(cS, 'pps_active') && cS.pps_active
                    K_kpps_hat_agg_raw = K_kpps_hat_agg_raw + sum(polS(ia).kpps_prime .* mass_dist_ia, 'all');
                    Bequest_kpps_generated_raw = Bequest_kpps_generated_raw + sum(polS(ia).kpps_prime .* mass_dist_ia, 'all') * prob_death;
                end
            end

            n_period = (1 + cS.n_ss)^cS.time_Step - 1;

            K_private_hat_agg_raw = K_k_hat_agg_raw + K_kpps_hat_agg_raw;
            K_private_hat_agg = K_private_hat_agg_raw / (1 + n_period);

            Bequest_generated_agg_raw = Bequest_k_generated_raw + Bequest_kpps_generated_raw;
            Bequest_generated_agg = Bequest_generated_agg_raw / (1 + n_period);

            C_agg = C_agg_raw;
        end


        % =========================================================================
        % ==           [核心修正] 转轨路径宏观量聚合器 v3.4
        % == 函数: get_path_aggregates
        % ==
        % == 目的:
        % ==   - [!!!] 将遗赠路径 (Bequest_gen_hat_path) 的输出修正为【总量】，
        % ==     与TPI求解器的迭代变量定义保持一致。
        % ==   - 资本路径 K_p_hat_path 保持为【总量】。
        % ==   - 消费 C_hat_path 保持为【人均】。
        % =========================================================================
        function aggrS = get_path_aggregates(Dist_path, Pol_path, cS, ss0, paramS)
            T = cS.T_sim;
            aD = cS.aD_new;

            L_path = zeros(1, T);
            C_hat_path = zeros(1, T);
            Bequest_gen_hat_path_raw = zeros(1, T);
            K_p_hat_end_of_period_path_raw = zeros(1, T);

            Z_path_abs = cS.Z_path_raw;
            total_pop_path = sum(Z_path_abs, 1);

            g_n_path = zeros(1, T);
            % 在完美的零冲击测试中, g_n_path(1) 应该等于 g_n_ss
            g_n_period_ss = (1+cS.n_ss)^cS.time_Step-1;
            g_n_path(1) = g_n_period_ss;
            g_n_path(2:T) = total_pop_path(2:T) ./ total_pop_path(1:T-1) - 1;

            for t = 1:T
                dist_t_abs = Dist_path(:,:,:,:,t);
                pol_t = Pol_path{t};
                L_t = 0; C_hat_t_raw = 0; Kp_hat_eop_t_raw = 0; Bq_gen_hat_t_raw = 0;

                for ia = 1:aD
                    mass_ia_abs = dist_t_abs(:,:,:,ia);
                    pol_ia_t = pol_t(ia);
                    if ia <= cS.aR_new
                        le_grid_slice = reshape(paramS.leGridV, [1,1,cS.nw_expanded]);
                        labor_supply_slice = cS.ageEffV_new(ia) .* le_grid_slice;
                        L_t = L_t + sum(labor_supply_slice .* mass_ia_abs, 'all');
                    end
                    C_hat_t_raw = C_hat_t_raw + sum(pol_ia_t.c .* mass_ia_abs, 'all');
                    k_prime_slice = pol_ia_t.k_prime;
                    if cS.pps_active, k_prime_slice = k_prime_slice + pol_ia_t.kpps_prime; end
                    Kp_hat_eop_t_raw = Kp_hat_eop_t_raw + sum(k_prime_slice .* mass_ia_abs, 'all');
                    prob_death_t = (1 - cS.s_pathV(ia, t));
                    Bq_gen_hat_t_raw = Bq_gen_hat_t_raw + sum(k_prime_slice .* mass_ia_abs, 'all') * prob_death_t;
                end

                L_path(t) = L_t;
                C_hat_path(t) = C_hat_t_raw / total_pop_path(t);
                K_p_hat_end_of_period_path_raw(t) = Kp_hat_eop_t_raw;
                Bequest_gen_hat_path_raw(t) = Bq_gen_hat_t_raw;
            end

            K_p_hat_path = zeros(1, T);
            K_p_hat_path(1) = ss0.K_private_hat;
            K_p_hat_path(2:T) = K_p_hat_end_of_period_path_raw(1:T-1) ./ (1 + g_n_path(2:T));

            % [核心修正] 遗赠路径是【总量】
            Bequest_gen_hat_path = zeros(1,T);
            Bequest_gen_hat_path(1) = ss0.Bequest_distributed_agg;
            % t期产生的遗赠总量，经过人口增长规范化后，成为t+1期的可用遗赠总量
            Bequest_gen_hat_path(2:T) = Bequest_gen_hat_path_raw(1:T-1) ./ (1 + g_n_path(2:T));

            aggrS = struct(...
                'L_path', L_path, ...
                'C_hat_path', C_hat_path, ...
                'K_p_hat_path', K_p_hat_path, ...
                'Bequest_gen_hat_path', Bequest_gen_hat_path ... % 返回【总量】遗赠
                );
        end


    end
end
