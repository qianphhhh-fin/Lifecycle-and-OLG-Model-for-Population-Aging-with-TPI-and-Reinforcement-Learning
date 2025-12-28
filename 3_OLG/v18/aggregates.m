% --- aggregates.m ---
classdef aggregates
    methods (Static)


        function [K_private_hat_agg, C_agg, Bequest_generated_agg, Total_Tax_Revenue] = get_aggregates(Dist, polS, cS, paramS)
            % =========================================================================
            % == 函数: get_aggregates
            % == 版本: [v3.3 - 稳态BGP聚合修正版]
            % ==
            % == 核心修改:
            % ==   - 恢复了在计算稳态资本和遗赠时，对 (1 + n_period) 的除法。
            % ==   - 理由: 在稳态(BGP)中，t期人口产生的期末总储蓄
            % ==     K_private_hat_agg_raw (去趋势A)，必须等于t+1期人口需要的
            % ==     期初总资本。而t+1期的总资本(去趋势A)是 t期的 (1+n) 倍。
            % ==     即 K_hat_{t+1} = K_hat_t * (1+n)。
            % ==     所以，为了让市场出清，必须有 K_hat_t = K_private_hat_agg_raw / (1+n)。
            % ==     这个修正确保了家庭部门的储蓄行为与BGP的投资需求一致。
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
            
            % [!!! 核心修正 !!!] 
            % 在稳态，t期产生的期末总储蓄 K_raw (去趋势A)，要供给 t+1 期人口使用。
            % t+1 期的资本存量 K_hat_{t+1} (去趋势A) 是 t 期的 K_hat_t * (1+n)。
            % 为了使市场出清 K_hat_t = K_raw / (1+n)。
            K_private_hat_agg_raw = K_k_hat_agg_raw + K_kpps_hat_agg_raw;
            K_private_hat_agg = K_private_hat_agg_raw / (1 + n_period);

            % [!!! 核心修正 !!!] 遗赠同理。
            Bequest_generated_agg_raw = Bequest_k_generated_raw + Bequest_kpps_generated_raw;
            Bequest_generated_agg = Bequest_generated_agg_raw / (1 + n_period);

            C_agg = C_agg_raw;
        end
         
        
        
        function aggrS = get_path_aggregates(Dist_path, Pol_path, cS, ss0, paramS)
            % =========================================================================
            % ==           [核心修正] 转轨路径宏观量聚合器 v3.7
            % == 函数: get_path_aggregates
            % ==
            % == 目的:
            % ==   - [!!!] 将宏观总量 K_p_hat_path 和 Bequest_gen_hat_path 定义为
            % ==     【总量】(仅技术去趋势)，与 L_path 的定义保持一致。
            % ==   - 依据会计恒等式：t-1期末产生的总资本/遗赠，直接成为 t 期初
            % ==     可用的总量。不再进行任何人口增长的除法调整。
            % ==   - 这确保了在n>0的BGP上，K 和 L 能以相同速率增长，从而稳定K/L比率和价格。
            % =========================================================================
            T = cS.T_sim;
            aD = cS.aD_new;

            L_path = zeros(1, T);
            C_hat_path_raw = zeros(1, T);
            Bequest_gen_hat_path_raw = zeros(1, T);
            K_p_hat_end_of_period_path_raw = zeros(1, T);

            Z_path_abs = cS.Z_path_raw;
            total_pop_path = sum(Z_path_abs, 1);

            for t = 1:T
                dist_t_abs = Dist_path(:,:,:,:,t);
                pol_t = Pol_path{t};
                L_t = 0; C_hat_t_raw_sum = 0; Kp_hat_eop_t_raw = 0; Bq_gen_hat_t_raw = 0;

                for ia = 1:aD
                    mass_ia_abs = dist_t_abs(:,:,:,ia);
                    pol_ia_t = pol_t(ia);
                    if ia <= cS.aR_new
                        le_grid_slice = reshape(paramS.leGridV, [1,1,cS.nw_expanded]);
                        labor_supply_slice = cS.ageEffV_new(ia) .* le_grid_slice;
                        L_t = L_t + sum(labor_supply_slice .* mass_ia_abs, 'all');
                    end
                    C_hat_t_raw_sum = C_hat_t_raw_sum + sum(pol_ia_t.c .* mass_ia_abs, 'all');
                    k_prime_slice = pol_ia_t.k_prime;
                    if cS.pps_active, k_prime_slice = k_prime_slice + pol_ia_t.kpps_prime; end
                    Kp_hat_eop_t_raw = Kp_hat_eop_t_raw + sum(k_prime_slice .* mass_ia_abs, 'all');
                    prob_death_t = (1 - cS.s_pathV(ia, t));
                    Bq_gen_hat_t_raw = Bq_gen_hat_t_raw + sum(k_prime_slice .* mass_ia_abs, 'all') * prob_death_t;
                end

                L_path(t) = L_t;
                C_hat_path_raw(t) = C_hat_t_raw_sum;
                K_p_hat_end_of_period_path_raw(t) = Kp_hat_eop_t_raw;
                Bequest_gen_hat_path_raw(t) = Bq_gen_hat_t_raw;
            end
            
            C_hat_path = C_hat_path_raw ./ total_pop_path;
            
            % [!!! 核心修正 !!!]
            % 定义 K_p_hat_path 和 Bequest_gen_hat_path 为【总量】
            % K_p_hat_path(t) 是 t 期初的资本总量。
            % K_p_hat_end_of_period_path_raw(t-1) 是 t-1 期末产生的资本总量。
            % 会计恒等式：t期初的资本 = t-1期末的储蓄。
            
            K_p_hat_path = zeros(1, T);
            Bequest_gen_hat_path = zeros(1, T);

            % t=1 的资本由 t=0 (稳态)的决策决定。
            % ss0.K_private_hat 是BGP规范化的稳态值。在t=1，总人口尚未增长，
            % 因此这个值就是t=1的资本总量。
            K_p_hat_path(1) = ss0.K_private_hat; 
            Bequest_gen_hat_path(1) = ss0.Bequest_distributed_agg;

            % t=2...T 时，t期初的资本总量 = t-1期末产生的资本总量
            K_p_hat_path(2:T) = K_p_hat_end_of_period_path_raw(1:T-1);
            Bequest_gen_hat_path(2:T) = Bequest_gen_hat_path_raw(1:T-1);

            aggrS = struct(...
                'L_path', L_path, ...
                'C_hat_path', C_hat_path, ...
                'K_p_hat_path', K_p_hat_path, ...
                'Bequest_gen_hat_path', Bequest_gen_hat_path ...
                );
        end    
    
    
    
    end
end
