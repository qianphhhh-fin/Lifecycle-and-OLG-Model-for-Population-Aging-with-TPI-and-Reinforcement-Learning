% --- aggregates.m ---
classdef aggregates
    methods (Static)


        function [K_private_hat_agg, C_agg, Bequest_generated_agg, Total_Tax_Revenue, L_hat_agg] = get_aggregates(Dist, polS, cS, paramS)
            % =========================================================================
            % == 函数: get_aggregates
            % == 版本: [v3.4 - 劳动聚合修正版]
            % ==
            % == 核心修改:
            % ==   - 新增了对总有效劳动供给 L_hat_agg 的聚合。
            % ==   - 聚合逻辑与 get_path_aggregates 中的逻辑完全一致，
            % ==     即对完整的(k,e,a)联合分布进行加权求和。
            % ==   - 函数现在返回5个值，包括 L_hat_agg。
            % =========================================================================

            K_k_hat_agg_raw = 0;
            K_kpps_hat_agg_raw = 0;
            C_agg_raw = 0;
            Bequest_k_generated_raw = 0;
            Bequest_kpps_generated_raw = 0;
            Total_Tax_Revenue = 0;
            L_hat_agg = 0; % [!!!] 初始化劳动聚合器

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
                
                % [!!! 核心修正: 聚合劳动 !!!]
                % 仅对工作年龄段(ia <= aR_new)的个体进行聚合
                if ia <= cS.aR_new
                    % 构建一个 (nk, nkpps, nw) 维度的劳动供给网格
                    labor_supply_slice = cS.ageEffV_new(ia) .* paramS.leGridV'; % shape: [1, 1, nw]
                    labor_supply_grid = repmat(reshape(labor_supply_slice, [1,1,cS.nw_expanded]), [cS.nk, cS.nkpps, 1]);
                    
                    % 将劳动供给网格与该(e,a)状态下的个体数量(质量)相乘并累加
                    L_hat_agg = L_hat_agg + sum(labor_supply_grid .* mass_dist_ia, 'all');
                end
            end

            n_period = (1 + cS.n_ss)^cS.time_Step - 1;
            
            K_private_hat_agg_raw = K_k_hat_agg_raw + K_kpps_hat_agg_raw;
            K_private_hat_agg = K_private_hat_agg_raw / (1 + n_period);

            Bequest_generated_agg_raw = Bequest_k_generated_raw + Bequest_kpps_generated_raw;
            Bequest_generated_agg = Bequest_generated_agg_raw / (1 + n_period);

            C_agg = C_agg_raw;
        end         
        
        
        function aggrS = get_path_aggregates(Dist_path, Pol_path, cS, ss0, paramS)
            % =========================================================================
            % ==           转轨路径宏观量聚合器 v3.8 (逻辑修正最终版)
            % == 函数: get_path_aggregates
            % ==
            % == 核心修正:
            % ==   - [!!!] 删除了所有跨期调整逻辑 (如 B(t) = B_raw(t-1))。
            % ==   - 此函数的职责被严格限定为：忠实地聚合并报告在 t=1...T
            % ==     每一期中，根据当期的分布和策略实际发生的宏观总量。
            % ==   - 所有跨期核算 (如 t-1 期的储蓄/遗赠决定 t 期的资本/继承)
            % ==     都移交给 TPI.m 中的主循环处理，以保证逻辑的唯一和清晰。
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
            
            % --- 计算期初资本路径 ---
            % K_p_hat_path(t) 是 t 期初的总资本量 (hat)
            K_p_hat_path = zeros(1, T);
            % t=1 的期初资本由 t=0 (稳态) 的储蓄决定
            K_p_hat_path(1) = ss0.K_private_hat;
            % t=2...T 的期初资本 = t-1...T-1 的期末储蓄
            K_p_hat_path(2:T) = K_p_hat_end_of_period_path_raw(1:T-1);
            
            % --- [!!! 核心修正 !!!] ---
            % Bequest_gen_hat_path(t) 就是 t 期产生的总遗赠量 (hat)
            % 不再进行任何时间移位。
            Bequest_gen_hat_path = Bequest_gen_hat_path_raw;
            
            % --- 计算人均消费路径 ---
            C_hat_path = C_hat_path_raw ./ total_pop_path;

            aggrS = struct(...
                'L_path', L_path, ...
                'C_hat_path', C_hat_path, ...
                'K_p_hat_path', K_p_hat_path, ...
                'Bequest_gen_hat_path', Bequest_gen_hat_path ...
                );
        end    
    
    
    end
end
