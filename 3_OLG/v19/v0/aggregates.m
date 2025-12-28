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
        
        
        function aggrS = get_path_aggregates(Dist_path, Pol_path, cS, ss0, paramS, dist0)
            % =========================================================================
            % ==           转轨路径宏观量聚合器 v4.1 (初始资本会计修正版)
            % == 函数: get_path_aggregates
            % ==
            % == 核心修正:
            % ==   - [BUG修复] 增加了输入参数 `dist0`，以修复变量无法识别的错误。
            % ==   - [!!! 关键会计修正 !!!] 正确计算t=1期的初始资本总量。
            % ==     通过将稳态的标准化资本(ss0.K_private_hat)乘以初始绝对总人口
            % ==     (sum(dist0, 'all')), 确保了从稳态到转轨的无缝衔接。
            % ==   - 输出的人均量现在基于正确的总量计算，保证了内部一致性。
            % =========================================================================
            T = cS.T_sim;
            aD = cS.aD_new;

            L_path = zeros(1, T);
            C_hat_path_raw = zeros(1, T);
            Bequest_gen_hat_path_raw = zeros(1, T);
            K_p_hat_end_of_period_path_raw = zeros(1, T);

            total_pop_path = sum(cS.Z_path_raw, 1);

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
            
            % --- 计算期初总资本路径 (Total Standardized Capital) ---
            K_p_hat_total_path = zeros(1, T);
            % [!!! 核心会计修正 !!!] t=1的期初资本 = 标准化稳态资本 * 初始绝对总人口
            K_p_hat_total_path(1) = ss0.K_private_hat * sum(dist0, 'all'); 
            % t>1的期初资本 = (t-1)的期末储蓄总量
            K_p_hat_total_path(2:T) = K_p_hat_end_of_period_path_raw(1:T-1);
            
            % --- 计算人均标准化路径 (Per-Capita Standardized Paths) ---
            % 人均消费 = 总消费 / 总人口
            C_hat_pc_path = C_hat_path_raw ./ total_pop_path;
            % 人均遗赠 = 总遗赠 / 总人口
            Bequest_gen_hat_pc_path = Bequest_gen_hat_path_raw ./ total_pop_path;

            aggrS = struct(...
                'L_path', L_path, ...
                'C_hat_pc_path', C_hat_pc_path, ...
                'K_p_hat_total_path', K_p_hat_total_path, ...
                'Bequest_gen_hat_pc_path', Bequest_gen_hat_pc_path ...
                );
        end
    
    
    end
end
