% --- aggregates.m ---
classdef aggregates
    methods (Static)
        function aggrS = get_aggregates(Dist, polS, cS, paramS, i_type)
            % =========================================================================
            % == 函数: get_aggregates
            % == 版本: [v8.0 - 会计统一最终版]
            % ==
            % == 核心修改:
            % ==   - [!!! 关键会计统一 !!!] 此函数现在是聚合的唯一来源。
            % ==   - 新增计算【期初】总资本存量 (K_p_hat_start_of_period)，
            % ==     它通过聚合状态变量 k 和 kpps 得到，代表了 K_t。
            % ==   - 原有的 K_p_hat_tplus1_raw 被重命名为 K_p_hat_end_of_period，
            % ==     它通过聚合策略变量 k_prime 和 kpps_prime 得到，代表了 K_{t+1}。
            % ==   - 移除了所有其他聚合函数，以消除逻辑不一致的根源。
            % =========================================================================
        
            K_k_hat_start_of_period = 0;
            K_kpps_hat_start_of_period = 0;
            K_k_hat_end_of_period = 0;
            K_kpps_hat_end_of_period = 0;
            C_agg_raw = 0;
            Bequest_generated_agg_raw = 0;
            Total_Tax_Revenue = 0;
            L_hat_agg = 0;

            % 为聚合构建状态变量网格
            k_grid_full = repmat(cS.kGridV, [1, cS.nkpps, cS.nw_expanded]);
            if cS.pps_active
                kpps_grid_full = repmat(reshape(cS.kppsGridV, [1, cS.nkpps]), [cS.nk, 1, cS.nw_expanded]);
            end
        
            for ia = 1:cS.aD_new
                mass_dist_ia = Dist(:,:,:,ia);
                if sum(mass_dist_ia, 'all') < 1e-30, continue; end
        
                pol_ia = polS(ia);

                % --- 1. 聚合【期初】资本存量 (K_t) ---
                K_k_hat_start_of_period = K_k_hat_start_of_period + sum(k_grid_full .* mass_dist_ia, 'all');
                if cS.pps_active
                    K_kpps_hat_start_of_period = K_kpps_hat_start_of_period + sum(kpps_grid_full .* mass_dist_ia, 'all');
                end

                % --- 2. 聚合当期流量 (消费, 税收, 劳动) ---
                C_agg_raw = C_agg_raw + sum(pol_ia.c .* mass_dist_ia, 'all');
                Total_Tax_Revenue = Total_Tax_Revenue + sum(pol_ia.tax_regular .* mass_dist_ia, 'all');
        
                if ia <= cS.aR_new
                    age_efficiency_for_type = cS.ageEff_by_type(i_type, ia);
                    labor_supply_slice = age_efficiency_for_type .* paramS.leGridV'; 
                    labor_supply_grid = repmat(reshape(labor_supply_slice, [1,1,cS.nw_expanded]), [cS.nk, cS.nkpps, 1]);
                    L_hat_agg = L_hat_agg + sum(labor_supply_grid .* mass_dist_ia, 'all');
                end
        
                % --- 3. 聚合【期末】资本存量 (K_{t+1}) ---
                K_k_hat_end_of_period = K_k_hat_end_of_period + sum(pol_ia.k_prime .* mass_dist_ia, 'all');
                if cS.pps_active
                     K_kpps_hat_end_of_period = K_kpps_hat_end_of_period + sum(pol_ia.kpps_prime .* mass_dist_ia, 'all');
                end
                
                % --- 4. 聚合当期产生的遗赠 (来自期末财富) ---
                k_prime_for_bequest = pol_ia.k_prime;
                if cS.pps_active
                    k_prime_for_bequest = k_prime_for_bequest + pol_ia.kpps_prime;
                end
                prob_death = (1 - cS.s_pathV(ia));
                Bequest_generated_agg_raw = Bequest_generated_agg_raw + sum(k_prime_for_bequest .* mass_dist_ia, 'all') * prob_death;
            end
            
            % --- 最终赋值 ---
            aggrS.K_start_of_period = K_k_hat_start_of_period + K_kpps_hat_start_of_period;
            aggrS.K_end_of_period = K_k_hat_end_of_period + K_kpps_hat_end_of_period;
            aggrS.C_agg = C_agg_raw;
            aggrS.L_hat = L_hat_agg;
            aggrS.Bequest_generated_agg = Bequest_generated_agg_raw;
            aggrS.Total_Tax_Revenue = Total_Tax_Revenue;
        end


function aggrS = get_path_aggregates(Dist_path, Pol_path, cS, ss0, paramS, dist0)
            % =========================================================================
            % ==           转轨路径宏观量聚合器 v5.1 (输出变量名修正版)
            % ==
            % == 核心修正:
            % ==   - [!!!] 将输出结构体中的人均遗赠字段名从 'Bequest_gen_hat_pc_path'
            % ==     修正为 'Bequest_gen_raw_pc_path'。
            % ==   - 这明确了其经济学含义是【未经BGP调整的人均原始遗赠产出】，
            % ==     并确保了与TPI核心引擎的变量名完全一致。
            % =========================================================================
            T = cS.T_sim;
            
            % --- 初始化输出路径 ---
            L_path = zeros(1, T);
            C_hat_path_raw = zeros(1, T);
            Bequest_gen_hat_path_raw = zeros(1, T);
            K_p_hat_total_path = zeros(1, T);

            total_pop_path = sum(cS.Z_path_raw, 1);
            
            % --- 构建用于聚合资本的网格 ---
            k_grid_full = repmat(cS.kGridV, [1, cS.nkpps, cS.nw_expanded]);
            if cS.pps_active
                kpps_grid_full = repmat(reshape(cS.kppsGridV, [1, cS.nkpps]), [cS.nk, 1, cS.nw_expanded]);
            end

            % --- 核心循环: 逐期聚合 ---
            for t = 1:T
                dist_t_abs = Dist_path(:,:,:,:,t);
                pol_t = Pol_path{t};
                
                % 1. [!!! 核心修正 !!!] 直接从期初分布聚合期初资本总量
                K_p_hat_t = 0;
                for ia = 1:cS.aD_new
                    dist_ia_t = dist_t_abs(:, :, :, ia);
                    K_p_hat_t = K_p_hat_t + sum(k_grid_full .* dist_ia_t, 'all');
                    if cS.pps_active
                         K_p_hat_t = K_p_hat_t + sum(kpps_grid_full .* dist_ia_t, 'all');
                    end
                end
                K_p_hat_total_path(t) = K_p_hat_t;

                % 2. 聚合当期的流量 (消费, 劳动, 期末产生的遗赠)
                L_t = 0; C_hat_t_raw_sum = 0; Bq_gen_hat_t_raw = 0;
                for ia = 1:cS.aD_new
                    mass_ia_abs = dist_t_abs(:,:,:,ia);
                    if sum(mass_ia_abs(:)) < 1e-30, continue; end
                    
                    pol_ia_t = pol_t(ia);
                    
                    % 聚合劳动
                    if ia <= cS.aR_new
                        le_grid_slice = reshape(paramS.leGridV, [1, 1, cS.nw_expanded]);
                        labor_supply_slice = cS.ageEffV_new(ia) .* le_grid_slice;
                        L_t = L_t + sum(labor_supply_slice .* mass_ia_abs, 'all');
                    end
                    
                    % 聚合消费
                    C_hat_t_raw_sum = C_hat_t_raw_sum + sum(pol_ia_t.c .* mass_ia_abs, 'all');
                    
                    % 聚合期末产生的遗赠
                    k_prime_slice = pol_ia_t.k_prime;
                    if cS.pps_active
                        k_prime_slice = k_prime_slice + pol_ia_t.kpps_prime; 
                    end
                    prob_death_t = (1 - cS.s_pathV(ia, t));
                    Bq_gen_hat_t_raw = Bq_gen_hat_t_raw + sum(k_prime_slice .* mass_ia_abs, 'all') * prob_death_t;
                end
                
                L_path(t) = L_t;
                C_hat_path_raw(t) = C_hat_t_raw_sum;
                Bequest_gen_hat_path_raw(t) = Bq_gen_hat_t_raw;
            end
            
            % --- 计算人均标准化路径 ---
            C_hat_pc_path = C_hat_path_raw ./ total_pop_path;
            
            % [!!! 核心修正: 统一变量名 !!!]
            Bequest_gen_raw_pc_path = zeros(1,T);
            non_zero_pop_mask = total_pop_path > 1e-9;
            Bequest_gen_raw_pc_path(non_zero_pop_mask) = Bequest_gen_hat_path_raw(non_zero_pop_mask) ./ total_pop_path(non_zero_pop_mask);


            % --- 打包输出结构体 ---
            aggrS = struct(...
                'L_path', L_path, ...
                'C_hat_pc_path', C_hat_pc_path, ...
                'K_p_hat_total_path', K_p_hat_total_path, ...
                'Bequest_gen_raw_pc_path', Bequest_gen_raw_pc_path ...
                );
end    
   
    
    end
end
