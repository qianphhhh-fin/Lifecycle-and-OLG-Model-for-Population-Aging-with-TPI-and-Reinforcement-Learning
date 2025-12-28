% --- aggregates.m ---
classdef aggregates
    methods (Static)


        function aggrS = get_aggregates_from_dist(Dist, cS)
            % =========================================================================
            % ==             【新增】从分布聚合当期资本的函数 v1.1
            % ==
            % == 功能:
            % ==   - 修正了原函数无法处理5维(含nTypes)分布的BUG。
            % ==   - 现在函数会正确地遍历所有家庭类型，并加总其资本存量。
            % ==   - 这确保了求解器看到的总资本供给与经济体中所有家庭的
            % ==     实际储蓄之和相一致。
            % =========================================================================
            K_p_hat_total = 0;

            % 构建用于聚合【当期】资本的网格 (状态变量 k, kpps)
            k_grid_full = repmat(cS.kGridV, [1, cS.nkpps, cS.nw_expanded]);
            if cS.pps_active
                kpps_grid_full = repmat(reshape(cS.kppsGridV, [1, cS.nkpps]), [cS.nk, 1, cS.nw_expanded]);
            end

            for i_type = 1:cS.nTypes
                Dist_type = Dist(:,:,:,:,i_type); % 取出当前类型的分布
                for ia = 1:cS.aD_new
                    mass_dist_ia = Dist_type(:,:,:,ia);
                    if sum(mass_dist_ia, 'all') < 1e-30, continue; end

                    K_p_hat_total = K_p_hat_total + sum(k_grid_full .* mass_dist_ia, 'all');
                    if cS.pps_active
                        K_p_hat_total = K_p_hat_total + sum(kpps_grid_full .* mass_dist_ia, 'all');
                    end
                end
            end

            aggrS.K_p_hat_total = K_p_hat_total;
        end
        
        
        function aggrS = get_aggregates(Dist, polS, cS, paramS)
            % =========================================================================
            % == 函数: get_aggregates
            % == 版本: [v8.2 - 聚合逻辑修正版]
            % ==
            % == 核心修改:
            % ==   - 新增输出 `L_hat_by_type`，用于在一次聚合中同时得到
            % ==     按类型划分的劳动供给。
            % ==   - 这避免了在求解器中对同一变量进行重复和可能不一致的计算，
            % ==     是消除求解误差的关键步骤。
            % =========================================================================

            K_k_hat_agg_raw = 0;
            K_kpps_hat_agg_raw = 0;
            C_agg_raw = 0;

            Bequest_generated_by_type_agg = zeros(cS.nTypes, 1);
            Total_Tax_Revenue = 0;
            L_hat_agg = 0;
            
            % [新增] 初始化按类型划分的劳动供给聚合器
            L_hat_by_type = zeros(cS.nTypes, 1);

            for ia = 1:cS.aD_new
                for i_type = 1:cS.nTypes
                    mass_dist_ia_itype = Dist(:,:,:,ia, i_type);
                    if sum(mass_dist_ia_itype, 'all') < 1e-30, continue; end

                    pol_ia_type = polS{i_type}(ia);

                    C_agg_raw = C_agg_raw + sum(pol_ia_type.c .* mass_dist_ia_itype, 'all');
                    Total_Tax_Revenue = Total_Tax_Revenue + sum(pol_ia_type.tax_regular .* mass_dist_ia_itype, 'all');

                    if ia <= cS.aR_new
                        age_eff_type = cS.ageEff_by_type(i_type, ia);
                        labor_supply_slice = age_eff_type .* paramS.leGridV';
                        labor_supply_grid = repmat(reshape(labor_supply_slice, [1,1,cS.nw_expanded]), [cS.nk, cS.nkpps, 1]);
                        
                        % [新增] 累加到特定类型的劳动供给中
                        labor_supply_this_group = sum(labor_supply_grid .* mass_dist_ia_itype, 'all');
                        L_hat_by_type(i_type) = L_hat_by_type(i_type) + labor_supply_this_group;
                    end

                    K_k_hat_agg_raw = K_k_hat_agg_raw + sum(pol_ia_type.k_prime .* mass_dist_ia_itype, 'all');
                    if cS.pps_active
                        K_kpps_hat_agg_raw = K_kpps_hat_agg_raw + sum(pol_ia_type.kpps_prime .* mass_dist_ia_itype, 'all');
                    end

                    k_prime_for_bequest = pol_ia_type.k_prime;
                    if cS.pps_active
                        k_prime_for_bequest = k_prime_for_bequest + pol_ia_type.kpps_prime;
                    end

                    prob_death = (1 - cS.s_pathV(ia));

                    bequest_this_group = sum(k_prime_for_bequest .* mass_dist_ia_itype, 'all') * prob_death;
                    Bequest_generated_by_type_agg(i_type) = Bequest_generated_by_type_agg(i_type) + bequest_this_group;

                end
            end
            L_hat_agg = sum(L_hat_by_type); % 总劳动供给现在是分类型加总的结果
            
            % --- 最终赋值 ---
            aggrS.K_p_hat_tplus1_raw = K_k_hat_agg_raw + K_kpps_hat_agg_raw;
            aggrS.C_agg = C_agg_raw;
            aggrS.L_hat = L_hat_agg;
            
            % [新增] 增加新的分类型输出
            aggrS.L_hat_by_type = L_hat_by_type;
            aggrS.Bequest_generated_by_type_agg = Bequest_generated_by_type_agg;
            aggrS.Bequest_generated_agg = sum(Bequest_generated_by_type_agg);

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
