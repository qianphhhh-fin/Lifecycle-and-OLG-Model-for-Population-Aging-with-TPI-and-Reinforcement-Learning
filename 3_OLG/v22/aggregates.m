% --- aggregates.m ---
classdef aggregates
    methods (Static)

        function aggrS = get_aggregates_from_dist(Dist, cS)
            % =========================================================================
            % ==             【新增】从分布聚合当期资本的函数 v1.0
            % ==
            % == 功能:
            % ==   - 仅接收一个分布(Dist)，聚合出该分布代表的【当期期初】总资本存量。
            % ==   - 这是最直接、最无歧义的资本存量定义。
            % =========================================================================
            K_p_hat_total = 0;

            % 构建用于聚合【当期】资本的网格 (状态变量 k, kpps)
            k_grid_full = repmat(cS.kGridV, [1, cS.nkpps, cS.nw_expanded]);
            if cS.pps_active
                kpps_grid_full = repmat(reshape(cS.kppsGridV, [1, cS.nkpps]), [cS.nk, 1, cS.nw_expanded]);
            end

            for ia = 1:cS.aD_new
                mass_dist_ia = Dist(:,:,:,ia);
                if sum(mass_dist_ia, 'all') < 1e-30, continue; end

                K_p_hat_total = K_p_hat_total + sum(k_grid_full .* mass_dist_ia, 'all');
                if cS.pps_active
                    K_p_hat_total = K_p_hat_total + sum(kpps_grid_full .* mass_dist_ia, 'all');
                end
            end

            aggrS.K_p_hat_total = K_p_hat_total;
        end

        function aggrS = get_aggregates(Dist, polS, cS, paramS)
            % =========================================================================
            % == 函数: get_aggregates
            % == 版本: [v7.0 - 职责明确版]
            % ==
            % == 核心逻辑:
            % ==   - 职责1: 聚合当期的流量变量 (C, L, Tax)。
            % ==   - 职责2: 聚合【下一期】的资本存量 K_t+1_raw (来自策略 k_prime)。
            % ==   - 职责3: 聚合当期死亡者产生的【期末】遗赠。
            % ==   - [!!!] 此函数自身【不】进行任何BGP调整。
            % =========================================================================

            K_k_hat_agg_raw = 0;
            K_kpps_hat_agg_raw = 0;
            C_agg_raw = 0;
            Bequest_generated_agg_raw = 0;
            Total_Tax_Revenue = 0;
            L_hat_agg = 0;

            for ia = 1:cS.aD_new
                mass_dist_ia = Dist(:,:,:,ia);
                if sum(mass_dist_ia, 'all') < 1e-30, continue; end

                pol_ia = polS(ia);

                % --- 聚合当期流量 (消费, 税收, 劳动) ---
                C_agg_raw = C_agg_raw + sum(pol_ia.c .* mass_dist_ia, 'all');
                Total_Tax_Revenue = Total_Tax_Revenue + sum(pol_ia.tax_regular .* mass_dist_ia, 'all');

                if ia <= cS.aR_new
                    labor_supply_slice = cS.ageEffV_new(ia) .* paramS.leGridV';
                    labor_supply_grid = repmat(reshape(labor_supply_slice, [1,1,cS.nw_expanded]), [cS.nk, cS.nkpps, 1]);
                    L_hat_agg = L_hat_agg + sum(labor_supply_grid .* mass_dist_ia, 'all');
                end

                % --- 聚合下一期初的资本存量 K_t+1_raw (来自策略变量 k_prime) ---
                K_k_hat_agg_raw = K_k_hat_agg_raw + sum(pol_ia.k_prime .* mass_dist_ia, 'all');
                if cS.pps_active
                    K_kpps_hat_agg_raw = K_kpps_hat_agg_raw + sum(pol_ia.kpps_prime .* mass_dist_ia, 'all');
                end

                % --- 聚合当期产生的遗赠 (来自策略变量 k_prime) ---
                k_prime_for_bequest = pol_ia.k_prime;
                if cS.pps_active
                    k_prime_for_bequest = k_prime_for_bequest + pol_ia.kpps_prime;
                end
                prob_death = (1 - cS.s_pathV(ia));
                Bequest_generated_agg_raw = Bequest_generated_agg_raw + sum(k_prime_for_bequest .* mass_dist_ia, 'all') * prob_death;
            end

            % --- 最终赋值 ---
            aggrS.K_p_hat_tplus1_raw = K_k_hat_agg_raw + K_kpps_hat_agg_raw;
            aggrS.C_agg = C_agg_raw;
            aggrS.L_hat = L_hat_agg;
            aggrS.Bequest_generated_agg = Bequest_generated_agg_raw;
            aggrS.Total_Tax_Revenue = Total_Tax_Revenue;
        end

        % 脚本: aggregates.m
        function aggrS = get_path_aggregates(Dist_path_h, Pol_path_h, cS, paramS)
            % =========================================================================
            % ==           转轨路径宏观量聚合器 v7.0 (时变退休年龄版)
            % ==
            % == 核心修改:
            % ==   - [!!!] 在聚合劳动供给时，不再使用固定的 cS.aR_new。
            % ==   - 在时间循环 t 内部，首先获取当期退休年龄索引 aR_t = cS.aR_new_path(t)。
            % ==   - 使用 aR_t 作为年龄循环的上界来判断一个年龄组是否提供劳动。
            % ==   - 这确保了总劳动供给的计算与家庭在VFI中面对的时变退休年龄完全一致。
            % =========================================================================
            T = cS.T_sim;
            nH = cS.nTypes;

            % --- 初始化输出路径 ---
            L_path = zeros(1, T);
            C_hat_total_path = zeros(1, T);
            K_p_hat_total_path = zeros(1, T);
            Bequest_gen_raw_total_path = zeros(1, T);
            Bequest_gen_raw_pc_path_h = zeros(nH, T); % 分类型人均

            total_pop_by_type_path = zeros(nH, T);
            for h = 1:nH
                total_pop_by_type_path(h,:) = sum(cS.Z_path_raw,1) * cS.type_weights(h);
            end

            k_grid_full = repmat(cS.kGridV, [1, cS.nkpps, cS.nw_expanded]);
            if cS.pps_active
                kpps_grid_full = repmat(reshape(cS.kppsGridV, [1, cS.nkpps]), [cS.nk, 1, cS.nw_expanded]);
            end

            % --- 核心循环: 逐期聚合 ---
            for t = 1:T
                % === 核心修改: 获取当期的退休年龄索引 ===
                aR_t = cS.aR_new_path(t);

                % --- 在每个时期内，循环所有家庭类型 h ---
                for h = 1:nH
                    cS_h = cS; % 获取特定类型的参数
                    cS_h.ageEffV_new = cS.ageEffV_new_h(:,h);

                    dist_t_h = Dist_path_h{h}(:,:,:,:,t);
                    pol_t_h = Pol_path_h{h}{t};

                    % 1. 直接从期初分布聚合期初资本总量
                    K_p_hat_t_h = 0;
                    for ia = 1:cS.aD_new
                        dist_ia_t_h = dist_t_h(:, :, :, ia);
                        K_p_hat_t_h = K_p_hat_t_h + sum(k_grid_full .* dist_ia_t_h, 'all');
                        if cS.pps_active
                            K_p_hat_t_h = K_p_hat_t_h + sum(kpps_grid_full .* dist_ia_t_h, 'all');
                        end
                    end
                    K_p_hat_total_path(t) = K_p_hat_total_path(t) + K_p_hat_t_h;

                    % 2. 聚合当期的流量 (消费, 劳动, 期末产生的遗赠)
                    L_t_h = 0; C_hat_t_h = 0; Bq_gen_hat_t_h = 0;
                    for ia = 1:cS.aD_new
                        mass_ia_abs_h = dist_t_h(:,:,:,ia);
                        if sum(mass_ia_abs_h(:)) < 1e-30, continue; end

                        pol_ia_t_h = pol_t_h(ia);

                        % === 核心修改: 使用 aR_t 判断是否提供劳动 ===
                        if ia <= aR_t
                            le_grid_slice = reshape(paramS.leGridV, [1, 1, cS.nw_expanded]);
                            labor_supply_slice = cS_h.ageEffV_new(ia) .* le_grid_slice;
                            L_t_h = L_t_h + sum(labor_supply_slice .* mass_ia_abs_h, 'all');
                        end

                        C_hat_t_h = C_hat_t_h + sum(pol_ia_t_h.c .* mass_ia_abs_h, 'all');

                        k_prime_slice = pol_ia_t_h.k_prime;
                        if cS.pps_active
                            k_prime_slice = k_prime_slice + pol_ia_t_h.kpps_prime;
                        end
                        prob_death_t = (1 - cS.s_pathV(ia, t));
                        Bq_gen_hat_t_h = Bq_gen_hat_t_h + sum(k_prime_slice .* mass_ia_abs_h, 'all') * prob_death_t;
                    end

                    % 累加到总量
                    L_path(t) = L_path(t) + L_t_h;
                    C_hat_total_path(t) = C_hat_total_path(t) + C_hat_t_h;
                    Bequest_gen_raw_total_path(t) = Bequest_gen_raw_total_path(t) + Bq_gen_hat_t_h;

                    % 计算分类型人均遗赠
                    pop_h_t = total_pop_by_type_path(h,t);
                    if pop_h_t > 1e-9
                        Bequest_gen_raw_pc_path_h(h, t) = Bq_gen_hat_t_h / pop_h_t;
                    end
                end % 结束 h 循环
            end % 结束 t 循环

            total_pop_path = sum(cS.Z_path_raw, 1);
            C_hat_pc_path = C_hat_total_path ./ total_pop_path;

            % --- 打包输出结构体 ---
            aggrS = struct(...
                'L_path', L_path, ...
                'C_hat_total_path', C_hat_total_path, ...
                'C_hat_pc_path', C_hat_pc_path, ...
                'K_p_hat_total_path', K_p_hat_total_path, ...
                'Bequest_gen_raw_total_path', Bequest_gen_raw_total_path, ...
                'Bequest_gen_raw_pc_path_h', Bequest_gen_raw_pc_path_h ...
                );
        end
        % function aggrS = get_path_aggregates(Dist_path_h, Pol_path_h, cS, paramS)
        %     % =========================================================================
        %     % ==           转轨路径宏观量聚合器 v6.0 (异质性家庭版)
        %     % ==
        %     % == 核心修改:
        %     % ==   - [!!!] 函数接口完全重构，以处理异质性家庭。
        %     % ==   - 输入变为单元数组: Dist_path_h, Pol_path_h。
        %     % ==   - 输出变为结构体，包含总量路径和分类型路径。
        %     % ==   - 内部逻辑增加了一个对家庭类型 h 的循环。
        %     % ==   - K_p_hat_total_path: 所有类型资本的总和。
        %     % ==   - L_path: 所有类型劳动供给的总和。
        %     % ==   - C_hat_total_path: 所有类型消费的总和。
        %     % ==   - Bequest_gen_raw_total_path: 所有类型产生遗赠的总和。
        %     % ==   - Bequest_gen_raw_pc_path_h: [nH x T]矩阵，各类型人均遗赠。
        %     % =========================================================================
        %     T = cS.T_sim;
        %     nH = cS.nTypes;
        % 
        %     % --- 初始化输出路径 ---
        %     L_path = zeros(1, T);
        %     C_hat_total_path = zeros(1, T);
        %     K_p_hat_total_path = zeros(1, T);
        %     Bequest_gen_raw_total_path = zeros(1, T);
        %     Bequest_gen_raw_pc_path_h = zeros(nH, T); % 分类型人均
        % 
        %     total_pop_by_type_path = zeros(nH, T);
        %     for h = 1:nH
        %         total_pop_by_type_path(h,:) = sum(cS.Z_path_raw,1) * cS.type_weights(h);
        %     end
        % 
        %     k_grid_full = repmat(cS.kGridV, [1, cS.nkpps, cS.nw_expanded]);
        %     if cS.pps_active
        %         kpps_grid_full = repmat(reshape(cS.kppsGridV, [1, cS.nkpps]), [cS.nk, 1, cS.nw_expanded]);
        %     end
        % 
        %     % --- 核心循环: 逐期聚合 ---
        %     for t = 1:T
        %         % --- 在每个时期内，循环所有家庭类型 h ---
        %         for h = 1:nH
        %             cS_h = cS; % 获取特定类型的参数
        %             cS_h.ageEffV_new = cS.ageEffV_new_h(:,h);
        % 
        %             dist_t_h = Dist_path_h{h}(:,:,:,:,t);
        %             pol_t_h = Pol_path_h{h}{t};
        % 
        %             % 1. 直接从期初分布聚合期初资本总量
        %             K_p_hat_t_h = 0;
        %             for ia = 1:cS.aD_new
        %                 dist_ia_t_h = dist_t_h(:, :, :, ia);
        %                 K_p_hat_t_h = K_p_hat_t_h + sum(k_grid_full .* dist_ia_t_h, 'all');
        %                 if cS.pps_active
        %                     K_p_hat_t_h = K_p_hat_t_h + sum(kpps_grid_full .* dist_ia_t_h, 'all');
        %                 end
        %             end
        %             K_p_hat_total_path(t) = K_p_hat_total_path(t) + K_p_hat_t_h;
        % 
        %             % 2. 聚合当期的流量 (消费, 劳动, 期末产生的遗赠)
        %             L_t_h = 0; C_hat_t_h = 0; Bq_gen_hat_t_h = 0;
        %             for ia = 1:cS.aD_new
        %                 mass_ia_abs_h = dist_t_h(:,:,:,ia);
        %                 if sum(mass_ia_abs_h(:)) < 1e-30, continue; end
        % 
        %                 pol_ia_t_h = pol_t_h(ia);
        % 
        %                 if ia <= cS.aR_new
        %                     le_grid_slice = reshape(paramS.leGridV, [1, 1, cS.nw_expanded]);
        %                     labor_supply_slice = cS_h.ageEffV_new(ia) .* le_grid_slice;
        %                     L_t_h = L_t_h + sum(labor_supply_slice .* mass_ia_abs_h, 'all');
        %                 end
        % 
        %                 C_hat_t_h = C_hat_t_h + sum(pol_ia_t_h.c .* mass_ia_abs_h, 'all');
        % 
        %                 k_prime_slice = pol_ia_t_h.k_prime;
        %                 if cS.pps_active
        %                     k_prime_slice = k_prime_slice + pol_ia_t_h.kpps_prime;
        %                 end
        %                 prob_death_t = (1 - cS.s_pathV(ia, t));
        %                 Bq_gen_hat_t_h = Bq_gen_hat_t_h + sum(k_prime_slice .* mass_ia_abs_h, 'all') * prob_death_t;
        %             end
        % 
        %             % 累加到总量
        %             L_path(t) = L_path(t) + L_t_h;
        %             C_hat_total_path(t) = C_hat_total_path(t) + C_hat_t_h;
        %             Bequest_gen_raw_total_path(t) = Bequest_gen_raw_total_path(t) + Bq_gen_hat_t_h;
        % 
        %             % 计算分类型人均遗赠
        %             pop_h_t = total_pop_by_type_path(h,t);
        %             if pop_h_t > 1e-9
        %                 Bequest_gen_raw_pc_path_h(h, t) = Bq_gen_hat_t_h / pop_h_t;
        %             end
        %         end % 结束 h 循环
        %     end % 结束 t 循环
        % 
        %     total_pop_path = sum(cS.Z_path_raw, 1);
        %     C_hat_pc_path = C_hat_total_path ./ total_pop_path;
        % 
        %     % --- 打包输出结构体 ---
        %     aggrS = struct(...
        %         'L_path', L_path, ...
        %         'C_hat_total_path', C_hat_total_path, ...
        %         'C_hat_pc_path', C_hat_pc_path, ...
        %         'K_p_hat_total_path', K_p_hat_total_path, ...
        %         'Bequest_gen_raw_total_path', Bequest_gen_raw_total_path, ...
        %         'Bequest_gen_raw_pc_path_h', Bequest_gen_raw_pc_path_h ...
        %         );
        % end



    end
end
