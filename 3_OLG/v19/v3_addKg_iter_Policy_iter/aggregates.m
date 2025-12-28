% --- aggregates.m ---
classdef aggregates
    methods (Static)
        % function [K_private_hat_agg, C_agg, Bequest_generated_agg, Total_Tax_Revenue, L_hat_agg] = get_aggregates(Dist, polS, cS, paramS)
        %     % =========================================================================
        %     % == 函数: get_aggregates
        %     % == 版本: [v4.0 - 资本聚合口径与TPI对齐修正版]
        %     % ==
        %     % == 核心修正:
        %     % ==   - [!!! 根本性BUG修复 !!!] 彻底改变资本聚合方式，与TPI的逻辑看齐。
        %     % ==   - 不再聚合策略变量 k'，而是直接聚合状态变量 k。
        %     % ==   - 通过构建与分布匹配的资本网格，直接对期初分布进行加权求和，
        %     % ==     从而得到【当期期初】的资本总量 K_t。
        %     % ==   - 移除了错误的BGP调整 "/(1+n_period)"，因其混淆了存量定义。
        %     % =========================================================================
        % 
        %     C_agg_raw = 0;
        %     Bequest_generated_agg_raw = 0;
        %     Total_Tax_Revenue = 0;
        %     L_hat_agg = 0;
        %     K_private_hat_agg = 0; % 初始化【当期】资本聚合器
        % 
        %     % --- 构建用于聚合【当期】资本的网格 ---
        %     % 这些网格代表了状态变量 (k, kpps)
        %     k_grid_full = repmat(cS.kGridV, [1, cS.nkpps, cS.nw_expanded]);
        %     if cS.pps_active
        %         kpps_grid_full = repmat(reshape(cS.kppsGridV, [1, cS.nkpps]), [cS.nk, 1, cS.nw_expanded]);
        %     end
        % 
        %     for ia = 1:cS.aD_new
        %         mass_dist_ia = Dist(:,:,:,ia);
        %         if sum(mass_dist_ia, 'all') < 1e-30, continue; end
        % 
        %         % --- [核心修正] 聚合【当期期初】资本存量 (K_t) ---
        %         K_private_hat_agg = K_private_hat_agg + sum(k_grid_full .* mass_dist_ia, 'all');
        %         if cS.pps_active
        %             K_private_hat_agg = K_private_hat_agg + sum(kpps_grid_full .* mass_dist_ia, 'all');
        %         end
        % 
        %         % --- 聚合当期流量 (消费, 税收, 劳动) ---
        %         C_agg_raw = C_agg_raw + sum(polS(ia).c .* mass_dist_ia, 'all');
        %         Total_Tax_Revenue = Total_Tax_Revenue + sum(polS(ia).tax_regular .* mass_dist_ia, 'all');
        % 
        %         if ia <= cS.aR_new
        %             labor_supply_slice = cS.ageEffV_new(ia) .* paramS.leGridV'; 
        %             labor_supply_grid = repmat(reshape(labor_supply_slice, [1,1,cS.nw_expanded]), [cS.nk, cS.nkpps, 1]);
        %             L_hat_agg = L_hat_agg + sum(labor_supply_grid .* mass_dist_ia, 'all');
        %         end
        % 
        %         % --- 聚合【当期期末产生】的遗赠 (Bequest_t) ---
        %         k_prime_slice = polS(ia).k_prime;
        %         if cS.pps_active
        %             k_prime_slice = k_prime_slice + polS(ia).kpps_prime;
        %         end
        %         prob_death = (1 - cS.s_pathV(ia));
        %         Bequest_generated_agg_raw = Bequest_generated_agg_raw + sum(k_prime_slice .* mass_dist_ia, 'all') * prob_death;
        %     end
        % 
        %     % --- 最终赋值，不进行任何BGP调整，因为所有量都是基于当期分布的真实聚合 ---
        %     C_agg = C_agg_raw;
        %     Bequest_generated_agg = Bequest_generated_agg_raw;
        % end

        function [K_private_hat_agg, C_agg, Bequest_generated_agg, Total_Tax_Revenue, L_hat_agg] = get_aggregates(Dist, polS, cS, paramS)
    % =========================================================================
    % == 函数: get_aggregates
    % == 版本: [v6.0 - 稳态专用版，逻辑修复]
    % ==
    % == 核心逻辑:
    % ==   - 此函数专为【稳态求解器】服务。
    % ==   - 它聚合【策略变量 k_prime】，代表所有家庭在期末意愿持有的总资产。
    % ==   - 这是下一期期初的【未经BGP调整】的资本总量 K_t+1_raw。
    % ==   - [!!!] 输出的 K_private_hat_agg 就是这个未经调整的 K_t+1_raw。
    % ==     BGP调整将在调用此函数的上层 (SS.m) 中进行。
    % =========================================================================

    K_k_hat_agg_raw = 0;
    K_kpps_hat_agg_raw = 0;
    C_agg_raw = 0;
    Bequest_generated_agg_raw = 0;
    Total_Tax_Revenue = 0;
    L_hat_agg = 0;

    for ia = 1:cS.aD_new
        mass_dist_ia = Dist(:,:,:,ia);
        C_agg_raw = C_agg_raw + sum(polS(ia).c .* mass_dist_ia, 'all');
        Total_Tax_Revenue = Total_Tax_Revenue + sum(polS(ia).tax_regular .* mass_dist_ia, 'all');

        % --- 聚合下一期初的资本存量 (来自策略变量 k_prime) ---
        K_k_hat_agg_raw = K_k_hat_agg_raw + sum(polS(ia).k_prime .* mass_dist_ia, 'all');
        if cS.pps_active
             K_kpps_hat_agg_raw = K_kpps_hat_agg_raw + sum(polS(ia).kpps_prime .* mass_dist_ia, 'all');
        end
        
        % --- 聚合当期产生的遗赠 (来自策略变量 k_prime) ---
        k_prime_for_bequest = polS(ia).k_prime;
        if cS.pps_active
            k_prime_for_bequest = k_prime_for_bequest + polS(ia).kpps_prime;
        end
        prob_death = (1 - cS.s_pathV(ia));
        Bequest_generated_agg_raw = Bequest_generated_agg_raw + sum(k_prime_for_bequest .* mass_dist_ia, 'all') * prob_death;

        % --- 聚合劳动 ---
        if ia <= cS.aR_new
            labor_supply_slice = cS.ageEffV_new(ia) .* paramS.leGridV';
            labor_supply_grid = repmat(reshape(labor_supply_slice, [1,1,cS.nw_expanded]), [cS.nk, cS.nkpps, 1]);
            L_hat_agg = L_hat_agg + sum(labor_supply_grid .* mass_dist_ia, 'all');
        end
    end
    
    % --- 最终赋值 ---
    % K_private_hat_agg 是未经调整的下一期资本总量
    K_private_hat_agg = K_k_hat_agg_raw + K_kpps_hat_agg_raw;
    C_agg = C_agg_raw;
    Bequest_generated_agg = Bequest_generated_agg_raw;
end        
        



% 脚本: aggregates.m

function aggrS = get_path_aggregates(Dist_path, Pol_path, cS, ss0, paramS, dist0)
            % =========================================================================
            % ==           转轨路径宏观量聚合器 v5.2 (总量输出修正版)
            % ==
            % == 核心修正:
            % ==   - 输出结构体中的遗赠字段明确为 'Bequest_gen_hat_path'，
            % ==     代表【未经BGP调整的有效单位遗赠总量】。
            % ==   - 移除了不再需要的人均量计算，使函数职责更单一。
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
                
                % 1. 直接从期初分布聚合期初资本总量
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
                    
                    % 聚合消费 (总量, hat单位)
                    C_hat_t_raw_sum = C_hat_t_raw_sum + sum(pol_ia_t.c .* mass_ia_abs, 'all');
                    
                    % 聚合期末产生的遗赠 (总量, hat单位)
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
            
            % --- 计算人均标准化消费路径 ---
            C_hat_pc_path = C_hat_path_raw ./ total_pop_path;

            % --- 打包输出结构体 ---
            aggrS = struct(...
                'L_path', L_path, ...
                'C_hat_pc_path', C_hat_pc_path, ...
                'K_p_hat_total_path', K_p_hat_total_path, ...
                'Bequest_gen_hat_path', Bequest_gen_hat_path_raw ... % [!!!] 输出总量
                );
end
    
    
    end
end
