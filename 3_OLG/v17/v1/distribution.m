% --- distribution.m ---
classdef distribution
    methods (Static)

        function Dist = get_ss_dist(polS, paramS, cS, Z_ss_norm, Bequest_Total_guess)
            % =========================================================================
            % == 函数: get_ss_dist
            % == 版本: [v3.0 - PPS 兼容版]
            % ==
            % == 核心修改:
            % ==   - 根据 cS.pps_active 状态，决定是处理一维(k)还是二维(k, kpps)的分布。
            % ==   - 引入新的迭代逻辑来处理二维资产的转移。
            % =========================================================================

            if isfield(cS, 'pps_active') && cS.pps_active
                Dist = distribution.get_ss_dist_pps(polS, paramS, cS, Z_ss_norm, Bequest_Total_guess);
            else
                Dist = distribution.get_ss_dist_no_pps(polS, paramS, cS, Z_ss_norm, Bequest_Total_guess);
            end
        end

        function Dist = get_ss_dist_pps(polS, paramS, cS, Z_ss_norm, Bequest_Total_guess)
            % --- 计算带PPS的稳态分布 ---
            nk = cS.nk;
            nkpps = cS.nkpps;
            nw = cS.nw_expanded;
            aD = cS.aD_new;

            mass_newborns = Z_ss_norm(1);
            if mass_newborns > 1e-9
                beq_per_newborn = Bequest_Total_guess / mass_newborns;
            else
                beq_per_newborn = 0;
            end

            % [遗赠分配逻辑] 新生儿 kpps=0, 所有遗赠进入常规资产 k
            [ik_lower, ~, w_upper] = utils.find_grid_and_weights(beq_per_newborn, cS.kGridV);
            w_lower = 1.0 - w_upper;
            ik_upper = min(nk, ik_lower + 1);
            ikpps_newborn = 1; % 新生儿PPS资产为0，在网格点1

            dist_newborn_cond_shape = zeros(nk, nkpps, nw);
            prob_mass_to_distribute = reshape(paramS.leProb1V, [1, 1, cS.nw]);
            if nw > cS.nw, prob_mass_to_distribute(1,1, (cS.nw+1):end) = 0; end % 冲击状态不给新生儿

            dist_newborn_cond_shape(ik_lower, ikpps_newborn, 1:cS.nw) = w_lower * prob_mass_to_distribute;
            dist_newborn_cond_shape(ik_upper, ikpps_newborn, 1:cS.nw) = dist_newborn_cond_shape(ik_upper, ikpps_newborn, 1:cS.nw) + w_upper * prob_mass_to_distribute;
            
            Dist_cond = zeros(nk, nkpps, nw, aD, 'double');
            Dist_cond(:, :, :, 1) = dist_newborn_cond_shape;
            
            % [分布迭代]
            for ia = 1:(aD - 1)
                dist_cond_ia_slice = Dist_cond(:, :, :, ia);
                if sum(dist_cond_ia_slice(:)) < 1e-30, continue; end
                
                dist_cond_ia_plus_1_next = zeros(nk, nkpps, nw, 'double');
                trans_mat_next_age = paramS.TrProbM_by_age{ia + 1};
                
                for ik = 1:nk
                    for ikpps = 1:nkpps
                        for ie = 1:nw
                            mass_start = dist_cond_ia_slice(ik, ikpps, ie);
                            if mass_start < 1e-30, continue; end
                            
                            % 获取下一期资产选择
                            k_prime = polS(ia).k_prime(ik, ikpps, ie);
                            kpps_prime = polS(ia).kpps_prime(ik, ikpps, ie);
                            
                            % 在两个资产维度上进行插值
                            [ik_p_lower, ~, w_k_upper] = utils.find_grid_and_weights(k_prime, cS.kGridV);
                            w_k_lower = 1.0 - w_k_upper;
                            ik_p_upper = min(nk, ik_p_lower + 1);
                            
                            [ikpps_p_lower, ~, w_kpps_upper] = utils.find_grid_and_weights(kpps_prime, cS.kppsGridV);
                            w_kpps_lower = 1.0 - w_kpps_upper;
                            ikpps_p_upper = min(nkpps, ikpps_p_lower + 1);

                            trans_probs_vec = trans_mat_next_age(ie, :);
                            for ie_next = 1:nw
                                prob_to_enext = trans_probs_vec(ie_next);
                                if prob_to_enext < 1e-12, continue; end
                                
                                mass_to_distribute = mass_start * prob_to_enext;
                                
                                % 将质量分配到四个相邻的网格点
                                mass_ll = mass_to_distribute * w_k_lower * w_kpps_lower;
                                mass_lu = mass_to_distribute * w_k_lower * w_kpps_upper;
                                mass_ul = mass_to_distribute * w_k_upper * w_kpps_lower;
                                mass_uu = mass_to_distribute * w_k_upper * w_kpps_upper;

                                dist_cond_ia_plus_1_next(ik_p_lower, ikpps_p_lower, ie_next) = dist_cond_ia_plus_1_next(ik_p_lower, ikpps_p_lower, ie_next) + mass_ll;
                                dist_cond_ia_plus_1_next(ik_p_lower, ikpps_p_upper, ie_next) = dist_cond_ia_plus_1_next(ik_p_lower, ikpps_p_upper, ie_next) + mass_lu;
                                dist_cond_ia_plus_1_next(ik_p_upper, ikpps_p_lower, ie_next) = dist_cond_ia_plus_1_next(ik_p_upper, ikpps_p_lower, ie_next) + mass_ul;
                                dist_cond_ia_plus_1_next(ik_p_upper, ikpps_p_upper, ie_next) = dist_cond_ia_plus_1_next(ik_p_upper, ikpps_p_upper, ie_next) + mass_uu;
                            end
                        end
                    end
                end
                sum_next_dist = sum(dist_cond_ia_plus_1_next(:));
                if sum_next_dist > 1e-9
                    Dist_cond(:, :, :, ia+1) = dist_cond_ia_plus_1_next / sum_next_dist;
                end
            end
            
            Dist = Dist_cond .* reshape(Z_ss_norm, [1, 1, 1, aD]);
        end

        function Dist = get_ss_dist_no_pps(polS, paramS, cS, Z_ss_norm, Bequest_Total_guess)
             % --- 这是您原始的 get_ss_dist 函数，重命名为 no_pps 版本 ---
            nk = cS.nk;
            nw = cS.nw_expanded;
            aD = cS.aD_new;

            mass_newborns = Z_ss_norm(1);

            if mass_newborns > 1e-9
                beq_per_newborn = Bequest_Total_guess / mass_newborns;
            else
                beq_per_newborn = 0;
            end

            [ik_lower, ~, w_upper] = utils.find_grid_and_weights(beq_per_newborn, cS.kGridV);
            w_lower = 1.0 - w_upper;
            ik_upper = min(nk, ik_lower + 1);

            dist_newborn_cond_shape = zeros(nk, 1, nw);
            prob_mass_to_distribute = reshape(paramS.leProb1V, [1, 1, cS.nw]);
            if nw > cS.nw, prob_mass_to_distribute(1,1, (cS.nw+1):end) = 0; end
            dist_newborn_cond_shape(ik_lower, 1, 1:cS.nw) = w_lower * prob_mass_to_distribute;
            dist_newborn_cond_shape(ik_upper, 1, 1:cS.nw) = dist_newborn_cond_shape(ik_upper, 1, 1:cS.nw) + w_upper * prob_mass_to_distribute;
            
            Dist_cond = zeros(nk, 1, nw, aD, 'double');
            Dist_cond(:, :, :, 1) = dist_newborn_cond_shape;
            
            for ia = 1:(aD - 1)
                dist_cond_ia_slice = Dist_cond(:, 1, :, ia);
                if sum(dist_cond_ia_slice(:)) < 1e-30, continue; end
                
                dist_cond_ia_plus_1_next = zeros(nk, 1, nw, 'double');
                trans_mat_next_age = paramS.TrProbM_by_age{ia + 1};
                
                for ik = 1:nk
                    for ie = 1:nw
                        mass_start = dist_cond_ia_slice(ik, 1, ie);
                        if mass_start < 1e-30, continue; end
                        
                        k_prime = polS(ia).k_prime(ik, 1, ie);
                        [ik_lower_p, ~, w_k_upper] = utils.find_grid_and_weights(k_prime, cS.kGridV);
                        w_k_lower = 1.0 - w_k_upper;
                        ik_upper_p = min(nk, ik_lower_p + 1);
                        
                        trans_probs_vec = trans_mat_next_age(ie, :);
                        for ie_next = 1:nw
                            prob_to_enext = trans_probs_vec(ie_next);
                            if prob_to_enext < 1e-12, continue; end
                            
                            mass_to_distribute = mass_start * prob_to_enext;
                            dist_cond_ia_plus_1_next(ik_lower_p, 1, ie_next) = dist_cond_ia_plus_1_next(ik_lower_p, 1, ie_next) + mass_to_distribute * w_k_lower;
                            dist_cond_ia_plus_1_next(ik_upper_p, 1, ie_next) = dist_cond_ia_plus_1_next(ik_upper_p, 1, ie_next) + mass_to_distribute * w_k_upper;
                        end
                    end
                end
                sum_next_dist = sum(dist_cond_ia_plus_1_next(:));
                if sum_next_dist > 1e-9
                    Dist_cond(:, :, :, ia+1) = dist_cond_ia_plus_1_next / sum_next_dist;
                end
            end
            
            Dist = Dist_cond .* reshape(Z_ss_norm, [1, 1, 1, aD]);
        end

    end
end