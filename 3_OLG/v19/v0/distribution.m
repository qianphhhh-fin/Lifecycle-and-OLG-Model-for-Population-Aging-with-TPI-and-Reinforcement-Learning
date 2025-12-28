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
            % =========================================================================
            % == 函数: get_ss_dist_pps (v1.1 - 冲击矩阵索引修正版)
            % == 核心修正:
            % ==   - [!!!] 修正了访问年龄相关冲击转移矩阵时的索引错误。
            % ==   - 原代码使用 {ia + 1}，现修正为 {ia}。
            % =========================================================================
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

            [ik_lower, ~, w_upper] = utils.find_grid_and_weights(beq_per_newborn, cS.kGridV);
            w_lower = 1.0 - w_upper;
            ik_upper = min(nk, ik_lower + 1);
            ikpps_newborn = 1;

            dist_newborn_cond_shape = zeros(nk, nkpps, nw);
            prob_mass_to_distribute = reshape(paramS.leProb1V, [1, 1, cS.nw]);
            if nw > cS.nw, prob_mass_to_distribute(1,1, (cS.nw+1):end) = 0; end

            dist_newborn_cond_shape(ik_lower, ikpps_newborn, 1:cS.nw) = w_lower * prob_mass_to_distribute;
            dist_newborn_cond_shape(ik_upper, ikpps_newborn, 1:cS.nw) = dist_newborn_cond_shape(ik_upper, ikpps_newborn, 1:cS.nw) + w_upper * prob_mass_to_distribute;

            Dist_cond = zeros(nk, nkpps, nw, aD, 'double');
            Dist_cond(:, :, :, 1) = dist_newborn_cond_shape;

            for ia = 1:(aD - 1)
                dist_cond_ia_slice = Dist_cond(:, :, :, ia);
                if sum(dist_cond_ia_slice(:)) < 1e-30, continue; end

                dist_cond_ia_plus_1_next = zeros(nk, nkpps, nw, 'double');
                % [!!! 核心修正 !!!] 使用正确的索引 {ia}
                trans_mat_next_age = paramS.TrProbM_by_age{ia};

                for ik = 1:nk
                    for ikpps = 1:nkpps
                        for ie = 1:nw
                            mass_start = dist_cond_ia_slice(ik, ikpps, ie);
                            if mass_start < 1e-30, continue; end

                            k_prime = polS(ia).k_prime(ik, ikpps, ie);
                            kpps_prime = polS(ia).kpps_prime(ik, ikpps, ie);

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
            % =========================================================================
            % == 函数: get_ss_dist_no_pps (v1.1 - 冲击矩阵索引修正版)
            % == 核心修正:
            % ==   - [!!!] 修正了访问年龄相关冲击转移矩阵时的索引错误。
            % ==   - 原代码使用 {ia + 1}，现修正为 {ia}。
            % =========================================================================
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
                % [!!! 核心修正 !!!] 使用正确的索引 {ia}
                trans_mat_next_age = paramS.TrProbM_by_age{ia};

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



        function Dist_path = simulate_dist_forward(dist0, Pol_path, cS, paramS, bequest_per_newborn_hat_path)
            % =========================================================================
            % == 函数: simulate_dist_forward (v3.9 - 新生儿会计修正版)
            % == 核心修正:
            % ==   - [!!!] 修正了最关键的会计错误。函数现在接收一个清晰定义的
            % ==     `bequest_per_newborn_hat_path` 路径，该路径代表【每个新生儿】
            % ==     在对应时期收到的【个体标准化资产】。
            % ==   - 删除了所有内部多余且错误的关于总量和人均量的转换计算。
            % ==   - 函数现在直接、干净地使用传入的个体资产量为新生儿进行插值。
            % =========================================================================
            T = cS.T_sim;
            nk = cS.nk;
            nkpps = cS.nkpps;
            nw = cS.nw_expanded;
            aD = cS.aD_new;

            Dist_path = zeros(nk, nkpps, nw, aD, T);
            Dist_path(:,:,:,:,1) = dist0;
            Z_path_abs = cS.Z_path_raw;

            for t = 1:(T-1)
                dist_t_abs = Dist_path(:,:,:,:,t);
                pol_t = Pol_path{t};
                s_path_t = cS.s_pathV(:, t);

                dist_survivors_tp1 = zeros(nk, nkpps, nw, aD);

                for ia = 1:(aD - 1)
                    dist_ia_t = dist_t_abs(:,:,:,ia);
                    if sum(dist_ia_t(:)) < 1e-30, continue; end

                    pol_ia_t = pol_t(ia);
                    trans_mat_next_age = paramS.TrProbM_by_age{ia};

                    for ik = 1:nk
                        for ikpps = 1:nkpps
                            for ie = 1:nw
                                mass_start_abs = dist_ia_t(ik, ikpps, ie);
                                if mass_start_abs < 1e-30, continue; end

                                mass_survived_abs = mass_start_abs * s_path_t(ia);
                                if mass_survived_abs < 1e-30, continue; end

                                k_prime_hat = pol_ia_t.k_prime(ik, ikpps, ie);
                                [ik_p_lower, ~, w_k_upper] = utils.find_grid_and_weights(k_prime_hat, cS.kGridV);
                                w_k_lower = 1.0 - w_k_upper;
                                ik_p_upper = min(nk, ik_p_lower + 1);

                                trans_probs_vec = trans_mat_next_age(ie, :);
                                for ie_next = 1:nw
                                    prob_to_enext = trans_probs_vec(ie_next);
                                    if prob_to_enext < 1e-12, continue; end

                                    mass_to_distribute_abs = mass_survived_abs * prob_to_enext;

                                    if cS.pps_active
                                        kpps_prime_hat = pol_ia_t.kpps_prime(ik, ikpps, ie);
                                        [ikpps_p_lower, ~, w_kpps_upper] = utils.find_grid_and_weights(kpps_prime_hat, cS.kppsGridV);
                                        w_kpps_lower = 1.0 - w_kpps_upper;
                                        ikpps_p_upper = min(nkpps, ikpps_p_lower + 1);

                                        dist_survivors_tp1(ik_p_lower, ikpps_p_lower, ie_next, ia+1) = dist_survivors_tp1(ik_p_lower, ikpps_p_lower, ie_next, ia+1) + mass_to_distribute_abs * w_k_lower * w_kpps_lower;
                                        dist_survivors_tp1(ik_p_lower, ikpps_p_upper, ie_next, ia+1) = dist_survivors_tp1(ik_p_lower, ikpps_p_upper, ie_next, ia+1) + mass_to_distribute_abs * w_k_lower * w_kpps_upper;
                                        dist_survivors_tp1(ik_p_upper, ikpps_p_lower, ie_next, ia+1) = dist_survivors_tp1(ik_p_upper, ikpps_p_lower, ie_next, ia+1) + mass_to_distribute_abs * w_k_upper * w_kpps_lower;
                                        dist_survivors_tp1(ik_p_upper, ikpps_p_upper, ie_next, ia+1) = dist_survivors_tp1(ik_p_upper, ikpps_p_upper, ie_next, ia+1) + mass_to_distribute_abs * w_k_upper * w_kpps_upper;
                                    else
                                        dist_survivors_tp1(ik_p_lower, 1, ie_next, ia+1) = dist_survivors_tp1(ik_p_lower, 1, ie_next, ia+1) + mass_to_distribute_abs * w_k_lower;
                                        dist_survivors_tp1(ik_p_upper, 1, ie_next, ia+1) = dist_survivors_tp1(ik_p_upper, 1, ie_next, ia+1) + mass_to_distribute_abs * w_k_upper;
                                    end
                                end
                            end
                        end
                    end

                    dist_newborns_tp1 = zeros(nk, nkpps, nw, aD);
                    mass_newborns_tp1_abs = Z_path_abs(1, t+1);

                    if mass_newborns_tp1_abs > 1e-9
                        % [!!! 核心会计修正: 直接使用传入的、干净的“每个新生儿收到的资产” !!!]
                        beq_per_newborn_hat_tp1 = bequest_per_newborn_hat_path(t+1); % t+1的新生儿接收t+1时期的值

                        [ik_lower, ~, w_upper] = utils.find_grid_and_weights(beq_per_newborn_hat_tp1, cS.kGridV);
                        w_lower = 1.0 - w_upper;
                        ik_upper = min(nk, ik_lower + 1);
                        ikpps_newborn = 1;

                        prob_mass_le = reshape(paramS.leProb1V, [1, 1, cS.nw]);

                        dist_newborns_tp1(ik_lower, ikpps_newborn, 1:cS.nw, 1) = w_lower * prob_mass_le * mass_newborns_tp1_abs;
                        dist_newborns_tp1(ik_upper, ikpps_newborn, 1:cS.nw, 1) = dist_newborns_tp1(ik_upper, ikpps_newborn, 1:cS.nw, 1) + w_upper * prob_mass_le * mass_newborns_tp1_abs;
                    end

                    Dist_path(:,:,:,:,t+1) = dist_survivors_tp1 + dist_newborns_tp1;
                end
            end
        end

    end
end
