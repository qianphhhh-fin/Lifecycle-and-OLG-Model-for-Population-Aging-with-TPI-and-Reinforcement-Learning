% --- distribution.m ---
classdef distribution
    methods (Static)

        function Dist = get_ss_dist(polS, paramS, cS, Z_ss_norm, Bequest_by_type_guess)
            % =========================================================================
            % == 函数: get_ss_dist
            % == 版本: [v4.1 - 按类型划分遗赠]
            % ==
            % == 核心修改:
            % ==   - 接收的遗赠参数从总量 `Bequest_Total_guess` 变为向量
            % ==     `Bequest_by_type_guess`。
            % ==   - 在循环中，将对应类型的遗赠总量 `Bequest_by_type_guess(i_type)`
            % ==     传递给子函数，用于计算该类型新生儿的初始资产。
            % =========================================================================

            Dist = zeros(cS.nk, cS.nkpps, cS.nw_expanded, cS.aD_new, cS.nTypes, 'double');

            for i_type = 1:cS.nTypes
                polS_type = polS{i_type};

                % [修改] 计算当前类型的绝对人口，并获取该类型产生的总遗赠
                Z_abs_type = Z_ss_norm * cS.type_weights(i_type);
                bequest_total_for_this_type = Bequest_by_type_guess(i_type);

                % [修改] 将特定类型的遗赠总量传递给子函数
                if isfield(cS, 'pps_active') && cS.pps_active
                    Dist_abs_type = distribution.get_ss_dist_pps(polS_type, paramS, cS, Z_abs_type, bequest_total_for_this_type);
                else
                    Dist_abs_type = distribution.get_ss_dist_no_pps(polS_type, paramS, cS, Z_abs_type, bequest_total_for_this_type);
                end

                Dist(:,:,:,:,i_type) = Dist_abs_type;
            end
        end

        function Dist_abs = get_ss_dist_pps(polS, paramS, cS, Z_abs, bequest_total_for_this_type)
            % =========================================================================
            % == 函数: get_ss_dist_pps
            % == 版本: [v2.1 - 按类型划分遗赠]
            % ==
            % == 核心修改:
            % ==   - 接收的遗赠参数变为 `bequest_total_for_this_type`。
            % ==   - 计算人均遗赠时，用该类型的总遗赠除以该类型的新生儿数量，
            % ==     实现了遗赠在类型内的封闭循环。
            % ==   - 输入 Z_abs 是绝对人口，输出 Dist_abs 也是绝对量分布。
            % =========================================================================
            nk = cS.nk;
            nkpps = cS.nkpps;
            nw = cS.nw_expanded;
            aD = cS.aD_new;

            mass_newborns_this_type = Z_abs(1);

            % [核心修改] 人均遗赠只在本类型内部分配
            if mass_newborns_this_type > 1e-9
                beq_per_newborn = bequest_total_for_this_type / mass_newborns_this_type;
            else
                beq_per_newborn = 0;
            end

            [ik_lower, ~, w_upper] = utils.find_grid_and_weights(beq_per_newborn, cS.kGridV);
            w_lower = 1.0 - w_upper;
            ik_upper = min(nk, ik_lower + 1);
            ikpps_newborn = 1;

            % --- 计算新生儿的绝对分布 ---
            dist_newborns_abs = zeros(nk, nkpps, nw);
            prob_mass_le = reshape(paramS.leProb1V, [1, 1, cS.nw]);
            if nw > cS.nw, prob_mass_le(1,1, (cS.nw+1):end) = 0; end

            dist_newborns_abs(ik_lower, ikpps_newborn, 1:cS.nw) = w_lower * prob_mass_le * mass_newborns_this_type;
            dist_newborns_abs(ik_upper, ikpps_newborn, 1:cS.nw) = dist_newborns_abs(ik_upper, ikpps_newborn, 1:cS.nw) + w_upper * prob_mass_le * mass_newborns_this_type;

            Dist_abs = zeros(nk, nkpps, nw, aD, 'double');
            Dist_abs(:, :, :, 1) = dist_newborns_abs;

            for ia = 1:(aD - 1)
                dist_abs_ia_slice = Dist_abs(:, :, :, ia);
                if sum(dist_abs_ia_slice(:)) < 1e-30, continue; end

                dist_abs_ia_plus_1_next = zeros(nk, nkpps, nw, 'double');
                trans_mat_next_age = paramS.TrProbM_by_age{ia};

                for ik = 1:nk
                    for ikpps = 1:nkpps
                        for ie = 1:nw
                            mass_start = dist_abs_ia_slice(ik, ikpps, ie);
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

                                dist_abs_ia_plus_1_next(ik_p_lower, ikpps_p_lower, ie_next) = dist_abs_ia_plus_1_next(ik_p_lower, ikpps_p_lower, ie_next) + mass_ll;
                                dist_abs_ia_plus_1_next(ik_p_lower, ikpps_p_upper, ie_next) = dist_abs_ia_plus_1_next(ik_p_lower, ikpps_p_upper, ie_next) + mass_lu;
                                dist_abs_ia_plus_1_next(ik_p_upper, ikpps_p_lower, ie_next) = dist_abs_ia_plus_1_next(ik_p_upper, ikpps_p_lower, ie_next) + mass_ul;
                                dist_abs_ia_plus_1_next(ik_p_upper, ikpps_p_upper, ie_next) = dist_abs_ia_plus_1_next(ik_p_upper, ikpps_p_upper, ie_next) + mass_uu;
                            end
                        end
                    end
                end
                Dist_abs(:, :, :, ia+1) = dist_abs_ia_plus_1_next;
            end
        end

        function Dist_abs = get_ss_dist_no_pps(polS, paramS, cS, Z_abs, bequest_total_for_this_type)
            % =========================================================================
            % == 函数: get_ss_dist_no_pps
            % == 版本: [v2.1 - 按类型划分遗赠]
            % == 核心修改:
            % ==   - 接收的遗赠参数变为 `bequest_total_for_this_type`。
            % ==   - 计算人均遗赠时，用该类型的总遗赠除以该类型的新生儿数量。
            % =========================================================================
            nk = cS.nk;
            nw = cS.nw_expanded;
            aD = cS.aD_new;

            mass_newborns_this_type = Z_abs(1);

            % [核心修改] 人均遗赠只在本类型内部分配
            if mass_newborns_this_type > 1e-9
                beq_per_newborn = bequest_total_for_this_type / mass_newborns_this_type;
            else
                beq_per_newborn = 0;
            end

            [ik_lower, ~, w_upper] = utils.find_grid_and_weights(beq_per_newborn, cS.kGridV);
            w_lower = 1.0 - w_upper;
            ik_upper = min(nk, ik_lower + 1);

            dist_newborns_abs = zeros(nk, 1, nw);
            prob_mass_le = reshape(paramS.leProb1V, [1, 1, cS.nw]);
            if nw > cS.nw, prob_mass_le(1,1, (cS.nw+1):end) = 0; end

            dist_newborns_abs(ik_lower, 1, 1:cS.nw) = w_lower * prob_mass_le * mass_newborns_this_type;
            dist_newborns_abs(ik_upper, 1, 1:cS.nw) = dist_newborns_abs(ik_upper, 1, 1:cS.nw) + w_upper * prob_mass_le * mass_newborns_this_type;

            Dist_abs = zeros(nk, 1, nw, aD, 'double');
            Dist_abs(:, :, :, 1) = dist_newborns_abs;

            for ia = 1:(aD - 1)
                dist_abs_ia_slice = Dist_abs(:, 1, :, ia);
                if sum(dist_abs_ia_slice(:)) < 1e-30, continue; end

                dist_abs_ia_plus_1_next = zeros(nk, 1, nw, 'double');
                trans_mat_next_age = paramS.TrProbM_by_age{ia};

                for ik = 1:nk
                    for ie = 1:nw
                        mass_start = dist_abs_ia_slice(ik, 1, ie);
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
                            dist_abs_ia_plus_1_next(ik_lower_p, 1, ie_next) = dist_abs_ia_plus_1_next(ik_lower_p, 1, ie_next) + mass_to_distribute * w_k_lower;
                            dist_abs_ia_plus_1_next(ik_upper_p, 1, ie_next) = dist_abs_ia_plus_1_next(ik_upper_p, 1, ie_next) + mass_to_distribute * w_k_upper;
                        end
                    end
                end
                Dist_abs(:, :, :, ia+1) = dist_abs_ia_plus_1_next;
            end
        end

        function Dist_path = simulate_dist_forward(dist0, Pol_path, cS, paramS, bequest_per_newborn_hat_path)
            % =========================================================================
            % == 函数: simulate_dist_forward (v4.1 - 新生儿会计与循环修正版)
            % == 核心修正:
            % ==   - [!!! 关键BUG修复 !!!] 将计算新生儿分布 (dist_newborns_tp1)
            % ==     的代码块从 for ia...end 循环中移出。新生儿的计算应该在所有
            % ==     幸存者都演进完毕后，独立地、一次性地进行。
            % ==   - [!!! 接口语义修正 !!!] 函数现在接收一个更清晰的输入
            % ==     `bequest_per_newborn_hat_path`，它代表【每个新生儿】在
            % ==     t+1 期初收到的【个体标准化资产】。这消除了函数内部所有
            % ==     易错的、关于总量和人均量的转换。
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

                % --- 步骤1: 演进所有幸存者到 t+1 期 ---
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
                end % End of ia loop

                % --- 步骤2: [!!! 核心修正 !!!] 在幸存者演进循环之外，独立计算 t+1 期的新生儿分布 ---
                dist_newborns_tp1 = zeros(nk, nkpps, nw, aD);
                mass_newborns_tp1_abs = Z_path_abs(1, t+1);

                if mass_newborns_tp1_abs > 1e-9
                    % 直接使用传入的、为 t+1 期新生儿准备的个体资产
                    beq_per_newborn_hat_tp1 = bequest_per_newborn_hat_path(t+1);

                    [ik_lower, ~, w_upper] = utils.find_grid_and_weights(beq_per_newborn_hat_tp1, cS.kGridV);
                    w_lower = 1.0 - w_upper;
                    ik_upper = min(nk, ik_lower + 1);
                    ikpps_newborn = 1;

                    prob_mass_le = reshape(paramS.leProb1V, [1, 1, cS.nw]);

                    dist_newborns_tp1(ik_lower, ikpps_newborn, 1:cS.nw, 1) = w_lower * prob_mass_le * mass_newborns_tp1_abs;
                    dist_newborns_tp1(ik_upper, ikpps_newborn, 1:cS.nw, 1) = dist_newborns_tp1(ik_upper, ikpps_newborn, 1:cS.nw, 1) + w_upper * prob_mass_le * mass_newborns_tp1_abs;
                end

                % --- 步骤3: 合并幸存者和新生儿，得到 t+1 期的完整分布 ---
                Dist_path(:,:,:,:,t+1) = dist_survivors_tp1 + dist_newborns_tp1;

            end % End of t loop
        end


    end
end
