% --- distribution.m ---
classdef distribution
    methods (Static)

               function Dist = get_ss_dist(polS, paramS, cS, Z_ss_norm, bq_per_capita)
            % =========================================================================
            % == 函数: get_ss_dist
            % == 版本: [v3.1 - 异质性框架接口修正版]
            % ==
            % == 核心修改:
            % ==   - 将遗赠输入参数名从 Bequest_Total_guess 修改为 bq_per_capita，
            % ==     以明确其经济含义是人均量。
            % =========================================================================

            if isfield(cS, 'pps_active') && cS.pps_active
                Dist = distribution.get_ss_dist_pps(polS, paramS, cS, Z_ss_norm, bq_per_capita);
            else
                Dist = distribution.get_ss_dist_no_pps(polS, paramS, cS, Z_ss_norm, bq_per_capita);
            end
        end

        function Dist = get_ss_dist_no_pps(polS, paramS, cS, Z_ss_norm, bq_per_capita)
    % =========================================================================
    % == 函数: get_ss_dist_no_pps (v1.3 - 生存率逻辑修正版)
    % == 核心修正:
    % ==   - [!!! 关键BUG修复 !!!] 与 get_ss_dist_pps 的修改一致，
    % ==     在计算稳态分布的代际演化时，引入了生存率 cS.s_pathV(ia)。
    % =========================================================================
    nk = cS.nk;
    nw = cS.nw_expanded;
    aD = cS.aD_new;

    mass_newborns = Z_ss_norm(1);

    if mass_newborns > 1e-9
        total_bequest_h = bq_per_capita * sum(Z_ss_norm);
        beq_per_newborn = total_bequest_h / mass_newborns;
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
        trans_mat_next_age = paramS.TrProbM_by_age{ia};

        for ik = 1:nk
            for ie = 1:nw
                mass_start = dist_cond_ia_slice(ik, 1, ie);
                if mass_start < 1e-30, continue; end
                
                % [!!! 核心修正: 引入生存率 !!!]
                mass_survived = mass_start * cS.s_pathV(ia);
                if mass_survived < 1e-30, continue; end

                k_prime = polS(ia).k_prime(ik, 1, ie);
                [ik_lower_p, ~, w_k_upper] = utils.find_grid_and_weights(k_prime, cS.kGridV);
                w_k_lower = 1.0 - w_k_upper;
                ik_upper_p = min(nk, ik_lower_p + 1);

                trans_probs_vec = trans_mat_next_age(ie, :);
                for ie_next = 1:nw
                    prob_to_enext = trans_probs_vec(ie_next);
                    if prob_to_enext < 1e-12, continue; end
                    
                    % 使用存活下来的质量进行分配
                    mass_to_distribute = mass_survived * prob_to_enext;
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

        function Dist = get_ss_dist_pps(polS, paramS, cS, Z_ss_norm, bq_per_capita)
    % =========================================================================
    % == 函数: get_ss_dist_pps (v1.3 - 生存率逻辑修正版)
    % == 核心修正:
    % ==   - [!!! 关键BUG修复 !!!] 在计算稳态分布的代际演化时，引入了生存率。
    % ==   - 旧逻辑: 错误地假设所有个体都能活到下一期，导致计算出的
    % ==     条件分布(Dist_cond)不正确。
    % ==   - 新逻辑: 在质量从 ia 演化到 ia+1 时，乘以了正确的生存率 s_pathV(ia)。
    % ==     这使得稳态分布的计算逻辑与转轨路径的前向模拟逻辑完全一致。
    % =========================================================================
    nk = cS.nk;
    nkpps = cS.nkpps;
    nw = cS.nw_expanded;
    aD = cS.aD_new;

    mass_newborns = Z_ss_norm(1);
    if mass_newborns > 1e-9
        total_bequest_h = bq_per_capita * sum(Z_ss_norm);
        beq_per_newborn = total_bequest_h / mass_newborns;
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
        trans_mat_next_age = paramS.TrProbM_by_age{ia};

        for ik = 1:nk
            for ikpps = 1:nkpps
                for ie = 1:nw
                    mass_start = dist_cond_ia_slice(ik, ikpps, ie);
                    if mass_start < 1e-30, continue; end
                    
                    % [!!! 核心修正: 引入生存率 !!!]
                    mass_survived = mass_start * cS.s_pathV(ia);
                    if mass_survived < 1e-30, continue; end

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

                        % 使用存活下来的质量进行分配
                        mass_to_distribute = mass_survived * prob_to_enext;

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

        % function Dist = get_ss_dist_pps(polS, paramS, cS, Z_ss_norm, bq_per_capita)
        %     % =========================================================================
        %     % == 函数: get_ss_dist_pps (v1.2 - 遗赠逻辑修正版)
        %     % == 核心修正:
        %     % ==   - [!!! 关键BUG修复 !!!] 修正了对遗赠的处理逻辑。
        %     % ==   - 此函数现在接收人均遗赠 bq_per_capita (与no_pps版一致)。
        %     % ==   - 内部计算新生儿收到的遗赠时，先将人均量乘以该类型总人口，
        %     % ==     得到该类型的总遗赠，再除以新生儿数量。
        %     % =========================================================================
        %     nk = cS.nk;
        %     nkpps = cS.nkpps;
        %     nw = cS.nw_expanded;
        %     aD = cS.aD_new;
        % 
        %     mass_newborns = Z_ss_norm(1);
        %     if mass_newborns > 1e-9
        %         % [!!! 核心修正 !!!]
        %         % bq_per_capita是该类型家庭的人均遗赠
        %         % 乘以该类型总人口 sum(Z_ss_norm)，得到该类型产生的总遗赠
        %         % 再除以该类型的新生儿数量 mass_newborns，得到每个新生儿收到的遗赠
        %         total_bequest_h = bq_per_capita * sum(Z_ss_norm);
        %         beq_per_newborn = total_bequest_h / mass_newborns;
        %     else
        %         beq_per_newborn = 0;
        %     end
        % 
        %     [ik_lower, ~, w_upper] = utils.find_grid_and_weights(beq_per_newborn, cS.kGridV);
        %     w_lower = 1.0 - w_upper;
        %     ik_upper = min(nk, ik_lower + 1);
        %     ikpps_newborn = 1;
        % 
        %     dist_newborn_cond_shape = zeros(nk, nkpps, nw);
        %     prob_mass_to_distribute = reshape(paramS.leProb1V, [1, 1, cS.nw]);
        %     if nw > cS.nw, prob_mass_to_distribute(1,1, (cS.nw+1):end) = 0; end
        % 
        %     dist_newborn_cond_shape(ik_lower, ikpps_newborn, 1:cS.nw) = w_lower * prob_mass_to_distribute;
        %     dist_newborn_cond_shape(ik_upper, ikpps_newborn, 1:cS.nw) = dist_newborn_cond_shape(ik_upper, ikpps_newborn, 1:cS.nw) + w_upper * prob_mass_to_distribute;
        % 
        %     Dist_cond = zeros(nk, nkpps, nw, aD, 'double');
        %     Dist_cond(:, :, :, 1) = dist_newborn_cond_shape;
        % 
        %     for ia = 1:(aD - 1)
        %         dist_cond_ia_slice = Dist_cond(:, :, :, ia);
        %         if sum(dist_cond_ia_slice(:)) < 1e-30, continue; end
        % 
        %         dist_cond_ia_plus_1_next = zeros(nk, nkpps, nw, 'double');
        %         trans_mat_next_age = paramS.TrProbM_by_age{ia};
        % 
        %         for ik = 1:nk
        %             for ikpps = 1:nkpps
        %                 for ie = 1:nw
        %                     mass_start = dist_cond_ia_slice(ik, ikpps, ie);
        %                     if mass_start < 1e-30, continue; end
        % 
        %                     k_prime = polS(ia).k_prime(ik, ikpps, ie);
        %                     kpps_prime = polS(ia).kpps_prime(ik, ikpps, ie);
        % 
        %                     [ik_p_lower, ~, w_k_upper] = utils.find_grid_and_weights(k_prime, cS.kGridV);
        %                     w_k_lower = 1.0 - w_k_upper;
        %                     ik_p_upper = min(nk, ik_p_lower + 1);
        % 
        %                     [ikpps_p_lower, ~, w_kpps_upper] = utils.find_grid_and_weights(kpps_prime, cS.kppsGridV);
        %                     w_kpps_lower = 1.0 - w_kpps_upper;
        %                     ikpps_p_upper = min(nkpps, ikpps_p_lower + 1);
        % 
        %                     trans_probs_vec = trans_mat_next_age(ie, :);
        %                     for ie_next = 1:nw
        %                         prob_to_enext = trans_probs_vec(ie_next);
        %                         if prob_to_enext < 1e-12, continue; end
        % 
        %                         mass_to_distribute = mass_start * prob_to_enext;
        % 
        %                         mass_ll = mass_to_distribute * w_k_lower * w_kpps_lower;
        %                         mass_lu = mass_to_distribute * w_k_lower * w_kpps_upper;
        %                         mass_ul = mass_to_distribute * w_k_upper * w_kpps_lower;
        %                         mass_uu = mass_to_distribute * w_k_upper * w_kpps_upper;
        % 
        %                         dist_cond_ia_plus_1_next(ik_p_lower, ikpps_p_lower, ie_next) = dist_cond_ia_plus_1_next(ik_p_lower, ikpps_p_lower, ie_next) + mass_ll;
        %                         dist_cond_ia_plus_1_next(ik_p_lower, ikpps_p_upper, ie_next) = dist_cond_ia_plus_1_next(ik_p_lower, ikpps_p_upper, ie_next) + mass_lu;
        %                         dist_cond_ia_plus_1_next(ik_p_upper, ikpps_p_lower, ie_next) = dist_cond_ia_plus_1_next(ik_p_upper, ikpps_p_lower, ie_next) + mass_ul;
        %                         dist_cond_ia_plus_1_next(ik_p_upper, ikpps_p_upper, ie_next) = dist_cond_ia_plus_1_next(ik_p_upper, ikpps_p_upper, ie_next) + mass_uu;
        %                     end
        %                 end
        %             end
        %         end
        %         sum_next_dist = sum(dist_cond_ia_plus_1_next(:));
        %         if sum_next_dist > 1e-9
        %             Dist_cond(:, :, :, ia+1) = dist_cond_ia_plus_1_next / sum_next_dist;
        %         end
        %     end
        % 
        %     Dist = Dist_cond .* reshape(Z_ss_norm, [1, 1, 1, aD]);
        % end        
        
        
        % function Dist = get_ss_dist_no_pps(polS, paramS, cS, Z_ss_norm, bq_per_capita)
        %     % =========================================================================
        %     % == 函数: get_ss_dist_no_pps (v1.2 - 接口一致性修正版)
        %     % == 核心修正:
        %     % ==   - 将输入参数名从 bequest_per_capita 修改为 bq_per_capita，
        %     % ==     与 get_ss_dist 和 get_ss_dist_pps 保持完全一致。
        %     % ==   - 内部逻辑无需修改。
        %     % =========================================================================
        %     nk = cS.nk;
        %     nw = cS.nw_expanded;
        %     aD = cS.aD_new;
        % 
        %     mass_newborns = Z_ss_norm(1);
        % 
        %     if mass_newborns > 1e-9
        %         total_bequest_h = bq_per_capita * sum(Z_ss_norm);
        %         beq_per_newborn = total_bequest_h / mass_newborns;
        %     else
        %         beq_per_newborn = 0;
        %     end
        % 
        %     [ik_lower, ~, w_upper] = utils.find_grid_and_weights(beq_per_newborn, cS.kGridV);
        %     w_lower = 1.0 - w_upper;
        %     ik_upper = min(nk, ik_lower + 1);
        % 
        %     dist_newborn_cond_shape = zeros(nk, 1, nw);
        %     prob_mass_to_distribute = reshape(paramS.leProb1V, [1, 1, cS.nw]);
        %     if nw > cS.nw, prob_mass_to_distribute(1,1, (cS.nw+1):end) = 0; end
        %     dist_newborn_cond_shape(ik_lower, 1, 1:cS.nw) = w_lower * prob_mass_to_distribute;
        %     dist_newborn_cond_shape(ik_upper, 1, 1:cS.nw) = dist_newborn_cond_shape(ik_upper, 1, 1:cS.nw) + w_upper * prob_mass_to_distribute;
        % 
        %     Dist_cond = zeros(nk, 1, nw, aD, 'double');
        %     Dist_cond(:, :, :, 1) = dist_newborn_cond_shape;
        % 
        %     for ia = 1:(aD - 1)
        %         dist_cond_ia_slice = Dist_cond(:, 1, :, ia);
        %         if sum(dist_cond_ia_slice(:)) < 1e-30, continue; end
        % 
        %         dist_cond_ia_plus_1_next = zeros(nk, 1, nw, 'double');
        %         trans_mat_next_age = paramS.TrProbM_by_age{ia};
        % 
        %         for ik = 1:nk
        %             for ie = 1:nw
        %                 mass_start = dist_cond_ia_slice(ik, 1, ie);
        %                 if mass_start < 1e-30, continue; end
        % 
        %                 k_prime = polS(ia).k_prime(ik, 1, ie);
        %                 [ik_lower_p, ~, w_k_upper] = utils.find_grid_and_weights(k_prime, cS.kGridV);
        %                 w_k_lower = 1.0 - w_k_upper;
        %                 ik_upper_p = min(nk, ik_lower_p + 1);
        % 
        %                 trans_probs_vec = trans_mat_next_age(ie, :);
        %                 for ie_next = 1:nw
        %                     prob_to_enext = trans_probs_vec(ie_next);
        %                     if prob_to_enext < 1e-12, continue; end
        % 
        %                     mass_to_distribute = mass_start * prob_to_enext;
        %                     dist_cond_ia_plus_1_next(ik_lower_p, 1, ie_next) = dist_cond_ia_plus_1_next(ik_lower_p, 1, ie_next) + mass_to_distribute * w_k_lower;
        %                     dist_cond_ia_plus_1_next(ik_upper_p, 1, ie_next) = dist_cond_ia_plus_1_next(ik_upper_p, 1, ie_next) + mass_to_distribute * w_k_upper;
        %                 end
        %             end
        %         end
        %         sum_next_dist = sum(dist_cond_ia_plus_1_next(:));
        %         if sum_next_dist > 1e-9
        %             Dist_cond(:, :, :, ia+1) = dist_cond_ia_plus_1_next / sum_next_dist;
        %         end
        %     end
        % 
        %     Dist = Dist_cond .* reshape(Z_ss_norm, [1, 1, 1, aD]);
        % end


        function Dist_path = simulate_dist_forward(dist0, Pol_path, cS, paramS, bequest_per_newborn_hat_path, h_idx)
            % =========================================================================
            % == 函数: simulate_dist_forward (v4.2 - 异质性新生儿修正版)
            % == 核心修正:
            % ==   - [!!! 关键BUG修复 !!!] 增加了一个输入参数 `h_idx`，用于
            % ==     正确计算【特定类型h】的新生儿绝对数量。
            % ==   - 新生儿质量现在计算为: 总新生儿 * 该类型权重。
            % ==     mass_newborns_tp1_abs = Z_path_abs(1, t+1) * cS.type_weights(h_idx);
            % ==   - 这修复了在异质性模型中新生儿被重复计算的根本性错误。
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
                
                % [!!! 核心修正 !!!] 计算特定类型 h 的新生儿质量
                mass_newborns_tp1_abs = Z_path_abs(1, t+1) * cS.type_weights(h_idx);

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
