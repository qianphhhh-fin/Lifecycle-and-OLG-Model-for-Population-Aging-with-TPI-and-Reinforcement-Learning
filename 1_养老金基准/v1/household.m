classdef household
    methods (Static)

        function [cS, paramS] = parameter_values()
            % =========================================================================
            % == 函数: parameter_values (版本 v5.2 - 对标 Cocco/Gomes 参数)
            % == 目的: 设置并返回使用Epstein-Zin效用函数的模型参数。
            % == 核心修改:
            % ==   1. 将 psi (IES) 从 1.5 修改为 0.5，使其小于1，以匹配
            % ==      Gomes(2020)文献中产生高风险资产持有的偏好结构。
            % ==   2. 将 rho (RRA) 和 beta (折现因子) 调整为与 life_cycle.m
            % ==      中更接近的值，以增强可比性。
            % =========================================================================

            % 1. 生命周期与人口结构
            tb_age = 18; tr_age = 61; td_age = 100;
            cS.tb = 1;
            cS.tr = tr_age - tb_age + 1;
            cS.td = td_age - tb_age + 1;
            cS.tn = cS.td;

            % 条件生存概率
            survprob_data = [0.99975, 0.99974, 0.99973, 0.99972, 0.99971, 0.99969, 0.99968, 0.99966, 0.99963, 0.99961, ...
                0.99958, 0.99954, 0.99951, 0.99947, 0.99943, 0.99938, 0.99934, 0.99928, 0.99922, 0.99916, ...
                0.99908, 0.999, 0.99891, 0.99881, 0.9987, 0.99857, 0.99843, 0.99828, 0.99812, 0.99794, ...
                0.99776, 0.99756, 0.99735, 0.99713, 0.9969, 0.99666, 0.9964, 0.99613, 0.99583, 0.99551, ...
                0.99515, 0.99476, 0.99432, 0.99383, 0.9933, 0.9927, 0.99205, 0.99133, 0.99053, 0.98961, ...
                0.98852, 0.98718, 0.98553, 0.98346, 0.98089, 0.97772, 0.97391, 0.96943, 0.96429, 0.95854, ...
                0.95221, 0.94537, 0.93805, 0.93027, 0.92202, 0.91327, 0.90393, 0.89389, 0.88304, 0.87126, ...
                0.85846, 0.84452, 0.82935, 0.81282, 0.79485, 0.77543, 0.75458, 0.7324, 0.70893, 0.68424, ...
                0.68424, 0.68424];
            cS.s_pathV = ones(1, cS.tn);
            len_surv = length(survprob_data);
            cS.s_pathV(1:min(cS.tn-1, len_surv)) = survprob_data(1:min(cS.tn-1, len_surv));

            % 2. 偏好 (Epstein-Zin) - [!!! 核心修改区域 !!!]
            cS.beta = 0.97;         % 主观折现因子 (Cocco: 0.97)
            cS.rho = 10.0;          % 相对风险规避系数 (Cocco: 10.0)
            cS.psi = 0.5;           % 跨期替代弹性 (Cocco: 0.5), [修改前为1.5]

            % 预计算EZ参数，与Cocco代码一致
            paramS.theta = (1.0 - cS.rho) / (1.0 - 1.0/cS.psi);
            paramS.psi_1 = 1.0 - 1.0/cS.psi;
            paramS.psi_2 = 1.0 / paramS.psi_1;

            % 3. 收入过程
            aa = (-2.170042 + 2.700381); b1 = 0.16818; b2 = -0.0323371 / 10; b3 = 0.0019704 / 100;
            paramS.g_t_V = zeros(1, cS.tn);
            for t = 1:cS.tr
                age = t + tb_age - 1;
                paramS.g_t_V(t) = (aa + b1 * age + b2 * age^2 + b3 * age^3);
            end
            paramS.sigma_z = sqrt(0.112572); % 注意：Cocco代码中似乎没有持久性收入冲击，这里保留你的设定

            % 4. 养老金体系
            cS.tau_y = 0.06;
            cS.Y_bar = 0.68212; % 与Cocco的ret_fac对齐
            cS.Q_max = 0.1;
            cS.tau_q = 0.03;

            % 5. 资产与回报率
            cS.R_f = 1.015;          % Cocco: r=1.015
            paramS.mu = 0.04;           % Cocco: mu=0.04
            paramS.sigma_epsilon = 0.2; % Cocco: sigr=0.2
            cS.R_p = 1.04; % 你的设定

            % 6. 数值求解网格设定
            cS.nW = 31;
            cS.nF = 31;
            cS.nZ = 5;
            cS.nEpsilon = 5;
            cS.w_min = 0.25; cS.w_max = 200; % 调整w_max以匹配Cocco的gcash范围
            cS.f_min = 0.01; cS.f_max = 50;
            cS.wGridV = exp(linspace(log(cS.w_min), log(cS.w_max), cS.nW))';
            cS.fGridV = exp(linspace(log(cS.f_min), log(cS.f_max), cS.nF))';
            cS.fGridV(1) = 0;cS.f_min = 0;

            cS.nC = 20;
            cS.nQ = 20;
            cS.nAlpha = 21; % 使用21个点以包含0和1
            cS.alphaGridV = linspace(0, 1, cS.nAlpha);
            cS.qGridV = linspace(0, cS.Q_max, cS.nQ);


            % 7. 离散化随机过程
            [z_nodes, z_trans_mat] = household.tauchen(cS.nZ, 0, paramS.sigma_z, 0, 2.857);
            paramS.zGridV = z_nodes;
            paramS.zProbs = z_trans_mat(1, :);

            [eps_nodes, eps_prob_mat] = household.tauchen(cS.nEpsilon, 0, paramS.sigma_epsilon, 0, 2.857);
            paramS.epsilonNodes = eps_nodes;
            paramS.epsilonProbs = eps_prob_mat(1, :);
            paramS.R_shock_V = cS.R_f + paramS.mu + paramS.epsilonNodes;
        end

        % ---------VFI----------------
        function [polS, valS] = VFI_PPS(paramS, cS)
            % =========================================================================
            % == 函数: VFI_PPS (版本 v8.2 - 修正终端条件)
            % == 目的: 使用Epstein-Zin递归效用函数求解模型。
            % == 核心修正:
            % ==   1. 修正了终期 T 的价值函数和策略函数。原始代码错误地将最终消费
            % ==      设定为仅有流动性财富 W_T，忽略了当期的养老金和年金支付。
            % ==   2. 修正后的最终消费 C_T = W_T + P_bar_T + Y_bar，这正确地
            % ==      反映了个体在生命终点消费所有可支配资源的经济学假设。
            % =========================================================================

            % 1. 初始化
            valS = -inf(cS.nW, cS.nF, cS.tn);
            polS_cell = cell(cS.tn, 1);
            w_prime_frac_grid = linspace(0, 1, cS.nC);

            % 2. 终期 T
            % [!!! 核心修正 !!!]
            % 正确计算终期的 Cash-on-Hand (CoH)
            [w_mat_T, f_mat_T] = ndgrid(cS.wGridV, cS.fGridV);
            % 终期的年金支付 P_bar 仍然由 F_T (代表 F_tr) 决定
            p_bar_mat_T = (1 - cS.tau_q) * f_mat_T / (cS.td - cS.tr + 1);
            % 终期的总消费等于所有可支配资源
            final_coh = w_mat_T + p_bar_mat_T + cS.Y_bar;

            % 根据终期消费设定终期价值函数和策略
            valS(:, :, cS.tn) = final_coh * (1 - cS.beta)^paramS.psi_2;
            polS_cell{cS.tn} = struct('c', final_coh, 'q_rate', 0, 'alpha', 0);

            % 3. 逆向归纳: 退休期 (t = T-1, ..., K)
            for t = cS.tn-1 : -1 : cS.tr
                age_idx = t;
                v_prime_interp = griddedInterpolant({cS.wGridV, cS.fGridV}, valS(:, :, t+1), 'linear','linear');

                [w_mat, f_mat] = ndgrid(cS.wGridV, cS.fGridV);
                % p_bar_mat 由当前状态 f_mat (代表 F_tr) 决定
                p_bar_mat = (1 - cS.tau_q) * f_mat / (cS.td - cS.tr + 1);
                CoH = w_mat + p_bar_mat + cS.Y_bar;

                w_prime_choices = CoH .* reshape(w_prime_frac_grid, [1,1,cS.nC]);
                alpha_choices = reshape(cS.alphaGridV, [1,1,1,cS.nAlpha]);

                c_choices = CoH - w_prime_choices;
                c_choices(c_choices <= 1e-6) = 1e-6;
                util_c_term = (1 - cS.beta) * (c_choices .^ paramS.psi_1);

                R_portfolio = (1-alpha_choices).*cS.R_f + alpha_choices.*reshape(paramS.R_shock_V,[1,1,1,1,cS.nEpsilon]);

                w_prime_next = reshape(w_prime_choices, [cS.nW, cS.nF, cS.nC, 1, 1]) .* R_portfolio;
                w_prime_next = max(min(w_prime_next, cS.w_max), cS.w_min);

                % 在退休期, F_{t+1} = F_t (即 F_tr 状态保持不变)
                f_prime_next = repmat(f_mat, [1, 1, cS.nC, cS.nAlpha, cS.nEpsilon]);

                v_prime_vals = v_prime_interp(w_prime_next, f_prime_next);
                v_prime_vals(v_prime_vals <= 0) = 1e-10;

                ev_term_inner = v_prime_vals.^(1-cS.rho);
                ev_term = sum(ev_term_inner .* reshape(paramS.epsilonProbs,1,1,1,1,cS.nEpsilon), 5);

                certainty_equiv_term = ev_term .^ (1/paramS.theta);

                total_value = (reshape(util_c_term, [cS.nW, cS.nF, cS.nC, 1]) + cS.beta * cS.s_pathV(age_idx) * certainty_equiv_term) .^ paramS.psi_2;

                S_retire = size(total_value);
                total_value_reshaped = reshape(total_value, [S_retire(1), S_retire(2), S_retire(3)*S_retire(4)]);
                [val_t_mat, best_idx] = max(real(total_value_reshaped), [], 3, 'omitnan');
                [best_c_idx, best_alpha_idx] = ind2sub([cS.nC, cS.nAlpha], best_idx);

                valS(:, :, t) = val_t_mat;

                opt_alpha = cS.alphaGridV(best_alpha_idx);
                opt_w_prime_frac = reshape(w_prime_frac_grid(best_c_idx), size(best_c_idx));
                w_prime_opt = CoH .* opt_w_prime_frac;
                c_opt = CoH - w_prime_opt;
                c_opt(c_opt <= 1e-6) = 1e-6;
                q_rate_opt = zeros(size(c_opt));

                polS_cell{t} = struct('c', c_opt, 'q_rate', q_rate_opt, 'alpha', opt_alpha);
            end

            v_k_interp = griddedInterpolant({cS.wGridV, cS.fGridV}, valS(:, :, cS.tr), 'linear','linear');

            % 4. 逆向归纳: 工作期 (t = K-1, ..., 1)
            z_growth_factor = exp(paramS.zGridV);

            for t = cS.tr-1 : -1 : cS.tb
                t
                age_idx = t;
                if t == cS.tr - 1
                    v_prime_interp_handle = @(w, f) v_k_interp(w, f);
                else
                    v_prime_interp = griddedInterpolant({cS.wGridV, cS.fGridV}, valS(:,:,t+1), 'linear','linear');
                    v_prime_interp_handle = @(w,f) v_prime_interp(w,f);
                end

                hat_Y_t = 1.0;
                G_t = exp(paramS.g_t_V(t+1) - paramS.g_t_V(t));

                [w_mat, f_mat] = ndgrid(cS.wGridV, cS.fGridV);

                q_rate_choices = reshape(cS.qGridV, [1,1,cS.nQ]);
                alpha_choices = reshape(cS.alphaGridV, [1,1,1,cS.nAlpha]);

                hat_Q_t = hat_Y_t .* q_rate_choices;
                CoH = w_mat + (1 - cS.tau_y) * (hat_Y_t - hat_Q_t);

                w_prime_choices = CoH .* reshape(w_prime_frac_grid, [1,1,1,cS.nC]);
                c_choices = CoH .* reshape(1 - w_prime_frac_grid, [1,1,1,cS.nC]);
                c_choices(c_choices <= 1e-6) = 1e-6;
                util_c_term = (1 - cS.beta) * (c_choices .^ paramS.psi_1);

                f_prime_base = (f_mat + hat_Q_t) * cS.R_p;

                R_portfolio = (1-alpha_choices).*cS.R_f + alpha_choices.*reshape(paramS.R_shock_V,1,1,1,1,cS.nEpsilon);

                w_invested = reshape(w_prime_choices, [cS.nW, cS.nF, cS.nQ, 1, cS.nC, 1, 1]) .* reshape(R_portfolio, [1, 1, 1, cS.nAlpha, 1, cS.nEpsilon, 1]);

                total_growth_denom = G_t .* reshape(z_growth_factor, [1,1,1,1,1,1,cS.nZ]);

                w_prime_next_norm = w_invested ./ total_growth_denom;
                w_prime_next_norm = max(min(w_prime_next_norm, cS.w_max), cS.w_min);

                f_prime_next_norm = reshape(f_prime_base, [cS.nW, cS.nF, cS.nQ, 1, 1, 1, 1]) ./ total_growth_denom;
                f_prime_next_norm = max(0, min(f_prime_next_norm, cS.f_max));

                f_prime_next_norm_expanded = repmat(f_prime_next_norm, [1,1,1,cS.nAlpha, cS.nC, cS.nEpsilon, 1]);

                v_prime_vals = v_prime_interp_handle(w_prime_next_norm, f_prime_next_norm_expanded);
                v_prime_vals(v_prime_vals <= 0) = 1e-10;

                v_prime_denorm = v_prime_vals .* total_growth_denom;
                ev_inner_term = v_prime_denorm.^(1-cS.rho);
                ev_over_eps = sum(ev_inner_term .* reshape(paramS.epsilonProbs,1,1,1,1,1,cS.nEpsilon,1), 6);
                ev_term = sum(ev_over_eps .* reshape(paramS.zProbs,1,1,1,1,1,1,cS.nZ), 7);
                certainty_equiv_term = ev_term .^ (1.0/paramS.theta);

                total_value = (reshape(util_c_term, [cS.nW, cS.nF, cS.nQ, 1, cS.nC]) + cS.beta * cS.s_pathV(age_idx) * certainty_equiv_term) .^ paramS.psi_2;

                S_work = size(total_value);
                total_value_reshaped = reshape(total_value, [S_work(1), S_work(2), S_work(3)*S_work(4)*S_work(5)]);
                [val_t, best_idx] = max(real(total_value_reshaped), [], 3, 'omitnan');
                [best_q_idx, best_alpha_idx, best_c_idx] = ind2sub([cS.nQ, cS.nAlpha, cS.nC], best_idx);

                valS(:, :, t) = val_t;

                opt_q_rate = cS.qGridV(best_q_idx);
                opt_alpha = cS.alphaGridV(best_alpha_idx);

                hat_Y_t_mat = repmat(hat_Y_t, cS.nW, cS.nF);
                Q_t_opt = hat_Y_t_mat .* opt_q_rate;
                CoH_opt = w_mat + (1 - cS.tau_y) * (hat_Y_t_mat - Q_t_opt);

                opt_w_prime_frac = reshape(w_prime_frac_grid(best_c_idx), size(best_c_idx));
                w_prime_opt = CoH_opt .* opt_w_prime_frac;
                c_opt = CoH_opt - w_prime_opt;

                polS_cell{t} = struct('c', c_opt, 'q_rate', opt_q_rate, 'alpha', opt_alpha);
            end
            polS = [polS_cell{:}];
        end


        function [polS, valS] = VFI_PPS_for(paramS, cS)
            % =========================================================================
            % == 函数: VFI_PPS_for (版本 v8.2 - 修正终端条件)
            % == 目的: 使用显式For循环实现EZ效用模型的求解，用于验证。
            % == 核心修正:
            % ==   1. 修正了终期 T 的价值函数和策略函数，以正确反映最终消费。
            % =========================================================================

            % 1. 初始化
            valS = -inf(cS.nW, cS.nF, cS.tn);
            polS_cell = cell(cS.tn, 1);
            w_prime_frac_grid = linspace(0, 1, cS.nC);

            % 2. 终期 T
            % [!!! 核心修正 !!!]
            [w_mat_T, f_mat_T] = ndgrid(cS.wGridV, cS.fGridV);
            p_bar_mat_T = (1 - cS.tau_q) * f_mat_T / (cS.td - cS.tr + 1);
            final_coh = w_mat_T + p_bar_mat_T + cS.Y_bar;

            valS(:, :, cS.tn) = final_coh * (1 - cS.beta)^paramS.psi_2;
            polS_cell{cS.tn} = struct('c', final_coh, 'q_rate', 0, 'alpha', 0);

            % 3. 逆向归纳: 退休期 (t = T-1, ..., K)
            for t = cS.tn-1 : -1 : cS.tr
                age_idx = t;
                v_prime_interp = griddedInterpolant({cS.wGridV, cS.fGridV}, valS(:, :, t+1), 'linear','linear');

                c_pol = zeros(cS.nW, cS.nF);
                alpha_pol = zeros(cS.nW, cS.nF);
                val_t_mat = zeros(cS.nW, cS.nF);

                for iff = 1:cS.nF
                    f_t = cS.fGridV(iff);
                    P_bar = (1 - cS.tau_q) * f_t / (cS.td - cS.tr + 1);

                    for iw = 1:cS.nW
                        w_t = cS.wGridV(iw);
                        max_v = -inf;
                        best_c = nan; best_alpha = nan;

                        CoH = w_t + cS.Y_bar + P_bar;
                        for iAlpha = 1:cS.nAlpha
                            alpha_t = cS.alphaGridV(iAlpha);
                            for iFrac = 1:cS.nC
                                w_prime_frac = w_prime_frac_grid(iFrac);
                                w_prime = CoH * w_prime_frac;
                                c_t = CoH - w_prime;
                                if c_t <= 1e-6, c_t = 1e-6; end
                                util_c_term = (1 - cS.beta)*(c_t^paramS.psi_1);

                                ev_term_inner = 0;
                                for iEps = 1:cS.nEpsilon
                                    R_port = (1 - alpha_t)*cS.R_f + alpha_t*paramS.R_shock_V(iEps);
                                    w_next = w_prime * R_port;
                                    w_next = max(min(w_next, cS.w_max), cS.w_min);

                                    % F_t 状态保持不变传递到下一期
                                    v_val = v_prime_interp(w_next, f_t);
                                    if v_val <= 0, v_val = 1e-10; end

                                    ev_term_inner = ev_term_inner + paramS.epsilonProbs(iEps) * v_val^(1-cS.rho);
                                end
                                certainty_equiv_term = ev_term_inner^(1/paramS.theta);

                                current_v = (util_c_term + cS.beta * cS.s_pathV(age_idx) * certainty_equiv_term)^paramS.psi_2;

                                if real(current_v) > max_v
                                    max_v = real(current_v);
                                    best_c = c_t;
                                    best_alpha = alpha_t;
                                end
                            end
                        end
                        val_t_mat(iw, iff) = max_v;
                        c_pol(iw, iff) = best_c;
                        alpha_pol(iw, iff) = best_alpha;
                        coh(iw, iff) = CoH;
                    end
                end
                valS(:, :, t) = val_t_mat;
                q_pol = zeros(size(c_pol));
                polS_cell{t} = struct('c', c_pol, 'q_rate', q_pol, 'alpha', alpha_pol);
            end

            v_k_interp = griddedInterpolant({cS.wGridV, cS.fGridV}, valS(:, :, cS.tr), 'linear','linear');

            % 4. 逆向归纳: 工作期 (t = K-1, ..., 1)
            for t = cS.tr-1 : -1 : cS.tb
                age_idx = t;
                if t == cS.tr - 1
                    v_prime_interp_handle = @(w,f) v_k_interp(w,f);
                else
                    v_prime_interp = griddedInterpolant({cS.wGridV, cS.fGridV}, valS(:,:,t+1), 'linear','linear');
                    v_prime_interp_handle = @(w,f) v_prime_interp(w,f);
                end

                c_pol = zeros(cS.nW, cS.nF);
                q_pol = zeros(cS.nW, cS.nF);
                alpha_pol = zeros(cS.nW, cS.nF);

                hat_Y_t = 1.0;
                G_t = exp(paramS.g_t_V(t+1) - paramS.g_t_V(t));

                for iff = 1:cS.nF
                    f_t = cS.fGridV(iff);
                    for iw = 1:cS.nW
                        w_t = cS.wGridV(iw);
                        max_v = -inf;
                        best_c = 0; best_q = 0; best_alpha = 0;

                        for iQ = 1:cS.nQ
                            q_rate = cS.qGridV(iQ);
                            Q_t = hat_Y_t * q_rate;
                            if Q_t > hat_Y_t, continue; end

                            CoH = w_t + (1 - cS.tau_y) * (hat_Y_t - Q_t);
                            for iAlpha = 1:cS.nAlpha
                                alpha_t = cS.alphaGridV(iAlpha);
                                for iFrac = 1:cS.nC
                                    w_prime_frac = w_prime_frac_grid(iFrac);
                                    w_prime = CoH * w_prime_frac;
                                    c_t = CoH - w_prime;
                                    if c_t <= 1e-6, c_t = 1e-6; end
                                    util_c_term = (1 - cS.beta)*(c_t^paramS.psi_1);

                                    f_invested = (f_t + Q_t) * cS.R_p;

                                    ev_term_inner = 0;
                                    for iZ = 1:cS.nZ
                                        z_prob = paramS.zProbs(iZ);
                                        total_growth_factor = G_t * exp(paramS.zGridV(iZ));

                                        f_next_norm = f_invested / total_growth_factor;
                                        f_next_norm = max(0, min(f_next_norm, cS.f_max));

                                        ev_cond_on_z = 0;
                                        for iEps = 1:cS.nEpsilon
                                            R_port = (1-alpha_t)*cS.R_f + alpha_t*paramS.R_shock_V(iEps);
                                            w_invested = w_prime * R_port;
                                            w_next_norm = w_invested / total_growth_factor;
                                            w_next_norm = max(min(w_next_norm, cS.w_max), cS.w_min);

                                            v_val = v_prime_interp_handle(w_next_norm, f_next_norm);
                                            if v_val <= 0, v_val = 1e-10; end

                                            ev_cond_on_z = ev_cond_on_z + paramS.epsilonProbs(iEps) * (v_val * total_growth_factor)^(1 - cS.rho);
                                        end
                                        ev_term_inner = ev_term_inner + z_prob * ev_cond_on_z;
                                    end

                                    certainty_equiv_term = ev_term_inner^(1/paramS.theta);
                                    current_v = (util_c_term + cS.beta * cS.s_pathV(age_idx) * certainty_equiv_term)^paramS.psi_2;

                                    if real(current_v) > max_v
                                        max_v = real(current_v);
                                        best_c = c_t;
                                        best_q = q_rate;
                                        best_alpha = alpha_t;
                                    end
                                end
                            end
                        end
                        valS(iw, iff, t) = max_v;
                        c_pol(iw, iff) = best_c;
                        q_pol(iw, iff) = best_q;
                        alpha_pol(iw, iff) = best_alpha;
                    end
                end
                polS_cell{t} = struct('c', c_pol, 'q_rate', q_pol, 'alpha', alpha_pol);
            end
            polS = [polS_cell{:}];
        end


        % -----------Cocco的模型------------
        function [polS, valS] = VFI_cocco_for(paramS, cS)
            % =========================================================================
            % == 函数: VFI_cocco_for
            % == 目的: 复刻 Cocco/Gomes (2020) 无养老金账户模型。
            % == 核心修改:
            % ==   1. 移除状态变量 F 和控制变量 q_rate。
            % ==   2. 状态空间简化为 W (流动性财富)。
            % ==   3. 价值函数 valS 变为二维 [nW, tn]。
            % ==   4. 插值对象 v_prime_interp 变为一维。
            % =========================================================================

            % 1. 初始化 (无 F 维度)
            valS = -inf(cS.nW, cS.tn);
            polS_cell = cell(cS.tn, 1);
            w_prime_frac_grid = linspace(0, 1, cS.nC);

            % 2. 终期 T
            % 终期消费等于所有可用的流动性财富 W 加上基础养老金 Y_bar
            final_coh = cS.wGridV + cS.Y_bar;
            valS(:, cS.tn) = final_coh * (1 - cS.beta)^paramS.psi_2;
            polS_cell{cS.tn} = struct('c', final_coh, 'alpha', 0);

            % 3. 逆向归纳: 退休期 (t = T-1, ..., K)
            for t = cS.tn-1 : -1 : cS.tr
                t
                age_idx = t;
                % 价值函数插值对象变为一维
                v_prime_interp = griddedInterpolant(cS.wGridV, valS(:, t+1), 'linear','linear');

                c_pol = zeros(cS.nW, 1);
                alpha_pol = zeros(cS.nW, 1);
                val_t_vec = zeros(cS.nW, 1);

                % 循环只针对财富状态 w
                for iw = 1:cS.nW
                    w_t = cS.wGridV(iw);
                    max_v = -inf;
                    best_c = nan; best_alpha = nan;

                    % CoH 只包括流动性财富和基础养老金
                    CoH = w_t + cS.Y_bar;

                    for iAlpha = 1:cS.nAlpha
                        alpha_t = cS.alphaGridV(iAlpha);
                        for iFrac = 1:cS.nC
                            w_prime_frac = w_prime_frac_grid(iFrac);
                            w_prime = CoH * w_prime_frac;
                            c_t = CoH - w_prime;
                            if c_t <= 1e-6, c_t = 1e-6; end
                            util_c_term = (1 - cS.beta)*(c_t^paramS.psi_1);

                            ev_term_inner = 0;
                            for iEps = 1:cS.nEpsilon
                                R_port = (1 - alpha_t)*cS.R_f + alpha_t*paramS.R_shock_V(iEps);
                                w_next = w_prime * R_port;
                                w_next = max(min(w_next, cS.w_max), cS.w_min);

                                % 插值调用变为一维
                                v_val = v_prime_interp(w_next);
                                if v_val <= 0, v_val = 1e-10; end

                                ev_term_inner = ev_term_inner + paramS.epsilonProbs(iEps) * v_val^(1-cS.rho);
                            end
                            certainty_equiv_term = ev_term_inner^(1/paramS.theta);

                            current_v = (util_c_term + cS.beta * cS.s_pathV(age_idx) * certainty_equiv_term)^paramS.psi_2;

                            if real(current_v) > max_v
                                max_v = real(current_v);
                                best_c = c_t;
                                best_alpha = alpha_t;
                            end
                        end
                    end
                    val_t_vec(iw) = max_v;
                    c_pol(iw) = best_c;
                    if best_c==CoH
                        best_alpha=nan;
                    end
                    alpha_pol(iw) = best_alpha;
                end
                valS(:, t) = val_t_vec;
                % 策略函数中不再有 q_rate
                polS_cell{t} = struct('c', c_pol, 'alpha', alpha_pol);
            end

            % 创建退休期第一年的价值函数插值器
            v_k_interp = griddedInterpolant(cS.wGridV, valS(:, cS.tr), 'linear','linear');

            % 4. 逆向归纳: 工作期 (t = K-1, ..., 1)
            for t = cS.tr-1 : -1 : cS.tb
                t
                age_idx = t;
                if t == cS.tr - 1
                    v_prime_interp_handle = @(w) v_k_interp(w);
                else
                    v_prime_interp = griddedInterpolant(cS.wGridV, valS(:,t+1), 'linear','linear');
                    v_prime_interp_handle = @(w) v_prime_interp(w);
                end

                c_pol = zeros(cS.nW, 1);
                alpha_pol = zeros(cS.nW, 1);

                hat_Y_t = 1.0;
                G_t = exp(paramS.g_t_V(t+1) - paramS.g_t_V(t));

                for iw = 1:cS.nW
                    w_t = cS.wGridV(iw);
                    max_v = -inf;
                    best_c = 0; best_alpha = 0;

                    % CoH 计算中移除养老金缴费 Q
                    CoH = w_t + (1 - cS.tau_y) * hat_Y_t;

                    for iAlpha = 1:cS.nAlpha
                        alpha_t = cS.alphaGridV(iAlpha);
                        for iFrac = 1:cS.nC
                            w_prime_frac = w_prime_frac_grid(iFrac);
                            w_prime = CoH * w_prime_frac;
                            c_t = CoH - w_prime;
                            if c_t <= 1e-6, c_t = 1e-6; end
                            util_c_term = (1 - cS.beta)*(c_t^paramS.psi_1);

                            ev_term_inner = 0;
                            for iZ = 1:cS.nZ
                                z_prob = paramS.zProbs(iZ);
                                total_growth_factor = G_t * exp(paramS.zGridV(iZ));

                                ev_cond_on_z = 0;
                                for iEps = 1:cS.nEpsilon
                                    R_port = (1-alpha_t)*cS.R_f + alpha_t*paramS.R_shock_V(iEps);
                                    w_invested = w_prime * R_port;
                                    w_next_norm = w_invested / total_growth_factor;
                                    w_next_norm = max(min(w_next_norm, cS.w_max), cS.w_min);

                                    % 插值调用变为一维，且不再需要 f_next_norm
                                    v_val = v_prime_interp_handle(w_next_norm);
                                    if v_val <= 0, v_val = 1e-10; end

                                    ev_cond_on_z = ev_cond_on_z + paramS.epsilonProbs(iEps) * (v_val * total_growth_factor)^(1 - cS.rho);
                                end
                                ev_term_inner = ev_term_inner + z_prob * ev_cond_on_z;
                            end

                            certainty_equiv_term = ev_term_inner^(1/paramS.theta);
                            current_v = (util_c_term + cS.beta * cS.s_pathV(age_idx) * certainty_equiv_term)^paramS.psi_2;

                            if real(current_v) > max_v
                                max_v = real(current_v);
                                best_c = c_t;
                                best_alpha = alpha_t;
                            end
                        end
                    end
                    valS(iw, t) = max_v;
                    c_pol(iw) = best_c;
                    alpha_pol(iw) = best_alpha;
                end
                polS_cell{t} = struct('c', c_pol, 'alpha', alpha_pol);
            end
            polS = [polS_cell{:}];
        end


        function [polS, valS] = VFI_cocco_raw(paramS, cS)
            % =========================================================================
            % == 函数: VFI_cocco_raw (版本 v1.2 - 全局优化验证版)
            % == 目的: 使用与 life_cycle.m 相同的经济环境（CoH状态变量），但
            % ==       采用数值上更稳健的全局消费搜索，以验证并寻找全局最优解。
            % == 核心修改:
            % ==   1. 彻底移除自适应局部消费搜索机制 (min_c, max_c, mpc)。
            % ==   2. 引入了基于储蓄率的全局消费网格搜索，确保覆盖所有可能的
            % ==      消费/储蓄分割，与 VFI_cocco_for 的方法论一致。
            % =========================================================================

            % 1. 初始化
            valS = -inf(cS.nW, cS.tn);
            polS_cell = cell(cS.tn, 1);

            % 状态变量网格直接定义为 Cash-on-Hand 网格
            cashGridV = cS.wGridV;

            % [!!! 核心方法论改变 !!!]
            % 定义一个全局的储蓄率网格，覆盖从0%到100%的可能性
            savings_frac_grid = linspace(0, 1, cS.nC);

            % 2. life_cycle.m 的经济环境和随机过程设定
            tb_age = 18;

            % Gauss-Hermite 节点和权重
            gh_grid = [-2.85697; -1.355626; 0.0; 1.355626; 2.85697];
            gh_weig = [0.0112574; 0.2220759; 0.5333333; 0.2220759; 0.0112574];
            n_gh = 5;

            % 收入冲击标准差
            smav = 0.1; % 持续性冲击
            smay = 0.1; % 暂时性冲击

            % 预计算冲击 realization
            persistent_shocks = exp(gh_grid * smav);
            transitory_shocks = exp(gh_grid * smay);
            return_shocks = cS.R_f + paramS.mu + gh_grid * paramS.sigma_epsilon;

            % 确定性收入年龄剖面
            aa = (-2.170042 + 2.700381); b1 = 0.16818; b2 = -0.0323371 / 10; b3 = 0.0019704 / 100;
            f_y = zeros(cS.tr, 1);
            for t = 1:cS.tr
                age = t + tb_age - 1;
                f_y(t) = exp(aa + b1*age + b2*age^2 + b3*age^3);
            end

            gy = zeros(cS.tr-1, 1);
            for t = 1:cS.tr-1
                gy(t) = f_y(t+1) / f_y(t);
            end

            % 3. 终期 T (采用逻辑一致的设定)
            final_coh_at_T = cashGridV + cS.Y_bar;
            valS(:, cS.tn) = final_coh_at_T * (1 - cS.beta)^paramS.psi_2;
            polS_cell{cS.tn} = struct('c', final_coh_at_T, 'alpha', 0);

            % 4. 逆向归纳: 退休期 (t = T-1, ..., K)
            for t = cS.tn-1 : -1 : cS.tr
                age_idx = t;
                % 下一期的 CoH_{t+1} = sav * R + Y_bar
                % 因此插值器的横坐标是定义在下一期CoH上的，即 cashGridV + Y_bar
                v_prime_interp = @(x) interp1(cashGridV + cS.Y_bar, valS(:, t+1), x, 'spline', 'extrap');

                for iCash = 1:cS.nW
                    coh_t = cashGridV(iCash);
                    max_v = -inf;
                    best_c = nan; best_alpha = nan;

                    for iAlpha = 1:cS.nAlpha
                        alpha_t = cS.alphaGridV(iAlpha);
                        % [!!! 核心修改: 不再使用动态局部网格，改用全局搜索 !!!]
                        for iFrac = 1:cS.nC
                            sav = coh_t * savings_frac_grid(iFrac);
                            c_t = coh_t - sav;

                            if c_t <= 1e-6, c_t = 1e-6; end

                            util_c_term = (1 - cS.beta)*(c_t^paramS.psi_1);

                            ev_term_inner = 0;
                            for iEps = 1:n_gh
                                R_port = (1-alpha_t)*cS.R_f + alpha_t*return_shocks(iEps);
                                coh_next = sav * R_port + cS.Y_bar;

                                v_val = v_prime_interp(coh_next);
                                if v_val <= 0, v_val = 1e-10; end
                                ev_term_inner = ev_term_inner + gh_weig(iEps) * v_val^(1-cS.rho);
                            end
                            certainty_equiv_term = ev_term_inner^(1/paramS.theta);
                            current_v = (util_c_term + cS.beta * cS.s_pathV(age_idx) * certainty_equiv_term)^paramS.psi_2;

                            if real(current_v) > max_v
                                max_v = real(current_v);
                                best_c = c_t;
                                best_alpha = alpha_t;
                            end
                        end
                    end
                    valS(iCash, t) = max_v;
                    polS_cell{t}.c(iCash,1) = best_c;
                    polS_cell{t}.alpha(iCash,1) = best_alpha;
                end
            end

            % 5. 逆向归纳: 工作期 (t = K-1, ..., 1)
            for t = cS.tr-1 : -1 : cS.tb
                age_idx = t;
                v_prime_interp = @(x) interp1(cashGridV, valS(:, t+1), x, 'spline', 'extrap');

                for iCash = 1:cS.nW
                    coh_t = cashGridV(iCash);
                    max_v = -inf;
                    best_c = nan; best_alpha = nan;

                    for iAlpha = 1:cS.nAlpha
                        alpha_t = cS.alphaGridV(iAlpha);
                        % [!!! 核心修改: 改用全局搜索 !!!]
                        for iFrac = 1:cS.nC
                            sav = coh_t * savings_frac_grid(iFrac);
                            c_t = coh_t - sav;
                            if c_t <= 1e-6, c_t = 1e-6; end

                            util_c_term = (1 - cS.beta)*(c_t^paramS.psi_1);

                            ev_term_inner = 0;
                            for iPersist = 1:n_gh % 持续性冲击
                                growth_factor = gy(t) * persistent_shocks(iPersist);

                                ev_cond_on_persist = 0;
                                for iTrans = 1:n_gh % 暂时性冲击
                                    y_trans_shock = transitory_shocks(iTrans);

                                    ev_cond_on_trans = 0;
                                    for iEps = 1:n_gh % 投资回报冲击
                                        R_port = (1-alpha_t)*cS.R_f + alpha_t*return_shocks(iEps);

                                        coh_next = (sav * R_port) / growth_factor + y_trans_shock;
                                        coh_next = max(min(coh_next, cashGridV(end)), cashGridV(1));

                                        v_val = v_prime_interp(coh_next);
                                        if v_val <= 0, v_val = 1e-10; end

                                        ev_cond_on_trans = ev_cond_on_trans + gh_weig(iEps) * (v_val * growth_factor)^(1-cS.rho);
                                    end
                                    ev_cond_on_persist = ev_cond_on_persist + gh_weig(iTrans) * ev_cond_on_trans;
                                end
                                ev_term_inner = ev_term_inner + gh_weig(iPersist) * ev_cond_on_persist;
                            end

                            certainty_equiv_term = ev_term_inner^(1/paramS.theta);
                            current_v = (util_c_term + cS.beta * cS.s_pathV(age_idx) * certainty_equiv_term)^paramS.psi_2;

                            if real(current_v) > max_v
                                max_v = real(current_v);
                                best_c = c_t;
                                best_alpha = alpha_t;
                            end
                        end
                    end
                    valS(iCash, t) = max_v;
                    polS_cell{t}.c(iCash,1) = best_c;
                    polS_cell{t}.alpha(iCash,1) = best_alpha;
                end
            end
            polS = [polS_cell{:}];
        end


        function test_cocco_comparison()
            % =========================================================================
            % == 函数: test_cocco_comparison
            % == 目的: 对比 VFI_cocco_for (W-based) 和 VFI_cocco_raw (CoH-based)
            % ==       的结果，验证两者在经济上的一致性。
            % =========================================================================
            fprintf('开始测试 Cocco/Gomes 模型复刻版本的等价性...\n');

            % 1. 获取参数并使用小网格
            [cS, paramS] = household.parameter_values();
            cS.nW = 51; cS.nC =21; cS.nAlpha = 51;
            % 为了测试，暂时重设wGridV
            cS.w_min = 0.25; cS.w_max = 200;
            cS.wGridV = exp(linspace(log(cS.w_min), log(cS.w_max), cS.nW))';
            cS.alphaGridV = linspace(0, 1, cS.nAlpha);

            fprintf('测试参数: nW=%d, nC=%d, nAlpha=%d\n', cS.nW, cS.nC, cS.nAlpha);
            fprintf('生命周期: tn=%d, tr=%d\n\n', cS.tn, cS.tr);




            % 3. 运行 Raw (CoH-based) 版本
            fprintf('--> 正在运行 CoH-based raw 版本 (VFI_cocco_raw)...\n');
            tic;
            [polS_raw, valS_raw] = household.VFI_cocco_raw(paramS, cS);
            time_raw = toc;
            fprintf('    CoH-based raw 版本完成，耗时: %.4f 秒\n\n', time_raw);


            % 2. 运行 For-Loop (W-based) 版本
            fprintf('--> 正在运行 W-based 版本 (VFI_cocco_for)...\n');
            tic;
            [polS_for, valS_for] = household.VFI_cocco_for(paramS, cS);
            time_for = toc;
            fprintf('    W-based 版本完成，耗时: %.4f 秒\n\n', time_for);
            % 4. 逐期比较结果
            fprintf('--> 正在逐期比较两个版本的结果...\n');

            num_periods = cS.tn - cS.tb;
            error_table = table('Size', [num_periods+1, 4], ...
                'VariableTypes', {'double', 'double', 'double', 'double'}, ...
                'VariableNames', {'t', 'Max_Val_Diff', 'Max_C_Diff', 'Max_Alpha_Diff'});

            row_idx = 1;
            w_grid = cS.wGridV;
            cash_grid_raw = cS.wGridV; % 在raw版本中，wGridV被用作cash网格

            for t = cS.tn : -1 : cS.tb
                % 获取当前 t 的策略和价值
                val_for_t = valS_for(:, t);
                c_for_t = polS_for(t).c;
                alpha_for_t = polS_for(t).alpha;

                val_raw_t = valS_raw(:, t);
                c_raw_t = polS_raw(t).c;
                alpha_raw_t = polS_raw(t).alpha;

                % 关键：将 raw (CoH-based) 的结果插值到 for (W-based) 的等效点上
                if t >= cS.tr % 退休期
                    % CoH = W + Y_bar
                    coh_equiv_grid = w_grid + cS.Y_bar;
                    coh_for_t = w_grid + cS.Y_bar;
                else % 工作期
                    % CoH = W + (1 - tau_y) * hat_Y, (hat_Y = 1)
                    coh_equiv_grid = w_grid + (1 - cS.tau_y);
                    coh_for_t = w_grid + (1-cS.tau_y);
                end

                % 使用线性插值进行转换
                val_raw_interp = interp1(cash_grid_raw, val_raw_t, coh_equiv_grid, 'linear', 'extrap');
                c_raw_interp = interp1(cash_grid_raw, c_raw_t, coh_equiv_grid, 'linear', 'extrap');
                alpha_raw_interp = interp1(cash_grid_raw, alpha_raw_t, coh_equiv_grid, 'linear', 'extrap');

                % 消费需要从绝对水平转换为相对CoH的比例才能比较
                % C_for 是从 W-based CoH 中消费的量
                % C_raw_interp 是从 W-based CoH 中消费的量 (经过插值)
                % 两者可以直接比较

                % 计算差异
                max_val_diff = max(abs(val_for_t - val_raw_interp));
                max_c_diff = max(abs(c_for_t - c_raw_interp));
                max_alpha_diff = max(abs(alpha_for_t - alpha_raw_interp));

                error_table(row_idx, :) = {t, max_val_diff, max_c_diff, max_alpha_diff};
                row_idx = row_idx + 1;
            end

            % 5. 显示结果
            error_table = sortrows(error_table, 't');
            fprintf('\n--> 逐期最大误差对比表 (W-based vs CoH-based):\n');
            disp(error_table);

            output_dir = 'debug';
            if ~exist(output_dir, 'dir'), mkdir(output_dir); end
            output_filename = fullfile(output_dir, 'cocco_comparison_errors.txt');
            try
                writetable(error_table, output_filename, 'Delimiter', '\t');
                fprintf('\n    误差对比表已成功保存至文件: %s\n\n', output_filename);
            catch ME
                fprintf('\n    保存文件时出错: %s\n', ME.message);
            end

            fprintf('--> 性能总结:\n');
            if time_for > 1e-9
                fprintf('    CoH-based raw 版本比 W-based for 版本慢 %.2f 倍。\n', time_raw / time_for);
            end
            fprintf('测试结束。\n');
        end

        % ----------------数值模拟----------------
        function resultsS = simulation(polS, paramS, cS)
            % =========================================================================
            % == 函数: simulation (版本 v2.1 - 增加更多结果变量)
            % == 目的: 基于 VFI_PPS 的求解结果，进行生命周期模拟。
            % == 核心修改:
            % ==   1. 增加存储 simY_absM, simCoH_absM, simPension_absM 矩阵。
            % ==   2. 在最终结果 resultsS 中，增加这些变量的均值以及财富收入比。
            % =========================================================================

            fprintf('开始进行生命周期模拟...\n');

            % 1. 初始化
            nsim = 2000; % 模拟个体数量
            rng(2024);   % 固定随机种子

            % 初始化存储矩阵
            simW_normM = zeros(cS.tn, nsim); % 归一化财富
            simF_normM = zeros(cS.tn, nsim); % 归一化养老金
            Y_levelM = ones(cS.tn, nsim);   % 收入水平（用于反归一化）

            simC_absM = zeros(cS.tn, nsim);      % 绝对消费
            simW_absM = zeros(cS.tn, nsim);      % 绝对财富
            simF_absM = zeros(cS.tn, nsim);      % 绝对养老金
            simAlpha_M = zeros(cS.tn, nsim);     % 风险资产配置比例
            simQ_rate_M = zeros(cS.tn, nsim);    % 养老金缴费率

            % [!!! 新增: 初始化更多结果矩阵 !!!]
            simY_absM = zeros(cS.tn, nsim);      % 绝对收入
            simCoH_absM = zeros(cS.tn, nsim);    % 绝对可支配资源
            simPension_absM = zeros(cS.tn, nsim);% 绝对养老金给付

            % 设置初始状态
            simW_normM(1, :) = cS.w_min; % 从最低财富开始
            simF_normM(1, :) = 0;       % 从零养老金开始
            simY_absM(1,:) = Y_levelM(1,:); % 初始收入水平

            % 2. 生成生命周期内的随机冲击
            rand_z = rand(cS.tn-1, nsim);
            rand_eps = rand(cS.tn-1, nsim);

            z_prob_cum = cumsum(paramS.zProbs);
            eps_prob_cum = cumsum(paramS.epsilonProbs);

            z_idx_M = zeros(cS.tn-1, nsim);
            eps_idx_M = zeros(cS.tn-1, nsim);
            for t_sim = 1:(cS.tn-1)
                for i_sim = 1:nsim
                    z_idx_M(t_sim, i_sim) = find(rand_z(t_sim, i_sim) <= z_prob_cum, 1, 'first');
                    eps_idx_M(t_sim, i_sim) = find(rand_eps(t_sim, i_sim) <= eps_prob_cum, 1, 'first');
                end
            end

            % 3. 逐期迭代模拟
            fprintf('--> 正在逐期模拟个体决策...\n');
            [f_grid, w_grid] = meshgrid(cS.fGridV, cS.wGridV);

            p_bar_norm_V = zeros(1, nsim);

            for t = 1:cS.tn-1
                w_t_normV = simW_normM(t, :);
                f_t_normV = simF_normM(t, :);

                w_t_normV = max(min(w_t_normV, cS.w_max), cS.w_min);
                f_t_normV = max(min(f_t_normV, cS.f_max), cS.f_min);

                if t >= cS.tr
                    f_for_interp_V = simF_normM(cS.tr, :);
                else
                    f_for_interp_V = f_t_normV;
                end
                c_t_normV = interp2(f_grid, w_grid, polS(t).c, f_for_interp_V, w_t_normV, 'linear');
                alpha_t_V = interp2(f_grid, w_grid, polS(t).alpha, f_for_interp_V, w_t_normV, 'linear');

                if t < cS.tr
                    q_rate_t_V = interp2(f_grid, w_grid, polS(t).q_rate, f_t_normV, w_t_normV, 'linear');
                else
                    q_rate_t_V = zeros(1, nsim);
                end
                q_rate_t_V =  0;

                Y_level_t_V = Y_levelM(t, :);

                if t < cS.tr % 工作期
                    hat_Y_t = 1.0;
                    hat_Q_t = hat_Y_t .* q_rate_t_V;
                    CoH_norm = w_t_normV + (1 - cS.tau_y) * (hat_Y_t - hat_Q_t);

                    % [!!! 新增: 记录绝对收入和CoH !!!]
                    simY_absM(t,:) = hat_Y_t .* Y_level_t_V;
                    simCoH_absM(t,:) = CoH_norm .* Y_level_t_V;

                    c_t_normV = min(c_t_normV, CoH_norm);
                    c_t_normV(c_t_normV < 0) = 0;

                    w_prime_norm = CoH_norm - c_t_normV;

                    R_shock_V = paramS.R_shock_V(eps_idx_M(t, :))';
                    R_portfolio_V = (1 - alpha_t_V) .* cS.R_f + alpha_t_V .* R_shock_V;

                    w_invested_abs = w_prime_norm .* Y_level_t_V;
                    w_next_abs = w_invested_abs .* R_portfolio_V;

                    f_invested_abs = (f_t_normV + hat_Q_t) .* Y_level_t_V;
                    f_next_abs = f_invested_abs * cS.R_p;

                    G_t = exp(paramS.g_t_V(t+1) - paramS.g_t_V(t));
                    z_shocks_V = paramS.zGridV(z_idx_M(t, :))';
                    total_growth_factor_V = G_t * exp(z_shocks_V);

                    Y_levelM(t+1, :) = Y_level_t_V .* total_growth_factor_V;
                    simY_absM(t+1,:) = Y_levelM(t+1,:); % 记录下一期的收入水平

                    simW_normM(t+1, :) = w_next_abs ./ Y_levelM(t+1, :);
                    simF_normM(t+1, :) = f_next_abs ./ Y_levelM(t+1, :);

                else % 退休期
                    if t == cS.tr
                        p_bar_norm_V = (1 - cS.tau_q) * f_t_normV / (cS.td - cS.tr + 1);
                    end

                    CoH_norm = w_t_normV + p_bar_norm_V + cS.Y_bar;

                    % [!!! 新增: 记录绝对收入、养老金和CoH !!!]
                    % 退休期的 "收入" 是指基础养老金
                    simY_absM(t, :) = cS.Y_bar .* Y_level_t_V;
                    simPension_absM(t, :) = p_bar_norm_V .* Y_level_t_V;
                    simCoH_absM(t, :) = CoH_norm .* Y_level_t_V;

                    c_t_normV = min(c_t_normV, CoH_norm);
                    c_t_normV(c_t_normV < 0) = 0;

                    w_prime_norm = CoH_norm - c_t_normV;

                    R_shock_V = paramS.R_shock_V(eps_idx_M(t, :))';
                    R_portfolio_V = (1 - alpha_t_V) .* cS.R_f + alpha_t_V .* R_shock_V;

                    w_invested_abs = w_prime_norm .* Y_level_t_V;
                    w_next_abs = w_invested_abs .* R_portfolio_V;

                    f_next_abs = zeros(1, nsim);

                    Y_levelM(t+1, :) = Y_levelM(t, :);
                    simY_absM(t+1, :) = cS.Y_bar .* Y_levelM(t+1,:); % 记录下一期的收入水平

                    simW_normM(t+1, :) = w_next_abs ./ Y_levelM(t+1, :);
                    simF_normM(t+1, :) = f_next_abs ./ Y_levelM(t+1, :);
                end

                simC_absM(t, :) = c_t_normV .* Y_level_t_V;
                simW_absM(t, :) = w_t_normV .* Y_level_t_V;
                simF_absM(t, :) = f_t_normV .* Y_level_t_V;
                simAlpha_M(t, :) = alpha_t_V;
                simQ_rate_M(t, :) = q_rate_t_V;
            end

            % 处理最后一期
            w_tn_abs = simW_normM(cS.tn, :) .* Y_levelM(cS.tn, :);
            simW_absM(cS.tn, :) = w_tn_abs;
            simF_absM(cS.tn, :) = simF_normM(cS.tn, :) .* Y_levelM(cS.tn, :);
            simC_absM(cS.tn, :) = w_tn_abs;
            simCoH_absM(cS.tn, :) = w_tn_abs;
            simPension_absM(cS.tn, :) = simPension_absM(cS.tn-1, :); % 假设最后一期pension和之前一样
            simAlpha_M(cS.tn, :) = 0;
            simQ_rate_M(cS.tn, :) = 0;

            % 4. 汇总并打包结果
            fprintf('--> 正在汇总模拟结果...\n');
            resultsS.mean_C = mean(simC_absM, 2, 'omitnan');
            resultsS.mean_W = mean(simW_absM, 2, 'omitnan');
            resultsS.mean_F = mean(simF_absM, 2, 'omitnan');
            resultsS.mean_Alpha = mean(simAlpha_M, 2, 'omitnan');
            resultsS.mean_Q_rate = mean(simQ_rate_M, 2, 'omitnan');
            % [!!! 新增: 汇总新增变量 !!!]
            resultsS.mean_Y = mean(simY_absM, 2, 'omitnan');
            resultsS.mean_CoH = mean(simCoH_absM, 2, 'omitnan');
            resultsS.mean_Pension = mean(simPension_absM, 2, 'omitnan');
            resultsS.mean_W_to_Y_ratio = resultsS.mean_W ./ resultsS.mean_Y;

            resultsS.ages = (1:cS.tn) + 18 - 1;

            % 5. 可视化
            fprintf('--> 正在生成可视化图表...\n');
            fig = figure('Name', 'Lifecycle Simulation Results', 'Position', [100, 100, 1200, 800]);

            % Panel 1: 财富、消费与收入
            subplot(2,2,1);
            plot(resultsS.ages, resultsS.mean_W, 'b-', 'LineWidth', 2);
            hold on;
            plot(resultsS.ages, resultsS.mean_C, 'r--', 'LineWidth', 2);
            plot(resultsS.ages, resultsS.mean_F, 'g-.', 'LineWidth', 2);
            plot(resultsS.ages, resultsS.mean_Y, 'k:', 'LineWidth', 2); % 新增收入曲线
            xline(cS.tr + 18 - 1, 'k:', 'Retirement', 'LineWidth', 1.5);
            hold off;
            title('Mean Lifecycle Profiles', 'Interpreter', 'latex');
            xlabel('Age', 'Interpreter', 'latex');
            ylabel('Value (absolute)', 'Interpreter', 'latex');
            legend({'Wealth ($\bar{W}$)', 'Consumption ($\bar{C}$)', ...
                'Pension Fund ($\bar{F}$)', 'Income ($\bar{Y}$)'}, ...
                'Location', 'northwest', 'Interpreter', 'latex');
            grid on;

            % Panel 2: 投资组合决策
            subplot(2,2,2);
            plot(resultsS.ages, resultsS.mean_Alpha, 'm-', 'LineWidth', 2);
            xline(cS.tr + 18 - 1, 'k:', 'Retirement', 'LineWidth', 1.5);
            title('Mean Risky Asset Allocation', 'Interpreter', 'latex');
            xlabel('Age', 'Interpreter', 'latex');
            ylabel('Risky Share ($\bar{\alpha}$)', 'Interpreter', 'latex');
            ylim([-0.05, 1.05]);
            grid on;

            % Panel 3: 养老金缴费决策
            subplot(2,2,3);
            plot(resultsS.ages, resultsS.mean_Q_rate, 'k-', 'LineWidth', 2);
            xline(cS.tr + 18 - 1, 'k:', 'Retirement', 'LineWidth', 1.5);
            title('Mean Pension Contribution Rate', 'Interpreter', 'latex');
            xlabel('Age', 'Interpreter', 'latex');
            ylabel('Contribution Rate ($\bar{q}$)', 'Interpreter', 'latex');
            ylim([-0.01, cS.Q_max + 0.01]);
            grid on;

            % Panel 4: 财富/收入 比率
            subplot(2,2,4);
            plot(resultsS.ages, resultsS.mean_W_to_Y_ratio, 'c-', 'LineWidth', 2);
            xline(cS.tr + 18 - 1, 'k:', 'Retirement', 'LineWidth', 1.5);
            title('Mean Wealth-to-Income Ratio', 'Interpreter', 'latex');
            xlabel('Age', 'Interpreter', 'latex');
            ylabel('$\bar{W} / \bar{Y}$', 'Interpreter', 'latex');
            grid on;

            sgtitle('Lifecycle Simulation Results (VFI method)', 'FontSize', 16, 'FontWeight', 'bold');

            % 保存图表
            output_dir = 'fig';
            if ~exist(output_dir, 'dir'), mkdir(output_dir); end
            output_filename = fullfile(output_dir, 'Lifecycle_Simulation_VFI.png');
            fprintf('--> 正在将模拟结果图表保存至: %s\n', output_filename);
            print(fig, output_filename, '-dpng', '-r300');

            fprintf('模拟结束。\n');
        end
        
                function resultsS = simulation_cocco(polS, paramS, cS)
            % =========================================================================
            % == 函数: simulation_cocco (版本 v1.0)
            % == 目的: 基于 VFI_cocco_for 的求解结果，进行无养老金账户的
            % ==       生命周期模拟，以对标 life_cycle.m 的经济环境。
            % == 核心修改:
            % ==   1. 移除了所有与 F (养老金账户) 相关的变量和逻辑。
            % ==   2. 状态变量只有 W (流动性财富)。插值函数变为一维。
            % ==   3. 调整了工作期和退休期的状态转移方程。
            % =========================================================================

            fprintf('开始进行 Cocco/Gomes 模型的生命周期模拟...\n');
            
            % 1. 初始化
            nsim = 2000;
            rng(2024);

            % 初始化存储矩阵 (无 F 维度)
            simW_normM = zeros(cS.tn, nsim); % 归一化财富
            Y_levelM = ones(cS.tn, nsim);   % 收入永久成分水平
            
            simC_absM = zeros(cS.tn, nsim);      % 绝对消费
            simW_absM = zeros(cS.tn, nsim);      % 绝对财富
            simAlpha_M = zeros(cS.tn, nsim);     % 风险资产配置比例
            
            simY_absM = zeros(cS.tn, nsim);      % 绝对收入 (包括暂时性冲击)
            simCoH_absM = zeros(cS.tn, nsim);    % 绝对可支配资源
            
            % 设置初始状态
            simW_normM(1, :) = cS.w_min; % 从最低财富开始
            Y_levelM(1,:) = exp(paramS.g_t_V(1)); % 初始化收入水平

            % 2. 生成生命周期内的随机冲击
            % 为了更贴近 Cocco/Gomes，这里使用正态分布随机数
            z_shocks_M = randn(cS.tr-1, nsim) * paramS.sigma_z; % 只有工作期有永久冲击
            eps_shocks_M = randn(cS.tn-1, nsim) * paramS.sigma_epsilon;

            % 为了兼容，我们仍然使用 Tauchen 方法的节点
            % [您可以选择切换回randn，或继续使用离散节点以保持和VFI一致]
            rand_z = rand(cS.tr-1, nsim);
            rand_eps = rand(cS.tn-1, nsim);
            z_prob_cum = cumsum(paramS.zProbs);
            eps_prob_cum = cumsum(paramS.epsilonProbs);
            z_idx_M = zeros(cS.tr-1, nsim);
            eps_idx_M = zeros(cS.tn-1, nsim);
            for t_sim = 1:(cS.tr-1)
                for i_sim = 1:nsim
                    z_idx_M(t_sim, i_sim) = find(rand_z(t_sim, i_sim) <= z_prob_cum, 1, 'first');
                end
            end
            for t_sim = 1:(cS.tn-1)
                 for i_sim = 1:nsim
                    eps_idx_M(t_sim, i_sim) = find(rand_eps(t_sim, i_sim) <= eps_prob_cum, 1, 'first');
                end
            end
            
            % 3. 逐期迭代模拟
            fprintf('--> 正在逐期模拟个体决策 (Cocco/Gomes model)...\n');
            
            for t = 1:cS.tn-1
                w_t_normV = simW_normM(t, :);
                w_t_normV = max(min(w_t_normV, cS.w_max), cS.w_min);

                % 插值函数变为一维 (interp1)
                c_t_normV = interp1(cS.wGridV, polS(t).c, w_t_normV, 'linear', 'extrap');
                alpha_t_V = interp1(cS.wGridV, polS(t).alpha, w_t_normV, 'linear', 'extrap');
                
                Y_level_t_V = Y_levelM(t, :);

                if t < cS.tr % 工作期
                    hat_Y_t = 1.0; % 归一化暂时性收入为1
                    % CoH = W + Y (无养老金缴费)
                    CoH_norm = w_t_normV + (1-cS.tau_y) * hat_Y_t;
                    
                    simY_absM(t,:) = hat_Y_t .* Y_level_t_V;
                    simCoH_absM(t,:) = CoH_norm .* Y_level_t_V;
                    
                    % 确保消费不超过CoH
                    c_t_normV = min(c_t_normV, CoH_norm);
                    c_t_normV(c_t_normV < 0) = 1e-6;
                    
                    w_prime_norm = CoH_norm - c_t_normV;
                    
                    % 使用离散冲击节点
                    R_shock_V = (cS.R_f + paramS.mu + paramS.epsilonNodes(eps_idx_M(t, :)))';
                    R_portfolio_V = (1 - alpha_t_V) .* cS.R_f + alpha_t_V .* R_shock_V;
                    
                    w_invested_abs = w_prime_norm .* Y_level_t_V;
                    w_next_abs = w_invested_abs .* R_portfolio_V;
                    
                    G_t = exp(paramS.g_t_V(t+1) - paramS.g_t_V(t));
                    z_shocks_V = paramS.zGridV(z_idx_M(t, :))';
                    total_growth_factor_V = G_t * exp(z_shocks_V);
                    
                    Y_levelM(t+1, :) = Y_level_t_V .* total_growth_factor_V;
                    
                    simW_normM(t+1, :) = w_next_abs ./ Y_levelM(t+1, :);

                else % 退休期
                    % CoH = W + Y_bar (基础养老金)
                    CoH_norm = w_t_normV + cS.Y_bar;
                    
                    % 退休期的 "收入" 是指基础养老金
                    simY_absM(t, :) = cS.Y_bar .* ones(1, nsim); 
                    simCoH_absM(t, :) = CoH_norm .* Y_level_t_V;

                    c_t_normV = min(c_t_normV, CoH_norm);
                    c_t_normV(c_t_normV < 0) = 1e-6;

                    w_prime_norm = CoH_norm - c_t_normV;
                    
                    R_shock_V = (cS.R_f + paramS.mu + paramS.epsilonNodes(eps_idx_M(t, :)))';
                    R_portfolio_V = (1 - alpha_t_V) .* cS.R_f + alpha_t_V .* R_shock_V;

                    w_invested_abs = w_prime_norm .* Y_level_t_V;
                    w_next_abs = w_invested_abs .* R_portfolio_V;
                    
                    % 收入永久成分不再增长
                    Y_levelM(t+1, :) = Y_levelM(t, :);

                    simW_normM(t+1, :) = w_next_abs ./ Y_levelM(t+1, :);
                end

                simC_absM(t, :) = c_t_normV .* Y_level_t_V;
                simW_absM(t, :) = w_t_normV .* Y_level_t_V;
                simAlpha_M(t, :) = alpha_t_V;
            end

            % 处理最后一期
            w_tn_abs = simW_normM(cS.tn, :) .* Y_levelM(cS.tn, :);
            simW_absM(cS.tn, :) = w_tn_abs;
            % 终期消费 = 剩余财富 + 最后一笔基础养老金
            simC_absM(cS.tn, :) = w_tn_abs + cS.Y_bar .* Y_levelM(cS.tn, :);
            simCoH_absM(cS.tn, :) = simC_absM(cS.tn, :);
            simY_absM(cS.tn, :) = cS.Y_bar .* Y_levelM(cS.tn, :);
            simAlpha_M(cS.tn, :) = 0;
            
            % 4. 汇总并打包结果
            fprintf('--> 正在汇总 Cocco/Gomes 模型模拟结果...\n');
            resultsS.mean_C = mean(simC_absM, 2, 'omitnan');
            resultsS.mean_W = mean(simW_absM, 2, 'omitnan');
            resultsS.mean_Alpha = mean(simAlpha_M, 2, 'omitnan');
            resultsS.mean_Y = mean(simY_absM, 2, 'omitnan');
            resultsS.mean_CoH = mean(simCoH_absM, 2, 'omitnan');
            % 计算财富收入比时，需要处理退休期收入为常数的情况
            valid_Y_idx = resultsS.mean_Y > 1e-6;
            resultsS.mean_W_to_Y_ratio = nan(cS.tn, 1);
            resultsS.mean_W_to_Y_ratio(valid_Y_idx) = resultsS.mean_W(valid_Y_idx) ./ resultsS.mean_Y(valid_Y_idx);
            
            resultsS.ages = (1:cS.tn) + 18 - 1;
            
            % 5. 可视化
            fprintf('--> 正在生成 Cocco/Gomes 模型可视化图表...\n');
            fig = figure('Name', 'Lifecycle Simulation (Cocco/Gomes Model)', 'Position', [100, 100, 1200, 800]);
            
            % Panel 1: 财富、消费与收入
            subplot(2,2,1);
            plot(resultsS.ages, resultsS.mean_W, 'b-', 'LineWidth', 2);
            hold on;
            plot(resultsS.ages, resultsS.mean_C, 'r--', 'LineWidth', 2);
            plot(resultsS.ages, resultsS.mean_Y, 'k:', 'LineWidth', 2);
            xline(cS.tr + 18 - 1, 'k:', 'Retirement', 'LineWidth', 1.5);
            hold off;
            title('Mean Lifecycle Profiles', 'Interpreter', 'latex');
            xlabel('Age', 'Interpreter', 'latex');
            ylabel('Value (absolute)', 'Interpreter', 'latex');
            legend({'Wealth ($\bar{W}$)', 'Consumption ($\bar{C}$)', 'Income ($\bar{Y}$)'}, ...
                'Location', 'northwest', 'Interpreter', 'latex');
            grid on;

            % Panel 2: 投资组合决策
            subplot(2,2,2);
            plot(resultsS.ages, resultsS.mean_Alpha, 'm-', 'LineWidth', 2);
            xline(cS.tr + 18 - 1, 'k:', 'Retirement', 'LineWidth', 1.5);
            title('Mean Risky Asset Allocation', 'Interpreter', 'latex');
            xlabel('Age', 'Interpreter', 'latex');
            ylabel('Risky Share ($\bar{\alpha}$)', 'Interpreter', 'latex');
            ylim([-0.05, 1.05]);
            grid on;
            
            % Panel 3: 财富/收入 比率
            subplot(2,2,3);
            plot(resultsS.ages, resultsS.mean_W_to_Y_ratio, 'c-', 'LineWidth', 2);
            xline(cS.tr + 18 - 1, 'k:', 'Retirement', 'LineWidth', 1.5);
            title('Mean Wealth-to-Income Ratio', 'Interpreter', 'latex');
            xlabel('Age', 'Interpreter', 'latex');
            ylabel('$\bar{W} / \bar{Y}$', 'Interpreter', 'latex');
            grid on;
            
            % (Panel 4 is left empty)
            subplot(2,2,4);
            axis off;

            sgtitle('Lifecycle Simulation Results (Cocco/Gomes Model via VFI\_cocco\_for)', 'FontSize', 16, 'FontWeight', 'bold');
            
            output_dir = 'fig';
            if ~exist(output_dir, 'dir'), mkdir(output_dir); end
            output_filename = fullfile(output_dir, 'Lifecycle_Simulation_Cocco.png');
            fprintf('--> 正在将模拟结果图表保存至: %s\n', output_filename);
            print(fig, output_filename, '-dpng', '-r300');

            fprintf('Cocco/Gomes 模型模拟结束。\n');
        end
        
        % ---------------测试函数-----------------
        function test_for_matrix()
            % =========================================================================
            % == 函数: test_for_matrix (版本 v6.0 - 误差报告版)
            % == 目的: 对比矩阵与循环版本结果，并将所有时期的最大误差系统性地
            % ==      记录在一个表格中，打印并保存为文件。
            % =========================================================================
            fprintf('开始测试 VFI (Epstein-Zin) 版本的矩阵与循环等价性...\n');

            % 1. 获取一套参数，并使用较小的网格进行快速测试
            [cS, paramS] = household.parameter_values();
            cS.nW = 8; cS.nF = 4; cS.nZ = 3; cS.nEpsilon = 3;
            cS.nC = 5; cS.nQ = 3; cS.nAlpha = 3;
            cS.wGridV = linspace(0.1, 10, cS.nW)';
            cS.fGridV = linspace(0, 10, cS.nF)';

            [z_nodes, z_trans_mat] = household.tauchen(cS.nZ, 0, paramS.sigma_z, 0, 2);
            paramS.zGridV = z_nodes;
            paramS.zProbs = z_trans_mat(1, :);

            [eps_nodes, eps_prob_mat] = household.tauchen(cS.nEpsilon, 0, paramS.sigma_epsilon, 0, 2);
            paramS.epsilonNodes = eps_nodes;
            paramS.epsilonProbs = eps_prob_mat(1, :);
            paramS.R_shock_V = cS.R_f + paramS.mu + paramS.epsilonNodes;

            cS.alphaGridV = linspace(0, 1, cS.nAlpha);
            cS.qGridV = linspace(0, cS.Q_max, cS.nQ);

            fprintf('测试参数: nW=%d, nF=%d, nZ=%d, nEpsilon=%d\n', cS.nW, cS.nF, cS.nZ, cS.nEpsilon);
            fprintf('生命周期: T=%d, K=%d\n\n', cS.tn, cS.tr);

            % 2. 运行 For-Loop 版本并计时
            fprintf('--> 正在运行 For-Loop 版本 (VFI_PPS_for)...\n');
            tic;
            [polS_for, valS_for] = household.VFI_PPS_for(paramS, cS);
            time_for = toc;
            fprintf('    For-Loop 版本完成，耗时: %.4f 秒\n\n', time_for);

            % 3. 运行矩阵化版本并计时
            fprintf('--> 正在运行矩阵化版本 (VFI_PPS)...\n');
            tic;
            [polS_mat, valS_mat] = household.VFI_PPS(paramS, cS);
            time_mat = toc;
            fprintf('    矩阵化版本完成，耗时: %.4f 秒\n\n', time_mat);

            % 4. 逐期比较结果并记录误差
            fprintf('--> 正在逐期比较两个版本的结果...\n');

            num_periods = cS.tn - cS.tb;
            error_table = table('Size', [num_periods, 5], ...
                'VariableTypes', {'double', 'double', 'double', 'double', 'double'}, ...
                'VariableNames', {'t', 'Max_Val_Diff', 'Max_C_Diff', 'Max_Q_Diff', 'Max_Alpha_Diff'});

            row_idx = 1;
            for t = cS.tn-1 : -1 : cS.tb
                % 计算值函数差异
                val_diff_t = abs(valS_for(:,:,t) - valS_mat(:,:,t));
                max_diff_t = max(val_diff_t(:), [], 'omitnan');

                % 计算策略函数差异
                c_diff = abs(polS_for(t).c - polS_mat(t).c);
                max_c_diff = max(c_diff(:), [], 'omitnan');

                q_diff = abs(polS_for(t).q_rate - polS_mat(t).q_rate);
                max_q_diff = max(q_diff(:), [], 'omitnan');

                a_diff = abs(polS_for(t).alpha - polS_mat(t).alpha);
                max_a_diff = max(a_diff(:), [], 'omitnan');

                % 将当前时期的所有最大误差存入表格
                error_table(row_idx, :) = {t, max_diff_t, max_c_diff, max_q_diff, max_a_diff};
                row_idx = row_idx + 1;
            end

            % 5. 整理、显示并保存结果
            error_table = sortrows(error_table, 't'); % 按时间 t 升序排列

            fprintf('\n--> 逐期最大误差对比表:\n');
            disp(error_table);

            output_filename = 'debug\comparison_errors.txt';
            try
                writetable(error_table, output_filename, 'Delimiter', '\t');
                fprintf('\n    误差对比表已成功保存至文件: %s\n\n', output_filename);
            catch ME
                fprintf('\n    保存文件时出错: %s\n', ME.message);
                fprintf('    请检查工作目录的写入权限。\n\n');
            end

            % 6. 性能总结
            fprintf('--> 性能总结:\n');
            if time_mat > 1e-9
                fprintf('    矩阵化版本比 For-Loop 版本快 %.2f 倍。\n', time_for / time_mat);
            end
            fprintf('测试结束。\n');
        end


        function test_baseline()
            % =========================================================================
            % == 函数: test_baseline (版本 v5.4 - 增加结果保存)
            % == 目的: 求解EZ模型，保存结果，可视化策略函数，并运行模拟。
            % == 核心修改:
            % ==   1. 在VFI求解后，新增了 `save` 命令，将求解结果 `polS` 和 `valS`
            % ==      保存到 `debug/vfi_results.mat` 文件中。
            % =========================================================================

            fprintf('开始运行基准模型(Epstein-Zin)求解与可视化...\n');

            output_dir_data = 'debug';
            if ~exist(output_dir_data, 'dir')
                mkdir(output_dir_data);
            end
            output_filename_data = fullfile(output_dir_data, 'vfi_results.mat');

            if ~exist(output_filename_data,'file')
                % 1. 获取基准参数并求解模型
                [cS, paramS] = household.parameter_values();
                fprintf('参数加载完毕。模型维度: nW=%d, nF=%d, nZ=%d, nEpsilon=%d\n', cS.nW, cS.nF, cS.nZ, cS.nEpsilon);
                fprintf('--> 正在求解 VFI 模型 (VFI_PPS)...\n');
                tic;
                [polS, valS] = household.VFI_PPS(paramS, cS);
                time_vfi = toc;
                fprintf('    模型求解完成，耗时: %.4f 秒\n\n', time_vfi);

                % [!!! 新增步骤: 保存求解结果 !!!]
                fprintf('--> 正在保存 VFI 求解结果...\n');

                try
                    save(output_filename_data, 'polS', 'valS', 'paramS', 'cS');
                    fprintf('    求解结果已成功保存至: %s\n\n', output_filename_data);
                catch ME
                    fprintf('    保存文件时出错: %s\n', ME.message);
                    fprintf('    请检查工作目录的写入权限。\n\n');
                end
            else
                load( output_filename_data)
            end

            % 2. 可视化准备
            % fprintf('--> 正在生成整合的策略函数可视化图表...\n');
            % tb_age = 18;
            % plot_ages = [25, 45, 60, 75];
            % plot_t = plot_ages - tb_age + 1;
            % [f_grid_plot, w_grid_plot] = meshgrid(cS.fGridV, cS.wGridV);
            %
            % % 3. 创建并绘制策略函数图表
            % fig = figure('Name', 'Lifecycle Model Results (Epstein-Zin)', 'Position', [50, 50, 1800, 1600]);
            % sgtitle('Lifecycle Model Policy Functions (Epstein-Zin, Normalized)', 'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'latex');
            %
            % for i = 1:length(plot_t)
            %     t = plot_t(i);
            %     age = plot_ages(i);
            %
            %     subplot_idx = i;
            %     subplot(4, length(plot_t), subplot_idx);
            %     val_slice = valS(:, :, t);
            %     surf(w_grid_plot, f_grid_plot, val_slice);
            %     title(sprintf('Age = %d', age), 'Interpreter', 'latex');
            %     view(135, 30);
            %     xlabel('$\tilde{W}$', 'Interpreter', 'latex');
            %     if i == 1, ylabel('$\tilde{F}$', 'Interpreter', 'latex'); end
            %     zlabel('$V(\tilde{W}, \tilde{F})$', 'Interpreter', 'latex');
            %
            %     subplot_idx = length(plot_t) + i;
            %     subplot(4, length(plot_t), subplot_idx);
            %     c_slice = polS(t).c;
            %     surf(w_grid_plot, f_grid_plot, c_slice);
            %     view(135, 30);
            %     xlabel('$\tilde{W}$', 'Interpreter', 'latex');
            %     if i == 1, ylabel('$\tilde{F}$', 'Interpreter', 'latex'); end
            %     zlabel('$\tilde{C}$', 'Interpreter', 'latex');
            %
            %     subplot_idx = 2 * length(plot_t) + i;
            %     subplot(4, length(plot_t), subplot_idx);
            %     if t < cS.tr
            %         q_slice = polS(t).q_rate;
            %         surf(w_grid_plot, f_grid_plot, q_slice);
            %         zlim([0, cS.Q_max]);
            %         view(135, 30);
            %         xlabel('$\tilde{W}$', 'Interpreter', 'latex');
            %         if i == 1, ylabel('$\tilde{F}$', 'Interpreter', 'latex'); end
            %         zlabel('$q$', 'Interpreter', 'latex');
            %     else
            %         text(0.5, 0.5, 'Retired', 'HorizontalAlignment', 'center', 'FontSize', 12, 'Interpreter', 'latex');
            %         axis off;
            %     end
            %
            %     subplot_idx = 3 * length(plot_t) + i;
            %     subplot(4, length(plot_t), subplot_idx);
            %     alpha_slice = polS(t).alpha;
            %     surf(w_grid_plot, f_grid_plot, alpha_slice);
            %     zlim([0, 1]);
            %     view(135, 30);
            %     xlabel('$\tilde{W}$', 'Interpreter', 'latex');
            %     if i == 1, ylabel('$\tilde{F}$', 'Interpreter', 'latex'); end
            %     zlabel('$\alpha$', 'Interpreter', 'latex');
            % end
            %
            % % 4. 保存策略函数图表
            % output_dir_fig = 'fig';
            % if ~exist(output_dir_fig, 'dir'), mkdir(output_dir_fig); end
            % output_filename_fig = fullfile(output_dir_fig, 'Lifecycle_Policies_Grid_EZ.png');
            % fprintf('--> 正在将策略函数图表保存为 %s (600 DPI)...\n', output_filename_fig);
            % print(fig, output_filename_fig, '-dpng', '-r600');
            % fprintf('    策略函数可视化图表生成并保存完毕。\n\n');

            % =============================================================
            % == 步骤: 运行并可视化生命周期模拟
            % =============================================================
            fprintf('--> 开始基于 VFI 结果进行生命周期模拟...\n');
            tic;
            household.simulation(polS, paramS, cS); % 调用模拟函数
            time_sim = toc;
            fprintf('    模拟与可视化完成，耗时: %.4f 秒\n\n', time_sim);

            fprintf('test_baseline 执行结束。\n');
        end

                function test_baseline_cocco()
            % =========================================================================
            % == 函数: test_baseline_cocco (版本 v1.0)
            % == 目的: 求解 Cocco/Gomes 模型，保存结果，并运行模拟。
            % ==       这是一个完整的复刻流程的执行入口。
            % == 核心调用:
            % ==   1. VFI求解: VFI_cocco_for
            % ==   2. 模拟与可视化: simulation_cocco
            % =========================================================================

            fprintf('开始运行 Cocco/Gomes 模型基准求解与模拟...\n');

            % 设定结果保存路径和文件名
            output_dir_data = 'debug';
            if ~exist(output_dir_data, 'dir')
                mkdir(output_dir_data); 
            end
            output_filename_data = fullfile(output_dir_data, 'vfi_results_cocco.mat');

            % 检查是否已有计算结果，避免重复计算
            if ~exist(output_filename_data, 'file')
                % 1. 获取基准参数并求解模型
                [cS, paramS] = household.parameter_values();
                fprintf('参数加载完毕。模型维度: nW=%d, nZ=%d, nEpsilon=%d\n', cS.nW, cS.nZ, cS.nEpsilon);
                
                % [!!! 核心调用: 求解无养老金账户的模型 !!!]
                fprintf('--> 正在求解 VFI 模型 (VFI_cocco_for)...\n');
                tic;
                [polS, valS] = household.VFI_cocco_for(paramS, cS);
                time_vfi = toc;
                fprintf('    模型求解完成，耗时: %.4f 秒\n\n', time_vfi);

                % 2. 保存求解结果
                fprintf('--> 正在保存 VFI 求解结果...\n');
                try
                    save(output_filename_data, 'polS', 'valS', 'paramS', 'cS');
                    fprintf('    求解结果已成功保存至: %s\n\n', output_filename_data);
                catch ME
                    fprintf('    保存文件时出错: %s\n', ME.message);
                end
            else
                fprintf('--> 检测到已有的 VFI 结果，直接加载: %s\n\n', output_filename_data);
                load(output_filename_data);
            end

            % (可视化策略函数的代码被注释掉，因为 simulation_cocco 会生成更直观的生命周期图)
            % 您可以根据需要取消注释来查看策略曲面
            %{
            fprintf('--> 正在生成策略函数可视化图表...\n');
            tb_age = 18;
            plot_ages = [25, 45, 60, 75];
            plot_t = plot_ages - tb_age + 1;
            
            fig = figure('Name', 'Cocco/Gomes Model Policy Functions', 'Position', [50, 50, 1600, 800]);
            sgtitle('Cocco/Gomes Policy Functions (Normalized)', 'FontSize', 16, 'FontWeight', 'bold');
            
            for i = 1:length(plot_t)
                t = plot_t(i);
                age = plot_ages(i);
            
                % Value Function
                subplot(3, length(plot_t), i);
                plot(cS.wGridV, valS(:, t), 'LineWidth', 2);
                title(sprintf('Value Function at Age %d', age));
                xlabel('Normalized Wealth ($\tilde{W}$)');
                ylabel('$V(\tilde{W})$');
                grid on;
            
                % Consumption Policy
                subplot(3, length(plot_t), i + length(plot_t));
                plot(cS.wGridV, polS(t).c, 'LineWidth', 2);
                title(sprintf('Consumption Policy at Age %d', age));
                xlabel('Normalized Wealth ($\tilde{W}$)');
                ylabel('$\tilde{C}$');
                grid on;
            
                % Alpha Policy
                subplot(3, length(plot_t), i + 2*length(plot_t));
                plot(cS.wGridV, polS(t).alpha, 'LineWidth', 2);
                title(sprintf('Risky Share Policy at Age %d', age));
                xlabel('Normalized Wealth ($\tilde{W}$)');
                ylabel('$\alpha$');
                ylim([-0.05, 1.05]);
                grid on;
            end
            %}

            % 3. 运行并可视化生命周期模拟
            % [!!! 核心调用: 运行无养老金账户的模拟 !!!]
            fprintf('--> 开始基于 VFI 结果进行生命周期模拟 (Cocco/Gomes)...\n');
            tic;
            household.simulation_cocco(polS, paramS, cS); % 调用 Cocco 专属模拟函数
            time_sim = toc;
            fprintf('    模拟与可视化完成，耗时: %.4f 秒\n\n', time_sim);

            fprintf('test_baseline_cocco 执行结束。\n');
        end
        % --- 辅助函数 ---
        function [y_grid_out, trProbM_out] = tauchen(N, rho, sigma, mu, m)
            if N == 1, y_grid_out = 0; trProbM_out = 1; return; end
            std_y = sqrt(sigma^2 / (1-rho^2)); y_max = m*std_y; y_min = -y_max;
            y = linspace(y_min, y_max, N); d = y(2)-y(1);
            trProbM_out = zeros(N,N);
            for j=1:N, for k=1:N
                    m_k = rho*y(j) + mu;
                    if k==1, trProbM_out(j,k) = normcdf((y(1)-m_k+d/2)/sigma);
                    elseif k==N, trProbM_out(j,k) = 1 - normcdf((y(N)-m_k-d/2)/sigma);
                    else, trProbM_out(j,k) = normcdf((y(k)-m_k+d/2)/sigma) - normcdf((y(k)-m_k-d/2)/sigma); end
            end, end
        y_grid_out = y(:);
        end

    end
end