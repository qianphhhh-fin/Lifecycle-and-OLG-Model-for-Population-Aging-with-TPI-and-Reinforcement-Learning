classdef utils
    methods (Static)

        function cS = set_parameters()
            % =================================================================
            % == 功能: 设置模型所有参数和网格 (V2.0 - 使用Tauchen)
            % == 核心: 所有随机过程统一使用 Tauchen 方法离散化。
            % =================================================================

            % --- 1. 生命周期与人口结构 ---
            cS.tb_age = 20;         % 工作起始年龄
            cS.tr_age = 66;         % 退休年龄
            cS.td_age = 100;        % 最高寿命
            cS.tb = 1;              % 模型起始期
            cS.tr = cS.tr_age - cS.tb_age + 1; % 模型退休期
            cS.tn = cS.td_age - cS.tb_age + 1; % 模型总期数

            % 条件生存概率
            survprob_data = [0.99845, 0.99839, 0.99833, 0.9983, 0.99827, 0.99826, 0.99824, 0.9982, 0.99813, 0.99804, 0.99795, 0.99785, 0.99776, 0.99766, 0.99755, 0.99743, 0.9973, 0.99718, 0.99707, 0.99696, 0.99685, 0.99672, 0.99656, 0.99635, 0.9961, 0.99579, 0.99543, 0.99504, 0.99463, 0.9942, 0.9937, 0.99311, 0.99245, 0.99172, 0.99091, 0.99005, 0.98911, 0.98803, 0.9868, 0.98545, 0.98409, 0.9827, 0.98123, 0.97961, 0.97786, 0.97603, 0.97414, 0.97207, 0.9697, 0.96699, 0.96393, 0.96055, 0.9569, 0.9531, 0.94921, 0.94508, 0.94057, 0.9357, 0.93031, 0.92424, 0.91717, 0.90922, 0.90089, 0.89282, 0.88503, 0.87622, 0.86576, 0.8544, 0.8423, 0.82942, 0.8154, 0.80002, 0.78404, 0.76842, 0.75382, 0.73996, 0.72464, 0.71057, 0.6961, 0.6809];
            cS.pi_pathV = ones(cS.tn, 1);
            len_surv = length(survprob_data);
            cS.pi_pathV(2:min(cS.tn, len_surv+1)) = survprob_data(1:min(cS.tn-1, len_surv));

            % --- 2. 偏好 (Epstein-Zin) ---
            cS.beta = 0.97;
            cS.gamma = 10.0;
            cS.psi = 0.5;
            cS.theta = (1.0 - cS.gamma) / (1.0 - 1.0/cS.psi);
            cS.psi_1 = 1.0 - 1.0/cS.psi;
            cS.psi_2 = 1.0 / cS.psi_1;

            % --- 3. 收入过程 ---
            aa = -2.170042 + 2.700381; b1 = 0.16818; b2 = -0.0323371/10; b3 = 0.0019704/100;
            g_t = zeros(cS.tr, 1);
            for t = 1:cS.tr
                age = t + cS.tb_age - 1;
                g_t(t) = aa + b1*age + b2*age^2 + b3*age^3;
            end
            cS.G_pathV = ones(cS.tn, 1);
            cS.G_pathV(1:cS.tr-1) = exp(g_t(2:end) - g_t(1:end-1));
            cS.sigma_u = 0.1;
            cS.sigma_z = 0.1;

            % --- 4. 养老金与税收 ---
            cS.tau_y = 0.06;
            cS.tau_q = 0.03;
            cS.Q_max = 0.1;
            cS.Y_bar_frac = 0.68212;

            % --- 5. 资产与回报率 ---
            cS.R_f = 1.015;
            cS.mu = 0.04;
            cS.sigma_eps = 0.2;
            cS.R_p = cS.R_f;

            % --- 6. 数值求解网格 ---
            cS.nW = 51;
            cS.nF = 31;
            cS.w_min = 0.25; cS.w_max = 200;
            cS.f_min = 0.0; cS.f_max = 100;
            cS.wGridV = exp(linspace(log(cS.w_min), log(cS.w_max), cS.nW))';
            cS.fGridV = exp(linspace(log(cS.f_min+1e-6), log(cS.f_max), cS.nF-1))';
            cS.fGridV = [0; cS.fGridV];
            cS.nC = 21;
            cS.nQ = 21;
            cS.nAlpha = 51;
            cS.alphaGridV = linspace(0, 1, cS.nAlpha)';
            cS.qGridV = linspace(0, cS.Q_max, cS.nQ)';
            cS.savingsFracGridV = linspace(0.001, 0.999, cS.nC)';
            
            % --- [!!! 核心修改: 使用 Tauchen 离散化所有随机冲击 !!!] ---
            cS.nShocks = 5; % 离散化节点数
            m_std = 3;      % 网格宽度为m_std个标准差

            % 对于i.i.d.冲击，rho=0
            [zNodes, zProbMat] = utils.tauchen(cS.nShocks, 0, cS.sigma_z, 0, m_std);
            [uNodes, uProbMat] = utils.tauchen(cS.nShocks, 0, cS.sigma_u, 0, m_std);
            [epsNodes, epsProbMat] = utils.tauchen(cS.nShocks, 0, cS.sigma_eps, 0, m_std);
            
            cS.zNodes = zNodes;
            cS.uNodes = uNodes;
            cS.epsNodes = epsNodes;
            
            % 因为是i.i.d. (rho=0), 转移矩阵的每一行都是相同的，即稳态分布
            cS.shockProbs = zProbMat(1, :)'; % 取第一行作为概率向量即可
            
            % 基于离散化的冲击节点计算风险回报率
            cS.R_shock_V = cS.R_f + cS.mu + cS.epsNodes;
            % --- [!!! 修正结束 !!!] ---
        end


        function [polS, valS] = value_function_iteration(cS)
            % =================================================================
            % == 功能: 逆向归纳求解模型 (V2.0 - 修正工作期逻辑)
            % == 核心:
            % ==  1. 工作期VFI严格遵循 E[max V] 原则，最大化在对u_t求期望内部。
            % ==  2. 工作期策略函数变为三维，依赖于(w, f, u)。
            % ==  3. 退休期逻辑不变，对每个F_K解1D问题。
            % =================================================================

            % --- 1. 初始化 ---
            valS = -inf(cS.nW, cS.nF, cS.tn);
            % 策略函数现在是cell数组，工作期将存放3D矩阵
            polS.c = cell(cS.tn, 1);
            polS.alpha = cell(cS.tn, 1);
            polS.q = cell(cS.tn, 1);

            % --- 2. 终期 T ---
            [w_grid_T, f_grid_T] = ndgrid(cS.wGridV, cS.fGridV);
            p_bar_T = (1 - cS.tau_q) * f_grid_T / (cS.tn - cS.tr + 1);
            final_c = w_grid_T + cS.Y_bar_frac + p_bar_T;
            valS(:, :, cS.tn) = final_c; % V_T = C_T
            
            polS.c{cS.tn} = final_c;
            polS.alpha{cS.tn} = zeros(cS.nW, cS.nF);
            polS.q{cS.tn} = zeros(cS.nW, cS.nF);

            % --- 3. 退休期 (t = T-1, ..., K) ---
            % 这部分逻辑是正确的，因为退休后没有暂时性收入冲击，无需改动
            fprintf('   Solving retirement period (t = %d to %d)...\n', cS.tn-1, cS.tr);
            val_retire_all_f = zeros(cS.nW, cS.nF, cS.tn - cS.tr + 1);
            val_retire_all_f(:, :, end) = valS(:, :, cS.tn);

            for t = cS.tn-1 : -1 : cS.tr
                fprintf('   t = %d\n', t);
                V_t_plus_1_interp_handles = cell(cS.nF, 1);
                for iF = 1:cS.nF
                    V_t_plus_1_interp_handles{iF} = griddedInterpolant(cS.wGridV, val_retire_all_f(:, iF, t - cS.tr + 2), 'linear', 'linear');
                end
                
                temp_val_t = zeros(cS.nW, cS.nF);
                temp_c_pol_t = zeros(cS.nW, cS.nF);
                temp_alpha_pol_t = zeros(cS.nW, cS.nF);

                parfor iF = 1:cS.nF
                    f_k = cS.fGridV(iF);
                    p_bar = (1 - cS.tau_q) * f_k / (cS.tn - cS.tr + 1);
                    v_interp = V_t_plus_1_interp_handles{iF};
                    
                    val_vec = zeros(cS.nW, 1);
                    c_pol_vec = zeros(cS.nW, 1);
                    alpha_pol_vec = zeros(cS.nW, 1);

                    for iW = 1:cS.nW
                        w_t = cS.wGridV(iW);
                        coh = w_t + cS.Y_bar_frac + p_bar;
                        
                        max_v = -inf;
                        best_c = coh; best_alpha = 0;

                        for iAlpha = 1:cS.nAlpha
                            alpha_t = cS.alphaGridV(iAlpha);
                            for iSav = 1:cS.nC
                                sav = coh * cS.savingsFracGridV(iSav);
                                c_t = coh - sav;

                                util_c_term = (1-cS.beta) * c_t^cS.psi_1;
                                
                                EV_inner = 0;
                                for iEps = 1:cS.nShocks
                                    R_port = (1-alpha_t)*cS.R_f + alpha_t*cS.R_shock_V(iEps);
                                    w_next = max(cS.w_min, min(cS.w_max, sav * R_port));
                                    
                                    v_prime = v_interp(w_next);
                                    EV_inner = EV_inner + cS.shockProbs(iEps) * v_prime^(1-cS.gamma);
                                end
                                
                                certainty_equiv_term = EV_inner^((1-1/cS.psi)/(1-cS.gamma));
                                current_v = (util_c_term + cS.beta * cS.pi_pathV(t+1) * certainty_equiv_term)^cS.psi_2;
                                
                                if real(current_v) > max_v
                                    max_v = real(current_v);
                                    best_c = c_t;
                                    best_alpha = alpha_t;
                                end
                            end
                        end
                        val_vec(iW) = max_v;
                        c_pol_vec(iW) = best_c;
                        alpha_pol_vec(iW) = best_alpha;
                    end
                    temp_val_t(:, iF) = val_vec;
                    temp_c_pol_t(:, iF) = c_pol_vec;
                    temp_alpha_pol_t(:, iF) = alpha_pol_vec;
                end
                val_retire_all_f(:, :, t - cS.tr + 1) = temp_val_t;
                polS.c{t} = temp_c_pol_t;
                polS.alpha{t} = temp_alpha_pol_t;
                polS.q{t} = zeros(cS.nW, cS.nF);
            end
            valS(:, :, cS.tr:cS.tn-1) = val_retire_all_f(:,:,1:end-1);


            % --- 4. 工作期 (t = K-1, ..., 1) ---
            fprintf('   Solving working period (t = %d to %d)...\n', cS.tr-1, cS.tb);
            for t = cS.tr-1 : -1 : cS.tb
                fprintf('   t = %d\n', t);
                V_t_plus_1_interp = griddedInterpolant({cS.wGridV, cS.fGridV}, valS(:, :, t+1), 'linear', 'linear');
                
                % 临时存储当前期的 V 和 policies
                val_t = zeros(cS.nW, cS.nF);
                c_pol_t = zeros(cS.nW, cS.nF, cS.nShocks);
                alpha_pol_t = zeros(cS.nW, cS.nF, cS.nShocks);
                q_pol_t = zeros(cS.nW, cS.nF, cS.nShocks);

                parfor iF = 1:cS.nF
                    f_t = cS.fGridV(iF);
                    
                    val_w_vec = zeros(cS.nW, 1);
                    c_pol_w_slice = zeros(cS.nW, cS.nShocks);
                    alpha_pol_w_slice = zeros(cS.nW, cS.nShocks);
                    q_pol_w_slice = zeros(cS.nW, cS.nShocks);

                    for iW = 1:cS.nW
                        w_t = cS.wGridV(iW);
                        
                        expected_max_v = 0;
                        
                        % 循环遍历当期冲击 u_t
                        for iU = 1:cS.nShocks
                            u_t_val = cS.uNodes(iU);
                            y_hat = exp(u_t_val);
                            
                            max_v_for_this_u = -inf;
                            best_c_inner = nan; best_alpha_inner = 0; best_q_rate_inner = 0;
                            
                            % 对给定的 u_t, 寻找最优决策
                            for iQ = 1:cS.nQ
                                q_rate_t = cS.qGridV(iQ);
                                q_t = q_rate_t * y_hat;
                                
                                coh = w_t + (1 - cS.tau_y)*(y_hat - q_t);
                                if coh <= 0, continue; end
                                
                                for iAlpha = 1:cS.nAlpha
                                    alpha_t = cS.alphaGridV(iAlpha);
                                    for iSav = 1:cS.nC
                                        sav = coh * cS.savingsFracGridV(iSav);
                                        c_t = coh - sav;
                                        if c_t <= 1e-6, continue; end
                                        
                                        util_c_term = (1-cS.beta) * c_t^cS.psi_1;
                                        
                                        % 计算期望未来价值 (Continuation Value)
                                        EV_inner = 0;
                                        for iZ = 1:cS.nShocks
                                            z_next = cS.zNodes(iZ);
                                            total_growth = cS.G_pathV(t) * exp(z_next);
                                            
                                            f_next_hat = (f_t + q_t) * cS.R_p / total_growth;
                                            f_next_hat = max(cS.f_min, min(cS.f_max, f_next_hat));
                                            
                                            EV_cond_on_z = 0;
                                            for iEps = 1:cS.nShocks
                                                R_port = (1-alpha_t)*cS.R_f + alpha_t*cS.R_shock_V(iEps);
                                                w_next_hat = sav * R_port / total_growth;
                                                w_next_hat = max(cS.w_min, min(cS.w_max, w_next_hat));
                                                
                                                v_prime_hat = V_t_plus_1_interp(w_next_hat, f_next_hat);
                                                v_prime_denorm = total_growth * v_prime_hat;
                                                EV_cond_on_z = EV_cond_on_z + cS.shockProbs(iEps) * v_prime_denorm^(1-cS.gamma);
                                            end
                                            EV_inner = EV_inner + cS.shockProbs(iZ) * EV_cond_on_z;
                                        end
                                        certainty_equiv_term = EV_inner^((1-1/cS.psi)/(1-cS.gamma));
                                        current_v = (util_c_term + cS.beta * cS.pi_pathV(t+1) * certainty_equiv_term)^cS.psi_2;
                                        
                                        if real(current_v) > max_v_for_this_u
                                            max_v_for_this_u = real(current_v);
                                            best_c_inner = c_t;
                                            best_alpha_inner = alpha_t;
                                            best_q_rate_inner = q_rate_t;
                                        end
                                    end
                                end
                            end
                            % 对 u_t 求期望
                            expected_max_v = expected_max_v + cS.shockProbs(iU) * max_v_for_this_u;
                            
                            % 存储依赖于 u_t 的策略
                            c_pol_w_slice(iW, iU) = best_c_inner;
                            alpha_pol_w_slice(iW, iU) = best_alpha_inner;
                            q_pol_w_slice(iW, iU) = best_q_rate_inner;
                        end
                        val_w_vec(iW) = expected_max_v;
                    end
                    % 收集 parfor 的结果
                    val_t(:, iF) = val_w_vec;
                    c_pol_t(:, iF, :) = c_pol_w_slice;
                    alpha_pol_t(:, iF, :) = alpha_pol_w_slice;
                    q_pol_t(:, iF, :) = q_pol_w_slice;
                end
                
                % 将当期结果存入全局变量
                valS(:, :, t) = val_t;
                polS.c{t} = c_pol_t;
                polS.alpha{t} = alpha_pol_t;
                polS.q{t} = q_pol_t;
            end
        end

        function [polS, valS] = value_function_iteration_matrix(cS)
            % =================================================================
            % == 功能: 逆向归纳求解模型 (V3.1 - 修正插值维度匹配问题)
            % == 核心:
            % ==  1. 采用 parfor + 矩阵化 的混合策略。
            % ==  2. [修正] 使用 repmat 扩展 f_next_hat_mat 的维度，使其
            % ==     与 w_next_hat_mat 的维度(6D)完全匹配，解决插值错误。
            % =================================================================

            % --- 1. 初始化 ---
            valS = -inf(cS.nW, cS.nF, cS.tn);
            polS.c = cell(cS.tn, 1);
            polS.alpha = cell(cS.tn, 1);
            polS.q = cell(cS.tn, 1);

            % --- 2. 终期 T ---
            [w_grid_T, f_grid_T] = ndgrid(cS.wGridV, cS.fGridV);
            p_bar_T = (1 - cS.tau_q) * f_grid_T / (cS.tn - cS.tr + 1);
            final_c = w_grid_T + cS.Y_bar_frac + p_bar_T;
            valS(1:cS.nW, 1:cS.nF, cS.tn) = final_c; % V_T = C_T
            
            polS.c{cS.tn} = final_c;
            polS.alpha{cS.tn} = zeros(cS.nW, cS.nF);
            polS.q{cS.tn} = zeros(cS.nW, cS.nF);

            % --- 3. 退休期 (t = T-1, ..., K) ---
            fprintf('   Solving retirement period (t = %d to %d)...\n', cS.tn-1, cS.tr);
            val_retire_all_f = zeros(cS.nW, cS.nF, cS.tn - cS.tr + 1);
            val_retire_all_f(:, :, end) = valS(1:cS.nW, 1:cS.nF, cS.tn);

            for t = cS.tn-1 : -1 : cS.tr
                fprintf('   t = %d\n', t);
                V_t_plus_1_interp_handles = cell(cS.nF, 1);
                for iF = 1:cS.nF
                    V_t_plus_1_interp_handles{iF} = griddedInterpolant(cS.wGridV, val_retire_all_f(:, iF, t - cS.tr + 2), 'linear', 'linear');
                end
                
                temp_val_t = zeros(cS.nW, cS.nF);
                temp_c_pol_t = zeros(cS.nW, cS.nF);
                temp_alpha_pol_t = zeros(cS.nW, cS.nF);

                parfor iF = 1:cS.nF
                    f_k = cS.fGridV(iF);
                    p_bar = (1 - cS.tau_q) * f_k / (cS.tn - cS.tr + 1);
                    v_interp = V_t_plus_1_interp_handles{iF};
                    
                    val_vec = zeros(cS.nW, 1);
                    c_pol_vec = zeros(cS.nW, 1);
                    alpha_pol_vec = zeros(cS.nW, 1);

                    for iW = 1:cS.nW
                        w_t = cS.wGridV(iW);
                        coh = w_t + cS.Y_bar_frac + p_bar;
                        
                        max_v = -inf;
                        best_c = coh; best_alpha = 0;

                        for iAlpha = 1:cS.nAlpha
                            alpha_t = cS.alphaGridV(iAlpha);
                            for iSav = 1:cS.nC
                                sav = coh * cS.savingsFracGridV(iSav);
                                c_t = coh - sav;

                                util_c_term = (1-cS.beta) * c_t^cS.psi_1;
                                
                                EV_inner = 0;
                                for iEps = 1:cS.nShocks
                                    R_port = (1-alpha_t)*cS.R_f + alpha_t*cS.R_shock_V(iEps);
                                    w_next = max(cS.w_min, min(cS.w_max, sav * R_port));
                                    
                                    v_prime = v_interp(w_next);
                                    EV_inner = EV_inner + cS.shockProbs(iEps) * v_prime^(1-cS.gamma);
                                end
                                
                                certainty_equiv_term = EV_inner^((1-1/cS.psi)/(1-cS.gamma));
                                current_v = (util_c_term + cS.beta * cS.pi_pathV(t+1) * certainty_equiv_term)^cS.psi_2;
                                
                                if real(current_v) > max_v
                                    max_v = real(current_v);
                                    best_c = c_t;
                                    best_alpha = alpha_t;
                                end
                            end
                        end
                        val_vec(iW) = max_v;
                        c_pol_vec(iW) = best_c;
                        alpha_pol_vec(iW) = best_alpha;
                    end
                    temp_val_t(:, iF) = val_vec;
                    temp_c_pol_t(:, iF) = c_pol_vec;
                    temp_alpha_pol_t(:, iF) = alpha_pol_vec;
                end
                val_retire_all_f(:, :, t - cS.tr + 1) = temp_val_t;
                polS.c{t} = temp_c_pol_t;
                polS.alpha{t} = temp_alpha_pol_t;
                polS.q{t} = zeros(cS.nW, cS.nF);
            end
            valS(1:cS.nW, 1:cS.nF, cS.tr:cS.tn-1) = val_retire_all_f(:,:,1:end-1);


            % --- 4. 工作期 (t = K-1, ..., 1) ---
            fprintf('   Solving working period (t = %d to %d)...\n', cS.tr-1, cS.tb);
            for t = cS.tr-1 : -1 : cS.tb
                fprintf('   t = %d\n', t);
                V_t_plus_1_interp = griddedInterpolant({cS.wGridV, cS.fGridV}, valS(1:cS.nW, 1:cS.nF, t+1), 'linear', 'linear');
                
                val_t = zeros(cS.nW, cS.nF);
                c_pol_t = zeros(cS.nW, cS.nF, cS.nShocks);
                alpha_pol_t = zeros(cS.nW, cS.nF, cS.nShocks);
                q_pol_t = zeros(cS.nW, cS.nF, cS.nShocks);

                parfor iF = 1:cS.nF
                    f_t = cS.fGridV(iF);
                    
                    val_w_vec = zeros(cS.nW, 1);
                    c_pol_w_slice = zeros(cS.nW, cS.nShocks);
                    alpha_pol_w_slice = zeros(cS.nW, cS.nShocks);
                    q_pol_w_slice = zeros(cS.nW, cS.nShocks);

                    % 循环遍历当期冲击 u_t, 在此循环内部进行矩阵化
                    for iU = 1:cS.nShocks
                        u_t_val = cS.uNodes(iU);
                        y_hat = exp(u_t_val);
                        
                        % --- 矩阵化开始 ---
                        % 1. 创建 W_t 和决策变量的网格
                        % Dims: [nW, nQ, nAlpha, nC]
                        [w_grid, q_grid, alpha_grid, sav_frac_grid] = ndgrid(cS.wGridV, cS.qGridV, cS.alphaGridV, cS.savingsFracGridV);
                        
                        % 2. 计算 coh, c_t, sav
                        q_t_mat = q_grid .* y_hat;
                        coh_mat = w_grid + (1 - cS.tau_y)*(y_hat - q_t_mat);
                        sav_mat = coh_mat .* sav_frac_grid;
                        c_t_mat = coh_mat - sav_mat;
                        
                        valid_choices = (c_t_mat > 1e-6) & (sav_mat >= 0);
                        
                        % 3. 计算期望未来价值 (Continuation Value)
                        prob_z = reshape(cS.shockProbs, [1, 1, 1, 1, cS.nShocks]);
                        prob_eps = reshape(cS.shockProbs, [1, 1, 1, 1, 1, cS.nShocks]);
                        
                        z_nodes_reshaped = reshape(cS.zNodes, [1, 1, 1, 1, cS.nShocks]);
                        eps_nodes_reshaped = reshape(cS.epsNodes, [1, 1, 1, 1, 1, cS.nShocks]);

                        total_growth = cS.G_pathV(t) .* exp(z_nodes_reshaped); % Dims: [1,1,1,1,nZ]
                        
                        f_next_hat_mat = (f_t + q_t_mat) .* cS.R_p ./ total_growth; % Dims: [nW,nQ,nAlpha,nC,nZ]
                        f_next_hat_mat = max(cS.f_min, min(cS.f_max, f_next_hat_mat));

                        R_port_mat = (1-alpha_grid).*cS.R_f + alpha_grid.*(cS.R_f + cS.mu + eps_nodes_reshaped); % Dims: [nW,nQ,nAlpha,nC,1,nEps]
                        w_next_hat_mat = sav_mat ./ total_growth .* R_port_mat; % Dims: [nW,nQ,nAlpha,nC,nZ,nEps]
                        w_next_hat_mat = max(cS.w_min, min(cS.w_max, w_next_hat_mat));
                        
                        % [!!! 核心修正 !!!] 扩展f_next_hat_mat的维度以匹配w_next_hat_mat
                        % f_next_hat不依赖于风险资产冲击(eps),因此我们将其复制扩展到第6个维度
                        f_next_hat_mat = repmat(f_next_hat_mat, 1, 1, 1, 1, 1, cS.nShocks);

                        v_prime_hat = V_t_plus_1_interp(w_next_hat_mat, f_next_hat_mat);
                        v_prime_denorm = total_growth .* v_prime_hat; % Dims: [nW,nQ,nAlpha,nC,nZ,nEps]
                        
                        % 4. 计算期望
                        EV_cond_on_z = sum(prob_eps .* v_prime_denorm.^(1-cS.gamma), 6); % Sum over dim 6 (nEps)
                        EV_inner = sum(prob_z .* EV_cond_on_z, 5); % Sum over dim 5 (nZ)

                        % 5. 计算当前总效用
                        util_c_term = (1-cS.beta) * c_t_mat.^cS.psi_1;
                        certainty_equiv_term = EV_inner.^((1-1/cS.psi)/(1-cS.gamma));
                        
                        current_v_mat = (util_c_term + cS.beta * cS.pi_pathV(t+1) * certainty_equiv_term).^cS.psi_2;
                        current_v_mat(~valid_choices) = -inf;
                        
                        % 6. 寻找最优决策
                        [max_v_for_w, max_idx_flat] = max(reshape(current_v_mat, cS.nW, []), [], 2);
                        
                        val_w_vec = val_w_vec + cS.shockProbs(iU) .* real(max_v_for_w);

                        [~, q_idx, alpha_idx, c_sav_idx] = ind2sub(size(q_grid(1,:,:,:)), max_idx_flat);
                        
                        best_q_rate = cS.qGridV(q_idx);
                        best_alpha = cS.alphaGridV(alpha_idx);
                        
                        best_q_abs = best_q_rate * y_hat;
                        best_coh = cS.wGridV + (1 - cS.tau_y)*(y_hat - best_q_abs);
                        best_sav = best_coh .* cS.savingsFracGridV(c_sav_idx);
                        best_c = best_coh - best_sav;

                        c_pol_w_slice(:, iU) = best_c;
                        alpha_pol_w_slice(:, iU) = best_alpha;
                        q_pol_w_slice(:, iU) = best_q_rate;
                    end
                    val_t(:, iF) = val_w_vec;
                    c_pol_t(:, iF, :) = c_pol_w_slice;
                    alpha_pol_t(:, iF, :) = alpha_pol_w_slice;
                    q_pol_t(:, iF, :) = q_pol_w_slice;
                end
                
                valS(1:cS.nW, 1:cS.nF, t) = val_t;
                polS.c{t} = c_pol_t;
                polS.alpha{t} = alpha_pol_t;
                polS.q{t} = q_pol_t;
            end
        end
        
        function [polS, valS] = value_function_iteration_matrix_nopps(cS)
            % =================================================================
            % == 功能: [Gomes复刻版 V2.1] 逆向归纳求解(无养老金账户)
            % == 核心:
            % ==  1. [修正] 遵循Gomes(2020)代码，将生存概率pi_pathV移入期望算子内部。
            % ==     这是与原文结果匹配的关键一步。
            % ==  2. 退休期由于问题简单，保留高效的 parfor 结构。
            % ==  3. 此版本是原版逻辑在无F状态下的直接、高效的体现。
            % =================================================================

            % --- 1. 初始化 (无F维度) ---
            valS = -inf(cS.nW, cS.tn);
            polS.c = cell(cS.tn, 1);
            polS.alpha = cell(cS.tn, 1);

            % --- 2. 终期 T ---
            valS(:, cS.tn) = (cS.wGridV+ cS.Y_bar_frac)*(1-cS.beta)^cS.psi_2; %  
            polS.c{cS.tn} = cS.wGridV + cS.Y_bar_frac;
            polS.alpha{cS.tn} = zeros(cS.nW, 1);

            % --- 3. 退休期 (t = T-1, ..., K) ---
            fprintf('   Solving retirement period (t = %d to %d) [No PPS]...\n', cS.tn-1, cS.tr);
            for t = cS.tn-1 : -1 : cS.tr
                fprintf('   t = %d\n', t);
                V_t_plus_1_interp = griddedInterpolant(cS.wGridV, valS(:, t+1), 'linear', 'linear');
                
                temp_val_t = zeros(cS.nW, 1);
                temp_c_pol_t = zeros(cS.nW, 1);
                temp_alpha_pol_t = zeros(cS.nW, 1);

                % 退休期问题简单，parfor over W 是最高效的方式
                for iW = 1:cS.nW
                    w_t = cS.wGridV(iW);
                    coh = w_t + cS.Y_bar_frac;
                    
                    max_v = -inf;
                    best_c = coh; best_alpha = 0;

                    for iAlpha = 1:cS.nAlpha
                        alpha_t = cS.alphaGridV(iAlpha);
                        for iSav = 1:cS.nC
                            sav = coh * cS.savingsFracGridV(iSav);
                            c_t = coh - sav;

                            util_c_term = (1-cS.beta) * c_t^cS.psi_1;
                            
                            EV_inner = 0;
                            for iEps = 1:cS.nShocks
                                R_port = (1-alpha_t)*cS.R_f + alpha_t*cS.R_shock_V(iEps);
                                w_next = max(cS.w_min, min(cS.w_max, sav * R_port));
                                v_prime = V_t_plus_1_interp(w_next);
                                % [!!! GOMES REPLICATION MODIFICATION !!!]
                                % 将生存概率 cS.pi_pathV(t+1) 移入期望内部，与Gomes代码一致
                                EV_inner = EV_inner + cS.shockProbs(iEps) * cS.pi_pathV(t+1) * v_prime^(1-cS.gamma);
                            end
                            
                            certainty_equiv_term = EV_inner^((1-1/cS.psi)/(1-cS.gamma));
                            % [!!! GOMES REPLICATION MODIFICATION !!!]
                            % beta 不再与 pi 相乘
                            current_v = (util_c_term + cS.beta * certainty_equiv_term)^cS.psi_2;
                            
                            if real(current_v) > max_v
                                max_v = real(current_v);
                                best_c = c_t;
                                best_alpha = alpha_t;
                            end
                        end
                    end
                    temp_val_t(iW) = max_v;
                    temp_c_pol_t(iW) = best_c;
                    temp_alpha_pol_t(iW) = best_alpha;
                end
                valS(:, t) = temp_val_t;
                polS.c{t} = temp_c_pol_t;
                polS.alpha{t} = temp_alpha_pol_t;
            end

            % --- 4. 工作期 (t = K-1, ..., 1) ---
            fprintf('   Solving working period (t = %d to %d) [No PPS, Matrixized]...\n', cS.tr-1, cS.tb);
            for t = cS.tr-1 : -1 : cS.tb
                fprintf('   t = %d\n', t);
                V_t_plus_1_interp = griddedInterpolant(cS.wGridV, valS(:, t+1), 'linear', 'linear');
                
                % --- 完全矩阵化开始 ---
                % 1. 创建 W_t, 冲击 u_t, 和决策变量的网格
                % Dims: [nW, nAlpha, nC, nU]
                [w_grid, alpha_grid, sav_frac_grid, u_grid] = ndgrid(cS.wGridV, cS.alphaGridV, cS.savingsFracGridV, cS.uNodes);

                % 2. 计算 coh, c_t, sav
                y_hat_mat = exp(u_grid);
                coh_mat = w_grid + y_hat_mat;
                sav_mat = coh_mat .* sav_frac_grid;
                c_t_mat = coh_mat - sav_mat;
                
                valid_choices = (c_t_mat > 1e-6) & (sav_mat >= 0);
                
                % 3. 计算期望未来价值 (Continuation Value)
                prob_z = reshape(cS.shockProbs, [1, 1, 1, 1, cS.nShocks]);
                prob_eps = reshape(cS.shockProbs, [1, 1, 1, 1, 1, cS.nShocks]);
                
                z_nodes_reshaped = reshape(cS.zNodes, [1, 1, 1, 1, cS.nShocks]);
                eps_nodes_reshaped = reshape(cS.epsNodes, [1, 1, 1, 1, 1, cS.nShocks]);

                total_growth = cS.G_pathV(t) .* exp(z_nodes_reshaped); % Dims: [1,1,1,1,nZ]
                
                R_port_mat = (1-alpha_grid).*cS.R_f + alpha_grid.*(cS.R_f + cS.mu + eps_nodes_reshaped); % Dims: [nW,nAlpha,nC,nU,1,nEps]
                w_next_hat_mat = sav_mat ./ total_growth .* R_port_mat; % Dims: [nW,nAlpha,nC,nU,nZ,nEps]
                w_next_hat_mat = max(cS.w_min, min(cS.w_max, w_next_hat_mat));

                v_prime_hat = V_t_plus_1_interp(w_next_hat_mat);
                v_prime_denorm = total_growth .* v_prime_hat; % Dims: [nW,nAlpha,nC,nU,nZ,nEps]
                
                % 4. 计算期望
                % [!!! GOMES REPLICATION MODIFICATION !!!]
                % 将生存概率 cS.pi_pathV(t+1) 移入期望内部，与Gomes代码一致
                v_prime_denorm_powered = v_prime_denorm.^(1-cS.gamma);
                EV_cond_on_z = sum(prob_eps .* v_prime_denorm_powered, 6); % Sum over dim 6 (nEps)
                EV_inner = sum(prob_z .* cS.pi_pathV(t+1) .* EV_cond_on_z, 5); % Sum over dim 5 (nZ)

                % 5. 计算当前总效用
                util_c_term = (1-cS.beta) * c_t_mat.^cS.psi_1;
                certainty_equiv_term = EV_inner.^((1-1/cS.psi)/(1-cS.gamma));
                
                % [!!! GOMES REPLICATION MODIFICATION !!!]
                % beta 不再与 pi 相乘
                current_v_mat = (util_c_term + cS.beta * certainty_equiv_term).^cS.psi_2;
                current_v_mat(~valid_choices) = -inf;
                
                % 6. 对每个 (W, U) 状态，寻找最优决策 (alpha, sav_frac)
                % Reshape to [nW*nU, nAlpha*nC]
                reshaped_v = reshape(permute(current_v_mat, [1, 4, 2, 3]), cS.nW * cS.nShocks, []);
                [max_v_for_wu, max_idx_flat] = max(reshaped_v, [], 2);
                
                % 7. 计算当期价值函数 (对 u 求期望)
                max_v_mat = reshape(real(max_v_for_wu), cS.nW, cS.nShocks);
                val_t = sum(max_v_mat .* reshape(cS.shockProbs, 1, []), 2);
                valS(:, t) = val_t;

                % 8. 恢复最优策略
                [alpha_idx, c_sav_idx] = ind2sub([cS.nAlpha, cS.nC], max_idx_flat);
                
                alpha_pol_t = cS.alphaGridV(alpha_idx);
                alpha_pol_t = reshape(alpha_pol_t, cS.nW, cS.nShocks);

                sav_frac_pol_t = cS.savingsFracGridV(c_sav_idx);
                sav_frac_pol_t = reshape(sav_frac_pol_t, cS.nW, cS.nShocks);

                [w_grid_pol, u_grid_pol] = ndgrid(cS.wGridV, cS.uNodes);
                y_hat_pol = exp(u_grid_pol);
                coh_pol = w_grid_pol + y_hat_pol;
                c_pol_t = coh_pol .* (1 - sav_frac_pol_t);

                polS.c{t} = c_pol_t;
                polS.alpha{t} = alpha_pol_t;
            end
        end        
        
        
        function simS = simulation_nopps(polS, cS)
            % =================================================================
            % == 功能: [Gomes复刻版] 进行生命周期模拟(无养老金账户)
            % == 核心:
            % ==  1. 完全在归一化空间中模拟，最后再转换为绝对值。
            % ==  2. 移除了所有与F和Q相关的变量和逻辑。
            % ==  3. 策略函数插值简化为只依赖于 W_hat (和 u_t)。
            % =================================================================
            
            nsim = 10000;
            rng(2025);

            % --- 初始化存储矩阵 (无F) ---
            simW_hatM = zeros(cS.tn, nsim);
            simC_hatM = zeros(cS.tn, nsim);
            simP_levelM = ones(cS.tn, nsim);
            simAlpha_M = zeros(cS.tn, nsim);
            
            % --- 生成随机冲击 ---
            cdf = cumsum(cS.shockProbs'); 
            cdf_reshaped = reshape(cdf, 1, 1, cS.nShocks);

            rand_draws_z = rand(cS.tr-1, nsim);
            [~, z_idx] = max(rand_draws_z <= cdf_reshaped, [], 3);
            z_shocks = cS.zNodes(z_idx);

            rand_draws_u = rand(cS.tr-1, nsim);
            [~, u_idx] = max(rand_draws_u <= cdf_reshaped, [], 3);
            u_shocks = cS.uNodes(u_idx);

            rand_draws_eps = rand(cS.tn-1, nsim);
            [~, eps_idx] = max(rand_draws_eps <= cdf_reshaped, [], 3);
            eps_shocks = cS.epsNodes(eps_idx);

            % --- 创建插值对象 (简化) ---
            pol_alpha_interp_work = cell(cS.tr-1, 1);
            pol_c_interp_work = cell(cS.tr-1, 1);
            for t = 1:cS.tr-1
                pol_alpha_interp_work{t} = griddedInterpolant({cS.wGridV, cS.uNodes}, polS.alpha{t}, 'linear', 'nearest');
                pol_c_interp_work{t} = griddedInterpolant({cS.wGridV, cS.uNodes}, polS.c{t}, 'linear', 'nearest');
            end
            
            pol_alpha_interp_ret = cell(cS.tn - cS.tr, 1);
            pol_c_interp_ret = cell(cS.tn - cS.tr, 1);
            for t = cS.tr : cS.tn-1
                pol_alpha_interp_ret{t - cS.tr + 1} = griddedInterpolant(cS.wGridV, polS.alpha{t}, 'linear', 'nearest');
                pol_c_interp_ret{t - cS.tr + 1} = griddedInterpolant(cS.wGridV, polS.c{t}, 'linear', 'nearest');
            end

            % --- 模拟 (在归一化空间中) ---
            for t = 1:cS.tn-1
                w_hat = simW_hatM(t, :);
                w_hat = max(cS.w_min, min(cS.w_max, w_hat));
                
                if t < cS.tr
                    % --- 工作期 ---
                    u_t_val = u_shocks(t,:);
                    y_hat = exp(u_t_val);

                    alpha_t = pol_alpha_interp_work{t}(w_hat, u_t_val);
                    c_hat = pol_c_interp_work{t}(w_hat, u_t_val);
                    
                    coh_hat = w_hat + y_hat;
                    c_hat = min(c_hat, coh_hat * 0.9999);
                    sav_hat = coh_hat - c_hat;
                    
                    R_port = (1-alpha_t).*cS.R_f + alpha_t.*(cS.R_f + cS.mu + eps_shocks(t,:));
                    total_growth_factor = cS.G_pathV(t) .* exp(z_shocks(t,:));
                    
                    simW_hatM(t+1,:) = sav_hat .* R_port ./ total_growth_factor;
                    simP_levelM(t+1,:) = simP_levelM(t,:) .* total_growth_factor;
                    
                else
                    % --- 退休期 ---
                    y_hat_retire = cS.Y_bar_frac;
                    coh_hat = w_hat + y_hat_retire;
                    
                    alpha_t = pol_alpha_interp_ret{t - cS.tr + 1}(w_hat);
                    c_hat = pol_c_interp_ret{t - cS.tr + 1}(w_hat);
                    
                    c_hat = min(c_hat, coh_hat*0.9999);
                    sav_hat = coh_hat - c_hat;
                    
                    R_port = (1-alpha_t).*cS.R_f + alpha_t.*(cS.R_f + cS.mu + eps_shocks(t,:));
                    
                    simW_hatM(t+1,:) = sav_hat .* R_port;
                    simP_levelM(t+1,:) = simP_levelM(t,:);
                end
                simC_hatM(t,:) = c_hat;
                simAlpha_M(t,:) = alpha_t;
            end
            
            % --- 最后一期 T ---
            w_hat_T = simW_hatM(cS.tn, :);
            y_hat_T = cS.Y_bar_frac;
            simC_hatM(cS.tn,:) = w_hat_T + y_hat_T;
            simAlpha_M(cS.tn,:) = 0;
            
            % --- 将归一化路径转换为绝对值 ---
            simW_absM = simW_hatM .* simP_levelM;
            simC_absM = simC_hatM .* simP_levelM;
            
            simY_absM = zeros(cS.tn, nsim);
            u_shocks_full = [u_shocks; zeros(cS.tn - (cS.tr-1), nsim)];
            for t = 1:cS.tr-1
                simY_absM(t,:) = exp(u_shocks_full(t,:)) .* simP_levelM(t,:);
            end
            simY_absM(cS.tr:cS.tn, :) = cS.Y_bar_frac .* simP_levelM(cS.tr:cS.tn, :);
            
            % --- 汇总结果 ---
            simS.mean_W = mean(simW_absM, 2);
            simS.mean_C = mean(simC_absM, 2);
            simS.mean_Y = mean(simY_absM, 2);
            simS.mean_Alpha = mean(simAlpha_M, 2, 'omitnan');
            valid_Y_idx = simS.mean_Y > 1e-6;
            simS.mean_W_to_Y = nan(cS.tn, 1);
            simS.mean_W_to_Y(valid_Y_idx) = simS.mean_W(valid_Y_idx) ./ simS.mean_Y(valid_Y_idx);
            
            % 添加Gomes模型特有的输出
            simS.mean_F = zeros(cS.tn, 1); % 无养老金
            simS.mean_Q_rate = zeros(cS.tn, 1); % 无缴费
        end

        function simS = simulation(polS, cS)
            % =================================================================
            % == 功能: 进行生命周期模拟 (V4.0 - 采用纯归一化模拟)
            % == 核心:
            % ==  1. [根本性修正] 整个模拟循环完全在归一化空间中进行。所有状态
            % ==     变量(W_hat, F_hat)和决策(C_hat, Q_hat)均为归一化值。
            % ==  2. 模拟结束后，才将归一化路径整体转换为绝对值用于绘图。
            % ==  3. 此方法彻底解决了VFI和模拟之间的不一致性，是理论上最正确、
            % ==     代码上最简洁的实现方式。
            % =================================================================
            
            nsim = 10000;
            rng(2025); % for reproducibility

            % --- 初始化存储矩阵 ---
            % _hatM 存储归一化值，_absM 存储最终的绝对值
            simW_hatM = zeros(cS.tn, nsim);
            simF_hatM = zeros(cS.tn, nsim);
            simC_hatM = zeros(cS.tn, nsim);
            simP_levelM = ones(cS.tn, nsim); % P_t的绝对水平
            simAlpha_M = zeros(cS.tn, nsim);
            simQ_rate_M = zeros(cS.tn, nsim);
            
            % --- 生成随机冲击 ---
            cdf = cumsum(cS.shockProbs'); 
            cdf_reshaped = reshape(cdf, 1, 1, cS.nShocks);

            rand_draws_z = rand(cS.tr-1, nsim);
            [~, z_idx] = max(rand_draws_z <= cdf_reshaped, [], 3);
            z_shocks = cS.zNodes(z_idx);

            rand_draws_u = rand(cS.tr-1, nsim);
            [~, u_idx] = max(rand_draws_u <= cdf_reshaped, [], 3);
            u_shocks = cS.uNodes(u_idx);

            rand_draws_eps = rand(cS.tn-1, nsim);
            [~, eps_idx] = max(rand_draws_eps <= cdf_reshaped, [], 3);
            eps_shocks = cS.epsNodes(eps_idx);

            % --- 创建插值对象 ---
            pol_q_interp_work = cell(cS.tr-1, 1);
            pol_alpha_interp_work = cell(cS.tr-1, 1);
            pol_c_interp_work = cell(cS.tr-1, 1);
            for t = 1:cS.tr-1
                pol_q_interp_work{t} = griddedInterpolant({cS.wGridV, cS.fGridV, cS.uNodes}, polS.q{t}, 'linear', 'nearest');
                pol_alpha_interp_work{t} = griddedInterpolant({cS.wGridV, cS.fGridV, cS.uNodes}, polS.alpha{t}, 'linear', 'nearest');
                pol_c_interp_work{t} = griddedInterpolant({cS.wGridV, cS.fGridV, cS.uNodes}, polS.c{t}, 'linear', 'nearest');
            end
            
            pol_alpha_interp_ret = cell(cS.tn - cS.tr, 1);
            pol_c_interp_ret = cell(cS.tn - cS.tr, 1);
            for t = cS.tr : cS.tn-1
                pol_alpha_interp_ret{t - cS.tr + 1} = griddedInterpolant({cS.wGridV, cS.fGridV}, polS.alpha{t}, 'linear', 'nearest');
                pol_c_interp_ret{t - cS.tr + 1} = griddedInterpolant({cS.wGridV, cS.fGridV}, polS.c{t}, 'linear', 'nearest');
            end

            % --- 模拟 (在归一化空间中) ---
            p_bar_hat = zeros(1, nsim);
            for t = 1:cS.tn-1
                w_hat = simW_hatM(t, :);
                f_hat = simF_hatM(t, :);
                
                % 确保状态变量在网格范围内
                w_hat = max(cS.w_min, min(cS.w_max, w_hat));
                f_hat = max(cS.f_min, min(cS.f_max, f_hat));
                
                if t < cS.tr
                    % --- 工作期 ---
                    u_t_val = u_shocks(t,:);
                    y_hat = exp(u_t_val);

                    q_rate_t = pol_q_interp_work{t}(w_hat, f_hat, u_t_val);
                    alpha_t = pol_alpha_interp_work{t}(w_hat, f_hat, u_t_val);
                    c_hat = pol_c_interp_work{t}(w_hat, f_hat, u_t_val);
                    
                    q_hat = q_rate_t .* y_hat;
                    coh_hat = w_hat + (1-cS.tau_y)*(y_hat - q_hat);
                    c_hat = min(c_hat, coh_hat * 0.9999); % 消费约束
                    sav_hat = coh_hat - c_hat;
                    
                    R_port = (1-alpha_t).*cS.R_f + alpha_t.*(cS.R_f + cS.mu + eps_shocks(t,:));
                    total_growth_factor = cS.G_pathV(t) .* exp(z_shocks(t,:));
                    
                    simW_hatM(t+1,:) = sav_hat .* R_port ./ total_growth_factor;
                    simF_hatM(t+1,:) = (f_hat + q_hat) * cS.R_p ./ total_growth_factor;
                    simP_levelM(t+1,:) = simP_levelM(t,:) .* total_growth_factor;
                    
                    simQ_rate_M(t,:) = q_rate_t;
                    
                else
                    % --- 退休期 ---
                    if t == cS.tr
                        p_bar_hat = (1 - cS.tau_q) * f_hat / (cS.tn - cS.tr + 1);
                    end
                    y_hat_retire = cS.Y_bar_frac;
                    coh_hat = w_hat + y_hat_retire + p_bar_hat;
                    
                    alpha_t = pol_alpha_interp_ret{t - cS.tr + 1}(w_hat, f_hat);
                    c_hat = pol_c_interp_ret{t - cS.tr + 1}(w_hat, f_hat);
                    
                    c_hat = min(c_hat, coh_hat*0.9999);
                    sav_hat = coh_hat - c_hat;
                    
                    R_port = (1-alpha_t).*cS.R_f + alpha_t.*(cS.R_f + cS.mu + eps_shocks(t,:));
                    
                    simW_hatM(t+1,:) = sav_hat .* R_port; % 退休后 G=1, z=0
                    simF_hatM(t+1,:) = f_hat;
                    simP_levelM(t+1,:) = simP_levelM(t,:);
                    
                    simQ_rate_M(t,:) = 0;
                end
                simC_hatM(t,:) = c_hat;
                simAlpha_M(t,:) = alpha_t;
            end
            
            % --- 最后一期 T ---
            w_hat_T = simW_hatM(cS.tn, :);
            y_hat_T = cS.Y_bar_frac;
            p_bar_hat_T = p_bar_hat;
            simC_hatM(cS.tn,:) = w_hat_T + y_hat_T + p_bar_hat_T;
            simAlpha_M(cS.tn,:) = 0;
            
            % --- 将归一化路径转换为绝对值 ---
            simW_absM = simW_hatM .* simP_levelM;
            simF_absM = simF_hatM .* simP_levelM;
            simC_absM = simC_hatM .* simP_levelM;
            
            simY_absM = zeros(cS.tn, nsim);
            u_shocks_full = [u_shocks; zeros(cS.tn - (cS.tr-1), nsim)];
            for t = 1:cS.tr-1
                simY_absM(t,:) = exp(u_shocks_full(t,:)) .* simP_levelM(t,:);
            end
            simY_absM(cS.tr:cS.tn, :) = cS.Y_bar_frac .* simP_levelM(cS.tr:cS.tn, :);
            
            % --- 汇总结果 ---
            simS.mean_W = mean(simW_absM, 2);
            simS.mean_F = mean(simF_absM, 2);
            simS.mean_C = mean(simC_absM, 2);
            simS.mean_Y = mean(simY_absM, 2);
            simS.mean_Alpha = mean(simAlpha_M, 2, 'omitnan');
            simS.mean_Q_rate = mean(simQ_rate_M, 2);
            valid_Y_idx = simS.mean_Y > 1e-6;
            simS.mean_W_to_Y = nan(cS.tn, 1);
            simS.mean_W_to_Y(valid_Y_idx) = simS.mean_W(valid_Y_idx) ./ simS.mean_Y(valid_Y_idx);
        end        
        
        function plot_results(simS, cS)
            % =================================================================
            % == 功能: 可视化模拟结果
            % =================================================================
            
            ages = (1:cS.tn) + cS.tb_age - 1;
            fig = figure('Name', 'Lifecycle Simulation Results (Extended Model)', 'Position', [100, 100, 1200, 800]);
            
            subplot(2,2,1);
            plot(ages, simS.mean_W, 'b-', 'LineWidth', 2, 'DisplayName', 'Wealth ($\bar{W}$)');
            hold on;
            plot(ages, simS.mean_C, 'r--', 'LineWidth', 2, 'DisplayName', 'Consumption ($\bar{C}$)');
            plot(ages, simS.mean_F, 'g-.', 'LineWidth', 2, 'DisplayName', 'Pension Fund ($\bar{F}$)');
            plot(ages, simS.mean_Y, 'k:', 'LineWidth', 2, 'DisplayName', 'Income ($\bar{Y}$)');
            xline(cS.tr_age, 'k:', 'Retirement', 'LineWidth', 1.5, 'HandleVisibility', 'off');
            hold off;
            title('Mean Lifecycle Profiles'); xlabel('Age'); ylabel('Value (absolute)');
            legend('Location', 'northwest', 'Interpreter', 'latex'); grid on; xlim([cS.tb_age, cS.td_age]);

            subplot(2,2,2);
            plot(ages, simS.mean_Alpha, 'm-', 'LineWidth', 2);
            xline(cS.tr_age, 'k:', 'Retirement', 'LineWidth', 1.5);
            title('Mean Risky Asset Allocation'); xlabel('Age'); ylabel('Risky Share ($\bar{\alpha}$)', 'Interpreter', 'latex');
            ylim([-0.05, 1.05]); grid on; xlim([cS.tb_age, cS.td_age]);

            subplot(2,2,3);
            plot(ages, simS.mean_Q_rate, 'k-', 'LineWidth', 2);
            xline(cS.tr_age, 'k:', 'Retirement', 'LineWidth', 1.5);
            title('Mean Pension Contribution Rate'); xlabel('Age'); ylabel('Contribution Rate ($\bar{q}$)', 'Interpreter', 'latex');
            ylim([-0.01, cS.Q_max + 0.01]); grid on; xlim([cS.tb_age, cS.td_age]);

            subplot(2,2,4);
            plot(ages, simS.mean_W_to_Y, 'c-', 'LineWidth', 2);
            xline(cS.tr_age, 'k:', 'Retirement', 'LineWidth', 1.5);
            title('Mean Wealth-to-Income Ratio'); xlabel('Age'); ylabel('$\bar{W} / \bar{Y}$', 'Interpreter', 'latex');
            grid on; xlim([cS.tb_age, cS.td_age]);
            
            sgtitle('Lifecycle Simulation Results (VFI method)', 'FontSize', 16, 'FontWeight', 'bold');
            
            output_dir = 'fig';
            if ~exist(output_dir, 'dir'), mkdir(output_dir); end
            output_filename = fullfile(output_dir, 'Lifecycle_Simulation_Extended_Model.png');
            print(fig, output_filename, '-dpng', '-r300');
        end
        
        % --- 辅助函数 ---
        function [y_grid_out, trProbM_out] = tauchen(N, rho, sigma, mu, m)
            % =================================================================
            % == 功能: Tauchen (1986) 方法离散化 AR(1) 过程
            % == y_t = mu + rho * y_{t-1} + eps_t, eps_t ~ N(0, sigma^2)
            % =================================================================
            if N == 1, y_grid_out = 0; trProbM_out = 1; return; end
            
            % 对于 i.i.d. 过程, rho=0, 过程的标准差就是冲击的标准差 sigma
            if rho < 1
                std_y = sqrt(sigma^2 / (1-rho^2));
            else % 随机游走过程
                std_y = sigma; 
            end

            y_max = m*std_y; y_min = -y_max;
            y_grid_out = linspace(y_min, y_max, N)';
            d = y_grid_out(2)-y_grid_out(1);
            
            trProbM_out = zeros(N,N);
            for j=1:N % today's state
                for k=1:N % tomorrow's state
                    m_k = (1-rho)*mu + rho*y_grid_out(j); % E[y_t+1 | y_t]
                    if k==1
                        trProbM_out(j,k) = normcdf((y_grid_out(1)-m_k+d/2)/sigma);
                    elseif k==N
                        trProbM_out(j,k) = 1 - normcdf((y_grid_out(N)-m_k-d/2)/sigma);
                    else
                        trProbM_out(j,k) = normcdf((y_grid_out(k)-m_k+d/2)/sigma) - normcdf((y_grid_out(k)-m_k-d/2)/sigma);
                    end
                end
            end
        end

    end
end