classdef main_transition_utils_bgp
% =========================================================================
% == 类说明: main_transition_utils_bgp (BGP过渡路径求解工具集 - 修正版)
% ==
% == [v2 修正] 修复了 Z_path 变量作用域错误，并使养老金计算逻辑
% ==          与稳态求解器保持一致。
% =========================================================================
    methods (Static)

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                     主求解器 (Master Solver)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function TransitionResults = solve_transition_path(K_p_path_guess, K_g_path_guess, L_path_guess, ...
                                                          Dist_old, valS_new, polS_new, ...
                                                          Z_path, A_path, theta_path, ...
                                                          cS, paramS, ...
                                                          MAX_ITER, TOL, LAMBDA)
            % 初始化
            K_p_path = K_p_path_guess;
            K_g_path = K_g_path_guess;
            L_path   = L_path_guess; % L_path现在是循环中的核心猜测变量之一
            
            fprintf('--- 启动过渡路径迭代求解器 (BGP版本) ---\n');
            fprintf('    最大迭代次数: %d, 收敛容忍度: %.1e, 松弛因子: %.2f\n\n', MAX_ITER, TOL, LAMBDA);
            
            % --- 主迭代循环 ---
            for it = 1:MAX_ITER
                % --- 步骤1: 获取价格路径 (基于当前的宏观路径猜测) ---
                [w_path, r_path, g_A_period_path] = main_transition_utils_bgp.get_price_path(...
                    K_p_path, K_g_path, L_path, A_path, cS);

                % --- 步骤2: 向后求解家庭问题 (Backward Induction) ---
                fprintf('   迭代 %d/%d: [1/2] 正在向后求解家庭问题...\n', it, MAX_ITER);

                % [修正] 将 L_path (当前猜测) 和 Z_path 一起传入，以正确计算养老金b_t
                polS_trans = main_transition_utils_bgp.solve_hh_problem_backward(...
                    w_path, r_path, g_A_period_path, theta_path, L_path, Z_path, ...
                    valS_new, polS_new, cS, paramS);
                
                % --- 步骤3: 向前模拟经济 (Forward Simulation) ---
                fprintf('   迭代 %d/%d: [2/2] 正在向前模拟经济聚合...\n', it, MAX_ITER);
                [K_p_path_new, K_g_path_new, L_path_new] = main_transition_utils_bgp.simulate_forward(...
                    Dist_old, polS_trans, Z_path, w_path, cS, paramS);

                % --- 步骤4: 检查收敛并更新猜测 ---
                error_Kp = max(abs(K_p_path_new - K_p_path));
                error_Kg = max(abs(K_g_path_new - K_g_path));
                error_L  = max(abs(L_path_new - L_path));
                max_error = max([error_Kp, error_Kg, error_L]);
                
                fprintf('       -> 最大误差: %.6f (Kp: %.6f, Kg: %.6f, L: %.6f)\n', ...
                        max_error, error_Kp, error_Kg, error_L);
                
                if max_error < TOL
                    fprintf('\n--- 收敛成功! ---\n');
                    TransitionResults.converged = true;
                    TransitionResults.iterations = it;
                    TransitionResults.final_error = max_error;
                    TransitionResults.K_p_path = K_p_path_new;
                    TransitionResults.K_g_path = K_g_path_new;
                    TransitionResults.L_path = L_path_new;
                    return;
                end
                
                % 使用松弛因子更新猜测
                K_p_path = (1 - LAMBDA) * K_p_path + LAMBDA * K_p_path_new;
                K_g_path = (1 - LAMBDA) * K_g_path + LAMBDA * K_g_path_new;
                L_path   = (1 - LAMBDA) * L_path   + LAMBDA * L_path_new;
            end
            
            fprintf('\n--- 达到最大迭代次数，求解未能收敛。---\n');
            TransitionResults.converged = false;
            TransitionResults.iterations = MAX_ITER;
            TransitionResults.final_error = max_error;
            TransitionResults.K_p_path = K_p_path; 
            TransitionResults.K_g_path = K_g_path;
            TransitionResults.L_path = L_path;
        end
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                         辅助函数 (Helpers)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [w_path, r_path, g_A_period_path] = get_price_path(K_p_path, K_g_path, L_path, A_path, cS)
            T = cS.T_sim;
            w_path = zeros(T, 1);
            r_path = zeros(T, 1);
            g_A_period_path = zeros(T, 1);

            for t = 1:T
                A_t_normalized = 1.0; 
                prices = main_steady_state_utils_bgp.get_prices_at_t(...
                    K_p_path(t), K_g_path(t), L_path(t), A_t_normalized, cS);
                w_path(t) = prices.w_t;
                r_path(t) = prices.r_mkt_t;
                
                if t < T
                    g_A_period_path(t) = (A_path(t+1) / A_path(t)) - 1;
                else
                    g_A_period_path(t) = (1 + cS.g_A_ss)^cS.time_Step - 1;
                end
            end
        end

        % [修正] 函数签名现在接收 L_path 和 Z_path
        function polS_trans = solve_hh_problem_backward(w_path, r_path, g_A_period_path, theta_path, L_path, Z_path, valS_new, polS_new, cS, paramS)
            T = cS.T_sim;
            valS_trans = -Inf(cS.nk, cS.nkpps, cS.nw_expanded, cS.aD_new, T);
            polS_trans_cell = cell(T, cS.aD_new);
            valS_trans(:,:,:,:,T) = valS_new;
            for ia = 1:cS.aD_new
                polS_trans_cell{T, ia} = polS_new(ia);
            end

            for t = (T-1) : -1 : 1
                M_vfi_t = struct('w_t', w_path(t), 'r_mkt_t', r_path(t));

                % [修正] 养老金计算现在使用当前猜测的 L_path(t) 和外生的 Z_path(:,t)
                % 这与稳态求解器的逻辑完全一致，更加严谨
                L_guess_t = L_path(t); % 使用当前的总劳动供给猜测值
                mass_retirees = sum(Z_path((cS.aR_new+1):end, t));
                M_vfi_t.b_t = (theta_path(t) * M_vfi_t.w_t * L_guess_t) / max(1e-9, mass_retirees);

                cS_t = cS;
                cS_t.g_A_ss = g_A_period_path(t) / cS_t.time_Step; 

                bV_payg_vfi_t = zeros(1, cS.aD_new);
                if cS_t.aR_new < cS_t.aD_new
                    bV_payg_vfi_t((cS_t.aR_new+1):cS_t.aD_new) = M_vfi_t.b_t;
                end

                for a_idx = cS.aD_new : -1 : 1
                    vPrime_kkppse_next = [];
                    if a_idx < cS.aD_new
                        vPrime_kkppse_next = valS_trans(:,:,:,a_idx+1,t+1);
                    end
                    
                    [valS_trans(:,:,:,a_idx,t), polS_trans_cell{t,a_idx}] = ...
                        main_steady_state_utils_bgp.HHSolutionByAge_unified_PARFOR(...
                            a_idx, vPrime_kkppse_next, M_vfi_t, bV_payg_vfi_t(a_idx), paramS, cS_t);
                end
            end
            
            polS_trans = cell2mat(polS_trans_cell);
        end
        

        function [K_p_path_new, K_g_path_new, L_path_new] = simulate_forward(Dist_old, polS_trans, Z_path, w_path, cS, paramS)
            T = cS.T_sim;
            K_p_path_new = zeros(T, 1);
            K_g_path_new = zeros(T, 1); % 公共资本的路径也在这里计算
            L_path_new = zeros(T, 1);
            
            Dist_current = Dist_old;
            

            % --- 修正公共资本演化逻辑 ---
            % 使用从 ss_old (初始稳态)中获取的公共资本作为路径的起点
            % ss_old 需要可以被此函数访问。一个简单的方法是修改主程序，将ss_old传入solve_transition_path，
            % 然后再传入这里。为简化，我们先假设一个合理的起点。
            global ss_old; % 声明使用全局变量 ss_old, 需要在主脚本中定义
            if isempty(ss_old)
                 error('全局变量 ss_old 未在主脚本中定义，无法获取公共资本起点');
            end
            K_g_path_new(1) = ss_old.K_public_hat; 
            K_g_path_new(1) = ss_old.K_public_hat; % 使用旧稳态的公共资本作为起点，ss_old 需要在主程序中作为全局或传入

            % 向前循环
            for t = 1:T
                % 1. 聚合当期劳动供给 L_t
                L_t_agg = 0;
                dist_t_all_ages = Dist_current;
                for ia = 1:cS.aR_new
                    mass_by_epsilon = squeeze(sum(dist_t_all_ages(:,:,:,ia), [1,2]));
                    if ~isempty(mass_by_epsilon)
                       L_t_agg = L_t_agg + sum( (cS.ageEffV_new(ia) * paramS.leGridV') .* mass_by_epsilon' , 'all');
                    end
                end
                L_path_new(t) = L_t_agg;

                if t == T, break; end

                % 2. 演化到下一期 (t+1)
                polS_t = polS_trans(t,:);
                Dist_next = zeros(size(Dist_current));
                
                % 获取下一期 (t+1) 的新生儿分布
                dist_newborn_next = zeros(cS.nk, cS.nkpps, cS.nw_expanded);
                dist_newborn_next(1, 1, 1:cS.nw) = paramS.leProb1V';
                Dist_next(:,:,:,1) = dist_newborn_next * Z_path(1, t+1);
                
                K_p_agg_next = 0;
                
                for ia = 1:(cS.aD_new-1)
                    dist_ia_t = squeeze(dist_t_all_ages(:,:,:,ia));
                    k_prime_idx_ia = squeeze(main_steady_state_utils_bgp.get_policy_index_matrix_unified(polS_t(ia), cS, 'k'));
                    kpps_prime_idx_ia = squeeze(main_steady_state_utils_bgp.get_policy_index_matrix_unified(polS_t(ia), cS, 'kpps'));
                    
                    dist_ia_plus_1_next = zeros(cS.nk, cS.nkpps, cS.nw_expanded);
                    trans_mat_next = paramS.TrProbM_by_age{ia + 1};

                    for ik = 1:cS.nk
                       for ikpps = 1:cS.nkpps
                           for ie = 1:cS.nw_expanded
                               mass = dist_ia_t(ik, ikpps, ie);
                               if mass < 1e-20, continue; end

                               ik_p = k_prime_idx_ia(ik, ikpps, ie);
                               ikpps_p = kpps_prime_idx_ia(ik, ikpps, ie);
                               
                               trans_probs = reshape(trans_mat_next(ie,:), [1, 1, cS.nw_expanded]);
                               dist_ia_plus_1_next(ik_p, ikpps_p, :) = dist_ia_plus_1_next(ik_p, ikpps_p, :) + mass * trans_probs;
                           end
                       end
                    end
                    
                    mass_at_ia = sum(dist_ia_t, 'all');
                    if mass_at_ia > 1e-12
                        rescale = Z_path(ia+1, t+1) / sum(dist_ia_plus_1_next,'all');
                        Dist_next(:,:,:,ia+1) = dist_ia_plus_1_next * rescale;
                    end
                end

                % 3. 聚合 t+1 的资本存量 K_{t+1}
                K_p_path_new(t+1) = sum(cS.kGridV(squeeze(Dist_next(1:cS.nk,:,:,:))),'all');
                
                % 更新分布
                Dist_current = Dist_next;
            end
            % 公共资本路径暂时保持线性插值的猜测（简化）
            K_g_path_new = K_g_path_new(1) + (K_g_path_new(end)-K_g_path_new(1))/(T-1) * (0:T-1)';

        end

    end
end