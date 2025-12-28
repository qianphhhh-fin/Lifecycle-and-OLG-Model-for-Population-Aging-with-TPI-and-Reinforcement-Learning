% --- START OF FILE main_olg_v12_utils_simplified.m ---
classdef main_olg_v12_utils_simplified
    methods (Static)
        function [ss, eq_found] = solve_steady_state_RA(cS)
            % 这是一个代表性消费者(RA)模型的稳态求解器
            % 它通过迭代资本存量 K 来寻找均衡点
            
            fprintf('\n--- [RA版] 开始求解稳态 ---\n');
            K_guess = 3.0; % 初始猜测值
            max_iter = 100; tol = 1e-7; damp = 0.5; eq_found = false;

            fprintf('%5s | %12s | %12s | %12s | %12s\n', 'Iter', 'K_guess', 'r_net(%)', 'K_supply', 'K_error');
            fprintf('%s\n', repmat('-', 65, 1));

            for iter = 1:max_iter
                % 步骤 1: 根据 K_guess 计算所有价格和宏观量
                L_ss = 1.0; 
                if K_guess <= 0, K_guess = 1e-6; end
                Y_guess = cS.A * (K_guess^cS.alpha) * (L_ss^(1-cS.alpha));
                r_rental_guess = cS.alpha * Y_guess / K_guess;
                w_guess = (1-cS.alpha) * Y_guess / L_ss;
                
                % 计算家庭的净回报率
                r_net_guess = (r_rental_guess * (1 - cS.tau_k)) - cS.ddk;
                
                % 步骤 2: 根据欧拉方程反解出家庭的消费
                % u'(C) = beta * (1+r_net) * u'(C) --> 1 = beta * (1+r_net)
                % 我们需要调整K，使得 r_net = (1/beta) - 1
                
                % 基于这个 r_net，计算对应的消费 C
                % 注意：消费 C 是由家庭的预算约束决定的
                T_l_guess = cS.tau_l * w_guess * L_ss;
                T_k_guess = cS.tau_k * r_rental_guess * K_guess;
                
                % 净可支配收入 NDI = (Y - ddk*K) - (T_l + T_k)
                NDI_guess = (Y_guess - cS.ddk * K_guess) - (T_l_guess + T_k_guess);
                
                % C_value = NDI / (1+tau_c)
                C_guess = NDI_guess / (1 + cS.tau_c);
                
                % 步骤 3: 根据市场出清条件，计算出与当前价格匹配的资本供给 K_supply
                % Y = C + I + G
                % I = ddk*K
                % G = T = T_l + T_k + T_c
                T_c_guess = cS.tau_c * C_guess;
                G_guess = T_l_guess + T_k_guess + T_c_guess;
                
                % 从 Y = C + I + G 反解出资本存量 K
                % A*K^alpha = C + ddk*K + G
                % A*K^alpha - ddk*K = C + G
                % 这是一个关于K的非线性方程，我们需要求解它
                f = @(k) cS.A * k.^cS.alpha - cS.ddk * k - (C_guess + G_guess);
                
                % 使用 fzero 求解使商品市场出清的 K
                % 我们需要一个合理的初始点，就用当前的 K_guess
                options = optimset('Display','off');
                [K_supply, ~, exitflag] = fzero(f, K_guess, options);
                if exitflag ~= 1
                   % 如果求解失败，可能是函数无解或范围问题，用一个简单的更新规则
                   K_supply = K_guess; % 保持不变以避免崩溃
                   warning('fzero failed to find a solution for K_supply.');
                end

                K_error = K_supply - K_guess;

                fprintf('%5d | %12.4f | %11.4f | %12.4f | %12.3e\n', iter, K_guess, r_net_guess*100, K_supply, K_error);

                if abs(K_error) < tol
                    fprintf('✅ RA稳态均衡收敛！\n');
                    eq_found = true;
                    break;
                end
                
                % 更新猜测
                K_guess = (1 - damp) * K_guess + damp * K_supply;
            end
            
            if ~eq_found
                warning('RA稳态未收敛'); 
                % 即使未收敛，也返回最后一次迭代的结果用于审计
            end

            % 步骤 4: 打包最终结果
            ss = struct();
            K_ss = K_guess;
            L_ss = 1.0;
            Y_ss = cS.A * (K_ss^cS.alpha) * (L_ss^(1-cS.alpha));
            r_rental_ss = cS.alpha * Y_ss / K_ss;
            w_ss = (1-cS.alpha) * Y_ss / L_ss;
            
            T_labor_ss = cS.tau_l * w_ss * L_ss;
            T_capital_ss = cS.tau_k * r_rental_ss * K_ss;
            NDI_ss = (Y_ss - cS.ddk * K_ss) - (T_labor_ss + T_capital_ss);
            C_value_ss = NDI_ss / (1 + cS.tau_c);
            T_consumption_ss = C_value_ss * cS.tau_c;
            T_total_ss = T_labor_ss + T_capital_ss + T_consumption_ss;
            G_ss = T_total_ss;
            
            ss.K = K_ss;
            ss.L = L_ss;
            ss.Y = Y_ss;
            ss.r_rental = r_rental_ss;
            ss.w = w_ss;
            ss.C_value = C_value_ss;
            ss.I_gross = cS.ddk * K_ss;
            ss.G = G_ss;
            ss.T_labor = T_labor_ss;
            ss.T_capital = T_capital_ss;
            ss.T_consumption = T_consumption_ss;
            ss.T_total = T_total_ss;
        end
    end
end
% --- END OF FILE main_olg_v12_utils_simplified.m ---