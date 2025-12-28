% --- firm.m ---
classdef firm
    methods (Static)
        function M_prices = get_prices_at_t(K_p_hat, K_g_hat, L_hat, cS)
            % 从宏观总量计算要素价格
            if K_p_hat <= 0, K_p_hat = 1e-8; end
            if L_hat <= 0, L_hat = 1e-8; end
            if K_g_hat <= 0, K_g_hat = 1e-8; end
            
            labor_exponent = 1 - cS.alpha - cS.gamma;
            Y_hat_t = (K_p_hat.^cS.alpha) .* (K_g_hat.^cS.gamma) .* (L_hat.^labor_exponent);
            
            MPK_p_period = cS.alpha .* Y_hat_t ./ K_p_hat;
            w_hat_t = labor_exponent .* Y_hat_t ./ L_hat;
            r_mkt_t = MPK_p_period - cS.ddk;
            
            M_prices = struct('Y_hat_t', Y_hat_t, 'w_hat_t', w_hat_t, 'r_mkt_t', r_mkt_t);
        end
    
         function M_prices = get_prices_from_KL(K_p_hat, L_hat, cS, params_ext)
            % =========================================================================
            % == 函数: get_prices_from_KL (新增)
            % == 版本: [v1.0]
            % ==
            % == 目的:
            % ==   - 一个专门的函数，用于根据给定的【资本】和【劳动】总量，
            % ==     计算所有均衡价格和产出。
            % ==   - 这对于求解资本市场不出清的初始时期尤其有用，因为它清晰地
            % ==     体现了“需求决定价格”的逻辑。
            % =========================================================================
            if K_p_hat <= 0, K_p_hat = 1e-8; end
            if L_hat <= 0, L_hat = 1e-8; end

            % 在BGP稳态下，公共资本K_g与产出Y保持一个固定比例
            g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
            n_period = (1 + cS.n_ss)^cS.time_Step - 1;
            g_total_period = (1 + g_A_period) * (1 + n_period) - 1;

            Kg_Y_ratio_factor = cS.I_g_to_Y_ratio_ss / (g_total_period + cS.ddk_g);

            % 从生产函数中反解出Y
            % Y = (K_p^alpha * L^(1-alpha-gamma) * (K_g/Y)^gamma)^ (1/(1-gamma))
            Y_hat_calc = (K_p_hat^cS.alpha * L_hat^(1 - cS.alpha - cS.gamma) * Kg_Y_ratio_factor^cS.gamma)^(1 / (1 - cS.gamma));
            
            % 根据Y计算K_g
            K_g_hat_calc = Kg_Y_ratio_factor * Y_hat_calc;

            % 调用现有函数计算价格
            M_prices = firm.get_prices_at_t(K_p_hat, K_g_hat_calc, L_hat, cS);
        end
    
    end
end