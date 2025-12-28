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
    end
end