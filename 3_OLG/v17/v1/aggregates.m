% --- aggregates.m ---
classdef aggregates
    methods (Static)
        function [K_private_hat_agg, C_agg, Bequest_generated_agg, Total_Tax_Revenue] = get_aggregates(Dist, polS, cS, paramS)
            % =========================================================================
            % == 函数: get_aggregates
            % == 版本: [v3.0 - PPS 聚合版]
            % ==
            % == 核心修改:
            % ==   - 将PPS资产 kpps_prime 加总到总私人资本 K_private_hat_agg 中。
            % ==   - 将PPS资产 kpps_prime 加总到产生的总遗赠中。
            % =========================================================================

            K_k_hat_agg_raw = 0; 
            K_kpps_hat_agg_raw = 0; % [新增] 追踪PPS资产
            C_agg_raw = 0;
            Bequest_k_generated_raw = 0;
            Bequest_kpps_generated_raw = 0; % [新增]
            Total_Tax_Revenue = 0;
            
            for ia = 1:cS.aD_new
                mass_dist_ia = Dist(:,:,:,ia);
                C_agg_raw = C_agg_raw + sum(polS(ia).c .* mass_dist_ia, 'all');
                Total_Tax_Revenue = Total_Tax_Revenue + sum(polS(ia).tax_regular .* mass_dist_ia, 'all');
                
                % [修改] 聚合两个资产
                K_k_hat_agg_raw = K_k_hat_agg_raw + sum(polS(ia).k_prime .* mass_dist_ia, 'all');
                if isfield(polS(ia), 'kpps_prime')
                    K_kpps_hat_agg_raw = K_kpps_hat_agg_raw + sum(polS(ia).kpps_prime .* mass_dist_ia, 'all');
                end
                
                prob_death = (1 - cS.s_pathV(ia));
                
                % [修改] 遗赠也包括两个账户的资产
                Bequest_k_generated_raw = Bequest_k_generated_raw + sum(polS(ia).k_prime .* mass_dist_ia, 'all') * prob_death;
                if isfield(polS(ia), 'kpps_prime')
                     Bequest_kpps_generated_raw = Bequest_kpps_generated_raw + sum(polS(ia).kpps_prime .* mass_dist_ia, 'all') * prob_death;
                end
            end
            
            n_period = (1 + cS.n_ss)^cS.time_Step - 1;

            % [修改] 总私人资本 = 常规资本 + PPS资本
            K_private_hat_agg_raw = K_k_hat_agg_raw + K_kpps_hat_agg_raw;
            K_private_hat_agg = K_private_hat_agg_raw / (1 + n_period);

            % [修改] 总遗赠 = 两个账户的遗赠之和
            Bequest_generated_agg_raw = Bequest_k_generated_raw + Bequest_kpps_generated_raw;
            Bequest_generated_agg = Bequest_generated_agg_raw / (1 + n_period);
            
            C_agg = C_agg_raw;
        end
    end
end