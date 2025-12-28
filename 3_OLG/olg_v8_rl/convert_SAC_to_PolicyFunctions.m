% --- START OF FILE convert_SAC_to_PolicyFunctions.m ---
function [cPolM_sac, kPolM_sac, cPpsPolM_choice_sac, valM_dummy_sac] = convert_SAC_to_PolicyFunctions(...
    sac_agent, cS, paramS_rl, M_fixed)
%CONVERT_SAC_TO_POLICYFUNCTIONS 将训练好的SAC Agent转换为VFI格式的策略函数
%   sac_agent: 训练好的rlSACAgent对象
%   cS: OLG模型参数
%   paramS_rl: RL相关参数 (包含leGridV)
%   M_fixed: 一个固定的宏观经济状态结构体 (与OLGEnv_v8_SAC.current_M格式相同)

fprintf('Converting SAC agent to policy function matrices for M_fixed...\n');

% 初始化输出矩阵
cPolM_sac  = zeros(cS.nk, cS.nkpps, cS.nw, cS.aD_new);
kPolM_sac  = zeros(cS.nk, cS.nkpps, cS.nw, cS.aD_new);
cPpsPolM_choice_sac = zeros(cS.nk, cS.nkpps, cS.nw, cS.aD_new);
valM_dummy_sac = NaN(cS.nk, cS.nkpps, cS.nw, cS.aD_new); % SAC不直接输出V

% 构建固定的bV_payg (与OLGEnv中逻辑一致)
bV_payg_fixed = zeros(1, cS.aD_new);
if cS.aR_new < cS.aD_new
    if isfield(M_fixed, 'b_payg_avg_retiree')
        bV_payg_fixed(cS.aR_new + 1 : cS.aD_new) = M_fixed.b_payg_avg_retiree;
    elseif isfield(M_fixed, 'b_payg_avg_for_obs') % Fallback if name differs
        bV_payg_fixed(cS.aR_new + 1 : cS.aD_new) = M_fixed.b_payg_avg_for_obs;
    else
        error('M_fixed needs b_payg_avg_retiree or b_payg_avg_for_obs field.');
    end
end
% 确保 M_fixed 中用于观测的 b_payg 字段存在
if ~isfield(M_fixed, 'b_payg_avg_for_obs')
    if isfield(M_fixed, 'b_payg_avg_retiree')
        M_fixed.b_payg_avg_for_obs = M_fixed.b_payg_avg_retiree;
    else
        error('M_fixed needs b_payg_avg_retiree or b_payg_avg_for_obs field.');
    end
end


% 定义观测归一化参数 (与 OLGEnv_v8_SAC 中一致)
% 更好的做法是将env的归一化参数传递过来，或者从env对象获取
obs_norm_min_conv = [cS.kMin, cS.kppsMin, 1, 1, ...
                     M_fixed.R_k_net_factor, M_fixed.w_gross, ... % 使用M_fixed的值作为范围的单点
                     M_fixed.TR_total, M_fixed.b_payg_avg_for_obs, ...
                     M_fixed.tau_l, M_fixed.theta_payg_actual];
obs_norm_max_conv = obs_norm_min_conv; % 因为M是固定的，所以min=max for M parts
obs_norm_max_conv(1:4) = [cS.kMax, cS.kppsMax, cS.aD_new, cS.nw]; % k,k_pps,age,eps有范围
obs_norm_range_conv = obs_norm_max_conv - obs_norm_min_conv;
obs_norm_range_conv(obs_norm_range_conv < 1e-6) = 1;


total_states_to_predict = cS.nk * cS.nkpps * cS.nw * cS.aD_new;
predicted_count = 0;

for a_idx = 1:cS.aD_new
    if mod(a_idx,2)==0, fprintf('  Converting age group %d/%d\n', a_idx, cS.aD_new); end
    for i_k = 1:cS.nk
        for i_kpps = 1:cS.nkpps
            for i_eps = 1:cS.nw
                current_k_val_conv = cS.kGridV(i_k);
                current_k_pps_val_conv = cS.kppsGridV(i_kpps);
                current_eps_idx_conv = i_eps;
                current_epsilon_val_conv = paramS_rl.leGridV(current_eps_idx_conv);

                % 1. 构建观测向量 (与OLGEnv_v8_SAC.getObservation一致)
                raw_obs_vec_conv = [current_k_val_conv, current_k_pps_val_conv, ...
                                    a_idx, current_eps_idx_conv, ...
                                    M_fixed.R_k_net_factor, M_fixed.w_gross, ...
                                    M_fixed.TR_total, M_fixed.b_payg_avg_for_obs, ...
                                    M_fixed.tau_l, M_fixed.theta_payg_actual];
                
                obs_conv_norm = (raw_obs_vec_conv - obs_norm_min_conv) ./ obs_norm_range_conv;
                obs_conv_norm = max(0, min(1, obs_conv_norm));
                obs_for_agent = {obs_conv_norm(:)'};

                % 2. 使用SAC agent的actor获取确定性行动 (比例)
                % 对于SAC，getAction通常会考虑探索。评估时，我们想要确定性行动。
                % 如果actor是确定性的(rlDeterministicActorRepresentation)，则直接用。
                % 如果是随机的，需要取分布的均值。
                % 对于MATLAB的SAC，直接调用getAction在评估模式下（或从actor表示中获取）
                % 通常，agent.UseExplorationPolicy = false; (但这不是标准属性)
                % 或者直接从actor网络获取输出：
                % action_proportions = predict(sac_agent.Actor.Representation.getModel(), obs_for_agent);
                % 简便方法:
                action_proportions_cell = getAction(sac_agent, obs_for_agent); % 确保agent处于评估模式
                action_proportions = action_proportions_cell{1};

                prop_pps_contrib_sac = max(0, min(1, action_proportions(1)));
                prop_non_pps_saving_sac = max(0, min(1, action_proportions(2)));

                % 3. 将比例行动转换回绝对值 (与OLGEnv_v8_SAC.step中的逻辑一致)
                actual_c_pps_sac = 0;
                max_permissible_cpps_sac = 0;
                if a_idx <= cS.aR_new && ...
                   cS.physAgeMap{a_idx}(1) <= cS.pps_contribution_age_max_idx && ...
                   cS.pps_active
                    
                    age_efficiency_sac = cS.ageEffV_new(a_idx);
                    current_gross_labor_income_sac = M_fixed.w_gross * age_efficiency_sac * current_epsilon_val_conv;
                    
                    if current_gross_labor_income_sac > 1e-6
                        max_cpps_by_frac_sac = current_gross_labor_income_sac * cS.pps_max_contrib_frac;
                        max_permissible_cpps_sac = min(cS.pps_annual_contrib_limit, max_cpps_by_frac_sac);
                        max_permissible_cpps_sac = max(0, max_permissible_cpps_sac);
                    end
                    actual_c_pps_sac = prop_pps_contrib_sac * max_permissible_cpps_sac;
                    actual_c_pps_sac = max(0, min(actual_c_pps_sac, max_permissible_cpps_sac));
                end
                cPpsPolM_choice_sac(i_k, i_kpps, i_eps, a_idx) = actual_c_pps_sac;

                % 计算资源
                paramS_hh_conv.tau_l = M_fixed.tau_l;
                paramS_hh_conv.theta_payg_actual_for_hh = M_fixed.theta_payg_actual;
                paramS_hh_conv.pps_tax_deferral_active = cS.pps_active;
                b_payg_this_age_conv = bV_payg_fixed(a_idx);

                [resources_after_pps_sac, ~, ~] = main_olg_v8_utils.HHIncome_Huggett(...
                    current_k_val_conv, M_fixed.R_k_net_factor, M_fixed.w_gross, ...
                    M_fixed.TR_total, b_payg_this_age_conv, actual_c_pps_sac, ...
                    a_idx, paramS_hh_conv, cS, current_epsilon_val_conv);

                % 计算k_prime和c
                consumption_floor_spending_sac = cS.cFloor * (1 + cS.tau_c);
                resources_for_kprime_c_above_floor_sac = resources_after_pps_sac - consumption_floor_spending_sac;
            
                actual_k_prime_sac = 0;
                if resources_for_kprime_c_above_floor_sac >= 0
                    actual_k_prime_sac = prop_non_pps_saving_sac * resources_for_kprime_c_above_floor_sac;
                    actual_k_prime_sac = max(cS.kMin, min(actual_k_prime_sac, resources_for_kprime_c_above_floor_sac));
                else
                    actual_k_prime_sac = cS.kMin;
                end
                actual_k_prime_sac = max(cS.kMin, min(actual_k_prime_sac, cS.kMax));
                kPolM_sac(i_k, i_kpps, i_eps, a_idx) = actual_k_prime_sac;

                consumption_expenditure_sac = resources_after_pps_sac - actual_k_prime_sac;
                current_c_sac = max(cS.cFloor, consumption_expenditure_sac / (1 + cS.tau_c));
                cPolM_sac(i_k, i_kpps, i_eps, a_idx) = current_c_sac;
                
                predicted_count = predicted_count + 1;
            end
        end
    end
end
fprintf('Conversion complete. Predicted %d states.\n', predicted_count);
end
% --- END OF FILE convert_SAC_to_PolicyFunctions.m ---