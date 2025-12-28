function [kHistM_out, kPpsHistM_out, cHistM_out, cppsHistM_out] = HHSimulation_olgm_naive(...
        eIdxM_group, R_k_net, w, TR, bV_payg, paramS_sim, cS_sim)
    % [带调试打印功能的版本]
    % 这是一个使用“幼稚策略”的确定性模拟器，用于与Python进行一致性检验。
    % 幼稚策略:
    % 1. PPS缴费 (c_pps) 永远为 0。
    % 2. 储蓄 (k') = 0.3 * (可用于消费和储蓄的总资源)。

    SAVE_FRAC = 0.30; % 储蓄率

    nSim = size(eIdxM_group, 1);
    aD = cS_sim.aD_new;

    kHistM_out = zeros(nSim, aD);
    kPpsHistM_out = zeros(nSim, aD);
    cHistM_out = zeros(nSim, aD);
    cppsHistM_out = zeros(nSim, aD);

    pps_return_factor = 1 + ((R_k_net - 1) + cS_sim.pps_return_rate_premium);
    k_next = zeros(nSim, 1);
    k_pps_next = zeros(nSim, 1);

    for a_idx = 1:aD
        k_now = k_next;
        k_pps_now = k_pps_next;
        kHistM_out(:, a_idx) = k_now;
        kPpsHistM_out(:, a_idx) = k_pps_now;

        c_pps_decision = zeros(nSim, 1);
        k_prime_decision = zeros(nSim, 1);
        c_decision = zeros(nSim, 1);

        for i_sim = 1:nSim
            % 为调试第一个个体设置一个标志
            if i_sim == 1 && a_idx <= 3
                paramS_sim.is_debugging = true;
            else
                paramS_sim.is_debugging = false;
            end
            
            epsilon_val = paramS_sim.leGridV(eIdxM_group(i_sim, a_idx));
            
            if paramS_sim.is_debugging
                fprintf('\n--- MATLAB DEBUG (i_sim=1, a_idx=%d) ---\n', a_idx);
                fprintf('State In: k_now=%.8f, k_pps_now=%.8f, eps_val=%.8f\n', k_now(i_sim), k_pps_now(i_sim), epsilon_val);
            end

            [resources, ~, non_capital_income] = main_olg_v8_utils.HHIncome_Huggett(...
                k_now(i_sim), R_k_net, w, TR, bV_payg(a_idx), ...
                c_pps_decision(i_sim), a_idx, paramS_sim, cS_sim, epsilon_val);
            
            if paramS_sim.is_debugging
                fprintf('HHIncome Out: resources=%.8f, non_cap_inc=%.8f\n', resources, non_capital_income);
            end

            k_prime_decision(i_sim) = SAVE_FRAC * resources;
            c_expend = resources - k_prime_decision(i_sim);
            c_decision(i_sim) = max(cS_sim.cFloor, c_expend / (1 + cS_sim.tau_c));
            
            if paramS_sim.is_debugging
                fprintf('Decision: k_prime=%.8f, c_decision=%.8f\n', k_prime_decision(i_sim), c_decision(i_sim));
                fprintf('-------------------------------------------\n');
            end
        end
        
        cHistM_out(:, a_idx) = c_decision;
        cppsHistM_out(:, a_idx) = c_pps_decision;

        if a_idx < aD
            kMin_double = double(cS_sim.kMin);
            kMax_double = double(cS_sim.kMax);
            k_next = max(kMin_double, min(kMax_double, k_prime_decision));

            pps_withdrawal = 0;
            if a_idx >= cS_sim.aR_new && cS_sim.pps_active
                pps_withdrawal = k_pps_now * cS_sim.pps_withdrawal_rate;
            end
            
            kppsMin_double = double(cS_sim.kppsMin);
            kppsMax_double = double(cS_sim.kppsMax);
            k_pps_next_unclamped = (k_pps_now + cppsHistM_out(:, a_idx) - pps_withdrawal) * pps_return_factor;
            k_pps_next = max(kppsMin_double, min(kppsMax_double, k_pps_next_unclamped));
        end
    end
end