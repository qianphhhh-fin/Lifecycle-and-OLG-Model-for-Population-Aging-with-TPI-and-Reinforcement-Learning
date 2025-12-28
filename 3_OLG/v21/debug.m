% --- 诊断脚本 v4.0：资本会计恒等式检验 ---

fprintf('\n\n==================== [TPI iter=1] Capital Accounting Identity Diagnosis ====================\n');

% --- 1. 检查 t=1 时的资本存量分割 ---
K_p_hat_t1 = aggr_supply.K_p_hat_total_path(1);
K_payg_hat_guess_t1 = K_payg_hat_total_path_guess(1);
K_firm_hat_implied_t1 = K_firm_hat_path(1); % 这是函数前面计算使用的值

fprintf('\n--- [Time t=1] Initial State Breakdown ---\n');
fprintf('%-25s : %.4f\n', 'Household Asset (Kp)', K_p_hat_t1);
fprintf('%-25s : %.4f\n', 'PAYG Fund Guess (Kpayg_g)', K_payg_hat_guess_t1);
fprintf('%-25s : %.4f\n', 'Implied Firm Capital', K_firm_hat_implied_t1);
fprintf('%-25s : %.4e\n', 'Accounting Check (Kp-Kpayg-Kfirm)', K_p_hat_t1 - K_payg_hat_guess_t1 - K_firm_hat_implied_t1);

% --- 2. 检查 t=2 时的资本存量目标 ---
K_p_hat_t2 = aggr_supply.K_p_hat_total_path(2);
K_payg_hat_target_t2 = K_payg_hat_total_target_path(2);
K_firm_hat_consistent_t2 = K_p_hat_t2 - K_payg_hat_target_t2; % 这是一个【应该】有的值

fprintf('\n--- [Time t=2] Evolved State Targets ---\n');
fprintf('%-25s : %.4f\n', 'Household Asset Target (Kp)', K_p_hat_t2);
fprintf('%-25s : %.4f\n', 'PAYG Fund Target (Kpayg_t)', K_payg_hat_target_t2);
fprintf('%-25s : %.4f\n', 'CONSISTENT Firm Capital', K_firm_hat_consistent_t2);


% --- 3. [核心] 检验利率的不一致性 ---
fprintf('\n--- [CORE] Interest Rate Inconsistency Test ---\n');

% 利率(t=2)是基于什么算出来的？是基于 K_firm(t=2) = Kp(t=2) - K_payg_GUESS(t=2)
K_payg_hat_guess_t2 = K_payg_hat_total_path_guess(2);
K_firm_hat_used_for_r_t2 = K_p_hat_t2 - K_payg_hat_guess_t2;
r_target_t2 = target_paths.r_path(2);
fprintf('Interest rate r(2) was calculated using K_firm = Kp(2) - Kpayg_guess(2)\n');
fprintf('  Kp(2)=%.4f, Kpayg_guess(2)=%.4f  =>  K_firm=%.4f  =>  r = %.4f\n', ...
    K_p_hat_t2, K_payg_hat_guess_t2, K_firm_hat_used_for_r_t2, r_target_t2);

% 一个完全一致的系统，它的利率(t=2)应该基于 CONSISTENT 的 K_firm(t=2)
K_g_hat_t2 = K_g_hat_total_path_guess(2);
L_hat_t2 = aggr_supply.L_path(2);
prices_consistent = firm.get_prices_at_t(K_firm_hat_consistent_t2, K_g_hat_t2, L_hat_t2, cS);
r_consistent_t2 = prices_consistent.r_mkt_t;
fprintf('A consistent r(2) should be based on K_firm = Kp(2) - Kpayg_TARGET(2)\n');
fprintf('  Kp(2)=%.4f, Kpayg_target(2)=%.4f  =>  K_firm=%.4f  =>  r = %.4f\n', ...
    K_p_hat_t2, K_payg_hat_target_t2, K_firm_hat_consistent_t2, r_consistent_t2);

fprintf('\nInterest Rate Gap (Target - Consistent) at t=2: %.4e\n', r_target_t2 - r_consistent_t2);

fprintf('\n================================== End Diagnosis ==================================\n\n');