% =========================================================================
% == 诊断脚本: debug_s_net_behavior.m
% == 目的: 诊断 S_net(K) 函数的行为，以理解 fzero 为何失败，并验证
% ==        最终修正版代码的行为是否符合经济学直觉。
% == 版本: 基于 debug_full_annuity.m 的最终修正版函数。
% =========================================================================
clear; close all;
addpath(pwd);
fprintf('=== 诊断脚本: 检查 S_net(K) 函数的行为 ===\n\n');

%% --- 1. 初始化 (与主脚本完全相同) ---
cS = ParameterValues_AiyagariStyle();
paramS = struct();

ngrid = 50; cS.nk = ngrid;
cS = generateGrids(cS);
[paramS.leLogGridV, paramS.leTrProbM, paramS.leStationaryDist_e] = EarningProcess(cS);
paramS.leGridV = exp(paramS.leLogGridV(:));


%% --- 2. 核心诊断循环 ---
fprintf('\n--- 2. 循环测试不同的 K_guess 值 ---\n');
K_test_vec = linspace(0.1, 20, 20); % 测试一系列K值
S_net_vec = zeros(size(K_test_vec));
K_model_vec = zeros(size(K_test_vec));
C_vec = zeros(size(K_test_vec));
Y_vec = zeros(size(K_test_vec));
I_vec = zeros(size(K_test_vec));
r_mkt_vec = zeros(size(K_test_vec));

for i = 1:length(K_test_vec)
    K_guess = K_test_vec(i);
    fprintf('   测试 K_guess = %.4f ...\n', K_guess);
        [S_net, ss_temp, ~, ~, K_model] = calculate_net_saving(K_guess, cS, paramS);
        S_net_vec(i) = S_net;
        K_model_vec(i) = K_model;
        C_vec(i) = ss_temp.C;
        Y_vec(i) = ss_temp.Y_from_production;
        I_vec(i) = ss_temp.I;
        r_mkt_vec(i) = ss_temp.r_mkt;
        fprintf('      S_net = %.4f, K_model = %.4f, r_mkt = %.4f\n', S_net, K_model, ss_temp.r_mkt);
end

%% --- 3. 结果可视化 ---
fprintf('\n--- 3. 可视化结果 ---\n');

figure('Name', 'S_net vs. K_guess Diagnostic (Final Corrected Version)');

% 图 1: 净储蓄 S_net vs K
subplot(2, 2, 1);
plot(K_test_vec, S_net_vec, '-o', 'LineWidth', 1.5, 'Color', [0, 0.4470, 0.7410]);
hold on;
plot(K_test_vec, zeros(size(K_test_vec)), 'k--');
xlabel('Aggregate Capital Guess (K_{guess})');
ylabel('Net Saving (S_{net})');
title('Net Saving (Y - C - I)');
grid on;
if any(S_net_vec > 0) && any(S_net_vec < 0)
    text(K_test_vec(end)*0.6, max(S_net_vec)*0.8, '✅ 符号变化', 'Color', 'g', 'FontWeight', 'bold');
else
    text(K_test_vec(end)*0.6, mean(S_net_vec), '❌ 无符号变化', 'Color', 'r', 'FontWeight', 'bold');
end


% 图 2: 资本市场出清 K_model vs K_guess
subplot(2, 2, 2);
plot(K_test_vec, K_model_vec, '-o', 'LineWidth', 1.5, 'Color', [0.8500, 0.3250, 0.0980]);
hold on;
plot(K_test_vec, K_test_vec, 'k--', 'DisplayName', '45-degree line');
xlabel('Aggregate Capital Guess (K_{guess})');
ylabel('Model Implied Capital (K_{model})');
title('Capital Market Clearing');
legend('Location','best');
grid on;

% 图 3: 宏观聚合量
subplot(2, 2, 3);
plot(K_test_vec, Y_vec, '-s', 'DisplayName', 'Y (Production)');
hold on;
plot(K_test_vec, C_vec, '-d', 'DisplayName', 'C (Consumption)');
plot(K_test_vec, I_vec, '-^', 'DisplayName', 'I (Investment)');
xlabel('Aggregate Capital Guess (K_{guess})');
ylabel('Macro Aggregates');
title('Y, C, I');
legend('Location','best');
grid on;

% 图 4: 市场利率
subplot(2, 2, 4);
plot(K_test_vec, r_mkt_vec*100, '-o', 'LineWidth', 1.5, 'Color', [0.4660, 0.6740, 0.1880]);
xlabel('Aggregate Capital Guess (K_{guess})');
ylabel('Market Interest Rate (%)');
title('Interest Rate');
grid on;

sgtitle('Diagnostic Plot for `calculate_net_saving` (Final Corrected Version)');

%% ========================================================================
%  == 辅助函数、求解器和验证函数
%  ========================================================================

% -------------------------------------------------------------------------
% --- 1. 参数和网格设定 ---
% -------------------------------------------------------------------------
function cS = ParameterValues_AiyagariStyle()
    fprintf('   设定 Aiyagari-Huggett 模型参数...\n');
    cS.beta = 0.96;      % 年度贴现因子
    cS.sigma = 2.0;      % 跨期替代弹性的倒数
    cS.cFloor = 1e-5;    % 消费最低值
    
    % 生产和税收
    cS.alpha = 0.36;     % 资本份额
    cS.ddk = 0.08;       % 年度折旧率
    cS.A = 1.0;          % TFP
    cS.tau_k = 0.10;     % 资本利得税
    cS.tau_l = 0.20;     % 劳动收入税
    cS.tau_c = 0.05;     % 消费税
    
    % 劳动冲击过程
    cS.nw = 5;           % 生产率冲击格点数
    
    % 政府
    cS.gov_exp_frac_Y = 0; % 为简化，假设G=T，并且G不由Y的分数决定，而是由税收内生决定
end

function cS = generateGrids(cS)
    cS.kMin = 0;
    cS.kMax = 50;
    power_k = 1.5; % 网格密度参数
    kGridV_temp = cS.kMin + (cS.kMax - cS.kMin) * (linspace(0, 1, cS.nk).^power_k);
    cS.kGridV = kGridV_temp(:);
end

function [leLogGridV, leTrProbM, leStationaryDist_e] = EarningProcess(cS)
    lePersistence = 0.95;
    leShockStd = 0.15;
    [leLogGridV_raw, leTrProbM] = tauchen(cS.nw, lePersistence, leShockStd, 0, 2.5);
    leLogGridV = leLogGridV_raw - mean(leLogGridV_raw);
    
    % 计算平稳分布
    [V, D] = eig(leTrProbM');
    [~, idx] = min(abs(diag(D) - 1));
    leStationaryDist_e = V(:,idx) ./ sum(V(:,idx));
end

function [y_grid_out, trProbM_out] = tauchen(N, rho, sigma, mu, m)
    std_y = sqrt(sigma^2 / (1-rho^2));
    y_max = m*std_y; y_min = -y_max;
    y = linspace(y_min, y_max, N);
    d = y(2)-y(1);
    trProbM_out = zeros(N,N);
    for j=1:N
        for k=1:N
            m_k = rho*y(j) + mu;
            if k==1, trProbM_out(j,k) = normcdf((y(1)-m_k+d/2)/sigma);
            elseif k==N, trProbM_out(j,k) = 1 - normcdf((y(N)-m_k-d/2)/sigma);
            else, trProbM_out(j,k) = normcdf((y(k)-m_k+d/2)/sigma) - normcdf((y(k)-m_k-d/2)/sigma); end
        end
    end
    y_grid_out = y(:);
end


% -------------------------------------------------------------------------
% --- 2. 稳态求解器 (高层) ---
% -------------------------------------------------------------------------
function [ss, eq_found, Dist, k_prime_policy] = solve_steady_state_iter_S(cS, paramS)
    excess_saving_wrapper = @(K_guess) fzero_wrapper_S(K_guess, cS, paramS);
    k_bracket=[0.1, 2]; options=optimset('TolX',1e-9);
    [K_eq,fval,exitflag]=fzero(excess_saving_wrapper,k_bracket,options); eq_found=(exitflag>0);
    if ~eq_found, warning('iter_S 未能找到均衡解'); ss=struct(); Dist=[]; k_prime_policy=[]; return; end
    fprintf('   iter_S 收敛！均衡资本 K_eq = %.8f, 最终净储蓄 S_net = %.3e\n', K_eq, fval);
    [~, ss, Dist, k_prime_policy] = calculate_net_saving(K_eq, cS, paramS);
end

function [ss, eq_found, Dist, k_prime_policy] = solve_steady_state_iter_K(cS, paramS)
    capital_market_wrapper = @(K_guess) fzero_wrapper_K(K_guess, cS, paramS);
    k_bracket=[0.1, 2]; options=optimset('TolX',1e-9,'Display','off');
    [K_eq,fval,exitflag]=fzero(capital_market_wrapper,k_bracket,options); eq_found=(exitflag>0);
    if ~eq_found, warning('iter_K 未能找到均衡解'); ss=struct(); Dist=[]; k_prime_policy=[]; return; end
    fprintf('   iter_K 收敛！均衡资本 K_eq = %.8f, 最终K误差 = %.3e\n', K_eq, fval);
    [~, ss, Dist, k_prime_policy] = calculate_net_saving(K_eq, cS, paramS);
end

function S_net_out = fzero_wrapper_S(K_guess, cS, paramS)
    [S_net_out] = calculate_net_saving(K_guess, cS, paramS);
end

function K_error_out = fzero_wrapper_K(K_guess, cS, paramS)
    [~,~,~,~,K_model]=calculate_net_saving(K_guess, cS, paramS);
    K_error_out = K_guess - K_model;
end


% -------------------------------------------------------------------------
% --- 3. 【核心计算引擎】 ---
% -------------------------------------------------------------------------
function [S_net, ss, Dist, k_prime_policy, K_model_out] = calculate_net_saving(K_guess, cS, paramS)
    % --- 1. 外循环给定 K_guess，计算价格和劳动供给 ---
    if K_guess <= 0, K_guess = 1e-8; end
    L_ss = sum(paramS.leStationaryDist_e .* paramS.leGridV); % 总劳动供给由外生冲击的平稳分布决定
    M_prices = get_prices_at_t(K_guess, L_ss, cS.A, cS);
    
    % --- 2. 求解家庭问题 (VFI) ---
    [k_prime_policy, ~] = HHSolution_VFI_Aiyagari(M_prices, cS, paramS);
    
    % --- 3. 求解稳态分布 ---
    k_prime_idx = get_policy_index_matrix(k_prime_policy, cS);
    Dist = solve_stationary_distribution_Aiyagari(k_prime_idx, cS, paramS);
    
    % --- 4. 聚合得到宏观量 ---
    [K_model_out, C_final, Tax_final] = aggregate_Aiyagari(Dist, k_prime_policy, M_prices, cS, paramS);
    
    % --- 5. 计算市场出清误差 ---
    G_final = Tax_final; % 政府预算平衡
    GDP = M_prices.Y;
    Depreciation = cS.ddk * K_guess;
    NDP = GDP - Depreciation;
    S_net = NDP - C_final - G_final;
    
    % --- 6. (可选) 打包所有宏观结果 ---
    if nargout > 1
        I_final = Depreciation; % 稳态下，投资等于折旧
        Y_exp = C_final + I_final + G_final;
        ss = struct();
        ss.K_physical = K_guess; ss.L = L_ss; ss.Y = Y_exp; ss.Y_from_production = GDP;
        ss.C = C_final; ss.I = I_final; ss.G = G_final; ss.w = M_prices.w; ss.r_mkt = M_prices.r_mkt;
    end
end

function M_prices = get_prices_at_t(K, L, A, cS)
    if K <= 0, K = 1e-8; end
    Y = A .* (K.^cS.alpha) .* (L.^(1-cS.alpha));
    MPK = cS.alpha .* Y ./ K;
    w_t = (1-cS.alpha) .* Y ./ L;
    r_mkt_t = MPK - cS.ddk;
    M_prices = struct('K', K, 'L', L, 'Y', Y, 'w', w_t, 'r_mkt', r_mkt_t);
end

function [kPolM, valM] = HHSolution_VFI_Aiyagari(M_prices, cS, paramS)
    % 使用标准的价值函数迭代法求解
    max_iter = 1000; tol = 1e-7; damp = 0.5;
    
    valM_old = zeros(cS.nk, cS.nw);
    kPolM = zeros(cS.nk, cS.nw);
    kPol_idx = zeros(cS.nk, cS.nw, 'uint16');
    
    r_net = M_prices.r_mkt * (1 - cS.tau_k);
    
    for iter = 1:max_iter
        % 1. 构建期望价值函数的插值器
        EV_interp = griddedInterpolant(cS.kGridV, valM_old, 'linear');
        
        % 2. 计算期望价值 E[V(k',e')]
        EV = EV_interp(cS.kGridV) * paramS.leTrProbM';
        
        valM_new = zeros(cS.nk, cS.nw);
        
        % 3. 逐个状态求解
        for ie = 1:cS.nw
            e_val = paramS.leGridV(ie);
            labor_income = M_prices.w * e_val;
            labor_tax = cS.tau_l * labor_income;
            
            for ik = 1:cS.nk
                k_now = cS.kGridV(ik);
                capital_income = M_prices.r_mkt * k_now;
                capital_tax = cS.tau_k * capital_income;
                
                cash_on_hand = (1+r_net)*k_now + labor_income * (1-cS.tau_l);
                % 预算约束: c(1+tau_c) + k' = (1+r_net)k + (1-tau_l)w*e
                
                % 向量化求解最优 k'
                c_expend_v = cash_on_hand - cS.kGridV';
                c_v = c_expend_v ./ (1 + cS.tau_c);
                
                util_v = -inf(size(c_v));
                valid_c_idx = c_v > cS.cFloor;
                util_v(valid_c_idx) = (c_v(valid_c_idx).^(1 - cS.sigma)) ./ (1 - cS.sigma);
                
                v_options = util_v + cS.beta .* EV(:, ie)';
                
                [max_val, idx_k_prime] = max(v_options);
                
                valM_new(ik, ie) = max_val;
                kPol_idx(ik, ie) = idx_k_prime;
                kPolM(ik, ie) = cS.kGridV(idx_k_prime);
            end
        end
        
        % 4. 检查收敛
        dist_v = max(abs(valM_new(:) - valM_old(:)));
        if dist_v < tol
            % fprintf('   VFI 在 %d 次迭代后收敛, 距离 = %.3e\n', iter, dist_v);
            break;
        end
        
        valM_old = (1-damp)*valM_old + damp*valM_new;
    end
    if iter == max_iter, warning('VFI 未在最大迭代次数内收敛'); end
    valM = valM_old;
end

function Dist = solve_stationary_distribution_Aiyagari(k_prime_idx, cS, paramS)
    % 迭代求解不变分布
    max_iter = 5000; tol = 1e-10;
    
    Dist_old = ones(cS.nk, cS.nw) ./ (cS.nk * cS.nw); % 从均匀分布开始
    
    for iter = 1:max_iter
        Dist_new = zeros(cS.nk, cS.nw);
        for ik = 1:cS.nk
            for ie = 1:cS.nw
                prob_mass = Dist_old(ik, ie);
                if prob_mass > 1e-20
                    ik_prime = k_prime_idx(ik, ie);
                    e_trans_probs = paramS.leTrProbM(ie, :);
                    Dist_new(ik_prime, :) = Dist_new(ik_prime, :) + prob_mass * e_trans_probs;
                end
            end
        end
        
        dist_d = max(abs(Dist_new(:) - Dist_old(:)));
        if dist_d < tol, break; end
        Dist_old = Dist_new;
    end
    if iter == max_iter, warning('分布迭代未在最大迭代次数内收敛'); end
    Dist = Dist_old;
end

function k_prime_idx=get_policy_index_matrix(kPolM, cS)
    k_prime_idx = zeros(cS.nk, cS.nw, 'uint16');
    for ie = 1:cS.nw
        for ik = 1:cS.nk
            [~, idx] = min(abs(cS.kGridV - kPolM(ik, ie)));
            k_prime_idx(ik, ie) = idx;
        end
    end
end

% -------------------------------------------------------------------------
% --- 4. 聚合与核算 ---
% -------------------------------------------------------------------------

function [K_agg, C_agg, Tax_agg] = aggregate_Aiyagari(Dist, k_prime_policy, M_prices, cS, paramS)
    K_agg = 0; C_agg = 0; Tax_agg = 0;
    
    for ie = 1:cS.nw
        for ik = 1:cS.nk
            mass = Dist(ik, ie);
            if mass < 1e-20, continue; end
            
            k_now = cS.kGridV(ik);
            e_val = paramS.leGridV(ie);
            k_prime = k_prime_policy(ik, ie);
            
            [c_val, tax_val] = backout_accounting_Aiyagari(k_now, k_prime, e_val, M_prices, cS);
            
            C_agg = C_agg + c_val * mass;
            Tax_agg = Tax_agg + tax_val * mass;
            K_agg = K_agg + k_prime * mass; % 下一期的资本存量
        end
    end
end

function [c_val, tax_total] = backout_accounting_Aiyagari(k_now, k_prime, e_val, M_prices, cS)
    labor_income = M_prices.w * e_val;
    capital_income = M_prices.r_mkt * k_now;
    
    labor_tax = cS.tau_l * labor_income;
    capital_tax = cS.tau_k * capital_income;
    
    % 从预算约束反推消费支出
    % budget: c(1+tau_c) + k' = (1+r)k + w*e - tau_k*r*k - tau_l*w*e
    % 也可以写成: c(1+tau_c) + k' = (1+r(1-tau_k))k + w*e(1-tau_l)
    r_net = M_prices.r_mkt * (1 - cS.tau_k);
    cash_available = (1 + r_net) * k_now + labor_income * (1 - cS.tau_l);
    
    c_expenditure = cash_available - k_prime;
    c_val = max(cS.cFloor, c_expenditure / (1 + cS.tau_c));
    
    consumption_tax = c_val * cS.tau_c;
    tax_total = labor_tax + capital_tax + consumption_tax;
end


% -------------------------------------------------------------------------
% --- 5. 验证函数 ---
% -------------------------------------------------------------------------

function validate_national_accounts(ss)
    fprintf('\n--- 国民账户一致性检验 ---\n');
    Y_exp=ss.Y; Y_prod=ss.Y_from_production;
    fprintf('   Y_prod vs Y_exp 残差: %.3e\n', Y_prod - Y_exp);
    
    Depreciation=ss.I;
    Y_ndp = Y_prod - Depreciation;
    NetLaborIncome = ss.w * ss.L * (1-cS.tau_l); % 税后
    NetCapitalIncome = ss.r_mkt * ss.K_physical * (1-cS.tau_k); % 税后
    
    % 国民收入核算: Y_ndp = C + I_net + G
    % 净投资 I_net = I - Depreciation = 0 in steady state.
    % Y_ndp = C + G.
    fprintf('   NDP vs (C+G) 残差: %.3e\n', Y_ndp - (ss.C + ss.G));
    
    % S_net = Y_ndp - C - G.
    NationalSaving_from_NDI = Y_ndp - ss.C - ss.G;
    NetInvestment = ss.I - Depreciation; % 稳态下净投资为0
    fprintf('   S_net vs I_net 残差: %.3e\n', NationalSaving_from_NDI - NetInvestment);
    
    if abs(Y_prod-Y_exp)<1e-7 && abs(NationalSaving_from_NDI-NetInvestment)<1e-7
        fprintf('   ✅ 宏观核算通过。\n');
    else
        fprintf('   ❌ 宏观核算失败。\n');
    end
end

function validate_micro_budget_constraints(ss, Dist, k_prime_policy, cS, paramS)
    M_sim = struct('r_mkt', ss.r_mkt, 'w', ss.w);
    max_abs_residual = 0;
    
    for ie = 1:cS.nw
        for ik = 1:cS.nk
            if Dist(ik, ie) < 1e-20, continue; end
            
            k_now = cS.kGridV(ik);
            e_val = paramS.leGridV(ie);
            k_prime = k_prime_policy(ik, ie);
            
            [c_val, tax_val] = backout_accounting_Aiyagari(k_now, k_prime, e_val, M_sim, cS);
            
            % 现金流入
            cash_in = (1 + M_sim.r_mkt) * k_now + M_sim.w * e_val;
            % 现金流出
            cash_out = c_val * (1 + cS.tau_c) + tax_val - cS.tau_c*c_val + k_prime; % tax_val已包含消费税，要加回来
            
            budget_residual = cash_in - (tax_val + c_val + k_prime);
            
            if abs(budget_residual) > max_abs_residual
                max_abs_residual = abs(budget_residual);
            end
        end
    end
    
    fprintf('   最大个体预算残差: %.3e\n', max_abs_residual);
    if max_abs_residual < 1e-9
        fprintf('   ✅ 微观审计通过。\n');
    else
        fprintf('   ❌ 微观审计失败。\n');
    end
end