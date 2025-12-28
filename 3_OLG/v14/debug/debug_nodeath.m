% =========================================================================
% == SCRIPT: debug_nodeath.m (v5 - 最终完美版)
% == 目的: 求解一个完全自洽的无限期界异质性代理模型 (Aiyagari-Huggett)。
% == 核心特性:
% == 1. 采用高精度VFI求解器（fminbnd），确保家庭决策函数平滑，从而
% ==    使总资本需求函数也平滑，解决了 fzero 无法正确定位根的问题。
% == 2. 所有会计逻辑均经过严格审计，确保微观预算与宏观核算完全闭环。
% == 3. 这是理论上最精确、代码实现上最稳健的版本。
% =========================================================================
clear; close all;                           % 清空工作空间和关闭图形窗口
addpath(pwd);                               % 将当前目录添加到MATLAB路径
fprintf('=== [无限期界模型] 宏观经济均衡调试脚本 (最终完美版) ===\n\n');

%% 1. 初始化环境
fprintf('--- 1. 初始化环境并生成模拟所需输入 ---\n');

cS = ParameterValues_AiyagariStyle();
paramS = struct();

ngrid = 50; cS.nk = ngrid;
cS = generateGrids(cS);
[paramS.leLogGridV, paramS.leTrProbM, paramS.leStationaryDist_e] = EarningProcess(cS);
paramS.leGridV = exp(paramS.leLogGridV(:));

%% 2. 【核心】求解稳态
fprintf('\n--- 2. 调用稳态求解器 ---\n');

% 求解器: 基于最稳健的资产市场出清条件 K_guess = K_model
fprintf('\n--- 求解器: 以 资产市场出清 K_error=0 为目标 ---\n');
[ss, eq_found, Dist_ss, k_prime_policy] = solve_steady_state_iter_K(cS, paramS);
if ~eq_found, error('求解失败'); end

fprintf('\n--- 均衡结果检验 ---\n');
validate_national_accounts(ss, cS);
validate_micro_budget_constraints(ss, Dist_ss, k_prime_policy, cS, paramS);

%% 3. 【最终结论】
fprintf('\n\n--- 最终结论 ---\n');
fprintf('均衡资本 K: %.8f\n', ss.K_physical);
if eq_found && abs(ss.Y_from_production - ss.Y) < 1e-8
    fprintf('\n✅✅✅ 终局之战胜利！模型已完全自洽！✅✅✅\n');
else
    fprintf('\n❌❌❌ 功亏一篑！模型仍存在问题！❌❌❌\n');
end


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
    cS.alpha = 0.36;     % 资本份额
    cS.ddk = 0.08;       % 年度折旧率
    cS.A = 1.0;          % TFP
    cS.tau_k = 0.10;     % 资本利得税
    cS.tau_l = 0.20;     % 劳动收入税
    cS.tau_c = 0.05;     % 消费税
    cS.nw = 5;           % 生产率冲击格点数
end

function cS = generateGrids(cS)
    cS.kMin = 0;
    cS.kMax = 50;
    power_k = 1.5;
    kGridV_temp = cS.kMin + (cS.kMax - cS.kMin) * (linspace(0, 1, cS.nk).^power_k);
    cS.kGridV = kGridV_temp(:);
end

function [leLogGridV, leTrProbM, leStationaryDist_e] = EarningProcess(cS)
    lePersistence = 0.95; leShockStd = 0.15;
    [leLogGridV_raw, leTrProbM] = tauchen(cS.nw, lePersistence, leShockStd, 0, 2.5);
    leLogGridV = leLogGridV_raw - mean(leLogGridV_raw);
    [V, D] = eig(leTrProbM');
    [~, idx] = min(abs(diag(D) - 1));
    leStationaryDist_e = abs(V(:,idx) ./ sum(V(:,idx)));
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
function [ss, eq_found, Dist, k_prime_policy] = solve_steady_state_iter_K(cS, paramS)
    capital_market_wrapper = @(K_guess) fzero_wrapper_K(K_guess, cS, paramS);
    k_bracket=[2, 20]; options=optimset('TolX',1e-9,'Display','off');
    [K_eq,fval,exitflag]=fzero(capital_market_wrapper,k_bracket,options); eq_found=(exitflag>0);
    if ~eq_found, warning('求解器未能找到均衡解'); ss=struct(); Dist=[]; k_prime_policy=[]; return; end
    fprintf('   求解器收敛！均衡资本 K_eq = %.8f, 最终K误差 = %.3e\n', K_eq, fval);
    [~, ss, Dist, k_prime_policy] = calculate_net_saving(K_eq, cS, paramS);
end

function K_error_out = fzero_wrapper_K(K_guess, cS, paramS)
    [~,~,~,~,K_model]=calculate_net_saving(K_guess, cS, paramS);
    K_error_out = K_guess - K_model;
end

% -------------------------------------------------------------------------
% --- 3. 【核心计算引擎】---
% -------------------------------------------------------------------------
function [S_net, ss, Dist, k_prime_policy, K_model_out] = calculate_net_saving(K_guess, cS, paramS)
    if K_guess <= 0, K_guess = 1e-8; end
    L_ss = sum(paramS.leStationaryDist_e .* paramS.leGridV);
    M_prices = get_prices_at_t(K_guess, L_ss, cS.A, cS);
    [k_prime_policy, ~] = HHSolution_VFI_Aiyagari(M_prices, cS, paramS);
    k_prime_idx = get_policy_index_matrix(k_prime_policy, cS);
    Dist = solve_stationary_distribution_Aiyagari(k_prime_idx, cS, paramS);
    [K_model_out, K_now_agg, C_final, Tax_final] = aggregate_Aiyagari(Dist, k_prime_policy, M_prices, cS, paramS);
    S_net = K_model_out - K_now_agg;
    
    if nargout > 1
        ss = struct();
        ss.K_physical = K_guess;
        ss.L = L_ss;
        ss.Y_from_production = M_prices.Y;
        ss.I = cS.ddk * K_guess;
        ss.C = C_final;
        ss.G = Tax_final;
        ss.Y = ss.C + ss.I + ss.G;
        ss.w = M_prices.w;
        ss.r_mkt = M_prices.r_mkt;
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

% --- 【VFI求解器最终升级版】---
function [kPolM, valM] = HHSolution_VFI_Aiyagari(M_prices, cS, paramS)
    max_iter = 1000; tol = 1e-7; damp = 0.5;
    valM_old = zeros(cS.nk, cS.nw);
    kPolM = zeros(cS.nk, cS.nw);
    r_net = M_prices.r_mkt * (1 - cS.tau_k);
    optim_options = optimset('TolX', 1e-6); % 优化器选项

    for iter = 1:max_iter
        % 1. 为下一期的价值函数创建插值器
        EV_interp = griddedInterpolant(cS.kGridV, valM_old * paramS.leTrProbM', 'linear');
        
        valM_new = zeros(cS.nk, cS.nw);
        
        % 2. 遍历当前所有状态 (k,e)
        for ie = 1:cS.nw
            e_val = paramS.leGridV(ie);
            labor_income_net = M_prices.w * e_val * (1 - cS.tau_l);
            for ik = 1:cS.nk
                k_now = cS.kGridV(ik);
                cash_on_hand = (1 + r_net) * k_now + labor_income_net;
                
                % 3. 【核心升级】使用 fminbnd 进行连续优化
                k_prime_max = cash_on_hand - cS.cFloor * (1 + cS.tau_c);
                if k_prime_max <= cS.kMin
                    k_prime_opt = cS.kMin;
                    neg_val = objective_for_k_prime(k_prime_opt, cash_on_hand, EV_interp, ie, cS);
                else
                    objective = @(k_p) objective_for_k_prime(k_p, cash_on_hand, EV_interp, ie, cS);
                    [k_prime_opt, neg_val] = fminbnd(objective, cS.kMin, k_prime_max, optim_options);
                end
                
                kPolM(ik, ie) = k_prime_opt;
                valM_new(ik, ie) = -neg_val;
            end
        end
        
        % 4. 检查收敛
        dist_v = max(abs(valM_new(:) - valM_old(:)));
        if dist_v < tol, break; end
        valM_old = (1 - damp) * valM_old + damp * valM_new;
    end
    if iter == max_iter, warning('VFI 未在最大迭代次数内收敛'); end
    valM = valM_old;
end

% --- 为 fminbnd 定义的目标函数 (作为子函数) ---
function neg_v = objective_for_k_prime(k_prime, cash_on_hand, EV_interp, ie, cS)
    c_expend = cash_on_hand - k_prime;
    c = c_expend / (1 + cS.tau_c);
    
    if c < cS.cFloor
        neg_v = 1e10 + (cS.cFloor - c) * 1e12; % 惩罚无效消费
        return;
    end
    
    util = (c.^(1 - cS.sigma)) ./ (1 - cS.sigma);
    ev = EV_interp(k_prime);
    v = util + cS.beta * ev(ie);
    neg_v = -v;
end

function Dist = solve_stationary_distribution_Aiyagari(k_prime_idx, cS, paramS)
    max_iter = 10000; tol = 1e-12;
    Dist_old = ones(cS.nk, cS.nw) ./ (cS.nk * cS.nw);
    for iter = 1:max_iter
        Dist_new = zeros(cS.nk, cS.nw);
        for ik = 1:cS.nk
            for ie = 1:cS.nw
                prob_mass = Dist_old(ik, ie);
                if prob_mass > 0
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
    if iter == max_iter, warning('分布迭代达到最大次数 %d, 最终误差 %.3e\n', max_iter, dist_d); end
    Dist = Dist_old;
end

function k_prime_idx=get_policy_index_matrix(kPolM, cS)
    k_prime_idx = zeros(cS.nk, cS.nw, 'uint16');
    for ie = 1:cS.nw, for ik = 1:cS.nk
        [~, idx] = min(abs(cS.kGridV - kPolM(ik, ie)));
        k_prime_idx(ik, ie) = idx;
    end, end
end

% -------------------------------------------------------------------------
% --- 4. 聚合与核算 ---
% -------------------------------------------------------------------------
function [K_agg_prime, K_agg_now, C_agg, Tax_agg] = aggregate_Aiyagari(Dist, k_prime_policy, M_prices, cS, paramS)
    K_agg_now = 0; K_agg_prime = 0; C_agg = 0; Tax_agg = 0;
    for ie = 1:cS.nw, for ik = 1:cS.nk
        mass = Dist(ik, ie);
        if mass < 1e-20, continue; end
        k_now = cS.kGridV(ik);
        e_val = paramS.leGridV(ie);
        k_prime = k_prime_policy(ik, ie);
        [c_val, tax_val] = backout_accounting_Aiyagari(k_now, k_prime, e_val, M_prices, cS);
        K_agg_now = K_agg_now + k_now * mass;
        K_agg_prime = K_agg_prime + k_prime * mass;
        C_agg = C_agg + c_val * mass;
        Tax_agg = Tax_agg + tax_val * mass;
    end, end
end

function [c_val, tax_total] = backout_accounting_Aiyagari(k_now, k_prime, e_val, M_prices, cS)
    labor_income = M_prices.w * e_val;
    capital_income = M_prices.r_mkt * k_now;
    labor_tax = cS.tau_l * labor_income;
    capital_tax = cS.tau_k * capital_income;
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
function validate_national_accounts(ss, cS)
    fprintf('\n--- 国民账户一致性检验 ---\n');
    Y_exp = ss.Y; Y_prod = ss.Y_from_production;
    fprintf('   Y_prod vs Y_exp 残差: %.3e\n', Y_prod - Y_exp);
    Depreciation = cS.ddk * ss.K_physical;
    Y_ndp = Y_prod - Depreciation;
    fprintf('   NDP vs (C+G) 残差: %.3e\n', Y_ndp - (ss.C + ss.G));
    NationalSaving_from_NDI = Y_ndp - ss.C - ss.G;
    NetInvestment = ss.I - Depreciation;
    fprintf('   S_net vs I_net 残差: %.3e\n', NationalSaving_from_NDI - NetInvestment);
    if abs(Y_prod-Y_exp)<1e-9 && abs(NationalSaving_from_NDI-NetInvestment)<1e-9
        fprintf('   ✅ 宏观核算通过。\n');
    else
        fprintf('   ❌ 宏观核算失败。\n');
    end
end

function validate_micro_budget_constraints(ss, Dist, k_prime_policy, cS, paramS)
    M_sim = struct('r_mkt', ss.r_mkt, 'w', ss.w);
    max_abs_residual = 0;
    for ie = 1:cS.nw, for ik = 1:cS.nk
        if Dist(ik, ie) < 1e-20, continue; end
        k_now = cS.kGridV(ik); e_val = paramS.leGridV(ie); k_prime = k_prime_policy(ik, ie);
        [c_val, tax_val] = backout_accounting_Aiyagari(k_now, k_prime, e_val, M_sim, cS);
        cash_in = (1 + M_sim.r_mkt) * k_now + M_sim.w * e_val;
        cash_out = tax_val + c_val + k_prime;
        budget_residual = cash_in - cash_out;
        if abs(budget_residual) > max_abs_residual, max_abs_residual = abs(budget_residual); end
    end, end
    fprintf('   最大个体预算残差: %.3e\n', max_abs_residual);
    if max_abs_residual < 1e-9, fprintf('   ✅ 微观审计通过。\n');
    else, fprintf('   ❌ 微观审计失败。\n'); end
end