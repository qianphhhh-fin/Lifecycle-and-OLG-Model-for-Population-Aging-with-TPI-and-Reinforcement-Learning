% =========================================================================
% == SCRIPT: Cohort Tracking Debugger
% == 目的：创建一个专门的调试环境，通过跟踪一个家庭的完整生命周期，
% ==       来验证 simulate_forward_with_policies 函数中的代际更替
% ==       （aging）和资产转移逻辑是否正确实现。
% =========================================================================
clear; close all;
addpath(pwd);
fprintf('=== 转型路径个体生命周期轨迹跟踪调试脚本 ===\n\n');

%% 1. 初始化环境 (与主调试脚本一致)
fprintf('--- 1. 初始化环境并生成模拟所需输入 ---\n');
cS = main_olg_v14_utils.ParameterValues_HuggettStyle();
paramS = struct();
 
ngrid = 30; cS.nk = ngrid; cS.nkpps = ngrid; cS.nkprime = ngrid; cS.npps = 5;
cS = main_olg_v14_utils.generateGrids(cS);
[Z_path, A_path, T_sim] = main_olg_v14_utils.load_exogenous_paths(cS);
cS = main_olg_v14_utils.calcaulte_theta_payg_path(cS, false);
cS.sim_years = cS.start_year:cS.time_Step:(cS.start_year + (T_sim-1)*cS.time_Step);
[paramS.leLogGridV, paramS.leTrProbM, paramS.leProb1V] = main_olg_v14_utils.EarningProcess_olgm(cS);
paramS.leGridV = exp(paramS.leLogGridV(:));
eIdxM = main_olg_v14_utils.LaborEndowSimulation_olgm_AgeGroup(cS, paramS);

% 求解初始稳态以获得初始分布
[ss_initial, eq_found, initial_dist_k, initial_dist_kpps] = main_olg_v14_utils.solve_steady_state_with_fund(Z_path(:,1), 0, cS, paramS, eIdxM);
if ~eq_found, error('初始稳态未求解成功'); end

% 猜测资本路径并计算价格路径
K_guess_path = linspace(ss_initial.K_total, ss_initial.K_total * 1.2, T_sim)';
w_path = zeros(T_sim, 1);
r_mkt_path = zeros(T_sim, 1);
b_path = zeros(T_sim, 1);
L_path = zeros(T_sim, 1);
for t = 1:T_sim
    [~, L_t] = main_olg_v14_utils.LaborSupply_Huggett(eIdxM, cS, paramS, Z_path(:,t));
    L_path(t) = L_t;
    M_prices_t = main_olg_v14_utils.get_prices_at_t(K_guess_path(t), L_t, A_path(t), cS);
    r_mkt_path(t) = M_prices_t.r_mkt_t;
    w_path(t) = M_prices_t.w_t;
    mass_workers_t = sum(Z_path(1:cS.aR_new, t));
    if mass_workers_t > 0, b_path(t) = cS.rho_prime_payg * (w_path(t) * L_t / mass_workers_t); else, b_path(t) = 0; end
end

% 加载预先计算好的策略函数
load('PolicyFunctions.mat');
fprintf('✅ 环境初始化完毕，所有模拟输入已准备就绪。\n\n');


%% 2. 运行带跟踪功能的模拟
fprintf('--- 2. 运行带有内部状态跟踪的前向模拟 ---\n');

% 调用一个特殊的、带跟踪功能的本地函数来执行模拟
[~, k_panel_history, kpps_panel_history] = simulate_forward_with_policies_and_track(...
    initial_dist_k, initial_dist_kpps, PolicyFunctions, T_sim, Z_path, A_path, ...
    w_path, r_mkt_path, b_path, L_path, cS, paramS, eIdxM);

fprintf('✅ 模拟完成，已捕获所有时期的资产面板历史。\n\n');


%% 3. 逐期跟踪单个家庭的生命周期
fprintf('--- 3. 跟踪单个家庭 (Household #1) 的生命周期轨迹 ---\n');

hh_to_track = 1; % 我们要跟踪的家庭编号
max_error = 0;
life_cycle_length = cS.aD_new; % 一个家庭的生命周期长度

% 存储被跟踪家庭的数据
tracked_k_path = zeros(life_cycle_length, 1);
tracked_kpps_path = zeros(life_cycle_length, 1);

fprintf('追踪家庭 #%d 的生命周期 (共 %d 个时期)...\n', hh_to_track, life_cycle_length);
fprintf('%s\n', repmat('-', 100, 1));
fprintf('%-5s | %-5s | %-12s | %-12s | %-15s | %-15s | %-12s\n', ...
    'Time', 'Age', 'K_start', 'KPPS_start', 'K_prime_made', 'KPPS_prime_made', 'Error(K)');
fprintf('%s\n', repmat('-', 100, 1));

% 我们跟踪一个从 t=1, a=1 开始的生命周期的个体
% 这个循环的长度是其生命周期长度
for t = 1:life_cycle_length-1
    a = t; % 在时期t，这个同生群的年龄是a
    
    % a. 从历史记录中获取 t 期、a 年龄的期初资产
    k_panel_t = k_panel_history{t};
    kpps_panel_t = kpps_panel_history{t};
    
    k_start_of_period = k_panel_t(hh_to_track, a);
    kpps_start_of_period = kpps_panel_t(hh_to_track, a);
    
    % 存储路径
    tracked_k_path(a) = k_start_of_period;
    tracked_kpps_path(a) = kpps_start_of_period;
    
    % b. 根据 t 期的策略函数，计算该家庭在 t 期做出的储蓄决策 k'
    Policies_t = PolicyFunctions(t, :);
    Policy_age = Policies_t{a};
    kPolInterp_cell = Policy_age{1};
    cPpsPolInterp_cell = Policy_age{2};
    
    ie = eIdxM(hh_to_track, a); % 获取该家庭在此时的生产力冲击状态
    kPolInterp = kPolInterp_cell{ie};
    cPpsPolInterp = cPpsPolInterp_cell{ie};
    
    if numel(kPolInterp.GridVectors) > 1
        k_prime_decision = kPolInterp(k_start_of_period, kpps_start_of_period);
        cpps_decision = cPpsPolInterp(k_start_of_period, kpps_start_of_period);
    else
        k_prime_decision = kPolInterp(k_start_of_period);
        cpps_decision = cPpsPolInterp(k_start_of_period);
    end
    k_prime_decision(isnan(k_prime_decision)) = cS.kMin;
    cpps_decision(isnan(cpps_decision)) = 0;

    % 计算kpps'
    pps_return_factor = 1 + r_mkt_path(t) + cS.pps_return_rate_premium;
    pps_withdrawal_gross = 0;
    if a > cS.aR_new && cS.pps_active, pps_withdrawal_gross = kpps_start_of_period .* cS.pps_withdrawal_rate; end
    k_pps_prime_decision = (kpps_start_of_period - pps_withdrawal_gross + cpps_decision) * pps_return_factor;
    k_pps_prime_decision = max(cS.kppsMin, min(cS.kppsMax, k_pps_prime_decision));

    % c. 从历史记录中获取 t+1 期、a+1 年龄的期初资产 (这是实际模拟的结果)
    k_panel_tp1 = k_panel_history{t+1};
    k_actual_next_period = k_panel_tp1(hh_to_track, a+1);
    
    % d. 比较决策与结果
    err = abs(k_prime_decision - k_actual_next_period);
    max_error = max(max_error, err);
    
    fprintf('%-5d | %-5d | %-12.4f | %-12.4f | %-15.4f | %-15.4f | %-12.3e\n', ...
        t, a, k_start_of_period, kpps_start_of_period, k_prime_decision, k_pps_prime_decision, err);
end
fprintf('%s\n', repmat('-', 100, 1));


%% 4. 结论
if max_error < 1e-12
    fprintf('\n✅ 验证通过！个体家庭的资产代际转移逻辑完全正确。\n');
    fprintf('   最大误差为 %.3e，在可接受的数值精度范围内。\n', max_error);
else
    fprintf('\n❌ 验证失败！个体家庭的资产代际转移存在错误。\n');
    fprintf('   最大误差为 %.3e，超过了数值精度范围。\n', max_error);
end

% 可视化被跟踪家庭的资产路径
figure('Name', '被跟踪家庭的生命周期资产路径');
plot(1:life_cycle_length, tracked_k_path, '-ob', 'DisplayName', '普通资产 (k)');
hold on;
plot(1:life_cycle_length, tracked_kpps_path, '-sr', 'DisplayName', '养老金资产 (kpps)');
xlabel('年龄组 (Age group)');
ylabel('资产水平');
title(sprintf('家庭 #%d 的生命周期资产积累路径', hh_to_track));
legend('show', 'Location', 'best');
grid on;


%% 本地函数：带跟踪功能的模拟器
function [Results, k_panel_history, kpps_panel_history] = simulate_forward_with_policies_and_track(initial_dist_k, initial_dist_kpps, PolicyFunctions, T_sim, Z_path, A_path, w_path, r_mkt_path, b_path, L_path, cS, paramS, eIdxM)
    % 这是一个本地函数，是 simulate_forward_with_policies 的一个特殊版本。
    % 它额外返回 k_panel_history 和 kpps_panel_history，以便我们进行检查。
    
    k_panel = initial_dist_k(:, 1:cS.aD_new);
    kpps_panel = initial_dist_kpps(:, 1:cS.aD_new);
    
    Results = struct();
    k_panel_history = cell(T_sim, 1);
    kpps_panel_history = cell(T_sim, 1);

    for t = 1:T_sim
        % *** 核心跟踪步骤：在循环开始时，存储当前状态 ***
        k_panel_history{t} = k_panel;
        kpps_panel_history{t} = kpps_panel;

        Z_t = Z_path(:, t);
        A_t = A_path(t);
        
        M_t = struct('w_t', w_path(t), 'r_mkt_t', r_mkt_path(t), 'r_net_period', r_mkt_path(t)*(1-cS.tau_k), 'b_t', b_path(t));
        M_t.current_t = t;
        cS_sim = cS;
        cS_sim.theta_t = cS.theta_path(t);
        if isfield(cS, 'sim_years'), cS_sim.pps_active = (cS.sim_years(t) >= cS_sim.pps_activation_year); else, cS_sim.pps_active = true; end
        Policies_t = PolicyFunctions(t, :);

        [kHist_t, kppsHist_t, ~, ~] = HHSimulation_olgm_transition(k_panel, kpps_panel, Policies_t, eIdxM, M_t, paramS, cS_sim);

        % 【正确的代际更替逻辑】
        k_prime_choices = kHist_t(:, 2:end);
        kpps_prime_choices = kppsHist_t(:, 2:end);
        k_panel_tp1 = zeros(size(k_panel));
        kpps_panel_tp1 = zeros(size(kpps_panel));
        k_panel_tp1(:, 2:cS.aD_new) = k_prime_choices(:, 1:cS.aD_new-1);
        kpps_panel_tp1(:, 2:cS.aD_new) = kpps_prime_choices(:, 1:cS.aD_new-1);
        k_panel = k_panel_tp1;
        kpps_panel = kpps_panel_tp1;
    end
    
    % 存储最后一次的面板，用于验证最后一步
    if T_sim > 0
        k_panel_history{T_sim + 1} = k_panel;
        kpps_panel_history{T_sim + 1} = kpps_panel;
    end
end

function [kHistM_out, kPpsHistM_out, cHistM_out, cppsHistM_out] = HHSimulation_olgm_transition(k_panel_in, kpps_panel_in, Policies_t, eIdxM_group, M_sim, paramS_sim, cS_sim)
    % =========================================================================
    % == 【函数已最终修正】HHSimulation_olgm_transition
    % == 核心修正：
    % == 1. 此函数现在【严格使用】传入的 k_panel_in 和 kpps_panel_in 作为
    % ==    每个家庭在每个年龄组的期初资产。
    % == 2. 不再初始化 kHistM_out 和 kPpsHistM_out 为零，而是用传入的
    % ==    面板数据作为历史记录的起点。
    % == 3. 循环中直接使用 k_panel_in(:, a_idx) 作为当前状态，而不是
    % ==    依赖于前一列的计算结果。
    % =========================================================================
    
    nSim = size(eIdxM_group, 1);
    aD = cS_sim.aD_new;
    
    % --- 1. 【核心修正】初始化输出矩阵，并将传入的面板作为历史记录的“期初”部分 ---
    kHistM_out = zeros(nSim, aD + 1);
    kPpsHistM_out = zeros(nSim, aD + 1);
    cHistM_out = zeros(nSim, aD);
    cppsHistM_out = zeros(nSim, aD);

    % kHistM_out 的前 aD 列代表 nSim 个家庭在 aD 个年龄组的期初资产
    kHistM_out(:, 1:aD) = k_panel_in;
    kPpsHistM_out(:, 1:aD) = kpps_panel_in;

    tr_per_capita_sim = 0; % 在转型路径中，意外遗赠通过其他机制处理

    % --- 2. 按年龄逐期模拟决策 ---
    for a_idx = 1:aD
        % --- 【核心修正】直接从传入的面板中获取当前状态 ---
        k_now = k_panel_in(:, a_idx);
        k_pps_now = kpps_panel_in(:, a_idx);

        % 从策略函数库中获取当前年龄的插值器
        Policy_age = Policies_t{a_idx};
        kPolInterp_cell = Policy_age{1};
        cPpsPolInterp_cell = Policy_age{2};

        k_next_decision = zeros(nSim, 1);
        cpps_decision = zeros(nSim, 1);
        
        % 遍历所有效率冲击状态，为每个家庭找到其决策
        for ie = 1:cS_sim.nw
            idx_sim = find(eIdxM_group(:, a_idx) == ie);
            if isempty(idx_sim), continue; end

            kPolInterp = kPolInterp_cell{ie};
            cPpsPolInterp = cPpsPolInterp_cell{ie};

            k_now_e = k_now(idx_sim);
            k_pps_now_e = k_pps_now(idx_sim);
            
            % 使用插值器得到决策
            if numel(kPolInterp.GridVectors) > 1
                k_next_decision(idx_sim) = kPolInterp(k_now_e, k_pps_now_e);
                cpps_decision(idx_sim) = cPpsPolInterp(k_now_e, k_pps_now_e);
            else
                k_next_decision(idx_sim) = kPolInterp(k_now_e);
                cpps_decision(idx_sim) = cPpsPolInterp(k_now_e);
            end
        end
        
        % 处理插值可能产生的NaN
        k_next_decision(isnan(k_next_decision)) = cS_sim.kMin;
        cpps_decision(isnan(cpps_decision)) = 0;
        
        % --- 3. 存储当期决策 k' 和 cpps ---
        cppsHistM_out(:, a_idx) = max(0, cpps_decision);
        
        % k' 决策被存储在历史矩阵的下一列
        kHistM_out(:, a_idx + 1) = max(cS_sim.kMin, k_next_decision);
        
        % --- 4. 反推消费 (此逻辑保持不变，它依赖于正确的 k_now 和 k_prime) ---
        tr_this_age = 0; if a_idx > 1, tr_this_age = tr_per_capita_sim; end
        if a_idx <= cS_sim.aR_new
            labor_income_gross = M_sim.w_t .* cS_sim.ageEffV_new(a_idx) .* paramS_sim.leGridV(eIdxM_group(:, a_idx));
            k_return = k_now .* (1 + M_sim.r_net_period);
            total_inflow = k_return + labor_income_gross + tr_this_age;
            payg_tax = cS_sim.theta_t .* labor_income_gross;
            labor_tax = cS_sim.tau_l .* max(0, labor_income_gross - cppsHistM_out(:, a_idx) - payg_tax);
            k_prime_outflow = kHistM_out(:, a_idx + 1);
            cpps_outflow = cppsHistM_out(:, a_idx);
            total_outflow_non_c = payg_tax + labor_tax + k_prime_outflow + cpps_outflow;
            c_expend = total_inflow - total_outflow_non_c;
        else
            k_return = k_now .* (1 + M_sim.r_net_period);
            b_age_val = M_sim.b_t;
            pps_withdrawal_gross = 0; if cS_sim.pps_active, pps_withdrawal_gross = k_pps_now .* cS_sim.pps_withdrawal_rate; end
            pps_withdrawal_net = pps_withdrawal_gross .* (1 - cS_sim.pps_tax_rate_withdrawal);
            total_inflow = k_return + b_age_val + tr_this_age + pps_withdrawal_net;
            k_prime_outflow = kHistM_out(:, a_idx + 1);
            c_expend = total_inflow - k_prime_outflow;
        end
        cHistM_out(:, a_idx) = max(cS_sim.cFloor, c_expend ./ (1 + cS_sim.tau_c));

        % --- 5. 更新PPS资产 (此逻辑保持不变) ---
        pps_return_factor = 1 + M_sim.r_mkt_t + cS_sim.pps_return_rate_premium;
        pps_withdrawal_gross_for_update = 0;
        if a_idx > cS_sim.aR_new && cS_sim.pps_active, pps_withdrawal_gross_for_update = k_pps_now .* cS_sim.pps_withdrawal_rate; end
        k_pps_after_withdrawal = k_pps_now - pps_withdrawal_gross_for_update;
        k_pps_next_unclamped = (k_pps_after_withdrawal + cppsHistM_out(:, a_idx)) * pps_return_factor;
        kPpsHistM_out(:, a_idx + 1) = max(cS_sim.kppsMin, min(cS_sim.kppsMax, k_pps_next_unclamped));
    end
end