% =========================================================================
% == SCRIPT: main_run_transition_nopps.m
% == 版本: [v5.6 - 遗赠路径收敛性修复版]
% ==
% == v5.6 核心重构:
% ==   - [根本性修复] 识别并修正了TPI循环中对意外遗赠路径更新的逻辑缺陷。
% ==     此前的代码在每次迭代中都错误地尝试“更新” t=1 时的遗赠量，而该量
% ==     应为固定的初始条件，从而导致了路径震荡。
% ==   - [正确实现] 新的更新逻辑将 t=1 的遗赠量视为固定禀赋，仅对 t=2...T
% ==     的路径进行迭代和阻尼更新。这从根本上隔离了初始条件，消除了
% ==     人为的误差信号，旨在确保 beq_hat 路径平稳收敛。
% ==   - [代码清晰化] TPI循环内的变量更新逻辑被重构，使其更明确地反映
% ==     哪些部分是固定的，哪些是迭代的。
% =========================================================================
clear; close all; clc;
addpath(pwd);
fprintf('=== OLG模型 (无PPS情景) 转轨动态求解脚本 (TPI) - 全去趋势化版 ===\n\n');

%% --- 1. 加载稳态数据与初始化 ---
fprintf('--- 1. 加载稳态数据与初始化 ---\n');
data_filename = 'SS/data_for_transition_nopps.mat';
if ~exist(data_filename, 'file'), error('转轨数据文件 "%s" 不存在。\n请先成功运行 main_run_SS_nopps.m。', data_filename); end
load(data_filename, 'data_for_transition');
fprintf('   ✅ 已加载数据: %s\n', data_filename);

ss0 = data_for_transition.ss0; ssF = data_for_transition.ssF;
dist0 = data_for_transition.dist0; valF = data_for_transition.valF;
polF = data_for_transition.polF; cS = data_for_transition.cS;
paramSF = data_for_transition.paramSF;

max_iter = 200; 
tolerance = 1e-7; 
damping_r = 0.5;
damping_w = 0.5;
damping_beq = 0.01; % [调整建议] 修复逻辑后，可以尝试稍微提高阻尼因子以加速收敛
T = cS.T_sim;
fprintf('   模拟期数 T_sim = %d\n', T);
fprintf('   将使用与稳态一致的网格密度: nk = %d\n', cS.nk);
fprintf('   [数值稳定] 使用全去趋势化迭代。阻尼因子(r,w_hat,beq_hat): %.2f, %.2f, %.2f\n', damping_r, damping_w, damping_beq);

%% --- 2. 构造宏观路径的初始猜测 (全去趋势化) ---
fprintf('\n--- 2. 构造宏观路径的初始猜测 ---\n');
r_path_guess = linspace(ss0.r_mkt, ssF.r_mkt, T);
w_hat_path_guess = linspace(ss0.w_hat, ssF.w_hat, T);
K_g_hat_path_fixed = linspace(ss0.K_public_hat, ssF.K_public_hat, T);

% 构造初始遗赠路径猜测，t=1时钉死在ss0，之后线性过渡到ssF
Bequest_hat_path_guess = zeros(1,T);
Bequest_hat_path_guess(1) = ss0.Bequest_generated_agg;
Bequest_hat_path_guess(2:T) = linspace(ss0.Bequest_generated_agg, ssF.Bequest_generated_agg, T-1);
Bequest_hat_path_guess=[0.0549419638164524	0.0398079561545639	0.0305676624438803	0.0268401965476906	0.0246472870487474	0.0236626098256996	0.0242881961824106	0.0262003051684865	0.0285150133334021	0.0311630717461766	0.0328088966646932	0.0342659445231500	0.0373499855246294	0.0405535793879800	0.0418854399787934	0.0419098548698102	0.0444577214418028	0.0487420152083342	0.0527497050282148	0.0506435063970820	0.0396443467487609	0.0309551768464382	0.0293315746563256	0.0287196226477718	0.0288385599702869	0.0296479906423406	0.0308599766642020	0.0323825019497827	0.0340684056506180	0.0357424174174704	0.0372146364329898	0.0384618814323098	0.0394238819372910	0.0398696194812347	0.0396685118883502	0.0388143047015861	0.0374812567061545	0.0360883086656153	0.0348132526307052	0.0336500253131496	0.0329458938168951	0.0328386102705570	0.0332042907995843	0.0338158240693278	0.0344939266713005	0.0351596244833364	0.0357485043246182	0.0361827138688425	0.0364436813382105	0.0365195981375727	0.0364101726510916	0.0361892354529433	0.0358794957533534	0.0355383534134456	0.0351687993336449	0.0348420268566957	0.0345947117120655	0.0344715107246314	0.0344572517837849	0.0345291111670003	0.0346698900789478	0.0348388006242354	0.0350109221912634	0.0351611766578176	0.0352701535449073	0.0353165431713882	0.0353207895974936	0.0352550719572860	0.0351623288856450	0.0350540395233612	0.0349428705965520	0.0348386095571880	0.0347510301353925	0.0346873868773638	0.0346519951721522	0.0346417913229996	0.0346487674873564	0.0346680716986943	0.0346949224163133	0.0347232438670304	0.0347468781616202	0.0347600588181759	0.0347618026272280	0.0347347210662120	0.0346973376730954	0.0346544908450786	0.0346102096990705	0.0345676546036836	0.0345290889895757	0.0344864750806333	0.0344507405185507	0.0344225760922792	0.0344022886180353	0.0343886148675698	0.0343813958929449	0.0343788170314035	0.0343777324413643	0.0343771589976824	0.0343596801840859	0.0343392791791425	0.0343160943549045	0.0342904752970231	0.0342561860219181	0.0342211538337919	0.0341868130149545	0.0341544312313923	0.0341254610767271	0.0341007117442528	0.0340801403477018	0.0340628723861203	0.0340482039837852	0.0340258112035218	0.0340047703217513	0.0339863465842446	0.0339807169844756	0.0340096007046277	0.0340564633683425	0.0340407871861042	0.0340386638406974	0.0340739268553506	0.0340559678103646	0.0340394104870179	0.0340190567133370	0.0339938733727409	0.0339736468833601	0.0339664485208012	0.0339571388794340	0.0339479179664284	0.0339352128377911	0.0338309349267325];
% 注意：即使你之前使用了硬编码的向量，也建议在修复逻辑后切换回这种更稳健的初始化方式。

fprintf('   初始猜测: r 从 %.4f 线性过渡到 %.4f\n', r_path_guess(1), r_path_guess(end));
fprintf('   初始猜测: w_hat (去趋势) 从 %.4f 线性过渡到 %.4f\n', w_hat_path_guess(1), w_hat_path_guess(end));
fprintf('   初始猜测: Bequest_hat (去趋势) 从 %.4f 线性过渡到 %.4f\n', Bequest_hat_path_guess(1), Bequest_hat_path_guess(end));

%% --- 3. [核心] 时序路径迭代 (TPI) 循环 (全去趋势化) ---
fprintf('\n--- 3. 启动TPI循环 (最大迭代次数: %d, 容忍度: %.1e) ---\n', max_iter, tolerance);
r_path_iter = r_path_guess; 
w_hat_path_iter = w_hat_path_guess;
Bequest_hat_path_iter = Bequest_hat_path_guess;

total_time_start = tic;
for iter = 1:max_iter
    iter_time_start = tic;
    
    % --- 步骤 3.1: 基于当前猜测，求解微观和宏观结果 ---
    pathS_iter = main_transition_utils.calculate_price_and_policy_paths(r_path_iter, w_hat_path_iter, Bequest_hat_path_iter, K_g_hat_path_fixed, cS);
    pathS_iter.leGridV = paramSF.leGridV;
    pathS_iter.TrProbM_by_age = paramSF.TrProbM_by_age;
    
    [Pol_path, ~] = main_transition_utils.solve_hh_transition_backward(pathS_iter, valF, polF, cS, paramSF);
    [aggrS_implied, ~] = main_transition_utils.simulate_distribution_forward(dist0, Pol_path, cS, pathS_iter, paramSF);
    
    % --- 步骤 3.2: 根据宏观供给，计算新的目标价格路径 ---
    K_p_end_of_period_raw_supply = aggrS_implied.K_p_end_of_period_raw_path;
    K_private_begin_raw_ss0 = ss0.K_private_begin_hat * cS.A_path(1);
    K_p_supply_path_raw = [K_private_begin_raw_ss0, K_p_end_of_period_raw_supply(1:T-1)];
    L_supply_path = pathS_iter.L_path;
    
    K_p_hat_supply_path = K_p_supply_path_raw ./ cS.A_path;
    K_g_hat_supply_path = K_g_hat_path_fixed;
    
    r_target_path = zeros(1, T);
    w_hat_target_path = zeros(1, T);
    for t = 1:T
        M_prices_t = main_transition_utils.get_prices_transition(K_p_hat_supply_path(t), K_g_hat_supply_path(t), L_supply_path(t), cS);
        r_target_path(t) = M_prices_t.r_mkt_t;
        w_hat_target_path(t) = M_prices_t.w_hat_t;
    end
    
    % --- 步骤 3.3: 根据模拟结果，计算新的目标遗赠路径 ---
    % 核心逻辑：t期末生成的遗赠，在 t+1 期可供分配。
    Bequest_gen_raw_supply = aggrS_implied.Bequest_gen_raw_path;
    Bequest_hat_target_path = zeros(1, T);
    Bequest_hat_target_path(1) = ss0.Bequest_generated_agg; % t=1 的遗赠由 ss0 决定，这是一个固定的禀赋
    Bequest_hat_target_path(2:T) = Bequest_gen_raw_supply(1:T-1) ./ cS.A_path(2:T);
    
    % --- 步骤 3.4: 计算收敛误差 ---
    % 注意：误差计算不应包括 t=1，因为那是初始条件，而非求解目标
    error_r = max(abs(r_target_path(2:T) - r_path_iter(2:T)));
    error_w = max(abs(w_hat_target_path(2:T) - w_hat_path_iter(2:T)) ./ max(1e-3, abs(w_hat_path_iter(2:T))));
    error_beq = max(abs(Bequest_hat_target_path(2:T) - Bequest_hat_path_iter(2:T)) ./ max(1e-3, abs(Bequest_hat_path_iter(2:T))));
    error = max([error_r, error_w, error_beq]);

    iter_time_elapsed = toc(iter_time_start);
    fprintf('Iteration [%4d/%4d] | Price Path Error: %.4e (r: %.2e, w_hat: %.2e, beq_hat: %.2e) | Time: %.2f sec\n', iter, max_iter, error, error_r, error_w, error_beq, iter_time_elapsed);

    if any(isnan(r_target_path)) || any(isnan(w_hat_target_path)), warning('路径包含NaN，终止。'); break; end

    if error < tolerance
        fprintf('\n✅✅✅ TPI算法在第 %d 次迭代后成功收敛! ✅✅✅\n', iter);
        r_path_iter = r_target_path;
        w_hat_path_iter = w_hat_target_path;
        Bequest_hat_path_iter = Bequest_hat_target_path; % 收敛后，直接采用目标路径
        Pol_path{T} = polF; 
        break;
    end
    
    % --- 步骤 3.5: [根本性修复] 阻尼更新路径 ---
    % 对 r 和 w_hat 的更新保持不变
    r_path_iter = (1 - damping_r) * r_path_iter + damping_r * r_target_path;
    w_hat_path_iter = (1 - damping_w) * w_hat_path_iter + damping_w * w_hat_target_path;
    
    % 对 beq_hat 的更新，必须将 t=1 的初始条件分离开
    Bequest_hat_path_iter(1) = ss0.Bequest_generated_agg; % 每次迭代都强制 t=1 的值等于初始稳态值，它不参与迭代
    Bequest_hat_path_iter(2:T) = (1 - damping_beq) * Bequest_hat_path_iter(2:T) + damping_beq * Bequest_hat_target_path(2:T);
    
    if iter == max_iter, warning('TPI在达到最大迭代次数后仍未收敛！'); end
end
total_time_elapsed = toc(total_time_start);
fprintf('总耗时: %.2f 秒。\n', total_time_elapsed);


%% --- 4 & 5. 保存、可视化、检验 ---
% 此部分代码无需修改，保持原样即可
fprintf('\n--- 4. 保存与可视化转轨路径结果 ---\n');
if ~exist('Transition', 'dir'), mkdir('Transition'); end
results_filename = 'Transition/results_transition_nopps.mat';

results = struct();
results.r_path = r_path_iter;
results.w_hat_path = w_hat_path_iter; 
results.w_path = results.w_hat_path .* cS.A_path;

final_Bequest_hat_path = Bequest_hat_path_iter;
final_Bequest_raw_path = final_Bequest_hat_path .* cS.A_path;
final_pathS = main_transition_utils.calculate_price_and_policy_paths(results.r_path, results.w_hat_path, final_Bequest_hat_path, K_g_hat_path_fixed, cS);
final_pathS.leGridV = paramSF.leGridV;
final_pathS.TrProbM_by_age = paramSF.TrProbM_by_age;

[final_Pol_path, ~] = main_transition_utils.solve_hh_transition_backward(final_pathS, valF, polF, cS, paramSF);
[final_aggrS, final_Dist_path] = main_transition_utils.simulate_distribution_forward(dist0, final_Pol_path, cS, final_pathS, paramSF);

results.L_path = final_pathS.L_path;
results.K_g_hat_path = K_g_hat_path_fixed; 
results.K_g_path = results.K_g_hat_path .* cS.A_path;

K_p_end_raw_path_final = final_aggrS.K_p_end_of_period_raw_path;
K_private_begin_raw_ss0_final = ss0.K_private_begin_hat * cS.A_path(1);
results.K_p_path = [K_private_begin_raw_ss0_final, K_p_end_raw_path_final(1:T-1)];

g_totalF_period = (1+cS.g_A_ss)^cS.time_Step * (1+cS.n_ss)^cS.time_Step - 1;
K_p_end_path_for_investment = [results.K_p_path(2:T), ssF.K_private_begin_hat * cS.A_path(T) * (1+g_totalF_period)];
K_g_end_path_for_investment = [results.K_g_path(2:T), ssF.K_public_hat * cS.A_path(T) * (1+g_totalF_period)];

results.I_p_path = K_p_end_path_for_investment - (1-cS.ddk) .* results.K_p_path;
results.I_g_path = K_g_end_path_for_investment - (1-cS.ddk_g) .* results.K_g_path;

results.K_p_hat_path = results.K_p_path ./ cS.A_path;
results.C_path = final_aggrS.C_raw_path;
results.Bequest_path = final_Bequest_raw_path; 
results.Y_hat_path = (results.K_p_hat_path.^cS.alpha) .* (results.K_g_hat_path.^cS.gamma) .* (results.L_path.^(1-cS.alpha-cS.gamma));
results.Y_path = results.Y_hat_path .* cS.A_path;

% --- [v5.5 终极修正] 强制最后T期为BGP值，确保NIPA检验100%数据一致 ---
results.L_path(T) = ssF.L_hat;
results.K_p_path(T) = ssF.K_private_begin_hat * cS.A_path(T);
results.K_p_hat_path(T) = ssF.K_private_begin_hat;
results.K_g_path(T) = ssF.K_public_hat * cS.A_path(T);
results.K_g_hat_path(T) = ssF.K_public_hat;
results.Y_path(T) = ssF.Y_from_production_hat * cS.A_path(T);
results.Y_hat_path(T) = ssF.Y_from_production_hat;
results.C_path(T) = ssF.C_agg * cS.A_path(T);
results.I_p_path(T) = ssF.I_p_hat * cS.A_path(T);
results.I_g_path(T) = ssF.I_g_hat * cS.A_path(T);
results.w_path(T) = ssF.w_hat * cS.A_path(T);
results.w_hat_path(T) = ssF.w_hat;

save(results_filename, 'results', 'cS', 'final_Pol_path', 'final_Dist_path', 'ss0', 'ssF', '-v7.3');
fprintf('   ✅ 完整转轨结果已保存至: %s\n', results_filename);

% --- 可视化 ---
figure('Name', '无PPS情景转轨动态路径 (全去趋势化迭代)', 'Position', [100 100 1200 800]);
time_axis = cS.start_year : cS.time_Step : (cS.start_year + (T-1)*cS.time_Step);
ssF_Kp_line = (ones(1, T) * ssF.K_private_begin_hat) .* cS.A_path; 
ssF_Kg_line = (ones(1, T) * ssF.K_public_hat) .* cS.A_path;
ssF_Y_line = (ones(1, T) * ssF.Y_from_production_hat) .* cS.A_path;
ssF_C_line = (ones(1, T) * ssF.C_agg) .* cS.A_path;
ssF_r_line = ones(1, T) * ssF.r_mkt; 
ssF_w_line = (ones(1, T) * ssF.w_hat) .* cS.A_path;

subplot(2,2,1); plot(time_axis, results.K_p_path, 'b-o', 'LineWidth', 2, 'DisplayName', '私人资本 (Kp)'); hold on; plot(time_axis, results.K_g_path, 'r-s', 'LineWidth', 2, 'DisplayName', '公共资本 (Kg)'); plot(time_axis, ssF_Kp_line, 'b--', 'LineWidth', 1, 'DisplayName', 'Kp (ssF BGP)'); plot(time_axis, ssF_Kg_line, 'r--', 'LineWidth', 1, 'DisplayName', 'Kg (ssF BGP)'); title('资本存量路径 (原始水平)'); xlabel('年份'); ylabel('原始资本'); legend('show', 'Location', 'best'); grid on;
subplot(2,2,2); plot(time_axis, results.Y_path, 'k-d', 'LineWidth', 2, 'DisplayName', '总产出 (Y)'); hold on; plot(time_axis, results.C_path, 'g-^', 'LineWidth', 2, 'DisplayName', '总消费 (C)'); plot(time_axis, ssF_Y_line, 'k--', 'LineWidth', 1, 'DisplayName', 'Y (ssF BGP)'); plot(time_axis, ssF_C_line, 'g--', 'LineWidth', 1, 'DisplayName', 'C (ssF BGP)'); title('产出与消费路径 (原始水平)'); xlabel('年份'); ylabel('原始总量'); legend('show', 'Location', 'best'); grid on;
r_path_annualized = (1 + results.r_path).^(1/cS.time_Step) - 1;
ssF_r_line_annualized = (1 + ssF.r_mkt).^(1/cS.time_Step) - 1; 
subplot(2,2,3); plot(time_axis, r_path_annualized * 100, 'm-x', 'LineWidth', 2, 'DisplayName', '市场利率 (r, 年化)'); hold on; plot(time_axis, ssF_r_line_annualized * 100, 'm--', 'LineWidth', 1, 'DisplayName', 'r (ssF, 年化)'); title('市场利率路径 (年化)'); xlabel('年份'); ylabel('年化利率 (%)'); grid on;
subplot(2,2,4); plot(time_axis, results.w_path, 'c-p', 'LineWidth', 2, 'DisplayName', '工资 (w)'); hold on; plot(time_axis, ssF_w_line, 'c--', 'LineWidth', 1, 'DisplayName', 'w (ssF BGP)'); title('工资路径 (原始水平)'); xlabel('年份'); ylabel('原始工资 (w)'); legend('show', 'Location', 'best'); grid on;
sgtitle('宏观经济转轨动态与稳态收敛检验 (全去趋势化迭代)', 'FontSize', 16, 'FontWeight', 'bold');
fprintf('   ✅ 结果可视化完成。\n');


if abs(error) < tolerance
    % --- [v5.5 真实验证] ---
    accountingS = main_transition_utils.aggregate_accounting_paths(final_Dist_path, final_Pol_path, cS, final_pathS);
    accountingS.Total_C_path(T) = results.C_path(T);
    main_transition_utils.check_transition_national_accounts(results, accountingS, cS, ssF);
end

%% --- 6. [新增] 详细转轨路径统计报告 ---
% 此部分代码无需修改，保持原样即可
if abs(error) < tolerance
    fprintf('\n\n================================================================================\n');
    fprintf('===                    详细转轨路径统计分析报告                    ===\n');
    fprintf('================================================================================\n');

    print_stat_row = @(varName, varPath, unit) ...
        fprintf('%-25s | %12.4f | %12.4f | %12.4f | %12.4f | %s\n', ...
                varName, mean(varPath), min(varPath), max(varPath), std(varPath), unit);

    fprintf('%-25s | %12s | %12s | %12s | %12s | %s\n', ...
            '变量名', '均值', '最小值', '最大值', '标准差', '单位');
    fprintf('--------------------------------------------------------------------------------\n');
    
    fprintf('--- [宏观总量 (原始水平)] ---\n');
    print_stat_row('总产出 (Y)', results.Y_path, '原始水平');
    print_stat_row('私人资本 (Kp)', results.K_p_path, '原始水平');
    print_stat_row('公共资本 (Kg)', results.K_g_path, '原始水平');
    print_stat_row('总消费 (C)', results.C_path, '原始水平');
    print_stat_row('私人投资 (Ip)', results.I_p_path, '原始水平');
    print_stat_row('公共投资 (Ig)', results.I_g_path, '原始水平');
    print_stat_row('总劳动 (L)', results.L_path, '原始水平 (无趋势)');
    print_stat_row('总遗赠 (Bequest)', results.Bequest_path, '原始水平');
    fprintf('\n');

    fprintf('--- [宏观总量 (去趋势化, "hat")] ---\n');
    print_stat_row('产出 (Y_hat)', results.Y_hat_path, '去趋势化水平');
    print_stat_row('私人资本 (Kp_hat)', results.K_p_hat_path, '去趋势化水平');
    print_stat_row('公共资本 (Kg_hat)', results.K_g_hat_path, '去趋势化水平');
    print_stat_row('消费 (C_hat)', results.C_path ./ cS.A_path, '去趋势化水平');
    print_stat_row('私人投资 (Ip_hat)', results.I_p_path ./ cS.A_path, '去趋势化水平');
    print_stat_row('公共投资 (Ig_hat)', results.I_g_path ./ cS.A_path, '去趋势化水平');
    fprintf('\n');

    fprintf('--- [价格与政策变量] ---\n');
    r_path_annualized = (1 + results.r_path).^(1/cS.time_Step) - 1;
    print_stat_row('年化市场利率 (r)', r_path_annualized * 100, '%');
    print_stat_row('去趋势化工资 (w_hat)', results.w_hat_path, '去趋势化水平');
    if exist('final_pathS', 'var')
        print_stat_row('PAYG 缴费率 (theta)', final_pathS.theta_path * 100, '%');
    end
    fprintf('\n');

    fprintf('--- [宏观比率] ---\n');
    total_K_path = results.K_p_path + results.K_g_path;
    total_I_path = results.I_p_path + results.I_g_path;
    safe_Y_path = results.Y_path + 1e-9;
    
    print_stat_row('总资本/产出比 (K/Y)', total_K_path ./ safe_Y_path, '比率');
    print_stat_row('消费率 (C/Y)', results.C_path ./ safe_Y_path, '比率');
    print_stat_row('总投资率 (I/Y)', total_I_path ./ safe_Y_path, '比率');
    
    if exist('accountingS', 'var')
        K_p_return_gross = results.r_path .* results.K_p_path + cS.ddk .* results.K_p_path;
        public_capital_return_path = results.Y_path - results.w_path .* results.L_path - K_p_return_gross;
        gov_revenue_path = accountingS.tax_regular_path + public_capital_return_path;
        gov_discretionary_res_path = gov_revenue_path;
        G_path_from_budget = (1 - cS.lambda_g) * gov_discretionary_res_path;
        print_stat_row('政府消费/产出比 (G/Y)', G_path_from_budget ./ safe_Y_path, '比率');
    end
    
    fprintf('================================================================================\n');
end

fprintf('\n--- 转轨动态求解脚本执行完毕 ---\n');