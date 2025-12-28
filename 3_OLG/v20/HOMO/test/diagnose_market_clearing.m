% =========================================================================
% == SCRIPT: diagnose_market_clearing.m
% == 版本: [v1.1 - 核心逻辑与并行修正版]
% ==
% == 目的:
% ==   1. 修正了 v1.0 中因传入 NaN 参数导致计算失败的根本问题。
% ==   2. 确保在分析一个市场的供求曲线时，其他市场的价格和转移支付
% ==      都固定在其全局均衡水平上，保证了分析的一致性。
% ==   3. 修正了劳动需求曲线的计算公式。
% ==   4. 为资本供给曲线和劳动供给曲线的计算启用了并行化(parfor)，
% ==      并添加了进度条，以大幅提升执行效率。
% =========================================================================

clear; close all;
addpath(pwd);
fprintf('=== 市场出清诊断脚本 (v1.1 - 核心逻辑与并行修正版) ===\n\n');

%% --- 1. 加载数据与创建输出文件夹 ---
fprintf('--- 1. 加载数据与创建输出文件夹 ---\n');
data_filename = 'SS/data_for_transition.mat';
if ~exist(data_filename, 'file'), error('稳态数据文件 "%s" 不存在。', data_filename); end
load(data_filename, 'data_for_transition');
fprintf('   ✅ 已加载数据: %s\n', data_filename);

output_folder = 'market_clearing_diagnosis';
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
    fprintf('   ✅ 已创建输出文件夹: %s\n', output_folder);
end

% 启动并行池
if isempty(gcp('nocreate')), parpool; end


%% --- 2. [核心] 分别对初始和终期稳态进行分析 ---
% --- 分析初始稳态 (ss0) ---
fprintf('\n\n==========================================================\n');
fprintf('###              开始诊断初始稳态 (ss0)              ###\n');
fprintf('==========================================================\n');
diagnose_one_steadystate(data_for_transition.ss0, ...
                         data_for_transition.cS, ...
                         data_for_transition.paramS0, ...
                         data_for_transition.cS.Z_path(:, 1), ...
                         'ss0_initial', output_folder);

% --- 分析终期稳态 (ssF) ---
fprintf('\n\n==========================================================\n');
fprintf('###              开始诊断终期稳态 (ssF)              ###\n');
fprintf('==========================================================\n');
diagnose_one_steadystate(data_for_transition.ssF, ...
                         data_for_transition.cS, ...
                         data_for_transition.paramSF, ...
                         data_for_transition.cS.Z_path(:, end), ...
                         'ssF_final', output_folder);

fprintf('\n\n--- ✅ 诊断脚本执行完毕 ---\n');


%% --- 辅助函数：执行单次稳态诊断的核心逻辑 ---
function diagnose_one_steadystate(ss, cS, paramS, Z_norm, file_label, output_folder)
    
    if isempty(ss) || ~isstruct(ss) || ~isfield(ss, 'r_mkt')
        fprintf('   [警告] 传入的稳态结构体为空或无效，跳过对【%s】的诊断。\n', file_label);
        return;
    end
    
    % --- 0. 准备均衡值和参数 ---
    r_eq = ss.r_mkt;
    w_eq = ss.w_hat;
    K_eq = ss.K_private_hat;
    L_eq = ss.L_hat;
    
    mass_total = sum(Z_norm);
    bq_eq = ss.Bequest_distributed_agg / mass_total;
    tr_eq = ss.TR_distributed_agg / mass_total;
    b_hat_eq = ss.b_hat;
    
    n_points = 21;
    perturb_pct = 0.5;

    % --- 1. 资本市场分析 (K_s(r) vs K_d(r)) ---
    fprintf('\n--- 正在分析资本市场 ---\n');
    r_vec = linspace(r_eq * (1 - perturb_pct), r_eq * (1 + perturb_pct), n_points);
    K_supply_vec_par = zeros(size(r_vec)); % Parfor 临时变量
    K_demand_vec = zeros(size(r_vec));

    % a. 计算资本供给曲线 K_s(r)
    fprintf('   计算资本供给曲线 K_s(r) (启用并行计算)...\n');
    tic;
    parfor_progress(n_points);
    parfor i = 1:n_points
        % [!!! 核心修正 !!!] 
        % 根据正确的函数签名，直接捕获第二个输出参数 `dist_i`
        [~, dist_i, ~, ~] = SS.run_micro_to_macro(r_vec(i), w_eq, bq_eq, tr_eq, b_hat_eq, Z_norm, cS, paramS, struct('A',1.0));
        
        % 直接从返回的分布中聚合当期资本存量
        aggrS_i = aggregates.get_aggregates_from_dist(dist_i, cS);
        K_supply_vec_par(i) = aggrS_i.K_p_hat_total;
        
        parfor_progress;
    end
    K_supply_vec = K_supply_vec_par; % 将结果从临时变量赋回
    parfor_progress(0);
    toc;

    % b. 计算资本需求曲线 K_d(r)
    g_A_period = (1 + cS.g_A_ss)^cS.time_Step - 1;
    n_period = (1 + cS.n_ss)^cS.time_Step - 1;
    g_total_period = (1 + g_A_period) * (1 + n_period) - 1;
    Kg_Y_ratio = cS.I_g_to_Y_ratio_ss / (g_total_period + cS.ddk_g);

    for i = 1:n_points
        target_r = r_vec(i);
        f_k_wrapper = @(K) (K > 0) .* (cS.alpha * ( (K.^cS.alpha * L_eq^(1-cS.alpha-cS.gamma) * Kg_Y_ratio^cS.gamma).^(1/(1-cS.gamma)) ) ./ K - (target_r + cS.ddk));
        search_interval_K = [1e-8, K_eq * 10];
        
        try
            K_demand_vec(i) = fzero(f_k_wrapper, search_interval_K);
        catch
            K_demand_vec(i) = NaN; 
        end
    end

    % c. 绘制资本市场图
    fig_K = figure('Name', ['Capital Market: ' file_label], 'Visible', 'on');
    plot(K_supply_vec, r_vec, 'b-o', 'LineWidth', 2, 'DisplayName', '资本供给 K_s(r)');
    hold on;
    plot(K_demand_vec, r_vec, 'r-d', 'LineWidth', 2, 'DisplayName', '资本需求 K_d(r)');
    plot(K_eq, r_eq, 'k*', 'MarkerSize', 12, 'LineWidth', 2, 'DisplayName', '均衡点');
    xlabel('资本量 (K)'); ylabel('利率 (r)');
    title(['资本市场出清 @ ' file_label]);
    legend('show', 'location', 'best'); grid on;
    
    filename_K = fullfile(output_folder, ['market_capital_' file_label '.png']);
    saveas(fig_K, filename_K);
    fprintf('   ✅ 资本市场诊断图已保存至: %s\n', filename_K);
    close(fig_K);


    % --- 2. 劳动市场分析 (L_s(w) vs L_d(w)) ---
    fprintf('\n--- 正在分析劳动市场 ---\n');
    w_vec = linspace(w_eq * (1 - perturb_pct), w_eq * (1 + perturb_pct), n_points);
    L_supply_vec_par = zeros(size(w_vec)); 
    L_demand_vec = zeros(size(w_vec));

    % a. 计算劳动供给曲线 L_s(w)
    fprintf('   计算劳动供给曲线 L_s(w) (启用并行计算)...\n');
    tic;
    parfor_progress(n_points);
    parfor i = 1:n_points
        % [!!! 核心修正 !!!] 
        % 捕获第一个输出 `macro_state`
        [macro_state, ~, ~, ~] = SS.run_micro_to_macro(r_eq, w_vec(i), bq_eq, tr_eq, b_hat_eq, Z_norm, cS, paramS, struct('A',1.0));
        L_supply_vec_par(i) = macro_state.L_hat;
        parfor_progress;
    end
    parfor_progress(0);
    L_supply_vec = L_supply_vec_par;
    toc;
    
    % b. 计算劳动需求曲线 L_d(w)
    labor_share = 1 - cS.alpha - cS.gamma;
    for i = 1:n_points
        target_w = w_vec(i);
        f_l_wrapper = @(L) (L > 0) .* (labor_share * ( (K_eq^cS.alpha * L.^(1-cS.alpha-cS.gamma) * Kg_Y_ratio^cS.gamma).^(1/(1-cS.gamma)) ) ./ L - target_w);
        search_interval_L = [1e-8, L_eq * 10];
        
        try
            L_demand_vec(i) = fzero(f_l_wrapper, search_interval_L);
        catch
            L_demand_vec(i) = NaN;
        end
    end

    % c. 绘制劳动市场图
    fig_L = figure('Name', ['Labor Market: ' file_label], 'Visible', 'on');
    plot(L_supply_vec, w_vec, 'b-o', 'LineWidth', 2, 'DisplayName', '劳动供给 L_s(w)');
    hold on;
    plot(L_demand_vec, w_vec, 'r-d', 'LineWidth', 2, 'DisplayName', '劳动需求 L_d(w)');
    plot(L_eq, w_eq, 'k*', 'MarkerSize', 12, 'LineWidth', 2, 'DisplayName', '均衡点');
    xlabel('劳动量 (L)'); ylabel('工资 (w)');
    title(['劳动市场出清 @ ' file_label]);
    legend('show', 'location', 'best'); grid on;

    filename_L = fullfile(output_folder, ['market_labor_' file_label '.png']);
    saveas(fig_L, filename_L);
    fprintf('   ✅ 劳动市场诊断图已保存至: %s\n', filename_L);
    close(fig_L);

end