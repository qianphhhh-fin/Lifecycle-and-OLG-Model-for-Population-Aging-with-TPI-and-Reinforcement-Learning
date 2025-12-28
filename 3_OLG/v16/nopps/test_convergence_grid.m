% =========================================================================
% == SCRIPT: main_run_convergence_test.m
% == 目的: 检验终期稳态(ssF)的国民账户闭合误差是否随网格密度增加而收敛
% =========================================================================
clear; close all; % ... (保留所有初始设置)

% --- [核心修改 1] 定义网格密度测试向量 ---
ngrid_vector = [40, 60, 80, 100, 150]; % 从粗到细的网格密度
results_table = table(); % 用于存储结果

fprintf('=== 开始收敛性检验 ===\n');

for i = 1:length(ngrid_vector)
    
    current_ngrid = ngrid_vector(i);
    fprintf('\n\n=========================================================\n');
    fprintf('=== 测试轮次 %d/%d: 网格密度 ngrid = %d ===\n', i, length(ngrid_vector), current_ngrid);
    fprintf('=========================================================\n');
    
    % --- [核心修改 2] 每次循环重新设定全局参数 ---
    
    % 1. 全局设置与初始化
    ngrid = current_ngrid;
    ngrid_pps = 1;
    cS = model_setup_utils_bgp.ParameterValues(); % 重新加载物理参数
    
    % 2. 生成对应密度的网格
    cS.nk = ngrid; 
    cS.nkpps = ngrid_pps; 
    cS.nkprime = ngrid; 
    cS.npps = ngrid_pps;
    cS = model_setup_utils_bgp.generateGrids(cS);
    fprintf('   ✅ 已为 ngrid=%d 重新生成全局网格。\n', current_ngrid);

    % 3. 保持其他所有参数不变 (养老金模式，年份等)
    cS.endogenous_theta_mode = true;
    cS.pps_active = false;
    cS.ss0_year = 2023;
    cS.start_year = 2023;
    if cS.endogenous_theta_mode
        cS.payg_replacement_rate = 0.40; % 保持与你报告中一致的参数
        cS.theta_max = 0.99;
    end
    
    % --- [核心修改 3] 只求解终期稳态(ssF)，并记录结果 ---
    
    % a. 生成外生路径 (这部分逻辑不变，但cS已更新)
    [Z_path, Z_path_raw, A_path, cS] = model_setup_utils_bgp.generate_exo_paths(cS, false);
    cS = model_setup_utils_bgp.calcaulte_theta_payg_path(cS, false);
    theta_path = cS.theta_path;

    % b. 为ssF设定环境 (逻辑不变)
    cSF = cS; 
    cSF.pps_active = false; 
    cSF.s_pathV = cS.s_pathV(:,end); 
    theta_for_ssF = []; % DB模式
    
    paramSF = struct();
    [paramSF.leGridV, paramSF.TrProbM_by_age, paramSF.leProb1V, cSF.nw_expanded] = model_setup_utils_bgp.EarningProcess_AgeDependent(cSF);
    paramSF.leLogGridV = log(paramSF.leGridV(1:cSF.nw));

    params_for_ssF = struct('Z', Z_path(:,end), 'A', 1.0, 'theta', theta_for_ssF, 'g_A_ss', cSF.g_A_ss, 'n_ss', cSF.n_ss);

    % c. 求解ssF (使用调整过容差的求解器)
    tic;
    fprintf('   ⚙️  启动统一稳态求解器 (求解 ssF, ngrid=%d)...\n', current_ngrid);
    [ssF_result, distF, polF, valF] = main_steady_state_utils_bgp.solve_steady_state_unified(cSF, paramSF, params_for_ssF, false, [], 'lsqnonlin'); % verbose=false, 减少屏幕输出
    solve_time = toc;

    if isempty(ssF_result)
        fprintf('   ❌ ngrid=%d 时求解失败！跳过此轮。\n', current_ngrid);
        continue;
    end
    
    % d. [关键] 计算并记录误差
    % 我们需要调用 display_national_accounts_unified 来获取误差结构体
    errorsF = main_steady_state_utils_bgp.display_national_accounts_unified(ssF_result, cSF, distF, polF, paramSF, false); % verbose=false
    
    % e. 将本轮结果存入table
    new_row = table(current_ngrid, errorsF.err_K_consistency, errorsF.err_NIPA, solve_time, ...
                    'VariableNames', {'GridPoints', 'Error_K_Consistency', 'Error_NIPA', 'SolveTime'});
    results_table = [results_table; new_row];
    
    fprintf('   ✅ ngrid=%d 测试完成。 NIPA误差 = %.4e\n', current_ngrid, errorsF.err_NIPA);
end

% --- [核心修改 4] 结果分析与可视化 ---
fprintf('\n\n===== 收敛性检验结果汇总 =====\n');
disp(results_table);

if height(results_table) > 1
    figure('Name', '离散化误差收敛性检验');
    
    % 使用双Y轴来展示误差和求解时间
    yyaxis left
    loglog(results_table.GridPoints, abs(results_table.Error_NIPA), 'bo-', 'LineWidth', 2, 'MarkerFaceColor', 'b');
    ylabel('国民账户闭合误差 |Y-C-I-G| (对数尺度)');
    ax = gca;
    ax.YAxis(1).Color = 'b';
    ax.YScale = 'log';
    ax.XScale = 'log';
    
    yyaxis right
    loglog(results_table.GridPoints, results_table.SolveTime, 'rs--', 'LineWidth', 1.5);
    ylabel('求解时间 (秒, 对数尺度)');
    ax.YAxis(2).Color = 'r';
    
    xlabel('资产网格点数 (对数尺度)');
    title('国民账户闭合误差 vs. 网格密度');
    legend('NIPA 误差', '求解时间', 'Location', 'best');
    grid on;
end