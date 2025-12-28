% =========================================================================
% == SCRIPT: analyze_equilibrium_difficulty.m
% == 版本: [v2.1 - 结果保存版]
% ==
% == 目的:
% ==   1. 延续 v2.0 的双稳态诊断功能。
% ==   2. [!!! 新增功能 !!!] 自动将生成的所有图表保存为 .png 文件到
% ==      一个名为 'error_analysis' 的子文件夹中。
% ==   3. 如果文件夹不存在，脚本会自动创建。
% ==
% == 如何使用:
% ==   1. 运行 main_run_SS.m 生成最新的稳态数据。
% ==   2. 直接运行此脚本。
% ==   3. 脚本执行完毕后，在当前目录下会生成一个 'error_analysis' 文件夹，
% ==      其中包含了四张命名清晰的诊断图。
% =========================================================================

clear; close all; clc;
addpath(pwd);
fprintf('=== 稳态均衡求解难度分析脚本 (v2.1 - 结果保存版) ===\n\n');

%% --- 1. 加载数据与创建输出文件夹 ---
fprintf('--- 1. 加载数据与创建输出文件夹 ---\n');
data_filename = 'SS/data_for_transition.mat';
if ~exist(data_filename, 'file'), error('稳态数据文件 "%s" 不存在。', data_filename); end
load(data_filename, 'data_for_transition');
fprintf('   ✅ 已加载数据: %s\n', data_filename);

% --- [新增] 创建输出文件夹 ---
output_folder = 'error_analysis';
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
    fprintf('   ✅ 已创建输出文件夹: %s\n', output_folder);
end

%% --- 2. [核心] 分别对初始和终期稳态进行分析 ---

% --- 分析初始稳态 (ss0) ---
fprintf('\n\n==========================================================\n');
fprintf('###              开始分析初始稳态 (ss0)              ###\n');
fprintf('==========================================================\n');
analyze_one_steadystate(data_for_transition.ss0, ...
                        data_for_transition.cS, ...
                        data_for_transition.paramS0, ...
                        data_for_transition.cS.Z_path(:, 1), ...
                        'ss0_initial', ... % 使用简短的文件名标签
                        output_folder);

% --- 分析终期稳态 (ssF) ---
fprintf('\n\n==========================================================\n');
fprintf('###              开始分析终期稳态 (ssF)              ###\n');
fprintf('==========================================================\n');
analyze_one_steadystate(data_for_transition.ssF, ...
                        data_for_transition.cS, ...
                        data_for_transition.paramSF, ...
                        data_for_transition.cS.Z_path(:, end), ...
                        'ssF_final', ... % 使用简短的文件名标签
                        output_folder);

fprintf('\n\n--- ✅ 分析脚本执行完毕 ---\n');


%% --- 辅助函数：执行单次稳态分析的核心逻辑 ---
function analyze_one_steadystate(ss, cS, paramS, Z_norm, file_label, output_folder)
    % --- 设定分析区间的参数 ---
    n_points_1D = 101; 
    n_points_2D = 41;  
    perturb_pct = 0.5; 
    
    if isempty(ss) || ~isstruct(ss) || ~isfield(ss, 'r_mkt')
        fprintf('   [警告] 传入的稳态结构体为空或无效，跳过对【%s】的分析。\n', file_label);
        return;
    end
    
    r_eq = ss.r_mkt;
    w_eq = ss.w_hat;
    mass_total = sum(Z_norm);
    bq_eq = ss.Bequest_distributed_agg / mass_total;
    tr_eq = ss.TR_distributed_agg / mass_total;
    b_hat_eq = ss.b_hat;

    fprintf('   分析中心点: r=%.4f, w=%.4f, bq=%.4f, tr=%.4f, b_hat=%.4f\n', r_eq, w_eq, bq_eq, tr_eq, b_hat_eq);
    fprintf('   分析区间: 将在均衡点附近 ±%.1f%% 的范围内进行扰动分析。\n', perturb_pct*100);

    % --- 一维敏感性分析：误差 vs. 利率 (r) ---
    fprintf('\n--- 正在进行一维敏感性分析 (Error vs. r) ---\n');
    r_vec = linspace(r_eq * (1 - perturb_pct), r_eq * (1 + perturb_pct), n_points_1D);
    r_vec = sort([r_vec, r_eq]); 
    n_points_1D = length(r_vec);
    errors_vs_r = zeros(6, n_points_1D);

    x_base = [NaN, w_eq, bq_eq, tr_eq, b_hat_eq]; 

    tic;
    parfor i = 1:n_points_1D
        x_guess_i = x_base;
        x_guess_i(1) = r_vec(i);
        F_error_i = SS.system_of_equations(x_guess_i, Z_norm, cS, paramS, struct('A',1.0));
        errors_vs_r(:, i) = F_error_i;
    end
    elapsed_time = toc;
    fprintf('   计算完成，耗时 %.2f 秒。\n', elapsed_time);

    fig1D = figure('Name', ['1D Sensitivity: ' file_label], 'Position', [100, 100, 1200, 800], 'Visible', 'off');
    error_labels = {'Capital Market (r)', 'Labor Market (w)', 'Bequest Market (BQ)', 'Gov Budget (TR)', 'PAYG System (b_{hat})', 'Capital Stock (K)'};
    colors = lines(6);

    for i = 1:6
        subplot(2, 3, i);
        plot(r_vec, errors_vs_r(i, :), 'LineWidth', 2, 'Color', colors(i,:));
        hold on;
        xline(r_eq, 'k--', 'LineWidth', 1.5);
        yline(0, 'k-', 'LineWidth', 0.5);
        title([error_labels{i} ' Error']);
        xlabel('Interest Rate (r)');
        ylabel('Error (Guess - Target)');
        grid on;
    end
    sgtitle(['1D Sensitivity: Market Clearing Errors vs. Interest Rate (r) @ ' file_label], 'FontSize', 14, 'FontWeight', 'bold');
    
    % [新增] 保存一维图
    filename_1D = fullfile(output_folder, ['sensitivity_1D_' file_label '.png']);
    saveas(fig1D, filename_1D);
    fprintf('   ✅ 一维敏感性图已保存至: %s\n', filename_1D);
    close(fig1D);

    % --- 二维误差曲面分析: 总误差 vs. (r, w) ---
    fprintf('\n--- 正在进行二维误差曲面分析 (Total Error vs. r, w) ---\n');
    r_grid_2D = linspace(r_eq * (1 - perturb_pct), r_eq * (1 + perturb_pct), n_points_2D);
    w_grid_2D = linspace(w_eq * (1 - perturb_pct), w_eq * (1 + perturb_pct), n_points_2D);
    [R_grid, W_grid] = meshgrid(r_grid_2D, w_grid_2D);
    Resnorm_surface = zeros(n_points_2D, n_points_2D);

    x_base_2D = [NaN, NaN, bq_eq, tr_eq, b_hat_eq];

    tic;
    x_base_2D_slice = x_base_2D;
    parfor i = 1:n_points_2D
        errors_col = zeros(n_points_2D, 1);
        for j = 1:n_points_2D
            x_guess_ij = x_base_2D_slice;
            x_guess_ij(1) = R_grid(j, i); 
            x_guess_ij(2) = W_grid(j, i);
            
            F_error_ij = SS.system_of_equations(x_guess_ij, Z_norm, cS, paramS, struct('A',1.0));
            errors_col(j) = sum(F_error_ij.^2);
        end
        Resnorm_surface(:, i) = errors_col;
    end
    elapsed_time = toc;
    fprintf('   计算完成，耗时 %.2f 秒。\n', elapsed_time);
    
    fig2D = figure('Name', ['2D Error Surface: ' file_label], 'Position', [200, 200, 1400, 700], 'Visible', 'off');
    
    subplot(1, 2, 1);
    surf(R_grid, W_grid, log10(Resnorm_surface));
    hold on;
    plot3(r_eq, w_eq, min(log10(Resnorm_surface(:))), 'r*', 'MarkerSize', 15, 'LineWidth', 3);
    xlabel('Interest Rate (r)'); ylabel('Wage (w)'); zlabel('log_{10}(Sum of Squared Errors)');
    title('Total Error Surface (Log Scale)'); colorbar; view(30, 45);

    subplot(1, 2, 2);
    contourf(R_grid, W_grid, log10(Resnorm_surface), 20);
    hold on;
    plot(r_eq, w_eq, 'r*', 'MarkerSize', 12, 'LineWidth', 2);
    xlabel('Interest Rate (r)'); ylabel('Wage (w)');
    title('Total Error Contour Map (Log Scale)'); colorbar; axis tight; grid on;

    sgtitle(['2D Error Surface: Total Error vs. (r, w) @ ' file_label], 'FontSize', 16, 'FontWeight', 'bold');
    
    % [新增] 保存二维图
    filename_2D = fullfile(output_folder, ['error_surface_2D_' file_label '.png']);
    saveas(fig2D, filename_2D);
    fprintf('   ✅ 二维误差曲面图已保存至: %s\n', filename_2D);
    close(fig2D);
end