% =========================================================================
% == SCRIPT: analyze_pps_structural_reform_welfare.m
% == 版本: [v1.1 - 缓存中间结果版]
% ==
% == 目的:
% ==   1. 比较'structural_reform'与'moderate'两种PPS发展路径下的福利差异。
% ==   2. 通过重构价值函数路径，计算各年龄、各类型家庭的补偿变化(CV)。
% ==   3. 生成LaTeX表格，以揭示结构性改革的福利代价在不同群体间的分配。
% ==
% == v1.1 核心修改:
% ==   - [!!!] 新增了缓存机制。脚本会检查是否存在包含已计算价值函数的
% ==     中间文件。如果存在，则跳过耗时的反向迭代计算，直接加载结果。
% ==   - 如果不存在中间文件，则执行计算并将结果保存，供后续运行使用。
% =========================================================================

clear; close all; clc;
addpath(pwd);
fprintf('=== PPS结构性改革福利分析脚本 (v1.1) ===\n\n');

%% --- 1. 用户设定与文件定义 ---
fprintf('--- 1. 用户设定与文件定义 ---\n');

% --- 分析设定 ---
reform_scenario = 'structural_reform';
base_scenario = 'moderate'; % 福利分析的基准情景

% --- 文件路径模板 ---
file_tpi_pattern = 'TRANS/TPI_results_het_pps_%s.mat';

% --- [核心修改] 中间结果缓存文件 ---
intermediate_results_file = 'results/welfare_analysis_pps_reform_intermediate.mat';

% --- 输出路径 ---
output_dir_tex = 'tex/tab/';
if ~exist(output_dir_tex, 'dir'), mkdir(output_dir_tex); end
output_tex_table_welfare = fullfile(output_dir_tex, 'tab4.6_0_welfare_pps_structural_reform.tex');

fprintf('   改革方案: %s\n', reform_scenario);
fprintf('   基准方案: %s\n', base_scenario);
fprintf('   LaTeX 表格将输出至: %s\n', output_tex_table_welfare);


%% --- 2. [数据加载] 加载所有情景的计算结果 ---
fprintf('\n--- 2. [数据加载] 循环加载方案结果 ---\n');

results_all = struct();
scenarios_to_load = {base_scenario, reform_scenario};

for i = 1:length(scenarios_to_load)
    scenario_name = scenarios_to_load{i};
    fprintf('   正在加载方案: [%s] ...\n', upper(scenario_name));

    file_tpi = sprintf(file_tpi_pattern, scenario_name);
    if ~exist(file_tpi, 'file')
        error('找不到文件: %s，请确认main_run_trans_pps.m已针对该情景运行。', file_tpi); 
    end

    results_all.(scenario_name) = load(file_tpi);
    fprintf('   ✅ 已加载: %s\n', file_tpi);
end
fprintf('--- 所有情景数据加载完毕 ---\n');


%% --- 3. [福利分析] 为各情景重构价值函数路径 (或从缓存加载) ---
fprintf('\n--- 3. [福利分析] 重构价值函数路径 ---\n');

skip_computation = false;
if exist(intermediate_results_file, 'file')
    fprintf('   发现中间结果文件: %s\n', intermediate_results_file);
    try
        loaded_data = load(intermediate_results_file, 'saved_val_paths');
        if isfield(loaded_data.saved_val_paths, base_scenario) && isfield(loaded_data.saved_val_paths, reform_scenario)
            fprintf('   ✅ 文件内容有效，正在加载已保存的价值函数路径...\n');
            results_all.(base_scenario).Val_path_h = loaded_data.saved_val_paths.(base_scenario);
            results_all.(reform_scenario).Val_path_h = loaded_data.saved_val_paths.(reform_scenario);
            skip_computation = true;
        else
            fprintf('   ⚠️ 文件内容不匹配当前分析的方案，将重新计算。\n');
        end
    catch ME
        fprintf('   ❌ 加载失败: %s。将重新计算。\n', ME.message);
    end
else
    fprintf('   未发现中间结果文件，将执行完整计算。\n');
end


if ~skip_computation
    for i_scen = 1:length(scenarios_to_load)
        scen_name = scenarios_to_load{i_scen};
        fprintf('   正在为方案 [%s] 运行反向迭代...\n', upper(scen_name));
        
        data = results_all.(scen_name);
        cS = data.cS;
        paramSF = data.results.ss_data.paramSF;
        valF_h = data.results.ss_data.valF_h;
        polF_h = data.results.ss_data.polF_h;
        iter_results = data.results.iter_results;
        
        Val_path_h = cell(cS.nTypes, 1);
        
        for h = 1:cS.nTypes
            pathS_h = get_pathS_for_backward(iter_results, h, cS, data.ss0);
            cS_h = get_cs_for_type(cS, h);
            [~, Val_path_h{h}] = household.backward_hh(pathS_h, cS_h, paramSF, valF_h{h}, polF_h{h});
        end
        
        results_all.(scen_name).Val_path_h = Val_path_h;
        fprintf('   ✅ 方案 [%s] 价值函数重构完成。\n', upper(scen_name));
    end

    % [核心修改] 保存计算结果到缓存文件
    fprintf('   正在保存中间计算结果至: %s ...\n', intermediate_results_file);
    saved_val_paths = struct();
    saved_val_paths.(base_scenario) = results_all.(base_scenario).Val_path_h;
    saved_val_paths.(reform_scenario) = results_all.(reform_scenario).Val_path_h;
    save(intermediate_results_file, 'saved_val_paths', '-v7.3');
    fprintf('   ✅ 中间结果保存成功。\n');
end

fprintf('--- 所有情景价值函数路径准备完毕 ---\n');


%% --- 4. [生成表格] 计算并生成福利分析(CV)表格 ---
fprintf('\n--- 4. [生成表格] 计算并生成福利分析(CV)表格 ---\n');

% --- 定义分析维度 ---
assessment_years_welfare = [2030, 2040, 2050, 2070, 2100];
type_labels = {'高收入', '中低收入'};
sector_labels = {'A. 城镇职工', 'B. 城乡居民'};

% --- 获取基准情景的数据 ---
data_base = results_all.(base_scenario);
data_reform = results_all.(reform_scenario);
cS_base = data_base.cS;
sigma = cS_base.sigma;
T_sim_base = cS_base.T_sim;

% --- 定义年龄组 ---
age_groups = {
    '青年 (25-29岁)', round((25 - cS_base.age1_orig)/cS_base.time_Step) + 1;
    '中年 (45-49岁)', round((45 - cS_base.age1_orig)/cS_base.time_Step) + 1;
    '临退 (60-64岁)', cS_base.aR_new_path(1);
    '老年 (70-74岁)', round((70 - cS_base.age1_orig)/cS_base.time_Step) + 1;
};

% --- 打开文件并写入LaTeX表头 ---
fileID_welfare = fopen(output_tex_table_welfare, 'w', 'n', 'UTF-8');
fprintf(fileID_welfare, '%% =======================================================\n');
fprintf(fileID_welfare, '%%  此文件由 analyze_pps_structural_reform_welfare.m 自动生成\n');
fprintf(fileID_welfare, '%% =======================================================\n');
fprintf(fileID_welfare, '\\begin{table}[h!]\n');
fprintf(fileID_welfare, '\\centering\n');
fprintf(fileID_welfare, '\\caption{结构性养老金改革的福利影响 (相对于温和改革情景的补偿变化 CV, \\%%)}\n');
fprintf(fileID_welfare, '\\label{tab:welfare_impact_pps_reform}\n');
fprintf(fileID_welfare, '\\begin{threeparttable}\n');

header_cols = ['l', repmat('c', 1, length(assessment_years_welfare))];
fprintf(fileID_welfare, '\\begin{tabular}{%s}\n', header_cols);
fprintf(fileID_welfare, '\\toprule\n');

header_line = '家庭类型与年龄组';
for year = assessment_years_welfare
    header_line = [header_line, sprintf(' & %d年', year)];
end
fprintf(fileID_welfare, '%s \\\\\n', header_line);
fprintf(fileID_welfare, '\\midrule\n');

% --- 核心CV计算与写入 ---
for i_sector = 1:length(sector_labels)
    fprintf(fileID_welfare, '\\multicolumn{%d}{l}{\\textbf{%s}} \\\\\n', 1 + length(assessment_years_welfare), sector_labels{i_sector});
    for i_type = 1:length(type_labels)
        h_idx = (i_sector - 1) * 2 + i_type;
        
        fprintf(fileID_welfare, '    \\hspace{1em}%s', type_labels{i_type});
        fprintf(fileID_welfare, repmat(' & ', 1, length(assessment_years_welfare)));
        fprintf(fileID_welfare, ' \\\\\n');

        for i_age = 1:size(age_groups, 1)
            age_label = age_groups{i_age, 1};
            a_idx = age_groups{i_age, 2};
            
            fprintf(fileID_welfare, '       \\hspace{2em}%s', age_label);

            for i_yr = 1:length(assessment_years_welfare)
                year = assessment_years_welfare(i_yr);
                t_idx = round((year - cS_base.start_year) / cS_base.time_Step) + 1;
                
                cv_str = '-';
                if t_idx <= T_sim_base
                    dist_base_slice = data_base.final_Dist_path_h{h_idx}(:, :, :, a_idx, t_idx);
                    val_base_slice  = data_base.Val_path_h{h_idx}(:, :, :, a_idx, t_idx);
                    val_reform_slice= data_reform.Val_path_h{h_idx}(:, :, :, a_idx, t_idx);
                    
                    mass_slice = sum(dist_base_slice, 'all');
                    
                    cv = NaN;
                    if mass_slice > 1e-12
                        EV_base = sum(val_base_slice .* dist_base_slice, 'all') / mass_slice;
                        EV_reform = sum(val_reform_slice .* dist_base_slice, 'all') / mass_slice;
                        
                        if EV_base < -1e-9 && EV_reform < -1e-9
                            cv = (EV_reform / EV_base)^(1 / (1 - sigma)) - 1;
                        end
                    end
                    cv_str = format_cv(cv);
                end
                fprintf(fileID_welfare, ' & %s', cv_str);
            end % year loop
            fprintf(fileID_welfare, ' \\\\\n');
        end % age group loop
    end % type loop
    if i_sector == 1, fprintf(fileID_welfare, '\\midrule\n'); end
end % sector loop

% --- 计算并写入总体平均 ---
fprintf(fileID_welfare, '\\midrule\n');
fprintf(fileID_welfare, '\\multicolumn{%d}{l}{\\textbf{C. 总体平均}} \\\\\n', 1 + length(assessment_years_welfare));
fprintf(fileID_welfare, '    人口加权平均 CV');

for i_yr = 1:length(assessment_years_welfare)
    year = assessment_years_welfare(i_yr);
    t_idx = round((year - cS_base.start_year) / cS_base.time_Step) + 1;
    
    cv_str = '-';
    if t_idx <= T_sim_base
        total_pop_t = sum(cS_base.Z_path_raw(:, t_idx), 'all');
        total_weighted_ev_base = 0;
        total_weighted_ev_reform = 0;
        
        for h_idx = 1:cS_base.nTypes
            for a_idx = 1:cS_base.aD_new
                dist_base_slice = data_base.final_Dist_path_h{h_idx}(:, :, :, a_idx, t_idx);
                val_base_slice  = data_base.Val_path_h{h_idx}(:, :, :, a_idx, t_idx);
                val_reform_slice= data_reform.Val_path_h{h_idx}(:, :, :, a_idx, t_idx);
                
                total_weighted_ev_base   = total_weighted_ev_base   + sum(val_base_slice .* dist_base_slice, 'all');
                total_weighted_ev_reform = total_weighted_ev_reform + sum(val_reform_slice .* dist_base_slice, 'all');
            end
        end
        
        EV_base_total_avg = total_weighted_ev_base / total_pop_t;
        EV_reform_total_avg = total_weighted_ev_reform / total_pop_t;
        
        cv_total = NaN;
        if EV_base_total_avg < -1e-9 && EV_reform_total_avg < -1e-9
             cv_total = (EV_reform_total_avg / EV_base_total_avg)^(1 / (1 - sigma)) - 1;
        end
        cv_str = format_cv(cv_total);
    end
    fprintf(fileID_welfare, ' & %s', cv_str);
end
fprintf(fileID_welfare, ' \\\\\n');


% --- 写入表尾 ---
fprintf(fileID_welfare, '\\bottomrule\n');
fprintf(fileID_welfare, '\\end{tabular}\n');
fprintf(fileID_welfare, '\\begin{tablenotes}[para,flushleft]\n');
fprintf(fileID_welfare, '  \\item 注：表中数值为补偿变化(CV)，衡量“结构性改革”情景相对于“温和改革”情景的福利差异。负值代表福利损失。例如，-0.5\\%%表示该群体愿意放弃其在温和改革情景下剩余生命周期消费的0.5\\%%，以避免进入结构性改革情景。\n');
fprintf(fileID_welfare, '\\end{tablenotes}\n');
fprintf(fileID_welfare, '\\end{threeparttable}\n');
fprintf(fileID_welfare, '\\end{table}\n');
fclose(fileID_welfare);
fprintf('   ✅ 成功生成 LaTeX 表格文件: %s\n', output_tex_table_welfare);

fprintf('\n--- 脚本执行完毕 ---\n');


%% --- 5. [生成图表] 比较两种情景下的宏观经济路径 ---
fprintf('\n--- 5. [生成图表] 比较两种情景下的宏观经济路径 ---\n');

time_axis = data_base.cS.start_year : data_base.cS.time_Step : (data_base.cS.start_year + (data_base.cS.T_sim-1)*data_base.cS.time_Step);
time_mask = (time_axis >= 2023) & (time_axis <= 2100);
time_axis_plot = time_axis(time_mask);


% --- 定义输出路径 ---
output_dir_fig = 'tex/fig/';
if ~exist(output_dir_fig, 'dir'), mkdir(output_dir_fig); end
output_fig_macro_comparison = fullfile(output_dir_fig, 'fig4.8-1_macro_comparison_pps_reforms.png');

% --- 准备绘图数据 ---
scenarios_to_plot = {base_scenario, reform_scenario};
plot_data = struct();

for i = 1:length(scenarios_to_plot)
    scen_name = scenarios_to_plot{i};
    data = results_all.(scen_name);
    
    % 利率 (年化, %)
    plot_data.(scen_name).r_path = ((1 + data.results.r_path).^(1/data.cS.time_Step) - 1) * 100;
    
    % 工资 (level)
    plot_data.(scen_name).w_path = data.results.w_path;

    % K/Y (私人资本/产出)
    plot_data.(scen_name).Kp_to_Y_path = data.results.K_p_path ./ data.results.Y_path;

    % C/Y (消费/产出)
    plot_data.(scen_name).C_to_Y_path = data.results.C_path ./ data.results.Y_path;
end

% --- 绘图风格设定 ---
fig_macro = figure('Name', '宏观经济路径比较：温和改革 vs. 结构性改革');
fig_macro.Units = 'centimeters';
fig_macro.Position = [15 15 14.99 5.5]; % [left bottom width height]

tcl_macro = tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');


% 灰度风格定义
style_base = {'-', 'Color', [0.5 0.5 0.5], 'LineWidth', 2.5};
style_reform = {'--', 'Color', [0.0 0.0 0.0], 'LineWidth', 2.5};
legend_labels = {'温和改革 (基准)', '结构性改革'};

% --- 子图1: 利率路径 ---
ax1 = nexttile;
hold(ax1, 'on');
plot(ax1, time_axis_plot, plot_data.(base_scenario).r_path(time_mask), style_base{:});
plot(ax1, time_axis_plot, plot_data.(reform_scenario).r_path(time_mask), style_reform{:});
hold(ax1, 'off');
grid(ax1, 'on'); box(ax1, 'on');
xlim(ax1, [2023, 2100]);
title(ax1, '(a) 年化真实利率(%)', 'FontSize', 10);
xlabel(ax1, '年份');
l1 = legend(ax1, legend_labels, 'position',[0.0765 0.1992 0.1684 0.1562],'IconColumnWidth',10,'Box','off');

% --- 子图2: 工资路径 ---
ax2 = nexttile;
hold(ax2, 'on');
plot(ax2, time_axis_plot, plot_data.(base_scenario).w_path(time_mask), style_base{:});
plot(ax2, time_axis_plot, plot_data.(reform_scenario).w_path(time_mask), style_reform{:});
hold(ax2, 'off');
grid(ax2, 'on'); box(ax2, 'on');
xlim(ax2, [2023, 2100]);
title(ax2, '(b) 工资水平', 'FontSize', 10);
xlabel(ax2, '年份');
l2 = legend(ax2, legend_labels, 'position',[0.3862 0.7083 0.1684 0.1562],'IconColumnWidth',10,'Box','off');

% --- 子图3: 宏观比率 (K/Y 和 C/Y) ---
ax3 = nexttile;
hold(ax3, 'on');
% 绘制 K/Y
p1 = plot(ax3, time_axis_plot, plot_data.(base_scenario).Kp_to_Y_path(time_mask), '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 2.5);
p2 = plot(ax3, time_axis_plot, plot_data.(reform_scenario).Kp_to_Y_path(time_mask), '--', 'Color', [0.0 0.0 0.0], 'LineWidth', 2.5);
% 绘制 C/Y
p3 = plot(ax3, time_axis_plot, plot_data.(base_scenario).C_to_Y_path(time_mask), '-.', 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
p4 = plot(ax3, time_axis_plot, plot_data.(reform_scenario).C_to_Y_path(time_mask), ':', 'Color', [0.0 0.0 0.0], 'LineWidth', 2);
hold(ax3, 'off');
grid(ax3, 'on'); box(ax3, 'on');
xlim(ax3, [2023, 2100]);
title(ax3, '(c) 宏观比率', 'FontSize', 10);
xlabel(ax3, '年份');
l3 = legend([p1 p2 p3 p4], {'K/Y (温和)', 'K/Y (结构性)', 'C/Y (温和)', 'C/Y (结构性)'}, 'position',[0.7013 0.5861 0.1455 0.2957],'IconColumnWidth',10,'Box','off');

%% --- 6. 保存新生成的图表 ---
fprintf('\n--- 6. 正在保存宏观比较图表 ---\n');
try
    exportgraphics(fig_macro, output_fig_macro_comparison, 'Resolution', 300);
    fprintf('   ✅ 成功保存图表至: %s\n', output_fig_macro_comparison);
catch ME
    fprintf('   ❌ 保存失败: %s\n', ME.message);
end

fprintf('\n--- 脚本执行完毕 ---\n');

%% --- [本地辅助函数] ---

function pathS_h = get_pathS_for_backward(iter_results, h, cS, ss0)
    % 为指定的家庭类型 h，打包后向迭代所需的所有时间路径
    pathS_h = struct();
    pathS_h.r_path = iter_results.r_path_iter;
    w_hat_path_converged = iter_results.w_hat_path_iter;
    pathS_h.w_hat_path = w_hat_path_converged;
    pathS_h.tr_per_hh_hat_path = iter_results.TR_hat_pc_path_iter;
    
    w_hat_path_for_b = [ss0.w_hat, w_hat_path_converged(1:end-1)];
    b_hat_formula_path_h_all = zeros(cS.nTypes, cS.T_sim);
    for t = 1:cS.T_sim
        b_hat_formula_path_h_all(:, t) = SS.calculate_formula_benefits(w_hat_path_for_b(t), cS);
    end
    adj_path_converged = iter_results.adj_path_iter;
    pathS_h.b_hat_path = adj_path_converged .* b_hat_formula_path_h_all(h, :);
    
    A_path_ext = [cS.A_path, cS.A_path(end) * (1 + ((1+cS.g_A_ss)^cS.time_Step-1))];
    g_A_path = A_path_ext(2:end)./A_path_ext(1:end-1)-1;
    pathS_h.g_A_path = g_A_path(1:cS.T_sim);
    pathS_h.theta_path = cS.theta_path_h(h, :);
end

function cS_h = get_cs_for_type(cS, h)
    % 为指定的家庭类型 h，准备一个专属的 cS 参数包
    cS_h = cS;
    cS_h.ageEffV_new = cS.ageEffV_new_h(:, h);
end

function str_out = format_cv(cv_value)
    % 辅助函数: 将CV数值格式化为带正负号和百分号的字符串
    if isnan(cv_value)
        str_out = '-';
        return;
    end
    
    if cv_value >= 0
        sign_char = '+';
    else
        sign_char = ''; % 负号会自动显示
    end
    str_out = sprintf('%s%.2f\\%%', sign_char, cv_value * 100);
end