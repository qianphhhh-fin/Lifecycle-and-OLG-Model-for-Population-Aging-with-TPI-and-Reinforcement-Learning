% =========================================================================
% == SCRIPT: welfare_analysis.m
% == 版本: [v3.0 - 异质性家庭分解版]
% ==
% == 目的:
% ==   - 基于一个包含四类异质性家庭的模型，进行福利分析。
% ==   - 核心任务 1: 生成一张LaTeX表格，按家庭类型和总体平均，
% ==     展示在不同截面年份下，所有年龄组加权平均的福利变化(CV)。
% ==   - 核心任务 2: 生成一张2x2的图表，分类型展示青年、中年、
% ==     临近退休三个年龄组的CV随时间的变化路径。
% ==
% == 前置条件:
% ==   - 已运行 `main_run_trans.m` 两次，分别生成了带与不带PPS的转轨结果。
% ==   - TRANS/TPI_results_het_nopps.mat & TPI_results_het_pps.mat
% ==   - TRANS/iter_results_het_nopps.mat & iter_results_het_pps.mat
% ==   - SS/data_for_het_transition_nopps.mat & SS/data_for_het_transition_pps.mat
% =========================================================================

clear; close all; clc;
addpath(pwd);
fprintf('=== OLG异质性模型福利分析脚本 (v3.0 - 类型分解版) ===\n\n');

%% --- 1. 设置与加载数据 ---
fprintf('--- 1. 设置与加载数据 ---\n');

% --- 文件路径定义 ---
file_path_nopps = 'TRANS/TPI_results_het_nopps.mat';
file_path_pps = 'TRANS/TPI_results_het_pps.mat';
iter_file_nopps = 'TRANS/iter_results_het_nopps.mat';
iter_file_pps = 'TRANS/iter_results_het_pps.mat';
ss_data_nopps_path = 'SS/data_for_het_transition_nopps.mat';
ss_data_pps_path = 'SS/data_for_het_transition_pps.mat';
output_tex_filename = 'tex/tab/tab4.1_welfare_analysis_by_type.tex';
output_fig_filename = 'tex/fig/fig4.3_welfare_analysis_detailed.png';

% --- 加载数据 ---
fprintf('   正在加载 "无PPS" (基准) 和 "有PPS" (改革) 情景数据...\n');
data_nopps = load(file_path_nopps);
iter_data_nopps = load(iter_file_nopps);
ss_data_nopps = load(ss_data_nopps_path);

data_pps = load(file_path_pps);
iter_data_pps = load(iter_file_pps);
ss_data_pps = load(ss_data_pps_path);
fprintf('   ✅ 数据加载完成。\n');


%% --- 2. [核心异质性修改] 重构价值函数路径 (分类型) ---
fprintf('\n--- 2. 重构价值函数路径 (分类型) ---\n');

nH = data_nopps.cS.nTypes;
Val_path_nopps_h = cell(nH, 1);
Val_path_pps_h = cell(nH, 1);

% --- a. 重构 "无PPS" 情景的价值函数路径 (循环所有类型) ---
fprintf('   正在为 "无PPS" 情景运行反向迭代...\n');
cS_nopps = data_nopps.cS;
paramSF_nopps = ss_data_nopps.data_for_transition.paramSF;
for h = 1:nH
    pathS_nopps_h = struct();
    pathS_nopps_h.r_path = data_nopps.results.r_path;
    pathS_nopps_h.w_hat_path = data_nopps.results.w_path ./ cS_nopps.A_path;
    pathS_nopps_h.b_hat_path = iter_data_nopps.b_hat_path_h_iter(h,:);
    pathS_nopps_h.tr_per_hh_hat_path = iter_data_nopps.TR_hat_pc_path_iter;
    A_path_ext_nopps = [cS_nopps.A_path, cS_nopps.A_path(end) * (1 + ((1+cS_nopps.g_A_ss)^cS_nopps.time_Step-1))];
    pathS_nopps_h.g_A_path = A_path_ext_nopps(2:end)./A_path_ext_nopps(1:end-1)-1;
    pathS_nopps_h.theta_path = cS_nopps.theta_path_h(h, :);
    
    cS_h_nopps = cS_nopps;
    cS_h_nopps.ageEffV_new = cS_nopps.ageEffV_new_h(:,h);

    valF_nopps_h = ss_data_nopps.data_for_transition.valF_h{h};
    polF_nopps_h = ss_data_nopps.data_for_transition.polF_h{h};
    
    [~, Val_path_nopps_h{h}] = household.backward_hh(pathS_nopps_h, cS_h_nopps, paramSF_nopps, valF_nopps_h, polF_nopps_h);
end

% --- b. 重构 "有PPS" 情景的价值函数路径 (循环所有类型) ---
fprintf('   正在为 "有PPS" 情景运行反向迭代...\n');
cS_pps = data_pps.cS;
paramSF_pps = ss_data_pps.data_for_transition.paramSF;
for h = 1:nH
    pathS_pps_h = struct();
    pathS_pps_h.r_path = data_pps.results.r_path;
    pathS_pps_h.w_hat_path = data_pps.results.w_path ./ cS_pps.A_path;
    pathS_pps_h.b_hat_path = iter_data_pps.b_hat_path_h_iter(h,:);
    pathS_pps_h.tr_per_hh_hat_path = iter_data_pps.TR_hat_pc_path_iter;
    A_path_ext_pps = [cS_pps.A_path, cS_pps.A_path(end) * (1 + ((1+cS_pps.g_A_ss)^cS_pps.time_Step-1))];
    pathS_pps_h.g_A_path = A_path_ext_pps(2:end)./A_path_ext_pps(1:end-1)-1;
    pathS_pps_h.theta_path = cS_pps.theta_path_h(h, :);

    cS_h_pps = cS_pps;
    cS_h_pps.ageEffV_new = cS_pps.ageEffV_new_h(:,h);

    valF_pps_h = ss_data_pps.data_for_transition.valF_h{h};
    polF_pps_h = ss_data_pps.data_for_transition.polF_h{h};

    [~, Val_path_pps_h{h}] = household.backward_hh(pathS_pps_h, cS_h_pps, paramSF_pps, valF_pps_h, polF_pps_h);
end

fprintf('   ✅ 价值函数路径重构完成。\n');


%% --- 3. [核心] 计算福利变化 (CV) ---
fprintf('\n--- 3. 计算各群体在特定年份的福利变化 ---\n');

% --- 定义分析维度 ---
assessment_years = [2023, 2030, 2040, 2050, 2070, 2100, 2150, 2200, 2250];
type_labels = {'高收入城镇职工'; '中低收入城镇职工'; '高收入居民'; '中低收入居民'};

% --- 定义绘图所需的代表性年龄组 (基于模型年龄索引 a_idx) ---
aR_idx = cS_nopps.aR_new;
age_groups_for_figure = {
    '青年 (25-29岁)', round((25 - cS_nopps.age1_orig)/cS_nopps.time_Step) + 1;
    '中年 (45-49岁)', round((45 - cS_nopps.age1_orig)/cS_nopps.time_Step) + 1;
    '临近退休 (60-64岁)', aR_idx
};
num_age_groups_fig = size(age_groups_for_figure, 1);

% --- 初始化结果存储 ---
num_years = length(assessment_years);
cv_table_data = zeros(nH, num_years); % 表格数据: 4种类型 x N个年份
cv_figure_data = zeros(num_age_groups_fig, num_years, nH); % 绘图数据: 3个年龄组 x N个年份 x 4种类型

% --- 核心计算循环 ---
sigma = cS_nopps.sigma;

for i_yr = 1:num_years
    year = assessment_years(i_yr);
    t_idx = round((year - cS_nopps.start_year) / cS_nopps.time_Step) + 1;
    
    fprintf('   正在评估年份: %d (t=%d)...\n', year, t_idx);
    
    if t_idx > cS_nopps.T_sim
        fprintf('      警告: 评估年份 %d (t=%d) 超出模拟范围，此列将留空。\n', year, t_idx);
        cv_table_data(:, i_yr) = NaN;
        cv_figure_data(:, i_yr, :) = NaN;
        continue;
    end
    
    for h = 1:nH
        % -- 为表格计算：所有年龄组的加权平均CV --
        total_mass_h_yr = 0;
        weighted_EV_nopps_h_yr = 0;
        weighted_EV_pps_h_yr = 0;

        for a_idx = 1:cS_nopps.aD_new
            dist_nopps_slice = data_nopps.final_Dist_path_h{h}(:, :, :, a_idx, t_idx);
            val_nopps_slice = Val_path_nopps_h{h}(:, :, :, a_idx, t_idx);
            mass_slice = sum(dist_nopps_slice, 'all');

            if mass_slice > 1e-12
                total_mass_h_yr = total_mass_h_yr + mass_slice;
                weighted_EV_nopps_h_yr = weighted_EV_nopps_h_yr + sum(val_nopps_slice .* dist_nopps_slice, 'all');
                
                dist_pps_slice = data_pps.final_Dist_path_h{h}(:, :, :, a_idx, t_idx);
                val_pps_slice = Val_path_pps_h{h}(:, :, :, a_idx, t_idx);
                weighted_EV_pps_h_yr = weighted_EV_pps_h_yr + sum(val_pps_slice .* dist_pps_slice, 'all');
            end
            
            % -- 为绘图计算：特定年龄组的CV --
            for i_age_grp = 1:num_age_groups_fig
                if a_idx == age_groups_for_figure{i_age_grp, 2}
                    EV_nopps_agegrp = sum(val_nopps_slice .* dist_nopps_slice, 'all') / max(1e-9, mass_slice);
                    
                    dist_pps_slice_agegrp = data_pps.final_Dist_path_h{h}(:, :, :, a_idx, t_idx);
                    val_pps_slice_agegrp = Val_path_pps_h{h}(:, :, :, a_idx, t_idx);
                    mass_pps_slice_agegrp = sum(dist_pps_slice_agegrp, 'all');
                    EV_pps_agegrp = sum(val_pps_slice_agegrp .* dist_pps_slice_agegrp, 'all') / max(1e-9, mass_pps_slice_agegrp);
                    
                    if EV_nopps_agegrp < 0 && EV_pps_agegrp < 0
                        cv_figure_data(i_age_grp, i_yr, h) = (EV_pps_agegrp / EV_nopps_agegrp)^(1 / (1 - sigma)) - 1;
                    else
                        cv_figure_data(i_age_grp, i_yr, h) = NaN;
                    end
                end
            end
        end

        EV_nopps_avg = weighted_EV_nopps_h_yr / max(1e-9, total_mass_h_yr);
        EV_pps_avg = weighted_EV_pps_h_yr / max(1e-9, total_mass_h_yr);
        
        if EV_nopps_avg < 0 && EV_pps_avg < 0
            cv_table_data(h, i_yr) = (EV_pps_avg / EV_nopps_avg)^(1 / (1 - sigma)) - 1;
        else
            cv_table_data(h, i_yr) = NaN;
        end
    end
end
fprintf('   ✅ 所有群体福利计算完成。\n');


%% --- 4. 生成 LaTeX 表格 ---
fprintf('\n--- 4. 生成 LaTeX 格式的福利分析表格 ---\n');

% --- 准备输出目录和文件 ---
[output_dir_tab, ~, ~] = fileparts(output_tex_filename);
if ~exist(output_dir_tab, 'dir'), mkdir(output_dir_tab); end
fileID = fopen(output_tex_filename, 'w', 'n', 'UTF-8');
if fileID == -1, error('无法创建或打开LaTex输出文件: %s', output_tex_filename); end

% --- 计算总体平均 ---
total_pop_by_year = sum(cS_nopps.Z_path_raw, 1);
cv_total_avg = zeros(1, num_years);
for i_yr = 1:num_years
    year = assessment_years(i_yr);
    t_idx = round((year - cS_nopps.start_year) / cS_nopps.time_Step) + 1;
    if t_idx > cS_nopps.T_sim, cv_total_avg(i_yr) = NaN; continue; end
    
    pop_weights_h = cS_nopps.type_weights;
    valid_mask = ~isnan(cv_table_data(:, i_yr));
    if any(valid_mask)
        cv_total_avg(i_yr) = sum(cv_table_data(valid_mask, i_yr) .* pop_weights_h(valid_mask)) / sum(pop_weights_h(valid_mask));
    else
        cv_total_avg(i_yr) = NaN;
    end
end


% --- 写入表格 ---
fprintf(fileID, '%% =======================================================\n');
fprintf(fileID, '%%  此文件由 welfare_analysis.m (v3.0) 自动生成\n');
fprintf(fileID, '%% =======================================================\n');
fprintf(fileID, '\\begin{table}[htbp]\n');
fprintf(fileID, '\\centering\n');
fprintf(fileID, '\\caption{引入个人养老金制度对不同类型家庭的福利影响 (截面分析)}\n');
fprintf(fileID, '\\label{tab4.1_welfare_analysis_by_type}\n');
fprintf(fileID, '\\begin{threeparttable}\n');

col_defs = ['l', repmat('c', 1, num_years)];
fprintf(fileID, '\\begin{tabular}{%s}\n', col_defs);
fprintf(fileID, '\\toprule\n');

fprintf(fileID, ' & \\multicolumn{%d}{c}{评估年份} \\\\\n', num_years);
fprintf(fileID, '\\cmidrule(lr){2-%d}\n', num_years + 1);

year_header_cells = arrayfun(@(y) num2str(y), assessment_years, 'UniformOutput', false);
year_header_line = strjoin(year_header_cells, ' & ');
fprintf(fileID, '家庭类型 & %s \\\\\n', year_header_line);
fprintf(fileID, '\\midrule\n');

% --- 循环写入各类型数据 ---
for h = 1:nH
    row_data_str = strjoin(arrayfun(@(cv) format_cv(cv), cv_table_data(h, :), 'UniformOutput', false), ' & ');
    fprintf(fileID, '%s & %s \\\\\n', type_labels{h}, row_data_str);
end

fprintf(fileID, '\\midrule\n');

% --- 写入总体平均 ---
total_avg_str = strjoin(arrayfun(@(cv) format_cv(cv), cv_total_avg, 'UniformOutput', false), ' & ');
fprintf(fileID, '\\textbf{总体平均} & %s \\\\\n', total_avg_str);

% --- 写入表格结尾和注释 ---
fprintf(fileID, '\\bottomrule\n');
fprintf(fileID, '\\end{tabular}\n');
fprintf(fileID, '\\begin{tablenotes}[para,flushleft]\n');
fprintf(fileID, '  \\footnotesize\n');
fprintf(fileID, '  \\item[注] 表中数值为“有PPS情景”相对于“无PPS情景”的补偿变化(CV)。CV值为正，代表该群体福利因改革而提高。例如，+1.00\\%%表示该群体愿意放弃其剩余生命周期消费的1\\%%以换取改革。每格数值为该类型家庭在该年份下，对所有年龄组的福利变化进行人口加权平均得到的结果。\n');
fprintf(fileID, '\\end{tablenotes}\n');
fprintf(fileID, '\\end{threeparttable}\n');
fprintf(fileID, '\\end{table}\n');
fclose(fileID);
fprintf('   ✅ 成功生成 LaTeX 表格文件: %s\n', output_tex_filename);

%% --- 5. 生成 2x2 福利变化图 ---
fprintf('\n--- 5. 生成 2x2 福利变化分解图 ---\n');

% Figure尺寸和灰度风格定义
fig = figure('Name', '各类型家庭内部福利变化分解', 'Position', [100 497 632 498]);
% [核心修改] 设置TileSpacing和Padding为'compact'以最小化间距
tcl = tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

% --- 定义灰度绘图风格 ---
style_young = {'-o', 'Color', [0 0 0], 'LineWidth', 2, 'MarkerSize', 5, 'MarkerFaceColor', [0.2 0.2 0.2]}; % 青年: 黑色实线 + 实心圆
style_middle = {'--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2, 'Marker', 's', 'MarkerSize', 5};             % 中年: 深灰虚线 + 空心方块
style_preretire = {':d', 'Color', [0.6 0.6 0.6], 'LineWidth', 2, 'MarkerSize', 5};                        % 临近退休: 浅灰点线 + 空心菱形
styles = {style_young, style_middle, style_preretire};

subplot_titles = {'(a) 高收入城镇职工', '(b) 中低收入城镇职工', ...
                  '(c) 高收入居民', '(d) 中低收入居民'};

for h = 1:nH
    ax = nexttile;
    hold(ax, 'on');
    
    for i_age_grp = 1:num_age_groups_fig
        plot(ax, assessment_years, cv_figure_data(i_age_grp, :, h) * 100, ...
            styles{i_age_grp}{:}, 'DisplayName', age_groups_for_figure{i_age_grp, 1});
    end
    
    yline(ax, 0, 'k--', 'LineWidth', 1, 'HandleVisibility','off');
    
    title(ax, subplot_titles{h}, 'FontName', 'SimSun', 'FontSize', 14);
    xlabel(ax, '评估年份', 'FontName', 'SimSun');
    if mod(h, 2) == 1 % 只在左侧的图表显示Y轴标签
        ylabel(ax, '补偿变化 (CV %)', 'FontName', 'SimSun');
    end
    grid(ax, 'on');
    box(ax, 'on');
    xlim(ax, [assessment_years(1), assessment_years(end)]);
    
    % 只在第一个子图中显示图例
    if h == 1
        legend(ax, 'show', 'Location', 'best', 'FontName', 'SimSun');
    end
end

% --- 保存图像 ---
[output_dir_fig, ~, ~] = fileparts(output_fig_filename);
if ~exist(output_dir_fig, 'dir'), mkdir(output_dir_fig); end
fprintf('   正在保存图像至: %s\n', output_fig_filename);
try
    exportgraphics(fig, output_fig_filename, 'Resolution', 300);
    fprintf('   ✅ 图像保存成功。\n');
catch ME
    fprintf('   ❌ 图像保存失败！错误信息: %s\n', ME.message);
end

fprintf('\n--- 福利分析脚本执行完毕 ---\n');

%% --- [辅助函数] ---
function str_out = format_cv(cv_value)
    % 辅助函数: 将CV数值格式化为带正负号和百分号的字符串
    if isnan(cv_value)
        str_out = '-';
        return;
    end
    
    if cv_value >= 0
        sign_char = '+';
    else
        sign_char = '';
    end
    str_out = sprintf('%s%.2f\\%%', sign_char, cv_value * 100);
end