% =========================================================================
% == SCRIPT: tab_flexiblelabor.m
% == 版本: [v1.4 - 表格拆分版]
% ==
% == 目的:
% ==   1. 生成两个独立的LaTeX表格:
% ==      - tex/tab/tab4.7_a.tex (劳动供给决策)
% ==      - tex/tab/tab4.7_b.tex (宏观变量对比)
% ==   2. 保持v1.3版本的所有计算逻辑和格式优化。
% =========================================================================

% --- 1. 初始化与设定 ---
clear; close all; clc;
addpath(pwd);
fprintf('=== 弹性/固定劳动供给对比分析脚本 (v1.4) ===\n\n');

% --- 分析设定 ---
scenarios = {'flexible', 'fixed'};

% --- [核心修改] 路径设定: 定义两个输出文件 ---
file_tpi_pattern = 'v1_flexible_labor/TRANS/TPI_results_het_nopps_%s.mat';
output_dir_tex = 'tex/tab/';
if ~exist(output_dir_tex, 'dir'), mkdir(output_dir_tex); end
output_tex_table_a = 'tex/tab/tab4.7_a_flexible_labor.tex';
output_tex_table_b = 'tex/tab/tab4.7_b_flexible_marco.tex';

% --- 表格参数 ---
assessment_years = [2030, 2040, 2050, 2060, 2070, 2100];

type_labels = {'高收入', '中低收入'};
sector_labels = {'A. 城镇职工', 'B. 城乡居民'};

macro_indicators = {'y_hat_pc', 'k_p_hat_pc', 'c_hat_pc', 'r', 'adj'};
macro_indicator_labels = {'有效人均产出 ($\hat{y}$)', '有效人均私人资本 ($\hat{k}_p$)', ...
                           '有效人均消费 ($\hat{c}$)', ...
                          '年化真实利率 (r, \%)', '养老金均衡调整因子 (adj)'};


%% --- 2. [数据加载] 加载所有情景的计算结果 ---
fprintf('--- 2. [数据加载] 循环加载各方案结果 ---\n');
results_all = struct();
for i = 1:length(scenarios)
    scenario_name = scenarios{i};
    fprintf('   正在加载方案: [%s] ...\n', upper(scenario_name));
    file_tpi = sprintf(file_tpi_pattern, scenario_name);
    if ~exist(file_tpi, 'file')
        error('找不到文件: %s，请确认main_run_trans.m已针对该情景运行。', file_tpi); 
    end
    results_all.(scenario_name) = load(file_tpi);
    fprintf('   ✅ 已加载: %s\n', file_tpi);
end
fprintf('--- 所有情景数据加载完毕 ---\n');


%% --- 3. [Panel A 计算] 平均劳动供给 ---
fprintf('\n--- 3. [Panel A 计算] 平均劳动供给决策 (l_choice) ---\n');
data_flex = results_all.flexible;
cS_flex = data_flex.cS;
paramS_flex = data_flex.results.ss_data.paramSF; 

age_groups = {
    '青年 (25-29岁)', round((25 - cS_flex.age1_orig)/cS_flex.time_Step) + 1;
    '中年 (45-49岁)', round((45 - cS_flex.age1_orig)/cS_flex.time_Step) + 1;
    '临退 (60-64岁)', cS_flex.aR_new_path(1); 
};

panel_a_data = nan(length(sector_labels) * length(type_labels) * size(age_groups, 1), length(assessment_years));

leGridV = paramS_flex.leGridV; 
[~, ~, le_grid_full] = ndgrid(cS_flex.kGridV, cS_flex.kppsGridV, leGridV); 

row_idx_a = 1;
for i_sector = 1:length(sector_labels)
    for i_type = 1:length(type_labels)
        h_idx = (i_sector - 1) * 2 + i_type;
        for i_age = 1:size(age_groups, 1)
            a_idx = age_groups{i_age, 2};
            
            age_eff_ia_h = cS_flex.ageEffV_new_h(a_idx, h_idx);
            if age_eff_ia_h < 1e-9, row_idx_a = row_idx_a + 1; continue; end
            
            denominator_grid = age_eff_ia_h .* le_grid_full;
            denominator_grid(abs(denominator_grid) < 1e-9) = 1e-9;

            for i_yr = 1:length(assessment_years)
                year = assessment_years(i_yr);
                t_idx = round((year - cS_flex.start_year) / cS_flex.time_Step) + 1;
                
                if t_idx > cS_flex.T_sim, continue; end
                
                dist_slice = data_flex.final_Dist_path_h{h_idx}(:, :, :, a_idx, t_idx);
                l_effective_slice = data_flex.final_Pol_path_h{h_idx}{t_idx}(a_idx).l;
                
                l_choice_slice = l_effective_slice ./ denominator_grid;
                
                mass_total = sum(dist_slice, 'all');
                if mass_total > 1e-12
                    weighted_l_choice = sum(dist_slice .* l_choice_slice, 'all');
                    avg_l_choice = weighted_l_choice / mass_total;
                    panel_a_data(row_idx_a, i_yr) = avg_l_choice;
                end
            end
            row_idx_a = row_idx_a + 1;
        end
    end
end
fprintf('   ✅ Panel A (平均劳动决策) 数据计算完成。\n');


%% --- 4. [Panel B 计算] 宏观变量(有效人均)对比 ---
fprintf('\n--- 4. [Panel B 计算] 宏观变量(有效人均)对比 ---\n');
panel_b_data = struct();

for i_scen = 1:length(scenarios)
    scen_name = scenarios{i_scen};
    data = results_all.(scen_name);
    results = data.results;
    cS = data.cS;
    
    total_pop_path = sum(cS.Z_path_raw, 1);
    Y_hat_pc_path = (results.Y_path ./ cS.A_path) ./ total_pop_path;
    K_p_hat_pc_path = (results.K_p_path ./ cS.A_path) ./ total_pop_path;
    C_hat_pc_path = (results.C_path ./ cS.A_path) ./ total_pop_path;
    
    aggr_supply = aggregates.get_path_aggregates(data.final_Dist_path_h, data.final_Pol_path_h, cS, results.ss_data.paramSF);
    L_hat_pc_path = aggr_supply.L_path ./ total_pop_path;

    r_path = results.iter_results.r_path_iter;
    adj_path = results.iter_results.adj_path_iter;
    
    paths = containers.Map(macro_indicators, {Y_hat_pc_path, K_p_hat_pc_path, C_hat_pc_path, r_path, adj_path});
    
    for i_ind = 1:length(macro_indicators)
        indicator = macro_indicators{i_ind};
        current_path = paths(indicator);
        for i_yr = 1:length(assessment_years)
            year = assessment_years(i_yr);
            t_idx = round((year - cS.start_year) / cS.time_Step) + 1;
            
            if t_idx > cS.T_sim, value = NaN;
            else
                value = current_path(t_idx);
                if strcmp(indicator, 'r'), value = ((1 + value)^(1/cS.time_Step) - 1) * 100; end
            end
            panel_b_data.(scen_name).(indicator)(i_yr) = value;
        end
    end
end
fprintf('   ✅ Panel B 数据计算完成。\n');


%% --- 5. [生成 LaTeX 表格] ---
fprintf('\n--- 5. [生成 LaTeX 表格] ---\n');

% --- 5a. 生成表格 a: 劳动供给决策 ---
fileID_a = fopen(output_tex_table_a, 'w', 'n', 'UTF-8');
fprintf(fileID_a, '%% =======================================================\n');
fprintf(fileID_a, '%%  此文件由 tab_flexiblelabor.m (v1.4) 自动生成\n');
fprintf(fileID_a, '%% =======================================================\n');
fprintf(fileID_a, '\\begin{table}[htbp]\n');
fprintf(fileID_a, '\\centering\n');
fprintf(fileID_a, '\\caption{弹性劳动供给模型下高收入家庭的平均劳动供给决策}\n');
fprintf(fileID_a, '\\label{tab4.7a:labor_decision}\n');
fprintf(fileID_a, '\\begin{threeparttable}\n');
fprintf(fileID_a, '\\begin{tabular}{lcccccc}\n');
fprintf(fileID_a, '\\toprule\n');
fprintf(fileID_a, '家庭类型与年龄组');
for year = assessment_years, fprintf(fileID_a, ' & %d', year); end
fprintf(fileID_a, ' \\\\\n');
fprintf(fileID_a, '\\midrule\n');

row_idx_a = 1;
for i_sector = 1:length(sector_labels)
    fprintf(fileID_a, '\\multicolumn{7}{l}{\\textbf{%s}} \\\\\n', sector_labels{i_sector});
    i_type = 1; % 只打印高收入组
    % for i_type = 1:length(type_labels)
    fprintf(fileID_a, '    \\hspace{1em}%s', type_labels{i_type});
    fprintf(fileID_a, repmat(' & ', 1, length(assessment_years)));
    fprintf(fileID_a, ' \\\\\n');
    for i_age = 1:size(age_groups, 1)
        fprintf(fileID_a, '        \\hspace{2em}%s', age_groups{i_age, 1});
        for i_yr = 1:length(assessment_years)
            val = panel_a_data(row_idx_a, i_yr);
            if isnan(val), fprintf(fileID_a, ' & -'); else, fprintf(fileID_a, ' & %.3f', val); end
        end
        fprintf(fileID_a, ' \\\\\n');
        row_idx_a = row_idx_a + 1;
    end
    % end
    row_idx_a = row_idx_a + size(age_groups, 1); % 跳过低收入组的数据行索引
    if i_sector < length(sector_labels), fprintf(fileID_a, '\\addlinespace\n'); end
    
end

fprintf(fileID_a, '\\bottomrule\n');
fprintf(fileID_a, '\\end{tabular}\n');
fprintf(fileID_a, '\\begin{tablenotes}[para,flushleft]\n');
fprintf(fileID_a, '  \\footnotesize\n');
fprintf(fileID_a, '  \\item[注] 表中数值为特定工作年龄群体将其时间禀赋用于市场工作的平均比例。为表格简洁，仅展示高收入群体决策；所有情景下，中低收入群体的平均劳动决策均为1.0。所有情景均基于无个人养老金账户(nopps)的设定。\n');
fprintf(fileID_a, '\\end{tablenotes}\n');
fprintf(fileID_a, '\\end{threeparttable}\n');
fprintf(fileID_a, '\\end{table}\n');
fclose(fileID_a);
fprintf('   ✅ 成功生成 LaTeX 表格文件: %s\n', output_tex_table_a);


% --- 5b. 生成表格 b: 宏观变量对比 ---
fileID_b = fopen(output_tex_table_b, 'w', 'n', 'UTF-8');
fprintf(fileID_b, '%% =======================================================\n');
fprintf(fileID_b, '%%  此文件由 tab_flexiblelabor.m (v1.4) 自动生成\n');
fprintf(fileID_b, '%% =======================================================\n');
fprintf(fileID_b, '\\begin{table}[htbp]\n');
fprintf(fileID_b, '\\centering\n');
fprintf(fileID_b, '\\caption{不同劳动供给模型下的关键宏观变量对比 (有效人均量)}\n');
fprintf(fileID_b, '\\label{tab4.7b:macro_comparison}\n');
fprintf(fileID_b, '\\begin{threeparttable}\n');
fprintf(fileID_b, '\\begin{tabular}{lcccccc}\n');
fprintf(fileID_b, '\\toprule\n');
fprintf(fileID_b, '宏观指标');
for year = assessment_years, fprintf(fileID_b, ' & %d', year); end
fprintf(fileID_b, ' \\\\\n');
fprintf(fileID_b, '\\midrule\n');

for i_ind = 1:length(macro_indicators)
    indicator = macro_indicators{i_ind};
    fprintf(fileID_b, '%s', macro_indicator_labels{i_ind});
    fprintf(fileID_b, repmat(' & ', 1, length(assessment_years)));
    fprintf(fileID_b, ' \\\\\n');
    
    for i_scen = 1:length(scenarios)
        scen_name = scenarios{i_scen};
        if strcmp(scen_name, 'flexible'), row_label = '弹性供给'; else, row_label = '固定供给'; end
        fprintf(fileID_b, '    \\hspace{1em}%s', row_label);
        
        for i_yr = 1:length(assessment_years)
            val = panel_b_data.(scen_name).(indicator)(i_yr);
             if isnan(val), fprintf(fileID_b, ' & -'); else, fprintf(fileID_b, ' & %.3f', val); end
        end
        fprintf(fileID_b, ' \\\\\n');
    end
    if i_ind < length(macro_indicators), fprintf(fileID_b, '\\addlinespace\n'); end
end

fprintf(fileID_b, '\\bottomrule\n');
fprintf(fileID_b, '\\end{tabular}\n');
fprintf(fileID_b, '\\begin{tablenotes}[para,flushleft]\n');
fprintf(fileID_b, '  \\footnotesize\n');
fprintf(fileID_b, '  \\item[注] 表中对比了两种模型下的关键宏观变量路径，所有总量均为剔除技术(A)和人口(N)双重影响后的有效人均量。固定供给模型假设所有工作年龄个体均提供1单位的劳动。所有情景均基于无个人养老金账户(nopps)的设定。\n');
fprintf(fileID_b, '\\end{tablenotes}\n');
fprintf(fileID_b, '\\end{threeparttable}\n');
fprintf(fileID_b, '\\end{table}\n');
fclose(fileID_b);
fprintf('   ✅ 成功生成 LaTeX 表格文件: %s\n', output_tex_table_b);

fprintf('\n--- 脚本执行完毕 ---\n');