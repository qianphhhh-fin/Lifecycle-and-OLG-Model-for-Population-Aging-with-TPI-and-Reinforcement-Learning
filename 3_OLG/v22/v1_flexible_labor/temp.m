% =========================================================================
% == SCRIPT: analyze_tfp_scenarios.m
% == 版本: [v1.0 - TFP情景分析与adj框架]
% ==
% == 目的:
% ==   1. 生成 表 4.3: 不同经济增长情景下的宏观经济与养老金可持续性
% ==      关键指标，核心为内生调整因子(adj)。
% ==   2. 生成 表 4.4: 不同经济增长情景对各类型家庭的福利影响
% ==      (CV, %)，以基准情景为比较对象。
% ==
% == 前置条件:
% ==   - 已成功运行 main_run_trans.m 并生成以下文件:
% ==     - TRANS/TPI_results_het_nopps_TFP_baseline.mat
% ==     - TRANS/TPI_results_het_nopps_TFP_optimistic.mat
% ==     - TRANS/TPI_results_het_nopps_TFP_pessimistic.mat
% =========================================================================

clear; close all; clc;
addpath(pwd);
fprintf('=== TFP冲击情景对比分析脚本 (adj框架 v1.0) ===\n\n');

%% --- 1. 用户设定与文件定义 ---
fprintf('--- 1. 用户设定与文件定义 ---\n');

% --- 分析设定 ---
scenarios = {'pessimistic', 'baseline', 'optimistic'};
base_scenario = 'baseline'; % 福利分析的基准情景

% --- 文件路径模板 ---
file_tpi_pattern = 'TRANS/TPI_results_het_nopps_TFP_%s.mat';

% --- 输出路径 ---
output_dir_tex = 'tex/tab/';
if ~exist(output_dir_tex, 'dir'), mkdir(output_dir_tex); end
output_tex_table_macro = fullfile(output_dir_tex, 'tab4.3_macro_sustainability.tex');
output_tex_table_welfare = fullfile(output_dir_tex, 'tab4.4_welfare_analysis.tex');

fprintf('   福利分析基准方案: %s\n', base_scenario);
fprintf('   LaTeX 表格将输出至: %s\n', output_dir_tex);


%% --- 2. [数据加载] 加载所有情景的计算结果 ---
fprintf('\n--- 2. [数据加载] 循环加载各方案结果 ---\n');

results_all = struct(); % 用于存储所有方案的核心计算结果

for i = 1:length(scenarios)
    scenario_name = scenarios{i};
    fprintf('   正在加载方案: [%s] ...\n', upper(scenario_name));

    % --- 加载数据 ---
    file_tpi = sprintf(file_tpi_pattern, scenario_name);
    if ~exist(file_tpi, 'file')
        error('找不到文件: %s，请确认main_run_trans.m已针对该情景运行。', file_tpi); 
    end

    results_all.(scenario_name) = load(file_tpi);
    fprintf('   ✅ 已加载: %s\n', file_tpi);
end
fprintf('--- 所有情景数据加载完毕 ---\n');




%% --- 3. [生成表4.3] 宏观经济与养老金可持续性 ---
fprintf('\n--- 3. [生成表4.3] 宏观经济与养老金可持续性 ---\n');

% --- [核心修正] 定义评估年份和指标 (增加 K/AL) ---
assessment_years_macro = [2030, 2050];
macro_indicators = {'r', 'w_hat', 'K_AL', 'adj'};
indicator_labels = {'年化真实利率 (r, \%)', '有效工资率 (w, 指数)', '资本-效率劳动比 (K/AL)', '均衡调整因子 (adj)'};
section_labels = {'A. 宏观经济与要素价格', 'B. PAYG养老金体系可持续性'};

% --- 打开文件并写入LaTeX表头 ---
fileID = fopen(output_tex_table_macro, 'w', 'n', 'UTF-8');
fprintf(fileID, '%% =======================================================\n');
fprintf(fileID, '%%  此文件由 analyze_tfp_scenarios.m (v1.4, K/AL版) 自动生成\n');
fprintf(fileID, '%% =======================================================\n');
fprintf(fileID, '\\begin{table}[htbp]\n');
fprintf(fileID, '\\centering\n');
fprintf(fileID, '\\caption{不同经济增长情景下的宏观经济与养老金可持续性关键指标}\n');
fprintf(fileID, '\\label{tab4.3:macro_sustainability_adj}\n');
fprintf(fileID, '\\begin{threeparttable}\n');
fprintf(fileID, '\\begin{tabular}{llccc}\n');
fprintf(fileID, '\\toprule\n');
fprintf(fileID, '指标 & 评估年份 & 悲观情景 & 基准情景 & 乐观情景 \\\\\n');
fprintf(fileID, '\\midrule\n');

% --- 循环提取数据并写入表格 ---
for i_sec = 1:length(section_labels)
    fprintf(fileID, '\\multicolumn{5}{l}{\\textbf{%s}} \\\\\n', section_labels{i_sec});
    
    current_indicators = {};
    if i_sec == 1, current_indicators = {'r', 'w_hat', 'K_AL'}; end
    if i_sec == 2, current_indicators = {'adj'}; end

    for i_ind = 1:length(current_indicators)
        indicator = current_indicators{i_ind};
        label_idx = find(strcmp(macro_indicators, indicator));
        
        fprintf(fileID, '%s', indicator_labels{label_idx});

        % --- [核心修正] 写入期初均衡数据 (从各自的 ss0 提取) ---
        fprintf(fileID, ' & 期初均衡');
        for i_scen = 1:length(scenarios)
            scen_name = scenarios{i_scen};
            data = results_all.(scen_name);
            ss0 = data.ss0;
            cS = data.cS;
            
            switch indicator
                case 'r'
                    value = ((1 + ss0.r_mkt)^(1/cS.time_Step) - 1) * 100;
                    fprintf(fileID, ' & %.3f\\%%', value);
                case 'w_hat'
                    value = ss0.w_hat;
                    fprintf(fileID, ' & %.3f', value);
                case 'K_AL'
                    % K_AL_ss = K_p_hat_ss / (A_ss * L_hat_ss) -> A_ss=1, K_p_hat / L_hat
                    value = ss0.K_private_hat / ss0.L_hat;
                    fprintf(fileID, ' & %.3f', value);
                case 'adj'
                    value = ss0.adj;
                    fprintf(fileID, ' & %.3f', value);
            end
        end
        fprintf(fileID, ' \\\\\n');
        
        % --- 写入转轨期数据 ---
        for i_yr = 1:length(assessment_years_macro)
            year = assessment_years_macro(i_yr);
            fprintf(fileID, ' & %d', year); % 评估年份列

            for i_scen = 1:length(scenarios)
                scen_name = scenarios{i_scen};
                data = results_all.(scen_name);
                t_idx = round((year - data.cS.start_year) / data.cS.time_Step) + 1;
                
                if t_idx > 0 && t_idx <= data.cS.T_sim
                    switch indicator
                        case 'r'
                            r_period = data.results.iter_results.r_path_iter(t_idx);
                            value = ((1 + r_period)^(1/data.cS.time_Step) - 1) * 100;
                            fprintf(fileID, ' & %.3f\\%%', value);
                        case 'w_hat'
                            value = data.results.iter_results.w_hat_path_iter(t_idx);
                            fprintf(fileID, ' & %.3f', value);
                        case 'K_AL'
                            aggr_supply = aggregates.get_path_aggregates(data.final_Dist_path_h, data.final_Pol_path_h, data.cS, data.results.ss_data.paramSF);
                            K_p_hat_t = aggr_supply.K_p_hat_total_path(t_idx);
                            L_hat_t = aggr_supply.L_path(t_idx);
                            % K/AL = K_hat*A*N / (A*L_hat*N) = K_hat / L_hat
                            value = K_p_hat_t / L_hat_t;
                            fprintf(fileID, ' & %.3f', value);
                        case 'adj'
                            value = data.results.iter_results.adj_path_iter(t_idx);
                            fprintf(fileID, ' & %.3f', value);
                    end
                else
                    fprintf(fileID, ' & -');
                end
            end
            fprintf(fileID, ' \\\\\n');
        end
        
        % --- 写入终期稳态数据 ---
        fprintf(fileID, ' & 终期稳态');
        for i_scen = 1:length(scenarios)
            scen_name = scenarios{i_scen};
            data = results_all.(scen_name);
            ssF = data.ssF;
            cS = data.cS;

            switch indicator
                case 'r'
                    value = ((1 + ssF.r_mkt)^(1/cS.time_Step) - 1) * 100;
                    fprintf(fileID, ' & %.3f\\%%', value);
                case 'w_hat'
                    value = ssF.w_hat;
                    fprintf(fileID, ' & %.3f', value);
                case 'K_AL'
                    value = ssF.K_private_hat / ssF.L_hat;
                    fprintf(fileID, ' & %.3f', value);
                case 'adj'
                    value = ssF.adj;
                    fprintf(fileID, ' & %.3f', value);
            end
        end
        fprintf(fileID, ' \\\\\n');
        
        if (i_sec == 1 && i_ind < length(current_indicators)) || (i_sec == 2 && i_ind < length(current_indicators))
            fprintf(fileID, '\\addlinespace\n');
        end
        
    end
end

% --- 写入表尾和注释 ---
fprintf(fileID, '\\bottomrule\n');
fprintf(fileID, '\\end{tabular}\n');
fprintf(fileID, '\\begin{tablenotes}\n');
fprintf(fileID, '  \\footnotesize\n');
fprintf(fileID, '  \\item[注] 资本-效率劳动比(K/AL)是私人资本存量与经技术调整后的总劳动供给之比。均衡调整因子(adj)是在保持缴费率固定的前提下，为保证PAYG体系收支平衡所必需的福利调整系数。adj=1.0表示福利可100\\%%兑现。有效工资率(w)是剔除技术进步影响的工资指数。所有情景均基于无个人养老金账户(nopps)的设定。\n');
fprintf(fileID, '\\end{tablenotes}\n');
fprintf(fileID, '\\end{threeparttable}\n');
fprintf(fileID, '\\end{table}\n');
fclose(fileID);
fprintf('   ✅ 成功生成 LaTeX 表格文件: %s\n', output_tex_table_macro);

%% --- 4. [福利分析] 为各情景重构价值函数路径 ---
fprintf('\n--- 4. [福利分析] 重构价值函数路径 ---\n');

for i_scen = 1:length(scenarios)
    scen_name = scenarios{i_scen};
    fprintf('   正在为方案 [%s] 运行反向迭代...\n', upper(scen_name));
    
    data = results_all.(scen_name);
    cS = data.cS;
    paramSF = data.results.ss_data.paramSF;
    valF_h = data.results.ss_data.valF_h;
    polF_h = data.results.ss_data.polF_h;
    iter_results = data.results.iter_results;
    
    Val_path_h = cell(cS.nTypes, 1);
    
    for h = 1:cS.nTypes
        % 获取向后迭代所需的所有路径
        pathS_h = get_pathS_for_backward(iter_results, h, cS, data.ss0);
        
        % 获取该类型特定的参数
        cS_h = get_cs_for_type(cS, h);
        
        % 运行后向迭代
        [~, Val_path_h{h}] = household.backward_hh(pathS_h, cS_h, paramSF, valF_h{h}, polF_h{h});
    end
    
    % 将计算得到的价值函数存回主结构体
    results_all.(scen_name).Val_path_h = Val_path_h;
    fprintf('   ✅ 方案 [%s] 价值函数重构完成。\n', upper(scen_name));
end
fprintf('--- 所有情景价值函数路径重构完毕 ---\n');


%% --- 5. [生成表4.4] 计算并生成福利分析(CV)表格 ---
fprintf('\n--- 5. [生成表4.4] 计算并生成福利分析(CV)表格 ---\n');

% --- 定义分析维度 ---
assessment_years_welfare = [2030, 2050, 2070, 2100];
% [核心修正] 明确定义改革情景的顺序，以匹配LaTex表头
reform_scenarios = {'optimistic', 'pessimistic'}; 
type_labels = {'高收入', '中低收入'};
sector_labels = {'A. 城镇职工', 'B. 城乡居民'};

% --- 获取基准情景的数据 ---
data_base = results_all.(base_scenario);
cS_base = data_base.cS;
sigma = cS_base.sigma;
T_sim_base = cS_base.T_sim;

% --- [核心修正] 在 cS_base 定义后，直接计算年龄组索引，避免使用 eval ---
age_groups = {
    '青年 (25-29岁)', round((25 - cS_base.age1_orig)/cS_base.time_Step) + 1;
    '中年 (45-49岁)', round((45 - cS_base.age1_orig)/cS_base.time_Step) + 1;
    '临退 (60-64岁)', cS_base.aR_new_path(1) % 使用初始退休年龄作为代表
    '老年 (70-74岁)', round((70 - cS_base.age1_orig)/cS_base.time_Step) + 1;
};

% --- 打开文件并写入LaTeX表头 ---
fileID_welfare = fopen(output_tex_table_welfare, 'w', 'n', 'UTF-8');
fprintf(fileID_welfare, '%% =======================================================\n');
fprintf(fileID_welfare, '%%  此文件由 analyze_tfp_scenarios.m (v1.2, 修正版) 自动生成\n');
fprintf(fileID_welfare, '%% =======================================================\n');
fprintf(fileID_welfare, '\\begin{table}[htbp]\n');
fprintf(fileID_welfare, '\\centering\n');
fprintf(fileID_welfare, '\\caption{不同经济增长情景对各类型家庭的福利影响 (相对于基准情景的补偿变化 CV, \\%%)}\n');
fprintf(fileID_welfare, '\\label{tab4.4:welfare_impact_tfp}\n');
fprintf(fileID_welfare, '\\begin{threeparttable}\n');

% [核心修正] 动态、安全地构建并写入表头
header_cols = 'l';
header_line_1 = {'家庭类型与年龄组'};
header_line_2 = {''};
cmidrule_parts = {};
num_data_cols = 0;

for i_yr = 1:length(assessment_years_welfare)
    header_cols = [header_cols, 'rr'];
    header_line_1{end+1} = sprintf('\\multicolumn{2}{c}{%d}', assessment_years_welfare(i_yr));
    header_line_2{end+1} = '乐观';
    header_line_2{end+1} = '悲观';
    cmidrule_parts{end+1} = sprintf('\\cmidrule(lr){%d-%d}', num_data_cols + 2, num_data_cols + 3);
    num_data_cols = num_data_cols + 2;
end

final_header_1 = strjoin(header_line_1, ' & ');
final_header_2 = strjoin(header_line_2, ' & ');
final_cmidrule = strjoin(cmidrule_parts, ' ');

fprintf(fileID_welfare, '\\begin{tabular}{%s}\n', header_cols);
fprintf(fileID_welfare, '\\toprule\n');
fprintf(fileID_welfare, '%s \\\\\n', final_header_1);
fprintf(fileID_welfare, '%s\n', final_cmidrule);
fprintf(fileID_welfare, '%s \\\\\n', final_header_2);
fprintf(fileID_welfare, '\\midrule\n');

% --- [核心逻辑修正] 核心CV计算与写入 ---
for i_sector = 1:length(sector_labels)
    fprintf(fileID_welfare, '\\multicolumn{%d}{l}{\\textbf{%s}} \\\\\n', num_data_cols + 1, sector_labels{i_sector});
    for i_type = 1:length(type_labels)
        h_idx = (i_sector - 1) * 2 + i_type;
        
        % 打印类型标签，并填充该行的空数据列
        fprintf(fileID_welfare, '    \\hspace{1em}%s', type_labels{i_type});
        fprintf(fileID_welfare, repmat(' & ', 1, num_data_cols));
        fprintf(fileID_welfare, ' \\\\\n');

        for i_age = 1:size(age_groups, 1)
            age_label = age_groups{i_age, 1};
            a_idx = age_groups{i_age, 2};
            
            fprintf(fileID_welfare, '       \\hspace{2em}%s', age_label);

            for i_yr = 1:length(assessment_years_welfare)
                year = assessment_years_welfare(i_yr);
                t_idx = round((year - cS_base.start_year) / cS_base.time_Step) + 1;

                if t_idx > T_sim_base
                    fprintf(fileID_welfare, ' & - & -');
                    continue;
                end
                
                % 按照 "乐观 vs. 基准", "悲观 vs. 基准" 的顺序循环
                for i_ref = 1:length(reform_scenarios)
                    reform_name = reform_scenarios{i_ref};
                    data_reform = results_all.(reform_name);
                    
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
                    fprintf(fileID_welfare, ' & %s', format_cv(cv));
                end % reform scenarios loop
            end % year loop
            fprintf(fileID_welfare, ' \\\\\n');
        end % age group loop
    end % type loop
    if i_sector == 1
         fprintf(fileID_welfare, '\\midrule\n');
    end
end % sector loop

% --- 计算并写入总体平均 ---
fprintf(fileID_welfare, '\\midrule\n');
fprintf(fileID_welfare, '\\multicolumn{%d}{l}{\\textbf{C. 总体平均}} \\\\\n', num_data_cols + 1);
fprintf(fileID_welfare, '    人口加权平均 CV');

for i_yr = 1:length(assessment_years_welfare)
    year = assessment_years_welfare(i_yr);
    t_idx = round((year - cS_base.start_year) / cS_base.time_Step) + 1;
    
    if t_idx > T_sim_base
        fprintf(fileID_welfare, ' & - & -');
        continue;
    end
    
    total_pop_t = sum(cS_base.Z_path_raw(:, t_idx), 'all');
    
    for i_ref = 1:length(reform_scenarios)
        reform_name = reform_scenarios{i_ref};
        data_reform = results_all.(reform_name);
        
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
        fprintf(fileID_welfare, ' & %s', format_cv(cv_total));
    end
end
fprintf(fileID_welfare, ' \\\\\n');


% --- 写入表尾 ---
fprintf(fileID_welfare, '\\bottomrule\n');
fprintf(fileID_welfare, '\\end{tabular}\n');
fprintf(fileID_welfare, '\\begin{tablenotes}[para,flushleft]\n');
fprintf(fileID_welfare, '  \\footnotesize\n');
% [核心修正] 更新注释以匹配简化的表头
fprintf(fileID_welfare, '  \\item[注] 表中数值为补偿变化(CV)，衡量“乐观”或“悲观”情景相对于“基准”情景的福利差异。正值代表福利改善。例如，+0.8\\%%表示该群体愿意放弃其在基准情景下剩余生命周期消费的0.8\\%%，以换取进入乐观增长情景。所有计算均基于固定缴费率和内生福利调整因子(adj)的框架。\n');
fprintf(fileID_welfare, '\\end{tablenotes}\n');
fprintf(fileID_welfare, '\\end{threeparttable}\n');
fprintf(fileID_welfare, '\\end{table}\n');
fclose(fileID_welfare);
fprintf('   ✅ 成功生成 LaTeX 表格文件: %s\n', output_tex_table_welfare);

fprintf('\n--- 脚本执行完毕 ---\n');

%% --- 7. [生成图4.7] 绘制三种情景的TFP增长率路径 ---
fprintf('\n--- 7. [生成图4.7] 绘制三种情景的TFP增长率路径 ---\n');

% --- 定义输出路径 ---
output_dir_fig = 'tex/fig/';
if ~exist(output_dir_fig, 'dir'), mkdir(output_dir_fig); end
output_fig_tfp = fullfile(output_dir_fig, 'fig4.7_tfp_path.png');

% --- 设定绘图参数 ---
% 从任一情景加载cS以获取时间轴参数
cS = results_all.baseline.cS; 
annual_years_vec = cS.start_year:cS.end_year;

% --- 设定灰度颜色和线型 ---
color_pessimistic = [0.65 0.65 0.65]; 
color_baseline = [0.0 0.0 0.0]; 
color_optimistic = [0.35 0.35 0.35];
colors = containers.Map({'pessimistic', 'baseline', 'optimistic'}, {color_pessimistic, color_baseline, color_optimistic});
line_styles = {':', '-', '--'}; 
line_styles_map = containers.Map({'pessimistic', 'baseline', 'optimistic'}, line_styles);
% [核心修正] 中文图例标签
legend_labels_map = containers.Map({'pessimistic', 'baseline', 'optimistic'}, {'悲观情景', '基准情景', '乐观情景'});

% --- [核心修正] 初始化图形，使用指定位置 ---
fig_tfp = figure('Name', 'TFP增长率路径情景', 'Position', [100 428 425 172]);
hold on;

% --- 循环生成并绘制路径 ---
for i = 1:length(scenarios)
    scen_name = scenarios{i};
    % 调用函数生成该情景的年度增长率路径
    [~, g_path_annual, ~] = utils.generate_tfp_path(cS, scen_name, false);
    
    plot(annual_years_vec, g_path_annual * 100, ...
        'Color', colors(scen_name), ...
        'LineStyle', line_styles_map(scen_name), ...
        'LineWidth', 2.5, ...
        'DisplayName', legend_labels_map(scen_name));
end

% --- [核心修正] 添加图表格式，移除 yline ---
xlabel('年份', 'FontName', 'SimSun');
ylabel('年化增长率 (%)', 'FontName', 'SimSun');
legend('show', 'Location', 'best', 'FontName', 'SimSun');
grid on;
box on;
xlim([cS.start_year, 2150]);

% --- 保存图像 ---
fprintf('   正在保存TFP路径图至: %s\n', output_fig_tfp);
try
    exportgraphics(fig_tfp, output_fig_tfp, 'Resolution', 300);
    fprintf('   ✅ 图像保存成功。\n');
catch ME
    fprintf('   ❌ 图像保存失败！错误信息: %s\n', ME.message);
end

%% --- [本地辅助函数] ---

function pathS_h = get_pathS_for_backward(iter_results, h, cS, ss0)
    % =====================================================================
    % == 函数: get_pathS_for_backward (本地辅助函数)
    % == 目的: 为指定的家庭类型 h，打包后向迭代所需的所有时间路径。
    % =====================================================================
    pathS_h = struct();
    pathS_h.r_path = iter_results.r_path_iter;
    w_hat_path_converged = iter_results.w_hat_path_iter;
    pathS_h.w_hat_path = w_hat_path_converged;
    pathS_h.tr_per_hh_hat_path = iter_results.TR_hat_pc_path_iter;
    
    % --- 重建养老金路径 (b_hat_path) ---
    % 养老金公式基于上一期的工资，所以需要构造一个错位的工资路径
    w_hat_path_for_b = [ss0.w_hat, w_hat_path_converged(1:end-1)];
    b_hat_formula_path_h_all = zeros(cS.nTypes, cS.T_sim);
    for t = 1:cS.T_sim
        % 调用SS命名空间下的静态函数来计算公式养老金
        b_hat_formula_path_h_all(:, t) = SS.calculate_formula_benefits(w_hat_path_for_b(t), cS);
    end
    adj_path_converged = iter_results.adj_path_iter;
    pathS_h.b_hat_path = adj_path_converged .* b_hat_formula_path_h_all(h, :);
    
    % --- 重建增长率和缴费率路径 ---
    A_path_ext = [cS.A_path, cS.A_path(end) * (1 + ((1+cS.g_A_ss)^cS.time_Step-1))];
    g_A_path = A_path_ext(2:end)./A_path_ext(1:end-1)-1;
    pathS_h.g_A_path = g_A_path(1:cS.T_sim);
    pathS_h.theta_path = cS.theta_path_h(h, :);
end


function cS_h = get_cs_for_type(cS, h)
    % =====================================================================
    % == 函数: get_cs_for_type (本地辅助函数)
    % == 目的: 为指定的家庭类型 h，准备一个专属的 cS 参数包。
    % =====================================================================
    cS_h = cS;
    cS_h.ageEffV_new = cS.ageEffV_new_h(:, h);
    % aR_new_path 已经是特定于此cS的，无需修改
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
    str_out = sprintf('%s%.1f\\%%', sign_char, cv_value * 100);
end