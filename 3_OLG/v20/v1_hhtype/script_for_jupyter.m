% =========================================================================
% == SCRIPT: welfare_analysis.m
% == 版本: [v2.1 - 价值函数重构简化版]
% ==
% == 目的:
% ==   - [核心修正] 简化了价值函数重构过程，移除了对 `final_aggr_supply`
% ==     和 `calculate_endogenous_TR_path` 的依赖，解决了 "无法识别的
% ==     字段名称" 错误。
% ==   - 现在直接使用 TPI 迭代收敛后的宏观价格和福利路径作为反向求解的输入。
% ==   - 分析逻辑和最终输出的表格结构保持不变。
% ==
% == 前置条件:
% ==   - 已运行 `fix_tpi_results_variables.m` 修正了数据文件。
% ==   - TRANS/TPI_results_main_nopps.mat
% ==   - TRANS/TPI_results_main_pps.mat
% ==   - SS/data_for_transition_nopps.mat
% ==   - SS/data_for_transition_pps.mat
% =========================================================================

clear; close all; clc;
addpath(pwd);
fprintf('=== OLG模型年龄组截面福利分析脚本 (v2.1) ===\n\n');

%% --- 1. 设置与加载数据 ---
fprintf('--- 1. 设置与加载数据 ---\n');

% --- 文件路径定义 ---
file_path_nopps = 'TRANS/TPI_results_main_nopps.mat';
file_path_pps = 'TRANS/TPI_results_main_pps.mat';
iter_file_nopps = 'TRANS/iter_results_main_nopps.mat';
iter_file_pps = 'TRANS/iter_results_main_pps.mat';
ss_data_nopps_path = 'SS/data_for_transition_nopps.mat';
ss_data_pps_path = 'SS/data_for_transition_pps.mat';
output_tex_filename = 'tex/tab/tab4.1_welfare_analysis_by_age_group.tex';

% --- 加载数据 ---
fprintf('   正在加载 "无PPS" (基准) 和 "有PPS" (改革) 情景数据...\n');
data_nopps = load(file_path_nopps);
iter_data_nopps = load(iter_file_nopps);
ss_data_nopps = load(ss_data_nopps_path);

data_pps = load(file_path_pps);
iter_data_pps = load(iter_file_pps);
ss_data_pps = load(ss_data_pps_path);
fprintf('   ✅ 数据加载完成。\n');


%% --- 2. [核心修正] 简化价值函数路径重构 ---
fprintf('\n--- 2. 重构价值函数路径 (简化版) ---\n');

% --- a. 重构 "无PPS" 情景的价值函数路径 ---
fprintf('   正在为 "无PPS" 情景运行反向迭代...\n');
cS_nopps = data_nopps.cS;
paramSF_nopps = ss_data_nopps.data_for_transition.paramSF;
pathS_nopps = struct();
pathS_nopps.r_path = data_nopps.results.r_path;
pathS_nopps.w_hat_path = data_nopps.results.w_path ./ cS_nopps.A_path;
% 直接使用 TPI 迭代收敛的路径作为输入
pathS_nopps.b_hat_path = iter_data_nopps.b_hat_path_iter;
pathS_nopps.tr_per_hh_hat_path = iter_data_nopps.TR_hat_pc_path_iter;
% 填充外生路径
A_path_ext_nopps = [cS_nopps.A_path, cS_nopps.A_path(end) * (1 + ((1+cS_nopps.g_A_ss)^cS_nopps.time_Step-1))];
pathS_nopps.g_A_path = A_path_ext_nopps(2:end)./A_path_ext_nopps(1:end-1)-1;
pathS_nopps.theta_path = cS_nopps.theta_path;
% 运行反向迭代
[~, Val_path_nopps] = household.backward_hh(pathS_nopps, cS_nopps, paramSF_nopps, ...
    ss_data_nopps.data_for_transition.valF, ss_data_nopps.data_for_transition.polF);

% --- b. 重构 "有PPS" 情景的价值函数路径 ---
fprintf('   正在为 "有PPS" 情景运行反向迭代...\n');
cS_pps = data_pps.cS;
paramSF_pps = ss_data_pps.data_for_transition.paramSF;
pathS_pps = struct();
pathS_pps.r_path = data_pps.results.r_path;
pathS_pps.w_hat_path = data_pps.results.w_path ./ cS_pps.A_path;
% 直接使用 TPI 迭代收敛的路径作为输入
pathS_pps.b_hat_path = iter_data_pps.b_hat_path_iter;
pathS_pps.tr_per_hh_hat_path = iter_data_pps.TR_hat_pc_path_iter;
% 填充外生路径
A_path_ext_pps = [cS_pps.A_path, cS_pps.A_path(end) * (1 + ((1+cS_pps.g_A_ss)^cS_pps.time_Step-1))];
pathS_pps.g_A_path = A_path_ext_pps(2:end)./A_path_ext_pps(1:end-1)-1;
pathS_pps.theta_path = cS_pps.theta_path;
% 运行反向迭代
[~, Val_path_pps] = household.backward_hh(pathS_pps, cS_pps, paramSF_pps, ...
    ss_data_pps.data_for_transition.valF, ss_data_pps.data_for_transition.polF);

fprintf('   ✅ 价值函数路径重构完成。\n');


%% --- 3. [核心] 计算特定截面下各年龄组的福利变化 (CV) ---
fprintf('\n--- 3. 计算各年龄组在特定年份的福利变化 ---\n');

% --- 定义分析维度 ---
assessment_years = [2023, 2030,2040, 2050,2070,2100,2150]; % 选择要进行截面分析的年份

% --- 定义代表性年龄组 (基于模型年龄索引 a_idx) ---
% 假设 cS.age1_orig = 20, cS.time_Step = 5
% a=1 -> 20-24岁, a=2 -> 25-29岁, etc.
aR_idx = cS_nopps.aR_new; % 退休年龄索引
working_age_groups = {
    '青年 (25-29岁)', round((25 - cS_nopps.age1_orig)/cS_nopps.time_Step) + 1;
    '中年 (45-49岁)', round((45 - cS_nopps.age1_orig)/cS_nopps.time_Step) + 1;
    '临近退休 (60-64岁)', aR_idx
};
retired_age_groups = {
    '刚退休 (65-69岁)', aR_idx + 1;
    '中年退休 (75-79岁)', aR_idx + 3;
    '高龄退休 (85-89岁)', aR_idx + 5
};
age_groups_to_analyze = [working_age_groups; retired_age_groups];
num_age_groups = size(age_groups_to_analyze, 1);

% --- 初始化结果存储 ---
% 使用 table 来存储，更易于阅读和导出
varTypes = ['string', repmat("string", 1, length(assessment_years))];
varNames = ['age_group', arrayfun(@(y) sprintf('年份_%d', y), assessment_years, 'UniformOutput', false)];
results_table = table('Size', [num_age_groups, length(varNames)], 'VariableTypes', varTypes, 'VariableNames', varNames);
results_table.age_group = age_groups_to_analyze(:,1);

% --- 核心计算循环 ---
sigma = cS_nopps.sigma;

for col_idx = 1:length(assessment_years)
    year = assessment_years(col_idx);
    t_idx = round((year - cS_nopps.start_year) / cS_nopps.time_Step) + 1;
    
    fprintf('   正在评估年份: %d (t=%d)...\n', year, t_idx);
    
    if t_idx > cS_nopps.T_sim
        fprintf('      警告: 评估年份超出模拟范围，此列将留空。\n');
        continue;
    end
    
    for row_idx = 1:num_age_groups
        age_label = age_groups_to_analyze{row_idx, 1};
        a_idx = age_groups_to_analyze{row_idx, 2};
        
        if a_idx > cS_nopps.aD_new
            results_table(row_idx, col_idx+1) = {"-"};
            continue;
        end
        
        % --- 计算期望剩余效用 (无PPS) ---
        dist_nopps = data_nopps.final_Dist_path(:, :, :, a_idx, t_idx);
        val_nopps = Val_path_nopps(:, :, :, a_idx, t_idx);
        mass_nopps = sum(dist_nopps, 'all');
        EV_nopps = sum(val_nopps .* dist_nopps, 'all') / max(1e-9, mass_nopps);
        
        % --- 计算期望剩余效用 (有PPS) ---
        dist_pps = data_pps.final_Dist_path(:, :, :, a_idx, t_idx);
        val_pps = Val_path_pps(:, :, :, a_idx, t_idx);
        mass_pps = sum(dist_pps, 'all');
        EV_pps = sum(val_pps .* dist_pps, 'all') / max(1e-9, mass_pps);
        
        % --- 计算 CV ---
        if EV_nopps >= 0 || EV_pps >= 0
            cv_str = 'N/A';
        else
            cv = (EV_pps / EV_nopps)^(1 / (1 - sigma)) - 1;
            cv_str = sprintf('%.2f\\%%', cv * 100);
        end
        results_table(row_idx, col_idx+1) = {cv_str};
    end
end
disp('福利分析结果:');
disp(results_table);
fprintf('   ✅ 所有年龄组福利计算完成。\n');

%% --- 4. 生成 LaTeX 表格 ---
fprintf('\n--- 4. 生成 LaTeX 格式的福利分析表格 ---\n');

% --- 准备输出目录和文件 ---
[output_dir, ~, ~] = fileparts(output_tex_filename);
if ~exist(output_dir, 'dir'), mkdir(output_dir); end
fileID = fopen(output_tex_filename, 'w', 'n', 'UTF-8');
if fileID == -1, error('无法创建或打开LaTex输出文件。'); end

% --- 写入表格头部 ---
fprintf(fileID, '%% =======================================================\n');
fprintf(fileID, '%%  此文件由 welfare_analysis.m (v2.3) 自动生成\n');
fprintf(fileID, '%% =======================================================\n');
fprintf(fileID, '\\begin{table}[htbp]\n');
fprintf(fileID, '\\centering\n');
fprintf(fileID, '\\caption{引入个人养老金制度对不同年龄组福利的影响 (截面分析)}\n');
fprintf(fileID, '\\label{tab4.1_welfare_analysis_by_age_group}\n');
fprintf(fileID, '\\begin{threeparttable}\n');

% --- 动态生成 tabular 列定义 ---
num_year_cols = length(assessment_years);
col_defs = ['l', repmat('c', 1, num_year_cols)];
fprintf(fileID, '\\begin{tabular}{%s}\n', col_defs);
fprintf(fileID, '\\toprule\n');

% --- [!!! 核心修正: 创建两级表头 !!!] ---
% --- 第一行表头: 使用 multicolumn ---
%  - 第一个单元格留空
%  - 第二个单元格横跨所有年份列
fprintf(fileID, ' & \\multicolumn{%d}{c}{评估年份} \\\\\n', num_year_cols);
%  - 在两行表头之间增加一条细线以改善视觉效果
fprintf(fileID, '\\cmidrule(lr){2-%d}\n', num_year_cols + 1);

% --- 第二行表头: 具体的年份 ---
year_header_cells = cellfun(@(y) num2str(y), num2cell(assessment_years), 'UniformOutput', false);
year_header_line = strjoin(year_header_cells, ' & ');
fprintf(fileID, '年龄组 & %s \\\\\n', year_header_line);
fprintf(fileID, '\\midrule\n');


% --- 插入工作期和退休期子标题 ---
fprintf(fileID, '\\multicolumn{%d}{l}{\\textit{工作期}} \\\\\n', num_year_cols + 1);

% --- 循环写入工作期数据 ---
for i = 1:size(working_age_groups, 1)
    temp_row_cell = cell(1, width(results_table));
    for j = 1:width(results_table)
        value = results_table{i, j};
        temp_row_cell{j} = convert_to_string(value);
    end
    row_data_str = strjoin(temp_row_cell, ' & ');
    fprintf(fileID, '%s \\\\\n', row_data_str);
end

fprintf(fileID, '\\midrule\n');
fprintf(fileID, '\\multicolumn{%d}{l}{\\textit{退休期}} \\\\\n', num_year_cols + 1);

% --- 循环写入退休期数据 ---
for i = 1:size(retired_age_groups, 1)
    row_idx = size(working_age_groups, 1) + i;
    temp_row_cell = cell(1, width(results_table));
    for j = 1:width(results_table)
        value = results_table{row_idx, j};
        temp_row_cell{j} = convert_to_string(value);
    end
    row_data_str = strjoin(temp_row_cell, ' & ');
    fprintf(fileID, '%s \\\\\n', row_data_str);
end

% --- 写入表格结尾和注释 ---
fprintf(fileID, '\\bottomrule\n');
fprintf(fileID, '\\end{tabular}\n');
fprintf(fileID, '\\begin{tablenotes}[para,flushleft]\n');
fprintf(fileID, '  \\footnotesize\n');
fprintf(fileID, '  \\item[注] 表中数值为“有PPS情景”相对于“无PPS情景”的补偿变化(CV)。\n');
fprintf(fileID, '  CV值为正，代表该年龄组在评估年份的剩余生命周期福利因引入PPS制度而提高。例如，+1.00\\%%表示该群体愿意放弃其剩余生命周期消费的1\\%%，以换取从“无PPS”制度切换到“有PPS”制度。计算基于各年龄组在该时刻的期望剩余终身效用。\n');
fprintf(fileID, '\\end{tablenotes}\n');
fprintf(fileID, '\\end{threeparttable}\n');
fprintf(fileID, '\\end{table}\n');

% --- 关闭文件 ---
fclose(fileID);
fprintf('   ✅ 成功生成 LaTeX 表格文件: %s\n', output_tex_filename);
fprintf('\n--- 福利分析脚本执行完毕 ---\n');


% --- [辅助函数] ---
function str_out = convert_to_string(value_in)
    % 一个稳健的函数，将不同类型的数据转换为字符串
    if isstring(value_in)
        str_out = char(value_in);
    elseif ischar(value_in)
        str_out = value_in;
    elseif isnumeric(value_in) || islogical(value_in)
        str_out = num2str(value_in);
    elseif iscell(value_in) % 如果输入是cell，递归处理
        str_out = convert_to_string(value_in{1});
    else
        str_out = 'Unsupported Type';
    end
end