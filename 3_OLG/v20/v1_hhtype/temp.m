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