% =========================================================================
% == SCRIPT: ss_comparison.m
% == 版本: [v2.2 - PAYG报告内容修正版]
% ==
% == 目的:
% ==   1. 加载单个转轨路径模拟的结果文件。
% ==   2. 提取并比较该路径的"初期稳态(ss0)"与"终期稳态(ssF)"。
% ==   3. 报告PAYG体系的规则化福利、总缴费和均衡调整因子(adj)。
% ==   4. 自动生成一个可直接用于论文的LaTeX表格 (Tab4.0_ss_comparison.tex)。
% ==
% == 前置条件:
% ==   - 已成功运行 main_run_trans.m 并生成以下文件:
% ==     - TRANS/TPI_results_het_nopps.mat
% =========================================================================

clear; close all; clc;
addpath(pwd);
fprintf('=== 转轨路径稳态比较与LaTeX表格生成脚本 (v2.2) ===\n\n');

%% --- 1. 加载数据 ---
fprintf('--- 1. 加载转轨路径的稳态数据 ---\n');

% --- 文件路径定义 ---
file_path = 'TRANS/TPI_results_het_nopps.mat';
output_tex_filename = 'tex/tab/Tab4.0_ss_comparison.tex';

% --- 加载数据 ---
try
    loaded_data = load(file_path, 'results', 'cS', 'ss0', 'ssF');
    results = loaded_data.results;
    ss0 = loaded_data.ss0;
    ssF = loaded_data.ssF;
    cS = loaded_data.cS;
    paramS_for_report = results.ss_data.paramSF;

    % --- 为核算报告重构稳态人口分布 (Z_ss_norm) 和总人口 ---
    Z0_total_abs = cS.Z_path_raw(:, 1);
    Z_ss_norm_0 = (Z0_total_abs / sum(Z0_total_abs)) * cS.type_weights';
    mass_total_0 = sum(Z_ss_norm_0, 'all');
    
    ZF_total_abs = cS.Z_path_raw(:, end);
    Z_ss_norm_F = (ZF_total_abs / sum(ZF_total_abs)) * cS.type_weights';
    mass_total_F = sum(Z_ss_norm_F, 'all');

    fprintf('   ✅ 已从 "%s" 加载转轨数据、初始稳态(ss0)和终期稳态(ssF)。\n', file_path);
catch ME
    error('加载转轨数据失败: %s', ME.message);
end

fprintf('   数据加载完成。\n');

%% --- 2. 生成国民经济核算报告 ---
fprintf('\n--- 2. 生成国民经济核算报告 ---\n');

% --- 为ss0和ssF分别生成报告结构体 (关闭命令行打印) ---
report_ss0 = SS.display_national_accounts(ss0, cS, paramS_for_report, Z_ss_norm_0, '', false, false);
report_ssF = SS.display_national_accounts(ssF, cS, paramS_for_report, Z_ss_norm_F, '', false, true);

% --- [核心修改] 手动计算规则化总福利，并添加到报告结构体中 ---
total_formula_benefit_ss0 = 0;
aR_idx_0 = cS.aR_new_path(1);
for h_idx = 1:cS.nTypes
    mass_retirees_h = sum(Z_ss_norm_0((aR_idx_0+1):end, h_idx));
    total_formula_benefit_ss0 = total_formula_benefit_ss0 + ss0.b_hat_formula_h(h_idx) * mass_retirees_h;
end
report_ss0.PAYG_formula_benefit = total_formula_benefit_ss0;

total_formula_benefit_ssF = 0;
aR_idx_F = cS.aR_new_path(end);
for h_idx = 1:cS.nTypes
    mass_retirees_h = sum(Z_ss_norm_F((aR_idx_F+1):end, h_idx));
    total_formula_benefit_ssF = total_formula_benefit_ssF + ssF.b_hat_formula_h(h_idx) * mass_retirees_h;
end
report_ssF.PAYG_formula_benefit = total_formula_benefit_ssF;

% --- 为报告结构体补充额外的比率 ---
report_ss0.Ig_Y_ratio = (report_ss0.I_g / report_ss0.Y) * 100;
report_ssF.Ig_Y_ratio = (report_ssF.I_g / report_ssF.Y) * 100;
report_ss0.Gc_Y_ratio = (report_ss0.G_c / report_ss0.Y) * 100;
report_ssF.Gc_Y_ratio = (report_ssF.G_c / report_ssF.Y) * 100;

fprintf('   ✅ 核算报告已在内存中生成并补充完毕。\n');


%% --- 3. [核心] 定义比较变量并计算差异 ---
fprintf('\n--- 3. 定义比较变量并计算表格内容 ---\n');

% --- [核心修改] 定义要比较的变量和格式 ---
vars_to_compare = {
    '\multicolumn{4}{l}{\textbf{A. 有效人均宏观总量}}', '', '', '', '';
    '产出 ($y$)', 'Y', '%.4f', 'per_capita', '';
    '私人资本 ($k_p$)', 'K_p', '%.4f', 'per_capita', '';
    '公共资本 ($k_g$)', 'K_g', '%.4f', 'per_capita', '';
    '私人消费 ($c$)', 'C', '%.4f', 'per_capita', '';
    '私人投资 ($i_p$)', 'I_p_acct', '%.4f', 'per_capita', '';
    '\multicolumn{4}{l}{\textbf{B. 核心宏观比率}}', '', '', '', '';
    '私人资本/产出比 ($K_p/Y$)', 'Kp_Y_ratio', '%.3f', 'diff', '';
    '消费/产出比 ($C/Y$)', 'C_Y_ratio', '%.2f', 'ratio_pct', '\%';
    '私人投资/产出比 ($I_p/Y$)', 'Ip_Y_ratio', '%.2f', 'ratio_pct', '\%';
    '公共投资/产出比 ($I_g/Y$)', 'Ig_Y_ratio', '%.2f', 'ratio_pct', '\%';
    '政府消费/产出比 ($G_c/Y$)', 'Gc_Y_ratio', '%.2f', 'ratio_pct', '\%';
    '\multicolumn{4}{l}{\textbf{C. 要素价格}}', '', '', '', '';
    '年化真实利率 ($r$)', 'r_annual', '%.2f', 'diff', '\%';
    '有效工资率 ($w$)', 'w', '%.4f', 'abs_change_pct', '';
    '\multicolumn{4}{l}{\textbf{D. PAYG养老金体系}}', '', '', '', '';
    '总缴费 (占Y比例)', 'PAYG_contrib', '%.2f', 'ratio_y', '\%';
    '规则化总福利 (占Y比例)', 'PAYG_formula_benefit', '%.2f', 'ratio_y', '\%';
    '均衡调整因子 (adj)', 'adj', '%.4f', 'diff_abs', '';
};


% --- 初始化存储结果的 cell 数组 ---
table_data = cell(size(vars_to_compare, 1), 4); % {Label, Val_ss0, Val_ssF, Change}

% --- 循环处理每个变量 ---
for i = 1:size(vars_to_compare, 1)
    label = vars_to_compare{i, 1};
    var_name = vars_to_compare{i, 2};
    
    if isempty(var_name), table_data{i, 1} = label; continue; end
    
    table_data{i, 1} = label;
    
    format_spec = vars_to_compare{i, 3};
    change_type = vars_to_compare{i, 4};
    unit = vars_to_compare{i, 5};
    
    % --- 提取数值 ---
    if strcmp(var_name, 'adj')
        val_ss0 = ss0.adj; val_ssF = ssF.adj;
    else
        if isfield(report_ss0, var_name), val_ss0 = report_ss0.(var_name); else, val_ss0 = NaN; end
        if isfield(report_ssF, var_name), val_ssF = report_ssF.(var_name); else, val_ssF = NaN; end
    end
    
    % --- 根据类型进行预处理 ---
    if strcmp(change_type, 'per_capita')
        val_ss0 = val_ss0 / mass_total_0;
        val_ssF = val_ssF / mass_total_F;
    elseif strcmp(change_type, 'ratio_y')
        if isfield(report_ss0, 'Y') && report_ss0.Y~=0, val_ss0 = (val_ss0 / report_ss0.Y)*100; end
        if isfield(report_ssF, 'Y') && report_ssF.Y~=0, val_ssF = (val_ssF / report_ssF.Y)*100; end
    end
    
    % --- 格式化数值列 ---
    if isnan(val_ss0), table_data{i, 2} = '-'; else, table_data{i, 2} = [sprintf(format_spec, val_ss0), unit]; end
    if isnan(val_ssF), table_data{i, 3} = '-'; else, table_data{i, 3} = [sprintf(format_spec, val_ssF), unit]; end

    % --- 计算并格式化变化列 ---
    change_val_str = '-';
    if ~isnan(val_ss0) && ~isnan(val_ssF)
        switch change_type
            case {'per_capita', 'abs_change_pct'} % 相对变化
                if abs(val_ss0) > 1e-9, change = (val_ssF / val_ss0 - 1) * 100; else, change = NaN; end
                if ~isnan(change), change_val_str = sprintf('%.2f\\%%', change); end
            case {'diff', 'ratio_y', 'ratio_pct'} % 百分点差异
                change = val_ssF - val_ss0;
                change_val_str = sprintf('%.2f pp', change);
            case 'diff_abs' % 绝对值差异
                change = val_ssF - val_ss0;
                change_val_str = sprintf('%+.4f', change);
        end
    end
    table_data{i, 4} = change_val_str;
end
fprintf('   ✅ 表格内容计算完成。\n');

%% --- 4. 生成 LaTeX 表格 ---
fprintf('\n--- 4. 生成 LaTeX 格式的稳态比较表格 ---\n');

% --- 准备输出目录和文件 ---
[output_dir, ~, ~] = fileparts(output_tex_filename);
if ~exist(output_dir, 'dir'), mkdir(output_dir); end
fileID = fopen(output_tex_filename, 'w', 'n', 'UTF-8');
if fileID == -1, error('无法创建或打开LaTex输出文件。'); end

% --- 写入表格框架 ---
fprintf(fileID, '%% =======================================================\n');
fprintf(fileID, '%%  此文件由 ss_comparison.m (v2.2) 自动生成\n');
fprintf(fileID, '%% =======================================================\n');
fprintf(fileID, '\\begin{table}[htbp]\n');
fprintf(fileID, '\\centering\n');
fprintf(fileID, '\\caption{转轨路径的初期与终期稳态宏观经济结果比较}\n');
fprintf(fileID, '\\label{tab4.0:ss_comparison}\n');
fprintf(fileID, '\\begin{threeparttable}\n');
fprintf(fileID, '\\begin{tabular}{lccc}\n');
fprintf(fileID, '\\toprule\n');

% [核心修改] 更新表头以包含动态年份
fprintf(fileID, '变量 & 起点 (%d年) & 终期稳态 (%d年) & 变化 \\\\\n', cS.start_year, cS.end_year);

fprintf(fileID, '\\midrule\n');

% --- 循环写入表格内容 ---
for i = 1:size(table_data, 1)
    if isempty(vars_to_compare{i, 2}) % 如果是子标题
        fprintf(fileID, '%s \\\\\n', table_data{i, 1});
    else
        row_str = strjoin(table_data(i,:), ' & ');
        fprintf(fileID, '%s \\\\\n', row_str);
    end
end

% --- 写入表格结尾和注释 ---
fprintf(fileID, '\\bottomrule\n');
fprintf(fileID, '\\end{tabular}\n');
fprintf(fileID, '\\begin{tablenotes}[para,flushleft]\n');
fprintf(fileID, '  \\footnotesize\n');
fprintf(fileID, '  \\item[注] “期初稳态”对应模型在转轨路径第一期(t=1)的人口与政策环境下的伪稳态解。“终期稳态”则对应路径最后一期(t=T)的BGP稳态解。\n');
fprintf(fileID, '  除特别注明外，所有流量和存量均为标准化值。财政与养老金流量已表示为占各自稳态产出(Y)的百分比。\n');
fprintf(fileID, '  均衡调整因子(adj)是在保持缴费率固定的前提下，为保证PAYG体系收支平衡所必需的福利调整系数。adj=1.0表示福利可100\\%%兑现。\n');
fprintf(fileID, '  “变化”列中，“\\%%”表示相对变化率，“pp”表示百分点差异。\n');
fprintf(fileID, '\\end{tablenotes}\n');
fprintf(fileID, '\\end{threeparttable}\n');
fprintf(fileID, '\\end{table}\n');

% --- 关闭文件 ---
fclose(fileID);
fprintf('   ✅ 成功生成 LaTeX 表格文件: %s\n', output_tex_filename);
fprintf('\n--- 稳态比较脚本执行完毕 ---\n');