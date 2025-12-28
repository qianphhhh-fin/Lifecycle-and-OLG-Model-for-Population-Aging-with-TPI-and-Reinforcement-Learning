%% 敏感性分析结果可视化与制表脚本 (已修改)
%
% 描述:
%   本脚本旨在自动读取所有已生成的敏感性分析结果(.mat文件)，
%   为每种参数类别生成一组具有经济学解释力的对比图表或LaTeX表格，
%   并严格按照博士论文的排版要求进行格式化。
%
% 流程:
%   1. 定义要进行可视化分析的参数列表，区分绘图与制表。
%   2. 对每个参数，自动在 'result/sensitivity/' 目录中查找相关结果文件。
%   3. 加载基准模型和各敏感性模型的结果，并为每个模型重新运行一次
%      标准化的蒙特卡洛模拟以获取完整的生命周期路径数据。
%   4. 对gamma参数，调用定制的绘图函数生成灰度图。
%   5. 对其他参数，动态发现存在的参数变动值，并调用统一的制表函数，
%      为三个核心决策变量(alpha, Q, C/W)分别生成符合规范的LaTeX表格。
%   6. 将生成的图像和表格保存到指定目录。
%
% 前置条件:
%   - result/vfi_results.mat (基准模型) 文件已存在。
%   - result/sensitivity/ 文件夹下所有敏感性分析的 .mat 文件已生成。

% --- 初始化环境 ---
clear;
close all;
clc;
fprintf('===== 敏感性分析结果可视化与制表脚本启动 =====\n\n');

%% --- 1. 全局配置 ---
% 定义要分析的参数及其配置
% 'param_name': 参数在cS结构体中的字段名
% 'latex_symbol': 在表格中显示的LaTeX符号
% 'is_table_param': 标记此参数是用于制表 (true)还是绘图 (false)
param_configs = {
    struct('param_name', 'gamma', 'is_table_param', false);
    struct('param_name', 'rp', 'sub_graph_name', '养老金收益率',   'is_table_param', true, 'latex_symbol', '\bar{R}_p');
    struct('param_name', 'sigr',  'sub_graph_name', '股票收益波动率','is_table_param', true, 'latex_symbol', '\sigma_r');
    struct('param_name', 'smav',  'sub_graph_name', '持久性收入冲击波动率','is_table_param', true, 'latex_symbol', '\sigma_z');
    struct('param_name', 'tr',    'sub_graph_name', '退休年龄','is_table_param', true, 'latex_symbol', 'T_R');
};

% --- 路径和制表/绘图常量 ---
BASE_RESULT_FILE = 'result/vfi_results.mat';
SENSITIVITY_DIR = 'result/sensitivity/';
IMAGE_OUTPUT_DIR = 'tex/fig/ch04/sensitivity/';
TABLE_OUTPUT_DIR = 'tex/tab/ch04/';
FIGURE_WIDTH_CM = 14.99; 

if ~exist(IMAGE_OUTPUT_DIR, 'dir'), mkdir(IMAGE_OUTPUT_DIR); end
if ~exist(TABLE_OUTPUT_DIR, 'dir'), mkdir(TABLE_OUTPUT_DIR); end

%% --- 2. 主循环：加载所有数据 ---
all_sensitivity_data = struct();

for i = 1:length(param_configs)
    config = param_configs{i};
    param_name = config.param_name;
    
    fprintf('--- 开始处理参数: %s ---\n', upper(param_name));
    
    % --- a. 发现相关文件 ---
    search_pattern = sprintf('vfi_results_%s_*.mat', param_name);
    sensitivity_files = dir(fullfile(SENSITIVITY_DIR, search_pattern));
    
    if isempty(sensitivity_files)
        warning('未找到参数 %s 的任何敏感性分析文件，跳过此参数。', param_name);
        continue;
    end
    
    file_list = {BASE_RESULT_FILE};
    for j = 1:length(sensitivity_files)
        file_list{end+1} = fullfile(SENSITIVITY_DIR, sensitivity_files(j).name);
    end
    
    % --- b. 加载数据并运行模拟 ---
    all_results_for_param = [];
    values = [];
    
    for j = 1:length(file_list)
        file_path = file_list{j};
        
        % 从文件名解析参数值
        if j == 1 % 基准模型
            temp_data = load(file_path, 'cS');
            current_value = getfield(temp_data.cS, param_name);
            label = sprintf('Baseline (%.3g)', current_value);
            is_baseline = true;
        else
            [~, fname, ~] = fileparts(file_path);
            token = regexp(fname, [param_name, '_([\d_p]+)'], 'tokens');
            value_str = token{1}{1};
            current_value = str2double(strrep(value_str, 'p', '.'));
            label = sprintf('%.4g', current_value);
            is_baseline = false;
        end
        
        fprintf('     - 模拟: %s = %s\n', param_name, label);
        
        metrics = run_simulation_for_file(file_path);
        
        values(j) = current_value;
        all_results_for_param(j).label = label;
        all_results_for_param(j).value = current_value;
        all_results_for_param(j).metrics = metrics;
        all_results_for_param(j).is_baseline = is_baseline;
        all_results_for_param(j).ages = (metrics.cS.tb : metrics.cS.td)';
    end
    
    % --- c. 排序并存储 ---
    [~, sort_idx] = sort(values);
    all_sensitivity_data.(param_name) = all_results_for_param(sort_idx);
end
fprintf('--- 所有参数数据加载与模拟完成 ---\n\n');

%% --- 3. 生成输出：绘图与制表 ---

% --- a. 为 gamma 参数生成灰度图 ---
fprintf('正在为 gamma 生成灰度图...\n');
plot_gamma_sensitivity(all_sensitivity_data.gamma, FIGURE_WIDTH_CM, IMAGE_OUTPUT_DIR);
fprintf('  ✅ 图表 sensitivity_gamma.png 已生成。\n\n');

% --- b. 为其他参数生成三张汇总表格 ---
fprintf('正在生成敏感性分析表格...\n');

% 准备制表所需的数据和配置
table_params_data = struct();
table_params_config = {};
for i = 1:length(param_configs)
    if param_configs{i}.is_table_param
        param_name = param_configs{i}.param_name;
        table_params_data.(param_name) = all_sensitivity_data.(param_name);
        table_params_config{end+1} = param_configs{i};
    end
end

% 定义三个指标的表格信息
metrics_to_tabulate = {
    struct('name', 'alpha', 'caption', '风险资产配置比例的敏感性分析', 'unit', '', 'file_suffix', 'alpha', 'is_q_table', false);
    struct('name', 'Q_real', 'caption', '个人养老金年缴费的敏感性分析', 'unit', ' (万元)', 'file_suffix', 'q', 'is_q_table', true);
    struct('name', 'c_ratio', 'caption', '消费/现金持有之比的敏感性分析', 'unit', '', 'file_suffix', 'cw', 'is_q_table', false);
};

% 调用制表函数
for i = 1:length(metrics_to_tabulate)
    generate_sensitivity_table(table_params_data, table_params_config, metrics_to_tabulate{i}, TABLE_OUTPUT_DIR);
    fprintf('  ✅ 表格 sensitivity_%s.tex 已生成。\n', metrics_to_tabulate{i}.file_suffix);
end


fprintf('\n===== 所有敏感性分析图表与表格已成功生成。 =====\n');


%% ========================================================================
%                     本地辅助函数：模拟、绘图、制表
% =========================================================================

function metrics = run_simulation_for_file(file_path)
    % 描述: 加载指定的.mat文件，严格按照标准流程运行模拟，并计算所有
    %       可能用到的绘图/制表指标，返回一个包含这些指标的结构体。
    
    data = load(file_path, 'vfi_results', 'cS');
    vfi_results = data.vfi_results;
    cS = data.cS;

    rng(42); 
    nsim = 10000; 
    C_policy = vfi_results.C_policy;
    A_policy = vfi_results.A_policy;
    Q_policy = vfi_results.Q_policy;
    simGPY = ones(cS.tn, nsim);       
    simY_norm = zeros(cS.tn, nsim);    
    simR = zeros(cS.tn, nsim);         
    simW_norm = zeros(cS.tn, nsim);    
    simF_norm = zeros(cS.tn, nsim);    
    simC_norm = zeros(cS.tn, nsim);    
    simQ_norm = zeros(cS.tn, nsim);    
    simA = zeros(cS.tn, nsim);         
    sim_norm_cash = zeros(cS.tn, nsim);

    work_periods = cS.tr - cS.tb;
    for i1 = 1:floor(nsim/2)
        z_shocks = randn(work_periods, 1) * cS.smav;
        u_shocks_log = randn(work_periods, 1) * cS.smay;
        r_shocks = randn(cS.tn, 1) * cS.sigr;
        simGPY(2:work_periods+1, i1) = cS.f_y(2:work_periods+1) ./ cS.f_y(1:work_periods) .* exp(z_shocks);
        simY_norm(1:work_periods, i1) = exp(u_shocks_log);
        simR(:, i1) = cS.rf + cS.mu + r_shocks;
        i2 = nsim/2 + i1;
        simGPY(2:work_periods+1, i2) = cS.f_y(2:work_periods+1) ./ cS.f_y(1:work_periods) .* exp(-z_shocks);
        simY_norm(1:work_periods, i2) = exp(-u_shocks_log);
        simR(:, i2) = cS.rf + cS.mu - r_shocks;
    end
    simY_norm(work_periods+1:cS.tn, :) = cS.ret_fac;

    for t = 1:cS.tn
        if t <= work_periods 
            norm_cash = simW_norm(t, :) + (1-cS.tau_y) * simY_norm(t, :);
        else 
            norm_pension_payout = (1-cS.tau_q) * simF_norm(t, :) * cS.pension_rate;
            norm_cash = simW_norm(t, :) + simY_norm(t, :) + norm_pension_payout;
        end
        sim_norm_cash(t, :) = norm_cash;
        
        Fc_interpolant = griddedInterpolant({cS.gcash, cS.gfund}, C_policy(:,:,t), 'spline', 'linear');
        Fa_interpolant = griddedInterpolant({cS.gcash, cS.gfund}, A_policy(:,:,t), 'spline', 'linear');
        Fq_interpolant = griddedInterpolant({cS.gcash, cS.gfund}, Q_policy(:,:,t), 'spline', 'linear');

        simC_norm(t, :) = Fc_interpolant(norm_cash, simF_norm(t, :));
        simA(t, :) = Fa_interpolant(norm_cash, simF_norm(t, :));
        simQ_norm(t, :) = Fq_interpolant(norm_cash, simF_norm(t, :));

        simA(t, isnan(simA(t,:))) = 0;
        simA(t, :) = max(min(simA(t, :), 1), 0);
        simC_norm(t, :) = max(min(simC_norm(t, :), norm_cash), 0);
        simQ_norm(t, :) = max(min(simQ_norm(t, :), cS.Q_max), 0);
        
        if t <= work_periods
            total_outflow = simC_norm(t, :) + (1-cS.tau_y)*simQ_norm(t, :);
            idx_exceed = total_outflow > norm_cash;
            if any(idx_exceed)
                scale_factor = norm_cash(idx_exceed) ./ total_outflow(idx_exceed) * 0.9999;
                simC_norm(t, idx_exceed) = simC_norm(t, idx_exceed) .* scale_factor;
                simQ_norm(t, idx_exceed) = simQ_norm(t, idx_exceed) .* scale_factor;
            end
            liquid_sav_norm = norm_cash - simC_norm(t, :) - (1-cS.tau_y)*simQ_norm(t, :);
        else 
            simQ_norm(t,:) = 0;
            simC_norm(t,:) = min(simC_norm(t,:), norm_cash * 0.9999);
            liquid_sav_norm = norm_cash - simC_norm(t, :);
        end

        if t < cS.tn
            simS_norm = simA(t, :) .* liquid_sav_norm;
            simB_norm = liquid_sav_norm - simS_norm;
            portfolio_return_next = simB_norm * cS.rf + simS_norm .* simR(t, :);
            simW_norm(t+1, :) = portfolio_return_next ./ simGPY(t+1, :);
            if t <= work_periods
                simF_norm(t+1, :) = ((simF_norm(t, :) + simQ_norm(t, :)) * cS.rp) ./ simGPY(t+1, :);
            else
                simF_norm(t+1, :) = simF_norm(t, :);
            end
        end
    end

    initial_income_yuan = 60000;
    simP_real = zeros(cS.tn, nsim);
    simP_real(1, :) = initial_income_yuan;
    for t = 1:(cS.tn - 1)
        simP_real(t+1, :) = simP_real(t, :) .* simGPY(t+1, :);
    end
    simQ_real = simQ_norm .* simP_real;

    sim_c_ratio = simC_norm ./ sim_norm_cash;
    sim_c_ratio(isinf(sim_c_ratio) | isnan(sim_c_ratio)) = NaN;

    metrics.cS = cS;
    metrics.med_alpha = nanmedian(simA, 2);
    metrics.med_c_ratio = nanmedian(sim_c_ratio, 2);
    metrics.med_Q_real = nanmedian(simQ_real, 2);
end

% =========================================================================

function plot_gamma_sensitivity(results, fig_width_cm, output_dir)
    % 为gamma参数生成灰度图
    fig = figure('Name', 'Gamma Sensitivity');
    set(fig, 'Units', 'centimeters', 'Position', [5, 5, fig_width_cm, 5]);
    
    % 使用灰度颜色 (从浅灰到黑)
    num_colors = length(results);
    colors = repmat(linspace(0.7, 0, num_colors)', 1, 3);
    
    retirement_age = results(1).metrics.cS.tr;
    tb = results(1).metrics.cS.tb;
    yuan_to_wanyuan = 1/10000;

    % 子图 (a): 风险资产配置比例
    ax1 = subplot(1, 3, 1);
    hold on;
    for i = 1:length(results)
        plot(results(i).ages, results(i).metrics.med_alpha, 'LineWidth', 1.5, 'Color', colors(i,:));
    end
    xline(retirement_age, '--', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
    hold off;
    title('(a) $\alpha_t$', 'Interpreter', 'latex'); xlabel('年龄'); 
    ylabel('比例'); 
    grid on; xlim([tb, 100]); ylim([0, 1]);
    
    % 子图 (b): 个人养老金缴费绝对值
    ax2 = subplot(1, 3, 2);
    work_periods = retirement_age - tb;
    hold on;
    for i = 1:length(results)
        plot(results(i).ages(1:work_periods), results(i).metrics.med_Q_real(1:work_periods) * yuan_to_wanyuan, 'LineWidth', 1.5, 'Color', colors(i,:));
    end
    hold off;
    title('(b) $Q_t$', 'Interpreter', 'latex'); xlabel('年龄'); 
    ylabel('绝对数额'); 
    grid on; xlim([tb, retirement_age]);
    
    % 子图 (c): 消费比例
    ax3 = subplot(1, 3, 3);
    hold on;
    for i = 1:length(results)
        plot(results(i).ages, results(i).metrics.med_c_ratio, 'LineWidth', 1.5, 'Color', colors(i,:), 'DisplayName', ['$\gamma = $ ', results(i).label]);
    end
    xline(retirement_age, '--', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
    hold off;
    title('(c) $C_t/W_t$', 'Interpreter', 'latex'); xlabel('年龄');
    ylabel('比例'); grid on; xlim([tb, 100]);
    
    lgd = legend(ax3, 'show','box','off','Location', 'best', 'FontSize', 7, 'Interpreter', 'latex');
    lgd.ItemTokenSize(1) = 10;
    
    output_filename = fullfile(output_dir, 'sensitivity_gamma.png');
    print(fig, output_filename, '-dpng', '-r400');
end

% =========================================================================

function generate_sensitivity_table(all_params_data, all_params_config, metric_info, output_dir)
    % 描述: 为指定指标(alpha, Q, or C/W)生成一个完整的LaTeX敏感性分析表格
    
    % --- 1. 定义年龄段 ---
    tb = all_params_data.rp(1).metrics.cS.tb;
    td = all_params_data.rp(1).metrics.cS.td;
    
    age_bins = {
        [tb, 29], [30, 39], [40, 49], [50, 60]
    };
    age_bin_labels = {'22-29', '30-39', '40-49', '50-60'};

    if ~metric_info.is_q_table
        age_bins = [age_bins, {[61, 70], [71, 80], [81, 90], [91, td]}];
        age_bin_labels = [age_bin_labels, {'61-70', '71-80', '81-90', '91-100'}];
    end

    % --- 2. 预计算所有格子的数据 ---
    processed_data = cell(length(all_params_config), 1);
    discovered_variations = cell(length(all_params_config), 1);

    for p_idx = 1:length(all_params_config)
        config = all_params_config{p_idx};
        param_name = config.param_name;
        param_results_all = all_params_data.(param_name);
        
        % 动态发现参数变动值 (非基准, 最多取前3个)
        param_results_variations = param_results_all(~[param_results_all.is_baseline]);
        num_variations = min(3, length(param_results_variations));
        if num_variations < 3
             warning('参数 %s 找到的变动少于3个，表格将只显示 %d 列', param_name, num_variations);
        end
        param_results = param_results_variations(1:num_variations);
        table_values = [param_results.value];
        discovered_variations{p_idx} = table_values;
        
        param_data = zeros(length(age_bins), num_variations);

        for v_idx = 1:num_variations
            result = param_results(v_idx);
            
            % 获取对应的指标序列
            switch metric_info.name
                case 'alpha'
                    metric_series = result.metrics.med_alpha;
                case 'Q_real'
                    metric_series = result.metrics.med_Q_real / 10000; % 单位万元
                case 'c_ratio'
                    metric_series = result.metrics.med_c_ratio;
            end
            
            % 按年龄段计算均值
            for b_idx = 1:length(age_bins)
                bin = age_bins{b_idx};
                age_indices = (result.ages >= bin(1) & result.ages <= bin(2));
                param_data(b_idx, v_idx) = nanmean(metric_series(age_indices));
            end
        end
        processed_data{p_idx} = param_data;
    end
    
    % --- 3. 写入LaTeX文件 ---
    filename = fullfile(output_dir, ['sensitivity_', metric_info.file_suffix, '.tex']);
    fid = fopen(filename, 'w', 'n', 'UTF-8');
    
    fprintf(fid, '%% Auto-generated by MATLAB script on %s\n\n', datestr(now));
    
    % 表格头部
    fprintf(fid, '\\begin{table}[h!]\n');
    fprintf(fid, '\\centering\n');
    fprintf(fid, '\\caption{%s}\n', metric_info.caption);
    fprintf(fid, '\\label{tab:sensitivity-%s}\n', metric_info.file_suffix);
    fprintf(fid, '\\small\n');
    fprintf(fid, '\\renewcommand{\\arraystretch}{1.1}\n');
    fprintf(fid, '\\makebox[\\textwidth][c]{\n');
    fprintf(fid, '\\begin{tabular}{cc}\n');
    
    % 写入两行子表
    write_table_row(fid, all_params_config(1:2), processed_data(1:2), discovered_variations(1:2), age_bin_labels,{'(a)','(b)'});
    fprintf(fid, '\\\\[0.8cm]\n'); % 两行之间的垂直间距
    write_table_row(fid, all_params_config(3:4), processed_data(3:4), discovered_variations(3:4), age_bin_labels, {'(c)','(d)'});

    % 表格尾部
    fprintf(fid, '\\end{tabular}\n');
    fprintf(fid, '}\n');
    fprintf(fid, '\\vspace{0.2cm}\n');
    fprintf(fid, '\\begin{tablenotes}\n');
    fprintf(fid, '\\small\n');
    fprintf(fid, '\\item 注：表中数值表示不同年龄段决策变量的均值%s。基准模型参数值未在表中汇报。\n', metric_info.unit);
    fprintf(fid, '\\end{tablenotes}\n');
    fprintf(fid, '\\end{table}\n');
    
    fclose(fid);
end

function write_table_row(fid, configs, data, variations, age_labels, prefix)
    % 辅助函数：向文件写入一行包含两个子表格的LaTeX代码

    % --- 子表1 ---
    conf1 = configs{1};
    data1 = data{1};
    vars1 = variations{1};
    col_format_1 = ['c|', repmat('c', 1, length(vars1))];
    header_1 = sprintf('%.4g & ', vars1); header_1 = header_1(1:end-2); % Remove trailing '& '

    fprintf(fid, '\\begin{tabular}{%s}\n', col_format_1);
    fprintf(fid, '\\multicolumn{%d}{c}{%s %s $%s$}\\\\\n', length(vars1)+1,prefix{1}, conf1.sub_graph_name, conf1.latex_symbol);
    fprintf(fid, '\\toprule\n');
    fprintf(fid, '年龄段 & %s \\\\\n', header_1);
    fprintf(fid, '\\midrule\n');
    for i = 1:length(age_labels)
        row_data_1 = sprintf('%.3f & ', data1(i,:)); row_data_1 = row_data_1(1:end-2);
        fprintf(fid, '%s & %s \\\\\n', age_labels{i}, row_data_1);
    end
    fprintf(fid, '\\bottomrule\n');
    fprintf(fid, '\\end{tabular}\n');
    
    fprintf(fid, '&\n'); % 两个子表之间的分隔符
    
    % --- 子表2 ---
    conf2 = configs{2};
    data2 = data{2};
    vars2 = variations{2};
    col_format_2 = ['c|', repmat('c', 1, length(vars2))];
    header_2 = sprintf('%.4g & ', vars2); header_2 = header_2(1:end-2);

    fprintf(fid, '\\begin{tabular}{%s}\n', col_format_2);
    fprintf(fid, '\\multicolumn{%d}{c}{%s %s $%s$}\\\\\n', length(vars2)+1, prefix{2}, conf2.sub_graph_name, conf2.latex_symbol);
    fprintf(fid, '\\toprule\n');
    fprintf(fid, '年龄段 & %s \\\\\n', header_2);
    fprintf(fid, '\\midrule\n');
    for i = 1:length(age_labels)
        row_data_2 = sprintf('%.3f & ', data2(i,:)); row_data_2 = row_data_2(1:end-2);
        fprintf(fid, '%s & %s \\\\\n', age_labels{i}, row_data_2);
    end
    fprintf(fid, '\\bottomrule\n');
    fprintf(fid, '\\end{tabular}\n');
end