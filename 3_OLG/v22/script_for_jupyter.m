% =========================================================================
% == SCRIPT: tab_pps.m
% == 版本: [v1.4 - 引入缓存机制]
% ==
% == 核心任务:
% ==   1. [新增] 引入缓存机制，将耗时计算结果保存，避免重复运行。
% ==   2. [重构] 将计算部分与图表生成部分分离。
% ==   3. [保留] 生成 表 4.5(现图4.8)、表 4.6 及附录福利表。
% =========================================================================

clear; close all; clc;
addpath(pwd);
fprintf('=== PPS改革情景对比分析脚本 (v1.4) ===\n\n');

%% --- 0. [!!!] 快速测试模式开关 ---
quick_test_mode = false;
if quick_test_mode
    fprintf('!!! 警告: 快速测试模式已开启。结果将不被缓存，且附录福利分析结果将不准确。!!!\n');
end
% =====================================================================

%% --- 1. 用户设定与文件定义 ---
fprintf('\n--- 1. 用户设定与文件定义 ---\n');

scenarios = {'None', 'conservative', 'moderate', 'ambitious'};
base_scenario_welfare = 'None';
scenarios_for_welfare_comp = {'conservative', 'moderate', 'ambitious'};

file_tpi_pattern = 'TRANS/PPS/TPI_results_het_pps_%s.mat';
file_tpi_nopps = 'TRANS/PPS/TPI_results_het_nopps.mat';

% --- 输出路径定义 ---
output_dir_tex = 'tex/tab/';
output_dir_fig = 'tex/fig/';
output_dir_results = 'results/'; % [新增] 缓存文件目录
if ~exist(output_dir_tex, 'dir'), mkdir(output_dir_tex); end
if ~exist(output_dir_fig, 'dir'), mkdir(output_dir_fig); end
if ~exist(output_dir_results, 'dir'), mkdir(output_dir_results); end

% --- 输出文件名定义 ---
output_tex_table_macro = fullfile(output_dir_tex, 'tab4.5_macro_pps.tex'); % 注意：此脚本不再生成此表
output_tex_table_adj = fullfile(output_dir_tex, 'tab4.6_adj_sustainability_pps.tex');
appendix_file_pattern = 'Appendix_tab1.1_welfare_analysis_pps_%s.tex';
output_fig_pps = fullfile(output_dir_fig, 'fig4.8_macro_pps.png');

% --- [新增] 缓存文件定义 ---
output_file_cache = fullfile(output_dir_results, 'tab_pps_processed_data.mat');
fprintf('   LaTeX 表格将输出至: %s\n', output_dir_tex);
fprintf('   图像将输出至: %s\n', output_dir_fig);
fprintf('   计算结果缓存文件为: %s\n', output_file_cache);


%% --- 2. [核心逻辑] 加载或计算处理后数据 ---
% 检查缓存文件是否存在。如果存在且不是快速测试模式，则加载；否则，执行完整计算。
if exist(output_file_cache, 'file') && ~quick_test_mode
    fprintf('\n--- 2. [缓存加载] 发现已处理数据, 直接加载... ---\n');
    load(output_file_cache, 'results_all', 'plot_data');
    fprintf('   ✅ 已加载: %s\n', output_file_cache);
else
    if quick_test_mode
        fprintf('\n--- 2. [重新计算] 快速测试模式开启, 将强制重新计算 (结果不保存) ---\n');
    else
        fprintf('\n--- 2. [重新计算] 未发现缓存文件, 开始执行完整计算... ---\n');
    end

    % --- 2.1. [数据加载] 加载所有情景的计算结果 ---
    fprintf('\n--- 2.1. [数据加载] 循环加载各方案结果 ---\n');
    results_all = struct();
    for i = 1:length(scenarios)
        scenario_name = scenarios{i};
        fprintf('   正在加载方案: [%s] ...\n', upper(scenario_name));
        if strcmp(scenario_name, 'None'), file_tpi = file_tpi_nopps;
        else, file_tpi = sprintf(file_tpi_pattern, scenario_name); end
        if ~exist(file_tpi, 'file'), error('找不到文件: %s', file_tpi); end
        results_all.(scenario_name) = load(file_tpi);
        fprintf('   ✅ 已加载: %s\n', file_tpi);
    end
    fprintf('--- 所有情景数据加载完毕 ---\n');

    % --- 2.2. [福利分析预处理] 为所有相关情景重构价值函数路径 ---
    fprintf('\n--- 2.2. [福利分析预处理] 重构价值函数路径 ---\n');
    scenarios_needing_vfi = [{base_scenario_welfare}, scenarios_for_welfare_comp];
    for i = 1:length(scenarios_needing_vfi)
        scen_name = scenarios_needing_vfi{i};
        fprintf('   正在为方案 [%s] 准备价值函数路径...\n', upper(scen_name));
        data = results_all.(scen_name);
        cS = data.cS;
        paramSF = data.results.ss_data.paramSF;
        valF_h = data.results.ss_data.valF_h;
        polF_h = data.results.ss_data.polF_h;
        iter_results = data.results.iter_results;

        if quick_test_mode
            fprintf('      >>> [快速模式] 跳过VFI, 创建占位符Val_path_h...\n');
            T = cS.T_sim;
            Val_path_h_placeholder = cell(cS.nTypes, 1);
            for h = 1:cS.nTypes
                Val_path_h_placeholder{h} = repmat(valF_h{h}, [1, 1, 1, 1, T]);
            end
            results_all.(scen_name).Val_path_h = Val_path_h_placeholder;
        else
            fprintf('      >>> [完整模式] 启动VFI反向迭代 (此过程可能较慢)...\n');
            Val_path_h = cell(cS.nTypes, 1);
            parfor h = 1:cS.nTypes % 可以考虑使用parfor加速
                pathS_h = get_pathS_for_backward(iter_results, h, cS, data.ss0);
                cS_h = get_cs_for_type(cS, h);
                [~, Val_path_h{h}] = household.backward_hh(pathS_h, cS_h, paramSF, valF_h{h}, polF_h{h});
            end
            results_all.(scen_name).Val_path_h = Val_path_h;
            fprintf('      >>> [完整模式] VFI反向迭代完成。\n');
        end
    end
    fprintf('--- 所有情景价值函数路径准备完毕 ---\n');
    
    % --- 2.3. [数据计算] 计算绘图所需数据 ---
    fprintf('\n--- 2.3. [数据计算] 为图4.8准备数据 ---\n');
    plot_data = struct();
    scenarios_for_fig_calc = {'conservative', 'moderate', 'ambitious'};
    
    % -- 计算恒定的PAYG/GDP比率 (只需计算一次) --
    data_ref = results_all.None;
    cS_ref = data_ref.cS;
    paramSF_ref = data_ref.results.ss_data.paramSF;
    t_idx_ref = 5; % 任选一个时间点
    aR_t_ref = cS_ref.aR_new_path(t_idx_ref);
    w_hat_t_ref = data_ref.results.iter_results.w_hat_path_iter(t_idx_ref);
    total_payg_contrib_ref = 0;
    for h = 1:cS_ref.nTypes
        L_h_t = TPI.get_labor_supply_path_for_type(h, {data_ref.final_Dist_path_h{h}(:,:,:,:,t_idx_ref)}, cS_ref, paramSF_ref, aR_t_ref);
        total_payg_contrib_ref = total_payg_contrib_ref + cS_ref.theta_path_h(h,t_idx_ref) * w_hat_t_ref * L_h_t;
    end
    Y_hat_t_ref = data_ref.results.Y_path(t_idx_ref) / cS_ref.A_path(t_idx_ref);
    plot_data.payg_gdp_ratio = total_payg_contrib_ref / Y_hat_t_ref;

    % -- 循环计算各PPS情景的PPS/GDP比率路径 --
    for i_scen = 1:length(scenarios_for_fig_calc)
        scen_name = scenarios_for_fig_calc{i_scen};
        data = results_all.(scen_name);
        cS = data.cS;
        paramSF = data.results.ss_data.paramSF;
        T = cS.T_sim;
        pps_gdp_ratio_path = zeros(1, T);

        for t = 1:T
            Y_hat_t = data.results.Y_path(t) / cS.A_path(t);
            aR_t = cS.aR_new_path(t);
            w_hat_t = data.results.iter_results.w_hat_path_iter(t);
            total_pps_contrib_t = 0;

            for h = 1:cS.nTypes
                cS_h = get_cs_for_type(cS, h);
                pol_t_h = data.final_Pol_path_h{h}{t};
                dist_t_h = data.final_Dist_path_h{h}(:,:,:,:,t);
                for ia = 1:aR_t
                    mass_ia_h = dist_t_h(:,:,:,ia);
                    if sum(mass_ia_h,'all') < 1e-30, continue; end
                    le_grid_slice = reshape(paramSF.leGridV, [1,1,cS.nw_expanded]);
                    labor_income_grid = w_hat_t * cS_h.ageEffV_new(ia) .* le_grid_slice;
                    labor_income_grid_full = repmat(labor_income_grid, [cS.nk, cS.nkpps, 1]);
                    pps_rate_grid = pol_t_h(ia).pps_contrib_rate;
                    total_pps_contrib_t = total_pps_contrib_t + sum(pps_rate_grid .* labor_income_grid_full .* mass_ia_h, 'all');
                end
            end

            if abs(Y_hat_t) > 1e-9, pps_gdp_ratio_path(t) = total_pps_contrib_t / Y_hat_t; end
        end
        plot_data.(scen_name).pps_gdp_ratio = pps_gdp_ratio_path;
    end
    fprintf('   ✅ 绘图数据计算完成。\n');

    % --- 2.4. [adj分析预处理] 计算各情景PPS给付总额 ---
    fprintf('\n--- 2.4. [adj分析预处理] 计算各情景PPS给付总额 ---\n');
    for i = 1:length(scenarios)
        scen_name = scenarios{i};
        data = results_all.(scen_name);
        cS = data.cS;
        T = cS.T_sim;
        total_pps_withdrawal_hat_path = zeros(1, T);
        if ~cS.pps_active
            results_all.(scen_name).total_pps_withdrawal_hat_path = total_pps_withdrawal_hat_path;
            fprintf('   方案 [%s] 无PPS, 给付额为0。\n', upper(scen_name));
            continue;
        end
        fprintf('   正在为方案 [%s] 计算PPS给付路径...\n', upper(scen_name));
        kpps_grid_vec = cS.kppsGridV;
        kpps_mat = repmat(reshape(kpps_grid_vec, [1, cS.nkpps]), [cS.nk, 1, cS.nw_expanded]);
        r_path = data.results.iter_results.r_path_iter;
        for t = 1:T
            aR_t = cS.aR_new_path(t);
            withdrawal_rate_t = cS.pps_withdrawal_rate;
            for h = 1:cS.nTypes
                dist_t_h = data.final_Dist_path_h{h}(:,:,:,:,t);
                for ia = (aR_t + 1) : cS.aD_new
                    dist_slice = dist_t_h(:, :, :, ia);
                    if sum(dist_slice, 'all') < 1e-30, continue; end
                    pps_asset_return_hat = kpps_mat * (1 + r_path(t));
                    pps_withdrawal_grid = withdrawal_rate_t * pps_asset_return_hat;
                    total_pps_withdrawal_hat_path(t) = total_pps_withdrawal_hat_path(t) + sum(pps_withdrawal_grid .* dist_slice, 'all');
                end
            end
        end
        results_all.(scen_name).total_pps_withdrawal_hat_path = total_pps_withdrawal_hat_path;
    end
    fprintf('--- 所有情景PPS给付路径计算完毕 ---\n');
    
    % --- 2.5. [保存缓存] ---
    if ~quick_test_mode
        fprintf('\n--- 2.5. [保存缓存] 所有计算完成, 保存结果到文件... ---\n');
        save(output_file_cache, 'results_all', 'plot_data', '-v7.3');
        fprintf('   ✅ 结果已保存至: %s\n', output_file_cache);
    end
end
fprintf('--- 数据准备完毕 ---\n');


%% --- 3. [生成图4.8] PPS与PAYG缴费占GDP比重路径 ---
fprintf('\n--- 3. [生成图4.8] PPS与PAYG缴费占GDP比重路径 ---\n');

% --- 1. 定义绘图参数 ---
scenarios_for_fig = {'conservative', 'moderate', 'ambitious'};
color_map = containers.Map(...
    {'conservative', 'moderate', 'ambitious'}, ...
    {[0.65 0.65 0.65], [0.0 0.0 0.0], [0.35 0.35 0.35]}); % 灰度: 浅, 黑, 深
line_style_map = containers.Map(...
    {'conservative', 'moderate', 'ambitious'}, ...
    {':', '-', '--'});
legend_label_map = containers.Map(...
    {'conservative', 'moderate', 'ambitious'}, ...
    {'保守方案', '中性方案', '激进方案'});

% --- 2. 绘图 ---
cS = results_all.moderate.cS; % 使用任一情景获取时间轴
T = cS.T_sim;
time_axis = cS.start_year : cS.time_Step : (cS.start_year + (T-1)*cS.time_Step);

fig_pps_gdp = figure('Name', 'PPS and PAYG Contribution Rates', 'Position', [100 428 425 172]);
hold on;

% -- 绘制恒定的PAYG比率虚线 --
yline(plot_data.payg_gdp_ratio * 100, ...
    'Color', [0.5 0.5 0.5], ...
    'LineStyle', ':', ...
    'LineWidth', 2, ...
    'DisplayName', 'PAYG缴费率');

% -- 循环绘制各PPS情景路径 --
for i = 1:length(scenarios_for_fig)
    scen_name = scenarios_for_fig{i};
    plot(time_axis, plot_data.(scen_name).pps_gdp_ratio * 100, ...
        'Color', color_map(scen_name), ...
        'LineStyle', line_style_map(scen_name), ...
        'LineWidth', 2.5, ...
        'DisplayName', legend_label_map(scen_name));
end

% --- 3. 添加图表格式 ---
xlabel('年份', 'FontName', 'SimSun');
ylabel('占GDP比重 (%)', 'FontName', 'SimSun');
legend('show', 'Location', 'best', 'FontName', 'SimSun');
grid on;
box on;
xlim([2023, 2100]);

% --- 4. 保存图像 ---
fprintf('   正在保存缴费率路径图至: %s\n', output_fig_pps);
try
    exportgraphics(fig_pps_gdp, output_fig_pps, 'Resolution', 300);
    fprintf('   ✅ 图像保存成功。\n');
catch ME
    fprintf('   ❌ 图像保存失败！错误信息: %s\n', ME.message);
end


%% --- 4. [生成附录表] 循环生成各情景福利分析(CV)表格 ---
fprintf('\n--- 4. [生成附录表] 循环生成各情景福利分析(CV)表格 ---\n');

% [!!! 核心修正: 在此逻辑块内明确定义所需变量 !!!]
assessment_years_welfare = [2023, 2030, 2040, 2050, 2070, 2100];
type_labels = {'高收入', '中低收入'};
sector_labels = {'A. 城镇职工', 'B. 城乡居民'};

data_base = results_all.(base_scenario_welfare);
cS_base = data_base.cS;
sigma = cS_base.sigma;
T_sim_base = cS_base.T_sim;

age_groups = {
    '青年 (25-29岁)', round((25 - cS_base.age1_orig)/cS_base.time_Step) + 1;
    '中年 (45-49岁)', round((45 - cS_base.age1_orig)/cS_base.time_Step) + 1;
    '临退 (60-64岁)', cS_base.aR_new_path(1);
    '老年 (70-74岁)', round((70 - cS_base.age1_orig)/cS_base.time_Step) + 1;
};

for i_reform_scen = 1:length(scenarios_for_welfare_comp)
    reform_scenario_name = scenarios_for_welfare_comp{i_reform_scen};
    data_reform = results_all.(reform_scenario_name);

    output_tex_table_welfare = fullfile(output_dir_tex, sprintf(appendix_file_pattern, reform_scenario_name));
    fprintf('   正在生成福利表: %s\n', output_tex_table_welfare);

    fileID_welfare = fopen(output_tex_table_welfare, 'w', 'n', 'UTF-8');
    fprintf(fileID_welfare, '%% =======================================================\n');
    fprintf(fileID_welfare, '%%  此文件由 tab_pps.m (v1.4) 自动生成\n');
    fprintf(fileID_welfare, '%% =======================================================\n');
    fprintf(fileID_welfare, '\\begin{table}[htbp]\n');
    fprintf(fileID_welfare, '\\centering\n');
    fprintf(fileID_welfare, '\\caption{引入个人养老金账户(PPS)的福利影响 (%s方案 vs. 无PPS情景, 补偿变化 CV, \\%%)}\n', reform_scenario_name);
    fprintf(fileID_welfare, '\\label{tab:welfare_impact_pps_%s}\n', reform_scenario_name);
    fprintf(fileID_welfare,'\\small\n');
    fprintf(fileID_welfare, '\\begin{threeparttable}\n');

    header_cols = ['l', repmat('c', 1, length(assessment_years_welfare))];
    header_line_1 = strjoin(cellfun(@(x) sprintf('%d', x), num2cell(assessment_years_welfare), 'UniformOutput', false), ' & ');
    fprintf(fileID_welfare, '\\begin{tabular}{%s}\n', header_cols);
    fprintf(fileID_welfare, '\\toprule\n');
    fprintf(fileID_welfare, '家庭类型与年龄组 & %s \\\\\n', header_line_1);
    fprintf(fileID_welfare, '\\midrule\n');

    for i_sector = 1:length(sector_labels)
        fprintf(fileID_welfare, '\\multicolumn{%d}{l}{\\textbf{%s}} \\\\\n', length(assessment_years_welfare) + 1, sector_labels{i_sector});
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
                    cv = NaN;
                    if t_idx > 0 && t_idx <= T_sim_base
                        dist_base_slice = data_base.final_Dist_path_h{h_idx}(:,:,:,a_idx,t_idx);
                        val_base_slice  = data_base.Val_path_h{h_idx}(:,:,:,a_idx,t_idx);
                        val_reform_slice= data_reform.Val_path_h{h_idx}(:,1,:,a_idx,t_idx);
                        mass_slice = sum(dist_base_slice, 'all');
                        if mass_slice > 1e-12
                            EV_base = sum(val_base_slice .* dist_base_slice, 'all') / mass_slice;
                            EV_reform = sum(val_reform_slice .* dist_base_slice, 'all') / mass_slice;
                            if EV_base < -1e-9 && EV_reform < -1e-9
                                cv = (EV_reform / EV_base)^(1 / (1 - sigma)) - 1;
                            end
                        end
                    end
                    fprintf(fileID_welfare, ' & %s', format_cv(cv));
                end
                fprintf(fileID_welfare, ' \\\\\n');
            end
        end
    end

    fprintf(fileID_welfare, '\\midrule\n');
    fprintf(fileID_welfare, '\\textbf{总体平均}');
    for i_yr = 1:length(assessment_years_welfare)
        year = assessment_years_welfare(i_yr);
        t_idx = round((year - cS_base.start_year) / cS_base.time_Step) + 1;
        cv_total = NaN;
        if t_idx > 0 && t_idx <= T_sim_base
            total_pop_t = sum(cS_base.Z_path_raw(:, t_idx));
            total_weighted_ev_base_all = 0;
            total_weighted_ev_reform_all = 0;
            for h_idx = 1:cS_base.nTypes
                for a_idx = 1:cS_base.aD_new
                    dist_base_slice = data_base.final_Dist_path_h{h_idx}(:,:,:,a_idx,t_idx);
                    val_base_slice  = data_base.Val_path_h{h_idx}(:,:,:,a_idx,t_idx);
                    val_reform_slice= data_reform.Val_path_h{h_idx}(:,1,:,a_idx,t_idx);
                    total_weighted_ev_base_all = total_weighted_ev_base_all + sum(val_base_slice .* dist_base_slice, 'all');
                    total_weighted_ev_reform_all = total_weighted_ev_reform_all + sum(val_reform_slice .* dist_base_slice, 'all');
                end
            end
            if total_pop_t > 1e-9
                EV_base_total_avg = total_weighted_ev_base_all / total_pop_t;
                EV_reform_total_avg = total_weighted_ev_reform_all / total_pop_t;
                if EV_base_total_avg < -1e-9 && EV_reform_total_avg < -1e-9
                     cv_total = (EV_reform_total_avg / EV_base_total_avg)^(1 / (1 - sigma)) - 1;
                end
            end
        end
        fprintf(fileID_welfare, ' & %s', format_cv(cv_total));
    end
    fprintf(fileID_welfare, ' \\\\\n');

    fprintf(fileID_welfare, '\\bottomrule\n');
    fprintf(fileID_welfare, '\\end{tabular}\n');
    fprintf(fileID_welfare, '\\begin{tablenotes}[para,flushleft]\n');
    fprintf(fileID_welfare, '  \\item 注：表中数值为补偿变化(CV)，衡量引入个人养老金账户相对于无个人养老金账户情景的福利差异。正值代表福利改善。\n');
    fprintf(fileID_welfare, '\\end{tablenotes}\n');
    fprintf(fileID_welfare, '\\end{threeparttable}\n');
    fprintf(fileID_welfare, '\\end{table}\n');
    fclose(fileID_welfare);
    fprintf('   ✅ 成功生成 LaTeX 表格文件: %s\n', output_tex_table_welfare);
end


%% --- 5. [生成表4.6] 计算并生成adj对比表格 ---
fprintf('\n--- 5. [生成表4.6] 计算并生成adj对比表格 ---\n');
assessment_years_adj = [2030, 2040, 2050, 2060, 2070, 2100];
scenarios_for_adj = {'None', 'conservative', 'moderate', 'ambitious'};
scenario_labels_adj = {'无PPS', '保守方案', '中性方案', '激进方案'};
fileID_adj = fopen(output_tex_table_adj, 'w', 'n', 'UTF-8');
fprintf(fileID_adj, '%% =======================================================\n');
fprintf(fileID_adj, '%%  此文件由 tab_pps.m (v1.4) 自动生成\n');
fprintf(fileID_adj, '%% =======================================================\n');
fprintf(fileID_adj, '\\begin{table}[htbp]\n');
fprintf(fileID_adj, '\\centering\n');
fprintf(fileID_adj, '\\caption{不同个人养老金(PPS)改革方案对退休人员综合福利替代水平的影响}\n');
fprintf(fileID_adj, '\\label{tab4.6:adj_welfare_pps}\n');
fprintf(fileID_adj, '\\small\n');
fprintf(fileID_adj, '\\begin{threeparttable}\n');

header_cols = 'l';
header_line_1 = {'年份'};
header_line_2_parts = {''};
cmidrule_parts = {};
current_col_idx = 1;
for i = 1:length(scenarios_for_adj)
    scen_label = scenario_labels_adj{i};
    if strcmp(scenarios_for_adj{i}, 'None')
        header_cols = [header_cols, 'c'];
        header_line_1{end+1} = sprintf('%s', scen_label);
        header_line_2_parts{end+1} = 'PAYG';
        cmidrule_parts{end+1} = sprintf('\\cmidrule(lr){%d-%d}', current_col_idx + 1, current_col_idx + 1);
        current_col_idx = current_col_idx + 1;
    else
        header_cols = [header_cols, 'rrr'];
        header_line_1{end+1} = sprintf('\\multicolumn{3}{c}{%s}', scen_label);
        header_line_2_parts = [header_line_2_parts, {'PAYG', '综合', '缴费率'}];
        cmidrule_parts{end+1} = sprintf('\\cmidrule(lr){%d-%d}', current_col_idx + 1, current_col_idx + 3);
        current_col_idx = current_col_idx + 3;
    end
end
final_header_1 = strjoin(header_line_1, ' & ');
final_header_2 = strjoin(header_line_2_parts, ' & ');
final_cmidrule = strjoin(cmidrule_parts, ' ');

fprintf(fileID_adj, '\\begin{tabular}{%s}\n', header_cols);
fprintf(fileID_adj, '\\toprule\n');
fprintf(fileID_adj, '%s \\\\\n', final_header_1);
fprintf(fileID_adj, '%s\n', final_cmidrule);
fprintf(fileID_adj, '%s \\\\\n', final_header_2);
fprintf(fileID_adj, '\\midrule\n');

for i_yr = 1:length(assessment_years_adj)
    year = assessment_years_adj(i_yr);
    fprintf(fileID_adj, '%d', year);
    for i_scen = 1:length(scenarios_for_adj)
        scen_name = scenarios_for_adj{i_scen};
        data = results_all.(scen_name);
        cS = data.cS;
        paramSF = data.results.ss_data.paramSF;
        t_idx = round((year - cS.start_year) / cS.time_Step) + 1;
        adj_payg = NaN;
        adj_comp = NaN;
        avg_pps_rate = NaN;
        if t_idx > 0 && t_idx <= cS.T_sim
            adj_payg = data.results.iter_results.adj_path_iter(t_idx);
            w_hat_for_b = data.ss0.w_hat;
            if t_idx > 1, w_hat_for_b = data.results.iter_results.w_hat_path_iter(t_idx-1); end
            aR_t = cS.aR_new_path(t_idx);
            b_hat_formula_h = SS.calculate_formula_benefits(w_hat_for_b, cS);
            total_formula_benefit_t = 0;
            for h = 1:cS.nTypes
                mass_retirees_h_t = sum(cS.Z_path_raw((aR_t + 1):end, t_idx)) * cS.type_weights(h);
                total_formula_benefit_t = total_formula_benefit_t + b_hat_formula_h(h) * mass_retirees_h_t;
            end
            if abs(total_formula_benefit_t) > 1e-9
                actual_payg_benefit_t = adj_payg * total_formula_benefit_t;
                pps_withdrawal_t = data.total_pps_withdrawal_hat_path(t_idx);
                adj_comp = (actual_payg_benefit_t + pps_withdrawal_t) / total_formula_benefit_t;
            end
            if ~strcmp(scen_name, 'None')
                w_hat_t = data.results.iter_results.w_hat_path_iter(t_idx);
                total_pps_contrib_t = 0;
                total_labor_income_workers_t = 0;
                for h = 1:cS.nTypes
                    cS_h = get_cs_for_type(cS, h);
                    dist_t_h = data.final_Dist_path_h{h}(:,:,:,:,t_idx);
                    pol_t_h = data.final_Pol_path_h{h}{t_idx};
                    for ia = 1:aR_t
                        mass_ia_h = dist_t_h(:,:,:,ia);
                        if sum(mass_ia_h, 'all') < 1e-30, continue; end
                        le_grid_slice = reshape(paramSF.leGridV, [1, 1, cS.nw_expanded]);
                        labor_income_grid = w_hat_t * cS_h.ageEffV_new(ia) .* le_grid_slice;
                        labor_income_grid_full = repmat(labor_income_grid, [cS.nk, cS.nkpps, 1]);
                        pps_rate_grid = pol_t_h(ia).pps_contrib_rate;
                        total_pps_contrib_t = total_pps_contrib_t + sum(pps_rate_grid .* labor_income_grid_full .* mass_ia_h, 'all');
                        total_labor_income_workers_t = total_labor_income_workers_t + sum(labor_income_grid_full .* mass_ia_h, 'all');
                    end
                end
                if abs(total_labor_income_workers_t) > 1e-9
                    avg_pps_rate = total_pps_contrib_t / total_labor_income_workers_t;
                else, avg_pps_rate = 0; end
            end
        end
        if strcmp(scen_name, 'None')
            fprintf(fileID_adj, ' & %.3f', adj_payg);
        else
            fprintf(fileID_adj, ' & %.3f & %.3f & %s', adj_payg, adj_comp, format_rate(avg_pps_rate));
        end
    end
    fprintf(fileID_adj, ' \\\\\n');
end

fprintf(fileID_adj, '\\bottomrule\n');
fprintf(fileID_adj, '\\end{tabular}\n');
fprintf(fileID_adj, '\\begin{tablenotes}[para,flushleft]\n');
fprintf(fileID_adj, '  \\item 注： “PAYG”列为在各情景下，为保证PAYG体系自身收支平衡所需的福利调整因子。“综合”列为一个模拟的综合福利兑现因子，其计算方式为 (实际PAYG福利 + PPS给付) / 公式化PAYG福利，衡量了引入PPS后，退休人员总养老金收入相对于原定PAYG福利水平的提升程度。“缴费率”为该年份所有工作年龄人口的个人养老金缴费总额占其劳动收入总额的加权平均比率。\n');
fprintf(fileID_adj, '\\end{tablenotes}\n');
fprintf(fileID_adj, '\\end{threeparttable}\n');
fprintf(fileID_adj, '\\end{table}\n');
fclose(fileID_adj);
fprintf('   ✅ 成功生成 LaTeX 表格文件: %s\n', output_tex_table_adj);

fprintf('\n--- 脚本执行完毕 ---\n');



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
    % =====================================================================
    % == 函数: get_cs_for_type (本地辅助函数)
    % == 版本: 1.1 - 修正PPS路径索引错误
    % == 目的: 为指定的家庭类型 h，准备一个专属的 cS 参数包。
    % =====================================================================
    cS_h = cS;
    cS_h.ageEffV_new = cS.ageEffV_new_h(:, h);
end

function str_out = format_cv(cv_value)
    % 辅助函数: 将CV数值格式化为带正负号和百分号的字符串
    if isnan(cv_value)
        str_out = '--';
        return;
    end

    if cv_value >= 0
        sign_char = '+';
    else
        sign_char = '';
    end
    str_out = sprintf('%s%.2f\\%%', sign_char, cv_value * 100);
end

function str_out = format_rate(rate_value)
    % 辅助函数: 将缴费率格式化为百分比字符串
    if isnan(rate_value)
        str_out = '--';
    else
        str_out = sprintf('%.2f\\%%', rate_value * 100);
    end
end