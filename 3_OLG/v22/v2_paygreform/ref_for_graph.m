% =========================================================================
% == SCRIPT: analyze_retirement_reforms.m
% ==
% == 版本: [v2.0 - 分析视角与可视化重构]
% ==
% == 目的:
% ==   1. [分析视角] 将福利分析的基准情景从'moderate'改为'mild'。
% ==   2. [年龄组调整] 福利分析聚焦于45-59岁的三个中年及临近退休年龄组。
% ==   3. [可视化重构] 图1改为直接展示退休年龄路径和PAYG收支缺口。
% ==   4. [风格统一] 所有图表统一为灰度风格。
% =========================================================================

clear; close all; clc;
addpath(pwd);
fprintf('=== 延迟退休改革多情景对比分析脚本 (v2.0 - 视角重构版) ===\n\n');

%% --- 1. 用户设定与文件定义 ---
fprintf('--- 1. 用户设定与文件定义 ---\n');

% --- 分析设定 ---
scenarios = {'mild','moderate', 'aggressive'};
% [核心修改] 将福利分析的基准改为'mild'
base_scenario_for_welfare = 'mild'; 
initial_fund_to_gdp_ratio = 0.0239;
year_max = 2100;

% --- 文件路径模板 (仅需TPI结果文件) ---
file_tpi_pattern = 'TRANS/TPI_results_het_nopps_Retire_%s.mat';

% --- 输出路径 ---
output_dir = 'tex/fig/';
if ~exist(output_dir, 'dir'), mkdir(output_dir); end
output_fig_policy_deficit = fullfile(output_dir, 'fig4.5_retire_reform_policy_deficit.png');
output_fig_welfare = fullfile(output_dir, 'fig4.6_retire_reform_welfare_impact.png');

fprintf('   基准方案 (福利): %s\n', base_scenario_for_welfare);
fprintf('   初始基金/GDP比例: %.2f%%\n', initial_fund_to_gdp_ratio * 100);

%% --- 2. [宏观分析] 加载数据并计算核心宏观路径 ---
fprintf('\n--- 2. [宏观分析] 循环加载各方案，计算核心宏观路径 ---\n');

results_all = struct(); % 用于存储所有方案的核心计算结果

for i = 1:length(scenarios)
    scenario_name = scenarios{i};
    fprintf('\n   正在处理方案: [%s] ...\n', upper(scenario_name));

    % --- a. 加载数据 ---
    file_tpi = sprintf(file_tpi_pattern, scenario_name);
    if ~exist(file_tpi, 'file'), error('找不到文件: %s', file_tpi); end

    data_tpi = load(file_tpi);
    cS = data_tpi.cS;
    
    % 从加载的数据中获取路径和参数
    Y_path_level = data_tpi.results.Y_path;
    w_path_level = data_tpi.results.w_path;
    dist_path_h = data_tpi.final_Dist_path_h;
    ss_data = data_tpi.results.ss_data;
    paramSF = ss_data.paramSF;

    % [核心修改] 存储绘图和后续分析所需的数据
    results_all.(scenario_name).data_tpi = data_tpi;
    results_all.(scenario_name).Y_path_level = Y_path_level;
    results_all.(scenario_name).aR_new_path = cS.aR_new_path; % 存储退休年龄索引路径

    % --- b. 重构收支路径以计算缺口 ---
    w_hat_path = w_path_level ./ cS.A_path; % 需要 hat 量

    [Total_Contributions_Level, Total_Formula_Benefits_Level] = ...
        recalculate_payg_flows_from_tpi_logic(w_hat_path, dist_path_h, cS, paramSF, ss_data.ss0);

    % === 计算PAYG年度缺口占GDP的比重 ===
    Deficit_Level = Total_Formula_Benefits_Level - Total_Contributions_Level;
    Deficit_to_GDP_path = Deficit_Level ./ Y_path_level;
    results_all.(scenario_name).Deficit_to_GDP_path = Deficit_to_GDP_path;
    
    fprintf('   ✅ 核心宏观路径处理完成。\n');
end

%% --- 3. [微观分析] 计算比较福利 (CV) ---
fprintf('\n--- 3. [微观分析] 计算各方案相对基准的福利变化 ---\n');

% --- a. 定义劳动负效用参数 ---
add_labor_disutility = false; 
labor_disutil_param_psi = 1.5; 

if add_labor_disutility
    fprintf('   *** 将在福利评估中加入劳动负效用 (ψ = %.2f) ***\n', labor_disutil_param_psi);
else
    fprintf('   *** 将在福利评估中【不】加入劳动负效用 ***\n');
end


% --- b. 重构基准方案 ('mild') 的价值和策略路径 ---
fprintf('   正在重构基准方案 [%s] 的路径...\n', upper(base_scenario_for_welfare));
base_results = results_all.(base_scenario_for_welfare);
cS_base = base_results.data_tpi.cS;
iter_results_base = base_results.data_tpi.results.iter_results;
ss_data_base = base_results.data_tpi.results.ss_data;
paramSF_base = ss_data_base.paramSF;
valF_base_h = ss_data_base.valF_h;
polF_base_h = ss_data_base.polF_h;

Val_path_base_h = cell(cS_base.nTypes, 1);
Pol_path_base_h = cell(cS_base.nTypes, 1);
for h = 1:cS_base.nTypes
    pathS_base_h = get_pathS_for_backward(base_results.data_tpi.results, iter_results_base, h, cS_base, ss_data_base.ss0);
    cS_h_base = get_cs_for_type(cS_base, h);
    [Pol_path_base_h{h}, Val_path_base_h{h}] = household.backward_hh(pathS_base_h, cS_h_base, paramSF_base, valF_base_h{h}, polF_base_h{h});
end

if add_labor_disutility
    Val_path_base_for_CV_h = cell(cS_base.nTypes, 1);
    for h = 1:cS_base.nTypes
        cS_h_base = get_cs_for_type(cS_base, h);
        Val_path_base_for_CV_h{h} = recalculate_vpath_with_labor_disutility(...
            Val_path_base_h{h}, Pol_path_base_h{h}, cS_h_base, paramSF_base, labor_disutil_param_psi);
    end
else
    Val_path_base_for_CV_h = Val_path_base_h;
end
fprintf('   ✅ 基准方案路径重构完成。\n');


% --- c. 循环计算各改革方案相对基准的CV ---
reform_scenarios_for_welfare = setdiff(scenarios, base_scenario_for_welfare);
cv_results = struct();

for i = 1:length(reform_scenarios_for_welfare)
    reform_name = reform_scenarios_for_welfare{i};
    fprintf('   正在计算方案 [%s] 相对于 [%s] 的福利...\n', upper(reform_name), upper(base_scenario_for_welfare));
    
    reform_results = results_all.(reform_name);
    cS_reform = reform_results.data_tpi.cS;
    iter_results_reform = reform_results.data_tpi.results.iter_results;
    ss_data_reform = reform_results.data_tpi.results.ss_data;
    paramSF_reform = ss_data_reform.paramSF;
    valF_reform_h = ss_data_reform.valF_h;
    polF_reform_h = ss_data_reform.polF_h;

    Val_path_reform_h = cell(cS_reform.nTypes, 1);
    Pol_path_reform_h = cell(cS_reform.nTypes, 1);
    for h = 1:cS_reform.nTypes
        pathS_reform_h = get_pathS_for_backward(reform_results.data_tpi.results, iter_results_reform, h, cS_reform, ss_data_reform.ss0);
        cS_h_reform = get_cs_for_type(cS_reform, h);
        [Pol_path_reform_h{h}, Val_path_reform_h{h}] = household.backward_hh(pathS_reform_h, cS_h_reform, paramSF_reform, valF_reform_h{h}, polF_reform_h{h});
    end

    if add_labor_disutility
        Val_path_reform_for_CV_h = cell(cS_reform.nTypes, 1);
        for h = 1:cS_reform.nTypes
            cS_h_reform = get_cs_for_type(cS_reform, h);
            Val_path_reform_for_CV_h{h} = recalculate_vpath_with_labor_disutility(...
                Val_path_reform_h{h}, Pol_path_reform_h{h}, cS_h_reform, paramSF_reform, labor_disutil_param_psi);
        end
    else
        Val_path_reform_for_CV_h = Val_path_reform_h;
    end
    
    % [核心修改] 调整福利分析的年龄组
    age_groups_for_figure = {
        '中年I (45-49岁)', round((45 - cS_base.age1_orig)/cS_base.time_Step) + 1;
        '中年II (50-54岁)', round((50 - cS_base.age1_orig)/cS_base.time_Step) + 1;
        '临近退休 (55-59岁)', round((55 - cS_base.age1_orig)/cS_base.time_Step) + 1
    };
    num_age_groups_fig = size(age_groups_for_figure, 1);
    
    T_max_comp = min(cS_base.T_sim, cS_reform.T_sim);
    
    cv_paths_all_types = zeros(num_age_groups_fig, T_max_comp);

    for i_age_grp = 1:num_age_groups_fig
        a_idx = age_groups_for_figure{i_age_grp, 2};
        if a_idx > cS_base.aD_new, continue; end

        for t = 1:T_max_comp
            total_mass_t = 0;
            weighted_EV_base_t = 0;
            weighted_EV_reform_t = 0;
            
            for h = 1:cS_base.nTypes
                dist_base_slice = base_results.data_tpi.final_Dist_path_h{h}(:, :, :, a_idx, t);
                mass_slice = sum(dist_base_slice, 'all');

                if mass_slice > 1e-12
                   total_mass_t = total_mass_t + mass_slice;
                   val_base_slice = Val_path_base_for_CV_h{h}(:, :, :, a_idx, t);
                   val_reform_slice = Val_path_reform_for_CV_h{h}(:, :, :, a_idx, t);
                   weighted_EV_base_t = weighted_EV_base_t + sum(val_base_slice .* dist_base_slice, 'all');
                   weighted_EV_reform_t = weighted_EV_reform_t + sum(val_reform_slice .* dist_base_slice, 'all');
                end
            end
            
            if total_mass_t < 1e-9
                cv_paths_all_types(i_age_grp, t) = 0;
                continue;
            end

            EV_base_avg = weighted_EV_base_t / total_mass_t;
            EV_reform_avg = weighted_EV_reform_t / total_mass_t;
            
            if EV_base_avg > -1e-9 || EV_reform_avg > -1e-9
                cv_paths_all_types(i_age_grp, t) = NaN;
                continue;
            end

            cv = (EV_reform_avg / EV_base_avg)^(1 / (1 - cS_base.sigma)) - 1;
            
            if ~isreal(cv)
                cv_paths_all_types(i_age_grp, t) = NaN;
            else
                cv_paths_all_types(i_age_grp, t) = cv;
            end
        end
    end
    cv_results.(reform_name) = cv_paths_all_types;
    fprintf('   ✅ 福利计算完成。\n');
end
fprintf('--- 福利分析全部完成 ---\n');


%% --- 4. 绘图 ---
fprintf('\n--- 4. 生成结果图表 ---\n');

% --- 设定灰度颜色和线型 ---
color_mild = [0.65 0.65 0.65]; 
color_moderate = [0.35 0.35 0.35]; 
color_aggressive = [0.0 0.0 0.0];
colors = containers.Map(scenarios, {color_mild, color_moderate, color_aggressive});
line_styles = {'--', '-', ':'}; 
line_styles_map = containers.Map(scenarios, line_styles);

% --- 图1: 政策路径、抚养比与PAYG收支缺口 ---
fig1 = figure('Name', '政策路径、抚养比与PAYG缺口', 'Position', [100 100 950 280]); 
tcl1 = tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact'); 

cS_plot = results_all.(scenarios{1}).data_tpi.cS;
time_axis_model = (cS_plot.start_year : cS_plot.time_Step : (cS_plot.start_year + (cS_plot.T_sim-1)*cS_plot.time_Step));
time_axis_annual = (cS_plot.start_year : 1 : min(year_max, time_axis_model(end)));
% -- Panel (a): 平均法定退休年龄路径 (精确复现年度路径) --
ax1a = nexttile;
hold(ax1a, 'on');

legend_policy = struct();
legend_policy.mild = '平缓(我国现行)';
legend_policy.moderate = '温和';
legend_policy.aggressive = '激进';

annual_years_vector = cS_plot.start_year : min(year_max, cS_plot.end_year);
for i = 1:length(scenarios)
    s = scenarios{i};
    % --- 根据方案名称设定改革参数 (与utils.retire_age_path保持一致) ---
    if strcmpi(s, 'mild')
        reform_start_year = 2025; initial_retire_age_male = 60; initial_retire_age_female = 55;
        target_retire_age = 67; annual_increment_years = 0.25;
    elseif strcmpi(s, 'moderate')
        reform_start_year = 2025; initial_retire_age_male = 60; initial_retire_age_female = 55;
        target_retire_age = 67; annual_increment_years = 1/3;
    elseif strcmpi(s, 'aggressive')
        reform_start_year = 2025; initial_retire_age_male = 60; initial_retire_age_female = 55;
        target_retire_age = 67; annual_increment_years = 0.5;
    end
    
    % [核心修改] 创建新的图例标签
    months_per_year = annual_increment_years * 12;
    legend_label = sprintf('%s(%.1d个月/年)', getfield(legend_policy,s), months_per_year);

    % --- 生成男女各自的年度路径 ---
    retire_age_male_annual_path = ones(size(annual_years_vector)) * initial_retire_age_male;
    retire_age_female_annual_path = ones(size(annual_years_vector)) * initial_retire_age_female;
    for idx = 2:length(annual_years_vector)
        year = annual_years_vector(idx);
        if year >= reform_start_year
            new_age_male = retire_age_male_annual_path(idx-1) + annual_increment_years;
            retire_age_male_annual_path(idx) = min(target_retire_age, new_age_male);
            new_age_female = retire_age_female_annual_path(idx-1) + annual_increment_years;
            retire_age_female_annual_path(idx) = min(target_retire_age, new_age_female);
        end
    end
    % --- 计算加权平均并绘图 ---
    avg_retire_age_annual_path = 0.5 * retire_age_male_annual_path + 0.5 * retire_age_female_annual_path;
    % [核心修改] 使用新的图例标签
    plot(ax1a, annual_years_vector, avg_retire_age_annual_path, 'LineWidth', 2.5, 'Color', colors(s), 'LineStyle', line_styles_map(s), 'DisplayName', legend_label);
end
hold(ax1a, 'off');
grid(ax1a, 'on'); box(ax1a, 'on');
xlim(ax1a, [annual_years_vector(1), annual_years_vector(end)]);
title(ax1a, '(a) 平均法定退休年龄', 'FontName', 'SimSun', 'FontSize', 12);
ylabel(ax1a, '物理年龄 (岁)', 'FontName', 'SimSun');
xlabel(ax1a, '年份', 'FontName', 'SimSun');
legend(ax1a, 'show', 'Location', 'southeast', 'FontName', 'SimSun');

% -- Panel (b): 老年抚养比路径 --
ax1b = nexttile;
hold(ax1b, 'on');
for i = 1:length(scenarios)
    s = scenarios{i};
    dependency_ratio_path = zeros(1, cS_plot.T_sim);
    for t = 1:cS_plot.T_sim
        aR_t = results_all.(s).aR_new_path(t); % 获取当期退休年龄索引
        working_pop_t = sum(cS_plot.Z_path_raw(1:aR_t, t));
        retired_pop_t = sum(cS_plot.Z_path_raw(aR_t+1:end, t));
        if working_pop_t > 1e-9
            dependency_ratio_path(t) = (retired_pop_t / working_pop_t) * 100;
        else
            dependency_ratio_path(t) = NaN;
        end
    end
    plot(ax1b, time_axis_model, dependency_ratio_path, 'LineWidth', 2.5, 'Color', colors(s), 'LineStyle', line_styles_map(s));
end
hold(ax1b, 'off');
grid(ax1b, 'on'); box(ax1b, 'on');
xlim(ax1b, [time_axis_annual(1), time_axis_annual(end)]);
title(ax1b, '(b) 老年抚养比', 'FontName', 'SimSun', 'FontSize', 12);
ylabel(ax1b, '退休人口 / 劳动人口 (%)', 'FontName', 'SimSun');
xlabel(ax1b, '年份', 'FontName', 'SimSun');

% -- Panel (c): PAYG年度收支缺口/GDP (年度插值) --
ax1c = nexttile;
hold(ax1c, 'on');

for i = 1:length(scenarios)
    s = scenarios{i};
    deficit_to_gdp_annual = interp1(time_axis_model, results_all.(s).Deficit_to_GDP_path, time_axis_annual, 'pchip');
    plot(ax1c, time_axis_annual, deficit_to_gdp_annual*100, 'LineWidth', 2.5, 'Color', colors(s), 'LineStyle', line_styles_map(s));
end
yline(ax1c, 0, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
hold(ax1c, 'off');
grid(ax1c, 'on'); box(ax1c, 'on');
xlim(ax1c, [time_axis_annual(1), time_axis_annual(end)]);
title(ax1c, '(c) 现收现付制年度收支缺口', 'FontName', 'SimSun', 'FontSize', 12);
ylabel(ax1c, '缺口占GDP比重 (%)', 'FontName', 'SimSun');
xlabel(ax1c, '年份', 'FontName', 'SimSun');

fprintf('   正在保存政策与缺口图...\n');
exportgraphics(fig1, output_fig_policy_deficit, 'Resolution', 300);


% --- 图2: 福利影响分析 ---
fig2 = figure('Name', '福利影响分析', 'Position', [461 646 637 276]);
tcl2 = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

cS_base = results_all.(base_scenario_for_welfare).data_tpi.cS;
T_max_comp_plot = size(cv_results.(reform_scenarios_for_welfare{1}), 2);
time_axis_welfare = cS_base.start_year : cS_base.time_Step : (cS_base.start_year + (T_max_comp_plot-1)*cS_base.time_Step);
time_mask_welfare = time_axis_welfare <= year_max;
time_axis_plot_welfare = time_axis_welfare(time_mask_welfare);

% --- 灰度绘图风格 ---
style_mid1 = {'-o', 'Color', [0.0 0.0 0.0], 'LineWidth', 2, 'MarkerSize', 4, 'MarkerFaceColor', [0.2 0.2 0.2]};
style_mid2 = {'--s', 'Color', [0.4 0.4 0.4], 'LineWidth', 2, 'MarkerSize', 4};
style_preretire = {':d', 'Color', [0.7 0.7 0.7], 'LineWidth', 2, 'MarkerSize', 4};
styles_age = {style_mid1, style_mid2, style_preretire};

% -- Panel (a): Moderate vs. Mild --
ax2a = nexttile;
hold(ax2a, 'on');
cv_data_mod = cv_results.moderate;
for i_age_grp = 1:num_age_groups_fig
    plot_data = cv_data_mod(i_age_grp, time_mask_welfare) * 100;
    plot(ax2a, time_axis_plot_welfare, plot_data, styles_age{i_age_grp}{:}, 'DisplayName', age_groups_for_figure{i_age_grp, 1});
end
yline(ax2a, 0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
hold(ax2a, 'off');
grid(ax2a, 'on'); box(ax2a, 'on');
xlim(ax2a, [time_axis_plot_welfare(1), time_axis_plot_welfare(end)]);
title(ax2a, '(a) 温和方案 vs. 平缓方案', 'FontName', 'SimSun', 'FontSize', 12);
ylabel(ax2a, ['相对于平缓方案的福利变化 (%)'], 'FontName', 'SimSun');
xlabel(ax2a, '年份', 'FontName', 'SimSun');
legend(ax2a, 'show', 'Location', 'best', 'FontName', 'SimSun');
ylim([-8,8])


% -- Panel (b): Aggressive vs. Mild --
ax2b = nexttile;
hold(ax2b, 'on');
cv_data_agg = cv_results.aggressive;
for i_age_grp = 1:num_age_groups_fig
    plot_data = cv_data_agg(i_age_grp, time_mask_welfare) * 100;
    plot(ax2b, time_axis_plot_welfare, plot_data, styles_age{i_age_grp}{:}, 'DisplayName', age_groups_for_figure{i_age_grp, 1});
end
yline(ax2b, 0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
hold(ax2b, 'off');
grid(ax2b, 'on'); box(ax2b, 'on');
xlim(ax2b, [time_axis_plot_welfare(1), time_axis_plot_welfare(end)]);
title(ax2b, '(b) 激进方案 vs. 平缓方案', 'FontName', 'SimSun', 'FontSize', 12);
xlabel(ax2b, '年份', 'FontName', 'SimSun');
ylim([-8,8])

fprintf('   正在保存福利影响图...\n');
exportgraphics(fig2, output_fig_welfare, 'Resolution', 300);

fprintf('\n--- 脚本执行完毕 ---\n');


%% --- [本地辅助函数] ---

function pathS_h = get_pathS_for_backward(results, iter_results, h, cS, ss0)
    % =====================================================================
    % == 函数: get_pathS_for_backward (v2.0 - 最终稳健版)
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
% == 目的: 为指定的家庭类型 h，准备一个专属的 cS 参数包。
% =====================================================================
cS_h = cS;
cS_h.ageEffV_new = cS.ageEffV_new_h(:, h);
end

function [Total_Contributions_Level, Total_Formula_Benefits_Level] = recalculate_payg_flows_from_tpi_logic(w_hat_path, Dist_path_h, cS, paramSF, ss0)
    % =========================================================================
    % == 函数: recalculate_payg_flows_from_tpi_logic
    % == 目的: 精确复现 TPI 中的养老金收支计算逻辑。
    % =========================================================================
    T = cS.T_sim;
    nH = cS.nTypes;
    total_pension_pot_hat_path = zeros(1, T);
    total_formula_benefits_hat_path = zeros(1, T);
    w_hat_path_for_formula = [ss0.w_hat, w_hat_path(1:T-1)];
    b_hat_formula_path_h = zeros(nH, T);
    for t = 1:T
        b_hat_formula_path_h(:, t) = SS.calculate_formula_benefits(w_hat_path_for_formula(t), cS);
    end
    for t = 1:T
        total_contrib_t = 0;
        total_formula_benefit_t = 0;
        aR_t = cS.aR_new_path(t);
        for h = 1:nH
            L_path_supply_h_t = get_labor_supply_for_type_local(h, {Dist_path_h{h}(:,:,:,:,t)}, cS, paramSF, aR_t);
            total_contrib_t = total_contrib_t + cS.theta_path_h(h,t) * w_hat_path(t) * L_path_supply_h_t;
            mass_retirees_path_h_t = sum(cS.Z_path_raw((aR_t + 1):end, t)) * cS.type_weights(h);
            total_formula_benefit_t = total_formula_benefit_t + b_hat_formula_path_h(h, t) * mass_retirees_path_h_t;
        end
        total_pension_pot_hat_path(t) = total_contrib_t;
        total_formula_benefits_hat_path(t) = total_formula_benefit_t;
    end
    Total_Contributions_Level = total_pension_pot_hat_path .* cS.A_path;
    Total_Formula_Benefits_Level = total_formula_benefits_hat_path .* cS.A_path;
end

function L_path_h = get_labor_supply_for_type_local(h, Dist_path_h, cS, paramS, aR_t_period)
% =========================================================================
% == 函数: get_labor_supply_for_type_local
% == 目的: 精确复现 TPI.get_labor_supply_path_for_type 的逻辑。
% =========================================================================
T = size(Dist_path_h{1}, 5);
L_path_h = zeros(1, T);
cS_h = cS;
cS_h.ageEffV_new = cS.ageEffV_new_h(:,h);
dist_path_for_type_h = Dist_path_h{1};
for t = 1:T
    L_t_h = 0;
    dist_t_abs_h = dist_path_for_type_h(:,:,:,:,t);
    if length(aR_t_period) > 1, aR_now = aR_t_period(t); else, aR_now = aR_t_period; end
    for ia = 1:aR_now
        mass_ia_abs_h = dist_t_abs_h(:,:,:,ia);
        if sum(mass_ia_abs_h(:)) < 1e-30, continue; end
        le_grid_slice = reshape(paramS.leGridV, [1, 1, cS.nw_expanded]);
        labor_supply_slice = cS_h.ageEffV_new(ia) .* le_grid_slice;
        L_t_h = L_t_h + sum(labor_supply_slice .* mass_ia_abs_h, 'all');
    end
    L_path_h(t) = L_t_h;
end
end

function Val_path_welfare = recalculate_vpath_with_labor_disutility(Val_path, Pol_path, cS, paramS, labor_disutil_param)
    % =====================================================================
    % == 函数: recalculate_vpath_with_labor_disutility
    % == 目的: 重新计算包含劳动负效用的“福利衡量用”价值函数路径。
    % =====================================================================
    T = cS.T_sim;
    Val_path_welfare = zeros(size(Val_path));
    Val_path_welfare(:,:,:,:,T) = Val_path(:,:,:,:,T);
    A_path_ext = [cS.A_path, cS.A_path(end) * (1 + ((1+cS.g_A_ss)^cS.time_Step-1))];
    g_A_path = A_path_ext(2:end)./A_path_ext(1:end-1)-1;

    for t = (T-1) : -1 : 1
        valS_t_welfare = -Inf(cS.nk, cS.nkpps, cS.nw_expanded, cS.aD_new);
        Vprime_tp1_welfare = Val_path_welfare(:,:,:,:,t+1);
        g_A_period = g_A_path(t);
        growth_factor_bgp = (1 + g_A_period);
        discount_factor_V_prime = cS.beta * growth_factor_bgp^(1 - cS.sigma);
        pol_t_struct_array = Pol_path{t}; 

        for a_idx = cS.aD_new : -1 : 1
            pol_age_t = pol_t_struct_array(a_idx);
            c_mat = pol_age_t.c;
            if isempty(c_mat), continue; end
            util_c_mat = (c_mat.^(1-cS.sigma))./(1-cS.sigma);
            labor_supply_mat = pol_age_t.l;
            nu = 1.0;
            psi = labor_disutil_param;
            labor_disutil_mat = psi * (labor_supply_mat.^(1+nu))./(1+nu);
            ev_on_policy_welfare = zeros(cS.nk, cS.nkpps, cS.nw_expanded);
            if a_idx < cS.aD_new
                vPrime_next_age_welfare = Vprime_tp1_welfare(:,:,:,a_idx+1);
                k_prime_mat = pol_age_t.k_prime;
                kpps_prime_mat = pol_age_t.kpps_prime;
                trans_mats = paramS.TrProbM_by_age;
                vPrime_interpolants_welfare = cell(cS.nw_expanded, 1);
                for ie_next = 1:cS.nw_expanded
                    if cS.nkpps > 1 && cS.pps_active
                         vPrime_interpolants_welfare{ie_next} = griddedInterpolant({cS.kGridV, cS.kppsGridV}, vPrime_next_age_welfare(:, :, ie_next), 'pchip', 'pchip');
                    else
                         vPrime_interpolants_welfare{ie_next} = griddedInterpolant(cS.kGridV, squeeze(vPrime_next_age_welfare(:, :, ie_next)), 'pchip', 'pchip');
                    end
                end
                for ie = 1:cS.nw_expanded
                    trans_prob_row = trans_mats{a_idx}(ie, :);
                    ev_slice = zeros(cS.nk, cS.nkpps);
                    for ie_next = 1:cS.nw_expanded
                        if trans_prob_row(ie_next) > 1e-9
                           if cS.nkpps > 1 && cS.pps_active
                               v_interp = vPrime_interpolants_welfare{ie_next}(k_prime_mat(:,:,ie), kpps_prime_mat(:,:,ie));
                           else
                               v_interp = vPrime_interpolants_welfare{ie_next}(k_prime_mat(:,1,ie));
                           end
                           ev_slice = ev_slice + trans_prob_row(ie_next) * v_interp;
                        end
                    end
                    ev_on_policy_welfare(:,:,ie) = ev_slice;
                end
            end
            survival_rate = cS.s_pathV(a_idx, t);
            k_prime_for_bequest = pol_age_t.k_prime;
            if cS.pps_active, k_prime_for_bequest = k_prime_for_bequest + pol_age_t.kpps_prime; end
            util_bequest_mat = utils.bequest_utility(k_prime_for_bequest, cS);
            Future_V_welfare_discounted = discount_factor_V_prime * (survival_rate * ev_on_policy_welfare + (1 - survival_rate) * util_bequest_mat);
            valS_t_welfare(:,:,:,a_idx) = util_c_mat - labor_disutil_mat + Future_V_welfare_discounted;
        end
        Val_path_welfare(:,:,:,:,t) = valS_t_welfare;
    end
end