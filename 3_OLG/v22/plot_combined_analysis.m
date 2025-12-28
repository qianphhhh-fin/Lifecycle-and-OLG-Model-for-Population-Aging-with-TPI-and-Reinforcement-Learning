% =========================================================================
% == SCRIPT: plot_combined_analysis.m
% == 版本: [v1.7 - PAYG单位修正与精确定位版]
% ==
% == 目的:
% ==   1. [核心修正] 纠正了图4.2(PAYG流量)的单位错误，通过直接使用
% ==      模型计算出的原始GDP比率进行绘图，移除了不必要的单位转换。
% ==   2. [格式化] 严格按照要求设置两个图窗的'Position'属性，
% ==      以获得精确的尺寸和位置。
% ==   3. 确保所有文本标签均为简体中文。
% ==
% == 前置条件:
% ==   - TRANS/TPI_results_het_nopps_Retire_mild.mat
% =========================================================================

clear; close all;
addpath(pwd);

fprintf('=== 异质性模型综合分析与绘图脚本 (v1.7 - PAYG单位修正版) ===\n\n');

%% --- 1. 用户设定与文件定义 ---
fprintf('--- 1. 用户设定与文件定义 ---\n');

% --- 真实世界单位转换参数 ---
real_gdp_2023 = 126.06; % 2023年中国GDP (万亿元)
annual_inflation_rate = 0.02; % 假设的年化温和通货膨胀率 (2%)

% --- 图4.1 & 4.2 通用设定 ---
year_max = 2100; % 控制图表展示的结束年份

% --- 图4.2 专属设定 ---
exo_rho = 0.40; % 外生目标养老金替代率 (相对于当期有效工资)

% --- 输入文件 ---
file_tpi_results = 'TRANS/TPI_results_het_nopps.mat';

% --- 输出文件 ---
output_dir = 'tex/fig';
output_fig_4_1 = fullfile(output_dir, 'fig4.1_aging_impact_dashboard_het.png');
output_fig_4_2 = fullfile(output_dir, 'fig4.2_exo_pension_fund_het.png');

if ~exist(file_tpi_results, 'file'), error('找不到TPI结果文件: %s', file_tpi_results); end
if ~exist(output_dir, 'dir'), mkdir(output_dir); fprintf('   已创建输出目录: %s\n', output_dir); end

%% --- 2. 数据加载与解包 ---
fprintf('--- 2. 数据加载与解包 ---\n');
fprintf('   正在加载模型转轨结果: %s\n', file_tpi_results);

data_loaded = load(file_tpi_results);
results = data_loaded.results;
cS = data_loaded.cS;
ss0 = data_loaded.ss0;
final_Dist_path_h = data_loaded.final_Dist_path_h;
final_Pol_path_h = data_loaded.final_Pol_path_h;
paramSF = results.ss_data.paramSF;

[~, ~, cS] = population.generate_Z_path(cS, false);
fprintf('   已调用 population.generate_Z_path 更新人口数据。\n');

T = cS.T_sim;
nH = cS.nTypes;
Y_model_base = results.Y_path(1); 

fprintf('   ✅ 数据加载完毕。模型模拟 %d 期, 包含 %d 类家庭。\n', T, nH);
fprintf('   [单位转换基准] 2023年GDP=%.2f万亿元, 年通胀率=%.1f%%\n', real_gdp_2023, annual_inflation_rate*100);

%% --- 3. [图4.1] 老龄化冲击下的宏观经济动态仪表盘 ---
fprintf('\n--- 3. [生成图4.1] 老龄化宏观影响仪表盘 ---\n');

% --- 3.1 准备绘图数据 ---
fprintf('   正在为图4.1准备数据 (含真实单位转换)...\n');
paths_4_1 = struct();

time_axis_full = (cS.start_year : cS.time_Step : (cS.start_year + (T-1)*cS.time_Step));
time_mask = time_axis_full <= year_max;
time_axis = time_axis_full(time_mask);
total_pop_path = sum(cS.Z_path_raw, 1) * 1000;
pop_path_masked = total_pop_path(time_mask);

K_real_agg = convertToRealWorldUnits(results.K_p_path(time_mask), time_axis, cS, Y_model_base, real_gdp_2023, annual_inflation_rate, 'trillion');
Y_real_agg = convertToRealWorldUnits(results.Y_path(time_mask), time_axis, cS, Y_model_base, real_gdp_2023, annual_inflation_rate, 'trillion');
paths_4_1.K_pc_real = (K_real_agg ./ pop_path_masked) * 1e8;
paths_4_1.Y_pc_real = (Y_real_agg ./ pop_path_masked) * 1e8;

aggr_supply = aggregates.get_path_aggregates(final_Dist_path_h, final_Pol_path_h, cS, paramSF);
L_path_hat = aggr_supply.L_path;
total_labor_income_level = results.w_path .* L_path_hat;
per_capita_labor_income_level = total_labor_income_level ./ total_pop_path;
paths_4_1.w_pc_real = convertToRealWorldUnits(per_capita_labor_income_level(time_mask), time_axis, cS, Y_model_base, real_gdp_2023, annual_inflation_rate, 'ten_thousand');

paths_4_1.r_path = results.r_path(time_mask);
mass_retirees_path = zeros(1, T);
for t = 1:T
    aR_t = cS.aR_new_path(t);
    mass_retirees_path(t) = sum(cS.Z_path_raw((aR_t+1):end, t));
end
paths_4_1.retiree_share = mass_retirees_path(time_mask) ./ (total_pop_path(time_mask)/1000);
fprintf('   ✅ 数据准备完毕。\n');

% --- 3.2 绘图 ---
fprintf('   正在生成图表 (1x3布局)...\n');
fig1 = figure('Name', '老龄化冲击下的宏观经济动态(真实单位)', 'Position', [100 679 787 241]);
tcl1 = tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'normal');

style_K_pc = {'k-', 'LineWidth', 2.5}; 
style_Y_pc = {'--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2.5}; 
style_w_pc = {'k-', 'LineWidth', 2.5};
style_r    = {'k-', 'LineWidth', 2.5};
style_bg_line  = {'-.', 'Color', [0.8 0.8 0.8], 'LineWidth', 2};
color_left_axis = 'black';
color_right_axis = [0.5 0.5 0.5];
font_size_title = 14;
font_size_axis  = 11;
yy_lim_bar = [15 80];
tick_start_year = ceil(time_axis(1)/50)*50;
xticks_vec = tick_start_year :30 : time_axis(end);
xticks_vec = unique([time_axis(1), xticks_vec]);

ax1 = nexttile;
hold(ax1, 'on'); set(ax1, 'FontName', 'Times New Roman', 'FontSize', font_size_axis);
yyaxis(ax1, 'right');
p_bg_line1 = plot(ax1, time_axis, paths_4_1.retiree_share * 100, style_bg_line{:}, 'DisplayName', ['退休人口占比' newline '(右轴, %)']);
set(ax1, 'YColor', color_right_axis, 'YLim', yy_lim_bar);
yyaxis(ax1, 'left');
p_K_pc = plot(ax1, time_axis, paths_4_1.K_pc_real, style_K_pc{:}, 'DisplayName', '人均私人资本');
p_Y_pc = plot(ax1, time_axis, paths_4_1.Y_pc_real, style_Y_pc{:}, 'DisplayName', '人均产出');
set(ax1, 'YColor', color_left_axis); ylabel(ax1, '万元', 'FontName', 'SimSun'); ylim(ax1, 'padded');
title(ax1, '(a) 人均资本与产出', 'FontSize', font_size_title, 'FontName', 'SimSun');
legend(ax1, [p_K_pc, p_Y_pc, p_bg_line1], 'Location', 'northwest', 'FontSize', 10, 'Box', 'off', 'FontName', 'SimSun');
grid(ax1, 'on'); xlim(ax1, [time_axis(1), time_axis(end)]); 
set(ax1, 'XTick', 2030:10:2100); % 根据图像调整刻度
xtickangle(ax1, 40);
hold(ax1, 'off');

ax2 = nexttile;
hold(ax2, 'on'); set(ax2, 'FontName', 'Times New Roman', 'FontSize', font_size_axis);
yyaxis(ax2, 'right');
p_bg_line2 = plot(ax2, time_axis, paths_4_1.retiree_share * 100, style_bg_line{:});
set(ax2, 'YColor', color_right_axis, 'YLim', yy_lim_bar);
yyaxis(ax2, 'left');
plot(ax2, time_axis, paths_4_1.w_pc_real, style_w_pc{:});
set(ax2, 'YColor', color_left_axis); ylabel(ax2, '万元', 'FontName', 'SimSun'); ylim(ax2, 'padded');
title(ax2, '(b) 人均劳动收入', 'FontSize', font_size_title, 'FontName', 'SimSun');
legend(ax2, p_bg_line2, ['退休人口占比' newline '(右轴, %)'], 'Location', 'southeast', 'FontSize', 10, 'Box', 'off', 'FontName', 'SimSun');
grid(ax2, 'on'); xlim(ax2, [time_axis(1), time_axis(end)]); 
set(ax2, 'XTick', 2030:10:2100); % 根据图像调整刻度
xtickangle(ax2, 40);
hold(ax2, 'off');

ax3 = nexttile;
hold(ax3, 'on'); set(ax3, 'FontName', 'Times New Roman', 'FontSize', font_size_axis);
yyaxis(ax3, 'right');
p_bg_line3 = plot(ax3, time_axis, paths_4_1.retiree_share * 100, style_bg_line{:});
set(ax3, 'YColor', color_right_axis, 'YLim', yy_lim_bar);
yyaxis(ax3, 'left');
r_annual = (1 + paths_4_1.r_path).^(1/cS.time_Step) - 1;
plot(ax3, time_axis, r_annual * 100, style_r{:});
set(ax3, 'YColor', color_left_axis); ylim(ax3, 'padded');
title(ax3, '(c) 年化真实利率 (%)', 'FontSize', font_size_title, 'FontName', 'SimSun');
legend(ax3, p_bg_line3, ['退休人口占比' newline '(右轴, %)'], 'Location', 'northeast', 'FontSize', 10, 'Box', 'off', 'FontName', 'SimSun');
grid(ax3, 'on'); xlim(ax3, [time_axis(1), time_axis(end)]); 
set(ax3, 'XTick', 2030:10:2100); % 根据图像调整刻度
xtickangle(ax3, 40);
hold(ax3, 'off');

fprintf('   正在保存图像至: %s\n', output_fig_4_1);
try
    exportgraphics(fig1, output_fig_4_1, 'Resolution', 600);
    fprintf('   ✅ 图4.1保存成功。\n');
catch ME
    fprintf('   ❌ 图像保存失败！错误信息: %s\n', ME.message);
end

%% --- 4. [图4.2] PAYG养老金收支流量分析 ---
fprintf('\n--- 4. [生成图4.2] PAYG收支流量分析 ---\n');

% --- 4.1 会计核算 (不进行真实世界单位转换) ---
fprintf('   正在进行PAYG收支流量的会计核算 (原始模型单位)...\n');
iter_results = results.iter_results;
adj_path = iter_results.adj_path_iter;
w_hat_path_iter = iter_results.w_hat_path_iter;
w_hat_path_for_b = [ss0.w_hat, w_hat_path_iter(1:T-1)];
b_hat_formula_path_h = zeros(nH, T);
for t = 1:T
    b_hat_formula_path_h(:, t) = SS.calculate_formula_benefits(w_hat_path_for_b(t), cS);
end
b_hat_path_h_actual = adj_path .* b_hat_formula_path_h;

Y_path_level = results.Y_path;

mass_retirees_path_h = zeros(nH, T);
for t = 1:T
    aR_t = cS.aR_new_path(t);
    mass_retirees_total_t = sum(cS.Z_path_raw((aR_t+1):end, t));
    for h = 1:nH
        mass_retirees_path_h(h, t) = mass_retirees_total_t * cS.type_weights(h);
    end
end

b_hat_path_exo = exo_rho * w_hat_path_iter;

Total_PAYG_Contrib_level = sum(b_hat_path_h_actual .* cS.A_path .* mass_retirees_path_h, 1);
Total_Promised_Benefit_level = sum(repmat(b_hat_path_exo, nH, 1) .* cS.A_path .* mass_retirees_path_h, 1);

% [核心修正] 直接计算原始比率，不进行任何单位转换
contrib_gdp_ratio = Total_PAYG_Contrib_level ./ Y_path_level;
benefit_gdp_ratio = Total_Promised_Benefit_level ./ Y_path_level;
fprintf('   ✅ 会计核算完成。\n');

% --- 4.2 准备绘图数据 ---
fprintf('   正在为图4.2准备数据...\n');
paths_4_2 = struct();
paths_4_2.contrib_gdp = contrib_gdp_ratio(time_mask);
paths_4_2.benefit_gdp = benefit_gdp_ratio(time_mask);
fprintf('   ✅ 数据准备完毕。\n');

% --- 4.3 绘图 ---
fprintf('   正在生成图表...\n');
fig2 = figure('Name', 'PAYG收支流量分析', 'Position', [461 407 550 243]);

ax = gca;
hold(ax, 'on');
style_contrib = {'k-', 'LineWidth', 2};
style_benefit = {'--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2};
color_deficit_area = [0.8 0.8 0.8];
set(ax, 'FontName', 'Times New Roman', 'FontSize', 11);

deficit_gdp = paths_4_2.benefit_gdp - paths_4_2.contrib_gdp;
p_area = area(ax, time_axis, deficit_gdp * 100, ...
    'FaceColor', color_deficit_area, 'FaceAlpha', 0.6, ...
    'EdgeColor', 'none', 'DisplayName', 'PAYG缺口');
p_contrib = plot(ax, time_axis, paths_4_2.contrib_gdp * 100, style_contrib{:}, 'DisplayName', 'PAYG缴费收入');
p_benefit = plot(ax, time_axis, paths_4_2.benefit_gdp * 100, style_benefit{:}, 'DisplayName', '承诺的福利支出');

ylim(ax, [0, max(paths_4_2.benefit_gdp)*100*1.1]);
ylabel(ax, '占GDP的百分比 (%)', 'FontName', 'SimSun');
xlabel(ax, '年份', 'FontName', 'SimSun');

title_obj = title(ax, 'PAYG收支流量占GDP比重', 'FontSize', 14, 'FontName', 'SimSun');
legend_obj = legend(ax, [p_contrib, p_benefit, p_area], 'Location', 'northwest', 'FontSize', 10, 'Box', 'off', 'FontName', 'SimSun');

grid(ax, 'on'); box(ax, 'on');
xlim(ax, [time_axis(1), time_axis(end)]);
set(ax, 'XTick', 2030:10:2100); % 根据图像调整刻度
xtickangle(ax, 40);
hold(ax, 'off');

fprintf('   ✅ 图表生成完毕。\n');

% --- 4.4 保存图像 ---
fprintf('   正在保存图像至: %s\n', output_fig_4_2);
try
    exportgraphics(fig2, output_fig_4_2, 'Resolution', 600);
    fprintf('   ✅ 图4.2保存成功。\n');
catch ME
    fprintf('   ❌ 图像保存失败！错误信息: %s\n', ME.message);
end

fprintf('\n--- 脚本执行完毕 ---\n');

%% --- [本地辅助函数] ---
function real_world_path = convertToRealWorldUnits(model_path, time_axis, cS, Y_model_base, real_gdp_base, inflation_rate, target_unit)
    scale_factor = real_gdp_base / Y_model_base;
    real_path = model_path * scale_factor;
    years_since_start = time_axis - cS.start_year;
    price_level_path = (1 + inflation_rate) .^ years_since_start;
    nominal_path = real_path .* price_level_path;
    
    switch lower(target_unit)
        case 'trillion' 
            real_world_path = nominal_path;
        case 'ten_thousand' 
            real_world_path = nominal_path * 1e8;
        otherwise
            error('未知的目标单位: %s。请使用 "trillion" 或 "ten_thousand"。', target_unit);
    end
end