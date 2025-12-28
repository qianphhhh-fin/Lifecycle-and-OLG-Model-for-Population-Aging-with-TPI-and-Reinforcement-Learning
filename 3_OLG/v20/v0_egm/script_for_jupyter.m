% =========================================================================
% == SCRIPT: exo_pension.m
% ==
% == 版本: [v2.1 - 1x2灰度精修版]
% ==
% == 目的:
% ==   - 采用1x2的横向布局，生成一幅紧凑、专业的灰度图。
% ==   - (a) 左图展示基金存量动态。
% ==   - (b) 右图展示PAYG收支流量动态。
% ==   - 严格遵循特定的尺寸、字体、刻度和颜色要求。
% ==
% == 前置条件:
% ==   - 'TRANS/TPI_results_main_nopps.mat'
% ==   - 'TRANS/iter_results_main_nopps.mat'
% =========================================================================

clear; close all;
addpath(pwd);

fprintf('--- 启动外生养老金基金模拟与绘图脚本 (v2.1) ---\n');

%% --- 1. 用户设定 ---
exo_rho = 0.40;
initial_fund_to_gdp_ratio = 0.055;
year_max = 2100;

%% --- 2. 定义文件路径与加载数据 ---
file_tpi_results = 'TRANS/TPI_results_main_nopps.mat';
file_iter_results = 'TRANS/iter_results_main_nopps.mat';
output_dir = 'tex/fig';
output_filename = fullfile(output_dir, 'fig4.2_exo_pension_fund.png');

if ~exist(file_tpi_results, 'file'), error('找不到TPI结果文件: %s', file_tpi_results); end
if ~exist(file_iter_results, 'file'), error('找不到TPI迭代文件: %s', file_iter_results); end
if ~exist(output_dir, 'dir'), mkdir(output_dir); fprintf('   已创建输出目录: %s\n', output_dir); end

fprintf('   正在加载基准路径数据...\n');
data_tpi = load(file_tpi_results);
data_iter = load(file_iter_results);

cS = data_tpi.cS;
results = data_tpi.results;
b_hat_path_iter = data_iter.b_hat_path_iter;
T = cS.T_sim;

fprintf('   ✅ 数据加载完毕。\n');

%% --- 3. [核心] 会计核算 ---
fprintf('   正在进行养老金缺口与基金动态的会计核算...\n');

r_path_level = results.r_path;
Y_path_level = results.Y_path;
w_hat_path = results.w_path ./ cS.A_path;
mass_retirees_path = sum(cS.Z_path_raw((cS.aR_new+1):end, :), 1);

b_hat_path_exo = exo_rho * w_hat_path;
b_hat_delta = b_hat_path_exo - b_hat_path_iter;

P_level_path = zeros(1, T);
P_level_path(1) = initial_fund_to_gdp_ratio * Y_path_level(1);
fprintf('   初始(t=1)基金余额设定为: %.4f (基于GDP的%.2f%%)\n', P_level_path(1), initial_fund_to_gdp_ratio*100);

for t = 1:(T-1)
    B_delta_total_level = b_hat_delta(t) * cS.A_path(t) * mass_retirees_path(t);
    P_level_path(t+1) = P_level_path(t) * (1 + r_path_level(t)) - B_delta_total_level;
end

Total_PAYG_Contrib_level = b_hat_path_iter .* cS.A_path .* mass_retirees_path;
Total_Promised_Benefit_level = b_hat_path_exo .* cS.A_path .* mass_retirees_path;

contrib_gdp_ratio = Total_PAYG_Contrib_level ./ Y_path_level;
benefit_gdp_ratio = Total_Promised_Benefit_level ./ Y_path_level;

fprintf('   ✅ 会计核算完成。\n');

%% --- 4. 准备绘图数据 ---
fprintf('   正在准备绘图数据...\n');
paths = struct();
paths.fund_level = P_level_path;
paths.retiree_share = mass_retirees_path ./ sum(cS.Z_path_raw, 1);
paths.contrib_gdp = contrib_gdp_ratio;
paths.benefit_gdp = benefit_gdp_ratio;

time_axis_full = (cS.start_year : cS.time_Step : (cS.start_year + (T-1)*cS.time_Step));
time_mask = time_axis_full <= year_max;
time_axis = time_axis_full(time_mask);

fields = fieldnames(paths);
for i = 1:length(fields)
    if size(paths.(fields{i}), 2) == T
        paths.(fields{i}) = paths.(fields{i})(time_mask);
    end
end
fprintf('   图表将展示到 %d 年。\n', year_max);

%% --- 5. [核心] 绘图 ---
fprintf('   正在生成图表...\n');
fig = figure('Name', '外生社保基金动态与流量分析', 'Position', [461, 603, 730, 319]);
tcl = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

% --- 统一定义灰度风格和格式 ---
font_size_title = 12;
font_size_axis  = 10;
color_axis_main = 'black';
color_axis_bg = [0.5 0.5 0.5];



% ================= PANEL (a): 基金存量动态 =================
ax1 = nexttile;
hold(ax1, 'on');

style_fund = {'k-', 'LineWidth', 2};
style_retiree_share = {'--', 'Color', [0.5 0.5 0.5], 'LineWidth', 2};
style_zero_line = {'--', 'Color', [0.7 0.7 0.7], 'LineWidth', 1};

set(ax1, 'FontName', 'Times New Roman', 'FontSize', font_size_axis);

yyaxis(ax1, 'right');
p_retiree = plot(ax1, time_axis, paths.retiree_share * 100, style_retiree_share{:}, 'DisplayName', '退休人口占比 (右轴, %)');
set(ax1, 'YColor', color_axis_bg, 'YLim', [15 70]);

yyaxis(ax1, 'left');
p_fund = plot(ax1, time_axis, paths.fund_level, style_fund{:}, 'DisplayName', '社保基金余额 (左轴)');
ylim([min(paths.fund_level)*0.9 max(paths.fund_level)*2])
yline(ax1, 0, style_zero_line{:});
set(ax1, 'YColor', color_axis_main);
ylabel(ax1, '基金余额 (万亿元)', 'FontName', 'SimSun');


title_obj1 = title(ax1, '(a) 基金存量动态', 'FontSize', font_size_title, 'FontName', 'SimSun');
legend_obj1 = legend(ax1, [p_fund, p_retiree], 'Location', 'south', 'FontSize', 9, 'Box', 'off', 'FontName', 'SimSun');

grid(ax1, 'on'); box(ax1, 'on');
xlim(ax1, [time_axis(1), time_axis(end)]);
set(ax1, 'XTick', time_axis);
xtickangle(ax1, 40);
hold(ax1, 'off');

% ================= PANEL (b): PAYG收支流量 =================
ax2 = nexttile;
hold(ax2, 'on');

style_contrib = {'k-', 'LineWidth', 2};
style_benefit = {'--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2};
color_deficit_area = [0.8 0.8 0.8]; % 浅灰色区域

set(ax2, 'FontName', 'Times New Roman', 'FontSize', font_size_axis);

deficit_gdp = paths.benefit_gdp - paths.contrib_gdp;
p_area = area(ax2, time_axis, deficit_gdp * 100, ...
    'FaceColor', color_deficit_area, 'FaceAlpha', 0.6, ...
    'EdgeColor', 'none', 'DisplayName', 'PAYG缺口');

p_contrib = plot(ax2, time_axis, paths.contrib_gdp * 100, style_contrib{:}, 'DisplayName', 'PAYG缴费收入');
p_benefit = plot(ax2, time_axis, paths.benefit_gdp * 100, style_benefit{:}, 'DisplayName', '承诺的福利支出');

ylim(ax2, [0, max(paths.benefit_gdp)*100*1.1]);
ylabel(ax2, '占GDP的百分比 (%)', 'FontName', 'SimSun');

title_obj2 = title(ax2, '(b) PAYG收支流量', 'FontSize', font_size_title, 'FontName', 'SimSun');
legend_obj2 = legend(ax2, [p_contrib, p_benefit, p_area], 'Location', 'northwest', 'FontSize', 9, 'Box', 'off', 'FontName', 'SimSun');

grid(ax2, 'on'); box(ax2, 'on');
xlim(ax2, [time_axis(1), time_axis(end)]);
set(ax2, 'XTick', time_axis);
xtickangle(ax2, 40);
hold(ax2, 'off');


fprintf('   ✅ 图表生成完毕。\n');

%% --- 6. 保存图像 ---
fprintf('   正在保存图像至: %s\n', output_filename);
try
    exportgraphics(fig, output_filename, 'Resolution', 600);
    fprintf('   ✅ 图像保存成功。\n');
catch ME
    fprintf('   ❌ 图像保存失败！错误信息: %s\n', ME.message);
end

fprintf('\n--- 脚本执行完毕 ---\n');