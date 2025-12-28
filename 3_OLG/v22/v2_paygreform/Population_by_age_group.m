% --- 1. 从 Excel 文件加载和准备所有年份的数据 ---

% 定义文件名和关键参数
filename = 'China_population_by_age_1950_2100.xlsx';

% 假设Excel中预测区间的变量名为 'High' 和 'Low'
% 如果实际名称不同 (例如 'Upper', 'Lower' 或中文)，请修改下面这两行
variant_median = 'Median';
variant_upper = 'High';
variant_lower = 'Low';

% 定义预测区间的起始年份 (与原脚本逻辑一致)
projection_start_year = 2020;

% 读取Excel数据
try
    data_table = readtable(filename);
catch
    error('无法读取Excel文件 "%s"。请确保文件与此脚本位于同一文件夹中，或文件在MATLAB的搜索路径下。', filename);
end

% 将人口单位从“个”转换为“百万”
data_table.Value = data_table.Value / 1e6;

% --- 为每个年龄段提取完整的时间序列数据 ---

% 定义一个帮助函数，用于按年龄和类型筛选数据
% 这会让后续代码更整洁
filter_data = @(age, variant) data_table(strcmp(data_table.Age, age) & strcmp(data_table.Variant, variant), :);

% 提取 0-14岁 的数据
median_0_14 = filter_data('0-14', variant_median);
upper_0_14  = filter_data('0-14', variant_upper);
lower_0_14  = filter_data('0-14', variant_lower);

% 提取 15-24岁 的数据
median_15_24 = filter_data('15-24', variant_median);
upper_15_24  = filter_data('15-24', variant_upper);
lower_15_24  = filter_data('15-24', variant_lower);

% 提取 25-64岁 的数据
median_25_64 = filter_data('25-64', variant_median);
upper_25_64  = filter_data('25-64', variant_upper);
lower_25_64  = filter_data('25-64', variant_lower);

% 提取 65岁以上 的数据
median_65_plus = filter_data('65+', variant_median);
upper_65_plus  = filter_data('65+', variant_upper);
lower_65_plus  = filter_data('65+', variant_lower);

% 检查数据是否成功提取
if isempty(median_0_14) || isempty(upper_0_14) || isempty(lower_0_14)
    error(['未能从Excel文件中成功提取所需数据。请检查：\n' ...
           '1. Excel文件中的列名是否为 Time, Age, Variant, Value。\n' ...
           '2. Variant列是否包含 %s, %s, %s 等字符串。\n' ...
           '3. Age列是否包含 0-14, 15-24, 25-64, 65+ 等字符串。'], ...
           variant_median, variant_upper, variant_lower);
end


% --- 2. 绘图设置 (此部分及之后保持完全不变) ---

% 创建一个新的图形窗口
fig = figure('Name', '中国人口年龄结构', 'Color', 'white');

% *** 设置输出图片尺寸为不超过14.99cm ***
set(fig, 'PaperUnits', 'centimeters');
% 设置宽度为14.99cm，高度按比例设为10cm，可自行调整
set(fig, 'PaperPosition', [0 0 14.99 5]);
set(fig,"Position",[680 641 628 237])

hold on; % 允许在同一坐标轴上绘制多个图形


% --- 3. 绘制95%预测区间 (半透明填充) ---

% 筛选出预测年份的数据
proj_years_idx = upper_0_14.Time >= projection_start_year;
proj_years = upper_0_14.Time(proj_years_idx)'; % 转置为行向量

% 定义填充函数简化代码
plot_fill = @(lower_data, upper_data, color) ...
    fill([proj_years, fliplr(proj_years)], ...
         [lower_data.Value(proj_years_idx)', fliplr(upper_data.Value(proj_years_idx)')], ...
         color, 'EdgeColor', 'none', 'FaceAlpha', 0.7);

% 绘制各个年龄段的填充区域
plot_fill(lower_0_14, upper_0_14, [1 0.85 0.85]);
plot_fill(lower_15_24, upper_15_24, [0.85 1 0.85]);
plot_fill(lower_25_64, upper_25_64, [0.85 0.9 1]);
plot_fill(lower_65_plus, upper_65_plus, [0.9 0.85 1]);


% --- 4. 绘制主要的人口趋势线 (直接绘制所有年份的数据，不再插值) ---
plot(median_0_14.Time, median_0_14.Value, 'Color', [0.85, 0.1, 0.1], 'LineWidth', 2.5);
plot(median_15_24.Time, median_15_24.Value, 'Color', [0.1, 0.6, 0.2], 'LineWidth', 2.5);
plot(median_25_64.Time, median_25_64.Value, 'Color', [0, 0.4470, 0.7410], 'LineWidth', 2.5);
plot(median_65_plus.Time, median_65_plus.Value, 'Color', [0.4940, 0.1840, 0.5560], 'LineWidth', 2.5);


% --- 5. 图表美化与设置 ---
ax = gca;
ax.FontSize = 12;
ax.XLim = [1950, 2100];
ax.YLim = [0, 900];
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.GridLineStyle = ':';
ax.GridAlpha = 0.5;
ax.Layer = 'top'; % 让网格线显示在线条上方
box on;

xlabel('年份', 'FontSize', 14);
ylabel('人口 (百万)', 'FontSize', 14);


% --- 6. 在曲线上直接添加文字标注 ---
% 注意：由于现在绘制的是实际数据，曲线位置可能微调，这些标注的坐标可能需要手动调整以达到最佳效果
t1=text(1.9857e+03,420.0599, '0 - 14', 'FontSize', 13, 'FontWeight', 'bold', 'Color', [0.85, 0.1, 0.1]);
t2=text(1.9717e+03,147.6048, '15 - 24', 'FontSize', 13, 'FontWeight', 'bold', 'Color', [0.1, 0.6, 0.2]);
t3=text(2010, 742.2156, '25 - 64', 'FontSize', 13, 'FontWeight', 'bold', 'Color', [0, 0.4470, 0.7410]);
t4=text(2.0072e+03,60.3892, '65+', 'FontSize', 13, 'FontWeight', 'bold', 'Color', [0.4940, 0.1840, 0.5560]);


% --- 7. 创建一个专门的图例来标注预测区间 ---
% 创建一个不可见的代理对象，只为显示图例
p_legend = patch(NaN, NaN, [0.85 0.85 0.85], 'FaceAlpha', 0.8);
legend(p_legend, '95% 预测区间', 'Location', 'northeast', 'FontSize', 13);
hold off;

% 提示：如果要保存图片，可以使用以下命令
% print('China_Population_Chart_Full_Data', '-dpng', '-r300'); % 保存为300 DPI的PNG图片