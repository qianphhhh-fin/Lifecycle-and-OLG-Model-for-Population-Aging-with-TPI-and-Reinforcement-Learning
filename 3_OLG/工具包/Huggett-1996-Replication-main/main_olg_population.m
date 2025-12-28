%% main_olg_population.m
% 人口结构变化模拟
% 模拟中国人口结构随时间变化的趋势
% 提取自main_olg.m

% 清除工作区
clear;
close all;

% 创建保存图片和数据的目录
if ~exist('fig', 'dir')
    mkdir('fig');
end
if ~exist('data', 'dir')
    mkdir('data');
end

% 设置随机种子以保证结果可重复
rng(42);

%% 模型参数初始化
% 模拟参数
periods = 50;  % 总模拟时期数，每个时期为5年
age_groups = 16;  % 16个五年期年龄组：22-101岁

% 人口参数
x = -0.02;  % 新生代人口增长率调整因子

% 存活概率：从一个年龄组到下一个年龄组的存活概率
% 提高高龄组存活率以增加老龄人口
beta_surv = [0.995, 0.99, 0.985, 0.98, 0.975, 0.97, 0.965, 0.96, ...
    0.95, 0.94, 0.92, 0.89, 0.85, 0.80, 0.75];  % 显著提高高龄存活率

% 退休参数
retirement_age = 8;  % 退休年龄组索引 (对应62岁开始)

% BGP稳态参数
bgp_tolerance = 0.001;  % 判断达到稳态的容差
bgp_window = 5;         % 连续几期满足条件才判定为稳态

%% 初始化人口分布
% 初始化状态变量 - 使用2023年中国人口结构数据
initial_pop = [76.20905535, 86.45596319, 113.8702459, 98.60198303, 86.64117824, 102.7890433, 112.0217783, 99.04620047, ...
              64.05142331, 66.93157492, 44.16815149, 25.40848066, 14.97325553, 6.872421945, 1.743059943, 0.216184341];

% 初始化人口矩阵，行为时期，列为年龄组
Z = zeros(periods+1, age_groups);
Z(1,:) = initial_pop;

%% 模拟人口演变
fprintf('开始模拟人口演变...\n');

% 记录每个时期的主要人口指标
dependency_ratio_history = zeros(periods, 1);  % 抚养比历史
retirees_ratio_history = zeros(periods, 1);    % 退休人口比例历史
young_ratio_history = zeros(periods, 1);       % 年轻工作者比例历史
old_worker_ratio_history = zeros(periods, 1);  % 老年工作者比例历史
total_population_history = zeros(periods, 1);  % 总人口历史
growth_rate_history = zeros(periods, 1);       % 人口增长率历史

% BGP稳态检测变量
pop_structure_ratios = zeros(periods, age_groups);  % 各年龄组占总人口比例
age_group_growth_rates = zeros(periods, age_groups); % 各年龄组增长率
bgp_reached = false;                              % 是否达到BGP稳态
bgp_period = 0;                                   % 达到BGP稳态的时期

% 运行人口演变模拟
for t = 1:periods
    fprintf('计算第 %d 个时期的人口演变...\n', t);
    
    % 更新人口分布
    Z_new = update_population(Z, t, x, beta_surv);
    Z(t+1,:) = Z_new;
    
    % 计算当期的人口相关指标
    dependency_ratio = calculate_dependency_ratio(Z(t,:), retirement_age);
    
    % 计算人口结构变化
    young_workers = sum(Z(t, 1:4));  % 22-41岁
    older_workers = sum(Z(t, 5:retirement_age-1));  % 42-61岁
    retirees = sum(Z(t, retirement_age:end));  % 62-101岁
    total_population = sum(Z(t,:));
    
    young_ratio = young_workers / total_population * 100;
    older_ratio = older_workers / total_population * 100;
    retirees_ratio = retirees / total_population * 100;
    
    % 计算人口增长率
    if t > 1
        growth_rate = (total_population - sum(Z(t-1,:))) / sum(Z(t-1,:)) * 100;
    else
        growth_rate = 0;
    end
    
    % 记录历史数据
    dependency_ratio_history(t) = dependency_ratio;
    retirees_ratio_history(t) = retirees_ratio;
    young_ratio_history(t) = young_ratio;
    old_worker_ratio_history(t) = older_ratio;
    total_population_history(t) = total_population;
    growth_rate_history(t) = growth_rate;
    
    % 计算各年龄组占总人口的比例
    pop_structure_ratios(t,:) = Z(t,:) / total_population;
    
    % 计算各年龄组人口增长率
    if t > 1
        for a = 1:age_groups
            age_group_growth_rates(t,a) = (Z(t,a) - Z(t-1,a)) / Z(t-1,a);
        end
    end
    
    fprintf('  时期 %d: 总人口=%.2f, 抚养比=%.2f, 退休人口比例=%.2f%%\n', ...
        t, total_population, dependency_ratio, retirees_ratio);
end

fprintf('人口演变模拟完成\n');

%% 检测BGP稳态
fprintf('\n开始检测BGP稳态...\n');

% 各年龄组人口结构比例稳定性检测
pop_structure_changes = zeros(periods-1, 1);
for t = 2:periods
    % 计算人口结构变化程度 (欧氏距离)
    pop_structure_changes(t-1) = norm(pop_structure_ratios(t,:) - pop_structure_ratios(t-1,:));
end

% 增长率稳定性检测
growth_rate_changes = zeros(periods-2, 1);
for t = 3:periods
    growth_rate_changes(t-2) = abs(growth_rate_history(t) - growth_rate_history(t-1));
end

% 抚养比稳定性检测
dependency_ratio_changes = zeros(periods-1, 1);
for t = 2:periods
    dependency_ratio_changes(t-1) = abs(dependency_ratio_history(t) - dependency_ratio_history(t-1));
end

% BGP稳态判断 - 连续bgp_window期满足条件
for t = bgp_window+2:periods
    % 检查结构变化、增长率变化和抚养比变化是否都在容差范围内
    structure_stable = all(pop_structure_changes(t-bgp_window:t-1) < bgp_tolerance);
    growth_stable = all(growth_rate_changes(t-bgp_window-1:t-2) < bgp_tolerance);
    dependency_stable = all(dependency_ratio_changes(t-bgp_window:t-1) < bgp_tolerance);
    
    if structure_stable && growth_stable && dependency_stable
        bgp_reached = true;
        bgp_period = t;
        break;
    end
end

if bgp_reached
    fprintf('检测到BGP稳态! 在第 %d 期达到稳态\n', bgp_period);
    fprintf('稳态特征:\n');
    fprintf('  - 人口增长率: %.4f%%\n', growth_rate_history(bgp_period));
    fprintf('  - 抚养比: %.4f\n', dependency_ratio_history(bgp_period));
    fprintf('  - 退休人口比例: %.2f%%\n', retirees_ratio_history(bgp_period));
else
    fprintf('在模拟期限内未检测到BGP稳态\n');
end

%% 绘制人口结构变化图表
% plot_population_dynamics(periods, Z, dependency_ratio_history, ...
%     retirees_ratio_history, young_ratio_history, old_worker_ratio_history, ...
%     total_population_history, growth_rate_history, retirement_age);

%% 绘制BGP稳态检测图表
plot_bgp_detection(periods, pop_structure_changes, growth_rate_changes, ...
    dependency_ratio_changes, bgp_tolerance, bgp_period, bgp_reached);

fprintf('人口结构变化图表已保存到fig文件夹\n');

%% 更新人口分布函数
function Z_new = update_population(Z, t, x, beta_surv)
    % 更新人口分布 - 16个年龄组
    Z_new = zeros(1, 16);
    
    % 最年轻组 (22-26岁) - 新生人口
    % 随时间变化的增长率：前期小负增长，后期大负增长
    growth_rate = x;
    if t < 6
        % 前6期较缓和的负增长
        growth_rate = -0.01 - 0.003 * t;  % 从-1%逐渐增大负增长率
    elseif t >= 6
        % 后期强化负增长
        growth_rate = -0.03 - 0.004 * min(t-6, 10);  % 最终达到-7%的负增长率
    end
    
    Z_new(1) = Z(t, 1) * (1 + growth_rate);
    
    % 其他年龄组 - 从上一组存活下来的人口
    for a = 2:16
        Z_new(a) = Z(t, a-1) * beta_surv(a-1);
    end
end

%% 计算抚养比
function dependency_ratio = calculate_dependency_ratio(Z, retirement_age)
    % 计算抚养比：退休人口/工作人口
    retired_population = sum(Z(retirement_age:end));
    working_population = sum(Z(1:retirement_age-1));
    
    if working_population > 0
        dependency_ratio = retired_population / working_population;
    else
        dependency_ratio = inf;  % 避免除零错误
    end
end

%% 绘制人口结构动态图表
function plot_population_dynamics(periods, Z, dependency_ratio_history, ...
    retirees_ratio_history, young_ratio_history, old_worker_ratio_history, ...
    total_population_history, growth_rate_history, retirement_age)
    
    % 设置中文显示支持
    set(0, 'DefaultAxesFontName', 'SimHei');
    set(0, 'DefaultTextFontName', 'SimHei');
    
    % 准备数据
    t = 1:periods;  % 时间轴
    
    % 1. 人口金字塔演变图 (选择几个关键时间点)
    key_periods = [1, 6, 11, 16, 21, 26];  % 选择展示的时期
    age_labels = {'22-26', '27-31', '32-36', '37-41', '42-46', '47-51', '52-56', '57-61', ...
        '62-66', '67-71', '72-76', '77-81', '82-86', '87-91', '92-96', '97-101'};
    
    figure('Position', [100, 100, 1200, 800]);
    subplot_cols = 3;
    subplot_rows = ceil(length(key_periods)/subplot_cols);
    
    for i = 1:length(key_periods)
        t_idx = key_periods(i);
        if t_idx <= periods
            subplot(subplot_rows, subplot_cols, i);
            barh(Z(t_idx,:), 'FaceColor', [0.4, 0.6, 0.9]);
            set(gca, 'YTick', 1:16, 'YTickLabel', age_labels);
            title(sprintf('第%d期人口分布', t_idx));
            xlabel('人口数量 (百万)');
            grid on;
            
            % 添加退休年龄线
            hold on;
            plot([0, max(Z(t_idx,:))*1.1], [retirement_age-0.5, retirement_age-0.5], 'r--');
            text(max(Z(t_idx,:))*0.7, retirement_age-0.3, '退休年龄', 'Color', 'r');
            hold off;
        end
    end
    
    sgtitle('中国人口金字塔演变', 'FontSize', 16);
    saveas(gcf, 'fig/population_pyramids.png');
    
    % 2. 人口结构比例变化
    figure('Position', [100, 100, 1200, 500]);
    
    % 人口比例堆积图
    subplot(1, 2, 1);
    area(t, [young_ratio_history, old_worker_ratio_history, retirees_ratio_history]);
    colormap summer;
    title('人口结构比例变化', 'FontSize', 14);
    xlabel('时期', 'FontSize', 12);
    ylabel('人口比例 (%)', 'FontSize', 12);
    legend('年轻工作者(22-41岁)', '老年工作者(42-61岁)', '退休者(62+岁)', 'Location', 'best');
    grid on;
    
    % 抚养比变化曲线
    subplot(1, 2, 2);
    yyaxis left;
    plot(t, dependency_ratio_history, 'r-', 'LineWidth', 2);
    ylabel('抚养比', 'Color', 'r', 'FontSize', 12);
    
    yyaxis right;
    plot(t, retirees_ratio_history, 'b-', 'LineWidth', 2);
    ylabel('退休人口比例 (%)', 'Color', 'b', 'FontSize', 12);
    
    title('抚养比与退休人口比例变化', 'FontSize', 14);
    xlabel('时期', 'FontSize', 12);
    legend('抚养比', '退休人口比例', 'Location', 'best');
    grid on;
    
    saveas(gcf, 'fig/population_structure.png');
    
    % 3. 人口总量和增长率变化
    figure('Position', [100, 100, 1000, 500]);
    
    yyaxis left;
    plot(t, total_population_history, 'k-', 'LineWidth', 2);
    ylabel('总人口 (百万)', 'Color', 'k', 'FontSize', 12);
    
    yyaxis right;
    plot(t, growth_rate_history, 'g-', 'LineWidth', 2);
    ylabel('人口增长率 (%)', 'Color', 'g', 'FontSize', 12);
    
    title('总人口与人口增长率变化', 'FontSize', 14);
    xlabel('时期', 'FontSize', 12);
    legend('总人口', '人口增长率', 'Location', 'best');
    grid on;
    
    saveas(gcf, 'fig/population_growth.png');
    
    % 4. 各年龄组人口变化热图
    figure('Position', [100, 100, 1200, 600]);
    
    % 创建热图数据
    heatmap_data = Z(1:periods,:)';
    
    % 绘制热图
    imagesc(t, 1:16, heatmap_data);
    colormap(jet);
    colorbar;
    title('各年龄组人口变化', 'FontSize', 14);
    xlabel('时期', 'FontSize', 12);
    ylabel('年龄组', 'FontSize', 12);
    set(gca, 'YTick', 1:16, 'YTickLabel', age_labels);
    
    % 添加退休年龄线
    hold on;
    plot([0, periods+1], [retirement_age-0.5, retirement_age-0.5], 'w--', 'LineWidth', 2);
    text(periods*0.05, retirement_age-0.3, '退休年龄', 'Color', 'w', 'FontSize', 12);
    hold off;
    
    saveas(gcf, 'fig/population_heatmap.png');
end 

%% 绘制BGP稳态检测图表
function plot_bgp_detection(periods, pop_structure_changes, growth_rate_changes, ...
    dependency_ratio_changes, bgp_tolerance, bgp_period, bgp_reached)
    
    % 设置中文显示支持
    set(0, 'DefaultAxesFontName', 'SimHei');
    set(0, 'DefaultTextFontName', 'SimHei');
    
    % 准备数据
    t1 = 2:periods;  % 结构变化和抚养比变化的时间轴
    t2 = 3:periods;  % 增长率变化的时间轴
    
    % 创建BGP稳态检测图表
    figure('Position', [100, 100, 1200, 800]);
    
    % 1. 人口结构变化
    subplot(3, 1, 1);
    plot(t1, pop_structure_changes, 'b-', 'LineWidth', 2);
    hold on;
    plot([2, periods], [bgp_tolerance, bgp_tolerance], 'r--', 'LineWidth', 1.5);
    if bgp_reached
        plot([bgp_period, bgp_period], [0, max(pop_structure_changes)*1.1], 'g--', 'LineWidth', 1.5);
        text(bgp_period+0.5, max(pop_structure_changes)*0.8, '稳态', 'Color', 'g', 'FontSize', 12);
    end
    hold off;
    title('人口结构变化程度', 'FontSize', 14);
    xlabel('时期', 'FontSize', 12);
    ylabel('人口结构变化', 'FontSize', 12);
    legend('结构变化', '稳态阈值', 'Location', 'best');
    grid on;
    
    % 2. 增长率变化
    subplot(3, 1, 2);
    plot(t2, growth_rate_changes, 'r-', 'LineWidth', 2);
    hold on;
    plot([3, periods], [bgp_tolerance, bgp_tolerance], 'r--', 'LineWidth', 1.5);
    if bgp_reached
        plot([bgp_period, bgp_period], [0, max(growth_rate_changes)*1.1], 'g--', 'LineWidth', 1.5);
    end
    hold off;
    title('人口增长率变化程度', 'FontSize', 14);
    xlabel('时期', 'FontSize', 12);
    ylabel('增长率变化', 'FontSize', 12);
    legend('增长率变化', '稳态阈值', 'Location', 'best');
    grid on;
    
    % 3. 抚养比变化
    subplot(3, 1, 3);
    plot(t1, dependency_ratio_changes, 'm-', 'LineWidth', 2);
    hold on;
    plot([2, periods], [bgp_tolerance, bgp_tolerance], 'r--', 'LineWidth', 1.5);
    if bgp_reached
        plot([bgp_period, bgp_period], [0, max(dependency_ratio_changes)*1.1], 'g--', 'LineWidth', 1.5);
    end
    hold off;
    title('抚养比变化程度', 'FontSize', 14);
    xlabel('时期', 'FontSize', 12);
    ylabel('抚养比变化', 'FontSize', 12);
    legend('抚养比变化', '稳态阈值', 'Location', 'best');
    grid on;
    
    % 总标题
    if bgp_reached
        sgtitle(sprintf('BGP稳态检测 - 第%d期达到稳态', bgp_period), 'FontSize', 16);
    else
        sgtitle('BGP稳态检测 - 未达到稳态', 'FontSize', 16);
    end
    
    saveas(gcf, 'fig/bgp_detection.png');
    
    % 如果达到稳态，绘制稳态结构图
    if bgp_reached
        % 稳态人口金字塔
        age_labels = {'22-26', '27-31', '32-36', '37-41', '42-46', '47-51', '52-56', '57-61', ...
            '62-66', '67-71', '72-76', '77-81', '82-86', '87-91', '92-96', '97-101'};
        
        figure('Position', [100, 100, 1000, 600]);
        barh(Z(bgp_period,:), 'FaceColor', [0.4, 0.7, 0.5]);
        set(gca, 'YTick', 1:16, 'YTickLabel', age_labels);
        title('BGP稳态人口结构', 'FontSize', 16);
        xlabel('人口数量 (百万)', 'FontSize', 14);
        ylabel('年龄组', 'FontSize', 14);
        grid on;
        
        saveas(gcf, 'fig/bgp_steady_state.png');
    end
end 