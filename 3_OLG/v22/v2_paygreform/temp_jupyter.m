% =========================================================================
% == SCRIPT: create_fig_structural_reform.m
% == 版本: [v1.2 - 尺寸与灰度修正版]
% ==
% == 目的:
% ==   1. 加载'structural_reform'情景下的转轨动态结果。
% ==   2. 对最终的分布路径进行后处理，聚合出分类型的PPS资本和总PAYG缴费。
% ==   3. 严格遵循'ref_for_graph.m'的风格，生成一个1x3的图，展示
% ==      从PAYG主导向PPS主导转变的关键宏观和政策路径。
% ==
% == v1.2 核心修改:
% ==   - [!!!] 将图形宽度严格设置为 14.99cm，以满足排版要求。
% ==   - [!!!] 将子图2 (PAYG) 的所有元素统一为灰度风格。
% =========================================================================

clear; close all; clc;
addpath(pwd);
fprintf('=== 结构性养老金改革(structural_reform)可视化脚本 (v1.2) ===\n\n');

%% --- 1. 用户设定与数据加载 ---
fprintf('--- 1. 设定与文件加载 ---\n');

% --- 输入文件 ---
input_file = 'TRANS/TPI_results_het_pps_structural_reform.mat';

% --- 输出路径 ---
output_dir = 'tex/fig/';
if ~exist(output_dir, 'dir'), mkdir(output_dir); end
output_filename = fullfile(output_dir, 'fig4.8_pps_structural_reform.png');

% --- 绘图设定 ---
plot_start_year = 2023;
plot_end_year = 2100;

% --- 加载数据 ---
fprintf('   正在加载文件: %s ...\n', input_file);
if ~exist(input_file, 'file'), error('找不到TPI结果文件: %s', input_file); end
data = load(input_file);
cS = data.cS;
results = data.results;
dist_path_h = data.final_Dist_path_h;
paramSF = data.results.ss_data.paramSF;
fprintf('   ✅ 数据加载完成。\n');


%% --- 2. [核心] 数据后处理与计算 ---
fprintf('\n--- 2. 正在进行后处理计算 ---\n');
T = cS.T_sim;
nH = cS.nTypes;
Y_path = results.Y_path;

% --- 2.1 聚合分类型的PPS资本路径 ---
fprintf('   2.1 聚合分类型PPS资本...\n');
[Kpps_hat_path_h, Kpps_path_h] = aggregate_pps_capital_by_type(dist_path_h, cS);
Kpps_to_GDP_path_h = Kpps_path_h ./ Y_path * 100; % 转换为百分比
fprintf('       ✅ 完成。\n');

% --- 2.2 聚合总PAYG缴费路径 ---
fprintf('   2.2 聚合总PAYG缴费...\n');
PAYG_Contrib_hat_path = aggregate_payg_contributions(dist_path_h, results.w_path, cS, paramSF);
PAYG_Contrib_path = PAYG_Contrib_hat_path .* cS.A_path;
PAYG_to_GDP_path = PAYG_Contrib_path ./ Y_path * 100; % 转换为百分比
fprintf('       ✅ 完成。\n');

% --- 2.3 计算宏观比率和价格路径 ---
fprintf('   2.3 计算宏观比率与价格路径...\n');
Kp_to_Y_path = results.K_p_path ./ Y_path;
C_to_Y_path = results.C_path ./ Y_path;
w_path = results.w_path;
r_annual_path = ((1 + results.r_path).^(1/cS.time_Step) - 1) * 100; % 年化利率
fprintf('       ✅ 完成。\n');

%% --- 3. [核心] 可视化 ---
fprintf('\n--- 3. 正在生成1x3图表 ---\n');

% --- 绘图风格设定 (遵循 ref_for_graph.m) ---
% [!!! 核心修改 !!!] 严格设置图形尺寸
fig = figure('Name', '结构性养老金改革路径分析');
fig.Units = 'centimeters';
fig.Position = [10 10 14.99 5.5]; % [left bottom width height]

tcl = tiledlayout(1, 3, 'Padding', 'tight', 'TileSpacing', 'compact');

time_axis = cS.start_year : cS.time_Step : (cS.start_year + (T-1)*cS.time_Step);
time_mask = (time_axis >= plot_start_year) & (time_axis <= plot_end_year);
time_axis_plot = time_axis(time_mask);

% --- 子图1: PPS 资本 与 缴费上限 ---
ax1 = nexttile;
hold(ax1, 'on');

bar_data = Kpps_to_GDP_path_h(:, time_mask)';
b1 = bar(ax1, time_axis_plot, bar_data, 'stacked');

colors_bar = [0.2 0.2 0.2; 0.45 0.45 0.45; 0.65 0.65 0.65; 0.85 0.85 0.85];
for i = 1:nH
    b1(i).FaceColor = colors_bar(i, :);
end

plot(ax1, time_axis_plot, cS.pps_fixed_path(time_mask) * 100, 'k-s', 'LineWidth', 0.3, 'MarkerFaceColor', 'k','MarkerSize',3);
ax1.YAxis.Color = 'k';

grid(ax1, 'on'); box(ax1, 'on');
xlim(ax1, [plot_start_year, plot_end_year]);
title(ax1, '（a）个人养老金路径（%）',  'FontSize', 9);
xlabel(ax1, '年份',  'FontSize', 9);
l1=legend(ax1, '高收城镇', '中低收城镇', '高收居民', '中低收居民', '缴费上限', 'Position',[0.0275 0.5513 0.1070 0.3654],...
    'IconColumnWidth',5,'Box','off','FontSize', 9);

% --- 子图2: PAYG 缴费 与 政策路径 ---
ax2 = nexttile;
hold(ax2, 'on');

% [!!! 核心修改 !!!] 统一为灰度图
p0=bar(ax2, time_axis_plot, PAYG_to_GDP_path(time_mask), 'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none', 'DisplayName','PAYG/GDP');
p1 = plot(ax2, time_axis_plot, mean(cS.theta_path_h(:, time_mask),1) * 100, '--', 'Color', [0.0 0.0 0.0], 'LineWidth', 1, 'DisplayName', 'PAYG费率');
p2 = plot(ax2, time_axis_plot, results.adj_path(time_mask), ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.5, 'DisplayName', '调整因子');

ax2.YAxis.Color = 'k';

grid(ax2, 'on'); box(ax2, 'on');
xlim(ax2, [plot_start_year, plot_end_year]);
title(ax2, '（b）PAYG路径（%）',  'FontSize', 9);
xlabel(ax2, '年份');
l2=legend(ax2, [p0, p1, p2], 'Position',[0.4267 0.7460 0.1764 0.1562],...
    'Box','off','IconColumnWidth',10,'FontSize', 9);


% --- 子图3: 关键宏观变量路径 ---
ax3 = nexttile;
hold(ax3, 'on');

yyaxis(ax3, 'left');
p_ky = plot(ax3, time_axis_plot, Kp_to_Y_path(time_mask), '-', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.5, 'DisplayName', 'K/Y（左）');
p_cy = plot(ax3, time_axis_plot, C_to_Y_path(time_mask), '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5, 'DisplayName', 'C/Y（左）');
ax3.YAxis(1).Color = 'k';

yyaxis(ax3, 'right');
% p_w = plot(ax3, time_axis_plot, w_path(time_mask), '-.', 'Color', [0.4 0.4 0.4], 'LineWidth', 2, 'DisplayName', '工资 (右轴)');
p_r = plot(ax3, time_axis_plot, r_annual_path(time_mask), ':', 'Color', [0.0 0.0 0.0], 'LineWidth', 1.5, 'DisplayName', '利率（右,%）');
ax3.YAxis(2).Color = 'k';

grid(ax3, 'on'); box(ax3, 'on');
xlim(ax3, [plot_start_year, plot_end_year]);
title(ax3, '（c）宏观经济路径', 'FontSize', 9);
xlabel(ax3, '年份',  'FontSize', 9);
l3=legend(ax3, 'show', 'Position',[0.8355 0.3860 0.1035 0.2260],...
    'IconColumnWidth',5,'Box','off', 'FontSize', 8);


%% --- 4. 保存图表 ---
fprintf('\n--- 4. 正在保存图表 ---\n');
try
    exportgraphics(fig, output_filename, 'Resolution', 300);
    fprintf('   ✅ 成功保存图表至: %s\n', output_filename);
catch ME
    fprintf('   ❌ 保存失败: %s\n', ME.message);
end

fprintf('\n--- 脚本执行完毕 ---\n');


%% --- [本地辅助函数] ---

function [Kpps_hat_path_h, Kpps_path_h] = aggregate_pps_capital_by_type(dist_path_h, cS)
    % 聚合每个类型家庭的PPS资本存量路径 (hat量 和 level量)
    nH = cS.nTypes;
    T = cS.T_sim;
    Kpps_hat_path_h = zeros(nH, T);

    kpps_grid_full = repmat(reshape(cS.kppsGridV, [1, cS.nkpps]), [cS.nk, 1, cS.nw_expanded]);

    for t = 1:T
        for h = 1:nH
            dist_t_h = dist_path_h{h}(:,:,:,:,t);
            kpps_hat_t_h = 0;
            for ia = 1:cS.aD_new
                dist_slice = dist_t_h(:,:,:,ia);
                kpps_hat_t_h = kpps_hat_t_h + sum(kpps_grid_full .* dist_slice, 'all');
            end
            Kpps_hat_path_h(h, t) = kpps_hat_t_h;
        end
    end
    Kpps_path_h = Kpps_hat_path_h .* cS.A_path;
end

function PAYG_Contrib_hat_path = aggregate_payg_contributions(dist_path_h, w_path, cS, paramSF)
    % 聚合经济体总的PAYG缴费路径 (hat量)
    nH = cS.nTypes;
    T = cS.T_sim;
    PAYG_Contrib_hat_path = zeros(1, T);
    w_hat_path = w_path ./ cS.A_path;

    for t = 1:T
        total_contrib_t_hat = 0;
        aR_t = cS.aR_new_path(t); % 使用时变退休年龄
        for h = 1:nH
            dist_t_h = dist_path_h{h}(:,:,:,:,t);
            cS_h = cS;
            cS_h.ageEffV_new = cS.ageEffV_new_h(:,h);
            
            % 聚合类型h在t期的劳动供给
            L_h_t = 0;
            for ia = 1:aR_t
                mass_ia = dist_t_h(:,:,:,ia);
                if sum(mass_ia(:)) < 1e-30, continue; end
                le_grid_slice = reshape(paramSF.leGridV, [1, 1, cS.nw_expanded]);
                labor_supply_slice = cS_h.ageEffV_new(ia) .* le_grid_slice;
                L_h_t = L_h_t + sum(labor_supply_slice .* mass_ia, 'all');
            end
            
            % 计算类型h在t期的总缴费
            total_contrib_t_hat = total_contrib_t_hat + cS.theta_path_h(h,t) * w_hat_path(t) * L_h_t;
        end
        PAYG_Contrib_hat_path(t) = total_contrib_t_hat;
    end
end