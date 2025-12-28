% =========================================================================
% == SCRIPT: lifecycle_type.m
% == 版本: [v1.0]
% ==
% == 目的:
% ==   - 绘制不同制度下各类型家庭的生命周期储蓄与消费剖面图。
% ==   - 核心任务: 生成一张2x2的图表 (fig4.4_lifecycle_type)，每个子图
% ==     对应一类家庭，展示其在2050年横截面的生命周期行为。
% ==   - 每个子图包含四条线:
% ==     1. 储蓄率 (无PPS) = k_prime / 到手收入
% ==     2. 消费率 (无PPS) = c / 到手收入
% ==     3. 储蓄率 (有PPS) = k_prime / 到手收入
% ==     4. 消费率 (有PPS) = c / 到手收入
% ==
% == 前置条件:
% ==   - 已运行 `main_run_trans.m` 两次，分别生成了带与不带PPS的转轨结果。
% ==   - TRANS/TPI_results_het_nopps.mat
% ==   - TRANS/TPI_results_het_pps.mat
% =========================================================================

clear; close all; clc;
addpath(pwd);
fprintf('=== OLG异质性模型生命周期剖面分析脚本 (v1.0) ===\n\n');

%% --- 1. 设置与加载数据 ---
fprintf('--- 1. 设置与加载数据 ---\n');

% --- 文件路径定义 ---
file_path_nopps = 'TRANS/TPI_results_het_nopps.mat';
file_path_pps = 'TRANS/TPI_results_het_pps.mat';
output_dir_fig = 'tex/fig/';
output_fig_filename = [output_dir_fig, 'fig4.4_lifecycle_type.png'];

% --- 检查并创建输出目录 ---
if ~exist(output_dir_fig, 'dir'), mkdir(output_dir_fig); end

% --- 加载数据 ---
fprintf('   正在加载 "无PPS" (基准) 和 "有PPS" (改革) 情景数据...\n');
if ~exist(file_path_nopps, 'file'), error('基准情景(nopps)的数据文件不存在: %s', file_path_nopps); end
if ~exist(file_path_pps, 'file'), error('改革情景(pps)的数据文件不存在: %s', file_path_pps); end
data_nopps = load(file_path_nopps);
data_pps = load(file_path_pps);
fprintf('   ✅ 数据加载完成。\n');


%% --- 2. [核心] 计算生命周期剖面 ---
fprintf('\n--- 2. 计算各类型家庭在2050年的生命周期剖面 ---\n');

% --- 统一定义与参数 ---
analysis_year = 2050;
type_labels = {'(a) 高收入城镇职工', '(b) 中低收入城镇职工', '(c) 高收入居民', '(d) 中低收入居民'};
scenarios = {'nopps', 'pps'};
loaded_data = {data_nopps, data_pps};

% --- 初始化结果存储容器 ---
% 维度: aD_new x nTypes x 2 (scenarios)
cS_ref = data_nopps.cS; % 使用一个cS作为参考
saving_rates_all = zeros(cS_ref.aD_new, cS_ref.nTypes, 2);
consumption_rates_all = zeros(cS_ref.aD_new, cS_ref.nTypes, 2);


% --- 遍历两种情景 (无PPS/有PPS) ---
for s = 1:2
    fprintf('   正在处理情景: %s ...\n', scenarios{s});
    
    % --- 为当前情景解包数据 ---
    data_scen = loaded_data{s};
    cS = data_scen.cS;
    final_Pol_path_h = data_scen.final_Pol_path_h;
    final_Dist_path_h = data_scen.final_Dist_path_h;

    % --- 定位到2050年对应的时间索引 t_idx ---
    t_idx = round((analysis_year - cS.start_year) / cS.time_Step) + 1;
    if t_idx > cS.T_sim
        warning('分析年份 %d 超出 %s 情景的模拟范围 T=%d', analysis_year, scenarios{s}, cS.T_sim);
        continue;
    end
    
    % --- 计算该时期的BGP增长因子 ---
    A_path_ext = [cS.A_path, cS.A_path(end) * (1 + ((1+cS.g_A_ss)^cS.time_Step-1))];
    g_A_period_path = A_path_ext(2:end)./A_path_ext(1:end-1)-1;
    growth_factor_bgp = 1 + g_A_period_path(t_idx);
    
    % --- 遍历四种家庭类型 ---
    for h = 1:cS.nTypes
        pol_t_h = final_Pol_path_h{h}{t_idx};
        dist_t_h_full = final_Dist_path_h{h}(:,:,:,:,t_idx);
        
        % --- 遍历所有年龄组 ---
        for ia = 1:cS.aD_new
            pol_age = pol_t_h(ia);
            dist_age = squeeze(dist_t_h_full(:, :, :, ia));
            
            total_mass_age = sum(dist_age, 'all');
            if total_mass_age < 1e-12, continue; end
            
            % 从策略中提取消费和储蓄决策 (hat量)
            c_mat = pol_age.c;
            k_prime_mat = pol_age.k_prime;
            
            % 根据预算约束反解出每个微观状态下的到手收入 (hat量)
            % 到手收入(hat) = (1+tau_c)*c(hat) + k_prime(hat) * (1+g_A)
            disposable_income_mat = (1 + cS.tau_c) * c_mat + k_prime_mat * growth_factor_bgp;
            
            % 聚合得到该年龄组的总量
            total_c_agg = sum(c_mat .* dist_age, 'all');
            total_k_prime_agg = sum(k_prime_mat .* dist_age, 'all');
            total_disposable_income_agg = sum(disposable_income_mat .* dist_age, 'all');
            
            % 计算平均比率
            if total_disposable_income_agg > 1e-9
                consumption_rates_all(ia, h, s) = total_c_agg / total_disposable_income_agg;
                saving_rates_all(ia, h, s) = total_k_prime_agg / total_disposable_income_agg;
            end
        end % 结束年龄循环
    end % 结束类型循环
end % 结束情景循环

fprintf('   ✅ 所有生命周期剖面计算完成。\n');


%% --- 3. 生成 2x2 图表 ---
fprintf('\n--- 3. 生成生命周期剖面图 (fig4.4) ---\n');

fig = figure('Name', '各类型家庭生命周期储蓄与消费剖面', 'Position', [100 495 706 500]);
tcl = tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

% --- 定义灰度绘图风格 ---
% 无PPS: 较深的颜色, 实线(储蓄)+虚线(消费)
style_sr_nopps = {'-o', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.5, 'MarkerSize', 3};
style_cr_nopps = {'--o', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.5, 'MarkerSize', 3};
% 有PPS: 较浅的颜色, 实线(储蓄)+虚线(消费)
style_sr_pps = {'-s', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5, 'MarkerSize', 3};
style_cr_pps = {'--s', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5, 'MarkerSize', 3};

% --- 准备X轴: 物理年龄 ---
phys_ages = cS_ref.age1_orig:cS_ref.time_Step:(cS_ref.age1_orig + (cS_ref.aD_new-1)*cS_ref.time_Step);

for h = 1:cS_ref.nTypes
    ax = nexttile;
    hold(ax, 'on');
    
    % --- 绘制 无PPS 情景 ---
    p1 = plot(ax, phys_ages, saving_rates_all(:, h, 1), style_sr_nopps{:}, 'DisplayName', '储蓄率 (无PPS)');
    p2 = plot(ax, phys_ages, consumption_rates_all(:, h, 1), style_cr_nopps{:}, 'DisplayName', '消费率 (无PPS)');
    
    % --- 绘制 有PPS 情景 ---
    p3 = plot(ax, phys_ages, saving_rates_all(:, h, 2), style_sr_pps{:}, 'DisplayName', '储蓄率 (有PPS)');
    p4 = plot(ax, phys_ages, consumption_rates_all(:, h, 2), style_cr_pps{:}, 'DisplayName', '消费率 (有PPS)');
    
    % --- 标注退休年龄 ---
    retire_line = xline(ax, cS_ref.ageRetire_orig, 'k-.', 'LineWidth', 1, 'DisplayName', '退休年龄', 'HandleVisibility', 'off');
    
    % --- 图形格式化 ---
    title(ax, type_labels{h}, 'FontName', 'SimSun', 'FontSize', 14);
    xlabel(ax, '年龄', 'FontName', 'SimSun');
    if mod(h, 2) == 1 % 只在左侧图表显示Y轴标签
        ylabel(ax, '占当期到手收入比例', 'FontName', 'SimSun');
    end
    grid(ax, 'on');
    box(ax, 'on');
    xlim(ax, [cS_ref.age1_orig, cS_ref.ageLast_orig - cS_ref.time_Step]);
    ylim(ax, [0, 1.0]); % 比例通常在0到1之间

    % --- 在第一个子图中显示图例 ---
    if h == 1
        lgd = legend(ax, [p1, p2, p3, p4], 'Location', 'North', 'FontName', 'SimSun');
        lgd.NumColumns = 2;
    end
end

% --- 保存图像 ---
fprintf('   正在保存图像至: %s\n', output_fig_filename);
try
    exportgraphics(fig, output_fig_filename, 'Resolution', 300);
    fprintf('   ✅ 图像保存成功。\n');
catch ME
    fprintf('   ❌ 图像保存失败！错误信息: %s\n', ME.message);
end

fprintf('\n--- 生命周期剖面分析脚本执行完毕 ---\n');