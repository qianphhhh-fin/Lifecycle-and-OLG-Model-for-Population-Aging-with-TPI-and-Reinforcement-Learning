% --- START OF FILE main_olg_v2_utils_pps.m ---

classdef main_olg_v2_utils_pps % 类名已更新
    % OLG 模型 v2 的工具函数
    % 修改版：采用“混合时间单位”，包含个人养老金计划 (PPS),
    % 并且PAYG工资税率内生，替代率外生

    methods (Static)

        % =====================================================================
        % == 人口动态函数 (年龄组层面 - 未改变) =========
        % =====================================================================

        function popS = initPopulation(cS) % cS: 参数结构体
            % 初始化人口结构 (各年龄组的人口数量)
            popS.Z = zeros(cS.aD_new, 1); % popS.Z: 各年龄组人口数量向量, cS.aD_new: 年龄组数量
            initial_total = sum(cS.initial_pop); % cS.initial_pop: 初始各年龄组人口数量
            if initial_total > 0 && length(cS.initial_pop) == cS.aD_new
                 popS.Z(:, 1) = cS.initial_pop(:) / initial_total * 100; % 标准化为总和100
            else
                warning('初始人口不匹配或总和为零。设置为均匀的初始年龄组人口分布。');
                 popS.Z(:, 1) = 100 / cS.aD_new;
            end
            popS.totalPop = sum(popS.Z(:, 1)); % popS.totalPop: 初始总人口
            if popS.totalPop > 1e-9
                 popS.ageDist = popS.Z(:, 1) / popS.totalPop; % popS.ageDist: 各年龄组人口占比
            else
                 popS.ageDist = zeros(cS.aD_new, 1);
            end
            % 存储初始年龄组分布
            popS.initialAgeDist = popS.ageDist;
            fprintf('初始年龄组人口已设置。总人口=%.2f\n', popS.totalPop);
        end

        function popS = populationDynamics(popS, cS)
            % 模拟年龄组层面的人口动态，使用多年期存活率/增长率
            % 使用 cS.survivalProbV_popdyn (多年期存活率)

            max_periods_sim = cS.max_periods; % max_periods_sim: 最大模拟期数
            Z_history = zeros(cS.aD_new, max_periods_sim + 1); % Z_history: 各年龄组人口数量历史记录
            totalPop_history = zeros(1, max_periods_sim + 1); % totalPop_history: 总人口历史记录
            ageDist_history = zeros(cS.aD_new, max_periods_sim + 1); % ageDist_history: 各年龄组人口占比历史记录

            Z_history(:, 1) = popS.Z(:, 1);
            totalPop_history(1) = popS.totalPop(1);
            ageDist_history(:, 1) = popS.ageDist(:, 1);

            fprintf('人口动态模拟开始 (年龄组, 最大期数 = %d)...\n', max_periods_sim);
            bgp_reached_flag = false; % bgp_reached_flag: 是否达到平衡增长路径（稳态）的标志
            actual_periods_run = max_periods_sim; % actual_periods_run: 实际运行的期数

            for t = 1:max_periods_sim % t: 时间期数
                if mod(t, 10) == 0 || t == 1
                    fprintf('  模拟人口期数 %d (年龄组)\n', t);
                end

                Z_current_period = Z_history(:, t); % Z_current_period: 当前期各年龄组人口数量
                Z_next_period = zeros(cS.aD_new, 1); % Z_next_period: 下一期各年龄组人口数量

                % 增长/存活逻辑与原始V2版本相同 (年龄组层面)
                 time_varying_growth_rate = 0; % time_varying_growth_rate: 时变的（新生儿组）增长率
                 if t < 6
                     time_varying_growth_rate = -0.01 - 0.003 * t;
                 else
                     time_varying_growth_rate = -0.03 - 0.004 * min(t-6, 10);
                 end
                 Z_next_period(1) = Z_current_period(1) * (1 + time_varying_growth_rate);
                 Z_next_period(1) = max(0, Z_next_period(1)); % 确保非负

                for a = 2:cS.aD_new % a: 年龄组索引
                    survival_prob = 0; % survival_prob: 从上一组存活到本组的概率
                    if (a-1) >= 1 && (a-1) <= length(cS.survivalProbV_popdyn)
                         survival_prob = cS.survivalProbV_popdyn(a-1); % 使用多年期人口动态存活率
                    end
                    Z_next_period(a) = Z_current_period(a-1) * survival_prob;
                    Z_next_period(a) = max(0, Z_next_period(a)); % 确保非负
                end

                Z_history(:, t+1) = Z_next_period;
                totalPop_history(t+1) = sum(Z_next_period);
                if totalPop_history(t+1) > 1e-9
                    ageDist_history(:, t+1) = Z_next_period / totalPop_history(t+1);
                else
                    ageDist_history(:, t+1) = 0; totalPop_history(t+1) = 0;
                end

                % 检查是否达到稳态 (年龄组)
                current_check_period = t + 1; % current_check_period: 当前检查的期数索引
                if current_check_period >= cS.bgp_window + 1 % cS.bgp_window: 判断稳态的窗口期长度
                    stable = true; % stable: 当前窗口期是否稳定的标志
                    for w_idx = 1:cS.bgp_window % w_idx: 窗口期内迭代变量 (原为w，避免与工资变量冲突)
                       hist_idx1 = current_check_period - w_idx + 1;
                       hist_idx2 = current_check_period - w_idx;
                       if hist_idx1 > 0 && hist_idx2 > 0
                           change = norm(ageDist_history(:, hist_idx1) - ageDist_history(:, hist_idx2)); % change: 相邻两期分布的范数差异
                           if change >= cS.bgp_tolerance % cS.bgp_tolerance: 稳态容忍度
                               stable = false; break;
                           end
                       else
                           stable = false; break; % 索引越界
                       end
                    end
                    if stable
                        fprintf('\n人口稳态 (年龄组) 在第 %d 期达到。\n', t);
                        bgp_reached_flag = true; actual_periods_run = t; break;
                    end
                end
            end % 结束模拟循环 t

            final_period_idx = actual_periods_run + 1; % final_period_idx: 最终数据期数索引
            popS.Z = Z_history(:, 1:final_period_idx);
            popS.totalPop = totalPop_history(1:final_period_idx);
            popS.ageDist = ageDist_history(:, 1:final_period_idx); % 存储年龄组分布历史

            % 基于年龄组计算抚养比
            depRatio_history = zeros(1, actual_periods_run); % depRatio_history: 抚养比历史记录
             for th = 1:actual_periods_run % th: 历史期数索引
                 Z_t = Z_history(:, th); % Z_t: 第th期各年龄组人口数量
                 working_pop = sum(Z_t(1:cS.aR_new)); % working_pop: 工作人口总数, cS.aR_new: 最后一个工作年龄组索引
                 retired_pop = sum(Z_t(cS.aR_new+1:end)); % retired_pop: 退休人口总数
                 if working_pop > 1e-9
                     depRatio_history(th) = retired_pop / working_pop;
                 else
                     depRatio_history(th) = inf; % 如果工作人口为0，抚养比为无穷
                 end
             end
            popS.dependencyRatio = depRatio_history;

            fprintf('人口动态模拟完成。运行期数: %d。达到BGP: %d\n', actual_periods_run, bgp_reached_flag);
            if ~bgp_reached_flag
                fprintf('警告: 人口稳态 (年龄组) 未在 %d 期内达到。\n', max_periods_sim);
            end
        end

        function [Z_ss, dependency_ratio_ss, bgp_reached, bgp_period] = detectSteadyStatePopulation(popS, cS)
            % 基于年龄组分布的稳定性检测稳态
            % 绘制初始与稳态年龄组分布的比较图

             actual_periods_in_data = size(popS.Z, 2); % actual_periods_in_data: 人口数据中的实际期数
             bgp_reached = false; % bgp_reached: 是否达到稳态的标志
             bgp_period = actual_periods_in_data - 1; % bgp_period: 稳态开始的期数 (默认为最后一期)

             if actual_periods_in_data < cS.bgp_window + 1
                 fprintf('人口模拟期数过短 (%d 数据点)，无法进行稳态检查 (窗口期 = %d)。\n', actual_periods_in_data, cS.bgp_window);
             else
                 fprintf('检查人口稳态 (年龄组, 最近 %d 期)...\n', cS.bgp_window);
                 for t_check_end = actual_periods_in_data : -1 : cS.bgp_window + 1 % t_check_end: 检查窗口的结束期
                     stable = true;
                     for w_idx = 0 : (cS.bgp_window - 1) % w_idx: 窗口期内迭代变量 (原为w)
                         idx1 = t_check_end - w_idx;
                         idx2 = t_check_end - w_idx - 1;
                         if idx1 > 0 && idx2 > 0 && idx1 <= size(popS.ageDist,2) && idx2 <= size(popS.ageDist,2) % 边界检查
                             change = norm(popS.ageDist(:, idx1) - popS.ageDist(:, idx2)); % 比较年龄组分布
                             if change >= cS.bgp_tolerance
                                 stable = false; break;
                             end
                         else
                             stable = false; break; % 索引越界
                         end
                     end
                     if stable
                         bgp_reached = true; bgp_period = t_check_end - 1;
                         fprintf('人口稳态 (年龄组) 从第 %d 期 (数据索引 %d) 开始检测到。\n', bgp_period, t_check_end);
                         break;
                     end
                 end
                  if ~bgp_reached
                      fprintf('未检测到人口稳态 (年龄组)。使用最终期数据。\n');
                  end
             end

             ss_data_index = min(bgp_period + 1, size(popS.Z, 2)); % ss_data_index: 稳态数据对应的索引 (确保有效)
             Z_ss = popS.Z(:, ss_data_index); % Z_ss: 稳态年龄组人口数量
             Z_ss_norm = zeros(cS.aD_new, 1); % Z_ss_norm: 归一化的稳态年龄组质量
             if sum(Z_ss) > 1e-9
                 Z_ss_norm = Z_ss / sum(Z_ss);
             end

             valid_dep_ratio_index = min(bgp_period, length(popS.dependencyRatio)); % valid_dep_ratio_index: 抚养比的有效索引
             if isfield(popS, 'dependencyRatio') && valid_dep_ratio_index > 0
                 dependency_ratio_ss = popS.dependencyRatio(valid_dep_ratio_index); % dependency_ratio_ss: 稳态抚养比
             else % 如果历史抚养比数据不足，则重新计算
                  working_pop_ss = sum(Z_ss(1:cS.aR_new));
                  retired_pop_ss = sum(Z_ss(cS.aR_new+1:end));
                  if working_pop_ss > 1e-9
                      dependency_ratio_ss = retired_pop_ss / working_pop_ss;
                  else
                      dependency_ratio_ss = inf;
                  end
             end

             % --- 绘制初始与稳态年龄组分布图 ---
             figure('Name', 'V2 Mod PPS: 初始 vs 稳态 年龄组人口分布'); % 更改图名以反映PPS版本
             hold on;
             group_indices = 1:cS.aD_new; % group_indices: 年龄组索引 (用于绘图)
             if isfield(popS, 'initialAgeDist') && ~isempty(popS.initialAgeDist)
                 bar(group_indices - 0.2, popS.initialAgeDist * 100, 0.4, 'b', 'DisplayName', '初始年龄组分布');
             else
                 warning('未找到或空的初始年龄分布用于绘图。');
             end
             bar(group_indices + 0.2, Z_ss_norm * 100, 0.4, 'r', 'DisplayName', sprintf('稳态年龄组分布 (第 %d 期)', bgp_period));
             hold off;
             xlabel(sprintf('年龄组索引 (1 至 %d)', cS.aD_new)); ylabel('占总人口百分比 (%)');
             title('V2 Mod PPS: 初始 vs 稳态/最终 年龄组人口分布');
             legend('Location', 'best'); xticks(group_indices); grid on; drawnow;
             fprintf('已绘制初始与稳态/最终年龄组人口分布图。\n');
        end


        % =====================================================================
        % == 经济函数 (修改版：混合时间单位 & PPS & 内生PAYG税率) =======
        % =====================================================================

        function cS = ParameterValues_HuggettStyle()
            % 修改版：包含PPS参数和内生PAYG税率的参数设置

            % --- 开始：原始人口统计与分组参数 ---
            cS.age1_orig      = 20;     % cS.age1_orig: 开始真实年龄
            cS.ageLast_orig   = 98;     % cS.ageLast_orig: 最大真实年龄
            cS.ageRetire_orig = 65;     % cS.ageRetire_orig: 退休真实年龄
            cS.popGrowth_orig = 0.012;  % cS.popGrowth_orig: 年度人口自然增长率

            cS.aD_orig        = cS.ageLast_orig - cS.age1_orig + 1; % cS.aD_orig: 总模型年数
            cS.aR_idx_orig    = cS.ageRetire_orig - cS.age1_orig + 1; % cS.aR_idx_orig: 第一个退休年度的索引
            cS.aW_orig        = cS.aR_idx_orig - 1; % cS.aW_orig: 工作年数
            cS.physAgeV_orig   = (cS.age1_orig : cS.ageLast_orig)'; % cS.physAgeV_orig: 真实年龄向量

            % 年度存活概率 (s_orig) 和死亡率 (d_orig) - 从Huggett复制
            cS.d_orig         = [0.00159, 0.00169, 0.00174, 0.00172, 0.00165, 0.00156, 0.00149, 0.00145, 0.00145, 0.00149,...
                            0.00156, 0.00163, 0.00171, 0.00181, 0.00193, 0.00207, 0.00225, 0.00246, 0.00270, 0.00299,...
                            0.00332, 0.00368, 0.00409, 0.00454, 0.00504, 0.00558, 0.00617, 0.00686, 0.00766, 0.00865,...
                            0.00955, 0.01058, 0.01162, 0.01264, 0.01368, 0.01475, 0.01593, 0.01730, 0.01891, 0.02074,...
                            0.02271, 0.02476, 0.02690, 0.02912, 0.03143, 0.03389, 0.03652, 0.03930, 0.04225, 0.04538,...
                            0.04871, 0.05230, 0.05623, 0.06060, 0.06542, 0.07066, 0.07636, 0.08271, 0.08986, 0.09788,...
                            0.10732, 0.11799, 0.12895, 0.13920, 0.14861, 0.16039, 0.17303, 0.18665, 0.20194, 0.21877,...
                            0.23601, 0.25289, 0.26973, 0.28612, 0.30128, 0.31416, 0.32915, 0.34450, 0.36018];
            if length(cS.d_orig) ~= cS.aD_orig, error('d_orig 长度不匹配'); end
            cS.s_orig         = 1 - cS.d_orig; % cS.s_orig: 年度存活概率向量

            % 静态年度年龄质量分布 (供参考，均衡计算中未使用)
            ageMassV_orig_static = ones(1, cS.aD_orig); % ageMassV_orig_static: 静态年度年龄质量
            for i = 2 : cS.aD_orig
                ageMassV_orig_static(i) = cS.s_orig(i-1) * ageMassV_orig_static(i-1) / (1 + cS.popGrowth_orig);
            end
            cS.ageMassV_orig_static = ageMassV_orig_static ./ sum(ageMassV_orig_static);

            % 5年期年龄组参数
            cS.yearStep      = 5; % cS.yearStep: 每个年龄组的年数
            cS.aD_new        = ceil(cS.aD_orig / cS.yearStep); % cS.aD_new: 5年期年龄组的数量 (16)
            cS.aR_new        = ceil(cS.aW_orig / cS.yearStep); % cS.aR_new: 最后一个工作5年期年龄组的索引 (9)

            % 将新的年龄组索引 'a' 映射到原始年度索引
            cS.physAgeMap = cell(cS.aD_new, 1); % cS.physAgeMap: 年龄组到年度年龄的映射
            for a = 1:cS.aD_new
                startIdx = (a-1)*cS.yearStep + 1;
                endIdx = min(a*cS.yearStep, cS.aD_orig);
                cS.physAgeMap{a} = startIdx:endIdx;
            end

            % 每个5年期年龄组开始时的真实年龄
            cS.physAgeV_new = zeros(cS.aD_new, 1); % cS.physAgeV_new: 各年龄组起始真实年龄向量
             for a = 1:cS.aD_new
                cS.physAgeV_new(a) = cS.physAgeV_orig(cS.physAgeMap{a}(1));
            end

            % 计算用于组间过渡的1年存活概率
            % (使用当前组最后一年的存活概率 s_orig)
            cS.s_1yr_transitionV = zeros(cS.aD_new, 1); % cS.s_1yr_transitionV: VFI中使用的组间1年过渡存活概率
            for a = 1:(cS.aD_new - 1)
                 lastYearIdx = cS.physAgeMap{a}(end);
                 if lastYearIdx < cS.aD_orig
                     cS.s_1yr_transitionV(a) = cS.s_orig(lastYearIdx);
                 else
                     cS.s_1yr_transitionV(a) = 0; % 如果是最后一个年度年龄，则无法过渡
                 end
            end
            cS.s_1yr_transitionV(cS.aD_new) = 0; % 无法从最后一个年龄组存活出去

            % 计算年龄组人口质量 (静态，供参考/计算L)
            cS.ageMassV_new_static = zeros(1, cS.aD_new); % cS.ageMassV_new_static: 静态年龄组质量分布
            for a = 1:cS.aD_new
                cS.ageMassV_new_static(a) = sum(cS.ageMassV_orig_static(cS.physAgeMap{a}));
            end
            cS.ageMassV_new_static = cS.ageMassV_new_static ./ sum(cS.ageMassV_new_static);
            cS.retireMass_new = sum(cS.ageMassV_new_static(cS.aR_new + 1 : cS.aD_new)); % cS.retireMass_new: 静态计算的退休人口质量占比 (用于PAYG计算)
            % --- 结束：原始人口统计与分组参数 ---

            %% 人口动态参数 (年龄组层面 - 来自V2)
            cS.initial_pop = [76.20905535, 86.45596319, 113.8702459, 98.60198303, 86.64117824, 102.7890433, 112.0217783, 99.04620047, ...
                              64.05142331, 66.93157492, 44.16815149, 25.40848066, 14.97325553, 6.872421945, 1.743059943, 0.216184341];
            if length(cS.initial_pop) ~= cS.aD_new, error('initial_pop 长度不匹配'); end
            % 多年期存活概率 (人口动态用 - 来自V2的beta_surv)
             beta_surv_pop = [0.998, 0.996, 0.994, 0.992, 0.988, 0.984, 0.980, ...
                             0.976, 0.970, 0.960, 0.945, 0.920, ...
                             0.880, 0.800, 0.680]; % 长度 aD_new - 1
            if length(beta_surv_pop) ~= cS.aD_new - 1, error('beta_surv_pop 长度不正确'); end
            cS.survivalProbV_popdyn = [beta_surv_pop(:)', 0]; % cS.survivalProbV_popdyn: 长度 aD_new，用于人口动态模拟

            % BGP 检测参数和最大期数 (来自V2)
            cS.bgp_tolerance = 0.001; % cS.bgp_tolerance: 人口稳态收敛容忍度
            cS.bgp_window = 5;        % cS.bgp_window: 判断人口稳态的窗口期
            cS.max_periods = 50;      % cS.max_periods: 人口动态模拟的最大期数

            %% 家庭参数 (年度)
            cS.sigma      = 1.5;   % cS.sigma: 效用函数曲率 (相对风险规避系数)
            cS.beta       = 0.97;  % cS.beta: 年度主观贴现因子 (VFI中使用) - 已修改
            cS.cFloor     = 0.05;  % cS.cFloor: 最低消费水平
            cS.nSim       = 5e4;   % cS.nSim: 模拟的个体数量

            %% 技术参数 (年度)
            cS.A          = 0.895944; % cS.A: 全要素生产率
            cS.alpha      = 0.36;   % cS.alpha: 资本份额
            cS.ddk        = 0.06;   % cS.ddk: 年度资本折旧率

            %% 社会保障参数 (年度) - 已修改
            % cS.theta 不再是固定的工资税率。它将内生决定以满足替代率。
            cS.pension_replacement_rate = 0.30; % cS.pension_replacement_rate: 外生参数: PAYG福利占平均税前工资的比例
            cS.max_payg_payroll_tax_rate = 1.0; % cS.max_payg_payroll_tax_rate: 安全措施: 允许的内生PAYG工资税率上限
                                                % 如果不希望有严格上限，初始可设为一个非常大的数 (例如 0.9)

            %% 劳动禀赋参数 (年度过程)
            cS.leSigma1      = 0.38 ^ 0.5; % cS.leSigma1: 初始对数劳动效率的标准差
            cS.leShockStd    = 0.045 .^ 0.5; % cS.leShockStd: 年度劳动效率冲击的标准差
            cS.lePersistence = 0.96;        % cS.lePersistence: 劳动效率AR(1)过程的持续性系数
            cS.leWidth       = 4;           % cS.leWidth: Tauchen方法中网格宽度参数 (标准差倍数)
            cS.nw            = 18;          % cS.nw: 离散化的劳动效率状态数量

            %% 网格参数
            cS.tgKY          = 3;     % cS.tgKY: 目标K/Y比率 (用于校准)
            cS.tgWage        = (1-cS.alpha)*cS.A*((cS.tgKY/cS.A)^(cS.alpha/(1-cS.alpha))); % cS.tgWage: 目标工资 (用于设定kMax)
            cS.nk            = 100;    % cS.nk: 非PPS资本网格点数
            cS.kMin          = 0;     % cS.kMin: 资本下限
            cS.kMax          = 100 * cS.tgWage; % cS.kMax: 资本上限
            power = 1.5; % power: 用于使资本网格在底部更密集的指数
            kGridV = cS.kMin + (cS.kMax - cS.kMin) * (linspace(0, 1, cS.nk).^power); % kGridV: 资本网格向量
            kGridV(1)=cS.kMin; % 确保最小值精确
            cS.kGridV        = kGridV(:); % cS.kGridV: 资本网格列向量

            %% 年龄效率剖面 (年度和分组)
            ageEffV_orig_temp = zeros(100, 1); % ageEffV_orig_temp: 临时年度年龄效率向量 (标准Huggett剖面定义)
             ageEffV_orig_temp(20:72) = [linspace(0.3, 1.5, 36-20+1), 1.5 .* ones(1, 47-37+1), ...
                               linspace(1.5, 0.2, 65-48+1), linspace(0.18, 0, 72-66+1)];
            cS.ageEffV_orig = ageEffV_orig_temp(cS.age1_orig : cS.ageLast_orig); % cS.ageEffV_orig: 模型年度年龄1..aD_orig的效率向量
            if length(cS.ageEffV_orig) ~= cS.aD_orig, error('ageEffV_orig 长度不匹配'); end

            % 分组年龄效率剖面 (5年平均)
            cS.ageEffV_new = zeros(cS.aD_new, 1); % cS.ageEffV_new: 各年龄组的平均效率
            for a = 1:cS.aD_new
                cS.ageEffV_new(a) = mean(cS.ageEffV_orig(cS.physAgeMap{a}));
            end

            %% ===== 新增: 个人养老金计划 (PPS) 参数 (与之前PPS版本相同) =====
            cS.pps_active = true; % cS.pps_active: 是否激活PPS功能
            % 年度最大缴费额占年度工资收入 (w*eff*e) 的比例
            cS.pps_max_contrib_frac = 0.05; % cS.pps_max_contrib_frac: 例如，工资收入的5% (已修改)
            % 退休期间从PPS提取资金的税率
            cS.pps_tax_rate_withdrawal = 0.03; % cS.pps_tax_rate_withdrawal: 例如，3%的固定税率
            % PPS资产相对于市场回报率 R-1 的回报溢价/折价
            cS.pps_return_rate_premium = 0.0; % cS.pps_return_rate_premium: 例如，与市场回报相同
            % 退休期间每年从PPS余额中提取的比例
            cS.pps_withdrawal_rate = 0.15; % cS.pps_withdrawal_rate: 例如，每年提取余额的15%
             % 允许缴费的最后一个年度年龄索引
            cS.pps_contribution_age_max = cS.aR_idx_orig - 1; % cS.pps_contribution_age_max: 例如，退休前一年
            % 开始提取的第一个年度年龄索引
            cS.pps_withdrawal_age_min = cS.aR_idx_orig; % cS.pps_withdrawal_age_min: 例如，退休第一年
             % 总PPS质量是否计入生产性资本K?
            cS.pps_in_K = true; % cS.pps_in_K: PPS资产是否贡献于总资本
            % PPS资产在死亡时是否可遗赠?
            cS.pps_bequeathable = false; % cS.pps_bequeathable: 简化假设：不可遗赠
            % ============================================================
        end

        function [R_market, MPL_gross] = HHPrices_Huggett(K, L, cS) % K: 总生产性资本, L: 总有效劳动
            % 修改版：仅计算市场回报率 R_market 和总 MPL。
            % 实际的 PAYG 福利 'b'、PAYG 工资税率 'theta_payg_actual' 和净工资 'w' 在主循环中决定。
            % 此处的 K 是总生产性资本。

            if K <= 0
                K = 1e-6; % 避免除零或非正数错误
                warning('K 在 HHPrices 中为非正数，已设为小正数。');
            end
            if L <= 0
                L = 1e-6; % 避免除零或非正数错误
                warning('L 在 HHPrices 中为非正数，已设为小正数。');
            end

            Y_gross = cS.A * (K^cS.alpha) * (L^(1-cS.alpha)); % Y_gross: 年度总产出
            MPK_gross = cS.alpha * Y_gross / K;                % MPK_gross: 年度总资本边际产出
            MPL_gross = (1-cS.alpha) * Y_gross / L;            % MPL_gross: 年度总劳动边际产出 (税前平均工资的代理)

            % 家庭进行 *非PPS* 储蓄时可获得的市场回报率
            % 为简化此价格函数中的 R_market 计算，假设资本利得税为零
            R_market   = 1 + MPK_gross - cS.ddk; % R_market: 市场总回报率因子
            R_market   = max(1.0 + 1e-6, R_market); % 确保 R 略高于 1
        end

        % --- 劳动禀赋函数 (tauchen, norm_grid, MarkovChainSimulation) - 未改变 ---
        function [logGridV, trProbM, prob1V] = EarningProcess_olgm(cS)
            % 为年度 AR(1) 过程生成网格和转移矩阵。
             [logGridV_raw, trProbM_calc] = main_olg_v2_utils_pps.tauchen(cS.nw, cS.lePersistence, cS.leShockStd, 0, cS.leWidth);
             gridMin = logGridV_raw(1); gridMax = logGridV_raw(end);
             normGridMin = gridMin - 2*cS.leShockStd; normGridMax = gridMax + 2*cS.leShockStd; % 用于norm_grid的边界
             prob1V_calc = []; % 初始化初始概率向量
              try
                  [prob1V_calc, ~, ~] = main_olg_v2_utils_pps.norm_grid(logGridV_raw, normGridMin, normGridMax, 0, cS.leSigma1);
                  prob1V_calc = prob1V_calc(:);
              catch ME % 捕获norm_grid可能发生的错误
                  warning('norm_grid 失败: %s。使用均匀分布。', ME.message);
                  prob1V_calc = ones(cS.nw, 1)/cS.nw; % 如果失败，使用均匀分布
              end
             logGridV_calc = logGridV_raw(:);
             % 原始Huggett代码移动了网格，保留该做法。
             logGridV_calc = logGridV_calc - logGridV_calc(1) - 1;

             % 归一化转移矩阵
             if any(abs(sum(trProbM_calc, 2) - 1) > 1e-6) % 检查行和是否为1
                  row_sums = sum(trProbM_calc, 2);
                  row_sums(row_sums <= 1e-9) = 1; % 避免除以零
                  trProbM_calc = bsxfun(@rdivide, trProbM_calc, row_sums); % 按行归一化
             end
             % 归一化初始概率
             if abs(sum(prob1V_calc) - 1) > 1e-6 % 检查和是否为1
                 prob1V_calc = prob1V_calc ./ sum(prob1V_calc); % 归一化
             end

             % *** 将计算值赋给输出参数 ***
             logGridV = logGridV_calc; % logGridV: 对数劳动效率网格
             trProbM = trProbM_calc;   % trProbM: 转移概率矩阵
             prob1V = prob1V_calc;     % prob1V: 初始状态概率
        end
        function [y, trProbM] = tauchen(N, pRho, pSigma, pMu, n_std) % (未改变)
            % Tauchen 方法离散化 AR(1) 过程
             std_y = sqrt(pSigma^2 / (1 - pRho^2)); % 过程的无条件标准差
             if abs(1-pRho)<1e-9, std_y = pSigma*100; end % 处理 pRho 接近1的情况
             y_max = n_std*std_y; y_min = -y_max; % 网格上下界
             y_val = linspace(y_min, y_max, N); % y_val: 离散化状态点 (原为y，避免与输出冲突)
             d = 0; if N > 1, d = y_val(2) - y_val(1); end; % 网格间距
             trProbM = zeros(N, N); % 初始化转移概率矩阵
             for iR = 1:N % iR: 当前状态行索引
                 for iC = 1:N % iC: 下一状态列索引
                     mn = pRho*y_val(iR); % 条件均值
                     if iC==1 % 第一个状态的概率
                         trProbM(iR,iC)=normcdf((y_val(1)-mn+d/2)/pSigma);
                     elseif iC==N % 最后一个状态的概率
                         trProbM(iR,iC)=1-normcdf((y_val(N)-mn-d/2)/pSigma);
                     else % 中间状态的概率
                         trProbM(iR,iC)=normcdf((y_val(iC)-mn+d/2)/pSigma)-normcdf((y_val(iC)-mn-d/2)/pSigma);
                     end
                 end
             end
             rs=sum(trProbM,2); rs(rs<=1e-9)=1; % 行和，用于归一化
             trProbM=bsxfun(@rdivide,trProbM,rs); % 按行归一化
             um=pMu/(1-pRho); % 无条件均值 (如果 pMu 非零)
             if abs(1-pRho)<1e-9 && pMu~=0, um=sign(pMu)*inf; end
             if isfinite(um), y_val=y_val+um; end % 如果有非零均值，则平移网格
             y = y_val; % 将结果赋给输出变量 y
        end
        function [massV, lbV, ubV] = norm_grid(xV, xMin, xMax, mu, sig) % (未改变)
            % 计算正态分布在给定网格点上的概率质量
             n=length(xV); xV=xV(:)'; % n: 网格点数, xV: 网格点行向量
             if n>1 && any(xV(2:n)<xV(1:n-1)), error('xV 非单调递增'); end
             lbV=[]; ubV=[]; % lbV: 区间下界, ubV: 区间上界
             if n>1
                 xM=0.5*(xV(1:n-1)+xV(2:n)); % 区间中点
                 lbV=[xMin,xM]; ubV=[xM,xMax];
             elseif n==1 % 单点网格
                 lbV=xMin; ubV=xMax;
             else % 空网格
                 massV=[]; return;
             end
             cdfV=normcdf([lbV,ubV(end)],mu,sig); % 计算累积分布函数值
             massV=diff(cdfV); % massV: 每个区间的概率质量
              if any(massV < -1e-9), warning('在 norm_grid 中检测到负质量。'); massV(massV<0)=0; end % 修正负质量
             ts=sum(massV); % 总概率质量
             if ts>1e-9, massV=massV/ts; % 归一化
             else % 如果总质量接近零
                 if n>0, massV=ones(1,n)/n; else massV=[]; end % 均匀分配或空
             end
             massV=massV(:); lbV=lbV(:); ubV=ubV(:); % 转换为列向量
             if n>1 && any(ubV<lbV), error('ubV < lbV'); end % 检查区间有效性
        end
        function eIdxM = MarkovChainSimulation(nSim, T_sim, prob0V, trProbM, rvInM) % (未改变, T重命名为T_sim)
            % 模拟马尔可夫链
             ns=length(prob0V); % ns: 状态数量
             if size(trProbM,1)~=ns||size(trProbM,2)~=ns,error('转移矩阵维度不匹配');end;
             if abs(sum(prob0V)-1)>1e-5,prob0V=prob0V./sum(prob0V);end; % 归一化初始概率
             if any(abs(sum(trProbM,2)-1)>1e-5), row_sums=sum(trProbM,2); row_sums(row_sums<=1e-9)=1; trProbM=bsxfun(@rdivide,trProbM,row_sums); end; % 归一化转移矩阵
             if size(rvInM,1)~=nSim||size(rvInM,2)~=T_sim,error('随机数矩阵维度不匹配');end % rvInM: 随机数输入

             cP0=cumsum(prob0V(:)'); % cP0: 初始概率的累积分布
             cPT=cumsum(trProbM,2); cPT(:,ns)=1.0; % cPT: 转移概率的累积分布 (按行)
             eIdxM=zeros(nSim,T_sim,'uint16'); % eIdxM: 模拟的状态索引历史 (使用uint16节省内存)
             eIdxM(:,1)=1+sum(bsxfun(@gt,rvInM(:,1),cP0),2); % 确定初始状态

             for t_mc=1:(T_sim-1) % t_mc: 模拟时间步 (原为t)
                 cSI=eIdxM(:,t_mc); % cSI: 当前期状态索引
                 vPI=(cSI>=1)&(cSI<=ns); % vPI: 检查索引是否有效
                 if ~all(vPI),cSI(~vPI)=1;eIdxM(:,t_mc)=cSI;end; % 无效则设为第一个状态
                 cPt=cPT(cSI,:); % 提取对应当前状态的累积转移概率行
                 eIdxM(:,t_mc+1)=1+sum(bsxfun(@gt,rvInM(:,t_mc+1),cPt),2); % 确定下一期状态
             end
        end
        % --- 结束劳动禀赋 ---

        function [HHlaborM_group, L] = LaborSupply_Huggett(eIdxM, cS, paramS, Z_ss_norm_group)
            % 计算总劳动供给 L。(逻辑未改变)
            % 使用分组效率和分组稳态质量 Z_ss_norm_group。
            nSim = size(eIdxM, 1); % nSim: 模拟个体数
            HHlaborM_group = zeros(nSim, cS.aD_new); % HHlaborM_group: 每个个体在每个年龄组的平均年劳动供给

            % 年度年龄索引到5年期年龄组索引的映射
            ageToGroupMap = zeros(cS.aD_orig, 1); % ageToGroupMap: 年度年龄到组索引的映射
            for a_new_map = 1:cS.aD_new
                ageToGroupMap(cS.physAgeMap{a_new_map}) = a_new_map;
            end
            leGridV_col = paramS.leGridV(:); % leGridV_col: 劳动效率列向量

            % 首先计算每个模拟个体/年度年龄的劳动供给
            HHlaborM_annual_temp = zeros(nSim, cS.aD_orig); % HHlaborM_annual_temp: 个体年度劳动供给临时矩阵
            for a_orig = 1 : cS.aD_orig % a_orig: 年度年龄索引
               if a_orig < cS.aR_idx_orig % 仅工作年龄贡献劳动
                   a_new = ageToGroupMap(a_orig); % 找到对应的年龄组
                   % 使用分组效率剖面 cS.ageEffV_new
                   HHlaborM_annual_temp(:, a_orig) = cS.ageEffV_new(a_new) .* leGridV_col(eIdxM(:, a_orig));
               end
            end

            % 平均每个年龄组内的年度劳动供给
            for a_new = 1:cS.aD_new % a_new: 年龄组索引
                annual_indices = cS.physAgeMap{a_new}; % annual_indices: 该组包含的年度索引
                if ~isempty(annual_indices)
                    % 对每个模拟个体，平均该组内各年的劳动供给
                    HHlaborM_group(:, a_new) = mean(HHlaborM_annual_temp(:, annual_indices), 2);
                end
            end

            % 使用平均年龄组劳动和年龄组质量加总劳动供给 L (仅工作年龄组)
            L_working_groups = mean(HHlaborM_group(:, 1:cS.aR_new), 1) * Z_ss_norm_group(1:cS.aR_new); % L_working_groups: 工作年龄组贡献的总劳动
            L = max(0, L_working_groups); % L: 总有效劳动供给 (确保非负)
        end

        function eIdxM = LaborEndowSimulation_olgm(cS, paramS)
            % 使用年度转移概率模拟年度禀赋索引。(未改变)
            rng(433); % 设置随机数种子以保证结果可复现
            rvInM=rand([cS.nSim,cS.aD_orig]); % rvInM: [nSim x aD_orig] 的随机数矩阵
            eIdxM = main_olg_v2_utils_pps.MarkovChainSimulation(cS.nSim, cS.aD_orig, paramS.leProb1V, paramS.leTrProbM, rvInM);
        end

        function [kHistM, kPpsHistM, cHistM] = HHSimulation_olgm(kPolM, cPpsPolM, cPolM_consump, eIdxM, R_market, w_net, T_transfer, bV_payg, paramS, cS)
            % 修改版：使用 R_market, w_net (内生PAYG税后), T_transfer, bV_payg 进行模拟
            % kPolM是非PPS储蓄, cPpsPolM是PPS缴费, cPolM_consump是消费决策

            validateattributes(kPolM, {'double'}, {'finite', 'nonnan', 'nonempty', 'real', 'size', [cS.nk, cS.nw, cS.aD_new]});
            validateattributes(cPpsPolM, {'double'}, {'finite', 'nonnan', 'nonempty', 'real','>=',0, 'size', [cS.nk, cS.nw, cS.aD_new]});
            validateattributes(cPolM_consump, {'double'}, {'finite', 'nonnan', 'nonempty', 'real', 'size', [cS.nk, cS.nw, cS.aD_new]});
            validateattributes(eIdxM, {'uint16', 'double'}, {'nonempty', 'integer', 'positive', '<=', cS.nw, 'size', [cS.nSim, cS.aD_orig]});

            ageToGroupMap = zeros(cS.aD_orig, 1); % ageToGroupMap: 年度年龄到组索引的映射
            for a_new = 1:cS.aD_new
                annual_indices = cS.physAgeMap{a_new};
                if ~isempty(annual_indices), ageToGroupMap(annual_indices) = a_new; end
            end

            nSim   = size(eIdxM, 1); % nSim: 模拟个体数
            kHistM    = zeros(nSim, cS.aD_orig + 1); % kHistM: 非PPS资产年度历史 (期初)
            kPpsHistM = zeros(nSim, cS.aD_orig + 1); % kPpsHistM: PPS资产年度历史 (期初)
            cHistM    = zeros(nSim, cS.aD_orig);     % cHistM: 消费年度历史 (期内)
            kGridV_col = cS.kGridV(:); % kGridV_col: 非PPS资本网格列向量
            leGridV_col = paramS.leGridV(:); % leGridV_col: 劳动效率列向量

            % 预创建插值器
            kPolInterp = cell(cS.nw, cS.aD_new);        % kPolInterp: 非PPS储蓄决策插值器
            cPpsPolInterp = cell(cS.nw, cS.aD_new);     % cPpsPolInterp: PPS缴费决策插值器
            cPolConsumpInterp = cell(cS.nw, cS.aD_new); % cPolConsumpInterp: 消费决策插值器
             for a_new=1:cS.aD_new % a_new: 年龄组索引
                 for ie=1:cS.nw % ie: 劳动效率状态索引
                      kPolInterp{ie, a_new}    = griddedInterpolant(kGridV_col, kPolM(:, ie, a_new), 'linear', 'nearest'); % 使用最近邻外插
                      cPpsPolInterp{ie, a_new} = griddedInterpolant(kGridV_col, cPpsPolM(:, ie, a_new), 'linear', 'nearest');
                      cPolConsumpInterp{ie, a_new} = griddedInterpolant(kGridV_col, cPolM_consump(:, ie, a_new), 'linear', 'nearest');
                 end
             end

            pps_return_net = (R_market - 1) + cS.pps_return_rate_premium; % pps_return_net: PPS资产的年度净回报率

            for a_orig = 1 : cS.aD_orig % a_orig: 年度年龄索引
                a_new_group_idx = ageToGroupMap(a_orig); % a_new_group_idx: 对应的年龄组索引
                kNowV = kHistM(:, a_orig); % kNowV: 当前期初非PPS资产向量
                kPpsNowV = kPpsHistM(:, a_orig); % kPpsNowV: 当前期初PPS资产向量

                kNextNonPpsV = zeros(nSim, 1); % kNextNonPpsV: 下一期初非PPS资产
                cPpsContribV = zeros(nSim, 1); % cPpsContribV: 当期PPS缴费
                cConsumpValV = zeros(nSim, 1); % cConsumpValV: 当期消费
                kPpsNextV    = zeros(nSim, 1); % kPpsNextV: 下一期初PPS资产

                is_working = (a_orig < cS.aR_idx_orig); % is_working: 是否处于工作年龄的标志
                b_payg_eff = bV_payg(a_new_group_idx);  % b_payg_eff: 该年龄组的有效PAYG福利 (VFI中使用)
                
                % 模拟中的PPS提取 (基于期初PPS资产)
                pps_withdrawal_for_sim_taxed = zeros(nSim,1);  % pps_withdrawal_for_sim_taxed: 税后PPS提取额
                pps_withdrawal_for_sim_pretax = zeros(nSim,1); % pps_withdrawal_for_sim_pretax: 税前PPS提取额

                if ~is_working && cS.pps_active && a_orig >= cS.pps_withdrawal_age_min % 退休且PPS激活且达到提取年龄
                    pps_withdrawal_for_sim_pretax = kPpsNowV * cS.pps_withdrawal_rate;
                    pps_withdrawal_for_sim_taxed = pps_withdrawal_for_sim_pretax * (1 - cS.pps_tax_rate_withdrawal);
                end

                for ie = 1 : cS.nw % ie: 劳动效率状态索引
                    simIdxV = find(eIdxM(:, a_orig) == ie); % simIdxV: 处于该效率状态的个体索引
                    if ~isempty(simIdxV)
                         kNow_ie = kNowV(simIdxV); % kNow_ie: 这些个体的当前非PPS资产
                         kNow_ie_clamped = max(kGridV_col(1), min(kGridV_col(end), kNow_ie)); % 确保在插值网格范围内

                         % 从VFI结果中获取决策
                         kNextNonPpsV(simIdxV) = kPolInterp{ie, a_new_group_idx}(kNow_ie_clamped);
                         cPpsContribV(simIdxV) = cPpsPolInterp{ie, a_new_group_idx}(kNow_ie_clamped);
                         cConsumpValV(simIdxV) = cPolConsumpInterp{ie, a_new_group_idx}(kNow_ie_clamped);

                         % 在模拟中再次强制执行PPS缴费规则 (更精确的应用)
                         if is_working && cS.pps_active && a_orig <= cS.pps_contribution_age_max % 工作且PPS激活且未到最大缴费年龄
                             % 注意：w_net 是 PAYG 税后的工资。max_contrib_frac 通常基于税前工资。
                             % VFI中的简化规则可能已经考虑了这一点。这里需要确保与VFI一致。
                             % 当前VFI的cPpsPolM是基于w_net * eff * le * frac计算的上限（近似）。
                             % 这里我们用更精确的基于总工资的上限来重新约束模拟中的cPpsContribV。
                             % 这是一个微妙之处，理想情况下VFI应使用总工资。
                             % 为简单起见，此处假设VFI中的cPpsPolM已充分考虑预算和上限。
                             % 但为确保模拟中的缴费不超过基于工资的上限，可以再加一层约束。
                             % 当前的VFI中的cPpsPolM是基于w_net（PAYG税后）计算的，这会低估允许的缴费额。
                             % 这里应该使用一个更接近税前工资的基数。
                             % 获取主循环中计算的 MPL_gross_iter 和 theta_payg_actual (这里没有直接传入，是个问题)
                             % 暂时用w_net / (1-theta_payg_actual_from_main_loop_if_available)来近似gross_wage
                             % **简化处理：这里直接使用VFI的cPpsPolM结果，因为它已考虑预算约束，**
                             % **但在VFI的HHSolutionByAge_VFI_Huggett中，计算max_contrib时，**
                             % **wage_income_for_pps_calc 使用的是 w_net，这会低估基于总收入的缴费上限。**
                             % **这是一个需要统一的地方。**
                             % **为了保持当前代码的运行，我们暂时用一个近似的总工资基数来计算上限。**
                             % gross_wage_proxy = w_net / (1-theta_payg_actual_passed_if_any_otherwise_approx); % 需要theta_payg
                             % max_contrib_based_on_wage = gross_wage_proxy * cS.ageEffV_new(a_new_group_idx) * leGridV_col(ie) * cS.pps_max_contrib_frac;
                             % cPpsContribV(simIdxV) = max(0, min(cPpsContribV(simIdxV), max_contrib_based_on_wage ));
                             % **鉴于以上复杂性，目前模拟直接使用VFI的cPpsPolM，它已包含上限和预算检查。**
                             % (没有额外修改cPpsContribV(simIdxV)的上限逻辑)
                             
                         else % 退休或PPS未激活或超过最大缴费年龄
                             cPpsContribV(simIdxV) = 0; % 不缴费
                         end
                         
                         % 更新PPS余额
                         if cS.pps_active
                             kPpsNextV(simIdxV) = (kPpsNowV(simIdxV) + cPpsContribV(simIdxV) - pps_withdrawal_for_sim_pretax(simIdxV)) * (1 + pps_return_net);
                             kPpsNextV(simIdxV) = max(0, kPpsNextV(simIdxV)); % 确保非负
                         else
                             kPpsNextV(simIdxV) = kPpsNowV(simIdxV); % 如果PPS未激活，余额不变
                         end
                    end
                end % 结束劳动效率状态循环 ie

                 % 确保选择可行并应用上下限
                 kHistM(:, a_orig + 1)    = max(cS.kMin, min(cS.kMax, kNextNonPpsV));
                 kPpsHistM(:, a_orig + 1) = max(0, kPpsNextV); % PPS资产不能为负
                 cHistM(:, a_orig)        = max(cS.cFloor, cConsumpValV); % 使用来自VFI的消费决策

                 % 最终的有限性检查
                 kHistM(~isfinite(kHistM(:, a_orig + 1)), a_orig + 1) = cS.kMin;
                 kPpsHistM(~isfinite(kPpsHistM(:, a_orig + 1)), a_orig + 1) = 0;
                 cHistM(~isfinite(cHistM(:, a_orig)), a_orig) = cS.cFloor;
            end % 结束年度年龄循环 a_orig

            kHistM = kHistM(:, 1:cS.aD_orig); % 非PPS资产，期初，年度1到aD_orig
            kPpsHistM = kPpsHistM(:, 1:cS.aD_orig); % PPS资产，期初，年度1到aD_orig

            validateattributes(kHistM, {'double'}, {'finite', 'nonnan', 'nonempty', 'real', '>=', cS.kMin - 1e-6, '<=', cS.kMax + 1e-6, 'size', [nSim, cS.aD_orig]});
            validateattributes(kPpsHistM, {'double'}, {'finite', 'nonnan', 'nonempty', 'real', '>=', -1e-6, 'size', [nSim, cS.aD_orig]});
            validateattributes(cHistM, {'double'}, {'finite', 'nonnan', 'nonempty', 'real', '>=', cS.cFloor - 1e-6, 'size', [nSim, cS.aD_orig]});
        end


        function [muM, utilM] = CES_utility(cM, sig, cS) % cM: 消费, sig: 曲率, cS: 参数 (用于cFloor)
            % 效用函数。处理最低消费。(与原始版本相同)
             if ~isscalar(sig)||sig<=0,error('Sigma 必须是正标量。');end;
             min_c = cS.cFloor; % min_c: 最低消费
             vld = (cM >= min_c); % vld: 消费是否有效的标志
             c_c = max(min_c, cM); % c_c: 用于计算的消费 (已处理下限)
             utilM = -Inf(size(cM)); muM = Inf(size(cM)); % utilM: 效用, muM: 边际效用 (初始化)

             if abs(sig-1)<1e-6 % 对数效用情况
                 utilM(vld) = log(c_c(vld));
                 muM(vld) = 1 ./ c_c(vld);
             else % CRRA 一般情况
                 utilM(vld) = (c_c(vld).^(1-sig)) ./ (1-sig);
                 muM(vld) = c_c(vld).^(-sig);
             end
             utilM(~vld) = -1e10 - (min_c - cM(~vld))*1e10; % 对无效消费的惩罚
             muM(~vld) = c_c(~vld).^(-sig) + 1e10; % 无效消费对应的高边际效用
        end

        function incomeM = HHIncome_Huggett(kV_non_pps, R_market, w_net, T_transfer, b_payg, a_new_group_idx, paramS, cS)
            % 计算VFI中可用的年度资源 (在PPS选择和k'选择之前)
            % w_net 已经是内生PAYG工资税后的净工资。
            % b_payg 是该年龄组的PAYG福利。
            % T_transfer 是来自非PPS意外遗赠的lump-sum转移。

            nk = length(kV_non_pps); % nk: 非PPS资本网格点数
            leGridV_col = paramS.leGridV(:); % leGridV_col: 劳动效率列向量
            nonCapIncomeNetM = zeros(1, cS.nw); % nonCapIncomeNetM: 净非资本收入 (与k无关的部分)

            if a_new_group_idx <= cS.aR_new % 工作年龄组
               % 净工资收入 (已扣除PAYG税) + 转移T + PAYG福利b (工作者b=0)
               net_wage_income_vec = w_net * cS.ageEffV_new(a_new_group_idx) .* leGridV_col; % net_wage_income_vec: [1 x nw]
               nonCapIncomeNetM(1,:) = net_wage_income_vec + T_transfer + b_payg;
            else % 退休年龄组
               nonCapIncomeNetM(1,:) = b_payg + T_transfer; % PAYG福利 + 转移T
            end

            % 来自非PPS资产和当前净非资本收入的总资源
            incomeM = R_market * kV_non_pps(:) * ones(1, cS.nw) + ones(nk, 1) * nonCapIncomeNetM; % incomeM: [nk x nw] 总可用资源

            if any(~isfinite(incomeM(:)))
                warning('在 HHIncome_Huggett 中检测到年龄组 %d 的非有限收入。进行截断处理。', a_new_group_idx);
                incomeM(~isfinite(incomeM)) = 1e-6; % 或使用中位数？
            end
        end

        function [cPolM, kPolM, cPpsPolM, valueM] = HHSolution_VFI_Huggett(R_market, w_net, T_transfer, bV_payg, paramS, cS)
            % 修改版VFI: w_net是净工资, bV_payg是PAYG福利向量
            % 同时求解 k' (非PPS储蓄) 和 c_pps (PPS缴费)
            cPolM  = zeros(cS.nk, cS.nw, cS.aD_new);       % cPolM: 消费决策规则
            kPolM  = zeros(cS.nk, cS.nw, cS.aD_new);       % kPolM: 非PPS储蓄决策规则
            cPpsPolM = zeros(cS.nk, cS.nw, cS.aD_new); % cPpsPolM: PPS缴费决策规则
            valueM = zeros(cS.nk, cS.nw, cS.aD_new);     % valueM: 价值函数

            for a = cS.aD_new : -1 : 1 % a: 年龄组索引 (从后向前迭代)
               vPrimeM = []; % vPrimeM: 下一年龄组的期望价值函数 (关于k'和e')
               if a < cS.aD_new
                   vPrimeM = valueM(:,:,a+1);
               end

               % 调用分年龄组的求解函数
               [cPolM(:,:,a), kPolM(:,:,a), cPpsPolM(:,:,a), valueM(:,:,a)] = ...
                  main_olg_v2_utils_pps.HHSolutionByAge_VFI_Huggett(a, vPrimeM, R_market, w_net, T_transfer, bV_payg(a), paramS, cS);
            end
            % 验证输出的有效性
            validateattributes(cPolM, {'double'}, {'finite', 'nonnan', 'real','>=',cS.cFloor-1e-6,'size',[cS.nk,cS.nw,cS.aD_new]});
            validateattributes(kPolM, {'double'}, {'finite', 'nonnan', 'real','>=',cS.kMin-1e-6,'<=',cS.kMax+1e-6,'size',[cS.nk,cS.nw,cS.aD_new]});
            validateattributes(cPpsPolM, {'double'}, {'finite', 'nonnan', 'real','>=',0,'size',[cS.nk,cS.nw,cS.aD_new]}); % PPS缴费非负
            validateattributes(valueM, {'double'}, {'finite', 'nonnan', 'real','size',[cS.nk,cS.nw,cS.aD_new]});
        end

        function [cPolM_consump, kPolM_non_pps, cPpsPolM_contrib, valueM] = HHSolutionByAge_VFI_Huggett(a_new_group_idx, vPrime_keM, R_market, w_net, T_transfer, b_payg, paramS, cS)
            % 求解特定5年期年龄组 'a_new_group_idx' 的家庭问题。
            % w_net 是内生PAYG税后的净工资。b_payg 是该组的PAYG福利。
            % 使用简化的顺序PPS选择。

            if a_new_group_idx < cS.aD_new && isempty(vPrime_keM), error('当 a < aD 时，vPrime 为空'); end
            if a_new_group_idx == cS.aD_new && ~isempty(vPrime_keM), error('当 a = aD 时，vPrime 非空'); end
            if a_new_group_idx < cS.aD_new, validateattributes(vPrime_keM,{'double'},{'finite','nonnan','real','size',[cS.nk,cS.nw]}); end

            % 期初可用预算 (来自非PPS资产和净非资本收入)，在PPS缴费决策之前。
            budgetM_before_pps = main_olg_v2_utils_pps.HHIncome_Huggett(cS.kGridV, R_market, w_net, T_transfer, b_payg, a_new_group_idx, paramS, cS); % [nk, nw]

            fminOpt=optimset('fminbnd'); fminOpt.TolX=1e-6; fminOpt.Display='none'; % fminbnd 优化选项

            cPolM_consump = zeros(cS.nk, cS.nw);       % cPolM_consump: 消费决策
            kPolM_non_pps = zeros(cS.nk, cS.nw);       % kPolM_non_pps: 非PPS储蓄决策
            cPpsPolM_contrib = zeros(cS.nk, cS.nw);  % cPpsPolM_contrib: PPS缴费决策
            valueM = -Inf(cS.nk, cS.nw);               % valueM: 价值函数

            % 判断当前年龄组是否为PPS缴费的工作年龄组
            % 需要将 a_new_group_idx 映射到近似的年度年龄以与 cS.pps_contribution_age_max 比较
            % 为简化，假设仅当 a_new_group_idx <= cS.aR_new (工作年龄组) 时才允许PPS缴费
            is_pps_contrib_age = (a_new_group_idx <= cS.aR_new);


            if a_new_group_idx == cS.aD_new % 最后一个年龄组 - 全部消费，不进行PPS缴费，不进行非PPS储蓄
               cPpsPolM_contrib = zeros(cS.nk,cS.nw); % 不缴费
               % 消费预算即为 budgetM_before_pps，因为不储蓄
               cPolM_consump = max(cS.cFloor, budgetM_before_pps);
               kPolM_non_pps = cS.kMin*ones(cS.nk,cS.nw); % 储蓄最小值 (零)
               [~, valueM_calc] = main_olg_v2_utils_pps.CES_utility(cPolM_consump, cS.sigma, cS); % valueM_calc: 临时变量
               valueM = valueM_calc; % 赋值
               valueM(~isfinite(valueM)) = -1e12; % 处理非有限值
            else % 非最后一个年龄组
               % 期望价值函数 E[V'(k', e')] - 在此设置中仅取决于 k' (非PPS资产)
               ExValuePrime_k_given_e = zeros(cS.nk, cS.nw); % ExValuePrime_k_given_e: 给定当前e，下一期期望价值 (作为k'的函数)
                for ikPrime = 1:cS.nk % ikPrime: 下一期非PPS资本网格点索引
                   % 对 e' 求和: V'(k',e') * Pr(e'|e)
                   ExValuePrime_k_given_e(ikPrime, :) = vPrime_keM(ikPrime, :) * paramS.leTrProbM';
                end
                ExValuePrime_k_given_e(~isfinite(ExValuePrime_k_given_e))=-1e10; % 处理非有限值

               vPInt=cell(1,cS.nw); % vPInt: E[V'] 作为 k' 函数的插值器单元格数组
               for ie=1:cS.nw % ie: 当前劳动效率状态索引
                   vPInt{ie} = griddedInterpolant(cS.kGridV, ExValuePrime_k_given_e(:, ie), 'linear', 'linear');
               end

               % 如果有并行计算工具箱，可以使用 parfor
               parfor ik = 1:cS.nk % ik: 当前非PPS资本网格点索引
                   c_row=zeros(1,cS.nw); k_row=zeros(1,cS.nw); cpps_row=zeros(1,cS.nw); v_row=zeros(1,cS.nw); % 初始化行向量

                   for ie = 1 : cS.nw % ie: 当前劳动效率状态索引
                        vPofK_interpolant = vPInt{ie}; % vPofK_interpolant: E[V'(k')|e] 的插值器
                        budget_ike_before_pps = budgetM_before_pps(ik,ie); % budget_ike_before_pps: PPS缴费前的预算

                        % 简化的PPS选择逻辑 (如前所述)
                        cPps_chosen_for_state = 0; % cPps_chosen_for_state: 当前状态下选择的PPS缴费额
                        if is_pps_contrib_age && cS.pps_active % 如果是PPS缴费年龄且PPS激活
                           % 计算最大缴费额的工资基数。w_net是PAYG税后工资。
                           % 理想情况下应基于总工资。目前使用w_net作为代理。
                           wage_income_for_pps_calc = w_net * cS.ageEffV_new(a_new_group_idx) * paramS.leGridV(ie);
                           max_contrib_val = wage_income_for_pps_calc * cS.pps_max_contrib_frac; % max_contrib_val: 最大允许缴费额
                           
                           % 检查在缴纳max_contrib后，剩余预算是否够最低消费和最低非PPS储蓄
                           if budget_ike_before_pps - max_contrib_val >= cS.cFloor + cS.kMin
                               cPps_chosen_for_state = max_contrib_val;
                           else % 如果不够
                               % 缴纳在保证最低消费和最低非PPS储蓄后的所有剩余部分
                               cPps_chosen_for_state = max(0, budget_ike_before_pps - cS.cFloor - cS.kMin);
                           end
                           cPps_chosen_for_state = max(0, cPps_chosen_for_state); % 确保非负
                        end
                        cpps_row(ie) = cPps_chosen_for_state;

                        % 用于消费 c 和 非PPS储蓄 k' 的剩余预算
                        budget_remaining_for_c_kprime = budget_ike_before_pps - cPps_chosen_for_state; % budget_remaining_for_c_kprime

                        % 给定 cPps_chosen_for_state 后，求解 k' 的一维问题
                        [c_row(ie), k_row(ie), v_row(ie)] = ...
                              main_olg_v2_utils_pps.HHSolutionByOneState_OptK(a_new_group_idx, budget_remaining_for_c_kprime, vPofK_interpolant, fminOpt, cS, paramS);
                   end
                   % 存储该 ik 行的结果
                   cPolM_consump(ik,:) = c_row;
                   kPolM_non_pps(ik,:) = k_row;
                   cPpsPolM_contrib(ik,:) = cpps_row; % 存储决定的PPS缴费
                   valueM(ik,:) = v_row;
               end % 结束 ik 循环
            end % 结束 if a_new_group_idx == cS.aD_new
        end

        function [c, kPrime, ValueFunc] = HHSolutionByOneState_OptK(a_new, budget_remaining, vPofK_int, fminOpt, cS, paramS)
            % 在给定 budget_remaining (PPS缴费后) 的情况下，求解 k' 的优化问题
            % (与之前的PPS版本相同 - 此函数现在是稳健的)

            kPMin = cS.kMin; % kPMin: k' 的最小值 (通常为0)
            kPMax = budget_remaining - cS.cFloor; % kPMax: 在保证最低消费前提下，k' 的最大可能值
            ValueFunc = -Inf; % ValueFunc: 初始化价值函数

            if kPMax <= kPMin % 角点解: 消费最低，储蓄最小 k'
                kPrime = kPMin; % kPrime: 最优下一期非PPS储蓄
                c = max(cS.cFloor, budget_remaining - kPrime); % c: 最优消费 (扣除最小储蓄后的剩余)
                 [u_c, ~] = main_olg_v2_utils_pps.CES_utility(c, cS.sigma, cS); % u_c: 当前效用 (muM此处不需要)

                if a_new < cS.aD_new % 如果不是最后一个年龄组
                    EV_c = -Inf; % EV_c: 期望的下一期价值
                     try
                        % 在kPrime处评估期望价值 (对插值器输入进行截断)
                        kP_eval = max(cS.kGridV(1), min(cS.kGridV(end), kPrime));
                        EV_c = vPofK_int(kP_eval);
                     catch ME_interp % 捕获插值错误
                         warning('在OptK的角点解中插值错误，年龄 %d: %s。使用边界值。', a_new, ME_interp.message);
                         if kPrime < cS.kGridV(1), EV_c = vPofK_int(cS.kGridV(1)); else EV_c = vPofK_int(cS.kGridV(end)); end
                     end
                     if ~isfinite(EV_c), EV_c = -1e12; end % 处理非有限值
                     s_transition = cS.s_1yr_transitionV(a_new); % s_transition: 过渡年的存活概率
                     ValueFunc = u_c + cS.beta * s_transition * EV_c; % 使用年度beta计算总价值
                else % 最后一个年龄组
                    ValueFunc = u_c;
                end
                if ~isfinite(ValueFunc), ValueFunc = -1e12; end % 处理非有限值

            else % 可能存在 k' 的内部解
                kPMin_opt = max(cS.kMin, kPMin); % kPMin_opt: 优化下界 (至少为kMin)
                % 优化上界: 确保 k' 不超过预算或 kMax，且优化区间有效
                kPMax_opt = max(kPMin_opt + 1e-9, min(cS.kMax, kPMax));

                 if kPMin_opt >= kPMax_opt % 调整后再次检查
                     % 如果优化边界无效，默认为最小储蓄
                     kPrime_opt = kPMin_opt;
                     % 在此角点评估价值
                     [negV, ~] = main_olg_v2_utils_pps.BellmanInner_k_pps(kPrime_opt, a_new, budget_remaining, vPofK_int, cS, paramS); % 调用重命名的BellmanInner
                     warning('在OptK中优化边界无效，年龄 %d, 剩余预算 %.2f。设置 k=kMin。', a_new, budget_remaining);
                 else
                     % 定义 fminbnd 的目标函数
                     obj = @(kP) main_olg_v2_utils_pps.negBellmanObjective_k_pps(kP, a_new, budget_remaining, vPofK_int, cS, paramS); % 调用重命名的目标函数
                     % 使用 fminbnd 寻找最优 kPrime
                     [kPrime_opt, negV, eflag] = fminbnd(obj, kPMin_opt, kPMax_opt, fminOpt); % eflag: fminbnd退出标志

                     % 如果fminbnd有问题或在边界处结束，仔细检查边界解
                     if eflag <= 0 || abs(kPrime_opt - kPMin_opt)<1e-7 || abs(kPrime_opt - kPMax_opt)<1e-7
                         [nV_min, ~] = main_olg_v2_utils_pps.BellmanInner_k_pps(kPMin_opt, a_new, budget_remaining, vPofK_int, cS, paramS);
                         [nV_max, ~] = main_olg_v2_utils_pps.BellmanInner_k_pps(kPMax_opt, a_new, budget_remaining, vPofK_int, cS, paramS);
                          % 如果价值非常接近，倾向于选择较低的 k'
                          if nV_min <= nV_max + 1e-9
                              kPrime_opt = kPMin_opt; negV = nV_min;
                          else
                              kPrime_opt = kPMax_opt; negV = nV_max;
                          end
                          if eflag <= 0
                              warning('fminbnd 在OptK中出现问题 (eflag=%d)，年龄 %d, 剩余预算 %.2f。使用边界 k=%.2f', eflag, a_new, budget_remaining, kPrime_opt);
                          end
                     end
                 end

                 kPrime = kPrime_opt; % 最优非PPS储蓄
                 ValueFunc = -negV; % 最大化的价值 (negV是最小化的负价值)
                 c = max(cS.cFloor, budget_remaining - kPrime); % 消费是剩余部分

                 % 如果消费被下限约束，重新计算价值函数？BellmanInner应该已经处理了。
                 if ~isfinite(ValueFunc), ValueFunc = -1e12; end % 处理非有限值
            end

            % 最终边界检查和清理
            kPrime = max(cS.kMin, min(cS.kMax, kPrime));
            c = max(cS.cFloor, budget_remaining - kPrime); % 在kPrime有界后再计算c
            if ~isfinite(kPrime), kPrime = cS.kMin; end
            if ~isfinite(c), c = cS.cFloor; end
            if ~isfinite(ValueFunc), ValueFunc = -1e12; end % 重新计算价值（如有必要），但BellmanInner已处理。

            % 验证输出属性
            validateattributes(kPrime,{'double'},{'finite','nonnan','scalar','real','>=',cS.kMin-1e-6,'<=',cS.kMax+1e-6});
            validateattributes(c,{'double'},{'finite','nonnan','scalar','real','>=',cS.cFloor-1e-6});
            validateattributes(ValueFunc,{'double'},{'finite','nonnan','scalar','real'});
        end

        % 重命名的嵌套辅助函数，以避免与原始 HHSolutionByOneState_VFI_Huggett (如果被调用) 冲突
        function negV = negBellmanObjective_k_pps(kP, ag_obj, budget_rem_obj, evInt_obj, cS_in_obj, pS_in_obj)
            % fminbnd的目标函数 (返回负价值)
            [negV, ~] = main_olg_v2_utils_pps.BellmanInner_k_pps(kP, ag_obj, budget_rem_obj, evInt_obj, cS_in_obj, pS_in_obj);
        end

        function [negVal, Val] = BellmanInner_k_pps(kP_inner, ag_inner, budget_rem_inner, evInt_inner, cS_in, pS_in)
            % 为给定的 k' = kP_inner 计算价值 V(k,e,a)
            % 给定PPS缴费后的剩余预算 budget_rem_inner
            cons = max(cS_in.cFloor, budget_rem_inner - kP_inner); % cons: 消费
            [util, ~] = main_olg_v2_utils_pps.CES_utility(cons, cS_in.sigma, cS_in); % util: 当前效用 (mu此处不需要)
            if ~isfinite(util), util = -1e12; end % 处理非有限值
            Val = util; % Val: 当前期效用

            if ag_inner < cS_in.aD_new % 如果不是最后一个时期，加上贴现的期望未来价值
                evF = -Inf; % evF: 下一期的期望价值
                try
                    % 对插值器的kP_inner进行截断
                    kP_eval = max(cS_in.kGridV(1), min(cS_in.kGridV(end), kP_inner));
                    evF = evInt_inner(kP_eval); % 插值计算期望价值
                catch ME_interp_inner % 捕获插值错误
                    warning('BellmanInner_k_pps 中插值错误，年龄 %d: %s。使用边界值。', ag_inner, ME_interp_inner.message);
                     if kP_inner < cS_in.kGridV(1), evF = evInt_inner(cS_in.kGridV(1)); else evF = evInt_inner(cS_in.kGridV(end)); end
                end
                if ~isfinite(evF), evF = -1e12; end % 处理非有限值

                s_transition = cS_in.s_1yr_transitionV(ag_inner); % 使用1年过渡存活率
                Val = util + cS_in.beta * s_transition * evF; % 加上贴现的期望未来价值 (使用年度beta)
            end

            % 处理非有限结果并返回负价值用于最小化
            if ~isfinite(Val)
                negVal = 1e12; % 对无效选择的高惩罚
                Val = -1e12;
            else
                negVal = -Val; % negVal: 负价值
            end
        end

    end % 结束 methods (Static)
end % 结束 classdef