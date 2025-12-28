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
            % 参数设定: 外生替代率, 内生PAYG税率, PPS
            % **核心参数调整为之前讨论的能产生较好行为的组合**

            % --- 人口统计与分组参数 (与你之前v2_pps版本一致) ---
            cS.age1_orig = 20; cS.ageLast_orig = 98; cS.ageRetire_orig = 65; cS.popGrowth_orig = 0.012;
            cS.aD_orig = cS.ageLast_orig - cS.age1_orig + 1; cS.aR_idx_orig = cS.ageRetire_orig - cS.age1_orig + 1;
            cS.aW_orig = cS.aR_idx_orig - 1; cS.physAgeV_orig = (cS.age1_orig : cS.ageLast_orig)';
            cS.d_orig = [0.00159,0.00169,0.00174,0.00172,0.00165,0.00156,0.00149,0.00145,0.00145,0.00149,0.00156,0.00163,0.00171,0.00181,0.00193,0.00207,0.00225,0.00246,0.00270,0.00299,0.00332,0.00368,0.00409,0.00454,0.00504,0.00558,0.00617,0.00686,0.00766,0.00865,0.00955,0.01058,0.01162,0.01264,0.01368,0.01475,0.01593,0.01730,0.01891,0.02074,0.02271,0.02476,0.02690,0.02912,0.03143,0.03389,0.03652,0.03930,0.04225,0.04538,0.04871,0.05230,0.05623,0.06060,0.06542,0.07066,0.07636,0.08271,0.08986,0.09788,0.10732,0.11799,0.12895,0.13920,0.14861,0.16039,0.17303,0.18665,0.20194,0.21877,0.23601,0.25289,0.26973,0.28612,0.30128,0.31416,0.32915,0.34450,0.36018];
            if length(cS.d_orig) ~= cS.aD_orig, error('d_orig 长度不匹配'); end
            cS.s_orig = 1 - cS.d_orig;
            % 静态人口质量计算 (如果HHPrices_Huggett备用逻辑需要)
            ageMassV_orig_static = ones(1, cS.aD_orig);
            for i = 2 : cS.aD_orig, ageMassV_orig_static(i) = cS.s_orig(i-1) * ageMassV_orig_static(i-1) / (1 + cS.popGrowth_orig); end
            if sum(ageMassV_orig_static)>0, cS.ageMassV_orig_static = ageMassV_orig_static ./ sum(ageMassV_orig_static); else, cS.ageMassV_orig_static = ones(1,cS.aD_orig)/cS.aD_orig; end;
            cS.yearStep = 5; cS.aD_new = ceil(cS.aD_orig / cS.yearStep); cS.aR_new = ceil(cS.aW_orig / cS.yearStep);
            cS.physAgeMap = cell(cS.aD_new, 1); for a = 1:cS.aD_new, startIdx = (a-1)*cS.yearStep + 1; endIdx = min(a*cS.yearStep, cS.aD_orig); cS.physAgeMap{a} = startIdx:endIdx; end
            cS.physAgeV_new = zeros(cS.aD_new, 1); for a = 1:cS.aD_new, cS.physAgeV_new(a) = cS.physAgeV_orig(cS.physAgeMap{a}(1)); end
            cS.s_1yr_transitionV = zeros(cS.aD_new, 1); for a = 1:(cS.aD_new - 1), lastYearIdx = cS.physAgeMap{a}(end); if lastYearIdx < cS.aD_orig, cS.s_1yr_transitionV(a) = cS.s_orig(lastYearIdx); else, cS.s_1yr_transitionV(a) = 0; end; end; cS.s_1yr_transitionV(cS.aD_new) = 0;
            % cS.retireMass_new (静态退休人口占比) 在这里不需要，因为b的计算会使用动态人口

            %% 人口动态参数
            cS.initial_pop = [76.2, 86.4, 113.8, 98.6, 86.6, 102.7, 112.0, 99.0, 64.0, 66.9, 44.1, 25.4, 14.9, 6.8, 1.7, 0.2];
            if length(cS.initial_pop) ~= cS.aD_new, cS.initial_pop = ones(cS.aD_new,1) * (100/cS.aD_new); warning('initial_pop mismatch, using uniform'); end
            beta_surv_pop = [0.998,0.996,0.994,0.992,0.988,0.984,0.980,0.976,0.970,0.960,0.945,0.920,0.880,0.800,0.680];
            if length(beta_surv_pop) ~= cS.aD_new - 1, error('beta_surv_pop length incorrect'); end
            cS.survivalProbV_popdyn = [beta_surv_pop(:)', 0];
            cS.bgp_tolerance = 0.001; cS.bgp_window = 5; cS.max_periods = 50;

            %% 家庭参数 (年度)
            cS.sigma      = 1.5;   % **使用之前能得到平滑结果的 sigma**
            cS.beta       = 1.011; % **使用之前能得到平滑结果的 beta**
            cS.cFloor     = 0.05;
            cS.nSim       = 5e4;

            %% 技术参数 (年度)
            cS.A          = 0.895944; cS.alpha = 0.36; cS.ddk = 0.06;

            %% 社会保障参数 (年度) - 外生替代率，内生PAYG税率
            cS.pension_replacement_rate = 0.30; % **例如，使用 0.30**
            cS.max_payg_payroll_tax_rate = 1.0; % PAYG税率上限

            %% 劳动禀赋参数
            cS.leSigma1 = 0.38^0.5; cS.leShockStd = 0.045^0.5; cS.lePersistence = 0.96; cS.leWidth = 4; cS.nw = 18;

            %% 网格参数 (kGridV代表非PPS资产的网格)
            cS.tgKY = 3; cS.tgWage = (1-cS.alpha)*cS.A*((cS.tgKY/cS.A)^(cS.alpha/(1-cS.alpha)));
            cS.nk = 50; % **如果之前nk=50时NoPPS能平滑，则保持；否则尝试增加到100或更高**
            cS.kMin = 0; cS.kMax = 100 * cS.tgWage;
            power = 1.5; kGridV = cS.kMin + (cS.kMax - cS.kMin) * (linspace(0, 1, cS.nk).^power);
            kGridV(1)=cS.kMin; cS.kGridV = kGridV(:);

            %% 年龄效率剖面
            ageEffV_orig_temp = zeros(100, 1);
            ageEffV_orig_temp(20:72) = [linspace(0.3,1.5,36-20+1), 1.5.*ones(1,47-37+1), linspace(1.5,0.2,65-48+1), linspace(0.18,0,72-66+1)];
            cS.ageEffV_orig = ageEffV_orig_temp(cS.age1_orig : cS.ageLast_orig);
            if length(cS.ageEffV_orig) ~= cS.aD_orig, error('ageEffV_orig length mismatch'); end
            cS.ageEffV_new = zeros(cS.aD_new, 1); for a = 1:cS.aD_new, cS.ageEffV_new(a) = mean(cS.ageEffV_orig(cS.physAgeMap{a})); end

            %% 个人养老金计划 (PPS) 参数
            cS.pps_active = true;
            cS.pps_max_contrib_frac = 0.1; % **设置为0，以便与NoPPS版本进行比较**
            cS.pps_tax_rate_withdrawal = 0.03;
            cS.pps_return_rate_premium = 0.0;
            cS.pps_withdrawal_rate = 0.15;
            cS.pps_contribution_age_max = cS.aR_idx_orig - 1;
            cS.pps_withdrawal_age_min = cS.aR_idx_orig;
            cS.pps_in_K = true;
            cS.pps_bequeathable = false;
        end

        function [R_market, MPL_gross] = HHPrices_Huggett(K_productive, L_total_eff, cS)
            % 计算市场回报率 R_market 和总MPL。
            % 这个版本不计算 PAYG 福利b 和净工资w，这些将在主循环中根据内生PAYG税率计算。
            % K_productive 是总生产性资本。

            if K_productive <= 0, K_productive=1e-6; warning('K_prod non-positive in HHPrices, set to small.'); end
            if L_total_eff <= 0, L_total_eff=1e-6; warning('L_eff non-positive in HHPrices, set to small.'); end

            Y_gross = cS.A * (K_productive^cS.alpha) * (L_total_eff^(1-cS.alpha));
            MPK_gross = cS.alpha * Y_gross / K_productive;
            MPL_gross = (1-cS.alpha) * Y_gross / L_total_eff;

            R_market   = 1 + MPK_gross - cS.ddk; % 假设无资本所得税影响市场总回报率的计算基础
            R_market   = max(1.0 + 1e-6, R_market);
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

         function incomeM_resources = HHIncome_Huggett(k_non_pps_grid, R_market_price, w_net_price, T_transfer_price, b_payg_price_group, a_new_idx_income, paramS_income, cS)
            % 计算VFI中可用的年度资源 (在PPS选择和非PPS储蓄k'选择之前)
            % k_non_pps_grid: 当前非PPS资产的网格点
            % R_market_price: 市场回报率因子 (通常是税前的，或一个基准回报率)
            % w_net_price: 净工资 (已扣除PAYG税和其他可能的劳动税)
            % T_transfer_price: 意外遗赠转移
            % b_payg_price_group: 该年龄组的PAYG养老金

            nk_income = length(k_non_pps_grid);
            leGridV_col_income = paramS_income.leGridV(:);
            nonCapIncomeNetM_income = zeros(1, cS.nw);

            if a_new_idx_income <= cS.aR_new % 工作年龄组
               net_wage_income_vec_income = w_net_price * cS.ageEffV_new(a_new_idx_income) .* leGridV_col_income;
               nonCapIncomeNetM_income(1,:) = net_wage_income_vec_income' + T_transfer_price + b_payg_price_group; % b对工作者通常为0
            else % 退休年龄组
               nonCapIncomeNetM_income(1,:) = b_payg_price_group + T_transfer_price;
            end

            % R_market_price 应该是家庭储蓄（非PPS）获得的实际回报率因子
            incomeM_resources = R_market_price * k_non_pps_grid(:) * ones(1, cS.nw) + ones(nk_income, 1) * nonCapIncomeNetM_income;

            if any(~isfinite(incomeM_resources(:))), warning('Non-finite income in HHIncome_Huggett group %d.',a_new_idx_income); incomeM_resources(~isfinite(incomeM_resources)) = 1e-6; end
         end

             function [cPolM, kPolM, cPpsPolM, valueM] = HHSolution_VFI_Huggett(R_vfi, w_vfi, T_vfi, bV_vfi, paramS_vfi, cS)
            % (与你上一轮提供的 main_olg_v2_utils_pps.m 中的此函数一致，它调用了下面的按年龄求解函数)
            cPolM  = zeros(cS.nk, cS.nw, cS.aD_new);
            kPolM  = zeros(cS.nk, cS.nw, cS.aD_new);
            cPpsPolM = zeros(cS.nk, cS.nw, cS.aD_new); % 初始化为0
            valueM = zeros(cS.nk, cS.nw, cS.aD_new);

            for a_vfi_loop_main = cS.aD_new : -1 : 1
               vPrimeM_vfi = [];
               if a_vfi_loop_main < cS.aD_new
                   vPrimeM_vfi = valueM(:,:,a_vfi_loop_main+1);
               end
               [cPolM(:,:,a_vfi_loop_main), kPolM(:,:,a_vfi_loop_main), cPpsPolM_temp, valueM(:,:,a_vfi_loop_main)] = ...
                  main_olg_v2_utils_pps.HHSolutionByAge_VFI_Huggett(a_vfi_loop_main, vPrimeM_vfi, ...
                  R_vfi, w_vfi, T_vfi, bV_vfi(a_vfi_loop_main), paramS_vfi, cS);
               
               % 只有当PPS实际可能贡献时才赋值 (以处理frac=0的情况)
               if cS.pps_active && cS.pps_max_contrib_frac > 1e-9
                   cPpsPolM(:,:,a_vfi_loop_main) = cPpsPolM_temp;
               end
            end
            validateattributes(cPolM, {'double'}, {'finite', 'nonnan', 'real','>=',cS.cFloor-1e-6,'size',[cS.nk,cS.nw,cS.aD_new]});
            validateattributes(kPolM, {'double'}, {'finite', 'nonnan', 'real','>=',cS.kMin-1e-6,'<=',cS.kMax+1e-6,'size',[cS.nk,cS.nw,cS.aD_new]});
            validateattributes(cPpsPolM, {'double'}, {'finite', 'nonnan', 'real','>=',0,'size',[cS.nk,cS.nw,cS.aD_new]});
            validateattributes(valueM, {'double'}, {'finite', 'nonnan', 'real','size',[cS.nk,cS.nw,cS.aD_new]});
        end

             function [cPolM_age, kPolM_age, cPpsPolM_age_actual, valueM_age_out] = HHSolutionByAge_VFI_Huggett(...
                a_idx, vPrime_keM_age, R_age, w_age, T_age, b_age_group, paramS_age, cS)
            %  并且当frac=0时cPpsPolM_age_actual会为0)
            if a_idx < cS.aD_new && isempty(vPrime_keM_age), error('VFI by Age: vPrime empty a<aD'); end
            if a_idx == cS.aD_new && ~isempty(vPrime_keM_age), error('VFI by Age: vPrime not empty a=aD'); end
            if a_idx < cS.aD_new, validateattributes(vPrime_keM_age,{'double'},{'finite','nonnan','real','size',[cS.nk,cS.nw]}); end

            budgetM_before_pps_age = main_olg_v2_utils_pps.HHIncome_Huggett(cS.kGridV, R_age, w_age, T_age, b_age_group, a_idx, paramS_age, cS);

            fminOpt_age=optimset('fminbnd'); fminOpt_age.TolX=1e-6; fminOpt_age.Display='none';

            cPolM_age = zeros(cS.nk, cS.nw);
            kPolM_age = zeros(cS.nk, cS.nw);
            cPpsPolM_age_actual = zeros(cS.nk, cS.nw); % 确保初始化为0
            valueM_age_out = -Inf(cS.nk, cS.nw);
            is_pps_contrib_age_val = (a_idx <= cS.aR_new);

            if a_idx == cS.aD_new
               cPolM_age = max(cS.cFloor, budgetM_before_pps_age);
               kPolM_age = cS.kMin*ones(cS.nk,cS.nw);
               [~, val_temp_age] = main_olg_v2_utils_pps.CES_utility(cPolM_age, cS.sigma, cS);
               valueM_age_out = val_temp_age;
               valueM_age_out(~isfinite(valueM_age_out)) = -1e12;
            else
               ExValuePrime_k_given_e_age = zeros(cS.nk, cS.nw);
                for ikPrime_age = 1:cS.nk, ExValuePrime_k_given_e_age(ikPrime_age, :) = vPrime_keM_age(ikPrime_age, :) * paramS_age.leTrProbM'; end
                ExValuePrime_k_given_e_age(~isfinite(ExValuePrime_k_given_e_age))=-1e10;
               vPInt_age=cell(1,cS.nw);
               for ie_age_interp=1:cS.nw % 使用与能收敛的简化模型相同的插值方法
                   vPInt_age{ie_age_interp} = griddedInterpolant(cS.kGridV, ExValuePrime_k_given_e_age(:, ie_age_interp), 'linear','linear');
               end
               % parfor ik_age_loop = 1:cS.nk % 可改为parfor
               parfor ik_age_loop = 1:cS.nk
                   c_row_age=zeros(1,cS.nw); k_row_age=zeros(1,cS.nw); cpps_row_for_this_k_state=zeros(1,cS.nw); v_row_age=zeros(1,cS.nw);
                   for ie_age_loop_inner = 1 : cS.nw
                        vPofK_interpolant_age = vPInt_age{ie_age_loop_inner}; budget_ike_before_pps_age = budgetM_before_pps_age(ik_age_loop,ie_age_loop_inner);
                        cPps_chosen_this_state_val = 0.0;
                        if is_pps_contrib_age_val && cS.pps_active && cS.pps_max_contrib_frac > 1e-9
                           wage_income_pps_limit_age = w_age * cS.ageEffV_new(a_idx) * paramS_age.leGridV(ie_age_loop_inner);
                           max_contrib_pps_age = wage_income_pps_limit_age * cS.pps_max_contrib_frac;
                           if budget_ike_before_pps_age - max_contrib_pps_age >= cS.cFloor + cS.kMin, cPps_chosen_this_state_val = max_contrib_pps_age;
                           else, cPps_chosen_this_state_val = min(max_contrib_pps_age, max(0, budget_ike_before_pps_age - cS.cFloor - cS.kMin )); end
                           cPps_chosen_this_state_val = max(0, cPps_chosen_this_state_val);
                        end
                        cpps_row_for_this_k_state(ie_age_loop_inner) = cPps_chosen_this_state_val;
                        budget_remaining_for_k_prime_and_c = budget_ike_before_pps_age - cPps_chosen_this_state_val;
                        [c_row_age(ie_age_loop_inner), k_row_age(ie_age_loop_inner), v_row_age(ie_age_loop_inner)] = ...
                              main_olg_v2_utils_pps.HHSolutionByOneState_OptK(a_idx, budget_remaining_for_k_prime_and_c, vPofK_interpolant_age, fminOpt_age, cS, paramS_age);
                   end
                   cPolM_age(ik_age_loop,:) = c_row_age; kPolM_age(ik_age_loop,:) = k_row_age;
                   cPpsPolM_age_actual(ik_age_loop,:) = cpps_row_for_this_k_state; valueM_age_out(ik_age_loop,:) = v_row_age;
               end
            end
        end

                function [c, kPrime, ValueFunc] = HHSolutionByOneState_OptK(a_new_osa_k, budget_remaining_osa_k, vPofK_int_osa_k, fminOpt_osa_k, cS, paramS_osa_k)
            % (与你之前v2_pps版本中的OptK逻辑相同)
            kPMin_osa_k = cS.kMin;
            kPMax_osa_k = budget_remaining_osa_k - cS.cFloor;
            ValueFunc = -Inf;

            % 将辅助函数定义为嵌套函数
            function negV_nested = negBellmanObjective_nested(kP_nested, ag_obj_nested, budget_obj_nested, evInt_obj_nested)
                [negV_nested, ~] = BellmanInner_nested(kP_nested, ag_obj_nested, budget_obj_nested, evInt_obj_nested);
            end

            function [negVal_nested, Val_nested] = BellmanInner_nested(kP_inner_nested, ag_inner_nested, budget_inner_nested, evInt_inner_nested)
                cons_nested = max(cS.cFloor, budget_inner_nested - kP_inner_nested); % 直接访问外部的cS
                [~, util_val_nested] = main_olg_v2_utils_pps.CES_utility(cons_nested, cS.sigma, cS); % 调用本类的CES_utility
                if ~isfinite(util_val_nested), util_val_nested = -1e12; end; Val_nested = util_val_nested;
                if ag_inner_nested < cS.aD_new
                    evF_val_nested = -Inf;
                    try 
                        kP_eval_nested = max(cS.kGridV(1),min(cS.kGridV(end),kP_inner_nested)); 
                        evF_val_nested=evInt_inner_nested(kP_eval_nested);
                    catch ME_bellman_nested
                        warning('BellmanInner_nested interp error age %d: %s', ag_inner_nested, ME_bellman_nested.message); 
                        if kP_inner_nested<cS.kGridV(1),evF_val_nested=evInt_inner_nested(cS.kGridV(1));else,evF_val_nested=evInt_inner_nested(cS.kGridV(end));end
                    end
                    if ~isfinite(evF_val_nested),evF_val_nested=-1e12;end
                    s_transition_nested = cS.s_1yr_transitionV(ag_inner_nested); 
                    Val_nested = util_val_nested + cS.beta*s_transition_nested*evF_val_nested;
                end
                if ~isfinite(Val_nested), negVal_nested=1e12; Val_nested=-1e12; else negVal_nested = -Val_nested; end
            end


            if kPMax_osa_k <= kPMin_osa_k
                kPrime = kPMin_osa_k; c = max(cS.cFloor, budget_remaining_osa_k - kPrime);
                 [~, u_c_osa_k] = main_olg_v2_utils_pps.CES_utility(c, cS.sigma, cS); 
                if a_new_osa_k < cS.aD_new
                    EV_c_osa_k = -Inf; try kP_eval_osa_k = max(cS.kGridV(1),min(cS.kGridV(end),kPrime)); EV_c_osa_k=vPofK_int_osa_k(kP_eval_osa_k); catch ME_osa_k; warning('OptK corner interp error age %d: %s',a_new_osa_k,ME_osa_k.message); if kPrime<cS.kGridV(1),EV_c_osa_k=vPofK_int_osa_k(cS.kGridV(1));else,EV_c_osa_k=vPofK_int_osa_k(cS.kGridV(end));end;end; if ~isfinite(EV_c_osa_k),EV_c_osa_k=-1e12;end
                    s_transition_osa_k = cS.s_1yr_transitionV(a_new_osa_k); ValueFunc = u_c_osa_k + cS.beta*s_transition_osa_k*EV_c_osa_k;
                else, ValueFunc = u_c_osa_k; end
                if ~isfinite(ValueFunc), ValueFunc = -1e12; end
            else
                kPMin_opt_osa_k = max(cS.kMin, kPMin_osa_k); kPMax_opt_osa_k = max(kPMin_opt_osa_k + 1e-9, min(cS.kMax, kPMax_osa_k));
                if kPMin_opt_osa_k >= kPMax_opt_osa_k
                    kPrime_opt_osa_k = kPMin_opt_osa_k; 
                    [negV_osa_k, ~] = BellmanInner_nested(kPrime_opt_osa_k, a_new_osa_k, budget_remaining_osa_k, vPofK_int_osa_k); % 调用嵌套函数
                    ValueFunc = -negV_osa_k; 
                    warning('OptK bounds invalid age %d',a_new_osa_k);
                else
                    obj_osa_k = @(kP_osa_k) negBellmanObjective_nested(kP_osa_k, a_new_osa_k, budget_remaining_osa_k, vPofK_int_osa_k); % 调用嵌套函数
                    [kPrime_opt_osa_k, negV_osa_k, eflag_osa_k] = fminbnd(obj_osa_k, kPMin_opt_osa_k, kPMax_opt_osa_k, fminOpt_osa_k);
                    if eflag_osa_k <= 0 || abs(kPrime_opt_osa_k-kPMin_opt_osa_k)<1e-7 || abs(kPrime_opt_osa_k-kPMax_opt_osa_k)<1e-7
                        [nV_min_osa_k,~]=BellmanInner_nested(kPMin_opt_osa_k,a_new_osa_k,budget_remaining_osa_k,vPofK_int_osa_k); 
                        [nV_max_osa_k,~]=BellmanInner_nested(kPMax_opt_osa_k,a_new_osa_k,budget_remaining_osa_k,vPofK_int_osa_k);
                        if nV_min_osa_k <= nV_max_osa_k + 1e-9, kPrime_opt_osa_k=kPMin_opt_osa_k; negV_osa_k=nV_min_osa_k; else, kPrime_opt_osa_k=kPMax_opt_osa_k; negV_osa_k=nV_max_osa_k; end
                        if eflag_osa_k <= 0, warning('fminbnd issue OptK age %d, eflag %d',a_new_osa_k,eflag_osa_k); end
                    end
                    ValueFunc = -negV_osa_k;
                end
                kPrime = kPrime_opt_osa_k; c = max(cS.cFloor, budget_remaining_osa_k - kPrime);
                if ~isfinite(ValueFunc), ValueFunc = -1e12; end
            end
            kPrime = max(cS.kMin, min(cS.kMax, kPrime)); c = max(cS.cFloor, budget_remaining_osa_k - kPrime);
            [~, u_final_c_osa_k] = main_olg_v2_utils_pps.CES_utility(c, cS.sigma, cS);
            if a_new_osa_k < cS.aD_new
                EV_A_final_osa_k = -Inf; try kP_eval_final_osa_k = max(cS.kGridV(1),min(cS.kGridV(end),kPrime)); EV_A_final_osa_k=vPofK_int_osa_k(kP_eval_final_osa_k); catch; if kPrime<cS.kGridV(1),EV_A_final_osa_k=vPofK_int_osa_k(cS.kGridV(1));else,EV_A_final_osa_k=vPofK_int_osa_k(cS.kGridV(end));end;end; if ~isfinite(EV_A_final_osa_k),EV_A_final_osa_k=-1e12;end
                s_transition_final_osa_k = cS.s_1yr_transitionV(a_new_osa_k); ValueFunc = u_final_c_osa_k + cS.beta*s_transition_final_osa_k*EV_A_final_osa_k;
            else, ValueFunc = u_final_c_osa_k; end
            if ~isfinite(kPrime), kPrime = cS.kMin; end; if ~isfinite(c), c = cS.cFloor; end; if ~isfinite(ValueFunc), ValueFunc = -1e12; end
            validateattributes(kPrime,{'double'},{'finite','nonnan','scalar','real','>=',cS.kMin-1e-6,'<=',cS.kMax+1e-6});
            validateattributes(c,{'double'},{'finite','nonnan','scalar','real','>=',cS.cFloor-1e-6});
            validateattributes(ValueFunc,{'double'},{'finite','nonnan','scalar','real'});
        end

        % 删除之前的静态方法 BellmanInner_k_pps 和 negBellmanObjective_k_pps
        % 因为它们现在是 HHSolutionByOneState_OptK 的嵌套函数了。
        % function negV = negBellmanObjective_k_pps(...) end
        % function [negVal, Val] = BellmanInner_k_pps(...) end

        function negV = negBellmanObjective_k_pps(kP, ag_obj, budget_rem_obj, evInt_obj, cS_in_obj, pS_in_obj)
            % (与之前v2_pps版本相同)
            [negV, ~] = main_olg_v2_utils_pps.BellmanInner_k_pps(kP, ag_obj, budget_rem_obj, evInt_obj, cS_in_obj, pS_in_obj);
        end

        function [negVal, Val] = BellmanInner_k_pps(kP_inner, ag_inner, budget_rem_inner, evInt_inner, cS_in, pS_in)
            % (与之前v2_pps版本相同)
            cons = max(cS_in.cFloor, budget_rem_inner - kP_inner);
            [util_val, ~] = main_olg_v2_utils_pps.CES_utility(cons, cS_in.sigma, cS_in); % 使用本类的方法
            if ~isfinite(util_val), util_val = -1e12; end; Val = util_val;
            if ag_inner < cS_in.aD_new
                evF_val = -Inf;
                try kP_eval_val = max(cS_in.kGridV(1),min(cS_in.kGridV(end),kP_inner)); evF_val=evInt_inner(kP_eval_val);
                catch ME_bellman_pps; warning('BellmanInner_k_pps interp error age %d: %s', ag_inner, ME_bellman_pps.message); if kP_inner<cS_in.kGridV(1),evF_val=evInt_inner(cS_in.kGridV(1));else,evF_val=evInt_inner(cS_in.kGridV(end));end; end
                if ~isfinite(evF_val),evF_val=-1e12;end
                s_transition_val = cS_in.s_1yr_transitionV(ag_inner); Val = util_val + cS_in.beta*s_transition_val*evF_val;
            end
            if ~isfinite(Val), negVal=1e12; Val=-1e12; else negVal = -Val; end
        end

    end % End methods (Static)
end % End classdef