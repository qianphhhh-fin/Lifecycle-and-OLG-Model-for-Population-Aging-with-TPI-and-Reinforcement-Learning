% --- 开始文件：main_olg_v8_utils.m ---

% =====================================================================
% === OLG 模型 V8 工具函数库: 内生PPS缴费决策的数值实现 ===
% =====================================================================
% 
% 理论基础：
% 本工具函数库实现了v8.tex理论模型中描述的所有核心算法和数值方法
% 包括值函数迭代(VFI)、家庭问题求解、一般均衡算法等关键组件
% 
% 主要功能模块：
% 1. 参数设定函数 - 校准模型所有参数
% 2. 人口动态函数 - 模拟人口老龄化过程  
% 3. 劳动禀赋过程 - Tauchen方法离散化AR(1)过程
% 4. 家庭效用函数 - CRRA效用和边际效用计算
% 5. VFI核心算法 - 内生PPS缴费的值函数迭代求解
% 6. 宏观经济函数 - 要素价格和市场出清条件
% 7. 一般均衡求解器 - 迭代求解K和tau_l的均衡值
% 
% V8模型核心创新：
% - HHSolutionByAge_VFI_Huggett_v8: 实现内生PPS缴费选择的VFI
% - 对每个状态(k,k_pps,ε,a)，个体选择最优的(c,k',c_pps)组合
% - PPS缴费受双重约束：收入比例上限和年度绝对上限
% - 通过离散网格搜索实现PPS缴费的优化选择
% =====================================================================

classdef main_olg_v8_utils % 定义 main_olg_v8_utils 类
    % OLG 模型 V8 的工具函数 (内生PPS缴费决策)
    % - 基于 V7/Baseline 版本修改，增加内生PPS缴费优化
    % - HHSolutionByAge_VFI_Huggett_v8 实现对 c_pps 和 k' 的联合优化

    methods (Static) % 定义静态方法

        % =====================================================================
        % == 模型参数设定函数 ==
        % =====================================================================
        % 对应v8.tex第4节"参数设定和参数校准"
        function cS = ParameterValues_HuggettStyle()
            % ParameterValues_HuggettStyle - 设置OLG模型V8的所有参数
            % 
            % 输出：
            %   cS - 包含所有模型参数的结构体
            % 
            % 参数分类：
            % 1. 人口动态参数：年龄结构、存活率、生育率
            % 2. 经济参数：偏好、生产技术、政府政策
            % 3. 养老金参数：PAYG和PPS制度设计
            % 4. PPS约束参数：缴费上限和提取规则（V8新增）

            % --- 人口结构基础参数（对应v8.tex人口动态设定） ---
            cS.age1_orig = 20;              % 模型起始年龄（岁）
            cS.ageLast_orig = 98;           % 模型终止年龄（岁）
            cS.ageRetire_orig = 65;         % 退休年龄（岁）
            cS.popGrowth_orig = 0.012;      % 原始人口增长率
            cS.aD_orig = cS.ageLast_orig - cS.age1_orig + 1;        % 年度年龄组数
            cS.aR_idx_orig = cS.ageRetire_orig - cS.age1_orig + 1;  % 退休年度年龄索引
            cS.aW_orig = cS.aR_idx_orig - 1;                        % 工作年数
            cS.physAgeV_orig = (cS.age1_orig : cS.ageLast_orig)';   % 年度年龄向量

            % --- 年度死亡率数据（基于中国生命表） ---
            % 对应v8.tex中的β_{surv,a-1,t-1}存活率参数
            cS.d_orig = [0.00159,0.00169,0.00174,0.00172,0.00165,0.00156,0.00149,0.00145,0.00145,0.00149,0.00156,0.00163,0.00171,0.00181,0.00193,0.00207,0.00225,0.00246,0.00270,0.00299,0.00332,0.00368,0.00409,0.00454,0.00504,0.00558,0.00617,0.00686,0.00766,0.00865,0.00955,0.01058,0.01162,0.01264,0.01368,0.01475,0.01593,0.01730,0.01891,0.02074,0.02271,0.02476,0.02690,0.02912,0.03143,0.03389,0.03652,0.03930,0.04225,0.04538,0.04871,0.05230,0.05623,0.06060,0.06542,0.07066,0.07636,0.08271,0.08986,0.09788,0.10732,0.11799,0.12895,0.13920,0.14861,0.16039,0.17303,0.18665,0.20194,0.21877,0.23601,0.25289,0.26973,0.28612,0.30128,0.31416,0.32915,0.34450,0.36018];
            if length(cS.d_orig) ~= cS.aD_orig
                error('年度死亡率数据 d_orig 长度与年龄跨度不匹配');
            end
            cS.s_orig = 1 - cS.d_orig;      % 年度存活率

            % --- 年龄组聚合参数（将年度年龄聚合为5年期年龄组） ---
            % 对应v8.tex模型中的16个5年期年龄组设定
            cS.yearStep = 5;                % 每个年龄组跨度（年）
            cS.aD_new = ceil(cS.aD_orig / cS.yearStep);    % 年龄组数量
            cS.aR_new = ceil(cS.aW_orig / cS.yearStep);    % 工作年龄组数量

            % 建立年度年龄到年龄组的映射关系
            cS.physAgeMap = cell(cS.aD_new, 1);
            for a = 1:cS.aD_new
                startIdx = (a-1)*cS.yearStep + 1;
                endIdx = min(a*cS.yearStep, cS.aD_orig);
                cS.physAgeMap{a} = startIdx:endIdx;
            end

            % 计算各年龄组代表性年龄
            cS.physAgeV_new = zeros(cS.aD_new, 1);
            for a = 1:cS.aD_new
                cS.physAgeV_new(a) = cS.physAgeV_orig(cS.physAgeMap{a}(1));
            end

            % 计算年龄组间转移存活率（用于VFI中的期望效用计算）
            cS.s_1yr_transitionV = zeros(cS.aD_new, 1);
            for a = 1:(cS.aD_new - 1)
                lastYearIdxInGroup = cS.physAgeMap{a}(end);
                if lastYearIdxInGroup < cS.aD_orig
                    cS.s_1yr_transitionV(a) = cS.s_orig(lastYearIdxInGroup);
                else
                    cS.s_1yr_transitionV(a) = 0;
                end
            end
            cS.s_1yr_transitionV(cS.aD_new) = 0;

            % --- 初始人口分布（基于2023年中国人口结构） ---
            cS.initial_pop = [76.2, 86.4, 113.8, 98.6, 86.6, 102.7, 112.0, 99.0, 64.0, 66.9, 44.1, 25.4, 14.9, 6.8, 1.7, 0.2];
            if length(cS.initial_pop) ~= cS.aD_new
                cS.initial_pop = ones(cS.aD_new,1) * (100/cS.aD_new);
                warning('initial_pop长度与年龄组数不匹配，已重设为均匀分布。');
            end

            % --- 年龄组间存活率（用于人口动态模拟） ---
            % 对应v8.tex公式中的β_{surv,a-1,t-1}参数
            beta_surv_pop = [0.998,0.996,0.994,0.992,0.988,0.984,0.980,0.976,0.970,0.960,0.945,0.920,0.880,0.800,0.680];
            if length(beta_surv_pop) ~= cS.aD_new - 1
                error('年龄组间存活率 beta_surv_pop 的长度对于 %d 个年龄组不正确。应为 %d。', cS.aD_new, cS.aD_new -1);
            end
            cS.survivalProbV_popdyn = [beta_surv_pop(:)', 0];

            % --- 人口动态收敛参数 ---
            cS.bgp_tolerance = 0.001;       % 稳态收敛容忍度
            cS.bgp_window = 5;              % 稳态检测窗口期
            cS.max_periods = 50;            % 最大模拟期数

            % --- 家庭偏好参数（对应v8.tex第2.2.1节） ---
            cS.sigma      = 1.5;            % 相对风险厌恶系数γ
            cS.beta       = 1.011;          % 主观贴现因子β_disc  
            cS.cFloor     = 0.05;           % 最低消费约束
            cS.nSim       = 1000;           % 蒙特卡洛模拟个体数

            % --- 生产技术参数（对应v8.tex第2.3节） ---
            cS.A          = 0.895944;       % 全要素生产率
            cS.alpha      = 0.36;           % 资本产出弹性
            cS.ddk        = 0.06;           % 资本折旧率δ

            % --- 政府财政参数（对应v8.tex第2.4节） ---
            cS.tau_k = 0.20;                % 资本所得税率
            cS.tau_c = 0.10;                % 消费税率
            cS.gov_exp_frac_Y = 0.15;       % 政府支出占GDP比例
            cS.gov_debt_frac_Y = 0.60;      % 政府债务占GDP比例

            % --- 劳动效率冲击过程参数 ---
            % 对应v8.tex中的ε_{a,t}随机过程
            cS.leSigma1 = 0.38^0.5;         % 初期效率分布标准差
            cS.leShockStd = 0.045^0.5;      % 效率冲击标准差
            cS.lePersistence = 0.96;        % AR(1)持续性参数
            cS.leWidth = 4;                 % Tauchen方法的标准差倍数
            cS.nw = 5;                      % 效率状态网格点数

            % --- 资产网格参数 ---
            % 非PPS资产网格（对应v8.tex中的k状态空间）
            cS.tgKY = 3;                    % 目标资本产出比
            cS.tgWage = (1-cS.alpha)*cS.A*((cS.tgKY/cS.A)^(cS.alpha/(1-cS.alpha)));
            cS.nk = 10;                     % 非PPS资产网格点数
            cS.kMin = 0;                    % 非PPS资产下界
            cS.kMax = 100 * cS.tgWage;      % 非PPS资产上界
            power = 1.5;                    % 网格密度参数
            kGridV = cS.kMin + (cS.kMax - cS.kMin) * (linspace(0, 1, cS.nk).^power);
            if cS.nk > 0
                kGridV(1)=cS.kMin;
            end
            cS.kGridV = kGridV(:);

            % --- 年龄效率剖面（对应v8.tex中的h_{a,t}） ---
            % 基于中国劳动者生产率-年龄关系校准
            ageEffV_orig_temp = zeros(100, 1);
            ageEffV_orig_temp(20:72) = [linspace(0.3,1.5,36-20+1), 1.5.*ones(1,47-37+1), linspace(1.5,0.2,65-48+1), linspace(0.18,0,72-66+1)];
            cS.ageEffV_orig = ageEffV_orig_temp(cS.age1_orig : cS.ageLast_orig);
            if length(cS.ageEffV_orig) ~= cS.aD_orig
                error('ageEffV_orig 年度年龄效率剖面长度不匹配');
            end
            % 计算年龄组平均效率
            cS.ageEffV_new = zeros(cS.aD_new, 1);
            for a = 1:cS.aD_new
                cS.ageEffV_new(a) = mean(cS.ageEffV_orig(cS.physAgeMap{a}));
            end

            % === V8模型核心：PPS制度参数设计 ===
            % 对应v8.tex第2.2.2节"PPS缴费约束"部分

            % --- PPS制度基础参数 ---
            cS.pps_active = true;                        % PPS制度激活标志
            cS.pps_tax_rate_withdrawal = 0.03;          % PPS提取阶段税率
            cS.pps_return_rate_premium = 0.02;          % PPS超额收益率
            cS.pps_withdrawal_rate = 0.15;              % 退休后年度提取比例
            cS.pps_in_K = true;                         % PPS资产是否计入生产性资本
            cS.pps_bequeathable = false;                % PPS资产是否可遗赠

            % --- V8模型关键创新：PPS缴费约束参数 ---
            % 实现v8.tex中描述的双重约束机制
            cS.pps_annual_contrib_limit = 1.5;          % PPS年度绝对缴费上限
            cS.pps_max_contrib_frac = 0.15;             % PPS缴费占劳动收入比例上限
            cS.pps_contribution_age_max_idx = cS.aR_idx_orig - 1;  % 最大缴费年度年龄
            cS.pps_withdrawal_age_min_idx = cS.aR_idx_orig;        % 最低提取年度年龄
            
            % V8核心：PPS缴费选择网格点数
            % 用于在VFI中离散化PPS缴费选择空间
            cS.n_pps_choice_grid_points = 10;           % PPS缴费选择网格点数

            % --- PPS资产网格（对应v8.tex中的k_pps状态空间） ---
            cS.nkpps = 10;                               % PPS资产网格点数
            cS.kppsMin = 0;                             % PPS资产下界
            cS.kppsMax = cS.kMax / 2;                   % PPS资产上界
            if cS.nkpps > 0
                cS.kppsMax = max(cS.kppsMax, 1e-3);
            end
            power_kpps = 1.5;                           % PPS资产网格密度参数
            if cS.nkpps > 1
                kppsGridV_temp = cS.kppsMin + (cS.kppsMax - cS.kppsMin) * (linspace(0, 1, cS.nkpps).^power_kpps);
                kppsGridV_temp(1) = cS.kppsMin;
            elseif cS.nkpps == 1
                kppsGridV_temp = cS.kppsMin;
            else
                kppsGridV_temp = [];
            end
            cS.kppsGridV = kppsGridV_temp(:);

            % --- 一般均衡求解参数 ---
            cS.max_iter_K_tau_l = 100;                 % K和tau_l迭代最大次数
            cS.tol_K_tau_l = 1e-4;                     % K和tau_l收敛容忍度
            cS.damp_K_v5 = 0.3;                        % K更新阻尼系数
            cS.damp_tau_l_v5 = 0.3;                    % tau_l更新阻尼系数
            cS.gbc_tol_for_internal_loop = 1e-3;       % 政府预算平衡容忍度

            % --- 收敛检测参数 ---
            cS.max_stagnation_iters = 10;              % 最大停滞迭代次数
            cS.min_norm_improvement_frac = 1e-3;       % 最小改进比例
            cS.max_tau_l_boundary_strikes = 5;         % tau_l边界冲击最大次数
            
            % --- PAYG税率约束参数 ---
            cS.tau_l_init_guess = 0.05;                % 所得税率初始猜测
            cS.tau_l_min = 0.00;                       % 所得税率下界
            cS.tau_l_max = 0.3;                        % 所得税率上界
            cS.max_total_labor_tax = 0.6;              % 总劳动税负上限
            cS.theta_payg_max = 0.35;                  % PAYG税率上限

            % --- 参考性PPS缴费时间表（V8中不再直接使用） ---
            % 保留用于与V7模型比较和参数参考
            cS.pps_fixed_contrib_schedule_frac = zeros(cS.aD_new, 1);
            num_working_age_groups = cS.aR_new;
            if num_working_age_groups > 0
                if num_working_age_groups == 1
                     cS.pps_fixed_contrib_schedule_frac(1:num_working_age_groups) = 0.05;
                elseif num_working_age_groups > 1
                    mid_point1 = ceil(num_working_age_groups / 3);
                    mid_point2 = ceil(2 * num_working_age_groups / 3);
                    if mid_point1 > 0
                        cS.pps_fixed_contrib_schedule_frac(1:mid_point1) = linspace(0.02, 0.06, mid_point1);
                    end
                    if mid_point2 > mid_point1
                        cS.pps_fixed_contrib_schedule_frac(mid_point1+1:mid_point2) = linspace(0.06, 0.10, mid_point2 - mid_point1);
                    end
                    if num_working_age_groups > mid_point2
                        cS.pps_fixed_contrib_schedule_frac(mid_point2+1:num_working_age_groups) = linspace(0.10, 0.04, num_working_age_groups - mid_point2);
                    end
                end
            end
            if cS.aR_new < cS.aD_new
                cS.pps_fixed_contrib_schedule_frac(cS.aR_new + 1 : cS.aD_new) = 0;
            end
            
            % --- 数值优化参数 ---
            cS.fminbnd_TolX = 1e-6;                    % fminbnd容忍度
            cS.fminbnd_Display = 'none';               % fminbnd显示设置

            fprintf('V8: 完整参数已设置完毕。\n');
            fprintf('    主要创新：内生PPS缴费选择，缴费网格点数=%d\n', cS.n_pps_choice_grid_points);
            fprintf('    PPS约束：收入比例上限=%.1f%%, 年度绝对上限=%.2f\n', ...
                cS.pps_max_contrib_frac*100, cS.pps_annual_contrib_limit);
        end

        % =====================================================================
        % == 人口动态模拟函数 ==
        % =====================================================================
        % 对应v8.tex第2.1节"人口动态设定"，实现人口老龄化过程模拟
        
        function popS = initPopulation(cS)
            % initPopulation - 初始化人口结构
            % 
            % 输入：
            %   cS - 参数结构体
            % 输出：
            %   popS - 包含初始人口分布的结构体
            % 
            % 功能：基于2023年中国人口结构设置16个年龄组的初始分布
            popS.Z = zeros(cS.aD_new, 1);
            initial_total = sum(cS.initial_pop);

            if initial_total > 0 && length(cS.initial_pop) == cS.aD_new
                popS.Z(:, 1) = cS.initial_pop(:) / initial_total * 100;
            else
                warning('初始人口数据不匹配或总和为零。将设置为均匀的初始年龄组人口分布。');
                popS.Z(:, 1) = 100 / cS.aD_new;
            end

            popS.totalPop = sum(popS.Z(:, 1));

            if popS.totalPop > 1e-9
                popS.ageDist = popS.Z(:, 1) / popS.totalPop;
            else
                popS.ageDist = zeros(cS.aD_new, 1);
            end
            popS.initialAgeDist = popS.ageDist;
            fprintf('初始年龄组人口已设置。总人口=%.2f (代表百分比基数)。\n', popS.totalPop);
        end

        function popS = populationDynamics(popS, cS)
            % populationDynamics - 模拟人口动态演进
            % 
            % 输入：
            %   popS - 初始人口结构
            %   cS - 参数结构体
            % 输出：
            %   popS - 包含人口演进历史的结构体
            % 
            % 功能：实现v8.tex公式Z_{a,t} = β_{surv,a-1,t-1}Z_{a-1,t-1}的人口转移
            max_periods_sim = cS.max_periods;
            Z_history = zeros(cS.aD_new, max_periods_sim + 1);
            totalPop_history = zeros(1, max_periods_sim + 1);
            ageDist_history = zeros(cS.aD_new, max_periods_sim + 1);

            Z_history(:, 1) = popS.Z(:, 1);
            totalPop_history(1) = popS.totalPop(1);
            ageDist_history(:, 1) = popS.ageDist(:, 1);

            fprintf('人口动态模拟开始 (年龄组, 最大期数 = %d)...\n', max_periods_sim);
            bgp_reached_flag = false;
            actual_periods_run = max_periods_sim;

            for t = 1:max_periods_sim
                if mod(t, 10) == 0 || t == 1
                    fprintf('  模拟人口期数 %d (年龄组)\n', t);
                end

                Z_current_period = Z_history(:, t);
                Z_next_period = zeros(cS.aD_new, 1);

                time_varying_growth_rate = 0;
                if t < 6
                    time_varying_growth_rate = -0.01 - 0.003 * t;
                else
                    time_varying_growth_rate = -0.03 - 0.004 * min(t-6, 10);
                end
                Z_next_period(1) = Z_current_period(1) * (1 + time_varying_growth_rate);
                Z_next_period(1) = max(0, Z_next_period(1));

                for a = 2:cS.aD_new
                    survival_prob_group = 0;
                    if (a-1) >= 1 && (a-1) <= length(cS.survivalProbV_popdyn)
                        survival_prob_group = cS.survivalProbV_popdyn(a-1);
                    end
                    Z_next_period(a) = Z_current_period(a-1) * survival_prob_group;
                    Z_next_period(a) = max(0, Z_next_period(a));
                end

                Z_history(:, t+1) = Z_next_period;
                totalPop_history(t+1) = sum(Z_next_period);
                if totalPop_history(t+1) > 1e-9
                    ageDist_history(:, t+1) = Z_next_period / totalPop_history(t+1);
                else
                    ageDist_history(:, t+1) = 0;
                    totalPop_history(t+1) = 0;
                end

                current_check_period_idx = t + 1;
                if current_check_period_idx >= cS.bgp_window + 1
                    stable = true;
                    for w_idx = 1:cS.bgp_window
                       hist_idx1 = current_check_period_idx - w_idx + 1;
                       hist_idx2 = current_check_period_idx - w_idx;
                       if hist_idx1 > 0 && hist_idx2 > 0 && hist_idx1 <= size(ageDist_history,2) && hist_idx2 <= size(ageDist_history,2)
                           change = norm(ageDist_history(:, hist_idx1) - ageDist_history(:, hist_idx2));
                           if change >= cS.bgp_tolerance
                               stable = false;
                               break;
                           end
                       else
                           stable = false; break;
                       end
                    end
                    if stable
                        fprintf('\n人口稳态 (年龄组) 在模拟期数 %d (对应历史数据索引 %d) 达到。\n', t, current_check_period_idx);
                        bgp_reached_flag = true;
                        actual_periods_run = t;
                        break;
                    end
                end
            end

            final_period_idx_to_store = min(actual_periods_run + 1, size(Z_history,2));
            popS.Z = Z_history(:, 1:final_period_idx_to_store);
            popS.totalPop = totalPop_history(1:final_period_idx_to_store);
            popS.ageDist = ageDist_history(:, 1:final_period_idx_to_store);

            depRatio_history = zeros(1, actual_periods_run);
            for th_loop = 1:actual_periods_run
                 Z_t_for_depratio = Z_history(:, th_loop + 1);
                 working_pop = sum(Z_t_for_depratio(1:cS.aR_new));
                 retired_pop = sum(Z_t_for_depratio(cS.aR_new+1:end));
                 if working_pop > 1e-9
                     depRatio_history(th_loop) = retired_pop / working_pop;
                 else
                     depRatio_history(th_loop) = inf;
                 end
            end
            popS.dependencyRatio = depRatio_history;

            fprintf('人口动态模拟完成。运行期数: %d。达到BGP: %d\n', actual_periods_run, bgp_reached_flag);
            if ~bgp_reached_flag
                fprintf('警告: 人口稳态 (年龄组) 未在 %d 期内达到。\n', max_periods_sim);
            end
        end

        function [Z_ss, dependency_ratio_ss, bgp_reached, bgp_period] = detectSteadyStatePopulation(popS, cS)
            % (Identical to main_olg_baseline_utils.m)
            actual_periods_in_data = size(popS.Z, 2);
            bgp_reached = false;
            bgp_period = actual_periods_in_data - 1;

            if actual_periods_in_data < cS.bgp_window + 1
                fprintf('人口模拟期数过短 (%d 数据点)，无法进行稳态检查 (窗口期 = %d)。\n', actual_periods_in_data, cS.bgp_window);
            else
                fprintf('检查人口稳态 (年龄组, 最近 %d 期)...\n', cS.bgp_window);
                for t_check_end_idx = actual_periods_in_data : -1 : cS.bgp_window + 1
                    stable = true;
                    for w_idx = 0 : (cS.bgp_window - 1)
                        idx1 = t_check_end_idx - w_idx;
                        idx2 = t_check_end_idx - w_idx - 1;
                        if idx1 > 0 && idx2 > 0 && idx1 <= size(popS.ageDist,2) && idx2 <= size(popS.ageDist,2)
                            change = norm(popS.ageDist(:, idx1) - popS.ageDist(:, idx2));
                            if change >= cS.bgp_tolerance
                                stable = false;
                                break;
                            end
                        else
                            stable = false; break;
                        end
                    end
                    if stable
                        bgp_reached = true;
                        bgp_period = t_check_end_idx - 1;
                        fprintf('人口稳态 (年龄组) 从模拟期数 %d (数据索引 %d) 开始检测到 (稳定窗口结束于此)。\n', bgp_period, t_check_end_idx);
                        break;
                    end
                end
                if ~bgp_reached
                     fprintf('未检测到人口稳态 (年龄组)。将使用最终期数据。\n');
                     bgp_period = actual_periods_in_data - 1;
                end
            end

            ss_data_index = min(bgp_period + 1, size(popS.Z, 2));
            Z_ss = popS.Z(:, ss_data_index);

            Z_ss_norm = zeros(cS.aD_new, 1);
            if sum(Z_ss) > 1e-9
                Z_ss_norm = Z_ss / sum(Z_ss);
            end
            
            dependency_ratio_ss = NaN; 
            if isfield(popS, 'dependencyRatio') && ~isempty(popS.dependencyRatio)
                valid_dep_ratio_index = min(max(1, bgp_period), length(popS.dependencyRatio)); % Ensure index is valid
                 if bgp_period == 0 && length(popS.dependencyRatio) >=1 % if initial period is SS, use first dep ratio
                    valid_dep_ratio_index = 1;
                 end
                 if valid_dep_ratio_index > 0 && valid_dep_ratio_index <= length(popS.dependencyRatio)
                    dependency_ratio_ss = popS.dependencyRatio(valid_dep_ratio_index);
                 end
            end
            if isnan(dependency_ratio_ss) % Fallback if not found or invalid
                 working_pop_ss = sum(Z_ss(1:cS.aR_new));
                 retired_pop_ss = sum(Z_ss(cS.aR_new+1:end));
                 if working_pop_ss > 1e-9
                     dependency_ratio_ss = retired_pop_ss / working_pop_ss;
                 else
                     dependency_ratio_ss = inf;
                 end
                 if (~isfield(popS, 'dependencyRatio') || isempty(popS.dependencyRatio)) && bgp_period > 0
                     warning('抚养比历史未找到或过短，已基于Z_ss重新计算。');
                 end
            end


            figure('Name', 'V8: 初始 vs 稳态/最终 年龄组人口分布');
            hold on;
            group_indices = 1:cS.aD_new;
            if isfield(popS, 'initialAgeDist') && ~isempty(popS.initialAgeDist)
                bar(group_indices - 0.2, popS.initialAgeDist * 100, 0.4, 'b', 'DisplayName', '初始年龄组分布');
            else
                warning('未找到或空的初始年龄分布用于绘图。');
            end
            bar(group_indices + 0.2, Z_ss_norm * 100, 0.4, 'r', 'DisplayName', sprintf('稳态年龄组分布 (模拟期 %d)', bgp_period));
            hold off;
            xlabel(sprintf('年龄组索引 (1 至 %d)', cS.aD_new));
            ylabel('占总人口百分比 (%)');
            title(sprintf('V8: 初始 vs 稳态/最终 年龄组人口分布 (稳态代表模拟期 t=%d)', bgp_period));
            legend('Location', 'best');
            xticks(group_indices);
            grid on;
            drawnow;
            fprintf('已绘制初始与稳态/最终年龄组人口分布图。\n');
        end

        % =====================================================================
        % == 劳动禀赋过程与收入分布函数 ==
        % =====================================================================
        % 对应v8.tex第2.2节中的劳动效率冲击ε_{a,t}的随机过程建模
        
        function [logGridV, trProbM, prob1V] = EarningProcess_olgm(cS)
            % EarningProcess_olgm - 生成劳动效率冲击的离散状态空间
            % 
            % 输入：
            %   cS - 参数结构体
            % 输出：
            %   logGridV - 劳动效率对数值网格
            %   trProbM - 状态转移概率矩阵
            %   prob1V - 初始分布概率向量
            % 
            % 功能：使用Tauchen方法将AR(1)劳动效率过程离散化为有限状态Markov链
            [logGridV_raw, trProbM_calc] = main_olg_v8_utils.tauchen(cS.nw, cS.lePersistence, cS.leShockStd, 0, cS.leWidth);
            gridMin_raw = logGridV_raw(1);
            gridMax_raw = logGridV_raw(end);
            normGridMin = gridMin_raw - 2*cS.leShockStd;
            normGridMax = gridMax_raw + 2*cS.leShockStd;
            prob1V_calc = [];
            try
                [prob1V_calc, ~, ~] = main_olg_v8_utils.norm_grid(logGridV_raw, normGridMin, normGridMax, 0, cS.leSigma1);
                prob1V_calc = prob1V_calc(:);
            catch ME
                warning('norm_grid 计算初始分布失败: %s。将使用均匀分布作为备用。', ME.message);
                prob1V_calc = ones(cS.nw, 1)/cS.nw;
            end

            logGridV_calc = logGridV_raw(:);
            logGridV_calc = logGridV_calc - logGridV_calc(1) - 1;

            if any(abs(sum(trProbM_calc, 2) - 1) > 1e-6)
                 row_sums = sum(trProbM_calc, 2);
                 row_sums(row_sums <= 1e-9) = 1;
                 trProbM_calc = bsxfun(@rdivide, trProbM_calc, row_sums);
            end

            if abs(sum(prob1V_calc) - 1) > 1e-6
                prob1V_calc = prob1V_calc ./ sum(prob1V_calc);
            end

            logGridV = logGridV_calc;
            trProbM = trProbM_calc;
            prob1V = prob1V_calc;
            fprintf('劳动禀赋过程参数已生成 (Tauchen & norm_grid)。\n');
        end

        function [y_grid_out, trProbM_out] = tauchen(N_states, persistence_rho, shock_sigma, mean_val_mu, num_std_dev_width)
            % (Identical to main_olg_baseline_utils.m)
            std_y_unconditional = sqrt(shock_sigma^2 / (1 - persistence_rho^2));
            if abs(1-persistence_rho) < 1e-9
                std_y_unconditional = shock_sigma * 100;
                warning('Tauchen: persistence_rho (%.4f) 接近1，std_y_unconditional可能不准确。', persistence_rho);
            end

            y_max_boundary = num_std_dev_width * std_y_unconditional;
            y_min_boundary = -y_max_boundary;
            y_grid_centered = linspace(y_min_boundary, y_max_boundary, N_states);

            step_size_d = 0;
            if N_states > 1
                step_size_d = y_grid_centered(2) - y_grid_centered(1);
            end

            trProbM_calc = zeros(N_states, N_states);
            for iRow = 1:N_states
                for iCol = 1:N_states
                    mean_next_y_conditional = persistence_rho * y_grid_centered(iRow);
                    if iCol == 1
                        trProbM_calc(iRow,iCol) = normcdf((y_grid_centered(1) - mean_next_y_conditional + step_size_d/2) / shock_sigma);
                    elseif iCol == N_states
                        trProbM_calc(iRow,iCol) = 1 - normcdf((y_grid_centered(N_states) - mean_next_y_conditional - step_size_d/2) / shock_sigma);
                    else
                        upper_bound_cdf = normcdf((y_grid_centered(iCol) - mean_next_y_conditional + step_size_d/2) / shock_sigma);
                        lower_bound_cdf = normcdf((y_grid_centered(iCol) - mean_next_y_conditional - step_size_d/2) / shock_sigma);
                        trProbM_calc(iRow,iCol) = upper_bound_cdf - lower_bound_cdf;
                    end
                end
            end

            row_sums_check = sum(trProbM_calc,2);
            row_sums_check(row_sums_check <= 1e-9) = 1;
            trProbM_out = bsxfun(@rdivide, trProbM_calc, row_sums_check);

            unconditional_mean_shift = mean_val_mu / (1-persistence_rho);
             if abs(1-persistence_rho) < 1e-9 && mean_val_mu ~= 0
                unconditional_mean_shift = sign(mean_val_mu) * inf;
                warning('Tauchen: rho=1 and mu non-zero, unconditional mean is ill-defined.');
            elseif abs(1-persistence_rho) < 1e-9 && mean_val_mu == 0
                unconditional_mean_shift = 0;
            end
            
            y_grid_out = y_grid_centered;
            if isfinite(unconditional_mean_shift)
                y_grid_out = y_grid_centered + unconditional_mean_shift;
            end
        end

        function [massV_out, lbV_out, ubV_out] = norm_grid(x_grid_points, overall_min_bound, overall_max_bound, dist_mean_mu, dist_std_sig)
            % (Identical to main_olg_baseline_utils.m)
            num_points = length(x_grid_points);
            x_grid_points_row = x_grid_points(:)';

            if num_points > 1 && any(x_grid_points_row(2:num_points) < x_grid_points_row(1:num_points-1))
                error('norm_grid: 输入的网格点 x_grid_points 必须单调递增。');
            end

            lower_bounds = [];
            upper_bounds = [];

            if num_points > 1
                mid_points = 0.5 * (x_grid_points_row(1:num_points-1) + x_grid_points_row(2:num_points));
                lower_bounds = [overall_min_bound, mid_points];
                upper_bounds = [mid_points, overall_max_bound];
            elseif num_points == 1
                lower_bounds = overall_min_bound;
                upper_bounds = overall_max_bound;
            else
                massV_out = [];
                lbV_out = [];
                ubV_out = [];
                return;
            end

            cdf_values_at_bounds = normcdf([lower_bounds, upper_bounds(end)], dist_mean_mu, dist_std_sig);
            massV_calc = diff(cdf_values_at_bounds);

            if any(massV_calc < -1e-9)
                warning('norm_grid: 检测到负概率质量。已将其设为0。');
                massV_calc(massV_calc < 0) = 0;
            end

            total_sum_mass = sum(massV_calc);
            if total_sum_mass > 1e-9
                massV_out = massV_calc / total_sum_mass;
            else
                if num_points > 0
                    massV_out = ones(1, num_points) / num_points;
                    warning('norm_grid: 总概率质量过小，已使用均匀分布。');
                else
                    massV_out = [];
                end
            end

            massV_out = massV_out(:);
            lbV_out = lower_bounds(:);
            ubV_out = upper_bounds(:);

            if num_points > 1 && any(ubV_out < lbV_out)
                error('norm_grid: 区间上界小于下界，发生错误。');
            end
        end

        function eIdxM_out = MarkovChainSimulation(num_simulations, num_periods_sim, initial_prob_dist_p0V, transition_matrix_trProbM, random_numbers_rvInM)
            % (Identical to main_olg_baseline_utils.m)
            num_states = length(initial_prob_dist_p0V);

            if size(transition_matrix_trProbM,1) ~= num_states || size(transition_matrix_trProbM,2) ~= num_states
                error('MarkovChainSimulation: 转移矩阵维度与初始分布长度不匹配。');
            end
            if abs(sum(initial_prob_dist_p0V) - 1) > 1e-5
                warning('MarkovChainSimulation: 初始分布 p0V 的和不为1，已重新归一化。');
                initial_prob_dist_p0V = initial_prob_dist_p0V ./ sum(initial_prob_dist_p0V);
            end
             if any(abs(sum(transition_matrix_trProbM, 2) - 1) > 1e-5)
                warning('MarkovChainSimulation: 转移矩阵 trProbM 的某些行和不为1，已重新归一化。');
                row_sums_tr = sum(transition_matrix_trProbM, 2);
                row_sums_tr(row_sums_tr <= 1e-9) = 1;
                transition_matrix_trProbM = bsxfun(@rdivide, transition_matrix_trProbM, row_sums_tr);
            end
            if size(random_numbers_rvInM,1) ~= num_simulations || size(random_numbers_rvInM,2) ~= num_periods_sim
                error('MarkovChainSimulation: 随机数矩阵维度与模拟参数不匹配。');
            end

            cumulative_initial_prob_cP0 = cumsum(initial_prob_dist_p0V(:)');
            cumulative_transition_prob_cPT = cumsum(transition_matrix_trProbM, 2);
            if num_states > 0
                cumulative_transition_prob_cPT(:, num_states) = 1.0;
            end

            eIdxM_out = zeros(num_simulations, num_periods_sim, 'uint16');
            if num_simulations > 0 && num_periods_sim > 0
                eIdxM_out(:, 1) = 1 + sum(bsxfun(@gt, random_numbers_rvInM(:,1), cumulative_initial_prob_cP0), 2);
            end

            for t_mc_loop = 1:(num_periods_sim - 1)
                current_state_indices_cSI = eIdxM_out(:, t_mc_loop);
                valid_indices_vPI = (current_state_indices_cSI >= 1) & (current_state_indices_cSI <= num_states);
                if ~all(valid_indices_vPI)
                    warning('MarkovChainSimulation: 在期 %d 检测到无效的当前状态索引。已重置为状态1。', t_mc_loop);
                    current_state_indices_cSI(~valid_indices_vPI) = 1;
                    eIdxM_out(:, t_mc_loop) = current_state_indices_cSI;
                end
                cPt_for_next_state = cumulative_transition_prob_cPT(current_state_indices_cSI, :);
                eIdxM_out(:, t_mc_loop + 1) = 1 + sum(bsxfun(@gt, random_numbers_rvInM(:, t_mc_loop + 1), cPt_for_next_state), 2);
            end
        end

        function eIdxM = LaborEndowSimulation_olgm(cS, paramS)
            % (Identical to main_olg_baseline_utils.m)
            rng(433);
            random_numbers_for_sim = rand([cS.nSim, cS.aD_orig]);
            eIdxM = main_olg_v8_utils.MarkovChainSimulation(cS.nSim, cS.aD_orig, ...
                                                          paramS.leProb1V, paramS.leTrProbM, ...
                                                          random_numbers_for_sim);
            fprintf('劳动禀赋路径已模拟 (%d 个体, %d 年度年龄)。\n', cS.nSim, cS.aD_orig);
        end

        function [HHlaborM_group, L_total_eff_pc] = LaborSupply_Huggett(eIdxM_annual, cS, paramS, Z_ss_norm_group)
            % (Identical to main_olg_baseline_utils.m)
            nSim = size(eIdxM_annual, 1);
            HHlaborM_group = zeros(nSim, cS.aD_new);
            ageToGroupMap_local = zeros(cS.aD_orig, 1);
            for a_new_map_idx = 1:cS.aD_new
                annual_indices_in_group = cS.physAgeMap{a_new_map_idx};
                if ~isempty(annual_indices_in_group)
                    ageToGroupMap_local(annual_indices_in_group) = a_new_map_idx;
                end
            end
            leGridV_col_local = paramS.leGridV(:);
            HHlaborM_annual_temp = zeros(nSim, cS.aD_orig);
            for a_orig_idx = 1 : cS.aD_orig
               if a_orig_idx < cS.aR_idx_orig
                   a_new_group_idx_current = ageToGroupMap_local(a_orig_idx);
                   if a_new_group_idx_current > 0 && a_new_group_idx_current <= cS.aR_new
                       HHlaborM_annual_temp(:, a_orig_idx) = cS.ageEffV_new(a_new_group_idx_current) .* leGridV_col_local(eIdxM_annual(:, a_orig_idx));
                   end
               end
            end
            for a_new_group_idx = 1:cS.aD_new
                annual_indices_for_this_group = cS.physAgeMap{a_new_group_idx};
                if ~isempty(annual_indices_for_this_group)
                    HHlaborM_group(:, a_new_group_idx) = mean(HHlaborM_annual_temp(:, annual_indices_for_this_group), 2);
                end
            end
            L_total_eff_pc_sum = 0;
            if cS.aR_new > 0
                mean_labor_per_working_group = mean(HHlaborM_group(:, 1:cS.aR_new), 1);
                L_total_eff_pc_sum = mean_labor_per_working_group * Z_ss_norm_group(1:cS.aR_new);
            end
            L_total_eff_pc = max(0, L_total_eff_pc_sum);
            fprintf('家庭劳动供给已计算。总体人均有效劳动供给 (L_eff_pc) = %.4f\n', L_total_eff_pc);
        end

        % =====================================================================
        % == 宏观经济与要素价格函数 ==
        % =====================================================================
        % 对应v8.tex第2.3节"生产部门"中的要素价格决定机制
        
        function [R_market_gross_factor, MPL_gross] = HHPrices_Huggett(K_productive, L_total_eff, cS)
            % HHPrices_Huggett - 根据边际生产力计算要素价格
            % 
            % 输入：
            %   K_productive - 总生产性资本（包含非PPS资本和PPS资本）
            %   L_total_eff - 总有效劳动供给
            %   cS - 模型参数
            % 输出：
            %   R_market_gross_factor - 市场毛回报因子(1+r_gross)
            %   MPL_gross - 劳动边际产品（毛工资率）
            % 
            % 功能：实现v8.tex公式：
            %   Y = A·K^α·L^(1-α)
            %   r_gross = α·(Y/K) - δ
            %   w_gross = (1-α)·(Y/L)
            if K_productive <= 0
                K_productive=1e-6;
                warning('HHPrices: K_productive 非正，已重置为 %.1e。', K_productive);
            end
            if L_total_eff <= 0
                L_total_eff=1e-6;
                warning('HHPrices: L_total_eff 非正，已重置为 %.1e。', L_total_eff);
            end
            Y_gross = cS.A * (K_productive^cS.alpha) * (L_total_eff^(1-cS.alpha));
            MPK_gross_val = cS.alpha * Y_gross / K_productive;
            MPL_gross = (1-cS.alpha) * Y_gross / L_total_eff;
            R_market_gross_factor = 1 + MPK_gross_val - cS.ddk;
            R_market_gross_factor = max(1.0 + 1e-6, R_market_gross_factor);
        end

        % =====================================================================
        % == 家庭效用函数 ==
        % =====================================================================
        % 对应v8.tex第2.2.1节"偏好与效用"中的CRRA效用函数
        
        function [muM, utilM] = CES_utility(cM_quantity, sigma_crra, cS_common)
            % CES_utility - 计算CRRA效用函数和边际效用
            % 
            % 输入：
            %   cM_quantity - 消费量（可为向量）
            %   sigma_crra - 相对风险厌恶系数γ
            %   cS_common - 包含最低消费约束的参数结构体
            % 输出：
            %   muM - 边际效用 ∂U/∂c
            %   utilM - 效用水平 U(c)
            % 
            % 功能：实现U(c) = (c^(1-γ))/(1-γ)，当γ=1时为对数效用
            if ~isscalar(sigma_crra) || sigma_crra <= 0
                error('CES_utility: sigma_crra 必须是正标量。');
            end
            min_c_quantity = cS_common.cFloor;
            is_valid_consumption = (cM_quantity >= min_c_quantity);
            c_adjusted_quantity = max(min_c_quantity, cM_quantity);
            utilM = -Inf(size(cM_quantity));
            muM   =  Inf(size(cM_quantity));
            if abs(sigma_crra - 1) < 1e-6
                utilM(is_valid_consumption) = log(c_adjusted_quantity(is_valid_consumption));
                muM(is_valid_consumption)   = 1 ./ c_adjusted_quantity(is_valid_consumption);
            else
                utilM(is_valid_consumption) = (c_adjusted_quantity(is_valid_consumption).^(1-sigma_crra)) ./ (1-sigma_crra);
                muM(is_valid_consumption)   = c_adjusted_quantity(is_valid_consumption).^(-sigma_crra);
            end
            utilM(~is_valid_consumption) = -1e10 - (min_c_quantity - cM_quantity(~is_valid_consumption))*1e10;
            if abs(sigma_crra - 1) < 1e-6
                 muM(~is_valid_consumption) = 1 ./ c_adjusted_quantity(~is_valid_consumption) + 1e10;
            else
                 muM(~is_valid_consumption) = c_adjusted_quantity(~is_valid_consumption).^(-sigma_crra) + 1e10;
            end
        end

        % =====================================================================
        % == V8 家庭问题核心函数：内生PPS缴费决策 ==
        % =====================================================================
        % 实现v8.tex第2.2.2节"预算约束与决策"中描述的家庭优化问题

        % --- 家庭收入计算函数 ---
        function [resources_for_c_and_k_prime, labor_income_gross_state, pps_deduction_actual_state] = HHIncome_Huggett(...
                k_now_val, R_k_net_factor, w_gross, ...
                TR_total, b_payg_val, c_pps_chosen_and_constrained_input_val, ...
                a_idx, paramS_hh, cS, epsilon_val)
            % HHIncome_Huggett - 计算家庭可支配资源
            % 
            % 输入：
            %   k_now_val - 当前非PPS资产
            %   R_k_net_factor - 税后资本回报因子(1+r_net)
            %   w_gross - 市场毛工资率
            %   TR_total - 总转移支付
            %   b_payg_val - PAYG养老金
            %   c_pps_chosen_and_constrained_input_val - 已选择的PPS缴费额
            %   a_idx - 年龄组索引
            %   paramS_hh - 家庭参数
            %   cS - 模型参数
            %   epsilon_val - 当前劳动效率
            % 输出：
            %   resources_for_c_and_k_prime - 可用于消费和储蓄的资源
            %   labor_income_gross_state - 税前劳动收入
            %   pps_deduction_actual_state - PPS税前扣除额
            % 
            % 功能：实现v8.tex中的预算约束计算，包括PPS税收递延处理

            labor_income_gross_state = 0;
            pps_deduction_actual_state = 0;
            non_capital_income = 0;
            % actual_pps_contribution_expenditure IS c_pps_chosen_and_constrained_input_val
            % No further rule-based calculation here, it's already determined upstream.

            if a_idx <= cS.aR_new % 如果是工作年龄组
                age_efficiency = cS.ageEffV_new(a_idx);
                labor_income_gross_state = w_gross * age_efficiency * epsilon_val;

                % PPS缴费 c_pps_chosen_and_constrained_input_val 已由VFI外层优化选择并施加约束
                actual_pps_contribution_expenditure = max(0, c_pps_chosen_and_constrained_input_val); % Ensure non-negative

                if isfield(paramS_hh, 'pps_tax_deferral_active') && paramS_hh.pps_tax_deferral_active
                    pps_deduction_actual_state = actual_pps_contribution_expenditure;
                else
                    pps_deduction_actual_state = 0;
                end

                labor_income_taxable_for_tau_l = labor_income_gross_state - pps_deduction_actual_state;
                labor_income_taxable_for_tau_l = max(0, labor_income_taxable_for_tau_l);

                income_tax_tau_l = labor_income_taxable_for_tau_l * paramS_hh.tau_l;
                payg_tax_theta = labor_income_gross_state * paramS_hh.theta_payg_actual_for_hh;

                labor_income_net_of_all_taxes = labor_income_gross_state - income_tax_tau_l - payg_tax_theta;
                non_capital_income = labor_income_net_of_all_taxes + TR_total + b_payg_val;
            else % 如果是退休年龄组
                actual_pps_contribution_expenditure = 0; % 退休期PPS缴费为0
                pps_deduction_actual_state = 0;
                non_capital_income = TR_total + b_payg_val;
            end

            capital_income_net_of_tax = (R_k_net_factor - 1) * k_now_val;
            resources_for_c_and_k_prime = k_now_val + capital_income_net_of_tax + non_capital_income - actual_pps_contribution_expenditure;

            if ~isfinite(resources_for_c_and_k_prime)
                resources_for_c_and_k_prime = -1e10;
                warning('HHIncome_Huggett (V8): 计算得到的资源为非有限值。');
            end
        end

        % --- V8 VFI 主函数：家庭生命周期优化求解器 ---
        function [cPolM_q, kPolM, cPpsPolM_choice, valM] = HHSolution_VFI_Huggett(...
                R_k_net_factor_vfi, w_gross_vfi, TR_total_vfi, bV_payg_vfi, paramS_vfi, cS)
            % HHSolution_VFI_Huggett - V8模型的值函数迭代主求解器
            % 
            % 输入：
            %   R_k_net_factor_vfi - 税后资本回报因子
            %   w_gross_vfi - 市场毛工资率
            %   TR_total_vfi - 总转移支付
            %   bV_payg_vfi - 各年龄组PAYG福利向量
            %   paramS_vfi - 家庭决策参数
            %   cS - 模型参数
            % 输出 (均为 4D 矩阵: nk x nkpps x nw x aD_new):
            %   cPolM_q - 最优消费策略
            %   kPolM - 最优非PPS储蓄策略
            %   cPpsPolM_choice - 最优PPS缴费选择策略 (V8核心创新)
            %   valM - 值函数
            % 
            % 功能：通过逆向迭代求解v8.tex第2.2.3节描述的动态规划问题
            %       实现内生PPS缴费选择的值函数迭代算法

            cPolM_q  = zeros(cS.nk, cS.nkpps, cS.nw, cS.aD_new);
            kPolM  = zeros(cS.nk, cS.nkpps, cS.nw, cS.aD_new);
            cPpsPolM_choice = zeros(cS.nk, cS.nkpps, cS.nw, cS.aD_new); % V8: Optimal PPS contribution
            valM = -Inf(cS.nk, cS.nkpps, cS.nw, cS.aD_new);

            % fprintf('  VFI V8 (HHSolution_VFI_Huggett): 开始逆向迭代...\n');
            total_age_groups_vfi = cS.aD_new;
            for a_idx = cS.aD_new : -1 : 1
                % if mod(cS.aD_new - a_idx + 1, 1) == 0 || a_idx == cS.aD_new || a_idx == 1
                %     progress_pct_vfi = (cS.aD_new - a_idx + 1) / total_age_groups_vfi * 100;
                %     fprintf('    VFI V8 年龄组 %2d / %2d (大约 %.0f%%)...\n', a_idx, total_age_groups_vfi, progress_pct_vfi);
                % end

                vPrime_kkppse_next = [];
                if a_idx < cS.aD_new
                    vPrime_kkppse_next = valM(:,:,:,a_idx+1);
                end
                eps_grid_for_vfi = paramS_vfi.leGridV;

                % 调用按年龄组求解的函数 (V8 version)
                [cPolM_q(:,:,:,a_idx), kPolM(:,:,:,a_idx), cPpsPolM_choice(:,:,:,a_idx), valM(:,:,:,a_idx)] = ...
                    main_olg_v8_utils.HHSolutionByAge_VFI_Huggett_v8(a_idx, vPrime_kkppse_next, ...
                    R_k_net_factor_vfi, w_gross_vfi, TR_total_vfi, bV_payg_vfi(a_idx), paramS_vfi, cS, eps_grid_for_vfi);
            end
            % fprintf('  VFI V8 (HHSolution_VFI_Huggett): 完成。\n');
        end

        % --- V8 VFI 核心算法：按年龄组求解内生PPS缴费选择 ---
        function [cPol_age_q, kPol_age, cPpsPol_age_choice, val_age] = HHSolutionByAge_VFI_Huggett_v8(...
                a_idx, vPrime_kkppse_next, R_k_net_factor_age, w_gross_age, TR_total_age, b_age_val, ...
                paramS_age, cS, epsilon_grid)
            % HHSolutionByAge_VFI_Huggett_v8 - V8模型的核心创新：内生PPS缴费决策VFI
            % 
            % 输入：
            %   a_idx - 当前年龄组索引
            %   vPrime_kkppse_next - 下一期值函数V_{a+1}(k',k'_pps,ε')
            %   R_k_net_factor_age - 税后资本回报因子
            %   w_gross_age - 市场毛工资率
            %   TR_total_age - 转移支付
            %   b_age_val - 当期PAYG福利
            %   paramS_age - 家庭决策参数
            %   cS - 模型参数
            %   epsilon_grid - 劳动效率网格
            % 输出：
            %   cPol_age_q - 该年龄组的最优消费策略
            %   kPol_age - 该年龄组的最优非PPS储蓄策略
            %   cPpsPol_age_choice - 该年龄组的最优PPS缴费策略（V8关键）
            %   val_age - 该年龄组的值函数
            %
            % 算法逻辑（实现v8.tex第3节描述的求解算法）：
            %   对于给定年龄组a和所有状态(k, k_pps, ε)：
            %   1. 离散化PPS缴费选择空间：构建从0到最大允许缴费额的网格
            %   2. 对每个可能的PPS缴费额c_pps_choice：
            %      a. 计算约束后的实际PPS缴费额
            %      b. 计算下一期PPS资产k'_pps
            %      c. 在给定c_pps下优化非PPS储蓄k'（内层优化）
            %      d. 计算该选择下的总价值（当期效用+贴现未来期望效用）
            %   3. 选择使总价值最大的c_pps_choice（外层优化）
            %   
            % 这种双层优化结构实现了v8.tex中的联合最优化：
            %   V_a(k,k_pps,ε) = max_{c_pps} max_{k'} [U(C) + β E[V_{a+1}(k',k'_pps,ε')]]

            cPol_age_q_init = zeros(cS.nk, cS.nkpps, cS.nw);
            kPol_age_init   = zeros(cS.nk, cS.nkpps, cS.nw);
            cPpsPol_age_choice_init = zeros(cS.nk, cS.nkpps, cS.nw); % V8: Optimal PPS contribution choice
            val_age_init    = -Inf(cS.nk, cS.nkpps, cS.nw);

            fminbnd_opts_iter = optimset('TolX', cS.fminbnd_TolX, 'Display', cS.fminbnd_Display);

            % --- 情况1: 最后一个年龄组 (矢量化计算) ---
            if a_idx == cS.aD_new
                % 矢量化最后期计算
                [K_grid, Kpps_grid, ~] = ndgrid(cS.kGridV, cS.kppsGridV, epsilon_grid);
                
                % 批量计算收入
                resources_batch = zeros(size(K_grid));
                for ie = 1:cS.nw
                    for ik = 1:cS.nk
                        for ikpps = 1:cS.nkpps
                            [resources_batch(ik, ikpps, ie), ~, ~] = main_olg_v8_utils.HHIncome_Huggett(...
                                K_grid(ik, ikpps, ie), R_k_net_factor_age, w_gross_age, TR_total_age, b_age_val, ...
                                0, a_idx, paramS_age, cS, epsilon_grid(ie));
                        end
                    end
                end
                
                % 加上PPS最终提取
                total_resources = resources_batch + Kpps_grid * (1 - cS.pps_tax_rate_withdrawal);
                
                % 矢量化消费和效用计算
                cPol_age_q_init = max(cS.cFloor, total_resources / (1 + cS.tau_c));
                kPol_age_init(:) = cS.kMin;
                cPpsPol_age_choice_init(:) = 0;
                
                % 批量计算效用
                for idx = 1:numel(cPol_age_q_init)
                    [~, val_age_init(idx)] = main_olg_v8_utils.CES_utility(cPol_age_q_init(idx), cS.sigma, cS);
                end
                
            else % --- 情况2: 非最后一个年龄组 (优化版) ---
                % 预计算期望值矩阵 (矢量化)
                EV_matrix = zeros(cS.nk, cS.nkpps, cS.nw);
                for ie_current = 1:cS.nw
                    % 矢量化期望值计算
                    transition_probs = paramS_age.leTrProbM(ie_current, :);
                    EV_slice = sum(vPrime_kkppse_next .* reshape(transition_probs, 1, 1, cS.nw), 3);
                    EV_matrix(:, :, ie_current) = EV_slice;
                end
                
                % 创建插值器
                EV_interpolants = cell(cS.nw, 1);
                for ie_current = 1:cS.nw
                    if cS.nk > 1 && cS.nkpps > 1
                        EV_interpolants{ie_current} = griddedInterpolant({cS.kGridV, cS.kppsGridV}, ...
                            EV_matrix(:, :, ie_current), 'linear', 'linear');
                    elseif cS.nk > 1
                        EV_interpolants{ie_current} = griddedInterpolant(cS.kGridV, ...
                            EV_matrix(:, 1, ie_current), 'linear', 'linear');
                    elseif cS.nkpps > 1
                        EV_interpolants{ie_current} = griddedInterpolant(cS.kppsGridV, ...
                            EV_matrix(1, :, ie_current)', 'linear', 'linear');
                    else
                        EV_interpolants{ie_current} = @(k_s, kp_s) EV_matrix(1, 1, ie_current);
                    end
                end
                
                % 预计算所有状态组合
                total_states = cS.nk * cS.nkpps * cS.nw;
                state_indices = cell(total_states, 1);
                state_params = cell(total_states, 1);
                
                state_idx = 1;
                for ik = 1:cS.nk
                    for ikpps = 1:cS.nkpps
                        for ie = 1:cS.nw
                            state_indices{state_idx} = [ik, ikpps, ie];
                            
                            % 预计算该状态的参数
                            k_state = cS.kGridV(ik);
                            k_pps_state = cS.kppsGridV(ikpps);
                            epsilon_state = epsilon_grid(ie);
                            
                            % 预计算PPS缴费资格和收入
                            model_age_group_start_year_idx = cS.physAgeMap{a_idx}(1);
                            is_pps_eligible = (a_idx <= cS.aR_new && ...
                                model_age_group_start_year_idx <= cS.pps_contribution_age_max_idx && ...
                                cS.pps_active);
                            
                            current_gross_labor_income = 0;
                            if is_pps_eligible
                                age_efficiency = cS.ageEffV_new(a_idx);
                                current_gross_labor_income = w_gross_age * age_efficiency * epsilon_state;
                            end
                            
                            % 预计算PPS缴费网格
                            cpps_grid = [0]; % 总是允许零缴费
                            if is_pps_eligible && current_gross_labor_income > 1e-6
                                max_cpps_by_frac = current_gross_labor_income * cS.pps_max_contrib_frac;
                                max_permissible_cpps = min(cS.pps_annual_contrib_limit, max_cpps_by_frac);
                                max_permissible_cpps = max(0, max_permissible_cpps);
                                
                                if max_permissible_cpps > 1e-6 && cS.n_pps_choice_grid_points > 1
                                    additional_choices = linspace(1e-6, max_permissible_cpps, ...
                                        cS.n_pps_choice_grid_points - 1);
                                    cpps_grid = unique([0, additional_choices]);
                                elseif max_permissible_cpps > 1e-6
                                    cpps_grid = unique([0, max_permissible_cpps]);
                                end
                            end
                            
                            state_params{state_idx} = struct(...
                                'k_state', k_state, ...
                                'k_pps_state', k_pps_state, ...
                                'epsilon_state', epsilon_state, ...
                                'ie', ie, ...
                                'cpps_grid', cpps_grid, ...
                                'is_pps_eligible', is_pps_eligible, ...
                                'current_gross_labor_income', current_gross_labor_income);
                            
                            state_idx = state_idx + 1;
                        end
                    end
                end
                
                % 多维并行计算
                results_cell = cell(total_states, 1);
                
                parfor state_idx_par = 1:total_states
                    % 获取状态参数
                    indices = state_indices{state_idx_par};
                    params = state_params{state_idx_par};
                    
                    % 本地化变量以避免broadcast开销
                    cS_local = cS;
                    paramS_age_local = paramS_age;
                    a_idx_local = a_idx;
                    R_k_net_factor_age_local = R_k_net_factor_age;
                    w_gross_age_local = w_gross_age;
                    TR_total_age_local = TR_total_age;
                    b_age_val_local = b_age_val;
                    fminbnd_opts_local = fminbnd_opts_iter;
                    EV_interpolants_local = EV_interpolants;
                    
                    % 状态优化
                    max_value = -Inf;
                    optimal_k_prime = cS_local.kMin;
                    optimal_c = cS_local.cFloor;
                    optimal_cpps = 0;
                    
                    % 对该状态的所有PPS缴费选择进行优化
                    for i_cpps = 1:length(params.cpps_grid)
                        c_pps_try = params.cpps_grid(i_cpps);
                        
                        % 计算k_pps_prime
                        pps_withdrawal = 0;
                        annual_age_check = cS_local.physAgeMap{a_idx_local}(1);
                        is_retired = (a_idx_local > cS_local.aR_new);
                        if is_retired && annual_age_check >= cS_local.pps_withdrawal_age_min_idx && cS_local.pps_active
                            pps_withdrawal = params.k_pps_state * cS_local.pps_withdrawal_rate;
                        end
                        
                        pps_return_factor = 1 + ((R_k_net_factor_age_local - 1) + cS_local.pps_return_rate_premium);
                        k_pps_prime = (params.k_pps_state + c_pps_try - pps_withdrawal) * pps_return_factor;
                        k_pps_prime = max(cS_local.kppsMin, min(cS_local.kppsMax, k_pps_prime));
                        
                        % 计算可用于c和k_prime的资源
                        [resources_for_c_k_prime, ~, ~] = main_olg_v8_utils.HHIncome_Huggett(...
                            params.k_state, R_k_net_factor_age_local, w_gross_age_local, ...
                            TR_total_age_local, b_age_val_local, c_pps_try, ...
                            a_idx_local, paramS_age_local, cS_local, params.epsilon_state);
                        
                        % 创建1D插值器优化k_prime
                        ev_func_1d = @(k_p) main_olg_v8_utils.CallInterpolator(...
                            EV_interpolants_local{params.ie}, k_p, k_pps_prime, cS_local);
                        
                        % 优化k_prime
                        [c_val, k_p_val, v_val] = main_olg_v8_utils.HHSolutionByOneState_OptK_Mod(...
                            a_idx_local, resources_for_c_k_prime, ev_func_1d, ...
                            fminbnd_opts_local, cS_local, paramS_age_local);
                        
                        % 更新最优解
                        if v_val > max_value
                            max_value = v_val;
                            optimal_k_prime = k_p_val;
                            optimal_c = c_val;
                            optimal_cpps = c_pps_try;
                        end
                    end
                    
                    % 存储结果
                    results_cell{state_idx_par} = struct(...
                        'indices', indices, ...
                        'optimal_c', optimal_c, ...
                        'optimal_k_prime', optimal_k_prime, ...
                        'optimal_cpps', optimal_cpps, ...
                        'max_value', max_value);
                end
                
                % 组装结果
                for state_idx = 1:total_states
                    result = results_cell{state_idx};
                    ik = result.indices(1);
                    ikpps = result.indices(2);
                    ie = result.indices(3);
                    
                    cPol_age_q_init(ik, ikpps, ie) = result.optimal_c;
                    kPol_age_init(ik, ikpps, ie) = result.optimal_k_prime;
                    cPpsPol_age_choice_init(ik, ikpps, ie) = result.optimal_cpps;
                    val_age_init(ik, ikpps, ie) = result.max_value;
                end
            end
            
            cPol_age_q = cPol_age_q_init;
            kPol_age = kPol_age_init;
            cPpsPol_age_choice = cPpsPol_age_choice_init;
            val_age = val_age_init;
        end

        % --- V8 VFI 优化版算法：多维并行化加速版本 ---
        function [cPol_age_q, kPol_age, cPpsPol_age_choice, val_age] = HHSolutionByAge_VFI_Huggett_v8_optimized(...
                a_idx, vPrime_kkppse_next, R_k_net_factor_age, w_gross_age, TR_total_age, b_age_val, ...
                paramS_age, cS, epsilon_grid)
            % HHSolutionByAge_VFI_Huggett_v8_optimized - V8模型的优化版VFI算法
            % 
            % 主要优化策略：
            % 1. 多维并行化：将3层嵌套循环展开为单层并行计算
            % 2. 预计算优化：避免重复计算状态参数和PPS网格
            % 3. 矢量化加速：使用MATLAB内置矩阵运算
            % 4. 内存访问优化：减少数据传输开销
            
            cPol_age_q_init = zeros(cS.nk, cS.nkpps, cS.nw);
            kPol_age_init   = zeros(cS.nk, cS.nkpps, cS.nw);
            cPpsPol_age_choice_init = zeros(cS.nk, cS.nkpps, cS.nw);
            val_age_init    = -Inf(cS.nk, cS.nkpps, cS.nw);

            fminbnd_opts_iter = optimset('TolX', cS.fminbnd_TolX, 'Display', cS.fminbnd_Display);

            % --- 情况1: 最后一个年龄组 (全矢量化) ---
            if a_idx == cS.aD_new
                % 预计算所有状态组合
                [K_mesh, Kpps_mesh, E_mesh] = ndgrid(cS.kGridV, cS.kppsGridV, epsilon_grid);
                
                % 批量计算所有状态的收入
                resources_all = zeros(size(K_mesh));
                parfor idx = 1:numel(K_mesh)
                    [i, j, k] = ind2sub(size(K_mesh), idx);
                    [resources_all(idx), ~, ~] = main_olg_v8_utils.HHIncome_Huggett(...
                        K_mesh(idx), R_k_net_factor_age, w_gross_age, TR_total_age, b_age_val, ...
                        0, a_idx, paramS_age, cS, E_mesh(idx));
                end
                
                % 加上PPS最终提取
                total_resources_all = resources_all + Kpps_mesh * (1 - cS.pps_tax_rate_withdrawal);
                
                % 矢量化消费和效用计算
                cPol_age_q_init = max(cS.cFloor, total_resources_all / (1 + cS.tau_c));
                kPol_age_init(:) = cS.kMin;
                cPpsPol_age_choice_init(:) = 0;
                
                % 批量计算效用
                parfor idx = 1:numel(cPol_age_q_init)
                    [~, val_age_init(idx)] = main_olg_v8_utils.CES_utility(cPol_age_q_init(idx), cS.sigma, cS);
                end
                
            else % --- 情况2: 非最后年龄组 (优化版并行计算) ---
                
                % 预计算期望值矩阵 (矢量化)
                EV_matrix = zeros(cS.nk, cS.nkpps, cS.nw);
                for ie_current = 1:cS.nw
                    transition_probs = paramS_age.leTrProbM(ie_current, :);
                    for ie_next = 1:cS.nw
                        EV_matrix(:, :, ie_current) = EV_matrix(:, :, ie_current) + ...
                            transition_probs(ie_next) * vPrime_kkppse_next(:, :, ie_next);
                    end
                end
                
                % 创建快速插值器
                EV_interpolants = cell(cS.nw, 1);
                for ie_current = 1:cS.nw
                    if cS.nk > 1 && cS.nkpps > 1
                        EV_interpolants{ie_current} = griddedInterpolant({cS.kGridV, cS.kppsGridV}, ...
                            EV_matrix(:, :, ie_current), 'linear', 'linear');
                    elseif cS.nk > 1
                        EV_interpolants{ie_current} = griddedInterpolant(cS.kGridV, ...
                            EV_matrix(:, 1, ie_current), 'linear', 'linear');
                    elseif cS.nkpps > 1
                        EV_interpolants{ie_current} = griddedInterpolant(cS.kppsGridV, ...
                            EV_matrix(1, :, ie_current)', 'linear', 'linear');
                    else
                        EV_interpolants{ie_current} = @(k_s, kp_s) EV_matrix(1, 1, ie_current);
                    end
                end
                
                % 展开状态空间为单层并行计算
                total_states = cS.nk * cS.nkpps * cS.nw;
                state_list = zeros(total_states, 6); % [ik, ikpps, ie, k_val, kpps_val, eps_val]
                
                idx = 1;
                for ik = 1:cS.nk
                    for ikpps = 1:cS.nkpps
                        for ie = 1:cS.nw
                            state_list(idx, :) = [ik, ikpps, ie, ...
                                cS.kGridV(ik), cS.kppsGridV(ikpps), epsilon_grid(ie)];
                            idx = idx + 1;
                        end
                    end
                end
                
                % 多维并行优化
                results = zeros(total_states, 4); % [optimal_c, optimal_k_prime, optimal_cpps, max_value]
                
                parfor state_idx = 1:total_states
                    % 获取状态信息
                    ik = state_list(state_idx, 1);
                    ikpps = state_list(state_idx, 2);
                    ie = state_list(state_idx, 3);
                    k_state = state_list(state_idx, 4);
                    k_pps_state = state_list(state_idx, 5);
                    epsilon_state = state_list(state_idx, 6);
                    
                    % 本地化参数以减少通信开销
                    cS_local = cS;
                    paramS_local = paramS_age;
                    a_idx_local = a_idx;
                    
                    % 预计算PPS缴费选择网格
                    cpps_grid = [0]; % 总是允许零缴费
                    
                    % 检查PPS缴费资格
                    model_age_group_start_year_idx = cS_local.physAgeMap{a_idx_local}(1);
                    is_pps_eligible = (a_idx_local <= cS_local.aR_new && ...
                        model_age_group_start_year_idx <= cS_local.pps_contribution_age_max_idx && ...
                        cS_local.pps_active);
                    
                    if is_pps_eligible
                        age_efficiency = cS_local.ageEffV_new(a_idx_local);
                        current_gross_labor_income = w_gross_age * age_efficiency * epsilon_state;
                        
                        if current_gross_labor_income > 1e-6
                            max_cpps_by_frac = current_gross_labor_income * cS_local.pps_max_contrib_frac;
                            max_permissible_cpps = min(cS_local.pps_annual_contrib_limit, max_cpps_by_frac);
                            max_permissible_cpps = max(0, max_permissible_cpps);
                            
                            if max_permissible_cpps > 1e-6 && cS_local.n_pps_choice_grid_points > 1
                                additional_choices = linspace(1e-6, max_permissible_cpps, ...
                                    cS_local.n_pps_choice_grid_points - 1);
                                cpps_grid = unique([0, additional_choices]);
                            elseif max_permissible_cpps > 1e-6
                                cpps_grid = unique([0, max_permissible_cpps]);
                            end
                        end
                    end
                    
                    % 对所有PPS选择进行优化
                    max_value = -Inf;
                    optimal_k_prime = cS_local.kMin;
                    optimal_c = cS_local.cFloor;
                    optimal_cpps = 0;
                    
                    for i_cpps = 1:length(cpps_grid)
                        c_pps_try = cpps_grid(i_cpps);
                        
                        % 计算k_pps_prime
                        pps_withdrawal = 0;
                        annual_age_check = cS_local.physAgeMap{a_idx_local}(1);
                        is_retired = (a_idx_local > cS_local.aR_new);
                        if is_retired && annual_age_check >= cS_local.pps_withdrawal_age_min_idx && cS_local.pps_active
                            pps_withdrawal = k_pps_state * cS_local.pps_withdrawal_rate;
                        end
                        
                        pps_return_factor = 1 + ((R_k_net_factor_age - 1) + cS_local.pps_return_rate_premium);
                        k_pps_prime = (k_pps_state + c_pps_try - pps_withdrawal) * pps_return_factor;
                        k_pps_prime = max(cS_local.kppsMin, min(cS_local.kppsMax, k_pps_prime));
                        
                        % 计算可用资源
                        [resources_for_c_k_prime, ~, ~] = main_olg_v8_utils.HHIncome_Huggett(...
                            k_state, R_k_net_factor_age, w_gross_age, TR_total_age, b_age_val, ...
                            c_pps_try, a_idx_local, paramS_local, cS_local, epsilon_state);
                        
                        % 1D优化k_prime
                        ev_func_1d = @(k_p) main_olg_v8_utils.CallInterpolator(...
                            EV_interpolants{ie}, k_p, k_pps_prime, cS_local);
                        
                        [c_val, k_p_val, v_val] = main_olg_v8_utils.HHSolutionByOneState_OptK_Mod(...
                            a_idx_local, resources_for_c_k_prime, ev_func_1d, ...
                            fminbnd_opts_iter, cS_local, paramS_local);
                        
                        % 更新最优解
                        if v_val > max_value
                            max_value = v_val;
                            optimal_k_prime = k_p_val;
                            optimal_c = c_val;
                            optimal_cpps = c_pps_try;
                        end
                    end
                    
                    % 存储结果
                    results(state_idx, :) = [optimal_c, optimal_k_prime, optimal_cpps, max_value];
                end
                
                % 重新组装结果到3D矩阵
                for state_idx = 1:total_states
                    ik = state_list(state_idx, 1);
                    ikpps = state_list(state_idx, 2);
                    ie = state_list(state_idx, 3);
                    
                    cPol_age_q_init(ik, ikpps, ie) = results(state_idx, 1);
                    kPol_age_init(ik, ikpps, ie) = results(state_idx, 2);
                    cPpsPol_age_choice_init(ik, ikpps, ie) = results(state_idx, 3);
                    val_age_init(ik, ikpps, ie) = results(state_idx, 4);
                end
            end
            
            % 返回结果
            cPol_age_q = cPol_age_q_init;
            kPol_age = kPol_age_init;
            cPpsPol_age_choice = cPpsPol_age_choice_init;
            val_age = val_age_init;
        end

        % --- HHSolutionByOneState_OptK_Mod (与V7/Baseline一致) ---
        function [c_quantity, kPrime_out, ValueFunc_out] = HHSolutionByOneState_OptK_Mod(...
                a_idx_current, budget_for_c_expend_and_kprime, EVprime_of_kprime_interp_1D, ...
                fminbnd_opts_in, cS, paramS_current_age) % paramS_current_age unused
            % (Identical to main_olg_baseline_utils.m)
            kPrime_min_bound = cS.kMin;
            kPrime_max_bound = budget_for_c_expend_and_kprime - (cS.cFloor * (1 + cS.tau_c));
            ValueFunc_out = -Inf;
            kPrime_out = cS.kMin; 
            c_quantity = cS.cFloor; 

            function negV_nested = negBellmanObjective_nested_k_only(kPrime_choice_nested)
                [negV_nested, ~] = BellmanInner_nested_k_only(kPrime_choice_nested);
            end

            function [negVal_out_nested, Val_out_nested] = BellmanInner_nested_k_only(kPrime_inner_nested)
                consumption_expenditure_nested = budget_for_c_expend_and_kprime - kPrime_inner_nested;
                current_c_quantity_nested = max(cS.cFloor, consumption_expenditure_nested / (1 + cS.tau_c) );
                [~, util_current_period_nested] = main_olg_v8_utils.CES_utility(current_c_quantity_nested, cS.sigma, cS);
                if ~isfinite(util_current_period_nested)
                    util_current_period_nested = -1e12;
                end
                Val_out_nested = util_current_period_nested;
                if a_idx_current < cS.aD_new
                    expected_future_value_nested = -Inf;
                    try
                        kPrime_eval_for_interp_nested = max(cS.kGridV(1), min(cS.kGridV(end), kPrime_inner_nested));
                        expected_future_value_nested = EVprime_of_kprime_interp_1D(kPrime_eval_for_interp_nested);
                    catch ME_interp_k_only
                        warning('BellmanInner_nested_k_only (V8): 插值错误 (年龄 %d, kP %.2e): %s. 将使用边界值。', ...
                                a_idx_current, kPrime_inner_nested, ME_interp_k_only.message);
                        if kPrime_inner_nested < cS.kGridV(1)
                            expected_future_value_nested = EVprime_of_kprime_interp_1D(cS.kGridV(1));
                        else
                            expected_future_value_nested = EVprime_of_kprime_interp_1D(cS.kGridV(end));
                        end
                    end
                    if ~isfinite(expected_future_value_nested)
                        expected_future_value_nested = -1e12;
                    end
                    s_transition_to_next_year_nested = cS.s_1yr_transitionV(a_idx_current);
                    Val_out_nested = util_current_period_nested + cS.beta * s_transition_to_next_year_nested * expected_future_value_nested;
                end
                if ~isfinite(Val_out_nested)
                    negVal_out_nested = 1e12;
                    Val_out_nested = -1e12;
                else
                    negVal_out_nested = -Val_out_nested;
                end
            end

            if kPrime_max_bound <= kPrime_min_bound
                kPrime_out = kPrime_min_bound;
                c_expenditure_at_corner = budget_for_c_expend_and_kprime - kPrime_out;
                c_quantity = max(cS.cFloor, c_expenditure_at_corner / (1 + cS.tau_c) );
                [~, utility_at_corner] = main_olg_v8_utils.CES_utility(c_quantity, cS.sigma, cS);
                if a_idx_current < cS.aD_new
                    EV_at_corner = -Inf;
                    try
                        kPrime_eval_corner = max(cS.kGridV(1),min(cS.kGridV(end),kPrime_out));
                        EV_at_corner = EVprime_of_kprime_interp_1D(kPrime_eval_corner);
                    catch
                        if kPrime_out<cS.kGridV(1), EV_at_corner=EVprime_of_kprime_interp_1D(cS.kGridV(1));
                        else, EV_at_corner=EVprime_of_kprime_interp_1D(cS.kGridV(end));end;
                    end
                    if ~isfinite(EV_at_corner), EV_at_corner=-1e12; end
                    s_transition_osa_k = cS.s_1yr_transitionV(a_idx_current);
                    ValueFunc_out = utility_at_corner + cS.beta*s_transition_osa_k*EV_at_corner;
                else
                    ValueFunc_out = utility_at_corner;
                end
                if ~isfinite(ValueFunc_out), ValueFunc_out = -1e12; end
            else
                kPrime_min_for_fminbnd = max(cS.kMin, kPrime_min_bound);
                kPrime_max_for_fminbnd = max(kPrime_min_for_fminbnd + 1e-9, min(cS.kMax, kPrime_max_bound));
                kPrime_optimal_fminbnd = kPrime_min_for_fminbnd; 
                negValue_fminbnd = Inf; % Initialize
                
                if kPrime_min_for_fminbnd >= kPrime_max_for_fminbnd
                    [negValue_fminbnd, ~] = BellmanInner_nested_k_only(kPrime_optimal_fminbnd);
                else
                    objective_func_for_fminbnd = @(kP_choice) negBellmanObjective_nested_k_only(kP_choice);
                    [kPrime_optimal_fminbnd, negValue_fminbnd, exitflag_fminbnd] = ...
                        fminbnd(objective_func_for_fminbnd, kPrime_min_for_fminbnd, kPrime_max_for_fminbnd, fminbnd_opts_in);
                    if exitflag_fminbnd <= 0 || ...
                       abs(kPrime_optimal_fminbnd - kPrime_min_for_fminbnd) < 1e-7 || ...
                       abs(kPrime_optimal_fminbnd - kPrime_max_for_fminbnd) < 1e-7
                        [negValue_at_min_bound, ~] = BellmanInner_nested_k_only(kPrime_min_for_fminbnd);
                        [negValue_at_max_bound, ~] = BellmanInner_nested_k_only(kPrime_max_for_fminbnd);
                        if negValue_at_min_bound <= negValue_at_max_bound + 1e-9
                            kPrime_optimal_fminbnd = kPrime_min_for_fminbnd;
                            negValue_fminbnd = negValue_at_min_bound;
                        else
                            kPrime_optimal_fminbnd = kPrime_max_for_fminbnd;
                            negValue_fminbnd = negValue_at_max_bound;
                        end
                    end
                end
                ValueFunc_out = -negValue_fminbnd;
                kPrime_out = kPrime_optimal_fminbnd;
                c_expenditure_optimal = budget_for_c_expend_and_kprime - kPrime_out;
                c_quantity = max(cS.cFloor, c_expenditure_optimal / (1 + cS.tau_c) );
                if ~isfinite(ValueFunc_out), ValueFunc_out = -1e12; end
            end

            kPrime_out = max(cS.kMin, min(cS.kMax, kPrime_out));
            c_expenditure_final_check = budget_for_c_expend_and_kprime - kPrime_out;
            c_quantity = max(cS.cFloor, c_expenditure_final_check / (1 + cS.tau_c) );
            
            [~, utility_final_check] = main_olg_v8_utils.CES_utility(c_quantity, cS.sigma, cS);
            if a_idx_current < cS.aD_new
                EV_final_check = -Inf;
                try
                    kPrime_eval_final_check = max(cS.kGridV(1),min(cS.kGridV(end),kPrime_out));
                    EV_final_check = EVprime_of_kprime_interp_1D(kPrime_eval_final_check);
                catch
                    if kPrime_out<cS.kGridV(1), EV_final_check=EVprime_of_kprime_interp_1D(cS.kGridV(1));
                    else, EV_final_check=EVprime_of_kprime_interp_1D(cS.kGridV(end));end;
                end
                if ~isfinite(EV_final_check), EV_final_check = -1e12; end
                s_transition_final_check = cS.s_1yr_transitionV(a_idx_current);
                ValueFunc_out = utility_final_check + cS.beta * s_transition_final_check * EV_final_check;
            else
                ValueFunc_out = utility_final_check;
            end

            if ~isfinite(kPrime_out), kPrime_out = cS.kMin; end
            if ~isfinite(c_quantity), c_quantity = cS.cFloor; end
            if ~isfinite(ValueFunc_out), ValueFunc_out = -1e12; end
        end

        % --- V8 家庭生命周期路径模拟器 ---
        function [kHistM_out, kPpsHistM_out, cHistM_out] = HHSimulation_olgm(...
                kPolM_4D_input, cPpsPolM_choice_4D_input, cPolM_consump_q_4D_input, eIdxM_annual_input, ...
                R_k_net_factor_hh_sim, w_gross_sim_price, TR_total_sim_transfer, bV_payg_sim_benefit, ...
                paramS_sim_household, cS_common_sim)
            % HHSimulation_olgm - 基于最优策略模拟家庭生命周期路径
            % 
            % 输入：
            %   kPolM_4D_input - 最优非PPS储蓄策略函数
            %   cPpsPolM_choice_4D_input - 最优PPS缴费选择策略（V8关键输入）
            %   cPolM_consump_q_4D_input - 最优消费策略函数
            %   eIdxM_annual_input - 劳动效率冲击路径
            %   R_k_net_factor_hh_sim - 税后资本回报因子
            %   w_gross_sim_price - 市场毛工资率
            %   TR_total_sim_transfer - 总转移支付
            %   bV_payg_sim_benefit - PAYG福利向量
            %   paramS_sim_household - 家庭参数
            %   cS_common_sim - 模型参数
            % 输出：
            %   kHistM_out - 非PPS资产生命周期路径矩阵
            %   kPpsHistM_out - PPS资产生命周期路径矩阵
            %   cHistM_out - 消费生命周期路径矩阵
            % 
            % 功能：使用VFI求解得到的策略函数，模拟所有家庭的生命周期决策
            %       包括内生PPS缴费选择和资产积累过程

            ageToGroupMap_sim = zeros(cS_common_sim.aD_orig,1);
            for a_map_idx_sim = 1:cS_common_sim.aD_new
                idx_map_sim = cS_common_sim.physAgeMap{a_map_idx_sim};
                if ~isempty(idx_map_sim)
                    ageToGroupMap_sim(idx_map_sim) = a_map_idx_sim;
                end
            end

            nSim_sim = size(eIdxM_annual_input,1);
            kHistM_out_temp    = zeros(nSim_sim, cS_common_sim.aD_orig + 1);
            kPpsHistM_out_temp = zeros(nSim_sim, cS_common_sim.aD_orig + 1);
            cHistM_out    = zeros(nSim_sim, cS_common_sim.aD_orig);

            leGridV_col_sim = paramS_sim_household.leGridV(:);

            kPolInterp_sim         = cell(cS_common_sim.nw, cS_common_sim.aD_new);
            cPpsPolInterp_choice_sim = cell(cS_common_sim.nw, cS_common_sim.aD_new); % V8
            cPolqInterp_sim        = cell(cS_common_sim.nw, cS_common_sim.aD_new);

            for ia_interp = 1:cS_common_sim.aD_new
                for ie_interp = 1:cS_common_sim.nw
                    if cS_common_sim.nk > 1 && cS_common_sim.nkpps > 1
                        kPolInterp_sim{ie_interp,ia_interp}    = griddedInterpolant({cS_common_sim.kGridV, cS_common_sim.kppsGridV}, squeeze(kPolM_4D_input(:,:,ie_interp,ia_interp)), 'linear', 'nearest');
                        cPpsPolInterp_choice_sim{ie_interp,ia_interp} = griddedInterpolant({cS_common_sim.kGridV, cS_common_sim.kppsGridV}, squeeze(cPpsPolM_choice_4D_input(:,:,ie_interp,ia_interp)), 'linear', 'nearest'); % V8
                        cPolqInterp_sim{ie_interp,ia_interp}   = griddedInterpolant({cS_common_sim.kGridV, cS_common_sim.kppsGridV}, squeeze(cPolM_consump_q_4D_input(:,:,ie_interp,ia_interp)), 'linear', 'nearest');
                    elseif cS_common_sim.nk > 1 && cS_common_sim.nkpps == 1
                        kPolInterp_sim{ie_interp,ia_interp}    = griddedInterpolant(cS_common_sim.kGridV, squeeze(kPolM_4D_input(:,1,ie_interp,ia_interp)), 'linear', 'nearest');
                        cPpsPolInterp_choice_sim{ie_interp,ia_interp} = griddedInterpolant(cS_common_sim.kGridV, squeeze(cPpsPolM_choice_4D_input(:,1,ie_interp,ia_interp)), 'linear', 'nearest'); % V8
                        cPolqInterp_sim{ie_interp,ia_interp}   = griddedInterpolant(cS_common_sim.kGridV, squeeze(cPolM_consump_q_4D_input(:,1,ie_interp,ia_interp)), 'linear', 'nearest');
                    elseif cS_common_sim.nk == 1 && cS_common_sim.nkpps > 1
                        kPolInterp_sim{ie_interp,ia_interp}    = griddedInterpolant(cS_common_sim.kppsGridV, squeeze(kPolM_4D_input(1,:,ie_interp,ia_interp))', 'linear', 'nearest');
                        cPpsPolInterp_choice_sim{ie_interp,ia_interp} = griddedInterpolant(cS_common_sim.kppsGridV, squeeze(cPpsPolM_choice_4D_input(1,:,ie_interp,ia_interp))', 'linear', 'nearest'); % V8
                        cPolqInterp_sim{ie_interp,ia_interp}   = griddedInterpolant(cS_common_sim.kppsGridV, squeeze(cPolM_consump_q_4D_input(1,:,ie_interp,ia_interp))', 'linear', 'nearest');
                    elseif cS_common_sim.nk == 1 && cS_common_sim.nkpps == 1
                        kPolInterp_sim{ie_interp,ia_interp}    = @(x,y) kPolM_4D_input(1,1,ie_interp,ia_interp);
                        cPpsPolInterp_choice_sim{ie_interp,ia_interp} = @(x,y) cPpsPolM_choice_4D_input(1,1,ie_interp,ia_interp); % V8
                        cPolqInterp_sim{ie_interp,ia_interp}   = @(x,y) cPolM_consump_q_4D_input(1,1,ie_interp,ia_interp);
                    else
                        error('HHSimulation_olgm (V8): nk 或 nkpps 为零，无法创建插值器。');
                    end
                end
            end

            pps_return_net_annual_factor_sim = 1 + ((R_k_net_factor_hh_sim - 1) + cS_common_sim.pps_return_rate_premium);

            for a_orig_loop_idx = 1:cS_common_sim.aD_orig
                a_new_group_idx_sim = ageToGroupMap_sim(a_orig_loop_idx);
                kNowV_annual_sim    = kHistM_out_temp(:, a_orig_loop_idx);
                kPpsNowV_annual_sim = kPpsHistM_out_temp(:, a_orig_loop_idx);
                kNextNonPpsV_from_policy = zeros(nSim_sim,1);
                cPpsDecisionFromPol_choice = zeros(nSim_sim,1); % V8: From optimal choice policy
                cConsumpValV_q_from_policy = zeros(nSim_sim,1);
                kPpsNextV_annual_sim     = zeros(nSim_sim,1);

                is_working_age_annual_sim = (a_orig_loop_idx < cS_common_sim.aR_idx_orig);
                is_pps_withdrawal_eligible_annual_sim = (~is_working_age_annual_sim && ...
                    cS_common_sim.pps_active && ...
                    a_orig_loop_idx >= cS_common_sim.pps_withdrawal_age_min_idx);
                pps_withdrawal_pretax_this_year_sim = zeros(nSim_sim,1);
                if is_pps_withdrawal_eligible_annual_sim
                    pps_withdrawal_pretax_this_year_sim = kPpsNowV_annual_sim * cS_common_sim.pps_withdrawal_rate;
                end
                
                actual_cpps_final_for_period_sim = zeros(nSim_sim,1); % This will store the cPps from policy

                for ie_sim_idx = 1 : cS_common_sim.nw
                    simIdxV_for_this_e = find(eIdxM_annual_input(:, a_orig_loop_idx) == ie_sim_idx);
                    if isempty(simIdxV_for_this_e), continue; end

                    kNow_clamped    = max(cS_common_sim.kGridV(1), min(cS_common_sim.kGridV(end), kNowV_annual_sim(simIdxV_for_this_e)));
                    kPpsNow_clamped = max(cS_common_sim.kppsGridV(1), min(cS_common_sim.kppsGridV(end), kPpsNowV_annual_sim(simIdxV_for_this_e)));

                    if cS_common_sim.nk > 1 && cS_common_sim.nkpps > 1
                        kNextNonPpsV_from_policy(simIdxV_for_this_e) = kPolInterp_sim{ie_sim_idx, a_new_group_idx_sim}(kNow_clamped, kPpsNow_clamped);
                        cPpsDecisionFromPol_choice(simIdxV_for_this_e) = cPpsPolInterp_choice_sim{ie_sim_idx, a_new_group_idx_sim}(kNow_clamped, kPpsNow_clamped); % V8
                        cConsumpValV_q_from_policy(simIdxV_for_this_e) = cPolqInterp_sim{ie_sim_idx, a_new_group_idx_sim}(kNow_clamped, kPpsNow_clamped);
                    elseif cS_common_sim.nk > 1
                        kNextNonPpsV_from_policy(simIdxV_for_this_e) = kPolInterp_sim{ie_sim_idx, a_new_group_idx_sim}(kNow_clamped);
                        cPpsDecisionFromPol_choice(simIdxV_for_this_e) = cPpsPolInterp_choice_sim{ie_sim_idx, a_new_group_idx_sim}(kNow_clamped); % V8
                        cConsumpValV_q_from_policy(simIdxV_for_this_e) = cPolqInterp_sim{ie_sim_idx, a_new_group_idx_sim}(kNow_clamped);
                    elseif cS_common_sim.nkpps > 1
                        kNextNonPpsV_from_policy(simIdxV_for_this_e) = kPolInterp_sim{ie_sim_idx, a_new_group_idx_sim}(kPpsNow_clamped);
                        cPpsDecisionFromPol_choice(simIdxV_for_this_e) = cPpsPolInterp_choice_sim{ie_sim_idx, a_new_group_idx_sim}(kPpsNow_clamped); % V8
                        cConsumpValV_q_from_policy(simIdxV_for_this_e) = cPolqInterp_sim{ie_sim_idx, a_new_group_idx_sim}(kPpsNow_clamped);
                    else
                        kNextNonPpsV_from_policy(simIdxV_for_this_e) = kPolInterp_sim{ie_sim_idx, a_new_group_idx_sim}(kNow_clamped, kPpsNow_clamped);
                        cPpsDecisionFromPol_choice(simIdxV_for_this_e) = cPpsPolInterp_choice_sim{ie_sim_idx, a_new_group_idx_sim}(kNow_clamped, kPpsNow_clamped); % V8
                        cConsumpValV_q_from_policy(simIdxV_for_this_e) = cPolqInterp_sim{ie_sim_idx, a_new_group_idx_sim}(kNow_clamped, kPpsNow_clamped);
                    end
                    
                    % In V8, cPpsDecisionFromPol_choice already reflects the optimal, constrained choice.
                    % No further rule application or constraint needed here for simulation, as it's part of the policy.
                    actual_cpps_final_for_period_sim(simIdxV_for_this_e) = cPpsDecisionFromPol_choice(simIdxV_for_this_e);
                end

                if cS_common_sim.pps_active
                    kPpsNextV_annual_sim = (kPpsNowV_annual_sim + actual_cpps_final_for_period_sim - pps_withdrawal_pretax_this_year_sim) * pps_return_net_annual_factor_sim;
                    kPpsNextV_annual_sim = max(cS_common_sim.kppsMin, min(cS_common_sim.kppsMax, kPpsNextV_annual_sim));
                else
                    kPpsNextV_annual_sim = kPpsNowV_annual_sim;
                end

                kHistM_out_temp(:, a_orig_loop_idx + 1) = max(cS_common_sim.kMin, min(cS_common_sim.kMax, kNextNonPpsV_from_policy));
                kPpsHistM_out_temp(:, a_orig_loop_idx + 1) = kPpsNextV_annual_sim;
                cHistM_out(:, a_orig_loop_idx) = max(cS_common_sim.cFloor, cConsumpValV_q_from_policy);
            end

            kHistM_out = kHistM_out_temp(:, 1:cS_common_sim.aD_orig);
            kPpsHistM_out = kPpsHistM_out_temp(:, 1:cS_common_sim.aD_orig);
        end


        % --- V8 一般均衡求解器：核心宏观经济算法 ---
        function [K_sol_out, tau_l_sol_out, gbc_res_final_out, converged_and_feasible_out, solution_details_out] = solve_K_tau_l_for_rho_prime(...
                rho_prime_payg_target_input, K_init_guess_input, cS_global, paramS_global_in, eIdxM_global_sim_paths)
            % solve_K_tau_l_for_rho_prime - V8模型的一般均衡求解器
            % 
            % 输入：
            %   rho_prime_payg_target_input - 目标PAYG替代率ρ'
            %   K_init_guess_input - 总资本的初始猜测值
            %   cS_global - 模型参数
            %   paramS_global_in - 派生参数
            %   eIdxM_global_sim_paths - 劳动效率路径
            % 输出：
            %   K_sol_out - 均衡总资本
            %   tau_l_sol_out - 均衡所得税率
            %   gbc_res_final_out - 政府预算残差
            %   converged_and_feasible_out - 收敛标志
            %   solution_details_out - 均衡详情
            % 
            % 功能：迭代求解满足以下条件的一般均衡：
            %   1. 资本市场出清：K_model = K_guess（家庭决策生成的资本=价格计算用资本）
            %   2. 政府预算平衡：一般预算收支平衡（通过内生调整tau_l）
            %   3. PAYG系统平衡：实现目标替代率（通过内生调整theta_payg）
            %   4. 家庭最优化：给定价格下家庭选择最优(c,k',c_pps)
            % 
            % 算法流程：
            %   - 外层迭代：调整K_guess和tau_l直至市场出清和预算平衡
            %   - 内层VFI：在给定价格下求解家庭的内生PPS缴费决策问题
            %   - 收敛检查：验证所有市场出清条件和约束条件
            K_current_guess = K_init_guess_input;
            tau_l_current_guess = cS_global.tau_l_init_guess;
            L_per_capita_global = paramS_global_in.L_per_capita;
            mass_workers_global = paramS_global_in.mass_workers_group;

            maxIter_ktl_loop = cS_global.max_iter_K_tau_l;
            tol_norm_ktl_loop = cS_global.tol_K_tau_l;
            dampK_ktl_loop = cS_global.damp_K_v5;
            damp_tau_l_ktl_loop = cS_global.damp_tau_l_v5;

            converged_and_feasible_out = false;
            K_sol_out = NaN; tau_l_sol_out = NaN; gbc_res_final_out = Inf;
            solution_details_out = struct();

            mass_retirees_global = sum(paramS_global_in.ageMassV(cS_global.aR_new + 1 : cS_global.aD_new));
            theta_payg_required_calc = 0;
            if mass_workers_global > 1e-9
                theta_payg_required_calc = rho_prime_payg_target_input * (mass_retirees_global / mass_workers_global);
            else
                if rho_prime_payg_target_input > 1e-9, theta_payg_required_calc = Inf;
                else, theta_payg_required_calc = 0;
                end
            end
            theta_payg_required_calc = max(0, theta_payg_required_calc);
            solution_details_out.theta_payg_required_before_cap = theta_payg_required_calc;

            if theta_payg_required_calc > cS_global.theta_payg_max + 1e-5
                if ~isfield(paramS_global_in, 'suppress_initial_theta_print') || ~paramS_global_in.suppress_initial_theta_print
                    fprintf('  solve_K_tau_l (V8): rho_prime_target=%.4f 导致理论theta_req=%.4f > theta_max=%.3f. 直接标记为不可行。\n', ...
                        rho_prime_payg_target_input, theta_payg_required_calc, cS_global.theta_payg_max);
                end
                converged_and_feasible_out = false;
                K_sol_out = K_init_guess_input; tau_l_sol_out = tau_l_current_guess; gbc_res_final_out = Inf;
                solution_details_out.theta_payg = min(theta_payg_required_calc, cS_global.theta_payg_max);
                solution_details_out.MPL_gross = NaN; solution_details_out.R_mkt_gross_factor = NaN; solution_details_out.b_payg = NaN;
                solution_details_out.T_bequest_Model = NaN; solution_details_out.C_model = NaN; solution_details_out.Y_model = NaN;
                solution_details_out.K_model_pps = NaN; solution_details_out.K_model_non_pps = NaN;
                return;
            end

            stagnation_counter_ktl = 0;
            prev_devNorm_ktl = Inf;
            tau_l_boundary_strike_count_ktl = 0;

            if ~isfield(paramS_global_in, 'suppress_inner_print_header') || ~paramS_global_in.suppress_inner_print_header
                fprintf('  solve_K_tau_l_for_rho_prime (V8): rho_prime_target=%.4f (理论theta_req=%.4f), K_init=%.2f, tau_l_init=%.3f\n', ...
                        rho_prime_payg_target_input, theta_payg_required_calc, K_current_guess, tau_l_current_guess);
                fprintf('  IterKTL | K_guess  | tau_l_gs | MPL_g    | theta_act| K_tot_mod| K_pps_mod| GBC_res  | K_dev    | tau_l_dev| Norm     | Improv   | Strikes\n');
                fprintf('  -------------------------------------------------------------------------------------------------------------------------------------------\n');
            end
            
            MPL_gross_iter_val = NaN; R_mkt_gross_factor_iter_val = NaN; theta_payg_actual_iter_val = NaN; b_payg_iter_val = NaN;
            T_bequest_model_iter_val = NaN; C_model_iter_val = NaN; Y_for_gbc_iter_val = NaN; gbc_residual_iter_val = Inf;
            K_model_from_sim_iter_val = K_current_guess; K_dev_from_sim_iter_val = Inf;
            K_model_pps_sim_iter_val = NaN; K_model_nonpps_sim_iter_val = NaN;

            for iter_ktl_idx = 1:maxIter_ktl_loop
                [R_mkt_gross_factor_iter_val, MPL_gross_iter_val] = main_olg_v8_utils.HHPrices_Huggett(K_current_guess, L_per_capita_global, cS_global);
                r_mkt_gross_iter_val = R_mkt_gross_factor_iter_val - 1;

                avg_worker_gross_wage_iter_val = 0;
                if mass_workers_global > 1e-9 && L_per_capita_global > 0 && MPL_gross_iter_val > 0
                    avg_worker_gross_wage_iter_val = (MPL_gross_iter_val * L_per_capita_global) / mass_workers_global;
                end
                b_payg_iter_val = rho_prime_payg_target_input * avg_worker_gross_wage_iter_val;
                b_payg_iter_val = max(0, b_payg_iter_val);

                theta_payg_actual_iter_val = theta_payg_required_calc;
                if (theta_payg_actual_iter_val + tau_l_current_guess) > cS_global.max_total_labor_tax
                    theta_payg_actual_iter_val = max(0, cS_global.max_total_labor_tax - tau_l_current_guess);
                end
                theta_payg_actual_iter_val = min(theta_payg_actual_iter_val, cS_global.theta_payg_max);
                theta_payg_actual_iter_val = max(0, theta_payg_actual_iter_val);

                r_k_net_hh_iter_val = r_mkt_gross_iter_val * (1 - cS_global.tau_k);
                R_k_net_hh_factor_iter_val = 1 + r_k_net_hh_iter_val;

                bV_payg_vec_iter_val = zeros(1, cS_global.aD_new);
                if cS_global.aR_new < cS_global.aD_new
                    bV_payg_vec_iter_val(cS_global.aR_new+1 : cS_global.aD_new) = b_payg_iter_val;
                end

                paramS_for_vfi_sim_iter = paramS_global_in;
                paramS_for_vfi_sim_iter.tau_l = tau_l_current_guess;
                paramS_for_vfi_sim_iter.theta_payg_actual_for_hh = theta_payg_actual_iter_val;
                paramS_for_vfi_sim_iter.pps_tax_deferral_active = cS_global.pps_active;

                TR_total_for_vfi_guess_val = 0.01 * MPL_gross_iter_val;
                if iter_ktl_idx > 1 && isfinite(T_bequest_model_iter_val)
                    TR_total_for_vfi_guess_val = T_bequest_model_iter_val;
                end
                max_vfi_tr_sub_iter = 5;
                tol_vfi_tr_sub_iter = 1e-3;
                cPolM_4D_from_vfi_final = []; kPolM_4D_from_vfi_final = []; 
                cPpsPolM_4D_choice_from_vfi_final = []; % V8
                TR_total_for_vfi_final_iter = TR_total_for_vfi_guess_val;

                for i_vfi_tr_sub_loop = 1:max_vfi_tr_sub_iter
                    [cPolM_vfi_temp, kPolM_vfi_temp, cPpsPolM_vfi_temp_choice, ~] = ... % V8 output
                        main_olg_v8_utils.HHSolution_VFI_Huggett(R_k_net_hh_factor_iter_val, MPL_gross_iter_val, ...
                        TR_total_for_vfi_guess_val, bV_payg_vec_iter_val, paramS_for_vfi_sim_iter, cS_global);
                    [kHistM_sim_for_bequest, ~, ~] = ...
                        main_olg_v8_utils.HHSimulation_olgm(kPolM_vfi_temp, cPpsPolM_vfi_temp_choice, cPolM_vfi_temp, ... % V8 input
                        eIdxM_global_sim_paths, R_k_net_hh_factor_iter_val, MPL_gross_iter_val, ...
                        TR_total_for_vfi_guess_val, bV_payg_vec_iter_val, paramS_for_vfi_sim_iter, cS_global);
                    
                    ageDeathMass_annual_iter_val = paramS_global_in.Z_ss_norm_annual(:) .* cS_global.d_orig(:);
                    mean_bequest_wealth_per_age_iter_val = mean(kHistM_sim_for_bequest, 1);
                    TotalBequests_pc_iter_val = sum(mean_bequest_wealth_per_age_iter_val(:) .* ageDeathMass_annual_iter_val(:));
                    T_bequest_model_new_iter_val = TotalBequests_pc_iter_val / (1 + paramS_global_in.popGrowthForDebt);
                    T_bequest_model_new_iter_val = max(0, T_bequest_model_new_iter_val);

                    T_bequest_model_iter_val = T_bequest_model_new_iter_val;
                    cPolM_4D_from_vfi_final = cPolM_vfi_temp;
                    kPolM_4D_from_vfi_final = kPolM_vfi_temp;
                    cPpsPolM_4D_choice_from_vfi_final = cPpsPolM_vfi_temp_choice; % V8
                    TR_total_for_vfi_final_iter = T_bequest_model_new_iter_val;

                    if abs(T_bequest_model_new_iter_val - TR_total_for_vfi_guess_val) < tol_vfi_tr_sub_iter || ...
                       i_vfi_tr_sub_loop == max_vfi_tr_sub_iter
                        break;
                    end
                    TR_total_for_vfi_guess_val = 0.5 * TR_total_for_vfi_guess_val + 0.5 * T_bequest_model_new_iter_val;
                end
                T_bequest_model_iter_val = TR_total_for_vfi_final_iter;

                [kHistM_non_pps_sim_iter_val, kPpsHistM_sim_iter_val, cHistM_sim_iter_val] = main_olg_v8_utils.HHSimulation_olgm(...
                    kPolM_4D_from_vfi_final, cPpsPolM_4D_choice_from_vfi_final, cPolM_4D_from_vfi_final, eIdxM_global_sim_paths, ... % V8
                    R_k_net_hh_factor_iter_val, MPL_gross_iter_val, TR_total_for_vfi_final_iter, bV_payg_vec_iter_val, paramS_for_vfi_sim_iter, cS_global);

                K_model_nonpps_sim_iter_val = mean(kHistM_non_pps_sim_iter_val, 1) * paramS_global_in.Z_ss_norm_annual;
                K_model_pps_sim_iter_val = 0;
                if cS_global.pps_active && cS_global.pps_in_K && (cS_global.pps_max_contrib_frac > 0 || cS_global.pps_annual_contrib_limit > 0)
                    K_model_pps_sim_iter_val = mean(kPpsHistM_sim_iter_val, 1) * paramS_global_in.Z_ss_norm_annual;
                    K_model_pps_sim_iter_val = max(0, K_model_pps_sim_iter_val);
                end
                K_model_from_sim_iter_val = K_model_nonpps_sim_iter_val + K_model_pps_sim_iter_val;
                K_model_from_sim_iter_val = max(1e-6, K_model_from_sim_iter_val);
                C_model_iter_val = mean(cHistM_sim_iter_val,1) * paramS_global_in.Z_ss_norm_annual;

                Y_for_gbc_iter_val = cS_global.A * (K_current_guess^cS_global.alpha) * (L_per_capita_global^(1-cS_global.alpha));
                G_val_target_iter_val = cS_global.gov_exp_frac_Y * Y_for_gbc_iter_val;
                B_val_target_iter_val = cS_global.gov_debt_frac_Y * Y_for_gbc_iter_val;

                gbc_residual_iter_val = main_olg_v8_utils.check_gbc_residual(K_current_guess, C_model_iter_val, Y_for_gbc_iter_val, ...
                    G_val_target_iter_val, B_val_target_iter_val, MPL_gross_iter_val, r_mkt_gross_iter_val, ...
                    theta_payg_actual_iter_val, tau_l_current_guess, ...
                    b_payg_iter_val, T_bequest_model_iter_val, 0, cS_global, paramS_global_in);

                K_dev_from_sim_iter_val = K_current_guess - K_model_from_sim_iter_val;
                tau_l_dev_raw_for_update = -gbc_residual_iter_val / (MPL_gross_iter_val * L_per_capita_global + 1e-9);
                current_devNorm_val = sqrt(K_dev_from_sim_iter_val^2 + (gbc_residual_iter_val)^2 );
                norm_improvement_val = prev_devNorm_ktl - current_devNorm_val;

                fprintf('  %7d | %8.4f | %8.4f | %8.4f | %8.4f | %8.4f | %8.4f | %8.2e | %8.4f | %8.4f | %8.2e | %.1e | %7d\n', ...
                    iter_ktl_idx, K_current_guess, tau_l_current_guess, MPL_gross_iter_val, theta_payg_actual_iter_val, ...
                    K_model_from_sim_iter_val, K_model_pps_sim_iter_val, gbc_residual_iter_val, ...
                    K_dev_from_sim_iter_val, tau_l_dev_raw_for_update, current_devNorm_val, norm_improvement_val, tau_l_boundary_strike_count_ktl);

                payg_fully_funded_by_actual_theta_check = (theta_payg_actual_iter_val >= theta_payg_required_calc - 1e-5);

                if current_devNorm_val < tol_norm_ktl_loop && ...
                   abs(gbc_residual_iter_val) < cS_global.gbc_tol_for_internal_loop && ...
                   payg_fully_funded_by_actual_theta_check
                    converged_and_feasible_out = true;
                    K_sol_out = K_model_from_sim_iter_val;
                    tau_l_sol_out = tau_l_current_guess;
                    gbc_res_final_out = gbc_residual_iter_val;
                    solution_details_out.R_mkt_gross_factor = R_mkt_gross_factor_iter_val;
                    solution_details_out.MPL_gross = MPL_gross_iter_val;
                    solution_details_out.theta_payg = theta_payg_actual_iter_val;
                    solution_details_out.b_payg = b_payg_iter_val;
                    solution_details_out.T_bequest_Model = T_bequest_model_iter_val;
                    solution_details_out.C_model = C_model_iter_val;
                    solution_details_out.Y_model = cS_global.A * (K_sol_out^cS_global.alpha) * (L_per_capita_global^(1-cS_global.alpha));
                    solution_details_out.K_model_pps = K_model_pps_sim_iter_val;
                    solution_details_out.K_model_non_pps = K_model_nonpps_sim_iter_val;
                    fprintf('  solve_K_tau_l (V8): K和tau_l成功收敛 (rho_prime_target=%.4f, 实际theta_act=%.4f).\n', ...
                        rho_prime_payg_target_input, theta_payg_actual_iter_val);
                    break;

                elseif current_devNorm_val < tol_norm_ktl_loop && ...
                       abs(gbc_residual_iter_val) < cS_global.gbc_tol_for_internal_loop && ...
                       ~payg_fully_funded_by_actual_theta_check
                    fprintf('  solve_K_tau_l (V8): K, tau_l, GBC收敛 (rho_prime=%.4f), 但实际theta_payg (%.4f) 因总税负上限低于理论需求 (%.4f)。标记为不可行。\n', ...
                        rho_prime_payg_target_input, theta_payg_actual_iter_val, theta_payg_required_calc);
                    converged_and_feasible_out = false;
                    K_sol_out = K_model_from_sim_iter_val; tau_l_sol_out = tau_l_current_guess; gbc_res_final_out = gbc_residual_iter_val;
                    solution_details_out.R_mkt_gross_factor = R_mkt_gross_factor_iter_val;
                    solution_details_out.MPL_gross = MPL_gross_iter_val;
                    solution_details_out.theta_payg = theta_payg_actual_iter_val;
                    solution_details_out.b_payg = b_payg_iter_val;
                    solution_details_out.T_bequest_Model = T_bequest_model_iter_val;
                    solution_details_out.C_model = C_model_iter_val;
                    solution_details_out.Y_model = Y_for_gbc_iter_val;
                    solution_details_out.K_model_pps = K_model_pps_sim_iter_val;
                    solution_details_out.K_model_non_pps = K_model_nonpps_sim_iter_val;
                    break;
                end

                K_guess_next_val = K_current_guess - dampK_ktl_loop * K_dev_from_sim_iter_val;
                K_current_guess = max(1e-3, K_guess_next_val);
                new_tau_l_unconstrained_val = tau_l_current_guess + damp_tau_l_ktl_loop * tau_l_dev_raw_for_update;
                tau_l_next_iter_constrained_val = max(cS_global.tau_l_min, min(cS_global.tau_l_max, new_tau_l_unconstrained_val));
                is_tau_l_at_boundary_now = ...
                    (abs(tau_l_next_iter_constrained_val - cS_global.tau_l_max) < 1e-7 && new_tau_l_unconstrained_val >= cS_global.tau_l_max) || ...
                    (abs(tau_l_next_iter_constrained_val - cS_global.tau_l_min) < 1e-7 && new_tau_l_unconstrained_val <= cS_global.tau_l_min);
                if is_tau_l_at_boundary_now && abs(gbc_residual_iter_val) > cS_global.gbc_tol_for_internal_loop
                    tau_l_boundary_strike_count_ktl = tau_l_boundary_strike_count_ktl + 1;
                else
                    tau_l_boundary_strike_count_ktl = 0;
                end
                tau_l_current_guess = tau_l_next_iter_constrained_val;

                if tau_l_boundary_strike_count_ktl >= cS_global.max_tau_l_boundary_strikes && abs(gbc_residual_iter_val) > cS_global.gbc_tol_for_internal_loop
                    fprintf('  警告 (V8): tau_l 在边界 (%.4f) 持续 %d 次迭代，且GBC (%.2e) 未平衡。为 rho_prime=%.4f 中止。\n', ...
                        tau_l_current_guess, tau_l_boundary_strike_count_ktl, gbc_residual_iter_val, rho_prime_payg_target_input);
                    converged_and_feasible_out = false;
                    K_sol_out = K_model_from_sim_iter_val; tau_l_sol_out = tau_l_current_guess; gbc_res_final_out = gbc_residual_iter_val;
                    solution_details_out.R_mkt_gross_factor = R_mkt_gross_factor_iter_val; solution_details_out.MPL_gross = MPL_gross_iter_val;
                    solution_details_out.theta_payg = theta_payg_actual_iter_val; solution_details_out.b_payg = b_payg_iter_val;
                    solution_details_out.T_bequest_Model = T_bequest_model_iter_val; solution_details_out.C_model = C_model_iter_val;
                    solution_details_out.Y_model = Y_for_gbc_iter_val;
                    solution_details_out.K_model_pps = K_model_pps_sim_iter_val; solution_details_out.K_model_non_pps = K_model_nonpps_sim_iter_val;
                    break;
                end

                if iter_ktl_idx > 1
                    if norm_improvement_val < (cS_global.min_norm_improvement_frac * prev_devNorm_ktl) && current_devNorm_val > tol_norm_ktl_loop
                        stagnation_counter_ktl = stagnation_counter_ktl + 1;
                    else
                        stagnation_counter_ktl = 0;
                    end
                end
                prev_devNorm_ktl = current_devNorm_val;

                if stagnation_counter_ktl >= cS_global.max_stagnation_iters && current_devNorm_val > tol_norm_ktl_loop
                    fprintf('  警告 (V8): 在 %d 次迭代后检测到范数停滞。范数: %.2e > 容忍度: %.1e。为 rho_prime=%.4f 中止。\n', ...
                        iter_ktl_idx, current_devNorm_val, tol_norm_ktl_loop, rho_prime_payg_target_input);
                    converged_and_feasible_out = false;
                    K_sol_out = K_model_from_sim_iter_val; tau_l_sol_out = tau_l_current_guess; gbc_res_final_out = gbc_residual_iter_val;
                    solution_details_out.R_mkt_gross_factor = R_mkt_gross_factor_iter_val; solution_details_out.MPL_gross = MPL_gross_iter_val;
                    solution_details_out.theta_payg = theta_payg_actual_iter_val; solution_details_out.b_payg = b_payg_iter_val;
                    solution_details_out.T_bequest_Model = T_bequest_model_iter_val; solution_details_out.C_model = C_model_iter_val;
                    solution_details_out.Y_model = Y_for_gbc_iter_val;
                    solution_details_out.K_model_pps = K_model_pps_sim_iter_val; solution_details_out.K_model_non_pps = K_model_nonpps_sim_iter_val;
                    break;
                end
            end

            if ~converged_and_feasible_out && iter_ktl_idx == maxIter_ktl_loop
                fprintf('  警告 (V8): K和tau_l迭代达到最大次数 (%d) 或在该次数内未达可行解 (rho_prime_target=%.4f).\n', maxIter_ktl_loop, rho_prime_payg_target_input);
                K_sol_out = K_model_from_sim_iter_val;
                tau_l_sol_out = tau_l_current_guess;
                gbc_res_final_out = gbc_residual_iter_val;
                if exist('MPL_gross_iter_val', 'var') && isfinite(MPL_gross_iter_val)
                    solution_details_out.R_mkt_gross_factor = R_mkt_gross_factor_iter_val;
                    solution_details_out.MPL_gross = MPL_gross_iter_val;
                    solution_details_out.theta_payg = theta_payg_actual_iter_val;
                    solution_details_out.b_payg = b_payg_iter_val;
                    solution_details_out.T_bequest_Model = T_bequest_model_iter_val;
                    solution_details_out.C_model = C_model_iter_val;
                    solution_details_out.Y_model = Y_for_gbc_iter_val;
                    solution_details_out.K_model_pps = K_model_pps_sim_iter_val;
                    solution_details_out.K_model_non_pps = K_model_nonpps_sim_iter_val;
                else
                    if ~isfield(solution_details_out, 'R_mkt_gross_factor'), solution_details_out.R_mkt_gross_factor = NaN; end
                    if ~isfield(solution_details_out, 'MPL_gross'), solution_details_out.MPL_gross = NaN; end
                    if ~isfield(solution_details_out, 'theta_payg'), solution_details_out.theta_payg = min(theta_payg_required_calc, cS_global.theta_payg_max); end
                    if ~isfield(solution_details_out, 'b_payg'), solution_details_out.b_payg = NaN; end
                    if ~isfield(solution_details_out, 'T_bequest_Model'), solution_details_out.T_bequest_Model = NaN; end
                    if ~isfield(solution_details_out, 'C_model'), solution_details_out.C_model = NaN; end
                    if ~isfield(solution_details_out, 'Y_model'), solution_details_out.Y_model = NaN; end
                    if ~isfield(solution_details_out, 'K_model_pps'), solution_details_out.K_model_pps = NaN; end
                    if ~isfield(solution_details_out, 'K_model_non_pps'), solution_details_out.K_model_non_pps = NaN; end
                end
            end

            if ~isfield(solution_details_out, 'theta_payg_required_before_cap') || ...
               (isfield(solution_details_out, 'theta_payg_required_before_cap') && isnan(solution_details_out.theta_payg_required_before_cap))
                 recalc_mass_retirees_final_sd = sum(paramS_global_in.ageMassV(cS_global.aR_new + 1 : cS_global.aD_new));
                 recalc_theta_req_final_sd = 0;
                 if mass_workers_global > 1e-9
                     recalc_theta_req_final_sd = rho_prime_payg_target_input * (recalc_mass_retirees_final_sd / mass_workers_global);
                 end
                 solution_details_out.theta_payg_required_before_cap = max(0, recalc_theta_req_final_sd);
            end
             if ~isfield(solution_details_out, 'K_model_pps')
                if exist('K_model_pps_sim_iter_val','var') && isfinite(K_model_pps_sim_iter_val)
                    solution_details_out.K_model_pps = K_model_pps_sim_iter_val;
                else
                    solution_details_out.K_model_pps = NaN;
                end
            end
             if ~isfield(solution_details_out, 'K_model_non_pps')
                 if exist('K_model_nonpps_sim_iter_val','var') && isfinite(K_model_nonpps_sim_iter_val)
                     solution_details_out.K_model_non_pps = K_model_nonpps_sim_iter_val;
                 else
                     solution_details_out.K_model_non_pps = NaN;
                 end
            end
            if converged_and_feasible_out && isfield(solution_details_out, 'MPL_gross') && isfinite(solution_details_out.MPL_gross)
                 solution_details_out.Y_model = cS_global.A * (K_sol_out^cS_global.alpha) * (L_per_capita_global^(1-cS_global.alpha));
            end
        end

        % --- check_gbc_residual (与V7/Baseline一致) ---
        function gbc_residual_out = check_gbc_residual(...
                K_val_market_input, C_val_model_input, Y_val_market_input, G_val_target_input, B_val_target_input, ...
                MPL_gross_val_input, r_mkt_gross_val_input, ...
                theta_payg_val_actual_input, tau_l_val_input, ...
                b_payg_val_per_retiree_input, T_bequest_val_pc_input, TR_gov_val_pc_input, ...
                cS_check, paramS_loc_check)
            % (Identical to main_olg_baseline_utils.m)
            L_per_capita_local_check = paramS_loc_check.L_per_capita;
            LaborTaxRev_general_part_calc = tau_l_val_input * MPL_gross_val_input * L_per_capita_local_check;
            CapitalTaxRev_calc = r_mkt_gross_val_input * K_val_market_input * cS_check.tau_k;
            ConsumptionTaxRev_calc = C_val_model_input * cS_check.tau_c;
            GeneralRevenue_calc = LaborTaxRev_general_part_calc + CapitalTaxRev_calc + ConsumptionTaxRev_calc;

            GovConsumption_calc = G_val_target_input;
            r_b_for_debt_service_calc = r_mkt_gross_val_input;
            DebtService_calc = (r_b_for_debt_service_calc - paramS_loc_check.popGrowthForDebt) * B_val_target_input;
            GovDirectTransfers_calc = TR_gov_val_pc_input;
            GeneralOutlays_calc = GovConsumption_calc + DebtService_calc + GovDirectTransfers_calc;
            gbc_residual_out = GeneralRevenue_calc - GeneralOutlays_calc;
        end

        % --- CallInterpolator (与V7/Baseline一致) ---
        function ev_val_out = CallInterpolator(interpolant_obj_input, k_val_input, k_pps_val_input, cS_local_interp)
            % (Identical to main_olg_baseline_utils.m)
            ev_val_out = -Inf;
            try
                if isa(interpolant_obj_input, 'griddedInterpolant')
                    if cS_local_interp.nk > 1 && cS_local_interp.nkpps > 1
                        ev_val_out = interpolant_obj_input(k_val_input, k_pps_val_input);
                    elseif cS_local_interp.nk > 1
                        ev_val_out = interpolant_obj_input(k_val_input);
                    elseif cS_local_interp.nkpps > 1
                        ev_val_out = interpolant_obj_input(k_pps_val_input);
                    else
                        if isscalar(interpolant_obj_input.Values)
                            ev_val_out = interpolant_obj_input.Values;
                        else
                            ev_val_out = interpolant_obj_input(k_val_input, k_pps_val_input);
                        end
                    end
                elseif isa(interpolant_obj_input, 'function_handle')
                    ev_val_out = interpolant_obj_input(k_val_input, k_pps_val_input);
                else
                    warning('CallInterpolator (V8): 未处理的插值器类型。');
                    ev_val_out = -1e11;
                end
            catch ME_call_interp_error
                warning('CallInterpolator (V8): 插值过程中发生错误: "%s"。将尝试限制输入值并重试。', ME_call_interp_error.message);
                k_clamped = k_val_input; 
                if cS_local_interp.nk > 0 && ~isempty(cS_local_interp.kGridV)
                    k_clamped = max(cS_local_interp.kGridV(1), min(cS_local_interp.kGridV(end), k_val_input));
                end
                k_pps_clamped = k_pps_val_input; 
                if cS_local_interp.nkpps > 0 && ~isempty(cS_local_interp.kppsGridV)
                    k_pps_clamped = max(cS_local_interp.kppsGridV(1), min(cS_local_interp.kppsGridV(end), k_pps_val_input));
                end
                try
                    if isa(interpolant_obj_input, 'griddedInterpolant')
                        if cS_local_interp.nk > 1 && cS_local_interp.nkpps > 1
                            ev_val_out = interpolant_obj_input(k_clamped, k_pps_clamped);
                        elseif cS_local_interp.nk > 1
                            ev_val_out = interpolant_obj_input(k_clamped);
                        elseif cS_local_interp.nkpps > 1
                            ev_val_out = interpolant_obj_input(k_pps_clamped);
                        else
                             if isscalar(interpolant_obj_input.Values)
                                ev_val_out = interpolant_obj_input.Values;
                            else
                                ev_val_out = interpolant_obj_input(k_clamped, k_pps_clamped);
                            end
                        end
                    elseif isa(interpolant_obj_input, 'function_handle')
                        ev_val_out = interpolant_obj_input(k_clamped, k_pps_clamped);
                    end
                catch
                    ev_val_out = -1e11;
                end
            end
            if ~isfinite(ev_val_out)
                ev_val_out = -1e12;
            end
        end

        % --- 性能比较测试函数 ---
        function compare_vfi_performance(cS)
            % compare_vfi_performance - 比较原版和优化版VFI的性能
            % 用于测试不同并行化策略的效果
            
            fprintf('\n=== VFI性能比较测试 ===\n');
            fprintf('网格参数：nk=%d, nkpps=%d, nw=%d, n_pps_choice_grid_points=%d\n', ...
                cS.nk, cS.nkpps, cS.nw, cS.n_pps_choice_grid_points);
            
            % 创建测试数据
            test_age = 5; % 中间年龄组
            vPrime_test = rand(cS.nk, cS.nkpps, cS.nw) * 100;
            R_k_test = 1.03;
            w_gross_test = 1.0;
            TR_total_test = 0.1;
            b_age_test = 0.05;
            
            % 创建测试参数
            paramS_test = struct();
            paramS_test.leTrProbM = eye(cS.nw) * 0.8 + (1-eye(cS.nw)) * 0.2/4; % 简单转移矩阵
            paramS_test.tau_l = 0.10; % 劳动所得税率
            paramS_test.theta_payg_actual_for_hh = 0.08; % PAYG税率
            epsilon_test = linspace(0.5, 1.5, cS.nw);
            
            fprintf('\n1. 测试原版VFI性能...\n');
            tic;
            [cPol_orig, kPol_orig, cPpsPol_orig, val_orig] = ...
                main_olg_v8_utils.HHSolutionByAge_VFI_Huggett_v8(...
                test_age, vPrime_test, R_k_test, w_gross_test, TR_total_test, b_age_test, ...
                paramS_test, cS, epsilon_test);
            time_orig = toc;
            fprintf('   原版完成时间: %.2f秒\n', time_orig);
            
            fprintf('\n2. 测试优化版VFI性能...\n');
            tic;
            [cPol_opt, kPol_opt, cPpsPol_opt, val_opt] = ...
                main_olg_v8_utils.HHSolutionByAge_VFI_Huggett_v8_optimized(...
                test_age, vPrime_test, R_k_test, w_gross_test, TR_total_test, b_age_test, ...
                paramS_test, cS, epsilon_test);
            time_opt = toc;
            fprintf('   优化版完成时间: %.2f秒\n', time_opt);
            
            % 计算性能改进
            speedup = time_orig / time_opt;
            fprintf('\n3. 性能比较结果：\n');
            fprintf('   加速比: %.2fx\n', speedup);
            if speedup > 1
                fprintf('   优化版比原版快 %.1f%%\n', (speedup-1)*100);
            else
                fprintf('   优化版比原版慢 %.1f%%\n', (1-speedup)*100);
            end
            
            % 验证结果一致性
            fprintf('\n4. 结果一致性检验：\n');
            max_c_diff = max(abs(cPol_orig(:) - cPol_opt(:)));
            max_k_diff = max(abs(kPol_orig(:) - kPol_opt(:)));
            max_cpps_diff = max(abs(cPpsPol_orig(:) - cPpsPol_opt(:)));
            max_val_diff = max(abs(val_orig(:) - val_opt(:)));
            
            fprintf('   消费策略最大差异: %.2e\n', max_c_diff);
            fprintf('   储蓄策略最大差异: %.2e\n', max_k_diff);
            fprintf('   PPS策略最大差异: %.2e\n', max_cpps_diff);
            fprintf('   值函数最大差异: %.2e\n', max_val_diff);
            
            tolerance = 1e-10;
            if max_c_diff < tolerance && max_k_diff < tolerance && ...
               max_cpps_diff < tolerance && max_val_diff < tolerance
                fprintf('   ✓ 结果一致性检验通过\n');
            else
                fprintf('   ✗ 结果存在显著差异，需要检查优化实现\n');
            end
            
            fprintf('\n=== 性能测试完成 ===\n');
        end
        
        % --- 进一步优化建议函数 ---
        function suggest_optimization_strategies(cS)
            % suggest_optimization_strategies - 分析并建议进一步的优化策略
            
            fprintf('\n=== VFI进一步优化建议 ===\n');
            
            % 计算当前计算复杂度
            states_per_age = cS.nk * cS.nkpps * cS.nw;
            total_pps_evaluations = states_per_age * cS.n_pps_choice_grid_points;
            total_ages = cS.aD_new;
            total_computations = total_pps_evaluations * total_ages;
            
            fprintf('当前计算规模分析：\n');
            fprintf('  每个年龄组状态数: %d\n', states_per_age);
            fprintf('  每个状态PPS选择数: %d\n', cS.n_pps_choice_grid_points);
            fprintf('  每个年龄组总优化次数: %d\n', total_pps_evaluations);
            fprintf('  总年龄组数: %d\n', total_ages);
            fprintf('  总计算量: %d次优化\n', total_computations);
            
            fprintf('\n优化策略建议：\n');
            
            % 建议1：调整网格密度
            fprintf('1. 网格参数调优：\n');
            if cS.n_pps_choice_grid_points > 5
                fprintf('   - 考虑减少PPS缴费网格点数(当前%d)到5-7个\n', cS.n_pps_choice_grid_points);
                savings1 = (1 - 5/cS.n_pps_choice_grid_points) * 100;
                fprintf('     预期加速: %.1f%%\n', savings1);
            end
            
            if cS.nk > 50
                fprintf('   - 考虑适度减少k网格密度(当前%d)\n', cS.nk);
            end
            
            if cS.nkpps > 30
                fprintf('   - 考虑适度减少k_pps网格密度(当前%d)\n', cS.nkpps);
            end
            
            % 建议2：算法改进
            fprintf('\n2. 算法改进策略：\n');
            fprintf('   - 实现自适应PPS网格：根据收入水平动态调整网格密度\n');
            fprintf('   - 使用插值加速：预计算关键点，插值得到中间值\n');
            fprintf('   - 实现多级并行：CPU核心数允许时使用嵌套parfor\n');
            
            % 建议3：数值优化
            fprintf('\n3. 数值计算优化：\n');
            fprintf('   - 预计算转移概率矩阵，避免重复计算\n');
            fprintf('   - 缓存中间结果，减少重复的收入计算\n');
            fprintf('   - 使用GPU加速矩阵运算（如果可用）\n');
            
            % 建议4：内存优化
            fprintf('\n4. 内存访问优化：\n');
            fprintf('   - 重组数据结构，提高缓存命中率\n');
            fprintf('   - 减少大矩阵的broadcast开销\n');
            fprintf('   - 使用单精度浮点数（如果精度允许）\n');
            
            % 计算理论最大加速比
            max_cores = feature('numcores');
            current_parallel_states = cS.nk;
            optimal_parallel_states = min(states_per_age, max_cores);
            theoretical_speedup = optimal_parallel_states / current_parallel_states;
            
            fprintf('\n5. 并行化潜力分析：\n');
            fprintf('   - 当前CPU核心数: %d\n', max_cores);
            fprintf('   - 当前并行状态数: %d\n', current_parallel_states);
            fprintf('   - 理论最优并行数: %d\n', optimal_parallel_states);
            fprintf('   - 理论最大加速比: %.2fx\n', theoretical_speedup);
            
            fprintf('\n=== 优化建议完成 ===\n');
        end

    end % 结束 Static 方法块
end % 结束 main_olg_v8_utils 类定义
% --- END OF FILE main_olg_v8_utils.m ---
