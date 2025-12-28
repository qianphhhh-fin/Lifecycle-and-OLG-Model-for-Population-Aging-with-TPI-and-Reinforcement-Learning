% --- 开始文件：main_olg_v8_utils.m ---

% =====================================================================
% === OLG 模型 V8 工具函数库: 内生PPS缴费决策的数值实现 (按年龄组模拟) ===
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
% 8. 家庭模拟器 - 按年龄组模拟生命周期路径（取代年度模拟）
% 
% V8模型核心创新：
% - HHSolutionByAge_VFI_Huggett_v8: 实现内生PPS缴费选择的VFI
% - 对每个状态(k,k_pps,ε,a)，个体选择最优的(c,k',c_pps)组合
% - PPS缴费受双重约束：收入比例上限和年度绝对上限
% - 通过离散网格搜索实现PPS缴费的优化选择
% - HHSimulation_olgm: 按年龄组模拟，简化计算过程并提高效率
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
            fprintf('V8: 开始设置参数...\n');
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
            
            % 🔧 计算年龄组死亡率（从年度死亡率聚合而来）
            cS.d_new = zeros(cS.aD_new, 1);
            for a = 1:cS.aD_new
                annual_indices_in_group = cS.physAgeMap{a};
                if ~isempty(annual_indices_in_group)
                    % 使用年龄组内年度死亡率的平均值
                    cS.d_new(a) = mean(cS.d_orig(annual_indices_in_group));
                else
                    cS.d_new(a) = 0.05; % 默认死亡率
                end
            end

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
            cS.beta       = 0.97;          % 主观贴现因子β_disc  
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
            cS.nw = 3;                      % 效率状态网格点数

            % --- 资产网格参数 ---
            % 非PPS资产网格（对应v8.tex中的k状态空间）
            cS.tgKY = 3;                    % 目标资本产出比
            cS.tgWage = (1-cS.alpha)*cS.A*((cS.tgKY/cS.A)^(cS.alpha/(1-cS.alpha)));
            cS.nk = 5;                     % 非PPS资产网格点数
            cS.kMin = 0;                    % 非PPS资产下界
            cS.kMax = 15 * cS.tgWage;      % 非PPS资产上界
            power = 1.5;                    % 网格密度参数
            kGridV = cS.kMin + (cS.kMax - cS.kMin) * (linspace(0, 1, cS.nk).^power); %
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
            cS.use_continuous_optimization = true; % <<<<<< 默认使用fmincon

            % --- PPS制度基础参数 ---
            cS.pps_active = true;                        % PPS制度激活标志
            cS.pps_tax_rate_withdrawal = 0.3;          % PPS提取阶段税率
            cS.pps_return_rate_premium = 0.08;          % <<<<<<< 修改这里，例如改为0.03 >>>>>>> PPS超额收益率 (相对于市场净回报r_k_net_hh)
            cS.pps_withdrawal_rate = 0.15;              % 退休后年度提取比例
            cS.pps_in_K = true;                         % PPS资产是否计入生产性资本
            cS.pps_bequeathable = true;                 % <<<<<< 修改这里为 true >>>>>>> PPS资产是否可遗赠

            % --- V8模型关键创新：PPS缴费约束参数 ---
            % 实现v8.tex中描述的双重约束机制
            cS.pps_contrib_limit = 9999;          % PPS年度绝对缴费上限
            cS.pps_max_contrib_frac = 1;             % PPS缴费占劳动收入比例上限
            cS.pps_contribution_age_max_idx = cS.aR_idx_orig - 1;  % 最大缴费年度年龄
            cS.pps_withdrawal_age_min_idx = cS.aR_idx_orig;        % 最低提取年度年龄
            
            % V8核心：PPS缴费选择网格点数
            % 用于在VFI中离散化PPS缴费选择空间
            cS.n_pps_choice_grid_points = 12;           % 进一步减少格点数，改善平滑度
            cS.power_pps_choice_grid = 1.3;             % 降低幂次参数，减少非线性程度

            % --- PPS资产网格（对应v8.tex中的k_pps状态空间） ---
            cS.nkpps = 5;                               % PPS资产网格点数
            cS.kppsMin = 0;                             % PPS资产下界
            cS.kppsMax = 2; % cS.kMax / 2;                   % PPS资产上界
            if cS.nkpps > 0
                cS.kppsMax = max(cS.kppsMax, 1e-3);
            end
            power_kpps = 1.5;                           % PPS资产网格密度参数
            if cS.nkpps > 1
                kppsGridV_temp = cS.kppsMin + (cS.kppsMax - cS.kppsMin) * (linspace(0, 1, cS.nkpps).^power_kpps); % .^power_kpps
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
            cS.damp_K_v5 = 0.1;                        % K更新阻尼系数
            cS.damp_tau_l_v5 = 0.1;                    % tau_l更新阻尼系数
            cS.gbc_tol_for_internal_loop = 1e-3;       % 政府预算平衡容忍度

            % --- 收敛检测参数 ---
            cS.max_stagnation_iters = 10;              % 最大停滞迭代次数
            cS.min_norm_improvement_frac = 1e-3;       % 最小改进比例
            cS.max_tau_l_boundary_strikes = 5;         % tau_l边界冲击最大次数
            
            % --- PAYG税率约束参数 ---
            cS.tau_l_init_guess = 0.1509;                % 所得税率初始猜测
            cS.tau_l_min = 0.00;                       % 所得税率下界
            cS.tau_l_max = 0.3;                        % 所得税率上界
            cS.max_total_labor_tax = 1;              % 总劳动税负上限
            cS.theta_payg_max = 1;                  % PAYG税率上限

            % --- 参考性PPS缴费时间表（V8中不再直接使用） ---
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

                    % --- 神经网络VFI近似参数 ---
        cS.use_NN_for_VFI = false; % <<<< 新增: 是否使用NN替代VFI内层循环
        cS.trained_nn_filename = 'trained_VFI_NN_v8.mat';
        cS.normalization_params_filename = 'normalization_params_VFI_NN_v8.mat';
        % 这些文件名应该与 train_NN_VFI.m 中保存的文件名一致

            fprintf('V8: 完整参数已设置完毕。\n');
            fprintf('    主要创新：内生PPS缴费选择，缴费网格点数=%d\n', cS.n_pps_choice_grid_points);
            fprintf('    PPS约束：收入比例上限=%.1f%%, 年度绝对上限=%.2f\n', ...
                cS.pps_max_contrib_frac*100, cS.pps_contrib_limit);
            fprintf('    PPS可遗赠 (cS.pps_bequeathable): %s\n', mat2str(cS.pps_bequeathable)); % <-- 新增的打印
        end

        function cS = generateGrids(cS)
            % generateGrids - 重新生成网格参数
            % 
            % 输入：
            %   cS - 参数结构体
            % 输出：
            %   cS - 更新了网格的参数结构体
            % 
            % 功能：根据当前的网格参数设置重新生成kGridV和kppsGridV
            
            % 重新生成非PPS资产网格
            power = 1.5; % 网格密度参数
            kGridV = cS.kMin + (cS.kMax - cS.kMin) * (linspace(0, 1, cS.nk).^power);
            if cS.nk > 0
                kGridV(1) = cS.kMin;
            end
            cS.kGridV = kGridV(:);
            
            % 重新生成PPS资产网格
            power_kpps = 1.5; % PPS资产网格密度参数
            if cS.nkpps > 1
                kppsGridV_temp = cS.kppsMin + (cS.kppsMax - cS.kppsMin) * (linspace(0, 1, cS.nkpps).^power_kpps);
                kppsGridV_temp(1) = cS.kppsMin;
            elseif cS.nkpps == 1
                kppsGridV_temp = cS.kppsMin;
            else
                kppsGridV_temp = [];
            end
            cS.kppsGridV = kppsGridV_temp(:);
            
            % 输出网格生成信息
            fprintf('网格已重新生成：nk=%d, nkpps=%d\n', cS.nk, cS.nkpps);
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
        
        % --- 劳动禀赋过程生成器 ---
        function [leLogGridV, leTrProbM, leProb1V] = EarningProcess_olgm(cS)
            % EarningProcess_olgm - 生成劳动效率冲击的离散Markov链
            %
            % 使用标准AR(1)过程: log(ε') = ρ*log(ε) + σ*u, u~N(0,1)
            % 
            % 输入：
            %   cS - 包含nw的参数结构体
            % 输出：
            %   leLogGridV - 对数劳动效率网格 (nw x 1)
            %   leTrProbM - 转移概率矩阵 (nw x nw)
            %   leProb1V - 初始分布概率 (nw x 1)

            fprintf('劳动禀赋过程参数已生成 (Tauchen & norm_grid)。\n');

            % 保守的AR(1)参数设置（确保合理的效率分布）
            lePersistence = 0.90;    % 持续性参数ρ（降低以减少极端值）
            leShockStd = 0.15;       % 创新标准差σ（降低以减少变异）
            Tauchen_q = 2.0;         % Tauchen网格范围参数（缩小范围）

            % 使用标准Tauchen方法离散化AR(1)过程
            % log(ε') = ρ*log(ε) + σ*u, u~N(0,1)
            mew = 0;                 % AR(1)过程的常数项（设为0）
            rho = lePersistence;     % 自相关系数
            sigma = leShockStd;      % 创新标准差
            znum = cS.nw;           % 状态数量

            % 调用标准Tauchen方法
            % [leLogGridV_raw, leTrProbM] = discretizeAR1_Tauchen(mew, rho, sigma, znum, Tauchen_q);
            [leLogGridV_raw, leTrProbM] = main_olg_v8_utils.tauchen(znum, rho, sigma, mew, Tauchen_q);


            % 标准化：使对数效率网格以0为中心
            % 这样确保效率分布相对对称
            leLogGridV = leLogGridV_raw - mean(leLogGridV_raw);

            % 进一步检查并调整极端值
            leGridV_test = exp(leLogGridV);
            efficiency_ratio = leGridV_test(end) / leGridV_test(1);

            % 如果比值仍然过大，进行额外压缩
            max_acceptable_ratio = 5.0; % 设定最大可接受比值为5
            if efficiency_ratio > max_acceptable_ratio
                compression_factor = log(max_acceptable_ratio) / log(efficiency_ratio);
                leLogGridV = leLogGridV * compression_factor;
                fprintf('效率分布压缩因子: %.3f\n', compression_factor);
            end

            % 计算平稳分布作为初始分布
            % 通过求解 π = π * P 得到
            try
                % 方法1：特征向量法
                [eigenvecs, eigenvals] = eig(leTrProbM');
                [~, idx] = min(abs(diag(eigenvals) - 1));
                leProb1V = eigenvecs(:, idx);
                leProb1V = leProb1V / sum(leProb1V);
                leProb1V = abs(leProb1V); % 确保为正数
            catch
                % 方法2：数值迭代法（备选）
                leProb1V = ones(cS.nw, 1) / cS.nw; % 初始均匀分布
                for iter = 1:1000
                    leProb1V_new = leTrProbM' * leProb1V;
                    if norm(leProb1V_new - leProb1V) < 1e-10
                        break;
                    end
                    leProb1V = leProb1V_new;
                end
            end

            % 确保概率和为1
            if sum(leProb1V) > 1e-9
                 leProb1V = leProb1V / sum(leProb1V);
            else
                 leProb1V = ones(cS.nw, 1) / cS.nw; % Fallback to uniform if sum is too small
                 warning('EarningProcess_olgm: 平稳分布概率和过小，已重置为均匀分布。');
            end


            % 最终验证和报告
            leGridV_final = exp(leLogGridV);
            efficiency_ratio_final = leGridV_final(end) / leGridV_final(1);
            mean_efficiency = sum(leGridV_final .* leProb1V);

            % 输出诊断信息
            fprintf('AR(1)参数: ρ=%.3f, σ=%.3f, q=%.1f\n', rho, sigma, Tauchen_q);
            fprintf('效率网格范围: [%.4f, %.4f], 比值=%.2f\n', ...
                leGridV_final(1), leGridV_final(end), efficiency_ratio_final);
            fprintf('平均效率: %.4f\n', mean_efficiency);
        end

        function [y_grid_out, trProbM_out] = tauchen(N_states, persistence_rho, shock_sigma, mean_val_mu, num_std_dev_width)
            % (Identical to main_olg_baseline_utils.m)
            std_y_unconditional = sqrt(shock_sigma^2 / (1 - persistence_rho^2));
            if abs(1-persistence_rho) < 1e-9
                std_y_unconditional = shock_sigma * 100; % A large number if persistence is 1
                % warning('Tauchen: persistence_rho (%.4f) 接近1，std_y_unconditional可能不准确。', persistence_rho);
            end

            y_max_boundary = num_std_dev_width * std_y_unconditional;
            y_min_boundary = -y_max_boundary;
            y_grid_centered = linspace(y_min_boundary, y_max_boundary, N_states);
             if N_states == 1 % Handle single state case
                y_grid_centered = 0; % Or some other appropriate single value
            end


            step_size_d = 0;
            if N_states > 1
                step_size_d = y_grid_centered(2) - y_grid_centered(1);
            end

            trProbM_calc = zeros(N_states, N_states);
            if N_states == 1
                trProbM_calc(1,1) = 1.0; % If only one state, it transitions to itself
            else
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
            end


            row_sums_check = sum(trProbM_calc,2);
            row_sums_check(row_sums_check <= 1e-9) = 1; % Avoid division by zero if a row sum is tiny
            trProbM_out = bsxfun(@rdivide, trProbM_calc, row_sums_check);

            unconditional_mean_shift = mean_val_mu / (1-persistence_rho);
             if abs(1-persistence_rho) < 1e-9 && mean_val_mu ~= 0
                unconditional_mean_shift = sign(mean_val_mu) * inf;
                % warning('Tauchen: rho=1 and mu non-zero, unconditional mean is ill-defined.');
            elseif abs(1-persistence_rho) < 1e-9 && mean_val_mu == 0
                unconditional_mean_shift = 0;
            end
            
            y_grid_out = y_grid_centered;
            if isfinite(unconditional_mean_shift)
                y_grid_out = y_grid_centered + unconditional_mean_shift;
            end
             y_grid_out = y_grid_out(:); % Ensure column vector
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
            % MarkovChainSimulation - 标准马尔可夫链模拟（保留原版本用于兼容性）
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
            if num_simulations > 0 && num_periods_sim > 0 && num_states > 0
                eIdxM_out(:, 1) = 1 + sum(bsxfun(@gt, random_numbers_rvInM(:,1), cumulative_initial_prob_cP0), 2);
            elseif num_states == 0 && num_simulations > 0 && num_periods_sim > 0
                 warning('MarkovChainSimulation: num_states is 0. Cannot simulate.');
                 return; % Return empty or error, as appropriate
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

        function eIdxM_group_out = MarkovChainSimulation_AgeGroup(num_simulations, cS, initial_prob_dist_p0V, transition_matrix_trProbM)
            % MarkovChainSimulation_AgeGroup - 直接模拟年龄组的马尔可夫链
            % 
            % 输入：
            %   num_simulations - 模拟个体数量
            %   cS - 包含年龄组信息的参数结构体
            %   initial_prob_dist_p0V - 初始状态分布
            %   transition_matrix_trProbM - 转移概率矩阵
            % 输出：
            %   eIdxM_group_out - 年龄组效率冲击索引矩阵 (nSim × aD_new)
            
            % 设置随机数种子确保可重复性
            rng(433);
            
            num_states = length(initial_prob_dist_p0V);
            num_age_groups = cS.aD_new;
            
            % 参数验证
            if size(transition_matrix_trProbM,1) ~= num_states || size(transition_matrix_trProbM,2) ~= num_states
                error('MarkovChainSimulation_AgeGroup: 转移矩阵维度与初始分布长度不匹配。');
            end
            
            % 归一化概率分布
            if abs(sum(initial_prob_dist_p0V) - 1) > 1e-5
                initial_prob_dist_p0V = initial_prob_dist_p0V ./ sum(initial_prob_dist_p0V);
            end
            if any(abs(sum(transition_matrix_trProbM, 2) - 1) > 1e-5)
                row_sums_tr = sum(transition_matrix_trProbM, 2);
                row_sums_tr(row_sums_tr <= 1e-9) = 1;
                transition_matrix_trProbM = bsxfun(@rdivide, transition_matrix_trProbM, row_sums_tr);
            end
            
            % 生成年龄组随机数
            random_numbers_group = rand(num_simulations, num_age_groups);
            
            % 计算累积概率
            cumulative_initial_prob_cP0 = cumsum(initial_prob_dist_p0V(:)');
            cumulative_transition_prob_cPT = cumsum(transition_matrix_trProbM, 2);
            if num_states > 0
                cumulative_transition_prob_cPT(:, num_states) = 1.0;
            end
            
            % 初始化结果矩阵
            eIdxM_group_out = zeros(num_simulations, num_age_groups, 'uint16');
            
            % 第一个年龄组：使用初始分布
            if num_simulations > 0 && num_age_groups > 0 && num_states > 0
                eIdxM_group_out(:, 1) = 1 + sum(bsxfun(@gt, random_numbers_group(:,1), cumulative_initial_prob_cP0), 2);
            end
            
            % 后续年龄组：使用转移概率
            for a_group = 2:num_age_groups
                current_state_indices = eIdxM_group_out(:, a_group - 1);
                
                % 验证状态索引有效性
                valid_indices = (current_state_indices >= 1) & (current_state_indices <= num_states);
                if ~all(valid_indices)
                    warning('MarkovChainSimulation_AgeGroup: 在年龄组 %d 检测到无效状态索引。已重置为状态1。', a_group-1);
                    current_state_indices(~valid_indices) = 1;
                    eIdxM_group_out(:, a_group - 1) = current_state_indices;
                end
                
                % 计算下一期状态
                cPt_for_next_state = cumulative_transition_prob_cPT(current_state_indices, :);
                eIdxM_group_out(:, a_group) = 1 + sum(bsxfun(@gt, random_numbers_group(:, a_group), cPt_for_next_state), 2);
            end
            
            fprintf('年龄组马尔可夫链模拟完成 (%d 个体, %d 年龄组)。\n', num_simulations, num_age_groups);
        end

        function eIdxM_group = LaborEndowSimulation_olgm(cS, paramS)
            % LaborEndowSimulation_olgm - 生成年龄组劳动禀赋路径
            % 🔧 已修改为年龄组版本，与Python版本保持一致
            rng(433);
            eIdxM_group = main_olg_v8_utils.MarkovChainSimulation_AgeGroup(cS.nSim, cS, ...
                                                          paramS.leProb1V, paramS.leTrProbM);
            fprintf('劳动禀赋路径已模拟 (%d 个体, %d 年龄组)。\n', cS.nSim, cS.aD_new);
        end

        function eIdxM_group = LaborEndowSimulation_olgm_AgeGroup(cS, paramS)
            % LaborEndowSimulation_olgm_AgeGroup - 直接生成年龄组的劳动禀赋路径
            % 
            % 输入：
            %   cS - 包含年龄组参数的结构体
            %   paramS - 包含马尔可夫链参数的结构体
            % 输出：
            %   eIdxM_group - 年龄组效率冲击索引矩阵 (nSim × aD_new)
            %
            % 功能：直接为年龄组生成效率冲击序列，避免年度到年龄组的转换
            
            eIdxM_group = main_olg_v8_utils.MarkovChainSimulation_AgeGroup(cS.nSim, cS, ...
                                                                          paramS.leProb1V, paramS.leTrProbM);
            fprintf('年龄组劳动禀赋路径已生成 (%d 个体, %d 年龄组)。\n', cS.nSim, cS.aD_new);
        end

        function [HHlaborM_group, L_total_eff_pc] = LaborSupply_Huggett(eIdxM_group, cS, paramS, Z_ss_norm_group)
            % LaborSupply_Huggett - 计算劳动供给（年龄组版本）
            % 🔧 完全使用年龄组逻辑，与Python版本保持一致
            % 
            % 输入参数：
            %   eIdxM_group - 年龄组效率冲击矩阵 (nSim × aD_new)
            %   cS - 模型参数
            %   paramS - 参数结构体
            %   Z_ss_norm_group - 年龄组人口分布
            %
            % 输出：
            %   HHlaborM_group - 年龄组劳动供给矩阵 (nSim × aD_new)
            %   L_total_eff_pc - 总体人均有效劳动供给

            nSim = size(eIdxM_group, 1);
            leGridV_col_local = paramS.leGridV(:);
            
            % 🔧 只处理年龄组数据，不再支持年度数据转换
            if size(eIdxM_group, 2) ~= cS.aD_new
                error('LaborSupply_Huggett: 效率冲击矩阵必须是年龄组格式 (%d × %d)', nSim, cS.aD_new);
            end
            
            fprintf('LaborSupply_Huggett: 处理年龄组效率冲击数据 (%d × %d)。\n', nSim, cS.aD_new);
            
            % 初始化输出矩阵
            HHlaborM_group = zeros(nSim, cS.aD_new);
            
            % 🔧 添加边界检查，确保索引不超出leGridV范围
            eIdxM_group_clipped = max(1, min(eIdxM_group, length(leGridV_col_local)));
            
            % 直接在年龄组层面计算劳动供给
            for a_group = 1:cS.aD_new
                if a_group <= cS.aR_new  % 工作年龄组
                    current_eIdx = eIdxM_group_clipped(:, a_group);
                    labor_eff_for_group = leGridV_col_local(current_eIdx);
                    HHlaborM_group(:, a_group) = cS.ageEffV_new(a_group) .* labor_eff_for_group;
                end
            end
            
            % 计算总体人均有效劳动供给
            L_total_eff_pc = 0;
            if cS.aR_new > 0
                mean_labor_per_working_group = mean(HHlaborM_group(:, 1:cS.aR_new), 1);
                L_total_eff_pc = mean_labor_per_working_group * Z_ss_norm_group(1:cS.aR_new);
            end
            L_total_eff_pc = max(0, L_total_eff_pc);
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
                % warning('HHPrices: K_productive 非正，已重置为 %.1e。', K_productive);
            end
            if L_total_eff <= 0
                L_total_eff=1e-6;
                % warning('HHPrices: L_total_eff 非正，已重置为 %.1e。', L_total_eff);
            end
            Y_gross = cS.A * (K_productive^cS.alpha) * (L_total_eff^(1-cS.alpha));
            MPK_gross_val = cS.alpha * Y_gross / K_productive;
            MPL_gross = (1-cS.alpha) * Y_gross / L_total_eff;
            R_market_gross_factor = 1 + MPK_gross_val - cS.ddk;
            R_market_gross_factor = max(1.0 + 1e-6, R_market_gross_factor); % Ensure R > 1
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
            muM   =  Inf(size(cM_quantity)); % High marginal utility if below floor
            
            if abs(sigma_crra - 1) < 1e-6 % Log utility (sigma = 1)
                utilM(is_valid_consumption) = log(c_adjusted_quantity(is_valid_consumption));
                muM(is_valid_consumption)   = 1 ./ c_adjusted_quantity(is_valid_consumption);
            else % CRRA utility (sigma ~= 1)
                utilM(is_valid_consumption) = (c_adjusted_quantity(is_valid_consumption).^(1-sigma_crra)) ./ (1-sigma_crra);
                muM(is_valid_consumption)   = c_adjusted_quantity(is_valid_consumption).^(-sigma_crra);
            end

            % Penalize consumption below floor
            % The penalty is applied to utility, and marginal utility is high due to c_adjusted_quantity being cFloor
            utilM(~is_valid_consumption) = -1e10 - (min_c_quantity - cM_quantity(~is_valid_consumption))*1e10;
            % For muM, if c_adjusted_quantity is used, it already gives high mu for c_adjusted_quantity = cFloor
            % If cM_quantity is used for muM below floor, it would be even higher or Inf.
            % The current muM calculation based on c_adjusted_quantity is standard.
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
                
                % 处理PPS提取税收
                pps_withdrawal_gross = 0;
                pps_withdrawal_tax = 0;
                pps_withdrawal_net = 0;
                
                % 检查是否有PPS提取
                if isfield(paramS_hh, 'current_pps_withdrawal') && paramS_hh.current_pps_withdrawal > 0
                    pps_withdrawal_gross = paramS_hh.current_pps_withdrawal;
                    
                    % 计算PPS提取税
                    pps_tax_rate_withdrawal = 0.15; % 默认提取税率
                    if isfield(paramS_hh, 'pps_tax_rate_withdrawal')
                        pps_tax_rate_withdrawal = paramS_hh.pps_tax_rate_withdrawal;
                    end
                    
                    pps_withdrawal_tax = pps_withdrawal_gross * pps_tax_rate_withdrawal;
                    pps_withdrawal_net = pps_withdrawal_gross - pps_withdrawal_tax;
                end
                
                non_capital_income = TR_total + b_payg_val + pps_withdrawal_net;
            end

            % 计算资本收入税后净值
            capital_income_gross = (R_k_net_factor - 1) * k_now_val / (1 - cS.tau_k);
            capital_income_tax = cS.tau_k * max(0, capital_income_gross);
            capital_income_net = capital_income_gross - capital_income_tax;
            
            resources_for_c_and_k_prime = k_now_val + capital_income_net + non_capital_income - actual_pps_contribution_expenditure;

            if ~isfinite(resources_for_c_and_k_prime)
                resources_for_c_and_k_prime = -1e10; % Or some other very small number to indicate infeasibility
                % warning('HHIncome_Huggett (V8): 计算得到的资源为非有限值。');
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
            %   cS - 模型参数 (应包含 use_continuous_optimization 标志)
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

            fprintf('  🔄 VFI V8 (HHSolution_VFI_Huggett): 开始逆向迭代...\n');
            total_age_groups_vfi = cS.aD_new;
            vfi_start_time = tic;
            
            % 初始化进度条变量
            progress_bar_length = 40;  % 进度条长度
            last_progress_shown = -1;  % 上次显示的进度百分比
            
            for a_idx = cS.aD_new : -1 : 1
                % 计算当前进度
                current_step = cS.aD_new - a_idx + 1;
                progress_pct = current_step / total_age_groups_vfi * 100;
                
                % 更新进度条（每5%或最后一个年龄组更新一次）
                if progress_pct - last_progress_shown >= 5 || a_idx == 1 || a_idx == cS.aD_new
                    % 计算进度条填充
                    filled_length = round(progress_pct / 100 * progress_bar_length);
                    bar_str = [repmat('█', 1, filled_length), repmat('░', 1, progress_bar_length - filled_length)];
                    
                    % 计算预估剩余时间
                    elapsed_time = toc(vfi_start_time);
                    if current_step > 1
                        estimated_total_time = elapsed_time * total_age_groups_vfi / current_step;
                        eta_seconds = estimated_total_time - elapsed_time;
                        eta_str = sprintf('ETA: %dm%.0fs', floor(eta_seconds/60), mod(eta_seconds, 60));
                    else
                        eta_str = 'ETA: --';
                    end
                    
                    % 显示进度条
                    fprintf('\r    📊 VFI进度: [%s] %3.0f%% (%2d/%2d) | %s', ...
                        bar_str, progress_pct, current_step, total_age_groups_vfi, eta_str);
                    
                    last_progress_shown = progress_pct;
                end

                vPrime_kkppse_next = [];
                if a_idx < cS.aD_new
                    vPrime_kkppse_next = valM(:,:,:,a_idx+1);
                end
                eps_grid_for_vfi = paramS_vfi.leGridV;

                % 调用按年龄组求解的函数 (V8 version)
                % 根据 cS.use_continuous_optimization 选择优化方法
                if ~isfield(cS, 'use_continuous_optimization')
                    warning('HHSolution_VFI_Huggett: cS.use_continuous_optimization 未设置，默认为 true (fmincon)。');
                    effective_use_continuous_optimization = true;
                else
                    effective_use_continuous_optimization = cS.use_continuous_optimization;
                end

                if effective_use_continuous_optimization
                    [cPolM_q(:,:,:,a_idx), kPolM(:,:,:,a_idx), cPpsPolM_choice(:,:,:,a_idx), valM(:,:,:,a_idx)] = ...
                        main_olg_v8_utils.HHSolutionByAge_VFI_Huggett_v8_fmincon(a_idx, vPrime_kkppse_next, ...
                        R_k_net_factor_vfi, w_gross_vfi, TR_total_vfi, bV_payg_vfi(a_idx), paramS_vfi, cS, eps_grid_for_vfi);
                else
                    % 使用离散选择的 PPS 优化 (nested_parallel 是其中一个实现)
                    [cPolM_q(:,:,:,a_idx), kPolM(:,:,:,a_idx), cPpsPolM_choice(:,:,:,a_idx), valM(:,:,:,a_idx)] = ...
                        main_olg_v8_utils.HHSolutionByAge_VFI_Huggett_v8_nested_parallel(a_idx, vPrime_kkppse_next, ...
                        R_k_net_factor_vfi, w_gross_vfi, TR_total_vfi, bV_payg_vfi(a_idx), paramS_vfi, cS, eps_grid_for_vfi);
                end

            end
            
            % 显示完成信息
            total_vfi_time = toc(vfi_start_time);
            fprintf('\n  ✅ VFI V8 (HHSolution_VFI_Huggett): 完成! 总耗时: %.2f秒\n', total_vfi_time);
        end

      

        % --- V8 家庭生命周期路径模拟器 (按年龄组模拟) ---
        function [kHistM_out, kPpsHistM_out, cHistM_out, cppsHistM_out] = HHSimulation_olgm(...
                kPolM_4D_input, cPpsPolM_choice_4D_input, cPolM_consump_q_4D_input, eIdxM_group, ...
                R_k_net_factor_hh_sim, w_gross_sim_price, TR_total_sim_transfer, bV_payg_sim_benefit, ...
                paramS_sim_household, cS_common_sim)
            % HHSimulation_olgm - 基于最优策略模拟家庭生命周期路径（年龄组版本）
            % 🔧 完全使用年龄组逻辑，与Python版本保持一致
            % 
            % 输入参数：
            %   eIdxM_group - 年龄组效率冲击矩阵 (nSim × aD_new)
            % 
            % 输出：
            %   kHistM_out - 非PPS资产路径 (nSim × aD_new)
            %   kPpsHistM_out - PPS资产路径 (nSim × aD_new)
            %   cHistM_out - 消费路径 (nSim × aD_new)
            %   cppsHistM_out - PPS缴费路径 (nSim × aD_new)
            
            nSim_sim = size(eIdxM_group, 1);
            
            % 🔧 只处理年龄组数据，不再支持年度数据转换
            if size(eIdxM_group, 2) ~= cS_common_sim.aD_new
                error('HHSimulation_olgm: 效率冲击矩阵必须是年龄组格式 (%d × %d)', nSim_sim, cS_common_sim.aD_new);
            end

            % 初始化结果矩阵（按年龄组）
            kHistM_out = zeros(nSim_sim, cS_common_sim.aD_new);
            kPpsHistM_out = zeros(nSim_sim, cS_common_sim.aD_new);
            cHistM_out = zeros(nSim_sim, cS_common_sim.aD_new);
            cppsHistM_out = zeros(nSim_sim, cS_common_sim.aD_new); % V8: PPS缴费路径

            leGridV_col_sim = paramS_sim_household.leGridV(:);

            % 创建策略函数插值器
            kPolInterp_sim = cell(cS_common_sim.nw, cS_common_sim.aD_new);
            cPpsPolInterp_choice_sim = cell(cS_common_sim.nw, cS_common_sim.aD_new); % V8
            cPolqInterp_sim = cell(cS_common_sim.nw, cS_common_sim.aD_new);

            for ia_interp = 1:cS_common_sim.aD_new
                for ie_interp = 1:cS_common_sim.nw
                    if cS_common_sim.nk > 1 && cS_common_sim.nkpps > 1
                        kPolInterp_sim{ie_interp,ia_interp} = griddedInterpolant({cS_common_sim.kGridV, cS_common_sim.kppsGridV}, squeeze(kPolM_4D_input(:,:,ie_interp,ia_interp)), 'spline', 'spline');
                        cPpsPolInterp_choice_sim{ie_interp,ia_interp} = griddedInterpolant({cS_common_sim.kGridV, cS_common_sim.kppsGridV}, squeeze(cPpsPolM_choice_4D_input(:,:,ie_interp,ia_interp)), 'spline', 'spline'); % V8
                        cPolqInterp_sim{ie_interp,ia_interp} = griddedInterpolant({cS_common_sim.kGridV, cS_common_sim.kppsGridV}, squeeze(cPolM_consump_q_4D_input(:,:,ie_interp,ia_interp)), 'spline', 'spline');
                    elseif cS_common_sim.nk > 1 && cS_common_sim.nkpps == 1 % nk > 1, nkpps = 1
                        kPolInterp_sim{ie_interp,ia_interp} = griddedInterpolant(cS_common_sim.kGridV, squeeze(kPolM_4D_input(:,1,ie_interp,ia_interp)), 'spline', 'spline');
                        cPpsPolInterp_choice_sim{ie_interp,ia_interp} = griddedInterpolant(cS_common_sim.kGridV, squeeze(cPpsPolM_choice_4D_input(:,1,ie_interp,ia_interp)), 'spline', 'spline'); % V8
                        cPolqInterp_sim{ie_interp,ia_interp} = griddedInterpolant(cS_common_sim.kGridV, squeeze(cPolM_consump_q_4D_input(:,1,ie_interp,ia_interp)), 'spline', 'spline');
                    elseif cS_common_sim.nk == 1 && cS_common_sim.nkpps > 1 % nk = 1, nkpps > 1
                        kPolInterp_sim{ie_interp,ia_interp} = griddedInterpolant(cS_common_sim.kppsGridV, squeeze(kPolM_4D_input(1,:,ie_interp,ia_interp))', 'spline', 'spline');
                        cPpsPolInterp_choice_sim{ie_interp,ia_interp} = griddedInterpolant(cS_common_sim.kppsGridV, squeeze(cPpsPolM_choice_4D_input(1,:,ie_interp,ia_interp))', 'spline', 'spline'); % V8
                        cPolqInterp_sim{ie_interp,ia_interp} = griddedInterpolant(cS_common_sim.kppsGridV, squeeze(cPolM_consump_q_4D_input(1,:,ie_interp,ia_interp))', 'spline', 'spline');
                    elseif cS_common_sim.nk == 1 && cS_common_sim.nkpps == 1 % nk = 1, nkpps = 1
                        kPolInterp_sim{ie_interp,ia_interp} = @(x,y) kPolM_4D_input(1,1,ie_interp,ia_interp);
                        cPpsPolInterp_choice_sim{ie_interp,ia_interp} = @(x,y) cPpsPolM_choice_4D_input(1,1,ie_interp,ia_interp); % V8
                        cPolqInterp_sim{ie_interp,ia_interp} = @(x,y) cPolM_consump_q_4D_input(1,1,ie_interp,ia_interp);
                    else
                        error('HHSimulation_olgm (V8): nk 或 nkpps 为零，无法创建插值器。');
                    end
                end
            end

            % PPS年度回报因子
            pps_return_net_factor_sim = 1 + ((R_k_net_factor_hh_sim - 1) + cS_common_sim.pps_return_rate_premium);

            % 按年龄组模拟
            for a_group_idx = 1:cS_common_sim.aD_new
                kNowV_group_sim = kHistM_out(:, a_group_idx);
                kPpsNowV_group_sim = kPpsHistM_out(:, a_group_idx);
                
                % 判断工作vs退休状态（基于年龄组）
                is_working_age_group = (a_group_idx <= cS_common_sim.aR_new);
                is_retired_group = ~is_working_age_group;
                
                % PPS提取判断（基于年龄组对应的第一个年度年龄）
                annual_age_check = cS_common_sim.physAgeMap{a_group_idx}(1);
                is_pps_withdrawal_eligible = (is_retired_group && ...
                    cS_common_sim.pps_active && ...
                    annual_age_check >= cS_common_sim.pps_withdrawal_age_min_idx);

                % 计算PPS提取
                pps_withdrawal_pretax_this_period = zeros(nSim_sim, 1);
                if is_pps_withdrawal_eligible
                    pps_withdrawal_pretax_this_period = kPpsNowV_group_sim * cS_common_sim.pps_withdrawal_rate;
                end
                
                % 为HHIncome_Huggett设置PPS提取参数
                paramS_hh_step = paramS_sim_household;
                paramS_hh_step.current_pps_withdrawal = 0;
                paramS_hh_step.pps_tax_rate_withdrawal = 0.15; % 默认提取税率
                if isfield(cS_common_sim, 'pps_tax_rate_withdrawal')
                    paramS_hh_step.pps_tax_rate_withdrawal = cS_common_sim.pps_tax_rate_withdrawal;
                end
                paramS_hh_step.tau_k = 0.2; % 默认资本所得税率
                if isfield(cS_common_sim, 'tau_k')
                    paramS_hh_step.tau_k = cS_common_sim.tau_k;
                end
                
                % 初始化决策变量
                kNextNonPpsV_from_policy = zeros(nSim_sim, 1);
                cPpsDecisionFromPol_choice = zeros(nSim_sim, 1); % V8: From optimal choice policy
                cConsumpValV_q_from_policy = zeros(nSim_sim, 1);

                % 按效率状态循环
                for ie_sim_idx = 1:cS_common_sim.nw
                    simIdxV_for_this_e = find(eIdxM_group(:, a_group_idx) == ie_sim_idx);
                    if isempty(simIdxV_for_this_e), continue; end

                    kNow_clamped = max(cS_common_sim.kGridV(1), min(cS_common_sim.kGridV(end), kNowV_group_sim(simIdxV_for_this_e)));
                    kPpsNow_clamped = max(cS_common_sim.kppsGridV(1), min(cS_common_sim.kppsGridV(end), kPpsNowV_group_sim(simIdxV_for_this_e)));

                    % 使用策略函数插值器
                    if cS_common_sim.nk > 1 && cS_common_sim.nkpps > 1
                        kNextNonPpsV_from_policy(simIdxV_for_this_e) = kPolInterp_sim{ie_sim_idx, a_group_idx}(kNow_clamped, kPpsNow_clamped);
                        cPpsDecisionFromPol_choice(simIdxV_for_this_e) = cPpsPolInterp_choice_sim{ie_sim_idx, a_group_idx}(kNow_clamped, kPpsNow_clamped); % V8
                        cConsumpValV_q_from_policy(simIdxV_for_this_e) = cPolqInterp_sim{ie_sim_idx, a_group_idx}(kNow_clamped, kPpsNow_clamped);
                    elseif cS_common_sim.nk > 1 && cS_common_sim.nkpps == 1
                        kNextNonPpsV_from_policy(simIdxV_for_this_e) = kPolInterp_sim{ie_sim_idx, a_group_idx}(kNow_clamped);
                        cPpsDecisionFromPol_choice(simIdxV_for_this_e) = cPpsPolInterp_choice_sim{ie_sim_idx, a_group_idx}(kNow_clamped); % V8
                        cConsumpValV_q_from_policy(simIdxV_for_this_e) = cPolqInterp_sim{ie_sim_idx, a_group_idx}(kNow_clamped);
                    elseif cS_common_sim.nk == 1 && cS_common_sim.nkpps > 1
                        kNextNonPpsV_from_policy(simIdxV_for_this_e) = kPolInterp_sim{ie_sim_idx, a_group_idx}(kPpsNow_clamped);
                        cPpsDecisionFromPol_choice(simIdxV_for_this_e) = cPpsPolInterp_choice_sim{ie_sim_idx, a_group_idx}(kPpsNow_clamped); % V8
                        cConsumpValV_q_from_policy(simIdxV_for_this_e) = cPolqInterp_sim{ie_sim_idx, a_group_idx}(kPpsNow_clamped);
                    else % nk=1, nkpps=1
                        kNextNonPpsV_from_policy(simIdxV_for_this_e) = kPolInterp_sim{ie_sim_idx, a_group_idx}(kNow_clamped(1), kPpsNow_clamped(1));
                        cPpsDecisionFromPol_choice(simIdxV_for_this_e) = cPpsPolInterp_choice_sim{ie_sim_idx, a_group_idx}(kNow_clamped(1), kPpsNow_clamped(1));
                        cConsumpValV_q_from_policy(simIdxV_for_this_e) = cPolqInterp_sim{ie_sim_idx, a_group_idx}(kNow_clamped(1), kPpsNow_clamped(1));
                    end
                    
                    % 为每个个体设置PPS提取金额
                    if is_pps_withdrawal_eligible
                        for idx_individual = 1:length(simIdxV_for_this_e)
                            sim_idx = simIdxV_for_this_e(idx_individual);
                            paramS_hh_step.current_pps_withdrawal = pps_withdrawal_pretax_this_period(sim_idx);
                        end
                    end
                end

                % 处理PPS缴费和账户更新
                actual_cpps_final_for_period_sim = cPpsDecisionFromPol_choice;
                
                if cS_common_sim.pps_active
                    % 检查是否符合PPS缴费资格（基于年龄组和对应的年度年龄）
                    can_contribute_pps_group = is_working_age_group && ...
                                               (annual_age_check <= cS_common_sim.pps_contribution_age_max_idx);
                    
                    if ~can_contribute_pps_group
                        actual_cpps_final_for_period_sim(:) = 0; % 如果不符合缴费资格，强制为0
                    end

                    % 更新PPS账户余额（按年龄组回报率）
                    kPpsNextV_group_sim = (kPpsNowV_group_sim + actual_cpps_final_for_period_sim - pps_withdrawal_pretax_this_period) * pps_return_net_factor_sim;
                    kPpsNextV_group_sim = max(cS_common_sim.kppsMin, min(cS_common_sim.kppsMax, kPpsNextV_group_sim));
                else
                    kPpsNextV_group_sim = zeros(nSim_sim, 1); % PPS未激活时账户余额为0
                end

                % 存储当期结果
                cHistM_out(:, a_group_idx) = max(cS_common_sim.cFloor, cConsumpValV_q_from_policy);
                cppsHistM_out(:, a_group_idx) = actual_cpps_final_for_period_sim; % V8: 保存PPS缴费路径
                
                % 更新下一期的资产状态
                if a_group_idx < cS_common_sim.aD_new
                    kHistM_out(:, a_group_idx + 1) = max(cS_common_sim.kMin, min(cS_common_sim.kMax, kNextNonPpsV_from_policy));
                    kPpsHistM_out(:, a_group_idx + 1) = kPpsNextV_group_sim;
                end
            end
        end

        % =====================================================================
        % == 一般均衡求解器：核心宏观经济算法 ===
        % =====================================================================
        
        function [K_sol_out, tau_l_sol_out, gbc_res_final_out, converged_and_feasible_out, solution_details_out] = solve_K_tau_l_for_rho_prime(...
                rho_prime_payg_target_input, K_init_guess_input, cS_global, paramS_global_in, eIdxM_global_sim_paths)
            % solve_K_tau_l_for_rho_prime - V8模型的一般均衡求解器
            % (MODIFIED FOR PPS BEQUESTS)
            
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
                fprintf('  IterKTL | K_guess  | tau_l_gs | MPL_g    | theta_act| K_tot_mod| K_pps_mod| GBC_res  | K_dev    | tau_l_dev| Norm     | Improv   | Strikes  | Time (s) |\n');
                fprintf('  -----------------------------------------------------------------------------------------------------------------------------------------------------------\n');
            end
            
            MPL_gross_iter_val = NaN; R_mkt_gross_factor_iter_val = NaN; theta_payg_actual_iter_val = NaN; b_payg_iter_val = NaN;
            T_bequest_model_iter_val = NaN; C_model_iter_val = NaN; Y_for_gbc_iter_val = NaN; gbc_residual_iter_val = Inf;
            K_model_from_sim_iter_val = K_current_guess; 
            K_dev_from_sim_iter_val = Inf;
            K_model_pps_sim_iter_val = NaN; 
            K_model_nonpps_sim_iter_val = NaN;

            iter_timer_start = tic; 

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

                TR_total_for_vfi_guess_val = 0.01 * MPL_gross_iter_val; % Initial guess for bequests
                if iter_ktl_idx > 1 && isfinite(T_bequest_model_iter_val)
                    TR_total_for_vfi_guess_val = T_bequest_model_iter_val;
                end
                max_vfi_tr_sub_iter = 5;
                tol_vfi_tr_sub_iter = 1e-3 * (MPL_gross_iter_val + 1e-9); % Scale tolerance
                cPolM_4D_from_vfi_final = []; kPolM_4D_from_vfi_final = []; 
                cPpsPolM_4D_choice_from_vfi_final = []; 
                TR_total_for_vfi_final_iter = TR_total_for_vfi_guess_val;

                for i_vfi_tr_sub_loop = 1:max_vfi_tr_sub_iter
                    [cPolM_vfi_temp, kPolM_vfi_temp, cPpsPolM_vfi_temp_choice, ~] = ... 
                        main_olg_v8_utils.HHSolution_VFI_Huggett(R_k_net_hh_factor_iter_val, MPL_gross_iter_val, ...
                        TR_total_for_vfi_guess_val, bV_payg_vec_iter_val, paramS_for_vfi_sim_iter, cS_global);
                    
                    % HHSimulation_olgm returns kHistM (non-PPS) and kPpsHistM (PPS) - 按年龄组模拟
                    [kHistM_sim_for_bequest, kPpsHistM_sim_for_bequest,cHistM_out, cppsHistM_out] = ...
                        main_olg_v8_utils.HHSimulation_olgm(kPolM_vfi_temp, cPpsPolM_vfi_temp_choice, cPolM_vfi_temp, ...
                        eIdxM_global_sim_paths, R_k_net_hh_factor_iter_val, MPL_gross_iter_val, ...
                        TR_total_for_vfi_guess_val, bV_payg_vec_iter_val, paramS_for_vfi_sim_iter, cS_global);
                    
                    % 🔧 计算年龄组死亡质量（使用年龄组死亡率）
                    ageDeathMass_group_iter_val = zeros(cS_global.aD_new, 1);
                    for a_death_idx = 1:cS_global.aD_new
                        group_mass = paramS_global_in.ageMassV(a_death_idx);
                        group_death_rate = cS_global.d_new(a_death_idx);  % 🔧 直接使用年龄组死亡率
                        ageDeathMass_group_iter_val(a_death_idx) = group_mass * group_death_rate;
                    end
                    
                    % Calculate non-PPS bequests (按年龄组)
                    mean_non_pps_bequest_wealth = mean(kHistM_sim_for_bequest, 1);
                    TotalNonPPSBequests_pc_iter_val = sum(mean_non_pps_bequest_wealth(:) .* ageDeathMass_group_iter_val(:));

                    % Calculate PPS bequests (if PPS is active and bequeathable) (按年龄组)
                    TotalPPSBequests_pc_iter_val = 0;
                    if cS_global.pps_active && cS_global.pps_bequeathable
                        if ~isempty(kPpsHistM_sim_for_bequest)
                            mean_pps_bequest_wealth = mean(kPpsHistM_sim_for_bequest, 1);
                            TotalPPSBequests_pc_iter_val = sum(mean_pps_bequest_wealth(:) .* ageDeathMass_group_iter_val(:));
                        end
                    end
                    
                    % Total bequests
                    TotalBequests_pc_iter_val = TotalNonPPSBequests_pc_iter_val + TotalPPSBequests_pc_iter_val;
                    
                    T_bequest_model_new_iter_val = TotalBequests_pc_iter_val / (1 + paramS_global_in.popGrowthForDebt);
                    T_bequest_model_new_iter_val = max(0, T_bequest_model_new_iter_val);

                    T_bequest_model_iter_val = T_bequest_model_new_iter_val; % Update for this iteration
                    cPolM_4D_from_vfi_final = cPolM_vfi_temp;
                    kPolM_4D_from_vfi_final = kPolM_vfi_temp;
                    cPpsPolM_4D_choice_from_vfi_final = cPpsPolM_vfi_temp_choice; 
                    TR_total_for_vfi_final_iter = T_bequest_model_new_iter_val; % This will be used for the main simulation

                    if abs(T_bequest_model_new_iter_val - TR_total_for_vfi_guess_val) < tol_vfi_tr_sub_iter || ...
                       i_vfi_tr_sub_loop == max_vfi_tr_sub_iter
                        break; % Converged or max iterations for bequest loop
                    end
                    TR_total_for_vfi_guess_val = 0.5 * TR_total_for_vfi_guess_val + 0.5 * T_bequest_model_new_iter_val; % Dampening
                end
                T_bequest_model_iter_val = TR_total_for_vfi_final_iter; % Final bequest value for this K, tau_l iteration

                [kHistM_non_pps_sim_iter_val, kPpsHistM_sim_iter_val, cHistM_sim_iter_val, ~] = main_olg_v8_utils.HHSimulation_olgm(...
                    kPolM_4D_from_vfi_final, cPpsPolM_4D_choice_from_vfi_final, cPolM_4D_from_vfi_final, eIdxM_global_sim_paths, ...
                    R_k_net_hh_factor_iter_val, MPL_gross_iter_val, TR_total_for_vfi_final_iter, bV_payg_vec_iter_val, paramS_for_vfi_sim_iter, cS_global);

                % 使用年龄组权重进行聚合计算
                K_model_nonpps_sim_iter_val = mean(kHistM_non_pps_sim_iter_val, 1) * paramS_global_in.ageMassV;
                K_model_pps_sim_iter_val = 0;
                if cS_global.pps_active && cS_global.pps_in_K && (cS_global.pps_max_contrib_frac > 0 || cS_global.pps_contrib_limit > 0)
                     if ~isempty(kPpsHistM_sim_iter_val)
                        K_model_pps_sim_iter_val = mean(kPpsHistM_sim_iter_val, 1) * paramS_global_in.ageMassV;
                        K_model_pps_sim_iter_val = max(0, K_model_pps_sim_iter_val);
                     end
                end
                K_model_from_sim_iter_val = K_model_nonpps_sim_iter_val + K_model_pps_sim_iter_val;
                K_model_from_sim_iter_val = max(1e-6, K_model_from_sim_iter_val);
                C_model_iter_val = mean(cHistM_sim_iter_val,1) * paramS_global_in.ageMassV;

                Y_for_gbc_iter_val = cS_global.A * (K_current_guess^cS_global.alpha) * (L_per_capita_global^(1-cS_global.alpha));
                G_val_target_iter_val = cS_global.gov_exp_frac_Y * Y_for_gbc_iter_val;
                B_val_target_iter_val = cS_global.gov_debt_frac_Y * Y_for_gbc_iter_val;

                gbc_residual_iter_val = main_olg_v8_utils.check_gbc_residual(K_current_guess, C_model_iter_val, Y_for_gbc_iter_val, ...
                    G_val_target_iter_val, B_val_target_iter_val, MPL_gross_iter_val, r_mkt_gross_iter_val, ...
                    theta_payg_actual_iter_val, tau_l_current_guess, ...
                    b_payg_iter_val, T_bequest_model_iter_val, 0, cS_global, paramS_global_in);

                K_dev_from_sim_iter_val = K_current_guess - K_model_from_sim_iter_val;
                tau_l_dev_raw_for_update = -gbc_residual_iter_val / (MPL_gross_iter_val * L_per_capita_global + 1e-9); % Avoid division by zero
                current_devNorm_val = sqrt(K_dev_from_sim_iter_val^2 + (gbc_residual_iter_val)^2 );
                norm_improvement_val = prev_devNorm_ktl - current_devNorm_val;
                
                elapsed_iter_time = toc(iter_timer_start); 
                
                fprintf('  %7d | %8.4f | %8.4f | %8.4f | %8.4f | %8.4f | %8.4f | %8.2e | %8.4f | %8.4f | %8.2e | %.1e | %7d | %8.2f |\n', ...
                    iter_ktl_idx, K_current_guess, tau_l_current_guess, MPL_gross_iter_val, theta_payg_actual_iter_val, ...
                    K_model_from_sim_iter_val, K_model_pps_sim_iter_val, gbc_residual_iter_val, ...
                    K_dev_from_sim_iter_val, tau_l_dev_raw_for_update, current_devNorm_val, norm_improvement_val, tau_l_boundary_strike_count_ktl, elapsed_iter_time);
                
                iter_timer_start = tic; 

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
                K_current_guess = max(1e-3, K_guess_next_val); % Ensure K remains positive
                new_tau_l_unconstrained_val = tau_l_current_guess + damp_tau_l_ktl_loop * tau_l_dev_raw_for_update;
                tau_l_next_iter_constrained_val = max(cS_global.tau_l_min, min(cS_global.tau_l_max, new_tau_l_unconstrained_val));
                
                is_tau_l_at_boundary_now = ...
                    (abs(tau_l_next_iter_constrained_val - cS_global.tau_l_max) < 1e-7 && new_tau_l_unconstrained_val >= cS_global.tau_l_max - 1e-7) || ...
                    (abs(tau_l_next_iter_constrained_val - cS_global.tau_l_min) < 1e-7 && new_tau_l_unconstrained_val <= cS_global.tau_l_min + 1e-7);

                if is_tau_l_at_boundary_now && abs(gbc_residual_iter_val) > cS_global.gbc_tol_for_internal_loop
                    tau_l_boundary_strike_count_ktl = tau_l_boundary_strike_count_ktl + 1;
                else
                    tau_l_boundary_strike_count_ktl = 0; % Reset if not at boundary or GBC is met
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
            end % End of K and tau_l iteration loop

            if ~converged_and_feasible_out && iter_ktl_idx == maxIter_ktl_loop
                fprintf('  警告 (V8): K和tau_l迭代达到最大次数 (%d) 或在该次数内未达可行解 (rho_prime_target=%.4f).\n', maxIter_ktl_loop, rho_prime_payg_target_input);
                K_sol_out = K_model_from_sim_iter_val; % Store last computed values
                tau_l_sol_out = tau_l_current_guess;
                gbc_res_final_out = gbc_residual_iter_val;
                % Ensure all solution_details fields are populated even if convergence fails
                if exist('MPL_gross_iter_val', 'var') && isfinite(MPL_gross_iter_val)
                    solution_details_out.R_mkt_gross_factor = R_mkt_gross_factor_iter_val;
                    solution_details_out.MPL_gross = MPL_gross_iter_val;
                    solution_details_out.theta_payg = theta_payg_actual_iter_val;
                    solution_details_out.b_payg = b_payg_iter_val;
                    solution_details_out.T_bequest_Model = T_bequest_model_iter_val;
                    solution_details_out.C_model = C_model_iter_val;
                    solution_details_out.Y_model = Y_for_gbc_iter_val; % Use Y based on K_current_guess for consistency if K_sol_out is from sim
                    solution_details_out.K_model_pps = K_model_pps_sim_iter_val;
                    solution_details_out.K_model_non_pps = K_model_nonpps_sim_iter_val;
                else % If loop didn't even run once or values are NaN
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

            % Ensure theta_payg_required_before_cap is always in details
            if ~isfield(solution_details_out, 'theta_payg_required_before_cap') || ...
               (isfield(solution_details_out, 'theta_payg_required_before_cap') && isnan(solution_details_out.theta_payg_required_before_cap))
                 recalc_mass_retirees_final_sd = sum(paramS_global_in.ageMassV(cS_global.aR_new + 1 : cS_global.aD_new));
                 recalc_theta_req_final_sd = 0;
                 if mass_workers_global > 1e-9
                     recalc_theta_req_final_sd = rho_prime_payg_target_input * (recalc_mass_retirees_final_sd / mass_workers_global);
                 end
                 solution_details_out.theta_payg_required_before_cap = max(0, recalc_theta_req_final_sd);
            end
             if ~isfield(solution_details_out, 'K_model_pps') % Ensure K_pps is in details
                if exist('K_model_pps_sim_iter_val','var') && isfinite(K_model_pps_sim_iter_val)
                    solution_details_out.K_model_pps = K_model_pps_sim_iter_val;
                else
                    solution_details_out.K_model_pps = NaN;
                end
            end
             if ~isfield(solution_details_out, 'K_model_non_pps') % Ensure K_non_pps is in details
                 if exist('K_model_nonpps_sim_iter_val','var') && isfinite(K_model_nonpps_sim_iter_val)
                     solution_details_out.K_model_non_pps = K_model_nonpps_sim_iter_val;
                 else
                     solution_details_out.K_model_non_pps = NaN;
                 end
            end
            % If converged, Y_model in details should be based on K_sol_out (which is K_model_from_sim)
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
            
            % Revenue from labor tax (general part, not PAYG)
            % This needs to consider that PPS contributions might be tax-deductible for tau_l
            % However, for a macro check, MPL * L is often used as the broad base.
            % A more precise calculation would require aggregate PPS deductions.
            % For now, using the simpler broad base, consistent with how tau_l_dev is often calculated.
            LaborTaxRev_general_part_calc = tau_l_val_input * MPL_gross_val_input * L_per_capita_local_check;
            % If PPS deductions were significant and affected the aggregate tax base for tau_l:
            % effective_labor_income_base_for_tau_l = MPL_gross_val_input * L_per_capita_local_check - Aggregate_PPS_Deductions_pc;
            % LaborTaxRev_general_part_calc = tau_l_val_input * effective_labor_income_base_for_tau_l;

            CapitalTaxRev_calc = r_mkt_gross_val_input * K_val_market_input * cS_check.tau_k;
            ConsumptionTaxRev_calc = C_val_model_input * cS_check.tau_c;
            GeneralRevenue_calc = LaborTaxRev_general_part_calc + CapitalTaxRev_calc + ConsumptionTaxRev_calc;

            GovConsumption_calc = G_val_target_input;
            r_b_for_debt_service_calc = r_mkt_gross_val_input; % Assuming gov debt pays market rate
            DebtService_calc = (r_b_for_debt_service_calc - paramS_loc_check.popGrowthForDebt) * B_val_target_input;
            GovDirectTransfers_calc = TR_gov_val_pc_input; % This is lump-sum from gov, usually zero in this model type
            
            GeneralOutlays_calc = GovConsumption_calc + DebtService_calc + GovDirectTransfers_calc;
            
            % Residual: Revenue - Outlays. Positive means surplus.
            gbc_residual_out = GeneralRevenue_calc - GeneralOutlays_calc;
        end

        % --- CallInterpolator (与V7/Baseline一致) ---
        function ev_val_out = CallInterpolator(interpolant_obj_input, k_val_input, k_pps_val_input, cS_local_interp)
            % (Identical to main_olg_baseline_utils.m)
            ev_val_out = -Inf; % Default for errors or out-of-bound
            try
                if isa(interpolant_obj_input, 'griddedInterpolant')
                    if cS_local_interp.nk > 1 && cS_local_interp.nkpps > 1
                        ev_val_out = interpolant_obj_input(k_val_input, k_pps_val_input);
                    elseif cS_local_interp.nk > 1 % nkpps must be 1
                        ev_val_out = interpolant_obj_input(k_val_input);
                    elseif cS_local_interp.nkpps > 1 % nk must be 1
                        ev_val_out = interpolant_obj_input(k_pps_val_input);
                    else % nk=1, nkpps=1
                        if isscalar(interpolant_obj_input.Values)
                            ev_val_out = interpolant_obj_input.Values; % The single value
                        else
                             % This case should ideally not happen if nk=1, nkpps=1,
                             % unless interpolant was built incorrectly.
                             % Try to evaluate as if it were 2D, though it might error or give an unexpected result.
                            ev_val_out = interpolant_obj_input(k_val_input, k_pps_val_input);
                        end
                    end
                elseif isa(interpolant_obj_input, 'function_handle') % For nk=1, nkpps=1 case with lambda
                    ev_val_out = interpolant_obj_input(k_val_input, k_pps_val_input);
                else
                    % warning('CallInterpolator (V8): 未处理的插值器类型。');
                    ev_val_out = -1e11; % Large negative value for error
                end
            catch ME_call_interp_error
                % warning('CallInterpolator (V8): 插值过程中发生错误: "%s"。将尝试限制输入值并重试。', ME_call_interp_error.message);
                % Clamp inputs to grid boundaries and retry
                k_clamped = k_val_input; 
                if cS_local_interp.nk > 0 && ~isempty(cS_local_interp.kGridV)
                    k_clamped = max(cS_local_interp.kGridV(1), min(cS_local_interp.kGridV(end), k_val_input));
                end
                
                k_pps_clamped = k_pps_val_input; 
                if cS_local_interp.nkpps > 0 && ~isempty(cS_local_interp.kppsGridV)
                    k_pps_clamped = max(cS_local_interp.kppsGridV(1), min(cS_local_interp.kppsGridV(end), k_pps_val_input));
                end
                
                try % Retry with clamped values
                    if isa(interpolant_obj_input, 'griddedInterpolant')
                        if cS_local_interp.nk > 1 && cS_local_interp.nkpps > 1
                            ev_val_out = interpolant_obj_input(k_clamped, k_pps_clamped);
                        elseif cS_local_interp.nk > 1
                            ev_val_out = interpolant_obj_input(k_clamped);
                        elseif cS_local_interp.nkpps > 1
                            ev_val_out = interpolant_obj_input(k_pps_clamped);
                        else % nk=1, nkpps=1
                             if isscalar(interpolant_obj_input.Values)
                                ev_val_out = interpolant_obj_input.Values;
                            else
                                ev_val_out = interpolant_obj_input(k_clamped, k_pps_clamped);
                            end
                        end
                    elseif isa(interpolant_obj_input, 'function_handle')
                        ev_val_out = interpolant_obj_input(k_clamped, k_pps_clamped);
                    end
                catch % If retry also fails
                    ev_val_out = -1e11; % Large negative value
                end
            end
            if ~isfinite(ev_val_out) % Final check for NaN or Inf
                ev_val_out = -1e12; % Ensure it's a very bad outcome
            end
        end
        
        % --- V8 VFI 连续优化版：基于fmincon的联合优化 ---
        function [cPol_age_q, kPol_age, cPpsPol_age_choice, val_age] = HHSolutionByAge_VFI_Huggett_v8_fmincon(...
                a_idx, vPrime_kkppse_next, R_k_net_factor_age, w_gross_age, TR_total_age, b_age_val, ...
                paramS_age, cS, epsilon_grid)
            % HHSolutionByAge_VFI_Huggett_v8_fmincon - 基于fmincon的连续优化VFI算法
            % (Code from your provided version, minor adjustments for robustness if any)
            
            cPol_age_q_init = zeros(cS.nk, cS.nkpps, cS.nw);
            kPol_age_init   = zeros(cS.nk, cS.nkpps, cS.nw);
            cPpsPol_age_choice_init = zeros(cS.nk, cS.nkpps, cS.nw);
            val_age_init    = -Inf(cS.nk, cS.nkpps, cS.nw);

            fmincon_opts = optimoptions('fmincon', ...
                'Display', 'none', ...
                'Algorithm', 'interior-point', ...
                'SpecifyObjectiveGradient', false, ... % 我们不提供解析梯度
                'FiniteDifferenceType', 'central', ...  % 使用中心差分
                'TolFun', 1e-7, ... % Slightly tighter
                'TolX', 1e-7, ...   % Slightly tighter
                'MaxIter', 500, ... % Increased iterations
                'MaxFunEvals', 2000); % Increased evals

            % fmincon_opts = optimoptions('fmincon', ...
            %     'Display', 'none', ...
            %     'Algorithm', 'sqp', ... % or 'interior-point' which can also use gradients
            %     'SpecifyObjectiveGradient', true, ... % <<<< TELL FMINCON TO USE GRADIENT >>>>
        
            %     'TolFun', 1e-7, ...
            %     'TolX', 1e-7, ...
            %     'MaxIter', 500, ...
            %     'MaxFunEvals', 2000);

            if a_idx == cS.aD_new
                [K_grid, Kpps_grid, Epsilon_ndgrid] = ndgrid(cS.kGridV, cS.kppsGridV, epsilon_grid);
                resources_batch = zeros(size(K_grid));
                % Can use parfor here if numel(K_grid) is large enough
                for i_nd = 1:numel(K_grid) 
                    [ik_nd, ikpps_nd, ie_nd] = ind2sub(size(K_grid), i_nd);
                    
                    % 为最后一期设置PPS提取参数
                    paramS_last_period = paramS_age;
                    paramS_last_period.current_pps_withdrawal = 0;
                    paramS_last_period.pps_tax_rate_withdrawal = cS.pps_tax_rate_withdrawal;
                    paramS_last_period.tau_k = cS.tau_k;
                    
                    % 如果是退休期且PPS激活，计算PPS提取
                    if cS.pps_active && a_idx > cS.aR_new
                        paramS_last_period.current_pps_withdrawal = Kpps_grid(ik_nd, ikpps_nd, ie_nd) * cS.pps_withdrawal_rate;
                    end
                    
                    [resources_batch(ik_nd, ikpps_nd, ie_nd), ~, ~] = main_olg_v8_utils.HHIncome_Huggett(...
                        K_grid(ik_nd, ikpps_nd, ie_nd), R_k_net_factor_age, w_gross_age, TR_total_age, b_age_val, ...
                        0, a_idx, paramS_last_period, cS, Epsilon_ndgrid(ik_nd, ikpps_nd, ie_nd));
                end
                total_resources = resources_batch;
                cPol_age_q_init = max(cS.cFloor, total_resources / (1 + cS.tau_c));
                kPol_age_init(:) = cS.kMin;
                cPpsPol_age_choice_init(:) = 0;
                util_temp_storage = zeros(size(cPol_age_q_init));
                % Can use parfor here
                for i_nd_util = 1:numel(cPol_age_q_init) 
                    [~, util_temp_storage(i_nd_util)] = main_olg_v8_utils.CES_utility(cPol_age_q_init(i_nd_util), cS.sigma, cS);
                end
                val_age_init = util_temp_storage;
            else 
                EV_matrix = zeros(cS.nk, cS.nkpps, cS.nw);
                for ie_current = 1:cS.nw
                    transition_probs = paramS_age.leTrProbM(ie_current, :);
                    EV_slice = sum(vPrime_kkppse_next .* reshape(transition_probs, 1, 1, cS.nw), 3);
                    EV_matrix(:, :, ie_current) = EV_slice;
                end
                
                EV_interpolants = cell(cS.nw, 1);
                for ie_current = 1:cS.nw
                    if cS.nk > 1 && cS.nkpps > 1
                        EV_interpolants{ie_current} = griddedInterpolant({cS.kGridV, cS.kppsGridV}, ...
                            EV_matrix(:, :, ie_current), 'spline', 'spline'); % Use spline for smoother derivatives
                    elseif cS.nk > 1 && cS.nkpps == 1
                        EV_interpolants{ie_current} = griddedInterpolant(cS.kGridV, ...
                            EV_matrix(:, 1, ie_current), 'spline', 'spline');
                    elseif cS.nk == 1 && cS.nkpps > 1
                         EV_interpolants{ie_current} = griddedInterpolant(cS.kppsGridV, ...
                            squeeze(EV_matrix(1, :, ie_current)), 'spline', 'spline');
                    else % nk=1, nkpps=1
                        EV_interpolants{ie_current} = @(k_s, kp_s) EV_matrix(1, 1, ie_current);
                    end
                end
                
                cS_local = cS; % Broadcast cS once
                paramS_local = paramS_age;
                a_idx_local = a_idx;
                R_k_net_factor_local = R_k_net_factor_age;
                w_gross_local = w_gross_age;
                TR_total_local = TR_total_age;
                b_age_val_local = b_age_val;
                EV_interpolants_local_bcast = EV_interpolants; % Broadcast once
                fmincon_opts_local_bcast = fmincon_opts;
                epsilon_grid_local_bcast = epsilon_grid;

                nk_const = cS_local.nk;
                nkpps_const = cS_local.nkpps;
                nw_const = cS_local.nw;
                
                % Temporary matrices for parfor output
                cPol_temp = zeros(nk_const, nkpps_const, nw_const);
                kPol_temp = zeros(nk_const, nkpps_const, nw_const);
                cPpsPol_temp = zeros(nk_const, nkpps_const, nw_const);
                val_temp = -Inf(nk_const, nkpps_const, nw_const);

                parfor ik = 1:nk_const
                    % Slices for this worker
                    cPol_slice_par = zeros(nkpps_const, nw_const);
                    kPol_slice_par = zeros(nkpps_const, nw_const);
                    cPpsPol_slice_par = zeros(nkpps_const, nw_const);
                    val_slice_par = -Inf(nkpps_const, nw_const);

                    for ikpps = 1:nkpps_const
                        for ie = 1:nw_const
                            k_state = cS_local.kGridV(ik);
                            k_pps_state = cS_local.kppsGridV(ikpps);
                            epsilon_state = epsilon_grid_local_bcast(ie);
                            
                            model_age_group_start_year_idx = cS_local.physAgeMap{a_idx_local}(1);
                            is_pps_eligible = (a_idx_local <= cS_local.aR_new && ...
                                model_age_group_start_year_idx <= cS_local.pps_contribution_age_max_idx && ...
                                cS_local.pps_active);
                            
                            max_permissible_cpps = 0;
                            if is_pps_eligible
                                age_efficiency = cS_local.ageEffV_new(a_idx_local);
                                current_gross_labor_income = w_gross_local * age_efficiency * epsilon_state;
                                if current_gross_labor_income > 1e-6
                                    max_cpps_by_frac = current_gross_labor_income * cS_local.pps_max_contrib_frac;
                                    max_permissible_cpps = min(cS_local.pps_contrib_limit, max_cpps_by_frac);
                                    max_permissible_cpps = max(0, max_permissible_cpps);
                                end
                            end
     % 无梯度版本                      
                            objective_function = @(x_prop) main_olg_v8_utils.fmincon_objective_helper_proportional(...
                                x_prop, k_state, k_pps_state, epsilon_state, a_idx_local, ie, ...
                                R_k_net_factor_local, w_gross_local, TR_total_local, b_age_val_local, ...
                                paramS_local, cS_local, EV_interpolants_local_bcast, max_permissible_cpps);
    % 有梯度版本
    % objective_function = @(x_prop) main_olg_v8_utils.fmincon_objective_helper_proportional_with_deri(...
    % x_prop, k_state, k_pps_state, epsilon_state, a_idx_local, ie, ...
    % R_k_net_factor_local, w_gross_local, TR_total_local, b_age_val_local, ...
    % paramS_local, cS_local, EV_interpolants_local_bcast, max_permissible_cpps);

                            lb_fmin = [0, 0];
                            ub_fmin = [1, 1];
                            x0_pps_prop_fmin = 0.5;
                            if max_permissible_cpps < 1e-9 
                                x0_pps_prop_fmin = 0;
                                ub_fmin(1) = 0; 
                            end
                            x0_fmin = [x0_pps_prop_fmin, 0.5]; % Start with 50% of max pps, 50% of spendable on k'

                            optimal_cpps_val = 0;
                            optimal_k_prime_val = cS_local.kMin;
                            optimal_c_val = cS_local.cFloor;
                            optimal_value_val = -Inf;


                                [x_opt_prop, fval, exitflag] = fmincon(objective_function, x0_fmin, [], [], [], [], ...
                                    lb_fmin, ub_fmin, [], fmincon_opts_local_bcast);
                                
                                    pps_prop_opt = x_opt_prop(1);
                                    k_prime_prop_opt = x_opt_prop(2);
                                    optimal_value_val = -fval;

                                    optimal_cpps_val = pps_prop_opt * max_permissible_cpps;
                                    optimal_cpps_val = max(0, min(optimal_cpps_val, max_permissible_cpps)); % Re-clamp for safety

                                    [resources_after_pps_opt, ~, ~] = main_olg_v8_utils.HHIncome_Huggett(...
                                        k_state, R_k_net_factor_local, w_gross_local, TR_total_local, b_age_val_local, ...
                                        optimal_cpps_val, a_idx_local, paramS_local, cS_local, epsilon_state);

                                    consumption_floor_spending_opt = cS_local.cFloor * (1 + cS_local.tau_c);
                                    resources_for_kprime_c_above_floor_opt = resources_after_pps_opt - consumption_floor_spending_opt;
                                    
                                    if resources_for_kprime_c_above_floor_opt >=0
                                        optimal_k_prime_val = k_prime_prop_opt * resources_for_kprime_c_above_floor_opt;
                                        optimal_k_prime_val = max(cS_local.kMin, min(optimal_k_prime_val, resources_for_kprime_c_above_floor_opt));
                                    else
                                        optimal_k_prime_val = cS_local.kMin;
                                    end
                                    optimal_k_prime_val = max(cS_local.kMin, min(optimal_k_prime_val, cS_local.kMax));

                                    consumption_expenditure_opt = resources_after_pps_opt - optimal_k_prime_val;
                                    optimal_c_val = max(cS_local.cFloor, consumption_expenditure_opt / (1 + cS_local.tau_c));


                            
                            val_slice_par(ikpps, ie) = optimal_value_val;
                            cPol_slice_par(ikpps, ie) = optimal_c_val;
                            kPol_slice_par(ikpps, ie) = optimal_k_prime_val;
                            cPpsPol_slice_par(ikpps, ie) = optimal_cpps_val;
                        end
                    end
                    % Assign slices from this worker to the temporary full matrices
                    val_temp(ik, :, :) = val_slice_par;
                    cPol_temp(ik, :, :) = cPol_slice_par;
                    kPol_temp(ik, :, :) = kPol_slice_par;
                    cPpsPol_temp(ik, :, :) = cPpsPol_slice_par;
                end % end parfor ik
                % Copy from temporary storage to output matrices
                val_age_init = val_temp;
                cPol_age_q_init = cPol_temp;
                kPol_age_init = kPol_temp;
                cPpsPol_age_choice_init = cPpsPol_temp;
            end 
            
            cPol_age_q = cPol_age_q_init;
            kPol_age = kPol_age_init;
            cPpsPol_age_choice = cPpsPol_age_choice_init;
            val_age = val_age_init;
        end

        % --- fmincon辅助函数：目标函数 (使用比例决策变量) ---
        function neg_value = fmincon_objective_helper_proportional(x_prop, k_state, k_pps_state, epsilon_state, a_idx, ie, ...
                R_k_net_factor_age, w_gross_age, TR_total_age, b_age_val, paramS_age, cS, EV_interpolants, max_permissible_cpps)
            % 该函数计算家庭最优化问题的负值目标函数，用于fmincon优化
            % 输入:
            %   x_prop - 决策变量向量 [pps比例, k'比例]
            %   k_state - 当前资产状态
            %   k_pps_state - 当前养老金账户状态
            %   epsilon_state - 当前劳动生产率冲击
            %   a_idx - 年龄指数
            %   ie - 劳动生产率状态指数
            %   R_k_net_factor_age - 净资本回报率因子
            %   w_gross_age - 总工资率
            %   TR_total_age - 总转移支付
            %   b_age_val - 养老金福利
            %   paramS_age - 年龄特定参数
            %   cS - 模型参数结构体
            %   EV_interpolants - 期望值函数插值器
            %   max_permissible_cpps - 最大允许的养老金缴费

            % 从比例决策变量提取实际决策
            pps_proportion = x_prop(1);
            k_prime_proportion = x_prop(2);

            % 计算实际养老金缴费并确保在允许范围内
            actual_c_pps = pps_proportion * max_permissible_cpps;
            actual_c_pps = max(0, min(actual_c_pps, max_permissible_cpps));



            % 计算养老金提取
            pps_withdrawal = 0;
            annual_age_check = cS.physAgeMap{a_idx}(1);
            is_retired = (a_idx > cS.aR_new); % 基于模型年龄组判断是否退休
            % 检查是否满足养老金提取条件
            if is_retired && annual_age_check >= cS.pps_withdrawal_age_min_idx && cS.pps_active
                pps_withdrawal = k_pps_state * cS.pps_withdrawal_rate;
            end

            % 为HHIncome_Huggett设置PPS提取参数
            paramS_obj = paramS_age;
            paramS_obj.current_pps_withdrawal = pps_withdrawal;
            paramS_obj.pps_tax_rate_withdrawal = cS.pps_tax_rate_withdrawal;
            paramS_obj.tau_k = cS.tau_k;

            % 重新计算资源，包括PPS提取税收
            [resources_after_pps_corrected, ~, ~] = main_olg_v8_utils.HHIncome_Huggett(...
                k_state, R_k_net_factor_age, w_gross_age, TR_total_age, b_age_val, ...
                actual_c_pps, a_idx, paramS_obj, cS, epsilon_state);

            % 使用修正后的资源重新计算消费和资产
            consumption_floor_spending_corrected = cS.cFloor * (1 + cS.tau_c);
            resources_for_kprime_c_above_floor_corrected = resources_after_pps_corrected - consumption_floor_spending_corrected;
            
            if resources_for_kprime_c_above_floor_corrected >= 0
                actual_k_prime = k_prime_proportion * resources_for_kprime_c_above_floor_corrected;
                actual_k_prime = max(cS.kMin, min(actual_k_prime, resources_for_kprime_c_above_floor_corrected));
                consumption_expenditure_corrected = resources_after_pps_corrected - actual_k_prime;
                current_c = max(cS.cFloor, consumption_expenditure_corrected / (1 + cS.tau_c));
            else
                actual_k_prime = cS.kMin;
                consumption_expenditure_corrected = resources_after_pps_corrected - actual_k_prime;
                current_c = max(cS.cFloor, consumption_expenditure_corrected / (1 + cS.tau_c));
            end
            
            actual_k_prime = max(cS.kMin, min(actual_k_prime, cS.kMax));

            % 计算养老金账户回报和下期养老金账户余额
            pps_return_factor = 1 + ((R_k_net_factor_age - 1) + cS.pps_return_rate_premium);
            k_pps_prime = (k_pps_state + actual_c_pps - pps_withdrawal) * pps_return_factor;
            k_pps_prime = max(cS.kppsMin, min(cS.kppsMax, k_pps_prime));

            % 计算当期效用
            [~, current_utility] = main_olg_v8_utils.CES_utility(current_c, cS.sigma, cS);
            if ~isfinite(current_utility)
                % 如果效用无穷，返回惩罚值
                neg_value = 1e12 + abs(current_c - cS.cFloor) * 1e10;
                return;
            end

            % 计算期望未来价值
            expected_future_value = -Inf;
            if a_idx < cS.aD_new
                try
                    % 确保评估点在网格范围内
                    k_prime_eval = max(cS.kGridV(1), min(cS.kGridV(end), actual_k_prime));
                    k_pps_prime_eval = max(cS.kppsGridV(1), min(cS.kppsGridV(end), k_pps_prime));
                    expected_future_value = main_olg_v8_utils.CallInterpolator(...
                        EV_interpolants{ie}, k_prime_eval, k_pps_prime_eval, cS);
                catch
                    fprintf('Error in CallInterpolator: %s\n', lasterr);
                    % 插值失败时设置惩罚值
                    expected_future_value = -1e11;
                end
            end
            if ~isfinite(expected_future_value)
                expected_future_value = -1e11;
            end

            % 计算生存概率和总价值
            s_transition = cS.s_1yr_transitionV(a_idx);
            total_value = current_utility + cS.beta * s_transition * expected_future_value;

            % 返回负值用于最小化问题
            if ~isfinite(total_value)
                neg_value = 1e12;
            else
                neg_value = -total_value;
            end
        end
        

    % --- V8 家庭问题核心函数：基于神经网络的策略预测 ---
    
    end % 结束 Static 方法块

end % 结束 main_olg_v8_utils 类定义