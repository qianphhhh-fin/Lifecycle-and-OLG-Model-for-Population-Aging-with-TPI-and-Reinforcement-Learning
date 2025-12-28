% --- START OF FILE main_olg_v7_utils.m ---
classdef main_olg_v7_utils % 定义 main_olg_v7_utils 类
    % OLG 模型 v7 的工具函数 (VFI包含k_pps作为状态变量, PPS所得税递延, PPS缴费规则化)
    % - HHSolutionByAge_VFI_Huggett_v7 中的 ik 循环启用了 parfor
    % - 所有依赖的v4/v6函数已内联至此文件，实现自包含。

    methods (Static) % 定义静态方法

        % =====================================================================
        % == 参数设定函数 ==
        % =====================================================================
        function cS = ParameterValues_HuggettStyle()
            % ParameterValues_HuggettStyle - 设置OLG模型的所有参数
            %
            % 输出:
            %   cS - 包含所有模型参数的结构体
            %
            % 说明:
            %   此函数首先加载一个基础版本的参数集 (源自v4模型)，
            %   然后根据v7模型的特定需求进行修改和添加。
            %   包括人口统计、家庭偏好、技术、社会保障、政府财政、
            %   劳动禀赋、网格以及个人养老金计划(PPS)等参数。

            % --- 首先, 获取 v4 版本的基准参数 ---
            % 人口统计与分组参数 (与之前v3/v4版本一致)
            cS.age1_orig = 20;                         % 初始经济活动年龄 (年度)
            cS.ageLast_orig = 98;                      % 最大存活年龄 (年度)
            cS.ageRetire_orig = 65;                    % 退休年龄 (年度)
            cS.popGrowth_orig = 0.012;                 % 年化人口增长率 (基准值)
            cS.aD_orig = cS.ageLast_orig - cS.age1_orig + 1; % 年度模型总年龄跨度
            cS.aR_idx_orig = cS.ageRetire_orig - cS.age1_orig + 1; % 年度模型退休年龄对应的索引
            cS.aW_orig = cS.aR_idx_orig - 1;           % 年度模型工作期长度
            cS.physAgeV_orig = (cS.age1_orig : cS.ageLast_orig)'; % 年度生理年龄向量

            % 年度死亡率 (来自外部数据)
            cS.d_orig = [0.00159,0.00169,0.00174,0.00172,0.00165,0.00156,0.00149,0.00145,0.00145,0.00149,0.00156,0.00163,0.00171,0.00181,0.00193,0.00207,0.00225,0.00246,0.00270,0.00299,0.00332,0.00368,0.00409,0.00454,0.00504,0.00558,0.00617,0.00686,0.00766,0.00865,0.00955,0.01058,0.01162,0.01264,0.01368,0.01475,0.01593,0.01730,0.01891,0.02074,0.02271,0.02476,0.02690,0.02912,0.03143,0.03389,0.03652,0.03930,0.04225,0.04538,0.04871,0.05230,0.05623,0.06060,0.06542,0.07066,0.07636,0.08271,0.08986,0.09788,0.10732,0.11799,0.12895,0.13920,0.14861,0.16039,0.17303,0.18665,0.20194,0.21877,0.23601,0.25289,0.26973,0.28612,0.30128,0.31416,0.32915,0.34450,0.36018];
            if length(cS.d_orig) ~= cS.aD_orig
                error('年度死亡率数据 d_orig 长度与年龄跨度不匹配');
            end
            cS.s_orig = 1 - cS.d_orig; % 年度存活率

            % 时间单位转换参数 ("混合时间单位"设置)
            cS.yearStep = 5; % 每个模型年龄组代表的年数
            cS.aD_new = ceil(cS.aD_orig / cS.yearStep); % 模型中的总年龄组数量
            cS.aR_new = ceil(cS.aW_orig / cS.yearStep); % 模型中工作年龄组的数量 (注意: 不是退休年龄组索引)

            % 建立年度年龄到模型年龄组的映射
            cS.physAgeMap = cell(cS.aD_new, 1); % 初始化映射单元格数组
            for a = 1:cS.aD_new
                startIdx = (a-1)*cS.yearStep + 1; % 该年龄组对应的年度年龄起始索引
                endIdx = min(a*cS.yearStep, cS.aD_orig); % 该年龄组对应的年度年龄结束索引
                cS.physAgeMap{a} = startIdx:endIdx; % 存储映射关系
            end

            % 模型年龄组对应的代表性生理年龄 (取每组的第一个年度年龄)
            cS.physAgeV_new = zeros(cS.aD_new, 1);
            for a = 1:cS.aD_new
                cS.physAgeV_new(a) = cS.physAgeV_orig(cS.physAgeMap{a}(1));
            end

            % 模型年龄组之间的年度化转移存活率 (用于VFI中的beta贴现)
            % s_1yr_transitionV(j) 指的是从第 j 个年龄组活到该年龄组最后一年，并再活一年进入下一个年龄组的概率
            % (近似为第j组最后一个年度年龄的存活率)
            cS.s_1yr_transitionV = zeros(cS.aD_new, 1);
            for a = 1:(cS.aD_new - 1) % 遍历到倒数第二个年龄组
                lastYearIdxInGroup = cS.physAgeMap{a}(end); % 获取当前年龄组最后一个年度年龄的索引
                if lastYearIdxInGroup < cS.aD_orig % 如果不是最后一个年度年龄
                    cS.s_1yr_transitionV(a) = cS.s_orig(lastYearIdxInGroup); % 使用该年度年龄的存活率
                else
                    cS.s_1yr_transitionV(a) = 0; % 如果是最后一个年度年龄，则无法转移到下一组
                end
            end
            cS.s_1yr_transitionV(cS.aD_new) = 0; % 最后一个年龄组的转移存活率为0

            % 人口动态模拟参数
            % 初始人口分布 (按年龄组, 数据来自某个基准校准, 总和约为1000)
            cS.initial_pop = [76.2, 86.4, 113.8, 98.6, 86.6, 102.7, 112.0, 99.0, 64.0, 66.9, 44.1, 25.4, 14.9, 6.8, 1.7, 0.2];
            if length(cS.initial_pop) ~= cS.aD_new % 如果提供的初始人口数据长度与年龄组数不匹配
                cS.initial_pop = ones(cS.aD_new,1) * (100/cS.aD_new); % 使用均匀分布
                warning('initial_pop长度与年龄组数不匹配，已重设为均匀分布。');
            end

            % 年龄组间的存活概率 (用于人口动态模拟，每 cS.yearStep 年的存活率)
            % beta_surv_pop(j) 表示从第 j 组存活到第 j+1 组的概率
            beta_surv_pop = [0.998,0.996,0.994,0.992,0.988,0.984,0.980,0.976,0.970,0.960,0.945,0.920,0.880,0.800,0.680];
            if length(beta_surv_pop) ~= cS.aD_new - 1 % 检查长度是否正确
                error('年龄组间存活率 beta_surv_pop 的长度对于 %d 个年龄组不正确。应为 %d。', cS.aD_new, cS.aD_new -1);
            end
            cS.survivalProbV_popdyn = [beta_surv_pop(:)', 0]; % 最后一个年龄组到下一个不存在的组的存活率为0

            cS.bgp_tolerance = 0.001; % 人口动态稳态判断的容忍度 (年龄结构变化范数)
            cS.bgp_window = 5;       % 人口动态稳态判断的窗口期 (期数)
            cS.max_periods = 50;     % 人口动态模拟的最大期数

            % 家庭参数 (年度化)
            cS.sigma      = 1.5;    % CRRA效用函数的风险厌恶系数
            cS.beta       = 1.011;  % 主观贴现因子 (年度)
            cS.cFloor     = 0.05;   % 最低消费量 (避免效用函数出问题)
            cS.nSim       = 1000;   % 微观模拟中的个体数量 (V7修改：减少以加快测试)

            % 技术参数 (年度化)
            cS.A          = 0.895944; % TFP 总产出因子
            cS.alpha      = 0.36;   % 资本产出弹性
            cS.ddk        = 0.06;   % 资本折旧率 (年度)

            % 社会保障参数 (年度化) - 外生PAYG税率theta (仅作为v4的基准设定，v7中内生)
            cS.theta      = 0.20;   % PAYG 工资税率 (外生基准值)

            % Heer (2020) 风格的政府财政参数 (年度化)
            cS.tau_k = 0.20;          % 资本所得税率 (对 r_gross_market - ddk 征收)
            cS.tau_c = 0.10;          % 消费税率
            cS.tau_l = 0.05;          % 额外的劳动所得税率 (在PAYG税率cS.theta之外)
            cS.gov_exp_frac_Y = 0.15; % 政府消费占GDP的比例
            cS.gov_debt_frac_Y = 0.60;% 政府债务占GDP的比例

            % 劳动禀赋参数 (基于年度模型)
            cS.leSigma1 = 0.38^0.5;     % 劳动效率初始分布的标准差 (对数值)
            cS.leShockStd = 0.045^0.5;  % 劳动效率冲击的标准差 (对数值)
            cS.lePersistence = 0.96;    % 劳动效率冲击的持续性 (AR1系数)
            cS.leWidth = 4;             % Tauchen方法离散化时的宽度 (标准差倍数)
            cS.nw = 5;                  % 劳动效率冲击状态的数量 (V7修改：减少以加快测试)

            % 资产网格参数 (kGridV 代表非PPS资产的网格)
            cS.tgKY = 3; % 目标资本产出比 (用于设定网格大致范围)
            cS.tgWage = (1-cS.alpha)*cS.A*((cS.tgKY/cS.A)^(cS.alpha/(1-cS.alpha))); % 基于目标K/Y的估算工资
            cS.nk = 30;                % 非PPS资产网格点数量
            cS.kMin = 0;               % 非PPS资产网格下限
            cS.kMax = 100 * cS.tgWage; % 非PPS资产网格上限 (基于估算工资的倍数)
            power = 1.5; % 网格点分布的幂指数 (>1表示点在下限附近更密集)
            kGridV = cS.kMin + (cS.kMax - cS.kMin) * (linspace(0, 1, cS.nk).^power);
            if cS.nk > 0
                kGridV(1)=cS.kMin; % 确保第一个点是下限
            end
            cS.kGridV = kGridV(:); % 存储为列向量

            % 年龄效率剖面 (年度数据)
            ageEffV_orig_temp = zeros(100, 1); % 临时向量，覆盖较宽年龄范围
            % 分段线性设定年龄效率 (20-36岁上升, 37-47岁平稳, 48-65岁下降, 66-72岁快速下降至0)
            ageEffV_orig_temp(20:72) = [linspace(0.3,1.5,36-20+1), 1.5.*ones(1,47-37+1), linspace(1.5,0.2,65-48+1), linspace(0.18,0,72-66+1)];
            cS.ageEffV_orig = ageEffV_orig_temp(cS.age1_orig : cS.ageLast_orig); % 截取模型相关的年度年龄段
            if length(cS.ageEffV_orig) ~= cS.aD_orig
                error('ageEffV_orig 年度年龄效率剖面长度不匹配');
            end
            % 将年度年龄效率转换为模型年龄组的平均效率
            cS.ageEffV_new = zeros(cS.aD_new, 1);
            for a = 1:cS.aD_new
                cS.ageEffV_new(a) = mean(cS.ageEffV_orig(cS.physAgeMap{a}));
            end

            % 个人养老金计划 (PPS) 参数 (v4 基准设定)
            cS.pps_active = true; % PPS计划是否激活
            cS.pps_max_contrib_frac = 0.1;    % PPS缴费占工资收入的比例上限 (法定，v7中可能被规则覆盖但仍作为法定上限)
            cS.pps_tax_rate_withdrawal = 0.03;% PPS提取时的税率
            cS.pps_return_rate_premium = 0.01; % PPS资产相对于市场净回报率的额外回报率
            cS.pps_withdrawal_rate = 0.15;    % 退休后每年从PPS余额中提取的比例
            cS.pps_contribution_age_max = cS.aR_idx_orig - 1; % 允许PPS缴费的最大年度年龄索引 (退休前一年)
            cS.pps_withdrawal_age_min = cS.aR_idx_orig;   % 允许PPS提取的最小年度年龄索引 (退休当年)
            cS.pps_in_K = true;               % PPS资产是否计入生产性资本K
            cS.pps_bequeathable = false;      % PPS资产是否可遗赠

            % --- V7 特有参数或对V4基准的修改 ---
            cS.tau_l_init_guess = 0.05; % 一般劳动所得税率 tau_l 的初始猜测值 (用于迭代求解)
            cS.tau_l_min = 0.00;        % tau_l 的下限
            cS.tau_l_max = 0.3;         % tau_l 的上限 (调整以适应模型)
            cS.max_total_labor_tax = 0.6; % 总劳动税负上限 (theta_payg + tau_l)

            cS.theta_payg_max = 0.35; % PAYG 税率的上限

            % PPS 参数 (V7 特化设定)
            % cS.pps_active 沿用v4设定 (true)
            cS.pps_annual_contrib_limit = 1.5; % PPS年度绝对缴费上限 (例如，金额单位)
            cS.pps_max_contrib_frac = 0.15;    % PPS缴费占工资收入的比例上限 (法定，高于v4)
            % cS.pps_in_K 沿用v4设定 (true)
            % cS.pps_bequeathable 沿用v4设定 (false)
            % cS.pps_tax_rate_withdrawal 沿用v4设定 (0.03)
            cS.pps_return_rate_premium = 0.01; % PPS资产回报率溢价 (略微提高)
            % cS.pps_withdrawal_rate 沿用v4设定 (0.15)
            % 年度年龄索引 (基于原始年度年龄)
            cS.pps_contribution_age_max_idx = cS.aR_idx_orig - 1; % 最大缴费年龄索引 (同v4: 退休前一年)
            cS.pps_withdrawal_age_min_idx = cS.aR_idx_orig;       % 最低提取年龄索引 (同v4: 退休当年)

            % PPS资产网格 (k_pps)
            cS.nkpps = 20;             % PPS资产网格点数量
            cS.kppsMin = 0;            % PPS资产网格下限
            cS.kppsMax = cS.kMax / 2;  % PPS资产网格上限 (示例：为非PPS资产上限的一半)
            if cS.nkpps > 0 % 确保kppsMax至少是一个小的正数（如果nkpps>0）
                cS.kppsMax = max(cS.kppsMax, 1e-3);
            end
            power_kpps = 1.5; % PPS资产网格分布幂指数
            if cS.nkpps > 1
                kppsGridV_temp = cS.kppsMin + (cS.kppsMax - cS.kppsMin) * (linspace(0, 1, cS.nkpps).^power_kpps);
                kppsGridV_temp(1) = cS.kppsMin; % 确保第一个点是下限
            elseif cS.nkpps == 1
                kppsGridV_temp = cS.kppsMin; % 如果只有一个点，则为下限 (或可设为平均值)
            else % nkpps = 0 (网格为空)
                kppsGridV_temp = [];
            end
            cS.kppsGridV = kppsGridV_temp(:); % 存储为列向量

            % 均衡求解器参数
            cS.max_iter_K_tau_l = 100; % K 和 tau_l 迭代求解的最大次数
            cS.tol_K_tau_l = 1e-4;     % K 和 tau_l 迭代求解的收敛容忍度 (范数)
            cS.damp_K_v5 = 0.3;        % K 迭代的阻尼系数
            cS.damp_tau_l_v5 = 0.3;    % tau_l 迭代的阻尼系数
            cS.gbc_tol_for_internal_loop = 1e-3; % 内层循环中政府预算约束的容忍度
            cS.gbc_tol_for_rho_search = 1e-2;    % 外层rho_prime搜索中GBC的容忍度 (未使用，因内层已严格平衡)

            cS.max_stagnation_iters = 10;        % 最大停滞迭代次数 (如果范数改善过小)
            cS.min_norm_improvement_frac = 1e-3; % 范数改善的最小比例 (低于此则认为停滞)
            cS.max_tau_l_boundary_strikes = 5;   % tau_l 触达边界且GBC未平衡的最大次数

            % V7 特有：PPS规则化缴费计划 (按年龄组占当期工资收入的比例)
            cS.pps_fixed_contrib_schedule_frac = zeros(cS.aD_new, 1); % 初始化为列向量
            
            num_working_age_groups = cS.aR_new; % 获取工作年龄组的数量
            if num_working_age_groups > 0 % 如果存在工作年龄组
                % 示例：分段线性设定缴费比例 (年轻时较低，中年较高，临近退休略降)
                if num_working_age_groups == 1 % 如果只有一个工作年龄组
                     cS.pps_fixed_contrib_schedule_frac(1:num_working_age_groups) = 0.05; % 示例值
                elseif num_working_age_groups > 1 % 如果有多个工作年龄组
                    mid_point1 = ceil(num_working_age_groups / 3); % 第一个三分点
                    mid_point2 = ceil(2 * num_working_age_groups / 3); % 第二个三分点
                    
                    if mid_point1 > 0 % 第一段 (早期工作)
                        cS.pps_fixed_contrib_schedule_frac(1:mid_point1) = linspace(0.02, 0.06, mid_point1);
                    end
                    if mid_point2 > mid_point1 % 第二段 (中期工作)
                        cS.pps_fixed_contrib_schedule_frac(mid_point1+1:mid_point2) = linspace(0.06, 0.10, mid_point2 - mid_point1);
                    end
                    if num_working_age_groups > mid_point2 % 第三段 (后期工作)
                        cS.pps_fixed_contrib_schedule_frac(mid_point2+1:num_working_age_groups) = linspace(0.10, 0.04, num_working_age_groups - mid_point2);
                    end
                end
            end
            % 确保退休年龄组的规则化缴费为0
            if cS.aR_new < cS.aD_new % 如果存在退休年龄组
                cS.pps_fixed_contrib_schedule_frac(cS.aR_new + 1 : cS.aD_new) = 0;
            end
            
            % fminbnd 优化器参数 (用于 HHSolutionByOneState_OptK_Mod)
            cS.fminbnd_TolX = 1e-6;      % fminbnd 的 TolX 参数
            cS.fminbnd_Display = 'none'; % fminbnd 的 Display 参数

            fprintf('V7: 完整参数已设置完毕。\n');
            fprintf('V7 PPS参数: 固定缴费比例计划(按年龄组)已设定。例如，第1组: %.3f, 第%d组: %.3f\n', ...
                cS.pps_fixed_contrib_schedule_frac(1), cS.aR_new, cS.pps_fixed_contrib_schedule_frac(cS.aR_new));
        end

        % =====================================================================
        % == 人口动态函数 (年龄组层面) ==
        % =====================================================================

        function popS = initPopulation(cS)
            % initPopulation - 初始化人口结构
            %
            % 输入:
            %   cS   - 参数结构体
            % 输出:
            %   popS - 包含初始人口信息和分布的结构体
            %
            % 说明:
            %   根据cS.initial_pop设置初始年龄组人口。
            %   如果cS.initial_pop无效或总和为零，则设置为均匀分布。
            %   计算初始总人口和年龄分布。

            popS.Z = zeros(cS.aD_new, 1); % 初始化年龄组人口数量向量 Z (列向量)
            initial_total = sum(cS.initial_pop); % 计算提供的初始人口总和

            if initial_total > 0 && length(cS.initial_pop) == cS.aD_new
                % 如果初始人口数据有效且长度匹配
                popS.Z(:, 1) = cS.initial_pop(:) / initial_total * 100; % 归一化为百分比，并存入第一列
            else
                % 如果数据无效或不匹配
                warning('初始人口数据不匹配或总和为零。将设置为均匀的初始年龄组人口分布。');
                popS.Z(:, 1) = 100 / cS.aD_new; % 设置为均匀分布
            end

            popS.totalPop = sum(popS.Z(:, 1)); % 计算总人口 (此时应为100)

            if popS.totalPop > 1e-9 % 避免除以零
                popS.ageDist = popS.Z(:, 1) / popS.totalPop; % 计算年龄分布 (应为归一化比例)
            else
                popS.ageDist = zeros(cS.aD_new, 1); % 如果总人口为零，则年龄分布也为零
            end
            popS.initialAgeDist = popS.ageDist; % 存储初始年龄分布

            fprintf('初始年龄组人口已设置。总人口=%.2f (代表百分比基数)。\n', popS.totalPop);
        end

        function popS = populationDynamics(popS, cS)
            % populationDynamics - 模拟人口随时间的动态演变
            %
            % 输入:
            %   popS - 包含初始人口信息的结构体
            %   cS   - 参数结构体
            % 输出:
            %   popS - 更新后的结构体，包含人口历史、总人口历史、年龄分布历史和抚养比历史
            %
            % 说明:
            %   模拟从初始人口分布开始，根据生育率(隐含在第一期人口增长中)和
            %   年龄组间的存活率(cS.survivalProbV_popdyn)，迭代计算未来各期的人口结构。
            %   模拟直到达到稳态(BGP)或达到最大模拟期数(cS.max_periods)。
            %   生育率部分采用随时间变化的设定。

            max_periods_sim = cS.max_periods; % 获取最大模拟期数
            % 初始化历史数据存储矩阵
            Z_history = zeros(cS.aD_new, max_periods_sim + 1);         % 年龄组人口数量历史
            totalPop_history = zeros(1, max_periods_sim + 1);         % 总人口历史
            ageDist_history = zeros(cS.aD_new, max_periods_sim + 1);  % 年龄分布历史

            % 存储初始期 (t=0, 对应历史数据索引1) 的数据
            Z_history(:, 1) = popS.Z(:, 1);
            totalPop_history(1) = popS.totalPop(1);
            ageDist_history(:, 1) = popS.ageDist(:, 1);

            fprintf('人口动态模拟开始 (年龄组, 最大期数 = %d)...\n', max_periods_sim);
            bgp_reached_flag = false; % 稳态是否达到的标志
            actual_periods_run = max_periods_sim; % 实际运行的期数

            for t = 1:max_periods_sim % 迭代模拟每一期
                if mod(t, 10) == 0 || t == 1 % 每10期或第一期打印进度
                    fprintf('  模拟人口期数 %d (年龄组)\n', t);
                end

                Z_current_period = Z_history(:, t); % 获取当前期的人口分布
                Z_next_period = zeros(cS.aD_new, 1); % 初始化下一期的人口分布

                % 计算下一期第一个年龄组的人口 (新生儿/年轻人口)
                % 这里使用一个随时间变化的人口增长率 (或可视为生育率变化的影响)
                time_varying_growth_rate = 0;
                if t < 6 % 前5期 (t=1 to 5)
                    time_varying_growth_rate = -0.01 - 0.003 * t; % 示例：增长率下降
                else % 后续期
                    time_varying_growth_rate = -0.03 - 0.004 * min(t-6, 10); % 示例：增长率进一步下降并趋于稳定
                end
                Z_next_period(1) = Z_current_period(1) * (1 + time_varying_growth_rate); % 基于上一期同龄组和增长率
                Z_next_period(1) = max(0, Z_next_period(1)); % 确保非负

                % 计算其他年龄组的人口 (基于上一期年轻一组的存活)
                for a = 2:cS.aD_new % 从第二个年龄组开始
                    survival_prob_group = 0; % 初始化年龄组间存活率
                    % cS.survivalProbV_popdyn(j) 是从第 j 组到第 j+1 组的存活率
                    if (a-1) >= 1 && (a-1) <= length(cS.survivalProbV_popdyn)
                        survival_prob_group = cS.survivalProbV_popdyn(a-1); % 获取从 a-1 组到 a 组的存活率
                    end
                    Z_next_period(a) = Z_current_period(a-1) * survival_prob_group; % 上一期 a-1 组人口 * 存活率
                    Z_next_period(a) = max(0, Z_next_period(a)); % 确保非负
                end

                % 存储下一期的人口数据 (在历史数据索引 t+1 处)
                Z_history(:, t+1) = Z_next_period;
                totalPop_history(t+1) = sum(Z_next_period);
                if totalPop_history(t+1) > 1e-9 % 避免除以零
                    ageDist_history(:, t+1) = Z_next_period / totalPop_history(t+1);
                else % 如果总人口趋近于零
                    ageDist_history(:, t+1) = 0; % 年龄分布设为零
                    totalPop_history(t+1) = 0;  % 总人口设为零
                end

                % 检查是否达到稳态 (BGP)
                current_check_period_idx = t + 1; % 当前检查的是历史数据中的这一点
                if current_check_period_idx >= cS.bgp_window + 1 % 确保有足够的数据进行窗口比较
                    stable = true; % 假设已稳定
                    for w_idx = 1:cS.bgp_window % 遍历窗口内的每两期间隔
                       hist_idx1 = current_check_period_idx - w_idx + 1; % 窗口中较新的一期
                       hist_idx2 = current_check_period_idx - w_idx;     % 窗口中较老的一期
                       if hist_idx1 > 0 && hist_idx2 > 0 && hist_idx1 <= size(ageDist_history,2) && hist_idx2 <= size(ageDist_history,2)
                           % 计算年龄分布向量的欧几里得范数差异
                           change = norm(ageDist_history(:, hist_idx1) - ageDist_history(:, hist_idx2));
                           if change >= cS.bgp_tolerance % 如果变化大于容忍度
                               stable = false; % 则未稳定
                               break; % 跳出窗口检查
                           end
                       else % 如果索引无效 (不太可能发生在此逻辑下，但作为保险)
                           stable = false; break;
                       end
                    end
                    if stable % 如果在整个窗口内都稳定
                        fprintf('\n人口稳态 (年龄组) 在模拟期数 %d (对应历史数据索引 %d) 达到。\n', t, current_check_period_idx);
                        bgp_reached_flag = true; % 标记已达到稳态
                        actual_periods_run = t; % 记录实际运行的模拟期数 (不含初始期)
                        break; % 跳出主模拟循环
                    end
                end
            end % 结束人口动态模拟循环

            % 截取实际运行期数的数据
            final_period_idx_to_store = min(actual_periods_run + 1, size(Z_history,2)); % +1 是因为包含初始期
            popS.Z = Z_history(:, 1:final_period_idx_to_store);
            popS.totalPop = totalPop_history(1:final_period_idx_to_store);
            popS.ageDist = ageDist_history(:, 1:final_period_idx_to_store);

            % 计算抚养比历史 (基于模拟的每一期, 不含初始期)
            depRatio_history = zeros(1, actual_periods_run);
            for th_loop = 1:actual_periods_run % 遍历模拟的每一期 (t=1 to actual_periods_run)
                 Z_t_for_depratio = Z_history(:, th_loop + 1); % 获取该期的人口分布 (th_loop+1 对应 Z_history 的索引)
                 working_pop = sum(Z_t_for_depratio(1:cS.aR_new)); % 工作年龄人口 (模型年龄组1到aR_new)
                 retired_pop = sum(Z_t_for_depratio(cS.aR_new+1:end)); % 退休年龄人口
                 if working_pop > 1e-9 % 避免除以零
                     depRatio_history(th_loop) = retired_pop / working_pop;
                 else
                     depRatio_history(th_loop) = inf; % 如果工作人口为零，抚养比为无穷大
                 end
            end
            popS.dependencyRatio = depRatio_history; % 存储抚养比历史

            fprintf('人口动态模拟完成。运行期数: %d。达到BGP: %d\n', actual_periods_run, bgp_reached_flag);
            if ~bgp_reached_flag
                fprintf('警告: 人口稳态 (年龄组) 未在 %d 期内达到。\n', max_periods_sim);
            end
        end

        function [Z_ss, dependency_ratio_ss, bgp_reached, bgp_period] = detectSteadyStatePopulation(popS, cS)
            % detectSteadyStatePopulation - 从已模拟的人口历史数据中检测稳态
            %
            % 输入:
            %   popS - 包含人口历史数据的结构体
            %   cS   - 参数结构体
            % 输出:
            %   Z_ss                - 稳态(或最终期)的年龄组人口数量 (未归一化)
            %   dependency_ratio_ss - 稳态(或最终期)的抚养比
            %   bgp_reached         - 是否检测到稳态的标志 (true/false)
            %   bgp_period          - 检测到稳态的模拟期数 (从0开始计数, 0代表初始期)
            %                         如果未达到稳态，则为最后一期
            %
            % 说明:
            %   此函数与 populationDynamics 中的稳态检测逻辑类似，但作用于已有的 popS.ageDist 历史数据。
            %   它从后向前检查，以找到最近的稳定窗口。
            %   同时绘制初始与稳态(或最终)年龄组人口分布图。

            actual_periods_in_data = size(popS.Z, 2); % 获取数据中的总期数 (包含初始期，所以数据点数比模拟期数多1)
            bgp_reached = false; % 初始化稳态标志
            bgp_period = actual_periods_in_data - 1; % 默认稳态期为最后一期 (模拟期数，从0开始)

            if actual_periods_in_data < cS.bgp_window + 1 % 如果数据点数不足以进行窗口检查
                fprintf('人口模拟期数过短 (%d 数据点)，无法进行稳态检查 (窗口期 = %d)。\n', actual_periods_in_data, cS.bgp_window);
            else
                fprintf('检查人口稳态 (年龄组, 最近 %d 期)...\n', cS.bgp_window);
                % 从数据的最后一期向前检查，t_check_end 是窗口的结束点在 popS.ageDist 中的索引
                for t_check_end_idx = actual_periods_in_data : -1 : cS.bgp_window + 1
                    stable = true; % 假设当前窗口是稳定的
                    for w_idx = 0 : (cS.bgp_window - 1) % 遍历窗口内的比较对
                        idx1 = t_check_end_idx - w_idx;       % 窗口中较新的一期
                        idx2 = t_check_end_idx - w_idx - 1;   % 窗口中较老的一期
                        if idx1 > 0 && idx2 > 0 && idx1 <= size(popS.ageDist,2) && idx2 <= size(popS.ageDist,2)
                            change = norm(popS.ageDist(:, idx1) - popS.ageDist(:, idx2));
                            if change >= cS.bgp_tolerance % 如果变化大于容忍度
                                stable = false; % 则不稳定
                                break; % 跳出窗口内比较
                            end
                        else % 索引无效 (不太可能)
                            stable = false; break;
                        end
                    end
                    if stable % 如果整个窗口都稳定
                        bgp_reached = true;
                        % t_check_end_idx 是稳定窗口的结束点 (数据索引)
                        % 稳态开始的模拟期数是这个窗口之前的那个点，即 idx2-1 (从0开始计数的模拟期)
                        % 或者说，稳定是从 t_check_end_idx - cS.bgp_window 这个数据点开始的
                        % 我们通常用稳定窗口的最后一点的前一点作为稳态的代表期
                        bgp_period = t_check_end_idx - 1 -1; % -1转为模拟期数, 再-1表示稳定窗口前的点
                        bgp_period = t_check_end_idx - cS.bgp_window ; % 稳态开始的期数 (0-indexed)
                        bgp_period = t_check_end_idx - 1; % 用窗口的最后一点作为稳态代表点 (数据索引转模拟期数)
                                                          % 这里把 bgp_period 理解为稳定状态所代表的模拟期数t
                                                          % 如果 t_check_end_idx 是数据索引, 则对应的模拟期数是 t_check_end_idx - 1
                        fprintf('人口稳态 (年龄组) 从模拟期数 %d (数据索引 %d) 开始检测到 (稳定窗口结束于此)。\n', bgp_period, t_check_end_idx);
                        break; % 找到最近的稳定窗口，跳出检查循环
                    end
                end
                if ~bgp_reached % 如果循环结束仍未检测到稳态
                     fprintf('未检测到人口稳态 (年龄组)。将使用最终期数据。\n');
                     bgp_period = actual_periods_in_data - 1; % 确保 bgp_period 是最后一期 (0-indexed sim period)
                end
            end

            % 获取稳态数据
            % ss_data_index 是 popS.Z 等历史数据中的列索引 (1-indexed)
            ss_data_index = min(bgp_period + 1, size(popS.Z, 2)); % +1 将0-indexed模拟期转为1-indexed数据列
            Z_ss = popS.Z(:, ss_data_index); % 未归一化的稳态人口数量

            % 计算归一化的稳态年龄分布 (用于绘图和部分计算)
            Z_ss_norm = zeros(cS.aD_new, 1);
            if sum(Z_ss) > 1e-9
                Z_ss_norm = Z_ss / sum(Z_ss);
            end

            % 获取稳态抚养比
            % valid_dep_ratio_index 对应 popS.dependencyRatio 中的索引 (1-indexed)
            % popS.dependencyRatio 存储的是模拟期数 t=1 到 t=actual_periods_run 的抚养比
            % 所以索引 i 对应模拟期数 i。
            % bgp_period 是0-indexed的模拟期数。
            valid_dep_ratio_index = min(bgp_period, length(popS.dependencyRatio)); % bgp_period如果是0, 则取第0个抚养比，不合理
            if bgp_period == 0 % 如果稳态是初始期 (不太可能，但作为边界情况)
                valid_dep_ratio_index = 1; % 尝试取第一个模拟期(t=1)的抚养比
            end
            
            if isfield(popS, 'dependencyRatio') && ~isempty(popS.dependencyRatio) && ...
               valid_dep_ratio_index > 0 && valid_dep_ratio_index <= length(popS.dependencyRatio)
                dependency_ratio_ss = popS.dependencyRatio(valid_dep_ratio_index);
            else % 如果抚养比历史数据有问题，则重新计算
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

            % 绘制初始与稳态/最终年龄组人口分布图
            figure('Name', 'V7: 初始 vs 稳态/最终 年龄组人口分布');
            hold on;
            group_indices = 1:cS.aD_new; % 年龄组的X轴坐标
            if isfield(popS, 'initialAgeDist') && ~isempty(popS.initialAgeDist)
                bar(group_indices - 0.2, popS.initialAgeDist * 100, 0.4, 'b', 'DisplayName', '初始年龄组分布');
            else
                warning('未找到或空的初始年龄分布用于绘图。');
            end
            bar(group_indices + 0.2, Z_ss_norm * 100, 0.4, 'r', 'DisplayName', sprintf('稳态年龄组分布 (模拟期 %d)', bgp_period));
            hold off;
            xlabel(sprintf('年龄组索引 (1 至 %d)', cS.aD_new));
            ylabel('占总人口百分比 (%)');
            title(sprintf('V7: 初始 vs 稳态/最终 年龄组人口分布 (稳态代表模拟期 t=%d)', bgp_period));
            legend('Location', 'best');
            xticks(group_indices);
            grid on;
            drawnow; % 立即显示图像
            fprintf('已绘制初始与稳态/最终年龄组人口分布图。\n');
        end


        % =====================================================================
        % == 劳动禀赋过程函数 ==
        % =====================================================================

        function [logGridV, trProbM, prob1V] = EarningProcess_olgm(cS)
            % EarningProcess_olgm - 生成劳动效率冲击的网格、转移矩阵和初始分布
            %
            % 输入:
            %   cS - 参数结构体，需要包含:
            %        nw, lePersistence, leShockStd, leWidth, leSigma1
            % 输出:
            %   logGridV - 劳动效率对数值的网格点 (列向量)
            %   trProbM  - 劳动效率状态间的转移概率矩阵 (P(i,j) = Pr(s' = j | s = i))
            %   prob1V   - 劳动效率初始状态的概率分布 (列向量)
            %
            % 说明:
            %   使用Tauchen方法生成离散化的AR(1)过程网格和转移矩阵。
            %   使用norm_grid生成初始分布。
            %   对数效率网格点进行了标准化处理。

            % 1. 使用Tauchen方法生成原始的对数效率网格和转移矩阵
            %    pMu (均值) 设为0，因为通常效率过程围绕某个均值，这里可以后续调整或通过ageEffV体现
            [logGridV_raw, trProbM_calc] = main_olg_v7_utils.tauchen(cS.nw, cS.lePersistence, cS.leShockStd, 0, cS.leWidth);

            % 2. 生成初始分布 (基于norm_grid)
            gridMin_raw = logGridV_raw(1); % 原始网格最小值
            gridMax_raw = logGridV_raw(end); % 原始网格最大值
            % 为norm_grid设定积分边界 (略微扩展以覆盖尾部)
            normGridMin = gridMin_raw - 2*cS.leShockStd;
            normGridMax = gridMax_raw + 2*cS.leShockStd;
            prob1V_calc = []; % 初始化初始分布
            try
                % cS.leSigma1 是初始对数效率分布的标准差
                [prob1V_calc, ~, ~] = main_olg_v7_utils.norm_grid(logGridV_raw, normGridMin, normGridMax, 0, cS.leSigma1);
                prob1V_calc = prob1V_calc(:); % 确保是列向量
            catch ME
                warning('norm_grid 计算初始分布失败: %s。将使用均匀分布作为备用。', ME.message);
                prob1V_calc = ones(cS.nw, 1)/cS.nw; % 备用：均匀分布
            end

            % 3. 标准化对数效率网格
            logGridV_calc = logGridV_raw(:); % 确保是列向量
            % 标准化：使得第一个点为-1 (或某个固定负值)，其他点相对平移
            % 这种标准化方式比较特殊，常见的是减去均值或使均值为0。这里的目的是什么？
            % 原始代码：logGridV_calc = logGridV_calc - logGridV_calc(1) - 1;
            % 考虑更常见的标准化：使其均值为0
            % logGridV_calc = logGridV_calc - mean(logGridV_calc);
            % 或者，如果希望exp(logGridV)均值为1，则 logGridV_calc = logGridV_calc - mean(logGridV_calc) - log(mean(exp(logGridV_calc - mean(logGridV_calc))))
            % 暂时保留原始代码的标准化方式，但需注意其含义。
            logGridV_calc = logGridV_calc - logGridV_calc(1) - 1; % 原始标准化

            % 4. 确保转移矩阵每行加总为1 (处理可能的数值误差)
            if any(abs(sum(trProbM_calc, 2) - 1) > 1e-6) % 如果行和显著不为1
                 row_sums = sum(trProbM_calc, 2);
                 row_sums(row_sums <= 1e-9) = 1; % 避免除以零
                 trProbM_calc = bsxfun(@rdivide, trProbM_calc, row_sums); % 按行归一化
            end

            % 5. 确保初始分布加总为1
            if abs(sum(prob1V_calc) - 1) > 1e-6 % 如果总和显著不为1
                prob1V_calc = prob1V_calc ./ sum(prob1V_calc); % 归一化
            end

            % 6. 赋值输出
            logGridV = logGridV_calc;
            trProbM = trProbM_calc;
            prob1V = prob1V_calc;
            fprintf('劳动禀赋过程参数已生成 (Tauchen & norm_grid)。\n');
        end

        function [y_grid_out, trProbM_out] = tauchen(N_states, persistence_rho, shock_sigma, mean_val_mu, num_std_dev_width)
            % tauchen - Tauchen (1986) 方法离散化 AR(1) 过程
            %
            % y_t = (1-rho)*mu + rho * y_{t-1} + epsilon_t,  epsilon_t ~ N(0, shock_sigma^2)
            %
            % 输入:
            %   N_states          - 离散状态的数量
            %   persistence_rho   - AR(1) 过程的持续性系数 (rho)
            %   shock_sigma       - 冲击的标准差 (sigma_epsilon)
            %   mean_val_mu       - AR(1) 过程的条件均值 (mu, 即 y 无条件均值若平稳)
            %                       (注意: 这里的mu是公式中的mu, 不是过程本身的无条件均值, 除非rho=0)
            %   num_std_dev_width - 网格宽度，以过程的无条件标准差的倍数表示
            %
            % 输出:
            %   y_grid_out        - 离散化的状态值网格 (1xN_states行向量)
            %   trProbM_out       - 状态间的转移概率矩阵 (N_states x N_states)
            %
            % 说明:
            %   计算AR(1)过程的无条件标准差 std_y。
            %   在 [-num_std_dev_width*std_y, num_std_dev_width*std_y] 区间内生成等距网格点。
            %   计算从每个当前状态到下一个所有状态的转移概率。
            %   最后将网格点加上过程的无条件均值 (如果 mu 和 rho 指定)。

            % 计算过程的无条件标准差 (如果过程是平稳的, |rho|<1)
            std_y_unconditional = sqrt(shock_sigma^2 / (1 - persistence_rho^2));
            if abs(1-persistence_rho) < 1e-9 % 如果 rho 接近 1 (单位根或爆炸过程)
                std_y_unconditional = shock_sigma * 100; % 使用一个较大的值作为近似
                warning('Tauchen: persistence_rho (%.4f) 接近1，std_y_unconditional可能不准确。', persistence_rho);
            end

            % 确定网格的上下界 (对称于0)
            y_max_boundary = num_std_dev_width * std_y_unconditional;
            y_min_boundary = -y_max_boundary;

            % 生成等距的网格点 (初始时中心为0)
            y_grid_centered = linspace(y_min_boundary, y_max_boundary, N_states);

            % 计算网格点之间的步长
            step_size_d = 0;
            if N_states > 1
                step_size_d = y_grid_centered(2) - y_grid_centered(1);
            end

            % 初始化转移概率矩阵
            trProbM_calc = zeros(N_states, N_states);

            % 计算转移概率
            for iRow = 1:N_states % 当前状态 y_grid_centered(iR)
                for iCol = 1:N_states % 下一状态 y_grid_centered(iC)
                    % 下一期 y' 的条件均值，给定当前 y = y_grid_centered(iR)
                    % 假设 AR(1) 过程为: y' = rho * y + epsilon
                    % 如果原始公式是 y' = (1-rho)*mu_uncond + rho*y + epsilon, 则这里 mean_next_y 应包含 (1-rho)*mu_uncond
                    % 但通常Tauchen先处理以0为中心的过程，最后平移。
                    % 这里 shock_sigma 是 epsilon 的标准差。
                    % mean_val_mu 在这里暂不使用，将在最后平移网格。
                    mean_next_y_conditional = persistence_rho * y_grid_centered(iRow);

                    % 对于每个目标区间 [y_grid_centered(iC) - d/2, y_grid_centered(iC) + d/2]
                    % 计算 Pr(y'落入此区间 | y_current)
                    if iCol == 1 % 第一个状态的区间是 (-inf, y_grid_centered(1) + d/2]
                        trProbM_calc(iRow,iCol) = normcdf((y_grid_centered(1) - mean_next_y_conditional + step_size_d/2) / shock_sigma);
                    elseif iCol == N_states % 最后一个状态的区间是 [y_grid_centered(N) - d/2, +inf)
                        trProbM_calc(iRow,iCol) = 1 - normcdf((y_grid_centered(N_states) - mean_next_y_conditional - step_size_d/2) / shock_sigma);
                    else % 中间状态的区间是 [y_grid_centered(iC) - d/2, y_grid_centered(iC) + d/2]
                        upper_bound_cdf = normcdf((y_grid_centered(iCol) - mean_next_y_conditional + step_size_d/2) / shock_sigma);
                        lower_bound_cdf = normcdf((y_grid_centered(iCol) - mean_next_y_conditional - step_size_d/2) / shock_sigma);
                        trProbM_calc(iRow,iCol) = upper_bound_cdf - lower_bound_cdf;
                    end
                end
            end

            % 确保转移矩阵每行加总为1 (处理数值误差)
            row_sums_check = sum(trProbM_calc,2);
            row_sums_check(row_sums_check <= 1e-9) = 1; % 避免除以零
            trProbM_out = bsxfun(@rdivide, trProbM_calc, row_sums_check);

            % 计算过程的无条件均值 (如果指定了mean_val_mu且过程平稳)
            % 无条件均值 E[y] = mu (如果过程是 y' = (1-rho)mu + rho*y + eps)
            unconditional_mean_shift = mean_val_mu / (1-persistence_rho);
            if abs(1-persistence_rho) < 1e-9 && mean_val_mu ~= 0 % rho=1且mu非零，均值是无穷或取决于初始条件
                unconditional_mean_shift = sign(mean_val_mu) * inf; % 或其他处理
                warning('Tauchen: rho=1 and mu non-zero, unconditional mean is ill-defined.');
            elseif abs(1-persistence_rho) < 1e-9 && mean_val_mu == 0 % rho=1且mu=0 (随机游走)
                unconditional_mean_shift = 0; % 网格保持中心为0 (或由y0决定)
            end

            % 将网格点平移，使其均值为 unconditional_mean_shift
            y_grid_out = y_grid_centered; % 默认为中心化的网格
            if isfinite(unconditional_mean_shift)
                y_grid_out = y_grid_centered + unconditional_mean_shift;
            end
        end

        function [massV_out, lbV_out, ubV_out] = norm_grid(x_grid_points, overall_min_bound, overall_max_bound, dist_mean_mu, dist_std_sig)
            % norm_grid - 计算给定网格点在正态分布下的概率质量
            %
            % 输入:
            %   x_grid_points     - 网格点向量 (1xN 或 Nx1)，必须单调递增
            %   overall_min_bound - 正态分布积分的下限 (-inf 如果需要整个范围)
            %   overall_max_bound - 正态分布积分的上限 (+inf 如果需要整个范围)
            %   dist_mean_mu      - 正态分布的均值 (mu)
            %   dist_std_sig      - 正态分布的标准差 (sigma)
            %
            % 输出:
            %   massV_out         - 每个网格点区间对应的概率质量 (Nx1 列向量)
            %   lbV_out           - 每个区间的下界 (Nx1 列向量)
            %   ubV_out           - 每个区间的上界 (Nx1 列向量)
            %
            % 说明:
            %   对于给定的网格点 x_grid_points = [x1, x2, ..., xN]，此函数将定义 N 个区间。
            %   区间的边界点是网格点之间的中点。
            %   第一个区间的下界是 overall_min_bound，最后一个区间的上界是 overall_max_bound。
            %   然后计算标准正态CDF来得到每个区间的概率质量。
            %   最终的 massV_out 会被归一化，使其总和为1。

            num_points = length(x_grid_points);
            x_grid_points_row = x_grid_points(:)'; % 确保是行向量

            % 检查网格点是否单调递增
            if num_points > 1 && any(x_grid_points_row(2:num_points) < x_grid_points_row(1:num_points-1))
                error('norm_grid: 输入的网格点 x_grid_points 必须单调递增。');
            end

            % 初始化区间边界向量
            lower_bounds = [];
            upper_bounds = [];

            if num_points > 1
                % 计算网格点之间的中点
                mid_points = 0.5 * (x_grid_points_row(1:num_points-1) + x_grid_points_row(2:num_points));
                % 定义区间边界
                lower_bounds = [overall_min_bound, mid_points];
                upper_bounds = [mid_points, overall_max_bound];
            elseif num_points == 1 % 如果只有一个网格点
                lower_bounds = overall_min_bound;
                upper_bounds = overall_max_bound;
            else % 如果没有网格点
                massV_out = [];
                lbV_out = [];
                ubV_out = [];
                return;
            end

            % 计算所有边界点的累积分布函数(CDF)值
            % 包含所有下界点和最后一个上界点
            cdf_values_at_bounds = normcdf([lower_bounds, upper_bounds(end)], dist_mean_mu, dist_std_sig);

            % 计算每个区间的概率质量 (相邻CDF值的差)
            massV_calc = diff(cdf_values_at_bounds);

            % 检查是否有负概率质量 (由于数值误差可能发生)
            if any(massV_calc < -1e-9) % 允许非常小的负值误差
                warning('norm_grid: 检测到负概率质量。已将其设为0。');
                massV_calc(massV_calc < 0) = 0; % 将负质量设为0
            end

            % 归一化概率质量，使其总和为1
            total_sum_mass = sum(massV_calc);
            if total_sum_mass > 1e-9 % 避免除以零
                massV_out = massV_calc / total_sum_mass;
            else % 如果总质量接近于零 (例如，网格范围远偏离分布)
                if num_points > 0
                    massV_out = ones(1, num_points) / num_points; % 使用均匀分布作为备用
                    warning('norm_grid: 总概率质量过小，已使用均匀分布。');
                else
                    massV_out = []; % 如果没有点，则为空
                end
            end

            % 确保输出为列向量
            massV_out = massV_out(:);
            lbV_out = lower_bounds(:);
            ubV_out = upper_bounds(:);

            % 最终检查边界有效性
            if num_points > 1 && any(ubV_out < lbV_out)
                error('norm_grid: 区间上界小于下界，发生错误。');
            end
        end

        function eIdxM_out = MarkovChainSimulation(num_simulations, num_periods_sim, initial_prob_dist_p0V, transition_matrix_trProbM, random_numbers_rvInM)
            % MarkovChainSimulation - 模拟马尔可夫链的路径
            %
            % 输入:
            %   num_simulations       - 模拟的路径数量 (个体数量)
            %   num_periods_sim       - 每条路径的时期长度 (时间长度)
            %   initial_prob_dist_p0V - 初始状态的概率分布 (行或列向量)
            %   transition_matrix_trProbM - 状态转移概率矩阵 (P(i,j) = Pr(s_t+1 = j | s_t = i))
            %   random_numbers_rvInM  - 用于模拟的随机数矩阵 (num_simulations x num_periods_sim)
            %                           每个元素应在 [0, 1) 区间内。
            %
            % 输出:
            %   eIdxM_out             - 模拟的状态索引矩阵 (num_simulations x num_periods_sim)
            %                           元素值为状态的索引 (1 到 num_states)。
            %
            % 说明:
            %   根据初始分布和转移矩阵，使用提供的随机数模拟多条马尔可夫链路径。
            %   状态索引用整数表示。

            num_states = length(initial_prob_dist_p0V); % 获取状态数量

            % 输入校验
            if size(transition_matrix_trProbM,1) ~= num_states || size(transition_matrix_trProbM,2) ~= num_states
                error('MarkovChainSimulation: 转移矩阵维度与初始分布长度不匹配。');
            end
            if abs(sum(initial_prob_dist_p0V) - 1) > 1e-5 % 检查初始分布和是否为1
                warning('MarkovChainSimulation: 初始分布 p0V 的和不为1，已重新归一化。');
                initial_prob_dist_p0V = initial_prob_dist_p0V ./ sum(initial_prob_dist_p0V);
            end
            if any(abs(sum(transition_matrix_trProbM, 2) - 1) > 1e-5) % 检查转移矩阵行和是否为1
                warning('MarkovChainSimulation: 转移矩阵 trProbM 的某些行和不为1，已重新归一化。');
                row_sums_tr = sum(transition_matrix_trProbM, 2);
                row_sums_tr(row_sums_tr <= 1e-9) = 1; % 避免除以零
                transition_matrix_trProbM = bsxfun(@rdivide, transition_matrix_trProbM, row_sums_tr);
            end
            if size(random_numbers_rvInM,1) ~= num_simulations || size(random_numbers_rvInM,2) ~= num_periods_sim
                error('MarkovChainSimulation: 随机数矩阵维度与模拟参数不匹配。');
            end

            % 计算累积概率分布
            cumulative_initial_prob_cP0 = cumsum(initial_prob_dist_p0V(:)'); % 初始分布的CDF (行向量)
            cumulative_transition_prob_cPT = cumsum(transition_matrix_trProbM, 2); % 转移矩阵每行的CDF
            % 确保最后一列的累积概率为1 (处理数值误差)
            if num_states > 0
                cumulative_transition_prob_cPT(:, num_states) = 1.0;
            end


            % 初始化输出的状态索引矩阵 (使用uint16节省内存，如果状态数不多)
            eIdxM_out = zeros(num_simulations, num_periods_sim, 'uint16');

            % 模拟第一个时期的状态
            % 对于每个模拟路径，比较其第一个随机数与cP0，确定初始状态
            % sum(bsxfun(@gt, rv, cdf_vector), 2) + 1 是一种高效的从CDF抽样的方法
            if num_simulations > 0 && num_periods_sim > 0
                eIdxM_out(:, 1) = 1 + sum(bsxfun(@gt, random_numbers_rvInM(:,1), cumulative_initial_prob_cP0), 2);
            end


            % 模拟后续时期的状态
            for t_mc_loop = 1:(num_periods_sim - 1) % 从第1期到倒数第二期，决定下一期状态
                current_state_indices_cSI = eIdxM_out(:, t_mc_loop); % 获取当前所有路径的状态索引

                % 校验当前状态索引的有效性 (理论上应该总是有效)
                valid_indices_vPI = (current_state_indices_cSI >= 1) & (current_state_indices_cSI <= num_states);
                if ~all(valid_indices_vPI) % 如果有无效索引
                    warning('MarkovChainSimulation: 在期 %d 检测到无效的当前状态索引。已重置为状态1。', t_mc_loop);
                    current_state_indices_cSI(~valid_indices_vPI) = 1; % 将无效索引重置为第一个状态
                    eIdxM_out(:, t_mc_loop) = current_state_indices_cSI; % 更新矩阵中的值
                end

                % 获取对应于当前状态的累积转移概率行
                % cPt_for_next_state 的维度是 num_simulations x num_states
                cPt_for_next_state = cumulative_transition_prob_cPT(current_state_indices_cSI, :);

                % 模拟下一期的状态
                eIdxM_out(:, t_mc_loop + 1) = 1 + sum(bsxfun(@gt, random_numbers_rvInM(:, t_mc_loop + 1), cPt_for_next_state), 2);
            end
        end

        function eIdxM = LaborEndowSimulation_olgm(cS, paramS)
            % LaborEndowSimulation_olgm - 模拟个体在整个生命周期内的劳动效率冲击路径
            %
            % 输入:
            %   cS     - 参数结构体，需要包含:
            %            nSim (模拟个体数), aD_orig (年度总年龄跨度)
            %   paramS - 参数结构体，需要包含:
            %            leProb1V (劳动效率初始状态分布)
            %            leTrProbM (劳动效率状态转移矩阵)
            %
            % 输出:
            %   eIdxM  - 模拟的劳动效率状态索引矩阵 (nSim x aD_orig)
            %            每一行代表一个个体，每一列代表一个年度年龄。
            %
            % 说明:
            %   为每个模拟个体，在每个年度年龄上，根据马尔可夫过程模拟其劳动效率状态。
            %   使用固定的随机数种子以保证结果可复现。

            rng(433); % 设置随机数生成器的种子，以保证结果可复现

            % 生成用于模拟的随机数矩阵
            % 维度：个体数量 (cS.nSim) x 年度年龄期数 (cS.aD_orig)
            random_numbers_for_sim = rand([cS.nSim, cS.aD_orig]);

            % 调用马尔可夫链模拟函数
            eIdxM = main_olg_v7_utils.MarkovChainSimulation(cS.nSim, cS.aD_orig, ...
                                                          paramS.leProb1V, paramS.leTrProbM, ...
                                                          random_numbers_for_sim);
            fprintf('劳动禀赋路径已模拟 (%d 个体, %d 年度年龄)。\n', cS.nSim, cS.aD_orig);
        end

        function [HHlaborM_group, L_total_eff_pc] = LaborSupply_Huggett(eIdxM_annual, cS, paramS, Z_ss_norm_group)
            % LaborSupply_Huggett - 计算家庭劳动供给 (模型年龄组层面和总体人均)
            %
            % 输入:
            %   eIdxM_annual    - 年度劳动效率状态索引矩阵 (nSim x aD_orig)
            %   cS              - 参数结构体，包含:
            %                     aD_orig, aD_new, aR_idx_orig, aR_new, physAgeMap, ageEffV_new
            %   paramS          - 参数结构体，包含:
            %                     leGridV (劳动效率网格)
            %   Z_ss_norm_group - 稳态的归一化年龄组人口分布 (列向量)
            %
            % 输出:
            %   HHlaborM_group  - 每个模拟个体在每个模型年龄组的平均劳动供给 (nSim x aD_new)
            %   L_total_eff_pc  - 总体人均有效劳动供给 (标量)
            %
            % 说明:
            %   劳动供给 = 年龄效率 * 随机效率冲击。
            %   首先计算每个个体在每个年度年龄的劳动供给。
            %   然后将其平均到模型年龄组层面。
            %   最后，使用稳态人口分布加权平均，得到总体人均有效劳动供给。
            %   只有工作年龄段的人提供劳动。

            nSim = size(eIdxM_annual, 1); % 获取模拟个体数量

            % 初始化每个个体在每个模型年龄组的劳动供给矩阵
            HHlaborM_group = zeros(nSim, cS.aD_new);

            % 建立年度年龄到模型年龄组的映射 (方便查找)
            ageToGroupMap_local = zeros(cS.aD_orig, 1);
            for a_new_map_idx = 1:cS.aD_new
                annual_indices_in_group = cS.physAgeMap{a_new_map_idx}; % 获取该年龄组包含的年度年龄索引
                if ~isempty(annual_indices_in_group)
                    ageToGroupMap_local(annual_indices_in_group) = a_new_map_idx; % 记录映射
                end
            end

            leGridV_col_local = paramS.leGridV(:); % 劳动效率冲击网格 (列向量)

            % 1. 计算每个个体在每个年度年龄的劳动供给
            HHlaborM_annual_temp = zeros(nSim, cS.aD_orig); % 临时存储年度劳动供给
            for a_orig_idx = 1 : cS.aD_orig % 遍历所有年度年龄
               % 仅对工作年龄段计算劳动供给 (年度模型退休年龄索引前)
               if a_orig_idx < cS.aR_idx_orig
                   a_new_group_idx_current = ageToGroupMap_local(a_orig_idx); % 获取对应的模型年龄组索引
                   if a_new_group_idx_current > 0 && a_new_group_idx_current <= cS.aR_new % 确保是工作年龄组
                       % 劳动供给 = 年龄效率(模型组) * 随机效率(年度个体)
                       % 注意: cS.ageEffV_new 是模型年龄组的平均效率
                       HHlaborM_annual_temp(:, a_orig_idx) = cS.ageEffV_new(a_new_group_idx_current) .* leGridV_col_local(eIdxM_annual(:, a_orig_idx));
                   end
               end
               % 非工作年龄，劳动供给保持为0 (已初始化)
            end

            % 2. 将年度劳动供给平均到模型年龄组
            for a_new_group_idx = 1:cS.aD_new % 遍历每个模型年龄组
                annual_indices_for_this_group = cS.physAgeMap{a_new_group_idx}; % 获取该组包含的年度年龄索引
                if ~isempty(annual_indices_for_this_group)
                    % 取这些年度劳动供给的均值，作为该模型年龄组的劳动供给
                    HHlaborM_group(:, a_new_group_idx) = mean(HHlaborM_annual_temp(:, annual_indices_for_this_group), 2);
                end
            end

            % 3. 计算总体人均有效劳动供给 L_total_eff_pc
            %   对每个工作年龄组，计算该组的平均劳动供给 (跨个体平均)
            %   然后用该组的人口占比加权
            L_total_eff_pc_sum = 0;
            if cS.aR_new > 0 % 确保有工作年龄组
                % 平均每个工作年龄组的劳动供给 (跨个体)
                mean_labor_per_working_group = mean(HHlaborM_group(:, 1:cS.aR_new), 1); % 结果是 1 x cS.aR_new 行向量
                % 与人口分布加权 (Z_ss_norm_group是列向量，转置使其匹配)
                L_total_eff_pc_sum = mean_labor_per_working_group * Z_ss_norm_group(1:cS.aR_new);
            end
            L_total_eff_pc = max(0, L_total_eff_pc_sum); % 确保非负

            fprintf('家庭劳动供给已计算。总体人均有效劳动供给 (L_eff_pc) = %.4f\n', L_total_eff_pc);
        end

        % =====================================================================
        % == 宏观经济函数 ==
        % =====================================================================

        function [R_market_gross_factor, MPL_gross] = HHPrices_Huggett(K_productive, L_total_eff, cS)
            % HHPrices_Huggett - 计算市场总回报率因子和总边际劳动产出 (MPL)
            %
            % 输入:
            %   K_productive  - 生产性资本总量
            %   L_total_eff   - 总有效劳动供给
            %   cS            - 参数结构体，需要包含: A, alpha, ddk
            %
            % 输出:
            %   R_market_gross_factor - 市场毛资本回报率因子 (1 + MPK_gross - ddk)
            %   MPL_gross             - 市场毛边际劳动产出
            %
            % 说明:
            %   假设柯布-道格拉斯生产函数 Y = A * K^alpha * L^(1-alpha)。
            %   MPK_gross 是资本的边际产出 (未扣除折旧)。
            %   MPL_gross 是劳动的边际产出。

            % 处理资本或劳动为非正的边界情况
            if K_productive <= 0
                K_productive=1e-6; % 设置为一个非常小的正数
                warning('HHPrices: K_productive 非正，已重置为 %.1e。', K_productive);
            end
            if L_total_eff <= 0
                L_total_eff=1e-6; % 设置为一个非常小的正数
                warning('HHPrices: L_total_eff 非正，已重置为 %.1e。', L_total_eff);
            end

            % 计算总产出 Y_gross
            Y_gross = cS.A * (K_productive^cS.alpha) * (L_total_eff^(1-cS.alpha));

            % 计算资本的边际产出 (MPK_gross, 未扣除折旧)
            MPK_gross_val = cS.alpha * Y_gross / K_productive;

            % 计算劳动的边际产出 (MPL_gross)
            MPL_gross = (1-cS.alpha) * Y_gross / L_total_eff;

            % 计算市场毛资本回报率因子
            % R_market_gross_factor = 1 + 净回报率 = 1 + (MPK_gross - 折旧率)
            R_market_gross_factor = 1 + MPK_gross_val - cS.ddk;
            % 确保回报率因子至少为1 (即净回报率非负，避免出现负利率导致的问题)
            R_market_gross_factor = max(1.0 + 1e-6, R_market_gross_factor); % 略大于1以防数值问题
        end

        % =====================================================================
        % == 家庭效用函数 ==
        % =====================================================================
        function [muM, utilM] = CES_utility(cM_quantity, sigma_crra, cS_common)
            % CES_utility - 计算CRRA效用函数值及其边际效用
            %
            % 输入:
            %   cM_quantity  - 消费量矩阵或向量
            %   sigma_crra   - CRRA风险厌恶系数 (sigma)
            %   cS_common    - 公共参数结构体，需要包含 cFloor (最低消费)
            %
            % 输出:
            %   muM          - 边际效用矩阵或向量 (c^(-sigma))
            %   utilM        - 效用值矩阵或向量
            %
            % 说明:
            %   U(c) = c^(1-sigma) / (1-sigma)  如果 sigma != 1
            %   U(c) = log(c)                   如果 sigma == 1
            %   消费量 c 会被限制在 cS_common.cFloor 以上。
            %   对于低于 cFloor 的消费，效用值会受到惩罚。

            if ~isscalar(sigma_crra) || sigma_crra <= 0
                error('CES_utility: sigma_crra 必须是正标量。');
            end

            min_c_quantity = cS_common.cFloor; % 获取最低消费量

            % 标记有效的消费量 (大于等于最低消费)
            is_valid_consumption = (cM_quantity >= min_c_quantity);
            % 将消费量限制在最低消费以上 (用于计算)
            c_adjusted_quantity = max(min_c_quantity, cM_quantity);

            % 初始化输出矩阵
            utilM = -Inf(size(cM_quantity)); % 默认效用为负无穷
            muM   =  Inf(size(cM_quantity)); % 默认边际效用为正无穷 (对应极低消费)

            % 计算效用和边际效用
            if abs(sigma_crra - 1) < 1e-6 % 对数效用情况 (sigma ≈ 1)
                utilM(is_valid_consumption) = log(c_adjusted_quantity(is_valid_consumption));
                muM(is_valid_consumption)   = 1 ./ c_adjusted_quantity(is_valid_consumption);
            else % CRRA效用情况 (sigma != 1)
                utilM(is_valid_consumption) = (c_adjusted_quantity(is_valid_consumption).^(1-sigma_crra)) ./ (1-sigma_crra);
                muM(is_valid_consumption)   = c_adjusted_quantity(is_valid_consumption).^(-sigma_crra);
            end

            % 对无效消费 (低于cFloor) 进行惩罚
            % 效用惩罚：一个非常大的负数，且随消费不足程度增加
            utilM(~is_valid_consumption) = -1e10 - (min_c_quantity - cM_quantity(~is_valid_consumption))*1e10;
            % 边际效用惩罚：在正常计算的边际效用基础上增加一个大数
            % (确保在优化时，低于cFloor的消费选项不被选择)
            if abs(sigma_crra - 1) < 1e-6
                 muM(~is_valid_consumption) = 1 ./ c_adjusted_quantity(~is_valid_consumption) + 1e10;
            else
                 muM(~is_valid_consumption) = c_adjusted_quantity(~is_valid_consumption).^(-sigma_crra) + 1e10;
            end
        end

        % =====================================================================
        % == V7 家庭问题求解相关函数 ==
        % =====================================================================

        % --- V7 家庭收入计算 ---
        function [resources_for_c_and_k_prime, labor_income_gross_state, pps_deduction_actual_state] = HHIncome_Huggett(...
                k_now_val, R_k_net_factor, w_gross, ...                     % 当前资产, 资本回报因子, 毛工资率
                TR_total, b_payg_val, c_pps_rule_based_input_val, ...     % 总转移支付, PAYG福利, 规则化PPS缴费输入
                a_idx, paramS_hh, cS, epsilon_val)                        % 年龄组索引, 家庭参数, 公共参数, 效率冲击值

            % HHIncome_Huggett (V7) - 计算家庭可用于消费和非PPS储蓄的总资源
            %
            % 输入:
            %   k_now_val        - 当前非PPS资产水平
            %   R_k_net_factor   - 家庭面临的税后资本净回报因子 (1 + r_net)
            %   w_gross          - 市场毛工资率
            %   TR_total         - 总的政府 lump-sum 转移支付 (人均，如遗赠分配)
            %   b_payg_val       - 当前年龄组的PAYG福利 (如果退休)
            %   c_pps_rule_based_input_val - 基于规则计算出的意向PPS缴费额 (可能还需受限)
            %   a_idx            - 当前模型年龄组索引
            %   paramS_hh        - 特定于家庭的参数结构体，需要:
            %                      tau_l (一般劳动所得税率)
            %                      theta_payg_actual_for_hh (家庭面临的实际PAYG税率)
            %                      pps_tax_deferral_active (PPS税收递延是否激活标志)
            %   cS               - 公共参数结构体，需要:
            %                      aR_new (工作年龄组数), ageEffV_new (年龄效率)
            %                      pps_max_contrib_frac (PPS缴费比例上限-法定)
            %                      pps_annual_contrib_limit (PPS年度绝对缴费上限)
            %   epsilon_val      - 当前个体的劳动效率冲击值 (非对数)
            %
            % 输出:
            %   resources_for_c_and_k_prime - 可用于消费和下一期非PPS储蓄的总资源
            %                                (已扣除实际发生的PPS缴费)
            %   labor_income_gross_state    - 当前状态下的税前劳动总收入
            %   pps_deduction_actual_state  - 实际发生的、符合税前扣除条件的PPS缴费额
            %
            % 说明:
            %   1. 计算税前劳动总收入。
            %   2. 根据规则和法定上限，确定实际的PPS缴费额。
            %   3. 如果PPS税收递延激活，从劳动收入中扣除PPS缴费，得到应税劳动收入。
            %   4. 计算一般劳动所得税 (tau_l) 和 PAYG 工资税 (theta)。
            %   5. 计算税后净劳动收入，加上转移支付和PAYG福利，得到非资本净收入。
            %   6. 计算税后资本收入。
            %   7. 总资源 = 当前资产 + 税后资本收入 + 非资本净收入 - 实际PPS缴费支出。

            labor_income_gross_state = 0;     % 初始化税前劳动总收入
            pps_deduction_actual_state = 0;   % 初始化实际PPS税前扣除额
            non_capital_income = 0;           % 初始化非资本净收入

            % --- A. 计算劳动相关收入和PPS缴费 (仅对工作年龄组) ---
            if a_idx <= cS.aR_new % 如果是工作年龄组
                % A1. 计算税前劳动总收入
                age_efficiency = cS.ageEffV_new(a_idx); % 获取该年龄组的平均年龄效率
                labor_income_gross_state = w_gross * age_efficiency * epsilon_val; % 税前总工资

                % A2. 确定实际的PPS缴费额 (基于规则输入并施加法定上限)
                %     c_pps_rule_based_input_val 是VFI外层根据年龄组规则算出的意向缴费
                c_pps_choice_val_limited_by_rule = max(0, c_pps_rule_based_input_val); % 确保非负
                
                % 法定上限1: 按收入比例
                max_pps_by_statutory_contrib_frac = labor_income_gross_state * cS.pps_max_contrib_frac;
                % 法定上限2: 年度绝对金额
                % (cS.pps_annual_contrib_limit)

                % 实际允许的PPS缴费是规则值和各法定上限的较小者
                actual_pps_contribution_expenditure = min(c_pps_choice_val_limited_by_rule, cS.pps_annual_contrib_limit);
                actual_pps_contribution_expenditure = min(actual_pps_contribution_expenditure, max_pps_by_statutory_contrib_frac);
                actual_pps_contribution_expenditure = max(0, actual_pps_contribution_expenditure); % 再次确保非负

                % A3. 确定税前可扣除的PPS缴费额
                if isfield(paramS_hh, 'pps_tax_deferral_active') && paramS_hh.pps_tax_deferral_active
                    % 如果税收递延激活，则实际发生的PPS缴费可以税前扣除
                    pps_deduction_actual_state = actual_pps_contribution_expenditure;
                else
                    pps_deduction_actual_state = 0; % 否则无税前扣除
                end

                % A4. 计算应纳一般所得税的劳动收入基数
                labor_income_taxable_for_tau_l = labor_income_gross_state - pps_deduction_actual_state;
                labor_income_taxable_for_tau_l = max(0, labor_income_taxable_for_tau_l); % 确保非负

                % A5. 计算各项劳动收入税
                income_tax_tau_l = labor_income_taxable_for_tau_l * paramS_hh.tau_l; % 一般劳动所得税
                payg_tax_theta = labor_income_gross_state * paramS_hh.theta_payg_actual_for_hh; % PAYG工资税 (对总工资征收)

                % A6. 计算税后净劳动收入
                labor_income_net_of_all_taxes = labor_income_gross_state - income_tax_tau_l - payg_tax_theta;

                % A7. 计算总的非资本净收入 (税后劳动收入 + 转移支付 + PAYG福利)
                %     注意：工作年龄组的 b_payg_val 通常为0
                non_capital_income = labor_income_net_of_all_taxes + TR_total + b_payg_val;
            else % 如果是退休年龄组
                % 退休者无劳动收入，只有转移支付和PAYG福利
                % PPS缴费也为0 (由VFI外层传入的 c_pps_rule_based_input_val 应为0)
                actual_pps_contribution_expenditure = 0; % 确保退休期PPS缴费为0
                pps_deduction_actual_state = 0;        % 无劳动收入，无扣除
                non_capital_income = TR_total + b_payg_val;
            end

            % --- B. 计算资本相关收入 ---
            % 家庭面临的税后资本净回报率 r_net = (R_k_net_factor - 1)
            capital_income_net_of_tax = (R_k_net_factor - 1) * k_now_val;

            % --- C. 计算可用于消费和下一期非PPS储蓄的总资源 ---
            % 总资源 = 当前非PPS资产 + 税后资本净收入 + 非资本净收入 - 实际PPS缴费支出
            % 注意: actual_pps_contribution_expenditure 才是实际从预算中扣除用于PPS的金额
            resources_for_c_and_k_prime = k_now_val + capital_income_net_of_tax + non_capital_income - actual_pps_contribution_expenditure;

            % 确保资源非无穷大 (例如，如果输入价格有问题)
            if ~isfinite(resources_for_c_and_k_prime)
                resources_for_c_and_k_prime = -1e10; % 返回一个非常小的值，表示预算极度受限
                warning('HHIncome_Huggett (V7): 计算得到的资源为非有限值。');
            end
        end

        % --- V7 VFI 主函数 (调用按年龄求解的函数) ---
        function [cPolM_q, kPolM, cPpsPolM_rule, valM] = HHSolution_VFI_Huggett(...
                R_k_net_factor_vfi, w_gross_vfi, TR_total_vfi, bV_payg_vfi, paramS_vfi, cS)
            % HHSolution_VFI_Huggett (V7) - 通过逆向迭代求解家庭的多期生命周期问题
            %
            % 输入:
            %   R_k_net_factor_vfi - 家庭面临的税后资本净回报因子 (1 + r_net)
            %   w_gross_vfi        - 市场毛工资率
            %   TR_total_vfi       - 总的政府 lump-sum 转移支付 (人均)
            %   bV_payg_vfi        - 各年龄组的PAYG福利向量 (1 x aD_new)
            %   paramS_vfi         - 特定于家庭的参数结构体 (传递给下一层)
            %   cS                 - 公共参数结构体
            %
            % 输出 (均为 4D 矩阵: nk x nkpps x nw x aD_new):
            %   cPolM_q            - 最优消费量策略
            %   kPolM              - 最优下一期非PPS资产策略
            %   cPpsPolM_rule      - 基于规则的实际PPS缴费额策略 (非选择变量，但记录在此)
            %   valM               - 价值函数
            %
            % 说明:
            %   从最后一个年龄组开始，逆向求解每个年龄组、每个状态下的最优决策。
            %   状态变量包括: 非PPS资产(k), PPS资产(k_pps), 劳动效率冲击(epsilon)。
            %   决策变量: 下一期非PPS资产(k')。PPS缴费(c_pps)在此版本中基于规则决定。
            %   消费(c)通过预算约束内生决定。

            % 初始化策略和价值函数矩阵
            cPolM_q  = zeros(cS.nk, cS.nkpps, cS.nw, cS.aD_new);    % 消费量
            kPolM  = zeros(cS.nk, cS.nkpps, cS.nw, cS.aD_new);      % 非PPS资产 k'
            cPpsPolM_rule = zeros(cS.nk, cS.nkpps, cS.nw, cS.aD_new);% 规则化PPS缴费 c_pps
            valM = -Inf(cS.nk, cS.nkpps, cS.nw, cS.aD_new);         % 价值函数 (初始化为负无穷)

            % fprintf('  VFI V7 (HHSolution_VFI_Huggett): 开始逆向迭代...\n');
            total_age_groups_vfi = cS.aD_new;
            for a_idx = cS.aD_new : -1 : 1 % 从最后一个年龄组向前迭代
                % 打印进度 (可选)
                % if mod(cS.aD_new - a_idx + 1, 1) == 0 || a_idx == cS.aD_new || a_idx == 1
                %     progress_pct_vfi = (cS.aD_new - a_idx + 1) / total_age_groups_vfi * 100;
                %     fprintf('    VFI V7 年龄组 %2d / %2d (大约 %.0f%%)...\n', a_idx, total_age_groups_vfi, progress_pct_vfi);
                % end

                vPrime_kkppse_next = []; % 初始化下一期的期望价值函数
                if a_idx < cS.aD_new % 如果不是最后一个年龄组
                    vPrime_kkppse_next = valM(:,:,:,a_idx+1); % 获取已计算的下一期价值函数
                end

                eps_grid_for_vfi = paramS_vfi.leGridV; % 获取劳动效率冲击网格

                % 调用按年龄组求解的函数
                [cPolM_q(:,:,:,a_idx), kPolM(:,:,:,a_idx), cPpsPolM_rule(:,:,:,a_idx), valM(:,:,:,a_idx)] = ...
                    main_olg_v7_utils.HHSolutionByAge_VFI_Huggett_v7(a_idx, vPrime_kkppse_next, ...
                    R_k_net_factor_vfi, w_gross_vfi, TR_total_vfi, bV_payg_vfi(a_idx), paramS_vfi, cS, eps_grid_for_vfi);
            end
            % fprintf('  VFI V7 (HHSolution_VFI_Huggett): 完成。\n');
        end

        % --- V7 VFI 按年龄组求解 (启用PARFOR) ---
        function [cPol_age_q, kPol_age, cPpsPol_age_rule, val_age] = HHSolutionByAge_VFI_Huggett_v7(...
                a_idx, vPrime_kkppse_next, R_k_net_factor_age, w_gross_age, TR_total_age, b_age_val, ...
                paramS_age, cS, epsilon_grid)
            % HHSolutionByAge_VFI_Huggett_v7 - 求解给定年龄组的家庭问题 (VFI的核心步骤)
            %
            % 输入: (与HHSolution_VFI_Huggett中的参数类似，但针对特定年龄组a_idx)
            %   a_idx                 - 当前年龄组索引
            %   vPrime_kkppse_next    - 下一期期望价值函数 (nk x nkpps x nw 矩阵)
            %   ... (其他价格和参数)
            %   epsilon_grid          - 劳动效率冲击网格 (1 x nw)
            %
            % 输出 (均为 3D 矩阵: nk x nkpps x nw，对应当前年龄组a_idx):
            %   cPol_age_q            - 最优消费量策略
            %   kPol_age              - 最优下一期非PPS资产策略
            %   cPpsPol_age_rule      - 基于规则的实际PPS缴费额
            %   val_age               - 当前期价值函数
            %
            % 说明:
            %   对于给定的年龄组a_idx和所有状态(k, k_pps, epsilon)：
            %   1. 确定本期基于规则的PPS缴费额 c_pps。
            %   2. 计算下一期的PPS资产 k_pps'。
            %   3. 构建用于优化 k' 的目标函数，该函数内嵌了Bellman方程。
            %      目标函数需要下一期期望价值函数的插值器 (基于k'和规则化的k_pps')。
            %   4. 使用fminbnd优化求解最优的k'。
            %   5. 根据最优k'和预算约束计算最优消费c。
            %   核心区别于V6：c_pps不是优化选择，而是基于年龄组规则外生给定(但受法定约束)。
            %   PARFOR用于加速对状态变量k的循环。

            % --- 初始化输出矩阵 (使用临时变量以适应parfor) ---
            cPol_age_q_init = zeros(cS.nk, cS.nkpps, cS.nw);
            kPol_age_init   = zeros(cS.nk, cS.nkpps, cS.nw);
            cPpsPol_age_rule_init = zeros(cS.nk, cS.nkpps, cS.nw); % 存储规则化后的PPS缴费
            val_age_init    = -Inf(cS.nk, cS.nkpps, cS.nw);    % 初始化为负无穷

            % --- 判断当前年龄组是否适用PPS缴费规则 ---
            % model_age_group_start_year_idx_v7 是当前年龄组对应的第一个年度生理年龄的索引 (相对于age1_orig)
            model_age_group_start_year_idx_v7 = cS.physAgeMap{a_idx}(1);
            is_group_eligible_for_pps_contrib_rule = ...
                (a_idx <= cS.aR_new && ... % 是工作年龄组
                 model_age_group_start_year_idx_v7 <= cS.pps_contribution_age_max_idx && ... % 年度年龄未超过最大缴费年龄
                 cS.pps_active && ...        % PPS计划已激活
                 (cS.pps_max_contrib_frac > 0 || cS.pps_annual_contrib_limit > 0) && ... % 法定缴费上限有效
                 cS.pps_fixed_contrib_schedule_frac(a_idx) > 0); % 且本年龄组的规则缴费比例大于0

            % --- fminbnd 优化器选项 ---
            fminbnd_opts_iter = optimset('TolX', cS.fminbnd_TolX, 'Display', cS.fminbnd_Display);

            % --- 情况1: 最后一个年龄组 (a_idx == cS.aD_new) ---
            if a_idx == cS.aD_new
                for ik_temp = 1:cS.nk % 遍历非PPS资产状态
                    for ikpps_temp = 1:cS.nkpps % 遍历PPS资产状态
                        k_now_val_last = cS.kGridV(ik_temp);       % 当前非PPS资产
                        k_pps_now_val_last = cS.kppsGridV(ikpps_temp); % 当前PPS资产
                        for ie_temp = 1:cS.nw % 遍历效率状态
                            % 最后一个时期，无PPS缴费
                            cpps_last_period = 0;
                            cPpsPol_age_rule_init(ik_temp, ikpps_temp, ie_temp) = cpps_last_period;

                            % 计算可用于消费的资源 (已扣除PPS缴费，这里为0)
                            [resources_non_pps_last, ~, ~] = main_olg_v7_utils.HHIncome_Huggett(...
                                k_now_val_last, R_k_net_factor_age, w_gross_age, TR_total_age, b_age_val, ...
                                cpps_last_period, a_idx, paramS_age, cS, epsilon_grid(ie_temp));

                            % 最后一个时期，全部PPS资产被提取 (税前)
                            pps_final_withdrawal_pretax = k_pps_now_val_last;
                            % 计算税后提取额
                            pps_final_withdrawal_net = pps_final_withdrawal_pretax * (1 - cS.pps_tax_rate_withdrawal);
                            % 总资源 = 非PPS资源 + 税后PPS提取额
                            total_resources_for_consumption_last_period = resources_non_pps_last + pps_final_withdrawal_net;

                            % 最优消费 (全部资源用于消费，扣除消费税)
                            cPol_age_q_init(ik_temp,ikpps_temp,ie_temp) = max(cS.cFloor, total_resources_for_consumption_last_period / (1 + cS.tau_c) );
                            % 下一期非PPS资产为最小值 (通常为0)
                            kPol_age_init(ik_temp,ikpps_temp,ie_temp) = cS.kMin;
                            % 计算效用
                            [~, val_age_init(ik_temp,ikpps_temp,ie_temp)] = main_olg_v7_utils.CES_utility(cPol_age_q_init(ik_temp,ikpps_temp,ie_temp), cS.sigma, cS);
                        end
                    end
                end
            else % --- 情况2: 非最后一个年龄组 ---
                % 为下一期期望价值函数创建插值器 (每个当前效率状态 ie_current 对应一个插值器)
                EV_interpolants = cell(cS.nw, 1); % 存储插值器的cell数组
                for ie_current = 1:cS.nw % 对每个当前的效率状态 ie_current
                    EV_for_interp_slice = zeros(cS.nk, cS.nkpps); % 准备用于插值的数据表 (k_next, k_pps_next)
                    % 计算 E[V'(k', k_pps', e') | e_current]
                    %   = sum_{e'} V'(k', k_pps', e') * P(e' | e_current)
                    for ik_next = 1:cS.nk % 遍历下一期非PPS资产状态 k'
                        for ikpps_next = 1:cS.nkpps % 遍历下一期PPS资产状态 k_pps'
                            expected_v_sum = 0;
                            for ie_next = 1:cS.nw % 遍历下一期效率状态 e'
                                expected_v_sum = expected_v_sum + vPrime_kkppse_next(ik_next, ikpps_next, ie_next) * paramS_age.leTrProbM(ie_current, ie_next);
                            end
                            EV_for_interp_slice(ik_next, ikpps_next) = expected_v_sum;
                        end
                    end

                    % 根据网格点数量创建合适的插值器
                    if cS.nk > 1 && cS.nkpps > 1 % 二维插值 (k', k_pps')
                        EV_interpolants{ie_current} = griddedInterpolant({cS.kGridV, cS.kppsGridV}, EV_for_interp_slice, 'linear', 'linear');
                    elseif cS.nk > 1 % 一维插值 (k' only, k_pps' is scalar)
                        EV_interpolants{ie_current} = griddedInterpolant(cS.kGridV, EV_for_interp_slice(:,1), 'linear', 'linear');
                    elseif cS.nkpps > 1 % 一维插值 (k_pps' only, k' is scalar)
                        EV_interpolants{ie_current} = griddedInterpolant(cS.kppsGridV, EV_for_interp_slice(1,:)', 'linear', 'linear');
                    else % 0维插值 (k' and k_pps' are both scalar, V' is a single value)
                        EV_interpolants{ie_current} = @(k_s, kp_s) EV_for_interp_slice(1,1); % 返回常数值的匿名函数
                    end
                end % 结束插值器创建循环

                % --- 使用 PARFOR 循环遍历当前非PPS资产状态 ik ---
                % 为parfor准备cell数组来存储每个ik切片的结果
                c_results_cell_par = cell(cS.nk,1);
                k_results_cell_par = cell(cS.nk,1);
                cpps_results_cell_par = cell(cS.nk,1);
                v_results_cell_par = cell(cS.nk,1);

                parfor ik_par = 1:cS.nk % PARFOR 循环开始
                    % 初始化用于存储当前ik下结果的临时数组 (parfor内部的切片变量)
                    c_slice_for_this_ik = zeros(cS.nkpps, cS.nw);
                    k_slice_for_this_ik = zeros(cS.nkpps, cS.nw);
                    cpps_slice_for_this_ik = zeros(cS.nkpps, cS.nw);
                    v_slice_for_this_ik = -Inf(cS.nkpps, cS.nw);

                    % 将部分外部变量引入parfor内部 (确保parfor能正确访问)
                    % 结构体通常会被整体广播，但为了明确，可以部分复制或传递
                    cS_in_par = cS; % 公共参数
                    paramS_age_in_par = paramS_age; % 年龄特定参数
                    epsilon_grid_in_par = epsilon_grid; % 效率网格
                    R_k_net_factor_age_in_par = R_k_net_factor_age; % 价格
                    w_gross_age_in_par = w_gross_age;
                    TR_total_age_in_par = TR_total_age;
                    b_age_val_in_par = b_age_val;
                    a_idx_in_par = a_idx; % 当前年龄组
                    is_group_eligible_for_pps_in_par = is_group_eligible_for_pps_contrib_rule; % PPS资格
                    fminbnd_opts_in_par = fminbnd_opts_iter; % 优化器选项
                    EV_interpolants_in_par = EV_interpolants; % 插值器cell数组 (会被广播)

                    % --- 内部串行循环遍历 k_pps 和 epsilon 状态 ---
                    for ikpps_srl = 1:cS_in_par.nkpps % 串行循环: 当前PPS资产状态
                        for ie_srl = 1 : cS_in_par.nw % 串行循环: 当前效率状态
                            % 获取当前状态变量值
                            k_state_srl = cS_in_par.kGridV(ik_par); % 当前非PPS资产 (来自parfor的ik_par)
                            k_pps_state_srl = cS_in_par.kppsGridV(ikpps_srl); % 当前PPS资产
                            epsilon_state_srl = epsilon_grid_in_par(ie_srl); % 当前效率值

                            % 1. 确定本期基于规则的实际PPS缴费额 c_pps_effective
                            cpps_effective_contribution_srl = 0; % 初始化
                            if is_group_eligible_for_pps_in_par % 如果该年龄组符合PPS缴费规则
                                % 获取该年龄组的年龄效率 和 规则缴费比例
                                age_efficiency_state_srl = cS_in_par.ageEffV_new(a_idx_in_par);
                                rule_contrib_rate_for_group = cS_in_par.pps_fixed_contrib_schedule_frac(a_idx_in_par);
                                % 计算税前劳动总收入
                                labor_income_gross_state_srl = w_gross_age_in_par * age_efficiency_state_srl * epsilon_state_srl;
                                % 根据规则计算意向缴费额
                                cpps_desired_by_rule_srl = labor_income_gross_state_srl * rule_contrib_rate_for_group;
                                % 施加法定上限1: 工资收入比例
                                max_pps_by_stat_contrib_frac_cap_srl = labor_income_gross_state_srl * cS_in_par.pps_max_contrib_frac;
                                % 施加法定上限2: 年度绝对金额
                                % (cS_in_par.pps_annual_contrib_limit)
                                % 最终有效缴费是规则值和各法定上限的较小者
                                cpps_effective_contribution_srl = min(cpps_desired_by_rule_srl, cS_in_par.pps_annual_contrib_limit);
                                cpps_effective_contribution_srl = min(cpps_effective_contribution_srl, max_pps_by_stat_contrib_frac_cap_srl);
                                cpps_effective_contribution_srl = max(0, cpps_effective_contribution_srl); % 确保非负
                            end
                            cpps_slice_for_this_ik(ikpps_srl, ie_srl) = cpps_effective_contribution_srl; % 存储规则化PPS缴费

                            % 2. 计算下一期的PPS资产 k_pps_prime (基于规则化缴费)
                            pps_withdrawal_pretax_this_period_srl = 0; % 初始化当期PPS提取额 (税前)
                            % 检查是否为退休期且符合提取条件
                            annual_age_check_srl = cS_in_par.physAgeMap{a_idx_in_par}(1); % 当前组的起始年度年龄索引
                            is_retired_group_check_srl = (a_idx_in_par > cS_in_par.aR_new); % 是否为退休年龄组
                            if is_retired_group_check_srl && annual_age_check_srl >= cS_in_par.pps_withdrawal_age_min_idx && cS_in_par.pps_active
                                pps_withdrawal_pretax_this_period_srl = k_pps_state_srl * cS_in_par.pps_withdrawal_rate;
                            end
                            % PPS资产演化方程
                            pps_return_factor_srl = 1 + ( (R_k_net_factor_age_in_par - 1) + cS_in_par.pps_return_rate_premium );
                            k_pps_prime_val_from_rule_srl = (k_pps_state_srl + cpps_effective_contribution_srl - pps_withdrawal_pretax_this_period_srl) * pps_return_factor_srl;
                            % 限制在PPS资产网格内
                            k_pps_prime_val_from_rule_srl = max(cS_in_par.kppsMin, min(cS_in_par.kppsMax, k_pps_prime_val_from_rule_srl));

                            % 3. 计算可用于消费和下一期非PPS储蓄的总资源
                            %    (已扣除实际发生的PPS缴费 cpps_effective_contribution_srl)
                            [resources_for_c_and_k_prime_state_srl, ~, ~] = main_olg_v7_utils.HHIncome_Huggett(...
                                k_state_srl, R_k_net_factor_age_in_par, w_gross_age_in_par, TR_total_age_in_par, b_age_val_in_par, ...
                                cpps_effective_contribution_srl, ... % 使用已确定的有效PPS缴费
                                a_idx_in_par, paramS_age_in_par, cS_in_par, epsilon_state_srl);

                            % 4. 创建用于优化 k_prime 的1D期望价值函数插值器
                            %    它将 k_prime作为输入，使用已确定的 k_pps_prime_val_from_rule_srl
                            ev_func_for_k_prime_opt_1d_srl = @(k_prime_choice_srl) main_olg_v7_utils.CallInterpolator(...
                                EV_interpolants_in_par{ie_srl}, k_prime_choice_srl, k_pps_prime_val_from_rule_srl, cS_in_par);

                            % 5. 优化求解最优的 k_prime (下一期非PPS资产)
                            [c_val_srl, k_p_val_srl, v_val_srl] = main_olg_v7_utils.HHSolutionByOneState_OptK_Mod(...
                                a_idx_in_par, resources_for_c_and_k_prime_state_srl, ... % 预算(已扣除c_pps)
                                ev_func_for_k_prime_opt_1d_srl, ...                      % 1D插值器 E[V'(k', fixed_k_pps')]
                                fminbnd_opts_in_par, cS_in_par, paramS_age_in_par);

                            % 存储结果到当前ik的切片中
                            c_slice_for_this_ik(ikpps_srl, ie_srl) = c_val_srl;
                            k_slice_for_this_ik(ikpps_srl, ie_srl) = k_p_val_srl;
                            v_slice_for_this_ik(ikpps_srl, ie_srl) = v_val_srl;
                        end % 结束串行 ie_srl 循环
                    end % 结束串行 ikpps_srl 循环

                    % 将当前ik_par计算得到的整个 (nkpps x nw) 切片结果存入对应的cell
                    c_results_cell_par{ik_par} = c_slice_for_this_ik;
                    k_results_cell_par{ik_par} = k_slice_for_this_ik;
                    cpps_results_cell_par{ik_par} = cpps_slice_for_this_ik;
                    v_results_cell_par{ik_par} = v_slice_for_this_ik;
                end % 结束 PARFOR ik_par 循环

                % --- 从cell数组组装最终的策略和价值函数矩阵 ---
                for ik_assemble = 1:cS.nk
                    if ~isempty(c_results_cell_par{ik_assemble}) % 检查cell是否为空
                        cPol_age_q_init(ik_assemble,:,:) = c_results_cell_par{ik_assemble};
                        kPol_age_init(ik_assemble,:,:) = k_results_cell_par{ik_assemble};
                        cPpsPol_age_rule_init(ik_assemble,:,:) = cpps_results_cell_par{ik_assemble};
                        val_age_init(ik_assemble,:,:) = v_results_cell_par{ik_assemble};
                    else
                        % 这种情况理论上不应发生，如果发生，说明parfor的某个迭代未能正确赋值
                        warning('HHSolutionByAge_VFI_Huggett_v7: 并行计算结果中发现空单元 (ik_assemble = %d, a_idx = %d)。策略可能不正确。', ik_assemble, a_idx);
                        % val_age_init 已经初始化为 -Inf，其他初始化为0，如果为空则保持这些初始值
                    end
                end
            end % 结束 if a_idx == cS.aD_new / else 条件

            % --- 将初始化后的策略矩阵赋值给函数的正式输出变量 ---
            cPol_age_q = cPol_age_q_init;
            kPol_age = kPol_age_init;
            cPpsPol_age_rule = cPpsPol_age_rule_init;
            val_age = val_age_init;
        end % 结束 HHSolutionByAge_VFI_Huggett_v7 函数


        % --- HHSolutionByOneState_OptK_Mod (优化k', 给定预算和1D E[V']) ---
        function [c_quantity, kPrime_out, ValueFunc_out] = HHSolutionByOneState_OptK_Mod(...
                a_idx_current, budget_for_c_expend_and_kprime, EVprime_of_kprime_interp_1D, ...
                fminbnd_opts_in, cS, paramS_current_age)
            % HHSolutionByOneState_OptK_Mod - 为单个状态优化下一期非PPS资产k'
            %
            % 输入:
            %   a_idx_current            - 当前年龄组索引
            %   budget_for_c_expend_and_kprime - 可用于消费支出和k'的总预算
            %                                   (此预算已扣除当期PPS缴费)
            %   EVprime_of_kprime_interp_1D - 下一期期望价值函数对k'的1D插值器
            %                                   E[V'(k', fixed_k_pps_prime)]
            %   fminbnd_opts_in          - fminbnd优化器的选项
            %   cS                       - 公共参数结构体
            %   paramS_current_age       - 当前年龄的特定参数 (未使用，但保留接口一致性)
            %
            % 输出:
            %   c_quantity    - 最优消费量
            %   kPrime_out    - 最优下一期非PPS资产 k'
            %   ValueFunc_out - 当前状态的价值函数 V(k, k_pps, e)
            %
            % 说明:
            %   给定总预算 (已扣除PPS缴费)，在 k' 的可行范围内，
            %   通过 fminbnd 寻找最大化 Bellman 方程右侧的 k'。
            %   Bellman 方程: V = U(c) + beta * s * E[V'(k', fixed_k_pps')]
            %   其中 c = (budget_for_c_expend_and_kprime - k') / (1 + tau_c)

            % --- 确定 k' 的优化范围 ---
            % k' 的下限是 cS.kMin (通常为0)
            kPrime_min_bound = cS.kMin;
            % k' 的上限是预算扣除最低消费支出后的剩余部分
            % 最低消费支出 = cS.cFloor * (1 + cS.tau_c)
            kPrime_max_bound = budget_for_c_expend_and_kprime - (cS.cFloor * (1 + cS.tau_c));

            ValueFunc_out = -Inf; % 初始化价值函数

            % --- 定义内部嵌套的Bellman目标函数 (用于fminbnd) ---
            function negV_nested = negBellmanObjective_nested_k_only(kPrime_choice_nested)
                % 输入: kPrime_choice_nested (标量, 当前尝试的k')
                % 输出: negV_nested (标量, 对应负的Bellman方程右侧值)
                [negV_nested, ~] = BellmanInner_nested_k_only(kPrime_choice_nested);
            end

            function [negVal_out_nested, Val_out_nested] = BellmanInner_nested_k_only(kPrime_inner_nested)
                % 输入: kPrime_inner_nested (标量, 当前尝试的k')
                % 输出: negVal_out_nested (负价值), Val_out_nested (价值)

                % 1. 计算消费支出和消费量
                consumption_expenditure_nested = budget_for_c_expend_and_kprime - kPrime_inner_nested;
                % 消费量 = 消费支出 / (1 + 消费税率)
                current_c_quantity_nested = max(cS.cFloor, consumption_expenditure_nested / (1 + cS.tau_c) );

                % 2. 计算当期效用 U(c)
                [~, util_current_period_nested] = main_olg_v7_utils.CES_utility(current_c_quantity_nested, cS.sigma, cS);
                if ~isfinite(util_current_period_nested) % 处理可能的无效效用
                    util_current_period_nested = -1e12; % 惩罚
                end
                Val_out_nested = util_current_period_nested; % 初始化总价值

                % 3. 如果不是最后一个时期，加上贴现的未来期望价值
                if a_idx_current < cS.aD_new
                    expected_future_value_nested = -Inf; % 初始化
                    try
                        % 确保 k' 在插值器的有效范围内
                        kPrime_eval_for_interp_nested = max(cS.kGridV(1), min(cS.kGridV(end), kPrime_inner_nested));
                        % 调用1D插值器获取 E[V'(k', fixed_k_pps')]
                        expected_future_value_nested = EVprime_of_kprime_interp_1D(kPrime_eval_for_interp_nested);
                    catch ME_interp_k_only
                        warning('BellmanInner_nested_k_only: 插值错误 (年龄 %d, kP %.2e): %s. 将使用边界值。', ...
                                a_idx_current, kPrime_inner_nested, ME_interp_k_only.message);
                        % 备用：如果插值失败，尝试使用网格边界值
                        if kPrime_inner_nested < cS.kGridV(1)
                            expected_future_value_nested = EVprime_of_kprime_interp_1D(cS.kGridV(1));
                        else
                            expected_future_value_nested = EVprime_of_kprime_interp_1D(cS.kGridV(end));
                        end
                    end
                    if ~isfinite(expected_future_value_nested) % 处理可能的无效期望价值
                        expected_future_value_nested = -1e12; % 惩罚
                    end

                    % 年龄组间的年度化转移存活率 (s_1yr_transitionV)
                    s_transition_to_next_year_nested = cS.s_1yr_transitionV(a_idx_current);
                    % 总价值 = U(c) + beta * s * E[V']
                    Val_out_nested = util_current_period_nested + cS.beta * s_transition_to_next_year_nested * expected_future_value_nested;
                end % 结束 if a_idx_current < cS.aD_new

                % 准备输出 (fminbnd最小化负价值)
                if ~isfinite(Val_out_nested)
                    negVal_out_nested = 1e12;  % 对应非常差的价值
                    Val_out_nested = -1e12;
                else
                    negVal_out_nested = -Val_out_nested; % fminbnd 需要最小化的目标
                end
            end % 结束 BellmanInner_nested_k_only

            % --- 主优化逻辑 ---
            if kPrime_max_bound <= kPrime_min_bound % 如果优化区间无效 (例如预算不足以支付最低消费)
                % 角点解：将所有可用预算（如果有）用于k'的下限，剩余的用于消费
                kPrime_out = kPrime_min_bound; % k' 取下限
                c_expenditure_at_corner = budget_for_c_expend_and_kprime - kPrime_out;
                c_quantity = max(cS.cFloor, c_expenditure_at_corner / (1 + cS.tau_c) ); % 计算对应消费

                % 计算此角点解的价值
                [~, utility_at_corner] = main_olg_v7_utils.CES_utility(c_quantity, cS.sigma, cS);
                if a_idx_current < cS.aD_new
                    EV_at_corner = -Inf;
                    try
                        kPrime_eval_corner = max(cS.kGridV(1),min(cS.kGridV(end),kPrime_out));
                        EV_at_corner = EVprime_of_kprime_interp_1D(kPrime_eval_corner);
                    catch % 备用插值
                        if kPrime_out<cS.kGridV(1), EV_at_corner=EVprime_of_kprime_interp_1D(cS.kGridV(1));
                        else, EV_at_corner=EVprime_of_kprime_interp_1D(cS.kGridV(end));end;
                    end
                    if ~isfinite(EV_at_corner), EV_at_corner=-1e12; end
                    s_transition_osa_k = cS.s_1yr_transitionV(a_idx_current);
                    ValueFunc_out = utility_at_corner + cS.beta*s_transition_osa_k*EV_at_corner;
                else % 如果是最后时期
                    ValueFunc_out = utility_at_corner;
                end
                if ~isfinite(ValueFunc_out), ValueFunc_out = -1e12; end
            else % 如果优化区间有效
                % 进一步约束优化区间在合理的资产范围内
                kPrime_min_for_fminbnd = max(cS.kMin, kPrime_min_bound); % 不低于全局kMin
                % 确保上限至少比下限大一点点，且不高于全局kMax
                kPrime_max_for_fminbnd = max(kPrime_min_for_fminbnd + 1e-9, min(cS.kMax, kPrime_max_bound));

                if kPrime_min_for_fminbnd >= kPrime_max_for_fminbnd % 如果调整后区间仍无效
                    kPrime_optimal_fminbnd = kPrime_min_for_fminbnd; % 取下限
                    [negValue_fminbnd, ~] = BellmanInner_nested_k_only(kPrime_optimal_fminbnd);
                    ValueFunc_out = -negValue_fminbnd;
                else % 区间有效，调用fminbnd
                    objective_func_for_fminbnd = @(kP_choice) negBellmanObjective_nested_k_only(kP_choice);
                    [kPrime_optimal_fminbnd, negValue_fminbnd, exitflag_fminbnd] = ...
                        fminbnd(objective_func_for_fminbnd, kPrime_min_for_fminbnd, kPrime_max_for_fminbnd, fminbnd_opts_in);

                    % 检查fminbnd的解是否在边界，或者是否成功找到内部解
                    if exitflag_fminbnd <= 0 || ... % fminbnd未成功，或解在边界
                       abs(kPrime_optimal_fminbnd - kPrime_min_for_fminbnd) < 1e-7 || ...
                       abs(kPrime_optimal_fminbnd - kPrime_max_for_fminbnd) < 1e-7
                        % 如果解在边界，比较两个边界点的价值
                        [negValue_at_min_bound, ~] = BellmanInner_nested_k_only(kPrime_min_for_fminbnd);
                        [negValue_at_max_bound, ~] = BellmanInner_nested_k_only(kPrime_max_for_fminbnd);
                        if negValue_at_min_bound <= negValue_at_max_bound + 1e-9 % 比较负价值 (越小越好)
                            kPrime_optimal_fminbnd = kPrime_min_for_fminbnd;
                            negValue_fminbnd = negValue_at_min_bound;
                        else
                            kPrime_optimal_fminbnd = kPrime_max_for_fminbnd;
                            negValue_fminbnd = negValue_at_max_bound;
                        end
                    end
                    ValueFunc_out = -negValue_fminbnd; % 存储最优价值
                end
                kPrime_out = kPrime_optimal_fminbnd; % 存储最优k'
                % 根据最优k'计算消费
                c_expenditure_optimal = budget_for_c_expend_and_kprime - kPrime_out;
                c_quantity = max(cS.cFloor, c_expenditure_optimal / (1 + cS.tau_c) );
                if ~isfinite(ValueFunc_out), ValueFunc_out = -1e12; end
            end % 结束主优化逻辑

            % --- 后处理和最终校验 ---
            % 确保k'在全局范围内
            kPrime_out = max(cS.kMin, min(cS.kMax, kPrime_out));
            % 重新计算最终消费量 (基于可能调整过的kPrime_out)
            c_expenditure_final_check = budget_for_c_expend_and_kprime - kPrime_out;
            c_quantity = max(cS.cFloor, c_expenditure_final_check / (1 + cS.tau_c) );

            % 重新计算最终价值函数值 (确保与最终的c_quantity和kPrime_out一致)
            [~, utility_final_check] = main_olg_v7_utils.CES_utility(c_quantity, cS.sigma, cS);
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

            % 确保输出值是有限的
            if ~isfinite(kPrime_out), kPrime_out = cS.kMin; end
            if ~isfinite(c_quantity), c_quantity = cS.cFloor; end
            if ~isfinite(ValueFunc_out), ValueFunc_out = -1e12; end
        end % 结束 HHSolutionByOneState_OptK_Mod

        % --- HHSimulation_olgm (家庭模拟，V7版本) ---
        function [kHistM_out, kPpsHistM_out, cHistM_out] = HHSimulation_olgm(...
                kPolM_4D_input, cPpsPolM_rule_4D_input, cPolM_consump_q_4D_input, eIdxM_annual_input, ...
                R_k_net_factor_hh_sim, w_gross_sim_price, TR_total_sim_transfer, bV_payg_sim_benefit, ...
                paramS_sim_household, cS_common_sim)
            % HHSimulation_olgm (V7) - 模拟家庭在生命周期内的资产和消费路径
            %
            % 输入:
            %   kPolM_4D_input          - 非PPS资产策略 (nk x nkpps x nw x aD_new)
            %   cPpsPolM_rule_4D_input  - 规则化PPS缴费策略 (nk x nkpps x nw x aD_new)
            %   cPolM_consump_q_4D_input- 消费量策略 (nk x nkpps x nw x aD_new)
            %   eIdxM_annual_input      - 年度劳动效率状态索引 (nSim x aD_orig)
            %   R_k_net_factor_hh_sim   - 家庭税后资本净回报因子
            %   w_gross_sim_price       - 市场毛工资率
            %   TR_total_sim_transfer   - 总转移支付 (人均)
            %   bV_payg_sim_benefit     - PAYG福利向量 (按年龄组)
            %   paramS_sim_household    - 家庭特定参数 (tau_l, theta_payg_actual_for_hh, etc.)
            %   cS_common_sim           - 公共参数
            %
            % 输出:
            %   kHistM_out    - 非PPS资产历史 (nSim x aD_orig)
            %   kPpsHistM_out - PPS资产历史 (nSim x aD_orig)
            %   cHistM_out    - 消费量历史 (nSim x aD_orig)
            %
            % 说明:
            %   根据给定的策略函数、价格、转移支付和效率路径，模拟个体的决策。
            %   首先创建策略函数的插值器。
            %   然后逐个年度年龄进行模拟。在每个年度年龄：
            %     1. 获取当前状态 (k, k_pps, e)。
            %     2. 使用插值器得到策略对应的 c_pps(规则), k', c。
            %     3. 根据规则和法定上限调整实际的 c_pps。
            %     4. 更新下一期的 k_pps 和 k。
            %     5. 记录消费 c。
            %   注意，策略函数是基于模型年龄组的，模拟是基于年度年龄的。

            % --- 1. 准备工作 ---
            % 建立年度年龄到模型年龄组的映射
            ageToGroupMap_sim = zeros(cS_common_sim.aD_orig,1);
            for a_map_idx_sim = 1:cS_common_sim.aD_new
                idx_map_sim = cS_common_sim.physAgeMap{a_map_idx_sim};
                if ~isempty(idx_map_sim)
                    ageToGroupMap_sim(idx_map_sim) = a_map_idx_sim;
                end
            end

            nSim_sim = size(eIdxM_annual_input,1); % 模拟个体数量
            % 初始化历史数据矩阵 (多存储一期用于k_next的初始赋值)
            kHistM_out_temp    = zeros(nSim_sim, cS_common_sim.aD_orig + 1);
            kPpsHistM_out_temp = zeros(nSim_sim, cS_common_sim.aD_orig + 1);
            cHistM_out    = zeros(nSim_sim, cS_common_sim.aD_orig); % 消费只记录aD_orig期

            leGridV_col_sim = paramS_sim_household.leGridV(:); % 效率冲击值网格

            % --- 2. 创建策略函数的插值器 ---
            % (每个效率状态ie 和 每个模型年龄组ia 对应一个插值器)
            kPolInterp_sim         = cell(cS_common_sim.nw, cS_common_sim.aD_new);
            cPpsPolInterp_rule_sim = cell(cS_common_sim.nw, cS_common_sim.aD_new);
            cPolqInterp_sim        = cell(cS_common_sim.nw, cS_common_sim.aD_new);

            for ia_interp = 1:cS_common_sim.aD_new % 遍历模型年龄组
                for ie_interp = 1:cS_common_sim.nw % 遍历效率状态
                    % 根据网格维度选择合适的插值方式
                    if cS_common_sim.nk > 1 && cS_common_sim.nkpps > 1 % 2D插值 (k, k_pps)
                        kPolInterp_sim{ie_interp,ia_interp}    = griddedInterpolant({cS_common_sim.kGridV, cS_common_sim.kppsGridV}, squeeze(kPolM_4D_input(:,:,ie_interp,ia_interp)), 'linear', 'nearest');
                        cPpsPolInterp_rule_sim{ie_interp,ia_interp} = griddedInterpolant({cS_common_sim.kGridV, cS_common_sim.kppsGridV}, squeeze(cPpsPolM_rule_4D_input(:,:,ie_interp,ia_interp)), 'linear', 'nearest');
                        cPolqInterp_sim{ie_interp,ia_interp}   = griddedInterpolant({cS_common_sim.kGridV, cS_common_sim.kppsGridV}, squeeze(cPolM_consump_q_4D_input(:,:,ie_interp,ia_interp)), 'linear', 'nearest');
                    elseif cS_common_sim.nk > 1 && cS_common_sim.nkpps == 1 % 1D插值 (k only)
                        kPolInterp_sim{ie_interp,ia_interp}    = griddedInterpolant(cS_common_sim.kGridV, squeeze(kPolM_4D_input(:,1,ie_interp,ia_interp)), 'linear', 'nearest');
                        cPpsPolInterp_rule_sim{ie_interp,ia_interp} = griddedInterpolant(cS_common_sim.kGridV, squeeze(cPpsPolM_rule_4D_input(:,1,ie_interp,ia_interp)), 'linear', 'nearest');
                        cPolqInterp_sim{ie_interp,ia_interp}   = griddedInterpolant(cS_common_sim.kGridV, squeeze(cPolM_consump_q_4D_input(:,1,ie_interp,ia_interp)), 'linear', 'nearest');
                    elseif cS_common_sim.nk == 1 && cS_common_sim.nkpps > 1 % 1D插值 (k_pps only)
                        kPolInterp_sim{ie_interp,ia_interp}    = griddedInterpolant(cS_common_sim.kppsGridV, squeeze(kPolM_4D_input(1,:,ie_interp,ia_interp))', 'linear', 'nearest');
                        cPpsPolInterp_rule_sim{ie_interp,ia_interp} = griddedInterpolant(cS_common_sim.kppsGridV, squeeze(cPpsPolM_rule_4D_input(1,:,ie_interp,ia_interp))', 'linear', 'nearest');
                        cPolqInterp_sim{ie_interp,ia_interp}   = griddedInterpolant(cS_common_sim.kppsGridV, squeeze(cPolM_consump_q_4D_input(1,:,ie_interp,ia_interp))', 'linear', 'nearest');
                    elseif cS_common_sim.nk == 1 && cS_common_sim.nkpps == 1 % 0D插值 (标量)
                        kPolInterp_sim{ie_interp,ia_interp}    = @(x,y) kPolM_4D_input(1,1,ie_interp,ia_interp); % 匿名函数返回标量值
                        cPpsPolInterp_rule_sim{ie_interp,ia_interp} = @(x,y) cPpsPolM_rule_4D_input(1,1,ie_interp,ia_interp);
                        cPolqInterp_sim{ie_interp,ia_interp}   = @(x,y) cPolM_consump_q_4D_input(1,1,ie_interp,ia_interp);
                    else
                        error('HHSimulation_olgm (V7): 非PPS资产网格nk或PPS资产网格nkpps为零，无法创建插值器。');
                    end
                end
            end % 结束插值器创建

            % PPS资产的年化净回报因子
            pps_return_net_annual_factor_sim = 1 + ((R_k_net_factor_hh_sim - 1) + cS_common_sim.pps_return_rate_premium);

            % --- 3. 按年度年龄进行生命周期模拟 ---
            for a_orig_loop_idx = 1:cS_common_sim.aD_orig % 遍历每个年度年龄
                a_new_group_idx_sim = ageToGroupMap_sim(a_orig_loop_idx); % 获取对应的模型年龄组

                % 获取当前年度的资产状态 (来自上一年度的决策结果)
                kNowV_annual_sim    = kHistM_out_temp(:, a_orig_loop_idx);
                kPpsNowV_annual_sim = kPpsHistM_out_temp(:, a_orig_loop_idx);

                % 初始化当期决策结果向量
                kNextNonPpsV_from_policy = zeros(nSim_sim,1);
                cPpsDecisionFromPol_rule = zeros(nSim_sim,1); % 来自策略的规则化PPS缴费 (可能还需约束)
                cConsumpValV_q_from_policy = zeros(nSim_sim,1);
                kPpsNextV_annual_sim     = zeros(nSim_sim,1); % 下一期PPS资产

                % 获取当前年度年龄 和 相关资格判断
                % current_annual_age_val_sim = cS_common_sim.age1_orig + a_orig_loop_idx - 1; % 当前生理年龄
                % model_age_group_start_year_idx_sim = cS_common_sim.physAgeMap{a_new_group_idx_sim}(1);

                is_working_age_annual_sim = (a_orig_loop_idx < cS_common_sim.aR_idx_orig); % 是否为工作年龄 (基于年度索引)
                % PPS缴费资格 (基于年度年龄索引 和 模型年龄组规则)
                is_pps_contrib_eligible_annual_sim = (is_working_age_annual_sim && ...
                    a_orig_loop_idx <= cS_common_sim.pps_contribution_age_max_idx && ... % 未超年度最大缴费年龄
                    cS_common_sim.pps_active && ...
                    cS_common_sim.pps_fixed_contrib_schedule_frac(a_new_group_idx_sim) > 0); % 且该组规则缴费>0

                % PPS提取资格 (基于年度年龄索引)
                is_pps_withdrawal_eligible_annual_sim = (~is_working_age_annual_sim && ...
                    cS_common_sim.pps_active && ...
                    a_orig_loop_idx >= cS_common_sim.pps_withdrawal_age_min_idx); % 已达年度最低提取年龄

                % 计算当期PPS提取额 (税前)
                pps_withdrawal_pretax_this_year_sim = zeros(nSim_sim,1);
                if is_pps_withdrawal_eligible_annual_sim
                    pps_withdrawal_pretax_this_year_sim = kPpsNowV_annual_sim * cS_common_sim.pps_withdrawal_rate;
                end

                % 存储最终实际发生的PPS缴费 (在所有约束施加后)
                actual_cpps_final_for_period_sim = zeros(nSim_sim,1);

                % --- 3a. 对每个效率状态分别插值获取策略 ---
                for ie_sim_idx = 1 : cS_common_sim.nw % 遍历效率状态
                    simIdxV_for_this_e = find(eIdxM_annual_input(:, a_orig_loop_idx) == ie_sim_idx); % 找到当前为此效率状态的个体
                    if isempty(simIdxV_for_this_e), continue; end % 如果没有个体为此状态，则跳过

                    % 获取这些个体的当前资产状态，并限制在网格范围内
                    kNow_clamped    = max(cS_common_sim.kGridV(1), min(cS_common_sim.kGridV(end), kNowV_annual_sim(simIdxV_for_this_e)));
                    kPpsNow_clamped = max(cS_common_sim.kppsGridV(1), min(cS_common_sim.kppsGridV(end), kPpsNowV_annual_sim(simIdxV_for_this_e)));

                    % 使用插值器获取策略值
                    if cS_common_sim.nk > 1 && cS_common_sim.nkpps > 1 % 2D 插值
                        kNextNonPpsV_from_policy(simIdxV_for_this_e) = kPolInterp_sim{ie_sim_idx, a_new_group_idx_sim}(kNow_clamped, kPpsNow_clamped);
                        cPpsDecisionFromPol_rule(simIdxV_for_this_e) = cPpsPolInterp_rule_sim{ie_sim_idx, a_new_group_idx_sim}(kNow_clamped, kPpsNow_clamped);
                        cConsumpValV_q_from_policy(simIdxV_for_this_e) = cPolqInterp_sim{ie_sim_idx, a_new_group_idx_sim}(kNow_clamped, kPpsNow_clamped);
                    elseif cS_common_sim.nk > 1 % 1D 插值 (k only)
                        kNextNonPpsV_from_policy(simIdxV_for_this_e) = kPolInterp_sim{ie_sim_idx, a_new_group_idx_sim}(kNow_clamped);
                        cPpsDecisionFromPol_rule(simIdxV_for_this_e) = cPpsPolInterp_rule_sim{ie_sim_idx, a_new_group_idx_sim}(kNow_clamped);
                        cConsumpValV_q_from_policy(simIdxV_for_this_e) = cPolqInterp_sim{ie_sim_idx, a_new_group_idx_sim}(kNow_clamped);
                    elseif cS_common_sim.nkpps > 1 % 1D 插值 (k_pps only)
                        kNextNonPpsV_from_policy(simIdxV_for_this_e) = kPolInterp_sim{ie_sim_idx, a_new_group_idx_sim}(kPpsNow_clamped);
                        cPpsDecisionFromPol_rule(simIdxV_for_this_e) = cPpsPolInterp_rule_sim{ie_sim_idx, a_new_group_idx_sim}(kPpsNow_clamped);
                        cConsumpValV_q_from_policy(simIdxV_for_this_e) = cPolqInterp_sim{ie_sim_idx, a_new_group_idx_sim}(kPpsNow_clamped);
                    else % 0D 插值 (标量)
                        kNextNonPpsV_from_policy(simIdxV_for_this_e) = kPolInterp_sim{ie_sim_idx, a_new_group_idx_sim}(kNow_clamped, kPpsNow_clamped); % 实际是调用匿名函数
                        cPpsDecisionFromPol_rule(simIdxV_for_this_e) = cPpsPolInterp_rule_sim{ie_sim_idx, a_new_group_idx_sim}(kNow_clamped, kPpsNow_clamped);
                        cConsumpValV_q_from_policy(simIdxV_for_this_e) = cPolqInterp_sim{ie_sim_idx, a_new_group_idx_sim}(kNow_clamped, kPpsNow_clamped);
                    end

                    % --- 3b. 根据规则和法定上限，确定这些个体实际的PPS缴费 ---
                    temp_actual_cpps_for_these_individuals = zeros(length(simIdxV_for_this_e),1);
                    if is_pps_contrib_eligible_annual_sim % 如果符合缴费资格
                        % 获取这些个体的税前劳动总收入
                        current_gross_labor_income_for_these_individuals = w_gross_sim_price * cS_common_sim.ageEffV_new(a_new_group_idx_sim) * leGridV_col_sim(ie_sim_idx);
                        % 从策略中获取的规则化意向缴费 (对这些个体是相同的值，因为k,k_pps对c_pps(rule)影响不大，主要看年龄组和e)
                        rule_based_cpps_from_policy_for_these = cPpsDecisionFromPol_rule(simIdxV_for_this_e);

                        % 法定上限1: 工资收入比例
                        max_cpps_by_income_frac_for_these = current_gross_labor_income_for_these_individuals * cS_common_sim.pps_max_contrib_frac;
                        % 法定上限2: 年度绝对金额 (cS_common_sim.pps_annual_contrib_limit)

                        % 最终有效缴费是规则值和各法定上限的较小者
                        temp_actual_cpps_for_these_individuals = min(rule_based_cpps_from_policy_for_these, cS_common_sim.pps_annual_contrib_limit);
                        temp_actual_cpps_for_these_individuals = min(temp_actual_cpps_for_these_individuals, max_cpps_by_income_frac_for_these);
                        temp_actual_cpps_for_these_individuals = max(0, temp_actual_cpps_for_these_individuals); % 确保非负
                    end
                    actual_cpps_final_for_period_sim(simIdxV_for_this_e) = temp_actual_cpps_for_these_individuals;
                end % 结束效率状态循环 ie_sim_idx

                % --- 3c. 更新下一期PPS资产 (所有个体) ---
                if cS_common_sim.pps_active
                    kPpsNextV_annual_sim = (kPpsNowV_annual_sim + actual_cpps_final_for_period_sim - pps_withdrawal_pretax_this_year_sim) * pps_return_net_annual_factor_sim;
                    % 限制在PPS资产网格范围内
                    kPpsNextV_annual_sim = max(cS_common_sim.kppsMin, min(cS_common_sim.kppsMax, kPpsNextV_annual_sim));
                else % 如果PPS未激活
                    kPpsNextV_annual_sim = kPpsNowV_annual_sim; % PPS资产不变
                end

                % --- 3d. 存储当期决策和下一期状态 ---
                % 下一期非PPS资产 (限制在网格内)
                kHistM_out_temp(:, a_orig_loop_idx + 1) = max(cS_common_sim.kMin, min(cS_common_sim.kMax, kNextNonPpsV_from_policy));
                % 下一期PPS资产
                kPpsHistM_out_temp(:, a_orig_loop_idx + 1) = kPpsNextV_annual_sim; % 已在上面约束过
                % 当期消费量 (确保不低于cFloor)
                cHistM_out(:, a_orig_loop_idx) = max(cS_common_sim.cFloor, cConsumpValV_q_from_policy);

            end % 结束年度年龄循环 a_orig_loop_idx

            % 截取最终的有效历史数据 (前aD_orig期)
            kHistM_out = kHistM_out_temp(:, 1:cS_common_sim.aD_orig);
            kPpsHistM_out = kPpsHistM_out_temp(:, 1:cS_common_sim.aD_orig);
            % cHistM_out 已经是 aD_orig 期

            % fprintf('家庭生命周期模拟完成 (V7)。\n');
        end


        % --- 核心均衡求解器 (solve_K_tau_l_for_rho_prime, V7版本) ---
        function [K_sol_out, tau_l_sol_out, gbc_res_final_out, converged_and_feasible_out, solution_details_out] = solve_K_tau_l_for_rho_prime(...
                rho_prime_payg_target_input, K_init_guess_input, cS_global, paramS_global_in, eIdxM_global_sim_paths)
            % solve_K_tau_l_for_rho_prime (V7) - 求解给定PAYG替代率下的均衡K和tau_l
            %
            % 输入:
            %   rho_prime_payg_target_input - 目标PAYG替代率 (b_payg / avg_worker_gross_wage)
            %   K_init_guess_input        - 总资本存量的初始猜测值
            %   cS_global                 - 全局公共参数结构体
            %   paramS_global_in          - 全局派生参数结构体 (包含人口、劳动供给等)
            %   eIdxM_global_sim_paths    - 预先模拟的个体劳动效率路径
            %
            % 输出:
            %   K_sol_out                 - 均衡的总资本存量
            %   tau_l_sol_out             - 均衡的一般劳动所得税率 (使GBC平衡)
            %   gbc_res_final_out         - 最终的政府一般预算残差
            %   converged_and_feasible_out- 均衡是否找到且可行的标志
            %   solution_details_out      - 包含均衡状态详细信息的结构体
            %
            % 说明:
            %   这是一个迭代过程，用于寻找一对 (K, tau_l) 使得：
            %     1. 资本市场出清：家庭最优决策汇总得到的总资本需求等于K。
            %     2. 政府一般预算平衡(GBC)：在tau_l调整下，政府非PAYG收支平衡 (TR_gov=0)。
            %     3. PAYG系统约束：实际PAYG税率不超过上限，且能支持目标替代率。
            %   迭代变量是 K_guess 和 tau_l_guess。
            %   内含VFI求解、微观模拟、价格更新、税率调整等步骤。

            % --- 1. 初始化迭代变量和参数 ---
            K_current_guess = K_init_guess_input;
            tau_l_current_guess = cS_global.tau_l_init_guess;
            L_per_capita_global = paramS_global_in.L_per_capita; % 人均有效劳动供给
            mass_workers_global = paramS_global_in.mass_workers_group; % 工人人口占比

            % 获取迭代控制参数
            maxIter_ktl_loop = cS_global.max_iter_K_tau_l;
            tol_norm_ktl_loop = cS_global.tol_K_tau_l;
            dampK_ktl_loop = cS_global.damp_K_v5;
            damp_tau_l_ktl_loop = cS_global.damp_tau_l_v5;

            % 初始化输出和状态变量
            converged_and_feasible_out = false;
            K_sol_out = NaN; tau_l_sol_out = NaN; gbc_res_final_out = Inf;
            solution_details_out = struct(); % 用于存储详细的均衡结果

            % --- 2. 计算目标替代率rho'对应的理论PAYG税率theta_req ---
            mass_retirees_global = sum(paramS_global_in.ageMassV(cS_global.aR_new + 1 : cS_global.aD_new));
            theta_payg_required_calc = 0;
            if mass_workers_global > 1e-9 % 避免除以零
                theta_payg_required_calc = rho_prime_payg_target_input * (mass_retirees_global / mass_workers_global);
            else % 如果没有工人
                if rho_prime_payg_target_input > 1e-9, theta_payg_required_calc = Inf; % 无法实现正替代率
                else, theta_payg_required_calc = 0; % 零替代率可以实现
                end
            end
            theta_payg_required_calc = max(0, theta_payg_required_calc); % 确保非负
            solution_details_out.theta_payg_required_before_cap = theta_payg_required_calc; % 存入细节

            % --- 2a. 预先检查理论theta_req是否已超过上限 ---
            if theta_payg_required_calc > cS_global.theta_payg_max + 1e-5 % 允许微小误差
                % 如果理论所需税率已超上限，则此rho'目标不可行
                if ~isfield(paramS_global_in, 'suppress_initial_theta_print') || ~paramS_global_in.suppress_initial_theta_print
                    fprintf('  solve_K_tau_l (V7): rho_prime_target=%.4f 导致理论theta_req=%.4f > theta_max=%.3f. 直接标记为不可行。\n', ...
                        rho_prime_payg_target_input, theta_payg_required_calc, cS_global.theta_payg_max);
                end
                converged_and_feasible_out = false;
                K_sol_out = K_init_guess_input; tau_l_sol_out = tau_l_current_guess; gbc_res_final_out = Inf;
                solution_details_out.theta_payg = min(theta_payg_required_calc, cS_global.theta_payg_max); % 记录受限的theta
                % 填充其他NaN细节以避免后续错误
                solution_details_out.MPL_gross = NaN; solution_details_out.R_mkt_gross_factor = NaN; solution_details_out.b_payg = NaN;
                solution_details_out.T_bequest_Model = NaN; solution_details_out.C_model = NaN; solution_details_out.Y_model = NaN;
                solution_details_out.K_model_pps = NaN; solution_details_out.K_model_non_pps = NaN;
                return; % 直接返回，无需迭代
            end

            % --- 3. 初始化迭代跟踪变量 ---
            stagnation_counter_ktl = 0;         % 停滞计数器
            prev_devNorm_ktl = Inf;             % 上一轮的偏差范数
            tau_l_boundary_strike_count_ktl = 0;% tau_l触界计数器

            % --- 3a. 打印迭代表头 (如果允许) ---
            if ~isfield(paramS_global_in, 'suppress_inner_print_header') || ~paramS_global_in.suppress_inner_print_header
                fprintf('  solve_K_tau_l_for_rho_prime (V7): rho_prime_target=%.4f (理论theta_req=%.4f), K_init=%.2f, tau_l_init=%.3f\n', ...
                        rho_prime_payg_target_input, theta_payg_required_calc, K_current_guess, tau_l_current_guess);
                fprintf('  IterKTL | K_guess  | tau_l_gs | MPL_g    | theta_act| K_tot_mod| K_pps_mod| GBC_res  | K_dev    | tau_l_dev| Norm     | Improv   | Strikes\n');
                fprintf('  -------------------------------------------------------------------------------------------------------------------------------------------\n');
            end

            % 初始化循环内会计算的变量，以确保在所有分支中都定义
            MPL_gross_iter_val = NaN; R_mkt_gross_factor_iter_val = NaN; theta_payg_actual_iter_val = NaN; b_payg_iter_val = NaN;
            T_bequest_model_iter_val = NaN; C_model_iter_val = NaN; Y_for_gbc_iter_val = NaN; gbc_residual_iter_val = Inf;
            K_model_from_sim_iter_val = K_current_guess; K_dev_from_sim_iter_val = Inf;
            K_model_pps_sim_iter_val = NaN; K_model_nonpps_sim_iter_val = NaN;

            % --- 4. 开始 K 和 tau_l 的迭代循环 ---
            for iter_ktl_idx = 1:maxIter_ktl_loop
                % --- 4a. 计算当前K_guess下的市场价格 ---
                [R_mkt_gross_factor_iter_val, MPL_gross_iter_val] = main_olg_v7_utils.HHPrices_Huggett(K_current_guess, L_per_capita_global, cS_global);
                r_mkt_gross_iter_val = R_mkt_gross_factor_iter_val - 1; % 市场毛资本回报率 (MPK-delta)

                % --- 4b. 计算PAYG福利b_payg (基于目标rho'和当前工资) ---
                avg_worker_gross_wage_iter_val = 0;
                if mass_workers_global > 1e-9 && L_per_capita_global > 0 && MPL_gross_iter_val > 0
                    avg_worker_gross_wage_iter_val = (MPL_gross_iter_val * L_per_capita_global) / mass_workers_global;
                end
                b_payg_iter_val = rho_prime_payg_target_input * avg_worker_gross_wage_iter_val;
                b_payg_iter_val = max(0, b_payg_iter_val); % 确保非负

                % --- 4c. 计算实际PAYG税率theta_actual (考虑上限和总税负约束) ---
                theta_payg_actual_iter_val = theta_payg_required_calc; % 从理论值开始
                % 约束1: PAYG税率本身不能超过cS_global.theta_payg_max (这一步在开始已检查过，这里是针对调整后的tau_l_guess)
                % 约束2: 总劳动税负 (theta + tau_l) 不能超过 cS_global.max_total_labor_tax
                if (theta_payg_actual_iter_val + tau_l_current_guess) > cS_global.max_total_labor_tax
                    theta_payg_actual_iter_val = max(0, cS_global.max_total_labor_tax - tau_l_current_guess);
                end
                 % 约束1再次应用，以防tau_l过高导致theta被压到低于0后，又被theta_max限制
                theta_payg_actual_iter_val = min(theta_payg_actual_iter_val, cS_global.theta_payg_max);
                theta_payg_actual_iter_val = max(0, theta_payg_actual_iter_val); % 确保theta非负

                % --- 4d. 计算家庭面临的税后回报率 ---
                r_k_net_hh_iter_val = r_mkt_gross_iter_val * (1 - cS_global.tau_k); % 家庭税后资本净回报率
                R_k_net_hh_factor_iter_val = 1 + r_k_net_hh_iter_val; % 因子

                % --- 4e. 构建PAYG福利向量 (按年龄组) ---
                bV_payg_vec_iter_val = zeros(1, cS_global.aD_new);
                if cS_global.aR_new < cS_global.aD_new % 如果存在退休年龄组
                    bV_payg_vec_iter_val(cS_global.aR_new+1 : cS_global.aD_new) = b_payg_iter_val;
                end

                % --- 4f. 设置当轮VFI和模拟所需的家庭参数 ---
                paramS_for_vfi_sim_iter = paramS_global_in; % 复制全局派生参数
                paramS_for_vfi_sim_iter.tau_l = tau_l_current_guess; % 使用当前猜测的tau_l
                paramS_for_vfi_sim_iter.theta_payg_actual_for_hh = theta_payg_actual_iter_val; % 家庭面临的实际PAYG税率
                paramS_for_vfi_sim_iter.pps_tax_deferral_active = cS_global.pps_active; % PPS税收递延状态

                % --- 4g. 迭代求解遗赠 T_bequest (作为总转移支付 TR_total 的一部分) ---
                %    TR_gov 在此模型中设定为0，所以 TR_total = T_bequest
                TR_total_for_vfi_guess_val = 0.01 * MPL_gross_iter_val; % 遗赠的初始猜测
                if iter_ktl_idx > 1 && isfinite(T_bequest_model_iter_val) % 如果不是第一轮K-tau_l迭代，用上一轮结果
                    TR_total_for_vfi_guess_val = T_bequest_model_iter_val;
                end

                max_vfi_tr_sub_iter = 5;    % 遗赠子迭代最大次数
                tol_vfi_tr_sub_iter = 1e-3; % 遗赠子迭代收敛容忍度
                % 初始化VFI策略矩阵，以确保在子循环外有定义
                cPolM_4D_from_vfi_final = []; kPolM_4D_from_vfi_final = []; cPpsPolM_4D_rule_from_vfi_final = [];

                for i_vfi_tr_sub_loop = 1:max_vfi_tr_sub_iter
                    % 调用VFI求解家庭问题 (给定TR_total_for_vfi_guess_val)
                    [cPolM_vfi_temp, kPolM_vfi_temp, cPpsPolM_vfi_temp_rule, ~] = ...
                        main_olg_v7_utils.HHSolution_VFI_Huggett(R_k_net_hh_factor_iter_val, MPL_gross_iter_val, ...
                        TR_total_for_vfi_guess_val, bV_payg_vec_iter_val, paramS_for_vfi_sim_iter, cS_global);

                    % 基于VFI策略，模拟家庭路径以计算模型内生的遗赠
                    [kHistM_sim_for_bequest, ~, ~] = ...
                        main_olg_v7_utils.HHSimulation_olgm(kPolM_vfi_temp, cPpsPolM_vfi_temp_rule, cPolM_vfi_temp, ...
                        eIdxM_global_sim_paths, R_k_net_hh_factor_iter_val, MPL_gross_iter_val, ...
                        TR_total_for_vfi_guess_val, bV_payg_vec_iter_val, paramS_for_vfi_sim_iter, cS_global);

                    % 计算总意外遗赠 (人均)
                    % k' (下一期资产) 乘以回报因子，才是死亡时可用于遗赠的财富
                    kprime_for_bequest_calc = zeros(cS_global.nSim, cS_global.aD_orig);
                    if cS_global.aD_orig > 1
                        % kHistM是期末资产，成为下一期的期初资产。若在期末死亡，则此资产被遗赠。
                        kprime_for_bequest_calc(:, 1:cS_global.aD_orig) = kHistM_sim_for_bequest(:, 1:cS_global.aD_orig);
                    end
                    % 死亡时财富 = k'(t+1) * R_net (如果用期末资产) 或 k(t) * R_net (如果用期初资产且死亡发生在期末消费后)
                    % 这里假设遗赠的是期末资产 kHistM (即下一期的期初资产)
                    % 不需要再乘以 R_k_net_hh_factor，因为kHistM本身就是期末价值
                    
                    % 修正：遗赠发生在本期决策后，下一期开始前。所以遗赠的是k'。
                    % 而kHistM(:,a)是第a期期末（即第a+1期期初）的资产。
                    % 如果个体在第a期结束时死亡，他遗赠的是kHistM(:,a)。
                    % 死亡率 d_orig(a_orig) 是指在年度年龄 a_orig 这一年内死亡的概率。
                    % 假设死亡发生在年末，消费和储蓄决策已完成后。
                    
                    ageDeathMass_annual_iter_val = paramS_global_in.Z_ss_norm_annual(:) .* cS_global.d_orig(:); % 各年度年龄死亡人口占比
                    % 遗赠的财富是 kHistM，它是在 period t 结束时持有的资产，即 period t+1 的期初资产
                    % 如果个体在 period t (年度年龄 a_orig) 死亡，遗赠的是 kHistM(:, a_orig_loop_idx)
                    mean_bequest_wealth_per_age_iter_val = mean(kHistM_sim_for_bequest, 1); % 按年度年龄平均
                    TotalBequests_pc_iter_val = sum(mean_bequest_wealth_per_age_iter_val(:) .* ageDeathMass_annual_iter_val(:));
                    % 遗赠在下一期分配，需除以人口增长因子 (1+n)
                    T_bequest_model_new_iter_val = TotalBequests_pc_iter_val / (1 + paramS_global_in.popGrowthForDebt);
                    T_bequest_model_new_iter_val = max(0, T_bequest_model_new_iter_val); % 确保非负

                    % 更新最终的VFI策略 (以备子迭代不收敛时使用最后一轮结果)
                    T_bequest_model_iter_val = T_bequest_model_new_iter_val;
                    cPolM_4D_from_vfi_final = cPolM_vfi_temp;
                    kPolM_4D_from_vfi_final = kPolM_vfi_temp;
                    cPpsPolM_4D_rule_from_vfi_final = cPpsPolM_vfi_temp_rule;

                    % 检查遗赠子迭代是否收敛
                    if abs(T_bequest_model_new_iter_val - TR_total_for_vfi_guess_val) < tol_vfi_tr_sub_iter || ...
                       i_vfi_tr_sub_loop == max_vfi_tr_sub_iter
                        TR_total_for_vfi_final_iter = T_bequest_model_new_iter_val; % TR_gov=0
                        break; % 跳出遗赠子迭代
                    end
                    % 更新遗赠猜测值 (阻尼迭代)
                    TR_total_for_vfi_guess_val = 0.5 * TR_total_for_vfi_guess_val + 0.5 * T_bequest_model_new_iter_val;
                end % 结束遗赠子迭代
                T_bequest_model_iter_val = TR_total_for_vfi_final_iter; % 最终确定的当轮遗赠/总转移

                % --- 4h. 使用最终确定的VFI策略和TR_total进行最终的微观模拟 ---
                [kHistM_non_pps_sim_iter_val, kPpsHistM_sim_iter_val, cHistM_sim_iter_val] = main_olg_v7_utils.HHSimulation_olgm(...
                    kPolM_4D_from_vfi_final, cPpsPolM_4D_rule_from_vfi_final, cPolM_4D_from_vfi_final, eIdxM_global_sim_paths, ...
                    R_k_net_hh_factor_iter_val, MPL_gross_iter_val, TR_total_for_vfi_final_iter, bV_payg_vec_iter_val, paramS_for_vfi_sim_iter, cS_global);

                % --- 4i. 汇总模型内生的宏观量 ---
                % 非PPS资本总量 (人均)
                K_model_nonpps_sim_iter_val = mean(kHistM_non_pps_sim_iter_val, 1) * paramS_global_in.Z_ss_norm_annual;
                % PPS资本总量 (人均)
                K_model_pps_sim_iter_val = 0;
                if cS_global.pps_active && cS_global.pps_in_K && (cS_global.pps_max_contrib_frac > 0 || cS_global.pps_annual_contrib_limit > 0)
                    K_model_pps_sim_iter_val = mean(kPpsHistM_sim_iter_val, 1) * paramS_global_in.Z_ss_norm_annual;
                    K_model_pps_sim_iter_val = max(0, K_model_pps_sim_iter_val);
                end
                % 总生产性资本 (模型模拟结果)
                K_model_from_sim_iter_val = K_model_nonpps_sim_iter_val + K_model_pps_sim_iter_val;
                K_model_from_sim_iter_val = max(1e-6, K_model_from_sim_iter_val); % 避免为零
                % 总消费量 (模型模拟结果, 人均)
                C_model_iter_val = mean(cHistM_sim_iter_val,1) * paramS_global_in.Z_ss_norm_annual;

                % --- 4j. 计算政府一般预算残差 (GBC) ---
                %   注意：Y_for_gbc 使用的是 K_current_guess, 而非 K_model_from_sim
                %   这是因为我们要看在当前猜测的K下，tau_l是否能平衡预算
                Y_for_gbc_iter_val = cS_global.A * (K_current_guess^cS_global.alpha) * (L_per_capita_global^(1-cS_global.alpha));
                G_val_target_iter_val = cS_global.gov_exp_frac_Y * Y_for_gbc_iter_val; % 目标政府支出
                B_val_target_iter_val = cS_global.gov_debt_frac_Y * Y_for_gbc_iter_val; % 目标政府债务

                gbc_residual_iter_val = main_olg_v7_utils.check_gbc_residual(K_current_guess, C_model_iter_val, Y_for_gbc_iter_val, ...
                    G_val_target_iter_val, B_val_target_iter_val, MPL_gross_iter_val, r_mkt_gross_iter_val, ...
                    theta_payg_actual_iter_val, tau_l_current_guess, ...
                    b_payg_iter_val, T_bequest_model_iter_val, 0, cS_global, paramS_global_in); % TR_gov=0

                % --- 4k. 计算偏差和更新规则 ---
                K_dev_from_sim_iter_val = K_current_guess - K_model_from_sim_iter_val; % 资本偏差
                % tau_l偏差 (通过GBC残差对总工资收入的比例来近似，总工资是tau_l的主要税基)
                % 如果GBC>0 (盈余)，则tau_l应调低 (dev<0)
                % 如果GBC<0 (赤字)，则tau_l应调高 (dev>0)
                % 所以是 -gbc_res / tax_base
                tau_l_dev_raw_for_update = -gbc_residual_iter_val / (MPL_gross_iter_val * L_per_capita_global + 1e-9); % 加1e-9避免除零

                % 综合偏差范数 (V6风格，这里GBC残差本身即是目标)
                current_devNorm_val = sqrt(K_dev_from_sim_iter_val^2 + (gbc_residual_iter_val)^2 );
                norm_improvement_val = prev_devNorm_ktl - current_devNorm_val;

                % 打印当轮K-tau_l迭代信息
                fprintf('  %7d | %8.4f | %8.4f | %8.4f | %8.4f | %8.4f | %8.4f | %8.2e | %8.4f | %8.4f | %8.2e | %.1e | %7d\n', ...
                    iter_ktl_idx, K_current_guess, tau_l_current_guess, MPL_gross_iter_val, theta_payg_actual_iter_val, ...
                    K_model_from_sim_iter_val, K_model_pps_sim_iter_val, gbc_residual_iter_val, ...
                    K_dev_from_sim_iter_val, tau_l_dev_raw_for_update, current_devNorm_val, norm_improvement_val, tau_l_boundary_strike_count_ktl);

                % --- 4l. 检查收敛和可行性 ---
                % PAYG系统是否能被实际theta覆盖 (即实际theta不因总税负上限而被压低到无法满足理论需求)
                payg_fully_funded_by_actual_theta_check = (theta_payg_actual_iter_val >= theta_payg_required_calc - 1e-5);

                if current_devNorm_val < tol_norm_ktl_loop && ...
                   abs(gbc_residual_iter_val) < cS_global.gbc_tol_for_internal_loop && ...
                   payg_fully_funded_by_actual_theta_check
                    % 条件1: K和GBC的综合偏差范数达标
                    % 条件2: GBC残差本身也足够小
                    % 条件3: PAYG系统能被实际theta资助
                    converged_and_feasible_out = true;
                    K_sol_out = K_model_from_sim_iter_val; % 使用模拟得到的K作为解
                    tau_l_sol_out = tau_l_current_guess;   % 使用当前tau_l作为解
                    gbc_res_final_out = gbc_residual_iter_val;
                    % 存储详细解信息
                    solution_details_out.R_mkt_gross_factor = R_mkt_gross_factor_iter_val;
                    solution_details_out.MPL_gross = MPL_gross_iter_val;
                    solution_details_out.theta_payg = theta_payg_actual_iter_val; % 实际使用的theta
                    solution_details_out.b_payg = b_payg_iter_val;
                    solution_details_out.T_bequest_Model = T_bequest_model_iter_val;
                    solution_details_out.C_model = C_model_iter_val;
                    % Y_model 应基于最终的 K_sol_out 计算
                    solution_details_out.Y_model = cS_global.A * (K_sol_out^cS_global.alpha) * (L_per_capita_global^(1-cS_global.alpha));
                    solution_details_out.K_model_pps = K_model_pps_sim_iter_val;
                    solution_details_out.K_model_non_pps = K_model_nonpps_sim_iter_val;

                    fprintf('  solve_K_tau_l (V7): K和tau_l成功收敛 (rho_prime_target=%.4f, 实际theta_act=%.4f).\n', ...
                        rho_prime_payg_target_input, theta_payg_actual_iter_val);
                    break; % 跳出K-tau_l迭代循环

                elseif current_devNorm_val < tol_norm_ktl_loop && ...
                       abs(gbc_residual_iter_val) < cS_global.gbc_tol_for_internal_loop && ...
                       ~payg_fully_funded_by_actual_theta_check
                    % K和GBC已收敛，但PAYG系统因总税负上限而无法完全满足理论theta_req
                    fprintf('  solve_K_tau_l (V7): K, tau_l, GBC收敛 (rho_prime=%.4f), 但实际theta_payg (%.4f) 因总税负上限低于理论需求 (%.4f)。标记为不可行。\n', ...
                        rho_prime_payg_target_input, theta_payg_actual_iter_val, theta_payg_required_calc);
                    converged_and_feasible_out = false; % 标记为不可行
                    K_sol_out = K_model_from_sim_iter_val; tau_l_sol_out = tau_l_current_guess; gbc_res_final_out = gbc_residual_iter_val;
                    % 仍然存储解的细节，供外层判断
                    solution_details_out.R_mkt_gross_factor = R_mkt_gross_factor_iter_val;
                    solution_details_out.MPL_gross = MPL_gross_iter_val;
                    solution_details_out.theta_payg = theta_payg_actual_iter_val;
                    solution_details_out.b_payg = b_payg_iter_val;
                    solution_details_out.T_bequest_Model = T_bequest_model_iter_val;
                    solution_details_out.C_model = C_model_iter_val;
                    solution_details_out.Y_model = Y_for_gbc_iter_val; % 使用基于K_guess的Y
                    solution_details_out.K_model_pps = K_model_pps_sim_iter_val;
                    solution_details_out.K_model_non_pps = K_model_nonpps_sim_iter_val;
                    break; % 跳出K-tau_l迭代循环
                end % 结束收敛检查

                % --- 4m. 更新 K_guess 和 tau_l_guess 进行下一轮迭代 ---
                % 更新 K_guess
                K_guess_next_val = K_current_guess - dampK_ktl_loop * K_dev_from_sim_iter_val;
                K_current_guess = max(1e-3, K_guess_next_val); % 确保K为正

                % 更新 tau_l_guess
                new_tau_l_unconstrained_val = tau_l_current_guess + damp_tau_l_ktl_loop * tau_l_dev_raw_for_update; % 注意符号
                tau_l_next_iter_constrained_val = max(cS_global.tau_l_min, min(cS_global.tau_l_max, new_tau_l_unconstrained_val));

                % 检查 tau_l 是否持续触达边界且GBC未平衡
                is_tau_l_at_boundary_now = ...
                    (abs(tau_l_next_iter_constrained_val - cS_global.tau_l_max) < 1e-7 && new_tau_l_unconstrained_val >= cS_global.tau_l_max) || ...
                    (abs(tau_l_next_iter_constrained_val - cS_global.tau_l_min) < 1e-7 && new_tau_l_unconstrained_val <= cS_global.tau_l_min);

                if is_tau_l_at_boundary_now && abs(gbc_residual_iter_val) > cS_global.gbc_tol_for_internal_loop
                    tau_l_boundary_strike_count_ktl = tau_l_boundary_strike_count_ktl + 1;
                else
                    tau_l_boundary_strike_count_ktl = 0; % 重置触界计数
                end
                tau_l_current_guess = tau_l_next_iter_constrained_val; % 更新tau_l猜测值

                % 如果tau_l持续触界且GBC不平衡，则中止
                if tau_l_boundary_strike_count_ktl >= cS_global.max_tau_l_boundary_strikes && abs(gbc_residual_iter_val) > cS_global.gbc_tol_for_internal_loop
                    fprintf('  警告 (V7): tau_l 在边界 (%.4f) 持续 %d 次迭代，且GBC (%.2e) 未平衡。为 rho_prime=%.4f 中止。\n', ...
                        tau_l_current_guess, tau_l_boundary_strike_count_ktl, gbc_residual_iter_val, rho_prime_payg_target_input);
                    converged_and_feasible_out = false;
                    K_sol_out = K_model_from_sim_iter_val; tau_l_sol_out = tau_l_current_guess; gbc_res_final_out = gbc_residual_iter_val;
                    solution_details_out.R_mkt_gross_factor = R_mkt_gross_factor_iter_val; solution_details_out.MPL_gross = MPL_gross_iter_val;
                    solution_details_out.theta_payg = theta_payg_actual_iter_val; solution_details_out.b_payg = b_payg_iter_val;
                    solution_details_out.T_bequest_Model = T_bequest_model_iter_val; solution_details_out.C_model = C_model_iter_val;
                    solution_details_out.Y_model = Y_for_gbc_iter_val;
                    solution_details_out.K_model_pps = K_model_pps_sim_iter_val; solution_details_out.K_model_non_pps = K_model_nonpps_sim_iter_val;
                    break; % 跳出K-tau_l迭代
                end

                % --- 4n. 检查迭代是否停滞 ---
                if iter_ktl_idx > 1 % 从第二次迭代开始检查
                    if norm_improvement_val < (cS_global.min_norm_improvement_frac * prev_devNorm_ktl) && current_devNorm_val > tol_norm_ktl_loop
                        stagnation_counter_ktl = stagnation_counter_ktl + 1;
                    else
                        stagnation_counter_ktl = 0; % 重置停滞计数
                    end
                end
                prev_devNorm_ktl = current_devNorm_val; % 更新上一轮偏差范数

                if stagnation_counter_ktl >= cS_global.max_stagnation_iters && current_devNorm_val > tol_norm_ktl_loop
                    fprintf('  警告 (V7): 在 %d 次迭代后检测到范数停滞。范数: %.2e > 容忍度: %.1e。为 rho_prime=%.4f 中止。\n', ...
                        iter_ktl_idx, current_devNorm_val, tol_norm_ktl_loop, rho_prime_payg_target_input);
                    converged_and_feasible_out = false;
                    K_sol_out = K_model_from_sim_iter_val; tau_l_sol_out = tau_l_current_guess; gbc_res_final_out = gbc_residual_iter_val;
                    solution_details_out.R_mkt_gross_factor = R_mkt_gross_factor_iter_val; solution_details_out.MPL_gross = MPL_gross_iter_val;
                    solution_details_out.theta_payg = theta_payg_actual_iter_val; solution_details_out.b_payg = b_payg_iter_val;
                    solution_details_out.T_bequest_Model = T_bequest_model_iter_val; solution_details_out.C_model = C_model_iter_val;
                    solution_details_out.Y_model = Y_for_gbc_iter_val;
                    solution_details_out.K_model_pps = K_model_pps_sim_iter_val; solution_details_out.K_model_non_pps = K_model_nonpps_sim_iter_val;
                    break; % 跳出K-tau_l迭代
                end
            end % 结束 K-tau_l 迭代循环 (iter_ktl_idx)

            % --- 5. 处理循环结束后的情况 (例如达到最大迭代次数但未收敛) ---
            if ~converged_and_feasible_out && iter_ktl_idx == maxIter_ktl_loop
                fprintf('  警告 (V7): K和tau_l迭代达到最大次数 (%d) 或在该次数内未达可行解 (rho_prime_target=%.4f).\n', maxIter_ktl_loop, rho_prime_payg_target_input);
                % 使用最后一轮的值作为输出 (尽管可能不是均衡解)
                K_sol_out = K_model_from_sim_iter_val;
                tau_l_sol_out = tau_l_current_guess;
                gbc_res_final_out = gbc_residual_iter_val;

                % 确保 solution_details_out 有值，即使是最后一轮不完美的迭代结果
                if exist('MPL_gross_iter_val', 'var') && isfinite(MPL_gross_iter_val) % 检查是否有有效的价格信息
                    solution_details_out.R_mkt_gross_factor = R_mkt_gross_factor_iter_val;
                    solution_details_out.MPL_gross = MPL_gross_iter_val;
                    solution_details_out.theta_payg = theta_payg_actual_iter_val;
                    solution_details_out.b_payg = b_payg_iter_val;
                    solution_details_out.T_bequest_Model = T_bequest_model_iter_val;
                    solution_details_out.C_model = C_model_iter_val;
                    solution_details_out.Y_model = Y_for_gbc_iter_val; % 基于K_guess的Y
                    solution_details_out.K_model_pps = K_model_pps_sim_iter_val;
                    solution_details_out.K_model_non_pps = K_model_nonpps_sim_iter_val;
                else % 如果连价格信息都没有 (例如第一轮就失败了)
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
            end % 结束最大迭代次数处理

            % --- 6. 确保 solution_details_out 中的一些关键字段有值 ---
            % (主要针对theta_payg_required_before_cap，因为其他字段在循环中应已赋值)
            if ~isfield(solution_details_out, 'theta_payg_required_before_cap') || ...
               (isfield(solution_details_out, 'theta_payg_required_before_cap') && isnan(solution_details_out.theta_payg_required_before_cap))
                 % 如果之前没有存入，重新计算一次 (虽然理论上在循环开始时已存)
                 recalc_mass_retirees_final_sd = sum(paramS_global_in.ageMassV(cS_global.aR_new + 1 : cS_global.aD_new));
                 recalc_theta_req_final_sd = 0;
                 if mass_workers_global > 1e-9
                     recalc_theta_req_final_sd = rho_prime_payg_target_input * (recalc_mass_retirees_final_sd / mass_workers_global);
                 end
                 solution_details_out.theta_payg_required_before_cap = max(0, recalc_theta_req_final_sd);
            end

            % 确保PPS资本字段存在
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
            
            % 如果成功收敛，确保 Y_model 是基于最终的 K_sol_out 计算的
            if converged_and_feasible_out && isfield(solution_details_out, 'MPL_gross') && isfinite(solution_details_out.MPL_gross)
                % K_sol_out 已经是 K_model_from_sim_iter_val 了
                 solution_details_out.Y_model = cS_global.A * (K_sol_out^cS_global.alpha) * (L_per_capita_global^(1-cS_global.alpha));
            end

        end % 结束 solve_K_tau_l_for_rho_prime (V7)


        % --- check_gbc_residual (政府一般预算残差检查, V7版本) ---
        function gbc_residual_out = check_gbc_residual(...
                K_val_market_input, C_val_model_input, Y_val_market_input, G_val_target_input, B_val_target_input, ...
                MPL_gross_val_input, r_mkt_gross_val_input, ...
                theta_payg_val_actual_input, tau_l_val_input, ...
                b_payg_val_per_retiree_input, T_bequest_val_pc_input, TR_gov_val_pc_input, ...
                cS_check, paramS_loc_check)
            % check_gbc_residual (V7) - 计算政府一般预算(非PAYG部分)的残差
            %
            % 输入:
            %   K_val_market_input        - 市场资本存量 (用于计算税收)
            %   C_val_model_input         - 模型总消费 (用于计算消费税)
            %   Y_val_market_input        - 市场总产出 (用于确定G和B的目标水平)
            %   G_val_target_input        - 目标政府消费支出
            %   B_val_target_input        - 目标政府债务水平
            %   MPL_gross_val_input       - 市场毛工资率
            %   r_mkt_gross_val_input     - 市场毛资本回报率 (MPK-delta)
            %   theta_payg_val_actual_input - 实际PAYG税率 (不在此函数中使用，因PAYG自平衡)
            %   tau_l_val_input           - 一般劳动所得税率
            %   b_payg_val_per_retiree_input- PAYG福利 (不在此函数中使用)
            %   T_bequest_val_pc_input    - 人均意外遗赠 (不在此函数中使用，因为TR_gov已包含)
            %   TR_gov_val_pc_input       - 政府直接转移支付 (除PAYG和遗赠分配外的，TR_gov=0目标)
            %   cS_check                  - 公共参数 (tau_k, tau_c)
            %   paramS_loc_check          - 局部参数 (L_per_capita, popGrowthForDebt)
            %
            % 输出:
            %   gbc_residual_out          - 政府一般预算残差 (收入 - 支出)
            %
            % 说明:
            %   一般预算不包括PAYG系统的收支，假设PAYG系统自身是平衡的（或其赤字/盈余不影响此预算）。
            %   TR_gov_val_pc_input 设定为0，因为模型目标是tau_l调整以使TR_gov=0。
            %   此函数用于检查在给定的tau_l下，非PAYG预算是否平衡。

            L_per_capita_local_check = paramS_loc_check.L_per_capita; % 人均有效劳动

            % --- 计算政府一般收入 ---
            % 1. 一般劳动所得税收入 (基于tau_l)
            %    税基是总劳动收入 MPL_gross * L_pc。
            %    如果tau_l的税基在HHIncome中已扣除PPS，则这里也应使用调整后的税基。
            %    但通常GBC检查使用宏观总量。这里假设tau_l对总劳动收入征收（或平均意义上）。
            LaborTaxRev_general_part_calc = tau_l_val_input * MPL_gross_val_input * L_per_capita_local_check;
            % 2. 资本所得税收入
            CapitalTaxRev_calc = r_mkt_gross_val_input * K_val_market_input * cS_check.tau_k;
            % 3. 消费税收入
            ConsumptionTaxRev_calc = C_val_model_input * cS_check.tau_c;
            % 总一般收入
            GeneralRevenue_calc = LaborTaxRev_general_part_calc + CapitalTaxRev_calc + ConsumptionTaxRev_calc;

            % --- 计算政府一般支出 ---
            % 1. 政府消费
            GovConsumption_calc = G_val_target_input;
            % 2. 政府债务的净利息支出
            %    假设政府以税后资本回报率 r_b = r_mkt_gross * (1-tau_k) 或直接 r_mkt_gross 发债
            %    这里使用 r_mkt_gross_val_input 作为政府债券利率的代理 (与Heer(2020)一致性)
            r_b_for_debt_service_calc = r_mkt_gross_val_input; % 假设政府借贷利率等于市场毛回报率
            % 债务动态: B' = (1+r_b)B - PrimarySurplus。稳态时 B'=(1+n)B。
            % (1+n)B = (1+r_b)B - PrimarySurplus  => PrimarySurplus = (r_b - n)B
            % 我们需要的是 DebtService = r_b * B (或 (r_b-n)B 如果是净增长部分)
            % 原始代码v6用的是 (r_b - popGrowth) * B，这代表为维持B/Y不变所需的基本盈余中的利息部分。
            % 如果直接是利息支出，则是 r_b * B。
            % 采用 (r_b - popGrowth) * B 作为净债务服务成本，以与债务稳定动态一致。
            DebtService_calc = (r_b_for_debt_service_calc - paramS_loc_check.popGrowthForDebt) * B_val_target_input;
            % 3. 政府直接转移支付 (TR_gov, 在此模型中目标为0)
            GovDirectTransfers_calc = TR_gov_val_pc_input; % 应为0

            % 总一般支出
            GeneralOutlays_calc = GovConsumption_calc + DebtService_calc + GovDirectTransfers_calc;

            % --- 计算GBC残差 (收入 - 支出) ---
            gbc_residual_out = GeneralRevenue_calc - GeneralOutlays_calc;
        end

        % --- CallInterpolator (安全调用插值函数, V7版本) ---
        function ev_val_out = CallInterpolator(interpolant_obj_input, k_val_input, k_pps_val_input, cS_local_interp)
            % CallInterpolator (V7) - 安全地调用griddedInterpolant或函数句柄
            %
            % 输入:
            %   interpolant_obj_input - griddedInterpolant对象 或 函数句柄
            %   k_val_input           - 非PPS资产插值点
            %   k_pps_val_input       - PPS资产插值点
            %   cS_local_interp       - 局部参数结构体 (含网格信息 nk, nkpps, kGridV, kppsGridV)
            %
            % 输出:
            %   ev_val_out            - 插值得到的期望价值，出错则为-Inf或-1e12
            %
            % 说明:
            %   根据插值对象的类型和网格维度，以合适的方式调用插值。
            %   包含错误处理和边界值限制，以增强鲁棒性。

            ev_val_out = -Inf; % 默认返回值 (表示极差的价值)

            try
                if isa(interpolant_obj_input, 'griddedInterpolant')
                    % --- 处理 griddedInterpolant 对象 ---
                    if cS_local_interp.nk > 1 && cS_local_interp.nkpps > 1 % 2D插值 (k, k_pps)
                        ev_val_out = interpolant_obj_input(k_val_input, k_pps_val_input);
                    elseif cS_local_interp.nk > 1 % 1D插值 (k only, k_pps是标量或不使用)
                        ev_val_out = interpolant_obj_input(k_val_input);
                    elseif cS_local_interp.nkpps > 1 % 1D插值 (k_pps only, k是标量或不使用)
                        ev_val_out = interpolant_obj_input(k_pps_val_input);
                    else % 0D插值 (k和k_pps都是标量，插值器只有一个值)
                        if isscalar(interpolant_obj_input.Values)
                            ev_val_out = interpolant_obj_input.Values; % 直接取值
                        else
                            % 理论上不应到此，但作为备用尝试调用 (可能失败)
                            ev_val_out = interpolant_obj_input(k_val_input, k_pps_val_input);
                        end
                    end
                elseif isa(interpolant_obj_input, 'function_handle')
                    % --- 处理函数句柄 (假设它能处理输入) ---
                    % 通常用于0D情况，或当插值器本身就是个复杂函数时
                    ev_val_out = interpolant_obj_input(k_val_input, k_pps_val_input);
                else
                    warning('CallInterpolator: 未处理的插值器类型。');
                    ev_val_out = -1e11; % 返回一个表示错误的较大负值
                end
            catch ME_call_interp_error % 如果插值调用出错
                warning('CallInterpolator: 插值过程中发生错误: "%s"。将尝试限制输入值并重试。', ME_call_interp_error.message);
                % 限制输入值在网格范围内
                k_clamped = max(cS_local_interp.kGridV(1), min(cS_local_interp.kGridV(end), k_val_input));
                k_pps_clamped = 0; % 默认值
                if cS_local_interp.nkpps > 0 && ~isempty(cS_local_interp.kppsGridV) % 确保PPS网格有效
                    k_pps_clamped = max(cS_local_interp.kppsGridV(1), min(cS_local_interp.kppsGridV(end), k_pps_val_input));
                end

                % 再次尝试插值 (使用限制后的值)
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
                                ev_val_out = interpolant_obj_input(k_clamped, k_pps_clamped); % 备用
                            end
                        end
                    elseif isa(interpolant_obj_input, 'function_handle')
                        ev_val_out = interpolant_obj_input(k_clamped, k_pps_clamped);
                    end
                catch % 如果再次失败
                    ev_val_out = -1e11; % 最终的错误返回值
                end
            end % 结束主try-catch块

            % 确保最终输出是有限的
            if ~isfinite(ev_val_out)
                ev_val_out = -1e12; % 使用一个标准的非常差的价值
            end
        end % 结束 CallInterpolator

    end % 结束 Static 方法块
end % 结束 main_olg_v7_utils 类定义
% --- END OF FILE main_olg_v7_utils.m ---
