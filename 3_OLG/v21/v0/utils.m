classdef utils
    methods (Static)

        function cS = ParameterValues()
    % =========================================================================
            % == 函数: ParameterValues
            % == 版本: [v_het_types - 异质性家庭类型版]
            % ==
            % == 核心修改:
            % ==   - [!!!] 引入家庭类型维度 (cS.nTypes) 和各类型人口权重 (cS.type_weights)。
            % ==   - 将部分核心偏好参数 (beta, phi_bequest) 改造为向量，以便为不同类型赋值。
            % ==   - 移除原有的通用年龄效率计算，改为调用新的、专门的 age_efficiency 函数。
            % ==   - 改造对 calcaulte_theta_payg_path 的调用，以接收区分城乡的缴费率路径。
            % =========================================================================

          % --- [新增] 异质性家庭类型定义 ---
            cS.nTypes = 4; % 1:高收入城镇职工, 2:中低收入城镇职工, 3:高收入居民, 4:中低收入居民
            
            % 校准依据: 国家统计局2023年数据，高收入定义为收入最高的20%群体
            pop_share_urban = 0.7; % 假设城镇人口总占比为70%
            pop_share_resident = 1 - pop_share_urban; % 农村人口总占比为30%
            
            pop_share_high_income_in_urban = 0.2; % 城镇人口中，高收入群体占20%
            pop_share_low_income_in_urban = 1 - pop_share_high_income_in_urban; % 城镇中低收入群体占80%
            
            pop_share_high_income_in_resident = 0.2; % 农村人口中，高收入群体占20%
            pop_share_low_income_in_resident = 1 - pop_share_high_income_in_resident; % 农村中低收入群体占80%

            cS.type_weights = [...
                pop_share_urban * pop_share_high_income_in_urban; ...         % 类型1: 高收入城镇职工
                pop_share_urban * pop_share_low_income_in_urban; ...      % 类型2: 中低收入城镇职工
                pop_share_resident * pop_share_high_income_in_resident; ...    % 类型3: 高收入居民
                pop_share_resident * pop_share_low_income_in_resident ... % 类型4: 中低收入居民
                ];
            % 确保权重和为1
            cS.type_weights = cS.type_weights / sum(cS.type_weights);
            % --- 核心时间与年龄参数 ---
            cS.time_Step = 5;                   % 模型中每个时期代表的年数 (例如，5年为一个年龄组)。
            cS.age1_orig = 20;                  % 个体进入模型的起始物理年龄 (例如，20岁)。
            cS.ageLast_orig = 98;               % 个体确定性死亡前的最大物理年龄 (例如，98岁)。
            cS.ageRetire_orig = 65;             % 个体的法定退休物理年龄 (例如，60岁)。

            % --- 行为参数 (家庭偏好) ---
            cS.beta = 0.98;                    % 主观折现因子 (年化)，反映了个体的耐心程度。模型中将根据time_Step进行调整。
            cS.sigma = 3;                       % 跨期替代弹性的倒数 (相对风险厌恶系数)，控制消费平滑的意愿。
            cS.phi_bequest = 0.5;                 % 遗赠效用权重，衡量“暖心式”(warm-glow)遗赠动机的强度。
            cS.cFloor = 0.005;                  % 最低消费水平，用于保证效用函数在消费趋近于零时有界。

            % --- 生产技术参数 ---
            cS.alpha = 0.3;                    % 产出中私人资本的份额 (资本产出弹性)。
            cS.gamma = 0.12;                    % 产出中公共资本的份额 (公共资本产出弹性)。
            cS.ddk = 1 - (1 - 0.03)^cS.time_Step; % 私人资本的折旧率 (模型期)，由年化率0.015转换而来。
            cS.ddk_g = cS.ddk;                  % 公共资本的折旧率 (模型期)，为简化设为与私人资本相同。
            cS.A = 1.0;                         % 全要素生产率(TFP)的初始基准水平，通常归一化为1。

            % --- 政府与政策参数 ---
            cS.tau_k = 0.05;                    % 资本利得税率。
            cS.tau_l = 0.06;                    % 劳动收入税率 (在PAYG缴费之外)。
            cS.tau_c = 0.03;                    % 消费税率。
            cS.G_c_to_Y_ratio_ss = 0.05;        % [新增] 稳态下，政府消费占GDP的目标比率。
            cS.I_g_to_Y_ratio_ss = 0.03;
            % --- 步骤 2.1b: 设定稳态下的目标政府债务/GDP比率 ---

                        % -----PAYG设定------
            cS.n_contrib_years = 35; % 平均缴费期
            cS.personal_account_factor_urban = 0.25; % 城镇职工PAYG中个人账户比例
            cS.resident_basic_benefit_ratio = 0.18; % 城乡居民能领到的基础养老金，社平工资的比例
            cS.resident_personal_factor = 0.05; % 城乡居民PAYG中个人账户比例
            cS.K_payg_to_Y_ratio_ss0 = 0.055; % 2023社保基金余额存量占GDP比重

            % --- 长期BGP(平衡增长路径)参数 ---
            cS.g_A_ss = 0.012;                  % 长期稳态下技术进步的年化增长率 (TFP growth rate)。
            cS.long_term_birth_rate_target_per_1000 = 3.5; % 长期稳态目标的年化粗出生率 (每千人)。
            cS.transition_period_years = 10;    % 外生路径从数据点过渡到长期稳态所需的年数。
            cS.curvature_param = 0.5;           % 过渡路径的曲率参数 (k<1为上凸收敛，k>1为下凸收敛)。
            cS.pop_data_last_year = 2050;       % 出生率和死亡率数据都将在此年份后保持恒定

            % --- 养老金体系参数 ---
            cS.pps_active = false;              % 控制是否激活私人养老金(PPS)模块的全局开关。
            cS.pps_contrib_rate = 0.03;         % PPS的强制缴费率 (占劳动收入的比例)。
            cS.pps_withdrawal_rate = 0.05;      % 退休后每个模型期强制提取PPS资产的比例。
            cS.pps_tax_rate_withdrawal = 0.03;  % 提取PPS资产时适用的税率。
            cS.theta_path = 0.06;               % 外生设定的PAYG养老金缴费率 (在DC模式下的基准值)。
            cS.pps_max = 0.1; % 最高缴费比例
            cS.pps_fixed = 0.1;% 固定pps缴费比例
            
            % 设定模拟的起始年份
            cS.ss0_year = 2023;
            cS.start_year = 2023;
            cS.end_year = 2100;

            % --- 收入过程参数 ---
            cS.nw = 1;   % 5                       % 劳动效率冲击离散状态的数量。
            % 年龄节点
            cS.young_shock_start_age = 26;
            cS.young_shock_end_age = 35;
            cS.old_shock_start_age = 65;
            cS.old_shock_end_age = 88;
            % 冲击年化概率峰值
            cS.p_shock_young_peak_annual = 0; %0.01;
            cS.p_shock_old_peak_annual = 0; %0.01;
            % 冲击支出比例 (kappa)
            cS.kappa_young = 0.9;
            cS.kappa_old = 0.9;

            cS = population.generate_mortality_path(cS); % 生成死亡率数据

            % [修改] 调用新的年龄效率函数
            cS = utils.age_efficiency(cS);
        end

        function cS = age_efficiency(cS)
            % =========================================================================
            % == 函数: age_efficiency
            % == 版本: [v1.1 - 基于2023年数据校准版]
            % ==
            % == 目的:
            % ==   - 根据2023年官方数据校准不同类型家庭的相对收入水平。
            % ==   - 高收入组(前20%)的平均收入约为全体平均的2.2倍，
            % ==     中低收入组(后80%)的平均收入约为全体平均的0.725倍。
            % =========================================================================

            % --- 1. 计算基准的物理年龄效率剖面 (与之前相同) ---
            ages_eff = (cS.age1_orig : cS.ageLast_orig)';
            beta0=-13.215788; beta1=1.349514; beta2=-0.043363; beta3=0.000585; beta4=-0.000003;
            log_age_eff_orig = beta0 + beta1*ages_eff + beta2*(ages_eff.^2) + beta3*(ages_eff.^3) + beta4*(ages_eff.^4);
            ageEffV_orig_unnormalized = exp(log_age_eff_orig);
            ageEffV_orig_unnormalized((cS.aR_idx_orig):end) = 0;

            mean_efficiency_working = mean(ageEffV_orig_unnormalized(1:(cS.aR_idx_orig - 1)));
            ageEffV_orig_normalized = ageEffV_orig_unnormalized / mean_efficiency_working;

            % --- 2. [核心校准] 定义各类型的效率缩放因子 ---
            % 依据: 高收入组(20%)收入是平均的2.2倍，中低收入组(80%)是平均的0.725倍
            % 推导: (2.2 * 0.2) + (X * 0.8) = 1  => 0.44 + 0.8X = 1 => 0.8X = 0.56 => X = 0.7
            % 为校准更精确，使用0.725，使得加权平均更接近1
            
            scale_high_income = 2.2; %;
            scale_low_income = 0.725; %; 

            type_scales = [
                scale_high_income;  % 1: 高收入城镇职工
                scale_low_income;   % 2: 中低收入城镇职工
                scale_high_income;  % 3: 高收入居民
                scale_low_income    % 4: 中低收入居民
                ];

            % --- 3. 生成特定于类型的物理年龄效率剖面 (与之前相同) ---
            ageEff_by_type_orig = ageEffV_orig_normalized .* type_scales';

            % --- 4. 将物理年龄剖面聚合到模型年龄组 (与之前相同) ---
            cS.ageEffV_new_h = zeros(cS.aD_new,cS.nTypes);
            for i_type = 1:cS.nTypes
                ageEff_type_i_orig = ageEff_by_type_orig(:, i_type);
                ageEff_type_i_new = zeros(cS.aD_new, 1);
                for a = 1:cS.aD_new
                    phys_ages_in_group = cS.physAgeMap{a};
                    model_ages_indices = phys_ages_in_group - cS.age1_orig + 1;
                    ageEff_type_i_new(a) = mean(ageEff_type_i_orig(model_ages_indices));
                end
                cS.ageEffV_new_h(:,i_type) = ageEff_type_i_new;
            end

            fprintf('✅ 已生成 %d 类家庭的年龄-效率剖面 (基于2023年数据校准)。\n', cS.nTypes);
        end


        function cS = generateGrids(cS, varargin)
            % [BGP修改] 支持可选的网格上限参数，实现自适应网格
            % [BGP修改] 修改默认值计算逻辑，适应稳态化模型
            % 新签名: generateGrids(cS, 'k_max', value, 'kpps_max', value)
            % 这将网格生成从固定参数变为动态服务

            % 使用 inputParser 处理可选参数
            p = inputParser;

            % [BGP修改] 简化默认值计算逻辑
            % 在稳态化求解中，这些默认值几乎不会被用到，因为求解器会根据猜测的
            % 标准化资本 k̂_guess 来动态地传入新的 k_max
            default_kMax = 10; % 直接设为合理的常数
            default_kppsMax = 2; % 直接设为合理的常数

            % 添加可选参数
            addParameter(p, 'k_max', default_kMax, @isnumeric);
            addParameter(p, 'kpps_max', default_kppsMax, @isnumeric);
            parse(p, varargin{:});

            % 获取最终的网格上限值
            final_k_max = p.Results.k_max;
            final_kpps_max = p.Results.kpps_max;

            % 设置网格参数
            cS.kMin = 0;
            cS.kMax = final_k_max;
            cS.kppsMin = 0;
            cS.kppsMax = final_kpps_max;

            % 使用非线性网格以更好地捕捉财富分布
            power_k = 5;

            % 生成常规资产网格
            if cS.nk > 1
                kGridV_temp = cS.kMin + (cS.kMax - cS.kMin) * (linspace(0, 1, cS.nk).^power_k);
            else
                kGridV_temp = cS.kMin;
            end
            cS.kGridV = kGridV_temp(:);

            % 生成PPS资产网格
                if cS.nkpps > 1
                    kppsGridV_temp = cS.kppsMin + (cS.kppsMax - cS.kppsMin) * (linspace(0, 1, cS.nkpps).^power_k);
                else
                    kppsGridV_temp = [cS.kppsMin];
                end
                cS.kppsGridV = kppsGridV_temp(:);
        end

        % =========================================================================
        % ==               最终修正版的 EarningProcess_AgeDependent
        % =========================================================================
        function [leGridV_expanded, TrProbM_by_age, leProb1V, nw_expanded] = EarningProcess_AgeDependent(cS)
            % =========================================================================
            % == 函数: EarningProcess_AgeDependent
            % == 版本: [v12 - 最终转置修正版]
            % ==
            % == 目的:
            % ==   纠正之前对 tauchen 函数输出的错误假设。
            % ==
            % == 核心修正:
            % ==   - [!!! 关键 !!!] 移除了对 tauchen 输出矩阵的转置操作。
            % ==     经过验证，tauchen 函数本身就生成了一个标准的“行随机”矩阵，
            % ==     额外的转置是错误的，并导致了行和不为1的警告。
            % ==     此修正确保了从源头开始就使用正确的转移概率矩阵。
            % =========================================================================

            % --- 1. 生成基础AR(1)过程 ---
            lePersistence = 0.9; leShockStd = 0.15^0.5; Tauchen_q = 2.0;
            [leLogGridV_normal, leTrProbM_from_tauchen] = utils.tauchen(cS.nw, lePersistence, leShockStd, 0, Tauchen_q);

            % [!!! 根本性修正 !!!]
            % tauchen函数返回的是一个标准的行随机矩阵 (行和为1)。
            % 无需转置。
            leTrProbM_normal = leTrProbM_from_tauchen;

            % 检查修正是否有效（可选但推荐）
            row_sums = sum(leTrProbM_normal, 2);
            if any(abs(row_sums - 1.0) > 1e-6)
                % 如果这里仍然报警，说明 tauchen 函数的内部实现有问题。
                warning('EarningProcess_AgeDependent:RowSumCheck', ...
                    '转移矩阵 leTrProbM_normal 的行和不为1，请检查 tauchen 函数的内部实现！');
            end

            % 现在基于正确的、行随机的 leTrProbM_normal 计算平稳分布
            % 注意：为了找到对应于特征值1的左特征向量（即平稳分布），我们需要转置矩阵
            [V, D] = eig(leTrProbM_normal');
            [~, c] = min(abs(diag(D)-1));
            % 确保 V(:,c) 是平稳分布向量
            leProb1V = abs(V(:,c)/sum(V(:,c))); % 新生儿的长期冲击分布

            % --- 2. [自动切换] 检测是否存在重大支出冲击 ---
            has_major_shocks = (isfield(cS, 'p_shock_young_peak_annual') && cS.p_shock_young_peak_annual > 0) || ...
                (isfield(cS, 'p_shock_old_peak_annual') && cS.p_shock_old_peak_annual > 0);

            if has_major_shocks
                % --- 路径 A: 生成包含重大支出冲击的扩展过程 ---
                % (此部分代码无需修改，因为它现在会接收到正确的行随机 leTrProbM_normal)
                fprintf('   正在生成【含】重大支出冲击的收入过程...\n');
                phys_ages = cS.age1_orig : cS.ageLast_orig;
                num_phys_ages = length(phys_ages);
                p_young_annual = zeros(num_phys_ages, 1);
                p_old_annual = zeros(num_phys_ages, 1);
                y_start_idx = find(phys_ages == cS.young_shock_start_age, 1);
                y_end_idx = find(phys_ages == cS.young_shock_end_age, 1);
                if ~isempty(y_start_idx) && ~isempty(y_end_idx)
                    p_young_annual(y_start_idx:y_end_idx) = cS.p_shock_young_peak_annual;
                end
                o_start_idx = find(phys_ages == cS.old_shock_start_age, 1);
                o_end_idx = find(phys_ages == cS.old_shock_end_age, 1);
                if ~isempty(o_start_idx) && ~isempty(o_end_idx)
                    p_old_annual(o_start_idx:o_end_idx) = cS.p_shock_old_peak_annual;
                end
                p_shock_path_model = zeros(cS.aD_new, 2);
                for a = 1:cS.aD_new
                    phys_age_indices_in_group = cS.physAgeMap{a};
                    p_no_shock_y = 1.0;
                    p_no_shock_o = 1.0;
                    for phys_idx = phys_age_indices_in_group
                        if phys_idx > 0 && phys_idx <= num_phys_ages
                            p_no_shock_y = p_no_shock_y * (1 - p_young_annual(phys_idx));
                            p_no_shock_o = p_no_shock_o * (1 - p_old_annual(phys_idx));
                        end
                    end
                    p_shock_path_model(a, 1) = 1 - p_no_shock_y;
                    p_shock_path_model(a, 2) = 1 - p_no_shock_o;
                end
                leGridV_expanded = [exp(leLogGridV_normal(:)); 0; 0];
                nw_expanded = cS.nw + 2;
                TrProbM_by_age = cell(cS.aD_new, 1);
                for a = 1:cS.aD_new
                    Tr_a = zeros(nw_expanded, nw_expanded);
                    p_y = p_shock_path_model(a, 1);
                    p_o = p_shock_path_model(a, 2);
                    p_total_shock = min(1.0, p_y + p_o);
                    p_y_norm = 0;
                    p_o_norm = 0;
                    if p_total_shock > 0
                        p_y_norm = p_y / p_total_shock;
                        p_o_norm = p_o / p_total_shock;
                    end
                    Tr_a(1:cS.nw, 1:cS.nw) = leTrProbM_normal * (1 - p_total_shock);
                    Tr_a(1:cS.nw, cS.nw + 1) = p_total_shock * p_y_norm;
                    Tr_a(1:cS.nw, cS.nw + 2) = p_total_shock * p_o_norm;
                    Tr_a(cS.nw + 1, 1:cS.nw) = leProb1V';
                    Tr_a(cS.nw + 2, 1:cS.nw) = leProb1V';
                    row_sums_a = sum(Tr_a, 2);
                    row_sums_a(row_sums_a == 0) = 1; % 防止除以零
                    Tr_a = Tr_a ./ row_sums_a;
                    TrProbM_by_age{a} = Tr_a;
                end
            else
                % --- 路径 B: 生成不含重大支出冲击的标准过程 ---
                fprintf('   检测到重大冲击参数为零，正在生成【无】重大支出冲击的标准收入过程...\n');
                nw_expanded = cS.nw;
                leGridV_expanded = exp(leLogGridV_normal(:));
                TrProbM_by_age = cell(cS.aD_new, 1);
                for a = 1:cS.aD_new
                    % 现在 leTrProbM_normal 已经是正确的行随机矩阵
                    TrProbM_by_age{a} = leTrProbM_normal;
                end
            end
        end

        function [y_grid_out, trProbM_out] = tauchen(N, rho, sigma, mu, m)
            % --- [新增] 处理 N=1 (同质代理人) 的边缘情况 ---
            if N == 1
                y_grid_out = 0;      % 唯一的网格点是对数均值0
                trProbM_out = 1;     % 转移概率为100%
                return;              % 直接返回，跳过后续计算
            end

            % [BGP不变] 标准的、用于离散化AR(1)过程的工具函数
            std_y = sqrt(sigma^2 / (1-rho^2)); y_max = m*std_y; y_min = -y_max;
            y = linspace(y_min, y_max, N); d = y(2)-y(1);
            trProbM_out = zeros(N,N);
            for j=1:N, for k=1:N
                    m_k = rho*y(j) + mu;
                    if k==1, trProbM_out(j,k) = normcdf((y(1)-m_k+d/2)/sigma);
                    elseif k==N, trProbM_out(j,k) = 1 - normcdf((y(N)-m_k-d/2)/sigma);
                    else, trProbM_out(j,k) = normcdf((y(k)-m_k+d/2)/sigma) - normcdf((y(k)-m_k-d/2)/sigma); end
            end, end
        y_grid_out = y(:);
        end

        function cS = calcaulte_theta_payg_path(cS, graph_flag)
            % =========================================================================
            % == 函数: calcaulte_theta_payg_path
            % == 版本: [v3.0 - 城乡分离路径版]
            % ==
            % == 核心修改:
            % ==   - 函数的最终输出从单一混合路径 cS.theta_path 修改为两条独立路径：
            % ==     cS.theta_path_urban 和 cS.theta_path_resident。
            % ==   - 这使得VFI可以根据家庭类型（城镇或居民）选择对应的缴费率。
            % ==   - 移除了原有的混合计算逻辑。
            % =========================================================================
            fprintf('正在构建【分离的城乡】养老金缴费率路径...\n');

            % --- 1. 定义关键参数的【起点】与【终点】 ---
            theta_urban_initial = 0.2;   % 初始值: 基于2023年 16%+8% 的法定费率
            theta_urban_final   = 0.2;  % 长期目标: 参考澳大利亚等国的成熟DC体系费率

            theta_resident_initial = 0.05; % 初始值: 对当前名义缴费的合理有效费率估算
            theta_resident_final   = 0.05; % 长期目标: 假设保基本的作用略有增强

            year_reach_final = 2050;
            fprintf('   [政策设定] 缴费率将从 %.2f%%(U)/%.1f%%(R) 在 %d 年线性过渡到 %.2f%%(U)/%.1f%%(R)\n', ...
                theta_urban_initial*100, theta_resident_initial*100, year_reach_final, theta_urban_final*100, theta_resident_final*100);

            % --- 2. 生成缴费率本身的年度路径 ---
            % (借用一个辅助函数来生成平滑过渡路径，实际应在utils中定义)
            smooth_transition = @(start_val, end_val, steps, ~) linspace(start_val, end_val, steps);

            theta_urban_path_annual_trans = smooth_transition(theta_urban_initial, theta_urban_final, year_reach_final - cS.start_year + 1, 1);
            theta_resident_path_annual_trans = smooth_transition(theta_resident_initial, theta_resident_final, year_reach_final - cS.start_year + 1, 1);

            annual_years_vector = cS.start_year:cS.end_year;
            full_theta_urban_path = ones(size(annual_years_vector)) * theta_urban_final;
            full_theta_resident_path = ones(size(annual_years_vector)) * theta_resident_final;

            trans_len = length(theta_urban_path_annual_trans);
            path_indices = 1:min(trans_len, length(full_theta_urban_path));
            full_theta_urban_path(path_indices) = theta_urban_path_annual_trans(path_indices);
            full_theta_resident_path(path_indices) = theta_resident_path_annual_trans(path_indices);

            % --- 3. 从年度路径中，提取出模型【期】所需的路径 ---
            model_year_indices = 1:(cS.end_year - cS.start_year + 1);
            model_year_indices_5yr = model_year_indices(1:cS.time_Step:end);

            T_sim = length(cS.start_year:cS.time_Step:cS.end_year);

            % [核心修改] 提取并存储两条独立的路径
            theta_path_urban = full_theta_urban_path(model_year_indices_5yr);
            theta_path_urban = theta_path_urban(1:min(T_sim, length(theta_path_urban)));

            theta_path_resident = full_theta_resident_path(model_year_indices_5yr);
            theta_path_resident = theta_path_resident(1:min(T_sim, length(theta_path_resident)));
            
            cS.theta_path_h = [theta_path_urban;
                theta_path_urban;...
                theta_path_resident;
                theta_path_resident];

            % --- 4. (可选) 可视化检查路径 ---
            if graph_flag
                T_plot = length(cS.theta_path_urban);
                time_axis = cS.start_year:cS.time_Step:(cS.start_year + cS.time_Step*(T_plot-1));
                figure('Name', 'Separated PAYG Contribution Rate Paths');
                plot(time_axis, cS.theta_path_urban, 'b-s', 'LineWidth', 2, 'MarkerFaceColor', 'b', 'DisplayName', '城镇职工缴费率 (θ_t^{urban})');
                hold on;
                plot(time_axis, cS.theta_path_resident, 'r-d', 'LineWidth', 2, 'MarkerFaceColor', 'r', 'DisplayName', '城乡居民缴费率 (θ_t^{resident})');
                title('模型PAYG缴费率路径 (城乡分离)');
                xlabel('年份');
                ylabel('缴费率 (θ_t)');
                legend('show', 'Location', 'best');
                grid on;
                ylim([0, max(cS.theta_path_urban)*1.2]);
            end
            fprintf('✅ 分离的城乡缴费率路径构建完成。\n');
        end


        % function [muM, utilM] = CES_utility(cM, sigma, cS)
        %     % [BGP修改] 此函数的作用是输入一个消费值，输出其效用
        %     % [BGP修改] 在标准化框架下，VFI求解器会向它传递标准化的消费 ĉ
        %     % [BGP修改] 这个函数会正确地计算 u(ĉ) = ĉ^(1-σ)/(1-σ)
        %     % [BGP修改] 无需知道绝对技术水平A，完全在标准化环境中操作
        %     if cM<=0
        %         c_adj = cS.cFloor;
        %         print('消费小于零')
        %     else
        %         c_adj = cM;
        %     end
        %     if abs(sigma - 1) < 1e-6
        %         utilM = log(c_adj);
        %         muM = 1./c_adj;
        %     else utilM = (c_adj.^(1-sigma))./(1-sigma);
        %         muM = c_adj.^(-sigma);
        %     end
        % end
        % 
        % function util_beq = bequest_utility(k_prime, cS)
        %     % [BGP修改] 此函数将被传递标准化的遗赠 k̂'，并正确计算其效用
        %     % [BGP修改] 在标准化框架下，遗赠效用函数处理标准化的遗赠值
        %     % [BGP修改] 无需知道绝对技术水平A，完全在标准化环境中操作
        %     if cS.phi_bequest <= 0
        %         util_beq = 0;
        %         return;
        %     end
        % 
        %     k_adj = max(cS.cFloor, k_prime);
        % 
        %     if abs(cS.sigma - 1) < 1e-6
        %         util_beq = cS.phi_bequest * log(k_adj);
        %     else
        %         util_beq = cS.phi_bequest * (k_adj.^(1-cS.sigma))./(1-cS.sigma);
        %     end
        % end
        % 

                function [muM, utilM] = CES_utility(cM, sigma, cS)
            % =========================================================================
            % == 函数: CES_utility (v2.0 - 向量化安全修正版)
            % == 核心修正:
            % ==   - [!!! 关键稳健性修正 !!!] 替换了非向量化的 `if cM<=0`
            % ==     判断。旧代码在处理向量输入时，仅判断第一个元素，
            % ==     当向量中包含非正值时会导致 NaN 或 Inf。
            % ==   - 新方法: 使用 `max(cM, cS.cFloor)` 对消费输入进行
            % ==     向量化处理，确保所有元素都大于等于消费下限，然后再
            % ==     代入效用函数计算。这保证了函数对任意向量输入的健壮性。
            % =========================================================================
            
            % 使用 max 函数进行向量化安全处理，确保消费不低于下限
            c_adj = max(cM, cS.cFloor);

            if abs(sigma - 1) < 1e-6
                utilM = log(c_adj);
                muM = 1./c_adj;
            else
                utilM = (c_adj.^(1-sigma))./(1-sigma);
                muM = c_adj.^(-sigma);
            end
        end

        function util_beq = bequest_utility(k_prime, cS)
            % =========================================================================
            % == 函数: bequest_utility (v2.0 - 向量化安全修正版)
            % == 核心修正:
            % ==   - [!!! 关键稳健性修正 !!!] 与CES_utility的修正类似，
            % ==     使用 `max(k_prime, cS.cFloor)` 对遗赠财富输入进行
            % ==     向量化处理，以防止非正值输入导致计算错误。
            % ==   - 这确保了函数在任何求解器中（VFI, EGM）的行为都
            % ==     是稳健和可预测的。
            % =========================================================================
            if cS.phi_bequest <= 0
                util_beq = zeros(size(k_prime)); % 确保输出维度正确
                return;
            end

            % 使用 max 函数进行向量化安全处理
            k_adj = max(k_prime, cS.cFloor);

            if abs(cS.sigma - 1) < 1e-6
                util_beq = cS.phi_bequest * log(k_adj);
            else
                util_beq = cS.phi_bequest * (k_adj.^(1-cS.sigma))./(1-cS.sigma);
            end
        end


        function A_path = generate_tfp_path(cS,graph_flag)
            % =========================================================================
            % == 函数: generate_tfp_path (v_curvature_control - 曲率可控最终版)
            % ==
            % == 核心修正:
            % ==   调用新的smooth_transition函数，并传入一个曲率参数，
            % ==   确保TFP增长率路径是平滑收敛的。
            % =========================================================================

            fprintf('   正在构建TFP路径 (v_curvature_control - 曲率可控最终版)...\n');
            model_sim_years = cS.start_year:cS.time_Step:cS.end_year;
            annual_years_vec = cS.start_year:cS.end_year;

            % --- 1. 加载数据和设定参数 ---
            tfp_data_pwt = readtable('data\PWT\china_pwt_data.xlsx');
            pwt_years = tfp_data_pwt.year(tfp_data_pwt.year >= 1979);
            pwt_tfp_level = tfp_data_pwt.rtfpna(tfp_data_pwt.year >= 1979);
            bai_projections = [2025, 0.0557; 2030, 0.0482; 2035, 0.0394; 2040, 0.0340; 2045, 0.0346; 2050, 0.0298];

            long_term_g_annual = cS.g_A_ss;
            fprintf('   TFP长期年化增长率将平滑收敛至: %.4f%% (曲率k=%.1f)\n', long_term_g_annual * 100, cS.curvature_param);

            data_end_year = bai_projections(end, 1);
            transition_end_year = data_end_year + cS.transition_period_years;

            % --- 2. 构建平滑的【年度增长率】路径 ---
            g_path_annual = zeros(size(annual_years_vec));
            g_anchor_years = [pwt_years(end); bai_projections(:,1)];
            g_anchor_rates = [0.058; bai_projections(:,2)];

            data_mask = annual_years_vec <= data_end_year;
            g_path_annual(data_mask) = interp1(g_anchor_years, g_anchor_rates, annual_years_vec(data_mask), 'pchip', 'extrap');
            last_data_g = g_path_annual(find(data_mask, 1, 'last'));

            transition_mask = annual_years_vec > data_end_year & annual_years_vec <= transition_end_year;
            if any(transition_mask)
                num_trans_years = sum(transition_mask);
                g_path_annual(transition_mask) = utils.smooth_transition(last_data_g, long_term_g_annual, num_trans_years, cS.curvature_param);
            end

            ss_mask = annual_years_vec > transition_end_year;
            g_path_annual(ss_mask) = long_term_g_annual;

            % --- 3. 从增长率路径构建水平路径 ---
            A_path_annual_level = zeros(size(annual_years_vec));

            % 找到起始年份的TFP水平作为锚点
            start_year_tfp_level_anchor = interp1([pwt_years; 2100], [pwt_tfp_level; pwt_tfp_level(end)*(1.02)^50], cS.start_year, 'pchip', 'extrap');

            start_year_idx = find(annual_years_vec == cS.start_year, 1);
            A_path_annual_level(start_year_idx) = start_year_tfp_level_anchor;

            % 向前积分
            for t_idx = (start_year_idx + 1):length(annual_years_vec)
                A_path_annual_level(t_idx) = A_path_annual_level(t_idx - 1) * (1 + g_path_annual(t_idx));
            end
            % 向后积分
            for t_idx = (start_year_idx - 1):-1:1
                A_path_annual_level(t_idx) = A_path_annual_level(t_idx + 1) / (1 + g_path_annual(t_idx + 1));
            end

            % --- 4. 采样并归一化 ---
            A_path = interp1(annual_years_vec, A_path_annual_level, model_sim_years, 'linear');
            if A_path(1) ~= 0 && ~isnan(A_path(1))
                A_path = A_path / A_path(1);
            else
                warning('初始TFP水平为0或NaN，无法归一化。');
            end
            fprintf('✅ TFP路径构建完成。\n');

            if graph_flag
                annual_years_vec = cS.start_year:cS.end_year;
                if length(model_sim_years) >= 2
                    A_path_annual_for_plot = interp1(model_sim_years, A_path, annual_years_vec, 'pchip');
                else
                    A_path_annual_for_plot = A_path;
                end
                utils.create_A_dynamics_plots(cS, A_path_annual_for_plot,...
                    annual_years_vec);
            end

        end

        function create_plots(cS, plot_years, Z_path, A_path)
            fprintf('   正在生成外生路径可视化图表...\n');
            T_sim = size(Z_path, 2);
            figure('Name', '外生路径可视化', 'Position', [100, 100, 1200, 500]);
            subplot(1, 2, 1);
            selected_ages = [1, 4, 8, 12, 16];
            age_labels = {};
            for i = 1:length(selected_ages)
                age_start = cS.age1_orig + (selected_ages(i) - 1) * cS.time_Step;
                age_end = age_start + cS.time_Step - 1;
                age_labels{i} = sprintf('%d-%d岁', age_start, age_end);
            end
            colors = lines(length(selected_ages));
            for i = 1:length(selected_ages)
                plot(plot_years, Z_path(selected_ages(i), :) * 100, 'Color', colors(i,:), 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 4);
                hold on;
            end
            xlabel('年份'); ylabel('人口占比 (%)'); title('人口分布路径 (Z\_path)');
            legend(age_labels, 'Location', 'best', 'FontSize', 8); grid on;
            xlim([plot_years(1), plot_years(end)]);
            subplot(1, 2, 2);
            plot(plot_years, A_path, 'b-', 'LineWidth', 2.5, 'Marker', 's', 'MarkerSize', 5, 'MarkerFaceColor', 'b');
            xlabel('年份'); ylabel('TFP水平 (相对于起始年份)'); title('全要素生产率路径 (A\_path)');
            grid on; xlim([plot_years(1), plot_years(end)]);
            avg_growth_rate = (A_path(end)/A_path(1))^(1/(T_sim*cS.time_Step)) - 1;
            text(0.05, 0.95, sprintf('平均年增长率: %.2f%%', avg_growth_rate*100), 'Units', 'normalized', 'FontSize', 10, 'BackgroundColor', 'white');
            sgtitle(sprintf('外生路径概览 (%d-%d年, %d个时期)', cS.start_year, cS.end_year, T_sim), 'FontSize', 14);
            fprintf('✅ 外生路径可视化完成。\n');
        end



        function pol_age_expanded = expand_policy_slice(pol_age_slice, nkpps)
            pol_age_expanded = pol_age_slice;
            fields = fieldnames(pol_age_expanded);
            for i = 1:length(fields)
                field_name = fields{i};
                pol_age_expanded.(field_name) = repmat(pol_age_slice.(field_name), [1, nkpps, 1]);
            end
        end

        function create_A_dynamics_plots(cS, A_path_annual, annual_years_vec)
            % =========================================================================
            % == 函数: create_exogenous_dynamics_plots
            % == 版本: [v_viz_growth_rate - 可视化增长率最终版]
            % ==
            % == 核心功能:
            % ==   1. 接收所有关键外生路径（人口、TFP、出生率、人口增长率）。
            % ==   2. 使用双Y轴绘制功能强大的“人口动态”和“技术动态”仪表盘。
            % ==   3. 智能处理年度和模型期数据，确保曲线对齐。
            % =========================================================================

            fprintf('   正在生成外生动态可视化图表 (v_viz_growth_rate)...\n');

            % --- 创建图形窗口 ---
            figure('Name', '外生路径核心动态可视化', 'Position', [100, 100, 1400, 600]);


            % --- 图 2: 技术动态 (TFP水平 vs 增长率) ---

            if ~isempty(A_path_annual)
                % 计算TFP年增长率 (年度数据)
                tfp_growth_rate_annual = (A_path_annual(2:end) ./ A_path_annual(1:end-1)) - 1;

                % 左轴: 绘制TFP水平 (归一化，年度数据)
                A_path_annual_normalized = A_path_annual / A_path_annual(1);
                yyaxis left
                plot(annual_years_vec, A_path_annual_normalized, 'k-s', 'LineWidth', 2.5, 'MarkerSize', 4, 'MarkerFaceColor', 'k');
                ylabel('TFP 水平 (相对于起始年份)', 'Color', 'k');
                ax_tech = gca;
                ax_tech.YAxis(1).Color = 'k';

                % 右轴: 绘制TFP年增长率 (年度数据)
                yyaxis right
                plot(annual_years_vec(2:end), tfp_growth_rate_annual * 100, 'g--', 'LineWidth', 2);
                hold on;
                yline(cS.g_A_ss * 100, 'g:', 'LineWidth', 1.5);
                hold off;
                ylabel('TFP 年增长率 (%)', 'Color', 'g');
                ax_tech.YAxis(2).Color = 'g';

                legend('TFP水平', 'TFP年增长率', sprintf('稳态增长率 (%.2f%%)', cS.g_A_ss*100), 'Location', 'best');
            else
                text(0.5, 0.5, 'TFP路径无法生成 (T_{sim} < 2)', 'HorizontalAlignment', 'center');
            end

            % 图形格式
            title('技术动态：TFP水平与增长率');
            xlabel('年份');
            grid on;
            xlim([cS.start_year, cS.end_year]);

            sgtitle(sprintf('外生路径概览 (%d-%d年)', cS.start_year, cS.end_year), 'FontSize', 16, 'FontWeight', 'bold');
            fprintf('✅ 外生动态可视化完成。\n');
        end

        function smooth_path = smooth_transition(start_val, end_val, num_steps, curvature)
            % =========================================================================
            % == 辅助函数: smooth_transition (v2.0 - 曲率可控版)
            % == 目的: 使用幂函数生成一个从start_val到end_val的、曲率可控的过渡路径。
            % ==
            % == 输入:
            % ==   start_val: 过渡路径的起始值
            % ==   end_val:   过渡路径的终点值
            % ==   num_steps: 过渡期的长度 (年份数)
            % ==   curvature: 曲率参数 'k'
            % ==              k > 1: 下凸式收敛 (二阶导>0, 变化率先慢后快)
            % ==              k = 1: 线性收敛 (等同于linspace)
            % ==              0 < k < 1: 上凸式收敛 (二阶导<0, 变化率先快后慢)
            % =========================================================================
            if nargin < 4
                curvature = 2; % 默认使用下凸式收敛
            end
            if num_steps <= 0
                smooth_path = [];
                return;
            end

            t = linspace(0, 1, num_steps);
            weights = t .^ curvature;
            smooth_path = (1 - weights) * start_val + weights * end_val;
        end


        function [idx_lower, idx_upper, w_upper] = find_grid_and_weights(value, gridV)
            % 功能: 在一个一维网格(gridV)中，找到给定值(value)所在的区间
            %       及其用于线性插值的权重。
            % 返回: 下界索引, 上界索引, 上界权重(w_upper)
            if value <= gridV(1)
                idx_lower = 1;
                idx_upper = 1;
                w_upper = 0.0;
                return;
            end
            if value >= gridV(end)
                idx_lower = length(gridV);
                idx_upper = length(gridV);
                w_upper = 0.0; % 或1.0，取决于如何处理边界外的值，这里设为0
                return;
            end
            idx_upper = find(gridV >= value, 1);
            idx_lower = idx_upper - 1;

            if isempty(idx_lower) || idx_lower == 0
                idx_lower = 1;
                idx_upper = 1;
                w_upper = 0;
                return;
            end

            grid_spacing = gridV(idx_upper) - gridV(idx_lower);
            if grid_spacing > 1e-9
                w_upper = (value - gridV(idx_lower)) / grid_spacing;
            else
                w_upper = 0.0;
            end
        end




        % =========================================================================
        % ==                     比较函数 (v2.0 - 格式修正版)
        % =========================================================================
        function compare_steady_states_from_files(file_initial, file_final)
            % =========================================================================
            % == 函数: compare_steady_states_from_files
            % == 版本: [v2.0 - 新报告格式兼容版]
            % ==
            % == 目的: 解析 v7.6 版本的报告，特别是区分会计投资和理论投资。
            % =========================================================================

            % --- 1. 定义要提取的变量 (关键词已更新) ---
            vars_to_extract = {
                % 宏观总量
                '国内生产总值 (Y)', 'Y', 'num';
                '总私人资本存量 (K_p)', 'K_p', 'num';
                '公共资本存量 (K_g)', 'K_g', 'num';
                '有效劳动需求 (L)', 'L', 'num';
                '私人消费 (C)', 'C', 'num';
                '私人总投资 (I_p, 会计值)', 'I_p_acct', 'num'; % [修改] 明确指定会计值
                'BGP理论投资需求', 'I_p_bgp', 'num'; % [新增] 提取BGP值
                '公共总投资 (I_g)', 'I_g', 'num';
                '政府消费 (G_c)', 'G_c', 'num';
                % 价格
                '真实利率 (r, 模型期)', 'r_annual', 'pct';
                '真实工资率 (w, 单位有效劳动)', 'w', 'num';
                % 核心比率
                '私人资本/产出比 (K_p/Y)', 'Kp_Y_ratio', 'num';
                '私人消费 (C)', 'C_Y_ratio', 'pct';
                '私人总投资 (I_p, 会计值)', 'Ip_Y_ratio', 'pct'; % [修改]
                % 检验
                '资源约束误差 (Y - 总支出)', 'RC_Error_pct', 'pct';
                '投资缺口', 'Invest_Gap_pct', 'pct';
                };

            % --- 2. 解析两个文件 ---
            fprintf('\n\n--- 正在启动稳态比较分析 ---\n');
            fprintf('   解析初始稳态报告: %s\n', file_initial);
            ss0_data = utils.parse_report_file(file_initial, vars_to_extract);

            fprintf('   解析终期稳态报告: %s\n', file_final);
            ssF_data = utils.parse_report_file(file_final, vars_to_extract);

            % --- 3. 生成并显示比较表格 ---
            fprintf('\n\n=================================================================================================\n');
            fprintf('###                             稳态结果横向比较 (Initial vs. Final)                          ###\n');
            fprintf('=================================================================================================\n');
            fprintf('%-30s | %20s | %20s | %15s\n', '变量', '初期稳态 (ss0)', '终期稳态 (ssF)', '变化 (%)');
            fprintf('%s\n', repmat('-', 1, 95));

            table_rows = {
                '--- 宏观总量 ---', '', '', '';
                '国内生产总值 (Y)', 'Y', 'num', '%.4f';
                '私人资本存量 (K_p)', 'K_p', 'num', '%.4f';
                '公共资本存量 (K_g)', 'K_g', 'num', '%.4f';
                '有效劳动需求 (L)', 'L', 'num', '%.4f';
                '私人消费 (C)', 'C', 'num', '%.4f';
                '私人投资 (会计值 I_p)', 'I_p_acct', 'num', '%.4f'; % [修改]
                '--- 价格 ---', '', '', '';
                '年化真实利率 (r)', 'r_annual', 'pct', '%.4f %%';
                '真实工资率 (w)', 'w', 'num', '%.4f';
                '--- 核心比率 ---', '', '', '';
                '私人资本/产出比 (K_p/Y)', 'Kp_Y_ratio', 'num', '%.2f';
                '消费/产出比 (C/Y)', 'C_Y_ratio', 'pct', '%.2f %%';
                '投资/产出比 (I_p/Y)', 'Ip_Y_ratio', 'pct', '%.2f %%'; % [修改]
                '--- 稳态检验 (占Y的百分比) ---', '', '', '';
                '投资缺口 (I_acct - I_bgp)', 'Invest_Gap_pct', 'pct', '%.4f %%';
                '资源约束误差', 'RC_Error_pct', 'pct', '%.4f %%';
                };

            for i = 1:size(table_rows, 1)
                label = table_rows{i, 1};
                var_name = table_rows{i, 2};

                if isempty(var_name), fprintf('%s\n', label); continue; end

                val0 = ss0_data.(var_name);
                valF = ssF_data.(var_name);
                format_spec = table_rows{i, 4};

                if val0 ~= 0 && ~isnan(val0) && ~isnan(valF)
                    if strcmp(table_rows{i,3}, 'pct_ratio') || (strcmp(table_rows{i,3}, 'pct') && ~contains(label, '利率'))
                        change_val = valF - val0;
                        change_str = sprintf('%15.2f pp', change_val); % pp = percentage points
                    else
                        change_pct = ((valF / val0) - 1) * 100;
                        change_str = sprintf('%15.2f %%', change_pct);
                    end
                else
                    change_str = sprintf('%15s', 'N/A');
                end

                if contains(format_spec, '%%')
                    fprintf(['%-30s | ' format_spec ' | ' format_spec ' | %s\n'], label, val0, valF, change_str);
                else
                    fprintf(['%-30s | ' format_spec ' | ' format_spec ' | %s\n'], label, val0, valF, change_str);
                end
            end
            fprintf('%s\n', repmat('-', 1, 95));
        end

        function data_struct = parse_report_file(filename, vars_def)
            % [核心修正] 更新正则表达式以适应新的、更明确的关键词
            data_struct = struct();
            try
                fid = fopen(filename, 'r', 'n', 'UTF-8');
                if fid == -1, error('无法打开文件: %s', filename); end
                content = fread(fid, '*char')';
                fclose(fid);
            catch ME
                error('读取文件 %s 时出错: %s', filename, ME.message);
            end

            for i = 1:size(vars_def, 1)
                keyword = vars_def{i, 1};
                var_name = vars_def{i, 2};
                data_type = vars_def{i, 3};

                % 使用更灵活的正则表达式来处理括号和附加词
                % 它会匹配 "关键词" 后面跟着可选的 "(...)"，然后是 "..." 和 ":"
                escaped_keyword = regexptranslate('escape', keyword);

                if strcmp(data_type, 'num')
                    pat = [escaped_keyword, '.*?:\s*([-\d\.]+)'];
                elseif strcmp(data_type, 'pct')
                    pat = [escaped_keyword, '.*?([-\d\.]+[eE]*[-\d+]*)\s*\%'];
                end

                tokens = regexp(content, pat, 'tokens');

                if ~isempty(tokens)
                    data_struct.(var_name) = str2double(tokens{1}{1});
                else
                    data_struct.(var_name) = NaN; % 如果找不到则标记为NaN
                end
            end
        end

        function report = check_k_grid_limit(Dist, kGridV, ss_or_t, grid_name)
    % = 'kpps'
    % =========================================================================
    % == 函数: check_k_grid_limit
    % == 目的: 检查分布在资产网格最后一个点的质量集中情况。
    % == 输入:
    % ==   Dist: 人口分布 (nk x nkpps x nw x na [x T])
    % ==   kGridV: 与第一个维度对应的资产网格向量
    % ==   ss_or_t: 当前是稳态('ss')还是特定时期t (整数)
    % ==   grid_name: 'k_private' 或 'kpps' 用于报告
    % =========================================================================
    
    total_mass = sum(Dist, 'all');
    if total_mass < 1e-9
        report = sprintf('   - [%s Grid Check t=%s]: Total mass is zero, skipping.\n', grid_name, num2str(ss_or_t));
        fprintf(report);
        return;
    end
    
    % 获取第一个维度的最后一个切片
    mass_on_last_point = sum(Dist(end, :, :, :), 'all');
    
    percentage_on_last_point = (mass_on_last_point / total_mass) * 100;
    
    k_max = kGridV(end);
    
    fprintf('   --- [诊断报告] %s 网格上限检查 (t=%s) ---\n', grid_name, num2str(ss_or_t));
    fprintf('      网格上限 (k_max): %.2f\n', k_max);
    fprintf('      最后一个网格点上的总人口质量: %.4e\n', mass_on_last_point);
    fprintf('      占总人口的百分比: %.9f %%\n', percentage_on_last_point);
    
    if percentage_on_last_point > 0.5
        fprintf('      >>> [警告] 超过 0.5%% 的人口挤在网格上限！\n');
        fprintf('          强烈建议将 k_max 提高到至少 %.2f 或更高。\n', k_max * 1.5);
    elseif percentage_on_last_point > 0.1
        fprintf('      >>> [注意] 有少量人口 (%.4f%%) 位于网格上限。\n', percentage_on_last_point);
        fprintf('          当前的 k_max 可能足够，但可以考虑适当提高以增加精度。\n');
    else
        fprintf('      >>> [✅ 通过] 网格上限设置合理。只有极少的人口质量在最后一个点。\n');
    end
    
    report = sprintf('      Percentage on last point for %s at t=%s: %.4f%%\n', grid_name, num2str(ss_or_t), percentage_on_last_point);
end

    end
end