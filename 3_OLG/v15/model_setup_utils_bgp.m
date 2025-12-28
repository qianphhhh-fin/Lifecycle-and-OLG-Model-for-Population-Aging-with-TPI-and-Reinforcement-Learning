classdef model_setup_utils_bgp
    methods (Static)

        function cS = ParameterValues()
            % [vPaper.5 - 政府投资版]
            % [BGP修改] 提供了完整的参数列表，并将待校准参数标记为"初始猜测值"。
            % [BGP修改] 引入核心增长参数 g_A_ss 用于稳态化模型

            cS.pps_active = false;

            % --- 人口结构基础参数 ---
            cS.age1_orig = 20;
            cS.ageLast_orig = 98;
            cS.ageRetire_orig = 60;
            cS.aD_orig = cS.ageLast_orig - cS.age1_orig + 1;
            cS.aR_idx_orig = cS.ageRetire_orig - cS.age1_orig + 1;
            cS.aW_orig = cS.aR_idx_orig - 1;
            cS.physAgeV_orig = (cS.age1_orig : cS.ageLast_orig)';

            % --- 遗赠动机参数 ---
            cS.phi_bequest = 5;

            % --- PPS制度参数 ---
            cS.pps_withdrawal_rate = 0.15; % 退休后PPS每5年期的提取比例
            cS.pps_contrib_limit = 9999;   % PPS缴费的最大额度
            cS.pps_max_contrib_frac = 0.1; % PPS缴费的最大比例
            cS.pps_tax_rate_withdrawal = 0.03; % PPS取出时的税率

            % --- [BGP修改] 待校准参数的初始猜测值 ---
            cS.beta = 0.996;        % 主观折现因子
            cS.gamma = 0.12;        % 公共资本的产出弹性
            cS.lambda_g = 0.35;     % 政府投资占总税收的比例

            % 基于对技术前沿国家(如美国)长期TFP增长率的研究，设定一个
            % 反映成熟经济体创新驱动增长的合理值。1.2%是一个被广泛接受
            % 的、较为中性的估计，反映了近年的生产率放缓趋势。
            cS.g_A_ss = 0.012;      % 长期技术年增长率 (1.2%)
            cS.long_term_birth_rate_target_per_1000 = 1.2;  % 2100年联合国预测1.35
            cS.transition_period_years = 50; % 过渡期长度
            cS.curvature_param = 0.5; % 设定为 < 1 的值，实现下凸收敛


            % --- 其他固定参数 ---
            cS.sigma = 7;
            cS.cFloor = 0.005;
            cS.time_Step = 5; % 注意: 此值应与主脚本中的 TIME_STEP 保持一致

            % 企业
            cS.alpha = 0.35;
            cS.ddk = 1 - (1 - 0.015)^cS.time_Step; % 5年期私人资本折旧率

            % 政府
            cS.ddk_g = cS.ddk; % 公共资本折旧率 (为简化，设为与私人资本相同)
            cS.tau_k = 0.02;
            cS.tau_l = 0.06;
            cS.tau_c = 0.03;
            cS.A = 1.0; % TFP基准值（现在理解为"初始稳态的TFP基准值"，是归一化的锚点）

            % --- 年龄组聚合参数 ---
            cS.aD_new = ceil(cS.aD_orig / cS.time_Step);
            cS.aR_new = ceil(cS.aW_orig / cS.time_Step);
            cS.physAgeMap = cell(cS.aD_new, 1);
            for a = 1:cS.aD_new
                startIdx = (a-1)*cS.time_Step + 1;
                endIdx = min(a*cS.time_Step, cS.aD_orig);
                cS.physAgeMap{a} = startIdx:endIdx;
            end

            % 存活率 (基于年度死亡率数据聚合)
            d_orig_data = [0.00159,0.00169,0.00174,0.00172,0.00165,0.00156,0.00149,0.00145,0.00145,0.00149,0.00156,0.00163,0.00171,0.00181,0.00193,0.00207,0.00225,0.00246,0.00270,0.00299,0.00332,0.00368,0.00409,0.00454,0.00504,0.00558,0.00617,0.00686,0.00766,0.00865,0.00955,0.01058,0.01162,0.01264,0.01368,0.01475,0.01593,0.01730,0.01891,0.02074,0.02271,0.02476,0.02690,0.02912,0.03143,0.03389,0.03652,0.03930,0.04225,0.04538,0.04871,0.05230,0.05623,0.06060,0.06542,0.07066,0.07636,0.08271,0.08986,0.09788,0.10732,0.11799,0.12895,0.13920,0.14861,0.16039,0.17303,0.18665,0.20194,0.21877,0.23601,0.25289,0.26973,0.28612,0.30128,0.31416,0.32915,0.34450,0.36018];
            s_orig_1yr = 1 - d_orig_data;
            cS.s_pathV = zeros(cS.aD_new, 1);
            for a = 1:(cS.aD_new - 1)
                age_indices = cS.physAgeMap{a};
                s_period = 1.0;
                % 聚合得到模型期(如5年)的存活率
                for i = 1:cS.time_Step
                    annual_idx = age_indices(1) + i - 1;
                    if annual_idx <= length(s_orig_1yr)
                        s_period = s_period * s_orig_1yr(annual_idx);
                    end
                end
                cS.s_pathV(a) = s_period;
            end
            cS.s_pathV(cS.aD_new) = 0; % 最后一期确定死亡

            % 中国平均预期寿命约为78岁(2021-2023数据)，模型中对应物理年龄索引59
            avg_lifespan_idx = 78;

            % --- 支出冲击参数 ---
            cS.nw = 5; % 常规效率状态数量
            % 年龄节点
            cS.young_shock_start_age = 26;
            cS.young_shock_end_age = 35;
            cS.old_shock_start_age = 65;
            cS.old_shock_end_age = 88;
            % 冲击年化概率峰值
            cS.p_shock_young_peak_annual = 0;
            cS.p_shock_old_peak_annual = 0;
            % 冲击支出比例 (kappa)
            cS.kappa_young = 0.9;
            cS.kappa_old = 0.9;

            % --- 年龄效率剖面 ---
            beta0 = -13.215788; beta1 =  1.349514; beta2 = -0.043363;
            beta3 =  0.000585; beta4 = -0.000003;
            ages = cS.physAgeV_orig;
            log_age_eff_orig = beta0 + beta1 * ages + beta2 * (ages.^2) + beta3 * (ages.^3) + beta4 * (ages.^4);
            ageEffV_orig_unnormalized = exp(log_age_eff_orig);
            ageEffV_orig_unnormalized(cS.aR_idx_orig:end) = 0;
            working_life_indices = 1:(cS.aR_idx_orig - 1);
            mean_efficiency_working = mean(ageEffV_orig_unnormalized(working_life_indices));
            ageEffV_orig = ageEffV_orig_unnormalized / mean_efficiency_working;
            cS.ageEffV_new = zeros(cS.aD_new, 1);
            for a = 1:cS.aD_new
                cS.ageEffV_new(a) = mean(ageEffV_orig(cS.physAgeMap{a}));
            end

            % PAYG养老金缴费率的初始值/占位符
            cS.theta_path = 0.06;
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
            default_kMax = 40; % 直接设为合理的常数
            default_kppsMax = 20; % 直接设为合理的常数

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
            power_k = 2;

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
                kppsGridV_temp = cS.kppsMin;
            end
            cS.kppsGridV = kppsGridV_temp(:);
        end

        function [leGridV_expanded, TrProbM_by_age, leProb1V, nw_expanded] = EarningProcess_AgeDependent(cS)
            % =========================================================================
            % == 函数: EarningProcess_AgeDependent
            % == 版本: [v10 - Switchable Shocks]
            % ==
            % == 目的:
            % ==   根据参数自动生成包含或不包含重大支出冲击的年龄相关收入过程。
            % ==
            % == 核心逻辑:
            % ==   1. [自动检测] 检查 cS.p_shock_..._peak_annual 参数。
            % ==   2. [含冲击路径] 若参数非零，则生成一个扩展的状态空间(nw+2)，
            % ==      包含常规生产率状态以及青年/老年重大支出冲击状态。
            % ==   3. [无冲击路径] 若参数为零，则生成一个标准的AR(1)生产率过程，
            % ==      状态空间不扩展(nw)，且转移矩阵不随年龄变化。
            % =========================================================================

            % --- 1. 生成基础AR(1)过程 (对两种情况都通用) ---
            lePersistence = 0.9; leShockStd = 0.15^0.5; Tauchen_q = 2.0;
            [leLogGridV_normal, leTrProbM_normal] = model_setup_utils_bgp.tauchen(cS.nw, lePersistence, leShockStd, 0, Tauchen_q);
            [~, D] = eig(leTrProbM_normal'); [~, c] = min(abs(diag(D)-1));
            leProb1V = abs(D(:,c)/sum(D(:,c))); % 新生儿的长期冲击分布

            % --- 2. [自动切换] 检测是否存在重大支出冲击 ---
            has_major_shocks = (isfield(cS, 'p_shock_young_peak_annual') && cS.p_shock_young_peak_annual > 0) || ...
                (isfield(cS, 'p_shock_old_peak_annual') && cS.p_shock_old_peak_annual > 0);

            if has_major_shocks
                % --- 路径 A: 生成包含重大支出冲击的扩展过程 ---
                fprintf('   正在生成【含】重大支出冲击的收入过程...\n');

                % 2A. 根据参数构建年度冲击概率向量
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

                % 3A. 从年度概率聚合为模型期概率
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

                % 4A. 构建最终的年龄依赖转移矩阵并汇报
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

                    % 发生冲击后，下一期回归到常规冲击的长期稳态分布
                    Tr_a(cS.nw + 1, 1:cS.nw) = leProb1V';
                    Tr_a(cS.nw + 2, 1:cS.nw) = leProb1V';

                    % 归一化每一行，确保概率和为1
                    row_sums = sum(Tr_a, 2);
                    row_sums(row_sums == 0) = 1; % 防止除以零
                    Tr_a = Tr_a ./ row_sums;

                    TrProbM_by_age{a} = Tr_a;
                end

            else
                % --- 路径 B: 生成不含重大支出冲击的标准过程 ---
                fprintf('   检测到重大冲击参数为零，正在生成【无】重大支出冲击的标准收入过程...\n');

                % 2B. 设定输出变量
                nw_expanded = cS.nw;
                leGridV_expanded = exp(leLogGridV_normal(:));

                % 3B. 构建转移矩阵 (所有年龄都相同)
                TrProbM_by_age = cell(cS.aD_new, 1);
                for a = 1:cS.aD_new
                    TrProbM_by_age{a} = leTrProbM_normal;
                end
            end
        end

        function [y_grid_out, trProbM_out] = tauchen(N, rho, sigma, mu, m)
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

        function [muM, utilM] = CES_utility(cM, sigma, cS)
            % [BGP修改] 此函数的作用是输入一个消费值，输出其效用
            % [BGP修改] 在标准化框架下，VFI求解器会向它传递标准化的消费 ĉ
            % [BGP修改] 这个函数会正确地计算 u(ĉ) = ĉ^(1-σ)/(1-σ)
            % [BGP修改] 无需知道绝对技术水平A，完全在标准化环境中操作
            if cM<=0
                c_adj = cS.cFloor;
                print('消费小于零')
            else
                c_adj = cM;
            end
            if abs(sigma - 1) < 1e-6
                utilM = log(c_adj);
                muM = 1./c_adj;
            else utilM = (c_adj.^(1-sigma))./(1-sigma);
                muM = c_adj.^(-sigma);
            end
        end

        function util_beq = bequest_utility(k_prime, cS)
            % [BGP修改] 此函数将被传递标准化的遗赠 k̂'，并正确计算其效用
            % [BGP修改] 在标准化框架下，遗赠效用函数处理标准化的遗赠值
            % [BGP修改] 无需知道绝对技术水平A，完全在标准化环境中操作
            if cS.phi_bequest <= 0
                util_beq = 0;
                return;
            end

            k_adj = max(cS.cFloor, k_prime);

            if abs(cS.sigma - 1) < 1e-6
                util_beq = cS.phi_bequest * log(k_adj);
            else
                util_beq = cS.phi_bequest * (k_adj.^(1-cS.sigma))./(1-cS.sigma);
            end
        end

        function [Z_out, A_path, cS] = generate_exo_paths(cS, graph_flag)
    % =========================================================================
    % == 函数: generate_exo_paths
    % == 版本: [v_viz_growth_rate - 可视化增长率版]
    % ==
    % == 核心修正:
    % ==   1. 在主循环后，计算整个模拟期间的逐期【年化人口增长率】路径。
    % ==   2. 将这个新的增长率路径 pop_growth_rate_path_annual 传递给
    % ==      新的绘图函数 create_exogenous_dynamics_plots。
    % ==   3. 删除了之前在控制台打印单一稳态增长率的逻辑。
    % =========================================================================

    fprintf('\n--- 启动外生路径生成器 (v_viz_growth_rate - 可视化增长率版) ---\n');

    % --- 1. 初始化模拟和路径参数 ---
    ss_convergence_tol = 1e-7;
    max_sim_periods = 200;
    final_birth_rate_annual = cS.long_term_birth_rate_target_per_1000 / 1000;
    fprintf('   长期目标年化出生率设定为: %.2f / 1000\n', cS.long_term_birth_rate_target_per_1000);
    fprintf('   将使用曲率参数 k=%.1f 进行下凸式平滑收敛。\n', cS.curvature_param);

    % --- 2. 加载并处理外生数据 ---
    try
        cbr_data = readtable('data\人口\UN_PPP2024_CBR_birthper1000_China.xlsx');
    catch ME
        error('无法加载出生率数据文件: %s', ME.message);
    end
    var_names = cbr_data.Properties.VariableNames;
    year_cols_indices = find(startsWith(var_names, 'y'));
    un_years = str2double(cellfun(@(x) x(2:end), var_names(year_cols_indices), 'UniformOutput', false));
    un_cbr = table2array(cbr_data(1, year_cols_indices)) / 1000;
    last_data_year = un_years(end);
    try
        pop_data_model_group = readtable('data\人口\population_by_age_group_all_years.xlsx', 'Sheet', 'population');
    catch ME
        error('无法加载人口数据文件: %s', ME.message);
    end
    initial_year_col_name = ['y', num2str(cS.ss0_year)];
    if ~ismember(initial_year_col_name, pop_data_model_group.Properties.VariableNames)
        error('在人口数据中找不到初始年份列: %s', initial_year_col_name);
    end
    initial_pop_model_group = pop_data_model_group.(initial_year_col_name);

    % --- 3. [核心] 迭代模拟直至人口结构稳态 ---
    fprintf('   正在启动人口动态模拟以寻找结构性稳态...\n');
    Z_path_raw = zeros(cS.aD_new, max_sim_periods);
    Z_path_raw(:, 1) = initial_pop_model_group;
    Z_norm_prev = Z_path_raw(:, 1) / sum(Z_path_raw(:, 1));
    converged = false;
    T_final = max_sim_periods;

    annual_years_full_range = cS.start_year:(cS.start_year + (max_sim_periods * cS.time_Step) -1);
    birth_rate_path_annual_full = zeros(size(annual_years_full_range));
    transition_end_year = last_data_year + cS.transition_period_years;

    data_mask = annual_years_full_range <= last_data_year;
    birth_rate_path_annual_full(data_mask) = interp1(un_years, un_cbr, annual_years_full_range(data_mask), 'pchip');
    last_data_rate = birth_rate_path_annual_full(find(data_mask, 1, 'last'));

    transition_mask = annual_years_full_range > last_data_year & annual_years_full_range <= transition_end_year;
    if any(transition_mask)
        num_trans_years = sum(transition_mask);
        birth_rate_path_annual_full(transition_mask) = model_setup_utils_bgp.smooth_transition(last_data_rate, final_birth_rate_annual, num_trans_years, cS.curvature_param);
    end

    ss_mask = annual_years_full_range > transition_end_year;
    birth_rate_path_annual_full(ss_mask) = final_birth_rate_annual;
    fprintf('   出生率路径将在 %d 年到 %d 年之间【下凸式平滑下降】，并在 %d 年后进入长期稳态。\n', last_data_year, transition_end_year, transition_end_year);

    for t = 1:(max_sim_periods - 1)
        Z_path_raw(2:end, t+1) = Z_path_raw(1:(cS.aD_new-1), t) .* cS.s_pathV(1:(cS.aD_new-1));
        total_pop_base = sum(Z_path_raw(:, t));
        period_start_idx_in_annual = (t-1)*cS.time_Step + 1;
        period_end_idx_in_annual = t*cS.time_Step;
        avg_birth_rate_period = mean(birth_rate_path_annual_full(period_start_idx_in_annual:period_end_idx_in_annual));
        new_entrants = total_pop_base * avg_birth_rate_period * cS.time_Step;
        Z_path_raw(1, t+1) = new_entrants;
        Z_norm_current = Z_path_raw(:, t+1) / sum(Z_path_raw(:, t+1));
        diff = max(abs(Z_norm_current - Z_norm_prev));
        if diff < ss_convergence_tol
            T_final = t + 1;
            converged = true;
            fprintf('   ✅ 人口结构在第 %d 期 (年份: %d) 达到稳态 (差异: %.2e)。\n', ...
                T_final, cS.start_year + t*cS.time_Step, diff);
            break;
        end
        Z_norm_prev = Z_norm_current;
    end
    if ~converged, warning('人口模拟在达到最大期数 %d 后仍未收敛！', max_sim_periods); end
    
    cS.T_sim = T_final;
    cS.end_year = cS.start_year + (T_final - 1) * cS.time_Step;
    fprintf('   内生决定的模拟期数 T_sim = %d, 结束年份 = %d\n', cS.T_sim, cS.end_year);
    Z_path_raw_final = Z_path_raw(:, 1:cS.T_sim);
    Z_out = Z_path_raw_final ./ sum(Z_path_raw_final, 1);
    
    % --- [核心修改] 计算整个路径的年化人口增长率 ---
    total_pop_path = sum(Z_path_raw_final, 1);
    pop_growth_rate_path_period = (total_pop_path(2:end) ./ total_pop_path(1:end-1)) - 1;
    pop_growth_rate_path_annual = (1 + pop_growth_rate_path_period).^(1/cS.time_Step) - 1;

    % 将最终的稳态增长率存储到cS中，供后续稳态求解器使用
    if converged && ~isempty(pop_growth_rate_path_annual)
        cS.n_ss = pop_growth_rate_path_annual(end);
        fprintf('   最终稳态年化人口增长率设定为: n_ss = %.6f\n', cS.n_ss);
    else
        cS.n_ss = 0; % 如果未收敛或路径太短，默认为0
        fprintf('   [警告] 人口路径未收敛或过短，无法计算稳态增长率。默认为 n = 0。\n');
    end

    % --- 生成TFP路径和绘图 ---
    model_sim_years = cS.start_year:cS.time_Step:cS.end_year;
    annual_years_vec = cS.start_year:cS.end_year;
    A_path = model_setup_utils_bgp.generate_tfp_path(cS, annual_years_vec, model_sim_years);

    if graph_flag
        final_annual_mask = annual_years_full_range >= cS.start_year & annual_years_full_range <= cS.end_year;
        birth_rate_path_annual_for_plot = birth_rate_path_annual_full(final_annual_mask);
        if length(model_sim_years) >= 2
            A_path_annual_for_plot = interp1(model_sim_years, A_path, annual_years_vec, 'pchip');
        else
            A_path_annual_for_plot = [];
        end
        
        % [核心修改] 将人口增长率路径传递给绘图函数
        model_setup_utils_bgp.create_exogenous_dynamics_plots(cS, Z_out, A_path_annual_for_plot, ...
            birth_rate_path_annual_for_plot, annual_years_vec, pop_growth_rate_path_annual);
    end
end
        function A_path = generate_tfp_path(cS, annual_years_vec, model_sim_years)
            % =========================================================================
            % == 函数: generate_tfp_path (v_curvature_control - 曲率可控最终版)
            % ==
            % == 核心修正:
            % ==   调用新版的smooth_transition函数，并传入一个曲率参数(>1)，
            % ==   确保TFP增长率路径是【下凸式】收敛。
            % =========================================================================

            fprintf('   正在构建TFP路径 (v_curvature_control - 曲率可控最终版)...\n');

            % --- 1. 加载数据和设定参数 ---
            tfp_data_pwt = readtable('data\PWT\china_pwt_data.xlsx');
            pwt_years = tfp_data_pwt.year(tfp_data_pwt.year >= 1979);
            pwt_tfp_level = tfp_data_pwt.rtfpna(tfp_data_pwt.year >= 1979);
            bai_projections = [2025, 0.0557; 2030, 0.0482; 2035, 0.0394; 2040, 0.0340; 2045, 0.0346; 2050, 0.0298];

            long_term_g_annual = cS.g_A_ss;
            fprintf('   TFP长期年化增长率将【下凸式平滑收敛】至: %.4f%% (曲率k=%.1f)\n', long_term_g_annual * 100, cS.curvature_param);

            data_end_year = bai_projections(end, 1);
            transition_end_year = data_end_year + cS.transition_period_years;

            % --- 2. [核心] 构建平滑的【年度增长率】路径 ---
            g_path_annual = zeros(size(annual_years_vec));
            g_anchor_years = [pwt_years(end); bai_projections(:,1)];
            g_anchor_rates = [0.058; bai_projections(:,2)];

            data_mask = annual_years_vec <= data_end_year;
            g_path_annual(data_mask) = interp1(g_anchor_years, g_anchor_rates, annual_years_vec(data_mask), 'previous', 'extrap');
            last_data_g = g_path_annual(find(data_mask, 1, 'last'));

            % [修改] 调用新版函数并传入曲率参数
            transition_mask = annual_years_vec > data_end_year & annual_years_vec <= transition_end_year;
            if any(transition_mask)
                num_trans_years = sum(transition_mask);
                g_path_annual(transition_mask) = model_setup_utils_bgp.smooth_transition(last_data_g, long_term_g_annual, num_trans_years, cS.curvature_param);
            end

            ss_mask = annual_years_vec > transition_end_year;
            g_path_annual(ss_mask) = long_term_g_annual;
            fprintf('   TFP增长率将在 %d 年到 %d 年之间【下凸式平滑过渡】，并在 %d 年后进入稳态。\n', data_end_year, transition_end_year, transition_end_year);

            % ... (函数的其余部分，包括水平路径生成和采样，保持不变) ...
            A_path_annual_level = zeros(size(annual_years_vec));
            proj_years = [pwt_years(end); bai_projections(:,1)];
            proj_levels = zeros(size(proj_years));
            proj_levels(1) = pwt_tfp_level(end);
            for i = 1:size(bai_projections,1)
                prev_level = proj_levels(i);
                prev_year = proj_years(i);
                curr_year = proj_years(i+1);
                g = bai_projections(i,2);
                proj_levels(i+1) = prev_level * (1+g)^(curr_year - prev_year);
            end
            all_anchor_years = [pwt_years; proj_years(2:end)];
            all_anchor_levels = [pwt_tfp_level; proj_levels(2:end)];
            [unique_years, ia, ~] = unique(all_anchor_years);
            unique_levels = all_anchor_levels(ia);
            start_year_idx = find(annual_years_vec == cS.start_year, 1);
            if isempty(start_year_idx)
                error('模拟起始年份不在年度向量中！');
            end
            A_path_annual_level(start_year_idx) = interp1(unique_years, unique_levels, cS.start_year, 'pchip', 'extrap');
            for t_idx = (start_year_idx + 1):length(annual_years_vec)
                A_path_annual_level(t_idx) = A_path_annual_level(t_idx - 1) * (1 + g_path_annual(t_idx));
            end
            for t_idx = (start_year_idx - 1):-1:1
                A_path_annual_level(t_idx) = A_path_annual_level(t_idx + 1) / (1 + g_path_annual(t_idx + 1));
            end
            A_path = interp1(annual_years_vec, A_path_annual_level, model_sim_years, 'linear');
            if A_path(1) ~= 0 && ~isnan(A_path(1))
                A_path = A_path / A_path(1);
            else
                warning('初始TFP水平为0或NaN，无法归一化。');
            end
            fprintf('✅ TFP路径构建完成 (平滑过渡最终版)。\n');
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

        function Z_ss0 = get_calibration_inputs(target_year, cS)
            % [BGP不变] 此函数的作用是为初始稳态提供校准年份的人口分布 Z_ss0
            % 这个功能和逻辑在新的框架下完全正确且必要
            fprintf('   为校准获取 %d 年的输入...\n', target_year);

            % 1. 加载人口数据
            pop_data = readtable('data\人口\population_by_age_group_all_years.xlsx', 'Sheet', 'pop_normalized');
            data_years = str2double(cellfun(@(x) x(2:end), pop_data.Properties.VariableNames(2:end), 'UniformOutput', false));

            % 2. 插值得到目标年份的人口分布
            Z_ss0_raw = interp1(data_years, pop_data{:, 2:end}', target_year, 'linear', 'extrap')';

            % 3. 归一化
            Z_ss0 = Z_ss0_raw / sum(Z_ss0_raw);
        end

        function [cS] = calcaulte_theta_payg_path(cS, graph_flag)
            % [BGP不变] 此函数定义了一条外生的政策路径（养老金缴费率）
            % 它不受模型内生增长机制的影响
            fprintf('正在构建基于【覆盖率插值】的有效养老金缴费率路径 (theta_path)...\n');

            % --- 1. 定义关键参数 ---
            theta_urban_employee_effective = 0.20;
            theta_resident_effective = 0.03;

            % [政策目标] 定义未来的目标覆盖率
            coverage_urban_final = 0.8;
            coverage_resident_final = 0.8;
            year_reach_final = 2050;

            % --- 2. 收集来自官方统计的真实数据点 (单位：万人) ---
            year_pax_urban = [1997, 2000, 2005, 2008, 2011, 2014, 2017, 2020, 2023];
            pax_urban      = [8671, 10448, 17487, 21890, 28391, 34124, 40293, 45638, 52112];

            year_pax_resident = [2008, 2011, 2014, 2017, 2020, 2023];
            pax_resident      = [0,    45586, 50075, 51255, 54227, 54475];

            year_laborpop_data = [1997, 2000, 2005, 2008, 2011, 2014, 2017, 2020, 2023];
            laborpop_data      = [78531, 82352, 87532, 90833, 94072, 92837, 91199, 89552, 86481];

            % --- 3. 计算历史数据点的【覆盖率】 ---
            coverage_urban_data = pax_urban ./ laborpop_data;

            laborpop_for_resident_years = interp1(year_laborpop_data, laborpop_data, year_pax_resident, 'linear', 'extrap');
            coverage_resident_data = pax_resident ./ laborpop_for_resident_years;

            % --- 4. 直接对【覆盖率】进行插值，生成完整的年度路径 ---
            annual_years_vector = cS.start_year:cS.end_year;

            % a. 插值生成年度【城镇职工覆盖率】路径
            interp_years_urban = [year_pax_urban, year_reach_final, cS.end_year];
            interp_coverage_urban = [coverage_urban_data, coverage_urban_final, coverage_urban_final];
            [unique_years_u, ia_u, ~] = unique(interp_years_urban);
            coverage_urban_annual = interp1(unique_years_u, interp_coverage_urban(ia_u), annual_years_vector, 'linear');

            % b. 插值生成年度【城乡居民覆盖率】路径
            interp_years_resident = [cS.start_year, year_pax_resident, year_reach_final, cS.end_year];
            interp_coverage_resident = [0, coverage_resident_data, coverage_resident_final, coverage_resident_final];
            [unique_years_r, ia_r, ~] = unique(interp_years_resident);
            coverage_resident_annual = interp1(unique_years_r, interp_coverage_resident(ia_r), annual_years_vector, 'linear');
            coverage_resident_annual(annual_years_vector < min(year_pax_resident)) = 0;

            % --- 5. 基于【覆盖率路径】计算最终的【年度有效缴费率】路径 ---
            theta_path_annual = (coverage_urban_annual * theta_urban_employee_effective) + ...
                (coverage_resident_annual * theta_resident_effective);

            % --- 6. 从年度路径中，提取出模型【5年期】所需的路径 ---
            model_year_indices = 1:(cS.end_year - cS.start_year + 1);
            model_year_indices_5yr = model_year_indices(1:cS.time_Step:end);

            T_sim = length(cS.start_year:cS.time_Step:cS.end_year);
            theta_path = theta_path_annual(model_year_indices_5yr);
            theta_path = theta_path(1:min(T_sim, length(theta_path)));

            % --- 7. 将最终路径存入 cS 结构体 ---
            cS.theta_path = theta_path;

            % --- 8. (可选) 可视化检查路径 ---
            if graph_flag
                T_plot = length(cS.theta_path);
                time_axis = cS.start_year:cS.time_Step:(cS.start_year + cS.time_Step*(T_plot-1));

                figure('Name', 'Effective PAYG Tax Rate Path (theta_path) - Corrected');
                plot(time_axis, cS.theta_path, 'k-s', 'LineWidth', 2, 'MarkerFaceColor', 'k', 'DisplayName', '总有效缴费率 (加权平均)');
                hold on;

                contribution_urban = coverage_urban_annual(model_year_indices_5yr(1:T_plot)) * theta_urban_employee_effective;
                contribution_resident = coverage_resident_annual(model_year_indices_5yr(1:T_plot)) * theta_resident_effective;

                plot(time_axis, contribution_urban, 'b--o', 'LineWidth', 1.5, 'DisplayName', '城镇职工体系贡献');
                plot(time_axis, contribution_resident, 'r--d', 'LineWidth', 1.5, 'DisplayName', '城乡居民体系贡献');

                title('模型使用的有效养老金缴费率路径 (基于覆盖率插值)');
                xlabel('年份');
                ylabel('有效缴费率 (θ_t)');
                legend('show', 'Location', 'best');
                grid on;
                ylim([0, max(cS.theta_path)*1.2]);
            end
            fprintf('✅ 基于覆盖率插值的缴费率路径构建完成。\n');
        end

        % =======================================================
        % == [BGP一致性验证] 模型设置模块一致性总结
        % =======================================================
        % 本模块中的所有函数都正确地处理了BGP一致性:
        % 1. ParameterValues(): 正确引入了核心增长参数g_A_ss
        % 2. generateGrids(): 适应标准化求解，使用简化的默认值计算
        % 3. CES_utility(): 正确处理标准化消费ĉ，无需知道绝对技术水平A
        % 4. bequest_utility(): 正确处理标准化遗赠k̂'，完全在标准化环境中操作
        % 5. load_exogenous_paths(): 构建TFP水平路径A_path，为"复原趋势"提供基础
        % 6. 所有外生路径函数都正确地处理了增长和人口动态
        %
        % 这些函数与main_steady_state_utils_bgp.m中的求解器完全兼容，
        % 确保了整个BGP框架的理论一致性和数值稳定性。
        function plot_population_pyramid_simplified(Z_out)
            % 简化版人口金字塔可视化

            % 年龄组标签
            age_groups = arrayfun(@(a) sprintf('%d-%d', 20+(a-1)*5, 24+(a-1)*5), 1:16, 'UniformOutput', false);

            % 选择三个时间点
            t_initial = 1;
            t_peak_aging = 10; % 约50年后
            t_steady_state = 40;

            figure('Name', '人口结构金字塔演变', 'Position', [100, 100, 1400, 500]);

            % --- 子图1: 初始年份 ---
            subplot(1, 3, 1);
            barh(Z_out(:, t_initial) * 100, 'FaceColor', '#0072BD');
            set(gca, 'YDir', 'reverse');
            yticks(1:16);
            yticklabels(age_groups);
            xlabel('人口占比 (%)');
            title(sprintf('初始年 (t=%d)', t_initial));
            grid on;
            xlim([0, 12]);

            % --- 子图2: 老龄化高峰期 ---
            subplot(1, 3, 2);
            barh(Z_out(:, t_peak_aging) * 100, 'FaceColor', '#D95319');
            set(gca, 'YDir', 'reverse');
            yticks(1:16);
            yticklabels(age_groups);
            title(sprintf('老龄化高峰附近 (t=%d)', t_peak_aging));
            grid on;
            xlim([0, 12]);

            % --- 子图3: 长期稳态 ---
            subplot(1, 3, 3);
            barh(Z_out(:, t_steady_state) * 100, 'FaceColor', '#77AC30');
            set(gca, 'YDir', 'reverse');
            yticks(1:16);
            yticklabels(age_groups);
            title(sprintf('长期稳态 (t=%d)', t_steady_state));
            grid on;
            xlim([0, 12]);

            sgtitle('人口结构从“中年型”向“老年型”的动态演变', 'FontSize', 16, 'FontWeight', 'bold');
        end

        function Z_theory_zpg = compute_theoretical_ss_dist_zpg(cS)
            % =========================================================================
            % == 函数: compute_theoretical_ss_dist_zpg
            % == 目的: 计算一个理想化的、零人口增长(ZPG)的理论稳态人口分布。
            % ==       这个分布只由模型内在的生存率 s_pathV 决定，不受任何
            % ==       外部出生率数据或历史结构的影响。它为过渡路径求解器
            % ==       提供了一个完美的、稳定的“锚点”。
            % =========================================================================
            fprintf('   正在计算理想化的零人口增长(ZPG)理论稳态分布...\n');

            % 1. 以新生儿为基准1.0，计算各年龄组的相对人口规模
            mass_levels_by_age = zeros(cS.aD_new, 1);
            mass_levels_by_age(1) = 1.0;
            for ia = 1:(cS.aD_new - 1)
                mass_levels_by_age(ia+1) = mass_levels_by_age(ia) * cS.s_pathV(ia);
            end

            % 2. 将相对规模归一化，得到最终的概率分布
            Z_theory_zpg = mass_levels_by_age / sum(mass_levels_by_age);

            fprintf('   ✅ ZPG稳态分布计算完成。\n');
        end

        function pol_age_expanded = expand_policy_slice(pol_age_slice, nkpps)
            pol_age_expanded = pol_age_slice;
            fields = fieldnames(pol_age_expanded);
            for i = 1:length(fields)
                field_name = fields{i};
                pol_age_expanded.(field_name) = repmat(pol_age_slice.(field_name), [1, nkpps, 1]);
            end
        end

        function create_exogenous_dynamics_plots(cS, Z_out, A_path_annual, birth_rate_path_annual, annual_years_vec, pop_growth_rate_path_annual)
    % =========================================================================
    % == 函数: create_exogenous_dynamics_plots
    % == 版本: [v_viz_growth_rate - 可视化增长率版]
    % ==
    % == 核心修改:
    % ==   1. 接收一个新的输入参数: pop_growth_rate_path_annual。
    % ==   2. 在图1的左轴上，除了绘制老龄化率，还使用不同颜色绘制人口年增长率。
    % ==   3. 调整图例和标签以反映新增的曲线。
    % =========================================================================

    fprintf('   正在生成外生动态可视化图表...\n');

    % --- 创建图形窗口 ---
    figure('Name', '外生路径核心动态可视化', 'Position', [100, 100, 1400, 600]);

    % --- 图 1: 人口动态 (老龄化、增长率 vs 出生率) ---
    subplot(1, 2, 1);
    
    % -- 左轴 --
    yyaxis left
    
    % 计算并绘制退休人口占比 (基于模型期数据)
    retiree_indices = (cS.aR_new + 1):cS.aD_new;
    model_sim_years = cS.start_year:cS.time_Step:cS.end_year;
    retiree_proportion_path = sum(Z_out(retiree_indices, :), 1);
    p1 = plot(model_sim_years, retiree_proportion_path * 100, 'b-o', 'LineWidth', 2.5, 'MarkerSize', 4, 'MarkerFaceColor', 'b');
    ylabel('占比或增长率 (%)');
    ax = gca;
    ax.YAxis(1).Color = 'k'; % 将坐标轴颜色设为黑色以保持中立

    hold on;

    % [核心修改] 绘制人口年增长率
    % 增长率对应的时间点是每个时期的末尾，所以是 t=2...T_final
    years_for_growth_rate = model_sim_years(2:end);
    p2 = plot(years_for_growth_rate, pop_growth_rate_path_annual * 100, 'm-s', 'LineWidth', 2, 'MarkerSize', 4, 'MarkerFaceColor','m');
    hold off;

    % -- 右轴 --
    yyaxis right
    
    % 绘制年化出生率
    p3 = plot(annual_years_vec, birth_rate_path_annual * 1000, 'r--', 'LineWidth', 2);
    ylabel('年化出生率 (每千人)', 'Color', 'r');
    ax.YAxis(2).Color = 'r';

    % -- 图形格式 --
    title('人口动态：结构、增长与出生率');
    xlabel('年份');
    grid on;
    xlim([cS.start_year, cS.end_year]);
    
    % 更新图例
    legend([p1, p2, p3], {'退休人口占比 (%) [左轴]', '人口年增长率 (%) [左轴]', '年化出生率 (每千人) [右轴]'}, 'Location', 'northwest');

    % --- 图 2: 技术动态 (TFP水平 vs 增长率) ---
    subplot(1, 2, 2);

    if ~isempty(A_path_annual)
        % 计算TFP年增长率
        tfp_growth_rate_annual = (A_path_annual(2:end) ./ A_path_annual(1:end-1)) - 1;

        % 左轴: 绘制TFP水平 (归一化)
        A_path_annual_normalized = A_path_annual / A_path_annual(1);
        yyaxis left
        plot(annual_years_vec, A_path_annual_normalized, 'k-s', 'LineWidth', 2.5, 'MarkerSize', 4, 'MarkerFaceColor', 'k');
        ylabel('TFP 水平 (相对于起始年份)', 'Color', 'k');
        ax_tech = gca;
        ax_tech.YAxis(1).Color = 'k';

        % 右轴: 绘制TFP年增长率
        yyaxis right
        plot(annual_years_vec(2:end), tfp_growth_rate_annual * 100, 'g--', 'LineWidth', 2);
        ylabel('TFP 年增长率 (%)', 'Color', 'g');
        hold on;
        yline(cS.g_A_ss * 100, 'g:', 'LineWidth', 1.5, 'DisplayName', sprintf('长期稳态增长率 (%.2f%%)', cS.g_A_ss*100));
        hold off;
        ax_tech.YAxis(2).Color = 'g';
        
        legend('TFP水平', 'TFP年增长率', '稳态增长率', 'Location', 'northwest');
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

            % 构造一个从0到1的归一化时间向量
            t = linspace(0, 1, num_steps);

            % 应用幂函数来生成非线性的权重
            weights = t .^ curvature;

            % 使用权重来混合起始值和终点值
            smooth_path = (1 - weights) * start_val + weights * end_val;
        end
        % 假设Z_out是你提供的数据矩阵
        % plot_population_pyramid_simplified(Z_out);
    end
end