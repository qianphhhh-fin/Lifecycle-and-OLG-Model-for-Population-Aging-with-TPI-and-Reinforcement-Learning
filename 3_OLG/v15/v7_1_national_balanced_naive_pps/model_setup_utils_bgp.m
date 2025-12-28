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

            % --- [BGP修改] 核心增长参数 ---
            cS.g_A_ss = 0.015;      % 长期技术年增长率，是稳态化改造的引擎

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

        function [Z_out, A_path] = load_exogenous_paths(cS)
    %=====================================================================================
    %== 函数说明: [v16 - NaN BUG最终修正版]
    %==
    %== 核心更新:
    %==   1. [致命BUG修正] 修正了 interp1 函数调用时的参数错误。之前由于
    %==      缺少关键的Y坐标向量(birth_rate_path_annual)，导致插值结果
    %==      为NaN，进而污染了整个人口路径。此版本已恢复正确的函数调用。
    %=====================================================================================

    % --- 0. 参数处理 ---    
    fprintf('--- [v16] 外生路径生成器 (NaN BUG最终修正版) ---\n');

    % --- 1. [全局设定] 加载并确定唯一的长期稳定出生率 (这部分不变) ---
    try
        cbr_data = readtable('data\人口\UN_PPP2024_CBR_birthper1000_China.xlsx');
    catch ME
        error('无法加载出生率数据文件: %s', ME.message);
    end
    var_names = cbr_data.Properties.VariableNames;
    year_cols_indices = find(startsWith(var_names, 'y'));
    un_years = str2double(cellfun(@(x) x(2:end), var_names(year_cols_indices), 'UniformOutput', false));
    un_cbr = table2array(cbr_data(1, year_cols_indices)) / 1000; % 千分率转为比率
    final_years_mask = un_years >= 2050 & un_years <= 2100;
    if ~any(final_years_mask)
        final_birth_rate_annual = un_cbr(end);
        warning('未找到2080-2100年的出生率数据，使用最后一个可用数据作为长期出生率。');
    else
        final_birth_rate_annual = mean(un_cbr(final_years_mask));
    end
    fprintf('   已确定全局长期年化出生率 (CBR): %.5f (%.2f / 1000)\n', final_birth_rate_annual, final_birth_rate_annual*1000);
    
    % --- 2. 模拟人口过渡路径 ---
    model_sim_years = cS.start_year:cS.time_Step:cS.end_year;
    fprintf('   正在基于【纯年龄组】和【正确总人口】逻辑模拟人口过渡路径...\n');
    
    % --- 步骤 2.1: 加载初始人口数据 (这部分不变) ---
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
    
    % --- 步骤 2.2: 准备出生率路径 (这部分不变) ---
    first_needed_birth_year = min(model_sim_years);
    last_needed_birth_year = max(model_sim_years);
    annual_birth_years_vec = first_needed_birth_year:last_needed_birth_year;
    birth_rate_path_annual = interp1(un_years, un_cbr, annual_birth_years_vec, 'linear', 'extrap');
    last_un_year_for_birth = un_years(end);
    constant_part_mask = annual_birth_years_vec > last_un_year_for_birth;
    birth_rate_path_annual(constant_part_mask) = final_birth_rate_annual;
    fprintf('   出生率路径在 %d 年后将固定为长期稳态值以保证收敛。\n', last_un_year_for_birth);

    % --- 步骤 2.3: 直接在年龄组层面进行迭代预测 ---
    Z_path_raw = zeros(cS.aD_new, cS.T_sim);
    Z_path_raw(:, 1) = initial_pop_model_group;
    
    for t = 1:(cS.T_sim - 1)
        % 步骤A: 现有所有年龄组的人口按存活率进入下一个年龄组
        Z_path_raw(2:end, t+1) = Z_path_raw(1:(cS.aD_new-1), t) .* cS.s_pathV(1:(cS.aD_new-1));

        % 步骤B: 直接使用【模型内】当前的总人口作为计算基数
        total_pop_base_for_birth = sum(Z_path_raw(:, t));
        
        % 步骤C: 计算新进入模型的人口 (第一年龄组)
        current_year = model_sim_years(t+1);
        
        % [!!!!! 致命BUG修正 !!!!!] 修正 interp1 调用参数
        % 正确语法: interp1(已知X, 已知Y, 查询X)
        birth_rate_hist = interp1(annual_birth_years_vec, birth_rate_path_annual, current_year, 'linear');
        
        % 新生儿数量 = 总人口基数 * 年化出生率 * 模型期长度
        new_entrants = total_pop_base_for_birth * birth_rate_hist * cS.time_Step;
        
        Z_path_raw(1, t+1) = new_entrants;
    end
    
    % --- 步骤 2.4: 归一化并生成TFP路径 ---
    Z_out = Z_path_raw ./ sum(Z_path_raw, 1);

    % 生成TFP路径
    annual_years_vec_for_tfp = cS.start_year:cS.end_year;
    A_path = model_setup_utils_bgp.generate_tfp_path(cS, annual_years_vec_for_tfp, model_sim_years);

    % --- 可视化部分 (保持不变) ---
    % fprintf('   正在生成补充人口动态可视化图表...\n');
    % figure('Name', '人口结构与总量动态演变', 'Position', [150, 150, 1200, 500]);
    % subplot(1, 2, 1);
    % retiree_indices = (cS.aR_new + 1):cS.aD_new;
    % retiree_proportion_path = sum(Z_out(retiree_indices, :), 1);
    % plot(model_sim_years, retiree_proportion_path * 100, 'r-o', 'LineWidth', 2, 'MarkerSize', 4, 'MarkerFaceColor', 'r');
    % title('退休年龄组人口占比变化 (老龄化趋势)');
    % xlabel('年份');
    % ylabel('退休人口占比 (%)');
    % grid on;
    % xlim([cS.start_year, cS.end_year]);
    % text(model_sim_years(1), retiree_proportion_path(1)*100, sprintf('  初值: %.2f%%', retiree_proportion_path(1)*100), 'VerticalAlignment', 'bottom');
    % text(model_sim_years(end), retiree_proportion_path(end)*100, sprintf('终值: %.2f%%  ', retiree_proportion_path(end)*100), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    % subplot(1, 2, 2);
    % total_population_path = sum(Z_path_raw, 1);
    % plot(model_sim_years, total_population_path, 'b-s', 'LineWidth', 2, 'MarkerSize', 4, 'MarkerFaceColor', 'b');
    % title('总人口（模型内）总量变化');
    % xlabel('年份');
    % ylabel('人口总量 (千人)');
    % grid on;
    % xlim([cS.start_year, cS.end_year]);
    % text(model_sim_years(1), total_population_path(1), sprintf('  初值: %.0f', total_population_path(1)), 'VerticalAlignment', 'bottom');
    % text(model_sim_years(end), total_population_path(end), sprintf('终值: %.0f  ', total_population_path(end)), 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
    % sgtitle('补充人口动态可视化：老龄化趋势与总量演变', 'FontSize', 16, 'FontWeight', 'bold');
    % fprintf('✅ 补充可视化完成。\n');
end        
        
function A_path = generate_tfp_path(cS, annual_years_vec, model_sim_years)
            fprintf('   正在基于历史和预测数据插值TFP路径 (A_path)...\n');
            tfp_data_pwt = readtable('data\PWT\china_pwt_data.xlsx');
            pwt_years = tfp_data_pwt.year(tfp_data_pwt.year >= 1979);
            pwt_tfp_level = tfp_data_pwt.rtfpna(tfp_data_pwt.year >= 1979);
            bai_projections = [2025, 0.0557; 2030, 0.0482; 2035, 0.0394; 2040, 0.0340; 2045, 0.0346; 2050, 0.0298];
            long_term_g = 0.015;
            proj_years = [pwt_years(end); bai_projections(:,1)];
            proj_tfp_level = zeros(size(proj_years));
            proj_tfp_level(1) = pwt_tfp_level(end);
            g_rate_path = [0.058; bai_projections(:,2); long_term_g];
            for i = 1:(length(proj_years)-1)
                start_yr = proj_years(i); end_yr = proj_years(i+1);
                g_rate = g_rate_path(i);
                proj_tfp_level(i+1) = proj_tfp_level(i) * (1 + g_rate)^(end_yr - start_yr);
            end
            final_anchor_year = cS.end_year + 20;
            final_anchor_level = proj_tfp_level(end) * (1 + long_term_g)^(final_anchor_year - proj_years(end));
            proj_years = [proj_years; final_anchor_year];
            proj_tfp_level = [proj_tfp_level; final_anchor_level];
            all_anchor_years = [pwt_years; proj_years(2:end)];
            all_anchor_levels = [pwt_tfp_level; proj_tfp_level(2:end)];
            [unique_years, ia, ~] = unique(all_anchor_years);
            unique_levels = all_anchor_levels(ia);
            A_path_annual_level = exp(interp1(unique_years, log(unique_levels), annual_years_vec, 'pchip', 'extrap'));
            A_path = interp1(annual_years_vec, A_path_annual_level, model_sim_years, 'linear');
            A_path = A_path / A_path(1);
            fprintf('✅ TFP路径构建完成。\n');
        end

        % --- 辅助函数：创建可视化图表 (从原函数中提取) ---
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

% 假设Z_out是你提供的数据矩阵
% plot_population_pyramid_simplified(Z_out);
    end
end