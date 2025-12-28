classdef model_setup_utils
    methods (Static)

        function cS = ParameterValues()
            % [vPaper.5 - 政府投资版]
            % [修改] 提供了完整的参数列表，并将待校准参数标记为“初始猜测值”。

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

            % --- [修改] 待校准参数的初始猜测值 ---
            cS.beta =0.996;        % 主观折现因子
            cS.gamma = 0.12;        % 公共资本的产出弹性
            cS.lambda_g = 0.35;    % 政府投资占总税收的比例

            % --- 其他固定参数 ---
            cS.sigma = 8;
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
            cS.gov_debt_frac_Y = 0;
            cS.A = 1.0; % TFP基准值

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
            cS.p_shock_young_peak_annual = 0.02;
            cS.p_shock_old_peak_annual = 0.02;
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
            % [修改] 支持可选的网格上限参数，实现自适应网格
            % 新签名: generateGrids(cS, 'k_max', value, 'kpps_max', value)
            % 这将网格生成从固定参数变为动态服务
            
            % 使用 inputParser 处理可选参数
            p = inputParser;
            
            % 计算默认值
            cS.tgKY = 3; 
            cS.tgWage = (1-cS.alpha)*cS.A*((cS.tgKY/cS.A)^(cS.alpha/(1-cS.alpha)));
            default_kMax = 50 * cS.tgWage;
            default_kppsMax = default_kMax / 2;
            
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
            power_k = 4;
            
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
            % [v.PlateauShock - 平台式冲击版]
            % 核心修改: 将老年冲击的概率生成逻辑，从“斜坡式”(linspace)
            %           修改为与青年冲击完全相同的“平台式”(固定概率)。

            % --- 1. 生成基础AR(1)过程 ---
            lePersistence = 0.9; leShockStd = 0.15^0.5; Tauchen_q = 2.0;
            [leLogGridV_normal, leTrProbM_normal] = model_setup_utils.tauchen(cS.nw, lePersistence, leShockStd, 0, Tauchen_q);
            [~, D] = eig(leTrProbM_normal'); [~, c] = min(abs(diag(D)-1));
            leProb1V = abs(D(:,c)/sum(D(:,c)));

            % --- 2. 根据参数构建年度冲击概率向量 ---
            phys_ages = cS.age1_orig : cS.ageLast_orig;
            num_phys_ages = length(phys_ages);
            p_young_annual = zeros(num_phys_ages, 1);
            p_old_annual = zeros(num_phys_ages, 1);

            % 处理青年冲击 (保持不变)
            y_start_idx = find(phys_ages == cS.young_shock_start_age, 1);
            y_end_idx = find(phys_ages == cS.young_shock_end_age, 1);
            if ~isempty(y_start_idx) && ~isempty(y_end_idx)
                p_young_annual(y_start_idx:y_end_idx) = cS.p_shock_young_peak_annual;
            end

            % [!!!!! 核心修改 !!!!!]
            % 处理老年冲击 (应用与青年冲击相同的“平台式”逻辑)
            o_start_idx = find(phys_ages == cS.old_shock_start_age, 1);
            o_end_idx = find(phys_ages == cS.old_shock_end_age, 1);
            if ~isempty(o_start_idx) && ~isempty(o_end_idx)
                p_old_annual(o_start_idx:o_end_idx) = cS.p_shock_old_peak_annual;
            end
            % [!!!!! 修改结束 !!!!!]

            % --- 3. 从年度概率聚合为模型期概率 ---
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

            % --- 4. 构建最终的年龄依赖转移矩阵并汇报 ---
            leGridV_expanded = [exp(leLogGridV_normal(:)); 0; 0];
            nw_expanded = cS.nw + 2;
            TrProbM_by_age = cell(cS.aD_new, 1);

            for a = 1:cS.aD_new
                Tr_a = zeros(nw_expanded, nw_expanded);
                p_y = p_shock_path_model(a, 1);
                p_o = p_shock_path_model(a, 2);
                p_total_shock = min(1.0, p_y + p_o);
                if p_total_shock > 0, p_y_norm = p_y/p_total_shock; p_o_norm = p_o/p_total_shock; else, p_y_norm=0; p_o_norm=0; end

                Tr_a(1:cS.nw, 1:cS.nw) = leTrProbM_normal * (1 - p_total_shock);
                Tr_a(1:cS.nw, cS.nw + 1) = p_total_shock * p_y_norm;
                Tr_a(1:cS.nw, cS.nw + 2) = p_total_shock * p_o_norm;
                Tr_a(cS.nw + 1, 1:cS.nw) = leProb1V';
                Tr_a(cS.nw + 2, 1:cS.nw) = leProb1V';

                Tr_a = Tr_a ./ sum(Tr_a, 2);
                TrProbM_by_age{a} = Tr_a;

                % fprintf('\n================== 模型年龄组 a = %d ==================\n', a);
                % fprintf('物理年龄范围: %d-%d岁\n', cS.age1_orig + (a-1)*cS.time_Step, cS.age1_orig + a*cS.time_Step -1);
                % fprintf('本期5年冲击概率: p(青年)=%.4f, p(老年)=%.4f, p(总)=%.4f\n', p_y, p_o, p_total_shock);
            end
        end

        function [y_grid_out, trProbM_out] = tauchen(N, rho, sigma, mu, m)
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
            c_adj = max(cS.cFloor, cM);
            if abs(sigma - 1) < 1e-6
                utilM = log(c_adj);
                muM = 1./c_adj;
            else utilM = (c_adj.^(1-sigma))./(1-sigma);
                muM = c_adj.^(-sigma);
            end
            % utilM(cM < cS.cFloor) = -1e10 - (cS.cFloor - cM(cM < cS.cFloor))*1e10;
        end

        function util_beq = bequest_utility(k_prime, cS)
            % =========================================================================
            % == FUNCTION: bequest_utility
            % == 目的: 计算给定遗赠水平 k_prime 的“暖光”效用。
            % == 形式: 与消费效用函数形式相同，但由强度参数 cS.phi_bequest 调节。
            % =========================================================================
            if cS.phi_bequest <= 0
                util_beq = 0;
                return;
            end

            % 我们使用与消费相同的 sigma，但引入 phi_bequest 作为权重
            % 同样使用 cFloor 来避免 log(0) 或负数次幂等数学问题
            k_adj = max(cS.cFloor, k_prime);

            if abs(cS.sigma - 1) < 1e-6
                util_beq = cS.phi_bequest * log(k_adj);
            else
                util_beq = cS.phi_bequest * (k_adj.^(1-cS.sigma))./(1-cS.sigma);
            end
        end


        function [Z_path, A_path] = load_exogenous_paths(cS, plot_flag)
            % =========================================================================
            % == FUNCTION: load_exogenous_paths (v4 - Self-Contained Final Version)
            % == 目的: 集中处理所有外部数据加载和路径生成。
            % == 核心功能:
            % == 1. 自动加载并处理联合国出生率预测数据。
            % == 2. 只读取初始年份的人口分布作为起点。
            % == 3. 基于模型内生存活率和外部出生率预测，向前模拟人口路径。
            % == 4. 基于PWT和预测数据，插值生成TFP路径。
            % == 5. 可选的可视化功能。
            % =========================================================================
            if nargin < 2, plot_flag = false; end % 默认不画图
            fprintf('--- [v4] 集中加载和处理所有外生路径 ---\n');

            % 定义模拟所需的年份向量 (每年) 和期数
            annual_years_vec = cS.start_year:cS.end_year;
            model_sim_years = cS.start_year:cS.time_Step:cS.end_year;
            T_sim = length(model_sim_years);

            %% 1. 加载并处理联合国出生率数据 (CBR)
            fprintf('   加载联合国出生率预测数据...\n');
            try
                cbr_data = readtable('data\人口\UN_PPP2024_CBR_birthper1000_China.xlsx');
            catch ME
                error('无法加载出生率数据文件。请确保文件路径正确: data\\人口\\UN_PPP2024_CBR_birthper1000_China.xlsx. 错误信息: %s', ME.message);
            end
            
            % 提取年份和数据
            var_names = cbr_data.Properties.VariableNames;
            year_cols_indices = find(startsWith(var_names, 'y')); % 年份列以'y'开头, e.g., 'y2024'
            
            if isempty(year_cols_indices)
                error('未找到年份列。请确保Excel文件包含以y开头的年份列（如y2024）。');
            end
            
            % 提取年份（去掉'y'前缀）
            un_years = str2double(cellfun(@(x) x(2:end), var_names(year_cols_indices), 'UniformOutput', false));
            
            % 提取出生率数据（假设只有一行中国数据）
            un_cbr_permille = table2array(cbr_data(1, year_cols_indices));
            
            % 从千分之转换为比率
            un_cbr = un_cbr_permille / 1000;

            % 使用插值为模拟期间的每一年生成出生率路径
            % [修正] 对于超出数据范围的年份，保持为最后一年的数据值，而不是外推
            birth_rate_path_annual = zeros(size(annual_years_vec));
            
            % 找出数据范围内和超出范围的年份
            max_data_year = max(un_years);
            min_data_year = min(un_years);
            
            for i = 1:length(annual_years_vec)
                year_i = annual_years_vec(i);
                if year_i <= max_data_year && year_i >= min_data_year
                    % 在数据范围内，使用线性插值
                    birth_rate_path_annual(i) = interp1(un_years, un_cbr, year_i, 'linear');
                elseif year_i > max_data_year
                    % 超出最大年份，使用最后一年的数据
                    birth_rate_path_annual(i) = un_cbr(end);
                else
                    % 小于最小年份，使用第一年的数据
                    birth_rate_path_annual(i) = un_cbr(1);
                end
            end
            
            fprintf('✅ 出生率路径处理完成。\n');

            %% 2. 生成模型内生的人口路径 (Z_path)
            fprintf('   正在基于初始分布和出生率模拟人口路径 (Z_path)...\n');
            
            % --- 步骤 2.1: 加载初始年份的人口分布 ---
            try
                pop_data = readtable('data\人口\population_by_age_group_all_years.xlsx', 'Sheet', 'population');
            catch ME
                error('无法加载人口分布数据文件。请确保文件路径正确: data\\人口\\population_by_age_group_all_years.xlsx. 错误信息: %s', ME.message);
            end
            
            % [修正] 使用ss0_year作为初始人口分布，因为这是实际有数据的年份
            initial_year_column_name = ['y', num2str(cS.ss0_year)];
            if ~ismember(initial_year_column_name, pop_data.Properties.VariableNames)
                error('在Excel文件中找不到指定的初始年份列: %s', initial_year_column_name);
            end
            initial_pop_dist = pop_data.(initial_year_column_name);
            
            % --- 步骤 2.2: 基于内生存活率和外生出生率，向前模拟年度人口分布 ---
            % [修正] 人口模拟需要从ss0_year开始，但annual_years_vec从start_year开始
            % 所以需要扩展一年来包含初始年份
            extended_years_vec = [cS.ss0_year, annual_years_vec];
            extended_birth_rate = [birth_rate_path_annual(1), birth_rate_path_annual]; % 假设ss0_year的出生率与start_year相同
            
            Pop_dist_annual = zeros(cS.aD_new, length(extended_years_vec));
            Pop_dist_annual(:, 1) = initial_pop_dist; % ss0_year的人口分布
            
            s_pathV_annual = cS.s_pathV.^(1/cS.time_Step);
            
            for t = 1:(length(extended_years_vec) - 1)
                total_pop_t = sum(Pop_dist_annual(:, t));
                newborns_t_plus_1 = total_pop_t * extended_birth_rate(t);
                
                Pop_dist_annual(1, t+1) = max(0, newborns_t_plus_1);
                Pop_dist_annual(2:end, t+1) = Pop_dist_annual(1:end-1, t) .* s_pathV_annual(1:end-1);
            end
            
            % 去掉初始年份，只保留过渡路径年份
            Pop_dist_annual = Pop_dist_annual(:, 2:end);

            % --- 步骤 2.3: 从年度路径中采样并归一化 ---
            Z_path_raw = interp1(annual_years_vec, Pop_dist_annual', model_sim_years, 'linear')';
            Z_path = Z_path_raw ./ sum(Z_path_raw, 1);
            
            fprintf('✅ 人口路径模拟完成 (%d-%d, %d个时期)。\n', cS.start_year, cS.end_year, T_sim);
            fprintf('   初始人口分布基于 %d 年数据。\n', cS.ss0_year);

            %% 3. 生成基于插值的TFP路径 (A_path)
            fprintf('   正在基于历史和预测数据插值TFP路径 (A_path)...\n');
            
            tfp_data_pwt = readtable('data\PWT\china_pwt_data.xlsx');
            pwt_years = tfp_data_pwt.year(tfp_data_pwt.year >= 1979);
            pwt_tfp_level = tfp_data_pwt.rtfpna(tfp_data_pwt.year >= 1979);
            
            bai_projections = [2025, 0.0557; 2030, 0.0482; 2035, 0.0394; 2040, 0.0340; 2045, 0.0346; 2050, 0.0298];
            long_term_g = 0.015;

            proj_years = [pwt_years(end); bai_projections(:,1)];
            proj_tfp_level = zeros(size(proj_years));
            proj_tfp_level(1) = pwt_tfp_level(end);
            
            g_rate_path = [0.058; bai_projections(:,2); long_term_g]; % 包含2020-2024的增长率

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
            
            fprintf('✅ TFP路径构建完成 (%d-%d)。\n', cS.start_year, cS.end_year);
            
            %% 4. 可选的可视化
            if plot_flag
                fprintf('   正在生成外生路径可视化图表...\n');
                
                % 创建年份向量用于绘图
                plot_years = model_sim_years;
                
                figure('Name', '外生路径可视化', 'Position', [100, 100, 1200, 500]);
                
                % 子图1: 人口分布路径 (Z_path)
                subplot(1, 2, 1);
                
                % 选择几个代表性年龄组进行可视化
                selected_ages = [1, 4, 8, 12, 16]; % 对应大约25岁、40岁、55岁、70岁、85岁
                age_labels = {};
                for i = 1:length(selected_ages)
                    age_start = cS.age1_orig + (selected_ages(i) - 1) * cS.time_Step;
                    age_end = age_start + cS.time_Step - 1;
                    age_labels{i} = sprintf('%d-%d岁', age_start, age_end);
                end
                
                % 绘制选定年龄组的人口占比变化
                colors = lines(length(selected_ages));
                for i = 1:length(selected_ages)
                    plot(plot_years, Z_path(selected_ages(i), :) * 100, 'Color', colors(i,:), 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 4);
                    hold on;
                end
                
                xlabel('年份');
                ylabel('人口占比 (%)');
                title('人口分布路径 (Z\_path)');
                legend(age_labels, 'Location', 'best', 'FontSize', 8);
                grid on;
                xlim([plot_years(1), plot_years(end)]);
                
                % 子图2: TFP路径 (A_path)
                subplot(1, 2, 2);
                plot(plot_years, A_path, 'b-', 'LineWidth', 2.5, 'Marker', 's', 'MarkerSize', 5, 'MarkerFaceColor', 'b');
                xlabel('年份');
                ylabel('TFP水平 (相对于起始年份)');
                title('全要素生产率路径 (A\_path)');
                grid on;
                xlim([plot_years(1), plot_years(end)]);
                
                % 添加一些统计信息
                avg_growth_rate = (A_path(end)/A_path(1))^(1/(T_sim*cS.time_Step)) - 1;
                text(0.05, 0.95, sprintf('平均年增长率: %.2f%%', avg_growth_rate*100), ...
                     'Units', 'normalized', 'FontSize', 10, 'BackgroundColor', 'white');
                
                sgtitle(sprintf('外生路径概览 (%d-%d年, %d个时期)', cS.start_year, cS.end_year, T_sim), 'FontSize', 14);
                
                fprintf('✅ 外生路径可视化完成。\n');
            end
        end


        function Z_ss0 = get_calibration_inputs(target_year, cS)
            % [新增] 获取校准年份所需的外生输入
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
            % [v14.4 - 已修正插值逻辑]
            % 核心修正：从插值【绝对人数】改为直接插值【覆盖率】，以保证路径的平滑和目标的准确达成。

            fprintf('正在构建基于【覆盖率插值】的有效养老金缴费率路径 (theta_path)...\n');

            % --- 1. 定义关键参数 ---
            theta_urban_employee_effective = 0.20; % 城镇职工体系的名义缴费率
            theta_resident_effective = 0.03;       % 城乡居民体系的名义缴费率 (近似值)

            % [政策目标] 定义未来的目标覆盖率
            coverage_urban_final = 0.8;    % 目标：城镇职工覆盖率达到劳动人口的60%
            coverage_resident_final = 0.8; % 目标：城乡居民覆盖率达到劳动人口的35%
            year_reach_final = 2050;        % 假设在2050年达到并维持这个目标覆盖率

            % --- 2. 收集来自官方统计的真实数据点 (单位：万人) ---
            year_pax_urban = [1997, 2000, 2005, 2008, 2011, 2014, 2017, 2020, 2023];
            pax_urban      = [8671, 10448, 17487, 21890, 28391, 34124, 40293, 45638, 52112];

            year_pax_resident = [2008, 2011, 2014, 2017, 2020, 2023];
            pax_resident      = [0,    45586, 50075, 51255, 54227, 54475];

            year_laborpop_data = [1997, 2000, 2005, 2008, 2011, 2014, 2017, 2020, 2023];
            laborpop_data      = [78531, 82352, 87532, 90833, 94072, 92837, 91199, 89552, 86481];

            % --- 3. [核心修正] 计算历史数据点的【覆盖率】 ---
            % a. 计算城镇职工历史覆盖率
            coverage_urban_data = pax_urban ./ laborpop_data;

            % b. 计算城乡居民历史覆盖率
            %    需要先插值得到对应年份的劳动人口
            laborpop_for_resident_years = interp1(year_laborpop_data, laborpop_data, year_pax_resident, 'linear', 'extrap');
            coverage_resident_data = pax_resident ./ laborpop_for_resident_years;

            % --- 4. [核心修正] 直接对【覆盖率】进行插值，生成完整的年度路径 ---
            annual_years_vector = cS.start_year:cS.end_year;

            % a. 插值生成年度【城镇职工覆盖率】路径
            interp_years_urban = [year_pax_urban, year_reach_final, cS.end_year];
            interp_coverage_urban = [coverage_urban_data, coverage_urban_final, coverage_urban_final]; % 达到目标后保持不变
            [unique_years_u, ia_u, ~] = unique(interp_years_urban);
            coverage_urban_annual = interp1(unique_years_u, interp_coverage_urban(ia_u), annual_years_vector, 'linear');

            % b. 插值生成年度【城乡居民覆盖率】路径
            interp_years_resident = [cS.start_year, year_pax_resident, year_reach_final, cS.end_year];
            interp_coverage_resident = [0, coverage_resident_data, coverage_resident_final, coverage_resident_final]; % 达到目标后保持不变
            [unique_years_r, ia_r, ~] = unique(interp_years_resident);
            coverage_resident_annual = interp1(unique_years_r, interp_coverage_resident(ia_r), annual_years_vector, 'linear');
            coverage_resident_annual(annual_years_vector < min(year_pax_resident)) = 0; % 确保制度开始前为0

            % --- 5. 基于【覆盖率路径】计算最终的【年度有效缴费率】路径 ---
            theta_path_annual = (coverage_urban_annual * theta_urban_employee_effective) + ...
                (coverage_resident_annual * theta_resident_effective);

            % --- 6. 从年度路径中，提取出模型【5年期】所需的路径 ---
            model_year_indices = 1:(cS.end_year - cS.start_year + 1);
            model_year_indices_5yr = model_year_indices(1:cS.time_Step:end);

            T_sim = length(cS.start_year:cS.time_Step:cS.end_year);
            theta_path = theta_path_annual(model_year_indices_5yr);
            theta_path = theta_path(1:min(T_sim, length(theta_path))); % 确保长度正确

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
                ylim([0, max(cS.theta_path)*1.2]); % 调整Y轴范围
            end
            fprintf('✅ 基于覆盖率插值的缴费率路径构建完成。\n');
        end
        


    end
end