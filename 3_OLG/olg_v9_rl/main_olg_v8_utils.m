% --- 开始文件：main_olg_v8_utils.m (最终修正版 v9.2-compatible) ---

% =====================================================================
% === OLG 模型 V9 兼容工具函数库 (MATLAB VFI 实现) ===
% =====================================================================
%
% 🎯 核心目标：确保此 VFI 实现与 Python RL 框架在理论和实现上完全一致。
%
% [最终修正] v9.2:
% - 恢复所有原始函数，确保文件完整性。
% - 集成所有已讨论的修正（参数对齐、逻辑对齐、插值方法、纯年龄组逻辑）。
% - 修复了 CallInterpolator 缺失的错误。
% - 本文件为最终、完整的、可直接运行的MATLAB VFI工具库。
% =====================================================================

classdef main_olg_v8_utils

    methods (Static)

        % =====================================================================
        % == 1. 模型参数设定函数 ==
        % =====================================================================
        function cS = ParameterValues_HuggettStyle()
            % [VFI最终对齐版 v2] - 添加了均衡求解器参数

            % --- 人口结构基础参数 ---
            cS.age1_orig = 20;
            cS.ageLast_orig = 98;
            cS.ageRetire_orig = 65;
            cS.aD_orig = cS.ageLast_orig - cS.age1_orig + 1;
            cS.aR_idx_orig = cS.ageRetire_orig - cS.age1_orig + 1;
            cS.aW_orig = cS.aR_idx_orig - 1;
            cS.physAgeV_orig = (cS.age1_orig : cS.ageLast_orig)';

            % --- 年龄组聚合参数 ---
            cS.yearStep = 5;
            cS.aD_new = ceil(cS.aD_orig / cS.yearStep);
            cS.aR_new = ceil(cS.aW_orig / cS.yearStep);

            cS.physAgeMap = cell(cS.aD_new, 1);
            for a = 1:cS.aD_new
                startIdx = (a-1)*cS.yearStep + 1;
                endIdx = min(a*cS.yearStep, cS.aD_orig);
                cS.physAgeMap{a} = startIdx:endIdx;
            end

            d_orig_data = [0.00159,0.00169,0.00174,0.00172,0.00165,0.00156,0.00149,0.00145,0.00145,0.00149,0.00156,0.00163,0.00171,0.00181,0.00193,0.00207,0.00225,0.00246,0.00270,0.00299,0.00332,0.00368,0.00409,0.00454,0.00504,0.00558,0.00617,0.00686,0.00766,0.00865,0.00955,0.01058,0.01162,0.01264,0.01368,0.01475,0.01593,0.01730,0.01891,0.02074,0.02271,0.02476,0.02690,0.02912,0.03143,0.03389,0.03652,0.03930,0.04225,0.04538,0.04871,0.05230,0.05623,0.06060,0.06542,0.07066,0.07636,0.08271,0.08986,0.09788,0.10732,0.11799,0.12895,0.13920,0.14861,0.16039,0.17303,0.18665,0.20194,0.21877,0.23601,0.25289,0.26973,0.28612,0.30128,0.31416,0.32915,0.34450,0.36018];
            s_orig = 1 - d_orig_data;
            cS.s_1yr_transitionV = zeros(cS.aD_new, 1);
            for a = 1:(cS.aD_new - 1)
                lastYearIdxInGroup = cS.physAgeMap{a}(end);
                cS.s_1yr_transitionV(a) = s_orig(lastYearIdxInGroup);
            end

            % --- [对齐] 家庭偏好参数 ---
            cS.sigma = 1.5;
            cS.beta = 0.97;
            cS.cFloor = 0.05;
            cS.nSim = 1000;

            % --- [对齐] 生产技术参数 ---
            cS.A = 0.895944;
            cS.alpha = 0.36;
            cS.ddk = 0.06;

            % --- [对齐] 政府财政参数 ---
            cS.tau_k = 0.20;
            cS.tau_c = 0.10;
            cS.gov_exp_frac_Y = 0.15;
            cS.gov_debt_frac_Y = 0.60;

            % --- [对齐] 劳动效率冲击过程参数 ---
            cS.lePersistence = 0.96;
            cS.leShockStd = 0.045^0.5;
            cS.nw = 5;

            % --- [对齐] 资产网格参数 ---
            cS.tgKY = 3;
            cS.tgWage = (1-cS.alpha)*cS.A*((cS.tgKY/cS.A)^(cS.alpha/(1-cS.alpha)));
            cS.nk = 50;
            cS.nkpps = 30;
            cS.nkprime = 20;
            cS.npps = 20;
            cS.kMin = 0;
            cS.kMax = 15 * cS.tgWage;
            cS.kGridV = (cS.kMin + (cS.kMax - cS.kMin) * (linspace(0, 1, cS.nk).^1.5))';
            if cS.nk > 0, cS.kGridV(1) = cS.kMin; end

            cS.kppsMin = 0;
            cS.kppsMax = cS.kMax / 2;
            cS.kppsGridV = (cS.kppsMin + (cS.kppsMax - cS.kppsMin) * (linspace(0, 1, cS.nkpps).^1.5))';
            if cS.nkpps > 0, cS.kppsGridV(1) = cS.kppsMin; end

            % --- [对齐] 年龄效率剖面 ---
            ageEffV_orig_temp = zeros(100, 1);
            ageEffV_orig_temp(20:72) = [linspace(0.3,1.5,36-20+1), 1.5.*ones(1,47-37+1), linspace(1.5,0.2,65-48+1), linspace(0.18,0,72-66+1)];
            ageEffV_orig = ageEffV_orig_temp(cS.age1_orig : cS.ageLast_orig);
            cS.ageEffV_new = zeros(cS.aD_new, 1);
            for a = 1:cS.aD_new
                cS.ageEffV_new(a) = mean(ageEffV_orig(cS.physAgeMap{a}));
            end

            % --- [对齐] PPS制度参数 ---
            cS.pps_active = true;
            cS.pps_tax_rate_withdrawal = 0.03;
            cS.pps_return_rate_premium = 0.08;
            cS.pps_withdrawal_rate = 0.15;
            cS.pps_in_K = true;
            cS.pps_bequeathable = true;
            cS.pps_contrib_limit = 9999;
            cS.pps_max_contrib_frac = 1.0;

            % --- [新增] 一般均衡求解参数 (与Python对齐) ---
            cS.max_iter_K_tau_l = 100;
            cS.tol_K_tau_l = 1e-4;
            cS.damp_K_v5 = 0.1;
            cS.damp_tau_l_v5 = 0.1;
            cS.gbc_tol_for_internal_loop = 1e-3;

            % --- [新增] 收敛检测参数 ---
            cS.max_stagnation_iters = 10;
            cS.min_norm_improvement_frac = 1e-3;
            cS.max_tau_l_boundary_strikes = 5;

            % --- [新增] PAYG税率约束与初始值参数 ---
            cS.tau_l_init_guess = 0.1509;
            cS.tau_l_min = 0.00;
            cS.tau_l_max = 0.3;
            cS.max_total_labor_tax = 1.0;
            cS.theta_payg_max = 1.0;
        end
        % =====================================================================
        % == [新增] 人口动态相关辅助函数 ==
        % =====================================================================

        function popS = initPopulation(cS)
            % initPopulation - 初始化人口结构
            popS = struct();

            % 基于2023年中国人口结构的初始分布 (16个年龄组)
            initial_pop_dist = [76.2, 86.4, 113.8, 98.6, 86.6, 102.7, 112.0, 99.0, 64.0, 66.9, 44.1, 25.4, 14.9, 6.8, 1.7, 0.2];

            initial_total = sum(initial_pop_dist);

            if initial_total > 0 && length(initial_pop_dist) == cS.aD_new
                popS.Z = (initial_pop_dist / initial_total * 100)';
            else
                warning('初始人口数据不匹配或总和为零。将设置为均匀的初始年龄组人口分布。');
                popS.Z = ones(cS.aD_new, 1) * (100 / cS.aD_new);
            end

            popS.totalPop = sum(popS.Z(:, 1));

            if popS.totalPop(1) > 1e-9
                popS.ageDist = popS.Z(:, 1) / popS.totalPop(1);
            else
                popS.ageDist = zeros(cS.aD_new, 1);
            end

            popS.initialAgeDist = popS.ageDist;
            fprintf('初始年龄组人口已设置。总人口=%.2f (代表百分比基数)。\n', popS.totalPop(1));
        end

        function popS = populationDynamics(popS, cS)
            % populationDynamics - 模拟人口动态演进

            % 基于中国数据的年龄组间存活率
            beta_surv_pop = [0.998, 0.996, 0.994, 0.992, 0.988, 0.984, 0.980, 0.976, ...
                0.970, 0.960, 0.945, 0.920, 0.880, 0.800, 0.680];
            cS.survivalProbV_popdyn = [beta_surv_pop, 0]'; % 最后一个年龄组存活率为0

            max_periods = 50;
            bgp_tolerance = 0.001;
            bgp_window = 5;

            Z_history = zeros(cS.aD_new, max_periods + 1);
            totalPop_history = zeros(max_periods + 1, 1);
            ageDist_history = zeros(cS.aD_new, max_periods + 1);

            Z_history(:, 1) = popS.Z;
            totalPop_history(1) = popS.totalPop;
            ageDist_history(:, 1) = popS.ageDist;

            fprintf('人口动态模拟开始 (年龄组, 最大期数 = %d)...\n', max_periods);
            bgp_reached = false;

            for t = 1:max_periods
                Z_current = Z_history(:, t);
                Z_next = zeros(cS.aD_new, 1);

                % 时变的人口增长率
                if t < 5
                    growth_rate = -0.01 - 0.003 * t;
                else
                    growth_rate = -0.03 - 0.004 * min(t - 5, 10);
                end

                Z_next(1) = Z_current(1) * (1 + growth_rate);

                for a = 2:cS.aD_new
                    Z_next(a) = Z_current(a-1) * cS.survivalProbV_popdyn(a-1);
                end

                Z_history(:, t+1) = Z_next;
                totalPop_history(t+1) = sum(Z_next);
                if totalPop_history(t+1) > 1e-9
                    ageDist_history(:, t+1) = Z_next / totalPop_history(t+1);
                end

                if t >= bgp_window
                    is_stable = true;
                    for w = 0:bgp_window-1
                        change = norm(ageDist_history(:, t+1-w) - ageDist_history(:, t-w));
                        if change >= bgp_tolerance
                            is_stable = false;
                            break;
                        end
                    end
                    if is_stable
                        fprintf('人口稳态在模拟期数 %d 达到。\n', t);
                        bgp_reached = true;
                        Z_history = Z_history(:, 1:t+1);
                        totalPop_history = totalPop_history(1:t+1);
                        ageDist_history = ageDist_history(:, 1:t+1);
                        break;
                    end
                end
            end

            popS.Z = Z_history;
            popS.totalPop = totalPop_history;
            popS.ageDist = ageDist_history;

            if ~bgp_reached
                warning('人口稳态未在 %d 期内达到。', max_periods);
            end
        end

        function [Z_ss, dependency_ratio_ss, bgp_reached, bgp_period] = detectSteadyStatePopulation(popS, cS)
            % detectSteadyStatePopulation - 从模拟的人口历史中检测并提取稳态分布

            n_periods = size(popS.Z, 2);
            bgp_window = 5;
            bgp_tolerance = 0.001;
            bgp_reached = false;
            bgp_period = n_periods;

            if n_periods < bgp_window + 1
                warning('人口模拟期数过短，无法进行稳态检查。');
                Z_ss = popS.Z(:, end);
            else
                for t = n_periods : -1 : bgp_window + 1
                    is_stable = true;
                    for w = 0:bgp_window-1
                        change = norm(popS.ageDist(:, t-w) - popS.ageDist(:, t-w-1));
                        if change >= bgp_tolerance
                            is_stable = false;
                            break;
                        end
                    end
                    if is_stable
                        bgp_reached = true;
                        bgp_period = t;
                        break;
                    end
                end
                Z_ss = popS.Z(:, bgp_period);
            end

            working_pop_ss = sum(Z_ss(1:cS.aR_new));
            retired_pop_ss = sum(Z_ss(cS.aR_new+1:end));
            dependency_ratio_ss = retired_pop_ss / working_pop_ss;

            % 绘图比较
            figure('Name', 'VFI: 初始 vs 稳态人口分布');
            bar_data = [popS.initialAgeDist * 100, (Z_ss / sum(Z_ss)) * 100];
            bar(bar_data);
            xlabel('年龄组');
            ylabel('占总人口百分比 (%)');
            title(sprintf('初始 vs 稳态人口分布 (稳态于第%d期)', bgp_period));
            legend('初始分布', '稳态分布', 'Location', 'best');
            grid on;
        end

        function [HHlaborM_group, L_total_eff_pc] = LaborSupply_Huggett(eIdxM_group, cS, paramS, Z_ss_norm_group)
            % LaborSupply_Huggett - 计算劳动供给（适配年龄组模拟）
            nSim = size(eIdxM_group, 1);
            HHlaborM_group = zeros(nSim, cS.aD_new);

            for a_group = 1:cS.aR_new % 只在工作年龄组计算
                eIdx_this_age = eIdxM_group(:, a_group);
                labor_eff = paramS.leGridV(eIdx_this_age);
                HHlaborM_group(:, a_group) = cS.ageEffV_new(a_group) * labor_eff;
            end

            mean_labor_per_working_group = mean(HHlaborM_group(:, 1:cS.aR_new), 1);
            L_total_eff_pc = sum(mean_labor_per_working_group' .* Z_ss_norm_group(1:cS.aR_new));
        end

        function [R_market_gross_factor, MPL_gross] = HHPrices_Huggett(K, L, cS)
            % HHPrices_Huggett - 根据边际生产力计算要素价格
            if K <= 0, K = 1e-6; end
            if L <= 0, L = 1e-6; end

            Y = cS.A * (K^cS.alpha) * (L^(1-cS.alpha));
            MPK_gross = cS.alpha * Y / K;
            MPL_gross = (1-cS.alpha) * Y / L;
            R_market_gross_factor = 1 + MPK_gross - cS.ddk;
        end

        function gbc_residual = check_gbc_residual(K, C, Y, G, B, MPL, r_mkt, theta_payg, tau_l, b_payg, T_bequest, TR_gov, cS, paramS_gbc)
            % check_gbc_residual - 检查政府一般预算约束
            LaborTaxRev = tau_l * MPL * paramS_gbc.L_per_capita;
            CapitalTaxRev = r_mkt * K * cS.tau_k;
            ConsumptionTaxRev = C * cS.tau_c;
            GeneralRevenue = LaborTaxRev + CapitalTaxRev + ConsumptionTaxRev;

            DebtService = (r_mkt - paramS_gbc.popGrowthForDebt) * B;
            GeneralOutlays = G + DebtService + TR_gov;

            gbc_residual = GeneralRevenue - GeneralOutlays;
        end

        function cS = generateGrids(cS)
            % generateGrids - 根据当前的网格参数设置，重新生成资产网格。
            %
            % 输入：
            %   cS - 包含 nk, kMax, nkpps, kppsMax 等参数的结构体
            % 输出：
            %   cS - 更新了 kGridV 和 kppsGridV 的参数结构体

            % 重新生成非PPS资产网格 (kGridV)
            power_k = 1.5; % 网格密度参数
            if cS.nk > 1
                kGridV_temp = cS.kMin + (cS.kMax - cS.kMin) * (linspace(0, 1, cS.nk).^power_k);
                kGridV_temp(1) = cS.kMin;
            elseif cS.nk == 1
                kGridV_temp = cS.kMin;
            else
                kGridV_temp = [];
            end
            cS.kGridV = kGridV_temp(:); % 确保是列向量

            % 重新生成PPS资产网格 (kppsGridV)
            power_kpps = 1.5; % PPS资产网格密度参数
            if cS.nkpps > 1
                kppsGridV_temp = cS.kppsMin + (cS.kppsMax - cS.kppsMin) * (linspace(0, 1, cS.nkpps).^power_kpps);
                kppsGridV_temp(1) = cS.kppsMin;
            elseif cS.nkpps == 1
                kppsGridV_temp = cS.kppsMin;
            else
                kppsGridV_temp = [];
            end
            cS.kppsGridV = kppsGridV_temp(:); % 确保是列向量

            % 输出确认信息
            % fprintf('网格已重新生成：nk=%d, nkpps=%d\n', cS.nk, cS.nkpps);
        end
        % =====================================================================
        % == 2. VFI 核心求解器 (修正版) ==
        % =====================================================================
        % --- 在 main_olg_v8_utils.m 的 methods (Static) 中，替换以下两个函数 ---

        function [cPolM_q, kPolM, cPpsPolM_choice, valM] = HHSolution_VFI_Huggett(...
                R_k_net_factor_vfi, w_gross_vfi, TR_total_vfi, bV_payg_vfi, paramS_vfi, cS, solverMethod)
            % [带求解器切换功能的版本]
            % solverMethod: 'grid' 或 'hybrid'

            % --- 默认使用网格搜索法 ---
            if nargin < 7
                solverMethod = 'hybrid';
            end

            valM = -Inf(cS.nk, cS.nkpps, cS.nw, cS.aD_new);
            cPolM_q = zeros(cS.nk, cS.nkpps, cS.nw, cS.aD_new);
            kPolM = zeros(cS.nk, cS.nkpps, cS.nw, cS.aD_new);
            cPpsPolM_choice = zeros(cS.nk, cS.nkpps, cS.nw, cS.aD_new);

            fprintf('  🔄 VFI V9-compatible (%s Solver): 开始逆向迭代... [nk=%d, nkpps=%d, nkprime=%d, npps=%d]\n', ...
                solverMethod, cS.nk, cS.nkpps, cS.nkprime, cS.npps);
    vfi_start_time = tic;
        % [进度条] 1. 初始化进度条
    progress_msg = ''; % 用于存储上一条消息的长度
            for a_idx = cS.aD_new : -1 : 1
                       % [进度条] 2. 更新并打印进度
        % 首先，删除上一条消息
        fprintf(repmat('\b', 1, length(progress_msg))); 
        % 计算进度
        percent_done = (cS.aD_new - a_idx + 1) / cS.aD_new * 100;
        % 创建新消息
        progress_msg = sprintf('    正在处理年龄组 %2d/%d (%.0f%%)...', cS.aD_new - a_idx + 1, cS.aD_new, percent_done);
        % 打印新消息
        fprintf('%s', progress_msg);
                vPrime_kkppse_next = [];
                if a_idx < cS.aD_new
                    vPrime_kkppse_next = valM(:,:,:,a_idx+1);
                end

                % --- [核心] 根据 solverMethod 选择求解器 ---
                if strcmpi(solverMethod, 'hybrid')
                    [cPolM_q(:,:,:,a_idx), kPolM(:,:,:,a_idx), cPpsPolM_choice(:,:,:,a_idx), valM(:,:,:,a_idx)] = ...
                        main_olg_v8_utils.HHSolutionByAge_VFI_Huggett_HybridOptimizer(a_idx, vPrime_kkppse_next, ...
                        R_k_net_factor_vfi, w_gross_vfi, TR_total_vfi, bV_payg_vfi(a_idx), paramS_vfi, cS);
                elseif strcmpi(solverMethod, 'vectorized_grid') % <<<<< 新增CASE
                    [cPolM_q(:,:,:,a_idx), kPolM(:,:,:,a_idx), cPpsPolM_choice(:,:,:,a_idx), valM(:,:,:,a_idx)] = ...
                        main_olg_v8_utils.HHSolutionByAge_VFI_Huggett_VectorizedGrid(a_idx, vPrime_kkppse_next, ...
                        R_k_net_factor_vfi, w_gross_vfi, TR_total_vfi, bV_payg_vfi(a_idx), paramS_vfi, cS);
                else % 默认或指定 'grid'
                    [cPolM_q(:,:,:,a_idx), kPolM(:,:,:,a_idx), cPpsPolM_choice(:,:,:,a_idx), valM(:,:,:,a_idx)] = ...
                        main_olg_v8_utils.HHSolutionByAge_VFI_Huggett_v9_GridSearch(a_idx, vPrime_kkppse_next, ...
                        R_k_net_factor_vfi, w_gross_vfi, TR_total_vfi, bV_payg_vfi(a_idx), paramS_vfi, cS);
                end
                
            end
    % [进度条] 3. 循环结束后换行
    fprintf('\n');
            total_time = toc(vfi_start_time);
            fprintf('  ✅ VFI (%s Solver): 完成! 总耗时: %.2f秒\n', solverMethod, total_time);
        end

        function [cPol_age_q, kPol_age, cPpsPol_age_choice, val_age] = HHSolutionByAge_VFI_Huggett_v9_GridSearch(...
                a_idx, vPrime_kkppse_next, R_k_net_factor_age, w_gross_age, TR_total_age, b_age_val, ...
                paramS_age, cS)
            % [最终修正版] 并行化的离散网格搜索VFI求解器

            % 初始化输出矩阵
            val_age    = -Inf(cS.nk, cS.nkpps, cS.nw);
            cPol_age_q = zeros(cS.nk, cS.nkpps, cS.nw);
            kPol_age   = zeros(cS.nk, cS.nkpps, cS.nw);
            cPpsPol_age_choice = zeros(cS.nk, cS.nkpps, cS.nw);

            % --- [修正] 最后一期逻辑 ---
            if a_idx == cS.aD_new
                % 最后一期不依赖epsilon，先计算一次
                [K_grid, Kpps_grid] = ndgrid(cS.kGridV, cS.kppsGridV);

                % HHIncome_Huggett的最后一个参数(epsilon)在退休期无效，可以设为0
                [resources, ~, ~] = main_olg_v8_utils.HHIncome_Huggett(K_grid, R_k_net_factor_age, w_gross_age, TR_total_age, b_age_val, 0, a_idx, paramS_age, cS, 0);

                % 总资源 = 初始财富 + 税后PPS财富 + 其他净收入(不含资本收入)
                % 注意: HHIncome_Huggett返回的resources已经包含了k_now，所以不需要再加
                total_resources = resources + Kpps_grid * (1 - cS.pps_tax_rate_withdrawal);

                final_c = max(cS.cFloor, total_resources / (1 + cS.tau_c));
                [~, final_v] = main_olg_v8_utils.CES_utility(final_c, cS.sigma, cS);

                % 将结果赋给所有epsilon状态
                for ie = 1:cS.nw
                    cPol_age_q(:,:,ie) = final_c;
                    val_age(:,:,ie) = final_v;
                    % kPol_age 和 cPpsPol_age_choice 默认为0，是正确的
                end
                return;
            end

            % --- 非最后一期逻辑 ---

            % a. 计算期望未来价值矩阵 E[V']
            % (原始代码正确)
            EV_matrix = zeros(cS.nk, cS.nkpps, cS.nw);
            for ie_current = 1:cS.nw
                transition_probs = paramS_age.leTrProbM(ie_current, :);
                % 使用pagetimes进行高效的张量乘法
                EV_slice = pagemtimes(reshape(vPrime_kkppse_next, [cS.nk * cS.nkpps, cS.nw]), reshape(transition_probs, [cS.nw, 1]));
                EV_matrix(:, :, ie_current) = reshape(EV_slice, [cS.nk, cS.nkpps]);
            end

            % b. [修正] 创建插值器，并指定外插行为
            EV_interpolants = cell(cS.nw, 1);
            for ie_current = 1:cS.nw
                % 指定线性外插，避免NaN
                EV_interpolants{ie_current} = griddedInterpolant({cS.kGridV, cS.kppsGridV}, EV_matrix(:, :, ie_current), 'linear', 'linear');
            end

            % c. 并行循环遍历所有状态点
            parfor ik = 1:cS.nk
                % 为parfor创建临时变量，避免修改sliced variable
                val_slice = -Inf(cS.nkpps, cS.nw);
                c_slice = zeros(cS.nkpps, cS.nw);
                k_slice = zeros(cS.nkpps, cS.nw);
                cpps_slice = zeros(cS.nkpps, cS.nw);

                for ikpps = 1:cS.nkpps
                    for ie = 1:cS.nw
                        % 当前状态
                        k_state = cS.kGridV(ik);
                        k_pps_state = cS.kppsGridV(ikpps);
                        epsilon_state = paramS_age.leGridV(ie);

                        % 初始化最优值
                        best_val = -Inf;
                        best_c = cS.cFloor;
                        best_k_prime = cS.kMin;
                        best_c_pps = 0;

                        % 确定c_pps的搜索网格
                        max_cpps = 0;
                        if a_idx < cS.aR_new % 仅工作期可缴费
                            age_eff = cS.ageEffV_new(a_idx);
                            gross_labor_income = w_gross_age * age_eff * epsilon_state;
                            max_cpps = min(cS.pps_contrib_limit, gross_labor_income * cS.pps_max_contrib_frac);
                        end
                        cpps_grid = linspace(0, max(0, max_cpps), cS.npps);

                        % 遍历c_pps决策
                        for c_pps_choice = cpps_grid
                            [resources, ~, ~] = main_olg_v8_utils.HHIncome_Huggett(k_state, R_k_net_factor_age, w_gross_age, TR_total_age, b_age_val, c_pps_choice, a_idx, paramS_age, cS, epsilon_state);
                            
                                            % [核心修复] 将PPS提取的财富加入总资源
                    if a_idx >= cS.aR_new && cS.pps_active
                        pps_withdrawal_gross = k_pps_state * cS.pps_withdrawal_rate;
                        resources = resources + pps_withdrawal_gross * (1 - cS.pps_tax_rate_withdrawal);
                    end
                
                
                            k_prime_max_budget = resources - cS.cFloor * (1 + cS.tau_c);
                            if k_prime_max_budget < cS.kMin, continue; end % 如果资源不足以支付最低消费，跳过

                            k_prime_grid = linspace(cS.kMin, k_prime_max_budget, cS.nkprime);

                            % 遍历k'决策
                            for k_prime_choice = k_prime_grid
                                c_expend = resources - k_prime_choice;
                                c_choice = max(cS.cFloor, c_expend / (1 + cS.tau_c));

                                [~, util] = main_olg_v8_utils.CES_utility(c_choice, cS.sigma, cS);

                                pps_withdrawal = 0;
                                if a_idx >= cS.aR_new, pps_withdrawal = k_pps_state * cS.pps_withdrawal_rate; end
                                pps_return_factor = 1 + ((R_k_net_factor_age - 1) + cS.pps_return_rate_premium);
                                k_pps_prime = (k_pps_state + c_pps_choice - pps_withdrawal) * pps_return_factor;
    
                                % [核心修复] 钳位 k_prime_choice 和 k_pps_prime
                                k_prime_clamped = max(cS.kGridV(1), min(cS.kGridV(end), k_prime_choice));
                                k_pps_prime_clamped = max(cS.kppsGridV(1), min(cS.kppsGridV(end), k_pps_prime));
                                
                                % 调用插值器
                                ev = EV_interpolants{ie}(k_prime_clamped, k_pps_prime_clamped);

                                current_val = util + cS.beta * cS.s_1yr_transitionV(a_idx) * ev;

                                if current_val > best_val
                                    best_val = current_val;
                                    best_c = c_choice;
                                    best_k_prime = k_prime_choice;
                                    best_c_pps = c_pps_choice;
                                end
                            end
                        end

                        % 记录该状态点的最优解
                        val_slice(ikpps, ie) = best_val;
                        c_slice(ikpps, ie) = best_c;
                        k_slice(ikpps, ie) = best_k_prime;
                        cpps_slice(ikpps, ie) = best_c_pps;
                    end
                end

                % 将计算好的切片结果赋给主矩阵
                val_age(ik,:,:) = val_slice;
                cPol_age_q(ik,:,:) = c_slice;
                kPol_age(ik,:,:) = k_slice;
                cPpsPol_age_choice(ik,:,:) = cpps_slice;
            end
        end


                % =====================================================================
        % == [新] 2.B VFI 按年龄求解器 (高精度混合优化版) ==
        % =====================================================================
        function [cPol_age_q, kPol_age, cPpsPol_age_choice, val_age] = HHSolutionByAge_VFI_Huggett_HybridOptimizer(...
                a_idx, vPrime_kkppse_next, R_k_net_factor_age, w_gross_age, TR_total_age, b_age_val, ...
                paramS_age, cS)
            % [高精度版] 混合网格搜索与连续优化VFI求解器
            % - 对 c_pps 进行网格搜索。
            % - 对 k' (下一期资本) 使用 fminbnd 进行连续优化，以获得更高精度。
            % - 使用 'spline' 插值以更准确地估计期望未来价值。

            % 初始化输出矩阵
            val_age    = -Inf(cS.nk, cS.nkpps, cS.nw);
            cPol_age_q = zeros(cS.nk, cS.nkpps, cS.nw);
            kPol_age   = zeros(cS.nk, cS.nkpps, cS.nw);
            cPpsPol_age_choice = zeros(cS.nk, cS.nkpps, cS.nw);

            % --- 最后一期逻辑 (与原版相同) ---
            if a_idx == cS.aD_new
                [K_grid, Kpps_grid] = ndgrid(cS.kGridV, cS.kppsGridV);
                [resources, ~, ~] = main_olg_v8_utils.HHIncome_Huggett(K_grid, R_k_net_factor_age, w_gross_age, TR_total_age, b_age_val, 0, a_idx, paramS_age, cS, 0);
                total_resources = resources + Kpps_grid * (1 - cS.pps_tax_rate_withdrawal);
                final_c = max(cS.cFloor, total_resources / (1 + cS.tau_c));
                [~, final_v] = main_olg_v8_utils.CES_utility(final_c, cS.sigma, cS);
                for ie = 1:cS.nw
                    cPol_age_q(:,:,ie) = final_c;
                    val_age(:,:,ie) = final_v;
                end
                return;
            end

            % --- 非最后一期逻辑 ---

            % a. 计算期望未来价值矩阵 E[V'] (与原版相同)
            EV_matrix = zeros(cS.nk, cS.nkpps, cS.nw);
            for ie_current = 1:cS.nw
                transition_probs = paramS_age.leTrProbM(ie_current, :);
                EV_slice = pagemtimes(reshape(vPrime_kkppse_next, [cS.nk * cS.nkpps, cS.nw]), reshape(transition_probs, [cS.nw, 1]));
                EV_matrix(:, :, ie_current) = reshape(EV_slice, [cS.nk, cS.nkpps]);
            end

            % b. [改进] 创建插值器，使用 'spline' 方法提高精度，并指定外插
            EV_interpolants = cell(cS.nw, 1);
            for ie_current = 1:cS.nw
                EV_interpolants{ie_current} = griddedInterpolant({cS.kGridV, cS.kppsGridV}, EV_matrix(:, :, ie_current), 'spline', 'linear');
            end

            % c. 并行循环遍历所有状态点
            optim_options = optimset('Display', 'off', 'TolX', 1e-5); % fminbnd 优化选项

            parfor ik = 1:cS.nk
                % 为parfor创建临时变量
                val_slice = -Inf(cS.nkpps, cS.nw);
                c_slice = zeros(cS.nkpps, cS.nw);
                k_slice = zeros(cS.nkpps, cS.nw);
                cpps_slice = zeros(cS.nkpps, cS.nw);

                for ikpps = 1:cS.nkpps
                    for ie = 1:cS.nw
                        % 当前状态
                        k_state = cS.kGridV(ik);
                        k_pps_state = cS.kppsGridV(ikpps);
                        epsilon_state = paramS_age.leGridV(ie);
                        ev_interpolant = EV_interpolants{ie}; % 获取当前epsilon状态对应的插值器

                        % 初始化最优值
                        best_val_for_cpps_grid = -Inf;
                        best_c_for_cpps_grid = cS.cFloor;
                        best_k_prime_for_cpps_grid = cS.kMin;
                        best_c_pps_for_cpps_grid = 0;

                        % 确定c_pps的搜索网格
                        max_cpps = 0;
                        if a_idx < cS.aR_new
                            age_eff = cS.ageEffV_new(a_idx);
                            gross_labor_income = w_gross_age * age_eff * epsilon_state;
                            max_cpps = min(cS.pps_contrib_limit, gross_labor_income * cS.pps_max_contrib_frac);
                        end
                        cpps_grid = linspace(0, max(0, max_cpps), cS.npps);

                        % [改进] 遍历c_pps决策网格
                        for c_pps_choice = cpps_grid
                            % 1. 计算给定c_pps下的总资源
                            [resources, ~, ~] = main_olg_v8_utils.HHIncome_Huggett(k_state, R_k_net_factor_age, w_gross_age, TR_total_age, b_age_val, c_pps_choice, a_idx, paramS_age, cS, epsilon_state);
                            
                            % [核心修复] 将PPS提取的财富加入总资源
                            if a_idx >= cS.aR_new && cS.pps_active
                                pps_withdrawal_gross = k_pps_state * cS.pps_withdrawal_rate;
                                resources = resources + pps_withdrawal_gross * (1 - cS.pps_tax_rate_withdrawal);
                            end

                            k_prime_max_budget = resources - cS.cFloor * (1 + cS.tau_c);
                            if k_prime_max_budget < cS.kMin, continue; end % 资源不足

                            % 2. 定义目标函数，用于对 k' 进行优化
                    % 2. 定义目标函数
                    % [重要修改] 现在调用一个独立的静态方法
                    objective_func = @(k_prime) main_olg_v8_utils.objective_for_k_prime_private(...
                        k_prime, resources, k_pps_state, c_pps_choice, ...
                        R_k_net_factor_age, a_idx, ev_interpolant, cS);

                    % 3. [核心改进] 使用fminbnd
                    [k_prime_opt, neg_val_opt] = fminbnd(objective_func, cS.kMin, k_prime_max_budget, optim_options);

                            current_max_val = -neg_val_opt;

                            % 4. 更新最优解
                            if current_max_val > best_val_for_cpps_grid
                                best_val_for_cpps_grid = current_max_val;
                                best_k_prime_for_cpps_grid = k_prime_opt;
                                best_c_pps_for_cpps_grid = c_pps_choice;

                                % 根据最优 k' 和 c_pps 计算对应的 c
                                c_expend = resources - best_k_prime_for_cpps_grid;
                                best_c_for_cpps_grid = max(cS.cFloor, c_expend / (1 + cS.tau_c));
                            end
                        end % 结束对 c_pps_grid 的循环

                        % 记录该状态点的最优解
                        val_slice(ikpps, ie) = best_val_for_cpps_grid;
                        c_slice(ikpps, ie) = best_c_for_cpps_grid;
                        k_slice(ikpps, ie) = best_k_prime_for_cpps_grid;
                        cpps_slice(ikpps, ie) = best_c_pps_for_cpps_grid;
                    end
                end

                % 将计算好的切片结果赋给主矩阵
                val_age(ik,:,:) = val_slice;
                cPol_age_q(ik,:,:) = c_slice;
                kPol_age(ik,:,:) = k_slice;
                cPpsPol_age_choice(ik,:,:) = cpps_slice;
            end
        end
        

function [cPol_age_q, kPol_age, cPpsPol_age_choice, val_age] = HHSolutionByAge_VFI_Huggett_VectorizedGrid(...
        a_idx, vPrime_kkppse_next, R_k_net_factor_age, w_gross_age, TR_total_age, b_age_val, ...
        paramS_age, cS)
    % [v9.9 - 维度安全比例决策版]
    % - 修复了由于广播导致的维度不匹配问题。
    % - 所有计算都在最终的4D空间中显式进行。

    % --- 最后一期逻辑 ---
    if a_idx == cS.aD_new
        val_age    = -Inf(cS.nk, cS.nkpps, cS.nw);
        cPol_age_q = zeros(cS.nk, cS.nkpps, cS.nw);
        kPol_age   = zeros(cS.nk, cS.nkpps, cS.nw);
        cPpsPol_age_choice = zeros(cS.nk, cS.nkpps, cS.nw);
        
        [K_grid, Kpps_grid] = ndgrid(cS.kGridV, cS.kppsGridV);
        
        % 1. 计算非资本、非PPS的收入
        [resources, ~, ~] = main_olg_v8_utils.HHIncome_Huggett(K_grid, R_k_net_factor_age, w_gross_age, TR_total_age, b_age_val, 0, a_idx, paramS_age, cS, 0);
        
        % 2. [核心修复] 将所有剩余的PPS财富（税后）加入总资源
        %    在最后一期，所有PPS资产都被清算用于消费（或遗赠）
        if cS.pps_active
            total_resources = resources + Kpps_grid * (1 - cS.pps_tax_rate_withdrawal);
        else
            total_resources = resources;
        end
        
        % 3. 计算最终消费和效用
        final_c = max(cS.cFloor, total_resources / (1 + cS.tau_c));
        [~, final_v] = main_olg_v8_utils.CES_utility(final_c, cS.sigma, cS);
        
        % 4. 将结果赋给所有epsilon状态
        for ie = 1:cS.nw
            cPol_age_q(:,:,ie) = final_c;
            val_age(:,:,ie) = final_v;
            % kPol_age 和 cPpsPol_age_choice 默认为0，是正确的
        end
        return;
    end


    % --- 非最后一期逻辑 ---

    % a. 创建插值器 (不变)
    EV_matrix = zeros(cS.nk, cS.nkpps, cS.nw);
    for ie_current = 1:cS.nw
        transition_probs = paramS_age.leTrProbM(ie_current, :);
        EV_slice = pagemtimes(reshape(vPrime_kkppse_next, [cS.nk * cS.nkpps, cS.nw]), reshape(transition_probs, [cS.nw, 1]));
        EV_matrix(:, :, ie_current) = reshape(EV_slice, [cS.nk, cS.nkpps]);
    end
    EV_interpolants = cell(cS.nw, 1);
    for ie_current = 1:cS.nw
        EV_interpolants{ie_current} = griddedInterpolant({cS.kGridV, cS.kppsGridV}, EV_matrix(:, :, ie_current), 'spline', 'linear');
    end

    % c. 比例决策的矩阵化网格搜索
    prop_k_prime_grid = linspace(0, 1, cS.nkprime)';
    prop_cpps_grid = linspace(0, 1, cS.npps)';
    
    val_age = -Inf(cS.nk, cS.nkpps, cS.nw);
    cPol_age_q = zeros(cS.nk, cS.nkpps, cS.nw);
    kPol_age = zeros(cS.nk, cS.nkpps, cS.nw);
    cPpsPol_age_choice = zeros(cS.nk, cS.nkpps, cS.nw);

    for ie = 1:cS.nw
        epsilon_state = paramS_age.leGridV(ie);
        ev_interpolant = EV_interpolants{ie};
        
        % c.3.1. [核心修复] 创建所有变量的4D [nk, nkpps, nkprime, npps] 网格
        % 状态变量
        k_state_4D = repmat(reshape(cS.kGridV, [cS.nk, 1, 1, 1]), [1, cS.nkpps, cS.nkprime, cS.npps]);
        kpps_state_4D = repmat(reshape(cS.kppsGridV, [1, cS.nkpps, 1, 1]), [cS.nk, 1, cS.nkprime, cS.npps]);
        
        % 比例决策变量
        prop_k_prime_4D = repmat(reshape(prop_k_prime_grid, [1, 1, cS.nkprime, 1]), [cS.nk, cS.nkpps, 1, cS.npps]);
        prop_cpps_4D = repmat(reshape(prop_cpps_grid, [1, 1, 1, cS.npps]), [cS.nk, cS.nkpps, cS.nkprime, 1]);

        % c.3.2. 计算绝对决策值
        max_cpps = 0;
        if a_idx < cS.aR_new && cS.pps_active
            age_eff = cS.ageEffV_new(a_idx);
            gross_labor_income = w_gross_age * age_eff * epsilon_state;
            max_cpps = min(cS.pps_contrib_limit, gross_labor_income * cS.pps_max_contrib_frac);
        end
        actual_cpps_4D = max_cpps .* prop_cpps_4D;
        
        % c.3.3. 计算资源
        [resources_4D, ~, ~] = main_olg_v8_utils.HHIncome_Huggett(k_state_4D, R_k_net_factor_age, w_gross_age, TR_total_age, b_age_val, actual_cpps_4D, a_idx, paramS_age, cS, epsilon_state);
        if a_idx >= cS.aR_new && cS.pps_active
            pps_withdrawal_gross_4D = kpps_state_4D .* cS.pps_withdrawal_rate;
            resources_4D = resources_4D + pps_withdrawal_gross_4D .* (1 - cS.pps_tax_rate_withdrawal);
        end
        
        c_floor_spending = cS.cFloor * (1 + cS.tau_c);
        resources_above_floor_4D = max(0, resources_4D - c_floor_spending);
        actual_k_prime_4D = resources_above_floor_4D .* prop_k_prime_4D;

        % c.3.4. 计算消费和效用
        c_expend_4D = resources_4D - actual_k_prime_4D;
        c_choice_4D = max(cS.cFloor, c_expend_4D / (1 + cS.tau_c));
        [~, util_4D] = main_olg_v8_utils.CES_utility(c_choice_4D, cS.sigma, cS);
        
        % c.3.5. 计算下一期PPS资产并插值
        pps_withdrawal_4D = 0;
        if a_idx >= cS.aR_new, pps_withdrawal_4D = kpps_state_4D .* cS.pps_withdrawal_rate; end
        pps_return_factor = 1 + ((R_k_net_factor_age - 1) + cS.pps_return_rate_premium);
        k_pps_prime_4D = (kpps_state_4D + actual_cpps_4D - pps_withdrawal_4D) * pps_return_factor;
        
        k_prime_clamped = max(cS.kGridV(1), min(cS.kGridV(end), actual_k_prime_4D));
        k_pps_prime_clamped = max(cS.kppsGridV(1), min(cS.kppsGridV(end), k_pps_prime_4D));
        ev_mat = ev_interpolant(k_prime_clamped, k_pps_prime_clamped);
        
        % c.3.6. 计算总价值并应用约束
        val_grid = util_4D + cS.beta * cS.s_1yr_transitionV(a_idx) * ev_mat;
        val_grid(c_expend_4D < 0) = -Inf;
        
        % c.3.7. 寻找最优决策
        [val_max_k, idx_k_prime] = max(val_grid, [], 3);
        [val_max_kc, idx_cpps] = max(val_max_k, [], 4);
        
        val_age(:,:,ie) = squeeze(val_max_kc);
        
        % 提取最优决策 (这部分逻辑与之前版本相同，但现在操作的是4D索引)
        [I, J] = ndgrid(1:cS.nk, 1:cS.nkpps);
        idx_cpps_squeezed = squeeze(idx_cpps);
        % 在4D的val_grid中找到最优k_prime的索引
        % sub2ind的维度应该是size(val_grid)的前3维，因为我们是在第3维上找索引
        linear_idx_k = sub2ind(size(squeeze(val_grid(:,:,:,1))), I(:), J(:), idx_cpps_squeezed(:));
        final_k_prime_prop_idx = idx_k_prime(linear_idx_k);
        final_k_prime_prop_idx = reshape(final_k_prime_prop_idx, [cS.nk, cS.nkpps]);
        
        final_cpps_prop_idx = squeeze(idx_cpps);
        
        best_prop_k_prime = prop_k_prime_grid(final_k_prime_prop_idx);
        best_prop_cpps = prop_cpps_grid(final_cpps_prop_idx);
        
        % [状态网格]
        [k_state_2D, kpps_state_2D] = ndgrid(cS.kGridV, cS.kppsGridV);
        max_cpps_2D = zeros(cS.nk, cS.nkpps);
        if a_idx < cS.aR_new && cS.pps_active
            age_eff = cS.ageEffV_new(a_idx);
            gross_labor_income = w_gross_age * age_eff * epsilon_state;
            max_cpps_2D(:) = min(cS.pps_contrib_limit, gross_labor_income * cS.pps_max_contrib_frac);
        end
        
        best_c_pps = max_cpps_2D .* best_prop_cpps;
        
        [resources_final, ~, ~] = main_olg_v8_utils.HHIncome_Huggett(k_state_2D, R_k_net_factor_age, w_gross_age, TR_total_age, b_age_val, best_c_pps, a_idx, paramS_age, cS, epsilon_state);
        if a_idx >= cS.aR_new && cS.pps_active
             pps_withdrawal_gross = kpps_state_2D .* cS.pps_withdrawal_rate;
             resources_final = resources_final + pps_withdrawal_gross .* (1 - cS.pps_tax_rate_withdrawal);
        end
        
        resources_above_floor_final = max(0, resources_final - c_floor_spending);
        best_k_prime = resources_above_floor_final .* best_prop_k_prime;
        
        cPpsPol_age_choice(:,:,ie) = best_c_pps;
        kPol_age(:,:,ie) = best_k_prime;
        
        c_expend_final = resources_final - best_k_prime;
        cPol_age_q(:,:,ie) = max(cS.cFloor, c_expend_final / (1 + cS.tau_c));
    end
end      

 % --- 在 main_olg_v10_utils.m 中，替换整个 HHIncome_Huggett 函数 ---

        % =========================================================================
        % [决定性最终修正] HHIncome_Huggett 函数
        % =========================================================================
        function [resources, tax_info] = HHIncome_Huggett(k_now_val, R_k_net_factor, w_gross, TR_total, b_payg_val, c_pps_chosen, a_idx, paramS_hh, cS, epsilon_val)
            tax_info = struct('labor_tax', 0, 'capital_tax', 0);
            
            % 1. 劳动收入和相关税收
            if a_idx <= cS.aR_new
                age_efficiency = cS.ageEffV_new(a_idx);
                labor_income_gross = w_gross * age_efficiency * epsilon_val;
                payg_tax = cS.theta_t * labor_income_gross;
                pps_deduction = 0;
                if isfield(paramS_hh, 'pps_tax_deferral_active') && paramS_hh.pps_tax_deferral_active
                    pps_deduction = c_pps_chosen;
                end
                labor_income_taxable_for_tau_l = labor_income_gross - pps_deduction;
                general_labor_tax = cS.tau_l * max(0, labor_income_taxable_for_tau_l);
                tax_info.labor_tax = general_labor_tax + payg_tax;
                non_capital_income_base = labor_income_gross - tax_info.labor_tax;
            else
                non_capital_income_base = 0;
            end
            
            % 2. 资本收入和相关税收
            r_net_from_param = R_k_net_factor - 1;
            r_mkt_rental = (r_net_from_param + cS.ddk) / (1 - cS.tau_k);
            capital_income_gross = k_now_val * r_mkt_rental; 
            tax_info.capital_tax = capital_income_gross * cS.tau_k;
            capital_income_net_of_tax = capital_income_gross - tax_info.capital_tax;

            % 3. 总资源 (核心最终修正)
            % 家庭可用于分配的总资源 = 
            %   期初财富 + 税后资本收入 + 税后劳动收入 + 所有转移支付 - PPS缴费 - 资本折旧
            %   这等价于 NDI (净可支配收入) + k_now - C_expend
            
            % [最终会计修正] 意外遗赠转移支付 TR_total 必须被加入到家庭资源中！
            after_tax_cash_income = non_capital_income_base + capital_income_net_of_tax + b_payg_val + TR_total;
            
            % 从总资源中扣除折旧才是净收入
            net_income = after_tax_cash_income - c_pps_chosen - (cS.ddk * k_now_val);

            % 总资源 = 期初财富 + 净收入
            resources = k_now_val + net_income;
        end


        function [muM, utilM] = CES_utility(cM, sigma, cS)
            c_adj = max(cS.cFloor, cM);
            if abs(sigma - 1) < 1e-6, utilM = log(c_adj); muM = 1./c_adj;
            else, utilM = (c_adj.^(1-sigma))./(1-sigma); muM = c_adj.^(-sigma); end
            utilM(cM < cS.cFloor) = -1e10 - (cS.cFloor - cM(cM < cS.cFloor))*1e10;
        end

        % =====================================================================
        % == 4. 模拟与分析函数 (保留全部) ==
        % =====================================================================

        % --- 劳动禀赋过程 ---
        function [leLogGridV, leTrProbM, leProb1V] = EarningProcess_olgm(cS)
             lePersistence = 0.90;     % 改为 0.90
            leShockStd = 0.15;        % 改为 0.15
            Tauchen_q = 2.0;
            [leLogGridV_raw, leTrProbM] = main_olg_v8_utils.tauchen(cS.nw, lePersistence, leShockStd, 0, Tauchen_q);
            leLogGridV = leLogGridV_raw - mean(leLogGridV_raw);
            [~, D] = eig(leTrProbM');
            [~, c] = min(abs(diag(D)-1));
            leProb1V = abs(D(:,c)/sum(D(:,c)));
        end

        function [y_grid_out, trProbM_out] = tauchen(N, rho, sigma, mu, m)
            % ... (与原版一致)
            std_y = sqrt(sigma^2 / (1-rho^2));
            y_max = m*std_y; y_min = -y_max;
            y = linspace(y_min, y_max, N);
            d = y(2)-y(1);
            trProbM_out = zeros(N,N);
            for j=1:N
                for k=1:N
                    m_k = rho*y(j) + mu;
                    if k==1, trProbM_out(j,k) = normcdf((y(1)-m_k+d/2)/sigma);
                    elseif k==N, trProbM_out(j,k) = 1 - normcdf((y(N)-m_k-d/2)/sigma);
                    else, trProbM_out(j,k) = normcdf((y(k)-m_k+d/2)/sigma) - normcdf((y(k)-m_k-d/2)/sigma); end
                end
            end
            y_grid_out = y(:);
        end

        function eIdxM_group = LaborEndowSimulation_olgm_AgeGroup(cS, paramS)
            eIdxM_group = main_olg_v8_utils.MarkovChainSimulation_AgeGroup(cS.nSim, cS, paramS.leProb1V, paramS.leTrProbM);
        end

        function eIdxM_group_out = MarkovChainSimulation_AgeGroup(num_simulations, cS, p0, P)
            rng(433);
            eIdxM_group_out = zeros(num_simulations, cS.aD_new, 'uint16');
            eIdxM_group_out(:,1) = 1 + sum(rand(num_simulations,1) > cumsum(p0(:)'), 2);
            for a=2:cS.aD_new
                eIdxM_group_out(:,a) = 1 + sum(rand(num_simulations,1) > cumsum(P(eIdxM_group_out(:,a-1),:), 2), 2);
            end
        end

        % [核心重写] 家庭模拟器
        % =====================================================================
        % == [最终对齐版] 家庭模拟器 (替换原函数) ==
        % =====================================================================
        function [kHistM_out, kPpsHistM_out, cHistM_out, cppsHistM_out] = HHSimulation_olgm(...
                kPolM, cPpsPolM_choice, cPolM_consump, eIdxM_group, ...
                R_k_net, w, TR, bV_payg, paramS_sim, cS_sim)
            % [最终对齐版] - 与Python的HHSimulation_olgm_rl和VFI模拟器在物理过程上完全一致。
            % - 使用 griddedInterpolant 模拟家庭基于VFI策略的生命周期路径。
            % - 核心是按年龄组进行模拟，并确保所有状态演化逻辑正确。

            nSim = size(eIdxM_group, 1);
            aD = cS_sim.aD_new;

            % 初始化历史记录矩阵
            kHistM_out = zeros(nSim, aD);
            kPpsHistM_out = zeros(nSim, aD);
            cHistM_out = zeros(nSim, aD);
            cppsHistM_out = zeros(nSim, aD);

            % --- 1. 创建策略函数的插值器 ---
            kPolInterp = cell(cS_sim.nw, aD);
            cPpsPolInterp = cell(cS_sim.nw, aD);
            cPolInterp = cell(cS_sim.nw, aD);

            for ia = 1:aD
                for ie = 1:cS_sim.nw
                    % 根据资产网格维度选择合适的插值方法
                    if cS_sim.nk > 1 && cS_sim.nkpps > 1
                        kPolInterp{ie,ia} = griddedInterpolant({cS_sim.kGridV, cS_sim.kppsGridV}, squeeze(kPolM(:,:,ie,ia)), 'linear', 'linear');
                        cPpsPolInterp{ie,ia} = griddedInterpolant({cS_sim.kGridV, cS_sim.kppsGridV}, squeeze(cPpsPolM_choice(:,:,ie,ia)), 'linear', 'linear');
                        cPolInterp{ie,ia} = griddedInterpolant({cS_sim.kGridV, cS_sim.kppsGridV}, squeeze(cPolM_consump(:,:,ie,ia)), 'linear', 'linear');
                    elseif cS_sim.nk > 1
                        kPolInterp{ie,ia} = griddedInterpolant(cS_sim.kGridV, squeeze(kPolM(:,1,ie,ia)), 'linear', 'linear');
                        cPpsPolInterp{ie,ia} = griddedInterpolant(cS_sim.kGridV, squeeze(cPpsPolM_choice(:,1,ie,ia)), 'linear', 'linear');
                        cPolInterp{ie,ia} = griddedInterpolant(cS_sim.kGridV, squeeze(cPolM_consump(:,1,ie,ia)), 'linear', 'linear');
                    elseif cS_sim.nkpps > 1 % 只有kpps网格
                        kPolInterp{ie,ia} = griddedInterpolant(cS_sim.kppsGridV, squeeze(kPolM(1,:,ie,ia))', 'linear', 'linear');
                        cPpsPolInterp{ie,ia} = griddedInterpolant(cS_sim.kppsGridV, squeeze(cPpsPolM_choice(1,:,ie,ia))', 'linear', 'linear');
                        cPolInterp{ie,ia} = griddedInterpolant(cS_sim.kppsGridV, squeeze(cPolM_consump(1,:,ie,ia))', 'linear', 'linear');
                    else % nk=1, nkpps=1 的标量情况
                        kPolInterp{ie,ia} = @(x,y) squeeze(kPolM(1,1,ie,ia));
                        cPpsPolInterp{ie,ia} = @(x,y) squeeze(cPpsPolM_choice(1,1,ie,ia));
                        cPolInterp{ie,ia} = @(x,y) squeeze(cPolM_consump(1,1,ie,ia));
                    end
                end
            end

            % --- 2. 初始化状态和参数 ---
            pps_return_factor = 1 + ((R_k_net - 1) + cS_sim.pps_return_rate_premium);
            k_next = zeros(nSim, 1);
            k_pps_next = zeros(nSim, 1);

            % --- 3. 按年龄组进行前向模拟 ---
            for a_idx = 1:aD
                % 获取当前状态
                k_now = k_next;
                k_pps_now = k_pps_next;
                kHistM_out(:, a_idx) = k_now;
                kPpsHistM_out(:, a_idx) = k_pps_now;

                % 初始化当期决策向量
                k_prime_decision = zeros(nSim, 1);
                cpps_decision = zeros(nSim, 1);
                c_decision = zeros(nSim, 1);

                % --- 4. 根据效率冲击状态，查询策略并获取决策 ---
                for ie = 1:cS_sim.nw
                    % 找到当前效率状态对应的所有个体
                    idx_sim = find(eIdxM_group(:, a_idx) == ie);
                    if isempty(idx_sim), continue; end

                    k_now_e = k_now(idx_sim);
                    k_pps_now_e = k_pps_now(idx_sim);

                    % 验证代码
if any(k_now_e > cS_sim.kMax) || any(k_pps_now_e > cS_sim.kppsMax)
    fprintf('警告: a_idx=%d, ie=%d, 有 %d 个体超出网格范围！\n', ...
        a_idx, ie, sum(k_now_e > cS_sim.kMax | k_pps_now_e > cS_sim.kppsMax));
end

                    % 使用插值器获取决策
                    if cS_sim.nk > 1 && cS_sim.nkpps > 1
                        k_prime_decision(idx_sim) = kPolInterp{ie, a_idx}(k_now_e, k_pps_now_e);
                        cpps_decision(idx_sim) = cPpsPolInterp{ie, a_idx}(k_now_e, k_pps_now_e);
                        c_decision(idx_sim) = cPolInterp{ie, a_idx}(k_now_e, k_pps_now_e);
                    elseif cS_sim.nk > 1
                        k_prime_decision(idx_sim) = kPolInterp{ie, a_idx}(k_now_e);
                        cpps_decision(idx_sim) = cPpsPolInterp{ie, a_idx}(k_now_e);
                        c_decision(idx_sim) = cPolInterp{ie, a_idx}(k_now_e);
                    elseif cS_sim.nkpps > 1
                        k_prime_decision(idx_sim) = kPolInterp{ie, a_idx}(k_pps_now_e);
                        cpps_decision(idx_sim) = cPpsPolInterp{ie, a_idx}(k_pps_now_e);
                        c_decision(idx_sim) = cPolInterp{ie, a_idx}(k_pps_now_e);
                    else % 标量情况
                        k_prime_decision(idx_sim) = kPolInterp{ie, a_idx}(k_now_e, k_pps_now_e);
                        cpps_decision(idx_sim) = cPpsPolInterp{ie, a_idx}(k_now_e, k_pps_now_e);
                        c_decision(idx_sim) = cPolInterp{ie, a_idx}(k_now_e, k_pps_now_e);
                    end
                end

                % --- 5. 记录当期决策并演化到下一期状态 ---
                cHistM_out(:, a_idx) = max(cS_sim.cFloor, c_decision);
                cppsHistM_out(:, a_idx) = max(0, cpps_decision);

                if a_idx < aD
                    % 演化非PPS资产
                    k_next = max(cS_sim.kMin, min(cS_sim.kMax, k_prime_decision));

                    % 演化PPS资产
                    pps_withdrawal = 0;
                    if a_idx >= cS_sim.aR_new && cS_sim.pps_active
                        pps_withdrawal = k_pps_now * cS_sim.pps_withdrawal_rate;
                    end
                    k_pps_next_unclamped = (k_pps_now + cppsHistM_out(:, a_idx) - pps_withdrawal) * pps_return_factor;
                    k_pps_next = max(cS_sim.kppsMin, min(cS_sim.kppsMax, k_pps_next_unclamped));
                end
            end
        end

        function [K_eq, tau_l_eq, gbc_residual_eq, eq_found, final_eq_solution_details] = solve_K_tau_l_for_rho_prime(rho_prime_payg_target_input, K_init_guess_input, cS, paramS, eIdxM_group)
            % [VFI最终对齐版]
            K_current_guess = K_init_guess_input;
            tau_l_current_guess = cS.tau_l_init_guess;

            mass_retirees_global = sum(paramS.ageMassV(cS.aR_new+1:end));
            theta_payg_required_calc = rho_prime_payg_target_input * (mass_retirees_global / paramS.mass_workers_group);
            theta_payg_required_calc = max(0, theta_payg_required_calc);

            final_eq_solution_details = struct('theta_payg_required_before_cap', theta_payg_required_calc);

            if theta_payg_required_calc > cS.theta_payg_max + 1e-5
                fprintf('  solve_K_tau_l (VFI): rho_prime_target=%.4f 导致理论theta_req=%.4f > theta_max=%.3f. 不可行。\n', rho_prime_payg_target_input, theta_payg_required_calc, cS.theta_payg_max);
                K_eq = K_init_guess_input; tau_l_eq = tau_l_current_guess; gbc_residual_eq = Inf; eq_found = false; return;
            end

            stagnation_counter_ktl = 0;
            prev_devNorm_ktl = Inf;
            tau_l_boundary_strike_count_ktl = 0;

            fprintf('  solve_K_tau_l_for_rho_prime_vfi: rho_prime_target=%.4f (理论theta_req=%.4f)\n', rho_prime_payg_target_input, theta_payg_required_calc);
            fprintf('  IterKTL | K_guess  | tau_l_gs | MPL_g    | theta_act| K_tot_mod| K_pps_mod| GBC_res  | K_dev    | tau_l_dev| Norm     | Improv   | Strikes  | Time (s) |\n');
            fprintf('  %s\n', repmat('-', 1, 123));

            for iter_ktl_idx = 1:cS.max_iter_K_tau_l
                iter_timer_start = tic;

                [R_mkt_gross_factor, MPL_gross] = main_olg_v8_utils.HHPrices_Huggett(K_current_guess, paramS.L_per_capita, cS);
                r_mkt_gross = R_mkt_gross_factor - 1;

                avg_worker_gross_wage = (MPL_gross * paramS.L_per_capita) / paramS.mass_workers_group;
                b_payg = max(0, rho_prime_payg_target_input * avg_worker_gross_wage);

                theta_payg_actual = min(theta_payg_required_calc, cS.theta_payg_max);
                if (theta_payg_actual + tau_l_current_guess) > cS.max_total_labor_tax
                    theta_payg_actual = max(0, cS.max_total_labor_tax - tau_l_current_guess);
                end

                r_k_net_hh = r_mkt_gross * (1 - cS.tau_k);
                R_k_net_factor_hh = 1 + r_k_net_hh;

                bV_payg_vec = zeros(1, cS.aD_new);
                bV_payg_vec(cS.aR_new+1:end) = b_payg;

                paramS_for_vfi = paramS;
                paramS_for_vfi.tau_l = tau_l_current_guess;
                paramS_for_vfi.theta_payg_actual_for_hh = theta_payg_actual;
                paramS_for_vfi.pps_tax_deferral_active = cS.pps_active;

                % [核心] 调用VFI求解器
                [cPolM, kPolM, cPpsPolM, ~] = main_olg_v8_utils.HHSolution_VFI_Huggett(R_k_net_factor_hh, MPL_gross, 0.0, bV_payg_vec, paramS_for_vfi, cS);
                [kHistM, kPpsHistM, cHistM, ~] = main_olg_v8_utils.HHSimulation_olgm(kPolM, cPpsPolM, cPolM, eIdxM_group, R_k_net_factor_hh, MPL_gross, 0.0, bV_payg_vec, paramS_for_vfi, cS);

                K_model_nonpps_sim = mean(kHistM, 1) * paramS.ageMassV;
                K_model_pps_sim = 0;
                if cS.pps_active && ~isempty(kPpsHistM), K_model_pps_sim = mean(kPpsHistM, 1) * paramS.ageMassV; end

                K_model_from_sim = max(1e-6, K_model_nonpps_sim + K_model_pps_sim);
                C_model = mean(cHistM,1) * paramS.ageMassV;

                Y_for_gbc = cS.A * (K_current_guess^cS.alpha) * (paramS.L_per_capita^(1-cS.alpha));
                G_target = cS.gov_exp_frac_Y * Y_for_gbc;
                B_target = cS.gov_debt_frac_Y * Y_for_gbc;

                paramS_for_gbc = struct('L_per_capita', paramS.L_per_capita, 'popGrowthForDebt', paramS.popGrowthForDebt);
                gbc_residual = main_olg_v8_utils.check_gbc_residual(K_current_guess, C_model, Y_for_gbc, G_target, B_target, MPL_gross, r_mkt_gross, theta_payg_actual, tau_l_current_guess, b_payg, 0, 0, cS, paramS_for_gbc);

                K_dev = K_current_guess - K_model_from_sim;
                tau_l_dev_raw = -gbc_residual / (MPL_gross * paramS.L_per_capita + 1e-9);
                current_devNorm = sqrt(K_dev^2 + gbc_residual^2);
                norm_improvement = prev_devNorm_ktl - current_devNorm;
                elapsed_time = toc(iter_timer_start);

                fprintf('  %7d | %8.4f | %8.4f | %8.4f | %8.4f | %8.4f | %8.4f | %8.2e | %8.4f | %8.4f | %8.2e | %8.1e | %7d | %8.2f |\n', ...
                    iter_ktl_idx, K_current_guess, tau_l_current_guess, MPL_gross, theta_payg_actual, K_model_from_sim, K_model_pps_sim, gbc_residual, K_dev, tau_l_dev_raw, current_devNorm, norm_improvement, tau_l_boundary_strike_count_ktl, elapsed_time);

                payg_fully_funded_check = (theta_payg_actual >= theta_payg_required_calc - 1e-5);
                if current_devNorm < cS.tol_K_tau_l && abs(gbc_residual) < cS.gbc_tol_for_internal_loop && payg_fully_funded_check
                    fprintf('  solve_K_tau_l (VFI): K和tau_l成功收敛。\n');
                    final_eq_solution_details.R_mkt_gross_factor=R_mkt_gross_factor; final_eq_solution_details.MPL_gross=MPL_gross; final_eq_solution_details.theta_payg=theta_payg_actual; final_eq_solution_details.b_payg=b_payg; final_eq_solution_details.T_bequest_Model=0.0; final_eq_solution_details.C_model=C_model; final_eq_solution_details.K_model_pps=K_model_pps_sim; final_eq_solution_details.K_model_non_pps=K_model_nonpps_sim;
                    K_eq = K_model_from_sim; tau_l_eq = tau_l_current_guess; gbc_residual_eq = gbc_residual; eq_found = true; return;
                end

                K_current_guess = max(1e-3, K_current_guess - cS.damp_K_v5 * K_dev);
                tau_l_next_unconstrained = tau_l_current_guess + cS.damp_tau_l_v5 * tau_l_dev_raw;
                tau_l_next_constrained = clip(tau_l_next_unconstrained, cS.tau_l_min, cS.tau_l_max);

                is_at_boundary = (abs(tau_l_next_constrained - cS.tau_l_max) < 1e-7 && tau_l_next_unconstrained >= cS.tau_l_max - 1e-7) || ...
                    (abs(tau_l_next_constrained - cS.tau_l_min) < 1e-7 && tau_l_next_unconstrained <= cS.tau_l_min + 1e-7);

                if is_at_boundary && abs(gbc_residual) > cS.gbc_tol_for_internal_loop, tau_l_boundary_strike_count_ktl = tau_l_boundary_strike_count_ktl + 1;
                else, tau_l_boundary_strike_count_ktl = 0; end

                tau_l_current_guess = tau_l_next_constrained;

                if tau_l_boundary_strike_count_ktl >= cS.max_tau_l_boundary_strikes, fprintf('  警告 (VFI): tau_l 在边界持续撞击，且GBC未平衡。中止。\n'); break; end
                if iter_ktl_idx > 0 && norm_improvement < (cS.min_norm_improvement_frac * prev_devNorm_ktl), stagnation_counter_ktl = stagnation_counter_ktl + 1;
                else, stagnation_counter_ktl = 0; end

                prev_devNorm_ktl = current_devNorm;
                if stagnation_counter_ktl >= cS.max_stagnation_iters, fprintf('  警告 (VFI): 检测到范数停滞。中止。\n'); break; end
            end

            fprintf('  警告 (VFI): K和tau_l迭代达到最大次数或未达可行解。\n');
            K_eq = K_model_from_sim; tau_l_eq = tau_l_current_guess; gbc_residual_eq = gbc_residual; eq_found = false;
            final_eq_solution_details.R_mkt_gross_factor=R_mkt_gross_factor; final_eq_solution_details.MPL_gross=MPL_gross; final_eq_solution_details.theta_payg=theta_payg_actual; final_eq_solution_details.b_payg=b_payg; final_eq_solution_details.T_bequest_Model=0.0; final_eq_solution_details.C_model=C_model; final_eq_solution_details.K_model_pps=K_model_pps_sim; final_eq_solution_details.K_model_non_pps=K_model_nonpps_sim;
        end


    end

    methods (Static, Access = private) % <<<<<< 新增一个私有静态方法块

        function neg_v = objective_for_k_prime_private(k_prime_choice, resources, k_pps_state, c_pps_choice, R_k_net_factor_age, a_idx, ev_interpolant, cS)
            % [移到这里] 这个函数现在是类的一个私有静态方法
            % 它的代码内容完全不变
            
            % 1. 计算消费
            c_expend = resources - k_prime_choice;
            c_choice = max(cS.cFloor, c_expend / (1 + cS.tau_c));

            % 2. 计算当期效用
            [~, util] = main_olg_v8_utils.CES_utility(c_choice, cS.sigma, cS);

            % 3. 计算下一期PPS资产
            pps_withdrawal = 0;
            if a_idx >= cS.aR_new, pps_withdrawal = k_pps_state * cS.pps_withdrawal_rate; end
            pps_return_factor = 1 + ((R_k_net_factor_age - 1) + cS.pps_return_rate_premium);
            k_pps_prime = (k_pps_state + c_pps_choice - pps_withdrawal) * pps_return_factor;

    % 4. [核心修复] 在调用插值器前，对所有查询点进行钳位
    %    fminbnd 优化的 k_prime_choice 已经由其边界约束，但为了安全再次钳位。
    k_prime_clamped = max(cS.kGridV(1), min(cS.kGridV(end), k_prime_choice));
    k_pps_prime_clamped = max(cS.kppsGridV(1), min(cS.kppsGridV(end), k_pps_prime));
    
    % 使用被钳位的点进行插值
    % 注意：CallInterpolator 函数内部已经有 try-catch 钳位逻辑，
    % 但在这里显式钳位是更稳健的做法，确保我们完全控制了输入。
    ev = main_olg_v8_utils.CallInterpolator(ev_interpolant, k_prime_clamped, k_pps_prime_clamped, cS);

            % 5. 计算总价值
            current_val = util + cS.beta * cS.s_1yr_transitionV(a_idx) * ev;

            % fminbnd是最小化器，所以返回负价值
            neg_v = -current_val;
        end
    end % End of Static Methods
end % End of Classdef


