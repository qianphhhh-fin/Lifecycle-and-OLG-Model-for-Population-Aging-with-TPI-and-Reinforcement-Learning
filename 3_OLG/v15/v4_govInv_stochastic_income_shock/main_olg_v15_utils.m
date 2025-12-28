% =========================================================================
% == CLASS: main_olg_v15_utils.m [vPaper.4 - 支出冲击模型]
% == 目的: 包含模型的所有参数、函数和求解器。
% == 核心修正:
% == 1. 模型从“收入冲击”转为“支出冲击”，冲击直接作用于预算约束。
% == 2. VFI求解器、反解函数、聚合函数均已适配支出冲击结构。
% == 3. 养老金b_t的计算被简化，因老年冲击不再直接影响养老金领取。
% == 4. 劳动供给L的求解采用内层迭代，确保市场出清和国民账户平衡。
% =========================================================================

classdef main_olg_v15_utils

    methods (Static)

        function cS = ParameterValues_ExpenditureShock()
            % [vPaper.5 - 政府投资版]
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

            % 其他基本参数
            cS.beta = 0.995; cS.sigma = 10; cS.cFloor = 0.005;
            cS.time_Step = 5;

            % 企业
            cS.alpha = 0.35; cS.ddk = 1 - (1 - 0.02)^cS.time_Step;

            % [!!!!! 核心修改 !!!!!]
            % --- 政府投资参数 ---
            % 公共资本的产出弹性
            cS.gamma = 0.1;
            % 政府投资占总税收的比例
            cS.lambda_g = 0.35;
            % 公共资本折旧率（为简化，设为与私人资本相同）
            cS.ddk_g = cS.ddk;

            cS.tau_k = 0.02; cS.tau_l = 0.06; cS.tau_c = 0.03;
            cS.gov_debt_frac_Y = 0; cS.A = 1.0;

            % --- 年龄组聚合参数 ---
            cS.aD_new = ceil(cS.aD_orig / cS.time_Step);
            cS.aR_new = ceil(cS.aW_orig / cS.time_Step);
            cS.physAgeMap = cell(cS.aD_new, 1);
            for a = 1:cS.aD_new
                startIdx = (a-1)*cS.time_Step + 1;
                endIdx = min(a*cS.time_Step, cS.aD_orig);
                cS.physAgeMap{a} = startIdx:endIdx;
            end

            % 存活率
            d_orig_data = [0.00159,0.00169,0.00174,0.00172,0.00165,0.00156,0.00149,0.00145,0.00145,0.00149,0.00156,0.00163,0.00171,0.00181,0.00193,0.00207,0.00225,0.00246,0.00270,0.00299,0.00332,0.00368,0.00409,0.00454,0.00504,0.00558,0.00617,0.00686,0.00766,0.00865,0.00955,0.01058,0.01162,0.01264,0.01368,0.01475,0.01593,0.01730,0.01891,0.02074,0.02271,0.02476,0.02690,0.02912,0.03143,0.03389,0.03652,0.03930,0.04225,0.04538,0.04871,0.05230,0.05623,0.06060,0.06542,0.07066,0.07636,0.08271,0.08986,0.09788,0.10732,0.11799,0.12895,0.13920,0.14861,0.16039,0.17303,0.18665,0.20194,0.21877,0.23601,0.25289,0.26973,0.28612,0.30128,0.31416,0.32915,0.34450,0.36018];
            s_orig_1yr = 1 - d_orig_data;
            cS.s_pathV = zeros(cS.aD_new, 1);
            for a = 1:(cS.aD_new - 1)
                age_indices = cS.physAgeMap{a};
                s_5yr = 1.0;
                for i = 1:length(age_indices)
                    s_5yr = s_5yr * s_orig_1yr(age_indices(i));
                end
                cS.s_pathV(a) = s_5yr;
            end
            cS.s_pathV(cS.aD_new) = 0;

            % 中国平均预期寿命约为78岁(2021-2023数据)，模型中对应物理年龄索引59
            avg_lifespan_idx = 78;

            % --- [核心变更] 支出冲击参数 ---
            cS.nw = 5; % 常规效率状态数量
            % 年龄节点
            cS.young_shock_start_age = 26;
            cS.young_shock_end_age = 35;
            cS.old_shock_start_age = 65;
            cS.old_shock_end_age = 88;
            % 冲击年化概率峰值
            cS.p_shock_young_peak_annual = 0.03;
            cS.p_shock_old_peak_annual = 0.03;
            % 冲击支出比例 (kappa)
            cS.kappa_young = 0.7;
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

            cS.theta_path = 0.06;
        end

        function cS = generateGrids(cS)
            cS.tgKY = 3; cS.tgWage = (1-cS.alpha)*cS.A*((cS.tgKY/cS.A)^(cS.alpha/(1-cS.alpha)));
            cS.kMin = 0; cS.kMax = 50 * cS.tgWage; cS.kppsMin = 0; cS.kppsMax = cS.kMax / 2;
            power_k = 4;
            if cS.nk > 1, kGridV_temp = cS.kMin + (cS.kMax - cS.kMin) * (linspace(0, 1, cS.nk).^power_k);
            else, kGridV_temp = cS.kMin; end
            cS.kGridV = kGridV_temp(:);
            if cS.nkpps > 1, kppsGridV_temp = cS.kppsMin + (cS.kppsMax - cS.kppsMin) * (linspace(0, 1, cS.nkpps).^power_k);
            else, kppsGridV_temp = cS.kppsMin; end
            cS.kppsGridV = kppsGridV_temp(:);
        end

        function [leGridV_expanded, TrProbM_by_age, leProb1V, nw_expanded] = EarningProcess_AgeDependent(cS)
            % [v.PlateauShock - 平台式冲击版]
            % 核心修改: 将老年冲击的概率生成逻辑，从“斜坡式”(linspace)
            %           修改为与青年冲击完全相同的“平台式”(固定概率)。

            % --- 1. 生成基础AR(1)过程 ---
            lePersistence = 0.9; leShockStd = 0.15^0.5; Tauchen_q = 2.0;
            [leLogGridV_normal, leTrProbM_normal] = main_olg_v15_utils.tauchen(cS.nw, lePersistence, leShockStd, 0, Tauchen_q);
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

        function [cPolM, kPolM, cPpsPolM, valM] = HHSolution_VFI(M_vfi, paramS_vfi, cS_vfi)
            valM = -Inf(cS_vfi.nk, cS_vfi.nkpps, cS_vfi.nw_expanded, cS_vfi.aD_new);
            cPolM = zeros(cS_vfi.nk, cS_vfi.nkpps, cS_vfi.nw_expanded, cS_vfi.aD_new);
            kPolM = zeros(cS_vfi.nk, cS_vfi.nkpps, cS_vfi.nw_expanded, cS_vfi.aD_new);
            cPpsPolM = zeros(cS_vfi.nk, cS_vfi.nkpps, cS_vfi.nw_expanded, cS_vfi.aD_new);
            bV_payg_vfi = zeros(1, cS_vfi.aD_new);
            if cS_vfi.aR_new < cS_vfi.aD_new
                bV_payg_vfi((cS_vfi.aR_new+1):cS_vfi.aD_new) = M_vfi.b_t;
            end
            for a_idx = cS_vfi.aD_new : -1 : 1
                vPrime_kkppse_next = [];
                if a_idx < cS_vfi.aD_new
                    vPrime_kkppse_next = valM(:,:,:,a_idx+1);
                end
                [cPolM(:,:,:,a_idx), kPolM(:,:,:,a_idx), cPpsPolM(:,:,:,a_idx), valM(:,:,:,a_idx)] = ...
                    main_olg_v15_utils.HHSolutionByAge_VFI_ExpenditureShock(a_idx, vPrime_kkppse_next, M_vfi, bV_payg_vfi(a_idx), paramS_vfi, cS_vfi);
            end
        end

        function [cPol_age_q, kPol_age, cPpsPol_age_choice, val_age] = HHSolutionByAge_VFI_ExpenditureShock(a_idx, vPrime_kkppse_next, M_age, b_age_val, paramS_age, cS)
            % [vPaper.4.9 - 暖光遗赠动机版]
            % 核心变更:
            % 1. 在效用函数中加入了“暖光”遗赠动机。
            % 2. 个体在决策时，会考虑其储蓄 k' 在其意外死亡时能带来的效用。
            %    Value = u(c) + β*s*EV' + β*(1-s)*u_beq(k')
            % 3. 在生命最后一期，所有资源都用于遗赠，并从中获得效用。

            val_age = -1e20 * ones(cS.nk, cS.nkpps, cS.nw_expanded);
            cPol_age_q = zeros(cS.nk, cS.nkpps, cS.nw_expanded);
            kPol_age = zeros(cS.nk, cS.nkpps, cS.nw_expanded);
            cPpsPol_age_choice = zeros(cS.nk, cS.nkpps, cS.nw_expanded);

            % --- [修改点 1] ---
            % 在最后一期，死亡是确定的 (s=0)。
            % 如果有遗赠动机，个体不会消费，而是将所有资源留作遗赠。
            if a_idx == cS.aD_new
                [K_grid, ~] = ndgrid(cS.kGridV, cS.kppsGridV);
                pretax_non_capital_income = b_age_val;
                capital_tax = cS.tau_k .* M_age.r_mkt_t .* K_grid;
                total_resources = K_grid .* (1 + M_age.r_mkt_t) + pretax_non_capital_income - capital_tax;

                if cS.phi_bequest > 0
                    % 有遗赠动机：不消费，全部储蓄（即遗赠）
                    k_prime_final = total_resources;
                    c_final = zeros(size(k_prime_final));

                    % 效用完全来自遗赠
                    util_final = main_olg_v15_utils.bequest_utility(k_prime_final, cS);

                    for ie = 1:cS.nw_expanded
                        cPol_age_q(:,:,ie) = c_final;
                        kPol_age(:,:,ie) = k_prime_final; % 记录遗赠额
                        val_age(:,:,ie) = util_final;
                    end

                else
                    % 无遗赠动机（同旧版）：消费掉所有资源
                    c_expend_final = total_resources;
                    final_c = c_expend_final / (1 + cS.tau_c);
                    [~, util_c] = main_olg_v15_utils.CES_utility(final_c, cS.sigma, cS);
                    for ie = 1:cS.nw_expanded
                        cPol_age_q(:,:,ie) = final_c; val_age(:,:,ie) = util_c;
                        kPol_age(:,:,ie) = 0;
                    end
                end
                cPpsPol_age_choice(:,:,:) = 0; % 最后一期没有养老金份额选择
                return;
            end

            effective_discount_factor = (cS.beta ^ cS.time_Step) * cS.s_pathV(a_idx);

            % --- [修改点 2] ---
            % 引入与遗赠效用相关的折现因子
            bequest_discount_factor = (cS.beta ^ cS.time_Step) * (1 - cS.s_pathV(a_idx));

            EV_matrix = zeros(cS.nk, cS.nkpps, cS.nw_expanded);
            if ~isempty(vPrime_kkppse_next)
                transition_matrix_next_age = paramS_age.TrProbM_by_age{a_idx + 1};
                vPrime_reshaped = reshape(vPrime_kkppse_next, [cS.nk * cS.nkpps, cS.nw_expanded]);
                for ie_current = 1:cS.nw_expanded
                    transition_probs = transition_matrix_next_age(ie_current, :);
                    EV_slice = vPrime_reshaped * transition_probs';
                    EV_matrix(:, :, ie_current) = reshape(EV_slice, [cS.nk, cS.nkpps]);
                end
            end
            EV_interpolants = cell(cS.nw_expanded, 1);
            is_pps_disabled = (cS.nkpps == 1);
            for ie_current=1:cS.nw_expanded
                if is_pps_disabled, EV_interpolants{ie_current}=griddedInterpolant(cS.kGridV,squeeze(EV_matrix(:,1,ie_current)),'pchip','none');
                else, EV_interpolants{ie_current}=griddedInterpolant({cS.kGridV, cS.kppsGridV},EV_matrix(:,:,ie_current),'pchip','none'); end
            end

            market_return_factor = 1 + M_age.r_mkt_t;
            for ie = 1:cS.nw_expanded, epsilon_state = paramS_age.leGridV(ie); ev_interpolant = EV_interpolants{ie};
                for ik = 1:cS.nk, for ikpps = 1:cS.nkpps
                        k_state = cS.kGridV(ik);
                        best_val = -1e20; best_k_prime = cS.kMin; best_c = 0;

                        capital_income = k_state * M_age.r_mkt_t;
                        labor_income_gross = 0; pension_benefit = 0;
                        if a_idx <= cS.aR_new
                            labor_income_gross = M_age.w_t*cS.ageEffV_new(a_idx)*epsilon_state;
                        else
                            pension_benefit = b_age_val;
                        end
                        payg_tax = cS.theta_t*labor_income_gross;
                        labor_tax = cS.tau_l*max(0, labor_income_gross - payg_tax);
                        capital_tax = cS.tau_k * capital_income;

                        net_cash_before_shock = k_state * market_return_factor + labor_income_gross + pension_benefit - (payg_tax + labor_tax + capital_tax);

                        shock_expenditure = 0;
                        if ie == cS.nw + 1
                            shock_expenditure = cS.kappa_young * net_cash_before_shock;
                        elseif ie == cS.nw + 2
                            shock_expenditure = cS.kappa_old * net_cash_before_shock;
                        end

                        net_cash_for_c_k_prime = net_cash_before_shock - shock_expenditure;

                        for k_prime_choice = cS.kGridV'
                            c_expend = net_cash_for_c_k_prime - k_prime_choice;

                            if c_expend <= 0
                                break;
                            end

                            c_choice = c_expend / (1 + cS.tau_c);
                            [~,util] = main_olg_v15_utils.CES_utility(c_choice, cS.sigma, cS);

                            if is_pps_disabled, ev = ev_interpolant(k_prime_choice); else, ev = ev_interpolant(k_prime_choice, 0); end
                            if isnan(ev), ev = -1e10; end

                            % --- [修改点 3] ---
                            % 在当期价值计算中加入遗赠效用
                            util_bequest = main_olg_v15_utils.bequest_utility(k_prime_choice, cS);

                            current_val = util + effective_discount_factor * ev + bequest_discount_factor * util_bequest;

                            if current_val > best_val, best_val = current_val; best_c = c_choice; best_k_prime = k_prime_choice; end
                        end

                        if isinf(best_val)||isnan(best_val), val_age(ik,ikpps,ie)=-1e20; kPol_age(ik,ikpps,ie)=cS.kMin; cPpsPol_age_choice(ik,ikpps,ie)=0; cPol_age_q(ik,ikpps,ie)=0;
                        else, val_age(ik,ikpps,ie)=best_val; kPol_age(ik,ikpps,ie)=best_k_prime; cPpsPol_age_choice(ik,ikpps,ie)=0; cPol_age_q(ik,ikpps,ie)=best_c; end
                end,end,end
        end

        function Dist=solve_steady_state_distribution(k_prime_idx,paramS,cS,Z_ss_norm)
            Dist = zeros(cS.nk, cS.nw_expanded, cS.aD_new);
            dist_newborn_ke = zeros(cS.nk, cS.nw_expanded);
            dist_newborn_ke(1, 1:cS.nw) = paramS.leProb1V';
            Dist(:, :, 1) = dist_newborn_ke * Z_ss_norm(1);

            for ia = 1:(cS.aD_new - 1)
                dist_ia_ke = Dist(:,:,ia);
                dist_ia_plus_1_ke = zeros(cS.nk, cS.nw_expanded);
                transition_matrix_next_age = paramS.TrProbM_by_age{ia + 1};
                for ik = 1:cS.nk, for ie = 1:cS.nw_expanded
                        mass_at_state = dist_ia_ke(ik, ie);
                        if mass_at_state < 1e-20, continue; end
                        ik_prime = k_prime_idx(ik, ie, ia);
                        transition_probs_e = transition_matrix_next_age(ie, :);
                        mass_surviving = mass_at_state * cS.s_pathV(ia);
                        dist_ia_plus_1_ke(ik_prime, :) = dist_ia_plus_1_ke(ik_prime, :) + mass_surviving * transition_probs_e;
                end, end
            Dist(:,:,ia+1) = dist_ia_plus_1_ke;
            end
            if abs(sum(Dist, 'all') - 1.0) > 1e-6
                warning('稳态联合分布总和(%.8f)不为1，可能存在会计不一致。', sum(Dist, 'all'));
            end
        end

        function [c_val,tax_val,shock_expenditure,payg_tax,labor_tax,capital_tax,consumption_tax,pension_benefit]=backout_accounting_expenditure_shock(k_now,k_prime,ia,ie,epsilon_val,M_sim,cS)
            % [vPaper.4.8 - 弹性消费下限版]
            % 核心变更: 移除了 c_val = max(cFloor, ...) 的强制修正。
            %           现在函数直接报告由预算约束反解出的实际消费值，
            %           即使它可能低于cFloor。

            if ~isfield(cS, 'theta_t'), cS.theta_t = cS.theta_path(1); end

            market_return_factor = 1 + M_sim.r_mkt_t;
            capital_income = k_now * M_sim.r_mkt_t;
            labor_income_gross = 0; pension_benefit = 0;
            if ia <= cS.aR_new
                labor_income_gross = M_sim.w_t * cS.ageEffV_new(ia) * epsilon_val;
            else
                pension_benefit = M_sim.b_t;
            end
            total_inflow = k_now * market_return_factor + labor_income_gross + pension_benefit;

            payg_tax = cS.theta_t * labor_income_gross;
            capital_tax = cS.tau_k * capital_income;
            labor_tax = cS.tau_l * max(0, labor_income_gross - payg_tax);
            net_cash_before_shock = total_inflow - (payg_tax + capital_tax + labor_tax);

            shock_expenditure = 0;
            if ie == cS.nw + 1
                shock_expenditure = cS.kappa_young * net_cash_before_shock;
            elseif ie == cS.nw + 2
                shock_expenditure = cS.kappa_old * net_cash_before_shock;
            end

            net_cash_for_c_k_prime = net_cash_before_shock - shock_expenditure;

            % [!!!!! 核心修改 !!!!!]
            c_expend_available = net_cash_for_c_k_prime - k_prime;

            % 我们不再用max(cFloor,...)来强制修正消费。
            % 我们相信上游的VFI和离散化步骤会保证c_expend_available是正的。
            c_val = c_expend_available / (1 + cS.tau_c);

            % 如果c_val最终为负或零，说明上游有更严重的问题，让它在后续计算中暴露出来。

            consumption_tax = c_val * cS.tau_c;
            tax_val = labor_tax + capital_tax + consumption_tax;

            % --- 预算平衡检验依然保留，作为最终防线 ---
            total_consumption_expenditure = c_val * (1 + cS.tau_c);
            total_tax_outflow = payg_tax + capital_tax + labor_tax;
            total_outflow = k_prime + total_consumption_expenditure + shock_expenditure + total_tax_outflow;
            budget_gap = total_inflow - total_outflow;

            if abs(budget_gap) > 1e-9
                % 现在，预算缺口应该非常接近于0，因为我们没有“丢弃”任何现金
                error('微观预算约束被违反！状态: (ia=%d, ie=%d, ik=%d), 预算缺口: %.3e', ...
                    ia, ie, find(abs(cS.kGridV - k_now) < 1e-10, 1), budget_gap);
            end
        end

        function [K_agg,C_utility_agg,Tax_agg,Bequest_tax_agg,L_agg,PensionIn_agg,PensionOut_agg,ShockExp_agg]=aggregate_expenditure_shock(Dist,k_prime_idx,M_sim,cS,paramS)
            if ~isfield(cS, 'theta_t'), cS.theta_t = cS.theta_path(1); end
            K_agg=0; C_utility_agg=0; Tax_agg=0; Bequest_tax_agg=0; L_agg=0; PensionIn_agg=0; PensionOut_agg=0; ShockExp_agg=0;

            for ia=1:cS.aD_new
                for ie=1:cS.nw_expanded
                    for ik=1:cS.nk
                        mass = Dist(ik,ie,ia);
                        if mass < 1e-20, continue; end
                        k_now=cS.kGridV(ik); epsilon_val=paramS.leGridV(ie);
                        idx_k_prime=k_prime_idx(ik,ie,ia); k_prime=cS.kGridV(idx_k_prime);

                        [c_val,tax_val,shock_exp,payg_tax,~,~,~,pension_benefit] = ...
                            main_olg_v15_utils.backout_accounting_expenditure_shock(k_now,k_prime,ia,ie,epsilon_val,M_sim,cS);

                        C_utility_agg = C_utility_agg + c_val * mass;
                        Tax_agg = Tax_agg + tax_val * mass;
                        ShockExp_agg = ShockExp_agg + shock_exp * mass;

                        prob_survive = cS.s_pathV(ia);
                        K_agg = K_agg + k_prime * mass * prob_survive;
                        prob_death = 1 - prob_survive;
                        Bequest_tax_agg = Bequest_tax_agg + k_prime * mass * prob_death;

                        if ia<=cS.aR_new, L_agg = L_agg + (cS.ageEffV_new(ia)*epsilon_val) * mass; end
                        PensionIn_agg = PensionIn_agg + payg_tax * mass;
                        PensionOut_agg = PensionOut_agg + pension_benefit * mass;
                    end
                end
            end
        end

        function display_national_accounts_gov_investment(ss, Dist, k_prime_idx, cS, paramS)
            fprintf('\n\n========================================================================\n');
            fprintf('===     国民经济核算详细报告 (政府投资模型版)     ===\n');
            fprintf('========================================================================\n');

            M_sim = struct('r_mkt_t', ss.r_mkt, 'w_t', ss.w, 'b_t', ss.b);

            [~, C_utility_agg, Tax_agg, Bequest_tax_agg, L_agg_check, PensionIn_agg, PensionOut_agg, ShockExp_agg] = ...
                main_olg_v15_utils.aggregate_expenditure_shock(Dist, k_prime_idx, M_sim, cS, paramS);

            % --- 宏观变量定义 ---
            Y_prod = ss.Y_from_production;
            K_p = ss.K_private;
            K_g = ss.K_public;

            Gov_Revenue_total = Tax_agg + Bequest_tax_agg;
            I_g_agg = cS.lambda_g * Gov_Revenue_total;         % 政府投资
            G_c_agg = (1 - cS.lambda_g) * Gov_Revenue_total; % 政府消费

            I_p_agg_gross = cS.ddk * K_p; % 私人部门总投资 = 折旧
            I_total_agg = I_p_agg_gross + I_g_agg; % 经济总投资

            C_total_agg = C_utility_agg + ShockExp_agg; % 总消费

            % 支出法GDP
            Y_exp_actual = C_total_agg + I_total_agg + G_c_agg;

            fprintf('--- A. 宏观产出与支出核算 ---\n');
            fprintf('   生产法 GDP (Y_prod):         %.6f\n', Y_prod);
            fprintf('   支出法 GDP (C+I_total+G_c):  %.6f\n', Y_exp_actual);
            fprintf('   ------------------------------------\n');
            fprintf('   核算误差 (Y_exp - Y_prod):     %.3e (此值应接近0)\n', Y_exp_actual - Y_prod);
            fprintf('   总消费 (C):                  %.6f (占GDP: %.2f%%)\n', C_total_agg, C_total_agg/Y_prod*100);
            fprintf('     - 常规消费 (有作用):       %.6f\n', C_utility_agg);
            fprintf('     - 重大冲击支出 (无作用):   %.6f\n', ShockExp_agg);
            fprintf('   总投资 (I_total = I_p+I_g):  %.6f (占GDP: %.2f%%)\n', I_total_agg, I_total_agg/Y_prod*100);
            fprintf('     - 私人投资 (I_p):          %.6f\n', I_p_agg_gross);
            fprintf('     - 政府投资 (I_g):          %.6f\n', I_g_agg);
            fprintf('   政府消费 (G_c):              %.6f (占GDP: %.2f%%)\n', G_c_agg, G_c_agg/Y_prod*100);

            fprintf('\n--- B. 政府与养老金体系核算 ---\n');
            fprintf('   政府总收入 (T):              %.6f\n', Gov_Revenue_total);
            fprintf('   政府总支出 (G_c+I_g):        %.6f\n', G_c_agg + I_g_agg);
            fprintf('   政府预算平衡 (T - G_c - I_g):%.3e (此值应接近0)\n', Gov_Revenue_total - (G_c_agg + I_g_agg));
            fprintf('   ------------------------------------\n');
            fprintf('   养老金体系平衡 (收入-支出):  %.3e (此值应接近0)\n', PensionIn_agg - PensionOut_agg);
            fprintf('   劳动供给核算误差 (L_agg-L_ss): %.3e\n', L_agg_check - ss.L);

            fprintf('\n--- C. 关键宏观比率 ---\n');
            fprintf('   私人资本产出比 (K_p/Y):      %.4f\n', K_p / Y_prod);
            fprintf('   公共资本产出比 (K_g/Y):      %.4f\n', K_g / Y_prod);
            fprintf('   总资本产出比 (K_total/Y):    %.4f\n', (K_p + K_g) / Y_prod);
            fprintf('\n========================================================================\n');
        end

        function plot_saving_rate_by_age(ss, Dist, k_prime_idx, cS, paramS)
            % =========================================================================
            % == FUNCTION: plot_saving_rate_by_age
            % == 目的: 计算并绘制出模型中各个年龄组的平均储蓄率。
            % == 储蓄率定义: S_rate = (Aggregate k_prime) / (Aggregate Net Cash for C/k')
            % =========================================================================

            fprintf('\n   正在计算各年龄组储蓄率...\n');

            % 初始化用于聚合的向量
            age_saving_agg = zeros(cS.aD_new, 1);
            age_disposable_resources_agg = zeros(cS.aD_new, 1);

            M_sim = struct('r_mkt_t', ss.r_mkt, 'w_t', ss.w, 'b_t', ss.b);

            % 遍历所有状态来计算每个年龄组的总储蓄和总可支配资源
            for ia = 1:cS.aD_new
                for ik = 1:cS.nk
                    for ie = 1:cS.nw_expanded
                        mass = Dist(ik, ie, ia);
                        if mass < 1e-20
                            continue;
                        end

                        % 获取当前状态和决策
                        k_now = cS.kGridV(ik);
                        epsilon_val = paramS.leGridV(ie);
                        idx_k_prime = k_prime_idx(ik, ie, ia);
                        k_prime = cS.kGridV(idx_k_prime);

                        % --- 为该代理人计算可用于消费和储蓄的净现金流 ---
                        % 这部分逻辑与 backout_accounting_... 函数中的预算约束计算一致
                        capital_income = k_now * M_sim.r_mkt_t;
                        labor_income_gross = 0;
                        pension_benefit = 0;
                        if ia <= cS.aR_new
                            labor_income_gross = M_sim.w_t * cS.ageEffV_new(ia) * epsilon_val;
                        else
                            pension_benefit = M_sim.b_t;
                        end

                        payg_tax = cS.theta_path(1) * labor_income_gross;
                        labor_tax = cS.tau_l * max(0, labor_income_gross - payg_tax);
                        capital_tax = cS.tau_k * capital_income;

                        net_cash_before_shock = k_now * (1 + M_sim.r_mkt_t) + labor_income_gross + pension_benefit - (payg_tax + labor_tax + capital_tax);

                        shock_expenditure = 0;
                        if ie == cS.nw + 1
                            shock_expenditure = cS.kappa_young * net_cash_before_shock;
                        elseif ie == cS.nw + 2
                            shock_expenditure = cS.kappa_old * net_cash_before_shock;
                        end

                        net_cash_for_c_k_prime = net_cash_before_shock - shock_expenditure;

                        % 聚合该代理人的储蓄和可支配资源
                        age_saving_agg(ia) = age_saving_agg(ia) + mass * k_prime;
                        age_disposable_resources_agg(ia) = age_disposable_resources_agg(ia) + mass * net_cash_for_c_k_prime;
                    end
                end
            end

            % 计算最终的储蓄率
            saving_rate_by_age = age_saving_agg ./ age_disposable_resources_agg;

            % 处理可能因分母为0导致的NaN（例如最后一期）
            saving_rate_by_age(isnan(saving_rate_by_age)) = 0;

            % --- 绘制图形 ---
            % 创建用于绘图的年龄向量（取年龄区间的中间点）
            age_vector_plot = cS.age1_orig + ((1:cS.aD_new) - 1) * cS.time_Step + floor(cS.time_Step / 2);

            figure('Name', '年龄-储蓄率剖面');
            plot(age_vector_plot, saving_rate_by_age * 100, 'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
            hold on;

            % 绘制一条垂线以标记退休年龄
            y_limits = ylim;
            plot([cS.ageRetire_orig, cS.ageRetire_orig], y_limits, 'k--', 'LineWidth', 1.5, 'DisplayName', sprintf('退休年龄 (%d岁)', cS.ageRetire_orig));

            hold off;
            grid on;
            title('各年龄组的平均储蓄率');
            xlabel('年龄 (岁)');
            ylabel('储蓄率 (%)');
            legend('储蓄率', '退休年龄', 'Location', 'best');
            xlim([cS.age1_orig, cS.ageLast_orig]);
        end

        function [ss, eq_found, Dist, k_prime_idx] = solve_steady_state_iter_K(Z_ss_norm, cS, paramS)
            % [vPaper.4.7 - 两阶段求解：fzero + 松弛法精炼]
            % 核心逻辑:
            % 1. 使用 fzero 快速求解。
            % 2. 如果结果不理想，则在 fzero 解附近启动带阻尼的松弛法进行精炼。

            % --- 第 1 阶段：使用 fzero 进行快速求解 ---
            capital_market_wrapper = @(K_guess) main_olg_v15_utils.fzero_wrapper_K(K_guess, Z_ss_norm, cS, paramS);
            k_bracket=[0.05, 15];
            options=optimset('TolX',1e-9,'Display','iter','MaxIter',500);

            fprintf('\n--- 阶段 1: 启动 fzero 求解器 ---\n');
            [K_eq, fval, exitflag] = fzero(capital_market_wrapper, k_bracket, options);
            fprintf('--- fzero 求解完成 ---\n');

            % --- 第 2 阶段：检查结果并按需启动松弛法精炼 ---
            refinement_tolerance = 1e-6;

            if (exitflag <= 0) || (abs(fval) > refinement_tolerance)
                fprintf('\n‼️ fzero 结果不理想 (exitflag=%d, error=%.3e > %.1e)。\n', exitflag, fval, refinement_tolerance);
                fprintf('--- 阶段 2: 在 K_eq = %.4f 附近启动松弛法精炼 ---\n', K_eq);

                % 松弛法参数
                lambda = 0.0001;  % 松弛因子 (0 < lambda <= 1), 较小值更稳定
                max_iter_refine = 10; % 最大迭代次数
                tol_refine = 1e-9;   % 收敛容忍度

                K_current = K_eq; % 从 fzero 的结果开始迭代

                fprintf('   松弛因子 λ = %.2f\n', lambda);
                fprintf('   Iter.      K_current        K_model_out         Error\n');
                fprintf('   --------------------------------------------------------\n');

                converged = false;
                for i = 1:max_iter_refine
                    % K_model = K_guess - error
                    K_model_out = K_current - capital_market_wrapper(K_current);

                    % 市场误差: 供给 - 需求 (与 fzero 的定义相反，但更直观)
                    current_error = K_model_out - K_current;

                    fprintf('   %4d    %12.6f    %12.6f    %15.4e\n', i, K_current, K_model_out, current_error);

                    if abs(current_error) < tol_refine
                        fprintf('   松弛法在 %d 次迭代后收敛！\n', i);
                        converged = true;
                        K_current = K_model_out; % 确保最终值为收敛点
                        break;
                    end

                    % 应用松弛法更新
                    K_current = (1 - lambda) * K_current + lambda * K_model_out;
                end

                if ~converged
                    warning('松弛法在 %d 次迭代后未能收敛。', max_iter_refine);
                end

                % 比较 fzero 的解和松弛法的解，取更优者
                final_error_relaxation = K_current - (K_current - capital_market_wrapper(K_current));
                if abs(final_error_relaxation) < abs(fval)
                    fprintf('✅ 松弛法找到了更优的解。\n');
                    fprintf('   原解: K=%.8f, 误差=%.3e\n', K_eq, fval);
                    fprintf('   新解: K=%.8f, 误差=%.3e\n', K_current, -final_error_relaxation); % 变回fzero误差定义
                    K_eq = K_current;
                    fval = -final_error_relaxation;
                    exitflag = 98; % 自定义flag，表示由松弛法得到
                else
                    fprintf('⚠️ 松弛法未能找到更优的解。仍采用 fzero 的原始结果。\n');
                end
            end

            eq_found = (exitflag > 0 || abs(fval) < refinement_tolerance);
            if ~eq_found, warning('iter_K 在两阶段求解后仍未能找到满意的均衡解'); ss=struct(); Dist=[]; k_prime_idx=[]; return; end

            fprintf('\n   最终求解结果: 均衡资本 K_eq = %.8f, 最终K误差 = %.3e\n', K_eq, fval);

            % --- 绘图与报告模块 (保持不变) ---
            % fprintf('\n   正在为诊断图计算函数值...\n');
            % plot_grid_points = 50;
            % k_plot_grid = linspace(k_bracket(1), k_bracket(2), plot_grid_points);
            % error_plot_values = zeros(size(k_plot_grid));
            % for i = 1:plot_grid_points
            %     error_plot_values(i) = capital_market_wrapper(k_plot_grid(i));
            % end
            %
            % figure('Name', '资本市场出清 - 求解器诊断图');
            % hold on;
            % plot(k_plot_grid, error_plot_values, 'b-', 'LineWidth', 2.5, 'DisplayName', '资本市场误差 (K_{demand} - K_{supply})');
            % line(k_bracket, [0 0], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1.0, 'DisplayName', '均衡线 (误差=0)');
            %
            % if exitflag == 98 % 如果是松弛法得到的解
            %     plot(K_eq, 0, 'ms', 'MarkerSize', 14, 'MarkerFaceColor', 'm', 'DisplayName', sprintf('最终均衡解 (松弛法) K_{eq} = %.4f', K_eq));
            % else
            %     plot(K_eq, 0, 'gp', 'MarkerSize', 16, 'MarkerFaceColor', 'g', 'DisplayName', sprintf('最终均衡解 (fzero) K_{eq} = %.4f', K_eq));
            % end
            %
            % plot(k_bracket, [capital_market_wrapper(k_bracket(1)), capital_market_wrapper(k_bracket(2))], 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'DisplayName', '初始求解区间');
            %
            % hold off;
            % grid on;
            % title('资本市场出清 - 求解器诊断图');
            % xlabel('资本存量 (K)');
            % ylabel('资本市场误差 (需求 - 供给)');
            % legend('show', 'Location', 'best');
            % fprintf('   诊断图已成功生成。\n\n');

            % --- 后续计算 ---
            [~, ~, ss, Dist, k_prime_idx] = main_olg_v15_utils.calculate_net_saving(K_eq, Z_ss_norm, cS, paramS);

            if any(Dist(:) > 0)
                mass_by_k = squeeze(sum(Dist, [2, 3]));
                max_occupied_ik = find(mass_by_k > 1e-20, 1, 'last');
                if isempty(max_occupied_ik), max_occupied_k_val = 0;
                else, max_occupied_k_val = cS.kGridV(max_occupied_ik); end
            else, max_occupied_k_val = 0; end

            fprintf('   资产网格上限 (kMax):                 %.4f\n', cS.kMax);
            fprintf('   稳态分布中达到的最大资本:          %.4f\n', max_occupied_k_val);
            if max_occupied_k_val >= cS.kMax * 0.99
                warning('最大占用资本已达到或非常接近网格上限kMax！可能需要调高kMax。');
            end
        end

        function [ss, eq_found, Dist, k_prime_idx] = solve_steady_state_iter_Kg(Z_ss_norm, cS, paramS)
            % [vPaper.5 - 政府投资版]
            % 核心逻辑:
            % 1. 使用 fsolve 求解一个包含两个未知数 (K_p, K_g) 和两个方程的系统。
            %    - Eq1: 私人资本市场出清 (K_p_supply = K_p_demand)
            %    - Eq2: 公共资本稳态 (I_g = ddk_g * K_g)

            % 定义匿名函数，作为 fsolve 的目标
            system_wrapper = @(x) main_olg_v15_utils.system_of_equations_Kg(x, Z_ss_norm, cS, paramS);

            % 为 K_p 和 K_g 提供初始猜测值
            k_p_guess_initial = 3.5; % 基于一般K/Y比率的合理猜测
            k_g_guess_initial = 1.0;
            x0 = [k_p_guess_initial, k_g_guess_initial];

            % 设置 fsolve 选项
            options = optimoptions('fsolve', 'Display', 'iter', 'TolFun', 1e-9, 'TolX', 1e-9, 'MaxIterations', 500);

            fprintf('\n--- 启动 fsolve 求解器 (求解 [K_p, K_g]) ---\n');
            [x_eq, fval, exitflag] = fsolve(system_wrapper, x0, options);
            fprintf('--- fsolve 求解完成 ---\n');

            eq_found = (exitflag > 0);
            if ~eq_found
                warning('fsolve 未能找到均衡解 (exitflag = %d)', exitflag);
                % ss=struct(); Dist=[]; k_prime_idx=[];
                % return;
            end

            K_p_eq = x_eq(1);
            K_g_eq = x_eq(2);

            fprintf('\n   最终求解结果: K_p = %.8f, K_g = %.8f\n', K_p_eq, K_g_eq);
            fprintf('   最终方程误差: error_Kp = %.3e, error_Kg = %.3e\n', fval(1), fval(2));

            % --- 求解成功后，进行最后一次计算以获得所有宏观变量 ---
            [~, ss, Dist, k_prime_idx] = main_olg_v15_utils.calculate_aggregates_for_Kg(K_p_eq, K_g_eq, Z_ss_norm, cS, paramS);

            % ... (此处可以保留原有的网格上限检查等)
        end

        function F_error = system_of_equations_Kg(x, Z_ss_norm, cS, paramS)
            % [vPaper.5 - 政府投资版]
            % 目的: 为 fsolve 提供一个误差向量 F_error = [error_Kp; error_Kg]

            K_p_guess = x(1);
            K_g_guess = x(2);

            % 调用核心计算函数，得到模型隐含的宏观聚合量
            [K_p_model, ss] = main_olg_v15_utils.calculate_aggregates_for_Kg(K_p_guess, K_g_guess, Z_ss_norm, cS, paramS);

            % 方程1的误差: 私人资本供给 - 私人资本需求
            % (fsolve 需要 error=0，所以我们定义为 guess - model_outcome)
            error_Kp = K_p_guess - K_p_model;

            % 方程2的误差: 公共投资 - 公共资本折旧
            Gov_Revenue_total = ss.Regular_tax + ss.Bequest_tax;
            I_g_model = cS.lambda_g * Gov_Revenue_total;
            Depreciation_g_model = cS.ddk_g * K_g_guess;
            error_Kg = I_g_model - Depreciation_g_model;

            F_error = [error_Kp; error_Kg];
        end

        function K_error_out = fzero_wrapper_K(K_guess, Z_ss_norm, cS, paramS)
            K_model = main_olg_v15_utils.calculate_net_saving(K_guess, Z_ss_norm, cS, paramS);
            K_error_out = K_guess - K_model;
        end

        function [K_p_model_out, ss, Dist, k_prime_idx] = calculate_aggregates_for_Kg(K_p_guess, K_g_guess, Z_ss_norm, cS, paramS)
            % [vPaper.5 - 政府投资版]
            % 核心修改:
            % 1. 接收 K_p 和 K_g 作为输入
            % 2. 调用已更新的 get_prices_at_t
            % 3. 返回完整的ss结构体，供误差计算和最终报告使用

            if K_p_guess <= 0, K_p_guess = 1e-8; end
            if K_g_guess <= 0, K_g_guess = 1e-8; end
            A_ss = cS.A; theta_ss = cS.theta_path(1);

            L_iter_tol = 1e-7; L_iter_max = 50; L_damping = 0.5;

            % 劳动供给的内层循环 (与之前逻辑相同)
            mean_e_by_age = zeros(cS.aD_new,1);
            e_dist_by_age = zeros(cS.aD_new, cS.nw_expanded);
            e_dist_by_age(1, 1:cS.nw) = paramS.leProb1V';
            mean_e_by_age(1) = e_dist_by_age(1, :) * paramS.leGridV(:);
            for ia = 1:(cS.aD_new - 1)
                e_dist_by_age(ia+1, :) = e_dist_by_age(ia, :) * paramS.TrProbM_by_age{ia+1};
                mean_e_by_age(ia+1) = e_dist_by_age(ia+1, :) * paramS.leGridV(:);
            end
            L_guess = 0;
            for ia = 1:cS.aR_new
                L_guess = L_guess + Z_ss_norm(ia) * cS.ageEffV_new(ia) * mean_e_by_age(ia);
            end

            for iter = 1:L_iter_max
                % [核心修改] 调用新的价格函数
                M_prices = main_olg_v15_utils.get_prices_at_t(K_p_guess, K_g_guess, L_guess, A_ss, cS);
                M_for_hh = M_prices;

                total_wage_bill = M_prices.w_t * L_guess;
                mass_retirees_ss = sum(Z_ss_norm((cS.aR_new+1):end));
                if mass_retirees_ss > 1e-9
                    M_for_hh.b_t = (theta_ss * total_wage_bill) / mass_retirees_ss;
                else, M_for_hh.b_t = 0; end

                cS_ss = cS; cS_ss.pps_active = false; cS_ss.theta_t = theta_ss;
                if ~cS_ss.pps_active, cS_ss.nkpps = 1; cS_ss.npps = 1; cS_ss = main_olg_v15_utils.generateGrids(cS_ss); end

                [~, kPolM, ~, ~] = main_olg_v15_utils.HHSolution_VFI(M_for_hh, paramS, cS_ss);
                k_prime_idx = main_olg_v15_utils.get_policy_index_matrix(kPolM, cS_ss);
                Dist = main_olg_v15_utils.solve_steady_state_distribution(k_prime_idx, paramS, cS_ss, Z_ss_norm);

                L_model = 0;
                for ia = 1:cS.aR_new, for ie = 1:cS.nw_expanded, for ik = 1:cS.nk
                            mass = Dist(ik, ie, ia);
                            if mass > 0
                                epsilon_val = paramS.leGridV(ie);
                                L_model = L_model + (cS.ageEffV_new(ia) * epsilon_val) * mass;
                            end
                end,end,end

        L_error = abs(L_model - L_guess);
        if L_error < L_iter_tol, break; end
        L_guess = L_damping * L_guess + (1 - L_damping) * L_model;
        if iter == L_iter_max, warning('劳动供给内层循环在 %d 次迭代后未收敛. L_error = %.3e', L_iter_max, L_error); end
            end

            % 聚合宏观变量
            [K_p_model_out, C_utility_final, Tax_final, Bequest_tax_final, ~, ~, ~, ~] = ...
                main_olg_v15_utils.aggregate_expenditure_shock(Dist, k_prime_idx, M_for_hh, cS_ss, paramS);

            % 填充完整的 ss 结构体，供外部函数使用
            ss = struct();
            ss.K_private = K_p_guess;
            ss.K_public = K_g_guess;
            ss.K_total = K_p_guess + K_g_guess;
            ss.L = L_guess;
            ss.Y_from_production = M_prices.Y_t;
            ss.w = M_prices.w_t;
            ss.r_mkt = M_prices.r_mkt_t;
            ss.b = M_for_hh.b_t;
            ss.Bequest_tax = Bequest_tax_final;
            ss.Regular_tax = Tax_final;
        end

        function k_prime_idx=get_policy_index_matrix(kPolM,cS)
            % [vPaper.4.1 - 已修正离散化规则]
            % 修正: 更改了离散化方法。不再寻找“最近”的网格点，
            %       而是寻找“可负担的最高”网格点，以确保预算约束在
            %       离散化后永不被违反，从而避免 c_expend < 0。
            k_prime_idx=zeros(cS.nk,cS.nw_expanded,cS.aD_new,'uint16');
            for ia=1:cS.aD_new
                for ie=1:cS.nw_expanded
                    for ik=1:cS.nk
                        % 获取连续空间下的最优储蓄决策
                        k_prime_continuous = kPolM(ik,1,ie,ia);

                        % 寻找所有小于等于该最优决策的网格点的索引
                        affordable_indices = find(cS.kGridV <= k_prime_continuous);

                        if isempty(affordable_indices)
                            % 如果没有任何网格点可负担（例如k_prime_continuous为负）
                            % 则强制选择最低的资产水平（通常是k=0，索引为1）
                            idx = 1;
                        else
                            % 在所有可负担的网格点中，选择最高的那个
                            idx = affordable_indices(end);
                        end
                        k_prime_idx(ik,ie,ia) = idx;
                    end
                end
            end
        end

        function M_prices = get_prices_at_t(K_p, K_g, L, A_t, cS)
            % [vPaper.5.1 - 政府投资修正版]
            % 核心修改: 将公共资本 K_g 作为外部性，修正生产函数，确保国民账户平衡。
            %           Y = (A * K_g^gamma) * K_p^alpha * L^(1-alpha)

            if K_p <= 0, K_p = 1e-8; end; if L <= 0, L = 1e-8; end; if K_g <= 0, K_g = 1e-8; end;

            % 有效的总要素生产率，受公共资本影响
            A_effective = A_t .* (K_g.^cS.gamma);

            % 生产函数，现在对私人要素是常数规模报酬
            Y_period = A_effective .* (K_p.^cS.alpha) .* (L.^(1-cS.alpha));

            % 私人资本的边际产出
            MPK_p_period = cS.alpha .* Y_period ./ K_p;

            % 工资 (现在占总产出的 1-alpha 份额)
            w_t = (1-cS.alpha) .* Y_period ./ L;

            % 私人资本的市场回报率
            r_mkt_t = MPK_p_period - cS.ddk;

            M_prices = struct('K_p', K_p, 'K_g', K_g, 'L_t', L, 'Y_t', Y_period, 'w_t', w_t, 'r_mkt_t', r_mkt_t);
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

        function ky_error = calibrate_ky_target_objective(p_common, cS, paramS, Z_ss_norm, ky_target)
            % =========================================================================
            % == FUNCTION: calibrate_ky_target_objective
            % == 目的: 这是一个用于fzero校准过程的目标函数包装器。
            % ==      它接收一个待校准的冲击概率(p_common)，运行整个模型求解流程，
            % ==      然后返回模型生成的K/Y与目标K/Y之间的误差。
            % =========================================================================

            fprintf('\n================== CALIBRATION STEP: Testing p_common = %.6f ==================\n', p_common);

            % --- 1. 创建参数的本地副本，以避免干扰外部fzero循环 ---
            cS_local = cS;
            paramS_local = paramS;

            % --- 2. 应用当前fzero迭代的冲击概率参数 ---
            cS_local.p_shock_young_peak_annual = p_common;
            cS_local.p_shock_old_peak_annual   = p_common;

            % --- 3. [!! 关键步骤 !!] 重新生成依赖于冲击概率的转移矩阵 ---
            % 因为冲击概率改变了，所以智能体的预期也必须改变，这会体现在年龄依赖的转移矩阵中。
            [paramS_local.leGridV, paramS_local.TrProbM_by_age, paramS_local.leProb1V, cS_local.nw_expanded] = ...
                main_olg_v15_utils.EarningProcess_AgeDependent(cS_local);
            % 关闭单个年龄组的详细概率输出，以保持校准过程的输出整洁
            close(findobj('type','figure','name','支出冲击概率路径')); % 可选：如果EarningProcess生成图形则关闭

            % --- 4. 调用稳态求解器，并处理可能的失败 ---
            % 在一个try-catch块中运行，以防止内部求解失败导致整个校准崩溃。
            ss = struct();
            eq_found = false;
            try
                % 为了让校准过程更快、更稳定，我们在校准步骤中不使用松弛法精炼，
                % 并且关闭其内部迭代的显示。
                % 我们通过修改solve_steady_state_iter_K使其接受一个选项来做到这一点。
                % (注意：此处的实现是直接调用fzero，与您代码中的solve_steady_state_iter_K逻辑一致)

                capital_market_wrapper = @(K_guess) main_olg_v15_utils.fzero_wrapper_K(K_guess, Z_ss_norm, cS_local, paramS_local);
                k_bracket=[0.05, 15];
                options_fzero_inner = optimset('TolX', 1e-7, 'Display', 'off'); % 关闭内部求解器的显示

                [K_eq, ~, exitflag] = fzero(capital_market_wrapper, k_bracket, options_fzero_inner);

                if exitflag > 0
                    eq_found = true;
                    % 如果找到了均衡K, 我们需要调用一次核心计算函数来得到完整的ss结构体
                    [~, ~, ss] = main_olg_v15_utils.calculate_net_saving(K_eq, Z_ss_norm, cS_local, paramS_local);
                end

            catch ME
                fprintf('   ⚠️ 内部求解器发生错误: %s\n', ME.message);
                eq_found = false;
            end

            % 如果内部求解失败，返回一个大的惩罚误差，引导fzero向另一个方向搜索
            if ~eq_found || isempty(fieldnames(ss))
                fprintf('   ⚠️ 内部K求解失败或未收敛，返回一个大的惩罚误差。\n');
                % 如果当前p_common太高导致K/Y过高，我们返回一个正误差；反之亦然
                % 假设K/Y与p_common正相关，这是一个安全的假设
                if p_common > 0.1 % 猜测一个比较高的概率值
                    ky_error = 100; % 返回大的正误差，告诉fzero要降低p
                else
                    ky_error = -100; % 返回大的负误差，告诉fzero要提高p
                end
                return;
            end

            % --- 5. 计算并返回误差 ---
            model_ky_ratio = ss.K_physical / ss.Y_from_production;
            ky_error = model_ky_ratio - ky_target; % 定义误差

            fprintf('>>> p_common = %.6f  -->  Model K/Y = %.4f  -->  Error = %.4f\n', p_common, model_ky_ratio, ky_error);
        end
    end
end
