% --- START OF FILE main_olg_v2_utils.m ---

classdef main_olg_v2_utils
    % 工具函数 for OLG 模型
    % 基于简化模型PAYG (外生theta, 内生b), 引入PPS决策 (新版)
    % - VFI over 5-year groups using ANNUAL beta/prices/survival
    % - Simulation and Aggregation using ANNUAL results/masses (tracks k and k_pps)
    % - Population dynamics still determine steady-state GROUP mass distribution

    methods (Static)

        % =====================================================================
        % == 人口动态函数 (与原始简化模型的 utils.m 相同) =====================
        % =====================================================================
        function popS = initPopulation(cS)
            popS.Z = zeros(cS.aD_new, 1);
            initial_total = sum(cS.initial_pop);
            if initial_total > 0 && length(cS.initial_pop) == cS.aD_new
                 popS.Z(:, 1) = cS.initial_pop(:) / initial_total * 100;
            else
                warning('Initial population mismatch or zero sum. Setting uniform.');
                 popS.Z(:, 1) = 100 / cS.aD_new;
            end
            popS.totalPop = sum(popS.Z(:, 1));
            if popS.totalPop > 1e-9, popS.ageDist = popS.Z(:, 1) / popS.totalPop;
            else, popS.ageDist = zeros(cS.aD_new, 1); end
            popS.initialAgeDist = popS.ageDist;
            fprintf('Initial GROUP population set. Total=%.2f\n', popS.totalPop);
        end

        function popS = populationDynamics(popS, cS)
            max_periods_sim = cS.max_periods;
            Z_history = zeros(cS.aD_new, max_periods_sim + 1);
            totalPop_history = zeros(1, max_periods_sim + 1);
            ageDist_history = zeros(cS.aD_new, max_periods_sim + 1);
            Z_history(:, 1) = popS.Z(:, 1);
            totalPop_history(1) = popS.totalPop(1);
            ageDist_history(:, 1) = popS.ageDist(:, 1);
            fprintf('Population dynamics starting (Groups, Max Periods = %d)...\n', max_periods_sim);
            bgp_reached_flag = false; actual_periods_run = max_periods_sim;
            for t = 1:max_periods_sim
                if mod(t, 10) == 0 || t == 1, fprintf('  Simulating population period %d (groups)\n', t); end
                Z_current_period = Z_history(:, t); Z_next_period = zeros(cS.aD_new, 1);
                 time_varying_growth_rate = 0; if t < 6, time_varying_growth_rate = -0.01 - 0.003 * t; else, time_varying_growth_rate = -0.03 - 0.004 * min(t-6, 10); end
                 Z_next_period(1) = Z_current_period(1) * (1 + time_varying_growth_rate); Z_next_period(1) = max(0, Z_next_period(1));
                for a = 2:cS.aD_new, survival_prob = 0; if (a-1) >= 1 && (a-1) <= length(cS.survivalProbV_popdyn), survival_prob = cS.survivalProbV_popdyn(a-1); end; Z_next_period(a) = Z_current_period(a-1) * survival_prob; Z_next_period(a) = max(0, Z_next_period(a)); end
                Z_history(:, t+1) = Z_next_period; totalPop_history(t+1) = sum(Z_next_period);
                if totalPop_history(t+1) > 1e-9, ageDist_history(:, t+1) = Z_next_period / totalPop_history(t+1); else, ageDist_history(:, t+1) = 0; totalPop_history(t+1) = 0; end
                current_check_period = t + 1;
                if current_check_period >= cS.bgp_window + 1, stable = true; for w = 1:cS.bgp_window, hist_idx1 = current_check_period - w + 1; hist_idx2 = current_check_period - w; if hist_idx1 > 0 && hist_idx2 > 0 && hist_idx1 <= size(ageDist_history,2) && hist_idx2 <= size(ageDist_history,2), change = norm(ageDist_history(:, hist_idx1) - ageDist_history(:, hist_idx2)); if change >= cS.bgp_tolerance, stable = false; break; end; else stable = false; break; end; end; if stable, fprintf('\nPopulation steady state (groups) reached at period %d.\n', t); bgp_reached_flag = true; actual_periods_run = t; break; end; end
            end
            final_period_idx = min(actual_periods_run + 1, size(Z_history,2)); popS.Z = Z_history(:, 1:final_period_idx); popS.totalPop = totalPop_history(1:final_period_idx); popS.ageDist = ageDist_history(:, 1:final_period_idx);
            depRatio_history = zeros(1, actual_periods_run);
             for th = 1:actual_periods_run, Z_t = Z_history(:, th); working_pop = sum(Z_t(1:cS.aR_new)); retired_pop = sum(Z_t(cS.aR_new+1:end)); if working_pop > 1e-9, depRatio_history(th) = retired_pop / working_pop; else depRatio_history(th) = inf; end; end
            popS.dependencyRatio = depRatio_history;
            fprintf('Population dynamics finished. Periods run: %d. BGP reached: %d\n', actual_periods_run, bgp_reached_flag);
            if ~bgp_reached_flag, fprintf('Warning: Population SS not reached.\n'); end
        end

        function [Z_ss, dependency_ratio_ss, bgp_reached, bgp_period] = detectSteadyStatePopulation(popS, cS)
            % (与之前版本相同)
             actual_periods_in_data = size(popS.Z, 2); bgp_reached = false; bgp_period = actual_periods_in_data - 1;
             if actual_periods_in_data < cS.bgp_window + 1, fprintf('Pop sim too short for SS check.\n');
             else
                 fprintf('Checking for population SS (groups, last %d periods)...\n', cS.bgp_window);
                 for t_check_end = actual_periods_in_data : -1 : cS.bgp_window + 1, stable = true; for w = 0 : (cS.bgp_window - 1), idx1 = t_check_end - w; idx2 = t_check_end - w - 1; if idx1 > 0 && idx2 > 0 && idx1 <= size(popS.ageDist,2) && idx2 <= size(popS.ageDist,2), change = norm(popS.ageDist(:, idx1) - popS.ageDist(:, idx2)); if change >= cS.bgp_tolerance, stable = false; break; end; else stable = false; break; end; end; if stable, bgp_reached = true; bgp_period = t_check_end - 1; fprintf('Pop SS detected period %d.\n', bgp_period); break; end; end
                  if ~bgp_reached, fprintf('Pop SS not detected.\n'); end
             end
             ss_data_index = min(bgp_period + 1, size(popS.Z, 2)); Z_ss = popS.Z(:, ss_data_index); Z_ss_norm = zeros(cS.aD_new, 1);
             if sum(Z_ss) > 1e-9, Z_ss_norm = Z_ss / sum(Z_ss); end
             valid_dep_ratio_index = min(bgp_period, length(popS.dependencyRatio));
            if isfield(popS, 'dependencyRatio') && ~isempty(popS.dependencyRatio) && valid_dep_ratio_index > 0 && valid_dep_ratio_index <= length(popS.dependencyRatio)
                 dependency_ratio_ss = popS.dependencyRatio(valid_dep_ratio_index);
            else, working_pop_ss = sum(Z_ss(1:cS.aR_new)); retired_pop_ss = sum(Z_ss(cS.aR_new+1:end)); if working_pop_ss > 1e-9, dependency_ratio_ss = retired_pop_ss / working_pop_ss; else dependency_ratio_ss = inf; end; if (~isfield(popS, 'dependencyRatio') || isempty(popS.dependencyRatio)) && bgp_period > 0, warning('Dependency ratio history not found/too short, re-calculated.');end;end
             figure('Name', 'main_olg_v2: Initial vs SS GROUP Population'); hold on; group_indices = 1:cS.aD_new;
             if isfield(popS, 'initialAgeDist') && ~isempty(popS.initialAgeDist), bar(group_indices - 0.2, popS.initialAgeDist * 100, 0.4, 'b', 'DisplayName', 'Initial'); else warning('Initial age dist not found.'); end
             bar(group_indices + 0.2, Z_ss_norm * 100, 0.4, 'r', 'DisplayName', sprintf('SS (P%d)', bgp_period)); hold off;
             xlabel(sprintf('Age Group (1-%d)', cS.aD_new)); ylabel('Percent of Pop (%)'); title('main_olg_v2: Initial vs SS GROUP Pop Dist'); legend('Location', 'best'); xticks(group_indices); grid on; drawnow;
             fprintf('Plotted initial vs SS GROUP pop dist.\n');
        end

        % =====================================================================
        % == 经济函数 (简化模型PAYG + PPS) ===================================
        % =====================================================================

        function cS = ParameterValues_HuggettStyle()
            % 参数设定: 简化模型PAYG (外生theta), 并加入PPS参数
            % (与你之前能收敛的简化模型参数一致，并加入PPS相关参数)

            %% 人口统计与分组参数
            cS.age1_orig = 20; cS.ageLast_orig = 98; cS.ageRetire_orig = 65; cS.popGrowth_orig = 0.012;
            cS.aD_orig = cS.ageLast_orig - cS.age1_orig + 1; cS.aR_idx_orig = cS.ageRetire_orig - cS.age1_orig + 1;
            cS.aW_orig = cS.aR_idx_orig - 1; cS.physAgeV_orig = (cS.age1_orig : cS.ageLast_orig)';
            cS.d_orig = [0.00159,0.00169,0.00174,0.00172,0.00165,0.00156,0.00149,0.00145,0.00145,0.00149,0.00156,0.00163,0.00171,0.00181,0.00193,0.00207,0.00225,0.00246,0.00270,0.00299,0.00332,0.00368,0.00409,0.00454,0.00504,0.00558,0.00617,0.00686,0.00766,0.00865,0.00955,0.01058,0.01162,0.01264,0.01368,0.01475,0.01593,0.01730,0.01891,0.02074,0.02271,0.02476,0.02690,0.02912,0.03143,0.03389,0.03652,0.03930,0.04225,0.04538,0.04871,0.05230,0.05623,0.06060,0.06542,0.07066,0.07636,0.08271,0.08986,0.09788,0.10732,0.11799,0.12895,0.13920,0.14861,0.16039,0.17303,0.18665,0.20194,0.21877,0.23601,0.25289,0.26973,0.28612,0.30128,0.31416,0.32915,0.34450,0.36018];
            if length(cS.d_orig) ~= cS.aD_orig, error('d_orig 长度不匹配'); end
            cS.s_orig = 1 - cS.d_orig;
            ageMassV_orig_static = ones(1, cS.aD_orig);
            for i = 2 : cS.aD_orig, ageMassV_orig_static(i) = cS.s_orig(i-1) * ageMassV_orig_static(i-1) / (1 + cS.popGrowth_orig); end
            if sum(ageMassV_orig_static) > 1e-9, cS.ageMassV_orig_static = ageMassV_orig_static ./ sum(ageMassV_orig_static); else, cS.ageMassV_orig_static(:) = 1/cS.aD_orig; end
            cS.yearStep = 5; cS.aD_new = ceil(cS.aD_orig / cS.yearStep); cS.aR_new = ceil(cS.aW_orig / cS.yearStep);
            cS.physAgeMap = cell(cS.aD_new, 1); for a = 1:cS.aD_new, startIdx = (a-1)*cS.yearStep + 1; endIdx = min(a*cS.yearStep, cS.aD_orig); cS.physAgeMap{a} = startIdx:endIdx; end
            cS.physAgeV_new = zeros(cS.aD_new, 1); for a = 1:cS.aD_new, cS.physAgeV_new(a) = cS.physAgeV_orig(cS.physAgeMap{a}(1)); end
            cS.s_1yr_transitionV = zeros(cS.aD_new, 1); for a = 1:(cS.aD_new - 1), lastYearIdx = cS.physAgeMap{a}(end); if lastYearIdx < cS.aD_orig, cS.s_1yr_transitionV(a) = cS.s_orig(lastYearIdx); else, cS.s_1yr_transitionV(a) = 0; end; end; cS.s_1yr_transitionV(cS.aD_new) = 0;
            cS.ageMassV_new_static = zeros(1, cS.aD_new); for a = 1:cS.aD_new, cS.ageMassV_new_static(a) = sum(cS.ageMassV_orig_static(cS.physAgeMap{a})); end
            if sum(cS.ageMassV_new_static) > 1e-9, cS.ageMassV_new_static = cS.ageMassV_new_static ./ sum(cS.ageMassV_new_static); else, cS.ageMassV_new_static(:) = 1/cS.aD_new; end
            cS.retireMass_new = sum(cS.ageMassV_new_static(cS.aR_new + 1 : cS.aD_new));

            %% 人口动态参数
            cS.initial_pop = [76.2, 86.4, 113.8, 98.6, 86.6, 102.7, 112.0, 99.0, 64.0, 66.9, 44.1, 25.4, 14.9, 6.8, 1.7, 0.2];
            if length(cS.initial_pop) ~= cS.aD_new, cS.initial_pop = ones(cS.aD_new,1) * (100/cS.aD_new); warning('initial_pop mismatch, using uniform'); end
            beta_surv_pop = [0.998,0.996,0.994,0.992,0.988,0.984,0.980,0.976,0.970,0.960,0.945,0.920,0.880,0.800,0.680];
            if length(beta_surv_pop) ~= cS.aD_new - 1, error('beta_surv_pop length incorrect'); end
            cS.survivalProbV_popdyn = [beta_surv_pop(:)', 0];
            cS.bgp_tolerance = 0.001; cS.bgp_window = 5; cS.max_periods = 50;

            %% 家庭参数 (年度)
            cS.sigma      = 1.5;
            cS.beta       = 1.011; % **使用能收敛的简化模型的beta**
            cS.cFloor     = 0.05;
            cS.nSim       = 5e4;

            %% 技术参数 (年度)
            cS.A          = 0.895944; cS.alpha = 0.36; cS.ddk = 0.06;

            %% 社会保障参数 (年度) - 简化模型设定
            cS.theta      = 0.20; % **使用能收敛的简化模型的theta**

            %% 劳动禀赋参数
            cS.leSigma1 = 0.38^0.5; cS.leShockStd = 0.045^0.5; cS.lePersistence = 0.96; cS.leWidth = 4; cS.nw = 18;

            %% 网格参数 (kGridV代表非PPS资产的网格)
            cS.tgKY = 3; cS.tgWage = (1-cS.alpha)*cS.A*((cS.tgKY/cS.A)^(cS.alpha/(1-cS.alpha)));
            cS.nk = 50; cS.kMin = 0; cS.kMax = 100 * cS.tgWage;
            power = 1.5; kGridV = cS.kMin + (cS.kMax - cS.kMin) * (linspace(0, 1, cS.nk).^power);
            kGridV(1)=cS.kMin; cS.kGridV = kGridV(:);

            %% 年龄效率剖面
            ageEffV_orig_temp = zeros(100, 1);
            ageEffV_orig_temp(20:72) = [linspace(0.3,1.5,36-20+1), 1.5.*ones(1,47-37+1), linspace(1.5,0.2,65-48+1), linspace(0.18,0,72-66+1)];
            cS.ageEffV_orig = ageEffV_orig_temp(cS.age1_orig : cS.ageLast_orig);
            if length(cS.ageEffV_orig) ~= cS.aD_orig, error('ageEffV_orig length mismatch'); end
            cS.ageEffV_new = zeros(cS.aD_new, 1); for a = 1:cS.aD_new, cS.ageEffV_new(a) = mean(cS.ageEffV_orig(cS.physAgeMap{a})); end

            %% 个人养老金计划 (PPS) 参数
            cS.pps_active = true; % **保持PPS激活，但下面的frac使其无效，以便对比**
            cS.pps_max_contrib_frac = 0.1; % **设置为0，模拟NoPPS的情况**
            cS.pps_tax_rate_withdrawal = 0.03;
            cS.pps_return_rate_premium = 0.0;
            cS.pps_withdrawal_rate = 0.15;
            cS.pps_contribution_age_max = cS.aR_idx_orig - 1;
            cS.pps_withdrawal_age_min = cS.aR_idx_orig;
            cS.pps_in_K = true;
            cS.pps_bequeathable = false;
        end

        function [Y, R, w_net, b_payg, MPL_gross_out] = HHPrices_Huggett(K_productive, L_total_eff, cS, paramS_in)
            % (与你上一轮修改后能输出合理替代率的版本相同)
            if K_productive <= 0, K_productive=1e-6; warning('K_prod non-positive in HHPrices, set to small.'); end
            if L_total_eff <= 0, L_total_eff=1e-6; warning('L_eff non-positive in HHPrices, set to small.'); end
            total_retirees_count_for_b = 0;
            if nargin < 4 || ~isfield(paramS_in, 'Z_ss_counts') || isempty(paramS_in.Z_ss_counts)
                warning('Z_ss_counts not provided to HHPrices_Huggett or is empty. Using static retireMass_new for b calculation, which might be inconsistent.');
                pop_normalizer = sum(cS.initial_pop(isfinite(cS.initial_pop))); if pop_normalizer == 0, pop_normalizer =100; else pop_normalizer = pop_normalizer / (sum(cS.initial_pop(isfinite(cS.initial_pop))>0) / cS.aD_new * 100) ; end
                total_retirees_count_for_b = cS.retireMass_new * pop_normalizer ;
                 if cS.retireMass_new < 1e-9 && sum(cS.initial_pop)>1e-9, total_retirees_count_for_b = 0; end
            else
                total_retirees_count_for_b = sum(paramS_in.Z_ss_counts(cS.aR_new+1:end));
            end
            Y   = cS.A * (K_productive^cS.alpha) * (L_total_eff^(1-cS.alpha));
            MPK_gross = cS.alpha * Y / K_productive; MPL_gross = (1-cS.alpha) * Y / L_total_eff; MPL_gross_out = MPL_gross;
            R   = 1 + (MPK_gross - cS.ddk); R   = max(1.0 + 1e-6, R);
            w_net = MPL_gross * (1 - cS.theta); w_net = max(0, w_net); % cS.theta是固定的PAYG贡献率
            b_payg = 0;
            if total_retirees_count_for_b > 1e-9 && L_total_eff > 0
                total_gross_wage_bill_payg = MPL_gross * L_total_eff; total_payg_revenue = cS.theta * total_gross_wage_bill_payg;
                b_payg = total_payg_revenue / total_retirees_count_for_b;
            end
            b_payg = max(0, b_payg);
        end

        function [logGridV, trProbM, prob1V] = EarningProcess_olgm(cS)
            % (与之前版本相同)
             [logGridV_raw, trProbM_calc] = main_olg_v2_utils.tauchen(cS.nw, cS.lePersistence, cS.leShockStd, 0, cS.leWidth);
             gridMin = logGridV_raw(1); gridMax = logGridV_raw(end); normGridMin = gridMin - 2*cS.leShockStd; normGridMax = gridMax + 2*cS.leShockStd; prob1V_calc = [];
              try [prob1V_calc, ~, ~] = main_olg_v2_utils.norm_grid(logGridV_raw, normGridMin, normGridMax, 0, cS.leSigma1); prob1V_calc = prob1V_calc(:);
              catch ME, warning('norm_grid failed: %s. Using uniform.', ME.message); prob1V_calc = ones(cS.nw, 1)/cS.nw; end
             logGridV_calc = logGridV_raw(:); logGridV_calc = logGridV_calc - logGridV_calc(1) - 1;
             if any(abs(sum(trProbM_calc, 2) - 1) > 1e-6), row_sums = sum(trProbM_calc, 2); row_sums(row_sums <= 1e-9) = 1; trProbM_calc = bsxfun(@rdivide, trProbM_calc, row_sums); end
             if abs(sum(prob1V_calc) - 1) > 1e-6, prob1V_calc = prob1V_calc ./ sum(prob1V_calc); end
             logGridV = logGridV_calc; trProbM = trProbM_calc; prob1V = prob1V_calc;
        end
        function [y, trProbM] = tauchen(N, pRho, pSigma, pMu, n_std) % (与之前版本相同)
             std_y = sqrt(pSigma^2 / (1 - pRho^2)); if abs(1-pRho)<1e-9, std_y = pSigma*100; end; y_max = n_std*std_y; y_min = -y_max; y_val = linspace(y_min, y_max, N); d = 0; if N > 1, d = y_val(2) - y_val(1); end; trProbM = zeros(N, N); for iR = 1:N, for iC = 1:N, mn = pRho*y_val(iR); if iC==1, trProbM(iR,iC)=normcdf((y_val(1)-mn+d/2)/pSigma); elseif iC==N, trProbM(iR,iC)=1-normcdf((y_val(N)-mn-d/2)/pSigma); else trProbM(iR,iC)=normcdf((y_val(iC)-mn+d/2)/pSigma)-normcdf((y_val(iC)-mn-d/2)/pSigma); end; end; end; rs=sum(trProbM,2); rs(rs<=1e-9)=1; trProbM=bsxfun(@rdivide,trProbM,rs); um=pMu/(1-pRho); if abs(1-pRho)<1e-9&&pMu~=0,um=sign(pMu)*inf; end; if isfinite(um), y_val=y_val+um; end; y = y_val;
        end
        function [massV, lbV, ubV] = norm_grid(xV, xMin, xMax, mu, sig) % (与之前版本相同)
             n=length(xV); xV=xV(:)'; if n>1&&any(xV(2:n)<xV(1:n-1)), error('xV non-decreasing'); end; lbV=[]; ubV=[]; if n>1, xM=0.5*(xV(1:n-1)+xV(2:n)); lbV=[xMin,xM]; ubV=[xM,xMax]; elseif n==1, lbV=xMin; ubV=xMax; else massV=[]; return; end; cdfV=normcdf([lbV,ubV(end)],mu,sig); massV=diff(cdfV); if any(massV < -1e-9), warning('Negative mass in norm_grid.'); massV(massV<0)=0; end; ts=sum(massV); if ts>1e-9, massV=massV/ts; else if n>0, massV=ones(1,n)/n; else massV=[]; end; end; massV=massV(:); lbV=lbV(:); ubV=ubV(:); if n>1&&any(ubV<lbV), error('ubV < lbV'); end
        end
        function eIdxM = MarkovChainSimulation(nSim, T_sim, prob0V, trProbM, rvInM) % (与之前版本相同)
             ns=length(prob0V); if size(trProbM,1)~=ns||size(trProbM,2)~=ns,error('TrMat size mismatch');end; if abs(sum(prob0V)-1)>1e-5,prob0V=prob0V./sum(prob0V);end; if any(abs(sum(trProbM,2)-1)>1e-5), row_sums=sum(trProbM,2); row_sums(row_sums<=1e-9)=1; trProbM=bsxfun(@rdivide,trProbM,row_sums); end; if size(rvInM,1)~=nSim||size(rvInM,2)~=T_sim,error('RV size mismatch');end; cP0=cumsum(prob0V(:)'); cPT=cumsum(trProbM,2); cPT(:,ns)=1.0; eIdxM=zeros(nSim,T_sim,'uint16'); eIdxM(:,1)=1+sum(bsxfun(@gt,rvInM(:,1),cP0),2); for t_mc=1:(T_sim-1), cSI=eIdxM(:,t_mc); vPI=(cSI>=1)&(cSI<=ns); if ~all(vPI),cSI(~vPI)=1;eIdxM(:,t_mc)=cSI;end; cPt=cPT(cSI,:); eIdxM(:,t_mc+1)=1+sum(bsxfun(@gt,rvInM(:,t_mc+1),cPt),2); end
        end
        function [HHlaborM_group, L] = LaborSupply_Huggett(eIdxM, cS, paramS, Z_ss_norm_group) % (与之前版本相同)
            nSim = size(eIdxM, 1); HHlaborM_group = zeros(nSim, cS.aD_new); ageToGroupMap = zeros(cS.aD_orig, 1); for a_new_map = 1:cS.aD_new, ageToGroupMap(cS.physAgeMap{a_new_map}) = a_new_map; end; leGridV_col = paramS.leGridV(:); HHlaborM_annual_temp = zeros(nSim, cS.aD_orig); for a_orig = 1 : cS.aD_orig, if a_orig < cS.aR_idx_orig, a_new = ageToGroupMap(a_orig); HHlaborM_annual_temp(:, a_orig) = cS.ageEffV_new(a_new) .* leGridV_col(eIdxM(:, a_orig)); end; end; for a_new = 1:cS.aD_new, annual_indices = cS.physAgeMap{a_new}; if ~isempty(annual_indices), HHlaborM_group(:, a_new) = mean(HHlaborM_annual_temp(:, annual_indices), 2); end; end; L_working_groups = mean(HHlaborM_group(:, 1:cS.aR_new), 1) * Z_ss_norm_group(1:cS.aR_new); L = max(0, L_working_groups);
        end
        function eIdxM = LaborEndowSimulation_olgm(cS, paramS) % (与之前版本相同)
            rng(433); rvInM=rand([cS.nSim,cS.aD_orig]);
            eIdxM = main_olg_v2_utils.MarkovChainSimulation(cS.nSim, cS.aD_orig, paramS.leProb1V, paramS.leTrProbM, rvInM);
        end

        function [kHistM, kPpsHistM, cHistM] = HHSimulation_olgm(kPolM, cPpsPolM, cPolM_consump, eIdxM, R_iter, w_iter, T_iter, bV_iter, paramS, cS)
            % 模拟家庭年度路径 (非PPS资产k, PPS资产k_pps, 消费c)
            % (与你之前v2_pps版本中的 HHSimulation_olgm 函数逻辑相同, 但插值方法改为'spline')
            validateattributes(kPolM, {'double'}, {'finite', 'nonnan', 'nonempty', 'real', 'size', [cS.nk, cS.nw, cS.aD_new]});
            validateattributes(cPpsPolM, {'double'}, {'finite', 'nonnan', 'nonempty', 'real','>=',0, 'size', [cS.nk, cS.nw, cS.aD_new]});
            validateattributes(cPolM_consump, {'double'}, {'finite', 'nonnan', 'nonempty', 'real', 'size', [cS.nk, cS.nw, cS.aD_new]});
            validateattributes(eIdxM, {'uint16', 'double'}, {'nonempty', 'integer', 'positive', '<=', cS.nw, 'size', [cS.nSim, cS.aD_orig]});

            ageToGroupMap = zeros(cS.aD_orig, 1);
            for a_new_idx = 1:cS.aD_new, annual_indices_map = cS.physAgeMap{a_new_idx}; if ~isempty(annual_indices_map), ageToGroupMap(annual_indices_map) = a_new_idx; end; end

            nSim_sim   = size(eIdxM, 1);
            kHistM_sim    = zeros(nSim_sim, cS.aD_orig + 1);
            kPpsHistM_sim = zeros(nSim_sim, cS.aD_orig + 1);
            cHistM_sim    = zeros(nSim_sim, cS.aD_orig);
            kGridV_col_sim = cS.kGridV(:);
            leGridV_col_sim = paramS.leGridV(:);

            kPolInterp_sim = cell(cS.nw, cS.aD_new);
            cPpsPolInterp_sim = cell(cS.nw, cS.aD_new);
            cPolConsumpInterp_sim = cell(cS.nw, cS.aD_new);
             for a_new_sim_interp=1:cS.aD_new
                 for ie_sim_interp=1:cS.nw
                      kPolInterp_sim{ie_sim_interp, a_new_sim_interp}    = griddedInterpolant(kGridV_col_sim, kPolM(:, ie_sim_interp, a_new_sim_interp),  'linear', 'linear');
                      cPpsPolInterp_sim{ie_sim_interp, a_new_sim_interp} = griddedInterpolant(kGridV_col_sim, cPpsPolM(:, ie_sim_interp, a_new_sim_interp),  'linear', 'linear');
                      cPolConsumpInterp_sim{ie_sim_interp, a_new_sim_interp} = griddedInterpolant(kGridV_col_sim, cPolM_consump(:, ie_sim_interp, a_new_sim_interp),  'linear', 'linear');
                 end
             end

            pps_net_annual_return_sim = (R_iter - 1) + cS.pps_return_rate_premium;

            for a_orig_sim_loop = 1 : cS.aD_orig
                a_new_group_idx_sim_loop = ageToGroupMap(a_orig_sim_loop);
                kNowV_sim_loop = kHistM_sim(:, a_orig_sim_loop);
                kPpsNowV_sim_loop = kPpsHistM_sim(:, a_orig_sim_loop);
                kNextNonPpsV_sim_loop = zeros(nSim_sim, 1);
                cPpsContribV_sim_loop = zeros(nSim_sim, 1);
                cConsumpValV_sim_loop = zeros(nSim_sim, 1);
                kPpsNextV_sim_loop    = zeros(nSim_sim, 1);
                is_working_sim_loop = (a_orig_sim_loop < cS.aR_idx_orig);
                pps_withdrawal_pretax_sim_loop = zeros(nSim_sim,1);
                if ~is_working_sim_loop && cS.pps_active && a_orig_sim_loop >= cS.pps_withdrawal_age_min
                    pps_withdrawal_pretax_sim_loop = kPpsNowV_sim_loop * cS.pps_withdrawal_rate;
                end
                for ie_sim_loop = 1 : cS.nw
                    simIdxV_loop = find(eIdxM(:, a_orig_sim_loop) == ie_sim_loop);
                    if ~isempty(simIdxV_loop)
                         kNow_ie_clamped_sim_loop = max(kGridV_col_sim(1), min(kGridV_col_sim(end), kNowV_sim_loop(simIdxV_loop)));
                         kNextNonPpsV_sim_loop(simIdxV_loop) = kPolInterp_sim{ie_sim_loop, a_new_group_idx_sim_loop}(kNow_ie_clamped_sim_loop);
                         cPpsContribV_sim_loop(simIdxV_loop) = cPpsPolInterp_sim{ie_sim_loop, a_new_group_idx_sim_loop}(kNow_ie_clamped_sim_loop);
                         cConsumpValV_sim_loop(simIdxV_loop) = cPolConsumpInterp_sim{ie_sim_loop, a_new_group_idx_sim_loop}(kNow_ie_clamped_sim_loop);
                         if is_working_sim_loop && cS.pps_active && a_orig_sim_loop <= cS.pps_contribution_age_max
                            cPpsContribV_sim_loop(simIdxV_loop) = max(0, cPpsContribV_sim_loop(simIdxV_loop));
                         else, cPpsContribV_sim_loop(simIdxV_loop) = 0; end
                         if cS.pps_active
                             kPpsNextV_sim_loop(simIdxV_loop) = (kPpsNowV_sim_loop(simIdxV_loop) + cPpsContribV_sim_loop(simIdxV_loop) - pps_withdrawal_pretax_sim_loop(simIdxV_loop)) * (1 + pps_net_annual_return_sim);
                             kPpsNextV_sim_loop(simIdxV_loop) = max(0, kPpsNextV_sim_loop(simIdxV_loop));
                         else, kPpsNextV_sim_loop(simIdxV_loop) = kPpsNowV_sim_loop(simIdxV_loop); end
                    end
                end
                 kHistM_sim(:, a_orig_sim_loop + 1)    = max(cS.kMin, min(cS.kMax, kNextNonPpsV_sim_loop));
                 kPpsHistM_sim(:, a_orig_sim_loop + 1) = max(0, kPpsNextV_sim_loop);
                 cHistM_sim(:, a_orig_sim_loop)        = max(cS.cFloor, cConsumpValV_sim_loop);
                 kHistM_sim(~isfinite(kHistM_sim(:, a_orig_sim_loop + 1)), a_orig_sim_loop + 1) = cS.kMin;
                 kPpsHistM_sim(~isfinite(kPpsHistM_sim(:, a_orig_sim_loop + 1)), a_orig_sim_loop + 1) = 0;
                 cHistM_sim(~isfinite(cHistM_sim(:, a_orig_sim_loop)), a_orig_sim_loop) = cS.cFloor;
            end
            kHistM = kHistM_sim(:, 1:cS.aD_orig); kPpsHistM = kPpsHistM_sim(:, 1:cS.aD_orig); cHistM = cHistM_sim;
            validateattributes(kHistM, {'double'}, {'finite', 'nonnan', 'nonempty', 'real', '>=', cS.kMin - 1e-6, '<=', cS.kMax + 1e-6, 'size', [nSim_sim, cS.aD_orig]});
            validateattributes(kPpsHistM, {'double'}, {'finite', 'nonnan', 'nonempty', 'real', '>=', -1e-6, 'size', [nSim_sim, cS.aD_orig]});
            validateattributes(cHistM, {'double'}, {'finite', 'nonnan', 'nonempty', 'real', '>=', cS.cFloor - 1e-6, 'size', [nSim_sim, cS.aD_orig]});
        end

        function [muM, utilM] = CES_utility(cM, sig, cS)
            % (与之前版本相同)
             if ~isscalar(sig)||sig<=0,error('Sigma must be positive scalar.');end; min_c = cS.cFloor; vld = (cM >= min_c); c_c = max(min_c, cM); utilM = -Inf(size(cM)); muM = Inf(size(cM)); if abs(sig-1)<1e-6, utilM(vld) = log(c_c(vld)); muM(vld) = 1 ./ c_c(vld); else, utilM(vld) = (c_c(vld).^(1-sig)) ./ (1-sig); muM(vld) = c_c(vld).^(-sig); end; utilM(~vld) = -1e10 - (min_c - cM(~vld))*1e10; muM(~vld) = c_c(~vld).^(-sig) + 1e10;
        end

        function incomeM = HHIncome_Huggett(k_non_pps_grid, R_price, w_price, T_price, b_price_group, a_new_idx, paramS_income, cS)
            % (与之前版本相同)
            nk_income = length(k_non_pps_grid); leGridV_col_income = paramS_income.leGridV(:); nonCapIncomeNetM_income = zeros(1, cS.nw);
            if a_new_idx <= cS.aR_new, net_wage_income_vec_income = w_price * cS.ageEffV_new(a_new_idx) .* leGridV_col_income; nonCapIncomeNetM_income(1,:) = net_wage_income_vec_income' + T_price + b_price_group;
            else, nonCapIncomeNetM_income(1,:) = b_price_group + T_price; end
            incomeM = R_price * k_non_pps_grid(:) * ones(1, cS.nw) + ones(nk_income, 1) * nonCapIncomeNetM_income;
            if any(~isfinite(incomeM(:))), warning('Non-finite income in HHIncome_Huggett group %d.',a_new_idx); incomeM(~isfinite(incomeM)) = 1e-6; end
        end

        function [cPolM, kPolM, cPpsPolM, valueM] = HHSolution_VFI_Huggett(R_vfi, w_vfi, T_vfi, bV_vfi, paramS_vfi, cS)
            % (与之前版本相同)
            cPolM  = zeros(cS.nk, cS.nw, cS.aD_new); kPolM  = zeros(cS.nk, cS.nw, cS.aD_new); cPpsPolM = zeros(cS.nk, cS.nw, cS.aD_new); valueM = zeros(cS.nk, cS.nw, cS.aD_new);
            for a_vfi_loop_main = cS.aD_new : -1 : 1
               vPrimeM_vfi = []; if a_vfi_loop_main < cS.aD_new, vPrimeM_vfi = valueM(:,:,a_vfi_loop_main+1); end
               [cPolM(:,:,a_vfi_loop_main), kPolM(:,:,a_vfi_loop_main), cPpsPolM(:,:,a_vfi_loop_main), valueM(:,:,a_vfi_loop_main)] = ...
                  main_olg_v2_utils.HHSolutionByAge_VFI_Huggett(a_vfi_loop_main, vPrimeM_vfi, R_vfi, w_vfi, T_vfi, bV_vfi(a_vfi_loop_main), paramS_vfi, cS);
            end
            validateattributes(cPolM, {'double'}, {'finite', 'nonnan', 'real','>=',cS.cFloor-1e-6,'size',[cS.nk,cS.nw,cS.aD_new]});
            validateattributes(kPolM, {'double'}, {'finite', 'nonnan', 'real','>=',cS.kMin-1e-6,'<=',cS.kMax+1e-6,'size',[cS.nk,cS.nw,cS.aD_new]});
            validateattributes(cPpsPolM, {'double'}, {'finite', 'nonnan', 'real','>=',0,'size',[cS.nk,cS.nw,cS.aD_new]});
            validateattributes(valueM, {'double'}, {'finite', 'nonnan', 'real','size',[cS.nk,cS.nw,cS.aD_new]});
        end

                function [cPolM_age, kPolM_age, cPpsPolM_age_actual, valueM_age_out] = HHSolutionByAge_VFI_Huggett(...
                a_idx, vPrime_keM_age, R_age, w_age, T_age, b_age_group, paramS_age, cS)
            % 分年龄组的VFI求解
            % 如果PPS不活跃或缴费分数为0，则cPpsPolM_age_actual将为0

            if a_idx < cS.aD_new && isempty(vPrime_keM_age), error('VFI by Age: vPrime empty a<aD'); end
            if a_idx == cS.aD_new && ~isempty(vPrime_keM_age), error('VFI by Age: vPrime not empty a=aD'); end
            if a_idx < cS.aD_new, validateattributes(vPrime_keM_age,{'double'},{'finite','nonnan','real','size',[cS.nk,cS.nw]}); end

            budgetM_before_any_saving_choice = main_olg_v2_utils.HHIncome_Huggett(cS.kGridV, R_age, w_age, T_age, b_age_group, a_idx, paramS_age, cS);

            fminOpt_age=optimset('fminbnd'); fminOpt_age.TolX=1e-6; fminOpt_age.Display='none';

            cPolM_age = zeros(cS.nk, cS.nw);
            kPolM_age = zeros(cS.nk, cS.nw); % 代表非PPS储蓄 k'
            cPpsPolM_age_actual = zeros(cS.nk, cS.nw); % 实际PPS缴费，强制初始化为0
            valueM_age_out = -Inf(cS.nk, cS.nw);

            is_pps_contrib_age_val = (a_idx <= cS.aR_new);

            if a_idx == cS.aD_new % 最后一个年龄组
               % cPpsPolM_age_actual 保持为0
               cPolM_age = max(cS.cFloor, budgetM_before_any_saving_choice);
               kPolM_age = cS.kMin*ones(cS.nk,cS.nw);
               [~, val_temp_age] = main_olg_v2_utils.CES_utility(cPolM_age, cS.sigma, cS);
               valueM_age_out = val_temp_age;
               valueM_age_out(~isfinite(valueM_age_out)) = -1e12;
            else
               ExValuePrime_k_given_e_age = zeros(cS.nk, cS.nw);
                for ikPrime_age = 1:cS.nk
                   ExValuePrime_k_given_e_age(ikPrime_age, :) = vPrime_keM_age(ikPrime_age, :) * paramS_age.leTrProbM';
                end
                ExValuePrime_k_given_e_age(~isfinite(ExValuePrime_k_given_e_age))=-1e10;

               vPInt_age=cell(1,cS.nw);
               for ie_age_interp=1:cS.nw
                   % 使用与能收敛的简化模型相同的插值方法
                   vPInt_age{ie_age_interp} = griddedInterpolant(cS.kGridV, ExValuePrime_k_given_e_age(:, ie_age_interp), 'linear','linear');
               end

               parfor ik_age_loop = 1:cS.nk
                   c_row_age=zeros(1,cS.nw); 
                   k_row_age=zeros(1,cS.nw); 
                   cpps_row_for_this_k_state=zeros(1,cS.nw); 
                   v_row_age=zeros(1,cS.nw);

                   for ie_age_loop_inner = 1 : cS.nw
                        vPofK_interpolant_age = vPInt_age{ie_age_loop_inner};
                        budget_ike_initial = budgetM_before_any_saving_choice(ik_age_loop,ie_age_loop_inner);
                        
                        cPps_chosen_this_state_val = 0.0; % 强制为0

                        if is_pps_contrib_age_val && cS.pps_active && cS.pps_max_contrib_frac > 1e-9 
                           wage_income_pps_limit_age = w_age * cS.ageEffV_new(a_idx) * paramS_age.leGridV(ie_age_loop_inner);
                           max_contrib_pps_age = wage_income_pps_limit_age * cS.pps_max_contrib_frac;
                           
                           if budget_ike_initial - max_contrib_pps_age >= cS.cFloor + cS.kMin 
                               cPps_chosen_this_state_val = max_contrib_pps_age;
                           else
                               cPps_chosen_this_state_val = min(max_contrib_pps_age, max(0, budget_ike_initial - cS.cFloor - cS.kMin ));
                           end
                           cPps_chosen_this_state_val = max(0, cPps_chosen_this_state_val);
                        end
                        cpps_row_for_this_k_state(ie_age_loop_inner) = cPps_chosen_this_state_val;
                        
                        budget_for_k_prime_and_c = budget_ike_initial - cPps_chosen_this_state_val;

                        [c_row_age(ie_age_loop_inner), k_row_age(ie_age_loop_inner), v_row_age(ie_age_loop_inner)] = ...
                              main_olg_v2_utils.HHSolutionByOneState_OptK(a_idx, budget_for_k_prime_and_c, vPofK_interpolant_age, fminOpt_age, cS, paramS_age);
                   end
                   cPolM_age(ik_age_loop,:) = c_row_age;
                   kPolM_age(ik_age_loop,:) = k_row_age; 
                   cPpsPolM_age_actual(ik_age_loop,:) = cpps_row_for_this_k_state; 
                   valueM_age_out(ik_age_loop,:) = v_row_age;
               end
            end
        end

        function [c, kPrime, ValueFunc] = HHSolutionByOneState_OptK(a_new_osa_k, budget_remaining_osa_k, vPofK_int_osa_k, fminOpt_osa_k, cS, paramS_osa_k)
            % 这个函数应该与原始简化模型utils.m中的单状态优化函数HHSolutionByOneState_VFI_Huggett的逻辑完全一致
            % 只是调用的辅助函数是本类中的
            kPMin_osa_k = cS.kMin;
            kPMax_osa_k = budget_remaining_osa_k - cS.cFloor;
            ValueFunc = -Inf;
            if kPMax_osa_k <= kPMin_osa_k
                kPrime = kPMin_osa_k; c = max(cS.cFloor, budget_remaining_osa_k - kPrime);
                 [~, u_c_osa_k] = main_olg_v2_utils.CES_utility(c, cS.sigma, cS); % 修改：获取第一个输出来匹配原始函数
                if a_new_osa_k < cS.aD_new
                    EV_c_osa_k = -Inf; try kP_eval_osa_k = max(cS.kGridV(1),min(cS.kGridV(end),kPrime)); EV_c_osa_k=vPofK_int_osa_k(kP_eval_osa_k); catch ME_osa_k; warning('OptK corner interp error age %d: %s',a_new_osa_k,ME_osa_k.message); if kPrime<cS.kGridV(1),EV_c_osa_k=vPofK_int_osa_k(cS.kGridV(1));else,EV_c_osa_k=vPofK_int_osa_k(cS.kGridV(end));end;end; if ~isfinite(EV_c_osa_k),EV_c_osa_k=-1e12;end
                    s_transition_osa_k = cS.s_1yr_transitionV(a_new_osa_k); ValueFunc = u_c_osa_k + cS.beta*s_transition_osa_k*EV_c_osa_k;
                else, ValueFunc = u_c_osa_k; end
                if ~isfinite(ValueFunc), ValueFunc = -1e12; end
            else
                kPMin_opt_osa_k = max(cS.kMin, kPMin_osa_k); kPMax_opt_osa_k = max(kPMin_opt_osa_k + 1e-9, min(cS.kMax, kPMax_osa_k));
                if kPMin_opt_osa_k >= kPMax_opt_osa_k
                    kPrime_opt_osa_k = kPMin_opt_osa_k; 
                    [negV_osa_k, ~] = main_olg_v2_utils.BellmanInner_NoPPS_Equivalent(kPrime_opt_osa_k, a_new_osa_k, budget_remaining_osa_k, vPofK_int_osa_k, cS, paramS_osa_k); % 调用一个等价的NoPPS BellmanInner
                    ValueFunc = -negV_osa_k; 
                    warning('OptK bounds invalid age %d',a_new_osa_k);
                else
                    obj_osa_k = @(kP_osa_k) main_olg_v2_utils.negBellmanObjective_NoPPS_Equivalent(kP_osa_k, a_new_osa_k, budget_remaining_osa_k, vPofK_int_osa_k, cS, paramS_osa_k); % 调用等价的NoPPS目标函数
                    [kPrime_opt_osa_k, negV_osa_k, eflag_osa_k] = fminbnd(obj_osa_k, kPMin_opt_osa_k, kPMax_opt_osa_k, fminOpt_osa_k);
                    if eflag_osa_k <= 0 || abs(kPrime_opt_osa_k-kPMin_opt_osa_k)<1e-7 || abs(kPrime_opt_osa_k-kPMax_opt_osa_k)<1e-7
                        [nV_min_osa_k,~]=main_olg_v2_utils.BellmanInner_NoPPS_Equivalent(kPMin_opt_osa_k,a_new_osa_k,budget_remaining_osa_k,vPofK_int_osa_k,cS,paramS_osa_k); 
                        [nV_max_osa_k,~]=main_olg_v2_utils.BellmanInner_NoPPS_Equivalent(kPMax_opt_osa_k,a_new_osa_k,budget_remaining_osa_k,vPofK_int_osa_k,cS,paramS_osa_k);
                        if nV_min_osa_k <= nV_max_osa_k + 1e-9, kPrime_opt_osa_k=kPMin_opt_osa_k; negV_osa_k=nV_min_osa_k; else, kPrime_opt_osa_k=kPMax_opt_osa_k; negV_osa_k=nV_max_osa_k; end
                        if eflag_osa_k <= 0, warning('fminbnd issue OptK age %d, eflag %d',a_new_osa_k,eflag_osa_k); end
                    end
                    ValueFunc = -negV_osa_k;
                end
                kPrime = kPrime_opt_osa_k; c = max(cS.cFloor, budget_remaining_osa_k - kPrime);
                if ~isfinite(ValueFunc), ValueFunc = -1e12; end
            end
            kPrime = max(cS.kMin, min(cS.kMax, kPrime)); c = max(cS.cFloor, budget_remaining_osa_k - kPrime);
            [~, u_final_c_osa_k] = main_olg_v2_utils.CES_utility(c, cS.sigma, cS); % 修改：获取第一个输出来匹配原始函数
            if a_new_osa_k < cS.aD_new
                EV_A_final_osa_k = -Inf; try kP_eval_final_osa_k = max(cS.kGridV(1),min(cS.kGridV(end),kPrime)); EV_A_final_osa_k=vPofK_int_osa_k(kP_eval_final_osa_k); catch; if kPrime<cS.kGridV(1),EV_A_final_osa_k=vPofK_int_osa_k(cS.kGridV(1));else,EV_A_final_osa_k=vPofK_int_osa_k(cS.kGridV(end));end;end; if ~isfinite(EV_A_final_osa_k),EV_A_final_osa_k=-1e12;end
                s_transition_final_osa_k = cS.s_1yr_transitionV(a_new_osa_k); ValueFunc = u_final_c_osa_k + cS.beta*s_transition_final_osa_k*EV_A_final_osa_k;
            else, ValueFunc = u_final_c_osa_k; end
            if ~isfinite(kPrime), kPrime = cS.kMin; end; if ~isfinite(c), c = cS.cFloor; end; if ~isfinite(ValueFunc), ValueFunc = -1e12; end
            validateattributes(kPrime,{'double'},{'finite','nonnan','scalar','real','>=',cS.kMin-1e-6,'<=',cS.kMax+1e-6});
            validateattributes(c,{'double'},{'finite','nonnan','scalar','real','>=',cS.cFloor-1e-6});
            validateattributes(ValueFunc,{'double'},{'finite','nonnan','scalar','real'});
        end

        % --- Nested Helper Functions for OptK (NoPPS Equivalent) ---
        function negV = negBellmanObjective_NoPPS_Equivalent(kP, ag_obj, budget_obj, evInt_obj, cS_in_obj, pS_in_obj)
            [negV, ~] = main_olg_v2_utils.BellmanInner_NoPPS_Equivalent(kP, ag_obj, budget_obj, evInt_obj, cS_in_obj, pS_in_obj);
        end

        function [negVal, Val] = BellmanInner_NoPPS_Equivalent(kP_inner, ag_inner, budget_inner, evInt_inner, cS_in, pS_in)
            % 这个函数的逻辑应该和原始简化模型utils.m中的BellmanInner完全一样
            cons = max(cS_in.cFloor, budget_inner - kP_inner);
            [~, util_val] = main_olg_v2_utils.CES_utility(cons, cS_in.sigma, cS_in); % 修改：获取第一个输出来匹配原始函数
            if ~isfinite(util_val), util_val = -1e12; end; Val = util_val;
            if ag_inner < cS_in.aD_new
                evF_val = -Inf;
                try kP_eval_val = max(cS_in.kGridV(1),min(cS_in.kGridV(end),kP_inner)); evF_val=evInt_inner(kP_eval_val);
                catch ME_bellman_no_pps; warning('BellmanInner_NoPPS_Equivalent interp error age %d: %s', ag_inner, ME_bellman_no_pps.message); if kP_inner<cS_in.kGridV(1),evF_val=evInt_inner(cS_in.kGridV(1));else,evF_val=evInt_inner(cS_in.kGridV(end));end; end
                if ~isfinite(evF_val),evF_val=-1e12;end
                s_transition_val = cS_in.s_1yr_transitionV(ag_inner); Val = util_val + cS_in.beta*s_transition_val*evF_val;
            end
            if ~isfinite(Val), negVal=1e12; Val=-1e12; else negVal = -Val; end
        end

    end % End methods (Static)
end % End classdef