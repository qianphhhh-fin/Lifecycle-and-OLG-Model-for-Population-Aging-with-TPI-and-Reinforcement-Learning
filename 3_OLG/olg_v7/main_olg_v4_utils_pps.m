% --- START OF FILE main_olg_v4_utils_pps.m ---

classdef main_olg_v4_utils_pps
    % OLG 模型 v4 的工具函数
    % 修改版：采用“混合时间单位”，包含个人养老金计划 (PPS),
    % PAYG工资税率外生 (theta), PAYG福利内生 (b_payg),
    % 并整合 Heer (2020) 风格的政府财政 (tau_k, tau_c, tau_l, G, B)

    methods (Static)

        % =====================================================================
        % == 人口动态函数 (年龄组层面 - 与V3一致) =========
        % =====================================================================

        function popS = initPopulation(cS) % cS: 参数结构体
            popS.Z = zeros(cS.aD_new, 1);
            initial_total = sum(cS.initial_pop);
            if initial_total > 0 && length(cS.initial_pop) == cS.aD_new
                 popS.Z(:, 1) = cS.initial_pop(:) / initial_total * 100;
            else
                warning('初始人口不匹配或总和为零。设置为均匀的初始年龄组人口分布。');
                 popS.Z(:, 1) = 100 / cS.aD_new;
            end
            popS.totalPop = sum(popS.Z(:, 1));
            if popS.totalPop > 1e-9
                 popS.ageDist = popS.Z(:, 1) / popS.totalPop;
            else
                 popS.ageDist = zeros(cS.aD_new, 1);
            end
            popS.initialAgeDist = popS.ageDist;
            fprintf('初始年龄组人口已设置。总人口=%.2f\n', popS.totalPop);
        end

        function popS = populationDynamics(popS, cS)
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
                    survival_prob = 0;
                    if (a-1) >= 1 && (a-1) <= length(cS.survivalProbV_popdyn)
                         survival_prob = cS.survivalProbV_popdyn(a-1);
                    end
                    Z_next_period(a) = Z_current_period(a-1) * survival_prob;
                    Z_next_period(a) = max(0, Z_next_period(a));
                end

                Z_history(:, t+1) = Z_next_period;
                totalPop_history(t+1) = sum(Z_next_period);
                if totalPop_history(t+1) > 1e-9
                    ageDist_history(:, t+1) = Z_next_period / totalPop_history(t+1);
                else
                    ageDist_history(:, t+1) = 0; totalPop_history(t+1) = 0;
                end

                current_check_period = t + 1;
                if current_check_period >= cS.bgp_window + 1
                    stable = true;
                    for w_idx = 1:cS.bgp_window
                       hist_idx1 = current_check_period - w_idx + 1;
                       hist_idx2 = current_check_period - w_idx;
                       if hist_idx1 > 0 && hist_idx2 > 0 && hist_idx1 <= size(ageDist_history,2) && hist_idx2 <= size(ageDist_history,2)
                           change = norm(ageDist_history(:, hist_idx1) - ageDist_history(:, hist_idx2));
                           if change >= cS.bgp_tolerance
                               stable = false; break;
                           end
                       else
                           stable = false; break;
                       end
                    end
                    if stable
                        fprintf('\n人口稳态 (年龄组) 在第 %d 期达到。\n', t);
                        bgp_reached_flag = true; actual_periods_run = t; break;
                    end
                end
            end

            final_period_idx = min(actual_periods_run + 1, size(Z_history,2));
            popS.Z = Z_history(:, 1:final_period_idx);
            popS.totalPop = totalPop_history(1:final_period_idx);
            popS.ageDist = ageDist_history(:, 1:final_period_idx);

            depRatio_history = zeros(1, actual_periods_run);
             for th = 1:actual_periods_run
                 Z_t = Z_history(:, th);
                 working_pop = sum(Z_t(1:cS.aR_new));
                 retired_pop = sum(Z_t(cS.aR_new+1:end));
                 if working_pop > 1e-9
                     depRatio_history(th) = retired_pop / working_pop;
                 else
                     depRatio_history(th) = inf;
                 end
             end
            popS.dependencyRatio = depRatio_history;

            fprintf('人口动态模拟完成。运行期数: %d。达到BGP: %d\n', actual_periods_run, bgp_reached_flag);
            if ~bgp_reached_flag
                fprintf('警告: 人口稳态 (年龄组) 未在 %d 期内达到。\n', max_periods_sim);
            end
        end

        function [Z_ss, dependency_ratio_ss, bgp_reached, bgp_period] = detectSteadyStatePopulation(popS, cS)
             actual_periods_in_data = size(popS.Z, 2);
             bgp_reached = false;
             bgp_period = actual_periods_in_data - 1;

             if actual_periods_in_data < cS.bgp_window + 1
                 fprintf('人口模拟期数过短 (%d 数据点)，无法进行稳态检查 (窗口期 = %d)。\n', actual_periods_in_data, cS.bgp_window);
             else
                 fprintf('检查人口稳态 (年龄组, 最近 %d 期)...\n', cS.bgp_window);
                 for t_check_end = actual_periods_in_data : -1 : cS.bgp_window + 1
                     stable = true;
                     for w_idx = 0 : (cS.bgp_window - 1)
                         idx1 = t_check_end - w_idx;
                         idx2 = t_check_end - w_idx - 1;
                         if idx1 > 0 && idx2 > 0 && idx1 <= size(popS.ageDist,2) && idx2 <= size(popS.ageDist,2)
                             change = norm(popS.ageDist(:, idx1) - popS.ageDist(:, idx2));
                             if change >= cS.bgp_tolerance
                                 stable = false; break;
                             end
                         else
                             stable = false; break;
                         end
                     end
                     if stable
                         bgp_reached = true; bgp_period = t_check_end - 1;
                         fprintf('人口稳态 (年龄组) 从第 %d 期 (数据索引 %d) 开始检测到。\n', bgp_period, t_check_end);
                         break;
                     end
                 end
                  if ~bgp_reached
                      fprintf('未检测到人口稳态 (年龄组)。使用最终期数据。\n');
                      bgp_period = actual_periods_in_data - 1; % Ensure bgp_period is last period if not reached
                  end
             end

             ss_data_index = min(bgp_period + 1, size(popS.Z, 2));
             Z_ss = popS.Z(:, ss_data_index);
             Z_ss_norm = zeros(cS.aD_new, 1);
             if sum(Z_ss) > 1e-9
                 Z_ss_norm = Z_ss / sum(Z_ss);
             end

             valid_dep_ratio_index = min(bgp_period, length(popS.dependencyRatio));
            if isfield(popS, 'dependencyRatio') && ~isempty(popS.dependencyRatio) && valid_dep_ratio_index > 0 && valid_dep_ratio_index <= length(popS.dependencyRatio)
                 dependency_ratio_ss = popS.dependencyRatio(valid_dep_ratio_index);
            else
                  working_pop_ss = sum(Z_ss(1:cS.aR_new));
                  retired_pop_ss = sum(Z_ss(cS.aR_new+1:end));
                  if working_pop_ss > 1e-9
                      dependency_ratio_ss = retired_pop_ss / working_pop_ss;
                  else
                      dependency_ratio_ss = inf;
                  end
                  if (~isfield(popS, 'dependencyRatio') || isempty(popS.dependencyRatio)) && bgp_period > 0
                      warning('抚养比历史未找到或过短，已重新计算。');
                  end
             end

             figure('Name', 'V4: 初始 vs 稳态 年龄组人口分布');
             hold on;
             group_indices = 1:cS.aD_new;
             if isfield(popS, 'initialAgeDist') && ~isempty(popS.initialAgeDist)
                 bar(group_indices - 0.2, popS.initialAgeDist * 100, 0.4, 'b', 'DisplayName', '初始年龄组分布');
             else
                 warning('未找到或空的初始年龄分布用于绘图。');
             end
             bar(group_indices + 0.2, Z_ss_norm * 100, 0.4, 'r', 'DisplayName', sprintf('稳态年龄组分布 (第 %d 期)', bgp_period));
             hold off;
             xlabel(sprintf('年龄组索引 (1 至 %d)', cS.aD_new)); ylabel('占总人口百分比 (%)');
             title('V4: 初始 vs 稳态/最终 年龄组人口分布');
             legend('Location', 'best'); xticks(group_indices); grid on; drawnow;
             fprintf('已绘制初始与稳态/最终年龄组人口分布图。\n');
        end


        % =====================================================================
        % == 经济函数 (修改版：混合时间单位 & PPS & 外生PAYG税率 & Heer GBC) =======
        % =====================================================================

        function cS = ParameterValues_HuggettStyle()
            % 参数设定: 外生 PAYG 税率 theta, PPS, Heer GBC

            % --- 人口统计与分组参数 (与之前v3版本一致) ---
            cS.age1_orig = 20; cS.ageLast_orig = 98; cS.ageRetire_orig = 65; cS.popGrowth_orig = 0.012; % Annual pop growth
            cS.aD_orig = cS.ageLast_orig - cS.age1_orig + 1; cS.aR_idx_orig = cS.ageRetire_orig - cS.age1_orig + 1;
            cS.aW_orig = cS.aR_idx_orig - 1; cS.physAgeV_orig = (cS.age1_orig : cS.ageLast_orig)';
            cS.d_orig = [0.00159,0.00169,0.00174,0.00172,0.00165,0.00156,0.00149,0.00145,0.00145,0.00149,0.00156,0.00163,0.00171,0.00181,0.00193,0.00207,0.00225,0.00246,0.00270,0.00299,0.00332,0.00368,0.00409,0.00454,0.00504,0.00558,0.00617,0.00686,0.00766,0.00865,0.00955,0.01058,0.01162,0.01264,0.01368,0.01475,0.01593,0.01730,0.01891,0.02074,0.02271,0.02476,0.02690,0.02912,0.03143,0.03389,0.03652,0.03930,0.04225,0.04538,0.04871,0.05230,0.05623,0.06060,0.06542,0.07066,0.07636,0.08271,0.08986,0.09788,0.10732,0.11799,0.12895,0.13920,0.14861,0.16039,0.17303,0.18665,0.20194,0.21877,0.23601,0.25289,0.26973,0.28612,0.30128,0.31416,0.32915,0.34450,0.36018];
            if length(cS.d_orig) ~= cS.aD_orig, error('d_orig 长度不匹配'); end
            cS.s_orig = 1 - cS.d_orig;
            cS.yearStep = 5; cS.aD_new = ceil(cS.aD_orig / cS.yearStep); cS.aR_new = ceil(cS.aW_orig / cS.yearStep);
            cS.physAgeMap = cell(cS.aD_new, 1); for a = 1:cS.aD_new, startIdx = (a-1)*cS.yearStep + 1; endIdx = min(a*cS.yearStep, cS.aD_orig); cS.physAgeMap{a} = startIdx:endIdx; end
            cS.physAgeV_new = zeros(cS.aD_new, 1); for a = 1:cS.aD_new, cS.physAgeV_new(a) = cS.physAgeV_orig(cS.physAgeMap{a}(1)); end
            cS.s_1yr_transitionV = zeros(cS.aD_new, 1); for a = 1:(cS.aD_new - 1), lastYearIdx = cS.physAgeMap{a}(end); if lastYearIdx < cS.aD_orig, cS.s_1yr_transitionV(a) = cS.s_orig(lastYearIdx); else, cS.s_1yr_transitionV(a) = 0; end; end; cS.s_1yr_transitionV(cS.aD_new) = 0;

            %% 人口动态参数
            cS.initial_pop = [76.2, 86.4, 113.8, 98.6, 86.6, 102.7, 112.0, 99.0, 64.0, 66.9, 44.1, 25.4, 14.9, 6.8, 1.7, 0.2]; % Sums to ~1000
            if length(cS.initial_pop) ~= cS.aD_new, cS.initial_pop = ones(cS.aD_new,1) * (100/cS.aD_new); warning('initial_pop mismatch, using uniform'); end
            % Survival probabilities for population dynamics (per cS.yearStep years)
            beta_surv_pop = [0.998,0.996,0.994,0.992,0.988,0.984,0.980,0.976,0.970,0.960,0.945,0.920,0.880,0.800,0.680]; % Survival from group j to j+1
            if length(beta_surv_pop) ~= cS.aD_new - 1, error('beta_surv_pop length incorrect for %d groups', cS.aD_new); end
            cS.survivalProbV_popdyn = [beta_surv_pop(:)', 0]; % Last group has 0 survival to next (non-existent) group
            cS.bgp_tolerance = 0.001; cS.bgp_window = 5; cS.max_periods = 50; % For population dynamics simulation

            %% 家庭参数 (年度)
            cS.sigma      = 1.5;
            cS.beta       = 1.011;
            cS.cFloor     = 0.05;
            cS.nSim       = 5e4; % Number of individuals in simulation

            %% 技术参数 (年度)
            cS.A          = 0.895944; cS.alpha = 0.36; cS.ddk = 0.06; % Depreciation rate

            %% 社会保障参数 (年度) - 外生 PAYG 税率 theta
            cS.theta      = 0.20; % Exogenous PAYG payroll tax rate

            %% Heer (2020) Style Government Finance Parameters (年度)
            cS.tau_k = 0.20;       % 资本所得税率 (on r_gross_market - ddk)
            cS.tau_c = 0.10;       % 消费税率
            cS.tau_l = 0.05;       % 额外的劳动所得税率 (on top of PAYG tax cS.theta)
            cS.gov_exp_frac_Y = 0.15; % 政府消费占GDP的比例
            cS.gov_debt_frac_Y = 0.60; % 政府债务占GDP的比例

            %% 劳动禀赋参数
            cS.leSigma1 = 0.38^0.5; cS.leShockStd = 0.045^0.5; cS.lePersistence = 0.96; cS.leWidth = 4; cS.nw = 18;

            %% 网格参数 (kGridV代表非PPS资产的网格)
            cS.tgKY = 3; cS.tgWage = (1-cS.alpha)*cS.A*((cS.tgKY/cS.A)^(cS.alpha/(1-cS.alpha)));
            cS.nk = 50;
            cS.kMin = 0; cS.kMax = 100 * cS.tgWage;
            power = 1.5; kGridV = cS.kMin + (cS.kMax - cS.kMin) * (linspace(0, 1, cS.nk).^power);
            kGridV(1)=cS.kMin; cS.kGridV = kGridV(:);

            %% 年龄效率剖面
            ageEffV_orig_temp = zeros(100, 1);
            ageEffV_orig_temp(20:72) = [linspace(0.3,1.5,36-20+1), 1.5.*ones(1,47-37+1), linspace(1.5,0.2,65-48+1), linspace(0.18,0,72-66+1)];
            cS.ageEffV_orig = ageEffV_orig_temp(cS.age1_orig : cS.ageLast_orig);
            if length(cS.ageEffV_orig) ~= cS.aD_orig, error('ageEffV_orig length mismatch'); end
            cS.ageEffV_new = zeros(cS.aD_new, 1); for a = 1:cS.aD_new, cS.ageEffV_new(a) = mean(cS.ageEffV_orig(cS.physAgeMap{a})); end

            %% 个人养老金计划 (PPS) 参数
            cS.pps_active = true;
            cS.pps_max_contrib_frac = 0.1; % Set to 0 to simulate no PPS initially, or e.g., 0.1 for active PPS
            cS.pps_tax_rate_withdrawal = 0.03; % Tax on PPS withdrawals
            cS.pps_return_rate_premium = 0.0;  % Premium over net market return for PPS assets
            cS.pps_withdrawal_rate = 0.15;     % Annual withdrawal rate from PPS balance during retirement
            cS.pps_contribution_age_max = cS.aR_idx_orig - 1; % Max annual age for PPS contributions
            cS.pps_withdrawal_age_min = cS.aR_idx_orig;   % Min annual age for PPS withdrawals
            cS.pps_in_K = true;         % Whether PPS assets count towards productive capital K
            cS.pps_bequeathable = false; % Whether PPS assets are bequeathable
        end

        function [R_market_gross_factor, MPL_gross] = HHPrices_Huggett(K_productive, L_total_eff, cS)
            % 计算市场总回报率因子 R_market_gross_factor 和总MPL。
            % R_market_gross_factor = 1 + (MPK_gross - cS.ddk)
            % MPL_gross
            if K_productive <= 0, K_productive=1e-6; warning('K_prod non-positive in HHPrices, set to small.'); end
            if L_total_eff <= 0, L_total_eff=1e-6; warning('L_eff non-positive in HHPrices, set to small.'); end

            Y_gross = cS.A * (K_productive^cS.alpha) * (L_total_eff^(1-cS.alpha));
            MPK_gross_val = cS.alpha * Y_gross / K_productive; % This is MPK before depreciation
            MPL_gross = (1-cS.alpha) * Y_gross / L_total_eff;

            R_market_gross_factor = 1 + MPK_gross_val - cS.ddk;
            R_market_gross_factor = max(1.0 + 1e-6, R_market_gross_factor); % Ensure R >= 1
        end

        % --- 劳动禀赋函数 (tauchen, norm_grid, MarkovChainSimulation) - 与V3一致 ---
        function [logGridV, trProbM, prob1V] = EarningProcess_olgm(cS)
             [logGridV_raw, trProbM_calc] = main_olg_v4_utils_pps.tauchen(cS.nw, cS.lePersistence, cS.leShockStd, 0, cS.leWidth);
             gridMin = logGridV_raw(1); gridMax = logGridV_raw(end);
             normGridMin = gridMin - 2*cS.leShockStd; normGridMax = gridMax + 2*cS.leShockStd;
             prob1V_calc = [];
              try
                  [prob1V_calc, ~, ~] = main_olg_v4_utils_pps.norm_grid(logGridV_raw, normGridMin, normGridMax, 0, cS.leSigma1);
                  prob1V_calc = prob1V_calc(:);
              catch ME
                  warning('norm_grid 失败: %s。使用均匀分布。', ME.message);
                  prob1V_calc = ones(cS.nw, 1)/cS.nw;
              end
             logGridV_calc = logGridV_raw(:);
             logGridV_calc = logGridV_calc - logGridV_calc(1) - 1; % Normalize

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
        end
        function [y, trProbM] = tauchen(N, pRho, pSigma, pMu, n_std) % (与V3一致)
             std_y = sqrt(pSigma^2 / (1 - pRho^2));
             if abs(1-pRho)<1e-9, std_y = pSigma*100; end
             y_max = n_std*std_y; y_min = -y_max;
             y_val = linspace(y_min, y_max, N);
             d = 0; if N > 1, d = y_val(2) - y_val(1); end;
             trProbM = zeros(N, N);
             for iR = 1:N
                 for iC = 1:N
                     mn = pRho*y_val(iR);
                     if iC==1
                         trProbM(iR,iC)=normcdf((y_val(1)-mn+d/2)/pSigma);
                     elseif iC==N
                         trProbM(iR,iC)=1-normcdf((y_val(N)-mn-d/2)/pSigma);
                     else
                         trProbM(iR,iC)=normcdf((y_val(iC)-mn+d/2)/pSigma)-normcdf((y_val(iC)-mn-d/2)/pSigma);
                     end
                 end
             end
             rs=sum(trProbM,2); rs(rs<=1e-9)=1;
             trProbM=bsxfun(@rdivide,trProbM,rs);
             um=pMu/(1-pRho);
             if abs(1-pRho)<1e-9 && pMu~=0, um=sign(pMu)*inf; end
             if isfinite(um), y_val=y_val+um; end
             y = y_val;
        end
        function [massV, lbV, ubV] = norm_grid(xV, xMin, xMax, mu, sig) % (与V3一致)
             n=length(xV); xV=xV(:)';
             if n>1 && any(xV(2:n)<xV(1:n-1)), error('xV 非单调递增'); end
             lbV=[]; ubV=[];
             if n>1
                 xM=0.5*(xV(1:n-1)+xV(2:n));
                 lbV=[xMin,xM]; ubV=[xM,xMax];
             elseif n==1
                 lbV=xMin; ubV=xMax;
             else
                 massV=[]; return;
             end
             cdfV=normcdf([lbV,ubV(end)],mu,sig);
             massV=diff(cdfV);
              if any(massV < -1e-9), warning('在 norm_grid 中检测到负质量。'); massV(massV<0)=0; end
             ts=sum(massV);
             if ts>1e-9, massV=massV/ts;
             else
                 if n>0, massV=ones(1,n)/n; else massV=[]; end
             end
             massV=massV(:); lbV=lbV(:); ubV=ubV(:);
             if n>1 && any(ubV<lbV), error('ubV < lbV'); end
        end
        function eIdxM = MarkovChainSimulation(nSim, T_sim, prob0V, trProbM, rvInM) % (与V3一致)
             ns=length(prob0V);
             if size(trProbM,1)~=ns||size(trProbM,2)~=ns,error('转移矩阵维度不匹配');end;
             if abs(sum(prob0V)-1)>1e-5,prob0V=prob0V./sum(prob0V);end;
             if any(abs(sum(trProbM,2)-1)>1e-5), row_sums=sum(trProbM,2); row_sums(row_sums<=1e-9)=1; trProbM=bsxfun(@rdivide,trProbM,row_sums); end;
             if size(rvInM,1)~=nSim||size(rvInM,2)~=T_sim,error('随机数矩阵维度不匹配');end

             cP0=cumsum(prob0V(:)');
             cPT=cumsum(trProbM,2); cPT(:,ns)=1.0;
             eIdxM=zeros(nSim,T_sim,'uint16');
             eIdxM(:,1)=1+sum(bsxfun(@gt,rvInM(:,1),cP0),2);

             for t_mc=1:(T_sim-1)
                 cSI=eIdxM(:,t_mc);
                 vPI=(cSI>=1)&(cSI<=ns);
                 if ~all(vPI),cSI(~vPI)=1;eIdxM(:,t_mc)=cSI;end;
                 cPt=cPT(cSI,:);
                 eIdxM(:,t_mc+1)=1+sum(bsxfun(@gt,rvInM(:,t_mc+1),cPt),2);
             end
        end
        % --- 结束劳动禀赋 ---

        function [HHlaborM_group, L] = LaborSupply_Huggett(eIdxM, cS, paramS, Z_ss_norm_group) % (与V3一致)
            nSim = size(eIdxM, 1);
            HHlaborM_group = zeros(nSim, cS.aD_new);
            ageToGroupMap = zeros(cS.aD_orig, 1);
            for a_new_map = 1:cS.aD_new
                ageToGroupMap(cS.physAgeMap{a_new_map}) = a_new_map;
            end
            leGridV_col = paramS.leGridV(:);
            HHlaborM_annual_temp = zeros(nSim, cS.aD_orig);
            for a_orig = 1 : cS.aD_orig
               if a_orig < cS.aR_idx_orig % Only working ages
                   a_new = ageToGroupMap(a_orig);
                   HHlaborM_annual_temp(:, a_orig) = cS.ageEffV_new(a_new) .* leGridV_col(eIdxM(:, a_orig));
               end
            end
            for a_new = 1:cS.aD_new
                annual_indices = cS.physAgeMap{a_new};
                if ~isempty(annual_indices)
                    HHlaborM_group(:, a_new) = mean(HHlaborM_annual_temp(:, annual_indices), 2);
                end
            end
            L_working_groups = mean(HHlaborM_group(:, 1:cS.aR_new), 1) * Z_ss_norm_group(1:cS.aR_new);
            L = max(0, L_working_groups);
        end

        function eIdxM = LaborEndowSimulation_olgm(cS, paramS) % (与V3一致)
            rng(433);
            rvInM=rand([cS.nSim,cS.aD_orig]);
            eIdxM = main_olg_v4_utils_pps.MarkovChainSimulation(cS.nSim, cS.aD_orig, paramS.leProb1V, paramS.leTrProbM, rvInM);
        end

        function [kHistM, kPpsHistM, cHistM] = HHSimulation_olgm(kPolM, cPpsPolM, cPolM_consump, eIdxM, R_k_net_factor_hh, w_net_hh, TR_total, bV_payg, paramS, cS) % (与V3一致)
            validateattributes(kPolM, {'double'}, {'finite', 'nonnan', 'nonempty', 'real', 'size', [cS.nk, cS.nw, cS.aD_new]});
            validateattributes(cPpsPolM, {'double'}, {'finite', 'nonnan', 'nonempty', 'real','>=',0, 'size', [cS.nk, cS.nw, cS.aD_new]});
            validateattributes(cPolM_consump, {'double'}, {'finite', 'nonnan', 'nonempty', 'real', 'size', [cS.nk, cS.nw, cS.aD_new]});
            validateattributes(eIdxM, {'uint16', 'double'}, {'nonempty', 'integer', 'positive', '<=', cS.nw, 'size', [cS.nSim, cS.aD_orig]});

            ageToGroupMap = zeros(cS.aD_orig, 1);
            for a_new = 1:cS.aD_new
                annual_indices = cS.physAgeMap{a_new};
                if ~isempty(annual_indices), ageToGroupMap(annual_indices) = a_new; end
            end

            nSim   = size(eIdxM, 1);
            kHistM    = zeros(nSim, cS.aD_orig + 1);
            kPpsHistM = zeros(nSim, cS.aD_orig + 1);
            cHistM    = zeros(nSim, cS.aD_orig); % Stores consumption quantity
            kGridV_col = cS.kGridV(:);

            kPolInterp = cell(cS.nw, cS.aD_new);
            cPpsPolInterp = cell(cS.nw, cS.aD_new);
            cPolConsumpInterp = cell(cS.nw, cS.aD_new);
             for a_new_interp=1:cS.aD_new
                 for ie_interp=1:cS.nw
                      kPolInterp{ie_interp, a_new_interp}    = griddedInterpolant(kGridV_col, kPolM(:, ie_interp, a_new_interp), 'linear', 'nearest');
                      cPpsPolInterp{ie_interp, a_new_interp} = griddedInterpolant(kGridV_col, cPpsPolM(:, ie_interp, a_new_interp), 'linear', 'nearest');
                      cPolConsumpInterp{ie_interp, a_new_interp} = griddedInterpolant(kGridV_col, cPolM_consump(:, ie_interp, a_new_interp), 'linear', 'nearest');
                 end
             end

            pps_return_net_annual = (R_k_net_factor_hh - 1) + cS.pps_return_rate_premium;

            for a_orig = 1 : cS.aD_orig
                a_new_group_idx = ageToGroupMap(a_orig);
                kNowV = kHistM(:, a_orig);
                kPpsNowV = kPpsHistM(:, a_orig);

                kNextNonPpsV = zeros(nSim, 1);
                cPpsContribV = zeros(nSim, 1);
                cConsumpValV = zeros(nSim, 1); % This will be consumption quantity
                kPpsNextV    = zeros(nSim, 1);

                is_working = (a_orig < cS.aR_idx_orig);

                pps_withdrawal_pretax_sim = zeros(nSim,1); % Pre-tax withdrawal from PPS account

                if ~is_working && cS.pps_active && a_orig >= cS.pps_withdrawal_age_min && cS.pps_max_contrib_frac > 1e-9 % PPS active and relevant
                    pps_withdrawal_pretax_sim = kPpsNowV * cS.pps_withdrawal_rate;
                end

                for ie = 1 : cS.nw
                    simIdxV = find(eIdxM(:, a_orig) == ie);
                    if ~isempty(simIdxV)
                         kNow_ie_clamped = max(kGridV_col(1), min(kGridV_col(end), kNowV(simIdxV)));

                         kNextNonPpsV(simIdxV) = kPolInterp{ie, a_new_group_idx}(kNow_ie_clamped);
                         cPpsContribV(simIdxV) = cPpsPolInterp{ie, a_new_group_idx}(kNow_ie_clamped);
                         cConsumpValV(simIdxV) = cPolConsumpInterp{ie, a_new_group_idx}(kNow_ie_clamped); % This is c_quantity

                         if is_working && cS.pps_active && a_orig <= cS.pps_contribution_age_max && cS.pps_max_contrib_frac > 1e-9
                             % VFI has already determined cPpsContribV based on budget and rules
                         else
                             cPpsContribV(simIdxV) = 0; % No contribution if not working age or PPS not active for contribution
                         end

                         if cS.pps_active && cS.pps_max_contrib_frac > 1e-9
                             kPpsNextV(simIdxV) = (kPpsNowV(simIdxV) + cPpsContribV(simIdxV) - pps_withdrawal_pretax_sim(simIdxV)) * (1 + pps_return_net_annual);
                             kPpsNextV(simIdxV) = max(0, kPpsNextV(simIdxV));
                         else
                             kPpsNextV(simIdxV) = kPpsNowV(simIdxV); % PPS balance remains if PPS not active for contributions/withdrawals this period
                         end
                    end
                end

                 kHistM(:, a_orig + 1)    = max(cS.kMin, min(cS.kMax, kNextNonPpsV));
                 kPpsHistM(:, a_orig + 1) = max(0, kPpsNextV);
                 cHistM(:, a_orig)        = max(cS.cFloor, cConsumpValV); % Store consumption quantity

                 kHistM(~isfinite(kHistM(:, a_orig + 1)), a_orig + 1) = cS.kMin;
                 kPpsHistM(~isfinite(kPpsHistM(:, a_orig + 1)), a_orig + 1) = 0;
                 cHistM(~isfinite(cHistM(:, a_orig)), a_orig) = cS.cFloor;
            end

            kHistM = kHistM(:, 1:cS.aD_orig);
            kPpsHistM = kPpsHistM(:, 1:cS.aD_orig);

            validateattributes(kHistM, {'double'}, {'finite', 'nonnan', 'nonempty', 'real', '>=', cS.kMin - 1e-6, '<=', cS.kMax + 1e-6, 'size', [nSim, cS.aD_orig]});
            validateattributes(kPpsHistM, {'double'}, {'finite', 'nonnan', 'nonempty', 'real', '>=', -1e-6, 'size', [nSim, cS.aD_orig]});
            validateattributes(cHistM, {'double'}, {'finite', 'nonnan', 'nonempty', 'real', '>=', cS.cFloor - 1e-6, 'size', [nSim, cS.aD_orig]});
        end


        function [muM, utilM] = CES_utility(cM_quantity, sig, cS) % (与V3一致)
             if ~isscalar(sig)||sig<=0,error('Sigma 必须是正标量。');end;
             min_c_quantity = cS.cFloor;
             vld = (cM_quantity >= min_c_quantity);
             c_c_quantity = max(min_c_quantity, cM_quantity);
             utilM = -Inf(size(cM_quantity)); muM = Inf(size(cM_quantity));

             if abs(sig-1)<1e-6
                 utilM(vld) = log(c_c_quantity(vld));
                 muM(vld) = 1 ./ c_c_quantity(vld);
             else
                 utilM(vld) = (c_c_quantity(vld).^(1-sig)) ./ (1-sig);
                 muM(vld) = c_c_quantity(vld).^(-sig);
             end
             utilM(~vld) = -1e10 - (min_c_quantity - cM_quantity(~vld))*1e10;
             muM(~vld) = c_c_quantity(~vld).^(-sig) + 1e10;
        end

         function incomeM_resources = HHIncome_Huggett(k_non_pps_grid, R_k_net_factor_hh_price, w_net_hh_price, TR_total_price, b_payg_price_group, a_new_idx_income, paramS_income, cS) % (与V3一致)
            nk_income = length(k_non_pps_grid);
            leGridV_col_income = paramS_income.leGridV(:);
            nonCapIncomeNetM_income = zeros(1, cS.nw);

            if a_new_idx_income <= cS.aR_new % 工作年龄组
               net_wage_income_vec_income = w_net_hh_price * cS.ageEffV_new(a_new_idx_income) .* leGridV_col_income;
               nonCapIncomeNetM_income(1,:) = net_wage_income_vec_income' + TR_total_price + b_payg_price_group;
            else % 退休年龄组
               nonCapIncomeNetM_income(1,:) = b_payg_price_group + TR_total_price;
            end

            % Total resources available before PPS contributions and before consumption tax on c
            incomeM_resources = R_k_net_factor_hh_price * k_non_pps_grid(:) * ones(1, cS.nw) + ones(nk_income, 1) * nonCapIncomeNetM_income;

            if any(~isfinite(incomeM_resources(:))), warning('Non-finite income in HHIncome_Huggett group %d.',a_new_idx_income); incomeM_resources(~isfinite(incomeM_resources)) = 1e-6; end
         end

        function [cPolM_quantity, kPolM, cPpsPolM, valueM] = HHSolution_VFI_Huggett(R_k_net_factor_vfi, w_net_hh_vfi, TR_total_vfi, bV_payg_vfi, paramS_vfi, cS) % (与V3一致)
            cPolM_quantity  = zeros(cS.nk, cS.nw, cS.aD_new); % Stores consumption quantity
            kPolM  = zeros(cS.nk, cS.nw, cS.aD_new);
            cPpsPolM = zeros(cS.nk, cS.nw, cS.aD_new);
            valueM = zeros(cS.nk, cS.nw, cS.aD_new);

            for a_vfi_loop_main = cS.aD_new : -1 : 1
               vPrimeM_vfi = [];
               if a_vfi_loop_main < cS.aD_new
                   vPrimeM_vfi = valueM(:,:,a_vfi_loop_main+1);
               end
               [cPolM_quantity(:,:,a_vfi_loop_main), kPolM(:,:,a_vfi_loop_main), cPpsPolM_temp, valueM(:,:,a_vfi_loop_main)] = ...
                  main_olg_v4_utils_pps.HHSolutionByAge_VFI_Huggett(a_vfi_loop_main, vPrimeM_vfi, ...
                  R_k_net_factor_vfi, w_net_hh_vfi, TR_total_vfi, bV_payg_vfi(a_vfi_loop_main), paramS_vfi, cS);

               if cS.pps_active && cS.pps_max_contrib_frac > 1e-9
                   cPpsPolM(:,:,a_vfi_loop_main) = cPpsPolM_temp;
               else
                   cPpsPolM(:,:,a_vfi_loop_main) = 0; % Ensure it's zero if PPS not meaningfully active
               end
            end
            % cPolM is quantity, so floor is cS.cFloor
            validateattributes(cPolM_quantity, {'double'}, {'finite', 'nonnan', 'real','>=',cS.cFloor-1e-6,'size',[cS.nk,cS.nw,cS.aD_new]});
            validateattributes(kPolM, {'double'}, {'finite', 'nonnan', 'real','>=',cS.kMin-1e-6,'<=',cS.kMax+1e-6,'size',[cS.nk,cS.nw,cS.aD_new]});
            validateattributes(cPpsPolM, {'double'}, {'finite', 'nonnan', 'real','>=',0,'size',[cS.nk,cS.nw,cS.aD_new]});
            validateattributes(valueM, {'double'}, {'finite', 'nonnan', 'real','size',[cS.nk,cS.nw,cS.aD_new]});
        end

        function [cPolM_age_quantity, kPolM_age, cPpsPolM_age_actual, valueM_age_out] = HHSolutionByAge_VFI_Huggett(...
            a_idx, vPrime_keM_age, R_k_net_factor_age, w_net_hh_age, TR_total_age, b_age_group, paramS_age, cS) % (与V3一致)

            if a_idx < cS.aD_new && isempty(vPrime_keM_age), error('VFI by Age: vPrime empty a<aD'); end
            if a_idx == cS.aD_new && ~isempty(vPrime_keM_age), error('VFI by Age: vPrime not empty a=aD'); end
            if a_idx < cS.aD_new, validateattributes(vPrime_keM_age,{'double'},{'finite','nonnan','real','size',[cS.nk,cS.nw]}); end

            % budgetM_before_pps_age is resources available before PPS choice and before c*tau_c tax
            budgetM_before_pps_age = main_olg_v4_utils_pps.HHIncome_Huggett(cS.kGridV, R_k_net_factor_age, w_net_hh_age, TR_total_age, b_age_group, a_idx, paramS_age, cS);

            fminOpt_age=optimset('fminbnd'); fminOpt_age.TolX=1e-6; fminOpt_age.Display='none';

            cPolM_age_quantity = zeros(cS.nk, cS.nw); % Stores consumption quantity
            kPolM_age = zeros(cS.nk, cS.nw);
            cPpsPolM_age_actual = zeros(cS.nk, cS.nw);
            valueM_age_out = -Inf(cS.nk, cS.nw);
            is_pps_contrib_age_val = (a_idx <= cS.aR_new && a_idx <= (cS.pps_contribution_age_max - cS.age1_orig)/cS.yearStep +1 ) ; % Check if in PPS contribution age range (group index)

            if a_idx == cS.aD_new % Last period
               cPolM_age_quantity = max(cS.cFloor, budgetM_before_pps_age / (1 + cS.tau_c) );
               kPolM_age = cS.kMin*ones(cS.nk,cS.nw);
               cPpsPolM_age_actual = zeros(cS.nk, cS.nw);
               [~, val_temp_age] = main_olg_v4_utils_pps.CES_utility(cPolM_age_quantity, cS.sigma, cS);
               valueM_age_out = val_temp_age;
               valueM_age_out(~isfinite(valueM_age_out)) = -1e12;
            else
               ExValuePrime_k_given_e_age = zeros(cS.nk, cS.nw);
                for ikPrime_age = 1:cS.nk, ExValuePrime_k_given_e_age(ikPrime_age, :) = vPrime_keM_age(ikPrime_age, :) * paramS_age.leTrProbM'; end
                ExValuePrime_k_given_e_age(~isfinite(ExValuePrime_k_given_e_age))=-1e10;
               vPInt_age=cell(1,cS.nw);
               for ie_age_interp=1:cS.nw
                   vPInt_age{ie_age_interp} = griddedInterpolant(cS.kGridV, ExValuePrime_k_given_e_age(:, ie_age_interp), 'linear','linear');
               end

               parfor ik_age_loop = 1:cS.nk
                   c_row_age_q=zeros(1,cS.nw); k_row_age=zeros(1,cS.nw); cpps_row_for_this_k_state=zeros(1,cS.nw); v_row_age=zeros(1,cS.nw);
                   for ie_age_loop_inner = 1 : cS.nw
                        vPofK_interpolant_age = vPInt_age{ie_age_loop_inner};
                        budget_ike_before_pps_age = budgetM_before_pps_age(ik_age_loop,ie_age_loop_inner);
                        cPps_chosen_this_state_val = 0.0;

                        if is_pps_contrib_age_val && cS.pps_active && cS.pps_max_contrib_frac > 1e-9
                           wage_income_for_pps_limit_age = w_net_hh_age * cS.ageEffV_new(a_idx) * paramS_age.leGridV(ie_age_loop_inner);
                           max_contrib_pps_based_on_wage = wage_income_for_pps_limit_age * cS.pps_max_contrib_frac;
                           
                           % Max cPps also limited by budget: budget_ike - (c_floor_cost + k_min)
                           min_expenditure_on_c_and_k_prime = (cS.cFloor * (1 + cS.tau_c)) + cS.kMin;
                           max_permissible_cPps_from_budget = budget_ike_before_pps_age - min_expenditure_on_c_and_k_prime;
                           
                           cPps_chosen_this_state_val = min(max_contrib_pps_based_on_wage, max(0, max_permissible_cPps_from_budget));
                           cPps_chosen_this_state_val = max(0, cPps_chosen_this_state_val); % Ensure non-negative
                        end
                        cpps_row_for_this_k_state(ie_age_loop_inner) = cPps_chosen_this_state_val;
                        budget_remaining_for_k_prime_and_c_expenditure = budget_ike_before_pps_age - cPps_chosen_this_state_val;

                        [c_row_age_q(ie_age_loop_inner), k_row_age(ie_age_loop_inner), v_row_age(ie_age_loop_inner)] = ...
                              main_olg_v4_utils_pps.HHSolutionByOneState_OptK(a_idx, budget_remaining_for_k_prime_and_c_expenditure, vPofK_interpolant_age, fminOpt_age, cS, paramS_age);
                   end
                   cPolM_age_quantity(ik_age_loop,:) = c_row_age_q; kPolM_age(ik_age_loop,:) = k_row_age;
                   cPpsPolM_age_actual(ik_age_loop,:) = cpps_row_for_this_k_state; valueM_age_out(ik_age_loop,:) = v_row_age;
               end
            end
        end

        function [c_quantity, kPrime, ValueFunc] = HHSolutionByOneState_OptK(a_new_osa_k, budget_remaining_expenditure_osa_k, vPofK_int_osa_k, fminOpt_osa_k, cS, paramS_osa_k) % (与V3一致)
            kPMin_osa_k = cS.kMin;
            kPMax_osa_k = budget_remaining_expenditure_osa_k - (cS.cFloor * (1 + cS.tau_c));
            ValueFunc = -Inf;

            function negV_nested = negBellmanObjective_nested(kP_nested, ag_obj_nested, budget_obj_nested, evInt_obj_nested)
                [negV_nested, ~] = BellmanInner_nested(kP_nested, ag_obj_nested, budget_obj_nested, evInt_obj_nested);
            end

            function [negVal_nested, Val_nested] = BellmanInner_nested(kP_inner_nested, ag_inner_nested, budget_exp_inner_nested, evInt_inner_nested)
                cons_quantity_nested = max(cS.cFloor, (budget_exp_inner_nested - kP_inner_nested) / (1 + cS.tau_c) );
                [~, util_val_nested] = main_olg_v4_utils_pps.CES_utility(cons_quantity_nested, cS.sigma, cS);
                if ~isfinite(util_val_nested), util_val_nested = -1e12; end; Val_nested = util_val_nested;
                if ag_inner_nested < cS.aD_new
                    evF_val_nested = -Inf;
                    try
                        kP_eval_nested = max(cS.kGridV(1),min(cS.kGridV(end),kP_inner_nested));
                        evF_val_nested=evInt_inner_nested(kP_eval_nested);
                    catch ME_bellman_nested
                        warning('BellmanInner_nested interp error age %d, kP %.2e: %s', ag_inner_nested, kP_inner_nested, ME_bellman_nested.message);
                        if kP_inner_nested<cS.kGridV(1),evF_val_nested=evInt_inner_nested(cS.kGridV(1));else,evF_val_nested=evInt_inner_nested(cS.kGridV(end));end
                    end
                    if ~isfinite(evF_val_nested),evF_val_nested=-1e12;end
                    s_transition_nested = cS.s_1yr_transitionV(ag_inner_nested);
                    Val_nested = util_val_nested + cS.beta*s_transition_nested*evF_val_nested;
                end
                if ~isfinite(Val_nested), negVal_nested=1e12; Val_nested=-1e12; else negVal_nested = -Val_nested; end
            end

            if kPMax_osa_k <= kPMin_osa_k
                kPrime = kPMin_osa_k;
                c_quantity = max(cS.cFloor, (budget_remaining_expenditure_osa_k - kPrime) / (1 + cS.tau_c) );
                 [~, u_c_osa_k] = main_olg_v4_utils_pps.CES_utility(c_quantity, cS.sigma, cS);
                if a_new_osa_k < cS.aD_new
                    EV_c_osa_k = -Inf; try kP_eval_osa_k = max(cS.kGridV(1),min(cS.kGridV(end),kPrime)); EV_c_osa_k=vPofK_int_osa_k(kP_eval_osa_k); catch ME_osa_k; warning('OptK corner interp error age %d: %s',a_new_osa_k,ME_osa_k.message); if kPrime<cS.kGridV(1),EV_c_osa_k=vPofK_int_osa_k(cS.kGridV(1));else,EV_c_osa_k=vPofK_int_osa_k(cS.kGridV(end));end;end; if ~isfinite(EV_c_osa_k),EV_c_osa_k=-1e12;end
                    s_transition_osa_k = cS.s_1yr_transitionV(a_new_osa_k); ValueFunc = u_c_osa_k + cS.beta*s_transition_osa_k*EV_c_osa_k;
                else, ValueFunc = u_c_osa_k; end
                if ~isfinite(ValueFunc), ValueFunc = -1e12; end
            else
                kPMin_opt_osa_k = max(cS.kMin, kPMin_osa_k);
                kPMax_opt_osa_k = max(kPMin_opt_osa_k + 1e-9, min(cS.kMax, kPMax_osa_k));

                if kPMin_opt_osa_k >= kPMax_opt_osa_k
                    kPrime_opt_osa_k = kPMin_opt_osa_k;
                    [negV_osa_k, ~] = BellmanInner_nested(kPrime_opt_osa_k, a_new_osa_k, budget_remaining_expenditure_osa_k, vPofK_int_osa_k);
                    ValueFunc = -negV_osa_k;
                else
                    obj_osa_k = @(kP_osa_k) negBellmanObjective_nested(kP_osa_k, a_new_osa_k, budget_remaining_expenditure_osa_k, vPofK_int_osa_k);
                    [kPrime_opt_osa_k, negV_osa_k, eflag_osa_k] = fminbnd(obj_osa_k, kPMin_opt_osa_k, kPMax_opt_osa_k, fminOpt_osa_k);

                    if eflag_osa_k <= 0 || abs(kPrime_opt_osa_k-kPMin_opt_osa_k)<1e-7 || abs(kPrime_opt_osa_k-kPMax_opt_osa_k)<1e-7
                        [nV_min_osa_k,~]=BellmanInner_nested(kPMin_opt_osa_k,a_new_osa_k,budget_remaining_expenditure_osa_k,vPofK_int_osa_k);
                        [nV_max_osa_k,~]=BellmanInner_nested(kPMax_opt_osa_k,a_new_osa_k,budget_remaining_expenditure_osa_k,vPofK_int_osa_k);
                        if nV_min_osa_k <= nV_max_osa_k + 1e-9
                             kPrime_opt_osa_k=kPMin_opt_osa_k; negV_osa_k=nV_min_osa_k;
                        else
                             kPrime_opt_osa_k=kPMax_opt_osa_k; negV_osa_k=nV_max_osa_k;
                        end
                        if eflag_osa_k <= 0 && ~(abs(kPMin_opt_osa_k - kPMax_opt_osa_k) < 1e-8) % Avoid warning if bounds are virtually identical
                           % warning('fminbnd issue OptK age %d, eflag %d, k_bounds [%.2e, %.2e], budget %.2e',a_new_osa_k,eflag_osa_k, kPMin_opt_osa_k, kPMax_opt_osa_k, budget_remaining_expenditure_osa_k);
                        end
                    end
                    ValueFunc = -negV_osa_k;
                end
                kPrime = kPrime_opt_osa_k;
                c_quantity = max(cS.cFloor, (budget_remaining_expenditure_osa_k - kPrime) / (1 + cS.tau_c) );
                if ~isfinite(ValueFunc), ValueFunc = -1e12; end
            end

            kPrime = max(cS.kMin, min(cS.kMax, kPrime));
            c_quantity = max(cS.cFloor, (budget_remaining_expenditure_osa_k - kPrime) / (1 + cS.tau_c) );

            [~, u_final_c_osa_k] = main_olg_v4_utils_pps.CES_utility(c_quantity, cS.sigma, cS);
            if a_new_osa_k < cS.aD_new
                EV_A_final_osa_k = -Inf;
                try
                    kP_eval_final_osa_k = max(cS.kGridV(1),min(cS.kGridV(end),kPrime));
                    EV_A_final_osa_k=vPofK_int_osa_k(kP_eval_final_osa_k);
                catch
                    if kPrime<cS.kGridV(1),EV_A_final_osa_k=vPofK_int_osa_k(cS.kGridV(1));else,EV_A_final_osa_k=vPofK_int_osa_k(cS.kGridV(end));end;
                end;
                if ~isfinite(EV_A_final_osa_k),EV_A_final_osa_k=-1e12;end
                s_transition_final_osa_k = cS.s_1yr_transitionV(a_new_osa_k);
                ValueFunc = u_final_c_osa_k + cS.beta*s_transition_final_osa_k*EV_A_final_osa_k;
            else
                ValueFunc = u_final_c_osa_k;
            end

            if ~isfinite(kPrime), kPrime = cS.kMin; end;
            if ~isfinite(c_quantity), c_quantity = cS.cFloor; end;
            if ~isfinite(ValueFunc), ValueFunc = -1e12; end

            validateattributes(kPrime,{'double'},{'finite','nonnan','scalar','real','>=',cS.kMin-1e-6,'<=',cS.kMax+1e-6});
            validateattributes(c_quantity,{'double'},{'finite','nonnan','scalar','real','>=',cS.cFloor-1e-6});
            validateattributes(ValueFunc,{'double'},{'finite','nonnan','scalar','real'});
        end

        % Unused static helper functions from V3, kept for reference or potential future use.
        % If used, they would also need cS.tau_c modifications similar to the nested Bellman.
        function negV = negBellmanObjective_k_pps(kP, ag_obj, budget_rem_obj, evInt_obj, cS_in_obj, pS_in_obj) % (与V3一致, 但未使用)
            [negV, ~] = main_olg_v4_utils_pps.BellmanInner_k_pps(kP, ag_obj, budget_rem_obj, evInt_obj, cS_in_obj, pS_in_obj);
        end

        function [negVal, Val] = BellmanInner_k_pps(kP_inner, ag_inner, budget_rem_inner, evInt_inner, cS_in, pS_in) % (与V3一致, 但未使用)
            % This version uses budget_rem_inner directly for consumption, assuming it's post-tax-on-c.
            % The main VFI path uses budget_remaining_expenditure which is pre-tax-on-c.
            cons = max(cS_in.cFloor, budget_rem_inner - kP_inner); % Assumes budget_rem_inner is for c_quantity + kP
            [util_val, ~] = main_olg_v4_utils_pps.CES_utility(cons, cS_in.sigma, cS_in);
            if ~isfinite(util_val), util_val = -1e12; end; Val = util_val;
            if ag_inner < cS_in.aD_new
                evF_val = -Inf;
                try kP_eval_val = max(cS_in.kGridV(1),min(cS_in.kGridV(end),kP_inner)); evF_val=evInt_inner(kP_eval_val);
                catch ME_bellman_pps; warning('BellmanInner_k_pps interp error age %d: %s', ag_inner, ME_bellman_pps.message); if kP_inner<cS_in.kGridV(1),evF_val=evInt_inner(cS_in.kGridV(1));else,evF_val=evInt_inner(cS_in.kGridV(end));end; end
                if ~isfinite(evF_val),evF_val=-1e12;end
                s_transition_val = cS_in.s_1yr_transitionV(ag_inner); Val = util_val + cS_in.beta*s_transition_val*evF_val;
            end
            if ~isfinite(Val), negVal=1e12; Val=-1e12; else negVal = -Val; end
        end

    end % End methods (Static)
end % End classdef

% --- END OF FILE main_olg_v4_utils_pps.m ---