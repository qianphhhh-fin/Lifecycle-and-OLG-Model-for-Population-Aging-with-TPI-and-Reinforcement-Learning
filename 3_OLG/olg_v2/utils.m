% --- START OF FILE utils.m ---

classdef utils
    % Utility functions for OLG model v2 (MODIFIED to Mixed Time Units)
    % - VFI over 5-year groups using ANNUAL beta/prices/survival
    % - Simulation and Aggregation using ANNUAL results/masses
    % - Population dynamics still determine steady-state GROUP mass distribution

    methods (Static)

        % =====================================================================
        % == Population Dynamics Functions (Group Level - UNCHANGED) =========
        % =====================================================================

        function popS = initPopulation(cS)
            % Initialize population structure (counts for groups)
            popS.Z = zeros(cS.aD_new, 1); % Use aD_new for groups
            initial_total = sum(cS.initial_pop); % initial_pop is per group
            if initial_total > 0 && length(cS.initial_pop) == cS.aD_new
                 popS.Z(:, 1) = cS.initial_pop(:) / initial_total * 100; % Scale to sum to 100
            else
                warning('Initial population mismatch or zero sum. Setting uniform initial GROUP population.');
                 popS.Z(:, 1) = 100 / cS.aD_new;
            end
            popS.totalPop = sum(popS.Z(:, 1));
            if popS.totalPop > 1e-9
                 popS.ageDist = popS.Z(:, 1) / popS.totalPop; % Group distribution
            else
                 popS.ageDist = zeros(cS.aD_new, 1);
            end
            % Store initial GROUP distribution
            popS.initialAgeDist = popS.ageDist;
            fprintf('Initial GROUP population set. Total=%.2f\n', popS.totalPop);
        end

        function popS = populationDynamics(popS, cS)
            % Simulate GROUP population dynamics using multi-year survival/growth
            % Uses cS.survivalProbV_popdyn (multi-year survival)

            max_periods_sim = cS.max_periods;
            Z_history = zeros(cS.aD_new, max_periods_sim + 1); % Group history
            totalPop_history = zeros(1, max_periods_sim + 1);
            ageDist_history = zeros(cS.aD_new, max_periods_sim + 1); % Group distribution history

            Z_history(:, 1) = popS.Z(:, 1);
            totalPop_history(1) = popS.totalPop(1);
            ageDist_history(:, 1) = popS.ageDist(:, 1);

            fprintf('Population dynamics simulation starting (Groups, Max Periods = %d)...\n', max_periods_sim);
            bgp_reached_flag = false;
            actual_periods_run = max_periods_sim;

            for t = 1:max_periods_sim
                if mod(t, 10) == 0 || t == 1
                    fprintf('  Simulating population period %d (groups)\n', t);
                end

                Z_current_period = Z_history(:, t);
                Z_next_period = zeros(cS.aD_new, 1);

                % Growth/survival logic remains the same as original V2 (group level)
                 time_varying_growth_rate = 0;
                 if t < 6, time_varying_growth_rate = -0.01 - 0.003 * t;
                 else, time_varying_growth_rate = -0.03 - 0.004 * min(t-6, 10); end
                 Z_next_period(1) = Z_current_period(1) * (1 + time_varying_growth_rate);
                 Z_next_period(1) = max(0, Z_next_period(1));

                for a = 2:cS.aD_new
                    survival_prob = 0;
                    if (a-1) >= 1 && (a-1) <= length(cS.survivalProbV_popdyn)
                         survival_prob = cS.survivalProbV_popdyn(a-1); % Use multi-year popdyn survival
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

                % Check for steady state convergence (groups)
                current_check_period = t + 1;
                if current_check_period >= cS.bgp_window + 1
                    stable = true;
                    for w = 1:cS.bgp_window
                       hist_idx1 = current_check_period - w + 1; hist_idx2 = current_check_period - w;
                       if hist_idx1 > 0 && hist_idx2 > 0
                           change = norm(ageDist_history(:, hist_idx1) - ageDist_history(:, hist_idx2));
                           if change >= cS.bgp_tolerance, stable = false; break; end
                       else stable = false; break; end
                    end
                    if stable
                        fprintf('\nPopulation steady state (groups) reached at period %d.\n', t);
                        bgp_reached_flag = true; actual_periods_run = t; break;
                    end
                end
            end % End simulation loop t

            final_period_idx = actual_periods_run + 1;
            popS.Z = Z_history(:, 1:final_period_idx);
            popS.totalPop = totalPop_history(1:final_period_idx);
            popS.ageDist = ageDist_history(:, 1:final_period_idx); % Group distribution history

            % Calculate dependency ratio based on groups
            depRatio_history = zeros(1, actual_periods_run);
             for th = 1:actual_periods_run
                 Z_t = Z_history(:, th);
                 % Use aR_new (index of last working group)
                 working_pop = sum(Z_t(1:cS.aR_new));
                 retired_pop = sum(Z_t(cS.aR_new+1:end));
                 if working_pop > 1e-9, depRatio_history(th) = retired_pop / working_pop; else depRatio_history(th) = inf; end
             end
            popS.dependencyRatio = depRatio_history;

            fprintf('Population dynamics simulation finished. Periods run: %d. BGP reached: %d\n', actual_periods_run, bgp_reached_flag);
            if ~bgp_reached_flag, fprintf('Warning: Population steady state (groups) not reached within %d periods.\n', max_periods_sim); end
        end

        function [Z_ss, dependency_ratio_ss, bgp_reached, bgp_period] = detectSteadyStatePopulation(popS, cS)
            % Detect steady state based on GROUP age distribution stability
            % PLOT initial vs steady-state GROUP distribution comparison.

             actual_periods_in_data = size(popS.Z, 2);
             bgp_reached = false;
             bgp_period = actual_periods_in_data - 1;

             if actual_periods_in_data < cS.bgp_window + 1
                 fprintf('Population simulation too short (%d data points) for SS check (window = %d).\n', actual_periods_in_data, cS.bgp_window);
             else
                 fprintf('Checking for population steady state (groups, last %d periods)...\n', cS.bgp_window);
                 for t_check_end = actual_periods_in_data : -1 : cS.bgp_window + 1
                     stable = true;
                     for w = 0 : (cS.bgp_window - 1)
                         idx1 = t_check_end - w; idx2 = t_check_end - w - 1;
                         if idx1 > 0 && idx2 > 0
                             change = norm(popS.ageDist(:, idx1) - popS.ageDist(:, idx2)); % Compare group distributions
                             if change >= cS.bgp_tolerance, stable = false; break; end
                         else stable = false; break; end
                     end
                     if stable
                         bgp_reached = true; bgp_period = t_check_end - 1;
                         fprintf('Population steady state (groups) detected starting period %d (data index %d).\n', bgp_period, t_check_end);
                         break;
                     end
                 end
                  if ~bgp_reached, fprintf('Population steady state (groups) not detected.\n'); end
             end

             ss_data_index = bgp_period + 1;
             Z_ss = popS.Z(:, ss_data_index); % Steady-state group counts
             Z_ss_norm = zeros(cS.aD_new, 1);
             if sum(Z_ss) > 1e-9, Z_ss_norm = Z_ss / sum(Z_ss); end % Normalized group mass

             if isfield(popS, 'dependencyRatio') && length(popS.dependencyRatio) >= bgp_period && bgp_period > 0
                 dependency_ratio_ss = popS.dependencyRatio(bgp_period);
             else
                  working_pop_ss = sum(Z_ss(1:cS.aR_new));
                  retired_pop_ss = sum(Z_ss(cS.aR_new+1:end));
                  if working_pop_ss > 1e-9, dependency_ratio_ss = retired_pop_ss / working_pop_ss; else dependency_ratio_ss = inf; end
             end

             % --- Plot Initial vs Steady-State GROUP Distribution ---
             figure('Name', 'V2 Mod: Initial vs Steady-State GROUP Population');
             hold on;
             group_indices = 1:cS.aD_new;
             if isfield(popS, 'initialAgeDist')
                 bar(group_indices - 0.2, popS.initialAgeDist * 100, 0.4, 'b', 'DisplayName', 'Initial Group Dist');
             else
                 warning('Initial age distribution not found for plotting.');
             end
             bar(group_indices + 0.2, Z_ss_norm * 100, 0.4, 'r', 'DisplayName', sprintf('SS Group Dist (Period %d)', bgp_period));
             hold off;
             xlabel(sprintf('Age Group Index (1 to %d)', cS.aD_new)); ylabel('Percent of Total Population (%)');
             title('V2 Mod: Initial vs Steady-State GROUP Population Distribution');
             legend('Location', 'best'); xticks(group_indices); grid on; drawnow;
             fprintf('Plotted initial vs steady-state/final GROUP population distribution.\n');
        end


        % =====================================================================
        % == Economic Functions (MODIFIED for Mixed Time Units) ==============
        % =====================================================================

        function cS = ParameterValues_HuggettStyle()
            % Combines V2 population dynamics setup with Huggett parameters/structure
            % Uses ANNUAL parameters for VFI/prices, GROUP parameters for pop dynamics

            %% Original Demographic Parameters (Annual - Copied from Huggett)
            cS.age1_orig      = 20;     % Start physical age
            cS.ageLast_orig   = 98;     % Max physical age
            cS.ageRetire_orig = 65;     % Retirement physical age
            cS.popGrowth_orig = 0.012;  % Annual population growth rate

            cS.aD_orig        = cS.ageLast_orig - cS.age1_orig + 1; % Total model years (79)
            cS.aR_idx_orig    = cS.ageRetire_orig - cS.age1_orig + 1; % Index of first retirement year (46)
            cS.aW_orig        = cS.aR_idx_orig - 1; % Number of working years (45)
            cS.physAgeV_orig   = (cS.age1_orig : cS.ageLast_orig)';

            % Annual Survival Probabilities (s_orig) & Death Rates (d_orig) - Copied from Huggett
            cS.d_orig         = [0.00159, 0.00169, 0.00174, 0.00172, 0.00165, 0.00156, 0.00149, 0.00145, 0.00145, 0.00149,...
                            0.00156, 0.00163, 0.00171, 0.00181, 0.00193, 0.00207, 0.00225, 0.00246, 0.00270, 0.00299,...
                            0.00332, 0.00368, 0.00409, 0.00454, 0.00504, 0.00558, 0.00617, 0.00686, 0.00766, 0.00865,...
                            0.00955, 0.01058, 0.01162, 0.01264, 0.01368, 0.01475, 0.01593, 0.01730, 0.01891, 0.02074,...
                            0.02271, 0.02476, 0.02690, 0.02912, 0.03143, 0.03389, 0.03652, 0.03930, 0.04225, 0.04538,...
                            0.04871, 0.05230, 0.05623, 0.06060, 0.06542, 0.07066, 0.07636, 0.08271, 0.08986, 0.09788,...
                            0.10732, 0.11799, 0.12895, 0.13920, 0.14861, 0.16039, 0.17303, 0.18665, 0.20194, 0.21877,...
                            0.23601, 0.25289, 0.26973, 0.28612, 0.30128, 0.31416, 0.32915, 0.34450, 0.36018];
            if length(cS.d_orig) ~= cS.aD_orig, error('Length mismatch d_orig'); end
            cS.s_orig         = 1 - cS.d_orig;

            % Static Annual Age Mass Distribution (For reference, not used in equilibrium calc)
            ageMassV_orig_static = ones(1, cS.aD_orig);
            for i = 2 : cS.aD_orig
                ageMassV_orig_static(i) = cS.s_orig(i-1) * ageMassV_orig_static(i-1) / (1 + cS.popGrowth_orig);
            end
            cS.ageMassV_orig_static = ageMassV_orig_static ./ sum(ageMassV_orig_static);


            %% Parameters for 5-Year Age Groups
            cS.yearStep      = 5;
            cS.aD_new        = ceil(cS.aD_orig / cS.yearStep); % Number of 5-year groups (16)
            cS.aR_new        = ceil(cS.aW_orig / cS.yearStep); % Index of last working 5-year group (9)

            % Map new group index 'a' to original annual indices
            cS.physAgeMap = cell(cS.aD_new, 1);
            for a = 1:cS.aD_new
                startIdx = (a-1)*cS.yearStep + 1;
                endIdx = min(a*cS.yearStep, cS.aD_orig);
                cS.physAgeMap{a} = startIdx:endIdx;
            end

            % Physical age corresponding to the start of each 5-year group
            cS.physAgeV_new = zeros(cS.aD_new, 1);
             for a = 1:cS.aD_new
                cS.physAgeV_new(a) = cS.physAgeV_orig(cS.physAgeMap{a}(1));
            end

            % Calculate 1-year survival probability for transition between groups
            % (Using survival prob s_orig of the last year in the current group)
            cS.s_1yr_transitionV = zeros(cS.aD_new, 1);
            for a = 1:(cS.aD_new - 1)
                 lastYearIdx = cS.physAgeMap{a}(end);
                 if lastYearIdx < cS.aD_orig
                     cS.s_1yr_transitionV(a) = cS.s_orig(lastYearIdx);
                 else
                     cS.s_1yr_transitionV(a) = 0;
                 end
            end
            cS.s_1yr_transitionV(cS.aD_new) = 0; % Cannot survive past the last group

            % Calculate GROUP population mass (static, for reference/ L calc)
            cS.ageMassV_new_static = zeros(1, cS.aD_new);
            for a = 1:cS.aD_new
                cS.ageMassV_new_static(a) = sum(cS.ageMassV_orig_static(cS.physAgeMap{a}));
            end
            cS.ageMassV_new_static = cS.ageMassV_new_static ./ sum(cS.ageMassV_new_static);
            cS.retireMass_new = sum(cS.ageMassV_new_static(cS.aR_new + 1 : cS.aD_new));


            %% Population Dynamics Parameters (Group Level - From V2)
            cS.initial_pop = [76.20905535, 86.45596319, 113.8702459, 98.60198303, 86.64117824, 102.7890433, 112.0217783, 99.04620047, ...
                              64.05142331, 66.93157492, 44.16815149, 25.40848066, 14.97325553, 6.872421945, 1.743059943, 0.216184341];
            if length(cS.initial_pop) ~= cS.aD_new, error('initial_pop length mismatch'); end

            % Multi-year survival probabilities (PopDyn - beta_surv from V2)
            beta_surv_pop = [0.998, 0.996, 0.994, 0.992, 0.988, 0.984, 0.980, ...
                             0.976, 0.970, 0.960, 0.945, 0.920, ...
                             0.880, 0.800, 0.680]; % Length aD_new - 1
            if length(beta_surv_pop) ~= cS.aD_new - 1, error('Length of beta_surv_pop incorrect'); end
            cS.survivalProbV_popdyn = [beta_surv_pop(:)', 0]; % Length aD_new

            % BGP Detection Parameters & Max Periods (From V2)
            cS.bgp_tolerance = 0.001;
            cS.bgp_window = 5;
            cS.max_periods = 50;

            %% Household Parameters (Annual)
            cS.sigma      = 1.5;   % Utility curvature
            cS.beta       = 1.011; % Annual discount factor (USED IN VFI)
            cS.cFloor     = 0.05;  % Consumption floor
            cS.nSim       = 5e4;   % Number of simulated individuals

            %% Technology Parameters (Annual)
            cS.A          = 0.895944;
            cS.alpha      = 0.36;
            cS.ddk        = 0.06;    % Annual depreciation rate (USED FOR PRICES)

            %% Social Security Parameters (Annual)
            cS.theta      = 0.2;     % Payroll tax rate for SS

            %% Labor Endowment Parameters (Annual Process)
            cS.leSigma1      = 0.38 ^ 0.5;
            cS.leShockStd    = 0.045 .^ 0.5;
            cS.lePersistence = 0.96;
            cS.leWidth       = 4;
            cS.nw            = 18;

            %% Grid Parameters
            cS.tgKY          = 3;
            cS.tgWage        = (1-cS.alpha)*cS.A*((cS.tgKY/cS.A)^(cS.alpha/(1-cS.alpha)));
            cS.nk            = 50;
            cS.kMin          = 0;
            cS.kMax          = 100 * cS.tgWage; % Adjusted max capital? Keep for now.
            kGridV           = linspace(cS.kMin, cS.kMax, cS.nk);
             % Make grid denser at the bottom (copied from V2)
            power = 1.5; % Adjust power for more/less density at bottom
            kGridV = cS.kMin + (cS.kMax - cS.kMin) * (linspace(0, 1, cS.nk).^power);
            kGridV(1)=cS.kMin; % Ensure minimum is exact
            cS.kGridV        = kGridV(:);

            %% Age Efficiency Profiles (Annual and Grouped)
            ageEffV_orig_temp = zeros(100, 1); % Use standard Huggett profile definition
             ageEffV_orig_temp(20:72) = [linspace(0.3, 1.5, 36-20+1), 1.5 .* ones(1, 47-37+1), ...
                               linspace(1.5, 0.2, 65-48+1), linspace(0.18, 0, 72-66+1)];
            cS.ageEffV_orig = ageEffV_orig_temp(cS.age1_orig : cS.ageLast_orig); % Vector for model annual ages 1..aD_orig
            if length(cS.ageEffV_orig) ~= cS.aD_orig, error('ageEffV_orig length mismatch'); end

            % Grouped Age Efficiency Profile (Average over 5 years)
            cS.ageEffV_new = zeros(cS.aD_new, 1);
            for a = 1:cS.aD_new
                cS.ageEffV_new(a) = mean(cS.ageEffV_orig(cS.physAgeMap{a}));
            end
        end

        function [Y, R, w, b, MPL_out] = HHPrices_Huggett(K, L, cS, paramS_in) % <--- 修改：增加 paramS_in 输入，增加 MPL_out 输出
            % Calculates ANNUAL prices R, w, b
            % MODIFIED to use dynamically simulated retiree counts from paramS_in.Z_ss_counts

            if K <= 0, K=1e-6; warning('K non-positive, setting to small value.');  end % Added K=1e-6
            if L <= 0, L=1e-6; warning('L non-positive, setting to small value.');  end % Added L=1e-6
            if nargin < 4 || ~isfield(paramS_in, 'Z_ss_counts') || isempty(paramS_in.Z_ss_counts)
                warning('Z_ss_counts not provided to HHPrices_Huggett or is empty. Using static retireMass_new for b calculation.');
                % 使用静态计算作为备用（如果paramS_in未正确传递或Z_ss_counts为空）
                total_retirees = cS.retireMass_new * (sum(cS.initial_pop)/100); % 近似总退休人口数 (如果initial_pop归一到100)
                                                                            % 这是一个粗略的备用，理想情况下应该报错或停止
                if cS.retireMass_new < 1e-9, total_retirees = 0; end
            else
                total_retirees = sum(paramS_in.Z_ss_counts(cS.aR_new+1:end)); % 使用动态稳态的退休人口数量
            end


            Y   = cS.A * (K^cS.alpha) * (L^(1-cS.alpha));
            MPK = cS.alpha * Y / K; % 修正MPK计算
            MPL = (1-cS.alpha) * Y / L; % 修正MPL计算
            MPL_out = MPL; % 输出税前MPL

            tau = 0; % 假设资本所得税为0（简化模型中）

            R   = 1 + (MPK - cS.ddk)*(1 - tau);
            w   = MPL*(1 - tau - cS.theta); % 净工资计算仍使用 cS.theta (固定贡献率)
            R   = max(1.0 + 1e-6, R);
            w   = max(0, w);

            % Social security benefit (Annualized)
            b = 0;
            if total_retirees > 1e-9 && L > 0 % 使用传入的 total_retirees
                total_wage_bill = MPL * L; % 税前总工资
                b = cS.theta * total_wage_bill / total_retirees; % 人均福利
            elseif L > 0 && total_retirees <= 1e-9
                % 如果没有退休人员，但有劳动者在缴费，理论上b可以为0，税收充公或形成盈余
                % 但当前模型是PAYG，所以b应该为0
                b = 0;
            end
            b = max(0, b);
        end

        function [logGridV, trProbM, prob1V] = EarningProcess_olgm(cS)
            % Generate grids and transition matrix for ANNUAL AR(1) process.
             [logGridV_raw, trProbM_calc] = utils.tauchen(cS.nw, cS.lePersistence, cS.leShockStd, 0, cS.leWidth);
             gridMin = logGridV_raw(1); gridMax = logGridV_raw(end);
             normGridMin = gridMin - 2*cS.leShockStd; normGridMax = gridMax + 2*cS.leShockStd;
             prob1V_calc = [];
              try
                  [prob1V_calc, ~, ~] = utils.norm_grid(logGridV_raw, normGridMin, normGridMax, 0, cS.leSigma1);
                  prob1V_calc = prob1V_calc(:);
              catch ME
                  warning('norm_grid failed: %s. Using uniform.', ME.message);
                  prob1V_calc = ones(cS.nw, 1)/cS.nw;
              end
             logGridV_calc = logGridV_raw(:);
             % Original Huggett code shifted grid, retain that.
             logGridV_calc = logGridV_calc - logGridV_calc(1) - 1;

             % Normalize transition matrix
             if any(abs(sum(trProbM_calc, 2) - 1) > 1e-6)
                  row_sums = sum(trProbM_calc, 2); row_sums(row_sums <= 1e-9) = 1;
                  trProbM_calc = bsxfun(@rdivide, trProbM_calc, row_sums);
             end
             % Normalize initial probabilities
             if abs(sum(prob1V_calc) - 1) > 1e-6
                 prob1V_calc = prob1V_calc ./ sum(prob1V_calc);
             end

             % *** Assign calculated values to output arguments ***
             logGridV = logGridV_calc;
             trProbM = trProbM_calc;
             prob1V = prob1V_calc;
             % *****************************************************
        end

        function [y, trProbM] = tauchen(N, pRho, pSigma, pMu, n_std) % (Unchanged)
             std_y = sqrt(pSigma^2 / (1 - pRho^2)); if abs(1-pRho)<1e-9, std_y = pSigma*100; end
             y_max = n_std*std_y; y_min = -y_max; y = linspace(y_min, y_max, N);
             d = 0; if N > 1, d = y(2) - y(1); end; trProbM = zeros(N, N);
             for iR = 1:N, for iC = 1:N, mn = pRho*y(iR); if iC==1, trProbM(iR,iC)=normcdf((y(1)-mn+d/2)/pSigma); elseif iC==N, trProbM(iR,iC)=1-normcdf((y(N)-mn-d/2)/pSigma); else trProbM(iR,iC)=normcdf((y(iC)-mn+d/2)/pSigma)-normcdf((y(iC)-mn-d/2)/pSigma); end; end; end
             rs=sum(trProbM,2); rs(rs<=1e-9)=1; trProbM=bsxfun(@rdivide,trProbM,rs);
             um=pMu/(1-pRho); if abs(1-pRho)<1e-9&&pMu~=0,um=sign(pMu)*inf; end; if isfinite(um), y=y+um; end
         end

        function [massV, lbV, ubV] = norm_grid(xV, xMin, xMax, mu, sig) % (Unchanged)
             n=length(xV); xV=xV(:)'; if n>1&&any(xV(2:n)<xV(1:n-1)), error('xV non-decreasing'); end
             lbV=[]; ubV=[]; if n>1, xM=0.5*(xV(1:n-1)+xV(2:n)); lbV=[xMin,xM]; ubV=[xM,xMax]; elseif n==1, lbV=xMin; ubV=xMax; else massV=[]; return; end
             cdfV=normcdf([lbV,ubV(end)],mu,sig); massV=diff(cdfV);
              if any(massV < -1e-9), warning('Negative mass detected in norm_grid.'); massV(massV<0)=0; end
             ts=sum(massV); if ts>1e-9, massV=massV/ts; else if n>0, massV=ones(1,n)/n; else massV=[]; end; end
             massV=massV(:); lbV=lbV(:); ubV=ubV(:); if n>1&&any(ubV<lbV), error('ubV < lbV'); end
        end

        function [HHlaborM_group, L] = LaborSupply_Huggett(eIdxM, cS, paramS, Z_ss_norm_group)
            % Calculates aggregate labor supply L.
            % Uses GROUPED efficiency and GROUPED steady-state mass Z_ss_norm_group.
            % Needs annual endowment simulation eIdxM to average labor within groups.

            nSim = size(eIdxM, 1);
            HHlaborM_group = zeros(nSim, cS.aD_new); % Average annual labor supply PER GROUP

            % Map annual age index to 5-year group index
            ageToGroupMap = zeros(cS.aD_orig, 1);
            for a_new_map = 1:cS.aD_new
                ageToGroupMap(cS.physAgeMap{a_new_map}) = a_new_map;
            end
            leGridV_col = paramS.leGridV(:);

            % Calculate labor for each simulant/annual age first
            HHlaborM_annual_temp = zeros(nSim, cS.aD_orig);
            for a_orig = 1 : cS.aD_orig
               if a_orig < cS.aR_idx_orig % Only working ages contribute
                   a_new = ageToGroupMap(a_orig); % Find the group
                   % Use GROUPED efficiency profile cS.ageEffV_new
                   HHlaborM_annual_temp(:, a_orig) = cS.ageEffV_new(a_new) .* leGridV_col(eIdxM(:, a_orig));
               end
            end

            % Average the annual labor supplies within each group
            for a_new = 1:cS.aD_new
                annual_indices = cS.physAgeMap{a_new};
                if ~isempty(annual_indices)
                    % Average labor supply over the years within the group for each simulant
                    HHlaborM_group(:, a_new) = mean(HHlaborM_annual_temp(:, annual_indices), 2);
                end
            end

            % Aggregate labor supply L using mean GROUP labor and GROUP masses
            L = mean(HHlaborM_group, 1) * Z_ss_norm_group(:);
            L = max(0, L); % Ensure non-negative
        end

        function eIdxM = MarkovChainSimulation(nSim, T, prob0V, trProbM, rvInM) % (Unchanged)
             ns=length(prob0V); if size(trProbM,1)~=ns||size(trProbM,2)~=ns,error('TrMat size mismatch');end; if abs(sum(prob0V)-1)>1e-5,prob0V=prob0V./sum(prob0V);end; if any(abs(sum(trProbM,2)-1)>1e-5), row_sums=sum(trProbM,2); row_sums(row_sums<=1e-9)=1; trProbM=bsxfun(@rdivide,trProbM,row_sums); end; if size(rvInM,1)~=nSim||size(rvInM,2)~=T,error('RV size mismatch');end
             cP0=cumsum(prob0V(:)'); cPT=cumsum(trProbM,2); cPT(:,ns)=1.0; eIdxM=zeros(nSim,T,'uint16'); eIdxM(:,1)=1+sum(bsxfun(@gt,rvInM(:,1),cP0),2); % Use uint16
             for t=1:(T-1), cSI=eIdxM(:,t); vPI=(cSI>=1)&(cSI<=ns); if ~all(vPI),cSI(~vPI)=1;eIdxM(:,t)=cSI;end; cPt=cPT(cSI,:); eIdxM(:,t+1)=1+sum(bsxfun(@gt,rvInM(:,t+1),cPt),2); end
        end

        function eIdxM = LaborEndowSimulation_olgm(cS, paramS)
            % Simulates ANNUAL endowment indices using ANNUAL transitions.
            rng(433); rvInM=rand([cS.nSim,cS.aD_orig]); % Use aD_orig
            eIdxM = utils.MarkovChainSimulation(cS.nSim, cS.aD_orig, paramS.leProb1V, paramS.leTrProbM, rvInM);
        end

                % =====================================================================
        % == HHSimulation_olgm (Corrected ageToGroupMap logic) ==============
        % =====================================================================
        function [kHistM, cHistM] = HHSimulation_olgm(kPolM, cPolM, eIdxM, cS)
            % Simulate ANNUAL capital and consumption paths using GROUP policies.

            validateattributes(kPolM, {'double'}, {'finite', 'nonnan', 'nonempty', 'real', ...
                               'size', [cS.nk, cS.nw, cS.aD_new]});
             validateattributes(cPolM, {'double'}, {'finite', 'nonnan', 'nonempty', 'real', ...
                               'size', [cS.nk, cS.nw, cS.aD_new]});
            validateattributes(eIdxM, {'uint16', 'double'}, {'nonempty', 'integer', 'positive', ...
                                '<=', cS.nw, 'size', [cS.nSim, cS.aD_orig]});

            % --- START CORRECTION ---
            % Map annual age index to 5-year group index correctly
            ageToGroupMap = zeros(cS.aD_orig, 1); % Initialize as numerical vector
            for a_new = 1:cS.aD_new
                annual_indices = cS.physAgeMap{a_new}; % Get annual indices for this group
                if ~isempty(annual_indices)
                    % Use standard parentheses () for indexing the numerical vector
                    ageToGroupMap(annual_indices) = a_new; % Assign group index a_new
                end
            end
            % --- END CORRECTION ---

            nSim   = size(eIdxM, 1);
            kHistM = zeros(nSim, cS.aD_orig + 1); % Store k at start of each year
            cHistM = zeros(nSim, cS.aD_orig);     % Store c during each year
            kGridV_col = cS.kGridV(:);

            % Pre-create interpolants for faster access
            kPolInterp = cell(cS.nw, cS.aD_new);
            cPolInterp = cell(cS.nw, cS.aD_new);
             for a_new=1:cS.aD_new
                 for ie=1:cS.nw
                      % Consider using 'pchip' for smoother interpolation if needed
                      kPolInterp{ie, a_new} = griddedInterpolant(kGridV_col, kPolM(:, ie, a_new), 'linear', 'none'); % Use 'none' to catch extrapolation issues
                      cPolInterp{ie, a_new} = griddedInterpolant(kGridV_col, cPolM(:, ie, a_new), 'linear', 'none');
                 end
             end

            for a_orig = 1 : cS.aD_orig
                a_new = ageToGroupMap(a_orig); % Get the relevant 5-year group index
                kNowV = kHistM(:, a_orig); % Capital at start of year a_orig

                kNextV = zeros(nSim, 1);
                cValV = zeros(nSim, 1);

                for ie = 1 : cS.nw
                    simIdxV = find(eIdxM(:, a_orig) == ie);
                    if ~isempty(simIdxV)
                         kNow_ie = kNowV(simIdxV);
                         % Clamp to grid for interpolation (required if 'none' extrapolation)
                         kNow_ie = max(kGridV_col(1), min(kGridV_col(end), kNow_ie));

                         try
                             kNextV(simIdxV) = kPolInterp{ie, a_new}(kNow_ie);
                             cValV(simIdxV) = cPolInterp{ie, a_new}(kNow_ie);
                             % Check for NaNs if 'none' extrapolation used
                             nan_k_idx = isnan(kNextV(simIdxV));
                             nan_c_idx = isnan(cValV(simIdxV));
                             if any(nan_k_idx) || any(nan_c_idx)
                                % Find original indices corresponding to NaNs for reporting
                                orig_nan_idx = simIdxV(nan_k_idx | nan_c_idx);
                                warning('NaN from interpolation detected in HHSimulation: age_orig=%d, group=%d, shock=%d, %d cases. Using nearest.', a_orig, a_new, ie, length(orig_nan_idx));

                                % Use nearest neighbor for NaNs
                                if any(nan_k_idx)
                                    kNextV(simIdxV(nan_k_idx)) = interp1(kGridV_col, kPolM(:, ie, a_new), kNow_ie(nan_k_idx), 'nearest', cS.kMin); % Extrapolate to kMin
                                end
                                if any(nan_c_idx)
                                     cValV(simIdxV(nan_c_idx)) = interp1(kGridV_col, cPolM(:, ie, a_new), kNow_ie(nan_c_idx), 'nearest', cS.cFloor); % Extrapolate to cFloor
                                end
                             end
                         catch ME_interp
                             warning('Interpolation error in HHSimulation: age_orig=%d, group=%d, shock=%d. %s. Using nearest.', a_orig, a_new, ie, ME_interp.message);
                             kNextV(simIdxV) = interp1(kGridV_col, kPolM(:, ie, a_new), kNow_ie, 'nearest', cS.kMin);
                             cValV(simIdxV) = interp1(kGridV_col, cPolM(:, ie, a_new), kNow_ie, 'nearest', cS.cFloor);
                         end
                    end
                end

                 kHistM(:, a_orig + 1) = max(cS.kMin, min(cS.kMax, kNextV));
                 cHistM(:, a_orig) = max(cS.cFloor, cValV);
                 % Ensure finiteness
                 kHistM(~isfinite(kHistM(:, a_orig + 1)), a_orig + 1) = cS.kMin;
                 cHistM(~isfinite(cHistM(:, a_orig)), a_orig) = cS.cFloor;
            end

            kHistM = kHistM(:, 1:cS.aD_orig); % Capital at START of years 1 to aD_orig

            validateattributes(kHistM, {'double'}, {'finite', 'nonnan', 'nonempty', 'real', ...
                               '>=', cS.kMin - 1e-6, '<=', cS.kMax + 1e-6, ...
                               'size', [nSim, cS.aD_orig]});
            validateattributes(cHistM, {'double'}, {'finite', 'nonnan', 'nonempty', 'real', ...
                                '>=', cS.cFloor - 1e-6, 'size', [nSim, cS.aD_orig]});
        end

        function [muM, utilM] = CES_utility(cM, sig, cS) % Pass cS for cFloor
             % Utility function. Handles cFloor. (Adapted from V2)
             if ~isscalar(sig)||sig<=0,error('Sigma must be positive scalar.');end;
             min_c = cS.cFloor;
             vld = (cM >= min_c);
             c_c = max(min_c, cM);
             utilM = -Inf(size(cM)); muM = Inf(size(cM));

             if abs(sig-1)<1e-6
                 utilM(vld) = log(c_c(vld));
                 muM(vld) = 1 ./ c_c(vld);
             else
                 utilM(vld) = (c_c(vld).^(1-sig)) ./ (1-sig);
                 muM(vld) = c_c(vld).^(-sig);
             end
             utilM(~vld) = -1e10 - (min_c - cM(~vld))*1e10; % Penalty
             muM(~vld) = c_c(~vld).^(-sig) + 1e10; % High MU
        end

        function incomeM = HHIncome_Huggett(kV, R, w, T, b, a_new, paramS, cS)
            % Calculates ANNUAL income for a given 5-year age group a_new
            % Uses GROUPED efficiency profile.

            nk = length(kV);
            leGridV_col = paramS.leGridV(:);

            if a_new <= cS.aR_new % Working age group
               nonCapIncomeV = cS.ageEffV_new(a_new) .* w .* leGridV_col + T; % [nw, 1]
            else % Retired age group
               nonCapIncomeV = b .* ones(cS.nw, 1) + T; % [nw, 1]
            end

            incomeM = R * kV(:) * ones(1, cS.nw) + ones(nk, 1) * nonCapIncomeV(:)'; % [nk, nw]

            if any(~isfinite(incomeM(:)))
                warning('Non-finite income detected age group %d. Clamping.', a_new);
                incomeM(~isfinite(incomeM)) = 1e-6; % Or median?
            end
        end

        function [cPolM, kPolM, valueM] = HHSolution_VFI_Huggett(R, w, T, bV_new, paramS, cS)
            % Solve HH problem via VFI over 5-year age groups, using ANNUAL prices/beta/survival.
            cPolM  = zeros(cS.nk, cS.nw, cS.aD_new);
            kPolM  = zeros(cS.nk, cS.nw, cS.aD_new);
            valueM = zeros(cS.nk, cS.nw, cS.aD_new); % Value at start of group 'a'

            for a = cS.aD_new : -1 : 1
               vPrimeM = []; % Value function for the *next* group (a+1)
               if a < cS.aD_new, vPrimeM = valueM(:,:,a+1); end

               [cPolM(:,:,a), kPolM(:,:,a), valueM(:,:,a)] = ...
                  utils.HHSolutionByAge_VFI_Huggett(a, vPrimeM, R, w, T, bV_new(a), paramS, cS);
            end

            validateattributes(cPolM, {'double'}, {'finite', 'nonnan', 'real','>=',cS.cFloor-1e-6,'size',[cS.nk,cS.nw,cS.aD_new]});
            validateattributes(kPolM, {'double'}, {'finite', 'nonnan', 'real','>=',cS.kMin-1e-6,'<=',cS.kMax+1e-6,'size',[cS.nk,cS.nw,cS.aD_new]});
            validateattributes(valueM, {'double'}, {'finite', 'nonnan', 'real','size',[cS.nk,cS.nw,cS.aD_new]});
        end

        function [cPolM, kPolM, valueM] = HHSolutionByAge_VFI_Huggett(a_new, vPrime_keM, R, w, T, b, paramS, cS)
            % Solves the HH problem for one specific 5-year age group 'a_new'.
            % Uses ANNUAL prices/beta/survival. Expectation uses ANNUAL transition matrix.

            if a_new < cS.aD_new && isempty(vPrime_keM), error('vPrime empty a<aD'); end
            if a_new == cS.aD_new && ~isempty(vPrime_keM), error('vPrime not empty a=aD'); end
            if a_new < cS.aD_new, validateattributes(vPrime_keM,{'double'},{'finite','nonnan','real','size',[cS.nk,cS.nw]}); end

            yM = utils.HHIncome_Huggett(cS.kGridV, R, w, T, b, a_new, paramS, cS); % Annual income

            fminOpt=optimset('fminbnd'); fminOpt.TolX=1e-6; fminOpt.Display='none';

            cPolM = zeros(cS.nk, cS.nw); kPolM = zeros(cS.nk, cS.nw); valueM = -Inf(cS.nk, cS.nw);

            if a_new == cS.aD_new % Last group
               cPolM = max(cS.cFloor, yM);
               kPolM = cS.kMin*ones(cS.nk,cS.nw);
               [~, valueM] = utils.CES_utility(cPolM, cS.sigma, cS);
               valueM(~isfinite(valueM)) = -1e12;
            else % Not the last group
               % Expected value function E[V'(k', e')] for start of group a+1
               % Uses ANNUAL transition matrix paramS.leTrProbM (as per Huggett example)
               ExValuePrime_k_given_e = zeros(cS.nk, cS.nw); % E[V'(k',e')|e]
                for ikPrime = 1:cS.nk
                   % Sum over e': V'(k',e') * Pr(e'|e)
                   ExValuePrime_k_given_e(ikPrime, :) = vPrime_keM(ikPrime, :) * paramS.leTrProbM';
                end
                ExValuePrime_k_given_e(~isfinite(ExValuePrime_k_given_e))=-1e10;

               vPInt=cell(1,cS.nw);
               for ie=1:cS.nw
                   try % Linear interp faster, PCHIP smoother
                      vPInt{ie} = griddedInterpolant(cS.kGridV, ExValuePrime_k_given_e(:, ie), 'linear', 'linear');
                   catch ME, error('E[V] interp fail a=%d,ie=%d:%s',a_new,ie,ME.message); end
               end

               results_c = zeros(cS.nk, cS.nw); results_k = zeros(cS.nk, cS.nw); results_v = zeros(cS.nk, cS.nw);

               parfor ik = 1:cS.nk % Use parfor if available
               % for ik = 1:cS.nk % Use for debugging
                   c_row=zeros(1,cS.nw); k_row=zeros(1,cS.nw); v_row=zeros(1,cS.nw);
                   for ie = 1 : cS.nw
                        vPofK = vPInt{ie}; % Interpolant E[V'(k')|e]
                        [c_row(ie), k_row(ie), v_row(ie)] = ...
                              utils.HHSolutionByOneState_VFI_Huggett(a_new, yM(ik,ie), vPofK, fminOpt, cS, paramS);
                   end
                   results_c(ik,:) = c_row; results_k(ik,:) = k_row; results_v(ik,:) = v_row;
               end
               cPolM = results_c; kPolM = results_k; valueM = results_v;
            end
        end

        function [c, kPrime, ValueFunc] = HHSolutionByOneState_VFI_Huggett(a_new, y, vPofK_int, fminOpt, cS, paramS)
            % Solves the optimization problem for a single state (k, e) at age group a_new.
            % Uses ANNUAL beta and 1-year transition survival probability.

            kPMin = cS.kMin;
            kPMax = y - cS.cFloor;
            ValueFunc = -Inf;

            if kPMax <= kPMin % Corner solution: Consume floor, save minimum
                kPrime = kPMin;
                c = max(cS.cFloor, y - kPrime);
                 [~, u_c] = utils.CES_utility(c, cS.sigma, cS);
                if a_new < cS.aD_new
                    EV_c = -Inf;
                     try
                        kP_eval = max(cS.kGridV(1), min(cS.kGridV(end), kPrime));
                        EV_c = vPofK_int(kP_eval);
                     catch ME_interp
                         warning('Interpolation error corner age %d: %s. Using boundary.', a_new, ME_interp.message);
                         if kPrime < cS.kGridV(1), EV_c = vPofK_int(cS.kGridV(1)); else EV_c = vPofK_int(cS.kGridV(end)); end
                     end
                     if ~isfinite(EV_c), EV_c = -1e12; end
                     s_transition = cS.s_1yr_transitionV(a_new); % Survival prob for the transition year
                     ValueFunc = u_c + cS.beta * s_transition * EV_c; % Use ANNUAL beta
                else % Last period
                    ValueFunc = u_c;
                end
                if ~isfinite(ValueFunc), ValueFunc = -1e12; end
            else % Interior solution possible
                kPMin_opt = max(cS.kMin, kPMin);
                kPMax_opt = max(kPMin_opt + 1e-9, min(cS.kMax, kPMax)); % Ensure feasible range within [kMin, kMax]

                 if kPMin_opt >= kPMax_opt % Check again after adjustments
                     kPrime_opt = kPMin_opt;
                     [negV, ~] = BellmanInner(kPrime_opt, a_new, y, vPofK_int, cS, paramS);
                     warning('Optimization bounds invalid age %d, y %.2f. Setting k=kMin.', a_new, y);
                 else
                    obj = @(kP) negBellmanObjective(kP, a_new, y, vPofK_int, cS, paramS);
                    [kPrime_opt, negV, eflag] = fminbnd(obj, kPMin_opt, kPMax_opt, fminOpt);
                     % Check boundary solutions carefully
                     if eflag <= 0 || abs(kPrime_opt - kPMin_opt)<1e-7 || abs(kPrime_opt - kPMax_opt)<1e-7
                         [nV_min, ~] = BellmanInner(kPMin_opt, a_new, y, vPofK_int, cS, paramS);
                         [nV_max, ~] = BellmanInner(kPMax_opt, a_new, y, vPofK_int, cS, paramS);
                          if nV_min <= nV_max + 1e-9 % Prefer lower k' if values are close
                              kPrime_opt = kPMin_opt; negV = nV_min;
                          else
                              kPrime_opt = kPMax_opt; negV = nV_max;
                          end
                          if eflag <= 0, warning('fminbnd issue (eflag=%d) age %d, y %.2f. Using boundary k=%.2f', eflag, a_new, y, kPrime_opt); end
                     end
                 end

                 kPrime = kPrime_opt;
                 ValueFunc = -negV;
                 c = max(cS.cFloor, y - kPrime);
                 % Recalculate ValueFunc if consumption was floored? BellmanInner should handle this.
                 if ~isfinite(ValueFunc), ValueFunc = -1e12; end
            end

            kPrime = max(cS.kMin, min(cS.kMax, kPrime)); % Final bounds check
            c = max(cS.cFloor, y - kPrime);
            if ~isfinite(kPrime), kPrime = cS.kMin; end
            if ~isfinite(c), c = cS.cFloor; end
            if ~isfinite(ValueFunc), ValueFunc = -1e12; end

            validateattributes(kPrime,{'double'},{'finite','nonnan','scalar','real','>=',cS.kMin-1e-6,'<=',cS.kMax+1e-6});
            validateattributes(c,{'double'},{'finite','nonnan','scalar','real','>=',cS.cFloor-1e-6});
            validateattributes(ValueFunc,{'double'},{'finite','nonnan','scalar','real'});

             % --- Nested Helper Functions ---
             function negV = negBellmanObjective(kP, ag_obj, y_obj, evInt_obj, cS_in_obj, pS_in_obj)
                  [negV, ~] = BellmanInner(kP, ag_obj, y_obj, evInt_obj, cS_in_obj, pS_in_obj);
             end

             function [negVal, Val] = BellmanInner(kP_inner, ag_inner, y_inner, evInt_inner, cS_in, pS_in)
                  % Calculates the value V(k,e,a) for a given choice k' = kP_inner
                  cons = max(cS_in.cFloor, y_inner - kP_inner);
                   [~, util] = utils.CES_utility(cons, cS_in.sigma, cS_in);
                   if ~isfinite(util), util = -1e12; end
                   Val = util;

                   if ag_inner < cS_in.aD_new % If not the last period
                       evF = -Inf;
                       try
                           kP_eval = max(cS_in.kGridV(1), min(cS_in.kGridV(end), kP_inner));
                           evF = evInt_inner(kP_eval);
                       catch ME_interp_inner
                           warning('Interpolation error BellmanInner age %d: %s. Using boundary.', ag_inner, ME_interp_inner.message);
                            if kP_inner < cS_in.kGridV(1), evF = evInt_inner(cS_in.kGridV(1)); else evF = evInt_inner(cS_in.kGridV(end)); end
                       end
                       if ~isfinite(evF), evF = -1e12; end

                       s_transition = cS_in.s_1yr_transitionV(ag_inner); % Use 1-yr transition survival
                       Val = util + cS_in.beta * s_transition * evF; % Use ANNUAL beta
                   end

                   if ~isfinite(Val), negVal = 1e12; Val = -1e12; else negVal = -Val; end
             end
        end % End HHSolutionByOneState

    end % End methods (Static)
end % End classdef

% --- END OF FILE utils.m ---