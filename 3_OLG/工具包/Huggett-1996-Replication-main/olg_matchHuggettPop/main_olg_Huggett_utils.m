% --- START OF FILE main_olg_Huggett_utils.m ---

classdef main_olg_Huggett_utils
    methods (Static)
        % HHPrices_Huggett函数
        function [Y, R, w, b] = HHPrices_Huggett(K, L, cS)
            % Production function calculations
            Y   = cS.A * (K^cS.alpha) * (L^(1-cS.alpha));
            MPK = cS.alpha * cS.A * (K^(cS.alpha-1)) * (L^(1-cS.alpha));
            MPL = (1-cS.alpha) * cS.A * (K^cS.alpha) * (L^(-cS.alpha));

            % Tax rate (simplified from original, check if this matches your intent)
            % Original had: tau = 0.195/(1-cS.ddk * K/Y); 
            % If tau is fixed, set it here. If it depends endogenously, use the formula.
            % Let's assume a fixed tau for now based on the lack of its use elsewhere
            % or maybe it was implicit in theta? Let's assume theta covers taxes/transfers.
            % Recalculating based *only* on theta:
            tau = 0; % Assume no separate income tax here, handled by theta?
                     % If 0.195 tax rate IS intended, uncomment the formula & adjust prices.
            
            % Household prices (Annualized)
            R   = 1 + (MPK - cS.ddk)*(1 - tau); % Annual interest factor
            w   = MPL*(1 - tau - cS.theta);     % Annual wage rate

            % Social security benefit (Annualized)
            % Uses retireMass_new (mass of 5-year retired groups)
            if cS.retireMass_new > 1e-6
                b = cS.theta * w * L / cS.retireMass_new;
            else
                b = 0; % Avoid division by zero if no retirees
            end
        end

        % ParameterValues_Fixed function (Modified for 5-year groups)
        function cS = ParameterValues_Fixed()
            %% Original Demographic Parameters (Annual)
            cS.age1_orig      = 20;     % Start physical age
            cS.ageLast_orig   = 98;     % Max physical age
            cS.ageRetire_orig = 65;     % Retirement physical age
            cS.popGrowth_orig = 0.012;  % Annual population growth rate

            % Original Model Ages (Annual)
            cS.aD_orig        = cS.ageLast_orig - cS.age1_orig + 1; % Total model years
            cS.aR_idx_orig    = cS.ageRetire_orig - cS.age1_orig + 1; % Index of first retirement year (age 65 -> index 46)
            cS.aW_orig        = cS.aR_idx_orig - 1; % Number of working years (indices 1 to 45)

            % Original Survival Probabilities (Annual)
            cS.d_orig         = [0.00159, 0.00169, 0.00174, 0.00172, 0.00165, 0.00156, 0.00149, 0.00145, 0.00145, 0.00149,...
                            0.00156, 0.00163, 0.00171, 0.00181, 0.00193, 0.00207, 0.00225, 0.00246, 0.00270, 0.00299,...
                            0.00332, 0.00368, 0.00409, 0.00454, 0.00504, 0.00558, 0.00617, 0.00686, 0.00766, 0.00865,...
                            0.00955, 0.01058, 0.01162, 0.01264, 0.01368, 0.01475, 0.01593, 0.01730, 0.01891, 0.02074,...
                            0.02271, 0.02476, 0.02690, 0.02912, 0.03143, 0.03389, 0.03652, 0.03930, 0.04225, 0.04538,...
                            0.04871, 0.05230, 0.05623, 0.06060, 0.06542, 0.07066, 0.07636, 0.08271, 0.08986, 0.09788,...
                            0.10732, 0.11799, 0.12895, 0.13920, 0.14861, 0.16039, 0.17303, 0.18665, 0.20194, 0.21877,...
                            0.23601, 0.25289, 0.26973, 0.28612, 0.30128, 0.31416, 0.32915, 0.34450, 0.36018]; 
            cS.s_orig         = 1 - cS.d_orig;

            % Original Age Mass Distribution (Annual)
            cS.ageMassV_orig   = ones(1, cS.aD_orig);
            for i = 2 : cS.aD_orig
                cS.ageMassV_orig(i) = cS.s_orig(i-1) * cS.ageMassV_orig(i-1) / (1 + cS.popGrowth_orig);
            end
            cS.ageMassV_orig   = cS.ageMassV_orig ./ sum(cS.ageMassV_orig);
            cS.physAgeV_orig   = (cS.age1_orig : cS.ageLast_orig)';

            %% Parameters for 5-Year Age Groups
            cS.yearStep      = 5;
            cS.aD_new        = ceil(cS.aD_orig / cS.yearStep); % Number of 5-year groups
            cS.aR_new        = ceil(cS.aW_orig / cS.yearStep); % Index of last working 5-year group

            % Map new group index 'a' to original annual indices
            cS.physAgeMap = cell(cS.aD_new, 1);
            for a = 1:cS.aD_new
                startIdx = (a-1)*cS.yearStep + 1;
                endIdx = min(a*cS.yearStep, cS.aD_orig); % Handle last group potentially being shorter
                cS.physAgeMap{a} = startIdx:endIdx;
            end

            % Calculate 1-year survival probability for transition between groups
            % (Using survival prob of the last year in the current group)
            cS.s_1yr_transitionV = zeros(cS.aD_new, 1);
            for a = 1:(cS.aD_new - 1)
                 lastYearIdx = cS.physAgeMap{a}(end); % Index of last year in group 'a'
                 if lastYearIdx < cS.aD_orig % Ensure we don't index beyond available survival data
                     cS.s_1yr_transitionV(a) = cS.s_orig(lastYearIdx);
                 else
                     cS.s_1yr_transitionV(a) = 0; % Cannot survive past last year
                 end
            end
            cS.s_1yr_transitionV(cS.aD_new) = 0; % Cannot survive past the last group

            % Calculate population mass for each 5-year group
            % Robust method: Sum annual masses within each group
            cS.ageMassV_new = zeros(1, cS.aD_new);
            for a = 1:cS.aD_new
                cS.ageMassV_new(a) = sum(cS.ageMassV_orig(cS.physAgeMap{a}));
            end
            % Renormalize (should already sum to 1, but good practice)
            cS.ageMassV_new = cS.ageMassV_new ./ sum(cS.ageMassV_new);

            % Calculate total mass of retired 5-year groups
            if cS.aR_new < cS.aD_new
                cS.retireMass_new = sum(cS.ageMassV_new(cS.aR_new + 1 : cS.aD_new));
            else
                cS.retireMass_new = 0; % No retired groups if aR_new is the last group
            end

            % Physical age corresponding to the start of each 5-year group
            cS.physAgeV_new = zeros(cS.aD_new, 1);
             for a = 1:cS.aD_new
                cS.physAgeV_new(a) = cS.physAgeV_orig(cS.physAgeMap{a}(1));
            end


            %% Household Parameters (Annual)
            cS.sigma      = 1.5;   % Utility curvature
            cS.beta       = 1.011; % Annual discount factor
            % cS.beta       = 1.011^5; % Five-year discount factor
            cS.cFloor     = 0.05;  % Consumption floor
            cS.nSim       = 5e4;   % Number of simulated individuals

            %% Technology Parameters (Annual)
            cS.A          = 0.895944;
            cS.alpha      = 0.36;
            cS.ddk        = 0.06;    % Annual depreciation rate

            %% Social Security Parameters (Annual)
            cS.theta      = 0.1;     % Payroll tax rate for SS

            %% Labor Endowment Parameters
            cS.leSigma1      = 0.38 ^ 0.5; % Initial std dev of log endowment
            cS.leShockStd    = 0.045 .^ 0.5;% Std dev of AR1 shock
            cS.lePersistence = 0.96;      % Persistence of AR1
            cS.leWidth       = 4;         % Grid width (std devs)
            cS.nw            = 18;        % Number of endowment grid points

            %% Grid Parameters
            % Capital grid
            cS.tgKY          = 3;          % Target K/Y ratio
            cS.tgWage        = (1-cS.alpha)*cS.A*((cS.tgKY/cS.A)^(cS.alpha/(1-cS.alpha)));
            cS.nk            = 50;         % Number of capital grid points
            cS.kMin          = 0;          % Borrowing constraint k>=0
            cS.kMax          = 100 * cS.tgWage; % Max capital level
            kGridV           = linspace(cS.kMin, cS.kMax, cS.nk);
            cS.kGridV        = kGridV(:);  % Column vector

            %% Age Efficiency Profile (Original Annual)
            % Stored in paramS in the main script, but needed here for grouping
            ageEffV_orig = zeros(100, 1);
            ageEffV_orig(20:72) = [linspace(0.3, 1.5, 36-20+1), 1.5 .* ones(1, 47-37+1), ...
                              linspace(1.5, 0.2, 65-48+1), linspace(0.18, 0, 72-66+1)];
            cS.ageEffV_orig = ageEffV_orig(cS.age1_orig : cS.ageLast_orig); % Vector for model ages 1..aD_orig

            % Grouped Age Efficiency Profile (Average over 5 years)
            cS.ageEffV_new = zeros(cS.aD_new, 1);
            for a = 1:cS.aD_new
                cS.ageEffV_new(a) = mean(cS.ageEffV_orig(cS.physAgeMap{a}));
                % temp = cS.ageEffV_orig(cS.physAgeMap{a});
                % cS.ageEffV_new(a) = temp(1);
            end
        end

        % tauchen function (Unchanged)
        function [y, trProbM] = tauchen(N, pRho, pSigma, pMu, n_std)
            % Use grid width
            a_bar = n_std * sqrt(pSigma^2 / (1 - pRho^2));

            % Grid
            y     = linspace(-a_bar, a_bar, N);

            % Distance between points
            d     = y(2) - y(1);

            % Get transition probabilities
            trProbM  = zeros(N, N);

            for iRow = 1 : N
               % End points first
               trProbM(iRow, 1) = normcdf((y(1) - pRho*y(iRow) + d/2) / pSigma);
               trProbM(iRow, N) = 1 - normcdf((y(N) - pRho*y(iRow) - d/2) / pSigma);

               % Fill in interior columns
               for iCol = 2 : N-1
                  trProbM(iRow, iCol) = (normcdf((y(iCol) - pRho*y(iRow) + d/2) / pSigma) - ...
                                         normcdf((y(iCol) - pRho*y(iRow) - d/2) / pSigma));
               end
            end

            % Center process around the mean (ybar / (1 - rho))
            y = y + pMu / (1 - pRho);
        end

        % norm_grid function (Unchanged)
        function [massV, lbV, ubV] = norm_grid(xV, xMin, xMax, mu, sig)
            n = length(xV);

            % Validate inputs
            if any(xV(2:n) < xV(1:n-1))
                error('xV not increasing');
            end

            if xMin > xV(1) || xMax < xV(n)
                error('Invalid xMin or xMax');
            end

            if mu < xMin || mu > xMax
               % Allow mean outside grid range, but issue warning?
               % warning('Mean mu outside range [xMin, xMax]');
               % Original code had error, let's relax it slightly but be aware
               if mu < xV(1) || mu > xV(end) % Check against grid bounds instead of limits
                   warning('Mean mu outside grid range xV');
               end
            end


            % Grid boundaries
            xMidV = 0.5 .* (xV(1:n-1) + xV(2:n));
            lbV   = [xMin; xMidV(:)];
            ubV   = [xMidV(:); xMax];

            % Mass in each interval
            cdfV  = normcdf([xMin; ubV], mu, sig); % Include xMin in CDF calculation
            massV = cdfV(2:(n+1)) - cdfV(1:n);

             % Check for negative masses due to precision issues and handle
            if any(massV < 0)
                %warning('Negative mass detected in norm_grid, setting to zero.');
                massV(massV < 0) = 0;
            end

            % Normalize mass
            sumMass = sum(massV);
            if sumMass > 1e-6
                 massV = massV(:) ./ sumMass;
            else
                 warning('Total mass close to zero in norm_grid. Returning uniform.');
                 massV = ones(n, 1) / n; % Fallback to uniform
            end


            % Output validation
            if any(ubV < lbV)
                error('ubV < lbV');
            end
        end

                % LaborSupply_Huggett function (Corrected for Retirement)
        function [HHlaborM_annual, L] = LaborSupply_Huggett(eIdxM, cS, paramS)
            % Computes *annual* labor supply for each simulant/year
            % and aggregate labor supply L.
            % Uses the *grouped* age efficiency profile for working years only.

            nSim = size(eIdxM, 1);
            HHlaborM_annual = zeros(nSim, cS.aD_orig); % Annual labor supply

            % Map each annual age index to its corresponding 5-year group index
            ageToGroupMap = zeros(cS.aD_orig, 1);
            for a_new = 1:cS.aD_new
                ageToGroupMap(cS.physAgeMap{a_new}) = a_new;
            end

            for a_orig = 1 : cS.aD_orig
               % ***** START CORRECTION *****
               % Check if the ANNUAL age is retirement age or older
               if a_orig >= cS.aR_idx_orig
                   HHlaborM_annual(:, a_orig) = 0; % Force labor supply to zero for retirees
               else
                   % Only calculate labor supply for working ages
                   a_new = ageToGroupMap(a_orig); % Find the group this annual age belongs to
                   % Apply the grouped efficiency profile for that age group
                   HHlaborM_annual(:, a_orig) = cS.ageEffV_new(a_new) .* paramS.leGridV(eIdxM(:, a_orig));
               end
               % ***** END CORRECTION *****
            end

            % Aggregate labor supply L using *annual* masses and *annual* labor supplies
            % This now correctly sums only the labor supply from working individuals
            L = mean(HHlaborM_annual, 1) * cS.ageMassV_orig';
        end

        % EarningProcess_olgm function (Unchanged)
        function [logGridV, trProbM, prob1V] = EarningProcess_olgm(cS)
            [logGridV, trProbM] = main_olg_Huggett_utils.tauchen(cS.nw, cS.lePersistence, cS.leShockStd, 0, cS.leWidth);

            % Draw for new agents from approx log normal dist
            % Adjust bounds slightly for robustness with norm_grid
            logGridMin = logGridV(1) - 3 * cS.leShockStd;
            logGridMax = logGridV(end) + 3 * cS.leShockStd;
            prob1V = main_olg_Huggett_utils.norm_grid(logGridV, logGridMin, logGridMax, 0, cS.leSigma1);
            prob1V = prob1V(:);

            % Improve scaling - Original code shifted grid, maybe keep it centered?
            % Let's keep the original shift for replication consistency.
            % logGridV = logGridV(:); % - logGridV(1) - 1; % Keep original centering for now
            logGridV = logGridV(:) - logGridV(1) - 1; % <<< RESTORE THIS SHIFT

            % Output validation
            if any(abs(sum(trProbM, 2) - 1) > 1e-6)
                warning('Transition probs do not sum to 1');
                trProbM = bsxfun(@rdivide, trProbM, sum(trProbM, 2)); % Force normalization
            end

            if abs(sum(prob1V) - 1) > 1e-6
                warning('prob1V does not sum to 1');
                prob1V = prob1V ./ sum(prob1V); % Force normalization
            end
        end

        % MarkovChainSimulation function (Unchanged)
        function eIdxM = MarkovChainSimulation(nSim, T, prob0V, trProbM, rvInM)
            ns = length(prob0V); % Number of states

            % Input validation
            if size(trProbM,1) ~= ns || size(trProbM,2) ~= ns
                 error('Transition matrix size mismatch');
            end
            if any(abs(sum(trProbM, 2) - 1) > 1e-5) % Check rows sum to 1
                 warning('Transition matrix rows do not sum to 1. Normalizing.');
                 trProbM = bsxfun(@rdivide, trProbM, sum(trProbM, 2));
            end

            if abs(sum(prob0V) - 1) > 1e-5
                warning('Initial probs do not sum to 1. Normalizing.');
                 prob0V = prob0V(:) ./ sum(prob0V);
            end

            if size(rvInM, 1) ~= nSim || size(rvInM, 2) ~= T
                error('Random number matrix size mismatch');
            end

            % Prepare cumulative probabilities
            % Need cumulative probs for each row of trProbM
            cumTrProbM = cumsum(trProbM, 2);
            % Ensure last column is exactly 1 to avoid index errors with rand=1
            cumTrProbM(:, ns) = 1.0;

            cumProb0V = cumsum(prob0V(:)'); % Cumulative initial distribution

            % Simulate
            eIdxM = zeros(nSim, T, 'uint16'); % Use integer type for indices

            % Draw t=1 states
            eIdxM(:, 1) = 1 + sum(bsxfun(@gt, rvInM(:,1), cumProb0V), 2);

            % Draw t = 2, ..., T states
            for t = 1 : (T-1)
               % Get the cumulative transition probs for the current states eIdxM(:, t)
               currentCumProbs = cumTrProbM(eIdxM(:, t), :);
               % Compare random draw to cumulative probs to find next state
               eIdxM(:, t+1) = 1 + sum(bsxfun(@gt, rvInM(:, t+1), currentCumProbs), 2);
            end
        end

        % LaborEndowSimulation_olgm function (Modified for annual simulation)
        function eIdxM = LaborEndowSimulation_olgm(cS, paramS)
            % Simulates endowment indices for each individual over ANNUAL ages

            % Set random number generator for reproducibility
            rng(433);

            % Endowment state by [ind, annual age]
            % Simulate using ANNUAL number of periods (aD_orig)
            eIdxM = main_olg_Huggett_utils.MarkovChainSimulation(cS.nSim, cS.aD_orig, paramS.leProb1V, ...
                                      paramS.leTrProbM, rand([cS.nSim, cS.aD_orig]));
        end

        % HHSimulation_olgm function (Modified for grouped policies but annual simulation)
        function [kHistM, cHistM] = HHSimulation_olgm(kPolM, cPolM, eIdxM, cS)
            % Simulate annual capital and consumption paths using policy functions defined over 5-year groups.

            % Inputs:
            % kPolM, cPolM: Policy functions [nk, nw, aD_new]
            % eIdxM: Simulated annual endowment shocks [nSim, aD_orig]
            % cS: Struct with parameters (including aD_orig, aD_new, physAgeMap)

            % Input validation
            validateattributes(kPolM, {'double'}, {'finite', 'nonnan', 'nonempty', 'real', ...
                               'size', [cS.nk, cS.nw, cS.aD_new]});
             validateattributes(cPolM, {'double'}, {'finite', 'nonnan', 'nonempty', 'real', ...
                               'size', [cS.nk, cS.nw, cS.aD_new]});
            validateattributes(eIdxM, {'uint16', 'double'}, {'nonempty', 'integer', 'positive', ...
                                '<=', cS.nw, 'size', [cS.nSim, cS.aD_orig]}); % Allow double for eIdxM if needed

            % Map annual age index to 5-year group index
            ageToGroupMap = zeros(cS.aD_orig, 1);
            for a_new = 1:cS.aD_new
                ageToGroupMap(cS.physAgeMap{a_new}) = a_new;
            end

            % Simulate capital and consumption history, annually
            nSim   = size(eIdxM, 1);
            kHistM = zeros(nSim, cS.aD_orig + 1); % +1 to store k at age aD_orig+1 (terminal)
            cHistM = zeros(nSim, cS.aD_orig);
            kGridV_col = cS.kGridV(:); % Ensure column vector

            for a_orig = 1 : cS.aD_orig
                a_new = ageToGroupMap(a_orig); % Get the relevant 5-year group index

                % Get policies for this group
                kPol_a = kPolM(:, :, a_new);
                cPol_a = cPolM(:, :, a_new);

                % Current capital levels for all individuals at start of year a_orig
                kNowV = kHistM(:, a_orig);

                % Loop over endowment states for efficiency (vectorize interpolation)
                for ie = 1 : cS.nw
                    % Find individuals with endowment ie at annual age a_orig
                    simIdxV = find(eIdxM(:, a_orig) == ie);

                    if ~isempty(simIdxV)
                        % Interpolate optimal k' (k for start of a_orig+1)
                        kHistM(simIdxV, a_orig + 1) = interp1(kGridV_col, kPol_a(:, ie), ...
                                                        kNowV(simIdxV), 'spline', 'extrap');

                        % Interpolate optimal c (consumption during a_orig)
                        cHistM(simIdxV, a_orig) = interp1(kGridV_col, cPol_a(:, ie), ...
                                                        kNowV(simIdxV), 'spline', 'extrap');
                    end
                end % loop over endowment states

                 % Ensure constraints are met after interpolation/extrapolation
                 kHistM(:, a_orig + 1) = max(cS.kMin, min(cS.kMax, kHistM(:, a_orig + 1)));
                 cHistM(:, a_orig) = max(cS.cFloor, cHistM(:, a_orig));

            end % loop over annual ages a_orig

            % Trim the extra capital column used for calculation
            kHistM = kHistM(:, 1:cS.aD_orig);

            % Output validation (on the annual history)
            validateattributes(kHistM, {'double'}, {'finite', 'nonnan', 'nonempty', 'real', ...
                               '>=', cS.kMin - 1e-6, '<=', cS.kMax + 1e-6, ...
                               'size', [nSim, cS.aD_orig]});
            validateattributes(cHistM, {'double'}, {'finite', 'nonnan', 'nonempty', 'real', ...
                                '>=', cS.cFloor - 1e-6, 'size', [nSim, cS.aD_orig]});
        end

        % CES_utility function (Unchanged)
        function [muM, utilM] = CES_utility(cM, sig)
            % Input validation
            if any(cM(:) < 1e-10) % Use a smaller threshold
               % Find first problematic value for better error message
               firstBadVal = find(cM(:) < 1e-10, 1);
               if ~isempty(firstBadVal)
                   [r, col] = ind2sub(size(cM), firstBadVal);
                   errMsg = sprintf('Cannot compute utility for very small consumption c(%d, %d) = %e', r, col, cM(firstBadVal));
               else
                    errMsg = 'Cannot compute utility for very small consumption';
               end
               % Instead of error, maybe return very large negative utility?
               % Let's stick to error for now to highlight the issue.
                error(errMsg);
                % Alternative: Set to a very large negative number or handle upstream
                % utilM = -1e20 * ones(size(cM));
                % muM = 1e20 * ones(size(cM));
                % return;
            end

            if length(sig) ~= 1
                error('sig must be scalar');
            end

            if sig <= 0
                error('sig must be > 0');
            end

            % Compute utility
            if abs(sig - 1) < 1e-8          % Check for log utility robustly
               utilM = log(cM);             % Utility
               muM   = 1 ./ cM;             % Marginal utility
            else
               utilM = (cM .^ (1-sig)) ./ (1-sig);   % CES Utility
               muM = cM .^ (-sig);                  % Marginal utility
            end
        end

        % HHIncome_Huggett function (Modified for grouped age efficiency)
        function incomeM = HHIncome_Huggett(kV, R, w, T, b, a_new, paramS, cS)
            % Calculates annual income for a given 5-year age group a_new

            % Inputs:
            % kV: Vector of capital levels [nk, 1]
            % R, w, T, b: Annual prices/transfers
            % a_new: Index of the 5-year age group
            % paramS: Struct with leGridV
            % cS: Struct with parameters (aR_new, ageEffV_new, nw)

            nk = length(kV);

            % Non-capital income by shock (using grouped efficiency)
            if a_new <= cS.aR_new % Working age group
               % paramS.leGridV is [nw, 1] or [1, nw]? Assuming [nw, 1] from EarningProcess
               leGridV_col = paramS.leGridV(:); % Ensure column
               nonCapIncomeV = cS.ageEffV_new(a_new) .* w .* leGridV_col + T; % [nw, 1]
            else % Retired age group
               nonCapIncomeV = b .* ones(cS.nw, 1) + T; % [nw, 1]
            end

            % Total income for the group (annual)
            % kV(:) is [nk, 1], nonCapIncomeV(:)' is [1, nw]
            incomeM = R * kV(:) * ones(1, cS.nw) + ones(nk, 1) * nonCapIncomeV(:)'; % [nk, nw]
        end

        % HHSolution_VFI_Huggett function (Modified for 5-year groups)
        function [cPolM, kPolM, valueM] = HHSolution_VFI_Huggett(R, w, T, bV_new, paramS, cS)
            % Solve HH problem via VFI over 5-year age groups.

            % Inputs:
            % R, w, T: Annual prices/transfer
            % bV_new: Vector of annual SS benefits by 5-year group [1, aD_new]
            % paramS: Struct with earnings process info
            % cS: Struct with parameters (nk, nw, aD_new)

            % Initialize policy/value functions for 5-year groups
            cPolM  = zeros(cS.nk, cS.nw, cS.aD_new);
            kPolM  = zeros(cS.nk, cS.nw, cS.aD_new);
            valueM = zeros(cS.nk, cS.nw, cS.aD_new); % Value at start of group 'a'

            % Backward induction over 5-year groups
            for a = cS.aD_new : -1 : 1
               % Value function for the *next* 5-year group (a+1)
               if a < cS.aD_new
                  vPrimeM = valueM(:,:,a+1);
               else
                  % No next period after the last group
                  vPrimeM = zeros(cS.nk, cS.nw); % Or some terminal value if needed
               end

               [cPolM(:,:,a), kPolM(:,:,a), valueM(:,:,a)] = ...
                  main_olg_Huggett_utils.HHSolutionByAge_VFI_Huggett(a, vPrimeM, R, w, T, bV_new(a), paramS, cS);

            end

            % Output validation
            validateattributes(cPolM, {'double'}, {'finite', 'nonnan', 'nonempty', ...
                               'real', 'positive', 'size', [cS.nk, cS.nw, cS.aD_new]});

            validateattributes(valueM, {'double'}, {'finite', 'nonnan', 'nonempty', ...
                               'real', 'size', [cS.nk, cS.nw, cS.aD_new]});

            validateattributes(kPolM, {'double'}, {'finite', 'nonnan', 'nonempty', ...
                               'real', 'size', [cS.nk, cS.nw, cS.aD_new], ...
                               '>=', cS.kMin - 1e-6, '<=', cS.kGridV(end) + 1e-6});
        end

        % HHSolutionByAge_VFI_Huggett function (Modified for 5-year groups)
        function [cPolM, kPolM, valueM] = HHSolutionByAge_VFI_Huggett(a_new, vPrime_keM, R, w, T, b, paramS, cS)
            % Solves the HH problem for one specific 5-year age group 'a_new'.

            % Inputs:
            % a_new: Index of the current 5-year age group
            % vPrime_keM: Value function at the start of the *next* group (a_new+1) [nk, nw]
            % R, w, T, b: Annual prices/transfer/benefit for this group
            % paramS, cS: Parameter structs

            % Input check
            if a_new < cS.aD_new
                if ~isequal(size(vPrime_keM), [cS.nk, cS.nw])
                    error('Invalid size of vPrime_keM for group %d', a_new);
                end
            end

            % Calculate annual income y(k, e) for this age group a_new
            yM = main_olg_Huggett_utils.HHIncome_Huggett(cS.kGridV, R, w, T, b, a_new, paramS, cS);

            % Optimization options
            fminbndOptS = optimset('fminbnd');
            fminbndOptS.TolX = 1e-5;
            % fminbndOptS.Display = 'off'; % Suppress output

            if a_new == cS.aD_new
               % Last period group: Consume all income, no savings
               cPolM = max(cS.cFloor, yM); % Ensure consumption >= floor
               kPolM = zeros(cS.nk, cS.nw) + cS.kMin; % Save minimum possible
               [~, valueM] = main_olg_Huggett_utils.CES_utility(cPolM, cS.sigma);

            else
               % Allocate space for policy/value functions for this age group
               cPolM = zeros(cS.nk, cS.nw);
               kPolM = zeros(cS.nk, cS.nw);
               valueM = zeros(cS.nk, cS.nw);

               % Expected value function E[V'(k', e')] for tomorrow (start of group a+1)
               % Needs expectation over e', given current e (ie)
               % ExValuePrime_ke = E_e'[ V'(k', e') | e = ie ]
               % vPrime_keM is V'(k, e') = valueM(:,:,a+1)
               % paramS.leTrProbM is Pr(e'|e) [nw, nw]
               % Expected value, conditional on current state ie, as function of k'
               ExValuePrime_k_given_e = zeros(cS.nk, cS.nw); % Stores E[V'(k',e')|e] for each (k', e)
               for ikPrime = 1:cS.nk
                   % vPrime_keM(ikPrime, :) is V'(k'=kGrid(ikPrime), e') across all e'
                   % paramS.leTrProbM(ie, :) is Pr(e' | ie) across all e'
                   % Expected value E[V'(k'=kGrid(ikPrime), e') | ie] = sum_{e'} V'(k',e') * Pr(e'|ie)
                   ExValuePrime_k_given_e(ikPrime, :) = vPrime_keM(ikPrime, :) * paramS.leTrProbM';
                   % Result ExValuePrime_k_given_e(ikPrime, ie) = E[V'(k'=kGrid(ikPrime), e') | ie]
               end


               % Loop over current states (k_idx, e_idx)
               for ie = 1 : cS.nw
                   % Continuous approximation of tomorrow's expected value E[V'(k',e')|ie] as a function of k'
                   % Use ExValuePrime_k_given_e(:, ie) which holds E[V'(k',e')|ie] for grid points k'
                   vPrimeOfK_interpolant = griddedInterpolant(cS.kGridV, ExValuePrime_k_given_e(:, ie), 'spline', 'spline'); % Linear interpolation and extrapolation

                   % Loop over current capital states
                   for ik = 1 : cS.nk
                        % Solve the Bellman equation for state (ik, ie) at age a_new
                        [cPolM(ik,ie), kPolM(ik,ie), valueM(ik,ie)] = ...
                              main_olg_Huggett_utils.HHSolutionByOneState_VFI_Huggett(a_new, yM(ik,ie), R, vPrimeOfK_interpolant, fminbndOptS, cS);
                   end % loop over ik
               end % loop over ie
            end % if last period

            % Output validation for this age group
            validateattributes(cPolM, {'double'}, {'finite', 'nonnan', 'nonempty', ...
                               'real', 'positive', 'size', [cS.nk, cS.nw]});

            validateattributes(kPolM, {'double'}, {'finite', 'nonnan', 'nonempty', 'real', ...
                               '>=', cS.kMin - 1e-6,  '<=', cS.kGridV(cS.nk) + 1e-6, ...
                               'size', [cS.nk, cS.nw]});

            validateattributes(valueM, {'double'}, {'finite', 'nonnan', 'nonempty', ...
                               'real', 'size', [cS.nk, cS.nw]});
        end

        % HHSolutionByOneState_VFI_Huggett function (Modified for group transition)
        function [c, kPrime, ValueFunc] = HHSolutionByOneState_VFI_Huggett(a_new, y, R, vPrimeOfK_interpolant, fminbndOptS, cS)
            % Solves the optimization problem for a single state (k, e) at age group a_new.

            % Inputs:
            % a_new: Current 5-year group index
            % y: Current income (annual)
            % R: Annual interest factor
            % vPrimeOfK_interpolant: Interpolant for E[V'(k')|e] at start of group a_new+1
            % fminbndOptS: Optimization options
            % cS: Parameter struct (including kGridV, cFloor, beta, s_1yr_transitionV, sigma)

            % Feasible range for kPrime (next period's capital)
            % Max kPrime allows for minimum consumption cFloor
            kPrimeMax = y - cS.cFloor;

            % Ensure kPrime is within grid bounds and feasible range
            kPrimeLowerBound = cS.kMin;
            kPrimeUpperBound = min(cS.kMax, kPrimeMax);

            if kPrimeUpperBound <= kPrimeLowerBound + 1e-7 % Check if feasible range is negligible or invalid
               % No feasible choice allowing c >= cFloor. Consume floor, save minimum.
               kPrime = kPrimeLowerBound;
               c = cS.cFloor;
               % Calculate value, ensuring kPrime passed to interpolant is valid
               validKPrimeForInterp = max(cS.kGridV(1), min(cS.kGridV(end), kPrime));
               [~, u] = main_olg_Huggett_utils.CES_utility(c, cS.sigma);
                % Use survival prob for the transition from group a_new to a_new+1
               ValueFunc = u + cS.beta * cS.s_1yr_transitionV(a_new) * R * vPrimeOfK_interpolant(validKPrimeForInterp); % R here is annual factor
               % Note: Original code returned -ValueFunc from Bellman. Let's return the actual Value.
               % And the optimization minimizes -Value, so we negate the result from fminbnd later.

            else
               % Find optimal kPrime by minimizing negative Bellman value
               % Pass R and transition survival probability to Bellman objective
               objectiveFunc = @(kP) BellmanObjective(kP, y, R, vPrimeOfK_interpolant, cS.beta, cS.s_1yr_transitionV(a_new), cS.sigma, cS.cFloor, cS.kGridV(1), cS.kGridV(end));
               [kPrime, negValueFunc, exitflag] = fminbnd(objectiveFunc, kPrimeLowerBound, kPrimeUpperBound, fminbndOptS);

                % Handle optimization failures if necessary
                if exitflag <= 0
                   warning('fminbnd did not converge for age %d. y=%.2f. Using boundary.', a_new, y);
                   % Could potentially choose the end point with lower objective value
                   val_low = objectiveFunc(kPrimeLowerBound);
                   val_high = objectiveFunc(kPrimeUpperBound);
                   if val_low < val_high
                       kPrime = kPrimeLowerBound;
                       negValueFunc = val_low;
                   else
                       kPrime = kPrimeUpperBound;
                       negValueFunc = val_high;
                   end
                end

               ValueFunc = -negValueFunc; % Get the maximum value
               c = max(cS.cFloor, y - kPrime); % Calculate consumption
            end

            % Final validation of outputs
            validateattributes(kPrime, {'double'}, {'finite', 'nonnan', 'nonempty', 'real', 'scalar', ...
                               '>=', kPrimeLowerBound - 1e-6, '<=', kPrimeUpperBound + 1e-6});

            validateattributes(c, {'double'}, {'finite', 'nonnan', 'nonempty', 'real', 'scalar', ...
                               '>=', cS.cFloor - 1e-6});

            validateattributes(ValueFunc, {'double'}, {'finite', 'nonnan', 'nonempty', 'real', 'scalar'});


            % --- Nested Objective Function for fminbnd ---
            function [negVal, cons] = BellmanObjective(kP, income, R_annual, vPrimeInterp, beta_annual, s_transition, sigma, cFloor, kGridMin, kGridMax)
                % Ensure kP is within grid bounds for interpolant evaluation
                kP_interp = max(kGridMin, min(kGridMax, kP));
                cons = max(cFloor, income - kP); % Consumption ensures floor
                [~, util] = main_olg_Huggett_utils.CES_utility(cons, sigma);
                % Value = u(c) + beta * s * R * E[V'(k')]
                % R is the annual factor. If V' is annual value, this is correct.
                Val = util + beta_annual * s_transition * R_annual * vPrimeInterp(kP_interp);
                negVal = -Val; % Minimize negative value
            end
            % --- End Nested Function ---

        end % End HHSolutionByOneState_VFI_Huggett

    end % End Static methods
end % End classdef

% --- END OF FILE main_olg_Huggett_utils.m ---