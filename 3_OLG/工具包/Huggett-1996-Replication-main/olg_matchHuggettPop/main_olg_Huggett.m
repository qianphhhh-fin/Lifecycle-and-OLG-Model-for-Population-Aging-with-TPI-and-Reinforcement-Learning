% --- START OF FILE main_olg_Huggett.m ---

% Huggett (1996) Replication with 5-Year Age Groups
% Multi-Period OLG Model with Earnings and Lifetime Uncertainty
% Based on Yanran Guo's code, Modified by AI Assistant
% Date: [Current Date]
% Uses main_olg_Huggett_utils class

%% Environment Setup
clc; clear; close all;
% addpath(genpath(pwd)); % If utils are in a subfolder

fprintf('Huggett (1996) OLG Model with 5-Year Age Groups\n');
fprintf('Initializing...\n');

%% Fixed Model Parameters (Now includes 5-year grouped parameters)
cS = main_olg_Huggett_utils.ParameterValues_Fixed();
fprintf('Parameters loaded. Number of 5-year groups: %d. Working groups: %d.\n', cS.aD_new, cS.aR_new);

%% Precalibrate Labor Endowment Process and Grid
% This part is independent of age grouping, defines the shock process
[paramS.leLogGridV, paramS.leTrProbM, paramS.leProb1V] = main_olg_Huggett_utils.EarningProcess_olgm(cS);
paramS.leGridV = exp(paramS.leLogGridV(:)); % Ensure column vector and take exp
fprintf('Labor endowment process calibrated (nw=%d).\n', cS.nw);

% Store original and grouped age efficiency profiles in paramS for clarity
paramS.ageEffV_orig = cS.ageEffV_orig;
paramS.ageEffV_new  = cS.ageEffV_new;


%% General Equilibrium Solution Steps
%{
Overall Strategy:
1. Precompute aggregate labor supply L (using grouped efficiency but annual simulation & mass).
2. Guess aggregate capital K and lump-sum transfers T.
3. Iterate:
    a. Compute prices (R, w, b) based on K, L. (b uses grouped retire mass)
    b. Solve HH problem (VFI) over 5-year groups using annual prices -> Policy functions cPolM, kPolM [nk, nw, aD_new].
    c. Simulate HH decisions *annually* using the group-based policies -> kHistM, cHistM [nSim, aD_orig].
    d. Compute model's aggregate capital K_model (from annual kHistM & annual mass).
    e. Compute model's aggregate bequests T_model (from annual kHistM, annual death rates & mass).
    f. Check deviation |K - K_model| + |T - T_model|. Update K, T. Repeat until convergence.
%}

%% Step 1: Precompute Aggregate Labor Supply L
fprintf('Simulating labor endowments for %d individuals over %d years...\n', cS.nSim, cS.aD_orig);
% Simulate annual endowment paths
eIdxM = main_olg_Huggett_utils.LaborEndowSimulation_olgm(cS, paramS); % [nSim, aD_orig]

fprintf('Calculating aggregate labor supply L...\n');
% ============= Calculate annual labor supply using grouped efficiency,
% then aggregate =================
% [~, L] = main_olg_Huggett_utils.LaborSupply_Huggett(eIdxM, cS, paramS);
% fprintf('Aggregate Labor Supply L = %.4f\n', L);

% --- START Calculate L using ANNUAL efficiency ---
HHlaborM_annual_L = zeros(cS.nSim, cS.aD_orig);
leGridV_col_L = paramS.leGridV(:); % Ensure column vector
for a_orig_L = 1 : cS.aD_orig
   if a_orig_L >= cS.aR_idx_orig
       HHlaborM_annual_L(:, a_orig_L) = 0; % Retirees have zero labor
   else
       % Use original ANNUAL efficiency profile stored in paramS
       HHlaborM_annual_L(:, a_orig_L) = paramS.ageEffV_orig(a_orig_L) .* leGridV_col_L(eIdxM(:, a_orig_L));
   end
end
L = mean(HHlaborM_annual_L, 1) * cS.ageMassV_orig'; % Aggregate using annual mass
fprintf('Aggregate Labor Supply L = %.4f\n', L); % This L should be closer to 9.2
% --- END Calculate L using ANNUAL efficiency ---


%% Step 2: Initial Guesses for K and T
% Use values from the original code as starting points, adjust if needed
KGuess = 6.3816; % Initial guess for aggregate capital
TGuess = 0.0806;  % Initial guess for aggregate accidental bequests (annual transfer T)

%% Step 3-6: Iteration to Find Equilibrium K* and T*

% Iteration parameters
maxIter  = 100;      % Maximum iterations
tolLevel = 1e-5;     % Convergence tolerance (slightly relaxed from 1e-6)
dampK    = 0.3;       % Damping factor for K update
dampT    = 0.3;       % Damping factor for T update
iter     = 0;        % Iteration counter
devNorm  = tolLevel + 1; % Ensure loop runs at least once

fprintf('\nStarting iteration to find equilibrium (K*, T*)...\n');
fprintf('Iter |   K_Guess  |   T_Guess  |   K_Model  |   T_Model  |  K_Dev   |  T_Dev   |   Norm   |\n');
fprintf('------------------------------------------------------------------------------------------\n');

while iter < maxIter && devNorm > tolLevel
    iter = iter + 1;

    % Step 3a: Compute Prices given KGuess, L
    [~, R, w, b_annual] = main_olg_Huggett_utils.HHPrices_Huggett(KGuess, L, cS);

    % Create vector of SS benefits by 5-year age group
    bV_new = zeros(1, cS.aD_new);
    if cS.aR_new < cS.aD_new % Check if there are retired groups
        bV_new(cS.aR_new + 1 : cS.aD_new) = b_annual;
    end

    % Step 3b: Solve HH Problem (VFI over 5-year groups)
    [cPolM, kPolM, ~] = main_olg_Huggett_utils.HHSolution_VFI_Huggett(R, w, TGuess, bV_new, paramS, cS);
    % cPolM, kPolM are size [nk, nw, aD_new]

    % Step 3c: Simulate HH Decisions Annually
    % Uses grouped policies kPolM, cPolM but annual simulation eIdxM
    [kHistM, ~] = main_olg_Huggett_utils.HHSimulation_olgm(kPolM, cPolM, eIdxM, cS);
    % kHistM is size [nSim, aD_orig]

    % Step 3d: Compute Model's Aggregate Capital K_model
    % Aggregate from annual simulation using annual population weights
    KModel = mean(kHistM, 1) * cS.ageMassV_orig'; % Weighted average over individuals and annual ages
    KModel = max(0.01, KModel); % Ensure K > 0

    % Step 3e: Compute Model's Aggregate Bequests T_model
    % Calculate expected capital holdings of those who die each year (annual)
    % Need capital at the *end* of the period (start of next) for bequests
    kprimeHistM = zeros(cS.nSim, cS.aD_orig);
     % Map annual age index to 5-year group index
     ageToGroupMap = zeros(cS.aD_orig, 1);
     for a_new_map = 1:cS.aD_new
         ageToGroupMap(cS.physAgeMap{a_new_map}) = a_new_map;
     end
     kGridV_col = cS.kGridV(:);

     for a_orig_beq = 1:cS.aD_orig
         a_new_beq = ageToGroupMap(a_orig_beq);
         kPol_a_beq = kPolM(:, :, a_new_beq); % Policy for this group
         kNowV_beq = kHistM(:, a_orig_beq);
         for ie_beq = 1:cS.nw
             simIdxV_beq = find(eIdxM(:, a_orig_beq) == ie_beq);
             if ~isempty(simIdxV_beq)
                  kprimeHistM(simIdxV_beq, a_orig_beq) = interp1(kGridV_col, kPol_a_beq(:, ie_beq), ...
                                                           kNowV_beq(simIdxV_beq), 'linear', 'extrap');
             end
         end
          kprimeHistM(:, a_orig_beq) = max(cS.kMin, min(cS.kMax, kprimeHistM(:, a_orig_beq)));
     end


    % Accidental bequests = Sum [ mass(a) * death_prob(a) * mean_k'(a) * R ] / (1+g)
    % Using original annual mass and death probabilities
    ageDeathMassV_orig = cS.ageMassV_orig .* cS.d_orig; % Mass of agents dying at end of each annual age
    % Mean k' * R for each annual age group
    meanKprimeR_orig = mean(kprimeHistM * R, 1); % Average across simulations for each age
    % Total value of bequests generated this period
    totalBequests = meanKprimeR_orig * ageDeathMassV_orig';

    % Distribute bequests lump-sum next period, adjusted for population growth
    % Division factor: (1+g) seems standard to make it per capita for the growing cohort of newborns.
    % Sticking to provided code: /(1-g)? Let's use (1+g) as it feels more standard. Revert if needed.
    TModel = totalBequests / (1 + cS.popGrowth_orig);
    TModel = max(0, TModel); % Ensure non-negative

    % Step 3f: Check Deviation and Update Guesses
    KDev = KGuess - KModel;
    TDev = TGuess - TModel;
    devV = [KDev; TDev];
    devNorm = norm(devV);

    fprintf('%4d | %10.4f | %10.4f | %10.4f | %10.4f | %9.4f | %9.4f | %9.5f |\n', ...
            iter, KGuess, TGuess, KModel, TModel, KDev, TDev, devNorm);

    % Check convergence
    if devNorm < tolLevel
        fprintf('\nEquilibrium Found!\n');
        break;
    end

    % Update guesses using damping
    KGuess = KGuess - dampK * KDev;
    TGuess = TGuess - dampT * TDev;

    % Ensure guesses remain non-negative
    KGuess = max(0.01, KGuess);
    TGuess = max(0, TGuess);

end % End equilibrium iteration loop

% Check if convergence failed
if iter == maxIter && devNorm >= tolLevel
    warning('Equilibrium search did not converge within %d iterations. Final norm: %e', maxIter, devNorm);
end

%% Final Equilibrium Results
K_eq = KGuess;
T_eq = TGuess;

% Compute final prices and policies at equilibrium K*, T*
[Y_eq, R_eq, w_eq, b_eq_annual] = main_olg_Huggett_utils.HHPrices_Huggett(K_eq, L, cS);
bV_eq_new = zeros(1, cS.aD_new);
if cS.aR_new < cS.aD_new
    bV_eq_new(cS.aR_new + 1 : cS.aD_new) = b_eq_annual;
end
[cPolM_eq, kPolM_eq, ~] = main_olg_Huggett_utils.HHSolution_VFI_Huggett(R_eq, w_eq, T_eq, bV_eq_new, paramS, cS);

fprintf('\n--- Final Equilibrium Results ---\n');
fprintf('Aggregate Capital (K*) = %.4f\n', K_eq);
fprintf('Aggregate Labor (L)    = %.4f\n', L);
fprintf('Aggregate Output (Y*)  = %.4f\n', Y_eq);
fprintf('Interest Rate (r*)     = %.4f (%.2f%% annual)\n', R_eq - 1, (R_eq-1)*100);
fprintf('Wage Rate (w*)         = %.4f\n', w_eq);
fprintf('SS Benefit (b*)        = %.4f\n', b_eq_annual);
fprintf('Bequests (T*)          = %.4f\n', T_eq);
fprintf('Capital/Output Ratio   = %.3f\n', K_eq / Y_eq);
fprintf('---------------------------------\n');


%% Figures (Plotting Policy Functions for a Representative Age Group)

% Choose a representative 5-year age group to plot
% e.g., Group 5 corresponds to ~ age 40-44 (Indices 21-25)
plot_a_new = 5;
if plot_a_new > cS.aD_new
    plot_a_new = cS.aD_new; % Ensure valid index
end
physAgeStart = cS.physAgeV_new(plot_a_new);
physAgeEnd = cS.physAgeV_orig(cS.physAgeMap{plot_a_new}(end));

fprintf('\nPlotting policy functions for age group %d (Approx. Physical Ages %d-%d)\n', ...
        plot_a_new, physAgeStart, physAgeEnd);

% Select representative endowment levels (low, medium, high)
ie_low = 1;
ie_med = round(cS.nw / 2);
ie_high = cS.nw;

% Figure 1: Optimal Consumption Policy c(k, e, a_new)
figure('Name', sprintf('Consumption Policy: Group %d (Age %d-%d)', plot_a_new, physAgeStart, physAgeEnd));
plot(cS.kGridV, cPolM_eq(:, ie_low, plot_a_new), '-k', 'LineWidth', 1.5);
hold on;
plot(cS.kGridV, cPolM_eq(:, ie_med, plot_a_new), '--b', 'LineWidth', 1.5);
plot(cS.kGridV, cPolM_eq(:, ie_high, plot_a_new), ':r', 'LineWidth', 1.5);
hold off;
xlabel('Capital (k)');
ylabel('Consumption (c)');
title(sprintf('Optimal Consumption c(k, e, a_{new}=%d)', plot_a_new));
legend(sprintf('Low Endowment (e=%d)', ie_low), ...
       sprintf('Medium Endowment (e=%d)', ie_med), ...
       sprintf('High Endowment (e=%d)', ie_high), ...
       'Location', 'northwest');
grid on;

% Figure 2: Optimal Savings Policy k'(k, e, a_new)
figure('Name', sprintf('Savings Policy: Group %d (Age %d-%d)', plot_a_new, physAgeStart, physAgeEnd));
plot(cS.kGridV, kPolM_eq(:, ie_low, plot_a_new), '-k', 'LineWidth', 1.5);
hold on;
plot(cS.kGridV, kPolM_eq(:, ie_med, plot_a_new), '--b', 'LineWidth', 1.5);
plot(cS.kGridV, kPolM_eq(:, ie_high, plot_a_new), ':r', 'LineWidth', 1.5);
% Add 45-degree line for reference (k'=k)
plot(cS.kGridV, cS.kGridV, '-g', 'LineWidth', 0.5);
hold off;
xlabel('Current Capital (k)');
ylabel('Next Period Capital (k'')');
title(sprintf('Optimal Savings k''(k, e, a_{new}=%d)', plot_a_new));
legend(sprintf('Low Endowment (e=%d)', ie_low), ...
       sprintf('Medium Endowment (e=%d)', ie_med), ...
       sprintf('High Endowment (e=%d)', ie_high), ...
       'k''=k Line', ...
       'Location', 'northwest');
grid on;
axis tight; % Adjust axis limits

fprintf('Simulation and plotting complete.\n');

% --- END OF FILE main_olg_Huggett.m ---