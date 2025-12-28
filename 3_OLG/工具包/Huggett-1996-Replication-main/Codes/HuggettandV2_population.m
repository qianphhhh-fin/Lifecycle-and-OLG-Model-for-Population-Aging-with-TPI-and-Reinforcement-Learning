% --- START OF FILE compare_population_distributions.m ---

% Script to calculate and compare the steady-state population distributions
% from the original Huggett model and the main_olg_v2 model.

clc;
clear;
close all;

fprintf('--- Comparing Huggett vs. V2 Population Distributions ---\n\n');

%% 1. Calculate Huggett Model Distribution

fprintf('Calculating Huggett distribution...\n');
try
    % Load Huggett parameters and calculate its ageMassV
    cS_huggett = HuggettUtils.ParameterValues_Fixed();
    ageMassV_huggett = cS_huggett.ageMassV(:); % Ensure column vector
    physAgeV_huggett = cS_huggett.physAgeV(:); % Physical ages (20-98)

    fprintf('Huggett distribution calculated (aD = %d, sum = %.6f).\n', ...
            length(ageMassV_huggett), sum(ageMassV_huggett));
catch ME_huggett
    error('Failed to calculate Huggett distribution. Ensure HuggettUtils.m is in the path. Error: %s', ME_huggett.message);
end

%% 2. Calculate V2 Model Distribution

fprintf('\nCalculating V2 distribution...\n');
% try
    % Load V2 parameters (this now includes aligned survival probs if you implemented it)
    cS_v2 = main_olg_v2_utils.ParameterValues_HuggettStyle();

    % Simulate V2 population dynamics to get the steady-state distribution
    % We need the 5-year survival probabilities for the simulation
    if ~isfield(cS_v2, 'survivalProbV_5yr')
         warning(['cS_v2 structure does not contain survivalProbV_5yr field needed for population dynamics.' ...
                  ' Attempting calculation based on annual equivalent. This might be inaccurate if alignment was not done.']);
         % Attempt reverse calculation (less accurate) if needed, or default
         if isfield(cS_v2,'survivalProbV_annual_equiv')
            cS_v2.survivalProbV_5yr = cS_v2.survivalProbV_annual_equiv .^ 5;
            cS_v2.survivalProbV_5yr(cS_v2.aD) = 0;
         else
             error('Cannot run v2 population dynamics without survivalProbV_5yr or survivalProbV_annual_equiv in cS_v2.');
         end
    end

    popS_v2 = main_olg_v2_utils.initPopulation(cS_v2);
    popS_v2 = main_olg_v2_utils.populationDynamics(popS_v2, cS_v2);
    [Z_ss_v2, ~, bgp_reached_v2, bgp_period_v2] = main_olg_v2_utils.detectSteadyStatePopulation(popS_v2, cS_v2);

    % Normalize to get the distribution
    Z_ss_norm_v2 = Z_ss_v2 / sum(Z_ss_v2);
    Z_ss_norm_v2 = Z_ss_norm_v2(:); % Ensure column vector

    % Approximate age groups for plotting v2 results
    age_group_centers_v2 = cS_v2.age1 + 2.5 + (0:(cS_v2.aD-1))*5; % Midpoint age for each 5-yr group

    fprintf('V2 distribution calculated (aD = %d, sum = %.6f, BGP reached: %d at period %d).\n', ...
            length(Z_ss_norm_v2), sum(Z_ss_norm_v2), bgp_reached_v2, bgp_period_v2);

% catch ME_v2
%     error('Failed to calculate V2 distribution. Ensure main_olg_v2_utils.m is in the path and parameters are correct. Error: %s', ME_v2.message);
% end

%% 3. Plot Comparison

fprintf('\nPlotting comparison...\n');

figure('Name', 'Population Distribution Comparison');

% Plot Huggett distribution (line plot is better for annual data)
plot(physAgeV_huggett, ageMassV_huggett * 100, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Huggett (Annual)');
hold on;

% Plot V2 distribution (bar chart is better for grouped data)
% We plot the bar at the center of the age group
bar(age_group_centers_v2, Z_ss_norm_v2 * 100, 0.5, 'r', 'DisplayName', 'V2 (5-Year Groups)'); % 0.5 controls bar width

hold off;
xlabel('Age / Age Group Center');
ylabel('Percent of Total Population (%)');
title('Comparison of Steady-State Population Distributions');
legend('Location', 'northeast');
grid on;
xlim([min(physAgeV_huggett(1), age_group_centers_v2(1))-5, max(physAgeV_huggett(end), age_group_centers_v2(end))+5]); % Adjust x-axis limits

fprintf('\nComparison complete. Check the generated plot.\n');

% --- END OF FILE compare_population_distributions.m ---