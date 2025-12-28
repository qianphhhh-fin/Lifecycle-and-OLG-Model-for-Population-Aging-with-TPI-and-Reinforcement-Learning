% Huggett (1996) Replication
% Multi-Period OLG Model with Earnings and Lifetime Uncertainty
% Yanran Guo
% 8/7/2018

%--------------------------------------------------------------------
%{
In this exercise, there are several points need to be noticed

(1). Grid for k
     Number of grids can be set arbitrarily
     kMin is set to be 0 due to borrowing constraint
     kMax is set by using function kgrid_olgm

(2). Labor endowment process
     This process is calibrated by using function EarningProcess_olgm
     - Use Tauchen (1986) method to produce finte state Markov approximation 
       of the AR(1) process
     - For the initial labor endowment, approximate a normal distribution
       on the grids provided by Tauchen method

(3). Age efficiency profile 
     When calibrating labor endowment process, function 'EarningProcess_olgm'
     generates pure labor endowment states and is net of the age efficiency
     profile. 
     So we need to generate a rough linear approximation of Huggett (1996), 
     by physical age levels

%}


%% Environment Setup
clc;                     % Clear screen
clear;                   % Clear memory
close all;               % Close open windows
addpath(genpath(pwd));   % Include subfolders of the current folder


%% Fixed Model Parameters
cS             = ParameterValues_Fixed;

% Precalibrate labor endowment process
[paramS.leLogGridV, paramS.leTrProbM, paramS.leProb1V] ...
               = EarningProcess_olgm(cS);
paramS.leGridV = exp(paramS.leLogGridV);

% Age efficiency profile
% A rough linear approximation of Huggett (1996), by physical age levels
ageEffV        = zeros(100, 1);
ageEffV(20:72) = [linspace(0.3, 1.5, 36-20+1), 1.5 .* ones(1, 47-37+1), ...
                  linspace(1.5, 0.2, 65-48+1), linspace(0.18, 0, 72-66+1)];
paramS.ageEffV = ageEffV(cS.age1 : cS.ageLast);


%% Solve HH Problem without Calibrating 
% Note that in this replication exercise, I solve HH problem with given
% prices. These prices (interest rate R, wage w, lump-sum transfer T, and 
% social security benefit b) can be calibrated later by having an extra
% outside loop. The calibration idea is
% Step 1. Precompute the aggregate labor supply, L. Since age efficiency is
%         fixed and there is no labor supply decision, L can be precalibrated 
%         without solving household problem
% Step 2. Guess total capital KGuess and lump-sum transfer of accidental bequests TGuess
% Step 3. Given KGuess and L, compute prices R, w, and b
% Step 4. Given prices R, w, b, and given TGuess, solve HH problem and 
%         get policy function for saving
% Step 5. Given solution to HH problem, compute the aggregate saving K and 
%         accidental bequests T
% Step 6. devV = [KGuess - K; TGuess - T]
%         If devV is above the tolerance level, update guess for K and T,
%         and go back to Step 2.


% Compute the agrregate labor supply
eIdxM        = LaborEndowSimulation_olgm(cS, paramS);
[~, L]       = LaborSupply_Huggett(eIdxM, cS, paramS);

% Assign values to K and T without calibrating these two
K            = 56.430264399364560;
T            = 0.902795718208179;

% Get the prices faced by household
[~, R, w, b] = HHPrices_Huggett(K, L, cS);
bV           = [zeros(1, cS.aR), ones(1, cS.aD - cS.aR) .* b];

% Backward Induction with Value Function Iteration
% HHSolution_VFI_olgm --> HHSolutionByAge_VFI_olgm --> HHSolutionByOneState_VFI_olgm
[cPolM, kPolM, valueM] ...
             = HHSolution_VFI_Huggett(R, w, T, bV, paramS, cS);


%% Figure
% Figure 1. 
% Policy function -- optimal consumption
figure;
plot(cS.kGridV, cPolM(:, 1, 23), 'k', 'LineWidth', 1);
hold on
plot(cS.kGridV, cPolM(:, 10, 23), 'k--', 'LineWidth', 1);
plot(cS.kGridV, cPolM(:, 18, 23), 'r', 'LineWidth', 1);
hold off
xlabel('wealth k');
ylabel('c(k,e,23)');
title('optimal consumption at model age 23')
legend('e=1', 'e=10', 'e=18', 'Location', 'northwest');

% Figure 2.
% Policy function -- optimal capital level
figure;
plot(cS.kGridV, kPolM(:, 1, 23), 'k', 'LineWidth', 1);
hold on
plot(cS.kGridV, kPolM(:, 10, 23), 'k--', 'LineWidth', 1);
plot(cS.kGridV, kPolM(:, 18, 23), 'r', 'LineWidth', 1);
hold off
xlabel('wealth k');
ylabel('kprime(k,e,23)');
title('optimal capital level at model age 23')
legend('e=1', 'e=10', 'e=18', 'Location', 'northwest');


%% Compute the Aggregate K and T from the Model and Compare with Guess Values
% Simulate capital histories from policy function indexed by [ind, age]
[kHistM, ~]   = HHSimulation_olgm(kPolM, cPolM, eIdxM, cS);

% Aggregate capital in the model economy
% Let the mass of households with age a be mu(a), then the aggregate capital is
% K = sum_(a=1 to aD) mu(a) x mean( kHistM(:,a) )
KModel        = mean(kHistM,1) * cS.ageMassV'; 
KModel        = max(0.01, KModel);

% Transfer and accidental bequests
kprimeHistM   = [kHistM(:,2:end), zeros(size(kHistM,1),1)];
ageDeathMassV = cS.ageMassV .* cS.d;
acci_bequest  = (mean(kprimeHistM * R, 1) * ageDeathMassV') / (1 - cS.popGrowth);

devV(1)       = K - KModel;
devV(2)       = T - acci_bequest;
fprintf('    Kdev: %5.4f   Tdev: %5.4f \n', devV(1), devV(2));

function [muM, utilM] = CES_utility(cM, sig)
%% Documentation
% CES utility function

%{
INPUTS:
(1). cM:    Consumption (matrix of any dimension)
(2). sig:   Sigma. Curvature parameter
(3). dbg:   Debugging parameter

OUTPUTS:
(1). muM:   Marginal utility
(2). utilM: Utility
%}


%%  Input Validation
% 1. Consumption cannot be too small
if any(cM(:) < 1e-8)
    error('Cannot compute utility for very small consumption');
end

% 2. Sigma must be a scalar
if length(sig) ~= 1
    error('sig must be scalar');
end

% 3. Sigma must be a positive number
if sig <= 0
    error('sig must be > 0');
end


%% Compute the Utility: log when sig=1; CES when sig>1
if sig == 1                            % Log utility
    utilM = log(cM);                    % Utility
    muM   = 1 ./ cM;                    % Marginal utility
else
    utilM = cM .^ (1-sig) ./ (1-sig);   % CES Utility
    muM = cM .^ (-sig);                 % Marginal utility
end


end

