function cS = ParameterValues_Fixed
%% Documentation
%  This function assigns parameter values to replicate Huggett (1996) model.


%% ************************ Demographics ************************ 
% Physical age
cS.age1      = 20;     % Everyone in the model economy starts at age = 20
cS.ageLast   = 98;     % Max age
cS.ageRetire = 65;     % HHs retire at age = 65
cS.popGrowth = 0.012;  % Population growth rate

% Model age
cS.aD        = cS.ageLast - cS.age1 + 1;
cS.aR        = cS.ageRetire - cS.age1 + 1;

% Survival probabilities
% Based on Jordan, C., Life contingencies, 2nd ed. (Society of Actuaries).
% Huggett (1996) dates it as the 1975 print, the following link is to the 1991 imprint. 
% But since it is still the 2nd edition (which first appeared 1967)
% Seems to be the right numbers. (the pg 342 of the pdf, 346 of book, appears to be the closest thing)
% https://vdocuments.mx/download/life-contingencies-chester-wallace-jordanpdf
% Conditional probability of death: physical age 20 to 98
cS.d         = [0.00159, 0.00169, 0.00174, 0.00172, 0.00165, 0.00156, 0.00149, 0.00145, 0.00145, 0.00149,...
                0.00156, 0.00163, 0.00171, 0.00181, 0.00193, 0.00207, 0.00225, 0.00246, 0.00270, 0.00299,...
                0.00332, 0.00368, 0.00409, 0.00454, 0.00504, 0.00558, 0.00617, 0.00686, 0.00766, 0.00865,...
                0.00955, 0.01058, 0.01162, 0.01264, 0.01368, 0.01475, 0.01593, 0.01730, 0.01891, 0.02074,...
                0.02271, 0.02476, 0.02690, 0.02912, 0.03143, 0.03389, 0.03652, 0.03930, 0.04225, 0.04538,...
                0.04871, 0.05230, 0.05623, 0.06060, 0.06542, 0.07066, 0.07636, 0.08271, 0.08986, 0.09788,...
                0.10732, 0.11799, 0.12895, 0.13920, 0.14861, 0.16039, 0.17303, 0.18665, 0.20194, 0.21877,...
                0.23601, 0.25289, 0.26973, 0.28612, 0.30128, 0.31416, 0.32915, 0.34450, 0.36018]; 
% Conditional survival probabilities. Act as a discount rate.
cS.s         = 1 - cS.d; 

% Mass of households by age
cS.ageMassV   = ones(1, cS.aD);
for i = 2 : length(cS.ageMassV)
    cS.ageMassV(i) = cS.s(i-1) * cS.ageMassV(i-1) / (1 + cS.popGrowth);
end
cS.ageMassV   = cS.ageMassV./sum(cS.ageMassV);

% Mass of retired households
cS.retireMass = sum(cS.ageMassV(cS.aR + 1 : end));

% Physical age for each model age
cS.physAgeV   = (cS.age1 : cS.ageLast)';


%% ****************************** Household *******************************
cS.sigma      = 1.5;   % Curvature of utility function
cS.beta       = 1.011; % Discount factor
cS.cFloor     = 0.05;  % Consumption floor. Introduced for numerical reasons.
cS.nSim       = 5e4;   % Number of individuals to simulate 


%% ****************************** Technology ****************************** 
cS.A          = 0.895944;
cS.alpha      = 0.36;
cS.ddk        = 0.06;


%% *************************** Social Security ****************************
cS.theta      = 0.1;


%% *************************** Labor Endowment ****************************
cS.leSigma1      = 0.38 ^ 0.5;
cS.leShockStd    = 0.045 .^ 0.5;
cS.lePersistence = 0.96;     
cS.leWidth       = 4;     % Number of standard deviations the grid is wide
cS.nw            = 18;    % Size of labor endowment grids


%% ******************************* Grids **********************************
% Capital grid
cS.tgKY          = 3;          % Targeted capital / output ratio is 3.
cS.tgWage        = (1-cS.alpha)*cS.A*((cS.tgKY/cS.A)^(cS.alpha/(1-cS.alpha)));
cS.nk            = 50;
cS.kMin          = 0;          % Due to the borrowing constraint: k>=0
cS.kMax 	     = 100 * cS.tgWage;
kGridV           = linspace(cS.kMin, cS.kMax, cS.nk); 
cS.kGridV        = kGridV(:);  % kGridV is a 1xnk row vector now, make it a nkx1 vector


end