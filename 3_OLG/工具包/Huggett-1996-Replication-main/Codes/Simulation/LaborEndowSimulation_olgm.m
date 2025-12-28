function eIdxM = LaborEndowSimulation_olgm(cS, paramS)
%% Documentation:
% This function simulates labor endowments for a cohort of household as it
% moves through the ages

% OUTPUT:
% eIdxM: labor endowment index by [ind, age]
%        It is a (nSim x age) matrix
%        nSim is the number of individuals we want to simulate
%        age  is the total ages in the model


% ******************************** Notice ********************************* 
% We simulate nSim individuals for each age
% eIdxM(999, 29): for age 29, there are 1000 simulated individuals. eIdxM(999, 29)
%                 shows the labor endomwnt index for the 999th individual
% eIdxM shows the labor endowment index, not the labor endowment itself!!!

% Example:
% My eIdxM is 3. It means that my labor endowment is the 3rd element in the
% labor endowment vector. It doesn't mean that my labor endowment is 3


%% Main

% Seed random number generator for repeatability
% It is important to use the same random numbers for every iteration over
% the guesses for calibration targeted parameter
% Otherwise simulated aggregates change a little bit every time, which will
% confuses Matlab equation solvers
rng(433);

% Endowment state by [ind, age]
eIdxM = MarkovChainSimulation(cS.nSim, cS.aD, paramS.leProb1V, ...
                              paramS.leTrProbM', rand([cS.nSim, cS.aD]));
                          
% Note:
% The transition matrix used in function 'MarkovChainSimulation' is 
% trProbM(s',s), which shows the prob of tomorrow's state being s' given
% today's state being s
% However, the transition matrix we computed using function 'tauchen' is
% trProbM(i,j) = Prob i --> j
% Hence in order to use the transition matrix computed by 'tauchen' as an
% input in 'MarkovChainSimulation', we have to transpose trProbM computed
% by 'tauchen'. 
% That is why the input in above code is [paramS.leTrProbM']


end