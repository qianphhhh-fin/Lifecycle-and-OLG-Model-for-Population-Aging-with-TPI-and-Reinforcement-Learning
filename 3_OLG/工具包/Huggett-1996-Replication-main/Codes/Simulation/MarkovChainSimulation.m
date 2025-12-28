function eIdxM = MarkovChainSimulation(nSim, T, prob0V, trProbM, rvInM)
%% Documentation:
% This function simulates history of a Markov Chain 

% INPUTS:
% (1). nSim:          number of individuals to simulate
% (2). T:             length of histories
%                     In a finite-period model (e.g. OLG), T is the total 
%                     age that we set households to live
% (3). prob0V:        prob of each state at date 1
% (4). trProbM(s',s): transition matrix, showing the prob of state being s'
%                     tomorrow given the state today is s
% (5). rvInM:         uniform random variables, by [ind, t]

% OUTPUT:
% eIdxM: labor endowment index by [ind, age]


% ******************************** Notice ********************************* 
% We simulate nSim individuals for each age
% eIdxM(999, 29): for age 29, there are 1000 simulated individuals. eIdxM(999, 29)
%                 shows the labor endomwnt index for the 999th individual
% eIdxM shows the labor endowment index, not the labor endowment itself!!!
% Example:
% My eIdxM is 3. It means that my labor endowment is the 3rd element in the
% labor endowment vector. It doesn't mean that my labor endowment is 3
% *************************************************************************


% Codes are adapted from Prof. Hendricks, which are adapted from CompEcon
% Toolbox (Prof. Hendricks says it seems to be wrong...)


%% Input Validation
ns = length(prob0V);

% Check the number of input
if nargin ~= 5
   error('Invalid nargin');
end

validateattributes(trProbM, {'double'}, {'finite', 'nonnan', 'nonempty', ...
                   'real', '>=', 0, '<=', 1, 'size', [ns, ns]})

% Check if probabilities sum to one
prSumV = sum(trProbM);
if max(abs( prSumV - 1 )) > 1e-5
    error('Probabilities do not sum to one');
end

validateattributes(prob0V(:), {'double'}, {'finite', 'nonnan', 'nonempty', ...
                   'real', '>=', 0, '<=', 1, 'size', [ns,1]})

if abs(sum(prob0V) - 1) > 1e-5
    error('Initial probs do not sum to 1');
end
               
validateattributes(rvInM, {'double', 'single'}, {'finite', 'nonnan', ...
                   'nonempty', 'real', '>=', 0, '<=', 1, 'size', [nSim, T]})


%% Preliminaries
% For each state, find cumulative probability distribution for next period
cumTrProbM = cumsum(trProbM);
cumTrProbM(ns, :) = 1;

% Need to transpose this for the formula below now by [s, s']
cumTrProbM = cumTrProbM';


%%  Iterate over dates
eIdxM = zeros([nSim, T]);

% Draw t=1
eIdxM(:, 1) = 1 + sum((rvInM(:,1) * ones(1, ns)) > (ones(nSim,1) * cumsum(prob0V(:)')), 2);

% For t = 2, ..., T
for t = 1 : (T-1)
   eIdxM(:, t+1) = 1 + sum((rvInM(:,t+1) * ones(1, ns)) > cumTrProbM(eIdxM(:,t), :), 2);
end


end