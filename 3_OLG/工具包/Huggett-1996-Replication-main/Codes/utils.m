function [Y, R, w, b] = HHPrices_Huggett(K, L, cS)
%% Documentation:
% This function is defined based on production function

% INPUTS
% (1). K:        aggregate capital in the economy
% (2). L:        aggregate labor supply
%                This term can be precalibrated without solving HH problem
% (3). cS.A:     productivity
% (4). cS.alpha: capital share
% (5). cS.ddk:   depreciation rate
% (6). cS.theta: social security tax
% (7). cS.retireMass: mass of retired households

% OUTPUTS
% (1). Y: production output
% (2). R: capital rental price faced by households
%         R = 1+ (MPK - depreciation)(1 - tax)
% (3). w: wage after tax
%         w = MPL(1 - tax - social_security)
% (4). b: social security benefit
%         (social_security x w x L) = b x all_retired_households

% ******************* Note ******************* 
% tax is calibrated in this function as 
% 0.195/(1-depreciation x K/Y)


%% Main
% Production
Y   = cS.A * (K^cS.alpha) * (L^(1-cS.alpha));
MPK = cS.alpha * cS.A * (K^(cS.alpha-1)) * (L^(1-cS.alpha));
MPL = (1-cS.alpha) * cS.A * (K^cS.alpha) * (L^(-cS.alpha));

% Tax
tau = 0.195/(1-cS.ddk * K/Y);

% Prices faced by households
R   = 1 + (MPK - cS.ddk)*(1 - tau);
w   = MPL*(1 - tau - cS.theta);

% Social security benefits
b   = cS.theta * w * L/cS.retireMass;


end
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
function [y, trProbM] = tauchen(N, pRho, pSigma, pMu, n_std)
%% Documentation
% Use Tauchen's (1986) method to produce finite state Markov approximation 
% of the AR(1) processes:
%           y_t = pMu + pRho y_{t-1} + epsion_t,
%            where epsilon_t ~ N (0, pSigma^2)

% INPUTS:
% (1). N:       Number of points in markov process
% (2). pRho :   Persistence parameter in AR(1) process
% (3). pSigma : Standard deviation of random component of AR(1) process
% (4). pMu :    optional(default=0.0)
%               Mean of AR(1) process
% (5). n_std :  optional(default=3)
%               # of std deviations to each side the processes should span

% OUTPUTS:
% (1). y :       array(dtype=float, ndim=1)
%                1d-Array of nodes in the state space
% (2). trProbM : array(dtype=float, ndim=2)
%                Matrix transition probabilities for Markov Process
%                trProbM(i,j) = Prob i -> j

% Adapted from Prof. Lutz Hendricks' code which is adapted from QuantEcon Julia code


%% Get discretized space
% Width of grid
a_bar = n_std * sqrt(pSigma^2 / (1 - pRho^2));

% Grid
y     = linspace(-a_bar, a_bar, N);

% Distance between points
d     = y(2) - y(1);


%% Get transition probabilities
trProbM  = zeros(N, N);

for iRow = 1 : N

   % Do end points first
   trProbM(iRow, 1) = normcdf((y(1) - pRho*y(iRow) + d/2) / pSigma);
   trProbM(iRow, N) = 1 - normcdf((y(N) - pRho*y(iRow) - d/2) / pSigma);

   % fill in the middle columns
   for iCol = 2 : N-1
      trProbM(iRow, iCol) = (normcdf((y(iCol) - pRho*y(iRow) + d/2) / pSigma) - ...
                             normcdf((y(iCol) - pRho*y(iRow) - d/2) / pSigma));
   end

end


% NOTE: 

% We need to shift this vector after finding probabilities because when
% finding the probabilites we use the function 'normcdf', which assumes its
% input argument is distributed N(0,1). After adding, the mean E(y) is no 
% longer 0. Hence we would be passing elements with the wrong distribution.

% It is ok to do after the fact because adding this constant to each term 
% effectively shifts the entire distribution. Because the normal distribution
% is symmetric and we just care about relative distances between points, 
% the probabilities will be the same.

% We could have shifted it before, but then we would need to evaluate the 
% cdf with a function that allows the distribution of input arguments to be 
% (pMu/(1 - pRho), 1) instead of (0, 1)


%% Center process around its mean (ybar / (1 - rho))
y = y + pMu / (1 - pRho); 


end
function [massV, lbV, ubV] = norm_grid(xV, xMin, xMax, mu, sig)
%% Documentation:

% Approximate a Normal distribution on a given grid
% Given a Normal distribution with parameters mu, sig and a grid of points
% (xV), this function returns the mass in the interval around each xV 
% implied by N(mu, sig^2) where the Normal is truncated at [xMin, xMax]
% This only works well if the grid is sufficiently tight!

% INPUTS:
% (1). xV:   a grid of points
% (2). xMin: lower bound of smallest interval
% (3). xMax: upper bound of largest  interval
% (4). mu:   mean of the normal distribution
% (5). sig:  standard deviation of the normal distribution
      
% OUTPUTS:
% (1). massV: mass on each grid point
%             Not a density!
%             DensityV(i) = massV(i) / (ubV(i) - lbV(i))
% (2). lbV:   interval lower bound
% (3). ubV:   interval upper bound

% Adapted from Prof. Lutz Hendricks' code


%% Input Validation

% Check-1 the number of inputs
if nargin ~= 5
   error([ mfilename, ': Invalid nargin' ]);
end

n = length(xV);

% Check-2 if xV is increasing
if any( xV(2:n) < xV(1:n-1) )
    warnmsg([ mfilename, ':  xV not increasing' ]);
    keyboard;
end

% Check-3 if xMin is the smallest in xV, xMax is the largest in xV
if xMin > xV(1)  ||  xMax < xV(n)
    warnmsg([ mfilename, ':  Invalid xMin or xMax' ]);
    keyboard;
end

% Check-4 if mu lies between xMin and xMax
if mu < xMin  ||  mu > xMax
    warnmsg([ mfilename, ':  Invalid mu' ]);
    keyboard;
end


%% Main
% Construct interval boundaries
% Symmetric around xV
xMidV = 0.5 .* (xV(1:n-1) + xV(2:n));
lbV   = [xMin; xMidV(:)];
ubV   = [xMidV(:); xMax];

% Find mass in each interval
cdfV  = normcdf( [xMin; ubV], mu, sig );
massV = cdfV(2:(n+1)) - cdfV(1:n);
massV = massV(:) ./ sum(massV);


%% Output check
validateattributes(massV, {'double'}, {'finite', 'nonnan', 'nonempty', 'real', 'size', [n,1], ...
                   '>=', 0, '<=', 1});

validateattributes(lbV, {'double'}, {'finite', 'nonnan', 'nonempty', 'real', 'size', [n,1], ...
                   '>=', xMin, '<=', xMax});

validateattributes(ubV, {'double'}, {'finite', 'nonnan', 'nonempty', 'real', 'size', [n,1], ...
                   '>=', xMin, '<=', xMax});

if any(ubV < lbV)
    error('ubV < lbV');
end
 

end
function [HHlaborM, L] = LaborSupply_Huggett(eIdxM, cS, paramS)
%% Documentation:
% This function computes labor supply by [ind, age] and aggregate labor supply
% in Huggett (1996) model
% Since there is random mortality, mass of households by age varies across age

% INPUT:
% eIdxM: labor endowment we have simulated for all ages and for each
%        age there are nSim simulated individuals

% OUTPUTS:
% (1). LSHistM: individual labor supply by [ind, age]
% (2). L:       aggregate labor supply 


%% Main:
% Individual labor supplies: efficiency * labor endowment shock
HHlaborM = zeros(cS.nSim, cS.aD);
for a = 1 : cS.aD
   HHlaborM(:, a) = paramS.ageEffV(a) .* paramS.leGridV(eIdxM(:,a));
end


% Aggregate labor supply
% Let the mass of households with age a be mu(a), then the aggregate labor 
% supply is
% L = sum_(a=1 to aD) mu(a) x mean( HHlaborM(:,a) )
L = mean(HHlaborM,1)* cS.ageMassV';


end
function [logGridV, trProbM, prob1V] = EarningProcess_olgm(cS)
%% Documentation:
% Calibrate labor endowment process
% This is net of the age efficiency profile!

% INPUTS: cS.
% (1). nw:            # of labor endowment states
% (2). lePersistence: persistence parameter in the AR(1) process
% (3). leShockStd:    std dev of random component in the AR(1) process
% (4). leWidth:       # of standard deviations the grid is wide
% (5). leSigma1:      standard deviation of normal distribution that labor
%                     endowments in the initial period follow

% OUTPUTS:
% (1). logGridV: log grid of endowment states
% (2). trProbM:  trProbM(i,j) = Prob i -> j
% (3). prob1V:   stationary and age 1 distribution


%% Main
[logGridV, trProbM] = tauchen(cS.nw, cS.lePersistence, cS.leShockStd, 0, cS.leWidth);

% New agents draw from an approximate log normal distribution
% On the grid defined by the AR(1)
prob1V              = norm_grid(logGridV, logGridV(1)-2, logGridV(end)+2, 0, cS.leSigma1);
prob1V              = prob1V(:);

% Improve scaling
logGridV            = logGridV(:) - logGridV(1) - 1;


%% Output Validation
validateattributes(trProbM, {'double'}, {'finite', 'nonnan', 'nonempty', 'real', '>=', 0, ...
                   '<', 1, 'size', [cS.nw, cS.nw]});

pSumV = sum(trProbM, 2);
if any(abs(pSumV - 1) > 1e-6)
    error('probs do not sum to 1');
end

validateattributes(prob1V, {'double'}, {'finite', 'nonnan', 'nonempty', 'real', '>=', 0, ...
                   'size', [cS.nw, 1]});
             
if abs(sum(prob1V) - 1) > 1e-6
    error('prob1V does not sum to 1');
end


end
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
function [kHistM, cHistM] = HHSimulation_olgm(kPolM, cPolM, eIdxM, cS)
%% Documentation:
% This function simulates a population of households
% The basic idea is
% (1). Populate a set of households
% (2). Households go through sequence of labor endowments given in eIdxM
% (3). Compute capital holdings of these households based on policy function kPolM

% INPUTS:
% (1). kPolM: k' policy function, by [ik, ie, a]
% (2). cPolM: policy function for consumption, by [ik, ie, a]
% (3). eIdxM: labor endowment index for each simulated individal

% OUTPUTS:
% (1). kHistM: capital stock histories for households by [ind, age]
% (2). cHistM: consumption histories for households by [ind, age]


%% Input Validation
validateattributes(kPolM, {'double'}, {'finite', 'nonnan', 'nonempty', 'real', ...
                   'size', [cS.nk,cS.nw,cS.aD]})

               
%% Simulate capital and consumption histories, age by age
nSim   = size(eIdxM, 1);
kHistM = zeros(nSim, cS.aD);
cHistM = zeros(nSim, cS.aD);

for a = 1 : cS.aD
    
   for ie = 1 : cS.nw
      % Find households with labor endowment ie at this age
      idxV = find(eIdxM(:,a) == ie);

      if ~isempty(idxV)
         if a < cS.aD
            % Find next period capital for each individual by interpolation
            kHistM(idxV, a+1) = interp1(cS.kGridV(:), kPolM(:,ie,a), ...
                                        kHistM(idxV, a), 'linear');
         end
         
         cHistM(idxV, a) = interp1(cS.kGridV(:), cPolM(:,ie,a), ...
                                   kHistM(idxV, a), 'linear');
      end % if idexV is not empty

   end % for ie

end % for a



%% Output Validation
validateattributes(kHistM, {'double'}, {'finite', 'nonnan', 'nonempty', 'real', ...
                   '>', cS.kGridV(1) - 1e-6})

               
end

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

function incomeM = HHIncome_Huggett(kV, R, w, T, b, a, paramS, cS)
%% Documentation:
% This function computes income of a household GIVEN HIS MODEL AGE
% Household Income = Non-capital income + Capital Income

% INPUTS
% (1). kV:      capital grids, it is a (nk x 1) vector
% (2). R and w: prices faced by households
%               R = 1+ (MPK - depreciation)(1 - tax)
%               w = MPL(1 - tax - social_security)
% (3). T:       transfer by accidental bequest
% (4). b:       social security benefit received when retired
% (4). a:       age

% OUTPUT
% incomeM: income of household at a given age
%          it is a (nk x nw) matrix


%% Non-capital income (by shock)
% nonCapIncomeV: (nw x 1) vector

% Non-capital income for each given age only depends on exogenous state
% variable e, the realization of labor endowment.

% Hence for each working HH, nonCapIncomeV is a (nw x 1) vector, each
% element corresponding to each labor endowment realization

% For retired HH, a>aR, they have the same non-capital income -- transfer
% In order to write in the same way as when they are working, we still
% write the non-capital income as a (nw x 1) vector. But in this case, each
% element in nonCapIncomeV takes the same value

if a <= cS.aR
   nonCapIncomeV = paramS.ageEffV(a) .* w .* paramS.leGridV + T;
else
   nonCapIncomeV = b .* ones([cS.nw, 1]) + T;
end


%% Total Income at a given age
% incomeM is a (nk x nw) matrix
%     k_1 + nw kinds of labor income
%     k_2 + nw kinds of labor income
%         ...      ...
%     k_nk + nw kinds of labor income

incomeM = R * kV(:) * ones([1, cS.nw]) + ones([length(kV),1]) * nonCapIncomeV(:)';


end

function [cPolM, kPolM, valueM] = HHSolution_VFI_Huggett(R, w, T, bV, paramS, cS)
%% Documentation:
% This function solves household problem in the multi-period OLG model with
% earning uncertainty using VALUE FUNCTION ITERATION

% INPUTS:
% (1). R, w:       prices faced by household, R=1+r
% (2). T:          lump sum transfers of accidental bequests
% (3). bV:         social security benefits
%                  It is a (aD x 1) vector: 
%                  for HH with age 1~aR, social security benefit is 0
%                  for HH with age aR+1~aD, social security benefit is b
% (4). paramS, cS: other parameters

% OUTPUTS:
% (1). cPolM:  policy function for consumption by each state [ik, ie, age]
% (2). kPolM:  policy function for saving by each state [ik, ie, age]
% (3). valueM: value function at each state [ik, ie, age]


%*********** Note ************
% cPolM, kPolM and valueM are 3-dimensional matrix
% 1-D: nk, each capital grid
% 2-D: nw, each realization of labor endowment
% 3-D: aD, each age


%% Main
% Initialize policy function by [ik, ie, age]
cPolM  = zeros(cS.nk, cS.nw, cS.aD);
kPolM  = zeros(cS.nk, cS.nw, cS.aD);
valueM = zeros(cS.nk, cS.nw, cS.aD);

% Backward induction

for a = cS.aD : -1 : 1

   % Next period value function
   if a < cS.aD
      vPrimeM = valueM(:,:,a+1);
   else
   % There is no next period, since a=aD is the last period
      vPrimeM = [];
   end

   [cPolM(:,:,a), kPolM(:,:,a), valueM(:,:,a)] = ...
      HHSolutionByAge_VFI_Huggett(a, vPrimeM, R, w, T, bV(a), paramS, cS);
  
end


%% Output Validation
validateattributes(cPolM, {'double'}, {'finite', 'nonnan', 'nonempty', ...
                   'real', 'positive', 'size', [cS.nk, cS.nw, cS.aD]})

validateattributes(valueM, {'double'}, {'finite', 'nonnan', 'nonempty', ...
                   'real', 'size', [cS.nk, cS.nw, cS.aD]})

validateattributes(kPolM, {'double'}, {'finite', 'nonnan', 'nonempty', ...
                   'real', 'size', [cS.nk, cS.nw, cS.aD], ...
                   '>=', cS.kMin - 1e-6, '<=', cS.kGridV(end) + 1e-6})


end

function [cPolM, kPolM, valueM] = HHSolutionByAge_VFI_Huggett(a, vPrime_keM, R, w, T, b, paramS, cS)
%% Documentation:
% This function solves the household problem for a GIVEN AGE using value 
% funtion iteration. 
% ***************************************************************
% We take the value function for age a+1 as given, unless a=aD *
% ***************************************************************

% INPUTS:
% (1). a:               current age
% (2). vPrimeM(ik, ie): value function for age a+1
%                       Ignored if a = aD
% (3). R, w:            prices for k and labor faced by households
% (4). T:               lump sum transfers of accidental bequests
% (5). b:               social security benefits at age a

% OUTPUTS:
% cPolM, kPolM: Policy functions, indexed by [ik, ie]
% valueM:       Value function, indexed by [ik, ie]

%*********** Note ************
% Since cPolM, kPolM and valueM here are policy functions and value function 
% for a given age, hence they are only 2-dimensional
% 1-D: nk, each capital grid
% 2-D: nw, each realization of labor endowment


%% Input check
if a < cS.aD
    if ~isequal(size(vPrime_keM), [cS.nk, cS.nw])
        error('Invalid size of cPrimeM');
    end
end


%% Main
% Income y at age a
% yM is a (nk x nw) matrix: row each grid for k; column each labor realization
yM               = HHIncome_Huggett(cS.kGridV, R, w, T, b, a, paramS, cS);

% Options for optimization
fminbndOptS      = optimset('fminbnd');
fminbndOptS.TolX = 1e-5;

if a == cS.aD
   % Eat all income and save nothing, since this is the last period
   cPolM       = yM;
   kPolM       = zeros(cS.nk, cS.nw);
   [~, valueM] = CES_utility(cPolM, cS.sigma);
   
else
   % Allocate space for policy functions and value function
   cPolM       = zeros(cS.nk, cS.nw);
   kPolM       = zeros(cS.nk, cS.nw);
   valueM      = zeros(cS.nk, cS.nw);

   % Loop over states [ik, ie]
   for ie = 1 : cS.nw
       
      % Expected value function, by kPrime grid point -- EV(k')
      % ExValuePrimeV is a (nk x 1) vector, for each capital, the expected 
      % value over all possible labor endowment realizations in next period
      ExValuePrimeV = zeros(cS.nk, 1);
      for ik = 1 : cS.nk
         ExValuePrimeV(ik) = paramS.leTrProbM(ie,:) * vPrime_keM(ik,:)';
      end

      % Continuous approximation of tomorrows EV(k')
      % vPrimeOfK is a function which approximate the discrete relationship
      % between parmaS.kGridV and ExValuePrimeV by using linear interpolation. 
      vPrimeOfK = griddedInterpolant(cS.kGridV, ExValuePrimeV, 'linear');

      % Loop over capital states
      for ik = 1 : cS.nk
            [cPolM(ik,ie), kPolM(ik,ie), valueM(ik,ie)] = ...
                  HHSolutionByOneState_VFI_Huggett(a, yM(ik,ie), R, vPrimeOfK, fminbndOptS, cS);
      end % for ik

   end % for ie

end % for the for-loop


%% Output Validation
validateattributes(cPolM, {'double'}, {'finite', 'nonnan', 'nonempty', ...
                   'real', 'positive', 'size', [cS.nk, cS.nw]})

validateattributes(kPolM, {'double'}, {'finite', 'nonnan', 'nonempty', 'real', ...
                   '>=', cS.kMin - 1e-6,  '<=', cS.kGridV(cS.nk) + 1e-6, ...
                   'size', [cS.nk, cS.nw]})

validateattributes(valueM, {'double'}, {'finite', 'nonnan', 'nonempty', ...
                   'real', 'size', [cS.nk, cS.nw]})


end

function [c, kPrime, ValueFunc] = HHSolutionByOneState_VFI_Huggett(a, y, R, vPrimeOfK, fminbndOptS, cS)
%% Documentation:
% This function ruturns policy functions (for consumption and saving) and
% value function when my current state is (ik, ie, a), by using
% VALUE FUNCTION ITERATION
%--------------------------------------------------------------------------

% Budget constraint is: k' = y - c

% INPUTS:
% (1). a:           today's age
% (2). y:           today's income, given my state is (ik, ie)
% (3). R:           price of renting capital faced by household
% (4). vPrimeOfK:   expected value next period
%                   It is a function of kPrime, as a griddedInterpolant
% (5). fminbndOptS: options for fminbnd (searching for the opt c)


%% Main
% Range of feasible kPrime
kPrimeMax = min(cS.kGridV(cS.nk), y - cS.cFloor);

if kPrimeMax <= cS.kGridV(1)
   % No feasible choice. Household gets c floor and saves nothing
   kPrime = cS.kGridV(1);
   
else
   % Find optimal kPrime
   [kPrime, ~, ~] = fminbnd(@Bellman, cS.kGridV(1), kPrimeMax, fminbndOptS);
   
end

[ValueFunc, c] = Bellman(kPrime);
ValueFunc      = -ValueFunc;


%% Output Validation
validateattributes(kPrime, {'double'}, {'finite', 'nonnan', 'nonempty', 'real', 'scalar', ...
                   '>=', cS.kGridV(1) - 1e-6, '<=', kPrimeMax + 1e-6})

validateattributes(c, {'double'}, {'finite', 'nonnan', 'nonempty', 'real', 'scalar', ...
                   '>=', cS.cFloor - 1e-6})

validateattributes(ValueFunc, {'double'}, {'finite', 'nonnan', 'nonempty', 'real', 'scalar'})


%% Nested: objective function

% RHS of Bellman x (-1)
% We should make it negative, because fminbnd finds the solution which gives
% the smallest value of objective function
% Hence fminbnd finds 
% c = argmin -(RHS of Bellman), which is equivalent to
% c = argmax RHS of Bellman

    function [Valfunc, c] = Bellman(kPrime)
        c       = max(cS.cFloor, y - kPrime);
        [~, u]  = CES_utility(c, cS.sigma);
        Valfunc = -(u + cS.beta *cS.s(a)* R * vPrimeOfK(kPrime));
    end


end