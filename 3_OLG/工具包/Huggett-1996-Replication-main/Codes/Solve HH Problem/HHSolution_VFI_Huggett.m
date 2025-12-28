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