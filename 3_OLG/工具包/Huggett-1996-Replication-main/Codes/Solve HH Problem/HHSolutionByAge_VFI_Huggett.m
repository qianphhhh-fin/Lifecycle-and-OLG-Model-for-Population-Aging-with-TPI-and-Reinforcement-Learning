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