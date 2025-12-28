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