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