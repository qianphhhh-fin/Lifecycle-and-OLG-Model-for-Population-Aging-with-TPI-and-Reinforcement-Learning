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