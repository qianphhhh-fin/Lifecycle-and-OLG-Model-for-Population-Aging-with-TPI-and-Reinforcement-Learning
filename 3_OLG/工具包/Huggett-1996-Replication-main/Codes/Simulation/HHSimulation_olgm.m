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