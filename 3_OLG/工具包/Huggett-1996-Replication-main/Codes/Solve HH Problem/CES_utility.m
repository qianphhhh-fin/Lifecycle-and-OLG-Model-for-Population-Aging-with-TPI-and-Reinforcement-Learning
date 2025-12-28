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