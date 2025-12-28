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