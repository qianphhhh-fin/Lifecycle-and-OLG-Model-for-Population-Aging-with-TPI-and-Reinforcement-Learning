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