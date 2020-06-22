function rP=sppq( vPPS, iFactor ) 

% Calculates the smoothed pitch perturbation quotient in % (sPPQ) 
% from the sequence of pitch periods (vPPS)
%
% Input parameters:
%   vPPS:      is the sequence of peak amplitude values.
%   iFactor:   is the smoothing factor
%
% Output parameters:
%   rP:        smoothed pitch perturbation quotient in % (sPPQ).

if nargin < 2
    iFactor = 5; 
end

try 
    rP = pq( vPPS, iFactor );
catch
    warning('rSPPQ not calculated, returning NaN');
    rP = NaN;
end