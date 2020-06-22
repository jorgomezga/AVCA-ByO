function rA=sapq(vPAS, iFactor ) 

% Calculates the smoothed amplitude perturbation quotient in % (sAPQ)
% from the sequence of pitch amplitudes (PAS)
% This sAPQ coincides with the perturbation quotient pq using the
% sequence of pitch amplitudes PAS
%
% Input parameters:
%   vPAS:      is the sequence of peak amplitude values.
%   iFactor:   is the smoothing factor
%
% Output parameters:
%   rA:        amplitude perturbation quotient in % (sAPQ).

if nargin < 2
    iFactor = 5; 
end

try 
    rA=pq( vPAS, iFactor );
catch 
    warning('sapq not calculated, returning NaN')
    rA = NaN;
end