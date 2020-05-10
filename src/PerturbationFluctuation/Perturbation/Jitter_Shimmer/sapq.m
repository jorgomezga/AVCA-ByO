function rA=sapq(vPAS, iFactor ) 

% Calculates the smoothed amplitude disturbance ratio in % (sAPQ)
% from the sequence of pitch amplitudes (PAS)
% This sAPQ coincides with the perturbation quotient pq using the
% sequence of pitch amplitudes PAS
%
% Input parameters:
%   vPAS:      is the sequence of peak amplitude values.
%   iFactor:   is the smoothing factor
%
% Output parameters:
%   rA:        amplitude disturbance ratio in% (sAPQ).

if nargin < 2
    iFactor = 5; 
end

rA=pq( vPAS, iFactor );