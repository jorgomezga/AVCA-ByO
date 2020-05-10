function rP=sppq( vPPS, iFactor ) 

% Calculates the disturbance ratio of the smoothed amplitude in % (sPPQ) 
% from the sequence of pitch periods (PPS)
% This sPPQ coincides with the disturbance quotient pq when
% the sequence of PPS pitch periods is used
%
% Input parameters:
%   vPPS:       is the sequence of peak amplitude values.
%   iFactor:   is the smoothing factor
%
% Output parameters:
%   rP:        amplitude disturbance ratio in % (sPPQ).

if nargin < 2
    iFactor = 5; 
end

rP=pq( vPPS, iFactor );