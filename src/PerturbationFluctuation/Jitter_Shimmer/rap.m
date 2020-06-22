function rP=rap( vPPS ) 

% Calculates the relative average disturbance of the pitch period in % (RAP).
% vPPS is the sequence of pitch periods.
% This RAP coincides with the disturbance ratio pq with a smoothing factor of
% Smoothing of 3 periods
%
% Input parameters:
%   vPPS is the sequence of peak amplitude values.
%
% Output parameters:
%   rP amplitude disturbance ratio in % (RAP).

rP=pq( vPPS, 3 );