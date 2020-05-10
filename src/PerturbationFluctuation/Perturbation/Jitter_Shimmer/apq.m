function rA=apq(vPAS) 

% Calculates the amplitude disturbance ratio in % (APQ).
% This APQ coincides with the perturbation quotient pq with a Smoothing factor of
% of 11 periods
%
% Input parameters:
%       vPAS is the sequence of peak amplitude values.
%
% Output parameters:
%       rA amplitude disturbance ratio in% (APQ).

rA=pq( vPAS, 11 );