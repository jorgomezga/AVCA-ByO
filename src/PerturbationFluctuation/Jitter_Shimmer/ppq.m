function rP=ppq( vPPS ) 

% Calculates the pitch perturbation quotient in % (PPQ).
% Input parameters:
%   vPPS is the sequence of peak amplitude values.
%
% Output parameters:
%   rP pitch perturbation quotient in % (ppq).

rP = pq( vPPS, 11 );