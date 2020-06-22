function rA=apq(vPAS) 

% Calculates the amplitude perturbation quotient % (APQ).
% This APQ coincides with the perturbation quotient pq with a Smoothing factor of
% of 11 periods
%
% Input parameters:
%       vPAS is the sequence of peak amplitude values.
%
% Output parameters:
%       rA amplitude perturbation quotient % (APQ).

try
    rA = pq( vPAS, 11 );
catch
    warning('apq not calculated, returning NaN')
    rA = NaN;
end