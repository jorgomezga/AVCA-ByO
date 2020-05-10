function [rJitta, rJitt]=jitter( vPPS )

% Calculates the jitter of the voice signal. Calculate its absolute value 
% in microseconds, and its relative value in % with respect
% to the average value of the period
%
% Input parameters:
% vPPS:     is the sequence of pitch periods in seconds of the signal.
%
% Output parameters:
% rJitta:   absolute value in microseconds
% rJitt:   relative value in % with respect to the average value of the period

iLongSec = length( vPPS );
if iLongSec < 2
    error( 'VPPS is too short!' ); 
end

rJittaAcum = 0;
for i=1:iLongSec-1
   rJittaAcum = rJittaAcum+abs(vPPS(i+1)-vPPS(i));
end

% The absolute jitter will be the average of the accumulated absolute jitter
% in microseconds
rJitta = 1e6*rJittaAcum/(iLongSec-1); 

% The average period of the signal in microseconds
rMedio=1e6*mean( vPPS );

% So the relative jitter will be
if rMedio~=0 
    rJitt=100*rJitta/rMedio;
else
    rJitt=0;
end 