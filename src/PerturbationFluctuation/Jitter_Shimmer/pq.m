function iPertQuot=pq( vSecuencia, iFactorProm )

% Calculates the disturbance ratio as a percentage of the signal vSequence,
% which can be a sequence of pitch values or peak amplitudes for example.
% Using the formula given by Kasuya.
%
% Input parameters:
%   vInput:       input sequence
%   iFactorProm:  is the averaging factor (must be odd).
%
% Output parameters:
%   iPertQuot:    disturbance ratio as a percentage of the signal

if nargin < 2
    iFactorProm = 5;
end

iLongitud=length( vSecuencia );

if iLongitud < 2
    error( 'The sequence of pitch periods is too short!' ); 
end
if iLongitud < iFactorProm 
    error( 'iLongitud < iFactorProm' ); 
end
if mod( iFactorProm, 2 ) == 0 
    error( 'iFactorProm must be odd' ); 
end

m=0.5*(iFactorProm-1);

rSuma=0;
for n=1:iLongitud-iFactorProm+1
    rDif=0;
    for r=1:iFactorProm
        rDif=rDif+vSecuencia( n+r-1 )-vSecuencia( n+m );
    end
    rDif=( 1/iFactorProm )*rDif;
    rDif=abs( rDif );
    rSuma=rSuma+rDif;
end
rNumerador=rSuma/( iLongitud-iFactorProm+1 );

suma2=0;
for n=1:iLongitud
    suma2=suma2+abs( vSecuencia(n) );
end
rDenominador=( 1/iLongitud )*suma2;

% Default value
iPertQuot = 0;
if rDenominador ~= 0 
    iPertQuot=100*( rNumerador/rDenominador );
end