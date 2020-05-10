function rCocPert=feijoo_pert( vSecuencia )

% Calculate the disturbance ratio proposed by Feijoo as a percentage, which is similar to Kasuya's PQ perturbation ratio.
%
% Input parameters:
%   vSecuencia:     sequence of amplitude values of the peaks of the cycles of the
%                   voice signal (the Feijoo coefficient A is obtained),
%                   or is the sequence of signal pitch periods
%                   (the coefficient P is obtained).
%
% Output parameters:
%   rCocPert:       disturbance ratio proposed by Feijoo in percent,

iLongitud=length( vSecuencia );
if iLongitud < 2, error( 'The signal is too short!' ); end

rMax=max( vSecuencia );

rSuma=0;
for i=1:iLongitud-1
   rSuma = rSuma + abs( vSecuencia(i+1)-vSecuencia(i) );
end

rCocPert = 100*(1/(iLongitud-1))*(1/rMax)*rSuma;