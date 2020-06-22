function [vJitta, vJitter]=jitter_env(vPPS, iNumPuntos)

% Computes iNumPoints jitter values (absolute jitter in microseconds, and relative
% jitter compared to the average pitch period).
% Each value corresponding to the absolute value of the difference between 
% adjacent pitch periods (values are repeated if necessary)
%
% Input parameters:
%   vPPS:       is the sequence of pitch periods in seconds of the signal.
%   iNumPuntos: is the number of desired jitter values
%
% Output parameters:
%   vJitta:    vector containing iNumPoints absolute jitter values in
%              microseconds
%   vJitter:   vector containing iNumPoints values of jitter relative to
%              average pitch period

if nargin < 2, iNumPuntos = 100; end 

iLongSec=length(vPPS);

if iLongSec < 2, error( 'The sequence of pitch periods is too short!' ); end
if iNumPuntos < 0, error( 'The number of points should be larger than 0' ); end

vJitta=zeros(1,iNumPuntos);
vJitter=zeros(1,iNumPuntos);

iIndiceIni=1;

iDesplazamiento=iLongSec/iNumPuntos;

% Mean period of the signal in microseconds
rFoMed=1e6*mean( vPPS );

for n=0:iNumPuntos-1
   indice=fix( iIndiceIni+n*iDesplazamiento );
   
   if ( n>0 ) && (indice==fix( iIndiceIni+(n-1)*iDesplazamiento ))
      vJitta( n+1 )=vJitta( n );
      vJitter( n+1 )=vJitter( n );
   else
      if indice+1 <= iLongSec
         vJitta( n+1 )=1e6*abs( vPPS(indice+1)-vPPS(indice) );

         if rFoMed ~= 0
            vJitter( n+1 )=100*vJitta(n+1)/rFoMed;
         else
            vJitter( n+1 )=0;
         end
      else
         vJitta( n+1 )=0;
         vJitter( n+1 )=0;
      end
   end
end