function [vShimdB, vShimmer]=shimmer_env( vPAS, iNumPuntos )

% Calculates iNumPoints shimmer values (absolute value in dB, and relative value
% with respect to the average amplitude of the peaks in %).
% each value corresponding to absolute value of the difference between 
% adjacent pitch amplitudes (repeated values if necessary)
%
% Input parameters:
%   vPAS:       is the sequence of amplitude peaks of the the voice signal.
%   iNumPuntos: is the number of desired shimmer values
%
% Output parameters
%   vShimdB:    vector containing iNumPoints shimmer values in absolute
%               values in dB
%   vShimmer:   vector containing iNumPoints relative shimmer values (in %) 
%               with respect to the average amplitude of the peaks 

if nargin < 2, iNumPuntos = 100; end

iLongitud=length( vPAS );

if iLongitud < 2, error( 'The sequence of amplitude pitch periods is too short!' ); end
if iNumPuntos < 0, error( 'The number of points should be larger than 0' ); end

vShimdB=zeros(1, iNumPuntos);
vShimmer=zeros(1, iNumPuntos);

iDesplazamiento=iLongitud/iNumPuntos;

% Mean amplitude of the peaks
rMedia=mean( vPAS );

iIndiceIni=1;

for n=0:iNumPuntos-1
   indice=fix(iIndiceIni+n*iDesplazamiento);
   
   if ( n>0 ) && ( indice==fix( iIndiceIni+(n-1)*iDesplazamiento ) )
      vShimdB(n+1)=vShimdB(n);
      vShimmer(n+1)=vShimmer(n);
   else
      if indice+1<=iLongitud
         if vPAS( indice+1 )~=0 && vPAS( indice+1 )~=0
            vShimdB( n+1 )=abs( 20*log10( vPAS( indice+1 )/vPAS( indice ) ) );
         else
            vShimdB( n+1 )=0;
         end
         if rMedia~=0
            vShimmer( n+1 )=100*abs( vPAS( indice+1 )-vPAS( indice ) )/rMedia;
         else
            vShimmer( n+1 )=0;
         end
      else
         vShimdB( n+1 )=0;
         vShimmer( n+1 )=0;
      end
   end
end