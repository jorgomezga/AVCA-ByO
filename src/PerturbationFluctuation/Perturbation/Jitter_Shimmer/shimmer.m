function [rShdB, rShim]=shimmer( vPAS )

% Calculates the shimmer of the signal. Calculate its absolute value
% in dB (rShdB), and its relative value in % (rShim)
%
% Input parameters:
% vPAS:     is the sequence of amplitudes of the peaks of the cycles of the voice signal.
%
% Output parameters:
% rShdB:     absolute value in dB
% rShim:     relative value in %

iLongitud=length( vPAS );
if iLongitud < 2, error( 'The pitch period sequence is too short!' ); end

rShdBAcum=0; 
rShimAcum=0;

for i=1:iLongitud-1
   if vPAS( i )~=0 && vPAS( i+ 1)~=0
      ShdB_i= abs( 20*log10(vPAS(i+1)/vPAS(i)) );
   else
      ShdB_i=0;
   end
   rShdBAcum = rShdBAcum + ShdB_i;

   shima_i= abs( vPAS(i+1)-vPAS(i) );
   rShimAcum = rShimAcum + shima_i;
end

% The shimmer in dB will be the average of the accumulated shimmer in dB
rShdB=rShdBAcum/(iLongitud-1);

% The average value of the amplitudes of the peaks is
rMed = mean( vPAS );

% Relative shimmer
if rMed~=0 
    rShimAcum=100*rShimAcum/rMed; 
else
    rShimAcum=0;
end

rShim=rShimAcum/(iLongitud-1);