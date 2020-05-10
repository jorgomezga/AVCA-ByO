function iF0=PitchCeps( vFrame, iFs, iT0 )

% Pitch estimation using the LPC method of Rabiner
%
% Inputs:
%   vFrame               = Input vFrame
%   iFs                  = Sampling frequency
%   iT0                  = Previous calculation of the pitch period
% Outputs:
%   vClippedFrame        = Clipped frame


iVentana=length( vFrame );

% If the previous segment was voiced (pitch period To > 0),
% the current segment pitch period is looked for, around the 30% of T0
% (as long as the range does not exceed 2 ms - 15 ms). If not, the
% pitch period is looked in the range between the values of 2 ms and 15 ms.
if iT0~=0
   iTmin=max( round(0.7*iT0), round(0.002*iFs) );
   iTmax=min( round(1.3*iT0), round(0.015*iFs) );
else
   iTmin=round(0.002*iFs);
   iTmax=round(0.015*iFs);
end

ENERGY_THRESH=20;

% Real cesptrum calculation
iN_FFT=2^( ceil( log2( iVentana )));
vFrame=vFrame.*hamming( iVentana );
vCepstrum = real( ifft( log( abs( fft( vFrame, iN_FFT ))), iN_FFT) );
vCepstrum = vCepstrum( 1:iN_FFT/2 );
VCepsAlta = vCepstrum( iTmin:iN_FFT/2 );

% Cepstral energy
vCepsEnergy = VCepsAlta.^2;

[iMax,ind_max]=max( vCepsEnergy(1:iTmax-iTmin ));
      
% Inicializamos el pitch a 0. Como si fuera no sonora 
iF0 = 0;
if ( iMax > ENERGY_THRESH*mean( vCepsEnergy ) )
      
   % Previous maximum values (of lower order and sufficient amplitude) are
   % looked  for. Only the last value (the first above the threshold) is
   % kept
   [max2,i_max2]=max( vCepsEnergy( 1:ind_max-iTmin ) );
   while max2>ENERGY_THRESH*mean( vCepsEnergy )
      ind_max=i_max2;
      [max2,i_max2]=max( vCepsEnergy(1:i_max2-iTmin) );
   end
		
   % i_max+Tmin-1 es el index in the vector of cepstrums
   iT=ind_max+iTmin-1-1;
   if iT~=0
      iF0=iFs/iT;
   end
   
end