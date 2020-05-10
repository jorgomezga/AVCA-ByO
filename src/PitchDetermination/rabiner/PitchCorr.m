function iF0 = PitchCorr( vFrame, iFs, iT0 )

% Pitch estimation using the LPC method of Rabiner
%
% Inputs:
%   vFrame               = Input vFrame
%   iFs                  = Sampling frequency
%   iT0                  = Previous calculation of the pitch period
% Outputs:
%   vClippedFrame        = Clipped frame

iLengthFrame = length( vFrame );

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

% vA are LPC coefficients and rG is the gain of the filter
% Autocorrelation is calcuated using the prediction error
[vA,rG]   = LPC2( vFrame,4, 0);
vErrorPred= filter(vA, rG, vFrame);

vAutocorr = xcorr( vErrorPred );
vAutocorr = vAutocorr(iLengthFrame+iTmin:2*iLengthFrame-1).^2;

[rMaximo1, iMaximo1]=max( vAutocorr( 1:iTmax-iTmin ) );

% Pitch initialized to 0, as it was unvoiced 
iF0=0;
if (rMaximo1 > ENERGY_THRESH*mean( vAutocorr )) 
      
   % Previous maximum values (of lower order and sufficient amplitude) are
   % looked  for. Only the last value (the first above the threshold) is
   % kept
   [~,iMaximo2]=max( vAutocorr(1:iMaximo1-iTmin) );
   while iMaximo2 > ENERGY_THRESH*mean( vAutocorr )
      iMaximo1=iMaximo2;
      [~, iMaximo2]=max( vAutocorr(1:iMaximo2-iTmin) );
   end
		
   iT=iMaximo1+iTmin-1;
   if iT~=0
      iF0=iFs/iT;
   end
   
end