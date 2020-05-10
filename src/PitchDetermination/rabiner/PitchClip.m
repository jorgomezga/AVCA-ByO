function iF0 = PitchClip( vFrame, iFs, iT0 )

% Pitch estimation using static clipping 
%
% Inputs:
%   vFrame            = Input vFrame
%   iFs               = Sampling frequency
%   iT0               = Previously calculated pitch period
% Outputs:
%   vClippedFrame     = Clipped frame

if nargin < 2, error( 'Not enough input parameters!' ); end
if nargin < 3, iT0 = 0; end

% If the previous segment was voiced (pitch period To > 0),
% the current segment pitch period is looked for, around the 30% of T0
% (as long as the range does not exceed 2 ms - 15 ms). If not, the
% pitch period is looked in the range between the values of 2 ms and 15 ms.
if iT0~=0
   iTmin=max( round(0.7*iT0), round(0.002*iFs) );
   iTmax=min( round(1.3*iT0), round(0.015*iFs) );
else
   iTmin=round( 0.002*iFs );
   iTmax=round( 0.015*iFs );
end

ENERGY_THRESH=20; 

iLengthFrame = length( vFrame );
iMeanFrame   = mean( vFrame );

% Whitening to cleanse the input signal and obtain a more clear peak of the
% maximum of the autocorrelation function
% The maximum and minimum signal amplitudes in the first and last thirds of the signal 
% are computed, using a 65% positive clipping level of the minimum value of
% the obtained maxima and a negative clipping level of the maximum value of
% the obtained minima
iMax1=max( vFrame(1:fix(iLengthFrame/3)) );
iMax2=max( vFrame(fix(2*iLengthFrame/3):iLengthFrame) );
rAmax=min( iMax1, iMax2 ); 
if rAmax <= iMeanFrame
   rAmax=max( vFrame );
end
rCLpos=0.65*( rAmax-iMeanFrame );
rTRmax=rCLpos+iMeanFrame;

iMin1=min( vFrame( 1:fix(iLengthFrame/3)) );
iMin2=min( vFrame( fix(2*iLengthFrame/3):iLengthFrame) );
rAmin=max( iMin1, iMin2 );
if rAmin >= iMeanFrame
   rAmin=min( vFrame );
end
rCLneg=0.65*( iMeanFrame-rAmin );
rTRmin=-rCLneg+iMeanFrame;

vSigClip=zeros( 1, iLengthFrame );
for i=1:iLengthFrame
   vSigClip(i)=iMeanFrame; 
end

for i=1:iLengthFrame   
   if vFrame(i) > rTRmax
      vSigClip(i)=vFrame(i);
   elseif vFrame(i) < rTRmin
      vSigClip(i)=vFrame(i);      
   end    
end
   
% Autocorrelation and elimination of repeated information
vAutocorr=xcorr( vSigClip );
vAutocorr=vAutocorr( iLengthFrame+iTmin:2*iLengthFrame-1 ).^2;
% Max of ACF
[rMax, iMax]=max( vAutocorr( 1:iTmax-iTmin ));
   
% iF0 initialized to 0, as if the segment was unvoiced
iF0=0;
if ( rMax >ENERGY_THRESH*mean( vAutocorr ) ) 
		
    iT=iMax+iTmin-1;
   if iT~=0
      iF0=iFs/iT;
   end   
end