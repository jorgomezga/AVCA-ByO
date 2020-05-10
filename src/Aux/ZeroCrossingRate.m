function iZCR = ZeroCrossingRate( vFrame )

% Calculate the zero-crossing rate of the frame
%
% Inputs:
%   vFrame               = Input vFrame
% Outputs:
%   iZCR                 = Zero-crossing rate

iZCR=0; 

iLengthFrame = length( vFrame );
if iLengthFrame > 0
   iZCR = 0.5*sum( abs( sign( vFrame(2:iLengthFrame) ) -...
       sign( vFrame( 1:iLengthFrame-1 ) ) ) );
   iZCR = iZCR/iLengthFrame;
end

return; 
