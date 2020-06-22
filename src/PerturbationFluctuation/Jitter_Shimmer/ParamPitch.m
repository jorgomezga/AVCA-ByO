function [vPPS, vPAS, vPitchCont, vAmpCont]=ParamPitch( vSignal, iFs, iInicio, iFinal, lMetodo)

% Calculates the pitch period and pitch amplitude sequences of the input voice
%
% Input parameters:
% "vSignal":     is a column vector containing the complete speech signal
% "iFs":         is the sample rate
% "iInicio":     is the first sample of the section to be analyzed
% "iFinal":      is the last sample of the section to be analyzed
% "iNumPuntos":  is the number of HNR values ​​that are returned
% "Method":      is the analysis method to be used: 0 (Boyanov), 1 (Kasuya).
%
% Output parameters:
% "vPPS":       sequence of signal pitch periods (cycle by cycle).
% "vPAS":       sequence of amplitudes of the peaks of the signal cycles.
%               (pitch amplitudes)
% "vPitchCont": sequence "contour of the pitch period" (of the same duration
%               that the section analyzed; the period in seconds of each cycle
%               of voice is repeated for all samples of the cycle)
% "vAmpCont":   sequence "contour of pitch amplitudes" (of the same duration
%               that the section analyzed; the amplitude of the peak of each cycle
%               of voice is repeated for all samples of the cycle)

if nargin < 5, lMetodo=0; end
if nargin < 4, iFinal=length( vSignal ); end
if nargin < 3, iInicio=1; end
if nargin < 2, error( 'Not enough input parameters!' ); end

[A, i_A, ~] = PicosMayores(vSignal, iFs, iInicio, iFinal, lMetodo);
iNumPicos   = length(A); 

% A contains the amplitude of the peaks that are obtained
vPAS        = A;

% In i_A are the positions, on the voice sequence s, of the peaks found,
% so calculating the distances between them will provide the sequence of periods
% of pitch in samples, and divided by fs in seconds. (As long as they are entirely
% sound). In addition, the contours of pitch and amplitudes are calculated.
% If there are silent segments in the considered voice segment, then, in
% those segments, no peaks will have been located, so the distance between the last peak before
% a silent segment and the first one afterwards, does not correspond to the value of a pitch period.
% We will consider invalid distances if they are greater than the maximum period allowed, 
% e.g 15 ms. In that case the contours take the arbitrary value 0.
vPPSmax  = 0.015;
vPPStemp = zeros(1,iNumPicos-1);

vPitchCont=zeros(1, iFinal-iInicio+1 );
vAmpCont=zeros(1, iFinal-iInicio+1 );

k=1;
for i=1:iNumPicos-1
   vPPStemp(i)=(i_A(i+1)-i_A(i))/iFs;
   if vPPStemp(i)<=vPPSmax
      % vPPS contains only the valid pitch periods
      vPPS(k)=vPPStemp(i);
      k=k+1;
      
      vPitchCont(i_A(i):i_A(i+1)-1)=vPPStemp(i);
      vAmpCont(i_A(i):i_A(i+1)-1)=vPAS(i);
   else
      vPitchCont(i_A(i):i_A(i+1)-1)=0;
      vAmpCont(i_A(i):i_A(i+1)-1)=0;
   end
end

% Finally extend the ends (first and last frame) of the contours calculated up to here, to
% occupy the entire interval from start to end.
vPitchCont( iInicio:i_A(1)-1 )=vPitchCont( i_A(1) );
vPitchCont( i_A(iNumPicos):iFinal )=vPitchCont( i_A(iNumPicos)-1 );

vAmpCont( iInicio:i_A(1)-1 )=vAmpCont( i_A(1) );
vAmpCont( i_A(iNumPicos):iFinal )=vAmpCont( i_A(iNumPicos)-1 );