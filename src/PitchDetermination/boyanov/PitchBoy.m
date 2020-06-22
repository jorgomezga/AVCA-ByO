function vTo=PitchBoy( vSignal, iFs, iInicio, iFinal )

% Calculate the pitch following Boyanov's recommendations.
%
% Input parameters:
%   vSignal: frame
%   iFs:     Sampling rate
%   iInicio: First sample of the section to be analyzed
%   iFinal:  The last sample of the section to be analyzed
%
% Output parameters:
%   vTo: the sequence of pitch period values (in samples) calculated using 
%        segments of length three pitch periods of the previous segment. 
%        The first segment lasts 34ms. There is no overlap of segments.
 
if nargin<4, iFinal=length(vSignal); end
if nargin<3, iInicio=1; end

% Check that the vector is of type column
if ~isvector( vSignal )
    error( 'vSignal is not a vector!' );
elseif size( vSignal, 2 ) ~= 1
    vSignal=vSignal';
end

vTo=[]; 

% We are going to segment the voice, from start to finish, in windows of
% size 3To, where To is the pitch period of the previous segment; the size of the
% first window will be 34 ms
iVentana=round(0.034*iFs);

iInfVentana=iInicio;
iSupVentana=iInfVentana + iVentana - 1;
rToAnt=0; 
n=1;

while iSupVentana<=iFinal
   
   % Frame using rectangular window
   vFrame = vSignal(iInfVentana:iSupVentana);
   
   rFo=PitchBoyanov( vFrame, iFs, rToAnt );
   
   % If sound
   if rFo~=0
       
      % vTo(n) is the pitch period in sample for the current segment. This value
      % will be used for the next segment
      vTo(n)=round(iFs/rFo);
      rToAnt=vTo(n);
      
      % Size of the next sample
      iVentana=3*vTo(n);
      n=n+1;
      
   % If the current segment is not sound, and at least one of the previous sound segment has been found,
   % the period of those segments is used to calculate the size of the next segment to
   % analyze
   elseif rToAnt~=0 
      
      vTo(n)=0;
      iVentana=3*rToAnt;
      n=n+1;
      
   % If not sound and no previous segments have been found, then the window 
   % size is kept (34ms)
   else
      vTo(n)=0;
      n=n+1;
   end
   
   % Next segment
   iInfVentana=iSupVentana+1;
   iSupVentana=iInfVentana+iVentana-1;      
end 

if nargout==0 
   figure; 
   plot(vTo/iFs);
   title('T0');
end