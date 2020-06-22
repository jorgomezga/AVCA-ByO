function [vHNR, rHNRmean]= CHNR( vSignal, iFs, iInicio, iFinal, iNumPuntos, iFo )

% Calculates the harmonic-to-noise ratio (HNR) of a voice signal, using
% the method of G. de Krom [1] based on the cepstrum.
%
% [1] Guus de Krom, "A cepstrum-based technique for determining a
% harmonics-to-noise ratio in speech signals ", J.Speech Hear.Res., vol.
% 36, no. 2, pp. 254-266, 1993.
%
% Input parameters
% vSignal:  Column vector containing the complete speech signal
% iFs:      the sampling rate
% i:        Start the first sample of the section to be analyzed
% iFinal:   the last sample of the section to be analyzed
% iNumPoints: the number of HNR values returned
% iFo:      Pitch frequency of the voice
%
% Output parameters
% rHNR:     Row vector of dimension "iNumPuntos" that contains the
%           instantaneous values of the HNR
% rHNRmean: the average value of the HNR

if nargin < 6, iFo = 0; end
if nargin < 5, iNumPuntos=100; end
if nargin < 4, iFinal=length( vSignal ); end
if nargin < 3, iInicio=1; end
if nargin < 2, error( 'Not enough input parameters!' ); end

% Check that the vector is of type column
if ~isvector( vSignal ) 
    error( 'vSignal is not a vector!' );
elseif size( vSignal, 2 ) ~= 1 
    vSignal=vSignal'; 
end

vHNR = zeros( 1, iNumPuntos ); 

% The HNR depends on the size of the chosen window. We must ensure that it
% has more than one signal period, so we either work with the
% longer period that is possible (15 ms) or we calculate the average period of the
% studied voice (Kasuya method; 40ms windows with 50% overlap).
if iFo == 0
    [rFo, ~] = PitchMedio(vSignal, iFs, iInicio, iFinal);
    iFo = iFs/rFo;
end

if iFo~=0
   iTamVent = 10*round( iFo );
else
   iTamVent = 5*round( 0.015*iFs );
end

% Window Size
iTamVent=2^( ceil( log2( iTamVent )));
iN = iFinal-iInicio+1; 

% The overlap is adjusted according to the provided value iNumPuntos
iResto = rem( iN, iNumPuntos-1 );
rSolape = (iTamVent - iResto) / (iNumPuntos-1);  
iAux = floor(iN/(iNumPuntos-1)); 
rAvance = iAux - rSolape; 

for n=0:iNumPuntos-1
    if n==0
        iIndice = iInicio;
    else
        iIndice = iInicio + fix( n*rAvance) - 1;
    end
    
    if (iIndice>=1) && ((iIndice+iTamVent-1)<=length( vSignal ))
        % vSegment is the voice without windowing (equivalent to rectangular windowing)
        vSegmento = vSignal( iIndice : (iIndice+iTamVent-1));
        vHNR(n+1) = CHNRi( vSegmento, iFs );
    end
end 

% Now we find the average of the values obtained to obtain the
% total HNR factor. The values zero are not accounted
rPromedio=0;
j=0;
for i=1:length(vHNR)
    if vHNR(i)~=0
        j=j+1;
        rPromedio=rPromedio+vHNR(i);
    end
end

rHNRmean=rPromedio/j;

if nargout == 0
   figure
   plot(vHNR);
   title('HNR');
end