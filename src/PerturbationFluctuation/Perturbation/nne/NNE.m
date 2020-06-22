function [vNNE, rNNEmean]=NNE( vSignal, iFs, iInicio, iFinal, iNumPuntos, iFo )

% Calculates the normalized noise energy (NNE) of a speech segment (according to
% the Kasuya'86 method)
%
% Input parameters
% vSignal:  Column vector containing the full speech signal
% iFs:      Sample Rate (Hz)
% iStart:   First sample of the signal to be analyzed
% iFinal:   Last sample of the signal to be analyzed
% iNumPuntos:   Number of NNE values returned
% iFo:      Pitch frequency of the voice
%
% Output parameters
% rNNE:     row vecotr of dimension "iNumPuntos" that contains the
%           instantaneous values of the NNE
% rNNEmean: average value of the NNE

if nargin < 6, iFo = 0; end
if nargin < 5, iNumPuntos=100; end
if nargin < 4, iFinal=length( vSignal ); end
if nargin < 3, iInicio=1; end
if nargin < 2, error( 'Not enough input parameters!' ); end

% Check that the vector is of type column 
if ~isvector(vSignal)
    error( 'El parámetro vSignal no es un vector columna' ); 
elseif size( vSignal, 2 ) ~= 1
    vSignal=vSignal'; 
end 

vNNE = zeros( 1, iNumPuntos ); 

% The NNE depends on the size of the chosen window. We must ensure that it
% has more than one signal period, so we either work with the
% longer period that is possible (15 ms) or we calculate the average period of the
% studied voice (Kasuya method; 40ms windows with 50% overlap).
if iFo == 0
    [rFo, ~] = PitchMedio(vSignal, iFs, iInicio, iFinal);
    iFo = iFs/rFo;
end

% With 7 pitch periods sufficient spectral resolution is achieved.
if iFo~=0
   iTamVent = 7*round( iFo );
else
   iTamVent = 7*round( 0.015*iFs );
end

iN = iFinal-iInicio+1; 

%%% To avoid errors
if iTamVent>iN
    % Adjusts window size to fit 3 windows to 50% overlap
    iTamVent = fix(2*(iFinal-iInicio)/(iNumPuntos+1));
end

if iNumPuntos==1
        vFrame = vSignal( iInicio:(1+iTamVent-1) );

        % Pitch calculation using the Boyanov Method
        iFoi=PitchBoyanov( vFrame, iFs, 0 );
        if iFoi~=0 
            % Only calculated if the segment is sound
            vFrame=vFrame-mean( vFrame );
            vFrame=vFrame.*hamming( iTamVent );

            vNNE=NNEi( vFrame, iFs, iFoi, [800;2500] );
        else
            % If it is not sound, then an arbitrary value is given
            vNNE=0;
        end
        
        return
end

% We need iNumPoints analysis windows of iTamVent size, so that we adjust 
% the overlap to fit the iN samples that we have.
iResto = rem( iN, iNumPuntos-1 );
rSolape = (iTamVent - iResto) / (iNumPuntos-1); 
iAux = floor(iN/(iNumPuntos-1)); % Number of windows without overlapping
rAvance = iAux - rSolape; 

for n=0:iNumPuntos-1
    if n==0
        iIndice = iInicio;
    else
        iIndice = iInicio + fix( n*rAvance ) - 1;
    end
    
    if iIndice+iTamVent-1 <= length( vSignal )
        % Voice frame considering rectangular frames
        vFrame = vSignal( iIndice:(iIndice+iTamVent-1) );

        % Calculate pitch using Boyanov's method. Using ToAnt=0 for a start
        iFoi=PitchBoyanov( vFrame, iFs, 0 );
        if iFoi~=0 % Si el segmento es sonoro se calcula su NNE
            vFrame=vFrame-mean( vFrame );
            vFrame=vFrame.*hamming( iTamVent );

           % vNNE(n+1)=NNEi( vFrame, iFs, iFoi, [800;2500] ); %valores óptimos para detección de patología
            vNNE(n+1)=NNEi( vFrame, iFs, iFoi, [1000; 5000] );
        else              
            % If it is not sound, then an arbitrary value is given
            vNNE(n+1)=0;
        end
    end
end

% Mean calculation, ignoring the 0 values
rPromedio=0;
j=0;
for i=1:length( vNNE )
    if vNNE(i)~=0
        j=j+1;
        rPromedio = rPromedio+vNNE(i);
    end
end

% Check that j>0
if j>0
    rNNEmean = rPromedio/j;
else
    rNNEmean = 0;
end

if nargout == 0       
    figure; plot( vNNE );
    title('NNE using Kasuya method');
end 