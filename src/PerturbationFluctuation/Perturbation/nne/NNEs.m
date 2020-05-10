function [mParametros, vNNEmean]=NNEs( vSignal, iFs, sTipo, iFrame, iSolape )

% Calculates the normalized noise energy (NNE) of a voice signal,
% using the Kasuya method [1].
%
% [1] H. Kasuya, S. Ogawa, K. Mashima and S. Ebihara. "Normalized noise energy
% as an acoustic measure to evaluate pathologic voice ", J.Acoust.Soc.Am.,
% vol. 80, no. 5, pp. 1329-1334, 1986.
%
% Input parameters:
%    vSignal: column vector containing the full speech signal
%    iFs:     sample rate in Hz (default 50,000)
%    sType[not implemented]:   any reasonable combination of:
%           'L' long term average
%           'r' RASTA filtering
%           'd' first temporary derivative
%           'D' second temporal derivative
%           'S' voiced / unvoiced detection
%   iFrame:     size of the analysis window (default power of 2 <30 ms)
%   iSolape:    overlap between windows (default iFrame / 2)
%
% Output parameters
%   mArray:     matrix with as many rows as windows and as many columns
%               as frequency bands to be calculated.
%   vNNEmean:   the average value of each mNNE column

if nargin<2, iFs=50000; end
if nargin<3, sTipo='M'; end
if nargin<4, iFrame=pow2(floor(log2(0.03*iFs))); end
if nargin<5, iSolape=floor(iFrame/2); end

% Check that the vector is of type column 
if ~isvector(vSignal)
    error( 'El parámetro vSignal no es un vector columna' ); 
elseif size( vSignal, 2 ) ~= 1
    vSignal=vSignal'; 
end 

% Number of needed windows
iLongitud = length( vSignal ); 
iNumVentanas = fix((iLongitud - iFrame + iSolape)/(iSolape ));

% Lower and upper frequencies (in Hz) of each of the bands of
% frequency on which the NNE is calculated.
mBandas = [1000 1500 2000 2500 3000 3500;
           3000 3000 3000 3000 4000 4000];
iNumBandas = size(mBandas,2);

iInicio=1;
iFinal=length( vSignal );

% mNNE is an array of iNumBandas rows for the different frequency ranges
mNNE = zeros( iNumBandas, iNumVentanas ); 

% The NNE depends on the size of the chosen window. We must ensure that it
% has more than one signal period, so we either work with the
% longer period that is possible (15 ms) or we calculate the average period of the
% studied voice (Kasuya method; 40ms windows with 50% overlap).
rFo = PitchMedio(vSignal, iFs, iInicio, iFinal);
if rFo~=0
   iTamVent = 7*round( iFs/rFo );     % 7 pitch periods
else
   iTamVent = 7*round( 0.015*iFs );   % 7 times the expected minimum
end
% 7 pitch periods will provide enough spectral resolution

iN = iFinal-iInicio+1; % Number of samples of analysis

% We need iNumWindows of iTamVent size, and adjust the overlap to fit the iNSamples 
iResto = rem( iN, iNumVentanas-1 );
rSolape = (iTamVent - iResto) / (iNumVentanas-1); % no se redondea para no acumular errores. 
iAux = floor(iN/(iNumVentanas-1)); % Tamaño que deberían tener las ventanas para caber sin solapar
rAvance = iAux - rSolape; % no se redondea para no acumular errores. 

for n=0:iNumVentanas-1
    if n==0
        iIndice = iInicio;
    else
        iIndice = iInicio + fix( n*rAvance ) - 1;
    end
    
    if iIndice+iTamVent-1 <= length( vSignal )
        % Unwindowed voice (rectangular window)
        vFrame = vSignal( iIndice:(iIndice+iTamVent-1) );

        % The pitch of the segment under study is calculated using the
        % Boyanov pitch method. ToAnt = 0
        iFoi=PitchBoyTem( vFrame, iFs, 0 );
        if iFoi~=0 
            % If sound
            vFrame=vFrame-mean( vFrame );
            vFrame=vFrame.*hamming( iTamVent );

            mNNE(:,n+1)=NNEi( vFrame, iFs, iFoi, mBandas);
        else
            % If silent
            mNNE(:,n+1)=0;
        end
    end
end

% We find the average of the values obtained to calculate the factor
% NNE total. Those zero are not take them into account
vPromedio=zeros(iNumBandas,1);
j=0;

for i=1:length( mNNE )  % OJO con el uso de length
    if mNNE(:,i)~=0
        j=j+1;
        vPromedio(:)=vPromedio(:)+mNNE(:,i);
    end
end

if j>0
    vNNEmean = vPromedio/j;
else
    vNNEmean = 0;
end

% Rows: windows; Columns: each frequency band.
mParametros = mNNE'; 

% % ------------------------------------------------------------------------
% % Post processing
% % ------------------------------------------------------------------------
% 
% % Hacer filtrado RASTA.
% if any(sTipo=='r')
%     iFrameRate=iFs/iFrame*iSolape ;
%     mParametros = FiltroRasta( mParametros, iFrameRate );
% end;
% 
% % Calcular las derivadas primera y segunda.
% if any(sTipo=='D') || any(sTipo=='d'),
%   mParametros=Derivate( mParametros, sTipo );
% end
% 
% % Eliminar segmentos sordos.
% if any(sTipo=='S')
%     mSignal=enframe(vSignal, iFrame, iSolape);
%     vVoiceType=Vowel( mSignal, iFs ); 
%     % Find the positions of the elements that are non-zero
%     mParametros=mParametros( (vVoiceType == 1), : );
%     if sum( vVoiceType ) == 0
%         disp('El fichero especificado no contiene segmentos sonoros...'); 
%     end
% else 
%     % Remove the first and the last window
%     mParametros=mParametros( 2:iNumVentanas-1, : );   
% end
% 
% % Promediar los parámetros a largo plazo.
% if any(sTipo=='L')
%    mParametros=mean( mParametros, 2 );
% end

% ----------------------------------------------------------------------
if nargout<1
   [iNumParams, iNumVentanas]=size( mParametros );
   vTime=((0:iNumVentanas-1)*iSolape+(iFrame-1)/2)/iFs;
   vCoefs=(1:iNumParams)-any(sTipo=='0')-any(sTipo=='e');
   figure; imagesc( vTime, vCoefs, mParametros );
   axis('xy');
   xlabel('Time (signal)');
   ylabel('NNE coefficient');
   vMap = (0:63)'/63;
   colormap([ vMap vMap vMap]);
   colorbar;
end