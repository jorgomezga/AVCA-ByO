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

% Con 7 periodos de pitch se consigue suficiente resolución espectral. 
if iFo~=0
   iTamVent = 7*round( iFo );
else
   iTamVent = 7*round( 0.015*iFs );
end

iN = iFinal-iInicio+1; % Número de muestras de análisis


%%% Condicion para evitar error
if iTamVent>iN
    % Se ajusta el tamaño de la ventana para que quepan 3 ventanas al 50%
    % de traslape
    iTamVent = fix(2*(iFinal-iInicio)/(iNumPuntos+1));
end

if iNumPuntos==1
        vFrame = vSignal( iInicio:(1+iTamVent-1) );

        % Se calcula el pitch del segmento en estudio usando el m�todo de
        % Boyanov para el dominio temporal. Se pasa ToAnt=0.
        iFoi=PitchBoyTem( vFrame, iFs, 0 );
        if iFoi~=0 % Si el segmento es sonoro se calcula su NNE
            vFrame=vFrame-mean( vFrame );
            vFrame=vFrame.*hamming( iTamVent );

            vNNE=NNEi( vFrame, iFs, iFoi, [800;2500] );
        else  % Si el segmento es sordo la NNE toma un valor arbitrario nulo
            vNNE=0;
        end
        
        return
end



% Necesitamos iNumPuntos ventanas de análisis de tamaño iTamVent, por lo
% que ajustamos el solapamiento para que quepan en las iN muestras que
% tenemos.
iResto = rem( iN, iNumPuntos-1 );
rSolape = (iTamVent - iResto) / (iNumPuntos-1); % no se redondea para no acumular errores. 
iAux = floor(iN/(iNumPuntos-1)); % Tamaño que deberían tener las ventanas para caber sin solapar
rAvance = iAux - rSolape; % no se redondea para no acumular errores. 

for n=0:iNumPuntos-1
    if n==0
        iIndice = iInicio;
    else
        iIndice = iInicio + fix( n*rAvance ) - 1;
    end
    
    if iIndice+iTamVent-1 <= length( vSignal )
        % vFrame es la voz sin enventanar (equivale a enventanado rectangular)
        vFrame = vSignal( iIndice:(iIndice+iTamVent-1) );

        % Se calcula el pitch del segmento en estudio usando el método de
        % Boyanov para el dominio temporal. Se pasa ToAnt=0.
        iFoi=PitchBoyTem( vFrame, iFs, 0 );
        if iFoi~=0 % Si el segmento es sonoro se calcula su NNE
            vFrame=vFrame-mean( vFrame );
            vFrame=vFrame.*hamming( iTamVent );

           % vNNE(n+1)=NNEi( vFrame, iFs, iFoi, [800;2500] ); %valores
           % óptimos para detección de patología
            vNNE(n+1)=NNEi( vFrame, iFs, iFoi, [1000; 5000] );
        else  % Si el segmento es sordo la NNE toma un valor arbitrario nulo
            vNNE(n+1)=0;
        end
    end
end

% Hallamos el promedio de los valores obtenidos para calcular el factor
% NNE total. Los que hayan salido cero no los tenemos en cuenta
rPromedio=0;
j=0;
for i=1:length( vNNE )
    if vNNE(i)~=0
        j=j+1;
        rPromedio = rPromedio+vNNE(i);
    end
end

% Comprobamos que j sea mayor que cero.
if j>0
    rNNEmean = rPromedio/j;
else
    rNNEmean = 0;
end

if nargout == 0       
    figure; plot( vNNE );
    title('NNE según el metodo de Kasuya');
end 