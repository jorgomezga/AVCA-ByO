function T=PitchKas( vSignal, iFs, iInicio, iFinal)

% Calculate pitch using the Kasuya pitch determination algorithm
%
% Input parameters
%   vSignal: vector column containing the complete speech signal
%   iFs:     the sampling frequency
%   iInicio: Staring sample of the section to be analyzed
%   iFinal:  the last sample of the section to be analyzed
%
% Output parameters
%   T:       the sequence of pitch period values (in samples) calculated on
%            segments of 40ms, with displacements of 20ms, as Kasuya and Feijoo do.

if nargin<2, iFs=11025; end
if nargin<3, iInicio=1; end
if nargin<4, iFinal=length( vSignal) ; end

% Check that the vector is of type column
if ~isvector(vSignal)
    error( 'El parámetro vSignal no es un vector columna' );
elseif size( vSignal, 2 ) ~= 1
    vSignal=vSignal';
end

% iFrame is the window size (40 ms) in samples
iFrame=fix(0.04*iFs);

% The displacement between windows is 20ms (overlap 50%).
iDesplazamiento=0.02*iFs;

% Valor inicial del indice con el que recorremos el vector de muestras. Cogemos ventanas
% centradas en el punto de an�lisis, pero el �ndice apunta al principio de la ventana.
iIndiceIni = iInicio-fix( iFrame/2 );

% El n�mero de muestras de an�lisis es:
iNumMuestras=iFinal-iInicio+1;

% Por lo que el n�mero de ventanas ser�:
iNumVent=ceil( iNumMuestras/iDesplazamiento );

% La matriz C contiene 3 filas y iNumVent columnas. Cada columna corresponde
% a un segmento, y contiene los 3 candidatos a ser su periodo fundamental en muestras.
C=zeros( 3, iNumVent );

for n=0:iNumVent-1
    
    iIndice=fix( iIndiceIni+n*iDesplazamiento);
    
    if (iIndice>=1) && ((iIndice+iFrame-1)<=length( vSignal )),
        
        % Segmento es la voz sin enventanar (equivale a enventanado rectangular)
        vFrame = vSignal( iIndice:(iIndice+iFrame-1) );
        
        C(:,n+1)=PitchSegKas( vFrame, iFs );
    end
end

% Obtenci�n de T a partir de C:

% Primero se obtiene la mediana de los primeros candidatos, contenidos en C1
C1=C(1,:);
Tmediana=median(C1);

% For each segment, if the first candidate is within 20% of the median,
% the first candidate is taken as the pitch period of the segment.
% If not, the same is done with the following candidates. If neither meets the condition,
% the first is taken as first candidate.
margen=round(0.2*Tmediana);
minEntorno=Tmediana-margen;
maxEntorno=Tmediana+margen;
T = zeros(1, iNumVent);
for n=1:iNumVent
    if (C(1,n)>minEntorno) && (C(1,n)<maxEntorno)
        T(n)=C(1,n);
    elseif (C(2,n)>minEntorno) && (C(2,n)<maxEntorno)
        T(n)=C(2,n);
    elseif (C(3,n)>minEntorno) && (C(3,n)<maxEntorno)
        T(n)=C(3,n);
    else
        T(n)=C(1,n);
    end
end