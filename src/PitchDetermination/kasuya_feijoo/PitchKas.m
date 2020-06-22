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
% Initial value of the index with which we will go through the sample vector. We take windows
% centered on the analysis point, but the index points to the beginning of the window.
iIndiceIni = iInicio-fix( iFrame/2 );

% Number of windows
iNumMuestras=iFinal-iInicio+1;
iNumVent=ceil( iNumMuestras/iDesplazamiento );

% Matrix C contains 3 rows and iNumVent columns. Each column corresponds
% to a segment, and contains the 3 candidates to be its fundamental period in samples.
C=zeros( 3, iNumVent );

for n=0:iNumVent-1
    
    iIndice=fix( iIndiceIni+n*iDesplazamiento);
    
    if (iIndice>=1) && ((iIndice+iFrame-1)<=length( vSignal ))
        
        % Rectangular windowed voice frame
        vFrame = vSignal( iIndice:(iIndice+iFrame-1) );
        
        C(:,n+1)=PitchSegKas( vFrame, iFs );
    end
end

% Obtaining T from C:

% First we obtain the median of the first candidates, contained in C1
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