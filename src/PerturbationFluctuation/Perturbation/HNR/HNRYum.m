function [vHNR, rHNR]=HNRYum( vSignal, iFs, iInicio, iFinal, iNumPuntos)

% Calculates the harmonic ratio to HNR noise in dB according to the Yumoto method
%   Yumoto, E., Gould, W. J., & Baer, T. (1982). Harmonics‐to‐noise ratio 
%   as an index of the degree of hoarseness. The journal of the Acoustical 
%   Society of America, 71(6), 1544-1550.
%
%   Yumoto, E., Sasaki, Y., & Okamura, H. (1984). Harmonics-to-noise ratio 
%   and psychophysical measurement of the degree of hoarseness. Journal of 
%   Speech, Language, and Hearing Research, 27(1), 2-6.
%
% Input parameters
%   vSignal:    Column vector containing the complete speech signal
%   iFs:        the sampling rate
%   iInicio:    Start the first sample of the section to be analyzed
%   iFinal:     the last sample of the section to be analyzed
%   iNumPoints: the number of HNR values returned
%
% Output parameters
%   vHNR:     Row vector of dimension "iNumPuntos" that contains the
%             instantaneous values of the HNR
%   rHNR:     the average value of the HNR

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

% The peaks of the section in a cycle-by-cycle analysis are obtained by the kasuya method,
% that is, with 40 ms windows and 50% overlap
try 
    [~, i_A, ~]=PicosMayores(vSignal, iFs, iInicio, iFinal, 1);
catch
    warning( 'It was not possible to calculate the HNR of the current signal, returning 0' )
    vHNR=zeros(1, iNumPuntos);
    rHNR=0;
    return
end

if isempty( i_A )
    vHNR=zeros(1, iNumPuntos);
    rHNR=0;
    return
end

% A matrix is obtained, whose rows contain the consecutive pitch periods of
% the voice signal, all of them with the duration of the maximum period, and filling in with
% zeros if necessary
C=SeparaCiclos(vSignal, iInicio, iFinal, i_A);

%HNRtotal=hnr_i_yum(C)

[n_periodos,periodoMax]=size( C );

iN=iFinal-iInicio+1;
iDesplazamiento=iN/iNumPuntos;

% Displacement in number of maximum pitch periods
desp_period=iDesplazamiento/periodoMax;

% The voice window size, on which each HNR value is to be calculated,
% in number of maximum pitch periods equals 60 (Boyanov recommends more than 50)
tam_vent_per=60;

% The initial index (in number of maximum periods) to go through the rows of C
% will be set to 1
indice_ini_per=1;

vHNR=zeros(1, iNumPuntos);

for n=0:iNumPuntos-1
   indice_per=fix(indice_ini_per+n*desp_period);
   
   if (n>0) && (indice_per==fix(indice_ini_per+(n-1)*desp_period))
      vHNR(n+1)=vHNR(n);
   else
      if (indice_per+tam_vent_per-1) <= n_periodos
          % The segment taken for the calculation of an HNR value is
         Ci=C(indice_per:indice_per+tam_vent_per-1, :);
         vHNR(n+1)=HNRiYumBoy( Ci );
      end
   end
end

% Now we find the average of the computed values to obtain the
% total HNR factor. The zero values are not included
rPromedio=0;
j=0;
for i=1:length(vHNR)
    if vHNR(i)~=0
        j=j+1;
        rPromedio=rPromedio+vHNR(i);
    end
end

rHNR=rPromedio/j;

if nargout == 0
   figure; plot(vHNR);
   title('HNR');
end
   