function vParametros=PLPi( vFrame, iFs, sTipo, iNumCoef )

% Calculates the Perceptual linear prediction coefficients of a single frame
%
% Simple use:
%    vParametros=PLPi(vSignal, iFs)	     % calculate PLPs with 12 coefs
%    vParametros=PLPi(vSignal, iFs, 'e') % include log energy
%
% Parámetros de entrada: 
%    vFrame  	 speech signal
%    iFs         sample rate in Hz 
%    sTipo       any sensible combination of the following:
%        'r'  RASTA filtering
%        'e'  include log energy
%        'E'  parameters are normalized with respect to the energy.
%
%    iNumCoef      number of PLP coefficients (default 12)
%
% Parámetros de salida:
%   'vParametros'   mel cepstrum output

if nargin<2, iFs=50000; end
if nargin<3, sTipo='M'; end
if nargin<4, iNumCoef=12; end

if any(sTipo=='r')
    dorasta=1;
else
    dorasta=0;
end

% Check that the vector is of type column
if ~isvector( vFrame ) 
    error( 'vFrame is not a vector!' );
elseif size( vFrame, 2 ) ~= 1 
    vFrame=vFrame'; 
end

[vParametros, ~, ~, ~, ~, ~] = rastaplp(vFrame, iFs, dorasta, iNumCoef-1);

% Coefficient normalization with respect to the 0th coefficient
if any(sTipo=='E')
   for i=2:iNumCoef
       % iNumCoef number of cepstral coefficients excluding 0'th coefficient 
      vParametros(i)=vParametros(i)/(vParametros(1)+eps);    
   end
end

if any(sTipo=='e')               
   rLogEn=LogEnergy( vFrame );
   vParametros=[rLogEn; vParametros];
end