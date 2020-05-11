function vParametros=MFCCi( vFrame, iFs, sTipo, iNumCoef, iNumFilters, iLowFreq, iHighFreq )

% Calculate the mel-frequency cepstrum coefficients of a signal frame
%
% Simple use:
%    vParametros=MFCCi(vFrame, iFs)	      % calculate mel cepstrum with 12 coefs
%    vParametros=MFCCi(vFrame, iFs, 'e0') % include log energy, 0th cepstral coef
%
% Input parameters:
%    vFrame	 speech signal
%    iFs         sample rate in Hz (default 11025)
%    sTipo       any sensible combination of the following:
%        '0'  include 0'th order cepstral coefficient
%        'e'  include log energy
%        'E'  parameters are normalized with respect to the energy.
%    iNumCoef      number of MFCC coefficients  (default 12)
%    iNumFilters   number of filters in filterbank (default floor(3*log(iFs)) )
%    iLowFreq      low end of the lowest filter as a fraction of iFs (default = 0)
%    iHighFreq     high end of highest filter as a fraction of iFs (default = 0.5)
%
% Output parameters:
%   'vParametros'   mel cepstrum output

if nargin<2, iFs=50000; end
if nargin<3, sTipo='M'; end
if nargin<4, iNumCoef=12; end
if nargin<5, iNumFilters=floor(3*log(iFs)); end
if nargin<7
    iHighFreq=0.5;
    if nargin<6, iLowFreq=0; end
end

% Check that the vector is of type column
if ~isvector( vFrame )
    error( 'vSignal is not a vector!' );
elseif size( vFrame, 2 ) ~= 1
    vFrame=vFrame';
end

iNFFT=pow2( floor(log2( length(vFrame) ) ) );

% Transform to the frequency domain
vEspectro = rfft(vFrame.');

% Build the triangular filter bank.
[m,a,b]=melbankm(iNumFilters, iNFFT, iFs, iLowFreq, iHighFreq);%, sTipo);

% Normalize the amplitude of each filter by its width.
if any(sTipo=='v')
    vAnchos = sum(m~=0, 2) * ones(1,b-a+1);
    m = m./vAnchos;
end

% Calculate the energy of each filter (integration in mel bands).
if any(sTipo=='p')
    % Filters operate in the power domain.
    pot = vEspectro( a:b ).*conj( vEspectro( a:b ) );
    umbral_pot = max( pot )*1E-6;
    y = max( m*pot', umbral_pot);
else
    % Filters operate in the amplitude domain.
    amp = abs( vEspectro( a:b ) );
    umbral_amp = max( amp )*1E-3;
    y = max( m*amp', umbral_amp);
end

% Calculate the logarithm of the energy of each filter.
y = log(y);

% Transform to the cepstrum domain.
vParametros = rdct( y );

% Normalize between the number of filters.
vParametros = vParametros / iNumFilters;

iNumCoef=iNumCoef+1;
if iNumFilters>iNumCoef
    vParametros( iNumCoef+1:end )=[];
elseif iNumFilters<iNumCoef
    vParametros=[vParametros; zeros( iNumCoef-iNumFilters, 1)];
end

% Coefficient normalization with respect to the 0th coefficient
if any(sTipo=='E')
    for i=2:iNumCoef  % iNumCoef number of cepstral coefficients excluding 0'th coefficient
        vParametros(i)=vParametros(i)/(vParametros(1)+eps);
    end
end

% include 0'th order cepstral coefficient
if ~any(sTipo=='0')
    vParametros(1)=[];
end

% include log energy
if any(sTipo=='e')
    rLogEn=LogEnergia( vFrame );
    vParametros=[rLogEn; vParametros];
end