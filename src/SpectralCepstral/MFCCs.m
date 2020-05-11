function mParametros=MFCCs( vSignal, iFs, sTipo, iNumCoef, iFrame,...
    iSolape, eOptions, iVerbosity )

% Calculate the mel cepstrum of a signal
%
% Simple use:
%    mParametros=MFCC(vSignal, iFs)	         % calculate mel cepstrum with 12 coefs, 256 sample frames
%    mParametros=MFCC(vSignal, iFs, 'e0dD')  % include log energy, 0th cepstral coef
%
% Input parameters:
%    vSignal  speech signal
%    iFs      sample rate in Hz (default 11025)
%    sTipo       any sensible combination of the following:
%        '0'  include 0'th order cepstral coefficient
%        'e'  include log energy
%        'E'  parameters are normalized with respect to the energy.
%
%    iNumCoef      number of cepstral coefficients excluding 0'th coefficient (default 12)
%    iFrame        length of frame (default power of 2 <30 ms))
%    iSolape       frame increment (default iFrame/2)
%    iNumFilters   number of filters in filterbank (default floor(3*log(iFs)) )
%    iLowFreq      low end of the lowest filter as a fraction of iFs (default = 0)
%    iHighFreq     high end of highest filter as a fraction of iFs (default = 0.5)
%    eOptions      structure including more parameters needed to parameterize.
%
% Output parameters
%   'mParametros'     mel cepstrum output: one frame per row

%% Check input parameteres
if nargin<2, iFs=11025; end
if nargin<3, sTipo='M'; end
if nargin<4, iNumCoef=12; end
if nargin<5, iFrame=pow2(floor(log2(0.03*iFs))); end
if nargin<6, iSolape=floor(iFrame/2); end
if nargin<7, eOptions=[]; end
if nargin<8, iVerbosity=0; end

%% Retrieve parameters
if ~isempty( eOptions ) && isfield( eOptions, 'iNumFilters' )
    iNumFilters = eOptions.iNumFilters;
else
    iNumFilters=floor(3*log(iFs));
end

if ~isempty( eOptions ) && isfield( eOptions, 'iLowFreq' )
    iLowFreq = eOptions.iLowFreq;
else
    iLowFreq = 0;
end

if ~isempty( eOptions ) && isfield( eOptions, 'iHighFreq' )
    iHighFreq = eOptions.iHighFreq;
else
    iHighFreq = 0.5;
end
        
if iLowFreq > iHighFreq 
    error( 'iLowFreq must be lower than iHighFreq' ); 
end

if iNumCoef<=0
    error( 'The number of MFCC coefficient should be positive and larger than 0' )
end

mSignal = enframe( vSignal, hamming( iFrame ), iSolape );

iNumVentanas = size( mSignal, 1 );
if iNumVentanas==0
    error( 'Empty Parameterization' )
end
   
%% MFCC parametrization on each one of the windows
if any( sTipo=='e' )
    mParametros  = zeros( iNumVentanas, iNumCoef+1 );
else
    mParametros  = zeros( iNumVentanas, iNumCoef );
end

for i=1:iNumVentanas
    if iVerbosity==1
        fprintf('            +.+Processing window %d of %d\n',i,iNumVentanas)
    end
    vFrame=mSignal( i, : ); 
    mParametros(i, :)=MFCCi( vFrame, iFs, sTipo, iNumCoef, iNumFilters, iLowFreq, iHighFreq  ); 
end


%% Plot
if nargout == 0
    
    % Convert to MFCCs very close to those genrated by feacalc    
    maxFreq = round(iFs/2);
    [mm,~] = melfcc( vSignal, iFs, 'maxfreq', maxFreq, 'numcep', iNumCoef,...
        'nbands', 22, 'fbtype', 'fcmel', 'dcttype', 1, 'usecmp', 1,...
        'wintime', iFrame/iFs, 'hoptime', iSolape/iFs, 'preemph', 0, 'dither', 1 );
    
    % .. then convert the cepstra back to audio (same options)
    [im,~] = invmelfcc(mm, iFs, 'maxfreq', maxFreq, 'numcep', iNumCoef,...
        'nbands', 22, 'fbtype', 'fcmel', 'dcttype', 1, 'usecmp', 1,...
        'wintime', iFrame/iFs, 'hoptime', iSolape/iFs, 'preemph', 0, 'dither', 1);
    
    % compare the spectrograms
    subplot(211)
    specgram(vSignal,512,iFs)
    caxis([-50 30])
    title('Original spectrogram')
    xlabel( '' ); 
    ylabel( 'Frequency (Hz)' ); 
    colorbar
    pretty_xyplot        
    
    subplot(212)
    specgram(im,512,iFs)
    caxis([-50 30])    
    title('Reconstructed spectrogram from MFCC coefficients')
    xlabel( 'Time (s)' ); 
    ylabel( 'Frequency (Hz)' ); 
    colorbar
    pretty_xyplot
    
    
    % Notice how spectral detail is blurred out e.g. the triangle hits around 6 kHz are broadened to a noise bank from 6-8 kHz.
    % save out the reconstruction
%     max(abs(im))
end