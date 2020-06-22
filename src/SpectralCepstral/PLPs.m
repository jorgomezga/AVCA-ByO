function mParametros=PLPs( vSignal, iFs, sTipo, iNumCoef, iFrame, iSolape, eOptions, iVerbosity  )

% Calculate the PLP coefficients of a signal
%
% Simple use:
%    mParametros=PLPs(vSignal, iFs)	        % calculate PLP
%    mParametros=PLPs(vSignal, iFs, 'edD')  % include log energy, 1st and 2nd derivatives
%
% Input parameters:
%    vSignal  speech signal
%    iFs      sample rate in Hz (default 11025)
%    sTipo    any sensible combination of the following:
%        'R'  rectangular window in time domain
%        'N'  Hanning window in time domain
%        'M'  Hamming window in time domain (default)
%        'e'  include log energy
%        'L'  long term average
%        'r'  RASTA filtering
%        'd'  include delta coefficients (dc/dt)
%        'D'  include delta-delta coefficients (d^2c/dt^2)
%        'S'  voiced/unvoiced detection
%
%    iNumCoef      number of PLP coefficients (default 12)
%    iFrame        length of frame (default power of 2 <30 ms))
%    iSolape       frame increment (default iFrame/2)
%
% Output parameters:
%    mParametros   PLP output: one frame per row

%% Parameter check
if nargin<2, iFs=11025; end
if nargin<3, sTipo=''; end
if nargin<4, iNumCoef=12; end
if nargin<5, iFrame=pow2(floor(log2(0.03*iFs))); end
if nargin<6, iSolape=floor(iFrame/2); end
if nargin<7, eOptions=[]; end
if nargin<8, iVerbosity=0; end

% Check that the vector is of type column
if ~isvector( vSignal )
    error( 'vSignal is not a vector!' );
elseif size( vSignal, 2 ) ~= 1
    vSignal=vSignal';
end

%% Framing plus windowing
mSignal = enframe( vSignal, hamming( iFrame ), iSolape );

iNumVentanas = size( mSignal, 1 );
if iNumVentanas==0
    error( 'Empty Parameterization' )
end

mParametros = rastaplp_byo(vSignal, iFs, 0, iNumCoef-1, iFrame, iSolape, sTipo );

% mParametros  = zeros( iNumVentanas, iNumCoef );
% for i=1:iNumVentanas
%     if iVerbosity==1
%         fprintf('            +.+Processing window %d of %d\n',i,iNumVentanas)
%     end
%     vFrame=mSignal( i, : );
%     mParametros(i, :)=PLPi( vFrame, iFs, sTipo, iNumCoef );
% end

%% Plots
if nargout == 0
    % https://www.ee.columbia.edu/~dpwe/resources/matlab/rastamat/mfccs.html
    % Calculate iNumCoef order PLP features without RASTA
    %     [~, spec2] = rastaplp(vSignal, iFs, 0, iNumCoef);
    
    % Calculate HTK-style PLPs
    maxFreq = round(iFs/2);    
    
    plp = melfcc( vSignal, iFs, 'lifterexp', -22, 'nbands', 20, ...
        'dcttype', 1, 'maxfreq', maxFreq, 'fbtype', 'htkmel', 'preemph', 0, 'dither', 1,...
        'modelorder', 12, 'usecmp',1, 'wintime', iFrame/iFs, 'hoptime', iSolape/iFs);
        
%     [mm,~] = melfcc( vSignal, iFs, 'maxfreq', maxFreq, 'numcep', iNumCoef,...
%         'nbands', 22, 'fbtype', 'fcmel', 'dcttype', 1, 'usecmp', 1,...
%         'wintime', iFrame/iFs, 'hoptime', iSolape/iFs, 'preemph', 0, 'dither', 1 );
    
    [im,~] = invmelfcc( plp, iFs, 'lifterexp', -22, 'nbands', 20, ...
        'dcttype', 1, 'maxfreq', maxFreq, 'fbtype', 'htkmel', ...
        'modelorder', 12, 'usecmp', 1, 'wintime', iFrame/iFs, 'hoptime', iSolape/iFs);
    
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
    title('Reconstructed spectrogram from PLP coefficients')
    xlabel( 'Time (s)' );
    ylabel( 'Frequency (Hz)' );
    colorbar
    pretty_xyplot
    
end