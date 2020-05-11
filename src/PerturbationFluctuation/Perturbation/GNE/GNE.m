function [vGNEValues, rGNEFactor]=GNE( vSignal, iFs, iInicio, iFinal,...
    iNumFrames, rFrameSize, iBW, iFShift )

% Calculates the GNE -Glottal-to-Noise Excitation ratio- according to
% the Michaelis algorithm [1]. The signal is resampled at 10 kHz,
% obtains the glottal waveform by reverse filtering and is divided into
% temporary windows. From each window the correlation between the
% Hilbert envelopes of different frequency bands is computed, being the GNE
% the maximum value.
%
% [1] D. Michaelis, T. Gramss and H. W. Strube, "Glottal-to-Noise Excitation
% ratio - a new measure for describing pathological voices ", Acustica / Acta
% Acoustic%, vol. 83, pp. 700-706, 1997.
%
% Input parameters:
%       vSenal: Vector that contains the samples of the voice signal whose
%               GNE is to be calculated.
%       iFs: Original sampling rate of the voice signal.
%       iStart: first sample of the section to be analyzed.
%       iFinal: last sample of the section to be analyzed.
%       iNumPtos: Number of temporary windows into which the signal is divided.
%               Indicates the number of calculated GNE values​ for each signal
%               voice. Must be a positive integer and sufficient to cover the
%               entire signal.
%       rFrameSize: Duration (in ms) of each analysis window.
%       iBW:    Bandwidth (Hz) of the frequency bands.
%               Recommended: 1000, 2000, 3000 Hz.
%       iFShift: Displacement (Hz) between the center frequencies of the
%               frequency bands. Recommended: 100, 200, 300 Hz.
% Output parameters:
%       vValoresGNE: row vector of length iNumPtos that stores the GNEs of
%           each time window.
%       rFactorGNE: Average value of the GNE parameter of the voice signal.

if nargin < 8, iFShift=100; end
if nargin < 7, iBW=3000;    end
if nargin < 6, rFrameSize=30e-3;         end
if nargin < 5, iNumFrames=20;            end
if nargin < 4, iFinal=length( vSignal ); end
if nargin < 3, iInicio=1;                end
if nargin < 2, error( 'Not enough input parameters!' ); end

if ((iInicio < 1) || (iFinal > length( vSignal )) || (iBW < 0) || (iFShift < 0) || (iFs < 0) )
   error('Wrong input arguments');
end

vGNEValues = zeros( 1, iNumFrames);
rGNEFactor = 0;

% Fixed sampling rate
iFsResample=10000;
% Frequency bands for the calculation of GNE (Hz).
vBand = [1000 3000];   

if iBW > (vBand(2)-vBand(1))
    disp(['Error: iBW has to be lower than the frequency margin: ', num2str(vBand(2)-vBand(1)), ' Hz']);
    return; 
end
if iBW > iFsResample/2
    disp(['Error: iBW has to be lower than the maximal frequency after resampling: ', num2str(iFsResample/2), ' Hz']);
    return; 
end

% Check that the vector is of type row
if ~isvector( vSignal ) 
    error( 'vSignal is not a vector!' );
elseif size( vSignal, 2 ) == 1 
    vSignal=vSignal'; 
end

vSignal=vSignal( iInicio:iFinal ); 

% ------------------------------------------------------------------------
% Resampling to 10 kHz.
% ------------------------------------------------------------------------
if iFs~=iFsResample
    vSignalResampled = resample( vSignal, iFsResample, iFs );
else
    vSignalResampled = vSignal;
end

% ------------------------------------------------------------------------
% Inverse filtering to obtain the fine grain structure of the signal 
% ------------------------------------------------------------------------
% Fixed parameters by [1].
iLPCOrder=13;           % Orden of the LPC fiter
rWindowLength = 30e-3;  % Length of the window (ms)
rIncrease     = 10e-3;  % increase between windows (ms).

vFineGrain = EstructuraFina( vSignalResampled, iFsResample, iLPCOrder, rWindowLength, rIncrease ); 

% ------------------------------------------------------------------------
% Temporal windowing of the fine grain structure
% ------------------------------------------------------------------------
% Size in samples of the window in temporal domain
iLongWindow = floor( iFsResample*rFrameSize ); 
if  ( iLongWindow + iNumFrames ) > length( vFineGrain )
    % Condition met when the length of the signal is not long enough to
    % fill a new frame
    iNumFrames= length( vFineGrain ) - iLongWindow ; 
end
% If the window is longer than the signal itself, we take only one window
if iLongWindow > length( vFineGrain )  
    iNumFrames = 1 ; 
    iLongWindow=length( vFineGrain ); 
end 

% We must calculate: the iDisplacement between windows after resampling
% (equal length to that of the fine structure) and for the
% Number of points or windows iNumPtos given, and the length of the window
% Original signal length v Fine structure, excitation
iLongSenal=length( vFineGrain );
    
if iNumFrames > 1
    % Displacement between temporary windows
    iDesplazamiento = floor ((iLongSenal-iLongWindow) / (iNumFrames-1));
    % Matrix that stores in each row a windowed frame of the
    % excitation or fine structure
    mVentanasEstFina=enframe( vFineGrain, iLongWindow, iDesplazamiento );   
    
else
    % Case of a single time frame: Taken centered in the middle of the signal    
    iInic=floor((iLongSenal-iLongWindow)/2);
    if iInic == 0
        iInic = 1; 
    end 
    mVentanasEstFina=vFineGrain( iInic:iInic+iLongWindow-1 );
end
 
% ------------------------------------------------- -----------------------
% Calculation of Hilbert envelopes.
% ------------------------------------------------- -----------------------
% The third index of the hypermatrix of dimension 3 indicates the number of
% temporary windows and for each one it contains an array whose rows are the
% Hilbert envelopes of each band (as many bands as rows).
hmEnvolventesHilbert = CalculaEnvolventesHilbert( mVentanasEstFina, iFsResample, iBW, iFShift, vBand );

% ------------------------------------------------- -----------------------
% Calculation of cross correlations.
% ------------------------------------------------- -----------------------
% Matrix with two columns with the indices of the pairs whose Hilbert envelopes
% have to be correlated.
mFutiles = CalculaParesHilbert( iLongWindow, iFsResample, iBW, iFShift, vBand );

% Number of pairs of frequency bands to be correlated.
iNumPares = length(mFutiles);

% Saves the maximums of the correlation of each pair.
vMaximos = zeros(iNumPares,1);

% For each temporal window
for i=1:iNumFrames
    for j=1:iNumPares
        % Michaelis reports the normalized correlation to 1 with deviation + -3 samples
        vCorr=xcorr( hmEnvolventesHilbert( mFutiles(j,1),:,i ), hmEnvolventesHilbert( mFutiles(j,2),:,i ), 3, 'coeff' );
        vMaximos(j) = max(vCorr);
    end
    % Vector row that stores the GNE of each temporary window (the maximum of
    % the maximums of the cross correlations).
    vGNEValues(1,i) = max(vMaximos);
end

% ------------------------------------------------------------------------
% GNE calculation
% ------------------------------------------------------------------------
% GNE of the voice signal as the average of the GNE of each time window.
% Those that are 0 are excluded.
rGNEFactor = mean( vGNEValues( vGNEValues~=0 ) ) ;

% Plotting of the correlation matrix
if nargout == 0
    iNumBandas=size( hmEnvolventesHilbert, 1 );  
    mCorr=zeros( iNumBandas, iNumBandas );   
    for p=1:iNumFrames
        for j=1:iNumBandas
            for i=1:iNumBandas
                vAux = xcorr( hmEnvolventesHilbert( j,:,p ), hmEnvolventesHilbert( i,:,p ), 3, 'coeff' );
                mCorr(i,j) = mCorr(i,j) + max( vAux ); 
            end
        end
    end
    mCorr=mCorr/iNumFrames; 
    mCorr=triu( mCorr );
    
    figure;    
    x = linspace( vBand(1)/1000, vBand(2)/1000, length(mCorr) );
    y = x;
    [X,Y] = meshgrid (x,y); % This generates the actual grid of x and y values.  
    axis( [vBand(1)/1000 vBand(2)/1000 vBand(1)/1000 vBand(2)/1000 0 1] )
    surf(X,Y,mCorr,'EdgeColor','none');
    xlabel( 'Frequency (kHz)' )
    ylabel( 'Frequency (kHz)' )
    shading interp 
    colorbar; 
    h = rotate3d;
    set(h,'ActionPostCallback',@align_axislabel)
    clc
end